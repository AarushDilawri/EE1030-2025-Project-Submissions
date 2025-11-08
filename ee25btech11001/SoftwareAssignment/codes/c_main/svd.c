#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../c_libs/likho.h"
#include "../c_libs/svd.h"
#include "../c_libs/padho.h"


static const double zero = 1e-12;
static const int    max = 1000000;

//Helper functions
double** creatematrix(int m, int n){
    double ** res = (double**)malloc(m*sizeof(double*));
    for(int i=0; i<m; i++){
        res[i] = (double*)calloc(n, sizeof(double));
    }
    return res;
}
void free_matrix(double** M, int m) {
    if (!M){
        return;
    }
    for(int i = 0; i < m; i++){
        free(M[i]);
    }
    free(M);
}
double* createvec(int n){
    double* res = (double*)calloc(n, sizeof(double));
    return res;
}
void freevec(double* v){ 
    if(v){
        free(v);
    }
}

double innerprod(int n, const double* v1, const double* v2){
    double sum = 0.0;
    for(int i=0;i<n;i++){
        sum += v1[i]*v2[i];
    } 
    return sum;
}

double norm(int n, const double* v){
    double sum = 0.0;
    for(int i=0;i<n;i++){
        sum += v[i]*v[i];
    } 
    return sqrt(sum);
}
double** transpose(double** A, int m, int n){
    double** B = creatematrix(n,m);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[j][i] = A[i][j];
        }
    } 
    return B;
}

//this fn returns the result in R
void multPQ(int k, double** P, double** Q, double** R){
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++){
            double s = 0.0;
            for(int r=0;r<k;r++) s += P[i][r] * Q[r][j];
            R[i][j] = s;
        }
    }
}

//multiply matrix and a vector
double* multAb(int m, int n, double** A, const double* b){
    double* v = createvec(m);
    for(int i=0;i<m;i++){
        v[i] = innerprod(n, A[i], b);
    } 
    return v;
}

//since calculating G is O(n3), calculating Gx will better 
double* formGx(int m, int n, double** A, double** AT, const double* x){
    double* t = multAb(m, n, A, x);    
    double* Gx = multAb(n, m, AT, t);  
    freevec(t);
    return Gx;
}

//Take the projections of previous columns along the new column just calculated and make them 0
//doing reorthogonalize once wasnt working for all images, so did twice
void reorthogonalize(int n, int i, double* r, double** Q){
    for(int c=0; c<2; c++){
        for(int j=0; j<=i; j++){
            double proj = 0.0;
            for(int s=0; s<n; s++){
                proj += Q[s][j] * r[s];
            }
            if (proj != 0.0){
                for(int s=0; s<n; s++){
                    r[s] -= proj * Q[s][j];
                }
            }
        }
    }
}


//size of Q : n x (k+1), 
void lanczos(double** X, double** XT, int m, int n, int k, double* alpha, double* beta, double** Q){
    double rootn = 1.0 / sqrt((double)n);
    for(int i=0; i<n; i++){
        Q[i][0] = rootn;
    }
    //In lancsoz algo, we have to initialize a random vector for the algo to begin, I did this one, I also tried {1,0,0...0},
    //but it is sparse , so it was causing inaccuracies 

    double* col = createvec(n);
    double* r   = createvec(n);

    for(int i=0; i<k; i++){ //I have alrdy found q1, so running the loop k times to find {q2, ... qk+1}
        for(int t=0; t<n; t++){
            col[t] = Q[t][i];
        }

        double* Aq = formGx(m, n, X, XT, col);
        alpha[i] = innerprod(n, col, Aq);

        //r = Aqi - alphai*qi - beta(i-1)*q(i-1) (derivation done in report)
        for(int t=0;t<n;t++){
            r[t] = Aq[t] - alpha[i]*col[t];
        } 
        //beta0 = 0
        if (i > 0){
            for(int t=0; t<n; t++) {
                r[t] -= beta[i-1] * Q[t][i-1];
            }
        }

        reorthogonalize(n,i,r,Q);

        freevec(Aq);

        beta[i] = norm(n, r);
        if (beta[i] < zero) {
            for(int t=0;t<n;t++){
                Q[t][i+1] = 0.0;
            } 
            continue;
        }
        for(int t=0;t<n;t++){
            Q[t][i+1] = r[t] / beta[i];
        } 
    }

    freevec(col);
    freevec(r);
}


//Perform wilk shift on the pxp block from [0][0] to [p-1][p-1]
//the last 2x2 block is [[a,d], [d,b]]
// mu = d - (b^2/(delta^2 + root(delta^2 + b^2))), delta =(a-d)/2
double wilkinson(int p, double** T){
    if (p < 2){
        return 0.0;
    } 
    double a = T[p-2][p-2];
    double b = T[p-2][p-1];
    double d = T[p-1][p-1];
    double delta = (a - d) / 2.0;
    double mu = d - (b*b) / (fabs(delta) + sqrt(delta*delta + b*b));
    return mu;
}

void findeigen(int n, int k, double** T, double** V, double** Tp, double** Qp, double** QpT ){
    int p = k;
    int iter = 0;

    double* v = createvec(k); //column starting from diagonal element
    double* w = createvec(k); //Householder coefficients
    double* buff = createvec(k); 

    while (p > 1 && iter < max){
        iter++;
        //if bidiagonal elems are close to 0, make them 0
        if (fabs(T[p-1][p-2]) < zero) {
            T[p-1][p-2] = 0.0;
            T[p-2][p-1] = 0.0;
            p--;
            continue;
        }

        double mu = wilkinson(p, T);

        
        for(int i=0;i<p;i++){
            for(int j=0;j<p;j++){
                Tp[i][j] = T[i][j];
            }
            Tp[i][i] = T[i][i] - mu;
        }
        // T = T -muI
        //Qp = I
        for(int l=0; l<p; l++){
            Qp[l][l]= 1.0;
        }


        // Repeated QR Householder on Tp
        // Tp becomes R
        for(int l=0; l<p-1; l++){
            int size = p - l;
            
            for(int i=0;i<size;i++){
                v[i] = Tp[l+i][l];
            } 
            double alpha = norm(size, v);
            if (alpha < 1e-15){
                continue;
            } 
            double sign; 
            if(v[0] >= 0){
                sign = -1.0;
            }else{
                sign = 1.0;
            }
            v[0] += sign * alpha;
            double mag = norm(size, v);
            if (mag < 1e-15){
                continue;
            }
            for(int i=0;i<size;i++){
                v[i] = v[i]/mag;
            }

            //Now, norm(v) = 1
            // H = I - 2vvt/vtv  
            //since vtv = 1,
            //H = I - 2vvt
            //W have to pre multiply this H with Tp (currently its Tp - mu I)
            //So, HTp = Tp - 2v(vt Tp)
            //So first calculate 2vt Tp
            // We are starting from j=l because its a loop, initially l=0 was the case, so at that time we updated the whole
            // matrix T, we had j=0 to j=p, so we found p householder coffs. stored in w.
            // for any general l, the size of the matrix T is size*size (size = p-l), so, 
            for(int j=l;j<p;j++){
                double ans = 0.0;
                for(int i=0;i<size;i++){
                    ans += v[i] * Tp[l+i][j];
                } 
                w[j] = 2.0 * ans;
            }
            
            //Now, w is a horz vector
            //Now, multiply v and w, and subtract it from T
            for(int i=0; i<size; i++){
                for(int j=l;j<p;j++){
                    Tp[l+i][j] -= v[i] * w[j];
                }
            }

            // Now Qp = I
            // update it by post multiplying H, so for l iterations on Qp, we multiplied Qp by l Householder matrices
            // thus completing 1 QR
            for(int i=0;i<p;i++){
                double ans = 0.0;
                for(int r=0;r<size;r++){
                    ans += Qp[i][l+r] * v[r];
                } 
                double z = 2.0 * ans;
                for(int r=0;r<size;r++){
                    Qp[i][l+r] -= z * v[r];
                } 
            }
        }

        QpT = transpose(Qp, p, p);

        // compute QpT = Tp * QpT
        for(int i=0;i<p;i++){
            for(int j=0;j<p;j++){
                double s = 0.0;
                for(int r=0;r<p;r++){
                    s += Tp[i][r] * QpT[r][j];
                } 
                QpT[i][j] = s; 
            }
        }

        // T = QpT
        for(int i=0;i<p;i++){
            for(int j=0;j<p;j++){
                T[i][j] = QpT[i][j];
            } 
        } 

        // add mu back to diagonal
        for(int i=0;i<p;i++){
            T[i][i] += mu;
        } 
        //Now T =RQp +muI 
        // Now, V needs to be updated, V can be broken into to sub matrices
        // V = [V1 V2], V2 has the last p columns of V
        int start_col = k - p;
        for(int row = 0; row < n; row++){
            for(int c=0;c<p;c++){
                double s = 0.0;
                for(int r=0;r<p;r++){
                    s += V[row][start_col + r] * Qp[r][c];
                } 
                buff[c] = s;
            }
            for(int c=0;c<p;c++){
                V[row][start_col + c] = buff[c];
            } 
        }
        
    }

    freevec(v);
    freevec(w);
    freevec(buff);
}

extern void likh(const char* basename, double** A, int width, int height);

void compute_svd(double** A, int m, int n, int k, const char* basename){
    double** og_A = creatematrix(m,n);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            og_A[i][j] = A[i][j];
        } 
    } 

    double** AT = transpose(A, m, n);

    double* alpha = createvec(k);
    double* beta  = createvec(k);
    double** Ql = creatematrix(n, k+1);

    lanczos(A, AT, m, n, k, alpha, beta, Ql);
    //Now, Ql is lanczos basis
    //alpha has diagonal elems of T
    //beta has side-wale diagonal elements of T

    double** T = creatematrix(k,k);
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++) T[i][j] = 0.0;
        T[i][i] = alpha[i];
        if (i < k-1){
            T[i][i+1] = beta[i];
            T[i+1][i] = beta[i];
        }
    }

    double** V = creatematrix(n, k); // V is matrix containing columns as eigenvectors of G
    for(int i=0;i<n;i++){
        for(int j=0;j<k;j++){
            V[i][j] = Ql[i][j];
        } 
    } 

    double** Tp  = creatematrix(k,k);
    double** Qp  = creatematrix(k,k);
    double** QpT = creatematrix(k,k);

    findeigen(n, k, T, V, Tp, Qp, QpT);
    
    double* vals = createvec(k);
    for(int i=0;i<k;i++){
        vals[i] = sqrt(fabs(T[i][i]));
    } 

    //Sort
    for(int i=0;i<k-1;i++){
        for(int j=i+1;j<k;j++){
            if (vals[j] > vals[i]){
                double tmp = vals[i];
                vals[i] = vals[j];
                vals[j] = tmp;
                for(int r=0;r<n;r++){
                    double tmp2 = V[r][i];
                    V[r][i] = V[r][j];
                    V[r][j] = tmp2;
                }
            }
        }
    }

    //U = (og_A*V)/vals
    double** U = creatematrix(m, k);
    double* col = createvec(n);
    for(int j=0;j<k;j++){
        for(int r=0;r<n;r++){
            col[r] = V[r][j];
        } 
        double s = vals[j];
        if (s < 1e-14){
            for(int i=0;i<m;i++){
                U[i][j] = 0.0;
            } 
            continue;
        }
        for(int i=0;i<m;i++){
            U[i][j] = innerprod(n, og_A[i], col) / s;
        } 
        double mag = 0.0;
        for(int i=0;i<m;i++){
            mag += U[i][j]*U[i][j];
        } 
        mag = sqrt(mag);
        if (mag > 0.0){
            for(int i=0;i<m;i++){
                U[i][j] /= mag;
            } 
        } 
    }
    freevec(col);

    double** finalA = creatematrix(m, n);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int t=0;t<k;t++){
                finalA[i][j] += U[i][t] * vals[t] * V[j][t];
            }
             
        }
    }

    double frob = 0.0;
    double frob_A = 0.0;
    for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
        double d = A[i][j] - finalA[i][j];
        frob += d*d;
        frob_A += A[i][j]*A[i][j];
    }
    }
    frob = sqrt(frob);
    double rel = frob / (sqrt(frob_A) + 1e-30);
    printf("Frobenius error: %e  relative: %e\n", frob, rel);


    //grayscale bana do
    unsigned char *apna_image = malloc(m * n);
    double chhota = 1e9, bada = -1e9;

    // find min and max val a pixel can hld
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = finalA[i][j];
            if (v < chhota) chhota = v;
            if (v > bada) bada = v;
        }
    }

    if (bada - chhota < zero)
        bada = chhota + 1.0;

    // normalize to 0–255
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double val = (finalA[i][j] - chhota) * 255.0 / (bada - chhota);
            apna_image[i * n + j] = (unsigned char)(val + 0.5);
        }
    }



    char out[300];
    snprintf(out, sizeof(out), "../figs/%s", basename);
    likh(out, finalA, n, m);

    

    //clean up
    free(apna_image);
    free_matrix(finalA, m);

    free(alpha);
    free(beta);
    free(vals);

    free_matrix(AT, n);
    free_matrix(og_A, m);
    free_matrix(Ql, n);
    free_matrix(T, k);

    free_matrix(Tp, k);
    free_matrix(Qp, k);
    free_matrix(QpT, k);

    free_matrix(V, n);
    free_matrix(U, m);
}



void compute_svd_rgb(double*** A_rgb, int m, int n, int k, const char* basename) {
    const char* RGB[3] = {"_R", "_G", "_B"};
    char clrname[3][300];
    //created seperate grayscale for every color
    for (int c = 0; c < 3; c++){
        snprintf(clrname[c], sizeof(clrname[c]), "%s%s", basename, RGB[c]);
        printf("\n--- Compressing %s channel ---\n", RGB[c]);
        compute_svd(A_rgb[c], m, n, k, clrname[c]);
    }
    //merged those seperate grayscales
    double*** comp_rgb= (double***)malloc(3 * sizeof(double**));
    for (int c  =0; c <3; c++){
        char path[300];
        snprintf(path, sizeof(path), "../../figs/%s%s.pgm",basename,   RGB[c]);
        int w,h;
        comp_rgb[c] = padh(path, &w, &h);
    }

    //save rgb
    likh_rgb(basename, comp_rgb,n, m);

    // cleanup
    for (int c = 0;   c < 3; c++) {
        for (int i= 0; i < m; i++)
            free(comp_rgb[c][i]);
        free(comp_rgb[c]);
    }
    free(comp_rgb);
}
