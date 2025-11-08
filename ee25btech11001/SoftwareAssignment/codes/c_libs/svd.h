#ifndef SVD_H
#define SVD_H

void compute_svd(double** A, int m, int n, int k, const char* basename);
void compute_svd_rgb(double*** A_rgb, int m, int n, int k, const char* basename);


#endif
