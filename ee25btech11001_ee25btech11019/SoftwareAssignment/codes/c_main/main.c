#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../c_libs/padho.h"
#include "../c_libs/svd.h"

int main() {
    char filename[256], outname[256], filepath[300];
    int k, mode;

    printf("Enter image filename (e.g. photo.jpg): ");
    scanf("%s", filename);
    printf("Enter value of k (compression rank): ");
    scanf("%d", &k);
    printf("Enter output base name: ");
    scanf("%s", outname);
    printf("Enter 1 for grayscale, 3 for RGB: ");
    scanf("%d", &mode);

    snprintf(filepath,sizeof(filepath), "../../figs/%s", filename);

    int w, h;

    if (mode == 1) {
        double** A = padh(filepath, &w, &h);
        compute_svd(A, h,w, k, outname);
        for (int i = 0; i < h; i++){
            free(A[i]);
        } 
        free(A);
    } else if (mode == 3) {
        double*** A_rgb = padh_rgb(filepath, &w, &h);
        compute_svd_rgb(A_rgb, h, w, k, outname);
        for (int c = 0; c< 3; c++) {
            for (int i =0; i< h; i++){
                free(A_rgb[c][i]);
            } 
            free(A_rgb[c]);
        }
        free(A_rgb);
    } else {
        printf("wrong\n");
    }

    printf("done\n");
    return 0;
}
