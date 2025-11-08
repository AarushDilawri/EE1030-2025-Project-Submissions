#define STB_IMAGE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#include <stdio.h>
#include <stdlib.h>
#include "../c_libs/padho.h"


double** padh(const char* f, int* w, int* h) {
    int ch; // channels
    unsigned char* img   = stbi_load(f, w, h, &ch, 0);
    if (!img) {
        fprintf(stderr, "image not loading %s\n", f);
        exit(1);
    }

    printf("Image read ho gayi %d x %d\n", *w, *h);

    double** A =  (double**)malloc((*h) * sizeof(double*));
    for (int i = 0; i < *h; i++)
        A[i] = (double*)malloc((*w) * sizeof(double));

    for (int y = 0; y <*h; y++) {
        for (int x = 0; x< *w; x++) {
            unsigned char r, g, b;
            if (ch >= 3) {
                r = img[(y *(*w) + x) * ch + 0];
                g =img[(y *(*w) +x) * ch +1];
                b = img[(y *(*w) + x) * ch + 2];
                A[y][x] = 0.299 *r + 0.587 *g + 0.114 *b; //formula
            } else {
                A[y][x] = img[y *(*w) + x];
            }
        }
    }

    stbi_image_free(img);
    printf("Image to matrix done\n");
    return A;
}
double*** padh_rgb(const char* f , int* w, int* h) {
    int ch;
    unsigned char* img = stbi_load(f, w, h, &ch, 3); 
    if (!img) {
        fprintf(stderr, "image not loading %s\n", f);
        exit(1);
    }

    printf("RGB img lo: %d x %d x %d\n", *w, *h, ch);

    // make a 3 channel array
    double*** A = (double***)malloc(3 * sizeof(double**));
    for (int c =0; c < 3; c++) {
        A[c] =(double**)malloc((*h) * sizeof(double*));
        for (int i = 0; i < *h; i++) {
            A[c][i] = (double*)malloc((*w)* sizeof(double));
        }
    }

    for (int y = 0; y <*h; y++) {
        for (int x =0; x < *w; x++) {
            for (int c =0; c < 3; c++) {
                A[c][y][x] = img[(y * (*w)+ x) * 3 +c];
            }
        }
    }

    stbi_image_free(img);
    return A;
}