#define STB_IMAGE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#include <stdio.h>
#include <stdlib.h>
#include "../c_libs/padho.h"


double** padh(const char* f, int* w, int* h) {
    int ch; // channels
    unsigned char* img = stbi_load(f, w, h, &ch, 0);
    if (!img) {
        fprintf(stderr, "Error: Cannot load image %s\n", f);
        exit(1);
    }

    printf("Image read successfully: %d x %d\n", *w, *h);

    double** A = (double**)malloc((*h) * sizeof(double*));
    for (int i = 0; i < *h; i++)
        A[i] = (double*)malloc((*w) * sizeof(double));

    for (int y = 0; y < *h; y++) {
        for (int x = 0; x < *w; x++) {
            unsigned char r, g, b;
            if (ch >= 3) {
                r = img[(y * (*w) + x) * ch + 0];
                g = img[(y * (*w) + x) * ch + 1];
                b = img[(y * (*w) + x) * ch + 2];
                A[y][x] = 0.299 * r + 0.587 * g + 0.114 * b; //formula
            } else {
                A[y][x] = img[y * (*w) + x];
            }
        }
    }

    stbi_image_free(img);
    printf("Converted image to grayscale matrix.\n");
    return A;
}
double*** padh_rgb(const char* f, int* w, int* h) {
    int ch;
    unsigned char* img = stbi_load(f, w, h, &ch, 3); // force 3 channels
    if (!img) {
        fprintf(stderr, "Error: Cannot load image %s\n", f);
        exit(1);
    }

    printf("RGB image loaded successfully: %d x %d x %d\n", *w, *h, ch);

    // Allocate 3-channel double array
    double*** A = (double***)malloc(3 * sizeof(double**));
    for (int c = 0; c < 3; c++) {
        A[c] = (double**)malloc((*h) * sizeof(double*));
        for (int i = 0; i < *h; i++) {
            A[c][i] = (double*)malloc((*w) * sizeof(double));
        }
    }

    for (int y = 0; y < *h; y++) {
        for (int x = 0; x < *w; x++) {
            for (int c = 0; c < 3; c++) {
                A[c][y][x] = img[(y * (*w) + x) * 3 + c];
            }
        }
    }

    stbi_image_free(img);
    return A;
}