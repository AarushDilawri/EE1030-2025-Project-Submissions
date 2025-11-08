#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image_write.h"
#include <stdio.h>
#include <stdlib.h>
#include "../c_libs/likho.h"


void likh(const char* basename, double** A, int w, int h) {
    unsigned char* out = (unsigned char*)malloc(w * h);
    double chhota = 1e9, bada = -1e9;

    // find range of pixels' values
    for (int i= 0; i < h; i++)
        for (int j = 0; j < w; j++) {
            if (A[i][j] <chhota) chhota = A[i][j];
            if (A[i][j] >bada) bada = A[i][j];
        }

    double range = (bada -chhota);
    if (range < 1e-12) range = 1.0;

    // normalize to 0,255
    for (int i = 0; i < h;i++)
        for (int j = 0; j < w; j++) {
            double v =(A[i][j] - chhota) / range * 255.0;
            out[i * w + j] = (unsigned char)(v + 0.5);
        }

    // Save 
    char pgm_bhai[300];
    snprintf(pgm_bhai, sizeof(pgm_bhai),"../../figs/%s.pgm", basename);
    FILE* f = fopen(pgm_bhai, "wb");
    if (!f){
        fprintf(stderr,  "image not opening %s for writing\n", pgm_bhai);
        free(out);
        return;
    }
    fprintf(f , "P5\n%d %d\n255\n", w, h);
    fwrite(out,1, w *h, f);
    fclose(f);
    printf("Saved PGM  %s\n", pgm_bhai);

    char jpg_bhai[300];
    snprintf(jpg_bhai, sizeof(jpg_bhai), "../../figs/%s.jpg",basename);
    stbi_write_jpg(jpg_bhai, w, h, 1, out, 90);
    printf("Saved JPEG  %s\n", jpg_bhai);
    free(out);
}

void likh_rgb(const char* basename, double*** A_rgb, int w, int h) {
    unsigned char* out = (unsigned char*)malloc(w *h * 3);
    double chhota = 1e9, bada = -1e9;

    // find pixel min,max 
    for (int c = 0; c < 3; c++)
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j++) {
                double v = A_rgb[c][i][j];
                if (v < chhota) chhota = v;
                if (v > bada) bada = v;
            }

    if (bada - chhota < 1e-12){
        bada =  chhota + 1.0;
    }

    // normalize 
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            for (int c = 0; c < 3; c++) {
                double v = (A_rgb[c][i][j] -chhota)/ (bada - chhota) * 255.0;
                if (v < 0) v = 0;
                if (v > 255) v = 255;
                out[(i * w + j) * 3 + c] =(unsigned char)(v+ 0.5);
            }
        }
    }

    char jpg_bhai[300];
    snprintf(jpg_bhai, sizeof(jpg_bhai), "../../figs/%s_rgb.jpg", basename);
    stbi_write_jpg(jpg_bhai, w, h, 3, out, 90);
    printf("Saved JPEG  %s\n", jpg_bhai);

    free(out);
}
 