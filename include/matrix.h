#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "gf2.h"

typedef struct{
    unsigned char *entries;
    int c;
    int r;
    int cnum;
    unsigned char **data;
}MAT_struct;

typedef MAT_struct BMAT[1];  //binary matrix

typedef struct
{
    gf2* entries;
    int c;
    int r;
    gf2 **data;
}gf2_MAT_struct;

typedef gf2_MAT_struct gf2_MAT[1];

#define b_mat_entry(mat,i,j) (*((mat)->data[i] + (j)))

void gf2_matrix_init(gf2_MAT mat, int rows, int cols, int m);
void gf2_matrix_print(const gf2_MAT mat);
void gf2_matrix_free(gf2_MAT mat);
gf2* gf2_mat_entry(gf2_MAT mat, int row, int col);
void generate_random_gf2_matrix(gf2_MAT mat);
void gf2_matrix_swap_rows(gf2_MAT mat, int row1, int row2);
int gf2_matrix_has_zero_rows(gf2_MAT mat);
void gf2_matrix_add_row(gf2_MAT mat, int row1, int row2);
void gf2_matrix_copy(gf2_MAT dst, gf2_MAT src);
int gf2_matrix_echelon(gf2_MAT mat_ech, gf2_MAT mat);


void bmatrix_init(BMAT mat, int rows, int cols);
void bmatrix_print(const BMAT mat);
void generate_random_bmatrix(BMAT mat);
void bmatrix_free(BMAT mat);
void bmatrix_swap_rows(BMAT mat, int row1, int row2);
int bmatrix_has_zero_rows(BMAT mat);
void bmatrix_add_row(BMAT mat, int row1, int row2);
void bmatrix_copy(BMAT dst, BMAT src);
int bmatrix_echelon(BMAT mat_ech, BMAT mat);


void test_matrix_operation();

#endif
