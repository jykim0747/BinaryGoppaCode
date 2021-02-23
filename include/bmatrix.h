#ifndef _BMATRIX_H_
#define _BMATRIX_H_

#include "gf2.h"

typedef struct{
    unsigned char *entries;
    int c;
    int r;
    int cnum;
    unsigned char **data;
}MAT_struct;

typedef MAT_struct BMAT[1];  //binary matrix


#define b_mat_entry(mat,i,j) (*((mat)->data[i] + (j)))


void bmatrix_init(BMAT mat, int rows, int cols);
void bmatrix_print(const BMAT mat);
void generate_random_bmatrix(BMAT mat);
void bmatrix_free(BMAT mat);
void bmatrix_swap_rows(BMAT mat, int row1, int row2);
int bmatrix_has_zero_rows(BMAT mat);
void bmatrix_add_row(BMAT mat, int row1, int row2);
void bmatrix_copy(BMAT dst, BMAT src);
int bmatrix_echelon(BMAT mat_ech, BMAT mat);

void bmatrix_set_gf2(BMAT mat, const gf2 src, int row);

void test_bmatrix_operation();

#endif
