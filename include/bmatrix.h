#ifndef _BMATRIX_H_
#define _BMATRIX_H_

// #include "gf2.h"
#include "gf2m.h"


typedef struct{
    unsigned char *entries;
    int c;
    int r;
    int cnum;
    unsigned char **data;
}MAT_struct;

typedef MAT_struct BMAT[1];  //binary matrix


#define b_mat_entry(mat,i,j) (*((mat)->data[i] + (j)))

void swap(int* a, int* b);

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
void bmatrix_generate_identity(BMAT mat);
void bmatrix_add_identity(BMAT dst, BMAT src);
int bmatrix_rank(BMAT mat);

int has_Identity_bmat(BMAT mat);
void make_Identity_bmat(BMAT mat, int *support);

void gf2m_to_bmat(BMAT mat, gf2m src, int m, int column);
int mat_concat_horizontal(BMAT dst, BMAT a, BMAT b);
int bmatrix_transpose(BMAT dst, BMAT src);
int bmatrix_mul(BMAT dst, BMAT src1, BMAT src2);
int bmatrix_inverse(BMAT dst, BMAT src);
void bmatrix_generate_inverse(BMAT dst, BMAT src);
int bmatrix_add(BMAT dst, BMAT src1, BMAT src2);
void bmatrix_slice(BMAT dst, BMAT src, int n);

void test_bmatrix_operation();

#endif
