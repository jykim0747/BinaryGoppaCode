#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "gf2.h"

typedef struct{
    int* entries;
    int c;
    int r;
    int **data;
}MAT_struct;

typedef MAT_struct MAT[1];  //binary로 다시 만들 예정

typedef struct
{
    gf2* entries;
    int c;
    int r;
    gf2 **data;
}gf2_MAT_struct;

typedef gf2_MAT_struct gf2_MAT[1];


void gf2_matrix_init(gf2_MAT mat, int rows, int cols, int m);
void gf2_matrix_print(const gf2_MAT mat);
void gf2_matrix_free(gf2_MAT mat);
gf2* gf2_mat_entry(gf2_MAT mat, int row, int col);
void generate_random_gf2_matrix(gf2_MAT mat);


void test_matrix_operation();

#endif
