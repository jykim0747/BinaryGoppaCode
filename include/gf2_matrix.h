#ifndef _gf2_MATRIX_H_
#define _gf2_MATRIX_H_

#include "gf2.h"

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
void gf2_matrix_swap_rows(gf2_MAT mat, int row1, int row2);
int gf2_matrix_has_zero_rows(gf2_MAT mat);
void gf2_matrix_add_row(gf2_MAT mat, int row1, int row2);
void gf2_matrix_mul_row(gf2_MAT mat, int row, gf2* src, gf2* mod);
void gf2_matrix_copy(gf2_MAT dst, gf2_MAT src);
void gf2_matrix_copy_row(gf2_MAT dst, gf2_MAT src, int row);
int gf2_matrix_echelon(gf2_MAT mat_ech, gf2_MAT mat, gf2* mod);
void gf2_matrix_set_zero(gf2_MAT mat);
void gf2_matrix_generate_identity(gf2_MAT mat);

void test_gf2_matrix_operation();

#endif
