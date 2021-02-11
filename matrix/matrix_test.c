#include "matrix.h"

void test_init_matrix()
{
    gf2_MAT A;
    int row = 3;
    int col = 3;
    int m = 3;

    gf2_matrix_init(A, row, col, m);
    generate_random_gf2_matrix(A);
    gf2_matrix_print(A);
    gf2_matrix_swap_rows(A, 0, 1);
    gf2_matrix_print(A);
    gf2_matrix_free(A);
}

void test_init_bmatrix()
{
    BMAT A;
    int row = 2;
    int col = 14;

    bmatrix_init(A, row, col);
    generate_random_bmatrix(A);
    bmatrix_print(A);
    bmatrix_swap_rows(A, 0, 1);
    bmatrix_print(A);
    bmatrix_free(A);
}

void test_matrix_operation()
{
    test_init_matrix();
    //test_init_bmatrix();

    
}