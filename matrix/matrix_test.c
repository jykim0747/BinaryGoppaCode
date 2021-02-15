#include "matrix.h"

void test_init_matrix()
{
    gf2_MAT A;
    int row = 2;
    int col = 14;
    int m = 3;

    gf2_matrix_init(A, row, col, m);
    generate_random_gf2_matrix(A);
    gf2_matrix_print(A);
    gf2_matrix_swap_rows(A, 0, 1);
    gf2_matrix_print(A);
    gf2_matrix_free(A);
}

void test_has_zero_matrix()
{
    gf2_MAT A;
    int row = 2;
    int col = 14;
    int m = 3;
    int res;

    gf2_matrix_init(A, row, col, m);
    generate_random_gf2_matrix(A);
    gf2_matrix_print(A);
    
    for(int i=0; i<A->c; ++i)
        gf2_set_zero(gf2_mat_entry(A, 0, i));

    gf2_matrix_print(A);
    res = gf2_matrix_has_zero_rows(A);
    printf("has zero %d\n", res);
    printf("add row \n");

    gf2_matrix_add_row(A, 0, 1);

    gf2_matrix_print(A);
    gf2_matrix_free(A);

}

void test_copy_matrix()
{
    gf2_MAT A, B;
    int row = 3;
    int col = 7;
    int m = 3;
    
    gf2_matrix_init(A, row, col, m);
    gf2_matrix_init(B, row, col, m);
    
    generate_random_gf2_matrix(A);
    gf2_matrix_print(A);
    gf2_matrix_copy(B, A);
    gf2_matrix_print(B);

    gf2_matrix_free(A);
    gf2_matrix_free(B);
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

void test_has_zero_bmatrix()
{
    BMAT A;
    int row = 2;
    int col = 14;
    int res;
    bmatrix_init(A, row, col);
    generate_random_bmatrix(A);
    bmatrix_print(A);
    b_mat_entry(A, 0, 0) = 0;
    b_mat_entry(A, 0, 1) = 0;
    bmatrix_print(A);
    res = bmatrix_has_zero_rows(A);
    printf("has zero %d\n", res);
    bmatrix_add_row(A, 0, 1);
    bmatrix_print(A);
    bmatrix_free(A);

}

void test_copy_bmatrix()
{
    BMAT A, B;
    int row = 3;
    int col = 7;
    
    bmatrix_init(A, row, col);
    bmatrix_init(B, row, col);
    
    generate_random_bmatrix(A);
    bmatrix_print(A);
    bmatrix_copy(B, A);
    bmatrix_print(B);

    bmatrix_free(A);
    bmatrix_free(B);
}

void test_echelon_form_bmatrix()
{
    BMAT A, B;
    int row = 5;
    int col = 7;
    
    bmatrix_init(A, row, col);
    bmatrix_init(B, row, col);
    
    generate_random_bmatrix(A);
    bmatrix_print(A);

    printf("echelon form\n");
    bmatrix_echelon(B, A);
    bmatrix_print(B);

    bmatrix_free(A);
    bmatrix_free(B);
}

void test_matrix_operation()
{
    //test_init_matrix();
    //test_has_zero_matrix();
    //test_copy_matrix();

    //test_init_bmatrix();
    //test_has_zero_bmatrix();
    //test_copy_bmatrix();
    test_echelon_form_bmatrix();

}