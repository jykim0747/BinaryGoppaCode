#include "matrix.h"

void test_init_matrix()
{
    MAT A;
    int row = 3;
    int col = 3;
    int m = 3;

    gf2_matrix_init(A, row, col, m);
    generate_random_gf2_matrix(A);
    gf2_matrix_print(A);
    gf2_matrix_free(A);

}

void test_matrix_operation()
{
    test_init_matrix();


    
}