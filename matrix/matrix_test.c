#include "bmatrix.h"
#include "gf2_matrix.h"

void test_init_matrix()
{
    gf2_MAT A, B;
    int row = 7;
    int col = 7;
    int m = 3;

    gf2_matrix_init(A, row, col, m);
    generate_random_gf2_matrix(A);
    gf2_matrix_print(A);
    gf2_matrix_swap_rows(A, 0, 1);
    gf2_matrix_print(A);
    gf2_matrix_free(A);

    gf2_matrix_init(B, row, col, m);
    gf2_matrix_generate_identity(B);
    gf2_matrix_print(B);
    gf2_matrix_free(B);
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

void test_echelon_form_matrix()
{
    gf2_MAT A, B;
    gf2 mod;
    int row = 5;
    int col = 7;
    int m = 3;
    
    gf2_init(&mod, m);
    gf2_generate_irreducible(&mod, m);

    gf2_matrix_init(A, row, col, m - 1);
    gf2_matrix_init(B, row, col, m - 1);
    
    generate_random_gf2_matrix(A);
    gf2_matrix_print(A);

    printf("echelon form\n");
    gf2_matrix_echelon(B, A, &mod);
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

void test_gf2_to_bmat()
{
    BMAT A;
    gf2 B;

     gf2_init(&B, 17);
    // gf2_set_index(&B, 4);
    // gf2_set_index(&B, 2);
    // gf2_set_index(&B, 1);
    // gf2_set_index(&B, 0);
    gf2_random_gen_fix(&B);

    printf("A = "); gf2_print(&B);

    bmatrix_init(A, B.deg, B.deg);
    bmatrix_set_gf2(A, B, 1);
    bmatrix_print(A);

    bmatrix_free(A);
}

/**
 * 세로
 */
void test_gf2m_to_bmat()
{
    BMAT A;
    gf2m B;
    int t = 3;
    int m = 5;

    gf2m_init(&B, t);
    gf2m_random_gen(&B, m);
    gf2m_fit_len(&B);

    printf("A = "); gf2m_print(&B);

    bmatrix_init(A, (t+1)*(m+1), 16);
    gf2m_to_bmat(A, B, m+1, 1);
    bmatrix_print(A);

    bmatrix_free(A);
}


void test_generate_identity_bmatrix()
{
    BMAT A;
    int row = 7;
    int col = 7;
    int rank;

    bmatrix_init(A, row, col);
    //bmatrix_generate_identity(A);
    generate_random_bmatrix(A);
    bmatrix_print(A);
    rank = bmatrix_rank(A);
    printf("rank = %d\n", rank);
    bmatrix_free(A);
}

void test_concat_horizontal()
{
    BMAT A, B, C;
    int row = 2;
    int col = 7;

    bmatrix_init(A, row, col);
    bmatrix_init(B, row, 20);
    bmatrix_init(C, row, A->c + B->c);

    generate_random_bmatrix(A); bmatrix_print(A);
    generate_random_bmatrix(B); bmatrix_print(B);

    printf("C = [A || B]");
    mat_concat_horizontal(C, A, B);
    bmatrix_print(C);

    bmatrix_free(A);
    bmatrix_free(B);
    bmatrix_free(C);
}

void test_matrix_transpose()
{
    BMAT A, B, ech;
    int row = 12;
    int col = 16;

    bmatrix_init(A, row, col);
    bmatrix_init(ech, row, col);
    bmatrix_init(B, A->c, A->r);

    generate_random_bmatrix(A); bmatrix_print(A);
    bmatrix_echelon(ech, A); bmatrix_print(ech);

    bmatrix_transpose(B, ech);
    printf("B = [A]^T \n");
    bmatrix_print(B);

    bmatrix_free(A);
    bmatrix_free(B);
    bmatrix_free(ech);
}

void test_matrix_mul()
{
    BMAT A, B, C;
    int row = 4;
    int col = 16;
    
    bmatrix_init(A, row, row);
    bmatrix_init(B, row, col);
    bmatrix_init(C, row, col);
    
    generate_random_bmatrix(A);
    bmatrix_print(A);

    generate_random_bmatrix(B);
    bmatrix_print(B);

    bmatrix_mul(C, A, B);
    printf("C =  A * B\n");
    bmatrix_print(C);


    bmatrix_free(A);
    bmatrix_free(B);

}

void test_matrix_inverse()
{
    BMAT A, B;
    int row = 10;
    int col = 10;
    
    bmatrix_init(A, row, col);
    bmatrix_init(B, col, col);
    
    bmatrix_generate_inverse(B, A);

    printf("A \n");
    bmatrix_print(A);

    printf("B =  inverse of A \n");
    bmatrix_print(B);

    bmatrix_free(A);
    bmatrix_free(B);

}

void test_matrix_slice()
{
    BMAT A, B;
    int row = 3;
    int col = 28;
    int n = 20;
    
    bmatrix_init(A, row, col);
    
    generate_random_bmatrix(A);

    printf("A \n");
    bmatrix_print(A);

    bmatrix_slice(B, A, n);
    printf("B =  %d bits of A \n", n);
    bmatrix_print(B);

    bmatrix_free(A);
    bmatrix_free(B);

}

void test_matrix_mat_2_gf2m()
{
    BMAT A;
    gf2m poly;
    int t = 3;
    int m = 4;
    int row = 12;
    int col = 1;
    
    bmatrix_init(A, row, col);    
    generate_random_bmatrix(A);

    printf("A \n");
    bmatrix_print(A);

    gf2m_init(&poly, t);
    bmatrix_to_gf2m(&poly, A, m, t);

    bmatrix_free(A);

}


void test_gf2_matrix_operation()
{
    test_init_matrix();
    //test_has_zero_matrix();
    //test_copy_matrix();
    //test_echelon_form_matrix();

}

void test_bmatrix_operation()
{
    //test_init_bmatrix();
    //test_has_zero_bmatrix();
    //test_copy_bmatrix();
    //test_echelon_form_bmatrix();
    // test_gf2_to_bmat();
    // test_generate_identity_bmatrix();
    //test_gf2m_to_bmat();
    //test_concat_horizontal();
    //test_matrix_transpose();
    //test_matrix_mul();
    //test_matrix_inverse();
    //test_matrix_slice();
    test_matrix_mat_2_gf2m();
}