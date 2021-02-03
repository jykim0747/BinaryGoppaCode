#include "gf2m.h"
#include <time.h>

void test_gf2m_init()
{
    gf2m A;
    int m = 3;
    int t = 5;

    gf2m_init(&A, t);
    gf2m_random_gen(&A, m);
    gf2m_print(&A);
}

void test_gf2m_add()
{
    gf2m A, B, C;
    int m = 3;
    int t = 5;
    
    gf2m_init(&A, t);
    gf2m_init(&B, t);

    gf2m_random_gen(&A, m);
    gf2m_random_gen(&B, m);

    printf("A :");  gf2m_print(&A);
    printf("B :");  gf2m_print(&B);
    gf2m_add(&C, &A, &B);
    printf("A + B :");  gf2m_print(&C);
}

void test_gf2m_mul()
{
    gf2m A, B, C;
    gf2 mod;

    int m = 3;
    int t = 5;
    
    gf2_init(&mod, m);
    gf2_set_index(&mod, 3);
    gf2_set_index(&mod, 1);
    gf2_set_index(&mod, 0);

    gf2m_init(&A, t);
    gf2m_init(&B, t);
    gf2m_init(&C, 1);

    gf2m_random_gen(&A, m);
    gf2m_random_gen(&B, m);

    printf("mod = ");       gf2_print(&mod);
    printf("A =");          gf2m_print(&A);
    printf("B =");          gf2m_print(&B);
    gf2m_mul(&C, &A, &B, &mod);
    printf("A * B =");  gf2m_print(&C);
}

void test_gf2m_monic()
{
    gf2m A, B;
    gf2 mod;
    gf2 xgcd;

    int m = 13;
    int t = 5;
    
    /* z^13 + z^12 + z^10 + z^9 + z^7 + z^4 + 1, irreducible */
    gf2_init(&mod, m);
    gf2_set_index(&mod, 13);
    gf2_set_index(&mod, 12);
    gf2_set_index(&mod, 10);
    gf2_set_index(&mod, 9);
    gf2_set_index(&mod, 7);
    gf2_set_index(&mod, 4);
    gf2_set_index(&mod, 0);

    gf2m_init(&A, t);
    gf2m_init(&B, 1);

    gf2m_random_gen(&A, m-1);

    printf("mod = ");       gf2_print(&mod);
    printf("A =");          gf2m_print(&A);

    xgcd = gf2m_monic(&B, &A, &mod);
    printf("monic of A = ");          gf2m_print(&B);
    printf("xcgd = ");      gf2_print(&xgcd);

}

void test_gf2m_math_operation(){
    //test_gf2m_init();
    //test_gf2m_add();
    //test_gf2m_mul();
    test_gf2m_monic();
    
}
