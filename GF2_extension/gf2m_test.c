#include "gf2m.h"
#include <time.h>

void test_gf2m_init()
{
    gf2m A;
    int m = 13;
    int t = 5;

    gf2m_init(&A, t);
    gf2m_random_gen(&A, m);
    gf2m_fit_len(&A);
    printf("len %d\n", A.deg);
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

void test_gf2m_mul_gf2_element()
{
    gf2m A, B, C, D;
    gf2 xgcd;
    gf2 mod;

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
    gf2m_init(&C, 1);
    gf2m_init(&D, 1);

    gf2m_random_gen(&A, m-1);

    printf("mod = ");       gf2_print(&mod);
    printf("A =");          gf2m_print(&A);

    xgcd = gf2m_monic(&B, &A, &mod);
    printf("monic of A = ");          gf2m_print(&B);
    printf("xcgd = ");      gf2_print(&xgcd);

    gf2m_mul_gf2_element(&C, &B, &A.term[t], &mod);

    printf("gf2m mul gf2 element = ");      gf2m_print(&C);

    gf2m_shift(&D, &C, 3);
    printf("shift << %d ", 3);      gf2m_print(&D);


}

void test_gf2m_division()
{
    gf2m A, B, Q, R;
    gf2 mod;

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
    gf2m_init(&B, t-1);
    gf2m_init(&Q, 1);
    gf2m_init(&R, 1);

    gf2m_random_gen(&A, m-1);
    gf2m_random_gen(&B, m-1);

    printf("mod = ");       gf2_print(&mod);
    printf("A =");          gf2m_print(&A);
    printf("B =");          gf2m_print(&B);
    gf2m_long_division(&Q, &R, &A, &B, &mod);
    printf("Q =");          gf2m_print(&Q);
    printf("R =");          gf2m_print(&R);

}

void test_gf2m_gcd() //확인 x
{
    gf2m A, B, gcd;
    gf2 mod;

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
    gf2m_init(&B, t-1);
    gf2m_init(&gcd, 1);

    gf2m_random_gen(&A, m-1);
    gf2m_random_gen(&B, m-1);

    printf("mod = ");       gf2_print(&mod);
    printf("A =");          gf2m_print(&A);
    printf("B =");          gf2m_print(&B);
    gf2m_gcd(&gcd, &A, &B, &mod);
    printf("gcd =");          gf2m_print(&gcd);

}

void test_gf2m_xgcd()
{
    gf2m A, B, ainv, binv, gcd;
    gf2 mod;

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
    gf2m_init(&B, t-1);
    gf2m_init(&gcd, 1);
    gf2m_init(&ainv, 1);
    gf2m_init(&binv, 1);

    gf2m_random_gen(&A, m-1);
    gf2m_random_gen(&B, m-1);

    printf("mod = ");       gf2_print(&mod);
    printf("A =");          gf2m_print(&A);
    printf("B =");          gf2m_print(&B);
    gf2m_xgcd(&gcd, &A, &B, &ainv, &binv, &mod);
    printf("gcd =");          gf2m_print(&gcd);
    printf("A Inv =");          gf2m_print(&ainv);
    printf("B Inv =");          gf2m_print(&binv);

}

void test_gf2m_square()
{
    gf2m A, B;
    gf2 mod;

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

    gf2m_square(&B, &A, &mod);

    printf("mod = ");               gf2_print(&mod);
    printf("A =");                  gf2m_print(&A);
    printf("square of A= ");        gf2m_print(&B);

}

void test_gf2m_powmod()
{
    gf2m A, B, gf2m_mod;
    gf2 mod;
    int index = 10;
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
    gf2m_init(&gf2m_mod, t-1);

    gf2m_random_gen(&A, m-1);
    gf2m_random_gen(&gf2m_mod, m-1);

    gf2m_repeated_squaremod(&B, &A, index, &gf2m_mod, &mod);

    printf("mod = ");               gf2_print(&mod);
    printf("gf2m_mod =");           gf2m_print(&gf2m_mod);
    printf("A =");                  gf2m_print(&A);
    printf("square of A^(2^%d)= ", index);        gf2m_print(&B);

}

void test_gf2m_squareroot()
{
    gf2m A, B, gf2m_mod;
    gf2 mod;
    int index = 10;
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
    gf2m_init(&gf2m_mod, t-1);

    gf2m_random_gen(&A, m-1);
    gf2m_random_gen(&gf2m_mod, m-1);

    //기약 다항식 생성 필요
    gf2m_square_root(&B, &A, &gf2m_mod, &mod);
    gf2m_repeated_squaremod(&B, &A, index, &gf2m_mod, &mod);

    printf("mod = ");               gf2_print(&mod);
    printf("gf2m_mod =");           gf2m_print(&gf2m_mod);
    printf("A =");                  gf2m_print(&A);
    printf("square of A^(2^%d)= ", index);        gf2m_print(&B);

}

void test_gf2_generate_irreducible_poly()
{
    gf2m A, gf2m_mod;
    gf2 mod;
    int m = 13;
    int t = 18; //큰 수는 행렬 연산에서 오래 걸림
    
    /* z^13+z^4+z^3+z^1+z^0, irreducible */
    gf2_init(&mod, m);
    gf2_set_index(&mod, 13);
    gf2_set_index(&mod, 4);
    gf2_set_index(&mod, 3);
    gf2_set_index(&mod, 1);
    gf2_set_index(&mod, 0);

    gf2m_init(&A, t);

    gf2m_generate_irreducible(&A, &mod, t);
    printf("irreducible polynomial : "); gf2m_print(&A);

}
void test_gf2m_is_irreducible()
{
    gf2m A;
    gf2 mod;
    int m = 13;
    int t = 8; //큰 수는 행렬 연산에서 오래 걸림
    int res;
    
    /* z^13+z^4+z^3+z^1+z^0, irreducible */
    gf2_init(&mod, m);
    gf2_set_index(&mod, 13);
    gf2_set_index(&mod, 4);
    gf2_set_index(&mod, 3);
    gf2_set_index(&mod, 1);
    gf2_set_index(&mod, 0);

    gf2m_init(&A, t);
    gf2m_random_gen(&A, m);
    printf("A = "); gf2m_print(&A);

    res = gf2m_is_irreducible(&A, &mod);
    if(res == IRREDUCIBLE)
        printf("is IRREDUCIBLE !\n");
    printf("is REDUCIBLE !\n");

}
void test_gf2m_math_operation(){
    //test_gf2m_init();
    //test_gf2m_add();
    //test_gf2m_mul();
    //test_gf2m_monic();
    //test_gf2m_mul_gf2_element();
    //test_gf2m_division();
    //test_gf2m_gcd();
    //test_gf2m_xgcd();
    //test_gf2m_square();
    //test_gf2m_powmod();
    //test_gf2_generate_irreducible_poly();
    test_gf2m_is_irreducible();

}
