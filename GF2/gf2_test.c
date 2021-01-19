#include "gf2.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void test_gf2_init();
void test_gf2_add();
void test_gf2_mul();

void test_gf2_math_operation();
void test_gf2_div();
void test_gf2_square();
void test_gf2_pow();
void test_gf2_mulmod();


void test_gf2_init()
{
    gf2 A, B;
    int t = 9;

    gf2_init(&A, t);
    gf2_init(&B, t);
    gf2_random_gen(&A);
    //A.binary[0] = 0x00;
    //A.binary[1] = 0x00;
    gf2_print(&A);
    
    printf(" fix %d \n", t);
    gf2_random_gen_fix(&B);
    gf2_print(&B);

}

void test_gf2_add()
{
    gf2 A, B, C;
    int t = 20;

    gf2_init(&A, t);
    gf2_random_gen(&A);

    gf2_init(&B, t/2);
    gf2_random_gen(&B);


    gf2_add(&C, &A, &B);

    printf("A = ");               gf2_print(&A);
    printf("B = ");               gf2_print(&B);
    printf("A + B = ");           gf2_print(&C);

}

void test_gf2_mul()
{
    gf2 A, B, C;
    int t = 20;
    int index = 7;

    gf2_init(&A, t);
    gf2_random_gen(&A);

    printf("A = ");               gf2_print(&A);
    gf2_mul_shift(&C, &A, index);
    printf("A << %d = ", index);  gf2_print(&C);

    gf2_init(&B, t/2);
    gf2_random_gen(&B);

    gf2_mul(&C, &B, &A);
    printf("A = ");               gf2_print(&A);
    printf("B = ");               gf2_print(&B);
    printf("A * B = ");           gf2_print(&C);

}

void test_gf2_div()
{
    gf2 A, B, Q, R;
    int m = 13;

    gf2_init(&A, m);
    gf2_random_gen(&A);

    gf2_init(&B, 8);
    gf2_random_gen(&B);

    gf2_init(&Q, 1);
    gf2_init(&R, 1);

    printf("A = ");                 gf2_print(&A);
    printf("B = ");                 gf2_print(&B);

    gf2_long_division(&Q, &R, &A, &B);

    printf("Q = ");                 gf2_print(&Q);
    printf("R = ");                 gf2_print(&R);

}

void test_gf2_square()
{
    gf2 A, B;
    int t = 20;

    gf2_init(&A, t);
    gf2_init(&B, 1);

    gf2_random_gen(&A);
    gf2_square(&B, &A);

    printf("A = ");                 gf2_print(&A);
    printf("A^2 = ");               gf2_print(&B);

}

void test_gf2_pow()
{
    gf2 A, B;
    int t = 13;
    int e = 6;

    gf2_init(&A, t);
    gf2_init(&B, 1);

    gf2_random_gen(&A);
    gf2_pow(&B, &A, e);

    printf("A = ");                 gf2_print(&A);
    printf("A^%d = ", e);               gf2_print(&B);

}

void test_gf2_mulmod()
{
    gf2 A, B, C, mod;
    int m = 13;

    gf2_init(&mod, m);
    gf2_random_gen_fix(&mod);

    gf2_init(&A, 5);
    gf2_random_gen_fix(&A);

    gf2_init(&B, 9);
    gf2_random_gen_fix(&B);

    gf2_init(&C, m);

    printf("mod = ");               gf2_print(&mod);
    printf("A = ");                 gf2_print(&A);
    printf("B = ");                 gf2_print(&B);

    gf2_mulmod(&C, &A, &B, &mod);

    printf("C = ");                 gf2_print(&C);

}


void test_gf2_powmod()
{
    gf2 A, C, mod;
    int e = 10;
    int m = 13;

    gf2_init(&mod, m);
    gf2_random_gen_fix(&mod);

    gf2_init(&A, 5);
    gf2_random_gen_fix(&A);

    gf2_init(&C, m);

    printf("mod = ");               gf2_print(&mod);
    printf("A = ");                 gf2_print(&A);

    gf2_powmod(&C, &A, e, &mod);

    printf("A^%d = ", e);                 gf2_print(&C);

}

void test_gf2_gcd()
{
    gf2 A, B, gcd;
    int m = 13;

    gf2_init(&A, m);
    gf2_random_gen(&A);

    gf2_init(&B, m);
    gf2_random_gen(&B);

    gf2_gcd(&gcd, &A, &B);

    printf("A = ");                 gf2_print(&A);
    printf("B = ");                 gf2_print(&B);
    printf("gcd = ");               gf2_print(&gcd);

}

void test_gf2_square_root()
{
    gf2 A, B, mod;
    int m = 13;

/* z^13 + z^12 + z^10 + z^9 + z^7 + z^4 + 1, irreducible */
    gf2_init(&mod, m);
    gf2_set_index(&mod, 13);
    gf2_set_index(&mod, 12);
    gf2_set_index(&mod, 10);
    gf2_set_index(&mod, 9);
    gf2_set_index(&mod, 7);
    gf2_set_index(&mod, 4);
    gf2_set_index(&mod, 0);

    gf2_init(&A, m-1);
    gf2_random_gen(&A);
    gf2_square_root(&B, &A, &mod);

    printf("mod = ");               gf2_print(&mod);
    printf("A = ");                 gf2_print(&A);
    printf("square root of A = ");  gf2_print(&B);

}

void test_gf2_math_operation()
{
    /*
    printf("Start add test \n");
    test_gf2_add();

    printf("Start mul test \n");
    test_gf2_mul();

    printf("Start division test \n");
    test_gf2_div();
    
    printf("Start square test \n");
    test_gf2_square();
    
    printf("Start pow test \n");
    test_gf2_pow();
    
    printf("Start mulmod test \n");
    test_gf2_mulmod();
    
    printf("Start powmod test \n");
    test_gf2_powmod();
    
    printf("Start gcd test \n");
    test_gf2_gcd();
    */
    printf("Start square root test \n");
    test_gf2_square_root();
}


int main(){

    srand((unsigned int)time(NULL));
    //test_gf2_init();
    test_gf2_math_operation();
  
    return 1;
}
