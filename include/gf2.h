#ifndef _GF2_H_
#define _GF2_H_

#include <stdio.h>
#include <stdlib.h>
#include "error.h"

#define size 500

/* GF2 element struct */
typedef struct 
{
    int deg;
    unsigned char binary[size];
}gf2;

#define binary_len(X, Y) {\
    for((Y)=0; (X)>0; (Y)++) (X) /= 2;\
}




void gf2_init(gf2* a, int t);
void gf2_print(gf2* a);
void gf2_print_pretty(gf2* a);
void gf2_random_gen(gf2* a);
void gf2_random_gen_fix(gf2* a);
void gf2_fit_len(gf2* a);
int gf2_is_zero(gf2* a);
int gf2_is_one(gf2* a);
void gf2_set_zero(gf2* a);
void gf2_set_one(gf2* a);
void gf2_copy(gf2* dst, gf2* a);
void gf2_set_index(gf2* a, int index);

void gf2_add(gf2* dst, gf2* a, gf2* b);
void gf2_addmod(gf2* dst, gf2* a, gf2* b, gf2* mod);
void gf2_mul_shift(gf2* dst, gf2* a, int index);
void gf2_mul_shcool(gf2* dst, gf2* a, gf2* b);
void gf2_mul(gf2* dst, gf2* a, gf2* b);
int gf2_mulmod(gf2* dst, gf2* a, gf2* b, gf2* mod);

int gf2_long_division(gf2* Q, gf2* R, gf2* A, gf2* B);
void gf2_square(gf2* dst, gf2* a);
void gf2_squaremod(gf2* dst, gf2* a, gf2* mod);
void gf2_repeated_squaremod(gf2* dst, gf2* a, int e, gf2* mod);

void gf2_left_to_right(gf2* dst, gf2* a, int e);
void gf2_pow(gf2* dst, gf2* a, int e);
void gf2_powmod(gf2* dst, gf2* a, int e, gf2* mod);

int gf2_gcd(gf2* gcd, gf2* a, gf2* b);
void gf2_square_root(gf2* dst, gf2* a, gf2* mod);

int gf2_is_irreducible(gf2* src);
void gf2_generate_irreducible(gf2* src, int degree);

void gf2_xgcd(gf2* gcd, gf2* x, gf2* y, gf2* a, gf2* b);
void gf2_diff(gf2* dst, gf2* src);
int gf2_berlekamp_factoring(gf2* src);

void test_gf2_init();
void test_gf2_math_operation();


#endif