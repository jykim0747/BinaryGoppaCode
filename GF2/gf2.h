#ifndef _GF2_H_
#define _GF2_H_

#define size 100

/* GF2 element struct */
typedef struct 
{
    int deg;
    unsigned char binary[size];
}gf2;

#define binary_len(X, Y) {\
    for((Y)=0; (X)>0; (Y)++) (X) /= 2;\
}

#define ZERO        0000
#define ONE         0001
#define NOT_ZERO    0002
#define NOT_ONE     0003

#define FAILURE     0020
#define SUCCESS     0021

void gf2_init(gf2* a, int t);
void gf2_print(gf2* a);
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

void gf2_left_to_right(gf2* dst, gf2* a, int e);
void gf2_pow(gf2* dst, gf2* a, int e);
void gf2_powmod(gf2* dst, gf2* a, int e, gf2* mod);

int fq_gcd(gf2* gcd, gf2* a, gf2* b);

#endif