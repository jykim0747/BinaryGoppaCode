#ifndef _GF2M_H_
#define _GF2M_H_

#include "gf2.h"

#define MAX_DEGREE 150

/* gf2 extension element struct */
typedef struct 
{
    int deg;                // 차수
    gf2 term[MAX_DEGREE];   // 항의 개수, 차수 범위 문제로 넉넉히 설정.
}gf2m;


void gf2m_init(gf2m* a, int t);
void gf2m_print(gf2m* src);

void gf2m_fit_len(gf2m* src);
void gf2m_random_gen(gf2m* src, int m);

void gf2m_set_zero(gf2m* src);
void gf2m_set_one(gf2m* src);
void gf2m_copy(gf2m* dst, gf2m* src);
void gf2m_set_index(gf2m* dst, gf2* src, int inedx);

int gf2m_is_one(gf2m* src);
int gf2m_is_zero(gf2m* src);

void gf2m_add(gf2m* dst, gf2m* src1, gf2m* src2);
void gf2m_mul_shcool(gf2m* dst, gf2m* a, gf2m* b, gf2* mod);
void gf2m_mul(gf2m* dst, gf2m* a, gf2m* b, gf2* mod);

gf2 gf2m_monic(gf2m* dst, gf2m* a, gf2* mod);
void gf2m_mul_gf2_element(gf2m* dst, gf2m* src, gf2* element, gf2* mod);
void gf2m_shift(gf2m* dst, gf2m* src, int shift);
int gf2m_long_division(gf2m* Q, gf2m* R, gf2m* A, gf2m* B, gf2* mod);
void gf2m_gcd(gf2m* gcd, gf2m* a, gf2m* b, gf2* mod);
void gf2m_xgcd(gf2m* gcd, gf2m* a, gf2m* b, gf2m* a_inv, gf2m* b_inv, gf2m* mod);

void test_gf2m_math_operation();

#endif