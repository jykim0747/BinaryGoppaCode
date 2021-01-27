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

void test_gf2m_math_operation();

#endif