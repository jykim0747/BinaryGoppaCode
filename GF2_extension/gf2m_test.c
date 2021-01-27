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

void test_gf2m_math_operation(){

    test_gf2m_init();
}
