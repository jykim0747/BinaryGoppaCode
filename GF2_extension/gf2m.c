#include "gf2m.h"


/////////////////////////////////////////////////////////////////////
static void gf2_print_noline(gf2* a)
{
    int q  = a->deg / 8;
    int qr = a->deg % 8;

    int i = q;
    int j;
    int tmp;

    if(gf2_is_zero(a) == ZERO){
        printf("0");
        return;
    }

    tmp = a->binary[q];
    for(j=qr; j>=0; j--){
        if((tmp>>j)&0x01)
            printf("z^%d+", 8*q+j);
    }

    for(i=q-1; i>=0; i--)
    {   
        tmp = a->binary[i];
        for(j=7; j>=0; j--){   
            if((tmp>>j)&0x01){
                printf("z^%d+", i*8+j);
            }
        }
    }
    printf("\b");
}
/////////////////////////////////////////////////////////////////////
void gf2m_print(gf2m* src)
{
    int i;
    for(i= src->deg; i>=0; i--)
    {
        printf("(");
        gf2_print_noline(&(src->term[i]));
        if(i == 0)
        {
            printf(")*X^%d", i);
        }
        else
        {
            printf(")*X^%d+", i);
        }
    }
    printf("\n");
}
/////////////////////////////////////////////////////////////////////
void gf2m_init(gf2m* a, int t)
{
    int i;
    for(i=0; i<MAX_DEGREE; ++i)
        gf2_set_zero(&a->term[i]);

    a->deg = t;
}

/////////////////////////////////////////////////////////////////////
void gf2m_fit_len(gf2m* src)
{
    int i;

    for(i=MAX_DEGREE - 1; i>=0; i--)
    {
        gf2_fit_len(&src->term[i]);
        if( gf2_is_zero(&src->term[i]) != ZERO)
        {
            src->deg = i;
            return ;
        }
    }
    src->deg = 0;
}
/////////////////////////////////////////////////////////////////////
/*
@   Randomly generating 
@   degree of GF2 : m
@   degree of GF2m : t
*/
void gf2m_random_gen(gf2m* src, int m)
{
    int i;
    for(i = src->deg - 1; i >= 0 ; i--)
    {
        gf2_init(&src->term[i], m);
        gf2_random_gen(&src->term[i]);
    }

    gf2_random_gen_fix(&src->term[src->deg]);
    gf2m_fit_len(src);

}
/////////////////////////////////////////////////////////////////////
void gf2m_set_zero(gf2m* src)
{
    int i;
    for(i=0; i<MAX_DEGREE; ++i)
        gf2_set_zero(&src->term[i]);
    src->deg = 0;   /* 0의 차수는 0으로 간주한다. */
}
/////////////////////////////////////////////////////////////////////
int gf2m_is_zero(gf2m* src)
{
    if((src->deg == 0) && (gf2_is_zero(&src->term[0]) == ZERO))
        return ZERO;
    return NOT_ZERO;
}
/////////////////////////////////////////////////////////////////////
void gf2m_set_one(gf2m* src)
{
    gf2m_set_zero(src);
    gf2_set_one(&src->term[0]);
    src->deg = 1;   /* 1의 차수는 상수항 1로 간주한다. */
}
/////////////////////////////////////////////////////////////////////
int gf2m_is_one(gf2m* src)
{
    if((src->deg == 1) && (gf2_is_one(&src->term[0]) == ONE))
        return ONE;
    return NOT_ONE;
}
/////////////////////////////////////////////////////////////////////
void gf2m_copy(gf2m* dst, gf2m* src)
{
    int i;
    
    dst->deg = src->deg;
    for(i = 0; i <=src->deg; ++i)
    {
        gf2_copy(&dst->term[i], &src->term[i]);
    }
}
/////////////////////////////////////////////////////////////////////
/*
@   dst->term[index] = src
*/
void gf2m_set_index(gf2m* dst, gf2* src, int index)
{
    gf2_copy(&dst->term[index], src);
}
/////////////////////////////////////////////////////////////////////
/*
@   dst = src1 + src2  
*/
void gf2m_add(gf2m* dst, gf2m* src1, gf2m* src2)
{
    int i;
    dst->deg = (src1->deg >= src2->deg) ? src1->deg : src2->deg;

    for(i = dst->deg; i>=0; i--)
    {
        gf2_add(&dst->term[i], &src1->term[i], &src2->term[i]);
    }
    gf2m_fit_len(dst);
}
/////////////////////////////////////////////////////////////////////
void gf2m_mul_shcool(gf2m* dst, gf2m* a, gf2m* b, gf2* mod)
{
    int i, j;
    gf2 gf2_tmp;
    gf2m mul_tmp, add_tmp, gf2m_tmp;
    
    gf2_init(&gf2_tmp, mod->deg);
    gf2m_init(&mul_tmp, a->deg + b->deg);
    gf2m_init(&add_tmp, 1);
    gf2m_init(&gf2m_tmp, 1);

    for(i = a->deg; i >= 0; i--)
    {
        for(j = b->deg; j >= 0; j--)
        {
            gf2_mulmod(&gf2_tmp, &(a->term[i]), &(b->term[j]), mod);
            gf2m_set_index(&mul_tmp, &gf2_tmp, i+j);
            gf2m_fit_len(&mul_tmp);
            gf2m_add(&gf2m_tmp, &add_tmp, &mul_tmp);
            gf2m_copy(&add_tmp, &gf2m_tmp);
            
            gf2m_set_zero(&mul_tmp);
            gf2_set_zero(&gf2_tmp);
        }
    }
    gf2m_copy(dst, &add_tmp);

}
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = a*b
*/
void gf2m_mul(gf2m* dst, gf2m* a, gf2m* b, gf2* mod)
{
    if(gf2m_is_one(a) == ONE)
    {
        gf2m_copy(dst, b);
        return;
    }
    if(gf2m_is_one(b) == ONE)
    {
        gf2m_copy(dst, a);
        return;
    }
    if((gf2m_is_zero(a) == ZERO) || (gf2m_is_zero(b) == ZERO))
    {
        gf2m_set_zero(dst);
        return;
    }
    gf2m_mul_shcool(dst, a, b, mod);
    gf2m_fit_len(dst);

}