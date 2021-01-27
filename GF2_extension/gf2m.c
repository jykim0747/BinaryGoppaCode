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
            printf(")*x^%d", i);
        }
        else
        {
            printf(")*x^%d+", i);
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