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
        if( gf2_is_zero(&src->term[i]) != ZERO)
        {
            gf2_fit_len(&src->term[i]);
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
    gf2_init(&src->term[src->deg], m);
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
/////////////////////////////////////////////////////////////////////
/*
@   dst : a의 monic 다항식
@     a : 입력값
@   mod : modular (이진유한체)
@   return : a의 최고차항의 역원
*/
gf2 gf2m_monic(gf2m* dst, gf2m* a, gf2* mod)
{
    gf2 gcd;
    gf2 inv;
    gf2 tmp;
    gf2 a_tmp;
    int i;

    gf2_init(&gcd, 1);
    gf2_init(&inv, 1);
    gf2_init(&tmp, 1);
    gf2_init(&a_tmp, 1);
    
    gf2_copy(&a_tmp, &(a->term[a->deg]));
    gf2_fit_len(&a_tmp);

    gf2_xgcd(&gcd, &inv, &tmp, &a_tmp, mod);

    for(i=a->deg; i>=0; i--)
    {
        gf2_mulmod(&(dst->term[i]), &inv, &(a->term[i]), mod);
    }

    return inv;
}
/////////////////////////////////////////////////////////////////////
/*
@   Q : Quotient
@   R : Remainder
@   A : Numerator
@   B : Denominator
@   A = QB + R
@   mod : mod gf2
*/
int gf2m_long_division(gf2m* Q, gf2m* R, gf2m* A, gf2m* B, gf2m* mod)
{
    gf2m Q_tmp;
    gf2m R_tmp;
    gf2m add_tmp;
    gf2m B_monic;
    gf2m tmp;
    gf2 inv;

    if(gf2m_is_zero(B) == ZERO)
    {
        printf("0으로 나눌 수 없습니다.\n");
        return FAILURE;
    }

    if(A->deg < B->deg)
    {
        gf2m_set_zero(Q);
        gf2m_copy(R, A);
        return SUCCESS;
    }
    if( (A->deg == 0) && (B->deg == 0))
    {
        gf2_long_division(&(Q->term[0]), &(R->term[0]), &(A->term[0]), &(B->term[0]));
        gf2m_fit_len(Q);
        gf2m_fit_len(R);
        return SUCCESS;
    }

    gf2m_init(&Q_tmp, A->deg);
    gf2m_init(&R_tmp, 1);
    gf2m_init(&add_tmp, 1);
    gf2m_init(&B_monic, 1);
    gf2m_init(&tmp, 1);
    
    gf2m_copy(&B_monic, B);
    gf2m_copy(R, A);
    
    //필요
    //inv = fqt_make_monic(&B_monic, B, mod);
    /*

    //printf("inv\n");
    //fq_print(&inv);
    //printf("monic\n");
    //fqt_print(&B_monic);

    fqt_fit_len(R);
    //printf("Rdeg %d, Bdeg %d \n", R->deg, B->deg);
    for(i=R->deg; i>=0; i--)//차수 이동후 덧셈, 곱셈 진행
    {
        len = R->deg - B->deg;
        if(len >= 0)
        {
            fqt_mul_coe(&tmp, &B_monic, &(R->term[i]), mod);
            fqt_set_index(&Q_tmp, &(R->term[i]), len);
            fqt_shift(&R_tmp, &tmp, len);
            fqt_add(&add_tmp, R, &R_tmp);
            fqt_copy(R, &add_tmp);
            
            fqt_set_zero(&R_tmp);
            fqt_set_zero(&tmp);
            fqt_fit_len(R);

        }
    }
    //monic으로 진행 했으므로 계산 후에 B와 Q를 변경시킨다.
    fqt_mul_coe(Q, &Q_tmp, &inv, mod);
    
    fqt_fit_len(Q);
    fqt_fit_len(R);

    fqt_clear(&Q_tmp);
    fqt_clear(&R_tmp);
    fqt_clear(&add_tmp);
    fqt_clear(&B_monic);
    fqt_clear(&tmp);
*/
    return SUCCESS;
    
}