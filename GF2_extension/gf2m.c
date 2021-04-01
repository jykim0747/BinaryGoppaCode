#include "gf2m.h"
#include "gf2_matrix.h"


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
        if((tmp>>j) & 0x01)
            printf("z^%d +", 8*q+j);
    }

    for(i=q-1; i>=0; i--)
    {   
        tmp = a->binary[i];
        for(j=7; j>=0; j--){   
            if((tmp>>j) & 0x01){
                printf("z^%d +", i*8+j);
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
    for(i = src->deg; i >= 0 ; i--)
    {
        gf2_init(&src->term[i], m);
        gf2_random_gen(&src->term[i]);
    }
    //gf2_init(&src->term[src->deg], m);
    ///gf2_random_gen_fix(&src->term[src->deg]);
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
    src->deg = 0;   /* 1의 차수는 상수항 0으로 간주한다. */
}
/////////////////////////////////////////////////////////////////////
int gf2m_is_one(gf2m* src)
{
    if((src->deg == 0) && (gf2_is_one(&src->term[0]) == ONE))
        return ONE;
    return NOT_ONE;
}
/////////////////////////////////////////////////////////////////////
void gf2m_copy(gf2m* dst, gf2m* src)
{
    int i;
    
    gf2m_set_zero(dst);
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
@   dst : dst = a*b (mod)
*/
void gf2m_mulmod(gf2m* dst, gf2m* a, gf2m* b, gf2m* gf2m_mod, gf2* mod)
{
    gf2m Q, tmp;
    gf2m_init(&Q, 1);
    gf2m_init(&tmp, 1);
    gf2m_mul(&tmp, a, b, mod);
    gf2m_long_division(&Q, dst, &tmp, gf2m_mod, mod);

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
    gf2m_fit_len(dst);

    return inv;
}
/////////////////////////////////////////////////////////////////////
/*
*   dst : src 모든 항에 element를 곱하는 함수
*/
void gf2m_mul_gf2_element(gf2m* dst, gf2m* src, gf2* element, gf2* mod)
{
    int i;
    dst->deg = src->deg;

    gf2_fit_len(element);
    for(i = 0; i<= src->deg; ++i)
    {
        gf2_mulmod(&(dst->term[i]), &(src->term[i]), element, mod);
    }
    gf2m_fit_len(dst);
}
/////////////////////////////////////////////////////////////////////
/* 
*   dst = src << shift
*/
void gf2m_shift(gf2m* dst, gf2m* src, int shift)
{
    int i;
    if(shift == 0)
    {
        gf2m_copy(dst, src);
        return;
    }

    dst->deg = src->deg + shift;

    for(i = src->deg; i >=0 ; i--)
    {
        gf2_copy(&(dst->term[i+shift]), &(src->term[i]));
    }
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
int gf2m_long_division(gf2m* Q, gf2m* R, gf2m* A, gf2m* B, gf2* mod)
{
    gf2m Q_tmp;
    gf2m R_tmp;
    gf2m add_tmp;
    gf2m B_monic;
    gf2m tmp;
    gf2 inv;
    int i;
    int len;

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
    inv = gf2m_monic(&B_monic, B, mod);

    gf2m_fit_len(R);
    for(i=R->deg; i>=0; i--)//차수 이동후 덧셈, 곱셈 진행
    {
        len = R->deg - B->deg;
        if(len >= 0)
        {
            gf2m_mul_gf2_element(&tmp, &B_monic, &(R->term[i]), mod);
            gf2m_set_index(&Q_tmp, &(R->term[i]), len);
            gf2m_shift(&R_tmp, &tmp, len);
            gf2m_add(&add_tmp, R, &R_tmp);
            gf2m_copy(R, &add_tmp);
            
            gf2m_set_zero(&R_tmp);
            gf2m_set_zero(&tmp);
            gf2m_fit_len(R);
        }
    }
    //monic으로 진행 했으므로 계산 후에 B와 Q를 변경시킨다.
    gf2m_mul_gf2_element(Q, &Q_tmp, &inv, mod);
    
    gf2m_fit_len(Q);
    gf2m_fit_len(R);

    return SUCCESS;
    
}
/////////////////////////////////////////////////////////////////////
/*
@   sage 다항식의 최대공약수. (z+1) 상수항을 1로 반환
@   gcd : a, b의 최대공약수
*/
void gf2m_gcd(gf2m* gcd, gf2m* a, gf2m* b, gf2* mod)
{
    gf2m Q, R;
    gf2m a_monic, b_monic;
    gf2m t0, t1, t2;

    gf2m_init(&R, 1);
    gf2m_init(&Q, 1);
    gf2m_init(&a_monic, a->deg);
    gf2m_init(&b_monic, b->deg);
    gf2m_init(&t0, a->deg);
    gf2m_init(&t1, b->deg);
    gf2m_init(&t2, a->deg);

    gf2m_monic(&a_monic, a, mod);
    gf2m_monic(&b_monic, b, mod);

    gf2m_copy(&t0, &a_monic);
    gf2m_copy(&t1, &b_monic);

    while(gf2m_is_zero(&t1) == NOT_ZERO)
    {
        gf2m_copy(&t2, &t0);
        gf2m_copy(&t0, &t1);
    
        gf2m_long_division(&Q, &R, &t2, &t1, mod);
        gf2m_copy(&t1, &R);
    }
    gf2m_fit_len(&t0);
    gf2m_copy(gcd, &t0);
}
/////////////////////////////////////////////////////////////////////
/*
@   gcd : a*a_inv + b*b_inv
@   (a,b)입력. b는 goppa poly(기약), (a_inv,b_inv)출력.
*/
void gf2m_xgcd(gf2m* gcd, gf2m* a, gf2m* b, gf2m* a_inv, gf2m* b_inv, gf2* mod)
{
    gf2m t0, t1, t2;
    gf2m v0, v1, v2;
    gf2m u0, u1, u2;
    gf2m R, Q;
    gf2m tmp;
    gf2 y, gcd_inv;
    gf2 gcd_tmp;

    int len = (a->deg <= b->deg) ? 2*b->deg : 2*a->deg;

    gf2m_init(&t0, len); gf2m_init(&t1, len); gf2m_init(&t2, len);
    gf2m_init(&v0, len); gf2m_init(&v1, len); gf2m_init(&v2, len);
    gf2m_init(&u0, len); gf2m_init(&u1, len); gf2m_init(&u2, len);
    gf2m_init(&tmp, len);

    gf2m_init(&Q, a->deg);
    gf2m_init(&R, a->deg);

    gf2m_copy(&t0, a);         //t0 = a
    gf2m_copy(&t1, b);         //t1 = b

    gf2m_set_one(&u0);         //u0 = 1
    gf2m_set_zero(&u1);        //u1 = 0

    gf2m_set_one(&v1);         //v1 = 1
    gf2m_set_zero(&v0);        //v0 = 0

    while(gf2m_is_zero(&t1) == NOT_ZERO)
    {
        gf2m_copy(&t2, &t0);             //t2 = t0
        gf2m_copy(&t0, &t1);             //t0 = t1

        gf2m_long_division(&Q, &R, &t2, &t1, mod); // R = t2%t1
        gf2m_copy(&t1, &R);              //t1 = R

        gf2m_set_zero(&Q);
        gf2m_set_zero(&R);

        gf2m_add(&tmp, &t2, &t1);        //t2 = t2 - t1
        
        gf2m_copy(&t2, &tmp);
        gf2m_long_division(&Q, &R, &t2, &t0, mod); // Q = t2/t0

        gf2m_copy(&u2, &u0);             //u2 = u0
        gf2m_copy(&v2, &v0);             //v2 = v0
        gf2m_copy(&u0, &u1);             //u0 = u1
        gf2m_copy(&v0, &v1);             //v0 = v1
        
        gf2m_set_zero(&tmp);
        gf2m_mul(&tmp, &Q, &u1, mod);    //u1 = u2 - qu1
        gf2m_add(&u1, &u2, &tmp);
        gf2m_set_zero(&tmp);
        gf2m_mul(&tmp, &Q, &v1, mod);    //v1 = v2 - qv1
        gf2m_add(&v1, &v2, &tmp);

        gf2m_set_zero(&Q);
        gf2m_set_zero(&tmp);
    }

    gf2m_fit_len(&u0);
    gf2m_fit_len(&v0);

    gf2m_copy(gcd, &t0);         //gcd = t0

    gf2_init(&gcd_inv, mod->deg-1); //gcd의 역원을 곱하여 반환
    gf2_xgcd(&gcd_tmp, &gcd_inv, &y, &gcd->term[0], mod);

    gf2m_mul_gf2_element(a_inv, &u0, &gcd_inv, mod);
    gf2m_mul_gf2_element(b_inv, &v0, &gcd_inv, mod);

    gf2m_fit_len(a_inv);    
    
}
/////////////////////////////////////////////////////////////////////
/*
@   dst = src^2
*/
void gf2m_square(gf2m* dst, gf2m* src, gf2* mod)
{
    int i;
    gf2 tmp;

    if( gf2m_is_zero(src) == ZERO)
    {
        gf2m_set_zero(dst);
        return ;
    } 

    if( gf2m_is_one(src) == ONE)
    {
        gf2m_set_one(dst);
        return ;
    }

    gf2_init(&tmp, 1);
    dst->deg = src->deg*2;

    for(i=src->deg; i>=0; i--)
    {
        if(gf2_is_zero(&src->term[i]) != ZERO)
        {
            gf2_squaremod(&tmp, &(src->term[i]), mod);
            gf2m_set_index(dst, &tmp, i*2);
        }
    }

    gf2m_fit_len(dst);
}
/////////////////////////////////////////////////////////////////////
/*
@   dst = src^2 (mod)
*/
void gf2m_squaremod(gf2m* dst, gf2m* src, gf2m* gf2m_mod, gf2* mod)
{
    gf2m Q,R;
    gf2m tmp;

    gf2m_init(&Q, 1);
    gf2m_init(&R, 1);
    gf2m_init(&tmp, 1);

    gf2m_square(&tmp, src, mod);
    gf2m_long_division(&Q, &R, &tmp, gf2m_mod, mod);
    gf2m_copy(dst, &R);

}
/////////////////////////////////////////////////////////////////////
void gf2m_left_to_right_mod(gf2m* dst, gf2m* a, gf2m* gf2m_mod, gf2* mod, int e)
{
    gf2m pow_tmp, pow_tmp2;
    gf2m Q, R;

    int tmp;
    int i,len;
    int e_tmp = e;

    gf2m_init(&pow_tmp, 1);
    gf2m_init(&pow_tmp2, 1);
    gf2m_init(&Q, a->deg);
    gf2m_init(&R, a->deg);

    tmp = e_tmp;
    binary_len(e_tmp, len);    //비트 길이 반환

    gf2m_set_one(&pow_tmp);
    for(i=len-1; i>=0; i--)
    {
        gf2m_set_zero(&pow_tmp2);
        gf2m_square(&pow_tmp2, &pow_tmp, mod);
        gf2m_fit_len(&pow_tmp2);
        gf2m_long_division(&Q, &R, &pow_tmp2, gf2m_mod, mod);
        gf2m_copy(&pow_tmp2, &R);
        if(((tmp>>i) & 0x1) == 1) //미확인
        {
            gf2m_set_zero(&pow_tmp);
            gf2m_mulmod(&pow_tmp, &pow_tmp2, a, gf2m_mod, mod);
        }
        else
        {
            gf2m_copy(&pow_tmp, &pow_tmp2);
        }
    }
    gf2m_copy(dst, &pow_tmp);
    gf2m_fit_len(dst);
}
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = arc^{2^e} mod (mod)
*/
void gf2m_repeated_squaremod(gf2m* dst, gf2m* src, int e, gf2m* gf2m_mod, gf2* mod)
{
    int i = 0;
    gf2m tmp, tmp2;
 
    if(e == 0)
    {
        gf2m_copy(dst, src);
        return;
    }

    gf2m_init(&tmp, 1);
    gf2m_init(&tmp2, 1);

    gf2m_copy(&tmp, src);
    while(i < e)
    {
        gf2m_squaremod(&tmp2, &tmp, gf2m_mod, mod);
        gf2m_copy(&tmp, &tmp2);
        i++;
    }
    gf2m_copy(dst, &tmp);
}
/////////////////////////////////////////////////////////////////////
void gf2m_powmod(gf2m* dst, gf2m* a, gf2m* gf2m_mod, gf2* mod, int e)
{
    gf2m Q, R;
    gf2m tmp;

    gf2m_init(&Q, a->deg);
    gf2m_init(&R, a->deg);
    gf2m_init(&tmp, a->deg);

    if(e == 0)
    {
        gf2m_set_one(dst);
    }
    if(e == 1)
    {
        gf2m_long_division(&Q, &R, a, gf2m_mod, mod);
        gf2m_copy(dst, &R);
    }
    else if( e == 2)
    {
        gf2m_square(&tmp, a, mod);
        gf2m_long_division(&Q, &R, &tmp, gf2m_mod, mod);
        gf2m_copy(dst, &R);
    }
    else
    {
        gf2m_left_to_right_mod(dst, a, gf2m_mod, mod, e);
    }
    
    gf2m_fit_len(dst);

}
/////////////////////////////////////////////////////////////////////
/*
@   gf2m 기약다항식 생성하는 함수
@   dst : 기약다항식
@     t : t차 ( gf2m )
*/
int gf2m_generate_irreducible(gf2m* dst, gf2* mod, int t)
{

    gf2 irre;   //t차 기약 다항식
    gf2m gf2m_mod;
    gf2m gf2m_tmp, gf2m_tmp2, gf2m_tmp3;
    gf2_MAT mat, ech_mat, mat_tmp;

    int i, j;
    int m = mod->deg - 1;
    int count = 0;
    int row = t;
    int col = t+1;

    gf2_init(&irre, t);
    gf2_generate_irreducible(&irre, t);
    
    gf2_matrix_init(&mat, row, col, m);
    gf2_matrix_init(&mat_tmp, row, col, m);
    gf2_matrix_init(&ech_mat, row, col, m);

    printf("irre ! "); gf2_print(&irre);
    
    gf2m_init(&gf2m_mod, t);
    gf2m_init(&gf2m_tmp, t - 1);
    gf2m_init(&gf2m_tmp2, t - 1);
    gf2m_init(&gf2m_tmp3, t - 1);

    /* gf2 to gf2m */
    {
        int qq = irre.deg / 8;
        int qr = irre.deg % 8;
        for(j = qr; j>=0; j--)
        {
            if((irre.binary[qq]>>j) & 0x01)
                gf2_set_one(&gf2m_mod.term[qq*8 + j]);
        }
        for(i = qq - 1; i>=0; i--)
            for(j = 7; j>=0; j--)
                if((irre.binary[i]>>j) & 0x01)
                    gf2_set_one(&gf2m_mod.term[i*8 + j]);
    }
    gf2m_fit_len(&gf2m_mod);

    do{
        gf2m_random_gen(&gf2m_tmp, m - 2);    //식 생성이 잘못됨
        gf2m_print(&gf2m_tmp);
        gf2m_set_one(&gf2m_tmp2);
        for(i=0; i<col; i++)
        {
            for(j=0; j<row; j++)
            {
                *gf2_mat_entry(mat, j, i) = (gf2m_tmp2.term[j]);
            }
            gf2m_mulmod(&gf2m_tmp3, &gf2m_tmp2, &gf2m_tmp, &gf2m_mod, mod);
            gf2m_copy(&gf2m_tmp2, &gf2m_tmp3);
        }
        gf2_matrix_echelon(ech_mat, mat, mod);

        count = 0;
        for(i=0; i<row; ++i)
        {
            if(gf2_is_one(gf2_mat_entry(mat, i, i)) == ONE)
                count ++;
        }
    }while(count != row);

    for(i=0; i<t; i++)
    {
        gf2_copy(&dst->term[i], gf2_mat_entry(ech_mat, i, col-1));
    }
    gf2_set_one(&dst->term[t]);
    dst->deg = t;
        
    gf2_matrix_free(ech_mat);
    gf2_matrix_free(mat);
    gf2_matrix_free(mat_tmp);

    return SUCCESS;
}
/////////////////////////////////////////////////////////////////////
/*
@   polynomial 제곱근 연산 방법. 사전 계산 X
@   dst : 제곱근 결과값
@     a : 입력값
@   gf2m_mod : modular
@   index : 지수 index
*/
void gf2m_square_root(gf2m* dst, gf2m* a, gf2m* gf2m_mod, gf2* mod)
{
    int i, index;
    gf2m Xi, Ri;
    gf2m Sq, Sq_root;
    gf2m Sq_odd, Sq_even;
    gf2m Sq_tmp, Sq_tmp2;

    gf2m X;
    gf2 gf2_X;
    gf2 gf2_tmp;

    gf2_init(&gf2_tmp, mod->deg-1);
    gf2_init(&gf2_X, mod->deg);
    gf2m_init(&X, 1);
    gf2_set_index(&gf2_X, 0);         //fq_X <- 'X'
    gf2m_set_index(&X, &gf2_X, 1);    //X <-fq_X
    gf2m_init(&Sq, 1);
    gf2m_init(&Sq_root, 1);
    gf2m_init(&Sq_tmp, 1);
    gf2m_init(&Sq_tmp, 2);
    gf2m_init(&Sq_odd, 1);
    gf2m_init(&Sq_even, 1);
    gf2m_init(&Xi, 1);

    gf2m_set_index(&Sq, &gf2_X, 1);  //(gf2m)X 생성.
    
    //Sq <- X^{2^{mt-1}} mod (gf2m_mod)
    index = (mod->deg  * gf2m_mod->deg) - 1;
    gf2m_repeated_squaremod(&Sq, &X, index, gf2m_mod, mod);
    
    for(i=0; i<=a->deg; i++)
    {
        gf2_powmod(&gf2_tmp, &(a->term[i]), (1<<(mod->deg-1)), mod);
        if((i % 2) == 0)
        {
            gf2m_set_index(&Sq_even, &gf2_tmp, i/2);
            gf2m_fit_len(&Sq_even);
        }
        else //odd
        {
            //Rx <- X^i*Sq mod (gf2m_mod)
            gf2m_powmod(&Xi, &X, gf2m_mod, mod, i/2);
            gf2m_mulmod(&Ri, &Sq, &Xi, gf2m_mod, mod);
            gf2m_mul_gf2_element(&Sq_odd, &Ri, &gf2_tmp, mod);
        }
        //누적 덧셈
        gf2m_add(&Sq_tmp, &Sq_even, &Sq_odd);
        gf2m_add(&Sq_root, &Sq_tmp2, &Sq_tmp);
        gf2m_copy(&Sq_tmp2, &Sq_root);

        gf2m_set_zero(&Sq_even);
        gf2m_set_zero(&Sq_odd);
        gf2_set_zero(&gf2_tmp);
    }
    gf2m_fit_len(&Sq_root);
    gf2m_copy(dst, &Sq_root);
}


/*
@   gf2m 미분
@   dst : 미분값
@   src : 입력값
*/
void gf2m_diff(gf2m* dst, gf2m* src)
{
    int i;

    for(i=1; i<src->deg; i+=2)
        gf2_copy(&dst->term[i-1], &src->term[i]);

    gf2m_fit_len(dst);
}
