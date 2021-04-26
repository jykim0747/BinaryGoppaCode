#include "gf2.h"
#include <stdio.h>
#include <stdlib.h>
#include "bmatrix.h"
/*
@   a : Fq element
@   t : degree of a (with constant) (ex. z^2+z+1 -> 3) maximum degree is 14 (z^13)
*/
void gf2_init(gf2* a, int t)
{
    int iter;

    a->deg = t + 1;    
    for(iter = 0; iter < size; ++iter)
        a->binary[iter] = 0;
}
/////////////////////////////////////////////////////////////////////
void gf2_print(gf2* a)
{
    int q  = a->deg / 8;
    int qr = a->deg % 8;

    int i = q;
    int j;
    int tmp;

    if(gf2_is_zero(a) == ZERO)
    {
        printf(" 0\n");
        return;
    }

    tmp = a->binary[q];
    for(j=qr; j>=0; j--)
    {
        if((tmp>>j)&0x01)
            printf("z^%d+", 8*q+j);
    }

    for(i=q-1; i>=0; i--)
    {   
        tmp = a->binary[i];
        for(j=7; j>=0; j--)
        {   
            if((tmp>>j)&0x01)
            {
                printf("z^%d+", i*8+j);
            }
        }
    }
    printf("\b\n");

}
/////////////////////////////////////////////////////////////////////
void gf2_print_pretty(gf2* a)
{
    int q  = a->deg / 8;
    int qr = a->deg % 8;

    int i = q;
    int j;
    int tmp;

    if(gf2_is_zero(a) == ZERO)
    {
        printf("0");
        return;
    }

    tmp = a->binary[q];
    for(j=qr; j>=0; j--)
    {
        if((tmp>>j)&0x01)
            printf("z^%d+", 8*q+j);
    }

    for(i=q-1; i>=0; i--)
    {   
        tmp = a->binary[i];
        for(j=7; j>=0; j--)
        {   
            if((tmp>>j)&0x01)
            {
                printf("z^%d+", i*8+j);
            }
        }
    }
    printf("\b");

}
/////////////////////////////////////////////////////////////////////
/*
@   fitting length
*/
void gf2_fit_len(gf2* a)
{
    int i,j;
    int q = a->deg / 8;

    for(j=7; j>=0; j--)
    {
        if(((a->binary[q])>>j) & 0x01)
        {
            a->deg = q * 8 + j + 1;
            return;
        }
    }
    for(i=q-1; i>=0; i--)
    {
        for(j=7; j>=0; j--)
        {
            if(((a->binary[i])>>j) & 0x01)
            {
                a->deg = i * 8 + j + 1;
                return;
            }
        }
    }
    a->deg = 1;
}
/////////////////////////////////////////////////////////////////////
/*
@   check zero
*/
int gf2_is_zero(gf2* a)
{
    if((a->binary[0] == 0) && (a->deg == 1))
        return ZERO;
    
    return NOT_ZERO;
}
/////////////////////////////////////////////////////////////////////
/*
@   check one
*/
int gf2_is_one(gf2* a)
{
    if((a->binary[0] == 1) && (a->deg == 1))
        return ONE;
    return NOT_ONE;
}
/////////////////////////////////////////////////////////////////////
/*
@   set zero
*/
void gf2_set_zero(gf2* a)
{
    int i;
    a->deg = 1;
    for(i=0; i<size; i++)
        a->binary[i] = 0;
}
/////////////////////////////////////////////////////////////////////
/*
@   set one
*/
void gf2_set_one(gf2* a)
{
    int i;
    a->deg = 1;
    for(i=0; i<size; i++)
        a->binary[i] = 0;

    a->binary[0] = 1;
}
/////////////////////////////////////////////////////////////////////
/*
@   Randomly generating less than the maximum order of a
*/
void gf2_random_gen(gf2* a)
{
    int q  = a->deg / 8;
    int qr = a->deg % 8;
    int i;

    for(i=q-1; i>=0; i--)
    {
        a->binary[i] = rand();
    }
    a->binary[q] = rand() % (1<<qr);

    gf2_fit_len(a);
}
/////////////////////////////////////////////////////////////////////
/*
@   Randomly generating by order of a
*/
void gf2_random_gen_fix(gf2* a)
{
    int q  = (a->deg - 1) / 8;
    int qr = (a->deg - 1) % 8;
    int i;

    for(i=q-1; i>=0; i--)
    {
        a->binary[i] = rand();
    }
    a->binary[q] = rand() & ((1<<qr)-1);
    a->binary[q] ^= 1<<qr;
}
/////////////////////////////////////////////////////////////////////
/*
@   copy
*/
void gf2_copy(gf2* dst, gf2* a)
{
    int i;
    
    dst->deg = a->deg;

    for(i=(a->deg); i>=0; i--)
    {
        dst->binary[i] = a->binary[i];
    }

}
/////////////////////////////////////////////////////////////////////
/*
@   set the index-th position of a to 1
*/
void gf2_set_index(gf2* a, int index)
{
    int idxq = (index / 8);
    int idxr = (index % 8);

    a->binary[idxq] ^= (1 << idxr);
}
/////////////////////////////////////////////////////////////////////
/*
@   dst = a + b
*/
void gf2_add(gf2* dst, gf2* a, gf2* b)
{
    int i;
    int len = (a->deg <= b->deg) ? b->deg : a->deg;

    for(i=0; i<= (len/8); i++)
        dst->binary[i] = a->binary[i] ^ b->binary[i];

    dst->deg = len;
    gf2_fit_len(dst);
}
/////////////////////////////////////////////////////////////////////
/*
@   dst = a + b mod (mod)
*/
void gf2_addmod(gf2* dst, gf2* a, gf2* b, gf2* mod)
{
    gf2 Q, R;
    gf2 tmp;

    gf2_init(&Q, a->deg);
    gf2_init(&R, a->deg);
    gf2_init(&tmp, a->deg);

    gf2_add(&tmp, a, b);
    gf2_long_division(&Q, &R, &tmp, mod);
    gf2_copy(dst, &R);

    gf2_fit_len(dst);

}
/////////////////////////////////////////////////////////////////////
/*
@   shift operation over binary field
@   dst = a << index
*/
void gf2_mul_shift(gf2* dst, gf2* a, int index)
{
    int i,j;
    unsigned char tmp;

    int aq = a->deg / 8;
    int idxq;
    int idxr; 

    if(index == 0)
    {    
        gf2_copy(dst, a);
        return;
    }

    dst->deg = a->deg + index;

    for(i = aq; i>=0; i--)
    {
        for(j=7; j>=0; j--)
        {
            tmp = a->binary[i];
            if(((tmp>>j) & 0x01 )==1)
            {
                idxq = (i*8+j+index) / 8;
                idxr = (j+index) % 8;
                dst->binary[idxq] ^= 1<< idxr;
            }
        }
    }
}
/////////////////////////////////////////////////////////////////////
/*
@   dst = a * b
*/
void gf2_mul_shcool(gf2* dst, gf2* a, gf2* b)
{
    int i,j;
    gf2 tmp;
    int aq = a->deg / 8;

    gf2_init(&tmp, a->deg + b->deg -1);
    gf2_init(dst, a->deg + b->deg -1);

    for(i=aq; i>=0; i--)
    {
        for(j=7; j>=0; j--)
        {
            if((a->binary[i]>>j) & 0x01)
            {
                gf2_mul_shift(&tmp, b, i*8 + j);
                gf2_add(dst, dst, &tmp);
                gf2_set_zero(&tmp);
            }
        }
    }
}
/////////////////////////////////////////////////////////////////////
/*
@   dst = a * b
*/
void gf2_mul(gf2* dst, gf2* a, gf2* b)
{
    if(gf2_is_one(a) == ONE)
    {
        gf2_copy(dst, b);
        return ;
    }
    if(gf2_is_one(b) == ONE)
    {
        gf2_copy(dst, a);
        return ;
    }
    if((gf2_is_zero(a) == ZERO) || (gf2_is_zero(b) == ZERO))
    {
        gf2_set_zero(dst);
        return ;
    }

    gf2_mul_shcool(dst, a, b);
    gf2_fit_len(dst);

}
/////////////////////////////////////////////////////////////////////
/*
@   dst = a * b mod (mod)
*/
int gf2_mulmod(gf2* dst, gf2* a, gf2* b, gf2* mod)
{
    gf2 Q;
    gf2 tmp;

    gf2_init(&Q, 1);

    gf2_mul(&tmp, a, b);
    gf2_long_division(&Q, dst, &tmp, mod);

    /* debug */
    /*
    printf("Q = ");                 gf2_print(&Q);
    printf("mod = ");                 gf2_print(mod);
    printf("dst = ");                 gf2_print(dst);
    printf("tmp = ");                 gf2_print(&tmp);
    */

    return SUCCESS;
}

/*
@   Q : Quotient
@   R : Remainder
@   A : Numerator
@   B : Denominator
@   A = QB + R
*/
int gf2_long_division(gf2* Q, gf2* R, gf2* A, gf2* B)
{
    gf2_fit_len(A);
    gf2_fit_len(B);

    if( gf2_is_zero(B) == ZERO)
    {
        printf("Cannot be divided by Zero \n");
        return FAILURE;
    }
    if( gf2_is_one(B) == ONE)
    {
        gf2_copy(Q, A);
        gf2_set_zero(R);
        return SUCCESS;
    }
    if( A->deg < B->deg)
    {
        gf2_set_zero(Q);
        gf2_copy(R, A);
        return SUCCESS;
    }
    else if(A ->deg == B->deg)
    {
        gf2_set_one(Q);
        gf2_copy(R, A);
        gf2_add(R, R, B);
        return SUCCESS;
    }
    
    gf2 R_tmp;
    gf2 B_tmp;
    int len = A->deg - B->deg;

    gf2_init(&R_tmp, 1);
    gf2_copy(&R_tmp, A);
    gf2_init(&B_tmp, 1);
    
    
    while(len >= 0)//len이 워드단위를 여러개 가질경우.
    {
        gf2_mul_shift(&B_tmp, B, len);
        gf2_add(&R_tmp, &R_tmp, &B_tmp);
        gf2_set_index(Q, len);
        gf2_set_zero(&B_tmp);
                    
        len = R_tmp.deg - B->deg;
    }
    gf2_copy(R, &R_tmp);

    gf2_fit_len(Q);
    gf2_fit_len(R);

    return SUCCESS;
}
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = a^2
*/
void gf2_square(gf2* dst, gf2* a)
{
    int i, j;
    int aq = (a->deg / 8);
    int dq, dr;

    if( gf2_is_one(a) == ONE)
    {
        gf2_set_one(dst);
        return;
    }

    gf2_init(dst, 2*(a->deg));
    for(i=aq; i>=0; i--)
    {
        for(j=7; j>=0; j--)
        {
            if((a->binary[i]>>j) & (0x01))
            {
                dq = (i*8+j)*2 / 8;
                dr = (i*8+j)*2 % 8;
                dst->binary[dq] ^= 1<< dr;
            }
        }
    }
    gf2_fit_len(dst);
}
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = a^2 mod (mod)
*/
void gf2_squaremod(gf2* dst, gf2* a, gf2* mod)
{
    gf2 tmp, Q;

    gf2_init(&Q, 1);
    gf2_square(&tmp, a);
    gf2_long_division(&Q, dst, &tmp, mod);
    
}
/////////////////////////////////////////////////////////////////////
/*
@   method of powering
*/
void gf2_left_to_right(gf2* dst, gf2* a, int e)
{
    int i, len, e_tmp = e;
    int tmp;
    gf2 pow_tmp;
    gf2 pow_tmp2;

    gf2_init(&pow_tmp, 1);
    gf2_init(&pow_tmp2, 1);

    binary_len(e_tmp, len);    //len : length of e bits

    gf2_set_one(&pow_tmp);
    for(i=len; i>=0; i--)
    {
        gf2_set_zero(&pow_tmp2);
        gf2_square(&pow_tmp2, &pow_tmp);

        tmp = e;
        if((tmp>>i) & 0x01)
        {
            gf2_mul(&pow_tmp, &pow_tmp2, a);
        }
        else
        {
            gf2_copy(&pow_tmp, &pow_tmp2);
        }   
    }
    gf2_copy(dst, &pow_tmp);
    gf2_fit_len(dst);
    
}
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = a^e
*/
void gf2_pow(gf2* dst, gf2* a, int e)
{
    if( e == 0)
    {
        gf2_set_one(dst);
        return;
    }
    else if( e == 1)
    {
        gf2_copy(dst, a);
        return;
    }
    else if( e == 2)
    {
        gf2_square(dst, a);
        return;
    }
    else
    {
        gf2_left_to_right(dst, a, e);
        gf2_fit_len(dst);
        return;
    }
}
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = a^e mod (mod)
*/
void gf2_powmod(gf2* dst, gf2* a, int e, gf2* mod)
{
    gf2 Q, tmp;

    gf2_init(&Q, 1);
    gf2_pow(&tmp, a, e);
    gf2_long_division(&Q, dst, &tmp, mod);

}
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = a^{2^e} mod (mod)
*/
void gf2_repeated_squaremod(gf2* dst, gf2* a, int e, gf2* mod)
{
    int i = 0;
    gf2 tmp, tmp2;
 
    if(e == 0)
    {
        gf2_copy(dst, a);
        return;
    }

    gf2_init(&tmp, 1);
    gf2_init(&tmp2, 1);

    gf2_copy(&tmp, a);
    while(i < e)
    {
        gf2_squaremod(&tmp2, &tmp, mod);
        gf2_copy(&tmp, &tmp2);
        i++;
    }
    gf2_copy(dst, &tmp);
}
/////////////////////////////////////////////////////////////////////
/*                                                        
@   gcd = gcd(A,B)    
@   if A = 0, then gcd(A, B) == B 
*/  
int gf2_gcd(gf2* gcd, gf2* a, gf2* b)
{  
    gf2 R, Q;
    gf2 t0, t1, t2;
    int res;
                                                          
    gf2_init(&R, 1);       
    gf2_init(&Q, 1);
    gf2_init(&t0, 1);                                      
    gf2_init(&t1, 1);
    gf2_init(&t2, 1);
                                                          
    gf2_copy(&t0, a);
    gf2_copy(&t1, b);
                                                          
    while(gf2_is_zero(&t1) != ZERO)
    {                                                     
        gf2_copy(&t2, &t0);
        gf2_copy(&t0, &t1);
        res = gf2_long_division (&Q, &R, &t2, &t1); 
        gf2_copy(&t1, &R);      
    }                                                     
    gf2_copy(gcd, &t0);                                    
                                          
    return res;                                           
}                                                         
/////////////////////////////////////////////////////////////////////
/*
@   dst : dst = a^{2^{m-1}} mod (mod)
*/
void gf2_square_root(gf2* dst, gf2* a, gf2* mod)
{
    int index;      /* index : 2^{m-1}의 m-1을 의미한다. */
    index = mod->deg - 2;

    gf2_repeated_squaremod(dst, a, index, mod);

}
/////////////////////////////////////////////////////////////////////
/*
@   return IRREDUCIBLE or REDUCIBLE
@   Berlekamp algorithm 사용
*/
int gf2_is_irreducible(gf2* src)
{
    int factors;

    if(src->deg == 0)
    {
        return REDUCIBLE;
    }

    factors = gf2_berlekamp_factoring(src);
    if(factors != 1)
        return REDUCIBLE;

    return IRREDUCIBLE;

}
/////////////////////////////////////////////////////////////////////
/*
@   generate irreducible polynomial GF2
*/
void gf2_generate_irreducible(gf2* src, int degree)
{
    gf2 tmp;
    int res = 1;

    gf2_init(&tmp, degree);
    
    while(res != IRREDUCIBLE)
    {
        gf2_random_gen_fix(&tmp);
        //z^6 + z^5 + z^4 + z^3 + z^2 + z + 1 = (z^3 + z + 1) * (z^3 + z^2 + 1) irre로 나옴
        // //z^6, z^6+z^2+z^0
        //  gf2_set_index(&tmp, 6);
        //  gf2_set_index(&tmp, 2);
        //  gf2_set_index(&tmp, 0);
        // gf2_set_index(&tmp, 3);
        // gf2_set_index(&tmp, 1);
        // gf2_set_index(&tmp, 0);
        // gf2_fit_len(&tmp);

        res = gf2_is_irreducible(&tmp);
    }
    gf2_copy(src, &tmp);
}
/////////////////////////////////////////////////////////////////////
/*
@   gcd : gcd = ax+by, a, b의 최대공약수 (출력값)
@     x : a의 역원 (출력값)
@     y : b의 역원 (출력값)
@     a : 입력 대상1 (입력값)
@     b : 입력 대상2 (입력값)
*/
void gf2_xgcd(gf2* gcd, gf2* x, gf2* y, gf2* a, gf2* b)
{
    gf2 t0, t1, t2;
    gf2 v0, v1, v2;
    gf2 u0, u1, u2;
    gf2 R , Q;
    gf2 tmp;

    int len = a->deg <= b->deg ? 2*(b->deg-1) : 2*(a->deg-1);

    gf2_init(&t0, 1);    gf2_init(&t1, 1);    gf2_init(&t2, 1);
    gf2_init(&v0, 1);    gf2_init(&v1, 1);    gf2_init(&v2, 1);
    gf2_init(&u0, 1);    gf2_init(&u1, 1);    gf2_init(&u2, 1);
    gf2_init(&tmp, len);

    gf2_init(&R, 1);
    gf2_init(&Q, 1);

    gf2_copy(&t0, a);
    gf2_copy(&t1, b);

    gf2_set_one(&v1);
    gf2_set_zero(&v0);

    gf2_set_one(&u0);
    gf2_set_zero(&u1);

    while(gf2_is_zero(&t1) != ZERO)
    {
        gf2_copy(&t2, &t0);
        gf2_copy(&t0, &t1);

        gf2_set_zero(&Q);
        gf2_long_division(&Q, &R, &t2, &t1);
        gf2_copy(&t1, &R);
        gf2_add(&tmp, &t2, &t1);
        gf2_copy(&t2, &tmp);
        gf2_set_zero(&Q);
        gf2_long_division(&Q, &R, &t2, &t0);
         
        gf2_copy(&u2, &u0);
        gf2_copy(&v2, &v0);
        gf2_copy(&u0, &u1);
        gf2_copy(&v0, &v1);

        gf2_mul(&tmp, &Q, &u1);

        gf2_add(&u1, &u2, &tmp);
        gf2_mul(&tmp, &Q, &v1);       //v1 = v2 - qv1
        gf2_add(&v1, &v2 ,&tmp);
        gf2_set_zero(&Q);
    }
    gf2_fit_len(&u0);
    gf2_fit_len(&v0);
    gf2_fit_len(&t0);

    gf2_copy(x, &u0);
    gf2_copy(y, &v0);
    gf2_copy(gcd, &t0);

}

/*
@   dst = difference of src
*/
void gf2_diff(gf2* dst, gf2* src)
{
    if(gf2_is_zero(src) == ZERO){
        gf2_set_zero(dst);
        return;
    }
    int i,j;
    int rq = src->deg / 8;
    int rr = src->deg % 8;

    dst->deg = src->deg - 1;
    for(j=rr; j>=0; j--){
        if((j%2==1) && (((src->binary[rq] >> (j)) & 0x01) == 1)){
            gf2_set_index(dst, 8*rq + j-1);
        }
    }

    for(i=rq-1; i>=0; i--){
        for(j=7; j>=0; j-=2){
            if(((src->binary[i] >> (j)) & 0x01) == 1){
                gf2_set_index(dst, i*8 + j-1);
            }
        }
    }
    gf2_fit_len(dst);

}

/*
@     src : 입력값
@     return : factoring 갯수. 없을 경우 1 반환
*/
int gf2_berlekamp_factoring(gf2* src)
{
    gf2 X;
    gf2 X_powering;
    gf2 prime;
    gf2 gcd;
    BMAT B, I;
    int rank;
    int i;
    int k = src->deg - 1;

    //미분 먼저.
    gf2_init(&prime, src->deg - 1);
    gf2_diff(&prime, src);
    gf2_gcd(&gcd, &prime, src); 
    if(gf2_is_one(&gcd) != ONE) //exists factor
        return REDUCIBLE;

    gf2_init(&X_powering, 1);
    gf2_init(&X, 1);
    gf2_set_index(&X, 1);   //X

    bmatrix_init(B, src->deg - 1, src->deg - 1);
    bmatrix_init(I, src->deg - 1, src->deg - 1);

    for(i=0; i<src->deg - 1; ++i){  //미분값 차수까지
        gf2_powmod(&X_powering, &X, 2*i, src);
        bmatrix_set_gf2(B, X_powering, i);
    }
    bmatrix_generate_identity(I);
    bmatrix_add_identity(B, I);

    rank = bmatrix_rank(B);
    bmatrix_free(B);
    bmatrix_free(I);

    return k - rank;
}

/*
@   num : number (int)
@   gf2 : gf2 struct. binary(char)
*/
gf2 numtogf2(int num)
{
    gf2 res;

    gf2_init(&res, 32);
    
    res.binary[0] = num & 0xff;
    res.binary[1] = (num >> 8) & 0xff;
    res.binary[2] = (num >> 16) & 0xff;
    res.binary[3] = (num >> 24) & 0xff;

    gf2_fit_len(&res);

    return res;
}


/*
@   gf2 : gf2 struct. binary(char)
@   num : number (int)
*/
int gf2tonum(gf2 src)
{
    int res;

    res = src.binary[0] ^ (src.binary[1] << 8) ^ (src.binary[2] << 16) ^ (src.binary[3] << 24);
    
    return res;
}