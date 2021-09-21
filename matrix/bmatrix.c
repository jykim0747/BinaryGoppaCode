#include "bmatrix.h"
#include "error.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void swap(int* a, int* b)
{
    int tmp;

    tmp = *a;
    *a = *b;
    *b = tmp;
}

void bmatrix_init(BMAT mat, int rows, int cols)
{
    if((rows == 0) || (cols == 0))
        return;

    mat->cnum = (cols + 7) / 8;

    if (rows != 0 && cols != 0) /* Allocate space for r*c small entries */
    {
        int i;
        mat->entries = (unsigned char *)calloc((rows * mat->cnum), sizeof(unsigned char));
        mat->data = (unsigned char **)malloc(rows * sizeof(unsigned char *)); /* Initialise rows */

        for (i = 0; i < rows; i++)
            mat->data[i] = mat->entries + i * mat->cnum;
    }
    else
    {
        mat->entries = NULL;
        mat->data = NULL;
    }

    mat->r = rows;
    mat->c = cols;
}

void bmatrix_print(const BMAT mat)
{
    int i, j, k;

    printf("A = matrix(GF(2),[\n");
    for (i = 0; i < mat->r; i++)
    {
        printf("[");
        for (j = 0; j < mat->cnum - 1; j++)
        {
            for (k = 7; k >= 0; k--)
            {
                printf("%x,", (b_mat_entry(mat, i, j) >> k) & 0x1);
            }
            if (j < mat->c - 1)
                printf(" ");
        }
        for (j = 0; j<= (mat->c + 7) % 8 ; j++)
        {
            printf("%x,", (b_mat_entry(mat, i, mat->cnum - 1) >> (7-j)) & 0x1);
        }
        printf("],\n");
    }
    printf("])\n");
}

void generate_random_bmatrix(BMAT mat)
{
    int i, j;
    int r = mat->c % 8;
    unsigned char mask[8] = {0xff, 0x80, 0xc0, 0xe0, 0xf0, 0xf8, 0xfc, 0xfe};
    for(i=0; i<mat->r; ++i){
        for(j=0; j<mat->cnum-1; ++j){
            b_mat_entry(mat, i, j) = rand();
        }
        b_mat_entry(mat, i, mat->cnum-1) = rand()%mask[r];
    }
}

void bmatrix_free(BMAT mat)
{
    if (mat->entries)
    {
        free(mat->entries);
        free(mat->data);   
    }
}

void bmatrix_swap_rows(BMAT mat, int row1, int row2)
{
    if (row1 != row2)
    {
        unsigned char * u;

        u = mat->data[row2];
        mat->data[row2] = mat->data[row1];
        mat->data[row1] = u; 
    }
}

void bmatrix_copy(BMAT dst, BMAT src)
{
    int i, j;

    if((src == NULL) || (dst == NULL))
    {
        return;
    }

    if((src->r != dst->r) || (src->c != dst->c))
    {
        return;
    }

    for(i=0; i<src->r; ++i)
        for(j=0; j<src->cnum; ++j)
            b_mat_entry(dst, i, j) = b_mat_entry(src, i, j);

}

int bmatrix_has_zero_rows(BMAT mat)
{
    int i, j;
    int count;

    for(i=0; i<mat->r; ++i){
        count = 0;
        for(j=0; j<mat->cnum; ++j){
            if(b_mat_entry(mat, i, j) == ZERO)
                count++;
        }
        if(count == mat->cnum)
            return ZERO;   
    }

    return NOT_ZERO;
}

/*
* mat->data[row1] = mat->data[row1] ^ mat->data[row2]
*/
void bmatrix_add_row(BMAT mat, int row1, int row2)
{
    int i;
    for(i=0; i<mat->cnum; ++i)
        b_mat_entry(mat, row1, i) ^= b_mat_entry(mat, row2, i);
}

int bmatrix_echelon(BMAT mat_ech, BMAT mat)
{
    BMAT tmp;
    int i, j, k, count;

    bmatrix_init(tmp, mat->r, mat->c);
    bmatrix_copy(tmp, mat);

    // if(bmatrix_has_zero_rows(tmp) == ZERO)
    // {
    //     printf("has zero rows\n");
    //     return FAILURE;
    // }

    for(i=0; i<mat->r; ++i)
    {
        int iq = i / 8;
        int ir = i % 8;
        count = 0;
        for(j = i; j<mat->r; j++)
        {
            if((b_mat_entry(tmp, j, iq)>>(7-ir)) & 0x01)
            {
                if(count == 0){
                    bmatrix_swap_rows(tmp, i, j);
                    count = 1;
                }
                else{
                    bmatrix_add_row(tmp, j, i);
                }
            }
        }
        for(k=0; k<i; k++){
            if((b_mat_entry(tmp, k, iq)>>(7-ir)) & 0x01){
                bmatrix_add_row(tmp, k, i);
            }
        }
    }
    bmatrix_copy(mat_ech, tmp);
    bmatrix_free(tmp);

    return SUCCESS;
}

void bmatrix_set_gf2(BMAT mat, const gf2 src, int row)
{
    int i, j;
    int qq = src.deg / 8;

    for(i=qq; i>=0; i--)
    {
        for(j=7; j>=0; j--)
        {
            if((src.binary[i]>>(7-j)) & 0x01){
                b_mat_entry(mat, row, i) ^= (1<<j);
            }
        }
    }
}

void bmatrix_generate_identity(BMAT mat)
{
    int i,j;
    int rq = mat->r / 8;
    int rr = mat->r % 8;

    if(mat->r != mat->c)
        return;

    for(i=0; i< rq; ++i)
    {
        for(j=0; j<8; ++j)
            b_mat_entry(mat, i*8 + j, i) ^= (1<<(7-j));
    }
    
    for(j=0; j< rr; ++j)
    {
        b_mat_entry(mat, rq*8 + j, rq) ^= (1<<(7-j));
    }

}

/*
*   src : identity bmatrix
*/
void bmatrix_add_identity(BMAT dst, BMAT src)
{
    int i, j;

    for(i=0; i<src->r; ++i)
        for(j=0; j<src->cnum; ++j)
            b_mat_entry(dst, i, j) ^= b_mat_entry(src, i, j);

}

/*
*   counting diagonal elements after Gaussian elimination
*   Sometimes it returns a different result
*   But I want full rank or not
*/
int bmatrix_rank(BMAT mat)
{
    int rank = 0;
    int i, j;
    BMAT tmp;

    int rq = mat->r / 8;
    int rr = mat->r % 8;

    bmatrix_init(tmp, mat->r, mat->c);
    bmatrix_echelon(tmp, mat);

    for(i=0; i< rq; ++i){
        for(j=0; j<8; ++j){
            if(((b_mat_entry(tmp, i*8 + j, i) >> (7-j)) & 0x01) == 1)
                rank++;
        }
    }
    
    for(j=0; j< rr; ++j){
        if(((b_mat_entry(tmp, rq*8 + j, rq) >> (7-j)) & 0x01) == 1)
            rank++;
    }

    bmatrix_free(tmp);

    return rank;
}

/*
*   mat = [I|M]인지 확인하는 함수.
*/
int has_Identity_bmat(BMAT mat)
{
    int i,j;
    int rq = mat->r / 8;
    int rr = mat->r % 8;
    int count = 0;

    for(i=0; i< rq; ++i)
    {
        for(j=0; j<8; ++j){
            if((b_mat_entry(mat, i*8 + j, i) >> (7-j)) & 0x01 )
                count++;
        }
    }
    
    for(j=0; j< rr; ++j)
    {
        if((b_mat_entry(mat, rq*8 + j, rq) >> (7-j)) & 0x01)
            count++;
    }

    if(count == mat->r)
        return IDENTITY;

    return NOT_IDENTITY;
}

void make_Identity_bmat(BMAT mat, int *support)
{
    int i,j;
    int iq, ir;

    for(i=0; i< mat->r; ++i)
    {
        iq = i / 8;
        ir = i % 8;
        if(!((b_mat_entry(mat, i, iq) >> (7-ir)) & 0x01) )
        {
            for(j=ir; j<mat->c; ++j){
                if((b_mat_entry(mat, i, j/8) >> (7-j%8) & 0x01)){
                    swap(support+i, support+j);
                    break;
                }
            }
        }
    }
}

/**
 * gf2m -> bmatrix 세로로 값을 넣는 함수
 */
void gf2m_to_bmat(BMAT mat, gf2m src, int m, int column)
{

    int cq = column / 8;
    int cr = column % 8;
    int i, j;
    int deg = m;
    int gf2;

    for (i = 0; i <= src.deg; ++i)
    {
        gf2 = gf2tonum(src.term[i]);
        for (j = 0; j < deg; ++j)
        {
            if (((gf2 >> j) & 0x01) == 0x1)
            {
                b_mat_entry(mat, i * deg + j, cq) ^= 1 << (7 - cr);
            }
        }
    }
}
/*
@   가로로 연접하는 함수
@   dst : dst <- [a||b]
@     a : 왼쪽
@     b : 오른쪽
*/
int mat_concat_horizontal(BMAT dst, BMAT a, BMAT b)
{
    int i, j, k;
    int r = a->c % 8;

    if (a->r != b->r)
        return FAILURE;

    if (r == 0){
        for (i = 0; i < a->r; ++i){
            for(k = 0; k<a->cnum; ++k){
                *(dst->data[i] + k) = *(a->data[i] + k);
            }
            for (j = 0; j < b->cnum; ++j){
                *(dst->data[i] + a->cnum + j) = *(b->data[i] + j);
            }
        }
    }
    else
    {
        for (i = 0; i < a->r; ++i){
            for(k = 0; k<a->cnum; ++k){
                *(dst->data[i] + k) = *(a->data[i] + k);
            }
            *(dst->data[i] + (a->cnum-1)) &= (unsigned char)~((1 << (8 - r)) - 1); //a 행렬에서 불필요한 부분 제거
            *(dst->data[i] + (a->cnum-1)) ^= *(b->data[i]) >> r;

            for (j = 0; j < b->cnum; ++j){//a,b의 결합 이후 부분. b의 나머지 부분
                *(dst->data[i] + (a->cnum) + j) = *(b->data[i]+j) << (8-r);
                *(dst->data[i] + (a->cnum) + j) ^= *(b->data[i]+j+1) >> r;
            }
        }
    }

    return SUCCESS;
}

int bmatrix_transpose(BMAT dst, BMAT src){

    int res = 0;
    int iter = 0;
    int j = 0;
    int k = 0;
    int cr = src->c & 0x07;
 
    if(src->r != dst->c) goto end;
    if(src->c != dst->r) goto end;

    for (iter = 0; iter < src->r; ++iter)
    {
        for (j = 0; j < src->cnum-1; ++j)
        {
            for (k = 0; k < 8; ++k)
            {
                b_mat_entry(dst, j*8 + k, iter/8) ^= (b_mat_entry(src, iter, j) << k & 0x80) >> (iter%8);
            }
        }

        cr = (cr == 0) ? 8 : cr;
        for (k = 0; k < cr; ++k)
        {
            b_mat_entry(dst, (src->cnum-1)*8+k, iter/8) ^= (b_mat_entry(src, iter, (src->cnum - 1)) << k & 0x80) >> (iter%8);
        }
    }

end:

    return res;
}

int bmatrix_mul(BMAT dst, BMAT src1, BMAT src2)
{
    int res = 0;
    int r = src1->c % 8;
    int i,j,k,l;

    if (src1->c != src2->r)
    {
        res = -2;
        goto end;
    }
    printf("r, src1 cnum: %d, %d\n", r, src1->cnum);
    for(i=0; i<src1->r; ++i){
        for(j=0; j<src1->cnum-1; ++j){
            for(k=0; k<8; ++k){
                if(((b_mat_entry(src1, i, j)>>(7-k)) & 0x01) == 1){
                for(l=src2->cnum-1; l>=0; l--){
                    b_mat_entry(dst, i, l) ^= b_mat_entry(src2, j*8 + k, l);
                }
                }
            }
        }
        for(k=0; k<r; ++k){
            if(((b_mat_entry(src1, i, src1->cnum-1)>>(7-k)) & 0x01) == 1){
            for(l=src2->cnum-1; l>=0; l--){
                b_mat_entry(dst, i, l) ^= b_mat_entry(src2, (src1->cnum-1)*8 + k, l);
            }
            }
        }
    }

end:
    return res;
}

/*
* dst = inverse of src
* SUCCESS : 역행렬 존재
* FAILURE : 역행렬 존재하지 않음
*/
int bmatrix_inverse(BMAT dst, BMAT src){
    int res = 0;
    BMAT tmp;
    BMAT I;
    int i, j, k, count;

    bmatrix_init(tmp, src->r, src->c);
    bmatrix_copy(tmp, src);

    if(bmatrix_has_zero_rows(tmp) == ZERO)
    {
        //printf("has zero rows\n");
        res = FAILURE;
        goto end;
    }

    bmatrix_init(I, src->r, src->c);
    bmatrix_generate_identity(I);

    for(i=0; i<src->r; ++i)
    {
        int iq = i / 8;
        int ir = i % 8;
        count = 0;
        for(j = i; j<src->r; j++)
        {
            if((b_mat_entry(tmp, j, iq)>>(7-ir)) & 0x01)
            {
                if(count == 0){
                    bmatrix_swap_rows(tmp, i, j);
                    bmatrix_swap_rows(I, i, j);
                    count = 1;
                }
                else{
                    bmatrix_add_row(tmp, j, i);
                    bmatrix_add_row(I, j, i);
                }
            }
        }
        for(k=0; k<i; k++){
            if((b_mat_entry(tmp, k, iq)>>(7-ir)) & 0x01){
                bmatrix_add_row(tmp, k, i);
                bmatrix_add_row(I, k, i);
            }
        }
    }

    if(has_Identity_bmat(tmp) == NOT_IDENTITY)
    {
        res = FAILURE;
        goto end;
    }    

    bmatrix_copy(dst, I);

end:
    bmatrix_free(tmp);
    bmatrix_free(I);

    return res;
}

void bmatrix_generate_inverse(BMAT dst, BMAT src){
    do{
        generate_random_bmatrix(src);
    }
    while(bmatrix_inverse(dst, src) == FAILURE);
}