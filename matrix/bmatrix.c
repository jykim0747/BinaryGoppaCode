#include "matrix.h"
#include "error.h"
#include "stdlib.h"

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
    int i, j, res;
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

    if(bmatrix_has_zero_rows(tmp) == ZERO)
    {
        printf("has zero rows\n");
        return FAILURE;
    }

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