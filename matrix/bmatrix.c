#include "matrix.h"
#include "error.h"
#include "stdlib.h"

void bmatrix_init(BMAT mat, int rows, int cols)
{
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