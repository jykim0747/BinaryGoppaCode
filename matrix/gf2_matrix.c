#include "matrix.h"
#include "error.h"
#include "stdlib.h"

void gf2_matrix_print(const gf2_MAT mat)
{
    int i, j;

    printf("A = matrix([\n");
    for (i = 0; i < mat->r; i++)
    {
        printf("[");
        for (j = 0; j < mat->c; j++)
        {
            gf2_print_pretty(gf2_mat_entry(mat, i, j));
            if (j + 1 < mat->c)
                printf(", ");
        }
        printf("],\n");
    }
    printf("])\n");
}

void gf2_matrix_init(gf2_MAT mat, int rows, int cols, int m)
{
    int i;
    if (rows != 0)
        mat->data = (gf2 **) malloc(rows * sizeof(gf2 *));
    else
        mat->data = NULL;

    if (rows != 0 && cols != 0)
    {
        mat->entries = (gf2 *) malloc((rows * cols) * sizeof(gf2));

        for (i = 0; i < rows * cols; ++i)
            gf2_init(mat->entries + i, m);

        for (i = 0; i < rows; i++)
            mat->data[i] = mat->entries + i * cols;
    }
    else
    {
        mat->entries = NULL;
        if (rows != 0)
        {
            for (i = 0; i < rows; i++)
                mat->data[i] = NULL;
        }
    }

    mat->r = rows;
    mat->c = cols;
}

gf2* gf2_mat_entry(gf2_MAT mat, int row, int col)
{
    return (mat->data[row] + col);
}

void generate_random_gf2_matrix(gf2_MAT mat)
{
    int i, j;

    for(i=0; i<mat->r; ++i){
        for(j=0; j<mat->c; ++j){
            gf2_random_gen(&mat->entries[i*mat->r + j]);
        }
    }

}

void gf2_matrix_free(gf2_MAT mat)
{
    if (mat->entries)
    {
        int i;

        for (i = 0; i < mat->r * mat->c; i++)
            gf2_set_zero(mat->entries + i);

        free(mat->entries);
        free(mat->data);
    } else if (mat->r != 0)
        free(mat->data);
}

void gf2_matrix_swap_rows(gf2_MAT mat, int row1, int row2)
{
    if (row1 != row2)
    {
        gf2 * u;

        u = mat->data[row2];
        mat->data[row2] = mat->data[row1];
        mat->data[row1] = u; 
    }
}