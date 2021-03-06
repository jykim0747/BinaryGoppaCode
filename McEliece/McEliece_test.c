#include "McEliece.h"

#include <stdio.h>
#include <string.h>


static void test_generate_paritiycheckmatrix()
{
    int i;
    int res = 0;
    Param ctx;

    memset(&ctx, 0x00, sizeof(Param));

    //2^m >= n && n > mt
    ctx.n = 16;
    ctx.t = 3;
    ctx.m = 5;

    gf2_init(&ctx.mod, ctx.m);
    gf2_generate_irreducible(&ctx.mod, ctx.m);
    gf2_fit_len(&ctx.mod);
    gf2_print(&ctx.mod);

    gf2m_init(&ctx.Goppa, ctx.t);
    gf2m_generate_irreducible(&ctx.Goppa, &ctx.mod, ctx.t);
    printf(" irreducible poly  = ");    gf2m_print(&ctx.Goppa);

    ctx.supportSet = (int*)malloc(sizeof(int) * (1<<ctx.m));
    for(i=0; i<(1<<ctx.m); ++i)
    {
        ctx.supportSet[i] = i;
    }
    Fisher_Yate(ctx.supportSet, 1<<ctx.m);
    
    while(get_paritycheck_matrix(&ctx) == FAILURE);
    
    printf(" paritycheck matrix \n");   bmatrix_print(ctx.paritycheckMatrix);

    bmatrix_free(ctx.paritycheckMatrix);
    if(ctx.supportSet)  free(ctx.supportSet);
}

void test_mceliece_operation(){

    test_generate_paritiycheckmatrix();

}