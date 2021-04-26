#include "McEliece.h"

#include <stdio.h>
#include <string.h>


static void test_generate_paritiycheckmatrix()
{
    int i;
    int res = 0;
    Param ctx;

    memset(&ctx, 0x00, sizeof(Param));

    ctx.n = 70;
    ctx.t = 5;
    ctx.m = 13;

    gf2_init(&ctx.mod, ctx.m);
    gf2_set_index(&ctx.mod, 13);
    gf2_set_index(&ctx.mod, 4);
    gf2_set_index(&ctx.mod, 3);
    gf2_set_index(&ctx.mod, 1);
    gf2_set_index(&ctx.mod, 0);
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
    
    res = get_paritycheck_matrix(ctx);
    
    printf(" !!!! %d \n", res);
    if(ctx.supportSet)  free(ctx.supportSet);
}

void test_mceliece_operation(){

    test_generate_paritiycheckmatrix();

}