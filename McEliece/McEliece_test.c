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

    generateSupportSet(&ctx);
    
    while(get_paritycheck_matrix(&ctx) == FAILURE);
    
    printf(" paritycheck matrix \n");   bmatrix_print(ctx.paritycheckMatrix);

    clearParam(&ctx);
}

static void test_generate_generatormatrix()
{
    int i;
    int res = 0;
    Param ctx;

    memset(&ctx, 0x00, sizeof(Param));

    //2^m >= n && n > mt
    ctx.n = 16;
    ctx.t = 3;
    ctx.m = 4;

    gf2_init(&ctx.mod, ctx.m);
    gf2_generate_irreducible(&ctx.mod, ctx.m);
    gf2_fit_len(&ctx.mod);
    gf2_print(&ctx.mod);

    gf2m_init(&ctx.Goppa, ctx.t);
    gf2m_generate_irreducible(&ctx.Goppa, &ctx.mod, ctx.t);
    printf(" irreducible poly  = ");    gf2m_print(&ctx.Goppa);

    generateSupportSet(&ctx);
    
    while(get_paritycheck_matrix(&ctx) == FAILURE);
    
    printf(" paritycheck matrix \n");   bmatrix_print(ctx.paritycheckMatrix);

    get_generator_matrix(&ctx);
    printf(" generator matrix \n");     bmatrix_print(ctx.generatorMatrix);


    validate(&ctx);
    clearParam(&ctx);
}

static void test_generate_key()
{
    int i;
    int res = 0;
    Param ctx;

    memset(&ctx, 0x00, sizeof(Param));

    //2^m >= n && n > mt
    ctx.n = 16;
    ctx.t = 3;
    ctx.m = 4;

    gf2_init(&ctx.mod, ctx.m);
    gf2_generate_irreducible(&ctx.mod, ctx.m);
    gf2_fit_len(&ctx.mod);
    gf2_print(&ctx.mod);

    gf2m_init(&ctx.Goppa, ctx.t);
    gf2m_generate_irreducible(&ctx.Goppa, &ctx.mod, ctx.t);
    printf(" irreducible poly  = ");    gf2m_print(&ctx.Goppa);

    generateSupportSet(&ctx);

    while(get_paritycheck_matrix(&ctx) == FAILURE);
    
    printf(" paritycheck matrix \n");   bmatrix_print(ctx.paritycheckMatrix);

    get_generator_matrix(&ctx);
    printf(" generator matrix \n");     bmatrix_print(ctx.generatorMatrix);

    get_generate_key(&ctx);
    
    clearParam(&ctx);
}

static void test_patterson_decoding()
{
    int i;
    int res = 0;
    BMAT message;
    BMAT enc;
    BMAT decrypted_message;
    Param ctx;
    gf2 tmp;

    memset(&ctx, 0x00, sizeof(Param));

    //2^m >= n && n > mt
    // ctx.n = 16;
    // ctx.t = 3;
    // ctx.m = 4;

    ctx.n = 42;
    ctx.t = 3;
    ctx.m = 13;

    gf2_init(&ctx.mod, ctx.m);
    //gf2_generate_irreducible(&ctx.mod, ctx.m);
    gf2_set_index(&ctx.mod, 13);
    gf2_set_index(&ctx.mod, 4);
    gf2_set_index(&ctx.mod, 3);
    gf2_set_index(&ctx.mod, 1);
    gf2_set_index(&ctx.mod, 0);
    
    gf2_fit_len(&ctx.mod);
    gf2_print(&ctx.mod);

    gf2m_init(&ctx.Goppa, ctx.t);
    gf2_init(&tmp, ctx.m);
    //gf2m_generate_irreducible(&ctx.Goppa, &ctx.mod, ctx.t);
    // irreducible poly
    //(z^0 )*X^3+(z^9 +z^7 +z^5 +z^4 +z^3 +z^2 +z^0 )*X^2+(z^10 +z^7 +z^5 +z^3 +z^2 +z^0 )*X^1+(z^9 +z^7 +z^5 +z^2 +z^1 )*X^0

    gf2_set_index(&tmp, 0);
    gf2m_set_index(&ctx.Goppa, &tmp, 3);

    gf2_set_index(&tmp, 9);
    gf2_set_index(&tmp, 7);
    gf2_set_index(&tmp, 5);
    gf2_set_index(&tmp, 4);
    gf2_set_index(&tmp, 3);
    gf2_set_index(&tmp, 2);
    gf2_fit_len(&tmp);
    gf2m_set_index(&ctx.Goppa, &tmp, 2);

    gf2_set_zero(&tmp);
    gf2_set_index(&tmp, 10);
    gf2_set_index(&tmp, 7);
    gf2_set_index(&tmp, 5);
    gf2_set_index(&tmp, 3);
    gf2_set_index(&tmp, 2);
    gf2_set_index(&tmp, 0);
    gf2_fit_len(&tmp);
    gf2m_set_index(&ctx.Goppa, &tmp, 1);

    gf2_set_zero(&tmp);
    gf2_set_index(&tmp, 9);
    gf2_set_index(&tmp, 7);
    gf2_set_index(&tmp, 5);
    gf2_set_index(&tmp, 2);
    gf2_set_index(&tmp, 1);
    gf2_fit_len(&tmp);
    gf2m_set_index(&ctx.Goppa, &tmp, 0);
    gf2m_fit_len(&ctx.Goppa);

    printf(" irreducible poly  = ");    gf2m_print(&ctx.Goppa);

    return ;
    generateSupportSet(&ctx);
    
    while(get_paritycheck_matrix(&ctx) == FAILURE);
    
    printf(" paritycheck matrix \n");   bmatrix_print(ctx.paritycheckMatrix);

    get_generator_matrix(&ctx);
    printf(" generator matrix \n");     bmatrix_print(ctx.generatorMatrix);


    bmatrix_init(message, 1, ctx.n - ctx.m * ctx.t);
    generate_random_bmatrix(message);

    printf("message \n");
    bmatrix_print(message);

    printf("key generation \n");
    get_generate_key(&ctx);

    printf("encryption \n");
    encryption(enc, message, &ctx);

    printf("decryption \n");
    decryption(decrypted_message, enc, &ctx);
    //추가 예정

    clearParam(&ctx);

}

//(z^12 +z^11 +z^10 +z^7 +z^5 +z^4 +z^3 +z^1 +z^0 )*X^5
//+(z^11 +z^10 +z^6 +z^5 +z^3 +z^1 +z^0 )*X^4
//+(z^13 +z^9 +z^8 +z^7 +z^6 +z^5 +z^2 +z^1 +z^0 )*X^3
//+(z^13 +z^11 +z^10 +z^9 +z^8 +z^6 +z^5 +z^4 +z^3 +z^1 )*X^2
//+(z^12 +z^10 +z^9 +z^8 +z^7 +z^5 +z^4 +z^2 +z^0 )*X^1
//+(z^12 +z^11 +z^10 +z^9 +z^8 +z^7 +z^6 +z^5 +z^3 +z^2 )*X^0
//roots : [(z^10 + z^9 + z^7 + z^6 + z^2 + 1, 1),
// (z^12 + z^11 + z^8 + z^6 + z^5 + z^4 + z^3 + z^2 + z + 1, 1)]
void test_find_root()
{
    gf2m A;
    gf2 test;
    gf2 mod;
    //int* support = NULL;
    int i;
    Param ctx;
    

    //2^m >= n && n > mt
    
    int t = 5;
    int m = 13;
    int n = 1<<m;
    
    memset(&ctx, 0x00, sizeof(Param));

    gf2_init(&mod, 13);
    gf2_set_index(&mod, 13);
    gf2_set_index(&mod, 4);
    gf2_set_index(&mod, 3);
    gf2_set_index(&mod, 1);
    gf2_set_index(&mod, 0);
    gf2_fit_len(&mod);
    
    printf("mod\n");
    gf2_print(&mod);
    printf("***********\n");

    gf2m_init(&A, t);
    //gf2m_random_gen(&A, m);


//(z^12 +z^11 +z^10 +z^7 +z^5 +z^4 +z^3 +z^1 +z^0 )*X^5
    gf2_init(&test, 13);
    gf2_set_index(&test, 12);
    gf2_set_index(&test, 11);
    gf2_set_index(&test, 10);
    gf2_set_index(&test, 7);
    gf2_set_index(&test, 5);
    gf2_set_index(&test, 4);
    gf2_set_index(&test, 3);
    gf2_set_index(&test, 1);
    gf2_set_index(&test, 0);
    gf2_fit_len(&test);
    gf2m_set_index(&A, &test, 5);

//+(z^11 +z^10 +z^6 +z^5 +z^3 +z^1 +z^0 )*X^4

    gf2_set_zero(&test);
    test.deg = 14;
    gf2_set_index(&test, 11);
    gf2_set_index(&test, 10);
    gf2_set_index(&test, 6);
    gf2_set_index(&test, 5);
    gf2_set_index(&test, 3);
    gf2_set_index(&test, 1);
    gf2_set_index(&test, 0);
    gf2_fit_len(&test);
    gf2m_set_index(&A, &test, 4);

//+(z^13 +z^9 +z^8 +z^7 +z^6 +z^5 +z^2 +z^1 +z^0 )*X^3
    gf2_set_zero(&test);
    test.deg = 14;
    gf2_set_index(&test, 13);
    gf2_set_index(&test, 9);
    gf2_set_index(&test, 8);
    gf2_set_index(&test, 7);
    gf2_set_index(&test, 6);
    gf2_set_index(&test, 5);
    gf2_set_index(&test, 2);
    gf2_set_index(&test, 1);
    gf2_set_index(&test, 0);
    gf2_fit_len(&test);
    gf2m_set_index(&A, &test, 3);

//+(z^13 +z^11 +z^10 +z^9 +z^8 +z^6 +z^5 +z^4 +z^3 +z^1 )*X^2

    gf2_set_zero(&test);
    test.deg = 14;
    gf2_set_index(&test, 13);
    gf2_set_index(&test, 11);
    gf2_set_index(&test, 10);
    gf2_set_index(&test, 9);
    gf2_set_index(&test, 8);
    gf2_set_index(&test, 6);
    gf2_set_index(&test, 5);
    gf2_set_index(&test, 4);
    gf2_set_index(&test, 3);
    gf2_set_index(&test, 1);
    gf2_fit_len(&test);
    gf2m_set_index(&A, &test, 2);

//+(z^12 +z^10 +z^9 +z^8 +z^7 +z^5 +z^4 +z^2 +z^0 )*X^1

    gf2_set_zero(&test);
    test.deg = 14;
    gf2_set_index(&test, 12);
    gf2_set_index(&test, 10);
    gf2_set_index(&test, 9);
    gf2_set_index(&test, 8);
    gf2_set_index(&test, 7);
    gf2_set_index(&test, 5);
    gf2_set_index(&test, 4);
    gf2_set_index(&test, 2);
    gf2_set_index(&test, 0);
    gf2_fit_len(&test);
    gf2m_set_index(&A, &test, 1);

    
//+(z^12 +z^11 +z^10 +z^9 +z^8 +z^7 +z^6 +z^5 +z^3 +z^2 )*X^0
    gf2_set_zero(&test);
    test.deg = 14;
    gf2_set_index(&test, 12);
    gf2_set_index(&test, 11);
    gf2_set_index(&test, 10);
    gf2_set_index(&test, 9);
    gf2_set_index(&test, 8);
    gf2_set_index(&test, 7);
    gf2_set_index(&test, 6);
    gf2_set_index(&test, 5);
    gf2_set_index(&test, 3);
    gf2_set_index(&test, 2);
    gf2_fit_len(&test);
    gf2m_set_index(&A, &test, 0);

    printf("A(X)\n");
    gf2m_print(&A);

    ctx.supportSet = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++)
    {
        *(ctx.supportSet + i) = i;
    }
    ctx.n = n;
    ctx.mod = mod;

    gf2 *tmp = NULL;
    tmp = (gf2*)calloc(n, sizeof(gf2));
    int sol_num = find_root(tmp, &A, &ctx);

    i = 0;
    printf("root_set %d개\n", sol_num);
    while (sol_num--)
    {
        gf2_print(&tmp[sol_num]);
    }

    if (ctx.supportSet) free(ctx.supportSet);
    if (tmp) free(tmp);
}

void test_mceliece_operation(){

    //test_generate_paritiycheckmatrix();
    //test_generate_generatormatrix();
    //test_generate_key();
    test_patterson_decoding();
    //test_find_root();
}