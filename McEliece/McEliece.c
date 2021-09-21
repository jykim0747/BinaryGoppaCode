#include "McEliece.h"

static void generateP(BMAT src, int num);

/*
@   len 길이를 가지는 set의 원소를 섞는 함수
@   set : 섞이는 집합. 갱신
@   len : set 길이(섞을 범위)
*/
void Fisher_Yate(int* set, int len)
{
    int i, cv;
    int tmp;
    for(i=len-1; i>=0; i--)
    {
        cv = rand() % len;
        tmp = *(set + cv);
        *(set + cv) = *(set + i);
        *(set + i) = tmp;
    }
}


void clearParam(Param* ctx){
    gf2_set_zero(&ctx->mod);
    gf2m_set_zero(&ctx->Goppa);
    if(ctx->supportSet) free(ctx->supportSet);
    bmatrix_free(ctx->paritycheckMatrix);
    bmatrix_free(ctx->generatorMatrix);

    bmatrix_free(ctx->key.P);
    bmatrix_free(ctx->key.S);
    bmatrix_free(ctx->key.SGP);
}

int get_paritycheck_matrix(Param* ctx)
{
    int i, ir, iq;
    int iter, jter;
    gf2m Lpoly;
    gf2 X;
    BMAT mat;
    BMAT mat_ech;
    int res;

    gf2_init(&X, 1);
    bmatrix_init(mat, ctx->t * ctx->m, ctx->n);

    for (i = 0; i < ctx->n; ++i)
    {

        /*
        Step 1: Compute a(i) in L = [a(0),a(1),a(2),...,a(n-1)]
        */

        X = numtogf2(ctx->supportSet[i]);

        /*
        Step 2: Set (X - a(i)) = Lpoly
        */

        gf2m_init(&Lpoly, 1);
        gf2_set_one(&Lpoly.term[1]);
        gf2m_set_index(&Lpoly, &X, 0);

        // printf("Lpoly (%d): ",i);   gf2m_print(&Lpoly);

        /*
        Step 3: Compute (X - a(i))^-1 in terms of g(X)
        A * Goppa_poly + L_inv * LPoly = val
        We use the extended Euclidean algorithm
        */
        gf2m inv;
        gf2m tmp;
        gf2m gcd;

        gf2m_init(&inv, 1);
        gf2m_init(&tmp, 1);
        gf2m_init(&gcd, 1);

        gf2m_xgcd(&gcd, &Lpoly, &ctx->Goppa, &inv, &tmp, &ctx->mod);

        /*
        Step 4: Generate the Parity-check matrix H from (X - a(i))^-1
        역원을 각 열에 대입
        */

        gf2m_to_bmat(mat, inv, ctx->m, i);
    }

    //echelon form 을 적용하여 항등 행렬 I를 가지는지 확인
    bmatrix_init(mat_ech, ctx->t * ctx->m, ctx->n);
    bmatrix_echelon(mat_ech, mat);

    if(bmatrix_has_zero_rows(mat_ech) == ZERO)
    {
        Fisher_Yate(ctx->supportSet, 1<<ctx->m);
        res = FAILURE;
        goto end;
    }

    //중복 코드 처리 필요
    if(has_Identity_bmat(mat_ech) == IDENTITY)
    {
        bmatrix_init(ctx->paritycheckMatrix, ctx->t * ctx->m, ctx->n);
        bmatrix_copy(ctx->paritycheckMatrix, mat);
        res = SUCCESS;
        goto end;
    }
    else{
        make_Identity_bmat(mat_ech, ctx->supportSet);
        res = FAILURE;
    }

end:
    bmatrix_free(mat_ech);
    bmatrix_free(mat);

    return res;
}

int get_generator_matrix(Param* ctx)
{
    BMAT mat = {0x0,};
    BMAT mat_ech = {0x0,};
    BMAT mat_tmp = {0x0,};
    BMAT mat_transpose = {0x0,};
    BMAT mat_I = {0x0,};
    BMAT* pm;
    int res;
    int i = 0;
    int k = 0;

    /* [k x n(=mt+k)], gmat = [gmat_tmp || gmat_I] */

    k = (ctx->paritycheckMatrix->c - ctx->paritycheckMatrix->r);

    pm = &ctx->paritycheckMatrix;
    bmatrix_init(mat, k, ctx->n - k);
    bmatrix_init(ctx->generatorMatrix, k, ctx->n);
    bmatrix_init(mat_ech, (*pm)->r, (*pm)->c);
    bmatrix_init(mat_transpose, (*pm)->c, (*pm)->r);
    
    bmatrix_echelon(mat_ech, *pm);
    bmatrix_transpose(mat_transpose, mat_ech);

    for(i = 0; i < k; ++i){
        mat->data[i] = mat_transpose->data[mat_transpose->c + i];
    }

    bmatrix_init(mat_I, k, k);
    bmatrix_generate_identity(mat_I);

    mat_concat_horizontal(ctx->generatorMatrix, mat, mat_I);
    
    bmatrix_free(mat);
    bmatrix_free(mat_tmp);
    bmatrix_free(mat_ech);
    bmatrix_free(mat_I);
    bmatrix_free(mat_transpose);

    return 0;
}

// G*H^T = 0
int validate(Param* ctx)
{
    BMAT z;
    BMAT pt;
    int k = (ctx->paritycheckMatrix->c - ctx->paritycheckMatrix->r);
    int n = ctx->n;

    bmatrix_init(pt, ctx->paritycheckMatrix->c, ctx->paritycheckMatrix->r);
    bmatrix_transpose(pt, ctx->paritycheckMatrix);

    bmatrix_init(z, k, n-k);
    bmatrix_mul(z, ctx->generatorMatrix, pt);

    bmatrix_free(pt);
    bmatrix_free(z);

    return 0;
}

int get_generate_key(Param* ctx){
    int res = 0;

    int k = ctx->generatorMatrix->r;
    int n = ctx->generatorMatrix->c;
    BMAT sInv;
    BMAT tmp;

    /*
        Generation of (k x k) random invertible matrix S
    */
    bmatrix_init(ctx->key.S, k, k);
    bmatrix_init(sInv, k, k);

    bmatrix_generate_inverse(sInv, ctx->key.S);
#ifdef DEBUG
    bmatrix_generate_identity(ctx->key.S);
#endif
    printf("S\n"); bmatrix_print(ctx->key.S);
    /*
        Generation of (n x n) random permutation matrix P
    */
    bmatrix_init(ctx->key.P, n, n);
    generateP(ctx->key.P, n);
#ifdef DEBUG
    bmatrix_generate_identity(ctx->key.P);
#endif
    printf("P\n"); bmatrix_print(ctx->key.P);

    /*
        Generation of G' = SGP
    */
    bmatrix_init(ctx->key.SGP, k, n);
    bmatrix_init(tmp, k, n);

    bmatrix_mul(tmp, ctx->key.S, ctx->generatorMatrix);
    bmatrix_mul(ctx->key.SGP, tmp, ctx->key.P);
    printf("SGP\n"); bmatrix_print(ctx->key.SGP);

    bmatrix_free(tmp);
    bmatrix_free(sInv);
    return res;
}

static void generateP(BMAT src, int num){

    int* set = NULL;
    int i = 0;
    int pq, pr;
    
    set = (int*)malloc(sizeof(int) * num);
    
    for(i=0; i<num; ++i){
        *(set + i) = i;
    }

    Fisher_Yate(set, num);

    for (i = 0; i < num; ++i)
    {
        pq = (*(set + i)) / (sizeof(unsigned char) * 8);
        pr = (*(set + i)) % (sizeof(unsigned char) * 8);
        b_mat_entry(src, i, pq) ^= (1 << pr);
    }

    if(set) free(set);
}