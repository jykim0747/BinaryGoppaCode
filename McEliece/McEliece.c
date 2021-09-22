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

int generateSupportSet(Param* ctx){
    int i = 0;
    int m = ctx->m;

    ctx->supportSet = (int*)malloc(sizeof(int) * (1<<m));
    for(i = 0; i < (1<<m); ++i)
    {
        ctx->supportSet[i] = i;
    }
    Fisher_Yate(ctx->supportSet, 1<<m);

    return 0;
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
        pq = *(set + i) / 8;
        pr = *(set + i) % 8;
        b_mat_entry(src, i, pq) ^= (1 << pr);
    }

    if(set) free(set);
}

/*
@   오류 생성 함수
@   error : 오류 벡터
@     num : 오류 개수
*/
void generateError(BMAT error, int num)
{
    int* error_index = NULL;
    int i;
    int q, r;

    error_index = (int*)malloc(sizeof(int) * error->c);

    for(i=0; i<error->c; i++)
        *(error_index+i) = i;

    Fisher_Yate(error_index, error->c);

    for(i=0; i<num; i++)
    {   
        q = *(error_index+i)/8;
        r = *(error_index+i)%8;
        b_mat_entry(error, 0, q) ^= (1 << r);
    }

    if(error_index) free(error_index);

}

/*
@   McEliece Encryption
@   dst : ciphertext
@   src : plaintext
@   ctx : context
*/
void encryption(BMAT dst, BMAT src, Param* ctx)
{
    BMAT errCode;

    bmatrix_init(errCode, 1, ctx->key.SGP->c);
    generateError(errCode, ctx->t);
    
    /*
        Generation of a codeword.
    */
    bmatrix_init(dst, 1, ctx->key.SGP->c);
    bmatrix_mul(dst, src, ctx->key.SGP);

    /*
        Codeword + Error vector
    */
    printf("Ciphertext = codeword + error\n");
    //mat_add(dst, dst, errCode); //구현 필요

    bmatrix_free(errCode);
}

/*
@   patterson decoding 알고리듬
@   src : 수신된 벡터
@   pmat : paritycheck matrix
@   support : support set
@   Goppa_poly : Goppa 다항식
@   mod : mod(fq)
*/
int patterson_decoding(BMAT dst, BMAT src, Param* ctx)
{
    int res = 0;
    int i;
    BMAT received_vec_tp;
    BMAT syndrome;
    
    gf2m synPoly;
    gf2m gcd, inv, tmp;
    gf2m X, synPolyInv;
    gf2m polyTmp;
    gf2m Ax, Bx, Dx;
    gf2m sigma, Ax2, Bx2;
    gf2 rootSet;

    /***********************************************
     * Step 1. Compute the syndrome: Hy^T
    ***********************************************/    

    bmatrix_init(received_vec_tp, src->c, src->r);
    bmatrix_init(syndrome, ctx->paritycheckMatrix->r, 1);
    
    bmatrix_mul(syndrome, ctx->paritycheckMatrix, received_vec_tp);
    bmatrix_free(received_vec_tp);

    //mat is zero 구현 필요
    /*
    if( mat_is_zero(syndrome) == ZERO ) //if syndrome is zero. -> return No-Error.
    {
        printf("No errors\n");
        res = SUCCESS;
        goto end;
    }
    */

    printf("syndrome \n");
    bmatrix_print(syndrome);

    gf2m_init(&synPoly, ctx->t);
    for(i = 0; i < syndrome->r; i++)
    {
        if(b_mat_entry(syndrome, i, 0) == 1)
        {
            synPoly.term[i/ctx->m].binary[0] ^= (1<<(i%ctx->m)); // 0 맞을까?
        }
    }
    printf("syndrome poly\n");
    gf2m_print(&synPoly);
    //matrix to polynomial

    /***********************************************
     * Step 2. Compute the inverse of syndrome
    ***********************************************/
    gf2m_init(&inv, 1);
    gf2m_init(&tmp, 1);
    gf2m_init(&gcd, 1);

    gf2m_xgcd(&gcd, &synPoly, &ctx->Goppa, &inv, &tmp, &ctx->mod);

    /***********************************************
     * Step 3. Get S(X)^-1 + X
    ***********************************************/

    gf2m_init(&X, 1);
    gf2_set_one(&X.term[1]);
    gf2m_add(&synPolyInv, &inv, &X);

    /***********************************************
     * Step 4. Get the square root T(X) of S(X)^-1 + X
    ***********************************************/

    gf2m_init(&polyTmp, ctx->t);
    gf2m_square_root(&polyTmp, &synPolyInv, &ctx->Goppa, & ctx->mod);

    /***********************************************
     * Step 5. Solve the key equation
     *  sigam(X) = a^2(X) + b^2(X)*X
     *  We use EEA
     *  b(X) * d(X) == a(X) mod g(X), where d(X) = sqrt(S(X)^-1 + X)
    ***********************************************/

    gf2m_init(&Ax, ctx->t/2);
    gf2m_init(&Bx, ctx->t/2);

    EEA_patterson(&polyTmp, &ctx->Goppa, &Ax, &Bx, &ctx->mod);

    /***********************************************
     * Step 6. Compute sigma poly. (error locator poly.)
     * sigma(X) = a(X)^2 + b(X)^2 * X
    ***********************************************/

    gf2m_init(&Ax2, 2*Ax.deg);
    gf2m_init(&Bx2, 2*Bx.deg);
    gf2m_init(&sigma, ctx->t);

    gf2m_square(&Ax2, &Ax, &ctx->mod);
    gf2m_square(&Bx2, &Bx, &ctx->mod);
    gf2m_mul(&Bx2, &Bx2, &X, &ctx->mod);
    gf2m_add(&sigma, &Ax2, &Bx2);

    /***********************************************
     * Step 7. Find all roots of sigma(X)
    ***********************************************/

   gf2_init(&rootSet, ctx->t);

end:

    return res;
}


/*
@   B*X == A mod Y
*/
int EEA_patterson(gf2m* x, gf2m* y, gf2m* a, gf2m* b, gf2m* mod)
{
    int t = y->deg;
    gf2m t0, t1, t2;
    gf2m v0, v1, v2;
    gf2m u0, u1, u2;
    gf2m R , Q;
    gf2m tmp;

    int len = (x->deg <= y->deg) ? 2 * x->deg : 2 * y->deg;

    gf2m_init(&t0, len); gf2m_init(&t1, len); gf2m_init(&t2, len);
    gf2m_init(&v0, len); gf2m_init(&v1, len); gf2m_init(&v2, len);
    gf2m_init(&u0, len); gf2m_init(&u1, len); gf2m_init(&u2, len);
    gf2m_init(&tmp, len);

    gf2m_init(&R, t);
    gf2m_init(&Q, t);

    gf2m_copy(&t0, x);         //t0 = x
    gf2m_copy(&t1, y);         //t1 = y

    gf2m_set_one(&u0);         //u0 = 1
    gf2m_set_zero(&u1);        //u1 = 0

    gf2m_set_one(&v1);         //v1 = 1
    gf2m_set_zero(&v0);        //v0 = 0

    while(t0.deg > (t/2))
    {
        gf2m_copy(&t2, &t0);              //t2 = t0
        gf2m_copy(&t0, &t1);              //t0 = t1

        gf2m_long_division(&Q, &R, &t2, &t1, mod);  // R = t2%t1
        gf2m_copy(&t1, &R);               //t1 = R

        gf2m_add(&tmp, &t2, &t1);         //t2 = t2 - t1
        gf2m_long_division(&Q, &R, &tmp, &t0, mod); // Q = t2/t0

        gf2m_copy(&u2, &u0);              //u2 = u0
        gf2m_copy(&v2, &v0);              //v2 = v0
        gf2m_copy(&u0, &u1);              //u0 = u1
        gf2m_copy(&v0, &v1);              //v0 = v1

        gf2m_mul(&tmp, &Q, &u1, mod);     //u1 = u2 - qu1
        gf2m_add(&u1, &u2, &tmp);
        gf2m_mul(&tmp, &Q, &v1, mod);     //v1 = v2 - qv1
        gf2m_add(&v1, &v2 ,&tmp);
        gf2m_set_zero(&Q);
        gf2m_set_zero(&tmp);
        
    }

    gf2m_fit_len(&u0);
    gf2m_fit_len(&v0);
    
    gf2m_copy(a, &t0);
    gf2m_copy(b, &u0);

    return SUCCESS;
}