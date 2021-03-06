#include "McEliece.h"

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