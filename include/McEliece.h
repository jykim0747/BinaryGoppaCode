#ifndef _MCELIECE_H_
#define _MCELIECE_H_

#include "bmatrix.h"
#include "gf2m.h"

typedef struct
{
    int n;
    int m;
    int t;
    gf2 mod;
    gf2m Goppa;
    int* supportSet;
    BMAT paritycheckMatrix;
    BMAT generatorMatrix;

} Param;

void Fisher_Yate(int* set, int len);
void clearParam(Param* ctx);
int get_paritycheck_matrix(Param* ctx);
int get_generator_matrix(Param* ctx);
int validate(Param* ctx);

void test_mceliece_operation();

#endif