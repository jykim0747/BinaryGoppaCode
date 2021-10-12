#ifndef _MCELIECE_H_
#define _MCELIECE_H_

#include "bmatrix.h"
#include "gf2m.h"

typedef struct
{
    BMAT S; //invertible matrix
    BMAT P; //permutation matrix
    BMAT SGP; // public key
} Key;

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
    Key key;
} Param;

void Fisher_Yate(int* set, int len);
void clearParam(Param* ctx);
int generateSupportSet(Param* ctx);
int get_paritycheck_matrix(Param* ctx);
int get_generator_matrix(Param* ctx);
int validate(Param* ctx);
int get_generate_key(Param* ctx);
void generateError(BMAT error, int num);

void encryption(BMAT dst, BMAT src, Param* ctx);
void decryption(BMAT dst, BMAT src, Param* ctx);
int patterson_decoding(BMAT dst, BMAT src, Param* ctx);
int EEA_patterson(gf2m* x, gf2m* y, gf2m* a, gf2m* b, gf2m* mod);

gf2 horner_method(gf2m* poly, gf2* src, gf2* mod);
int find_root(gf2* root_set, gf2m* poly, Param* ctx);

void test_mceliece_operation();

#endif