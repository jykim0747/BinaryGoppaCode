#include "gf2.h"
#include "gf2m.h"
#include "bmatrix.h"
#include "gf2_matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){

    srand((unsigned int)time(NULL));
    //test_gf2_init();
    test_gf2_math_operation();
    //test_gf2m_math_operation();

    //test_gf2_matrix_operation();
    //test_bmatrix_operation();

    return 1;
}