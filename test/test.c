#include "gf2.h"
#include "gf2m.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){

    srand((unsigned int)time(NULL));
    //test_gf2_init();
    //test_gf2_math_operation();
  
    test_gf2m_math_operation();

    return 1;
}