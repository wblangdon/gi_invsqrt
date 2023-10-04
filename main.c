/*WBL 4 October 2023 "$Revision: 1.5 $" */

//WBL  4 Oct 2023 Simplify gi_cbrt/invsqrt/ r1.5

//Compile: gcc -c main.c; gcc gi_invsqrt.o main.o -lm

//Use:     ./a.out 100 (expect 0.1)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "invsqrt.h"


int main(int argc, char *argv[]){
  const double input  = (argc>1 && argv[1])? atof(argv[1]) : 1.0;
  const double output = invsqrt(input);

  printf("%g = invsqrt(%g)\n",output,input);

  return 0;
}
