# gi_invsqrt

reciprocal square root

Assumes use of Linux tcsh scripts and a recent GNU gcc compiler.

C Code derived from:

  W. B. Langdon. Genetic Improvement of Data gives double precision invsqrt.
  In  7th edition of GI @ GECCO Companion 2019, pages 1709-1714, Prague, Czech Republic, 2019. ACM 
  http://www.cs.ucl.ac.uk/staff/W.Langdon/ftp/papers/langdon_2019_GI7.pdf
  https://doi.org/10.1145/3319619.3326800

The evolved code is accurate to double precision

GGGP project at http://www.cs.ucl.ac.uk/staff/W.Langdon/gggp

Based on gi_cbrt.tar.gz README.txt 12 March 2020
gi_cbrt.tar.gz contains full results and code for several functions:
only invsqrt is included here


Compiling
~~~~~~~~~

gcc -c invsqrt.c

to compile example:

gcc -c main.c; gcc invsqrt.o main.o -lm


Testing
~~~~~~~
./a.out 4
0.5 = invsqrt(4)

