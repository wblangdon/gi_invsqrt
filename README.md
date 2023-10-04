# gi_invsqrt

W.Langdon cs.ucl.ac.uk 26 June 2018 based on opencv_ssbse2016/gp/README.txt r1.9 $Revision: 1.10 $

reciprocal square root

Assumes use of Linux tcsh scripts and a recent GNU gcc compiler.

C Code derived from

W. B. Langdon. Genetic Improvement of Data gives double precision invsqrt.
In  7th edition of GI @ GECCO Companion 2019, pages 1709-1714, Prague, Czech Republic, 2019. ACM 
http://www.cs.ucl.ac.uk/staff/W.Langdon/ftp/papers/langdon_2019_GI7.pdf
doi:10.1145/3319619.3326800

The evolved code is accurate to double precision

GGGP project at http://www.cs.ucl.ac.uk/staff/W.Langdon/gggp

We do not include copies of the GNU C library since 
1) you may already have it installed and 2) it is large.
We downloaded glibc-2.27.tar.gz from www.gnu.org/s/libc
Similarly the C version of CMA-ES may be found via
github.com/cma-es/c-cmaes/archive/master.zip
documentation at http://cma.gforge.inria.fr/cmaes_sourcecode_page.html#C

Based on gi_cbrt.tar.gz README.txt 12 March 2020
gi_cbrt.tar.gz conatins results and code for several functions:
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

