/* WBL 3 March 2019 $Revision: 1.22 $
based on gi_cbrt.c r1.20 and e_sqrt.c r1.9
glibc-2.29/sysdeps/powerpc/fpu/e_sqrt.c as 2.27 except #include <fenv.h>
*/
/* Double-precision floating point square root.
   Copyright (C) 1997-2018 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */

/*Modifications:
WBL  4 Oct 2023 Rename invsqrt
           Strip away evolutionary framework to leave only actual invsqrt
WBL  9 Mar 2019 Back to r1.9 make sure loglike returns 0 only if diff is 0
WBL  6 Mar 2019 Remove much debug before trying CMA-ES
WBL  4 Mar 2019 Convert for use as GI version of inverse square root
           based on glibc-2.29/sysdeps/powerpc/fpu/e_sqrt.c with e_sqrt.c r1.9
	   for invsqrt no need to divide in Newton-Raphson so can avoid sy
WBL 26 Jun 2018 For release kit, consolidate different versions
WBL 21 May 2018 Reverted back to r1.10 for testing mode,
           add read_double,print_double,tweak_double
WBL 18 May 2018 hack from sqrt to cbrt
 */

#ifndef NDEBUG
#include <stdio.h>
#include <assert.h>
#include <float.h>
#endif

#include <math.h>
#include <fenv.h>
//#include <fenv_libc.h> powerpc specific
#include <inttypes.h>
#include <stdint.h>
//#include <sysdep.h>
//#include <ldsodefs.h>
#include "math_private.h"

#include "invsqrt.h"
#include "t_invsqrt.c"

#ifndef _ARCH_PPCSQ

//#define __builtin_fma(a, b, c) (((a)*(b))+(c))
//__builtin_fma really just a placeholder so easy to see code changes
#define __builtin_fma(a, c) ((a)+(c))

/* from sysdeps/powerpc/fpu/fenv_libc.h, ie powerpc specific
This operation (i) sets the appropriate FPSCR bits for its
   parameter, (ii) converts sNaN to the corresponding qNaN, and (iii)
   otherwise passes its parameter through unchanged (in particular, -0
   and +0 stay as they were).  The `obvious' way to do this is optimised
   out by gcc.
#define f_wash(x) \
   ({ double d; asm volatile ("fmul %0,%1,%2" \
			      : "=f"(d) \
			      : "f" (x), "f"((float)1.0)); d; })
*/
#define f_wash(x) (x)

static const ieee_float_shape_type a_nan = {.word = 0x7fc00000 };
static const ieee_float_shape_type a_inf = {.word = 0x7f800000 };
static const float two108 = 3.245185536584267269e+32;
static const float two54 = 18014398509481984.0;

/* treat comments with caution many are from powerpc sqrt

   The method is based on a description in
   Computation of elementary functions on the IBM RISC System/6000 processor,
   P. W. Markstein, IBM J. Res. Develop, 34(1) 1990.
   Basically, it consists of two interleaved Newton-Raphson approximations,
   one to find the actual square root, and one to find its reciprocal
   without the expense of a division operation.   The tricky bit here
   is the use of the POWER/PowerPC multiply-add operation to get the
   required accuracy with high speed.

   The argument reduction works by a combination of table lookup to
   obtain the initial guesses, and some careful modification of the
   generated guesses (which mostly runs on the integer unit, while the
   Newton-Raphson is running on the FPU).  */

//WBL 3 March 2019 Replace table driven cbrt by table driven invsqrt based on e_sqrt.c's __slow_ieee754_sqrt
//double table_ieee754_invsqrt (const double x)
double invsqrt(const double x)
{

  const float inf = a_inf.value;

  if (x > 0)
    {
      /* schedule the EXTRACT_WORDS to get separation between the store
	 and the load.  */
      ieee_double_shape_type ew_u;
      ieee_double_shape_type iw_u;
      ew_u.value = (x);
      if (x != inf)
	{
	  /* Variables named starting with 's' exist in the
	     argument-reduced space, so that 2 > sx >= 0.5,
	     1.41... > sg >= 0.70.., 0.70.. >= sy > 0.35... .
	     Variables named ending with 'i' are integer versions of
	     floating-point values.  */
	  double sx;	/* The value of which we're trying to find the
			   square root.  */
	  double sg, g;	/* Guess of the square root of x.  */
	  double sd;/*d; * Difference between the square of the guess and x.  */
	  //double sy;	/* Estimate of 1/2g (overestimated by 1ulp).  */
	  //double sy2;	/* 2*sy */
	  //double e;	/* Difference between y*g and 1/2 (se = e * fsy).  */
	  double shx;	/* == sx * fsg */
	  double fsg;	/* sg*fsg == g.  */
	  /* fenv_t fe;	   Saved floating-point environment (stores rounding
			   mode and whether the inexact exception is
			   enabled).  */
	  uint32_t xi0, xi1, sxi, fsgi;
	  //const float *t_sqrt;

	  //fe = fegetenv_register ();
	  /* complete the EXTRACT_WORDS (xi0,xi1,x) operation.  */
	  xi0 = ew_u.parts.msw;
	  xi1 = ew_u.parts.lsw;
	  //relax_fenv_state ();
	  sxi = (xi0 & 0x3fffffff) | 0x3fe00000;
	  /* schedule the INSERT_WORDS (sx, sxi, xi1) to get separation
	     between the store and the load.  */
	  iw_u.parts.msw = sxi;
	  iw_u.parts.lsw = xi1;
	  const int idx = (xi0 >> (52 - 32 - 8 - 1) & 0x3fe)/2;
	  assert(idx >= 0);
	  assert(idx < 512);
	  sg = __t_invsqrt[idx];
	  //sy = t_sqrt[1];
	  //printf("Table values %g %g\n",sg);
	  //assert(sg==sg_);
	  /* complete the INSERT_WORDS (sx, sxi, xi1) operation.  */
	  sx = iw_u.value;

	  /* Here we have three Newton-Raphson iterations each of a
	     division and a square root and the remainder of the
	     argument reduction, all interleaved.   */
	  sd = __builtin_fma (sg, -sx*sg*sg*sg)/2;
	  //printf("Normalised target %g initial error %g\n",sx,sd);
	  const uint32_t tmp = (xi0 + 0x40000000) >> 1 & 0x7ff00000;
	  fsgi = 0x7fe00000 - tmp; /*invert exponent*/

	  //printf("xi0 0x%08x, sxi 0x%08x, __t_cbrt[%3d] %g, \n",
	  //	         xi0,        sxi,          idx, sg);
	  //printf("tmp 0x%08x, fsgi 0x%08x\n",
	  //	  tmp,        fsgi);

	  //sy2 = sy + sy;
	  sg = __builtin_fma (sd, sg);	/* 16-bit approximation to
						   sqrt(sx). */

	  /* schedule the INSERT_WORDS (fsg, fsgi, 0) to get separation
	     between the store and the load.  */
	  INSERT_WORDS (fsg, fsgi, 0);
	  iw_u.parts.msw = fsgi;
	  iw_u.parts.lsw = (0);
	  //e = -__builtin_fma (sy, sg, -almost_half);
	  sd = __builtin_fma (sg, -sx*sg*sg*sg)/2;
	  //printf("Second estimate %g 2nd error %g\n",sg,sd);
	  if ((xi0 & 0x7ff00000) == 0)
	    goto denorm;
	  //sy = __builtin_fma (e, sy2, sy);
	  sg = __builtin_fma (sd, sg);	/* 32-bit approximation to
						   sqrt(sx).  */
	  //sy2 = sy + sy;
	  /* complete the INSERT_WORDS (fsg, fsgi, 0) operation.  */
	  fsg = iw_u.value;
	  //e = -__builtin_fma (sy, sg, -almost_half);
	  sd = __builtin_fma (sg, -sx*sg*sg*sg)/2;
	  //printf("3rd estimate %g 3rd error %g\n",sg,sd);
	  //sy = __builtin_fma (e, sy2, sy);
	  shx = sx * fsg;
	  sg = __builtin_fma (sd, sg);	/* 64-bit approximation to
						   sqrt(sx), but perhaps
						   rounded incorrectly.  */
	  //printf("fsg %g last normalised estimate %g shx %g\n",fsg,sg,shx);

	  //sy2 = sy + sy;
	  g = sg * fsg;
	  //e = -__builtin_fma (sy, sg, -almost_half);
	  //do we need this? d = __builtin_fma (g, -shx);
	  //sy = __builtin_fma (e, sy2, sy);
	  //printf("shx %g last estimate %g error %g\n",shx,g,d);
	  //fesetenv_register (fe);
	  const double ans = g;//__builtin_fma (d, g);
	  //printf("table_ieee754_invsqrt(%g) returns %g\n",x,ans);
	  return ans;
	denorm:
	  //assert(0); //we are not seeing this...
	  /* For denormalised numbers, we normalise, calculate the
	     square root, and return an adjusted result.  */
	  //fesetenv_register (fe);
	  return invsqrt(x * two108) * two54;
	}
    }
  else if (x < 0)
    {
      /* For some reason, some PowerPC32 processors don't implement
	 FE_INVALID_SQRT.  */
#ifdef FE_INVALID_SQRT
      __feraiseexcept (FE_INVALID_SQRT);

      fenv_union_t u = { .fenv = fegetenv_register () };
      if ((u.l & FE_INVALID) == 0)
#endif
	//__feraiseexcept (FE_INVALID);
	feraiseexcept (FE_INVALID);
      return a_nan.value;
    }
  return f_wash (x);
}
#endif /* _ARCH_PPCSQ  */
