/* Copyright (C) 1999, 2000, 2001 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdlib.h>
#include <stdio.h>

#include "../config.h"
#include <check.h>

#include "matrices.h"

/* Basic operations: allocation, deallocation, etcetera. */

evectmatrix create_evectmatrix(int N, int c, int p,
			       int localN, int Nstart, int allocN)
{
     evectmatrix X;
 
     CHECK(localN <= N && allocN >= localN && Nstart < N,
	   "invalid N arguments");
    
     X.N = N;
     X.localN = localN;
     X.Nstart = Nstart;
     X.allocN = allocN;
     X.c = c;
     
     X.n = localN * c;
     X.alloc_p = X.p = p;
     
     if (allocN > 0) {
	  CHK_MALLOC(X.data, scalar, allocN * c * p);
     }
     else
	  X.data = NULL;

     return X;
}

void destroy_evectmatrix(evectmatrix X)
{
     free(X.data);
}

sqmatrix create_sqmatrix(int p)
{
     sqmatrix X;

     X.alloc_p = X.p = p;
     if (p > 0) {
	  CHK_MALLOC(X.data, scalar, p * p);
     }
     else
	  X.data = (scalar*) NULL;
     return X;
}

void destroy_sqmatrix(sqmatrix X)
{
     free(X.data);
}

/***********************************************************************/

/* a few general matrix operations for diagonal matrices; these
   will be used both by evectmatrix and sqmatrix routines: */

/* compute diag = diagonal elements of Xt * Y, where X and Y are n x p. */
void matrix_XtY_diag(scalar *X, scalar *Y, int n, int p, scalar *diag)
{
     int i, j;

     for (j = 0; j < p; ++j) {
	  ASSIGN_ZERO(diag[j]);
     }
     
     for (i = 0; i < n; ++i)
	  for (j = 0; j < p; ++j) {
	       ACCUMULATE_SUM_CONJ_MULT(diag[j], X[i*p+j], Y[i*p+j]);
	  }
}

/* compute diag = diagonal elements of Re[Xt * Y], where X and Y are n x p. */
void matrix_XtY_diag_real(scalar *X, scalar *Y, int n, int p, real *diag)
{
     int i, j;

     for (j = 0; j < p; ++j) {
	  diag[j] = 0;
     }
     
     for (i = 0; i < n; ++i)
	  for (j = 0; j < p; ++j) {
	       diag[j] += (SCALAR_RE(X[i*p+j]) * SCALAR_RE(Y[i*p+j]) + 
			   SCALAR_IM(X[i*p+j]) * SCALAR_IM(Y[i*p+j]));
	  }
}

/* compute diag = diagonal elements of Xt * X, where X is n x p. */
void matrix_XtX_diag_real(scalar *X, int n, int p, real *diag)
{
     int i, j;

     for (j = 0; j < p; ++j) {
	  diag[j] = 0;
     }
     
     for (i = 0; i < n; ++i)
	  for (j = 0; j < p; ++j) {
	       ACCUMULATE_SUM_SQ(diag[j], X[i*p+j]);
	  }
}

/* compute X += a * Y * diag(diag), where X and Y are n x p */
void matrix_XpaY_diag(scalar *X, real a, scalar *Y, 
		      scalar *diag, int n, int p)
{
     int i, j;

     for (i = 0; i < n; ++i) {
	  for (j = 0; j < p; ++j) {
	       scalar c;
	       ASSIGN_MULT(c, Y[i*p+j], diag[j]);
	       ASSIGN_SCALAR(X[i*p+j],
			     SCALAR_RE(X[i*p+j]) + a * SCALAR_RE(c),
			     SCALAR_IM(X[i*p+j]) + a * SCALAR_IM(c));
	  }
     }
}

/* compute X += a * Y * diag(diag), where X and Y are n x p and diag is real */
void matrix_XpaY_diag_real(scalar *X, real a, scalar *Y, 
			   real *diag, int n, int p)
{
     int i, j;

     for (i = 0; i < n; ++i) {
	  for (j = 0; j < p; ++j) {
	       real d = a * diag[j];
	       ASSIGN_SCALAR(X[i*p+j],
			     SCALAR_RE(X[i*p+j]) + d * SCALAR_RE(Y[i*p+j]),
			     SCALAR_IM(X[i*p+j]) + d * SCALAR_IM(Y[i*p+j]));
	  }
     }
}

/* compute X = X * diag1 + Y * diag2, where X and Y are n x p and 
   diag1 and diag2 are real diagonal matrices */
void matrix_X_diag_real_pY_diag_real(scalar *X, real *diag1,
				     scalar *Y, real *diag2, int n, int p)
{
          int i, j;

	  for (i = 0; i < n; ++i) {
	       for (j = 0; j < p; ++j) {
		    real d1 = diag1[j], d2 = diag2[j];
		    ASSIGN_SCALAR(X[i*p+j],
				  d1 * SCALAR_RE(X[i*p+j]) + 
				  d2 * SCALAR_RE(Y[i*p+j]),
				  d1 * SCALAR_IM(X[i*p+j]) +
				  d2 * SCALAR_IM(Y[i*p+j]));
	       }
	  }
}

/* compute Re [ trace A * diag(diag) ], where A is p by p. */
real matrix_re_trace_A_diag_real(scalar *A, real *diag, int p)
{
     real trace = 0.0;
     int i;
     for (i = 0; i < p; ++i)
	  trace += SCALAR_RE(A[i*(p+1)]) * diag[i];
     return trace;
}

scalar matrix_diag_trace(scalar *diag, int p)
{
     scalar trace = SCALAR_INIT_ZERO;
     int i;
     for (i = 0; i < p; ++i) {
	  ACCUMULATE_SUM(trace, diag[i]);
     }
     return trace;
}

real matrix_diag_real_trace(real *diag, int p)
{
     real trace = 0.0;
     int i;
     for (i = 0; i < p; ++i)
	  trace += diag[i];
     return trace;
}
