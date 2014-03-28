/* Copyright (C) 1999-2012, Massachusetts Institute of Technology.
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

#ifndef MATRICES_H
#define MATRICES_H

#include "scalar.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct {
     int N, localN, Nstart, allocN;
     int c;
     int n, p, alloc_p;
     scalar *data;
} evectmatrix;

typedef struct {
     int p, alloc_p;
     scalar *data;
} sqmatrix;

/* try to keep track of flops, at least from evectmatrix multiplications */
extern double evectmatrix_flops;

/* general creation/destruction operations: */

extern evectmatrix create_evectmatrix(int N, int c, int p,
				      int localN, int Nstart, int allocN);
extern void destroy_evectmatrix(evectmatrix X);
extern sqmatrix create_sqmatrix(int p);
extern void destroy_sqmatrix(sqmatrix X);

/* diagonal matrix utils: */

extern void matrix_XtY_diag(scalar *X, scalar *Y, int n, int p, scalar *diag);
extern void matrix_XtY_diag_real(scalar *X, scalar *Y, int n, int p, 
				 real *diag);
extern void matrix_XtX_diag_real(scalar *X, int n, int p, real *diag);
extern void matrix_XpaY_diag(scalar *X, real a, scalar *Y,
			     scalar *diag, int n, int p);
extern void matrix_XpaY_diag_real(scalar *X, real a, scalar *Y,
				  real *diag, int n, int p);
extern void matrix_X_diag_real_pY_diag_real(scalar *X, real *diag1,
					    scalar *Y, real *diag2,
					    int n, int p);
extern real matrix_re_trace_A_diag_real(scalar *A, real *diag, int p);
extern scalar matrix_diag_trace(scalar *diag, int p);
extern real matrix_diag_real_trace(real *diag, int p);

/* evectmatrix operations, defined in evectmatrix.c: */

extern void evectmatrix_resize(evectmatrix *A, int p, short preserve_data);
extern void evectmatrix_copy(evectmatrix X, evectmatrix Y);
extern void evectmatrix_copy_slice(evectmatrix X, evectmatrix Y,
				   int ix, int iy, int p);
extern void evectmatrix_aXpbY(real a, evectmatrix X, real b, evectmatrix Y);
extern void evectmatrix_aXpbYS_sub(real a, evectmatrix X, 
				   real b, evectmatrix Y,
				   sqmatrix S, int Soffset, short sdagger);
extern void evectmatrix_XeYS(evectmatrix X, evectmatrix Y,
			     sqmatrix S, short sherm);
extern void evectmatrix_XpaYS(evectmatrix X, real a, evectmatrix Y,
			      sqmatrix S, short sdagger);
extern void evectmatrix_XtX(sqmatrix U, evectmatrix X, sqmatrix S);
extern void evectmatrix_XtY(sqmatrix U, evectmatrix X, evectmatrix Y,
			    sqmatrix S);
extern void evectmatrix_XtY_slice(sqmatrix U, evectmatrix X, evectmatrix Y,
				  int ix, int iy, int p, sqmatrix S);
extern void evectmatrixXtY_sub(sqmatrix U, int Uoffset,
			       evectmatrix X, evectmatrix Y, sqmatrix S);
extern void evectmatrix_XtY_diag(evectmatrix X, evectmatrix Y, scalar *diag,
				 scalar *scratch_diag);
extern void evectmatrix_XtY_diag_real(evectmatrix X, evectmatrix Y,
				      real *diag, real *scratch_diag);
extern void evectmatrix_XtX_diag_real(evectmatrix X, real *diag,
				      real *scratch_diag);
extern scalar evectmatrix_traceXtY(evectmatrix X, evectmatrix Y);

/* sqmatrix operations, defined in sqmatrix.c: */

extern void sqmatrix_assert_hermitian(sqmatrix A);
extern void sqmatrix_copy(sqmatrix A, sqmatrix B);
extern void sqmatrix_resize(sqmatrix *A, int p, short preserve_data);
extern void sqmatrix_copy_upper2full(sqmatrix F, sqmatrix U);
extern void sqmatrix_symmetrize(sqmatrix Asym, sqmatrix A);
extern scalar sqmatrix_trace(sqmatrix U);
extern scalar sqmatrix_traceAtB(sqmatrix A, sqmatrix B);
extern void sqmatrix_AeBC(sqmatrix A, sqmatrix B, short bdagger,
			  sqmatrix C, short cdagger);
extern void sqmatrix_ApaBC(sqmatrix A, real a, sqmatrix B, short bdagger,
			   sqmatrix C, short cdagger);
extern void sqmatrix_ApaB(sqmatrix A, real a, sqmatrix B);
extern void sqmatrix_aApbB(real a, sqmatrix A, real b, sqmatrix B);
extern int sqmatrix_invert(sqmatrix U, short positive_definite,
			    sqmatrix Work);
extern void sqmatrix_eigensolve(sqmatrix U, real *eigenvals, sqmatrix W);
extern void sqmatrix_sqrt(sqmatrix Usqrt, sqmatrix U, sqmatrix W);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* MATRICES_H */
