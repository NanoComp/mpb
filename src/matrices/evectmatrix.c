/* Copyright (C) 1999, 2000 Massachusetts Institute of Technology.
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

#include <stdio.h>
#include <stdlib.h>

#include "../config.h"
#include <mpiglue.h>

#include <check.h>

#include "matrices.h"
#include "blasglue.h"

/* Operations on evectmatrix blocks:
       X + a Y, X * S, X + a Y * S, Xt * X, Xt * Y, trace(Xt * Y), etc.
   (X, Y: evectmatrix, S: sqmatrix) */

/* X = Y */
void evectmatrix_copy(evectmatrix X, evectmatrix Y)
{
     CHECK(X.n == Y.n && X.p == Y.p, "arrays not conformant");

     blasglue_copy(X.n * X.p, Y.data, 1, X.data, 1);
}

/* compute X = a*X + b*Y; X and Y may be equal. */
void evectmatrix_aXpbY(real a, evectmatrix X, real b, evectmatrix Y)
{
     CHECK(X.n == Y.n && X.p == Y.p, "arrays not conformant");
     
     if (a != 1.0)
	  blasglue_scal(X.n * X.p, a, X.data, 1);

     blasglue_axpy(X.n * X.p, b, Y.data, 1, X.data, 1);
}

/* Compute X = a*X + b*Y*S.  Instead of using the entire S matrix, however,
   we use only a Y.p x Y.p submatrix, beginning at the element indexed
   by Soffset.  If sdagger != 0, then the adjoint of the submatrix is
   used instead of the submatrix. */
void evectmatrix_aXpbYS_sub(real a, evectmatrix X, real b, evectmatrix Y,
			    sqmatrix S, int Soffset, short sdagger)
{
     if (S.p == 0)  /* we treat the S.p == 0 case as if S were the identity */
	  evectmatrix_aXpbY(a, X, b, Y);
     else {
	  CHECK(X.n == Y.n && X.p == Y.p && X.p <= S.p,
		"arrays not conformant");
	  CHECK(Soffset + (Y.p-1)*S.p + Y.p <= S.p*S.p,
		"submatrix exceeds matrix bounds");
	  blasglue_gemm('N', sdagger ? 'C' : 'N', X.n, X.p, X.p,
			b, Y.data, Y.p, S.data + Soffset, S.p,
			a, X.data, X.p);
     }
}

/* compute X = YS.  If sherm != 0, then S is assumed to be Hermitian.
   This can be used to make the multiplication more efficient. */
void evectmatrix_XeYS(evectmatrix X, evectmatrix Y, sqmatrix S, short sherm)
{
     CHECK(S.p == 0 || S.p == Y.p, "arrays not conformant");
     evectmatrix_aXpbYS_sub(0.0, X, 1.0, Y, S, 0, sherm);
}

/* compute X += a Y * S.  If sdagger != 0, then St is used instead of S. */
void evectmatrix_XpaYS(evectmatrix X, real a, evectmatrix Y,
		       sqmatrix S, short sdagger)
{
     CHECK(S.p == 0 || S.p == Y.p, "arrays not conformant");
     evectmatrix_aXpbYS_sub(1.0, X, a, Y, S, 0, sdagger);
}

/* compute U = adjoint(X) * X */
void evectmatrix_XtX(sqmatrix U, evectmatrix X)
{
     CHECK(X.p == U.p, "matrices not conformant");
     
/*
     blasglue_gemm('C', 'N', X.p, X.p, X.n,
		   1.0, X.data, X.p, X.data, X.p, 0.0, U.data, U.p);
*/

     /* take advantage of the fact that U is Hermitian and only write
	out the upper triangle of the matrix */
     blasglue_herk('U', 'C', X.p, X.n, 1.0, X.data, X.p, 0.0, U.data, U.p);

     /* Now, copy the conjugate of the upper half onto the lower half of U */
     {
	  int i, j;

	  for (i = 0; i < U.p; ++i)
	       for (j = i + 1; j < U.p; ++j) {
		    ASSIGN_CONJ(U.data[j * U.p + i], U.data[i * U.p + j]);
	       }
     }

     MPI_Allreduce(U.data, U.data, U.p * U.p * SCALAR_NUMVALS,
		   SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
}

/* compute U = adjoint(X) * Y */
void evectmatrix_XtY(sqmatrix U, evectmatrix X, evectmatrix Y)
{
     CHECK(X.p == Y.p && X.n == Y.n && X.p == U.p, "matrices not conformant");
     
     blasglue_gemm('C', 'N', X.p, X.p, X.n,
		   1.0, X.data, X.p, Y.data, Y.p, 0.0, U.data, U.p);

     MPI_Allreduce(U.data, U.data, U.p * U.p * SCALAR_NUMVALS,
		   SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
}

/* Compute adjoint(X) * Y, storing the result in U at an offset Uoffset
   with the matrix (i.e. as a submatrix within U). */
void evectmatrixXtY_sub(sqmatrix U, int Uoffset, evectmatrix X, evectmatrix Y)
{
     int i;

     CHECK(X.p == Y.p && X.n == Y.n && U.p >= Y.p, "matrices not conformant");
     CHECK(Uoffset + (Y.p-1)*U.p + Y.p <= U.p*U.p,
	   "submatrix exceeds matrix bounds");

     blasglue_gemm('C', 'N', X.p, X.p, X.n,
                   1.0, X.data, X.p, Y.data, Y.p, 0.0, U.data + Uoffset, U.p);

     for (i = 0; i < Y.p; ++i) {
	  MPI_Allreduce(U.data + Uoffset + i*U.p, U.data + Uoffset + i*U.p,
			Y.p * SCALAR_NUMVALS,
			SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
     }
}

/* compute only the diagonal elements of XtY */
void evectmatrix_XtY_diag(evectmatrix X, evectmatrix Y, scalar *diag)
{
     matrix_XtY_diag(X.data, Y.data, X.n, X.p, diag);
     MPI_Allreduce(diag, diag, X.p, SCALAR_MPI_TYPE * SCALAR_NUMVALS,
		   MPI_SUM, MPI_COMM_WORLD);
}

void evectmatrix_XtY_diag_real(evectmatrix X, evectmatrix Y, real *diag)
{
     matrix_XtY_diag_real(X.data, Y.data, X.n, X.p, diag);
     MPI_Allreduce(diag, diag, X.p, SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
}

/* compute only the diagonal elements of XtX */
void evectmatrix_XtX_diag_real(evectmatrix X, real *diag)
{
     matrix_XtX_diag_real(X.data, X.n, X.p, diag);
     MPI_Allreduce(diag, diag, X.p, SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
}

/* compute trace(adjoint(X) * Y) */
scalar evectmatrix_traceXtY(evectmatrix X, evectmatrix Y)
{
     scalar trace;

     CHECK(X.p == Y.p && X.n == Y.n, "matrices not conformant");
     
     trace = blasglue_dotc(X.n * X.p, X.data, 1, Y.data, 1);

     MPI_Allreduce(&trace, &trace, SCALAR_NUMVALS, SCALAR_MPI_TYPE,
		   MPI_SUM, MPI_COMM_WORLD);

     return trace;
}
