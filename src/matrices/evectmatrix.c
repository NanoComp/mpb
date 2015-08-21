/* Copyright (C) 1999-2014 Massachusetts Institute of Technology.
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
#include <string.h>

#include "config.h"
#include <mpiglue.h>

#include <check.h>

#include "matrices.h"
#include "blasglue.h"

double evectmatrix_flops = 0;

/* Operations on evectmatrix blocks:
       X + a Y, X * S, X + a Y * S, Xt * X, Xt * Y, trace(Xt * Y), etc.
   (X, Y: evectmatrix, S: sqmatrix) */

/* X = Y */
void evectmatrix_copy(evectmatrix X, evectmatrix Y)
{
     CHECK(X.n == Y.n && X.p == Y.p, "arrays not conformant");

     blasglue_copy(X.n * X.p, Y.data, 1, X.data, 1);
}

/* set p selected columns of X to those in Y, starting at ix and iy.  */
void evectmatrix_copy_slice(evectmatrix X, evectmatrix Y, 
			    int ix, int iy, int p)
{
     CHECK(ix + p <= X.p && iy + p <= Y.p && ix >= 0 && iy >= 0 && X.n == Y.n,
	   "invalid arguments to evectmatrix_copy_slice");

     if (ix == 0 && iy == 0 && p == X.p && p == Y.p)
	  evectmatrix_copy(X, Y);
     else if (p == 1)
	  blasglue_copy(X.n, Y.data + iy, Y.p, X.data + ix, X.p);
     else {
	  int i;
	  for (i = 0; i < X.n; ++i)
	       blasglue_copy(p, Y.data + iy + i * Y.p, 1,
			     X.data + ix + i * X.p, 1);
     }
}

/* Resize A from its current size to an nxp matrix, assuming that
   A was initially allocated to hold at least this big a matrix.
   If preserve_data is nonzero, copies the existing data in A (or
   a subset of it, if the matrix is shrinking) to the corresponding
   entries of the resized matrix. */
void evectmatrix_resize(evectmatrix *A, int p, short preserve_data)
{
     CHECK(p <= A->alloc_p, "tried to resize beyond allocated limit");

     if (preserve_data) {
	  int i, j;
	  
	  if (p < A->p) {
	       for (i = 0; i < A->n; ++i)
		    for (j = 0; j < p; ++j)
			 A->data[i*p + j] = A->data[i*A->p + j];
	  }
	  else {
	       for (i = A->n-1; i >= 0; --i)
		    for (j = A->p-1; j >= 0; --j)
			 A->data[i*p + j] = A->data[i*A->p + j];
	  }
     }

     A->p = p;
}

/* compute X = a*X + b*Y; X and Y may be equal. */
void evectmatrix_aXpbY(real a, evectmatrix X, real b, evectmatrix Y)
{
     CHECK(X.n == Y.n && X.p == Y.p, "arrays not conformant");
     
     if (a != 1.0)
	  blasglue_rscal(X.n * X.p, a, X.data, 1);

     blasglue_axpy(X.n * X.p, b, Y.data, 1, X.data, 1);
     evectmatrix_flops += X.N * X.c * X.p * 3;
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
	  evectmatrix_flops += X.N * X.c * X.p * (3 + 2 * X.p);
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

/* compute U = adjoint(X) * X, with S a scratch matrix. */
void evectmatrix_XtX(sqmatrix U, evectmatrix X, sqmatrix S)
{
     CHECK(X.p == U.p && U.p <= S.alloc_p, "matrices not conformant");
     
/*
     blasglue_gemm('C', 'N', X.p, X.p, X.n,
		   1.0, X.data, X.p, X.data, X.p, 0.0, S.data, U.p);
*/

     /* take advantage of the fact that U is Hermitian and only write
	out the upper triangle of the matrix */
     memset(S.data, 0, sizeof(scalar) * (U.p * U.p));
     blasglue_herk('U', 'C', X.p, X.n, 1.0, X.data, X.p, 0.0, S.data, U.p);
     evectmatrix_flops += X.N * X.c * X.p * (X.p - 1);

     /* Now, copy the conjugate of the upper half onto the lower half of S */
     {
	  int i, j;

	  for (i = 0; i < U.p; ++i)
	       for (j = i + 1; j < U.p; ++j) {
		    ASSIGN_CONJ(S.data[j * U.p + i], S.data[i * U.p + j]);
	       }
     }

     mpi_allreduce(S.data, U.data, U.p * U.p * SCALAR_NUMVALS,
		   real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
}

/* Dot p selected columns of X with those in Y, starting at ix and iy.
   Stores the result in U, with S a scratch matrix. */
void evectmatrix_XtY_slice(sqmatrix U, evectmatrix X, evectmatrix Y,
			   int ix, int iy, int p, sqmatrix S)
{
     CHECK(ix + p <= X.p && iy + p <= Y.p && ix >= 0 && iy >= 0 && X.n == Y.n
           && p == U.p && p <= S.alloc_p, "invalid arguments to XtY_slice");

     memset(S.data, 0, sizeof(scalar) * (U.p * U.p));
     blasglue_gemm('C', 'N', p, p, X.n,
                   1.0, X.data + ix, X.p, Y.data + iy, Y.p, 0.0, S.data, U.p);
     evectmatrix_flops += X.N * X.c * p * (2*p);

     mpi_allreduce(S.data, U.data, U.p * U.p * SCALAR_NUMVALS,
                   real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
}

/* compute U = adjoint(X) * Y, with S a scratch matrix. */
void evectmatrix_XtY(sqmatrix U, evectmatrix X, evectmatrix Y, sqmatrix S)
{
     CHECK(X.p == Y.p, "matrices not conformant");
     
     evectmatrix_XtY_slice(U, X, Y, 0, 0, X.p, S);
}

/* Compute adjoint(X) * Y, storing the result in U at an offset
   Uoffset with the matrix (i.e. as a submatrix within U).  S is a
   scratch matrix (at least Y.p by Y.p). */
void evectmatrixXtY_sub(sqmatrix U, int Uoffset, evectmatrix X, evectmatrix Y,
			sqmatrix S)
{
     int i;

     CHECK(X.p == Y.p && X.n == Y.n && U.p >= Y.p, "matrices not conformant");
     CHECK(Uoffset + (Y.p-1)*U.p + Y.p <= U.p*U.p,
	   "submatrix exceeds matrix bounds");
     CHECK(Y.p <= S.alloc_p, "scratch matrix too small");
     
     memset(S.data, 0, sizeof(scalar) * (Y.p * Y.p));
     blasglue_gemm('C', 'N', X.p, X.p, X.n,
		   1.0, X.data, X.p, Y.data, Y.p, 0.0, S.data, Y.p);
     evectmatrix_flops += X.N * X.c * X.p * (2*X.p);

     for (i = 0; i < Y.p; ++i) {
	  mpi_allreduce(S.data + i*Y.p, U.data + Uoffset + i*U.p, 
			Y.p * SCALAR_NUMVALS,
			real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
     }
}

/* Compute only the diagonal elements of XtY, storing in diag
   (with scratch_diag a scratch array of the same size as diag). */
void evectmatrix_XtY_diag(evectmatrix X, evectmatrix Y, scalar *diag,
			  scalar *scratch_diag)
{
     matrix_XtY_diag(X.data, Y.data, X.n, X.p, scratch_diag);
     evectmatrix_flops += X.N * X.c * X.p * 2;
     mpi_allreduce(scratch_diag, diag, X.p * SCALAR_NUMVALS, 
		   real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
}

/* As above, but only compute real parts of diagonal. */
void evectmatrix_XtY_diag_real(evectmatrix X, evectmatrix Y, real *diag,
			       real *scratch_diag)
{
     matrix_XtY_diag_real(X.data, Y.data, X.n, X.p, scratch_diag);
     evectmatrix_flops += X.N * X.c * X.p * (2*X.p);
     mpi_allreduce(scratch_diag, diag, X.p,
		   real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
}

/* As above, but compute only the diagonal elements of XtX. */
void evectmatrix_XtX_diag_real(evectmatrix X, real *diag, real *scratch_diag)
{
     matrix_XtX_diag_real(X.data, X.n, X.p, scratch_diag);
     evectmatrix_flops += X.N * X.c * X.p * (2*X.p);
     mpi_allreduce(scratch_diag, diag, X.p,
		   real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
}

/* compute trace(adjoint(X) * Y) */
scalar evectmatrix_traceXtY(evectmatrix X, evectmatrix Y)
{
     scalar trace, trace_scratch;

     CHECK(X.p == Y.p && X.n == Y.n, "matrices not conformant");
     
     trace_scratch = blasglue_dotc(X.n * X.p, X.data, 1, Y.data, 1);
     evectmatrix_flops += X.N * X.c * X.p * (2*X.p) + X.p;

     mpi_allreduce(&trace_scratch, &trace, SCALAR_NUMVALS,
		   real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);

     return trace;
}
