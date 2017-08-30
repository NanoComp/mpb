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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "config.h"
#include <check.h>

#include "matrices.h"
#include "blasglue.h"

/* Simple operations on sqmatrices.  Also, some not-so-simple operations,
   like inversion and eigenvalue decomposition. */

#ifdef DEBUG
static double max2(double a, double b) { return a > b ? a : b; }
#endif

void sqmatrix_assert_hermitian(sqmatrix A)
{
#ifdef DEBUG
     double err = 0, maxsq = 0;
     int i, j, p = A.p;

     for (i = 0; i < p; ++i)
	  for (j = 0; j < p; ++j)
	       maxsq = max2(maxsq, SCALAR_NORMSQR(A.data[i*p + j]));
     for (i = 0; i < p; ++i) {
	  err = max2(err, (SCALAR_IM(A.data[i*p + i]) *
			   SCALAR_IM(A.data[i*p + i])) / maxsq);
	  for (j = i + 1; j < p; ++j) {
	       scalar x;
	       ASSIGN_CONJ(x, A.data[i*p + j]);
	       ACCUMULATE_DIFF(x, A.data[j*p + i]);
	       err = max2(err, SCALAR_NORMSQR(x) / maxsq);
	  }
     }
     CHECK(err < 1e-10, "sqmatrix_assert_hermitian failed");
#else
     (void) A; /* unused */
#endif
}

/* A = B */
void sqmatrix_copy(sqmatrix A, sqmatrix B)
{
     CHECK("A.p == B.p", "arrays not conformant"); 
     blasglue_copy(A.p * A.p, B.data, 1, A.data, 1);
}

/* Resize A from its current size to a pxp matrix, assuming that
   A was initially allocated to hold at least this big a matrix.
   If preserve_data is nonzero, copies the existing data in A (or
   a subset of it, if the matrix is shrinking) to the corresponding
   entries of the resized matrix. */
void sqmatrix_resize(sqmatrix *A, int p, short preserve_data)
{
     CHECK(p <= A->alloc_p, "tried to resize beyond allocated limit");

     if (preserve_data) {
	  int i, j;

	  if (p < A->p) {
	       for (i = 0; i < p; ++i)
		    for (j = 0; j < p; ++j)
			 A->data[i*p + j] = A->data[i*A->p + j];
	  }
	  else {
	       for (i = A->p-1; i >= 0; --i)
		    for (j = A->p-1; j >= 0; --j)
			 A->data[i*p + j] = A->data[i*A->p + j];
	  }
     }

     A->p = p;
}

/* U contains the upper triangle of a Hermitian matrix; we copy this
   to F and also fill in the lower triangle with the adjoint of the upper. */
void sqmatrix_copy_upper2full(sqmatrix F, sqmatrix U)
{
     int i, j;

     CHECK(F.p == U.p, "arrays not conformant");
     for (i = 0; i < U.p; ++i) {
	  for (j = 0; j < i; ++j) {
	       ASSIGN_CONJ(F.data[i*U.p + j], U.data[j*U.p + i]);
	  }
	  for (; j < U.p; ++j)
	       F.data[i*U.p + j] = U.data[i*U.p + j];
     }
     sqmatrix_assert_hermitian(F);
}

/* Asym = (A + adjoint(A)) / 2.  Asym is thus Hermitian. */
void sqmatrix_symmetrize(sqmatrix Asym, sqmatrix A)
{
     int i, j;

     CHECK(Asym.p == A.p, "arrays not conformant");

     for (i = 0; i < A.p; ++i)
	  for (j = 0; j < A.p; ++j) {
	       int ij = i * A.p + j, ji = j * A.p + i;
	       ASSIGN_SCALAR(Asym.data[ij],
			     0.5 * (SCALAR_RE(A.data[ij]) +
				    SCALAR_RE(A.data[ji])),
			     0.5 * (SCALAR_IM(A.data[ij]) -
				    SCALAR_IM(A.data[ji])));
	  }
     sqmatrix_assert_hermitian(Asym);
}

/* trace(U) */
scalar sqmatrix_trace(sqmatrix U)
{
     int i;
     scalar trace = SCALAR_INIT_ZERO;

     for (i = 0; i < U.p; ++i)
	  ACCUMULATE_SUM(trace, U.data[i*U.p + i]);

     return trace;
}

/* compute trace(adjoint(A) * B) */
scalar sqmatrix_traceAtB(sqmatrix A, sqmatrix B)
{
     scalar trace;

     CHECK(A.p == B.p, "matrices not conformant");

     trace = blasglue_dotc(A.p * A.p, A.data, 1, B.data, 1);

     return trace;
}

/* A = B * C.  If bdagger != 0, then adjoint(B) is used; similarly for C.
   A must be distinct from B and C.   Note that since the matrices
   are stored in row-major order, the most efficient operation should
   be B * adjoint(C), assuming the BLAS is sane.  i.e. if C is hermitian,
   you should use cdagger = 1.  Conversely, the worst operation is
   probably adjoint(B) * C. */
void sqmatrix_AeBC(sqmatrix A, sqmatrix B, short bdagger,
		   sqmatrix C, short cdagger)
{
     CHECK(A.p == B.p && A.p == C.p, "matrices not conformant");

     blasglue_gemm(bdagger ? 'C' : 'N', cdagger ? 'C' : 'N', A.p, A.p, A.p,
                   1.0, B.data, B.p, C.data, C.p, 0.0, A.data, A.p);
}

/* A += a B * C.  bdagger, cdagger are as for sqmatrix_AeBC, above. */
void sqmatrix_ApaBC(sqmatrix A, real a, sqmatrix B, short bdagger,
		    sqmatrix C, short cdagger)
{
     CHECK(A.p == B.p && A.p == C.p, "matrices not conformant");

     blasglue_gemm(bdagger ? 'C' : 'N', cdagger ? 'C' : 'N', A.p, A.p, A.p,
                   a, B.data, B.p, C.data, C.p, 1.0, A.data, A.p);
}

/* A += a B */
void sqmatrix_ApaB(sqmatrix A, real a, sqmatrix B)
{
     CHECK(A.p == B.p, "matrices not conformant");

     blasglue_axpy(A.p * A.p, a, B.data, 1, A.data, 1);
}

/* compute A = a*A + b*B; A and B may be equal. */
void sqmatrix_aApbB(real a, sqmatrix A, real b, sqmatrix B)
{
     CHECK(A.p == B.p, "arrays not conformant");

     if (a != 1.0)
          blasglue_rscal(A.p * A.p, a, A.data, 1);

     blasglue_axpy(A.p * A.p, b, B.data, 1, A.data, 1);
}

/* U <- 1/U.  U must be Hermitian and, if positive_definite != 0,
   positive-definite (e.g. U = Yt*Y).  Work is a scratch matrix.
   Returns 1 on success, 0 if failure (e.g. matrix singular) */
int sqmatrix_invert(sqmatrix U, short positive_definite,
		     sqmatrix Work)
{
     int i, j;

     sqmatrix_assert_hermitian(U);
     if (positive_definite) {
	  /* factorize U: */
	  if (!lapackglue_potrf('U', U.p, U.data, U.p)) return 0;

	  /* QUESTION: would it be more efficient to stop here,
	     returning the Cholesky factorization of U?  This
	     could then be used to multiply by 1/U without
	     ever calculating the inverse explicitly.  It
	     would probably be more numerically stable, but
	     how do the computational costs compare? */

	  /* Compute 1/U (upper half only) */
	  if (!lapackglue_potri('U', U.p, U.data, U.p)) return 0;
     }
     else {
	  int *ipiv;
	  CHK_MALLOC(ipiv, int, U.p);

	  CHECK(Work.p * Work.p >= U.p, "scratch matrix is too small");

	  /* factorize U: */
	  if (!lapackglue_hetrf('U', U.p, U.data, U.p,
				ipiv, Work.data, Work.p * Work.p)) return 0;
	  /* Compute 1/U (upper half only) */
	  if (!lapackglue_hetri('U', U.p, U.data, U.p, ipiv, Work.data))
	       return 0;

	  free(ipiv);
     }

     /* Now, copy the conjugate of the upper half
	onto the lower half of U */
     for (i = 0; i < U.p; ++i)
	  for (j = i + 1; j < U.p; ++j) {
	       ASSIGN_CONJ(U.data[j * U.p + i], U.data[i * U.p + j]);
	  }

     return 1;
}

/* U <- eigenvectors of Ux=lambda B x, while B is overwritten (by its
   Cholesky factors).  U and B must be Hermitian, and B must be
   positive-definite; if B==NULL then it is taken to be the
   identity. eigenvals <- eigenvalues.  W is a work array.  The
   columns of adjoint(U') are the eigenvectors, so that U ==
   adjoint(U') D U', with D = diag(eigenvals).

   The eigenvalues are returned in ascending order. */
void sqmatrix_gen_eigensolve(sqmatrix U, sqmatrix B, real *eigenvals, sqmatrix W)
{
     real *work;
     scalar *morework;
     int nwork;

     sqmatrix_assert_hermitian(U);
     CHK_MALLOC(work, real, 3*U.p - 2);
     if (W.p * W.p >= 3 * U.p - 1) {
         morework = W.data;
         nwork = W.p * W.p;
     }
     else {
         CHK_MALLOC(morework, scalar, 3 * U.p - 1);
         nwork = 3 * U.p - 1;
     }
     if (B.data) {
         CHECK(B.p == U.p, "mismatched matrix sizes in sqmatrix_eigensolve");
         sqmatrix_assert_hermitian(B);
         lapackglue_hegv(1, 'V', 'U', U.p, U.data, B.p, B.data, U.p, eigenvals,
                         morework, nwork, work);
     }
     else {
         lapackglue_heev('V', 'U', U.p, U.data, U.p, eigenvals,
                         morework, nwork, work);
     }

     if (morework != W.data) free(morework);
     free(work);
}

void sqmatrix_eigensolve(sqmatrix U, real *eigenvals, sqmatrix W)
{
    sqmatrix B;
    B.data = NULL;
    sqmatrix_gen_eigensolve(U, B, eigenvals, W);
}

/* Compute eigenvalues of a general (non-Hermitian) matrix A.
   Does not compute the eigenvectors. */
void sqmatrix_eigenvalues(sqmatrix A, scalar_complex *eigvals)
{
    sqmatrix B; /* make a copy of A, since geev overwrites array */
    scalar *work, work1;
    real *rwork;
    int lwork;
    B = create_sqmatrix(A.p);
    sqmatrix_copy(B, A);
    CHK_MALLOC(rwork, real, 2*A.p);
    lapackglue_geev('N','N', A.p, B.data, A.p, eigvals, NULL,1,NULL,1,
                    &work1, -1, rwork);
    lwork = (int) (SCALAR_RE(work1) + 0.5);
    CHK_MALLOC(work, scalar, lwork);
    lapackglue_geev('N','N', A.p, B.data, A.p, eigvals, NULL,1,NULL,1,
                    work, lwork, rwork);
    free(work);
    free(rwork);
    destroy_sqmatrix(B);
}

/* Compute Usqrt <- sqrt(U), where U must be Hermitian positive-definite.
   W is a work array, and must be the same size as U.  Both U and
   W are overwritten. */
void sqmatrix_sqrt(sqmatrix Usqrt, sqmatrix U, sqmatrix W)
{
     real *eigenvals;

     sqmatrix_assert_hermitian(U);
     CHECK(Usqrt.p == U.p && U.p == W.p, "matrices not conformant");

     CHK_MALLOC(eigenvals, real, U.p);

     sqmatrix_eigensolve(U, eigenvals, W);

     {
	  int i;

	  /* Compute W = diag(sqrt(eigenvals)) * U; i.e. the rows of W
	     become the rows of U times sqrt(corresponding eigenvalue) */
	  for (i = 0; i < U.p; ++i) {
	       CHECK(eigenvals[i] > 0, "non-positive eigenvalue");

	       blasglue_copy(U.p, U.data + i*U.p, 1, W.data + i*U.p, 1);
	       blasglue_rscal(U.p, sqrt(eigenvals[i]), W.data + i*U.p, 1);
	  }
     }

     free(eigenvals);

     /* compute Usqrt = Ut * W */
     sqmatrix_AeBC(Usqrt, U, 1, W, 0);
}
