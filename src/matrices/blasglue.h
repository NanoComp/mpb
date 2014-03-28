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

#ifndef BLASGLUE_H
#define BLASGLUE_H

#include "scalar.h"

extern void blasglue_axpy(int n, real a, scalar *x, int incx,
			  scalar *y, int incy);
extern void blasglue_scal(int n, scalar a, scalar *x, int incx);
extern void blasglue_rscal(int n, real a, scalar *x, int incx);
extern void blasglue_copy(int n, scalar *x, int incx, scalar *y, int incy);
extern scalar blasglue_dotc(int n, scalar *x, int incx, scalar *y, int incy);
void blasglue_gemm(char transa, char transb, int m, int n, int k,
                   real a, scalar *A, int fdA, scalar *B, int fdB,
                   real b, scalar *C, int fdC);
extern void blasglue_herk(char uplo, char trans, int n, int k,
			  real a, scalar *A, int fdA,
			  real b, scalar *C, int fdC);
extern int lapackglue_potrf(char uplo, int n, scalar *A, int fdA);
extern int lapackglue_potri(char uplo, int n, scalar *A, int fdA);
extern int lapackglue_hetrf(char uplo, int n, scalar *A, int fdA,
			     int *ipiv, scalar *work, int lwork);
extern int lapackglue_hetri(char uplo, int n, scalar *A, int fdA,
			     int *ipiv, scalar *work);
extern void lapackglue_heev(char jobz, char uplo, int n, scalar *A, int fdA,
			    real *w, scalar *work, int lwork, real *rwork);
extern void lapackglue_syev(char jobz, char uplo, int n, real *A, int fdA,
			    real *w, real *work, int lwork);

#endif /* BLASGLUE_H */
