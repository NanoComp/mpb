#ifndef BLASGLUE_H
#define BLASGLUE_H

#include "scalar.h"

extern void blasglue_axpy(int n, real a, scalar *x, int incx,
			  scalar *y, int incy);
extern void blasglue_scal(int n, real a, scalar *x, int incx);
extern void blasglue_copy(int n, scalar *x, int incx, scalar *y, int incy);
extern scalar blasglue_dotc(int n, scalar *x, int incx, scalar *y, int incy);
void blasglue_gemm(char transa, char transb, int m, int n, int k,
                   real a, scalar *A, int fdA, scalar *B, int fdB,
                   real b, scalar *C, int fdC);
extern void blasglue_herk(char uplo, char trans, int n, int k,
			  real a, scalar *A, int fdA,
			  real b, scalar *C, int fdC);
extern void lapackglue_potrf(char uplo, int n, scalar *A, int fdA);
extern void lapackglue_potri(char uplo, int n, scalar *A, int fdA);
extern void lapackglue_heev(char jobz, char uplo, int n, scalar *A, int fdA,
			    real *w, scalar *work, int lwork, real *rwork);

#endif /* BLASGLUE_H */
