/* Glue code to make the interface to BLAS routines more C-like.

   These routines take care of translating between C and Fortran
   argument conventions, including array formats--the C code can
   call these routines using ordinary row-major convention, and
   the arguments will be translated to Fortran's column-major
   format automatically.  (No data movement is required for this.)

   Note that, where in Fortran you pass the leading dimension ("ld")
   of each array to the routines, in C (here) we pass the final
   dimension ("fd") of each array.

   This code also automatically selects the right version of the BLAS
   routines, depending upon which data type is defined in scalar.h. */

/* This file also includes similar glue code for some LAPACK routines. */

/*************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <fortranize.h>
#include <check.h>

#include "blasglue.h"
#include "scalar.h"

/*************************************************************************/

/* Define a macro F(x,X) that works similarly to the FORTRANIZE
   macro except that it appends an appropriate BLAS prefix (c,z,s,d)
   to the routine name depending upon the type defined in scalar.h */

#ifdef SCALAR_COMPLEX
#  ifdef SCALAR_SINGLE_PREC
#    define F(x,X) FORTRANIZE(c##x, C##X)
#  else
#    define F(x,X) FORTRANIZE(z##x, Z##X)
#  endif
#else
#  ifdef SCALAR_SINGLE_PREC
#    define F(x,X) FORTRANIZE(s##x, S##X)
#  else
#    define F(x,X) FORTRANIZE(d##x, D##X)
#  endif
#endif

/*************************************************************************/

/* Prototypes for BLAS and LAPACK functions.  Note that we need to
   wrap these in extern "C" if this is compiled under C++, or all
   hell will break loose.  (i.e. C++'s name munging will conflict
   with Fortran's.) */

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

extern void F(axpy,AXPY) (int *, scalar *, scalar *, int *, scalar *, int *);
extern void F(scal,SCAL) (int *, scalar *, scalar *, int *);
extern void F(copy,COPY) (int *, scalar *, int *, scalar *, int *);
extern scalar F(dotc,DOTC) (int *, scalar *, int *, scalar *, int *);
extern scalar F(dot,DOT) (int *, scalar *, int *, scalar *, int *);
extern void F(gemm,GEMM) (char *, char *, int *, int *, int *,
			  scalar *, scalar *, int *, scalar *, int *,
			  scalar *, scalar *, int *);
extern void F(herk,HERK) (char *, char *, int *, int *,
			  scalar *, scalar *, int *,
			  scalar *, scalar *, int *);
extern void F(syrk,SYRK) (char *, char *, int *, int *,
			  scalar *, scalar *, int *,
			  scalar *, scalar *, int *);
extern void F(potrf,POTRF) (char *, int *, scalar *, int *, int *);
extern void F(potri,POTRI) (char *, int *, scalar *, int *, int *);
extern void F(heev,HEEV) (char *, char *, int *, scalar *, int *, real *,
			  scalar *, int *, real *, int *);
extern void F(syev,SYEV) (char *, char *, int *, scalar *, int *, real *,
			  scalar *, int *, int *);

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

/*************************************************************************/

void blasglue_axpy(int n, real a, scalar *x, int incx, scalar *y, int incy)
{
     scalar alpha;

     ASSIGN_REAL(alpha, a);

     F(axpy,AXPY) (&n, &alpha, x, &incx, y, &incy);
}

void blasglue_scal(int n, real a, scalar *x, int incx)
{
     scalar alpha;

     ASSIGN_REAL(alpha, a);

     F(scal,SCAL) (&n, &alpha, x, &incx);
}

void blasglue_copy(int n, scalar *x, int incx, scalar *y, int incy)
{
     F(copy,COPY) (&n, x, &incx, y, &incy);
}

scalar blasglue_dotc(int n, scalar *x, int incx, scalar *y, int incy)
{
#ifdef SCALAR_COMPLEX
     return (F(dotc,DOTC) (&n, x, &incx, y, &incy));
#else
     return (F(dot,DOT) (&n, x, &incx, y, &incy));
#endif
}

void blasglue_gemm(char transa, char transb, int m, int n, int k,
		   real a, scalar *A, int fdA, scalar *B, int fdB,
		   real b, scalar *C, int fdC)
{
     scalar alpha, beta;

     CHECK(A != C && B != C, "gemm output array must be distinct");
     
     ASSIGN_REAL(alpha,a);
     ASSIGN_REAL(beta,b);

     F(gemm,GEMM) (&transb, &transa, &n, &m, &k,
		   &alpha, B, &fdB, A, &fdA, &beta, C, &fdC);
}

void blasglue_herk(char uplo, char trans, int n, int k,
		   real a, scalar *A, int fdA,
		   real b, scalar *C, int fdC)
{
     scalar alpha, beta;
     
     CHECK(A != C, "herk output array must be distinct");
     
     ASSIGN_REAL(alpha,a);
     ASSIGN_REAL(beta,b);

     uplo = uplo == 'U' ? 'L' : 'U';
     trans = (trans == 'C' || trans == 'T') ? 'N' : 'C';

#ifdef SCALAR_COMPLEX
     F(herk,HERK) (&uplo, &trans, &n, &k,
		   &alpha, A, &fdA, &beta, C, &fdC);
#else
     F(syrk,SYRK) (&uplo, &trans, &n, &k,
		   &alpha, A, &fdA, &beta, C, &fdC);
#endif
}

/*************************************************************************/

#ifndef NO_LAPACK

void lapackglue_potrf(char uplo, int n, scalar *A, int fdA)
{
     int info;

     uplo = uplo == 'U' ? 'L' : 'U';

     F(potrf,POTRF) (&uplo, &n, A, &fdA, &info);

     CHECK(info >= 0, "invalid argument in potrf");
     CHECK(info <= 0, "non positive-definite matrix in potrf");
}

void lapackglue_potri(char uplo, int n, scalar *A, int fdA)
{
     int info;

     uplo = uplo == 'U' ? 'L' : 'U';

     F(potri,POTRI) (&uplo, &n, A, &fdA, &info);

     CHECK(info >= 0, "invalid argument in potri");
     CHECK(info <= 0, "zero diagonal element (singular matrix) in potri");
}

void lapackglue_heev(char jobz, char uplo, int n, scalar *A, int fdA, 
		     real *w, scalar *work, int lwork, real *rwork)
{
     int info;

     uplo = uplo == 'U' ? 'L' : 'U';

#ifdef SCALAR_COMPLEX
     F(heev,HEEV) (&jobz, &uplo, &n, A, &fdA, w, work, &lwork, rwork, &info);
#else
     F(syev,SYEV) (&jobz, &uplo, &n, A, &fdA, w, work, &lwork, &info);
#endif

     CHECK(info >= 0, "invalid argument in heev");
     CHECK(info <= 0, "failure to converge in heev");
}

#endif
