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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "../src/config.h"
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>

static sqmatrix A, Ainv;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

extern void Aop(evectmatrix Xin, evectmatrix Xout, void *data,
		int is_current_eigenvector, evectmatrix Work);
extern void Ainvop(evectmatrix Xin, evectmatrix Xout, void *data,
		   evectmatrix Y, real *eigenvals);
extern void Cop(evectmatrix Xin, evectmatrix Xout, void *data,
		evectmatrix Y, real *eigenvals, sqmatrix YtY);
extern void printmat(scalar *A, int m, int n, int ldn);
extern void printmat_matlab(scalar *A, int m, int n);

extern void debug_check_memory_leaks(void);

real norm_diff(scalar *a, scalar *b, int n)
{
     real bmag = 0.0, diffmag = 0.0;
     int i;
     for (i = 0; i < n; ++i) {
          scalar d;
          ASSIGN_SCALAR(d,
                        SCALAR_RE(b[i]) - SCALAR_RE(a[i]),
                        SCALAR_IM(b[i]) - SCALAR_IM(a[i]));
          bmag += SCALAR_NORMSQR(b[i]);
          diffmag += SCALAR_NORMSQR(d);
     }
     return sqrt(diffmag / bmag);
}

#define NWORK 3

int main(int argc, char **argv)
{
     int i, j, n = 0, p;
     sqmatrix X, U, YtY;
     evectmatrix Y, Y2, Ystart, W[NWORK];
     real *eigvals, *eigvals_dense, sum = 0.0;
     int num_iters;

     if (argc >= 2)
	  n = atoi(argv[1]);

     srand(argc >= 3 ? atoi(argv[2]) : time(NULL));

     CHECK(n > 0, "illegal argument\nSyntax: eigs_test <n> [<seed>]");

     X = create_sqmatrix(n);
     A = create_sqmatrix(n);
     Ainv = create_sqmatrix(n);
     U = create_sqmatrix(n);

     /* fill X with random data */
     for (i = 0; i < n*n; ++i)
	  ASSIGN_SCALAR(X.data[i],
			rand() * 1.0 / RAND_MAX,
			rand() * 1.0 / RAND_MAX);

     /* assign A = adjoint(X) * X to get a Hermitian matrix: */
     sqmatrix_AeBC(A, X, 1, X, 0);

     /* increase diagonal elements of A so that our preconditioner
	has a chance of being useful: */
     for (i = 0; i < n; ++i)
       ASSIGN_SCALAR(A.data[i * n + i],
		     n * SCALAR_RE(A.data[i * n + i]),
		     n * SCALAR_IM(A.data[i * n + i]));

     sqmatrix_copy(U, A);

     sqmatrix_copy(Ainv, A);
     sqmatrix_invert(Ainv);
     
     if (n <= 10) {
	  printf("Solving for eigenvalues of %d x %d matrix: \nA = ", n, n);
	  printmat_matlab(A.data, n, n);
     }

     CHK_MALLOC(eigvals_dense, real, n);

     sqmatrix_eigensolve(U, eigvals_dense, X);

     /* The eigenvectors are actually the columns of U'.  Assign U = U': */
     for (i = 0; i < n; ++i)
	  for (j = i + 1; j < n; ++j) {
	       scalar dummy;
	       dummy = U.data[i*n + j];
	       U.data[i*n + j] = U.data[j*n + i];
	       U.data[j*n + i] = dummy;
	  }
     for (i = 0; i < n * n; ++i)
	  ASSIGN_CONJ(U.data[i], U.data[i]);

     p = MIN(MIN(5, MAX(n/4, 2)), n);
     printf("\nSolving for %d eigenvals out of %d.\n", p, n);

     printf("\nSolved A by dense eigensolver.\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals_dense[i];
       printf("  %f", eigvals_dense[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);
     printf("\nEigenvectors are (by column): \n");
     printmat(U.data, n, p, n);

     YtY = create_sqmatrix(p);

     Y = create_evectmatrix(n, 1, p, n, 0, n);
     Y2 = create_evectmatrix(n, 1, p, n, 0, n);
     Ystart = create_evectmatrix(n, 1, p, n, 0, n);
     for (i = 0; i < NWORK; ++i)
	  W[i] = create_evectmatrix(n, 1, p, n, 0, n);
     CHK_MALLOC(eigvals, real, p);

     for (i = 0; i < n*p; ++i)
	  ASSIGN_REAL(Ystart.data[i], rand() * 1.0 / RAND_MAX);

     /* Check inverse Ainvop: */
     Aop(Ystart, Y, NULL, 0, Y2);
     Ainvop(Y, Y2, NULL, Ystart, NULL);
     printf("\nDifference |Y - (1/A)*A*Y| / |Y| = %g\n",
	    norm_diff(Ystart.data, Y2.data, Y.n * Y.p));

     evectmatrix_copy(Y, Ystart);
     eigensolver(Y, eigvals, Aop, NULL, Cop, NULL, NULL, NULL,
		 W, NWORK, 1e-10,&num_iters, EIGS_DEFAULT_FLAGS);
     printf("\nSolved for eigenvectors after %d iterations.\n", num_iters);
     printf("\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals[i];
       printf("  %f", eigvals[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);

     /* Change phase of eigenvectors to match those solved for previously: */
     for (i = 0; i < p; ++i) {
	  scalar phase;

	  ASSIGN_DIV(phase, U.data[i], Y.data[i]);

	  for (j = 0; j < n; ++j) {
	       ASSIGN_MULT(Y.data[j*p + i], Y.data[j*p + i], phase);
	  }
     }

     printf("Eigenvectors are (by column): \n");
     printmat(Y.data, n, p, p);
     evectmatrix_XtX(YtY, Y);
     printf("adjoint(Y) * Y:\n");
     printmat(YtY.data, p, p, p);

     printf("\nSolving with exact inverse preconditioner...\n");
     evectmatrix_copy(Y, Ystart);
     eigensolver(Y, eigvals, Aop, NULL, Ainvop, NULL, NULL, NULL,
		 W, NWORK, 1e-10, &num_iters, EIGS_DEFAULT_FLAGS);
     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals[i];
       printf("  %f", eigvals[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);

     printf("\nSolving without conjugate-gradient...\n");
     evectmatrix_copy(Y, Ystart);
     eigensolver(Y, eigvals, Aop, NULL, Cop, NULL, NULL, NULL,
		 W, NWORK - 1, 1e-10, &num_iters, EIGS_DEFAULT_FLAGS);
     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals[i];
       printf("  %f", eigvals[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);

     printf("\nSolving without preconditioning...\n");
     evectmatrix_copy(Y, Ystart);
     eigensolver(Y, eigvals, Aop, NULL, NULL, NULL, NULL, NULL,
		 W, NWORK, 1e-10, &num_iters, EIGS_DEFAULT_FLAGS);
     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals[i];
       printf("  %f", eigvals[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);

     printf("\nSolving without conjugate-gradient or preconditioning...\n");
     evectmatrix_copy(Y, Ystart);
     eigensolver(Y, eigvals, Aop, NULL, NULL, NULL, NULL, NULL,
		 W, NWORK - 1, 1e-10, &num_iters, EIGS_DEFAULT_FLAGS);
     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals[i];
       printf("  %f", eigvals[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);

     destroy_sqmatrix(A);
     destroy_sqmatrix(Ainv);
     destroy_sqmatrix(X);
     destroy_sqmatrix(U);
     destroy_sqmatrix(YtY);     
     destroy_evectmatrix(Y);
     destroy_evectmatrix(Y2);
     destroy_evectmatrix(Ystart);
     for (i = 0; i < NWORK; ++i)
	  destroy_evectmatrix(W[i]);

     free(eigvals);
     free(eigvals_dense);

     debug_check_memory_leaks();

     return EXIT_SUCCESS;
}

void Aop(evectmatrix Xin, evectmatrix Xout, void *data,
	 int is_current_eigenvector, evectmatrix Work)
{
     CHECK(A.p == Xin.n && A.p == Xout.n && Xin.p == Xout.p,
	   "matrices not conformant");

     blasglue_gemm('N', 'N', Xout.n, Xout.p, Xin.n,
		   1.0, A.data, A.p, Xin.data, Xin.p, 0.0, Xout.data, Xout.p);
}

void Ainvop(evectmatrix Xin, evectmatrix Xout, void *data,
	    evectmatrix Y, real *eigenvals)
{
     CHECK(Ainv.p == Xin.n && Ainv.p == Xout.n && Xin.p == Xout.p,
	   "matrices not conformant");

     blasglue_gemm('N', 'N', Xout.n, Xout.p, Xin.n,
		   1.0, Ainv.data, Ainv.p,
		   Xin.data, Xin.p, 0.0, Xout.data, Xout.p);
}

void Cop_old(evectmatrix Xin, evectmatrix Xout, void *data,
	 evectmatrix Y, real *eigenvals, sqmatrix YtY)
{
     int in, ip;

     CHECK(A.p == Xin.n && A.p == Xout.n && Xin.p == Xout.p,
           "matrices not conformant");

     evectmatrix_XeYS(Xout, Xin, YtY, 1);

     for (in = 0; in < Xout.n; ++in) {
	  real diag;
	  
	  diag = SCALAR_NORMSQR(A.data[in * A.p + in]);
	  diag = (diag == 0.0) ? 1.0 : 1.0 / sqrt(diag);
	  
	  for (ip = 0; ip < Xout.p; ++ip) {
	       scalar xin = Xout.data[in * Xout.p + ip];
	       ASSIGN_SCALAR(Xout.data[in * Xout.p + ip],
			     diag * SCALAR_RE(xin),
			     diag * SCALAR_IM(xin));
	  }
     }
}

void Cop(evectmatrix Xin, evectmatrix Xout, void *data,
	  evectmatrix Y, real *eigenvals, sqmatrix YtY)
{
     int in, ip;

     CHECK(A.p == Xin.n && A.p == Xout.n && Xin.p == Xout.p,
           "matrices not conformant");

     evectmatrix_XeYS(Xout, Xin, YtY, 1);

     for (in = 0; in < Xout.n; ++in) {
	  scalar diag = A.data[in * A.p + in];
	  
	  for (ip = 0; ip < Xout.p; ++ip) {
	       scalar scale;
	       ASSIGN_SCALAR(scale,
			     SCALAR_RE(diag) - 0*eigenvals[ip],
			     SCALAR_IM(diag));
	       ASSIGN_DIV(Xout.data[in * Xout.p + ip],
			  Xout.data[in * Xout.p + ip],
			  scale);
	  }
     }
}

void printmat(scalar *A, int m, int n, int ldn)
{
  int i, j;

  for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) {
#ifdef SCALAR_COMPLEX
	 printf("  (%6.3f,%6.3f)", A[i*ldn + j].re, A[i*ldn + j].im);
#else
	 printf("  %6.3f", A[i*ldn + j]);
#endif

	 if (j > 7) {
	   printf("  ...");
	   break;
	 }
       }
       printf("\n");
       if (i > 7) {
	 printf("   ...\n");
	 break;
       }
  }
}

void printmat_matlab(scalar *A, int m, int n)
{
  int i, j;

  printf("[");
  for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) {
#ifdef SCALAR_COMPLEX
	 printf("  %g+%gi", A[i*n + j].re, A[i*n + j].im);
#else
	 printf("  %g", A[i*n + j]);
#endif
       }
    printf(";\n");
  }
  printf("]\n");
}
