#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>

static sqmatrix A;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

extern void Aop(evectmatrix Xin, evectmatrix Xout);
extern void Cop(evectmatrix Xin, evectmatrix Xout);
extern void printmat(scalar *A, int m, int n);
extern void printmat_matlab(scalar *A, int m, int n);

extern void debug_check_memory_leaks(void);

int main(int argc, char **argv)
{
     int i, j, n = 0, p;
     sqmatrix X, U;
     evectmatrix Y, W[3];
     real *eigvals, *eigvals_dense, sum = 0.0;
     int num_iters;

     srand(time(NULL));

     if (argc >= 2)
	  n = atoi(argv[1]);

     CHECK(n > 0, "illegal argument\nSyntax: eigs_test <n>");

     X = create_sqmatrix(n);
     A = create_sqmatrix(n);
     U = create_sqmatrix(n);

     /* fill X with random data */
     for (i = 0; i < n*n; ++i)
	  ASSIGN_REAL(X.data[i], rand() * 1.0 / RAND_MAX);

     /* assign A = adjoint(X) * X to get a Hermitian matrix: */
     sqmatrix_AeBC(A, X, 1, X, 0);
     sqmatrix_AeBC(U, X, 1, X, 0); /* also assign U */
     
     if (n <= 20) {
	  printf("Solving for eigenvalues of %d x %d matrix: \nA = ", n, n);
	  printmat_matlab(A.data, n, n);
     }

     eigvals_dense = (real*) malloc(sizeof(real) * n);
     CHECK(eigvals_dense, "out of memory");

     sqmatrix_eigensolve(U, eigvals_dense, X);

     p = MIN(MIN(5, MAX(n/4, 2)), n);
     printf("\nSolving for %d eigenvals out of %d.\n", p, n);

     printf("\nSolved A by dense eigensolver.\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals_dense[i];
       printf("  %f", eigvals_dense[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);
     printf("\nEigenvectors are (by row): \n");
     printmat(U.data, p, n);

     Y = create_evectmatrix(n, 1, p, n, 0, n);
     for (i = 0; i < 3; ++i)
	  W[i] = create_evectmatrix(n, 1, p, n, 0, n);
     eigvals = (real*) malloc(sizeof(real) * p);
     CHECK(eigvals, "out of memory");

     for (i = 0; i < n*p; ++i)
	  ASSIGN_REAL(Y.data[i], rand() * 1.0 / RAND_MAX);

     eigensolver(Y, eigvals, Aop, Cop, W, 3, 1e-10, &num_iters);

     printf("\nSolved for eigenvectors after %d iterations.\n", num_iters);
     printf("\nEigenvalues = ");
     for (sum = 0.0, i = 0; i < p; ++i) {
       sum += eigvals[i];
       printf("  %f", eigvals[i]);
     }
     printf("\nEigenvalue sum = %f\n", sum);
     printf("Eigenvectors are (by column): \n");
     printmat(Y.data, n, p);
     
     destroy_sqmatrix(A);
     destroy_sqmatrix(X);
     destroy_sqmatrix(U);
     destroy_evectmatrix(Y);
     for (i = 0; i < 3; ++i)
	  destroy_evectmatrix(W[i]);

     free(eigvals);
     free(eigvals_dense);

     debug_check_memory_leaks();

     return EXIT_SUCCESS;
}

void Aop(evectmatrix Xin, evectmatrix Xout)
{
     CHECK(A.p == Xin.n && A.p == Xout.n && Xin.p == Xout.p,
	   "matrices not conformant");

     blasglue_gemm('N', 'N', Xout.n, Xout.p, Xin.n,
		   1.0, A.data, A.p, Xin.data, Xin.p, 0.0, Xout.data, Xout.p);
}

void Cop(evectmatrix Xin, evectmatrix Xout)
{
     int in, ip;

     CHECK(A.p == Xin.n && A.p == Xout.n && Xin.p == Xout.p,
           "matrices not conformant");

#ifdef NO_PRECONDITIONER
     evectmatrix_copy(Xout, Xin);
#else
     for (in = 0; in < Xin.n; ++in)
	  for (ip = 0; ip < Xin.p; ++ip) {
	       scalar xin = Xin.data[in * Xin.p + ip];
	       real diag;
	       
	       diag = SCALAR_NORMSQR(A.data[in * A.p + in]);
	       diag = (diag == 0.0) ? 1.0 : 1.0 / sqrt(diag);
	       
	       ASSIGN_SCALAR(Xout.data[in * Xout.p + ip],
			     diag * SCALAR_RE(xin),
			     diag * SCALAR_IM(xin));
	       }
#endif
}

void printmat(scalar *A, int m, int n)
{
  int i, j;

  for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) {
#ifdef SCALAR_COMPLEX
	 printf("  (%6.3f,%6.3f)", A[i*n + j].re, A[i*n + j].im);
#else
	 printf("  %6.3f", A[i*n + j]);
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
	 printf("  %f", SCALAR_RE(A[i*n + j]));
       }
    printf(";\n");
  }
  printf("]\n");
}
