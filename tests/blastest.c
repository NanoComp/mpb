#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <blasglue.h>
#include <check.h>

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
    }
    printf("\n");
  }
}

void printmat_matlab(scalar *A, int m, int n)
{
  int i, j;

  printf("[");
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      printf("  %6.3f", SCALAR_RE(A[i*n + j]));
    }
    printf(";\n");
  }
  printf("]\n");
}

int main(int argc, char **argv)
{
  const int N = 4;
  int i,j;
  scalar A[] = { 3.3, 6.2, 7.1, 9.1,
		 -2.3, 3.6, 0.3, 9.7,
		 6.7, -0.1, 1.1, 4.8,
		 8.4, 7.7, 5.9, -1.8 };
  scalar B[] = { 1.1, 2.2, 3.3, 4.4,
		 8.8, 7.7, 6.6, 5.5,
		 6.1, 8.2, 9.7, 3.6,
		 6.3, 2.9, 5.5, 8.1 };
  scalar C[16], D[16], E[16];
  real eigvals[4], wrk[20];

  printf("A = "); printmat_matlab(A,N,N);
  printf("B = "); printmat_matlab(B,N,N);

  blasglue_gemm('N', 'N', N, N, N, 1.0, A, N, B, N, 0.0, C, N);
  printf("\nC = A * B\n");
  printmat(C,N,N);

  blasglue_gemm('N', 'N', N, N, N, 1.0, C, N, B, N, 0.0, D, N);
  printf("\nC * B\n");
  printmat(D,N,N);

  blasglue_herk('U', 'C', N, N, 1.0, A, N, 0.0, D, N);
  /* Now, copy the conjugate of the upper half
     onto the lower half of D */
  for (i = 0; i < N; ++i)
    for (j = i + 1; j < N; ++j) {
      ASSIGN_CONJ(D[j * N + i], D[i * N + j]);
    }
  printf("\nD = transpose(A) * A\n");
  printmat(D,N,N);

  lapackglue_potrf('U', N, D, N);
  lapackglue_potri('U', N, D, N);
  /* Now, copy the conjugate of the upper half
     onto the lower half of D */
  for (i = 0; i < N; ++i)
    for (j = i + 1; j < N; ++j) {
      ASSIGN_CONJ(D[j * N + i], D[i * N + j]);
    }
  printf("\ninverse(D)\n");
  printmat(D,N,N);

  /* D = At * A, again */
  blasglue_herk('U', 'C', N, N, 1.0, A, N, 0.0, D, N);
  lapackglue_heev('V', 'U', N, D, N, eigvals, E, 16, wrk);
  printf("\neigenvals of D: ");
  for (i = 0; i < 4; ++i) printf("  %6.3f", eigvals[i]);
  printf("\neigenvects of D: \n");
  printmat(D,N,N);

  /* Compute E = diag(sqrt(eigenvals)) * D; i.e. the rows of E
     become the rows of D times sqrt(corresponding eigenvalue) */
  for (i = 0; i < N; ++i) {
    CHECK(eigvals[i] > 0, "non-positive eigenvalue");
    
    blasglue_copy(N, D + i*N, 1, E + i*N, 1);
    blasglue_scal(N, sqrt(eigvals[i]), E + i*N, 1);
  }

  /* compute C = adjoint(D) * E == sqrt (At * A) */
  blasglue_gemm('C', 'N', N, N, N, 1.0, D, N, E, N, 0.0, C, N);
  printf("\nsqrtm(D)\n");
  printmat(C,N,N);

  blasglue_gemm('C', 'N', N, N, N, 1.0, E, N, E, N, 0.0, C, N);
  printf("\nsqrtm(D) * sqrtm(D)\n");
  printmat(C,N,N);

  return EXIT_SUCCESS;
}
