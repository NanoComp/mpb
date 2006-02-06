/* Copyright (C) 1999, 2000, 2001, 2002, Massachusetts Institute of Technology.
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
#include <math.h>

#include "config.h"
#include <blasglue.h>
#include <check.h>

extern void debug_check_memory_leaks(void);  

void printmat(scalar *A, int m, int n)
{
  int i, j;

  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
#ifdef SCALAR_COMPLEX
      printf("  (%6.3f,%6.3f)", (double)A[i*n + j].re, (double)A[i*n + j].im);
#else
      printf("  %6.3f", (double)A[i*n + j]);
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
#ifdef SCALAR_COMPLEX
         printf("  %g+%gi", (double)A[i*n + j].re, (double)A[i*n + j].im);
#else
         printf("  %g", (double)A[i*n + j]);
#endif
    }
    printf(";\n");
  }
  printf("]\n");
}

int main(void)
{
  const int N = 4;
  int i,j;
#ifndef SCALAR_COMPLEX
  scalar A[] = { 3.3, 6.2, 7.1, 9.1,
                 -2.3, 3.6, 0.3, 9.7,
                 6.7, -0.1, 1.1, 4.8,
                 8.4, 7.7, 5.9, -1.8 };
  scalar B[] = { 1.1, 2.2, 3.3, 4.4,
		 8.8, 7.7, 6.6, 5.5,
		 6.1, 8.2, 9.7, 3.6,
		 6.3, 2.9, 5.5, 8.1 };
#else
  scalar A[] = { {3.3, 6.2} , {7.1, 9.1}, {2.3, 8.2}, {-3.2, 6.6},
		 {-2.3, 3.6}, {0.3, 9.7}, {1.9,-4.9}, {7.1, 7.1},
		 {6.7, -0.1}, {1.1, 4.8}, {-9.7, 3.7}, {-0.01, -0.2},
		 {8.4, 7.7}, {5.9, -1.8}, {8.8, 9.9}, {0.0, 0.1} };
  scalar B[] = { {1.1, 2.2}, {3.3, 4.4}, {1.2, 2.3}, {3.4, 4.5},
		 {8.8, 7.7}, {6.6, 5.5}, {3.5, 7.2}, {-0.3, 6.1},
		 {6.1, 8.2}, {9.7, 3.6}, {-5.1, 6.1}, {2.3, 8.1},
		 {6.3, 2.9}, {5.5, 8.1}, {8.5, 6.7}, {9.0, 2.4} };
#endif
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
  printf("\nD = A' * A\n");
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

  /* Compute eigenvectors and eigenvalues: */
  lapackglue_heev('V', 'U', N, D, N, eigvals, E, 16, wrk);
  /* Choose a deterministic phase for each row/eigenvector: */
  for (i = 0; i < N; ++i) {
       scalar phase;
       real len;
       for (j = 0; (len = sqrt(SCALAR_NORMSQR(D[i*N + j]))) < 1e-6; ++j)
	    ;
       /* phase to make D[i*N+j] purely real: */
       ASSIGN_SCALAR(phase, SCALAR_RE(D[i*N+j])/len, -SCALAR_IM(D[i*N+j])/len);
       ASSIGN_MULT(D[i*N+j], D[i*N+j], phase);
       if (SCALAR_RE(D[i*N+j]) < 0) { /* pick deterministic (positive) sign */
	    ASSIGN_SCALAR(phase, -SCALAR_RE(phase), -SCALAR_IM(phase));
	    ASSIGN_SCALAR(D[i*N+j], -SCALAR_RE(D[i*N+j]),-SCALAR_IM(D[i*N+j]));
       }
       for (j = j + 1; j < N; ++j)
	    ASSIGN_MULT(D[i*N + j], D[i*N + j], phase);
  }
  printf("\n[v,d] = eig(D);\n");
  printf("\ndiag(d)\n  ");
  for (i = 0; i < 4; ++i) printf("  %6.3f", (double)eigvals[i]);
  printf("\nv'\n");
  printmat(D,N,N);
  blasglue_gemm('C', 'N', N, N, N, 1.0, D, N, D, N, 0.0, C, N);
  printf("\nv * v'\n");
  printmat(C,N,N);

  /* Compute E = diag(sqrt(eigenvals)) * D; i.e. the rows of E
     become the rows of D times sqrt(corresponding eigenvalue) */
  for (i = 0; i < N; ++i) {
    CHECK(eigvals[i] > 0, "non-positive eigenvalue");
    
    blasglue_copy(N, D + i*N, 1, E + i*N, 1);
    blasglue_rscal(N, sqrt(eigvals[i]), E + i*N, 1);
  }

  /* compute C = adjoint(D) * E == sqrt (At * A) */
  blasglue_gemm('C', 'N', N, N, N, 1.0, D, N, E, N, 0.0, C, N);
  printf("\nsqrtm(D)\n");
  printmat(C,N,N);

  blasglue_gemm('C', 'N', N, N, N, 1.0, E, N, E, N, 0.0, C, N);
  printf("\nsqrtm(D) * sqrtm(D)\n");
  printmat(C,N,N);

  debug_check_memory_leaks();

  return EXIT_SUCCESS;
}
