#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include <config.h>
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>
#include <maxwell.h>

#define NX 8
#define NY 1
#define NZ 1
#define NUM_BANDS 8

#define NWORK 3

#define EPS_LOW 1.0
#define EPS_HIGH 9.0
#define EPS_HIGH_X 0.25

/* placeholder routine for when we don't want to do preconditioning */
void identity_op(evectmatrix Xin, evectmatrix Xout, void *data)
{
     evectmatrix_copy(Xout, Xin);
}

int main(int argc, char **argv)
{
     maxwell_data *mdata;
     int local_N, N_start, alloc_N;
     real G1[3] = {1,0,0}, G2[3] = {0,100,0}, G3[3] = {0,0,100};
     real kvector[3] = {0.5,0,0};
     evectmatrix H, Hstart, W[NWORK];
     real *eigvals;
     int i, j, k;
     int num_iters;

     printf("Syntax: maxwell_test [<kx>] [<seed>]\n");

     if (argc > 1)
	  kvector[0] = atof(argv[1]);

     srand(argc > 2 ? atoi(argv[2]) : time(NULL));
     
     printf("Creating Maxwell data...\n");
     mdata = create_maxwell_data(NX, NY, NZ, &local_N, &N_start, &alloc_N,
				 NUM_BANDS, NUM_BANDS);
     CHECK(mdata, "NULL mdata");

     printf("Setting k vector to (%g, %g, %g)...\n",
	    kvector[0], kvector[1], kvector[2]);
     update_maxwell_data_k(mdata, kvector, G1, G2, G3);

     printf("Initializing dielectric...\n");
     /* set up dielectric function (a simple Bragg mirror) */
     for (i = 0; i < NX; ++i) {
	  real epsinv = (i < EPS_HIGH_X * NX) ? (1/EPS_HIGH) : (1/EPS_LOW);

	  for (j = 0; j < NY; ++j)
	       for (k = 0; k < NZ; ++k) {
		    int ijk = k + NZ * (j + NY * i); /* row-major index */
		    mdata->eps_inv[ijk].m00 = epsinv;
		    mdata->eps_inv[ijk].m11 = epsinv;
		    mdata->eps_inv[ijk].m22 = epsinv;
		    mdata->eps_inv[ijk].m01 = 0.0;
		    mdata->eps_inv[ijk].m02 = 0.0;
		    mdata->eps_inv[ijk].m12 = 0.0;
	       }
     }
		    

     printf("Allocating fields...\n");
     H = create_evectmatrix(NX * NY * NZ, 2, NUM_BANDS,
			    local_N, N_start, alloc_N);
     Hstart = create_evectmatrix(NX * NY * NZ, 2, NUM_BANDS,
				 local_N, N_start, alloc_N);
     for (i = 0; i < NWORK; ++i)
	  W[i] = create_evectmatrix(NX * NY * NZ, 2, NUM_BANDS,
				    local_N, N_start, alloc_N);

     eigvals = (real*) malloc(sizeof(real) * NUM_BANDS);
     CHECK(eigvals, "out of memory");

     printf("Initializing fields...\n");
     for (i = 0; i < H.n * H.p; ++i)
          ASSIGN_REAL(Hstart.data[i], rand() * 1.0 / RAND_MAX);

     printf("Solving for eigenvectors...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 maxwell_operator, (void *) mdata,
		 identity_op, NULL,
		 W, NWORK, 1e-4, &num_iters);

     printf("\nSolved for eigenvectors after %d iterations.\n", num_iters);
     printf("\nEigenvalues (sqrt) = \n");
     for (i = 0; i < NUM_BANDS; ++i) {
	  printf("   %f (%f)\n", eigvals[i], sqrt(eigvals[i]));
     }
     printf("\n");

#if 0
     printf("Eigenvectors:\n");
     for (i = 0; i < (local_N > 8 ? 8 : local_N); ++i) {
	  for (j = 0; j < (H.p > 2 ? 2 : H.p); ++j)
	       printf("   [ (%g,%g) (%g,%g) ]",
		      H.data[(i*2) * H.p + j].re, H.data[(i*2) * H.p + j].im,
		      H.data[(i*2+1)*H.p + j].re, H.data[(i*2+1)*H.p + j].im);
	  printf("\n");
     }
#endif
     
     destroy_evectmatrix(H);
     destroy_evectmatrix(Hstart);
     for (i = 0; i < NWORK; ++i)
          destroy_evectmatrix(W[i]);

     destroy_maxwell_data(mdata);

     free(eigvals);

     debug_check_memory_leaks();

     return EXIT_SUCCESS;
}
