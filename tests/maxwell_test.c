#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#ifdef HAVE_GETOPT
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#elif defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif
#endif

#include <config.h>
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>
#include <maxwell.h>

#define NX 32
#define NY 1
#define NZ 1
#define NUM_BANDS 8

#define NWORK 3

#define KX 0.5
#define EPS_LOW 2.25
#define EPS_HIGH 12.25
#define EPS_HIGH_X 0.30

#ifdef ENABLE_PROF
#  define PROF_ITERS 10
#else
#  define PROF_ITERS 1
#endif

/*************************************************************************/

/* routines for analytic calculation of Bragg bands: */

static const TWOPI = 6.2831853071795864769252867665590057683943388;

/* We have an analytic expression for k as a function of omega
   for Bragg mirrors.  This will have to be numerically inverted
   to find omega as a function of k.

   n1 and n2 are the indices of the two dielectrics, and f1 and f2
   are their thicknesses as a fraction of the lattice constant
   (we should have f1 + f2 == 1). */
real analytic_bragg_k(real omega, real n1, real f1, real n2, real f2)
{
     real phase1, phase2, c1, s1, c2, s2, b;

     CHECK(fabs(f1 + f2 - 1) < 1e-6, "invalid params to analytic_bragg_k");

     phase1 = TWOPI * n1 * f1 * omega;
     phase2 = TWOPI * n2 * f2 * omega;
     c1 = cos(phase1); s1 = sin(phase1);
     c2 = cos(phase2); s2 = sin(phase2);

     b = c1*c2 - 0.5 * (n1/n2 + n2/n1) * s1*s2;

     if (fabs(b) > 1)
	  return (-1.0);

     return fabs(atan2(sqrt(1-b*b), b) / TWOPI);
}

/* Solve for Bragg omega for the given k and other parameters,
   using omega_guess as a starting guess. 

   We can't do anything clever like Newton's method or even an
   ordinary bisection search because there are regions of omega
   in which analytic_bragg_k is not defined (i.e. in the band gap). */
real bragg_omega(real omega_guess,
		 real k,
		 real n1, real f1, real n2, real f2,
		 real tolerance)
{
     real omega_guess_low = omega_guess - 0.2, 
	  omega_guess_high = omega_guess + 0.2;
     real k_cur;
     real k_best = -1.0, omega_best = 0.0;
     real tol;

     if (omega_guess_low < 0.0)
	  omega_guess_low = 0.0;

     for (tol = (omega_guess_high - omega_guess_low) / 4.0;
	  tol > tolerance;
	  tol *= 0.25) {
	  for (omega_guess = omega_guess_low + tol;
	       omega_guess < omega_guess_high;
	       omega_guess += tol) {
	       k_cur = analytic_bragg_k(omega_guess, n1, f1, n2, f2);
	       if (fabs(k_cur - k) < fabs(k_best - k)) {
		    k_best = k_cur;
		    omega_best = omega_guess;
	       }
	  }

	  CHECK(k_best > 0.0, "No valid omega values in guess range!");

	  omega_guess_low = omega_best - tol;
	  omega_guess_high = omega_best + tol;
     }

     return omega_best;
}


/*************************************************************************/

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

/*************************************************************************/

void usage(void)
{
     printf("Syntax: maxwell_test [options]\n"
	    "Options:\n"
            "   -h           Print this help\n"
	    "   -s <seed>    Set random seed\n"
	    "   -k <kx>      Set kx wavevector component [dflt. = %f]\n"
	    "   -b <n>       Compute n bands [default = %d]\n"
	    "   -n <index>   Specify high-dielectric index [dflt. = %f]\n"
	    "   -f <f>       Specify high-index fill fraction [dflt. = %f]\n"
	    "   -x <nx>      Use nx points in x direction [dflt. = %d]\n"
	    "   -y <ny>      Use ny points in y direction [dflt. = %d]\n"
	    "   -z <nz>      Use nz points in z direction [dflt. = %d]\n"
	    "   -e           Solve for TE polarization only.\n"
	    "   -m           Solve for TM polarization only.\n",
	    KX, NUM_BANDS, sqrt(EPS_HIGH), EPS_HIGH_X, NX, NY, NZ);
}

/*************************************************************************/

int main(int argc, char **argv)
{
     maxwell_data *mdata;
     int local_N, N_start, alloc_N;
     real G1[3] = {1,0,0}, G2[3] = {0,100,0}, G3[3] = {0,0,100};
     real kvector[3] = {KX,0,0};
     evectmatrix H, Hstart, W[NWORK];
     real *eigvals;
     int i, j, k, iters;
     int num_iters;
     polarization_t polarization = NO_POLARIZATION;
     int nx = NX, ny = NY, nz = NZ;
     int num_bands = NUM_BANDS;
     real high_eps = EPS_HIGH;
     real high_index_fill = EPS_HIGH_X;

     srand(time(NULL));

#ifdef HAVE_GETOPT
     {
          extern char *optarg;
          extern int optind;
          int c;

          while ((c = getopt(argc, argv, "hs:k:b:n:f:x:y:z:em")) != -1)
	       switch (c) {
		   case 'h':
			usage();
			exit(EXIT_SUCCESS);
			break;
		   case 's':
			srand(atoi(optarg));
			break;	
		   case 'k':
			kvector[0] = atof(optarg);
			break;
		   case 'b':
			num_bands = atoi(optarg);
			CHECK(num_bands > 0, "num_bands must be positive");
			break;
		   case 'n':
			high_eps = atof(optarg);
			CHECK(high_eps > 0.0, "index must be positive");
			high_eps = high_eps * high_eps;
			break;
		   case 'f':
			high_index_fill = atof(optarg);
			CHECK(high_index_fill > 0.0, "fill must be positive");
			break;
		   case 'x':
			nx = atoi(optarg);
			CHECK(nx > 0, "x size must be positive");
			break;
		   case 'y':
			ny = atoi(optarg);
			CHECK(ny > 0, "y size must be positive");
			break;
		   case 'z':
			nz = atoi(optarg);
			CHECK(nz > 0, "z size must be positive");
			break;
		   case 'e':
			polarization = TE_POLARIZATION;
			break;
		   case 'm':
			polarization = TM_POLARIZATION;
			break;
		   default:
			usage();
			exit(EXIT_FAILURE);
	       }

	  if (argc != optind) {
	       usage();
	       exit(EXIT_FAILURE);
	  }
     }     
#endif

     printf("Creating Maxwell data...\n");
     mdata = create_maxwell_data(nx, ny, nz, &local_N, &N_start, &alloc_N,
				 num_bands, num_bands);
     CHECK(mdata, "NULL mdata");

     set_maxwell_data_polarization(mdata, polarization);

     printf("Setting k vector to (%g, %g, %g)...\n",
	    kvector[0], kvector[1], kvector[2]);
     update_maxwell_data_k(mdata, kvector, G1, G2, G3);

     printf("Initializing dielectric...\n");
     /* set up dielectric function (a simple Bragg mirror) */
     for (i = 0; i < nx; ++i) {
	  real epsinv = (i < high_index_fill*nx) ? (1/high_eps) : (1/EPS_LOW);

	  for (j = 0; j < ny; ++j)
	       for (k = 0; k < nz; ++k) {
		    int ijk = k + nz * (j + ny * i); /* row-major index */
		    mdata->eps_inv[ijk].m00 = epsinv;
		    mdata->eps_inv[ijk].m11 = epsinv;
		    mdata->eps_inv[ijk].m22 = epsinv;
		    mdata->eps_inv[ijk].m01 = 0.0;
		    mdata->eps_inv[ijk].m02 = 0.0;
		    mdata->eps_inv[ijk].m12 = 0.0;
	       }
     }

     printf("Allocating fields...\n");
     H = create_evectmatrix(nx * ny * nz, 2, num_bands,
			    local_N, N_start, alloc_N);
     Hstart = create_evectmatrix(nx * ny * nz, 2, num_bands,
				 local_N, N_start, alloc_N);
     for (i = 0; i < NWORK; ++i)
	  W[i] = create_evectmatrix(nx * ny * nz, 2, num_bands,
				    local_N, N_start, alloc_N);

     eigvals = (real*) malloc(sizeof(real) * num_bands);
     CHECK(eigvals, "out of memory");

     for (iters = 0; iters < PROF_ITERS; ++iters) {

     printf("Initializing fields...\n");
     for (i = 0; i < H.n * H.p; ++i)
          ASSIGN_REAL(Hstart.data[i], rand() * 1.0 / RAND_MAX);

     /*****************************************/
     printf("\nSolving for eigenvectors with preconditioning...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 maxwell_operator, (void *) mdata,
		 maxwell_preconditioner, (void *) mdata, NULL,
		 maxwell_constraint, (void *) mdata,
		 W, NWORK, 1e-4, &num_iters);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("Eigenvalues (sqrt) = \n");
     for (i = 0; i < num_bands; ++i) {
	  printf("   %f (%f)\n", eigvals[i], sqrt(eigvals[i]));
     }
     printf("\n");

     printf("Corresponding analytic frequencies are:\n");
     for (i = 0;  i < num_bands; ++i) {
          printf("   %f\n",  bragg_omega(sqrt(eigvals[i]),
					 kvector[0],
					 sqrt(high_eps),
					 high_index_fill,
					 sqrt(EPS_LOW),
					 1.0 - high_index_fill,
					 1.0e-7));
     }
     printf("\n");

     }

#ifndef ENABLE_PROF
     /*****************************************/

     printf("\nSolving for eigenvectors without preconditioning...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 maxwell_operator, (void *) mdata,
		 NULL, NULL, NULL,
		 maxwell_constraint, (void *) mdata,
		 W, NWORK, 1e-4, &num_iters);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("Eigenvalues (sqrt) = \n");
     for (i = 0; i < num_bands; ++i) {
	  printf("   %f (%f)\n", eigvals[i], sqrt(eigvals[i]));
     }
     printf("\n");

     /*****************************************/
     
     printf("\nSolving for eigenvectors without conj. grad...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 maxwell_operator, (void *) mdata,
		 maxwell_preconditioner, (void *) mdata, NULL,
		 maxwell_constraint, (void *) mdata,
		 W, NWORK - 1, 1e-4, &num_iters);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("Eigenvalues (sqrt) = \n");
     for (i = 0; i < num_bands; ++i) {
	  printf("   %f (%f)\n", eigvals[i], sqrt(eigvals[i]));
     }
     printf("\n");

     /*****************************************/
     printf("\nSolving for eigenvectors without precond. or conj. grad...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 maxwell_operator, (void *) mdata,
		 NULL, NULL, NULL,
		 maxwell_constraint, (void *) mdata,
		 W, NWORK - 1, 1e-4, &num_iters);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("Eigenvalues (sqrt) = \n");
     for (i = 0; i < num_bands; ++i) {
	  printf("   %f (%f)\n", eigvals[i], sqrt(eigvals[i]));
     }
     printf("\n");

     /*****************************************/
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
