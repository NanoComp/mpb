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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "../src/config.h"
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>
#include <maxwell.h>

#if defined(HAVE_GETOPT_H)
#  include <getopt.h>
#endif
#if defined(HAVE_UNISTD_H)
#  include <unistd.h>
#endif

#define NX 32
#define NY 1
#define NZ 1
#define NUM_BANDS 8
#define NUM_FFT_BANDS 5

#define NWORK 3

#define KX 0.5
#define EPS_LOW 1.00
#define EPS_HIGH 9.00
#define EPS_HIGH_X 0.25

#define ERROR_TOL 1e-4

#ifdef ENABLE_PROF
#  define PROF_ITERS 10
#else
#  define PROF_ITERS 1
#endif

#define MESH_SIZE 7

/*************************************************************************/

typedef struct {
     real eps_high, eps_low, eps_high_x;
} epsilon_data;

#define INVERSION_SYM 1

static void epsilon(symmetric_matrix *eps, symmetric_matrix *eps_inv,
		    real r[3], void *edata_v)
{
     epsilon_data *edata = (epsilon_data *) edata_v;
     real eps_val;

#if INVERSION_SYM
     if (fabs(r[0]) < 0.5*edata->eps_high_x 
	 || fabs(r[0]-1.0) < 0.5*edata->eps_high_x)
	  eps_val = edata->eps_high;
#else
     if ((r[0] < edata->eps_high_x && r[0] >= 0.0) ||
	 (r[0] >= 1.0 && r[0] - 1.0 < edata->eps_high_x))
	  eps_val = edata->eps_high;
#endif
     else
	  eps_val = edata->eps_low;
     eps->m00 = eps->m11 = eps->m22 = eps_val;
     eps_inv->m00 = eps_inv->m11 = eps_inv->m22 = 1.0 / eps_val;
#ifdef WITH_HERMITIAN_EPSILON
     CASSIGN_ZERO(eps->m01);
     CASSIGN_ZERO(eps->m02);
     CASSIGN_ZERO(eps->m12);
     CASSIGN_ZERO(eps_inv->m01);
     CASSIGN_ZERO(eps_inv->m02);
     CASSIGN_ZERO(eps_inv->m12);
#else
     eps->m01 = eps->m02 = eps->m12 = 0.0;
     eps_inv->m01 = eps_inv->m02 = eps_inv->m12 = 0.0;
#endif
}

/*************************************************************************/

/* routines for analytic calculation of Bragg bands: */

static const double TWOPI = 6.2831853071795864769252867665590057683943388;

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

     for (tol = (omega_guess_high - omega_guess_low) / 10.0;
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
	    "   -m           Solve for TM polarization only.\n"
	    "   -t <freq>    Set target frequency [dflt. none].\n"
	    "   -c <tol>     Set convergence tolerance [dflt. %e].\n"
	    "   -g <NMESH>   Set mesh size [dflt. %d].\n"
	    "   -1           Stop after first computation.\n"
	    "   -p           Use simple preconditioner.\n"
	    "   -E <err>     Exit with error if the error exceeds <err>\n"
	    "   -v           Verbose output.\n",
	    KX, NUM_BANDS, sqrt(EPS_HIGH), EPS_HIGH_X, NX, NY, NZ,
	    ERROR_TOL, MESH_SIZE);
}

/*************************************************************************/

int main(int argc, char **argv)
{
     maxwell_data *mdata;
     maxwell_target_data *mtdata = NULL;
     int local_N, N_start, alloc_N;
     real R[3][3] = { {1,0,0}, {0,0.01,0}, {0,0,0.01} };
     real G[3][3] = { {1,0,0}, {0,100,0}, {0,0,100} };
     real kvector[3] = {KX,0,0};
     evectmatrix H, Hstart, W[NWORK];
     real *eigvals;
     int i, iters;
     int num_iters;
     int parity = NO_PARITY;
     int nx = NX, ny = NY, nz = NZ;
     int num_bands = NUM_BANDS;
     real target_freq = 0.0;
     int do_target = 0;
     evectoperator op;
     evectpreconditioner pre_op;
     void *op_data, *pre_op_data;
     real error_tol = ERROR_TOL;
     int mesh_size = MESH_SIZE, mesh[3];
     epsilon_data ed;
     int stop1 = 0;
     int verbose = 0;
     int which_preconditioner = 2;
     double max_err = 1e20;

     srand(time(NULL));

     ed.eps_high = EPS_HIGH;
     ed.eps_low = EPS_LOW;
     ed.eps_high_x = EPS_HIGH_X;

#ifdef HAVE_GETOPT
     {
          extern char *optarg;
          extern int optind;
          int c;

          while ((c = getopt(argc, argv, "hs:k:b:n:f:x:y:z:emt:c:g:1pvE:"))
		 != -1)
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
			ed.eps_high = atof(optarg);
			CHECK(ed.eps_high > 0.0, "index must be positive");
			ed.eps_high = ed.eps_high * ed.eps_high;
			break;
		   case 'f':
			ed.eps_high_x = atof(optarg);
			CHECK(ed.eps_high_x > 0.0, "fill must be positive");
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
			parity = EVEN_Z_PARITY;
			break;
		   case 'm':
			parity = ODD_Z_PARITY;
			break;
		   case 't':
			target_freq = fabs(atof(optarg));
			do_target = 1;
			break;
		   case 'E':
			max_err = fabs(atof(optarg));
			CHECK(max_err > 0, "maximum error must be positive");
			break;
		   case 'c':
			error_tol = fabs(atof(optarg));
			break;
		   case 'g':
			mesh_size = atoi(optarg);
			CHECK(mesh_size > 0, "mesh size must be positive");
			break;
		   case '1':
			stop1 = 1;
			break;
		   case 'p':
			which_preconditioner = 1;
			break;
		   case 'v':
			verbose = 1;
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

#ifdef ENABLE_PROF
     stop1 = 1;
#endif

     mesh[0] = mesh[1] = mesh[2] = mesh_size;

     printf("Creating Maxwell data...\n");
     mdata = create_maxwell_data(nx, ny, nz, &local_N, &N_start, &alloc_N,
				 num_bands, NUM_FFT_BANDS);
     CHECK(mdata, "NULL mdata");

     set_maxwell_data_parity(mdata, parity);

     printf("Setting k vector to (%g, %g, %g)...\n",
	    kvector[0], kvector[1], kvector[2]);
     update_maxwell_data_k(mdata, kvector, G[0], G[1], G[2]);

     printf("Initializing dielectric...\n");
     /* set up dielectric function (a simple Bragg mirror) */
     set_maxwell_dielectric(mdata, mesh, R, G, epsilon, &ed);

     if (verbose && ny == 1 && nz == 1) {
	  printf("dielectric function:\n");
	  for (i = 0; i < nx; ++i) {
	       if (mdata->eps_inv[i].m00 == mdata->eps_inv[i].m11)
		    printf("  eps(%g) = %g\n", i * 1.0 / nx, 
			   1.0/mdata->eps_inv[i].m00);
	  
	       else
		    printf("  eps(%g) = x: %g OR y: %g\n", i * 1.0 / nx, 
			   1.0/mdata->eps_inv[i].m00,
			   1.0/mdata->eps_inv[i].m11);
	  }
	  printf("\n");
     }

     printf("Allocating fields...\n");
     H = create_evectmatrix(nx * ny * nz, 2, num_bands,
			    local_N, N_start, alloc_N);
     Hstart = create_evectmatrix(nx * ny * nz, 2, num_bands,
				 local_N, N_start, alloc_N);
     for (i = 0; i < NWORK; ++i)
	  W[i] = create_evectmatrix(nx * ny * nz, 2, num_bands,
				    local_N, N_start, alloc_N);

     CHK_MALLOC(eigvals, real, num_bands);

     for (iters = 0; iters < PROF_ITERS; ++iters) {

     printf("Initializing fields...\n");
     for (i = 0; i < H.n * H.p; ++i)
          ASSIGN_REAL(Hstart.data[i], rand() * 1.0 / RAND_MAX);

     /*****************************************/
     if (do_target) {
	  printf("\nSolving for eigenvectors close to %f...\n", target_freq);
	  mtdata = create_maxwell_target_data(mdata, target_freq);
	  op = maxwell_target_operator;
	  if (which_preconditioner == 1)
	       pre_op = maxwell_target_preconditioner;
	  else
	       pre_op = maxwell_target_preconditioner2;
	  op_data = (void *) mtdata;
	  pre_op_data = (void *) mtdata;
     }
     else {
	  op = maxwell_operator;
	  if (which_preconditioner == 1)
	       pre_op = maxwell_preconditioner;
	  else
	       pre_op = maxwell_preconditioner2;
	  op_data = (void *) mdata;
	  pre_op_data = (void *) mdata;
     }

     /*****************************************/
     printf("\nSolving for eigenvectors with preconditioning...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 op, op_data,
		 pre_op, pre_op_data,
		 maxwell_parity_constraint, (void *) mdata,
		 W, NWORK, error_tol, &num_iters, EIGS_DEFAULT_FLAGS);

     if (do_target)
	  eigensolver_get_eigenvals(H, eigvals, maxwell_operator, mdata,
				    W[0], W[1]);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("%15s%15s%15s%15s\n","eigenval", "frequency", "exact freq.", 
	    "error");
     for (i = 0; i < num_bands; ++i) {
	  double err;
	  real freq = sqrt(eigvals[i]);
	  real exact_freq = bragg_omega(freq, kvector[0], sqrt(ed.eps_high),
					ed.eps_high_x, sqrt(ed.eps_low),
					1.0 - ed.eps_high_x, 1.0e-7);
	  printf("%15f%15f%15f%15e\n", eigvals[i], freq, exact_freq,
		 err = fabs(freq - exact_freq) / exact_freq);
	  CHECK(err <= max_err, "error exceeds tolerance");
     }
     printf("\n");

     }

     if (!stop1) {

     /*****************************************/

     printf("\nSolving for eigenvectors without preconditioning...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 op, op_data,
		 NULL, NULL,
		 maxwell_parity_constraint, (void *) mdata,
		 W, NWORK, error_tol, &num_iters, EIGS_DEFAULT_FLAGS);

     if (do_target)
	  eigensolver_get_eigenvals(H, eigvals, maxwell_operator, mdata,
				    W[0], W[1]);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("%15s%15s%15s%15s\n","eigenval", "frequency", "exact freq.", 
	    "error");
     for (i = 0; i < num_bands; ++i) {
	  double err;
	  real freq = sqrt(eigvals[i]);
	  real exact_freq = bragg_omega(freq, kvector[0], sqrt(ed.eps_high),
					ed.eps_high_x, sqrt(ed.eps_low),
					1.0 - ed.eps_high_x, 1.0e-7);
	  printf("%15f%15f%15f%15e\n", eigvals[i], freq, exact_freq,
		 err = fabs(freq - exact_freq) / exact_freq);
	  CHECK(err <= max_err, "error exceeds tolerance");
     }
     printf("\n");

     /*****************************************/
     
     printf("\nSolving for eigenvectors without conj. grad...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 op, op_data,
		 pre_op, pre_op_data,
		 maxwell_parity_constraint, (void *) mdata,
		 W, NWORK - 1, error_tol, &num_iters, EIGS_DEFAULT_FLAGS);

     if (do_target)
	  eigensolver_get_eigenvals(H, eigvals, maxwell_operator, mdata,
				    W[0], W[1]);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("%15s%15s%15s%15s\n","eigenval", "frequency", "exact freq.", 
	    "error");
     for (i = 0; i < num_bands; ++i) {
	  double err;
	  real freq = sqrt(eigvals[i]);
	  real exact_freq = bragg_omega(freq, kvector[0], sqrt(ed.eps_high),
					ed.eps_high_x, sqrt(ed.eps_low),
					1.0 - ed.eps_high_x, 1.0e-7);
	  printf("%15f%15f%15f%15e\n", eigvals[i], freq, exact_freq,
		 err = fabs(freq - exact_freq) / exact_freq);
	  CHECK(err <= max_err, "error exceeds tolerance");
     }
     printf("\n");

     /*****************************************/
     printf("\nSolving for eigenvectors without precond. or conj. grad...\n");
     evectmatrix_copy(H, Hstart);
     eigensolver(H, eigvals,
		 op, op_data,
		 NULL, NULL,
		 maxwell_parity_constraint, (void *) mdata,
		 W, NWORK - 1, error_tol, &num_iters, EIGS_DEFAULT_FLAGS);

     if (do_target)
	  eigensolver_get_eigenvals(H, eigvals, maxwell_operator, mdata,
				    W[0], W[1]);

     printf("Solved for eigenvectors after %d iterations.\n", num_iters);
     printf("%15s%15s%15s%15s\n","eigenval", "frequency", "exact freq.", 
	    "error");
     for (i = 0; i < num_bands; ++i) {
	  double err;
	  real freq = sqrt(eigvals[i]);
	  real exact_freq = bragg_omega(freq, kvector[0], sqrt(ed.eps_high),
					ed.eps_high_x, sqrt(ed.eps_low),
					1.0 - ed.eps_high_x, 1.0e-7);
	  printf("%15f%15f%15f%15e\n", eigvals[i], freq, exact_freq,
		 err = fabs(freq - exact_freq) / exact_freq);
	  CHECK(err <= max_err, "error exceeds tolerance");
     }
     printf("\n");

     /*****************************************/

     }
     
     destroy_evectmatrix(H);
     destroy_evectmatrix(Hstart);
     for (i = 0; i < NWORK; ++i)
          destroy_evectmatrix(W[i]);

     destroy_maxwell_target_data(mtdata);
     destroy_maxwell_data(mdata);

     free(eigvals);

     debug_check_memory_leaks();

     return EXIT_SUCCESS;
}
