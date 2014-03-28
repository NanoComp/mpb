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

/**************************************************************************/

/* Here, we define the external functions callable from Guile, as defined
   by mpb.scm. */

/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <stddef.h>

/* Header files for my eigensolver routines: */
#include "config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>
#include <maxwell.h>

/* Header file for the ctl-file (Guile) interface; automatically
   generated from mpb.scm */
#include <ctl-io.h>

/* Routines from libctl/utils/geom.c: */
#include <ctlgeom.h>

/* shared declarations for the files in mpb-ctl: */
#include "mpb.h"

/*  matrix/field "smobs" (Scheme objects) */
#include "matrix-smob.h"
#include "field-smob.h"

#if defined(DEBUG) && defined(HAVE_FEENABLEEXCEPT)
#  ifndef _GNU_SOURCE
#    define _GNU_SOURCE 1
#  endif
#  include <fenv.h>
#  if !HAVE_DECL_FEENABLEEXCEPT
int feenableexcept (int EXCEPTS);
#  endif
#endif

/**************************************************************************/

/* The following are hook functions called from main() when
   starting the program and just before exiting.   We use
   them to initialize MPI and OpenMP */

#if defined(SCALAR_SINGLE_PREC)
#  define FFTW(x) fftwf_ ## x
#elif defined(SCALAR_LONG_DOUBLE_PREC)
#  define FFTW(x) fftwl_ ## x
#else
#  define FFTW(x) fftw_ ## x
#endif

#ifdef USE_OPENMP
#  include <omp.h>
#  include <fftw3.h>
#endif

#if defined(HAVE_MPI) && (defined(HAVE_LIBFFTW3F_MPI) || defined(HAVE_LIBFFTW3L_MPI) || defined(HAVE_LIBFFTW3_MPI))
#  define HAVE_FFTW3_MPI
#  include <fftw3-mpi.h>
#endif

void ctl_start_hook(int *argc, char ***argv)
{
     MPI_Init(argc, argv);

#ifdef USE_OPENMP
     {
	  char *senthread = getenv("OMP_NUM_THREADS");
	  int i, nthread = senthread ? (atoi(senthread) > 0
					? atoi(senthread) : 1) : 1;
	  for (i = 0; i < *argc; ++i)
	       if (!strncasecmp("--nthread=", (*argv)[i], 10)) {
		    int j;
		    CHECK((nthread=atoi((*argv)[i]+10)) > 0,
			  "invalid argument for --nthread=...");
		    *argc -= 1;
		    for (j = i; j < *argc; ++j)
			 (*argv)[j] = (*argv)[j+1];
	       }
	  omp_set_num_threads(nthread);
	  CHECK(FFTW(init_threads)(), "error initializing threaded FFTW");
	  FFTW(plan_with_nthreads)(nthread);
     }
#endif

#ifdef HAVE_FFTW3_MPI
     FFTW(mpi_init)();
#endif

#ifdef HAVE_LIBCTL_QUIET
     {
	  extern int libctl_quiet;
	  libctl_quiet = !mpi_is_master();
     }
#endif

     /* disable: Guile garbage collector sometimes triggers an fp exception */
#if defined(DEBUG) && defined(HAVE_FEENABLEEXCEPT) && 0
     feenableexcept(FE_INVALID | FE_OVERFLOW); /* crash on NaN/overflow */
#endif
}

void ctl_stop_hook(void)
{
#ifdef HAVE_FFTW3_MPI
     FFTW(mpi_cleanup)();
#endif
     MPI_Finalize();
}

/* The following is a hook function called from main() when initializing
   Guile, which can export any additional symbols to Guile: */
void ctl_export_hook(void)
{
     register_matrix_smobs();
     register_field_smobs();
}

/**************************************************************************/

/* Some Guile-callable functions so that ctl files can know something
   about MPI. */

boolean mpi_is_masterp(void)
{
     return mpi_is_master();
}

boolean using_mpip(void)
{
#ifdef HAVE_MPI
     return 1;
#else
     return 0;
#endif
}

integer mpi_num_procs(void)
{
     int num_procs;
     MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
     return num_procs;
}

integer mpi_proc_index(void)
{
     int proc_num;
     MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);
     return proc_num;
}

number mpi_max(number num)
{
     double x = num, xmax;
     mpi_allreduce(&x, &xmax, 1, double, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     return xmax;
}

/**************************************************************************/

/* a couple of utilities to convert libctl data types to the data
   types of the eigensolver & maxwell routines: */

void vector3_to_arr(real arr[3], vector3 v)
{
     arr[0] = v.x;
     arr[1] = v.y;
     arr[2] = v.z;
}

void matrix3x3_to_arr(real arr[3][3], matrix3x3 m)
{
     vector3_to_arr(arr[0], m.c0);
     vector3_to_arr(arr[1], m.c1);
     vector3_to_arr(arr[2], m.c2);
}

/**************************************************************************/

/* global variables for retaining data about the eigenvectors between
   calls from Guile: */

int nwork_alloc = 0;

maxwell_data *mdata = NULL;
maxwell_target_data *mtdata = NULL;
evectmatrix H, W[MAX_NWORK], Hblock;

vector3 cur_kvector;
scalar_complex *curfield = NULL;
int curfield_band;
char curfield_type = '-';

void curfield_reset(void) { curfield = NULL; curfield_type = '-'; }

/* R[i]/G[i] are lattice/reciprocal-lattice vectors */
real R[3][3], G[3][3];
matrix3x3 Rm, Gm; /* same thing, but matrix3x3 */
real Vol; /* computational cell volume = |det Rm| */

/* index of current kpoint, for labeling output */
int kpoint_index = 0;

/**************************************************************************/

scalar_complex cnumber2cscalar(cnumber c)
{
     scalar_complex cs;
     CASSIGN_SCALAR(cs, cnumber_re(c), cnumber_im(c));
     return cs;
}

cnumber cscalar2cnumber(scalar_complex cs)
{
     return make_cnumber(CSCALAR_RE(cs), CSCALAR_IM(cs));
}

cvector3 cscalar32cvector3(const scalar_complex *cs)
{
     cvector3 v;
     v.x = cscalar2cnumber(cs[0]);
     v.y = cscalar2cnumber(cs[1]);
     v.z = cscalar2cnumber(cs[2]);
     return v;
}

void cvector32cscalar3(scalar_complex *cs, cvector3 v)
{
     cs[0] = cnumber2cscalar(v.x);
     cs[1] = cnumber2cscalar(v.y);
     cs[2] = cnumber2cscalar(v.z);
}

/**************************************************************************/

/* initialize the field to random numbers; should only be called
   after init-params.  (Guile-callable.) */
void randomize_fields(void)
{
     int i;

     if (!mdata)
	  return;
     mpi_one_printf("Initializing fields to random numbers...\n");
     for (i = 0; i < H.n * H.p; ++i) {
	  ASSIGN_SCALAR(H.data[i], rand() * 1.0 / RAND_MAX,
			rand() * 1.0 / RAND_MAX);
     }
}

/**************************************************************************/

/* Guile-callable functions for getting/setting the kpoint index. */

integer get_kpoint_index(void)
{
     return kpoint_index;
}

void set_kpoint_index(integer i)
{
     kpoint_index = i;
}

/**************************************************************************/

/* return a string describing the current parity, used for frequency
   and filename prefixes */
const char *parity_string(maxwell_data *d)
{
     static char s[128];
     strcpy(s, "");
     if (d->parity & EVEN_Z_PARITY)
	  strcat(s, (d->nz == 1) ? "te" : "zeven");
     else if (d->parity & ODD_Z_PARITY)
	  strcat(s, (d->nz == 1) ? "tm" : "zodd");
     if (d->parity & EVEN_Y_PARITY)
	  strcat(s, "yeven");
     else if (d->parity & ODD_Y_PARITY)
	  strcat(s, "yodd");
     return s;
}

/* Set the current parity to solve for. (init-params should have
   already been called.  (Guile-callable; see mpb.scm.) 

   p >= 0 means a bitwise OR of the various parity constants from
   maxwell.h (NO_PARITY, EVEN_Z_PARITY, etcetera).

   p = -1 means the parity of the previous call, 
       or NO_PARITY if this is the first call */

void set_parity(integer p)
{
     static int last_p = -2;  /* initialize to some non-value */

     if (!mdata) {
	  mpi_one_fprintf(stderr,
		  "init-params must be called before set-parity!\n");
	  return;
     }

     if (p == -1)
	  p = last_p < 0 ? NO_PARITY : last_p;

     set_maxwell_data_parity(mdata, p);
     CHECK(mdata->parity == p, "k vector incompatible with parity");
     mpi_one_printf("Solving for band polarization: %s.\n",
		    parity_string(mdata));

     last_p = p;
     set_kpoint_index(0);  /* reset index */
}

/**************************************************************************/

/* Guile-callable function: init-params, which initializes any data
   that we need for the eigenvalue calculation.  When this function
   is called, the input variables (the geometry, etcetera) have already
   been read into the global variables defined in ctl-io.h.  
   
   p is the parity to use for the coming calculation, although
   this can be changed by calling set-parity.  p is interpreted
   in the same way as for set-parity.

   If reset_fields is false, then any fields from a previous run are
   retained if they are of the same dimensions.  Otherwise, new
   fields are allocated and initialized to random numbers. */
void init_params(integer p, boolean reset_fields)
{
     int i, local_N, N_start, alloc_N;
     int nx, ny, nz;
     int have_old_fields = 0;
     int block_size;
     
     /* Output a bunch of stuff so that the user can see what we're
	doing and what we've read in. */
     
     mpi_one_printf("init-params: initializing eigensolver data\n");
#ifndef SCALAR_COMPLEX
     mpi_one_printf("  -- assuming INVERSION SYMMETRY in the geometry.\n");
#endif
     
     mpi_one_printf("Computing %d bands with %e tolerance.\n",
		    num_bands, tolerance);
     if (target_freq != 0.0)
	  mpi_one_printf("Target frequency is %g\n", target_freq);
     
     get_grid_size_n(&nx, &ny, &nz);

     {
	  int true_rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);
	  if (true_rank < dimensions)
	       dimensions = true_rank;
	  else if (true_rank > dimensions) {
	       mpi_one_fprintf(stderr, 
			       "WARNING: rank of grid is > dimensions.\n"
			       "         setting extra grid dims. to 1.\n");
	       /* force extra dims to be 1: */
	       if (dimensions <= 2)
		    nz = 1;
	       if (dimensions <= 1)
		    ny = 1;
	  }
     }

     mpi_one_printf("Working in %d dimensions.\n", dimensions);
     mpi_one_printf("Grid size is %d x %d x %d.\n", nx, ny, nz);

     if (eigensolver_block_size != 0 && eigensolver_block_size < num_bands) {
	  block_size = eigensolver_block_size;
	  if (block_size < 0) {
	       /* Guess a block_size near -block_size, chosen so that
		  all blocks are nearly equal in size: */
	       block_size = (num_bands - block_size - 1) / (-block_size);
	       block_size = (num_bands + block_size - 1) / block_size;
	  }
	  mpi_one_printf("Solving for %d bands at a time.\n", block_size);
     }
     else
	  block_size = num_bands;

     if (mdata) {  /* need to clean up from previous init_params call */
	  if (nx == mdata->nx && ny == mdata->ny && nz == mdata->nz &&
	      block_size == Hblock.alloc_p && num_bands == H.p &&
	      eigensolver_nwork == nwork_alloc)
	       have_old_fields = 1; /* don't need to reallocate */
	  else {
	       destroy_evectmatrix(H);
	       for (i = 0; i < nwork_alloc; ++i)
		    destroy_evectmatrix(W[i]);
	       if (Hblock.data != H.data)
		    destroy_evectmatrix(Hblock);
	  }
	  destroy_maxwell_target_data(mtdata); mtdata = NULL;
	  destroy_maxwell_data(mdata); mdata = NULL;
	  curfield_reset();
     }
     else
	  srand(time(NULL)); /* init random seed for field initialization */
   
     if (deterministicp) {  /* check input variable "deterministic?" */
	  /* seed should be the same for each run, although
	     it should be different for each process: */
	  int rank;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  srand(314159 * (rank + 1));
     }

     mpi_one_printf("Creating Maxwell data...\n");
     mdata = create_maxwell_data(nx, ny, nz, &local_N, &N_start, &alloc_N,
                                 block_size, NUM_FFT_BANDS);
     CHECK(mdata, "NULL mdata");

     if (target_freq != 0.0)
	  mtdata = create_maxwell_target_data(mdata, target_freq);
     else
	  mtdata = NULL;

     if (!have_old_fields) {
	  mpi_one_printf("Allocating fields...\n");
	  H = create_evectmatrix(nx * ny * nz, 2, num_bands,
				 local_N, N_start, alloc_N);
	  nwork_alloc = eigensolver_nwork;
	  for (i = 0; i < nwork_alloc; ++i)
	       W[i] = create_evectmatrix(nx * ny * nz, 2, block_size,
					 local_N, N_start, alloc_N);
	  if (block_size < num_bands)
	       Hblock = create_evectmatrix(nx * ny * nz, 2, block_size,
					   local_N, N_start, alloc_N);
	  else
	       Hblock = H;
     }

     init_epsilon();

     mpi_one_printf("%d k-points:\n", k_points.num_items);
     for (i = 0; i < k_points.num_items; ++i)
	  mpi_one_printf("     (%g,%g,%g)\n", k_points.items[i].x,
			 k_points.items[i].y, k_points.items[i].z);

     set_parity(p);
     if (!have_old_fields || reset_fields)
	  randomize_fields();

     {
	  int ierr = check_maxwell_dielectric(mdata, negative_epsilon_okp);
	  if (ierr == 1)
	       mpi_one_fprintf(stderr,
			   "ERROR: non positive-definite dielectric tensor\n");
	  else if (ierr == 2)
	       mpi_one_fprintf(stderr, 
		       "ERROR: dielectric tensor must not couple xy "
		       "plane with z direction for 2D TE/TM calculations\n");
	  CHECK(!ierr, "invalid dielectric function\n");
     }

     evectmatrix_flops = eigensolver_flops; /* reset, if changed */
}

/**************************************************************************/

/* When we are solving for a few bands at a time, we solve for the
   upper bands by "deflation"--by continually orthogonalizing them
   against the already-computed lower bands.  (This constraint
   commutes with the eigen-operator, of course, so all is well.) */

typedef struct {
     evectmatrix Y;  /* the vectors to orthogonalize against; Y must
			itself be normalized (Yt Y = 1) */
     int p;  /* the number of columns of Y to orthogonalize against */
     scalar *S;  /* a matrix for storing the dot products; should have
		    at least p * X.p elements (see below for X) */
     scalar *S2; /* a scratch matrix the same size as S */
} deflation_data;

static void deflation_constraint(evectmatrix X, void *data)
{
     deflation_data *d = (deflation_data *) data;

     CHECK(X.n == d->Y.n && d->Y.p >= d->p, "invalid dimensions");

     /* (Sigh...call the BLAS functions directly since we are not
	using all the columns of Y...evectmatrix is not set up for
	this case.) */

     /* compute S = Xt Y (i.e. all the dot products): */
     blasglue_gemm('C', 'N', X.p, d->p, X.n,
		   1.0, X.data, X.p, d->Y.data, d->Y.p, 0.0, d->S2, d->p);
     mpi_allreduce(d->S2, d->S, d->p * X.p * SCALAR_NUMVALS,
		   real, SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);

     /* compute X = X - Y*St = (1 - Y Yt) X */
     blasglue_gemm('N', 'C', X.n, X.p, d->p,
		   -1.0, d->Y.data, d->Y.p, d->S, d->p,
		   1.0, X.data, X.p);
}

/**************************************************************************/

/* Solve for the bands at a given k point.
   Must only be called after init_params! */
void solve_kpoint(vector3 kvector)
{
     int i, total_iters = 0, ib, ib0;
     real *eigvals;
     real k[3];
     int flags;
     deflation_data deflation;
     int prev_parity;

     /* if we get too close to singular k==0 point, just set k=0
	to exploit our special handling of this k */
     if (vector3_norm(kvector) < 1e-10)
	  kvector.x = kvector.y = kvector.z = 0;

     mpi_one_printf("solve_kpoint (%g,%g,%g):\n",
		    kvector.x, kvector.y, kvector.z);
     
     curfield_reset();

     if (num_bands == 0) {
	  mpi_one_printf("  num-bands is zero, not solving for any bands\n");
	  return;
     }

     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before solve-kpoint!\n");
	  return;
     }

     /* if this is the first k point, print out a header line for
	for the frequency grep data: */
     if (!kpoint_index && mpi_is_master()) {
	  printf("%sfreqs:, k index, k1, k2, k3, kmag/2pi",
		 parity_string(mdata));
	  for (i = 0; i < num_bands; ++i)
	       printf(", %s%sband %d",
		      parity_string(mdata),
		      mdata->parity == NO_PARITY ? "" : " ",
		      i + 1);
	  printf("\n");
     }

     prev_parity = mdata->parity;
     cur_kvector = kvector;
     vector3_to_arr(k, kvector);
     update_maxwell_data_k(mdata, k, G[0], G[1], G[2]);
     CHECK(mdata->parity == prev_parity,
	   "k vector is incompatible with specified parity");

     CHK_MALLOC(eigvals, real, num_bands);

     flags = eigensolver_flags; /* ctl file input variable */
     if (verbose)
	  flags |= EIGS_VERBOSE;

     /* constant (zero frequency) bands at k=0 are handled specially,
        so remove them from the solutions for the eigensolver: */
     if (mdata->zero_k && !mtdata) {
	  int in, ip;
	  ib0 = maxwell_zero_k_num_const_bands(H, mdata);
	  for (in = 0; in < H.n; ++in)
	       for (ip = 0; ip < H.p - ib0; ++ip)
		    H.data[in * H.p + ip] = H.data[in * H.p + ip + ib0];
	  evectmatrix_resize(&H, H.p - ib0, 1);
     }
     else
	  ib0 = 0; /* solve for all bands */

     /* Set up deflation data: */
     if (H.data != Hblock.data) {
	  deflation.Y = H;
	  deflation.p = 0;
	  CHK_MALLOC(deflation.S, scalar, H.p * Hblock.p);
	  CHK_MALLOC(deflation.S2, scalar, H.p * Hblock.p);
     }

     for (ib = ib0; ib < num_bands; ib += Hblock.alloc_p) {
	  evectconstraint_chain *constraints;
	  int num_iters;

	  /* don't solve for too many bands if the block size doesn't divide
	     the number of bands: */
	  if (ib + mdata->num_bands > num_bands) {
	       maxwell_set_num_bands(mdata, num_bands - ib);
	       for (i = 0; i < nwork_alloc; ++i)
		    evectmatrix_resize(&W[i], num_bands - ib, 0);
	       evectmatrix_resize(&Hblock, num_bands - ib, 0);
	  }

	  mpi_one_printf("Solving for bands %d to %d...\n",
			 ib + 1, ib + Hblock.p);

	  constraints = NULL;
	  constraints = evect_add_constraint(constraints,
					     maxwell_parity_constraint,
					     (void *) mdata);

	  if (mdata->zero_k)
	       constraints = evect_add_constraint(constraints,
						  maxwell_zero_k_constraint,
						  (void *) mdata);

	  if (Hblock.data != H.data) {  /* initialize fields of block from H */
	       int in, ip;
	       for (in = 0; in < Hblock.n; ++in)
		    for (ip = 0; ip < Hblock.p; ++ip)
			 Hblock.data[in * Hblock.p + ip] =
			      H.data[in * H.p + ip + (ib-ib0)];
	       deflation.p = ib-ib0;
	       if (deflation.p > 0)
		    constraints = evect_add_constraint(constraints,
						       deflation_constraint,
						       &deflation);
	  }

	  if (mtdata) {  /* solving for bands near a target frequency */
	       if (eigensolver_davidsonp)
		    eigensolver_davidson(
			 Hblock, eigvals + ib,
			 maxwell_target_operator, (void *) mtdata,
			 simple_preconditionerp ? 
			 maxwell_target_preconditioner :
			 maxwell_target_preconditioner2,
			 (void *) mtdata,
			 evectconstraint_chain_func,
			 (void *) constraints,
			 W, nwork_alloc, tolerance, &num_iters, flags, 0.0);
	       else
		    eigensolver(Hblock, eigvals + ib,
				maxwell_target_operator, (void *) mtdata,
				simple_preconditionerp ? 
				maxwell_target_preconditioner :
				maxwell_target_preconditioner2,
				(void *) mtdata,
				evectconstraint_chain_func,
				(void *) constraints,
				W, nwork_alloc, tolerance, &num_iters, flags);

	       /* now, diagonalize the real Maxwell operator in the
		  solution subspace to get the true eigenvalues and
		  eigenvectors: */
	       CHECK(nwork_alloc >= 2, "not enough workspace");
	       eigensolver_get_eigenvals(Hblock, eigvals + ib,
					 maxwell_operator,mdata, W[0],W[1]);
	  }
	  else {
	       if (eigensolver_davidsonp)
		    eigensolver_davidson(
			 Hblock, eigvals + ib,
			 maxwell_operator, (void *) mdata,
			 simple_preconditionerp ?
			 maxwell_preconditioner :
			 maxwell_preconditioner2,
			 (void *) mdata,
			 evectconstraint_chain_func,
			 (void *) constraints,
			 W, nwork_alloc, tolerance, &num_iters, flags, 0.0);
	       else
		    eigensolver(Hblock, eigvals + ib,
				maxwell_operator, (void *) mdata,
				simple_preconditionerp ?
				maxwell_preconditioner :
				maxwell_preconditioner2,
				(void *) mdata,
				evectconstraint_chain_func,
				(void *) constraints,
				W, nwork_alloc, tolerance, &num_iters, flags);
	  }
	  
	  if (Hblock.data != H.data) {  /* save solutions of current block */
	       int in, ip;
	       for (in = 0; in < Hblock.n; ++in)
		    for (ip = 0; ip < Hblock.p; ++ip)
			 H.data[in * H.p + ip + (ib-ib0)] =
			      Hblock.data[in * Hblock.p + ip];
	  }

	  evect_destroy_constraints(constraints);
	  
	  mpi_one_printf("Finished solving for bands %d to %d after "
			 "%d iterations.\n", ib + 1, ib + Hblock.p, num_iters);
	  total_iters += num_iters * Hblock.p;
     }

     if (num_bands - ib0 > Hblock.alloc_p)
	  mpi_one_printf("Finished k-point with %g mean iterations/band.\n",
			 total_iters * 1.0 / num_bands);

     /* Manually put in constant (zero-frequency) solutions for k=0: */
     if (mdata->zero_k && !mtdata) {
	  int in, ip;
	  evectmatrix_resize(&H, H.alloc_p, 1);
	  for (in = 0; in < H.n; ++in)
	       for (ip = H.p - ib0 - 1; ip >= 0; --ip)
		    H.data[in * H.p + ip + ib0] = H.data[in * H.p + ip];
	  maxwell_zero_k_set_const_bands(H, mdata);
	  for (ib = 0; ib < ib0; ++ib)
	       eigvals[ib] = 0;
     }

     /* Reset scratch matrix sizes: */
     evectmatrix_resize(&Hblock, Hblock.alloc_p, 0);
     for (i = 0; i < nwork_alloc; ++i)
	  evectmatrix_resize(&W[i], W[i].alloc_p, 0);
     maxwell_set_num_bands(mdata, Hblock.alloc_p);

     /* Destroy deflation data: */
     if (H.data != Hblock.data) {
	  free(deflation.S2);
	  free(deflation.S);
     }

     if (num_write_output_vars > 1) {
	  /* clean up from prev. call */
	  free(freqs.items);
	  free(parity);
     }

     CHK_MALLOC(parity, char, strlen(parity_string(mdata)) + 1);
     parity = strcpy(parity, parity_string(mdata));

     iterations = total_iters; /* iterations output variable */

     /* create freqs array for storing frequencies in a Guile list */
     freqs.num_items = num_bands;
     CHK_MALLOC(freqs.items, number, freqs.num_items);
     
     set_kpoint_index(kpoint_index + 1);

     mpi_one_printf("%sfreqs:, %d, %g, %g, %g, %g",
		    parity,
		    kpoint_index, (double)k[0], (double)k[1], (double)k[2],
		    vector3_norm(matrix3x3_vector3_mult(Gm, kvector)));
     for (i = 0; i < num_bands; ++i) {
	  freqs.items[i] =
	       negative_epsilon_okp ? eigvals[i] : sqrt(eigvals[i]);
	  mpi_one_printf(", %g", freqs.items[i]);
     }
     mpi_one_printf("\n");

     eigensolver_flops = evectmatrix_flops;

     free(eigvals);
}

/**************************************************************************/

/* Return a list of the z/y parities, one for each band. */

number_list compute_zparities(void)
{
     number_list z_parity;
     z_parity.num_items = num_bands;
     z_parity.items = maxwell_zparity(H, mdata);
     return z_parity;
}

number_list compute_yparities(void)
{
     number_list y_parity;
     y_parity.num_items = num_bands;
     y_parity.items = maxwell_yparity(H, mdata);
     return y_parity;
}

/**************************************************************************/

/* Compute the group velocity dw/dk in the given direction d (where
   the length of d is ignored).  d is in the reciprocal lattice basis.
   Should only be called after solve_kpoint.  Returns a list of the
   group velocities, one for each band, in units of c. */
number_list compute_group_velocity_component(vector3 d)
{
     number_list group_v;
     real *gv_scratch;
     real u[3];
     int i, ib;

     group_v.num_items = 0;  group_v.items = (number *) NULL;

     curfield_reset(); /* has the side effect of overwriting curfield scratch */

     if (!mdata) {
	  mpi_one_fprintf(stderr, "init-params must be called first!\n");
	  return group_v;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr, "solve-kpoint must be called first!\n");
	  return group_v;
     }

     /* convert d to unit vector in Cartesian coords: */
     d = unit_vector3(matrix3x3_vector3_mult(Gm, d));
     u[0] = d.x; u[1] = d.y; u[2] = d.z;

     group_v.num_items = num_bands;
     CHK_MALLOC(group_v.items, number, group_v.num_items);
     CHK_MALLOC(gv_scratch, real, group_v.num_items * 2);
     
     /* now, compute group_v.items = diag Re <H| curl 1/eps i u x |H>: */

     /* ...we have to do this in blocks of eigensolver_block_size since
	the work matrix W[0] may not have enough space to do it all at once. */
     
     for (ib = 0; ib < num_bands; ib += Hblock.alloc_p) {
	  if (ib + mdata->num_bands > num_bands) {
	       maxwell_set_num_bands(mdata, num_bands - ib);
	       evectmatrix_resize(&W[0], num_bands - ib, 0);
	       evectmatrix_resize(&Hblock, num_bands - ib, 0);
	  }
	  if (Hblock.data != H.data) {  /* initialize fields of block from H */
	       int in, ip;
	       for (in = 0; in < Hblock.n; ++in)
		    for (ip = 0; ip < Hblock.p; ++ip)
			 Hblock.data[in * Hblock.p + ip] =
			      H.data[in * H.p + ip + ib];
	  }

	  maxwell_ucross_op(Hblock, W[0], mdata, u);
	  evectmatrix_XtY_diag_real(Hblock, W[0], gv_scratch,
				    gv_scratch + group_v.num_items);
	  {
	       int ip;
	       for (ip = 0; ip < Hblock.p; ++ip)
		    group_v.items[ib + ip] = gv_scratch[ip];
	  }
     }

     free(gv_scratch);

     /* Reset scratch matrix sizes: */
     evectmatrix_resize(&Hblock, Hblock.alloc_p, 0);
     evectmatrix_resize(&W[0], W[0].alloc_p, 0);
     maxwell_set_num_bands(mdata, Hblock.alloc_p);

     /* The group velocity is given by:

	grad_k(omega)*d = grad_k(omega^2)*d / 2*omega
	   = grad_k(<H|maxwell_op|H>)*d / 2*omega
	   = Re <H| curl 1/eps i u x |H> / omega
        
        Note that our k is in units of 2*Pi/a, and omega is in
        units of 2*Pi*c/a, so the result will be in units of c. */
     for (i = 0; i < num_bands; ++i) {
	  if (freqs.items[i] == 0)  /* v is undefined in this case */
	       group_v.items[i] = 0.0;  /* just set to zero */
	  else
	       group_v.items[i] /= 
		    negative_epsilon_okp ? sqrt(fabs(freqs.items[i]))
		    : freqs.items[i];
     }
     
     return group_v;
}

/* as above, but only computes for given band */
number compute_1_group_velocity_component(vector3 d, integer b)
{
     real u[3];
     int ib = b - 1;
     real group_v, scratch;

     curfield_reset();

     if (!mdata) {
	  mpi_one_fprintf(stderr, "init-params must be called first!\n");
	  return group_v;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr, "solve-kpoint must be called first!\n");
	  return group_v;
     }

     /* convert d to unit vector in Cartesian coords: */
     d = unit_vector3(matrix3x3_vector3_mult(Gm, d));
     u[0] = d.x; u[1] = d.y; u[2] = d.z;

     evectmatrix_resize(&W[0], 1, 0);
     CHECK(nwork_alloc > 1, "eigensolver-nwork is too small");
     evectmatrix_resize(&W[1], 1, 0);

     {  /* initialize fields of block from H */
	  int in;
	  scalar *data = W[1].data;
	  for (in = 0; in < W[1].n; ++in)
	       data[in] = H.data[in * H.p + ib];
     }

     maxwell_ucross_op(W[1], W[0], mdata, u);
     evectmatrix_XtY_diag_real(W[1], W[0], &group_v, &scratch);

     /* Reset scratch matrix sizes: */
     evectmatrix_resize(&W[1], W[1].alloc_p, 0);
     evectmatrix_resize(&W[0], W[0].alloc_p, 0);

     if (freqs.items[ib] == 0)  /* v is undefined in this case */
	  group_v = 0.0;  /* just set to zero */
     else
	  group_v /= negative_epsilon_okp ? sqrt(fabs(freqs.items[ib]))
	       : freqs.items[ib];

     return group_v;
}

/* returns group velocity for band b, in Cartesian coordinates */
vector3 compute_1_group_velocity(integer b)
{
     vector3 v, d;
     matrix3x3 RmT = matrix3x3_transpose(Rm);
     d.x = 1; d.y = d.z = 0;
     v.x = compute_1_group_velocity_component(matrix3x3_vector3_mult(RmT,d),b);
     d.y = 1; d.x = d.z = 0;
     v.y = compute_1_group_velocity_component(matrix3x3_vector3_mult(RmT,d),b);
     d.z = 1; d.y = d.x = 0;
     v.z = compute_1_group_velocity_component(matrix3x3_vector3_mult(RmT,d),b);
     return v;
}

/* as above, but returns "group velocity" given by gradient of
   frequency with respect to k in reciprocal coords ... this is useful
   for band optimization. */
vector3 compute_1_group_velocity_reciprocal(integer b)
{
     return matrix3x3_vector3_mult(matrix3x3_transpose(Gm), 
				   compute_1_group_velocity(b));
}

/**************************************************************************/
