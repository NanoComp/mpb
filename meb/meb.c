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

/**************************************************************************/

/* Here, we define the external functions callable from Guile, as defined
   by meb.scm. */

/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <stddef.h>

/* GNU Guile library header file: */
#include <guile/gh.h>

/* Header files for my eigensolver routines: */
#include "../src/config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>

/* Header file for the ctl-file (Guile) interface; automatically
   generated from meb.scm */
#include <ctl-io.h>

/* Routines from libctl/utils/geom.c: */
#include <ctlgeom.h>

/* Functions for the elastic eigenproblem, analogous to src/maxwell */
#include "elastic.h"

/* shared declarations for the files in meb-ctl: */
#include "meb.h"

/**************************************************************************/

/* force use of Guile garbage-collector when we run low */
static void *malloc_hook(size_t n)
{
     void *p = malloc(n);
     if (!p && n) { /* garbage-collect, then try again */
          scm_gc();
          p = malloc(n);
     }
     CHECK(p || !n, "out of memory");
     return p;
}

/* The following are hook functions called from main() when
   starting the program and just before exiting.   We use
   them to initialize MPI. */
#include <fenv.h>
void ctl_start_hook(int *argc, char ***argv)
{
     my_malloc_hook = malloc_hook;
     MPI_Init(argc, argv);
     feenableexcept(FE_INVALID | FE_OVERFLOW);
}

void ctl_stop_hook(void)
{
     MPI_Finalize();
}

/* The following is a hook function called from main() when initializing
   Guile, which can export any additional symbols to Guile: */
void ctl_export_hook(void)
{
     /* register_matrix_smobs(); */
     /* register_field_smobs(); */
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

elastic_data *edata = NULL;
evectmatrix V, W[MAX_NWORK], Vblock;

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

     if (!edata)
	  return;
     mpi_one_printf("Initializing fields to random numbers...\n");
     for (i = 0; i < V.n * V.p; ++i) {
	  ASSIGN_SCALAR(V.data[i], rand() * 1.0 / RAND_MAX,
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
const char *parity_string(elastic_data *d)
{
     static char s[128];
     (void) d; /* unused, for now */
     strcpy(s, "");
     return s;
}

/* Set the current parity to solve for. (init-params should have
   already been called.  (Guile-callable; see meb.scm.) 

   p >= 0 means a bitwise OR of the various parity constants from
   elastic.h (NO_PARITY, EVEN_Z_PARITY, etcetera).

   p = -1 means the parity of the previous call, 
       or NO_PARITY if this is the first call */

void set_parity(integer p)
{
     static int last_p = -2;  /* initialize to some non-value */

     if (!edata) {
	  mpi_one_fprintf(stderr,
		  "init-params must be called before set-parity!\n");
	  return;
     }

     if (p == -1)
	  p = last_p < 0 ? NO_PARITY : last_p;

     set_elastic_data_parity(edata, p);
     CHECK(edata->parity == p, "k vector incompatible with parity");
     mpi_one_printf("Solving for band polarization: %s.\n",
		    parity_string(edata));

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

     if (edata) {  /* need to clean up from previous init_params call */
	  if (nx == edata->nx && ny == edata->ny && nz == edata->nz &&
	      block_size == Vblock.alloc_p && num_bands == V.p &&
	      eigensolver_nwork == nwork_alloc)
	       have_old_fields = 1; /* don't need to reallocate */
	  else {
	       destroy_evectmatrix(V);
	       for (i = 0; i < nwork_alloc; ++i)
		    destroy_evectmatrix(W[i]);
	       if (Vblock.data != V.data)
		    destroy_evectmatrix(Vblock);
	  }
	  destroy_elastic_data(edata); edata = NULL;
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

     mpi_one_printf("Creating Elastic data...\n");
     edata = create_elastic_data(nx, ny, nz, &local_N, &N_start, &alloc_N,
                                 block_size, NUM_FFT_BANDS);
     CHECK(edata, "NULL edata");

     if (!have_old_fields) {
	  mpi_one_printf("Allocating fields...\n");
	  V = create_evectmatrix(nx * ny * nz, 3, num_bands,
				 local_N, N_start, alloc_N);
	  nwork_alloc = eigensolver_nwork;
	  for (i = 0; i < nwork_alloc; ++i)
	       W[i] = create_evectmatrix(nx * ny * nz, 3, block_size,
					 local_N, N_start, alloc_N);
	  if (block_size < num_bands)
	       Vblock = create_evectmatrix(nx * ny * nz, 3, block_size,
					   local_N, N_start, alloc_N);
	  else
	       Vblock = V;
     }

     init_elastic_materials();

     mpi_one_printf("%d k-points:\n", k_points.num_items);
     for (i = 0; i < k_points.num_items; ++i)
	  mpi_one_printf("     (%g,%g,%g)\n", k_points.items[i].x,
			 k_points.items[i].y, k_points.items[i].z);

     set_parity(p);
     if (!have_old_fields || reset_fields)
	  randomize_fields();

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

     mpi_one_printf("solve_kpoint (%g,%g,%g):\n",
		    kvector.x, kvector.y, kvector.z);
     
     curfield_reset();

     if (num_bands == 0) {
	  mpi_one_printf("  num-bands is zero, not solving for any bands\n");
	  return;
     }

     if (!edata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before solve-kpoint!\n");
	  return;
     }

     /* if this is the first k point, print out a header line for
	for the frequency grep data: */
     if (!kpoint_index && mpi_is_master()) {
	  printf("%sfreqs:, k index, k1, k2, k3, kmag/2pi",
		 parity_string(edata));
	  for (i = 0; i < num_bands; ++i)
	       printf(", %s%sband %d",
		      parity_string(edata),
		      edata->parity == NO_PARITY ? "" : " ",
		      i + 1);
	  printf("\n");
     }

     prev_parity = edata->parity;
     cur_kvector = kvector;
     vector3_to_arr(k, kvector);
     update_elastic_data_k(edata, k, G[0], G[1], G[2]);
     CHECK(edata->parity == prev_parity,
	   "k vector is incompatible with specified parity");

     CHK_MALLOC(eigvals, real, num_bands);

     flags = eigensolver_flags; /* ctl file input variable */
     if (verbose)
	  flags |= EIGS_VERBOSE;

     /* constant-u (zero frequency) bands at k=0 are handled
        specially, so we should remove them from the solutions for the
        eigensolver: */
     if (edata->zero_k) {
          int in, ip;
          ib0 = elastic_zero_k_num_const_bands(V, edata);
	  elastic_zero_k_set_const_bands(V, edata);
	  for (ib = 0; ib < ib0; ++ib)
               eigvals[ib] = 0;
     }
     else
	  ib0 = 0; /* solve for all bands */

     /* Set up deflation data: */
     if (V.data != Vblock.data) {
	  deflation.Y = V;
	  deflation.p = ib0;
	  CHK_MALLOC(deflation.S, scalar, V.p * Vblock.p);
	  CHK_MALLOC(deflation.S2, scalar, V.p * Vblock.p);
     }

     for (ib = ib0; ib < num_bands; ib += Vblock.alloc_p) {
	  evectconstraint_chain *constraints;
	  int num_iters;

	  /* don't solve for too many bands if the block size doesn't divide
	     the number of bands: */
	  if (ib + edata->num_bands > num_bands) {
	       elastic_set_num_bands(edata, num_bands - ib);
	       for (i = 0; i < nwork_alloc; ++i)
		    evectmatrix_resize(&W[i], num_bands - ib, 0);
	       evectmatrix_resize(&Vblock, num_bands - ib, 0);
	  }

	  mpi_one_printf("Solving for bands %d to %d...\n",
			 ib + 1, ib + Vblock.p);

	  constraints = NULL;
	  constraints = evect_add_constraint(constraints,
					     elastic_parity_constraint,
					     (void *) edata);

	  if (Vblock.data != V.data) {  /* initialize fields of block from V */
	       int in, ip;
	       for (in = 0; in < Vblock.n; ++in)
		    for (ip = 0; ip < Vblock.p; ++ip)
			 Vblock.data[in * Vblock.p + ip] =
			      V.data[in * V.p + ip + ib];
	       deflation.p = ib;
	       if (deflation.p > 0)
		    constraints = evect_add_constraint(constraints,
						       deflation_constraint,
						       &deflation);
	  }

	  if (eigensolver_davidsonp)
	       eigensolver_davidson(
		    Vblock, eigvals + ib,
		    elastic_operator, (void *) edata,
		    NULL, (void *) edata,
		    evectconstraint_chain_func,
		    (void *) constraints,
		    W, nwork_alloc, tolerance, &num_iters, flags, 0.0);
	  else
	       eigensolver(Vblock, eigvals + ib,
			   elastic_operator, (void *) edata,
			   NULL, (void *) edata,
			   evectconstraint_chain_func,
			   (void *) constraints,
			   W, nwork_alloc, tolerance, &num_iters, flags);
	  
	  if (Vblock.data != V.data) {  /* save solutions of current block */
	       int in, ip;
	       for (in = 0; in < Vblock.n; ++in)
		    for (ip = 0; ip < Vblock.p; ++ip)
			 V.data[in * V.p + ip + ib] =
			      Vblock.data[in * Vblock.p + ip];
	  }

	  evect_destroy_constraints(constraints);
	  
	  mpi_one_printf("Finished solving for bands %d to %d after "
			 "%d iterations.\n", ib + 1, ib + Vblock.p, num_iters);
	  total_iters += num_iters * Vblock.p;
     }

     if (num_bands - ib0 > Vblock.alloc_p)
	  mpi_one_printf("Finished k-point with %g mean iterations/band.\n",
			 total_iters * 1.0 / num_bands);

     /* Reset scratch matrix sizes: */
     evectmatrix_resize(&Vblock, Vblock.alloc_p, 0);
     for (i = 0; i < nwork_alloc; ++i)
	  evectmatrix_resize(&W[i], W[i].alloc_p, 0);
     elastic_set_num_bands(edata, Vblock.alloc_p);

     /* Destroy deflation data: */
     if (V.data != Vblock.data) {
	  free(deflation.S2);
	  free(deflation.S);
     }

     if (num_write_output_vars > 1) {
	  /* clean up from prev. call */
	  free(freqs.items);
	  free(parity);
     }

     CHK_MALLOC(parity, char, strlen(parity_string(edata)) + 1);
     parity = strcpy(parity, parity_string(edata));

     iterations = total_iters; /* iterations output variable */

     /* create freqs array for storing frequencies in a Guile list */
     freqs.num_items = num_bands;
     CHK_MALLOC(freqs.items, number, freqs.num_items);
     
     set_kpoint_index(kpoint_index + 1);

     mpi_one_printf("%sfreqs:, %d, %g, %g, %g, %g",
		    parity,
		    kpoint_index, k[0], k[1], k[2],
		    vector3_norm(matrix3x3_vector3_mult(Gm, kvector)));
     for (i = 0; i < num_bands; ++i) {
	  freqs.items[i] = sqrt(eigvals[i]);
	  mpi_one_printf(", %g", freqs.items[i]);
     }
     mpi_one_printf("\n");

     eigensolver_flops = evectmatrix_flops;

     free(eigvals);
}

/**************************************************************************/


