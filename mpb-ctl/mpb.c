/* Copyright (C) 1999, 2000, 2001 Massachusetts Institute of Technology.
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

/* GNU Guile library header file: */
#include <guile/gh.h>

/* Header files for my eigensolver routines: */
#include "../src/config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <matrixio.h>
#include <eigensolver.h>
#include <maxwell.h>

/* Header file for the ctl-file (Guile) interface; automatically
   generated from mpb.scm */
#include <ctl-io.h>

/* Routines from libctl/utils/geom.c: */
#include <ctlgeom.h>

/* shared declarations for the files in mpb-ctl: */
#include "mpb.h"

/* this integer flag is defined by main.c from libctl, and is
   set when the user runs the program with --verbose */
extern int verbose;

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))

#define TWOPI 6.2831853071795864769252867665590057683943388

/**************************************************************************/

/* The following are hook functions called from main() when
   starting the program and just before exiting.   We use
   them to initialize MPI. */

void ctl_start_hook(int *argc, char ***argv)
{
     MPI_Init(argc, argv);
}

void ctl_stop_hook(void)
{
     MPI_Finalize();
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

static void vector3_to_arr(real arr[3], vector3 v)
{
     arr[0] = v.x;
     arr[1] = v.y;
     arr[2] = v.z;
}

static void matrix3x3_to_arr(real arr[3][3], matrix3x3 m)
{
     vector3_to_arr(arr[0], m.c0);
     vector3_to_arr(arr[1], m.c1);
     vector3_to_arr(arr[2], m.c2);
}

/**************************************************************************/

#define MAX_NWORK 10
static int nwork_alloc = 0;

#define NUM_FFT_BANDS 20 /* max number of bands to FFT at a time */

/* global variables for retaining data about the eigenvectors between
   calls from Guile: */

static maxwell_data *mdata = NULL;
static maxwell_target_data *mtdata = NULL;
static evectmatrix H, W[MAX_NWORK], Hblock;

static vector3 cur_kvector;
static scalar_complex *curfield = NULL;
static int curfield_band;
static char curfield_type = '-';

static void curfield_reset(void) { curfield = NULL; curfield_type = '-'; }

/* R[i]/G[i] are lattice/reciprocal-lattice vectors */
static real R[3][3], G[3][3];
static matrix3x3 Rm, Gm; /* same thing, but matrix3x3 */

/* index of current kpoint, for labeling output */
static int kpoint_index = 0;

#define USE_GEOMETRY_TREE 1
geom_box_tree geometry_tree = NULL; /* recursive tree of geometry objects
				       for fast searching */

/**************************************************************************/

typedef struct {
     maxwell_dielectric_function eps_file_func;
     void *eps_file_func_data;
} epsilon_func_data;

/* Given a position r in the basis of the lattice vectors, return the
   corresponding dielectric tensor and its inverse.  Should be
   called from within init_params (or after init_params), so that the
   geometry input variables will have been read in (for use in libgeom).

   This function is passed to set_maxwell_dielectric to initialize
   the dielectric tensor array for eigenvector calculations. */

static void epsilon_func(symmetric_matrix *eps, symmetric_matrix *eps_inv,
			 real r[3], void *edata)
{
     epsilon_func_data *d = (epsilon_func_data *) edata;
     material_type material;
     vector3 p;
     boolean inobject;

     /* p needs to be in the lattice *unit* vector basis, while r is
	in the lattice vector basis.  Also, shift origin to the center
        of the grid. */
     p.x = (r[0] - 0.5) * geometry_lattice.size.x;
     p.y = dimensions <= 1 ? 0 : (r[1] - 0.5) * geometry_lattice.size.y;
     p.z = dimensions <= 2 ? 0 : (r[2] - 0.5) * geometry_lattice.size.z;

     /* call search routine from libctl/utils/libgeom/geom.c: */
#if USE_GEOMETRY_TREE
     material = material_of_point_in_tree_inobject(p, geometry_tree, 
						   &inobject);
#else
     material = material_of_point_inobject(p, &inobject);
#endif

#if defined(DEBUG_GEOMETRY_TREE) && USE_GEOMETRY_TREE
     {
	  material_type m2 = material_of_point_inobject(p, &inobject);
	  CHECK(m2.which_subclass == material.which_subclass &&
		m2.subclass.dielectric_data ==
		material.subclass.dielectric_data,
		"material_of_point & material_of_point_in_tree don't agree!");
     }
#endif

     if (material.which_subclass == MATERIAL_TYPE_SELF) {
	  material = default_material;
	  inobject = 0;  /* treat as a "nothing" object */
     }

     /* if we aren't in any geometric object and we have an epsilon
	file, use that. */
     if (!inobject && d->eps_file_func) {
	  d->eps_file_func(eps, eps_inv, r, d->eps_file_func_data);
     }
     else {
	  boolean destroy_material = 0;

	  while (material.which_subclass == MATERIAL_FUNCTION) {
	       material_type m;
	       SCM mo;
	       /* material_func is a Scheme function, taking a position
		  vector and returning a material at that point: */
	       mo = gh_call1(material.subclass.
			     material_function_data->material_func,
			     ctl_convert_vector3_to_scm(p));
	       material_type_input(mo, &m);
	       if (destroy_material)
		    material_type_destroy(material);
	       material = m;
	       destroy_material = 1;
	  }
	       
	  switch (material.which_subclass) {
	      case DIELECTRIC:
	      {
		   real eps_val = material.subclass.dielectric_data->epsilon;
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
		   break;
	      }
	      case DIELECTRIC_ANISOTROPIC:
	      {
		   dielectric_anisotropic *d =
			material.subclass.dielectric_anisotropic_data;
		   eps->m00 = d->epsilon_diag.x;
		   eps->m11 = d->epsilon_diag.y;
		   eps->m22 = d->epsilon_diag.z;
#ifdef WITH_HERMITIAN_EPSILON
		   eps->m01.re = d->epsilon_offdiag.x;
		   eps->m02.re = d->epsilon_offdiag.y;
		   eps->m12.re = d->epsilon_offdiag.z;
		   eps->m01.im = d->epsilon_offdiag_imag.x;
		   eps->m02.im = d->epsilon_offdiag_imag.y;
		   eps->m12.im = d->epsilon_offdiag_imag.z;
#else
		   eps->m01 = d->epsilon_offdiag.x;
		   eps->m02 = d->epsilon_offdiag.y;
		   eps->m12 = d->epsilon_offdiag.z;
		   CHECK(d->epsilon_offdiag_imag.x == 0.0 &&
			 d->epsilon_offdiag_imag.y == 0.0 &&
			 d->epsilon_offdiag_imag.z == 0.0,
			 "epsilon-offdiag-imag is only supported when MPB is configured --with-hermitian-epsilon");
#endif
		   maxwell_sym_matrix_invert(eps_inv, eps);
		   break;
	      }
	      case MATERIAL_FUNCTION:
		   CHECK(0, "invalid use of material-function");
		   break;
	      case MATERIAL_TYPE_SELF:
		   CHECK(0, "invalid use of material-type");
		   break;
	  }
	  if (destroy_material)
	       material_type_destroy(material);
     }
}

/**************************************************************************/

/* initialize the field to random numbers; should only be called
   after init-params.  (Guile-callable.) */
void randomize_fields(void)
{
     int i;

     if (!mdata) {
	  mpi_one_fprintf(stderr,
		  "init-params must be called before randomize-fields!\n");
	  return;
     }
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
     int mesh[3];
     int have_old_fields = 0;
     int tree_depth, tree_nobjects;
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

     mpi_one_printf("Mesh size is %d.\n", mesh_size);
     mesh[0] = mesh_size;
     mesh[1] = (dimensions > 1) ? mesh_size : 1;
     mesh[2] = (dimensions > 2) ? mesh_size : 1;

     Rm.c0 = vector3_scale(nx == 1 ? 1 : geometry_lattice.size.x, 
			   geometry_lattice.basis.c0);
     Rm.c1 = vector3_scale(ny == 1 ? 1 : geometry_lattice.size.y,
			   geometry_lattice.basis.c1);
     Rm.c2 = vector3_scale(nz == 1 ? 1 : geometry_lattice.size.z,
			   geometry_lattice.basis.c2);
     mpi_one_printf("Lattice vectors:\n");
     mpi_one_printf("     (%g, %g, %g)\n", Rm.c0.x, Rm.c0.y, Rm.c0.z);  
     mpi_one_printf("     (%g, %g, %g)\n", Rm.c1.x, Rm.c1.y, Rm.c1.z);
     mpi_one_printf("     (%g, %g, %g)\n", Rm.c2.x, Rm.c2.y, Rm.c2.z);
  
     Gm = matrix3x3_inverse(matrix3x3_transpose(Rm));
     mpi_one_printf("Reciprocal lattice vectors (/ 2 pi):\n");
     mpi_one_printf("     (%g, %g, %g)\n", Gm.c0.x, Gm.c0.y, Gm.c0.z);  
     mpi_one_printf("     (%g, %g, %g)\n", Gm.c1.x, Gm.c1.y, Gm.c1.z);
     mpi_one_printf("     (%g, %g, %g)\n", Gm.c2.x, Gm.c2.y, Gm.c2.z);
     
     if (eigensolver_nwork > MAX_NWORK) {
	  mpi_one_printf("(Reducing nwork = %d to maximum: %d.)\n",
		 eigensolver_nwork, MAX_NWORK);
	  eigensolver_nwork = MAX_NWORK;
     }

     matrix3x3_to_arr(R, Rm);
     matrix3x3_to_arr(G, Gm);

     /* we must do this to correct for a non-orthogonal lattice basis: */
     geom_fix_objects();

     mpi_one_printf("Geometric objects:\n");
     if (mpi_is_master())
	  for (i = 0; i < geometry.num_items; ++i) {
	       display_geometric_object_info(5, geometry.items[i]);
	       
	       if (geometry.items[i].material.which_subclass == DIELECTRIC)
		    printf("%*sdielectric constant epsilon = %g\n",
			   5 + 5, "",
			   geometry.items[i].material.
			   subclass.dielectric_data->epsilon);
	  }

     destroy_geom_box_tree(geometry_tree);  /* destroy any tree from
					       previous runs */
     geometry_tree =  create_geom_box_tree();
     if (verbose && mpi_is_master()) {
	  printf("Geometry object bounding box tree:\n");
	  display_geom_box_tree(5, geometry_tree);
     }
     geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
     mpi_one_printf("Geometric object tree has depth %d and %d object nodes"
	    " (vs. %d actual objects)\n",
	    tree_depth, tree_nobjects, geometry.num_items);

     mpi_one_printf("%d k-points:\n", k_points.num_items);
     for (i = 0; i < k_points.num_items; ++i)
	  mpi_one_printf("     (%g,%g,%g)\n", k_points.items[i].x,
			 k_points.items[i].y, k_points.items[i].z);
     
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

     mpi_one_printf("Initializing dielectric function...\n");
     {
	  epsilon_func_data d;
	  get_epsilon_file_func(epsilon_input_file,
				&d.eps_file_func, &d.eps_file_func_data);
	  set_maxwell_dielectric(mdata, mesh, R, G, epsilon_func, &d);
	  destroy_epsilon_file_func_data(d.eps_file_func_data);
     }

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
	  printf("%sfreqs:, k index, kx, ky, kz, kmag/2pi",
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
		    kpoint_index, k[0], k[1], k[2],
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
     number *gv_scratch;
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
     CHK_MALLOC(gv_scratch, number, group_v.num_items);
     
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
	  evectmatrix_XtY_diag_real(Hblock, W[0],
				    group_v.items + ib, gv_scratch);
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

/**************************************************************************/

/* The following routines take the eigenvectors computed by solve-kpoint
   and compute the field (D, H, or E) in position space for one of the bands.
   This field is stored in the global curfield (actually an alias for
   mdata->fft_data, since the latter is unused and big enough).  This
   field can then be manipulated with subsequent "*-field-*" functions
   below.  You can also get the scalar field, epsilon.

   All of these functions are designed to be called by the user
   via Guile. */

void get_dfield(int which_band)
{
     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-dfield!\n");
	  return;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr,
			  "solve-kpoint must be called before get-dfield!\n");
	  return;
     }
     if (which_band < 1 || which_band > H.p) {
	  mpi_one_fprintf(stderr,
			  "must have 1 <= band index <= num_bands (%d)\n",H.p);
	  return;
     }

     curfield = (scalar_complex *) mdata->fft_data;
     curfield_band = which_band;
     curfield_type = 'd';
     maxwell_compute_d_from_H(mdata, H, curfield, which_band - 1, 1);

     /* Here, we correct for the fact that compute_d_from_H actually
	computes just (k+G) x H, whereas the actual D field is
	i/omega i(k+G) x H...so, there is an added factor of -1/omega. */
     {
	  int i, N;
	  double scale;
	  N = mdata->fft_output_size;

	  if (freqs.items[which_band - 1] != 0.0) {
	       scale = -1.0 / freqs.items[which_band - 1];
	  }
	  else
	       scale = -1.0; /* arbitrary */

	  for (i = 0; i < 3*N; ++i) {
	       curfield[i].re *= scale;
	       curfield[i].im *= scale;
	  }
     }
}

void get_hfield(integer which_band)
{
     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-hfield!\n");
	  return;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr,
			  "solve-kpoint must be called before get-hfield!\n");
	  return;
     }
     if (which_band < 1 || which_band > H.p) {
	  mpi_one_fprintf(stderr,
			  "must have 1 <= band index <= num_bands (%d)\n",H.p);
	  return;
     }

     curfield = (scalar_complex *) mdata->fft_data;
     curfield_band = which_band;
     curfield_type = 'h';
     maxwell_compute_h_from_H(mdata, H, curfield, which_band - 1, 1);
}

void get_efield_from_dfield(void)
{
     if (!curfield || curfield_type != 'd') {
	  mpi_one_fprintf(stderr, "get-dfield must be called before "
		  "get-efield-from-dfield!\n");
	  return;
     }
     CHECK(mdata, "unexpected NULL mdata");
     maxwell_compute_e_from_d(mdata, curfield, 1);
     curfield_type = 'e';
}

void get_efield(integer which_band)
{
     get_dfield(which_band);
     get_efield_from_dfield();
}

/* Extract the mean epsilon from the effective inverse dielectric tensor,
   which contains two eigenvalues that correspond to the mean epsilon,
   and one which corresponds to the harmonic mean. */
static real mean_epsilon(const symmetric_matrix *eps_inv)
{
     real eps_eigs[3];
     maxwell_sym_matrix_eigs(eps_eigs, eps_inv);
     /* the harmonic mean should be the largest eigenvalue (smallest
	epsilon), so we'll ignore it and average the other two: */
     return 2.0 / (eps_eigs[0] + eps_eigs[1]);
}

/* get the dielectric function, and compute some statistics */
void get_epsilon(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *epsilon;
     real eps_mean = 0, eps_inv_mean = 0, eps_high = -1e20, eps_low = 1e20;
     int fill_count = 0;

     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-epsilon!\n");
	  return;
     }

     curfield = (scalar_complex *) mdata->fft_data;
     epsilon = (real *) curfield;
     curfield_band = 0;
     curfield_type = 'n';

     /* get epsilon.  Recall that we actually have an inverse
	dielectric tensor at each point; define an average index by
	the inverse of the average eigenvalue of the 1/eps tensor.
	i.e. 3/(trace 1/eps). */

     N = mdata->fft_output_size;
     last_dim = mdata->last_dim;
     last_dim_stored =
	  mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

     for (i = 0; i < N; ++i) {
          epsilon[i] = mean_epsilon(mdata->eps_inv + i);
	  if (epsilon[i] < eps_low)
	       eps_low = epsilon[i];
	  if (epsilon[i] > eps_high)
	       eps_high = epsilon[i];
	  eps_mean += epsilon[i];
	  eps_inv_mean += 1/epsilon[i];
	  if (epsilon[i] > 1.0001)
	       ++fill_count;
#ifndef SCALAR_COMPLEX
	  /* most points need to be counted twice, by rfftw output symmetry: */
	  {
	       int last_index;
#  ifdef HAVE_MPI
	       if (nz == 1) /* 2d calculation: 1st dim. is truncated one */
		    last_index = i / nx + local_y_start;
	       else
		    last_index = i % last_dim_stored;
#  else
	       last_index = i % last_dim_stored;
#  endif
	       if (last_index != 0 && 2*last_index != last_dim) {
		    eps_mean += epsilon[i];
		    eps_inv_mean += 1/epsilon[i];
		    if (epsilon[i] > 1.0001)
			 ++fill_count;
	       }
	  }
#endif
     }

     mpi_allreduce_1(&eps_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce_1(&eps_inv_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce_1(&eps_low, real, SCALAR_MPI_TYPE,
		     MPI_MIN, MPI_COMM_WORLD);
     mpi_allreduce_1(&eps_high, real, SCALAR_MPI_TYPE,
		     MPI_MAX, MPI_COMM_WORLD);
     mpi_allreduce_1(&fill_count, int, MPI_INT,
                   MPI_SUM, MPI_COMM_WORLD);
     N = mdata->nx * mdata->ny * mdata->nz;
     eps_mean /= N;
     eps_inv_mean = N/eps_inv_mean;

     mpi_one_printf("epsilon: %g-%g, mean %g, harm. mean %g, "
		    "%g%% > 1, %g%% \"fill\"\n",
		    eps_low, eps_high, eps_mean, eps_inv_mean,
		    (100.0 * fill_count) / N, 
		    100.0 * (eps_mean-eps_low) / (eps_high-eps_low));
}

/* get the specified component of the dielectric tensor,
   or the inverse tensor if inv != 0 */
void get_epsilon_tensor(int c1, int c2, int imag, int inv)
{
     int i, N;
     real *epsilon;
     int conj = 0, offset = 0;

     curfield_type = '-'; /* only used internally, for now */
     epsilon = (real *) mdata->fft_data;
     N = mdata->fft_output_size;

     switch (c1 * 3 + c2) {
	 case 0:
	      offset = offsetof(symmetric_matrix, m00);
	      break;
	 case 1:
	      offset = offsetof(symmetric_matrix, m01);
	      break;
	 case 2:
	      offset = offsetof(symmetric_matrix, m02);
	      break;
	 case 3:
	      offset = offsetof(symmetric_matrix, m01); /* = conj(m10) */
	      conj = imag;
	      break;
	 case 4:
	      offset = offsetof(symmetric_matrix, m11);
	      break;
	 case 5:
	      offset = offsetof(symmetric_matrix, m12);
	      break;
	 case 6:
	      offset = offsetof(symmetric_matrix, m02); /* = conj(m20) */
	      conj = imag;
	      break;
	 case 7:
	      offset = offsetof(symmetric_matrix, m12); /* = conj(m21) */
	      conj = imag;
	      break;
	 case 8:
	      offset = offsetof(symmetric_matrix, m22);
	      break;
     }

#ifdef WITH_HERMITIAN_EPSILON
     if (c1 != c2 && imag)
	  offset += offsetof(scalar_complex, im);
#endif

     for (i = 0; i < N; ++i) {
	  if (inv) {
	       epsilon[i] = 
		    *((real *) (((char *) &mdata->eps_inv[i]) + offset));
	  }
	  else {
	       symmetric_matrix eps;
	       maxwell_sym_matrix_invert(&eps, &mdata->eps_inv[i]);
	       epsilon[i] = *((real *) (((char *) &eps) + offset));
	  }
	  if (conj)
	       epsilon[i] = -epsilon[i];
     }     
}

/**************************************************************************/

/* Replace curfield (either d or h) with the scalar energy density function,
   normalized to one.  While we're at it, compute some statistics about
   the relative strength of different field components.  Also return
   the integral of the energy density, which we used to normalize it,
   in case the user needs the unnormalized version. */
number_list compute_field_energy(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real comp_sum2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, comp_sum[6];
     real energy_sum = 0.0, normalization;
     real *energy_density = (real *) curfield;
     number_list retval = { 0, 0 };

     if (!curfield || !strchr("dh", curfield_type)) {
	  mpi_one_fprintf(stderr, "The D or H field must be loaded first.\n");
	  return retval;
     }

     N = mdata->fft_output_size;
     last_dim = mdata->last_dim;
     last_dim_stored =
	  mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

     for (i = 0; i < N; ++i) {
	  scalar_complex field[3];
	  real
	       comp_sqr0,comp_sqr1,comp_sqr2,comp_sqr3,comp_sqr4,comp_sqr5;

	  /* energy is either |curfield|^2 or |curfield|^2 / epsilon,
	     depending upon whether it is H or D. */
	  if (curfield_type == 'h') {
	       field[0] =   curfield[3*i];
	       field[1] = curfield[3*i+1];
	       field[2] = curfield[3*i+2];
	  }
	  else
	       assign_symmatrix_vector(field, mdata->eps_inv[i], curfield+3*i);

	  comp_sum2[0] += comp_sqr0 = field[0].re *   curfield[3*i].re;
	  comp_sum2[1] += comp_sqr1 = field[0].im *   curfield[3*i].im;
	  comp_sum2[2] += comp_sqr2 = field[1].re * curfield[3*i+1].re;
	  comp_sum2[3] += comp_sqr3 = field[1].im * curfield[3*i+1].im;
	  comp_sum2[4] += comp_sqr4 = field[2].re * curfield[3*i+2].re;
	  comp_sum2[5] += comp_sqr5 = field[2].im * curfield[3*i+2].im;

	  /* Note: here, we write to energy_density[i]; this is
	     safe, even though energy_density is aliased to curfield,
	     since energy_density[i] is guaranteed to come at or before
	     curfield[i] (which we are now done with). */

	  energy_sum += energy_density[i] = 
	       comp_sqr0+comp_sqr1+comp_sqr2+comp_sqr3+comp_sqr4+comp_sqr5;
#ifndef SCALAR_COMPLEX
	  /* most points need to be counted twice, by rfftw output symmetry: */
	  {
	       int last_index;
#  ifdef HAVE_MPI
	       if (nz == 1) /* 2d calculation: 1st dim. is truncated one */
		    last_index = i / nx + local_y_start;
	       else
		    last_index = i % last_dim_stored;
#  else
	       last_index = i % last_dim_stored;
#  endif
	       if (last_index != 0 && 2*last_index != last_dim) {
		    energy_sum += energy_density[i];
		    comp_sum2[0] += comp_sqr0;
		    comp_sum2[1] += comp_sqr1;
		    comp_sum2[2] += comp_sqr2;
		    comp_sum2[3] += comp_sqr3;
		    comp_sum2[4] += comp_sqr4;
		    comp_sum2[5] += comp_sqr5;
	       }
	  }
#endif
     }

     mpi_allreduce_1(&energy_sum, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce(comp_sum2, comp_sum, 6, real, SCALAR_MPI_TYPE,
                   MPI_SUM, MPI_COMM_WORLD);

     normalization = 1.0 / (energy_sum == 0 ? 1 : energy_sum);
     for (i = 0; i < N; ++i)
	  energy_density[i] *= normalization;

     mpi_one_printf("%c-energy-components:, %d, %d",
	    curfield_type, kpoint_index, curfield_band);
     for (i = 0; i < 6; ++i) {
	  comp_sum[i] *= normalization;
	  if (i % 2 == 1)
	       mpi_one_printf(", %g", comp_sum[i] + comp_sum[i-1]);
     }
     mpi_one_printf("\n");

     /* remember that we now have energy density; denoted by capital D/H */
     curfield_type = toupper(curfield_type);

     /* The return value is a list of 7 items: the total energy,
	followed by the 6 elements of the comp_sum array (the fraction
	of the energy in the real/imag. parts of each field component). */

     retval.num_items = 7;
     CHK_MALLOC(retval.items, number, retval.num_items);

     /* Return the total energy divided by nx*ny*nz to account for the
	scaling of the Fourier transform, so that the integral is
	consistent with the integral in frequency-domain. */
     retval.items[0] = energy_sum/H.N;

     for (i = 0; i < 6; ++i)
	  retval.items[i+1] = comp_sum[i];

     return retval;
}

/**************************************************************************/

/* Fix the phase of the current field (e/h/d) to a canonical value.
   Also changes the phase of the corresponding eigenvector by the
   same amount, so that future calculations will have a consistent
   phase.

   The following procedure is used, derived from a suggestion by Doug
   Allan of Corning: First, choose the phase to maximize the sum of
   the squares of the real parts of the components.  This doesn't fix
   the overall sign, though.  That is done (after incorporating the
   above phase) by: (1) find the largest absolute value of the real
   part, (2) find the point with the greatest spatial array index that
   has |real part| at least half of the largest value, and (3) make
   that point positive.

   In the case of inversion symmetry, on the other hand, the overall phase
   is already fixed, to within a sign, by the choice to make the Fourier
   transform purely real.  So, in that case we simply pick a sign, in
   a manner similar to (2) and (3) above. */
void fix_field_phase(void)
{
     int i, N;
     real sq_sum2[2] = {0,0}, sq_sum[2], maxabs = 0.0;
     int maxabs_index = 0, maxabs_sign = 1;
     double theta;
     scalar phase;

     if (!curfield || !strchr("dhe", curfield_type)) {
          mpi_one_fprintf(stderr, "The D/H/E field must be loaded first.\n");
          return;
     }
     N = mdata->fft_output_size * 3;

#ifdef SCALAR_COMPLEX
     /* Compute the phase that maximizes the sum of the squares of
	the real parts of the components.  Equivalently, maximize
	the real part of the sum of the squares. */
     for (i = 0; i < N; ++i) {
	  real a,b;
	  a = curfield[i].re; b = curfield[i].im;
	  sq_sum2[0] += a*a - b*b;
	  sq_sum2[1] += 2*a*b;
     }
     mpi_allreduce(sq_sum2, sq_sum, 2, real, SCALAR_MPI_TYPE,
                   MPI_SUM, MPI_COMM_WORLD);
     /* compute the phase = exp(i*theta) maximizing the real part of
	the sum of the squares.  i.e., maximize:
	    cos(2*theta)*sq_sum[0] - sin(2*theta)*sq_sum[1] */
     theta = 0.5 * atan2(-sq_sum[1], sq_sum[0]);
     phase.re = cos(theta);
     phase.im = sin(theta);
#else /* ! SCALAR_COMPLEX */
     phase = 1;
#endif /* ! SCALAR_COMPLEX */

     /* Next, fix the overall sign.  We do this by first computing the
	maximum |real part| of the jmax component (after multiplying
	by phase), and then finding the last spatial index at which
	|real part| is at least half of this value.  The sign is then
	chosen to make the real part positive at that point. 

        (Note that we can't just make the point of maximum |real part|
         positive, as that would be ambiguous in the common case of an
         oscillating field within the unit cell.)

        In the case of inversion symmetry (!SCALAR_COMPLEX), we work with
        (real part - imag part) instead of (real part), to insure that we
        have something that is nonzero somewhere. */

     for (i = 0; i < N; ++i) {
#ifdef SCALAR_COMPLEX
	  real r = fabs(curfield[i].re * phase.re - curfield[i].im * phase.im);
#else
	  real r = fabs(curfield[i].re - curfield[i].im);
#endif
	  if (r > maxabs)
	       maxabs = r;
     }
     mpi_allreduce_1(&maxabs, real, SCALAR_MPI_TYPE,
		     MPI_MAX, MPI_COMM_WORLD);
     for (i = N - 1; i >= 0; --i) {
#ifdef SCALAR_COMPLEX
	  real r = curfield[i].re * phase.re - curfield[i].im * phase.im;
#else
	  real r = curfield[i].re - curfield[i].im;
#endif
	  if (fabs(r) >= 0.5 * maxabs) {
	       maxabs_index = i;
	       maxabs_sign = r < 0 ? -1 : 1;
	       break;
	  }
     }
     if (i >= 0)  /* convert index to global index in distributed array: */
	  maxabs_index += mdata->local_y_start * mdata->nx * mdata->nz;
     {
	  /* compute maximum index and corresponding sign over all the 
	     processors, using the MPI_MAXLOC reduction operation: */
	  struct twoint_struct {int i; int s;} x;
	  x.i = maxabs_index; x.s = maxabs_sign;
	  mpi_allreduce_1(&x, struct twoint_struct, MPI_2INT,
			  MPI_MAXLOC, MPI_COMM_WORLD);
	  maxabs_index = x.i; maxabs_sign = x.s;
     }
     ASSIGN_SCALAR(phase,
		   SCALAR_RE(phase)*maxabs_sign, SCALAR_IM(phase)*maxabs_sign);

     mpi_one_printf("Fixing %c-field (band %d) phase by %g + %gi; "
		    "max ampl. = %g\n", curfield_type, curfield_band,
		    SCALAR_RE(phase), SCALAR_IM(phase), maxabs);

     /* Now, multiply everything by this phase, *including* the
	stored "raw" eigenvector in H, so that any future fields
	that we compute will have a consistent phase: */
     for (i = 0; i < N; ++i) {
	  real a,b;
	  a = curfield[i].re; b = curfield[i].im;
	  curfield[i].re = a*SCALAR_RE(phase) - b*SCALAR_IM(phase);
	  curfield[i].im = a*SCALAR_IM(phase) + b*SCALAR_RE(phase);
     }
     for (i = 0; i < H.n; ++i) {
          ASSIGN_MULT(H.data[i*H.p + curfield_band - 1], 
		      H.data[i*H.p + curfield_band - 1], phase);
     }
}

/**************************************************************************/

/* Functions to return epsilon, fields, energies, etcetera, at a specified
   point, linearly interpolating if necessary. */

static real get_val(int ix, int iy, int iz,
		    real *data, int stride, int conjugate)
{
     int nx = mdata->nx, ny = mdata->ny, nz = mdata->nz;
     
#ifndef SCALAR_COMPLEX
     {
	  int nlast = mdata->last_dim_size / 2;
	  if ((nz > 1 ? iz : (ny > 1 ? iy : ix)) >= nlast) {
	       ix = ix ? nx - ix : ix;
	       iy = iy ? ny - iy : iy;
	       iz = iz ? nz - iz : iz;
	       conjugate = conjugate ? 1 : 0;
	  }
	  if (nz > 1) nz = nlast; else if (ny > 1) ny = nlast; else nx = nlast;
     }
#else
     conjugate = 0;
#endif

#ifdef HAVE_MPI
     CHECK(0, "not yet implemented!");
#else
     if (conjugate)
	  return -data[(((ix * ny) + iy) * nz + iz) * stride];
     else
	  return data[(((ix * ny) + iy) * nz + iz) * stride];
#endif
}

static real interp_val(vector3 p, real *data, int stride, int conjugate)
{
     double ipart;
     real rx, ry, rz, dx, dy, dz;
     int x, y, z, x2, y2, z2;
     int nx = mdata->nx, ny = mdata->ny, nz = mdata->nz;

     rx = modf(p.x/geometry_lattice.size.x + 0.5, &ipart); if (rx < 0) rx += 1;
     ry = modf(p.y/geometry_lattice.size.y + 0.5, &ipart); if (ry < 0) ry += 1;
     rz = modf(p.z/geometry_lattice.size.z + 0.5, &ipart); if (rz < 0) rz += 1;

     /* get the point corresponding to r in the grid: */
     x = rx * nx;
     y = ry * ny;
     z = rz * nz;

     /* get the difference between (x,y,z) and the actual point */
     dx = rx * nx - x;
     dy = ry * ny - y;
     dz = rz * nz - z;

     /* get the other closest point in the grid, with periodic boundaries: */
     x2 = (nx + (dx >= 0.0 ? x + 1 : x - 1)) % nx;
     y2 = (ny + (dy >= 0.0 ? y + 1 : y - 1)) % ny;
     z2 = (nz + (dz >= 0.0 ? z + 1 : z - 1)) % nz;

     /* take abs(d{xyz}) to get weights for {xyz} and {xyz}2: */
     dx = fabs(dx);
     dy = fabs(dy);
     dz = fabs(dz);

#define D(x,y,z) (get_val(x,y,z, data, stride, conjugate))

     return(((D(x,y,z)*(1.0-dx) + D(x2,y,z)*dx) * (1.0-dy) +
             (D(x,y2,z)*(1.0-dx) + D(x2,y2,z)*dx) * dy) * (1.0-dz) +
            ((D(x,y,z2)*(1.0-dx) + D(x2,y,z2)*dx) * (1.0-dy) +
             (D(x,y2,z2)*(1.0-dx) + D(x2,y2,z2)*dx) * dy) * dz);

#undef D
}

static scalar_complex interp_cval(vector3 p, real *data, int stride,
				  int conjugate)
{
     scalar_complex cval;
     cval.re = interp_val(p, data, stride, conjugate);
     cval.im = interp_val(p, data + 1, stride, !conjugate);
     return cval;
}

static symmetric_matrix interp_eps_inv(vector3 p)
{
     int stride = sizeof(symmetric_matrix) / sizeof(real);
     symmetric_matrix eps_inv;

     eps_inv.m00 = interp_val(p, &mdata->eps_inv->m00, stride, 0);
     eps_inv.m11 = interp_val(p, &mdata->eps_inv->m11, stride, 0);
     eps_inv.m22 = interp_val(p, &mdata->eps_inv->m22, stride, 0);
#ifdef WITH_HERMITIAN_EPSILON
     eps_inv.m01 = interp_cval(p, &mdata->eps_inv->m01.re, stride, 0);
     eps_inv.m02 = interp_cval(p, &mdata->eps_inv->m02.re, stride, 0);
     eps_inv.m12 = interp_cval(p, &mdata->eps_inv->m12.re, stride, 0);
#else
     eps_inv.m01 = interp_val(p, &mdata->eps_inv->m01, stride, 0);
     eps_inv.m02 = interp_val(p, &mdata->eps_inv->m02, stride, 0);
     eps_inv.m12 = interp_val(p, &mdata->eps_inv->m12, stride, 0);
#endif
     return eps_inv;
}

number get_epsilon_point(vector3 p)
{
     symmetric_matrix eps_inv;
     eps_inv = interp_eps_inv(p);
     return mean_epsilon(&eps_inv);
}

matrix3x3 get_epsilon_inverse_tensor_point(vector3 p)
{
     matrix3x3 eps_inv3x3;
     symmetric_matrix eps_inv;
     eps_inv = interp_eps_inv(p);

     eps_inv3x3.c0.x = eps_inv.m00;
     eps_inv3x3.c1.y = eps_inv.m11;
     eps_inv3x3.c2.z = eps_inv.m22;
#ifdef WITH_HERMITIAN_EPSILON /* only return real part; see below for im */
     eps_inv3x3.c0.y = eps_inv3x3.c1.x = CSCALAR_RE(eps_inv.m01);
     eps_inv3x3.c0.z = eps_inv3x3.c2.x = CSCALAR_RE(eps_inv.m02);
     eps_inv3x3.c1.z = eps_inv3x3.c2.y = CSCALAR_RE(eps_inv.m12);
#else
     eps_inv3x3.c0.y = eps_inv3x3.c1.x = eps_inv.m01;
     eps_inv3x3.c0.z = eps_inv3x3.c2.x = eps_inv.m02;
     eps_inv3x3.c1.z = eps_inv3x3.c2.y = eps_inv.m12;
#endif
     return eps_inv3x3;
}

matrix3x3 get_epsilon_inverse_tensor_im_point(vector3 p)
{
     matrix3x3 eps_inv3x3 = {{0,0,0},{0,0,0},{0,0,0}};
     symmetric_matrix eps_inv;
     eps_inv = interp_eps_inv(p);

#ifdef WITH_HERMITIAN_EPSILON
     eps_inv3x3.c0.y = -(eps_inv3x3.c1.x = CSCALAR_IM(eps_inv.m01));
     eps_inv3x3.c0.z = -(eps_inv3x3.c2.x = CSCALAR_IM(eps_inv.m02));
     eps_inv3x3.c1.z = -(eps_inv3x3.c2.y = CSCALAR_IM(eps_inv.m12));
#endif
     return eps_inv3x3;
}

number get_energy_point(vector3 p)
{
     CHECK(curfield && strchr("DH", curfield_type),
	   "compute-field-energy must be called before get-energy-point");
     return interp_val(p, (real *) curfield, 1, 0);
}

vector3 get_bloch_field_re_point(vector3 p)
{
     scalar_complex field[3];
     int conjugate;
     vector3 F;

     CHECK(curfield && strchr("dhe", curfield_type),
	   "field must be must be loaded before get-*field*-point");
     conjugate = curfield_type != 'h'; /* non pseudo-vectors flip */
     field[0] = interp_cval(p, &curfield[0].re, 6, conjugate);
     field[1] = interp_cval(p, &curfield[1].re, 6, conjugate);
     field[2] = interp_cval(p, &curfield[2].re, 6, conjugate);
     F.x = field[0].re; F.y = field[1].re; F.z = field[2].re;
     return F;
}

vector3 get_bloch_field_im_point(vector3 p)
{
     scalar_complex field[3];
     int conjugate;
     vector3 F;

     CHECK(curfield && strchr("dhe", curfield_type),
	   "field must be must be loaded before get-*field*-point");
     conjugate = curfield_type != 'h'; /* non pseudo-vectors flip */
     field[0] = interp_cval(p, &curfield[0].re, 6, conjugate);
     field[1] = interp_cval(p, &curfield[1].re, 6, conjugate);
     field[2] = interp_cval(p, &curfield[2].re, 6, conjugate);
     F.x = field[0].im; F.y = field[1].im; F.z = field[2].im;
     return F;
}

static void get_field_point_reim(vector3 p, vector3 *re, vector3 *im)
{
     scalar_complex field[3], phase;
     double phase_phi;
     int conjugate;

     CHECK(curfield && strchr("dhe", curfield_type),
	   "field must be must be loaded before get-*field*-point");
     conjugate = curfield_type != 'h'; /* non pseudo-vectors flip */
     field[0] = interp_cval(p, &curfield[0].re, 6, conjugate);
     field[1] = interp_cval(p, &curfield[1].re, 6, conjugate);
     field[2] = interp_cval(p, &curfield[2].re, 6, conjugate);

     phase_phi = TWOPI * (cur_kvector.x * (p.x/geometry_lattice.size.x+0.5) +
			  cur_kvector.y * (p.y/geometry_lattice.size.y+0.5) +
			  cur_kvector.z * (p.z/geometry_lattice.size.z+0.5));
     CASSIGN_SCALAR(phase, cos(phase_phi), sin(phase_phi));
     CASSIGN_MULT(field[0], field[0], phase);
     CASSIGN_MULT(field[1], field[1], phase);
     CASSIGN_MULT(field[2], field[2], phase);

     re->x = field[0].re; re->y = field[1].re; re->z = field[2].re;
     im->x = field[0].im; im->y = field[1].im; im->z = field[2].im;
}

vector3 get_field_re_point(vector3 p)
{
     vector3 Fre, Fim;
     get_field_point_reim(p, &Fre, &Fim);
     return Fre;
}

vector3 get_field_im_point(vector3 p)
{
     vector3 Fre, Fim;
     get_field_point_reim(p, &Fre, &Fim);
     return Fim;
}

/**************************************************************************/

/* compute the fraction of the field energy that is located in the
   given range of dielectric constants: */
number compute_energy_in_dielectric(number eps_low, number eps_high)
{
     int N, i, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *energy = (real *) curfield;
     real epsilon, energy_sum = 0.0;

     if (!curfield || !strchr("DH", curfield_type)) {
          mpi_one_fprintf(stderr, "The D or H energy density must be loaded first.\n");
          return 0.0;
     }

     N = mdata->fft_output_size;
     last_dim = mdata->last_dim;
     last_dim_stored =
	  mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

     for (i = 0; i < N; ++i) {
	  epsilon = mean_epsilon(mdata->eps_inv +i);
	  if (epsilon >= eps_low && epsilon <= eps_high) {
	       energy_sum += energy[i];
#ifndef SCALAR_COMPLEX
	       /* most points are counted twice, by rfftw output symmetry: */
	       {
		    int last_index;
#  ifdef HAVE_MPI
		    if (nz == 1) /* 2d: 1st dim. is truncated one */
			 last_index = i / nx + local_y_start;
		    else
			 last_index = i % last_dim_stored;
#  else
		    last_index = i % last_dim_stored;
#  endif
		    if (last_index != 0 && 2*last_index != last_dim)
			 energy_sum += energy[i];
	       }
#endif
	  }
     }
     mpi_allreduce_1(&energy_sum, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     return energy_sum;
}

/**************************************************************************/

/* Prepend the prefix to the fname, and (if parity_suffix is true)
   append a parity specifier (if any) (e.g. ".te"), returning a new
   string, which should be deallocated with free().  fname or prefix
   may be NULL, in which case they are treated as the empty string. */
static char *fix_fname(const char *fname, const char *prefix,
		       maxwell_data *d, int parity_suffix)
{
     char *s;
     CHK_MALLOC(s, char,
		(fname ? strlen(fname) : 0) + 
		(prefix ? strlen(prefix) : 0) + 20);
     strcpy(s, prefix ? prefix : "");
     strcat(s, fname ? fname : "");
     if (parity_suffix && d->parity != NO_PARITY) {
	  /* assumes parity suffix is less than 20 characters;
	     currently it is less than 12 */
	  strcat(s, ".");
	  strcat(s, parity_string(d));
     }
     return s;
}

static void output_scalarfield(real *vals,
			       const int dims[3],
			       const int local_dims[3],
			       const int start[3],
			       matrixio_id file_id,
			       const char *dataname,
			       int last_dim_index,
			       int last_dim_start, int last_dim_size,
			       int first_dim_start, int first_dim_size,
			       int write_start0_special)
{
     matrixio_id data_id = -1;

     fieldio_write_real_vals(vals, 3, dims, 
			     local_dims, start, file_id, 0, 
			     dataname, &data_id);
     
#ifndef SCALAR_COMPLEX
     {
	  int start_new[3], local_dims_new[3];

	  start_new[0] = start[0];
	  start_new[1] = start[1];
	  start_new[2] = start[2];
	  local_dims_new[0] = local_dims[0];
	  local_dims_new[1] = local_dims[1];
	  local_dims_new[2] = local_dims[2];

	  maxwell_scalarfield_otherhalf(mdata, vals);
	  start_new[last_dim_index] = last_dim_start;
	  local_dims_new[last_dim_index] = last_dim_size;
	  start_new[0] = first_dim_start;
	  local_dims_new[0] = first_dim_size;
	  if (write_start0_special) {
	       /* The conjugated array half may be discontiguous.
		  First, write the part not containing start_new[0], and
		  then write the start_new[0] slab. */
	       fieldio_write_real_vals(vals +
				       local_dims_new[1] * local_dims_new[2],
				       3, dims, local_dims_new, start_new, 
				       file_id, 1, dataname, &data_id);
	       local_dims_new[0] = 1;
	       start_new[0] = 0;
	       fieldio_write_real_vals(vals, 3, dims, 
				       local_dims_new, start_new,
				       file_id, 1, dataname, &data_id);
	  }
	  else {
	       fieldio_write_real_vals(vals, 3, dims, 
				       local_dims_new, start_new, 
				       file_id, 1, dataname, &data_id);
	  }
     }
#endif

     if (data_id >= 0)
	  matrixio_close_dataset(data_id);
}

/* given the field in curfield, store it to HDF (or whatever) using
   the matrixio (fieldio) routines.  Allow the component to be specified
   (which_component 0/1/2 = x/y/z, -1 = all) for vector fields. 
   Also allow the user to specify a prefix string for the filename. */
void output_field_to_file(integer which_component, string filename_prefix)
{
     char fname[100], *fname2, description[100];
     int dims[3], local_dims[3], start[3] = {0,0,0};
     matrixio_id file_id = -1;
     int attr_dims[2] = {3, 3};
     real output_k[3]; /* kvector in reciprocal lattice basis */
     real output_R[3][3];

     /* where to put "otherhalf" block of output, only used for real scalars */
     int last_dim_index = 0;
     int last_dim_start = 0, last_dim_size = 0;
     int first_dim_start = 0, first_dim_size = 0;
     int write_start0_special = 0;

     if (!curfield) {
	  mpi_one_fprintf(stderr, 
		  "fields, energy dens., or epsilon must be loaded first.\n");
	  return;
     }
     
#ifdef HAVE_MPI
     /* The first two dimensions (x and y) of the position-space fields
	are transposed when we use MPI, so we need to transpose everything. */
     dims[0] = mdata->ny;
     local_dims[1] = dims[1] = mdata->nx;
     local_dims[2] = dims[2] = mdata->nz;
     local_dims[0] = mdata->local_ny;
     start[0] = mdata->local_y_start;
#  ifndef SCALAR_COMPLEX
     /* Ugh, hairy.  See also maxwell_vectorfield_otherhalf. */
     if (dims[2] == 1) {
	  last_dim_index = 0;
	  first_dim_size = local_dims[0];
	  first_dim_start = dims[0] - (start[0] + local_dims[0] - 1);

	  if (start[0] == 0)
	       --first_dim_size; /* DC frequency is not in other half */
	  if (start[0] + local_dims[0] == mdata->last_dim_size / 2 &&
	      dims[0] % 2 == 0) {
	       --first_dim_size; /* Nyquist frequency is not in other half */
	       ++first_dim_start;
	  }

	  last_dim_start = first_dim_start;
	  last_dim_size = first_dim_size;
     }
     else {
	  last_dim_index = 2;
	  local_dims[last_dim_index] = mdata->last_dim_size / 2;
	  if (start[0] == 0) {
	       first_dim_size = local_dims[0] - 1;
	       first_dim_start = dims[0] - first_dim_size;
	       write_start0_special = 1;
	  }
	  else {
	       first_dim_start = dims[0] - (start[0] + local_dims[0] - 1);
	       first_dim_size = local_dims[0];
	  }
	  last_dim_start = local_dims[last_dim_index];
	  last_dim_size = dims[last_dim_index] - local_dims[last_dim_index];
     }
#  endif /* ! SCALAR_COMPLEX */
     output_k[0] = R[1][0]*mdata->current_k[0] + R[1][1]*mdata->current_k[1]
	  + R[1][2]*mdata->current_k[2];
     output_k[1] = R[0][0]*mdata->current_k[0] + R[0][1]*mdata->current_k[1]
	  + R[0][2]*mdata->current_k[2];
     output_k[2] = R[2][0]*mdata->current_k[0] + R[2][1]*mdata->current_k[1]
	  + R[2][2]*mdata->current_k[2];
     output_R[0][0]=R[1][0]; output_R[0][1]=R[1][1]; output_R[0][2]=R[1][2];
     output_R[1][0]=R[0][0]; output_R[1][1]=R[0][1]; output_R[1][2]=R[0][2];
     output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
#else /* ! HAVE_MPI */
     dims[0] = mdata->nx;
     local_dims[1] = dims[1] = mdata->ny;
     local_dims[2] = dims[2] = mdata->nz;
     local_dims[0] = mdata->local_nx;
#  ifndef SCALAR_COMPLEX
     last_dim_index = dims[2] == 1 ? (dims[1] == 1 ? 0 : 1) : 2;
     local_dims[last_dim_index] = mdata->last_dim_size / 2;
     last_dim_start = local_dims[last_dim_index];
     last_dim_size = dims[last_dim_index] - local_dims[last_dim_index];
     first_dim_start = last_dim_index ? 0 : last_dim_start;
     first_dim_size = last_dim_index ? local_dims[0] : last_dim_size;
#  endif
     start[0] = mdata->local_x_start;
     output_k[0] = R[0][0]*mdata->current_k[0] + R[0][1]*mdata->current_k[1]
	  + R[0][2]*mdata->current_k[2];
     output_k[1] = R[1][0]*mdata->current_k[0] + R[1][1]*mdata->current_k[1]
	  + R[1][2]*mdata->current_k[2];
     output_k[2] = R[2][0]*mdata->current_k[0] + R[2][1]*mdata->current_k[1]
	  + R[2][2]*mdata->current_k[2];
     output_R[0][0]=R[0][0]; output_R[0][1]=R[0][1]; output_R[0][2]=R[0][2];
     output_R[1][0]=R[1][0]; output_R[1][1]=R[1][1]; output_R[1][2]=R[1][2];
     output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
#endif /* ! HAVE_MPI */
     
     if (strchr("dhe", curfield_type)) { /* outputting vector field */
	  matrixio_id data_id[6] = {-1,-1,-1,-1,-1,-1};
	  int i;

	  sprintf(fname, "%c.k%02d.b%02d",
		  curfield_type, kpoint_index, curfield_band);
	  if (which_component >= 0) {
	       char comp_str[] = ".x";
	       comp_str[1] = 'x' + which_component;
	       strcat(fname, comp_str);
	  }
	  sprintf(description, "%c field, kpoint %d, band %d, freq=%g",
		  curfield_type, kpoint_index, curfield_band, 
		  freqs.items[curfield_band - 1]);
	  fname2 = fix_fname(fname, filename_prefix, mdata, 1);
	  mpi_one_printf("Outputting fields to %s...\n", fname2);
	  file_id = matrixio_create(fname2);
	  free(fname2);
	  fieldio_write_complex_field(curfield, 3, dims, local_dims, start,
				      which_component, output_k,
				      file_id, 0, data_id);

#ifndef SCALAR_COMPLEX
	  /* Here's where it gets hairy. */
	  maxwell_vectorfield_otherhalf(mdata, curfield,
					output_k[0], output_k[1], output_k[2]);
	  start[last_dim_index] = last_dim_start;
	  local_dims[last_dim_index] = last_dim_size;
	  start[0] = first_dim_start;
	  local_dims[0] = first_dim_size;
	  if (write_start0_special) {
	       /* The conjugated array half may be discontiguous.
		  First, write the part not containing start[0], and
		  then write the start[0] slab. */
	       fieldio_write_complex_field(curfield + 
					   3 * local_dims[1] * local_dims[2],
					   3, dims, local_dims, start,
					   which_component, NULL, 
					   file_id, 1, data_id);
	       local_dims[0] = 1;
	       start[0] = 0;
	       fieldio_write_complex_field(curfield, 3, dims,local_dims,start,
					   which_component, NULL,
					   file_id, 1, data_id);
	  }
	  else {
	       fieldio_write_complex_field(curfield, 3, dims,local_dims,start,
					   which_component, NULL,
					   file_id, 1, data_id);
	  }
#endif

	  for (i = 0; i < 6; ++i)
	       if (data_id[i] >= 0)
		    matrixio_close_dataset(data_id[i]);
	  matrixio_write_data_attr(file_id, "Bloch wavevector",
				   output_k, 1, attr_dims);
     }
     else if (strchr("DHn", curfield_type)) { /* scalar field */
	  if (curfield_type == 'n') {
	       sprintf(fname, "epsilon");
	       sprintf(description, "dielectric function, epsilon");
	  }
	  else {
	       sprintf(fname, "%cpwr.k%02d.b%02d",
		       tolower(curfield_type), kpoint_index, curfield_band);
	       sprintf(description,
		       "%c field energy density, kpoint %d, band %d, freq=%g",
		       curfield_type, kpoint_index, curfield_band, 
		       freqs.items[curfield_band - 1]);
	  }
	  fname2 = fix_fname(fname, filename_prefix, mdata, 
			     /* no parity suffix for epsilon: */
			     curfield_type != 'n');
	  mpi_one_printf("Outputting %s...\n", fname2);
	  file_id = matrixio_create(fname2);
	  free(fname2);

	  output_scalarfield((real *) curfield, dims, 
			     local_dims, start, file_id, "data",
			     last_dim_index, last_dim_start, last_dim_size,
			     first_dim_start, first_dim_size,
			     write_start0_special);

	  if (curfield_type == 'n') {
	       int c1, c2, inv;
	       char dataname[100];

	       for (inv = 0; inv < 2; ++inv)
		    for (c1 = 0; c1 < 3; ++c1)
			 for (c2 = c1; c2 < 3; ++c2) {
			      get_epsilon_tensor(c1,c2, 0, inv);
			      sprintf(dataname, "%s.%c%c",
				      inv ? "epsilon_inverse" : "epsilon",
				      c1 + 'x', c2 + 'x');
			      output_scalarfield((real *) curfield, dims,
						 local_dims, start,
						 file_id, dataname,
						 last_dim_index,
						 last_dim_start, last_dim_size,
						 first_dim_start,
						 first_dim_size,
						 write_start0_special);
#if defined(WITH_HERMITIAN_EPSILON)
			      if (c1 != c2) {
				   get_epsilon_tensor(c1,c2, 1, inv);
				   strcat(dataname, ".i");
#ifndef SCALAR_COMPLEX /* scalarfield_otherhalf isn't right */
				   strcat(dataname, ".screwy");
#endif
				   output_scalarfield((real *) curfield, dims,
						      local_dims, start,
						      file_id, dataname,
						      last_dim_index,
						      last_dim_start,
						      last_dim_size,
						      first_dim_start,
						      first_dim_size,
						      write_start0_special);
			      }
#endif
			 }
	  }

     }
     else
	  mpi_one_fprintf(stderr, "unknown field type!\n");

     if (file_id >= 0) {
	  matrixio_write_data_attr(file_id, "lattice vectors",
				   &output_R[0][0], 2, attr_dims);
	  matrixio_write_string_attr(file_id, "description", description);

	  matrixio_close(file_id);
     }

     /* We have destroyed curfield (by multiplying it by phases,
	and/or reorganizing in the case of real-amplitude fields). */
     curfield_reset();
}

/**************************************************************************/

/* For curfield an energy density, compute the fraction of the energy
   that resides inside the given list of geometric objects.   Later
   objects in the list have precedence, just like the ordinary
   geometry list. */
number compute_energy_in_object_list(geometric_object_list objects)
{
     int i, j, k, n1, n2, n3, n_other, n_last, rank, last_dim;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
     real s1, s2, s3, c1, c2, c3;
     real *energy = (real *) curfield;
     real energy_sum = 0;

     if (!curfield || !strchr("DH", curfield_type)) {
          mpi_one_fprintf(stderr, "The D or H energy density must be loaded first.\n");
          return 0.0;
     }

     for (i = 0; i < objects.num_items; ++i)
	  geom_fix_object(objects.items[i]);

     n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;
     n_other = mdata->other_dims;
     n_last = mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     last_dim = mdata->last_dim;
     rank = (n3 == 1) ? (n2 == 1 ? 1 : 2) : 3;

     s1 = geometry_lattice.size.x / n1;
     s2 = geometry_lattice.size.y / n2;
     s3 = geometry_lattice.size.z / n3;
     c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5;
     c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
     c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

     /* Here we have different loops over the coordinates, depending
	upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and "index" describing the corresponding index in 
	the curfield array.

        This was all stolen from maxwell_eps.c...it would be better
        if we didn't have to cut and paste, sigh. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI
     
     for (i = 0; i < n1; ++i)
	  for (j = 0; j < n2; ++j)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j, k2 = k;
	  int index = ((i * n2 + j) * n3 + k);

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j + local_y_start, k2 = k;
	  int index = ((j * n1 + i) * n3 + k);

#  endif /* HAVE_MPI */

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

     for (i = 0; i < n_other; ++i)
	  for (j = 0; j < n_last; ++j)
     {
	  int index = i * n_last + j;
	  int i2, j2, k2;
	  switch (rank) {
	      case 2: i2 = i; j2 = j; k2 = 0; break;
	      case 3: i2 = i / n2; j2 = i % n2; k2 = j; break;
	      default: i2 = j; j2 = k2 = 0;  break;
	  }

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = mdata->last_dim_size / 2;
     else
	  local_n3 = 1;
     
     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < local_n3; ++k)
     {
#         define i2 i
	  int j2 = j + local_y_start;
#         define k2 k
	  int index = ((j * n1 + i) * local_n3 + k);

#  endif /* HAVE_MPI */

#endif /* not SCALAR_COMPLEX */

	  {
	       vector3 p;
	       int n;
	       p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;
	       for (n = objects.num_items - 1; n >= 0; --n)
		    if (point_in_periodic_fixed_objectp(p, objects.items[n])) {
			 if (objects.items[n].material.which_subclass
			     == MATERIAL_TYPE_SELF)
			      break; /* treat as a "nothing" object */
			 energy_sum += energy[index];
#ifndef SCALAR_COMPLEX
			 {
			      int last_index;
#  ifdef HAVE_MPI
			      if (n3 == 1)
				   last_index = j + local_y_start;
			      else
				   last_index = k;
#  else
			      last_index = j;
#  endif
			      if (last_index != 0 && 2*last_index != last_dim)
				   energy_sum += energy[index];
			 }
#endif
			 break;
		    }
	  }
     }

     mpi_allreduce_1(&energy_sum, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     return energy_sum;
}

/**************************************************************************/

#ifndef HAVE_GH_CALL4
static SCM gh_call4 (SCM proc, SCM arg1, SCM arg2, SCM arg3, SCM arg4)
{
     return scm_apply (proc, arg1, 
		       scm_cons(arg2, 
				scm_cons2 (arg3, arg4, scm_listofnull)));
}
#endif


/* Compute the integral of f(energy/field, epsilon, r) over the cell. */
number compute_field_integral(function f)
{
     int i, j, k, n1, n2, n3, n_other, n_last, rank, last_dim;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
     real s1, s2, s3, c1, c2, c3;
     int integrate_energy, hfield;
     real *energy = (real *) curfield;
     real integral = 0;

     if (!curfield || !strchr("dheDH", curfield_type)) {
          mpi_one_fprintf(stderr, "The D or H energy/field must be loaded first.\n");
          return 0.0;
     }

     integrate_energy = strchr("DH", curfield_type) != NULL;
     hfield = strchr("h", curfield_type) != NULL;

     n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;
     n_other = mdata->other_dims;
     n_last = mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     last_dim = mdata->last_dim;
     rank = (n3 == 1) ? (n2 == 1 ? 1 : 2) : 3;

     s1 = geometry_lattice.size.x / n1;
     s2 = geometry_lattice.size.y / n2;
     s3 = geometry_lattice.size.z / n3;
     c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5;
     c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
     c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

     /* Here we have different loops over the coordinates, depending
	upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and "index" describing the corresponding index in 
	the curfield array.

        This was all stolen from maxwell_eps.c...it would be better
        if we didn't have to cut and paste, sigh. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI
     
     for (i = 0; i < n1; ++i)
	  for (j = 0; j < n2; ++j)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j, k2 = k;
	  int index = ((i * n2 + j) * n3 + k);

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j + local_y_start, k2 = k;
	  int index = ((j * n1 + i) * n3 + k);

#  endif /* HAVE_MPI */

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

     for (i = 0; i < n_other; ++i)
	  for (j = 0; j < n_last; ++j)
     {
	  int index = i * n_last + j;
	  int i2, j2, k2;
	  switch (rank) {
	      case 2: i2 = i; j2 = j; k2 = 0; break;
	      case 3: i2 = i / n2; j2 = i % n2; k2 = j; break;
	      default: i2 = j; j2 = k2 = 0;  break;
	  }

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = mdata->last_dim_size / 2;
     else
	  local_n3 = 1;
     
     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < local_n3; ++k)
     {
#         define i2 i
	  int j2 = j + local_y_start;
#         define k2 k
	  int index = ((j * n1 + i) * local_n3 + k);

#  endif /* HAVE_MPI */

#endif /* not SCALAR_COMPLEX */

	  {
	       real epsilon;
	       vector3 p;

	       epsilon = mean_epsilon(mdata->eps_inv + index);
	       
	       p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;
	       if (integrate_energy)
		    integral += 
			 ctl_convert_number_to_c(
			      gh_call3(f,
				     ctl_convert_number_to_scm(energy[index]),
				       ctl_convert_number_to_scm(epsilon),
				       ctl_convert_vector3_to_scm(p)));	  
	       else {
		    vector3 Fre, Fim;
		    double phase_phi;
		    scalar_complex phase;

		    phase_phi = TWOPI * 
			 (cur_kvector.x * (p.x/geometry_lattice.size.x+0.5) +
			  cur_kvector.y * (p.y/geometry_lattice.size.y+0.5) +
			  cur_kvector.z * (p.z/geometry_lattice.size.z+0.5));
		    CASSIGN_SCALAR(phase, cos(phase_phi), sin(phase_phi));
		    CASSIGN_MULT_RE(Fre.x, curfield[3*index+0], phase);
		    CASSIGN_MULT_IM(Fim.x, curfield[3*index+0], phase);
		    CASSIGN_MULT_RE(Fre.y, curfield[3*index+1], phase);
		    CASSIGN_MULT_IM(Fim.y, curfield[3*index+1], phase);
		    CASSIGN_MULT_RE(Fre.z, curfield[3*index+2], phase);
		    CASSIGN_MULT_IM(Fim.z, curfield[3*index+2], phase);
		    integral +=
                         ctl_convert_number_to_c(
                              gh_call4(f,
				       ctl_convert_vector3_to_scm(Fre),
				       ctl_convert_vector3_to_scm(Fim),
				       ctl_convert_number_to_scm(epsilon),
                                       ctl_convert_vector3_to_scm(p)));
	       }

#ifndef SCALAR_COMPLEX
	       {
		    int last_index;
#  ifdef HAVE_MPI
		    if (n3 == 1)
			 last_index = j + local_y_start;
		    else
			 last_index = k;
#  else
		    last_index = j;
#  endif
		    
		    if (last_index != 0 && 2*last_index != last_dim) {
			 int i2c, j2c, k2c;
			 i2c = i2 ? (n1 - i2) : 0;
			 j2c = j2 ? (n2 - j2) : 0;
			 k2c = k2 ? (n3 - k2) : 0;
			 p.x = i2c * s1 - c1; 
			 p.y = j2c * s2 - c2; 
			 p.z = k2c * s3 - c3;
			 if (integrate_energy)
			      integral += 
				   ctl_convert_number_to_c(
					gh_call3(f,
				      ctl_convert_number_to_scm(energy[index]),
				      ctl_convert_number_to_scm(epsilon),
				      ctl_convert_vector3_to_scm(p)));
			 else {
			      vector3 Fre, Fim;
			      double phase_phi;
			      scalar_complex phase, Fx, Fy, Fz;
			      
			      Fx = curfield[3*index+0];
			      Fy = curfield[3*index+1];
			      Fz = curfield[3*index+2];
			      /* conjugate and flip field */
			      if (hfield) { /* pseudo-vectors aren't flipped */
				   Fx.im= -Fx.im; Fy.im= -Fy.im; Fz.im= -Fz.im;
			      }
			      else {
				   Fx.re= -Fx.re; Fy.re= -Fy.re; Fz.re= -Fz.re;
			      }

			      phase_phi = TWOPI * 
				   (cur_kvector.x 
				    * (p.x/geometry_lattice.size.x+0.5) +
				    cur_kvector.y 
				    * (p.y/geometry_lattice.size.y+0.5) +
				    cur_kvector.z 
				    * (p.z/geometry_lattice.size.z+0.5));
			      CASSIGN_SCALAR(phase, 
					     cos(phase_phi), sin(phase_phi));
			      CASSIGN_MULT_RE(Fre.x, Fx, phase);
			      CASSIGN_MULT_IM(Fim.x, Fx, phase);
			      CASSIGN_MULT_RE(Fre.y, Fy, phase);
			      CASSIGN_MULT_IM(Fim.y, Fy, phase);
			      CASSIGN_MULT_RE(Fre.z, Fz, phase);
			      CASSIGN_MULT_IM(Fim.z, Fz, phase);

			      integral +=
				   ctl_convert_number_to_c(
					gh_call4(f,
					    ctl_convert_vector3_to_scm(Fre),
				            ctl_convert_vector3_to_scm(Fim),
				            ctl_convert_number_to_scm(epsilon),
                                            ctl_convert_vector3_to_scm(p)));
			 }
		    }
	       }
#endif
	  }
     }

     mpi_allreduce_1(&integral, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     return integral;
}

/**************************************************************************/
