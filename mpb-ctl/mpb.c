/* Copyright (C) 1999 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/**************************************************************************/

/* Here, we define the external functions callable from Guile, as defined
   by photon.scm. */

/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* GNU Guile library header file: */
#include <guile/gh.h>

/* Header files for my eigensolver routines: */
#include <config.h>
#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <eigensolver.h>
#include <maxwell.h>

/* Header file for the ctl-file (Guile) interface; automatically
   generated from photon.scm */
#include <ctl-io.h>

/* Routines from libctl/utils/libgeom: */
#include <geom.h>

/**************************************************************************/

/* function to display a little information about a geometric object to
   convince the user that we've read it in correctly. */
static void display_object_info(geometric_object obj)
{
     printf("     object center = (%g,%g,%g), epsilon = %g\n",
	    obj.center.x, obj.center.y, obj.center.z,
	    obj.material.epsilon);    
     switch (obj.which_subclass) {
	 case CYLINDER:
	      printf("          cylinder with height %g, axis (%g, %g, %g)\n",
		     obj.subclass.cylinder_data->height,
		     obj.subclass.cylinder_data->axis.x,
		     obj.subclass.cylinder_data->axis.y,
		     obj.subclass.cylinder_data->axis.z);
	      break;
	 case SPHERE:
	      printf("          sphere with radius %g\n",
		     obj.subclass.sphere_data->radius);
	      break;
	 case BLOCK:
	      printf("          %s with size (%g,%g,%g)\n",
		     obj.subclass.block_data->which_subclass == BLOCK_SELF 
		     ? "block" : "ellipsoid",
		     obj.subclass.block_data->size.x,
		     obj.subclass.block_data->size.y,
		     obj.subclass.block_data->size.z);
	      printf("          projection matrix: %10.6f%10.6f%10.6f\n"
		     "                             %10.6f%10.6f%10.6f\n"
		     "                             %10.6f%10.6f%10.6f\n",
		     obj.subclass.block_data->projection_matrix.c0.x,
		     obj.subclass.block_data->projection_matrix.c1.x,
		     obj.subclass.block_data->projection_matrix.c2.x,
		     obj.subclass.block_data->projection_matrix.c0.y,
		     obj.subclass.block_data->projection_matrix.c1.y,
		     obj.subclass.block_data->projection_matrix.c2.y,
		     obj.subclass.block_data->projection_matrix.c0.z,
		     obj.subclass.block_data->projection_matrix.c1.z,
		     obj.subclass.block_data->projection_matrix.c2.z);
	      break;
	 case GEOMETRIC_OBJECT_SELF:
	      printf("          generic geometric object\n");
	      break;
	 default:
	      printf("          UNKNOWN OBJECT TYPE!\n");
     }
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

#define NWORK 3
#define NUM_FFT_BANDS 20 /* max number of bands to FFT at a time */

/* global variables for retaining data about the eigenvectors between
   calls from Guile: */

static maxwell_data *mdata = NULL;
static maxwell_target_data *mtdata = NULL;
static evectmatrix H, W[NWORK];

/* R[i]/G[i] are lattice/reciprocal-lattice vectors */
static real R[3][3], G[3][3];
static matrix3x3 Rm, Gm; /* same thing, but matrix3x3 */

/**************************************************************************/

/* Given a position r in the basis of the lattice vectors, return the
   corresponding dielectric constant.  edata is ignored.  Should be
   called from within init_params (or after init_params), so that the
   geometry input variables will have been read in (for use in libgeom).

   This function is passed to set_maxwell_dielectric to initialize
   the dielectric tensor array for eigenvector calculations. */

static real epsilon_func(real r[3], void *edata)
{
     material_type material;
     vector3 p;

     /* p needs to be in the lattice *unit* vector basis, while r is
	in the lattice vector basis; convert: */
     p.x = r[0] * geometry_lattice.size.x;
     p.y = r[1] * geometry_lattice.size.y;
     p.z = r[2] * geometry_lattice.size.z;
     
     material = material_of_point(p);  /* from libctl/utils/libgeom/geom.c */
     return material.epsilon;
}

/**************************************************************************/

/* Guile-callable function: init-params, which initializes any data
   that we need for the eigenvalue calculation.  When this function
   is called, the input variables (the geometry, etcetera) have already
   been read into the global variables defined in ctl-io.h.  */
void init_params(void)
{
     int i, local_N, N_start, alloc_N;
     int nx, ny, nz;
     int mesh[3];
     
     /* Output a bunch of stuff so that the user can see what we're
	doing and what we've read in. */
     
     printf("@@@@@ init-params: initializing eigensolver data\n");
     
     printf("Computing %d bands with %e tolerance.\n", num_bands, tolerance);
     if (target_freq != 0.0)
	  printf("Target frequency is %g\n", target_freq);
     
     nx = grid_size.x;
     ny = grid_size.y;
     nz = grid_size.z;

     {
	  int true_rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);
	  if (true_rank < dimensions)
	       dimensions = true_rank;
     }

     printf("Working in %d dimensions.\n", dimensions);

     printf("Grid size is %d x %d x %d.\n", nx, ny, nz);

     printf("Mesh size is %d.\n", mesh_size);
     mesh[0] = mesh_size;
     mesh[1] = (dimensions > 1) ? mesh_size : 1;
     mesh[2] = (dimensions > 2) ? mesh_size : 1;

     Rm.c0 = vector3_scale(geometry_lattice.size.x, geometry_lattice.basis.c0);
     Rm.c1 = vector3_scale(geometry_lattice.size.y, geometry_lattice.basis.c1);
     Rm.c2 = vector3_scale(geometry_lattice.size.z, geometry_lattice.basis.c2);
     printf("Lattice vectors:\n");
     printf("     (%g, %g, %g)\n", Rm.c0.x, Rm.c0.y, Rm.c0.z);  
     printf("     (%g, %g, %g)\n", Rm.c1.x, Rm.c1.y, Rm.c1.z);
     printf("     (%g, %g, %g)\n", Rm.c2.x, Rm.c2.y, Rm.c2.z);
  
     Gm = matrix3x3_inverse(Rm);
     printf("Reciprocal lattice vectors:\n");
     printf("     (%g, %g, %g)\n", Gm.c0.x, Gm.c0.y, Gm.c0.z);  
     printf("     (%g, %g, %g)\n", Gm.c1.x, Gm.c1.y, Gm.c1.z);
     printf("     (%g, %g, %g)\n", Gm.c2.x, Gm.c2.y, Gm.c2.z);
     
     matrix3x3_to_arr(R, Rm);
     matrix3x3_to_arr(G, Gm);

     printf("\nGeometric objects:\n");
     for (i = 0; i < geometry.num_items; ++i)
	  display_object_info(geometry.items[i]);

     printf("\n%d k-points:\n", k_points.num_items);
     for (i = 0; i < k_points.num_items; ++i)
	  printf("     (%g,%g,%g)\n", k_points.items[i].x,
		 k_points.items[i].y, k_points.items[i].z);
     
     if (mdata) {  /* need to clean up from previous init_params call */
	  destroy_evectmatrix(H);
	  for (i = 0; i < NWORK; ++i)
	       destroy_evectmatrix(W[i]);
	  destroy_maxwell_target_data(mtdata); mtdata = NULL;
	  destroy_maxwell_data(mdata); mdata = NULL;
     }
     else
	  srand(time(NULL)); /* init random seed for field initialization */
   
     printf("Creating Maxwell data...\n");
     mdata = create_maxwell_data(nx, ny, nz, &local_N, &N_start, &alloc_N,
                                 num_bands, NUM_FFT_BANDS);
     CHECK(mdata, "NULL mdata");

     /* can also be TE_POLARIZATION or TM_POLARIZATION; what is the
	best way to set this? */
     set_maxwell_data_polarization(mdata, NO_POLARIZATION);

     printf("Initializing dielectric function...\n");
     set_maxwell_dielectric(mdata, mesh, R, epsilon_func, NULL);

     if (target_freq != 0.0)
	  mtdata = create_maxwell_target_data(mdata, target_freq);
     else
	  mtdata = NULL;

     printf("Allocating fields...\n");
     H = create_evectmatrix(nx * ny * nz, 2, num_bands,
                            local_N, N_start, alloc_N);
     for (i = 0; i < NWORK; ++i)
          W[i] = create_evectmatrix(nx * ny * nz, 2, num_bands,
                                    local_N, N_start, alloc_N);

     printf("Initializing fields to random numbers...\n");
     for (i = 0; i < H.n * H.p; ++i)
          ASSIGN_REAL(H.data[i], rand() * 1.0 / RAND_MAX);

     printf("Stuff for grepping:\n");
     printf("sumfrq:, k index, kx, ky, kz, kmag/2pi");
     for (i = 0; i < num_bands; ++i)
	  printf(", band %d", i + 1);
     printf("\n");

     printf("@@@@@\n");
}

/* Solve for the bands at a given k point.
   Must only be called after init_params! */
void solve_kpoint(vector3 kvector)
{
     static int call_count = 0;
     int i, num_iters;
     real *eigvals;
     real k[3];

     printf("@@@@@ solve_kpoint (%g,%g,%g):\n",
	    kvector.x, kvector.y, kvector.z);

     vector3_to_arr(k, kvector);
     update_maxwell_data_k(mdata, k, G[0], G[1], G[2]);

     eigvals = (real*) malloc(sizeof(real) * num_bands);
     CHECK(eigvals, "out of memory");

     printf("Solving for bands...\n");

     if (mtdata) {  /* solving for bands near a target frequency */
	  eigensolver(H, eigvals,
		      maxwell_target_operator, (void *) mtdata,
		      maxwell_target_preconditioner, (void *) mtdata, NULL,
		      maxwell_constraint, (void *) mdata,
		      W, NWORK, tolerance, &num_iters);
	  /* now, diagonalize the real Maxwell operator in the
	     solution subspace to get the true eigenvalues and
	     eigenvectors: */
	  eigensolver_get_eigenvals(H, eigvals, maxwell_operator, mdata,
                                    W[0], W[1]);
     }
     else
	  eigensolver(H, eigvals,
		      maxwell_operator, (void *) mdata,
		      maxwell_preconditioner, (void *) mdata, NULL,
		      maxwell_constraint, (void *) mdata,
		      W, NWORK, tolerance, &num_iters);

     printf("Finished solving for bands after %d iterations.\n", num_iters);
     
     if (num_write_output_vars > 1)
	  destroy_output_vars(); /* we are required by libctl to call this
				    in order to deallocate output vars from
				    previous calls. */

     /* create freqs array for storing frequencies in a Guile list */
     freqs.num_items = num_bands;
     freqs.items = (number *) malloc(freqs.num_items * sizeof(number));
     CHECK(freqs.items, "out of memory");
     
     printf("sumfrq:, %d, %g, %g, %g, %g",
	    ++call_count, k[0], k[1], k[2],
	    vector3_norm(matrix3x3_vector3_mult(Gm, kvector)));
     for (i = 0; i < num_bands; ++i) {
	  freqs.items[i] = sqrt(eigvals[i]);
	  printf(", %g", freqs.items[i]);
     }
     printf("\n");

     free(eigvals);

     printf("@@@@@\n");
}

