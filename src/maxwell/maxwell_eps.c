/* Copyright (C) 1999 Massachusetts Institute of Technology.
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
#include <math.h>

#include "../config.h"
#include <check.h>

#include "maxwell.h"

/* Set Vinv = inverse of V, where both V and Vinv are symmetric matrices. */
static void sym_matrix_invert(symmetric_matrix *Vinv, symmetric_matrix V)
{
     double detinv;
     
     /* compute the determinant: */
     detinv = V.m00*V.m11*V.m22 - V.m02*V.m11*V.m02 +
	      2.0 * V.m01*V.m12*V.m02 -
	      V.m01*V.m01*V.m22 - V.m12*V.m12*V.m00;

     /* don't bother to check for singular matrices, as that shouldn't
	be possible in the context in which this is used */
     detinv = 1.0/detinv;
     
     Vinv->m00 = detinv * (V.m11*V.m22 - V.m12*V.m12);
     Vinv->m11 = detinv * (V.m00*V.m22 - V.m02*V.m02);
     Vinv->m22 = detinv * (V.m11*V.m00 - V.m01*V.m01);
     
     Vinv->m02 = detinv * (V.m01*V.m12 - V.m11*V.m02);
     Vinv->m01 = -detinv * (V.m01*V.m22 - V.m12*V.m02);
     Vinv->m12 = -detinv * (V.m00*V.m12 - V.m01*V.m02);
}

#define K_PI 3.141592653589793238462643383279502884197
#define MAX_MOMENT_MESH 12 /* max # of moment-mesh vectors */
#define SMALL 1.0e-6
#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))

/* A function to set up the mesh given the grid dimensions, mesh size,
   and lattice/reciprocal vectors.  (Any mesh sizes < 1 are taken to
   be 1.)  The returned values are:

   mesh_center: the center of the mesh, relative to the integer
                mesh coordinates; e.g. the mesh_center for a 3x3x3
                mesh is the point (1,1,1). 
   mesh_prod: the product of the mesh sizes.
   moment_mesh: an array of size_moment_mesh vectors, in lattice
                coordinates, of a spherically-symmetric mesh of
                points centered on the origin, designed to be
		used for averaging the first moment of epsilon at 
		a grid point (for finding the local surface normal). */
static void get_mesh(int nx, int ny, int nz, const int mesh_size[3],
		     real R[3][3], real G[3][3],
		     real mesh_center[3], int *mesh_prod,
		     real moment_mesh[MAX_MOMENT_MESH][3],
		     int *size_moment_mesh)
{
     int i,j;
     real min_diam = 1e20;
     real mesh_total[3] = { 0, 0, 0 };
     int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);
	  
     *mesh_prod = 1;
     for (i = 0; i < 3; ++i) {
	  int ms = MAX2(mesh_size[i], 1);
	  mesh_center[i] = (ms - 1) * 0.5;
	  *mesh_prod *= ms;
     }

     /* Set the mesh to compute the "dipole moment" of epsilon over,
	in Cartesian coordinates (used to compute the normal vector
	at surfaces).  Ideally, a spherically-symmetrical distribution
	of points on the radius-0.5 sphere. */
     if (rank == 1) {
	  /* one-dimensional: just two points (really, normal
	     vector is always along x): */
	  *size_moment_mesh = 2;
	  moment_mesh[0][0] = -0.5;
	  moment_mesh[0][1] = 0.0;
	  moment_mesh[0][2] = 0.0;
	  moment_mesh[1][0] = 0.5;
	  moment_mesh[1][1] = 0.0;
	  moment_mesh[1][2] = 0.0;
     }
     else if (rank == 2) {
	  /* two-dimensional: 6 points at 60-degree intervals: */
	  *size_moment_mesh = 6;
	  for (i = 0; i < *size_moment_mesh; ++i) {
	       moment_mesh[i][0] = cos(i * K_PI / 3.0) * 0.5;
	       moment_mesh[i][1] = sin(i * K_PI / 3.0) * 0.5;
	       moment_mesh[i][2] = 0.0;
	  }
     }
     else {
	  /* three-dimensional: use the vertices of an icosahedron,
	     which are at angles of acos(1/sqrt(5)) from each other
	     (roughly 63 degrees).  The vertices are at points (0, +/-
	     g, +/- 1), and cyclic shifts thereof, where g is the
	     golden ratio, scaled by 0.5/sqrt(2-g) onto the
	     unit-diameter sphere (see also the
	     comp.graphics.algorithms FAQ). */
	  real golden = (sqrt(5.0) - 1.0) * 0.5;  /* the golden ratio */
	  *size_moment_mesh = 12;
	  for (i = 0; i < *size_moment_mesh; ++i) {
	       int shift = i / 4;
	       moment_mesh[i][(shift + 0) % 3] = 0.0;
	       moment_mesh[i][(shift + 1) % 3] = 
		    golden * (1 - 2 * (i & 1)) * 0.5 / sqrt(2 - golden);
	       moment_mesh[i][(shift + 2) % 3] = 
		    1.0 * (1 - (i & 2)) * 0.5 / sqrt(2 - golden);
	  }
     }

     CHECK(*size_moment_mesh <= MAX_MOMENT_MESH, "yikes, moment mesh too big");

     /* scale the moment-mesh vectors to the diameter of the smallest
	grid direction: */

     /* first, find the minimum distance between grid points along the
	lattice directions: */
     for (i = 0; i < rank; ++i) {
	  real ri = R[i][0] * R[i][0] + R[i][1] * R[i][1] + R[i][2] * R[i][2];
	  ri = sqrt(ri) / (i == 0 ? nx : (i == 1 ? ny : nz));
	  min_diam = MIN2(min_diam, ri);
     }
     
     /* scale moment_mesh to this diameter: */
     for (i = 0; i < *size_moment_mesh; ++i) {
	  real len = 0;
	  for (j = 0; j < 3; ++j) {
	       moment_mesh[i][j] *= min_diam;
	       len += moment_mesh[i][j] * moment_mesh[i][j];
	       mesh_total[j] += moment_mesh[i][j];
	  }
	  CHECK(fabs(len - min_diam*min_diam/4) < SMALL,
		"bug in get_mesh: moment_mesh vector is wrong length");
     }
     CHECK(fabs(mesh_total[0]) + fabs(mesh_total[1]) + fabs(mesh_total[2])
	   < SMALL, "bug in get_mesh: moment_mesh does not average to zero");

     /* Now, convert the moment_mesh vectors to lattice/grid coordinates;
        to do this, we multiply by the G matrix (inverse of R transposed) */
     for (i = 0; i < *size_moment_mesh; ++i) {
	  real v[3];
	  for (j = 0; j < 3; ++j)    /* make a copy of moment_mesh[i] */
	       v[j] = moment_mesh[i][j];
	  for (j = 0; j < 3; ++j)
	       moment_mesh[i][j] = G[j][0]*v[0] + G[j][1]*v[1] + G[j][2]*v[2];
     }
}

/* The following function initializes the dielectric tensor md->eps_inv,
   using the dielectric function epsilon(r, epsilon_data).

   epsilon is averaged over a rectangular mesh spanning the space between
   grid points; the size of the mesh is given by mesh_size.

   R1, R2, and R2 are the spatial lattice vectors, and are used to convert
   the discretization grid into spatial coordinates (with the origin at
   the (0,0,0) grid element.

   In most places, the dielectric tensor is equal to 1/eps * identity,
   but at dielectric interfaces it varies depending upon the polarization
   of the field (for faster convergence).  In particular, it depends upon
   the direction of the field relative to the surface normal vector, so
   we must compute the latter.  The surface normal is approximated by
   the "dipole moment" of the dielectric function over a spherical mesh.

   Implementation note: md->eps_inv is chosen to have dimensions matching
   the output of the FFT.  Thus, its dimensions depend upon whether we are
   doing a real or complex and serial or parallel FFT. */

void set_maxwell_dielectric(maxwell_data *md,
			    const int mesh_size[3],
			    real R[3][3], real G[3][3],
			    dielectric_function epsilon,
			    void *epsilon_data)
{
     real s1, s2, s3, m1, m2, m3;  /* grid/mesh steps */
     real mesh_center[3];
     real moment_mesh[MAX_MOMENT_MESH][3];
     real eps_inv_total = 0.0;
     int i, j, k;
     int mesh_prod;
     int size_moment_mesh = 0;
     int nx, ny, nz;
#ifdef HAVE_MPI
     int local_ny, local_y_start;
#endif

     nx = md->nx; ny = md->ny; nz = md->nz;

     get_mesh(nx, ny, nz, mesh_size, R, G,
	      mesh_center, &mesh_prod, moment_mesh, &size_moment_mesh);

#ifdef DEBUG
     printf("Using moment mesh:\n");
     for (i = 0; i < size_moment_mesh; ++i)
	  printf("   (%g, %g, %g)\n",
		 moment_mesh[i][0], moment_mesh[i][1], moment_mesh[i][2]);
#endif

     s1 = 1.0 / nx;
     s2 = 1.0 / ny;
     s3 = 1.0 / nz;
     m1 = s1 / MAX2(1, mesh_size[0] - 1);
     m2 = s2 / MAX2(1, mesh_size[1] - 1);
     m3 = s3 / MAX2(1, mesh_size[2] - 1);

     /* Here we have different loops over the coordinates, depending
	upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and eps_index describing the corresponding index in 
	the array md->eps_inv[]. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI
     
     for (i = 0; i < nx; ++i)
	  for (j = 0; j < ny; ++j)
	       for (k = 0; k < nz; ++k)
     {
#         define i2 i
#         define j2 j
#         define k2 k
	  int eps_index = ((i * ny + j) * nz + k);

#  else /* HAVE_MPI */

     local_ny = md->fft_output_size / (nx * nz);
     local_y_start = md->fft_output_N_start / (nx * nz);

     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_ny; ++j)
          for (i = 0; i < nx; ++i)
	       for (k = 0; k < nz; ++k)
     {
#         define i2 i
	  int j2 = j + local_y_start;
#         define k2 k
	  int eps_index = ((j * nx + i) * nz + k);

#  endif

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

#  else /* HAVE_MPI */

#    error not yet implemented!

#  endif

#endif /* not SCALAR_COMPLEX */

	  int mi, mj, mk;
	  real eps_mean = 0.0, eps_inv_mean = 0.0, norm_len;
	  real norm0, norm1, norm2;
	  short means_different_p;
	  
	  for (mi = 0; mi < mesh_size[0]; ++mi)
	       for (mj = 0; mj < mesh_size[1]; ++mj)
		    for (mk = 0; mk < mesh_size[2]; ++mk) {
			 real r[3], eps;			 
			 r[0] = i2 * s1 + (mi - mesh_center[0]) * m1;
			 r[1] = j2 * s2 + (mj - mesh_center[1]) * m2;
			 r[2] = k2 * s3 + (mk - mesh_center[2]) * m3;
			 eps = epsilon(r, epsilon_data);
			 eps_mean += eps;
			 eps_inv_mean += 1.0 / eps;
		    }
     
	  /* compute inverses of epsilon means: */
	  eps_mean = mesh_prod / eps_mean;
	  eps_inv_mean = eps_inv_mean / mesh_prod;

	  means_different_p = fabs(eps_mean - eps_inv_mean) > SMALL;

	  if (means_different_p) {
	       real moment0 = 0, moment1 = 0, moment2 = 0;

	       for (mi = 0; mi < size_moment_mesh; ++mi) {
		    real r[3], eps;
		    r[0] = i2 * s1 + moment_mesh[mi][0];
		    r[1] = j2 * s2 + moment_mesh[mi][1];
		    r[2] = k2 * s3 + moment_mesh[mi][2];
		    eps = epsilon(r, epsilon_data);
		    moment0 += eps * moment_mesh[mi][0];
		    moment1 += eps * moment_mesh[mi][1];
		    moment2 += eps * moment_mesh[mi][2];
	       }
	       /* need to convert moment from lattice to cartesian coords: */
	       norm0 = R[0][0]*moment0 + R[1][0]*moment1 + R[2][0]*moment2;
	       norm1 = R[0][1]*moment0 + R[1][1]*moment1 + R[2][1]*moment2;
	       norm2 = R[0][2]*moment0 + R[1][2]*moment1 + R[2][2]*moment2;
	  
	       norm_len = sqrt(norm0*norm0 + norm1*norm1 + norm2*norm2);
	  }
	  
	  if (means_different_p && norm_len > SMALL) {
	       norm0 /= norm_len;
	       norm1 /= norm_len;
	       norm2 /= norm_len;
	       
	       /* compute effective inverse dielectric tensor: */
	       
	       md->eps_inv[eps_index].m00 = 
		    (eps_inv_mean-eps_mean) * norm0*norm0 + eps_mean;
	       md->eps_inv[eps_index].m11 = 
		    (eps_inv_mean-eps_mean) * norm1*norm1 + eps_mean;
	       md->eps_inv[eps_index].m22 = 
		    (eps_inv_mean-eps_mean) * norm2*norm2 + eps_mean;
	       md->eps_inv[eps_index].m01 = 
		    (eps_inv_mean-eps_mean) * norm0*norm1;
	       md->eps_inv[eps_index].m02 = 
		    (eps_inv_mean-eps_mean) * norm0*norm2;
	       md->eps_inv[eps_index].m12 = 
		    (eps_inv_mean-eps_mean) * norm1*norm2;
	  }
	  else { /* undetermined normal vector and/or constant eps */
	       md->eps_inv[eps_index].m00 = eps_mean;
	       md->eps_inv[eps_index].m11 = eps_mean;
	       md->eps_inv[eps_index].m22 = eps_mean;
	       md->eps_inv[eps_index].m01 = 0.0;
	       md->eps_inv[eps_index].m02 = 0.0;
	       md->eps_inv[eps_index].m12 = 0.0;
	  }
	  
	  eps_inv_total += (md->eps_inv[eps_index].m00 + 
			    md->eps_inv[eps_index].m11 + 
			    md->eps_inv[eps_index].m22);
     }  /* end of loop body */
     
     md->eps_inv_mean = eps_inv_total / (3 * md->fft_output_size);
}
