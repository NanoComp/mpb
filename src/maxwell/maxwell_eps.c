/* Copyright (C) 1999, 2000 Massachusetts Institute of Technology.
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
extern void maxwell_sym_matrix_invert(symmetric_matrix *Vinv, 
				      const symmetric_matrix *V)
{
     real m00 = V->m00, m11 = V->m11, m22 = V->m22,
	  m01 = V->m01, m02 = V->m02, m12 = V->m12;

     if (m01 == 0.0 && m02 == 0.0 && m12 == 0.0) {
	  /* optimize common case of a diagonal matrix: */
	  Vinv->m00 = 1.0 / m00;
	  Vinv->m11 = 1.0 / m11;
	  Vinv->m22 = 1.0 / m22;
	  Vinv->m01 = Vinv->m02 = Vinv->m12 = 0.0;
     }
     else {
	  double detinv;
	  
	  /* compute the determinant: */
	  detinv = m00*m11*m22 - m02*m11*m02 + 2.0 * m01*m12*m02 -
	       m01*m01*m22 - m12*m12*m00;
	  
	  CHECK(detinv != 0.0, "singular 3x3 matrix");
	  
	  detinv = 1.0/detinv;
	  
	  Vinv->m00 = detinv * (m11*m22 - m12*m12);
	  Vinv->m11 = detinv * (m00*m22 - m02*m02);
	  Vinv->m22 = detinv * (m11*m00 - m01*m01);
	  
	  Vinv->m02 = detinv * (m01*m12 - m11*m02);
	  Vinv->m01 = detinv * (m12*m02 - m01*m22);
	  Vinv->m12 = detinv * (m01*m02 - m00*m12);
     }
}

#define K_PI 3.141592653589793238462643383279502884197
#define SMALL 1.0e-6
#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))

#define MAX_MOMENT_MESH 12 /* max # of moment-mesh vectors */
#define MOMENT_MESH_R 0.5

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
	of points on the radius MOMENT_MESH_R sphere. */
     if (rank == 1) {
	  /* one-dimensional: just two points (really, normal
	     vector is always along x): */
	  *size_moment_mesh = 2;
	  moment_mesh[0][0] = -MOMENT_MESH_R;
	  moment_mesh[0][1] = 0.0;
	  moment_mesh[0][2] = 0.0;
	  moment_mesh[1][0] = MOMENT_MESH_R;
	  moment_mesh[1][1] = 0.0;
	  moment_mesh[1][2] = 0.0;
     }
     else if (rank == 2) {
	  /* two-dimensional: 6 points at 60-degree intervals: */
	  *size_moment_mesh = 6;
	  for (i = 0; i < *size_moment_mesh; ++i) {
	       moment_mesh[i][0] = cos(i * K_PI / 3.0) * MOMENT_MESH_R;
	       moment_mesh[i][1] = sin(i * K_PI / 3.0) * MOMENT_MESH_R;
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
	       moment_mesh[i][(shift + 1) % 3] = MOMENT_MESH_R *
		    golden * (1 - 2 * (i & 1)) / sqrt(2 - golden);
	       moment_mesh[i][(shift + 2) % 3] = MOMENT_MESH_R *
		    1.0 * (1 - (i & 2)) / sqrt(2 - golden);
	  }
     }

     CHECK(*size_moment_mesh <= MAX_MOMENT_MESH, "yikes, moment mesh too big");

     /* scale the moment-mesh vectors so that the sphere has a
	diameter of 2*MOMENT_MESH_R times the diameter of the
	smallest grid direction: */

     /* first, find the minimum distance between grid points along the
	lattice directions: */
     for (i = 0; i < rank; ++i) {
	  real ri = R[i][0] * R[i][0] + R[i][1] * R[i][1] + R[i][2] * R[i][2];
	  ri = sqrt(ri) / (i == 0 ? nx : (i == 1 ? ny : nz));
	  min_diam = MIN2(min_diam, ri);
     }
     
     /* scale moment_mesh by this diameter: */
     for (i = 0; i < *size_moment_mesh; ++i) {
	  real len = 0;
	  for (j = 0; j < 3; ++j) {
	       moment_mesh[i][j] *= min_diam;
	       len += moment_mesh[i][j] * moment_mesh[i][j];
	       mesh_total[j] += moment_mesh[i][j];
	  }
	  CHECK(fabs(len - min_diam*min_diam*(MOMENT_MESH_R*MOMENT_MESH_R)) 
		< SMALL,
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
   using the dielectric function epsilon(&eps, &eps_inv, r, epsilon_data).

   epsilon is averaged over a rectangular mesh spanning the space between
   grid points; the size of the mesh is given by mesh_size.

   R1, R2, and R2 are the spatial lattice vectors, and are used to convert
   the discretization grid into spatial coordinates (with the origin at
   the (0,0,0) grid element.

   In most places, the dielectric tensor is equal to eps_inv, but at
   dielectric interfaces it varies depending upon the polarization of
   the field (for faster convergence).  In particular, it depends upon
   the direction of the field relative to the surface normal vector,
   so we must compute the latter.  The surface normal is approximated
   by the "dipole moment" of the dielectric function over a spherical
   mesh.

   Implementation note: md->eps_inv is chosen to have dimensions matching
   the output of the FFT.  Thus, its dimensions depend upon whether we are
   doing a real or complex and serial or parallel FFT. */

void set_maxwell_dielectric(maxwell_data *md,
			    const int mesh_size[3],
			    real R[3][3], real G[3][3],
			    maxwell_dielectric_function epsilon,
			    void *epsilon_data)
{
     real s1, s2, s3, m1, m2, m3;  /* grid/mesh steps */
     real mesh_center[3];
     real moment_mesh[MAX_MOMENT_MESH][3];
     real eps_inv_total = 0.0;
     int i, j, k;
     int mesh_prod;
     real mesh_prod_inv;
     int size_moment_mesh = 0;
     int nx, ny, nz;
#ifdef HAVE_MPI
     int local_ny, local_y_start;
#endif

     nx = md->nx; ny = md->ny; nz = md->nz;

     get_mesh(nx, ny, nz, mesh_size, R, G,
	      mesh_center, &mesh_prod, moment_mesh, &size_moment_mesh);
     mesh_prod_inv = 1.0 / mesh_prod;

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
	  symmetric_matrix eps_mean = {0,0,0,0,0,0}, 
			   eps_inv_mean = {0,0,0,0,0,0}, eps_mean_inv;
	  real norm_len;
	  real norm0, norm1, norm2;
	  short means_different_p, diag_eps_p;
	  
	  for (mi = 0; mi < mesh_size[0]; ++mi)
	       for (mj = 0; mj < mesh_size[1]; ++mj)
		    for (mk = 0; mk < mesh_size[2]; ++mk) {
			 real r[3];
			 symmetric_matrix eps, eps_inv;
			 r[0] = i2 * s1 + (mi - mesh_center[0]) * m1;
			 r[1] = j2 * s2 + (mj - mesh_center[1]) * m2;
			 r[2] = k2 * s3 + (mk - mesh_center[2]) * m3;
			 epsilon(&eps, &eps_inv, r, epsilon_data);
			 eps_mean.m00 += eps.m00;
			 eps_mean.m11 += eps.m11;
			 eps_mean.m22 += eps.m22;
			 eps_mean.m01 += eps.m01;
			 eps_mean.m02 += eps.m02;
			 eps_mean.m12 += eps.m12;
			 eps_inv_mean.m00 += eps_inv.m00;
			 eps_inv_mean.m11 += eps_inv.m11;
			 eps_inv_mean.m22 += eps_inv.m22;
			 eps_inv_mean.m01 += eps_inv.m01;
			 eps_inv_mean.m02 += eps_inv.m02;
			 eps_inv_mean.m12 += eps_inv.m12;
		    }
     
	  diag_eps_p =
	       eps_mean.m01==0.0 && eps_mean.m02==0.0 && eps_mean.m12==0.0;
	  if (diag_eps_p) { /* handle the common case of diagonal matrices: */
	       eps_mean_inv.m00 = mesh_prod / eps_mean.m00;
	       eps_mean_inv.m11 = mesh_prod / eps_mean.m11;
	       eps_mean_inv.m22 = mesh_prod / eps_mean.m22;
	       eps_mean_inv.m01 = eps_mean_inv.m02 = eps_mean_inv.m12 = 0.0;
	       eps_inv_mean.m00 *= mesh_prod_inv;
	       eps_inv_mean.m11 *= mesh_prod_inv;
	       eps_inv_mean.m22 *= mesh_prod_inv;

	       means_different_p = 
		    fabs(eps_mean_inv.m00 - eps_inv_mean.m00) > SMALL ||
		    fabs(eps_mean_inv.m11 - eps_inv_mean.m11) > SMALL ||
		    fabs(eps_mean_inv.m22 - eps_inv_mean.m22) > SMALL;
	  }
	  else {
	       eps_mean.m00 *= mesh_prod_inv;
	       eps_mean.m11 *= mesh_prod_inv;
	       eps_mean.m22 *= mesh_prod_inv;
	       eps_mean.m01 *= mesh_prod_inv;
	       eps_mean.m02 *= mesh_prod_inv;
	       eps_mean.m12 *= mesh_prod_inv;
	       maxwell_sym_matrix_invert(&eps_mean_inv, &eps_mean);
	       eps_inv_mean.m00 *= mesh_prod_inv;
	       eps_inv_mean.m11 *= mesh_prod_inv;
	       eps_inv_mean.m22 *= mesh_prod_inv;
	       eps_inv_mean.m01 *= mesh_prod_inv;
	       eps_inv_mean.m02 *= mesh_prod_inv;
	       eps_inv_mean.m12 *= mesh_prod_inv;
	       
	       means_different_p = 
		    fabs(eps_mean_inv.m00 - eps_inv_mean.m00) > SMALL ||
		    fabs(eps_mean_inv.m11 - eps_inv_mean.m11) > SMALL ||
		    fabs(eps_mean_inv.m22 - eps_inv_mean.m22) > SMALL ||
		    fabs(eps_mean_inv.m01 - eps_inv_mean.m01) > SMALL ||
		    fabs(eps_mean_inv.m02 - eps_inv_mean.m02) > SMALL ||
		    fabs(eps_mean_inv.m12 - eps_inv_mean.m12) > SMALL;
	  }

	  /* if the two averaging methods yielded different results,
	     which usually happens if epsilon is not constant, then
	     we need to find the normal vector to the dielectric interface: */
	  if (means_different_p) {
	       real moment0 = 0, moment1 = 0, moment2 = 0;

	       for (mi = 0; mi < size_moment_mesh; ++mi) {
		    real r[3], eps_trace;
		    symmetric_matrix eps, eps_inv;
		    r[0] = i2 * s1 + moment_mesh[mi][0];
		    r[1] = j2 * s2 + moment_mesh[mi][1];
		    r[2] = k2 * s3 + moment_mesh[mi][2];
		    epsilon(&eps, &eps_inv, r, epsilon_data);
		    eps_trace = eps.m00 + eps.m11 + eps.m22;
		    moment0 += eps_trace * moment_mesh[mi][0];
		    moment1 += eps_trace * moment_mesh[mi][1];
		    moment2 += eps_trace * moment_mesh[mi][2];
	       }
	       /* need to convert moment from lattice to cartesian coords: */
	       norm0 = R[0][0]*moment0 + R[1][0]*moment1 + R[2][0]*moment2;
	       norm1 = R[0][1]*moment0 + R[1][1]*moment1 + R[2][1]*moment2;
	       norm2 = R[0][2]*moment0 + R[1][2]*moment1 + R[2][2]*moment2;
	  
	       norm_len = sqrt(norm0*norm0 + norm1*norm1 + norm2*norm2);
	  }
	  
	  if (means_different_p && norm_len > SMALL) {
	       real x0, x1, x2;

	       norm_len = 1.0/norm_len;
	       norm0 *= norm_len;
	       norm1 *= norm_len;
	       norm2 *= norm_len;

	       /* Compute the effective inverse dielectric tensor.
		  We define this as:
                    1/2 ( {eps_inv_mean, P} + {eps_mean_inv, 1-P} )
	          where P is the projection matrix onto the normal direction
                  (P = norm ^ norm), and {a,b} is the anti-commutator ab+ba.
		   = 1/2 {eps_inv_mean - eps_mean_inv, P} + eps_mean_inv
		   = 1/2 (n_i x_j + x_i n_j) + (eps_mean_inv)_ij
		  where n_k is the kth component of the normal vector and 
		     x_i = (eps_inv_mean - eps_mean_inv)_ik n_k  
		  Note the implied summations (Einstein notation).

		  Note that the resulting matrix is symmetric, and we get just
		  eps_inv_mean if eps_inv_mean == eps_mean_inv, as desired.

	          Note that P is idempotent, so for scalar epsilon this
	          is just eps_inv_mean * P + eps_mean_inv * (1-P)
                        = (1/eps_inv_mean * P + eps_mean * (1-P)) ^ (-1),
	          which corresponds to the expression in the Meade paper. */
	       
	       x0 = (eps_inv_mean.m00 - eps_mean_inv.m00) * norm0;
	       x1 = (eps_inv_mean.m11 - eps_mean_inv.m11) * norm1;
	       x2 = (eps_inv_mean.m22 - eps_mean_inv.m22) * norm2;
	       if (diag_eps_p) {
		    md->eps_inv[eps_index].m01 = 0.5*(x0*norm1 + x1*norm0);
		    md->eps_inv[eps_index].m02 = 0.5*(x0*norm2 + x2*norm0);
		    md->eps_inv[eps_index].m12 = 0.5*(x1*norm2 + x2*norm1);
	       }
	       else {
		    x0 += ((eps_inv_mean.m01 - eps_mean_inv.m01) * norm1 + 
			   (eps_inv_mean.m02 - eps_mean_inv.m02) * norm2);
		    x1 += ((eps_inv_mean.m01 - eps_mean_inv.m01) * norm0 + 
			   (eps_inv_mean.m12 - eps_mean_inv.m12) * norm2);
		    x2 += ((eps_inv_mean.m02 - eps_mean_inv.m02) * norm0 +
			   (eps_inv_mean.m12 - eps_mean_inv.m12) * norm1);

		    md->eps_inv[eps_index].m01 = (0.5*(x0*norm1 + x1*norm0) 
						  + eps_mean_inv.m01);
		    md->eps_inv[eps_index].m02 = (0.5*(x0*norm2 + x2*norm0) 
						  + eps_mean_inv.m02);
		    md->eps_inv[eps_index].m12 = (0.5*(x1*norm2 + x2*norm1) 
						  + eps_mean_inv.m12);
	       }
	       md->eps_inv[eps_index].m00 = x0*norm0 + eps_mean_inv.m00;
	       md->eps_inv[eps_index].m11 = x1*norm1 + eps_mean_inv.m11;
	       md->eps_inv[eps_index].m22 = x2*norm2 + eps_mean_inv.m22;
	  }
	  else { /* undetermined normal vector and/or constant eps */
	       md->eps_inv[eps_index] = eps_mean_inv;
	  }
	  
	  eps_inv_total += (md->eps_inv[eps_index].m00 + 
			    md->eps_inv[eps_index].m11 + 
			    md->eps_inv[eps_index].m22);
     }  /* end of loop body */
     
     md->eps_inv_mean = eps_inv_total / (3 * md->fft_output_size);
}
