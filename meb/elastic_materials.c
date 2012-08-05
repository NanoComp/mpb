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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "config.h"
#include <check.h>
#include <mpiglue.h>
#include <mpi_utils.h>

#include "elastic.h"

#define MAX2(a,b) ((a) > (b) ? (a) : (b))

/**************************************************************************/

/* The following function initializes the elastic materials arrays:
                 rho, sqrt_rhoinv, rhoct2, rhocl2
   in md using the function mat(&rho, &ct, &cl, r, mat_data).

   The quantities are averaged over a rectangular mesh spanning the space
   between grid points; the size of the mesh is given by mesh_size.

   R[0..2] are the spatial lattice vectors, and are used to convert
   the discretization grid into spatial coordinates (with the origin
   at the (0,0,0) grid element.

   In the future, we may expand this to compute effective material
   tensors instead of scalar averages, similar to MPB.

   Implementation note: rho etc. are chosen to have dimensions matching
   the output of the FFT.  Thus, their dimensions depend upon whether we are
   doing a real or complex and serial or parallel FFT. */

void set_elastic_materials(elastic_data *md,
			    const int mesh_size[3],
			    real R[3][3], real G[3][3],
			    elastic_material_function mat,
			    void *mat_data)
{
     real s1, s2, s3, m1, m2, m3;  /* grid/mesh steps */
     real mesh_center[3];
     int i, j, k;
     int mesh_prod;
     real mesh_prod_inv;
     int n1, n2, n3;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
#ifndef SCALAR_COMPLEX
     int n_other, n_last, rank;
#endif

     n1 = md->nx; n2 = md->ny; n3 = md->nz;

     mesh_prod = 1;
     for (i = 0; i < 3; ++i) {
	  int ms = MAX2(mesh_size[i], 1);
	  mesh_center[i] = (ms - 1) * 0.5;
	  mesh_prod *= ms;
     }
     mesh_prod_inv = 1.0 / mesh_prod;

     s1 = 1.0 / n1;
     s2 = 1.0 / n2;
     s3 = 1.0 / n3;
     m1 = s1 / MAX2(1, mesh_size[0]);
     m2 = s2 / MAX2(1, mesh_size[1]);
     m3 = s3 / MAX2(1, mesh_size[2]);

     /* Here we have different loops over the coordinates, depending
	upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and mat_index describing the corresponding index in 
	the array md->rho[]. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI
     
     for (i = 0; i < n1; ++i)
	  for (j = 0; j < n2; ++j)
	       for (k = 0; k < n3; ++k)
     {
#         define i2 i
#         define j2 j
#         define k2 k
	  int mat_index = ((i * n2 + j) * n3 + k);

#  else /* HAVE_MPI */

     local_n2 = md->local_ny;
     local_y_start = md->local_y_start;

     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < n3; ++k)
     {
#         define i2 i
	  int j2 = j + local_y_start;
#         define k2 k
	  int mat_index = ((j * n1 + i) * n3 + k);

#  endif /* HAVE_MPI */

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

     n_other = md->other_dims;
     n_last = md->last_dim_size / 2;
     rank = (n3 == 1) ? (n2 == 1 ? 1 : 2) : 3;

     for (i = 0; i < n_other; ++i)
	  for (j = 0; j < n_last; ++j)
     {
	  int mat_index = i * n_last + j;
	  int i2, j2, k2;
	  switch (rank) {
	      case 2: i2 = i; j2 = j; k2 = 0; break;
	      case 3: i2 = i / n2; j2 = i % n2; k2 = j; break;
	      default: i2 = j; j2 = k2 = 0;  break;
	  }

#  else /* HAVE_MPI */

     local_n2 = md->local_ny;
     local_y_start = md->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = md->last_dim_size / 2;
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
	  int mat_index = ((j * n1 + i) * local_n3 + k);

#  endif  /* HAVE_MPI */

#endif /* not SCALAR_COMPLEX */

     {
	  int mi, mj, mk;
	  real rho_mean = 0, sqrt_rhoinv_mean = 0, rhoct2_mean = 0, rhocl2_mean = 0;
	  
	  for (mi = 0; mi < mesh_size[0]; ++mi)
	       for (mj = 0; mj < mesh_size[1]; ++mj)
		    for (mk = 0; mk < mesh_size[2]; ++mk) {
			 real r[3];
			 real rho, cl, ct;
			 r[0] = i2 * s1 + (mi - mesh_center[0]) * m1;
			 r[1] = j2 * s2 + (mj - mesh_center[1]) * m2;
			 r[2] = k2 * s3 + (mk - mesh_center[2]) * m3;
			 mat(&rho, &ct, &cl, r, mat_data);
			 rho_mean += rho;
			 sqrt_rhoinv_mean += 1.0 / sqrt(rho);
			 rhoct2_mean += rho * ct*ct;
			 rhocl2_mean += rho * cl*cl;
		    }
     
	  md->rho[mat_index] = rho_mean * mesh_prod_inv;
	  md->sqrt_rhoinv[mat_index] = sqrt_rhoinv_mean * mesh_prod_inv;
	  md->rhoct2[mat_index] = rhoct2_mean * mesh_prod_inv;
	  md->rhocl2[mat_index] = rhocl2_mean * mesh_prod_inv;
     }}  /* end of loop body */
}
