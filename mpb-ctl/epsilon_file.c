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

/* This file defines functions to initialize the dielectric function
   by reading epsilon values from a file, rather than using the
   geometry.  Actually, we would like to use the geometry in addition
   to the epsilon file, for added flexibility.  So, we return an epsilon
   function that can be used when no geometric objects are found. */

/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../src/config.h"
#include <check.h>
#include <matrices.h>
#include <matrixio.h>
#include <maxwell.h>

#include "mpb.h"

typedef struct {
     int nx, ny, nz;
     real *data;
} epsilon_file_data;

static void epsilon_file_func(symmetric_matrix *eps, symmetric_matrix *eps_inv,
			      real r[3], void *edata)
{
     epsilon_file_data *d = (epsilon_file_data *) edata;
     real rx, ry, rz, dx, dy, dz;
     int x, y, z, x2, y2, z2;
     real eps_val;

     /* make sure r is positive: */
     rx = r[0] >= 0.0 ? r[0] : (r[0] + (1 + (int) (-r[0])));
     ry = r[1] >= 0.0 ? r[1] : (r[1] + (1 + (int) (-r[1])));
     rz = r[2] >= 0.0 ? r[2] : (r[2] + (1 + (int) (-r[2])));

     /* make sure r is in [0,1) */
     rx = rx < 1.0 ? rx : rx - ((int) rx);
     ry = ry < 1.0 ? ry : ry - ((int) ry);
     rz = rz < 1.0 ? rz : rz - ((int) rz);

     /* get the point corresponding to r in the epsilon array grid: */
     x = rx * d->nx;
     y = ry * d->ny;
     z = rz * d->nz;

     /* get the difference between (x,y,z) and the actual point */
     dx = rx * d->nx - x;
     dy = ry * d->ny - y;
     dz = rz * d->nz - z;
     
     /* get the other closest point in the grid, with periodic boundaries: */
     x2 = (d->nx + (dx >= 0.0 ? x + 1 : x - 1)) % d->nx;
     y2 = (d->ny + (dy >= 0.0 ? y + 1 : y - 1)) % d->ny;
     z2 = (d->nz + (dz >= 0.0 ? z + 1 : z - 1)) % d->nz;

     /* take abs(d{xyz}) to get weights for {xyz} and {xyz}2: */
     dx = fabs(dx);
     dy = fabs(dy);
     dz = fabs(dz);

     /* define a macro to give us epsilon(x,y,z) on the grid,
	in row-major order (the order used by HDF5): */
#define EPS(x,y,z) (d->data[((x)*d->ny + (y))*d->nz + (z)])

     /* compute the effective epsilon by linear interpolation: */
     eps_val = (((EPS(x,y,z)*(1.0-dx) + EPS(x2,y,z)*dx) * (1.0-dy) +
		 (EPS(x,y2,z)*(1.0-dx) + EPS(x2,y2,z)*dx) * dy) * (1.0-dz) +
		((EPS(x,y,z2)*(1.0-dx) + EPS(x2,y,z2)*dx) * (1.0-dy) +
		 (EPS(x,y2,z2)*(1.0-dx) + EPS(x2,y2,z2)*dx) * dy) * dz);

#undef EPS

     eps->m00 = eps->m11 = eps->m22 = eps_val;
     eps->m01 = eps->m02 = eps->m12 = 0.0;
     eps_inv->m00 = eps_inv->m11 = eps_inv->m22 = 1.0 / eps_val;
     eps_inv->m01 = eps_inv->m02 = eps_inv->m12 = 0.0;
}

void get_epsilon_file_func(const char *fname,
			   maxwell_dielectric_function *func,
			   void **func_data)
{
     if (fname && fname[0]) {
	  matrixio_id file_id;
	  epsilon_file_data *d;
	  int rank = 3, dims[3];

	  printf("Using background dielectric from file \"%s\"...\n", fname);

	  CHK_MALLOC(d, epsilon_file_data, 1);
	  
	  file_id = matrixio_open(fname, 1);
	  d->data = matrixio_read_real_data(file_id, NULL, &rank, dims,
					    0,0,0, NULL);
	  CHECK(d->data, "couldn't find dataset in dielectric file");
	  matrixio_close(file_id);
	  
	  d->nx = rank >= 1 ? dims[0] : 1;
	  d->ny = rank >= 2 ? dims[1] : 1;
	  d->nz = rank >= 3 ? dims[2] : 1;

	  printf("    ...read %dx%dx%d dielectric function\n",
		 d->nx, d->ny, d->nz);

	  *func = epsilon_file_func;
	  *func_data = (void*) d;
     }
     else {
	  *func = NULL;
	  *func_data = NULL;
     }
}

void destroy_epsilon_file_func_data(void *func_data)
{
     epsilon_file_data *d = (epsilon_file_data *) func_data;
     if (d) {
	  free(d->data);
	  free(d);
     }
}
