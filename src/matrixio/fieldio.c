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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <check.h>
#include <matrices.h>

#include "matrixio.h"

#define TWOPI 6.2831853071795864769252867665590057683943388
#define MAX2(a,b) ((a) > (b) ? (a) : (b))

/* note that kvector here is given in the cartesian basis */
void fieldio_write_complex_field(scalar_complex *field,
				 int rank,
				 const int dims[3],
				 int local_nx, int local_x_start,
				 const int copies[3],
				 const real kvector[3],
				 real R[3][3],
				 const char *fname,
				 const char *description)
{
     int i, j, k, component, ri_part;
     int cpies[3];
     int total_dims[3], local_dims[3], start[3] = {0,0,0}, localN;
     real s[3]; /* the step size between grid points dotted with k */
     char *fname2, name[] = "imaginary part of x component";
     matrixio_id file_ids[3][2], data_ids[3][2];
     real lastphase = 0.0;

     for (i = 0; i < 3; ++i) {
	  cpies[i] = MAX2(copies[i], 1); /* make sure copies are non-zero */
	  total_dims[i] = dims[i] * cpies[i];
	  local_dims[i] = dims[i];
     }
     rank = total_dims[2] == 1 ? (total_dims[1] == 1 ? 1 : 2) : 3;
     local_dims[0] = local_nx;
     localN = local_dims[0] * local_dims[1] * local_dims[2];

     for (i = 0; i < 3; ++i) {
	  s[i] = 0;
	  for (j = 0; j < 3; ++j)
	       s[i] += R[i][j] * kvector[j];
	  s[i] *= TWOPI / dims[i];
     }

     /* allocate string for holding fname + ".x.r" or similar */
     fname2 = (char *) malloc(strlen(fname) + 4 + 1);
     CHECK(fname2, "out of memory!");
     strcpy(fname2, fname);
     strcat(fname2, ".x.r");

     /* create output files & data sets for writing */
     for (component = 0; component < 3; ++component)
	  for (ri_part = 0; ri_part < 2; ++ri_part) {
	       fname2[strlen(fname) + 1] = 'x' + component;
	       fname2[strlen(fname) + 3] = ri_part ? 'i' : 'r';
	       sprintf(name, "%s part of %c component",
		       ri_part ? "imag." : "real", 'x' + component);
	       file_ids[component][ri_part] =
		    matrixio_create(fname2);
	       data_ids[component][ri_part] =
		    matrixio_create_dataset(file_ids[component][ri_part],
					    name, description,
					    rank, total_dims);
	  }
     free(fname2);

     /* Now, multiply field by exp(i k*r): */
     for (i = 0; i < local_nx; ++i)
	  for (j = 0; j < dims[1]; ++j) 
	       for (k = 0; k < dims[2]; ++k) {
		    int ijk = ((i * dims[1] + j) * dims[2] + k) * 3;
		    real phase = s[0]*(i+local_x_start) + s[1]*j + s[2]*k;
		    double c = cos(phase), s = sin(phase);
		    real re, im;
		    
		    re = field[ijk].re; im = field[ijk].im;
		    field[ijk].re = re * c - im * s;
		    field[ijk].im = im * c + re * s;

		    re = field[ijk+1].re; im = field[ijk+1].im;
		    field[ijk+1].re = re * c - im * s;
		    field[ijk+1].im = im * c + re * s;

		    re = field[ijk+2].re; im = field[ijk+2].im;
		    field[ijk+2].re = re * c - im * s;
		    field[ijk+2].im = im * c + re * s;
	       }

     /* loop over copies, multiplying by phases and writing out hyperslabs: */
     for (i = 0; i < cpies[0]; ++i)
	  for (j = 0; j < cpies[1]; ++j)
	       for (k = 0; k < cpies[2]; ++k) {
		    int n;
		    real phase = 
			 s[0]*i*dims[0] + s[1]*j*dims[1] + s[2]*k*dims[2];
		    double 
			 c = cos(phase - lastphase), 
			 s = sin(phase - lastphase);

		    /* keep track of previous phase, so that we
		       can subtract it from the exponent to undo the
		       phase change of the previous loop iteration */
		    lastphase = phase;

		    /* multiply by overall phase of this copy: */
		    if (phase != 0.0)
			 for (n = 0; n < localN; ++n) {
			      int ijk = n * 3;
			      real re, im;
			      
			      re = field[ijk].re; im = field[ijk].im;
			      field[ijk].re = re * c - im * s;
			      field[ijk].im = im * c + re * s;
			      
			      re = field[ijk+1].re; im = field[ijk+1].im;
			      field[ijk+1].re = re * c - im * s;
			      field[ijk+1].im = im * c + re * s;
			      
			      re = field[ijk+2].re; im = field[ijk+2].im;
			      field[ijk+2].re = re * c - im * s;
			      field[ijk+2].im = im * c + re * s;
			 }

		    /* start[] is the beginning of this hyperslab: */
		    start[0] = local_x_start + i * dims[0];
		    start[1] = j * dims[1];
		    start[2] = k * dims[2];

		    /* now, write out the hyperslab for each component
		       and for both the real and imaginary parts: */
		    for (component = 0; component < 3; ++component)
			 for (ri_part = 0; ri_part < 2; ++ri_part)
			      matrixio_write_real_data(
				   data_ids[component][ri_part],
				   local_dims, start, 6,
				   ri_part ? &field[component].im
				   : &field[component].re);
	       }

     /* close data sets and files */
     for (component = 0; component < 3; ++component)
	  for (ri_part = 0; ri_part < 2; ++ri_part) {
	       matrixio_close_dataset(data_ids[component][ri_part]);
	       matrixio_close(file_ids[component][ri_part]);
	  }
}

void fieldio_write_real_vals(real *vals,
			     int rank,
			     const int dims[3],
			     int local_nx, int local_x_start,
			     const int copies[3],
			     const char *fname,
			     const char *description)
{
     int i, j, k;
     int cpies[3], total_dims[3], local_dims[3], start[3] = {0,0,0};
     matrixio_id file_id, data_id;

     for (i = 0; i < 3; ++i) {
	  cpies[i] = MAX2(copies[i], 1); /* make sure copies are non-zero */
	  total_dims[i] = dims[i] * cpies[i];
	  local_dims[i] = dims[i];
     }
     rank = total_dims[2] == 1 ? (total_dims[1] == 1 ? 1 : 2) : 3;
     local_dims[0] = local_nx;

     /* create output file & data set for writing */
     file_id = matrixio_create(fname);
     data_id = matrixio_create_dataset(file_id, description, NULL,
				       rank, total_dims);

     /* loop over copies, writing out hyperslabs: */
     for (i = 0; i < cpies[0]; ++i)
	  for (j = 0; j < cpies[1]; ++j)
	       for (k = 0; k < cpies[2]; ++k) {
		    /* start[] is the beginning of this hyperslab: */
		    start[0] = local_x_start + i * dims[0];
		    start[1] = j * dims[1];
		    start[2] = k * dims[2];
		    matrixio_write_real_data(data_id,local_dims,start,1,vals);
       }

     /* close data set and file */
     matrixio_close_dataset(data_id);
     matrixio_close(file_id);
}
