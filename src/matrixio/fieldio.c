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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../config.h"

#include <check.h>
#include <matrices.h>

#include "matrixio.h"

#define TWOPI 6.2831853071795864769252867665590057683943388
#define MAX2(a,b) ((a) > (b) ? (a) : (b))

/* note that kvector here is given in the cartesian basis */
void fieldio_write_complex_field(scalar_complex *field,
				 int rank,
				 const int dims[3],
				 int which_component,
				 int local_nx, int local_x_start,
				 const int copies[3],
				 const real kvector[3],
				 real R[3][3],
				 matrixio_id file_id)
{
     int i, j, k, component, ri_part;
     int cpies[3];
     int total_dims[3], local_dims[3], start[3] = {0,0,0}, localN;
     real s[3]; /* the step size between grid points dotted with k */
     char name[] = "x.i";
     matrixio_id data_ids[3][2];
     scalar_complex *phasex, *phasey, *phasez;
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

     /* create data sets for writing */
     for (component = 0; component < 3; ++component)
	  if (component == which_component || which_component < 0)
	       for (ri_part = 0; ri_part < 2; ++ri_part) {
		    name[0] = 'x' + component;
		    name[2] = ri_part ? 'i' : 'r';
		    data_ids[component][ri_part] =
			 matrixio_create_dataset(file_id, name, NULL,
						 rank, total_dims);
	       }

     /* cache exp(ikx) along each of the directions, for speed */
     CHK_MALLOC(phasex, scalar_complex, local_nx);
     CHK_MALLOC(phasey, scalar_complex, dims[1]);
     CHK_MALLOC(phasez, scalar_complex, dims[2]);
     for (i = 0; i < local_nx; ++i) {
	  real phase = s[0] * (i + local_x_start);
	  phasex[i].re = cos(phase);
	  phasex[i].im = sin(phase);
     }
     for (j = 0; j < dims[1]; ++j) {
	  real phase = s[1] * j;
	  phasey[j].re = cos(phase);
	  phasey[j].im = sin(phase);
     }
     for (k = 0; k < dims[2]; ++k) {
	  real phase = s[2] * k;
	  phasez[k].re = cos(phase);
	  phasez[k].im = sin(phase);
     }

     /* Now, multiply field by exp(i k*r): */
     for (i = 0; i < local_nx; ++i) {
	  scalar_complex px = phasex[i];

	  for (j = 0; j < dims[1]; ++j) {
	       scalar_complex py;
	       real re = phasey[j].re, im = phasey[j].im;
	       py.re = px.re * re - px.im * im;
	       py.im = px.re * im + px.im * re;
	       
	       for (k = 0; k < dims[2]; ++k) {
		    int ijk = ((i * dims[1] + j) * dims[2] + k) * 3;
		    real p_re, p_im;
		    real re = phasez[k].re, im = phasez[k].im;

		    p_re = py.re * re - py.im * im;
		    p_im = py.re * im + py.im * re;
		    
		    re = field[ijk].re; im = field[ijk].im;
		    field[ijk].re = re * p_re - im * p_im;
		    field[ijk].im = im * p_re + re * p_im;

		    re = field[ijk+1].re; im = field[ijk+1].im;
		    field[ijk+1].re = re * p_re - im * p_im;
		    field[ijk+1].im = im * p_re + re * p_im;

		    re = field[ijk+2].re; im = field[ijk+2].im;
		    field[ijk+2].re = re * p_re - im * p_im;
		    field[ijk+2].im = im * p_re + re * p_im;
	       }
	  }
     }

     free(phasez);
     free(phasey);
     free(phasex);

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
			 if (component == which_component ||
			     which_component < 0)
			      for (ri_part = 0; ri_part < 2; ++ri_part)
				   matrixio_write_real_data(
					data_ids[component][ri_part],
					local_dims, start, 6,
					ri_part ? &field[component].im
					: &field[component].re);
	       }

     /* close data sets */
     for (component = 0; component < 3; ++component)
	  if (component == which_component || which_component < 0)
	       for (ri_part = 0; ri_part < 2; ++ri_part) {
		    matrixio_close_dataset(data_ids[component][ri_part]);
	       }
}

void fieldio_write_real_vals(real *vals,
			     int rank,
			     const int dims[3],
			     int local_nx, int local_x_start,
			     const int copies[3],
			     matrixio_id file_id)
{
     int i, j, k;
     int cpies[3], total_dims[3], local_dims[3], start[3] = {0,0,0};
     matrixio_id data_id;

     for (i = 0; i < 3; ++i) {
	  cpies[i] = MAX2(copies[i], 1); /* make sure copies are non-zero */
	  total_dims[i] = dims[i] * cpies[i];
	  local_dims[i] = dims[i];
     }
     rank = total_dims[2] == 1 ? (total_dims[1] == 1 ? 1 : 2) : 3;
     local_dims[0] = local_nx;

     /* create output file & data set for writing */
     data_id = matrixio_create_dataset(file_id, "data", NULL,
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
}
