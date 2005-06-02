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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "config.h"

#include <check.h>
#include <matrices.h>

#include "matrixio.h"

#define TWOPI 6.2831853071795864769252867665590057683943388
#define MAX2(a,b) ((a) > (b) ? (a) : (b))

/* note that kvector here is given in the reciprocal basis 
   ...data_id should be of length at 2*num_components */
void fieldio_write_complex_field(scalar_complex *field,
				 int rank,
				 const int dims[3],
				 const int local_dims[3],
				 const int start[3],
				 int which_component, int num_components,
				 const real *kvector,
				 matrixio_id file_id,
				 int append,
				 matrixio_id data_id[])
{
     int i, j, k, component, ri_part;

     rank = dims[2] == 1 ? (dims[1] == 1 ? 1 : 2) : 3;

     if (kvector) {
	  real s[3]; /* the step size between grid points dotted with k */
	  scalar_complex *phasex, *phasey, *phasez;

	  for (i = 0; i < 3; ++i)
	       s[i] = TWOPI * kvector[i] / dims[i];
	  
	  /* cache exp(ikx) along each of the directions, for speed */
	  CHK_MALLOC(phasex, scalar_complex, local_dims[0]);
	  CHK_MALLOC(phasey, scalar_complex, local_dims[1]);
	  CHK_MALLOC(phasez, scalar_complex, local_dims[2]);
	  for (i = 0; i < local_dims[0]; ++i) {
	       real phase = s[0] * (i + start[0]);
	       phasex[i].re = cos(phase);
	       phasex[i].im = sin(phase);
	  }
	  for (j = 0; j < local_dims[1]; ++j) {
	       real phase = s[1] * (j + start[1]);
	       phasey[j].re = cos(phase);
	       phasey[j].im = sin(phase);
	  }
	  for (k = 0; k < local_dims[2]; ++k) {
	       real phase = s[2] * (k + start[2]);
	       phasez[k].re = cos(phase);
	       phasez[k].im = sin(phase);
	  }

	  /* Now, multiply field by exp(i k*r): */
	  for (i = 0; i < local_dims[0]; ++i) {
	       scalar_complex px = phasex[i];
	       
	       for (j = 0; j < local_dims[1]; ++j) {
		    scalar_complex py;
		    real re = phasey[j].re, im = phasey[j].im;
		    py.re = px.re * re - px.im * im;
		    py.im = px.re * im + px.im * re;
		    
		    for (k = 0; k < local_dims[2]; ++k) {
			 int ijk = ((i*local_dims[1] + j)*local_dims[2] + k)*3;
			 real p_re, p_im;
			 real re = phasez[k].re, im = phasez[k].im;
			 
			 p_re = py.re * re - py.im * im;
			 p_im = py.re * im + py.im * re;
			 
			 for (component = 0; component < 3; ++component) {
			      int ijkc = ijk + component;
			      re = field[ijkc].re; im = field[ijkc].im;
			      field[ijkc].re = re * p_re - im * p_im;
			      field[ijkc].im = im * p_re + re * p_im;
			 }
		    }
	       }
	  }
	  
	  free(phasez);
	  free(phasey);
	  free(phasex);
     }

     /* write hyperslabs for each field component: */
     for (component = 0; component < num_components; ++component)
	  if (component == which_component ||
	      which_component < 0)
	       for (ri_part = 0; ri_part < 2; ++ri_part) {
		    char name[] = "x.i";
		    name[0] = (num_components == 1 ? 'c' : 'x') + component;
		    name[2] = ri_part ? 'i' : 'r';

		    if (!append)
			 data_id[component*2 + ri_part] =
			      matrixio_create_dataset(file_id, name, NULL,
						      rank, dims);
		    
		    matrixio_write_real_data(
			 data_id[component*2 + ri_part], local_dims, start, 
			 2 * num_components,
			 ri_part ? &field[component].im
			 : &field[component].re);
	       }
}

void fieldio_write_real_vals(real *vals,
			     int rank,
			     const int dims[3],
			     const int local_dims[3],
			     const int start[3],
			     matrixio_id file_id,
			     int append,
			     const char *dataname,
			     matrixio_id *data_id)
{
     rank = dims[2] == 1 ? (dims[1] == 1 ? 1 : 2) : 3;

     if (!append || *data_id < 0)
	  *data_id = matrixio_create_dataset(file_id, dataname, 
					     NULL, rank,dims);

     matrixio_write_real_data(*data_id,local_dims,start,1,vals);
}
