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
#include <math.h>

#include "../config.h"
#include <check.h>

#include <mpiglue.h>
#include "elastic.h"

/**************************************************************************/

/* function to call z and y parity constraints, if necessary */
void elastic_parity_constraint(evectmatrix X, void *data)
{
     elastic_data *d = (elastic_data *) data;

     CHECK(d, "null elastic data pointer!");
     CHECK(X.c == 3, "fields don't have 3 components!");

     CHECK(d->parity == NO_PARITY, "parity not yet supported in MEB");

/* TODO
     if (d->parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
	  elastic_zparity_constraint(X, data);
     if (d->parity & (EVEN_Y_PARITY | ODD_Y_PARITY))
	  elastic_yparity_constraint(X, data);
*/
}

/**************************************************************************/


/* to fix problems with slow convergence for k ~ 0, manually "put in"
   the k = 0 solution: first three bands have u = constant
   (v = sqrt(rho)) and higher bands are orthogonal.  Note that
   if certain symmetries are imposed, however, some of these constant
   bands may be excluded. */

/* return the number of constant (zero-frequency) bands: */
int elastic_zero_k_num_const_bands(evectmatrix X, elastic_data *d)
{
     int num_const_bands;
     
     CHECK(d, "null elastic data pointer!");
     CHECK(X.c == 3, "fields don't have 3 components!");
     
     num_const_bands = 3;
     if (d->parity & ODD_Z_PARITY)
	  num_const_bands -= 1;
     if (d->parity & ODD_Y_PARITY)
	  num_const_bands -= 1;

     if (num_const_bands > X.p)
	  num_const_bands = X.p;
     
     return num_const_bands;
}

#define MIN2(a,b) ((a) < (b) ? (a) : (b))

void elastic_zero_k_set_const_bands(evectmatrix X, elastic_data *d)
{
     int i, j, b, cur_band_start, num_const_bands;
     scalar *field = (scalar *) d->fft_data;
     scalar_complex *cfield = (scalar_complex *) d->fft_data;
     double scale = 0;
     int comp[3] = {0, 1, 2};
     
     CHECK(d, "null elastic data pointer!");
     CHECK(X.c == 3, "fields don't have 3 components!");

     if (X.p < 1)
	  return;

     num_const_bands = elastic_zero_k_num_const_bands(X, d);
     CHECK(num_const_bands <= 3 && num_const_bands <= X.p, "bug: too many const bands");
     if (d->parity & ODD_Y_PARITY)
	  comp[1] = 2;

     /* compute normalization scale factor: */
     for (i = 0; i < d->fft_output_size; ++i)
	  scale += d->rho[i];
     scale *= X.N;
     scale = 1.0 / sqrt(scale);

     /* set bands to sqrt(rho) in each direction, doing num_fft_bands
	at a time */
     for (cur_band_start = 0; cur_band_start < num_const_bands;
          cur_band_start += d->num_fft_bands) {
          int cur_num_bands = MIN2(d->num_fft_bands, 
				   num_const_bands - cur_band_start);

	  /* set bands to sqrt(rho) */
	  for (i = 0; i < d->fft_output_size; ++i) {
	       real sqrt_rho = 1.0 / d->sqrt_rhoinv[i];
	       for (b = 0; b < cur_num_bands; ++b) {
		    int ic = comp[b + cur_band_start];
		    int ib = 3 * (i * cur_num_bands + b);
		    CASSIGN_SCALAR(cfield[ib + ic], sqrt_rho * scale, 0);
		    CASSIGN_SCALAR(cfield[ib + (ic+1)%3], 0, 0);
		    CASSIGN_SCALAR(cfield[ib + (ic+2)%3], 0, 0);
	       }
	  }
	  
	  /* convert to Fourier space via FFT */
	  elastic_compute_fft(-1, d, field,
			      cur_num_bands*3, cur_num_bands*3, 1);

	  /* assign X to FFT (in field == cfield) */
	  for (i = 0; i < d->other_dims; ++i)
	       for (j = 0; j < d->last_dim; ++j) {
		    int ij = i * d->last_dim + j;
		    int ij2 = i * d->last_dim_size + j;
		    int k;

		    for (b = 0; b < cur_num_bands; ++b)
			 for (k = 0; k < 3; ++k)
			      X.data[(ij * 3 + k) * X.p + b + cur_band_start]
				   = field[3 * (ij2*cur_num_bands + b) + k];
	       }
     }
}

/**************************************************************************/

