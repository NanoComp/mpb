/* Copyright (C) 1999-2014 Massachusetts Institute of Technology.
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

#include "elastic.h"

/* This file is has too many #ifdef's...blech. */

#define MIN2(a,b) ((a) < (b) ? (a) : (b))
#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define TWOPI 6.2831853071795864769252867665590057683943388

elastic_data *create_elastic_data(int nx, int ny, int nz,
				  int *local_N, int *N_start, int *alloc_N,
				  int num_bands,
				  int max_fft_bands)
{
     int n[3], rank = (nz == 1) ? (ny == 1 ? 1 : 2) : 3;
     elastic_data *d = 0;
     int fft_data_size;

     n[0] = nx;
     n[1] = ny;
     n[2] = nz;

#ifndef HAVE_FFTW
#  error Non-FFTW FFTs are not currently supported.
#endif
     
#ifdef HAVE_FFTW
     CHECK(sizeof(fftw_real) == sizeof(real),
	   "floating-point type is inconsistent with FFTW!");
#endif

     CHK_MALLOC(d, elastic_data, 1);

     d->nx = nx;
     d->ny = ny;
     d->nz = nz;
     
     d->max_fft_bands = MIN2(num_bands, max_fft_bands);
     elastic_set_num_bands(d, num_bands);

     d->current_k[0] = d->current_k[1] = d->current_k[2] = 0.0;
     d->zero_k = 1;
     d->parity = NO_PARITY;

     d->last_dim_size = d->last_dim = n[rank - 1];

     /* ----------------------------------------------------- */
#ifndef HAVE_MPI 
     d->local_nx = nx; d->local_ny = ny;
     d->local_x_start = d->local_y_start = 0;
     *local_N = *alloc_N = nx * ny * nz;
     *N_start = 0;
     d->other_dims = *local_N / d->last_dim;

     d->fft_data = 0;  /* initialize it here for use in specific planner? */

#  ifdef HAVE_FFTW
#    ifdef SCALAR_COMPLEX
     d->fft_output_size = fft_data_size = nx * ny * nz;
     d->plan = fftwnd_create_plan_specific(rank, n, FFTW_BACKWARD,
					   FFTW_ESTIMATE | FFTW_IN_PLACE,
					   (fftw_complex*) d->fft_data,
					   3 * d->num_fft_bands,
					   (fftw_complex*) d->fft_data,
					   3 * d->num_fft_bands);
     d->iplan = fftwnd_create_plan_specific(rank, n, FFTW_FORWARD,
					    FFTW_ESTIMATE | FFTW_IN_PLACE,
					    (fftw_complex*) d->fft_data,
					    3 * d->num_fft_bands,
					    (fftw_complex*) d->fft_data,
					    3 * d->num_fft_bands);
#    else /* not SCALAR_COMPLEX */
     d->last_dim_size = 2 * (d->last_dim / 2 + 1);
     d->fft_output_size = (fft_data_size = d->other_dims * d->last_dim_size)/2;
     d->plan = rfftwnd_create_plan_specific(rank, n, FFTW_COMPLEX_TO_REAL,
					    FFTW_ESTIMATE | FFTW_IN_PLACE,
					    (fftw_real*) d->fft_data,
					    3 * d->num_fft_bands,
					    (fftw_real*) d->fft_data,
					    3 * d->num_fft_bands);
     d->iplan = rfftwnd_create_plan_specific(rank, n, FFTW_REAL_TO_COMPLEX,
					     FFTW_ESTIMATE | FFTW_IN_PLACE,
					     (fftw_real*) d->fft_data,
					     3 * d->num_fft_bands,
					     (fftw_real*) d->fft_data,
					     3 * d->num_fft_bands);
#    endif /* not SCALAR_COMPLEX */
#  endif /* HAVE_FFTW */

#else /* HAVE_MPI */
     /* ----------------------------------------------------- */

#  ifdef HAVE_FFTW

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

#    ifdef SCALAR_COMPLEX
     d->iplan = fftwnd_mpi_create_plan(mpb_comm, rank, n,
				       FFTW_FORWARD,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);
     {
	  int nt[3]; /* transposed dimensions for reverse FFT */
	  nt[0] = n[1]; nt[1] = n[0]; nt[2] = n[2]; 
	  d->plan = fftwnd_mpi_create_plan(mpb_comm, rank, nt,
					   FFTW_BACKWARD,
					   FFTW_ESTIMATE | FFTW_IN_PLACE);
     }

     fftwnd_mpi_local_sizes(d->iplan, &d->local_nx, &d->local_x_start,
			    &d->local_ny, &d->local_y_start,
			    &fft_data_size);
     
     d->fft_output_size = nx * d->local_ny * nz;

#    else /* not SCALAR_COMPLEX */

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

     d->iplan = rfftwnd_mpi_create_plan(mpb_comm, rank, n,
					FFTW_REAL_TO_COMPLEX,
					FFTW_ESTIMATE | FFTW_IN_PLACE);

     /* Unlike fftwnd_mpi, we do *not* pass transposed dimensions for
	the reverse transform here--we always pass the dimensions of the
	original real array, and rfftwnd_mpi assumes that if one
	transform is transposed, then the other is as well. */
     d->plan = rfftwnd_mpi_create_plan(mpb_comm, rank, n,
				       FFTW_COMPLEX_TO_REAL,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);

     rfftwnd_mpi_local_sizes(d->iplan, &d->local_nx, &d->local_x_start,
			     &d->local_ny, &d->local_y_start,
			     &fft_data_size);

     d->last_dim_size = 2 * (d->last_dim / 2 + 1);
     if (rank == 2)
	  d->fft_output_size = nx * d->local_ny * nz;
     else
	  d->fft_output_size = nx * d->local_ny * (d->last_dim_size / 2);

#    endif /* not SCALAR_COMPLEX */
     
     *local_N = d->local_nx * ny * nz;
     *N_start = d->local_x_start * ny * nz;
     *alloc_N = *local_N;
     d->other_dims = *local_N / d->last_dim;

#  endif /* HAVE_FFTW */

#endif /* HAVE_MPI */
     /* ----------------------------------------------------- */

#ifdef HAVE_FFTW
     CHECK(d->plan && d->iplan, "FFTW plan creation failed");
#endif

     CHK_MALLOC(d->rho, real, d->fft_output_size);
     CHK_MALLOC(d->sqrt_rhoinv, real, d->fft_output_size);
     CHK_MALLOC(d->rhoct2, real, d->fft_output_size);
     CHK_MALLOC(d->rhocl2, real, d->fft_output_size);

     /* A scratch output array is required because the "ordinary" arrays
	are not in a cartesian basis (or even a constant basis). */
     fft_data_size *= d->max_fft_bands;
     CHK_MALLOC(d->fft_data, scalar, 3 * fft_data_size);

     /* Scratch array(s) for operator.f ...this is way more than
	is necessary, so it will probably change in the future */
     CHK_MALLOC(d->crap, scalar_complex, nx*ny*nz * 21);

     d->local_N = *local_N;
     d->N_start = *N_start;
     d->alloc_N = *alloc_N;
     d->N = nx * ny * nz;

     return d;
}

void destroy_elastic_data(elastic_data *d)
{
     if (d) {

#ifdef HAVE_FFTW
#  ifdef HAVE_MPI
#    ifdef SCALAR_COMPLEX
	  fftwnd_mpi_destroy_plan(d->plan);
	  fftwnd_mpi_destroy_plan(d->iplan);
#    else /* not SCALAR_COMPLEX */
	  rfftwnd_mpi_destroy_plan(d->plan);
	  rfftwnd_mpi_destroy_plan(d->iplan);
#    endif /* not SCALAR_COMPLEX */
#  else /* not HAVE_MPI */
#    ifdef SCALAR_COMPLEX
	  fftwnd_destroy_plan(d->plan);
	  fftwnd_destroy_plan(d->iplan);
#    else /* not SCALAR_COMPLEX */
	  rfftwnd_destroy_plan(d->plan);
	  rfftwnd_destroy_plan(d->iplan);
#    endif /* not SCALAR_COMPLEX */
#  endif /* not HAVE_MPI */
#endif /* HAVE FFTW */

	  free(d->rho);
	  free(d->sqrt_rhoinv);
	  free(d->rhoct2);
	  free(d->rhocl2);
	  free(d->fft_data);
	  free(d->crap);

	  free(d);
     }
}

void elastic_set_num_bands(elastic_data *d, int num_bands)
{
     d->num_bands = num_bands;
     d->num_fft_bands = MIN2(num_bands, d->max_fft_bands);
}

/* Set the current k point for the Elastic solver.  k is given in the
   basis of the reciprocal lattice vectors, G1, G2, and G3. */
void update_elastic_data_k(elastic_data *d, real k[3],
			   real G1[3], real G2[3], real G3[3])
{
     real kx, ky, kz;
     
     d->G[0][0] = G1[0];
     d->G[0][1] = G1[1];
     d->G[0][2] = G1[2];
     d->G[1][0] = G2[0];
     d->G[1][1] = G2[1];
     d->G[1][2] = G2[2];
     d->G[2][0] = G3[0];
     d->G[2][1] = G3[1];
     d->G[2][2] = G3[2];

     kx = G1[0]*k[0] + G2[0]*k[1] + G3[0]*k[2];
     ky = G1[1]*k[0] + G2[1]*k[1] + G3[1]*k[2];
     kz = G1[2]*k[0] + G2[2]*k[1] + G3[2]*k[2];

     d->current_k[0] = kx;
     d->current_k[1] = ky;
     d->current_k[2] = kz;

     d->zero_k = kx == 0.0 && ky == 0.0 && kz == 0.0;

     /* make sure current parity is still valid: */
     set_elastic_data_parity(d, d->parity);
}

void set_elastic_data_parity(elastic_data *d, int parity)
{
     if ((parity & EVEN_Z_PARITY) && (parity & ODD_Z_PARITY))
	  parity &= ~(EVEN_Z_PARITY | ODD_Z_PARITY);
     if (d->current_k[2] != 0.0)
	  parity &= ~(EVEN_Z_PARITY | ODD_Z_PARITY);
     if ((parity & EVEN_Y_PARITY) && (parity & ODD_Y_PARITY))
	  parity &= ~(EVEN_Y_PARITY | ODD_Y_PARITY);
     if (d->current_k[1] != 0.0)
	  parity &= ~(EVEN_Y_PARITY | ODD_Y_PARITY);
     d->parity = parity;
}
