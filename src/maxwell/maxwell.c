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

/* This file is has too many #ifdef's...blech. */

#define MIN2(a,b) ((a) < (b) ? (a) : (b))

maxwell_data *create_maxwell_data(int nx, int ny, int nz,
				  int *local_N, int *N_start, int *alloc_N,
				  int num_bands,
				  int num_fft_bands)
{
     int n[3], rank = (nz == 1) ? (ny == 1 ? 1 : 2) : 3;
     maxwell_data *d = 0;
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

     d = (maxwell_data *) malloc(sizeof(maxwell_data));
     CHECK(d, "out of memory");

     d->nx = nx;
     d->ny = ny;
     d->nz = nz;
     
     d->num_fft_bands = MIN2(num_bands, num_fft_bands);

     d->current_k[0] = d->current_k[1] = d->current_k[2] = 0.0;
     d->polarization = NO_POLARIZATION;

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
     d->plan = fftwnd_create_plan_specific(rank, n, FFTW_FORWARD,
					   FFTW_ESTIMATE | FFTW_IN_PLACE,
					   (FFTW_COMPLEX*) d->fft_data,
					   3 * num_fft_bands,
					   (FFTW_COMPLEX*) d->fft_data,
					   3 * num_fft_bands);
     d->iplan = fftwnd_create_plan_specific(rank, n, FFTW_BACKWARD,
					    FFTW_ESTIMATE | FFTW_IN_PLACE,
					    (FFTW_COMPLEX*) d->fft_data,
					    3 * num_fft_bands,
					    (FFTW_COMPLEX*) d->fft_data,
					    3 * num_fft_bands);
#    else /* not SCALAR_COMPLEX */
     d->last_dim_size = 2 * (d->last_dim / 2 + 1);
     d->fft_output_size = fft_data_size = d->other_dims * d->last_dim_size;
     d->plan = rfftwnd_create_plan(rank, n, FFTW_FORWARD,
				   FFTW_ESTIMATE | FFTW_IN_PLACE,
				   REAL_TO_COMPLEX);
     d->plan = rfftwnd_create_plan(rank, n, FFTW_BACKWARD,
				   FFTW_ESTIMATE | FFTW_IN_PLACE,
				   COMPLEX_TO_REAL);
#    endif /* not SCALAR_COMPLEX */
#  endif /* HAVE_FFTW */

#else /* HAVE_MPI */
     /* ----------------------------------------------------- */

#  ifdef HAVE_FFTW

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

#    ifdef SCALAR_COMPLEX
     d->plan = fftwnd_mpi_create_plan(MPI_COMM_WORLD, rank, n,
				      FFTW_FORWARD,
				      FFTW_ESTIMATE | FFTW_IN_PLACE);
     d->iplan = fftwnd_mpi_create_plan(MPI_COMM_WORLD, rank, n,
				       FFTW_BACKWARD,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);

     fftwnd_mpi_local_sizes(plan, &d->local_nx, &d->local_x_start,
			    &d->local_ny, &d->local_y_start,
			    &fft_data_size);
     
     d->fft_output_size = nx * d->local_ny * nz;

#    else /* not SCALAR_COMPLEX */

     CHECK(rank > 2, "rank <= 2 MPI computations must use SCALAR_COMPLEX");

#  error rfftw MPI transforms not yet supported

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

     d->eps_inv = (symmetric_matrix*) malloc(sizeof(symmetric_matrix)
	                                     * d->fft_output_size);
     CHECK(d->eps_inv, "out of memory");

     /* a scratch output array is required because the "ordinary" arrays
	are not in a cartesian basis (or even a constant basis). */
     d->fft_data = (scalar*) malloc(sizeof(scalar) * 3
				    * num_fft_bands * fft_data_size);
     CHECK(d->fft_data, "out of memory");

     d->k_plus_G = (k_data*) malloc(sizeof(k_data) * *local_N);
     d->k_plus_G_normsqr = (real*) malloc(sizeof(real) * *local_N);
     CHECK(d->k_plus_G && d->k_plus_G_normsqr, "out of memory");

     d->eps_inv_mean = 1.0;

     d->local_N = *local_N;
     d->N_start = *N_start;
     d->alloc_N = *alloc_N;
     d->num_bands = num_bands;
     d->N = nx * ny * nz;

     return d;
}

void destroy_maxwell_data(maxwell_data *d)
{
     if (d) {

#ifdef HAVE_FFTW
#  ifdef HAVE_MPI
	  fftwnd_mpi_destroy_plan(d->plan);
	  fftwnd_mpi_destroy_plan(d->iplan);
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

	  free(d->eps_inv);
	  free(d->fft_data);
	  free(d->k_plus_G);
	  free(d->k_plus_G_normsqr);

	  free(d);
     }
}

/* compute a = b x c */
static void compute_cross(real *a0, real *a1, real *a2,
			  real b0, real b1, real b2,
			  real c0, real c1, real c2)
{
     *a0 = b1 * c2 - b2 * c1;
     *a1 = b2 * c0 - b0 * c2;
     *a2 = b0 * c1 - b1 * c0;
}

/* Set the current k point for the Maxwell solver.  k is given in the
   basis of the reciprocal lattice vectors, G1, G2, and G3. */
void update_maxwell_data_k(maxwell_data *d, real k[3],
			   real G1[3], real G2[3], real G3[3])
{
     int nx = d->nx, ny = d->ny, nz = d->nz;
     int cx = d->nx/2, cy = d->ny/2, cz = d->nz/2;
     k_data *kpG = d->k_plus_G;
     real *kpGn2 = d->k_plus_G_normsqr;
     int x, y, z;
     real kx, ky, kz;

     kx = G1[0]*k[0] + G2[0]*k[1] + G3[0]*k[2];
     ky = G1[1]*k[0] + G2[1]*k[1] + G3[1]*k[2];
     kz = G1[2]*k[0] + G2[2]*k[1] + G3[2]*k[2];

     if (kx == 0.0 && ky == 0.0 && kz == 0.0) {
	  printf("detected zero k\n");
	  kx = 1e-5;
     }

     d->current_k[0] = kx;
     d->current_k[1] = ky;
     d->current_k[2] = kz;

     /* make sure current polarization is still valid: */
     set_maxwell_data_polarization(d, d->polarization);

     for (x = d->local_x_start; x < d->local_x_start + d->local_nx; ++x) {
	  int kxi = (x > cx) ? (x - nx) : x;
	  for (y = 0; y < ny; ++y) {
	       int kyi = (y > cy) ? (y - ny) : y;
	       for (z = 0; z < nz; ++z, kpG++, kpGn2++) {
		    int kzi = (z > cz) ? (z - nz) : z;
		    real kpGx, kpGy, kpGz, a, b, c, leninv;

		    /* Compute k+G: */
		    kpGx = kx + G1[0]*kxi + G2[0]*kyi + G3[0]*kzi;
		    kpGy = ky + G1[1]*kxi + G2[1]*kyi + G3[1]*kzi;
		    kpGz = kz + G1[2]*kxi + G2[2]*kyi + G3[2]*kzi;

		    a = kpGx*kpGx + kpGy*kpGy + kpGz*kpGz;
		    kpG->kmag = sqrt(a);
		    *kpGn2 = a;
		    
		    /* Now, compute the two normal vectors: */

		    if (kpGx == 0.0 && kpGy == 0.0) {
			 /* just put n in the x direction if k+G is in z: */
			 kpG->nx = 1.0;
			 kpG->ny = 0.0;
			 kpG->nz = 0.0;
		    }
		    else {
			 /* otherwise, let n = z x (k+G), normalized: */
			 compute_cross(&a, &b, &c,
				       0.0, 0.0, 1.0,
				       kpGx, kpGy, kpGz);
			 leninv = 1.0 / sqrt(a*a + b*b + c*c);
			 kpG->nx = a * leninv;
			 kpG->ny = b * leninv;
			 kpG->nz = c * leninv;
		    }

		    /* m = n x (k+G), normalized */
		    compute_cross(&a, &b, &c,
				  kpG->nx, kpG->ny, kpG->nz,
				  kpGx, kpGy, kpGz);
		    leninv = 1.0 / sqrt(a*a + b*b + c*c);
		    kpG->mx = a * leninv;
		    kpG->my = b * leninv;
		    kpG->mz = c * leninv;

#ifdef DEBUG
#define DOT(u0,u1,u2,v0,v1,v2) ((u0)*(v0) + (u1)*(v1) + (u2)*(v2))

		    /* check orthogonality */
		    CHECK(fabs(DOT(kpGx, kpGy, kpGz,
				   kpG->nx, kpG->ny, kpG->nz)) < 1e-6,
			  "vectors not orthogonal!");
		    CHECK(fabs(DOT(kpGx, kpGy, kpGz,
				   kpG->mx, kpG->my, kpG->mz)) < 1e-6,
			  "vectors not orthogonal!");
		    CHECK(fabs(DOT(kpG->mx, kpG->my, kpG->mz,
				   kpG->nx, kpG->ny, kpG->nz)) < 1e-6,
			  "vectors not orthogonal!");

		    /* check normalization */
		    CHECK(fabs(DOT(kpG->nx, kpG->ny, kpG->nz,
				   kpG->nx, kpG->ny, kpG->nz) - 1.0) < 1e-6,
			  "vectors not unit vectors!");
		    CHECK(fabs(DOT(kpG->mx, kpG->my, kpG->mz,
				   kpG->mx, kpG->my, kpG->mz) - 1.0) < 1e-6,
			  "vectors not unit vectors!");
#endif
	       }
	  }
     }
}

void set_maxwell_data_polarization(maxwell_data *d,
				   polarization_t polarization)
{
     if (d->current_k[2] != 0.0 || d->nz != 1)
	  polarization = NO_POLARIZATION;
     d->polarization = polarization;
}

maxwell_target_data *create_maxwell_target_data(maxwell_data *md, 
						real target_frequency)
{
     maxwell_target_data *d;

     d = (maxwell_target_data *) malloc(sizeof(maxwell_target_data));
     CHECK(d, "out of memory");

     d->d = md;
     d->target_frequency = target_frequency;

     d->T = create_evectmatrix(md->N, 2, md->num_bands, 
			       md->local_N, md->N_start, md->alloc_N);

     return d;
}

void destroy_maxwell_target_data(maxwell_target_data *d)
{
     if (d) {
	  destroy_evectmatrix(d->T);
	  free(d);
     }
}
