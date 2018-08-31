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

#include "imaxwell.h"
#include "check.h"

/* This file is has too many #ifdef's...blech. */

#define MIN2(a,b) ((a) < (b) ? (a) : (b))
#define MAX2(a,b) ((a) > (b) ? (a) : (b))

maxwell_data *create_maxwell_data(int nx, int ny, int nz,
				  int *local_N, int *N_start, int *alloc_N,
				  int num_bands,
				  int max_fft_bands)
{
     int n[3], rank = (nz == 1) ? (ny == 1 ? 1 : 2) : 3;
     maxwell_data *d = 0;
     int fft_data_size;

     n[0] = nx;
     n[1] = ny;
     n[2] = nz;

#if !defined(HAVE_FFTW) && !defined(HAVE_FFTW3)
#  error Non-FFTW FFTs are not currently supported.
#endif
     

#if defined(HAVE_FFTW)
     CHECK(sizeof(fftw_real) == sizeof(real),
	   "floating-point type is inconsistent with FFTW!");
#endif

     CHK_MALLOC(d, maxwell_data, 1);

     d->nx = nx;
     d->ny = ny;
     d->nz = nz;
     
     d->max_fft_bands = MIN2(num_bands, max_fft_bands);
     maxwell_set_num_bands(d, num_bands);

     d->current_k[0] = d->current_k[1] = d->current_k[2] = 0.0;
     d->parity = NO_PARITY;

     d->last_dim_size = d->last_dim = n[rank - 1];

     /* ----------------------------------------------------- */
     d->nplans = 1;
#ifndef HAVE_MPI 
     d->local_nx = nx; d->local_ny = ny;
     d->local_x_start = d->local_y_start = 0;
     *local_N = *alloc_N = nx * ny * nz;
     *N_start = 0;
     d->other_dims = *local_N / d->last_dim;

     d->fft_data = 0;  /* initialize it here for use in specific planner? */

#  if defined(HAVE_FFTW3)
     d->nplans = 0; /* plans will be created as needed */
#    ifdef SCALAR_COMPLEX
     d->fft_output_size = fft_data_size = nx * ny * nz;
#    else
     d->last_dim_size = 2 * (d->last_dim / 2 + 1);
     d->fft_output_size = (fft_data_size = d->other_dims * d->last_dim_size)/2;
#    endif

#  elif defined(HAVE_FFTW)
#    ifdef SCALAR_COMPLEX
     d->fft_output_size = fft_data_size = nx * ny * nz;
     d->plans[0] = fftwnd_create_plan_specific(rank, n, FFTW_BACKWARD,
					   FFTW_ESTIMATE | FFTW_IN_PLACE,
					   (fftw_complex*) d->fft_data,
					   3 * d->num_fft_bands,
					   (fftw_complex*) d->fft_data,
					   3 * d->num_fft_bands);
     d->iplans[0] = fftwnd_create_plan_specific(rank, n, FFTW_FORWARD,
					    FFTW_ESTIMATE | FFTW_IN_PLACE,
					    (fftw_complex*) d->fft_data,
					    3 * d->num_fft_bands,
					    (fftw_complex*) d->fft_data,
					    3 * d->num_fft_bands);
#    else /* not SCALAR_COMPLEX */
     d->last_dim_size = 2 * (d->last_dim / 2 + 1);
     d->fft_output_size = (fft_data_size = d->other_dims * d->last_dim_size)/2;
     d->plans[0] = rfftwnd_create_plan_specific(rank, n, FFTW_COMPLEX_TO_REAL,
					    FFTW_ESTIMATE | FFTW_IN_PLACE,
					    (fftw_real*) d->fft_data,
					    3 * d->num_fft_bands,
					    (fftw_real*) d->fft_data,
					    3 * d->num_fft_bands);
     d->iplans[0] = rfftwnd_create_plan_specific(rank, n, FFTW_REAL_TO_COMPLEX,
					     FFTW_ESTIMATE | FFTW_IN_PLACE,
					     (fftw_real*) d->fft_data,
					     3 * d->num_fft_bands,
					     (fftw_real*) d->fft_data,
					     3 * d->num_fft_bands);
#    endif /* not SCALAR_COMPLEX */
#  endif /* HAVE_FFTW */

#else /* HAVE_MPI */
     /* ----------------------------------------------------- */

#  if defined(HAVE_FFTW3)
{
     int i;
     ptrdiff_t np[3], local_nx, local_ny, local_x_start, local_y_start;

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

     d->nplans = 0; /* plans will be created as needed */

     for (i = 0; i < rank; ++i) np[i] = n[i];
     
#    ifndef SCALAR_COMPLEX
     d->last_dim_size = 2 * (np[rank-1] = d->last_dim / 2 + 1);
#    endif

     fft_data_size = *alloc_N 
	  = FFTW(mpi_local_size_transposed)(rank, np, mpb_comm,
					    &local_nx, &local_x_start,
					    &local_ny, &local_y_start);
#    ifndef SCALAR_COMPLEX
     fft_data_size = (*alloc_N *= 2); // convert to # of real scalars
#    endif

     d->local_nx = local_nx;
     d->local_x_start = local_x_start;
     d->local_ny = local_ny;
     d->local_y_start = local_y_start;

     d->fft_output_size = nx * d->local_ny * (rank==3 ? np[2] : nz);
     *local_N = d->local_nx * ny * nz;
     *N_start = d->local_x_start * ny * nz;
     d->other_dims = *local_N / d->last_dim;
}
#  elif defined(HAVE_FFTW)

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

#    ifdef SCALAR_COMPLEX
     d->iplans[0] = fftwnd_mpi_create_plan(mpb_comm, rank, n,
				       FFTW_FORWARD,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);
     {
	  int nt[3]; /* transposed dimensions for reverse FFT */
	  nt[0] = n[1]; nt[1] = n[0]; nt[2] = n[2]; 
	  d->plans[0] = fftwnd_mpi_create_plan(mpb_comm, rank, nt,
					   FFTW_BACKWARD,
					   FFTW_ESTIMATE | FFTW_IN_PLACE);
     }

     fftwnd_mpi_local_sizes(d->iplans[0], &d->local_nx, &d->local_x_start,
			    &d->local_ny, &d->local_y_start,
			    &fft_data_size);
     
     d->fft_output_size = nx * d->local_ny * nz;

#    else /* not SCALAR_COMPLEX */

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

     d->iplans[0] = rfftwnd_mpi_create_plan(mpb_comm, rank, n,
					FFTW_REAL_TO_COMPLEX,
					FFTW_ESTIMATE | FFTW_IN_PLACE);

     /* Unlike fftwnd_mpi, we do *not* pass transposed dimensions for
	the reverse transform here--we always pass the dimensions of the
	original real array, and rfftwnd_mpi assumes that if one
	transform is transposed, then the other is as well. */
     d->plans[0] = rfftwnd_mpi_create_plan(mpb_comm, rank, n,
				       FFTW_COMPLEX_TO_REAL,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);

     rfftwnd_mpi_local_sizes(d->iplans[0], &d->local_nx, &d->local_x_start,
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
     CHECK(d->plans[0] && d->iplans[0], "FFTW plan creation failed");
#endif

     CHK_MALLOC(d->eps_inv, symmetric_matrix, d->fft_output_size);
     d->mu_inv = NULL;

     /* A scratch output array is required because the "ordinary" arrays
	are not in a cartesian basis (or even a constant basis). */
     fft_data_size *= d->max_fft_bands;
#if defined(HAVE_FFTW3)
     d->fft_data = (scalar *) FFTW(malloc)(sizeof(scalar) * 3 * fft_data_size);
     CHECK(d->fft_data, "out of memory!");
     d->fft_data2 = d->fft_data; /* works in-place */
#else     
     CHK_MALLOC(d->fft_data, scalar, 3 * fft_data_size);
     d->fft_data2 = d->fft_data; /* works in-place */
#endif

     CHK_MALLOC(d->k_plus_G, k_data, *local_N);
     CHK_MALLOC(d->k_plus_G_normsqr, real, *local_N);

     d->eps_inv_mean = 1.0;
     d->mu_inv_mean = 1.0;

     d->local_N = *local_N;
     d->N_start = *N_start;
     d->alloc_N = *alloc_N;
     d->N = nx * ny * nz;

     return d;
}

void destroy_maxwell_data(maxwell_data *d)
{
     if (d) {
	  int i;

	  for (i = 0; i < d->nplans; ++i) {
#if defined(HAVE_FFTW3)
	       FFTW(destroy_plan)((fftplan) (d->plans[i]));
	       FFTW(destroy_plan)((fftplan) (d->iplans[i]));
#elif defined(HAVE_FFTW)
#  ifdef HAVE_MPI
#    ifdef SCALAR_COMPLEX
	       fftwnd_mpi_destroy_plan((fftplan) (d->plans[i]));
	       fftwnd_mpi_destroy_plan((fftplan) (d->iplans[i]));
#    else /* not SCALAR_COMPLEX */
	       rfftwnd_mpi_destroy_plan((fftplan) (d->plans[i]));
	       rfftwnd_mpi_destroy_plan((fftplan) (d->iplans[i]));
#    endif /* not SCALAR_COMPLEX */
#  else /* not HAVE_MPI */
#    ifdef SCALAR_COMPLEX
	       fftwnd_destroy_plan((fftplan) (d->plans[i]));
	       fftwnd_destroy_plan((fftplan) (d->iplans[i]));
#    else /* not SCALAR_COMPLEX */
	       rfftwnd_destroy_plan((fftplan) (d->plans[i]));
	       rfftwnd_destroy_plan((fftplan) (d->iplans[i]));
#    endif /* not SCALAR_COMPLEX */
#  endif /* not HAVE_MPI */
#endif /* HAVE FFTW */
	  }

	  free(d->eps_inv);
          if (d->mu_inv) free(d->mu_inv);
#if defined(HAVE_FFTW3)
	  FFTW(free)(d->fft_data);
	  if (d->fft_data2 != d->fft_data)
	       FFTW(free)(d->fft_data2);
#else
	  free(d->fft_data);
#endif
	  free(d->k_plus_G);
	  free(d->k_plus_G_normsqr);

	  free(d);
     }
}

void maxwell_set_num_bands(maxwell_data *d, int num_bands)
{
     d->num_bands = num_bands;
     d->num_fft_bands = MIN2(num_bands, d->max_fft_bands);
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

typedef struct {
    double amp;
    int rank;
} amp_and_rank;

void maxwell_dominant_planewave(maxwell_data *d, evectmatrix H, int band, double kdom[3])
{
    double max_amp = 0;
    int max_i = 0, i = 0;

    CHECK(d != NULL, "maxwell_data is NULL");
    CHECK(1 <= band && band <= H.p, "band out of range");

    /* find biggest planewave component in band 'band' of H on current process: */
    for (i = 0; i < H.localN; ++i) {
        double amp = SCALAR_NORMSQR(H.data[(i*2+0)*H.p + band-1]) + SCALAR_NORMSQR(H.data[(i*2+1)*H.p + band-1]);
        if (amp > max_amp) {
            max_amp = amp;
            max_i = i;
        }
    }

    k_data cur_k = d->k_plus_G[max_i];

    /* set kdom to cur_k.kmag * (cross product of cur_k.m and cur_k.n) */
    compute_cross(&kdom[0], &kdom[1], &kdom[2],
                  cur_k.mx, cur_k.my, cur_k.mz,
                  cur_k.nx, cur_k.ny, cur_k.nz);

    kdom[0] *= cur_k.kmag;
    kdom[1] *= cur_k.kmag;
    kdom[2] *= cur_k.kmag;

#ifdef HAVE_MPI
    /* find kdom of biggest max_amp across all processes */
    amp_and_rank local_res = {max_amp, my_global_rank()};
    amp_and_rank global_res = {0, 0};

    MPI_Allreduce(&local_res, &global_res, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpb_comm);
    MPI_Bcast(kdom, 3, MPI_DOUBLE, global_res.rank, mpb_comm);
#endif
}

/* Set the current k point for the Maxwell solver.  k is given in the
   basis of the reciprocal lattice vectors, G1, G2, and G3. */
void update_maxwell_data_k(maxwell_data *d, real k[3],
			   real G1[3], real G2[3], real G3[3])
{
     int nx = d->nx, ny = d->ny, nz = d->nz;
     int cx = MAX2(1,d->nx/2), cy = MAX2(1,d->ny/2), cz = MAX2(1,d->nz/2);
     k_data *kpG = d->k_plus_G;
     real *kpGn2 = d->k_plus_G_normsqr;
     int x, y, z;
     real kx, ky, kz;

     kx = G1[0]*k[0] + G2[0]*k[1] + G3[0]*k[2];
     ky = G1[1]*k[0] + G2[1]*k[1] + G3[1]*k[2];
     kz = G1[2]*k[0] + G2[2]*k[1] + G3[2]*k[2];

     d->zero_k = kx == 0.0 && ky == 0.0 && kz == 0.0;

     d->current_k[0] = kx;
     d->current_k[1] = ky;
     d->current_k[2] = kz;

     /* make sure current parity is still valid: */
     set_maxwell_data_parity(d, d->parity);

     for (x = d->local_x_start; x < d->local_x_start + d->local_nx; ++x) {
	  int kxi = (x >= cx) ? (x - nx) : x;
	  for (y = 0; y < ny; ++y) {
	       int kyi = (y >= cy) ? (y - ny) : y;
	       for (z = 0; z < nz; ++z, kpG++, kpGn2++) {
		    int kzi = (z >= cz) ? (z - nz) : z;
		    real kpGx, kpGy, kpGz, a, b, c, leninv;

		    /* Compute k+G (noting that G is negative because
		       of the choice of sign in the FFTW Fourier transform): */
		    kpGx = kx - (G1[0]*kxi + G2[0]*kyi + G3[0]*kzi);
		    kpGy = ky - (G1[1]*kxi + G2[1]*kyi + G3[1]*kzi);
		    kpGz = kz - (G1[2]*kxi + G2[2]*kyi + G3[2]*kzi);

		    a = kpGx*kpGx + kpGy*kpGy + kpGz*kpGz;
		    kpG->kmag = sqrt(a);
		    *kpGn2 = a;
		    
		    /* Now, compute the two normal vectors: */
		    /* (Note that we choose them so that m has odd/even
		       parity in z/y, and n is even/odd in z/y.) */

		    if (a == 0) {
			 kpG->nx = 0.0; kpG->ny = 1.0; kpG->nz = 0.0;
			 kpG->mx = 0.0; kpG->my = 0.0; kpG->mz = 1.0;
		    }
		    else {
			 if (kpGx == 0.0 && kpGy == 0.0) {
			      /* put n in the y direction if k+G is in z: */
			      kpG->nx = 0.0;
			      kpG->ny = 1.0;
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
		    }

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

void set_maxwell_data_parity(maxwell_data *d, int parity)
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

maxwell_target_data *create_maxwell_target_data(maxwell_data *md, 
						real target_frequency)
{
     maxwell_target_data *d;

     CHK_MALLOC(d, maxwell_target_data, 1);

     d->d = md;
     d->target_frequency = target_frequency;

     return d;
}

void destroy_maxwell_target_data(maxwell_target_data *d)
{
     if (d) {
	  free(d);
     }
}
