#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <config.h>
#include <check.h>

#include "maxwell.h"

/* This file is has too many #ifdef's...blech. */

maxwell_data *create_maxwell_data(int nx, int ny, int nz,
				  int *local_N, int *alloc_N,
				  int num_bands,
				  int num_fft_bands)
{
     int n[3] = {nx, ny, nz}, rank = (nz == 1) ? 2 : 3;
     maxwell_data *d = 0;
     int fft_data_size;

#ifndef HAVE_FFTW
#  error Non-FFTW FFTs are not supported.
#endif
     
#ifdef HAVE_FFTW
     CHECK(sizeof(fftw_real) == sizeof(real),
	   "floating-point type is inconsistent with FFTW!");
#endif

     d = malloc(sizeof(maxwell_data));
     CHECK(d, "out of memory");

     d->nx = nx;
     d->ny = ny;
     d->nz = nz;
     
     d->num_fft_bands = num_fft_bands;

     d->last_dim = (rank == 2) ? ny : nz;

     /* ----------------------------------------------------- */
#ifndef HAVE_MPI 
     d->local_nx = nx; d->local_ny = ny;
     d->local_x_start = d->local_y_start = 0;
     *local_N = *alloc_N = nx * ny * nz;
     d->other_dims = *local_N / d->last_dim;

     d->eps_inv = (real*) malloc(sizeof(real) * 6 * *local_N);
     CHECK(d->eps_inv, "out of memory");

#  ifdef HAVE_FFTW
#    ifdef SCALAR_COMPLEX
     fft_data_size = nx * ny * nz;
     d->plan = fftwnd_create_plan_specific(rank, n, FFTW_FORWARD,
					   FFTW_MEASURE | FFTW_IN_PLACE,
					   (FFTW_COMPLEX*) d->fft_data,
					   3 * num_fft_bands,
					   (FFTW_COMPLEX*) d->fft_data,
					   3 * num_fft_bands);
     d->iplan = fftwnd_create_plan_specific(rank, n, FFTW_BACKWARD,
					    FFTW_MEASURE | FFTW_IN_PLACE,
					    (FFTW_COMPLEX*) d->fft_data,
					    3 * num_fft_bands,
					    (FFTW_COMPLEX*) d->fft_data,
					    3 * num_fft_bands);
#    else /* not SCALAR_COMPLEX */
     fft_data_size = d->other_dims * 2 * (d->last_dim / 2 + 1);
     d->plan = rfftwnd_create_plan(rank, n, FFTW_FORWARD,
				   FFTW_MEASURE | FFTW_IN_PLACE,
				   REAL_TO_COMPLEX);
     d->plan = rfftwnd_create_plan(rank, n, FFTW_BACKWARD,
				   FFTW_MEASURE | FFTW_IN_PLACE,
				   COMPLEX_TO_REAL);
#    endif /* not SCALAR_COMPLEX */
#  endif /* HAVE_FFTW */

#else /* HAVE_MPI */
     /*----------------------------------------------------- */

#  ifdef HAVE_FFTW

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

#    ifdef SCALAR_COMPLEX
     d->plan = fftwnd_mpi_create_plan(MPI_COMM_WORLD, rank, n,
				      FFTW_FORWARD,
				      FFTW_MEASURE | FFTW_IN_PLACE);
     d->iplan = fftwnd_mpi_create_plan(MPI_COMM_WORLD, rank, n,
				       FFTW_BACKWARD,
				       FFTW_MEASURE | FFTW_IN_PLACE);

     fftwnd_mpi_local_sizes(plan, &d->local_nx, &d->local_x_start,
			    &d->local_ny, &d->local_y_start,
			    &fft_data_size);
#    else /* not SCALAR_COMPLEX */

     CHECK(rank > 2, "rank <= 2 MPI computations must use SCALAR_COMPLEX");

#  error rfftw MPI transforms not yet supported

#    endif /* not SCALAR_COMPLEX */
     
     *local_N = *alloc_N = = d->local_nx * ny * nz;
     d->other_dims = *local_N / d->last_dim;

     d->eps_inv = (real*) malloc(sizeof(real) * EPS_MATRIX_COMPONENTS
				 * d->local_ny * nx * nz);
     CHECK(d->eps_inv, "out of memory");

#  endif /* HAVE_FFTW */

#endif /* HAVE_MPI */
     /* ----------------------------------------------------- */

#ifdef HAVE_FFTW
     CHECK(d->plan && d->iplan, "FFTW plan creation failed");
#endif

#ifdef SCALAR_COMPLEX
     d->fft_data = NULL; /* transforms are done in place */ 
#else
     /* real-complex transforms need a scratch array for rfftw's padding */
     d->fft_data = (scalar*) malloc(sizeof(scalar) * 3
				    * num_fft_bands * fft_data_size);
#endif

     d->k_plus_G = (real*) malloc(sizeof(real) * 3 * *local_N);
     d->k_plus_G_normsqr_inv = (real*) malloc(sizeof(real) * *local_N);
     d->eps_inv_mean = (real*) malloc(sizeof(real) * num_bands);
     CHECK(d->k_plus_G && d->k_plus_G_normsqr_inv && d->eps_inv_mean,
	   "out of memory");

     return d;
}

void update_maxwell_data_k(maxwell_data *d, real k[3],
			   real G1[3], real G2[3], real G3[3])
{
     int nx = d->nx, ny = d->ny, nz = d->nz;
     int cx = d->nx/2, cy = d->ny/2, cz = d->nz/2;
     real *kpG = d->k_plus_G, *kpGni = d->k_plus_G_normsqr_inv;
     int x, y, z;

     for (x = d->local_x_start; x < d->local_x_start + d->local_nx;
	  ++x, kpG += ny*nz*3, kpGni += ny*nz)
	  for (y = 0; y < ny; ++y, kpG += nz*3, kpGni += nz)
	       for (z = 0; z < nz; ++z, kpG += 3, kpGni += 1) {
		    int kxi, kyi, kzi;
		    
		    kxi = (x > cx) ? (x - nx) : x;
		    kyi = (y > cy) ? (y - ny) : y;
		    kzi = (z > cz) ? (z - nz) : z;

		    kpG[0] = k[0] + G1[0]*kxi + G2[0]*kyi + G3[0]*kzi;
		    kpG[1] = k[1] + G1[1]*kxi + G2[1]*kyi + G3[1]*kzi;
		    kpG[2] = k[2] + G1[2]*kxi + G2[2]*kyi + G3[2]*kzi;

		    *kpGni = 1.0 / (kpG[0]*kpG[0] + kpG[1]*kpG[1] 
				    + kpG[2]*kpG[2]);
	       }
}

#define MIN2(a,b) ((a) < (b) ? (a) : (b))

/* assign a = b x c (cross product). a may be the same as b or c. */
static void assign_cross(scalar *a, const real *b, const scalar *c)
{
     real b0 = b[0], b1 = b[1], b2 = b[2];
     scalar c0 = c[0], c1 = c[1], c2 = c[2];

     ASSIGN_SCALAR(a[0],
		   b1 * SCALAR_RE(c2) - b1 * SCALAR_RE(c1),
		   b1 * SCALAR_IM(c2) - b1 * SCALAR_IM(c1));
     ASSIGN_SCALAR(a[1],
		   b2 * SCALAR_RE(c[3]) - b[3] * SCALAR_RE(c2),
		   b2 * SCALAR_IM(c[3]) - b[3] * SCALAR_IM(c2));
     ASSIGN_SCALAR(a[2],
		   b0 * SCALAR_RE(c1) - b1 * SCALAR_RE(c0),
		   b0 * SCALAR_IM(c1) - b1 * SCALAR_IM(c0));
}

static void compute_fft(int dir, maxwell_data *d, scalar *array, 
			int howmany, int stride, int dist)
{
#ifdef HAVE_FFTW

#  ifdef SCALAR_COMPLEX

#    ifndef HAVE_MPI

     fftwnd(dir > 0 ? d->plan : d->iplan,
	    howmany,
	    (fftw_complex *) array, stride, dist,
	    0, 0, 0);

#    else /* HAVE_MPI */

     CHECK(dist == howmany && stride == 1,
	   "weird strides and dists don't work with fftwnd_mpi");

     fftwnd_mpi(dir > 0 ? d->plan : d->iplan,
		howmany,
		(fftw_complex *) array);

#    endif /* HAVE_MPI */

#  else /* not SCALAR_COMPLEX */

#    ifndef HAVE_MPI

     if (dir > 0)
	  rfftwnd_real_to_complex(d->plan,
				  howmany,
				  (fftw_real *) array, stride, dist,
				  0, 0, 0);
     else
	  rfftwnd_complex_to_real(d->iplan,
				  howmany,
				  (fftw_complex *) array, stride, dist,
				  0, 0, 0);

#    else /* HAVE_MPI */

#      error rfftwnd_mpi not supported yet.

#    endif /* HAVE_MPI */

#  endif /* not SCALAR_COMPLEX */

#else /* not HAVE_FFTW */
#  error only FFTW ffts are supported
#endif /* not HAVE_FFTW */
}

/* assigns newv = matrix * oldv.  matrix is symmetric and so is stored
   in "packed" format. */
static void assign_matrix_vector(scalar_complex *newv, const real *matrix,
				 const scalar_complex *oldv)
{
     real m00 = matrix[0], m01 = matrix[1], m02 = matrix[2],
	                   m11 = matrix[3], m12 = matrix[4],
	                                    m22 = matrix[5];
     scalar_complex v0 = oldv[0], v1 = oldv[1], v2 = oldv[2];

     newv[0].re = m00 * v0.re + m01 * v1.re + m02 * v2.re;
     newv[0].im = m00 * v0.im + m01 * v1.im + m02 * v2.im;

     newv[1].re = m01 * v0.re + m11 * v1.re + m12 * v2.re;
     newv[1].im = m01 * v0.im + m11 * v1.im + m12 * v2.im;

     newv[2].re = m02 * v0.re + m12 * v1.re + m22 * v2.re;     
     newv[2].im = m02 * v0.im + m12 * v1.im + m22 * v2.im;     
}

/* Compute Xout = curl(1/epsilon * curl(Xin)) */
void maxwell_operator(evectmatrix Xin, evectmatrix Xout, void *data)
{
     maxwell_data *d = (maxwell_data *) data;
     int cur_band_start;
     scalar *fft_data;
     scalar_complex *cdata;
     real *kpG;
     int i, b;
     
     CHECK(d, "null maxwell data pointer!");
     CHECK(Xin.c == 3, "fields don't have 3 components!");

     fft_data = d->fft_data;
     kpG = d->k_plus_G;

#ifdef SCALAR_COMPLEX

     /* first, compute Xout = curl(Xin): */
     for (i = 0; i < Xin.localN; ++i)
	  for (b = 0; b < 3 * Xin.p; ++b)
	       assign_cross(&Xout.data[i * 3 * Xin.p + b], 
			    &kpG[i], &Xin.data[i * 3 * Xin.p + b]);

     /* Now, multiply Xout by 1/epsilon using FFT/IFFT pair: */

     compute_fft(+1, d, Xout.data, Xout.p * Xout.c, Xout.p * Xout.c, 1);

     cdata = (scalar_complex *) Xout.data;
     for (i = 0; i < Xout.localN; ++i)
	  for (b = 0; b < Xout.p; ++b)
	       assign_matrix_vector(&cdata[3 * (i * Xout.p + b)],
				    &d->eps_inv[EPS_MATRIX_COMPONENTS * i],
				    &cdata[3 * (i * Xout.p + b)]);

     compute_fft(-1, d, Xout.data, Xout.p * Xout.c, Xout.p * Xout.c, 1);

     /* finally, compute Xout = curl(Xout): */
     for (i = 0; i < Xout.localN; ++i)
	  for (b = 0; b < 3 * Xout.p; ++b)
	       assign_cross(&Xout.data[i * 3 * Xout.p + b], 
			    &kpG[i], &Xout.data[i * 3 * Xout.p + b]);

#else /* not SCALAR_COMPLEX */

     /* compute Maxwell operator on Xin bands, num_fft_bands at a time: */

     for (cur_band_start = 0; cur_band_start < Xin.p;
	  cur_band_start += d->num_fft_bands) {
	  int num_fft_bands = MIN2(d->num_fft_bands, Xin.p - cur_band_start);

#error Not supported yet!

     }

#endif
     
}

