#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <config.h>
#include <check.h>

#include "maxwell.h"

/* This file is has too many #ifdef's...blech. */

maxwell_data *create_maxwell_data(int nx, int ny, int nz,
				  int *local_N, int *N_start, int *alloc_N,
				  int num_bands,
				  int num_fft_bands)
{
     int n[3] = {nx, ny, nz}, rank = (nz == 1) ? (ny == 1 ? 1 : 2) : 3;
     maxwell_data *d = 0;
     int fft_data_size;

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
     
     d->num_fft_bands = num_fft_bands;

     d->last_dim = n[rank - 1];

     /* ----------------------------------------------------- */
#ifndef HAVE_MPI 
     d->local_nx = nx; d->local_ny = ny;
     d->local_x_start = d->local_y_start = 0;
     *local_N = *alloc_N = nx * ny * nz;
     *N_start = 0;
     d->other_dims = *local_N / d->last_dim;

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
     d->fft_output_size = fft_data_size = 
	  d->other_dims * 2 * (d->last_dim / 2 + 1);
     d->plan = rfftwnd_create_plan(rank, n, FFTW_FORWARD,
				   FFTW_ESTIMATE | FFTW_IN_PLACE,
				   REAL_TO_COMPLEX);
     d->plan = rfftwnd_create_plan(rank, n, FFTW_BACKWARD,
				   FFTW_ESTIMATE | FFTW_IN_PLACE,
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
     d->k_plus_G_normsqr_inv = (real*) malloc(sizeof(real) * *local_N);
     d->eps_inv_mean = (real*) malloc(sizeof(real) * num_bands);
     CHECK(d->k_plus_G && d->k_plus_G_normsqr_inv && d->eps_inv_mean,
	   "out of memory");

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
	  free(d->k_plus_G_normsqr_inv);
	  free(d->eps_inv_mean);

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

void update_maxwell_data_k(maxwell_data *d, real k[3],
			   real G1[3], real G2[3], real G3[3])
{
     int nx = d->nx, ny = d->ny, nz = d->nz;
     int cx = d->nx/2, cy = d->ny/2, cz = d->nz/2;
     k_data *kpG = d->k_plus_G;
     real *kpGni = d->k_plus_G_normsqr_inv;
     int x, y, z;

     for (x = d->local_x_start; x < d->local_x_start + d->local_nx; ++x) {
	  int kxi = (x > cx) ? (x - nx) : x;
	  for (y = 0; y < ny; ++y) {
	       int kyi = (y > cy) ? (y - ny) : y;
	       for (z = 0; z < nz; ++z, kpG++, kpGni++) {
		    int kzi = (z > cz) ? (z - nz) : z;
		    real a, b, c, leninv;

		    /* Compute k+G: */
		    a = kpG->kx = k[0] + G1[0]*kxi + G2[0]*kyi + G3[0]*kzi;
		    b = kpG->ky = k[1] + G1[1]*kxi + G2[1]*kyi + G3[1]*kzi;
		    c = kpG->kz = k[2] + G1[2]*kxi + G2[2]*kyi + G3[2]*kzi;

		    a = a*a + b*b + c*c;
		    kpG->kmag = sqrt(a);
		    *kpGni = 1.0 / a;
		    
		    /* Now, compute the two normal vectors: */

		    if (kpG->kx == 0.0 && kpG->ky == 0.0) {
			 /* just put n in the x direction if k+G is in z: */
			 kpG->nx = 1.0;
			 kpG->ny = 0.0;
			 kpG->nz = 0.0;
		    }
		    else {
			 /* otherwise, let n = z x (k+G), normalized: */
			 compute_cross(&a, &b, &c,
				       0.0, 0.0, 1.0,
				       kpG->kx, kpG->ky, kpG->kz);
			 leninv = 1.0 / sqrt(a*a + b*b + c*c);
			 kpG->nx = a * leninv;
			 kpG->ny = b * leninv;
			 kpG->nz = c * leninv;
		    }

		    /* m = n x (k+G), normalized */
		    compute_cross(&a, &b, &c,
				  kpG->nx, kpG->ny, kpG->nz,
				  kpG->kx, kpG->ky, kpG->kz);
		    leninv = 1.0 / sqrt(a*a + b*b + c*c);
		    kpG->mx = a * leninv;
		    kpG->my = b * leninv;
		    kpG->mz = c * leninv;
	       }
	  }
     }
}

/* assign a = scale * k x v (cross product), going from transverse to
   cartesian coordinates.
  
   Here, a = (a[0],a[1],a[2]) and k = (k.kx,k.ky,k.kz) are in
   cartesian coordinates.  (v[0],v[1]) is in the transverse basis of
   k.m and k.n. */
static void assign_cross_t2c(scalar *a, const k_data k, const scalar *v,
			     real scale)
{
     scalar v0 = v[0], v1 = v[1];

     /* Note that k x m = |k| n, k x n = - |k| m.  Therefore,
        k x v = k x (v0 m + v1 n) = (v0 n - v1 m) * |k|. */

     scale *= k.kmag;

     ASSIGN_SCALAR(a[0],
		   (SCALAR_RE(v0)*k.nx - SCALAR_RE(v1)*k.mx) * scale,
		   (SCALAR_IM(v0)*k.nx - SCALAR_IM(v1)*k.mx) * scale);
     ASSIGN_SCALAR(a[1],
		   (SCALAR_RE(v0)*k.ny - SCALAR_RE(v1)*k.my) * scale,
		   (SCALAR_IM(v0)*k.ny - SCALAR_IM(v1)*k.my) * scale);
     ASSIGN_SCALAR(a[2],
		   (SCALAR_RE(v0)*k.nz - SCALAR_RE(v1)*k.mz) * scale,
		   (SCALAR_IM(v0)*k.nz - SCALAR_IM(v1)*k.mz) * scale);

#ifdef DEBUG
     scale = SCALAR_NORMSQR(a[0])+SCALAR_NORMSQR(a[1])+SCALAR_NORMSQR(a[2]);
     CHECK(!BADNUM(scale), "yikes, crazy number!");
#endif
}

/* assign v = k x a (cross product), going from cartesian to
   transverse coordinates.
  
   Here, a = (a[0],a[1],a[2]) and k = (k.kx,k.ky,k.kz) are in
   cartesian coordinates.  (v[0],v[1]) is in the transverse basis of
   k.m and k.n. */
static void assign_cross_c2t(scalar *v, const k_data k, const scalar *a)
{
     scalar a0 = a[0], a1 = a[1], a2 = a[2];
     scalar at0, at1;

     /* First, compute at0 = a*m and at1 = a*n.  (Components of a that
	are parallel to k are killed anyway by the cross product.) */

     ASSIGN_SCALAR(at0,
          SCALAR_RE(a0)*k.mx + SCALAR_RE(a1)*k.my + SCALAR_RE(a2)*k.mz,
	  SCALAR_IM(a0)*k.mx + SCALAR_IM(a1)*k.my + SCALAR_IM(a2)*k.mz);
     ASSIGN_SCALAR(at1,
          SCALAR_RE(a0)*k.nx + SCALAR_RE(a1)*k.ny + SCALAR_RE(a2)*k.nz,
	  SCALAR_IM(a0)*k.nx + SCALAR_IM(a1)*k.ny + SCALAR_IM(a2)*k.nz);

     /* Now, k x a = k x (at0*m + at1*n) = (at0*n - at1*m) * |k|. */

     ASSIGN_SCALAR(v[0],
		   - k.kmag * SCALAR_RE(at1),
		   - k.kmag * SCALAR_IM(at1));
     ASSIGN_SCALAR(v[1],
		   k.kmag * SCALAR_RE(at0),
		   k.kmag * SCALAR_IM(at0));

#ifdef DEBUG
     {
	  real dummy = SCALAR_NORMSQR(v[0]) + SCALAR_NORMSQR(v[1]);
	  CHECK(!BADNUM(dummy), "yikes, crazy number!");
     }
#endif
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
static void assign_matrix_vector(scalar_complex *newv,
				 const symmetric_matrix matrix,
				 const scalar_complex *oldv)
{
     scalar_complex v0 = oldv[0], v1 = oldv[1], v2 = oldv[2];

     newv[0].re = matrix.m00 * v0.re + matrix.m01 * v1.re + matrix.m02 * v2.re;
     newv[0].im = matrix.m00 * v0.im + matrix.m01 * v1.im + matrix.m02 * v2.im;

     newv[1].re = matrix.m01 * v0.re + matrix.m11 * v1.re + matrix.m12 * v2.re;
     newv[1].im = matrix.m01 * v0.im + matrix.m11 * v1.im + matrix.m12 * v2.im;

     newv[2].re = matrix.m02 * v0.re + matrix.m12 * v1.re + matrix.m22 * v2.re;
     newv[2].im = matrix.m02 * v0.im + matrix.m12 * v1.im + matrix.m22 * v2.im;

#ifdef DEBUG
     {
	  real dummy;
	  dummy = SCALAR_NORMSQR(newv[0]) + SCALAR_NORMSQR(newv[1])
	       + SCALAR_NORMSQR(newv[2]);
	  CHECK(!BADNUM(dummy), "yikes, crazy number!");
     }
#endif
}

#define MIN2(a,b) ((a) < (b) ? (a) : (b))

/* Compute Xout = curl(1/epsilon * curl(Xin)) */
void maxwell_operator(evectmatrix Xin, evectmatrix Xout, void *data)
{
     maxwell_data *d = (maxwell_data *) data;
     int cur_band_start;
     scalar *fft_data;
     scalar_complex *cdata;
     real scale;
     int i, b;
     
     CHECK(d, "null maxwell data pointer!");
     CHECK(Xin.c == 2, "fields don't have 2 components!");
     CHECK(d->num_fft_bands == Xin.p, "not yet implemented...");

     fft_data = d->fft_data;

     scale = -1.0 / Xout.N;  /* scale factor to normalize FFT; 
				negative sign comes from 2 i's from curls */

#ifdef SCALAR_COMPLEX

     /* first, compute fft_data = curl(Xin): */
     for (i = 0; i < Xin.localN; ++i) {
	  k_data cur_k = d->k_plus_G[i];
	  for (b = 0; b < Xin.p; ++b)
	       assign_cross_t2c(&fft_data[3 * (i * d->num_fft_bands + b)], 
				cur_k, &Xin.data[2 * (i * Xin.p + b)],
				scale);
     }

     /* Now, multiply fft_data by 1/epsilon using FFT/IFFT pair: */

     compute_fft(+1, d, fft_data, d->num_fft_bands*3, d->num_fft_bands*3, 1);

     cdata = (scalar_complex *) fft_data;
     for (i = 0; i < d->fft_output_size; ++i) {
	  symmetric_matrix eps_inv = d->eps_inv[i];
	  for (b = 0; b < d->num_fft_bands; ++b)
	       assign_matrix_vector(&cdata[3 * (i * Xout.p + b)],
				    eps_inv,
				    &cdata[3 * (i * Xout.p + b)]);
     }

     compute_fft(-1, d, fft_data, d->num_fft_bands*3, d->num_fft_bands*3, 1);

     /* finally, compute Xout = curl(fft_data): */
     for (i = 0; i < Xout.localN; ++i) {
	  k_data cur_k = d->k_plus_G[i];
	  for (b = 0; b < Xout.p; ++b)
	       assign_cross_c2t(&Xout.data[2 * (i * Xout.p + b)], 
				cur_k, &fft_data[3 * (i * Xout.p + b)]);
     }

#else /* not SCALAR_COMPLEX */

     /* compute Maxwell operator on Xin bands, num_fft_bands at a time: */

     for (cur_band_start = 0; cur_band_start < Xin.p;
	  cur_band_start += d->num_fft_bands) {
	  int num_fft_bands = MIN2(d->num_fft_bands, Xin.p - cur_band_start);

#error Not supported yet!

     }

#endif
     
}

