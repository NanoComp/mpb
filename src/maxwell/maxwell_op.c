#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <config.h>
#include <check.h>

#include "maxwell.h"

/* assign a = scale * k x v (cross product), going from transverse to
   cartesian coordinates.
  
   Here, a = (a[0],a[1],a[2]) and k = (k.kx,k.ky,k.kz) are in
   cartesian coordinates.  (v[0],v[vstride]) is in the transverse basis of
   k.m and k.n. */
static void assign_cross_t2c(scalar *a, const k_data k,
			     const scalar *v, int vstride,
			     real scale)
{
     scalar v0 = v[0], v1 = v[vstride];

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
   cartesian coordinates.  (v[0],v[vstride]) is in the transverse basis of
   k.m and k.n. */
static void assign_cross_c2t(scalar *v, int vstride,
			     const k_data k, const scalar *a)
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
     ASSIGN_SCALAR(v[vstride],
		   k.kmag * SCALAR_RE(at0),
		   k.kmag * SCALAR_IM(at0));

#ifdef DEBUG
     {
	  real dummy = SCALAR_NORMSQR(v[0]) + SCALAR_NORMSQR(v[vstride]);
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
void maxwell_operator(evectmatrix Xin, evectmatrix Xout, void *data,
		      int is_current_eigenvector)
{
     maxwell_data *d = (maxwell_data *) data;
     int cur_band_start;
     scalar *fft_data;
     scalar_complex *cdata;
     real scale;
     int i, b;
     real *e_field_sums = 0;
     
     CHECK(d, "null maxwell data pointer!");
     CHECK(Xin.c == 2, "fields don't have 2 components!");
     CHECK(d->num_fft_bands == Xin.p, "not yet implemented...");

     fft_data = d->fft_data;

     if (is_current_eigenvector) {
	  e_field_sums = (real *) malloc(sizeof(real) * d->num_fft_bands);
	  CHECK(e_field_sums, "out of memory!");
	  for (b = 0; b < Xin.p; ++b)
	       d->eps_inv_mean[b] = 0.0;
     }

     scale = -1.0 / Xout.N;  /* scale factor to normalize FFT; 
				negative sign comes from 2 i's from curls */

#ifdef SCALAR_COMPLEX

     /********************************************/

     /* first, compute fft_data = curl(Xin): */
     for (i = 0; i < Xin.localN; ++i) {
	  k_data cur_k = d->k_plus_G[i];
	  for (b = 0; b < Xin.p; ++b)
	       assign_cross_t2c(&fft_data[3 * (i * d->num_fft_bands + b)], 
				cur_k, &Xin.data[i * 2 * Xin.p + b], Xin.p,
				scale);
     }

     /********************************************/

     /* Now, multiply fft_data by 1/epsilon using FFT/IFFT pair.
        While we're at it, compute avererage 1/epsilon for each band
	(if is_current_eigenvector is true). */

     compute_fft(+1, d, fft_data, d->num_fft_bands*3, d->num_fft_bands*3, 1);

     if (is_current_eigenvector)
	  for (b = 0; b < d->num_fft_bands; ++b)
	       e_field_sums[b] = 0.0;

     cdata = (scalar_complex *) fft_data;

     if (is_current_eigenvector)
	  for (i = 0; i < d->fft_output_size; ++i) {
	       symmetric_matrix eps_inv = d->eps_inv[i];
	       for (b = 0; b < d->num_fft_bands; ++b) {
		    int ib = 3 * (i * d->num_fft_bands + b);
		    scalar_complex prod0, prod1, prod2;
		    
		    prod0 = cdata[ib];
		    prod1 = cdata[ib + 1];
		    prod2 = cdata[ib + 2];
		    e_field_sums[b] += prod0.re * prod0.re +
			               prod0.im * prod0.im + 
		                       prod1.re * prod1.re +
                                       prod1.im * prod1.im +
		                       prod2.re * prod2.re +
                                       prod2.im * prod2.im;
		    
		    assign_matrix_vector(&cdata[ib], eps_inv, &cdata[ib]);
		    
		    d->eps_inv_mean[b] +=
			 prod0.re*cdata[ib].re   + prod0.im*cdata[ib].im   +
			 prod1.re*cdata[ib+1].re + prod1.im*cdata[ib+1].im +
			 prod2.re*cdata[ib+2].re + prod2.im*cdata[ib+2].im;
	       }
	  }
     else
          for (i = 0; i < d->fft_output_size; ++i) {
               symmetric_matrix eps_inv = d->eps_inv[i];
               for (b = 0; b < d->num_fft_bands; ++b) {
                    int ib = 3 * (i * d->num_fft_bands + b);
		    assign_matrix_vector(&cdata[ib], eps_inv, &cdata[ib]);
	       }
	  }

     if (is_current_eigenvector)
	  for (b = 0; b < d->num_fft_bands; ++b)
	       d->eps_inv_mean[b] /= e_field_sums[b];

     compute_fft(-1, d, fft_data, d->num_fft_bands*3, d->num_fft_bands*3, 1);

     /********************************************/

     /* finally, compute Xout = curl(fft_data): */
     for (i = 0; i < Xout.localN; ++i) {
	  k_data cur_k = d->k_plus_G[i];
	  for (b = 0; b < Xout.p; ++b)
	       assign_cross_c2t(&Xout.data[i * 2 * Xout.p + b], Xout.p,
				cur_k, &fft_data[3 * (i * Xout.p + b)]);
     }

     /********************************************/

#else /* not SCALAR_COMPLEX */

     /* compute Maxwell operator on Xin bands, num_fft_bands at a time: */

     for (cur_band_start = 0; cur_band_start < Xin.p;
	  cur_band_start += d->num_fft_bands) {
	  int num_fft_bands = MIN2(d->num_fft_bands, Xin.p - cur_band_start);

#error Not supported yet!

     }

#endif

     if (is_current_eigenvector)
	  free(e_field_sums);
}

