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
#include <math.h>

#include "../config.h"
#include <check.h>

#include "maxwell.h"

/**************************************************************************/

/* assign a = v going from transverse to cartesian coordinates.  
   Here, a = (a[0],a[1],a[2]) is in cartesian coordinates.
   (v[0],v[vstride]) is in the transverse basis of k.m and k.n. */
static void assign_t2c(scalar *a, const k_data k,
		       const scalar *v, int vstride)
{
     scalar v0 = v[0], v1 = v[vstride];

     ASSIGN_SCALAR(a[0],
		   SCALAR_RE(v0)*k.mx + SCALAR_RE(v1)*k.nx,
		   SCALAR_IM(v0)*k.mx + SCALAR_IM(v1)*k.nx);
     ASSIGN_SCALAR(a[1],
		   SCALAR_RE(v0)*k.my + SCALAR_RE(v1)*k.ny,
		   SCALAR_IM(v0)*k.my + SCALAR_IM(v1)*k.ny);
     ASSIGN_SCALAR(a[2],
		   SCALAR_RE(v0)*k.mz + SCALAR_RE(v1)*k.nz,
		   SCALAR_IM(v0)*k.mz + SCALAR_IM(v1)*k.nz);
}

/* assign a = k x v (cross product), going from transverse to
   cartesian coordinates.
  
   Here, a = (a[0],a[1],a[2]) and k = (k.kx,k.ky,k.kz) are in
   cartesian coordinates.  (v[0],v[vstride]) is in the transverse basis of
   k.m and k.n. */
static void assign_cross_t2c(scalar *a, const k_data k,
			     const scalar *v, int vstride)
{
     scalar v0 = v[0], v1 = v[vstride];

     /* Note that k x m = |k| n, k x n = - |k| m.  Therefore,
        k x v = k x (v0 m + v1 n) = (v0 n - v1 m) * |k|. */

     ASSIGN_SCALAR(a[0],
		   (SCALAR_RE(v0)*k.nx - SCALAR_RE(v1)*k.mx) * k.kmag,
		   (SCALAR_IM(v0)*k.nx - SCALAR_IM(v1)*k.mx) * k.kmag);
     ASSIGN_SCALAR(a[1],
		   (SCALAR_RE(v0)*k.ny - SCALAR_RE(v1)*k.my) * k.kmag,
		   (SCALAR_IM(v0)*k.ny - SCALAR_IM(v1)*k.my) * k.kmag);
     ASSIGN_SCALAR(a[2],
		   (SCALAR_RE(v0)*k.nz - SCALAR_RE(v1)*k.mz) * k.kmag,
		   (SCALAR_IM(v0)*k.nz - SCALAR_IM(v1)*k.mz) * k.kmag);

#ifdef DEBUG
     {
	  real num;
	  num = SCALAR_NORMSQR(a[0])+SCALAR_NORMSQR(a[1])+SCALAR_NORMSQR(a[2]);
	  CHECK(!BADNUM(num), "yikes, crazy number!");
     }
#endif
}

/* assign v = scale * k x a (cross product), going from cartesian to
   transverse coordinates.
  
   Here, a = (a[0],a[1],a[2]) and k = (k.kx,k.ky,k.kz) are in
   cartesian coordinates.  (v[0],v[vstride]) is in the transverse basis of
   k.m and k.n. */
static void assign_cross_c2t(scalar *v, int vstride,
			     const k_data k, const scalar *a,
			     real scale)
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

     scale *= k.kmag;  /* combine scale factor and |k|*/
     ASSIGN_SCALAR(v[0],
		   - scale * SCALAR_RE(at1),
		   - scale * SCALAR_IM(at1));
     ASSIGN_SCALAR(v[vstride],
		   scale * SCALAR_RE(at0),
		   scale * SCALAR_IM(at0));

#ifdef DEBUG
     {
	  real dummy = SCALAR_NORMSQR(v[0]) + SCALAR_NORMSQR(v[vstride]);
	  CHECK(!BADNUM(dummy), "yikes, crazy number!");
     }
#endif
}

/* compute a = u x v, where a and u are in cartesian coordinates and
   v is in transverse coordinates. */
static void assign_ucross_t2c(scalar *a, const real u[3], const k_data k,
			     const scalar *v, int vstride)
{
     scalar v0 = v[0], v1 = v[vstride];
     real vx_r, vy_r, vz_r;
#ifdef SCALAR_COMPLEX
     real vx_i, vy_i, vz_i;
#endif

     /* Note that v = (vx,vy,vz) = (v0 m + v1 n). */

     vx_r = SCALAR_RE(v0)*k.mx + SCALAR_RE(v1)*k.nx;
     vy_r = SCALAR_RE(v0)*k.my + SCALAR_RE(v1)*k.ny;
     vz_r = SCALAR_RE(v0)*k.mz + SCALAR_RE(v1)*k.nz;

#ifdef SCALAR_COMPLEX
     vx_i = SCALAR_IM(v0)*k.mx + SCALAR_IM(v1)*k.nx;
     vy_i = SCALAR_IM(v0)*k.my + SCALAR_IM(v1)*k.ny;
     vz_i = SCALAR_IM(v0)*k.mz + SCALAR_IM(v1)*k.nz;
#endif

     ASSIGN_SCALAR(a[0],
		   u[1] * vz_r - u[2] * vy_r,
		   u[1] * vz_i - u[2] * vy_i);
     ASSIGN_SCALAR(a[1],
		   u[2] * vx_r - u[0] * vz_r,
		   u[2] * vx_i - u[0] * vz_i);
     ASSIGN_SCALAR(a[2],
		   u[0] * vy_r - u[1] * vx_r,
		   u[0] * vy_i - u[1] * vx_i);
}

/**************************************************************************/

void maxwell_compute_fft(int dir, maxwell_data *d, scalar *array, 
			 int howmany, int stride, int dist)
{
#ifdef HAVE_FFTW

#  ifdef SCALAR_COMPLEX

#    ifndef HAVE_MPI

     fftwnd(dir < 0 ? d->plan : d->iplan,
	    howmany,
	    (fftw_complex *) array, stride, dist,
	    0, 0, 0);

#    else /* HAVE_MPI */

     CHECK(stride == howmany && dist == 1,
	   "weird strides and dists don't work with fftwnd_mpi");

     fftwnd_mpi(dir < 0 ? d->plan : d->iplan,
		howmany,
		(fftw_complex *) array, (fftw_complex *) NULL,
		FFTW_TRANSPOSED_ORDER);

#    endif /* HAVE_MPI */

#  else /* not SCALAR_COMPLEX */

#    ifndef HAVE_MPI

     if (dir > 0)
	  rfftwnd_real_to_complex(d->iplan,
				  howmany,
				  (fftw_real *) array, stride, dist,
				  0, 0, 0);
     else
	  rfftwnd_complex_to_real(d->plan,
				  howmany,
				  (fftw_complex *) array, stride, dist,
				  0, 0, 0);

#    else /* HAVE_MPI */

     CHECK(stride == howmany && dist == 1,
	   "weird strides and dists don't work with rfftwnd_mpi");
     
     rfftwnd_mpi(dir < 0 ? d->plan : d->iplan,
		 howmany, array, (scalar *) NULL,
		 FFTW_TRANSPOSED_ORDER);

#    endif /* HAVE_MPI */

#  endif /* not SCALAR_COMPLEX */

#else /* not HAVE_FFTW */
#  error only FFTW ffts are supported
#endif /* not HAVE_FFTW */
}

/**************************************************************************/

/* assigns newv = matrix * oldv.  matrix is symmetric and so is stored
   in "packed" format. */
void assign_symmatrix_vector(scalar_complex *newv,
			     const symmetric_matrix matrix,
			     const scalar_complex *oldv)
{
     scalar_complex v0 = oldv[0], v1 = oldv[1], v2 = oldv[2];

#if defined(WITH_HERMITIAN_EPSILON)
     newv[0].re = matrix.m00 * v0.re;
     newv[0].im = matrix.m00 * v0.im;
     CACCUMULATE_SUM_MULT(newv[0], matrix.m01, v1);
     CACCUMULATE_SUM_MULT(newv[0], matrix.m02, v2);

     newv[1].re = matrix.m11 * v1.re;
     newv[1].im = matrix.m11 * v1.im;
     CACCUMULATE_SUM_CONJ_MULT(newv[1], matrix.m01, v0);
     CACCUMULATE_SUM_MULT(newv[1], matrix.m12, v2);

     newv[2].re = matrix.m22 * v2.re;
     newv[2].im = matrix.m22 * v2.im;
     CACCUMULATE_SUM_CONJ_MULT(newv[2], matrix.m02, v0);
     CACCUMULATE_SUM_CONJ_MULT(newv[2], matrix.m12, v1);
#else
     newv[0].re = matrix.m00 * v0.re + matrix.m01 * v1.re + matrix.m02 * v2.re;
     newv[0].im = matrix.m00 * v0.im + matrix.m01 * v1.im + matrix.m02 * v2.im;

     newv[1].re = matrix.m01 * v0.re + matrix.m11 * v1.re + matrix.m12 * v2.re;
     newv[1].im = matrix.m01 * v0.im + matrix.m11 * v1.im + matrix.m12 * v2.im;

     newv[2].re = matrix.m02 * v0.re + matrix.m12 * v1.re + matrix.m22 * v2.re;
     newv[2].im = matrix.m02 * v0.im + matrix.m12 * v1.im + matrix.m22 * v2.im;
#endif

#ifdef DEBUG
     {
	  real dummy;
	  dummy = CSCALAR_NORMSQR(newv[0]) + CSCALAR_NORMSQR(newv[1])
	       + CSCALAR_NORMSQR(newv[2]);
	  CHECK(!BADNUM(dummy), "yikes, crazy number!");
     }
#endif
}

/* compute the D field in position space from Hin, which holds the H
   field in Fourier space, for the specified bands; this amounts to
   taking the curl and then Fourier transforming.  The output array,
   dfield, is fft_output_size x cur_num_bands x 3, where
   fft_output_size is the local spatial indices and 3 is the field
   components. 

   Note: actually, this computes just (k+G) x H, whereas the actual D
   field is i/omega i(k+G) x H...so, we are really computing -omega*D,
   here. */
void maxwell_compute_d_from_H(maxwell_data *d, evectmatrix Hin, 
			      scalar_complex *dfield,
			      int cur_band_start, int cur_num_bands)
{
     scalar *fft_data = (scalar *) dfield;
     int i, j, b;

     CHECK(Hin.c == 2, "fields don't have 2 components!");
     CHECK(d, "null maxwell data pointer!");
     CHECK(dfield, "null field output data!");
     CHECK(cur_band_start >= 0 && cur_band_start + cur_num_bands <= Hin.p,
	   "invalid range of bands for computing fields");

     /* first, compute fft_data = curl(Hin) (really (k+G) x H) : */
     for (i = 0; i < d->other_dims; ++i)
	  for (j = 0; j < d->last_dim; ++j) {
	       int ij = i * d->last_dim + j;
	       int ij2 = i * d->last_dim_size + j;
	       k_data cur_k = d->k_plus_G[ij];
	       
	       for (b = 0; b < cur_num_bands; ++b)
		    assign_cross_t2c(&fft_data[3 * (ij2*cur_num_bands 
						    + b)], 
				     cur_k, 
				     &Hin.data[ij * 2 * Hin.p + 
					      b + cur_band_start],
				     Hin.p);
	  }

     /* now, convert to position space via FFT: */
     maxwell_compute_fft(+1, d, fft_data,
			 cur_num_bands*3, cur_num_bands*3, 1);
}

/* Compute E (output in dfield) from D (input in dfield); this amounts
   to just dividing by the dielectric tensor.  dfield is in position
   space and corresponds to the output from maxwell_compute_d_from_H,
   above. */
void maxwell_compute_e_from_d(maxwell_data *d,
			      scalar_complex *dfield,
			      int cur_num_bands)
{
     int i, b;

     CHECK(d, "null maxwell data pointer!");
     CHECK(dfield, "null field input/output data!");

     for (i = 0; i < d->fft_output_size; ++i) {
	  symmetric_matrix eps_inv = d->eps_inv[i];
	  for (b = 0; b < cur_num_bands; ++b) {
	       int ib = 3 * (i * cur_num_bands + b);
	       assign_symmatrix_vector(&dfield[ib], eps_inv, &dfield[ib]);
	  }
     }	  
}

/* Compute the magnetic (H) field in Fourier space from the electric
   field (e) in position space; this amouns to Fourier transforming and
   then taking the curl.  Also, multiply by scale.  Other
   parameters are as in compute_d_from_H. 

   Note: we actually compute (k+G) x E, whereas the actual H field
   is -i/omega i(k+G) x E...so, we are actually computing omega*H, here. */
void maxwell_compute_H_from_e(maxwell_data *d, evectmatrix Hout, 
			      scalar_complex *efield,
			      int cur_band_start, int cur_num_bands,
			      real scale)
{
     scalar *fft_data = (scalar *) efield;
     int i, j, b;

     CHECK(Hout.c == 2, "fields don't have 2 components!");
     CHECK(d, "null maxwell data pointer!");
     CHECK(efield, "null field output data!");
     CHECK(cur_band_start >= 0 && cur_band_start + cur_num_bands <= Hout.p,
	   "invalid range of bands for computing fields");

     /* convert back to Fourier space */
     maxwell_compute_fft(-1, d, fft_data,
			 cur_num_bands*3, cur_num_bands*3, 1);
     
     /* then, compute Hout = curl(fft_data) (* scale factor): */
     
     for (i = 0; i < d->other_dims; ++i)
	  for (j = 0; j < d->last_dim; ++j) {
	       int ij = i * d->last_dim + j;
	       int ij2 = i * d->last_dim_size + j;
	       k_data cur_k = d->k_plus_G[ij];
	       
	       for (b = 0; b < cur_num_bands; ++b)
		    assign_cross_c2t(&Hout.data[ij * 2 * Hout.p + 
					       b + cur_band_start],
				     Hout.p, cur_k, 
				     &fft_data[3 * (ij2*cur_num_bands + b)],
				     scale);
	  }
}


/* Compute H field in position space from Hin.  Parameters and output
   formats are the same as for compute_d_from_H, above. */
void maxwell_compute_h_from_H(maxwell_data *d, evectmatrix Hin, 
			      scalar_complex *hfield,
			      int cur_band_start, int cur_num_bands)
{
     scalar *fft_data = (scalar *) hfield;
     int i, j, b;

     CHECK(Hin.c == 2, "fields don't have 2 components!");
     CHECK(d, "null maxwell data pointer!");
     CHECK(hfield, "null field output data!");
     CHECK(cur_band_start >= 0 && cur_band_start + cur_num_bands <= Hin.p,
	   "invalid range of bands for computing fields");

     /* first, compute fft_data = Hin, with the vector field converted 
	from transverse to cartesian basis: */
     for (i = 0; i < d->other_dims; ++i)
	  for (j = 0; j < d->last_dim; ++j) {
	       int ij = i * d->last_dim + j;
	       int ij2 = i * d->last_dim_size + j;
               k_data cur_k = d->k_plus_G[ij];
	       
	       for (b = 0; b < cur_num_bands; ++b)
		    assign_t2c(&fft_data[3 * (ij2*cur_num_bands 
					      + b)], 
			       cur_k,
			       &Hin.data[ij * 2 * Hin.p + 
					b + cur_band_start],
			       Hin.p);
	  }

     /* now, convert to position space via FFT: */
     maxwell_compute_fft(+1, d, fft_data,
			 cur_num_bands*3, cur_num_bands*3, 1);
}

/**************************************************************************/

/* The following functions take a complex or real vector field
   in position space, as output by rfftwnd (or rfftwnd_mpi), and
   compute the "other half" of the array.

   That is, rfftwnd outputs only half of the logical FFT output,
   since the other half is redundant (see the FFTW manual).  This
   is fine for computation, but for visualization/output we
   want the whole array, redundant or not.  So, we output the array
   in two stages, first outputting the array we are given, then
   using the functions below to compute the other half and outputting
   that.

   Given an array A(i,j,k), the redundant half is given by the
   following identity for transforms of real data:

          A(nx-i,ny-j,nz-k) = A(i,j,k)*

   where nx-i/ny-j/nz-k are interpreted modulo nx/ny/nz.  (This
   means that zero coordinates are handled specially: nx-0 = 0.)

   Note that actually, the other "half" is actually slightly less
   than half of the array.  Note also that the other half, in the
   case of distributed MPI transforms, is not necessarily contiguous,
   due to special handling of zero coordinates.

   There is an additional complication.  The array with the symmetry
   above may have been multiplied by exp(ikx) phases to get its Bloch
   state.  In this case, we must use the identity:

          A(R-x)*exp(ik(R-x)) = ( A(x) exp(ikx) exp(-ikR) )*

   where R is a lattice vector.  That is, we not only must conjugate
   and reverse the order, but we also may need to multiply by exp(-ikR)
   before conjugating.  Unfortunately, R depends upon where we are
   in the array, because of the fact that the indices are interpreted
   modulo n (i.e. the zero indices aren't reordered).  e.g. for
   the point (nx-i,ny-j,nz-k), R = Rx*(i!=0) + Ry*(j!=0) + Rz*(k!=0).

   Another complication is that, for 2d rfftwnd_mpi transforms, the
   truncated dimension (in the transformed, transposed array) is
   the *first* dimension rather than the last.

   This code is a little too subtle for my tastes; real FFTs are a pain. */

#define TWOPI 6.2831853071795864769252867665590057683943388

/* This function takes a complex vector field and replaces it with its
   other "half."  phase{x,y,z} is the phase k*R{x,y,z}, in "units" of
   2*pi.  (Equivalently, phase{x,y,z} is the k vector in the
   reciprocal lattice basis.) */
void maxwell_vectorfield_otherhalf(maxwell_data *d, scalar_complex *field,
				   real phasex, real phasey, real phasez)
{
#ifndef SCALAR_COMPLEX
     int i, j, jmin = 1;
     int rank, n_other, n_last, n_last_stored, n_last_new, nx, ny, nz, nxmax;
#  ifdef HAVE_MPI
     int local_x_start;
#  endif
     scalar_complex pz, pxz, pyz, pxyz;

     nxmax = nx = d->nx; ny = d->ny; nz = d->nz;
     n_other = d->other_dims;
     n_last = d->last_dim;
     n_last_stored = d->last_dim_size / 2;
     n_last_new = n_last - n_last_stored; /* < n_last_stored always */
     rank = (nz == 1) ? (ny == 1 ? 1 : 2) : 3;

#  ifdef HAVE_MPI
     local_x_start = d->local_y_start;
     CHECK(rank == 2 || rank == 3, "unsupported rfftwnd_mpi rank!");
     if (rank == 2) {
	  n_other = nx;
	  n_last_new = ny = d->local_ny;
	  if (local_x_start == 0)
	       --n_last_new; /* DC frequency should not be in other half */
	  else
	       jmin = 0;
	  if (local_x_start + ny == n_last_stored && n_last % 2 == 0)
	       --n_last_new; /* Nyquist freq. should not be in other half */
	  n_last_stored = ny;
     }
     else { /* rank == 3 */
	  ny = nx;
	  nx = d->local_ny;
	  nxmax = local_x_start ? nx - 1 : nx;
	  n_other = nx * ny;
     }
#  endif /* HAVE_MPI */

     /* compute p = exp(i*phase) factors: */
     phasex *= -TWOPI; phasey *= -TWOPI; phasez *= -TWOPI;
     switch (rank) { /* treat z as the last/truncated dimension always */
	 case 3: break;
#  if defined(HAVE_MPI) && ! defined(SCALAR_COMPLEX)
	 case 2: phasez = phasex; phasex = phasey; phasey = 0; break;
#  else
	 case 2: phasez = phasey; phasey = 0; break;
#  endif
	 case 1: phasez = phasex; phasex = phasey = 0; break;
     }
     CASSIGN_SCALAR(pz, cos(phasez), sin(phasez));
     phasex += phasez;
     CASSIGN_SCALAR(pxz, cos(phasex), sin(phasex));
     phasex += phasey;
     CASSIGN_SCALAR(pxyz, cos(phasex), sin(phasex));
     phasey += phasez;
     CASSIGN_SCALAR(pyz, cos(phasey), sin(phasey));

/* convenience macros to copy vectors, vectors times phases, 
   and conjugated vectors: */
#  define ASSIGN_V(f,k,f2,k2) { f[3*(k)+0] = f2[3*(k2)+0];               \
                                f[3*(k)+1] = f2[3*(k2)+1];               \
                                f[3*(k)+2] = f2[3*(k2)+2]; }
#  define ASSIGN_VP(f,k,f2,k2,p) { CASSIGN_MULT(f[3*(k)+0], f2[3*(k2)+0], p); \
                                   CASSIGN_MULT(f[3*(k)+1], f2[3*(k2)+1], p); \
                                   CASSIGN_MULT(f[3*(k)+2], f2[3*(k2)+2], p); }
#  define ASSIGN_CV(f,k,f2,k2) { CASSIGN_CONJ(f[3*(k)+0], f2[3*(k2)+0]); \
                                 CASSIGN_CONJ(f[3*(k)+1], f2[3*(k2)+1]); \
                                 CASSIGN_CONJ(f[3*(k)+2], f2[3*(k2)+2]); }

     /* First, swap the order of elements and multiply by exp(ikR)
        phase factors.  We have to be careful here not to double-swap
        any element pair; this is prevented by never swapping with a
        "conjugated" point that is earlier in the array.  */

     if (rank == 3) {
	  int ix, iy;
	  for (ix = 0; 2*ix <= nxmax; ++ix) {
	       int xdiff, ixc;
#  ifdef HAVE_MPI
	       if (local_x_start == 0) {
		    xdiff = ix != 0; ixc = (nx - ix) % nx;
	       }
	       else {
		    xdiff = 1; ixc = nx-1 - ix;
	       }
#  else
	       xdiff = ix != 0; ixc = (nx - ix) % nx;
#  endif
	       for (iy = 0; iy < ny; ++iy) {
		    int ydiff = iy != 0;
		    int i = ix * ny + iy, ic = ixc * ny + (ny - iy) % ny, jmax;
		    if (ic < i)
			 continue;
		    jmax = n_last_new;
		    if (ic == i)
			 jmax = (jmax + 1) / 2;
		    for (j = 1; j <= jmax; ++j) {
			 int jc = n_last_new + 1 - j;
			 int ij = i*n_last_stored + j;
			 int ijc = ic*n_last_stored + jc;
			 scalar_complex f_tmp[3];
			 switch (xdiff*2 + ydiff) { /* pick exp(-ikR) phase */
			     case 3: /* xdiff && ydiff */
				  ASSIGN_VP(f_tmp, 0, field, ijc, pxyz);
				  ASSIGN_VP(field, ijc, field, ij, pxyz);
				  ASSIGN_V(field, ij, f_tmp, 0);
				  break;
			     case 2: /* xdiff && !ydiff */
				  ASSIGN_VP(f_tmp, 0, field, ijc, pxz);
				  ASSIGN_VP(field, ijc, field, ij, pxz);
				  ASSIGN_V(field, ij, f_tmp, 0);
				  break;
			     case 1: /* !xdiff && ydiff */
				  ASSIGN_VP(f_tmp, 0, field, ijc, pyz);
				  ASSIGN_VP(field, ijc, field, ij, pyz);
				  ASSIGN_V(field, ij, f_tmp, 0);
				  break;
			     case 0: /* !xdiff && !ydiff */
				  ASSIGN_VP(f_tmp, 0, field, ijc, pz);
				  ASSIGN_VP(field, ijc, field, ij, pz);
				  ASSIGN_V(field, ij, f_tmp, 0);
				  break;
			 }
		    }
	       }
	  }

	  /* Next, conjugate, and remove the holes from the array
	     corresponding to the DC and Nyquist frequencies (which were in
	     the first half already): */
	  for (i = 0; i < n_other; ++i)
	       for (j = 1; j < n_last_new + 1; ++j) {
		    int ij = i*n_last_stored + j, ijnew = i*n_last_new + j-1;
		    ASSIGN_CV(field, ijnew, field, ij);
	       }
     }
     else /* if (rank <= 2) */ {
	  int i;
	  if (rank == 1) /* (note that 1d MPI transforms are not allowed) */
	       nx = 1; /* x dimension is handled by j (last dimension) loop */

#  ifdef HAVE_MPI
	  for (i = 0; i < nx; ++i)
#  else
	  for (i = 0; 2*i <= nx; ++i)
#  endif
	  {
	       int xdiff = i != 0, ic = (nx - i) % nx;
	       int jmax = n_last_new + (jmin - 1);
#  ifndef HAVE_MPI
	       if (ic == i)
		    jmax = (jmax + 1) / 2;
#  endif
	       for (j = jmin; j <= jmax; ++j) {
		    scalar_complex f_tmp[3];
#  ifdef HAVE_MPI
		    int jc = jmax + jmin - j;
		    int ij = j * nx + i;
		    int ijc = jc * nx + ic;
		    if (ijc < ij)
			 break;
#  else  /* ! HAVE_MPI */
		    int jc = n_last_new + 1 - j;
		    int ij = i*n_last_stored + j;
		    int ijc = ic*n_last_stored + jc;
#  endif /* ! HAVE_MPI */
		    if (xdiff) {
			 ASSIGN_VP(f_tmp, 0, field, ijc, pxz);
			 ASSIGN_VP(field, ijc, field, ij, pxz);
			 ASSIGN_V(field, ij, f_tmp, 0);
		    }
		    else {
			 ASSIGN_VP(f_tmp, 0, field, ijc, pz);
			 ASSIGN_VP(field, ijc, field, ij, pz);
			 ASSIGN_V(field, ij, f_tmp, 0);
		    }
	       }
	  }

	  /* Next, conjugate, and remove the holes from the array
	     corresponding to the DC and Nyquist frequencies (which were in
	     the first half already): */
	  for (i = 0; i < nx; ++i)
	       for (j = jmin; j < n_last_new + jmin; ++j) {
#  ifdef HAVE_MPI
		    int ij = j*nx + i, ijnew = (j-jmin)*nx + i;
#  else
		    int ij = i*n_last_stored + j, ijnew = i*n_last_new + j-1;
#  endif
		    ASSIGN_CV(field, ijnew, field, ij);
	       }
     }

#  undef ASSIGN_V
#  undef ASSIGN_VP
#  undef ASSIGN_CV

#endif /* ! SCALAR_COMPLEX */
}

/* Similar to vectorfield_otherhalf, above, except that it operates on
   a real scalar field, which is assumed to have come from e.g. the
   absolute values of a complex field (and thus no phase factors or
   conjugations are necessary). */
void maxwell_scalarfield_otherhalf(maxwell_data *d, real *field)
{
#ifndef SCALAR_COMPLEX
     int i, j, jmin = 1;
     int rank, n_other, n_last, n_last_stored, n_last_new, nx, ny, nz, nxmax;
#  ifdef HAVE_MPI
     int local_x_start;
#  endif

     nxmax = nx = d->nx; ny = d->ny; nz = d->nz;
     n_other = d->other_dims;
     n_last = d->last_dim;
     n_last_stored = d->last_dim_size / 2;
     n_last_new = n_last - n_last_stored; /* < n_last_stored always */
     rank = (nz == 1) ? (ny == 1 ? 1 : 2) : 3;

#  ifdef HAVE_MPI
     local_x_start = d->local_y_start;
     CHECK(rank == 2 || rank == 3, "unsupported rfftwnd_mpi rank!");
     if (rank == 2) {
	  n_other = nx;
	  n_last_new = ny = d->local_ny;
	  if (local_x_start == 0)
	       --n_last_new; /* DC frequency should not be in other half */
	  else
	       jmin = 0;
	  if (local_x_start + ny == n_last_stored && n_last % 2 == 0)
	       --n_last_new; /* Nyquist freq. should not be in other half */
	  n_last_stored = ny;
     }
     else { /* rank == 3 */
	  ny = nx;
	  nx = d->local_ny;
	  nxmax = local_x_start ? nx - 1 : nx;
	  n_other = nx * ny;
     }
#  endif /* HAVE_MPI */

     /* First, swap the order of elements and multiply by exp(ikR)
        phase factors.  We have to be careful here not to double-swap
        any element pair; this is prevented by never swapping with a
        "conjugated" point that is earlier in the array.  */

     if (rank == 3) {
	  int ix, iy;
	  for (ix = 0; 2*ix <= nxmax; ++ix) {
	       int ixc;
#  ifdef HAVE_MPI
	       if (local_x_start == 0)
		    ixc = (nx - ix) % nx;
	       else
		    ixc = nx-1 - ix;
#  else
	       ixc = (nx - ix) % nx;
#  endif
	       for (iy = 0; iy < ny; ++iy) {
		    int i = ix * ny + iy, ic = ixc * ny + (ny - iy) % ny, jmax;
		    if (ic < i)
			 continue;
		    jmax = n_last_new;
		    if (ic == i)
			 jmax = (jmax + 1) / 2;
		    for (j = 1; j <= jmax; ++j) {
			 int jc = n_last_new + 1 - j;
			 int ij = i*n_last_stored + j;
			 int ijc = ic*n_last_stored + jc;
			 real f_tmp;
			 f_tmp = field[ijc];
			 field[ijc] = field[ij];
			 field[ij] = f_tmp;
		    }
	       }
	  }

	  /* Next, conjugate, and remove the holes from the array
	     corresponding to the DC and Nyquist frequencies (which were in
	     the first half already): */
	  for (i = 0; i < n_other; ++i)
	       for (j = 1; j < n_last_new + 1; ++j) {
		    int ij = i*n_last_stored + j, ijnew = i*n_last_new + j-1;
		    field[ijnew] = field[ij];
	       }
     }
     else /* if (rank <= 2) */ {
	  int i;
	  if (rank == 1) /* (note that 1d MPI transforms are not allowed) */
	       nx = 1; /* x dimension is handled by j (last dimension) loop */

#  ifdef HAVE_MPI
	  for (i = 0; i < nx; ++i)
#  else
	  for (i = 0; 2*i <= nx; ++i)
#  endif
	  {
	       int ic = (nx - i) % nx;
	       int jmax = n_last_new + (jmin - 1);
#  ifndef HAVE_MPI
	       if (ic == i)
		    jmax = (jmax + 1) / 2;
#  endif
	       for (j = jmin; j <= jmax; ++j) {
		    real f_tmp;
#  ifdef HAVE_MPI
		    int jc = jmax + jmin - j;
		    int ij = j * nx + i;
		    int ijc = jc * nx + ic;
		    if (ijc < ij)
			 break;
#  else  /* ! HAVE_MPI */
		    int jc = n_last_new + 1 - j;
		    int ij = i*n_last_stored + j;
		    int ijc = ic*n_last_stored + jc;
#  endif /* ! HAVE_MPI */
		    f_tmp = field[ijc];
		    field[ijc] = field[ij];
		    field[ij] = f_tmp;
	       }
	  }

	  /* Next, remove the holes from the array corresponding to
	     the DC and Nyquist frequencies (which were in the first
	     half already): */
	  for (i = 0; i < nx; ++i)
	       for (j = jmin; j < n_last_new + jmin; ++j) {
#  ifdef HAVE_MPI
		    int ij = j*nx + i, ijnew = (j-jmin)*nx + i;
#  else
		    int ij = i*n_last_stored + j, ijnew = i*n_last_new + j-1;
#  endif
		    field[ijnew] = field[ij];
	       }
     }
#endif /* ! SCALAR_COMPLEX */
}

/**************************************************************************/

#define MIN2(a,b) ((a) < (b) ? (a) : (b))

/* Compute Xout = curl(1/epsilon * curl(Xin)) */
void maxwell_operator(evectmatrix Xin, evectmatrix Xout, void *data,
		      int is_current_eigenvector, evectmatrix Work)
{
     maxwell_data *d = (maxwell_data *) data;
     int cur_band_start;
     scalar_complex *cdata;
     real scale;
     
     CHECK(d, "null maxwell data pointer!");
     CHECK(Xin.c == 2, "fields don't have 2 components!");

     (void) is_current_eigenvector;  /* unused */
     (void) Work;

     cdata = (scalar_complex *) d->fft_data;
     scale = -1.0 / Xout.N;  /* scale factor to normalize FFT; 
				negative sign comes from 2 i's from curls */

     /* compute the operator, num_fft_bands at a time: */
     for (cur_band_start = 0; cur_band_start < Xin.p; 
	  cur_band_start += d->num_fft_bands) {
	  int cur_num_bands = MIN2(d->num_fft_bands, Xin.p - cur_band_start);

	  maxwell_compute_d_from_H(d, Xin, cdata,
				   cur_band_start, cur_num_bands);
	  maxwell_compute_e_from_d(d, cdata, cur_num_bands);
	  maxwell_compute_H_from_e(d, Xout, cdata,
				   cur_band_start, cur_num_bands, scale);
     }
}

/* Compute the operation Xout = (M - w^2) Xin, where M is the Maxwell
   operator and w is the target frequency.  This shifts the eigenvalue
   spectrum so that the smallest magnitude eigenvalues are those
   nearest to w.  However, there are negative eigenvalues (the
   operator is not positive-definite), and the smallest eigenvectors
   (not taking the absolute value) are the same as those of M. */
void maxwell_target_operator1(evectmatrix Xin, evectmatrix Xout, void *data,
			     int is_current_eigenvector, evectmatrix Work)
{
     maxwell_target_data *d = (maxwell_target_data *) data;
     real omega_sqr = d->target_frequency * d->target_frequency;

     maxwell_operator(Xin, Xout, d->d, is_current_eigenvector, Work);
     evectmatrix_aXpbY(1.0, Xout, -omega_sqr, Xin);
}

/* Compute the operation Xout = (M - w^2)^2 Xin, where M is the
   Maxwell operator and w is the target frequency.  This shifts the
   eigenvalue spectrum so that the smallest eigenvalues are those
   nearest to w. */
void maxwell_target_operator(evectmatrix Xin, evectmatrix Xout, void *data,
			     int is_current_eigenvector, evectmatrix Work)
{
     if (Xin.n != 0)
	  CHECK(Work.data && Work.data != Xin.data && Work.data != Xout.data,
		"maxwell_target_operator must have distinct workspace!");

     maxwell_target_operator1(Xin, Work, data, is_current_eigenvector, Xout);

     /* N.B. maxwell_operator(), and thus maxwell_target_operator1(),
	doesn't actually need the workspace, so we can safely pass
	Work here for the scratch parameter: */
     maxwell_target_operator1(Work, Xout, data, is_current_eigenvector, Work);
}

/* Compute the operation Xout = curl 1/epsilon * i u x Xin, which 
   is useful operation in computing the group velocity (derivative
   of the maxwell operator).  u is a vector in cartesian coordinates. */
void maxwell_ucross_op(evectmatrix Xin, evectmatrix Xout,
		       maxwell_data *d, const real u[3])
{
     scalar *fft_data;
     scalar_complex *cdata;
     real scale;
     int cur_band_start;
     int i, j, b;

     CHECK(d, "null maxwell data pointer!");
     CHECK(Xin.c == 2, "fields don't have 2 components!");

     cdata = (scalar_complex *) (fft_data = d->fft_data);
     scale = -1.0 / Xout.N;  /* scale factor to normalize FFT;
                                negative sign comes from 2 i's from curls */

     /* compute the operator, num_fft_bands at a time: */
     for (cur_band_start = 0; cur_band_start < Xin.p;
          cur_band_start += d->num_fft_bands) {
          int cur_num_bands = MIN2(d->num_fft_bands, Xin.p - cur_band_start);
	  
	  /* first, compute fft_data = u x Xin: */
	  for (i = 0; i < d->other_dims; ++i)
	       for (j = 0; j < d->last_dim; ++j) {
		    int ij = i * d->last_dim + j;
		    int ij2 = i * d->last_dim_size + j;
		    k_data cur_k = d->k_plus_G[ij];
		    
		    for (b = 0; b < cur_num_bands; ++b)
			 assign_ucross_t2c(&fft_data[3 * (ij2*cur_num_bands
							  + b)], 
					   u, cur_k, 
					   &Xin.data[ij * 2 * Xin.p + 
						    b + cur_band_start],
					   Xin.p);
	       }
	  
	  /* now, convert to position space via FFT: */
	  maxwell_compute_fft(+1, d, fft_data,
			      cur_num_bands*3, cur_num_bands*3, 1);
	  
          maxwell_compute_e_from_d(d, cdata, cur_num_bands);
          maxwell_compute_H_from_e(d, Xout, cdata,
                                   cur_band_start, cur_num_bands, scale);
     }
}
