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

#include "../src/config.h"
#include <check.h>

#include "elastic.h"

/**************************************************************************/

void elastic_compute_fft(int dir, elastic_data *d, scalar *array, 
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

#define ELASTICFFT F77_FUNC(elasticfft,ELASTICFFT)

void ELASTICFFT(int *dir, elastic_data *d, scalar *array,
	   int *howmany, int *stride, int *dist)
{
     elastic_compute_fft(*dir, d, array, *howmany, *stride, *dist);
}

/**************************************************************************/

/* Compute u (displacement field) from v; this amounts to dividing v by
   sqrt(rho).  Input and output are in field. */
void elastic_compute_u_from_v(elastic_data *d,
			      scalar_complex *field,
			      int cur_num_bands)
{
     int i, b;

     CHECK(d, "null elastic data pointer!");
     CHECK(field, "null field input/output data!");

     for (i = 0; i < d->fft_output_size; ++i) {
	  real sqrt_rho_inv = d->sqrt_rhoinv[i];
	  for (b = 0; b < cur_num_bands; ++b) {
	       int ib = 3 * (i * cur_num_bands + b);
	       int k;
	       for (k = 0; k < 3; ++k)
		    CASSIGN_SCALAR(field[ib + k],
				   CSCALAR_RE(field[ib + k]) * sqrt_rho_inv,
				   CSCALAR_IM(field[ib + k]) * sqrt_rho_inv)
	  }
     }	  
}

/* Compute v field in position space from Vin.  */
void elastic_compute_v_from_V(elastic_data *d, evectmatrix Vin, 
			      scalar_complex *vfield,
			      int cur_band_start, int cur_num_bands)
{
     scalar *fft_data = (scalar *) vfield;
     int i, j, b;

     CHECK(Vin.c == 3, "fields don't have 3 components!");
     CHECK(d, "null elastic data pointer!");
     CHECK(vfield, "null field output data!");
     CHECK(cur_num_bands <= d->num_fft_bands, "too many bands to FFT at once");
     CHECK(cur_band_start >= 0 && cur_band_start + cur_num_bands <= Vin.p,
	   "invalid range of bands for computing fields");

     /* first, compute fft_data = Vin, with the vector field converted 
	from transverse to cartesian basis: */
     for (i = 0; i < d->other_dims; ++i)
	  for (j = 0; j < d->last_dim; ++j) {
	       int ij = i * d->last_dim + j;
	       int ij2 = i * d->last_dim_size + j;
	       int k;

	       for (b = 0; b < cur_num_bands; ++b) 
		    for (k = 0; k < 3; ++k)
			 fft_data[3 * (ij2*cur_num_bands + b) + k] =
			      Vin.data[(ij * 3 + k) * Vin.p + 
				      b + cur_band_start];
	  }

     /* now, convert to position space via FFT: */
     elastic_compute_fft(+1, d, fft_data,
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
void elastic_vectorfield_otherhalf(elastic_data *d, scalar_complex *field,
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
void elastic_scalarfield_otherhalf(elastic_data *d, real *field)
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

#define AV F77_FUNC(av,AV)
extern void AV(int *p, int *ldp,
	       scalar *v, scalar *Av,
	       elastic_data *edata,
	       real *sqrtrhoinv, real *rhoct2, real *rhocl2,
	       real *vkx, real *vky, real *vkz,
	       int *Ni, int *Nj, int *Nk,
	       real *b1x, real *b1y, real *b1z,
	       real *b2x, real *b2y, real *b2z,
	       real *b3x, real *b3y, real *b3z,
	       scalar_complex *Datax,
	       scalar_complex *Datay,
	       scalar_complex *Dataz,
	       scalar_complex *Data1,
	       scalar_complex *Data2,
	       scalar_complex *Data3,
	       scalar_complex *Data4,
	       scalar_complex *Data5,
	       scalar_complex *Data6,
	       scalar_complex *Data7,
	       scalar_complex *Data8,
	       scalar_complex *Data9,
	       scalar_complex *Data10,
	       scalar_complex *Data11,
	       scalar_complex *Data12,
	       scalar_complex *Data13,
	       scalar_complex *Data14,
	       scalar_complex *Data15,
	       scalar_complex *Data16,
	       scalar_complex *Data17,
	       scalar_complex *Data18);

/* Compute Xout = OP * Xin, where OP is the elastic-band eigen-operator
   defined in operator.f. */
void elastic_operator(evectmatrix Xin, evectmatrix Xout, void *data,
		      int is_current_eigenvector, evectmatrix Work)
{
     elastic_data *d = (elastic_data *) data;
     int cur_band_start;
     scalar_complex *cdata;
     
     CHECK(d, "null elastic data pointer!");
     CHECK(Xin.c == 3, "fields don't have 3 components!");

     (void) is_current_eigenvector;  /* unused */
     (void) Work; /* unused */
     cdata = (scalar_complex *) d->fft_data;

     /* compute the operator, num_fft_bands at a time: */
     for (cur_band_start = 0; cur_band_start < Xin.p; 
	  cur_band_start += d->num_fft_bands) {
	  int cur_num_bands = MIN2(d->num_fft_bands, Xin.p - cur_band_start);

	  AV(&cur_num_bands, &Xin.p,
	     Xin.data + cur_band_start, Xout.data + cur_band_start,
	     d,
	     d->sqrt_rhoinv, d->rhoct2, d->rhocl2,
	     &d->current_k[0], &d->current_k[1], &d->current_k[2],
	     &d->nx, &d->ny, &d->nz,
	     &d->G[0][0], &d->G[0][1], &d->G[0][2],
	     &d->G[1][0], &d->G[1][1], &d->G[1][2],
	     &d->G[2][0], &d->G[2][1], &d->G[2][2],

	     /* FIXME: in the future, it should be possible to just
		use the Work array and the fft_data array, and to
	        also use the FFTW plans (and above compute_v_from_V etc) */
	     d->crap + Xin.N * 0,
	     d->crap + Xin.N * 1,
	     d->crap + Xin.N * 2,
	     d->crap + Xin.N * 3,
	     d->crap + Xin.N * 4,
	     d->crap + Xin.N * 5,
	     d->crap + Xin.N * 6,
	     d->crap + Xin.N * 7,
	     d->crap + Xin.N * 8,
	     d->crap + Xin.N * 9,
	     d->crap + Xin.N * 10,
	     d->crap + Xin.N * 11,
	     d->crap + Xin.N * 12,
	     d->crap + Xin.N * 13,
	     d->crap + Xin.N * 14,
	     d->crap + Xin.N * 15,
	     d->crap + Xin.N * 16,
	     d->crap + Xin.N * 17,
	     d->crap + Xin.N * 18,
	     d->crap + Xin.N * 19,
	     d->crap + Xin.N * 20);
     }
}

/**************************************************************************/

#define PRECOND F77_FUNC(precond,PRECOND)
extern void PRECOND(int *p, int *ldp,
	       scalar *v, scalar *Av,
	       elastic_data *edata,
	       real *sqrtrhoinv, real *rhoct2, real *rhocl2,
	       real *vkx, real *vky, real *vkz,
	       int *Ni, int *Nj, int *Nk,
	       real *b1x, real *b1y, real *b1z,
	       real *b2x, real *b2y, real *b2z,
	       real *b3x, real *b3y, real *b3z,
	       scalar_complex *Datax,
	       scalar_complex *Datay,
	       scalar_complex *Dataz,
	       scalar_complex *Data1,
	       scalar_complex *Data2,
	       scalar_complex *Data3,
	       scalar_complex *Data4,
	       scalar_complex *Data5,
	       scalar_complex *Data6,
	       scalar_complex *Data7,
	       scalar_complex *Data8,
	       scalar_complex *Data9,
	       scalar_complex *Data10,
	       scalar_complex *Data11,
	       scalar_complex *Data12,
	       scalar_complex *Data13,
	       scalar_complex *Data14,
	       scalar_complex *Data15,
	       scalar_complex *Data16,
	       scalar_complex *Data17,
	       scalar_complex *Data18);

/* Compute Xout = OP * Xin, where OP is the elastic-band eigen-operator preconditioner
   defined in operator.f. */
void elastic_preconditioner(evectmatrix Xin, evectmatrix Xout, void *data,
			    evectmatrix Y, real *eigenvals,
			    sqmatrix YtY)
{
     elastic_data *d = (elastic_data *) data;
     int cur_band_start;
     scalar_complex *cdata;
     
     CHECK(d, "null elastic data pointer!");
     CHECK(Xin.c == 3, "fields don't have 3 components!");

     (void) Y; /* unused */
     (void) eigenvals;
     (void) YtY;
     cdata = (scalar_complex *) d->fft_data;

     /* compute the operator, num_fft_bands at a time: */
     for (cur_band_start = 0; cur_band_start < Xin.p; 
	  cur_band_start += d->num_fft_bands) {
	  int cur_num_bands = MIN2(d->num_fft_bands, Xin.p - cur_band_start);

	  PRECOND(&cur_num_bands, &Xin.p,
		  Xin.data + cur_band_start, Xout.data + cur_band_start,
		  d,
		  d->sqrt_rhoinv, d->rhoct2, d->rhocl2,
		  &d->current_k[0], &d->current_k[1], &d->current_k[2],
		  &d->nx, &d->ny, &d->nz,
		  &d->G[0][0], &d->G[0][1], &d->G[0][2],
		  &d->G[1][0], &d->G[1][1], &d->G[1][2],
		  &d->G[2][0], &d->G[2][1], &d->G[2][2],
		  
		  /* FIXME: in the future, it should be possible to just
		     use the Work array and the fft_data array, and to
		     also use the FFTW plans (and above compute_v_from_V etc) */
		  d->crap + Xin.N * 0,
		  d->crap + Xin.N * 1,
		  d->crap + Xin.N * 2,
		  d->crap + Xin.N * 3,
		  d->crap + Xin.N * 4,
		  d->crap + Xin.N * 5,
		  d->crap + Xin.N * 6,
		  d->crap + Xin.N * 7,
		  d->crap + Xin.N * 8,
		  d->crap + Xin.N * 9,
		  d->crap + Xin.N * 10,
		  d->crap + Xin.N * 11,
		  d->crap + Xin.N * 12,
		  d->crap + Xin.N * 13,
		  d->crap + Xin.N * 14,
		  d->crap + Xin.N * 15,
		  d->crap + Xin.N * 16,
		  d->crap + Xin.N * 17,
		  d->crap + Xin.N * 18,
		  d->crap + Xin.N * 19,
		  d->crap + Xin.N * 20);
     }
}

