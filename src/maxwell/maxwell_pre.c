/* Copyright (C) 1999-2020 Massachusetts Institute of Technology.
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

#include <mpiglue.h>
#include "maxwell.h"

#define PRECOND_SUBTR_EIGS 0

#define PRECOND_MIN_DENOM 5e-3
#define MIN2(a, b) ((a) < (b) ? (a) : (b))
#define MAX2(a, b) ((a) > (b) ? (a) : (b))

/* used to be: MAX2(x, PRECOND_MIN_DENOM) but now zero k-point is
   handled specially */
#define FIX_DENOM(x) ((x) == 0 ? 1.0 : (x))

void maxwell_simple_precondition(evectmatrix X, void *data, real *eigenvals) {
  maxwell_data *d = (maxwell_data *)data;
  int i, c, b;
  real *kpGn2 = d->k_plus_G_normsqr;

#if !PRECOND_SUBTR_EIGS
  (void)eigenvals; /* unused */
#endif

  for (i = 0; i < X.localN; ++i) {
    for (c = 0; c < X.c; ++c) {
      for (b = 0; b < X.p; ++b) {
        int index = (i * X.c + c) * X.p + b;
        real scale = kpGn2[i] * d->eps_inv_mean;

#if PRECOND_SUBTR_EIGS
        if (eigenvals) {
          scale -= eigenvals[b];
          scale = 1.0 / (scale + copysign(PRECOND_MIN_DENOM, scale));
        }
        else
#else
        { scale = 1.0 / FIX_DENOM(scale); }
#endif

          ASSIGN_SCALAR(X.data[index], scale * SCALAR_RE(X.data[index]),
                        scale * SCALAR_IM(X.data[index]));
      }
    }
  }
}

void maxwell_preconditioner(evectmatrix Xin, evectmatrix Xout, void *data, evectmatrix Y,
                            real *eigenvals, sqmatrix YtY) {
  (void)Y; /* unused */
  evectmatrix_XeYS(Xout, Xin, YtY, 1);
  maxwell_simple_precondition(Xout, data, eigenvals);
}

void maxwell_target_preconditioner(evectmatrix Xin, evectmatrix Xout, void *data, evectmatrix Y,
                                   real *eigenvals, sqmatrix YtY) {
  maxwell_target_data *td = (maxwell_target_data *)data;
  maxwell_data *d = td->d;
#if PRECOND_SUBTR_EIGS
  real omega_sqr = td->target_frequency * td->target_frequency;
#endif
  int i, c, b;
  real *kpGn2 = d->k_plus_G_normsqr;

  (void)Y; /* unused */
#if !PRECOND_SUBTR_EIGS
  (void)eigenvals; /* unused */
#endif

  evectmatrix_XeYS(Xout, Xin, YtY, 1);

  for (i = 0; i < Xout.localN; ++i) {
    for (c = 0; c < Xout.c; ++c) {
      for (b = 0; b < Xout.p; ++b) {
        int index = (i * Xout.c + c) * Xout.p + b;
        real scale = kpGn2[i] * d->eps_inv_mean;

#if PRECOND_SUBTR_EIGS
        scale -= omega_sqr;
#endif

        scale = scale * scale;

#if PRECOND_SUBTR_EIGS
        if (eigenvals) {
          scale -= eigenvals[b];
          scale = 1.0 / (scale + copysign(PRECOND_MIN_DENOM, scale));
        }
        else
#else
        { scale = 1.0 / FIX_DENOM(scale); }
#endif

          ASSIGN_SCALAR(Xout.data[index], scale * SCALAR_RE(Xout.data[index]),
                        scale * SCALAR_IM(Xout.data[index]));
      }
    }
  }
}

/**************************************************************************/

/* Fancy preconditioners */

/* Compute 'a' where v = k x a (cross product).

   Here, a = (a[0],a[1],a[2]) and k = (k.kx,k.ky,k.kz) are in
   cartesian coordinates.  (v[0],v[vstride]) is in the transverse basis of
   k.m and k.n.

   We can't compute 'a' exactly, since there is no way to find the
   component of a parallel to k.  So, we only compute the transverse
   component of 'a'--this is the main approximation in our preconditioner.
*/
static void assign_crossinv_t2c(scalar *a, const k_data k, const scalar *v, int vstride) {
  /* k x v = k x (k x a) = (k*a)k - k^2 a
           = -(a_transverse) * k^2
     ==> a_transverse = -1/k^2 * k x v     */

  /* Thus, we just do the same thing as assign_cross_t2c
     in maxwell_op.c, except that we divide by -k^2: */

  scalar v0 = v[0], v1 = v[vstride];
  real kmag_inv = -1.0 / FIX_DENOM(k.kmag);

  ASSIGN_SCALAR(a[0], (SCALAR_RE(v0) * k.nx - SCALAR_RE(v1) * k.mx) * kmag_inv,
                (SCALAR_IM(v0) * k.nx - SCALAR_IM(v1) * k.mx) * kmag_inv);
  ASSIGN_SCALAR(a[1], (SCALAR_RE(v0) * k.ny - SCALAR_RE(v1) * k.my) * kmag_inv,
                (SCALAR_IM(v0) * k.ny - SCALAR_IM(v1) * k.my) * kmag_inv);
  ASSIGN_SCALAR(a[2], (SCALAR_RE(v0) * k.nz - SCALAR_RE(v1) * k.mz) * kmag_inv,
                (SCALAR_IM(v0) * k.nz - SCALAR_IM(v1) * k.mz) * kmag_inv);
}

/* Compute 'v' * scale, where a = k x v, going from cartesian to transverse
   coordinates.  Since v is tranvserse to k, we can compute this inverse
   exactly. */
static void assign_crossinv_c2t(scalar *v, int vstride, const k_data k, const scalar *a,
                                real scale) {
  /* As in assign_crossinv_t2c above, we find:

     v = v_transverse = -1/k^2 * k x a

     So, we do the same thing as in assign_cross_c2t of maxwell_op.c,
     with the additional -1/k^2 factor. */

  scalar a0 = a[0], a1 = a[1], a2 = a[2];
  scalar at0, at1;

  ASSIGN_SCALAR(at0, SCALAR_RE(a0) * k.mx + SCALAR_RE(a1) * k.my + SCALAR_RE(a2) * k.mz,
                SCALAR_IM(a0) * k.mx + SCALAR_IM(a1) * k.my + SCALAR_IM(a2) * k.mz);
  ASSIGN_SCALAR(at1, SCALAR_RE(a0) * k.nx + SCALAR_RE(a1) * k.ny + SCALAR_RE(a2) * k.nz,
                SCALAR_IM(a0) * k.nx + SCALAR_IM(a1) * k.ny + SCALAR_IM(a2) * k.nz);

  /* combine scale factor and k * (-1/k^2) */
  scale = -scale / FIX_DENOM(k.kmag);

  ASSIGN_SCALAR(v[0], -scale * SCALAR_RE(at1), -scale * SCALAR_IM(at1));
  ASSIGN_SCALAR(v[vstride], scale * SCALAR_RE(at0), scale * SCALAR_IM(at0));
}

/* Fancy preconditioner.  This is very similar to maxwell_op, except that
   the steps are (approximately) inverted: */

void maxwell_preconditioner2(evectmatrix Xin, evectmatrix Xout, void *data, evectmatrix Y,
                             real *eigenvals, sqmatrix YtY) {
  maxwell_data *d = (maxwell_data *)data;
  int cur_band_start;
  scalar *fft_data, *fft_data2;
  scalar_complex *cdata;
  real scale;
  int i, j, b;

  (void)Y;         /* unused */
  (void)eigenvals; /* unused */

  CHECK(d, "null maxwell data pointer!");
  CHECK(Xin.c == 2, "fields don't have 2 components!");

  if (Xout.data != Xin.data) evectmatrix_XeYS(Xout, Xin, YtY, 1);

  fft_data = d->fft_data;
  fft_data2 = d->fft_data2;
  cdata = (scalar_complex *)fft_data;

  scale = -1.0 / Xout.N; /* scale factor to normalize FFT;
                            negative sign comes from 2 i's from curls */

  for (cur_band_start = 0; cur_band_start < Xout.p; cur_band_start += d->num_fft_bands) {
    int cur_num_bands = MIN2(d->num_fft_bands, Xout.p - cur_band_start);

    /********************************************/
    /* Compute approx. inverse of curl (inverse cross product with k): */

    for (i = 0; i < d->other_dims; ++i)
      for (j = 0; j < d->last_dim; ++j) {
        int ij = i * d->last_dim + j;
        int ij2 = i * d->last_dim_size + j;
        k_data cur_k = d->k_plus_G[ij];

        for (b = 0; b < cur_num_bands; ++b)
          assign_crossinv_t2c(&fft_data2[3 * (ij2 * cur_num_bands + b)], cur_k,
                              &Xout.data[ij * 2 * Xout.p + b + cur_band_start], Xout.p);
      }

    /********************************************/
    /* Multiply by epsilon: */

    /* convert to position space via FFT: */
    maxwell_compute_fft(+1, d, fft_data2, fft_data, cur_num_bands * 3, cur_num_bands * 3, 1);

    /* multiply by epsilon in position space.  Don't bother to
       invert the whole epsilon-inverse tensor; just take the
       inverse of the average epsilon-inverse (= trace / 3). */
    for (i = 0; i < d->fft_output_size; ++i) {
      symmetric_matrix eps_inv = d->eps_inv[i];
      real eps = 3.0 / (eps_inv.m00 + eps_inv.m11 + eps_inv.m22);
      for (b = 0; b < cur_num_bands; ++b) {
        int ib = 3 * (i * cur_num_bands + b);
        cdata[ib].re *= eps;
        cdata[ib].im *= eps;
        cdata[ib + 1].re *= eps;
        cdata[ib + 1].im *= eps;
        cdata[ib + 2].re *= eps;
        cdata[ib + 2].im *= eps;
      }
    }

    /* convert back to Fourier space */
    maxwell_compute_fft(-1, d, fft_data, fft_data2, cur_num_bands * 3, cur_num_bands * 3, 1);

    /********************************************/
    /* Finally, do second inverse curl (inverse cross product with k): */

    for (i = 0; i < d->other_dims; ++i)
      for (j = 0; j < d->last_dim; ++j) {
        int ij = i * d->last_dim + j;
        int ij2 = i * d->last_dim_size + j;
        k_data cur_k = d->k_plus_G[ij];

        for (b = 0; b < cur_num_bands; ++b)
          assign_crossinv_c2t(&Xout.data[ij * 2 * Xout.p + b + cur_band_start], Xout.p, cur_k,
                              &fft_data2[3 * (ij2 * cur_num_bands + b)], scale);
      }

    /********************************************/

  } /* end of cur_band_start loop */
}

void maxwell_target_preconditioner2(evectmatrix Xin, evectmatrix Xout, void *data, evectmatrix Y,
                                    real *eigenvals, sqmatrix YtY) {
  maxwell_target_data *d = (maxwell_target_data *)data;
  maxwell_preconditioner2(Xin, Xout, d->d, Y, eigenvals, YtY);
  maxwell_preconditioner2(Xout, Xout, d->d, Y, eigenvals, YtY);
}
