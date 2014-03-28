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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stddef.h>

#include "config.h"
#include <check.h>
#include <matrixio.h>
#include <mpiglue.h>
#include <mpi_utils.h>
#include <maxwell.h>

#include <ctl-io.h>
#include <ctlgeom.h>

#include "mpb.h"
#include "field-smob.h"

/**************************************************************************/

/* The following routines take the eigenvectors computed by solve-kpoint
   and compute the field (D, H, or E) in position space for one of the bands.
   This field is stored in the global curfield (actually an alias for
   mdata->fft_data, since the latter is unused and big enough).  This
   field can then be manipulated with subsequent "*-field-*" functions
   below.  You can also get the scalar field, epsilon.

   All of these functions are designed to be called by the user
   via Guile. */

void get_dfield(int which_band)
{
     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-dfield!\n");
	  return;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr,
			  "solve-kpoint must be called before get-dfield!\n");
	  return;
     }
     if (which_band < 1 || which_band > H.p) {
	  mpi_one_fprintf(stderr,
			  "must have 1 <= band index <= num_bands (%d)\n",H.p);
	  return;
     }

     curfield = (scalar_complex *) mdata->fft_data;
     curfield_band = which_band;
     curfield_type = 'd';
     maxwell_compute_d_from_H(mdata, H, curfield, which_band - 1, 1);

     /* Here, we correct for the fact that compute_d_from_H actually
	computes just (k+G) x H, whereas the actual D field is
	i/omega i(k+G) x H...so, there is an added factor of -1/omega. 

        We also divide by the cell volume so that the integral of |H|^2
        or of D*E is unity.  (From the eigensolver + FFT, they are
        initially normalized to sum to nx*ny*nz.) */
     {
	  int i, N;
	  double scale;
	  N = mdata->fft_output_size;

	  if (freqs.items[which_band - 1] != 0.0) {
	       scale = -1.0 / freqs.items[which_band - 1];
	  }
	  else
	       scale = -1.0; /* arbitrary */

	  scale /= sqrt(Vol);

	  for (i = 0; i < 3*N; ++i) {
	       curfield[i].re *= scale;
	       curfield[i].im *= scale;
	  }
     }
}

void get_hfield(integer which_band)
{
     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-hfield!\n");
	  return;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr,
			  "solve-kpoint must be called before get-hfield!\n");
	  return;
     }
     if (which_band < 1 || which_band > H.p) {
	  mpi_one_fprintf(stderr,
			  "must have 1 <= band index <= num_bands (%d)\n",H.p);
	  return;
     }

     curfield = (scalar_complex *) mdata->fft_data;
     curfield_band = which_band;
     curfield_type = 'h';
     maxwell_compute_h_from_H(mdata, H, curfield, which_band - 1, 1);

     /* Divide by the cell volume so that the integral of |H|^2
        or of D*E is unity.  (From the eigensolver + FFT, they are
        initially normalized to sum to nx*ny*nz.) */

     {
	  int i, N;
	  double scale;
	  N = mdata->fft_output_size;

	  scale = 1.0 / sqrt(Vol);
	  for (i = 0; i < 3*N; ++i) {
	       curfield[i].re *= scale;
	       curfield[i].im *= scale;
	  }
     }
}

void get_efield_from_dfield(void)
{
     if (!curfield || curfield_type != 'd') {
	  mpi_one_fprintf(stderr, "get-dfield must be called before "
		  "get-efield-from-dfield!\n");
	  return;
     }
     CHECK(mdata, "unexpected NULL mdata");
     maxwell_compute_e_from_d(mdata, curfield, 1);
     curfield_type = 'e';
}

void get_efield(integer which_band)
{
     get_dfield(which_band);
     get_efield_from_dfield();
}

/* Extract the mean epsilon from the effective inverse dielectric tensor,
   which contains two eigenvalues that correspond to the mean epsilon,
   and one which corresponds to the harmonic mean. */
real mean_epsilon_from_matrix(const symmetric_matrix *eps_inv)
{
     real eps_eigs[3];
     maxwell_sym_matrix_eigs(eps_eigs, eps_inv);
     /* the harmonic mean should be the largest eigenvalue (smallest
	epsilon), so we'll ignore it and average the other two: */
     return 2.0 / (eps_eigs[0] + eps_eigs[1]);
}

/* get the dielectric function, and compute some statistics */
void get_epsilon(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *epsilon;
     real eps_mean = 0, eps_inv_mean = 0, eps_high = -1e20, eps_low = 1e20;
     int fill_count = 0;

     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-epsilon!\n");
	  return;
     }

     curfield = (scalar_complex *) mdata->fft_data;
     epsilon = (real *) curfield;
     curfield_band = 0;
     curfield_type = 'n';

     /* get epsilon.  Recall that we actually have an inverse
	dielectric tensor at each point; define an average index by
	the inverse of the average eigenvalue of the 1/eps tensor.
	i.e. 3/(trace 1/eps). */

     N = mdata->fft_output_size;
     last_dim = mdata->last_dim;
     last_dim_stored =
	  mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

     for (i = 0; i < N; ++i) {
          epsilon[i] = mean_epsilon_from_matrix(mdata->eps_inv + i);
	  if (epsilon[i] < eps_low)
	       eps_low = epsilon[i];
	  if (epsilon[i] > eps_high)
	       eps_high = epsilon[i];
	  eps_mean += epsilon[i];
	  eps_inv_mean += 1/epsilon[i];
	  if (epsilon[i] > 1.0001)
	       ++fill_count;
#ifndef SCALAR_COMPLEX
	  /* most points need to be counted twice, by rfftw output symmetry: */
	  {
	       int last_index;
#  ifdef HAVE_MPI
	       if (nz == 1) /* 2d calculation: 1st dim. is truncated one */
		    last_index = i / nx + local_y_start;
	       else
		    last_index = i % last_dim_stored;
#  else
	       last_index = i % last_dim_stored;
#  endif
	       if (last_index != 0 && 2*last_index != last_dim) {
		    eps_mean += epsilon[i];
		    eps_inv_mean += 1/epsilon[i];
		    if (epsilon[i] > 1.0001)
			 ++fill_count;
	       }
	  }
#endif
     }

     mpi_allreduce_1(&eps_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce_1(&eps_inv_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce_1(&eps_low, real, SCALAR_MPI_TYPE,
		     MPI_MIN, MPI_COMM_WORLD);
     mpi_allreduce_1(&eps_high, real, SCALAR_MPI_TYPE,
		     MPI_MAX, MPI_COMM_WORLD);
     mpi_allreduce_1(&fill_count, int, MPI_INT,
                   MPI_SUM, MPI_COMM_WORLD);
     N = mdata->nx * mdata->ny * mdata->nz;
     eps_mean /= N;
     eps_inv_mean = N/eps_inv_mean;

     mpi_one_printf("epsilon: %g-%g, mean %g, harm. mean %g, "
		    "%g%% > 1, %g%% \"fill\"\n",
		    eps_low, eps_high, eps_mean, eps_inv_mean,
		    (100.0 * fill_count) / N, 
		    eps_high == eps_low ? 100.0 :
		    100.0 * (eps_mean-eps_low) / (eps_high-eps_low));
}

/* get the specified component of the dielectric tensor,
   or the inverse tensor if inv != 0 */
void get_epsilon_tensor(int c1, int c2, int imag, int inv)
{
     int i, N;
     real *epsilon;
     int conj = 0, offset = 0;

     curfield_type = '-'; /* only used internally, for now */
     epsilon = (real *) mdata->fft_data;
     N = mdata->fft_output_size;

     switch (c1 * 3 + c2) {
	 case 0:
	      offset = offsetof(symmetric_matrix, m00);
	      break;
	 case 1:
	      offset = offsetof(symmetric_matrix, m01);
	      break;
	 case 2:
	      offset = offsetof(symmetric_matrix, m02);
	      break;
	 case 3:
	      offset = offsetof(symmetric_matrix, m01); /* = conj(m10) */
	      conj = imag;
	      break;
	 case 4:
	      offset = offsetof(symmetric_matrix, m11);
	      break;
	 case 5:
	      offset = offsetof(symmetric_matrix, m12);
	      break;
	 case 6:
	      offset = offsetof(symmetric_matrix, m02); /* = conj(m20) */
	      conj = imag;
	      break;
	 case 7:
	      offset = offsetof(symmetric_matrix, m12); /* = conj(m21) */
	      conj = imag;
	      break;
	 case 8:
	      offset = offsetof(symmetric_matrix, m22);
	      break;
     }

#ifdef WITH_HERMITIAN_EPSILON
     if (c1 != c2 && imag)
	  offset += offsetof(scalar_complex, im);
#endif

     for (i = 0; i < N; ++i) {
	  if (inv) {
	       epsilon[i] = 
		    *((real *) (((char *) &mdata->eps_inv[i]) + offset));
	  }
	  else {
	       symmetric_matrix eps;
	       maxwell_sym_matrix_invert(&eps, &mdata->eps_inv[i]);
	       epsilon[i] = *((real *) (((char *) &eps) + offset));
	  }
	  if (conj)
	       epsilon[i] = -epsilon[i];
     }     
}

/**************************************************************************/

/* internal function for compute_field_energy, below */
double compute_field_energy_internal(real comp_sum[6])
{
    int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real comp_sum2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
     real energy_sum = 0.0;
     real *energy_density = (real *) curfield;

     N = mdata->fft_output_size;
     last_dim = mdata->last_dim;
     last_dim_stored =
	  mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

     for (i = 0; i < N; ++i) {
	  scalar_complex field[3];
	  real
	       comp_sqr0,comp_sqr1,comp_sqr2,comp_sqr3,comp_sqr4,comp_sqr5;

	  /* energy is either |curfield|^2 or |curfield|^2 / epsilon,
	     depending upon whether it is H or D. */
	  if (curfield_type == 'd') 
	       assign_symmatrix_vector(field, mdata->eps_inv[i], curfield+3*i);
	  else {
	       field[0] =   curfield[3*i];
	       field[1] = curfield[3*i+1];
	       field[2] = curfield[3*i+2];
	  }

	  comp_sum2[0] += comp_sqr0 = field[0].re *   curfield[3*i].re;
	  comp_sum2[1] += comp_sqr1 = field[0].im *   curfield[3*i].im;
	  comp_sum2[2] += comp_sqr2 = field[1].re * curfield[3*i+1].re;
	  comp_sum2[3] += comp_sqr3 = field[1].im * curfield[3*i+1].im;
	  comp_sum2[4] += comp_sqr4 = field[2].re * curfield[3*i+2].re;
	  comp_sum2[5] += comp_sqr5 = field[2].im * curfield[3*i+2].im;

	  /* Note: here, we write to energy_density[i]; this is
	     safe, even though energy_density is aliased to curfield,
	     since energy_density[i] is guaranteed to come at or before
	     curfield[i] (which we are now done with). */

	  energy_sum += energy_density[i] = 
	       comp_sqr0+comp_sqr1+comp_sqr2+comp_sqr3+comp_sqr4+comp_sqr5;
#ifndef SCALAR_COMPLEX
	  /* most points need to be counted twice, by rfftw output symmetry: */
	  {
	       int last_index;
#  ifdef HAVE_MPI
	       if (nz == 1) /* 2d calculation: 1st dim. is truncated one */
		    last_index = i / nx + local_y_start;
	       else
		    last_index = i % last_dim_stored;
#  else
	       last_index = i % last_dim_stored;
#  endif
	       if (last_index != 0 && 2*last_index != last_dim) {
		    energy_sum += energy_density[i];
		    comp_sum2[0] += comp_sqr0;
		    comp_sum2[1] += comp_sqr1;
		    comp_sum2[2] += comp_sqr2;
		    comp_sum2[3] += comp_sqr3;
		    comp_sum2[4] += comp_sqr4;
		    comp_sum2[5] += comp_sqr5;
	       }
	  }
#endif
     }

     mpi_allreduce_1(&energy_sum, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce(comp_sum2, comp_sum, 6, real, SCALAR_MPI_TYPE,
                   MPI_SUM, MPI_COMM_WORLD);

     /* remember that we now have energy density; denoted by capital D/H */
     curfield_type = toupper(curfield_type);

     return energy_sum;
}

/* Replace curfield (either d or h) with the scalar energy density function,
   normalized to one.  While we're at it, compute some statistics about
   the relative strength of different field components.  Also return
   the integral of the energy density, which should be unity. */
number_list compute_field_energy(void)
{
     int i;
     real energy_sum, comp_sum[6];
     number_list retval = { 0, 0 };

     if (!curfield || !strchr("dh", curfield_type)) {
	  mpi_one_fprintf(stderr, "The D or H field must be loaded first.\n");
	  return retval;
     }

     energy_sum = compute_field_energy_internal(comp_sum);

     mpi_one_printf("%c-energy-components:, %d, %d",
	    curfield_type, kpoint_index, curfield_band);
     for (i = 0; i < 6; ++i) {
	  comp_sum[i] /= (energy_sum == 0 ? 1 : energy_sum);
	  if (i % 2 == 1)
	       mpi_one_printf(", %g", comp_sum[i] + comp_sum[i-1]);
     }
     mpi_one_printf("\n");

      /* The return value is a list of 7 items: the total energy,
	followed by the 6 elements of the comp_sum array (the fraction
	of the energy in the real/imag. parts of each field component). */

     retval.num_items = 7;
     CHK_MALLOC(retval.items, number, retval.num_items);

     retval.items[0] = energy_sum * Vol / H.N;

     for (i = 0; i < 6; ++i)
	  retval.items[i+1] = comp_sum[i];

     return retval;
}

/* compute |F|^2 for the current field, whether or not this is an
   energy density */
void compute_field_squared(void)
{
     real comp_sum[6]; /* unused */

     if (!curfield || !strchr("dhecv", curfield_type)) {
          mpi_one_fprintf(stderr, "A vector field must be loaded first.\n");
     }

     curfield_type = 'c';  /* force it to just square the field */
     compute_field_energy_internal(comp_sum);
     curfield_type = 'R'; /* generic real scalar field */
}

/**************************************************************************/

/* replace the current field with its scalar divergence; only works
   for Bloch fields */

void compute_field_divergence(void)
{
     int i, j, N;
     scalar *field = (scalar *) curfield;
     scalar *field2 = mdata->fft_data == mdata->fft_data2 ? field : (field == mdata->fft_data ? mdata->fft_data2 : mdata->fft_data);
     real scale;

     if (!curfield || !strchr("dhec", curfield_type)) {
          mpi_one_fprintf(stderr, "A Bloch-periodic field must be loaded.\n");
          return;
     }

     /* convert back to Fourier space */
     maxwell_compute_fft(-1, mdata, field, field2, 3, 3, 1);

     /* compute (k+G) dot field */
     for (i = 0; i < mdata->other_dims; ++i)
          for (j = 0; j < mdata->last_dim; ++j) {
               int ij = i * mdata->last_dim_size + j;
	       k_data cur_k = mdata->k_plus_G[ij];
	       /* k+G = |k+G| (m x n) */
	       real kx = cur_k.kmag * (cur_k.my*cur_k.nz-cur_k.mz*cur_k.ny);
	       real ky = cur_k.kmag * (cur_k.mz*cur_k.nx-cur_k.mx*cur_k.nz);
	       real kz = cur_k.kmag * (cur_k.mx*cur_k.ny-cur_k.my*cur_k.nz);
	       ASSIGN_SCALAR(field2[ij],
			     SCALAR_RE(field2[3*ij+0]) * kx +
			     SCALAR_RE(field2[3*ij+1]) * ky +
			     SCALAR_RE(field2[3*ij+2]) * kz,
			     SCALAR_IM(field2[3*ij+0]) * kx +
			     SCALAR_IM(field2[3*ij+1]) * ky +
			     SCALAR_IM(field2[3*ij+2]) * kz);
	  }

     /* convert scalar field back to position space */
     maxwell_compute_fft(+1, mdata, field2, field, 1, 1, 1);

     /* multiply by i (from divergence) and normalization (from FFT)
        and 2*pi (from k+G) */
     scale = TWOPI / H.N;
     N = mdata->fft_output_size;
     for (i = 0; i < N; ++i) {
	  CASSIGN_SCALAR(curfield[i],
			 -CSCALAR_IM(curfield[i]) * scale,
			 CSCALAR_RE(curfield[i]) * scale);
     }

     curfield_type = 'C'; /* complex (Bloch) scalar field */
}

/**************************************************************************/

/* Fix the phase of the current field (e/h/d) to a canonical value.
   Also changes the phase of the corresponding eigenvector by the
   same amount, so that future calculations will have a consistent
   phase.

   The following procedure is used, derived from a suggestion by Doug
   Allan of Corning: First, choose the phase to maximize the sum of
   the squares of the real parts of the components.  This doesn't fix
   the overall sign, though.  That is done (after incorporating the
   above phase) by: (1) find the largest absolute value of the real
   part, (2) find the point with the greatest spatial array index that
   has |real part| at least half of the largest value, and (3) make
   that point positive.

   In the case of inversion symmetry, on the other hand, the overall phase
   is already fixed, to within a sign, by the choice to make the Fourier
   transform purely real.  So, in that case we simply pick a sign, in
   a manner similar to (2) and (3) above. */
void fix_field_phase(void)
{
     int i, N;
     real sq_sum2[2] = {0,0}, sq_sum[2], maxabs = 0.0;
     int maxabs_index = 0, maxabs_sign = 1;
     double theta;
     scalar phase;

     if (!curfield || !strchr("dhecv", curfield_type)) {
          mpi_one_fprintf(stderr, "The D/H/E field must be loaded first.\n");
          return;
     }
     N = mdata->fft_output_size * 3;

#ifdef SCALAR_COMPLEX
     /* Compute the phase that maximizes the sum of the squares of
	the real parts of the components.  Equivalently, maximize
	the real part of the sum of the squares. */
     for (i = 0; i < N; ++i) {
	  real a,b;
	  a = curfield[i].re; b = curfield[i].im;
	  sq_sum2[0] += a*a - b*b;
	  sq_sum2[1] += 2*a*b;
     }
     mpi_allreduce(sq_sum2, sq_sum, 2, real, SCALAR_MPI_TYPE,
                   MPI_SUM, MPI_COMM_WORLD);
     /* compute the phase = exp(i*theta) maximizing the real part of
	the sum of the squares.  i.e., maximize:
	    cos(2*theta)*sq_sum[0] - sin(2*theta)*sq_sum[1] */
     theta = 0.5 * atan2(-sq_sum[1], sq_sum[0]);
     phase.re = cos(theta);
     phase.im = sin(theta);
#else /* ! SCALAR_COMPLEX */
     phase = 1;
#endif /* ! SCALAR_COMPLEX */

     /* Next, fix the overall sign.  We do this by first computing the
	maximum |real part| of the jmax component (after multiplying
	by phase), and then finding the last spatial index at which
	|real part| is at least half of this value.  The sign is then
	chosen to make the real part positive at that point. 

        (Note that we can't just make the point of maximum |real part|
         positive, as that would be ambiguous in the common case of an
         oscillating field within the unit cell.)

        In the case of inversion symmetry (!SCALAR_COMPLEX), we work with
        (real part - imag part) instead of (real part), to insure that we
        have something that is nonzero somewhere. */

     for (i = 0; i < N; ++i) {
#ifdef SCALAR_COMPLEX
	  real r = fabs(curfield[i].re * phase.re - curfield[i].im * phase.im);
#else
	  real r = fabs(curfield[i].re - curfield[i].im);
#endif
	  if (r > maxabs)
	       maxabs = r;
     }
     mpi_allreduce_1(&maxabs, real, SCALAR_MPI_TYPE,
		     MPI_MAX, MPI_COMM_WORLD);
     for (i = N - 1; i >= 0; --i) {
#ifdef SCALAR_COMPLEX
	  real r = curfield[i].re * phase.re - curfield[i].im * phase.im;
#else
	  real r = curfield[i].re - curfield[i].im;
#endif
	  if (fabs(r) >= 0.5 * maxabs) {
	       maxabs_index = i;
	       maxabs_sign = r < 0 ? -1 : 1;
	       break;
	  }
     }
     if (i >= 0)  /* convert index to global index in distributed array: */
	  maxabs_index += mdata->local_y_start * mdata->nx * mdata->nz;
     {
	  /* compute maximum index and corresponding sign over all the 
	     processors, using the MPI_MAXLOC reduction operation: */
	  struct twoint_struct {int i; int s;} x;
	  x.i = maxabs_index; x.s = maxabs_sign;
	  mpi_allreduce_1(&x, struct twoint_struct, MPI_2INT,
			  MPI_MAXLOC, MPI_COMM_WORLD);
	  maxabs_index = x.i; maxabs_sign = x.s;
     }
     ASSIGN_SCALAR(phase,
		   SCALAR_RE(phase)*maxabs_sign, SCALAR_IM(phase)*maxabs_sign);

     mpi_one_printf("Fixing %c-field (band %d) phase by %g + %gi; "
		    "max ampl. = %g\n", curfield_type, curfield_band,
		    SCALAR_RE(phase), SCALAR_IM(phase), maxabs);

     /* Now, multiply everything by this phase, *including* the
	stored "raw" eigenvector in H, so that any future fields
	that we compute will have a consistent phase: */
     for (i = 0; i < N; ++i) {
	  real a,b;
	  a = curfield[i].re; b = curfield[i].im;
	  curfield[i].re = a*SCALAR_RE(phase) - b*SCALAR_IM(phase);
	  curfield[i].im = a*SCALAR_IM(phase) + b*SCALAR_RE(phase);
     }
     for (i = 0; i < H.n; ++i) {
          ASSIGN_MULT(H.data[i*H.p + curfield_band - 1], 
		      H.data[i*H.p + curfield_band - 1], phase);
     }
}

/**************************************************************************/

/* Functions to return epsilon, fields, energies, etcetera, at a specified
   point, linearly interpolating if necessary. */

static real get_val(int ix, int iy, int iz,
		    int nx, int ny, int nz, int last_dim_size,
		    real *data, int stride, int conjugate)
{
#ifndef SCALAR_COMPLEX
     {
	  int nlast = last_dim_size / 2;
	  if ((nz > 1 ? iz : (ny > 1 ? iy : ix)) >= nlast) {
	       ix = ix ? nx - ix : ix;
	       iy = iy ? ny - iy : iy;
	       iz = iz ? nz - iz : iz;
	       conjugate = conjugate ? 1 : 0;
	  }
	  else
	       conjugate = 0;
	  if (nz > 1) nz = nlast; else if (ny > 1) ny = nlast; else nx = nlast;
     }
#else
     conjugate = 0;
#endif

#ifdef HAVE_MPI
     CHECK(0, "get-*-point not yet implemented for MPI!");
#else
     if (conjugate)
	  return -data[(((ix * ny) + iy) * nz + iz) * stride];
     else
	  return data[(((ix * ny) + iy) * nz + iz) * stride];
#endif
}

static real interp_val(vector3 p, int nx, int ny, int nz, int last_dim_size,
		       real *data, int stride, int conjugate)
{
     double ipart;
     real rx, ry, rz, dx, dy, dz;
     int x, y, z, x2, y2, z2;

     rx = modf(p.x/geometry_lattice.size.x + 0.5, &ipart); if (rx < 0) rx += 1;
     ry = modf(p.y/geometry_lattice.size.y + 0.5, &ipart); if (ry < 0) ry += 1;
     rz = modf(p.z/geometry_lattice.size.z + 0.5, &ipart); if (rz < 0) rz += 1;

     /* get the point corresponding to r in the grid: */
     x = rx * nx;
     y = ry * ny;
     z = rz * nz;

     /* get the difference between (x,y,z) and the actual point */
     dx = rx * nx - x;
     dy = ry * ny - y;
     dz = rz * nz - z;

     /* get the other closest point in the grid, with periodic boundaries: */
     x2 = (nx + (dx >= 0.0 ? x + 1 : x - 1)) % nx;
     y2 = (ny + (dy >= 0.0 ? y + 1 : y - 1)) % ny;
     z2 = (nz + (dz >= 0.0 ? z + 1 : z - 1)) % nz;

     /* take abs(d{xyz}) to get weights for {xyz} and {xyz}2: */
     dx = fabs(dx);
     dy = fabs(dy);
     dz = fabs(dz);

#define D(x,y,z) (get_val(x,y,z,nx,ny,nz,last_dim_size, data,stride,conjugate))

     return(((D(x,y,z)*(1.0-dx) + D(x2,y,z)*dx) * (1.0-dy) +
             (D(x,y2,z)*(1.0-dx) + D(x2,y2,z)*dx) * dy) * (1.0-dz) +
            ((D(x,y,z2)*(1.0-dx) + D(x2,y,z2)*dx) * (1.0-dy) +
             (D(x,y2,z2)*(1.0-dx) + D(x2,y2,z2)*dx) * dy) * dz);

#undef D
}

static scalar_complex interp_cval(vector3 p, 
				  int nx, int ny, int nz, int last_dim_size,
				  real *data, int stride)
{
     scalar_complex cval;
     cval.re = interp_val(p, nx,ny,nz,last_dim_size, data, stride, 0);
     cval.im = interp_val(p, nx,ny,nz,last_dim_size,data + 1, stride, 1);
     return cval;
}

#define f_interp_val(p,f,data,stride,conj) interp_val(p,f->nx,f->ny,f->nz,f->last_dim_size,data,stride,conj)
#define f_interp_cval(p,f,data,stride) interp_cval(p,f->nx,f->ny,f->nz,f->last_dim_size,data,stride)

static symmetric_matrix interp_eps_inv(vector3 p)
{
     int stride = sizeof(symmetric_matrix) / sizeof(real);
     symmetric_matrix eps_inv;

     eps_inv.m00 = f_interp_val(p, mdata, &mdata->eps_inv->m00, stride, 0);
     eps_inv.m11 = f_interp_val(p, mdata, &mdata->eps_inv->m11, stride, 0);
     eps_inv.m22 = f_interp_val(p, mdata, &mdata->eps_inv->m22, stride, 0);
#ifdef WITH_HERMITIAN_EPSILON
     eps_inv.m01 = f_interp_cval(p, mdata, &mdata->eps_inv->m01.re, stride);
     eps_inv.m02 = f_interp_cval(p, mdata, &mdata->eps_inv->m02.re, stride);
     eps_inv.m12 = f_interp_cval(p, mdata, &mdata->eps_inv->m12.re, stride);
#else
     eps_inv.m01 = f_interp_val(p, mdata, &mdata->eps_inv->m01, stride, 0);
     eps_inv.m02 = f_interp_val(p, mdata, &mdata->eps_inv->m02, stride, 0);
     eps_inv.m12 = f_interp_val(p, mdata, &mdata->eps_inv->m12, stride, 0);
#endif
     return eps_inv;
}

number get_epsilon_point(vector3 p)
{
     symmetric_matrix eps_inv;
     eps_inv = interp_eps_inv(p);
     return mean_epsilon_from_matrix(&eps_inv);
}

cmatrix3x3 get_epsilon_inverse_tensor_point(vector3 p)
{
     symmetric_matrix eps_inv;
     eps_inv = interp_eps_inv(p);

#ifdef WITH_HERMITIAN_EPSILON
     return make_hermitian_cmatrix3x3(eps_inv.m00,eps_inv.m11,eps_inv.m22,
				      cscalar2cnumber(eps_inv.m01),
				      cscalar2cnumber(eps_inv.m02),
				      cscalar2cnumber(eps_inv.m12));
#else
     return make_hermitian_cmatrix3x3(eps_inv.m00,eps_inv.m11,eps_inv.m22,
				      make_cnumber(eps_inv.m01,0),
				      make_cnumber(eps_inv.m02,0),
				      make_cnumber(eps_inv.m12,0));
#endif
}

number get_energy_point(vector3 p)
{
     CHECK(curfield && strchr("DHR", curfield_type),
	   "compute-field-energy must be called before get-energy-point");
     return f_interp_val(p, mdata, (real *) curfield, 1, 0);
}

cvector3 get_bloch_field_point(vector3 p)
{
     scalar_complex field[3];
     cvector3 F;

     CHECK(curfield && strchr("dhecv", curfield_type),
	   "field must be must be loaded before get-*field*-point");
     field[0] = f_interp_cval(p, mdata, &curfield[0].re, 6);
     field[1] = f_interp_cval(p, mdata, &curfield[1].re, 6);
     field[2] = f_interp_cval(p, mdata, &curfield[2].re, 6);
     F.x = cscalar2cnumber(field[0]);
     F.y = cscalar2cnumber(field[1]);
     F.z = cscalar2cnumber(field[2]);
     return F;
}

cvector3 get_field_point(vector3 p)
{
     scalar_complex field[3], phase;
     cvector3 F;

     CHECK(curfield && strchr("dhecv", curfield_type),
	   "field must be must be loaded before get-*field*-point");
     field[0] = f_interp_cval(p, mdata, &curfield[0].re, 6);
     field[1] = f_interp_cval(p, mdata, &curfield[1].re, 6);
     field[2] = f_interp_cval(p, mdata, &curfield[2].re, 6);

     if (curfield_type != 'v') {
	  double phase_phi = TWOPI * 
	       (cur_kvector.x * (p.x/geometry_lattice.size.x) +
		cur_kvector.y * (p.y/geometry_lattice.size.y) +
		cur_kvector.z * (p.z/geometry_lattice.size.z));
	  CASSIGN_SCALAR(phase, cos(phase_phi), sin(phase_phi));
	  CASSIGN_MULT(field[0], field[0], phase);
	  CASSIGN_MULT(field[1], field[1], phase);
	  CASSIGN_MULT(field[2], field[2], phase);
     }

     F.x = cscalar2cnumber(field[0]);
     F.y = cscalar2cnumber(field[1]);
     F.z = cscalar2cnumber(field[2]);
     return F;
}

cnumber get_bloch_cscalar_point(vector3 p)
{
     CHECK(curfield && strchr("C", curfield_type),
	   "cscalar must be must be loaded before get-*cscalar*-point");
     
     return cscalar2cnumber(f_interp_cval(p, mdata, &curfield[0].re, 2));
}

cnumber get_cscalar_point(vector3 p)
{
     scalar_complex s;

     CHECK(curfield && strchr("C", curfield_type),
	   "cscalar must be must be loaded before get-*cscalar*-point");
     
     s = f_interp_cval(p, mdata, &curfield[0].re, 2);

     if (curfield_type == 'C') {
	  scalar_complex phase;
	  double phase_phi = TWOPI * 
	       (cur_kvector.x * (p.x/geometry_lattice.size.x) +
		cur_kvector.y * (p.y/geometry_lattice.size.y) +
		cur_kvector.z * (p.z/geometry_lattice.size.z));
	  CASSIGN_SCALAR(phase, cos(phase_phi), sin(phase_phi));
	  CASSIGN_MULT(s, s, phase);
     }

     return cscalar2cnumber(s);
}

number rscalar_field_get_point(SCM fo, vector3 p)
{
     field_smob *f = assert_field_smob(fo);
     CHECK(f->type == RSCALAR_FIELD_SMOB, 
	   "invalid argument to rscalar-field-get-point");
     return f_interp_val(p, f, f->f.rs, 1, 0);
}

cvector3 cvector_field_get_point_bloch(SCM fo, vector3 p)
{
     scalar_complex field[3];
     cvector3 F;
     field_smob *f = assert_field_smob(fo);
     CHECK(f->type == CVECTOR_FIELD_SMOB, 
	   "invalid argument to cvector-field-get-point");
     field[0] = f_interp_cval(p, f, &f->f.cv[0].re, 6);
     field[1] = f_interp_cval(p, f, &f->f.cv[1].re, 6);
     field[2] = f_interp_cval(p, f, &f->f.cv[2].re, 6);
     F.x = cscalar2cnumber(field[0]);
     F.y = cscalar2cnumber(field[1]);
     F.z = cscalar2cnumber(field[2]);
     return F;
}

cvector3 cvector_field_get_point(SCM fo, vector3 p)
{
     scalar_complex field[3];
     cvector3 F;
     field_smob *f = assert_field_smob(fo);
     CHECK(f->type == CVECTOR_FIELD_SMOB, 
	   "invalid argument to cvector-field-get-point");

     field[0] = f_interp_cval(p, f, &f->f.cv[0].re, 6);
     field[1] = f_interp_cval(p, f, &f->f.cv[1].re, 6);
     field[2] = f_interp_cval(p, f, &f->f.cv[2].re, 6);
     
     if (f->type_char != 'v') { /* v fields have no kvector */
	  scalar_complex phase;
	  double phase_phi = TWOPI * 
	       (cur_kvector.x * (p.x/geometry_lattice.size.x) +
		cur_kvector.y * (p.y/geometry_lattice.size.y) +
		cur_kvector.z * (p.z/geometry_lattice.size.z));
	  CASSIGN_SCALAR(phase, cos(phase_phi), sin(phase_phi));
	  CASSIGN_MULT(field[0], field[0], phase);
	  CASSIGN_MULT(field[1], field[1], phase);
	  CASSIGN_MULT(field[2], field[2], phase);
     }

     F.x = cscalar2cnumber(field[0]);
     F.y = cscalar2cnumber(field[1]);
     F.z = cscalar2cnumber(field[2]);
     return F;
}

cnumber cscalar_field_get_point_bloch(SCM fo, vector3 p)
{
     field_smob *f = assert_field_smob(fo);
     CHECK(f->type == CSCALAR_FIELD_SMOB, 
	   "invalid argument to cscalar-field-get-point-bloch");
     return cscalar2cnumber(f_interp_cval(p, f, &f->f.cv[0].re, 2));
}

cnumber cscalar_field_get_point(SCM fo, vector3 p)
{
     scalar_complex s;
     field_smob *f = assert_field_smob(fo);
     CHECK(f->type == CSCALAR_FIELD_SMOB, 
	   "invalid argument to cscalar-field-get-point");

     s = f_interp_cval(p, f, &f->f.cv[0].re, 2);
     
     if (f->type_char == 'C') { /* have kvector */
	  scalar_complex phase;
	  double phase_phi = TWOPI * 
	       (cur_kvector.x * (p.x/geometry_lattice.size.x) +
		cur_kvector.y * (p.y/geometry_lattice.size.y) +
		cur_kvector.z * (p.z/geometry_lattice.size.z));
	  CASSIGN_SCALAR(phase, cos(phase_phi), sin(phase_phi));
	  CASSIGN_MULT(s, s, phase);
     }

     return cscalar2cnumber(s);
}

/**************************************************************************/

/* compute the fraction of the field energy that is located in the
   given range of dielectric constants: */
number compute_energy_in_dielectric(number eps_low, number eps_high)
{
     int N, i, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *energy = (real *) curfield;
     real epsilon, energy_sum = 0.0;

     if (!curfield || !strchr("DHR", curfield_type)) {
          mpi_one_fprintf(stderr, "The D or H energy density must be loaded first.\n");
          return 0.0;
     }

     N = mdata->fft_output_size;
     last_dim = mdata->last_dim;
     last_dim_stored =
	  mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

     for (i = 0; i < N; ++i) {
	  epsilon = mean_epsilon_from_matrix(mdata->eps_inv +i);
	  if (epsilon >= eps_low && epsilon <= eps_high) {
	       energy_sum += energy[i];
#ifndef SCALAR_COMPLEX
	       /* most points are counted twice, by rfftw output symmetry: */
	       {
		    int last_index;
#  ifdef HAVE_MPI
		    if (nz == 1) /* 2d: 1st dim. is truncated one */
			 last_index = i / nx + local_y_start;
		    else
			 last_index = i % last_dim_stored;
#  else
		    last_index = i % last_dim_stored;
#  endif
		    if (last_index != 0 && 2*last_index != last_dim)
			 energy_sum += energy[i];
	       }
#endif
	  }
     }
     mpi_allreduce_1(&energy_sum, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     energy_sum *= Vol / H.N;
     return energy_sum;
}

/**************************************************************************/

/* Prepend the prefix to the fname, and (if parity_suffix is true)
   append a parity specifier (if any) (e.g. ".te"), returning a new
   string, which should be deallocated with free().  fname or prefix
   may be NULL, in which case they are treated as the empty string. */
static char *fix_fname(const char *fname, const char *prefix,
		       maxwell_data *d, int parity_suffix)
{
     char *s;
     CHK_MALLOC(s, char,
		(fname ? strlen(fname) : 0) + 
		(prefix ? strlen(prefix) : 0) + 20);
     strcpy(s, prefix ? prefix : "");
     strcat(s, fname ? fname : "");
     if (parity_suffix && d->parity != NO_PARITY) {
	  /* assumes parity suffix is less than 20 characters;
	     currently it is less than 12 */
	  strcat(s, ".");
	  strcat(s, parity_string(d));
     }
     return s;
}

static void output_scalarfield(real *vals,
			       const int dims[3],
			       const int local_dims[3],
			       const int start[3],
			       matrixio_id file_id,
			       const char *dataname,
			       int last_dim_index,
			       int last_dim_start, int last_dim_size,
			       int first_dim_start, int first_dim_size,
			       int write_start0_special)
{
     matrixio_id data_id = {-1, 1};

     fieldio_write_real_vals(vals, 3, dims, 
			     local_dims, start, file_id, 0, 
			     dataname, &data_id);
     
#ifndef SCALAR_COMPLEX
     {
	  int start_new[3], local_dims_new[3];

	  start_new[0] = start[0];
	  start_new[1] = start[1];
	  start_new[2] = start[2];
	  local_dims_new[0] = local_dims[0];
	  local_dims_new[1] = local_dims[1];
	  local_dims_new[2] = local_dims[2];

	  maxwell_scalarfield_otherhalf(mdata, vals);
	  start_new[last_dim_index] = last_dim_start;
	  local_dims_new[last_dim_index] = last_dim_size;
	  start_new[0] = first_dim_start;
	  local_dims_new[0] = first_dim_size;
	  if (write_start0_special) {
	       /* The conjugated array half may be discontiguous.
		  First, write the part not containing start_new[0], and
		  then write the start_new[0] slab. */
	       fieldio_write_real_vals(vals +
				       local_dims_new[1] * local_dims_new[2],
				       3, dims, local_dims_new, start_new, 
				       file_id, 1, dataname, &data_id);
	       local_dims_new[0] = 1;
	       start_new[0] = 0;
	       fieldio_write_real_vals(vals, 3, dims, 
				       local_dims_new, start_new,
				       file_id, 1, dataname, &data_id);
	  }
	  else {
	       fieldio_write_real_vals(vals, 3, dims, 
				       local_dims_new, start_new, 
				       file_id, 1, dataname, &data_id);
	  }
     }
#endif

     if (data_id.id >= 0)
	  matrixio_close_dataset(data_id);
}

/* given the field in curfield, store it to HDF (or whatever) using
   the matrixio (fieldio) routines.  Allow the component to be specified
   (which_component 0/1/2 = x/y/z, -1 = all) for vector fields. 
   Also allow the user to specify a prefix string for the filename. */
void output_field_to_file(integer which_component, string filename_prefix)
{
     char fname[100], *fname2, description[100];
     int dims[3], local_dims[3], start[3] = {0,0,0};
     matrixio_id file_id = {-1,1};
     int attr_dims[2] = {3, 3};
     real output_k[3]; /* kvector in reciprocal lattice basis */
     real output_R[3][3];

     /* where to put "otherhalf" block of output, only used for real scalars */
     int last_dim_index = 0;
     int last_dim_start = 0, last_dim_size = 0;
     int first_dim_start = 0, first_dim_size = 0;
     int write_start0_special = 0;

     if (!curfield) {
	  mpi_one_fprintf(stderr, 
		  "fields, energy dens., or epsilon must be loaded first.\n");
	  return;
     }
     
#ifdef HAVE_MPI
     /* The first two dimensions (x and y) of the position-space fields
	are transposed when we use MPI, so we need to transpose everything. */
     dims[0] = mdata->ny;
     local_dims[1] = dims[1] = mdata->nx;
     local_dims[2] = dims[2] = mdata->nz;
     local_dims[0] = mdata->local_ny;
     start[0] = mdata->local_y_start;
#  ifndef SCALAR_COMPLEX
     /* Ugh, hairy.  See also maxwell_vectorfield_otherhalf. */
     if (dims[2] == 1) {
	  last_dim_index = 0;
	  first_dim_size = local_dims[0];
	  first_dim_start = dims[0] - (start[0] + local_dims[0] - 1);

	  if (start[0] == 0)
	       --first_dim_size; /* DC frequency is not in other half */
	  if (start[0] + local_dims[0] == mdata->last_dim_size / 2 &&
	      dims[0] % 2 == 0) {
	       --first_dim_size; /* Nyquist frequency is not in other half */
	       ++first_dim_start;
	  }

	  last_dim_start = first_dim_start;
	  last_dim_size = first_dim_size;
     }
     else {
	  last_dim_index = 2;
	  local_dims[last_dim_index] = mdata->last_dim_size / 2;
	  if (start[0] == 0) {
	       first_dim_size = local_dims[0] - 1;
	       first_dim_start = dims[0] - first_dim_size;
	       write_start0_special = 1;
	  }
	  else {
	       first_dim_start = dims[0] - (start[0] + local_dims[0] - 1);
	       first_dim_size = local_dims[0];
	  }
	  last_dim_start = local_dims[last_dim_index];
	  last_dim_size = dims[last_dim_index] - local_dims[last_dim_index];
     }
#  endif /* ! SCALAR_COMPLEX */
     output_k[0] = R[1][0]*mdata->current_k[0] + R[1][1]*mdata->current_k[1]
	  + R[1][2]*mdata->current_k[2];
     output_k[1] = R[0][0]*mdata->current_k[0] + R[0][1]*mdata->current_k[1]
	  + R[0][2]*mdata->current_k[2];
     output_k[2] = R[2][0]*mdata->current_k[0] + R[2][1]*mdata->current_k[1]
	  + R[2][2]*mdata->current_k[2];
     output_R[0][0]=R[1][0]; output_R[0][1]=R[1][1]; output_R[0][2]=R[1][2];
     output_R[1][0]=R[0][0]; output_R[1][1]=R[0][1]; output_R[1][2]=R[0][2];
     output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
#else /* ! HAVE_MPI */
     dims[0] = mdata->nx;
     local_dims[1] = dims[1] = mdata->ny;
     local_dims[2] = dims[2] = mdata->nz;
     local_dims[0] = mdata->local_nx;
#  ifndef SCALAR_COMPLEX
     last_dim_index = dims[2] == 1 ? (dims[1] == 1 ? 0 : 1) : 2;
     local_dims[last_dim_index] = mdata->last_dim_size / 2;
     last_dim_start = local_dims[last_dim_index];
     last_dim_size = dims[last_dim_index] - local_dims[last_dim_index];
     first_dim_start = last_dim_index ? 0 : last_dim_start;
     first_dim_size = last_dim_index ? local_dims[0] : last_dim_size;
#  endif
     start[0] = mdata->local_x_start;
     output_k[0] = R[0][0]*mdata->current_k[0] + R[0][1]*mdata->current_k[1]
	  + R[0][2]*mdata->current_k[2];
     output_k[1] = R[1][0]*mdata->current_k[0] + R[1][1]*mdata->current_k[1]
	  + R[1][2]*mdata->current_k[2];
     output_k[2] = R[2][0]*mdata->current_k[0] + R[2][1]*mdata->current_k[1]
	  + R[2][2]*mdata->current_k[2];
     output_R[0][0]=R[0][0]; output_R[0][1]=R[0][1]; output_R[0][2]=R[0][2];
     output_R[1][0]=R[1][0]; output_R[1][1]=R[1][1]; output_R[1][2]=R[1][2];
     output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
#endif /* ! HAVE_MPI */

     if (strchr("Rv", curfield_type)) /* generic scalar/vector field */
	  output_k[0] = output_k[1] = output_k[2] = 0.0; /* don't know k */
     
     if (strchr("dhecv", curfield_type)) { /* outputting vector field */
	  matrixio_id data_id[6] = {{-1,1},{-1,1},{-1,1},{-1,1},{-1,1},{-1,1}};
	  int i;

	  sprintf(fname, "%c.k%02d.b%02d",
		  curfield_type, kpoint_index, curfield_band);
	  if (which_component >= 0) {
	       char comp_str[] = ".x";
	       comp_str[1] = 'x' + which_component;
	       strcat(fname, comp_str);
	  }
	  sprintf(description, "%c field, kpoint %d, band %d, freq=%g",
		  curfield_type, kpoint_index, curfield_band, 
		  freqs.items[curfield_band - 1]);
	  fname2 = fix_fname(fname, filename_prefix, mdata, 1);
	  mpi_one_printf("Outputting fields to %s...\n", fname2);
	  file_id = matrixio_create(fname2);
	  free(fname2);
	  fieldio_write_complex_field(curfield, 3, dims, local_dims, start,
				      which_component, 3, output_k,
				      file_id, 0, data_id);

#ifndef SCALAR_COMPLEX
	  /* Here's where it gets hairy. */
	  maxwell_vectorfield_otherhalf(mdata, curfield,
					output_k[0], output_k[1], output_k[2]);
	  start[last_dim_index] = last_dim_start;
	  local_dims[last_dim_index] = last_dim_size;
	  start[0] = first_dim_start;
	  local_dims[0] = first_dim_size;
	  if (write_start0_special) {
	       /* The conjugated array half may be discontiguous.
		  First, write the part not containing start[0], and
		  then write the start[0] slab. */
	       fieldio_write_complex_field(curfield + 
					   3 * local_dims[1] * local_dims[2],
					   3, dims, local_dims, start,
					   which_component, 3, NULL, 
					   file_id, 1, data_id);
	       local_dims[0] = 1;
	       start[0] = 0;
	       fieldio_write_complex_field(curfield, 3, dims,local_dims,start,
					   which_component, 3, NULL,
					   file_id, 1, data_id);
	  }
	  else {
	       fieldio_write_complex_field(curfield, 3, dims,local_dims,start,
					   which_component, 3, NULL,
					   file_id, 1, data_id);
	  }
#endif

	  for (i = 0; i < 6; ++i)
	       if (data_id[i].id >= 0)
		    matrixio_close_dataset(data_id[i]);
	  matrixio_write_data_attr(file_id, "Bloch wavevector",
				   output_k, 1, attr_dims);
     }
     else if (strchr("C", curfield_type)) { /* outputting cmplx scalar field */
	  matrixio_id data_id[2] = {{-1,1},{-1,1}};
	  int i;

	  sprintf(fname, "%c.k%02d.b%02d",
		  curfield_type, kpoint_index, curfield_band);
	  sprintf(description, "%c field, kpoint %d, band %d, freq=%g",
		  curfield_type, kpoint_index, curfield_band, 
		  freqs.items[curfield_band - 1]);
	  fname2 = fix_fname(fname, filename_prefix, mdata, 1);
	  mpi_one_printf("Outputting complex scalar field to %s...\n", fname2);
	  file_id = matrixio_create(fname2);
	  free(fname2);
	  fieldio_write_complex_field(curfield, 3, dims, local_dims, start,
				      which_component, 1, output_k,
				      file_id, 0, data_id);

#ifndef SCALAR_COMPLEX
	  /* Here's where it gets hairy. */
	  maxwell_cscalarfield_otherhalf(mdata, curfield,
					output_k[0], output_k[1], output_k[2]);
	  start[last_dim_index] = last_dim_start;
	  local_dims[last_dim_index] = last_dim_size;
	  start[0] = first_dim_start;
	  local_dims[0] = first_dim_size;
	  if (write_start0_special) {
	       /* The conjugated array half may be discontiguous.
		  First, write the part not containing start[0], and
		  then write the start[0] slab. */
	       fieldio_write_complex_field(curfield + 
					   local_dims[1] * local_dims[2],
					   3, dims, local_dims, start,
					   which_component, 1, NULL, 
					   file_id, 1, data_id);
	       local_dims[0] = 1;
	       start[0] = 0;
	       fieldio_write_complex_field(curfield, 3, dims,local_dims,start,
					   which_component, 1, NULL,
					   file_id, 1, data_id);
	  }
	  else {
	       fieldio_write_complex_field(curfield, 3, dims,local_dims,start,
					   which_component, 1, NULL,
					   file_id, 1, data_id);
	  }
#endif

	  for (i = 0; i < 2; ++i)
	       if (data_id[i].id >= 0)
		    matrixio_close_dataset(data_id[i]);
	  matrixio_write_data_attr(file_id, "Bloch wavevector",
				   output_k, 1, attr_dims);
     }
     else if (strchr("DHnR", curfield_type)) { /* scalar field */
	  if (curfield_type == 'n') {
	       sprintf(fname, "epsilon");
	       sprintf(description, "dielectric function, epsilon");
	  }
	  else {
	       sprintf(fname, "%cpwr.k%02d.b%02d",
		       tolower(curfield_type), kpoint_index, curfield_band);
	       sprintf(description,
		       "%c field energy density, kpoint %d, band %d, freq=%g",
		       curfield_type, kpoint_index, curfield_band, 
		       freqs.items[curfield_band - 1]);
	  }
	  fname2 = fix_fname(fname, filename_prefix, mdata, 
			     /* no parity suffix for epsilon: */
			     curfield_type != 'n');
	  mpi_one_printf("Outputting %s...\n", fname2);
	  file_id = matrixio_create(fname2);
	  free(fname2);

	  output_scalarfield((real *) curfield, dims, 
			     local_dims, start, file_id, "data",
			     last_dim_index, last_dim_start, last_dim_size,
			     first_dim_start, first_dim_size,
			     write_start0_special);

	  if (curfield_type == 'n') {
	       int c1, c2, inv;
	       char dataname[100];

	       for (inv = 0; inv < 2; ++inv)
		    for (c1 = 0; c1 < 3; ++c1)
			 for (c2 = c1; c2 < 3; ++c2) {
			      get_epsilon_tensor(c1,c2, 0, inv);
			      sprintf(dataname, "%s.%c%c",
				      inv ? "epsilon_inverse" : "epsilon",
				      c1 + 'x', c2 + 'x');
			      output_scalarfield((real *) curfield, dims,
						 local_dims, start,
						 file_id, dataname,
						 last_dim_index,
						 last_dim_start, last_dim_size,
						 first_dim_start,
						 first_dim_size,
						 write_start0_special);
#if defined(WITH_HERMITIAN_EPSILON)
			      if (c1 != c2) {
				   get_epsilon_tensor(c1,c2, 1, inv);
				   strcat(dataname, ".i");
#ifndef SCALAR_COMPLEX /* scalarfield_otherhalf isn't right */
				   strcat(dataname, ".screwy");
#endif
				   output_scalarfield((real *) curfield, dims,
						      local_dims, start,
						      file_id, dataname,
						      last_dim_index,
						      last_dim_start,
						      last_dim_size,
						      first_dim_start,
						      first_dim_size,
						      write_start0_special);
			      }
#endif
			 }
	  }

     }
     else
	  mpi_one_fprintf(stderr, "unknown field type!\n");

     if (file_id.id >= 0) {
	  matrixio_write_data_attr(file_id, "lattice vectors",
				   &output_R[0][0], 2, attr_dims);
	  matrixio_write_string_attr(file_id, "description", description);

	  matrixio_close(file_id);
     }

     /* We have destroyed curfield (by multiplying it by phases,
	and/or reorganizing in the case of real-amplitude fields). */
     curfield_reset();
}

/**************************************************************************/

/* For curfield an energy density, compute the fraction of the energy
   that resides inside the given list of geometric objects.   Later
   objects in the list have precedence, just like the ordinary
   geometry list. */
number compute_energy_in_object_list(geometric_object_list objects)
{
     int i, j, k, n1, n2, n3, n_other, n_last, rank, last_dim;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
     real s1, s2, s3, c1, c2, c3;
     real *energy = (real *) curfield;
     real energy_sum = 0;

     if (!curfield || !strchr("DHR", curfield_type)) {
          mpi_one_fprintf(stderr, "The D or H energy density must be loaded first.\n");
          return 0.0;
     }

     for (i = 0; i < objects.num_items; ++i)
	  geom_fix_object(objects.items[i]);

     n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;
     n_other = mdata->other_dims;
     n_last = mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     last_dim = mdata->last_dim;
     rank = (n3 == 1) ? (n2 == 1 ? 1 : 2) : 3;

     s1 = geometry_lattice.size.x / n1;
     s2 = geometry_lattice.size.y / n2;
     s3 = geometry_lattice.size.z / n3;
     c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5;
     c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
     c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

     /* Here we have different loops over the coordinates, depending
	upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and "index" describing the corresponding index in 
	the curfield array.

        This was all stolen from maxwell_eps.c...it would be better
        if we didn't have to cut and paste, sigh. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI
     
     for (i = 0; i < n1; ++i)
	  for (j = 0; j < n2; ++j)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j, k2 = k;
	  int index = ((i * n2 + j) * n3 + k);

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j + local_y_start, k2 = k;
	  int index = ((j * n1 + i) * n3 + k);

#  endif /* HAVE_MPI */

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

     for (i = 0; i < n_other; ++i)
	  for (j = 0; j < n_last; ++j)
     {
	  int index = i * n_last + j;
	  int i2, j2, k2;
	  switch (rank) {
	      case 2: i2 = i; j2 = j; k2 = 0; break;
	      case 3: i2 = i / n2; j2 = i % n2; k2 = j; break;
	      default: i2 = j; j2 = k2 = 0;  break;
	  }

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = mdata->last_dim_size / 2;
     else
	  local_n3 = 1;
     
     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < local_n3; ++k)
     {
#         define i2 i
	  int j2 = j + local_y_start;
#         define k2 k
	  int index = ((j * n1 + i) * local_n3 + k);

#  endif /* HAVE_MPI */

#endif /* not SCALAR_COMPLEX */

	  {
	       vector3 p;
	       int n;
	       p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;
	       for (n = objects.num_items - 1; n >= 0; --n)
		    if (point_in_periodic_fixed_objectp(p, objects.items[n])) {
			 if (objects.items[n].material.which_subclass
			     == MATERIAL_TYPE_SELF)
			      break; /* treat as a "nothing" object */
			 energy_sum += energy[index];
#ifndef SCALAR_COMPLEX
			 {
			      int last_index;
#  ifdef HAVE_MPI
			      if (n3 == 1)
				   last_index = j + local_y_start;
			      else
				   last_index = k;
#  else
			      last_index = j;
#  endif
			      if (last_index != 0 && 2*last_index != last_dim)
				   energy_sum += energy[index];
			 }
#endif
			 break;
		    }
	  }
     }

     mpi_allreduce_1(&energy_sum, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     energy_sum *= Vol / H.N;
     return energy_sum;
}

/**************************************************************************/

/* Compute the integral of f(energy/field, epsilon, r) over the cell. */
cnumber compute_field_integral(function f)
{
     int i, j, k, n1, n2, n3, n_other, n_last, rank, last_dim;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
     real s1, s2, s3, c1, c2, c3;
     int integrate_energy;
     real *energy = (real *) curfield;
     cnumber integral = {0,0};
     vector3 kvector = {0,0,0};

     if (!curfield || !strchr("dheDHRcv", curfield_type)) {
          mpi_one_fprintf(stderr, "The D or H energy/field must be loaded first.\n");
          return integral;
     }
     if (curfield_type != 'v')
	  kvector = cur_kvector;

     integrate_energy = strchr("DHR", curfield_type) != NULL;

     n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;
     n_other = mdata->other_dims;
     n_last = mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     last_dim = mdata->last_dim;
     rank = (n3 == 1) ? (n2 == 1 ? 1 : 2) : 3;

     s1 = geometry_lattice.size.x / n1;
     s2 = geometry_lattice.size.y / n2;
     s3 = geometry_lattice.size.z / n3;
     c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5;
     c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
     c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

     /* Here we have different loops over the coordinates, depending
	upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and "index" describing the corresponding index in 
	the curfield array.

        This was all stolen from maxwell_eps.c...it would be better
        if we didn't have to cut and paste, sigh. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI
     
     for (i = 0; i < n1; ++i)
	  for (j = 0; j < n2; ++j)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j, k2 = k;
	  int index = ((i * n2 + j) * n3 + k);

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < n3; ++k)
     {
	  int i2 = i, j2 = j + local_y_start, k2 = k;
	  int index = ((j * n1 + i) * n3 + k);

#  endif /* HAVE_MPI */

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

     for (i = 0; i < n_other; ++i)
	  for (j = 0; j < n_last; ++j)
     {
	  int index = i * n_last + j;
	  int i2, j2, k2;
	  switch (rank) {
	      case 2: i2 = i; j2 = j; k2 = 0; break;
	      case 3: i2 = i / n2; j2 = i % n2; k2 = j; break;
	      default: i2 = j; j2 = k2 = 0;  break;
	  }

#  else /* HAVE_MPI */

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = mdata->last_dim_size / 2;
     else
	  local_n3 = 1;
     
     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_n2; ++j)
          for (i = 0; i < n1; ++i)
	       for (k = 0; k < local_n3; ++k)
     {
#         define i2 i
	  int j2 = j + local_y_start;
#         define k2 k
	  int index = ((j * n1 + i) * local_n3 + k);

#  endif /* HAVE_MPI */

#endif /* not SCALAR_COMPLEX */

	  {
	       real epsilon;
	       vector3 p;

	       epsilon = mean_epsilon_from_matrix(mdata->eps_inv + index);
	       
	       p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;
	       if (integrate_energy) {
		    integral.re +=
			 ctl_convert_number_to_c(
			      gh_call3(f,
				     ctl_convert_number_to_scm(energy[index]),
				       ctl_convert_number_to_scm(epsilon),
				       ctl_convert_vector3_to_scm(p)));	  
	       }
	       else {
		    cvector3 F;
		    double phase_phi;
		    scalar_complex phase;
		    cnumber integrand;

		    phase_phi = TWOPI * 
			 (kvector.x * (p.x/geometry_lattice.size.x) +
			  kvector.y * (p.y/geometry_lattice.size.y) +
			  kvector.z * (p.z/geometry_lattice.size.z));
		    CASSIGN_SCALAR(phase, cos(phase_phi), sin(phase_phi));
		    CASSIGN_MULT_RE(F.x.re, curfield[3*index+0], phase);
		    CASSIGN_MULT_IM(F.x.im, curfield[3*index+0], phase);
		    CASSIGN_MULT_RE(F.y.re, curfield[3*index+1], phase);
		    CASSIGN_MULT_IM(F.y.im, curfield[3*index+1], phase);
		    CASSIGN_MULT_RE(F.z.re, curfield[3*index+2], phase);
		    CASSIGN_MULT_IM(F.z.im, curfield[3*index+2], phase);
		    integrand =
			 ctl_convert_cnumber_to_c(
                              gh_call3(f,
				       ctl_convert_cvector3_to_scm(F),
				       ctl_convert_number_to_scm(epsilon),
                                       ctl_convert_vector3_to_scm(p)));
		    integral.re += integrand.re;
		    integral.im += integrand.im;
	       }

#ifndef SCALAR_COMPLEX
	       {
		    int last_index;
#  ifdef HAVE_MPI
		    if (n3 == 1)
			 last_index = j + local_y_start;
		    else
			 last_index = k;
#  else
		    last_index = j;
#  endif
		    
		    if (last_index != 0 && 2*last_index != last_dim) {
			 int i2c, j2c, k2c;
			 i2c = i2 ? (n1 - i2) : 0;
			 j2c = j2 ? (n2 - j2) : 0;
			 k2c = k2 ? (n3 - k2) : 0;
			 p.x = i2c * s1 - c1; 
			 p.y = j2c * s2 - c2; 
			 p.z = k2c * s3 - c3;
			 if (integrate_energy)
			      integral.re += 
				   ctl_convert_number_to_c(
					gh_call3(f,
				      ctl_convert_number_to_scm(energy[index]),
				      ctl_convert_number_to_scm(epsilon),
				      ctl_convert_vector3_to_scm(p)));
			 else {
			      cvector3 F;
			      double phase_phi;
			      scalar_complex phase, Fx, Fy, Fz;
			      cnumber integrand;
			      
			      Fx = curfield[3*index+0];
			      Fy = curfield[3*index+1];
			      Fz = curfield[3*index+2];
			      Fx.im= -Fx.im; Fy.im= -Fy.im; Fz.im= -Fz.im;

			      phase_phi = TWOPI * 
				   (kvector.x 
				    * (p.x/geometry_lattice.size.x) +
				    kvector.y 
				    * (p.y/geometry_lattice.size.y) +
				    kvector.z 
				    * (p.z/geometry_lattice.size.z));
			      CASSIGN_SCALAR(phase, 
					     cos(phase_phi), sin(phase_phi));
			      CASSIGN_MULT_RE(F.x.re, Fx, phase);
			      CASSIGN_MULT_IM(F.x.im, Fx, phase);
			      CASSIGN_MULT_RE(F.y.re, Fy, phase);
			      CASSIGN_MULT_IM(F.y.im, Fy, phase);
			      CASSIGN_MULT_RE(F.z.re, Fz, phase);
			      CASSIGN_MULT_IM(F.z.im, Fz, phase);

			      integrand =
				   ctl_convert_cnumber_to_c(
					gh_call3(f,
					    ctl_convert_cvector3_to_scm(F),
				            ctl_convert_number_to_scm(epsilon),
                                            ctl_convert_vector3_to_scm(p)));
			      integral.re += integrand.re;
			      integral.im += integrand.im;
			 }
		    }
	       }
#endif
	  }
     }

     integral.re *= Vol / H.N;
     integral.im *= Vol / H.N;
     {
	  cnumber integral_sum;
	  mpi_allreduce(&integral, &integral_sum, 2, number, 
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  return integral_sum;
     }
}

number compute_energy_integral(function f)
{
     if (!curfield || !strchr("DHR", curfield_type)) {
          mpi_one_fprintf(stderr, "The D or H energy density must be loaded first.\n");
          return 0.0;
     }

     return cnumber_re(compute_field_integral(f));
}

/**************************************************************************/
