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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stddef.h>

#include "../src/config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <matrixio.h>

#include <ctl-io.h>
#include <ctlgeom.h>

#include "meb.h"

/**************************************************************************/

void get_ufield(integer which_band)
{
     if (!edata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-ufield!\n");
	  return;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr,
			  "solve-kpoint must be called before get-ufield!\n");
	  return;
     }
     if (which_band < 1 || which_band > V.p) {
	  mpi_one_fprintf(stderr,
			  "must have 1 <= band index <= num_bands (%d)\n",V.p);
	  return;
     }

     curfield = (scalar_complex *) edata->fft_data;
     curfield_band = which_band;
     curfield_type = 'u';
     elastic_compute_v_from_V(edata, V, curfield, which_band - 1, 1);
     elastic_compute_u_from_v(edata, curfield, 1);

     /* Divide by the cell volume so that the integral of rho*|U|^2
        is unity.  (From the eigensolver + FFT, it is
        initially normalized to sum to nx*ny*nz.) */

     {
	  int i, N;
	  double scale;
	  N = edata->fft_output_size;

	  scale = 1.0 / sqrt(Vol);
	  for (i = 0; i < 3*N; ++i) {
	       curfield[i].re *= scale;
	       curfield[i].im *= scale;
	  }
     }
}

void get_vfield(integer which_band)
{
     if (!edata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-vfield!\n");
	  return;
     }
     if (!kpoint_index) {
	  mpi_one_fprintf(stderr,
			  "solve-kpoint must be called before get-vfield!\n");
	  return;
     }
     if (which_band < 1 || which_band > V.p) {
	  mpi_one_fprintf(stderr,
			  "must have 1 <= band index <= num_bands (%d)\n",V.p);
	  return;
     }

     curfield = (scalar_complex *) edata->fft_data;
     curfield_band = which_band;
     curfield_type = 'v';
     elastic_compute_v_from_V(edata, V, curfield, which_band - 1, 1);

     /* Divide by the cell volume so that the integral of |V|^2
        is unity.  (From the eigensolver + FFT, it is
        initially normalized to sum to nx*ny*nz.) */

     {
	  int i, N;
	  double scale;
	  N = edata->fft_output_size;

	  scale = 1.0 / sqrt(Vol);
	  for (i = 0; i < 3*N; ++i) {
	       curfield[i].re *= scale;
	       curfield[i].im *= scale;
	  }
     }
}

/* get the density function, and compute some statistics */
void get_rho(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *rho;
     real rho_mean = 0, rho_inv_mean = 0, rho_high = -1e20, rho_low = 1e20;
     int fill_count = 0;

     if (!edata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-rho!\n");
	  return;
     }

     curfield = (scalar_complex *) edata->fft_data;
     rho = (real *) curfield;
     curfield_band = 0;
     curfield_type = 'r';

     N = edata->fft_output_size;
     last_dim = edata->last_dim;
     last_dim_stored =
	  edata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = edata->nx; nz = edata->nz; local_y_start = edata->local_y_start;

     for (i = 0; i < N; ++i) {
          rho[i] = edata->rho[i]; /* mean_rho(edata->rho + i); */
	  if (rho[i] < rho_low)
	       rho_low = rho[i];
	  if (rho[i] > rho_high)
	       rho_high = rho[i];
	  rho_mean += rho[i];
	  rho_inv_mean += 1/rho[i];
	  if (rho[i] > 1.0001)
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
		    rho_mean += rho[i];
		    rho_inv_mean += 1/rho[i];
		    if (rho[i] > 1.0001)
			 ++fill_count;
	       }
	  }
#endif
     }

     mpi_allreduce_1(&rho_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce_1(&rho_inv_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce_1(&rho_low, real, SCALAR_MPI_TYPE,
		     MPI_MIN, MPI_COMM_WORLD);
     mpi_allreduce_1(&rho_high, real, SCALAR_MPI_TYPE,
		     MPI_MAX, MPI_COMM_WORLD);
     mpi_allreduce_1(&fill_count, int, MPI_INT,
                   MPI_SUM, MPI_COMM_WORLD);
     N = edata->nx * edata->ny * edata->nz;
     rho_mean /= N;
     rho_inv_mean = N/rho_inv_mean;

     mpi_one_printf("rho: %g-%g, mean %g, harm. mean %g, "
		    "%g%% > 1, %g%% \"fill\"\n",
		    rho_low, rho_high, rho_mean, rho_inv_mean,
		    (100.0 * fill_count) / N, 
		    rho_high == rho_low ? 100.0 :
		    100.0 * (rho_mean-rho_low) / (rho_high-rho_low));
}

/* get the transverse velocity function ct */
void get_ct(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *ct;

     if (!edata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-ct!\n");
	  return;
     }

     curfield = (scalar_complex *) edata->fft_data;
     ct = (real *) curfield;
     curfield_band = 0;
     curfield_type = 't';

     N = edata->fft_output_size;
     last_dim = edata->last_dim;
     last_dim_stored =
	  edata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = edata->nx; nz = edata->nz; local_y_start = edata->local_y_start;

     for (i = 0; i < N; ++i) {
          ct[i] = sqrt(edata->rhoct2[i] / edata->rho[i]);
     }
}

/* get the transverse velocity funclion cl */
void get_cl(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *cl;

     if (!edata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-cl!\n");
	  return;
     }

     curfield = (scalar_complex *) edata->fft_data;
     cl = (real *) curfield;
     curfield_band = 0;
     curfield_type = 'l';

     N = edata->fft_output_size;
     last_dim = edata->last_dim;
     last_dim_stored =
	  edata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = edata->nx; nz = edata->nz; local_y_start = edata->local_y_start;

     for (i = 0; i < N; ++i) {
          cl[i] = sqrt(edata->rhocl2[i] / edata->rho[i]);
     }
}

/**************************************************************************/

/* Replace curfield (either u or v) with the scalar energy density function,
   normalized to one.  While we're at it, compute some statistics about
   the relative strength of different field components.  Also return
   the integral of the energy density, which should be unity. */
number_list compute_field_energy(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real comp_sum2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, comp_sum[6];
     real energy_sum = 0.0;
     real *energy_density = (real *) curfield;
     number_list retval = { 0, 0 };

     if (!curfield || !strchr("uv", curfield_type)) {
	  mpi_one_fprintf(stderr, "The U or V field must be loaded first.\n");
	  return retval;
     }

     N = edata->fft_output_size;
     last_dim = edata->last_dim;
     last_dim_stored =
	  edata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = edata->nx; nz = edata->nz; local_y_start = edata->local_y_start;

     for (i = 0; i < N; ++i) {
	  scalar_complex field[3];
	  real
	       comp_sqr0,comp_sqr1,comp_sqr2,comp_sqr3,comp_sqr4,comp_sqr5;

	  /* energy is either |u|^2 * rho or |v|^2. */
	  field[0] =   curfield[3*i];
	  field[1] = curfield[3*i+1];
	  field[2] = curfield[3*i+2];
	  
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
	       (comp_sqr0+comp_sqr1+comp_sqr2+comp_sqr3+comp_sqr4+comp_sqr5)
	       * (curfield_type == 'u' ? edata->rho[i] : 1.0);
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

     mpi_one_printf("%c-energy-components:, %d, %d",
	    curfield_type, kpoint_index, curfield_band);
     for (i = 0; i < 6; ++i) {
	  comp_sum[i] /= (energy_sum == 0 ? 1 : energy_sum);
	  if (i % 2 == 1)
	       mpi_one_printf(", %g", comp_sum[i] + comp_sum[i-1]);
     }
     mpi_one_printf("\n");


     /* remember that we now have energy density; denoted by capital U/V */
     curfield_type = toupper(curfield_type);

     /* The return value is a list of 7 items: the total energy,
	followed by the 6 elements of the comp_sum array (the fraction
	of the energy in the real/imag. parts of each field component). */

     retval.num_items = 7;
     CHK_MALLOC(retval.items, number, retval.num_items);

     retval.items[0] = energy_sum * Vol / V.N;

     for (i = 0; i < 6; ++i)
	  retval.items[i+1] = comp_sum[i];

     return retval;
}

/**************************************************************************/

/* Fix the phase of the current field (u/v) to a canonical value.
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

     if (!curfield || !strchr("uv", curfield_type)) {
          mpi_one_fprintf(stderr, "The U/V field must be loaded first.\n");
          return;
     }
     N = edata->fft_output_size * 3;

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
	  maxabs_index += edata->local_y_start * edata->nx * edata->nz;
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
     for (i = 0; i < V.n; ++i) {
          ASSIGN_MULT(V.data[i*V.p + curfield_band - 1], 
		      V.data[i*V.p + curfield_band - 1], phase);
     }
}

/**************************************************************************/

/* compute the fraction of the field energy that is located in the
   given range of densities: */
number compute_energy_in_rho(number rho_low, number rho_high)
{
     int N, i, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *energy = (real *) curfield;
     real rho, energy_sum = 0.0;

     if (!curfield || !strchr("UV", curfield_type)) {
          mpi_one_fprintf(stderr, "The U or V energy density must be loaded first.\n");
          return 0.0;
     }

     N = edata->fft_output_size;
     last_dim = edata->last_dim;
     last_dim_stored =
	  edata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = edata->nx; nz = edata->nz; local_y_start = edata->local_y_start;

     for (i = 0; i < N; ++i) {
	  rho = edata->rho[i];
	  if (rho >= rho_low && rho <= rho_high) {
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
     energy_sum *= Vol / V.N;
     return energy_sum;
}

/**************************************************************************/

/* Prepend the prefix to the fname, and (if parity_suffix is true)
   append a parity specifier (if any) (e.g. ".te"), returning a new
   string, which should be deallocated with free().  fname or prefix
   may be NULL, in which case they are treated as the empty string. */
static char *fix_fname(const char *fname, const char *prefix,
		       elastic_data *d, int parity_suffix)
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
     matrixio_id data_id = -1;

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

	  elastic_scalarfield_otherhalf(edata, vals);
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

     if (data_id >= 0)
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
     matrixio_id file_id = -1;
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
		  "fields, energy dens., or rho must be loaded first.\n");
	  return;
     }
     
#ifdef HAVE_MPI
     /* The first two dimensions (x and y) of the position-space fields
	are transposed when we use MPI, so we need to transpose everything. */
     dims[0] = edata->ny;
     local_dims[1] = dims[1] = edata->nx;
     local_dims[2] = dims[2] = edata->nz;
     local_dims[0] = edata->local_ny;
     start[0] = edata->local_y_start;
#  ifndef SCALAR_COMPLEX
     /* Ugh, hairy.  See also elastic_vectorfield_otherhalf. */
     if (dims[2] == 1) {
	  last_dim_index = 0;
	  first_dim_size = local_dims[0];
	  first_dim_start = dims[0] - (start[0] + local_dims[0] - 1);

	  if (start[0] == 0)
	       --first_dim_size; /* DC frequency is not in other half */
	  if (start[0] + local_dims[0] == edata->last_dim_size / 2 &&
	      dims[0] % 2 == 0) {
	       --first_dim_size; /* Nyquist frequency is not in other half */
	       ++first_dim_start;
	  }

	  last_dim_start = first_dim_start;
	  last_dim_size = first_dim_size;
     }
     else {
	  last_dim_index = 2;
	  local_dims[last_dim_index] = edata->last_dim_size / 2;
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
     output_k[0] = R[1][0]*edata->current_k[0] + R[1][1]*edata->current_k[1]
	  + R[1][2]*edata->current_k[2];
     output_k[1] = R[0][0]*edata->current_k[0] + R[0][1]*edata->current_k[1]
	  + R[0][2]*edata->current_k[2];
     output_k[2] = R[2][0]*edata->current_k[0] + R[2][1]*edata->current_k[1]
	  + R[2][2]*edata->current_k[2];
     output_R[0][0]=R[1][0]; output_R[0][1]=R[1][1]; output_R[0][2]=R[1][2];
     output_R[1][0]=R[0][0]; output_R[1][1]=R[0][1]; output_R[1][2]=R[0][2];
     output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
#else /* ! HAVE_MPI */
     dims[0] = edata->nx;
     local_dims[1] = dims[1] = edata->ny;
     local_dims[2] = dims[2] = edata->nz;
     local_dims[0] = edata->local_nx;
#  ifndef SCALAR_COMPLEX
     last_dim_index = dims[2] == 1 ? (dims[1] == 1 ? 0 : 1) : 2;
     local_dims[last_dim_index] = edata->last_dim_size / 2;
     last_dim_start = local_dims[last_dim_index];
     last_dim_size = dims[last_dim_index] - local_dims[last_dim_index];
     first_dim_start = last_dim_index ? 0 : last_dim_start;
     first_dim_size = last_dim_index ? local_dims[0] : last_dim_size;
#  endif
     start[0] = edata->local_x_start;
     output_k[0] = R[0][0]*edata->current_k[0] + R[0][1]*edata->current_k[1]
	  + R[0][2]*edata->current_k[2];
     output_k[1] = R[1][0]*edata->current_k[0] + R[1][1]*edata->current_k[1]
	  + R[1][2]*edata->current_k[2];
     output_k[2] = R[2][0]*edata->current_k[0] + R[2][1]*edata->current_k[1]
	  + R[2][2]*edata->current_k[2];
     output_R[0][0]=R[0][0]; output_R[0][1]=R[0][1]; output_R[0][2]=R[0][2];
     output_R[1][0]=R[1][0]; output_R[1][1]=R[1][1]; output_R[1][2]=R[1][2];
     output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
#endif /* ! HAVE_MPI */

     if (strchr("uv", curfield_type)) { /* outputting vector field */
	  matrixio_id data_id[6] = {-1,-1,-1,-1,-1,-1};
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
	  fname2 = fix_fname(fname, filename_prefix, edata, 1);
	  mpi_one_printf("Outputting fields to %s...\n", fname2);
	  file_id = matrixio_create(fname2);
	  free(fname2);
	  fieldio_write_complex_field(curfield, 3, dims, local_dims, start,
				      which_component, 3, output_k,
				      file_id, 0, data_id);

#ifndef SCALAR_COMPLEX
	  /* Here's where it gets hairy. */
	  elastic_vectorfield_otherhalf(edata, curfield,
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
	       if (data_id[i] >= 0)
		    matrixio_close_dataset(data_id[i]);
	  matrixio_write_data_attr(file_id, "Bloch wavevector",
				   output_k, 1, attr_dims);
     }
     else if (strchr("UVrtlR", curfield_type)) { /* scalar field */
	  if (curfield_type == 'r') {
	       sprintf(fname, "rho");
	       sprintf(description, "dielectric function, rho");
	  }
	  else if (curfield_type == 't') {
	       sprintf(fname, "ct");
	       sprintf(description, "transverse velocity, ct");
	  }
	  else if (curfield_type == 'l') {
	       sprintf(fname, "cl");
	       sprintf(description, "transverse velocity, cl");
	  }
	  else {
	       sprintf(fname, "%cpwr.k%02d.b%02d",
		       tolower(curfield_type), kpoint_index, curfield_band);
	       sprintf(description,
		       "%c field energy density, kpoint %d, band %d, freq=%g",
		       curfield_type, kpoint_index, curfield_band, 
		       freqs.items[curfield_band - 1]);
	  }
	  fname2 = fix_fname(fname, filename_prefix, edata, 
			     /* no parity suffix for material constants: */
			     !strchr("rtl", curfield_type));
	  mpi_one_printf("Outputting %s...\n", fname2);
	  file_id = matrixio_create(fname2);
	  free(fname2);

	  output_scalarfield((real *) curfield, dims, 
			     local_dims, start, file_id, "data",
			     last_dim_index, last_dim_start, last_dim_size,
			     first_dim_start, first_dim_size,
			     write_start0_special);
     }
     else
	  mpi_one_fprintf(stderr, "unknown field type!\n");

     if (file_id >= 0) {
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

     if (!curfield || !strchr("UV", curfield_type)) {
          mpi_one_fprintf(stderr, "The U or V energy density must be loaded first.\n");
          return 0.0;
     }

     for (i = 0; i < objects.num_items; ++i)
	  geom_fix_object(objects.items[i]);

     n1 = edata->nx; n2 = edata->ny; n3 = edata->nz;
     n_other = edata->other_dims;
     n_last = edata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     last_dim = edata->last_dim;
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

        This was all stolen from elastic_eps.c...it would be better
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

     local_n2 = edata->local_ny;
     local_y_start = edata->local_y_start;

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

     local_n2 = edata->local_ny;
     local_y_start = edata->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = edata->last_dim_size / 2;
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
     energy_sum *= Vol / V.N;
     return energy_sum;
}

/**************************************************************************/

/* Compute the integral of f(energy/field, rho, r) over the cell. */
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

     if (!curfield || !strchr("uvUV", curfield_type)) {
          mpi_one_fprintf(stderr, "The U or V energy/field must be loaded first.\n");
          return integral;
     }

     integrate_energy = strchr("UV", curfield_type) != NULL;

     n1 = edata->nx; n2 = edata->ny; n3 = edata->nz;
     n_other = edata->other_dims;
     n_last = edata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     last_dim = edata->last_dim;
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

        This was all stolen from elastic_eps.c...it would be better
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

     local_n2 = edata->local_ny;
     local_y_start = edata->local_y_start;

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

     local_n2 = edata->local_ny;
     local_y_start = edata->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = edata->last_dim_size / 2;
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
	       real rho;
	       vector3 p;

	       rho = edata->rho[index]; /* mean_rho(edata->rho_inv + index); */
	       
	       p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;
	       if (integrate_energy) {
		    integral.re +=
			 ctl_convert_number_to_c(
			      gh_call3(f,
				     ctl_convert_number_to_scm(energy[index]),
				       ctl_convert_number_to_scm(rho),
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
				       ctl_convert_number_to_scm(rho),
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
				      ctl_convert_number_to_scm(rho),
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
				            ctl_convert_number_to_scm(rho),
                                            ctl_convert_vector3_to_scm(p)));
			      integral.re += integrand.re;
			      integral.im += integrand.im;
			 }
		    }
	       }
#endif
	  }
     }

     integral.re *= Vol / V.N;
     integral.im *= Vol / V.N;
     {
	  cnumber integral_sum;
	  mpi_allreduce(&integral, &integral_sum, 2, number, 
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  return integral_sum;
     }
}

number compute_energy_integral(function f)
{
     if (!curfield || !strchr("UV", curfield_type)) {
          mpi_one_fprintf(stderr, "The U or V energy density must be loaded first.\n");
          return 0.0;
     }

     return cnumber_re(compute_field_integral(f));
}

/**************************************************************************/
