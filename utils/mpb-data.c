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
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <ctlgeom.h>

#include "config.h"

#include "check.h"
#include "matrixio.h"

int verbose = 0;

/* A macro to set x = fractional part of x input, xi = integer part,
   with 0 <= x < 1.0.   Note that we need the second test (if x >= 1.0)
   below, because x may start out as -0 or -1e-23 or something so that
   it is < 0 but x + 1.0 == 1.0, thanks to the wonders of floating point.
   (This has actually happened, on an Alpha.) */
#define MODF_POSITIVE(x, xi)                                                                       \
  {                                                                                                \
    x = modf(x, &xi);                                                                              \
    if (x < 0) {                                                                                   \
      x += 1.0;                                                                                    \
      if (x >= 1.0)                                                                                \
        x = 0;                                                                                     \
      else                                                                                         \
        xi -= 1.0;                                                                                 \
    }                                                                                              \
  }

#define ADJ_POINT(i1, i2, nx, dx, xi, xi2)                                                         \
  {                                                                                                \
    if (dx >= 0.0) {                                                                               \
      i2 = i1 + 1;                                                                                 \
      if (i2 >= nx) {                                                                              \
        i2 -= nx;                                                                                  \
        xi2 = xi + 1.0;                                                                            \
      }                                                                                            \
      else                                                                                         \
        xi2 = xi;                                                                                  \
    }                                                                                              \
    else {                                                                                         \
      i2 = i1 - 1;                                                                                 \
      if (i2 < 0) {                                                                                \
        i2 += nx;                                                                                  \
        xi2 = xi - 1.0;                                                                            \
      }                                                                                            \
      else                                                                                         \
        xi2 = xi;                                                                                  \
      dx = -dx;                                                                                    \
    }                                                                                              \
  }

void add_cmplx_times_phase(real *sum_re, real *sum_im, real d_re, real d_im, double ix, double iy,
                           double iz, real *s, real scale_by) {
  static real phase = 0.0, p_re = 1.0, p_im = 0.0;
  real new_phase;

  new_phase = ix * s[0] + iy * s[1] + iz * s[2];
  if (new_phase != phase) {
    phase = new_phase;
    p_re = cos(phase);
    p_im = sin(phase);
  }
  *sum_re += (d_re * p_re - d_im * p_im) * scale_by;
  *sum_im += (d_re * p_im + d_im * p_re) * scale_by;
}

#define TWOPI 6.2831853071795864769252867665590057683943388

#define MAX2(a, b) ((a) >= (b) ? (a) : (b))
#define MIN2(a, b) ((a) < (b) ? (a) : (b))

void map_data(real *d_in_re, real *d_in_im, int n_in[3], real *d_out_re, real *d_out_im,
              int n_out[3], matrix3x3 coord_map, real *kvector, short pick_nearest,
              short transpose) {
  int i, j, k;
  real s[3]; /* phase difference per cell in each lattice direction */
  real min_out_re = 1e20, max_out_re = -1e20, min_out_im = 1e20, max_out_im = -1e20;
  real shiftx, shifty, shiftz;

  CHECK(d_in_re && d_out_re, "invalid arguments");
  CHECK((d_out_im && d_in_im) || (!d_out_im && !d_in_im),
        "both input and output must be real or complex");

  coord_map.c0 = vector3_scale(1.0 / n_out[0], coord_map.c0);
  coord_map.c1 = vector3_scale(1.0 / n_out[1], coord_map.c1);
  coord_map.c2 = vector3_scale(1.0 / n_out[2], coord_map.c2);

  for (i = 0; i < 3; ++i) {
    if (kvector)
      s[i] = kvector[i] * TWOPI;
    else
      s[i] = 0;
  }

  /* Compute shift so that the origin of the output cell
     is mapped to the origin of the original primitive cell: */
  shiftx = 0.5 - (coord_map.c0.x * 0.5 * n_out[0] + coord_map.c1.x * 0.5 * n_out[1] +
                  coord_map.c2.x * 0.5 * n_out[2]);
  shifty = 0.5 - (coord_map.c0.y * 0.5 * n_out[0] + coord_map.c1.y * 0.5 * n_out[1] +
                  coord_map.c2.y * 0.5 * n_out[2]);
  shiftz = 0.5 - (coord_map.c0.z * 0.5 * n_out[0] + coord_map.c1.z * 0.5 * n_out[1] +
                  coord_map.c2.z * 0.5 * n_out[2]);

  for (i = 0; i < n_out[0]; ++i)
    for (j = 0; j < n_out[1]; ++j)
      for (k = 0; k < n_out[2]; ++k) {
        real x, y, z;
        double xi, yi, zi, xi2, yi2, zi2;
        double dx, dy, dz, mdx, mdy, mdz;
        int i1, j1, k1, i2, j2, k2;
        int ijk;

        if (transpose)
          ijk = (j * n_out[0] + i) * n_out[2] + k;
        else
          ijk = (i * n_out[1] + j) * n_out[2] + k;

        /* find the point corresponding to d_out[i,j,k] in
           the input array, and also find the next-nearest
           points. */
        x = coord_map.c0.x * i + coord_map.c1.x * j + coord_map.c2.x * k + shiftx;
        y = coord_map.c0.y * i + coord_map.c1.y * j + coord_map.c2.y * k + shifty;
        z = coord_map.c0.z * i + coord_map.c1.z * j + coord_map.c2.z * k + shiftz;
        MODF_POSITIVE(x, xi);
        MODF_POSITIVE(y, yi);
        MODF_POSITIVE(z, zi);
        i1 = x * n_in[0];
        j1 = y * n_in[1];
        k1 = z * n_in[2];
        dx = x * n_in[0] - i1;
        dy = y * n_in[1] - j1;
        dz = z * n_in[2] - k1;
        ADJ_POINT(i1, i2, n_in[0], dx, xi, xi2);
        ADJ_POINT(j1, j2, n_in[1], dy, yi, yi2);
        ADJ_POINT(k1, k2, n_in[2], dz, zi, zi2);

        /* dx, mdx, etcetera, are the weights for the various
           points in the input data, which we use for linearly
           interpolating to get the output point. */
        if (pick_nearest) {
          /* don't interpolate */
          dx = dx <= 0.5 ? 0.0 : 1.0;
          dy = dy <= 0.5 ? 0.0 : 1.0;
          dz = dz <= 0.5 ? 0.0 : 1.0;
        }
        mdx = 1.0 - dx;
        mdy = 1.0 - dy;
        mdz = 1.0 - dz;

        /* Now, linearly interpolate the input to get the
           output.  If the input/output are complex, we
           also need to multiply by the appropriate phase
           factor, depending upon which unit cell we are in. */

#define IN_INDEX(i, j, k) ((i * n_in[1] + j) * n_in[2] + k)
        if (d_out_im) {
          d_out_re[ijk] = 0.0;
          d_out_im[ijk] = 0.0;
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i1, j1, k1)],
                                d_in_im[IN_INDEX(i1, j1, k1)], xi, yi, zi, s, mdx * mdy * mdz);
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i1, j1, k2)],
                                d_in_im[IN_INDEX(i1, j1, k2)], xi, yi, zi2, s, mdx * mdy * dz);
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i1, j2, k1)],
                                d_in_im[IN_INDEX(i1, j2, k1)], xi, yi2, zi, s, mdx * dy * mdz);
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i1, j2, k2)],
                                d_in_im[IN_INDEX(i1, j2, k2)], xi, yi2, zi2, s, mdx * dy * dz);
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i2, j1, k1)],
                                d_in_im[IN_INDEX(i2, j1, k1)], xi2, yi, zi, s, dx * mdy * mdz);
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i2, j1, k2)],
                                d_in_im[IN_INDEX(i2, j1, k2)], xi2, yi, zi2, s, dx * mdy * dz);
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i2, j2, k1)],
                                d_in_im[IN_INDEX(i2, j2, k1)], xi2, yi2, zi, s, dx * dy * mdz);
          add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk, d_in_re[IN_INDEX(i2, j2, k2)],
                                d_in_im[IN_INDEX(i2, j2, k2)], xi2, yi2, zi2, s, dx * dy * dz);
          min_out_im = MIN2(min_out_im, d_out_im[ijk]);
          max_out_im = MAX2(max_out_im, d_out_im[ijk]);
        }
        else {
          d_out_re[ijk] = d_in_re[IN_INDEX(i1, j1, k1)] * mdx * mdy * mdz +
                          d_in_re[IN_INDEX(i1, j1, k2)] * mdx * mdy * dz +
                          d_in_re[IN_INDEX(i1, j2, k1)] * mdx * dy * mdz +
                          d_in_re[IN_INDEX(i1, j2, k2)] * mdx * dy * dz +
                          d_in_re[IN_INDEX(i2, j1, k1)] * dx * mdy * mdz +
                          d_in_re[IN_INDEX(i2, j1, k2)] * dx * mdy * dz +
                          d_in_re[IN_INDEX(i2, j2, k1)] * dx * dy * mdz +
                          d_in_re[IN_INDEX(i2, j2, k2)] * dx * dy * dz;
        }
        min_out_re = MIN2(min_out_re, d_out_re[ijk]);
        max_out_re = MAX2(max_out_re, d_out_re[ijk]);
#undef IN_INDEX
      }

  if (verbose) {
    printf("real part range: %g .. %g\n", min_out_re, max_out_re);
    if (d_out_im) printf("imag part range: %g .. %g\n", min_out_im, max_out_im);
  }
}

void handle_dataset(matrixio_id in_file, matrixio_id out_file, const char *name_re,
                    const char *name_im, matrix3x3 Rout, matrix3x3 coord_map, real *kvector,
                    double resolution, scalar_complex scaleby, real multiply_size[3],
                    int pick_nearest, int transpose) {
  real *d_in_re = NULL, *d_in_im = NULL, *d_out_re = NULL, *d_out_im = NULL;
  int in_dims[3] = {1, 1, 1}, out_dims[3] = {1, 1, 1}, out_dims2[3], rank = 3;
  int i, N;
  int start[3] = {0, 0, 0};
  matrixio_id data_id;
  char out_name[1000];

  d_in_re = matrixio_read_real_data(in_file, name_re, &rank, in_dims, 0, 0, 0, NULL);
  if (!d_in_re) goto done;

  if (verbose) printf("Found dataset %s...\n", name_re);

  if (name_im) {
    d_in_im = matrixio_read_real_data(in_file, name_im, &rank, out_dims, 0, 0, 0, NULL);
    if (!d_in_im) {
      fprintf(stderr, "mpb-data: found %s dataset but not %s\n", name_re, name_im);
      goto done;
    }

    for (i = 0; i < 3; ++i) {
      CHECK(out_dims[i] == in_dims[i], "re/im datasets must have same size!");
    }

    if (verbose) printf("   and imaginary part dataset %s...\n", name_im);
  }

  if (verbose)
    printf("Input data is rank %d, size %dx%dx%d.\n", rank, in_dims[0], in_dims[1], in_dims[2]);

  if (resolution > 0) {
    out_dims[0] = vector3_norm(Rout.c0) * resolution + 0.5;
    out_dims[1] = vector3_norm(Rout.c1) * resolution + 0.5;
    out_dims[2] = vector3_norm(Rout.c2) * resolution + 0.5;
  }
  else {
    for (i = 0; i < 3; ++i)
      out_dims[i] = in_dims[i] * multiply_size[i];
  }
  for (i = rank; i < 3; ++i)
    out_dims[i] = 1;
  for (N = 1, i = 0; i < 3; ++i)
    N *= (out_dims[i] = MAX2(out_dims[i], 1));

  if (transpose) {
    out_dims2[0] = out_dims[1];
    out_dims2[1] = out_dims[0];
    out_dims2[2] = out_dims[2];
  }
  else {
    out_dims2[0] = out_dims[0];
    out_dims2[1] = out_dims[1];
    out_dims2[2] = out_dims[2];
  }

  if (verbose) printf("Output data %dx%dx%d.\n", out_dims2[0], out_dims2[1], out_dims2[2]);

  CHK_MALLOC(d_out_re, real, N);
  if (d_in_im) { CHK_MALLOC(d_out_im, real, N); }

  map_data(d_in_re, d_in_im, in_dims, d_out_re, d_out_im, out_dims, coord_map, kvector,
           pick_nearest, transpose);

  if (d_out_im) { /* multiply * scaleby for complex data */
    for (i = 0; i < N; ++i) {
      scalar_complex d;
      CASSIGN_SCALAR(d, d_out_re[i], d_out_im[i]);
      CASSIGN_MULT(d, scaleby, d);
      d_out_re[i] = CSCALAR_RE(d);
      d_out_im[i] = CSCALAR_IM(d);
    }
  }

  strcpy(out_name, name_re);
  if (out_file.id == in_file.id) strcat(out_name, "-new");
  if (verbose) printf("Writing dataset to %s...\n", out_name);
  data_id = matrixio_create_dataset(out_file, out_name, "", rank, out_dims2);
  matrixio_write_real_data(data_id, out_dims2, start, 1, d_out_re);
  matrixio_close_dataset(data_id);

  if (d_out_im) {
    strcpy(out_name, name_im);
    if (out_file.id == in_file.id) strcat(out_name, "-new");
    if (verbose) printf("Writing dataset to %s...\n", out_name);
    data_id = matrixio_create_dataset(out_file, out_name, "", rank, out_dims2);
    matrixio_write_real_data(data_id, out_dims2, start, 1, d_out_im);
    matrixio_close_dataset(data_id);
  }

  if (verbose) printf("Successfully wrote out data.\n");

done:
  free(d_in_re);
  free(d_in_im);
  free(d_out_re);
  free(d_out_im);
}

void handle_cvector_dataset(matrixio_id in_file, matrixio_id out_file, matrix3x3 Rout,
                            matrix3x3 coord_map, matrix3x3 cart_map, real *kvector,
                            double resolution, scalar_complex scaleby, real multiply_size[3],
                            int pick_nearest, int transpose) {
  real *d_in[3][2] = {{0, 0}, {0, 0}, {0, 0}};
  int in_dims[3] = {1, 1, 1}, out_dims[3] = {1, 1, 1}, out_dims2[3], rank = 3;
  int i, N, dim, ri;
  int start[3] = {0, 0, 0};
  matrixio_id data_id;

  for (dim = 0; dim < 3; ++dim)
    for (ri = 0; ri < 2; ++ri) {
      char nam[] = "x.r";
      int dims[3] = {1, 1, 1}, rnk = 3;

      nam[0] = 'x' + dim;
      nam[2] = ri ? 'i' : 'r';
      d_in[dim][ri] = matrixio_read_real_data(in_file, nam, &rnk, dims, 0, 0, 0, NULL);
      if (!d_in[dim][ri]) goto bad;
      if (!dim && !ri) {
        rank = rnk;
        for (i = 0; i < 3; ++i)
          in_dims[i] = dims[i];
      }
      else {
        if (rank != rnk || in_dims[0] != dims[0] || in_dims[1] != dims[1] || in_dims[2] != dims[2])
          goto bad;
      }
    }

  if (verbose) printf("Found complex vector dataset...\n");

  if (verbose)
    printf("Input data is rank %d, size %dx%dx%d.\n", rank, in_dims[0], in_dims[1], in_dims[2]);

  /* rotate vector field according to cart_map */
  if (verbose)
    printf("Rotating vectors by matrix [ %10f%10f%10f\n"
           "                             %10f%10f%10f\n"
           "                             %10f%10f%10f ]\n",
           cart_map.c0.x, cart_map.c1.x, cart_map.c2.x, cart_map.c0.y, cart_map.c1.y, cart_map.c2.y,
           cart_map.c0.z, cart_map.c1.z, cart_map.c2.z);
  N = in_dims[0] * in_dims[1] * in_dims[2];
  for (ri = 0; ri < 2; ++ri)
    for (i = 0; i < N; ++i) {
      vector3 v;
      v.x = d_in[0][ri][i];
      v.y = d_in[1][ri][i];
      v.z = d_in[2][ri][i];
      v = matrix3x3_vector3_mult(cart_map, v);
      d_in[0][ri][i] = v.x;
      d_in[1][ri][i] = v.y;
      d_in[2][ri][i] = v.z;
    }

  if (resolution > 0) {
    out_dims[0] = vector3_norm(Rout.c0) * resolution + 0.5;
    out_dims[1] = vector3_norm(Rout.c1) * resolution + 0.5;
    out_dims[2] = vector3_norm(Rout.c2) * resolution + 0.5;
  }
  else {
    for (i = 0; i < 3; ++i)
      out_dims[i] = in_dims[i] * multiply_size[i];
  }
  for (i = rank; i < 3; ++i)
    out_dims[i] = 1;
  for (N = 1, i = 0; i < 3; ++i)
    N *= (out_dims[i] = MAX2(out_dims[i], 1));

  if (transpose) {
    out_dims2[0] = out_dims[1];
    out_dims2[1] = out_dims[0];
    out_dims2[2] = out_dims[2];
  }
  else {
    out_dims2[0] = out_dims[0];
    out_dims2[1] = out_dims[1];
    out_dims2[2] = out_dims[2];
  }

  if (verbose) printf("Output data %dx%dx%d.\n", out_dims2[0], out_dims2[1], out_dims2[2]);

  for (dim = 0; dim < 3; ++dim) {
    real *d_out_re, *d_out_im;
    char nam[] = "x.r-new";

    CHK_MALLOC(d_out_re, real, N);
    CHK_MALLOC(d_out_im, real, N);

    map_data(d_in[dim][0], d_in[dim][1], in_dims, d_out_re, d_out_im, out_dims, coord_map, kvector,
             pick_nearest, transpose);

    for (i = 0; i < N; ++i) { /* multiply * scaleby */
      scalar_complex d;
      CASSIGN_SCALAR(d, d_out_re[i], d_out_im[i]);
      CASSIGN_MULT(d, scaleby, d);
      d_out_re[i] = CSCALAR_RE(d);
      d_out_im[i] = CSCALAR_IM(d);
    }

    nam[0] = 'x' + dim;
    if (out_file.id != in_file.id) nam[3] = 0;
    if (verbose) printf("Writing dataset to %s...\n", nam);
    data_id = matrixio_create_dataset(out_file, nam, "", rank, out_dims2);
    matrixio_write_real_data(data_id, out_dims2, start, 1, d_out_re);
    matrixio_close_dataset(data_id);

    nam[2] = 'i';
    if (verbose) printf("Writing dataset to %s...\n", nam);
    data_id = matrixio_create_dataset(out_file, nam, "", rank, out_dims2);
    matrixio_write_real_data(data_id, out_dims2, start, 1, d_out_im);
    matrixio_close_dataset(data_id);

    if (verbose) printf("Successfully wrote out data.\n");
  }

  for (dim = 0; dim < 3; ++dim)
    for (ri = 0; ri < 2; ++ri)
      free(d_in[dim][ri]);
  return;

bad:
  for (dim = 0; dim < 3; ++dim)
    for (ri = 0; ri < 2; ++ri)
      free(d_in[dim][ri]);
  /* try individual datasets */
  for (dim = 0; dim < 3; ++dim) {
    char namr[] = "x.r";
    char nami[] = "x.i";

    namr[0] = 'x' + dim;
    nami[0] = 'x' + dim;
    handle_dataset(in_file, out_file, namr, nami, Rout, coord_map, kvector, resolution, scaleby,
                   multiply_size, pick_nearest, transpose);

    namr[1] = 0;
    handle_dataset(in_file, out_file, namr, NULL, Rout, coord_map, kvector, resolution, scaleby,
                   multiply_size, pick_nearest, transpose);
  }
}

void handle_file(const char *fname, const char *out_fname, const char *data_name, int rectify,
                 int have_ve, vector3 ve, double resolution, scalar_complex scaleby,
                 real multiply_size[3], int pick_nearest, int transpose) {
  matrixio_id in_file, out_file;
  real *R, *kvector, *copies;
  int dims[2], rank;
  matrix3x3 Rin = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, Rout, coord_map;
  matrix3x3 cart_map = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
#define NUM_DATANAMES 13
  char datanames[NUM_DATANAMES][30] = {"data",
                                       "epsilon.xx",
                                       "epsilon.xy",
                                       "epsilon.xz",
                                       "epsilon.yy",
                                       "epsilon.yz",
                                       "epsilon.zz",
                                       "epsilon_inverse.xx",
                                       "epsilon_inverse.xy",
                                       "epsilon_inverse.xz",
                                       "epsilon_inverse.yy",
                                       "epsilon_inverse.yz",
                                       "epsilon_inverse.zz"};
  int i;

  if (verbose) printf("Reading file %s...\n", fname);
  in_file = matrixio_open(fname, out_fname != NULL);

  if (data_name && !data_name[0]) data_name = NULL;

  R = matrixio_read_data_attr(in_file, "lattice vectors", &rank, 2, dims);
  if (R && rank == 2 && dims[0] == 3 && dims[1] == 3) {
    Rin.c0.x = R[0 * 3 + 0];
    Rin.c0.y = R[0 * 3 + 1];
    Rin.c0.z = R[0 * 3 + 2];
    Rin.c1.x = R[1 * 3 + 0];
    Rin.c1.y = R[1 * 3 + 1];
    Rin.c1.z = R[1 * 3 + 2];
    Rin.c2.x = R[2 * 3 + 0];
    Rin.c2.y = R[2 * 3 + 1];
    Rin.c2.z = R[2 * 3 + 2];
    if (verbose) printf("Read lattice vectors.\n");
  }
  free(R);

  kvector = matrixio_read_data_attr(in_file, "Bloch wavevector", &rank, 1, dims);
  if (rank != 1 || dims[0] != 3) {
    free(kvector);
    kvector = NULL;
  }
  else if (verbose)
    printf("Read Bloch wavevector (%g, %g, %g)\n", kvector[0], kvector[1], kvector[2]);

  copies = matrixio_read_data_attr(in_file, "lattice copies", &rank, 1, dims);
  if (copies && rank == 1 && dims[0] == 3) {
    Rin.c0 = vector3_scale(copies[0], Rin.c0);
    Rin.c1 = vector3_scale(copies[1], Rin.c1);
    Rin.c2 = vector3_scale(copies[2], Rin.c2);
    if (kvector) {
      kvector[0] *= copies[0];
      kvector[1] *= copies[1];
      kvector[2] *= copies[2];
    }
    if (verbose) printf("Read lattice copies (%g, %g, %g)\n", copies[0], copies[1], copies[2]);
  }
  free(copies);

  if (verbose)
    printf("Input lattice = (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n", Rin.c0.x, Rin.c0.y, Rin.c0.z,
           Rin.c1.x, Rin.c1.y, Rin.c1.z, Rin.c2.x, Rin.c2.y, Rin.c2.z);

  Rout = Rin;

  if (rectify) {
    double V;

    /* Orthogonalize the output lattice vectors.  If have_ve
       is true, then the first new lattice vector should be in
       the direction of the ve unit vector; otherwise, the first
       new lattice vector is the first original lattice vector.
       Note that we do this in such a way as to preserve the
       volume of the unit cell, and so that our first vector
       (in the direction of ve) smoothly interpolates between
       the original lattice vectors. */

    if (have_ve)
      ve = unit_vector3(ve);
    else
      ve = unit_vector3(Rout.c0);

    /* First, compute c0 in the direction of ve by smoothly
       interpolating the old c0/c1/c2 (formula is slightly tricky): */
    V = vector3_dot(vector3_cross(Rout.c0, Rout.c1), Rout.c2);
    Rout.c1 = vector3_minus(Rout.c1, Rout.c0);
    Rout.c2 = vector3_minus(Rout.c2, Rout.c0);
    Rout.c0 = vector3_scale(V / vector3_dot(vector3_cross(Rout.c1, Rout.c2), ve), ve);

    /* Now, orthogonalize c1 and c2: */
    Rout.c1 = vector3_minus(Rout.c1, vector3_scale(vector3_dot(ve, Rout.c1), ve));
    Rout.c2 = vector3_minus(Rout.c2, vector3_scale(vector3_dot(ve, Rout.c2), ve));
    Rout.c2 = vector3_minus(
        Rout.c2,
        vector3_scale(vector3_dot(Rout.c1, Rout.c2) / vector3_dot(Rout.c1, Rout.c1), Rout.c1));

    cart_map.c0 = unit_vector3(Rout.c0);
    cart_map.c1 = unit_vector3(Rout.c1);
    cart_map.c2 = unit_vector3(Rout.c2);
    cart_map = matrix3x3_inverse(cart_map);
  }

  if (transpose) { /* swap first two rows of cart_map */
    vector3 v;
    cart_map = matrix3x3_transpose(cart_map);
    v = cart_map.c0;
    cart_map.c0 = cart_map.c1;
    cart_map.c1 = v;
    cart_map = matrix3x3_transpose(cart_map);
  }

  Rout.c0 = vector3_scale(multiply_size[0], Rout.c0);
  Rout.c1 = vector3_scale(multiply_size[1], Rout.c1);
  Rout.c2 = vector3_scale(multiply_size[2], Rout.c2);

  if (verbose)
    printf("Output lattice = (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n", Rout.c0.x, Rout.c0.y, Rout.c0.z,
           Rout.c1.x, Rout.c1.y, Rout.c1.z, Rout.c2.x, Rout.c2.y, Rout.c2.z);

  coord_map = matrix3x3_mult(matrix3x3_inverse(Rin), Rout);

  if (out_fname) {
    if (verbose) printf("Creating output file %s...\n", out_fname);
    out_file = matrixio_create(out_fname);
  }
  else {
    if (verbose) printf("Writing output datasets to input file %s...\n", fname);
    out_file = in_file;
  }

  for (i = 0; i < NUM_DATANAMES; ++i) {
    const char *dname = datanames[i];
    char name_re[300], name_im[300];

    if (data_name) dname = data_name;

    strcpy(name_re, dname);
    handle_dataset(in_file, out_file, name_re, NULL, Rout, coord_map, kvector, resolution, scaleby,
                   multiply_size, pick_nearest, transpose);

    sprintf(name_re, "%s.r", dname);
    sprintf(name_im, "%s.i", dname);
    handle_dataset(in_file, out_file, name_re, name_im, Rout, coord_map, kvector, resolution,
                   scaleby, multiply_size, pick_nearest, transpose);

    if (data_name) break;
  }

  /* handle complex vector fields x.{ri}, y.{ri}, z.{ri} */
  handle_cvector_dataset(in_file, out_file, Rout, coord_map, cart_map, kvector, resolution, scaleby,
                         multiply_size, pick_nearest, transpose);

  free(kvector);

  matrixio_close(in_file);
  if (out_file.id != in_file.id) matrixio_close(out_file);
}

void usage(FILE *f) {
  fprintf(f, "Usage: mpb-data [options] [<filenames>]\n"
             "Options:\n"
             "         -h : this help message\n"
             "         -V : print version number and copyright\n"
             "         -v : verbose output\n"
             "  -o <file> : output to <file> (first input file only)\n"
             "         -r : output rectangular cell\n"
             " -e <x,y,z> : as -r, but first axis of cell is along <x,y,z>\n"
             "     -n <n> : output resolution of n grid points per a\n"
             "    -x <mx>\n"
             "    -y <my>\n"
             "    -z <mx> : output mx/my/mz periods in the x/y/z directions\n"
             " -P <angle> : multiply phase shift of <angle> degrees\n"
             "     -m <s> : same as -x <s> -y <s> -z <s>\n"
             "         -T : transpose first two dimensions (x & y) of data\n"
             "         -p : pixellized output (no grid interpolation)\n"
             "  -d <name> : use dataset <name> in the input files (default: all mpb datasets)\n"
             "           -- you can also specify a dataset via <filename>:<name>\n");
}

/* given an fname of the form <filename>:<data_name>, return a pointer
   to a newly-allocated string containing <filename>, and point data_name
   to the position of <data_name> in fname.  The user must free() the
<filename> string. */
static char *split_fname(char *fname, char **data_name) {
  int fname_len;
  char *colon, *filename;

  fname_len = strlen(fname);
  colon = strchr(fname, ':');
  if (colon) {
    int colon_len = strlen(colon);
    filename = (char *)malloc(sizeof(char) * (fname_len - colon_len + 1));
    CHECK(filename, "out of memory");
    strncpy(filename, fname, fname_len - colon_len + 1);
    filename[fname_len - colon_len] = 0;
    *data_name = colon + 1;
  }
  else { /* treat as if ":" were at the end of fname */
    filename = (char *)malloc(sizeof(char) * (fname_len + 1));
    CHECK(filename, "out of memory");
    strcpy(filename, fname);
    *data_name = fname + fname_len;
  }
  return filename;
}

int main(int argc, char **argv) {
  char *out_fname = NULL, *data_name = NULL;
  int rectify = 0, have_ve = 0;
  double phaseangle = 0;
  double resolution = 0;
  vector3 ve = {1, 0, 0};
  real multiply_size[3] = {1, 1, 1};
  int pick_nearest = 0, transpose = 0;
  int ifile, c;
  extern char *optarg;
  extern int optind;
  scalar_complex scaleby = {1, 0}, phase;

  while ((c = getopt(argc, argv, "hVvo:x:y:z:m:d:n:prTe:P:")) != -1)
    switch (c) {
      case 'h': usage(stdout); return EXIT_SUCCESS;
      case 'V':
        printf("mpb-data " PACKAGE_VERSION " by Steven G. Johnson.\n"
               "Copyright (C) 1999-2014 Massachusetts Institute of Technology.\n"
               "This is free software, and you are welcome to redistribute it under the\n"
               "terms of the GNU General Public License (GPL).  mpb-data comes with\n"
               "ABSOLUTELY NO WARRANTY; see the GPL for more details.\n");
        return EXIT_SUCCESS;
      case 'v': verbose = 1; break;
      case 'o':
        free(out_fname);
        out_fname = (char *)malloc(sizeof(char) * (strlen(optarg) + 1));
        CHECK(out_fname, "out of memory");
        strcpy(out_fname, optarg);
        break;
      case 'd':
        free(data_name);
        data_name = (char *)malloc(sizeof(char) * (strlen(optarg) + 1));
        CHECK(data_name, "out of memory");
        strcpy(data_name, optarg);
        break;
      case 'x': multiply_size[0] = atof(optarg); break;
      case 'y': multiply_size[1] = atof(optarg); break;
      case 'z': multiply_size[2] = atof(optarg); break;
      case 'm':
        multiply_size[0] = atof(optarg);
        multiply_size[1] = atof(optarg);
        multiply_size[2] = atof(optarg);
        break;
      case 'n':
        resolution = atof(optarg);
        CHECK(resolution > 0, "invalid resolution for -n (must be positive)");
        break;
      case 'P': phaseangle = atof(optarg); break;
      case 'p': pick_nearest = 1; break;
      case 'T': transpose = 1; break;
      case 'e':
        have_ve = 1;
        if (3 != sscanf(optarg, "%lf,%lf,%lf", &ve.x, &ve.y, &ve.z)) {
          fprintf(stderr, "Invalid -e argument \"%s\"\n", optarg);
          usage(stderr);
          return EXIT_FAILURE;
        }
        rectify = 1;
        break;
      case 'r': rectify = 1; break;
      default:
        fprintf(stderr, "Invalid argument -%c\n", c);
        usage(stderr);
        return EXIT_FAILURE;
    }
  if (optind == argc) { /* no parameters left */
    usage(stderr);
    return EXIT_FAILURE;
  }

  CASSIGN_SCALAR(phase, cos(TWOPI * phaseangle / 360.0), sin(TWOPI * phaseangle / 360.0));
  CASSIGN_MULT(scaleby, scaleby, phase);

  for (ifile = optind; ifile < argc; ++ifile) {
    char *dname, *h5_fname;
    h5_fname = split_fname(argv[ifile], &dname);
    if (!dname[0]) dname = data_name;

    handle_file(h5_fname, out_fname, dname, rectify, have_ve, ve, resolution, scaleby,
                multiply_size, pick_nearest, transpose);

    if (out_fname) free(out_fname);
    out_fname = NULL;
    free(h5_fname);
  }
  free(data_name);

  return EXIT_SUCCESS;
}
