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
#include <sphere-quad.h>
#include <mpiglue.h>
#include <mpi_utils.h>

#include "maxwell.h"
#include "xyz_loop.h"

/**************************************************************************/

/* Lapack eigenvalue functions */
#ifdef SCALAR_SINGLE_PREC
#define HEEV F77_FUNC(cheev, CHEEV)
#define SYEV F77_FUNC(ssyev, SSYEV)
#else
#define HEEV F77_FUNC(zheev, ZHEEV)
#define SYEV F77_FUNC(dsyev, DSYEV)
#endif
extern void HEEV(char *, char *, int *, scalar_complex *, int *, real *, scalar_complex *, int *,
                 real *, int *);
extern void SYEV(char *, char *, int *, real *, int *, real *, real *, int *, int *);

/* compute the 3 real eigenvalues of the matrix V */
void maxwell_sym_matrix_eigs(real eigs[3], const symmetric_matrix *V) {
  int n = 3, nw = 9, info;
#if defined(WITH_HERMITIAN_EPSILON)
  scalar_complex Vm[3][3], W[9];
  real W2[9];
  CASSIGN_SCALAR(Vm[0][0], V->m00, 0);
  CASSIGN_SCALAR(Vm[1][1], V->m11, 0);
  CASSIGN_SCALAR(Vm[2][2], V->m22, 0);
  Vm[0][1] = V->m01;
  CASSIGN_CONJ(Vm[1][0], V->m01);
  Vm[0][2] = V->m02;
  CASSIGN_CONJ(Vm[2][0], V->m02);
  Vm[1][2] = V->m12;
  CASSIGN_CONJ(Vm[2][1], V->m12);
  HEEV("V", "U", &n, &Vm[0][0], &n, eigs, W, &nw, W2, &info);
#else
  real Vm[3][3], W[9];
  Vm[0][0] = V->m00;
  Vm[1][1] = V->m11;
  Vm[2][2] = V->m22;
  Vm[0][1] = Vm[1][0] = V->m01;
  Vm[0][2] = Vm[2][0] = V->m02;
  Vm[1][2] = Vm[2][1] = V->m12;
  SYEV("V", "U", &n, &Vm[0][0], &n, eigs, W, &nw, &info);
#endif
  CHECK(info >= 0, "invalid argument in heev");
  CHECK(info <= 0, "failure to converge in heev");
}

/* Set Vinv = inverse of V, where both V and Vinv are real-symmetric
   (or possibly complex-Hermitian) matrices. */
void maxwell_sym_matrix_invert(symmetric_matrix *Vinv, const symmetric_matrix *V) {
  real m00 = V->m00, m11 = V->m11, m22 = V->m22;

#if defined(WITH_HERMITIAN_EPSILON)
  scalar_complex m01 = V->m01, m02 = V->m02, m12 = V->m12;

  if (m01.re == 0.0 && m01.im == 0.0 && m02.re == 0.0 && m02.im == 0.0 && m12.re == 0.0 &&
      m12.im == 0.0) {
    /* optimize common case of a diagonal matrix: */
    Vinv->m00 = 1.0 / m00;
    Vinv->m11 = 1.0 / m11;
    Vinv->m22 = 1.0 / m22;
    CASSIGN_ZERO(Vinv->m01);
    CASSIGN_ZERO(Vinv->m02);
    CASSIGN_ZERO(Vinv->m12);
  }
  else {
    double detinv;

    /* compute the determinant: */
    detinv = m00 * m11 * m22 - m11 * CSCALAR_NORMSQR(m02) - CSCALAR_NORMSQR(m01) * m22 -
             CSCALAR_NORMSQR(m12) * m00 +
             2.0 * ((m01.re * m12.re - m01.im * m12.im) * m02.re +
                    (m01.re * m12.im + m01.im * m12.re) * m02.im);

    CHECK(detinv != 0.0, "singular 3x3 matrix");

    detinv = 1.0 / detinv;

    Vinv->m00 = detinv * (m11 * m22 - CSCALAR_NORMSQR(m12));
    Vinv->m11 = detinv * (m00 * m22 - CSCALAR_NORMSQR(m02));
    Vinv->m22 = detinv * (m11 * m00 - CSCALAR_NORMSQR(m01));

    CASSIGN_SCALAR(Vinv->m02, detinv * (m01.re * m12.re - m01.im * m12.im - m11 * m02.re),
                   -detinv * (-m01.re * m12.im - m01.im * m12.re + m11 * m02.im));

    CASSIGN_SCALAR(Vinv->m01, detinv * (m12.re * m02.re + m12.im * m02.im - m22 * m01.re),
                   -detinv * (m12.im * m02.re - m12.re * m02.im + m22 * m01.im));

    CASSIGN_SCALAR(Vinv->m12, detinv * (m01.re * m02.re + m01.im * m02.im - m00 * m12.re),
                   -detinv * (m01.im * m02.re - m01.re * m02.im + m00 * m12.im));
  }

#else  /* real matrix */
  real m01 = V->m01, m02 = V->m02, m12 = V->m12;

  if (m01 == 0.0 && m02 == 0.0 && m12 == 0.0) {
    /* optimize common case of a diagonal matrix: */
    Vinv->m00 = 1.0 / m00;
    Vinv->m11 = 1.0 / m11;
    Vinv->m22 = 1.0 / m22;
    Vinv->m01 = Vinv->m02 = Vinv->m12 = 0.0;
  }
  else {
    double detinv;

    /* compute the determinant: */
    detinv = m00 * m11 * m22 - m02 * m11 * m02 + 2.0 * m01 * m12 * m02 - m01 * m01 * m22 -
             m12 * m12 * m00;

    CHECK(detinv != 0.0, "singular 3x3 matrix");

    detinv = 1.0 / detinv;

    Vinv->m00 = detinv * (m11 * m22 - m12 * m12);
    Vinv->m11 = detinv * (m00 * m22 - m02 * m02);
    Vinv->m22 = detinv * (m11 * m00 - m01 * m01);

    Vinv->m02 = detinv * (m01 * m12 - m11 * m02);
    Vinv->m01 = detinv * (m12 * m02 - m01 * m22);
    Vinv->m12 = detinv * (m01 * m02 - m00 * m12);
  }
#endif /* real matrix */
}

/* Returns whether or not V is positive-definite. */
int maxwell_sym_matrix_positive_definite(symmetric_matrix *V) {
  real det2, det3;
  real m00 = V->m00, m11 = V->m11, m22 = V->m22;

#if defined(WITH_HERMITIAN_EPSILON)
  scalar_complex m01 = V->m01, m02 = V->m02, m12 = V->m12;

  det2 = m00 * m11 - CSCALAR_NORMSQR(m01);
  det3 = det2 * m22 - m11 * CSCALAR_NORMSQR(m02) - CSCALAR_NORMSQR(m12) * m00 +
         2.0 * ((m01.re * m12.re - m01.im * m12.im) * m02.re +
                (m01.re * m12.im + m01.im * m12.re) * m02.im);
#else  /* real matrix */
  real m01 = V->m01, m02 = V->m02, m12 = V->m12;

  det2 = m00 * m11 - m01 * m01;
  det3 = det2 * m22 - m02 * m11 * m02 + 2.0 * m01 * m12 * m02 - m12 * m12 * m00;
#endif /* real matrix */

  return (m00 > 0.0 && det2 > 0.0 && det3 > 0.0);
}

#define EQ(x1, x2) (fabs((x1) - (x2)) < tol)
static int sym_matrix_eq(symmetric_matrix V1, symmetric_matrix V2, double tol) {
  if (!EQ(V1.m00, V2.m00) || !EQ(V1.m11, V2.m11) || !EQ(V1.m22, V2.m22)) return 0;
#if defined(WITH_HERMITIAN_EPSILON)
  return (EQ(V1.m01.re, V2.m01.re) && EQ(V1.m01.im, V2.m01.im) && EQ(V1.m02.re, V2.m02.re) &&
          EQ(V1.m02.im, V2.m02.im) && EQ(V1.m12.re, V2.m12.re) && EQ(V1.m12.im, V2.m12.im));
#else
  return (EQ(V1.m01, V2.m01) && EQ(V1.m02, V2.m02) && EQ(V1.m12, V2.m12));
#endif
}

/* rotate A by a unitary (real) rotation matrix R:
      RAR = transpose(R) * A * R
*/
void maxwell_sym_matrix_rotate(symmetric_matrix *RAR, const symmetric_matrix *A_, double R[3][3]) {
  int i, j;
  double A[3][3], AR[3][3];
  A[0][0] = A_->m00;
  A[1][1] = A_->m11;
  A[2][2] = A_->m22;
  A[0][1] = A[1][0] = ESCALAR_RE(A_->m01);
  A[0][2] = A[2][0] = ESCALAR_RE(A_->m02);
  A[1][2] = A[2][1] = ESCALAR_RE(A_->m12);
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      AR[i][j] = A[i][0] * R[0][j] + A[i][1] * R[1][j] + A[i][2] * R[2][j];
  for (i = 0; i < 3; ++i)
    for (j = i; j < 3; ++j)
      A[i][j] = R[0][i] * AR[0][j] + R[1][i] * AR[1][j] + R[2][i] * AR[2][j];
  RAR->m00 = A[0][0];
  RAR->m11 = A[1][1];
  RAR->m22 = A[2][2];
  ESCALAR_RE(RAR->m01) = A[0][1];
  ESCALAR_RE(RAR->m02) = A[0][2];
  ESCALAR_RE(RAR->m12) = A[1][2];
#if defined(WITH_HERMITIAN_EPSILON)
  A[0][0] = A[1][1] = A[2][2] = 0;
  A[1][0] = -(A[0][1] = ESCALAR_IM(A_->m01));
  A[2][0] = -(A[0][2] = ESCALAR_IM(A_->m02));
  A[2][1] = -(A[1][2] = ESCALAR_IM(A_->m12));
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      AR[i][j] = A[i][0] * R[0][j] + A[i][1] * R[1][j] + A[i][2] * R[2][j];
  for (i = 0; i < 3; ++i)
    for (j = i; j < 3; ++j)
      A[i][j] = R[0][i] * AR[0][j] + R[1][i] * AR[1][j] + R[2][i] * AR[2][j];
  ESCALAR_IM(RAR->m01) = A[0][1];
  ESCALAR_IM(RAR->m02) = A[0][2];
  ESCALAR_IM(RAR->m12) = A[1][2];
#endif
}

/**************************************************************************/

int check_maxwell_dielectric(maxwell_data *d, int negative_epsilon_okp) {
  int i, require_2d;

  require_2d = d->nz == 1 && (d->parity & (EVEN_Z_PARITY | ODD_Z_PARITY));

  for (i = 0; i < d->fft_output_size; ++i) {
    if (!negative_epsilon_okp && !maxwell_sym_matrix_positive_definite(d->eps_inv + i)) return 1;
    if (require_2d) {
#if defined(WITH_HERMITIAN_EPSILON)
      if (d->eps_inv[i].m02.re != 0.0 || d->eps_inv[i].m02.im != 0.0 ||
          d->eps_inv[i].m12.re != 0.0 || d->eps_inv[i].m12.im != 0.0)
        return 2;
#else  /* real matrix */
      if (d->eps_inv[i].m02 != 0.0 || d->eps_inv[i].m12 != 0.0) return 2;
#endif /* real matrix */
    }
  }
  return 0;
}

/**************************************************************************/

#define K_PI 3.141592653589793238462643383279502884197
#define SMALL 1.0e-6
#define MAX2(a, b) ((a) > (b) ? (a) : (b))
#define MIN2(a, b) ((a) < (b) ? (a) : (b))

#define MAX_MOMENT_MESH NQUAD /* max # of moment-mesh vectors */
#define MOMENT_MESH_R 0.5

/* A function to set up the mesh given the grid dimensions, mesh size,
   and lattice/reciprocal vectors.  (Any mesh sizes < 1 are taken to
   be 1.)  The returned values are:

   mesh_center: the center of the mesh, relative to the integer
                mesh coordinates; e.g. the mesh_center for a 3x3x3
                mesh is the point (1,1,1).
   mesh_prod: the product of the mesh sizes.
   moment_mesh: an array of size_moment_mesh vectors, in lattice
                coordinates, of a "spherically-symmetric" mesh of
                points centered on the origin, designed to be
                used for averaging the first moment of epsilon at
                a grid point (for finding the local surface normal).
   moment_mesh_weights: an array of size_moment_mesh weights to multiply
                        the integrand values by.  */
static void get_mesh(int nx, int ny, int nz, const int mesh_size[3], real R[3][3], real G[3][3],
                     real mesh_center[3], int *mesh_prod, real moment_mesh[MAX_MOMENT_MESH][3],
                     real moment_mesh_weights[MAX_MOMENT_MESH], int *size_moment_mesh) {
  int i, j;
  real min_diam = 1e20;
  real mesh_total[3] = {0, 0, 0};
  int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);
  real weight_sum = 0.0;

  *mesh_prod = 1;
  for (i = 0; i < 3; ++i) {
    int ms = MAX2(mesh_size[i], 1);
    mesh_center[i] = (ms - 1) * 0.5;
    *mesh_prod *= ms;
  }

  *size_moment_mesh = num_sphere_quad[rank - 1];
  for (i = 0; i < num_sphere_quad[rank - 1]; ++i) {
    for (j = 0; j < 3; ++j)
      moment_mesh[i][j] = sphere_quad[rank - 1][i][j] * MOMENT_MESH_R;
    moment_mesh_weights[i] = sphere_quad[rank - 1][i][3];
  }

  CHECK(*size_moment_mesh <= MAX_MOMENT_MESH, "yikes, moment mesh too big");

  for (i = 0; i < *size_moment_mesh; ++i)
    weight_sum += moment_mesh_weights[i];
  CHECK(fabs(weight_sum - 1.0) < SMALL, "bug, incorrect moment weights");

  /* scale the moment-mesh vectors so that the sphere has a
     diameter of 2*MOMENT_MESH_R times the diameter of the
     smallest grid direction: */

  /* first, find the minimum distance between grid points along the
     lattice directions (should we use the maximum instead?): */
  for (i = 0; i < rank; ++i) {
    real ri = R[i][0] * R[i][0] + R[i][1] * R[i][1] + R[i][2] * R[i][2];
    ri = sqrt(ri) / (i == 0 ? nx : (i == 1 ? ny : nz));
    min_diam = MIN2(min_diam, ri);
  }

  /* scale moment_mesh by this diameter: */
  for (i = 0; i < *size_moment_mesh; ++i) {
    real len = 0;
    for (j = 0; j < 3; ++j) {
      moment_mesh[i][j] *= min_diam;
      len += moment_mesh[i][j] * moment_mesh[i][j];
      mesh_total[j] += moment_mesh[i][j];
    }
    CHECK(fabs(len - min_diam * min_diam * (MOMENT_MESH_R * MOMENT_MESH_R)) < SMALL,
          "bug in get_mesh: moment_mesh vector is wrong length");
  }
  CHECK(fabs(mesh_total[0]) + fabs(mesh_total[1]) + fabs(mesh_total[2]) < SMALL,
        "bug in get_mesh: moment_mesh does not average to zero");

  /* Now, convert the moment_mesh vectors to lattice/grid coordinates;
     to do this, we multiply by the G matrix (inverse of R transposed) */
  for (i = 0; i < *size_moment_mesh; ++i) {
    real v[3];
    for (j = 0; j < 3; ++j) /* make a copy of moment_mesh[i] */
      v[j] = moment_mesh[i][j];
    for (j = 0; j < 3; ++j)
      moment_mesh[i][j] = G[j][0] * v[0] + G[j][1] * v[1] + G[j][2] * v[2];
  }
}

/**************************************************************************/

/* The following function initializes the dielectric tensor md->eps_inv,
   using the dielectric function epsilon(&eps, &eps_inv, r, epsilon_data).

   epsilon is averaged over a rectangular mesh spanning the space between
   grid points; the size of the mesh is given by mesh_size.

   R[0..2] are the spatial lattice vectors, and are used to convert
   the discretization grid into spatial coordinates (with the origin
   at the (0,0,0) grid element.

   In most places, the dielectric tensor is equal to eps_inv, but at
   dielectric interfaces it varies depending upon the polarization of
   the field (for faster convergence).  In particular, it depends upon
   the direction of the field relative to the surface normal vector,
   so we must compute the latter.  The surface normal is approximated
   by the "dipole moment" of the dielectric function over a spherical
   mesh.

   Implementation note: md->eps_inv is chosen to have dimensions matching
   the output of the FFT.  Thus, its dimensions depend upon whether we are
   doing a real or complex and serial or parallel FFT. */

void set_maxwell_dielectric(maxwell_data *md, const int mesh_size[3], real R[3][3], real G[3][3],
                            maxwell_dielectric_function epsilon,
                            maxwell_dielectric_mean_function mepsilon, void *epsilon_data) {
  real s1, s2, s3, m1, m2, m3; /* grid/mesh steps */
  real mesh_center[3];
  real moment_mesh[MAX_MOMENT_MESH][3];
  real moment_mesh_weights[MAX_MOMENT_MESH];
  real eps_inv_total = 0.0;
  int i, j, k;
  int mesh_prod;
  real mesh_prod_inv;
  int size_moment_mesh = 0;
  int n1, n2, n3;
#ifdef HAVE_MPI
  int local_n2, local_y_start, local_n3;
#endif
#ifndef SCALAR_COMPLEX
  int n_other, n_last, rank;
#endif

  n1 = md->nx;
  n2 = md->ny;
  n3 = md->nz;

  get_mesh(n1, n2, n3, mesh_size, R, G, mesh_center, &mesh_prod, moment_mesh, moment_mesh_weights,
           &size_moment_mesh);
  mesh_prod_inv = 1.0 / mesh_prod;

  s1 = 1.0 / n1;
  s2 = 1.0 / n2;
  s3 = 1.0 / n3;
  m1 = s1 / MAX2(1, mesh_size[0]);
  m2 = s2 / MAX2(1, mesh_size[1]);
  m3 = s3 / MAX2(1, mesh_size[2]);

  LOOP_XYZ(md) {
    int mi, mj, mk;
#ifdef WITH_HERMITIAN_EPSILON
    symmetric_matrix eps_mean, eps_inv_mean, eps_mean_inv;
#else
    symmetric_matrix eps_mean, eps_inv_mean, eps_mean_inv;
#endif
    real norm_len;
    real norm0, norm1, norm2;
    short means_different_p, diag_eps_p;

    {
      real r[3], normal[3];
      r[0] = i1 * s1;
      r[1] = i2 * s2;
      r[2] = i3 * s3;
      if (mepsilon &&
          mepsilon(&eps_mean, &eps_inv_mean, normal, s1, s2, s3, mesh_prod_inv, r, epsilon_data)) {

        maxwell_sym_matrix_invert(md->eps_inv + xyz_index, &eps_mean);
        goto got_eps_inv;

        norm0 = R[0][0] * normal[0] + R[1][0] * normal[1] + R[2][0] * normal[2];
        norm1 = R[0][1] * normal[0] + R[1][1] * normal[1] + R[2][1] * normal[2];
        norm2 = R[0][2] * normal[0] + R[1][2] * normal[1] + R[2][2] * normal[2];
        means_different_p = 1;
        diag_eps_p = DIAG_SYMMETRIC_MATRIX(eps_mean);
        maxwell_sym_matrix_invert(&eps_mean_inv, &eps_mean);

#if !defined(SCALAR_COMPLEX) && 0 /* check inversion symmetry */
        {
          symmetric_matrix eps_mean2, eps_inv_mean2;
          real normal2[3], r2[3], nc[3];
          r2[0] = n1 == 0 ? r[0] : 1.0 - r[0];
          r2[1] = n2 == 0 ? r[1] : 1.0 - r[1];
          r2[2] = n3 == 0 ? r[2] : 1.0 - r[2];
          CHECK(mepsilon(&eps_mean2, &eps_inv_mean2, normal2, s1, s2, s3, mesh_prod_inv, r2,
                         epsilon_data),
                "mepsilon symmetry is broken");
          CHECK(sym_matrix_eq(eps_mean, eps_mean2, 1e-10) &&
                    sym_matrix_eq(eps_inv_mean, eps_inv_mean2, 1e-10),
                "inversion symmetry is broken");
          nc[0] = normal[1] * normal2[2] - normal[2] * normal2[1];
          nc[1] = normal[2] * normal2[0] - normal[0] * normal2[2];
          nc[2] = normal[0] * normal2[1] - normal[1] * normal2[0];
          CHECK(sqrt(nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2]) < 1e-6,
                "normal-vector symmetry is broken");
        }
#endif

        goto got_mean;
      }
    }

    eps_mean.m00 = eps_mean.m11 = eps_mean.m22 = eps_inv_mean.m00 = eps_inv_mean.m11 =
        eps_inv_mean.m22 = 0.0;
    ASSIGN_ESCALAR(eps_mean.m01, 0, 0);
    ASSIGN_ESCALAR(eps_mean.m02, 0, 0);
    ASSIGN_ESCALAR(eps_mean.m12, 0, 0);
    ASSIGN_ESCALAR(eps_inv_mean.m01, 0, 0);
    ASSIGN_ESCALAR(eps_inv_mean.m02, 0, 0);
    ASSIGN_ESCALAR(eps_inv_mean.m12, 0, 0);

    for (mi = 0; mi < mesh_size[0]; ++mi)
      for (mj = 0; mj < mesh_size[1]; ++mj)
        for (mk = 0; mk < mesh_size[2]; ++mk) {
          real r[3];
          symmetric_matrix eps, eps_inv;
          r[0] = i1 * s1 + (mi - mesh_center[0]) * m1;
          r[1] = i2 * s2 + (mj - mesh_center[1]) * m2;
          r[2] = i3 * s3 + (mk - mesh_center[2]) * m3;
          epsilon(&eps, &eps_inv, r, epsilon_data);
          eps_mean.m00 += eps.m00;
          eps_mean.m11 += eps.m11;
          eps_mean.m22 += eps.m22;
          eps_inv_mean.m00 += eps_inv.m00;
          eps_inv_mean.m11 += eps_inv.m11;
          eps_inv_mean.m22 += eps_inv.m22;
#ifdef WITH_HERMITIAN_EPSILON
          CACCUMULATE_SUM(eps_mean.m01, eps.m01);
          CACCUMULATE_SUM(eps_mean.m02, eps.m02);
          CACCUMULATE_SUM(eps_mean.m12, eps.m12);
          CACCUMULATE_SUM(eps_inv_mean.m01, eps_inv.m01);
          CACCUMULATE_SUM(eps_inv_mean.m02, eps_inv.m02);
          CACCUMULATE_SUM(eps_inv_mean.m12, eps_inv.m12);
#else
          eps_mean.m01 += eps.m01;
          eps_mean.m02 += eps.m02;
          eps_mean.m12 += eps.m12;
          eps_inv_mean.m01 += eps_inv.m01;
          eps_inv_mean.m02 += eps_inv.m02;
          eps_inv_mean.m12 += eps_inv.m12;
#endif
        }

    diag_eps_p = DIAG_SYMMETRIC_MATRIX(eps_mean);
    if (diag_eps_p) { /* handle the common case of diagonal matrices: */
      eps_mean_inv.m00 = mesh_prod / eps_mean.m00;
      eps_mean_inv.m11 = mesh_prod / eps_mean.m11;
      eps_mean_inv.m22 = mesh_prod / eps_mean.m22;
#ifdef WITH_HERMITIAN_EPSILON
      CASSIGN_ZERO(eps_mean_inv.m01);
      CASSIGN_ZERO(eps_mean_inv.m02);
      CASSIGN_ZERO(eps_mean_inv.m12);
#else
      eps_mean_inv.m01 = eps_mean_inv.m02 = eps_mean_inv.m12 = 0.0;
#endif
      eps_inv_mean.m00 *= mesh_prod_inv;
      eps_inv_mean.m11 *= mesh_prod_inv;
      eps_inv_mean.m22 *= mesh_prod_inv;

      means_different_p = fabs(eps_mean_inv.m00 - eps_inv_mean.m00) > SMALL ||
                          fabs(eps_mean_inv.m11 - eps_inv_mean.m11) > SMALL ||
                          fabs(eps_mean_inv.m22 - eps_inv_mean.m22) > SMALL;
    }
    else {
      eps_inv_mean.m00 *= mesh_prod_inv;
      eps_inv_mean.m11 *= mesh_prod_inv;
      eps_inv_mean.m22 *= mesh_prod_inv;
      eps_mean.m00 *= mesh_prod_inv;
      eps_mean.m11 *= mesh_prod_inv;
      eps_mean.m22 *= mesh_prod_inv;
#ifdef WITH_HERMITIAN_EPSILON
      eps_mean.m01.re *= mesh_prod_inv;
      eps_mean.m01.im *= mesh_prod_inv;
      eps_mean.m02.re *= mesh_prod_inv;
      eps_mean.m02.im *= mesh_prod_inv;
      eps_mean.m12.re *= mesh_prod_inv;
      eps_mean.m12.im *= mesh_prod_inv;
      eps_inv_mean.m01.re *= mesh_prod_inv;
      eps_inv_mean.m01.im *= mesh_prod_inv;
      eps_inv_mean.m02.re *= mesh_prod_inv;
      eps_inv_mean.m02.im *= mesh_prod_inv;
      eps_inv_mean.m12.re *= mesh_prod_inv;
      eps_inv_mean.m12.im *= mesh_prod_inv;
#else
      eps_mean.m01 *= mesh_prod_inv;
      eps_mean.m02 *= mesh_prod_inv;
      eps_mean.m12 *= mesh_prod_inv;
      eps_inv_mean.m01 *= mesh_prod_inv;
      eps_inv_mean.m02 *= mesh_prod_inv;
      eps_inv_mean.m12 *= mesh_prod_inv;
#endif
      maxwell_sym_matrix_invert(&eps_mean_inv, &eps_mean);

      means_different_p = fabs(eps_mean_inv.m00 - eps_inv_mean.m00) > SMALL ||
                          fabs(eps_mean_inv.m11 - eps_inv_mean.m11) > SMALL ||
                          fabs(eps_mean_inv.m22 - eps_inv_mean.m22) > SMALL;
#ifdef WITH_HERMITIAN_EPSILON
      means_different_p = means_different_p ||
                          fabs(eps_mean_inv.m01.re - eps_inv_mean.m01.re) > SMALL ||
                          fabs(eps_mean_inv.m02.re - eps_inv_mean.m02.re) > SMALL ||
                          fabs(eps_mean_inv.m12.re - eps_inv_mean.m12.re) > SMALL ||
                          fabs(eps_mean_inv.m01.im - eps_inv_mean.m01.im) > SMALL ||
                          fabs(eps_mean_inv.m02.im - eps_inv_mean.m02.im) > SMALL ||
                          fabs(eps_mean_inv.m12.im - eps_inv_mean.m12.im) > SMALL;
#else
      means_different_p = means_different_p || fabs(eps_mean_inv.m01 - eps_inv_mean.m01) > SMALL ||
                          fabs(eps_mean_inv.m02 - eps_inv_mean.m02) > SMALL ||
                          fabs(eps_mean_inv.m12 - eps_inv_mean.m12) > SMALL;
#endif
    }

    /* if the two averaging methods yielded different results,
       which usually happens if epsilon is not constant, then
       we need to find the normal vector to the dielectric interface: */
    if (means_different_p) {
      real moment0 = 0, moment1 = 0, moment2 = 0;

      for (mi = 0; mi < size_moment_mesh; ++mi) {
        real r[3], eps_trace;
        symmetric_matrix eps, eps_inv;
        r[0] = i1 * s1 + moment_mesh[mi][0];
        r[1] = i2 * s2 + moment_mesh[mi][1];
        r[2] = i3 * s3 + moment_mesh[mi][2];
        epsilon(&eps, &eps_inv, r, epsilon_data);
        eps_trace = eps.m00 + eps.m11 + eps.m22;
        eps_trace *= moment_mesh_weights[mi];
        moment0 += eps_trace * moment_mesh[mi][0];
        moment1 += eps_trace * moment_mesh[mi][1];
        moment2 += eps_trace * moment_mesh[mi][2];
      }

      /* need to convert moment from lattice to cartesian coords: */
      norm0 = R[0][0] * moment0 + R[1][0] * moment1 + R[2][0] * moment2;
      norm1 = R[0][1] * moment0 + R[1][1] * moment1 + R[2][1] * moment2;
      norm2 = R[0][2] * moment0 + R[1][2] * moment1 + R[2][2] * moment2;

    got_mean:

      norm_len = sqrt(norm0 * norm0 + norm1 * norm1 + norm2 * norm2);
    }

    if (means_different_p && norm_len > SMALL) {
      real x0, x1, x2;

      norm_len = 1.0 / norm_len;
      norm0 *= norm_len;
      norm1 *= norm_len;
      norm2 *= norm_len;

      /* Compute the effective inverse dielectric tensor.
         We define this as:
            1/2 ( {eps_inv_mean, P} + {eps_mean_inv, 1-P} )
         where P is the projection matrix onto the normal direction
         (P = norm ^ norm), and {a,b} is the anti-commutator ab+ba.
          = 1/2 {eps_inv_mean - eps_mean_inv, P} + eps_mean_inv
          = 1/2 (n_i conj(x_j) + x_i n_j) + (eps_mean_inv)_ij
         where n_k is the kth component of the normal vector and
            x_i = (eps_inv_mean - eps_mean_inv)_ik n_k
         Note the implied summations (Einstein notation).

         Note that the resulting matrix is symmetric, and we get just
         eps_inv_mean if eps_inv_mean == eps_mean_inv, as desired.

         Note that P is idempotent, so for scalar epsilon this
         is just eps_inv_mean * P + eps_mean_inv * (1-P)
               = (1/eps_inv_mean * P + eps_mean * (1-P)) ^ (-1),
         which corresponds to the expression in the Meade paper. */

      x0 = (eps_inv_mean.m00 - eps_mean_inv.m00) * norm0;
      x1 = (eps_inv_mean.m11 - eps_mean_inv.m11) * norm1;
      x2 = (eps_inv_mean.m22 - eps_mean_inv.m22) * norm2;
      if (diag_eps_p) {
#ifdef WITH_HERMITIAN_EPSILON
        md->eps_inv[xyz_index].m01.re = 0.5 * (x0 * norm1 + x1 * norm0);
        md->eps_inv[xyz_index].m01.im = 0.0;
        md->eps_inv[xyz_index].m02.re = 0.5 * (x0 * norm2 + x2 * norm0);
        md->eps_inv[xyz_index].m02.im = 0.0;
        md->eps_inv[xyz_index].m12.re = 0.5 * (x1 * norm2 + x2 * norm1);
        md->eps_inv[xyz_index].m12.im = 0.0;
#else
        md->eps_inv[xyz_index].m01 = 0.5 * (x0 * norm1 + x1 * norm0);
        md->eps_inv[xyz_index].m02 = 0.5 * (x0 * norm2 + x2 * norm0);
        md->eps_inv[xyz_index].m12 = 0.5 * (x1 * norm2 + x2 * norm1);
#endif
      }
      else {
#ifdef WITH_HERMITIAN_EPSILON
        real x0i, x1i, x2i;
        x0 += ((eps_inv_mean.m01.re - eps_mean_inv.m01.re) * norm1 +
               (eps_inv_mean.m02.re - eps_mean_inv.m02.re) * norm2);
        x1 += ((eps_inv_mean.m01.re - eps_mean_inv.m01.re) * norm0 +
               (eps_inv_mean.m12.re - eps_mean_inv.m12.re) * norm2);
        x2 += ((eps_inv_mean.m02.re - eps_mean_inv.m02.re) * norm0 +
               (eps_inv_mean.m12.re - eps_mean_inv.m12.re) * norm1);
        x0i = ((eps_inv_mean.m01.im - eps_mean_inv.m01.im) * norm1 +
               (eps_inv_mean.m02.im - eps_mean_inv.m02.im) * norm2);
        x1i = (-(eps_inv_mean.m01.im - eps_mean_inv.m01.im) * norm0 +
               (eps_inv_mean.m12.im - eps_mean_inv.m12.im) * norm2);
        x2i = -((eps_inv_mean.m02.im - eps_mean_inv.m02.im) * norm0 +
                (eps_inv_mean.m12.im - eps_mean_inv.m12.im) * norm1);

        md->eps_inv[xyz_index].m01.re = (0.5 * (x0 * norm1 + x1 * norm0) + eps_mean_inv.m01.re);
        md->eps_inv[xyz_index].m02.re = (0.5 * (x0 * norm2 + x2 * norm0) + eps_mean_inv.m02.re);
        md->eps_inv[xyz_index].m12.re = (0.5 * (x1 * norm2 + x2 * norm1) + eps_mean_inv.m12.re);
        md->eps_inv[xyz_index].m01.im = (0.5 * (x0i * norm1 - x1i * norm0) + eps_mean_inv.m01.im);
        md->eps_inv[xyz_index].m02.im = (0.5 * (x0i * norm2 - x2i * norm0) + eps_mean_inv.m02.im);
        md->eps_inv[xyz_index].m12.im = (0.5 * (x1i * norm2 - x2i * norm1) + eps_mean_inv.m12.im);
#else
        x0 += ((eps_inv_mean.m01 - eps_mean_inv.m01) * norm1 +
               (eps_inv_mean.m02 - eps_mean_inv.m02) * norm2);
        x1 += ((eps_inv_mean.m01 - eps_mean_inv.m01) * norm0 +
               (eps_inv_mean.m12 - eps_mean_inv.m12) * norm2);
        x2 += ((eps_inv_mean.m02 - eps_mean_inv.m02) * norm0 +
               (eps_inv_mean.m12 - eps_mean_inv.m12) * norm1);

        md->eps_inv[xyz_index].m01 = (0.5 * (x0 * norm1 + x1 * norm0) + eps_mean_inv.m01);
        md->eps_inv[xyz_index].m02 = (0.5 * (x0 * norm2 + x2 * norm0) + eps_mean_inv.m02);
        md->eps_inv[xyz_index].m12 = (0.5 * (x1 * norm2 + x2 * norm1) + eps_mean_inv.m12);
#endif
      }
      md->eps_inv[xyz_index].m00 = x0 * norm0 + eps_mean_inv.m00;
      md->eps_inv[xyz_index].m11 = x1 * norm1 + eps_mean_inv.m11;
      md->eps_inv[xyz_index].m22 = x2 * norm2 + eps_mean_inv.m22;
    }
    else { /* undetermined normal vector and/or constant eps */
      md->eps_inv[xyz_index] = eps_mean_inv;
    }
  got_eps_inv:

    eps_inv_total +=
        (md->eps_inv[xyz_index].m00 + md->eps_inv[xyz_index].m11 + md->eps_inv[xyz_index].m22);
  }
}
} /* end of loop body */

mpi_allreduce_1(&eps_inv_total, real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
n1 = md->fft_output_size;
mpi_allreduce_1(&n1, int, MPI_INT, MPI_SUM, mpb_comm);
md->eps_inv_mean = eps_inv_total / (3 * n1);
}

void set_maxwell_mu(maxwell_data *md, const int mesh_size[3], real R[3][3], real G[3][3],
                    maxwell_dielectric_function mu, maxwell_dielectric_mean_function mmu,
                    void *mu_data) {
  symmetric_matrix *eps_inv = md->eps_inv;
  real eps_inv_mean = md->eps_inv_mean;
  if (md->mu_inv == NULL) { CHK_MALLOC(md->mu_inv, symmetric_matrix, md->fft_output_size); }
  /* just re-use code to set epsilon, but initialize mu_inv instead */
  md->eps_inv = md->mu_inv;
  set_maxwell_dielectric(md, mesh_size, R, G, mu, mmu, mu_data);
  md->eps_inv = eps_inv;
  md->mu_inv_mean = md->eps_inv_mean;
  md->eps_inv_mean = eps_inv_mean;
}
