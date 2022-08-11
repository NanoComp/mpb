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
#  define HEEV F77_FUNC(cheev,CHEEV)
#  define SYEV F77_FUNC(ssyev,SSYEV)
#else
#  define HEEV F77_FUNC(zheev,ZHEEV)
#  define SYEV F77_FUNC(dsyev,DSYEV)
#endif
extern void HEEV(char *, char *, int *, scalar_complex *, int *, real *,
		 scalar_complex *, int *, real *, int *);
extern void SYEV(char *, char *, int *, real *, int *, real *,
		 real *, int *, int *);

/* compute the 3 real eigenvalues of the matrix V */
void maxwell_sym_matrix_eigs(real eigs[3], const symmetric_matrix *V)
{
     int n = 3, nw = 9, info;
#if defined(WITH_HERMITIAN_EPSILON)
     scalar_complex Vm[3][3], W[9];
     real W2[9];
     CASSIGN_SCALAR(Vm[0][0], V->m00, 0);
     CASSIGN_SCALAR(Vm[1][1], V->m11, 0);
     CASSIGN_SCALAR(Vm[2][2], V->m22, 0);
     Vm[0][1] = V->m01; CASSIGN_CONJ(Vm[1][0], V->m01);
     Vm[0][2] = V->m02; CASSIGN_CONJ(Vm[2][0], V->m02);
     Vm[1][2] = V->m12; CASSIGN_CONJ(Vm[2][1], V->m12);
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
void maxwell_sym_matrix_invert(symmetric_matrix *Vinv,
			       const symmetric_matrix *V)
{
     real m00 = V->m00, m11 = V->m11, m22 = V->m22;

#if defined(WITH_HERMITIAN_EPSILON)
     scalar_complex m01 = V->m01, m02 = V->m02, m12 = V->m12;

     if (m01.re == 0.0 && m01.im == 0.0 &&
	 m02.re == 0.0 && m02.im == 0.0 &&
	 m12.re == 0.0 && m12.im == 0.0) {
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
	  detinv = m00*m11*m22 - m11*CSCALAR_NORMSQR(m02) -
	       CSCALAR_NORMSQR(m01)*m22 - CSCALAR_NORMSQR(m12)*m00 +
	       2.0 * ((m01.re * m12.re - m01.im * m12.im) * m02.re +
		      (m01.re * m12.im + m01.im * m12.re) * m02.im);

	  CHECK(detinv != 0.0, "singular 3x3 matrix");

	  detinv = 1.0/detinv;

	  Vinv->m00 = detinv * (m11*m22 - CSCALAR_NORMSQR(m12));
	  Vinv->m11 = detinv * (m00*m22 - CSCALAR_NORMSQR(m02));
	  Vinv->m22 = detinv * (m11*m00 - CSCALAR_NORMSQR(m01));

	  CASSIGN_SCALAR(Vinv->m02,
			 detinv * (m01.re*m12.re-m01.im*m12.im - m11*m02.re),
			 -detinv*(-m01.re*m12.im-m01.im*m12.re + m11*m02.im));

	  CASSIGN_SCALAR(Vinv->m01,
			 detinv * (m12.re*m02.re+m12.im*m02.im - m22*m01.re),
			 -detinv * (m12.im*m02.re-m12.re*m02.im + m22*m01.im));

	  CASSIGN_SCALAR(Vinv->m12,
			 detinv * (m01.re*m02.re+m01.im*m02.im - m00*m12.re),
			 -detinv * (m01.im*m02.re-m01.re*m02.im + m00*m12.im));
     }

#else /* real matrix */
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
	  detinv = m00*m11*m22 - m02*m11*m02 + 2.0 * m01*m12*m02 -
	       m01*m01*m22 - m12*m12*m00;

	  CHECK(detinv != 0.0, "singular 3x3 matrix");

	  detinv = 1.0/detinv;

	  Vinv->m00 = detinv * (m11*m22 - m12*m12);
	  Vinv->m11 = detinv * (m00*m22 - m02*m02);
	  Vinv->m22 = detinv * (m11*m00 - m01*m01);

	  Vinv->m02 = detinv * (m01*m12 - m11*m02);
	  Vinv->m01 = detinv * (m12*m02 - m01*m22);
	  Vinv->m12 = detinv * (m01*m02 - m00*m12);
     }
#endif /* real matrix */
}

/* Returns whether or not V is positive-definite. */
int maxwell_sym_matrix_positive_definite(symmetric_matrix *V)
{
     real det2, det3;
     real m00 = V->m00, m11 = V->m11, m22 = V->m22;

#if defined(WITH_HERMITIAN_EPSILON)
     scalar_complex m01 = V->m01, m02 = V->m02, m12 = V->m12;

     det2 = m00*m11 - CSCALAR_NORMSQR(m01);
     det3 = det2*m22 - m11*CSCALAR_NORMSQR(m02) - CSCALAR_NORMSQR(m12)*m00 +
	  2.0 * ((m01.re * m12.re - m01.im * m12.im) * m02.re +
		 (m01.re * m12.im + m01.im * m12.re) * m02.im);
#else /* real matrix */
     real m01 = V->m01, m02 = V->m02, m12 = V->m12;

     det2 = m00*m11 - m01*m01;
     det3 = det2*m22 - m02*m11*m02 + 2.0 * m01*m12*m02 - m12*m12*m00;
#endif /* real matrix */

     return (m00 > 0.0 && det2 > 0.0 && det3 > 0.0);
}

#define EQ(x1,x2) (fabs((x1) - (x2)) < tol)
static int sym_matrix_eq(symmetric_matrix V1, symmetric_matrix V2, double tol)
{
     if (!EQ(V1.m00,V2.m00) || !EQ(V1.m11,V2.m11) || !EQ(V1.m22,V2.m22))
	  return 0;
#if defined(WITH_HERMITIAN_EPSILON)
     return(EQ(V1.m01.re,V2.m01.re) && EQ(V1.m01.im,V2.m01.im) &&
	    EQ(V1.m02.re,V2.m02.re) && EQ(V1.m02.im,V2.m02.im) &&
	    EQ(V1.m12.re,V2.m12.re) && EQ(V1.m12.im,V2.m12.im));
#else
     return(EQ(V1.m01,V2.m01) && EQ(V1.m02,V2.m02) && EQ(V1.m12,V2.m12));
#endif
}

/* rotate A by a unitary (real) rotation matrix R:
      RAR = transpose(R) * A * R
*/
void maxwell_sym_matrix_rotate(symmetric_matrix *RAR,
			       const symmetric_matrix *A_,
			       double R[3][3])
{
     int i,j;
     double A[3][3], AR[3][3];
     A[0][0] = A_->m00;
     A[1][1] = A_->m11;
     A[2][2] = A_->m22;
     A[0][1] = A[1][0] = ESCALAR_RE(A_->m01);
     A[0][2] = A[2][0] = ESCALAR_RE(A_->m02);
     A[1][2] = A[2][1] = ESCALAR_RE(A_->m12);
     for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j)
	  AR[i][j] = A[i][0]*R[0][j] + A[i][1]*R[1][j] + A[i][2]*R[2][j];
     for (i = 0; i < 3; ++i) for (j = i; j < 3; ++j)
	  A[i][j] = R[0][i]*AR[0][j] + R[1][i]*AR[1][j] + R[2][i]*AR[2][j];
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
     for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j)
	  AR[i][j] = A[i][0]*R[0][j] + A[i][1]*R[1][j] + A[i][2]*R[2][j];
     for (i = 0; i < 3; ++i) for (j = i; j < 3; ++j)
	  A[i][j] = R[0][i]*AR[0][j] + R[1][i]*AR[1][j] + R[2][i]*AR[2][j];
     ESCALAR_IM(RAR->m01) = A[0][1];
     ESCALAR_IM(RAR->m02) = A[0][2];
     ESCALAR_IM(RAR->m12) = A[1][2];
#endif
}

/* compute a rotation matrix Rot that rotates to a coordinate system whose first
   axis lies along the (unit-length) vector [n0, n1, n2]                          */
void maxwell_rotation_matrix(double Rot[3][3], double n0, double n1, double n2)
{
   Rot[0][0] = n0;
   Rot[1][0] = n1;
   Rot[2][0] = n2;
   if (fabs(n0) > 1e-2 || fabs(n1) > 1e-2) { /* (z x n) */
        Rot[0][2] = n1;
        Rot[1][2] = -n0;
        Rot[2][2] = 0;
   }
   else { /* n is ~ parallel to z direction, use (x x n) instead */
        Rot[0][2] = 0;
        Rot[1][2] = -n2;
        Rot[2][2] = n1;
   }
   { /* normalize second column */
        double s = Rot[0][2]*Rot[0][2]+Rot[1][2]*Rot[1][2]+Rot[2][2]*Rot[2][2];
        s = 1.0 / sqrt(s);
        Rot[0][2] *= s;
        Rot[1][2] *= s;
        Rot[2][2] *= s;
   }
   /* 1st column is 2nd column x 0th column */
   Rot[0][1] = Rot[1][2] * Rot[2][0] - Rot[2][2] * Rot[1][0];
   Rot[1][1] = Rot[2][2] * Rot[0][0] - Rot[0][2] * Rot[2][0];
   Rot[2][1] = Rot[0][2] * Rot[1][0] - Rot[1][2] * Rot[0][0];
}

/**************************************************************************/

int check_maxwell_dielectric(maxwell_data *d,
			     int negative_epsilon_okp)
{
     int i, require_2d;

     require_2d = d->nz == 1 && (d->parity & (EVEN_Z_PARITY | ODD_Z_PARITY));

     for (i = 0; i < d->fft_output_size; ++i) {
	  if (!negative_epsilon_okp &&
	      !maxwell_sym_matrix_positive_definite(d->eps_inv + i))
	       return 1;
	  if (require_2d) {
#if defined(WITH_HERMITIAN_EPSILON)
	       if (d->eps_inv[i].m02.re != 0.0 ||
		   d->eps_inv[i].m02.im != 0.0 ||
		   d->eps_inv[i].m12.re != 0.0 ||
		   d->eps_inv[i].m12.im != 0.0)
		    return 2;
#else /* real matrix */
	       if (d->eps_inv[i].m02 != 0.0 || d->eps_inv[i].m12 != 0.0)
		    return 2;
#endif /* real matrix */
	  }
     }
     return 0;
}

/**************************************************************************/

#define K_PI 3.141592653589793238462643383279502884197
#define SMALL 1.0e-6
#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))

#define MAX_MOMENT_MESH NQUAD /* max # of moment-mesh vectors */
#define MOMENT_MESH_R_1D 0.5                /* diagonal/2 of a D-dimensional voxel */
#define MOMENT_MESH_R_2D 0.7071067811865475 /* sqrt(1/2) */
#define MOMENT_MESH_R_3D 0.8660254037844386 /* sqrt(3/4) */

/* Function to set up a voxel averaging mesh given the mesh size (any mesh
   sizes < 1 are taken to be 1). The returned values are:

   mesh_center: the center of the mesh, relative to the integer
                mesh coordinates; e.g. the mesh_center for a 3x3x3
                mesh is the point (1,1,1).
   mesh_prod: the product of the mesh sizes.                             */
static void get_mesh(const int mesh_size[3], real mesh_center[3], int *mesh_prod)
{
     int i;

     *mesh_prod = 1;
     for (i = 0; i < 3; ++i) {
	    int ms = MAX2(mesh_size[i], 1);
	    mesh_center[i] = (ms - 1) * 0.5;
	    *mesh_prod *= ms;
     }
}

/* Function to set up a spherical mesh for moment averaging (to compute normal vectors)
   given the grid dimensions nx, ny, nz and and the lattice & reciprocal vectors R & G.
   The returned values are:

   moment_mesh: an array of size_moment_mesh vectors, in lattice
                coordinates, of a "spherically-symmetric" mesh of
                points centered on the origin, designed to be
                used for averaging the first moment of epsilon at
                a grid point (for finding the local surface normal).
   moment_mesh_weights: an array of size_moment_mesh weights to multiply
                the integrand values by.
   size_moment_mesh: number of integration points (2 in 1D, NUMQUAD2=12
                in 2D, NUMQUAD3=50 in 3D)                                    */
static void get_moment_mesh(int nx, int ny, int nz,
                  real R[3][3], real G[3][3],
                  real moment_mesh[MAX_MOMENT_MESH][3],
                  real moment_mesh_weights[MAX_MOMENT_MESH],
                  int *size_moment_mesh)
{
     int i,j;
     real max_diam = 0;
     real mesh_total[3] = { 0, 0, 0 };
     int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);
     real weight_sum = 0.0;
     real moment_mesh_r = rank == 1 ? MOMENT_MESH_R_1D :
                         (rank == 2 ? MOMENT_MESH_R_2D : MOMENT_MESH_R_3D);

     *size_moment_mesh = num_sphere_quad[rank-1];
     for (i = 0; i < num_sphere_quad[rank-1]; ++i) {
	  for (j = 0; j < 3; ++j)
	       moment_mesh[i][j] = sphere_quad[rank-1][i][j] * moment_mesh_r;
	  moment_mesh_weights[i] = sphere_quad[rank-1][i][3];
     }

     CHECK(*size_moment_mesh <= MAX_MOMENT_MESH, "yikes, moment mesh too big");

     for (i = 0; i < *size_moment_mesh; ++i)
	  weight_sum += moment_mesh_weights[i];
     CHECK(fabs(weight_sum - 1.0) < SMALL, "bug, incorrect moment weights");

     /* scale the moment-mesh vectors so that the sphere has a diameter of
        2*moment_mesh_r times the diameter of the largest grid direction
        (to ensure the spherical mesh encloses the entire voxel):          */

     /* first, find the maximum distance between grid points along the
        lattice directions: */
     for (i = 0; i < rank; ++i) {
	  real ri = R[i][0] * R[i][0] + R[i][1] * R[i][1] + R[i][2] * R[i][2];
	  ri = sqrt(ri) / (i == 0 ? nx : (i == 1 ? ny : nz));
	  max_diam = MAX2(max_diam, ri);
     }

     /* scale moment_mesh by this diameter: */
     for (i = 0; i < *size_moment_mesh; ++i) {
	  real len = 0;
	  for (j = 0; j < 3; ++j) {
	       moment_mesh[i][j] *= max_diam;
	       len += moment_mesh[i][j] * moment_mesh[i][j];
	       mesh_total[j] += moment_mesh[i][j];
	  }
	  CHECK(fabs(len - max_diam*max_diam*(moment_mesh_r*moment_mesh_r))
		< SMALL,
		"bug in get_moment_mesh: moment_mesh vector is wrong length");
     }
     CHECK(fabs(mesh_total[0]) + fabs(mesh_total[1]) + fabs(mesh_total[2])
	   < SMALL, "bug in get_moment_mesh: moment_mesh does not average to zero");

     /* Now, convert the moment_mesh vectors to lattice/grid coordinates;
        to do this, we multiply by the G matrix (inverse of R transposed) */
     for (i = 0; i < *size_moment_mesh; ++i) {
	  real v[3];
	  for (j = 0; j < 3; ++j)    /* make a copy of moment_mesh[i] */
	       v[j] = moment_mesh[i][j];
	  for (j = 0; j < 3; ++j)
	       moment_mesh[i][j] = G[j][0]*v[0] + G[j][1]*v[1] + G[j][2]*v[2];
     }
}

static symmetric_matrix average_eps_inv_over_mesh(real r[3],
               maxwell_dielectric_function epsilon,
               const int mesh_size[3],
               real mesh_center[3],
               real m1, real m2, real m3,
               real mesh_prod_inv,
               void *epsilon_data) {
   symmetric_matrix eps_inv_mean;
   int mi, mj, mk;

   eps_inv_mean.m00 = eps_inv_mean.m11 = eps_inv_mean.m22 = 0.0;
   ASSIGN_ESCALAR(eps_inv_mean.m01, 0, 0);
   ASSIGN_ESCALAR(eps_inv_mean.m02, 0, 0);
   ASSIGN_ESCALAR(eps_inv_mean.m12, 0, 0);

   /* sum up eps_inv over mesh points */
   for (mi = 0; mi < mesh_size[0]; ++mi) {
      for (mj = 0; mj < mesh_size[1]; ++mj) {
         for (mk = 0; mk < mesh_size[2]; ++mk) {
            real r_mesh[3];
            symmetric_matrix eps, eps_inv;
            r_mesh[0] = r[0] + (mi - mesh_center[0]) * m1;
            r_mesh[1] = r[1] + (mj - mesh_center[1]) * m2;
            r_mesh[2] = r[2] + (mk - mesh_center[2]) * m3;
            epsilon(&eps, &eps_inv, r_mesh, epsilon_data);
            eps_inv_mean.m00 += eps_inv.m00;
            eps_inv_mean.m11 += eps_inv.m11;
            eps_inv_mean.m22 += eps_inv.m22;
            EACCUMULATE_SUM(eps_inv_mean.m01, eps_inv.m01);
            EACCUMULATE_SUM(eps_inv_mean.m02, eps_inv.m02);
            EACCUMULATE_SUM(eps_inv_mean.m12, eps_inv.m12);
         }
      }
   }

   /* rescale to get average */
   eps_inv_mean.m00 *= mesh_prod_inv;
   eps_inv_mean.m11 *= mesh_prod_inv;
   eps_inv_mean.m22 *= mesh_prod_inv;
   /* if not a diagonal matrix, rescale off-diagonal entries as well */
   if (!DIAG_SYMMETRIC_MATRIX(eps_inv_mean)) {
#ifdef WITH_HERMITIAN_EPSILON
      eps_inv_mean.m01.re *= mesh_prod_inv;
      eps_inv_mean.m01.im *= mesh_prod_inv;
      eps_inv_mean.m02.re *= mesh_prod_inv;
      eps_inv_mean.m02.im *= mesh_prod_inv;
      eps_inv_mean.m12.re *= mesh_prod_inv;
      eps_inv_mean.m12.im *= mesh_prod_inv;
#else
      eps_inv_mean.m01 *= mesh_prod_inv;
      eps_inv_mean.m02 *= mesh_prod_inv;
      eps_inv_mean.m12 *= mesh_prod_inv;
#endif
   }

   return eps_inv_mean;
}

/* Function to detect whether a voxel centered at `r` with sides `s1`, `s2`, & `s3`
   intersects a material interface (returns 1) or not (returns 0).
   Intersection is based on whether or not the permittivity tensor evaluated at the
   voxel corners deviates significantly (`> SMALL`) from its value at the voxel center.
   The approach is analogous to that taken in `mean_epsilon_func(...)`.                 */
short detect_interface_via_corner_check(real r[3],
            maxwell_dielectric_function epsilon,
            int rank,
            real s1, real s2, real s3,
            void *epsilon_data) {
   short is_interface;
   int i;
   const int num_corners[3] = {2, 4, 8};
   const real corners[3][8][3] = {
            { {-0.5,0,0}, {0.5,0,0},
            {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0} },
            { {-0.5,-0.5,0}, {0.5,0.5,0}, {-0.5,0.5,0}, {0.5,-0.5,0},
            {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0} },
            { {0.5,0.5,0.5}, {0.5,0.5,-0.5}, {0.5,-0.5,0.5}, {0.5,-0.5,-0.5},
            {-0.5,0.5,0.5}, {-0.5,0.5,-0.5}, {-0.5,-0.5,0.5}, {-0.5,-0.5,-0.5} } };
   symmetric_matrix eps, eps_inv, eps_corner, eps_inv_corner;
   real r_corner[3];

   epsilon(&eps, &eps_inv, r, epsilon_data); /* eps at voxel center */
   for (i = 0; i < num_corners[rank - 1]; ++i) {
      r_corner[0] = r[0] + corners[rank - 1][i][0] * s1;
      r_corner[1] = r[1] + corners[rank - 1][i][1] * s2;
      r_corner[2] = r[2] + corners[rank - 1][i][2] * s3;
      epsilon(&eps_corner, &eps_inv_corner, r_corner, epsilon_data);

      /* define the presence of an interface from whether an epsilon component at
         a voxel corner varies by more than SMALL from the voxel center's epsilon */
      is_interface = fabs(eps.m00 - eps_corner.m00) > SMALL ||
                     fabs(eps.m11 - eps_corner.m11) > SMALL ||
                     fabs(eps.m22 - eps_corner.m22) > SMALL ||
#ifdef WITH_HERMITIAN_EPSILON
                     fabs(eps.m01.re - eps_corner.m01.re) > SMALL ||
                     fabs(eps.m01.im - eps_corner.m01.im) > SMALL ||
                     fabs(eps.m02.re - eps_corner.m02.re) > SMALL ||
                     fabs(eps.m02.im - eps_corner.m02.im) > SMALL ||
                     fabs(eps.m12.re - eps_corner.m12.re) > SMALL ||
                     fabs(eps.m12.im - eps_corner.m12.im) > SMALL;
#else
                     fabs(eps.m01 - eps_corner.m01) > SMALL ||
                     fabs(eps.m02 - eps_corner.m02) > SMALL ||
                     fabs(eps.m12 - eps_corner.m12) > SMALL;
#endif
      if (is_interface)
         return is_interface;
   }

   return is_interface;
}

/**************************************************************************/

/* The following function initializes the dielectric tensor md->eps_inv,
   using the dielectric function epsilon(&eps, &eps_inv, r, epsilon_data).

   epsilon is (Kottke) averaged over a rectangular mesh spanning the space
   between grid points; the size of the mesh is given by mesh_size.

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

void set_maxwell_dielectric(maxwell_data *md,
			    const int mesh_size[3],
			    real R[3][3], real G[3][3],
			    maxwell_dielectric_function epsilon,
			    maxwell_dielectric_mean_function mepsilon,
			    void *epsilon_data)
{
     real s1, s2, s3, m1, m2, m3;  /* grid/mesh steps */
     real mesh_center[3];
     real moment_mesh[MAX_MOMENT_MESH][3];
     real moment_mesh_weights[MAX_MOMENT_MESH];
     real eps_inv_total = 0.0;
     int mesh_prod;
     real mesh_prod_inv;
     int size_moment_mesh = 0;
     int n1 = md->nx, n2 = md->ny, n3 = md->nz;
     int rank = n3 > 1 ? 3 : (n2 > 1 ? 2 : 1);
     short is_interface;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
#ifndef SCALAR_COMPLEX
     int n_other, n_last;
#endif

     /* integration mesh for checking whether voxel intersects an interface */
     get_mesh(mesh_size, mesh_center, &mesh_prod);
     mesh_prod_inv = 1.0 / mesh_prod;

     s1 = 1.0 / n1;
     s2 = 1.0 / n2;
     s3 = 1.0 / n3;
     m1 = s1 / MAX2(1, mesh_size[0]);
     m2 = s2 / MAX2(1, mesh_size[1]);
     m3 = s3 / MAX2(1, mesh_size[2]);

     /* spherical integration mesh for computing normal vectors */
     get_moment_mesh(n1, n2, n3, R, G, moment_mesh, moment_mesh_weights, &size_moment_mesh);

     LOOP_XYZ(md) {
	     int mi, mj, mk;
	     symmetric_matrix eps_mean, eps_inv_mean;
	     real norm_len;
	     real norm0, norm1, norm2;
        real r[3], normal[3];

	     r[0] = i1 * s1;
	     r[1] = i2 * s2;
	     r[2] = i3 * s3;

       if (epsilon && mesh_prod == 1) {
         epsilon(&eps_mean, &eps_inv_mean, r, epsilon_data);
         md->eps_inv[xyz_index] = eps_inv_mean;
         goto got_eps_inv;
       }

	     if (mepsilon && mepsilon(&eps_mean, &eps_inv_mean, normal,
					s1, s2, s3, mesh_prod_inv,
					r, epsilon_data)) {

		     maxwell_sym_matrix_invert(md->eps_inv + xyz_index, &eps_mean);

		     goto got_eps_inv;

#if !defined(SCALAR_COMPLEX) && 0 /* check inversion symmetry */
		    {
			 symmetric_matrix eps_mean2, eps_inv_mean2;
			 real normal2[3], r2[3], nc[3];
			 r2[0] = n1 == 0 ? r[0] : 1.0 - r[0];
			 r2[1] = n2 == 0 ? r[1] : 1.0 - r[1];
			 r2[2] = n3 == 0 ? r[2] : 1.0 - r[2];
			 CHECK(mepsilon(&eps_mean2, &eps_inv_mean2, normal2,
					s1, s2, s3, mesh_prod_inv,
					r2, epsilon_data),
			       "mepsilon symmetry is broken");
			 CHECK(sym_matrix_eq(eps_mean,eps_mean2,1e-10) &&
			       sym_matrix_eq(eps_inv_mean,eps_inv_mean2,1e-10),
			       "inversion symmetry is broken");
			 nc[0] = normal[1]*normal2[2] - normal[2]*normal2[1];
			 nc[1] = normal[2]*normal2[0] - normal[0]*normal2[2];
			 nc[2] = normal[0]*normal2[1] - normal[1]*normal2[0];
			 CHECK(sqrt(nc[0]*nc[0]+nc[1]*nc[1]+nc[2]*nc[2])<1e-6,
			       "normal-vector symmetry is broken");
		    }
#endif
     }

     /* detect whether voxel intersects an interface */
     is_interface = detect_interface_via_corner_check(r, epsilon, rank,
                                                      s1, s2, s3, epsilon_data);

     /* if an interface was detected, then we need to find the normal vector to
        the dielectric interface for Kottke averaging: */
     if (is_interface) {
	       real moment0 = 0, moment1 = 0, moment2 = 0;

	       for (mi = 0; mi < size_moment_mesh; ++mi) {
		       real r_mesh[3], eps_trace;
		       symmetric_matrix eps, eps_inv;
		       r_mesh[0] = r[0] + moment_mesh[mi][0];
		       r_mesh[1] = r[1] + moment_mesh[mi][1];
		       r_mesh[2] = r[2] + moment_mesh[mi][2];
		       epsilon(&eps, &eps_inv, r_mesh, epsilon_data);
		       eps_trace = eps.m00 + eps.m11 + eps.m22;
		       eps_trace *= moment_mesh_weights[mi];
		       moment0 += eps_trace * moment_mesh[mi][0];
		       moment1 += eps_trace * moment_mesh[mi][1];
		       moment2 += eps_trace * moment_mesh[mi][2];
	       }

	       /* need to convert moment from lattice to cartesian coords: */
	       norm0 = R[0][0]*moment0 + R[1][0]*moment1 + R[2][0]*moment2;
	       norm1 = R[0][1]*moment0 + R[1][1]*moment1 + R[2][1]*moment2;
	       norm2 = R[0][2]*moment0 + R[1][2]*moment1 + R[2][2]*moment2;

	       norm_len = sqrt(norm0*norm0 + norm1*norm1 + norm2*norm2);
	  }

      /* === Kottke averaging of interface voxels === */
      if (is_interface && norm_len > SMALL) {
           double Rot[3][3];
           symmetric_matrix tau;

           /* --- compute mean[τ(ε)] --- */
           /* τ(ε) is defined by [Kottke PRE, Eq. (4)]:
                       ( -1/ε₁₁   ε₁₂/ε₁₁         ε₁₃/ε₁₁        )
                τ(ε) = ( ε₂₁/ε₁₁  ε₂₂-ε₂₁ε₁₂/ε₁₁  ε₂₃-ε₂₁ε₁₃/ε₁₁ )
                       ( ε₃₁/ε₁₁  ε₃₂-ε₃₁ε₁₂/ε₁₁  ε₃₃-ε₃₁ε₁₃/ε₁₁ )
              where subscripts refer to an orthogonal coordinate system where 1-hat
              is orthogonal to the interface, and 2-hat & 3-hat are parallel.
              Thus, at each position in the averaging step, we need to rotate the
              cartesian ε appropriately (would be good to do without allocations,
              but don't bother to at the moment).                                    */

           /* compute rotation matrix to interface coordinate system */
           norm_len = 1.0/norm_len;
           norm0 *= norm_len;
           norm1 *= norm_len;
           norm2 *= norm_len;
           maxwell_rotation_matrix(Rot, norm0, norm1, norm2);

           /* sum over voxel mesh points, accumulating mesh_prod*mean[τ(ε)] */
           tau.m00 = tau.m11 = tau.m22 = 0.0;
           ASSIGN_ESCALAR(tau.m01, 0.0, 0.0);
           ASSIGN_ESCALAR(tau.m02, 0.0, 0.0);
           ASSIGN_ESCALAR(tau.m12, 0.0, 0.0);
           for (mi = 0; mi < mesh_size[0]; ++mi) {
              for (mj = 0; mj < mesh_size[1]; ++mj) {
                 for (mk = 0; mk < mesh_size[2]; ++mk) {
                    real r_mesh[3];
                    symmetric_matrix eps, eps_inv, teps;
                    r_mesh[0] = r[0] + (mi - mesh_center[0]) * m1;
                    r_mesh[1] = r[1] + (mj - mesh_center[1]) * m2;
                    r_mesh[2] = r[2] + (mk - mesh_center[2]) * m3;
                    epsilon(&eps, &eps_inv, r_mesh, epsilon_data); /* tensor in cartesian system */

                    /* rotate epsilon tensor to interface coordinate system */
                    maxwell_sym_matrix_rotate(&teps, &eps, Rot);

                    /* integrate τ(ε) */
                    tau.m00 += -1.0/teps.m00;                                                /* -1/ε₁₁         */
                    tau.m11 += teps.m11 - ESCALAR_NORMSQR(teps.m01)/teps.m00;                /* ε₂₂-ε₂₁ε₁₂/ε₁₁ */
                    tau.m22 += teps.m22 - ESCALAR_NORMSQR(teps.m02)/teps.m00;                /* ε₃₃-ε₃₁ε₁₃/ε₁₁ */

#ifdef WITH_HERMITIAN_EPSILON
                    CACCUMULATE_SCALAR(tau.m01, teps.m01.re/teps.m00, teps.m01.im/teps.m00); /* ε₁₂/ε₁₁        */
                    CACCUMULATE_SCALAR(tau.m02, teps.m02.re/teps.m00, teps.m02.im/teps.m00); /* ε₁₃/ε₁₁        */
                    CACCUMULATE_SCALAR(tau.m12,                                              /* ε₂₃-ε₂₁ε₁₃/ε₁₁ */
                        teps.m12.re - CSCALAR_MULT_CONJ_RE(teps.m02, teps.m01)/teps.m00,
                        teps.m12.im - CSCALAR_MULT_CONJ_IM(teps.m02, teps.m01)/teps.m00 );

#else
                    tau.m01 += teps.m01/teps.m00;                                            /* ε₁₂/ε₁₁        */
                    tau.m02 += teps.m02/teps.m00;                                            /* ε₁₃/ε₁₁        */
                    tau.m12 += teps.m12 - teps.m01*teps.m02/teps.m00;                        /* ε₂₃-ε₂₁ε₁₃/ε₁₁ */
#endif
                 }
              }
           }

           /* normalize τ-summation to get mean */
           tau.m00 *= mesh_prod_inv;
           tau.m11 *= mesh_prod_inv;
           tau.m22 *= mesh_prod_inv;
#ifdef WITH_HERMITIAN_EPSILON
           tau.m01.re *= mesh_prod_inv;
           tau.m02.re *= mesh_prod_inv;
           tau.m12.re *= mesh_prod_inv;
           tau.m01.im *= mesh_prod_inv;
           tau.m02.im *= mesh_prod_inv;
           tau.m12.im *= mesh_prod_inv;
#else
           tau.m01 *= mesh_prod_inv;
           tau.m02 *= mesh_prod_inv;
           tau.m12 *= mesh_prod_inv;
#endif

           /* --- compute τ⁻¹[mean(τ(ε))] (i.e. the Kottke-averaged permittivity) --- */
           /* τ⁻¹(τ) is defined by [Kottke PRE, Eq. (23)]:
                         ( -1/τ₁₁    -τ₁₂/τ₁₁        -τ₁₃/τ₁₁       )
                τ⁻¹(τ) = ( -τ₂₁/τ₁₁  τ₂₂-τ₂₁τ₁₂/τ₁₁  τ₂₃-τ₂₁τ₁₃/τ₁₁ ) = meanᴷ[ε]
                         ( -τ₃₁/τ₁₁  τ₃₂-τ₃₁τ₁₂/τ₁₁  τ₃₃-τ₃₁τ₁₃/τ₁₁ )
           */
           eps_mean.m00 = -1/tau.m00;
           eps_mean.m11 = tau.m11 - ESCALAR_NORMSQR(tau.m01)/tau.m00;
           eps_mean.m22 = tau.m22 - ESCALAR_NORMSQR(tau.m02)/tau.m00;

           ASSIGN_ESCALAR(eps_mean.m01,
                    -ESCALAR_RE(tau.m01)/tau.m00, -ESCALAR_IM(tau.m01)/tau.m00);
           ASSIGN_ESCALAR(eps_mean.m02,
                    -ESCALAR_RE(tau.m02)/tau.m00, -ESCALAR_IM(tau.m02)/tau.m00);
           ASSIGN_ESCALAR(eps_mean.m12,
                    ESCALAR_RE(tau.m12) - ESCALAR_MULT_CONJ_RE(tau.m02, tau.m01)/tau.m00,
                    ESCALAR_IM(tau.m12) - ESCALAR_MULT_CONJ_IM(tau.m02, tau.m01)/tau.m00 );

           /* --- rotate eps_mean (i.e. τ⁻¹) back to the cartesian coordinate system --- */
#define SWAP(a,b) { double xxx = a; a = b; b = xxx; }  /* invert via tranposition */
           SWAP(Rot[0][1], Rot[1][0]);
           SWAP(Rot[0][2], Rot[2][0]);
           SWAP(Rot[2][1], Rot[1][2]);
           maxwell_sym_matrix_rotate(&eps_mean, &eps_mean, Rot);
#undef SWAP

           /* invert eps_mean to get the Kottke averaged inverse permittivity */
            maxwell_sym_matrix_invert(md->eps_inv + xyz_index, &eps_mean);

	  }
	  else { /* undetermined normal vector and/or constant eps */
        if (is_interface) { /* if normal is nearly zero but an interface was detected
                               fall back to ordinary averaging (e.g. for smooth epsilon) */
           eps_inv_mean = average_eps_inv_over_mesh(r, epsilon, mesh_size, mesh_center,
                                       m1, m2, m3, mesh_prod_inv, epsilon_data);
        }
        else {              /* otherwise, assume constant epsilon in voxel */
           epsilon(&eps_mean, &eps_inv_mean, r, epsilon_data);
        }
        md->eps_inv[xyz_index] = eps_inv_mean;
	  }
     got_eps_inv:

	  eps_inv_total += (md->eps_inv[xyz_index].m00 +
			    md->eps_inv[xyz_index].m11 +
			    md->eps_inv[xyz_index].m22);
     }}}  /* end of loop body */

     mpi_allreduce_1(&eps_inv_total, real, SCALAR_MPI_TYPE,
		     MPI_SUM, mpb_comm);
     n1 = md->fft_output_size;
     mpi_allreduce_1(&n1, int, MPI_INT, MPI_SUM, mpb_comm);
     md->eps_inv_mean = eps_inv_total / (3 * n1);
}

void set_maxwell_mu(maxwell_data *md,
                    const int mesh_size[3],
                    real R[3][3], real G[3][3],
                    maxwell_dielectric_function mu,
                    maxwell_dielectric_mean_function mmu,
                    void *mu_data) {
    symmetric_matrix *eps_inv = md->eps_inv;
    real eps_inv_mean = md->eps_inv_mean;
    if (md->mu_inv == NULL) {
        CHK_MALLOC(md->mu_inv, symmetric_matrix, md->fft_output_size);
    }
    /* just re-use code to set epsilon, but initialize mu_inv instead */
    md->eps_inv = md->mu_inv;
    set_maxwell_dielectric(md, mesh_size, R, G, mu, mmu, mu_data);
    md->eps_inv = eps_inv;
    md->mu_inv_mean = md->eps_inv_mean;
    md->eps_inv_mean = eps_inv_mean;
}
