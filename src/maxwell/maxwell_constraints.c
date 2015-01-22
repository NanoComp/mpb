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
#include <math.h>

#include "config.h"
#include <check.h>

#include <mpiglue.h>
#include "maxwell.h"

/**************************************************************************/

/* function to call z and y parity constraints, if necessary */
void maxwell_parity_constraint(evectmatrix X, void *data)
{
     maxwell_data *d = (maxwell_data *) data;

     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");

     if (d->parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
	  maxwell_zparity_constraint(X, data);
     if (d->parity & (EVEN_Y_PARITY | ODD_Y_PARITY))
	  maxwell_yparity_constraint(X, data);
}

/**************************************************************************/

/* In 3d dielectric structures having a z=0 mirror plane (symmetric
   under z -> -z), the states will exhibit an analogue of TM and TE
   polarizations (when k has no z-component to break the symmetry).

   In this case, the states can be classified as "even" or "odd" with
   respect to mirror-flips through z=0.  This is referred to as the
   "parity" of the state, where even is parity +1 and odd is parity -1
   (the state is an eigenvector of the mirror flip operator with this
   eigenvalue).  Even/odd states are the analogues of TE/TM states,
   respectively (and in the mirror plane itself they are truly TE/TM
   polarized).

   Note that the magnetic field is a pseudo-vector, so the mirror
   operation acts specially on it.  Also, because of the way the
   m,n transverse basis for H is chosen, the basis vectors transform
   in a very simple way (just flip sign). */

/* Project X to its even or odd component, so that we can solve
   for only one parity of states (the projection operator, like the
   mirror flip operator, commutes with the Maxwell operator, so this
   projection should not slow convergence).  */
void maxwell_zparity_constraint(evectmatrix X, void *data)
{
     maxwell_data *d = (maxwell_data *) data;
     int i, j, b, nxy, nz;
     int zparity = ((d->parity & EVEN_Z_PARITY) ? +1 :
		    ((d->parity & ODD_Z_PARITY) ? -1 : 0));
     
     if (zparity == 0)
	  return;

     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");

     if (d->nz > 1) {
	  nxy = d->other_dims;
	  nz = d->last_dim;
     }
     else {  /* common case (2d system): even/odd == TE/TM */
	  nxy = d->other_dims * d->last_dim;
	  if (zparity == +1)
	       for (i = 0; i < nxy; ++i) 
		    for (b = 0; b < X.p; ++b) {
			 ASSIGN_ZERO(X.data[(i * X.c + 1) * X.p + b]);
		    }
	  else if (zparity == -1)
	       for (i = 0; i < nxy; ++i) 
		    for (b = 0; b < X.p; ++b) {
			 ASSIGN_ZERO(X.data[(i * X.c) * X.p + b]);
		    }
	  return;
     }

     for (i = 0; i < nxy; ++i) {
	  for (j = 0; 2*j <= nz; ++j) {
	       int ij = i * nz + j; 
	       int ij2 = i * nz + (j > 0 ? nz - j : 0);
	       for (b = 0; b < X.p; ++b) {
		    scalar u,v, u2,v2;
		    u = X.data[(ij * 2) * X.p + b];
		    v = X.data[(ij * 2 + 1) * X.p + b];
		    u2 = X.data[(ij2 * 2) * X.p + b];
		    v2 = X.data[(ij2 * 2 + 1) * X.p + b];
		    ASSIGN_SCALAR(X.data[(ij * 2) * X.p + b],
				  0.5*(SCALAR_RE(u) + zparity*SCALAR_RE(u2)),
				  0.5*(SCALAR_IM(u) + zparity*SCALAR_IM(u2)));
		    ASSIGN_SCALAR(X.data[(ij * 2 + 1) * X.p + b],
				  0.5*(SCALAR_RE(v) - zparity*SCALAR_RE(v2)),
				  0.5*(SCALAR_IM(v) - zparity*SCALAR_IM(v2)));
		    ASSIGN_SCALAR(X.data[(ij2 * 2) * X.p + b],
				  0.5*(SCALAR_RE(u2) + zparity*SCALAR_RE(u)),
				  0.5*(SCALAR_IM(u2) + zparity*SCALAR_IM(u)));
		    ASSIGN_SCALAR(X.data[(ij2 * 2 + 1) * X.p + b],
				  0.5*(SCALAR_RE(v2) - zparity*SCALAR_RE(v)),
				  0.5*(SCALAR_IM(v2) - zparity*SCALAR_IM(v)));
	       }
	  }
     }
}

/* Compute the parity of all of the states in X, returning an array
   of the parities (which the caller should deallocate with free).
   The parity of an arbitrary state is defined as the expectation value
   of the mirror flip operator, and will be +1/-1 for even/odd eigenstates
   and something in between for everything else.  Assumes that the
   columns of X are normalized to 1. */
double *maxwell_zparity(evectmatrix X, maxwell_data *d)
{
     int i, j, b, nxy, nz;
     double *zparity, *zp_scratch, *norm_scratch;

     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");

     CHK_MALLOC(zparity, double, X.p);
     CHK_MALLOC(zp_scratch, double, X.p);
     for (b = 0; b < X.p; ++b)
	  zp_scratch[b] = 0.0;
     CHK_MALLOC(norm_scratch, double, X.p);
     for (b = 0; b < X.p; ++b)
	  norm_scratch[b] = 0.0;

     if (d->nz > 1) {
	  nxy = d->other_dims;
	  nz = d->last_dim;
     }
     else {
	  nxy = d->other_dims * d->last_dim;
	  nz = 1;
     }

     for (i = 0; i < nxy; ++i)
	  for (j = 0; 2*j <= nz; ++j) {
	       int ij = i * nz + j; 
	       int ij2 = i * nz + (j > 0 ? nz - j : 0);
	       for (b = 0; b < X.p; ++b) {
		    scalar u,v, u2,v2;
		    u = X.data[(ij * 2) * X.p + b];
		    v = X.data[(ij * 2 + 1) * X.p + b];
		    u2 = X.data[(ij2 * 2) * X.p + b];
		    v2 = X.data[(ij2 * 2 + 1) * X.p + b];
		    zp_scratch[b] += (ij == ij2 ? 1.0 : 2.0) *
			 (SCALAR_RE(u) * SCALAR_RE(u2) +
			  SCALAR_IM(u) * SCALAR_IM(u2) -
			  SCALAR_RE(v) * SCALAR_RE(v2) -
			  SCALAR_IM(v) * SCALAR_IM(v2));
		    norm_scratch[b] += (ij == ij2 ? 1.0 : 2.0) *
                         (SCALAR_RE(u) * SCALAR_RE(u) +
                          SCALAR_IM(u) * SCALAR_IM(u) +
                          SCALAR_RE(v) * SCALAR_RE(v) +
                          SCALAR_IM(v) * SCALAR_IM(v));
	       }
	  }

     mpi_allreduce(zp_scratch, zparity, X.p,
		   double, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce(norm_scratch, zp_scratch, X.p,
		   double, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     
     for (b = 0; b < X.p; ++b)
         zparity[b] /= zp_scratch[b];
     free(zp_scratch);
     free(norm_scratch);
     
     return zparity;
}

/**************************************************************************/

/* Similar to the zparity functions above, but for the y -> -y mirror flip. */

/* Project X to its even or odd component, so that we can solve
   for only one parity of states (the projection operator, like the
   mirror flip operator, commutes with the Maxwell operator, so this
   projection should not slow convergence).  */
void maxwell_yparity_constraint(evectmatrix X, void *data)
{
     maxwell_data *d = (maxwell_data *) data;
     int i, j, k, b, nx, ny, nz;
     int yparity = ((d->parity & EVEN_Y_PARITY) ? +1 :
		    ((d->parity & ODD_Y_PARITY) ? -1 : 0));

     if (yparity == 0)
	  return;

     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");

     nx = d->local_nx;
     ny = d->ny;
     nz = d->nz;

     for (i = 0; i < nx; ++i) {
	  for (j = 0; 2*j <= ny; ++j) {
	       int ij = i * ny + j; 
	       int ij2 = i * ny + (j > 0 ? ny - j : 0);
	       for (k = 0; k < nz; ++k) {
		    int ijk = ij * nz + k;
		    int ijk2 = ij2 * nz + k;
		    for (b = 0; b < X.p; ++b) {
			 scalar u,v, u2,v2;
			 u = X.data[(ijk * 2) * X.p + b];
			 v = X.data[(ijk * 2 + 1) * X.p + b];
			 u2 = X.data[(ijk2 * 2) * X.p + b];
			 v2 = X.data[(ijk2 * 2 + 1) * X.p + b];
			 ASSIGN_SCALAR(X.data[(ijk * 2) * X.p + b],
				  0.5*(SCALAR_RE(u) - yparity*SCALAR_RE(u2)),
				  0.5*(SCALAR_IM(u) - yparity*SCALAR_IM(u2)));
			 ASSIGN_SCALAR(X.data[(ijk * 2 + 1) * X.p + b],
				  0.5*(SCALAR_RE(v) + yparity*SCALAR_RE(v2)),
				  0.5*(SCALAR_IM(v) + yparity*SCALAR_IM(v2)));
			 ASSIGN_SCALAR(X.data[(ijk2 * 2) * X.p + b],
				  0.5*(SCALAR_RE(u2) - yparity*SCALAR_RE(u)),
				  0.5*(SCALAR_IM(u2) - yparity*SCALAR_IM(u)));
			 ASSIGN_SCALAR(X.data[(ijk2 * 2 + 1) * X.p + b],
				  0.5*(SCALAR_RE(v2) + yparity*SCALAR_RE(v)),
				  0.5*(SCALAR_IM(v2) + yparity*SCALAR_IM(v)));
		    }
	       }
	  }
     }
}

/* Compute the parity of all of the states in X, returning an array
   of the parities (which the caller should deallocate with free).
   The parity of an arbitrary state is defined as the expectation value
   of the mirror flip operator, and will be +1/-1 for even/odd eigenstates
   and something in between for everything else.  Assumes that the
   columns of X are normalized to 1. */
double *maxwell_yparity(evectmatrix X, maxwell_data *d)
{
     int i, j, k, b, nx, ny, nz;
     double *yparity, *yp_scratch, *norm_scratch;

     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");

     CHK_MALLOC(yparity, double, X.p);
     CHK_MALLOC(yp_scratch, double, X.p);
     for (b = 0; b < X.p; ++b)
	  yp_scratch[b] = 0.0;
     CHK_MALLOC(norm_scratch, double, X.p);
     for (b = 0; b < X.p; ++b)
	  norm_scratch[b] = 0.0;

     nx = d->local_nx;
     ny = d->ny;
     nz = d->nz;

     for (i = 0; i < nx; ++i) {
	  for (j = 0; 2*j <= ny; ++j) {
	       int ij = i * ny + j; 
	       int ij2 = i * ny + (j > 0 ? ny - j : 0);
	       for (k = 0; k < nz; ++k) {
		    int ijk = ij * nz + k;
		    int ijk2 = ij2 * nz + k;
		    for (b = 0; b < X.p; ++b) {
			 scalar u,v, u2,v2;
			 u = X.data[(ijk * 2) * X.p + b];
			 v = X.data[(ijk * 2 + 1) * X.p + b];
			 u2 = X.data[(ijk2 * 2) * X.p + b];
			 v2 = X.data[(ijk2 * 2 + 1) * X.p + b];
			 yp_scratch[b] += (ijk == ijk2 ? 1.0 : 2.0) *
			      (SCALAR_RE(v) * SCALAR_RE(v2) +
			       SCALAR_IM(v) * SCALAR_IM(v2) -
			       SCALAR_RE(u) * SCALAR_RE(u2) -
			       SCALAR_IM(u) * SCALAR_IM(u2));
			 norm_scratch[b] += (ijk == ijk2 ? 1.0 : 2.0) *
			      (SCALAR_RE(v) * SCALAR_RE(v) +
			       SCALAR_IM(v) * SCALAR_IM(v) +
			       SCALAR_RE(u) * SCALAR_RE(u) +
			       SCALAR_IM(u) * SCALAR_IM(u));
		    }
	       }
	  }
     }

     mpi_allreduce(yp_scratch, yparity, X.p,
		   double, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     mpi_allreduce(norm_scratch, yp_scratch, X.p,
		   double, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     
     for (b = 0; b < X.p; ++b)
         yparity[b] /= yp_scratch[b];
     free(yp_scratch);
     free(norm_scratch);
     
     return yparity;
}

/**************************************************************************/

/* to fix problems with slow convergence for k ~ 0, manually "put in"
   the k = 0 solution: first two bands are constant and higher bands are
   orthogonal.  Note that in the TE/TM case, only one band is constant. 
   Also note that, in Fourier space, a constant field corresponds to
   1 in the DC component and 0 elsewhere. */

/* return the number of constant (zero-frequency) bands: */
int maxwell_zero_k_num_const_bands(evectmatrix X, maxwell_data *d)
{
     int num_const_bands, m_band = 1, n_band = 1;
     
     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");
     
     if (d->parity & (ODD_Z_PARITY | EVEN_Y_PARITY))
	  m_band = 0;
     if (d->parity & (ODD_Y_PARITY | EVEN_Z_PARITY))
	  n_band = 0;
     
     num_const_bands = m_band + n_band;
     
     if (num_const_bands > X.p)
	  num_const_bands = X.p;
     
     return num_const_bands;
}

void maxwell_zero_k_set_const_bands(evectmatrix X, maxwell_data *d)
{
     int i, j, num_const_bands, m_band = 1, n_band = 1;
     
     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");

     if (X.p < 1)
	  return;

     num_const_bands = maxwell_zero_k_num_const_bands(X, d);

     /* Initialize num_const_bands to zero: */
     for (i = 0; i < X.n; ++i) 
	  for (j = 0; j < num_const_bands; ++j) {
	       ASSIGN_ZERO(X.data[i * X.p + j]);
	  }
     
     if (X.Nstart > 0)
	  return;  /* DC frequency is not on this process */
		      
     /* Set DC components to 1 (in two parities) for num_const_bands: */

     if (d->parity & (ODD_Z_PARITY | EVEN_Y_PARITY))
	  m_band = 0;
     if (d->parity & (ODD_Y_PARITY | EVEN_Z_PARITY))
	  n_band = 0;

     if (m_band) {
	  ASSIGN_SCALAR(X.data[0], 1.0, 0.0);
	  ASSIGN_SCALAR(X.data[X.p], 0.0, 0.0);
     }
     if (n_band && (!m_band || X.p >= 2)) {
	  ASSIGN_SCALAR(X.data[m_band], 0.0, 0.0);
	  ASSIGN_SCALAR(X.data[X.p + m_band], 1.0, 0.0);
     }
}

/* during eigensolution (for upper bands), their DC components are
   constrained to be zero */
void maxwell_zero_k_constraint(evectmatrix X, void *data)
{
     int j;

     if (X.Nstart > 0)
	  return;  /* DC frequency is not on this process */
		      
     for (j = 0; j < X.p; ++j) {
	  ASSIGN_ZERO(X.data[j]);
	  ASSIGN_ZERO(X.data[X.p + j]);
     }
     (void)data; /* avoid warning about unused parameter */
}

/**************************************************************************/

