/* Copyright (C) 1999 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../config.h"
#include <check.h>

#include "maxwell.h"

#define PRECOND_SUBTR_EIGS 0

#define PRECOND_MIN_DENOM 0.01

void maxwell_preconditioner(evectmatrix Xin, evectmatrix Xout, void *data,
			    evectmatrix Y, real *eigenvals)
{
     maxwell_data *d = (maxwell_data *) data;
     int i, c, b;
     real *kpGn2 = d->k_plus_G_normsqr;

     for (i = 0; i < Xout.localN; ++i) {
	  for (c = 0; c < Xout.c; ++c) {
	       for (b = 0; b < Xout.p; ++b) {
		    int index = (i * Xout.c + c) * Xout.p + b;
		    real scale = kpGn2[i] * d->eps_inv_mean;

#if PRECOND_SUBTR_EIGS
		    scale -= eigenvals[b];
		    scale = 1.0 / (scale + copysign(PRECOND_MIN_DENOM, scale));
#else
		    scale = 1.0 / (scale + PRECOND_MIN_DENOM);
#endif

		    ASSIGN_SCALAR(Xout.data[index],
				  scale * SCALAR_RE(Xin.data[index]),
				  scale * SCALAR_IM(Xin.data[index]));
	       }
	  }
     }
}

void maxwell_target_preconditioner(evectmatrix Xin, evectmatrix Xout, 
				   void *data,
				   evectmatrix Y, real *eigenvals)
{
     maxwell_target_data *td = (maxwell_target_data *) data;
     maxwell_data *d = td->d;
     real omega_sqr = td->target_frequency * td->target_frequency;
     int i, c, b;
     real *kpGn2 = d->k_plus_G_normsqr;

     for (i = 0; i < Xout.localN; ++i) {
	  for (c = 0; c < Xout.c; ++c) {
	       for (b = 0; b < Xout.p; ++b) {
		    int index = (i * Xout.c + c) * Xout.p + b;
		    real scale = kpGn2[i] * d->eps_inv_mean - omega_sqr;

		    scale = scale * scale;

#if PRECOND_SUBTR_EIGS
		    scale -= eigenvals[b];
		    scale = 1.0 / (scale + copysign(PRECOND_MIN_DENOM, scale));
#else
		    scale = 1.0 / (scale + PRECOND_MIN_DENOM);
#endif

		    ASSIGN_SCALAR(Xout.data[index],
				  scale * SCALAR_RE(Xin.data[index]),
				  scale * SCALAR_IM(Xin.data[index]));
	       }
	  }
     }
}

void maxwell_constraint(evectmatrix X, void *data)
{
     maxwell_data *d = (maxwell_data *) data;
     int i, j;

     CHECK(d, "null maxwell data pointer!");
     CHECK(X.c == 2, "fields don't have 2 components!");

     /* Enforce the polarization of the current eigenvector.
        Note that this constraint is preserved by the maxwell
        operator (assuming it's two-dimensional), so it shouldn't
        screw up convergence. */
     
     if (d->polarization == TE_POLARIZATION)
	  for (i = 0; i < X.localN; ++i) 
	       for (j = 0; j < X.p; ++j) {
		    ASSIGN_ZERO(X.data[(i * X.c + 1) * X.p + j]);
	       }
     else if (d->polarization == TM_POLARIZATION)
	  for (i = 0; i < X.localN; ++i) 
	       for (j = 0; j < X.p; ++j) {
		    ASSIGN_ZERO(X.data[(i * X.c) * X.p + j]);
	       }
}
