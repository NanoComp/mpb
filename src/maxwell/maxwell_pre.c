#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <config.h>
#include <check.h>

#include "maxwell.h"

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
		    real scale = kpGn2[i] * d->eps_inv_mean[b];
/*         real scale = kpGn2[i] * d->eps_inv_mean[b] - eigenvals[b]; */

		    scale = 1.0 / (scale + 0.01);
/*		    scale = 1.0 / (scale + copysign(0.01, scale));  */
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
