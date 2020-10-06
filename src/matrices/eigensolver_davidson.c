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

/* This file contains an alternative eigensolver, currently experimental,
   based on the Davidson method (a preconditioned variant of Lanczos):

   M. Crouzeix, B. Philippe, and M. Sadkane, "The Davidson Method,"
   SIAM J. Sci. Comput. 15, no. 1, pp. 62-76 (January 1994). */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <scalar.h>
#include <matrices.h>
#include <blasglue.h>

#include "eigensolver.h"
#include "verbosity.h"

extern void eigensolver_get_eigenvals_aux(evectmatrix Y, real *eigenvals,
                                          evectoperator A, void *Adata,
                                          evectmatrix Work1, evectmatrix Work2,
                                          sqmatrix U, sqmatrix Usqrt,
                                          sqmatrix Uwork);

#define STRINGIZEx(x) #x /* a hack so that we can stringize macro values */
#define STRINGIZE(x) STRINGIZEx(x)

/**************************************************************************/

#define EIGENSOLVER_MAX_ITERATIONS 100000
#define FEEDBACK_TIME 4.0 /* elapsed time before we print progress feedback */

/**************************************************************************/

void eigensolver_davidson(evectmatrix Y, real *eigenvals,
			  evectoperator A, void *Adata,
			  evectpreconditioner K, void *Kdata,
			  evectconstraint constraint, void *constraint_data,
			  evectmatrix Work[], int nWork,
			  real tolerance, int *num_iterations,
			  int flags,
			  real target)
{
     int nbasis, q;
     evectmatrix *AV, *V;
     sqmatrix VAV, S, Swork, U, S2, S3, I;
     mpiglue_clock_t prev_feedback_time;
     int iteration = 0, ibasis = 0;
     real *eigenvals2, prev_E = 0;

     prev_feedback_time = MPIGLUE_CLOCK;

#ifdef DEBUG
     flags |= EIGS_VERBOSE;
#endif

     CHECK(nWork >= 4, "not enough workspace");

     nbasis = nWork / 2;
     V = Work;
     AV = Work + nbasis;

     q = Y.p * nbasis;
     VAV = create_sqmatrix(q);
     S = create_sqmatrix(q);
     Swork = create_sqmatrix(q);

     sqmatrix_resize(&VAV, 0, 0);
     sqmatrix_resize(&S, 0, 0);
     sqmatrix_resize(&Swork, 0, 0);

     CHK_MALLOC(eigenvals2, real, q);

     U = create_sqmatrix(Y.p);
     S2 = create_sqmatrix(Y.p);
     S3 = create_sqmatrix(Y.p);

     I = create_sqmatrix(0);

     if (constraint)
	  constraint(Y, constraint_data);

     evectmatrix_XtX(U, Y, S3);
     CHECK(sqmatrix_invert(U, 1, S3), "singular YtY at start");
     sqmatrix_sqrt(S2, U, S3); /* S2 = 1/sqrt(Yt*Y) */
     evectmatrix_XeYS(V[0], Y, S2, 1); /* V[0] = orthonormalize Y */

     do {
	  real E;
	  int itarget, i;

	  A(V[ibasis], AV[ibasis], Adata, 0, Y);

	  q = Y.p * (ibasis + 1);
	  sqmatrix_resize(&VAV, q, 1);
	  sqmatrix_resize(&S, q, 0);
	  sqmatrix_resize(&Swork, q, 0);

	  for (i = 0; i <= ibasis; ++i) {
	       evectmatrixXtY_sub(VAV, Y.p * (q * i + ibasis),
				  V[i], AV[ibasis], S3);
	  }
	  sqmatrix_copy_upper2full(S, VAV);

	  sqmatrix_eigensolve(S, eigenvals2, Swork);

	  /* find index itarget of start of "window" around the
	     target frequency : */
	  if (target == 0.0) /* not attempting targeted eigensolver */
	       itarget = 0;
	  else {
	       /* note that this technique seems to have convergence trouble */
	       for (itarget = 0; itarget + Y.p < q &&
			 fabs(target - eigenvals2[itarget]) >
			 fabs(target - eigenvals2[itarget + Y.p]); ++itarget)
		    ;
	  }

	  for (E = 0.0, i = 0; i < Y.p; ++i) {
	       E += (eigenvals[i] = eigenvals2[itarget + i]);
	  }
	  mpi_assert_equal(E);

	  /* compute Y = best eigenvectors */
	  for (i = 0; i <= ibasis; ++i) {
	       evectmatrix_aXpbYS_sub(i ? 1.0 : 0.0, Y,
				      1.0, V[i],
				      S, itarget * q + Y.p * i, 1);
	  }

	  if (iteration > 0 && mpi_is_master() &&
	      ((flags & EIGS_VERBOSE) ||
	       MPIGLUE_CLOCK_DIFF(MPIGLUE_CLOCK, prev_feedback_time)
	       > FEEDBACK_TIME)) {
		if (mpb_verbosity >= 1) {
			printf("    iteration %4d: "
				"trace = %0.16g (%g%% change)\n", iteration, E,
				200.0 * fabs(E - prev_E) / (fabs(E) + fabs(prev_E)));
			fflush(stdout); /* make sure output appears */
	       }
               prev_feedback_time = MPIGLUE_CLOCK; /* reset feedback clock */
          }

	  if (iteration > 0 &&
              fabs(E - prev_E) < tolerance * 0.5 * (fabs(E) +
						    fabs(prev_E) + 1e-7))
               break; /* convergence!  hooray! */

	  /* compute new directions from residual & update basis: */
	  {
	       int ibasis2 = (ibasis + 1) % nbasis;

	       /* compute V[ibasis2] = AY */
#if 1
	       for (i = 0; i <= ibasis; ++i) {
		    evectmatrix_aXpbYS_sub(i ? 1.0 : 0.0, V[ibasis2],
					   1.0, AV[i],
					   S, itarget * q + Y.p * i, 1);
	       }
#else
	       A(Y, V[ibasis2], Adata, 1, Y);
#endif

	       /* handle restart case: */
	       if (ibasis2 == 0) {
		    evectmatrix_copy(AV[0], V[0]);
		    evectmatrix_copy(V[0], Y);
		    sqmatrix_resize(&VAV, Y.p, 0);
		    evectmatrix_XtY(VAV, V[0], AV[0], S3);
		    ibasis2 = 1;
		    evectmatrix_copy(V[ibasis2], AV[0]);
	       }

	       /* V[ibasis2] = residual = AY - Y * eigenvals */
	       matrix_XpaY_diag_real(V[ibasis2].data,
				     -1.0, Y.data,
				     eigenvals, Y.n, Y.p);

	       /* AV[ibasis2] = precondition V[ibasis2]: */
	       if (K != NULL)
		    K(V[ibasis2], AV[ibasis2], Kdata, Y, eigenvals, I);
	       else
		    evectmatrix_copy(AV[ibasis2], V[ibasis2]);

	       /* project by the constraints, if any: */
	       if (constraint)
		    constraint(AV[ibasis2], constraint_data);

	       /* orthogonalize against previous V: */
	       for (i = 0; i < ibasis2; ++i) {
		    evectmatrix_XtY(U, V[i], AV[ibasis2], S3);
		    evectmatrix_XpaYS(AV[ibasis2], -1.0, V[i], U, 0);
	       }

	       /* orthonormalize within itself: */
	       evectmatrix_XtX(U, AV[ibasis2], S3);
	       CHECK(sqmatrix_invert(U, 1, S3), "non-independent AV subspace");
	       sqmatrix_sqrt(S2, U, S3);
	       evectmatrix_XeYS(V[ibasis2], AV[ibasis2], S2, 1);

	       ibasis = ibasis2;
	  }

	  prev_E = E;
     } while (++iteration < EIGENSOLVER_MAX_ITERATIONS);

     CHECK(iteration < EIGENSOLVER_MAX_ITERATIONS,
           "failure to converge after "
           STRINGIZE(EIGENSOLVER_MAX_ITERATIONS)
           " iterations");

     evectmatrix_XtX(U, Y, S3);
     CHECK(sqmatrix_invert(U, 1, S3), "singular YtY at end");
     eigensolver_get_eigenvals_aux(Y, eigenvals, A, Adata,
				   V[0], AV[0], U, S3, S2);

     free(eigenvals2);

     destroy_sqmatrix(VAV);
     destroy_sqmatrix(S);
     destroy_sqmatrix(Swork);
     destroy_sqmatrix(U);
     destroy_sqmatrix(S2);
     destroy_sqmatrix(S3);
     destroy_sqmatrix(I);

     *num_iterations = iteration;
}
