/* Copyright (C) 1999, 2000 Massachusetts Institute of Technology.
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

#include "../config.h"
#include <mpiglue.h>
#include <check.h>
#include <scalar.h>
#include <matrices.h>

#include "eigensolver.h"

static void eigensolver_get_eigenvals_aux(evectmatrix Y, real *eigenvals,
					  evectoperator A, void *Adata,
					  evectmatrix Work1, evectmatrix Work2,
					  sqmatrix U, sqmatrix Usqrt,
					  sqmatrix Uwork);

/* Check if we have converged yet, by seeing if fractional change
   in eigenvals (or their trace E) does not exceed the tolerance. */
static int check_converged(real E, real *eigenvals,
			   real prev_E, real *prev_eigenvals,
			   int num_bands,
			   real tolerance)
{
     if (!eigenvals || !prev_eigenvals)
	return(fabs(E - prev_E) < tolerance * 0.5 * (fabs(E) + fabs(prev_E)));
     else {
	  int i;

	  for (i = 0; i < num_bands; ++i) {
	  if (fabs(eigenvals[i] - prev_eigenvals[i]) >
	      tolerance * 0.5 * (fabs(eigenvals[i]) + fabs(prev_eigenvals[i])))
	       return 0;
	  }
	  
	  return 1;
     }
}

#define STRINGIZEx(x) #x /* a hack so that we can stringize macro values */
#define STRINGIZE(x) STRINGIZEx(x)
#define EIGENSOLVER_MAX_ITERATIONS 10000
#define FEEDBACK_TIME 4.0 /* elapsed time before we print progress feedback */

/* Number of iterations after which to reset conjugate gradient
   direction to steepest descent.  (Picked after some experimentation.
   Is there a better basis?  Should this change with the problem
   size?) */
#define CG_RESET_ITERS 100

/* maximum fractional change in trace for conjugate-gradient to be used;
   we don't want to use cg until we're near the minimum, so that the
   functional approximates a quadratic form. */
#define DELAY_CG_THRESH 0.1

/* Preconditioned eigensolver.  Finds the lowest Y.p eigenvectors
   and eigenvalues of the operator A, and returns them in Y and
   eigenvals (which should be an array of length Y.p).  
   Work[] contains workspace: nWork matrices.

   C is a preconditioner, and should be an approximate inverse for
   A - eigenvals[band].  If it is NULL then no preconditioning is done.

   Cdata_update is called to update Cdata when Y undergoes a
   rotation (multiplied by a unitary transformation).  (It
   can be NULL if you don't need to do anything to Cdata when
   there is a rotation.)

   constraint is a function which applies some other constraint to
   Y, or NULL if none.  Note that this can easily screw up
   convergence, so be careful.

   See eigensolver.h for more information on the form of 
   A, C, and Cdata_update.

   Initially, Y should contain a guess for the eigenvectors (or
   random data if you don't have a reasonable guess).

   nWork must be >= 2.  If nWork >= 3, preconditioned conjugate-
   gradient minimization is used; otherwise, we use preconditioned
   steepest-descent.  Currently, there is no advantage to using
   nWork > 3.

   tolerance is the convergence parameter.  Upon exit, 
   num_iterations holds the number of iterations that were required.

   flags indicate eigensolver options: different ways to do the
   iterative solver.  It should normally be EIGS_DEFAULT_FLAGS.

   NOTE: A and C are assumed to be linear operators. */

void eigensolver(evectmatrix Y, real *eigenvals,
		 evectoperator A, void *Adata,
		 evectpreconditioner C, void *Cdata,
 		 evectpreconditioner_data_updater Cdata_update,
		 evectconstraint constraint, void *constraint_data,
		 evectmatrix Work[], int nWork,
		 real tolerance, int *num_iterations,
		 int flags)
{
     evectmatrix G, D, X;
     sqmatrix U, YtAYU, S1, S2;
     short usingConjugateGradient = 0;
     real E, prev_E = 0.0;
     real y_norm;
     real *prev_eigenvals = NULL;
     real dE = 0.0, d2E, traceGtX, prev_traceGtX = 0.0;
     real lambda = 0.0, prev_lambda = 0.001, prev_change = -1e20;
     int i, iteration = 0;
     mpiglue_clock_t prev_feedback_time;
     
     prev_feedback_time = MPIGLUE_CLOCK;

     CHECK(nWork >= 2, "not enough workspace");
     G = Work[0];
     X = Work[1];
     
#ifdef DEBUG
     flags |= EIGS_VERBOSE;
#endif

     /* analytic-multimin only works when we are diagonalizing each step: */
     if ((flags & EIGS_ANALYTIC_MULTIMIN) &&
	 !(flags & EIGS_DIAGONALIZE_EACH_STEP) &&
	  (flags & EIGS_VERBOSE))
	  fprintf(stderr, "eigensolver: "
		  "must use diagonalize-each-step with analytic-multimin\n");
	  
     /* Only use CG when we have enough workspace.  Also, don't use it
	with EIGS_ANALYTIC_MULTIMIN since that screws up the CG condition */
     usingConjugateGradient = nWork >= 3 &&
	  !((flags & EIGS_ANALYTIC_MULTIMIN) &&
	    (flags & EIGS_DIAGONALIZE_EACH_STEP));
     if (usingConjugateGradient) {
	  int i;
	  D = Work[2];
	  /* we must initialize D to zero even though we multiply
	     it by zero (initial gamma) later on...otherwise, D
	     might contain NaN values */
	  for (i = 0; i < D.n * D.p; ++i)
	       ASSIGN_ZERO(D.data[i]);
     }
     else
	  D = X;

     /* U = 1/(Yt Y), YtAYU = Yt A Y U in the loop below,
	except when we are doing EIGS_DIAGONALIZE_EACH_STEP, in which
        case U = identity and YtAYU = diag(eigenvals), so we don't
        store them. */
     U = create_sqmatrix(Y.p);
     YtAYU = create_sqmatrix(Y.p);

     /* scratch matrices */
     S1 = create_sqmatrix(Y.p);
     S2 = create_sqmatrix(Y.p);

     if (flags & EIGS_DIAGONALIZE_EACH_STEP &&
	 flags & EIGS_CONVERGE_EACH_EIGENVALUE) {
	  CHK_MALLOC(prev_eigenvals, real, Y.p);
	  for (i = 0; i < Y.p; ++i)
	       prev_eigenvals[i] = 0.0;
     }

     for (i = 0; i < Y.p; ++i)
	  eigenvals[i] = 0.0;

     /* notation: for a matrix Z, Zt denotes adjoint(Z) */

     if (flags & EIGS_NORMALIZE_FIRST_STEP &&
         !(flags & EIGS_DIAGONALIZE_EACH_STEP)) {
	  if (constraint)
	       constraint(Y, constraint_data);
	  evectmatrix_XtX(U, Y);
	  sqmatrix_invert(U);
	  sqmatrix_sqrt(S1, U, S2); /* S1 = 1/sqrt(Yt*Y) */
	  evectmatrix_XeYS(X, Y, S1, 1);
	  evectmatrix_copy(Y, X);
     }
     
     /* The following loop performs an unconstrained minimization of
	the functional:

	E(Y) = trace [ Yt*A*Y / (Yt*Y) ]
	
	At the end, Y / sqrt(Yt*Y) will be the lowest eigenvectors.
	
	This is equivalent to minimizing trace[Zt*A*Z] under the
	constraint Zt*Z == 1.  (Z = Y / sqrt(Yt*Y)).  */
     
     do {
	  if (constraint)
	       constraint(Y, constraint_data);

	  evectmatrix_XtX(U, Y);
	  y_norm = sqrt(Y.p / SCALAR_RE(sqmatrix_trace(U)));
	  sqmatrix_invert(U);

	  if (flags & EIGS_DIAGONALIZE_EACH_STEP) {
	       /* First, orthonormalize: */
	       sqmatrix_sqrt(S1, U, S2); /* S1 = 1/sqrt(Yt*Y) */
	       evectmatrix_XeYS(X, Y, S1, 1);

	       /* Note that we use YtAYU as a scratch matrix below;
		  we also don't compute the actual value of Yt A Y U at
		  the end of the block because it is just diag(eigenvals) */
	       
	       /* Now, compute eigenvectors: */
	       A(X, Y, Adata, 1); /* Y = AX = A (Y/sqrt(U)) */
	       evectmatrix_XtY(S2, X, Y);
	       sqmatrix_eigensolve(S2, eigenvals, YtAYU);
	       
	       /* Compute G, Y, and D in new basis: */
	       evectmatrix_XeYS(G, Y, S2, 1);
	       evectmatrix_XeYS(Y, X, S2, 1);
	       if (usingConjugateGradient) {
		    /* mult. by S1 (not unitary) breaks conjugacy (?) */
		    sqmatrix_AeBC(YtAYU, S1, 0, S2, 1);
		    evectmatrix_copy(X, D);
		    evectmatrix_XeYS(D, X, YtAYU, 0);
	       }
	       y_norm = 1.0; /* Yt Y is now the identity matrix */
	       
	       if (Cdata_update)
		    Cdata_update(Cdata, S2);
	       
	       /* skip this since YtG should be diag(eigenvals): */
	       /* evectmatrix_XtY(YtAYU, Y, G); */
	       /* (also don't bother to compute U = 1/(YtY) = identity) */
	  }
	  else {
	       A(Y, X, Adata, 1); /* X = AY */
	       evectmatrix_XeYS(G, X, U, 1); /* note that U = adjoint(U) */
	       evectmatrix_XtY(YtAYU, Y, G);
	  }

	  if (flags & EIGS_DIAGONALIZE_EACH_STEP)
	       E = matrix_diag_real_trace(eigenvals, Y.p);
	  else
	       E = SCALAR_RE(sqmatrix_trace(YtAYU));

	  if (flags & EIGS_DIAGONALIZE_EACH_STEP) {
	       real E_check = 0.0;
	       for (i = 0; i < Y.p; ++i)
		    E_check += eigenvals[i];
	       CHECK(fabs(E_check - E) < 
		     (fabs(E) + fabs(E_check)) * 0.5 * tolerance,
		     "eigenvalue sum does not match trace!");
	  }

	  CHECK(!BADNUM(E), "crazy number detected in trace!!\n");

	  if (prev_change > E - prev_E) {
	       if (flags & EIGS_VERBOSE)
		    printf("    trace decreased more on it. %d than before\n",
			   iteration);
	  }
	  if (E > prev_E && iteration > 0) {
	       if (flags & EIGS_VERBOSE)
		    printf("    trace increased on it. %d!\n", iteration);
	  }

	  if (iteration > 0 && check_converged(E, eigenvals, 
					       prev_E, prev_eigenvals, Y.p,
					       tolerance))
	       break; /* convergence!  hooray! */

#if 0  /* yikes, don't do this: we need the eigenvals array below! */
	  if (flags & EIGS_DIAGONALIZE_EACH_STEP) {
	       /* initial eigenvalues are completely bogus and shouldn't be
		  used in preconditioner; kill them: */
	       if (iteration == 0 || 
		   fabs(E - prev_E) > 0.1 * 0.5 * (fabs(E) + fabs(prev_E)))
		    for (i = 0; i < Y.p; ++i)
			 eigenvals[i] = 0.0;
	  }
#endif

	  /* disable conjugate-gradient until we are close to the minimum */
	  if ((flags & EIGS_DELAY_CG) &&
	      (iteration == 0 ||
	       fabs(E - prev_E) > DELAY_CG_THRESH * 
	       0.5 * (fabs(E) + fabs(prev_E)))) {
	       prev_traceGtX = 0.0;	       
	  }

	  if ((flags & EIGS_VERBOSE) ||
	      MPIGLUE_CLOCK_DIFF(MPIGLUE_CLOCK, prev_feedback_time)
	      > FEEDBACK_TIME) {
	       printf("    iteration %4d: "
		      "trace = %g (%g%% change)\n", iteration + 1, E,
		      200.0 * fabs(E - prev_E) / (fabs(E) + fabs(prev_E)));
	       fflush(stdout); /* make sure output appears */
	       prev_feedback_time = MPIGLUE_CLOCK; /* reset feedback clock */
	  }

	  /* Compute gradient of functional G = (1 - Y U Yt) A Y U: */
	  if (flags & EIGS_DIAGONALIZE_EACH_STEP) {
	       /* G = G - Y * (YtAYU * U), optimized since
		  U is the identity and YtAYU = diag(eigenvals) */
	       matrix_XpaY_diag_real(G.data, -1.0, Y.data, eigenvals, G.n,G.p);
	  }
	  else {
	       sqmatrix_AeBC(S2, U, 0, YtAYU, 0);
	       evectmatrix_XpaYS(G, -1.0, Y, S2);
	  }
	  
	  if (C != NULL)
	       C(G, X, Cdata, Y, eigenvals);  /* X = precondition(G) */
	  else
	       evectmatrix_copy(X, G);  /* X = G if no preconditioner */

	  if (flags & EIGS_PROJECT_PRECONDITIONING) {
	       /* Operate projection P = (1 - Y U Yt) on X: */
	       evectmatrix_XtY(S1, Y, X);  /* S1 = Yt X */

	       if (flags & EIGS_DIAGONALIZE_EACH_STEP) {
		    /* U is the identity matrix */
		    evectmatrix_XpaYS(X, -1.0, Y, S1);
	       }
	       else {
		    sqmatrix_AeBC(S2, U, 0, S1, 0);
		    evectmatrix_XpaYS(X, -1.0, Y, S2);
	       }
	  }

	  traceGtX = 2.0 * SCALAR_RE(evectmatrix_traceXtY(G,X));
	  
	  /* In conjugate-gradient, the minimization direction D is
	     a combination of X with the previous search directions.
	     Otherwise, we just use D = X.
	     
	     We must also compute the derivative dE of E along the
	     search direction.  This is given by 2 Re[trace(Gt*D)]. */
	  
	  if (usingConjugateGradient) {
	       real gamma;
	       
	       if (prev_traceGtX == 0.0)
		    gamma = 0.0;
	       else
		    gamma = traceGtX / prev_traceGtX;

	       if ((flags & EIGS_RESET_CG) &&
		   (iteration + 1) % CG_RESET_ITERS == 0) {
		    /* periodically forget previous search directions,
		       and just juse D = X */
		    gamma = 0.0;
		    if (flags & EIGS_VERBOSE)
			 printf("    resetting CG direction...\n");
	       }
	       
	       evectmatrix_aXpbY(gamma, D, 1.0, X);
	       
	       if (gamma != 0.0)
		    dE = 2.0 * SCALAR_RE(evectmatrix_traceXtY(G,D));
	       else
		    dE = traceGtX;
	  }
	  else {
	       dE = traceGtX;
	       D = X;
	  }

	  /*** Minimization of trace along D direction ***/

	  if ((flags & EIGS_ANALYTIC_MULTIMIN) && 
	      (flags & EIGS_DIAGONALIZE_EACH_STEP)) {
	       real *YtAX_diag, *XtAX_diag, *XtX_diag;
	       real x_norm;

	       /* For this minimization method, we write
		      Y' = cos(lambda) Y + sin(lambda) X,
		  where X is the projection of D perpendicular to Y
		  (Yt X = 0).  (We also normalize the traces of Y and
		  X to equal p.)  Here, lambda is a pxp matrix,
		  diag(lambda_diag).  Thus, we are minimizing over a
		  p-dimensional space.  We then approximate our
		  function E(lambda) by dropping terms of O(lambda^3)
		  and higher, and thus can minimize analytically.
		  This is actually exact for the case of p=1 (1
		  eigenvector).

		  This stuff only works if we are diagonalizing on
	          each step!! */
	       
	       /* make a copy of D, if we are using conj. gradient.
		  (otherwise, X already equals D).  This is necessary
		  to avoid screwing up the conjugacy condition when
		  we project to space orthogonal to Y, below */
	       if (usingConjugateGradient)
		    evectmatrix_copy(X, D);

	       /* make sure Yt X = 0 */
	       if (usingConjugateGradient ||
		   (!(flags & EIGS_PROJECT_PRECONDITIONING) && C != NULL)) {
		    /* Operate projection P = (1 - Y U Yt) on D: */
		    evectmatrix_XtY(S1, Y, X);  /* S1 = Yt X */
		    
		    /* U is the identity matrix */
		    evectmatrix_XpaYS(X, -1.0, Y, S1);
	       }

	       /* Now, we need some 1d arrays of length p to store the
		  diagonals of our matrices.  We'll just use the
		  storage of S1, S2, and YtAYU for this (they hold pxp
		  scalar values), since we aren't using them for
		  anything else. */
	       YtAX_diag = (real *) S1.data;
	       XtAX_diag = (real *) S2.data;
	       XtX_diag = (real *) YtAYU.data;

	       A(X, G, Adata, 0); /* G = A X */

	       evectmatrix_XtY_diag_real(Y, G, YtAX_diag);
	       evectmatrix_XtY_diag_real(X, G, XtAX_diag);
	       evectmatrix_XtX_diag_real(X, XtX_diag);
	       x_norm = sqrt(X.p / matrix_diag_real_trace(XtX_diag, Y.p));

	       /* now compute the values of diag(lambda).  Actually,
		  compute cos(lambda) * y_norm and sin(lambda) * x_norm,
		  and store these values in XtAX_diag and XtX_diag,
	          for combining Y and X in the next step. */
	       for (i = 0; i < Y.p; ++i) {
		    lambda = 0.5 * atan2(-2.0 * YtAX_diag[i] * x_norm * y_norm,
				    (XtAX_diag[i] -
				     XtX_diag[i] * eigenvals[i]*y_norm*y_norm)
				    * x_norm*x_norm);
		    XtAX_diag[i] = cos(lambda) * y_norm;
		    XtX_diag[i] = sin(lambda) * x_norm;
	       }
	       
	       /* Y = Y * cos(lambda)*y_norm + X * sin(lambda)*x_norm */
               matrix_X_diag_real_pY_diag_real(Y.data, XtAX_diag,
					       X.data, XtX_diag, Y.n, Y.p);
	  }
	  else if (flags & EIGS_ANALYTIC_LINMIN) {
	       real traceYtAXU, traceUYtAYUXtX, traceXtAXU;
	       real x_norm;

	       /* For this minimization method, we write
		      Y' = cos(lambda) Y + sin(lambda) X,
		  where X is the projection of D perpendicular to Y
		  (Yt X = 0).  (We also normalize the traces of Y and
		  X to equal p.) We then approximate our function
		  E(lambda) by dropping terms of O(lambda^3) and
		  higher, and thus can minimize analytically.  This is
		  actually exact for the case of p=1 (1 eigenvector). */
	       
	       /* make a copy of D, if we are using conj. gradient.
		  (otherwise, X already equals D).  This is necessary
		  to avoid screwing up the conjugacy condition when
		  we project to space orthogonal to Y, below */
	       if (usingConjugateGradient)
		    evectmatrix_copy(X, D);

	       /* make sure Yt X = 0 */
	       if (usingConjugateGradient ||
		   (!(flags & EIGS_PROJECT_PRECONDITIONING) && C != NULL)) {
		    /* Operate projection P = (1 - Y U Yt) on D: */
		    evectmatrix_XtY(S1, Y, X);  /* S1 = Yt X */
		    
		    if (flags & EIGS_DIAGONALIZE_EACH_STEP) {
			 /* U is the identity matrix */
			 evectmatrix_XpaYS(X, -1.0, Y, S1);
		    }
		    else {
			 sqmatrix_AeBC(S2, U, 0, S1, 0);
			 evectmatrix_XpaYS(X, -1.0, Y, S2);
		    }
	       }

	       A(X, G, Adata, 0); /* G = A X */
	       
	       if (flags & EIGS_DIAGONALIZE_EACH_STEP) {
		    /* U is the identity matrix */
		    traceYtAXU = SCALAR_RE(evectmatrix_traceXtY(Y, G));
		    traceXtAXU = SCALAR_RE(evectmatrix_traceXtY(X, G));

		    evectmatrix_XtY(S1, X, X); /* S1 = Xt X */
		    traceUYtAYUXtX = 
			 matrix_re_trace_A_diag_real(S1.data, eigenvals, S1.p);

		    x_norm = sqrt(X.p / SCALAR_RE(sqmatrix_trace(S1)));
	       }
	       else {
		    evectmatrix_XtY(S1, Y, G); /* S1 = Yt A X */
		    traceYtAXU = SCALAR_RE(sqmatrix_traceAtB(U, S1));

		    evectmatrix_XtY(S1, X, G); /* S1 = Xt A X */
		    traceXtAXU = SCALAR_RE(sqmatrix_traceAtB(U, S1));

		    evectmatrix_XtY(S1, X, X); /* S1 = Xt X */
		    /* S2 = YtAYU Xt X: */
		    sqmatrix_AeBC(S2, YtAYU, 0, S1, 1);
		    traceUYtAYUXtX = SCALAR_RE(sqmatrix_traceAtB(U, S2));

		    x_norm = sqrt(X.p / SCALAR_RE(sqmatrix_trace(S1)));
	       }

	       lambda = 0.5 * atan2(-2.0 * traceYtAXU * x_norm * y_norm,
				    (traceXtAXU -
				     traceUYtAYUXtX * y_norm*y_norm)
				    * x_norm*x_norm);
	       if ((flags & EIGS_VERBOSE) && fabs(lambda) > 0.3)
		    printf("    large lambda = %g\n", lambda);
	       evectmatrix_aXpbY(cos(lambda)*y_norm, Y,
				 sin(lambda)*x_norm, X);
	  }
	  else {
	       real E2;

	       /* For this minimization method, we use the value
		  of the function and its derivative at the current point,
		  along with its value at a second point, to construct
		  a quadratic approximation for the function which
		  we can then minimize. */

	       /* Now, let's evaluate the functional at a point slightly
		  along the current direction, where "slightly" means
		  half the previous stepsize: */
	       evectmatrix_aXpbY(1.0, Y, prev_lambda*0.5, D);
	       evectmatrix_XtX(U, Y);
	       sqmatrix_invert(U);
	       A(Y, G, Adata, 0); /* G = AY */
	       evectmatrix_XtY(S1, Y, G);
	       E2 = SCALAR_RE(sqmatrix_traceAtB(U, S1));
	       
	       /* Minimizing Y + lambda * D: */
	       
	       /* At this point, we know the value of the function at Y (E),
		  the derivative (dE), and the value at a second point Y' (E2).
		  We fit this data to a quadratic and use that to predict the
		  minimum along the direction D: */
	       
	       d2E = 2.0 * (E2 - E - dE * prev_lambda*0.5) /
		    (prev_lambda*prev_lambda*0.25);
	       
	       /* Actually, we'll model things by a cos() curve, with
		  d2E being an approx. for the 2nd derivative, since
		  this is more well-behaved for small d2E: */
	       
	       lambda = 0.5 * atan2(-dE, 0.5*d2E);

	       if ((flags & EIGS_VERBOSE) && fabs(lambda) > 3.14159*0.2)
		    printf("    it. %d: large lambda/pi = %g\n", 
			   iteration + 1, lambda/3.14159);
	       
	       /* Compute new Y.  Note that we have already shifted Y. */
	       evectmatrix_aXpbY(1.0, Y, lambda - prev_lambda*0.5, D);
	  }

	  prev_traceGtX = traceGtX;
	  prev_lambda = lambda;
	  prev_change = E - prev_E;
	  prev_E = E;

	  if (prev_eigenvals)
	       for (i = 0; i < Y.p; ++i)
		    prev_eigenvals[i] = eigenvals[i];
     } while (++iteration < EIGENSOLVER_MAX_ITERATIONS);
     
     CHECK(iteration < EIGENSOLVER_MAX_ITERATIONS,
	   "failure to converge after "
	   STRINGIZE(EIGENSOLVER_MAX_ITERATIONS)
	   " iterations");

     if (!(flags & EIGS_DIAGONALIZE_EACH_STEP)) {
	  /* Now that we've converged, we need to find the actual eigenvectors
	     and eigenvalues. */

	  eigensolver_get_eigenvals_aux(Y, eigenvals, A, Adata,
					X, G, U, S1, S2);
     }
     
     *num_iterations = iteration;

     destroy_sqmatrix(U);
     destroy_sqmatrix(S1);
     destroy_sqmatrix(S2);
     destroy_sqmatrix(YtAYU);

     if (prev_eigenvals)
	  free(prev_eigenvals);
}


static void eigensolver_get_eigenvals_aux(evectmatrix Y, real *eigenvals,
					  evectoperator A, void *Adata,
					  evectmatrix Work1, evectmatrix Work2,
					  sqmatrix U, sqmatrix Usqrt,
					  sqmatrix Uwork)
{
     sqmatrix_sqrt(Usqrt, U, Uwork); /* Usqrt = 1/sqrt(Yt*Y) */
     evectmatrix_XeYS(Work1, Y, Usqrt, 1);

     A(Work1, Work2, Adata, 1);
     evectmatrix_XtY(U, Work1, Work2);

     sqmatrix_eigensolve(U, eigenvals, Uwork);
     evectmatrix_XeYS(Y, Work1, U, 1);
}

void eigensolver_get_eigenvals(evectmatrix Y, real *eigenvals,
			       evectoperator A, void *Adata,
			       evectmatrix Work1, evectmatrix Work2)
{
     sqmatrix U, Usqrt, Uwork;
 
     U = create_sqmatrix(Y.p);
     Usqrt = create_sqmatrix(Y.p);
     Uwork = create_sqmatrix(Y.p);

     evectmatrix_XtX(U, Y);
     sqmatrix_invert(U);

     eigensolver_get_eigenvals_aux(Y, eigenvals, A, Adata, Work1, Work2,
				   U, Usqrt, Uwork);

     destroy_sqmatrix(U);
     destroy_sqmatrix(Usqrt);
     destroy_sqmatrix(Uwork);
}

/* Subroutines for chaining constraints, to make it easy to pass
   multiple constraint functions to the eigensolver: */

evectconstraint_chain *evect_add_constraint(evectconstraint_chain *constraints,
					    evectconstraint C,
					    void *constraint_data)
{
     evectconstraint_chain *new_constraints;

     CHK_MALLOC(new_constraints, evectconstraint_chain, 1);

     new_constraints->C = C;
     new_constraints->constraint_data = constraint_data;
     new_constraints->next = constraints;
     return new_constraints;
}

void evect_destroy_constraints(evectconstraint_chain *constraints)
{
     while (constraints) {
	  evectconstraint_chain *cur_constraint = constraints;
	  constraints = constraints->next;
	  free(cur_constraint);
     }
}

void evectconstraint_chain_func(evectmatrix X, void *data)
{
     evectconstraint_chain *constraints = (evectconstraint_chain *) data;

     while (constraints) {
	  if (constraints->C)
	       constraints->C(X, constraints->constraint_data);
          constraints = constraints->next;
     }
}
