#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <config.h>
#include <check.h>
#include <scalar.h>
#include <matrices.h>

#include "eigensolver.h"

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

#define NORMALIZE_FIRST_STEP 0
#define DIAGONALIZE_EACH_STEP 1
#define PROJECT_PRECONDITIONING 0
#define CONVERGE_EACH_EIGENVALUE 1

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

   NOTE: A and C are assumed to be linear operators. */

void eigensolver(evectmatrix Y, real *eigenvals,
		 evectoperator A, void *Adata,
		 evectpreconditioner C, void *Cdata,
 		 evectpreconditioner_data_updater Cdata_update,
		 evectconstraint constraint, void *constraint_data,
		 evectmatrix Work[], int nWork,
		 real tolerance, int *num_iterations)
{
     evectmatrix G, D, X;
     sqmatrix U, Usqrt, YtAYU, UYtAYU;
     short usingConjugateGradient = 0;
     real E, E2, prev_E = 0.0;
     real *prev_eigenvals = NULL;
     real dE = 0.0, d2E, traceGtX, prev_traceGtX = 0.0;
     real lambda, prev_lambda = 0.001;
     int i, iteration = 0;
     
     CHECK(nWork >= 2, "not enough workspace");
     G = Work[0];
     X = Work[1];
     
     usingConjugateGradient = nWork >= 3;
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

     /* the following matrices are named after the original uses I put
	them to, but beware that I occasionally use them as scratch
	matrices for other purposes */
     U = create_sqmatrix(Y.p);
     Usqrt = create_sqmatrix(Y.p);
     YtAYU = create_sqmatrix(Y.p);
     UYtAYU = create_sqmatrix(Y.p);
     
#if DIAGONALIZE_EACH_STEP && CONVERGE_EACH_EIGENVALUE
     prev_eigenvals = (real*) malloc(sizeof(real) * Y.p);
     CHECK(prev_eigenvals, "out of memory");
     for (i = 0; i < Y.p; ++i)
	  prev_eigenvals[i] = 0.0;
#endif

     for (i = 0; i < Y.p; ++i)
	  eigenvals[i] = 0.0;

     /* notation: for a matrix Z, Zt denotes adjoint(Z) */

#if NORMALIZE_FIRST_STEP && ! DIAGONALIZE_EACH_STEP
     if (constraint)
	  constraint(Y, constraint_data);
     evectmatrix_XtX(U, Y);
     sqmatrix_invert(U);
     sqmatrix_sqrt(Usqrt, U, UYtAYU); /* Usqrt = 1/sqrt(Yt*Y) */
     evectmatrix_XeYS(X, Y, Usqrt, 1);
     evectmatrix_copy(Y, X);
#endif
     
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
	  sqmatrix_invert(U);

#if DIAGONALIZE_EACH_STEP
	  
	  /* First, orthonormalize: */
	  sqmatrix_sqrt(Usqrt, U, UYtAYU); /* Usqrt = 1/sqrt(Yt*Y) */
	  evectmatrix_XeYS(X, Y, Usqrt, 1);

	  /* Now, compute eigenvectors: */
	  A(X, Y, Adata, 1); /* Y = AX = A (Y/sqrt(U)) */
	  evectmatrix_XtY(YtAYU, X, Y);
	  sqmatrix_eigensolve(YtAYU, eigenvals, UYtAYU);

	  /* Compute G, Y, and D in new basis: */
	  evectmatrix_XeYS(G, Y, YtAYU, 1);
	  evectmatrix_XeYS(Y, X, YtAYU, 1);
	  if (usingConjugateGradient) {
	       /* mult. by Usqrt (not unitary) breaks conjugacy (?) */
	       sqmatrix_AeBC(UYtAYU, Usqrt, 0, YtAYU, 1);
	       evectmatrix_copy(X, D);
	       evectmatrix_XeYS(D, X, UYtAYU, 0);
	  }  

	  if (Cdata_update)
	       Cdata_update(Cdata, YtAYU);
     
	  /* skip this since YtG should be diag(eigenvals)? */
	  evectmatrix_XtY(YtAYU, Y, G);

#else

	  A(Y, X, Adata, 1); /* X = AY */
	  evectmatrix_XeYS(G, X, U, 1); /* note that U = adjoint(U) */
	  evectmatrix_XtY(YtAYU, Y, G);

#endif

	  E = SCALAR_RE(sqmatrix_trace(YtAYU));

#if DIAGONALIZE_EACH_STEP
	  E2 = 0.0;
	  for (i = 0; i < Y.p; ++i)
	       E2 += eigenvals[i];
	  CHECK(fabs(E2 - E) < (fabs(E) + fabs(E2)) * 0.5 * tolerance,
		"eigenvalue sum does not match trace!");
#endif

	  CHECK(!BADNUM(E), "crazy number detected in trace!!\n");

	  if (iteration > 0 && check_converged(E, eigenvals, 
					       prev_E, prev_eigenvals, Y.p,
					       tolerance))
	       break; /* convergence!  hooray! */

#if DIAGONALIZE_EACH_STEP
	  /* initial eigenvalues are completely bogus and shouldn't be
	     used in preconditioner; kill them: */
	  if (iteration == 0 || 
	      fabs(E - prev_E) > 0.1 * 0.5 * (fabs(E) + fabs(prev_E)))
	       for (i = 0; i < Y.p; ++i)
		    eigenvals[i] = 0.0;
#endif
	  
#ifdef DEBUG
	  if (iteration % 10 == 0)
	       printf("it. %4d: tr. = %g (%g%% change)\n", iteration, E,
		      200.0 * fabs(E - prev_E) / (fabs(E) + fabs(prev_E)));
#endif

#if DIAGONALIZE_EACH_STEP /* U is the identity matrix */
	  /* optimize this since YtAYU should be diagonal? */
	  evectmatrix_XpaYS(G, -1.0, Y, YtAYU); /* G is now the gradient of
						    the functional */
#else
	  sqmatrix_AeBC(UYtAYU, U, 0, YtAYU, 0);
	  evectmatrix_XpaYS(G, -1.0, Y, UYtAYU); /* G is now the gradient of
						    the functional */
#endif
	  
	  if (C)
	       C(G, X, Cdata, Y, eigenvals);  /* X = precondition(G) */
	  else
	       evectmatrix_copy(X, G);  /* X = G if no preconditioner */

#if PROJECT_PRECONDITIONING  /* Operate projection P = (1 - Y U Yt) on X: */
	  evectmatrix_XtY(Usqrt, Y, X);  /* Usqrt = Yt X */
	  
#  if DIAGONALIZE_EACH_STEP /* U is the identity matrix */
          evectmatrix_XpaYS(X, -1.0, Y, Usqrt);
#  else
          sqmatrix_AeBC(UYtAYU, U, 0, Usqrt, 0);
          evectmatrix_XpaYS(X, -1.0, Y, UYtAYU);
#  endif
#endif /* PROJECT_PRECONDITIONING */

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
	  
	  /* Now, let's evaluate the functional at a point slightly
	     along the current direction, where "slightly" means
	     half the previous stepsize: */
	  evectmatrix_aXpbY(1.0, Y, prev_lambda*0.5, D);
	  evectmatrix_XtX(U, Y);
	  sqmatrix_invert(U);
	  A(Y, G, Adata, 0); /* G = AY */
	  evectmatrix_XtY(YtAYU, Y, G);
	  E2 = SCALAR_RE(sqmatrix_traceAtB(U, YtAYU));
	  
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
	  
	  /* Compute new Y.  Note that we have already shifted Y slightly. */
	  evectmatrix_aXpbY(1.0, Y, lambda - prev_lambda*0.5, D);
	  
	  prev_traceGtX = traceGtX;
	  prev_lambda = lambda;
	  prev_E = E;

	  if (prev_eigenvals)
	       for (i = 0; i < Y.p; ++i)
		    prev_eigenvals[i] = eigenvals[i];
     } while (++iteration < EIGENSOLVER_MAX_ITERATIONS);
     
     CHECK(iteration < EIGENSOLVER_MAX_ITERATIONS,
	   "failure to converge after "
	   STRINGIZE(EIGENSOLVER_MAX_ITERATIONS)
	   " iterations");

#if !DIAGONALIZE_EACH_STEP
     /* Now that we've converged, we need to find the actual eigenvectors
	and eigenvalues. */

     sqmatrix_sqrt(Usqrt, U, UYtAYU); /* Usqrt = 1/sqrt(Yt*Y) */
     evectmatrix_XeYS(X, Y, Usqrt, 1);
     
     A(X, G, Adata, 1);
     evectmatrix_XtY(U, X, G);
     
     sqmatrix_eigensolve(U, eigenvals, YtAYU);
     evectmatrix_XeYS(Y, X, U, 1);
#endif
     
     *num_iterations = iteration;

     destroy_sqmatrix(U);
     destroy_sqmatrix(Usqrt);
     destroy_sqmatrix(YtAYU);
     destroy_sqmatrix(UYtAYU);

     if (prev_eigenvals)
	  free(prev_eigenvals);
}

