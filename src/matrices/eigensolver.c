#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <config.h>
#include <check.h>
#include <scalar.h>
#include <matrices.h>

#include "eigensolver.h"

#define STRINGIZEx(x) #x /* a hack so that we can stringize macro values */
#define STRINGIZE(x) STRINGIZEx(x)
#define EIGENSOLVER_MAX_ITERATIONS 10000

/* Preconditioned eigensolver.  Finds the lowest Y.p eigenvectors
   and eigenvalues of the operator A, and returns them in Y and
   eigenvals (which should be an array of length Y.p).  C is a
   preconditioner, and should be an approximate inverse for A.
   Work[] contains workspace: nWork matrices.

   Initially, Y should contain a guess for the eigenvectors.

   nWork must be >= 2.  If nWork >= 3, preconditioned conjugate-
   gradient minimization is used; otherwise, we use preconditioned
   steepest-descent.  Currently, there is no advantage to using
   nWork > 3.

   tolerance is the convergence parameter.  Upon exit, 
   num_iterations holds the number of iterations that were required.

   NOTE: A and C are assumed to be linear operators. */

void eigensolver(evectmatrix Y, real *eigenvals,
		 evectoperator A, void *Adata,
		 evectoperator C, void *Cdata,
		 evectmatrix Work[], int nWork,
		 real tolerance, int *num_iterations)
{
     evectmatrix G, D, X;
     sqmatrix U, YtAYU, UYtAYU;
     short usingConjugateGradient = 0;
     real E, E2, prev_E = 0.0;
     real dE = 0.0, d2E, traceGtX, prev_traceGtX = 0.0;
     real lambda, prev_lambda = 0.001;
     int iteration = 0;
     
     CHECK(nWork >= 2, "not enough workspace");
     G = Work[0];
     X = Work[1];
     
     usingConjugateGradient = nWork >= 3;
     if (usingConjugateGradient)
	  D = Work[2];
     else
	  D = X;
     
     U = create_sqmatrix(Y.p);
     YtAYU = create_sqmatrix(Y.p);
     UYtAYU = create_sqmatrix(Y.p);
     
     /* notation: for a matrix Z, Zt denotes adjoint(Z) */
     
     /* The following loop performs an unconstrained minimization of
	the functional:

	E(Y) = trace [ Yt*A*Y / (Yt*Y) ]
	
	At the end, Y / sqrt(Yt*Y) will be the lowest eigenvectors.
	
	This is equivalent to minimizing trace[Zt*A*Z] under the
	constraint Zt*Z == 1.  (Z = Y / sqrt(Yt*Y)).  */
     
     do {
	  evectmatrix_XtX(U, Y);
	  sqmatrix_invert(U);
	  A(Y, X, Adata); /* X = AY */
	  evectmatrix_XeYS(G, X, U, 1); /* note that U = adjoint(U) */

	  evectmatrix_XtY(YtAYU, Y, G);
	  E = SCALAR_RE(sqmatrix_trace(YtAYU));

	  if (fabs(E - prev_E) < tolerance * 0.5 * (fabs(E) + fabs(prev_E)))
	       break; /* convergence!  hooray! */
	  
	  sqmatrix_AeBC(UYtAYU, U, 0, YtAYU, 0);
	  evectmatrix_XpaYS(G, -1.0, Y, UYtAYU); /* G is now the gradient of
						    the functional */
	  
	  C(G, X, Cdata);  /* X = precondition(G) */

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
	  A(Y, G, Adata); /* G = AY */
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
     } while (++iteration < EIGENSOLVER_MAX_ITERATIONS);
     
     CHECK(iteration < EIGENSOLVER_MAX_ITERATIONS,
	   "failure to converge after "
	   STRINGIZE(EIGENSOLVER_MAX_ITERATIONS)
	   " iterations");

     /* Now that we've converged, we need to find the actual eigenvectors
	and eigenvalues. */

     {
       sqmatrix Usqrt = create_sqmatrix(U.p);
       int i;

       sqmatrix_sqrt(Usqrt, U, UYtAYU); /* Usqrt = 1/sqrt(Yt*Y) */
       evectmatrix_XeYS(X, Y, Usqrt, 1);

#if 1
       evectmatrix_XtY(U, X, X);
       A(X, G, Adata);
       evectmatrix_XtY(U, X, G);
#else
       sqmatrix_AeBC(UYtAYU, Usqrt, 0, YtAYU, 0);
       sqmatrix_invert(Usqrt);
       sqmatrix_AeBC(U, UYtAYU, 0, Usqrt, 1); /* U == 1/sqrt(Yt*Y) * Yt *
					              A * Y * 1/sqrt(Yt*Y) */
#endif

       sqmatrix_eigensolve(U, eigenvals, YtAYU);
       evectmatrix_XeYS(Y, X, U, 1);

       destroy_sqmatrix(Usqrt);
     }
     
     *num_iterations = iteration;

     destroy_sqmatrix(U);
     destroy_sqmatrix(YtAYU);
     destroy_sqmatrix(UYtAYU);
}
