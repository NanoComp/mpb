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
#include <blasglue.h>

#include "eigensolver.h"

extern void eigensolver_get_eigenvals_aux(evectmatrix Y, real *eigenvals,
                                          evectoperator A, void *Adata,
                                          evectmatrix Work1, evectmatrix Work2,
                                          sqmatrix U, sqmatrix Usqrt,
                                          sqmatrix Uwork);

#define STRINGIZEx(x) #x /* a hack so that we can stringize macro values */
#define STRINGIZE(x) STRINGIZEx(x)

#define K_PI 3.141592653589793238462643383279502884197
#define MIN2(a,b) ((a) < (b) ? (a) : (b))
#define MAX2(a,b) ((a) > (b) ? (a) : (b))

/**************************************************************************/

#define EIGENSOLVER_MAX_ITERATIONS 10000
#define FEEDBACK_TIME 4.0 /* elapsed time before we print progress feedback */

/* Number of iterations after which to reset conjugate gradient
   direction to steepest descent.  (Picked after some experimentation.
   Is there a better basis?  Should this change with the problem
   size?) */
#define CG_RESET_ITERS 70

/**************************************************************************/

typedef struct {
     real d_norm;
     sqmatrix YtAY, DtAD, symYtAD, YtY, DtD, symYtD, S1, S2, S3;
} trace_func_data;

static real trace_func(real *trace_deriv,
		       real theta, void *data)
{
     real trace;
     trace_func_data *d = (trace_func_data *) data;

     {
	  real c = cos(theta), s = sin(theta) / d->d_norm;
	  sqmatrix_copy(d->S1, d->YtY);
	  sqmatrix_aApbB(c*c, d->S1, s*s, d->DtD);
	  sqmatrix_ApaB(d->S1, 2*s*c, d->symYtD);
	  sqmatrix_invert(d->S1);
	  
	  sqmatrix_copy(d->S2, d->YtAY);
	  sqmatrix_aApbB(c*c, d->S2, s*s, d->DtAD);
	  sqmatrix_ApaB(d->S2, 2*s*c, d->symYtAD);
	  
	  trace = SCALAR_RE(sqmatrix_traceAtB(d->S2, d->S1));
     }

     if (trace_deriv) {
	  real c2 = cos(2*theta), s2 = sin(2*theta);
	  
	  sqmatrix_copy(d->S3, d->YtAY);
	  sqmatrix_ApaB(d->S3, -1.0/(d->d_norm*d->d_norm), d->DtAD);
	  sqmatrix_aApbB(-0.5 * s2, d->S3, c2/d->d_norm, d->symYtAD);
	  
	  *trace_deriv = SCALAR_RE(sqmatrix_traceAtB(d->S1, d->S3));
	  
	  sqmatrix_AeBC(d->S3, d->S1, 0, d->S2, 1);
	  sqmatrix_AeBC(d->S2, d->S3, 0, d->S1, 1);
	  
	  sqmatrix_copy(d->S3, d->YtY);
	  sqmatrix_ApaB(d->S3, -1.0/(d->d_norm*d->d_norm), d->DtD);
	  sqmatrix_aApbB(-0.5 * s2, d->S3, c2/d->d_norm, d->symYtD);
	  
	  *trace_deriv -= SCALAR_RE(sqmatrix_traceAtB(d->S2, d->S3));
	  *trace_deriv *= 2;
     }

     return trace;
}

/**************************************************************************/

static real linmin(real *improvement,
		   real xmin, real f_xmin, real df_xmin, real xmax, real x0,
		   real tolerance,
		   real (*func) (real *df, real x, void *data), void *data)
{
     real x_prev, x, dx, f, df, f_x0, df_x0, f_xmax, df_xmax, s, f_xstart;
     int eval_count = 0, is_xstart;

     CHECK(df_xmin * (x0 - xmin) < 0.0, "linmin: bad initial guess!");

     s = xmax > xmin ? 1.0 : -1.0;
     CHECK(x0*s < xmax*s && x0*s > xmin*s,
	   "linmin: initial guess out of range");

     /* first, bracket the minimum: 
        (the following works, but is not very smart...fix it!) */

     do {
	  real xmin2 = xmin, f_xmin2 = f_xmin, df_xmin2 = df_xmin;
	  dx = (x0 - xmin) * 2.0;
	  for (x = xmin + dx; x*s <= xmax*s; x += dx) {
	       f = func(&df, x, data);
	       ++eval_count;
	       if (df * (x - xmin) > 0.0)
		    break;
	       else {
		    xmin2 = x; f_xmin2 = f; df_xmin2 = df;
	       }
	  }
	  if (x*s <= xmax*s) {
	       xmin = xmin2; f_xmin = f_xmin2; df_xmin = df_xmin2;
	       xmax = x; f_xmax = f; df_xmax = df;
	       break;
	  }
	  x0 = 0.5 * (x0 + xmin);
     } while (fabs(x0 - xmin) > tolerance * (fabs(x0) + tolerance));
     CHECK(fabs(x0 - xmin) > tolerance * (fabs(x0) + tolerance),
	   "linmin: failed to bracket minimum!");

     if (x0*s <= xmin*s || x0*s >= xmax*s)
	  x0 = 0.5 * (xmin + xmax);

     /* Now, find the root of the derivative by Ridder's Method.

	Replace this with a more robust algorithm at some point so
	that we don't accidentally converge to a maximum?  (If there
	is one between xmin and xmax.)  */

     if (xmin > xmax) {
	  x = xmin; f = f_xmin; df = df_xmin;
	  xmin = xmax; f_xmin = f_xmax; df_xmin = df_xmax;
	  xmax = x; f_xmax = f; df_xmax = df;
     }
     
     is_xstart = 1;
     x_prev = x0;
     while (1) {
	  f_x0 = func(&df_x0, x0, data);
	  ++eval_count;
	  if (is_xstart) {
	       f_xstart = f_x0;
	       is_xstart = 0;
	  }

	  if (df_x0 == 0)
	       break;
	  if (df_xmin == 0) {
	       x0 = xmin;
	       break;
	  }
	  if (df_xmax == 0) {
	       x0 = xmax;
	       break;
	  }

	  x = x0 + (x0 - xmin) * (df_xmin > df_xmax ? 1 : -1) *
	       df_x0 / sqrt(df_x0*df_x0 - df_xmin*df_xmax);

	  if (MAX2(fabs(x - x_prev), MIN2(fabs(x - xmin), fabs(x - xmax))) <
	      tolerance * (fabs(x) + tolerance)) {
	       x0 = x;
	       break;
	  }

	  f = func(&df, x, data);
	  ++eval_count;

	  if (df * df_x0 > 0 || (df - df_x0) * (x - x0) < 0) {
	       if (x < x0) {
		    if (df_xmin*df > 0 || (df_xmin - df)*(xmin - x) < 0) {
			 xmin = x0; f_xmin = f_x0; df_xmin = df_x0;
		    }
		    else {
			 xmax = x; f_xmax = f; df_xmax = df;
		    }
	       }
	       else if (df_xmin*df_x0 > 0 || (df_xmin - df_x0)*(xmin-x0) < 0) {
		    xmin = x; f_xmin = f; df_xmin = df;
	       }
	       else {
		    xmax = x0; f_xmax = f_x0; df_xmax = df_x0;
	       }
	  }
	  else {
	       if (x < x0) {
		    xmin = x; f_xmin = f; df_xmin = df;
		    xmax = x0; f_xmax = f_x0; df_xmax = df_x0;
	       }
	       else {
		    xmin = x0; f_xmin = f_x0; df_xmin = df_x0;
		    xmax = x; f_xmax = f; df_xmax = df;
	       }
	  }

	  x0 = 0.5 * (xmin + xmax);
	  x_prev = x;
     }

     f_x0 = func(NULL, x0, data);
     ++eval_count;
     *improvement = (f_xstart - f_x0) * 2.0 /
	  (fabs(f_xstart) + fabs(f_x0) + tolerance);

#if 0
     printf("linmin: from %g to %g (%g%% improvement) in %d evaluations.\n",
	    f_xstart, f_x0, *improvement * 100, eval_count);
#endif

     return x0;
}

/**************************************************************************/

void eigensolver(evectmatrix Y, real *eigenvals,
                 evectoperator A, void *Adata,
                 evectpreconditioner K, void *Kdata,
                 evectconstraint constraint, void *constraint_data,
                 evectmatrix Work[], int nWork,
                 real tolerance, int *num_iterations,
                 int flags)
{
     evectmatrix G, D, X, prev_G;
     short usingConjugateGradient = 0, use_polak_ribiere = 0,
	   use_linmin = 1;
     real E, prev_E = 0.0;
     real traceGtX, prev_traceGtX = 0.0;
     real theta, prev_theta = 0.5;
     int i, iteration = 0;
     mpiglue_clock_t prev_feedback_time;
     sqmatrix YtAYU, DtAD, symYtAD, YtY, U, DtD, symYtD, S1, S2, S3;
     trace_func_data tfd;

     prev_feedback_time = MPIGLUE_CLOCK;
     
#ifdef DEBUG
     flags |= EIGS_VERBOSE;
#endif

     CHECK(nWork >= 2, "not enough workspace");
     G = Work[0];
     X = Work[1];

     usingConjugateGradient = nWork >= 3;
     if (usingConjugateGradient) {
          D = Work[2];
          for (i = 0; i < D.n * D.p; ++i)
               ASSIGN_ZERO(D.data[i]);
     }
     else
          D = X;

     use_polak_ribiere = nWork >= 4;
     if (use_polak_ribiere) {
          prev_G = Work[3];
          for (i = 0; i < Y.n * Y.p; ++i)
               ASSIGN_ZERO(prev_G.data[i]);
     }
     else
          prev_G = G;
     
     YtAYU = create_sqmatrix(Y.p);  /* holds Yt A Y */
     DtAD = create_sqmatrix(Y.p);  /* holds Dt A D */
     symYtAD = create_sqmatrix(Y.p);  /* holds (Yt A D + Dt A Y) / 2 */
     YtY = create_sqmatrix(Y.p);  /* holds Yt Y */
     U = create_sqmatrix(Y.p);  /* holds 1 / (Yt Y) */
     DtD = create_sqmatrix(Y.p);  /* holds Dt D */
     symYtD = create_sqmatrix(Y.p);  /* holds (Yt D + Dt Y) / 2 */

     /* Notation note: "t" represents a dagger superscript, so
	Yt represents adjoint(Y), or Y' in MATLAB syntax. */

     /* scratch matrices: */
     S1 = create_sqmatrix(Y.p);
     S2 = create_sqmatrix(Y.p);
     S3 = create_sqmatrix(Y.p);

     tfd.YtAY = S1; tfd.DtAD = DtAD; tfd.symYtAD = symYtAD;
     tfd.YtY = YtY; tfd.DtD = DtD; tfd.symYtD = symYtD;
     tfd.S1 = YtAYU; tfd.S2 = S2; tfd.S3 = S3;

     for (i = 0; i < Y.p; ++i)
          eigenvals[i] = 0.0;

     if (constraint)
	  constraint(Y, constraint_data);

     do {
	  real y_norm, d_norm;

	  evectmatrix_XtX(YtY, Y);

	  y_norm = sqrt(SCALAR_RE(sqmatrix_trace(YtY)) / Y.p);
	  blasglue_scal(Y.p * Y.n, 1/y_norm, Y.data, 1);
	  blasglue_scal(Y.p * Y.p, 1/(y_norm*y_norm), YtY.data, 1);

	  sqmatrix_copy(U, YtY);

	  sqmatrix_invert(U);

	  A(Y, X, Adata, 1, G); /* X = AY; G is scratch */
	  evectmatrix_XeYS(G, X, U, 1); /* G = AYU; note that U is Hermitian */
	  evectmatrix_XtY(YtAYU, Y, G);

	  E = SCALAR_RE(sqmatrix_trace(YtAYU));
	  CHECK(!BADNUM(E), "crazy number detected in trace!!\n");

	  if (iteration > 0 &&
              fabs(E - prev_E) < tolerance * 0.5 * (E + prev_E + 1e-7))
               break; /* convergence!  hooray! */
	  
	  if ((flags & EIGS_VERBOSE) ||
              MPIGLUE_CLOCK_DIFF(MPIGLUE_CLOCK, prev_feedback_time)
              > FEEDBACK_TIME) {
               printf("    iteration %4d: "
                      "trace = %g (%g%% change)\n", iteration + 1, E,
                      200.0 * fabs(E - prev_E) / (fabs(E) + fabs(prev_E)));
               fflush(stdout); /* make sure output appears */
               prev_feedback_time = MPIGLUE_CLOCK; /* reset feedback clock */
          }

	  /* Compute gradient of functional: G = (1 - Y U Yt) A Y U */
	  sqmatrix_AeBC(S1, U, 0, YtAYU, 0);
	  evectmatrix_XpaYS(G, -1.0, Y, S1);

	  /* set X = precondition(G): */
	  if (K != NULL)
	       K(G, X, Kdata, 
		 Y, NULL /* no eigenvals; not diagonalized */, 
		 YtY);  
	  else
               evectmatrix_copy(X, G);  /* preconditioner is identity */
	  
	  if (flags & EIGS_PROJECT_PRECONDITIONING) {
               /* Operate projection P = (1 - Y U Yt) on X: */
               evectmatrix_XtY(symYtD, Y, X);  /* symYtD = Yt X */
	       sqmatrix_AeBC(S1, U, 0, symYtD, 0);
	       evectmatrix_XpaYS(X, -1.0, Y, S1);
          }

	  /* In conjugate-gradient, the minimization direction D is
             a combination of X with the previous search directions.
             Otherwise, we just have D = X. */

	  traceGtX = SCALAR_RE(evectmatrix_traceXtY(G, X));
	  if (usingConjugateGradient) {
               real gamma, gamma_numerator;

               if (use_polak_ribiere) {
                    /* assign G = G - prev_G and copy prev_G = G in the
                       same loop.  We can't use the BLAS routines because
                       we would then need an extra n x p array. */
                    for (i = 0; i < Y.n * Y.p; ++i) {
                         scalar g = G.data[i];
                         ACCUMULATE_DIFF(G.data[i], prev_G.data[i]);
                         prev_G.data[i] = g;
                    }
                    gamma_numerator = SCALAR_RE(evectmatrix_traceXtY(G, X));
               }
               else /* otherwise, use Fletcher-Reeves (ignore prev_G) */
                    gamma_numerator = traceGtX;

               if (prev_traceGtX == 0.0)
                    gamma = 0.0;
               else
                    gamma = gamma_numerator / prev_traceGtX;

	       if ((flags & EIGS_RESET_CG) &&
                   (iteration + 1) % CG_RESET_ITERS == 0) {
		    /* periodically forget previous search directions,
		       and just juse D = X */
                    gamma = 0.0;
                    if (flags & EIGS_VERBOSE)
                         printf("    resetting CG direction...\n");
               }

               evectmatrix_aXpbY(gamma, D, 1.0, X);
          }

	  /* Minimize the trace along Y + lamba*D: */

	  if (use_linmin) {
	       real dE, d2E;
	       real d_norm2;
	       real improvement;

	  linmin:

	       A(D, G, Adata, 0, X); /* G = A D; X is scratch */
	       evectmatrix_XtX(DtD, D);
	       d_norm = sqrt(d_norm2 = SCALAR_RE(sqmatrix_trace(DtD)) / Y.p);
	       tfd.d_norm = d_norm;
	       evectmatrix_XtY(DtAD, D, G);
	       
	       evectmatrix_XtY(S1, Y, D);
	       sqmatrix_symmetrize(symYtD, S1);

	       evectmatrix_XtY(S1, Y, G);
	       sqmatrix_symmetrize(symYtAD, S1);

	       sqmatrix_AeBC(S1, U, 0, symYtD, 1);
	       dE = 2.0 * (SCALAR_RE(sqmatrix_traceAtB(U, symYtAD)) -
			   SCALAR_RE(sqmatrix_traceAtB(YtAYU, S1))) / d_norm;

	       sqmatrix_copy(S2, DtD);
	       sqmatrix_ApaBC(S2, -4.0, symYtD, 0, S1, 0);
	       sqmatrix_AeBC(S3, symYtAD, 0, S1, 0);
	       sqmatrix_AeBC(S1, U, 0, S2, 1);
	       d2E = 2.0 * (SCALAR_RE(sqmatrix_traceAtB(U, DtAD)) -
			    SCALAR_RE(sqmatrix_traceAtB(YtAYU, S1)) -
			    4.0 * SCALAR_RE(sqmatrix_traceAtB(U, S3))) 
		    / d_norm2;
	       
	       /* this is just Newton-Raphson to find a root of
		  the first derivative: */
	       theta = -dE/d2E;

	       if (d2E < 0) {
		    if (flags & EIGS_VERBOSE)
			 printf("    near maximum in trace\n");
		    theta = dE > 0 ? -fabs(prev_theta) : fabs(prev_theta);
	       }
	       else if (-0.5*dE*theta > 2.0 * fabs(E-prev_E)) {
		    if (flags & EIGS_VERBOSE)
			 printf("    large trace change predicted (%g%%)\n",
				-0.5*dE*theta/E * 100.0);
	       }
	       if (fabs(theta) >= K_PI) {
		    if (flags & EIGS_VERBOSE)
			 printf("    large theta (%g)\n", theta);
		    theta = dE > 0 ? -fabs(prev_theta) : fabs(prev_theta);
	       }

	       /* Set S1 to YtAYU * YtY = YtAY for use in linmin.
		  (tfd.YtAY == S1). */
	       sqmatrix_AeBC(S1, YtAYU, 0, YtY, 1);

	       theta = linmin(&improvement,
			      0, E, dE, dE > 0 ? -K_PI : K_PI, theta,
			      tolerance, trace_func, &tfd);

	       if (improvement > 0 && improvement < tolerance &&
		   !(flags & EIGS_ALWAYS_EXACT_LINMIN)) {
		    use_linmin = 0;
		    if (flags & EIGS_VERBOSE)
			 printf("    switching to approximate "
				"line minimization\n");
	       }

	       /* Shift Y to new location minimizing the trace along D: */
	       evectmatrix_aXpbY(cos(theta), Y, sin(theta) / d_norm, D);
	  }
	  else { /* approximate line minimization */
	       real dE, E2, d2E, t;

	       /* Here, we do an approximate line minimization along D
		  by evaluating dE (the derivative) at the current point,
		  and the trace E2 at a second point, and then approximating
		  the second derivative d2E by finite differences.  Then,
		  we use one step of Newton's method on the derivative. */

	       d_norm = sqrt(SCALAR_RE(evectmatrix_traceXtY(D,D)) / Y.p);

	       /* dE = 2 * tr Gt D.  (Use prev_G instead of G so that
		  it works even when we are using Polak-Ribiere.) */
	       dE = 2.0 * SCALAR_RE(evectmatrix_traceXtY(prev_G, D)) / d_norm;

	       /* shift Y by prev_theta along D, in the downhill direction: */
	       t = dE < 0 ? -fabs(prev_theta) : fabs(prev_theta);
	       evectmatrix_aXpbY(1.0, Y, t / d_norm, D);

	       evectmatrix_XtX(U, Y);
	       sqmatrix_invert(U);  /* U = 1 / (Yt Y) */
	       A(Y, X, Adata, 1, G); /* X = AY; G is scratch */
	       evectmatrix_XtY(S1, Y, X);  /* S1 = Yt A Y */
	       
	       E2 = SCALAR_RE(sqmatrix_traceAtB(S1, U));

	       /* Get finite-difference approximation for the 2nd derivative
		  of the trace.  Equivalently, fit to a quadratic of the
		  form:  E(theta) = E + dE theta + 1/2 d2E theta^2 */
	       d2E = (E2 - E - dE * t) / (0.5 * t * t);

	       theta = -dE/d2E;

	       /* If the 2nd derivative is negative, or a big shift
		  in the trace is predicted (compared to the previous
		  iteration), then this approximate line minimization is
		  probably not very good; switch back to the exact
		  line minimization.  Hopefully, we won't have to
		  abort like this very often, as it wastes operations. */
	       if (d2E < 0 || -0.5*dE*theta > 10.0 * fabs(E-prev_E)) {
		    if (flags & EIGS_VERBOSE)
			 printf("    switching back to exact "
				"line minimization\n");
		    use_linmin = 1;
		    evectmatrix_aXpbY(1.0, Y, -t / d_norm, D);
		    goto linmin;
	       }

	       /* Shift Y to new location, hopefully minimizing the trace: */
               evectmatrix_aXpbY(1.0, Y, (theta - t) / d_norm, D);
	  }

	  if (constraint)
               constraint(Y, constraint_data);

	  prev_traceGtX = traceGtX;
          prev_theta = theta;
          prev_E = E;
     } while (++iteration < EIGENSOLVER_MAX_ITERATIONS);

     CHECK(iteration < EIGENSOLVER_MAX_ITERATIONS,
           "failure to converge after "
           STRINGIZE(EIGENSOLVER_MAX_ITERATIONS)
           " iterations");

     eigensolver_get_eigenvals_aux(Y, eigenvals, A, Adata,
				   X, G, U, S1, S2);

     *num_iterations = iteration;
     
     destroy_sqmatrix(S3);
     destroy_sqmatrix(S2);
     destroy_sqmatrix(S1);
     destroy_sqmatrix(symYtD);
     destroy_sqmatrix(DtD);
     destroy_sqmatrix(U);
     destroy_sqmatrix(YtY);
     destroy_sqmatrix(symYtAD);
     destroy_sqmatrix(DtAD);
     destroy_sqmatrix(YtAYU);
}
