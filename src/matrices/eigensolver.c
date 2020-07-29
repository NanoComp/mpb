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
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <scalar.h>
#include <matrices.h>
#include <blasglue.h>

#include "eigensolver.h"
#include "linmin.h"
#include "verbosity.h"

extern void eigensolver_get_eigenvals_aux(evectmatrix Y, real *eigenvals, evectoperator A,
                                          void *Adata, evectmatrix Work1, evectmatrix Work2,
                                          sqmatrix U, sqmatrix Usqrt, sqmatrix Uwork);

#define STRINGIZEx(x) #x /* a hack so that we can stringize macro values */
#define STRINGIZE(x) STRINGIZEx(x)

#define K_PI 3.141592653589793238462643383279502884197
#define MIN2(a, b) ((a) < (b) ? (a) : (b))
#define MAX2(a, b) ((a) > (b) ? (a) : (b))

#if defined(SCALAR_LONG_DOUBLE_PREC)
#define fabs fabsl
#define cos cosl
#define sin sinl
#define sqrt sqrtl
#define atan atanl
#define atan2 atan2l
#endif

/* Evalutate op, and set t to the elapsed time (in seconds). */
#define TIME_OP(t, op)                                                                             \
  {                                                                                                \
    mpiglue_clock_t xxx_time_op_start_time = MPIGLUE_CLOCK;                                        \
    { op; }                                                                                        \
    (t) = MPIGLUE_CLOCK_DIFF(MPIGLUE_CLOCK, xxx_time_op_start_time);                               \
  }

/**************************************************************************/

#define EIGENSOLVER_MAX_ITERATIONS 100000
#define FEEDBACK_TIME 4.0 /* elapsed time before we print progress feedback */

/* Number of iterations after which to reset conjugate gradient
   direction to steepest descent.  (Picked after some experimentation.
   Is there a better basis?  Should this change with the problem
   size?) */
#define CG_RESET_ITERS 70

/* Threshold for trace(1/YtBY) = trace(U) before we reorthogonalize: */
#define EIGS_TRACE_U_THRESHOLD 1e8

/**************************************************************************/

/* estimated times/iteration for different iteration schemes, based
   on the measure times for various operations and the operation counts: */

#define EXACT_LINMIN_TIME(t_AZ, t_KZ, t_ZtW, t_ZS, t_ZtZ, t_linmin)                                \
  ((t_AZ)*2 + (t_KZ) + (t_ZtW)*4 + (t_ZS)*2 + (t_ZtZ)*2 + (t_linmin))

#define APPROX_LINMIN_TIME(t_AZ, t_KZ, t_ZtW, t_ZS, t_ZtZ)                                         \
  ((t_AZ)*2 + (t_KZ) + (t_ZtW)*2 + (t_ZS)*2 + (t_ZtZ)*2)

/* Guess for the convergence slowdown factor due to the approximate
   line minimization.  It is probably best to be conservative, as the
   exact line minimization is more reliable and we only want to
   abandon it if there is a big speed gain. */
#define APPROX_LINMIN_SLOWDOWN_GUESS 2.0

/* We also don't want to use the approximate line minimization if
   the exact line minimization makes a big difference in the value
   of the trace that's achieved (i.e. if one step of Newton's method
   on the trace derivative does not do a good job).  The following
   is the maximum improvement by the exact line minimization (over
   one step of Newton) at which we'll allow the use of approximate line
   minimization. */
#define APPROX_LINMIN_IMPROVEMENT_THRESHOLD 0.05

/**************************************************************************/

typedef struct {
  sqmatrix YtAY, DtAD, symYtAD, YtBY, DtBD, symYtBD, S1, S2, S3;
  real lag, d_lag, trace_YtLY, trace_DtLD, trace_YtLD;
} trace_func_data;

static linmin_real trace_func(linmin_real theta, linmin_real *trace_deriv, void *data) {
  linmin_real trace;
  trace_func_data *d = (trace_func_data *)data;
  linmin_real c = cos(theta), s = sin(theta);

  {
    sqmatrix_copy(d->S1, d->YtBY);
    sqmatrix_aApbB(c * c, d->S1, s * s, d->DtBD);
    sqmatrix_ApaB(d->S1, 2 * s * c, d->symYtBD);
    if (!sqmatrix_invert(d->S1, 1, d->S2)) {
      /* if c or s is small, we sometimes run into numerical
         difficulties, and it is better to use the Taylor expansion */
      if (c < 1e-4 * s && c != 0) {
        sqmatrix_copy(d->S1, d->DtBD);
        CHECK(sqmatrix_invert(d->S1, 1, d->S2), "singular DtBD!");
        sqmatrix_AeBC(d->S2, d->S1, 0, d->symYtBD, 1);
        sqmatrix_AeBC(d->S3, d->S2, 0, d->S1, 1);
        sqmatrix_aApbB(1 / (s * s), d->S1, -2 * c / (s * s * s), d->S3);
      }
      else if (s < 1e-4 * c && s != 0) {
        sqmatrix_copy(d->S1, d->YtBY);
        CHECK(sqmatrix_invert(d->S1, 1, d->S2), "singular DtBD!");
        sqmatrix_AeBC(d->S2, d->S1, 0, d->symYtBD, 1);
        sqmatrix_AeBC(d->S3, d->S2, 0, d->S1, 1);
        sqmatrix_aApbB(1 / (c * c), d->S1, -2 * s / (c * c * c), d->S3);
      }
      else {
        CHECK(0, "inexplicable singularity in linmin trace_func");
      }
    }

    sqmatrix_copy(d->S2, d->YtAY);
    sqmatrix_aApbB(c * c, d->S2, s * s, d->DtAD);
    sqmatrix_ApaB(d->S2, 2 * s * c, d->symYtAD);

    trace = SCALAR_RE(sqmatrix_traceAtB(d->S2, d->S1)) +
            (c * c * d->trace_YtLY + s * s * d->trace_DtLD + 2 * s * c * d->trace_YtLD) *
                (c * d->lag + s * d->d_lag);
  }

  if (trace_deriv) {
    linmin_real c2 = cos(2 * theta), s2 = sin(2 * theta);

    sqmatrix_copy(d->S3, d->YtAY);
    sqmatrix_ApaB(d->S3, -1.0, d->DtAD);
    sqmatrix_aApbB(-0.5 * s2, d->S3, c2, d->symYtAD);

    *trace_deriv = SCALAR_RE(sqmatrix_traceAtB(d->S1, d->S3));

    sqmatrix_AeBC(d->S3, d->S1, 0, d->S2, 1);
    sqmatrix_AeBC(d->S2, d->S3, 0, d->S1, 1);

    sqmatrix_copy(d->S3, d->YtBY);
    sqmatrix_ApaB(d->S3, -1.0, d->DtBD);
    sqmatrix_aApbB(-0.5 * s2, d->S3, c2, d->symYtBD);

    *trace_deriv -= SCALAR_RE(sqmatrix_traceAtB(d->S2, d->S3));
    *trace_deriv *= 2;

    *trace_deriv += (-s2 * d->trace_YtLY + s2 * d->trace_DtLD + 2 * c2 * d->trace_YtLD) *
                    (c * d->lag + s * d->d_lag);
    *trace_deriv += (c * c * d->trace_YtLY + s * s * d->trace_DtLD + 2 * s * c * d->trace_YtLD) *
                    (-s * d->lag + c * d->d_lag);
  }

  return trace;
}

/**************************************************************************/

#define EIG_HISTORY_SIZE 5

/* find generalized eigenvectors Y of (A,B) by minimizing Rayleigh quotient

        tr [ Yt A Y / (Yt B Y) ] + lag * tr [ Yt L Y ]

   where lag is a Lagrange multiplier and L is a Hermitian operator
   implementing some constraint tr [ Yt L Y ] = 0 on the eigenvectors
   (if L is not NULL).

   Constraints that commute with A and B (and L) are specified via the
   "constraint" argument, which gives the projection operator for
   the constraint(s).
*/

void eigensolver_lagrange(evectmatrix Y, real *eigenvals, evectoperator A, void *Adata,
                          evectoperator B, void *Bdata, evectpreconditioner K, void *Kdata,
                          evectconstraint constraint, void *constraint_data, evectoperator L,
                          void *Ldata, real *lag, evectmatrix Work[], int nWork, real tolerance,
                          int *num_iterations, int flags) {
  real convergence_history[EIG_HISTORY_SIZE];
  evectmatrix G, D, X, BY, prev_G, BD;
  real g_lag = 0, d_lag = 0, prev_g_lag = 0;
  short usingConjugateGradient = 0, use_polak_ribiere = 0, use_linmin = 1;
  real E, prev_E = 0.0;
  real d_scale = 1.0;
  real traceGtX, prev_traceGtX = 0.0;
  real theta, prev_theta = 0.5;
  int i, iteration = 0, num_emergency_restarts = 0;
  mpiglue_clock_t prev_feedback_time;
  real time_AZ, time_KZ = 0, time_ZtZ, time_ZtW, time_ZS, time_linmin = 0;
  real linmin_improvement = 0;
  sqmatrix YtAYU, DtAD, symYtAD, YtBY, U, DtBD, symYtBD, S1, S2, S3;
  trace_func_data tfd;

  prev_feedback_time = MPIGLUE_CLOCK;

#ifdef DEBUG
  flags |= EIGS_VERBOSE;
#endif

  CHECK(nWork >= 2, "not enough workspace");
  G = Work[0];
  X = Work[1];

  if (B) {
    CHECK(nWork >= 3, "not enough workspace for generalized problem");
    BY = Work[2];
  }
  else
    BY = Y;

  usingConjugateGradient = nWork >= 3 + (B != NULL);
  if (usingConjugateGradient) {
    D = Work[2 + (B != NULL)];
    for (i = 0; i < D.n * D.p; ++i)
      ASSIGN_ZERO(D.data[i]);
  }
  else
    D = X;

  BD = B ? BY : D; /* storage for B*D (re-use B*Y storage) */

  use_polak_ribiere = nWork >= 4 + (B != NULL);
  if (use_polak_ribiere) {
    prev_G = Work[3 + (B != NULL)];
    for (i = 0; i < Y.n * Y.p; ++i)
      ASSIGN_ZERO(prev_G.data[i]);
    if (flags & EIGS_ORTHOGONAL_PRECONDITIONER) /* see below */
      fprintf(stderr, "WARNING: Polak-Ribiere may not work with the "
                      "orthogonal-preconditioner option.\n");
  }
  else
    prev_G = G;

  YtAYU = create_sqmatrix(Y.p);   /* holds Yt A Y */
  DtAD = create_sqmatrix(Y.p);    /* holds Dt A D */
  symYtAD = create_sqmatrix(Y.p); /* holds (Yt A D + Dt A Y) / 2 */
  YtBY = create_sqmatrix(Y.p);    /* holds Yt B Y */
  U = create_sqmatrix(Y.p);       /* holds 1 / (Yt Y) */
  DtBD = create_sqmatrix(Y.p);    /* holds Dt B D */
  symYtBD = create_sqmatrix(Y.p); /* holds (Yt B D + Dt B Y) / 2 */

  /* Notation note: "t" represents a dagger superscript, so
     Yt represents adjoint(Y), or Y' in MATLAB syntax. */

  /* scratch matrices: */
  S1 = create_sqmatrix(Y.p);
  S2 = create_sqmatrix(Y.p);
  S3 = create_sqmatrix(Y.p);

  tfd.YtAY = S1;
  tfd.DtAD = DtAD;
  tfd.symYtAD = symYtAD;
  tfd.YtBY = YtBY;
  tfd.DtBD = DtBD;
  tfd.symYtBD = symYtBD;
  tfd.S1 = YtAYU;
  tfd.S2 = S2;
  tfd.S3 = S3;

restartY:

  if (flags & EIGS_ORTHONORMALIZE_FIRST_STEP) {
    if (B) {
      B(Y, BY, Bdata, 1, G); /* B*Y; G is scratch */
      evectmatrix_XtY(U, Y, BY, S2);
    }
    else {
      evectmatrix_XtX(U, Y, S2);
    }
    sqmatrix_assert_hermitian(U);
    CHECK(sqmatrix_invert(U, 1, S2), "non-independent initial Y");
    sqmatrix_sqrt(S1, U, S2);      /* S1 = 1/sqrt(Yt*Y) */
    evectmatrix_XeYS(G, Y, S1, 1); /* G = orthonormalize Y */
    evectmatrix_copy(Y, G);
  }

  for (i = 0; i < Y.p; ++i)
    eigenvals[i] = 0.0;

  for (i = 0; i < EIG_HISTORY_SIZE; ++i)
    convergence_history[i] = 10000.0;

  if (constraint) constraint(Y, constraint_data);

  do {
    real y_norm, gamma_numerator = 0;

    if (flags & EIGS_FORCE_APPROX_LINMIN) use_linmin = 0;

    if (B) {
      B(Y, BY, Bdata, 1, G); /* B*Y; G is scratch */
      TIME_OP(time_ZtZ, evectmatrix_XtY(YtBY, Y, BY, S2));
    }
    else {
      TIME_OP(time_ZtZ, evectmatrix_XtX(YtBY, Y, S2));
    }
    sqmatrix_assert_hermitian(YtBY);

    y_norm = sqrt(SCALAR_RE(sqmatrix_trace(YtBY)) / Y.p);
    blasglue_rscal(Y.p * Y.n, 1 / y_norm, Y.data, 1);
    if (B) blasglue_rscal(Y.p * Y.n, 1 / y_norm, BY.data, 1);
    blasglue_rscal(Y.p * Y.p, 1 / (y_norm * y_norm), YtBY.data, 1);

    sqmatrix_copy(U, YtBY);
    if (!sqmatrix_invert(U, 1, S2)) { /* non-independent Y columns */
      /* emergency restart with random Y */
      CHECK(iteration + 10 * ++num_emergency_restarts < EIGENSOLVER_MAX_ITERATIONS,
            "too many emergency restarts");
      if (mpb_verbosity >= 2) {
        mpi_one_printf("    emergency randomization of Y on iter. %d\n", iteration);
      }
      for (i = 0; i < Y.p * Y.n; ++i)
        ASSIGN_SCALAR(Y.data[i], rand() * 1.0 / RAND_MAX - 0.5, rand() * 1.0 / RAND_MAX - 0.5);
      goto restartY;
    }

    /* If trace(1/YtBY) gets big, it means that the columns
       of Y are becoming nearly parallel.  This sometimes happens,
       especially in the targeted eigensolver, because the
       preconditioner pushes all the columns towards the ground
       state.  If it gets too big, it seems to be a good idea
       to re-orthogonalize, resetting conjugate-gradient, as
       otherwise we start to encounter numerical problems. */
    if (flags & EIGS_REORTHOGONALIZE) {
      real traceU = SCALAR_RE(sqmatrix_trace(U));
      mpi_assert_equal(traceU);
      if (traceU > EIGS_TRACE_U_THRESHOLD * U.p) {
        if (mpb_verbosity >= 2) { mpi_one_printf("    re-orthonormalizing Y\n"); }
        sqmatrix_sqrt(S1, U, S2);      /* S1 = 1/sqrt(Yt*Y) */
        evectmatrix_XeYS(G, Y, S1, 1); /* G = orthonormalize Y */
        evectmatrix_copy(Y, G);
        prev_traceGtX = 0.0;
        if (B) {
          B(Y, BY, Bdata, 1, G); /* B*Y; G is scratch */
          evectmatrix_XtY(YtBY, Y, BY, S2);
        }
        else
          evectmatrix_XtX(YtBY, Y, S2);
        y_norm = sqrt(SCALAR_RE(sqmatrix_trace(YtBY)) / Y.p);
        blasglue_rscal(Y.p * Y.n, 1 / y_norm, Y.data, 1);
        if (B) blasglue_rscal(Y.p * Y.n, 1 / y_norm, BY.data, 1);
        blasglue_rscal(Y.p * Y.p, 1 / (y_norm * y_norm), YtBY.data, 1);
        sqmatrix_copy(U, YtBY);
        CHECK(sqmatrix_invert(U, 1, S2), "non-independent Y after re-orthogonalization");
      }
    }

    TIME_OP(time_AZ, A(Y, X, Adata, 1, G)); /* X = AY; G is scratch */

#ifdef DEBUG
    evectmatrix_XtY(S1, Y, X, S2);
    sqmatrix_assert_hermitian(S1);
#endif

    /* G = AYU; note that U is Hermitian: */
    TIME_OP(time_ZS, evectmatrix_XeYS(G, X, U, 1));

    TIME_OP(time_ZtW, evectmatrix_XtY(YtAYU, Y, G, S2));
    E = SCALAR_RE(sqmatrix_trace(YtAYU));
    CHECK(!BADNUM(E), "crazy number detected in trace!!\n");
    mpi_assert_equal(E);

    if (L) {
      /* X = LY, no scratch */
      L(Y, X, Ldata, 1, X);
      g_lag = tfd.trace_YtLY = SCALAR_RE(evectmatrix_traceXtY(Y, X));
      E += *lag * g_lag;
    }

    convergence_history[iteration % EIG_HISTORY_SIZE] =
        200.0 * fabs(E - prev_E) / (fabs(E) + fabs(prev_E));

    if (iteration > 0 && mpi_is_master() &&
        ((flags & EIGS_VERBOSE) ||
         MPIGLUE_CLOCK_DIFF(MPIGLUE_CLOCK, prev_feedback_time) > FEEDBACK_TIME)) {
      if (mpb_verbosity >= 2) {
        mpi_one_printf("    iteration %4d: "
                       "trace = %0.16g (%g%% change)\n",
                       iteration, (double)E,
                       (double)convergence_history[iteration % EIG_HISTORY_SIZE]);
      }
      if (flags & EIGS_VERBOSE) debug_output_malloc_count();
      fflush(stdout);                     /* make sure output appears */
      prev_feedback_time = MPIGLUE_CLOCK; /* reset feedback clock */
    }

    if (iteration > 0 && fabs(E - prev_E) < tolerance * 0.5 * (E + prev_E + 1e-7))
      break; /* convergence!  hooray! */

    /* Compute gradient of functional: G = (1 - BY U Yt) A Y U */
    sqmatrix_AeBC(S1, U, 0, YtAYU, 0);
    evectmatrix_XpaYS(G, -1.0, BY, S1, 1);

    if (L) { /* include Lagrange gradient; note X = LY from above */
      evectmatrix_aXpbY(1.0, G, *lag, X);
    }

    /* set X = precondition(G): */
    if (K != NULL) {
      TIME_OP(time_KZ, K(G, X, Kdata, Y, NULL, YtBY));
      /* Note: we passed NULL for eigenvals since we haven't
         diagonalized YAY (nor are the Y's orthonormal). */
    }
    else
      evectmatrix_copy(X, G); /* preconditioner is identity */

    /* We have to apply the constraint here, in case it doesn't
       commute with the preconditioner. */
    if (constraint) constraint(X, constraint_data);

    if (flags & EIGS_PROJECT_PRECONDITIONING) {
      /* Operate projection P = (1 - BY U Yt) on X: */
      evectmatrix_XtY(symYtBD, Y, X, S2); /* symYtBD = Yt X */
      sqmatrix_AeBC(S1, U, 0, symYtBD, 0);
      evectmatrix_XpaYS(X, -1.0, BY, S1, 0);
    }

    /* Now, for the case of EIGS_ORTHOGONAL_PRECONDITIONER, we
       need to use G as scratch space in order to avoid the need
       for an extra column bundle.  Before that, we need to do
       any computations that we need with G.  (Yes, we're
       playing tricksy games here, but isn't it fun?) */

    mpi_assert_equal(traceGtX = SCALAR_RE(evectmatrix_traceXtY(G, X)) + g_lag * g_lag);
    if (usingConjugateGradient) {
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

        {
          real g = g_lag;
          g_lag -= prev_g_lag;
          prev_g_lag = g;
        }
        gamma_numerator += g_lag * prev_g_lag;
      }
      else /* otherwise, use Fletcher-Reeves (ignore prev_G) */
        gamma_numerator = traceGtX;
      mpi_assert_equal(gamma_numerator);
    }

    /* The motivation for the following code came from a trick I
       noticed in Sleijpen and Van der Vorst, "A Jacobi-Davidson
       iteration method for linear eigenvalue problems," SIAM
       J. Matrix Anal. Appl. 17, 401-425 (April 1996).  (The
       motivation in our case comes from the fact that if you
       look at the Hessian matrix of the problem, it has a
       projection operator just as in the above reference, and
       so we should the same technique to invert it.)  So far,
       though, the hoped-for savings haven't materialized; maybe
       we need a better preconditioner first. */
    if (flags & EIGS_ORTHOGONAL_PRECONDITIONER) {
      real traceGtX_delta; /* change in traceGtX when we update X */

      /* set G = precondition(Y): */
      if (K != NULL)
        K(Y, G, Kdata, Y, NULL, YtBY);
      else
        evectmatrix_copy(G, Y); /* preconditioner is identity */

      /* let X = KG - KY S3t, where S3 is chosen so that YtX = 0:
         S3 = (YtKG)t / (YtKY).  Recall that, at this point,
         X holds KG and G holds KY.  K is assumed Hermitian. */
      evectmatrix_XtY(S1, Y, G, S2);
      CHECK(sqmatrix_invert(S1, 0, S2), "singular YtKY"); /* S1 = 1 / (YtKY) */
      evectmatrix_XtY(S2, X, Y, S3);                      /* S2 = GtKY = (YtKG)t */
      sqmatrix_AeBC(S3, S2, 0, S1, 1);
      evectmatrix_XpaYS(X, -1.0, G, S3, 1);

      /* Update traceGtX and gamma_numerator.  The update
         for gamma_numerator isn't really right in the case
         of Polak-Ribiere; it amounts to doing a weird combination
         of P-R and Fletcher-Reeves...what will happen?  (To
         do the right thing, I think we would need an extra
         column bundle.)  */
      traceGtX_delta = -SCALAR_RE(sqmatrix_traceAtB(S3, S2));
      traceGtX += traceGtX_delta;
      if (usingConjugateGradient) gamma_numerator += traceGtX_delta;
    }

    /* In conjugate-gradient, the minimization direction D is
       a combination of X with the previous search directions.
       Otherwise, we just have D = X. */

    if (usingConjugateGradient) {
      real gamma;

      if (prev_traceGtX == 0.0)
        gamma = 0.0;
      else
        gamma = gamma_numerator / prev_traceGtX;

      if ((flags & EIGS_DYNAMIC_RESET_CG) &&
          2.0 * convergence_history[iteration % EIG_HISTORY_SIZE] >=
              convergence_history[(iteration + 1) % EIG_HISTORY_SIZE]) {
        gamma = 0.0;
        if (flags & EIGS_VERBOSE && mpb_verbosity >= 2)
          mpi_one_printf("    dynamically resetting CG direction...\n");
        for (i = 1; i < EIG_HISTORY_SIZE; ++i)
          convergence_history[(iteration + i) % EIG_HISTORY_SIZE] = 10000.0;
      }

      if ((flags & EIGS_RESET_CG) && (iteration + 1) % CG_RESET_ITERS == 0) {
        /* periodically forget previous search directions,
           and just juse D = X */
        gamma = 0.0;
        if (flags & EIGS_VERBOSE && mpb_verbosity >= 2)
          mpi_one_printf("    resetting CG direction...\n");
      }

      mpi_assert_equal(gamma * d_scale);
      evectmatrix_aXpbY(gamma * d_scale, D, 1.0, X);
      d_lag = gamma * d_scale * d_lag + g_lag;
    }

    d_scale = 1.0;

    /* Minimize the trace along Y + lambda*D: */

    if (!use_linmin) {
      real dE, E2, d2E, t, d_norm;

      /* Here, we do an approximate line minimization along D
         by evaluating dE (the derivative) at the current point,
         and the trace E2 at a second point, and then approximating
         the second derivative d2E by finite differences.  Then,
         we use one step of Newton's method on the derivative.
         This has the advantage of requiring two fewer O(np^2)
         matrix multiplications compared to the exact linmin. */

      if (B) B(D, BD, Bdata, 0, BD); /* B*Y; no scratch */

      d_norm = sqrt(SCALAR_RE(evectmatrix_traceXtY(BD, D)) / Y.p);
      mpi_assert_equal(d_norm);

      /* dE = 2 * tr Gt D.  (Use prev_G instead of G so that
         it works even when we are using Polak-Ribiere.) */
      dE = 2.0 * SCALAR_RE(evectmatrix_traceXtY(prev_G, D)) / d_norm;

      /* shift Y by prev_theta along D, in the downhill direction: */
      t = dE < 0 ? -fabs(prev_theta) : fabs(prev_theta);
      evectmatrix_aXpbY(1.0, Y, t / d_norm, D);

      if (B) {
        B(Y, BY, Bdata, 1, G); /* B*Y; G is scratch */
        evectmatrix_XtY(U, Y, BY, S2);
      }
      else
        evectmatrix_XtX(U, Y, S2);
      CHECK(sqmatrix_invert(U, 1, S2), "singular YtBY"); /* U = 1 / (Yt B Y) */
      A(Y, G, Adata, 1, X);                              /* G = AY; X is scratch */
      evectmatrix_XtY(S1, Y, G, S2);                     /* S1 = Yt A Y */

      E2 = SCALAR_RE(sqmatrix_traceAtB(S1, U));

      if (L) {
        *lag += (t / d_norm) * d_lag;
        L(Y, X, Ldata, 1, X);
        E2 += *lag * SCALAR_RE(evectmatrix_traceXtY(Y, X));
      }

      mpi_assert_equal(E2);

      /* Get finite-difference approximation for the 2nd derivative
         of the trace.  Equivalently, fit to a quadratic of the
         form:  E(theta) = E + dE theta + 1/2 d2E theta^2 */
      d2E = (E2 - E - dE * t) / (0.5 * t * t);

      theta = -dE / d2E;

      /* If the 2nd derivative is negative, or a big shift
         in the trace is predicted (compared to the previous
         iteration), then this approximate line minimization is
         probably not very good; switch back to the exact
         line minimization.  Hopefully, we won't have to
         abort like this very often, as it wastes operations. */
      if (d2E < 0 || -0.5 * dE * theta > 20.0 * fabs(E - prev_E)) {
        if (flags & EIGS_FORCE_APPROX_LINMIN) {
          if (flags & EIGS_VERBOSE && mpb_verbosity >= 2)
            mpi_one_printf("    using previous stepsize\n");
        }
        else {
          if (flags & EIGS_VERBOSE && mpb_verbosity >= 2)
            mpi_one_printf("    switching back to exact "
                           "line minimization\n");
          use_linmin = 1;
          evectmatrix_aXpbY(1.0, Y, -t / d_norm, D);
          if (L) *lag -= (t / d_norm) * d_lag;
          prev_theta = atan(prev_theta); /* convert to angle */
          /* don't do this again: */
          flags |= EIGS_FORCE_EXACT_LINMIN;
        }
      }
      else {
        /* Shift Y by theta, hopefully minimizing the trace: */
        evectmatrix_aXpbY(1.0, Y, (theta - t) / d_norm, D);
        if (L) *lag += ((theta - t) / d_norm) * d_lag;
      }
    }
    if (use_linmin) {
      real dE, d2E;

      if (B) B(D, BD, Bdata, 0, G); /* B*Y; G is scratch */

      d_scale = sqrt(SCALAR_RE(evectmatrix_traceXtY(BD, D)) / Y.p);
      mpi_assert_equal(d_scale);
      blasglue_rscal(Y.p * Y.n, 1 / d_scale, D.data, 1);
      if (B) blasglue_rscal(Y.p * Y.n, 1 / d_scale, BD.data, 1);

      A(D, G, Adata, 0, X); /* G = A D; X is scratch */
      if (B)
        evectmatrix_XtY(DtBD, D, BD, S2);
      else
        evectmatrix_XtX(DtBD, D, S2);
      evectmatrix_XtY(DtAD, D, G, S2);
      sqmatrix_assert_hermitian(DtBD);
      sqmatrix_assert_hermitian(DtAD);

      evectmatrix_XtY(S1, Y, BD, S2);
      sqmatrix_symmetrize(symYtBD, S1);

      evectmatrix_XtY(S1, Y, G, S2);
      sqmatrix_symmetrize(symYtAD, S1);

      sqmatrix_AeBC(S1, U, 0, symYtBD, 1);
      dE = 2.0 *
           (SCALAR_RE(sqmatrix_traceAtB(U, symYtAD)) - SCALAR_RE(sqmatrix_traceAtB(YtAYU, S1)));

      sqmatrix_copy(S2, DtBD);
      sqmatrix_ApaBC(S2, -4.0, symYtBD, 0, S1, 0);
      sqmatrix_AeBC(S3, symYtAD, 0, S1, 0);
      sqmatrix_AeBC(S1, U, 0, S2, 1);
      d2E = 2.0 * (SCALAR_RE(sqmatrix_traceAtB(U, DtAD)) - SCALAR_RE(sqmatrix_traceAtB(YtAYU, S1)) -
                   4.0 * SCALAR_RE(sqmatrix_traceAtB(U, S3)));

      if (L) {
        d_lag *= 1 / d_scale;
        tfd.d_lag = d_lag;
        tfd.lag = *lag;
        /* note: tfd.trace_YtLY was set above */
        L(D, X, Ldata, 0, X);
        tfd.trace_DtLD = SCALAR_RE(evectmatrix_traceXtY(D, X));
        tfd.trace_YtLD = SCALAR_RE(evectmatrix_traceXtY(Y, X));
        dE += tfd.lag * 2.0 * tfd.trace_YtLD + tfd.d_lag * tfd.trace_YtLY;
        d2E += tfd.lag * 2.0 * tfd.trace_DtLD + tfd.d_lag * 4.0 * tfd.trace_YtLD;
      }
      else {
        tfd.d_lag = tfd.lag = tfd.trace_YtLY = tfd.trace_DtLD = tfd.trace_YtLD = 0;
      }

      /* this is just Newton-Raphson to find a root of
         the first derivative: */
      theta = -dE / d2E;

      if (d2E < 0) {
        if (flags & EIGS_VERBOSE && mpb_verbosity >= 2)
          mpi_one_printf("    near maximum in trace\n");
        theta = dE > 0 ? -fabs(prev_theta) : fabs(prev_theta);
      }
      else if (-0.5 * dE * theta > 2.0 * fabs(E - prev_E)) {
        if (flags & EIGS_VERBOSE && mpb_verbosity >= 2)
          mpi_one_printf("    large trace change predicted "
                         "(%g%%)\n",
                         (double)(-0.5 * dE * theta / E * 100.0));
      }
      if (fabs(theta) >= K_PI) {
        if (flags & EIGS_VERBOSE && mpb_verbosity >= 2)
          mpi_one_printf("    large theta (%g)\n", (double)theta);
        theta = dE > 0 ? -fabs(prev_theta) : fabs(prev_theta);
      }

      /* Set S1 to YtAYU * YtBY = YtAY for use in linmin.
         (tfd.YtAY == S1). */
      sqmatrix_AeBC(S1, YtAYU, 0, YtBY, 1);
      sqmatrix_assert_hermitian(S1);

      mpi_assert_equal(theta);
      {
        linmin_real new_E, new_dE;
        TIME_OP(time_linmin,
                theta = linmin(&new_E, &new_dE, theta, E, dE, 0.1, MIN2(tolerance, 1e-6), 1e-14, 0,
                               dE > 0 ? -K_PI : K_PI, trace_func, &tfd, flags & EIGS_VERBOSE));
        linmin_improvement = fabs(E - new_E) * 2.0 / fabs(E + new_E);
      }
      mpi_assert_equal(theta);

      CHECK(fabs(theta) <= K_PI, "converged out of bounds!");

      /* Shift Y to new location minimizing the trace along D: */
      evectmatrix_aXpbY(cos(theta), Y, sin(theta), D);
      if (L) *lag = *lag * cos(theta) + d_lag * sin(theta);
    }

    /* In exact arithmetic, we don't need to do this, but in practice
       it is probably a good idea to keep errors from adding up and
       eventually violating the constraints. */
    if (constraint) constraint(Y, constraint_data);

    prev_traceGtX = traceGtX;
    prev_theta = theta;
    prev_E = E;

    /* Finally, we use the times for the various operations to
       help us pick an algorithm for the next iteration: */
    {
      real t_exact, t_approx;
      t_exact = EXACT_LINMIN_TIME(time_AZ, time_KZ, time_ZtW, time_ZS, time_ZtZ, time_linmin);
      t_approx = APPROX_LINMIN_TIME(time_AZ, time_KZ, time_ZtW, time_ZS, time_ZtZ);
      if (flags & EIGS_PROJECT_PRECONDITIONING) {
        t_exact += time_ZtW + time_ZS;
        t_approx += time_ZtW + time_ZS;
      }

      /* Sum the times over the processors so that all the
         processors compare the same, average times. */
      mpi_allreduce_1(&t_exact, real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
      mpi_allreduce_1(&t_approx, real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);

      if (!(flags & EIGS_FORCE_EXACT_LINMIN) && linmin_improvement > 0 &&
          linmin_improvement <= APPROX_LINMIN_IMPROVEMENT_THRESHOLD &&
          t_exact > t_approx * APPROX_LINMIN_SLOWDOWN_GUESS) {
        if ((flags & EIGS_VERBOSE) && use_linmin && mpb_verbosity >= 2)
          mpi_one_printf("    switching to approximate "
                         "line minimization (decrease time by %g%%)\n",
                         (double)((t_exact - t_approx) * 100.0 / t_exact));
        use_linmin = 0;
      }
      else if (!(flags & EIGS_FORCE_APPROX_LINMIN)) {
        if ((flags & EIGS_VERBOSE) && !use_linmin && mpb_verbosity >= 2)
          mpi_one_printf("    switching back to exact "
                         "line minimization\n");
        use_linmin = 1;
        prev_theta = atan(prev_theta); /* convert to angle */
      }
    }
  } while (++iteration < EIGENSOLVER_MAX_ITERATIONS);

  CHECK(iteration < EIGENSOLVER_MAX_ITERATIONS,
        "failure to converge after " STRINGIZE(EIGENSOLVER_MAX_ITERATIONS) " iterations");

  if (B) {
    B(Y, BY, Bdata, 1, G); /* B*Y; G is scratch */
    evectmatrix_XtY(U, Y, BY, S2);
  }
  else
    evectmatrix_XtX(U, Y, S2);
  CHECK(sqmatrix_invert(U, 1, S2), "singular YtBY at end");
  eigensolver_get_eigenvals_aux(Y, eigenvals, A, Adata, X, G, U, S1, S2);

  *num_iterations = iteration;

  destroy_sqmatrix(S3);
  destroy_sqmatrix(S2);
  destroy_sqmatrix(S1);
  destroy_sqmatrix(symYtBD);
  destroy_sqmatrix(DtBD);
  destroy_sqmatrix(U);
  destroy_sqmatrix(YtBY);
  destroy_sqmatrix(symYtAD);
  destroy_sqmatrix(DtAD);
  destroy_sqmatrix(YtAYU);
}

void eigensolver(evectmatrix Y, real *eigenvals, evectoperator A, void *Adata, evectoperator B,
                 void *Bdata, evectpreconditioner K, void *Kdata, evectconstraint constraint,
                 void *constraint_data, evectmatrix Work[], int nWork, real tolerance,
                 int *num_iterations, int flags) {
  eigensolver_lagrange(Y, eigenvals, A, Adata, B, Bdata, K, Kdata, constraint, constraint_data, 0,
                       0, 0, Work, nWork, tolerance, num_iterations, flags);
}
