/* Copyright (C) 1996 Jorge J. More'.
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

/* The routines in this file were taken from the MINPACK-2 package by
   Jorge J. More', specifically its line search subroutines in the
   MINPACK-2/csrch, in the files dcsrch.f and dcstep.f.  MINPACK-2 is
   software for the solution of systems of nonlinear equations,
   nonlinear least squares problems, and minimization problems.  Web
   pages and ftp sites for MINPACK-2 can be found at:

        http://www.mcs.anl.gov/~more/minpack-2/minpack-2.html
        ftp://info.mcs.anl.gov/pub/MINPACK-2/
        http://www-fp.mcs.anl.gov/otc/minpack/summary.html

   It implements the line search algorithm described in:

   Jorge J. More and David J. Thuente, "Line search algorithms with
   guaranteed sufficient decrease," ACM Trans. on Mathematical
   Software, vol. 20, no. 3, pp. 286-307 (September 1994).

   The original code was under the copyright and license listed below,
   but Jorge J. More' graciously granted me permission to distribute
   under the terms of the GNU General Public License.

   Original copyright and license statement:

    * This program discloses material protectable under copyright laws of
    * the United States. Permission to copy and modify this software and
    * its documentation for internal research use is hereby granted,
    * provided that this notice is retained thereon and on all copies or
    * modifications.  The University of Chicago makes no representations
    * as to the suitability and operability of this software for any
    * purpose.  It is provided "as is" without express or implied
    * warranty.
    *
    * Use of this software for commercial purposes is expressly
    * prohibited without contacting.
    *
    * Jorge J. More'
    * Mathematics and Computer Science Division
    * Argonne National Laboratory
    * 9700 S. Cass Ave.
    * Argonne, Illinois 60439-4844
    * e-mail: more@mcs.anl.gov
    *
    * Argonne National Laboratory with facilities in the states of
    * Illinois and Idaho, is owned by The United States Government, and
    * operated by the University of Chicago under provision of a contract
    * with the Department of Energy.
*/

/* minpack2-linmin.f -- translated by f2c (version 19991025).
   C code cleaned up by Steven G. Johnson <stevenj@alum.mit.edu>.
*/

#include <math.h>
#include <string.h>

#include "config.h"

#include "linmin.h"
#define double linmin_real

/* Definitions so that we don't need -lf2c or f2c.h: */

typedef double doublereal;
typedef int integer;
typedef int logical;
typedef int ftnlen;

#ifndef HAVE_STRNCMP

/* provide a strncmp replacement if the system does not provide one: */
static int strncmp(const char *s1, const char *s2, size_t n) {
  size_t i;
  for (i = 0; i < n && s1[i] && s2[i] && s1[i] == s2[i]; ++i)
    ;
  if (i >= n)
    return 0;
  else
    return (s1[i] - s2[i]);
}

#endif /* ! HAVE_STRNCMP */

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#if defined(SCALAR_LONG_DOUBLE_PREC)
#define ABS(x) fabsl(x)
#define sqrt sqrtl
#else
#define ABS(x) fabs(x)
#endif
#define s_cmp(s1, s2, len1, len2) strncmp(s1, s2, MIN(len1, len2))
#define s_copy(s1, s2, len1, len2) strcpy(s1, s2)
#define TRUE_ 1
#define FALSE_ 0

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int dcsrch(doublereal *stp, doublereal *f, doublereal *g, doublereal *ftol,
                            doublereal *gtol, doublereal *xtol, char *task, doublereal *stpmin,
                            doublereal *stpmax, integer *isave, doublereal *dsave) {
  /* System generated locals */
  doublereal d__1;

  /* Builtin functions */

  /* Local variables */
  integer stage;
  doublereal finit, ginit, width, ftest, gtest, stmin, stmax, width1, fm, gm, fx, fy, gx, gy;
  logical brackt;
  extern /* Subroutine */ int dcstep(doublereal *, doublereal *, doublereal *, doublereal *,
                                     doublereal *, doublereal *, doublereal *, doublereal *,
                                     doublereal *, logical *, doublereal *, doublereal *);
  doublereal fxm, fym, gxm, gym, stx, sty;

  /*     ********** */

  /*     Subroutine dcsrch */

  /*     This subroutine finds a step that satisfies a sufficient */
  /*     decrease condition and a curvature condition. */

  /*     Each call of the subroutine updates an interval with */
  /*     endpoints stx and sty. The interval is initially chosen */
  /*     so that it contains a minimizer of the modified function */

  /*           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0). */

  /*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
  /*     interval is chosen so that it contains a minimizer of f. */

  /*     The algorithm is designed to find a step that satisfies */
  /*     the sufficient decrease condition */

  /*           f(stp) <= f(0) + ftol*stp*f'(0), */

  /*     and the curvature condition */

  /*           ABS(f'(stp)) <= gtol*ABS(f'(0)). */

  /*     If ftol is less than gtol and if, for example, the function */
  /*     is bounded below, then there is always a step which satisfies */
  /*     both conditions. */

  /*     If no step can be found that satisfies both conditions, then */
  /*     the algorithm stops with a warning. In this case stp only */
  /*     satisfies the sufficient decrease condition. */

  /*     A typical invocation of dcsrch has the following outline: */

  /*     Evaluate the function at stp = 0.0d0; store in f. */
  /*     Evaluate the gradient at stp = 0.0d0; store in g. */
  /*     Choose a starting step stp. */

  /*     task = 'START' */
  /*  10 continue */
  /*        call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax, */
  /*    +               isave,dsave) */
  /*        if (task .eq. 'FG') then */
  /*           Evaluate the function and the gradient at stp */
  /*           go to 10 */
  /*           end if */

  /*     NOTE: The user must not alter work arrays between calls. */

  /*     The subroutine statement is */

  /*       subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, */
  /*                         task,isave,dsave) */
  /*     where */

  /*       stp is a double precision variable. */
  /*         On entry stp is the current estimate of a satisfactory */
  /*            step. On initial entry, a positive initial estimate */
  /*            must be provided. */
  /*         On exit stp is the current estimate of a satisfactory step */
  /*            if task = 'FG'. If task = 'CONV' then stp satisfies */
  /*            the sufficient decrease and curvature condition. */

  /*       f is a double precision variable. */
  /*         On initial entry f is the value of the function at 0. */
  /*            On subsequent entries f is the value of the */
  /*            function at stp. */
  /*         On exit f is the value of the function at stp. */

  /*       g is a double precision variable. */
  /*         On initial entry g is the derivative of the function at 0. */
  /*            On subsequent entries g is the derivative of the */
  /*            function at stp. */
  /*         On exit g is the derivative of the function at stp. */

  /*       ftol is a double precision variable. */
  /*         On entry ftol specifies a nonnegative tolerance for the */
  /*            sufficient decrease condition. */
  /*         On exit ftol is unchanged. */

  /*       gtol is a double precision variable. */
  /*         On entry gtol specifies a nonnegative tolerance for the */
  /*            curvature condition. */
  /*         On exit gtol is unchanged. */

  /*       xtol is a double precision variable. */
  /*         On entry xtol specifies a nonnegative relative tolerance */
  /*            for an acceptable step. The subroutine exits with a */
  /*            warning if the relative difference between sty and stx */
  /*            is less than xtol. */
  /*         On exit xtol is unchanged. */

  /*       task is a character variable of length at least 60. */
  /*         On initial entry task must be set to 'START'. */
  /*         On exit task indicates the required action: */

  /*            If task(1:2) = 'FG' then evaluate the function and */
  /*            derivative at stp and call dcsrch again. */

  /*            If task(1:4) = 'CONV' then the search is successful. */

  /*            If task(1:4) = 'WARN' then the subroutine is not able */
  /*            to satisfy the convergence conditions. The exit value of */
  /*            stp contains the best point found during the search. */

  /*            If task(1:5) = 'ERROR' then there is an error in the */
  /*            input arguments. */

  /*         On exit with convergence, a warning or an error, the */
  /*            variable task contains additional information. */

  /*       stpmin is a double precision variable. */
  /*         On entry stpmin is a nonnegative lower bound for the step. */
  /*         On exit stpmin is unchanged. */

  /*       stpmax is a double precision variable. */
  /*         On entry stpmax is a nonnegative upper bound for the step. */
  /*         On exit stpmax is unchanged. */

  /*       isave is an integer work array of dimension 2. */

  /*       dsave is a double precision work array of dimension 13. */

  /*     Subprograms called */

  /*       MINPACK-2 ... dcstep */

  /*     MINPACK-1 Project. June 1983. */
  /*     Argonne National Laboratory. */
  /*     Jorge J. More' and David J. Thuente. */

  /*     MINPACK-2 Project. November 1993. */
  /*     Argonne National Laboratory and University of Minnesota. */
  /*     Brett M. Averick, Richard G. Carter, and Jorge J. More'. */

  /*     ********** */
  /*     Initialization block. */
  /* Parameter adjustments */
  --dsave;
  --isave;

  /* Function Body */
  if (s_cmp(task, "START", (ftnlen)5, (ftnlen)5) == 0) {
    /*        Check the input arguments for errors. */
    if (*stp < *stpmin) { s_copy(task, "ERROR: STP .LT. STPMIN", task_len, (ftnlen)22); }
    if (*stp > *stpmax) { s_copy(task, "ERROR: STP .GT. STPMAX", task_len, (ftnlen)22); }
    if (*g >= 0.) { s_copy(task, "ERROR: INITIAL G .GE. ZERO", task_len, (ftnlen)26); }
    if (*ftol < 0.) { s_copy(task, "ERROR: FTOL .LT. ZERO", task_len, (ftnlen)21); }
    if (*gtol < 0.) { s_copy(task, "ERROR: GTOL .LT. ZERO", task_len, (ftnlen)21); }
    if (*xtol < 0.) { s_copy(task, "ERROR: XTOL .LT. ZERO", task_len, (ftnlen)21); }
    if (*stpmin < 0.) { s_copy(task, "ERROR: STPMIN .LT. ZERO", task_len, (ftnlen)23); }
    if (*stpmax < *stpmin) { s_copy(task, "ERROR: STPMAX .LT. STPMIN", task_len, (ftnlen)25); }
    /*        Exit if there are errors on input. */
    if (s_cmp(task, "ERROR", (ftnlen)5, (ftnlen)5) == 0) { return 0; }
    /*        Initialize local variables. */
    brackt = FALSE_;
    stage = 1;
    finit = *f;
    ginit = *g;
    gtest = *ftol * ginit;
    width = *stpmax - *stpmin;
    width1 = width / .5;
    /*        The variables stx, fx, gx contain the values of the step, */
    /*        function, and derivative at the best step. */
    /*        The variables sty, fy, gy contain the value of the step, */
    /*        function, and derivative at sty. */
    /*        The variables stp, f, g contain the values of the step, */
    /*        function, and derivative at stp. */
    stx = 0.;
    fx = finit;
    gx = ginit;
    sty = 0.;
    fy = finit;
    gy = ginit;
    stmin = 0.;
    stmax = *stp + *stp * 4.;
    s_copy(task, "FG", task_len, (ftnlen)2);
    goto L10;
  }
  else {
    /*        Restore local variables. */
    if (isave[1] == 1) { brackt = TRUE_; }
    else {
      brackt = FALSE_;
    }
    stage = isave[2];
    ginit = dsave[1];
    gtest = dsave[2];
    gx = dsave[3];
    gy = dsave[4];
    finit = dsave[5];
    fx = dsave[6];
    fy = dsave[7];
    stx = dsave[8];
    sty = dsave[9];
    stmin = dsave[10];
    stmax = dsave[11];
    width = dsave[12];
    width1 = dsave[13];
  }
  /*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
  /*     algorithm enters the second stage. */
  ftest = finit + *stp * gtest;
  if (stage == 1 && *f <= ftest && *g >= 0.) { stage = 2; }
  /*     Test for warnings. */
  if (brackt && (*stp <= stmin || *stp >= stmax)) {
    s_copy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS", task_len, (ftnlen)41);
  }
  if (brackt && stmax - stmin <= *xtol * stmax) {
    s_copy(task, "WARNING: XTOL TEST SATISFIED", task_len, (ftnlen)28);
  }
  if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
    s_copy(task, "WARNING: STP = STPMAX", task_len, (ftnlen)21);
  }
  if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
    s_copy(task, "WARNING: STP = STPMIN", task_len, (ftnlen)21);
  }
  /*     Test for convergence. */
  if (*f <= ftest && ABS(*g) <= *gtol * (-ginit)) {
    s_copy(task, "CONVERGENCE", task_len, (ftnlen)11);
  }
  /*     Test for termination. */
  if (s_cmp(task, "WARN", (ftnlen)4, (ftnlen)4) == 0 ||
      s_cmp(task, "CONV", (ftnlen)4, (ftnlen)4) == 0) {
    goto L10;
  }
  /*     A modified function is used to predict the step during the */
  /*     first stage if a lower function value has been obtained but */
  /*     the decrease is not sufficient. */
  if (stage == 1 && *f <= fx && *f > ftest) {
    /*        Define the modified function and derivative values. */
    fm = *f - *stp * gtest;
    fxm = fx - stx * gtest;
    fym = fy - sty * gtest;
    gm = *g - gtest;
    gxm = gx - gtest;
    gym = gy - gtest;
    /*        Call dcstep to update stx, sty, and to compute the new step. */
    dcstep(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &stmin, &stmax);
    /*        Reset the function and derivative values for f. */
    fx = fxm + stx * gtest;
    fy = fym + sty * gtest;
    gx = gxm + gtest;
    gy = gym + gtest;
  }
  else {
    /*       Call dcstep to update stx, sty, and to compute the new step. */
    dcstep(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &stmax);
  }
  /*     Decide if a bisection step is needed. */
  if (brackt) {
    if ((d__1 = sty - stx, ABS(d__1)) >= width1 * .66) { *stp = stx + (sty - stx) * .5; }
    width1 = width;
    width = (d__1 = sty - stx, ABS(d__1));
  }
  /*     Set the minimum and maximum steps allowed for stp. */
  if (brackt) {
    stmin = MIN(stx, sty);
    stmax = MAX(stx, sty);
  }
  else {
    stmin = *stp + (*stp - stx) * 1.1;
    stmax = *stp + (*stp - stx) * 4.;
  }
  /*     Force the step to be within the bounds stpmax and stpmin. */
  *stp = MAX(*stp, *stpmin);
  *stp = MIN(*stp, *stpmax);
  /*     If further progress is not possible, let stp be the best */
  /*     point obtained during the search. */
  if ((brackt && (*stp <= stmin || *stp >= stmax)) || (brackt && stmax - stmin <= *xtol * stmax)) {
    *stp = stx;
  }
  /*     Obtain another function and derivative. */
  s_copy(task, "FG", task_len, (ftnlen)2);
L10:
  /*     Save local variables. */
  if (brackt) { isave[1] = 1; }
  else {
    isave[1] = 0;
  }
  isave[2] = stage;
  dsave[1] = ginit;
  dsave[2] = gtest;
  dsave[3] = gx;
  dsave[4] = gy;
  dsave[5] = finit;
  dsave[6] = fx;
  dsave[7] = fy;
  dsave[8] = stx;
  dsave[9] = sty;
  dsave[10] = stmin;
  dsave[11] = stmax;
  dsave[12] = width;
  dsave[13] = width1;
  return 0;
} /* dcsrch */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int dcstep(doublereal *stx, doublereal *fx, doublereal *dx, doublereal *sty,
                            doublereal *fy, doublereal *dy, doublereal *stp, doublereal *fp,
                            doublereal *dp, logical *brackt, doublereal *stpmin,
                            doublereal *stpmax) {
  /* System generated locals */
  doublereal d__1, d__2, d__3;

  /* Local variables */
  doublereal sgnd, stpc, stpf, stpq, p, q, gamma, r__, s, theta;

  /*     ********** */

  /*     Subroutine dcstep */

  /*     This subroutine computes a safeguarded step for a search */
  /*     procedure and updates an interval that contains a step that */
  /*     satisfies a sufficient decrease and a curvature condition. */

  /*     The parameter stx contains the step with the least function */
  /*     value. If brackt is set to .true. then a minimizer has */
  /*     been bracketed in an interval with endpoints stx and sty. */
  /*     The parameter stp contains the current step. */
  /*     The subroutine assumes that if brackt is set to .true. then */

  /*           MIN(stx,sty) < stp < MAX(stx,sty), */

  /*     and that the derivative at stx is negative in the direction */
  /*     of the step. */

  /*     The subroutine statement is */

  /*       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, */
  /*                         stpmin,stpmax) */

  /*     where */

  /*       stx is a double precision variable. */
  /*         On entry stx is the best step obtained so far and is an */
  /*            endpoint of the interval that contains the minimizer. */
  /*         On exit stx is the updated best step. */

  /*       fx is a double precision variable. */
  /*         On entry fx is the function at stx. */
  /*         On exit fx is the function at stx. */

  /*       dx is a double precision variable. */
  /*         On entry dx is the derivative of the function at */
  /*            stx. The derivative must be negative in the direction of */
  /*            the step, that is, dx and stp - stx must have opposite */
  /*            signs. */
  /*         On exit dx is the derivative of the function at stx. */

  /*       sty is a double precision variable. */
  /*         On entry sty is the second endpoint of the interval that */
  /*            contains the minimizer. */
  /*         On exit sty is the updated endpoint of the interval that */
  /*            contains the minimizer. */

  /*       fy is a double precision variable. */
  /*         On entry fy is the function at sty. */
  /*         On exit fy is the function at sty. */

  /*       dy is a double precision variable. */
  /*         On entry dy is the derivative of the function at sty. */
  /*         On exit dy is the derivative of the function at the exit sty. */

  /*       stp is a double precision variable. */
  /*         On entry stp is the current step. If brackt is set to .true. */
  /*            then on input stp must be between stx and sty. */
  /*         On exit stp is a new trial step. */

  /*       fp is a double precision variable. */
  /*         On entry fp is the function at stp */
  /*         On exit fp is unchanged. */

  /*       dp is a double precision variable. */
  /*         On entry dp is the the derivative of the function at stp. */
  /*         On exit dp is unchanged. */

  /*       brackt is an logical variable. */
  /*         On entry brackt specifies if a minimizer has been bracketed. */
  /*            Initially brackt must be set to .false. */
  /*         On exit brackt specifies if a minimizer has been bracketed. */
  /*            When a minimizer is bracketed brackt is set to .true. */

  /*       stpmin is a double precision variable. */
  /*         On entry stpmin is a lower bound for the step. */
  /*         On exit stpmin is unchanged. */

  /*       stpmax is a double precision variable. */
  /*         On entry stpmax is an upper bound for the step. */
  /*         On exit stpmax is unchanged. */

  /*     MINPACK-1 Project. June 1983 */
  /*     Argonne National Laboratory. */
  /*     Jorge J. More' and David J. Thuente. */

  /*     MINPACK-2 Project. November 1993. */
  /*     Argonne National Laboratory and University of Minnesota. */
  /*     Brett M. Averick and Jorge J. More'. */

  /*     ********** */
  sgnd = *dp * (*dx / ABS(*dx));
  /*     First case: A higher function value. The minimum is bracketed. */
  /*     If the cubic step is closer to stx than the quadratic step, the */
  /*     cubic step is taken, otherwise the average of the cubic and */
  /*     quadratic steps is taken. */
  if (*fp > *fx) {
    theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
    /* Computing MAX */
    d__1 = ABS(theta), d__2 = ABS(*dx), d__1 = MAX(d__1, d__2), d__2 = ABS(*dp);
    s = MAX(d__1, d__2);
    /* Computing 2nd power */
    d__1 = theta / s;
    gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
    if (*stp < *stx) { gamma = -gamma; }
    p = gamma - *dx + theta;
    q = gamma - *dx + gamma + *dp;
    r__ = p / q;
    stpc = *stx + r__ * (*stp - *stx);
    stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp - *stx);
    if ((d__1 = stpc - *stx, ABS(d__1)) < (d__2 = stpq - *stx, ABS(d__2))) { stpf = stpc; }
    else {
      stpf = stpc + (stpq - stpc) / 2.;
    }
    *brackt = TRUE_;
    /*     Second case: A lower function value and derivatives of opposite */
    /*     sign. The minimum is bracketed. If the cubic step is farther from */
    /*     stp than the secant step, the cubic step is taken, otherwise the */
    /*     secant step is taken. */
  }
  else if (sgnd < 0.) {
    theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
    /* Computing MAX */
    d__1 = ABS(theta), d__2 = ABS(*dx), d__1 = MAX(d__1, d__2), d__2 = ABS(*dp);
    s = MAX(d__1, d__2);
    /* Computing 2nd power */
    d__1 = theta / s;
    gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
    if (*stp > *stx) { gamma = -gamma; }
    p = gamma - *dp + theta;
    q = gamma - *dp + gamma + *dx;
    r__ = p / q;
    stpc = *stp + r__ * (*stx - *stp);
    stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
    if ((d__1 = stpc - *stp, ABS(d__1)) > (d__2 = stpq - *stp, ABS(d__2))) { stpf = stpc; }
    else {
      stpf = stpq;
    }
    *brackt = TRUE_;
    /*     Third case: A lower function value, derivatives of the same sign, */
    /*     and the magnitude of the derivative decreases. */
  }
  else if (ABS(*dp) < ABS(*dx)) {
    /*        The cubic step is computed only if the cubic tends to infinity */
    /*        in the direction of the step or if the minimum of the cubic */
    /*        is beyond stp. Otherwise the cubic step is defined to be the */
    /*        secant step. */
    theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
    /* Computing MAX */
    d__1 = ABS(theta), d__2 = ABS(*dx), d__1 = MAX(d__1, d__2), d__2 = ABS(*dp);
    s = MAX(d__1, d__2);
    /*        The case gamma = 0 only arises if the cubic does not tend */
    /*        to infinity in the direction of the step. */
    /* Computing MAX */
    /* Computing 2nd power */
    d__3 = theta / s;
    d__1 = 0., d__2 = d__3 * d__3 - *dx / s * (*dp / s);
    gamma = s * sqrt((MAX(d__1, d__2)));
    if (*stp > *stx) { gamma = -gamma; }
    p = gamma - *dp + theta;
    q = gamma + (*dx - *dp) + gamma;
    r__ = p / q;
    if (r__ < 0. && gamma != 0.) { stpc = *stp + r__ * (*stx - *stp); }
    else if (*stp > *stx) {
      stpc = *stpmax;
    }
    else {
      stpc = *stpmin;
    }
    stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
    if (*brackt) {
      /*           A minimizer has been bracketed. If the cubic step is */
      /*           closer to stp than the secant step, the cubic step is */
      /*           taken, otherwise the secant step is taken. */
      if ((d__1 = stpc - *stp, ABS(d__1)) < (d__2 = stpq - *stp, ABS(d__2))) { stpf = stpc; }
      else {
        stpf = stpq;
      }
      if (*stp > *stx) {
        /* Computing MIN */
        d__1 = *stp + (*sty - *stp) * .66;
        stpf = MIN(d__1, stpf);
      }
      else {
        /* Computing MAX */
        d__1 = *stp + (*sty - *stp) * .66;
        stpf = MAX(d__1, stpf);
      }
    }
    else {
      /*           A minimizer has not been bracketed. If the cubic step is */
      /*           farther from stp than the secant step, the cubic step is */
      /*           taken, otherwise the secant step is taken. */
      if ((d__1 = stpc - *stp, ABS(d__1)) > (d__2 = stpq - *stp, ABS(d__2))) { stpf = stpc; }
      else {
        stpf = stpq;
      }
      stpf = MIN(*stpmax, stpf);
      stpf = MAX(*stpmin, stpf);
    }
    /*     Fourth case: A lower function value, derivatives of the same sign, */
    /*     and the magnitude of the derivative does not decrease. If the */
    /*     minimum is not bracketed, the step is either stpmin or stpmax, */
    /*     otherwise the cubic step is taken. */
  }
  else {
    if (*brackt) {
      theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
      /* Computing MAX */
      d__1 = ABS(theta), d__2 = ABS(*dy), d__1 = MAX(d__1, d__2), d__2 = ABS(*dp);
      s = MAX(d__1, d__2);
      /* Computing 2nd power */
      d__1 = theta / s;
      gamma = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
      if (*stp > *sty) { gamma = -gamma; }
      p = gamma - *dp + theta;
      q = gamma - *dp + gamma + *dy;
      r__ = p / q;
      stpc = *stp + r__ * (*sty - *stp);
      stpf = stpc;
    }
    else if (*stp > *stx) {
      stpf = *stpmax;
    }
    else {
      stpf = *stpmin;
    }
  }
  /*     Update the interval which contains a minimizer. */
  if (*fp > *fx) {
    *sty = *stp;
    *fy = *fp;
    *dy = *dp;
  }
  else {
    if (sgnd < 0.) {
      *sty = *stx;
      *fy = *fx;
      *dy = *dx;
    }
    *stx = *stp;
    *fx = *fp;
    *dx = *dp;
  }
  /*     Compute the new step. */
  *stp = stpf;
  return 0;
} /* dcstep */
