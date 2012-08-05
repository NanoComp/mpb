/* Copyright (C) 1999-2012, Massachusetts Institute of Technology.
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

/* Nice C wrappers around minpack2-linmin.c functions. */

#include <stdlib.h>
#include <stdio.h>

#include "config.h"
#include <mpi_utils.h>
#include <check.h>

#include "linmin.h"

#define double linmin_real

extern int dcsrch(double *stp, double *f, double *g,
		   double *ftol, double *gtol, double *xtol,
		   char *task,
		   double *stpmin, double *stpmax, int *isave,
		   double *dsave);

/* Our line minimization algorithm is designed only for f(x) for x >= 0.
   If the minimum lies for negative x, we use reverse_func below to
   transform f(x) into f(-x). */
typedef struct {
     linmin_func f;
     void *f_data;
} linmin_func_data;
double reverse_func(double x, double *deriv, void *data)
{
     linmin_func_data *d = (linmin_func_data *) data;
     double val;
     val = d->f(-x, deriv, d->f_data);
     *deriv = -*deriv;
     return val;
}

double linmin(double *converged_f, double *converged_df,
	      double x_guess, double f_0, double df_0,
	      double f_tol, double df_tol, double x_tol,
	      double x_min, double x_max,
	      linmin_func f, void *f_data, int verbose)
{
     if (df_0 > 0) {  /* minimum lies for negative x; transform to f(-x) */
	  linmin_func_data d;
	  double x;
	  d.f = f;
	  d.f_data = f_data;
	  x = linmin(converged_f, converged_df, 
		     -x_guess, f_0, -df_0, f_tol, df_tol, x_tol,
		     -x_min, -x_max, reverse_func, &d, verbose);
	  *converged_df = -*converged_df;
	  return(-x);
     }
     else if (df_0 == 0) { /* already at minimum! */
	  *converged_f = f_0;
          *converged_df = df_0;
          return 0;
     }
     else {
	  char task[300] = "START";
	  int isave[2];
	  double dsave[13], x, f_x, df_x;
	  int iters = 0;

	  x = x_guess;
	  f_x = f_0; df_x = df_0; /* initially, pass in f and df at x=0 */
	  dcsrch(&x, &f_x, &df_x, &f_tol, &df_tol, &x_tol,
		 task, &x_min, &x_max, isave, dsave);
	  
	  while (*task == 'F') {
	       f_x = f(x, &df_x, f_data);
	       mpi_assert_equal(x);
	       mpi_assert_equal(f_x);
	       ++iters;
	       dcsrch(&x, &f_x, &df_x, &f_tol, &df_tol, &x_tol,
		      task, &x_min, &x_max, isave, dsave);
	  }

	  if (*task != 'C') {  /* not converged; warning or error */
	       if (verbose || *task == 'E')
		    mpi_one_fprintf(stderr, "linmin: %s\n", task);
	       CHECK(*task != 'E', "linmin failure");
	  }

	  if (verbose)
	       mpi_one_printf("    linmin: converged after %d iterations.\n",
			      iters);

	  *converged_f = f_x;
	  *converged_df = df_x;
	  return x;
     }
}
