/* Copyright (C) 1999-2014 Massachusetts Institute of Technology.
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

/* optimization routines for material grids */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>

#include "mpb.h"

#ifdef HAVE_NLOPT_H
#  include <nlopt.h>
#endif

/**************************************************************************/
/* optimization of band gaps as a function of the material grid */

typedef struct {
     boolean do_min;
     vector3_list ks;
     int b1, b2;
     int ngrids;
     material_grid *grids;
     int iter, unsolved;
     double *f1s, *f2s; /* arrays of length ks.num_items for freqs */
     double *work; /* work array of length ntot */
} maxgap_func_data;

/* the constraint is either an upper bound for band b1
   or a lower bound for band b2 */
typedef enum { BAND1_CONSTRAINT = 0, BAND2_CONSTRAINT = 1 } band_constraint_kind;

typedef struct {
     maxgap_func_data *d;
     int ik; /* index of k point for this constraint (in d->ks) */
     band_constraint_kind kind;
} band_constraint_data;

static double band_constraint(int n, const double *u, double *grad, void *data)
{
     band_constraint_data *cdata = (band_constraint_data *) data;
     maxgap_func_data *d = cdata->d;
     int ik = cdata->ik;
     int kind = cdata->kind;
     double *work = d->work;
     double val = 0;

     /* Strictly speaking, we should call material_grids_set here.  However
	we rely on an implementation detail of our MMA code: it always
	evaluates the objective function before evaluating the constraints,
	and hence we can set the material_grids once in the objective. */

     /* We will typically have more than one band per k-point
	(typically 2 bands), and we don't need to call solve_kpoint
	more than once per band.  Here we exploit the fact that our
	MMA code always calls all the constraints at once (in
	sequence); it never changes u in between one constraint & the next. */
     if (!vector3_equal(cur_kvector, d->ks.items[ik]) || d->unsolved) {
	  randomize_fields();
	  solve_kpoint(d->ks.items[ik]);
     }
     d->unsolved = 0;

     if (grad) memset(work, 0, sizeof(double) * (n-2));
     if (kind == BAND1_CONSTRAINT) {
	  if (grad) {
	       material_grids_addgradient(work, 1.0, d->b1, 
					  d->grids, d->ngrids);
	       grad[n-1] = -1;
	       grad[n-2] = 0;
	  }
	  val = (d->f1s[ik] = freqs.items[d->b1-1]) - u[n-1];
     }
     else {
	  if (grad) {
	       material_grids_addgradient(work, -1.0, d->b2, 
					  d->grids, d->ngrids);
	       grad[n-1] = 0;
	       grad[n-2] = 1;
	  }
	  val = u[n-2] - (d->f2s[ik] = freqs.items[d->b2-1]);
     }
     if (grad) /* gradient w.r.t. epsilon needs to be summed over processes */
	  mpi_allreduce(work, grad, n-2, double, MPI_DOUBLE, 
			MPI_SUM, MPI_COMM_WORLD);

     return val;
}

static double maxgap_func(int n, const double *u, double *grad, void *data)
{
     maxgap_func_data *d = (maxgap_func_data *) data;
     double gap, f1 = u[n-1], f2 = u[n-2];

     /* set the material grids, for use in the constraint functions
	and also for outputting in verbose mode */
     material_grids_set(u, d->grids, d->ngrids);
     reset_epsilon();
     d->iter++;
     d->unsolved = 1;

     gap = (f2 - f1) * 2.0 / (f1 + f2);
     
     if (grad) {
	  memset(grad, 0, sizeof(double) * (n-2));
	  grad[n-1] = 2.0 * ((f1 + f2) - (f1 - f2)) / ((f1+f2)*(f1+f2));
	  grad[n-2] = 2.0 * (-(f1 + f2) - (f1 - f2)) / ((f1+f2)*(f1+f2));
	  if (d->do_min) {
	       grad[n-1] = -grad[n-1];
	       grad[n-2] = -grad[n-2];
	  }
     }

     mpi_one_printf("material-grid-%sgap:, %d, %g, %g, %0.15g\n", 
		    d->do_min ? "min" : "max", d->iter, f1, f2, gap);
     
     if (verbose) {
	  char prefix[256];
	  get_epsilon();
	  snprintf(prefix, 256, "%sgap-%04d-", 
		   d->do_min ? "min" : "max", d->iter);
	  output_field_to_file(-1, prefix);
     }

     return d->do_min ? gap : -gap;
}

static number material_grids_maxmin_gap(boolean do_min,
					vector3_list kpoints, 
					integer band1, integer band2,
					number func_tol, number eps_tol,
					integer maxeval, number maxtime)
{
     maxgap_func_data d;
     int i, n;
     double *u, *lb, *ub, *u_tol, func_min;
     band_constraint_data *cdata;
     int have_uprod;

     CHECK(band1>0 && band1 <= num_bands && band2>0 && band2 <= num_bands,
	   "invalid band numbers in material-grid-maxgap");
     d.ks = kpoints;
     d.b1 = band1; d.b2 = band2;
     d.grids = get_material_grids(geometry, &d.ngrids);
     d.iter = 0;
     d.unsolved = 1;
     d.do_min = do_min;
     d.f1s = (double *) malloc(sizeof(double) * kpoints.num_items*2);
     d.f2s = d.f1s + kpoints.num_items;

     n = material_grids_ntot(d.grids, d.ngrids) + 2;
     u = (double *) malloc(sizeof(double) * n * 5);
     lb = u + n; ub = lb + n; u_tol = ub + n; d.work = u_tol + n;

     material_grids_get(u, d.grids, d.ngrids);
     u[n-1] = 0; /* band1 max */
     u[n-2] = HUGE_VAL; /* band2 min */

     cdata = (band_constraint_data*) malloc(sizeof(band_constraint_data)
					    * kpoints.num_items*2);
     for (i = 0; i < kpoints.num_items; ++i) {
	  band_constraint_kind kind;
	  for (kind = BAND1_CONSTRAINT; kind <= BAND2_CONSTRAINT; ++kind) {
	       cdata[2*i + kind].d = &d;
	       cdata[2*i + kind].ik = i;
	       cdata[2*i + kind].kind = kind;

	       /* compute initial band min/max */
	       band_constraint(n, u, NULL, &cdata[2*i + kind]);
	       if (kind == BAND1_CONSTRAINT && d.f1s[i] > u[n-1])
		    u[n-1] = d.f1s[i];
	       else if (kind == BAND2_CONSTRAINT && d.f2s[i] < u[n-2])
		    u[n-2] = d.f2s[i];
	  }
     }
     u[n-1] *= 1.001; u[n-2] /= 1.001; /* ensure feasibility of initial u */

     for (i = 0; i < d.ngrids && d.grids[i].material_grid_kind != U_PROD; ++i);
     have_uprod = i < d.ngrids;
     for (i = 0; i < n-2; ++i) {
	  ub[i] = 1;
	  u_tol[i] = eps_tol;
	  /* bound u slightly about 0 for uprod grids, as when u=0
	     the gradient is problematic (especially for multiple u's = 0 */
	  lb[i] = have_uprod ? 1e-4 : 0;
	  if (u[i] < lb[i]) u[i] = lb[i];
     }
     u_tol[n-1] = u_tol[n-2] = 0;
     lb[n-1] = lb[n-2] = 0;
     ub[n-1] = ub[n-2] = HUGE_VAL;

#if defined(HAVE_NLOPT_H) && defined(HAVE_NLOPT)
 {
     nlopt_result res;
     extern int mma_verbose;
     mma_verbose = kpoints.num_items*2;
     res = nlopt_minimize_constrained(
	  NLOPT_LD_MMA, n, maxgap_func, &d,
	  kpoints.num_items*2, band_constraint, 
	  cdata, sizeof(band_constraint_data),
	  lb, ub, u, &func_min,
	  -HUGE_VAL, func_tol,0, 0,u_tol, maxeval,maxtime);
     CHECK(res > 0, "failure of nlopt_minimize");
 }
#else
     CHECK(0, "nlopt library is required for material-grid-maxgap");
#endif

     maxgap_func(n, u, NULL, &d);

     /* recompute bands and get actual gap size */

     u[n-1] = 0; /* band1 max */
     u[n-2] = HUGE_VAL; /* band2 min */
     for (i = 0; i < kpoints.num_items; ++i) {
	  band_constraint_kind kind;
	  for (kind = BAND1_CONSTRAINT; kind <= BAND2_CONSTRAINT; ++kind) {
	       band_constraint(n, u, NULL, &cdata[2*i + kind]);
	       if (kind == BAND1_CONSTRAINT && d.f1s[i] > u[n-1])
		    u[n-1] = d.f1s[i];
	       else if (kind == BAND2_CONSTRAINT && d.f2s[i] < u[n-2])
		    u[n-2] = d.f2s[i];
	  }
     }

     func_min = (u[n-2] - u[n-1]) * 2.0 / (u[n-1] + u[n-2]);
     mpi_one_printf("material-grid-%sgap:, %d, %g, %g, %0.15g\n",
                    d.do_min ? "min" : "max", d.iter+1, 
		    u[n-1], u[n-2], func_min);
     func_min = d.do_min ? func_min : -func_min;

     free(cdata);
     free(u);
     free(d.grids);
     free(d.f1s);

     return(do_min ? func_min : -func_min);
}

number material_grids_maxgap(vector3_list kpoints, 
			     integer band1, integer band2,
			     number func_tol, number eps_tol,
			     integer maxeval, number maxtime)
{
     return material_grids_maxmin_gap(0, kpoints, band1, band2,
				      func_tol, eps_tol, maxeval, maxtime);
}

number material_grids_mingap(vector3_list kpoints, 
			     integer band1, integer band2,
			     number func_tol, number eps_tol,
			     integer maxeval, number maxtime)
{
     return material_grids_maxmin_gap(1, kpoints, band1, band2,
				      func_tol, eps_tol, maxeval, maxtime);
}

/**************************************************************************/
