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

#ifndef LINMIN_H
#define LINMIN_H

typedef double (*linmin_func) (double x, double *deriv, void *data);

extern double linmin(double *converged_f, double *converged_df,
		     double x_guess, double f_0, double df_0,
		     double f_tol, double df_tol, double x_tol,
		     double x_min, double x_max,
		     linmin_func f, void *f_data, int verbose);

#endif /* LINMIN_H */
