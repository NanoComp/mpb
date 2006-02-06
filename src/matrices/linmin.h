/* Copyright (C) 1999, 2000, 2001, 2002, Massachusetts Institute of Technology.
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

#if defined(SCALAR_LONG_DOUBLE_PREC)
typedef long double linmin_real;
#else
typedef double linmin_real;
#endif

typedef linmin_real (*linmin_func) (linmin_real x, linmin_real *deriv, void *data);

extern linmin_real linmin(linmin_real *converged_f, linmin_real *converged_df,
		     linmin_real x_guess, linmin_real f_0, linmin_real df_0,
		     linmin_real f_tol, linmin_real df_tol, linmin_real x_tol,
		     linmin_real x_min, linmin_real x_max,
		     linmin_func f, void *f_data, int verbose);

#endif /* LINMIN_H */
