/* Copyright (C) 1999 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <matrices.h>

typedef void (*evectoperator) (evectmatrix Xin, evectmatrix Xout,
			       void *data,
			       int is_current_eigenvector);

typedef void (*evectpreconditioner) (evectmatrix Xin, evectmatrix Xout,
				     void *data,
				     evectmatrix Y, real *eigenvals);

typedef void (*evectpreconditioner_data_updater) (void *data,
						  sqmatrix Yrotation);

typedef void (*evectconstraint) (evectmatrix X, void *data);

extern void eigensolver(evectmatrix Y, real *eigenvals,
			evectoperator A, void *Adata,
			evectpreconditioner C, void *Cdata,
			evectpreconditioner_data_updater Cdata_update,
			evectconstraint constraint, void *constraint_data,
			evectmatrix Work[], int nWork,
			real tolerance, int *num_iterations);

extern void eigensolver_get_eigenvals(evectmatrix Y, real *eigenvals,
				      evectoperator A, void *Adata,
				      evectmatrix Work1, evectmatrix Work2);

#endif /* EIGENSOLVER_H */
