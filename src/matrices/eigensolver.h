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
			real tolerance, int *num_iterations,
			int flags);

extern void eigensolver_get_eigenvals(evectmatrix Y, real *eigenvals,
				      evectoperator A, void *Adata,
				      evectmatrix Work1, evectmatrix Work2);

/* eigensolver option flags, designed to be combined with a bitwise or ('|');
   each flag should set exactly one bit. */
#define EIGS_NORMALIZE_FIRST_STEP (1<<0)
#define EIGS_DIAGONALIZE_EACH_STEP (1<<1)
#define EIGS_PROJECT_PRECONDITIONING (1<<2)
#define EIGS_CONVERGE_EACH_EIGENVALUE (1<<3)
#define EIGS_VERBOSE (1<<4)

/* default flags: what we think works best most of the time: */
#define EIGS_DEFAULT_FLAGS (EIGS_DIAGONALIZE_EACH_STEP | EIGS_CONVERGE_EACH_EIGENVALUE)

#endif /* EIGENSOLVER_H */
