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
#define EIGS_RESET_CG (1<<5)
#define EIGS_ANALYTIC_LINMIN (1<<6)
#define EIGS_ANALYTIC_MULTIMIN (1<<7)
#define EIGS_DELAY_CG (1<<8)

/* default flags: what we think works best most of the time: */
#define EIGS_DEFAULT_FLAGS (EIGS_DIAGONALIZE_EACH_STEP | EIGS_CONVERGE_EACH_EIGENVALUE | EIGS_RESET_CG | EIGS_ANALYTIC_LINMIN | EIGS_DELAY_CG)

typedef struct evectconstraint_chain_struct {
     evectconstraint C;
     void *constraint_data;
     struct evectconstraint_chain_struct *next;
} evectconstraint_chain;

extern evectconstraint_chain *evect_add_constraint(evectconstraint_chain 
						   *constraints,
						   evectconstraint C,
						   void *constraint_data);
extern void evect_destroy_constraints(evectconstraint_chain *constraints);
extern void evectconstraint_chain_func(evectmatrix X, void *data);

#endif /* EIGENSOLVER_H */
