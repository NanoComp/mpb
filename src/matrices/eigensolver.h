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
