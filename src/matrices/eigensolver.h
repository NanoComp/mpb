#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <matrices.h>

typedef void (*evectoperator) (evectmatrix Xin, evectmatrix Xout);

extern void eigensolver(evectmatrix Y, real *eigenvals,
			evectoperator A, evectoperator C,
			evectmatrix Work[], int nWork,
			real tolerance, int *num_iterations);

#endif /* EIGENSOLVER_H */
