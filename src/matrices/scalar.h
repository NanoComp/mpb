#ifndef SCALAR_H
#define SCALAR_H

#ifdef SCALAR_SINGLE_PREC
typedef float real;
#define SCALAR_MPI_TYPE MPI_FLOAT
#else
typedef double real;
#define SCALAR_MPI_TYPE MPI_DOUBLE
#endif

/************************** scalars are complex **************************/
#ifdef SCALAR_COMPLEX

typedef struct {
     real re, im;
} scalar;

#define SCALAR_NUMVALS 2

#define ASSIGN_CONJ(a, b) { (a).re = (b).re; (a).im = -(b).im; }

#define SCALAR_INIT_ZERO { 0.0, 0.0 }
#define ASSIGN_ZERO(a) { (a).re = 0.0; (a).im = 0.0; }

#define ACCUMULATE_SUM(sum, a) { (sum).re += (a).re; (sum).im += (a).im; }

#define SCALAR_RE(a) ((a).re)
#define SCALAR_IM(a) ((a).im)

#define ASSIGN_REAL(a, c) { (a).re = (c); (a).im = 0.0; }
#define ASSIGN_SCALAR(a, re, im) { (a).re = (re); (a).im = (im); }

#define SCALAR_NORMSQR(a) ((a).re * (a).re + (a).im * (a).im)

/*************************** scalars are real ****************************/
#else /* scalars are real */

typedef real scalar;

#define SCALAR_NUMVALS 1

#define ASSIGN_CONJ(a, b) (a) = (b);

#define SCALAR_INIT_ZERO 0.0
#define ASSIGN_ZERO(a) (a) = 0.0;

#define ACCUMULATE_SUM(sum, a) (sum) += (a);

#define SCALAR_RE(a) (a)
#define SCALAR_IM(a) 0.0

#define ASSIGN_REAL(a, c) (a) = (c);
#define ASSIGN_SCALAR(a, re, im) (a) = (re);

#define SCALAR_NORMSQR(a) ((a) * (a))

#endif

#endif /* SCALAR_H */
