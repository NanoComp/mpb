#ifndef SCALAR_H
#define SCALAR_H

#ifdef SCALAR_SINGLE_PREC
typedef float real;
#define SCALAR_MPI_TYPE MPI_FLOAT
#else
typedef double real;
#define SCALAR_MPI_TYPE MPI_DOUBLE
#endif

typedef struct {
     real re, im;
} scalar_complex;

/************************** scalars are complex **************************/
#ifdef SCALAR_COMPLEX

typedef scalar_complex scalar;

#define SCALAR_NUMVALS 2

#define ASSIGN_CONJ(a, b) { (a).re = (b).re; (a).im = -(b).im; }

#define SCALAR_INIT_ZERO { 0.0, 0.0 }
#define ASSIGN_ZERO(a) { (a).re = 0.0; (a).im = 0.0; }

#define ACCUMULATE_SUM(sum, a) { (sum).re += (a).re; (sum).im += (a).im; }

#define SCALAR_RE(a) ((a).re)
#define SCALAR_IM(a) ((a).im)

#define ASSIGN_REAL(a, c) { (a).re = (c); (a).im = 0.0; }
#define ASSIGN_SCALAR(a, real, imag) { (a).re = (real); (a).im = (imag); }

#define SCALAR_NORMSQR(a) ((a).re * (a).re + (a).im * (a).im)

/* a = b * c */
#define ASSIGN_MULT(a, b, c) { \
     real bbbb_re = (b).re, bbbb_im = (b).im; \
     real cccc_re = (c).re, cccc_im = (c).im; \
     ASSIGN_SCALAR(a, bbbb_re * cccc_re - bbbb_im * cccc_im, \
		   bbbb_re * cccc_im + bbbb_im * cccc_re); \
}

/* a = b / c = b * conj(c) / |c|^2 */
#define ASSIGN_DIV(a, b, c) { \
     scalar aaaa_tmp; real aaaa_tmp_norm; \
     ASSIGN_CONJ(aaaa_tmp, c); \
     aaaa_tmp_norm = 1.0 / SCALAR_NORMSQR(aaaa_tmp); \
     ASSIGN_MULT(aaaa_tmp, b, aaaa_tmp); \
     ASSIGN_SCALAR(a, aaaa_tmp.re*aaaa_tmp_norm, aaaa_tmp.im*aaaa_tmp_norm); \
}

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
#define ASSIGN_SCALAR(a, real, imag) (a) = (real);

#define SCALAR_NORMSQR(a) ((a) * (a))

#define ASSIGN_MULT(a, b, c) (a) = (b) * (c);
#define ASSIGN_DIV(a, b, c) (a) = (b) / (c);

#endif

#endif /* SCALAR_H */
