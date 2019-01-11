/* Copyright (C) 1999-2019 Massachusetts Institute of Technology.
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

#ifndef SCALAR_H
#define SCALAR_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if defined(SCALAR_SINGLE_PREC)
typedef float real;
#define SCALAR_MPI_TYPE MPI_FLOAT
#elif defined(SCALAR_LONG_DOUBLE_PREC)
typedef long double real;
#define SCALAR_MPI_TYPE MPI_LONG_DOUBLE
#else
typedef double real;
#define SCALAR_MPI_TYPE MPI_DOUBLE
#endif

/********************* complex types and operations **********************/
typedef struct {
     real re, im;
} scalar_complex;

#define CSCALAR_NUMVALS (2)

#define CSCALAR_INIT_ZERO { 0.0, 0.0 }

#define CSCALAR_RE(a) ((a).re)
#define CSCALAR_IM(a) ((a).im)

#define CASSIGN_SCALAR(a, real, imag) { (a).re = (real); (a).im = (imag); }
#define CACCUMULATE_SCALAR(a, real, imag) { (a).re +=(real); (a).im +=(imag); }
#define CACCUMULATE_DIFF_SCALAR(a, real, imag) { (a).re -=(real); (a).im -=(imag); }

#define CASSIGN_ZERO(a) CASSIGN_SCALAR(a, 0.0, 0.0);
#define CASSIGN_REAL(a, c) CASSIGN_SCALAR(a, c, 0.0)
#define CASSIGN_CONJ(a, b) CASSIGN_SCALAR(a, CSCALAR_RE(b), -CSCALAR_IM(b))

#define CSCALAR_NORMSQR(a) ((a).re * (a).re + (a).im * (a).im)

/* a = b * c */
#define CASSIGN_MULT(a, b, c) { \
     real bbbb_re = (b).re, bbbb_im = (b).im; \
     real cccc_re = (c).re, cccc_im = (c).im; \
     CASSIGN_SCALAR(a, bbbb_re * cccc_re - bbbb_im * cccc_im, \
		    bbbb_re * cccc_im + bbbb_im * cccc_re); \
}

/* a = b / c = b * conj(c) / |c|^2 */
#define CASSIGN_DIV(a, b, c) { \
     scalar_complex aaaa_tmp; real aaaa_tmp_norm; \
     CASSIGN_CONJ(aaaa_tmp, c); \
     aaaa_tmp_norm = 1.0 / CSCALAR_NORMSQR(aaaa_tmp); \
     CASSIGN_MULT(aaaa_tmp, b, aaaa_tmp); \
     CASSIGN_SCALAR(a, aaaa_tmp.re*aaaa_tmp_norm, aaaa_tmp.im*aaaa_tmp_norm); \
}

/* Re (b * c) */
#define CSCALAR_MULT_RE(b, c) ((b).re * (c).re - (b).im * (c).im)

/* IM (b * c) */
#define CSCALAR_MULT_IM(b, c) ((b).re * (c).im + (b).im * (c).re)

/* Re (b * conj(c)) */
#define CSCALAR_MULT_CONJ_RE(b, c) ((b).re * (c).re + (b).im * (c).im)

/* Im (b * conj(c)) */
#define CSCALAR_MULT_CONJ_IM(b, c) ((b).im * (c).re - (b).re * (c).im)

/* a = Re (b * c) */
#define CASSIGN_MULT_RE(a, b, c) { \
     real bbbb_re = (b).re, bbbb_im = (b).im; \
     real cccc_re = (c).re, cccc_im = (c).im; \
     (a) = bbbb_re * cccc_re - bbbb_im * cccc_im; \
}

/* a = Im (b * c) */
#define CASSIGN_MULT_IM(a, b, c) { \
     real bbbb_re = (b).re, bbbb_im = (b).im; \
     real cccc_re = (c).re, cccc_im = (c).im; \
     (a) = bbbb_re * cccc_im + bbbb_im * cccc_re; \
}

/* a += b * c */
#define CACCUMULATE_SUM_MULT(a, b, c) { \
     real bbbb_re = (b).re, bbbb_im = (b).im; \
     real cccc_re = (c).re, cccc_im = (c).im; \
     CACCUMULATE_SCALAR(a, bbbb_re * cccc_re - bbbb_im * cccc_im, \
		       bbbb_re * cccc_im + bbbb_im * cccc_re); \
}

/* a += conj(b) * c */
#define CACCUMULATE_SUM_CONJ_MULT(a, b, c) { \
     real bbbb_re = (b).re, bbbb_im = (b).im; \
     real cccc_re = (c).re, cccc_im = (c).im; \
     CACCUMULATE_SCALAR(a, bbbb_re * cccc_re + bbbb_im * cccc_im, \
		       bbbb_re * cccc_im - bbbb_im * cccc_re); \
}

/* a += |b|^2 */
#define CACCUMULATE_SUM_SQ(a, b) { \
     real bbbb_re = (b).re, bbbb_im = (b).im; \
     (a) += bbbb_re * bbbb_re + bbbb_im * bbbb_im; \
}

#define CACCUMULATE_SUM(sum, a) CACCUMULATE_SCALAR(sum,CSCALAR_RE(a),CSCALAR_IM(a))
#define CACCUMULATE_DIFF(sum, a) CACCUMULATE_DIFF_SCALAR(sum,CSCALAR_RE(a),CSCALAR_IM(a))

/************************** scalars are complex **************************/
#ifdef SCALAR_COMPLEX

typedef scalar_complex scalar;

#define SCALAR_NUMVALS CSCALAR_NUMVALS

#define SCALAR_INIT_ZERO CSCALAR_INIT_ZERO

#define SCALAR_RE(a) CSCALAR_RE(a)
#define SCALAR_IM(a) CSCALAR_IM(a)

#define ASSIGN_SCALAR(a, real, imag) CASSIGN_SCALAR(a, real, imag)
#define ACCUMULATE_SCALAR(a, real, imag) CACCUMULATE_SCALAR(a, real, imag)
#define ACCUMULATE_DIFF_SCALAR(a, real, imag) CACCUMULATE_DIFF_SCALAR(a, real, imag)

#define SCALAR_NORMSQR(a) CSCALAR_NORMSQR(a)

/* a = b * c */
#define ASSIGN_MULT(a, b, c) CASSIGN_MULT(a, b, c)

/* a = b / c = b * conj(c) / |c|^2 */
#define ASSIGN_DIV(a, b, c) CASSIGN_DIV(a, b, c)

#define SCALAR_MULT_RE(b, c) CSCALAR_MULT_RE(b, c) /* Re (b * c) */
#define SCALAR_MULT_IM(b, c) CSCALAR_MULT_IM(b, c) /* Im (b * c) */

#define SCALAR_MULT_CONJ_RE(b, c) CSCALAR_MULT_CONJ_RE(b, c)
#define SCALAR_MULT_CONJ_IM(b, c) CSCALAR_MULT_CONJ_IM(b, c)

/* a = Re (b * c) */
#define ASSIGN_MULT_RE(a, b, c) CASSIGN_MULT_RE(a, b, c)

/* a = Im (b * c) */
#define ASSIGN_MULT_IM(a, b, c) CASSIGN_MULT_IM(a, b, c)

/* a += b * c */
#define ACCUMULATE_SUM_MULT(a, b, c) CACCUMULATE_SUM_MULT(a, b, c)

/* a += conj(b) * c */
#define ACCUMULATE_SUM_CONJ_MULT(a, b, c) CACCUMULATE_SUM_CONJ_MULT(a, b, c)

/* a += |b|^2 */
#define ACCUMULATE_SUM_SQ(a, b) CACCUMULATE_SUM_SQ(a, b)

/*************************** scalars are real ****************************/
#else /* scalars are real */

typedef real scalar;

#define SCALAR_NUMVALS 1

#define SCALAR_INIT_ZERO 0.0

#define SCALAR_RE(a) (a)
#define SCALAR_IM(a) (0.0)

#define ASSIGN_SCALAR(a, real, imag) (a) = (real);
#define ACCUMULATE_SCALAR(a, real, imag) (a) += (real);
#define ACCUMULATE_DIFF_SCALAR(a, real, imag) (a) -= (real);

#define SCALAR_NORMSQR(a) ((a) * (a))

#define ASSIGN_MULT(a, b, c) (a) = (b) * (c);
#define ASSIGN_DIV(a, b, c) (a) = (b) / (c);

#define SCALAR_MULT_RE(b, c) ((b) * (c))
#define SCALAR_MULT_IM(b, c) (0.0)

#define SCALAR_MULT_CONJ_RE(b, c) ((b) * (c))
#define SCALAR_MULT_CONJ_IM(b, c) (0.0)

#define ASSIGN_MULT_RE(a, b, c) (a) = (b) * (c);
#define ASSIGN_MULT_IM(a, b, c) (a) = 0.0;

#define ACCUMULATE_SUM_MULT(a, b, c) (a) += (b) * (c);
#define ACCUMULATE_SUM_CONJ_MULT(a, b, c) (a) += (b) * (c);
#define ACCUMULATE_SUM_SQ(a, b) { real bbbb = (b); (a) += bbbb * bbbb; }

#endif /* scalars are real */

#define ASSIGN_ZERO(a) ASSIGN_SCALAR(a, 0.0, 0.0);
#define ASSIGN_REAL(a, c) ASSIGN_SCALAR(a, c, 0.0)
#define ASSIGN_CONJ(a, b) ASSIGN_SCALAR(a, SCALAR_RE(b), -SCALAR_IM(b))
#define ACCUMULATE_SUM(sum, a) ACCUMULATE_SCALAR(sum,SCALAR_RE(a),SCALAR_IM(a))
#define ACCUMULATE_DIFF(sum, a) ACCUMULATE_DIFF_SCALAR(sum,SCALAR_RE(a),SCALAR_IM(a))

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* SCALAR_H */
