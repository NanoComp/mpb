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

#ifndef ELASTIC_H
#define ELASTIC_H

#include <scalar.h>
#include <matrices.h>

#if defined(HAVE_LIBFFTW)
#  include <fftw.h>
#  include <rfftw.h>
#  ifdef HAVE_MPI
#    include <fftw_mpi.h>
#    include <rfftw_mpi.h>
#  endif
#elif defined(HAVE_LIBDFFTW)
#  include <dfftw.h>
#  include <drfftw.h>
#  ifdef HAVE_MPI
#    include <dfftw_mpi.h>
#    include <drfftw_mpi.h>
#  endif
#elif defined(HAVE_LIBSFFTW)
#  include <sfftw.h>
#  include <srfftw.h>
#  ifdef HAVE_MPI
#    include <sfftw_mpi.h>
#    include <srfftw_mpi.h>
#  endif
#endif

#if defined(HAVE_LIBFFTW) || defined(HAVE_LIBDFFTW) || defined(HAVE_LIBSFFTW)
#  define HAVE_FFTW 1
#endif

/* Data structure to hold the upper triangle of a symmetric real matrix
   or possibly a Hermitian complex matrix (e.g. the dielectric tensor). */
typedef struct {
#ifdef WITH_HERMITIAN_EPSILON
     real m00, m11, m22;
     scalar_complex m01, m02, m12;
#else
     real m00, m01, m02,
               m11, m12,
                    m22;
#endif
} symmetric_matrix;

#ifdef WITH_HERMITIAN_EPSILON
#  define DIAG_SYMMETRIC_MATRIX(m) ((m).m01.re == 0.0 && (m).m01.im == 0.0 && \
				    (m).m02.re == 0.0 && (m).m02.im == 0.0 && \
				    (m).m12.re == 0.0 && (m).m12.im == 0.0)
#else
#  define DIAG_SYMMETRIC_MATRIX(m) ((m).m01 == 0.0 && \
				    (m).m02 == 0.0 && \
				    (m).m12 == 0.0)
#endif

#define NO_PARITY (0)
#define EVEN_Z_PARITY (1<<0)
#define ODD_Z_PARITY (1<<1)
#define EVEN_Y_PARITY (1<<2)
#define ODD_Y_PARITY (1<<3)

typedef struct {
     int nx, ny, nz;
     int local_nx, local_ny;
     int local_x_start, local_y_start;
     int last_dim, last_dim_size, other_dims;

     int num_bands;
     int N, local_N, N_start, alloc_N;

     int fft_output_size;

     int max_fft_bands, num_fft_bands;

     real G[3][3]; /* G vectors (in cartesian) */
     real current_k[3];  /* (in cartesian basis) */
     int zero_k;
     int parity;

#ifdef HAVE_FFTW
#  ifdef HAVE_MPI
#    ifdef SCALAR_COMPLEX
     fftwnd_mpi_plan plan, iplan;
#    else
     rfftwnd_mpi_plan plan, iplan;
#    endif
#  else
#    ifdef SCALAR_COMPLEX
     fftwnd_plan plan, iplan;
#    else
     rfftwnd_plan plan, iplan;
#    endif
#  endif
#endif

     scalar *fft_data;
     scalar_complex *crap;

     real *rho, *sqrt_rhoinv, *rhoct2, *rhocl2;
} elastic_data;

extern elastic_data *create_elastic_data(int nx, int ny, int nz,
					 int *local_N, int *N_start,
					 int *alloc_N,
					 int num_bands,
					 int num_fft_bands);
extern void destroy_elastic_data(elastic_data *d);

extern void elastic_set_num_bands(elastic_data *d, int num_bands);

extern void update_elastic_data_k(elastic_data *d, real k[3],
				  real G1[3], real G2[3], real G3[3]);

extern void set_elastic_data_parity(elastic_data *d, int parity);

typedef void (*elastic_material_function) (real *rho, real *ct, real *cl,
					   real r[3], void *material_data);

extern void set_elastic_materials(elastic_data *md,
				  const int mesh_size[3],
				  real R[3][3], real G[3][3],
				  elastic_material_function mat,
				  void *material_data);

extern void elastic_sym_matrix_eigs(real eigs[3], const symmetric_matrix *V);
extern void elastic_sym_matrix_invert(symmetric_matrix *Vinv,
                                      const symmetric_matrix *V);

extern void elastic_compute_v_from_V(elastic_data *d, evectmatrix Xin,
				     scalar_complex *dfield,
				     int cur_band_start, int cur_num_bands);
extern void elastic_compute_u_from_v(elastic_data *d,
				     scalar_complex *field,
				     int cur_num_bands);

extern void elastic_vectorfield_otherhalf(elastic_data *d,
					  scalar_complex *field,
					  real phasex,real phasey,real phasez);
extern void elastic_scalarfield_otherhalf(elastic_data *d, real *field);

void assign_symmatrix_vector(scalar_complex *newv,
                             const symmetric_matrix matrix,
                             const scalar_complex *oldv);

extern void elastic_compute_fft(int dir, elastic_data *d, scalar *array,
				int howmany, int stride, int dist);

extern void elastic_operator(evectmatrix Xin, evectmatrix Xout, void *data,
			     int is_current_eigenvector, evectmatrix Work);
extern void elastic_preconditioner(evectmatrix Xin, evectmatrix Xout,
				   void *data,
				   evectmatrix Y, real *eigenvals,
				   sqmatrix YtY);

extern int elastic_zero_k_num_const_bands(evectmatrix X, elastic_data *d);
extern void elastic_zero_k_set_const_bands(evectmatrix X, elastic_data *d);

extern void elastic_parity_constraint(evectmatrix X, void *data);
extern void elastic_zparity_constraint(evectmatrix X, void *data);
extern void elastic_yparity_constraint(evectmatrix X, void *data);
extern real *elastic_zparity(evectmatrix X, elastic_data *d);
extern real *elastic_yparity(evectmatrix X, elastic_data *d);

#endif /* ELASTIC_H */
