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

#ifndef MAXWELL_H
#define MAXWELL_H

#include "scalar.h"
#include "matrices.h"

/* This data structure is designed to hold k+G related data at a given
   point.  kmag is the length of the k+G vector.  The m and n vectors are
   orthonormal vectors orthogonal to (kx,ky,kz).  These are used
   as the basis for the H vector (to maintain transversality). */
typedef struct {
     real kmag;
     real mx, my, mz;
     real nx, ny, nz;
} k_data;


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

     real current_k[3];  /* (in cartesian basis) */
     int parity;

     void *plan, *iplan;
     scalar *fft_data;
     
     int zero_k;  /* non-zero if k is zero (handled specially) */
     k_data *k_plus_G;
     real *k_plus_G_normsqr;

     symmetric_matrix *eps_inv;
     real eps_inv_mean;
} maxwell_data;

extern maxwell_data *create_maxwell_data(int nx, int ny, int nz,
					 int *local_N, int *N_start,
					 int *alloc_N,
					 int num_bands,
					 int num_fft_bands);
extern void destroy_maxwell_data(maxwell_data *d);

extern void maxwell_set_num_bands(maxwell_data *d, int num_bands);

extern void update_maxwell_data_k(maxwell_data *d, real k[3],
				  real G1[3], real G2[3], real G3[3]);

extern void set_maxwell_data_parity(maxwell_data *d, int parity);

typedef void (*maxwell_dielectric_function) (symmetric_matrix *eps,
					     symmetric_matrix *eps_inv,
					     real r[3], void *epsilon_data);

extern void set_maxwell_dielectric(maxwell_data *md,
				   const int mesh_size[3],
				   real R[3][3], real G[3][3],
				   maxwell_dielectric_function epsilon,
				   void *epsilon_data);

extern void maxwell_sym_matrix_eigs(real eigs[3], const symmetric_matrix *V);
extern void maxwell_sym_matrix_invert(symmetric_matrix *Vinv,
                                      const symmetric_matrix *V);

extern void maxwell_compute_fft(int dir, maxwell_data *d, scalar *array,
				int howmany, int stride, int dist);
extern void maxwell_compute_d_from_H(maxwell_data *d, evectmatrix Xin,
				     scalar_complex *dfield,
				     int cur_band_start, int cur_num_bands);
extern void maxwell_compute_h_from_H(maxwell_data *d, evectmatrix Hin,
				     scalar_complex *hfield,
				     int cur_band_start, int cur_num_bands);
extern void maxwell_compute_e_from_d(maxwell_data *d,
				     scalar_complex *dfield,
				     int cur_num_bands);

extern void maxwell_vectorfield_otherhalf(maxwell_data *d,
					  scalar_complex *field,
					  real phasex,real phasey,real phasez);
extern void maxwell_cscalarfield_otherhalf(maxwell_data *d, 
					   scalar_complex *field,
					   real phasex, real phasey, 
					   real phasez);
extern void maxwell_scalarfield_otherhalf(maxwell_data *d, real *field);

void assign_symmatrix_vector(scalar_complex *newv,
                             const symmetric_matrix matrix,
                             const scalar_complex *oldv);

extern void maxwell_operator(evectmatrix Xin, evectmatrix Xout, void *data,
			     int is_current_eigenvector, evectmatrix Work);
extern void maxwell_simple_precondition(evectmatrix X,
					void *data, real *eigenvals);
extern void maxwell_preconditioner(evectmatrix Xin, evectmatrix Xout,
				   void *data,
				   evectmatrix Y, real *eigenvals,
				   sqmatrix YtY);
extern void maxwell_preconditioner2(evectmatrix Xin, evectmatrix Xout,
				    void *data,
				    evectmatrix Y, real *eigenvals,
				    sqmatrix YtY);

extern void maxwell_ucross_op(evectmatrix Xin, evectmatrix Xout,
			      maxwell_data *d, const real u[3]);

extern void maxwell_parity_constraint(evectmatrix X, void *data);
extern void maxwell_zparity_constraint(evectmatrix X, void *data);
extern void maxwell_yparity_constraint(evectmatrix X, void *data);

extern int maxwell_zero_k_num_const_bands(evectmatrix X, maxwell_data *d);
extern void maxwell_zero_k_set_const_bands(evectmatrix X, maxwell_data *d);
extern void maxwell_zero_k_constraint(evectmatrix X, void *data);

extern real *maxwell_zparity(evectmatrix X, maxwell_data *d);
extern real *maxwell_yparity(evectmatrix X, maxwell_data *d);

typedef struct {
     maxwell_data *d;
     real target_frequency;
} maxwell_target_data;

extern maxwell_target_data *create_maxwell_target_data(maxwell_data *d,
						       real target_frequency);
extern void destroy_maxwell_target_data(maxwell_target_data *d);
extern void maxwell_target_operator1(evectmatrix Xin, evectmatrix Xout,
				     void *data,
				     int is_current_eigenvector,
				     evectmatrix Work);
extern void maxwell_target_operator(evectmatrix Xin, evectmatrix Xout,
				    void *data, int is_current_eigenvector,
				    evectmatrix Work);
extern void maxwell_target_preconditioner(evectmatrix Xin, evectmatrix Xout,
					  void *data,
					  evectmatrix Y, real *eigenvals,
					  sqmatrix YtY);
extern void maxwell_target_preconditioner2(evectmatrix Xin, evectmatrix Xout,
					   void *data,
					   evectmatrix Y, real *eigenvals,
					   sqmatrix YtY);

extern void spherical_quadrature_points(real *x, real *y, real *z,
					real *weight, int num_sq_pts);

extern int check_maxwell_dielectric(maxwell_data *d,
				    int negative_epsilon_okp);

#endif /* MAXWELL_H */
