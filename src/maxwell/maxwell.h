#ifndef MAXWELL_H
#define MAXWELL_H

#include <scalar.h>
#include <matrices.h>

#ifdef HAVE_FFTW
#  include <fftw.h>
#  include <rfftw.h>
#  ifdef HAVE_MPI
#    include <fftwnd_mpi.h>
#  endif
#endif

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
   (e.g. the dielectric tensor). */
typedef struct {
     real m00, m01, m02,
               m11, m12,
                    m22;
} symmetric_matrix;

typedef struct {
     int nx, ny, nz;
     int local_nx, local_ny;
     int local_x_start, local_y_start;
     int last_dim, other_dims;

     int fft_output_size;

     int num_fft_bands;

#ifdef HAVE_FFTW
#  ifdef HAVE_MPI
     fftwnd_mpi_plan plan, iplan;
#  else
#    ifdef SCALAR_COMPLEX
     fftwnd_plan plan, iplan;
#    else
     rfftwnd_plan plan, iplan;
#    endif
#  endif
#endif

     scalar *fft_data;
     
     k_data *k_plus_G;
     real *k_plus_G_normsqr;

     symmetric_matrix *eps_inv;
     real *eps_inv_mean;
} maxwell_data;

extern maxwell_data *create_maxwell_data(int nx, int ny, int nz,
					 int *local_N, int *N_start,
					 int *alloc_N,
					 int num_bands,
					 int num_fft_bands);
extern void destroy_maxwell_data(maxwell_data *d);

extern void update_maxwell_data_k(maxwell_data *d, real k[3],
				  real G1[3], real G2[3], real G3[3]);

extern void maxwell_operator(evectmatrix Xin, evectmatrix Xout, void *data,
			     int is_current_eigenvector);
extern void maxwell_preconditioner(evectmatrix Xin, evectmatrix Xout,
				   void *data,
				   evectmatrix Y, real *eigenvals);
#endif /* MAXWELL_H */
