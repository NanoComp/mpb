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

typedef struct {
     int nx, ny, nz;
     int local_nx, local_ny;
     int local_x_start, local_y_start;

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
     
     real *k_plus_G;

     real *eps_inv;

     real *k_plus_G_normsqr_inv;
     real *eps_inv_mean;
} maxwell_data;

#define EPS_MATRIX_COMPONENTS 6  /* see below */

/* eps_inv is an array of 3x3 real matrices.  However, since these matrices
   are guaranteed to be symmetric, we need to store only 6 components,
   whose offsets in the array are indicated as follows:

                          [ 0 1 2 ]
                          [ x 3 4 ]
                          [ x x 5 ]

   An "x" here indicates a matrix element that is not stored. */

extern maxwell_data *create_maxwell_data(int nx, int ny, int nz,
					 int *local_N, int *alloc_N,
					 epsilon_function eps,
					 int num_bands,
					 int num_fft_bands,
					 real k[3]);

extern void update_maxwell_data_k(maxwell_data *d, real k[3],
				  real G1[3], real G2[3], real G3[3]);


#endif /* MAXWELL_H */
