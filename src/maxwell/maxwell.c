#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <config.h>
#include <check.h>

#include "maxwell.h"

/* This file is has too many #ifdef's...blech. */

#define MIN2(a,b) ((a) < (b) ? (a) : (b))

maxwell_data *create_maxwell_data(int nx, int ny, int nz,
				  int *local_N, int *N_start, int *alloc_N,
				  int num_bands,
				  int num_fft_bands)
{
     int n[3] = {nx, ny, nz}, rank = (nz == 1) ? (ny == 1 ? 1 : 2) : 3, i;
     maxwell_data *d = 0;
     int fft_data_size;

#ifndef HAVE_FFTW
#  error Non-FFTW FFTs are not currently supported.
#endif
     
#ifdef HAVE_FFTW
     CHECK(sizeof(fftw_real) == sizeof(real),
	   "floating-point type is inconsistent with FFTW!");
#endif

     d = (maxwell_data *) malloc(sizeof(maxwell_data));
     CHECK(d, "out of memory");

     d->nx = nx;
     d->ny = ny;
     d->nz = nz;
     
     d->num_fft_bands = MIN2(num_bands, num_fft_bands);

     d->current_k[0] = d->current_k[1] = d->current_k[2] = 0.0;
     d->polarization = NO_POLARIZATION;

     d->last_dim_size = d->last_dim = n[rank - 1];

     /* ----------------------------------------------------- */
#ifndef HAVE_MPI 
     d->local_nx = nx; d->local_ny = ny;
     d->local_x_start = d->local_y_start = 0;
     *local_N = *alloc_N = nx * ny * nz;
     *N_start = 0;
     d->other_dims = *local_N / d->last_dim;

     d->fft_data = 0;  /* initialize it here for use in specific planner? */

#  ifdef HAVE_FFTW
#    ifdef SCALAR_COMPLEX
     d->fft_output_size = fft_data_size = nx * ny * nz;
     d->plan = fftwnd_create_plan_specific(rank, n, FFTW_FORWARD,
					   FFTW_ESTIMATE | FFTW_IN_PLACE,
					   (FFTW_COMPLEX*) d->fft_data,
					   3 * num_fft_bands,
					   (FFTW_COMPLEX*) d->fft_data,
					   3 * num_fft_bands);
     d->iplan = fftwnd_create_plan_specific(rank, n, FFTW_BACKWARD,
					    FFTW_ESTIMATE | FFTW_IN_PLACE,
					    (FFTW_COMPLEX*) d->fft_data,
					    3 * num_fft_bands,
					    (FFTW_COMPLEX*) d->fft_data,
					    3 * num_fft_bands);
#    else /* not SCALAR_COMPLEX */
     d->last_dim_size = 2 * (d->last_dim / 2 + 1);
     d->fft_output_size = fft_data_size = d->other_dims * d->last_dim_size;
     d->plan = rfftwnd_create_plan(rank, n, FFTW_FORWARD,
				   FFTW_ESTIMATE | FFTW_IN_PLACE,
				   REAL_TO_COMPLEX);
     d->plan = rfftwnd_create_plan(rank, n, FFTW_BACKWARD,
				   FFTW_ESTIMATE | FFTW_IN_PLACE,
				   COMPLEX_TO_REAL);
#    endif /* not SCALAR_COMPLEX */
#  endif /* HAVE_FFTW */

#else /* HAVE_MPI */
     /* ----------------------------------------------------- */

#  ifdef HAVE_FFTW

     CHECK(rank > 1, "rank < 2 MPI computations are not supported");

#    ifdef SCALAR_COMPLEX
     d->plan = fftwnd_mpi_create_plan(MPI_COMM_WORLD, rank, n,
				      FFTW_FORWARD,
				      FFTW_ESTIMATE | FFTW_IN_PLACE);
     d->iplan = fftwnd_mpi_create_plan(MPI_COMM_WORLD, rank, n,
				       FFTW_BACKWARD,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);

     fftwnd_mpi_local_sizes(plan, &d->local_nx, &d->local_x_start,
			    &d->local_ny, &d->local_y_start,
			    &fft_data_size);
     
     d->fft_output_size = nx * d->local_ny * nz;

#    else /* not SCALAR_COMPLEX */

     CHECK(rank > 2, "rank <= 2 MPI computations must use SCALAR_COMPLEX");

#  error rfftw MPI transforms not yet supported

#    endif /* not SCALAR_COMPLEX */
     
     *local_N = d->local_nx * ny * nz;
     *N_start = d->local_x_start * ny * nz;
     *alloc_N = *local_N;
     d->other_dims = *local_N / d->last_dim;

#  endif /* HAVE_FFTW */

#endif /* HAVE_MPI */
     /* ----------------------------------------------------- */

#ifdef HAVE_FFTW
     CHECK(d->plan && d->iplan, "FFTW plan creation failed");
#endif

     d->eps_inv = (symmetric_matrix*) malloc(sizeof(symmetric_matrix)
	                                     * d->fft_output_size);
     CHECK(d->eps_inv, "out of memory");

     /* a scratch output array is required because the "ordinary" arrays
	are not in a cartesian basis (or even a constant basis). */
     d->fft_data = (scalar*) malloc(sizeof(scalar) * 3
				    * num_fft_bands * fft_data_size);
     CHECK(d->fft_data, "out of memory");

     d->k_plus_G = (k_data*) malloc(sizeof(k_data) * *local_N);
     d->k_plus_G_normsqr = (real*) malloc(sizeof(real) * *local_N);
     CHECK(d->k_plus_G && d->k_plus_G_normsqr, "out of memory");

     d->eps_inv_mean = 1.0;

     d->local_N = *local_N;
     d->N_start = *N_start;
     d->alloc_N = *alloc_N;
     d->num_bands = num_bands;
     d->N = nx * ny * nz;

     return d;
}

void destroy_maxwell_data(maxwell_data *d)
{
     if (d) {

#ifdef HAVE_FFTW
#  ifdef HAVE_MPI
	  fftwnd_mpi_destroy_plan(d->plan);
	  fftwnd_mpi_destroy_plan(d->iplan);
#  else /* not HAVE_MPI */
#    ifdef SCALAR_COMPLEX
	  fftwnd_destroy_plan(d->plan);
	  fftwnd_destroy_plan(d->iplan);
#    else /* not SCALAR_COMPLEX */
	  rfftwnd_destroy_plan(d->plan);
	  rfftwnd_destroy_plan(d->iplan);
#    endif /* not SCALAR_COMPLEX */
#  endif /* not HAVE_MPI */
#endif /* HAVE FFTW */

	  free(d->eps_inv);
	  free(d->fft_data);
	  free(d->k_plus_G);
	  free(d->k_plus_G_normsqr);

	  free(d);
     }
}

/* compute a = b x c */
static void compute_cross(real *a0, real *a1, real *a2,
			  real b0, real b1, real b2,
			  real c0, real c1, real c2)
{
     *a0 = b1 * c2 - b2 * c1;
     *a1 = b2 * c0 - b0 * c2;
     *a2 = b0 * c1 - b1 * c0;
}

void update_maxwell_data_k(maxwell_data *d, real k[3],
			   real G1[3], real G2[3], real G3[3])
{
     int nx = d->nx, ny = d->ny, nz = d->nz;
     int cx = d->nx/2, cy = d->ny/2, cz = d->nz/2;
     k_data *kpG = d->k_plus_G;
     real *kpGn2 = d->k_plus_G_normsqr;
     int x, y, z;
     real kx = k[0], ky = k[1], kz = k[2];

     if (kx == 0.0 && ky == 0.0 && kz == 0.0) {
	  printf("detected zero k\n");
	  kx = 1e-5;
     }

     d->current_k[0] = kx;
     d->current_k[1] = ky;
     d->current_k[2] = kz;

     /* make sure current polarization is still valid: */
     set_maxwell_data_polarization(d, d->polarization);

     for (x = d->local_x_start; x < d->local_x_start + d->local_nx; ++x) {
	  int kxi = (x > cx) ? (x - nx) : x;
	  for (y = 0; y < ny; ++y) {
	       int kyi = (y > cy) ? (y - ny) : y;
	       for (z = 0; z < nz; ++z, kpG++, kpGn2++) {
		    int kzi = (z > cz) ? (z - nz) : z;
		    real kpGx, kpGy, kpGz, a, b, c, leninv;

		    /* Compute k+G: */
		    kpGx = kx + G1[0]*kxi + G2[0]*kyi + G3[0]*kzi;
		    kpGy = ky + G1[1]*kxi + G2[1]*kyi + G3[1]*kzi;
		    kpGz = kz + G1[2]*kxi + G2[2]*kyi + G3[2]*kzi;

		    a = kpGx*kpGx + kpGy*kpGy + kpGz*kpGz;
		    kpG->kmag = sqrt(a);
		    *kpGn2 = a;
		    
		    /* Now, compute the two normal vectors: */

		    if (kpGx == 0.0 && kpGy == 0.0) {
			 /* just put n in the x direction if k+G is in z: */
			 kpG->nx = 1.0;
			 kpG->ny = 0.0;
			 kpG->nz = 0.0;
		    }
		    else {
			 /* otherwise, let n = z x (k+G), normalized: */
			 compute_cross(&a, &b, &c,
				       0.0, 0.0, 1.0,
				       kpGx, kpGy, kpGz);
			 leninv = 1.0 / sqrt(a*a + b*b + c*c);
			 kpG->nx = a * leninv;
			 kpG->ny = b * leninv;
			 kpG->nz = c * leninv;
		    }

		    /* m = n x (k+G), normalized */
		    compute_cross(&a, &b, &c,
				  kpG->nx, kpG->ny, kpG->nz,
				  kpGx, kpGy, kpGz);
		    leninv = 1.0 / sqrt(a*a + b*b + c*c);
		    kpG->mx = a * leninv;
		    kpG->my = b * leninv;
		    kpG->mz = c * leninv;

#ifdef DEBUG
#define DOT(u0,u1,u2,v0,v1,v2) ((u0)*(v0) + (u1)*(v1) + (u2)*(v2))

		    /* check orthogonality */
		    CHECK(fabs(DOT(kpGx, kpGy, kpGz,
				   kpG->nx, kpG->ny, kpG->nz)) < 1e-6,
			  "vectors not orthogonal!");
		    CHECK(fabs(DOT(kpGx, kpGy, kpGz,
				   kpG->mx, kpG->my, kpG->mz)) < 1e-6,
			  "vectors not orthogonal!");
		    CHECK(fabs(DOT(kpG->mx, kpG->my, kpG->mz,
				   kpG->nx, kpG->ny, kpG->nz)) < 1e-6,
			  "vectors not orthogonal!");

		    /* check normalization */
		    CHECK(fabs(DOT(kpG->nx, kpG->ny, kpG->nz,
				   kpG->nx, kpG->ny, kpG->nz) - 1.0) < 1e-6,
			  "vectors not unit vectors!");
		    CHECK(fabs(DOT(kpG->mx, kpG->my, kpG->mz,
				   kpG->mx, kpG->my, kpG->mz) - 1.0) < 1e-6,
			  "vectors not unit vectors!");
#endif
	       }
	  }
     }
}

void set_maxwell_data_polarization(maxwell_data *d,
				   polarization_t polarization)
{
     if (d->current_k[2] != 0.0 || d->nz != 1)
	  polarization = NO_POLARIZATION;
     d->polarization = polarization;
}

maxwell_target_data *create_maxwell_target_data(maxwell_data *md, 
						real target_frequency)
{
     maxwell_target_data *d;

     d = (maxwell_target_data *) malloc(sizeof(maxwell_target_data));
     CHECK(d, "out of memory");

     d->d = md;
     d->target_frequency = target_frequency;

     d->T = create_evectmatrix(md->N, 2, md->num_bands, 
			       md->local_N, md->N_start, md->alloc_N);

     return d;
}

void destroy_maxwell_target_data(maxwell_target_data *d)
{
     if (d) {
	  destroy_evectmatrix(d->T);
	  free(d);
     }
}

/* Set Vinv = inverse of V, where both V and Vinv are symmetric matrices. */
static void sym_matrix_invert(symmetric_matrix *Vinv, symmetric_matrix V)
{
     double detinv;
     
     /* compute the determinant: */
     detinv = V.m00*V.m11*V.m22 - V.m02*V.m11*V.m02 +
	      2.0 * V.m01*V.m12*V.m02 -
	      V.m01*V.m01*V.m22 - V.m12*V.m12*V.m00;

     /* don't bother to check for singular matrices, as that shouldn't
	be possible in the context in which this is used */
     detinv = 1.0/detinv;
     
     Vinv->m00 = detinv * (V.m11*V.m22 - V.m12*V.m12);
     Vinv->m11 = detinv * (V.m00*V.m22 - V.m02*V.m02);
     Vinv->m22 = detinv * (V.m11*V.m00 - V.m01*V.m01);
     
     Vinv->m02 = detinv * (V.m01*V.m12 - V.m11*V.m02);
     Vinv->m01 = -detinv * (V.m01*V.m22 - V.m12*V.m02);
     Vinv->m12 = -detinv * (V.m00*V.m12 - V.m01*V.m02);
}

#define SMALL 1.0e-6

/* The following function initializes the dielectric tensor md->eps_inv,
   using the dielectric function epsilon(r, epsilon_data).

   epsilon is averaged over a rectangular mesh spanning the space between
   grid points; the size of the mesh is given by mesh_size.

   R1, R2, and R2 are the spatial lattice vectors, and are used to convert
   the discretization grid into spatial coordinates (with the origin at
   the (0,0,0) grid element.

   In most places, the dielectric tensor is equal to 1/eps * identity,
   but at dielectric interfaces it varies depending upon the polarization
   of the field (for faster convergence).  In particular, it depends upon
   the direction of the field relative to the surface normal vector, so
   we must compute the latter.  The surface normal is approximated by
   the "dipole moment" of the dielectric function over the mesh.

   Implementation note: md->eps_inv is chosen to have dimensions matching
   the output of the FFT.  Thus, its dimensions depend upon whether we are
   doing a real or complex and serial or parallel FFT. */

void set_maxwell_dielectric(maxwell_data *md,
			    int mesh_size[3],
			    real R[3][3],
			    dielectric_function epsilon,
			    void *epsilon_data)
{
     real s1, s2, s3, m1, m2, m3;  /* grid/mesh steps */
     real mesh_center[3];
     int i, j, k;
     int mesh_prod = mesh_size[0] * mesh_size[1] * mesh_size[2];
     real eps_inv_total = 0.0;
     int nx, ny, nz, local_ny, local_y_start, local_x_end;

     nx = md->nx; ny = md->ny; nz = md->nz;

     {
	  int mesh_divisions[3];

	  for (i = 0; i < 3; ++i) {
	       mesh_divisions[i] = mesh_size[i] <= 1 ? 1 : mesh_size[i] - 1;
	       mesh_center[i] = (mesh_size[i] - 1) * 0.5;
	  }

	       s1 = 1.0 / nx;
	       s2 = 1.0 / ny;
	       s3 = 1.0 / nz;
	       m1 = s1 / mesh_divisions[0];
	       m2 = s2 / mesh_divisions[1];
	       m3 = s3 / mesh_divisions[2];
     }


     /* Here we have different loops over the coordinates, depending
	upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and eps_index describing the corresponding index in 
	the array md->eps_inv[]. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI

     for (i = 0; i < nx; ++i)
	  for (j = 0; j < ny; ++j)
	       for (k = 0; k < nz; ++k)
     {
#         define i2 i
#         define j2 j
#         define k2 k
	  int eps_index = ((i * ny + j) * nz + k);

#  else /* HAVE_MPI */

     local_ny = md->fft_output_size / (nx * nz);
     local_y_start = md->fft_output_N_start / (nx * nz);

     /* first two dimensions are transposed in MPI output: */
     for (j = 0; j < local_ny; ++j)
          for (i = 0; i < nx; ++i)
	       for (k = 0; k < nz; ++k)
     {
#         define i2 i
	  int j2 = j + local_y_start;
#         define k2 k
	  int eps_index = ((j * nx + i) * nz + k);

#  endif

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

#  else /* HAVE_MPI */

#    error not yet implemented!

#  endif

#endif /* not SCALAR_COMPLEX */

	  int mi, mj, mk;
	  real eps_mean = 0.0, eps_inv_mean = 0.0, norm_len;
	  real norm0, norm1, norm2, moment0=0.0, moment1=0.0, moment2=0.0;
	  real del0, del1, del2;
	  short means_different_p;
	  
	  for (mi = 0; mi < mesh_size[0]; ++mi)
	       for (mj = 0; mj < mesh_size[1]; ++mj)
		    for (mk = 0; mk < mesh_size[2]; ++mk) {
			 real r[3], eps;
			 
			 del0 = (mi - mesh_center[0]) * m1;
			 r[0] = i2 * s1 + del0;
			 del1 = (mj - mesh_center[1]) * m2;
			 r[1] = j2 * s2 + del1;
			 del2 = (mk - mesh_center[2]) * m3;
			 r[2] = k2 * s3 + del2;
			 
			 eps = epsilon(r, epsilon_data);
			 
			 eps_mean += eps;
			 eps_inv_mean += 1.0 / eps;
			 
			 moment0 += eps * del0;
			 moment1 += eps * del1;
			 moment2 += eps * del2;
		    }
     
	  eps_mean = eps_mean / mesh_prod;
	  eps_inv_mean = mesh_prod / eps_inv_mean;

	  means_different_p = fabs(eps_mean - eps_inv_mean) > SMALL;

	  if (means_different_p) {
	       /* need to convert moment from lattice to cartesian coords: */
	       norm0 = R[0][0]*moment0 + R[1][0]*moment1 + R[2][0]*moment2;
	       norm1 = R[0][1]*moment0 + R[1][1]*moment1 + R[2][1]*moment2;
	       norm2 = R[0][2]*moment0 + R[1][2]*moment1 + R[2][2]*moment2;
	  
	       norm_len = sqrt(norm0*norm0 + norm1*norm1 + norm2*norm2);
	  }
	  
	  if (means_different_p && norm_len > SMALL) {
	       symmetric_matrix eps;
	       
	       norm0 /= norm_len;
	       norm1 /= norm_len;
	       norm2 /= norm_len;
	       
	       /* compute effective dielectric tensor: */
	       
	       eps.m00 = (eps_inv_mean-eps_mean) * norm0*norm0 + eps_mean;
	       eps.m11 = (eps_inv_mean-eps_mean) * norm1*norm1 + eps_mean;
	       eps.m22 = (eps_inv_mean-eps_mean) * norm2*norm2 + eps_mean;
	       eps.m01 = (eps_inv_mean-eps_mean) * norm0*norm1;
	       eps.m02 = (eps_inv_mean-eps_mean) * norm0*norm2;
	       eps.m12 = (eps_inv_mean-eps_mean) * norm1*norm2;

	       sym_matrix_invert(&md->eps_inv[eps_index], eps);
	  }
	  else { /* undetermined normal vector and/or constant eps */
	       md->eps_inv[eps_index].m00 = 1.0 / eps_mean;
	       md->eps_inv[eps_index].m11 = 1.0 / eps_mean;
	       md->eps_inv[eps_index].m22 = 1.0 / eps_mean;
	       md->eps_inv[eps_index].m01 = 0.0;
	       md->eps_inv[eps_index].m02 = 0.0;
	       md->eps_inv[eps_index].m12 = 0.0;
	  }
	  
	  eps_inv_total += (md->eps_inv[eps_index].m00 + 
			    md->eps_inv[eps_index].m11 + 
			    md->eps_inv[eps_index].m22);
     }  /* end of loop body */
     
     md->eps_inv_mean = eps_inv_total / (3 * md->fft_output_size);
}
