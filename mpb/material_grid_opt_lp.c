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

/* optimization using SDP */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>

#include <scalar.h>
#include <blasglue.h>
#include <time.h>

#include "mpb.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MOSEK_H
#include <mosek.h>    /* Include the MOSEK definition file.  */
#else
#define MSKAPI
#define MSKCONST const
#endif

/* FIXME: we don't require cblas header elsewhere — change to
          use Fortran API via blasglue.c? */
#ifdef HAVE_CBLAS
#include <cblas.h>
#else
typedef enum {CblasRowMajor=101, CblasColMajor=102} CBLAS_LAYOUT;
typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;
typedef enum {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
static void cblas_zgemv(const CBLAS_LAYOUT __Order, const CBLAS_TRANSPOSE __TransA, const int __M, const int __N, const void *__alpha, const void *__A, const int __lda, const void *__X, const int __incX, const void *__beta, void *__Y, const int __incY) { CHECK(0, "missing cblas.h"); }
static void cblas_dgemv(const CBLAS_LAYOUT __Order, const CBLAS_TRANSPOSE __TransA, const int __M, const int __N, const double __alpha, const double *__A, const int __lda, const double *__X, const int __incX, const double __beta, double *__Y, const int __incY) { CHECK(0, "missing cblas.h"); }
static void cblas_zdotc_sub(const int __N, const void *__X, const int __incX, const void *__Y, const int __incY, void *__dotc) { CHECK(0, "missing cblas.h"); }
static double cblas_ddot(const int __N, const double *__X, const int __incX, const double *__Y, const int __incY) { CHECK(0, "missing cblas.h"); }
static void cblas_dsymv(const CBLAS_LAYOUT __Order, const CBLAS_UPLO __Uplo, const int __N, const double __alpha, const double *__A, const int __lda, const double *__X, const int __incX, const double __beta, double *__Y, const int __incY) { CHECK(0, "missing cblas.h"); }
#endif

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))
typedef int bool;
#define true 1
#define false 0
#define NUM_THREADS 4

extern void debug_check_memory_leaks(void);

/**************************************************************************/

/* Note: the usage of this func is very limited!
   modify the epsilon with random perturbations; we assume meshsize = 1,
   hence, eps and eps_inv tensor are diagonal. should only be called after
   init-params.  (Guile-callable.) */
void randomize_epsilon(number delta)
{
  int i, N, ngrids;
  double allowance, eps_mean, eps_min, eps_max, eps_diff, eps_mid, eps_val;
  material_grid *mg;
  symmetric_matrix eps_inv;
  real m00, m11, m22, m01, m02, m12;

     if (!mdata)
	  return;
     mpi_one_printf("Randomizing epsilon...\n");

     mg = get_material_grids(geometry, &ngrids);
     eps_max = mg->epsilon_max;
     eps_min = mg->epsilon_min;
     eps_diff = eps_max - eps_min;
     eps_mid = (eps_max + eps_min)/2;

     N = mdata->fft_output_size;
     mpi_one_printf("total number of epsilon N = %d\n", N);
     allowance = delta*N*eps_diff;
     while(allowance>0){
       i = floor(rand()%N);

       eps_inv = mdata->eps_inv[i];

       m01 = eps_inv.m01; m02 = eps_inv.m02; m12 = eps_inv.m12;
       CHECK(m01 == 0.0 && m02 == 0.0 && m12 == 0.0, "epsilon tensor not diagnal");

       m00 = eps_inv.m00; m11 = eps_inv.m11; m22 = eps_inv.m22;
       eps_mean = 3.0/(m00+m11+m22);
       eps_val = (eps_mean > eps_mid) ? eps_min : eps_max;
       allowance =  allowance - fabs(eps_val - eps_mean);

       if(allowance>0)
	 {
	   mpi_one_printf("updating element %d, allowane left %f\n", i, allowance);
	   mpi_one_printf("original eps_mean = %f, update to eps_val %f\n", eps_mean, eps_val);
	   mdata->eps_inv[i].m00 = 1.0/eps_val;
	   mdata->eps_inv[i].m11 = 1.0/eps_val;
	   mdata->eps_inv[i].m22 = 1.0/eps_val;

	 }



     }

}
/**************************************************************************/

void randomize_meshgrid(number delta)
{
  int i, ngrids, ntot;
  double allowance, mg_mean, mg_min, mg_max, mg_diff, mg_mid, mg_val;
  material_grid *grids;
  double *u;
  mpi_one_printf("Randomizing meshgrid...\n");

  grids = get_material_grids(geometry, &ngrids);
  ntot = material_grids_ntot(grids, ngrids);
  u = (double *) malloc(sizeof(double) * ntot);
  material_grids_get(u, grids, ngrids);

     mg_max = 1.0;
     mg_min = 0.0;
     mg_diff = mg_max - mg_min;
     mg_mid = (mg_max + mg_min)/2;

     allowance = delta*ntot*mg_diff;
     while(allowance>0){
       i = floor(rand()%ntot);

       mg_val = (u[i] > mg_mid) ? mg_min : mg_max;
       allowance =  allowance - fabs(mg_val - u[i]);

       if(allowance>0)
	 u[i] = mg_val;
     }
	/* update u and epsilon */
        /* mpi_one_printf("Allowance left : %g % \n", allowance/(delta*ntot*mg_diff)*100.0); */
	material_grids_set(u, grids, ngrids);
	reset_epsilon();

}

/**************************************************************************/
/* dynamically expanding array struct */

typedef struct {
  double *data;
  int used;
  int size;
} double_array;


static void initArray(double_array *a, int initialSize) {
  a->data = (double *)malloc(initialSize * sizeof(double));
  a->used = 0;
  a->size = initialSize;
}

static void freeArray(double_array *a) {
  free(a->data);
  a->data = NULL;
  a->used = a->size = 0;
}

void print_solsta(int solsta, char *solstastr)
{
#ifdef HAVE_MOSEK_H
  switch(solsta)
    {
    /* case MSK_SOL_STA_BEGIN: */
    /*   sprintf (solstastr, "MSK_SOL_STA_BEGIN = $%d", 0);  */
    /* case MSK_SOL_STA_END: */
    /*   sprintf (solstastr, "MSK_SOL_STA_END = %d", 16);  */

    case MSK_SOL_STA_UNKNOWN:
      sprintf (solstastr, "MSK_SOL_STA_UNKNOWN = %d", 0);
      break;
    case MSK_SOL_STA_OPTIMAL:
      sprintf (solstastr, "MSK_SOL_STA_OPTIMAL = %d", 1);
      break;
    case MSK_SOL_STA_PRIM_FEAS:
      sprintf (solstastr, "MSK_SOL_STA_PRIM_FEAS = %d", 2);
      break;
    case MSK_SOL_STA_DUAL_FEAS:
      sprintf (solstastr, "MSK_SOL_STA_DUAL_FEAS = %d", 3);
      break;
    case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      sprintf (solstastr, "MSK_SOL_STA_PRIM_AND_DUAL_FEAS = %d", 4);
      break;
    case MSK_SOL_STA_PRIM_INFEAS_CER:
      sprintf (solstastr, "MSK_SOL_STA_PRIM_INFEAS_CER = %d", 5);
      break;
    case MSK_SOL_STA_DUAL_INFEAS_CER:
      sprintf (solstastr, "MSK_SOL_STA_DUAL_INFEAS_CER = %d", 6);
      break;
    /* case : */
    /*   sprintf (solstastr, " = %d", 7);  */
    case MSK_SOL_STA_NEAR_OPTIMAL:
      sprintf (solstastr, "MSK_SOL_STA_NEAR_OPTIMAL = %d", 8);
      break;
    case MSK_SOL_STA_NEAR_PRIM_FEAS:
      sprintf (solstastr, "MSK_SOL_STA_NEAR_PRIM_FEAS = %d", 9);
      break;
    case MSK_SOL_STA_NEAR_DUAL_FEAS:
      sprintf (solstastr, "MSK_SOL_STA_NEAR_DUAL_FEAS = %d", 10);
      break;
    case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
      sprintf (solstastr, "MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS = %d", 11);
      break;
    case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
      sprintf (solstastr, "MSK_SOL_STA_NEAR_PRIM_INFEAS_CER = %d", 12);
      break;
    case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      sprintf (solstastr, "MSK_SOL_STA_NEAR_DUAL_INFEAS_CER = %d", 13);
      break;
    case MSK_SOL_STA_INTEGER_OPTIMAL:
      sprintf (solstastr, "MSK_SOL_STA_INTEGER_OPTIMAL = %d", 14);
      break;
    case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
      sprintf (solstastr, "MSK_SOL_STA_NEAR_INTEGER_OPTIMAL = %d", 15);
      break;
    default:
      sprintf (solstastr, "Other solution status = %d", solsta);
      break;
    }
#else
	sprintf (solstastr, "Solution status = %d", solsta);
#endif
}

/*  returns all combinations of "n" number of coefficients "c", where c */
/*  adds up to "k" */
/*  used in SDP2LPcoef.m */
/*  total number of ways to write k = c_1 + c_2 + ... + c_n for non-negative */
/* c_i: (^{n+k-1}_k)? */
static double_array *sumCombination( int k, int n)
 {

   int i,row,col,len;
   double_array *c,  *ctmpt;

   c = (double_array *)malloc(sizeof(double_array));

   if(n==1)
     {
       initArray(c, 1);
       c->data[0] = 1.0;
       c->used = 1;
       return c;
     }
   else if(n==2)
     {
       initArray(c, 2*(k+1));
       for(i=0; i<=k; ++i)
	 { c->data[2*i] = k-i; c->data[2*i+1] =i; }
       c->used = 2*(k+1);
       return c;
     }
   else if(n>2)
     {
       initArray(c, 0);
       ctmpt = (double_array *)malloc(sizeof(double_array));

       for(i=0; i<=k; ++i)
   	 {

   	   /* initArray(ctmpt, 2*(k-i+1)); */
   	   ctmpt = sumCombination( k-i, n-1);
	   len = ctmpt->size/(n-1);

   	   if (c->used + ctmpt->used + len > c->size)
   	     {
   	       c->size = c->used + ctmpt->used + len;
   	       c->data = (double *)realloc(c->data, c->size * sizeof(double));
   	     }

   	   for(row=0; row<len; ++row)
   	     {
   	       for(col=0; col<n-1 ; ++col)
   		 c->data[c->used + row*n + col] = ctmpt->data[row*(n-1)+col];
   	       c->data[c->used + row*n + n-1] = i;
   	     }
   	   c->used = c->used + ctmpt->used + len ;

   	 }
       freeArray(ctmpt);
       return c;
     }

 }

/**************************************************************************/

static void printmatrix(scalar *A, int m, int n)
{
  int i, j;

  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
#ifdef SCALAR_COMPLEX
      mpi_one_printf("  (%6.6f,%6.6f)", (double)A[i*n + j].re, (double)A[i*n + j].im);
#else
       mpi_one_printf("  %6.6f", (double)A[i*n + j]);
#endif
    }
     mpi_one_printf("\n");
  }
}

/* print subset of matrix */
static void printmatrix_sub(scalar *A, int m, int n, int lda)
{
  int i, j;

  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
#ifdef SCALAR_COMPLEX
      mpi_one_printf("  (%6.6f,%6.6f)", (double)A[i*lda + j].re, (double)A[i*lda + j].im);
#else
      mpi_one_printf("  %6.6f", (double)A[i*lda + j]);
#endif
    }
    mpi_one_printf("\n");
  }
}

static void printmatrix_double(double *A, int m, int n)
{
  int i, j;

  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j)
      mpi_one_printf("  %6.6f", (double)A[i*n + j]);
    mpi_one_printf("\n");
  }
}


/**************************************************************************/
static void MSKAPI printstr(void *handle,
                            MSKCONST char str[])
{
    mpi_one_printf("%s",str);
} /* printstr */

/**************************************************************************/
/* res = v1'*v2 */
static scalar_complex cvector_inproduct(scalar_complex *cv1, double *v2, int len)
{
    int i;
    scalar_complex res;
    CASSIGN_SCALAR(res,0.0,0.0);

    for (i = 0; i < len; ++i) {
        res.re += cv1[i].re*v2[i];
        res.im += cv1[i].im*v2[i];
    }
    return res;

}

static scalar_complex cvector_inproductT(scalar_complex *cv1, double *v2, int len, int stride)
{
    int i;
    scalar_complex res;
    CASSIGN_SCALAR(res,0.0,0.0);

    for (i = 0; i < len; ++i) {
        res.re += cv1[i*stride].re*v2[i];
        res.im += cv1[i*stride].im*v2[i];
    }
    return res;

}

static void cvector_product(scalar_complex *cv1, double *v2, int len, int stride1, scalar_complex *cv2, int stride2)
{
    int i;

    for (i = 0; i < len; ++i) {
        cv2[i*stride2].re = cv1[i*stride1].re*v2[i];
        cv2[i*stride2].im = cv1[i*stride1].im*v2[i];
    }

}


static void get_Dfield(int which_band,scalar_complex *cfield,int cur_num_bands,real scale)
{
    int i, N;
    maxwell_compute_d_from_H(mdata, H, cfield, which_band-1, cur_num_bands);
    N = mdata->fft_output_size;
    /*Dfield/Efield in position space, hence dimension = fft_output_size*3 = (nx*ny*nz)*3*/
    for (i = 0; i < 3*N; ++i) {
        cfield[i].re *= scale;
        cfield[i].im *= scale;
    }
}

static void get_Efield_from_Dfield1(const scalar_complex *dfield, scalar_complex *efield, int cur_num_bands)
{

    int i,b;

    /*Dfield/Efield in position space, hence dimension = fft_output_size*3*/
    int N = mdata->fft_output_size;

    for (i = 0; i < N; ++i) {
        symmetric_matrix eps_inv = mdata->eps_inv[i];
        for (b = 0; b < cur_num_bands; ++b) {
            int ib = 3 * (i * cur_num_bands + b);
            assign_symmatrix_vector(&efield[ib], eps_inv, &dfield[ib]);
        }
    }

}

static void get_Efield_from_Dfield(scalar_complex *cfield, int cur_num_bands)
{
    maxwell_compute_e_from_d(mdata, cfield, cur_num_bands);
}

static void get_Efield(int which_band,scalar_complex *cfield,int cur_num_bands, real scale)
{
    get_Dfield(which_band,cfield,cur_num_bands,scale);
    get_Efield_from_Dfield(cfield,cur_num_bands);
}


static void deps_du(double *v, double scalegrad, const material_grid *grids, int ngrids)
{
    int i, j, k, n1, n2, n3, n_other, n_last, rank, last_dim;
#ifdef HAVE_MPI
    int local_n2, local_y_start, local_n3;
#endif
    real s1, s2, s3, c1, c2, c3;

    int ntot = material_grids_ntot(grids, ngrids);

    n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;
    n_other = mdata->other_dims;
    n_last = mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
    last_dim = mdata->last_dim;
    rank = (n3 == 1) ? (n2 == 1 ? 1 : 2) : 3;

    s1 = geometry_lattice.size.x / n1;
    s2 = geometry_lattice.size.y / n2;
    s3 = geometry_lattice.size.z / n3;
    c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5;
    c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
    c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

    /* Here we have different loops over the coordinates, depending
    upon whether we are using complex or real and serial or
        parallel transforms.  Each loop must define, in its body,
        variables (i2,j2,k2) describing the coordinate of the current
        point, and "index" describing the corresponding index in
    the curfield array.

        This was all stolen from fields.c...it would be better
        if we didn't have to cut and paste, sigh. */

#ifdef SCALAR_COMPLEX

#  ifndef HAVE_MPI

    for (i = 0; i < n1; ++i)
        for (j = 0; j < n2; ++j)
            for (k = 0; k < n3; ++k)
            {
                int i2 = i, j2 = j, k2 = k;
                int index = ((i * n2 + j) * n3 + k);

#  else /* HAVE_MPI */

    local_n2 = mdata->local_ny;
    local_y_start = mdata->local_y_start;

    /* first two dimensions are transposed in MPI output: */
    for (j = 0; j < local_n2; ++j)
        for (i = 0; i < n1; ++i)
            for (k = 0; k < n3; ++k)
            {
                int i2 = i, j2 = j + local_y_start, k2 = k;
                int index = ((j * n1 + i) * n3 + k);

#  endif /* HAVE_MPI */

#else /* not SCALAR_COMPLEX */

#  ifndef HAVE_MPI

    for (i = 0; i < n_other; ++i)
        for (j = 0; j < n_last; ++j)
        {
            int index = i * n_last + j;
            int i2, j2, k2;
            switch (rank) {
            case 2: i2 = i; j2 = j; k2 = 0; break;
            case 3: i2 = i / n2; j2 = i % n2; k2 = j; break;
            default: i2 = j; j2 = k2 = 0;  break;
            }

#  else /* HAVE_MPI */

    local_n2 = mdata->local_ny;
    local_y_start = mdata->local_y_start;

    /* For a real->complex transform, the last dimension is cut in
    half.  For a 2d transform, this is taken into account in local_ny
    already, but for a 3d transform we must compute the new n3: */
    if (n3 > 1)
        local_n3 = mdata->last_dim_size / 2;
    else
        local_n3 = 1;

    /* first two dimensions are transposed in MPI output: */
    for (j = 0; j < local_n2; ++j)
        for (i = 0; i < n1; ++i)
            for (k = 0; k < local_n3; ++k)
            {
#         define i2 i
                int j2 = j + local_y_start;
#         define k2 k
                int index = ((j * n1 + i) * local_n3 + k);

#  endif /* HAVE_MPI */

#endif /* not SCALAR_COMPLEX */

                {
                    vector3 p;

                    p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;

                    material_grids_addgradient_point(
                                v + index*ntot, p, scalegrad, grids,ngrids);

#ifndef SCALAR_COMPLEX
                    {
                        int last_index;
#  ifdef HAVE_MPI
                        if (n3 == 1)
                            last_index = j + local_y_start;
                        else
                            last_index = k;
#  else
                        last_index = j;
#  endif

                        if (last_index != 0 && 2*last_index != last_dim) {
                            int i2c, j2c, k2c;
                            i2c = i2 ? (n1 - i2) : 0;
                            j2c = j2 ? (n2 - j2) : 0;
                            k2c = k2 ? (n3 - k2) : 0;
                            p.x = i2c * s1 - c1;
                            p.y = j2c * s2 - c2;
                            p.z = k2c * s3 - c3;

                            material_grids_addgradient_point(
                                        v+index*ntot, p, scalegrad, grids,ngrids);
                        }
                    }
#endif /* !SCALAR_COMPLEX */

                }

            }
}
/*************/
static scalar_complex compute_fields_energy(scalar_complex *field1, scalar_complex *field2, bool update)
{
    int i, last_dim, last_dim_stored, nx, nz, local_y_start;
    scalar_complex energy_sum;
    scalar_complex *product = (scalar_complex *) field1;

    last_dim = mdata->last_dim;
    last_dim_stored =
            mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
    nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

    energy_sum.re = 0;
    energy_sum.im = 0;
    for (i = 0; i < mdata->fft_output_size; ++i) {
        scalar_complex field[3];
        real comp_sqr = 0;
        real comp_sqri = 0;

        field[0] =   field1[3*i];
        field[1] = field1[3*i+1];
        field[2] = field1[3*i+2];

        comp_sqr += field[0].re *   field2[3*i].re + field[0].im *   field2[3*i].im;
        comp_sqr += field[1].re * field2[3*i+1].re + field[1].im * field2[3*i+1].im;
        comp_sqr += field[2].re * field2[3*i+2].re + field[2].im * field2[3*i+2].im;

        comp_sqri += field[0].re *   field2[3*i].im - field[0].im *   field2[3*i].re;
        comp_sqri += field[1].re * field2[3*i+1].im - field[1].im * field2[3*i+1].re;
        comp_sqri += field[2].re * field2[3*i+2].im - field[2].im * field2[3*i+2].re;

        /* Note: here, we write to product[i]; this is
         safe, even though product is aliased to field1,
         since product[i] is guaranteed to come at or before
         field1[i] (which we are now done with). */

        energy_sum.re += comp_sqr;
        energy_sum.im += comp_sqri;
        if (update == true)
            CASSIGN_SCALAR(product[i],comp_sqr,comp_sqri);



#ifndef SCALAR_COMPLEX
        /* most points need to be counted twice, by rfftw output symmetry: */
        {
            int last_index;
#  ifdef HAVE_MPI
            if (nz == 1) /* 2d calculation: 1st dim. is truncated one */
                last_index = i / nx + local_y_start;
            else
                last_index = i % last_dim_stored;
#  else
            last_index = i % last_dim_stored;
#  endif /*HAVE_MPI*/
            if (last_index != 0 && 2*last_index != last_dim) {
                energy_sum.re += comp_sqr;
                energy_sum.im -= comp_sqri;
            }
        }
#endif
    }

    mpi_allreduce_1(&energy_sum.re, real, SCALAR_MPI_TYPE,
                    MPI_SUM, MPI_COMM_WORLD);
    mpi_allreduce_1(&energy_sum.im, real, SCALAR_MPI_TYPE,
                    MPI_SUM, MPI_COMM_WORLD);
    return energy_sum;
}

/********************/


// returns transposed Asp. Key is in the way Asp is updated (stride = final count, Asp starts from the offset (=current count))
static void material_grids_SPt(scalar_complex *Asp, const double *depsdu, const double *u, real scalegrad, int ntot, int band1, int band2, int stride)
{

    int i,ui;
    scalar_complex *field1,*field2,A0, Aisum, Aitmpt, fieldsum1;
    int cur_num_bands = 1;

    CHECK(band1 <= num_bands && band2 <= num_bands, "reducedA0 called for uncomputed band\n");
    field1 = (scalar_complex *) malloc(sizeof(scalar_complex) * mdata->fft_output_size*3);
    field2 = (scalar_complex *) malloc(sizeof(scalar_complex) * mdata->fft_output_size*3);
    /* Ai = (scalar_complex *) malloc(sizeof(scalar_complex) * ntot); */

    /* compute A0: Dfield_1'*Efield_2 */
    if (band1)
        get_Dfield(band1, field1, cur_num_bands, 1.0);

    if (band2!=band1)
        get_Efield(band2,field2, cur_num_bands, 1.0);
    else
        get_Efield_from_Dfield1(field1, field2, cur_num_bands);

    A0 = compute_fields_energy(field1,field2,false);
    /* A0 has been reduced (mpi_allreduce_1) in the compute_field_energy routine */

    scalegrad *= 1.0/H.N;
    CASSIGN_SCALAR(A0,A0.re*scalegrad,A0.im*scalegrad);

    /* compute Ai: - Efield_1'*depsdu*Efield_2 */
    if (band1!=band2)
        get_Efield_from_Dfield(field1, cur_num_bands);
    else
    {
        free(field1);
        field1 = field2;
    }

    fieldsum1 = compute_fields_energy(field1,field2,true);
    scalegrad *= -1.0;

    for (ui = 0; ui < ntot; ++ui)
    {
        CASSIGN_SCALAR(Aitmpt, 0.0,0.0);

        for (i = 0; i < mdata->fft_output_size; ++i)
        {
            Aitmpt.re += depsdu[i*ntot+ui]*field1[i].re;
            Aitmpt.im += depsdu[i*ntot+ui]*field1[i].im;
        }
        CASSIGN_SCALAR(Asp[ui*stride],Aitmpt.re*scalegrad,Aitmpt.im*scalegrad);

    /* field1, field2, and despdu (of size fft_output_size = nx*local_y*nz) only have the local block of data, and their products
       should be summed up (mpi_reduce) to account for global info. local_ny ~= ny/mpi_comm_size */

	mpi_allreduce_1(&Asp[ui*stride].re, real, SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
	mpi_allreduce_1(&Asp[ui*stride].im, real, SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);

    }

    Aisum = cvector_inproductT(Asp, u, ntot, stride);
    Asp[ntot*stride].re =  A0.re - Aisum.re;
    Asp[ntot*stride].im =  A0.im - Aisum.im;

    free(field2);
}
// returns transposed Asp. Key is in the way Asp(Ai) is updated (stride = final count, Asp starts from the offset (=current count))
static void material_grids_SPt_blas(scalar_complex *Ai, scalar_complex *depsdu_sc, scalar_complex *u_sc, real scalegrad, int ntot, int band1, int band2, int stride)
{

  int M,N;
  scalar_complex alpha, beta;
  scalar_complex *field1,*field2, A0, Aisum, foo;
  int cur_num_bands = 1;

    CHECK(band1 <= num_bands && band2 <= num_bands, "reducedA0 called for uncomputed band\n");
    field1 = (scalar_complex *) malloc(sizeof(scalar_complex) * mdata->fft_output_size*3);
    field2 = (scalar_complex *) malloc(sizeof(scalar_complex) * mdata->fft_output_size*3);

    /* compute A0: Dfield_1'*Efield_2 */
    if (band1)
        get_Dfield(band1, field1, cur_num_bands, 1.0);

    if (band2!=band1)
        get_Efield(band2,field2, cur_num_bands, 1.0);
    else
        get_Efield_from_Dfield1(field1, field2, cur_num_bands);

    A0 = compute_fields_energy(field1,field2,false);
    /* A0 has been reduced (mpi_allreduce_1) in the compute_field_energy routine */

    scalegrad *= 1.0/H.N;
    CASSIGN_SCALAR(A0,A0.re*scalegrad,A0.im*scalegrad);

    /* mpi_one_printf("(band1, band2) = (%d, %d), A0 = %g + i %g\n", band1, band2, A0.re, A0.im); */

    /* compute Ai: - Efield_1'*depsdu*Efield_2 */
    if (band1!=band2)
        get_Efield_from_Dfield(field1, cur_num_bands);
    else
    {
        free(field1);
        field1 = field2;
    }

    /* updata field1 = Efield_1.*Efield_2*/
    foo = compute_fields_energy(field1,field2,true);
    scalegrad *= -1.0;

    M = mdata->fft_output_size; N = ntot;
    CASSIGN_SCALAR(alpha, scalegrad, 0.0);
    CASSIGN_SCALAR(beta, 0.0, 0.0);

    cblas_zgemv(CblasRowMajor, CblasTrans, M, N, &alpha, depsdu_sc, N, field1, 1, &beta, Ai, stride);

    /* field1, field2, and despdu (of size fft_output_size = nx*local_y*nz) only have the local block of data, and their products
       should be summed up (mpi_reduce) to account for global info. local_ny ~= ny/mpi_comm_size */
    int i;
    for (i=0; i<ntot; i++)
      {
	mpi_allreduce_1(&Ai[i*stride].re, real, SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
	mpi_allreduce_1(&Ai[i*stride].im, real, SCALAR_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
      }

    cblas_zdotc_sub(ntot, Ai, stride, u_sc, 1, &Aisum);


    Ai[ntot*stride].re =  A0.re - Aisum.re;
    Ai[ntot*stride].im =  A0.im - Aisum.im;

    free(field2);
}

static void SDP2LP(double_array *Alin, scalar_complex *Avec, double_array *avec, int spdim, int n)
{

  int ridx,cidx,count,jn,ja;
  double *Amat,*vectmpt, *vecnorm;

  int stride = spdim*(spdim+1)/2;
  /* size of reduced matrices, converted from complex to real */
  int matsize = spdim*2;
  /* number of approximating vectors */
  int na = avec->used/matsize;

  Amat = (double *) calloc(matsize*matsize,sizeof(double));
  vectmpt = (double *) calloc(matsize,sizeof(double));
  vecnorm = (double *) malloc(sizeof(double)*na);

  for(ja = 0; ja < na; ++ja)
    vecnorm[ja] = cblas_ddot(matsize, avec->data+ja*matsize, 1, avec->data+ja*matsize, 1);

  for (jn = 0; jn < n; ++jn)
    {
      count = jn*stride;
      for(ridx=0; ridx<spdim; ++ridx)
	{
	  for(cidx=0; cidx<=ridx; ++cidx)
	    {
	      Amat[ridx*matsize+cidx] = Avec[count].re;
	      Amat[(ridx+spdim)*matsize+cidx+spdim] = Avec[count].re;
	      Amat[(ridx+spdim)*matsize+cidx] = Avec[count].im;
	      if (ridx!=cidx){
		Amat[(cidx+spdim)*matsize+ridx] = -Avec[count].im;
	      }
	      ++count;
	    }
	}

      for(ja = 0; ja < na; ++ja)
	{
	  cblas_dsymv(CblasRowMajor,CblasLower,matsize,1/vecnorm[ja],Amat,matsize,avec->data+ja*matsize,1,0,vectmpt,1);
	  Alin->data[Alin->used+ja*n+jn] = cblas_ddot(matsize, vectmpt, 1, avec->data+ja*matsize, 1);
	  }
    }
  /* mpi_one_printf("Alin = \n"); */
  /* printmatrix_double(Alin->data+Alin->used,na,n); */

  Alin->used = Alin->used + na*n;
  free(Amat);
  free(vectmpt);
  free(vecnorm);
}

static int SDP2LP_v1(double_array *Alin, scalar_complex *Avec, int spdim, int n, int kk)
{

  int ridx,cidx,count,jn,ja;
  double *Amat,*vectmpt, *vecnorm;
  double_array *avec;

  int stride = spdim*(spdim+1)/2;
  /* size of reduced matrices, converted from complex to real */
  int matsize = spdim*2;

  avec = sumCombination(kk, matsize);

  /* number of approximating vectors */
  int na = avec->used/matsize;

  if (Alin->size < n*na + Alin->used){
    Alin->size = n*na + Alin->used;
    Alin->data = (double *)realloc(Alin->data, Alin->size * sizeof(double));
  }

  Amat = (double *) calloc(matsize*matsize,sizeof(double));
  vectmpt = (double *) calloc(matsize,sizeof(double));
  vecnorm = (double *) malloc(sizeof(double)*na);

  for(ja = 0; ja < na; ++ja)
    vecnorm[ja] = cblas_ddot(matsize, avec->data+ja*matsize, 1, avec->data+ja*matsize, 1);

  for (jn = 0; jn < n; ++jn)
    {
      count = jn*stride;
      for(ridx=0; ridx<spdim; ++ridx)
	{
	  for(cidx=0; cidx<=ridx; ++cidx)
	    {
	      Amat[ridx*matsize+cidx] = Avec[count].re;
	      Amat[(ridx+spdim)*matsize+cidx+spdim] = Avec[count].re;
	      Amat[(ridx+spdim)*matsize+cidx] = Avec[count].im;
	      if (ridx!=cidx){
		Amat[(cidx+spdim)*matsize+ridx] = -Avec[count].im;
	      }
	      ++count;
	    }
	}

      for(ja = 0; ja < na; ++ja)
	{
	  cblas_dsymv(CblasRowMajor,CblasLower,matsize,1/vecnorm[ja],Amat,matsize,avec->data+ja*matsize,1,0,vectmpt,1);
	  Alin->data[Alin->used+ja*n+jn] = cblas_ddot(matsize, vectmpt, 1, avec->data+ja*matsize, 1);
	  }
    }

  Alin->used = Alin->used + na*n;
  free(Amat);
  free(vectmpt);
  free(vecnorm);
  freeArray(avec);

  return na;
}



static void Ad2Ad(double *Amat, const double *Avec,int spdim, int n)
{

  int ridx,cidx,count,jn;

  int stride = spdim*(spdim+1)/2;
  /* size of reduced matrices, converted from complex to real */
  int matsize = spdim*2;
  int matsize2 = matsize*matsize;

  memset(Amat, 0, sizeof(double)*matsize2);
  for (jn = 0; jn < n; ++jn)
    {
      count = jn*stride;
      for(ridx=0; ridx<spdim; ++ridx)
	{
	  for(cidx=0; cidx<=ridx; ++cidx)
	    {
	      Amat[matsize2*jn + ridx*matsize+cidx] = Avec[count];
	      Amat[matsize2*jn + (ridx+spdim)*matsize+cidx+spdim] = Avec[count];
	      Amat[matsize2*jn + (ridx+spdim)*matsize+cidx] = Avec[count+1];
	      if (ridx!=cidx){
		Amat[matsize2*jn + (cidx+spdim)*matsize+ridx] = -Avec[count+1];
	      }
	      count = count+2;
	    }
	}

    }

}

/**************************************************************************/

/* a quick implementation of the matlab "find" function;
   "find" returns the indices of v that are strictly positive;
   v is not expected to be long; n is the length of v;
   pos = 1, first occurrence; pos > 1, first pos occurrences;
   pos = -1, last occurrence;
 */
static int find(const double *v, double scalev, double shiftval, int pos, const int n)
{
    int i,idx = -1;
    double val;
    for (i = 0; i < n; ++i)
    {
        val = v[i]*scalev+shiftval;
        if(val>0) {
            idx = i;
            if (pos ==1)
                break;
        }
    }
    if(idx>=0 || (idx < 0 && pos ==-1)){ /* idx = -1; */
        /* mpi_one_printf("find return case 1: %d\n", idx); */
        return idx;
    }
    else{ /* can't find anything positive, return the last index */
        /* mpi_one_printf("find return case 2: %d\n", n); */
        return n;
    }
}


static double  detectFluctuation(double *obj,int irun,int maxflucfreq,double usum,double utol)
{
  int i,j;
  double errs[maxflucfreq],errsum;
  int backshift = MIN2((irun+1)/2,maxflucfreq);

 if (backshift>0)
   {
     for (j=1; j<=backshift; ++j)
       {
	 errsum = 0.0;
	 for (i=0; i<j; ++i)
	   {
	     errs[i] = fabs((obj[irun-i] - obj[irun-i-j])/obj[irun-i]);
	     errsum += errs[i];
	   }
	  mpi_one_printf("errsum = %+1.6e, usum = %+1.6e\n",errsum,usum);
	 if(errsum/j <= utol)
	   {  mpi_one_printf("fluctuation in gap detected, with flucfreq = %d\n",j);
	     /* to break the while loop */
	     usum = utol/10;
	     break;
	   }
       }
   }
 return usum;
}


/*******************************************************************************/
number run_matgrid_optgap_lp_DCG(vector3_list kpoints,
			     integer band1, integer band2,
			     integer maxrun, number utol,
			     number low_tol, number upp_tol,
			     integer kk, integer DCG, char *title)
{

    int k, ntot, n, ngrids,*nl, *nu,nk,ridx,cidx,rband,cband,count,nltmpt,nutmpt,stride,maxspdim;
    double *u,*eigenvalues,*depsdu,lambda_l,lambda_u,usum,usum_new;
    material_grid *grids;
    double lowtol = (double) low_tol; /* determines the size of lower subspace */
    double upptol = (double) upp_tol; /* determines the size of upper subspace */
    double scale;
    int irun; /* ,maxrun; */
    double *gap, *obj;
    int maxflucfreq = 5;
    double_array *bvec, *cvec, *Alin;
    int *Nl, *Nu, NB, NC, blockL, blockH;
    scalar_complex *Altemp, *Autemp, *depsdu_sc, *u_sc, *y_sc;
    double *Ak;
    /* double *Altemp1, *Autemp1; */
    double *Akmat, *eigvals, *foo;
    int MSK_FLAG, DCG_FLAG = 0, maxDCG=3, DCGrun = 0;
    double tolShift = 1e-4;
    char prefix[256];
 /**************/
#ifndef HAVE_MOSEK
    CHECK(0, "mosek is required for material_grids_optgap\n");
#else
    MSKint32t i, j;
    MSKrescodee  r;
    MSKenv_t     env = NULL;
    MSKtask_t    task = NULL;

    MSKint32t NUMCON;
    MSKint32t NUMVAR;

    MSKint32t  *acols, asub[2];
    double  aval[2];
    double *y;

 /**************/

    nk = kpoints.num_items;

    reset_epsilon();
    grids = get_material_grids(geometry, &ngrids);
    ntot = material_grids_ntot(grids, ngrids);

    /* u = [ubar_1, ubar_2, ... , ubar_ntot, theta, lambda_l, lambda_u] */
    n = ntot + 3;
    mpi_one_printf("number of decision variable n = %d\n",n);
    u = (double *) malloc(sizeof(double) * n);

    /* storing the dimensions of the lower and upper subspaces for each k */
    nl = (int *) malloc(sizeof(int) * nk * 2);
    nu = nl + nk;
    /* storing the number of the lower and upper approximating vectors for each k */
    Nl = (int *) malloc(sizeof(int) * nk * 2);
    Nu = Nl + nk;

    gap = (double *) malloc(sizeof(double) * (maxrun));
    obj = (double *) malloc(sizeof(double) * (maxrun));

    eigenvalues = (double *) malloc(sizeof(double) * num_bands);
    blockL = band1*(band1+1)/2 * n;
    blockH = (num_bands-band1)*(num_bands-band1+1)/2 * n;
    Altemp = (scalar_complex *) calloc(blockL * nk, sizeof(scalar_complex));
    Autemp = (scalar_complex *) calloc(blockH * nk, sizeof(scalar_complex));
    Ak = (double *) malloc(sizeof(double) * MAX2(band1*(band1+1), (num_bands-band1)*(num_bands-band1+1)));
    maxspdim = MAX2(band1, num_bands-band1);
    Akmat = (double *) malloc(sizeof(double) * maxspdim*maxspdim*4);
    eigvals = (double *) calloc(maxspdim*2, sizeof(double));
    foo = (double *) calloc(maxspdim*maxspdim*4, sizeof(double));


    /* n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz; */

    depsdu = (double *) malloc(sizeof(double) * ntot*mdata->fft_output_size);
    depsdu_sc = (scalar_complex *) calloc(ntot*mdata->fft_output_size, sizeof(scalar_complex));
    u_sc = (scalar_complex *) calloc(ntot, sizeof(scalar_complex));

    deps_du(depsdu, 1.0, grids, ngrids);
    for(j=0; j<ntot*mdata->fft_output_size; ++j)
      depsdu_sc[j].re = depsdu[j];

    /* bvec = (double_array *)malloc(sizeof(double_array)); */
    /* cvec = (double_array *)malloc(sizeof(double_array)); */
    Alin = (double_array *)malloc(sizeof(double_array));

    acols = (MSKint32t *) malloc(sizeof(MSKint32t)*n);
    for (j=0; j<n; ++j)
      acols[j] = j;


    r = MSK_makeenv(&env,NULL);

    usum = ntot; irun = 0;
    while (irun < maxrun && usum >= utol)
      {
	/* set to 1 to allow optimization */
	MSK_FLAG = 1;
	DCGrun = 0;

	NB = 0; NC = 0;

	initArray(Alin,2*(kk+1));
        lambda_l = 0.0;
        lambda_u = 100.0;

        material_grids_get(u, grids, ngrids);
	for(j=0; j<ntot; ++j)
	  u_sc[j].re = u[j];

	for (k = 0; k < nk; ++k) {
            randomize_fields();
            solve_kpoint(kpoints.items[k]);

	    for (j = 0; j < num_bands; ++j)
	      eigenvalues[j]  = freqs.items[j]*freqs.items[j];

            /* lambda_l and lambda_u updated at each k, only used outside the "for" loop for k */
            lambda_l = MAX2(eigenvalues[band1-1], lambda_l);
            lambda_u = MIN2(eigenvalues[band2-1], lambda_u);

            /* find(const double *v, double scalev, double shiftval, int pos, const int n) */
            nltmpt = (band1-1) - find(eigenvalues,-1,eigenvalues[band1-1]*(1-lowtol),-1,num_bands);
            nutmpt = find(eigenvalues,1,-eigenvalues[band2-1]*(1+upptol),1,num_bands) - (band2-1);

            nl[k] = nltmpt>=1 ? nltmpt : nltmpt+1;//(band1>1 ? nltmpt+1 : nltmpt );
            nu[k] = nutmpt>=1 ? nutmpt : nutmpt+1;//(num_bands-band2 ? nutmpt+1 : nutmpt);

            mpi_one_printf("nl[%d] = %d, nu[%d] = %d\n",k,nl[k],k,nu[k]);

	    scale = 1;
            count = 0;
            stride = (nl[k]+1)*nl[k]/2;
            for (ridx = 0; ridx < nl[k]; ++ridx)
            {
                rband = band1-nl[k]+1+ridx;
                for (cidx = 0; cidx <= ridx; ++cidx)
                {
                    cband = band1-nl[k]+1+cidx;

                    material_grids_SPt_blas(Altemp+blockL*k+count, depsdu_sc, u_sc, -scale, ntot, rband, cband, stride);
                    CASSIGN_SCALAR(Altemp[blockL*k+count+(n-2)*stride],(rband==cband)? 1:0.0, 0.0);
                    CASSIGN_SCALAR(Altemp[blockL*k+count+(n-1)*stride], 0.0, 0.0);
                    ++count;
                }
            }

	    count = 0;
            stride = (nu[k]+1)*nu[k]/2;
            for (ridx = 0; ridx < nu[k]; ++ridx)
            {
                rband = band2+ridx;
                for (cidx = 0; cidx <= ridx; ++cidx)
                {
                    cband = band2+cidx;

		    material_grids_SPt_blas(Autemp+blockH*k+count, depsdu_sc, u_sc, scale, ntot, rband, cband, stride);
                    CASSIGN_SCALAR(Autemp[blockH*k+count+stride*(n-2)], 0.0, 0.0);
                    CASSIGN_SCALAR(Autemp[blockH*k+count+stride*(n-1)],(rband==cband)? -1:0.0, 0.0);
                    ++count;
                }
            }

	    bvec = sumCombination(kk, 2*nl[k]);
	    cvec = sumCombination(kk, 2*nu[k]);
	    Nl[k] = bvec->used/nl[k]/2; NB += Nl[k];
	    Nu[k] = cvec->used/nu[k]/2; NC += Nu[k];

	    if (Alin->size < Alin->used + n*(Nl[k]+Nu[k])){
	      Alin->size = Alin->used + n*(Nl[k]+Nu[k]);
	      Alin->data = (double *)realloc(Alin->data, Alin->size * sizeof(double));
	    }

	    SDP2LP(Alin, Altemp+blockL*k, bvec, nl[k], n);
	    freeArray(bvec);
	    SDP2LP(Alin, Autemp+blockH*k, cvec, nu[k], n);
	    freeArray(cvec);

	}
	gap[irun] = 2*(lambda_u-lambda_l)/(lambda_l+lambda_u);
	mpi_one_printf("Before minimization, %d, gap is %0.15g \n",irun+1, gap[irun]);

	real *Altemp1 = (real *) Altemp;
	real *Autemp1 = (real *) Autemp;

	/*************************************************/
	/* OPTIMIZATION */
	NUMCON = NB+NC+1+ntot;
	NUMVAR = n;
	if ( r==MSK_RES_OK )
	  {
	    r = MSK_maketask(env,NUMCON,0,&task);
	    MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
	    r = MSK_putintparam(task,MSK_IPAR_NUM_THREADS,NUM_THREADS);

	    /* append NUMCON empty constraints */
	    if ( r == MSK_RES_OK )
	      r = MSK_appendcons(task,NUMCON);
	    /* append NUMVAR variables initially fixed at zero */
	    if ( r == MSK_RES_OK )
	      r = MSK_appendvars(task,NUMVAR);

	    /* set linear term c_j in the objective */
	    if ( r ==MSK_RES_OK )
	      r = MSK_putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE);
	    if ( r ==MSK_RES_OK )
	      r = MSK_putcj(task,NUMVAR-2,-2.0);
	    if ( r ==MSK_RES_OK )
	      r = MSK_putcj(task,NUMVAR-1,2.0);

	    /* set bounds on variables */
	    for (j=0; j<NUMVAR && r==MSK_RES_OK; ++j)
	      r = MSK_putvarbound( task,j,MSK_BK_LO,0.0,MSK_INFINITY);

	    /* set the bounds and row i of constraint */
	    for (i=0; i<NB+NC && r==MSK_RES_OK; ++i)
	      {
		r = MSK_putconbound(task, i, MSK_BK_LO, 0.0, MSK_INFINITY);
		r = MSK_putarow(task, i, n, acols, Alin->data+i*n);
	      }

	    asub[1] = ntot; aval[0] = 1.0; aval[1] = -1.0;
	    for (i=0; i<ntot && r==MSK_RES_OK; ++i)
	      {
		r = MSK_putconbound(task, NB+NC+i, MSK_BK_UP, -MSK_INFINITY, 0.0);
		asub[0]  = i;
		r = MSK_putarow(task, NB+NC+i, 2, asub, aval);
	      }

	    if ( r ==MSK_RES_OK )
	      r = MSK_putconbound(task, NUMCON-1, MSK_BK_FX, 1.0, 1.0);
	    asub[0] = n-2; asub[1] = n-1; aval[0] = 1.0; aval[1] = 1.0;
	    if ( r ==MSK_RES_OK )
	      r = MSK_putarow(task, NUMCON-1, 2, asub, aval);

	  }


	while (MSK_FLAG == 1 || DCG_FLAG == 1)
	  {
	    if ( r==MSK_RES_OK )
	      {
		/* set to zero to prevent multiple optimization */
		MSK_FLAG = 0;
		DCG_FLAG = 0;

		MSKrescodee trmcode;

		mpi_one_printf("start msk optimitization\n");

		/* Run optimizer */
		/* r = MSK_optimizetrm(task,&trmcode); */


		if (r == MSK_RES_OK)
		  r = MSK_optimizetrm(task,&trmcode);

		mpi_one_printf("msk optimitization done\n");
		MSK_solutionsummary (task,MSK_STREAM_MSG);

		if ( r==MSK_RES_OK )
		  {
		    MSKsolstae solsta;
		    r = MSK_getsolsta (task,MSK_SOL_ITR,&solsta);

		    switch(solsta)
		      {
		      case MSK_SOL_STA_OPTIMAL:
		      case MSK_SOL_STA_NEAR_OPTIMAL:
			y   = (double*) MSK_calloctask(task,n,sizeof(MSKrealt));

			MSK_getxx(task, MSK_SOL_ITR, y);
			MSK_getprimalobj (task, MSK_SOL_ITR, obj+irun);

			usum = 0.0;
			usum_new = 0.0;
			for(j=0; j<ntot; ++j){
			  usum += fabs(y[j]/y[ntot]-u_sc[j].re);
			  usum_new += fabs(y[j]/y[ntot]-u[j]);
			  u[j] = (double) y[j]/y[ntot];
			}
			usum /= ntot;
			usum_new /= ntot;

			/* if no sufficient improvement, stop DCG */
			if (usum_new < tolShift)
			  DCGrun = maxDCG;

			mpi_one_printf("After maximization, %d, objective is %0.15g, change_in_u = %g\n",irun+1, obj[irun], usum);

			break;
		      case MSK_SOL_STA_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_PRIM_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
			mpi_one_printf("Primal or dual infeasibility certificate found.\n");
			break;

		      case MSK_SOL_STA_UNKNOWN:
			{
			char symname[MSK_MAX_STR_LEN];
			char desc[MSK_MAX_STR_LEN];

			/* The solutions status is unknown. The termination code
			   indicating why the optimizer terminated prematurely. */

			mpi_one_printf("The solution status is unknown.\n");
			if ( r!=MSK_RES_OK )
			  {
			    /* A system failure e.g. out of space. */
			    MSK_getcodedesc(r,symname,desc);
			     mpi_one_printf("  Response code: %s\n",symname);
			  }
			else
			  {
			    /* No system failure e.g. an iteration limit is reached.  */
			    MSK_getcodedesc(trmcode,symname,desc);
			     mpi_one_printf("  Termination code: %s\n",symname);
			  }
			break;
			}
		      default:
			mpi_one_printf("Other solution status.");
			break;
		      }
		  }
		else
		   mpi_one_printf("Error while optimizing.\n");
	      }

	    if(DCG > 0 && DCGrun < maxDCG)
	      {
		++DCGrun;
		double alpha = 1.0, beta = 0.0;
		NB = 0; NC = 0;

		int idx;
		for (k = 0; k<nk; ++k)
		  {
		    /* lower subspace */
		    cblas_dgemv (CblasRowMajor,CblasTrans, n, nl[k]*(nl[k]+1), alpha, Altemp1+k*blockL*2, nl[k]*(nl[k]+1), y, 1, beta, Ak, 1);
		    Ad2Ad(Akmat, Ak, nl[k], 1);
		    /* Az2Ad(Akmat, Ak,nl[k], 1); */
		    lapackglue_syev('V', 'L', 2*nl[k] , Akmat, 2*nl[k], eigvals, foo, 6*nl[k]);


		    /* forth argument pos = 1 --> first occurance */
		    idx = find(eigvals, 1.0, tolShift, 1, 2*nl[k]);

		    if (idx>0)
		      {
			DCG_FLAG = 1;
			bvec = (double_array *)malloc(sizeof(double_array));
			initArray(bvec, 2*nl[k]*idx);
			bvec->used = 2*nl[k]*idx;

			for (j=0; j<idx; ++j)
			  for (i=0; i<2*nl[k]; ++i)
			    bvec->data[2*nl[k]*j+i] = Akmat[2*nl[k]*i+j];

			NB += idx;
			if (Alin->size < Alin->used + n*idx){
			  Alin->size = Alin->used + n*idx;
			  Alin->data = (double *)realloc(Alin->data, Alin->size * sizeof(double));
			}
			SDP2LP(Alin, Altemp+k*blockL, bvec, nl[k], n);
			free(bvec);
		      } /* end of "if idx>0" */

		    /* upper subspace */
		    /* cblas_zgemv (CblasRowMajor,CblasTrans, n, nu[k]*(nu[k]+1)/2, &alpha, Autemp+k*blockH*2, nu[k]*(nu[k]+1), y_sc, 1, &beta, Ak, 1); */
		    /* Az2Ad(Akmat, Ak,nu[k], 1); */
		    cblas_dgemv (CblasRowMajor,CblasTrans, n, nu[k]*(nu[k]+1), alpha, Autemp1+k*blockH*2, nu[k]*(nu[k]+1), y, 1, beta, Ak, 1);
		    Ad2Ad(Akmat, Ak,nu[k], 1);
		    lapackglue_syev('V', 'L', 2*nu[k] , Akmat, 2*nu[k], eigvals, foo, 6*nu[k]);

		    /* forth argument pos = 1 --> first occurance */
		    idx = find(eigvals, 1.0, tolShift, 1, 2*nu[k]);

		    if (idx>0)
		      {
			DCG_FLAG = 1;
			cvec = (double_array *)malloc(sizeof(double_array));
			initArray(cvec, 2*nu[k]*idx);
			cvec->used = 2*nu[k]*idx;

			for (j=0; j<idx; ++j)
			  for (i=0; i<2*nu[k]; ++i)
			    cvec->data[2*nu[k]*j+i] = Akmat[2*nu[k]*i+j];

			NC += idx;
			if (Alin->size < Alin->used + n*idx){
			  Alin->size = Alin->used + n*idx;
			  Alin->data = (double *)realloc(Alin->data, Alin->size * sizeof(double));
			}
			SDP2LP(Alin, Autemp+k*blockH, cvec, nu[k], n);
			free(cvec);
		      } /* end of "if idx>0" */

		  }

	    	/*************************************************/
		/* OPTIMIZATION */
		if ( DCG_FLAG && r==MSK_RES_OK )
		  {
		    if (r == MSK_RES_OK)
		      r = MSK_getnumcon(task,&NUMCON);

		    /* append NUMCON empty constraints */
		    if ( r == MSK_RES_OK )
		      r = MSK_appendcons(task,NB+NC);

		    /* set the bounds and row i of constraint */
		    for (i=0; i<NB+NC && r==MSK_RES_OK; ++i)
		      {
			r = MSK_putconbound(task, NUMCON+i, MSK_BK_LO, 0.0, MSK_INFINITY);
			r = MSK_putarow(task, NUMCON+i, n, acols, Alin->data+(NUMCON-ntot-1+i)*n);
		      }
		    /* if (r == MSK_RES_OK) */
		      /* r = MSK_putintparam(task,MSK_IPAR_OPTIMIZER,MSK_OPTIMIZER_FREE_SIMPLEX); */
		  }
		NUMCON += NB+NC;
		mpi_one_printf("DCG run = %d\n", DCG);
	      } /* only if DCG is enabled */

	  } /* only if MSG_FLAG or DCG_FLAG */


	/* update u and epsilon */
	material_grids_set(u, grids, ngrids);
	reset_epsilon();

	/* detect flunctuation */
	usum = detectFluctuation(gap,irun,maxflucfreq,usum,utol);
	MSK_freetask(task,y);


	MSK_deletetask(&task);
	freeArray(Alin);


       	/* output epsilon file at each iteration */

       get_epsilon();
       snprintf(prefix, 256, "%s%04d-", title, irun+1);
       output_field_to_file(-1, prefix);
       strcat(prefix,"grid");
       save_material_grid(*grids, prefix);

	++irun;

      }

    /*************************************************/

    free(nl);
    free(Nl);
    free(u);
    free(u_sc);
    free(eigenvalues);
    free(depsdu);
    free(depsdu_sc);
    free(Altemp);
    free(Autemp);

    free(gap);
    free(obj);
    free(acols);
    free(Ak);
    free(Akmat);
    free(eigvals);
    free(foo);
    MSK_deleteenv(&env);


#endif
    return 1;
}

/* /\*******************************************************************************\/ */
/* number run_matgrid_optgap_lp(vector3_list kpoints, */
/* 			     integer band1, integer band2, */
/* 			     integer maxrun, number utol, */
/* 			     number low_tol, number upp_tol, */
/* 			     integer kk, char *title) */
/* { */
/*   int k, ntot, n, ngrids,*nl, *nu,nk,ridx,cidx,rband,cband,count,nltmpt,nutmpt,stride,irun; */
/*     double *u,*eigenvalues,*depsdu,lambda_l,lambda_u,usum,scale; */
/*     material_grid *grids; */
/*     double lowtol = (double) low_tol; /\* determines the size of lower subspace *\/ */
/*     double upptol = (double) upp_tol; /\* determines the size of upper subspace *\/ */
/*     double *gap, *obj; */
/*     int maxflucfreq = 5; */
/*     double_array *bvec, *cvec, *Alin; */
/*     int *Nl, *Nu, NB, NC; */
/*     scalar_complex *Altemp, *Autemp, *depsdu_sc, *u_sc; */
/*     char prefix[256]; */
/*  /\**************\/ */
/* #ifndef HAVE_MOSEK */
/*     CHECK(0, "mosek is required for material_grids_optgap\n"); */
/* #else */
/*     MSKint32t i, j; */
/*     MSKrescodee  r; */
/*     MSKenv_t     env = NULL; */
/*     MSKtask_t    task = NULL; */

/*     MSKint32t NUMCON; */
/*     MSKint32t NUMVAR; */

/*     MSKint32t  *acols, asub[2]; */
/*     double  aval[2]; */
/*     double *y; */

/*  /\**************\/ */

/*     nk = kpoints.num_items; */
/*     grids = get_material_grids(geometry, &ngrids); */
/*     ntot = material_grids_ntot(grids, ngrids); */

/*     /\* u = [ubar_1, ubar_2, ... , ubar_ntot, theta, lambda_l, lambda_u] *\/ */
/*     n = ntot + 3; */
/*     mpi_one_printf("number of decision variable n = %d\n",n); */
/*     u = (double *) malloc(sizeof(double) * n); */
/*     y = (double *) malloc(sizeof(double) * n); */


/*     /\* storing the dimensions of the lower and upper subspaces for each k *\/ */
/*     nl = (int *) malloc(sizeof(int) * nk * 2); */
/*     nu = nl + nk; */
/*     /\* storing the number of the lower and upper approximating vectors for each k *\/ */
/*     Nl = (int *) malloc(sizeof(int) * nk * 2); */
/*     Nu = Nl + nk; */

/*     gap = (double *) malloc(sizeof(double) * (maxrun)); */
/*     obj = (double *) malloc(sizeof(double) * (maxrun)); */

/*     eigenvalues = (double *) malloc(sizeof(double) * num_bands); */
/*     Altemp = (scalar_complex *) calloc(band1*(band1+1)/2*n, sizeof(scalar_complex)); */
/*     Autemp = (scalar_complex *) calloc((num_bands-band1)*(num_bands-band1+1)/2*n, sizeof(scalar_complex)); */


/*     /\* n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz; *\/ */

/*     depsdu = (double *) malloc(sizeof(double) * ntot*mdata->fft_output_size); */
/*     depsdu_sc = (scalar_complex *) calloc(ntot*mdata->fft_output_size, sizeof(scalar_complex)); */
/*     u_sc = (scalar_complex *) calloc(ntot, sizeof(scalar_complex)); */

/*     deps_du(depsdu, 1.0, grids, ngrids); */
/*     for(j=0; j<ntot*mdata->fft_output_size; j++) */
/*       depsdu_sc[j].re = depsdu[j]; */

/*     Alin = (double_array *)malloc(sizeof(double_array)); */
/*     initArray(Alin,2); */

/*     acols = (MSKint32t *) malloc(sizeof(MSKint32t)*n); */
/*     for (j=0; j<n; j++) */
/*       acols[j] = j; */

/*    /\* OPTIMIZATION *\/ */
/*     r = MSK_makeenv(&env,NULL); */
/*     NUMVAR = n; */

/*     usum = ntot; irun = 0; */
/*     while (irun < maxrun && usum >= utol) */
/*       { */

/* 	if ( r==MSK_RES_OK ) */
/* 	  { */
/* 	    r = MSK_maketask(env,0,NUMVAR,&task); */
/* 	    MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); */
/* 	  } */

/* 	/\* append NUMVAR variables initially fixed at zero *\/ */
/* 	if ( r == MSK_RES_OK ) */
/* 	  r = MSK_appendvars(task,NUMVAR); */
/* 	/\* set bounds on variables *\/ */
/* 	for (j=0; j<NUMVAR && r==MSK_RES_OK; ++j) */
/* 	  r = MSK_putvarbound( task,j,MSK_BK_LO,0.0,MSK_INFINITY); */



/* 	NB = 0; NC = 0; */

/*         lambda_l = 0.0; */
/*         lambda_u = 100.0; */

/*         material_grids_get(u, grids, ngrids); */
/* 	for(j=0; j<ntot; j++) */
/* 	  u_sc[j].re = u[j]; */

/* 	for (k = 0; k < nk; ++k) { */
/*             randomize_fields(); */
/*             solve_kpoint(kpoints.items[k]); */

/* 	    int rank;  */
/* 	    MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
/* 	    /\* printf("k = %d, r = %d, fft_output_size = %d, nx = %d, ny = %d, nz = %d\n", k, rank, mdata->fft_output_size, mdata->nx, mdata->ny, mdata->nz); *\/ */


/* 	    for (j = 0; j < num_bands; ++j) */
/* 	      eigenvalues[j]  = freqs.items[j]*freqs.items[j]; */

/*             /\* lambda_l and lambda_u updated at each k, only used outside the "for" loop for k *\/ */
/*             lambda_l = MAX2(eigenvalues[band1-1], lambda_l); */
/*             lambda_u = MIN2(eigenvalues[band2-1], lambda_u); */

/*             /\* find(const double *v, double scalev, double shiftval, int pos, const int n) *\/ */
/*             nltmpt = (band1-1) - find(eigenvalues,-1,eigenvalues[band1-1]*(1-lowtol),-1,num_bands); */
/*             nutmpt = find(eigenvalues,1,-eigenvalues[band2-1]*(1+upptol),1,num_bands) - (band2-1); */

/*             nl[k] = nltmpt>=1 ? nltmpt : nltmpt+1;//(band1>1 ? nltmpt+1 : nltmpt ); */
/*             nu[k] = nutmpt>=1 ? nutmpt : nutmpt+1;//(num_bands-band2 ? nutmpt+1 : nutmpt); */

/*             mpi_one_printf("nl[%d] = %d, nu[%d] = %d\n",k,nl[k],k,nu[k]); */

/* 	    scale = 1; */
/*             count = 0; */
/*             stride = (nl[k]+1)*nl[k]/2; */
/*             for (ridx = 0; ridx < nl[k]; ++ridx) */
/*             { */
/*                 rband = band1-nl[k]+1+ridx; */
/*                 for (cidx = 0; cidx <= ridx; ++cidx) */
/*                 { */
/*                     cband = band1-nl[k]+1+cidx; */

/*                     material_grids_SPt_blas(Altemp+count, depsdu_sc, u_sc, -scale, ntot, rband, cband, stride); */
/*                     CASSIGN_SCALAR(Altemp[count+(n-2)*stride],(rband==cband)? 1:0.0, 0.0); */
/*                     CASSIGN_SCALAR(Altemp[count+(n-1)*stride], 0.0, 0.0); */
/*                     ++count; */
/*                 } */
/*             } */

/* 	    count = 0; */
/*             stride = (nu[k]+1)*nu[k]/2; */
/*             for (ridx = 0; ridx < nu[k]; ++ridx) */
/*             { */
/*                 rband = band2+ridx; */
/*                 for (cidx = 0; cidx <= ridx; ++cidx) */
/*                 { */
/*                     cband = band2+cidx; */

/* 		    material_grids_SPt_blas(Autemp+count, depsdu_sc, u_sc, scale, ntot, rband, cband, stride); */
/*                     CASSIGN_SCALAR(Autemp[count+stride*(n-2)], 0.0, 0.0); */
/*                     CASSIGN_SCALAR(Autemp[count+stride*(n-1)],(rband==cband)? -1:0.0, 0.0); */
/*                     ++count; */
/*                 } */
/*             } */

/* 	    bvec = sumCombination(kk, 2*nl[k]); */
/* 	    cvec = sumCombination(kk, 2*nu[k]); */
/* 	    Nl[k] = bvec->used/nl[k]/2; */
/* 	    Nu[k] = cvec->used/nu[k]/2; */

/* 	    Alin->used = 0; */
/* 	    if (Alin->size < n*(Nl[k]+Nu[k])){ */
/* 	      Alin->size = n*(Nl[k]+Nu[k]); */
/* 	      Alin->data = (double *)realloc(Alin->data, Alin->size * sizeof(double)); */
/* 	    } */

/* 	    SDP2LP(Alin, Altemp, bvec, nl[k], n); */
/* 	    freeArray(bvec); */
/* 	    SDP2LP(Alin, Autemp, cvec, nu[k], n); */
/* 	    freeArray(cvec); */

/* 	    if ( r == MSK_RES_OK ) */
/* 	      r = MSK_appendcons(task,Nl[k]+Nu[k]); */
/* 	    /\* set row i of constraint *\/ */
/* 	    for (i=0; i<Nl[k]+Nu[k] && r==MSK_RES_OK; ++i) */
/* 		r = MSK_putarow(task, NB+NC+i, n, acols, Alin->data+i*n); */

/* 	    NB += Nl[k];NC += Nu[k]; */
/* 	} */
/* 	gap[irun] = 2*(lambda_u-lambda_l)/(lambda_l+lambda_u); */
/* 	mpi_one_printf("Before maximization, %d, gap is %0.15g \n",irun+1, gap[irun]); */




/* 	if ( r==MSK_RES_OK ) */
/* 	  { */


/* 	    /\* set linear term c_j in the objective *\/ */
/* 	    if ( r ==MSK_RES_OK ) */
/* 	      r = MSK_putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE); */
/* 	    if ( r ==MSK_RES_OK ) */
/* 	      r = MSK_putcj(task,NUMVAR-2,-2.0); */
/* 	    if ( r ==MSK_RES_OK ) */
/* 	      r = MSK_putcj(task,NUMVAR-1,2.0); */

/* 	    /\* set the bounds and row i of constraint *\/ */
/* 	    for (i=0; i<NB+NC && r==MSK_RES_OK; ++i) */
/* 	      { */
/* 	    	r = MSK_putconbound(task, i, MSK_BK_LO, 0.0, MSK_INFINITY); */
/* 	    	/\* r = MSK_putarow(task, i, n, acols, Alin->data+i*n); *\/ */
/* 	      } */

/* 	    /\* append 1+ntot empty constraints *\/ */
/* 	    if ( r == MSK_RES_OK ) */
/* 	      r = MSK_appendcons(task,1+ntot); */

/* 	    asub[1] = ntot; aval[0] = 1.0; aval[1] = -1.0; */
/* 	    for (i=0; i<ntot && r==MSK_RES_OK; ++i) */
/* 	      { */
/* 		r = MSK_putconbound(task, NB+NC+i, MSK_BK_UP, -MSK_INFINITY, 0.0); */
/* 		asub[0]  = i; */
/* 		r = MSK_putarow(task, NB+NC+i, 2, asub, aval); */
/* 	      } */

/* 	    if ( r ==MSK_RES_OK ) */
/* 	      r = MSK_putconbound(task, NB+NC+ntot, MSK_BK_FX, 1.0, 1.0); */
/* 	    asub[0] = n-2; asub[1] = n-1; aval[0] = 1.0; aval[1] = 1.0; */
/* 	    if ( r ==MSK_RES_OK ) */
/* 	      r = MSK_putarow(task, NB+NC+ntot, 2, asub, aval); */

/* 	  } */
/* 	mpi_one_printf("Starting msk\n"); */
/* 	if ( r==MSK_RES_OK ) */
/* 	  { */
/* 	     mpi_one_printf("MSK_RES_OK, Starting msk\n"); */


/* 	    r = MSK_putintparam(task,MSK_IPAR_NUM_THREADS,NUM_THREADS); */

/* 	    MSKrescodee trmcode; */
/* 	    r = MSK_optimizetrm(task,&trmcode); */

/* 	     mpi_one_printf("msk optimization done\n"); */
/* 	    MSK_solutionsummary (task,MSK_STREAM_MSG); */


/* 	    if ( r==MSK_RES_OK ) */
/* 	      { */

/* 		  mpi_one_printf("check msk solution\n"); */

/* 		MSKsolstae solsta; */
/* 		r = MSK_getsolsta (task,MSK_SOL_ITR,&solsta); */

/* 		switch(solsta) */
/* 		  { */
/* 		  case MSK_SOL_STA_OPTIMAL: */
/* 		  case MSK_SOL_STA_NEAR_OPTIMAL: */
/* 		    /\* y   = (double*) MSK_calloctask(task,n,sizeof(MSKrealt)); *\/ */

/* 		    MSK_getxx(task, MSK_SOL_ITR, y); */
/* 		    MSK_getprimalobj (task, MSK_SOL_ITR, obj+irun); */

/* 		    usum = 0.0; */
/* 		    for(j=0; j<ntot; j++){ */
/* 		      usum += fabs(y[j]/y[ntot]-u[j]); */
/* 		      u[j] = (double) y[j]/y[ntot]; */
/* 		    } */
/* 		    usum /= ntot; */
/* 		     mpi_one_printf("After maximization, %d, objective is %0.15g, change_in_u = %g\n",irun+1, obj[irun], usum); */

/* 		    /\* detect flunctuation *\/ */
/* 		    usum = detectFluctuation(gap,irun,maxflucfreq,usum,utol); */

/* 		    /\* update u and epsilon *\/ */
/* 		    material_grids_set(u, grids, ngrids); */
/* 		    reset_epsilon(); */

/* 		    /\* MSK_freetask(task,y); *\/ */

/* 		    break; */
/* 		  case MSK_SOL_STA_DUAL_INFEAS_CER: */
/* 		  case MSK_SOL_STA_PRIM_INFEAS_CER: */
/* 		  case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: */
/* 		  case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER: */
/* 		     mpi_one_printf("Primal or dual infeasibility certificate found.\n"); */
/* 		    break; */

/* 		  case MSK_SOL_STA_UNKNOWN: */
/* 		    { */
/* 		      char symname[MSK_MAX_STR_LEN]; */
/* 		      char desc[MSK_MAX_STR_LEN]; */

/* 		      /\* The solutions status is unknown. The termination code */
/* 			 indicating why the optimizer terminated prematurely. *\/ */

/* 		       mpi_one_printf("The solution status is unknown.\n"); */
/* 		      if ( r!=MSK_RES_OK ) */
/* 			{ */
/* 			  /\* A system failure e.g. out of space. *\/ */
/* 			  MSK_getcodedesc(r,symname,desc); */
/* 			   mpi_one_printf("  Response code: %s\n",symname); */
/* 			} */
/* 		      else */
/* 			{ */
/* 			  /\* No system failure e.g. an iteration limit is reached.  *\/ */
/* 			  MSK_getcodedesc(trmcode,symname,desc); */
/* 			   mpi_one_printf("  Termination code: %s\n",symname); */
/* 			} */
/* 		      break; */
/* 		    } */
/* 		  default: */
/* 		     mpi_one_printf("Other solution status."); */
/* 		    break; */
/* 		  } */
/* 	      } */
/* 	    else */
/* 	       mpi_one_printf("Error while optimizing.\n"); */
/* 	  } */

/*        	/\* output epsilon file at each iteration *\/ */
/* 	get_epsilon(); */
/* 	/\* snprintf(prefix, 256, "%s%d-%04d-", title, band1, irun); *\/ */
/* 	snprintf(prefix, 256, "%s%04d-", title, irun+1); */
/* 	output_field_to_file(-1, prefix); */
/* 	strcat(prefix,"grid"); */
/* 	save_material_grid(*grids, prefix); */

/* 	MSK_deletetask(&task); */

/* 	++irun; */
/*       } */

/*     /\*************************************************\/ */

/*     free(nl); */
/*     free(Nl); */
/*     free(u); */
/*     free(u_sc); */
/*     free(y); */
/*     free(eigenvalues); */
/*     free(depsdu); */
/*     free(depsdu_sc); */
/*     free(Altemp); */
/*     free(Autemp); */
/*     freeArray(Alin); */
/*     free(gap); */
/*     free(obj); */
/*     free(acols); */
/*     MSK_deleteenv(&env); */


/* #endif */
/*     return 1; */
/* } */


/*******************************************************************************/
number run_matgrid_optgap_lp(vector3_list kpoints,
			     integer band1, integer band2,
			     integer maxrun, number utol,
			     number low_tol, number upp_tol,
			     integer kk, char *title)
{
  int k, ntot, n, ngrids,*nl, *nu,nk,ridx,cidx,rband,cband,count,nltmpt,nutmpt,stride,irun;
    double *u,*eigenvalues,*depsdu,lambda_l,lambda_u,usum,scale;
    material_grid *grids;
    double lowtol = (double) low_tol; /* determines the size of lower subspace */
    double upptol = (double) upp_tol; /* determines the size of upper subspace */
    double *gap, *obj;
    int maxflucfreq = 5;
    double_array *bvec, *cvec, *Alin;
    int *Nl, *Nu, NB, NC;
    scalar_complex *Altemp, *Autemp, *depsdu_sc, *u_sc;
    char prefix[256];
    char solstastr[30];
 /**************/
#ifndef HAVE_MOSEK
    CHECK(0, "mosek is required for material_grids_optgap\n");
#else
    MSKint32t i, j;
    MSKrescodee  r;
    MSKenv_t     env = NULL;
    MSKtask_t    task = NULL;

    MSKint32t NUMCON;
    MSKint32t NUMVAR;

    MSKint32t  *acols, asub[2];
    double  aval[2];
    double *y;

 /**************/
    int numtasks = NUM_THREADS;

#ifdef HAVE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    /* int optRank = numtasks-1; */
#endif

    nk = kpoints.num_items;
    grids = get_material_grids(geometry, &ngrids);
    /* synchronize_material_grid(grids); */

    ntot = material_grids_ntot(grids, ngrids);

    /* u = [ubar_1, ubar_2, ... , ubar_ntot, theta, lambda_l, lambda_u] */
    n = ntot + 3;
    mpi_one_printf("number of decision variable n = %d\n",n);
    u = (double *) malloc(sizeof(double) * n);
    y = (double *) malloc(sizeof(double) * n);


    /* storing the dimensions of the lower and upper subspaces for each k */
    nl = (int *) malloc(sizeof(int) * nk * 2);
    nu = nl + nk;
    /* storing the number of the lower and upper approximating vectors for each k */
    Nl = (int *) malloc(sizeof(int) * nk * 2);
    Nu = Nl + nk;

    gap = (double *) malloc(sizeof(double) * (maxrun));
    obj = (double *) malloc(sizeof(double) * (maxrun));

    eigenvalues = (double *) malloc(sizeof(double) * num_bands);
    Altemp = (scalar_complex *) calloc(band1*(band1+1)/2*n, sizeof(scalar_complex));
    Autemp = (scalar_complex *) calloc((num_bands-band1)*(num_bands-band1+1)/2*n, sizeof(scalar_complex));

    depsdu = (double *) malloc(sizeof(double) * ntot*mdata->fft_output_size);
    depsdu_sc = (scalar_complex *) calloc(ntot*mdata->fft_output_size, sizeof(scalar_complex));
    u_sc = (scalar_complex *) calloc(ntot, sizeof(scalar_complex));

    deps_du(depsdu, 1.0, grids, ngrids);
    for(j=0; j<ntot*mdata->fft_output_size; j++)
      depsdu_sc[j].re = depsdu[j];

    Alin = (double_array *)malloc(sizeof(double_array));
    initArray(Alin,2);

    acols = (MSKint32t *) malloc(sizeof(MSKint32t)*n);
    for (j=0; j<n; j++)
      acols[j] = j;

   /* OPTIMIZATION */
    r = MSK_makeenv(&env,NULL);
    NUMVAR = n;

    usum = ntot; irun = 0;
    while (irun < maxrun && usum >= utol)
      {
	r = MSK_maketask(env,0,NUMVAR,&task);
	MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

#ifdef HAVE_MPI
	if (mpi_is_master())
#endif
	if (r==MSK_RES_OK)
	  {
	      r = MSK_appendvars(task,NUMVAR);
	    for (j=0; j<NUMVAR && r==MSK_RES_OK; ++j)
	      r = MSK_putvarbound( task,j,MSK_BK_LO,0.0,MSK_INFINITY);
	  }

	NB = 0; NC = 0;
	Alin->used = 0;

        lambda_l = 0.0;
        lambda_u = 100.0;

        material_grids_get(u, grids, ngrids);
	for(j=0; j<ntot; j++)
	  u_sc[j].re = u[j];

	for (k = 0; k < nk; ++k) {
            randomize_fields();
            solve_kpoint(kpoints.items[k]);

	    for (j = 0; j < num_bands; ++j)
	      eigenvalues[j]  = freqs.items[j]*freqs.items[j];

            /* lambda_l and lambda_u updated at each k, only used outside the "for" loop for k */
            lambda_l = MAX2(eigenvalues[band1-1], lambda_l);
            lambda_u = MIN2(eigenvalues[band2-1], lambda_u);

            /* find(const double *v, double scalev, double shiftval, int pos, const int n) */
            nltmpt = (band1-1) - find(eigenvalues,-1,eigenvalues[band1-1]*(1-lowtol),-1,num_bands);
            nutmpt = find(eigenvalues,1,-eigenvalues[band2-1]*(1+upptol),1,num_bands) - (band2-1);

            nl[k] = nltmpt>=1 ? nltmpt : nltmpt+1;//(band1>1 ? nltmpt+1 : nltmpt );
            nu[k] = nutmpt>=1 ? nutmpt : nutmpt+1;//(num_bands-band2 ? nutmpt+1 : nutmpt);

            mpi_one_printf("nl[%d] = %d, nu[%d] = %d\n",k,nl[k],k,nu[k]);

	    scale = 1;

	    compute_subspace(band1, 'l', nl[k], -scale, Altemp, depsdu_sc, u_sc, ntot);
	    compute_subspace(band2, 'u', nu[k], scale, Autemp, depsdu_sc, u_sc, ntot);

	    /* Alin->used = 0; */
	    Nl[k] = SDP2LP_v1(Alin, Altemp, nl[k], n, kk);
	    Nu[k] = SDP2LP_v1(Alin, Autemp, nu[k], n, kk);

	    NB += Nl[k];NC += Nu[k];

	}
	gap[irun] = 2*(lambda_u-lambda_l)/(lambda_l+lambda_u);
	mpi_one_printf("Before maximization %d, NB = %d, NC = %d, gap is %0.15g \n",irun+1, NB, NC, gap[irun]);


	/* MPI_Bcast(Alin->data, (NB+NC)*n, MPI_DOUBLE, 0, MPI_COMM_WORLD); */

	/*************************************************/
#ifdef HAVE_MPI
	printf("About to check rank, rank = %d, status = %d\n", rank, r);
	if (mpi_is_master())
	  {
#endif

	/* OPTIMIZATION */
	if ( r==MSK_RES_OK )
	  {
	    r = MSK_appendcons(task,NB+NC);
	    for (i=0; i<NB+NC && r==MSK_RES_OK; ++i)
	      {
		r = MSK_putarow(task, i, n, acols, Alin->data+i*n);
		r = MSK_putconbound(task, i, MSK_BK_LO, 0.0, MSK_INFINITY);
	      }

	    /* set linear term c_j in the objective */
	    if ( r ==MSK_RES_OK )
	      r = MSK_putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE);
	    if ( r ==MSK_RES_OK )
	      r = MSK_putcj(task,NUMVAR-2,-2.0);
	    if ( r ==MSK_RES_OK )
	      r = MSK_putcj(task,NUMVAR-1,2.0);


	    /* append 1+ntot empty constraints */
	    if ( r == MSK_RES_OK )
	      r = MSK_appendcons(task,1+ntot);

	    asub[1] = ntot; aval[0] = 1.0; aval[1] = -1.0;
	    for (i=0; i<ntot && r==MSK_RES_OK; ++i)
	      {
		r = MSK_putconbound(task, NB+NC+i, MSK_BK_UP, -MSK_INFINITY, 0.0);
		asub[0]  = i;
		r = MSK_putarow(task, NB+NC+i, 2, asub, aval);
	      }

	    if ( r ==MSK_RES_OK )
	      r = MSK_putconbound(task, NB+NC+ntot, MSK_BK_FX, 1.0, 1.0);
	    asub[0] = n-2; asub[1] = n-1; aval[0] = 1.0; aval[1] = 1.0;
	    if ( r ==MSK_RES_OK )
	      r = MSK_putarow(task, NB+NC+ntot, 2, asub, aval);

	    printf("Starting msk, status = %d\n", r);

	    printf("MSK_RES_OK, Starting msk\n");
	    r = MSK_putintparam(task,MSK_IPAR_NUM_THREADS, numtasks);

	    MSKrescodee trmcode;
	    r = MSK_optimizetrm(task,&trmcode);

	    printf("msk optimization done\n");
	    MSK_solutionsummary (task,MSK_STREAM_MSG);

	    if ( r==MSK_RES_OK )
	      {

		MSKsolstae solsta;
		r = MSK_getsolsta (task,MSK_SOL_BAS,&solsta);

		switch(solsta)
		  {
		  case MSK_SOL_STA_OPTIMAL:
		  case MSK_SOL_STA_NEAR_OPTIMAL:
		    MSK_getxx(task, MSK_SOL_BAS, y);
		    MSK_getprimalobj (task, MSK_SOL_BAS, obj+irun);

		    usum = 0.0;
		    for(j=0; j<ntot; j++){
		      usum += fabs(y[j]/y[ntot]-u[j]);
		      u[j] = (double) y[j]/y[ntot];
		    }
		    usum /= ntot;
		    printf("After maximization, %d, objective is %0.15g, change_in_u = %g\n",irun+1, obj[irun], usum);

		    /* detect flunctuation */
		    usum = detectFluctuation(gap,irun,maxflucfreq,usum,utol);

		    break;
		  case MSK_SOL_STA_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_PRIM_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
		      case MSK_SOL_STA_UNKNOWN:
		      default:
			print_solsta(solsta, solstastr);
			printf("%s\n", solstastr);
			break;
		  } /* end of switch */
	      }
	    else
	      printf("Error while optimizing.\n");
	  }
#ifdef HAVE_MPI
	  }
	mpi_one_printf("broadcasting usum\n");
	MPI_Bcast(&usum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	/* update u and epsilon */
	material_grids_set(u, grids, ngrids);
	synchronize_material_grid(grids);
	reset_epsilon();

       	/* output epsilon file at each iteration */
	get_epsilon();
	/* snprintf(prefix, 256, "%s%d-%04d-", title, band1, irun); */
	snprintf(prefix, 256, "%s%04d-", title, irun+1);
	output_field_to_file(-1, prefix);
	strcat(prefix,"grid");
	save_material_grid(*grids, prefix);

	MSK_deletetask(&task);

	++irun;
      }

    /*************************************************/

    free(nl);
    free(Nl);
    free(u);
    free(u_sc);
    free(y);
    free(eigenvalues);
    free(depsdu);
    free(depsdu_sc);
    free(Altemp);
    free(Autemp);
    freeArray(Alin);
    free(gap);
    free(obj);
    free(acols);
    MSK_deleteenv(&env);


#endif
    return 1;
}


int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

/*******************************************************************************/
number run_matgrid_optgap_lp_fa(vector3_list kpoints,
			     integer band1, integer band2,
			     integer maxrun, number utol,
			     number sp_tol, integer kk,
			     number delta, char *title)
{
  int k, ntot, n, ngrids,nl, nu,nk,ridx,cidx,rband,cband,count,nltmpt,nutmpt,stride,irun;
    double *u,*eigenvalues,*depsdu,lambda_l,lambda_u,usum,scale;
    material_grid *grids;
    double lowtol = (double) sp_tol; /* determines the size of lower subspace */
    double upptol = (double) sp_tol; /* determines the size of upper subspace */
    double *gap, *obj, obj_ij;

    int maxflucfreq = 5;
    double_array *bvec, *cvec, *Blin, *Clin;
    int Nl, Nu, NB, NC, ccount;
    scalar_complex *Altemp, *Autemp, *depsdu_sc, *u_sc;
    char prefix[256];
    char solstastr[30];

    struct timeval tbegin, tend, tdiff;

 /**************/
#ifndef HAVE_MOSEK
    CHECK(0, "mosek is required for material_grids_optgap\n");
#else
    MSKint32t i, j, ii, jj;
    MSKrescodee  r;
    MSKenv_t     env = NULL;
    MSKtask_t    task = NULL;
    MSKtask_t    taskij = NULL;
    MSKrescodee trmcode;

    MSKint32t NUMCON;
    MSKint32t NUMVAR;

    MSKint32t *ayidx_ij;
    double *y, *xx,*delf_ij, lhs;
    double min_obj, knorm;
 /**************/
    nk = kpoints.num_items;

    reset_epsilon();
    grids = get_material_grids(geometry, &ngrids);
    ntot = material_grids_ntot(grids, ngrids);
    mpi_one_printf("delta = %f, delta*ntot = %f\n", delta, delta*ntot);

    /* u = [ubar_1, ubar_2, ... , ubar_ntot, theta] */
    n = ntot + 1;
    mpi_one_printf("number of decision variable n = %d\n",n);
    u = (double *) malloc(sizeof(double) * n);

    gap = (double *) malloc(sizeof(double) * (maxrun));
    obj = (double *) malloc(sizeof(double) * (maxrun));

    eigenvalues = (double *) malloc(sizeof(double) * num_bands);
    Altemp = (scalar_complex *) calloc(band1*(band1+1)/2*n, sizeof(scalar_complex));
    Autemp = (scalar_complex *) calloc((num_bands-band1)*(num_bands-band1+1)/2*n, sizeof(scalar_complex));

    depsdu = (double *) malloc(sizeof(double) * ntot*mdata->fft_output_size);
    depsdu_sc = (scalar_complex *) calloc(ntot*mdata->fft_output_size, sizeof(scalar_complex));
    u_sc = (scalar_complex *) calloc(ntot, sizeof(scalar_complex));

    deps_du(depsdu, 1.0, grids, ngrids);
    for(j=0; j<ntot*mdata->fft_output_size; ++j)
      depsdu_sc[j].re = depsdu[j];
    free(depsdu);

    /* decision variables are [ybar, qbar, theta] */
    int numconij = 3*ntot+2;
    int numvarij = 2*ntot+1;
    ayidx_ij = (MSKint32t *) calloc(n, sizeof(MSKint32t));
    for (i=0; i<n; ++i)
      ayidx_ij[i] = i;
    xx = (double *) calloc(numvarij, sizeof(double));
    y = (double *) calloc(numconij, sizeof(double));
    delf_ij = (double *) calloc(ntot,sizeof(double));


    Blin = (double_array *)malloc(sizeof(double_array));
    initArray(Blin,2);

    Clin = (double_array *)malloc(sizeof(double_array));
    initArray(Clin,2);

   /* OPTIMIZATION */
    r = MSK_makeenv(&env,NULL);
    if (r == MSK_RES_OK)
      r = MSK_initenv(env);

    NUMVAR = ntot+1;

    usum = ntot; irun = 0;
    while (irun < maxrun && usum >= utol)
      {
	NB = 0; NC = 0;
	Blin->used = 0;
	Clin->used = 0;

        lambda_l = 0.0;
        lambda_u = 100.0;
	min_obj = 100;

        material_grids_get(u, grids, ngrids);
	for(j=0; j<ntot; ++j)
	  u_sc[j].re = u[j];

	for (k = 0; k < nk; ++k) {
            randomize_fields();
            solve_kpoint(kpoints.items[k]);

	    for (j = 0; j < num_bands; ++j)
	      eigenvalues[j]  = freqs.items[j]*freqs.items[j];

            /* lambda_l and lambda_u updated at each k, only used outside the "for" loop for k */
            lambda_l = MAX2(eigenvalues[band1-1], lambda_l);
            lambda_u = MIN2(eigenvalues[band2-1], lambda_u);

            /* find(const double *v, double scalev, double shiftval, int pos, const int n) */
            nltmpt = (band1-1) - find(eigenvalues,-1,eigenvalues[band1-1]*(1-lowtol),-1,num_bands);
            nutmpt = find(eigenvalues,1,-eigenvalues[band2-1]*(1+upptol),1,num_bands) - (band2-1);

            nl = nltmpt>=1 ? nltmpt : nltmpt+1;
            nu = nutmpt>=1 ? nutmpt : nutmpt+1;

            mpi_one_printf("nl[%d] = %d, nu[%d] = %d\n",k,nl,k,nu);

	    scale = 1.0;
	    knorm = kpoints.items[k].x*kpoints.items[k].x + kpoints.items[k].y*kpoints.items[k].y + kpoints.items[k].z*kpoints.items[k].z;
	    Nl = 0;
	    if (knorm>1e-3) {
	      compute_subspace_fa(band1, 'l', nl, scale, Altemp, depsdu_sc, u_sc,ntot);
	      Nl = SDP2LP_v1(Blin, Altemp, nl, n, kk);
	      NB += Nl;
	    }

	    compute_subspace_fa(band2, 'u', nu, scale, Autemp, depsdu_sc, u_sc,ntot);
	    Nu = SDP2LP_v1(Clin, Autemp, nu, n, kk);
	    NC += Nu;

	    mpi_one_printf("Nl[%d] = %d, Nu[%d] = %d\n",k,Nl,k,Nu);
	}
	gap[irun] = 2*(lambda_u-lambda_l)/(lambda_l+lambda_u);
	mpi_one_printf("Before maximization, %d, gap is %0.15g \n",irun+1, gap[irun]);

	/*************************************************/
	/* Outer Optimiztaion Task */
	if ( r== MSK_RES_OK )
	    r = MSK_maketask(env,NB*NC,NUMVAR,&task);
	if ( r == MSK_RES_OK )
	    MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
	/* append NUMVAR variables */
	if ( r == MSK_RES_OK )
	  r = MSK_appendvars(task,NUMVAR);
	/* set bounds on variables */
	for (j=0; j<ntot && r==MSK_RES_OK; ++j)
	  r = MSK_putvarbound(task,j,MSK_BK_RA,0.0,1.0);
	if ( r == MSK_RES_OK )
	  r = MSK_putvarbound(task,ntot,MSK_BK_FR,-MSK_INFINITY,+MSK_INFINITY);
	/* set objectives */
	if ( r == MSK_RES_OK )
	  r = MSK_putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE);
	/* decision variable of the outer task */
	/*   [u_1, u_2, ..., u_ntot, t] */
	if ( r == MSK_RES_OK )
	  r = MSK_putcj(task,ntot,1);

	/*************************************************/
	/* Inner tasks: (i,j) subproblem */

	/* create a subtask with numvarij variables, and numconij constraints */
	r = MSK_maketask(env,numconij,numvarij,&taskij);
	/* MSK_linkfunctotaskstream(taskij,MSK_STREAM_LOG,NULL,printstr); */
	/* append numvarij variables */
	if ( r == MSK_RES_OK )
	  r = MSK_appendvars(taskij,numvarij);
	/* set lower bounds on variables, all >= 0 */
	for (i=0; i<numvarij && r==MSK_RES_OK; ++i)
	  r = MSK_putvarbound(taskij,i,MSK_BK_LO,0.0,MSK_INFINITY);

	if ( r == MSK_RES_OK )
	  r = MSK_appendcons(taskij,numconij);
	/* set row i of constraint */
	for (i=0; i<ntot && r==MSK_RES_OK; ++i)
	  {
	    /* ybar - x \theta <= qbar */
	    r = MSK_putaij(taskij, i, i, 1.0);
	    r = MSK_putaij(taskij, i, ntot+i, -1.0);
	    r = MSK_putaij(taskij, i, 2*ntot, -u[i]);
	    r = MSK_putconbound(taskij, i, MSK_BK_UP, -MSK_INFINITY, 0.0);

	    /*  ybar - x \theta >= - qbar */
	    r = MSK_putaij(taskij, ntot+i, i, -1.0);
	    r = MSK_putaij(taskij, ntot+i, ntot+i, -1.0);
	    r = MSK_putaij(taskij, ntot+i, 2*ntot, +u[i]);
	    r = MSK_putconbound(taskij, ntot+i, MSK_BK_UP, -MSK_INFINITY, 0);

	    /* ybar - e \theta <= 0 */
	    r = MSK_putaij(taskij, 2*ntot+i, i, 1.0);
	    r = MSK_putaij(taskij, 2*ntot+i, 2*ntot, -1.0);
	    r = MSK_putconbound(taskij, 2*ntot+i, MSK_BK_UP, -MSK_INFINITY, 0.0);
	  }

	/* e'x qbar - delta x ntot x \theta < = 0 */
	for (i=0; i<ntot && r==MSK_RES_OK; ++i)
	    r = MSK_putaij(taskij, 3*ntot, ntot+i, 1.0);
	r = MSK_putaij(taskij, 3*ntot, 2*ntot, -delta*ntot);

	r = MSK_putconbound(taskij, 3*ntot, MSK_BK_UP, -MSK_INFINITY, 0.0);
	r = MSK_putconbound(taskij, 3*ntot+1, MSK_BK_FX, 1, 1);

	if ( r ==MSK_RES_OK )
	  r = MSK_putobjsense(taskij,MSK_OBJECTIVE_SENSE_MINIMIZE);

	if (r == MSK_RES_OK) {
	  r = MSK_putintparam(taskij,MSK_IPAR_NUM_THREADS,1);
	  r = MSK_putintparam(taskij,MSK_IPAR_OPTIMIZER,MSK_OPTIMIZER_FREE_SIMPLEX);
	}

	/* append emptry constraint */
	if ( r == MSK_RES_OK )
	  r = MSK_appendcons(task,NB*NC);

	ii = 0; jj = 0;
	mpi_one_printf("NB = %d, NC = %d \n", NB, NC);

	for(ii=0; ii<NB; ++ii)
	  {
	    /* /\* time(&tbegin); *\/ */
	    /* gettimeofday(&tbegin, NULL); */
	    mpi_one_printf("msk subproblem (%d, %d)\n", ii, jj);

	    for(jj=0; jj<NC; ++jj)
	      {
		ccount =  ii*NC+jj;

		for(j=0; j<ntot && r ==MSK_RES_OK; ++j)
		  r = MSK_putaij(taskij, 3*ntot+1, j, (Clin->data[jj*n+j] + Blin->data[ii*n+j]));
		r = MSK_putaij(taskij, 3*ntot+1, 2*ntot, (Clin->data[jj*n+ntot] + Blin->data[ii*n+ntot]));

		for(j = 0; j<ntot && r ==MSK_RES_OK; ++j)
		  r = MSK_putcj(taskij, j, 2*(Clin->data[jj*n+j] - Blin->data[ii*n+j]));
		r = MSK_putcj(taskij, 2*ntot, 2*(Clin->data[jj*n+ntot] - Blin->data[ii*n+ntot]));

		/* ==== */

		r = MSK_optimizetrm(taskij,&trmcode);
		/* MSK_solutionsummary (taskij,MSK_STREAM_MSG); */
		/* ==== */

		if ( r==MSK_RES_OK )
		  {

		    MSKsolstae solsta;
		    r = MSK_getsolsta (taskij,MSK_SOL_BAS,&solsta);

		    switch(solsta)
		      {
		      case MSK_SOL_STA_OPTIMAL:
		      case MSK_SOL_STA_NEAR_OPTIMAL:
			MSK_getxx(taskij, MSK_SOL_BAS, xx);
			MSK_getyslice(taskij, MSK_SOL_BAS, 0, 2*ntot, y);
			MSK_getprimalobj (taskij, MSK_SOL_BAS, &obj_ij);

			for(j = 0; j<ntot && r ==MSK_RES_OK; ++j)
			  delf_ij[j] = (y[j] - y[ntot+j])*xx[2*ntot];
			delf_ij[ntot] = -1.0;

			if(r == MSK_RES_OK)
			  r = MSK_putarow(task, ccount, n, ayidx_ij, delf_ij);

			/* lhs constant (lowerbound) = delf_ij'*uhat - ~f_ij(xhat)  */
			lhs = cblas_ddot(ntot,delf_ij,1,u,1);
			lhs -= obj_ij;

			min_obj = MIN2(min_obj,obj_ij);

			if(r == MSK_RES_OK)
			  r = MSK_putconbound(task, ccount, MSK_BK_LO, lhs, MSK_INFINITY);

			break;
		      default:
			print_solsta(solsta, solstastr);
			mpi_one_printf("%s\n", solstastr);
			break;
		      }
		  }
		else
		  mpi_one_printf("Error while optimizing.\n");
	      }
	    /* /\* time(&tend); *\/ */
	    /* /\* mpi_one_printf("The loop used %f seconds.\n", difftime(tend, tbegin)); *\/ */
	    /* gettimeofday(&tend, NULL); */
	    /* timeval_subtract(&tdiff, &tend, &tbegin); */
	    /* mpi_one_printf("The loop used %ld.%06ld seconds\n", tdiff.tv_sec, tdiff.tv_usec);	     */

	    mpi_one_printf("min_obj = %f\n", min_obj);
	  }

	/*************************************************/
	mpi_one_printf("number of constraints in outer task = %d\n", NB*NC);

	/* OPTIMIZATION of Outer Task*/
	mpi_one_printf("Starting outer task\n");
	if ( r==MSK_RES_OK )
	  {
	    r = MSK_putintparam(task,MSK_IPAR_OPTIMIZER,MSK_OPTIMIZER_FREE_SIMPLEX);
	    r = MSK_putintparam(task,MSK_IPAR_NUM_THREADS,1);

	    /* MSKrescodee trmcode; */
	    r = MSK_optimizetrm(task,&trmcode);

	    /* MSK_solutionsummary (task,MSK_STREAM_MSG); */

	    if ( r==MSK_RES_OK )
	      {

		  mpi_one_printf("check msk solution\n");

		/* MSKsolstae solsta; */
		  MSKsolstae solsta;
		r = MSK_getsolsta (task,MSK_SOL_BAS,&solsta);

		switch(solsta)
		  {
		  case MSK_SOL_STA_OPTIMAL:
		  case MSK_SOL_STA_NEAR_OPTIMAL:
		    MSK_getxx(task, MSK_SOL_BAS, xx);
		    mpi_one_printf("t = %g\n", xx[ntot]);
		    usum = 0.0;
		    for(j=0; j<ntot; ++j){
		      usum += fabs(xx[j]-u[j]);

		      u[j] = xx[j];
		    }
		    usum /= ntot;

		    MSK_getprimalobj (task, MSK_SOL_BAS, obj+irun);
		    mpi_one_printf("After maximization, %d, objective is %0.15g, change_in_u = %g\n",irun+1, obj[irun], usum);

		    /* detect flunctuation */
		    usum = detectFluctuation(gap,irun,maxflucfreq,usum,utol);
		    /* update u and epsilon */
		    material_grids_set(u, grids, ngrids);
		    reset_epsilon();

		    break;
		  default:
		     print_solsta(solsta, solstastr);
		     mpi_one_printf("%s\n", solstastr);
		    break;
		  }
	      }
	    else
	       mpi_one_printf("Error while optimizing.\n");
	  }

       	/* output epsilon file at each iteration */
	get_epsilon();
	snprintf(prefix, 256, "%s%04d-", title, irun+1);
	output_field_to_file(-1, prefix);
	strcat(prefix,"grid");
	save_material_grid(*grids, prefix);

	mpi_one_printf("about to delete taskij\n");
	MSK_deletetask(&taskij);
	mpi_one_printf("about to delete task\n");
	MSK_deletetask(&task);

	++irun;
      }

    /*************************************************/

    free(u);
    free(u_sc);
    free(eigenvalues);
    free(depsdu_sc);
    free(Altemp);
    free(Autemp);
    free(gap);
    free(obj);
    free(ayidx_ij);
    free(xx);
    free(y);
    free(delf_ij);

    MSK_deleteenv(&env);
#endif
    return 1;
}

void compute_subspace(int band, char lu, int nsp, double scale, scalar_complex *Atemp, scalar_complex *depsdu_sc, scalar_complex *u_sc, int ntot)
{
  int stride = (nsp+1)*nsp/2;
  int count = 0;
  int ridx, rband, cidx, cband;
  int n = ntot+3;

  for (ridx = 0; ridx < nsp; ++ridx)
    {
      if ( lu == 'l')
	rband = band-nsp+1+ridx;
      else
	rband = band+ridx;

      for (cidx = 0; cidx <= ridx; ++cidx)
	{

	    if ( lu == 'l')
	      {
		cband = band-nsp+1+cidx;
		/* CASSIGN_SCALAR(Atemp[count+(ntot+1)*stride],(rband==cband)? 1:0.0, 0.0); */
		/* CASSIGN_SCALAR(Atemp[count+(ntot+2)*stride], 0.0, 0.0); */
	      }
	    else
	      {
		cband = band+cidx;
		/* CASSIGN_SCALAR(Atemp[count+(ntot+1)*stride], 0.0, 0.0); */
		/* CASSIGN_SCALAR(Atemp[count+(ntot+2)*stride],(rband==cband)? -1:0.0, 0.0); */
	      }

	    material_grids_SPt_blas(Atemp+count, depsdu_sc, u_sc, scale, ntot, rband, cband, stride);

	    if ( lu == 'l')
	      {

		CASSIGN_SCALAR(Atemp[count+(n-2)*stride],(rband==cband)? 1:0.0, 0.0);
		CASSIGN_SCALAR(Atemp[count+(n-1)*stride], 0.0, 0.0);
	      }
	    else
	      {

		CASSIGN_SCALAR(Atemp[count+(n-2)*stride], 0.0, 0.0);
		CASSIGN_SCALAR(Atemp[count+(n-1)*stride],(rband==cband)? -1:0.0, 0.0);
	      }

	    ++count;
	}
    }

}

void compute_subspace_fa(int band, char lu, int nsp, double scale, scalar_complex *Atemp, scalar_complex *depsdu_sc, scalar_complex *u_sc, int ntot)
{
  int stride = (nsp+1)*nsp/2;
  int count = 0;
  int ridx, rband, cidx, cband;

  for (ridx = 0; ridx < nsp; ++ridx)
    {
      if ( lu == 'l')
	rband = band-nsp+1+ridx;
      else
	rband = band+ridx;

      for (cidx = 0; cidx <= ridx; ++cidx)
	{
	  if ( lu == 'l')
	    cband = band-nsp+1+cidx;
	  else
	    cband = band+cidx;

	    material_grids_SPt_blas(Atemp+count, depsdu_sc, u_sc, scale, ntot, rband, cband, stride);
	  ++count;
	}
    }

}

void findFail(bool *optStat, int totalStat, int totalStatRound, int *failList, int *totalFail, int offset)
{
  int i, failCnt = 0, failCntRound = 0;

  for(i=0; i<totalStat; ++i)
    {
      if(optStat[i] != 1)
	{
	  failList[failCnt] = i + offset;
	  ++failCnt;
	}
    }

  for(i=totalStat; i<totalStatRound; ++i)
    {
      if(optStat[i] != 1)
	{
	  ++failCntRound;
	}
    }
  printf("failCnt + failCntRound = %d + %d = %d, totalFail = %d\n", failCnt, failCntRound, failCnt+failCntRound, *totalFail);

  CHECK( (failCnt+failCntRound)==*totalFail, "total number of fail is wrong\n");
  totalFail = failCnt;

}

/*******************************************************************************/
number run_matgrid_optgap_lp_fa_mpi(vector3_list kpoints,
			     integer band1, integer band2,
			     integer maxrun, number utol, number sp_tol,
			     integer kk, number delta, integer optRank, char *title)
{


  int i, j, k, ii, jj, ntot, n, ngrids,nl,nu,nk,ridx,cidx,rband,cband,count,nltmpt,nutmpt,stride,irun;
    double *u,*eigenvalues,*depsdu,lambda_l,lambda_u,usum,scale;
    material_grid *grids;
    double lowtol = (double) sp_tol;
    double upptol = (double) sp_tol;
    double *gap, *obj, obj_ij;
    int maxflucfreq = 5;
    /* double_array *bvec, *cvec; */
    double_array *Blin, *Clin;
    int Nl, Nu, NB, NC, ccount;
    scalar_complex *Altemp, *Autemp, *depsdu_sc, *u_sc;
    char prefix[256];
    char solstastr[30];
    struct timeval tbegin, tend, tdiff;


 /**************/
#ifndef HAVE_MOSEK
    CHECK(0, "mosek is required for material_grids_optgap\n");
#else
    MSKrescodee  r;
    MSKenv_t     env = NULL;
    MSKtask_t    task = NULL;
    MSKtask_t    taskij = NULL;
    MSKrescodee trmcode;

    MSKint32t NUMCON;
    MSKint32t NUMVAR;

    MSKint32t *ayidx_ij;
    double *y, *xx, *yopt, *delf_ij, *delf, *lhs_ij, *lhs, *Csub;

    int numtasks, rank, trueCount_ij, trueCount, offset, NCsub, failCount, NCround;
    int *optStat_ij, *optStat;
    int *failList;
    /* double min_obj_ij, min_obj, knorm; */
    double knorm;
    struct {
        double val;
        int   rank;
    } min_obj_ij, min_obj;
 /**************/

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    grids = get_material_grids(geometry, &ngrids);
    ntot = material_grids_ntot(grids, ngrids);
    mpi_one_printf("delta = %f, delta*ntot = %f\n", delta, delta*ntot);
    nk = kpoints.num_items;
  /* u = [ubar_1, ubar_2, ... , ubar_ntot, theta] */
    n = ntot + 1;
    mpi_one_printf("number of decision variable n = %d\n",n);
    u = (double *) malloc(sizeof(double) * n);

    gap = (double *) malloc(sizeof(double) * (maxrun));
    obj = (double *) malloc(sizeof(double) * (maxrun));

    eigenvalues = (double *) malloc(sizeof(double) * num_bands);
    Altemp = (scalar_complex *) calloc(band1*(band1+1)/2*n, sizeof(scalar_complex));
    Autemp = (scalar_complex *) calloc((num_bands-band1)*(num_bands-band1+1)/2*n, sizeof(scalar_complex));

    depsdu = (double *) malloc(sizeof(double) * ntot*mdata->fft_output_size);
    depsdu_sc = (scalar_complex *) calloc(ntot*mdata->fft_output_size, sizeof(scalar_complex));
    u_sc = (scalar_complex *) calloc(ntot, sizeof(scalar_complex));

    deps_du(depsdu, 1.0, grids, ngrids);
    for(j=0; j<ntot*mdata->fft_output_size; ++j)
      depsdu_sc[j].re = depsdu[j];
    free(depsdu);

    Blin = (double_array *)malloc(sizeof(double_array));
    initArray(Blin,2);
    Clin = (double_array *)malloc(sizeof(double_array));
    initArray(Clin,2);

    ayidx_ij = (MSKint32t *) calloc(n, sizeof(MSKint32t));
    for (i=0; i<n; ++i)
	  ayidx_ij[i] = i;
    /* decision variables are [ybar, qbar, theta] */
    int numconij = 3*ntot+2;
    int numvarij = 2*ntot+1;
    xx = (double *) calloc(numvarij, sizeof(double));
    yopt = (double *) calloc(ntot, sizeof(double)); /* save the worst y of fa */
    y = (double *) calloc(numconij, sizeof(double));


   /* OPTIMIZATION */
    r = MSK_makeenv(&env,NULL);
    if (r == MSK_RES_OK)
      r = MSK_initenv(env);

    NUMVAR = ntot+1;

    usum = ntot; irun = 0;
    while (irun < maxrun && usum >= utol)
      {

	NB = 0; NC = 0;
	Blin->used = 0;
	Clin->used = 0;

        lambda_l = 0.0;
        lambda_u = 100.0;

        material_grids_get(u, grids, ngrids);
	for(j=0; j<ntot; ++j)
	  u_sc[j].re = u[j];


	for (k = 0; k < nk; ++k) {
            randomize_fields();
            solve_kpoint(kpoints.items[k]);

	    for (j = 0; j < num_bands; ++j)
	      eigenvalues[j]  = freqs.items[j]*freqs.items[j];

	    /* lambda_l and lambda_u updated at each k, only used outside the "for" loop for k */
	    lambda_l = MAX2(eigenvalues[band1-1], lambda_l);
	    lambda_u = MIN2(eigenvalues[band2-1], lambda_u);

	    /* find(const double *v, double scalev, double shiftval, int pos, const int n) */
	    nltmpt = (band1-1) - find(eigenvalues,-1,eigenvalues[band1-1]*(1-lowtol),-1,num_bands);
	    nutmpt = find(eigenvalues,1,-eigenvalues[band2-1]*(1+upptol),1,num_bands) - (band2-1);

	    nl = nltmpt>=1 ? nltmpt : nltmpt+1;
	    nu = nutmpt>=1 ? nutmpt : nutmpt+1;

	    /* mpi_one_printf("nl[%d] = %d, nu[%d] = %d\n",k,nl,k,nu); */

	    /* compute subspaces and their linear approximations */
	    scale = 1.0;
	    knorm = kpoints.items[k].x*kpoints.items[k].x + kpoints.items[k].y*kpoints.items[k].y + kpoints.items[k].z*kpoints.items[k].z;
	    Nl = 0;
	    if (knorm>1e-3) {
	      compute_subspace_fa(band1, 'l', nl, scale, Altemp, depsdu_sc, u_sc,ntot);
	      Nl = SDP2LP_v1(Blin, Altemp, nl, n, kk);
	      NB += Nl;
	    }

	    compute_subspace_fa(band2, 'u', nu, scale, Autemp, depsdu_sc, u_sc,ntot);
	    Nu = SDP2LP_v1(Clin, Autemp, nu, n, kk);
	    NC += Nu;

	    /* mpi_one_printf("Nl[%d] = %d, Nu[%d] = %d\n",k,Nl,k,Nu); */
	}
	gap[irun] = 2*(lambda_u-lambda_l)/(lambda_l+lambda_u);

	mpi_one_printf("Before maximization, %d, U = %g, L = %g, gap is %0.15g \n",irun+1, lambda_u, lambda_l, gap[irun]);


	if(rank==optRank)
	  {
	    /* Outer Optimiztaion Task */
	    if ( r== MSK_RES_OK )
	      r = MSK_maketask(env,NB*NC,NUMVAR,&task);
	    if ( r == MSK_RES_OK )
	      MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
	    /* append NUMVAR variables */
	    if ( r == MSK_RES_OK )
	      r = MSK_appendvars(task,NUMVAR);
	    /* set bounds on variables */
	    for (j=0; j<ntot && r==MSK_RES_OK; ++j)
	      r = MSK_putvarbound(task,j,MSK_BK_RA,0.0,1.0);
	    if ( r == MSK_RES_OK )
	      r = MSK_putvarbound(task,ntot,MSK_BK_FR,-MSK_INFINITY,+MSK_INFINITY);
	    /* set objectives */
	    if ( r == MSK_RES_OK )
	      r = MSK_putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE);
	    /* decision variable of the outer task */
	    /*   [u_1, u_2, ..., u_ntot, t] */
	    if ( r == MSK_RES_OK )
	      r = MSK_putcj(task,ntot,1);
	  }
	/*************************************************/
	/* Inner tasks: (i,j) subproblem */

	/* create a subtask with numvarij variables, and numconij constraints */
	r = MSK_maketask(env,numconij,numvarij,&taskij);
	/* MSK_linkfunctotaskstream(taskij,MSK_STREAM_LOG,NULL,printstr); */

	/* append numvarij variables */
	if ( r == MSK_RES_OK )
	  r = MSK_appendvars(taskij,numvarij);
	/* set lower bounds on variables, all >= 0 */
	for (i=0; i<numvarij && r==MSK_RES_OK; ++i)
	  r = MSK_putvarbound(taskij,i,MSK_BK_LO,0.0,MSK_INFINITY);

	if ( r == MSK_RES_OK )
	  r = MSK_appendcons(taskij,numconij);
	/* set row i of constraint */
	for (i=0; i<ntot && r==MSK_RES_OK; ++i)
	  {
	    /* ybar - x \theta <= qbar */
	    r = MSK_putaij(taskij, i, i, 1.0);
	    r = MSK_putaij(taskij, i, ntot+i, -1.0);
	    r = MSK_putaij(taskij, i, 2*ntot, -u[i]);
	    r = MSK_putconbound(taskij, i, MSK_BK_UP, -MSK_INFINITY, 0.0);

	    /*  ybar - x \theta >= - qbar */
	    r = MSK_putaij(taskij, ntot+i, i, -1.0);
	    r = MSK_putaij(taskij, ntot+i, ntot+i, -1.0);
	    r = MSK_putaij(taskij, ntot+i, 2*ntot, +u[i]);
	    r = MSK_putconbound(taskij, ntot+i, MSK_BK_UP, -MSK_INFINITY, 0);

	    /* ybar - e \theta <= 0 */
	    r = MSK_putaij(taskij, 2*ntot+i, i, 1.0);
	    r = MSK_putaij(taskij, 2*ntot+i, 2*ntot, -1.0);
	    r = MSK_putconbound(taskij, 2*ntot+i, MSK_BK_UP, -MSK_INFINITY, 0.0);
	  }

	/* e'x qbar - delta x ntot x \theta < = 0 */
	for (i=0; i<ntot && r==MSK_RES_OK; ++i)
	    r = MSK_putaij(taskij, 3*ntot, ntot+i, 1.0);
	r = MSK_putaij(taskij, 3*ntot, 2*ntot, -delta*ntot);

	r = MSK_putconbound(taskij, 3*ntot, MSK_BK_UP, -MSK_INFINITY, 0.0);
	r = MSK_putconbound(taskij, 3*ntot+1, MSK_BK_FX, 1, 1);

	if ( r ==MSK_RES_OK )
	  r = MSK_putobjsense(taskij,MSK_OBJECTIVE_SENSE_MINIMIZE);

	if (r == MSK_RES_OK) {
	  r = MSK_putintparam(taskij,MSK_IPAR_NUM_THREADS,1);
	  r = MSK_putintparam(taskij,MSK_IPAR_OPTIMIZER,MSK_OPTIMIZER_FREE_SIMPLEX);
	}

	/* append emptry constraint */
	if(rank==optRank)
	  {
	    if ( r == MSK_RES_OK )
	      r = MSK_appendcons(task,NB*NC);
	  }

	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	NCsub = ceil(NC*1.0/numtasks);
	NCround = NCsub*numtasks;
	printf("rank = %d, NB = %d, NC = %d, NCround = %d, NCsub = %d \n", rank, NB, NC, NCround, NCsub);

	Csub = (double *) malloc(NCsub*n*sizeof(double));
	delf_ij = (double *) malloc(NCsub*n*sizeof(double));
	lhs_ij = (double *) malloc(NCsub*sizeof(double));
	optStat_ij = (int *)calloc(sizeof(int),NCsub);
	delf = (double *) malloc(NCround*n*sizeof(double));
	lhs = (double *) malloc(NCround*sizeof(double));
	optStat = (int *)malloc(NCround*sizeof(int));

	if (NC < NCround && Clin->size < NCround*n){
	      mpi_one_printf("Expand Clin->data with zeros\n");
	      Clin->size = NCround*n;
	      Clin->data = (double *)realloc(Clin->data, NCround*n*sizeof(double));
	    }

	mpi_one_printf("scatter Clin->data\n");
	MPI_Scatter(Clin->data, NCsub*n, MPI_DOUBLE, Csub, NCsub*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	mpi_one_printf("broadcast Bin->data\n");
	MPI_Bcast(Blin->data, NB*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	offset = 0;
	min_obj_ij.rank = rank;
	for(ii=0; ii<NB; ++ii)
	  {
	    mpi_one_printf("msk subproblem (%d, %d)\n", ii, NC);
	    min_obj_ij.val = 100;
	    trueCount_ij = 0; /* trueCount = 0; */
	    for(jj=0; jj<NCsub; ++jj)
	      {
		for(j=0; j<ntot && r ==MSK_RES_OK; ++j)
		    r = MSK_putaij(taskij, 3*ntot+1, j, (Csub[jj*n+j] + Blin->data[ii*n+j]));
		r = MSK_putaij(taskij, 3*ntot+1, 2*ntot, (Csub[jj*n+ntot] + Blin->data[ii*n+ntot]));

		for(j = 0; j<ntot && r ==MSK_RES_OK; ++j)
		  r = MSK_putcj(taskij, j, 2*(Csub[jj*n+j] - Blin->data[ii*n+j]));
		r = MSK_putcj(taskij, 2*ntot, 2*(Csub[jj*n+ntot] - Blin->data[ii*n+ntot]));

		r = MSK_optimizetrm(taskij,&trmcode);
		/* MSK_solutionsummary (taskij,MSK_STREAM_MSG); */

		if ( r==MSK_RES_OK )
		  {
		    MSKsolstae solsta;

		    /* r = MSK_getsolsta (taskij,MSK_SOL_ITR,&solsta); */
		    r = MSK_getsolsta(taskij,MSK_SOL_BAS,&solsta);

		    switch(solsta)
		      {
		      case MSK_SOL_STA_OPTIMAL:
		      case MSK_SOL_STA_NEAR_OPTIMAL:
			/* MSK_getxx(taskij, MSK_SOL_ITR, xx); */
			/* MSK_getyslice(taskij, MSK_SOL_ITR, 0, 2*ntot, y); */
			/* MSK_getprimalobj (taskij, MSK_SOL_ITR, &obj_ij); */
			MSK_getxx(taskij, MSK_SOL_BAS, xx);
			MSK_getyslice(taskij, MSK_SOL_BAS, 0, 2*ntot, y);
			MSK_getprimalobj (taskij, MSK_SOL_BAS, &obj_ij);

			for(j = 0; j<ntot && r ==MSK_RES_OK; ++j)
			  delf_ij[jj*n+j] = (y[j] - y[ntot+j])*xx[2*ntot];
			delf_ij[jj*n+ntot] = -1.0; /* coefficient for t */

			/* lhs constant (lowerbound) = delf_ij'*uhat - ~f_ij(xhat)  */
			lhs_ij[jj] = cblas_ddot(ntot,delf_ij+jj*n,1,u,1);
			lhs_ij[jj] = lhs_ij[jj] - obj_ij;

			/* printf("rank = %d, (%d,%d), obj_ij = %f\n", rank, ii,jj, obj_ij); */
			if (rank*NCsub+jj < NC)
			  {
			    if (obj_ij < min_obj_ij.val)
			      {
				for(j=0; j<ntot; ++j)
				  yopt[j] = xx[j]/xx[2*ntot];
				min_obj_ij.val = obj_ij;
			     }
				/* min_obj_ij = MIN2(min_obj_ij,obj_ij); */
			  }

			optStat_ij[jj] = 1;
			++trueCount_ij;

			break;
		      case MSK_SOL_STA_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_PRIM_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
		      case MSK_SOL_STA_UNKNOWN:
		      default:
			optStat_ij[jj] = 0;
			/* print_solsta(solsta, solstastr); */
			/* printf("rank %d, %s\n", rank, solstastr); */
			break;
		      }
		  }
		else
		  mpi_one_printf("Error while optimizing.\n");
	      } /* end of jj < NCsub */

	    MPI_Allgather(delf_ij, NCsub*n, MPI_DOUBLE, delf, NCsub*n, MPI_DOUBLE, MPI_COMM_WORLD);
	    MPI_Allgather(lhs_ij, NCsub, MPI_DOUBLE, lhs, NCsub, MPI_DOUBLE, MPI_COMM_WORLD);
	    MPI_Allgather(optStat_ij, NCsub, MPI_INT, optStat, NCsub, MPI_INT, MPI_COMM_WORLD);
	    MPI_Allreduce(&trueCount_ij, &trueCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&min_obj_ij, &min_obj, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
	    /* MPI_Allreduce(&min_obj_ij, &min_obj, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); */
    	    mpi_one_printf("min_obj = %f from rank = %d\n", min_obj.val, min_obj.rank);
	    MPI_Bcast(yopt, ntot, MPI_DOUBLE, min_obj.rank, MPI_COMM_WORLD);




	    if(rank==optRank)
	      {
		for(jj=0; jj<NC & r==MSK_RES_OK; ++jj)
		  {
		    r = MSK_putarow(task, jj+offset, n, ayidx_ij, delf+jj*n);
		    r = MSK_putconbound(task, jj+offset, MSK_BK_LO, lhs[jj], MSK_INFINITY);
		  }

		/* failList stores the indices of the failed opt instances */
		failCount = MAX2(NCround-trueCount,0);
		if(failCount>0)
		  {
		    printf("trueCount = %d, failCount = %d > 0, so we are going to remove some constraints\n",trueCount, failCount);

		    failList = (int *)malloc(failCount*sizeof(int));
		    findFail(optStat, NC, NCround, failList, &failCount, offset);
		    MSK_removecons(task, failCount, failList);
		    free(failList);
		  }

		/* index of the current last constraint*/
		offset += NC-failCount;
	      }

	  } /* end of NB */
	free(Csub);
	free(delf_ij);
	free(lhs_ij);
	free(optStat_ij);

	free(delf);
	free(lhs);
	free(optStat);

	/*************************************************/
	if(r==MSK_RES_OK)
	  {
	    if(rank==optRank)
	      {
		printf("rank = %d, number of constraints in outer task = %d\n", rank, offset);

		/* OPTIMIZATION of Outer Task*/
		printf("Rank %d starting outer task\n", rank);

		if ( r==MSK_RES_OK )
		  r = MSK_putintparam(task,MSK_IPAR_OPTIMIZER,MSK_OPTIMIZER_FREE_SIMPLEX);
		else
		  printf("rank %d, error while optimizing with code %d .\n", rank, r);

		if ( r==MSK_RES_OK )
		  r = MSK_putintparam(task,MSK_IPAR_NUM_THREADS,numtasks);
		else
		  printf("rank %d, error while optimizing with code %d .\n", rank, r);

		if ( r==MSK_RES_OK )
		  /* MSKrescodee trmcode; */
		  r = MSK_optimizetrm(task,&trmcode);
		else
		  printf("rank %d, error while optimizing with code %d .\n", rank, r);

		MSK_solutionsummary (task,MSK_STREAM_MSG);

		if ( r==MSK_RES_OK )
		  {
		    printf("rank %d check msk solution\n", rank);

		    MSKsolstae solsta;
		    r = MSK_getsolsta (task,MSK_SOL_BAS,&solsta);

		    printf("rank %d got msk solution\n", rank);

		    switch(solsta)
		      {
		      case MSK_SOL_STA_OPTIMAL:
		      case MSK_SOL_STA_NEAR_OPTIMAL:
			printf("rank %d optimal solution\n", rank);

			MSK_getxx(task, MSK_SOL_BAS, xx);
			printf("rank %d, t = %g\n", rank, xx[ntot]);
			usum = 0.0;

			for(j=0; j<ntot; ++j){
			  usum += fabs(xx[j]-u[j]);
			  u[j] = xx[j];
			}
			usum /= ntot;

			MSK_getprimalobj (task, MSK_SOL_BAS, obj+irun);
			printf("On rank %d, after maximization, %d, objective is %0.15g, change_in_u = %g\n",rank, irun+1, obj[irun], usum);

			/* detect flunctuation */
			usum = detectFluctuation(gap,irun,maxflucfreq,usum,utol);

			break;
		      case MSK_SOL_STA_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_PRIM_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
		      case MSK_SOL_STA_UNKNOWN:
		      default:
			print_solsta(solsta, solstastr);
			printf("%s\n", solstastr);
			break;
		      } /* end of switch */
		  }
		 else
		   printf("rank %d, error while optimizing with code %d .\n", rank, r);

		printf("about to delete task\n");
		MSK_deletetask(&task);
	      } /* end of rank == optRank */
 	  } /* end of if r==MSK_RES_OK */
	else
	  printf("rank %d, error while optimizing with code %d .\n", rank, r);

	mpi_one_printf("about to delete taskij\n");
	MSK_deletetask(&taskij);

	mpi_one_printf("broadcasting u\n");
	MPI_Bcast(u, ntot, MPI_DOUBLE, optRank, MPI_COMM_WORLD);
	mpi_one_printf("broadcasting usum\n");
	MPI_Bcast(&usum, 1, MPI_DOUBLE, optRank, MPI_COMM_WORLD);


	/* temporarily update grid with optimal y and output the .h5 file */
	material_grids_set(yopt, grids, ngrids);
	mpi_one_printf("outputting grid-y\n");
	snprintf(prefix, 256, "%s%04d-grid-y", title, irun+1);
	save_material_grid(*grids, prefix);

	/* update u and epsilon */
	material_grids_set(u, grids, ngrids);
	/* synchronize_material_grid(grids); */
	reset_epsilon();
       	/* output epsilon file at each iteration */
	mpi_one_printf("getting eps\n");
	get_epsilon();
	mpi_one_printf("outputting eps \n");
	snprintf(prefix, 256, "%s%04d-", title, irun+1);
	output_field_to_file(-1, prefix);
	mpi_one_printf("outputting mg u\n");
	strcat(prefix,"grid");
	save_material_grid(*grids, prefix);

	++irun;
      } /* end of while irun < maxrun loop */

    /*************************************************/

    free(u);
    free(u_sc);
    free(eigenvalues);
    free(depsdu_sc);
    free(Altemp);
    free(Autemp);
    free(gap);
    free(obj);
    free(ayidx_ij);
    free(Blin);
    free(Clin);

    free(xx);
    free(y);

    MSK_deleteenv(&env);
#endif
    return 1;

}
