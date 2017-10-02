/* Copyright (C) 1999-2014 Massachusetts Institute of Technology.
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "config.h"

#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <matrixio.h>

#include "mpb.h"
#include "matrixio.h"


void BtH_overlap(scalar *BtH, int band_start, int n_bands,
                 int iG1_min, int iG2_min, int iG3_min,
                 int iG1_max, int iG2_max, int iG3_max)
{
    int final_band = band_start + n_bands, n = 0, ib, ibb, i1, i2, i3, ixt, iyt, c;
    int local_x_start = mdata->local_x_start, local_nx = mdata->local_nx;
    int nx = mdata->nx, ny = mdata->ny, nz = mdata->nz;  /* # of G-vecs used in calculation along xyz */
    int nG = (iG1_max-iG1_min+1)*(iG2_max-iG2_min+1)*(iG3_max-iG3_min+1); /* # of G vecs req'ed */
    int *ix, *iy, *iz;
    scalar *BtH_local;
	 CHK_MALLOC(ix,int,nG);
    CHK_MALLOC(iy,int,nG);
    CHK_MALLOC(iz,int,nG);

    CHECK(iG1_min <= iG1_max && iG2_min <= iG2_max && iG3_min <= iG3_max,
          "req'ed G vectors must be incrementing (min<=max)");
    CHECK( ((iG1_min > -nx/2) && (iG1_max <= nx/2))       ||
           ((iG1_min == iG1_max) && (iG1_max == nx-1) && (nx-1 == 0)),
           "iG1 out of bounds");
    CHECK( ((iG2_min > -ny/2) && (iG2_max <= ny/2))       ||
           ((iG2_min == iG2_max) && (iG2_max == ny-1) && (ny-1 == 0)),
           "iG2 out of bounds");
    CHECK( ((iG3_min > -nz/2) && (iG3_max <= nz/2))       ||
           ((iG3_min == iG3_max) && (iG3_max == nz-1) && (nz-1 == 0)),
           "iG3 out of bounds");

    CHECK(mdata, "init-params must be called before shift_overlap");
    CHECK(band_start + n_bands <= num_bands, "not enough bands in shift_overlap");


    /* -- construct lists for the _specific_ G vectors we actually want to pick out --
     * We adopt the convention that i{1,2,3} refers to G-vectors (sign-flip incorporated)
     *        G(i1,i2,i3) = i1*G1 + i2*G2 + i3*G3
     * On the other hand, i{x,y,z} refer to their associated array positions in
     * the evectmatrices H and Hblock.                                                  */
    for (i1 = iG1_min; i1 <= iG1_max; ++i1) {
       ixt = i1 <= 0 ? -i1 : nx - i1;
       for (i2 = iG2_min; i2 <= iG2_max; ++i2) {
          iyt = i2 <= 0 ? -i2 : ny - i2;
          for (i3 = iG3_min; i3 <= iG3_max; ++i3) {
             ix[n] = ixt;
             iy[n] = iyt;
             iz[n] = i3 <= 0 ? -i3 : nz - i3;
             ++n;
          }
       }
    }

    mpi_one_printf("\nSDOS calculation info:\n");
    mpi_one_printf("   iG1_min = %3d, iG2_min = %3d, iG3_min = %3d\n", iG1_min, iG2_min, iG3_min);
    mpi_one_printf("   iG1_max = %3d, iG2_max = %3d, iG3_max = %3d\n", iG1_max, iG2_max, iG3_max);
    mpi_one_printf("   nx =      %3d, ny =      %3d, nz =      %3d\n", nx, ny, nz);
    mpi_one_printf("   nG = %d, nG_avail = %d, n_bands = %d\n\n", nG, nx*ny*nz,n_bands);

    /* temporary BtH array since MPI summation needs to be done out-of-place */
    BtH_local = (scalar*) calloc(nG*n_bands, sizeof(scalar));
    CHECK(BtH_local || nG*n_bands == 0, "out of memory");

    /* ...we have to do this in blocks of eigensolver_block_size since
     * the work matrix W[0] may not have enough space to do it at once. */
    for (ib = band_start; ib < final_band; ib += Hblock.alloc_p) {
       if (ib + Hblock.alloc_p > final_band) {
           maxwell_set_num_bands(mdata, final_band - ib);
           evectmatrix_resize(&Hblock, final_band - ib, 0);
       }

       /* Beware the notation: H is really B and Hblock is blocks of H */
       maxwell_compute_H_from_B(mdata, H, Hblock,
                                (scalar_complex *) mdata->fft_data,
                                ib, 0, Hblock.p);

       /* Now we calculate a matrix BtH = Tr(H*adjoint(B)) in the G basis:
        * the trace works over the polarization subspace (c-indices)      */
       for (n = 0; n < nG; ++n) {     /* req'd G-vecs, indices from i{x,y,z} */
         int ixn = ix[n] - local_x_start, iyn = iy[n], izn = iz[n];
         if (ixn >= 0 && ixn < local_nx) {
            for (ibb = ib; ibb < ib + Hblock.p; ++ibb) { /* bands in the block */
               scalar polsum;
               ASSIGN_ZERO(polsum);
               for (c = 0; c < 2; ++c) {    /* polarization (transverse basis) */
                  ACCUMULATE_SUM_CONJ_MULT(polsum, /* a += conj(b) * c for (a,b,c) call */
                     H.data[(((ixn*ny+iyn)*nz+izn)*2+c)*H.p           +ibb   ],
                     Hblock.data[(((ixn*ny+iyn)*nz+izn)*2+c)*Hblock.p +ibb-ib] );
               }
               BtH_local[n*n_bands+ibb-band_start] = polsum; /* row-major */
            }
         }
      }
    }

    /* sum BtH_local arrays from all processes into BtH */
    mpi_allreduce(BtH_local, BtH, nG*n_bands * SCALAR_NUMVALS,
                  real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
    free(BtH_local);

    /* Reset matrix sizes and free memory */
    evectmatrix_resize(&Hblock, Hblock.alloc_p, 0);
    maxwell_set_num_bands(mdata, Hblock.alloc_p);
    free(ix);
    free(iy);
    free(iz);
}


/* round to nearest integer; for conversion of float vector3
 * entries to integer values (only to be safe, in case of
 * floating-point mischief */
int myround(float x){
    assert(x >= INT_MIN - 0.5);
	assert(x <= INT_MAX + 0.5);
    if (x >= 0)
       return (int) (x+0.5);
    return (x-0.5); /* else (for negative values) */
}


/* Prepend the prefix to the fname, and (if parity_suffix is true)
 * append a parity specifier (if any) (e.g. ".te"), returning a new
 * string, which should be deallocated with free().  fname or prefix
 * may be NULL, in which case they are treated as the empty string. */
/* copied from fields.c */
static char *fix_fname(const char *fname, const char *prefix,
		       maxwell_data *d, int parity_suffix)
{
     char *s;
     CHK_MALLOC(s, char,
		(fname ? strlen(fname) : 0) +
		(prefix ? strlen(prefix) : 0) + 20);
     strcpy(s, prefix ? prefix : "");
     strcat(s, fname ? fname : "");
     if (parity_suffix && d->parity != NO_PARITY) {
	  /* assumes parity suffix is less than 20 characters;
  	   * currently it is less than 12 */
	  strcat(s, ".");
	  strcat(s, parity_string(d));
     }
     return s;
}


/* calculate the spectral density of states (sdos) and write it
 * to an .h5 file. The variable band_start follows the Guile
 * convention and starts at 1 (conversion to 0-index'ing inside) */
void get_sdos(number freq_min, number freq_max, integer freq_num,
              number eta,
              integer band_start, integer n_bands,
              vector3 iG_min, vector3 iG_max,
              char *saveprefix)
{
    int i,b,n;
	int iG1_min = myround(iG_min.x), iG2_min = myround(iG_min.y), iG3_min = myround(iG_min.z);
	int iG1_max = myround(iG_max.x), iG2_max = myround(iG_max.y), iG3_max = myround(iG_max.z);
    int nG1 = iG1_max-iG1_min+1,     nG2 = iG2_max-iG2_min+1,     nG3 = iG3_max-iG3_min+1;
    int nG = nG1 * nG2 * nG3; /* total # of G vecs req'ed */
    real iGspan[6] = {iG1_min, iG1_max, iG2_min, iG2_max, iG3_min, iG3_max}, /* requires C99 (or >) */
         freqspan[3] = {freq_min, freq_max, freq_num};
    scalar *BtH, ctemp;
    real *spanfreqs, *spanfreqs2, *freqs2re, *freqs2im, *sdos,
         df, fpref, npref = 2*Vol/3.141592653589793;
	int iodims0[1] = {freq_num*nG}, iodims1[1] = {3}, iodims2[1] = {6},
        iodims3[1] = {3}, iodims4[1] = {1}, iostart[1] = {0}; /* for .h5 write */
    matrixio_id file_id, data_id;
    char *savename;

    /* allocate arrays */
    CHK_MALLOC(BtH, scalar, nG*n_bands);
    CHK_MALLOC(spanfreqs,  real, freq_num);
    CHK_MALLOC(spanfreqs2, real, freq_num);
    CHK_MALLOC(freqs2re,   real, n_bands);
    CHK_MALLOC(freqs2im,   real, n_bands);
    CHK_MALLOC(sdos,       real, freq_num*nG);
    CHK_MALLOC(savename,   char, 256); /* assume savename less than 256 chars */

    /* in Guile interface, band_start is considered to start at 1; in C it is zero-index */
    band_start -= 1;

    /* create a frequency array that spans freq_min to freq_max in freq_num steps */
    spanfreqs[0] = freq_min;
    df = (freq_max-freq_min)/(freq_num-1);
    for (i = 1; i < freq_num; ++i)
       spanfreqs[i] = spanfreqs[i-1]+df;
    for (i = 0; i < freq_num; ++i)
       spanfreqs2[i] = pow(spanfreqs[i],2);

    /* get the squared eigenfrequencies (real and imag), starting from band_start */
    for (b = 0; b < n_bands; ++b){
	   freqs2re[b] = pow(freqs.items[b+band_start],2) - pow(eta,2);
       freqs2im[b] = -2*freqs.items[b+band_start]*eta;
    }

    /* compute the overlap between G-components of bands */
    BtH_overlap(BtH, band_start, n_bands, iG1_min, iG2_min, iG3_min,
                iG1_max, iG2_max, iG3_max);

    /* print G-sum of BtH (per band) shows how much of each band the
     * req'd G-vectors effectively cover (= 1 is all of it) */
    mpi_one_printf("G-summed BtH of bands in SDOS (sum over req'd G-vectors)\n   (band, |sum|^2) = ");
    for (b = 0; b < n_bands; ++b) {
       ASSIGN_ZERO(ctemp);
       for (n = 0; n < nG; ++n) {
          ACCUMULATE_SUM(ctemp, BtH[n*n_bands+b]);
       }
       mpi_one_printf("(%d, %.3g)%s", b+1, SCALAR_NORMSQR(ctemp), (b == (n_bands-1)) ? "\n" : ", ");
    }

    /* calculate the sdos by summation of terms */
    for (i = 0; i < freq_num; ++i) {      /* frequency loop */
       fpref = npref*spanfreqs[i];
       for (n = 0; n < nG; ++n) {         /* G-vector loop */
          sdos[i*nG+n] = 0;               /* init = 0 before adding anything */
          for (b = 0; b < n_bands; ++b) { /* band loop */
             ASSIGN_SCALAR(ctemp, /* = (omegan-i*eta)^2 - omega^2 (= denom) */
                           freqs2re[b]-spanfreqs2[i],
                           freqs2im[b]);
             ASSIGN_DIV(ctemp, *(BtH+n*n_bands+b), ctemp); /* = BtH/denom (= fraction) */
             sdos[i*nG+n] += SCALAR_IM(ctemp);  /* += Im(fraction) | sum over included bands */
          }
          sdos[i*nG+n] *= fpref; /* multiply by prefactor */
       }
    }

    /* construct a savename from saveprefix and misc meta info */
    sprintf(savename, "sdos.k%d", kpoint_index);
    savename = fix_fname(savename, saveprefix, mdata, 1);

    if (mpi_is_master()) {
       /* write contents of sdos to hdf5 format; so far just as one long 1d array to */
       file_id = matrixio_create_serial(savename);

       data_id = matrixio_create_dataset(file_id, "sdos", "remember to unfold", 1, iodims0);
       matrixio_write_real_data(data_id, iodims0, iostart, 1, sdos);

       data_id = matrixio_create_dataset(file_id, "freqspan", /* write some meta data as well */
                                "freq_min, freq_max, freq_num", 1, iodims1);
       matrixio_write_real_data(data_id, iodims1, iostart, 1, freqspan);

       data_id = matrixio_create_dataset(file_id, "iGspan",
                    "iG1_min, iG1_max, iG2_min, iG2_max, iG3_min, iG3_max", 1, iodims2);
       matrixio_write_real_data(data_id, iodims2, iostart, 1, iGspan);

       data_id = matrixio_create_dataset(file_id, "kpoint", "kx, ky, kz", 1, iodims3);
       matrixio_write_real_data(data_id, iodims3, iostart, 1, mdata->current_k);

       data_id = matrixio_create_dataset(file_id, "eta", "imag part in omegan - i*eta", 1, iodims4);
       matrixio_write_real_data(data_id, iodims4, iostart, 1, (real *) &eta);

       matrixio_close_dataset(data_id);
       matrixio_close(file_id);
   }

    /* free up allocated memory from heap */
    free(BtH);
    free(spanfreqs);
    free(spanfreqs2);
    free(freqs2re);
    free(freqs2im);
    free(sdos);
}
