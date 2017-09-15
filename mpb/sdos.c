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


void BtH_overlap(scalar_complex *BtH, int band_start, int n_bands, 
                 int iG1_start, int iG2_start, int iG3_start,
                 int iG1_end,   int iG2_end,   int iG3_end)
{
    int final_band = band_start + n_bands, n = 0, ib, ibb, i1, i2, i3, ixt, iyt, c;
    int nx = mdata->nx, ny = mdata->ny, nz = mdata->nz;  /* # of G-vecs used in calculation along xyz */
    int nG = (iG1_end-iG1_start+1)*(iG2_end-iG2_start+1)*(iG3_end-iG3_start+1); /* # of G vecs req'ed */
    int *ix, *iy, *iz; 
    scalar_complex polsum;
	CHK_MALLOC(ix,int,nG);
    CHK_MALLOC(iy,int,nG);
    CHK_MALLOC(iz,int,nG);

    CHECK(iG1_start <= iG1_end && iG2_start <= iG2_end && iG3_start <= iG3_end,
          "req'ed G vectors must be incrementing (start<=end)"); 
    CHECK( ((iG1_start > -nx/2) && (iG1_end <= nx/2))       || 
           ((iG1_start == iG1_end) && (iG1_end == nx-1) && (nx-1 == 0)),
           "iG1 out of bounds");
    CHECK( ((iG2_start > -ny/2) && (iG2_end <= ny/2))       || 
           ((iG2_start == iG2_end) && (iG2_end == ny-1) && (ny-1 == 0)),
           "iG2 out of bounds");
    CHECK( ((iG3_start > -nz/2) && (iG3_end <= nz/2))       || 
           ((iG3_start == iG3_end) && (iG3_end == nz-1) && (nz-1 == 0)),
           "iG3 out of bounds");

    CHECK(mdata, "init-params must be called before shift_overlap");
    CHECK(band_start + n_bands <= num_bands, "not enough bands in shift_overlap");
 

    /* -- construct lists for the _specific_ G vectors we actually want to pick out -- 
     * We adopt the convention that i{1,2,3} refers to G-vectors (sign-flip incorporated)
     *        G(i1,i2,i3) = i1*G1 + i2*G2 + i3*G3                  
     * On the other hand, i{x,y,z} refer to their associated array positions in  
     * the evectmatrices H and Hblock.                                                  */
    for (i1 = iG1_start; i1 <= iG1_end; ++i1) {
       ixt = i1 <= 0 ? -i1 : nx - i1;
       for (i2 = iG2_start; i2 <= iG2_end; ++i2) {
          iyt = i2 <= 0 ? -i2 : ny - i2;
          for (i3 = iG3_start; i3 <= iG3_end; ++i3) {
             ix[n] = ixt;
             iy[n] = iyt;
             iz[n] = i3 <= 0 ? -i3 : nz - i3;
             ++n; 
          }
       }
    }

    printf("\n");
    printf("   iG1_start = %3d, iG2_start = %3d, iG3_start = %3d\n", iG1_start, iG2_start, iG3_start);
    printf("   iG1_end =   %3d, iG2_end =   %3d, iG3_end =   %3d\n", iG1_end, iG2_end, iG3_end);
    printf("   nx =        %3d, ny =        %3d, nz =        %3d\n", nx, ny, nz);
    printf("   nG = %d, nG_avail = %d, n_bands = %d\n\n", nG, nx*ny*nz,n_bands);

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
          for (ibb = ib; ibb < ib + Hblock.p; ++ibb) { /* bands in the block */
             CASSIGN_ZERO(polsum);
             for (c = 0; c < 2; ++c) {    /* polarization (transverse basis) */
                CACCUMULATE_SUM_CONJ_MULT( polsum, /* a += conj(b) * c for (a,b,c) call */
                            H.data[(((ix[n]*ny+iy[n])*nz+iz[n])*2+c)*H.p+ibb],
	                   Hblock.data[(((ix[n]*ny+iy[n])*nz+iz[n])*2+c)*Hblock.p+ibb-ib] );
             }
             BtH[n*n_bands+ibb-band_start] = polsum; /* row-major */ 
          }
       }
    }
    

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


/* calculate the spectral density of states (sdos) and write it 
 * to an .h5 file */
void get_sdos(number freq_min, number freq_max, integer freq_num, 
                     number eta,
                     integer band_start, integer n_bands, 
                     vector3 iG_start, vector3 iG_end) 
{
    int i,b,n;
	int iG1_start = myround(iG_start.x), iG2_start = myround(iG_start.y), iG3_start = myround(iG_start.z); 
	int iG1_end = myround(iG_end.x),     iG2_end = myround(iG_end.y),     iG3_end = myround(iG_end.z); 
    int nG1 = iG1_end-iG1_start+1,       nG2 = iG2_end-iG2_start+1,       nG3 = iG3_end-iG3_start+1; 
    int nG = nG1 * nG2 * nG3; /* total # of G vecs req'ed */
    scalar_complex *BtH, ctemp;
    real npref = 2*Vol/3.141592653589793, *spanfreqs, *spanfreqs2, *freqs2, df, fpref;
	real *sdos; 
	int iodims[1] = {freq_num*nG}, iostart[1] = {0}; /* for .h5 write */
    matrixio_id file_id, data_id; 
    
    /* allocate arrays */ 
    CHK_MALLOC(BtH, scalar_complex, nG*n_bands); 
    CHK_MALLOC(spanfreqs,  real, freq_num);
    CHK_MALLOC(spanfreqs2, real, freq_num);
    CHK_MALLOC(freqs2,     real, n_bands);
    CHK_MALLOC(sdos,       real, freq_num*nG);
 
    /* create a frequency array that spans freq_min to freq_max in freq_num steps */
    spanfreqs[0] = freq_min;
    df = (freq_max-freq_min)/(freq_num-1);
    for (i = 1; i < freq_num; ++i) 
       spanfreqs[i] = spanfreqs[i-1]+df;
    for (i = 0; i < freq_num; ++i) 
       spanfreqs2[i] = pow(spanfreqs[i],2);
    
    
    /* get the squared eigenfrequencies starting from band_start */
    for (b = 0; b < n_bands; ++b) 
	   freqs2[b] = pow(freqs.items[b+band_start],2);

    /* compute the overlap between G-components of bands */
    BtH_overlap(BtH, band_start, n_bands, iG1_start, iG2_start, iG3_start,
                iG1_end, iG2_end, iG3_end);
     
    /* for (b = 0; b < n_bands; ++b) {
       CASSIGN_ZERO(ctemp);
       for (n = 0; n < nG; ++n) {
          CACCUMULATE_SUM(ctemp,BtH[n*n_bands+b]);
       }
       printf("   G-summed BtH at band %d = %.5g+i%.5g\n", b+1, CSCALAR_RE(ctemp), CSCALAR_IM(ctemp));
    } */

    /* calculate the sdos by summation of terms */
    for (i = 0; i < freq_num; ++i) {      /* frequency loop */
       fpref = npref*spanfreqs[i];
       for (n = 0; n < nG; ++n) {         /* G-vector loop */
          sdos[i*nG+n] = 0;               /* init = 0 before adding anything */
          for (b = 0; b < n_bands; ++b) { /* band loop */
             CASSIGN_SCALAR(ctemp, freqs2[b]-spanfreqs2[i], -eta); /* = omegan^2 - omega^2 - i*eta (= denom) */
             CASSIGN_DIV(ctemp, *(BtH+n*n_bands+b), ctemp); /* = BtH/denom (= fraction) */
             sdos[i*nG+n] += CSCALAR_IM(ctemp);  /* += Im(fraction) | sum over included bands */
          }
          sdos[i*nG+n] *= fpref; /* multiply by prefactor */
       }
    }

    /* write contents of sdos to hdf5 format; so just as one long 1d array to */ 
    file_id = matrixio_create("temporary_save_name");
    data_id = matrixio_create_dataset(file_id,"sdos",NULL, 1, iodims);
    matrixio_write_real_data(data_id, iodims, iostart, 1, sdos); 
    matrixio_close_dataset(data_id);
    matrixio_close(file_id);
}
