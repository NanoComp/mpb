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

#include "config.h"

#include <check.h>
#include <blasglue.h>
#include <matrices.h>
#include <matrixio.h>

#include "mpb.h"

sqmatrix shift_overlap(int band_start, int n_bands, int s1, int s2, int s3)
{
    sqmatrix U, S1, S2;
    int final_band = band_start + n_bands, ib;

    CHECK(mdata, "init-params must be called before shift_overlap");
    CHECK(band_start + n_bands <= num_bands, "not enough bands in shift_overlap");

    U = create_sqmatrix(n_bands);

    /* ...we have to do this in blocks of eigensolver_block_size since
       the work matrix W[0] may not have enough space to do it at once. */
    S1 = create_sqmatrix(n_bands);
    S2 = create_sqmatrix(n_bands);

    for (ib = band_start; ib < final_band; ib += Hblock.alloc_p) {
        if (ib + Hblock.alloc_p > final_band) {
            maxwell_set_num_bands(mdata, final_band - ib);
            evectmatrix_resize(&Hblock, final_band - ib, 0);
        }
        maxwell_compute_H_from_shifted_B(mdata, H, Hblock,
                                         (scalar_complex *) mdata->fft_data,
                                         s1, s2, s3,
                                         ib, 0, Hblock.p);

        evectmatrix_XtY_slice2(U, H, Hblock, band_start, 0, n_bands, Hblock.p,
                               ib-band_start, S1, S2);
    }

    /* Reset scratch matrix sizes: */
    evectmatrix_resize(&Hblock, Hblock.alloc_p, 0);
    maxwell_set_num_bands(mdata, Hblock.alloc_p);

    destroy_sqmatrix(S2);
    destroy_sqmatrix(S1);
    return U;
}

/* complex argument */
static double scalar_carg(scalar_complex z)
{
    return atan2(CSCALAR_IM(z), CSCALAR_RE(z));
}

number_list bott_indices(integer band_start, integer n_bands)
{
    sqmatrix U, W, S1, S2, Un, Wn;
    scalar_complex *eigvals;
    number_list bott;
    int i, n;

    CHECK(mdata, "init-params must be called before shift_overlap");
    CHECK(mdata->nz == 1, "Bott index is not defined (yet) for 3d");

    U = shift_overlap(band_start-1, n_bands, 0, 1, 0);
    W = shift_overlap(band_start-1, n_bands, 1, 0, 0);

    /* allocate scratch arrays */
    S1 = create_sqmatrix(n_bands);
    S2 = create_sqmatrix(n_bands);
    Un = create_sqmatrix(n_bands);
    Wn = create_sqmatrix(n_bands);
    CHK_MALLOC(eigvals, scalar_complex, n_bands);

    bott.num_items = n_bands;
    CHK_MALLOC(bott.items, number, n_bands);
    for (n = 1; n <= n_bands; ++n) { /* compute n-th bott index */
        sqmatrix_copy(Un, U);
        sqmatrix_copy(Wn, W);
        sqmatrix_resize(&Un, n, 1);
        sqmatrix_resize(&Wn, n, 1);
        sqmatrix_resize(&S1, n, 0);
        sqmatrix_resize(&S2, n, 0);
        sqmatrix_AeBC(S1, Wn, 0, Un, 0);
        sqmatrix_AeBC(S2, Wn, 1, Un, 1);
        sqmatrix_AeBC(Un, S1, 0, S2, 0);
        sqmatrix_eigenvalues(Un, eigvals);
        for (i = 0; i < n; ++i)
            bott.items[n-1] += scalar_carg(eigvals[i]);
        bott.items[i] /= 6.28318530717958647692528676655900; /* twopi */
    }
    destroy_sqmatrix(Wn);
    destroy_sqmatrix(Un);
    destroy_sqmatrix(S2);
    destroy_sqmatrix(S1);
    destroy_sqmatrix(W);
    destroy_sqmatrix(U);
    free(eigvals);
    return bott;
}
