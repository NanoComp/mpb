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
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>

#include "mpb.h"
#include <ctl-io.h>
#include <xyz_loop.h>
#include <maxwell.h>



/* compute and return adjoint(cv1)*cv2*/
static cnumber cvector3_cdot(cvector3 cv1, cvector3 cv2)
{
    cnumber dest;
    dest.re = (cv1.x.re*cv2.x.re + cv1.y.re*cv2.y.re + cv1.z.re*cv2.z.re +
               cv1.x.im*cv2.x.im + cv1.y.im*cv2.y.im + cv1.z.im*cv2.z.im);
    dest.im = (cv1.x.re*cv2.x.im + cv1.y.re*cv2.y.im + cv1.z.re*cv2.z.im -
               cv1.x.im*cv2.x.re - cv1.y.im*cv2.y.re - cv1.z.im*cv2.z.re);
    return dest;
}


/* Right now, this assumes that the user loads the hfield via (get-hfield band) beforehand. 
   Probably, this ought to be implemented instead by just computing this internally, here,
   using maxwell_compute_h_from_H; this would also be more natural for eventually supporting
   magnetic materials. Similarly, it would generalize to D and E fields. */
cnumber transformed_overlap(matrix3x3 W, vector3 w)
{
    #ifdef HAVE_MPI
        int local_n2, local_y_start;
    #endif
    int n1, n2, n3;
    real s1, s2, s3, c1, c2, c3;
    cnumber integral = {0,0};
    matrix3x3 invW;
    number detW;

    if (!curfield || !strchr("hb", curfield_type)) {
        mpi_one_fprintf(stderr, "The real-space H (or B) field must be loaded first.\n");
        return integral;
    }
    #ifdef HAVE_MPI
        CHECK(0, "transformed_overlap(..) not yet implemented for MPI!");
    #endif

    /* should also check and error if SCALAR_COMPLEX not defined, or if mu_inv != NULL */

    n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;

    s1 = geometry_lattice.size.x / n1;
    s2 = geometry_lattice.size.y / n2;
    s3 = geometry_lattice.size.z / n3;
    c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5;
    c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
    c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

    invW = matrix3x3_inverse(W);
    detW = matrix3x3_determinant(W);

    LOOP_XYZ(mdata) {
        vector3 p, pt;
        cvector3 F, Ft;
        cnumber integrand;

        p.x = i1 * s1 - c1; 
        p.y = i2 * s2 - c2; 
        p.z = i3 * s3 - c3;

        /* Field value at current point p */
        F.x = cscalar2cnumber(curfield[3*xyz_index+0]); /* since F.x is a cnumber and curfield is a */
        F.y = cscalar2cnumber(curfield[3*xyz_index+1]); /* scalar_complex; we have to convert here  */
        F.z = cscalar2cnumber(curfield[3*xyz_index+2]); 

        /* First, obtain new transformed coordinate pt = invW*p - w */
        pt = vector3_minus(matrix3x3_vector3_mult(invW, p), w);

        /* Next, obtain field value at transformed coordinate pt; interpolation is 
           needed to ensure generality in the case of fractional translations. 
           Unfortunately, that is _NOT_ compatible with MPI, since get_val is not
           implemented for MPI */
        Ft = get_bloch_field_point(pt); 

        /* Transform the components of Ft by W; this is a bit tedious to do, since 
           there are no matrix3x3*cvector3 routines, nor any CASSIGN_CVECTOR_RE.
           Instead, we do manually what is done in matrix3x3_vector3_mult, twice  */
        Ft.x.re = W.c0.x * Ft.x.re + W.c1.x * Ft.y.re + W.c2.x * Ft.z.re;
        Ft.y.re = W.c0.y * Ft.x.re + W.c1.y * Ft.y.re + W.c2.y * Ft.z.re;
        Ft.z.re = W.c0.z * Ft.x.re + W.c1.z * Ft.y.re + W.c2.z * Ft.z.re;
        Ft.x.im = W.c0.x * Ft.x.im + W.c1.x * Ft.y.im + W.c2.x * Ft.z.im;
        Ft.y.im = W.c0.y * Ft.x.im + W.c1.y * Ft.y.im + W.c2.y * Ft.z.im;
        Ft.z.im = W.c0.z * Ft.x.im + W.c1.z * Ft.y.im + W.c2.z * Ft.z.im;

        integrand = cvector3_cdot(F, Ft);
        integral.re += integrand.re;
        integral.im += integrand.im;
    }}}

    integral.re *= Vol / H.N;
    integral.im *= Vol / H.N;

    { /* C89 requires new declarations to exist only in blocks; hence the {...} here */
    /* this isn't really needed, since we don't support MPI anyway; still works
       in serial though, so might as well keep it around */
    cnumber integral_sum;
    mpi_allreduce(&integral, &integral_sum, 2, number,
                  MPI_DOUBLE, MPI_SUM, mpb_comm);
    integral_sum.re *= detW;
    integral_sum.im *= detW;

    return integral_sum;
    }
}

