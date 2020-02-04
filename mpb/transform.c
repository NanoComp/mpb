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



/* compute and return adjoint(cv1)*cv2 */
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
    vector3 kvector = cur_kvector;
    matrix3x3 invW;
    number detW;

    if (!curfield || !strchr("hb", curfield_type)) {
        mpi_one_fprintf(stderr, "The real-space H (or B) field must be loaded first.\n");
        return integral;
    }
    #ifdef HAVE_MPI
        CHECK(0, "transformed_overlap(..) not yet implemented for MPI!");
    #endif
    /* TODO: should also check and error if SCALAR_COMPLEX not defined, or if mu_inv != NULL */
    
    /* prepare before looping ... */
    n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;

    s1 = geometry_lattice.size.x / n1; /* pixel spacings */
    s2 = geometry_lattice.size.y / n2;
    s3 = geometry_lattice.size.z / n3;
    c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5; /* offsets (negative) */
    c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
    c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

    invW = matrix3x3_inverse(W);
    detW = matrix3x3_determinant(W);
    if (fabs(fabs(detW) - 1.0) > 1e-12) {
        mpi_one_fprintf(stderr, "Valid symmetry operations {W|w} must have |det(W)| = 1.0\n");
        return integral;
    }

    kvector.x *= TWOPI/geometry_lattice.size.x; /* hoist these rescalings outside the  */
    kvector.y *= TWOPI/geometry_lattice.size.y; /* loop; might be that licm takes care */
    kvector.z *= TWOPI/geometry_lattice.size.z; /* of it - but better to be sure       */


    /* Loop over coordinates (introduces int vars i1, i2, i3, & xyz_index) */
    LOOP_XYZ(mdata) { /* implies two opening braces '{{' */
        vector3 p, pt;
        cvector3 F, Ft;
        cnumber integrand;
        double deltaphi;
        scalar_complex phase;

        /* Current lattice position in loop */
        p.x = i1 * s1 - c1;
        p.y = i2 * s2 - c2;
        p.z = i3 * s3 - c3;

        /* Bloch field value at current point p (exludes exp(ikr) factor) */
        F.x = cscalar2cnumber(curfield[3*xyz_index+0]); /* since F.x is a cnumber and curfield is a */
        F.y = cscalar2cnumber(curfield[3*xyz_index+1]); /* scalar_complex; we have to convert here  */
        F.z = cscalar2cnumber(curfield[3*xyz_index+2]); 

        /* Transformed coordinate pt = invW*p - w */
        pt = vector3_minus(matrix3x3_vector3_mult(invW, p), w);

        /* Field value at transformed coordinate pt; interpolation is needed to ensure
           generality in the case of fractional translations. Unfortunately, this
           precludes compatibility with MPI, since get_val is not implemented for MPI */
        Ft = get_bloch_field_point(pt); /* excludes exp(ikr) factor */

        /* Transform the vector components of Ft by W; a bit tedious to do, since 
           there are no matrix3x3*cvector3 routines, nor any CASSIGN_CVECTOR_RE/IM.
           Instead, we do manually what is done in matrix3x3_vector3_mult, twice  */
        Ft.x.re = W.c0.x*Ft.x.re + W.c1.x*Ft.y.re + W.c2.x*Ft.z.re;
        Ft.y.re = W.c0.y*Ft.x.re + W.c1.y*Ft.y.re + W.c2.y*Ft.z.re;
        Ft.z.re = W.c0.z*Ft.x.re + W.c1.z*Ft.y.re + W.c2.z*Ft.z.re;
        Ft.x.im = W.c0.x*Ft.x.im + W.c1.x*Ft.y.im + W.c2.x*Ft.z.im;
        Ft.y.im = W.c0.y*Ft.x.im + W.c1.y*Ft.y.im + W.c2.y*Ft.z.im;
        Ft.z.im = W.c0.z*Ft.x.im + W.c1.z*Ft.y.im + W.c2.z*Ft.z.im;

        /* Inner product of H and OH (for O = {W|w}), in Bloch form */
        integrand = cvector3_cdot(F, Ft);

        /* So far, we have excluded the Bloch phases; they must be included, however.
           It saves two trigonometric operations if we do them just jointly for p and pt.
           Note that rescaling by TWOPI/geometry_lattice.xyz is hoisted outside loop      */
        deltaphi = (kvector.x*(-p.x+pt.x) + kvector.y*(-p.y+pt.y) + kvector.z*(-p.z+pt.z));
        CASSIGN_SCALAR(phase, cos(deltaphi), sin(deltaphi));

        /* Add integrand-contribution to integral, keeping Bloch phases in mind */
        integral.re += integrand.re*phase.re - integrand.im*phase.im;
        integral.im += integrand.re*phase.im + integrand.im*phase.re;
    }}}

    integral.re *= Vol / H.N;
    integral.im *= Vol / H.N;

    /* C89 requires new type declarations to occur only at the start of a blocks; 
       hence {...} here */
    { 
    /* Also, we don't really need to do this, since we don't support MPI anyway; 
       still works in serial though, so might as well keep it around */
    cnumber integral_sum;
    mpi_allreduce(&integral, &integral_sum, 2, number,
                  MPI_DOUBLE, MPI_SUM, mpb_comm);

    integral_sum.re *= detW; /* H & B are pseudovectors => transform includes det(W) */
    integral_sum.im *= detW;

    return integral_sum;
    }
}
