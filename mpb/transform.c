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



/* Right now, this assumes that the user loads the bfield or dfield via (get-hfield band) or
   (get-field band) beforehand. 
   Instead, this could be implemented by just computing this here, e.g. by having us calling
   get_hfield or get_bfield; that seems to achieve the same (i.e. setting curfield(_type)) */
cnumber transformed_overlap(matrix3x3 W, vector3 w)
{
    /* TODO: Could we just call get_dfield or get_hfield here instead of having the user do it? 
             then the input would have a choice between e or d type input */
    int n1, n2, n3;
    real s1, s2, s3, c1, c2, c3;
    cnumber integral = {0,0};
    number detW;
    vector3 kvector = cur_kvector;
    matrix3x3 invW;

    if (!curfield || !strchr("db", curfield_type)) {
        mpi_one_fprintf(stderr, "a real-space H or D-field must be loaded beforehand.\n");
        return integral;
    }
    #ifdef HAVE_MPI
        CHECK(0, "transformed_overlap(..) not yet implemented for MPI!"); /* TODO: Remove if PR#112 is merged */
    #endif
    
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
        mpi_one_fprintf(stderr, "valid symmetry operations {W|w} must have |det(W)| = 1.0\n");
        return integral;
    }

    kvector.x *= TWOPI/geometry_lattice.size.x; /* hoist these rescalings outside the  */
    kvector.y *= TWOPI/geometry_lattice.size.y; /* loop; might be that licm takes care */
    kvector.z *= TWOPI/geometry_lattice.size.z; /* of it - but better to be sure       */


    /* Loop over coordinates (introduces int vars i1, i2, i3, xyz_index) */
    LOOP_XYZ(mdata) { /* implies two opening braces '{{' */
        vector3 p, pt;
        scalar_complex F[3], Ft[3], integrand, phase;
        double deltaphi;

        /* Current lattice position in loop */
        p.x = i1 * s1 - c1;
        p.y = i2 * s2 - c2;
        p.z = i3 * s3 - c3;

        /* Transformed coordinate pt = invW*p - w */
        pt = vector3_minus(matrix3x3_vector3_mult(invW, p), w);

        /* Bloch field value at transformed coordinate pt: interpolation is needed to ensure
           generality in the case of fractional translations. Unfortunately, this currently
           precludes compatibility with MPI, since get_val is not implemented for MPI */
        /* TODO: Remove disclaimer above if PR#112 is merged */
        get_bloch_field_point_(Ft, pt); /* sets Ft to F at p [excludes exp(ikr) factor] */

        /* Transform the vector components of Ft by W; a bit tedious to do, since 
           there are no matrix3x3*cvector3 routines, nor any CASSIGN_CVECTOR_RE/IM.
           Instead, we do manually what is done in matrix3x3_vector3_mult, twice  */
        Ft[0].re = W.c0.x*Ft[0].re + W.c1.x*Ft[1].re + W.c2.x*Ft[2].re;
        Ft[0].im = W.c0.x*Ft[0].im + W.c1.x*Ft[1].im + W.c2.x*Ft[2].im;
        Ft[1].re = W.c0.y*Ft[0].re + W.c1.y*Ft[1].re + W.c2.y*Ft[2].re;
        Ft[1].im = W.c0.y*Ft[0].im + W.c1.y*Ft[1].im + W.c2.y*Ft[2].im;
        Ft[2].re = W.c0.z*Ft[0].re + W.c1.z*Ft[1].re + W.c2.z*Ft[2].re;
        Ft[2].im = W.c0.z*Ft[0].im + W.c1.z*Ft[1].im + W.c2.z*Ft[2].im;
        
        /* Get the Bloch field value at current point p [exludes exp(ikr) factor].
           We multiply the input field F (either B or D-field) with μ⁻¹ or ε⁻¹ to get
           H- or E-fields, as the relevant overlap is ⟨F|Ft⟩ = (H|Bt) or (E|Dt), with
           t-postscript denoting a field transformed by {W|w}. Here, we essentially
           adapt some boiler-plate code from compute_field_energy_internal in fields.c   */
        if (curfield_type == 'd') {
            assign_symmatrix_vector(F, mdata->eps_inv[xyz_index], curfield+3*xyz_index);
        }
        else if (curfield_type == 'b' && mdata->mu_inv != NULL) {
            assign_symmatrix_vector(F, mdata->mu_inv[xyz_index],  curfield+3*xyz_index);
        }
        else {
            F[0] = curfield[3*xyz_index];
            F[1] = curfield[3*xyz_index+1];
            F[2] = curfield[3*xyz_index+2];
        }
        
        /* Inner product of F and Ft={W|w}F in Bloch form */
        CASSIGN_CONJ_MULT(integrand, F[0], Ft[0]);          /* add adjoint(F)*Ft to integrand */
        CACCUMULATE_SUM_CONJ_MULT(integrand, F[1], Ft[1]);
        CACCUMULATE_SUM_CONJ_MULT(integrand, F[2], Ft[2]);

        /* So far, we have excluded the Bloch phases; they must be included, however.
           It saves two trigonometric operations if we do them just jointly for p and pt.
           Note that rescaling by TWOPI/geometry_lattice.xyz is hoisted outside loop      */
        deltaphi = (kvector.x*(pt.x-p.x) + kvector.y*(pt.y-p.y) + kvector.z*(pt.z-p.z));
        CASSIGN_SCALAR(phase, cos(deltaphi), sin(deltaphi));

        /* Add integrand-contribution to integral, keeping Bloch phases in mind */
        integral.re += CSCALAR_MULT_RE(integrand, phase);
        integral.im += CSCALAR_MULT_IM(integrand, phase); 
    }}}

    integral.re *= Vol / H.N;
    integral.im *= Vol / H.N;

    /* C89 requires new type declarations to occur only at the start of a blocks; 
       hence {...} here */
    { 
    /* Also, we don't really need to do this, since we don't support MPI anyway; 
       still works in serial though, so might as well keep it around */
    /* TODO: Remove disclaimer above if PR#112 is merged */
    cnumber integral_sum;
    mpi_allreduce(&integral, &integral_sum, 2, number,
                  MPI_DOUBLE, MPI_SUM, mpb_comm);

    if (curfield_type == 'b') {  /* H & B are pseudovectors => transform includes det(W) */
        integral_sum.re *= detW; 
        integral_sum.im *= detW;
    }

    return integral_sum;
    }
}
