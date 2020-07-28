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



/* If curfield is set to the real-space D-field (of band-index i), computes the overlap 
         ∫ Eᵢ(r){W|w}Dᵢ(r) dr  =  ∫ Eᵢ(r)(WDᵢ)({W|w}⁻¹r) dr,
   for a symmetry operation {W|w} with point-group part W and translation part w; each 
   specified in the lattice basis. The vector fields Eᵢ and Dᵢ include Bloch phases.
   If instead curfield is set to the real-space B-field, the overlap
         ∫ Hᵢ(r){W|v}Bᵢ(r) dr  =  det(W) ∫ Hᵢ(r)(WBᵢ)({W|w}⁻¹r) dr,
   is computed instead. Note that a factor det(W) is then included since B & H are
   pseudovectors. As a result, the computed symmetry expectation values are independent
   of whether the D- or B-field is used.
   No other choices for curfield are allowed: to set curfield to the real-space B- or D-
   field call get-bfield/dfield in Scheme, or get_bfield/dfield in C.                   
   Usually, it will be more convenient to use the accessor compute_symmetry(i, W, w) 
   which defaults to the B-field (since μ = 1 usually) and instead takes a band-index i .*/
cnumber transformed_overlap(matrix3x3 W, vector3 w)
{
    int n1, n2, n3;
    real s1, s2, s3, c1, c2, c3;
    cnumber integral = {0,0}, integral_sum;
    number detW;
    vector3 kvector = cur_kvector;
    matrix3x3 invW, Wc;

    if (!curfield || !strchr("db", curfield_type)) {
        mpi_one_fprintf(stderr, "a real-space H- or D-field must be loaded beforehand\n");
        return integral;
    }
    /* TODO: Is special-casing for #ifndef SCALAR_COMPLEX needed? */
    /* TODO: Doesn't seem to work with mpbi */
    
    /* prepare before looping ... */
    n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;

    s1 = geometry_lattice.size.x / n1; /* pixel spacings */
    s2 = geometry_lattice.size.y / n2;
    s3 = geometry_lattice.size.z / n3;
    c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5; /* offsets (negative) */
    c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
    c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

    detW = matrix3x3_determinant(W);
    if (fabs(fabs(detW) - 1.0) > 1e-12) {
        mpi_one_fprintf(stderr, "valid symmetry operations {W|w} must have |det(W)| = 1.0\n");
        return integral;
    }
    invW = matrix3x3_inverse(W);
    /* W is specified in the lattice basis, but field *vectors* evaluated in real-space are
       in a Cartesian basis: when we transform the vector components, we must account for this
       difference. We transform W to a Cartesian basis Wc=RWR⁻¹ (R=geometry_lattice.basis, a
       matrix w/ columns of Cartesian basis vectors) and then use Wc to transform vector fields */
    Wc = matrix3x3_mult(matrix3x3_mult(geometry_lattice.basis, W),
                        matrix3x3_inverse(geometry_lattice.basis));

    kvector.x *= TWOPI/geometry_lattice.size.x; /* hoist these rescalings outside the  */
    kvector.y *= TWOPI/geometry_lattice.size.y; /* loop; might be that licm takes care */
    kvector.z *= TWOPI/geometry_lattice.size.z; /* of it - but better to be sure       */

    /* Loop over coordinates (introduces int vars i1, i2, i3, xyz_index) */
    LOOP_XYZ(mdata) { /* implies two opening braces '{{' */
        vector3 p, pt;
        scalar_complex F[3], Ft[3], Ftemp[3], integrand, phase;
        double deltaphi;

        /* Current lattice coordinate */
        p.x = i1 * s1 - c1;
        p.y = i2 * s2 - c2;
        p.z = i3 * s3 - c3;

        /* Transformed coordinate pt = {W|w}⁻¹p = W⁻¹(p-w) since {W|w}⁻¹={W⁻¹|-W⁻¹w} */
        pt = matrix3x3_vector3_mult(invW, vector3_minus(p, w));

        /* Bloch field value at transformed coordinate pt: interpolation is needed to ensure
           generality in the case of fractional translations.                                */
        get_bloch_field_point_(Ftemp, pt); /* assigns Ftemp to field at p [excludes exp(ikr) factor] */

        /* Transform the vector components of Ftemp by W to obtain Ft; a bit repetitious but 
           we just write out the matrix-product manually here, for real and imag components.
           Since F, Ftemp, and Ft are in a Cartesian basis, we use Wc instead of W.          */
        Ft[0].re = Wc.c0.x*Ftemp[0].re + Wc.c1.x*Ftemp[1].re + Wc.c2.x*Ftemp[2].re;
        Ft[0].im = Wc.c0.x*Ftemp[0].im + Wc.c1.x*Ftemp[1].im + Wc.c2.x*Ftemp[2].im;
        Ft[1].re = Wc.c0.y*Ftemp[0].re + Wc.c1.y*Ftemp[1].re + Wc.c2.y*Ftemp[2].re;
        Ft[1].im = Wc.c0.y*Ftemp[0].im + Wc.c1.y*Ftemp[1].im + Wc.c2.y*Ftemp[2].im;
        Ft[2].re = Wc.c0.z*Ftemp[0].re + Wc.c1.z*Ftemp[1].re + Wc.c2.z*Ftemp[2].re;
        Ft[2].im = Wc.c0.z*Ftemp[0].im + Wc.c1.z*Ftemp[1].im + Wc.c2.z*Ftemp[2].im;

        /* Get the Bloch field value at current point p [exludes exp(ikr) factor].
           We multiply the input field F (either B or D-field) with μ⁻¹ or ε⁻¹ to get
           H- or E-fields, as the relevant overlap is ⟨F|Ft⟩ = (H|Bt) or (E|Dt), with
           t-postscript denoting a field transformed by {W|w}. Here, we essentially
           adapt some boiler-plate code from compute_field_energy_internal in fields.c       */
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
        CASSIGN_CONJ_MULT(integrand, F[0], Ft[0]);  /* add adjoint(F)*Ft to integrand */
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

    mpi_allreduce(&integral, &integral_sum, 2, number,
                  MPI_DOUBLE, MPI_SUM, mpb_comm);

    if (curfield_type == 'b') {  /* H & B are pseudovectors => transform includes det(W) */
        integral_sum.re *= detW;
        integral_sum.im *= detW;
    }

    return integral_sum;
}

cnumber compute_symmetry(int which_band, matrix3x3 W, vector3 w)
{
    cnumber symval;
    get_bfield(which_band);
    symval = transformed_overlap(W, w);

    return symval;
}

cnumber_list compute_symmetries(matrix3x3 W, vector3 w)
{
    cnumber_list symvals;
    int ib;
    symvals.num_items = num_bands;
    CHK_MALLOC(symvals.items, cnumber, num_bands);

    for (ib = 0; ib < num_bands; ++ib){
        symvals.items[ib] = compute_symmetry(ib+1, W, w);
    }
    return symvals;
}