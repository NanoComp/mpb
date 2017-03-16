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

#include "matrix-smob.h"

#include "mpb.h"

#ifndef HAVE_SCM_MAKE_SMOB_TYPE
/* null mark function, for smobs containing no SCM objects */
static SCM mark_null(SCM obj) { (void) obj; return SCM_BOOL_F; }
#endif

/*************************************************************************/

long scm_tc16_smob_evectmatrix = 0;

static SCM evectmatrix_p(SCM obj)
{
     return ctl_convert_boolean_to_scm(EVECTMATRIX_P(obj));
}

static int print_evectmatrix(SCM obj, SCM port, scm_print_state *pstate)
{
     char buf[256];
     evectmatrix *pm = EVECTMATRIX(obj);
     (void) pstate; /* unused argument */

     scm_puts("#<evectmatrix ", port);
     sprintf(buf, "(%dx%d)x%d", pm->N, pm->c, pm->p);
     scm_puts(buf, port);
#ifdef SCALAR_COMPLEX
     scm_puts(" complex", port);
#else
     scm_puts(" real", port);
#endif
     if (pm->localN < pm->N) {
	  sprintf(buf, ", (%dx%d)x%d local", pm->localN, pm->c, pm->p);
	  scm_puts(buf, port);
     }
     scm_putc('>', port);
     return 1;
}

static size_t free_evectmatrix(SCM obj)
{
     evectmatrix *pm = EVECTMATRIX(obj);
     destroy_evectmatrix(*pm);
     free(pm);
     return 0;
}

#define mark_evectmatrix mark_null

/*************************************************************************/

long scm_tc16_smob_sqmatrix = 0;

static SCM sqmatrix_p(SCM obj)
{
     return ctl_convert_boolean_to_scm(SQMATRIX_P(obj));
}

static int print_sqmatrix(SCM obj, SCM port, scm_print_state *pstate)
{
     char buf[256];
     sqmatrix *pm = SQMATRIX(obj);
     (void) pstate; /* unused argument */

     scm_puts("#<sqmatrix ", port);
     sprintf(buf, "%dx%d", pm->p, pm->p);
     scm_puts(buf, port);
#ifdef SCALAR_COMPLEX
     scm_puts(" complex", port);
#else
     scm_puts(" real", port);
#endif
     if (pm->alloc_p > pm->p) {
	  sprintf(buf, ", %dx%d alloc", pm->alloc_p, pm->alloc_p);
	  scm_puts(buf, port);
     }
     scm_putc('>', port);
     return 1;
}

static size_t free_sqmatrix(SCM obj)
{
     sqmatrix *pm = SQMATRIX(obj);
     destroy_sqmatrix(*pm);
     free(pm);
     return 0;
}

#define mark_sqmatrix mark_null

/*************************************************************************/

/* return a Scheme object *copy* of m */
SCM evectmatrix2scm(evectmatrix m)
{
     SCM obj;
     evectmatrix *mp;
     CHK_MALLOC(mp, evectmatrix, 1);
     *mp = create_evectmatrix(m.N, m.c, m.p, m.localN, m.Nstart, m.allocN);
     evectmatrix_copy(*mp, m);
     NEWCELL_SMOB(obj, evectmatrix, mp);
     return obj;
}

/* return a Scheme object *copy* of the given columns of m */
SCM evectmatrix_slice2scm(evectmatrix m, int p_start, int p)
{
     SCM obj;
     evectmatrix *mp;
     CHECK(p_start >= 0 && p_start + p <= m.p && p >= 0,
	   "invalid arguments in evectmatrix_slice2scm");
     CHK_MALLOC(mp, evectmatrix, 1);
     *mp = create_evectmatrix(m.N, m.c, p, m.localN, m.Nstart, m.allocN);
     evectmatrix_copy_slice(*mp, m, 0, p_start, p);
     NEWCELL_SMOB(obj, evectmatrix, mp);
     return obj;
}

SCM sqmatrix2scm(sqmatrix m)
{
     SCM obj;
     sqmatrix *mp;
     CHK_MALLOC(mp, sqmatrix, 1);
     *mp = create_sqmatrix(m.p);
     sqmatrix_copy(*mp, m);
     NEWCELL_SMOB(obj, sqmatrix, mp);
     return obj;
}

/*************************************************************************/

void register_matrix_smobs(void)
{
#ifdef HAVE_SCM_MAKE_SMOB_TYPE
     scm_tc16_smob_evectmatrix = scm_make_smob_type("evectmatrix", 0);
     scm_set_smob_free(scm_tc16_smob_evectmatrix, free_evectmatrix);
     scm_set_smob_print(scm_tc16_smob_evectmatrix, print_evectmatrix);

     scm_tc16_smob_sqmatrix = scm_make_smob_type("sqmatrix", 0);
     scm_set_smob_free(scm_tc16_smob_sqmatrix, free_sqmatrix);
     scm_set_smob_print(scm_tc16_smob_sqmatrix, print_sqmatrix);
#else /* old way to register smobs */
     MAKE_SMOBFUNS(evectmatrix);
     MAKE_SMOBFUNS(sqmatrix);
     REGISTER_SMOBFUNS(evectmatrix);
     REGISTER_SMOBFUNS(sqmatrix);
#endif     

     gh_new_procedure("evectmatrix?", evectmatrix_p, 1, 0, 0);
     gh_new_procedure("sqmatrix?", sqmatrix_p, 1, 0, 0);
}

/*************************************************************************/

sqmatrix *assert_sqmatrix_smob(SCM mo)
{
     sqmatrix *m = SAFE_SQMATRIX(mo);
     CHECK(m, "wrong type argument: expecting sqmatrix");
     return m;
}

evectmatrix *assert_evectmatrix_smob(SCM mo)
{
     evectmatrix *m = SAFE_EVECTMATRIX(mo);
     CHECK(m, "wrong type argument: expecting evectmatrix");
     return m;
}

/*************************************************************************/

integer sqmatrix_size(SCM mo)
{
     sqmatrix *m = assert_sqmatrix_smob(mo);
     return m->p;
}

cnumber sqmatrix_ref(SCM mo, integer i, integer j)
{
     sqmatrix *m = assert_sqmatrix_smob(mo);
     cnumber c;
     CHECK(m && i >= 0 && j >= 0 && i < m->p && j < m->p,
	   "invalid arguments to sqmatrix-ref");
     c.re = SCALAR_RE(m->data[i * m->p + j]);
     c.im = SCALAR_IM(m->data[i * m->p + j]);
     scm_remember_upto_here_1(mo);
     return c;
}

SCM sqmatrix_mult(SCM Ao, SCM Bo)
{
    sqmatrix *A = assert_sqmatrix_smob(Ao);
    sqmatrix *B = assert_sqmatrix_smob(Bo);
    sqmatrix C;
    SCM obj;
    CHECK(A->p == B->p, "only equal-size matrices can be multiplied");
    C = create_sqmatrix(A->p);
    sqmatrix_AeBC(C, *A, 0, *B, 0);
    obj = sqmatrix2scm(C);
    destroy_sqmatrix(C);
    scm_remember_upto_here_2(Ao, Bo);
    return obj;
}

SCM sqmatrix_diagm(cnumber_list d)
{
    int i, p = d.num_items;
    SCM obj;
    sqmatrix D = create_sqmatrix(p);
    for (i = 0; i < p*p; ++i)
        CASSIGN_ZERO(D.data[i]);
    for (i = 0; i < p; ++i)
        CASSIGN_SCALAR(D.data[i*p+i], d.items[i].re, d.items[i].im);
    obj = sqmatrix2scm(D);
    destroy_sqmatrix(D);
    return obj;
}

cnumber_list sqmatrix_eigvals(SCM Ao)
{
    sqmatrix *A = assert_sqmatrix_smob(Ao);
    cnumber_list eigvals;
    eigvals.num_items = A->p;
    CHK_MALLOC(eigvals.items, cnumber, A->p);
    sqmatrix_eigenvalues(*A, (scalar_complex *) eigvals.items);
    scm_remember_upto_here_1(Ao);
    return eigvals;
}

/*************************************************************************/

SCM get_eigenvectors(integer b_start, integer num_bands)
{
     CHECK(mdata, "init-params must be called before get-eigenvectors");

     return evectmatrix_slice2scm(H, b_start - 1, num_bands);
}

void set_eigenvectors(SCM mo, integer b_start)
{
     evectmatrix *m = assert_evectmatrix_smob(mo);
     CHECK(mdata, "init-params must be called before set-eigenvectors");

     evectmatrix_copy_slice(H, *m, b_start - 1, 0, m->p);
     curfield_reset();
     scm_remember_upto_here_1(mo);
}

SCM dot_eigenvectors(SCM mo, integer b_start)
{
     evectmatrix *m = assert_evectmatrix_smob(mo);
     sqmatrix U;
     SCM obj;
     int final_band = b_start-1 + m->p;

     CHECK(mdata, "init-params must be called before dot-eigenvectors");
     CHECK(final_band <= num_bands, "not enough bands in dot-eigenvectors");

     U = create_sqmatrix(m->p);
     if (mdata->mu_inv == NULL) {
         sqmatrix S = create_sqmatrix(m->p);
         evectmatrix_XtY_slice(U, *m, H, 0, b_start - 1, m->p, S);
         destroy_sqmatrix(S);
     }
     else {
         /* ...we have to do this in blocks of eigensolver_block_size since
            the work matrix W[0] may not have enough space to do it at once. */
         int ib;
         sqmatrix S1 = create_sqmatrix(m->p);
         sqmatrix S2 = create_sqmatrix(m->p);

         for (ib = b_start-1; ib < final_band; ib += Hblock.alloc_p) {
             if (ib + Hblock.alloc_p > final_band) {
                 maxwell_set_num_bands(mdata, final_band - ib);
                 evectmatrix_resize(&Hblock, final_band - ib, 0);
             }
             maxwell_compute_H_from_B(mdata, H, Hblock,
                                      (scalar_complex *) mdata->fft_data,
                                      ib, 0, Hblock.p);
             
             evectmatrix_XtY_slice2(U, *m, Hblock, 0, 0, m->p, Hblock.p,
                                    ib-(b_start-1), S1, S2);
         }

         /* Reset scratch matrix sizes: */
         evectmatrix_resize(&Hblock, Hblock.alloc_p, 0);
         maxwell_set_num_bands(mdata, Hblock.alloc_p);

         destroy_sqmatrix(S2);
         destroy_sqmatrix(S1);
     }
     obj = sqmatrix2scm(U);
     destroy_sqmatrix(U);
     scm_remember_upto_here_1(mo);
     return obj;
}

void scale_eigenvector(integer b, cnumber scale)
{
     scalar s;

     CHECK(mdata, "init-params must be called before scale-eigenvector");
     CHECK(b > 0 && b <= H.p, "invalid band number in scale-eigenvector");

#ifndef SCALAR_COMPLEX
     CHECK(fabs(cnumber_im(scale) * cnumber_re(scale)) < 1e-14,
	   "scale-eigenvector must be called with real argument in mpbi");
#endif     
     ASSIGN_SCALAR(s, cnumber_re(scale), cnumber_im(scale));
     blasglue_scal(H.n, s, H.data + b-1, H.p);
     curfield_reset();
}

void output_eigenvectors(SCM mo, char *filename)
{
     evectmatrix *m = assert_evectmatrix_smob(mo);
     evectmatrixio_writeall_raw(filename, *m);
     curfield_reset();
     scm_remember_upto_here_1(mo);
}

SCM input_eigenvectors(char *filename, integer num_bands)
{
     SCM mo = get_eigenvectors(1, num_bands);
     {
	  evectmatrix *m = assert_evectmatrix_smob(mo);
	  evectmatrixio_readall_raw(filename, *m);
     }
     scm_remember_upto_here_1(mo);     
     return mo;
}

void save_eigenvectors(char *filename)
{
     CHECK(mdata, "init-params must be called before save-eigenvectors");
     printf("Saving eigenvectors to \"%s\"...\n", filename);
     evectmatrixio_writeall_raw(filename, H);
}

void load_eigenvectors(char *filename)
{
     CHECK(mdata, "init-params must be called before load-eigenvectors");
     printf("Loading eigenvectors from \"%s\"...\n", filename);
     evectmatrixio_readall_raw(filename, H);
     curfield_reset();
}

/*************************************************************************/
