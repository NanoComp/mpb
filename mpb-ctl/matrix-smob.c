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

#include <stdio.h>
#include <stdlib.h>

#include "../src/config.h"

#include <check.h>
#include <blasglue.h>
#include <matrices.h>

#include "matrix-smob.h"

#include "mpb.h"

/* null mark function, for smobs containing no SCM objects */
static SCM mark_null(SCM obj) { (void) obj; return SCM_BOOL_F; }

/*************************************************************************/

long scm_tc16_smob_evectmatrix = 0;

static SCM evectmatrix_p(SCM obj)
{
     return gh_bool2scm(EVECTMATRIX_P(obj));
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

#define mark_evectmatrix mark_null

static size_t free_evectmatrix(SCM obj)
{
     evectmatrix *pm = EVECTMATRIX(obj);
     destroy_evectmatrix(*pm);
     return 0;
}

/*************************************************************************/

long scm_tc16_smob_sqmatrix = 0;

static SCM sqmatrix_p(SCM obj)
{
     return gh_bool2scm(SQMATRIX_P(obj));
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

#define mark_sqmatrix mark_null

static size_t free_sqmatrix(SCM obj)
{
     sqmatrix *pm = SQMATRIX(obj);
     destroy_sqmatrix(*pm);
     free(pm);
     return 0;
}

/*************************************************************************/

/* return a Scheme object *copy* of m */
SCM evectmatrix2scm(evectmatrix m)
{
     SCM obj;
     evectmatrix *mp;
     mp = (evectmatrix *) malloc(sizeof(evectmatrix));
     CHECK(mp, "out of memory");
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
     mp = (evectmatrix *) malloc(sizeof(evectmatrix));
     CHECK(mp, "out of memory");
     *mp = create_evectmatrix(m.N, m.c, p, m.localN, m.Nstart, m.allocN);
     evectmatrix_copy_slice(*mp, m, 0, p_start, p);
     NEWCELL_SMOB(obj, evectmatrix, mp);
     return obj;
}

SCM sqmatrix2scm(sqmatrix m)
{
     SCM obj;
     sqmatrix *mp;
     mp = (sqmatrix *) malloc(sizeof(sqmatrix));
     CHECK(mp, "out of memory");
     *mp = create_sqmatrix(m.p);
     sqmatrix_copy(*mp, m);
     NEWCELL_SMOB(obj, sqmatrix, mp);
     return obj;
}

/*************************************************************************/

MAKE_SMOBFUNS(evectmatrix);
MAKE_SMOBFUNS(sqmatrix);

void register_matrix_smobs(void)
{
     REGISTER_SMOBFUNS(evectmatrix);
     REGISTER_SMOBFUNS(sqmatrix);
     gh_new_procedure("evectmatrix?", evectmatrix_p, 1, 0, 0);
     gh_new_procedure("sqmatrix?", sqmatrix_p, 1, 0, 0);
}

/*************************************************************************/

integer sqmatrix_size(SCM mo)
{
     sqmatrix *m = SAFE_SQMATRIX(mo);
     CHECK(m, "invalid argument to sqmatrix-size");
     return m->p;
}

cnumber sqmatrix_ref(SCM mo, integer i, integer j)
{
     sqmatrix *m = SAFE_SQMATRIX(mo);
     cnumber c;
     CHECK(m && i >= 0 && j >= 0 && i < m->p && j < m->p,
	   "invalid arguments to sqmatrix-ref");
     c.re = SCALAR_RE(m->data[i * m->p + j]);
     c.im = SCALAR_IM(m->data[i * m->p + j]);
     return c;
}

/*************************************************************************/

SCM get_eigenvectors(integer b_start, integer num_bands)
{
     CHECK(mdata, "init-params must be called before get-eigenvectors");

     return evectmatrix_slice2scm(H, b_start - 1, num_bands);
}

void set_eigenvectors(SCM mo, integer b_start)
{
     evectmatrix *m = SAFE_EVECTMATRIX(mo);
     CHECK(m, "invalid argument to evectmatrix-size");
     CHECK(mdata, "init-params must be called before set-eigenvectors");

     evectmatrix_copy_slice(H, *m, b_start - 1, 0, m->p);
     curfield_reset();
}

SCM dot_eigenvectors(SCM mo, integer b_start)
{
     evectmatrix *m = SAFE_EVECTMATRIX(mo);
     sqmatrix U, S;
     SCM obj;

     CHECK(m, "invalid argument to evectmatrix-size");
     CHECK(mdata, "init-params must be called before dot-eigenvectors");

     U = create_sqmatrix(m->p);
     S = create_sqmatrix(m->p);
     evectmatrix_XtY_slice(U, *m, H, 0, b_start - 1, m->p, S);
     destroy_sqmatrix(S);
     obj = sqmatrix2scm(U);
     destroy_sqmatrix(U);
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

/*************************************************************************/
