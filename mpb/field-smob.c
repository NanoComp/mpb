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

#include "config.h"

#include <check.h>
#include <mpiglue.h>
#include <mpi_utils.h>

#include "field-smob.h"

#include "mpb.h"

#ifndef HAVE_SCM_MAKE_SMOB_TYPE
/* null mark function, for smobs containing no SCM objects */
static SCM mark_null(SCM obj) { (void) obj; return SCM_BOOL_F; }
#endif

/*************************************************************************/

long scm_tc16_smob_field_smob = 0;

static SCM field_p(SCM obj)
{
     return ctl_convert_boolean_to_scm(FIELD_P(obj));
}

static SCM rscalar_field_p(SCM obj)
{
     return ctl_convert_boolean_to_scm(RSCALAR_FIELD_P(obj));
}

static SCM cscalar_field_p(SCM obj)
{
     return ctl_convert_boolean_to_scm(CSCALAR_FIELD_P(obj));
}

static SCM cvector_field_p(SCM obj)
{
     return ctl_convert_boolean_to_scm(CVECTOR_FIELD_P(obj));
}

static int print_field_smob(SCM obj, SCM port, scm_print_state *pstate)
{
     char buf[256];
     field_smob *pf = FIELD(obj);
     (void) pstate; /* unused argument */

     scm_puts("#<field ", port);
     sprintf(buf, "%dx%dx%d", pf->nx, pf->ny, pf->nz);
     scm_puts(buf, port);
     switch (pf->type) {
	 case RSCALAR_FIELD_SMOB:
	      scm_puts(" real scalar field", port);
	      break;
	 case CSCALAR_FIELD_SMOB:
	      scm_puts(" complex scalar field", port);
	      break;
	 case CVECTOR_FIELD_SMOB:
	      scm_puts(" complex vector field", port);
	      break;
     }
     if (pf->local_ny < pf->ny) {
	  sprintf(buf, ", y=%d-%d local",
		  pf->local_y_start, pf->local_y_start + pf->local_ny - 1);
	  scm_puts(buf, port);
     }
     scm_putc('>', port);
     return 1;
}

static size_t free_field_smob(SCM obj)
{
     field_smob *pf = FIELD(obj);
     free(pf->f.rs);
     free(pf);
     return 0;
}

#define mark_field_smob mark_null

SCM field2scm(field_smob *pf)
{
     SCM obj;
     NEWCELL_SMOB(obj, field_smob, pf);
     return obj;
}

/*************************************************************************/

void register_field_smobs(void)
{
#ifdef HAVE_SCM_MAKE_SMOB_TYPE
     scm_tc16_smob_field_smob = scm_make_smob_type("field", 0);
     scm_set_smob_free(scm_tc16_smob_field_smob, free_field_smob);
     scm_set_smob_print(scm_tc16_smob_field_smob, print_field_smob);
#else /* old way to register smobs */
     MAKE_SMOBFUNS(field_smob);
     REGISTER_SMOBFUNS(field_smob);
#endif

     gh_new_procedure("field?", field_p, 1, 0, 0);
     gh_new_procedure("rscalar-field?", rscalar_field_p, 1, 0, 0);
     gh_new_procedure("cscalar-field?", cscalar_field_p, 1, 0, 0);
     gh_new_procedure("cvector-field?", cvector_field_p, 1, 0, 0);
}

/*************************************************************************/

static field_smob curfield_smob;

field_smob *update_curfield_smob(void)
{
     CHECK(mdata, "init-params must be called before manipulating fields");
     curfield_smob.nx = mdata->nx;
     curfield_smob.ny = mdata->ny;
     curfield_smob.nz = mdata->nz;
     curfield_smob.N = mdata->fft_output_size;
     curfield_smob.local_ny = mdata->local_ny;
     curfield_smob.local_y_start = mdata->local_y_start;
     curfield_smob.last_dim = mdata->last_dim;
     curfield_smob.last_dim_size = mdata->last_dim_size;
     curfield_smob.other_dims = mdata->other_dims;
     curfield_smob.type_char = curfield_type;
     if (strchr("dhbecv", curfield_type)) { /* complex vector field */
	  curfield_smob.type = CVECTOR_FIELD_SMOB;
	  curfield_smob.f.cv = curfield;
     }
     else if (strchr("DHBnR", curfield_type)) { /* real scalar field */
	  curfield_smob.type = RSCALAR_FIELD_SMOB;
	  curfield_smob.f.rs = (real *) curfield;
     }
     else if (strchr("C", curfield_type)) { /* complex scalar field */
	  curfield_smob.type = CSCALAR_FIELD_SMOB;
	  curfield_smob.f.cs = curfield;
     }
     else {
	  curfield_smob.type = RSCALAR_FIELD_SMOB; /* arbitrary */
	  curfield_smob.f.rs = (real *) curfield;
	  if (!curfield_smob.f.rs)
	       curfield_smob.f.rs = (real *) mdata->fft_data;
	  return 0;
     }
     return &curfield_smob;
}

static void update_curfield(field_smob *pf)
{
     if (pf == &curfield_smob) {
	  curfield_type = curfield_smob.type_char;
	  curfield = curfield_smob.f.cv;
     }
}

boolean cur_fieldp(SCM obj)
{
     if (SCM_NIMP(obj) && SCM_SYMBOLP(obj)) {
	  char *s = ctl_symbol2newstr(obj);
	  int ret = !strcmp(s, "cur-field");
	  free(s);
	  return ret;
     }
     return 0;
}

/*************************************************************************/

field_smob *assert_field_smob(SCM fo)
{
     field_smob *f = SAFE_FIELD(fo);
     CHECK(f, "wrong type argument: expecting field");
     return f;
}

/*************************************************************************/

SCM rscalar_field_make(SCM f0)
{
     int i;
     field_smob *pf;
     field_smob *pf0 = assert_field_smob(f0);
     CHK_MALLOC(pf, field_smob, 1);
     *pf = *pf0;
     pf->type = RSCALAR_FIELD_SMOB;
     pf->type_char = 'R';
     CHK_MALLOC(pf->f.rs, real, pf->N);
     for (i = 0; i < pf->N; ++i)
	  pf->f.rs[i] = 0.0;
     scm_remember_upto_here_1(f0);
     return field2scm(pf);
}

SCM cscalar_field_make(SCM f0)
{
     int i;
     field_smob *pf;
     field_smob *pf0 = assert_field_smob(f0);
     CHK_MALLOC(pf, field_smob, 1);
     *pf = *pf0;
     pf->type = CSCALAR_FIELD_SMOB;
     pf->type_char = 'C';
     CHK_MALLOC(pf->f.cs, scalar_complex, pf->N);
     for (i = 0; i < pf->N; ++i)
	  CASSIGN_ZERO(pf->f.cs[i]);
     scm_remember_upto_here_1(f0);
     return field2scm(pf);
}

SCM cvector_field_make(SCM f0)
{
     int i;
     field_smob *pf;
     field_smob *pf0 = assert_field_smob(f0);
     CHECK(mdata, "init-params must be called before cvector-field-make");
     CHK_MALLOC(pf, field_smob, 1);
     *pf = *pf0;
     pf->type = CVECTOR_FIELD_SMOB;
     pf->type_char = 'c';
     CHK_MALLOC(pf->f.cv, scalar_complex, 3 * pf->N);
     for (i = 0; i < pf->N * 3; ++i)
	  CASSIGN_ZERO(pf->f.cv[i]);
     scm_remember_upto_here_1(f0);
     return field2scm(pf);
}

void cvector_field_nonblochB(SCM f)
{
     field_smob *pf = assert_field_smob(f);
     pf->type_char = 'v';
     update_curfield(pf);
     scm_remember_upto_here_1(f);
}

SCM field_make(SCM f0)
{
     field_smob *pf0 = assert_field_smob(f0);
     switch (pf0->type) {
	 case RSCALAR_FIELD_SMOB:
	      return rscalar_field_make(f0);
	 case CSCALAR_FIELD_SMOB:
	      return cscalar_field_make(f0);
	 case CVECTOR_FIELD_SMOB:
	      return cvector_field_make(f0);
     }
     scm_remember_upto_here_1(f0);
     return SCM_UNDEFINED;
}

static boolean fields_conform(field_smob *f1, field_smob *f2)
{
#define EQF(field) (f1->field == f2->field)
     return (EQF(nx) && EQF(ny) && EQF(nz) &&
	     EQF(N) && EQF(local_ny) && EQF(local_y_start) &&
	     EQF(last_dim) && EQF(last_dim_size) && EQF(other_dims));
#undef EQF
}

boolean fields_conformp(SCM f1o, SCM f2o)
{
     field_smob *f1 = assert_field_smob(f1o);
     field_smob *f2 = assert_field_smob(f2o);
     boolean ret = fields_conform(f1, f2);
     scm_remember_upto_here_1(f1o);
     scm_remember_upto_here_1(f2o);
     return ret;
}

static void field_set(field_smob *fd, field_smob *fs)
{
     int i;

     CHECK(fd->type == fs->type && fields_conform(fd, fs),
	   "fields for field-set! must conform");
     switch (fs->type) {
         case RSCALAR_FIELD_SMOB:
	      CHECK(fs->type_char != '-', "must load field for field-set!");
	      for (i = 0; i < fs->N; ++i)
		   fd->f.rs[i] = fs->f.rs[i];
	      break;
         case CSCALAR_FIELD_SMOB:
	      CHECK(fs->type_char != '-', "must load field for field-set!");
	      for (i = 0; i < fs->N; ++i)
		   fd->f.cs[i] = fs->f.cs[i];
	      break;
         case CVECTOR_FIELD_SMOB:
	      CHECK(fs->type_char != '-', "must load field for field-set!");
	      for (i = 0; i < fs->N * 3; ++i)
		   fd->f.cv[i] = fs->f.cv[i];
	      break;
     }
     fd->type_char = fs->type_char;
     update_curfield(fd);
}

void field_setB(SCM dest, SCM src)
{
     field_smob *fd = assert_field_smob(dest);
     field_smob *fs = assert_field_smob(src);
     field_set(fd, fs);
     scm_remember_upto_here_1(dest);
     scm_remember_upto_here_1(src);
}

void field_load(SCM src)
{
     field_smob *fs = assert_field_smob(src);
     CHECK(mdata, "init-params must be called before field-load");
     update_curfield_smob();
     CHECK(fields_conform(fs, &curfield_smob),
	   "argument for field-load must conform to current size");
     curfield_smob.type = fs->type;
     field_set(&curfield_smob, fs);
     scm_remember_upto_here_1(src);
}

void field_mapLB(SCM dest, function f, SCM_list src)
{
     field_smob *pd = assert_field_smob(dest);
     field_smob **ps;
     int i, j;
     CHK_MALLOC(ps, field_smob *, src.num_items);
     for (j = 0; j < src.num_items; ++j) {
	  ps[j] = assert_field_smob(src.items[j]);
	  CHECK(fields_conform(pd, ps[j]),
		"fields for field-map! must conform");
     }
     for (i = 0; i < pd->N; ++i) {
	  list arg_list = SCM_EOL;
	  SCM result;
	  for (j = src.num_items - 1; j >= 0; --j) {
	       SCM item = SCM_EOL;
	       switch (ps[j]->type) {
		   case RSCALAR_FIELD_SMOB:
			item = ctl_convert_number_to_scm(ps[j]->f.rs[i]);
			break;
		   case CSCALAR_FIELD_SMOB:
			item = cnumber2scm(cscalar2cnumber(ps[j]->f.cs[i]));
			break;
		   case CVECTOR_FIELD_SMOB:
			item =
			     cvector32scm(cscalar32cvector3(ps[j]->f.cv+3*i));
			break;
	       }
	       arg_list = gh_cons(item, arg_list);
	  }
	  result = gh_apply(f, arg_list);
	  switch (pd->type) {
	      case RSCALAR_FIELD_SMOB:
		   pd->f.rs[i] = ctl_convert_number_to_c(result);
		   break;
	      case CSCALAR_FIELD_SMOB:
		   pd->f.cs[i] = cnumber2cscalar(scm2cnumber(result));
		   break;
	      case CVECTOR_FIELD_SMOB:
		   cvector32cscalar3(pd->f.cv+3*i, scm2cvector3(result));
		   break;
	  }
     }
     if (src.num_items == 1 && ps[0]->type == pd->type)
	  pd->type_char = ps[0]->type_char;
     else if (src.num_items > 1)
	  switch (pd->type) {
	      case RSCALAR_FIELD_SMOB:
		   pd->type_char = 'R';
		   break;
	      case CSCALAR_FIELD_SMOB:
		   pd->type_char = 'C';
		   break;
	      case CVECTOR_FIELD_SMOB:
		   pd->type_char = 'c';
		   break;
	  }
     scm_remember_upto_here_1(src);
     free(ps);
     update_curfield(pd);
     scm_remember_upto_here_1(dest);
}

/*************************************************************************/

static cvector3 cvector3_conj(cvector3 c)
{
     cvector3 cc;
     cc.x = cnumber_conj(c.x);
     cc.y = cnumber_conj(c.y);
     cc.z = cnumber_conj(c.z);
     return cc;
}

/* Compute the integral of f(r, {fields}) over the cell. */
cnumber integrate_fieldL(function f, SCM_list fields)
{
     int i, j, k, n1, n2, n3, n_other, n_last, rank, last_dim;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
     real s1, s2, s3, c1, c2, c3;
     int ifield;
     field_smob **pf;
     cnumber integral = {0,0};

     CHK_MALLOC(pf, field_smob *, fields.num_items);
     for (ifield = 0; ifield < fields.num_items; ++ifield) {
          pf[ifield] = assert_field_smob(fields.items[ifield]);
          CHECK(fields_conform(pf[0], pf[ifield]),
                "fields for integrate-fields must conform");
     }

     if (fields.num_items > 0) {
	  n1 = pf[0]->nx; n2 = pf[0]->ny; n3 = pf[0]->nz;
	  n_other = pf[0]->other_dims;
	  n_last = pf[0]->last_dim_size
	       / (sizeof(scalar_complex)/sizeof(scalar));
	  last_dim = pf[0]->last_dim;
     }
     else {
	  n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;
	  n_other = mdata->other_dims;
	  n_last = mdata->last_dim_size
	       / (sizeof(scalar_complex)/sizeof(scalar));
	  last_dim = mdata->last_dim;
     }
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

        This was all stolen from maxwell_eps.c...it would be better
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

     if (fields.num_items > 0) {
	  local_n2 = pf[0]->local_ny;
	  local_y_start = pf[0]->local_y_start;
     }
     else {
	  local_n2 = mdata->local_ny;
	  local_y_start = mdata->local_y_start;
     }

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

     if (fields.num_items > 0) {
	  local_n2 = pf[0]->local_ny;
	  local_y_start = pf[0]->local_y_start;
     }
     else {
	  local_n2 = mdata->local_ny;
	  local_y_start = mdata->local_y_start;
     }

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1) {
	  if (fields.num_items > 0)
	       local_n3 = pf[0]->last_dim_size / 2;
	  else
	       local_n3 = mdata->last_dim_size / 2;
     }
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
	       list arg_list = SCM_EOL;
	       cnumber integrand;
	       vector3 p;

	       p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;

	       for (ifield = fields.num_items - 1; ifield >= 0; --ifield) {
		    SCM item = SCM_EOL;
		    switch (pf[ifield]->type) {
			case RSCALAR_FIELD_SMOB:
			     item = ctl_convert_number_to_scm(pf[ifield]->f.rs[index]);
			     break;
			case CSCALAR_FIELD_SMOB:
			     item = cnumber2scm(cscalar2cnumber(
				  pf[ifield]->f.cs[index]));
			     break;
			case CVECTOR_FIELD_SMOB:
                        item = cvector32scm(cscalar32cvector3(
			     pf[ifield]->f.cv+3*index));
                        break;
		    }
		    arg_list = gh_cons(item, arg_list);
	       }
	       arg_list = gh_cons(vector32scm(p), arg_list);
	       integrand = ctl_convert_cnumber_to_c(gh_apply(f, arg_list));
	       integral.re += integrand.re;
	       integral.im += integrand.im;

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
			 arg_list = SCM_EOL;
			 for (ifield = fields.num_items - 1;
			      ifield >= 0; --ifield) {
			      SCM item = SCM_UNDEFINED;
			      switch (pf[ifield]->type) {
				  case RSCALAR_FIELD_SMOB:
				       item = ctl_convert_number_to_scm(
					    pf[ifield]->f.rs[index]);
				       break;
				  case CSCALAR_FIELD_SMOB:
				       item = cnumber2scm(cscalar2cnumber(
					    pf[ifield]->f.cs[index]));
				       break;
				  case CVECTOR_FIELD_SMOB:
				       item = cvector32scm(
					    cvector3_conj(cscalar32cvector3(
						 pf[ifield]->f.cv+3*index)));
				       break;
			      }
			      arg_list = gh_cons(item, arg_list);
			 }
			 arg_list = gh_cons(vector32scm(p), arg_list);
			 integrand =
			      ctl_convert_cnumber_to_c(gh_apply(f, arg_list));
			 integral.re += integrand.re;
			 integral.im += integrand.im;
		    }
	       }
#endif
	  }
     }

     free(pf);

     integral.re *= Vol / (n1 * n2 * n3);
     integral.im *= Vol / (n1 * n2 * n3);
     {
	  cnumber integral_sum;
	  mpi_allreduce(&integral, &integral_sum, 2, number,
			MPI_DOUBLE, MPI_SUM, mpb_comm);
	  return integral_sum;
     }
}
