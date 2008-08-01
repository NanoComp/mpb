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
#include <math.h>

#include "config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>

#include "mpb.h"

/**************************************************************************/

/* For Guile 1.6, to access this with reasonable efficiency requires
   some ugly code from the guts of libguile/unif.c.  In Guile 1.8,
   they provided a documented way (scm_array_get_handle) to do this,
   but in this case you are also required to call scm_array_handle_release,
   via material_grid_array_release.  In our code, you can only have
   one material_grid array pointer at a time. */
#ifdef HAVE_SCM_ARRAY_GET_HANDLE
static scm_t_array_handle cur_material_grid_array_handle;
#endif
static double *material_grid_array(const material_grid *g)
{
#ifdef HAVE_SCM_ARRAY_GET_HANDLE
     scm_array_get_handle(g->matgrid, &cur_material_grid_array_handle);
     return (double *) scm_array_handle_uniform_writable_elements(
	  &cur_material_grid_array_handle);
#else
     CHECK(SCM_ARRAYP(g->matgrid), "bug: matgrid is not an array");
     return (double *) SCM_CELL_WORD_1(SCM_ARRAY_V(g->matgrid));
#endif
}

static void material_grid_array_release(const material_grid *g)
{
#ifdef HAVE_SCM_ARRAY_GET_HANDLE
     (void) g;
     scm_array_handle_release(&cur_material_grid_array_handle);
#else
     (void) g;
#endif
}

/* Get the interpolated value at p from the material grid g.
   p.x/p.y/p.z must be in (-1,2).  This involves a bit more Guile
   internals than I would like, ripped out of scm_uniform_vector_ref
   in libguile/unif.c, but the alternative is a lot of overhead given
   that we know for certain that the material grid is a uniform 3d
   array of g->size double-precision values. */
real material_grid_val(vector3 p, const material_grid *g)
{
     real val;
     CHECK(SCM_ARRAYP(g->matgrid), "bug: matgrid is not an array");
     val = linear_interpolate(p.x, p.y, p.z, material_grid_array(g),
			      g->size.x, g->size.y, g->size.z, 1);
     material_grid_array_release(g);
     return val;
}

/* Returns true if m is a material grid and has the same epsilon min/max
   as mg */
static int compatible_matgrids(const material_grid *mg,
			       const material_type *m)
{
     return (m->which_subclass == MATERIAL_GRID &&
	     m->subclass.material_grid_data->epsilon_min == mg->epsilon_min &&
	     m->subclass.material_grid_data->epsilon_max == mg->epsilon_max);
}

static int matgrid_val_count = 0; /* cache for gradient calculation */
double matgrid_val(vector3 p, geom_box_tree tp, int oi,
		   const material_grid *mg)
{
     double uprod = 1.0, umin = 1.0, usum = 0.0, u;
     matgrid_val_count = 0;
     CHECK(sizeof(real) == sizeof(double), "material grids require double precision");
     if (tp) {
	  do {
	       u = material_grid_val(
		    to_geom_box_coords(p, &tp->objects[oi]),
		    tp->objects[oi].o->material
		    .subclass.material_grid_data);
	       if (u < umin) umin = u;
	       uprod *= u;
	       usum += u; ++matgrid_val_count;
	       tp = geom_tree_search_next(p, tp, &oi);
	  } while (tp &&
		   compatible_matgrids(mg, &tp->objects[oi].o->material));
     }
     if (!tp && compatible_matgrids(mg, &default_material)) {
	  p.x = no_size_x ? 0 : p.x / geometry_lattice.size.x;
	  p.y = no_size_y ? 0 : p.y / geometry_lattice.size.y;
	  p.z = no_size_z ? 0 : p.z / geometry_lattice.size.z;
	  u = material_grid_val(p, 
				default_material.subclass.material_grid_data);
	  if (u < umin) umin = u;
	  uprod *= u;
	  usum += u; ++matgrid_val_count;
     }
     return (mg->material_grid_kind == U_MIN ? umin
	     : (mg->material_grid_kind == U_PROD ? uprod 
		: usum / matgrid_val_count));
}

/**************************************************************************/

material_grid *get_material_grids(geometric_object_list g, int *ngrids)
{
     int i, nalloc = 0;
     material_grid *grids = 0;
     *ngrids = 0;
     for (i = 0; i < g.num_items; ++i)
	  if (g.items[i].material.which_subclass == MATERIAL_GRID) {
	       int j;
	       for (j = 0; j < *ngrids; ++j)
		    if (material_grid_equal(&grids[j],
					    g.items[i].material.subclass
					    .material_grid_data))
			 break;
	       if (j < *ngrids) continue;
	       if (j >= nalloc) {
		    nalloc = nalloc * 2 + 1;
		    grids = realloc(grids, sizeof(material_grid) * nalloc);
	       }
	       grids[j] = *g.items[i].material.subclass.material_grid_data;
	       ++*ngrids;
	  }
     if (default_material.which_subclass == MATERIAL_GRID) {
	  int j;
	  for (j = 0; j < *ngrids; ++j)
	       if (material_grid_equal(&grids[j],
				       default_material.subclass
				       .material_grid_data))
		    break;
	  if (j == *ngrids) {
	       if (j >= nalloc) {
		    nalloc = nalloc * 2 + 1;
		    grids = realloc(grids, sizeof(material_grid) * nalloc);
	       }
	       grids[j] = *default_material.subclass.material_grid_data;
	       ++*ngrids;
	  }
     }
     return grids;
}

int material_grids_ntot(const material_grid *grids, int ngrids)
{
     int i, ntot = 0;
     for (i = 0; i < ngrids; ++i)
	  ntot += grids[i].size.x * grids[i].size.y * grids[i].size.z;
     return ntot;
}

void material_grids_set(const double *u, material_grid *grids, int ngrids)
{
     int i, j = 0;
     CHECK(sizeof(real) == sizeof(double), "material grids require double precision");
     for (i = 0; i < ngrids; ++i) {
          int ntot = grids[i].size.x * grids[i].size.y * grids[i].size.z;
	  double *a = material_grid_array(&grids[i]);
	  int k;
	  for (k = 0; k < ntot; ++k)
	       a[k] = u[j + k];
	  material_grid_array_release(&grids[i]);
	  j += ntot;
     }
}

void material_grids_get(double *u, const material_grid *grids, int ngrids)
{
     int i, j = 0;
     CHECK(sizeof(real) == sizeof(double), "material grids require double precision");
     for (i = 0; i < ngrids; ++i) {
          int ntot = grids[i].size.x * grids[i].size.y * grids[i].size.z;
	  double *a = material_grid_array(&grids[i]);
	  int k;
	  for (k = 0; k < ntot; ++k)
	       u[j + k] = a[k];
	  material_grid_array_release(&grids[i]);
	  j += ntot;
     }
}

/**************************************************************************/
/* The addgradient function adds to v the gradient, scaled by
   scalegrad, of the frequency of the given band, with respect to
   changes in the material grid values.  This requires that
   solve_kpoint has already been called to solve for the fields (with
   enough bands).   (Note that the band index starts at 1!!)

   By perturbation theory, the change in frequency for a small change
   deps in epsilon is (-omega/2) times the integral of deps |E|^2,
   where E is normalized so that integral eps |E|^2 = 1 (the default
   normalization in MPB).  Thus for a particular "pixel" in the
   material grid with value u, that component of the gradient is:
        |E|^2 * dV * (eps_max-eps_min) * interpolation_weight
   where |E|^2 is the field at that point, dV is the volume of the
   |E|^2 voxel, and interpolation_weight is the weight of that grid
   pixel that contributes to the |E|^2 voxel in the linear interpolation.

   For U_MIN: Where multiple grids overlap, only those grids that contribute
   the minimum u contribute, and for other grids the gradient is zero.
   This unfortunately makes the gradient only piecewise continuous.

   For U_PROD: The gradient is multiplied by the product of u's from 
   overlapping grids, divided by the u from the current grid.  This
   unfortunately makes the gradient zero when two or more u's are zero,
   stalling convergence, although we try to avoid this by making the
   minimum u = 1e-4 instead of 0.

   For U_SUM: The gradient is divided by the number of overlapping grids.
   This doesn't have the property that u=0 in one grid makes the total
   u=0, unfortunately, which is desirable if u=0 indicates "drilled holes".
*/

/* add the weights from linear_interpolate (see the linear_interpolate
   function in epsilon_file.c) to data ... this has to be changed if
   linear_interpolate is changed!! ...also multiply by scaleby
   etc. for different gradient types */
static void add_interpolate_weights(real rx, real ry, real rz, real *data, 
				    int nx, int ny, int nz, int stride,
				    double scaleby,
				    const real *udata, 
				    int ukind, double uval)
{
     int x, y, z, x2, y2, z2;
     real dx, dy, dz, u;

     /* mirror boundary conditions for r just beyond the boundary */
     if (rx < 0.0) rx = -rx; else if (rx > 1.0) rx = 1.0 - rx;
     if (ry < 0.0) ry = -ry; else if (ry > 1.0) ry = 1.0 - ry;
     if (rz < 0.0) rz = -rz; else if (rz > 1.0) rz = 1.0 - rz;

     /* get the point corresponding to r in the epsilon array grid: */
     x = rx * nx; if (x == nx) --x;
     y = ry * ny; if (y == ny) --y;
     z = rz * nz; if (z == nz) --z;

     /* get the difference between (x,y,z) and the actual point
        ... we shift by 0.5 to center the data points in the pixels */
     dx = rx * nx - x - 0.5;
     dy = ry * ny - y - 0.5;
     dz = rz * nz - z - 0.5;

     /* get the other closest point in the grid, with mirror boundaries: */
     x2 = (dx >= 0.0 ? x + 1 : x - 1);
     if (x2 < 0) x2++; else if (x2 == nx) x2--;
     y2 = (dy >= 0.0 ? y + 1 : y - 1);
     if (y2 < 0) y2++; else if (y2 == ny) y2--;
     z2 = (dz >= 0.0 ? z + 1 : z - 1);
     if (z2 < 0) z2++; else if (z2 == nz) z2--;

     /* take abs(d{xyz}) to get weights for {xyz} and {xyz}2: */
     dx = fabs(dx);
     dy = fabs(dy);
     dz = fabs(dz);

     /* define a macro to give us data(x,y,z) on the grid,
	in row-major order (the order used by HDF5): */
#define D(x,y,z) (data[(((x)*ny + (y))*nz + (z)) * stride])
#define U(x,y,z) (udata[(((x)*ny + (y))*nz + (z)) * stride])
     
     u = (((U(x,y,z)*(1.0-dx) + U(x2,y,z)*dx) * (1.0-dy) +
	   (U(x,y2,z)*(1.0-dx) + U(x2,y2,z)*dx) * dy) * (1.0-dz) +
	  ((U(x,y,z2)*(1.0-dx) + U(x2,y,z2)*dx) * (1.0-dy) +
	   (U(x,y2,z2)*(1.0-dx) + U(x2,y2,z2)*dx) * dy) * dz);

     if (ukind == U_MIN && u != uval) return;
     if (ukind == U_PROD) scaleby *= uval / u;
     
     D(x,y,z) += (1.0-dx) * (1.0-dy) * (1.0-dz) * scaleby;
     D(x2,y,z) += dx * (1.0-dy) * (1.0-dz) * scaleby;
     D(x,y2,z) += (1.0-dx) * dy * (1.0-dz) * scaleby;
     D(x2,y2,z) += dx * dy * (1.0-dz) * scaleby;
     D(x,y,z2) += (1.0-dx) * (1.0-dy) * dz * scaleby;
     D(x2,y,z2) += dx * (1.0-dy) * dz * scaleby;
     D(x,y2,z2) += (1.0-dx) * dy * dz * scaleby;
     D(x2,y2,z2) += dx * dy * dz * scaleby;

#undef D
}

static void material_grids_addgradient_point(double *v, 
					     vector3 p, double scalegrad,
					     const material_grid *grids, 
					     int ngrids)
{
     geom_box_tree tp;
     int oi, i;
     material_grid *mg;
     double uval;
     int kind;
     
     tp = geom_tree_search(p, geometry_tree, &oi);
     if (tp && tp->objects[oi].o->material.which_subclass == MATERIAL_GRID)
          mg = tp->objects[oi].o->material.subclass.material_grid_data;
     else if (!tp && default_material.which_subclass == MATERIAL_GRID)
	  mg = default_material.subclass.material_grid_data;
     else
          return; /* no material grids at this point */

     uval = matgrid_val(p, tp, oi, mg);
     scalegrad *= (mg->epsilon_max - mg->epsilon_min);
     if ((kind = mg->material_grid_kind) == U_SUM)
	  scalegrad /= matgrid_val_count;

     if (tp) {
	  do {
	       vector3 pb = to_geom_box_coords(p, &tp->objects[oi]);
	       vector3 sz = tp->objects[oi].o->material
		    .subclass.material_grid_data->size;
	       double *vcur = v, *ucur;
	       for (i = 0; i < ngrids; ++i) {
		    if (material_grid_equal(grids+i,
					    tp->objects[oi].o->material
					    .subclass.material_grid_data))
			 break;
		    else
			 vcur += (int) (grids[i].size.x * grids[i].size.y 
					* grids[i].size.z);
	       }
	       CHECK(i < ngrids, "bug in material_grid_gradient_point");
	       ucur = material_grid_array(grids+i);
	       add_interpolate_weights(pb.x, pb.y, pb.z, 
				       vcur, sz.x, sz.y, sz.z, 1, scalegrad,
				       ucur, kind, uval);
	       material_grid_array_release(grids+i);
	       tp = geom_tree_search_next(p, tp, &oi);
	  } while (tp &&
		   compatible_matgrids(mg, &tp->objects[oi].o->material));
     }
     if (!tp && compatible_matgrids(mg, &default_material)) {
	  vector3 pb;
	  vector3 sz = default_material.subclass.material_grid_data->size;
	  double *vcur = v, *ucur;
	  for (i = 0; i < ngrids; ++i) {
	       if (material_grid_equal(grids+i, default_material
				       .subclass.material_grid_data))
		    break;
	       else
		    vcur += (int) (grids[i].size.x * grids[i].size.y 
				   * grids[i].size.z);
	  }
	  CHECK(i < ngrids, "bug in material_grid_gradient_point");
	  pb.x = no_size_x ? 0 : p.x / geometry_lattice.size.x;
	  pb.y = no_size_y ? 0 : p.y / geometry_lattice.size.y;
	  pb.z = no_size_z ? 0 : p.z / geometry_lattice.size.z;
	  ucur = material_grid_array(grids+i);
	  add_interpolate_weights(pb.x, pb.y, pb.z, 
				  vcur, sz.x, sz.y, sz.z, 1, scalegrad,
				  ucur, kind, uval);
	  material_grid_array_release(grids+i);
     }
}

void material_grids_addgradient(double *v,
				double scalegrad, int band,
				const material_grid *grids, int ngrids)
{
     int i, j, k, n1, n2, n3, n_other, n_last, rank, last_dim;
#ifdef HAVE_MPI
     int local_n2, local_y_start, local_n3;
#endif
     real s1, s2, s3, c1, c2, c3;
     real *Esqr;

     CHECK(band <= num_bands, "addgradient called for uncomputed band");
     if (band) {
	  scalegrad *= -freqs.items[band - 1]/2;
	  get_efield(band);
     }
     compute_field_squared();
     Esqr = (real *) curfield;
     scalegrad *= Vol / H.N;

     n1 = mdata->nx; n2 = mdata->ny; n3 = mdata->nz;
     n_other = mdata->other_dims;
     n_last = mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     last_dim = mdata->last_dim;
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

        This was all stolen from fields.c...it would be better
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

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

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

     local_n2 = mdata->local_ny;
     local_y_start = mdata->local_y_start;

     /* For a real->complex transform, the last dimension is cut in
	half.  For a 2d transform, this is taken into account in local_ny
	already, but for a 3d transform we must compute the new n3: */
     if (n3 > 1)
	  local_n3 = mdata->last_dim_size / 2;
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
	       vector3 p;

	       p.x = i2 * s1 - c1; p.y = j2 * s2 - c2; p.z = k2 * s3 - c3;

	       material_grids_addgradient_point(
		    v, p, Esqr[index]*scalegrad, grids,ngrids);

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
			 
			 material_grids_addgradient_point(
			      v, p, Esqr[index]*scalegrad, grids,ngrids);
		    }
	       }
#endif /* !SCALAR_COMPLEX */

	  }

     }
}

/**************************************************************************/
/* some routines mainly for debugging */

void print_material_grids_gradient(integer band)
{
     int ngrids;
     material_grid *grids = get_material_grids(geometry, &ngrids);
     int i, ntot = material_grids_ntot(grids, ngrids);
     double *grad = (double *) malloc(sizeof(double) * ntot);
     for (i = 0; i < ntot; ++i) grad[i] = 0;
     material_grids_addgradient(grad, 1.0, band, grids, ngrids);
     for (i = 0; i < ntot; ++i)
	  mpi_one_printf(", %g", grad[i]);
     free(grad);
     free(grids);
}

number material_grids_approx_gradient(vector3 kpoint, integer band, 
				      integer iu, number du)
{
     int ngrids;
     material_grid *grids = get_material_grids(geometry, &ngrids);
     int i, ntot = material_grids_ntot(grids, ngrids);
     double *u = (double *) malloc(sizeof(double) * ntot);
     double f0, f1, dfdu;
     solve_kpoint(kpoint);
     f0 = freqs.items[band-1];
     for (i = 0; i < ntot; ++i) u[i] = 0;
     material_grids_addgradient(u, 1.0, band, grids, ngrids);
     dfdu = u[iu];
     material_grids_get(u, grids, ngrids);
     u[iu] += du;
     material_grids_set(u, grids, ngrids);
     reset_epsilon();
     solve_kpoint(kpoint);
     f1 = freqs.items[band-1];
     u[iu] -= du;
     material_grids_set(u, grids, ngrids);
     reset_epsilon();
     mpi_one_printf("approxgrad: ntot=%d, u[%d] = %g -> f_%d = %g, u += %g -> f_%d = %g; df/du = %g vs. analytic %g\n", ntot, iu, u[iu], band, f0, du, band, f1, (f1-f0)/du, dfdu);
     free(u);
     free(grids);
     return (f1 - f0) / du;
}

/**************************************************************************/

#ifdef HAVE_NLOPT_H
#  include <nlopt.h>
#endif

typedef struct {
     boolean do_min;
     vector3_list ks;
     int b1, b2;
     int ngrids;
     material_grid *grids;
     int iter;
     SCM field1, field2;
     double *work;
     int ntot;
} maxgap_func_data;

/* SCM to pass to nlopt for optimization */
static double maxgap_func(int n, const double *u, double *grad, void *data)
{
     maxgap_func_data *d = (maxgap_func_data *) data;
     int i;
     double f1 = 0, f2 = HUGE_VAL, gap, scale, *work = d->work;
     char prefix[256];
     SCM curfield = gh_lookup("cur-field");

     material_grids_set(u, d->grids, d->ngrids);
     reset_epsilon();

     for (i = 0; i < d->ks.num_items; ++i) {
	  solve_kpoint(d->ks.items[i]);
	  if (freqs.items[d->b1 - 1] > f1) {
	       f1 = freqs.items[d->b1- 1];
	       get_efield(d->b1);
	       if (d->field1 == SCM_EOL) d->field1 = field_make(curfield);
	       field_setB(d->field1, curfield);
	  }
	  if (freqs.items[d->b2 - 1] < f2) {
	       f2 = freqs.items[d->b2- 1];
	       get_efield(d->b2);
	       if (d->field2 == SCM_EOL) d->field2 = field_make(curfield);
	       field_setB(d->field2, curfield);
	  }
     }

     gap = (f2 - f1) * 2.0 / (f1 + f2);
     for (i = 0; i < n; ++i) work[i] = 0;

     /* d(-gap)/df1 * (-f1/2)*/
     scale = -f1 * ((f1 + f2) - (f1 - f2)) / ((f1+f2)*(f1+f2));
     if (d->do_min) scale = -scale;
     field_load(d->field1);
     material_grids_addgradient(work, scale, 0, d->grids, d->ngrids);

     /* d(-gap)/df2 * (-f2/2)*/
     scale = -f2 * (-(f1 + f2) - (f1 - f2)) / ((f1+f2)*(f1+f2));
     if (d->do_min) scale = -scale;
     field_load(d->field2);
     material_grids_addgradient(work, scale, 0, d->grids, d->ngrids);

     mpi_allreduce(work, grad, d->ntot, double, MPI_DOUBLE, 
		   MPI_SUM, MPI_COMM_WORLD);

     mpi_one_printf("material-grid-%sgap:, %d, %g, %g, %0.15g\n", 
		    d->do_min ? "min" : "max", ++d->iter, f1, f2, gap);
     
     if (verbose) {
	  get_epsilon();
	  snprintf(prefix, 256, "%sgap-%04d-", 
		   d->do_min ? "min" : "max", d->iter);
	  output_field_to_file(-1, prefix);
     }

     return d->do_min ? gap : -gap;
}

static number material_grids_maxmin_gap(boolean do_min,
					vector3_list kpoints, 
					integer band1, integer band2,
					number func_tol, number eps_tol,
					integer maxeval, number maxtime)
{
     maxgap_func_data d;
     int i, ntot;
     double *u, *lb, *ub, *u_tol, *work, func_min;
     int have_uprod;

     CHECK(band1>0 && band1 <= num_bands && band2>0 && band2 <= num_bands,
	   "invalid band numbers in material-grid-maxgap");
     d.ks = kpoints;
     d.b1 = band1; d.b2 = band2;
     d.grids = get_material_grids(geometry, &d.ngrids);
     d.iter = 0;
     d.field1 = d.field2 = SCM_EOL;
     d.do_min = do_min;

     ntot = material_grids_ntot(d.grids, d.ngrids);
     u = (double *) malloc(sizeof(double) * 5 * ntot);
     lb = u + ntot; ub = lb + ntot; u_tol = ub + ntot; work = u_tol + ntot;
     material_grids_get(u, d.grids, d.ngrids);
     for (i = 0; i < d.ngrids && d.grids[i].material_grid_kind != U_PROD; ++i);
     have_uprod = i < d.ngrids;
     for (i = 0; i < ntot; ++i) {
	  ub[i] = 1;
	  u_tol[i] = eps_tol;
	  /* bound u slightly about 0 for uprod grids, as when u=0
	     the gradient is problematic (especially for multiple u's = 0 */
	  lb[i] = have_uprod ? 1e-4 : 0;
	  if (u[i] < lb[i]) u[i] = lb[i];
     }

     d.work = work;
     d.ntot = ntot;

#if defined(HAVE_NLOPT_H) && defined(HAVE_NLOPT)
 {
     nlopt_result res;
     res = nlopt_minimize(NLOPT_LD_MMA, ntot, maxgap_func, &d,
			  lb, ub, u, &func_min,
			  -HUGE_VAL, func_tol, 0, 0, u_tol, maxeval, maxtime);
     CHECK(res > 0, "failure of nlopt_minimize");
 }
#else
     CHECK(0, "nlopt library is required for material-grid-maxgap");
#endif

     d.field1 = d.field2 = SCM_EOL;

     material_grids_set(u, d.grids, d.ngrids);
     free(u);
     free(d.grids);

     return -func_min;
}

number material_grids_maxgap(vector3_list kpoints, 
			     integer band1, integer band2,
			     number func_tol, number eps_tol,
			     integer maxeval, number maxtime)
{
     return material_grids_maxmin_gap(0, kpoints, band1, band2,
				      func_tol, eps_tol, maxeval, maxtime);
}

number material_grids_mingap(vector3_list kpoints, 
			     integer band1, integer band2,
			     number func_tol, number eps_tol,
			     integer maxeval, number maxtime)
{
     return material_grids_maxmin_gap(1, kpoints, band1, band2,
				      func_tol, eps_tol, maxeval, maxtime);
}

/**************************************************************************/
