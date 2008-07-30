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
     if (p.x >= 1.0) p.x -= 1.0; else if (p.x < 0.0) p.x += 1.0;
     if (p.y >= 1.0) p.y -= 1.0; else if (p.y < 0.0) p.y += 1.0;
     if (p.z >= 1.0) p.z -= 1.0; else if (p.z < 0.0) p.z += 1.0;
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

double matgrid_valprod(vector3 p, geom_box_tree tp, int oi,
		       const material_grid *mg)
{
     double u = 1.0;
     if (tp) {
	  do {
	       u *= material_grid_val(
		    to_geom_box_coords(p, &tp->objects[oi]),
		    tp->objects[oi].o->material
		    .subclass.material_grid_data);
	       tp = geom_tree_search_next(p, tp, &oi);
	  } while (tp &&
		   compatible_matgrids(mg, &tp->objects[oi].o->material));
     }
     if (!tp && compatible_matgrids(mg, &default_material)) {
	  p.x = no_size_x ? 0 : p.x / geometry_lattice.size.x + 0.5;
	  p.y = no_size_y ? 0 : p.y / geometry_lattice.size.y + 0.5;
	  p.z = no_size_z ? 0 : p.z / geometry_lattice.size.z + 0.5;
	  u *= material_grid_val(p, 
				 default_material.subclass.material_grid_data);
     }
     return u;
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

int material_grids_ntot(material_grid *grids, int ngrids)
{
     int i, ntot = 0;
     for (i = 0; i < ngrids; ++i)
	  ntot += grids[i].size.x * grids[i].size.y * grids[i].size.z;
     return ntot;
}

void material_grids_set(const double *u, material_grid *grids, int ngrids)
{
     int i, j = 0;
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
        |E|^2 * dV * (product of u's at that point) * (eps_max-eps_min) / u
   where |E|^2 is the field at that point and dV is the volume of the
   pixel.

   This whole process is complicated by the fact that a given material grid
   pixel may appear via multiple objects in the structure (although actually
   this is pretty easy to handle).
*/

static void material_grid_addgradient(double *v, double scalegrad,
				      const material_grid *g,
				      const geometric_object *o,
				      int o_is_default)
{
     int ix,iy,iz, nx, ny, nz;

     /* First, compute the volume element dV for grid pixels in this object. */
     {
	  vector3 p0 = {0,0,0}, e1 = {1,0,0}, e2 = {0,1,0}, e3 = {0,0,1}; 
	  p0 = from_geom_object_coords(p0, *o);
	  e1 = vector3_minus(from_geom_object_coords(e1, *o), p0);
	  e2 = vector3_minus(from_geom_object_coords(e2, *o), p0);
	  e3 = vector3_minus(from_geom_object_coords(e3, *o), p0);

	  /* if we are not a 3d cell, require that at least one
	     of the e vectors be along the "flat" dimension(s), and
	     set these vectors to unit vectors */
	  if (no_size_x) {
	       if (e1.y == 0 && e1.z == 0) e1.x = 1;
	       else if (e2.y == 0 && e2.z == 0) e2.x = 1;
	       else if (e3.y == 0 && e3.z == 0) e3.x = 1;
	       else CHECK(0, "invalid object for material_grid");
	  }
	  if (no_size_y) {
	       if (e1.x == 0 && e1.z == 0) e1.y = 1;
	       else if (e2.x == 0 && e2.z == 0) e2.y = 1;
	       else if (e3.x == 0 && e3.z == 0) e3.y = 1;
	       else CHECK(0, "invalid object for material_grid");
	  }
	  if (no_size_z) {
	       if (e1.y == 0 && e1.x == 0) e1.z = 1;
	       else if (e2.y == 0 && e2.x == 0) e2.z = 1;
	       else if (e3.y == 0 && e3.x == 0) e3.z = 1;
	       else CHECK(0, "invalid object for material_grid");
	  }

	  scalegrad *= fabs(vector3_dot(e1, vector3_cross(e2, e3)))
	       / (g->size.x * g->size.y * g->size.z);
     }

     scalegrad *= (g->epsilon_max - g->epsilon_min);

     nx = g->size.x; ny = g->size.y; nz = g->size.z;
     for (ix = 0; ix < nx; ++ix)
     for (iy = 0; iy < ny; ++iy)
     for (iz = 0; iz < nz; ++iz) {
	  double u = 1.0;
	  geom_box_tree tp;
	  int oi;
	  vector3 p;
	  int found_o = 0;
	  double Esqr;

	  p.x = (ix+0.5) / g->size.x;
	  p.y = (iy+0.5) / g->size.y;
	  p.z = (iz+0.5) / g->size.z;
	  p = from_geom_object_coords(p, *o);

	  tp = geom_tree_search(p, geometry_tree, &oi);
	  if ((!o_is_default && !tp) ||
	      (tp && !compatible_matgrids(g, &tp->objects[oi].o->material)))
	       continue;

	  Esqr = get_energy_point(p);

	  /* multiply values of all grids *except* from the current
	     object, and make sure that we find the current object somewhere */
	  if (tp) {
	       do {
		    if (tp->objects[oi].o == o)
			 found_o = 1;
		    else {
			 u *= material_grid_val(
			      to_geom_box_coords(p, &tp->objects[oi]),
			      tp->objects[oi].o->material
			      .subclass.material_grid_data);
		    }
		    tp = geom_tree_search_next(p, tp, &oi);
	       } while (tp &&
			compatible_matgrids(g, &tp->objects[oi].o->material));
	  }
	  if (!tp && compatible_matgrids(g, &default_material)) {
	       if (o_is_default)
		    found_o = 1;
	       else {
		    p.x = no_size_x ? 0 : p.x / geometry_lattice.size.x + 0.5;
		    p.y = no_size_y ? 0 : p.y / geometry_lattice.size.y + 0.5;
		    p.z = no_size_z ? 0 : p.z / geometry_lattice.size.z + 0.5;
		    u *= material_grid_val(
			 p, default_material.subclass.material_grid_data);
	       }
	  }
	  if (!found_o) continue;

	  v[(ix * ny + iy) * nz + iz] += scalegrad * u * Esqr;
     }
}

void material_grids_addgradient(double *v, double scalegrad, int band,
				const material_grid *grids, int ngrids)
{
     int i;
     CHECK(band <= num_bands, "addgradient called for uncomputed band");
     scalegrad *= -freqs.items[band - 1]/2;
     get_efield(band);
     compute_field_squared();
     for (i = 0; i < ngrids; ++i) {
	  int ntot = grids[i].size.x * grids[i].size.y * grids[i].size.z;
	  int j;

	  for (j = 0; j < geometry.num_items; ++j)
	       if (geometry.items[j].material.which_subclass==MATERIAL_GRID &&
		   material_grid_equal(&grids[i], geometry.items[j].material
				       .subclass.material_grid_data))
		    material_grid_addgradient(v, scalegrad, &grids[i],
					      &geometry.items[j], 0);

	  if (default_material.which_subclass == MATERIAL_GRID &&
	      material_grid_equal(&grids[i], default_material
				  .subclass.material_grid_data)) {
	       vector3 cen = {0,0,0};
	       geometric_object o;
	       o = make_block(default_material, cen,
			      geometry_lattice.basis1,
			      geometry_lattice.basis2,
			      geometry_lattice.basis3,
			      geometry_lattice.size);
	       material_grid_addgradient(v, scalegrad, &grids[i], &o, 1);
	       geometric_object_destroy(o);
	  }
	  
	  v += ntot;
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
