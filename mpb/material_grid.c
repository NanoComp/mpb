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

/* Get the interpolated value at p from the material grid g.
   p.x/p.y/p.z must be in (-1,2).  This involves a bit more Guile
   internals than I would like, ripped out of scm_uniform_vector_ref
   in libguile/unif.c, but the alternative is a lot of overhead given
   that we know for certain that the material grid is a uniform 3d
   array of g->size double-precision values. */
real material_grid_val(vector3 p, material_grid *g)
{
     if (p.x >= 1.0) p.x -= 1.0; else if (p.x < 0.0) p.x += 1.0;
     if (p.y >= 1.0) p.y -= 1.0; else if (p.y < 0.0) p.y += 1.0;
     if (p.z >= 1.0) p.z -= 1.0; else if (p.z < 0.0) p.z += 1.0;
     CHECK(SCM_ARRAYP(g->matgrid), "bug: matgrid is not an array");
     return linear_interpolate(p.x, p.y, p.z,
			       (double *) SCM_CELL_WORD_1(SCM_ARRAY_V(
							       g->matgrid)),
			       g->size.x, g->size.y, g->size.z, 1);
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

