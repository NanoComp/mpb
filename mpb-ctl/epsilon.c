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

#include "../src/config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <maxwell.h>

#include <ctl-io.h>
#include <ctlgeom.h>

#include "mpb.h"

#define USE_GEOMETRY_TREE 1
static geom_box_tree geometry_tree = NULL; /* recursive tree of geometry 
					      objects for fast searching */

/**************************************************************************/

typedef struct {
     maxwell_dielectric_function eps_file_func;
     void *eps_file_func_data;
} epsilon_func_data;

/* Given a position r in the basis of the lattice vectors, return the
   corresponding dielectric tensor and its inverse.  Should be
   called from within init_params (or after init_params), so that the
   geometry input variables will have been read in (for use in libgeom).

   This function is passed to set_maxwell_dielectric to initialize
   the dielectric tensor array for eigenvector calculations. */

static void epsilon_func(symmetric_matrix *eps, symmetric_matrix *eps_inv,
			 real r[3], void *edata)
{
     epsilon_func_data *d = (epsilon_func_data *) edata;
     material_type material;
     vector3 p;
     boolean inobject;

     /* p needs to be in the lattice *unit* vector basis, while r is
	in the lattice vector basis.  Also, shift origin to the center
        of the grid. */
     p.x = (r[0] - 0.5) * geometry_lattice.size.x;
     p.y = dimensions <= 1 ? 0 : (r[1] - 0.5) * geometry_lattice.size.y;
     p.z = dimensions <= 2 ? 0 : (r[2] - 0.5) * geometry_lattice.size.z;

     /* call search routine from libctl/utils/libgeom/geom.c: */
#if USE_GEOMETRY_TREE
     material = material_of_point_in_tree_inobject(p, geometry_tree, 
						   &inobject);
#else
     material = material_of_point_inobject(p, &inobject);
#endif

#if defined(DEBUG_GEOMETRY_TREE) && USE_GEOMETRY_TREE
     {
	  material_type m2 = material_of_point_inobject(p, &inobject);
	  CHECK(m2.which_subclass == material.which_subclass &&
		m2.subclass.dielectric_data ==
		material.subclass.dielectric_data,
		"material_of_point & material_of_point_in_tree don't agree!");
     }
#endif

     if (material.which_subclass == MATERIAL_TYPE_SELF) {
	  material = default_material;
	  inobject = 0;  /* treat as a "nothing" object */
     }

     /* if we aren't in any geometric object and we have an epsilon
	file, use that. */
     if (!inobject && d->eps_file_func) {
	  d->eps_file_func(eps, eps_inv, r, d->eps_file_func_data);
     }
     else {
	  boolean destroy_material = 0;

	  while (material.which_subclass == MATERIAL_FUNCTION) {
	       material_type m;
	       SCM mo;
	       /* material_func is a Scheme function, taking a position
		  vector and returning a material at that point: */
	       mo = gh_call1(material.subclass.
			     material_function_data->material_func,
			     ctl_convert_vector3_to_scm(p));
	       material_type_input(mo, &m);
	       if (destroy_material)
		    material_type_destroy(material);
	       material = m;
	       destroy_material = 1;
	  }
	       
	  switch (material.which_subclass) {
	      case DIELECTRIC:
	      {
		   real eps_val = material.subclass.dielectric_data->epsilon;
		   eps->m00 = eps->m11 = eps->m22 = eps_val;
		   eps_inv->m00 = eps_inv->m11 = eps_inv->m22 = 1.0 / eps_val;
#ifdef WITH_HERMITIAN_EPSILON
		   CASSIGN_ZERO(eps->m01);
		   CASSIGN_ZERO(eps->m02);
		   CASSIGN_ZERO(eps->m12);
		   CASSIGN_ZERO(eps_inv->m01);
		   CASSIGN_ZERO(eps_inv->m02);
		   CASSIGN_ZERO(eps_inv->m12);
#else
		   eps->m01 = eps->m02 = eps->m12 = 0.0;
		   eps_inv->m01 = eps_inv->m02 = eps_inv->m12 = 0.0;
#endif
		   break;
	      }
	      case DIELECTRIC_ANISOTROPIC:
	      {
		   dielectric_anisotropic *d =
			material.subclass.dielectric_anisotropic_data;
		   eps->m00 = d->epsilon_diag.x;
		   eps->m11 = d->epsilon_diag.y;
		   eps->m22 = d->epsilon_diag.z;
#ifdef WITH_HERMITIAN_EPSILON
		   CASSIGN_SCALAR(eps->m01, d->epsilon_offdiag.x.re,
				  d->epsilon_offdiag.x.im +
				  d->epsilon_offdiag_imag.x);
		   CASSIGN_SCALAR(eps->m02, d->epsilon_offdiag.y.re,
				  d->epsilon_offdiag.y.im +
				  d->epsilon_offdiag_imag.y);
		   CASSIGN_SCALAR(eps->m12, d->epsilon_offdiag.z.re,
				  d->epsilon_offdiag.z.im +
				  d->epsilon_offdiag_imag.z);
#else
		   eps->m01 = d->epsilon_offdiag.x.re;
		   eps->m02 = d->epsilon_offdiag.y.re;
		   eps->m12 = d->epsilon_offdiag.z.re;
		   CHECK(vector3_norm(vector3_plus(
			     cvector3_im(d->epsilon_offdiag),
			     d->epsilon_offdiag_imag)) == 0.0,
			 "imaginary epsilon-offdiag is only supported when MPB is configured --with-hermitian-eps");
#endif
		   maxwell_sym_matrix_invert(eps_inv, eps);
		   break;
	      }
	      case MATERIAL_FUNCTION:
		   CHECK(0, "invalid use of material-function");
		   break;
	      case MATERIAL_TYPE_SELF:
		   CHECK(0, "invalid use of material-type");
		   break;
	  }
	  if (destroy_material)
	       material_type_destroy(material);
     }
}

/**************************************************************************/

/* Initialize the dielectric function of the global mdata structure,
   along with other geometry data.  Should be called from init-params,
   or in general when global input vars have been loaded and mdata
   allocated. */
void init_epsilon(void)
{
     int mesh[3], i;
#if USE_GEOMETRY_TREE
     int tree_depth, tree_nobjects;
#endif
     number no_size; 

     no_size = 2.0 / ctl_get_number("infinity");

     mpi_one_printf("Mesh size is %d.\n", mesh_size);
     mesh[0] = mesh_size;
     mesh[1] = (dimensions > 1) ? mesh_size : 1;
     mesh[2] = (dimensions > 2) ? mesh_size : 1;

     Rm.c0 = vector3_scale(geometry_lattice.size.x <= no_size ? 
			   1 : geometry_lattice.size.x, 
			   geometry_lattice.basis.c0);
     Rm.c1 = vector3_scale(geometry_lattice.size.y <= no_size ? 
			   1 : geometry_lattice.size.y, 
			   geometry_lattice.basis.c1);
     Rm.c2 = vector3_scale(geometry_lattice.size.z <= no_size ? 
			   1 : geometry_lattice.size.z, 
			   geometry_lattice.basis.c2);
     mpi_one_printf("Lattice vectors:\n");
     mpi_one_printf("     (%g, %g, %g)\n", Rm.c0.x, Rm.c0.y, Rm.c0.z);  
     mpi_one_printf("     (%g, %g, %g)\n", Rm.c1.x, Rm.c1.y, Rm.c1.z);
     mpi_one_printf("     (%g, %g, %g)\n", Rm.c2.x, Rm.c2.y, Rm.c2.z);
     Vol = fabs(matrix3x3_determinant(Rm));
     mpi_one_printf("Cell volume = %g\n", Vol);
  
     Gm = matrix3x3_inverse(matrix3x3_transpose(Rm));
     mpi_one_printf("Reciprocal lattice vectors (/ 2 pi):\n");
     mpi_one_printf("     (%g, %g, %g)\n", Gm.c0.x, Gm.c0.y, Gm.c0.z);  
     mpi_one_printf("     (%g, %g, %g)\n", Gm.c1.x, Gm.c1.y, Gm.c1.z);
     mpi_one_printf("     (%g, %g, %g)\n", Gm.c2.x, Gm.c2.y, Gm.c2.z);
     
     if (eigensolver_nwork > MAX_NWORK) {
	  mpi_one_printf("(Reducing nwork = %d to maximum: %d.)\n",
		 eigensolver_nwork, MAX_NWORK);
	  eigensolver_nwork = MAX_NWORK;
     }

     matrix3x3_to_arr(R, Rm);
     matrix3x3_to_arr(G, Gm);

     /* we must do this to correct for a non-orthogonal lattice basis: */
     geom_fix_objects();

     mpi_one_printf("Geometric objects:\n");
     if (mpi_is_master())
	  for (i = 0; i < geometry.num_items; ++i) {
	       display_geometric_object_info(5, geometry.items[i]);
	       
	       if (geometry.items[i].material.which_subclass == DIELECTRIC)
		    printf("%*sdielectric constant epsilon = %g\n",
			   5 + 5, "",
			   geometry.items[i].material.
			   subclass.dielectric_data->epsilon);
	  }

#if USE_GEOMETRY_TREE
     destroy_geom_box_tree(geometry_tree);  /* destroy any tree from
					       previous runs */
     geometry_tree =  create_geom_box_tree();
     if (verbose && mpi_is_master()) {
	  printf("Geometry object bounding box tree:\n");
	  display_geom_box_tree(5, geometry_tree);
     }
     geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
     mpi_one_printf("Geometric object tree has depth %d and %d object nodes"
	    " (vs. %d actual objects)\n",
	    tree_depth, tree_nobjects, geometry.num_items);
#endif

     mpi_one_printf("Initializing dielectric function...\n");
     {
	  epsilon_func_data d;
	  get_epsilon_file_func(epsilon_input_file,
				&d.eps_file_func, &d.eps_file_func_data);
	  set_maxwell_dielectric(mdata, mesh, R, G, epsilon_func, &d);
	  destroy_epsilon_file_func_data(d.eps_file_func_data);
     }
}
