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

int no_size_x = 0, no_size_y = 0, no_size_z = 0;

geom_box_tree geometry_tree = NULL; /* recursive tree of geometry 
				       objects for fast searching */

/**************************************************************************/

typedef struct {
     maxwell_dielectric_function eps_file_func;
     void *eps_file_func_data;
} epsilon_func_data;

static material_type make_dielectric(double epsilon)
{
     material_type m;
     m.which_subclass = DIELECTRIC;
     CHK_MALLOC(m.subclass.dielectric_data, dielectric, 1);
     m.subclass.dielectric_data->epsilon = epsilon;
     return m;
}

static void material_eps(material_type material,
			 symmetric_matrix *eps, symmetric_matrix *eps_inv)
{
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
	 case MATERIAL_GRID:
	      CHECK(0, "invalid use of material-grid");
	      break;
	 case MATERIAL_FUNCTION:
	      CHECK(0, "invalid use of material-function");
	      break;
	 case MATERIAL_TYPE_SELF:
	      CHECK(0, "invalid use of material-type");
	      break;
     }
}

/* Given a position r in the basis of the lattice vectors, return the
   corresponding dielectric tensor and its inverse.  Should be
   called from within init_params (or after init_params), so that the
   geometry input variables will have been read in (for use in libgeom).

   This function is passed to set_maxwell_dielectric to initialize
   the dielectric tensor array for eigenvector calculations. */

static void epsilon_func(symmetric_matrix *eps, symmetric_matrix *eps_inv,
			 const real r[3], void *edata)
{
     epsilon_func_data *d = (epsilon_func_data *) edata;
     geom_box_tree tp;
     int oi;
     material_type material;
     vector3 p;
     boolean inobject;

     /* p needs to be in the lattice *unit* vector basis, while r is
	in the lattice vector basis.  Also, shift origin to the center
        of the grid. */
     p.x = no_size_x ? 0 : (r[0] - 0.5) * geometry_lattice.size.x;
     p.y = no_size_y ? 0 : (r[1] - 0.5) * geometry_lattice.size.y;
     p.z = no_size_z ? 0 : (r[2] - 0.5) * geometry_lattice.size.z;

     /* call search routine from libctl/utils/libgeom/geom.c: 
        (we have to use the lower-level geom_tree_search to
         support material-grid types, which have funny semantics) */
     tp = geom_tree_search(shift_to_unit_cell(p), geometry_tree, &oi);
     if (tp) {
	  inobject = 1;
	  material = tp->objects[oi].o->material;
     }
     else {
	  inobject = 0;
	  material = default_material;
     }

#ifdef DEBUG_GEOMETRY_TREE
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
	  tp = 0; inobject = 0;  /* treat as a "nothing" object */
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

	  /* For a material grid, we interpolate the point (in "object"
	     coordinates) into the grid.  More than that, however,
	     we check if the same point intersects the *same* material grid
	     from multiple objects -- if so, we take the product of
	     the interpolated grid values. */
	  if (material.which_subclass == MATERIAL_GRID) {
	       material_type mat_eps;
	       mat_eps = make_dielectric(
		    matgrid_valprod(p, tp, oi, 
				    material.subclass.material_grid_data)
		    * (material.subclass.material_grid_data->epsilon_max -
		       material.subclass.material_grid_data->epsilon_min) +
		    material.subclass.material_grid_data->epsilon_min);
	       if (destroy_material)
		    material_type_destroy(material);
	       material = mat_eps;
	       destroy_material = 1;
	  }
	       
	  material_eps(material, eps, eps_inv);
	  if (destroy_material)
	       material_type_destroy(material);
     }
}

static int variable_material(int which_subclass)
{
     return (which_subclass == MATERIAL_GRID ||
	     which_subclass == MATERIAL_FUNCTION);
}

static int mean_epsilon_func(symmetric_matrix *meps, 
			     symmetric_matrix *meps_inv,
			     real n[3],
			     real d1, real d2, real d3, real tol,
			     const real r[3], void *edata)
{
     epsilon_func_data *d = (epsilon_func_data *) edata;
     vector3 p;
     const geometric_object *o1 = 0, *o2 = 0;
     vector3 shiftby1, shiftby2, normal;
     geom_box pixel;
     double fill;
     material_type mat1, mat2;
     int id1 = -1, id2 = -1;
     int i;
     const int num_neighbors[3] = { 3, 5, 9 };
     const int neighbors[3][9][3] = { 
	  { {0,0,0}, {-1,0,0}, {1,0,0}, 
	    {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0} },
	  { {0,0,0},
	    {-1,-1,0}, {1,1,0}, {-1,1,0}, {1,-1,0},
	    {0,0,0},{0,0,0},{0,0,0},{0,0,0} },
	  { {0,0,0},
	    {1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
	    {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1} }
     };

     /* p needs to be in the lattice *unit* vector basis, while r is
        in the lattice vector basis.  Also, shift origin to the center
        of the grid. */
     p.x = no_size_x ? 0 : (r[0] - 0.5) * geometry_lattice.size.x;
     p.y = no_size_y ? 0 : (r[1] - 0.5) * geometry_lattice.size.y;
     p.z = no_size_z ? 0 : (r[2] - 0.5) * geometry_lattice.size.z;
     d1 *= no_size_x ? 0 : geometry_lattice.size.x * 0.5;
     d2 *= no_size_y ? 0 : geometry_lattice.size.y * 0.5;
     d3 *= no_size_z ? 0 : geometry_lattice.size.z * 0.5;

     for (i = 0; i < num_neighbors[dimensions - 1]; ++i) {
	  const geometric_object *o;
	  material_type mat;
	  vector3 q, z, shiftby;
	  int id;
	  q.x = p.x + neighbors[dimensions - 1][i][0] * d1;
	  q.y = p.y + neighbors[dimensions - 1][i][1] * d2;
	  q.z = p.z + neighbors[dimensions - 1][i][2] * d3;
	  z = shift_to_unit_cell(q);
	  o = object_of_point_in_tree(z, geometry_tree, &shiftby, &id);
	  shiftby = vector3_plus(shiftby, vector3_minus(q, z));
	  if ((id == id1 && vector3_equal(shiftby, shiftby1)) ||
	      (id == id2 && vector3_equal(shiftby, shiftby2)))
	       continue;
	  mat = (o && o->material.which_subclass != MATERIAL_TYPE_SELF)
	       ? o->material : default_material;
	  if (id1 == -1) {
	       o1 = o;
	       shiftby1 = shiftby;
	       id1 = id;
	       mat1 = mat;
	  }
	  else if (id2 == -1 || ((id >= id1 && id >= id2) &&
				 (id1 == id2 
				  || material_type_equal(&mat1,&mat2)))) {
	       o2 = o;
	       shiftby2 = shiftby;
	       id2 = id;
	       mat2 = mat;
	  }
	  else if (!(id1 < id2 && 
		     (id1 == id || material_type_equal(&mat1,&mat))) &&
		   !(id2 < id1 &&
		     (id2 == id || material_type_equal(&mat2,&mat))))
	       return 0; /* too many nearby objects for analysis */
     }

     CHECK(id1 > -1, "bug in object_of_point_in_tree?");
     if (id2 == -1) { /* only one nearby object/material */
	  id2 = id1;
	  o2 = o1;
	  mat2 = mat1;
	  shiftby2 = shiftby1;
     }

     if ((o1 && variable_material(o1->material.which_subclass)) ||
	 (o2 && variable_material(o2->material.which_subclass)) ||
	 ((variable_material(default_material.which_subclass)
	   || d->eps_file_func)
	  && (!o1 || !o2 ||
	      o1->material.which_subclass == MATERIAL_TYPE_SELF ||
	      o2->material.which_subclass == MATERIAL_TYPE_SELF)))
	  return 0; /* arbitrary material functions are non-analyzable */
	      
     material_eps(mat1, meps, meps_inv);

     /* check for trivial case of only one object/material */
     if (id1 == id2 || material_type_equal(&mat1, &mat2)) { 
	  n[0] = n[1] = n[2] = 0;
	  return 1;
     }

     if (id1 > id2)
	  normal = normal_to_fixed_object(vector3_minus(p, shiftby1), *o1);
     else
	  normal = normal_to_fixed_object(vector3_minus(p, shiftby2), *o2);

     n[0] = no_size_x ? 0 : normal.x / geometry_lattice.size.x;
     n[1] = no_size_y ? 0 : normal.y / geometry_lattice.size.y;
     n[2] = no_size_z ? 0 : normal.z / geometry_lattice.size.z;

     pixel.low.x = p.x - d1;
     pixel.high.x = p.x + d1;
     pixel.low.y = p.y - d2;
     pixel.high.y = p.y + d2;
     pixel.low.z = p.z - d3;
     pixel.high.z = p.z + d3;

     tol = tol > 0.01 ? 0.01 : tol;
     if (id1 > id2) {
	  pixel.low = vector3_minus(pixel.low, shiftby1);
	  pixel.high = vector3_minus(pixel.high, shiftby1);
	  fill = box_overlap_with_object(pixel, *o1, tol, 100/tol);
     }
     else {
	  pixel.low = vector3_minus(pixel.low, shiftby2);
	  pixel.high = vector3_minus(pixel.high, shiftby2);
	  fill = 1 - box_overlap_with_object(pixel, *o2, tol, 100/tol);
     }

     {
	  symmetric_matrix eps2, epsinv2;
	  material_eps(mat2, &eps2, &epsinv2);

	  meps->m00 = fill * (meps->m00 - eps2.m00) + eps2.m00;
	  meps->m11 = fill * (meps->m11 - eps2.m11) + eps2.m11;
	  meps->m22 = fill * (meps->m22 - eps2.m22) + eps2.m22;
#ifdef WITH_HERMITIAN_EPSILON
	  CASSIGN_SCALAR(meps->m01, 
			 fill * (CSCALAR_RE(meps->m01) -
				 CSCALAR_RE(eps2.m01)) + CSCALAR_RE(eps2.m01),
			 fill * (CSCALAR_IM(meps->m01) -
				 CSCALAR_IM(eps2.m01)) + CSCALAR_IM(eps2.m01));
	  CASSIGN_SCALAR(meps->m02, 
			 fill * (CSCALAR_RE(meps->m02) -
				 CSCALAR_RE(eps2.m02)) + CSCALAR_RE(eps2.m02),
			 fill * (CSCALAR_IM(meps->m02) -
				 CSCALAR_IM(eps2.m02)) + CSCALAR_IM(eps2.m02));
	  CASSIGN_SCALAR(meps->m12, 
			 fill * (CSCALAR_RE(meps->m12) -
				 CSCALAR_RE(eps2.m12)) + CSCALAR_RE(eps2.m12),
			 fill * (CSCALAR_IM(meps->m12) -
				 CSCALAR_IM(eps2.m12)) + CSCALAR_IM(eps2.m12));
#else
	  meps->m01 = fill * (meps->m01 - eps2.m01) + eps2.m01;
	  meps->m02 = fill * (meps->m02 - eps2.m02) + eps2.m02;
	  meps->m12 = fill * (meps->m12 - eps2.m12) + eps2.m12;
#endif

	  meps_inv->m00 = fill * (meps_inv->m00 - epsinv2.m00) + epsinv2.m00;
	  meps_inv->m11 = fill * (meps_inv->m11 - epsinv2.m11) + epsinv2.m11;
	  meps_inv->m22 = fill * (meps_inv->m22 - epsinv2.m22) + epsinv2.m22;
#ifdef WITH_HERMITIAN_EPSILON
	  CASSIGN_SCALAR(meps_inv->m01, 
			 fill * (CSCALAR_RE(meps_inv->m01) -
			   CSCALAR_RE(epsinv2.m01)) + CSCALAR_RE(epsinv2.m01),
			 fill * (CSCALAR_IM(meps_inv->m01) -
			   CSCALAR_IM(epsinv2.m01)) + CSCALAR_IM(epsinv2.m01));
	  CASSIGN_SCALAR(meps_inv->m02, 
			 fill * (CSCALAR_RE(meps_inv->m02) -
			   CSCALAR_RE(epsinv2.m02)) + CSCALAR_RE(epsinv2.m02),
			 fill * (CSCALAR_IM(meps_inv->m02) -
			   CSCALAR_IM(epsinv2.m02)) + CSCALAR_IM(epsinv2.m02));
	  CASSIGN_SCALAR(meps_inv->m12, 
			 fill * (CSCALAR_RE(meps_inv->m12) -
			   CSCALAR_RE(epsinv2.m12)) + CSCALAR_RE(epsinv2.m12),
			 fill * (CSCALAR_IM(meps_inv->m12) -
			   CSCALAR_IM(epsinv2.m12)) + CSCALAR_IM(epsinv2.m12));
#else
	  meps_inv->m01 = fill * (meps_inv->m01 - epsinv2.m01) + epsinv2.m01;
	  meps_inv->m02 = fill * (meps_inv->m02 - epsinv2.m02) + epsinv2.m02;
	  meps_inv->m12 = fill * (meps_inv->m12 - epsinv2.m12) + epsinv2.m12;
#endif

     }

     return 1;
}

/**************************************************************************/

/* Initialize the dielectric function of the global mdata structure,
   along with other geometry data.  Should be called from init-params,
   or in general when global input vars have been loaded and mdata
   allocated. */
void init_epsilon(void)
{
     int mesh[3], i;
     int tree_depth, tree_nobjects;
     number no_size; 

     no_size = 2.0 / ctl_get_number("infinity");

     mpi_one_printf("Mesh size is %d.\n", mesh_size);
     mesh[0] = mesh_size;
     mesh[1] = (dimensions > 1) ? mesh_size : 1;
     mesh[2] = (dimensions > 2) ? mesh_size : 1;

     no_size_x = geometry_lattice.size.x <= no_size;
     no_size_y = geometry_lattice.size.y <= no_size || dimensions < 2;
     no_size_z = geometry_lattice.size.z <= no_size || dimensions < 3;

     Rm.c0 = vector3_scale(no_size_x ? 1 : geometry_lattice.size.x, 
			   geometry_lattice.basis.c0);
     Rm.c1 = vector3_scale(no_size_y ? 1 : geometry_lattice.size.y, 
			   geometry_lattice.basis.c1);
     Rm.c2 = vector3_scale(no_size_z ? 1 : geometry_lattice.size.z, 
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

     destroy_geom_box_tree(geometry_tree);  /* destroy any tree from
					       previous runs */
     {
	  geom_box b0;
	  b0.low = vector3_plus(geometry_center,
				vector3_scale(-0.5, geometry_lattice.size));
	  b0.high = vector3_plus(geometry_center,
				 vector3_scale(0.5, geometry_lattice.size));
	  /* pad tree boundaries to allow for sub-pixel averaging */
	  b0.low.x -= geometry_lattice.size.x / mdata->nx;
	  b0.low.y -= geometry_lattice.size.y / mdata->ny;
	  b0.low.z -= geometry_lattice.size.z / mdata->nz;
	  b0.high.x += geometry_lattice.size.x / mdata->nx;
	  b0.high.y += geometry_lattice.size.y / mdata->ny;
	  b0.high.z += geometry_lattice.size.z / mdata->nz;
	  geometry_tree = create_geom_box_tree0(geometry, b0);
     }
     if (verbose && mpi_is_master()) {
	  printf("Geometry object bounding box tree:\n");
	  display_geom_box_tree(5, geometry_tree);
     }
     geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
     mpi_one_printf("Geometric object tree has depth %d and %d object nodes"
	    " (vs. %d actual objects)\n",
	    tree_depth, tree_nobjects, geometry.num_items);

     mpi_one_printf("Initializing dielectric function...\n");
     {
	  epsilon_func_data d;
	  get_epsilon_file_func(epsilon_input_file,
				&d.eps_file_func, &d.eps_file_func_data);
	  set_maxwell_dielectric(mdata, mesh, R, G, 
				 epsilon_func, mean_epsilon_func, &d);
	  destroy_epsilon_file_func_data(d.eps_file_func_data);
     }
}
