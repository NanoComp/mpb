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

int no_size_x = 0, no_size_y = 0, no_size_z = 0;
const int negative_mu_okp = 0; /* mu must always be > 0 */

geom_box_tree geometry_tree = NULL; /* recursive tree of geometry
                                       objects for fast searching */

/**************************************************************************/

typedef struct {
  maxwell_dielectric_function epsilon_file_func;
  void *epsilon_file_func_data;
  maxwell_dielectric_function mu_file_func;
  void *mu_file_func_data;
} medium_func_data;

static material_type make_medium(double epsilon, double mu) {
  material_type m;
  m.which_subclass = MEDIUM;
  CHK_MALLOC(m.subclass.medium_data, medium, 1);
  m.subclass.medium_data->epsilon = epsilon;
  m.subclass.medium_data->mu = mu;
  return m;
}

static material_type make_epsilon(double epsilon) { return make_medium(epsilon, 1.0); }

static material_type make_mu(double mu) { return make_medium(1.0, mu); }

static int variable_material(int which_subclass) {
  return (which_subclass == MATERIAL_GRID || which_subclass == MATERIAL_FUNCTION);
}

/**************************************************************************/

#define epsilon_CURFIELD_TYPE 'n'
#define mu_CURFIELD_TYPE 'm'

#include "epsilon.c"
#include "mu.c"

/**************************************************************************/

static int material_has_mu(material_type m) {
  if (m.which_subclass == MEDIUM) return m.subclass.medium_data->mu != 1;
  if (m.which_subclass == MEDIUM_ANISOTROPIC) {
    medium_anisotropic *a = m.subclass.medium_anisotropic_data;
    return (a->mu_diag.x != 1 || a->mu_diag.y != 1 || a->mu_diag.z != 1 ||
            a->mu_offdiag.x.re != 0 || a->mu_offdiag.y.re != 0 || a->mu_offdiag.z.re != 0 ||
            a->mu_offdiag.x.im + a->mu_offdiag_imag.x != 0 ||
            a->mu_offdiag.y.im + a->mu_offdiag_imag.y != 0 ||
            a->mu_offdiag.z.im + a->mu_offdiag_imag.z != 0);
  }
  if (m.which_subclass == MATERIAL_GRID)
    return (m.subclass.material_grid_data->mu_min != 1 ||
            m.subclass.material_grid_data->mu_max != 1);
  return 0;
}

/* return true if we could potentially have mu != 1 */
static int has_mu(medium_func_data *d) {
  int i;
  if (d->mu_file_func || force_mup || material_has_mu(default_material)) return 1;
  for (i = 0; i < geometry.num_items; ++i)
    if (material_has_mu(geometry.items[i].material)) return 1;
  return 0;
}

/**************************************************************************/

void reset_epsilon(void) {
  medium_func_data d;
  int mesh[3];

  mesh[0] = mesh_size;
  mesh[1] = (dimensions > 1) ? mesh_size : 1;
  mesh[2] = (dimensions > 2) ? mesh_size : 1;

  get_epsilon_file_func(epsilon_input_file, &d.epsilon_file_func, &d.epsilon_file_func_data);
  get_epsilon_file_func(mu_input_file, &d.mu_file_func, &d.mu_file_func_data);
  mpi_one_printf("Initializing epsilon function...\n");
  set_maxwell_dielectric(mdata, mesh, R, G, epsilon_func, mean_epsilon_func, &d);
  if (has_mu(&d)) {
    mpi_one_printf("Initializing mu function...\n");
    set_maxwell_mu(mdata, mesh, R, G, mu_func, mean_mu_func, &d);
  }
  destroy_epsilon_file_func_data(d.epsilon_file_func_data);
  destroy_epsilon_file_func_data(d.mu_file_func_data);
}

/* Initialize the dielectric function of the global mdata structure,
   along with other geometry data.  Should be called from init-params,
   or in general when global input vars have been loaded and mdata
   allocated. */
void init_epsilon(void) {
  int i;
  int tree_depth, tree_nobjects;
  number no_size;

  no_size = 2.0 / ctl_get_number("infinity");

  mpi_one_printf("Mesh size is %d.\n", mesh_size);

  no_size_x = geometry_lattice.size.x <= no_size;
  no_size_y = geometry_lattice.size.y <= no_size || dimensions < 2;
  no_size_z = geometry_lattice.size.z <= no_size || dimensions < 3;

  Rm.c0 = vector3_scale(no_size_x ? 1 : geometry_lattice.size.x, geometry_lattice.basis.c0);
  Rm.c1 = vector3_scale(no_size_y ? 1 : geometry_lattice.size.y, geometry_lattice.basis.c1);
  Rm.c2 = vector3_scale(no_size_z ? 1 : geometry_lattice.size.z, geometry_lattice.basis.c2);
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
    mpi_one_printf("(Reducing nwork = %d to maximum: %d.)\n", eigensolver_nwork, MAX_NWORK);
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

      if (geometry.items[i].material.which_subclass == MEDIUM)
        printf("%*sepsilon = %g, mu = %g\n", 5 + 5, "",
               geometry.items[i].material.subclass.medium_data->epsilon,
               geometry.items[i].material.subclass.medium_data->mu);
    }

  destroy_geom_box_tree(geometry_tree); /* destroy any tree from
                                           previous runs */
  {
    geom_box b0;
    b0.low = vector3_plus(geometry_center, vector3_scale(-0.5, geometry_lattice.size));
    b0.high = vector3_plus(geometry_center, vector3_scale(0.5, geometry_lattice.size));
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

  reset_epsilon();
}
