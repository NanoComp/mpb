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

#ifndef MPB_H
#define MPB_H

#include <string.h>
#include <maxwell.h>
#include <ctl-io.h>
#include <ctlgeom.h>

/* this integer flag is defined by main.c from libctl, and is
   set when the user runs the program with --verbose */
extern int verbose;

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))

#define TWOPI 6.2831853071795864769252867665590057683943388

/**************************************************************************/

extern void get_epsilon_file_func(const char *fname,
				  maxwell_dielectric_function *func,
				  void **func_data);
extern void destroy_epsilon_file_func_data(void *func_data);

extern real linear_interpolate(real rx, real ry, real rz,
			       real *data, int nx, int ny, int nz, int stride);

/**************************************************************************/

/* global variables for retaining data about the eigenvectors between
   calls from Guile: */

#define MAX_NWORK 10
extern int nwork_alloc;

#define NUM_FFT_BANDS 20 /* max number of bands to FFT at a time */

extern maxwell_data *mdata;
extern maxwell_target_data *mtdata;
extern evectmatrix H, W[MAX_NWORK], Hblock;

extern vector3 cur_kvector;
extern scalar_complex *curfield;
extern int curfield_band;
extern char curfield_type;

extern void curfield_reset(void);

/* R[i]/G[i] are lattice/reciprocal-lattice vectors */
extern real R[3][3], G[3][3];
extern matrix3x3 Rm, Gm; /* same thing, but matrix3x3 */
extern real Vol; /* computational cell volume = |det Rm| */

/* index of current kpoint, for labeling output */
extern int kpoint_index;

/* in fields.c */
extern void compute_field_squared(void);
void get_efield(integer which_band);
real mean_medium_from_matrix(const symmetric_matrix *eps_inv);
void get_bloch_field_point_(scalar_complex *field, vector3 p);

/**************************************************************************/

extern void vector3_to_arr(real arr[3], vector3 v);
extern void matrix3x3_to_arr(real arr[3][3], matrix3x3 m);
extern scalar_complex cnumber2cscalar(cnumber c);
extern cnumber cscalar2cnumber(scalar_complex cs);
extern cvector3 cscalar32cvector3(const scalar_complex *cs);
extern void cvector32cscalar3(scalar_complex *cs, cvector3 v);

/**************************************************************************/
/* in epsilon.c */

extern int no_size_x, no_size_y, no_size_z;
extern geom_box_tree geometry_tree;
extern void reset_epsilon(void);
extern void init_epsilon(void);

/**************************************************************************/
/* material_grid.c */

typedef enum { U_MIN = 0, U_PROD = 1, U_SUM = 2 } material_grid_kinds;
extern real material_grid_val(vector3 p, const material_grid *g);
extern double matgrid_val(vector3 p, geom_box_tree tp, int oi,
			  const material_grid *mg);
material_grid *get_material_grids(geometric_object_list g, int *ngrids);
int material_grids_ntot(const material_grid *grids, int ngrids);
void material_grids_set(const double *u, material_grid *grids, int ngrids);
void material_grids_get(double *u, const material_grid *grids, int ngrids);
void material_grids_addgradient(double *v,
				double scalegrad, int band,
				const material_grid *grids, int ngrids);

/**************************************************************************/

extern const char *parity_string(maxwell_data *d);

#endif /* MPB_H */
