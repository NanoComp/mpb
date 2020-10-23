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

/* THIS FILE IS TO BE #included IN medium.c.  ALSO USED TO GENERATE mu.c */

/**************************************************************************/

static void material_epsilon(material_type material,
			 symmetric_matrix *eps, symmetric_matrix *eps_inv)
{
     switch (material.which_subclass) {
	 case MEDIUM:
	 {
	      real eps_val = material.subclass.medium_data->epsilon;
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
	 case MEDIUM_ANISOTROPIC:
	 {
              medium_anisotropic *d =
		   material.subclass.medium_anisotropic_data;
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
     medium_func_data *d = (medium_func_data *) edata;
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
     tp = geom_tree_search(p = shift_to_unit_cell(p), geometry_tree, &oi);
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
		m2.subclass.medium_data ==
		material.subclass.medium_data,
		"material_of_point & material_of_point_in_tree don't agree!");
     }
#endif

     if (material.which_subclass == MATERIAL_TYPE_SELF) {
	  material = default_material;
	  tp = 0; inobject = 0;  /* treat as a "nothing" object */
     }

     /* if we aren't in any geometric object and we have an epsilon
	file, use that. */
     if (!inobject && d->epsilon_file_func) {
	  d->epsilon_file_func(eps, eps_inv, r, d->epsilon_file_func_data);
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
	       mat_eps = make_epsilon(
		    matgrid_val(p, tp, oi, 
				material.subclass.material_grid_data)
		    * (material.subclass.material_grid_data->epsilon_max -
		       material.subclass.material_grid_data->epsilon_min) +
		    material.subclass.material_grid_data->epsilon_min);
	       if (destroy_material)
		    material_type_destroy(material);
	       material = mat_eps;
	       destroy_material = 1;
	  }
	       
	  material_epsilon(material, eps, eps_inv);
	  if (destroy_material)
	       material_type_destroy(material);
     }
}

static int mean_epsilon_func(symmetric_matrix *meps, 
			     symmetric_matrix *meps_inv,
			     real n[3],
			     real d1, real d2, real d3, real tol,
			     const real r[3], void *edata)
{
     medium_func_data *d = (medium_func_data *) edata;
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

#if 0 /* no averaging */
     epsilon_func(meps, meps_inv, r, edata);
     n[0] = n[1] = n[2] = 0;
     return 1;
#endif

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
	   || d->epsilon_file_func)
	  && (!o1 || !o2 ||
	      o1->material.which_subclass == MATERIAL_TYPE_SELF ||
	      o2->material.which_subclass == MATERIAL_TYPE_SELF)))
	  return 0; /* arbitrary material functions are non-analyzable */
	      
     material_epsilon(mat1, meps, meps_inv);

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

    /* a neighbor point may have zero measure in the voxel (e.g. material terminating
       on a voxel corner); just treat such "boundary points" as homogeneous */
     if (fabs(fill) < tol) {          /* mat1 fills voxel (meps & meps_inv already set) */
        n[0] = n[1] = n[2] = 0;
        return 1;
     }
     else if (fabs(1 - fill) < tol) { /* mat2 fills voxel */
        material_epsilon(mat2, meps, meps_inv);
        n[0] = n[1] = n[2] = 0;
        return 1;
     }

     {
	  symmetric_matrix eps2, epsinv2;
	  symmetric_matrix eps1, delta;
	  double Rot[3][3], norm, n0, n1, n2;
	  material_epsilon(mat2, &eps2, &epsinv2);
	  eps1 = *meps;

	  /* make Cartesian orthonormal frame relative to interface */
	  n0 = R[0][0] * n[0] + R[1][0] * n[1] + R[2][0] * n[2];
	  n1 = R[0][1] * n[0] + R[1][1] * n[1] + R[2][1] * n[2];
	  n2 = R[0][2] * n[0] + R[1][2] * n[1] + R[2][2] * n[2];
	  norm = sqrt(n0*n0 + n1*n1 + n2*n2);
	  if (norm == 0.0)
	       return 0;
	  norm = 1.0 / norm;
	  Rot[0][0] = n0 = n0 * norm;
	  Rot[1][0] = n1 = n1 * norm;
	  Rot[2][0] = n2 = n2 * norm;
	  if (fabs(n0) > 1e-2 || fabs(n1) > 1e-2) { /* (z x n) */
	       Rot[0][2] = n1;
	       Rot[1][2] = -n0;
	       Rot[2][2] = 0;
	  }
	  else { /* n is ~ parallel to z direction, use (x x n) instead */
	       Rot[0][2] = 0;
	       Rot[1][2] = -n2;
	       Rot[2][2] = n1;
	  }
	  { /* normalize second column */
	       double s = Rot[0][2]*Rot[0][2]+Rot[1][2]*Rot[1][2]+Rot[2][2]*Rot[2][2];
	       s = 1.0 / sqrt(s);
	       Rot[0][2] *= s;
	       Rot[1][2] *= s;
	       Rot[2][2] *= s;
	  }
	  /* 1st column is 2nd column x 0th column */
	  Rot[0][1] = Rot[1][2] * Rot[2][0] - Rot[2][2] * Rot[1][0];
	  Rot[1][1] = Rot[2][2] * Rot[0][0] - Rot[0][2] * Rot[2][0];
	  Rot[2][1] = Rot[0][2] * Rot[1][0] - Rot[1][2] * Rot[0][0];

	  /* rotate epsilon tensors to surface parallel/perpendicular axes */
	  maxwell_sym_matrix_rotate(&eps1, &eps1, Rot);
	  maxwell_sym_matrix_rotate(&eps2, &eps2, Rot);

#define AVG (fill * (EXPR(eps1)) + (1-fill) * (EXPR(eps2)))

#define EXPR(eps) (-1 / eps.m00)
	  delta.m00 = AVG;
#undef EXPR
#define EXPR(eps) (eps.m11 - ESCALAR_NORMSQR(eps.m01) / eps.m00)
	  delta.m11 = AVG;
#undef EXPR
#define EXPR(eps) (eps.m22 - ESCALAR_NORMSQR(eps.m02) / eps.m00)
	  delta.m22 = AVG;
#undef EXPR

#define EXPR(eps) (ESCALAR_RE(eps.m01) / eps.m00)
	  ESCALAR_RE(delta.m01) = AVG;
#undef EXPR
#define EXPR(eps) (ESCALAR_RE(eps.m02) / eps.m00)
	  ESCALAR_RE(delta.m02) = AVG;
#undef EXPR
#define EXPR(eps) (ESCALAR_RE(eps.m12) - ESCALAR_MULT_CONJ_RE(eps.m02, eps.m01) / eps.m00)
	  ESCALAR_RE(delta.m12) = AVG;
#undef EXPR

#ifdef WITH_HERMITIAN_EPSILON
#  define EXPR(eps) (ESCALAR_IM(eps.m01) / eps.m00)
	  ESCALAR_IM(delta.m01) = AVG;
#  undef EXPR
#  define EXPR(eps) (ESCALAR_IM(eps.m02) / eps.m00)
	  ESCALAR_IM(delta.m02) = AVG;
#  undef EXPR
#  define EXPR(eps) (ESCALAR_IM(eps.m12) - ESCALAR_MULT_CONJ_IM(eps.m02, eps.m01) / eps.m00)
	  ESCALAR_IM(delta.m12) = AVG;
#  undef EXPR
#endif /* WITH_HERMITIAN_EPSILON */

	  meps->m00 = -1/delta.m00;
	  meps->m11 = delta.m11 - ESCALAR_NORMSQR(delta.m01) / delta.m00;
	  meps->m22 = delta.m22 - ESCALAR_NORMSQR(delta.m02) / delta.m00;
	  ASSIGN_ESCALAR(meps->m01, -ESCALAR_RE(delta.m01)/delta.m00,
			-ESCALAR_IM(delta.m01)/delta.m00);
	  ASSIGN_ESCALAR(meps->m02, -ESCALAR_RE(delta.m02)/delta.m00,
			-ESCALAR_IM(delta.m02)/delta.m00);
	  ASSIGN_ESCALAR(meps->m12, 
			ESCALAR_RE(delta.m12) 
			- ESCALAR_MULT_CONJ_RE(delta.m02, delta.m01)/delta.m00,
			ESCALAR_IM(delta.m12) 
			- ESCALAR_MULT_CONJ_IM(delta.m02, delta.m01)/delta.m00);

#define SWAP(a,b) { double xxx = a; a = b; b = xxx; }	  
	  /* invert rotation matrix = transpose */
	  SWAP(Rot[0][1], Rot[1][0]);
	  SWAP(Rot[0][2], Rot[2][0]);
	  SWAP(Rot[2][1], Rot[1][2]);
	  maxwell_sym_matrix_rotate(meps, meps, Rot); /* rotate back */
#undef SWAP

#  ifdef DEBUG
	  CHECK(negative_epsilon_okp 
		|| maxwell_sym_matrix_positive_definite(meps),
		"negative mean epsilon from Kottke algorithm");
#  endif
     }

     return 1;
}

/****************************************************************************/
/* epsilon output functions (see also fields.c) */

/* get the epsilon function, and compute some statistics */
void get_epsilon(void)
{
     int i, N, last_dim, last_dim_stored, nx, nz, local_y_start;
     real *epsilon;
     real eps_mean = 0, eps_inv_mean = 0, eps_high = -1e20, eps_low = 1e20;
     int fill_count = 0;

     if (!mdata) {
	  mpi_one_fprintf(stderr,
			  "init-params must be called before get-epsilon!\n");
	  return;
     }

     curfield = (scalar_complex *) mdata->fft_data;
     epsilon = (real *) curfield;
     curfield_band = 0;
     curfield_type = epsilon_CURFIELD_TYPE;

     /* get epsilon.  Recall that we actually have an inverse
	dielectric tensor at each point; define an average index by
	the inverse of the average eigenvalue of the 1/eps tensor.
	i.e. 3/(trace 1/eps). */

     N = mdata->fft_output_size;
     last_dim = mdata->last_dim;
     last_dim_stored =
	  mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
     nx = mdata->nx; nz = mdata->nz; local_y_start = mdata->local_y_start;

     for (i = 0; i < N; ++i) {
          if (mdata->eps_inv == NULL)
              epsilon[i] = 1.0;
          else
              epsilon[i] = mean_medium_from_matrix(mdata->eps_inv + i);
	  if (epsilon[i] < eps_low)
	       eps_low = epsilon[i];
	  if (epsilon[i] > eps_high)
	       eps_high = epsilon[i];
	  eps_mean += epsilon[i];
	  eps_inv_mean += 1/epsilon[i];
	  if (epsilon[i] > 1.0001)
	       ++fill_count;
#ifndef SCALAR_COMPLEX
	  /* most points need to be counted twice, by rfftw output symmetry: */
	  {
	       int last_index;
#  ifdef HAVE_MPI
	       if (nz == 1) /* 2d calculation: 1st dim. is truncated one */
		    last_index = i / nx + local_y_start;
	       else
		    last_index = i % last_dim_stored;
#  else
	       last_index = i % last_dim_stored;
#  endif
	       if (last_index != 0 && 2*last_index != last_dim) {
		    eps_mean += epsilon[i];
		    eps_inv_mean += 1/epsilon[i];
		    if (epsilon[i] > 1.0001)
			 ++fill_count;
	       }
	  }
#endif
     }

     mpi_allreduce_1(&eps_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, mpb_comm);
     mpi_allreduce_1(&eps_inv_mean, real, SCALAR_MPI_TYPE,
		     MPI_SUM, mpb_comm);
     mpi_allreduce_1(&eps_low, real, SCALAR_MPI_TYPE,
		     MPI_MIN, mpb_comm);
     mpi_allreduce_1(&eps_high, real, SCALAR_MPI_TYPE,
		     MPI_MAX, mpb_comm);
     mpi_allreduce_1(&fill_count, int, MPI_INT,
                   MPI_SUM, mpb_comm);
     N = mdata->nx * mdata->ny * mdata->nz;
     eps_mean /= N;
     eps_inv_mean = N/eps_inv_mean;

     mpi_one_printf("epsilon: %g-%g, mean %g, harm. mean %g, "
		    "%g%% > 1, %g%% \"fill\"\n",
		    eps_low, eps_high, eps_mean, eps_inv_mean,
		    (100.0 * fill_count) / N, 
		    eps_high == eps_low ? 100.0 :
		    100.0 * (eps_mean-eps_low) / (eps_high-eps_low));
}

/* get the specified component of the dielectric tensor,
   or the inverse tensor if inv != 0 */
void get_epsilon_tensor(int c1, int c2, int imag, int inv)
{
     int i, N;
     real *epsilon;
     int conj = 0, offset = 0;

     curfield_type = '-'; /* only used internally, for now */
     epsilon = (real *) mdata->fft_data;
     N = mdata->fft_output_size;

     switch (c1 * 3 + c2) {
	 case 0:
	      offset = offsetof(symmetric_matrix, m00);
	      break;
	 case 1:
	      offset = offsetof(symmetric_matrix, m01);
	      break;
	 case 2:
	      offset = offsetof(symmetric_matrix, m02);
	      break;
	 case 3:
	      offset = offsetof(symmetric_matrix, m01); /* = conj(m10) */
	      conj = imag;
	      break;
	 case 4:
	      offset = offsetof(symmetric_matrix, m11);
	      break;
	 case 5:
	      offset = offsetof(symmetric_matrix, m12);
	      break;
	 case 6:
	      offset = offsetof(symmetric_matrix, m02); /* = conj(m20) */
	      conj = imag;
	      break;
	 case 7:
	      offset = offsetof(symmetric_matrix, m12); /* = conj(m21) */
	      conj = imag;
	      break;
	 case 8:
	      offset = offsetof(symmetric_matrix, m22);
	      break;
     }

#ifdef WITH_HERMITIAN_EPSILON
     if (c1 != c2 && imag)
	  offset += offsetof(scalar_complex, im);
#endif

     for (i = 0; i < N; ++i) {
	  if (inv) {
	       epsilon[i] = 
		    *((real *) (((char *) &mdata->eps_inv[i]) + offset));
	  }
	  else {
	       symmetric_matrix eps;
	       maxwell_sym_matrix_invert(&eps, &mdata->eps_inv[i]);
	       epsilon[i] = *((real *) (((char *) &eps) + offset));
	  }
	  if (conj)
	       epsilon[i] = -epsilon[i];
     }     
}

