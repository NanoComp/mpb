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

/* matrixio.c: This file layers a "matrixio" abstraction on top of the
   HDF5 binary i/o interface.  This abstraction should make HDF5 much
   easier to use for our purposes, and could also conceivably allow
   us to replace HDF5 with some other file format (e.g. HDF4). */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "config.h"

#include <check.h>
#include <scalar.h>
#include <mpi_utils.h>

#include "matrixio.h"

#include <mpiglue.h>

/*****************************************************************************/

/* If we have the H5Pset_fapl_mpio function (which is available if HDF5 was
   compiled for MPI), then we can perform collective file i/o operations
   (e.g. all processes call H5Fcreate at the same time to create one file).

   If we don't, however, then we deal with it by having one process
   work with the file at a time: "exclusive" access.  The following
   macro helps us select different bits of code depending upon whether
   this is the case. */

#ifdef HAVE_H5PSET_MPI /* old name for this routine */
#  define H5Pset_fapl_mpio H5Pset_mpi
#  ifndef HAVE_H5PSET_FAPL_MPIO
#    define HAVE_H5PSET_FAPL_MPIO 1
#  endif
#endif

#ifdef HAVE_H5PSET_FAPL_MPIO
#  define IF_EXCLUSIVE(yes,no) no
#else
#  define IF_EXCLUSIVE(yes,no) yes
#endif

/*****************************************************************************/

/* Normally, HDF5 prints out all sorts of error messages, e.g. if a dataset
   can't be found, in addition to returning an error code.  The following
   macro can be wrapped around code to temporarily suppress error messages. */

#define SUPPRESS_HDF5_ERRORS(statements) { \
     H5E_auto_t xxxxx_err_func; \
     void *xxxxx_err_func_data; \
     H5Eget_auto(&xxxxx_err_func, &xxxxx_err_func_data); \
     H5Eset_auto(NULL, NULL); \
     { statements; } \
     H5Eset_auto(xxxxx_err_func, xxxxx_err_func_data); \
}

/*****************************************************************************/

/* Wrappers to write/read an attribute attached to id.  HDF5 attributes
   can *not* be attached to files, in which case we'll write/read it
   as an ordinary dataset.  Ugh. */

static void write_attr(matrixio_id id, matrixio_id_ type_id,
		       matrixio_id_ space_id,
		       const char *name, const void *val)
{
#if defined(HAVE_HDF5)
     hid_t attr_id;

#ifndef HAVE_H5PSET_FAPL_MPIO
     if (!mpi_is_master() && id.parallel)
	  return; /* only one process should add attributes */
#else
     /* otherwise, the operations must be performed collectively */
#endif
     
     if (H5I_FILE == H5Iget_type(id.id)) {
          attr_id = H5Dcreate(id.id, name, type_id, space_id, H5P_DEFAULT);
          CHECK(attr_id >= 0, "error creating HDF attr");
          H5Dwrite(attr_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
          H5Dclose(attr_id);
     }
     else {
	  attr_id = H5Acreate(id.id, name, type_id, space_id, H5P_DEFAULT);
	  CHECK(attr_id >= 0, "error creating HDF attr");
	  H5Awrite(attr_id, type_id, val);
	  H5Aclose(attr_id);
     }
#endif
}

static matrixio_id open_attr(matrixio_id id, matrixio_id_ *type_id,
			     matrixio_id_ *space_id, const char *name)
{
     matrixio_id attr_id;
     attr_id.parallel = id.parallel;
     attr_id.id = -1;
#if defined(HAVE_HDF5)
     if (H5I_FILE == H5Iget_type(id.id)) {
          SUPPRESS_HDF5_ERRORS(attr_id.id = H5Dopen(id.id, name));
	  if (attr_id.id >= 0) {
	       *type_id = H5Dget_type(attr_id.id);
	       *space_id = H5Dget_space(attr_id.id);
	  }
     }
     else {
          SUPPRESS_HDF5_ERRORS(attr_id.id = H5Aopen_name(id.id, name));
	  if (attr_id.id >= 0) {
	       *type_id = H5Aget_type(attr_id.id);
	       *space_id = H5Aget_space(attr_id.id);
	  }
     }

#endif
     return attr_id;
}

static void read_attr(matrixio_id id, matrixio_id attr_id,
		      matrixio_id_ mem_type_id, void *val)
{
#if defined(HAVE_HDF5)
     if (H5I_FILE == H5Iget_type(id.id)) {
	  H5Dread(attr_id.id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
     }
     else {
	  H5Aread(attr_id.id, mem_type_id, val);
     }
#endif
}

static void close_attr(matrixio_id id, matrixio_id attr_id)
{
#if defined(HAVE_HDF5)
     if (H5I_FILE == H5Iget_type(id.id)) {
          H5Dclose(attr_id.id);
     }
     else {
          H5Aclose(attr_id.id);
     }
#endif
}

/*****************************************************************************/

void matrixio_write_string_attr(matrixio_id id, const char *name,
				const char *val)
{
#if defined(HAVE_HDF5)
     hid_t type_id;
     hid_t space_id;

     if (!val || !name || !name[0] || !val[0])
	  return; /* don't try to create empty attributes */
     
     type_id = H5Tcopy(H5T_C_S1);
     H5Tset_size(type_id, strlen(val) + 1);
     space_id = H5Screate(H5S_SCALAR);
     write_attr(id, type_id, space_id, name, val);
     H5Sclose(space_id);
     H5Tclose(type_id);
#endif
}

void matrixio_write_data_attr(matrixio_id id, const char *name,
			      const real *val, int rank, const int *dims)
{
#if defined(HAVE_HDF5)
     hid_t type_id;
     hid_t space_id;
     hsize_t *space_dims;
     int i;

     if (!val || !name || !name[0] || rank < 0 || !dims)
	  return; /* don't try to create empty attributes */
     
#if defined(SCALAR_SINGLE_PREC)
     type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
     type_id = H5T_NATIVE_LDOUBLE;
#else
     type_id = H5T_NATIVE_DOUBLE;
#endif

     if (rank > 0) {
	  CHK_MALLOC(space_dims, hsize_t, rank);
	  for (i = 0; i < rank; ++i)
	       space_dims[i] = dims[i];
	  space_id = H5Screate_simple(rank, space_dims, NULL);
	  free(space_dims);
     }
     else {
	  space_id = H5Screate(H5S_SCALAR);
     }

     write_attr(id, type_id, space_id, name, val);
     H5Sclose(space_id);
#endif
}

char *matrixio_read_string_attr(matrixio_id id, const char *name)
{
#if defined(HAVE_HDF5)
     matrixio_id attr_id;
     hid_t type_id;
     hid_t space_id;
     int len;
     char *s = NULL;

     if (!name || !name[0])
	  return NULL; /* don't try to read empty-named attributes */
     
     attr_id = open_attr(id, &type_id, &space_id, name);
     if (attr_id.id < 0)
	  return NULL;

     if (H5Sget_simple_extent_npoints(space_id) == 1) {
	  len = H5Tget_size(type_id);
	  H5Tclose(type_id);
	  
	  type_id = H5Tcopy(H5T_C_S1);
	  H5Tset_size(type_id, len);
	  
	  CHK_MALLOC(s, char, len);
	  read_attr(id, attr_id, type_id, s);
     }

     H5Tclose(type_id);
     H5Sclose(space_id);
     close_attr(id, attr_id);

     return s;
#else
     return NULL;
#endif
}

real *matrixio_read_data_attr(matrixio_id id, const char *name,
			      int *rank, int max_rank, int *dims)
{
#if defined(HAVE_HDF5)
     matrixio_id attr_id;
     hid_t type_id, mem_type_id, space_id;
     real *d = NULL;

     if (!name || !name[0] || max_rank < 0 || !dims)
	  return NULL; /* don't try to create empty attributes */
     
#if defined(SCALAR_SINGLE_PREC)
     mem_type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
     mem_type_id = H5T_NATIVE_LDOUBLE;
#else
     mem_type_id = H5T_NATIVE_DOUBLE;
#endif

     attr_id = open_attr(id, &type_id, &space_id, name);
     if (attr_id.id < 0)
          return NULL;

     *rank = H5Sget_simple_extent_ndims(space_id);
     if (*rank <= max_rank) {
	  if (*rank > 0) {
	       int i;
	       hsize_t *space_dims, *maxdims;
	       CHK_MALLOC(space_dims, hsize_t, *rank);
	       CHK_MALLOC(maxdims, hsize_t, *rank);
	       H5Sget_simple_extent_dims(space_id, space_dims, maxdims);
	       for (i = 0; i < *rank; ++i)
		    dims[i] = space_dims[i];
	       free(maxdims);
	       free(space_dims);
	  }
	  CHK_MALLOC(d, real, H5Sget_simple_extent_npoints(space_id));
          read_attr(id, attr_id, mem_type_id, d);
     }

     H5Tclose(type_id);
     H5Sclose(space_id);
     close_attr(id, attr_id);

     return d;
#else
     return NULL;
#endif
}

/*****************************************************************************/

#define FNAME_SUFFIX ".h5"  /* standard HDF5 filename suffix */

static char *add_fname_suffix(const char *fname)
{
     int oldlen = strlen(fname);
     int suflen = strlen(FNAME_SUFFIX);
     char *new_fname;

     CHECK(fname, "null filename!");

     CHK_MALLOC(new_fname, char, oldlen + suflen + 1);

     strcpy(new_fname, fname);

     /* only add suffix if it is not already there: */
     if (strstr(new_fname, FNAME_SUFFIX) != new_fname + oldlen - suflen)
	  strcat(new_fname, FNAME_SUFFIX);

     return new_fname;
}

/*****************************************************************************/

#ifndef HAVE_H5PSET_FAPL_MPIO
static int matrixio_critical_section_tag = 0;
#endif

static matrixio_id matrixio_create_(const char *fname, int parallel)
{
#if defined(HAVE_HDF5)
     char *new_fname;
     matrixio_id id;
     hid_t access_props;

     access_props = H5Pcreate (H5P_FILE_ACCESS);
     
#  if defined(HAVE_MPI) && defined(HAVE_H5PSET_FAPL_MPIO)
     if (parallel)
	  CHECK(H5Pset_fapl_mpio(access_props, MPI_COMM_WORLD, MPI_INFO_NULL)
		>= 0, "error initializing MPI file access");
#  endif

     new_fname = add_fname_suffix(fname);

#  ifdef HAVE_H5PSET_FAPL_MPIO
     id.id = H5Fcreate(new_fname, H5F_ACC_TRUNC, H5P_DEFAULT, access_props);
#  else
     if (parallel) mpi_begin_critical_section(matrixio_critical_section_tag);
     if (mpi_is_master() || !parallel)
	  id.id = H5Fcreate(new_fname, H5F_ACC_TRUNC,H5P_DEFAULT,access_props);
     else
	  id.id = H5Fopen(new_fname, H5F_ACC_RDWR, access_props);
#  endif
     id.parallel = parallel;

     CHECK(id.id >= 0, "error creating HDF output file");

     free(new_fname);

     H5Pclose(access_props);

     return id;
#else
     mpi_one_fprintf(stderr,
		     "matrixio: cannot output \"%s\" (compiled without HDF)\n",
		     fname);
     {
	  matrixio_id id = {0,0};
	  return id;
     }
#endif
}

matrixio_id matrixio_create(const char *fname) {
     return matrixio_create_(fname, 1);
}

matrixio_id matrixio_create_serial(const char *fname) {
     return matrixio_create_(fname, 0);
}

static matrixio_id matrixio_open_(const char *fname, int read_only, int parallel)
{
#if defined(HAVE_HDF5)
     char *new_fname;
     matrixio_id id;
     hid_t access_props;

     access_props = H5Pcreate (H5P_FILE_ACCESS);
     
#  if defined(HAVE_MPI) && defined(HAVE_H5PSET_FAPL_MPIO)
     if (parallel)
	  H5Pset_fapl_mpio(access_props, MPI_COMM_WORLD, MPI_INFO_NULL);
#  endif

     new_fname = add_fname_suffix(fname);

     IF_EXCLUSIVE(if (parallel) mpi_begin_critical_section(matrixio_critical_section_tag),0);

     if (read_only)
	  id.id = H5Fopen(new_fname, H5F_ACC_RDONLY, access_props);
     else
	  id.id = H5Fopen(new_fname, H5F_ACC_RDWR, access_props);
     id.parallel = parallel;
     CHECK(id.id >= 0, "error opening HDF input file");

     free(new_fname);

     H5Pclose(access_props);

     return id;
#else
     CHECK(0, "no matrixio implementation is linked");
     {
	  matrixio_id id = {0,0};
	  return id;
     }
#endif
}

void matrixio_close(matrixio_id id)
{
#if defined(HAVE_HDF5)
     CHECK(H5Fclose(id.id) >= 0, "error closing HDF file");
     IF_EXCLUSIVE(if (id.parallel) mpi_end_critical_section(matrixio_critical_section_tag++),0);
#endif
}

matrixio_id matrixio_open(const char *fname, int read_only) {
     return matrixio_open_(fname, read_only, 1);
}

matrixio_id matrixio_open_serial(const char *fname, int read_only) {
     return matrixio_open_(fname, read_only, 0);
}

/*****************************************************************************/

matrixio_id matrixio_create_sub(matrixio_id id, 
				const char *name, const char *description)
{
     matrixio_id sub_id;
     sub_id.id = 0;
     sub_id.parallel = id.parallel;
#if defined(HAVE_HDF5)

#  ifdef HAVE_H5PSET_FAPL_MPIO /* H5Gcreate is collective */
     sub_id.id = H5Gcreate(id.id, name, 0 /* ==> default size */ );
     matrixio_write_string_attr(sub_id, "description", description);
#  else
     /* when running a parallel job, only the master process creates the
	group.  It flushes the group to disk and then the other processes
	open the group.  Is this the right thing to do, or is the
        H5Gcreate function parallel-aware? */

     if (mpi_is_master() || !id.parallel) {
	  sub_id.id = H5Gcreate(id.id, name, 0 /* ==> default size */ );
	  matrixio_write_string_attr(sub_id, "description", description);
	  
	  H5Fflush(sub_id.id, H5F_SCOPE_GLOBAL);

	  IF_EXCLUSIVE(0,if (id.parallel) MPI_Barrier(MPI_COMM_WORLD));
     }
     else {
	  IF_EXCLUSIVE(0,if (id.parallel) MPI_Barrier(MPI_COMM_WORLD));

	  sub_id.id = H5Gopen(id.id, name);
     }
#  endif
#endif
     return sub_id;
}

void matrixio_close_sub(matrixio_id id)
{
#if defined(HAVE_HDF5)
     CHECK(H5Gclose(id.id) >= 0, "error closing HDF group");
#endif
}

/*****************************************************************************/

matrixio_id matrixio_open_dataset(matrixio_id id,
				  const char *name,
				  int rank, const int *dims)
{
     matrixio_id data_id;
     data_id.id = 0;
     data_id.parallel = id.parallel;
#if defined(HAVE_HDF5)
 {
     int i, rank_copy;
     hid_t space_id;
     hsize_t *dims_copy, *maxdims;

     CHECK((data_id.id = H5Dopen(id.id, name)) >= 0, "error in H5Dopen");

     CHECK((space_id = H5Dget_space(data_id.id)) >= 0,
	   "error in H5Dget_space");

     rank_copy = H5Sget_simple_extent_ndims(space_id);
     CHECK(rank == rank_copy, "rank in HDF5 file doesn't match expected rank");
     
     CHK_MALLOC(dims_copy, hsize_t, rank);
     CHK_MALLOC(maxdims, hsize_t, rank);
     H5Sget_simple_extent_dims(space_id, dims_copy, maxdims);
     free(maxdims);
     for (i = 0; i < rank; ++i) {
	  CHECK(dims_copy[i] == dims[i],
		"array size in HDF5 file doesn't match expected size");
     }
     free(dims_copy);

     H5Sclose(space_id);
 }
#endif
     return data_id;
}

/*****************************************************************************/

matrixio_id matrixio_create_dataset(matrixio_id id,
				    const char *name, const char *description,
				    int rank, const int *dims)
{
     matrixio_id data_id;
     data_id.id = 0;
     data_id.parallel = id.parallel;
#if defined(HAVE_HDF5)
 {
     int i;
     hid_t space_id, type_id;
     hsize_t *dims_copy;

     /* delete pre-existing datasets, or we'll have an error; I think
        we can only do this on the master process. (?) */
     if (matrixio_dataset_exists(id, name)) {
#  ifdef HAVE_H5PSET_FAPL_MPIO /* H5Gunlink is collective */
	  matrixio_dataset_delete(id, name);
#  else
	  if (mpi_is_master() || !id.parallel) {
	       matrixio_dataset_delete(id, name);
	       H5Fflush(id.id, H5F_SCOPE_GLOBAL);
	  }
	  IF_EXCLUSIVE(0,if (id.parallel) MPI_Barrier(MPI_COMM_WORLD));
#  endif
     }

     CHECK(rank > 0, "non-positive rank");

     CHK_MALLOC(dims_copy, hsize_t, rank);
     for (i = 0; i < rank; ++i)
          dims_copy[i] = dims[i];

     space_id = H5Screate_simple(rank, dims_copy, NULL);

     free(dims_copy);

#if defined(SCALAR_SINGLE_PREC)
     type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
     type_id = H5T_NATIVE_LDOUBLE;
#else
     type_id = H5T_NATIVE_DOUBLE;
#endif
     
     /* Create the dataset.  Note that, on parallel machines, H5Dcreate
	should do the right thing; it is supposedly a collective operation. */
     IF_EXCLUSIVE(
	  if (mpi_is_master() || !id.parallel)
	       data_id.id = H5Dcreate(id.id,name,type_id,space_id,H5P_DEFAULT);
	  else
	       data_id.id = H5Dopen(id.id, name),
	  data_id.id = H5Dcreate(id.id, name, type_id, space_id, H5P_DEFAULT));

     H5Sclose(space_id);  /* the dataset should have its own copy now */
     
     matrixio_write_string_attr(data_id, "description", description);
 }
#endif
     return data_id;
}

void matrixio_close_dataset(matrixio_id data_id)
{
#if defined(HAVE_HDF5)
     CHECK(H5Dclose(data_id.id) >= 0, "error closing HDF dataset");
#endif
}

int matrixio_dataset_exists(matrixio_id id, const char *name)
{
#if defined(HAVE_HDF5)
     hid_t data_id;
     SUPPRESS_HDF5_ERRORS(data_id = H5Dopen(id.id, name));
     if (data_id >= 0)
	  H5Dclose(data_id);
     return (data_id >= 0);
#else
     return 0;
#endif
}

void matrixio_dataset_delete(matrixio_id id, const char *name)
{
#if defined(HAVE_HDF5)
     H5Gunlink(id.id, name);
#endif
}

/*****************************************************************************/

void matrixio_write_real_data(matrixio_id data_id,
			      const int *local_dims, const int *local_start,
			      int stride,
			      real *data)
{
#if defined(HAVE_HDF5)
     int rank;
     hsize_t *dims, *maxdims;
     hid_t space_id, type_id, mem_space_id;
     start_t *start;
     hsize_t *strides, *count, count_prod;
     int i;
     real *data_copy;
     int data_copy_stride = 1, free_data_copy = 0, do_write = 1;

     /*******************************************************************/
     /* Get dimensions of dataset */
     
     space_id = H5Dget_space(data_id.id);

     rank = H5Sget_simple_extent_ndims(space_id);
     
     CHK_MALLOC(dims, hsize_t, rank);
     CHK_MALLOC(maxdims, hsize_t, rank);

     H5Sget_simple_extent_dims(space_id, dims, maxdims);

     free(maxdims);

#if defined(SCALAR_SINGLE_PREC)
     type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
     type_id = H5T_NATIVE_LDOUBLE;
#else
     type_id = H5T_NATIVE_DOUBLE;
#endif

     /*******************************************************************/
     /* if stride > 1, make a contiguous copy; hdf5 is much faster
	in this case. */

     if (stride > 1) {
	  int N = 1;
	  for (i = 0; i < rank; ++i)
	       N *= local_dims[i];
	  CHK_MALLOC(data_copy, real, N);
	  if (data_copy) {
	       free_data_copy = 1;
	       for (i = 0; i < (N & 3); ++i)
		    data_copy[i] = data[i * stride];
	       for (; i < N; i += 4) {
		    real d0 = data[i * stride];
		    real d1 = data[(i + 1) * stride];
		    real d2 = data[(i + 2) * stride];
		    real d3 = data[(i + 3) * stride];
		    data_copy[i] = d0;
		    data_copy[i+1] = d1;
		    data_copy[i+2] = d2;
		    data_copy[i+3] = d3;
	       }
	       CHECK(i == N, "bug in matrixio copy routine");
	  }
	  else {
	       data_copy = data;
	       data_copy_stride = stride;
	  }
     }
     else
	  data_copy = data;

     /*******************************************************************/
     /* Before we can write the data to the data set, we must define
	the dimensions and "selections" of the arrays to be read & written: */

     CHK_MALLOC(start, start_t, rank);
     CHK_MALLOC(strides, hsize_t, rank);
     CHK_MALLOC(count, hsize_t, rank);

     count_prod = 1;
     for (i = 0; i < rank; ++i) {
	  start[i] = local_start[i];
	  count[i] = local_dims[i];
	  strides[i] = 1;
	  count_prod *= count[i];
     }

     if (count_prod > 0) {
	  H5Sselect_hyperslab(space_id, H5S_SELECT_SET,
			      start, NULL, count, NULL);

	  for (i = 0; i < rank; ++i)
	       start[i] = 0;
	  strides[rank - 1] = data_copy_stride;
	  count[rank - 1] *= data_copy_stride;
	  mem_space_id = H5Screate_simple(rank, count, NULL);
	  count[rank - 1] = local_dims[rank - 1];
	  H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
			      start, data_copy_stride <= 1 ? NULL : strides,
			      count, NULL);
     }
     else { /* this can happen on leftover processes in MPI */
	  H5Sselect_none(space_id);
	  mem_space_id = H5Scopy(space_id); /* can't create an empty space */
	  H5Sselect_none(mem_space_id);
	  do_write = 0;  /* HDF5 complains about empty dataspaces otherwise */
     }

     /*******************************************************************/
     /* Write the data, then free all the stuff we've allocated. */

     if (do_write)
	  H5Dwrite(data_id.id, type_id, mem_space_id, space_id, H5P_DEFAULT, 
		   data_copy);

     if (free_data_copy)
	  free(data_copy);
     H5Sclose(mem_space_id);
     free(count);
     free(strides);
     free(start);
     free(dims);
     H5Sclose(space_id);
#endif
}

#if defined(HAVE_HDF5)
/* check if the given name is a dataset in group_id, and if so set d
   to point to a char** with a copy of name. */
static herr_t find_dataset(hid_t group_id, const char *name, void *d)
{
     char **dname = (char **) d;
     H5G_stat_t info;

     H5Gget_objinfo(group_id, name, 1, &info);
     if (info.type == H5G_DATASET) {
	  CHK_MALLOC(*dname, char, strlen(name) + 1);
	  strcpy(*dname, name);
	  return 1;
     }
     return 0;
}
#endif

/*****************************************************************************/

/* Read real data from the file/group 'id', from the dataset 'name'.

   If name is NULL, reads from the first dataset in 'id'.

   If data is non-NULL, then data must have dimensions given in rank
   and dims (* stride); actually, what is read in is the hyperslab given by the
   local_dim0* parameters.  The dataset is read into 'data' with the
   given 'stride'.  Returns the data pointer.

   If data is NULL, then upon output rank and dims point to the size
   of the array, and a pointer to the (malloc'ed) data is returned.
   On input, *rank should point to the maximum allowed rank (e.g. the
   length of the dims array)!  The local_dim* and stride parameters
   are ignored here.

   Returns NULL if the dataset could not be found in id. */
real *matrixio_read_real_data(matrixio_id id,
			      const char *name,
			      int *rank, int *dims,
			      int local_dim0, int local_dim0_start,
			      int stride,
			      real *data)
{
#if defined(HAVE_HDF5)
     hid_t space_id, type_id, data_id, mem_space_id;
     hsize_t *dims_copy, *maxdims;
     char *dname;
     int i;

     CHECK(*rank > 0, "non-positive rank");

     /*******************************************************************/
     /* Open the data set and check the dimensions: */

     if (name) {
	  CHK_MALLOC(dname, char, strlen(name) + 1);
	  strcpy(dname, name);
     }
     else {
	  if (H5Giterate(id.id, "/", NULL, find_dataset, &dname) < 0)
	       return NULL;
     }
     SUPPRESS_HDF5_ERRORS(data_id = H5Dopen(id.id, dname));
     free(dname);
     if (data_id < 0)
	  return NULL;

     CHECK((space_id = H5Dget_space(data_id)) >= 0,
	   "error in H5Dget_space");

     {
	  int filerank = H5Sget_simple_extent_ndims(space_id);

	  if (data) {
	       CHECK(*rank == filerank,
		     "rank in HDF5 file doesn't match expected rank");
	  }
	  else {
	       CHECK(*rank >= filerank,
		     "rank in HDF5 file is too big");
	       *rank = filerank;
	  }
     }
     
     CHK_MALLOC(dims_copy, hsize_t, *rank);
     CHK_MALLOC(maxdims, hsize_t, *rank);

     H5Sget_simple_extent_dims(space_id, dims_copy, maxdims);
     free(maxdims);

     if (data)
	  for (i = 0; i < *rank; ++i) {
	       CHECK(dims_copy[i] == dims[i],
		     "array size in HDF5 file doesn't match expected size");
	  }
     else
	  for (i = 0; i < *rank; ++i)
	       dims[i] = dims_copy[i];

#if defined(SCALAR_SINGLE_PREC)
     type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
     type_id = H5T_NATIVE_LDOUBLE;
#else
     type_id = H5T_NATIVE_DOUBLE;
#endif

     /*******************************************************************/
     /* Before we can read the data from the data set, we must define
	the dimensions and "selections" of the arrays to be read & written: */

     if (data) {
	  start_t *start;
	  hsize_t *strides, *count;

	  CHK_MALLOC(start, start_t, *rank);
	  CHK_MALLOC(strides, hsize_t, *rank);
	  CHK_MALLOC(count, hsize_t, *rank);
	  
	  for (i = 0; i < *rank; ++i) {
	       start[i] = 0;
	       strides[i] = 1;
	       count[i] = dims[i];
	  }
	  
	  dims_copy[0] = local_dim0;
	  dims_copy[*rank - 1] *= stride;
	  start[0] = 0;
	  strides[*rank - 1] = stride;
	  count[0] = local_dim0;
	  mem_space_id = H5Screate_simple(*rank, dims_copy, NULL);
	  H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
			      start, strides, count, NULL);
	  
	  start[0] = local_dim0_start;
	  count[0] = local_dim0;
	  H5Sselect_hyperslab(space_id, H5S_SELECT_SET,
			      start, NULL, count, NULL);

	  free(count);
	  free(strides);
	  free(start);
     }
     else {
	  int N = 1;
	  for (i = 0; i < *rank; ++i)
	       N *= dims[i];
	  CHK_MALLOC(data, real, N);

	  mem_space_id = H5S_ALL;
	  H5Sclose(space_id);
	  space_id = H5S_ALL;
     }

     /*******************************************************************/
     /* Read the data, then free all the H5 identifiers. */

     CHECK(H5Dread(data_id, type_id, mem_space_id, space_id, H5P_DEFAULT, 
		   data) >= 0,
	   "error reading HDF5 dataset");

     if (mem_space_id != H5S_ALL)
	  H5Sclose(mem_space_id);
     free(dims_copy);
     if (space_id != H5S_ALL)
	  H5Sclose(space_id);
     H5Dclose(data_id);

     return data;
#else
     CHECK(0, "no matrixio implementation is linked");
     return NULL;
#endif
}
