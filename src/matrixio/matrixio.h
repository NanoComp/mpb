/* Copyright (C) 1999, 2000 Massachusetts Institute of Technology.
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

#ifndef MATRIXIO_H
#define MATRIXIO_H

#include <matrices.h>

#if defined(HAVE_HDF5)
#  include <hdf5.h>
typedef hid_t matrixio_id;
#else /* no HDF */
typedef int matrixio_id; /* dummy */
#endif

extern matrixio_id matrixio_create(const char *fname);
extern matrixio_id matrixio_open(const char *fname);
extern void matrixio_close(matrixio_id id);

extern matrixio_id matrixio_create_sub(matrixio_id id,
                                const char *name, const char *description);
extern void matrixio_close_sub(matrixio_id id);

extern matrixio_id matrixio_create_dataset(matrixio_id id,
                                    const char *name, const char *description,
                                    int rank, const int *dims);
extern void matrixio_close_dataset(matrixio_id data_id);

extern void matrixio_write_real_data(matrixio_id data_id,
                              const int *local_dims, const int *local_start,
                              int stride,
                              real *data);
extern real *matrixio_read_real_data(matrixio_id id,
				     const char *name,
				     int *rank, int *dims,
				     int local_dim0, int local_dim0_start,
				     int stride,
				     real *data);

extern void matrixio_write_string_attr(matrixio_id id, const char *name,
				       const char *val);
extern void matrixio_write_data_attr(matrixio_id id, const char *name,
				     const real *val, int rank, 
				     const int *dims);
extern char *matrixio_read_string_attr(matrixio_id id, const char *name);
extern real *matrixio_read_data_attr(matrixio_id id, const char *name,
				     int *rank, int max_rank, int *dims);

extern void evectmatrixio_writeall_raw(const char *filename, evectmatrix a);
extern void evectmatrixio_readall_raw(const char *filename, evectmatrix a);

extern void fieldio_write_complex_field(scalar_complex *field,
					int rank,
					const int dims[3],
					int which_component,
					int local_nx, int local_x_start,
					const int copies[3],
					const real kvector[3],
					real R[3][3],
					matrixio_id file_id);
void fieldio_write_real_vals(real *vals,
                             int rank,
                             const int dims[3],
                             int local_nx, int local_x_start,
                             const int copies[3],
			     matrixio_id file_id);

#endif /* MATRIXIO_H */
