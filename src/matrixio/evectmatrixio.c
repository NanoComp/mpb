/* Copyright (C) 1999-2020 Massachusetts Institute of Technology.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "config.h"

#include <check.h>
#include <matrices.h>

#include "matrixio.h"

void evectmatrixio_writeall_raw(const char *filename, evectmatrix a) {
  int dims[4], start[4] = {0, 0, 0, 0};
  const int rank = 4;
  matrixio_id file_id, data_id;

  dims[0] = a.N;
  dims[1] = a.c;
  dims[2] = a.p;
  dims[3] = SCALAR_NUMVALS;

  start[0] = a.Nstart;

  file_id = matrixio_create(filename);
  data_id = matrixio_create_dataset(file_id, "rawdata", NULL, rank, dims);

  dims[0] = a.localN;
  matrixio_write_real_data(data_id, dims, start, 1, (real *)a.data);

  matrixio_close_dataset(data_id);
  matrixio_close(file_id);
}

void evectmatrixio_readall_raw(const char *filename, evectmatrix a) {
  int rank = 4, dims[4];
  matrixio_id file_id;

  dims[0] = a.N;
  dims[1] = a.c;
  dims[2] = a.p;
  dims[3] = SCALAR_NUMVALS;

  file_id = matrixio_open(filename, 1);

  CHECK(matrixio_read_real_data(file_id, "rawdata", &rank, dims, a.localN, a.Nstart, 1,
                                (real *)a.data),
        "error reading data set in file");

  matrixio_close(file_id);
}
