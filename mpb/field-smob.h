/* Copyright (C) 1999-2012, Massachusetts Institute of Technology.
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

#ifndef FIELD_SMOB_H
#define FIELD_SMOB_H

#include "ctl-io.h"
#include "my-smob.h"

#include <matrices.h>

extern long scm_tc16_smob_field_smob;

typedef enum { 
     RSCALAR_FIELD_SMOB, CVECTOR_FIELD_SMOB, CSCALAR_FIELD_SMOB
} field_smob_type;

typedef struct {
     field_smob_type type;
     char type_char;
     union {
	  real *rs;
	  scalar_complex *cv;
	  scalar_complex *cs;
     } f;
     int nx, ny, nz, N;
     int local_ny, local_y_start;
     int last_dim, last_dim_size, other_dims;
} field_smob;

#define FIELD_P(X) T_SMOB_P(field_smob, X)
#define FIELD(X) T_SMOB(field_smob, X)
#define SAFE_FIELD(X) (cur_fieldp(X) ? update_curfield_smob() : SAFE_T_SMOB(field_smob, X))

#define RSCALAR_FIELD_P(X) (FIELD_P(X) && ((FIELD(X))->type == RSCALAR_FIELD_SMOB))
#define CSCALAR_FIELD_P(X) (FIELD_P(X) && ((FIELD(X))->type == CSCALAR_FIELD_SMOB))
#define CVECTOR_FIELD_P(X) (FIELD_P(X) && ((FIELD(X))->type == CVECTOR_FIELD_SMOB))

extern field_smob *update_curfield_smob(void);
extern void register_field_smobs(void);
extern field_smob *assert_field_smob(SCM fo);

#endif /* FIELD_SMOB_H */
