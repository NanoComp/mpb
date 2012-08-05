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

#ifndef MATRIX_SMOB_H
#define MATRIX_SMOB_H

#include "ctl-io.h"
#include "my-smob.h"

#include <matrices.h>

extern long scm_tc16_smob_evectmatrix;
extern long scm_tc16_smob_sqmatrix;

#define EVECTMATRIX_P(X) T_SMOB_P(evectmatrix, X)
#define EVECTMATRIX(X) T_SMOB(evectmatrix, X)
#define SAFE_EVECTMATRIX(X) SAFE_T_SMOB(evectmatrix, X)

#define SQMATRIX_P(X) T_SMOB_P(sqmatrix, X)
#define SQMATRIX(X) T_SMOB(sqmatrix, X)
#define SAFE_SQMATRIX(X) SAFE_T_SMOB(sqmatrix, X)

extern SCM evectmatrix2scm(evectmatrix m);
extern SCM evectmatrix_sub2scm(evectmatrix m, int p_start, int p_end);
extern SCM sqmatrix2scm(sqmatrix m);

extern integer sqmatrix_size(object mo);
extern cnumber sqmatrix_ref(object mo, integer i, integer j);

extern void register_matrix_smobs(void);
extern sqmatrix *assert_sqmatrix_smob(SCM mo);
extern evectmatrix *assert_evectmatrix_smob(SCM mo);

#endif /* MATRIX_SMOB_H */
