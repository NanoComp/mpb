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

#ifndef _MPB_VERBOSITY_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
A global "verbosity" level.
   0: minimal output
   1: a little
   2: a lot (default)
   3: debugging
*/
extern int mpb_verbosity;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _MPB_VERBOSITY_H */
