; Copyright (C) 1999 Massachusetts Institute of Technology.
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
;  This program is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with this program; if not, write to the Free Software
;  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

(define-class material-type no-parent
  (define-property epsilon 'number no-default positive?)
  (define-property conductivity 'number (make-default 0.0)))

; use the solid geometry classes, variables, etcetera in libgeom:
; (one specifications file can include another specifications file)
(include "../../libctl/utils/libgeom/geom.scm")

; ****************************************************************

; Add some predefined variables, for convenience:

(define vacuum (make material-type (epsilon 1.0)))
(define air vacuum)

(define infinity 1.0e20) ; big number for infinite dimensions of objects

(set! default-material air)

; ****************************************************************

(define-input-var k-points '() (make-list-type 'vector3))

(define-input-var num-bands 1 'integer)
(define-input-var tolerance 1.0e-4 'number positive?)
(define-input-var target-freq 0.0 'number positive?)

(define-input-var grid-size (vector3 16 16 16) 'vector3
  (lambda (v) (vector-for-all? v (lambda (x) (and (> x 0) (integer? x))))))
(define-input-var mesh-size 7 'integer positive?)

(define-output-var freqs (make-list-type 'number))

; ****************************************************************

; (init-params) initializes the geometry, etcetera, and does everything
; else that's needed to get ready for an eigenvalue calculation.
; This should be called after the input variables are changed.
(define-external-function init-params true false no-return-value)

; (solve-kpoint kpoint) solves for the specified bands at the given k point.
; Requires that (init-params) has been called, and does not re-read the
; input variables, but does write the output vars.
(define-external-function solve-kpoint false true no-return-value 'vector3)

; define a (run) function to do a vanilla calculation (solve for the
; bands at all the k points).
(define (run)
  (set! interactive false)  ; don't be interactive if we call (run)
  (init-params)
  (map solve-kpoint k-points))

; it would be nice to have this routine later:
; (define-external-function energy-in-object false false
;   'number 'geometric-object)
