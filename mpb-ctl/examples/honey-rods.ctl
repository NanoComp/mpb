; A honeycomb lattice of dielectric rods in air.  (This structure has
; a complete (overlapping TE/TM) band gap.)  A honeycomb lattice is really
; just a triangular lattice with two rods per unit cell, so we just
; take the lattice, k-points, etcetera from tri-rods.ctl.

(set! num-bands 8)

(define-param r 0.15) ; the rod radius
(define-param eps 12) ; the rod dielectric constant

; triangular lattice:
(set! geometry-lattice (make lattice
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))

; Two rods per unit cell, at the correct positions to form a honeycomb
; lattice, and arranged to have inversion symmetry:
(set! geometry
      (list (make cylinder
	      (center (/ 6) (/ 6) 0) (radius r) (height infinity)
	      (material (make dielectric (epsilon eps))))
	    (make cylinder
	      (center (/ -6) (/ -6) 0) (radius r) (height infinity)
	      (material (make dielectric (epsilon eps))))))

; The k-points list, for the Brillouin zone of a triangular lattice:
(set! k-points (list (vector3 0 0 0)          ; Gamma
                     (vector3 0 0.5 0)        ; M
                     (vector3 (/ -3) (/ 3) 0) ; K
                     (vector3 0 0 0)))        ; Gamma
(define-param k-interp 4) ; number of k-points to interpolate
(set! k-points (interpolate k-interp k-points))

(define-param res 64) ; the grid resolution
(set! grid-size (vector3 res res 1))

(run-tm)
(run-te)

; Since there is a complete gap, we could instead see it just by using:
; (run)
; The gap is between bands 12 and 13 in this case.  (Note that there is
; a false gap between bands 2 and 3, which disappears as you increase the
; k-point resolution.)
