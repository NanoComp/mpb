; A line-defect waveguide in a 2d triangular lattice of dielectric
; rods (c.f. tri-rods.ctl), formed by a row of missing rods along the
; "x" direction.  (Here, "x" and "y" refer to the first and second
; basis directions.)  This structure supports a single guided band
; within the band gap, much like the analogous waveguide in a square
; lattice of rods (see "Photonic Crystals" by Joannopoulos et al.).

(define-param supercell-y 7) ; the (odd) number of lateral supercell periods

(set! geometry-lattice (make lattice 
			 (basis1 (/ (sqrt 3) 2) 0.5)
			 (basis2 (/ (sqrt 3) 2) -0.5)
			 (size 1 supercell-y 1)))

(define-param eps 12) ; the dielectric constant of the rods
(define-param r 0.2) ; the rod radius in the bulk crystal

(set! geometry (list (make cylinder
		       (center 0 0 0) (radius r) (height infinity)
		       (material (make dielectric (epsilon eps))))))

(set! geometry 
      (append
       ; duplicate the bulk crystal rods over the supercell:
       (geometric-objects-lattice-duplicates geometry)

       ; add a rod of air, to erase a row of rods and form a waveguide:
       (list 
	(make cylinder (center 0) (radius r) (height infinity)
	      (material air)))))

(define Gamma (vector3 0 0 0))
(define K' (lattice->reciprocal (vector3 0.5 0 0))) ; edge of Brillouin zone.
(set! k-points (interpolate 4 (list Gamma K')))

; the bigger the supercell, the more bands you need to compute to get
; to the defect modes (the lowest band is "folded" supercell-y times):
(define-param extra-bands 5) ; number of extra bands to compute above the gap
(set! num-bands (+ supercell-y extra-bands))

(define-param res 32) ; the resolution (grid points/a)
(set! grid-size (vector3 res (* res supercell-y) 1))

; Compute the TM modes, outputting the Ez field in the *middle* of the
; band.  (In general, the guided mode in such an air defect may have
; exited the gap by the time it reaches the edge of the Brillouin
; zone at K'.)
(run-tm 
 (output-at-kpoint (list-ref k-points (quotient (length k-points) 2))
		   fix-efield-phase output-efield-z))
