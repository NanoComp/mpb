; A line defect in a 2d triangular lattice of air holes, formed by a
; row of missing holes along the "x" direction.  (Here, "x" and "y"
; refer to the first and second basis directions.)

(define-param supercell-y 7) ; the (odd) number of lateral supercell periods

(set! geometry-lattice (make lattice 
			 (basis1 (/ (sqrt 3) 2) 0.5)
			 (basis2 (/ (sqrt 3) 2) -0.5)
			 (size 1 supercell-y 1)))

(set! default-material (make dielectric (epsilon 12)))

(set! geometry (list (make cylinder
		       (center 0 0 0) (radius 0.3) (height infinity)
		       (material air))))

(define-param rd 0.2) ; the defect rod radius
(set! geometry 
      (append
       (geometric-objects-lattice-duplicates geometry)
       (list (make cylinder (center 0) (radius rd) (height infinity)
		   (material default-material)))))

(define Gamma (vector3 0 0 0))
(define K' (lattice->reciprocal (vector3 0.5 0 0))) ; edge of Brillouin zone.
(set! k-points (interpolate 4 (list Gamma K')))

; the bigger the supercell, the more bands you need to compute to get
; to the defect modes.  
(define-param extra-bands 5) ; number of extra bands to compute above the gap
(set! num-bands (+ supercell-y extra-bands))

(define-param res 32) ; the resolution (grid points/a)
(set! grid-size (vector3 res (* res supercell-y) 1))

; Compute the TE modes, outputting the Hz field at the K' point:
(run-te (output-at-kpoint K' fix-hfield-phase output-hfield-z))
