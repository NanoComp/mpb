; Compute the bands at the X point for a quarter-wave stack Bragg
; mirror (this is the point that defines the band gap edges).

; the high and low indices:
(define-param n-lo 1.0)
(define-param n-hi 3.0)

(define-param w-hi (/ n-lo (+ n-hi n-lo))) ; a quarter-wave stack

(set! geometry-lattice (make lattice (size 1 no-size no-size))) ; 1d cell
(set! default-material (make dielectric (index n-lo)))
(set! geometry 
      (list
       (make cylinder 
	 (material (make dielectric (index n-hi)))
	 (center 0 0 0) (axis 1 0 0) 
	 (radius infinity) (height w-hi))))

(define-param kx 0.5)
(set! k-points (list (vector3 kx 0 0)))

(set-param! resolution 32)
(set-param! num-bands 8)

(run-tm output-hfield-y) ; note that TM and TE bands are degenerate, so we only need TM
