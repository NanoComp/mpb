; Compute band structure for a square lattice of dielectric rods
; in air.

; Define various parameters with define-param so that they are
; settable from the command-line (with mpb <param>=<value>):
(define-param r 0.2) ; radius of the rods
(define-param eps 11.56) ; dielectric constant
(define-param k-interp 4) ; number of k points to interpolate

(define Si (make dielectric (epsilon eps)))

(set! geometry-lattice (make lattice (size 1 1 no-size))) ; 2d cell

(set! geometry 
      (list
       (make cylinder 
	 (material Si) 
	 (center 0 0) (radius r) (height infinity))))
      
(define Gamma (vector3 0 0 0))
(define X (vector3 0.5 0 0))
(define M (vector3 0.5 0.5 0))
(set! k-points (interpolate k-interp (list Gamma X M Gamma)))

(set-param! resolution 32)
(set-param! num-bands 8)

; Compute the TE and TM bands.  Wrap in the (begin-time message ...)
; construct from libctl so that we report the total elapsed time:
(begin-time
 "total time for both TE and TM bands: "
 (run-te)
 (run-tm))

(display-eigensolver-stats)
