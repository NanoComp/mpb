; Compute the bands at the X point for a quarter-wave stack Bragg
; mirror (this is the point that defines the band gap edges).

; Using define-param, we can set the x resolution on the command line
; with e.g. 'mpb nx=32 bragg.ctl'; defaults to 16 grid points.
(define-param nx 16)

; the high and low indices:
(define-param n-lo 1.0)
(define-param n-hi 3.0)

(define-param w-hi (/ n-lo (+ n-hi n-lo))) ; a quarter-wave stack

(set! default-material (make dielectric (index n-lo)))
(set! geometry 
      (list
       (make cylinder 
	 (material (make dielectric (index n-hi)))
	 (center 0 0 0) (axis 1 0 0) 
	 (radius infinity) (height w-hi))))

(define-param kx 0.5)
(set! k-points (list (vector3 kx 0 0)))

(set! grid-size (vector3 nx 1 1))
(set! num-bands 8)

(run-tm) ; note that TM and TE bands are degenerate, so we only need TM
