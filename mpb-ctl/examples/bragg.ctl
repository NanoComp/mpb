; Compute the bands at the X point for a quarter-wave stack Bragg
; mirror (this is the point that defines the band gap edges).

; Using define-param, we can set the x resolution on the command line
; with e.g. 'mpb nx=32 bragg.ctl'; defaults to 16 grid points.
(define-param nx 16)

(set! geometry 
      (list
       (make cylinder 
	 (material (make dielectric (epsilon 9.0)))
	 (center (vector3 0)) (axis (vector3 1)) 
	 (radius infinity) (height 0.25))))

(define-param kx 0.5)
(set! k-points (list (vector3 kx 0 0)))

(set! grid-size (vector3 nx 1 1))
(set! num-bands 8)

(run-tm) ; note that TM and TE bands are degenerate, so we only need TM
