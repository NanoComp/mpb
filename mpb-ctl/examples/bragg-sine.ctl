; Compute the band structure for a Bragg mirror consisting of a
; sinusoidally-varying dielectric index.

; The index will vary sinusoidally between index-min and index-max:
(define-param index-min 1)
(define-param index-max 3)

(define pi (* 4 (atan 1))) ; 3.14159...

; Define a function of position p (in the lattice basis) that returns
; the material at that position.  In this case, we use the function:
;        index-min + 0.5 * (index-max - index-min)
;                        * (1 + cos(2*pi*x))
; This is periodic, and also has inversion symmetry.
(define (eps-func p)
  (make dielectric
    (index (+ index-min (* 0.5 (- index-max index-min)
			   (+ 1 (cos (* 2 pi (vector3-x p)))))))))

; We'll just make it the default material, so that it goes everywhere.
(set! default-material (make material-function (material-func eps-func)))

(set! k-points (interpolate 9 (list (vector3 0 0 0) (vector3 0.5 0 0))))

(define-param nx 32)
(set! grid-size (vector3 nx 1 1))

(set! num-bands 8)

; the TM and TE bands are degenerate, so we only need TM:
(run-tm)
