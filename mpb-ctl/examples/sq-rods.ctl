(define-param norun false)

(define Si (make material-type (epsilon 12.0)))

(set! geometry 
      (list
       (make cylinder 
	 (material Si) 
	 (center (vector3 0))
	 (radius 0.2) (height infinity))))
      
(define Gamma (vector3 0 0 0))
(define X (vector3 0.5 0 0))
(define M (vector3 0.5 0.5 0))

(set! k-points (interpolate 4 (list Gamma X M Gamma)))

(set! dimensions 2)  ; causes (run) to solve for TE/TM bands.
(set! grid-size (vector3 32 32 1))

(set! num-bands 5)

(if (not norun)
    (run))
