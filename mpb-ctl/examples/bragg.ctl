(define-param nx 16)
(define-param target-frequency 0.0)

(set! geometry 
      (list
       (make cylinder 
	 (material (make material-type (epsilon 9.0)))
	 (center (vector3 0)) (axis (vector3 1)) 
	 (radius infinity) (height 0.25))))

(define X (vector3 0.5 0 0))
(set! k-points (list (vector3 0.5)))

(set! dimensions 2)
(set! grid-size (vector3 nx 1 1))

(set! num-bands 8)
(set! target-freq target-frequency)

(run)
