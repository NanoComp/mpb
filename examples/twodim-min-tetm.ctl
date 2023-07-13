(define-param tri? false)
(set-param! mesh-size 1)
(define-param mg-res 32)
(set-param! resolution 64)
(define-param eps 12)
(set-param! tolerance 1e-13)

(if (not deterministic?)
    (let ((time (gettimeofday)))
      (set! *random-state*
	    (seed->random-state (+ (car time)
				   (cdr time))))))
(define-param r 0.1)
(define-param rand-init? true)
(define-param kind U-SUM)
(define mg (make material-grid
	     (epsilon-min 1)
	     (epsilon-max eps)
	     (material-grid-kind kind)
	     (size mg-res mg-res 1)
	     (matgrid-init
	      (if rand-init?
		  (lambda (x y z) (random:uniform))
		  (lambda (x y z) 
		    (if (< (+ (* x x) (* y y)) (* r r))
			0.0 1.0))))))

(if tri?
    (set! geometry-lattice (make lattice (size 1 1 no-size)
				 (basis1 (/ (sqrt 3) 2) 0.5)
				 (basis2 (/ (sqrt 3) 2) -0.5)))
    (set! geometry-lattice (make lattice (size 1 1 no-size))))

(set! geometry 
      (list 
       (make block (e1 1 0) (e2 0 1) 
	     (center 0) (size 1 1) ( material mg))))

(define-param band 1)
(set! num-bands 1)
(define-param gridfile false)
(if gridfile
    (load-material-grid! mg gridfile (vector3 supercell supercell 1)))
(define-param noise 0)
(if (> noise 0)
    (randomize-material-grid! mg noise))

(init-params NO-PARITY true)
(set! interactive? false)

(define-param epsfile false)
(if epsfile
    (material-grids-match-epsilon-file! epsfile 1e-4))

(define-param gaptol 0)
(define-param epstol 0)
(define-param maxeval 40)

(define-param kx 0.1234)
(define-param ky 0.3456)

(material-grids-min-tetm-gap
 (vector3 kx ky) band
 gaptol epstol maxeval 0)

(save-material-grid mg (string-append (get-filename-prefix) "grid"))

(set! k-points (list (vector3 kx ky)))
(set! num-bands 2)
(run)
