; This file contains the Scheme commands from the user tutorial section
; of the manual.  It is meant to be run interactively by the user.

; *** Our First Band Structure ***

(display "********** Square lattice of rods in air **********\n")

(set! num-bands 8)
(set! k-points (list (vector3 0 0 0)     ; Gamma
                     (vector3 0.5 0 0)   ; X
                     (vector3 0.5 0.5 0) ; M
                     (vector3 0 0 0)))   ; Gamma
(set! k-points (interpolate 4 k-points))
(set! geometry (list (make cylinder 
                       (center 0 0 0) (radius 0.2) (height infinity)
                       (material (make dielectric (epsilon 12))))))
(set! grid-size (vector3 32 32 1))

(display "********** Square lattice of rods: TE bands**********\n")
(run-te)

(display "********** Square lattice of rods: TM bands **********\n")
(run-tm)

(display "********** Square lattice of rods: TM, w/efield **********\n")
(run-tm output-efield-z)

(display "********** Square lattice of rods: TE, w/hfield & dpwr **********\n")
(run-te (output-at-kpoint (vector3 0.5 0 0) output-hfield-z output-dpwr))

; *** Bands of a Triangular Lattice ***

(display "********** Triangular lattice of rods in air **********\n")

(set! geometry-lattice (make lattice
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))

(set! k-points (list (vector3 0 0 0)          ; Gamma
                     (vector3 0 0.5 0)        ; M
                     (vector3 (/ -3) (/ 3) 0) ; K
                     (vector3 0 0 0)))        ; Gamma
(set! k-points (interpolate 4 k-points))

(run-tm)

; *** Maximizing the First TM Gap ***

(display "********** Maximizing the first TM gap **********\n")

(define (first-tm-gap r)
  (set! geometry (list (make cylinder
                         (center 0 0 0) (radius r) (height infinity)
                         (material (make dielectric (epsilon 12))))))
  (run-tm)
  (if (null? gap-list)
      0.0  ; no gaps, return 0
      (caar gap-list))) ; the first element of the first list in gap-list

(set! num-bands 2)
(set! mesh-size 7) ; increase from default value of 3

(define result (maximize first-tm-gap 0.1 0.1 0.5))
(display-many "radius at maximum: " (max-arg result) "\n")
(display-many "gap size at maximum: " (max-val result) "\n")

(set! mesh-size 3) ; reset to default value of 3

; *** Finding a Point-defect State ***

(display "********** 5x5 point defect **********\n")

(set! geometry-lattice (make lattice (size 5 5 1)))
(set! geometry (geometric-objects-lattice-duplicates geometry))
(set! geometry (append geometry 
                       (list (make cylinder (center 0 0 0) 
                                   (radius 0.2) (height infinity)
                                   (material air)))))
(set! grid-size (vector3 (* 16 5) (* 16 5) 1))
(set! k-points (list (vector3 0.5 0.5 0)))

(set! num-bands 50)
(run-tm)

(output-efield-z 25)
(get-dfield 25)  ; compute the D field for band 25
(compute-field-energy)  ; compute the energy density from D
(display-many
 "energy in cylinder: "
 (compute-energy-in-objects (make cylinder (center 0 0 0)
                                  (radius 1.0) (height infinity)
                                  (material air)))
 "\n")

(display "********** 5x5 point defect, targeted solver **********\n")
(set! num-bands 1)  ; only need to compute a single band, now!
(set! target-freq (/ (+ 0.2812 0.4174) 2))
(set! tolerance 1e-8)
(run-tm)

; *** Tuning the Point-defect Mode ***

(display "********** Tuning the 5x5 point defect **********\n")

(define old-geometry geometry) ; save the 5x5 grid with a missing rod
(define (rootfun eps)
  ; add the cylinder of epsilon = eps to the old geometry:
  (set! geometry (append old-geometry
                         (list (make cylinder (center 0 0 0)
                                     (radius 0.2) (height infinity)
                                     (material (make dielectric
                                                 (epsilon eps)))))))
  (run-tm)  ; solve for the mode (using the targeted solver)
  (display-many "epsilon = " eps " gives freq. = " (list-ref freqs 0) "\n")
  (- (list-ref freqs 0) 0.314159))  ; return 1st band freq. - 0.314159

(define rooteps (find-root rootfun 0.01 1 12))
(display-many "root (value of epsilon) is at: " rooteps "\n")

(define rootval (rootfun rooteps))
(display-many "root function at " rooteps " = " rootval "\n")
