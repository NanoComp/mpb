; A triangular lattice of dielectric rods in air.  (This structure has
; a band-gap for TM fields.)  This file is used in the "Data Analysis
; Tutorial" section of the MPB manual.

(set! num-bands 8)

(set! geometry-lattice (make lattice
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))
(set! geometry (list (make cylinder
                       (center 0 0 0) (radius 0.2) (height infinity)
                       (material (make dielectric (epsilon 12))))))

(set! k-points (list (vector3 0 0 0)          ; Gamma
                     (vector3 0 0.5 0)        ; M
                     (vector3 (/ -3) (/ 3) 0) ; K
                     (vector3 0 0 0)))        ; Gamma
(set! k-points (interpolate 4 k-points))

(set! grid-size (vector3 32 32 1))

(run-tm (output-at-kpoint (vector3 (/ -3) (/ 3) 0) output-efield-z))
(run-te)

