; 2d system: triangular lattice of air holes in dielectric
; This structure has a complete band gap (i.e. a gap in both TE and TM
; simultaneously) for a hole radius of 0.45a and a dielectric constant of
; 12.   (See, e.g., the book "Photonic Crystals" by Joannopoulos et al.)

; first, define the lattice vectors and k-points for a triangular lattice:

(set! geometry-lattice (make lattice (size 1 1 no-size)
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))

(define-param kz 0) ; use non-zero kz to consider vertical propagation

(set! k-points (list (vector3 0 0 kz)          ; Gamma
                     (vector3 0 0.5 kz)        ; M
                     (vector3 (/ -3) (/ 3) kz) ; K
                     (vector3 0 0 kz)))        ; Gamma
(set! k-points (interpolate 4 k-points))

; Now, define the geometry, etcetera:

(define-param eps 12) ; the dielectric constant of the background
(define-param r 0.45) ; the hole radius

(set! default-material (make dielectric (epsilon eps)))
(set! geometry (list (make cylinder (center 0) (material air)
			   (radius r) (height infinity))))

(set-param! resolution 32)
(set-param! num-bands 8)

(if (= kz 0)
    (begin
      (run-te)
      (run-tm))
    (run)) ; if kz != 0 there are no purely te and tm bands
