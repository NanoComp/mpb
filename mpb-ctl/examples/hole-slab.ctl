; Photonic crystal slab consisting of a triangular lattice of air
; holes in a finite-thickness dielectric slab.  See the paper:
; S. G. Johnson, S. Fan, P. R. Villeneuve, J. D. Joannopoulos,
; L. A. Kolodziejski, "Guided modes in photonic crystal slabs,"
; PRB 60, 5751 (August 1999).

(define-param h 0.5) ; the thickness of the slab
(define-param eps 12.0) ; the dielectric constant of the slab
(define-param r 0.3) ; the radius of the holes
(define-param supercell-h 4) ; height of the supercell

; triangular lattice with vertical supercell:
(set! geometry-lattice (make lattice (size 1 1 supercell-h)
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))

(set! geometry
      (list (make block (material (make dielectric (epsilon eps)))
		  (center 0) (size 1 1 h))
	    (make cylinder (material air)
		  (center 0) (radius r) (height supercell-h))))

; 1st Brillouin zone of a triangular lattice:
(define Gamma (vector3 0 0 0))
(define M (vector3 0 0.5 0))
(define K (vector3 (/ -3) (/ 3) 0))

(define-param only-K false) ; run with only-K=true to only do this k-point
(if only-K
    (set! k-points (list K))
    (set! k-points (interpolate 4 (list Gamma M K Gamma))))

(set! grid-size (vector3 32 32 32))
(set! num-bands 20)

; Define a band function to output the magnetic fields only at the K point:
(define (output-hfield-z-only-at-K band)
  (if (vector3= current-k K)
      (output-hfield-z band)))

(run output-hfield-z-only-at-K)

(display-eigensolver-stats)
