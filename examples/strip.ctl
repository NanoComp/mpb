; Compute modes of a rectangular Si strip waveguide on top of oxide. 
; Note that you should only pay attention, here, to the guided modes,
; which are the modes whose frequency falls under the light line --
; that is, frequency < beta / 1.45, where 1.45 is the SiO2 index.

; Since there's no special lengthscale here, I'll just
; use microns.  In general, if you use units of x, the frequencies
; output are equivalent to x/lambda; so here, the freqeuncies will be
; output as um/lambda, e.g. 1.5um would correspond to the frequency
; 1/1.5 = 0.6667.

; define parameters which can easily be changed on the command-line
; (e.g.: mpb w=0.5 buried-strip.ctl)
(define-param w 0.3) ; Si width (um)
(define-param h 0.25) ; Si height (um)

(define Si (make dielectric (index 3.45)))
(define SiO2 (make dielectric (index 1.45)))

; Define the computational cell.  We'll make x the propagation direction.
; the other cell sizes should be big enough so that the boundaries are
; far away from the mode field.
(define-param sc-y 2) ; supercell width (um)
(define-param sc-z 2) ; supercell height (um)
(set! geometry-lattice (make lattice (size no-size sc-y sc-z)))

; define the 2d blocks for the strip and substrate
(set! geometry
      (list
       (make block (size infinity infinity (* 0.5 (- sc-z h)))
	     (center 0 0 (* 0.25 (+ sc-z h))) (material SiO2))
       (make block (size infinity w h) (center 0 0 0) (material Si))))


; The k (i.e. beta, i.e. wavenumber) points to look at, in units of 2*pi/um.
; We'll look at num-k points from k-min to k-max.
(define-param num-k 9)
(define-param k-min 0.1)
(define-param k-max 3.0)
(set! k-points (interpolate num-k (list (vector3 k-min) (vector3 k-max))))

(set-param! resolution 32) ; pixels/um

; Increase this to see more modes.  (The guided ones are the ones below the
; light line, i.e. those with frequencies < kmag / 1.45, where kmag
; is the corresponding column in the output if you grep for "freqs:".)
(set-param! num-bands 4)

(set-param! filename-prefix "strip-") ; use this prefix for output files

; compute num-bands lowest frequencies as a function of k. Also display
; "parities", i.e. whether the mode is symmetric or anti-symmetric
; through the y=0 and z=0 planes.
(run display-yparities display-zparities)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Above, we outputted the dispersion relation: frequency (omega) as a
; function of wavevector kx (beta).  Alternatively, you can compute
; beta for a given omega -- for example, you might want to find the
; modes and wavevectors at a fixed wavelength of 1.55 microns.  You
; can do that using the find-k function:

(define-param omega (/ 1 1.55)) ; frequency corresponding to 1.55um

; Output the x component of the Poynting vector for num-bands bands at omega
(find-k NO-PARITY omega 1 num-bands (vector3 1) 1e-3
	(* omega 3.45) (* omega 0.1) (* omega 4) output-poynting-x
	display-yparities display-group-velocities)

