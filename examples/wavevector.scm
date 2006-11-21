; Code to output the data needed for a wavevector diagram at a grid of
; k points, suitable for plotting with a contour-plot program.

(define (kgrid kx-min kx-max ky-min ky-max nkx nky)
  (map (lambda (kx)
	 (interpolate nky (list (vector3 kx ky-min) (vector3 kx ky-max))))
       (interpolate nkx (list kx-min kx-max))))

; output a bunch of lines of the form:
; kgrid:, band#, kx, frequencies at kys...
; Frequencies above the light cone omega > c |k| / n-lightcone are
; excluded (multiplied by -1); set n-lightcone = 0 to disable this.
(define (wavevector-diagram kgrid parity n-lightcone)
  (map (lambda (kylist)
	 (set! k-points kylist)
	 (run-parity parity true)
	 (map
	  (lambda (band)
	    (print "kgrid:, " band ", " (vector3-x (car kylist)))
	    (map (lambda (freqs k) 
		   (print ", " 
			  (* (if (and (positive? n-lightcone)
				      (> (list-ref freqs (- band 1))
					 (* n-lightcone
					    (vector3-norm 
					     (reciprocal->cartesian k)))))
				 -1
				 1)
			     (list-ref freqs (- band 1)))))
		 all-freqs k-points)
	    (print "\n"))
	  (arith-sequence 1 1 num-bands)))
	 kgrid))

  