; A simple function to compute the density of states (DOS) by summing
; Gaussian smoothing functions around each computed frequency point.
; This scheme was suggested by Xavier Gonze and implemented in ABINIT
; (http://www.mapr.ucl.ac.be/Fr/PCPM/ABINIT/), with further
; suggestions by Doug Allan of Corning.  (There are also more
; sophisticated algorithms employing the group velocity, but we have
; not implemented them yet.)

; To apply it to output, say, the density of states at 100 points
; points in the frequency range 0 to 1 you would do the following
; (*after* using (run) to compute the frequencies at a bunch of
; k-points filling the irreducible Brillouin zone):
;
; (include "dos.scm")
; (print-dos 0 1 100)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Given a list of frequencies, freqs, and a broadening width df,
; return a function (dos f) that returns the density of states at a
; frequency f.  The integral of dos equals the number of elements of
; freqs.  See also broaden, below, to compute a default df.
(define ((broaden-df freqs df) f)
  (/ (fold-left
      + 0 (map (lambda (f0) (exp (- (sqr (/ (- f f0) df))))) freqs))
     (* 2 df (sqrt (atan 1)))))

; Return the median difference between consecutive numbers in the list nums.
(define (median-diff nums)
  (let ((snums (sort nums <)) (n (- (length nums) 1)))
    (let ((sdiff (sort (map (lambda (x y) (- y x))
			    (reverse (cdr (reverse snums)))
			    (cdr snums)) <)))
      (* 0.5 (+ (list-ref sdiff (quotient n 2))
		(list-ref sdiff (- (quotient (+ n 1) 2) 1)))))))

; Like broaden-df, above, but use a default df given by (median-diff freqs).
(define (broaden freqs)
  (broaden-df freqs (median-diff freqs)))

; Return the broaden (dos) function applied to all-freqs, the list of
; all frequencies computed from the last run.  (Actually, all-freqs
; is the list of lists of frequencies at each k-point.)
(define (all-freqs-broaden)
  (broaden (fold-left append '() all-freqs)))

; Output the DOS from all-freqs, at num-freq points from freq-min to freq-max.
(define (print-dos freq-min freq-max num-freq)
  (let ((dos (all-freqs-broaden)))
    (map (lambda (f) (print "dos:, " f ", " (dos f) "\n"))
	 (interpolate (max 0 (- num-freq 2)) (list freq-min freq-max)))))
