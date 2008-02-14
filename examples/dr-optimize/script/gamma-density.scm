
(load "utl.scm")
(load "gsl.scm")

;; Read the mean and variance from standard input.
(let ((mean (read))
      (variance (read)))
  ;; For now, I'll just hardwire the range of the plot in,
  ;; as well as the number of interpolating points.
  (let ((range-lo 10.0)
	(range-hi 1600)
	(bin-size 5)			; That which is used to make
					; the data histograms; has nothing
					; to do with this algorithm, just
					; rescales the output.
	(x-count 800))
    ;; One bad thing that happens is that all the output populations come out
    ;; the same, giving variance 0.  In this case, we artificially arrange for
    ;; the standard deviation to be half the bin-size!  Note that this is the
    ;; bin-size used for histogramming, and the idea is to make the gamma
    ;; distribution look like the histogram distribution.
    (or (< 0.0
	   variance)
	(set! variance (/ (* bin-size
			     bin-size)
			  4.0)))
    ;; Another thing that can happen is that the mean is 0.  In this case we
    ;; artificially change the mean to the center of the lowest histogram bin.
    (or (< 0.0
	   mean)
	(set! mean (/ bin-size
		      2.0)))
    (let ((x-range-width (- range-hi
			    range-lo)))
      (let ((x-delta (/ x-range-width
			x-count)))
	(do ((x-index 0
		      (+ x-index
			 1))
	     (x-value range-lo
		      (+ x-value
			 x-delta)))
	    ((= x-index
		x-count)
	     #f)
	  (display x-value)
	  (display "\t")
	  (display (* bin-size
		      (gamma-density x-value mean variance)))
	  (newline))))))
    
