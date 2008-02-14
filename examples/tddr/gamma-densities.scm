
(load "utl.scm")
(load "gsl.scm")

;; Read the mean and variance from standard input.
(let ((mean (read))
      (variance (read)))
  ;; For now, I'll just hardwire the range of the plot in,
  ;; as well as the number of interpolating points.
  (let ((range-lo 10.0)
	(range-hi 250)
	(bin-size 5)			; That which is used to make
					; the data histograms; has nothing
					; to do with this algorithm, just
					; rescales the output.
	(x-count 500))
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
    
