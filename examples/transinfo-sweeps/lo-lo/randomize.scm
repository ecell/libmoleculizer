
(load "utl.scm")
(load "gsl.scm")

(define (gamma-randomize-real positive-real
			      randomize-fraction)
  (let ((desired-sd (* positive-real
		       randomize-fraction)))
    (let ((desired-variance (* desired-sd
			       desired-sd)))
      (sample-gamma positive-real
		    desired-variance))))

;;; Produces integers whose mean is roughly positive-real.
(define (gamma-randomize-integer positive-real
				 randomize-fraction)
  (floor (+ (gamma-randomize-real positive-real
				  randomize-fraction)
	    0.5)))

;;; Copies header (which gives a sort of "variable name" to each column)
;;; from the input-values-file to the output-values-file.  Then randomizes
;;; the values on the first following row in input-values-file to generate
;;; all the remaining rows of the output-values-file
(define (randomize-values input-values-file
			  output-values-file
			  randomize-function
			  randomize-count
			  randomize-fraction)
  ;; Function to randomize one row.
  ;; 
  ;; Since the first column of the values file now gives the output
  ;; directories, we use it to generate names for output
  ;; directories for the randomized version.
  (define (randomize-values-vector values-vector
				   randomize-index)
    (let ((result-vector (make-vector)))
      ;; Generate a directory name by adding the randomize-index as a
      ;; suffix.  Note here we're assuming that the output directory reads
      ;; as a symbol.  Would it be better to put quoted strings in the
      ;; values file?
      (vector-push! result-vector
		    (let ((dir-name-port (open-output-string)))
		      (display (vector-ref values-vector
					   0)
			       dir-name-port)
		      (display "_"
			       dir-name-port)
		      (display randomize-index
			       dir-name-port)
		      (close-port dir-name-port)))
      ;; Randomize each of the (remaining) values on this line of the values
      ;; file into the corresponding position in the current output line.
      (do ((value-index 1
			(+ value-index
			   1)))
	  ((<= (vector-length values-vector)
	       value-index)
	   result-vector)
	(vector-push! result-vector
		      (randomize-function (vector-ref values-vector
						      value-index)
					  randomize-fraction)))))

  (let ((input-values-port (open-input-file input-values-file))
	(output-values-port (open-output-file output-values-file)))
    ;; Copy the header line from the input port to the output
    ;; port.
    (display (get-line input-values-port)
	     output-values-port)
    (newline output-values-port)

    ;; Read a single vector of values from the input port.  We randomize
    ;; these many times to generate the body of the output.
    (let ((values-vector (get-line-vector input-values-port)))
      ;; We're done with the input file.
      (close-port input-values-port)
      ;; Do the many randomizations, writing each as a tab-delimited line
      ;; in the output file.
      (do ((randomize-index 0
			    (+ randomize-index
			       1)))
	  ((= randomize-index
	      randomize-count)
	   #f)
	(put-line-vector (randomize-values-vector values-vector
						  randomize-index)
			 output-values-port)))
    (close-port output-values-port)))


;;; For now, this will only work with populations, not a mixture of
;;; populations and rates.
(define (randomize-pops input-pops-file
			output-pops-file
			randomize-count
			randomize-fraction)
  (randomize-values input-pops-file
		    output-pops-file
		    gamma-randomize-integer
		    randomize-count
		    randomize-fraction))

;;; For now, this will only work with rates, not a mixture of rates and
;;; populations.
(define (randomize-rates input-rates-file
			 output-rates-file
			 randomize-count
			 randomize-fraction)
  (randomize-values input-rates-file
		    output-rates-file
		    gamma-randomize-real
		    randomize-count
		    randomize-fraction))
