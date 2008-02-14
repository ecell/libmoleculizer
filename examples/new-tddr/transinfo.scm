
(load "utl.scm")
(load "gsl.scm")
(load "bayes2.scm")

(define (riemann-sum f
		     x-lo
		     x-delta
		     x-count)
  (do ((the-sum 0.0
		(+ the-sum
		   (* x-delta
		      (f x))))
       (x x-lo
	  (+ x
	     x-delta))
       (x-index 0
		(+ x-index
		   1)))
      ((= x-index
	  x-count)
       the-sum)))

;;; Here, f should be a function of the index into the weights-vector,
;;; which is taken to be a scaled probability distribution on its
;;; index.
(define (expected-value-discrete weights-vector
				 f)
  (do ((scaled-sum 0.0
		   (+ scaled-sum
		      (* (vector-ref weights-vector
				     weight-index)
			 (f weight-index))))
       (total-weight 0.0
		     (+ total-weight
			(vector-ref weights-vector
				    weight-index)))
       (weight-index 0
		     (+ weight-index
			1)))
      ((= weight-index
	  (vector-length weights-vector))
       (/ scaled-sum
	  total-weight))))

;;; For calculating dose/response transinformation using gamma distributions
;;; fitted to conditional responses to a finite list of doses.
(define (gamma-transinfo dose-weights-vector
			 mean-responses-vector
			 response-variances-vector
			 x-lo
			 x-hi
			 x-count)
  ;; Just to make invocation easier.
  (define (gd x
	      dose-index)
    (gamma-density x
		   (vector-ref mean-responses-vector
			       dose-index)
		   (vector-ref response-variances-vector
			       dose-index)))
  ;; The predictive response density at a particular value x of the
  ;; response variable.
  (define (predictive-density x)
    (expected-value-discrete dose-weights-vector
			     (lambda (dose-index)
			       (gd x
				   dose-index))))
  ;; Takess the Kullback-Leibler distance from the response distribution
  ;; conditional on dose to the predictive distribution.
  (define (kld dose-index)
    (let ((x-range-width (- x-hi
			    x-lo)))
      (let ((x-delta (/ x-range-width
			x-count)))
	(riemann-sum (lambda (x)
		       (let ((x-density (gd x
					    dose-index)))
			 (* x-density
			    (log (/ x-density
				    (predictive-density x))))))
		     x-lo
		     x-delta
		     x-count))))
  ;; And now the transinformation is the expected value of the
  ;; Kullback-Leibler distance.
  (expected-value-discrete dose-weights-vector
			   kld))

;;; Note that this function "stubs off" the weight vector, assuming
;;; that all weights are 1.  I will want to bring the weight vector
;;; forward for information capacity calculations.
(define (gamma-transinfo-file stats-file
			      x-lo
			      x-hi
			      x-count)
  (let ((dose-weights-vector (make-vector))
	(mean-responses-vector (make-vector))
	(response-variances-vector (make-vector))
	(stats-port (open-input-file stats-file)))
    (do ((stats-finished #f
			 stats-finished))
	(stats-finished
	 (begin (close-port stats-port)
		(gamma-transinfo dose-weights-vector
				 mean-responses-vector
				 response-variances-vector
				 x-lo
				 x-hi
				 x-count)))
      ;; We don't use the dump header here.
      (let ((dmp-header (read stats-port)))
	(if (eof-object? dmp-header)
	    (set! stats-finished
		  #t)
	    (let ((mean (read stats-port)))
	      (and (eof-object? mean)
		   (thread-error
		    "transinfo.scm: mean missing at end of file."))
	      (let ((variance (read stats-port)))
		(and (eof-object? variance)
		     (thread-error
		      "transinfo.scm: variance missing at end of file."))
		;; For now, all the weights are equal.
		(vector-push! dose-weights-vector 1.0)
		(vector-push! mean-responses-vector mean)
		(vector-push! response-variances-vector variance))))))))

;;; Get the stats file name from command line options, and takes the
;;; transinformation as described above.
(let ((options (get-options)))
  (let ((stats-file-pair (string-map-ref options
					 "-stats-file"))
	(transinfo-file-pair (string-map-ref options
					     "-transinfo-file")))
    (or stats-file-pair
	(thread-error
	 "transinfo.scm: no -stats-file option on command line."))
    (or transinfo-file-pair
	(thread-error
	 "transinfo.scm: no -transinfo-file option on command line."))
    ;; These numbers will generally have to come from eyeballing
    ;; plots, I expect.  That's how I got these.
    (let ((transinfo (let ((x-lo 10.0)
			   (x-hi 300)
			   (x-count 600))
		       (* bits-per-nit
			  (gamma-transinfo-file (cdr stats-file-pair)
						x-lo
						x-hi
						x-count))))
	  (transinfo-port (open-output-file (cdr transinfo-file-pair))))
      ;; Write message to console.
      (display "The transinformation is ")
      (display transinfo)
      (display " bits.")
      (newline)
      ;; Record transinformation for later use.
      (display transinfo
	       transinfo-port)
      (newline transinfo-port)
      (close-port transinfo-port))))
    
