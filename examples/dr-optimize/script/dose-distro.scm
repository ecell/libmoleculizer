
(load "utl.scm")
(load "gsl.scm")

;; Read the stats file, converting pertinent columns to vectors, and
;; getting the low and hi doses for purposes of tabulation.
(define (make-stats-cuts stats-port
			 return-the-values)
  (let ((dose-vector (make-vector))
	(dose-column-index 0)
	(response-vector (make-vector))
	(response-column-index 1)
	(variance-vector (make-vector))
	(variance-column-index 2)
	(dose-lo #f)
	(dose-hi #f))
    ;; Do I have a canned routine to convert a tab-delimited file into
    ;; cut vectors?
    (do ((stats-vector (get-line-vector stats-port)
		       (get-line-vector stats-port)))
	((not stats-vector)
	 (return-the-values dose-vector
			    response-vector
			    variance-vector
			    dose-lo
			    dose-hi))
      (let ((dose (vector-ref stats-vector
			      dose-column-index))
	    (response (vector-ref stats-vector
				  response-column-index))
	    (variance (vector-ref stats-vector
				  variance-column-index)))
	;; Push the statistics onto their vectors.
	(vector-push! dose-vector
		      dose)
	(vector-push! response-vector
		      response)
	(vector-push! variance-vector
		      variance)
	;; Adjust the maximum and minimum doses.
	(and (or (not dose-lo)
		 (< dose
		    dose-lo))
	     (set! dose-lo
		   dose))
	(and (or (not dose-hi)
		 (> dose
		    dose-hi))
	     (set! dose-hi
		   dose))))))

;;; Evidently, I need a basic routine to do this, which is very similar
;;; to the GSL interface in C.
(define (make-stat-interpolator dose-vector
				stat-vector)
  (let ((dose-count (vector-length dose-vector)))
    (let ((stat-interpolator (make-interpolator dose-count)))
      ;; Not really necessary to insert these in order, since
      ;; init! now sorts as needed.
      (do ((dose-index 0
		       (+ dose-index
			  1)))
	  ((= dose-index
	      dose-count)
	   (begin
	     (interpolator-init! stat-interpolator)
	     stat-interpolator))
	(interpolator-set-point! stat-interpolator
				 dose-index
				 (vector-ref dose-vector
					     dose-index)
				 (vector-ref stat-vector
					     dose-index))))))

;; Tablulates function to output port.  Doesn't close port.
(define (tabulate-function output-port
			   function
			   arg-lo
			   arg-hi
			   arg-count)
  (let ((arg-delta (/ (- arg-hi
			 arg-lo)
		      arg-count)))
    (do ((arg-index 0
		    (+ arg-index
		       1))
	 (arg arg-lo
	      (+ arg
		 arg-delta)))
	((= arg-index
	    arg-count)
	 output-port)
      (put-line-vector (vector arg
			       (function arg))
		       output-port))))

;; Tabulates, but with logarithmically spaced arguments.  This is to prevent
;; there being too little coverage of the low range in semi-log plots.  This
;; should make the plotted points equally spaced on the horizontal,
;; logarithmic axis.
(define (log-tabulate-function output-port
			       function
			       arg-lo
			       arg-hi
			       arg-count)
  ;; I seem to be having trouble with getting out of range, probably due to
  ;; (exp (log x)) not being x.
  (let ((log-arg-delta (* 0.999 (/ (- (log arg-hi)
				      (log arg-lo))
				   arg-count))))
    (do ((arg-index 0
		    (+ arg-index
		       1))
	 ;; I seem to be having trouble with getting out of range, probably
	 ;; due to (exp (log x)) not being x.
	 (log-arg (* 1.001 (log arg-lo))
		  (+ log-arg
		     log-arg-delta)))
	((= arg-index
	    arg-count)
	 output-port)
      (put-line-vector (let ((arg (exp log-arg)))
			 (display "arg = ")
			 (display arg)
			 (newline)
			 (vector arg
				 (function arg)))
		       output-port))))

(define (tabulate-interpolated-stats stats-file
				     dr-file
				     var-file
				     deriv-file
				     tabulate-point-count)
  ;; Convert the stats file into cut vectors, and get the domain
  ;; interval.
  (let ((dose-vector #f)
	(response-vector #f)
	(variance-vector #f)
	(dose-lo #f)
	(dose-hi #f))
    (let ((stats-port (open-input-file stats-file)))
      (make-stats-cuts stats-port
		       (lambda (d-vector
				r-vector
				v-vector
				lo
				hi)
			 (set! dose-vector
			       d-vector)
			 (set! response-vector
			       r-vector)
			 (set! variance-vector
			       v-vector)
			 (set! dose-lo
			       lo)
			 (set! dose-hi
			       hi)))
      (close-port stats-port))

    ;; Construct the response interpolator, which we use to tabulate both
    ;; the dose-response curve and its derivative.
    (let ((dr-interpolator (make-stat-interpolator dose-vector
						   response-vector)))

      ;; Tabulate the Bezier-interpolated dose-response curve.
      (let ((dr-port (open-output-file dr-file)))
	(tabulate-function dr-port
			   (lambda (dose)
			     (interpolator-eval dr-interpolator
						dose))
			   dose-lo
			   dose-hi
			   tabulate-point-count)
	(close-port dr-port))

      ;; Tabulate the derivative of the dose-response curve.
      (let ((deriv-port (open-output-file deriv-file)))
	(tabulate-function deriv-port
			   (lambda (dose)
			     (interpolator-eval-derivative dr-interpolator
							   dose))
			   dose-lo
			   dose-hi
			   tabulate-point-count)
	(close-port deriv-port)))

    ;; Construct the variance interpolator, and do its tabulation.
    (let ((variance-interpolator (make-stat-interpolator dose-vector
							 variance-vector)))
      (let ((var-port (open-output-file var-file)))
	(tabulate-function var-port
			   (lambda (dose)
			     (interpolator-eval variance-interpolator
						dose))
			   dose-lo
			   dose-hi
			   tabulate-point-count)
	(close-port var-port)))))

(define (log-tabulate-interpolated-stats stats-file
					 dr-file
					 var-file
					 deriv-file
					 tabulate-point-count)
  ;; Convert the stats file into cut vectors, and get the domain
  ;; interval.
  (let ((dose-vector #f)
	(response-vector #f)
	(variance-vector #f)
	(dose-lo #f)
	(dose-hi #f))
    (let ((stats-port (open-input-file stats-file)))
      (make-stats-cuts stats-port
		       (lambda (d-vector
				r-vector
				v-vector
				lo
				hi)
			 (set! dose-vector
			       d-vector)
			 (set! response-vector
			       r-vector)
			 (set! variance-vector
			       v-vector)
			 (set! dose-lo
			       lo)
			 (set! dose-hi
			       hi)))
      (close-port stats-port))

    ;; Construct the response interpolator, which we use to tabulate both
    ;; the dose-response curve and its derivative.
    (let ((dr-interpolator (make-stat-interpolator dose-vector
						   response-vector)))

      ;; Tabulate the Bezier-interpolated dose-response curve.
      (let ((dr-port (open-output-file dr-file)))
	(log-tabulate-function dr-port
			   (lambda (dose)
			     (interpolator-eval dr-interpolator
						dose))
			   dose-lo
			   dose-hi
			   tabulate-point-count)
	(close-port dr-port))

      ;; Tabulate the derivative of the dose-response curve.
      (let ((deriv-port (open-output-file deriv-file)))
	(log-tabulate-function deriv-port
			   (lambda (dose)
			     (interpolator-eval-derivative dr-interpolator
							   dose))
			   dose-lo
			   dose-hi
			   tabulate-point-count)
	(close-port deriv-port)))

    ;; Construct the variance interpolator, and do its tabulation.
    (let ((variance-interpolator (make-stat-interpolator dose-vector
							 variance-vector)))
      (let ((var-port (open-output-file var-file)))
	(log-tabulate-function var-port
			   (lambda (dose)
			     (interpolator-eval variance-interpolator
						dose))
			   dose-lo
			   dose-hi
			   tabulate-point-count)
	(close-port var-port)))))

(define (doit)
  (let ((options (get-options))
	(tabulate-point-count 100))
    (let (;; The file where experimental/simulation statistics are found.
	  (stats-file-pair (string-map-ref options
					   "-stats-file"))
	  ;; Output file for tabulated, interpolated dose-response.
	  (dr-file-pair (string-map-ref options
					"-dr-file"))
	  ;; Output file for tabulated, interpolated dose variances.
	  (var-file-pair (string-map-ref options
					 "-var-file"))
	  ;; Output file for tabulated derivative of dose-response.
	  (deriv-file-pair (string-map-ref options
					   "-deriv-file")))

      ;; Check that we really were given all the filenames we need.
      (or (pair? stats-file-pair)
	  (thread-error
	   "dose-distro.scm: no -stats-file option on command line."))
      (or (pair? dr-file-pair)
	  (thread-error
	   "dose-distro.scm: no -dr-file option on command line."))
      (or (pair? var-file-pair)
	  (thread-error
	   "dose-distro.scm: no -var-file option on command line."))
      (or (pair? deriv-file-pair)
	  (thread-error
	   "dose-distro.scm: no -deriv-file option on command line."))

      (log-tabulate-interpolated-stats (cdr stats-file-pair)
				       (cdr dr-file-pair)
				       (cdr var-file-pair)
				       (cdr deriv-file-pair)
				       tabulate-point-count))))

(doit)