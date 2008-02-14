
(load "utl.scm")
(load "gsl.scm")

;;; Transinformation calculations using continuous dose distribution
;;; given by the (normalized) derivative of the dose/response function.

;;; I will probably have to make some adjustments so that all our evaluations
;;; of interpolated functions fall strictly within the domain of those
;;; functions which, as seen earlier, is not always easy.

(define (expected-value weight-function
			value-function
			arg-lo
			arg-delta
			arg-count)
  ;; Being a little anal here to avoid evaluating weight-function
  ;; more than once inside the loop.
  (let ((weighted-sum 0.0)
	(total-weight 0.0))
    (do ((arg-index 0
		    (+ arg-index
		       1))
	 (arg arg-lo
	      (+ arg
		 arg-delta)))
	((= arg-index
	    arg-count)
	 (/ weighted-sum
	    total-weight))
      (let ((weight-factor (* (weight-function arg)
			      arg-delta)))
	(set! weighted-sum
	      (+ weighted-sum
		 (* (value-function arg)
		    weight-factor)))
	(set! total-weight
	      (+ total-weight
		 weight-factor))))))
       
	  
(define (caching-transinfo dose-density-function
			   mean-response-function
			   response-variance-function
			   log-dose-lo
			   log-dose-delta
			   dose-count
			   response-lo
			   response-delta
			   response-count)

  ;; The density function of the response distribution conditional on the
  ;; given dose.
  (define (conditional-density response
			       log-dose)
    ;; Having lots of trouble with bad arguments to the gamma density
    ;; function.  These tests are to see if they're appearing at the level of
    ;; mean and variance, or if they're coming in the calculation of the gamma
    ;; parameters from the mean and the variance.
    (or (< 0.0
	   response)
	(thread-error 
	 "conditional-density: response is not strictly positive."))
    (let ((mean-response (mean-response-function log-dose))
	  (response-variance (response-variance-function log-dose)))
      (or (< 0.0
	     mean-response)
	  (let ((msg-port (open-output-string))
		(dose (exp log-dose)))
	    (display "conditional-density: mean response "
		     msg-port)
	    (display mean-response
		     msg-port)
	    (display " to dose "
		     msg-port)
	    (display dose
		     msg-port)
	    (display " is not strictly positive."
		     msg-port)
	    (thread-error (close-port msg-port))))
      (or (< 0.0
	     response-variance)
	  (let ((msg-port (open-output-string))
		(dose (exp log-dose)))
	    (display "conditional-density: variance "
		     msg-port)
	    (display response-variance
		     msg-port)
	    (display " of response to dose "
		     msg-port)
	    (display dose
		     msg-port)
	    (display " is not strictly positive."
		     msg-port)
	    (thread-error (close-port msg-port))))
      (gamma-density response
		     mean-response
		     response-variance)))

  ;; Calculation of the predictive density function.  We want to cache values
  ;; returned by this function for use in the
  (define (predictive-density response)
    (expected-value dose-density-function
		    (lambda (log-dose)
		      (conditional-density response
					   log-dose))
		    log-dose-lo
		    log-dose-delta
		    dose-count))

  ;; Tabulate the predictive density function, a minor acceleration.
  ;; 
  ;; A new twist here is to use an interpolator, rather than a vector, to hold
  ;; the cache, so that I use generic integration routine.
  (display "Caching predictive density.")
  (newline)
  (define fast-predictive-density
    (let ((cache (make-interpolator response-count)))
      (do ((response-index 0
			   (+ response-index
			      1))
	   (response response-lo
		     (+ response
			response-delta)))
	  ((= response-index
	      response-count)
	   (begin 
	     (interpolator-init! cache)
	     (lambda (response)
	       (interpolator-eval cache
				  response))))
	(interpolator-set-point! cache
				 response-index
				 response
				 (predictive-density response)))))

  ;; Compute the Kullback-Leibler distance from the response distribution
  ;; given the dose to the predictive response distribution.
  ;;
  ;; One inefficiency here is computing the (gamma) conditional-density twice
  ;; every time through the loop.  But this lets us use a generic expected
  ;; value.
  (define (kld log-dose)
    (expected-value (lambda (response)
		      (conditional-density response
					   log-dose))
		    (lambda (response)
		      (let ((cd (conditional-density response
						     log-dose))
			    (pd (fast-predictive-density response)))
			;; These densities sometimes round off to 0.
			(if (and (< 0.0
				    cd)
				 (< 0.0
				    pd))
			    (- (log cd)
			       (log pd))
			    0.0)))
		    response-lo
		    response-delta
		    response-count))

  ;; Now the transinformation is the expected value with respect to dose of
  ;; the Kullback-Leibler distance above.
  ;;
  ;; Must not forget to convert from nits to bits.
  (display "Taking expected value of Kullback-Leibler distance.")
  (newline)
  (/ (expected-value dose-density-function
		     kld
		     log-dose-lo
		     log-dose-delta
		     dose-count)
     (log 2)))

;; Read the stats file, converting pertinent columns to vectors, and
;; getting the low and hi doses for purposes of tabulation.
;;
;; The bin-width is the width of the histogram bins. I use it to generate
;; a fake mean response and response-variance when either of these is
;; legitimately 0.
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

(define (caching-transinfo-from-port stats-port
				     log-dose-lo
				     log-dose-delta
				     dose-count
				     response-lo
				     response-delta
				     response-count
				     bin-width)
  ;; Cut the columns of the file as vectors.
  (display "Cutting columns of stats file.")
  (newline)
  (let ((dose-vector #f)
	(response-vector #f)
	(variance-vector #f)
	(dose-lo #f)
	(dose-hi #f))
    (make-stats-cuts stats-port
		     (lambda (dv
			      rv
			      vv
			      dl
			      dh)
		       (set! dose-vector
			     dv)
		       (set! response-vector
			     rv)
		       (set! variance-vector
			     vv)
		       (set! dose-lo
			     dl)
		       (set! dose-hi
			     dh)))
    ;; Construct interpolators for the functions given by the columns.
    (display "Constructing interpolators.")
    (newline)
    (let ((response-interpolator (make-stat-interpolator dose-vector
							 response-vector))
	  (variance-interpolator (make-stat-interpolator dose-vector
							 variance-vector)))
      (let ((min-response (/ bin-width
			     2.0))
	    (min-variance (/ (* bin-width
				bin-width)
			     4.0)))
	(let ((dose-density-function
	       (lambda (log-dose)
		 ;; For the (non-normalized) weighting function on doses, we
		 ;; use the derivative of the response function, as we have
		 ;; discussed in the seminar.
		 ;; 
		 ;; Sometimes the derivative of the interpolated dose/response
		 ;; curve is not strictly positive; the dose/response curve
		 ;; is sometimes not monotone in the regions where it is nearly
		 ;; horizontal.
		 ;;
		 ;; It's okay for this to be flat 0.
		 (let ((dose (exp log-dose)))
		   ;; Here is where the Jacobian of the exponential function
		   ;; comes in.
		   (let ((dose-density
			  (* dose
			     (interpolator-eval-derivative
			      response-interpolator
			      dose))))
		     (if (< 0.0
			    dose-density)
			 dose-density
			 0.0)))))
	      (mean-response-function
	       (lambda (log-dose)
		 ;; Sometimes responses are flat 0.
		 ;; 
		 ;; Sometimes, when responses are nearly 0, the interpolating
		 ;; function will go slightly negative. E.g. before starting
		 ;; the "main ascent" on the dose-response curve.
		 (let ((response (interpolator-eval response-interpolator
						    (exp log-dose))))
		   (if (< min-response
			  response)
		       response
		       min-response))))
	      (response-variance-function
	       (lambda (log-dose)
		 ;; Sometimes variance of the response distribution is 0,
		 ;; typically when the responses are all flat 0.
		 (let ((variance (interpolator-eval variance-interpolator
						    (exp log-dose))))
		   (if (< min-variance
			  variance)
		       variance
		       min-variance)))))
	  ;; Calculate the transinformation.
	  (display "Calculating transinformation.")
	  (newline)
	  (caching-transinfo dose-density-function
			     mean-response-function
			     response-variance-function
			     log-dose-lo
			     log-dose-delta
			     dose-count
			     response-lo
			     response-delta
			     response-count))))))

;;; What this file does when loaded.
(define (doit)
    ;; Get the command-line options; basically input and output files.
    (let ((options (get-options)))
      (let (;; The file where experimental/simulation statistics are found.
	    (stats-file-pair (string-map-ref options
					     "-stats-file"))
	    ;; Output file for transinformation.
	    (transinfo-file-pair (string-map-ref options
						 "-transinfo-file")))

	;; Check that we really were given all the filenames we need.
	(or (pair? stats-file-pair)
	    (thread-error
	     "dose-distro.scm: no -stats-file option on command line."))
	(or (pair? transinfo-file-pair)
	    (thread-error
	     "dose-distro.scm: no -transinfo-file option on command line."))

	;; Set up the wired-in options.
	;;
	;; It will probably be necessary to treat doses logarithmically, if
	;; only to hold down the size of the integral.
	;;
	;; As I saw in "dose-distro.scm", using logarithmic doses causes
	;; problems with getting somewhat outside the range of the
	;; interpolated functions.  Here, using dose-inset-epsilon, we retreat
	;; to a slightly smaller range of doses.
	(let ((dose-inset-epsilon 1.0)
	      (dose-lo 25)
	      (dose-hi 102400)
	      (dose-count 200)
	      (response-lo 10.0)
	      (response-hi 1600)
	      (response-count 400)
	      ;; This is the histogramming bin width.  I use it to generate
	      ;; a fake response (bin-width/2) and fake response variance
	      ;; when either of these turns out to be 0, as they sometimes
	      ;; legitimately do.
	      (bin-width 5))
	  (let ((log-dose-lo (+ (log dose-lo)
				dose-inset-epsilon))
		(log-dose-hi (- (log dose-hi)
				dose-inset-epsilon)))
	    (let ((log-dose-delta (/ (- log-dose-hi
					log-dose-lo)
				     dose-count))
		  (response-delta (/ (- response-hi
					response-lo)
				     response-count)))
	      (let ((stats-port (open-input-file (cdr stats-file-pair)))
		    (transinfo-port (open-output-file (cdr transinfo-file-pair))))
		(let ((transinfo (caching-transinfo-from-port stats-port
							      log-dose-lo
							      log-dose-delta
							      dose-count
							      response-lo
							      response-delta
							      response-count
							      bin-width)))
		  ;; Write message to console.
		  (display "The transinformation is ")
		  (display transinfo)
		  (display " bits.")
		  (newline)

		  ;; Record the transinformation in the output file for later
		  ;; use.
		  (display transinfo
			   transinfo-port)
		  (newline transinfo-port))
		(close-port stats-port)
		(close-port transinfo-port))))))))


(doit)