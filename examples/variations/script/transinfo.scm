;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Nu - a C++ friendly Scheme byte-code compiler.
;;; Copyright (C) 2004  Walter Lawrence (Larry) Lok.
;;;
;;; This program is free software; you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation; either version 2 of the License, or
;;; (at your option) any later version.

;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.

;;; You should have received a copy of the GNU General Public License
;;; along with this program; if not, write to the Free Software
;;; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
;;;    
;;; Contact information:
;;;   Larry Lok, Research Fellow          Voice: 510-981-8740
;;;   The Molecular Sciences Institute      Fax: 510-647-0699
;;;   2168 Shattuck Ave.                  Email: lok@molsci.org
;;;   Berkeley, CA 94704
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(load "utl.scm")
(load "gsl.scm")
(load "bayes2.scm")

;;; This version of transinfo.scm tries out cacheing of the predictive density
;;; in a vector, so that it is not recomputed over and over.  But this means
;;; that I don't get to use generic Riemann sum.  I could, if I replaced the
;;; very simple cache vector used here with an interpolator.

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
  (define x-delta
    (let ((x-range-width (- x-hi
			    x-lo)))
      (/ x-range-width
	 x-count)))

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

  ;; Rather than do this sum over and over again, we cache these
  ;; values.
  (define predictive-density-cache
    (let ((cache (make-vector x-count)))
      (do ((cache-index 0
			(+ cache-index
			   1))
	   (x x-lo
	      (+ x
		 x-delta)))
	  ((= cache-index
	      x-count)
	   cache)
	(vector-set! cache
		     cache-index
		     (predictive-density x)))))
		 
  ;; Takess the Kullback-Leibler distance from the response distribution
  ;; conditional on dose to the predictive distribution.
  (define (kld dose-index)
    ;; We can't do a generic Riemann sum here, since I want to try
    ;; cacheing of the predictive density.
    (do ((the-sum 0.0
		  (+ the-sum
		     (* x-delta
			(let ((x-density (gd x
					     dose-index))
			      (pred-density (vector-ref
					     predictive-density-cache
					     x-index)))
			  ;; Sometimes the predictive density will
			  ;; underflow giving zero numerator but
			  ;; non-zero denominator.  To avoid this, we
			  ;; just return 0 when either density is 0.
			  (if (and (< 0.0
				      x-density)
				   (< 0.0
				      pred-density))
			      (* x-density
				 (log (/ x-density
					 pred-density)))
			      0.0)))))
	 (x x-lo
	    (+ x
	       x-delta))
	 (x-index 0
		  (+ x-index
		     1)))
	((= x-index
	    x-count)
	 the-sum)))

  ;; And now the transinformation is the expected value of the
  ;; Kullback-Leibler distance.
  (expected-value-discrete dose-weights-vector
			   kld))

;;; Note that this function "stubs off" the weight vector, assuming
;;; that all weights are 1.  I will want to bring the weight vector
;;; forward for information capacity calculations.
;;;
;;; Note that the bin-size is the size of the bins in histogramming.  This
;;; default for variance (when sample variance comes out 0) is what I used in
;;; gamma-density.scm as a way of making the gamma distribution look as much
;;; as possible like the histogram distribution.
(define (gamma-transinfo-file stats-file
			      x-lo
			      x-hi
			      x-count
			      bin-size)
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
      ;;
      ;; Why am I not using get-tab-vector here?
      (let ((dose (read stats-port)))
	(if (eof-object? dose)
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
		;; Read and ignore the standard error.
		(let ((standard-error (read stats-port)))
		  (and (eof-object? standard-error)
		       (thread-error
			"transinfo.scm: standard-error missing at end of file."))
		  ;; Default the variance, if 0, to make the standard
		  ;; deviation equal to half the histogram bin-size.
		  (or (< 0.0
			 variance)
		      (set! variance (/ (* bin-size
					   bin-size)
					4.0)))
		  ;; Default the mean, if 0, to the middle of the lowest
		  ;; histogram bin.
		  (or (< 0.0
			 mean)
		      (set! mean (/ bin-size
				    2.0)))
		  ;; For now, all the weights are equal.
		  (vector-push! dose-weights-vector 1.0)
		  (vector-push! mean-responses-vector mean)
		  (vector-push! response-variances-vector variance)))))))))

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
    ;;
    ;; Here, bin-size is the size of the histogram bins, and I use it here to
    ;; make the gamma distribution look like the histogram distribution as in
    ;; gamma-density.scm.
    (let ((transinfo (let ((x-lo 10.0)
			   (x-hi 1600)
			   (x-count 800)
			   (bin-size 5))
		       (* bits-per-nit
			  (gamma-transinfo-file (cdr stats-file-pair)
						x-lo
						x-hi
						x-count
						bin-size))))
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
    
