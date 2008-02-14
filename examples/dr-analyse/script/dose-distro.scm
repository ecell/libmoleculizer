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

;; Read the stats file, converting pertinent columns to vectors, and
;; getting the low and hi doses for purposes of tabulation.
(define (make-stats-cuts stats-port
			 return-the-values)
  (let ((log-dose-vector (make-vector))
	(dose-column-index 0)
	(response-vector (make-vector))
	(response-column-index 1)
	(variance-vector (make-vector))
	(variance-column-index 2)
	(log-dose-lo #f)
	(log-dose-hi #f))
    ;; Do I have a canned routine to convert a tab-delimited file into
    ;; cut vectors?
    (do ((stats-vector (get-line-vector stats-port)
		       (get-line-vector stats-port)))
	((not stats-vector)
	 (return-the-values log-dose-vector
			    response-vector
			    variance-vector
			    log-dose-lo
			    log-dose-hi))
      (let ((log-dose (log (vector-ref stats-vector
				       dose-column-index)))
	    (response (vector-ref stats-vector
				  response-column-index))
	    (variance (vector-ref stats-vector
				  variance-column-index)))
	;; Push the statistics onto their vectors.
	(vector-push! log-dose-vector
		      log-dose)
	(vector-push! response-vector
		      response)
	(vector-push! variance-vector
		      variance)
	;; Adjust the maximum and minimum doses.
	(and (or (not log-dose-lo)
		 (< log-dose
		    log-dose-lo))
	     (set! log-dose-lo
		   log-dose))
	(and (or (not log-dose-hi)
		 (> log-dose
		    log-dose-hi))
	     (set! log-dose-hi
		   log-dose))))))

;;; Evidently, I need a basic routine to do this, which is very similar
;;; to the GSL interface in C.
(define (make-stat-interpolator log-dose-vector
				stat-vector)
  (let ((dose-count (vector-length log-dose-vector)))
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
				 (vector-ref log-dose-vector
					     dose-index)
				 (vector-ref stat-vector
					     dose-index))))))

(define (log-tabulate-function output-port
			       function
			       log-arg-lo
			       log-arg-hi
			       arg-count)
  ;; I seem to be having trouble with getting out of range, probably due to
  ;; (exp (log x)) not being x.
  (let ((log-arg-delta (/ (- log-arg-hi
			     log-arg-lo)
			  arg-count)))
    (do ((arg-index 0
		    (+ arg-index
		       1))
	 ;; I seem to be having trouble with getting out of range, probably
	 ;; due to (exp (log x)) not being x.
	 (log-arg log-arg-lo
		  (+ log-arg
		     log-arg-delta)))
	((= arg-index
	    arg-count)
	 output-port)
      (put-line-vector (vector (exp log-arg)
			       (function log-arg))
		       output-port))))


(define (log-tabulate-interpolated-stats stats-file
					 dr-file
					 var-file
					 deriv-file
					 tabulate-point-count)
  ;; Convert the stats file into cut vectors, and get the domain
  ;; interval.
  (let ((log-dose-vector #f)
	(response-vector #f)
	(variance-vector #f)
	(log-dose-lo #f)
	(log-dose-hi #f))
    (let ((stats-port (open-input-file stats-file)))
      (make-stats-cuts stats-port
		       (lambda (ld-vector
				r-vector
				v-vector
				ld-lo
				ld-hi)
			 (set! log-dose-vector
			       ld-vector)
			 (set! response-vector
			       r-vector)
			 (set! variance-vector
			       v-vector)
			 (set! log-dose-lo
			       ld-lo)
			 (set! log-dose-hi
			       ld-hi)))
      (close-port stats-port))

    ;; Construct the response interpolator, which we use to tabulate both
    ;; the dose-response curve and its derivative.
    (let ((dr-interpolator (make-stat-interpolator log-dose-vector
						   response-vector)))

      ;; Tabulate the Bezier-interpolated dose-response curve.
      (let ((dr-port (open-output-file dr-file)))
	(log-tabulate-function dr-port
			   (lambda (log-dose)
			     (interpolator-eval dr-interpolator
						log-dose))
			   log-dose-lo
			   log-dose-hi
			   tabulate-point-count)
	(close-port dr-port))

      ;; Tabulate the derivative of the dose-response curve.
      (let ((deriv-port (open-output-file deriv-file)))
	(log-tabulate-function deriv-port
			   (lambda (log-dose)
			     (interpolator-eval-derivative dr-interpolator
							   log-dose))
			   log-dose-lo
			   log-dose-hi
			   tabulate-point-count)
	(close-port deriv-port)))

    ;; Construct the variance interpolator, and do its tabulation.
    (let ((variance-interpolator (make-stat-interpolator log-dose-vector
							 variance-vector)))
      (let ((var-port (open-output-file var-file)))
	(log-tabulate-function var-port
			   (lambda (log-dose)
			     (interpolator-eval variance-interpolator
						log-dose))
			   log-dose-lo
			   log-dose-hi
			   tabulate-point-count)
	(close-port var-port)))))

(define (doit)
  (let ((options (get-options))
	(tabulate-point-count 200))
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