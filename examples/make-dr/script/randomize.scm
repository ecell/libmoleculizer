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

(define (gamma-randomize-real positive-real
			      randomize-fraction)
  (let ((desired-sd (* positive-real
		       randomize-fraction)))
    (let ((desired-variance (* desired-sd
			       desired-sd)))
      (sample-gamma positive-real
		    desired-variance))))

(define (make-gamma-real-randomizer randomize-fraction)
  (lambda (positive-real)
    (gamma-randomize-real positive-real
			  randomize-fraction)))

(define (make-correlated-real-randomizer randomize-fraction)
  (let ((desired-variance (* randomize-fraction
			     randomize-fraction)))
    (let ((multiplier (sample-gamma 1
				    desired-variance)))
      (lambda (positive-real)
	(* positive-real
	   multiplier)))))

;;; Produces integers whose mean is roughly positive-real.
(define (gamma-randomize-integer positive-real
				 randomize-fraction)
  (floor (+ (gamma-randomize-real positive-real
				  randomize-fraction)
	    0.5)))

(define (make-gamma-integer-randomizer randomize-fraction)
  (lambda (positive-real)
    (gamma-randomize-integer positive-real
			     randomize-fraction)))

(define (make-correlated-integer-randomizer randomize-fraction)
  (let ((desired-variance (* randomize-fraction
			     randomize-fraction)))
    (let ((multiplier (sample-gamma 1
				    desired-variance)))
      (lambda (positive-real)
	(floor (+ (* positive-real
		     multiplier)
		  0.5))))))

;;; Copies header (which gives a sort of "variable name" to each column)
;;; from the input-values-file to the output-values-file.  Then randomizes
;;; the values on the first following row in input-values-file to generate
;;; all the remaining rows of the output-values-file
(define (randomize-values input-values-file
			  output-values-file
			  make-randomizer-function
			  randomize-count
			  randomize-fraction)
  ;; Function to randomize one row.
  ;;
  ;; When we are doing correlated randomization, we want each member of a row
  ;; to be multiplied by the same "randomizing factor," so we create a new
  ;; randomizer for each row.
  (define (randomize-values-vector values-vector
				   randomize-index)
    (let ((result-vector (make-vector))
	  (randomizer (make-randomizer-function randomize-fraction)))
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
		      (randomizer (vector-ref values-vector
					      value-index))))))

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
      ;; Emit the nominal value as the first randomization.
      (put-line-vector values-vector
		       output-values-port)
      ;; Do the many randomizations, writing each as a tab-delimited line
      ;; in the output file.
      (do ((randomize-index 1
			    (+ randomize-index
			       1)))
	  ((<= randomize-count
	       randomize-index)
	   #f)
	(put-line-vector (randomize-values-vector values-vector
						  randomize-index)
			 output-values-port)))
    (close-port output-values-port)))


;;; For now, this will only work with populations, not a mixture of
;;; populations and rates.
;;;
;;; Note that this is doing CORRELATED randomization, which is probably
;;; only suitable for populations (or concentrations).
(define (randomize-pops input-pops-file
			output-pops-file
			randomize-count
			randomize-fraction)
  (randomize-values input-pops-file
		    output-pops-file
		    make-correlated-integer-randomizer
		    randomize-count
		    randomize-fraction))

;;; For now, this will only work with rates, not a mixture of rates and
;;; populations.
;;;
;;; Note that this is doing CORRELATED randomization, which is probably
;;; only suitable for populations or concentrations.
(define (randomize-rates input-rates-file
			 output-rates-file
			 randomize-count
			 randomize-fraction)
  (randomize-values input-rates-file
		    output-rates-file
		    make-correlated-real-randomizer
		    randomize-count
		    randomize-fraction))
