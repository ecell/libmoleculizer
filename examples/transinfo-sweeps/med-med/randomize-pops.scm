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

(load "randomize.scm")

(let ((options (get-options)))
  (let ((input-file-pair (string-map-ref options
					 "-input"))
	(output-file-pair (string-map-ref options
					  "-output"))
	(randomize-count-pair (string-map-ref options
					      "-count"))
	(randomize-fraction-pair (string-map-ref options
						 "-fraction")))
    (or input-file-pair
	(thread-error "randomize-pops.scm: no -input option given."))
    (or output-file-pair
	(thread-error "randomize-pops.scm: no -output option given."))
    (or randomize-count-pair
	(thread-error "randomize-pops.scm: no -count option given."))
    (or randomize-fraction-pair
	(thread-error "randomize-pops.scm: no -fraction option given."))
    
    ;; The numerical arguments have to be read to get numbers.
    (let ((count (let ((count-port (open-input-string
				    (cdr randomize-count-pair))))
		   (let ((count (read count-port)))
		     (close-port count-port)
		     count)))
	  (fraction (let ((fraction-port (open-input-string
					  (cdr randomize-fraction-pair))))
		      (let ((fraction (read fraction-port)))
			(close-port fraction-port)
			fraction))))
      (randomize-pops (cdr input-file-pair)
		      (cdr output-file-pair)
		      count
		      fraction))))

