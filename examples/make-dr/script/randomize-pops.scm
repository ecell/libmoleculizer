
(load "script/randomize.scm")

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

