(use-modules (gdb)
	     (gdb iterator)
	     )

(define (try-car pair) (if pair (car pair) #f))

(define* 
  (metis-options-fn optarr #:optional (force #f))
  (define enabled?
    (and
      (equal? (value-type optarr) (type-pointer (lookup-type "idx_t")))
      (or
	force
	(let ((opt (try-car (lookup-symbol "options"))))
	  (and opt (equal? (symbol-value opt #:frame (selected-frame)) optarr))))))
  (define (children _)
    (define (opt i-str ty-str)
      (define i-sym-str (string-join (list "METIS_OPTION_" i-str) ""))
      (define i (symbol-value (car (lookup-symbol i-sym-str))))
      (define ty (if ty-str 
		   (let ((sym (try-car (lookup-symbol ty-str)))) 
			      (if sym (symbol-type sym) #f))
		   #f))
      (define indexed (value-subscript optarr i))
      (cons 
	i-str
	(if (= (value->integer indexed) -1)
	  "<unset>"
	  (if ty (value-cast indexed ty) indexed))))

    (make-list-iterator
      (list
	(opt "PTYPE" "mptype_et")
	(opt "OBJTYPE" "mobjtype_et")
	(opt "CTYPE" "mctype_et")
	(opt "IPTYPE" "miptype_et")
	(opt "RTYPE" "mrtype_et")
	(opt "CONTIG" "bool")
	(opt "DBGLVL" #f)
	(opt "NIPARTS" #f)
	(opt "NITER" #f)
	(opt "NCUTS" #f)
	(opt "SEED" #f)
	(opt "ONDISK" "bool")
	(opt "MINCONN" #f)
	(opt "CONTIG" "bool")
	(opt "COMPRESS" "bool")
	(opt "CCORDER" "bool")
	(opt "PFACTOR" #f)
	(opt "NSEPS" #f)
	(opt "UFACTOR" #f)
	(opt "NUMBERING" #f)
	(opt "DROPEDGES" #f)
	(opt "NO2HOP" "bool"))))
  (if enabled? (make-pretty-printer-worker "map" (lambda (_) #f) children) #f))

(define metis-options
  (make-pretty-printer "metis-options" 
		       (lambda (self o) (metis-options-fn o #f))))

(set-pretty-printers! (cons metis-options
			(pretty-printers)))
