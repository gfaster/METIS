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
      (define i (symbol-value (car (lookup-symbol i-str))))
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
	(opt "METIS_OPTION_PTYPE" "mptype_et")
	(opt "METIS_OPTION_OBJTYPE" "mobjtype_et")
	(opt "METIS_OPTION_CTYPE" "mctype_et")
	(opt "METIS_OPTION_IPTYPE" "miptype_et")
	(opt "METIS_OPTION_RTYPE" "mrtype_et")
	(opt "METIS_OPTION_CONTIG" "bool")
	(opt "METIS_OPTION_DBGLVL" #f)
	(opt "METIS_OPTION_NIPARTS" #f)
	(opt "METIS_OPTION_NITER" #f)
	(opt "METIS_OPTION_NCUTS" #f)
	(opt "METIS_OPTION_SEED" #f)
	(opt "METIS_OPTION_ONDISK" "bool")
	(opt "METIS_OPTION_MINCONN" #f)
	(opt "METIS_OPTION_CONTIG" "bool")
	(opt "METIS_OPTION_COMPRESS" "bool")
	(opt "METIS_OPTION_CCORDER" "bool")
	(opt "METIS_OPTION_PFACTOR" #f)
	(opt "METIS_OPTION_NSEPS" #f)
	(opt "METIS_OPTION_UFACTOR" #f)
	(opt "METIS_OPTION_NUMBERING" #f)
	(opt "METIS_OPTION_DROPEDGES" #f)
	(opt "METIS_OPTION_NO2HOP" "bool")
	(opt "METIS_OPTION_TWOHOP" "bool")
	(opt "METIS_OPTION_FAST" "bool"))))
  (if enabled? (make-pretty-printer-worker "map" (lambda (_) #f) children) #f))

(define metis-options
  (make-pretty-printer "metis-options" 
		       (lambda (self o) (metis-options-fn o #f))))

(set-pretty-printers! (cons metis-options
			(pretty-printers)))
