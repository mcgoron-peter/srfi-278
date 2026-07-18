;;; SPDX-FileCopyrightText: 2026 Peter McGoron
;;; SPDX-License-Identifier: MIT

(define (nan? obj)
  (and (number? obj) (r7rs:nan? obj)))

(define (exact-integer? obj)
  (and (integer? obj) (exact? obj)))

(define (exact-zero? z)
  (and (zero? z) (exact? z)))

(define (strictly-real? obj)
  (and (complex? obj)
       (exact-zero? (imag-part obj))))

(define (strictly-rational? obj)
  (and (rational? obj)
       (exact-zero? (imag-part obj))))

(define (strictly-integer? obj)
  (and (integer? obj)
       (exact-zero? (imag-part obj))))

(define (conjugate z)
  ;; Return exact value given possibly exact arguments.
  ;; For implementations with only inexact complex numbers.
  (if (strictly-real? z)
      z
      (make-rectangular (real-part z)
                        (- (imag-part z)))))

(cond-expand
  ((or (library (srfi 144))
       (library (srfi 276)))
   (define (make-flonum* x k)
     (cond
       ((eqv? k +inf.0) +inf.0)
       ((eqv? k -inf.0) (* (sign x) 0.0))
       (else (make-flonum (flonum x) (exact k)))))
   (define (cssqs z)
     ;; NOTE: We don't have access to the overflow and underflow flags
     ;; here, so we cannot implement this part of Kahan's algorithm
     ;; verbatim. We need to do a trick.
     (let* ((x (real-part z))
            (y (imag-part z))
            (x^2 (square x))
            (y^2 (square y))
            (rho (+ x^2 y^2)))
       (if (and (or (nan? rho) (infinite? rho))
                (or (infinite? x) (infinite? y)))
           (values +inf.0 0)
           (let ((underflowed? (or (< x^2 fl-least)
                                   (< y^2 fl-least)))
                 (overflowed? (or (infinite? rho)
                                  (infinite? x^2)
                                  (infinite? y^2))))
             (if (or overflowed?
                     (and underflowed? (< rho (/ fl-least fl-epsilon))))
                 (let ((k (flexponent (flonum (max (abs x) (abs y))))))
                   ;; NOTE: This isn't scalb, this is ldexp
                   (values (+ (square (make-flonum* x (- k)))
                              (square (make-flonum* y (- k))))
                           k))
                 (values rho 0))))))
   ;; These are defined on all values, and work around a contradiction
   ;; in the R7RS-Small.
   (define (even*? k)
     (and (integer? k) (not (infinite? k)) (even? k)))
   (define (odd*? k)
     (and (integer? k) (not (infinite? k)) (odd? k)))
   (define (csqrt z)
     (let*-values (((x) (real-part z))
                   ((y) (imag-part z))
                   ((rho k) (cssqs z))
                   ((rho) (if (not (nan? x))
                              (+ (make-flonum* (abs x)
                                               (- k))
                                 (sqrt rho))))
                   ((rho) (if (even*? k)
                              (+ rho rho)
                              rho))
                   ((k) (if (odd*? k)
                            (/ (- k 1) 2)
                            (- (/ k 2) 1)))
                   ((rho) (make-flonum* (sqrt rho) k))
                   ((zeta) rho)
                   ((eta) y)
                   ((eta) (if (and (not (zero? rho)) (not (infinite? y)))
                              (/ eta rho 2.0)
                              eta)))
       (if (and (not (zero? rho)) (not (infinite? y))
                (negative? x))
           (make-rectangular (abs eta) (* rho (sign y)))
           (make-rectangular zeta eta)))))
  (else (define csqrt sqrt)))

;;; TODO: cacos, casin?

(define (sinh z)
  ;; The exponential formula does not keep the sign of zero.
  (let ((value (* -i (sin (* +i z)))))
    (if (strictly-real? z)
        (real-part value)
        value)))

(define (cosh z)
  (let ((value (* (cos (* +i z)))))
    (if (strictly-real? z)
        (real-part value)
        value)))

(cond-expand
  ((or (library (srfi 144))
       (library (srfi 276)))
   (define (casin z)
     (let ((x (real-part z))
           (s:1-z (csqrt (- 1 z)))
           (s:1+z (csqrt (+ 1 z))))
       (make-rectangular (atan x
                               (real-part (* s:1-z s:1+z)))
                         (flasinh (imag-part (* (conjugate s:1-z)
                                                s:1+z)))))))
  (else (define casin asin)))

(define (asinh z)
  (* -i (casin (* +i z))))

(define (%acosh z)
  (let* ((x (real-part z))
         (y (imag-part z))
         (sqrt:z-1 (csqrt (- z 1)))
         (sqrt:z+1 (csqrt (+ z 1))))
    (make-rectangular (flasinh (real-part (* (conjugate sqrt:z+1)
                                             sqrt:z-1)))
                      (* 2 (atan (imag-part sqrt:z-1)
                                 (real-part sqrt:z+1))))))

(define (acosh z)
  (if (eqv? z 1)
      0
      (%acosh z)))

;;; These are constants defined in Kahan's paper.
;;; They don't have better names, AFAIK.

(define theta (/ (sqrt fl-greatest) 4))
(define rho (/ theta))

(define (sign x)
  (cond
    ((and (zero? x) (not (eqv? x 0.0))) -1)
    ((negative? x) -1)
    (else 1)))

(define (%atanh z)
  (let* ((z (* (sign (real-part z)) (conjugate z)))
         (x (real-part z))
         (y (imag-part z)))
    (cond
      ((or (> x theta) (> (abs y) theta))
       (make-rectangular (real-part (/ z))
                         (* fl-pi/2 (sign y))))
      ((and (= x 1.0) (zero? y))
       (make-rectangular -inf.0
                         (* (sign y) .7853981633974483)))
      ((= x 1.0)
       (let ((ay+rho (+ (abs y) rho)))
         (make-rectangular (log (/ (sqrt (sqrt (+ 4.0 (square y))))
                                   (sqrt ay+rho)))
                           (* (/ (+ fl-pi/2
                                    (atan ay+rho 2.0))
                                 2.0)
                              (sign y)))))
      (else
       (let ((y+rho (square (+ (abs y) rho))))
         (make-rectangular (/ (fllog1+ (/ (* 4.0 x)
                                          (+ (square (- 1.0 x))
                                             y+rho)))
                              4.0)
                           (/ (angle (+ (* (- 1.0 x) (+ 1.0 x))
                                        (- y+rho)
                                        (* +2.0i y)))
                              2.0)))))))

(define (atanh z)
  (cond
    ((and (strictly-real? z) (eqv? (abs z) 1))
     (error 'atanh "atanh has a singularity at 1 and -1"))
    ((and (strictly-real? z) (< -1.0 z 1.0))
     (flatanh (inexact z)))
    (else (* (sign (real-part z)) (conjugate (%atanh z))))))

(define tanh-overflow-treshold
  (real-part (/ (asinh fl-greatest) 4)))

(define (tanh z)
  (let ((x (real-part z))
        (y (imag-part z)))
    (if (> (abs x) tanh-overflow-treshold)
        (if (strictly-real? z)
            (* 1.0 (sign x))
            (make-rectangular (* 1.0 (sign x))
                              (* 0.0 (sign y))))
        (let* ((t (tan y))
               (beta (+ 1.0 (square t)))
               (s (if (eqv? x 0)
                      0.0              ; Avoid divide-by-exact-zero errors
                      (sinh x)))
               (rho (sqrt (+ 1.0 (square s)))))
          (if (infinite? t)
              (make-rectangular (/ rho s) (/ t))
              (let ((ret (if (strictly-real? z)
                             (* beta rho s)
                             (make-rectangular (* beta rho s)
                                               t))))
                (/ ret (+ 1.0 (* beta (square s))))))))))

(define (round-away x)
  (if (negative? x)
      (truncate (- x 1/2))
      (truncate (+ x 1/2))))
