;;; SPDX-FileCopyrightText: 2026 Peter McGoron
;;; SPDX-License-Identifier: MIT

(define signed-imaginary-zero?
  ;; True only when the sign of *imaginary* zero is distinguished.
  ;; Gauche, CHICKEN, for example, don't, even when the distinguish
  ;; the sign of real zero.
  (not (eqv? (make-rectangular 0.0 0.0)
             (make-rectangular 0.0 -0.0))))

(define needs-strict-definition?
  (cond
    ((not (real? 0.0+0.0i)) #f)
    (signed-imaginary-zero? #t)
    ((exact? (imag-part 0.0)) #t)
    (else #f)))

(define (nan? obj)
  (and (number? obj) (r7rs:nan? obj)))

(define (exact-integer? obj)
  (and (integer? obj) (exact? obj)))

(define strictly-real?
  (if needs-strict-definition?
      (lambda (obj)
        (and (complex? obj)
             (zero? (imag-part obj))
             (exact? (imag-part obj))))
      real?))

(define strictly-rational?
  (if needs-strict-definition?
      (lambda (obj)
        (and (rational? obj)
             (exact? (imag-part obj))))
      rational?))

(define strictly-integer?
  (if needs-strict-definition?
      (lambda (obj)
        (and (integer? obj)
             (exact? (imag-part obj))))
      integer?))

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
           (let ((underflowed? (or (< x^2 fl-least-normal)
                                   (< y^2 fl-least-normal)))
                 (overflowed? (or (infinite? rho)
                                  (infinite? x^2)
                                  (infinite? y^2))))
             (if (or overflowed?
                     (and underflowed? (< rho (/ fl-least-normal
                                                 fl-epsilon))))
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
                   ((eta) (if (and (not (zero? rho)) (not (infinite? eta)))
                              (/ eta rho 2.0)
                              eta)))
       (if (and (not (zero? rho)) (negative? x))
           (make-rectangular (abs eta) (* rho (sign y)))
           (make-rectangular zeta eta)))))
  (else (define csqrt sqrt)))

(define (sinh z)
  (if (strictly-real? z)
      (flsinh (flonum z))
      (make-rectangular (* (flsinh (real-part z))
                           (cos (imag-part z)))
                        (* (flcosh (real-part z))
                           (sin (imag-part z))))))

(define (cosh z)
  (if (strictly-real? z)
      (flcosh (flonum z))
      (let ((x (flonum (real-part z)))
            (y (imag-part z)))
        (make-rectangular (* (flcosh x) (cos y))
                          (* (flsinh x) (sin y))))))

(cond-expand
  ((or (library (srfi 144))
       (library (srfi 276)))
   (define (casin z)
     (let ((x (real-part z))
           (s:1-z (csqrt (- 1 z)))
           (s:1+z (csqrt (+ 1 z))))
       (make-rectangular (atan x (real-part (* s:1-z s:1+z)))
                         (flasinh (imag-part (* (conjugate s:1-z)
                                                s:1+z)))))))
  (else (define casin asin)))

;;; Kahan's algorithm can cope with all unsigned zeros, or all signed
;;; zeroes, but some Scheme implementations have signed real zeroes
;;; but unsigned imaginary zeroes. This implementation of asinh causes
;;; the *clockwise* direction to be chosen for any signed zero input,
;;; because of these multiplications by +i storing the signed zero in
;;; the imaginary part, removing the sign.
;;;
;;; Not only is that the opposite convention of the unsigned zero case,
;;; it also erases sign information from the real part! This hack fixes
;;; this.

(cond-expand
  (chicken-6
   (define (*-i z)
     (make-rectangular (imag-part z)
                       (- (real-part z))))
   (define (*+i z)
     (make-rectangular (- (imag-part z))
                       (real-part z))))
  (else (define (*-i z) (* -i z))
        (define (*+i z) (* +i z))))

(define asinh
  (if signed-imaginary-zero?
      (lambda (z)
        (if (strictly-real? z)
            (flasinh (flonum z))
            (*-i (casin (*+i z)))))
      (lambda (z)
        (let ((w (* -i (casin (* +i z))))
              (x (real-part z))
              (y (imag-part z)))
          (cond
            ((positive? y)                        ; First or second
             (make-rectangular (* (sign x)        ; quadrant. The CCW rule
                                  (real-part w))  ; was applied, meaning
                               (imag-part w)))    ; that we should flip
                                                  ; the sign.

            ((zero? y) w)
            (else                                 ; Third or fourth
             (make-rectangular (* (sign x)        ; quadrant. We might
                                  -1              ; need to flip the sign,
                                  (real-part w))  ; but in the opposite
                               (imag-part w)))    ; scenarios.
          )))))

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
    ((and (exact? x) (zero? x)) 1)
    ((and (zero? x) (not (eqv? x 0.0))) -1)
    ((negative? x) -1)
    (else 1)))

(cond-expand
  ;; Work around a bug(?) in CHICKEN 6 where
  ;;   (* 1 +inf.0+1.0i) => +inf.0+nan.0i
  (chicken-6
   (define (c* x z)
     (make-rectangular (* x (real-part z))
                       (* x (imag-part z)))))
  (else (define c* *)))

(define (%atanh z)
  (let* ((z (c* (sign (real-part z)) (conjugate z)))
         (x (real-part z))
         (y (imag-part z)))
    (cond
      ((or (> x theta) (> (abs y) theta))
       (make-rectangular (real-part (/ z))
                         (* fl-pi/2 (sign y))))
      ((and (= x 1.0) (zero? y))
       (make-rectangular -inf.0
                         (* (sign y) fl-pi/4)))
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
                                        (make-rectangular
                                         0.0
                                         (* 2.0 y))))
                              2.0)))))))

(define (atanh z)
  (cond
    ((and (strictly-real? z) (eqv? (abs z) 1))
     (error 'atanh "atanh has a singularity at 1 and -1"))
    ((and (strictly-real? z) (< -1.0 z 1.0))
     (flatanh (flonum z)))
    (else (c* (sign (real-part z)) (conjugate (%atanh z))))))

(define tanh-overflow-treshold
  (real-part (/ (flasinh fl-greatest) 4)))

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
