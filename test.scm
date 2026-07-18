;;; SPDX-FileCopyrightText: 2026 Peter McGoron
;;; SPDX-License-Identifier: MIT

(import (except (scheme base) exact-integer?)
        (scheme complex)
        (rename (scheme inexact) (nan? r7rs:nan?))
        (srfi 64)
        (srfi 278))

(test-begin "SRFI 278")

(test-group "exact-integer?"
  (test-assert (exact-integer? 1))
  (test-assert (exact-integer? 0))
  (test-assert (not (exact-integer? 0.0)))
  (test-assert (not (exact-integer? "0.0"))))

(test-group "nan?"
  (test-assert (nan? +nan.0))
  (test-assert (not (nan? +inf.0)))
  (test-assert (not (nan? "NaN"))))

(test-group "round-away"
  (test-approximate 4.0
                    (round-away 3.5)
                    1e-6)
  (test-approximate 3.0
                    (round-away 2.5)
                    1e-6)
  (test-approximate 3.0
                    (round-away 2.6)
                    1e-6)
  (test-approximate 2.0
                    (round-away 2.4)
                    1e-6)
  (test-eqv 3 (round-away 5/2))
  (test-approximate -4.0
                    (round-away -3.5)
                    1e-6)
  (test-approximate -3.0
                    (round-away -2.5)
                    1e-6))

(test-group "conjugate"
  (test-eqv 1 (conjugate 1))
  (test-eqv -1 (conjugate -1))
  (test-eqv 1-2i (conjugate 1+2i))
  (test-eqv -0.0 (imag-part (conjugate 1.0+0.0i)))
  (test-eqv +i (conjugate -i)))

(test-group "sinh"
  ;; TODO: Better way to test?
  (test-assert (zero? (sinh 0)))
  (test-eqv 0.0 (sinh 0.0))
  (test-eqv -0.0 (sinh -0.0))
  (test-eqv +inf.0 (sinh +inf.0))
  (test-eqv -inf.0 (sinh -inf.0))
  (test-approximate (/ (- (exp 1) (exp -1)) 2)
                    (sinh 1)
                    1e-6))

(test-group "cosh"
  (test-eqv 1.0 (cosh 0))
  (test-eqv +inf.0 (cosh +inf.0))
  (test-eqv +inf.0 (cosh -inf.0))
  (test-approximate (/ (+ (exp 1) (exp -1)) 2)
                    (cosh 1)
                    1e-6))

(test-group "tanh"
  (test-eqv "(tanh 0.0)" 0.0 (tanh 0.0))
  (test-eqv "(tanh -0.0)" -0.0 (tanh -0.0))
  (test-eqv "(tanh +inf.0)" 1.0 (tanh +inf.0))
  (test-eqv "(tanh -inf.0)" -1.0 (tanh -inf.0))
  (test-approximate (/ (sinh 10) (cosh 10))
                    (tanh 10)
                    1e-6))

(define (naive-atanh z)
  (/ (- (log (+ 1 z)) (log (- 1 z)))
     2))

(define (test-approximate/special expect actual error)
  (cond
    ((eqv? expect +0.0)
     (test-assert (or (eqv? actual 0.0)
                      (eqv? actual 0))))
    ((eqv? expect -0.0)
     (test-assert (eqv? actual -0.0)))
    ((infinite? expect)
     (test-eqv expect actual))
    (else (test-approximate expect actual error))))

(define (test-real-and-imag z-expect z-actual)
  (test-approximate/special (real-part z-expect)
                            (real-part z-actual)
                            1e-6)
  (test-approximate/special (imag-part z-expect)
                            (imag-part z-actual)
                            1e-6))

(test-group "atanh"
  (test-group "(atanh 1.0+0.0i)"
    (test-real-and-imag -inf.0+.7853981633974483i
                        (atanh 1.0+0.0i)))
  (test-group "(atanh 1.0-0.0i)"
    (test-real-and-imag -inf.0-.7853981633974483i
                        (atanh 1.0-0.0i)))
  (test-group "(atanh -1.0+0.0i)"
    (test-real-and-imag +inf.0+.7853981633974483i
                        (atanh -1.0+0.0i)))
  (test-group "(atanh -1.0-0.0i)"
    (test-real-and-imag +inf.0-.7853981633974483i
                        (atanh -1.0-0.0i)))
  (test-group "(atanh 2.0+0.0i)"
    (test-real-and-imag .5493061443340549+1.5707963267948966i
                        (atanh 2.0+0.0i)))
  (test-group "(atanh 2.0-0.0i)"
    (test-real-and-imag .5493061443340549-1.5707963267948966i
                        (atanh 2.0-0.0i)))
  (test-group "(atanh 0.0+1.0i)"
    (test-real-and-imag 0.0+0.7853981633974483i
                        (atanh 0.0+1.0i)))
  (test-group "(atanh 0.0-1.0i)"
    (test-real-and-imag 0.0-.7853981633974483i
                        (atanh 0.0-1.0i)))
  (test-approximate "(atanh 0.5)"
                    (naive-atanh 0.5)
                    (atanh 0.5)
                    1e-6)
  (test-approximate "(atanh 0.99)"
                    (naive-atanh 0.99)
                    (atanh 0.99)
                    1e-6)
  (test-assert "(atanh 0)" (zero? (atanh 0)))
  (test-group "(atanh 2.0)"
    (test-real-and-imag (naive-atanh 2.0)
                        (atanh 2.0)))
  (test-group "(atanh 5+10i)"
    (test-real-and-imag (naive-atanh 5+10i)
                        (atanh 5+10i))))

(test-group "acosh"
  (test-group "(acosh 1.0+0.0i)"
    (test-real-and-imag 0.0+0.0i
                        (acosh 1.0+0.0i)))
  (test-group "(acosh 1.0-0.0i)"
    (test-real-and-imag 0.0-0.0i
                        (acosh 1.0-0.0i)))
  (test-group "(acosh 0.0+0.0i)"
    (test-real-and-imag 0.+1.5707963267948966i
                        (acosh 0.0+0.0i)))
  (test-group "(acosh 0.0-0.0i)"
    (test-real-and-imag 0.-1.5707963267948966i
                        (acosh 0.0-0.0i)))
  (test-group "(acosh -1+0.0i)"
    (test-real-and-imag 0.0+3.141592653589793i
                        (acosh -1+0.0i)))
  (test-group "(acosh -1-0.0i)"
    (test-real-and-imag 0.0-3.141592653589793i
                        (acosh -1-0.0i)))
  (test-group "(acosh 3+4i)"
    (test-real-and-imag 2.305509031243477+.9368124611557198i
                        (acosh 3+4i))))

(define (naive-asinh z)
  (log (+ z (sqrt (+ 1 (square z))))))

(test-group "asinh"
  (test-group "(asinh 0.0+2.0i)"
    (test-real-and-imag 1.3169578969248166+1.5707963267948966i
                        (asinh 0.0+2.0i)))
  (test-group "(asinh -0.0+2.0i)"
    (test-real-and-imag -1.3169578969248166+1.5707963267948966i
                        (asinh -0.0+2.0i)))
  (test-group "(asinh 0.0-2.0i)"
    (test-real-and-imag 1.3169578969248166-1.5707963267948966i
                        (asinh 0.0-2.0i)))
  (test-group "(asinh -0.0-2.0i)"
    (test-real-and-imag -1.3169578969248166-1.5707963267948966i
                        (asinh -0.0-2.0i)))
  (test-group "(asinh +i)"
    (test-real-and-imag +1.5707963267948966i
                        (asinh +i)))
  (test-group "(asinh 5)"
    (test-approximate (naive-asinh 5)
                      (asinh 5)
                      1e-6))
  (test-group "(asinh 1+2i)"
    (test-real-and-imag (naive-asinh 1+2i)
                        (asinh 1+2i))))

(test-end "SRFI 278")
