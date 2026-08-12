;;; SPDX-FileCopyrightText: 2026 Peter McGoron
;;; SPDX-License-Identifier: MIT

(define-library (srfi 278)
  (import (except (scheme base) exact-integer? real? rational? integer?
                                rationalize)
          (prefix (only (scheme base) real? rational? integer?)
                  r5rs:)
          (scheme write)
          (rename (scheme inexact)
                  (nan? r7rs:nan?))
          (scheme complex)
          (scheme case-lambda))
  (export nan? exact-integer?
          imaginary? real? rational? integer?
          sinh cosh tanh asinh acosh atanh
          conjugate rationalize round-away)
  (cond-expand
    ((or chicken (library (srfi 144)))  ; TODO: fix CHICKEN here
     (import (only (srfi 144)
                   flonum
                   fl-greatest
                   fl-epsilon
                   fl-least
                   flnormalized?
                   fladjacent
                   fl-pi/2
                   fl-pi/4
                   make-flonum
                   flexponent
                   flasinh
                   flsinh
                   flcosh
                   flatanh
                   fllog1+))
     (begin
       (define fl-least-normal
         (do ((candidate fl-least (* 2.0 candidate)))
             ((flnormalized? candidate) candidate)))))
    ;; If you don't have SRFI 144, you have to define the following
    ;; here:
    ;;
    ;; flonum (which is probably just `inexact`)
    ;; fl-greatest
    ;; fl-least-normal (not fl-least; the smallest normal number)
    ;; fl-epsilon
    ;; fl-pi/2
    ;; fl-pi/4
    ;; make-flonum (aka ldexp)
    ;; flexponent (aka logb)
    ;; fllog1+
    ;; flasinh
    ;; flsinh
    ;; flcosh
    ;; flatanh
    ;;
    ;; If your inexact real type is a IEEE 754 format number, then you
    ;; can use the SRFI 144 sample implementation.
    ;;
    ;; The other imports are used to implement Kahan's complex number
    ;; versions of Scheme's built-in procedures, like sqrt, asin, etc.
    ;; If neither 276 or 144 are available, then the implementation will
    ;; fall back to the standard procedures. Note that if the
    ;; implementation's procedures are inaccurate, it will affect the
    ;; results of the procedures here.
    ;;
    ;; I would be very interested in any Schemes using non-standard
    ;; floating point formats.
)
  (include "278.scm"))
