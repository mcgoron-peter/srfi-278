;;; SPDX-FileCopyrightText: 2026 Peter McGoron
;;; SPDX-License-Identifier: MIT

(define-library (srfi 278)
  (import (except (scheme base) exact-integer?)
          (scheme write)
          (rename (scheme inexact)
                  (nan? r7rs:nan?))
          (scheme complex))
  (export nan? exact-integer?
          strictly-real? strictly-rational? strictly-integer?
          sinh cosh tanh asinh acosh atanh
          conjugate
          round-away)
  (cond-expand
    ((library (srfi 276))
     (import (rename (only (srfi 276)
                           :flonum
                           :greatest :least :epsilon :pi/2
                           :expoennt
                           :make-flonum
                           :asin
                           :atanh
                           :log1+)
                     (:flonum flonum)
                     (:greatest fl-greatest)
                     (:least fl-least)
                     (:epsilon fl-epsilon)
                     (:pi/2 fl-pi/2)
                     (:make-flonum make-flonum)
                     (:exponent flexponent)
                     (:asinh flasinh)
                     (:atanh flatanh)
                     (:log1+ fllog1+))))
    ((library (srfi 144))
     (import (only (srfi 144)
                   flonum
                   fl-greatest
                   fl-least
                   fl-epsilon
                   fl-pi/2
                   make-flonum
                   flexponent
                   flasinh
                   flatanh
                   fllog1+)))
    ;; If you don't have SRFI 144, you have to define the following
    ;; here:
    ;; 
    ;; fl-greatest
    ;; fl-pi/2
    ;; fllog1+
    ;; flasinh
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
