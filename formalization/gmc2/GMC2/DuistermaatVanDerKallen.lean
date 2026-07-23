import Mathlib.Algebra.Polynomial.Laurent

/-!
# Duistermaat--van der Kallen

This is the sole imported mathematical theorem in the formalization.  The
one-variable statement below is the specialization used by GMC(2).  A future
mathlib import of the multivariate theorem can replace this declaration
without changing downstream files.
-/

namespace GMC2

open LaurentPolynomial

/-- In one dimension, `0 ∈ conv(support f)` is equivalent to the support
containing an exponent on each weak side of zero. -/
def SupportStraddlesZero {R : Type*} [Semiring R]
    (f : LaurentPolynomial R) : Prop :=
  (∃ i ∈ f.coeff.support, i ≤ 0) ∧ (∃ j ∈ f.coeff.support, 0 ≤ j)

/-- Imported Duistermaat--van der Kallen constant-term theorem, in exactly
the form consumed by the lower-face proof.

Reference: J. J. Duistermaat and W. van der Kallen,
*Constant terms in powers of a Laurent polynomial* (1998).
-/
axiom duistermaat_van_der_kallen
    {K : Type*} [Field K] [CharZero K] (f : LaurentPolynomial K)
    (hconv : SupportStraddlesZero f) :
    ∃ r : ℕ, 0 < r ∧ (f ^ r).coeff 0 ≠ 0

end GMC2
