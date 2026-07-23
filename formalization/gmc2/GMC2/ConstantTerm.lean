import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Algebra.Polynomial.Laurent

/-!
# Constant terms and Frobenius

This file proves that coefficient extraction at exponent zero commutes with
the characteristic-`p` Frobenius on Laurent polynomials.
-/

namespace GMC2

open LaurentPolynomial

/-- Angular constant-term extraction. -/
def constantTerm {R : Type*} [Semiring R] (f : LaurentPolynomial R) : R :=
  f.coeff 0

theorem constantTerm_frobenius
    {R : Type*} [CommRing R] {p : ℕ} [CharP R p]
    (hp : 0 < p) (f : LaurentPolynomial R) :
    constantTerm (f ^ p) = constantTerm f ^ p := by
  classical
  rw [show f = ∑ k ∈ f.coeff.support, AddMonoidAlgebra.single k (f.coeff k) by
    simpa [Finsupp.sum] using (AddMonoidAlgebra.sum_coeff_single f).symm]
  rw [Finset.sum_pow_char]
  simp only [constantTerm, AddMonoidAlgebra.single_pow,
    AddMonoidAlgebra.coeff_sum, Finsupp.single_apply]
  rw [Finset.sum_eq_single 0]
  · simp
  · intro b hb hb0
    simp [hb0, hp.ne']
  · simp

end GMC2
