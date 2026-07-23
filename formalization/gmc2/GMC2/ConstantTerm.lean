import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Algebra.CharP.Algebra
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
    (hp : p.Prime) (f : LaurentPolynomial R) :
    constantTerm (f ^ p) = constantTerm f ^ p := by
  classical
  letI : Fact p.Prime := ⟨hp⟩
  letI : CharP (LaurentPolynomial R) p :=
    charP_of_injective_ringHom (p := p) (f := LaurentPolynomial.C) (by
      intro x y h
      simpa using congrArg (fun g : LaurentPolynomial R ↦ g.coeff 0) h)
  rw [show f = ∑ k ∈ f.coeff.support, AddMonoidAlgebra.single k (f.coeff k) by
    simpa [Finsupp.sum] using (AddMonoidAlgebra.sum_coeff_single f).symm]
  rw [sum_pow_char]
  simp only [constantTerm, AddMonoidAlgebra.single_pow, AddMonoidAlgebra.coeff_sum,
    Finset.sum_apply, Finsupp.single_apply]
  change
    (∑ i ∈ f.coeff.support, if p • i = 0 then f.coeff i ^ p else 0) =
      (∑ i ∈ f.coeff.support, if i = 0 then f.coeff i else 0) ^ p
  simp [nsmul_eq_mul, hp.ne_zero]

end GMC2
