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

/-- Angular constant-term extraction as an additive homomorphism. -/
noncomputable def constantTermHom {R : Type*} [Semiring R] :
    LaurentPolynomial R →+ R :=
  (Finsupp.applyAddHom 0).comp AddMonoidAlgebra.coeffAddEquiv.toAddMonoidHom

/-- Angular constant-term extraction. -/
noncomputable def constantTerm {R : Type*} [Semiring R]
    (f : LaurentPolynomial R) : R :=
  constantTermHom f

@[simp] theorem constantTermHom_single
    {R : Type*} [Semiring R] (n : ℤ) (a : R) :
    constantTermHom (AddMonoidAlgebra.single n a) =
      if n = 0 then a else 0 := by
  change (Finsupp.single n a) 0 = if n = 0 then a else 0
  exact Finsupp.single_apply

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
  change
    constantTermHom
      (∑ k ∈ f.coeff.support,
        AddMonoidAlgebra.single k (f.coeff k) ^ p) =
      constantTermHom
        (∑ k ∈ f.coeff.support,
          AddMonoidAlgebra.single k (f.coeff k)) ^ p
  rw [map_sum, map_sum]
  simp only [AddMonoidAlgebra.single_pow, constantTermHom_single]
  simp [hp.ne_zero]

end GMC2
