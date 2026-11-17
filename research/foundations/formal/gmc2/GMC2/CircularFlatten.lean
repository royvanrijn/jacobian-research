import GMC2.GaussianModel
import GMC2.WeightedInitial
import Mathlib.Algebra.MonoidAlgebra.MapDomain

/-!
# Flattening circular polynomials

The Laurent-polynomial-of-polynomials representation is canonically the
monoid algebra on angular/radial bidegrees.  This file records that
equivalence and the projection which forgets radial degree.
-/

namespace GMC2

open Polynomial LaurentPolynomial

/-- Curry/uncurry identifies a circular polynomial with a finitely
supported expression on `ℤ × ℕ`. -/
noncomputable def circularBigradedEquiv
    (R : Type*) [CommRing R] :
    CircularPolynomial R ≃+* Bigraded R :=
  (AddMonoidAlgebra.mapRingEquiv ℤ (Polynomial.toFinsuppIso R)).trans
    AddMonoidAlgebra.curryRingEquiv.symm

@[simp] theorem circularBigradedEquiv_coeff
    {R : Type*} [CommRing R] (P : CircularPolynomial R) (k : ℤ) (j : ℕ) :
    (circularBigradedEquiv R P).coeff (k, j) = (P.coeff k).coeff j := by
  rfl

/-- Projection from a bidegree to its angular coordinate. -/
def angularProjection : (ℤ × ℕ) →+ ℤ where
  toFun x := x.1
  map_zero' := rfl
  map_add' _ _ := rfl

/-- Forget radial degree, summing coefficients with the same angular
weight. -/
noncomputable def forgetRadial
    (R : Type*) [CommRing R] : Bigraded R →+* LaurentPolynomial R :=
  AddMonoidAlgebra.mapDomainRingHom R angularProjection

theorem forgetRadial_coeff
    {R : Type*} [CommRing R] (f : Bigraded R) (k : ℤ) :
    (forgetRadial R f).coeff k =
      ∑ x ∈ f.coeff.support with x.1 = k, f.coeff x := by
  classical
  rw [forgetRadial, AddMonoidAlgebra.mapDomainRingHom_apply,
    AddMonoidAlgebra.coeff_mapDomain]
  simpa [angularProjection] using
    (Finsupp.mapDomain_apply_eq_sum (fun x : ℤ × ℕ ↦ x.1) f.coeff
      (a := (k, 0)))

end GMC2
