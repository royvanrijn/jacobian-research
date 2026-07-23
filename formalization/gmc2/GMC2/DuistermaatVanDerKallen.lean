import Mathlib.Algebra.Polynomial.Laurent
import Mathlib.Algebra.CharP.Algebra
import Mathlib.RingTheory.Localization.FractionRing

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

/-- Domain-valued form obtained from the imported field theorem by passing
to the fraction field. -/
theorem duistermaat_van_der_kallen_domain
    {A : Type*} [CommRing A] [IsDomain A] [CharZero A]
    (f : LaurentPolynomial A) (hconv : SupportStraddlesZero f) :
    ∃ r : ℕ, 0 < r ∧ (f ^ r).coeff 0 ≠ 0 := by
  letI : CharZero (FractionRing A) := IsFractionRing.charZero A
  let ι : A →+* FractionRing A := algebraMap A (FractionRing A)
  let f' : LaurentPolynomial (FractionRing A) :=
    AddMonoidAlgebra.mapRingHom ℤ ι f
  have hconv' : SupportStraddlesZero f' := by
    rcases hconv with ⟨⟨a, ha, ha0⟩, ⟨b, hb, hb0⟩⟩
    refine ⟨⟨a, ?_, ha0⟩, ⟨b, ?_, hb0⟩⟩
    · apply Finsupp.mem_support_iff.mpr
      change (AddMonoidAlgebra.mapRingHom ℤ ι f).coeff a ≠ 0
      rw [AddMonoidAlgebra.coeff_mapRingHom]
      intro h
      apply Finsupp.mem_support_iff.mp ha
      apply IsFractionRing.injective A (FractionRing A)
      simpa using h
    · apply Finsupp.mem_support_iff.mpr
      change (AddMonoidAlgebra.mapRingHom ℤ ι f).coeff b ≠ 0
      rw [AddMonoidAlgebra.coeff_mapRingHom]
      intro h
      apply Finsupp.mem_support_iff.mp hb
      apply IsFractionRing.injective A (FractionRing A)
      simpa using h
  obtain ⟨r, hr, hcoeff⟩ := duistermaat_van_der_kallen f' hconv'
  refine ⟨r, hr, ?_⟩
  intro hzero
  apply hcoeff
  change ((AddMonoidAlgebra.mapRingHom ℤ ι f) ^ r).coeff 0 = 0
  rw [← map_pow, AddMonoidAlgebra.coeff_mapRingHom, hzero, map_zero]

end GMC2
