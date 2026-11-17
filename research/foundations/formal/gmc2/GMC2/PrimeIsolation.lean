import GMC2.FactorialQuotient
import Mathlib.Algebra.CharP.Basic
import Mathlib.Algebra.Polynomial.Basic

/-!
# Prime-dilated factorial isolation

This is the arithmetic heart of the lower-face proof, stated independently
of Laurent polynomials.  It says that after Frobenius dilation, division by
the lowest factorial kills every higher radial degree modulo `p`.
-/

namespace GMC2

open Polynomial

/-- The normalized factorial sum attached to the Frobenius dilation of `f`.
The term indexed by `j` represents radial degree `jp`. -/
noncomputable def primeDilatedMoment
    {R : Type*} [CommRing R] (p n : ℕ) (f : R[X]) : R :=
  ∑ j ∈ f.support,
    (factorialQuotient (n * p) (j * p) : R) * f.coeff j ^ p

theorem cast_factorialQuotient_mul_eq_zero
    {R : Type*} [CommRing R] {p n j : ℕ} [CharP R p]
    (hp : p.Prime) (hnj : n < j) :
    (factorialQuotient (n * p) (j * p) : R) = 0 := by
  rw [CharP.cast_eq_zero_iff R p]
  exact prime_dvd_factorialQuotient_mul hp.pos hnj

/-- Factorial isolation: only the lowest nonzero coefficient survives. -/
theorem primeDilatedMoment_eq_lowestCoeff_pow
    {R : Type*} [CommRing R] {p n : ℕ} [CharP R p]
    (hp : p.Prime) (f : R[X])
    (hlow : ∀ j < n, f.coeff j = 0) (hc : f.coeff n ≠ 0) :
    primeDilatedMoment p n f = f.coeff n ^ p := by
  classical
  have hnmem : n ∈ f.support := Polynomial.mem_support_iff.mpr hc
  rw [primeDilatedMoment, Finset.sum_eq_single n]
  · rw [factorialQuotient, Nat.div_self (Nat.factorial_pos _)]
    simp
  · intro j hj hne
    rcases lt_or_gt_of_ne hne with hjn | hnj
    · simp [hlow j hjn, hp.ne_zero]
    · simp [cast_factorialQuotient_mul_eq_zero hp hnj]
  · exact fun h ↦ (h hnmem).elim

/-- The contradiction form used by the GMC argument. -/
theorem lowestCoeff_eq_zero_of_primeDilatedMoment_eq_zero
    {K : Type*} [Field K] {p n : ℕ} [CharP K p]
    (hp : p.Prime) (f : K[X])
    (hlow : ∀ j < n, f.coeff j = 0)
    (hmoment : primeDilatedMoment p n f = 0) :
    f.coeff n = 0 := by
  by_contra hc
  have hpow : f.coeff n ^ p = 0 := by
    rw [← primeDilatedMoment_eq_lowestCoeff_pow hp f hlow hc, hmoment]
  exact (pow_ne_zero p hc) hpow

end GMC2
