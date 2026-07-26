import GMC2.GaussianModel
import GMC2.PrimeIsolation
import Mathlib.Algebra.CharP.Lemmas
import Mathlib.Tactic

/-!
# Normalized factorial moments

The characteristic-zero moment is divided by the lowest factorial before
specialization.  This module makes that cancellation and the subsequent
Frobenius calculation explicit.
-/

namespace GMC2

open Polynomial

/-- The factorial functional normalized by `N!`; it is used only when all
degrees below `N` vanish. -/
noncomputable def normalizedFactorialMomentHom
    {R : Type*} [CommRing R] (N : ℕ) : R[X] →ₗ[R] R :=
  Polynomial.lsum fun j ↦
    LinearMap.mulLeft R (factorialQuotient N j : R)

noncomputable def normalizedFactorialMoment
    {R : Type*} [CommRing R] (N : ℕ) (f : R[X]) : R :=
  normalizedFactorialMomentHom N f

theorem normalizedFactorialMoment_apply
    {R : Type*} [CommRing R] (N : ℕ) (f : R[X]) :
    normalizedFactorialMoment N f =
      f.sum fun j a ↦ (factorialQuotient N j : R) * a := by
  rfl

@[simp] theorem normalizedFactorialMoment_add
    {R : Type*} [CommRing R] (N : ℕ) (f g : R[X]) :
    normalizedFactorialMoment N (f + g) =
      normalizedFactorialMoment N f + normalizedFactorialMoment N g :=
  (normalizedFactorialMomentHom N).map_add f g

@[simp] theorem normalizedFactorialMoment_sum
    {R : Type*} [CommRing R] {ι : Type*}
    (N : ℕ) (s : Finset ι) (f : ι → R[X]) :
    normalizedFactorialMoment N (∑ i ∈ s, f i) =
      ∑ i ∈ s, normalizedFactorialMoment N (f i) := by
  simp [normalizedFactorialMoment]

@[simp] theorem normalizedFactorialMoment_monomial
    {R : Type*} [CommRing R] (N j : ℕ) (a : R) :
    normalizedFactorialMoment N (Polynomial.monomial j a) =
      (factorialQuotient N j : R) * a := by
  simp [normalizedFactorialMoment, normalizedFactorialMomentHom]

theorem factorialFunctional_eq_factorial_mul_normalized
    {R : Type*} [CommRing R] (N : ℕ) (f : R[X])
    (hlow : ∀ j < N, f.coeff j = 0) :
    factorialFunctional f =
      (Nat.factorial N : R) * normalizedFactorialMoment N f := by
  classical
  rw [factorialFunctional, normalizedFactorialMoment_apply, Polynomial.sum_def,
    Polynomial.sum_def, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j hj
  have hjN : N ≤ j := by
    by_contra h
    exact (Polynomial.mem_support_iff.mp hj)
      (hlow j (Nat.lt_of_not_ge h))
  have hcast :
      (Nat.factorial N : R) * (factorialQuotient N j : R) =
        (Nat.factorial j : R) := by
    rw [← Nat.cast_mul]
    congr 1
    exact factorialQuotient_mul hjN
  rw [← hcast]
  ring

theorem normalizedFactorialMoment_eq_zero_of_factorialFunctional_eq_zero
    {R : Type*} [CommRing R] [IsDomain R] [CharZero R]
    (N : ℕ) (f : R[X])
    (hlow : ∀ j < N, f.coeff j = 0)
    (hmoment : factorialFunctional f = 0) :
    normalizedFactorialMoment N f = 0 := by
  have hfac : (Nat.factorial N : R) ≠ 0 := by
    exact_mod_cast Nat.factorial_ne_zero N
  have hmul :
      (Nat.factorial N : R) * normalizedFactorialMoment N f = 0 := by
    rw [← factorialFunctional_eq_factorial_mul_normalized N f hlow, hmoment]
  exact (mul_eq_zero.mp hmul).resolve_left hfac

theorem normalizedFactorialMoment_map
    {R S : Type*} [CommRing R] [CommRing S]
    (φ : R →+* S) (N : ℕ) (f : R[X]) :
    φ (normalizedFactorialMoment N f) =
      normalizedFactorialMoment N (f.map φ) := by
  induction f using Polynomial.induction_on' with
  | add f g hf hg =>
      simp [hf, hg]
  | monomial j a =>
      simp

/-- In prime characteristic, the normalized moment of a Frobenius power
is exactly the finite sum used by `primeDilatedMoment`. -/
theorem normalizedFactorialMoment_frobenius
    {K : Type*} [CommRing K] {p n : ℕ} [CharP K p]
    (hp : p.Prime) (f : K[X]) :
    normalizedFactorialMoment (n * p) (f ^ p) =
      primeDilatedMoment p n f := by
  classical
  letI : Fact p.Prime := ⟨hp⟩
  have hrepr :
      f = ∑ j ∈ f.support, Polynomial.monomial j (f.coeff j) := by
    rw [← Polynomial.sum_def]
    exact (Polynomial.sum_monomial_eq f).symm
  conv_lhs => rw [hrepr, sum_pow_char]
  rw [normalizedFactorialMoment_sum, primeDilatedMoment]
  apply Finset.sum_congr rfl
  intro j hj
  simp [normalizedFactorialMoment_monomial]

end GMC2
