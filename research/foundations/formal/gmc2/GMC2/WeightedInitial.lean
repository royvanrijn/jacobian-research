import Mathlib.Algebra.MonoidAlgebra.Basic
import Mathlib.Data.Finsupp.Basic
import Mathlib.Tactic

/-!
# Weighted initial forms of finitely supported bigraded expressions

This module packages the associated-graded calculation used by the lower
face argument.  The first coordinate is the integral angular weight and the
second is the natural radial degree.
-/

namespace GMC2

/-- A flattened circular polynomial: finitely supported in angular weight
and radial degree. -/
abbrev Bigraded (R : Type*) [Semiring R] :=
  AddMonoidAlgebra R (ℤ × ℕ)

/-- Rational weight attached to an angular/radial bidegree. -/
def bidegreeWeight (θ : ℚ) (x : ℤ × ℕ) : ℚ :=
  (x.2 : ℚ) - θ * (x.1 : ℚ)

@[simp] theorem bidegreeWeight_add (θ : ℚ) (x y : ℤ × ℕ) :
    bidegreeWeight θ (x + y) =
      bidegreeWeight θ x + bidegreeWeight θ y := by
  simp [bidegreeWeight]
  ring

/-- Every monomial of `f` has weight at least `ρ`. -/
def HasLowerBidegreeWeight {R : Type*} [Semiring R]
    (θ ρ : ℚ) (f : Bigraded R) : Prop :=
  ∀ x ∈ f.coeff.support, ρ ≤ bidegreeWeight θ x

/-- The exact weight-`ρ` part of a bigraded expression. -/
noncomputable def weightedInitial {R : Type*} [Semiring R]
    (θ ρ : ℚ) (f : Bigraded R) : Bigraded R :=
  ⟨f.coeff.filter fun x ↦ bidegreeWeight θ x = ρ⟩

@[simp] theorem weightedInitial_coeff {R : Type*} [Semiring R]
    (θ ρ : ℚ) (f : Bigraded R) (x : ℤ × ℕ) :
    (weightedInitial θ ρ f).coeff x =
      if bidegreeWeight θ x = ρ then f.coeff x else 0 := by
  rfl

theorem weightedInitial_mul {R : Type*} [CommSemiring R]
    (θ ρ σ : ℚ) (f g : Bigraded R)
    (hf : HasLowerBidegreeWeight θ ρ f)
    (hg : HasLowerBidegreeWeight θ σ g) :
    weightedInitial θ (ρ + σ) (f * g) =
      weightedInitial θ ρ f * weightedInitial θ σ g := by
  classical
  ext x
  simp only [weightedInitial_coeff, AddMonoidAlgebra.coeff_mul]
  simp only [weightedInitial, Finsupp.sum, Finsupp.support_filter,
    Finset.sum_filter, Finsupp.filter_apply]
  by_cases hx : bidegreeWeight θ x = ρ + σ
  · rw [if_pos hx]
    apply Finset.sum_congr rfl
    intro u hu
    by_cases huEq : bidegreeWeight θ u = ρ
    · simp only [huEq, if_true]
      apply Finset.sum_congr rfl
      intro v hv
      by_cases huv : u + v = x
      · have hsum :
            bidegreeWeight θ u + bidegreeWeight θ v = ρ + σ := by
          rw [← bidegreeWeight_add, huv, hx]
        have hvLower := hg v hv
        have hvEq : bidegreeWeight θ v = σ := by linarith
        simp [huv, hvEq]
      · simp [huv]
    · rw [if_neg huEq]
      apply Finset.sum_eq_zero
      intro v hv
      by_cases huv : u + v = x
      · have hsum :
            bidegreeWeight θ u + bidegreeWeight θ v = ρ + σ := by
          rw [← bidegreeWeight_add, huv, hx]
        have huLower := hf u hu
        have hvLower := hg v hv
        have huStrict : ρ < bidegreeWeight θ u :=
          lt_of_le_of_ne huLower (Ne.symm huEq)
        have : False := by linarith
        exact this.elim
      · simp [huv]
  · rw [if_neg hx]
    symm
    apply Finset.sum_eq_zero
    intro u hu
    by_cases huEq : bidegreeWeight θ u = ρ
    · simp only [huEq, if_true]
      apply Finset.sum_eq_zero
      intro v hv
      by_cases hvEq : bidegreeWeight θ v = σ
      · simp only [hvEq, if_true]
        by_cases huv : u + v = x
        · exfalso
          apply hx
          rw [← huv, bidegreeWeight_add, huEq, hvEq]
        · simp [huv]
      · simp [hvEq]
    · simp [huEq]

theorem hasLowerBidegreeWeight_mul {R : Type*} [CommSemiring R]
    (θ ρ σ : ℚ) (f g : Bigraded R)
    (hf : HasLowerBidegreeWeight θ ρ f)
    (hg : HasLowerBidegreeWeight θ σ g) :
    HasLowerBidegreeWeight θ (ρ + σ) (f * g) := by
  classical
  intro x hx
  have hx' := AddMonoidAlgebra.support_coeff_mul_subset f g hx
  rcases Finset.mem_add.mp hx' with ⟨u, hu, v, hv, huv⟩
  rw [← huv, bidegreeWeight_add]
  exact add_le_add (hf u hu) (hg v hv)

theorem hasLowerBidegreeWeight_pow {R : Type*} [CommSemiring R]
    (θ ρ : ℚ) (f : Bigraded R)
    (hf : HasLowerBidegreeWeight θ ρ f) (m : ℕ) :
    HasLowerBidegreeWeight θ ((m : ℚ) * ρ) (f ^ m) := by
  induction m with
  | zero =>
      intro x hx
      have hx0 : x = 0 := by
        by_contra hne
        exact (Finsupp.mem_support_iff.mp hx) (by
          simp [AddMonoidAlgebra.one_def, hne])
      subst x
      simp [bidegreeWeight]
  | succ m ih =>
      rw [pow_succ]
      convert hasLowerBidegreeWeight_mul θ ((m : ℚ) * ρ) ρ (f ^ m) f ih hf using 1
      push_cast
      ring

theorem weightedInitial_pow {R : Type*} [CommSemiring R]
    (θ ρ : ℚ) (f : Bigraded R)
    (hf : HasLowerBidegreeWeight θ ρ f) (m : ℕ) :
    weightedInitial θ ((m : ℚ) * ρ) (f ^ m) =
      weightedInitial θ ρ f ^ m := by
  induction m with
  | zero =>
      ext x
      by_cases hx : x = 0
      · subst x
        simp [weightedInitial_coeff, bidegreeWeight, AddMonoidAlgebra.one_def]
      · simp [weightedInitial_coeff, bidegreeWeight, AddMonoidAlgebra.one_def, hx]
  | succ m ih =>
      rw [pow_succ, pow_succ]
      have hm := hasLowerBidegreeWeight_pow θ ρ f hf m
      rw [show ((m + 1 : ℕ) : ℚ) * ρ = (m : ℚ) * ρ + ρ by
        push_cast
        ring]
      rw [weightedInitial_mul θ ((m : ℚ) * ρ) ρ (f ^ m) f hm hf, ih]

end GMC2
