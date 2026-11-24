import Mathlib.Algebra.Algebra.Hom

/-!
# Mathieu kernels of linear functionals

The definitions use positive powers and one cutoff per fixed marker, as in
the quadratic-complement dimension classification.
-/

namespace MathieuQuadratic

variable {K A B : Type*} [Field K] [CommRing A] [Algebra K A]

/-- Every positive power of the phase is in the kernel. -/
def Pure (trace : A →ₗ[K] K) (phase : A) : Prop :=
  ∀ n : ℕ, 0 < n → trace (phase ^ n) = 0

/-- Each fixed multiplier has its own eventual vanishing threshold. -/
def Safe (trace : A →ₗ[K] K) (phase : A) : Prop :=
  ∀ marker : A, ∃ cutoff : ℕ, ∀ n ≥ cutoff, trace (marker * phase ^ n) = 0

/-- The commutative Mathieu property for the kernel of a functional. -/
def IsMathieu (trace : A →ₗ[K] K) : Prop :=
  ∀ phase : A, Pure trace phase → Safe trace phase

/-- One marker has nonzero moment at every positive order. -/
def AllOrderWitness (trace : A →ₗ[K] K) (phase marker : A) : Prop :=
  Pure trace phase ∧ ∀ n : ℕ, 0 < n → trace (marker * phase ^ n) ≠ 0

theorem not_isMathieu_of_allOrderWitness {trace : A →ₗ[K] K} {phase marker : A}
    (h : AllOrderWitness trace phase marker) : ¬ IsMathieu trace := by
  intro mathieu
  obtain ⟨cutoff, hcutoff⟩ := mathieu phase h.1 marker
  exact h.2 (cutoff + 1) (Nat.zero_lt_succ _) (hcutoff _ (Nat.le_succ _))

/-- Unit-times-power exclusion in a specified scalar extension. -/
def NonpowerUpToUnits [CommRing B] (phase : B) : Prop :=
  ∀ (u : Bˣ) (g : B) (e : ℕ), 2 ≤ e → phase ≠ ↑u * g ^ e

/-- The deduction of the sharp dimension criterion from the two geometric
inputs in Sections 3–5 and 8 of the classification. For the absolute
criterion, the specified extension is to the algebraic closure. -/
theorem sharp_dimension_of_inputs [CommRing B] [Algebra K B] (numVars : ℕ)
    (trace : A →ₗ[K] K) (intoExtension : A →ₐ[K] B)
    (low : numVars ≤ 2 → IsMathieu trace)
    (high : 3 ≤ numVars → ∃ phase marker,
      NonpowerUpToUnits (intoExtension phase) ∧ AllOrderWitness trace phase marker) :
    (IsMathieu trace ↔ numVars ≤ 2) ∧
      ((¬ ∃ phase marker, NonpowerUpToUnits (intoExtension phase) ∧
          AllOrderWitness trace phase marker) ↔ numVars ≤ 2) := by
  constructor
  · refine ⟨fun h ↦ ?_, low⟩
    by_contra hn
    obtain ⟨phase, marker, _, witness⟩ := high (Nat.succ_le_of_lt (Nat.lt_of_not_ge hn))
    exact not_isMathieu_of_allOrderWitness witness h
  · constructor
    · intro absent
      by_contra hn
      exact absent (high (Nat.succ_le_of_lt (Nat.lt_of_not_ge hn)))
    · rintro hn ⟨phase, marker, _, witness⟩
      exact not_isMathieu_of_allOrderWitness witness (low hn)

end MathieuQuadratic
