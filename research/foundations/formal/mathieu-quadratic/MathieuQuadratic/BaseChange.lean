import MathieuQuadratic.Basic
import Mathlib.RingTheory.TensorProduct.Basic

/-!
# Scalar extension of a fixed phase

This is Section 2 of the sharp quadratic-complement dimension classification.
-/

open scoped TensorProduct
open TensorProduct

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CommRing A] [Algebra K A]
variable (L : Type*) [Field L] [Algebra K L]

/-- The extended functional with target identified with the extension field. -/
noncomputable def extendTrace (trace : A →ₗ[K] K) : L ⊗[K] A →ₗ[L] L :=
  (AlgebraTensorModule.rid K L L).toLinearMap.comp (trace.baseChange L)

@[simp] theorem extendTrace_tmul (trace : A →ₗ[K] K) (c : L) (a : A) :
    extendTrace L trace (c ⊗ₜ[K] a) = algebraMap K L (trace a) * c := by
  simp [extendTrace, Algebra.smul_def]

theorem pure_includeRight_iff (trace : A →ₗ[K] K) (phase : A) :
    Pure (extendTrace L trace)
      (Algebra.TensorProduct.includeRight (A := L) phase : L ⊗[K] A) ↔
      Pure trace phase := by
  let i : A →ₐ[K] L ⊗[K] A := Algebra.TensorProduct.includeRight
  change Pure (extendTrace L trace) (i phase) ↔ Pure trace phase
  have eq (n : ℕ) : extendTrace L trace (i phase ^ n) =
      algebraMap K L (trace (phase ^ n)) := by
    rw [← map_pow, show i (phase ^ n) = 1 ⊗ₜ[K] (phase ^ n) from rfl,
      extendTrace_tmul, mul_one]
  simp only [Pure, eq, map_eq_zero_iff _ (RingHom.injective (algebraMap K L))]

theorem safe_includeRight_iff (trace : A →ₗ[K] K) (phase : A) :
    Safe (extendTrace L trace)
      (Algebra.TensorProduct.includeRight (A := L) phase : L ⊗[K] A) ↔
      Safe trace phase := by
  let i : A →ₐ[K] L ⊗[K] A := Algebra.TensorProduct.includeRight
  change Safe (extendTrace L trace) (i phase) ↔ Safe trace phase
  constructor
  · intro extended marker
    obtain ⟨cutoff, hcutoff⟩ := extended (i marker)
    refine ⟨cutoff, fun n hn ↦ ?_⟩
    have h := hcutoff n hn
    rwa [← map_pow, ← map_mul, show i (marker * phase ^ n) = 1 ⊗ₜ[K] (marker * phase ^ n)
      from rfl, extendTrace_tmul, mul_one,
      map_eq_zero_iff _ (RingHom.injective (algebraMap K L))] at h
  · intro safe marker
    induction marker using TensorProduct.induction_on with
    | zero => exact ⟨0, by simp⟩
    | tmul c a =>
        obtain ⟨cutoff, hcutoff⟩ := safe a
        refine ⟨cutoff, fun n hn ↦ ?_⟩
        rw [← map_pow, show i (phase ^ n) = 1 ⊗ₜ[K] (phase ^ n) from rfl,
          Algebra.TensorProduct.tmul_mul_tmul, mul_one, extendTrace_tmul,
          hcutoff n hn, map_zero, zero_mul]
    | add a b ha hb =>
        obtain ⟨cutoffA, hcutoffA⟩ := ha
        obtain ⟨cutoffB, hcutoffB⟩ := hb
        refine ⟨max cutoffA cutoffB, fun n hn ↦ ?_⟩
        rw [add_mul, map_add, hcutoffA n ((le_max_left _ _).trans hn),
          hcutoffB n ((le_max_right _ _).trans hn), add_zero]

theorem allOrderWitness_includeRight_iff (trace : A →ₗ[K] K) (phase marker : A) :
    AllOrderWitness (extendTrace L trace)
      (Algebra.TensorProduct.includeRight (A := L) phase : L ⊗[K] A)
      (Algebra.TensorProduct.includeRight (A := L) marker) ↔
    AllOrderWitness trace phase marker := by
  let i : A →ₐ[K] L ⊗[K] A := Algebra.TensorProduct.includeRight
  change AllOrderWitness (extendTrace L trace) (i phase) (i marker) ↔ _
  have eq (k : ℕ) : extendTrace L trace (i marker * i phase ^ k) =
      algebraMap K L (trace (marker * phase ^ k)) := by
    rw [← map_pow, ← map_mul, show i (marker * phase ^ k) = 1 ⊗ₜ[K] (marker * phase ^ k)
      from rfl, extendTrace_tmul, mul_one]
  have hp : Pure (extendTrace L trace) (i phase) ↔ Pure trace phase :=
    pure_includeRight_iff L trace phase
  simp only [AllOrderWitness, hp, eq, ne_eq,
    map_eq_zero_iff _ (RingHom.injective (algebraMap K L))]

theorem unsafe_includeRight_iff (trace : A →ₗ[K] K) (phase : A) :
    (¬ Safe (extendTrace L trace)
      (Algebra.TensorProduct.includeRight (A := L) phase : L ⊗[K] A)) ↔
      ¬ Safe trace phase :=
  not_congr (safe_includeRight_iff L trace phase)

theorem isMathieu_of_extended (trace : A →ₗ[K] K)
    (h : IsMathieu (extendTrace L trace)) : IsMathieu trace := by
  intro phase hp
  exact (safe_includeRight_iff L trace phase).1
    (h _ ((pure_includeRight_iff L trace phase).2 hp))

end MathieuQuadratic
