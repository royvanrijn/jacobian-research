import MathieuQuadratic.Counterexample
import MathieuQuadratic.Newton

/-!
# The sharp dimension deduction

A low-dimensional model identifies the extended trace with Laurent
constant coefficient. A high-dimensional model identifies the trace
with a normalized functional on the diagonal norm localization.
-/

open scoped TensorProduct

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

/-- The Laurent necessity input in the two ranks needed for the classification. -/
def LowRankLaurentTheorem : Prop :=
  ∀ n : ℕ, n = 1 ∨ n = 2 →
    NewtonLaurentNecessity ℂ (Fin n)

omit [CharZero K] in
private theorem transported_witness {L B : Type*} [Field L] [Algebra K L] [CommRing B] [Algebra K B]
    (coords : A ≃ₐ[K] B) (trace : A →ₗ[K] K) (f g : B)
    (hn : NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L) f))
    (hw : AllOrderWitness (trace.comp coords.symm.toLinearMap) f g) :
    ∃ phase marker,
      NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L) phase) ∧
        AllOrderWitness trace phase marker := by
  refine ⟨coords.symm f, coords.symm g, ?_, ?_, ?_⟩
  · intro unit root exponent he eq
    let lift := Algebra.TensorProduct.congr (AlgEquiv.refl : L ≃ₐ[L] L) coords
    have hl : lift (Algebra.TensorProduct.includeRight (coords.symm f)) =
        Algebra.TensorProduct.includeRight f := by
      simp [lift, Algebra.TensorProduct.includeRight_apply]
    have hx := congrArg lift eq
    rw [hl, map_mul, map_pow] at hx
    exact hn (Units.map lift.toMonoidHom unit) (lift root) exponent he hx
  · intro n hn
    rw [← map_pow]
    exact hw.1 n hn
  · intro n hn
    rw [← map_pow, ← map_mul]
    exact hw.2 n hn

/-- Sharp classification with the trace comparison inputs and the Laurent
necessity theorem stated explicitly. -/
theorem sharp_classification_of_models {L : Type*} [Field L] [CharZero L] [Algebra K L]
    [IsAlgClosed L] (numVars : ℕ) (hpos : 0 < numVars) (trace : A →ₗ[K] K)
    (dvk : LowRankLaurentTheorem)
    (split : numVars ≤ 2 →
      ∃ coordinates : L ⊗[K] A ≃ₐ[L] AddMonoidAlgebra L (Fin numVars → ℤ),
        ∀ f, extendTrace L trace f = laurentCT (coordinates f))
    (diagonal : 3 ≤ numVars →
      ∃ (a : Fin 2 → K) (b : Fin (numVars - 2) → K) (e : A ≃ₐ[K] NormComplement a b),
        a 0 ≠ 0 ∧ ComplementTrace a b (trace.comp e.symm.toLinearMap)) :
    (IsMathieu trace ↔ numVars ≤ 2) ∧
      ((¬ ∃ phase marker, NonpowerUpToUnits
        (Algebra.TensorProduct.includeRight (R := K) (A := L) phase) ∧
        AllOrderWitness trace phase marker) ↔ numVars ≤ 2) := by
  apply sharp_dimension_of_inputs numVars trace (Algebra.TensorProduct.includeRight (R := K) (A := L))
  · intro h
    obtain ⟨coords, ht⟩ := split h
    exact isMathieu_of_split_laurent trace coords ht (laurentSeparation_of_complex (laurentSeparation_of_newton (dvk numVars (by omega))))
  · intro h
    obtain ⟨a, b, coords, ha, ht⟩ := diagonal h
    obtain ⟨hn, hw, _⟩ := uniform_counterexample (L := L) a b ht (by omega) ha
    exact transported_witness coords trace _ _ hn hw

end MathieuQuadratic
