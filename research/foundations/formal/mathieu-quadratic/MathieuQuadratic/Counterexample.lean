import MathieuQuadratic.BaseChange
import MathieuQuadratic.Nonpower
import MathieuQuadratic.Rotations
import MathieuQuadratic.Scalar

/-!
# The uniform counterexample on the actual norm localization

The trace inputs are normalization, invariance under the scalar two, and
annihilation of the infinitesimal rotations. They isolate the geometric
Reynolds-functional input to the uniform construction.
-/

open scoped BigOperators TensorProduct

namespace MathieuQuadratic

variable {K : Type*} [Field K] [CharZero K] {M : ℕ}
variable (a : Fin 2 → K) (b : Fin M → K)

/-- The normalized trace properties used by the uniform construction. -/
structure ComplementTrace (trace : NormComplement a b →ₗ[K] K) : Prop where
  normalized : trace 1 = 1
  scalar : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f
  rotations : AnnihilatesRotations a b trace

noncomputable def complementPhase : NormComplement a b :=
  normalizerPhase (normUnit a b) (height a b)

omit [CharZero K] in
theorem complementPhase_eq :
    complementPhase a b =
      ((normUnit a b : NormComplement a b) + 1) *
        ((normUnit a b : NormComplement a b) -
          ((normUnit a b : NormComplement a b) + 1) ^ 2 *
            algebraMap (QuadPolynomial K M) (NormComplement a b) (heightNorm b)) ^ 2 *
        (↑((normUnit a b)⁻¹) : NormComplement a b) ^ 3 :=
  normalizerPhase_eq_regular _ _

theorem ComplementTrace.normalizer {trace : NormComplement a b →ₗ[K] K}
    (h : ComplementTrace a b trace) (hM : 0 < M) :
    NormalizerTrace M trace (normUnit a b) (height a b) :=
  ⟨trace_weight_of_scalar a b trace h.scalar,
    heightMoment_of_recurrence trace _ M hM h.normalized
      (height_recurrence_of_rotations a b trace h.rotations)⟩

variable {L : Type*} [Field L] [Algebra K L] [CharZero L] [IsAlgClosed L]

/-- The uniform phase and fixed marker on the diagonal quadratic
complement, with exact mixed moments and absolute unit-power exclusion.
The trace must satisfy the three stated Reynolds identities. -/
theorem uniform_counterexample {trace : NormComplement a b →ₗ[K] K}
    (h : ComplementTrace a b trace) (hM : 0 < M) (ha : a 0 ≠ 0) :
    let P := complementPhase a b
    let F := uniformPhase M P
    let H := (normUnit a b : NormComplement a b) ^ M
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L) F) ∧
      AllOrderWitness trace F H ∧
      ∀ m : ℕ, 0 < m → trace (H * F ^ m) =
        algebraMap ℚ K (∑ r ∈ Finset.range (m + 1), (m.choose r : ℚ) * betaMoment M (M * m + r)) :=
  ⟨uniformPhase_absolute_nonpower a b ha,
    uniformPhase_allOrderWitness (h.normalizer a b hM) hM⟩

theorem complement_not_mathieu {trace : NormComplement a b →ₗ[K] K}
    (h : ComplementTrace a b trace) (hM : 0 < M) : ¬ IsMathieu trace :=
  not_isMathieu_of_allOrderWitness ((uniformPhase_allOrderWitness (h.normalizer a b hM) hM).1)

/-- The ternary specialization uses P itself with the single marker q. -/
theorem direct_ternary_counterexample (b : Fin 1 → K) {trace : NormComplement a b →ₗ[K] K}
    (h : ComplementTrace a b trace) (ha : a 0 ≠ 0) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L) (complementPhase a b)) ∧
      AllOrderWitness trace (complementPhase a b) (normUnit a b : NormComplement a b) ∧
      ∀ m : ℕ, 0 < m → trace ((normUnit a b : NormComplement a b) * complementPhase a b ^ m) =
        algebraMap ℚ K (betaMoment 1 m) :=
  ⟨directPhase_absolute_nonpower a b ha, normalizer_direct_ternary (h.normalizer a b Nat.one_pos)⟩

end MathieuQuadratic
