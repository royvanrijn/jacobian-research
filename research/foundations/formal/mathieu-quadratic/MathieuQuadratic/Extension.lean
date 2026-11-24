import MathieuQuadratic.Counterexample
import MathieuQuadratic.LocalizationBaseChange

/-!
# The same normalizer phase and marker after extension
-/

open scoped TensorProduct BigOperators
open Algebra.TensorProduct

namespace MathieuQuadratic

variable {K L : Type*} [Field K] [Field L] [Algebra K L] {M : ℕ}
variable (a : Fin 2 → K) (b : Fin M → K)

private theorem map_heightNorm :
    MvPolynomial.map (algebraMap K L) (heightNorm b) =
      heightNorm (algebraMap K L ∘ b) := by
  simp [heightNorm]

theorem map_quadNorm :
    MvPolynomial.map (algebraMap K L) (quadNorm a b) =
      quadNorm (algebraMap K L ∘ a) (algebraMap K L ∘ b) := by
  rw [quadNorm, map_add, map_heightNorm, quadNorm]
  simp [planeNorm, Function.comp_def]

local notation "A" => NormComplement a b
local notation "aL" => algebraMap K L ∘ a
local notation "bL" => algebraMap K L ∘ b
local notation "B" => NormComplement aL bL

private noncomputable def identifyExtendedNorm :
    Localization.Away (MvPolynomial.map (algebraMap K L) (quadNorm a b)) ≃ₐ[L] B :=
  IsLocalization.algEquivOfAlgEquiv
    (M := Submonoid.powers (MvPolynomial.map (algebraMap K L) (quadNorm a b)))
    (T := Submonoid.powers (quadNorm aL bL))
    (Localization.Away (MvPolynomial.map (algebraMap K L) (quadNorm a b))) B
    (AlgEquiv.refl : QuadPolynomial L M ≃ₐ[L] QuadPolynomial L M)
    (by
      rw [Submonoid.map_powers]
      change Submonoid.powers (MvPolynomial.map (algebraMap K L) (quadNorm a b)) =
        Submonoid.powers (quadNorm aL bL)
      rw [map_quadNorm])

/-- The norm localization over the larger field is the base change of the
original norm localization. -/
noncomputable def complementBaseChange : L ⊗[K] A ≃ₐ[L] B :=
  (localizationBaseChange (L := L) (quadNorm a b)).trans (identifyExtendedNorm a b)

private theorem complementBaseChange_polynomial (f : QuadPolynomial K M) :
    complementBaseChange (L := L) a b (includeRight (algebraMap (QuadPolynomial K M) A f)) =
      algebraMap (QuadPolynomial L M) B (MvPolynomial.map (algebraMap K L) f) := by
  rw [complementBaseChange, AlgEquiv.trans_apply, localizationBaseChange_includeRight,
    identifyExtendedNorm, IsLocalization.algEquivOfAlgEquiv_eq]
  rfl

private theorem complementBaseChange_normUnit :
    Units.map ((complementBaseChange (L := L) a b).toMonoidHom.comp includeRight.toMonoidHom)
      (normUnit a b) = normUnit aL bL := by
  apply Units.ext
  change complementBaseChange (L := L) a b (includeRight (normUnit a b : A)) = (normUnit aL bL : B)
  rw [normUnit_val, complementBaseChange_polynomial, map_quadNorm, normUnit_val]

private theorem complementBaseChange_height :
    complementBaseChange (L := L) a b (includeRight (height a b)) = height aL bL := by
  simp only [height, map_mul]
  rw [complementBaseChange_polynomial, map_heightNorm]
  have h := congrArg (fun u : Bˣ ↦ ((↑(u⁻¹) : B))) (complementBaseChange_normUnit (L := L) a b)
  change _ * _ = _ * _
  rw [← h, ← map_inv]
  rfl

theorem complementBaseChange_phase :
    complementBaseChange (L := L) a b (includeRight (complementPhase a b)) =
      complementPhase aL bL := by
  have hn := congrArg (fun u : Bˣ ↦ (u : B)) (complementBaseChange_normUnit (L := L) a b)
  have hi := congrArg (fun u : Bˣ ↦ (↑(u⁻¹) : B)) (complementBaseChange_normUnit (L := L) a b)
  rw [← map_inv] at hi
  simp only [complementPhase, normalizerPhase, map_mul, map_add, map_one, map_pow, map_sub,
    complementBaseChange_height]
  rw [show complementBaseChange (L := L) a b (includeRight (↑((normUnit a b)⁻¹) : A)) =
    (↑((normUnit aL bL)⁻¹) : B) from hi,
    show complementBaseChange (L := L) a b (includeRight (normUnit a b : A)) =
    (normUnit aL bL : B) from hn]

theorem complementBaseChange_uniformPhase :
    complementBaseChange (L := L) a b (includeRight (uniformPhase M (complementPhase a b))) =
      uniformPhase M (complementPhase aL bL) := by
  simp only [uniformPhase, map_mul, map_add, map_one, map_pow, complementBaseChange_phase]

theorem complementBaseChange_marker :
    complementBaseChange (L := L) a b (includeRight ((normUnit a b : A) ^ M)) =
      (normUnit aL bL : B) ^ M := by
  simp only [map_pow, normUnit_val, complementBaseChange_polynomial, map_quadNorm]

end MathieuQuadratic
