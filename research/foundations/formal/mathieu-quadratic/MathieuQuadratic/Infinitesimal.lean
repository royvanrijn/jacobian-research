import MathieuQuadratic.Rotations
import Mathlib.Algebra.DualNumber

/-!
# First-order rotations on the norm localization

The polynomial rotation has an actual pullback over the dual numbers.
Its degree-one part on every numerator with norm-power denominator is
the localized derivative used in the height recurrence.
-/

open scoped DualNumber
open TrivSqZeroExt

namespace MathieuQuadratic

variable {K : Type*} [Field K] {M : ℕ}
variable (a : Fin 2 → K) (b : Fin M → K)

local notation "R" => QuadPolynomial K M
local notation "A" => NormComplement a b

private noncomputable def firstOrderPolynomial (i : Fin 2) (j : Fin M) :
    R →ₐ[K] A[ε] where
  toFun f := ⟨algebraMap R A f, algebraMap R A (rotation a b i j f)⟩
  map_zero' := by apply TrivSqZeroExt.ext <;> simp
  map_one' := by
    apply TrivSqZeroExt.ext
    · exact map_one _
    · simp
  map_add' f g := by
    apply TrivSqZeroExt.ext
    · exact map_add _ _ _
    · change algebraMap R A (rotation a b i j (f + g)) = _
      simp only [map_add, snd_add]
      rfl
  map_mul' f g := by
    apply TrivSqZeroExt.ext
    · exact map_mul _ _ _
    · change algebraMap R A (rotation a b i j (f * g)) = _
      rw [Derivation.leibniz, DualNumber.snd_mul]
      simp only [smul_eq_mul, map_add, map_mul]
      change _ = algebraMap R A f * algebraMap R A (rotation a b i j g) +
        algebraMap R A (rotation a b i j f) * algebraMap R A g
      ring
  commutes' c := by
    apply TrivSqZeroExt.ext
    · exact (IsScalarTower.algebraMap_apply K R A c).symm
    · change algebraMap R A (rotation a b i j (algebraMap K R c)) = 0
      rw [Derivation.map_algebraMap, map_zero]

private theorem firstOrderPolynomial_norm (i : Fin 2) (j : Fin M) :
    firstOrderPolynomial a b i j (quadNorm a b) = inl (normUnit a b : A) := by
  change (⟨algebraMap R A (quadNorm a b), algebraMap R A (rotation a b i j (quadNorm a b))⟩ :
    A[ε]) = _
  rw [rotation_quadNorm, map_zero, normUnit_val]
  rfl

private theorem firstOrderPolynomial_unit (i : Fin 2) (j : Fin M) :
    IsUnit (firstOrderPolynomial a b i j (quadNorm a b)) := by
  rw [firstOrderPolynomial_norm]
  exact (normUnit a b).isUnit.map (TrivSqZeroExt.inlHom A A)

/-- The first-order pullback on the regular coordinate ring. -/
noncomputable def firstOrderRotation (i : Fin 2) (j : Fin M) : A →ₐ[K] A[ε] :=
  IsLocalization.Away.liftAlgHom (quadNorm a b) (firstOrderPolynomial_unit a b i j)

private theorem firstOrderRotation_polynomial (i : Fin 2) (j : Fin M) (f : R) :
    firstOrderRotation a b i j (algebraMap R A f) = firstOrderPolynomial a b i j f :=
  IsLocalization.Away.lift_eq (quadNorm a b) (firstOrderPolynomial_unit a b i j) f

theorem firstOrderRotation_norm (i : Fin 2) (j : Fin M) :
    firstOrderRotation a b i j (normUnit a b : A) = inl (normUnit a b : A) := by
  rw [normUnit_val, firstOrderRotation_polynomial, firstOrderPolynomial_norm, normUnit_val]

theorem firstOrderRotation_inv_norm (i : Fin 2) (j : Fin M) :
    firstOrderRotation a b i j (↑((normUnit a b)⁻¹) : A) = inl (↑((normUnit a b)⁻¹) : A) := by
  let ψ := firstOrderRotation a b i j
  let U := normUnit a b
  refine (U.isUnit.map ψ).mul_left_cancel ?_
  rw [← map_mul, Units.mul_inv, map_one, firstOrderRotation_norm, inl_mul_inl, Units.mul_inv,
    inl_one]

theorem firstOrderRotation_fraction (i : Fin 2) (j : Fin M) (f : R) (n : ℕ) :
    (firstOrderRotation a b i j (algebraMap R A f * (↑((normUnit a b)⁻¹) : A) ^ n)).snd =
      algebraMap R A (rotation a b i j f) * (↑((normUnit a b)⁻¹) : A) ^ n := by
  rw [map_mul, map_pow, firstOrderRotation_inv_norm, firstOrderRotation_polynomial]
  change ((firstOrderPolynomial a b i j) f *
    (TrivSqZeroExt.inlHom A A) (↑((normUnit a b)⁻¹) : A) ^ n).snd = _
  rw [← map_pow, DualNumber.snd_mul]
  change algebraMap R A f * 0 +
    algebraMap R A (rotation a b i j f) * (↑((normUnit a b)⁻¹) : A) ^ n = _
  rw [mul_zero, zero_add]

theorem annihilatesRotations_of_firstOrder_invariance (trace : A →ₗ[K] K)
    (h : ∀ (i : Fin 2) (j : Fin M) (f : A), trace (firstOrderRotation a b i j f).snd = 0) :
    AnnihilatesRotations a b trace := by
  intro i j f n
  rw [← firstOrderRotation_fraction]
  exact h i j _

end MathieuQuadratic
