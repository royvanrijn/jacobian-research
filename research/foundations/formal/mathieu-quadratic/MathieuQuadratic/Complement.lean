import MathieuQuadratic.Uniform
import Mathlib.Algebra.MvPolynomial.PDeriv
import Mathlib.RingTheory.Localization.Away.Basic

/-!
# The diagonal quadratic complement

All coordinates and their norm are defined over the given field. The
coordinate ring inverts the single polynomial q.
-/

open scoped BigOperators

namespace MathieuQuadratic

variable {K : Type*} [Field K] {M : ℕ}

abbrev QuadPolynomial (K : Type*) [Field K] (M : ℕ) :=
  MvPolynomial (Fin 2 ⊕ Fin M) K

noncomputable def planeNorm (a : Fin 2 → K) : QuadPolynomial K M :=
  ∑ i, MvPolynomial.C (a i) * MvPolynomial.X (Sum.inl i) ^ 2

noncomputable def heightNorm (b : Fin M → K) : QuadPolynomial K M :=
  ∑ i, MvPolynomial.C (b i) * MvPolynomial.X (Sum.inr i) ^ 2

noncomputable def quadNorm (a : Fin 2 → K) (b : Fin M → K) : QuadPolynomial K M :=
  planeNorm a + heightNorm b

abbrev NormComplement (a : Fin 2 → K) (b : Fin M → K) :=
  Localization.Away (quadNorm a b)

noncomputable def normUnit (a : Fin 2 → K) (b : Fin M → K) : (NormComplement a b)ˣ :=
  (IsLocalization.Away.algebraMap_isUnit (S := NormComplement a b) (quadNorm a b)).unit

@[simp] theorem normUnit_val (a : Fin 2 → K) (b : Fin M → K) :
    (normUnit a b : NormComplement a b) =
      algebraMap (QuadPolynomial K M) (NormComplement a b) (quadNorm a b) :=
  IsUnit.unit_spec _

noncomputable def height (a : Fin 2 → K) (b : Fin M → K) : NormComplement a b :=
  algebraMap (QuadPolynomial K M) (NormComplement a b) (heightNorm b) *
    (↑((normUnit a b)⁻¹) : NormComplement a b)

end MathieuQuadratic
