import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import MathieuQuadratic.QuadraticCoordinates
import MathieuQuadratic.TraceUniqueness

/-!
# Independence of the quadratic trace from diagonal coordinates
-/

open MvPolynomial
open scoped Matrix

namespace MathieuQuadratic

variable {K I J T : Type*} [Field K] [CharZero K] [Infinite K]
  [Fintype I] [Fintype J] [Fintype T]

omit [CharZero K] [Infinite K] in
private theorem coordinateMatrix_equiv_inverse (e : (I → K) ≃ₗ[K] (J → K)) [DecidableEq I] :
    coordinateMatrix e.symm.toLinearMap * coordinateMatrix e.toLinearMap = 1 := by
  apply Matrix.mulVec_injective
  funext v
  rw [← Matrix.mulVec_mulVec, Matrix.one_mulVec, coordinateMatrix_mulVec, coordinateMatrix_mulVec,
    LinearEquiv.coe_coe, LinearEquiv.coe_coe, LinearEquiv.symm_apply_apply]

section Diagonal

variable [DecidableEq I] [DecidableEq J] [Nonempty I] [Nonempty J]
variable (d : I → Kˣ) (r : J → Kˣ)
variable (e : (I → K) ≃ₗ[K] (J → K))

omit [CharZero K] [DecidableEq I] [DecidableEq J] [Nonempty I] [Nonempty J] in
private theorem diagonal_norm_change (he : ∀ v, diagonalNorm r (e v) = diagonalNorm d v) :
    linearCoordinateEquiv e (diagonalNormPolynomial d) = diagonalNormPolynomial r := by
  apply MvPolynomial.funext
  intro v
  rw [eval_linearCoordinateEquiv, eval_diagonalNormPolynomial, eval_diagonalNormPolynomial,
    ← he (e.symm v), LinearEquiv.apply_symm_apply]

/-- The localized change of diagonal coordinates. -/
noncomputable def diagonalLocalizationEquiv (he : ∀ v, diagonalNorm r (e v) = diagonalNorm d v) :
    Localization.Away (diagonalNormPolynomial d) ≃ₐ[K]
      Localization.Away (diagonalNormPolynomial r) :=
  localizationCoordinateEquiv e _ _ (diagonal_norm_change d r e he)

/-- A linear isometry between diagonal forms preserves the normalized trace
on the two localizations. -/
theorem algebraicTrace_diagonalLocalizationEquiv (he : ∀ v, diagonalNorm r (e v) = diagonalNorm d v) :
    (algebraicTrace r).comp (diagonalLocalizationEquiv d r e he).toLinearMap = algebraicTrace d := by
  let A := Localization.Away (diagonalNormPolynomial d)
  let U := coordinateMatrix e.symm.toLinearMap
  have hn : ∀ v, diagonalNorm d (U *ᵥ v) = diagonalNorm r v := by
    intro v
    rw [coordinateMatrix_mulVec, LinearEquiv.coe_coe, ← he (e.symm v), LinearEquiv.apply_symm_apply]
  have him := inverseDiagonal_conjugate_between d r U (coordinateMatrix e.toLinearMap)
    (coordinateMatrix_equiv_inverse e) (diagonalConjugate_of_preserves_between d r U hn)
  have hc : Fintype.card I = Fintype.card J := by
    simpa [Module.finrank_fintype_fun_eq_card] using LinearEquiv.finrank_eq e
  apply LinearMap.ext
  intro f
  let φ := diagonalLocalizationEquiv d r e he
  let s := IsLocalization.Away.sec (diagonalNormPolynomial d) f
  have hs := IsLocalization.Away.sec_spec (diagonalNormPolynomial d) f
  change algebraicTrace r (φ f) = algebraicTrace d f
  apply algebraicTrace_fraction_between r d U hc him f (φ f) s.1 s.2 hs
  have hg := congrArg φ hs
  change φ (f * algebraMap (MvPolynomial I K) A (diagonalNormPolynomial d ^ s.2)) =
      φ (algebraMap (MvPolynomial I K) A s.1) at hg
  dsimp only [φ, diagonalLocalizationEquiv] at hg
  rw [map_mul, localizationCoordinateEquiv_polynomial, localizationCoordinateEquiv_polynomial,
    map_pow, diagonal_norm_change d r e he, linearCoordinateEquiv_matrixPullback] at hg
  exact hg

end Diagonal

namespace QuadraticDiagonalChart

variable {q : QuadraticForm K (I → K)}
variable [DecidableEq J] [Nonempty J]

/-- The trace transported to the original coordinates of the form. -/
noncomputable def trace (chart : QuadraticDiagonalChart q J) :
    Localization.Away (quadraticCoordinatePolynomial q) →ₗ[K] K :=
  (algebraicTrace chart.weights).comp chart.localization.toLinearMap

@[simp] theorem trace_one (chart : QuadraticDiagonalChart q J) : chart.trace 1 = 1 := by
  simp [trace]

theorem trace_localizationAt (chart : QuadraticDiagonalChart q J)
    (r : MvPolynomial J K) (hr : r = diagonalNormPolynomial chart.weights) :
    (algebraicTraceAt chart.weights r hr).comp (chart.localizationAt r hr).toLinearMap = chart.trace := by
  subst r
  rfl

/-- The normalized trace transported by any two diagonal charts is the same. -/
theorem trace_eq [DecidableEq T] [Nonempty T] (chart : QuadraticDiagonalChart q J)
    (other : QuadraticDiagonalChart q T) : chart.trace = other.trace := by
  let e := chart.coordinates.symm.trans other.coordinates
  have he (v : J → K) : diagonalNorm other.weights (e v) = diagonalNorm chart.weights v := by
    change diagonalNorm other.weights (other.coordinates (chart.coordinates.symm v)) = _
    rw [other.preserves, ← chart.preserves (chart.coordinates.symm v), LinearEquiv.apply_symm_apply]
  let h := diagonalLocalizationEquiv chart.weights other.weights e he
  have ht : h.toAlgHom.comp chart.localization.toAlgHom = other.localization.toAlgHom := by
    apply IsLocalization.algHom_ext (Submonoid.powers (quadraticCoordinatePolynomial q))
    apply AlgHom.ext
    intro p
    change h (chart.localization (algebraMap _ _ p)) = other.localization (algebraMap _ _ p)
    dsimp only [h, diagonalLocalizationEquiv, localization, localizationAt]
    rw [localizationCoordinateEquiv_polynomial, localizationCoordinateEquiv_polynomial,
      localizationCoordinateEquiv_polynomial, linearCoordinateEquiv_comp]
    have heq : chart.coordinates.trans e = other.coordinates := by
      ext x i
      simp [e]
    rw [heq]
  apply LinearMap.ext
  intro f
  have happ := LinearMap.congr_fun (algebraicTrace_diagonalLocalizationEquiv
    chart.weights other.weights e he) (chart.localization f)
  change algebraicTrace other.weights (h (chart.localization f)) =
    algebraicTrace chart.weights (chart.localization f) at happ
  rw [show h (chart.localization f) = other.localization f from AlgHom.congr_fun ht f] at happ
  exact happ.symm

end QuadraticDiagonalChart

/-- The normalized trace on the coordinate ring of a nondegenerate
quadratic form of positive dimension. -/
noncomputable def quadraticTrace (q : QuadraticForm K (I → K))
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft) (hpos : 0 < Fintype.card I) :
    Localization.Away (quadraticCoordinatePolynomial q) →ₗ[K] K := by
  letI : Nonempty (Fin (Fintype.card I)) := ⟨⟨0, hpos⟩⟩
  exact (Classical.choice (exists_fin_diagonal_chart q hn (Fintype.card I) rfl)).trace

/-- Every isometric diagonal chart computes the canonical quadratic trace. -/
theorem quadraticTrace_eq_chart {q : QuadraticForm K (I → K)}
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft) (hpos : 0 < Fintype.card I)
    [DecidableEq J] [Nonempty J] (chart : QuadraticDiagonalChart q J) :
    quadraticTrace q hn hpos = chart.trace := by
  letI : Nonempty (Fin (Fintype.card I)) := ⟨⟨0, hpos⟩⟩
  exact (Classical.choice (exists_fin_diagonal_chart q hn (Fintype.card I) rfl)).trace_eq chart

@[simp] theorem quadraticTrace_one {q : QuadraticForm K (I → K)}
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft) (hpos : 0 < Fintype.card I) :
    quadraticTrace q hn hpos 1 = 1 := by
  letI : Nonempty (Fin (Fintype.card I)) := ⟨⟨0, hpos⟩⟩
  exact (Classical.choice (exists_fin_diagonal_chart q hn (Fintype.card I) rfl)).trace_one

section GeometricNormalization

open scoped TensorProduct

variable {q : QuadraticForm K (I → K)} [DecidableEq J] [Nonempty J]

/-- Express the scalar extension of a functional in a diagonal chart
over an algebraic closure. -/
noncomputable def QuadraticDiagonalChart.geometricLift (chart : QuadraticDiagonalChart q J)
    (τ : Localization.Away (quadraticCoordinatePolynomial q) →ₗ[K] K) :
    Localization.Away (diagonalNormPolynomial (extendedDiagonal (L := AlgebraicClosure K) chart.weights))
      →ₗ[AlgebraicClosure K] AlgebraicClosure K :=
  (extendTrace (AlgebraicClosure K) (τ.comp chart.localization.symm.toLinearMap)).comp
    (diagonalComplementBaseChange chart.weights).symm.toLinearMap

/-- The canonical trace extends to the normalized geometric point trace
in every diagonal chart. -/
theorem quadraticTrace_geometricLift_eq (chart : QuadraticDiagonalChart q J)
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft) (hpos : 0 < Fintype.card I) :
    chart.geometricLift (quadraticTrace q hn hpos) =
      algebraicTrace (extendedDiagonal (L := AlgebraicClosure K) chart.weights) := by
  rw [quadraticTrace_eq_chart hn hpos chart]
  have he : chart.trace.comp chart.localization.symm.toLinearMap = algebraicTrace chart.weights := by
    apply LinearMap.ext
    intro f
    change algebraicTrace chart.weights (chart.localization (chart.localization.symm f)) = _
    rw [AlgEquiv.apply_symm_apply]
  rw [QuadraticDiagonalChart.geometricLift, he, ← algebraicTrace_baseChange]
  apply LinearMap.ext
  intro f
  simp only [LinearMap.comp_apply, AlgEquiv.toLinearMap_apply, AlgEquiv.apply_symm_apply]

/-- A normalized ground-field functional invariant under the geometric
points over the algebraic closure is the canonical trace. -/
theorem quadraticTrace_eq_of_geometricInvariant (chart : QuadraticDiagonalChart q J)
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft) (hpos : 0 < Fintype.card I)
    (τ : Localization.Away (quadraticCoordinatePolynomial q) →ₗ[K] K) (hτ : τ 1 = 1)
    (invariant : ∀ g f, chart.geometricLift τ
      (orthogonalScalarAction (extendedDiagonal (L := AlgebraicClosure K) chart.weights) g f) =
      chart.geometricLift τ f) :
    τ = quadraticTrace q hn hpos := by
  let L := AlgebraicClosure K
  let w := extendedDiagonal (L := L) chart.weights
  have h1 : chart.geometricLift τ 1 = 1 := by
    rw [QuadraticDiagonalChart.geometricLift, LinearMap.comp_apply, AlgEquiv.toLinearMap_apply, map_one]
    change extendTrace L (τ.comp chart.localization.symm.toLinearMap) (1 ⊗ₜ[K] 1) = 1
    rw [extendTrace_tmul, LinearMap.comp_apply, AlgEquiv.toLinearMap_apply, map_one, hτ, map_one, mul_one]
  have hgeom := (existsUnique_geometricPointTrace w).unique
    ⟨h1, invariant⟩ (algebraicTrace_normalized_invariant w)
  rw [← quadraticTrace_geometricLift_eq chart hn hpos] at hgeom
  apply LinearMap.ext
  intro f
  apply (RingHom.injective (algebraMap K L))
  let E := diagonalComplementBaseChange (L := L) chart.weights
  have he := LinearMap.congr_fun hgeom (E (1 ⊗ₜ[K] chart.localization f))
  simpa only [QuadraticDiagonalChart.geometricLift, LinearMap.comp_apply, AlgEquiv.toLinearMap_apply,
    AlgEquiv.symm_apply_apply, extendTrace_tmul, mul_one, E] using he

end GeometricNormalization

end MathieuQuadratic
