import MathieuQuadratic.Diagonal
import MathieuQuadratic.GeometricAction
import Mathlib.Algebra.MvPolynomial.Funext

/-!
# Independence of the coordinate model

A linear change of variables preserving the quadratic form induces
an isomorphism of its norm localizations. The high-dimensional chart
in the draft is obtained by the rank-two orthogonal splitting.
-/

open scoped BigOperators Matrix
open MvPolynomial

namespace MathieuQuadratic

variable {K I J : Type*} [Field K] [CharZero K] [Infinite K]
  [Fintype I] [Fintype J]

private noncomputable def standardVector {D : Type*} (i : D) : D → K := by
  classical
  exact fun j ↦ if i = j then 1 else 0

private noncomputable def linearSubstitution (e : (J → K) →ₗ[K] (I → K)) :
    MvPolynomial I K →ₐ[K] MvPolynomial J K :=
  aeval fun i ↦ ∑ j, C (e (standardVector j) i) * X j

omit [CharZero K] [Infinite K] [Fintype I] in
private theorem eval_linearSubstitution (e : (J → K) →ₗ[K] (I → K))
    (x : J → K) (p : MvPolynomial I K) :
    eval x (linearSubstitution e p) = eval (e x) p := by
  classical
  have he : (aeval x : MvPolynomial J K →ₐ[K] K).comp (linearSubstitution e) =
      aeval (e x) := by
    apply MvPolynomial.algHom_ext
    intro i
    simp only [AlgHom.comp_apply, linearSubstitution, aeval_X, map_sum, map_mul, aeval_C,
      Algebra.algebraMap_self, RingHom.id_apply]
    have hel := congrArg (fun v : I → K ↦ v i) (e.pi_apply_eq_sum_univ x)
    unfold standardVector
    simpa [mul_comm, smul_eq_mul] using hel.symm
  exact AlgHom.congr_fun he p

/-- A linear coordinate equivalence gives an isomorphism of the
polynomial coordinate rings, by substitution along its inverse. -/
noncomputable def linearCoordinateEquiv (e : (I → K) ≃ₗ[K] (J → K)) :
    MvPolynomial I K ≃ₐ[K] MvPolynomial J K := by
  refine AlgEquiv.ofAlgHom (linearSubstitution e.symm.toLinearMap)
    (linearSubstitution e.toLinearMap) ?_ ?_
  · apply AlgHom.ext
    intro p
    apply MvPolynomial.funext
    intro x
    simp [eval_linearSubstitution]
  · apply AlgHom.ext
    intro p
    apply MvPolynomial.funext
    intro x
    simp [eval_linearSubstitution]

omit [CharZero K] in
theorem eval_linearCoordinateEquiv (e : (I → K) ≃ₗ[K] (J → K))
    (x : J → K) (p : MvPolynomial I K) :
    eval x (linearCoordinateEquiv e p) = eval (e.symm x) p :=
  eval_linearSubstitution e.symm.toLinearMap x p

/-- The coordinate polynomial associated with a quadratic form.
Polarization supplies the entries of its symmetric matrix. -/
noncomputable def quadraticCoordinatePolynomial (q : QuadraticForm K (I → K)) :
    MvPolynomial I K :=
  letI : Invertible (2 : K) := invertibleOfNonzero two_ne_zero
  ∑ i, ∑ j, C (QuadraticMap.associated (R := K) q (standardVector i) (standardVector j)) *
    X i * X j

omit [Infinite K] in
theorem eval_quadraticCoordinatePolynomial (q : QuadraticForm K (I → K)) (v : I → K) :
    eval v (quadraticCoordinatePolynomial q) = q v := by
  classical
  letI : Invertible (2 : K) := invertibleOfNonzero two_ne_zero
  let B := QuadraticMap.associated (R := K) q
  rw [← QuadraticMap.associated_eq_self_apply K q v]
  change eval v (∑ i, ∑ j, C (B (standardVector i) (standardVector j)) * X i * X j) = B v v
  have hl : B v v = ∑ i, v i * B (standardVector i) v := by
    unfold standardVector
    simpa [smul_eq_mul] using congrArg (fun f : (I → K) →ₗ[K] K ↦ f v)
      (B.pi_apply_eq_sum_univ v)
  have hr (i : I) : B (standardVector i) v =
      ∑ j, v j * B (standardVector i) (standardVector j) := by
    let f := B (standardVector i)
    change f v = ∑ j, v j * f (standardVector j)
    unfold standardVector
    simpa [smul_eq_mul] using f.pi_apply_eq_sum_univ v
  simp only [map_sum, map_mul, eval_C, eval_X]
  simp_rw [hl, hr, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ Finset.sum_congr rfl fun j _ ↦ ?_
  ring

/-- Equivalent quadratic forms have equivalent norm localizations. -/
noncomputable def localizationOfIsometry
    {q : QuadraticForm K (I → K)} {q' : QuadraticForm K (J → K)} (e : q.IsometryEquiv q') :
    Localization.Away (quadraticCoordinatePolynomial q) ≃ₐ[K]
      Localization.Away (quadraticCoordinatePolynomial q') := by
  have hp : linearCoordinateEquiv e.toLinearEquiv (quadraticCoordinatePolynomial q) =
      quadraticCoordinatePolynomial q' := by
    apply MvPolynomial.funext
    intro v
    rw [eval_linearCoordinateEquiv, eval_quadraticCoordinatePolynomial,
      eval_quadraticCoordinatePolynomial]
    exact e.symm.map_app v
  exact IsLocalization.algEquivOfAlgEquiv (M := Submonoid.powers (quadraticCoordinatePolynomial q))
    (T := Submonoid.powers (quadraticCoordinatePolynomial q'))
    _ _ (linearCoordinateEquiv e.toLinearEquiv) (by rw [Submonoid.map_powers, hp])

omit [CharZero K] [Infinite K] [Fintype I] [Fintype J] in
/-- The matrix whose columns are the coordinate images under a linear map. -/
noncomputable def coordinateMatrix (e : (J → K) →ₗ[K] (I → K)) : Matrix I J K :=
  fun i j ↦ e (standardVector j) i

omit [CharZero K] [Infinite K] [Fintype I] in
theorem coordinateMatrix_mulVec (e : (J → K) →ₗ[K] (I → K)) (v : J → K) :
    coordinateMatrix e *ᵥ v = e v := by
  classical
  funext i
  have h := congrArg (fun w : I → K ↦ w i) (e.pi_apply_eq_sum_univ v)
  change (∑ j, e (standardVector j) i * v j) = e v i
  unfold standardVector
  simpa [mul_comm, smul_eq_mul] using h.symm

omit [CharZero K] in
theorem linearCoordinateEquiv_matrixPullback (e : (I → K) ≃ₗ[K] (J → K)) (p : MvPolynomial I K) :
    linearCoordinateEquiv e p = matrixPolynomialPullback (coordinateMatrix e.symm.toLinearMap) p := rfl

omit [CharZero K] in
theorem linearCoordinateEquiv_comp {T : Type*} [Fintype T]
    (e : (I → K) ≃ₗ[K] (J → K)) (g : (J → K) ≃ₗ[K] (T → K)) (p : MvPolynomial I K) :
    linearCoordinateEquiv g (linearCoordinateEquiv e p) = linearCoordinateEquiv (e.trans g) p := by
  apply MvPolynomial.funext
  intro x
  rw [eval_linearCoordinateEquiv, eval_linearCoordinateEquiv, eval_linearCoordinateEquiv]
  rfl

omit [CharZero K] in
/-- A polynomial coordinate equivalence carrying one norm to another
induces the corresponding equivalence on localizations. -/
noncomputable def localizationCoordinateEquiv (e : (I → K) ≃ₗ[K] (J → K))
    (p : MvPolynomial I K) (r : MvPolynomial J K) (h : linearCoordinateEquiv e p = r) :
    Localization.Away p ≃ₐ[K] Localization.Away r :=
  IsLocalization.algEquivOfAlgEquiv _ _ (linearCoordinateEquiv e)
    (by rw [Submonoid.map_powers, h])

omit [CharZero K] in
theorem localizationCoordinateEquiv_polynomial (e : (I → K) ≃ₗ[K] (J → K))
    (p : MvPolynomial I K) (r : MvPolynomial J K) (h : linearCoordinateEquiv e p = r)
    (s : MvPolynomial I K) :
    localizationCoordinateEquiv e p r h (algebraMap _ (Localization.Away p) s) =
      algebraMap _ (Localization.Away r) (linearCoordinateEquiv e s) := by
  rw [localizationCoordinateEquiv, IsLocalization.algEquivOfAlgEquiv_eq]

/-- A diagonal chart records the linear isometry before localization. -/
structure QuadraticDiagonalChart (q : QuadraticForm K (I → K)) (J : Type*) [Fintype J] where
  weights : J → Kˣ
  coordinates : (I → K) ≃ₗ[K] (J → K)
  preserves : ∀ v, diagonalNorm weights (coordinates v) = q v

namespace QuadraticDiagonalChart

variable {q : QuadraticForm K (I → K)} (chart : QuadraticDiagonalChart q J)

theorem norm_polynomial :
    linearCoordinateEquiv chart.coordinates (quadraticCoordinatePolynomial q) =
      diagonalNormPolynomial chart.weights := by
  apply MvPolynomial.funext
  intro v
  rw [eval_linearCoordinateEquiv, eval_quadraticCoordinatePolynomial, eval_diagonalNormPolynomial,
    ← chart.preserves (chart.coordinates.symm v), LinearEquiv.apply_symm_apply]

/-- The localized chart with its norm in any given polynomial presentation. -/
noncomputable def localizationAt (r : MvPolynomial J K)
    (hr : r = diagonalNormPolynomial chart.weights) :
    Localization.Away (quadraticCoordinatePolynomial q) ≃ₐ[K] Localization.Away r :=
  localizationCoordinateEquiv chart.coordinates _ r (chart.norm_polynomial.trans hr.symm)

noncomputable def localization :
    Localization.Away (quadraticCoordinatePolynomial q) ≃ₐ[K]
      Localization.Away (diagonalNormPolynomial chart.weights) :=
  chart.localizationAt (diagonalNormPolynomial chart.weights) rfl

end QuadraticDiagonalChart

omit [Infinite K] in
/-- A nondegenerate form admits a diagonal chart indexed by its dimension. -/
theorem exists_fin_diagonal_chart (q : QuadraticForm K (I → K))
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft)
    (N : ℕ) (hN : Fintype.card I = N) :
    Nonempty (QuadraticDiagonalChart q (Fin N)) := by
  have hdim : Module.finrank K (I → K) = N := by
    simpa [Module.finrank_fintype_fun_eq_card] using hN
  clear hN
  subst N
  letI : Invertible (2 : K) := invertibleOfNonzero two_ne_zero
  obtain ⟨d, ⟨e⟩⟩ := q.equivalent_weightedSumSquares_units_of_nondegenerate' hn
  refine ⟨⟨d, e.toLinearEquiv, fun v ↦ ?_⟩⟩
  change diagonalNorm d (e v) = q v
  simpa only [QuadraticMap.weightedSumSquares_apply, diagonalNorm, diagonalPair,
    Units.smul_def, smul_eq_mul, mul_assoc] using e.map_app v

/-- Diagonalization and transport give a norm localization with one
coordinate per dimension. -/
theorem exists_fin_diagonal_localization (q : QuadraticForm K (I → K))
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft)
    (N : ℕ) (hN : Fintype.card I = N) :
    ∃ d : Fin N → Kˣ, Nonempty (Localization.Away (quadraticCoordinatePolynomial q) ≃ₐ[K]
      Localization.Away (diagonalNormPolynomial d)) := by
  obtain ⟨chart⟩ := exists_fin_diagonal_chart q hn N hN
  exact ⟨chart.weights, ⟨chart.localization⟩⟩

/-- The higher-dimensional coordinate ring in the draft admits the
displayed rank-two diagonal chart, by a norm-preserving isomorphism. -/
theorem exists_diagonal_localization
    (q : QuadraticForm K (I → K))
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft)
    (M : ℕ) (hM : Fintype.card I = M + 2) :
    ∃ (a : Fin 2 → Kˣ) (b : Fin M → Kˣ)
      (e : Localization.Away (quadraticCoordinatePolynomial q) ≃ₐ[K]
        NormComplement (fun i ↦ (a i : K)) (fun j ↦ (b j : K))),
      e (algebraMap (MvPolynomial I K) _ (quadraticCoordinatePolynomial q)) =
        (normUnit (fun i ↦ (a i : K)) (fun j ↦ (b j : K)) :
          NormComplement (fun i ↦ (a i : K)) (fun j ↦ (b j : K))) := by
  have hdim : Module.finrank K (I → K) = M + 2 := by simpa [Module.finrank_fintype_fun_eq_card] using hM
  obtain ⟨a, b, ⟨iso⟩⟩ := exists_plane_height_diagonal q hn M hdim
  let p := quadNorm (fun i ↦ (a i : K)) (fun j ↦ (b j : K))
  have hp : linearCoordinateEquiv iso.toLinearEquiv (quadraticCoordinatePolynomial q) = p := by
    apply MvPolynomial.funext
    intro v
    rw [eval_linearCoordinateEquiv, eval_quadraticCoordinatePolynomial]
    change q (iso.symm v) = _
    rw [eval_quadNorm]
    exact iso.symm.map_app v
  let e : Localization.Away (quadraticCoordinatePolynomial q) ≃ₐ[K] Localization.Away p :=
    IsLocalization.algEquivOfAlgEquiv _ _ (linearCoordinateEquiv iso.toLinearEquiv)
      (by rw [Submonoid.map_powers, hp])
  refine ⟨a, b, e, ?_⟩
  dsimp only [e]
  rw [IsLocalization.algEquivOfAlgEquiv_eq, hp, normUnit_val]

end MathieuQuadratic
