import MathieuQuadratic.Classification
import MathieuQuadratic.Extension
import MathieuQuadratic.LaurentBaseChange
import MathieuQuadratic.OrthogonalTrace
import MathieuQuadratic.QuadraticCoordinates
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure

/-!
# Statements of The Mathieu property for quadratic complements

This file gives the manuscript statements on the localized coordinate
rings. The mixed-moment proposition, the ternary example and the explicit
witness use geometric trace invariance. The dimension theorem also uses
the low-dimensional trace comparisons and the complex Laurent theorem.
-/

open scoped BigOperators TensorProduct Interval
open Algebra.TensorProduct

namespace MathieuQuadratic

variable {K : Type*} [Field K] [CharZero K]

/-- Termwise integration of a polynomial from zero to one over K. -/
noncomputable def polynomialIntegral : Polynomial K →ₗ[K] K :=
  Polynomial.lsum (fun d ↦ LinearMap.mulRight K (1 / ((d : K) + 1)))

omit [CharZero K] in
@[simp] theorem polynomialIntegral_monomial (d : ℕ) (c : K) :
    polynomialIntegral (Polynomial.monomial d c) = c / ((d : K) + 1) := by
  rw [polynomialIntegral, Polynomial.lsum_apply, Polynomial.sum_monomial_index]
  · simp [LinearMap.mulRight_apply, div_eq_mul_inv]
  · exact map_zero _

section Diagonal

variable {M : ℕ} (a : Fin 2 → K) (b : Fin M → K)
variable {trace : NormComplement a b →ₗ[K] K}

/-- The draft's height integral, with its polynomial integrand over K. -/
theorem draft_height_integral (h : GeometricTrace a b trace) (hM : 0 < M)
    (g : Polynomial K) :
    trace (Polynomial.aeval (height a b) g) =
      (M : K) * polynomialIntegral ((Polynomial.X : Polynomial K) ^ (M - 1) *
        g.comp (Polynomial.X ^ 2)) := by
  induction g using Polynomial.induction_on' with
  | add f g hf hg =>
    rw [map_add, map_add, Polynomial.add_comp, mul_add, map_add, mul_add, hf, hg]
  | monomial d c =>
    have hd : M - 1 + 2 * d + 1 = M + 2 * d := by omega
    have hcoeff : (↑(M - 1 + 2 * d) : K) + 1 = (↑(M + 2 * d) : K) := by
      exact_mod_cast hd
    rw [Polynomial.aeval_monomial, ← Algebra.smul_def, map_smul, smul_eq_mul,
      (h.complementTrace.normalizer a b hM).height,
      ← Polynomial.C_mul_X_pow_eq_monomial, Polynomial.mul_comp, Polynomial.C_comp,
      Polynomial.pow_comp, Polynomial.X_comp, ← pow_mul, mul_left_comm, ← pow_add,
      Polynomial.C_mul_X_pow_eq_monomial, polynomialIntegral_monomial, hcoeff]
    ring

/-- The product phase in the manuscript's moment proposition. -/
noncomputable def complementProductPhase (α β γ : ℕ) : NormComplement a b :=
  productPhase (normUnit a b) (height a b) α β γ

/-- The geometric mixed-moment formula, including negative-index vanishing. -/
theorem draft_mixed_moments (h : GeometricTrace a b trace) (hM : 0 < M)
    (α β γ m j : ℕ) (ℓ : ℤ) (hℓ : (α : ℤ) * m - (γ : ℤ) * m ≤ ℓ) :
    trace ((↑(normUnit a b ^ ℓ) : NormComplement a b) * height a b ^ j *
      complementProductPhase a b α β γ ^ m) =
      if ℓ ≤ (α * m : ℕ) then algebraMap ℚ K
        (mixedBetaMoment M (γ * m) j *
          ((Ring.choose ((β * m : ℕ) - (M + 2 * j : ℕ) : ℤ)
            ((α * m : ℕ) - ℓ : ℤ).toNat : ℤ) : ℚ))
      else 0 :=
  product_mixedMoment trace _ _ M α β γ m j hM ℓ hℓ
    (h.complementTrace.normalizer a b hM).weight (h.complementTrace.normalizer a b hM).height

/-- The first discrepancy in the geometric mixed-moment formula. -/
theorem draft_first_error (h : GeometricTrace a b trace) (hM : 0 < M)
    (α β γ m j : ℕ) :
    trace ((↑(normUnit a b ^ ((α : ℤ) * m - (γ : ℤ) * m - 1)) : NormComplement a b) *
      height a b ^ j * complementProductPhase a b α β γ ^ m) =
      algebraMap ℚ K (mixedBetaMoment M (γ * m) j *
        ((Ring.choose ((β * m : ℕ) - (M + 2 * j : ℕ) : ℤ) (γ * m + 1) : ℤ) : ℚ) +
          (M : ℚ) * (-2) ^ (γ * m) / (γ * m + 1 : ℕ)) :=
  product_mixedMoment_first_error trace _ _ M α β γ m j hM
    (h.complementTrace.normalizer a b hM).weight (h.complementTrace.normalizer a b hM).height

/-- The pure and marked moments for the manuscript's two witnesses. -/
theorem draft_product_witness (h : GeometricTrace a b trace) (hM : 0 < M)
    (γ : ℕ) (hγ : M ≤ γ) :
    let F := complementProductPhase a b M M γ
    let H := (normUnit a b : NormComplement a b) ^ M
    AllOrderWitness trace F H ∧ ∀ m : ℕ, 0 < m →
      trace (H * F ^ m) = algebraMap ℚ K (mixedBetaMoment M (γ * m) 0) := by
  have ht := h.complementTrace.normalizer a b hM
  have hp (m : ℕ) (hm : 0 < m) :=
    product_uniform_moments trace _ _ M γ m hM hγ hm ht.weight ht.height
  exact ⟨⟨fun m hm ↦ (hp m hm).1, fun m hm ↦ (hp m hm).2.trans_ne
      (by
        apply (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
        exact (mixedBetaMoment_pos M (γ * m) 0 hM).ne')⟩,
    fun m hm ↦ (hp m hm).2⟩

/-- The explicit nonpower phase and prescribed multiplier of the main theorem. -/
theorem draft_explicit_counterexample (h : GeometricTrace a b trace) (hM : 0 < M)
    (ha : a 0 ≠ 0) (hb : ∀ i, b i ≠ 0) :
    let F := complementProductPhase a b M M (M + 1)
    let H := (normUnit a b : NormComplement a b) ^ M
    NonpowerUpToUnits (includeRight (R := K) (A := AlgebraicClosure K) F) ∧
      AllOrderWitness trace F H ∧
      ∀ m : ℕ, 0 < m → trace (H * F ^ m) = algebraMap ℚ K (mixedBetaMoment M ((M + 1) * m) 0) := by
  let i : Fin M := ⟨0, hM⟩
  exact ⟨productPhase_absolute_nonpower a b ha i (hb i) M M (M + 1) (by simp),
    draft_product_witness a b h hM (M + 1) (by omega)⟩

/-- The same product phase and marker in the scalar extension. -/
theorem draft_same_witness_after_extension (L : Type*) [Field L] [Algebra K L]
    (h : GeometricTrace a b trace) (hM : 0 < M) (γ : ℕ) (hγ : M ≤ γ) :
    AllOrderWitness ((extendTrace L trace).comp (complementBaseChange (L := L) a b).symm.toLinearMap)
      (complementProductPhase (algebraMap K L ∘ a) (algebraMap K L ∘ b) M M γ)
      ((normUnit (algebraMap K L ∘ a) (algebraMap K L ∘ b) :
        NormComplement (algebraMap K L ∘ a) (algebraMap K L ∘ b)) ^ M) := by
  rw [complementProductPhase, ← complementBaseChange_productPhase a b,
    ← complementBaseChange_marker a b]
  obtain ⟨hp, hm⟩ := (allOrderWitness_includeRight_iff L trace _ _).2
    (draft_product_witness a b h hM γ hγ).1
  constructor
  · intro n hn
    rw [LinearMap.comp_apply, AlgEquiv.toLinearMap_apply, ← map_pow, AlgEquiv.symm_apply_apply]
    exact hp n hn
  · intro n hn
    rw [LinearMap.comp_apply, AlgEquiv.toLinearMap_apply, map_mul, ← map_pow,
      AlgEquiv.symm_apply_apply, AlgEquiv.symm_apply_apply]
    exact hm n hn

end Diagonal

/-- The lower-pole ternary witness in the opening example. -/
theorem draft_ternary_phase (a : Fin 2 → K) (b : Fin 1 → K)
    {trace : NormComplement a b →ₗ[K] K} (h : GeometricTrace a b trace)
    (ha : a 0 ≠ 0) (hb : b 0 ≠ 0) :
    let F₀ := complementProductPhase a b 1 1 1
    NonpowerUpToUnits (includeRight (R := K) (A := AlgebraicClosure K) F₀) ∧
      AllOrderWitness trace F₀ (normUnit a b : NormComplement a b) := by
  exact ⟨productPhase_absolute_nonpower a b ha 0 hb 1 1 1 (by decide),
    by simpa using (draft_product_witness a b h Nat.one_pos 1 le_rfl).1⟩

/-- The opening moments for x²+y²+z² and F₀ = F₁,₁,₁. -/
theorem draft_opening
    {trace : NormComplement (fun _ : Fin 2 ↦ (1 : K)) (fun _ : Fin 1 ↦ (1 : K)) →ₗ[K] K}
    (h : GeometricTrace (fun _ : Fin 2 ↦ (1 : K)) (fun _ : Fin 1 ↦ (1 : K)) trace)
    (m : ℕ) (hm : 0 < m) :
    let a : Fin 2 → K := fun _ ↦ 1
    let b : Fin 1 → K := fun _ ↦ 1
    let F₀ := complementProductPhase a b 1 1 1
    let q : NormComplement a b := normUnit a b
    trace (F₀ ^ m) = 0 ∧ trace (q * F₀ ^ m) = algebraMap ℚ K (mixedBetaMoment 1 m 0) ∧
      trace (q * F₀ ^ m) ≠ 0 ∧
      (mixedBetaMoment 1 m 0 : ℝ) = ∫ x in (0 : ℝ)..1, (1 - x ^ 2) ^ m := by
  obtain ⟨hw, moments⟩ := draft_product_witness
    (fun _ : Fin 2 ↦ (1 : K)) (fun _ : Fin 1 ↦ (1 : K)) h Nat.one_pos 1 le_rfl
  exact ⟨hw.1 m hm, by simpa using moments m hm, by simpa using hw.2 m hm, by
    simpa using mixedBetaMoment_cast_eq_integral 1 m 0 Nat.one_pos⟩

/-- The first marked moment over the rationals. -/
theorem draft_opening_first_beta : mixedBetaMoment 1 1 0 = 2 / 3 := by
  norm_num [mixedBetaMoment, powerBetaMoment, powerBetaAntiderivative, powerBetaQuotient,
    powerBetaCoefficient, Finset.sum_range_succ]


section Classification

variable {A : Type*} [CommRing A] [Algebra K A]

/-- The higher-dimensional trace comparison also identifies the norm
unit. This keeps the fixed multiplier in the statement of the draft. -/
structure DiagonalTraceModel (M : ℕ) (q : Aˣ) (trace : A →ₗ[K] K) where
  plane : Fin 2 → K
  heightCoeffs : Fin M → K
  coordinates : A ≃ₐ[K] NormComplement plane heightCoeffs
  norm_eq : coordinates (q : A) = (normUnit plane heightCoeffs : NormComplement plane heightCoeffs)
  plane_nonzero : plane 0 ≠ 0
  height_nonzero : ∀ i, heightCoeffs i ≠ 0
  geometric : GeometricTrace plane heightCoeffs (trace.comp coordinates.symm.toLinearMap)

theorem DiagonalTraceModel.counterexample {M : ℕ} {q : Aˣ} {trace : A →ₗ[K] K}
    (model : DiagonalTraceModel M q trace) (hM : 0 < M) :
    let F := model.coordinates.symm (complementProductPhase model.plane model.heightCoeffs M M (M + 1))
    NonpowerUpToUnits (includeRight (R := K) (A := AlgebraicClosure K) F) ∧
      AllOrderWitness trace F ((q : A) ^ M) := by
  obtain ⟨hn, hw, _⟩ := draft_explicit_counterexample
    model.plane model.heightCoeffs model.geometric hM model.plane_nonzero model.height_nonzero
  have he : model.coordinates.symm (normUnit model.plane model.heightCoeffs : NormComplement
      model.plane model.heightCoeffs) = (q : A) := by
    rw [← model.norm_eq, AlgEquiv.symm_apply_apply]
  have hp := transport_nonpower_and_witness model.coordinates trace _ _ hn hw
  have hm : model.coordinates.symm ((normUnit model.plane model.heightCoeffs :
      NormComplement model.plane model.heightCoeffs) ^ M) = (q : A) ^ M := by
    rw [map_pow, he]
  rw [hm] at hp
  exact hp

/-- Theorem 1 for a norm coordinate ring with its trace comparisons. The
large-rank conclusion retains the prescribed norm-power marker. -/
theorem draft_dimension_classification (N : ℕ) (hN : 0 < N)
    (q : Aˣ) (trace : A →ₗ[K] K) (dvk : LowRankLaurentTheorem)
    (split : N ≤ 2 →
      ∃ coordinates : AlgebraicClosure K ⊗[K] A ≃ₐ[AlgebraicClosure K]
          AddMonoidAlgebra (AlgebraicClosure K) (Fin N → ℤ),
        ∀ f, extendTrace (AlgebraicClosure K) trace f = laurentCT (coordinates f))
    (diagonal : 3 ≤ N → Nonempty (DiagonalTraceModel (N - 2) q trace)) :
    (IsMathieu trace ↔ N ≤ 2) ∧
      (3 ≤ N → ∃ F : A, NonpowerUpToUnits (includeRight (R := K) (A := AlgebraicClosure K) F) ∧
        AllOrderWitness trace F ((q : A) ^ (N - 2))) := by
  have high (hr : 3 ≤ N) :
      ∃ F : A, NonpowerUpToUnits (includeRight (R := K) (A := AlgebraicClosure K) F) ∧
        AllOrderWitness trace F ((q : A) ^ (N - 2)) := by
    obtain ⟨model⟩ := diagonal hr
    exact ⟨_, model.counterexample (by omega)⟩
  refine ⟨⟨?_, ?_⟩, high⟩
  · intro h
    by_contra hn
    obtain ⟨_, _, witness⟩ := high (by omega)
    exact not_isMathieu_of_allOrderWitness witness h
  · intro hn
    obtain ⟨coordinates, ht⟩ := split hn
    exact isMathieu_of_split_laurent trace coordinates ht
      (laurentSeparation_of_complex (laurentSeparation_of_newton (dvk N (by omega))))

/-- The scalar-extension corollary of the draft with a single set of
trace comparisons over K. The splitting field is allowed to differ from
the extension field. -/
theorem draft_mathieu_iff_after_extension (L : Type*) [Field L] [Algebra K L]
    (N : ℕ) (hN : 0 < N) (q : Aˣ) (trace : A →ₗ[K] K) (dvk : LowRankLaurentTheorem)
    (split : N ≤ 2 →
      ∃ coordinates : AlgebraicClosure K ⊗[K] A ≃ₐ[AlgebraicClosure K]
          AddMonoidAlgebra (AlgebraicClosure K) (Fin N → ℤ),
        ∀ f, extendTrace (AlgebraicClosure K) trace f = laurentCT (coordinates f))
    (diagonal : 3 ≤ N → Nonempty (DiagonalTraceModel (N - 2) q trace)) :
    (IsMathieu (extendTrace L trace) ↔ IsMathieu trace) ∧
      (3 ≤ N → ∃ F : A,
        AllOrderWitness trace F ((q : A) ^ (N - 2)) ∧
        AllOrderWitness (extendTrace L trace) (includeRight (A := L) F)
          (includeRight (A := L) ((q : A) ^ (N - 2)))) := by
  obtain ⟨base, high⟩ := draft_dimension_classification N hN q trace dvk split diagonal
  refine ⟨⟨isMathieu_of_extended L trace, ?_⟩, ?_⟩
  · intro h
    have hn := base.1 h
    obtain ⟨coordinates, ht⟩ := split hn
    exact isMathieu_after_extension_of_split_laurent trace coordinates ht
      (laurentSeparation_of_newton (dvk N (by omega)))
  · intro hn
    obtain ⟨F, _, hw⟩ := high hn
    exact ⟨F, hw, (allOrderWitness_includeRight_iff L trace _ _).2 hw⟩

end Classification

section ArbitraryQuadraticForm

open MvPolynomial

variable {I : Type*} [Fintype I] [Infinite K]

/-- The unit represented by the norm on its quadratic complement. -/
noncomputable def quadraticNormUnit (q : QuadraticForm K (I → K)) :
    (Localization.Away (quadraticCoordinatePolynomial q))ˣ :=
  (IsLocalization.Away.algebraMap_isUnit (quadraticCoordinatePolynomial q)).unit

/-- Diagonalization supplies the coordinates and the prescribed multiplier
for an arbitrary nondegenerate form of dimension M + 2. -/
theorem exists_DiagonalTraceModel_of_form (q : QuadraticForm K (I → K))
    (hn : (QuadraticMap.associated (R := K) q).SeparatingLeft)
    (M : ℕ) (hdim : Fintype.card I = M + 2)
    (trace : Localization.Away (quadraticCoordinatePolynomial q) →ₗ[K] K)
    (compatible : ∀ (a : Fin 2 → Kˣ) (b : Fin M → Kˣ)
      (e : Localization.Away (quadraticCoordinatePolynomial q) ≃ₐ[K]
        NormComplement (fun i ↦ (a i : K)) (fun j ↦ (b j : K))),
      e (quadraticNormUnit q : Localization.Away (quadraticCoordinatePolynomial q)) =
        (normUnit (fun i ↦ (a i : K)) (fun j ↦ (b j : K)) :
          NormComplement (fun i ↦ (a i : K)) (fun j ↦ (b j : K))) →
      GeometricTrace (fun i ↦ (a i : K)) (fun j ↦ (b j : K)) (trace.comp e.symm.toLinearMap)) :
    Nonempty (DiagonalTraceModel M (quadraticNormUnit q) trace) := by
  obtain ⟨a, b, e, he⟩ := exists_diagonal_localization q hn M hdim
  have hunit : e (quadraticNormUnit q : Localization.Away (quadraticCoordinatePolynomial q)) =
      (normUnit (fun i ↦ (a i : K)) (fun j ↦ (b j : K)) :
        NormComplement (fun i ↦ (a i : K)) (fun j ↦ (b j : K))) := by
    have hq : (quadraticNormUnit q : Localization.Away (quadraticCoordinatePolynomial q)) =
        algebraMap (MvPolynomial I K) _ (quadraticCoordinatePolynomial q) := IsUnit.unit_spec _
    rw [hq, he]
  exact ⟨⟨fun i ↦ a i, fun j ↦ b j, e, hunit, (a 0).ne_zero, fun i ↦ (b i).ne_zero,
    compatible a b e hunit⟩⟩

end ArbitraryQuadraticForm

end MathieuQuadratic
