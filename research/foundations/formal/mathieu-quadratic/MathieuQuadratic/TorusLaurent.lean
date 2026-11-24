import MathieuQuadratic.Laurent
import Mathlib.Analysis.Fourier.AddCircle
import Mathlib.MeasureTheory.Integral.Prod

/-!
# The torus integral of a planar Laurent polynomial

The normalized torus mean gives the Laurent constant term. The geometric
series under the same integral gives the generating function in Section 4
of Duistermaat–van der Kallen.
-/

open MeasureTheory

namespace MathieuQuadratic

local instance : Fact (0 < (1 : ℝ)) := ⟨zero_lt_one⟩

/-- The compact rank-two torus, in additive coordinates of period one. -/
abbrev LaurentTorus := AddCircle (1 : ℝ) × AddCircle (1 : ℝ)

/-- Normalized Haar measure on the compact rank-two torus. -/
noncomputable def laurentTorusHaar : Measure LaurentTorus :=
  AddCircle.haarAddCircle.prod AddCircle.haarAddCircle

instance : IsProbabilityMeasure laurentTorusHaar := by
  unfold laurentTorusHaar
  infer_instance

private noncomputable def torusCharacter : Multiplicative (Fin 2 → ℤ) →* C(LaurentTorus, ℂ) where
  toFun n := (fourier (n.toAdd 0)).comp ContinuousMap.fst *
    (fourier (n.toAdd 1)).comp ContinuousMap.snd
  map_one' := by
    ext p
    change fourier 0 p.1 * fourier 0 p.2 = 1
    rw [fourier_zero, fourier_zero, mul_one]
  map_mul' n m := by
    ext p
    change fourier (n.toAdd 0 + m.toAdd 0) p.1 * fourier (n.toAdd 1 + m.toAdd 1) p.2 =
      (fourier (n.toAdd 0) p.1 * fourier (n.toAdd 1) p.2) *
      (fourier (m.toAdd 0) p.1 * fourier (m.toAdd 1) p.2)
    rw [fourier_add, fourier_add]
    ring

local notation "Laurent₂" => AddMonoidAlgebra ℂ (Fin 2 → ℤ)

/-- Restrict a complex Laurent polynomial to the compact torus. -/
noncomputable def torusEvaluate : Laurent₂ →ₐ[ℂ] C(LaurentTorus, ℂ) :=
  AddMonoidAlgebra.lift ℂ _ _ torusCharacter

private theorem integral_fourier_circle (n : ℤ) :
    (∫ x, fourier n x ∂AddCircle.haarAddCircle (T := 1)) = if n = 0 then 1 else 0 := by
  have h := congrFun (@fourierCoeff_fourier 1 _ n) 0
  rw [fourierCoeff, Pi.single_apply] at h
  simpa only [neg_zero, fourier_zero, one_smul, eq_comm] using h

/-- The normalized torus integral extracts the constant coefficient. -/
theorem integral_torusEvaluate (f : Laurent₂) :
    (∫ x, torusEvaluate f x ∂laurentTorusHaar) = laurentCT f := by
  classical
  induction f using AddMonoidAlgebra.induction_linear with
  | zero => simp
  | add f g hf hg =>
    rw [map_add]
    change (∫ x, torusEvaluate f x + torusEvaluate g x ∂laurentTorusHaar) = _
    rw [integral_add, hf, hg, map_add]
    all_goals
      apply Continuous.integrable_of_hasCompactSupport
      · exact ContinuousMap.continuous _
      · exact HasCompactSupport.of_compactSpace _
  | single n c =>
    rw [torusEvaluate, AddMonoidAlgebra.lift_single]
    change (∫ x, c * (fourier (n 0) x.1 * fourier (n 1) x.2) ∂laurentTorusHaar) = _
    rw [integral_const_mul, laurentTorusHaar, integral_prod_mul,
      integral_fourier_circle, integral_fourier_circle]
    change c * _ = (Finsupp.single n c) 0
    rw [Finsupp.single_apply]
    have hn : n = 0 ↔ n 0 = 0 ∧ n 1 = 0 := by
      constructor
      · rintro rfl
        exact ⟨rfl, rfl⟩
      · rintro ⟨h₀, h₁⟩
        funext i
        fin_cases i <;> assumption
    by_cases h₀ : n 0 = 0 <;> by_cases h₁ : n 1 = 0 <;> simp [hn, h₀, h₁]

/-- Near zero, the moment generating function equals the torus integral
of the Laurent polynomial divided by its geometric-series denominator. -/
theorem hasSum_torusLaurent_resolvent (f : Laurent₂) (t : ℂ)
    (ht : ‖t‖ * ‖torusEvaluate f‖ < 1) :
    HasSum (fun n : ℕ ↦ t ^ n * laurentCT (f ^ (n + 1)))
      (∫ x, torusEvaluate f x / (1 - t * torusEvaluate f x) ∂laurentTorusHaar) := by
  let F := torusEvaluate f
  let r := ‖t‖ * ‖F‖
  have hr : 0 ≤ r := mul_nonneg (norm_nonneg _) (norm_nonneg _)
  have hb (x : LaurentTorus) : ‖F x‖ ≤ ‖F‖ := ContinuousMap.norm_coe_le_norm _ _
  let a (n : ℕ) (x : LaurentTorus) := F x * (t * F x) ^ n
  have ha (n : ℕ) : Integrable (a n) laurentTorusHaar := by
    apply Continuous.integrable_of_hasCompactSupport
    · dsimp [a]
      fun_prop
    · exact HasCompactSupport.of_compactSpace _
  have hbound (n : ℕ) (x : LaurentTorus) : ‖a n x‖ ≤ ‖F‖ * r ^ n := by
    simp only [a, r, norm_mul, norm_pow]
    gcongr <;> exact hb x
  have hs : Summable (fun n ↦ ∫ x, ‖a n x‖ ∂laurentTorusHaar) := by
    apply ((summable_geometric_of_lt_one hr ht).mul_left ‖F‖).of_nonneg_of_le
    · exact fun n ↦ integral_nonneg (fun _ ↦ norm_nonneg _)
    · intro n
      calc
        _ ≤ ∫ _ : LaurentTorus, ‖F‖ * r ^ n ∂laurentTorusHaar :=
          integral_mono (ha n).norm (integrable_const _) (hbound n)
        _ = _ := by simp
  have hseries := hasSum_integral_of_summable_integral_norm ha hs
  have hsum (x : LaurentTorus) : (∑' n, a n x) = F x / (1 - t * F x) := by
    have htx : ‖t * F x‖ < 1 := by
      rw [norm_mul]
      exact (mul_le_mul_of_nonneg_left (hb x) (norm_nonneg _)).trans_lt ht
    simpa only [a, div_eq_mul_inv] using ((hasSum_geometric_of_norm_lt_one htx).mul_left (F x)).tsum_eq
  have hcoeff (n : ℕ) :
      (∫ x, a n x ∂laurentTorusHaar) = t ^ n * laurentCT (f ^ (n + 1)) := by
    have hx (x : LaurentTorus) : a n x = t ^ n * torusEvaluate (f ^ (n + 1)) x := by
      rw [map_pow, ContinuousMap.pow_apply, pow_succ]
      simp only [a, F, mul_pow]
      ring
    simp only [hx, integral_const_mul, integral_torusEvaluate]
  simpa only [hcoeff, hsum, F] using hseries

/-- Purity makes the resolvent integral vanish throughout its initial
disk of convergence. -/
theorem integral_torusLaurent_resolvent_eq_zero (f : Laurent₂) (hp : Pure laurentCT f)
    (t : ℂ) (ht : ‖t‖ * ‖torusEvaluate f‖ < 1) :
    (∫ x, torusEvaluate f x / (1 - t * torusEvaluate f x) ∂laurentTorusHaar) = 0 := by
  have h := (hasSum_torusLaurent_resolvent f t ht).tsum_eq
  rw [← h]
  simp only [hp _ (Nat.succ_pos _), mul_zero, tsum_zero]

end MathieuQuadratic
