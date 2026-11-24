import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Algebra.BigOperators.Field

open scoped BigOperators InnerProductSpace Classical

namespace QuadraticNorms

private theorem complex_star_mul_self (z : ℂ) : star z * z = ((‖z‖ ^ 2 : ℝ) : ℂ) := by
  change (starRingEnd ℂ) z * z = _
  rw [mul_comm, Complex.mul_conj, Complex.normSq_eq_norm_sq]

section Frame

variable {ι E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]

/-- Synthesis after analysis in the finite character frame. -/
noncomputable def bankFrameOperator (s : Finset ι) (v : ι → E) (x : E) : E :=
  ∑ i ∈ s, ⟪v i, x⟫_ℂ • v i

private theorem bank_frame_synthesis_sq_bound (s : Finset ι) (v : ι → E) {δ : ℝ} (hδ : 0 ≤ δ)
    (hdiag : ∀ i ∈ s, ‖v i‖ = 1)
    (hoff : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → ‖⟪v i, v j⟫_ℂ‖ ≤ δ) (a : ι → ℂ) :
    ‖∑ i ∈ s, a i • v i‖ ^ 2 ≤ (1 + (s.card : ℝ) * δ) * ∑ i ∈ s, ‖a i‖ ^ 2 := by
  let y := ∑ i ∈ s, a i • v i
  let b := ∑ i ∈ s, ‖a i‖ ^ 2
  let d (i j : ι) : ℂ := if i = j then 1 else 0
  let G (i j : ι) := ⟪v i, v j⟫_ℂ - d i j
  have hG (i : ι) (hi : i ∈ s) (j : ι) (hj : j ∈ s) : ‖G i j‖ ≤ δ := by
    by_cases hh : i = j
    · subst j
      simpa [G, d, inner_self_eq_norm_sq_to_K, hdiag i hi] using hδ
    · simpa only [G, d, if_neg hh, sub_zero] using hoff i hi j hj hh
  have hd : (∑ i ∈ s, ∑ j ∈ s, star (a i) * a j * d i j) = (b : ℂ) := by
    rw [show b = ∑ i ∈ s, ‖a i‖ ^ 2 from rfl, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    simp only [d, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq, if_pos hi]
    exact complex_star_mul_self _
  have hinner : ((‖y‖ ^ 2 : ℝ) : ℂ) =
      (b : ℂ) + ∑ i ∈ s, ∑ j ∈ s, star (a i) * a j * G i j := by
    have he : ⟪y, y⟫_ℂ = ((‖y‖ ^ 2 : ℝ) : ℂ) :=
      (inner_self_eq_norm_sq_to_K (𝕜 := ℂ) y).trans (Complex.ofReal_pow _ _).symm
    rw [← he]
    simp only [y, sum_inner, inner_sum, inner_smul_left, inner_smul_right, starRingEnd_apply,
      Finset.mul_sum]
    rw [Finset.sum_comm, ← hd, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl (fun j _ => by dsimp [G]; ring)
  have hb : 0 ≤ b := Finset.sum_nonneg (fun _ _ => sq_nonneg _)
  have herr : ‖∑ i ∈ s, ∑ j ∈ s, star (a i) * a j * G i j‖ ≤ δ * (∑ i ∈ s, ‖a i‖) ^ 2 := by
    refine (norm_sum_le _ _).trans ?_
    calc
      _ ≤ ∑ i ∈ s, ∑ j ∈ s, ‖a i‖ * ‖a j‖ * δ :=
        Finset.sum_le_sum (fun i hi => (norm_sum_le _ _).trans
          (Finset.sum_le_sum (fun j hj => by
            rw [norm_mul, norm_mul, norm_star]
            exact mul_le_mul_of_nonneg_left (hG i hi j hj) (by positivity))))
      _ = ∑ i ∈ s, (‖a i‖ * (∑ j ∈ s, ‖a j‖)) * δ := by
        simp only [Finset.mul_sum, Finset.sum_mul]
      _ = _ := by
        rw [← Finset.sum_mul, ← Finset.sum_mul, pow_two]
        ring
  have hsum : (∑ i ∈ s, ‖a i‖) ^ 2 ≤ (s.card : ℝ) * b := by
    have h := Finset.sum_mul_sq_le_sq_mul_sq s (fun i => ‖a i‖) (fun _ => (1 : ℝ))
    simpa only [mul_one, one_pow, Finset.sum_const, nsmul_eq_mul, mul_comm, b] using h
  calc
    ‖y‖ ^ 2 = ‖((‖y‖ ^ 2 : ℝ) : ℂ)‖ := by rw [Complex.norm_real, Real.norm_of_nonneg (by positivity)]
    _ ≤ b + δ * (∑ i ∈ s, ‖a i‖) ^ 2 := by
      rw [hinner]
      refine (norm_add_le ..).trans (add_le_add ?_ herr)
      exact le_of_eq (by rw [Complex.norm_real, Real.norm_of_nonneg hb])
    _ ≤ b + δ * ((s.card : ℝ) * b) := add_le_add le_rfl (mul_le_mul_of_nonneg_left hsum hδ)
    _ = _ := by ring

/-- A frame with unit diagonal and off-diagonal Gram entries bounded
by δ has operator norm at most 1 plus δ times the frame size. -/
theorem bank_frame_operator_norm_le (s : Finset ι) (v : ι → E) {δ : ℝ} (hδ : 0 ≤ δ)
    (hdiag : ∀ i ∈ s, ‖v i‖ = 1)
    (hoff : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → ‖⟪v i, v j⟫_ℂ‖ ≤ δ) (x : E) :
    ‖bankFrameOperator s v x‖ ≤ (1 + (s.card : ℝ) * δ) * ‖x‖ := by
  let a i := ⟪v i, x⟫_ℂ
  let b := ∑ i ∈ s, ‖a i‖ ^ 2
  let z := bankFrameOperator s v x
  let C := 1 + (s.card : ℝ) * δ
  have hC : 0 ≤ C := by dsimp [C]; positivity
  have hb : 0 ≤ b := Finset.sum_nonneg (fun _ _ => sq_nonneg _)
  have hinner : ⟪z, x⟫_ℂ = (b : ℂ) := by
    rw [show b = ∑ i ∈ s, ‖a i‖ ^ 2 from rfl, Complex.ofReal_sum]
    simp only [z, bankFrameOperator, sum_inner, inner_smul_left]
    exact Finset.sum_congr rfl (fun i _ => complex_star_mul_self (a i))
  have hSynth : ‖z‖ ^ 2 ≤ C * b := bank_frame_synthesis_sq_bound s v hδ hdiag hoff a
  have hNorm : b ≤ ‖z‖ * ‖x‖ := by
    rw [← Real.norm_of_nonneg hb, ← Complex.norm_real, ← hinner]
    exact norm_inner_le_norm z x
  have hAnalysis : b ≤ C * ‖x‖ ^ 2 := by
    by_cases hb₀ : b = 0
    · rw [hb₀]
      positivity
    · have hh : b ^ 2 ≤ (‖z‖ * ‖x‖) ^ 2 := pow_le_pow_left₀ hb hNorm 2
      rw [mul_pow] at hh
      have hbound := hh.trans (mul_le_mul_of_nonneg_right hSynth (sq_nonneg ‖x‖))
      apply (mul_le_mul_iff_of_pos_left (lt_of_le_of_ne hb (Ne.symm hb₀))).mp
      nlinarith only [hbound]
  apply le_of_sq_le_sq _ (by positivity)
  exact (hSynth.trans (mul_le_mul_of_nonneg_left hAnalysis hC)).trans_eq (by ring)

end Frame

section WeightedCell

variable {X ι : Type*} [Fintype X]

/-- Put weighted functions in ordinary Euclidean space by multiplying
by the square roots of the weights. -/
private noncomputable def bankCellVector (w : X → ℝ) (f : X → ℂ) : EuclideanSpace ℂ X :=
  WithLp.toLp 2 (fun x => (Real.sqrt (w x) : ℂ) * f x)

private theorem bank_cell_vector_norm_sq (w : X → ℝ) (hw : ∀ x, 0 ≤ w x) (f : X → ℂ) :
    ‖bankCellVector w f‖ ^ 2 = ∑ x, w x * ‖f x‖ ^ 2 := by
  rw [EuclideanSpace.norm_sq_eq]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  change ‖(Real.sqrt (w x) : ℂ) * f x‖ ^ 2 = _
  rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg (Real.sqrt_nonneg _),
    mul_pow, Real.sq_sqrt (hw x)]

/-- The Gram entry for functions on a cell with its conditional weights. -/
noncomputable def bankCellInner (w : X → ℝ) (f g : X → ℂ) : ℂ :=
  ∑ x, (w x : ℂ) * (star (f x) * g x)

private theorem bank_cell_vector_inner (w : X → ℝ) (hw : ∀ x, 0 ≤ w x) (f g : X → ℂ) :
    ⟪bankCellVector w f, bankCellVector w g⟫_ℂ = bankCellInner w f g := by
  rw [PiLp.inner_apply, bankCellInner]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  simp only [bankCellVector, PiLp.toLp_apply, RCLike.inner_apply, map_mul, Complex.conj_ofReal]
  change (Real.sqrt (w x) : ℂ) * g x * ((Real.sqrt (w x) : ℂ) * star (f x)) =
    (w x : ℂ) * (star (f x) * g x)
  have hs : (Real.sqrt (w x) : ℂ) ^ 2 = (w x : ℂ) := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt (hw x)]
  calc
    _ = (Real.sqrt (w x) : ℂ) ^ 2 * (star (f x) * g x) := by ring
    _ = _ := by rw [hs]

/-- The character kernel on one cell with conditional probability weights. -/
noncomputable def bankCellKernel (s : Finset ι) (χ : ι → X → ℂ) (d x : X) : ℂ :=
  ∑ i ∈ s, star (χ i d) * χ i x

noncomputable def bankCellApply (w : X → ℝ) (s : Finset ι) (χ : ι → X → ℂ) (f : X → ℂ) (x : X) : ℂ :=
  ∑ d, (w d : ℂ) * bankCellKernel s χ d x * f d

private theorem bank_cell_apply_eq (w : X → ℝ) (s : Finset ι) (χ : ι → X → ℂ)
    (f : X → ℂ) (x : X) :
    bankCellApply w s χ f x = ∑ i ∈ s, bankCellInner w (χ i) f * χ i x := by
  unfold bankCellApply bankCellKernel bankCellInner
  simp only [Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun d _ => ?_))
  ring

/-- The character kernel on a retained cell has the asserted L²
energy bound from its off-diagonal weighted Gram estimate. -/
theorem bank_cell_apply_energy_le (w : X → ℝ) (hw : ∀ x, 0 ≤ w x) (hprob : ∑ x, w x = 1)
    (s : Finset ι) (χ : ι → X → ℂ) (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1)
    {δ : ℝ} (hδ : 0 ≤ δ) (hoff : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → ‖bankCellInner w (χ i) (χ j)‖ ≤ δ)
    (f : X → ℂ) :
    (∑ x, w x * ‖bankCellApply w s χ f x‖ ^ 2) ≤
      (1 + (s.card : ℝ) * δ) ^ 2 * ∑ x, w x * ‖f x‖ ^ 2 := by
  let v i := bankCellVector w (χ i)
  have hdiag (i : ι) (hi : i ∈ s) : ‖v i‖ = 1 := by
    have h : ‖v i‖ ^ 2 = (1 : ℝ) ^ 2 := by
      rw [bank_cell_vector_norm_sq w hw]
      simp only [hunit i hi, one_pow, mul_one, hprob]
    exact (pow_left_inj₀ (norm_nonneg _) (by norm_num) (by decide : (2 : ℕ) ≠ 0)).mp h
  have hcross (i : ι) (hi : i ∈ s) (j : ι) (hj : j ∈ s) (hne : i ≠ j) : ‖⟪v i, v j⟫_ℂ‖ ≤ δ := by
    rw [bank_cell_vector_inner w hw]
    exact hoff i hi j hj hne
  have hvec : bankCellVector w (bankCellApply w s χ f) =
      bankFrameOperator s v (bankCellVector w f) := by
    apply WithLp.ofLp_injective 2
    funext x
    simp only [bankFrameOperator, WithLp.ofLp_sum, Finset.sum_apply, WithLp.ofLp_smul,
      Pi.smul_apply, smul_eq_mul]
    simp only [v, bank_cell_vector_inner w hw]
    simp only [bankCellVector, WithLp.ofLp_toLp, bank_cell_apply_eq, Finset.mul_sum]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  rw [← bank_cell_vector_norm_sq w hw, ← bank_cell_vector_norm_sq w hw, hvec]
  have hn := bank_frame_operator_norm_le s v hδ hdiag hcross (bankCellVector w f)
  exact (pow_le_pow_left₀ (norm_nonneg _) hn 2).trans_eq (mul_pow ..)

private theorem bank_cell_inner_mul_right (w : X → ℝ) (f g : X → ℂ) (ζ : ℂ) :
    bankCellInner w f (fun x => ζ * g x) = ζ * bankCellInner w f g := by
  rw [bankCellInner, bankCellInner, Finset.mul_sum]
  exact Finset.sum_congr rfl (fun _ _ => by ring)

private theorem bank_cell_inner_unit_self (w : X → ℝ) (hprob : ∑ x, w x = 1) (f : X → ℂ)
    (hf : ∀ x, ‖f x‖ = 1) : bankCellInner w f f = 1 := by
  unfold bankCellInner
  simp only [complex_star_mul_self, hf, one_pow, Complex.ofReal_one, mul_one, ← Complex.ofReal_sum, hprob]

/-- The kernel reproduces the phase of a frame representative,
with error at most the frame size times the off-diagonal bound. -/
theorem bank_cell_reproduces_phase (w : X → ℝ) (hprob : ∑ x, w x = 1)
    (s : Finset ι) (χ : ι → X → ℂ) (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1)
    {δ : ℝ} (hδ : 0 ≤ δ) {j : ι} (hj : j ∈ s)
    (hoff : ∀ i ∈ s, i ≠ j → ‖bankCellInner w (χ i) (χ j)‖ ≤ δ)
    {ζ : ℂ} (hζ : ‖ζ‖ = 1) (x : X) :
    ‖bankCellApply w s χ (fun d => ζ * χ j d) x - ζ * χ j x‖ ≤ (s.card : ℝ) * δ := by
  have hsole : (∑ i ∈ s, if i = j then ζ * χ i x else 0) = ζ * χ j x := by
    rw [Finset.sum_ite_eq', if_pos hj]
  have he : bankCellApply w s χ (fun d => ζ * χ j d) x - ζ * χ j x =
      ∑ i ∈ s, ζ * (bankCellInner w (χ i) (χ j) - (if i = j then 1 else 0)) * χ i x := by
    rw [bank_cell_apply_eq, ← hsole, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [bank_cell_inner_mul_right]
    split <;> ring
  rw [he]
  refine (norm_sum_le _ _).trans ?_
  calc
    _ ≤ ∑ _i ∈ s, δ := by
      refine Finset.sum_le_sum (fun i hi => ?_)
      rw [norm_mul, norm_mul, hζ, hunit i hi x, one_mul, mul_one]
      by_cases hij : i = j
      · subst i
        rw [if_pos rfl, bank_cell_inner_unit_self w hprob _ (hunit j hj), sub_self, norm_zero]
        exact hδ
      · rw [if_neg hij, sub_zero]
        exact hoff i hi hij
    _ = _ := by rw [Finset.sum_const, nsmul_eq_mul]

omit [Fintype X] in
/-- The kernel bound before dividing by a retained cell mass. -/
theorem bank_cell_kernel_norm_le (s : Finset ι) (χ : ι → X → ℂ) (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1)
    (d x : X) : ‖bankCellKernel s χ d x‖ ≤ s.card := by
  unfold bankCellKernel
  refine (norm_sum_le _ _).trans_eq ?_
  simp_rw [norm_mul, norm_star]
  calc
    _ = ∑ _i ∈ s, (1 : ℝ) := Finset.sum_congr rfl (fun i hi => by rw [hunit i hi d, hunit i hi x, one_mul])
    _ = _ := by rw [Finset.sum_const, nsmul_eq_mul, mul_one]

end WeightedCell

section RetainedCells

variable {X ι J : Type*} [Fintype X] [DecidableEq J]

/-- The mass of a logarithmic product cell. -/
noncomputable def bankLogCellMass (w : X → ℝ) (cell : X → J) (j : J) : ℝ :=
  ∑ d ∈ Finset.univ.filter (fun d => cell d = j), w d

private noncomputable def bankConditionalWeight (w : X → ℝ) (cell : X → J) (j : J) (d : X) : ℝ :=
  if cell d = j then w d / bankLogCellMass w cell j else 0

private theorem bank_conditional_weight_nonneg (w : X → ℝ) (hw : ∀ d, 0 ≤ w d)
    (cell : X → J) (j : J) (hν : 0 < bankLogCellMass w cell j) (d : X) :
    0 ≤ bankConditionalWeight w cell j d := by
  unfold bankConditionalWeight
  split
  · exact div_nonneg (hw d) hν.le
  · exact le_rfl

private theorem bank_conditional_weight_sum (w : X → ℝ) (cell : X → J) (j : J)
    (hν : 0 < bankLogCellMass w cell j) :
    ∑ d, bankConditionalWeight w cell j d = 1 := by
  unfold bankConditionalWeight
  rw [← Finset.sum_filter, ← Finset.sum_div]
  exact div_self hν.ne'

private theorem bank_conditional_weight_inner (w : X → ℝ) (cell : X → J) (j : J) (f g : X → ℂ) :
    bankCellInner (bankConditionalWeight w cell j) f g =
      (∑ d ∈ Finset.univ.filter (fun d => cell d = j), (w d : ℂ) * (star (f d) * g d)) /
        (bankLogCellMass w cell j : ℂ) := by
  rw [bankCellInner, Finset.sum_div, Finset.sum_filter]
  refine Finset.sum_congr rfl (fun d _ => ?_)
  by_cases hd : cell d = j
  · rw [bankConditionalWeight, if_pos hd, if_pos hd, Complex.ofReal_div]
    ring
  · simp [bankConditionalWeight, hd]

/-- The progression kernel, with low-mass cells discarded. -/
noncomputable def bankProgressionKernel (w : X → ℝ) (cell : X → J) (ξ : ℝ)
    (s : Finset ι) (χ : ι → X → ℂ) (d x : X) : ℂ :=
  if cell d = cell x ∧ ξ ≤ bankLogCellMass w cell (cell x) then
    bankCellKernel s χ d x / (bankLogCellMass w cell (cell x) : ℂ)
  else 0

noncomputable def bankProgressionApply (w : X → ℝ) (cell : X → J) (ξ : ℝ)
    (s : Finset ι) (χ : ι → X → ℂ) (f : X → ℂ) (x : X) : ℂ :=
  ∑ d, (w d : ℂ) * bankProgressionKernel w cell ξ s χ d x * f d

private theorem bank_progression_apply_retained (w : X → ℝ) (cell : X → J) (ξ : ℝ)
    (s : Finset ι) (χ : ι → X → ℂ) (f : X → ℂ) (x : X) (hx : ξ ≤ bankLogCellMass w cell (cell x)) :
    bankProgressionApply w cell ξ s χ f x = bankCellApply (bankConditionalWeight w cell (cell x)) s χ f x := by
  rw [bankProgressionApply, bankCellApply]
  refine Finset.sum_congr rfl (fun d _ => ?_)
  by_cases hd : cell d = cell x
  · simp only [bankProgressionKernel, bankConditionalWeight, hd, hx, and_self, if_true, Complex.ofReal_div]
    ring
  · simp [bankProgressionKernel, bankConditionalWeight, hd]

theorem bank_progression_kernel_norm_le (w : X → ℝ) (cell : X → J) {ξ : ℝ} (hξ : 0 < ξ)
    (s : Finset ι) (χ : ι → X → ℂ) (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1) (d x : X) :
    ‖bankProgressionKernel w cell ξ s χ d x‖ ≤ (s.card : ℝ) / ξ := by
  unfold bankProgressionKernel
  split
  · rename_i h
    have hν : 0 < bankLogCellMass w cell (cell x) := hξ.trans_le h.2
    rw [norm_div, Complex.norm_real, Real.norm_of_nonneg hν.le]
    exact (div_le_div_of_nonneg_right (bank_cell_kernel_norm_le s χ hunit d x) hν.le).trans
      (div_le_div_of_nonneg_left (by positivity) hξ h.2)
  · rw [norm_zero]
    exact div_nonneg (by positivity) hξ.le

theorem bank_progression_kernel_mean_norm_le (w : X → ℝ) (hw : ∀ d, 0 ≤ w d)
    (cell : X → J) {ξ : ℝ} (hξ : 0 < ξ) (s : Finset ι) (χ : ι → X → ℂ)
    (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1) (x : X) :
    (∑ d, w d * ‖bankProgressionKernel w cell ξ s χ d x‖) ≤ s.card := by
  let ν := bankLogCellMass w cell (cell x)
  by_cases hx : ξ ≤ ν
  · have hν : 0 < ν := hξ.trans_le hx
    have hpoint (d : X) :
        w d * ‖bankProgressionKernel w cell ξ s χ d x‖ ≤
          if cell d = cell x then w d * ((s.card : ℝ) / ν) else 0 := by
      by_cases hd : cell d = cell x
      · have hret : cell d = cell x ∧ ξ ≤ bankLogCellMass w cell (cell x) := ⟨hd, hx⟩
        rw [bankProgressionKernel, if_pos hret, if_pos hd, norm_div, Complex.norm_real, Real.norm_of_nonneg hν.le]
        exact mul_le_mul_of_nonneg_left
          (div_le_div_of_nonneg_right (bank_cell_kernel_norm_le s χ hunit d x) hν.le) (hw d)
      · simp [bankProgressionKernel, hd]
    refine (Finset.sum_le_sum (fun d _ => hpoint d)).trans_eq ?_
    rw [← Finset.sum_filter, ← Finset.sum_mul]
    change ν * ((s.card : ℝ) / ν) = _
    exact mul_div_cancel₀ _ hν.ne'
  · have hn (d : X) : bankProgressionKernel w cell ξ s χ d x = 0 := by
      unfold bankProgressionKernel
      exact if_neg (fun h => hx h.2)
    simp [hn]

private theorem bank_conditional_gram_bound (w : X → ℝ) (cell : X → J) {ξ η : ℝ} (hξ : 0 < ξ)
    (hη : 0 ≤ η) {c : J} (hc : ξ ≤ bankLogCellMass w cell c) {f g : X → ℂ}
    (hb : ‖∑ d ∈ Finset.univ.filter (fun d => cell d = c), (w d : ℂ) * (star (f d) * g d)‖ ≤ η) :
    ‖bankCellInner (bankConditionalWeight w cell c) f g‖ ≤ η / ξ := by
  have hν := hξ.trans_le hc
  rw [bank_conditional_weight_inner, norm_div, Complex.norm_real, Real.norm_of_nonneg hν.le]
  exact (div_le_div_of_nonneg_right hb hν.le).trans (div_le_div_of_nonneg_left hη hξ hc)

/-- The retained-cell kernel reproduces a character branch up to the
weighted off-diagonal error divided by the cell-mass threshold. -/
theorem bank_progression_reproduces_phase (w : X → ℝ) (cell : X → J) {ξ η : ℝ}
    (hξ : 0 < ξ) (hη : 0 ≤ η) (s : Finset ι) (χ : ι → X → ℂ)
    (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1) (x : X) (hx : ξ ≤ bankLogCellMass w cell (cell x))
    {j : ι} (hj : j ∈ s) {ζ : ℂ} (hζ : ‖ζ‖ = 1)
    (hoff : ∀ i ∈ s, i ≠ j →
      ‖∑ d ∈ Finset.univ.filter (fun d => cell d = cell x), (w d : ℂ) * (star (χ i d) * χ j d)‖ ≤ η) :
    ‖bankProgressionApply w cell ξ s χ (fun d => ζ * χ j d) x - ζ * χ j x‖ ≤
      (s.card : ℝ) * (η / ξ) := by
  rw [bank_progression_apply_retained w cell ξ s χ _ x hx]
  have hprob := bank_conditional_weight_sum w cell (cell x) (hξ.trans_le hx)
  exact bank_cell_reproduces_phase _ hprob s χ hunit (by positivity) hj
    (fun i hi hij => bank_conditional_gram_bound w cell hξ hη hx (hoff i hi hij)) hζ x

private theorem bank_conditional_summation (w : X → ℝ) (cell : X → J) (j : J)
    (hν : 0 < bankLogCellMass w cell j) (g : X → ℝ) :
    (∑ x ∈ Finset.univ.filter (fun x => cell x = j), w x * g x) =
      bankLogCellMass w cell j * ∑ x, bankConditionalWeight w cell j x * g x := by
  rw [Finset.sum_filter, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  by_cases hx : cell x = j
  · simp only [if_pos hx, bankConditionalWeight]
    field_simp
  · simp [bankConditionalWeight, hx]

/-- The full kernel energy bound is obtained by summing the conditional
energy bounds on the orthogonal retained cells. -/
theorem bank_progression_apply_energy_le (w : X → ℝ) (hw : ∀ d, 0 ≤ w d)
    (cell : X → J) {ξ η : ℝ} (hξ : 0 < ξ) (hη : 0 ≤ η)
    (s : Finset ι) (χ : ι → X → ℂ) (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1)
    (hoff : ∀ c, ξ ≤ bankLogCellMass w cell c → ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      ‖∑ d ∈ Finset.univ.filter (fun d => cell d = c), (w d : ℂ) * (star (χ i d) * χ j d)‖ ≤ η)
    (f : X → ℂ) :
    (∑ x, w x * ‖bankProgressionApply w cell ξ s χ f x‖ ^ 2) ≤
      (1 + (s.card : ℝ) * (η / ξ)) ^ 2 * ∑ x, w x * ‖f x‖ ^ 2 := by
  let C := (1 + (s.card : ℝ) * (η / ξ)) ^ 2
  let cells := Finset.univ.image cell
  let A := bankProgressionApply w cell ξ s χ f
  have hmap (x : X) (_hx : x ∈ Finset.univ) : cell x ∈ cells := Finset.mem_image_of_mem cell (Finset.mem_univ x)
  have hcell (c : J) :
      (∑ x ∈ Finset.univ.filter (fun x => cell x = c), w x * ‖A x‖ ^ 2) ≤
        C * ∑ x ∈ Finset.univ.filter (fun x => cell x = c), w x * ‖f x‖ ^ 2 := by
    by_cases hc : ξ ≤ bankLogCellMass w cell c
    · have hν := hξ.trans_le hc
      let W := bankConditionalWeight w cell c
      have hW := bank_conditional_weight_nonneg w hw cell c hν
      have hprob := bank_conditional_weight_sum w cell c hν
      have hlocal := bank_cell_apply_energy_le W hW hprob s χ hunit (div_nonneg hη hξ.le)
        (fun i hi j hj hij => bank_conditional_gram_bound w cell hξ hη hc (hoff c hc i hi j hj hij)) f
      have heq : (∑ x, W x * ‖A x‖ ^ 2) = ∑ x, W x * ‖bankCellApply W s χ f x‖ ^ 2 := by
        refine Finset.sum_congr rfl (fun x _ => ?_)
        by_cases hx : cell x = c
        · have hret : ξ ≤ bankLogCellMass w cell (cell x) := by rwa [hx]
          have happly := bank_progression_apply_retained w cell ξ s χ f x hret
          rw [hx] at happly
          change W x * ‖bankProgressionApply w cell ξ s χ f x‖ ^ 2 = _
          rw [happly]
        · simp [W, bankConditionalWeight, hx]
      rw [bank_conditional_summation w cell c hν, bank_conditional_summation w cell c hν,
        show (∑ x, bankConditionalWeight w cell c x * ‖A x‖ ^ 2) =
          ∑ x, W x * ‖bankCellApply W s χ f x‖ ^ 2 from heq]
      exact (mul_le_mul_of_nonneg_left hlocal hν.le).trans_eq (by unfold C; ring)
    · have hz : (∑ x ∈ Finset.univ.filter (fun x => cell x = c), w x * ‖A x‖ ^ 2) = 0 := by
        refine Finset.sum_eq_zero (fun x hx => ?_)
        have hx' := (Finset.mem_filter.mp hx).2
        have hA : A x = 0 := by
          unfold A bankProgressionApply
          refine Finset.sum_eq_zero (fun d _ => ?_)
          have hret : ¬ (cell d = cell x ∧ ξ ≤ bankLogCellMass w cell (cell x)) :=
            fun h => hc (hx' ▸ h.2)
          rw [bankProgressionKernel, if_neg hret, mul_zero, zero_mul]
        simp [hA]
      rw [hz]
      exact mul_nonneg (sq_nonneg _) (Finset.sum_nonneg (fun x _ => mul_nonneg (hw x) (sq_nonneg _)))
  calc
    _ = ∑ c ∈ cells, ∑ x ∈ Finset.univ.filter (fun x => cell x = c), w x * ‖A x‖ ^ 2 :=
      (Finset.sum_fiberwise_of_maps_to hmap (fun x => w x * ‖A x‖ ^ 2)).symm
    _ ≤ ∑ c ∈ cells, C * ∑ x ∈ Finset.univ.filter (fun x => cell x = c), w x * ‖f x‖ ^ 2 :=
      Finset.sum_le_sum (fun c _ => hcell c)
    _ = _ := by rw [← Finset.mul_sum, Finset.sum_fiberwise_of_maps_to hmap]

/-- In the accuracy range used by the graph comparison, the retained
frame operator has norm at most two, hence energy cost at most four. -/
theorem bank_progression_frame_energy_le_four (w : X → ℝ) (hw : ∀ d, 0 ≤ w d)
    (cell : X → J) {ξ η : ℝ} (hξ : 0 < ξ) (hη : 0 ≤ η)
    (s : Finset ι) (χ : ι → X → ℂ) (hunit : ∀ i ∈ s, ∀ x, ‖χ i x‖ = 1)
    (hacc : (s.card : ℝ) * (η / ξ) ≤ 1)
    (hoff : ∀ c, ξ ≤ bankLogCellMass w cell c → ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      ‖∑ d ∈ Finset.univ.filter (fun d => cell d = c), (w d : ℂ) * (star (χ i d) * χ j d)‖ ≤ η)
    (f : X → ℂ) :
    (∑ x, w x * ‖bankProgressionApply w cell ξ s χ f x‖ ^ 2) ≤
      4 * ∑ x, w x * ‖f x‖ ^ 2 := by
  refine (bank_progression_apply_energy_le w hw cell hξ hη s χ hunit hoff f).trans ?_
  apply mul_le_mul_of_nonneg_right _ (Finset.sum_nonneg (fun x _ => mul_nonneg (hw x) (sq_nonneg _)))
  have h := pow_le_pow_left₀ (by positivity : 0 ≤ 1 + (s.card : ℝ) * (η / ξ))
    (by linarith only [hacc] : 1 + (s.card : ℝ) * (η / ξ) ≤ 2) 2
  norm_num at h
  exact h

end RetainedCells

end QuadraticNorms
