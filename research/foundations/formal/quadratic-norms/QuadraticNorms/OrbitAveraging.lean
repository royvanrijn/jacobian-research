import Mathlib.Analysis.Real.Sqrt
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Module.Pi

open scoped BigOperators

namespace QuadraticNorms

variable {Y : Type*} [Fintype Y]

/-- The real inner product weighted by a probability mass function on a finite space. -/
noncomputable def weightedPair (π : Y → ℝ) (f g : Y → ℝ) : ℝ :=
  ∑ y, π y * (f y * g y)

private theorem weightedPair_comm (π f g : Y → ℝ) :
    weightedPair π f g = weightedPair π g f := by simp [weightedPair, mul_comm (f _) (g _)]

private theorem weightedPair_add_left (π f g h : Y → ℝ) :
    weightedPair π (f + g) h = weightedPair π f h + weightedPair π g h := by
  simp only [weightedPair, Pi.add_apply, add_mul, mul_add, Finset.sum_add_distrib]

private theorem weightedPair_add_right (π f g h : Y → ℝ) :
    weightedPair π f (g + h) = weightedPair π f g + weightedPair π f h := by
  simp only [weightedPair, Pi.add_apply, mul_add, Finset.sum_add_distrib]

theorem weightedPair_smul_left (π f g : Y → ℝ) (r : ℝ) :
    weightedPair π (r • f) g = r * weightedPair π f g := by
  rw [weightedPair, weightedPair, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro y _
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

private theorem weightedPair_smul_right (π f g : Y → ℝ) (r : ℝ) :
    weightedPair π f (r • g) = r * weightedPair π f g := by
  rw [weightedPair_comm π f, weightedPair_smul_left, weightedPair_comm π g f]

private theorem weightedPair_self_nonneg (π f : Y → ℝ) (hπ : ∀ y, 0 ≤ π y) :
    0 ≤ weightedPair π f f := by
  apply Finset.sum_nonneg
  exact fun y _ => mul_nonneg (hπ y) (mul_self_nonneg (f y))

theorem weightedPair_cauchySchwarz (π f g : Y → ℝ) (hπ : ∀ y, 0 ≤ π y) :
    weightedPair π f g ^ 2 ≤ weightedPair π f f * weightedPair π g g := by
  let wf := fun y => Real.sqrt (π y) * f y
  let wg := fun y => Real.sqrt (π y) * g y
  have hprod (a b : Y → ℝ) (y : Y) :
      (Real.sqrt (π y) * a y) * (Real.sqrt (π y) * b y) = π y * (a y * b y) := by
    rw [show (Real.sqrt (π y) * a y) * (Real.sqrt (π y) * b y) =
      Real.sqrt (π y) ^ 2 * (a y * b y) by ring, Real.sq_sqrt (hπ y)]
  simpa only [weightedPair, wf, wg, hprod, pow_two] using
    (Finset.sum_mul_sq_le_sq_mul_sq Finset.univ wf wg)

private theorem weightedPair_abs_le_of_norm_le {π f g : Y → ℝ} (hπ : ∀ y, 0 ≤ π y)
    {Λ : ℝ} (hΛ : 0 ≤ Λ)
    (hspectral : weightedPair π g g ≤ Λ ^ 2 * weightedPair π f f) :
    |weightedPair π f g| ≤ Λ * weightedPair π f f := by
  have hnorm := weightedPair_self_nonneg π f hπ
  have hsq := (weightedPair_cauchySchwarz π f g hπ).trans
    (mul_le_mul_of_nonneg_left hspectral hnorm)
  have hsq' : |weightedPair π f g| ^ 2 ≤ (Λ * weightedPair π f f) ^ 2 := by
    rw [sq_abs, mul_pow]
    exact hsq.trans_eq (by ring)
  exact le_of_sq_le_sq hsq' (mul_nonneg hΛ hnorm)

/-- The indicator of the selected orbit. -/
noncomputable def orbitIndicator [DecidableEq Y] (Ω : Finset Y) : Y → ℝ :=
  fun y => if y ∈ Ω then 1 else 0

theorem weightedPair_indicator_left [DecidableEq Y] (π g : Y → ℝ) (Ω : Finset Y) :
    weightedPair π (orbitIndicator Ω) g = ∑ y ∈ Ω, π y * g y := by
  unfold weightedPair orbitIndicator
  have hfilter : Finset.univ.filter (fun y => y ∈ Ω) = Ω := by ext y; simp
  rw [← hfilter, Finset.sum_filter]
  apply Finset.sum_congr rfl
  intro y _
  by_cases hy : y ∈ Ω <;> simp [hy]


/-- The spectral orbit estimate, using constancy of the image on the orbit.
For equivariant operators this constancy follows from invariance of its indicator. -/
theorem spectral_orbit_deviation [DecidableEq Y] {π : Y → ℝ}
    (hπ : ∀ y, 0 ≤ π y) (hprob : ∑ y, π y = 1)
    (Ω : Finset Y) (y₀ : Y) (T : (Y → ℝ) →ₗ[ℝ] (Y → ℝ))
    {d Λ : ℝ} (hd : 0 < d) (hΛ : 0 ≤ Λ)
    (hμ : 0 < ∑ y ∈ Ω, π y)
    (hunit : T 1 = d • 1)
    (hself : ∀ f g, weightedPair π (T f) g = weightedPair π f (T g))
    (hconstant : ∀ y ∈ Ω, T (orbitIndicator Ω) y = T (orbitIndicator Ω) y₀)
    (hspectral :
      let h := orbitIndicator Ω - (∑ y ∈ Ω, π y) • (1 : Y → ℝ)
      weightedPair π (T h) (T h) ≤ Λ ^ 2 * weightedPair π h h) :
    |T (orbitIndicator Ω) y₀ / d - ∑ y ∈ Ω, π y| ≤
      Λ / d * (1 - ∑ y ∈ Ω, π y) := by
  let μ := ∑ y ∈ Ω, π y
  let f := orbitIndicator Ω
  let h := f - μ • (1 : Y → ℝ)
  change 0 < μ at hμ
  change weightedPair π (T h) (T h) ≤ Λ ^ 2 * weightedPair π h h at hspectral
  change |(T f) y₀ / d - μ| ≤ Λ / d * (1 - μ)
  have hone : weightedPair π 1 1 = 1 := by simpa [weightedPair]
  have hf1 : weightedPair π f 1 = μ := by
    rw [weightedPair_indicator_left]; simp [μ]
  have hff : weightedPair π f f = μ := by
    rw [weightedPair_indicator_left]
    exact Finset.sum_congr rfl (by intro y hy; simp [f, orbitIndicator, hy])
  have hf : f = h + μ • (1 : Y → ℝ) := by simp [h]
  have hh1 : weightedPair π h 1 = 0 := by
    have he := congrArg (fun k => weightedPair π k 1) hf
    rw [weightedPair_add_left, weightedPair_smul_left, hf1, hone] at he
    linarith
  have h1h : weightedPair π 1 h = 0 := (weightedPair_comm ..).trans hh1
  have hhh : weightedPair π h h = μ * (1 - μ) := by
    have he := congrArg (fun k => weightedPair π k k) hf
    rw [hff, weightedPair_add_left, weightedPair_add_right, weightedPair_smul_right,
      weightedPair_smul_left, weightedPair_add_right, weightedPair_smul_right,
      hh1, h1h, hone] at he
    linarith
  have hTf1 : weightedPair π (T h) 1 = 0 := by
    rw [hself, hunit, weightedPair_smul_right, hh1, mul_zero]
  have h1Tf : weightedPair π 1 (T h) = 0 := (weightedPair_comm ..).trans hTf1
  have havg : weightedPair π f (T f) = μ * (T f) y₀ := by
    rw [weightedPair_indicator_left, Finset.sum_congr rfl (fun y hy => by rw [hconstant y hy]),
      Finset.sum_mul]
  have hmain : μ * (T f) y₀ = d * μ ^ 2 + weightedPair π h (T h) := by
    calc
      μ * (T f) y₀ = weightedPair π f (T f) := havg.symm
      _ = _ := by
        rw [hf, map_add, map_smul, hunit, smul_smul,
          weightedPair_add_left, weightedPair_add_right, weightedPair_smul_right,
          weightedPair_smul_left, weightedPair_add_right, weightedPair_smul_right,
          hh1, h1Tf, hone]
        ring
  have hbound := weightedPair_abs_le_of_norm_le hπ hΛ hspectral
  rw [hhh] at hbound
  have hpos : 0 < d * μ := mul_pos hd hμ
  have herr : T f y₀ / d - μ = weightedPair π h (T h) / (d * μ) := by
    field_simp
    linarith
  rw [herr, abs_div, abs_of_pos hpos, div_le_iff₀ hpos]
  have hcancel : Λ / d * (1 - μ) * (d * μ) = Λ * (μ * (1 - μ)) := by
    field_simp
  rwa [hcancel]

/-- The orbit estimate with a second normalized eigenfunction.
Subtracting both projections leaves squared norm μ − μ² − b². -/
theorem spectral_orbit_two_eigenfunctions [DecidableEq Y] {π : Y → ℝ}
    (hπ : ∀ y, 0 ≤ π y) (hprob : ∑ y, π y = 1)
    (Ω : Finset Y) (y₀ : Y) (T : (Y → ℝ) →ₗ[ℝ] (Y → ℝ)) (ψ : Y → ℝ)
    {d eigen Λ : ℝ} (hΛ : 0 ≤ Λ)
    (hunit : T 1 = d • 1) (hchar : T ψ = eigen • ψ)
    (horth : weightedPair π 1 ψ = 0) (hnorm : weightedPair π ψ ψ = 1)
    (hself : ∀ f g, weightedPair π (T f) g = weightedPair π f (T g))
    (hconstant : ∀ y ∈ Ω, T (orbitIndicator Ω) y = T (orbitIndicator Ω) y₀)
    (hspectral :
      let f := orbitIndicator Ω
      let h := f - (∑ y ∈ Ω, π y) • (1 : Y → ℝ) - weightedPair π f ψ • ψ
      weightedPair π (T h) (T h) ≤ Λ ^ 2 * weightedPair π h h) :
    let μ := ∑ y ∈ Ω, π y
    let b := weightedPair π (orbitIndicator Ω) ψ
    |μ * T (orbitIndicator Ω) y₀ - (d * μ ^ 2 + eigen * b ^ 2)| ≤ Λ * (μ - μ ^ 2 - b ^ 2) := by
  let μ := ∑ y ∈ Ω, π y
  let f := orbitIndicator Ω
  let b := weightedPair π f ψ
  let h := f - μ • (1 : Y → ℝ) - b • ψ
  change weightedPair π (T h) (T h) ≤ Λ ^ 2 * weightedPair π h h at hspectral
  have hone : weightedPair π 1 1 = 1 := by simpa [weightedPair]
  have horth' : weightedPair π ψ 1 = 0 := (weightedPair_comm ..).trans horth
  have hf1 : weightedPair π f 1 = μ := by rw [weightedPair_indicator_left]; simp [μ]
  have hfψ : weightedPair π f ψ = b := rfl
  have hff : weightedPair π f f = μ := by
    rw [weightedPair_indicator_left]
    exact Finset.sum_congr rfl (by intro y hy; simp [f, orbitIndicator, hy])
  have hf : f = h + μ • (1 : Y → ℝ) + b • ψ := by dsimp [h]; abel
  have hh1 : weightedPair π h 1 = 0 := by
    have he := congrArg (fun g => weightedPair π g 1) hf
    simp only [weightedPair_add_left, weightedPair_smul_left, hf1, hone, horth', mul_one, mul_zero] at he
    linarith only [he]
  have hhψ : weightedPair π h ψ = 0 := by
    have he := congrArg (fun g => weightedPair π g ψ) hf
    simp only [weightedPair_add_left, weightedPair_smul_left, hfψ, hnorm, horth, mul_one, mul_zero] at he
    linarith only [he]
  have h1h : weightedPair π 1 h = 0 := (weightedPair_comm ..).trans hh1
  have hψh : weightedPair π ψ h = 0 := (weightedPair_comm ..).trans hhψ
  have hhh : weightedPair π h h = μ - μ ^ 2 - b ^ 2 := by
    have he := congrArg (fun g => weightedPair π g g) hf
    simp only [weightedPair_add_left, weightedPair_add_right, weightedPair_smul_left, weightedPair_smul_right,
      hff, hh1, hhψ, h1h, hψh, horth, horth', hone, hnorm, mul_one, mul_zero, add_zero, zero_add] at he
    nlinarith only [he]
  have hTh1 : weightedPair π (T h) 1 = 0 := by rw [hself, hunit, weightedPair_smul_right, hh1, mul_zero]
  have hThψ : weightedPair π (T h) ψ = 0 := by rw [hself, hchar, weightedPair_smul_right, hhψ, mul_zero]
  have h1Th : weightedPair π 1 (T h) = 0 := (weightedPair_comm ..).trans hTh1
  have hψTh : weightedPair π ψ (T h) = 0 := (weightedPair_comm ..).trans hThψ
  have havg : weightedPair π f (T f) = μ * (T f) y₀ := by
    rw [weightedPair_indicator_left, Finset.sum_congr rfl (fun y hy => by rw [hconstant y hy]),
      Finset.sum_mul]
  have hmain : μ * (T f) y₀ = d * μ ^ 2 + eigen * b ^ 2 + weightedPair π h (T h) := by
    calc
      _ = weightedPair π f (T f) := havg.symm
      _ = _ := by
        rw [hf]
        simp only [map_add, map_smul, hunit, hchar, smul_smul, weightedPair_add_left, weightedPair_add_right,
          weightedPair_smul_left, weightedPair_smul_right, hh1, hhψ, h1Th, hψTh,
          horth, horth', hone, hnorm, mul_one, mul_zero, add_zero, zero_add]
        ring
  have hbound := weightedPair_abs_le_of_norm_le hπ hΛ hspectral
  rw [hhh] at hbound
  change |μ * (T f) y₀ - (d * μ ^ 2 + eigen * b ^ 2)| ≤ _
  rwa [hmain, add_sub_cancel_left]

omit [Fintype Y] in
/-- Equivariance supplies the orbit-constancy hypothesis of the spectral estimate. -/
theorem orbit_indicator_image_constant {G : Type*} [Group G] [MulAction G Y] [DecidableEq Y]
    (Ω : Finset Y) (y₀ : Y) (T : (Y → ℝ) →ₗ[ℝ] (Y → ℝ)) (K : Subgroup G)
    (horbit : ∀ y, y ∈ Ω ↔ ∃ k ∈ K, k • y₀ = y)
    (hequiv : ∀ (k : G) f y, T (fun z => f (k • z)) y = T f (k • y)) :
    ∀ y ∈ Ω, T (orbitIndicator Ω) y = T (orbitIndicator Ω) y₀ := by
  intro y hy
  obtain ⟨k, hk, rfl⟩ := (horbit y).mp hy
  rw [← hequiv]
  suffices (fun z => orbitIndicator Ω (k • z)) = orbitIndicator Ω by rw [this]
  funext z
  have hb : k • z ∈ Ω ↔ z ∈ Ω := by
    constructor
    · intro hz
      obtain ⟨a, ha, hh⟩ := (horbit _).mp hz
      exact (horbit _).mpr ⟨k⁻¹ * a, K.mul_mem (K.inv_mem hk) ha, by
        rw [mul_smul, hh, inv_smul_smul]⟩
    · intro hz
      obtain ⟨a, ha, rfl⟩ := (horbit _).mp hz
      exact (horbit _).mpr ⟨k * a, K.mul_mem hk ha, mul_smul ..⟩
  simp [orbitIndicator, hb]

theorem weightedPair_sub_left (π f g h : Y → ℝ) :
    weightedPair π (f - g) h = weightedPair π f h - weightedPair π g h := by
  simp only [weightedPair, Pi.sub_apply, sub_mul, mul_sub, Finset.sum_sub_distrib]

/-- Detailed balance makes a finite kernel self-adjoint for the weighted pairing. -/
theorem weightedPair_kernel_self_adjoint (π f g : Y → ℝ) (K : Y → Y → ℝ)
    (hbalance : ∀ x y, π x * K x y = π y * K y x) :
    weightedPair π (fun x => ∑ y, K x y * f y) g =
      weightedPair π f (fun x => ∑ y, K x y * g y) := by
  unfold weightedPair
  calc
    _ = ∑ x, ∑ y, π x * K x y * (f y * g x) := by
      simp only [Finset.sum_mul, Finset.mul_sum]
      refine Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => ?_
      ring
    _ = ∑ y, ∑ x, π y * K y x * (f y * g x) := by
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => by rw [hbalance x y]
    _ = _ := by
      simp only [Finset.mul_sum]
      refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => ?_
      ring

/-- A nonnegative reversible kernel with row sum c has weighted norm at most c. -/
theorem weightedPair_nonnegative_kernel_le (π f : Y → ℝ) (K : Y → Y → ℝ) {c : ℝ}
    (hπ : ∀ x, 0 ≤ π x) (hK : ∀ x y, 0 ≤ K x y)
    (hbalance : ∀ x y, π x * K x y = π y * K y x) (hrow : ∀ x, ∑ y, K x y = c) :
    weightedPair π (fun x => ∑ y, K x y * f y) (fun x => ∑ y, K x y * f y) ≤
      c ^ 2 * weightedPair π f f := by
  have hcauchy (x : Y) : (∑ y, K x y * f y) ^ 2 ≤ c * ∑ y, K x y * (f y * f y) := by
    have h := weightedPair_cauchySchwarz (K x) 1 f (hK x)
    simp only [weightedPair, Pi.one_apply, one_mul, mul_one, hrow] at h
    exact h
  have hstat (y : Y) : ∑ x, π x * K x y = π y * c := by
    simp_rw [hbalance]
    rw [← Finset.mul_sum, hrow]
  unfold weightedPair
  calc
    _ = ∑ x, π x * (∑ y, K x y * f y) ^ 2 := by simp only [pow_two]
    _ ≤ ∑ x, π x * (c * ∑ y, K x y * (f y * f y)) :=
      Finset.sum_le_sum fun x _ => mul_le_mul_of_nonneg_left (hcauchy x) (hπ x)
    _ = c * ∑ y, (∑ x, π x * K x y) * (f y * f y) := by
      simp only [Finset.mul_sum, Finset.sum_mul]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => ?_
      ring
    _ = _ := by
      simp only [hstat, Finset.mul_sum]
      refine Finset.sum_congr rfl fun y _ => ?_
      ring

end QuadraticNorms
