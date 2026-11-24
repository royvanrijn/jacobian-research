import Mathlib.Analysis.Fourier.ZMod
import Mathlib.NumberTheory.DirichletCharacter.Orthogonality
import Mathlib.NumberTheory.MulChar.Lemmas
import Mathlib.Data.Int.GCD

open scoped BigOperators

namespace QuadraticNorms

/-- Squared mass in each fiber of a finite weighted family. -/
noncomputable def pointwiseFiberEnergy {ι α : Type*} [Fintype α] [DecidableEq α]
    (s : Finset ι) (r : ι → α) (F : ι → ℂ) : ℝ :=
  ∑ a, ‖∑ i ∈ s.filter (fun i => r i = a), F i‖ ^ 2

/-- The energy of a weighted integer family after reduction modulo q. -/
noncomputable def pointwiseResidueEnergy {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (q : ℕ) [NeZero q] : ℝ :=
  pointwiseFiberEnergy s (fun i => (n i : ZMod q)) F

private theorem complex_sq_norm (z : ℂ) : ((‖z‖ ^ 2 : ℝ) : ℂ) = z * star z := by
  rw [← Complex.normSq_eq_norm_sq, ← Complex.mul_conj, starRingEnd_apply]

private theorem finite_twisted_energy_expand {ι κ : Type*} [Fintype κ]
    (s : Finset ι) (F : ι → ℂ) (ψ : κ → ι → ℂ) :
    ((∑ k, ‖∑ i ∈ s, F i * ψ k i‖ ^ 2 : ℝ) : ℂ) =
      ∑ i ∈ s, ∑ j ∈ s, F i * star (F j) * ∑ k, ψ k i * star (ψ k j) := by
  rw [Complex.ofReal_sum]
  simp only [complex_sq_norm, ← starRingEnd_apply, map_sum, map_mul,
    Finset.sum_mul, Finset.mul_sum]
  conv_rhs => rw [Finset.sum_comm (s := s) (t := s)]
  rw [Finset.sum_comm (s := Finset.univ) (t := s)]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.sum_comm (s := Finset.univ) (t := s)]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
  ring

/-- Expanding the squared fiber sums counts exactly the pairs with equal residues. -/
theorem pointwise_fiber_energy_eq_collisions {ι α : Type*} [Fintype α] [DecidableEq α]
    (s : Finset ι) (r : ι → α) (F : ι → ℂ) :
    (pointwiseFiberEnergy s r F : ℂ) =
      ∑ i ∈ s, ∑ j ∈ s, if r i = r j then F i * star (F j) else 0 := by
  have hsum (a : α) : ∑ i ∈ s.filter (fun i => r i = a), F i =
      ∑ i ∈ s, F i * (if r i = a then 1 else 0 : ℂ) := by
    simp only [Finset.sum_filter, mul_ite, mul_one, mul_zero]
  simp only [pointwiseFiberEnergy, hsum]
  rw [finite_twisted_energy_expand s F (fun a i => if r i = a then 1 else 0)]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.sum_eq_single_of_mem (r i) (Finset.mem_univ _)]
  · by_cases h : r i = r j
    · simp [h]
    · simp [h, Ne.symm h]
  · intro a _ ha
    simp [Ne.symm ha]

/-- The additive phase e(θn) used in the reflected comparison. -/
noncomputable def pointwiseAdditivePhase (θ : ℝ) (n : ℤ) : ℂ :=
  Complex.exp (2 * (Real.pi : ℂ) * Complex.I * (θ : ℂ) * (n : ℂ))

/-- The Fourier sum of a finite weighted integer family. -/
noncomputable def pointwiseFourierSum {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (θ : ℝ) : ℂ :=
  ∑ i ∈ s, F i * pointwiseAdditivePhase θ (n i)

private theorem pointwise_additive_phase_add (α β : ℝ) (n : ℤ) :
    pointwiseAdditivePhase (α + β) n =
      pointwiseAdditivePhase α n * pointwiseAdditivePhase β n := by
  unfold pointwiseAdditivePhase
  rw [← Complex.exp_add, Complex.ofReal_add]
  congr 1
  ring

private theorem pointwise_additive_phase_grid {q : ℕ} [NeZero q] (a : ZMod q) (n : ℤ) :
    pointwiseAdditivePhase ((a.val : ℝ) / q) n = ZMod.stdAddChar ((n : ZMod q) * a) := by
  have ha : ((n * (a.val : ℤ) : ℤ) : ZMod q) = (n : ZMod q) * a := by
    rw [Int.cast_mul, Int.cast_natCast, ZMod.natCast_zmod_val]
  rw [← ha, ZMod.stdAddChar_coe, pointwiseAdditivePhase,
    Complex.ofReal_div, Complex.ofReal_natCast, Complex.ofReal_natCast, Int.cast_mul, Int.cast_natCast]
  congr 1
  ring

private theorem pointwise_additive_character_orthogonality {q : ℕ} [NeZero q]
    (x y : ZMod q) :
    ∑ a : ZMod q, ZMod.stdAddChar (x * a) * star (ZMod.stdAddChar (y * a)) =
      if x = y then (q : ℂ) else 0 := by
  calc
    _ = ∑ a : ZMod q, ZMod.stdAddChar ((x - y) * a) := by
      refine Finset.sum_congr rfl fun a _ => ?_
      have hy : ‖ZMod.stdAddChar (y * a)‖ = 1 := by
        rw [ZMod.stdAddChar_apply]
        exact Circle.norm_coe _
      rw [sub_mul, AddChar.map_sub_eq_div, div_eq_mul_inv, Complex.inv_eq_conj hy,
        starRingEnd_apply]
    _ = _ := by
      by_cases h : x = y
      · simp [h, ZMod.card]
      · rw [if_neg h]
        exact AddChar.sum_eq_zero_of_ne_one (ZMod.isPrimitive_stdAddChar q (sub_ne_zero.mpr h))

/-- Parseval on the translated target grid j/q + β, including all q grid points. -/
theorem pointwise_fourier_grid_parseval {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (q : ℕ) [NeZero q] (β : ℝ) :
    ∑ a : ZMod q, ‖pointwiseFourierSum s n F ((a.val : ℝ) / q + β)‖ ^ 2 =
      (q : ℝ) * pointwiseResidueEnergy s n
        (fun i => F i * pointwiseAdditivePhase β (n i)) q := by
  let G : ι → ℂ := fun i => F i * pointwiseAdditivePhase β (n i)
  have hfourier (a : ZMod q) : pointwiseFourierSum s n F ((a.val : ℝ) / q + β) =
      ∑ i ∈ s, G i * ZMod.stdAddChar ((n i : ZMod q) * a) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [pointwise_additive_phase_add, pointwise_additive_phase_grid]
    dsimp [G]
    ring
  simp only [hfourier]
  apply Complex.ofReal_injective
  rw [finite_twisted_energy_expand, Complex.ofReal_mul, Complex.ofReal_natCast,
    pointwiseResidueEnergy, pointwise_fiber_energy_eq_collisions, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [pointwise_additive_character_orthogonality]
  split
  · dsimp [G]
    ring
  · rw [mul_zero, mul_zero]

private theorem pointwise_dilated_cast_eq_iff {q A : ℕ} [NeZero q] (n m : ℤ) :
    (((A : ℤ) * n : ℤ) : ZMod q) = (((A : ℤ) * m : ℤ) : ZMod q) ↔
      (n : ZMod (q / q.gcd A)) = (m : ZMod (q / q.gcd A)) := by
  rw [ZMod.intCast_eq_intCast_iff, ZMod.intCast_eq_intCast_iff]
  constructor
  · intro h
    simpa only [Int.gcd_natCast_natCast, Int.natCast_ediv] using
      h.cancel_left_div_gcd (by exact_mod_cast NeZero.pos q)
  · intro h
    have hA : (q.gcd A : ℤ) * (A / q.gcd A : ℕ) = (A : ℤ) := by
      exact_mod_cast Nat.mul_div_cancel' (Nat.gcd_dvd_right q A)
    have hq : (q.gcd A : ℤ) * (q / q.gcd A : ℕ) = (q : ℤ) := by
      exact_mod_cast Nat.mul_div_cancel' (Nat.gcd_dvd_left q A)
    have hh := (h.mul_left (c := ((A / q.gcd A : ℕ) : ℤ))).mul_left'
        (c := (q.gcd A : ℤ))
    simpa only [← mul_assoc, hA, hq] using hh

/-- Dilating all labels by A replaces q by q/gcd(q,A) in the energy. -/
theorem pointwise_residue_energy_dilate {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (q : ℕ) [NeZero q] (A : ℕ) :
    letI : NeZero (q / q.gcd A) := ⟨Nat.ne_of_gt
      (Nat.div_pos (Nat.le_of_dvd (NeZero.pos q) (Nat.gcd_dvd_left q A))
        (Nat.gcd_pos_of_pos_left _ (NeZero.pos q)))⟩
    pointwiseResidueEnergy s (fun i => (A : ℤ) * n i) F q =
      pointwiseResidueEnergy s n F (q / q.gcd A) := by
  apply Complex.ofReal_injective
  simp only [pointwiseResidueEnergy, pointwise_fiber_energy_eq_collisions,
    pointwise_dilated_cast_eq_iff]

private theorem pointwise_multiplicative_character_orthogonality {q : ℕ} [NeZero q]
    (x y : ZMod q) (hy : IsUnit y) :
    ∑ χ : DirichletCharacter ℂ q, χ x * star (χ y) =
      if x = y then (q.totient : ℂ) else 0 := by
  have hinv : Ring.inverse y = y⁻¹ :=
    left_inv_eq_right_inv (Ring.inverse_mul_cancel _ hy) (ZMod.mul_inv_of_unit y hy)
  simpa only [MulChar.star_apply', MulChar.inv_apply, hinv, mul_comm, eq_comm] using
    DirichletCharacter.sum_char_inv_mul_char_eq ℂ hy x

/-- The normalized character energy of the labels divisible by d,
with the characters evaluated on the integer cofactors. -/
noncomputable def pointwiseGcdCharacterEnergy {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (q d : ℕ) : ℝ :=
  if h : 0 < q / d then
    letI : NeZero (q / d) := ⟨Nat.ne_of_gt h⟩
    (∑ χ : DirichletCharacter ℂ (q / d),
      ‖∑ i ∈ s.filter (fun i => (d : ℤ) ∣ n i), F i * χ ((n i / (d : ℤ) : ℤ) : ZMod (q / d))‖ ^ 2) /
      ((q / d).totient : ℝ)
  else 0

private theorem pointwise_cofactor_unit_iff {q d : ℕ} (hd : 0 < d) (hdq : d ∣ q)
    (a : ℤ) (hda : (d : ℤ) ∣ a) :
    IsUnit ((a / (d : ℤ) : ℤ) : ZMod (q / d)) ↔ a.gcd (q : ℤ) = d := by
  have hmul : a.gcd (q : ℤ) = d * (a / (d : ℤ)).gcd ((q : ℤ) / (d : ℤ)) := by
    calc
      _ = ((d : ℤ) * (a / (d : ℤ))).gcd ((d : ℤ) * ((q : ℤ) / (d : ℤ))) := by
        rw [Int.mul_ediv_cancel' hda, Int.mul_ediv_cancel' (Int.natCast_dvd_natCast.mpr hdq)]
      _ = _ := by rw [Int.gcd_mul_left, Int.natAbs_natCast]
  rw [ZMod.coe_int_isUnit_iff_isCoprime, Int.isCoprime_iff_gcd_eq_one,
    Int.gcd_comm ((q / d : ℕ) : ℤ), Int.natCast_ediv, hmul]
  constructor
  · intro h
    rw [h, mul_one]
  · intro h
    apply Nat.eq_of_mul_eq_mul_left hd
    simpa using h

private theorem pointwise_cofactor_cast_eq_iff {q d : ℕ} (hd : 0 < d) (hdq : d ∣ q)
    (a b : ℤ) (ha : (d : ℤ) ∣ a) (hb : (d : ℤ) ∣ b) :
    ((a / (d : ℤ) : ℤ) : ZMod (q / d)) = ((b / (d : ℤ) : ℤ) : ZMod (q / d)) ↔
      (a : ZMod q) = (b : ZMod q) := by
  rw [ZMod.intCast_eq_intCast_iff, ZMod.intCast_eq_intCast_iff, Int.natCast_ediv]
  have hdq' := Int.natCast_dvd_natCast.mpr hdq
  constructor
  · exact fun h => h.of_div ha hb hdq'
  · intro h
    have hd' : (d : ℤ) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hd
    apply (Int.ModEq.mul_left_cancel_iff' hd').mp
    simpa only [Int.mul_ediv_cancel' ha, Int.mul_ediv_cancel' hb,
      Int.mul_ediv_cancel' hdq'] using h

private theorem pointwise_gcd_eq_of_cast_eq {q : ℕ} {a b : ℤ}
    (h : (a : ZMod q) = (b : ZMod q)) : a.gcd (q : ℤ) = b.gcd (q : ℤ) := by
  obtain ⟨k, hk⟩ := Int.modEq_iff_dvd.mp ((ZMod.intCast_eq_intCast_iff a b q).mp h)
  have hb : b = a + (q : ℤ) * k := by omega
  conv_rhs => rw [hb, Int.gcd_add_mul_left_left]

private theorem pointwise_gcd_character_energy_collisions {ι : Type*}
    (s : Finset ι) (n : ι → ℤ) (F : ι → ℂ) (q : ℕ) [NeZero q]
    {d : ℕ} (hd : d ∈ q.divisors) :
    (pointwiseGcdCharacterEnergy s n F q d : ℂ) =
      ∑ i ∈ s.filter (fun i => (n i).gcd (q : ℤ) = d),
        ∑ j ∈ s, if (n i : ZMod q) = (n j : ZMod q) then F i * star (F j) else 0 := by
  classical
  have hdpos := Nat.pos_of_mem_divisors hd
  have hdq := (Nat.mem_divisors.mp hd).1
  have hr : 0 < q / d := Nat.div_pos (Nat.le_of_dvd (NeZero.pos q) hdq) hdpos
  letI : NeZero (q / d) := ⟨Nat.ne_of_gt hr⟩
  let t := s.filter (fun i => (n i).gcd (q : ℤ) = d)
  have hdvd (i : ι) (hi : i ∈ t) : (d : ℤ) ∣ n i := by
    have hg := (Finset.mem_filter.mp hi).2
    rw [← hg]
    exact Int.gcd_dvd_left _ _
  have hchar (χ : DirichletCharacter ℂ (q / d)) :
      ∑ i ∈ s.filter (fun i => (d : ℤ) ∣ n i), F i * χ ((n i / (d : ℤ) : ℤ) : ZMod (q / d)) =
        ∑ i ∈ t, F i * χ ((n i / (d : ℤ) : ℤ) : ZMod (q / d)) := by
    symm
    apply Finset.sum_subset
    · intro i hi
      exact Finset.mem_filter.mpr ⟨(Finset.mem_filter.mp hi).1, hdvd i hi⟩
    · intro i hi hnot
      have hn : ¬ (n i).gcd (q : ℤ) = d := fun hg =>
        hnot (Finset.mem_filter.mpr ⟨(Finset.mem_filter.mp hi).1, hg⟩)
      rw [MulChar.map_nonunit χ
        (mt (pointwise_cofactor_unit_iff hdpos hdq (n i) (Finset.mem_filter.mp hi).2).mp hn),
        mul_zero]
  unfold pointwiseGcdCharacterEnergy
  rw [dif_pos hr, Complex.ofReal_div, Complex.ofReal_natCast]
  simp only [hchar]
  rw [finite_twisted_energy_expand]
  have hφ : ((q / d).totient : ℂ) ≠ 0 := by
    exact_mod_cast Nat.ne_of_gt (Nat.totient_pos.mpr hr)
  apply (div_eq_iff hφ).mpr
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Finset.sum_mul]
  calc
    _ = ∑ j ∈ t, (if (n i : ZMod q) = (n j : ZMod q) then F i * star (F j) else 0) *
        ((q / d).totient : ℂ) := by
      refine Finset.sum_congr rfl fun j hj => ?_
      have hu : IsUnit ((n j / (d : ℤ) : ℤ) : ZMod (q / d)) :=
        (pointwise_cofactor_unit_iff hdpos hdq (n j) (hdvd j hj)).mpr (Finset.mem_filter.mp hj).2
      rw [pointwise_multiplicative_character_orthogonality _ _ hu,
        if_congr (pointwise_cofactor_cast_eq_iff hdpos hdq (n i) (n j) (hdvd i hi) (hdvd j hj))
          rfl rfl]
      split <;> simp only [zero_mul, mul_zero]
    _ = _ := by
      apply Finset.sum_subset (Finset.filter_subset _ _)
      intro j hj hnot
      have hn : (n i : ZMod q) ≠ (n j : ZMod q) := fun he =>
        hnot (Finset.mem_filter.mpr
          ⟨hj, (pointwise_gcd_eq_of_cast_eq he).symm.trans (Finset.mem_filter.mp hi).2⟩)
      rw [if_neg hn, zero_mul]

/-- Multiplicative orthogonality summed over every gcd profile of q.
The character zeros restrict the cofactor to units, including when q/d is one. -/
theorem pointwise_residue_energy_eq_gcd_character_energy {ι : Type*} (s : Finset ι)
    (n : ι → ℤ) (F : ι → ℂ) (q : ℕ) [NeZero q] :
    pointwiseResidueEnergy s n F q =
      ∑ d ∈ q.divisors, pointwiseGcdCharacterEnergy s n F q d := by
  apply Complex.ofReal_injective
  rw [Complex.ofReal_sum, pointwiseResidueEnergy, pointwise_fiber_energy_eq_collisions]
  symm
  calc
    _ = ∑ d ∈ q.divisors, ∑ i ∈ s.filter (fun i => (n i).gcd (q : ℤ) = d),
          ∑ j ∈ s, if (n i : ZMod q) = (n j : ZMod q) then F i * star (F j) else 0 := by
      exact Finset.sum_congr rfl fun _ hd => pointwise_gcd_character_energy_collisions s n F q hd
    _ = _ := by
      apply Finset.sum_fiberwise_of_maps_to
      intro i _
      exact Nat.mem_divisors.mpr
        ⟨Int.natCast_dvd_natCast.mp (Int.gcd_dvd_right _ _), NeZero.ne q⟩

private theorem pointwise_sq_norm_sum_le {ι : Type*} (s : Finset ι) (F : ι → ℂ)
    (M : ℕ) (hM : s.card ≤ M) :
    ‖∑ i ∈ s, F i‖ ^ 2 ≤ (M : ℝ) * ∑ i ∈ s, ‖F i‖ ^ 2 := by
  calc
    _ ≤ (∑ i ∈ s, ‖F i‖) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) (norm_sum_le s F) 2
    _ ≤ (s.card : ℝ) * ∑ i ∈ s, ‖F i‖ ^ 2 := by
      simpa using Finset.sum_mul_sq_le_sq_mul_sq s (fun _ => (1 : ℝ)) (fun i => ‖F i‖)
    _ ≤ _ := by gcongr

/-- A bound on the largest fiber gives a squared-sum bound for the energy. -/
theorem pointwise_fiber_energy_le_sq_norm_sum {ι α : Type*} [Fintype α] [DecidableEq α]
    (s : Finset ι) (r : ι → α) (F : ι → ℂ) (M : ℕ)
    (hM : ∀ a, (s.filter fun i => r i = a).card ≤ M) :
    pointwiseFiberEnergy s r F ≤ (M : ℝ) * ∑ i ∈ s, ‖F i‖ ^ 2 := by
  unfold pointwiseFiberEnergy
  calc
    _ ≤ ∑ a, (M : ℝ) * ∑ i ∈ s.filter (fun i => r i = a), ‖F i‖ ^ 2 :=
      Finset.sum_le_sum fun a _ => pointwise_sq_norm_sum_le _ F M (hM a)
    _ = _ := by rw [← Finset.mul_sum, Finset.sum_fiberwise]

private theorem pointwise_fiber_energy_le_refinement {ι α γ : Type*} [Fintype α]
    [Fintype γ] [DecidableEq α] [DecidableEq γ] (s : Finset ι) (r : ι → γ)
    (c : γ → α) (F : ι → ℂ) (M : ℕ)
    (hM : ∀ a, (Finset.univ.filter fun b => c b = a).card ≤ M) :
    pointwiseFiberEnergy s (fun i => c (r i)) F ≤ (M : ℝ) * pointwiseFiberEnergy s r F := by
  let G (b : γ) := ∑ i ∈ s.filter (fun i => r i = b), F i
  have hg (a : α) : ∑ i ∈ s.filter (fun i => c (r i) = a), F i =
      ∑ b ∈ Finset.univ.filter (fun b => c b = a), G b := by
    dsimp only [G]
    rw [Finset.sum_fiberwise_eq_sum_filter]
    simp
  simpa only [pointwiseFiberEnergy, hg, G] using
    pointwise_fiber_energy_le_sq_norm_sum Finset.univ c G M hM

private theorem pointwise_mod_reduction_card_le {q Q : ℕ} [NeZero q] [NeZero Q]
    (hqQ : q ∣ Q) (a : ZMod q) :
    (Finset.univ.filter fun b => ZMod.castHom hqQ (ZMod q) b = a).card ≤ Q / q := by
  rw [← Finset.card_range (Q / q)]
  apply Finset.card_le_card_of_injOn (fun b : ZMod Q => b.val / q)
  · intro b _
    simp only [Finset.coe_range]
    exact (Nat.div_lt_iff_lt_mul (NeZero.pos q)).mpr
      (by simpa only [Nat.div_mul_cancel hqQ] using b.val_lt)
  · intro b hb c hc h
    change b.val / q = c.val / q at h
    apply ZMod.val_injective Q
    have hcast (x : ZMod Q) : (ZMod.castHom hqQ (ZMod q) x).val = x.val % q := by
      rw [ZMod.castHom_apply, ZMod.cast_eq_val, ZMod.val_natCast]
    have hr : b.val % q = c.val % q := by
      rw [← hcast b, ← hcast c, (Finset.mem_filter.mp hb).2, (Finset.mem_filter.mp hc).2]
    calc
      b.val = b.val % q + q * (b.val / q) := (Nat.mod_add_div _ _).symm
      _ = c.val % q + q * (c.val / q) := by rw [hr, h]
      _ = c.val := Nat.mod_add_div _ _

/-- Refinement from q to a multiple Q loses at most the index Q/q in energy. -/
theorem pointwise_residue_energy_le_refinement {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) {q Q : ℕ} [NeZero q] [NeZero Q] (hqQ : q ∣ Q) :
    pointwiseResidueEnergy s n F q ≤ (Q / q : ℕ) * pointwiseResidueEnergy s n F Q := by
  simpa only [pointwiseResidueEnergy, ZMod.castHom_apply, ZMod.cast_intCast hqQ] using
    pointwise_fiber_energy_le_refinement s (fun i => (n i : ZMod Q))
      (ZMod.castHom hqQ (ZMod q)) F (Q / q) (pointwise_mod_reduction_card_le hqQ)

/-- Refining to lcm(q,L) costs at most L in the manuscript's energy estimate. -/
theorem pointwise_residue_energy_le_lcm {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (q L : ℕ) [NeZero q] [NeZero L] :
    letI : NeZero (q.lcm L) := ⟨Nat.ne_of_gt (Nat.lcm_pos (NeZero.pos q) (NeZero.pos L))⟩
    pointwiseResidueEnergy s n F q ≤ (L : ℝ) * pointwiseResidueEnergy s n F (q.lcm L) := by
  letI : NeZero (q.lcm L) := ⟨Nat.ne_of_gt (Nat.lcm_pos (NeZero.pos q) (NeZero.pos L))⟩
  refine (pointwise_residue_energy_le_refinement s n F (Nat.dvd_lcm_left q L)).trans ?_
  have hle : q.lcm L / q ≤ L := Nat.div_le_of_le_mul
    (Nat.le_of_dvd (Nat.mul_pos (NeZero.pos q) (NeZero.pos L))
      (Nat.lcm_dvd (Nat.dvd_mul_right _ _) (Nat.dvd_mul_left _ _)))
  apply mul_le_mul_of_nonneg_right (by exact_mod_cast hle)
  exact Finset.sum_nonneg fun _ _ => sq_nonneg _

/-- A residue class meets a finite subset of a length-H integer interval
in at most 1 + H/q points. The subset may be a discarded tail. -/
theorem pointwise_residue_energy_le_interval (s : Finset ℤ) (F : ℤ → ℂ) (lo : ℤ)
    (H q : ℕ) [NeZero q] (hs : ∀ i ∈ s, lo ≤ i ∧ i ≤ lo + (H : ℤ)) :
    pointwiseResidueEnergy s id F q ≤ (1 + (H : ℝ) / q) * ∑ i ∈ s, ‖F i‖ ^ 2 := by
  have hcard (a : ZMod q) :
      (s.filter fun (i : ℤ) => (i : ZMod q) = a).card ≤ H / q + 1 := by
    let shift (i : ℤ) := (i - lo).toNat
    have hshift (i : ℤ) (hi : i ∈ s) : (shift i : ℤ) = i - lo :=
      Int.toNat_of_nonneg (sub_nonneg.mpr (hs i hi).1)
    rw [← Finset.card_range (H / q + 1)]
    apply Finset.card_le_card_of_injOn (fun i => shift i / q)
    · intro i hi
      have his := (Finset.mem_filter.mp hi).1
      have hle : shift i ≤ H := by
        apply Int.ofNat_le.mp
        rw [hshift i his]
        linarith [(hs i his).2]
      simpa only [Finset.coe_range, Set.mem_Iio, Nat.lt_succ_iff] using
        Nat.div_le_div_right hle
    · intro i hi j hj h
      obtain ⟨his, hia⟩ := Finset.mem_filter.mp hi
      obtain ⟨hjs, hja⟩ := Finset.mem_filter.mp hj
      have hcast (k : ℤ) (hk : k ∈ s) : (shift k : ZMod q) = (k : ZMod q) - (lo : ZMod q) := by
        rw [← Int.cast_natCast, hshift k hk, Int.cast_sub]
      have he : shift i = shift j := (Nat.residueClassesEquiv q).injective
        (Prod.ext (by change (shift i : ZMod q) = (shift j : ZMod q); rw [hcast i his, hcast j hjs,
          hia, hja]) h)
      have hz := congrArg (fun k : ℕ => (k : ℤ)) he
      rw [hshift i his, hshift j hjs] at hz
      exact sub_left_inj.mp hz
  refine (pointwise_fiber_energy_le_sq_norm_sum s (fun i => (i : ZMod q)) F
    (H / q + 1) hcard).trans ?_
  apply mul_le_mul_of_nonneg_right _ (Finset.sum_nonneg fun _ _ => sq_nonneg _)
  push_cast
  linarith [Nat.cast_div_le (α := ℝ) (m := H) (n := q)]

private theorem pointwise_character_energy_congr (s : Finset ℤ) (F : ℤ → ℂ)
    {r q : ℕ} [NeZero r] [NeZero q] (h : r = q) :
    (∑ χ : DirichletCharacter ℂ r, ‖∑ i ∈ s, F i * χ (i : ZMod r)‖ ^ 2) /
        (r.totient : ℝ) =
      (∑ χ : DirichletCharacter ℂ q, ‖∑ i ∈ s, F i * χ (i : ZMod q)‖ ^ 2) /
        (q.totient : ℝ) := by
  subst r
  rfl

/-- The normalized squared cofactor character sums over any selection of
characters are bounded by the interval occupancy times the squared weight. -/
theorem pointwise_character_energy_le_interval (s : Finset ℤ) (F : ℤ → ℂ) (lo : ℤ)
    (H q : ℕ) [NeZero q] (C : Finset (DirichletCharacter ℂ q))
    (hs : ∀ i ∈ s, lo ≤ i ∧ i ≤ lo + (H : ℤ)) :
    (∑ χ ∈ C, ‖∑ i ∈ s, F i * χ (i : ZMod q)‖ ^ 2) / (q.totient : ℝ) ≤
      (1 + (H : ℝ) / q) * ∑ i ∈ s, ‖F i‖ ^ 2 := by
  have h₁ : 1 ∈ q.divisors := Nat.mem_divisors.mpr ⟨one_dvd q, NeZero.ne q⟩
  have hbound : pointwiseGcdCharacterEnergy s id F q 1 ≤ pointwiseResidueEnergy s id F q := by
    rw [pointwise_residue_energy_eq_gcd_character_energy]
    apply Finset.single_le_sum _ h₁
    intro d _
    unfold pointwiseGcdCharacterEnergy
    split
    · exact div_nonneg (Finset.sum_nonneg fun _ _ => sq_nonneg _) (Nat.cast_nonneg _)
    · exact le_rfl
  have hone : pointwiseGcdCharacterEnergy s id F q 1 =
      (∑ χ : DirichletCharacter ℂ q, ‖∑ i ∈ s, F i * χ (i : ZMod q)‖ ^ 2) /
        (q.totient : ℝ) := by
    have hq : 0 < q / 1 := by simpa using NeZero.pos q
    letI : NeZero (q / 1) := ⟨Nat.ne_of_gt hq⟩
    rw [pointwiseGcdCharacterEnergy, dif_pos hq]
    simp only [id_eq, Nat.cast_one, Int.ediv_one, one_dvd, Finset.sum_filter, if_true]
    exact pointwise_character_energy_congr s F (Nat.div_one q)
  calc
    _ ≤ pointwiseGcdCharacterEnergy s id F q 1 := by
      rw [hone]
      apply div_le_div_of_nonneg_right _ (Nat.cast_nonneg _)
      exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) fun _ _ _ => sq_nonneg _
    _ ≤ pointwiseResidueEnergy s id F q := hbound
    _ ≤ _ := pointwise_residue_energy_le_interval s F lo H q hs

theorem pointwise_additive_phase_continuous (n : ℤ) :
    Continuous (fun θ => pointwiseAdditivePhase θ n) := by
  apply Complex.continuous_exp.comp
  exact (continuous_const.mul Complex.continuous_ofReal).mul continuous_const

theorem pointwise_fourier_sum_continuous {ι : Type*} (s : Finset ι) (n : ι → ℤ) (F : ι → ℂ) :
    Continuous (pointwiseFourierSum s n F) := by
  apply continuous_finsetSum
  intro i _
  exact continuous_const.mul (pointwise_additive_phase_continuous (n i))

/-- A uniform translated residue-energy estimate integrates over the entire
grid with only the interval length and the Parseval normalization q. -/
theorem pointwise_fourier_grid_integral_le {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (q : ℕ) [NeZero q] {lo hi B : ℝ} (hlo : lo ≤ hi)
    (hB : ∀ β ∈ Set.Icc lo hi, pointwiseResidueEnergy s n
      (fun i => F i * pointwiseAdditivePhase β (n i)) q ≤ B / (q : ℝ)) :
    (∫ β in lo..hi, ∑ a : ZMod q, ‖pointwiseFourierSum s n F ((a.val : ℝ) / q + β)‖ ^ 2) ≤
      (hi - lo) * B := by
  have hcont : Continuous (fun β : ℝ =>
      ∑ a : ZMod q, ‖pointwiseFourierSum s n F ((a.val : ℝ) / q + β)‖ ^ 2) :=
    continuous_finsetSum _ fun _ _ =>
      (((pointwise_fourier_sum_continuous s n F).comp (continuous_const.add continuous_id)).norm).pow 2
  calc
    _ ≤ ∫ β in lo..hi, B := by
      apply intervalIntegral.integral_mono_on hlo (hcont.intervalIntegrable _ _)
        (continuous_const.intervalIntegrable _ _)
      intro β hβ
      rw [pointwise_fourier_grid_parseval]
      have hq : (q : ℝ) ≠ 0 := by exact_mod_cast NeZero.ne q
      calc
        _ ≤ (q : ℝ) * (B / q) := mul_le_mul_of_nonneg_left (hB β hβ) (Nat.cast_nonneg _)
        _ = B := by field_simp
    _ = _ := by rw [intervalIntegral.integral_const, smul_eq_mul]

end QuadraticNorms
