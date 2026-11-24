import Mathlib.GroupTheory.Perm.DomMulAct
import Mathlib.GroupTheory.GroupAction.Quotient
import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Algebra.BigOperators.Field

open scoped BigOperators Classical

namespace QuadraticNorms

private theorem bank_mem_permutation_orbit_iff {α β : Type*} [Fintype α] [DecidableEq β]
    (x y : α → β) :
    y ∈ MulAction.orbit (Equiv.Perm α)ᵈᵐᵃ x ↔
      ∀ b, Fintype.card {i // y i = b} = Fintype.card {i // x i = b} := by
  rw [MulAction.mem_orbit_iff]
  constructor
  · rintro ⟨g, rfl⟩ b
    let e : Equiv.Perm α := DomMulAct.mk.symm g
    exact Fintype.card_congr (Equiv.subtypeEquiv e (fun _ => Iff.rfl))
  · intro h
    let e := Equiv.ofFiberEquiv (f := y) (g := x) (fun b => Fintype.equivOfCardEq (h b))
    refine ⟨DomMulAct.mk e, funext (fun i => ?_)⟩
    simpa only [DomMulAct.smul_apply, Equiv.symm_apply_apply, Equiv.Perm.smul_def, e] using
      Equiv.ofFiberEquiv_map (fun b => Fintype.equivOfCardEq (h b)) i

private theorem bank_permutation_fiber_card {α β : Type*} [Fintype α] [DecidableEq α]
    [Fintype β] [DecidableEq β] (x : α → β) :
    Fintype.card {y : α → β // ∀ b, Fintype.card {i // y i = b} = Fintype.card {i // x i = b}} *
      (∏ b : β, (Fintype.card {i // x i = b}).factorial) = (Fintype.card α).factorial := by
  let G := (Equiv.Perm α)ᵈᵐᵃ
  letI : Fintype G := Fintype.ofEquiv (Equiv.Perm α) DomMulAct.mk
  letI : Fintype (MulAction.stabilizer G x) := Fintype.ofFinite _
  have horb : Fintype.card (MulAction.orbit G x) =
      Fintype.card {y : α → β // ∀ b, Fintype.card {i // y i = b} = Fintype.card {i // x i = b}} :=
    Fintype.card_congr (Equiv.subtypeEquiv (Equiv.refl _) (bank_mem_permutation_orbit_iff x))
  have hstab : Fintype.card (MulAction.stabilizer G x) =
      ∏ b : β, (Fintype.card {i // x i = b}).factorial := by
    let e : {g : Equiv.Perm α // x ∘ g = x} ≃ (MulAction.stabilizer G x) :=
      Equiv.subtypeEquiv DomMulAct.mk (fun g => by
        change x ∘ g = x ↔ DomMulAct.mk g ∈ MulAction.stabilizer G x
        rw [DomMulAct.mem_stabilizer_iff, Equiv.symm_apply_apply])
    exact (Fintype.card_congr e.symm).trans (DomMulAct.stabilizer_card x)
  have hgroup : Fintype.card G = (Fintype.card α).factorial :=
    (Fintype.card_congr DomMulAct.mk.symm).trans Fintype.card_perm
  have h := MulAction.card_orbit_mul_card_stabilizer_eq_card_group G x
  rw [horb, hstab, hgroup] at h
  exact h

/-- One ordered draw for every slot of every finite bank. -/
abbrev PrimeBankDraws {G : Type*} (P : G → Finset ℕ) (k : G → ℕ) :=
  ∀ g : G, Fin (k g) → ↥(P g)

/-- The numerical product of the draws, with repetitions included. -/
def primeBankDrawProduct {G : Type*} [Fintype G] {P : G → Finset ℕ} {k : G → ℕ}
    (z : PrimeBankDraws P k) : ℕ :=
  ∏ g, ∏ j : Fin (k g), ((z g j) : ℕ)

section Count

variable {G : Type*} [Fintype G] [DecidableEq G] {P : G → Finset ℕ} {k : G → ℕ}

omit [DecidableEq G] in
private theorem bank_draw_factorization (hp : ∀ g p, p ∈ P g → p.Prime)
    (hdis : Pairwise fun g g' => Disjoint (P g) (P g')) (z : PrimeBankDraws P k)
    (g : G) (p : ↥(P g)) :
    (primeBankDrawProduct z).factorization (p : ℕ) = Fintype.card {j : Fin (k g) // z g j = p} := by
  have hnon (a : G) : (∏ j : Fin (k a), (z a j : ℕ)) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr (fun j _ => (hp a _ (z a j).property).ne_zero)
  rw [primeBankDrawProduct, Nat.factorization_prod_apply (fun a _ => hnon a),
    Finset.sum_eq_single g]
  · rw [Nat.factorization_prod_apply (fun j _ => (hp g _ (z g j).property).ne_zero),
      Fintype.card_subtype, Finset.card_eq_sum_ones, Finset.sum_filter]
    refine Finset.sum_congr rfl (fun j _ => ?_)
    rw [(hp g _ (z g j).property).factorization, Finsupp.single_apply]
    by_cases hj : z g j = p
    · simp [hj]
    · have hv : (z g j : ℕ) ≠ (p : ℕ) := fun h => hj (Subtype.ext h)
      simp [hj, hv]
  · intro a _ hne
    rw [Nat.factorization_prod_apply (fun j _ => (hp a _ (z a j).property).ne_zero)]
    refine Finset.sum_eq_zero (fun j _ => ?_)
    have hn : (z a j : ℕ) ≠ (p : ℕ) := by
      intro h
      exact Finset.disjoint_left.mp (hdis hne) (z a j).property (h.symm ▸ p.property)
    rw [(hp a _ (z a j).property).factorization, Finsupp.single_apply, if_neg hn]
  · simp

omit [DecidableEq G] in
private theorem bank_draw_product_eq_iff (hp : ∀ g p, p ∈ P g → p.Prime)
    (hdis : Pairwise fun g g' => Disjoint (P g) (P g')) (x y : PrimeBankDraws P k) :
    primeBankDrawProduct x = primeBankDrawProduct y ↔
      ∀ g (p : ↥(P g)), Fintype.card {j : Fin (k g) // x g j = p} =
        Fintype.card {j : Fin (k g) // y g j = p} := by
  constructor
  · intro h g p
    rw [← bank_draw_factorization hp hdis x g p, ← bank_draw_factorization hp hdis y g p, h]
  · intro h
    unfold primeBankDrawProduct
    refine Finset.prod_congr rfl (fun g _ => ?_)
    obtain ⟨e, he⟩ := MulAction.mem_orbit_iff.mp ((bank_mem_permutation_orbit_iff (y g) (x g)).mpr (h g))
    exact Fintype.prod_equiv (DomMulAct.mk.symm e : Equiv.Perm (Fin (k g)))
      (fun j => (x g j : ℕ)) (fun j => (y g j : ℕ)) (fun j => by
        have hh := congrFun he j
        rw [DomMulAct.smul_apply, Equiv.Perm.smul_def] at hh
        exact congrArg Subtype.val hh.symm)

private theorem bank_draw_product_fiber_card (hp : ∀ g p, p ∈ P g → p.Prime)
    (hdis : Pairwise fun g g' => Disjoint (P g) (P g')) {n : ℕ}
    (hn : ∃ z : PrimeBankDraws P k, primeBankDrawProduct z = n) :
    (Fintype.card {z : PrimeBankDraws P k // primeBankDrawProduct z = n} : ℝ) =
      ∏ g, ((k g).factorial : ℝ) / ∏ p : ↥(P g), ((n.factorization p).factorial : ℝ) := by
  obtain ⟨x, hx⟩ := hn
  let q (g : G) (z : Fin (k g) → ↥(P g)) :=
    ∀ p : ↥(P g), Fintype.card {j // z j = p} = Fintype.card {j // x g j = p}
  have heq : Fintype.card {z : PrimeBankDraws P k // primeBankDrawProduct z = n} =
      Fintype.card (∀ g, {z // q g z}) := by
    apply Fintype.card_congr
    refine (Equiv.subtypeEquiv (Equiv.refl _) (fun z => ?_)).trans Equiv.subtypePiEquivPi
    rw [← hx]
    exact bank_draw_product_eq_iff hp hdis z x
  rw [heq, Fintype.card_pi, Nat.cast_prod]
  refine Finset.prod_congr rfl (fun g _ => ?_)
  have hfreq (p : ↥(P g)) : Fintype.card {j // x g j = p} = n.factorization (p : ℕ) := by
    rw [← bank_draw_factorization hp hdis x g p, hx]
  have hf : Fintype.card {z // q g z} *
      (∏ p : ↥(P g), (Fintype.card {j // x g j = p}).factorial) = (k g).factorial := by
    simpa only [Fintype.card_fin] using bank_permutation_fiber_card (x g)
  rw [show (∏ p : ↥(P g), (Fintype.card {j // x g j = p}).factorial) =
      ∏ p : ↥(P g), (n.factorization p).factorial from Finset.prod_congr rfl (fun p _ => by rw [hfreq])] at hf
  apply (eq_div_iff (by positivity : (∏ p : ↥(P g), ((n.factorization p).factorial : ℝ)) ≠ 0)).mpr
  exact_mod_cast hf

end Count

/-- Reciprocal mass of a prime bank. -/
noncomputable def primeBankDrawMass (P : Finset ℕ) : ℝ := ∑ p : ↥P, (1 : ℝ) / (p : ℕ)

/-- The product probability weight of one ordered draw. -/
noncomputable def primeBankDrawWeight {G : Type*} [Fintype G] {P : G → Finset ℕ} {k : G → ℕ}
    (z : PrimeBankDraws P k) : ℝ :=
  ∏ g, ∏ j : Fin (k g), 1 / (primeBankDrawMass (P g) * (z g j : ℕ))

/-- The mass of all ordered draws having numerical product n. -/
noncomputable def primeBankProductMass {G : Type*} [Fintype G] [DecidableEq G]
    (P : G → Finset ℕ) (k : G → ℕ) (n : ℕ) : ℝ :=
  ∑ z : PrimeBankDraws P k, if primeBankDrawProduct z = n then primeBankDrawWeight z else 0

private theorem bank_draw_weight_product {G : Type*} [Fintype G] {P : G → Finset ℕ} {k : G → ℕ}
    (z : PrimeBankDraws P k) :
    primeBankDrawWeight z = 1 / ((∏ g, primeBankDrawMass (P g) ^ k g) * (primeBankDrawProduct z : ℝ)) := by
  unfold primeBankDrawWeight primeBankDrawProduct
  simp_rw [Finset.prod_div_distrib, Finset.prod_const_one, Finset.prod_mul_distrib,
    Finset.prod_const, Finset.card_univ, Fintype.card_fin, Nat.cast_prod]

theorem prime_bank_product_mass_eq_zero_of_unrepresented {G : Type*} [Fintype G] [DecidableEq G]
    (P : G → Finset ℕ) (k : G → ℕ) {n : ℕ} (h : ¬ ∃ z : PrimeBankDraws P k, primeBankDrawProduct z = n) :
    primeBankProductMass P k n = 0 := by
  unfold primeBankProductMass
  exact Finset.sum_eq_zero (fun z _ => if_neg (fun hz => h ⟨z, hz⟩))

/-- The exact draw-product law, including factorials for repeated
labels in every bank. -/
theorem prime_bank_product_mass_exact {G : Type*} [Fintype G] [DecidableEq G]
    (P : G → Finset ℕ) (k : G → ℕ) (hp : ∀ g p, p ∈ P g → p.Prime)
    (hdis : Pairwise fun g g' => Disjoint (P g) (P g'))
    (hV : ∀ g, 0 < primeBankDrawMass (P g))
    {n : ℕ} (hn : 0 < n) (h : ∃ z : PrimeBankDraws P k, primeBankDrawProduct z = n) :
    (n : ℝ) * primeBankProductMass P k n =
      ∏ g, ((k g).factorial : ℝ) /
        (primeBankDrawMass (P g) ^ k g * ∏ p : ↥(P g), ((n.factorization p).factorial : ℝ)) := by
  let W := ∏ g, primeBankDrawMass (P g) ^ k g
  have hW : 0 < W := Finset.prod_pos (fun g _ => pow_pos (hV g) _)
  have hnR : 0 < (n : ℝ) := by exact_mod_cast hn
  have hconst : primeBankProductMass P k n =
      (Fintype.card {z : PrimeBankDraws P k // primeBankDrawProduct z = n} : ℝ) / (W * n) := by
    unfold primeBankProductMass
    rw [← Finset.sum_filter, Fintype.card_subtype]
    calc
      _ = ∑ _z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z = n),
          1 / (W * n) := by
        refine Finset.sum_congr rfl (fun z hz => ?_)
        rw [bank_draw_weight_product, (Finset.mem_filter.mp hz).2]
      _ = _ := by simp only [Finset.sum_const, nsmul_eq_mul, mul_one_div]
  rw [hconst, bank_draw_product_fiber_card hp hdis h]
  rw [← mul_div_assoc, mul_comm W, ← div_div, mul_div_cancel_left₀ _ hnR.ne']
  unfold W
  rw [← Finset.prod_div_distrib]
  exact Finset.prod_congr rfl (fun g _ => by rw [div_div, mul_comm])

section WeightBounds

variable {G : Type*} [Fintype G] [DecidableEq G] {P : G → Finset ℕ} {k : G → ℕ}

omit [DecidableEq G] in
private theorem prime_bank_draw_weight_nonneg (z : PrimeBankDraws P k) :
    0 ≤ primeBankDrawWeight z := by
  unfold primeBankDrawWeight primeBankDrawMass
  positivity

/-- The ordered draw law is normalized whenever every bank has
positive reciprocal mass. -/
theorem prime_bank_draw_weight_sum (hV : ∀ g, 0 < primeBankDrawMass (P g)) :
    ∑ z : PrimeBankDraws P k, primeBankDrawWeight z = 1 := by
  have hlocal (g : G) :
      (∑ z : Fin (k g) → ↥(P g), ∏ j : Fin (k g), 1 / (primeBankDrawMass (P g) * (z j : ℕ))) = 1 := by
    have hOne : (∑ p : ↥(P g), 1 / (primeBankDrawMass (P g) * (p : ℕ))) = 1 := by
      calc
        _ = (∑ p : ↥(P g), ((1 : ℝ) / (p : ℕ))) / primeBankDrawMass (P g) := by
          rw [Finset.sum_div]
          exact Finset.sum_congr rfl (fun _ _ => by rw [div_div, mul_comm])
        _ = _ := div_self (hV g).ne'
    simpa only [Fintype.piFinset_univ, hOne, Finset.prod_const_one] using
      Finset.sum_prod_piFinset (Finset.univ : Finset ↥(P g))
        (fun (_ : Fin (k g)) p => 1 / (primeBankDrawMass (P g) * (p : ℕ)))
  simpa only [Fintype.piFinset_univ, hlocal, Finset.prod_const_one, primeBankDrawWeight] using
    (Finset.prod_univ_sum (fun g : G => (Finset.univ : Finset (Fin (k g) → ↥(P g))))
      (fun g z => ∏ j : Fin (k g), 1 / (primeBankDrawMass (P g) * (z j : ℕ)))).symm

private theorem prime_bank_product_mass_mul_le (hp : ∀ g p, p ∈ P g → p.Prime)
    (hdis : Pairwise fun g g' => Disjoint (P g) (P g')) {v : ℝ} (hv : 0 < v)
    (hV : ∀ g, v ≤ primeBankDrawMass (P g)) {m : ℕ} (hm : ∀ g, k g ≤ m)
    (n : ℕ) :
    (n : ℝ) * primeBankProductMass P k n ≤
      ((m.factorial : ℝ) * (max 1 v⁻¹) ^ m) ^ Fintype.card G := by
  let a := max (1 : ℝ) v⁻¹
  have ha₁ : 1 ≤ a := le_max_left ..
  by_cases hn : n = 0
  · rw [hn, Nat.cast_zero, zero_mul]
    exact pow_nonneg (mul_nonneg (by positivity) (pow_nonneg (zero_le_one.trans ha₁) _)) _
  by_cases hz : ∃ z : PrimeBankDraws P k, primeBankDrawProduct z = n
  · have hV0 g : 0 < primeBankDrawMass (P g) := hv.trans_le (hV g)
    rw [prime_bank_product_mass_exact P k hp hdis hV0 (Nat.pos_of_ne_zero hn) hz]
    rw [← Finset.card_univ, ← Finset.prod_const]
    apply Finset.prod_le_prod (fun g _ =>
      div_nonneg (by positivity) (mul_nonneg (pow_nonneg (hV0 g).le _) (by positivity)))
    intro g _
    let V := primeBankDrawMass (P g)
    let F := ∏ p : ↥(P g), ((n.factorization p).factorial : ℝ)
    have hF : 1 ≤ F := Finset.one_le_prod (fun p _ => by
      exact_mod_cast Nat.succ_le_of_lt (Nat.factorial_pos (n.factorization (p : ℕ))))
    have hInv : V⁻¹ ≤ a := ((inv_le_inv₀ (hV0 g) hv).mpr (hV g)).trans (le_max_right ..)
    have hFac : ((k g).factorial : ℝ) ≤ (m.factorial : ℝ) := by exact_mod_cast Nat.factorial_le (hm g)
    change ((k g).factorial : ℝ) / (V ^ k g * F) ≤ (m.factorial : ℝ) * a ^ m
    calc
      _ ≤ ((k g).factorial : ℝ) / V ^ k g := by
        apply div_le_div_of_nonneg_left (by positivity) (pow_pos (hV0 g) _)
        exact le_mul_of_one_le_right (pow_nonneg (hV0 g).le _) hF
      _ = ((k g).factorial : ℝ) * V⁻¹ ^ k g := by rw [inv_pow, div_eq_mul_inv]
      _ ≤ _ := mul_le_mul hFac
        ((pow_le_pow_left₀ (inv_nonneg.mpr (hV0 g).le) hInv _).trans (pow_le_pow_right₀ ha₁ (hm g)))
        (pow_nonneg (inv_nonneg.mpr (hV0 g).le) _) (by positivity)
  · rw [prime_bank_product_mass_eq_zero_of_unrepresented P k hz, mul_zero]
    exact pow_nonneg (mul_nonneg (by positivity) (pow_nonneg (zero_le_one.trans ha₁) _)) _

/-- The product law is bounded by L to the explicit graph exponent,
divided by the numerical product. -/
theorem prime_bank_product_mass_log_bound (hp : ∀ g p, p ∈ P g → p.Prime)
    (hdis : Pairwise fun g g' => Disjoint (P g) (P g')) {v : ℝ} (hv : 0 < v)
    (hV : ∀ g, v ≤ primeBankDrawMass (P g)) {m : ℕ} (hm : ∀ g, k g ≤ m)
    {L C : ℝ} (hL : 1 ≤ L) (hC : (Fintype.card G : ℝ) ≤ C * Real.log L) {n : ℕ} (hn : 0 < n) :
    primeBankProductMass P k n ≤
      L ^ (C * (Real.log m.factorial + m * Real.log (max 1 v⁻¹))) / (n : ℝ) := by
  let a := max (1 : ℝ) v⁻¹
  let H := (m.factorial : ℝ) * a ^ m
  have ha₁ : 1 ≤ a := le_max_left ..
  have ha : 0 < a := zero_lt_one.trans_le ha₁
  have hF : (0 : ℝ) < m.factorial := by exact_mod_cast Nat.factorial_pos m
  have hF₁ : (1 : ℝ) ≤ m.factorial := by exact_mod_cast Nat.succ_le_of_lt (Nat.factorial_pos m)
  have hH : 0 < H := mul_pos hF (pow_pos ha _)
  have hH₁ : 1 ≤ H := one_le_mul_of_one_le_of_one_le hF₁ (one_le_pow₀ ha₁)
  have hb : H ^ Fintype.card G ≤ L ^ (C * (Real.log m.factorial + m * Real.log a)) := by
    have he : Real.log H = Real.log m.factorial + m * Real.log a := by
      rw [show H = (m.factorial : ℝ) * a ^ m from rfl,
        Real.log_mul hF.ne' (pow_ne_zero _ ha.ne'), Real.log_pow]
    rw [← Real.exp_log (pow_pos hH _), Real.log_pow, Real.rpow_def_of_pos (zero_lt_one.trans_le hL)]
    apply Real.exp_le_exp.mpr
    refine (mul_le_mul_of_nonneg_right hC (Real.log_nonneg hH₁)).trans_eq ?_
    rw [he]
    ring
  apply (le_div_iff₀ (show (0 : ℝ) < n by exact_mod_cast hn)).mpr
  rw [mul_comm (primeBankProductMass P k n)]
  exact (prime_bank_product_mass_mul_le hp hdis hv hV hm n).trans hb

/-- Coefficient obtained by summing a tuple test over one product. -/
noncomputable def primeBankProductCoefficient (f : PrimeBankDraws P k → ℂ) (n : ℕ) : ℂ :=
  ∑ z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z = n),
    (primeBankDrawWeight z : ℂ) * f z

/-- Cauchy–Schwarz on a product fiber, with its probability mass as
the exact factor in the squared-norm bound. -/
theorem prime_bank_product_coefficient_sq_le (f : PrimeBankDraws P k → ℂ) (n : ℕ) :
    ‖primeBankProductCoefficient f n‖ ^ 2 ≤
      primeBankProductMass P k n *
        ∑ z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z = n),
          primeBankDrawWeight z * ‖f z‖ ^ 2 := by
  let t := Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z = n)
  let w (z : PrimeBankDraws P k) := primeBankDrawWeight z
  have hw z : 0 ≤ w z := prime_bank_draw_weight_nonneg z
  have h₁ : ‖primeBankProductCoefficient f n‖ ≤ ∑ z ∈ t, w z * ‖f z‖ := by
    unfold primeBankProductCoefficient
    refine (norm_sum_le _ _).trans_eq ?_
    exact Finset.sum_congr rfl (fun z _ => by rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg (hw z)])
  have h₂ := Finset.sum_mul_sq_le_sq_mul_sq t (fun z => Real.sqrt (w z))
    (fun z => Real.sqrt (w z) * ‖f z‖)
  simp only [← mul_assoc, ← pow_two (Real.sqrt (w _)), mul_pow, Real.sq_sqrt (hw _)] at h₂
  rw [primeBankProductMass, ← Finset.sum_filter]
  exact (pow_le_pow_left₀ (norm_nonneg _) h₁ 2).trans h₂

/-- A uniform mass bound on the product window bounds the pushed
forward energy by the tuple energy. -/
theorem prime_bank_product_coefficient_energy_le (f : PrimeBankDraws P k → ℂ)
    (s : Finset ℕ) {C : ℝ} (hC : 0 ≤ C) (hmass : ∀ n ∈ s, primeBankProductMass P k n ≤ C) :
    (∑ n ∈ s, ‖primeBankProductCoefficient f n‖ ^ 2) ≤
      C * ∑ z : PrimeBankDraws P k, primeBankDrawWeight z * ‖f z‖ ^ 2 := by
  let a (z : PrimeBankDraws P k) := primeBankDrawWeight z * ‖f z‖ ^ 2
  have ha z : 0 ≤ a z := mul_nonneg (prime_bank_draw_weight_nonneg z) (sq_nonneg _)
  calc
    _ ≤ ∑ n ∈ s, C * ∑ z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z = n),
        a z := by
      exact Finset.sum_le_sum (fun n hn => (prime_bank_product_coefficient_sq_le f n).trans
        (mul_le_mul_of_nonneg_right (hmass n hn) (Finset.sum_nonneg (fun z _ => ha z))))
    _ = C * ∑ z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z ∈ s), a z := by
      rw [← Finset.mul_sum, Finset.sum_fiberwise_eq_sum_filter]
    _ ≤ _ := mul_le_mul_of_nonneg_left (Finset.sum_le_sum_of_subset_of_nonneg
      (Finset.filter_subset _ _) (fun z _ _ => ha z)) hC

/-- The tuple-to-product energy bound on a product range whose
numerical products are at least U, with the explicit graph exponent. -/
theorem prime_bank_product_dyadic_energy_le (hp : ∀ g p, p ∈ P g → p.Prime)
    (hdis : Pairwise fun g g' => Disjoint (P g) (P g')) {v : ℝ} (hv : 0 < v)
    (hV : ∀ g, v ≤ primeBankDrawMass (P g)) {m : ℕ} (hm : ∀ g, k g ≤ m)
    {L C U : ℝ} (hL : 1 ≤ L) (hC : (Fintype.card G : ℝ) ≤ C * Real.log L)
    (hU : 0 < U) (f : PrimeBankDraws P k → ℂ) (s : Finset ℕ) (hs : ∀ n ∈ s, U ≤ (n : ℝ)) :
    (∑ n ∈ s, ‖primeBankProductCoefficient f n‖ ^ 2) ≤
      (L ^ (C * (Real.log m.factorial + m * Real.log (max 1 v⁻¹))) / U) *
        ∑ z : PrimeBankDraws P k, primeBankDrawWeight z * ‖f z‖ ^ 2 := by
  have hL₀ : 0 < L := zero_lt_one.trans_le hL
  apply prime_bank_product_coefficient_energy_le f s (by positivity)
  intro n hn
  have hn₀ : 0 < n := by exact_mod_cast hU.trans_le (hs n hn)
  exact (prime_bank_product_mass_log_bound hp hdis hv hV hm hL hC hn₀).trans
    (div_le_div_of_nonneg_left (by positivity) hU (hs n hn))

/-- Conditional mean of a tuple test on a product fiber. Empty
fibers have zero mean. -/
noncomputable def primeBankConditionalMean (f : PrimeBankDraws P k → ℂ) (n : ℕ) : ℂ :=
  primeBankProductCoefficient f n / (primeBankProductMass P k n : ℂ)

/-- Conditional averaging from ordered tuples to numerical products
is a contraction for the probability-weighted squared norm. -/
theorem prime_bank_product_conditional_mean_energy_le (f : PrimeBankDraws P k → ℂ) (s : Finset ℕ) :
    (∑ n ∈ s, primeBankProductMass P k n * ‖primeBankConditionalMean f n‖ ^ 2) ≤
      ∑ z : PrimeBankDraws P k, primeBankDrawWeight z * ‖f z‖ ^ 2 := by
  let a (z : PrimeBankDraws P k) := primeBankDrawWeight z * ‖f z‖ ^ 2
  have ha z : 0 ≤ a z := mul_nonneg (prime_bank_draw_weight_nonneg z) (sq_nonneg _)
  have hpoint (n) : (primeBankProductMass P k n) * ‖primeBankConditionalMean f n‖ ^ 2 ≤
      ∑ z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z = n), a z := by
    let p := primeBankProductMass P k n
    have hnonneg : 0 ≤ p := by
      unfold p primeBankProductMass
      exact Finset.sum_nonneg (fun z _ => by split; exact prime_bank_draw_weight_nonneg z; exact le_rfl)
    change p * ‖primeBankConditionalMean f n‖ ^ 2 ≤ _
    by_cases hp : p = 0
    · rw [hp, zero_mul]
      exact Finset.sum_nonneg (fun z _ => ha z)
    · have hp₀ : 0 < p := lt_of_le_of_ne hnonneg (Ne.symm hp)
      unfold primeBankConditionalMean
      rw [norm_div, Complex.norm_real, Real.norm_of_nonneg hnonneg, div_pow]
      change p * (‖primeBankProductCoefficient f n‖ ^ 2 / p ^ 2) ≤ _
      have he : p * (‖primeBankProductCoefficient f n‖ ^ 2 / p ^ 2) =
          ‖primeBankProductCoefficient f n‖ ^ 2 / p := by field_simp
      rw [he, div_le_iff₀ hp₀, mul_comm]
      exact prime_bank_product_coefficient_sq_le f n
  calc
    _ ≤ ∑ n ∈ s, ∑ z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z = n), a z :=
      Finset.sum_le_sum (fun n _ => hpoint n)
    _ = ∑ z ∈ Finset.univ.filter (fun z : PrimeBankDraws P k => primeBankDrawProduct z ∈ s), a z :=
      Finset.sum_fiberwise_eq_sum_filter Finset.univ s primeBankDrawProduct a
    _ ≤ _ := Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) (fun z _ _ => ha z)

end WeightBounds

end QuadraticNorms
