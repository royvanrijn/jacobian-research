import Mathlib.Algebra.BigOperators.Field
import Mathlib.Tactic.DeriveFintype
import Mathlib.Data.Nat.ChineseRemainder
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Analysis.SpecialFunctions.Pow.Real

open scoped BigOperators

namespace QuadraticNorms

/-- Exact valuation at a prime, read from the canonical representative.
At modulus p^k the predicate is residue-invariant for exponents below k. -/
abbrev BankExactValuation {n : ℕ} (p : ℕ) (z : ZMod n) (e : ℕ) : Prop :=
  p ^ e ∣ z.val ∧ ¬ p ^ (e + 1) ∣ z.val

private theorem filter_precomp_equiv_card {α β : Type*} [Fintype α] [Fintype β]
    (g : α ≃ β) (P : β → Prop) [DecidablePred P] :
    (Finset.univ.filter (fun a => P (g a))).card = (Finset.univ.filter P).card := by
  refine Finset.card_bij' (fun a _ => g a) (fun b _ => g.symm b) ?_ ?_ ?_ ?_
  · intro a ha
    simpa only [Finset.mem_filter, Finset.mem_univ, true_and] using ha
  · intro b hb
    simpa only [Finset.mem_filter, Finset.mem_univ, true_and, Equiv.apply_symm_apply] using hb
  · exact fun a _ => g.symm_apply_apply a
  · exact fun b _ => g.apply_symm_apply b

private theorem zmod_val_divisible_card {n d : ℕ} [NeZero n] (hd : d ∣ n) :
    (Finset.univ.filter (fun x : ZMod n => d ∣ x.val)).card = n / d := by
  have hn := Nat.pos_of_ne_zero (NeZero.ne n)
  have hd₀ := Nat.pos_of_dvd_of_pos hd hn
  rw [← Finset.card_range (n / d)]
  refine Finset.card_bij' (fun x _ => x.val / d) (fun i _ => ((d * i : ℕ) : ZMod n)) ?_ ?_ ?_ ?_
  · intro x _
    apply Finset.mem_range.mpr
    apply (Nat.div_lt_iff_lt_mul hd₀).mpr
    rw [Nat.div_mul_cancel hd]
    exact ZMod.val_lt x
  · intro i hi
    have hil : d * i < n := by
      rw [← Nat.mul_div_cancel' hd]
      exact Nat.mul_lt_mul_of_pos_left (Finset.mem_range.mp hi) hd₀
    exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, by rw [ZMod.val_natCast_of_lt hil]; exact dvd_mul_right ..⟩
  · intro x hx
    rw [Nat.mul_div_cancel' (Finset.mem_filter.mp hx).2, ZMod.natCast_zmod_val]
  · intro i hi
    have hil : d * i < n := by
      rw [← Nat.mul_div_cancel' hd]
      exact Nat.mul_lt_mul_of_pos_left (Finset.mem_range.mp hi) hd₀
    rw [ZMod.val_natCast_of_lt hil, Nat.mul_div_cancel_left _ hd₀]

private theorem prime_power_divisible_card {p k j : ℕ} (hp : 0 < p) [NeZero (p ^ k)]
    (hj : j ≤ k) :
    (Finset.univ.filter (fun x : ZMod (p ^ k) => p ^ j ∣ x.val)).card = p ^ (k - j) := by
  rw [zmod_val_divisible_card (Nat.pow_dvd_pow p hj), Nat.pow_div hj hp]

private theorem bank_exact_valuation_card {p k e : ℕ} (hp : 0 < p) [NeZero (p ^ k)]
    (he : e < k) :
    (Finset.univ.filter (fun x : ZMod (p ^ k) => BankExactValuation p x e)).card =
      p ^ (k - e) - p ^ (k - (e + 1)) := by
  let S := Finset.univ.filter (fun x : ZMod (p ^ k) => p ^ e ∣ x.val)
  let T := Finset.univ.filter (fun x : ZMod (p ^ k) => p ^ (e + 1) ∣ x.val)
  have heq : Finset.univ.filter (fun x : ZMod (p ^ k) => BankExactValuation p x e) = S \ T := by
    ext x
    simp [S, T]
  have hsub : T ⊆ S := fun x hx => Finset.mem_filter.mpr
    ⟨Finset.mem_univ x, (Nat.pow_dvd_pow p (Nat.le_succ e)).trans (Finset.mem_filter.mp hx).2⟩
  rw [heq, Finset.card_sdiff_of_subset hsub]
  change (Finset.univ.filter (fun x : ZMod (p ^ k) => p ^ e ∣ x.val)).card -
    (Finset.univ.filter (fun x : ZMod (p ^ k) => p ^ (e + 1) ∣ x.val)).card = _
  rw [prime_power_divisible_card hp he.le, prime_power_divisible_card hp he]

private theorem bank_exact_valuation_zero_iff {p n : ℕ} (x : ZMod n) :
    BankExactValuation p x 0 ↔ ¬ p ∣ x.val := by
  simp only [BankExactValuation, pow_zero, zero_add, pow_one, one_dvd, true_and]

private theorem bank_exact_valuation_pos_dvd {p n e : ℕ} {x : ZMod n} (he : 0 < e)
    (hx : BankExactValuation p x e) : p ∣ x.val :=
  (dvd_pow_self p he.ne').trans hx.1

private theorem bank_cast_eq_zero_iff {p k : ℕ} [NeZero (p ^ k)] (hk : 0 < k)
    (z : ZMod (p ^ k)) :
    ZMod.castHom (dvd_pow_self p hk.ne') (ZMod p) z = 0 ↔ p ∣ z.val := by
  let f := ZMod.castHom (dvd_pow_self p hk.ne') (ZMod p)
  change f z = 0 ↔ _
  have heq : f z = (z.val : ZMod p) := by rw [← map_natCast f z.val, ZMod.natCast_zmod_val]
  rw [heq, ZMod.natCast_eq_zero_iff]

private theorem no_common_bank_prime {p k : ℕ} [NeZero (p ^ k)] (hk : 0 < k)
    (u : (ZMod (p ^ k))ˣ) {d : ZMod (p ^ k)} (hd : ¬ p ∣ d.val) (x : ZMod (p ^ k)) :
    ¬ (p ∣ x.val ∧ p ∣ ((u : ZMod (p ^ k)) * x + d).val) := by
  let f := ZMod.castHom (dvd_pow_self p hk.ne') (ZMod p)
  have hf (y : ZMod (p ^ k)) : f y = 0 ↔ p ∣ y.val := bank_cast_eq_zero_iff hk y
  rintro ⟨hx, hy⟩
  have h := (hf _).mpr hy
  rw [map_add, map_mul, (hf x).mpr hx, mul_zero, zero_add] at h
  exact hd ((hf d).mp h)

/-- The local law at a single bank prime. -/
noncomputable def bankPrimeValuationWeight (p e f : ℕ) : ℝ :=
  if e = 0 then
    if f = 0 then 1 - 2 / (p : ℝ) else (1 - 1 / (p : ℝ)) / (p : ℝ) ^ f
  else
    if f = 0 then (1 - 1 / (p : ℝ)) / (p : ℝ) ^ e else 0

/-- The local cell in the normalized line parameter x, with second
coordinate u*x+d. -/
def bankValuationCell {n : ℕ} [NeZero n] (p : ℕ) (u : (ZMod n)ˣ) (d : ZMod n)
    (e f : ℕ) : Finset (ZMod n) :=
  Finset.univ.filter (fun x => BankExactValuation p x e ∧ BankExactValuation p ((u : ZMod n) * x + d) f)

private theorem bank_valuation_cell_card_zero_zero {p k : ℕ} (hp : 0 < p) (hk : 0 < k)
    [NeZero (p ^ k)] (u : (ZMod (p ^ k))ˣ) {d : ZMod (p ^ k)} (hd : ¬ p ∣ d.val) :
    (bankValuationCell p u d 0 0).card = p ^ k - 2 * p ^ (k - 1) := by
  let g := (Units.mulLeft u).trans (Equiv.addRight d)
  let S := Finset.univ.filter (fun x : ZMod (p ^ k) => p ∣ x.val)
  let T := Finset.univ.filter (fun x : ZMod (p ^ k) => p ∣ (g x).val)
  have hdis : Disjoint S T := by
    rw [Finset.disjoint_left]
    exact fun x hx hy => no_common_bank_prime hk u hd x
      ⟨(Finset.mem_filter.mp hx).2, (Finset.mem_filter.mp hy).2⟩
  have heq : bankValuationCell p u d 0 0 = Finset.univ \ (S ∪ T) := by
    ext x
    simp only [bankValuationCell, bank_exact_valuation_zero_iff, Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.mem_sdiff, Finset.mem_union, not_or, S, T, g]
    rfl
  have hS : S.card = p ^ (k - 1) := by
    simpa only [pow_one, S] using prime_power_divisible_card (p := p) hp hk
  have hT : T.card = S.card := filter_precomp_equiv_card g (fun x => p ∣ x.val)
  rw [heq, Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_union_of_disjoint hdis,
    hT, hS, Finset.card_univ, ZMod.card, two_mul]

private theorem bank_valuation_cell_card_pos_zero {p k e : ℕ} (hp : 0 < p) (he : 0 < e)
    (hek : e < k) [NeZero (p ^ k)] (u : (ZMod (p ^ k))ˣ) {d : ZMod (p ^ k)} (hd : ¬ p ∣ d.val) :
    (bankValuationCell p u d e 0).card = p ^ (k - e) - p ^ (k - (e + 1)) := by
  have heq : bankValuationCell p u d e 0 =
      Finset.univ.filter (fun x : ZMod (p ^ k) => BankExactValuation p x e) := by
    ext x
    simp only [bankValuationCell, bank_exact_valuation_zero_iff, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨And.left, fun hx => ⟨hx, fun hy =>
      no_common_bank_prime (by omega) u hd x ⟨bank_exact_valuation_pos_dvd he hx, hy⟩⟩⟩
  rw [heq, bank_exact_valuation_card hp hek]

private theorem bank_valuation_cell_card_zero_pos {p k f : ℕ} (hp : 0 < p) (hf : 0 < f)
    (hfk : f < k) [NeZero (p ^ k)] (u : (ZMod (p ^ k))ˣ) {d : ZMod (p ^ k)} (hd : ¬ p ∣ d.val) :
    (bankValuationCell p u d 0 f).card = p ^ (k - f) - p ^ (k - (f + 1)) := by
  let g := (Units.mulLeft u).trans (Equiv.addRight d)
  have heq : bankValuationCell p u d 0 f =
      Finset.univ.filter (fun x : ZMod (p ^ k) => BankExactValuation p (g x) f) := by
    ext x
    simp only [bankValuationCell, bank_exact_valuation_zero_iff, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨And.right, fun hy => ⟨fun hx =>
      no_common_bank_prime (by omega) u hd x ⟨hx, bank_exact_valuation_pos_dvd hf hy⟩, hy⟩⟩
  rw [heq, filter_precomp_equiv_card g (fun y => BankExactValuation p y f), bank_exact_valuation_card hp hfk]

private theorem prime_power_ratio {p j k : ℕ} (hp : 0 < p) (hj : j ≤ k) :
    ((p ^ (k - j) : ℕ) : ℝ) / (p : ℝ) ^ k = 1 / (p : ℝ) ^ j := by
  have hpR : 0 < (p : ℝ) := by exact_mod_cast hp
  rw [Nat.cast_pow]
  field_simp
  rw [← pow_add, Nat.sub_add_cancel hj]

private theorem normalized_single_valuation {p k e : ℕ} (hp : 0 < p) (he : e < k) :
    ((p ^ (k - e) - p ^ (k - (e + 1)) : ℕ) : ℝ) / (p : ℝ) ^ k =
      (1 - 1 / (p : ℝ)) / (p : ℝ) ^ e := by
  have hpR : 0 < (p : ℝ) := by exact_mod_cast hp
  have hle : p ^ (k - (e + 1)) ≤ p ^ (k - e) := Nat.pow_le_pow_right hp (by omega)
  rw [Nat.cast_sub hle, sub_div, prime_power_ratio hp he.le, prime_power_ratio hp he, pow_succ]
  field_simp

/-- The exact normalized record count at modulus p^k. The two requested
exponents are smaller than k, so neither is a saturated valuation. -/
theorem bank_prime_valuation_law {p k e f : ℕ} (hp : p.Prime) (he : e < k) (hf : f < k)
    [NeZero (p ^ k)] (u : (ZMod (p ^ k))ˣ) (d : ZMod (p ^ k)) (hd : ¬ p ∣ d.val) :
    ((bankValuationCell p u d e f).card : ℝ) / (p : ℝ) ^ k = bankPrimeValuationWeight p e f := by
  have hk : 0 < k := by omega
  by_cases he₀ : e = 0
  · subst e
    by_cases hf₀ : f = 0
    · subst f
      have hp₀ : 0 < (p : ℝ) := by exact_mod_cast hp.pos
      have hle : 2 * p ^ (k - 1) ≤ p ^ k := by
        calc
          _ ≤ p * p ^ (k - 1) := Nat.mul_le_mul_right _ hp.two_le
          _ = p ^ k := by rw [← pow_succ', Nat.sub_add_cancel hk]
      rw [bank_valuation_cell_card_zero_zero hp.pos hk u hd, Nat.cast_sub hle, sub_div,
        Nat.cast_pow, div_self (by positivity : (p : ℝ) ^ k ≠ 0), Nat.cast_mul, Nat.cast_ofNat,
        mul_div_assoc, prime_power_ratio hp.pos hk]
      simp [bankPrimeValuationWeight, div_eq_mul_inv]
    · rw [bank_valuation_cell_card_zero_pos hp.pos (Nat.pos_of_ne_zero hf₀) hf u hd,
        normalized_single_valuation hp.pos hf]
      simp [bankPrimeValuationWeight, hf₀]
  · by_cases hf₀ : f = 0
    · subst f
      rw [bank_valuation_cell_card_pos_zero hp.pos (Nat.pos_of_ne_zero he₀) he u hd,
        normalized_single_valuation hp.pos he]
      simp [bankPrimeValuationWeight, he₀]
    · have hempty : bankValuationCell p u d e f = ∅ := by
        apply Finset.eq_empty_iff_forall_notMem.mpr
        intro x hx
        obtain ⟨ha, hb⟩ := (Finset.mem_filter.mp hx).2
        exact no_common_bank_prime hk u hd x ⟨bank_exact_valuation_pos_dvd (Nat.pos_of_ne_zero he₀) ha,
          bank_exact_valuation_pos_dvd (Nat.pos_of_ne_zero hf₀) hb⟩
      simp [hempty, bankPrimeValuationWeight, he₀, hf₀]

/-- The valuation law on the fixed-shift line. -/
theorem bank_prime_shift_valuation_law {p k e f : ℕ} (hp : p.Prime) (he : e < k) (hf : f < k)
    [NeZero (p ^ k)] (d : ZMod (p ^ k)) (hd : ¬ p ∣ d.val) :
    ((Finset.univ.filter (fun x : ZMod (p ^ k) =>
      BankExactValuation p x e ∧ BankExactValuation p (x + d) f)).card : ℝ) / (p : ℝ) ^ k =
      bankPrimeValuationWeight p e f := by
  simpa only [bankValuationCell, Units.val_one, one_mul] using bank_prime_valuation_law hp he hf 1 d hd

/-- The valuation law on the reflected sum line with unit
coefficients and a target prime to the bank prime. -/
theorem bank_prime_sum_valuation_law {p k e f : ℕ} (hp : p.Prime) (he : e < k) (hf : f < k)
    [NeZero (p ^ k)] (A B : (ZMod (p ^ k))ˣ) (N : ZMod (p ^ k)) (hN : ¬ p ∣ N.val) :
    ((Finset.univ.filter (fun z : ZMod (p ^ k) × ZMod (p ^ k) =>
      (A : ZMod (p ^ k)) * z.1 + (B : ZMod (p ^ k)) * z.2 = N ∧
        BankExactValuation p z.1 e ∧ BankExactValuation p z.2 f)).card : ℝ) / (p : ℝ) ^ k =
      bankPrimeValuationWeight p e f := by
  let U := -A * B⁻¹
  let D := (↑B⁻¹ : ZMod (p ^ k)) * N
  let g (x : ZMod (p ^ k)) := (U : ZMod (p ^ k)) * x + D
  have hg (x : ZMod (p ^ k)) : g x = (↑B⁻¹ : ZMod (p ^ k)) * (N - A * x) := by
    simp only [g, U, D, Units.val_mul, Units.val_neg]
    ring
  have hline (x y : ZMod (p ^ k)) : A * x + B * y = N ↔ y = g x := by
    rw [hg]
    constructor
    · intro h
      calc
        y = (↑B⁻¹ : ZMod (p ^ k)) * (B * y) := by rw [← mul_assoc, Units.inv_mul, one_mul]
        _ = _ := congrArg ((↑B⁻¹ : ZMod (p ^ k)) * ·) (eq_sub_of_add_eq' h)
    · intro h
      apply eq_sub_iff_add_eq'.mp
      rw [h, ← mul_assoc, Units.mul_inv, one_mul]
  have hD : ¬ p ∣ D.val := by
    have hk : 0 < k := by omega
    let reduce := ZMod.castHom (dvd_pow_self p hk.ne') (ZMod p)
    intro h
    have hzero : reduce D = 0 := (bank_cast_eq_zero_iff hk D).mpr h
    change reduce ((↑B⁻¹ : ZMod (p ^ k)) * N) = 0 at hzero
    rw [map_mul] at hzero
    have hunit : IsUnit (reduce (↑B⁻¹ : ZMod (p ^ k))) := (B⁻¹).isUnit.map reduce
    have hz : reduce N = 0 := hunit.mul_left_cancel (hzero.trans (mul_zero _).symm)
    exact hN ((bank_cast_eq_zero_iff hk N).mp hz)
  have hcard : (Finset.univ.filter (fun z : ZMod (p ^ k) × ZMod (p ^ k) =>
      (A : ZMod (p ^ k)) * z.1 + (B : ZMod (p ^ k)) * z.2 = N ∧
        BankExactValuation p z.1 e ∧ BankExactValuation p z.2 f)).card =
      (bankValuationCell p U D e f).card := by
    refine Finset.card_bij' (fun z _ => z.1) (fun x _ => (x, g x)) ?_ ?_ ?_ (fun _ _ => rfl)
    · intro z hz
      obtain ⟨hl, hleft, hright⟩ := (Finset.mem_filter.mp hz).2
      have heq := (hline z.1 z.2).mp hl
      have hv : BankExactValuation p (g z.1) f := by rw [← heq]; exact hright
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, hleft, hv⟩
    · intro x hx
      have hv := (Finset.mem_filter.mp hx).2
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, (hline x (g x)).mpr rfl, hv⟩
    · intro z hz
      exact Prod.ext rfl ((hline z.1 z.2).mp (Finset.mem_filter.mp hz).2.1).symm
  rw [hcard]
  exact bank_prime_valuation_law hp he hf U D hD

/-- The three divisibility outcomes at a bank prime. -/
inductive BankPrimeOutcome
  | neither
  | left
  | right
  deriving DecidableEq, Fintype

/-- The distribution after forgetting the exact positive valuation. -/
noncomputable def bankPrimeOutcomeWeight (r : ℝ) : BankPrimeOutcome → ℝ
  | .neither => 1 - 2 * r
  | .left => r
  | .right => r

def bankPrimeOutcome {n : ℕ} (p : ℕ) (u : (ZMod n)ˣ) (d x : ZMod n) : BankPrimeOutcome :=
  if p ∣ x.val then .left else if p ∣ ((u : ZMod n) * x + d).val then .right else .neither

def bankOutcomeCell {n : ℕ} [NeZero n] (p : ℕ) (u : (ZMod n)ˣ) (d : ZMod n)
    (s : BankPrimeOutcome) : Finset (ZMod n) :=
  Finset.univ.filter (fun x => bankPrimeOutcome p u d x = s)

theorem bank_prime_outcome_law {p k : ℕ} (hp : p.Prime) (hk : 0 < k) [NeZero (p ^ k)]
    (u : (ZMod (p ^ k))ˣ) (d : ZMod (p ^ k)) (hd : ¬ p ∣ d.val) (s : BankPrimeOutcome) :
    ((bankOutcomeCell p u d s).card : ℝ) / (p : ℝ) ^ k =
      bankPrimeOutcomeWeight (1 / (p : ℝ)) s := by
  cases s with
  | neither =>
    have heq : bankOutcomeCell p u d .neither = bankValuationCell p u d 0 0 := by
      ext x
      simp only [bankOutcomeCell, bankValuationCell, bank_exact_valuation_zero_iff,
        Finset.mem_filter, Finset.mem_univ, true_and, bankPrimeOutcome]
      by_cases h : p ∣ x.val <;> by_cases h' : p ∣ ((u : ZMod (p ^ k)) * x + d).val <;> simp [h, h']
    rw [heq, bank_prime_valuation_law hp hk hk u d hd]
    simp [bankPrimeOutcomeWeight, bankPrimeValuationWeight, div_eq_mul_inv]
  | left =>
    have heq : bankOutcomeCell p u d .left = Finset.univ.filter (fun x => p ∣ x.val) := by
      ext x
      simp only [bankOutcomeCell, Finset.mem_filter, Finset.mem_univ, true_and, bankPrimeOutcome]
      by_cases h : p ∣ x.val <;> by_cases h' : p ∣ ((u : ZMod (p ^ k)) * x + d).val <;> simp [h, h']
    have h : (Finset.univ.filter (fun x : ZMod (p ^ k) => p ∣ x.val)).card = p ^ (k - 1) := by
      simpa only [pow_one] using prime_power_divisible_card hp.pos hk
    rw [heq, h, prime_power_ratio hp.pos hk]
    simp [bankPrimeOutcomeWeight]
  | right =>
    let g := (Units.mulLeft u).trans (Equiv.addRight d)
    have heq : bankOutcomeCell p u d .right = Finset.univ.filter (fun x => p ∣ (g x).val) := by
      ext x
      simp only [bankOutcomeCell, Finset.mem_filter, Finset.mem_univ, true_and, bankPrimeOutcome]
      change (if p ∣ x.val then BankPrimeOutcome.left else
        if p ∣ (g x).val then .right else .neither) = .right ↔ p ∣ (g x).val
      by_cases h : p ∣ x.val
      · have h' : ¬ p ∣ (g x).val := fun h' => no_common_bank_prime hk u hd x ⟨h, h'⟩
        simp [h, h']
      · by_cases h' : p ∣ (g x).val <;> simp [h, h']
    rw [heq, filter_precomp_equiv_card g (fun x => p ∣ x.val)]
    have h : (Finset.univ.filter (fun x : ZMod (p ^ k) => p ∣ x.val)).card = p ^ (k - 1) := by
      simpa only [pow_one] using prime_power_divisible_card hp.pos hk
    rw [h, prime_power_ratio hp.pos hk]
    simp [bankPrimeOutcomeWeight]

section CRT

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Reduction at the i-th modulus from the product modulus. -/
def bankCRTCoordinate (m : ι → ℕ) (x : ZMod (∏ i, m i)) (i : ι) : ZMod (m i) :=
  ZMod.castHom (Finset.dvd_prod_of_mem m (Finset.mem_univ i)) (ZMod (m i)) x

theorem bank_crt_coordinate_bijective (m : ι → ℕ) (hm : ∀ i, m i ≠ 0)
    (hcop : Pairwise (fun i j => (m i).Coprime (m j))) :
    Function.Bijective (bankCRTCoordinate m) := by
  letI : ∀ i, NeZero (m i) := fun i => ⟨hm i⟩
  letI : NeZero (∏ i, m i) := ⟨Finset.prod_ne_zero_iff.mpr (fun i _ => hm i)⟩
  apply (Fintype.bijective_iff_surjective_and_card _).mpr
  constructor
  · intro x
    let a := Nat.chineseRemainderOfFinset (fun i => (x i).val) m Finset.univ
      (fun i _ => hm i) (fun i _ j _ h => hcop h)
    refine ⟨(a.val : ZMod (∏ i, m i)), funext (fun i => ?_)⟩
    change (ZMod.castHom _ (ZMod (m i))) (a.val : ZMod (∏ i, m i)) = _
    rw [map_natCast]
    exact ((ZMod.natCast_eq_natCast_iff _ _ _).mpr (a.property i (Finset.mem_univ i))).trans
      (ZMod.natCast_zmod_val (x i))
  · rw [ZMod.card, Fintype.card_pi]
    simp only [ZMod.card]

/-- The simultaneous CRT cell for a family of local events. -/
noncomputable def bankCRTCell (m : ι → ℕ) [NeZero (∏ i, m i)]
    (cells : ∀ i, Finset (ZMod (m i))) : Finset (ZMod (∏ i, m i)) := by
  classical
  exact Finset.univ.filter (fun x => ∀ i, bankCRTCoordinate m x i ∈ cells i)

theorem bank_crt_cell_card (m : ι → ℕ) [NeZero (∏ i, m i)] (hm : ∀ i, m i ≠ 0)
    (hcop : Pairwise (fun i j => (m i).Coprime (m j))) (cells : ∀ i, Finset (ZMod (m i))) :
    (bankCRTCell m cells).card = ∏ i, (cells i).card := by
  classical
  let g := Equiv.ofBijective (bankCRTCoordinate m) (bank_crt_coordinate_bijective m hm hcop)
  have h := Fintype.card_piFinset cells
  rw [← h]
  refine Finset.card_bij' (fun x _ => g x) (fun y _ => g.symm y) ?_ ?_ ?_ ?_
  · intro x hx
    exact Fintype.mem_piFinset.mpr (Finset.mem_filter.mp hx).2
  · intro y hy
    apply Finset.mem_filter.mpr
    refine ⟨Finset.mem_univ _, ?_⟩
    change ∀ i, g (g.symm y) i ∈ cells i
    rw [Equiv.apply_symm_apply]
    exact Fintype.mem_piFinset.mp hy
  · exact fun x _ => g.symm_apply_apply x
  · exact fun y _ => g.apply_symm_apply y

theorem bank_crt_cell_proportion (m : ι → ℕ) [NeZero (∏ i, m i)] (hm : ∀ i, m i ≠ 0)
    (hcop : Pairwise (fun i j => (m i).Coprime (m j))) (cells : ∀ i, Finset (ZMod (m i))) :
    ((bankCRTCell m cells).card : ℝ) / (∏ i, (m i : ℝ)) =
      ∏ i, ((cells i).card : ℝ) / (m i : ℝ) := by
  rw [bank_crt_cell_card m hm hcop, Nat.cast_prod, Finset.prod_div_distrib]

omit [Fintype ι] [DecidableEq ι] in
private theorem bank_prime_power_pairwise {p : ι → ℕ} (hp : ∀ i, (p i).Prime)
    (hDistinct : Function.Injective p) (k : ι → ℕ) :
    Pairwise (fun i j => (p i ^ k i).Coprime (p j ^ k j)) :=
  fun i j h => ((Nat.coprime_primes (hp i) (hp j)).mpr (fun he => h (hDistinct he))).pow _ _

/-- The exact product law for prime valuations on the joint residue
space. Every local exponent is taken below the corresponding modulus
exponent. -/
theorem bank_joint_prime_valuation_law (p k e f : ι → ℕ)
    (hp : ∀ i, (p i).Prime) (hDistinct : Function.Injective p)
    (he : ∀ i, e i < k i) (hf : ∀ i, f i < k i)
    (u : ∀ i, (ZMod (p i ^ k i))ˣ) (d : ∀ i, ZMod (p i ^ k i))
    (hd : ∀ i, ¬ p i ∣ (d i).val) :
    letI : ∀ i, NeZero (p i ^ k i) := fun i => ⟨pow_ne_zero _ (hp i).ne_zero⟩
    letI : NeZero (∏ i, p i ^ k i) := ⟨Finset.prod_ne_zero_iff.mpr
      (fun i _ => pow_ne_zero _ (hp i).ne_zero)⟩
    ((bankCRTCell (fun i => p i ^ k i) (fun i => bankValuationCell (p i) (u i) (d i) (e i) (f i))).card : ℝ) /
        (∏ i, (p i : ℝ) ^ k i) = ∏ i, bankPrimeValuationWeight (p i) (e i) (f i) := by
  let m := fun i => p i ^ k i
  letI : ∀ i, NeZero (m i) := fun i => ⟨pow_ne_zero _ (hp i).ne_zero⟩
  letI : NeZero (∏ i, m i) := ⟨Finset.prod_ne_zero_iff.mpr (fun i _ => NeZero.ne (m i))⟩
  change _ / (∏ i, (p i : ℝ) ^ k i) = _
  simp_rw [← Nat.cast_pow]
  rw [bank_crt_cell_proportion m (fun i => NeZero.ne (m i)) (bank_prime_power_pairwise hp hDistinct k)]
  apply Finset.prod_congr rfl
  intro i _
  rw [Nat.cast_pow]
  exact bank_prime_valuation_law (hp i) (he i) (hf i) (u i) (d i) (hd i)

/-- The product law after keeping only the two divisibility indicators. -/
theorem bank_joint_prime_outcome_law (p k : ι → ℕ) (s : ι → BankPrimeOutcome)
    (hp : ∀ i, (p i).Prime) (hDistinct : Function.Injective p) (hk : ∀ i, 0 < k i)
    (u : ∀ i, (ZMod (p i ^ k i))ˣ) (d : ∀ i, ZMod (p i ^ k i))
    (hd : ∀ i, ¬ p i ∣ (d i).val) :
    letI : ∀ i, NeZero (p i ^ k i) := fun i => ⟨pow_ne_zero _ (hp i).ne_zero⟩
    letI : NeZero (∏ i, p i ^ k i) := ⟨Finset.prod_ne_zero_iff.mpr
      (fun i _ => pow_ne_zero _ (hp i).ne_zero)⟩
    ((bankCRTCell (fun i => p i ^ k i) (fun i => bankOutcomeCell (p i) (u i) (d i) (s i))).card : ℝ) /
        (∏ i, (p i : ℝ) ^ k i) = ∏ i, bankPrimeOutcomeWeight (1 / (p i : ℝ)) (s i) := by
  let m := fun i => p i ^ k i
  letI : ∀ i, NeZero (m i) := fun i => ⟨pow_ne_zero _ (hp i).ne_zero⟩
  letI : NeZero (∏ i, m i) := ⟨Finset.prod_ne_zero_iff.mpr (fun i _ => NeZero.ne (m i))⟩
  change _ / (∏ i, (p i : ℝ) ^ k i) = _
  simp_rw [← Nat.cast_pow]
  rw [bank_crt_cell_proportion m (fun i => NeZero.ne (m i)) (bank_prime_power_pairwise hp hDistinct k)]
  apply Finset.prod_congr rfl
  intro i _
  rw [Nat.cast_pow]
  exact bank_prime_outcome_law (hp i) (hk i) (u i) (d i) (hd i) (s i)

/-- The product distribution gives the exact proportion of every event
in the finite record space, by grouping the actual residues by record. -/
theorem bank_joint_outcome_event_proportion (p k : ι → ℕ)
    (hp : ∀ i, (p i).Prime) (hDistinct : Function.Injective p) (hk : ∀ i, 0 < k i)
    (u : ∀ i, (ZMod (p i ^ k i))ˣ) (d : ∀ i, ZMod (p i ^ k i))
    (hd : ∀ i, ¬ p i ∣ (d i).val) (event : (ι → BankPrimeOutcome) → Prop) [DecidablePred event] :
    letI : NeZero (∏ i, p i ^ k i) := ⟨Finset.prod_ne_zero_iff.mpr
      (fun i _ => pow_ne_zero _ (hp i).ne_zero)⟩
    ((Finset.univ.filter (fun x : ZMod (∏ i, p i ^ k i) =>
      event (fun i => bankPrimeOutcome (p i) (u i) (d i)
        (bankCRTCoordinate (fun j => p j ^ k j) x i)))).card : ℝ) / (∏ i, (p i : ℝ) ^ k i) =
      ∑ s, if event s then ∏ i, bankPrimeOutcomeWeight (1 / (p i : ℝ)) (s i) else 0 := by
  classical
  let m := fun i => p i ^ k i
  letI : ∀ i, NeZero (m i) := fun i => ⟨pow_ne_zero _ (hp i).ne_zero⟩
  letI : NeZero (∏ i, m i) := ⟨Finset.prod_ne_zero_iff.mpr (fun i _ => NeZero.ne (m i))⟩
  let record (x : ZMod (∏ i, m i)) (i : ι) :=
    bankPrimeOutcome (p i) (u i) (d i) (bankCRTCoordinate m x i)
  let region := Finset.univ.filter (fun x => event (record x))
  have hcard : region.card = ∑ s, (region.filter (fun x => record x = s)).card :=
    Finset.card_eq_sum_card_fiberwise (fun _ _ => Finset.mem_univ _)
  change (region.card : ℝ) / _ = _
  rw [hcard, Nat.cast_sum, Finset.sum_div]
  refine Finset.sum_congr rfl (fun s _ => ?_)
  let cells := fun i => bankOutcomeCell (p i) (u i) (d i) (s i)
  have hmem (x : ZMod (∏ i, m i)) : x ∈ bankCRTCell m cells ↔ record x = s := by
    simp only [bankCRTCell, bankOutcomeCell, record, cells, Finset.mem_filter,
      Finset.mem_univ, true_and, funext_iff]
  by_cases h : event s
  · rw [if_pos h]
    have heq : region.filter (fun x => record x = s) = bankCRTCell m cells := by
      ext x
      rw [hmem x]
      simp only [region, Finset.mem_filter, Finset.mem_univ, true_and]
      exact ⟨And.right, fun hx => ⟨hx.symm ▸ h, hx⟩⟩
    rw [heq]
    exact bank_joint_prime_outcome_law p k s hp hDistinct hk u d hd
  · rw [if_neg h]
    have heq : region.filter (fun x => record x = s) = ∅ := by
      apply Finset.eq_empty_iff_forall_notMem.mpr
      intro x hx
      have ⟨hr, he⟩ := Finset.mem_filter.mp hx
      exact h (he ▸ (Finset.mem_filter.mp hr).2)
    rw [heq, Finset.card_empty, Nat.cast_zero, zero_div]

end CRT

end QuadraticNorms
