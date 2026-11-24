import QuadraticNorms.NormRootCounts
import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.Data.Nat.Sqrt

open scoped BigOperators

namespace QuadraticNorms

private theorem isSquare_iff_even_factorization {n : ℕ} (hn : n ≠ 0) :
    IsSquare n ↔ ∀ p, Even (n.factorization p) := by
  constructor
  · rintro ⟨d, rfl⟩ p
    rw [← pow_two, Nat.factorization_pow]
    change Even (2 * d.factorization p)
    exact even_two_mul _
  · intro h
    refine ⟨∏ p ∈ n.primeFactors, p ^ (n.factorization p / 2), ?_⟩
    rw [← pow_two, ← Finset.prod_pow]
    calc
      n = ∏ p ∈ n.primeFactors, p ^ n.factorization p :=
        Nat.prod_primeFactors_pow_factorization hn
      _ = _ := by
        refine Finset.prod_congr rfl (fun p _ => ?_)
        rw [← pow_mul]
        congr 1
        have := Nat.even_iff.mp (h p)
        omega

private theorem isSquare_mul_iff_of_coprime {m n : ℕ} (hm : m ≠ 0) (hn : n ≠ 0)
    (hc : m.Coprime n) : IsSquare (m * n) ↔ IsSquare m ∧ IsSquare n := by
  rw [isSquare_iff_even_factorization (mul_ne_zero hm hn), isSquare_iff_even_factorization hm,
    isSquare_iff_even_factorization hn]
  have hzero (p : ℕ) : m.factorization p = 0 ∨ n.factorization p = 0 := by
    by_contra! h
    have hm' : p ∈ m.primeFactors := by
      rw [← Nat.support_factorization, Finsupp.mem_support_iff]
      exact h.1
    have hn' : p ∈ n.primeFactors := by
      rw [← Nat.support_factorization, Finsupp.mem_support_iff]
      exact h.2
    exact Finset.disjoint_left.mp hc.disjoint_primeFactors hm' hn'
  constructor
  · intro h
    have hadd (p : ℕ) := h p
    simp only [Nat.factorization_mul_apply_of_coprime hc] at hadd
    constructor <;> intro p <;> rcases hzero p with hzero | hzero
    · simp [hzero]
    · simpa only [hzero, add_zero] using hadd p
    · simpa only [hzero, zero_add] using hadd p
    · simp [hzero]
  · rintro ⟨hm', hn'⟩ p
    rw [Nat.factorization_mul_apply_of_coprime hc]
    exact (hm' p).add (hn' p)

noncomputable def squareIndicator : ArithmeticFunction ℕ := by
  classical
  exact ⟨fun n => if n ≠ 0 ∧ IsSquare n then 1 else 0, by simp⟩

theorem squareIndicator_apply (n : ℕ) [Decidable (n ≠ 0 ∧ IsSquare n)] :
    squareIndicator n = if n ≠ 0 ∧ IsSquare n then 1 else 0 := by
  classical
  by_cases h : n ≠ 0 ∧ IsSquare n <;> simp [squareIndicator, h]

theorem squareIndicator_multiplicative : squareIndicator.IsMultiplicative := by
  classical
  refine ⟨by rw [squareIndicator_apply]; simp, fun {m n} hc => ?_⟩
  by_cases hm : m = 0
  · simp [hm]
  by_cases hn : n = 0
  · simp [hn]
  simp only [squareIndicator_apply, isSquare_mul_iff_of_coprime hm hn hc]
  by_cases hms : IsSquare m <;> by_cases hns : IsSquare n <;> simp [hm, hn, hms, hns]

/-- The divisor envelope obtained by scaling primitive representations
through the possible square divisors of the norm. -/
noncomputable def normDivisorEnvelope (root : ArithmeticFunction ℕ) : ArithmeticFunction ℕ :=
  squareIndicator * root

theorem norm_divisor_envelope_prime_power {root : ArithmeticFunction ℕ} {p : ℕ} (hp : p.Prime) (k : ℕ) :
    normDivisorEnvelope root (p ^ k) =
      ∑ j ∈ Finset.range (k + 1), if (k - j) % 2 = 0 then root (p ^ j) else 0 := by
  classical
  rw [normDivisorEnvelope, ArithmeticFunction.mul_apply,
    Nat.sum_divisorsAntidiagonal' (fun a b => squareIndicator a * root b),
    Nat.sum_divisors_prime_pow hp]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  rw [Nat.pow_div (by simpa using Finset.mem_range.mp hj) hp.pos, squareIndicator_apply]
  have hSq (a : ℕ) : IsSquare (p ^ a) ↔ a % 2 = 0 := by
    rw [← Nat.even_iff, isSquare_iff_even_factorization (pow_ne_zero _ hp.ne_zero)]
    constructor
    · intro h
      simpa only [hp.factorization_pow, Finsupp.single_eq_same] using h p
    · intro ha q
      rw [hp.factorization_pow]
      by_cases hq : q = p
      · simpa [hq] using ha
      · simp [Finsupp.single_eq_of_ne hq]
  by_cases hparity : (k - j) % 2 = 0 <;> simp [hSq, hp.ne_zero, hparity]

theorem norm_divisor_envelope_prime_power_step {root : ArithmeticFunction ℕ} {p : ℕ} (hp : p.Prime) (k : ℕ) :
    normDivisorEnvelope root (p ^ (k + 2)) =
      normDivisorEnvelope root (p ^ k) + root (p ^ (k + 2)) := by
  rw [norm_divisor_envelope_prime_power hp (k + 2), norm_divisor_envelope_prime_power hp k,
    Finset.sum_range_succ, Finset.sum_range_succ]
  have he : (k + 2 - (k + 1)) % 2 ≠ 0 := by omega
  simp only [he, if_false, Nat.sub_self, Nat.zero_mod, if_true, add_zero]
  congr 1
  refine Finset.sum_congr rfl (fun j hj => ?_)
  have hj' := Finset.mem_range.mp hj
  rw [show (k + 2 - j) % 2 = (k - j) % 2 by omega]

private theorem binary_form_zero {b c x y : ℤ} (hdisc : 0 < 4 * c - b ^ 2)
    (h : form b c x y = 0) : x = 0 ∧ y = 0 := by
  have hf : (2 * x + b * y) ^ 2 + (4 * c - b ^ 2) * y ^ 2 = 0 := by
    unfold form at h
    linear_combination 4 * h
  have hy : y = 0 := by
    by_contra hne
    linarith only [hf, sq_nonneg (2 * x + b * y), mul_pos hdisc (sq_pos_of_ne_zero hne)]
  rw [hy] at hf
  have hx : x = 0 := by
    by_contra hne
    nlinarith only [hf, sq_pos_of_ne_zero hne]
  exact ⟨hx, hy⟩

/-- For discriminant magnitude greater than four, a root at a
positive norm identifies at most one opposite pair of vectors. -/
theorem binary_form_eq_or_neg_of_root_congruence {b c : ℤ} (hdisc : 4 < 4 * c - b ^ 2)
    {n : ℕ} (hn : 0 < n) {v w : ℤ × ℤ} {t : ZMod n}
    (hv : form b c v.1 v.2 = n) (hw : form b c w.1 w.2 = n)
    (hvt : (v.1 : ZMod n) + (v.2 : ZMod n) * t = 0)
    (hwt : (w.1 : ZMod n) + (w.2 : ZMod n) * t = 0) :
    w = v ∨ w = (-v.1, -v.2) := by
  let d := v.1 * w.2 - w.1 * v.2
  let g := 2 * v.1 * w.1 + b * (v.1 * w.2 + v.2 * w.1) + 2 * c * v.2 * w.2
  have hgram : g ^ 2 + (4 * c - b ^ 2) * d ^ 2 = 4 * (n : ℤ) ^ 2 := by
    calc
      _ = 4 * form b c v.1 v.2 * form b c w.1 w.2 := by dsimp [g, d, form]; ring
      _ = _ := by rw [hv, hw]; ring
  have hdvd : (n : ℤ) ∣ d := by
    apply (ZMod.intCast_zmod_eq_zero_iff_dvd ..).mp
    calc
      (d : ZMod n) = ((v.1 : ZMod n) + (v.2 : ZMod n) * t) * (w.2 : ZMod n) -
          ((w.1 : ZMod n) + (w.2 : ZMod n) * t) * (v.2 : ZMod n) := by
        dsimp [d]
        push_cast
        ring
      _ = _ := by rw [hvt, hwt]; ring
  have hnZ : (0 : ℤ) < n := by exact_mod_cast hn
  have hd : d = 0 := by
    by_contra hd0
    obtain ⟨k, hk⟩ := hdvd
    have hk0 : k ≠ 0 := by rintro rfl; simp at hk; exact hd0 hk
    have hk1 : 1 ≤ k ^ 2 := by have := sq_pos_of_ne_zero hk0; omega
    have hdSq : (n : ℤ) ^ 2 ≤ d ^ 2 := by
      rw [hk, mul_pow]
      exact le_mul_of_one_le_right (sq_nonneg _) hk1
    have hcoeff := mul_le_mul_of_nonneg_right (show 5 ≤ 4 * c - b ^ 2 by omega) (sq_nonneg d)
    nlinarith only [hgram, hdSq, hcoeff, sq_nonneg g, sq_pos_of_pos hnZ]
  have hg : g = 2 * (n : ℤ) ∨ g = -(2 * (n : ℤ)) := by
    apply sq_eq_sq_iff_eq_or_eq_neg.mp
    rw [hd] at hgram
    nlinarith only [hgram]
  rcases hg with hg | hg
  · have hform : form b c (w.1 - v.1) (w.2 - v.2) = 0 := by
      dsimp [g, form] at hv hw hg ⊢
      linear_combination hv + hw - hg
    have h := binary_form_zero (by omega : 0 < 4 * c - b ^ 2) hform
    exact Or.inl (Prod.ext (sub_eq_zero.mp h.1) (sub_eq_zero.mp h.2))
  · have hform : form b c (w.1 + v.1) (w.2 + v.2) = 0 := by
      dsimp [g, form] at hv hw hg ⊢
      linear_combination hv + hw + hg
    have h := binary_form_zero (by omega : 0 < 4 * c - b ^ 2) hform
    exact Or.inr (Prod.ext (eq_neg_of_add_eq_zero_left h.1) (eq_neg_of_add_eq_zero_left h.2))

/-- Each root accounts for at most two primitive signed representations. -/
theorem binary_primitive_norm_card_le {b c : ℤ} (hdisc : 4 < 4 * c - b ^ 2)
    {n : ℕ} (hn : 0 < n) (s : Finset (ℤ × ℤ))
    (hs : ∀ v ∈ s, Int.gcd v.1 v.2 = 1 ∧ form b c v.1 v.2 = n) :
    s.card ≤ 2 * normRootCount b c n := by
  classical
  letI : NeZero n := ⟨hn.ne'⟩
  let residue (v : ℤ × ℤ) := primitiveNormResidue n v.1 v.2
  have hspec (v : ℤ × ℤ) (hv : v ∈ s) :
      residue v ∈ normRootResidues b c n ∧ (v.1 : ZMod n) + (v.2 : ZMod n) * residue v = 0 := by
    have h := primitiveNormResidue_spec (hs v hv).1 (hs v hv).2
    exact ⟨mem_normRootResidues.mpr h.1, h.2⟩
  rw [norm_root_count_apply]
  apply Finset.card_le_mul_card_image_of_maps_to (f := residue) (s := s) (t := normRootResidues b c n)
    (fun v hv => (hspec v hv).1) 2
  intro t _
  let fiber := s.filter (fun v => residue v = t)
  by_cases hne : fiber.Nonempty
  · obtain ⟨v, hv⟩ := hne
    obtain ⟨hrep, ht⟩ := Finset.mem_filter.mp hv
    have hvt := (hspec v hrep).2
    rw [ht] at hvt
    have hsub : fiber ⊆ {v, (-v.1, -v.2)} := by
      intro w hw
      obtain ⟨hwrep, hwt⟩ := Finset.mem_filter.mp hw
      have hwcong := (hspec w hwrep).2
      rw [hwt] at hwcong
      simpa only [Finset.mem_insert, Finset.mem_singleton] using
        binary_form_eq_or_neg_of_root_congruence hdisc hn (hs v hrep).2 (hs w hwrep).2 hvt hwcong
    exact (Finset.card_le_card hsub).trans (by simpa using Finset.card_insert_le v {(-v.1, -v.2)})
  · change fiber.card ≤ 2
    rw [Finset.not_nonempty_iff_eq_empty.mp hne]
    simp

/-- The full signed count is bounded by the convolution of the root count
with square divisibility. The supplied finite fibers enumerate actual coordinates. -/
theorem binary_norm_fiber_card_le {b c : ℤ} (hdisc : 4 < 4 * c - b ^ 2)
    (fiber : ℕ → Finset (ℤ × ℤ))
    (hfiber : ∀ n v, v ∈ fiber n ↔ form b c v.1 v.2 = n) {n : ℕ} (hn : 0 < n) :
    (fiber n).card ≤ 2 * normDivisorEnvelope (normRootCount b c) n := by
  classical
  let factors := n.divisorsAntidiagonal.filter (fun a => IsSquare a.1)
  let primitive (a : ℕ × ℕ) := (fiber a.2).filter fun v => Int.gcd v.1 v.2 = 1
  let scale (d : ℕ) (v : ℤ × ℤ) := ((d : ℤ) * v.1, (d : ℤ) * v.2)
  let reconstruct (t : Σ _ : ℕ × ℕ, ℤ × ℤ) := scale t.1.1.sqrt t.2
  have hsubset : fiber n ⊆ (factors.sigma primitive).image reconstruct := by
    intro v hv
    let d := Int.gcd v.1 v.2
    have hN := (hfiber n v).mp hv
    have hd : 0 < d := by
      by_contra! hd
      obtain ⟨h1, h2⟩ := Int.gcd_eq_zero_iff.mp (show d = 0 by omega)
      simp [form, h1, h2] at hN
      omega
    obtain ⟨x, hx⟩ := Int.gcd_dvd_left v.1 v.2
    obtain ⟨y, hy⟩ := Int.gcd_dvd_right v.1 v.2
    let w : ℤ × ℤ := (x, y)
    have hscale : scale d w = v := Prod.ext hx.symm hy.symm
    have hgcd : Int.gcd w.1 w.2 = 1 := by
      apply Nat.eq_of_mul_eq_mul_left hd
      calc
        _ = Int.gcd ((d : ℤ) * w.1) ((d : ℤ) * w.2) := by rw [Int.gcd_mul_left, Int.natAbs_natCast]
        _ = d := by change Int.gcd (scale d w).1 (scale d w).2 = d; rw [hscale]
        _ = _ := (mul_one _).symm
    have hprod : (d : ℤ) ^ 2 * form b c w.1 w.2 = n := by
      calc
        _ = form b c (scale d w).1 (scale d w).2 := by dsimp [scale, form]; ring
        _ = _ := by rw [hscale, hN]
    have hw₀ : 0 ≤ form b c w.1 w.2 := by
      by_contra! hneg
      have ht := mul_neg_of_pos_of_neg (sq_pos_of_pos (show 0 < (d : ℤ) by exact_mod_cast hd)) hneg
      omega
    let a := (form b c w.1 w.2).toNat
    have ha : (a : ℤ) = form b c w.1 w.2 := Int.toNat_of_nonneg hw₀
    rw [← ha] at hprod
    have hpair : (d ^ 2, a) ∈ factors :=
      Finset.mem_filter.mpr
        ⟨Nat.mem_divisorsAntidiagonal.mpr ⟨by exact_mod_cast hprod, hn.ne'⟩, ⟨d, pow_two d⟩⟩
    refine Finset.mem_image.mpr ⟨⟨(d ^ 2, a), w⟩, Finset.mem_sigma.mpr
      ⟨hpair, Finset.mem_filter.mpr ⟨(hfiber a w).mpr ha.symm, hgcd⟩⟩, ?_⟩
    dsimp only [reconstruct]
    rw [pow_two, Nat.sqrt_eq, hscale]
  calc
    _ ≤ ((factors.sigma primitive).image reconstruct).card := Finset.card_le_card hsubset
    _ ≤ (factors.sigma primitive).card := Finset.card_image_le
    _ = ∑ a ∈ factors, (primitive a).card := Finset.card_sigma factors primitive
    _ ≤ ∑ a ∈ factors, 2 * normRootCount b c a.2 := by
      refine Finset.sum_le_sum fun a ha => binary_primitive_norm_card_le hdisc
        (Nat.pos_of_ne_zero (Nat.right_ne_zero_of_mem_divisorsAntidiagonal (Finset.mem_filter.mp ha).1))
        (primitive a) ?_
      intro v hv
      obtain ⟨h, hg⟩ := Finset.mem_filter.mp hv
      exact ⟨hg, (hfiber a.2 v).mp h⟩
    _ = 2 * normDivisorEnvelope (normRootCount b c) n := by
      rw [← Finset.mul_sum, normDivisorEnvelope, ArithmeticFunction.mul_apply]
      congr 1
      change ∑ a ∈ n.divisorsAntidiagonal.filter _, _ = _
      rw [Finset.sum_filter]
      refine Finset.sum_congr rfl fun a ha => ?_
      rw [squareIndicator_apply]
      have ha₀ := Nat.left_ne_zero_of_mem_divisorsAntidiagonal ha
      by_cases hs : IsSquare a.1 <;> simp [hs, ha₀]

end QuadraticNorms
