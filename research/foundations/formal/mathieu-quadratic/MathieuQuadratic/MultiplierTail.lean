import MathieuQuadratic.SharpMultiplier
import Mathlib.FieldTheory.RatFunc.AsPolynomial

/-!
# The rational test for a fixed block expansion

The rational function is presented by a common denominator. At a
nonnegative integer the denominator is a product of positive rational
rising factorials, so the test reduces to an ordinary polynomial.
-/

open Polynomial
open scoped BigOperators

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

private noncomputable def tailDenominatorTerm (M γ j : ℕ) : ℚ[X] :=
  (ascPochhammer ℚ j).comp (C ((M : ℚ) / 2 + 1) + C (γ : ℚ) * X)

private theorem tailDenominatorTerm_eval (M γ j m : ℕ) :
    (tailDenominatorTerm M γ j).eval (m : ℚ) =
      (ascPochhammer ℚ j).eval ((M : ℚ) / 2 + (γ * m : ℕ) + 1) := by
  rw [tailDenominatorTerm, eval_comp]
  congr 1
  simp
  ring

private noncomputable def tailNumeratorTerm (M : ℕ) (ℓ : ℤ) (j : ℕ) : ℚ[X] :=
  if (M + 2 * j : ℕ) ≤ ℓ then
    let n := (ℓ - (M + 2 * j : ℕ)).toNat
    C ((ascPochhammer ℚ j).eval ((M : ℚ) / 2) / n.factorial) *
      (descPochhammer ℚ n).comp (C (M : ℚ) * X - C (↑(M + 2 * j) : ℚ))
  else 0

private theorem tailNumeratorTerm_eval (M : ℕ) (ℓ : ℤ) (j m : ℕ) (hm : M + 2 * j ≤ M * m) :
    (tailNumeratorTerm M ℓ j).eval (m : ℚ) =
      if (M + 2 * j : ℕ) ≤ ℓ then
        (ascPochhammer ℚ j).eval ((M : ℚ) / 2) *
          ((Ring.choose ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ)
            (ℓ - (M + 2 * j : ℕ)).toNat : ℤ) : ℚ)
      else 0 := by
  unfold tailNumeratorTerm
  split_ifs
  · let n := (ℓ - (M + 2 * j : ℕ)).toNat
    let u := M * m - (M + 2 * j)
    have hu : (u : ℤ) = ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ) := Int.ofNat_sub hm
    rw [← hu, Ring.choose_natCast, Int.cast_natCast, Nat.cast_choose_eq_descPochhammer_div,
      eval_mul, eval_C, eval_comp, eval_sub, eval_mul, eval_C, eval_X, eval_C,
      show (M : ℚ) * m - (↑(M + 2 * j) : ℚ) = ↑u by
        dsimp [u]; rw [Nat.cast_sub hm, Nat.cast_mul]]
    ring
  · simp

/-- The polynomial common denominator in the draft's tail test. -/
noncomputable def multiplierTailDenominator (M γ : ℕ) (support : Finset (ℤ × ℕ)) : K[X] :=
  ∏ p ∈ support, (tailDenominatorTerm M γ p.2).map (algebraMap ℚ K)

/-- The numerator after putting the terms of R_h over a common denominator. -/
noncomputable def multiplierTailNumerator (M γ : ℕ) (support : Finset (ℤ × ℕ))
    (c : ℤ × ℕ → K) : K[X] :=
  ∑ p ∈ support, C (c p) * (tailNumeratorTerm M p.1 p.2).map (algebraMap ℚ K) *
    ∏ t ∈ support.erase p, (tailDenominatorTerm M γ t.2).map (algebraMap ℚ K)

/-- The rational test in the draft, written over the common denominator. -/
noncomputable def multiplierTailFunction (M γ : ℕ) (support : Finset (ℤ × ℕ))
    (c : ℤ × ℕ → K) : RatFunc K :=
  algebraMap K[X] (RatFunc K) (multiplierTailNumerator M γ support c) /
    algebraMap K[X] (RatFunc K) (multiplierTailDenominator M γ support)

private theorem eval_map_at_nat (p : ℚ[X]) (m : ℕ) :
    (p.map (algebraMap ℚ K)).eval (m : K) = algebraMap ℚ K (p.eval (m : ℚ)) := by
  rw [eval_map, ← map_natCast (algebraMap ℚ K) m]
  simp

private theorem tailDenominatorTerm_eval_ne (M γ j m : ℕ) :
    ((tailDenominatorTerm M γ j).map (algebraMap ℚ K)).eval (m : K) ≠ 0 := by
  rw [eval_map_at_nat, tailDenominatorTerm_eval]
  exact (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
    (ascPochhammer_pos _ _ (by positivity)).ne'

theorem multiplierTailDenominator_eval_ne (M γ : ℕ) (support : Finset (ℤ × ℕ)) (m : ℕ) :
    (multiplierTailDenominator (K := K) M γ support).eval (m : K) ≠ 0 := by
  rw [multiplierTailDenominator, eval_prod]
  exact Finset.prod_ne_zero_iff.2 (fun p _ ↦ tailDenominatorTerm_eval_ne M γ p.2 m)

theorem multiplierTailFunction_eq_zero_iff (M γ : ℕ) (support : Finset (ℤ × ℕ))
    (c : ℤ × ℕ → K) :
    multiplierTailFunction M γ support c = 0 ↔ multiplierTailNumerator M γ support c = 0 := by
  have hd : multiplierTailDenominator (K := K) M γ support ≠ 0 := by
    intro h
    exact multiplierTailDenominator_eval_ne (K := K) M γ support 0 (by rw [h, eval_zero])
  have hi := IsFractionRing.injective K[X] (RatFunc K)
  rw [multiplierTailFunction, div_eq_zero_iff,
    map_eq_zero_iff _ hi, map_eq_zero_iff _ hi, or_iff_left hd]

/-- The common-denominator expression equals the termwise rational sum. -/
theorem multiplierTailFunction_eq_sum (M γ : ℕ) (support : Finset (ℤ × ℕ))
    (c : ℤ × ℕ → K) :
    multiplierTailFunction M γ support c =
      ∑ p ∈ support,
        algebraMap K[X] (RatFunc K) (C (c p) * (tailNumeratorTerm M p.1 p.2).map (algebraMap ℚ K)) /
          algebraMap K[X] (RatFunc K) ((tailDenominatorTerm M γ p.2).map (algebraMap ℚ K)) := by
  rw [multiplierTailFunction, multiplierTailNumerator, map_sum, Finset.sum_div]
  refine Finset.sum_congr rfl fun p hp ↦ ?_
  rw [multiplierTailDenominator, ← Finset.prod_erase_mul _ _ hp, map_mul, map_mul,
    map_mul, map_prod]
  have hd (t : ℤ × ℕ) :
      algebraMap K[X] (RatFunc K) ((tailDenominatorTerm M γ t.2).map (algebraMap ℚ K)) ≠ 0 := by
    intro h
    have hpoly : (tailDenominatorTerm M γ t.2).map (algebraMap ℚ K) = 0 :=
      (IsFractionRing.injective K[X] (RatFunc K)) (h.trans (map_zero _).symm)
    exact tailDenominatorTerm_eval_ne (K := K) M γ t.2 0 (by rw [hpoly, eval_zero])
  rw [mul_comm (∏ t ∈ support.erase p, _), mul_div_mul_right _ _
    (Finset.prod_ne_zero_iff.2 (fun t _ ↦ hd t))]

private theorem betaMoment_cast_ne_zero (M k j : ℕ) (hM : 0 < M) :
    algebraMap ℚ K (mixedBetaMoment M k j) ≠ 0 :=
  (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
    (mixedBetaMoment_pos M k j hM).ne'

/-- The finite Laurent expansion version of the tail formula, with its
constant beta factor and the denominator displayed separately. -/
theorem product_multiplier_tail_formula (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M γ : ℕ) (hM : 0 < M) (hγ : M < γ)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K))
    (support : Finset (ℤ × ℕ)) (c : ℤ × ℕ → K) (m : ℕ)
    (hm : ∀ p ∈ support, p.1.natAbs + M + 2 * p.2 ≤ m) :
    trace ((∑ p ∈ support, c p • ((↑(U ^ p.1) : A) * S ^ p.2)) *
        productPhase U S M M γ ^ m) =
      algebraMap ℚ K (mixedBetaMoment M (γ * m) 0) *
        ((multiplierTailNumerator M γ support c).eval (m : K) /
          (multiplierTailDenominator M γ support).eval (m : K)) := by
  have hd := multiplierTailDenominator_eval_ne (K := K) M γ support m
  simp only [Finset.sum_mul, map_sum, smul_mul_assoc, map_smul, smul_eq_mul,
    multiplierTailNumerator, multiplierTailDenominator, eval_finsetSum, eval_mul, eval_C, eval_prod, Finset.sum_div,
    Finset.mul_sum]
  refine Finset.sum_congr rfl fun p hp ↦ ?_
  rw [product_tail_monomial trace U S hM weight height p.1 p.2 m hγ (hm p hp), map_mul,
    eval_map_at_nat, tailNumeratorTerm_eval M p.1 p.2 m (by
      have hh := hm p hp
      have hprod := Nat.mul_le_mul_right m (show 1 ≤ M by omega)
      rw [one_mul] at hprod
      omega)]
  rw [← Finset.prod_erase_mul support
    (fun t ↦ ((tailDenominatorTerm M γ t.2).map (algebraMap ℚ K)).eval (m : K)) hp]
  let d := ((tailDenominatorTerm M γ p.2).map (algebraMap ℚ K)).eval (m : K)
  let e := ∏ t ∈ support.erase p, ((tailDenominatorTerm M γ t.2).map (algebraMap ℚ K)).eval (m : K)
  have hd' : d ≠ 0 := tailDenominatorTerm_eval_ne M γ p.2 m
  have he : e ≠ 0 := Finset.prod_ne_zero_iff.2 (fun t _ ↦ tailDenominatorTerm_eval_ne M γ t.2 m)
  have hde : d = algebraMap ℚ K ((ascPochhammer ℚ p.2).eval
      ((M : ℚ) / 2 + (γ * m : ℕ) + 1)) := by
    dsimp [d]
    rw [eval_map_at_nat, tailDenominatorTerm_eval]
  split_ifs
  · simp only [map_mul, map_div₀]
    rw [← hde]
    let b := algebraMap ℚ K (mixedBetaMoment M (γ * m) 0)
    let n := algebraMap ℚ K ((ascPochhammer ℚ p.2).eval ((M : ℚ) / 2))
    let v := algebraMap ℚ K ((Ring.choose ((M * m : ℕ) - (M + 2 * p.2 : ℕ) : ℤ)
      (p.1 - (M + 2 * p.2 : ℕ)).toNat : ℤ) : ℚ)
    change c p * (b * (n / d * v)) = b * (c p * (n * v) * e / (e * d))
    field_simp [hd', he]
  · simp

/-- A polynomial over a characteristic-zero field vanishes at all large
integers exactly when it is zero. -/
private theorem eventually_eval_zero_iff (p : K[X]) :
    (∃ m₀ : ℕ, ∀ m ≥ m₀, p.eval (m : K) = 0) ↔ p = 0 := by
  constructor
  · rintro ⟨m₀, hm⟩
    apply p.eq_of_infinite_eval_eq 0
    apply (Set.infinite_range_of_injective (f := fun m : ℕ ↦ ((m + m₀ : ℕ) : K)) (by
      intro a b hab
      exact Nat.add_right_cancel (Nat.cast_injective hab))).mono
    rintro _ ⟨m, rfl⟩
    simpa using hm (m + m₀) (by omega)
  · rintro rfl
    exact ⟨0, by simp⟩

/-- Eventual vanishing of the marked moments is equivalent to the
rational function test, for a given finite block expansion. -/
theorem product_multiplier_tail_iff (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M γ : ℕ) (hM : 0 < M) (hγ : M < γ)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K))
    (support : Finset (ℤ × ℕ)) (c : ℤ × ℕ → K) :
    (∃ m₀ : ℕ, ∀ m ≥ m₀,
        trace ((∑ p ∈ support, c p • ((↑(U ^ p.1) : A) * S ^ p.2)) *
          productPhase U S M M γ ^ m) = 0) ↔
      multiplierTailFunction M γ support c = 0 := by
  let cutoff := support.sup (fun p ↦ p.1.natAbs + M + 2 * p.2)
  have hz (m : ℕ) (hm : cutoff ≤ m) :
      trace ((∑ p ∈ support, c p • ((↑(U ^ p.1) : A) * S ^ p.2)) *
          productPhase U S M M γ ^ m) = 0 ↔
        (multiplierTailNumerator M γ support c).eval (m : K) = 0 := by
    rw [product_multiplier_tail_formula trace U S M γ hM hγ weight height support c m
      (fun p hp ↦ (Finset.le_sup hp).trans hm), mul_eq_zero, div_eq_zero_iff,
      or_iff_right (betaMoment_cast_ne_zero M (γ * m) 0 hM),
      or_iff_left (multiplierTailDenominator_eval_ne M γ support m)]
  rw [multiplierTailFunction_eq_zero_iff, ← eventually_eval_zero_iff]
  constructor <;> rintro ⟨m₀, hm⟩ <;> refine ⟨max cutoff m₀, fun m h ↦ ?_⟩
  · exact (hz m ((le_max_left _ _).trans h)).1 (hm m ((le_max_right _ _).trans h))
  · exact (hz m ((le_max_left _ _).trans h)).2 (hm m ((le_max_right _ _).trans h))

/-- A nonzero rational test gives nonzero marked moments at every
sufficiently large order. -/
theorem product_multiplier_tail_eventually_ne_zero (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M γ : ℕ) (hM : 0 < M) (hγ : M < γ)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K))
    (support : Finset (ℤ × ℕ)) (c : ℤ × ℕ → K)
    (hn : multiplierTailFunction M γ support c ≠ 0) :
    ∃ m₀ : ℕ, ∀ m ≥ m₀,
      trace ((∑ p ∈ support, c p • ((↑(U ^ p.1) : A) * S ^ p.2)) *
        productPhase U S M M γ ^ m) ≠ 0 := by
  let p := multiplierTailNumerator M γ support c
  have hp : p ≠ 0 := mt (multiplierTailFunction_eq_zero_iff M γ support c).2 hn
  have hroots := (Polynomial.finite_setOf_isRoot hp).preimage
    (Nat.cast_injective (R := K)).injOn
  obtain ⟨r, hr⟩ := hroots.bddAbove
  let cutoff := support.sup (fun q ↦ q.1.natAbs + M + 2 * q.2)
  refine ⟨max cutoff (r + 1), fun m hm ↦ ?_⟩
  have hcut (q : ℤ × ℕ) (hq : q ∈ support) : q.1.natAbs + M + 2 * q.2 ≤ m :=
    (Finset.le_sup hq).trans ((le_max_left _ _).trans hm)
  have hmne : p.eval (m : K) ≠ 0 := by
    intro hz
    have hle : m ≤ r := hr (show (m : K) ∈ {x | p.IsRoot x} from hz)
    have hgt := (le_max_right cutoff (r + 1)).trans hm
    omega
  rw [product_multiplier_tail_formula trace U S M γ hM hγ weight height support c m hcut]
  exact mul_ne_zero (betaMoment_cast_ne_zero M (γ * m) 0 hM)
    (div_ne_zero hmne (multiplierTailDenominator_eval_ne M γ support m))

/-- The ternary cancellation example has a nonzero first moment and
vanishes from the third order. -/
theorem draft_ternary_tail_example (a : Fin 2 → K) (b : Fin 1 → K)
    {trace : NormComplement a b →ₗ[K] K} (h : GeometricTrace a b trace) :
    let U := normUnit a b
    let S := height a b
    let F := complementProductPhase a b 1 1 2
    let H := (4 : K) • ((U : NormComplement a b) ^ 4 * S) +
      (15 : K) • ((U : NormComplement a b) ^ 3 * S) - (U : NormComplement a b)
    trace (H * F) ≠ 0 ∧ ∀ m ≥ 3, trace (H * F ^ m) = 0 := by
  intro U S F H
  have ht := h.complementTrace.normalizer a b Nat.one_pos
  have hmark (m : ℕ) (hm : 0 < m) :
      trace ((U : NormComplement a b) * F ^ m) = algebraMap ℚ K (mixedBetaMoment 1 (2 * m) 0) := by
    simpa only [pow_one] using (draft_product_witness a b h Nat.one_pos 2 (by omega)).2 m hm
  have hterm (m ℓ : ℕ) :
      trace ((U : NormComplement a b) ^ ℓ * S * F ^ m) =
        if ℓ ≤ m then algebraMap ℚ K (mixedBetaMoment 1 (2 * m) 1 *
          ((Ring.choose ((m : ℤ) - 3) (m - ℓ) : ℤ) : ℚ)) else 0 := by
    have he := product_mixedMoment trace U S 1 1 1 2 m 1 Nat.one_pos ℓ (by omega) ht.weight ht.height
    by_cases hℓ : ℓ ≤ m
    · rw [if_pos hℓ]
      have hℓ' : (ℓ : ℤ) ≤ (1 * m : ℕ) := by omega
      rw [if_pos hℓ', one_mul, ← Int.ofNat_sub hℓ, Int.toNat_natCast] at he
      simpa [F, complementProductPhase] using he
    · rw [if_neg hℓ]
      have hℓ' : ¬ (ℓ : ℤ) ≤ (1 * m : ℕ) := by omega
      rw [if_neg hℓ'] at he
      simpa [F, complementProductPhase] using he
  have hsum (m : ℕ) :
      trace (H * F ^ m) = 4 * trace ((U : NormComplement a b) ^ 4 * S * F ^ m) +
        15 * trace ((U : NormComplement a b) ^ 3 * S * F ^ m) -
          trace ((U : NormComplement a b) * F ^ m) := by
    simp [H, sub_mul, add_mul, map_smul, smul_eq_mul]
  constructor
  · have hfirst : trace (H * F) = -(algebraMap ℚ K (mixedBetaMoment 1 2 0)) := by
      rw [← pow_one F, hsum, hterm 1 4, hterm 1 3, hmark 1 (by decide)]
      norm_num
    rw [hfirst]
    exact neg_ne_zero.mpr (betaMoment_cast_ne_zero 1 2 0 Nat.one_pos)
  · intro m hm
    let C₀ := mixedBetaMoment 1 (2 * m) 0
    let C₁ := mixedBetaMoment 1 (2 * m) 1
    have hn : ((m : ℤ) - 3) = (↑(m - 3) : ℤ) := (Int.ofNat_sub hm).symm
    have h3 : trace ((U : NormComplement a b) ^ 3 * S * F ^ m) = algebraMap ℚ K C₁ := by
      rw [hterm, if_pos hm, hn, Ring.choose_natCast, Nat.choose_self]
      simp [C₁]
    have h4 : trace ((U : NormComplement a b) ^ 4 * S * F ^ m) =
        algebraMap ℚ K (C₁ * (↑(m - 3) : ℚ)) := by
      rw [hterm]
      by_cases h4 : 4 ≤ m
      · rw [if_pos h4, hn, Ring.choose_natCast, Nat.choose_symm_of_eq_add
          (show m - 3 = m - 4 + 1 by omega), Nat.choose_one_right, Int.cast_natCast]
      · simp [h4, show m - 3 = 0 by omega]
    have hbeta := mixedBetaMoment_succ_height 1 (2 * m) 0 Nat.one_pos
    change (↑(1 + 2 * 0 : ℕ) : ℚ) * C₀ =
      (↑(1 + 2 * 0 + 2 * (2 * m) + 2 : ℕ) : ℚ) * C₁ at hbeta
    push_cast at hbeta
    ring_nf at hbeta
    rw [hsum, h4, h3, hmark m (by omega)]
    rw [← map_ofNat (algebraMap ℚ K) 4, ← map_ofNat (algebraMap ℚ K) 15, ← map_mul,
      ← map_mul, ← map_add, ← map_sub]
    apply (map_eq_zero_iff _ (RingHom.injective (algebraMap ℚ K))).2
    change 4 * (C₁ * (↑(m - 3) : ℚ)) + 15 * C₁ - C₀ = 0
    rw [Nat.cast_sub hm]
    linear_combination -hbeta

end MathieuQuadratic
