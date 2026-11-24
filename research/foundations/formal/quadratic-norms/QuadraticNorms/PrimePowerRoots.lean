import QuadraticNorms.NormRoots
import QuadraticNorms.PrimeCores
import Mathlib.Algebra.Field.ZMod

namespace QuadraticNorms

private def rootPolynomial (b c t : ℤ) : ℤ := t ^ 2 - b * t + c

theorem HasNormRoot.of_dvd {b c : ℤ} {m n : ℕ} (h : HasNormRoot b c n)
    (hm : m ∣ n) : HasNormRoot b c m := by
  obtain ⟨t, ht⟩ := h
  let f := ZMod.castHom hm (ZMod m)
  exact ⟨f t, by simpa only [map_add, map_sub, map_mul, map_pow, map_intCast, map_zero]
    using congrArg f ht⟩

theorem HasNormRoot.mul {b c : ℤ} {m n : ℕ} (hmn : m.Coprime n)
    (hm : HasNormRoot b c m) (hn : HasNormRoot b c n) : HasNormRoot b c (m * n) := by
  obtain ⟨x, hx⟩ := hm
  obtain ⟨y, hy⟩ := hn
  refine ⟨(ZMod.chineseRemainder hmn).symm (x, y), ?_⟩
  apply (ZMod.chineseRemainder hmn).injective
  simp only [map_add, map_sub, map_mul, map_pow, map_intCast,
    map_zero, RingEquiv.apply_symm_apply]
  exact Prod.ext hx hy

theorem hasNormRoot_iff_primePowers {b c : ℤ} {n : ℕ} (hn : 0 < n) :
    HasNormRoot b c n ↔
      ∀ p k : ℕ, p.Prime → 0 < k → p ^ k ∣ n → HasNormRoot b c (p ^ k) := by
  refine ⟨fun h _ _ _ _ hdiv => h.of_dvd hdiv, ?_⟩
  refine Nat.recOnPosPrimePosCoprime (motive := fun n => 0 < n →
    (∀ p k : ℕ, p.Prime → 0 < k → p ^ k ∣ n → HasNormRoot b c (p ^ k)) →
      HasNormRoot b c n) ?_ ?_ ?_ ?_ n hn
  · intro p k hp hk _ h
    exact h p k hp (by omega) dvd_rfl
  · omega
  · intro _ _
    exact ⟨0, Subsingleton.elim ..⟩
  · intro m n hm hn hcoprime ih₁ ih₂ _ h
    exact (ih₁ (by omega) (fun p k hp hk hd => h p k hp hk
      (hd.trans (dvd_mul_right ..)))).mul hcoprime
      (ih₂ (by omega) (fun p k hp hk hd => h p k hp hk
        (hd.trans (dvd_mul_left ..))))

private theorem rootPolynomial_dvd_iff {b c : ℤ} {n : ℕ} [NeZero n] :
    HasNormRoot b c n ↔ ∃ t : ℤ, (n : ℤ) ∣ rootPolynomial b c t := by
  refine ⟨?_, ?_⟩
  · rintro ⟨t, ht⟩
    refine ⟨t.val, (ZMod.intCast_zmod_eq_zero_iff_dvd ..).mp ?_⟩
    simpa [rootPolynomial] using ht
  · rintro ⟨t, ht⟩
    exact ⟨t, by simpa [rootPolynomial] using (ZMod.intCast_zmod_eq_zero_iff_dvd ..).mpr ht⟩

/-- An unramified root modulo a prime lifts to every positive power. -/
theorem HasNormRoot.of_unramified {b c : ℤ} {p : ℕ} (hp : p.Prime)
    (hdisc : ¬ (p : ℤ) ∣ b ^ 2 - 4 * c) (h : HasNormRoot b c p) (k : ℕ) :
    HasNormRoot b c (p ^ k) := by
  letI : Fact p.Prime := ⟨hp⟩
  obtain ⟨t, ht⟩ := rootPolynomial_dvd_iff.mp h
  have hderiv (a : ℤ) (ha : (p : ℤ) ∣ rootPolynomial b c a) :
      ((2 * a - b : ℤ) : ZMod p) ≠ 0 := by
    intro hd
    have ha' := (ZMod.intCast_zmod_eq_zero_iff_dvd ..).mpr ha
    apply hdisc
    apply (ZMod.intCast_zmod_eq_zero_iff_dvd ..).mp
    push_cast at *
    unfold rootPolynomial at ha
    have ha'' : (a : ZMod p) ^ 2 - (b : ZMod p) * (a : ZMod p) + (c : ZMod p) = 0 := by
      exact_mod_cast ha'
    linear_combination hd * (2 * (a : ZMod p) - (b : ZMod p)) - 4 * ha''
  cases k with
  | zero =>
    rw [pow_zero]
    exact ⟨0, Subsingleton.elim ..⟩
  | succ k =>
    have hlift : ∀ k : ℕ, ∃ a : ℤ, (p : ℤ) ^ (k + 1) ∣ rootPolynomial b c a := by
      intro k
      induction k with
      | zero => exact ⟨t, by simpa using ht⟩
      | succ k ih =>
        obtain ⟨a, j, hj⟩ := ih
        have ha : (p : ℤ) ∣ rootPolynomial b c a := by
          exact (dvd_pow_self _ (by omega : k + 1 ≠ 0)).trans ⟨j, hj⟩
        let r : ZMod p := -(j : ZMod p) / ((2 * a - b : ℤ) : ZMod p)
        have hr : (j : ZMod p) + (r.val : ZMod p) * ((2 * a - b : ℤ) : ZMod p) = 0 := by
          rw [ZMod.natCast_val, ZMod.cast_id]
          dsimp only [r]
          rw [div_mul_cancel₀ _ (hderiv a ha), add_neg_cancel]
        obtain ⟨w, hw⟩ := (ZMod.intCast_zmod_eq_zero_iff_dvd (j + r.val * (2 * a - b)) p).mp
          (by simpa only [Int.cast_add, Int.cast_mul, Int.cast_natCast] using hr)
        refine ⟨a + (p : ℤ) ^ (k + 1) * r.val, w + (p : ℤ) ^ k * (r.val : ℤ) ^ 2, ?_⟩
        generalize hd : (p : ℤ) ^ k = d
        simp only [pow_succ, hd] at hj ⊢
        unfold rootPolynomial at hj ⊢
        linear_combination hj + (d * p) * hw
    letI : NeZero (p ^ (k + 1)) := ⟨pow_ne_zero _ hp.ne_zero⟩
    apply rootPolynomial_dvd_iff.mpr
    obtain ⟨a, ha⟩ := hlift k
    exact ⟨a, by exact_mod_cast ha⟩

/-- A root modulo an odd ramified prime has no lift to its square
when the discriminant is minus that prime. -/
theorem not_hasNormRoot_ramified_sq {b c : ℤ} {p : ℕ} (hp : p.Prime)
    (hdisc : 4 * c - b ^ 2 = (p : ℤ)) :
    ¬ HasNormRoot b c (p ^ 2) := by
  letI : NeZero (p ^ 2) := ⟨pow_ne_zero _ hp.ne_zero⟩
  intro h
  obtain ⟨a, ha⟩ := rootPolynomial_dvd_iff.mp h
  have ha' : (p : ℤ) ^ 2 ∣ rootPolynomial b c a := by exact_mod_cast ha
  have hsplit : (2 * a - b) ^ 2 = 4 * rootPolynomial b c a - p := by
    unfold rootPolynomial
    linear_combination -hdisc
  have hpz : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have hlin : (p : ℤ) ∣ 2 * a - b := by
    apply hpz.dvd_of_dvd_pow (n := 2)
    rw [hsplit]
    exact Int.dvd_sub
      (((dvd_pow_self _ (by decide : 2 ≠ 0)).trans ha').trans (dvd_mul_left ..)) dvd_rfl
  have hdv : (p : ℤ) ^ 2 ∣ (p : ℤ) := by
    have heq : 4 * rootPolynomial b c a - (2 * a - b) ^ 2 = (p : ℤ) := by
      rw [hsplit]
      ring
    have hsub : (p : ℤ) ^ 2 ∣ 4 * rootPolynomial b c a - (2 * a - b) ^ 2 :=
      Int.dvd_sub (ha'.trans (dvd_mul_left ..)) (pow_dvd_pow_of_dvd hlin 2)
    simpa only [heq] using hsub
  have hle : p ^ 2 ≤ p := Nat.le_of_dvd hp.pos (by exact_mod_cast hdv)
  nlinarith [hp.two_le]

end QuadraticNorms
