import QuadraticNorms.PrimePowerRoots
import QuadraticNorms.ExactSumsets
import Mathlib.NumberTheory.LegendreSymbol.QuadraticReciprocity

namespace QuadraticNorms

private theorem hasNormRoot_iff_square_discriminant {b c : ℤ} {p : ℕ}
    (hp : p.Prime) (hp₂ : p ≠ 2) :
    HasNormRoot b c p ↔ IsSquare ((b ^ 2 - 4 * c : ℤ) : ZMod p) := by
  letI : Fact p.Prime := ⟨hp⟩
  have h₂ : (2 : ZMod p) ≠ 0 := by
    simpa only [Nat.cast_ofNat] using (@ZMod.prime_ne_zero p 2 ⟨hp⟩ ⟨Nat.prime_two⟩ hp₂)
  constructor
  · rintro ⟨t, ht⟩
    refine ⟨2 * t - (b : ZMod p), ?_⟩
    push_cast
    linear_combination -4 * ht
  · rintro ⟨s, hs⟩
    refine ⟨(s + (b : ZMod p)) / 2, ?_⟩
    push_cast at hs
    field_simp
    linear_combination -hs

private theorem legendre_negative_prime {p q : ℕ} [Fact p.Prime] [Fact q.Prime]
    (hp₂ : p ≠ 2) (hq₄ : q % 4 = 3) :
    legendreSym p (-(q : ℤ)) = legendreSym q p := by
  have hq₂ : q ≠ 2 := by omega
  have hpodd := (Nat.Prime.mod_two_eq_one_iff_ne_two (Fact.out : p.Prime)).mpr hp₂
  rw [legendreSym.at_neg hp₂, legendreSym.quadratic_reciprocity' hp₂ hq₂,
    ZMod.χ₄_eq_neg_one_pow hpodd, mul_comm (p / 2) (q / 2), pow_mul,
    ZMod.neg_one_pow_div_two_of_three_mod_four hq₄]

private theorem square_negative_prime_iff {p q : ℕ}
    (hp : p.Prime) (hq : q.Prime) (hp₂ : p ≠ 2) (hq₄ : q % 4 = 3) (hpq : p ≠ q) :
    IsSquare ((-(q : ℤ)) : ZMod p) ↔ IsSquare (p : ZMod q) := by
  letI : Fact p.Prime := ⟨hp⟩
  letI : Fact q.Prime := ⟨hq⟩
  have hq' : ((-(q : ℤ)) : ZMod p) ≠ 0 := by
    simpa using ZMod.prime_ne_zero p q hpq
  have h₁ : legendreSym p (-(q : ℤ)) = 1 ↔ _ :=
    legendreSym.eq_one_iff p (a := -(q : ℤ)) (by simpa only [Int.cast_neg, Int.cast_natCast] using hq')
  have h₂ := legendreSym.eq_one_iff' q (ZMod.prime_ne_zero q p hpq.symm)
  rw [legendre_negative_prime hp₂ hq₄] at h₁
  simpa only [Int.cast_neg, Int.cast_natCast] using h₁.symm.trans h₂

private theorem sqrt3_residues (p : ℕ) : IsSquare (p : ZMod 3) ↔ p % 3 = 0 ∨ p % 3 = 1 := by
  have table : ∀ x : ZMod 3, IsSquare x ↔ x.val = 0 ∨ x.val = 1 := by
    decide
  simpa only [ZMod.val_natCast] using table p

private theorem sqrt7_residues (p : ℕ) : IsSquare (p : ZMod 7) ↔
    p % 7 = 0 ∨ p % 7 = 1 ∨ p % 7 = 2 ∨ p % 7 = 4 := by
  have table : ∀ x : ZMod 7, IsSquare x ↔
      x.val = 0 ∨ x.val = 1 ∨ x.val = 2 ∨ x.val = 4 := by
    decide
  simpa only [ZMod.val_natCast] using table p

private theorem sqrt11_residues (p : ℕ) : IsSquare (p : ZMod 11) ↔
    p % 11 = 0 ∨ p % 11 = 1 ∨ p % 11 = 3 ∨ p % 11 = 4 ∨ p % 11 = 5 ∨ p % 11 = 9 := by
  have table : ∀ x : ZMod 11, IsSquare x ↔
      x.val = 0 ∨ x.val = 1 ∨ x.val = 3 ∨ x.val = 4 ∨ x.val = 5 ∨ x.val = 9 := by
    decide
  simpa only [ZMod.val_natCast] using table p

private theorem prime_mod_ne_zero {p q : ℕ} (hp : p.Prime) (hq : q.Prime)
    (hne : p ≠ q) : p % q ≠ 0 := by
  intro h
  obtain h₁ | h₂ := (Nat.dvd_prime hp).mp (Nat.dvd_of_mod_eq_zero h)
  · exact hq.ne_one h₁
  · exact hne h₂.symm

private theorem root3_prime {p : ℕ} (hp : p.Prime) (hne : p ≠ 3) :
    HasNormRoot (-1) 1 p ↔ p % 3 = 1 := by
  by_cases h₂ : p = 2
  · subst p
    unfold HasNormRoot
    decide
  rw [hasNormRoot_iff_square_discriminant hp h₂]
  rw [show (-1 : ℤ) ^ 2 - 4 * 1 = -(3 : ℤ) by decide]
  have hs := square_negative_prime_iff hp (by decide : Nat.Prime 3) h₂ (by decide) hne
  rw [sqrt3_residues] at hs
  simpa [prime_mod_ne_zero hp (by decide) hne] using hs

private theorem root7_prime {p : ℕ} (hp : p.Prime) (hne : p ≠ 7) :
    HasNormRoot 1 2 p ↔ p % 7 = 1 ∨ p % 7 = 2 ∨ p % 7 = 4 := by
  by_cases h₂ : p = 2
  · subst p
    unfold HasNormRoot
    decide
  rw [hasNormRoot_iff_square_discriminant hp h₂]
  rw [show (1 : ℤ) ^ 2 - 4 * 2 = -(7 : ℤ) by decide]
  have hs := square_negative_prime_iff hp (by decide : Nat.Prime 7) h₂ (by decide) hne
  rw [sqrt7_residues] at hs
  simpa [prime_mod_ne_zero hp (by decide) hne] using hs

private theorem root8_prime {p : ℕ} (hp : p.Prime) (hne : p ≠ 2) :
    HasNormRoot 0 2 p ↔ p % 8 = 1 ∨ p % 8 = 3 := by
  letI : Fact p.Prime := ⟨hp⟩
  unfold HasNormRoot
  simp only [Int.cast_zero, zero_mul, sub_zero, Int.cast_ofNat]
  rw [← ZMod.exists_sq_eq_neg_two_iff hne]
  constructor
  · rintro ⟨t, ht⟩
    exact ⟨t, by linear_combination -ht⟩
  · rintro ⟨t, ht⟩
    exact ⟨t, by linear_combination -ht⟩

theorem hasNormRoot11_prime_iff {p : ℕ} (hp : p.Prime) (hne : p ≠ 11) :
    HasNormRoot 1 3 p ↔ p % 11 = 1 ∨ p % 11 = 3 ∨ p % 11 = 4 ∨
      p % 11 = 5 ∨ p % 11 = 9 := by
  by_cases h₂ : p = 2
  · subst p
    unfold HasNormRoot
    decide
  rw [hasNormRoot_iff_square_discriminant hp h₂]
  rw [show (1 : ℤ) ^ 2 - 4 * 3 = -(11 : ℤ) by decide]
  have hs := square_negative_prime_iff hp (by decide : Nat.Prime 11) h₂ (by decide) hne
  rw [sqrt11_residues] at hs
  simpa [prime_mod_ne_zero hp (by decide) hne] using hs

private theorem prime_not_dvd_neg_prime {p r : ℕ} (hp : p.Prime) (hr : r.Prime) (h : p ≠ r) :
    ¬ (p : ℤ) ∣ -(r : ℤ) := by
  intro hd
  have hd' : p ∣ r := by
    have := Int.dvd_neg.mpr hd
    simpa only [neg_neg, Int.natCast_dvd_natCast] using this
  exact h ((Nat.dvd_prime hr).mp hd' |>.resolve_left hp.ne_one)

private theorem prime_not_dvd_neg_eight {p : ℕ} (hp : p.Prime) (h : p ≠ 2) :
    ¬ (p : ℤ) ∣ (-8 : ℤ) := by
  intro hd
  have hd' : p ∣ 2 ^ 3 := Int.natCast_dvd_natCast.mp
    (Int.dvd_neg.mpr hd : (p : ℤ) ∣ (↑(2 ^ 3 : ℕ) : ℤ))
  exact h ((Nat.dvd_prime Nat.prime_two).mp (hp.dvd_of_dvd_pow hd') |>.resolve_left hp.ne_one)

/-- A simple root at each core prime lifts to every power. The ramified
prime has a root at exponent one only. -/
private theorem primitive_iff_primeCore_or_dilate {b c : ℤ} {r : ℕ}
    {allowed : ℕ → Prop} (hp : r.Prime) (hr : ¬ allowed r)
    (hdisc : 0 < 4 * c - b ^ 2) (hdisc' : 4 * c - b ^ 2 ≤ 11)
    (hprime : ∀ p, p.Prime → p ≠ r → (HasNormRoot b c p ↔ allowed p))
    (hgood : ∀ p, p.Prime → p ≠ r → ¬ (p : ℤ) ∣ b ^ 2 - 4 * c)
    (hrRoot : HasNormRoot b c r) (hrSquare : ¬ HasNormRoot b c (r ^ 2))
    (n : ℕ) : PrimitiveValues b c n ↔ PrimeCore allowed n ∨ Dilate r (PrimeCore allowed) n := by
  rw [primitiveValues_iff_normRoot hdisc hdisc']
  constructor
  · rintro ⟨hn, hroot⟩
    have hcore (m : ℕ) (hm : m ∣ n) (hrm : ¬ r ∣ m) : PrimeCore allowed m := by
      refine ⟨Nat.pos_of_dvd_of_pos hm hn, fun p hpp hpm => ?_⟩
      have hpne : p ≠ r := fun he => hrm (he ▸ hpm)
      exact (hprime p hpp hpne).mp (hroot.of_dvd (hpm.trans hm))
    by_cases hrd : r ∣ n
    · obtain ⟨m, rfl⟩ := hrd
      refine Or.inr ⟨m, hcore m (dvd_mul_left ..) ?_, rfl⟩
      intro hrm
      exact hrSquare (hroot.of_dvd (by
        rw [pow_two]
        exact Nat.mul_dvd_mul_left _ hrm))
    · exact Or.inl (hcore n dvd_rfl hrd)
  · have hcoreRoot {m : ℕ} (hm : PrimeCore allowed m) : HasNormRoot b c m := by
      apply (hasNormRoot_iff_primePowers hm.1).mpr
      intro p k hpp hk hpk
      have hpa := hm.2 p hpp ((dvd_pow_self _ hk.ne').trans hpk)
      have hpne : p ≠ r := fun he => hr (he ▸ hpa)
      exact HasNormRoot.of_unramified hpp (hgood p hpp hpne)
        ((hprime p hpp hpne).mpr hpa) k
    rintro (h | ⟨m, hm, rfl⟩)
    · exact ⟨h.1, hcoreRoot h⟩
    · refine ⟨Nat.mul_pos hp.pos hm.1, HasNormRoot.mul ?_ hrRoot (hcoreRoot hm)⟩
      exact hp.coprime_iff_not_dvd.mpr (fun h => hr (hm.2 r hp h))

theorem P3_iff_core (n : ℕ) : P3 n ↔ C3 n ∨ Dilate 3 C3 n :=
  primitive_iff_primeCore_or_dilate (by decide) (by decide) (by decide) (by decide)
    (fun _ => root3_prime) (fun _ hp hne => prime_not_dvd_neg_prime hp (by decide) hne)
    (by unfold HasNormRoot; decide)
    (not_hasNormRoot_ramified_sq (by decide) (by decide)) n

theorem P7_iff_core (n : ℕ) : P7 n ↔ C7 n ∨ Dilate 7 C7 n :=
  primitive_iff_primeCore_or_dilate (by decide) (by decide) (by decide) (by decide)
    (fun _ => root7_prime) (fun _ hp hne => prime_not_dvd_neg_prime hp (by decide) hne)
    (by unfold HasNormRoot; decide)
    (not_hasNormRoot_ramified_sq (by decide) (by decide)) n

theorem P8_iff_core (n : ℕ) : P8 n ↔ C8 n ∨ Dilate 2 C8 n :=
  primitive_iff_primeCore_or_dilate (by decide) (by decide) (by decide) (by decide)
    (fun _ => root8_prime) (fun _ => prime_not_dvd_neg_eight)
    (by unfold HasNormRoot; decide)
    (by unfold HasNormRoot; decide) n

theorem P11_iff_core (n : ℕ) : P11 n ↔ C11 n ∨ Dilate 11 C11 n :=
  primitive_iff_primeCore_or_dilate (by decide) (by decide) (by decide) (by decide)
    (fun _ => hasNormRoot11_prime_iff) (fun _ hp hne => prime_not_dvd_neg_prime hp (by decide) hne)
    (by unfold HasNormRoot; decide)
    (not_hasNormRoot_ramified_sq (by decide) (by decide)) n

theorem restoration3 : PrimitiveRestoration C3 P3 3 :=
  ⟨fun n hn => (P3_iff_core n).mpr (Or.inl hn),
    fun n hn => (P3_iff_core _).mpr (Or.inr ⟨n, hn, rfl⟩)⟩

theorem restoration7 : PrimitiveRestoration C7 P7 7 :=
  ⟨fun n hn => (P7_iff_core n).mpr (Or.inl hn),
    fun n hn => (P7_iff_core _).mpr (Or.inr ⟨n, hn, rfl⟩)⟩

theorem restoration8 : PrimitiveRestoration C8 P8 2 :=
  ⟨fun n hn => (P8_iff_core n).mpr (Or.inl hn),
    fun n hn => (P8_iff_core _).mpr (Or.inr ⟨n, hn, rfl⟩)⟩

theorem restoration11 : PrimitiveRestoration C11 P11 11 :=
  ⟨fun n hn => (P11_iff_core n).mpr (Or.inl hn),
    fun n hn => (P11_iff_core _).mpr (Or.inr ⟨n, hn, rfl⟩)⟩

end QuadraticNorms
