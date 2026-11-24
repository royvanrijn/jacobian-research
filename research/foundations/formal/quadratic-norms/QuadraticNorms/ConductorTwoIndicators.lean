import QuadraticNorms.OrdinaryIndicators

namespace QuadraticNorms

private theorem eisenstein_even_second_coordinate (x y : ℤ) :
    ∃ u v : ℤ, (2 : ℤ) ∣ v ∧ form (-1) 1 u v = form (-1) 1 x y := by
  by_cases hy : (2 : ℤ) ∣ y
  · exact ⟨x, y, hy, rfl⟩
  by_cases hx : (2 : ℤ) ∣ x
  · exact ⟨y, x, hx, by unfold form; ring⟩
  · exact ⟨x, x - y, by omega, by unfold form; ring⟩

/-- For x²+3y² the conductor-two restriction leaves all ordinary
Eisenstein values available. -/
theorem ordinaryThree_iff_eisenstein {n : ℕ} : OrdinaryValues 0 3 n ↔ OrdinaryValues (-1) 1 n := by
  constructor
  · rintro ⟨hn, x, y, hxy⟩
    refine ⟨hn, x + y, 2 * y, ?_⟩
    rw [← hxy]
    unfold form
    ring
  · rintro ⟨hn, x, y, hxy⟩
    obtain ⟨u, _, ⟨v, rfl⟩, h⟩ := eisenstein_even_second_coordinate x y
    refine ⟨hn, u - v, v, ?_⟩
    rw [← hxy, ← h]
    unfold form
    ring

private theorem gaussian_even_second_coordinate_of_odd {x y : ℤ} {n : ℕ}
    (h : form 0 1 x y = n) (hodd : n % 2 = 1) :
    ∃ u v : ℤ, (2 : ℤ) ∣ v ∧ form 0 1 u v = n := by
  by_cases hy : (2 : ℤ) ∣ y
  · exact ⟨x, y, hy, h⟩
  refine ⟨y, x, ?_, by rw [← h]; unfold form; ring⟩
  have hm := (Int.mod_modEq x 2).pow 2 |>.add ((Int.mod_modEq y 2).pow 2)
  have heq : x ^ 2 + y ^ 2 = (n : ℤ) := by simpa [form] using h
  rw [heq] at hm
  have hymod : y % 2 = 1 := by omega
  by_contra! hx
  have hxmod : x % 2 = 1 := by omega
  rw [hxmod, hymod] at hm
  change (0 : ℤ) = (n : ℤ) % 2 at hm
  omega

theorem ordinaryFour_iff_odd_or_four_dilate {n : ℕ} :
    OrdinaryValues 0 4 n ↔
      (n % 2 = 1 ∧ OrdinaryValues 0 1 n) ∨ Dilate 4 (OrdinaryValues 0 1) n := by
  constructor
  · rintro ⟨hn, x, y, hform⟩
    have hxy : form 0 1 x (2 * y) = n := by
      rw [← hform]
      unfold form
      ring
    by_cases ho : n % 2 = 1
    · exact Or.inl ⟨ho, hn, x, 2 * y, hxy⟩
    have h_even : (n : ℤ) % 2 = 0 := by omega
    have hmod := (Int.mod_modEq x 2).pow 2
    have hxmod : x % 2 = 0 := by
      by_contra! hh
      have hxone : x % 2 = 1 := by omega
      rw [hxone] at hmod
      have hquad : x ^ 2 + 4 * y ^ 2 = (n : ℤ) := by simpa [form] using hform
      change (1 : ℤ) = x ^ 2 % 2 at hmod
      have heq : (n : ℤ) % 2 = x ^ 2 % 2 := by omega
      omega
    obtain ⟨u, rfl⟩ := Int.dvd_of_emod_eq_zero hxmod
    let m := (u ^ 2 + y ^ 2).toNat
    have hm : (m : ℤ) = u ^ 2 + y ^ 2 := Int.toNat_of_nonneg (by positivity)
    have heq : (n : ℤ) = 4 * (m : ℤ) := by
      rw [← hform, hm]
      unfold form
      ring
    exact Or.inr ⟨m, ⟨by omega, u, y, by simpa [form] using hm.symm⟩,
      by exact_mod_cast heq.symm⟩
  · rintro (⟨ho, hn, x, y, hxy⟩ | ⟨m, ⟨hm, x, y, hxy⟩, rfl⟩)
    · obtain ⟨u, _, ⟨v, rfl⟩, h⟩ := gaussian_even_second_coordinate_of_odd hxy ho
      refine ⟨hn, u, v, ?_⟩
      rw [← h]
      unfold form
      ring
    · refine ⟨by omega, 2 * x, y, ?_⟩
      push_cast
      rw [← hxy]
      unfold form
      ring


theorem ordinaryThree_iff {n : ℕ} :
    OrdinaryValues 0 3 n ↔ 0 < n ∧ ∀ p, p.Prime → p % 3 = 2 → Even (padicValNat p n) := by
  rw [ordinaryThree_iff_eisenstein, ordinary_iff_even_inert_valuations
    (r := 3) (allowed := fun p => p % 3 = 1) (by decide) P3_iff_core]
  refine and_congr_right (fun _ => ⟨?_, ?_⟩)
  · exact fun h p hp hmod => h p hp (by omega) (by omega)
  · intro h p hp hpne hr
    apply h p hp
    have hpmod : p % 3 ≠ 0 := by
      intro h
      exact hpne (((Nat.dvd_prime hp).mp (Nat.dvd_of_mod_eq_zero h)).resolve_left (by decide)).symm
    omega

private theorem valuation_eq_at_four_mul {p m : ℕ} (hp : p.Prime) (hpne : p ≠ 2)
    (hm : 0 < m) : padicValNat p (4 * m) = padicValNat p m := by
  have hp4 : ¬ p ∣ 4 := by
    intro hd
    have hdiv : p ∣ 2 ^ 2 := hd
    exact hpne ((Nat.dvd_prime Nat.prime_two).mp (hp.dvd_of_dvd_pow hdiv) |>.resolve_left hp.ne_one)
  rw [← Nat.factorization_def _ hp, ← Nat.factorization_def _ hp,
    Nat.factorization_mul (by decide : 4 ≠ 0) hm.ne', Finsupp.add_apply,
    Nat.factorization_eq_zero_of_not_dvd hp4, zero_add]

private theorem unused_two_valuation_iff {n : ℕ} (hn : 0 < n) :
    padicValNat 2 n ≠ 1 ↔ n % 2 = 1 ∨ 4 ∣ n := by
  rw [← Nat.factorization_def _ Nat.prime_two]
  have htwo := Nat.prime_two.dvd_iff_one_le_factorization hn.ne'
  have hfour := Nat.prime_two.pow_dvd_iff_le_factorization (k := 2) hn.ne'
  change 4 ∣ n ↔ _ at hfour
  rw [hfour]
  omega

theorem ordinaryFour_iff {n : ℕ} :
    OrdinaryValues 0 4 n ↔ 0 < n ∧ padicValNat 2 n ≠ 1 ∧
      ∀ p, p.Prime → p % 4 = 3 → Even (padicValNat p n) := by
  rw [ordinaryFour_iff_odd_or_four_dilate]
  constructor
  · rintro (⟨hodd, hG⟩ | ⟨m, hG, rfl⟩)
    · have ⟨hn, hp⟩ := ordinaryGaussian_iff.mp hG
      exact ⟨hn, (unused_two_valuation_iff hn).mpr (Or.inl hodd), hp⟩
    · have ⟨hm, h⟩ := ordinaryGaussian_iff.mp hG
      have hn : 0 < 4 * m := by omega
      refine ⟨hn, (unused_two_valuation_iff hn).mpr (Or.inr (dvd_mul_right ..)), ?_⟩
      intro p hp hmod
      rw [valuation_eq_at_four_mul hp (by omega) hm]
      exact h p hp hmod
  · rintro ⟨hn, hv, h⟩
    rcases (unused_two_valuation_iff hn).mp hv with hodd | ⟨m, rfl⟩
    · exact Or.inl ⟨hodd, ordinaryGaussian_iff.mpr ⟨hn, h⟩⟩
    · have hm : 0 < m := by omega
      refine Or.inr ⟨m, ordinaryGaussian_iff.mpr ⟨hm, ?_⟩, rfl⟩
      intro p hp hmod
      rw [← valuation_eq_at_four_mul hp (by omega) hm]
      exact h p hp hmod

private theorem normSeven_ordinary_iff {n : ℕ} :
    OrdinaryValues 1 2 n ↔ 0 < n ∧
      ∀ p, p.Prime → (p % 7 = 3 ∨ p % 7 = 5 ∨ p % 7 = 6) → Even (padicValNat p n) := by
  rw [ordinary_iff_even_inert_valuations (r := 7) (allowed := fun p => p % 7 = 1 ∨ p % 7 = 2 ∨ p % 7 = 4)
    (by decide) P7_iff_core]
  refine and_congr_right (fun _ => ⟨?_, ?_⟩)
  · exact fun h p hp hmod => h p hp (by omega) (by omega)
  · intro h p hp hpne hr
    apply h p hp
    have hpmod : p % 7 ≠ 0 := by
      intro h
      exact hpne (((Nat.dvd_prime hp).mp (Nat.dvd_of_mod_eq_zero h)).resolve_left (by decide)).symm
    omega

theorem ordinarySeven_iff {n : ℕ} :
    OrdinaryValues 0 7 n ↔ 0 < n ∧ padicValNat 2 n ≠ 1 ∧
      ∀ p, p.Prime → (p % 7 = 3 ∨ p % 7 = 5 ∨ p % 7 = 6) → Even (padicValNat p n) := by
  constructor
  · rintro ⟨hn, x, y, hxy⟩
    have hbase : OrdinaryValues 1 2 n := by
      refine ⟨hn, x - y, 2 * y, ?_⟩
      rw [← hxy]
      unfold form
      ring
    refine ⟨hn, (unused_two_valuation_iff hn).mpr ?_, (normSeven_ordinary_iff.mp hbase).2⟩
    have hresx : x % 4 = 0 ∨ x % 4 = 1 ∨ x % 4 = 2 ∨ x % 4 = 3 := by omega
    have hresy : y % 4 = 0 ∨ y % 4 = 1 ∨ y % 4 = 2 ∨ y % 4 = 3 := by omega
    have hmod := ((Int.mod_modEq x 4).pow 2).add ((Int.ModEq.refl 7).mul ((Int.mod_modEq y 4).pow 2))
    have hform : x ^ 2 + 7 * y ^ 2 = (n : ℤ) := by simpa [form] using hxy
    rw [hform] at hmod
    change ((x % 4) ^ 2 + 7 * (y % 4) ^ 2) % 4 = (n : ℤ) % 4 at hmod
    rcases hresx with hx | hx | hx | hx <;> rcases hresy with hy | hy | hy | hy <;>
      rw [hx, hy] at hmod <;> norm_num at hmod <;> omega
  · rintro ⟨hn, hv, h⟩
    rcases (unused_two_valuation_iff hn).mp hv with hodd | ⟨m, rfl⟩
    · obtain ⟨_, x, y, hxy⟩ := normSeven_ordinary_iff.mpr ⟨hn, h⟩
      have hform : x ^ 2 + x * y + 2 * y ^ 2 = (n : ℤ) := by simpa [form] using hxy
      have hyeven : (2 : ℤ) ∣ y := by
        have hx := Int.mod_modEq x 2
        have hy := Int.mod_modEq y 2
        have hmod := ((hx.pow 2).add (hx.mul hy)).add ((Int.ModEq.refl 2).mul (hy.pow 2))
        rw [hform] at hmod
        by_contra! hny
        have hyone : y % 2 = 1 := by omega
        have hxzero : x % 2 = 0 ∨ x % 2 = 1 := by omega
        change ((x % 2) ^ 2 + (x % 2) * (y % 2) + 2 * (y % 2) ^ 2) % 2 = (n : ℤ) % 2 at hmod
        rcases hxzero with hxzero | hxzero <;> rw [hxzero, hyone] at hmod <;>
          norm_num at hmod <;> omega
      obtain ⟨v, rfl⟩ := hyeven
      refine ⟨hn, x + v, v, ?_⟩
      rw [← hxy]
      unfold form
      ring
    · have hm : 0 < m := by omega
      have hmvalues : OrdinaryValues 1 2 m := by
        refine normSeven_ordinary_iff.mpr ⟨hm, ?_⟩
        intro p hp hmod
        rw [← valuation_eq_at_four_mul hp (by omega) hm]
        exact h p hp hmod
      obtain ⟨_, x, y, hxy⟩ := hmvalues
      refine ⟨hn, 2 * x + y, y, ?_⟩
      push_cast
      rw [← hxy]
      unfold form
      ring

/-- At a shift congruent to two modulo four, the ordinary x²+4y²
indicator has identically zero pair count. -/
theorem ordinaryFour_no_pairs_at_shift_two_mod_four {a b : ℕ} {h : ℤ}
    (hh : h % 4 = 2) (heq : (a : ℤ) + h = b) (ha : OrdinaryValues 0 4 a) :
    ¬ OrdinaryValues 0 4 b := by
  have residues {n : ℕ} (hn : OrdinaryValues 0 4 n) : n % 4 = 0 ∨ n % 4 = 1 := by
    obtain ⟨_, x, y, hx⟩ := hn
    have hxy : x ^ 2 + 4 * y ^ 2 = (n : ℤ) := by simpa [form] using hx
    have hxmod := (Int.mod_modEq x 4).pow 2
    have hmod : (x % 4) ^ 2 % 4 = (n : ℤ) % 4 := by
      change (x % 4) ^ 2 % 4 = x ^ 2 % 4 at hxmod
      omega
    have hcases : x % 4 = 0 ∨ x % 4 = 1 ∨ x % 4 = 2 ∨ x % 4 = 3 := by omega
    rcases hcases with hh | hh | hh | hh <;> rw [hh] at hmod <;> norm_num at hmod <;> omega
  intro hb
  have := residues ha
  have := residues hb
  omega

end QuadraticNorms
