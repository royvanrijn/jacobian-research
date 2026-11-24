import QuadraticNorms.PrimitiveClassification
import Mathlib.NumberTheory.SumTwoSquares

namespace QuadraticNorms

/-- Positive ordinary values with no gcd restriction on the coordinates. -/
def OrdinaryValues (b c : ℤ) (n : ℕ) : Prop :=
  0 < n ∧ ∃ x y : ℤ, form b c x y = n

theorem ordinaryValues_iff_square_times_primitive {b c : ℤ} {n : ℕ} :
    OrdinaryValues b c n ↔ ∃ g m : ℕ, 0 < g ∧ PrimitiveValues b c m ∧ n = g ^ 2 * m := by
  constructor
  · rintro ⟨hn, x, y, hform⟩
    have hpos : 0 < Int.gcd x y := by
      by_contra! h
      obtain ⟨rfl, rfl⟩ := Int.gcd_eq_zero_iff.mp (by omega : Int.gcd x y = 0)
      simp [form] at hform
      omega
    obtain ⟨g, u, v, hg, hgcd, hxu, hyv⟩ := Int.exists_gcd_one' hpos
    have heq : (n : ℤ) = (g : ℤ) ^ 2 * form b c u v := by
      rw [← hform, hxu, hyv]
      unfold form
      ring
    have hvalue : 0 < form b c u v := by
      apply pos_of_mul_pos_right (a := (g : ℤ) ^ 2) ?_ (sq_nonneg _)
      rw [← heq]
      exact_mod_cast hn
    let m := (form b c u v).toNat
    have hcast : (m : ℤ) = form b c u v := Int.toNat_of_nonneg hvalue.le
    refine ⟨g, m, hg, ⟨by omega, u, v, hgcd, hcast.symm⟩, ?_⟩
    rw [← hcast] at heq
    exact_mod_cast heq
  · rintro ⟨g, m, hg, ⟨hm, x, y, _, hform⟩, rfl⟩
    refine ⟨Nat.mul_pos (pow_pos hg 2) hm, x * g, y * g, ?_⟩
    rw [Int.natCast_mul, Int.natCast_pow]
    calc
      _ = (g : ℤ) ^ 2 * form b c x y := by unfold form; ring
      _ = _ := by rw [hform]

/-- The inert-prime valuation criterion used by the Gaussian indicator. -/
theorem ordinaryGaussian_iff {n : ℕ} :
    OrdinaryValues 0 1 n ↔ 0 < n ∧ ∀ p, p.Prime → p % 4 = 3 → Even (padicValNat p n) := by
  constructor
  · rintro ⟨hn, x, y, hform⟩
    have hxy : n = x.natAbs ^ 2 + y.natAbs ^ 2 := by
      have hr : (n : ℤ) = |x| ^ 2 + |y| ^ 2 := by simpa [form] using hform.symm
      rw [← Int.natCast_natAbs, ← Int.natCast_natAbs] at hr
      exact_mod_cast hr
    refine ⟨hn, fun p hp hmod => ?_⟩
    by_cases hd : p ∣ n
    · exact (Nat.eq_sq_add_sq_iff.mp ⟨_, _, hxy⟩) p
        (Nat.mem_primeFactors.mpr ⟨hp, hd, hn.ne'⟩) hmod
    · rw [← Nat.factorization_def _ hp, Nat.factorization_eq_zero_of_not_dvd hd]
      exact Even.zero
  · rintro ⟨hn, hv⟩
    obtain ⟨x, y, hxy⟩ := Nat.eq_sq_add_sq_iff.mpr
      (fun p hp hmod => hv p (Nat.prime_of_mem_primeFactors hp) hmod)
    refine ⟨hn, x, y, ?_⟩
    have h : (n : ℤ) = (x : ℤ) ^ 2 + (y : ℤ) ^ 2 := by exact_mod_cast hxy
    simpa [form] using h.symm

private theorem squarefree_in_primitiveCore {m r : ℕ} {allowed : ℕ → Prop}
    (hm : 0 < m) (hs : Squarefree m) (hr : r.Prime)
    (h : ∀ p, p.Prime → p ∣ m → p ≠ r → allowed p) :
    PrimeCore allowed m ∨ Dilate r (PrimeCore allowed) m := by
  by_cases hd : r ∣ m
  · obtain ⟨k, rfl⟩ := hd
    refine Or.inr ⟨k, ⟨Nat.pos_of_mul_pos_left hm, ?_⟩, rfl⟩
    intro p hp hpk
    apply h p hp (hpk.trans (dvd_mul_left ..))
    intro he
    subst p
    exact hr.not_isUnit (hs r (Nat.mul_dvd_mul_left r hpk))
  · exact Or.inl ⟨hm, fun p hp hpm => h p hp hpm (fun he => hd (he ▸ hpm))⟩

theorem ordinary_iff_even_inert_valuations {b c : ℤ} {r : ℕ} {allowed : ℕ → Prop}
    (hr : r.Prime)
    (hP : ∀ n, PrimitiveValues b c n ↔ PrimeCore allowed n ∨ Dilate r (PrimeCore allowed) n)
    {n : ℕ} :
    OrdinaryValues b c n ↔ 0 < n ∧ ∀ p, p.Prime → p ≠ r → ¬ allowed p → Even (padicValNat p n) := by
  constructor
  · intro h
    have hn := h.1
    obtain ⟨g, m, hg, hm, heq⟩ := ordinaryValues_iff_square_times_primitive.mp h
    have hmpos := hm.1
    refine ⟨hn, fun p hp hpne hmod => ?_⟩
    have hnot (k : ℕ) (hk : PrimeCore allowed k) : ¬ p ∣ k := fun hd => hmod (hk.2 p hp hd)
    have hpnot : ¬ p ∣ m := by
      rcases (hP m).mp hm with hcore | ⟨k, hk, rfl⟩
      · exact hnot _ hcore
      · intro hd
        rcases hp.dvd_mul.mp hd with hd₂ | hdk
        · exact hpne ((Nat.dvd_prime hr).mp hd₂ |>.resolve_left hp.ne_one)
        · exact hnot _ hk hdk
    rw [← Nat.factorization_def _ hp, heq, Nat.factorization_mul (pow_ne_zero _ hg.ne') hmpos.ne',
      Nat.factorization_pow, Finsupp.add_apply, Finsupp.smul_apply,
      Nat.factorization_eq_zero_of_not_dvd hpnot, add_zero, nsmul_eq_mul]
    exact ⟨g.factorization p, by omega⟩
  · rintro ⟨hn, hv⟩
    obtain ⟨m, g, hm, hg, heq, hs⟩ := Nat.sq_mul_squarefree_of_pos hn
    refine ordinaryValues_iff_square_times_primitive.mpr ⟨g, m, hg, (hP m).mpr ?_, heq.symm⟩
    have hsupport : ∀ p, p.Prime → p ∣ m → p ≠ r → allowed p := by
      intro p hp hpd hpne
      by_contra! h
      have hev := hv p hp hpne h
      rw [← Nat.factorization_def _ hp, ← heq, Nat.factorization_mul (pow_ne_zero _ hg.ne') hm.ne',
        Nat.factorization_pow, Finsupp.add_apply, Finsupp.smul_apply,
        Nat.factorization_eq_one_of_squarefree hs hp hpd, nsmul_eq_mul] at hev
      obtain ⟨j, hj⟩ := hev
      omega
    exact squarefree_in_primitiveCore hm hs hr hsupport

/-- The inert-prime valuation criterion used by the ordinary
x²+2y² indicator. The ramified prime two has unrestricted ordinary
valuation because even powers may occur in the removed square factor. -/
theorem ordinaryTwo_iff {n : ℕ} :
    OrdinaryValues 0 2 n ↔ 0 < n ∧
      ∀ p, p.Prime → (p % 8 = 5 ∨ p % 8 = 7) → Even (padicValNat p n) := by
  rw [ordinary_iff_even_inert_valuations (r := 2)
    (allowed := fun p => p % 8 = 1 ∨ p % 8 = 3) (by decide) P8_iff_core]
  refine and_congr_right (fun _ => ⟨?_, ?_⟩)
  · exact fun h p hp hmod => h p hp (by omega) (by omega)
  · intro h p hp hpne hn
    apply h p hp
    have hpodd := hp.eq_two_or_odd.resolve_left hpne
    have hmod := Nat.mod_mod_of_dvd p (by decide : 2 ∣ 8)
    rw [hpodd] at hmod
    omega

end QuadraticNorms
