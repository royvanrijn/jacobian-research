import QuadraticNorms.MinusElevenPrimitiveCounts
import Mathlib.NumberTheory.ArithmeticFunction.Defs

open scoped BigOperators

namespace QuadraticNorms

private theorem mem_minusEleven_roots {n : ℕ} [NeZero n] {t : ZMod n} :
    t ∈ minusElevenRoots n ↔ t ^ 2 - t + 3 = 0 := by
  simpa only [minusElevenRoots, Int.cast_one, one_mul, Int.cast_ofNat] using
    (mem_normRootResidues (b := 1) (c := 3) (t := t))

/-- The root count at positive moduli, extended by zero at modulus zero. -/
noncomputable def minusElevenRootCount : ArithmeticFunction ℕ :=
  normRootCount 1 3

theorem minusEleven_root_count_apply {n : ℕ} [NeZero n] :
    minusElevenRootCount n = (minusElevenRoots n).card :=
  norm_root_count_apply 1 3 n

theorem minusEleven_root_count_of_no_root {n : ℕ} (hn : 0 < n) (h : ¬ HasNormRoot 1 3 n) :
    minusElevenRootCount n = 0 := norm_root_count_of_no_root hn h

@[simp] theorem minusEleven_root_count_one : minusElevenRootCount 1 = 1 :=
  norm_root_count_one 1 3

theorem minusEleven_root_count_multiplicative : minusElevenRootCount.IsMultiplicative :=
  norm_root_count_multiplicative 1 3

theorem minusEleven_unramified_root_count_le_two {p k : ℕ} (hp : p.Prime) (hne : p ≠ 11)
    (hk : 0 < k) : minusElevenRootCount (p ^ k) ≤ 2 := by
  apply norm_root_count_unramified_le_two 1 3 hp _ hk
  intro h
  have hd : p ∣ 11 := Int.natCast_dvd_natCast.mp (by norm_num at h ⊢; exact h)
  exact hne ((Nat.dvd_prime (by decide : Nat.Prime 11)).mp hd |>.resolve_left hp.ne_one)

theorem minusEleven_ramified_root_count_le_one : minusElevenRootCount 11 ≤ 1 := by
  rw [minusEleven_root_count_apply]
  apply Finset.card_le_one.mpr
  have h : ∀ a : ZMod 11, a ^ 2 - a + 3 = 0 → a = 6 := by decide
  exact fun a ha b hb => (h a (mem_minusEleven_roots.mp ha)).trans (h b (mem_minusEleven_roots.mp hb)).symm

theorem minusEleven_ramified_root_count_eq_zero {k : ℕ} (hk : 2 ≤ k) :
    minusElevenRootCount (11 ^ k) = 0 := by
  refine minusEleven_root_count_of_no_root (by positivity) (fun h => ?_)
  exact not_hasNormRoot_ramified_sq (b := 1) (c := 3) (by decide) (by decide)
    (h.of_dvd (pow_dvd_pow 11 hk))

theorem minusEleven_nonsplit_root_count_eq_zero {p k : ℕ} (hp : p.Prime) (hne : p ≠ 11)
    (hs : ¬ MinusElevenSplitPrime p) (hk : 0 < k) : minusElevenRootCount (p ^ k) = 0 := by
  refine minusEleven_root_count_of_no_root (pow_pos hp.pos _) (fun h => hs ⟨hp, ?_⟩)
  have hh := (hasNormRoot11_prime_iff hp hne).mp (h.of_dvd (dvd_pow_self p hk.ne'))
  simpa only [List.mem_cons, List.not_mem_nil, or_false] using hh

end QuadraticNorms
