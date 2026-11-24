import QuadraticNorms.PrimeCores
import Mathlib.NumberTheory.ArithmeticFunction.Moebius
import Mathlib.Data.ZMod.Basic

open scoped BigOperators ArithmeticFunction.Moebius

namespace QuadraticNorms

/-- The residue pairs counted by Γ_N in the subgroup-sum theorem. -/
def residuePairs (m A B n : ℕ) (U V : Finset (Fin m)) : Finset (Fin m × Fin m) :=
  (U ×ˢ V).filter fun uv => (A * uv.1.val + B * uv.2.val) % m = n % m

/-- The exact congruence obstruction in the subgroup-sum theorem. -/
theorem residuePairs_nonempty_of_representation {m A B n x y : ℕ}
    {U V : Finset (Fin m)} (hm : 0 < m) (h : A * x + B * y = n)
    (hx : (⟨x % m, Nat.mod_lt x hm⟩ : Fin m) ∈ U)
    (hy : (⟨y % m, Nat.mod_lt y hm⟩ : Fin m) ∈ V) :
    (residuePairs m A B n U V).Nonempty := by
  refine ⟨(⟨x % m, Nat.mod_lt x hm⟩, ⟨y % m, Nat.mod_lt y hm⟩), ?_⟩
  refine Finset.mem_filter.mpr ⟨Finset.mem_product.mpr ⟨hx, hy⟩, ?_⟩
  exact (((Nat.ModEq.refl A).mul (Nat.mod_modEq x m)).add
    ((Nat.ModEq.refl B).mul (Nat.mod_modEq y m))).trans
    (show A * x + B * y ≡ n [MOD m] by rw [h])

theorem C12_sum_mod_twelve {n : ℕ} (h : InSumset C12 C12 n) : n % 12 = 2 := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have := C12_mod_twelve ha
  have := C12_mod_twelve hb
  omega

/-- The good-prime local count when both summands require units.
Its one or two forbidden residues are zero and the target. -/
theorem goodPrime_unit_pair_count (p : ℕ) [NeZero p] (N : ZMod p) :
    (Finset.univ.filter fun a : ZMod p => a ≠ 0 ∧ N - a ≠ 0).card =
      p - (if N = 0 then 1 else 2) := by
  classical
  have heq : (Finset.univ.filter fun a : ZMod p => a ≠ 0 ∧ N - a ≠ 0) =
      Finset.univ \ ({0, N} : Finset (ZMod p)) := by
    ext a
    simp [sub_ne_zero, eq_comm, and_comm]
  rw [heq, Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_univ, ZMod.card]
  by_cases h : N = 0
  · simp [h]
  · simp [h, Ne.symm h]

noncomputable section

open scoped Classical

/-- The numerical ordered pairs whose two summands are coprime,
indexed by the first coordinate. -/
def coprimePairFinset (allowed : ℕ → Prop) (n : ℕ) : Finset ℕ := by
  classical
  exact (pairFinset (PrimeCore allowed) (PrimeCore allowed) n).filter
    (fun a => a.gcd (n - a) = 1)

private theorem card_divisiblePairs {allowed : ℕ → Prop} {n d : ℕ}
    (hn : 0 < n) (hd : d ∣ n) (hcore : PrimeCore allowed d) :
    ((pairFinset (PrimeCore allowed) (PrimeCore allowed) n).filter
      (fun a => d ∣ a.gcd (n - a))).card =
    pairCount (PrimeCore allowed) (PrimeCore allowed) (n / d) := by
  classical
  have hdpos := Nat.pos_of_dvd_of_pos hd hn
  unfold pairCount
  refine (Finset.card_bij' (fun a _ => d * a) (fun a _ => a / d) ?_ ?_ ?_ ?_).symm
  · intro a ha
    obtain ⟨hal, ha₁, ha₂⟩ := mem_pairFinset.mp ha
    have heq : n - d * a = d * (n / d - a) := by
      rw [Nat.mul_sub_left_distrib, Nat.mul_div_cancel' hd]
    exact Finset.mem_filter.mpr ⟨mem_pairFinset.mpr
      ⟨by rw [← Nat.mul_div_cancel' hd]; exact Nat.mul_le_mul_left _ hal,
        primeCore_mul hcore ha₁, heq ▸ primeCore_mul hcore ha₂⟩,
      by rw [heq, Nat.gcd_mul_left]; exact dvd_mul_right _ _⟩
  · intro a ha
    obtain ⟨hab, hgcd⟩ := Finset.mem_filter.mp ha
    obtain ⟨hal, ha₁, ha₂⟩ := mem_pairFinset.mp hab
    have hda := hgcd.trans (Nat.gcd_dvd_left a (n - a))
    have hdb := hgcd.trans (Nat.gcd_dvd_right a (n - a))
    have heq : (n - a) / d = n / d - a / d := by
      apply Nat.eq_of_mul_eq_mul_left hdpos
      rw [Nat.mul_div_cancel' hdb, Nat.mul_sub_left_distrib,
        Nat.mul_div_cancel' hd, Nat.mul_div_cancel' hda]
    exact mem_pairFinset.mpr
      ⟨Nat.div_le_div_right hal, primeCore_of_dvd ha₁ (Nat.div_dvd_of_dvd hda),
        heq ▸ primeCore_of_dvd ha₂ (Nat.div_dvd_of_dvd hdb)⟩
  · intro a _
    exact Nat.mul_div_cancel_left a hdpos
  · intro a ha
    have hgcd := (Finset.mem_filter.mp ha).2
    exact Nat.mul_div_cancel' (hgcd.trans (Nat.gcd_dvd_left a (n - a)))

/-- Exact Möbius inversion for the coprime numerical pair count.
Specializing to C₁₂ gives equation (coprime-inversion). -/
theorem coprime_pair_moebius (allowed : ℕ → Prop) {n : ℕ} (hn : 0 < n) :
    ((coprimePairFinset allowed n).card : ℤ) =
      ∑ d ∈ n.divisors.filter (PrimeCore allowed),
        μ d * (pairCount (PrimeCore allowed) (PrimeCore allowed) (n / d) : ℤ) := by
  classical
  let pairs := pairFinset (PrimeCore allowed) (PrimeCore allowed) n
  let divisors := n.divisors.filter (PrimeCore allowed)
  change _ = ∑ d ∈ divisors, _
  symm
  calc
    (∑ d ∈ divisors, μ d * (pairCount (PrimeCore allowed) (PrimeCore allowed) (n / d) : ℤ)) =
        ∑ d ∈ divisors, ∑ a ∈ pairs, if d ∣ a.gcd (n - a) then μ d else 0 := by
      apply Finset.sum_congr rfl
      intro d hd
      obtain ⟨hdn, hcore⟩ := Finset.mem_filter.mp hd
      rw [← card_divisiblePairs hn (Nat.dvd_of_mem_divisors hdn) hcore, ← Finset.sum_filter]
      simp [pairs, mul_comm]
    _ = ∑ a ∈ pairs, ∑ d ∈ divisors, if d ∣ a.gcd (n - a) then μ d else 0 :=
      Finset.sum_comm
    _ = ∑ a ∈ pairs, if a.gcd (n - a) = 1 then (1 : ℤ) else 0 := by
      apply Finset.sum_congr rfl
      intro a ha
      obtain ⟨hal, ha₁, _⟩ := mem_pairFinset.mp ha
      let g := a.gcd (n - a)
      have hgpos : 0 < g := Nat.gcd_pos_of_pos_left _ ha₁.1
      have hga : g ∣ a := Nat.gcd_dvd_left ..
      have hgn : g ∣ n := by
        rw [← Nat.add_sub_of_le hal]
        exact Nat.dvd_add hga (Nat.gcd_dvd_right ..)
      have hfilter : divisors.filter (fun d => d ∣ g) = g.divisors := by
        ext d
        simp only [divisors, Finset.mem_filter, Nat.mem_divisors]
        exact ⟨fun h => ⟨h.2, hgpos.ne'⟩, fun h =>
          ⟨⟨⟨h.1.trans hgn, hn.ne'⟩, primeCore_of_dvd ha₁ (h.1.trans hga)⟩, h.1⟩⟩
      rw [← Finset.sum_filter, hfilter, ← ArithmeticFunction.coe_mul_zeta_apply,
        ArithmeticFunction.moebius_mul_coe_zeta, ArithmeticFunction.one_apply]
    _ = ((coprimePairFinset allowed n).card : ℤ) := by
      rw [← Finset.sum_filter]
      simp [coprimePairFinset, pairs]

end
end QuadraticNorms
