import QuadraticNorms.NormRoots
import Mathlib.NumberTheory.ArithmeticFunction.Defs
import Mathlib.Algebra.Field.ZMod

namespace QuadraticNorms

/-- The residues solving the polynomial of the binary norm form. -/
noncomputable def normRootResidues (b c : ℤ) (n : ℕ) [NeZero n] : Finset (ZMod n) := by
  classical
  exact Finset.univ.filter (fun t : ZMod n => t ^ 2 - (b : ZMod n) * t + (c : ZMod n) = 0)

theorem mem_normRootResidues {b c : ℤ} {n : ℕ} [NeZero n] {t : ZMod n} :
    t ∈ normRootResidues b c n ↔ t ^ 2 - (b : ZMod n) * t + (c : ZMod n) = 0 := by
  simp [normRootResidues]

/-- The polynomial root count, extended by zero at modulus zero. -/
noncomputable def normRootCount (b c : ℤ) : ArithmeticFunction ℕ :=
  ⟨fun n => if h : n = 0 then 0 else (letI : NeZero n := ⟨h⟩; (normRootResidues b c n).card),
    by simp⟩

theorem norm_root_count_apply (b c : ℤ) (n : ℕ) [NeZero n] :
    normRootCount b c n = (normRootResidues b c n).card := by
  change (if h : n = 0 then _ else _) = _
  rw [dif_neg (NeZero.ne n)]

theorem norm_root_count_of_no_root {b c : ℤ} {n : ℕ} (hn : 0 < n) (h : ¬ HasNormRoot b c n) :
    normRootCount b c n = 0 := by
  letI : NeZero n := ⟨hn.ne'⟩
  rw [norm_root_count_apply, Finset.card_eq_zero]
  apply Finset.eq_empty_iff_forall_notMem.mpr
  exact fun t ht => h ⟨t, mem_normRootResidues.mp ht⟩

@[simp] theorem norm_root_count_one (b c : ℤ) : normRootCount b c 1 = 1 := by
  rw [norm_root_count_apply]
  have h : normRootResidues b c 1 = Finset.univ := by
    ext t
    simp only [mem_normRootResidues, Finset.mem_univ, iff_true]
    exact Subsingleton.elim ..
  rw [h, Finset.card_univ, ZMod.card]

/-- The root count is multiplicative by the Chinese remainder theorem. -/
theorem norm_root_count_multiplicative (b c : ℤ) : (normRootCount b c).IsMultiplicative := by
  refine ⟨norm_root_count_one b c, fun {m n} hcoprime => ?_⟩
  by_cases hm : m = 0
  · simp [hm]
  by_cases hn : n = 0
  · simp [hn]
  letI : NeZero m := ⟨hm⟩
  letI : NeZero n := ⟨hn⟩
  letI : NeZero (m * n) := ⟨mul_ne_zero hm hn⟩
  let f := ZMod.chineseRemainder hcoprime
  have hpair (t : ZMod (m * n)) : t ∈ normRootResidues b c (m * n) ↔
      f t ∈ normRootResidues b c m ×ˢ normRootResidues b c n := by
    rw [mem_normRootResidues, Finset.mem_product, mem_normRootResidues, mem_normRootResidues]
    constructor
    · intro ht
      have h := congrArg f ht
      simp only [map_add, map_sub, map_pow, map_mul, map_intCast, map_zero] at h
      exact ⟨congrArg Prod.fst h, congrArg Prod.snd h⟩
    · intro h
      apply f.injective
      simp only [map_add, map_sub, map_pow, map_mul, map_intCast, map_zero]
      exact Prod.ext h.1 h.2
  rw [norm_root_count_apply, norm_root_count_apply, norm_root_count_apply, ← Finset.card_product]
  refine Finset.card_bij' (fun t _ => f t) (fun x _ => f.symm x)
    (fun t ht => (hpair t).mp ht) (fun x hx => (hpair (f.symm x)).mpr (by simpa only [RingEquiv.apply_symm_apply]))
    (fun t _ => f.symm_apply_apply t) (fun x _ => f.apply_symm_apply x)

private theorem norm_root_prime_card_le_two (b c : ℤ) {p : ℕ} (hp : p.Prime) :
    letI : NeZero p := ⟨hp.ne_zero⟩
    (normRootResidues b c p).card ≤ 2 := by
  letI : Fact p.Prime := ⟨hp⟩
  letI : NeZero p := ⟨hp.ne_zero⟩
  by_cases hne : (normRootResidues b c p).Nonempty
  · obtain ⟨a, ha⟩ := hne
    have hsub : normRootResidues b c p ⊆ {a, (b : ZMod p) - a} := by
      intro u hu
      have ha0 := mem_normRootResidues.mp ha
      have hu0 := mem_normRootResidues.mp hu
      have hprod : (u - a) * (u + a - (b : ZMod p)) = 0 := by linear_combination hu0 - ha0
      rcases mul_eq_zero.mp hprod with h | h
      · exact Finset.mem_insert.mpr (Or.inl (sub_eq_zero.mp h))
      · exact Finset.mem_insert_of_mem (Finset.mem_singleton.mpr (eq_sub_of_add_eq (sub_eq_zero.mp h)))
    exact (Finset.card_le_card hsub).trans (by simpa using Finset.card_insert_le a {(b : ZMod p) - a})
  · rw [Finset.not_nonempty_iff_eq_empty.mp hne]
    simp

/-- At a prime outside the discriminant, reduction injects the
prime-power roots into the at most two roots modulo p. -/
theorem norm_root_count_unramified_le_two (b c : ℤ) {p k : ℕ} (hp : p.Prime)
    (hne : ¬ (p : ℤ) ∣ 4 * c - b ^ 2) (hk : 0 < k) : normRootCount b c (p ^ k) ≤ 2 := by
  letI : Fact p.Prime := ⟨hp⟩
  letI : NeZero (p ^ k) := ⟨pow_ne_zero _ hp.ne_zero⟩
  let f := ZMod.castHom (dvd_pow_self p hk.ne') (ZMod p)
  have hmap {a : ZMod (p ^ k)} (ha : a ∈ normRootResidues b c (p ^ k)) :
      f a ∈ normRootResidues b c p := by
    apply mem_normRootResidues.mpr
    simpa only [map_add, map_sub, map_pow, map_mul, map_intCast, map_zero] using
      congrArg f (mem_normRootResidues.mp ha)
  have hInj : Set.InjOn f ↑(normRootResidues b c (p ^ k)) := by
    intro a ha u hu hau
    have ha' := mem_normRootResidues.mp ha
    have hu' := mem_normRootResidues.mp hu
    have hfa := mem_normRootResidues.mp (hmap ha)
    let v := a + u - (b : ZMod (p ^ k))
    have hv : f v ≠ 0 := by
      intro hzero
      have hdouble : 2 * f a - (b : ZMod p) = 0 := by
        simpa only [v, map_add, map_sub, map_intCast, ← hau, two_mul] using hzero
      have hdisc : ((4 * c - b ^ 2 : ℤ) : ZMod p) = 0 := by
        push_cast
        linear_combination 4 * hfa - (2 * f a - (b : ZMod p)) * hdouble
      exact hne ((ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp hdisc)
    have hvunit : IsUnit v := by
      rw [← ZMod.natCast_zmod_val v]
      apply (ZMod.isUnit_natCast_iff_not_dvd_pow hp hk).mpr
      intro hd
      apply hv
      rw [← ZMod.natCast_zmod_val v, map_natCast]
      exact (ZMod.natCast_eq_zero_iff _ _).mpr hd
    have hmul : (a - u) * v = 0 := by dsimp [v]; linear_combination ha' - hu'
    exact sub_eq_zero.mp (hvunit.mul_right_cancel (hmul.trans (zero_mul v).symm))
  rw [norm_root_count_apply]
  exact (Finset.card_le_card_of_injOn f (fun _ h => hmap h) hInj).trans
    (norm_root_prime_card_le_two b c hp)

end QuadraticNorms
