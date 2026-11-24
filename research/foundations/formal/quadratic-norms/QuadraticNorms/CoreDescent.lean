import QuadraticNorms.PrimitiveClassification
import QuadraticNorms.PackedCertificates
import Mathlib.Tactic.NormNum.Prime

namespace QuadraticNorms

/-- A target with no prime divisor in the split-prime core. -/
def IsReduced (C : ℕ → Prop) (n : ℕ) : Prop :=
  ∀ p : ℕ, p.Prime → p ∣ n → ¬ C p

private theorem inSumset_dilate_one {P : ℕ → Prop} {n : ℕ} :
    InSumset P (Dilate 1 P) n ↔ InSumset P P n := by
  simp [InSumset, Dilate]

private theorem coreSum_of_prime_descent' {allowed eligible : ℕ → Prop} {k : ℕ}
    (hreduced : ∀ n, eligible n → IsReduced (PrimeCore allowed) n →
      InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) n)
    (hstep : ∀ p n, p.Prime → PrimeCore allowed p → p ∣ n → eligible n →
      eligible (n / p) ∨ InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) n)
    {n : ℕ} (hnpos : 0 < n) (hn : eligible n) :
    InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) n := by
  apply coreSum_of_prime_descent _ _ hnpos hn
  · intro m hm h
    exact hreduced m hm (fun p hp hd hc => h p hp hd ((primeCore_prime hp).mp hc))
  · exact fun p m hp ha hd hm => hstep p m hp ((primeCore_prime hp).mpr ha) hd hm

/-- Eisenstein descent works for the ordinary and mixed target classes. -/
theorem C3_core_coverage {k r : ℕ}
    (hcover : ∀ n, 0 < n → n % 6 = r → IsReduced C3 n → InSumset C3 (Dilate k C3) n)
    {n : ℕ} (hn : 0 < n) (hclass : n % 6 = r) : InSumset C3 (Dilate k C3) n := by
  apply coreSum_of_prime_descent'
    (eligible := fun n => 0 < n ∧ n % 6 = r) (fun n h hred => hcover n h.1 h.2 hred)
    _ hn ⟨hn, hclass⟩
  intro p m hp hcore hpdiv ⟨hm, hclass⟩
  have hquot : 0 < m / p := Nat.div_pos (Nat.le_of_dvd hm hpdiv) hp.pos
  have hpone := C3_mod_six hcore
  have heq : p * (m / p) = m := Nat.mul_div_cancel' hpdiv
  have hmod : p * (m / p) ≡ m / p [MOD 6] := by
    simpa using (show p ≡ 1 [MOD 6] from hpone).mul (Nat.ModEq.refl (m / p))
  rw [heq] at hmod
  exact Or.inl ⟨hquot, (show m % 6 = (m / p) % 6 from hmod).symm.trans hclass⟩

/-- Dividing an even target by an odd prime factor keeps the
quotient even and at least two. A forbidden divisor of the quotient
would divide the original target. -/
private theorem eligible_even_quotient {p n m : ℕ} (hp : p.Prime) (hodd : p % 2 = 1)
    (hpd : p ∣ n) (hn : 2 ≤ n) (heven : n % 2 = 0) (havoid : ¬ m ∣ n) :
    2 ≤ n / p ∧ (n / p) % 2 = 0 ∧ ¬ m ∣ n / p := by
  have hq : 0 < n / p := Nat.div_pos (Nat.le_of_dvd (by omega) hpd) hp.pos
  have heq : p * (n / p) = n := Nat.mul_div_cancel' hpd
  have hmod : p * (n / p) ≡ n / p [MOD 2] := by
    simpa using (show p ≡ 1 [MOD 2] from hodd).mul (Nat.ModEq.refl (n / p))
  rw [heq] at hmod
  have hqeven : (n / p) % 2 = 0 := (show n % 2 = (n / p) % 2 from hmod).symm.trans heven
  exact ⟨by omega, hqeven, fun h => havoid (h.trans (Nat.div_dvd_of_dvd hpd))⟩

theorem C8_ordinary_core_coverage
    (hcover : ∀ n, 2 ≤ n → n % 2 = 0 → ¬ 8 ∣ n → IsReduced C8 n → InSumset C8 C8 n)
    {n : ℕ} (hn : 2 ≤ n) (heven : n % 2 = 0) (h₈ : ¬ 8 ∣ n) : InSumset C8 C8 n := by
  apply inSumset_dilate_one.mp
  apply coreSum_of_prime_descent' (eligible := fun n => 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 8 ∣ n)
    _ _ (by omega) ⟨hn, heven, h₈⟩
  · exact fun m ⟨h₂, he, h₈⟩ hred => inSumset_dilate_one.mpr (hcover m h₂ he h₈ hred)
  · intro p m hp hcore hd ⟨h₂, he, h₈⟩
    have := C8_mod_eight hcore
    exact Or.inl (eligible_even_quotient hp (by omega) hd h₂ he h₈)

theorem C8_mixed_core_coverage
    (hcover : ∀ n, 3 ≤ n → n % 2 = 1 → IsReduced C8 n → InSumset C8 (Dilate 2 C8) n)
    (hprime : ∀ p, p.Prime → C8 p → InSumset C8 (Dilate 2 C8) p)
    {n : ℕ} (hn : 3 ≤ n) (hodd : n % 2 = 1) : InSumset C8 (Dilate 2 C8) n := by
  apply coreSum_of_prime_descent' (eligible := fun n => 3 ≤ n ∧ n % 2 = 1)
    (fun m ⟨hm, ho⟩ hr => hcover m hm ho hr) _ (by omega) ⟨hn, hodd⟩
  intro p m hp hcore hd ⟨hm, hodd⟩
  have hpmod := C8_mod_eight hcore
  have hpodd : p % 2 = 1 := by omega
  have heq : p * (m / p) = m := Nat.mul_div_cancel' hd
  have hmod : p * (m / p) ≡ m / p [MOD 2] := by
    simpa using (show p ≡ 1 [MOD 2] from hpodd).mul (Nat.ModEq.refl (m / p))
  rw [heq] at hmod
  have hqodd : (m / p) % 2 = 1 := (show m % 2 = (m / p) % 2 from hmod).symm.trans hodd
  by_cases hq : m / p = 1
  · rw [hq, mul_one] at heq
    exact Or.inr (heq ▸ hprime p hp hcore)
  · have hqge : 3 ≤ m / p := by
      have hqrem := Nat.mod_add_div (m / p) 2
      omega
    exact Or.inl ⟨hqge, hqodd⟩

theorem C11_core_coverage
    (hcover : ∀ n, 2 ≤ n → n % 2 = 0 → ¬ 11 ∣ n → IsReduced C11 n → InSumset C11 C11 n)
    {n : ℕ} (hn : 2 ≤ n) (heven : n % 2 = 0) (h₁₁ : ¬ 11 ∣ n) : InSumset C11 C11 n := by
  apply inSumset_dilate_one.mp
  apply coreSum_of_prime_descent' (eligible := fun n => 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 11 ∣ n)
    _ _ (by omega) ⟨hn, heven, h₁₁⟩
  · exact fun m ⟨h₂, he, h₁₁⟩ hred => inSumset_dilate_one.mpr (hcover m h₂ he h₁₁ hred)
  · intro p m hp hcore hd ⟨h₂, he, h₁₁⟩
    have hpne : p ≠ 2 := by
      intro hp₂
      subst p
      have h := hcore.2 2 Nat.prime_two dvd_rfl
      contradiction
    exact Or.inl (eligible_even_quotient hp (hp.eq_two_or_odd.resolve_left hpne) hd h₂ he h₁₁)

private theorem C7_scale {g n : ℕ} (hg : C7 g) (h : InSumset C7 C7 n) :
    InSumset C7 C7 (g * n) :=
  h.scale (fun _ ha => primeCore_mul hg ha) (fun _ hb => primeCore_mul hg hb)

private theorem C7_product_patches :
    InSumset C7 C7 (11 * 11) ∧ InSumset C7 C7 (11 * 29) ∧ InSumset C7 C7 (29 * 29) := by
  have hc₂ : C7 2 := (primeCore_prime (by decide : Nat.Prime 2)).mpr (by decide)
  have hc₄ := primeCore_mul hc₂ hc₂
  have hc₈ := primeCore_mul hc₂ hc₄
  have hc₃₂ := primeCore_mul hc₄ hc₈
  have hc₁₁₃ : C7 113 := (primeCore_prime (by decide : Nat.Prime 113)).mpr (by decide)
  have hc₃₁₇ : C7 317 := (primeCore_prime (by norm_num : Nat.Prime 317)).mpr (by decide)
  have hc₈₀₉ : C7 809 := (primeCore_prime (by norm_num : Nat.Prime 809)).mpr (by decide)
  exact ⟨⟨8, 113, hc₈, hc₁₁₃, by decide⟩, ⟨2, 317, hc₂, hc₃₁₇, by decide⟩,
    ⟨32, 809, hc₃₂, hc₈₀₉, by decide⟩⟩

/-- The −7 exceptional quotients are split primes themselves. If
neither quotient can be covered by descending, their three products
are covered by the explicit pairs checked above. -/
theorem C7_core_coverage
    (hcover : ∀ n, 2 ≤ n → ¬ 7 ∣ n → n ≠ 11 → n ≠ 29 → IsReduced C7 n → InSumset C7 C7 n)
    (hprime : ∀ p, p.Prime → C7 p → p ≠ 11 → p ≠ 29 → InSumset C7 C7 p)
    {n : ℕ} (hn : 2 ≤ n) (h₇ : ¬ 7 ∣ n) (h₁₁ : n ≠ 11) (h₂₉ : n ≠ 29) :
    InSumset C7 C7 n := by
  induction n using Nat.strong_induction_on with
  | h n ih =>
    by_cases hred : IsReduced C7 n
    · exact hcover n hn h₇ h₁₁ h₂₉ hred
    unfold IsReduced at hred
    push Not at hred
    obtain ⟨p, hp, hd, hcore⟩ := hred
    obtain ⟨m, rfl⟩ := hd
    have hmp : 0 < m := by
      by_contra! h
      simp [show m = 0 by omega] at hn
    by_cases hm : m = 1
    · subst m
      simp only [mul_one] at *
      exact hprime p hp hcore h₁₁ h₂₉
    · have hm₂ : 2 ≤ m := by omega
      have hltp : p < p * m := by
        simpa using Nat.mul_lt_mul_of_pos_left (show 1 < m by omega) hp.pos
      have hmlt : m < p * m := by
        simpa using Nat.mul_lt_mul_of_pos_right hp.one_lt hmp
      have hm₇ : ¬ 7 ∣ m := fun h => h₇ (h.trans (dvd_mul_left ..))
      by_cases hsmall : m = 11 ∨ m = 29
      · have hc : C7 m := hsmall.elim
          (fun e => e ▸ (primeCore_prime (by decide : Nat.Prime 11)).mpr (by decide))
          (fun e => e ▸ (primeCore_prime (by decide : Nat.Prime 29)).mpr (by decide))
        by_cases hsp : p = 11 ∨ p = 29
        · rcases hsp with rfl | rfl <;> rcases hsmall with rfl | rfl
          · exact C7_product_patches.1
          · exact C7_product_patches.2.1
          · simpa [Nat.mul_comm] using C7_product_patches.2.1
          · exact C7_product_patches.2.2
        · rw [Nat.mul_comm]
          exact C7_scale hc (ih p hltp hp.two_le
            (fun h => h₇ (h.trans (dvd_mul_right ..))) (fun h => hsp (Or.inl h))
            (fun h => hsp (Or.inr h)))
      · exact C7_scale hcore (ih m hmlt hm₂ hm₇
          (fun h => hsmall (Or.inl h)) (fun h => hsmall (Or.inr h)))

/-- Coverage statements required on reduced targets and the two
split-prime boundaries in the four-discriminant result. -/
structure ReducedCoreCoverage : Prop where
  eisensteinOrdinary : ∀ n, 0 < n → n % 6 = 2 → IsReduced C3 n → InSumset C3 C3 n
  eisensteinMixed : ∀ n, 0 < n → n % 6 = 4 → IsReduced C3 n → InSumset C3 (Dilate 3 C3) n
  seven : ∀ n, 2 ≤ n → ¬ 7 ∣ n → n ≠ 11 → n ≠ 29 → IsReduced C7 n → InSumset C7 C7 n
  sevenPrimes : ∀ p, p.Prime → C7 p → p ≠ 11 → p ≠ 29 → InSumset C7 C7 p
  eightOrdinary : ∀ n, 2 ≤ n → n % 2 = 0 → ¬ 8 ∣ n → IsReduced C8 n → InSumset C8 C8 n
  eightMixed : ∀ n, 3 ≤ n → n % 2 = 1 → IsReduced C8 n → InSumset C8 (Dilate 2 C8) n
  eightPrimes : ∀ p, p.Prime → C8 p → InSumset C8 (Dilate 2 C8) p
  eleven : ∀ n, 2 ≤ n → n % 2 = 0 → ¬ 11 ∣ n → IsReduced C11 n → InSumset C11 C11 n

/-- The four rows of the exact overview, reduced to the finite and
infinite coverage obligations on the split-prime cores. -/
theorem exact_overview_of_reduced_coverage (h : ReducedCoreCoverage) (n : ℕ) :
    (InSumset P3 P3 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ n % 18 ≠ 0 ∧ n % 18 ≠ 12) ∧
    (InSumset P7 P7 n ↔ 2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203) ∧
    (InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n) ∧
    (InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n) := by
  refine ⟨P3_sum_classification restoration3 ?_ ?_ n,
    P7_sum_classification restoration7 (fun _ h₂ h₇ h₁₁ h₂₉ =>
      C7_core_coverage h.seven h.sevenPrimes h₂ h₇ h₁₁ h₂₉) n,
    P8_sum_classification restoration8 (fun _ => C8_ordinary_core_coverage h.eightOrdinary)
      (fun _ => C8_mixed_core_coverage h.eightMixed h.eightPrimes) n,
    P11_sum_classification restoration11 (fun _ => C11_core_coverage h.eleven) n⟩
  · intro m hm
    exact inSumset_dilate_one.mp
      (C3_core_coverage (fun a ha hr hred => inSumset_dilate_one.mpr
        (h.eisensteinOrdinary a ha hr hred)) (by omega) hm)
  · intro m hm
    exact C3_core_coverage h.eisensteinMixed (by omega) hm

private theorem C3_ordinary_necessary {n : ℕ} (h : InSumset C3 C3 n) : n % 6 = 2 := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have := C3_mod_six ha
  have := C3_mod_six hb
  omega

private theorem C3_mixed_necessary {n : ℕ} (h : InSumset C3 (Dilate 3 C3) n) : n % 6 = 4 := by
  obtain ⟨a, _, ha, ⟨b, hb, rfl⟩, rfl⟩ := h
  have := C3_mod_six ha
  have := C3_mod_six hb
  omega

private theorem C7_ordinary_local_necessary {n : ℕ} (h : InSumset C7 C7 n) :
    2 ≤ n ∧ ¬ 7 ∣ n := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have hapos := ha.1
  have hbpos := hb.1
  refine ⟨by omega, fun hd => ?_⟩
  have ⟨ha₇, _⟩ := P7_ramified_pair (restoration7.core a ha) (restoration7.core b hb) hd
  have hbad := ha.2 7 (by decide) ha₇
  contradiction

private def finiteC7 (n : ℕ) : Bool :=
  primeCoreTest (fun p => decide (p % 7 = 1 ∨ p % 7 = 2 ∨ p % 7 = 4)) n

private theorem finiteC7_spec {n : ℕ} : finiteC7 n = true ↔ C7 n := by
  simpa [finiteC7, C7] using
    (primeCoreTest_spec (allowed := fun p => decide (p % 7 = 1 ∨ p % 7 = 2 ∨ p % 7 = 4)) (n := n))

private theorem C7_ordinary_excludes_small : ¬ InSumset C7 C7 11 ∧ ¬ InSumset C7 C7 29 := by
  have excludes {n : ℕ}
      (htest : ((List.range (n + 1)).any (fun a => finiteC7 a && finiteC7 (n - a))) = false) :
      ¬ InSumset C7 C7 n := by
    rintro ⟨a, b, ha, hb, rfl⟩
    have htrue : ((List.range (a + b + 1)).any (fun t => finiteC7 t && finiteC7 (a + b - t))) = true := by
      apply List.any_eq_true.mpr
      refine ⟨a, List.mem_range.mpr (by omega), ?_⟩
      rw [Bool.and_eq_true]
      refine ⟨finiteC7_spec.mpr ha, ?_⟩
      simpa using finiteC7_spec.mpr hb
    rw [htest] at htrue
    contradiction
  exact ⟨excludes (by decide +kernel), excludes (by decide +kernel)⟩

theorem C8_ordinary_necessary {n : ℕ} (h : InSumset C8 C8 n) :
    2 ≤ n ∧ n % 2 = 0 ∧ ¬ 8 ∣ n := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have := ha.1
  have := hb.1
  have := C8_mod_eight ha
  have := C8_mod_eight hb
  omega

theorem C8_mixed_necessary {n : ℕ} (h : InSumset C8 (Dilate 2 C8) n) :
    3 ≤ n ∧ n % 2 = 1 := by
  obtain ⟨a, _, ha, ⟨b, hb, rfl⟩, rfl⟩ := h
  have := ha.1
  have := hb.1
  have := C8_mod_eight ha
  omega

private theorem C11_ordinary_necessary {n : ℕ} (h : InSumset C11 C11 n) :
    2 ≤ n ∧ n % 2 = 0 ∧ ¬ 11 ∣ n := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have := ha.1
  have := hb.1
  have hpa := restoration11.core a ha
  have hpb := restoration11.core b hb
  have := P11_odd hpa
  have := P11_odd hpb
  refine ⟨by omega, by omega, fun hd => ?_⟩
  have ⟨ha₁₁, _⟩ := P11_ramified_pair hpa hpb hd
  have hbad := ha.2 11 (by decide) ha₁₁
  contradiction

/-- All core sumset identities follow from coverage on the reduced
targets and the split primes that cannot be reached by descent. -/
theorem core_sumsets_of_reduced_coverage (h : ReducedCoreCoverage) (n : ℕ) :
    (InSumset C3 C3 n ↔ n % 6 = 2) ∧
    (InSumset C3 (Dilate 3 C3) n ↔ n % 6 = 4) ∧
    (InSumset C7 C7 n ↔ 2 ≤ n ∧ ¬ 7 ∣ n ∧ n ≠ 11 ∧ n ≠ 29) ∧
    (InSumset C8 C8 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 8 ∣ n) ∧
    (InSumset C8 (Dilate 2 C8) n ↔ 3 ≤ n ∧ n % 2 = 1) ∧
    (InSumset C11 C11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 11 ∣ n) := by
  refine ⟨⟨C3_ordinary_necessary, ?_⟩, ⟨C3_mixed_necessary, ?_⟩, ⟨?_, ?_⟩,
    ⟨C8_ordinary_necessary, ?_⟩, ⟨C8_mixed_necessary, ?_⟩, ⟨C11_ordinary_necessary, ?_⟩⟩
  · intro hn
    exact inSumset_dilate_one.mp (C3_core_coverage
      (fun m hm hr hred => inSumset_dilate_one.mpr (h.eisensteinOrdinary m hm hr hred)) (by omega) hn)
  · intro hn
    exact C3_core_coverage h.eisensteinMixed (by omega) hn
  · intro hn
    obtain ⟨h₂, h₇⟩ := C7_ordinary_local_necessary hn
    exact ⟨h₂, h₇, fun he => C7_ordinary_excludes_small.1 (he ▸ hn),
      fun he => C7_ordinary_excludes_small.2 (he ▸ hn)⟩
  · exact fun ⟨h₂, h₇, h₁₁, h₂₉⟩ => C7_core_coverage h.seven h.sevenPrimes h₂ h₇ h₁₁ h₂₉
  · exact fun ⟨h₂, he, h₈⟩ => C8_ordinary_core_coverage h.eightOrdinary h₂ he h₈
  · exact fun ⟨h₃, ho⟩ => C8_mixed_core_coverage h.eightMixed h.eightPrimes h₃ ho
  · exact fun ⟨h₂, he, h₁₁⟩ => C11_core_coverage h.eleven h₂ he h₁₁

end QuadraticNorms
