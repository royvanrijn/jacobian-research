import QuadraticNorms.PackedCertificates
import QuadraticNorms.PrimeCores
import Mathlib.Tactic.Ring

namespace QuadraticNorms

/-- If an integer 1 modulo six has a forbidden divisor, then one
such divisor is at most its square root. -/
theorem C3_of_small_prime_test {n bound : ℕ} (hn : n % 6 = 1) (hbound : n < bound ^ 2)
    (h : ∀ p < bound, p.Prime → p % 6 = 5 → ¬ p ∣ n) : C3 n := by
  refine ⟨by omega, ?_⟩
  intro p hp hpn
  by_contra hgood
  have hpmod : p % 6 = 5 := by
    have hcop : p ≠ 2 := by intro he; subst p; omega
    have hpodd := hp.eq_two_or_odd.resolve_left hcop
    have hp₃ : p % 3 ≠ 0 := by
      intro he
      have hp₃ := Nat.dvd_of_mod_eq_zero he
      rcases (Nat.dvd_prime hp).mp hp₃ with h | rfl <;> omega
    have h₂ := Nat.mod_mod_of_dvd p (by decide : 2 ∣ 6)
    have h₃ := Nat.mod_mod_of_dvd p (by decide : 3 ∣ 6)
    rw [hpodd] at h₂
    omega
  have hpbound : bound ≤ p := by
    by_contra! hlt
    exact h p hlt hp hpmod hpn
  obtain ⟨m, rfl⟩ := hpn
  have hmlt : m < bound := by
    by_contra! hle
    have := Nat.mul_le_mul hpbound hle
    nlinarith
  have hmpos : 0 < m := by
    by_contra! hm
    simp [show m = 0 by omega] at hn
  have hm_mod : m % 6 = 5 := by
    have := (show p ≡ 5 [MOD 6] from hpmod).mul (Nat.ModEq.refl m)
    change p * m % 6 = (5 * m) % 6 at this
    omega
  have hmcore : ¬ C3 m := fun hm => by have := C3_mod_six hm; omega
  unfold C3 PrimeCore at hmcore
  push Not at hmcore
  obtain ⟨q, hq, hqm, hqbad⟩ := hmcore hmpos
  have hqn := hqm.trans (dvd_mul_left m p)
  have hqlt := (Nat.le_of_dvd hmpos hqm).trans_lt hmlt
  have hqmod : q % 6 = 5 := by
    have hqne : q ≠ 2 := by
      intro he
      subst q
      have := Nat.dvd_trans hqm (dvd_mul_left m p)
      omega
    have hqodd := hq.eq_two_or_odd.resolve_left hqne
    have hq₃ : q % 3 ≠ 0 := by
      intro he
      have hq₃ := Nat.dvd_of_mod_eq_zero he
      rcases (Nat.dvd_prime hq).mp hq₃ with h | rfl <;> omega
    have h₂ := Nat.mod_mod_of_dvd q (by decide : 2 ∣ 6)
    have h₃ := Nat.mod_mod_of_dvd q (by decide : 3 ∣ 6)
    rw [hqodd] at h₂
    omega
  exact h q hqlt hq hqmod hqn

private theorem inertIndex_iff {p k : ℕ} (hp : p % 6 = 5) :
    p ∣ 6 * k + 1 ↔ ∃ j, k = (5 * p - 1) / 6 + p * j := by
  have hstart : 6 * ((5 * p - 1) / 6) + 1 = 5 * p := by omega
  constructor
  · rintro ⟨a, ha⟩
    have hm := (show p ≡ 5 [MOD 6] from hp).mul (Nat.ModEq.refl a)
    rw [← ha] at hm
    change (6 * k + 1) % 6 = (5 * a) % 6 at hm
    have hamod : a % 6 = 5 := by omega
    have hquot : a = 5 + 6 * (a / 6) := by omega
    refine ⟨a / 6, ?_⟩
    conv at ha => rhs; rw [hquot, Nat.mul_add, ← Nat.mul_left_comm]
    omega
  · rintro ⟨j, rfl⟩
    exact ⟨5 + 6 * j, by
      calc
        6 * ((5 * p - 1) / 6 + p * j) + 1 =
            (6 * ((5 * p - 1) / 6) + 1) + 6 * p * j := by ring
        _ = p * (5 + 6 * j) := by rw [hstart]; ring⟩

private def forbiddenPrimeMask (H p : ℕ) : ℕ :=
  progressionMask p (H / p + 1) <<< ((5 * p - 1) / 6)

private theorem forbiddenPrimeMask_bit {H p k : ℕ} (hp : p % 6 = 5) (hk : k ≤ H) :
    (forbiddenPrimeMask H p).testBit k = true ↔ p ∣ 6 * k + 1 := by
  rw [forbiddenPrimeMask, Nat.testBit_shiftLeft, Bool.and_eq_true, decide_eq_true_eq,
    progressionMask_bit, inertIndex_iff hp]
  constructor
  · rintro ⟨hk', j, _, heq⟩
    exact ⟨j, by omega⟩
  · rintro ⟨j, heq⟩
    have hp' : 0 < p := by omega
    have hj : j ≤ H / p :=
      (Nat.le_div_iff_mul_le hp').mpr (by rw [Nat.mul_comm]; omega)
    exact ⟨by omega, j, Nat.lt_succ_of_le hj, by omega⟩

private def inertPrimesBelow (bound : ℕ) : List ℕ :=
  (List.range bound).filter (fun p => decide (p.Prime ∧ p % 6 = 5))

private theorem mem_inertPrimesBelow {bound p : ℕ} :
    p ∈ inertPrimesBelow bound ↔ p < bound ∧ p.Prime ∧ p % 6 = 5 := by
  simp [inertPrimesBelow]

private def badIndexMask (H : ℕ) (primes : List ℕ) : ℕ :=
  primes.foldr (fun p acc => forbiddenPrimeMask H p ||| acc) 0

private theorem badIndexMask_bit {H k : ℕ} (hk : k ≤ H) (primes : List ℕ)
    (hp : ∀ p ∈ primes, p % 6 = 5) :
    (badIndexMask H primes).testBit k = true ↔ ∃ p ∈ primes, p ∣ 6 * k + 1 := by
  induction primes with
  | nil => simp [badIndexMask]
  | cons p primes ih =>
    change (forbiddenPrimeMask H p ||| badIndexMask H primes).testBit k = true ↔ _
    rw [Nat.testBit_or, Bool.or_eq_true, forbiddenPrimeMask_bit (hp _ (by simp)) hk,
      ih (fun q hq => hp q (List.mem_cons_of_mem _ hq))]
    simp only [List.mem_cons, exists_eq_or_imp]

/-- Bit j of the Eisenstein membership mask encodes the value 6j+1.
Only genuine inert primes below the chosen square-root bound are sieved. -/
def eisensteinCoreMask (H bound : ℕ) : ℕ :=
  let all := 2 ^ (H + 1) - 1
  all &&& (all ^^^ badIndexMask H (inertPrimesBelow bound))

theorem eisensteinCoreMask_bit {H bound j : ℕ} (hj : j ≤ H) (hbound : 6 * j + 1 < bound ^ 2) :
    (eisensteinCoreMask H bound).testBit j = true ↔ C3 (6 * j + 1) := by
  unfold eisensteinCoreMask
  rw [Nat.testBit_and, Nat.testBit_xor, Nat.testBit_two_pow_sub_one, decide_eq_true (by omega),
    Bool.true_and, Bool.true_xor]
  have hnot (b : Bool) : (!b) = true ↔ b ≠ true := by cases b <;> decide
  rw [hnot]
  have hbad := badIndexMask_bit hj (inertPrimesBelow bound)
    (fun p hp => (mem_inertPrimesBelow.mp hp).2.2)
  constructor
  · intro h
    apply C3_of_small_prime_test (bound := bound) (by omega) hbound
    intro p hl hp hmod hd
    exact h (hbad.mpr ⟨p, mem_inertPrimesBelow.mpr ⟨hl, hp, hmod⟩, hd⟩)
  · intro h hset
    obtain ⟨p, hp, hpd⟩ := hbad.mp hset
    have ⟨_, hprime, hmod⟩ := mem_inertPrimesBelow.mp hp
    have := h.2 p hprime hpd
    omega

/-- The anchor index is i. In the ordinary scan the partner index
must be at least i, while the mixed scan also permits partner zero. -/
def eisensteinWitnessMask (scale : ℕ) (lower : ℕ → ℕ) (anchors : List ℕ) (values : ℕ) : ℕ :=
  anchors.foldr (fun i acc =>
    if values.testBit i then ((values >>> lower i) <<< (scale * i + lower i)) ||| acc else acc) 0

theorem eisensteinWitnessMask_bit {scale : ℕ} {lower : ℕ → ℕ} {anchors : List ℕ} {values k : ℕ} :
    (eisensteinWitnessMask scale lower anchors values).testBit k = true ↔
    ∃ i ∈ anchors, values.testBit i = true ∧ scale * i + lower i ≤ k ∧
      values.testBit (k - scale * i) = true := by
  induction anchors with
  | nil => simp [eisensteinWitnessMask]
  | cons i anchors ih =>
    change (if values.testBit i then
      ((values >>> lower i) <<< (scale * i + lower i)) |||
        eisensteinWitnessMask scale lower anchors values
      else eisensteinWitnessMask scale lower anchors values).testBit k = true ↔ _
    by_cases h : values.testBit i = true
    · rw [if_pos h, Nat.testBit_or, Nat.testBit_shiftLeft, Nat.testBit_shiftRight,
        Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq, ih]
      simp only [List.mem_cons, exists_eq_or_imp, h, true_and]
      have heq : ∀ hk : scale * i + lower i ≤ k,
          lower i + (k - (scale * i + lower i)) = k - scale * i := fun hk => by omega
      aesop
    · rw [if_neg h, ih]
      simp [List.mem_cons, exists_eq_or_imp, h]

/-- Pack the target interval into a single bitwise inclusion test. -/
def eisensteinScanCheck (scale H short : ℕ) (lower : ℕ → ℕ) (values : ℕ) : Bool :=
  let requested := 2 ^ (H + 1) - 1
  (requested &&& eisensteinWitnessMask scale lower (List.range short) values) == requested

theorem eisensteinScanCheck_sound {scale H short values k : ℕ} {lower : ℕ → ℕ}
    (h : eisensteinScanCheck scale H short lower values = true) (hk : k ≤ H) :
    ∃ i < short, values.testBit i = true ∧ scale * i + lower i ≤ k ∧
      values.testBit (k - scale * i) = true := by
  have hcheck : (2 ^ (H + 1) - 1) &&& eisensteinWitnessMask scale lower (List.range short) values =
      2 ^ (H + 1) - 1 := by simpa [eisensteinScanCheck] using h
  have hbit := congrArg (fun bits : ℕ => bits.testBit k) hcheck
  rw [Nat.testBit_and, Nat.testBit_two_pow_sub_one, decide_eq_true (by omega), Bool.true_and] at hbit
  simpa only [List.mem_range] using eisensteinWitnessMask_bit.mp hbit

private def eisensteinFullMask : ℕ := eisensteinCoreMask 16666666 10000

/-- The numerical test in Proposition finite-certificates, ordinary row. -/
def eisensteinOrdinaryCheck : Bool :=
  eisensteinScanCheck 1 16666666 36 id eisensteinFullMask

/-- The numerical test in Proposition finite-certificates, mixed row. -/
def eisensteinMixedCheck : Bool :=
  eisensteinScanCheck 3 16666665 34 (fun _ => 0) eisensteinFullMask

theorem eisensteinOrdinaryCertificate_sound (h : eisensteinOrdinaryCheck = true) {k : ℕ}
    (hk : 6 * k + 2 < 100000000) :
    ∃ i ≤ 35, 2 * i ≤ k ∧ C3 (6 * i + 1) ∧ C3 (6 * (k - i) + 1) := by
  obtain ⟨i, hi, hbit, hil, hkbit⟩ := eisensteinScanCheck_sound
    (scale := 1) (lower := id) (short := 36) (values := eisensteinFullMask)
    (by exact h) (by omega : k ≤ 16666666)
  simp only [one_mul, id_eq] at hil hkbit
  have hin : C3 (6 * i + 1) :=
    (eisensteinCoreMask_bit (by omega : i ≤ 16666666) (by omega)).mp hbit
  have hkn : C3 (6 * (k - i) + 1) :=
    (eisensteinCoreMask_bit (by omega : k - i ≤ 16666666) (by omega)).mp hkbit
  exact ⟨i, by omega, by omega, hin, hkn⟩

theorem eisensteinMixedCertificate_sound (h : eisensteinMixedCheck = true) {k : ℕ}
    (hk : 6 * k + 4 < 100000000) :
    ∃ j ≤ 33, 3 * j ≤ k ∧ C3 (6 * j + 1) ∧ C3 (6 * (k - 3 * j) + 1) := by
  obtain ⟨j, hj, hbit, hjl, hkbit⟩ := eisensteinScanCheck_sound
    (scale := 3) (lower := fun _ => 0) (short := 34) (values := eisensteinFullMask)
    (by exact h) (by omega : k ≤ 16666665)
  have hjn : C3 (6 * j + 1) :=
    (eisensteinCoreMask_bit (by omega : j ≤ 16666666) (by omega)).mp hbit
  have hkn : C3 (6 * (k - 3 * j) + 1) :=
    (eisensteinCoreMask_bit (by omega : k - 3 * j ≤ 16666666) (by omega)).mp hkbit
  exact ⟨j, by omega, by omega, hjn, hkn⟩

end QuadraticNorms
