import QuadraticNorms.PackedCertificates
import QuadraticNorms.PrimeCores

namespace QuadraticNorms

private theorem primeCore_of_small_test {allowed : ℕ → Prop} {n bound : ℕ}
    (hone : allowed 1) (hmul : ∀ a b, allowed a → allowed b → allowed (a * b))
    (hcancel : ∀ a b, allowed (a * b) → allowed b → allowed a)
    (hn : 0 < n) (hres : allowed n) (hbound : n < bound ^ 2)
    (htest : ∀ p < bound, p.Prime → ¬ allowed p → ¬ p ∣ n) :
    PrimeCore allowed n := by
  refine ⟨hn, ?_⟩
  intro p hp hpn
  by_contra! hbad
  have hpbound : bound ≤ p := by
    by_contra! hpbound
    exact htest p hpbound hp hbad hpn
  obtain ⟨m, rfl⟩ := hpn
  have hmp : 0 < m := Nat.pos_of_mul_pos_left hn
  have hmlt : m < bound := by
    by_contra! hle
    have hcontr := Nat.mul_le_mul hpbound hle
    rw [← pow_two] at hcontr
    omega
  apply hbad (hcancel p m hres ?_)
  apply primeCore_induction hone (fun _ _ h => h) hmul
  refine ⟨hmp, fun q hq hqm => ?_⟩
  by_contra! hqbad
  exact htest q ((Nat.le_of_dvd hmp hqm).trans_lt hmlt) hq hqbad (hqm.trans (dvd_mul_left ..))

private def periodicMask (period bound : ℕ) (residues : List ℕ) : ℕ :=
  packedShiftUnion residues (progressionMask period (bound / period + 1))

private theorem periodicMask_bit {period bound n : ℕ} {residues : List ℕ}
    (hn : n ≤ bound) (hres : ∀ r ∈ residues, r < period) :
    (periodicMask period bound residues).testBit n = true ↔ n % period ∈ residues := by
  rw [periodicMask, packedShiftUnion_bit]
  constructor
  · rintro ⟨r, hr, hrn, hbit⟩
    obtain ⟨j, _, hj⟩ := progressionMask_bit .. |>.mp hbit
    have heq : n = r + period * j := by omega
    rw [heq, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt (hres r hr)]
    exact hr
  · intro h
    have hmod := Nat.mod_add_div n period
    exact ⟨n % period, h, Nat.mod_le _ _, (progressionMask_bit ..).mpr
      ⟨n / period, Nat.lt_succ_of_le (Nat.div_le_div_right hn), by omega⟩⟩

private def badMultiplesMask (bound : ℕ) (primes : List ℕ) : ℕ :=
  primes.foldr (fun p acc => progressionMask p (bound / p + 1) ||| acc) 0

private theorem badMultiplesMask_bit {bound n : ℕ} (hn : n ≤ bound) (primes : List ℕ)
    (hpositive : ∀ p ∈ primes, 0 < p) :
    (badMultiplesMask bound primes).testBit n = true ↔ ∃ p ∈ primes, p ∣ n := by
  induction primes with
  | nil => simp [badMultiplesMask]
  | cons p primes ih =>
    change (progressionMask p (bound / p + 1) ||| badMultiplesMask bound primes).testBit n = true ↔ _
    rw [Nat.testBit_or, Bool.or_eq_true, progressionMask_bit, ih (fun q hq => hpositive q (by simp [hq]))]
    simp only [List.mem_cons, exists_eq_or_imp]
    have hequiv : (∃ j < bound / p + 1, n = p * j) ↔ p ∣ n := by
      constructor
      · rintro ⟨j, _, rfl⟩
        exact dvd_mul_right ..
      · rintro ⟨j, rfl⟩
        have hle := (Nat.le_div_iff_mul_le (hpositive p (by simp))).mpr (by simpa [Nat.mul_comm] using hn)
        exact ⟨j, Nat.lt_succ_of_le hle, rfl⟩
    rw [hequiv]

/-- Sieve values on prescribed residue classes by removing multiples
of the forbidden primes below the square-root bound. -/
def periodicCoreMask (period bound primeBound : ℕ) (residues : List ℕ) : ℕ :=
  let all := 2 ^ (bound + 1) - 1
  let excluded := (List.range primeBound).filter
    (fun p => decide (p.Prime ∧ p % period ∉ residues))
  let candidate := periodicMask period bound residues
  (all &&& candidate) &&& (all ^^^ badMultiplesMask bound excluded)

theorem periodicCoreMask_bit {period bound primeBound n : ℕ} {residues : List ℕ}
    (hn : 0 < n) (hle : n ≤ bound) (hsqrt : n < primeBound ^ 2)
    (hres : ∀ r ∈ residues, r < period)
    (hone : 1 % period ∈ residues)
    (hmul : ∀ a b, a % period ∈ residues → b % period ∈ residues → (a * b) % period ∈ residues)
    (hcancel : ∀ a b, (a * b) % period ∈ residues → b % period ∈ residues →
      a % period ∈ residues) :
    (periodicCoreMask period bound primeBound residues).testBit n = true ↔
      PrimeCore (fun p => p % period ∈ residues) n := by
  unfold periodicCoreMask
  simp only [Nat.testBit_and, Nat.testBit_two_pow_sub_one, decide_eq_true (by omega : n < bound + 1),
    Bool.true_and, Nat.testBit_xor, Bool.true_xor, Bool.and_eq_true]
  rw [periodicMask_bit hle hres]
  let bad := (List.range primeBound).filter (fun p => decide (p.Prime ∧ p % period ∉ residues))
  have hmem : ∀ p, p ∈ bad ↔ p < primeBound ∧ p.Prime ∧ p % period ∉ residues := fun _ => by simp [bad]
  have hbad := badMultiplesMask_bit hle bad (fun p hp => (hmem p |>.mp hp).2.1.pos)
  have notbit (b : Bool) : (!b) = true ↔ ¬ b = true := by cases b <;> decide
  rw [notbit]
  constructor
  · rintro ⟨hnres, hno⟩
    exact primeCore_of_small_test hone hmul hcancel hn hnres hsqrt
      (fun p hp hpp hpbad hpd => hno (hbad.mpr ⟨p, (hmem p).mpr ⟨hp, hpp, hpbad⟩, hpd⟩))
  · intro hc
    refine ⟨primeCore_induction hone (fun _ _ h => h) hmul hc, fun hset => ?_⟩
    obtain ⟨p, hp, hpd⟩ := hbad.mp hset
    obtain ⟨_, hprime, hmod⟩ := (hmem p).mp hp
    exact hmod (hc.2 p hprime hpd)

theorem periodicCoreMask_zero {period bound primeBound : ℕ} {residues : List ℕ}
    (hres : ∀ r ∈ residues, r < period) (hzero : 0 ∉ residues) :
    (periodicCoreMask period bound primeBound residues).testBit 0 = false := by
  have hz : (periodicMask period bound residues).testBit 0 = false := by
    apply Bool.eq_false_iff.mpr
    intro h
    exact hzero (by simpa only [Nat.zero_mod] using (periodicMask_bit (Nat.zero_le bound) hres).mp h)
  simp only [periodicCoreMask, Nat.testBit_and, hz, Bool.and_false, Bool.false_and]

end QuadraticNorms
