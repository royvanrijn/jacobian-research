import QuadraticNorms.Definitions
import Mathlib.Data.Nat.Factorization.Induction
import Mathlib.Data.Nat.ModEq

namespace QuadraticNorms

theorem primeCore_prime {allowed : ℕ → Prop} {p : ℕ} (hp : p.Prime) :
    PrimeCore allowed p ↔ allowed p := by
  constructor
  · exact fun h => h.2 p hp dvd_rfl
  · intro h
    refine ⟨hp.pos, ?_⟩
    intro q hq hqp
    rcases (Nat.dvd_prime hp).mp hqp with rfl | rfl
    · exact (hq.ne_one rfl).elim
    · exact h

/-- A property closed under multiplication holds on the core once
it holds for one and the allowed primes. -/
theorem primeCore_induction {allowed P : ℕ → Prop}
    (hone : P 1) (hprime : ∀ p, p.Prime → allowed p → P p)
    (hmul : ∀ a b, P a → P b → P (a * b))
    {n : ℕ} (hn : PrimeCore allowed n) : P n := by
  refine Nat.recOnMul (motive := fun k => PrimeCore allowed k → P k) ?_ ?_ ?_ ?_ n hn
  · exact fun h => (Nat.lt_irrefl _ h.1).elim
  · exact fun _ => hone
  · exact fun p hp h => hprime p hp ((primeCore_prime hp).mp h)
  · intro a b hPa hPb hab
    exact hmul a b (hPa (primeCore_of_dvd hab (dvd_mul_right a b)))
      (hPb (primeCore_of_dvd hab (dvd_mul_left b a)))

theorem C3_mod_six {n : ℕ} (hn : C3 n) : n % 6 = 1 := by
  apply primeCore_induction (P := fun n => n ≡ 1 [MOD 6]) (by decide) ?_ ?_ hn
  · intro p hp h
    have hpodd := hp.eq_two_or_odd.resolve_left (by omega)
    change p % 6 = 1
    have h₃ := Nat.mod_mod_of_dvd p (by decide : 3 ∣ 6)
    have h₂ := Nat.mod_mod_of_dvd p (by decide : 2 ∣ 6)
    rw [h] at h₃
    rw [hpodd] at h₂
    omega
  · intro a b ha hb
    simpa using ha.mul hb

theorem C12_mod_twelve {n : ℕ} (hn : C12 n) : n % 12 = 1 := by
  apply primeCore_induction (P := fun n => n ≡ 1 [MOD 12]) (by decide) ?_ ?_ hn
  · exact fun _ _ h => h
  · intro a b ha hb
    simpa using ha.mul hb

theorem C8_mod_eight {n : ℕ} (hn : C8 n) : n % 8 = 1 ∨ n % 8 = 3 := by
  apply primeCore_induction (P := fun n => n ≡ 1 [MOD 8] ∨ n ≡ 3 [MOD 8])
    (by decide) ?_ ?_ hn
  · exact fun _ _ h => h
  · intro a b ha hb
    rcases ha with ha | ha <;> rcases hb with hb | hb
    · exact Or.inl (by simpa using ha.mul hb)
    · exact Or.inr (by simpa using ha.mul hb)
    · exact Or.inr (by simpa using ha.mul hb)
    · exact Or.inl ((ha.mul hb).trans (by decide))

/-- Closure under scaling in either the ordinary or mixed norm problem. -/
theorem coreSum_scale {allowed : ℕ → Prop} {k g n : ℕ}
    (hg : PrimeCore allowed g)
    (h : InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) n) :
    InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) (g * n) := by
  apply InSumset.scale (fun _ ha => primeCore_mul hg ha) _ h
  rintro _ ⟨a, ha, rfl⟩
  exact ⟨g * a, primeCore_mul hg ha, by ac_rfl⟩

/-- Strong induction for core sumsets using the exact split-prime
quotient and the cases in which there is no split prime factor. -/
theorem coreSum_of_prime_descent {allowed eligible : ℕ → Prop} {k : ℕ}
    (hreduced : ∀ n, eligible n →
      (∀ p, p.Prime → p ∣ n → ¬ allowed p) →
      InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) n)
    (hstep : ∀ p n, p.Prime → allowed p → p ∣ n → eligible n →
      eligible (n / p) ∨ InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) n)
    {n : ℕ} (hnpos : 0 < n) (hn : eligible n) :
    InSumset (PrimeCore allowed) (Dilate k (PrimeCore allowed)) n := by
  induction n using Nat.strong_induction_on with
  | h n ih =>
    by_cases hex : ∃ p, p.Prime ∧ p ∣ n ∧ allowed p
    · obtain ⟨p, hp, hd, ha⟩ := hex
      rcases hstep p n hp ha hd hn with hsmall | hdone
      · have hpos := Nat.div_pos (Nat.le_of_dvd hnpos hd) hp.pos
        have hlt := Nat.div_lt_self hnpos hp.one_lt
        have hscaled := coreSum_scale ((primeCore_prime hp).mpr ha)
          (ih _ hlt hpos hsmall)
        simpa [Nat.mul_div_cancel' hd] using hscaled
      · exact hdone
    · exact hreduced n hn (fun p hp hd ha => hex ⟨p, hp, hd, ha⟩)

end QuadraticNorms
