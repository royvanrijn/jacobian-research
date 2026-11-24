import Mathlib.Data.Int.GCD
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

namespace QuadraticNorms

/-- The integral binary quadratic form x² + bxy + cy². -/
def form (b c x y : ℤ) : ℤ := x ^ 2 + b * x * y + c * y ^ 2

/-- Positive numerical values at primitive integer vectors. -/
def PrimitiveValues (b c : ℤ) (n : ℕ) : Prop :=
  0 < n ∧ ∃ x y : ℤ, Int.gcd x y = 1 ∧ form b c x y = n

def P3 := PrimitiveValues (-1) 1
def P7 := PrimitiveValues 1 2
def P8 := PrimitiveValues 0 2
def P11 := PrimitiveValues 1 3

/-- The positive integers with no prime divisor outside the allowed primes. -/
def PrimeCore (allowed : ℕ → Prop) (n : ℕ) : Prop :=
  0 < n ∧ ∀ p : ℕ, p.Prime → p ∣ n → allowed p

def C3 := PrimeCore (fun p => p % 3 = 1)
def C7 := PrimeCore (fun p => p % 7 = 1 ∨ p % 7 = 2 ∨ p % 7 = 4)
def C8 := PrimeCore (fun p => p % 8 = 1 ∨ p % 8 = 3)
def C11 := PrimeCore (fun p => p % 11 = 1 ∨ p % 11 = 3 ∨ p % 11 = 4 ∨
  p % 11 = 5 ∨ p % 11 = 9)
def C12 := PrimeCore (fun p => p % 12 = 1)

theorem primeCore_one (allowed : ℕ → Prop) : PrimeCore allowed 1 := by
  refine ⟨Nat.one_pos, ?_⟩
  intro p hp hd
  exact (hp.ne_one (Nat.dvd_one.mp hd)).elim

theorem primeCore_mul {allowed : ℕ → Prop} {a b : ℕ}
    (ha : PrimeCore allowed a) (hb : PrimeCore allowed b) :
    PrimeCore allowed (a * b) := by
  refine ⟨Nat.mul_pos ha.1 hb.1, ?_⟩
  intro p hp hd
  rcases hp.dvd_mul.mp hd with h | h
  · exact ha.2 p hp h
  · exact hb.2 p hp h

theorem primeCore_of_dvd {allowed : ℕ → Prop} {a b : ℕ}
    (hb : PrimeCore allowed b) (hab : a ∣ b) : PrimeCore allowed a := by
  refine ⟨Nat.pos_of_dvd_of_pos hab hb.1, ?_⟩
  exact fun p hp hpa => hb.2 p hp (hpa.trans hab)

/-- Dilation by a fixed coefficient. -/
def Dilate (k : ℕ) (P : ℕ → Prop) (n : ℕ) : Prop :=
  ∃ a, P a ∧ k * a = n

/-- The binary sumset of two predicates on numerical values. -/
def InSumset (P Q : ℕ → Prop) (n : ℕ) : Prop :=
  ∃ a b, P a ∧ Q b ∧ a + b = n

/-- Ordered value pairs are indexed by their first coordinate, with the
second coordinate determined by the target. -/
def pairFinset (P Q : ℕ → Prop) [DecidablePred P] [DecidablePred Q] (n : ℕ) :
    Finset ℕ :=
  (Finset.range (n + 1)).filter (fun a => P a ∧ Q (n - a))

def pairCount (P Q : ℕ → Prop) [DecidablePred P] [DecidablePred Q] (n : ℕ) : ℕ :=
  (pairFinset P Q n).card

theorem mem_pairFinset {P Q : ℕ → Prop} [DecidablePred P] [DecidablePred Q]
    {n a : ℕ} :
    a ∈ pairFinset P Q n ↔ a ≤ n ∧ P a ∧ Q (n - a) := by
  simp [pairFinset, Nat.lt_succ_iff]

theorem pairCount_pos_iff {P Q : ℕ → Prop} [DecidablePred P] [DecidablePred Q]
    {n : ℕ} :
    0 < pairCount P Q n ↔ InSumset P Q n := by
  rw [pairCount, Finset.card_pos]
  constructor
  · rintro ⟨a, ha⟩
    obtain ⟨hle, hP, hQ⟩ := mem_pairFinset.mp ha
    exact ⟨a, n - a, hP, hQ, by omega⟩
  · rintro ⟨a, b, hP, hQ, hab⟩
    exact ⟨a, mem_pairFinset.mpr ⟨by omega, hP, by simpa [← hab] using hQ⟩⟩

theorem InSumset.scale {P Q : ℕ → Prop} {g n : ℕ}
    (hP : ∀ a, P a → P (g * a)) (hQ : ∀ b, Q b → Q (g * b))
    (h : InSumset P Q n) : InSumset P Q (g * n) := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  exact ⟨g * a, g * b, hP a ha, hQ b hb, (Nat.mul_add ..).symm⟩

end QuadraticNorms
