import Mathlib.Data.Finset.Sort
import Mathlib.Data.Nat.Factorization.Defs

open scoped BigOperators

namespace QuadraticNorms

private def prefixPower (primes exponents : List ℕ) : ℕ :=
  (primes.zipWith (· ^ ·) exponents).prod

private def divisorProduct (exponents : List ℕ) : ℕ :=
  (exponents.map (· + 1)).prod

/-- Maximize the divisor product on prefixes of the given primes.
Searching all positive exponent orders avoids an exponent rearrangement lemma. -/
def splitDivisorEnvelope (cap : ℕ) : List ℕ → ℕ → ℕ
  | [], _ => 1
  | p :: rest, X =>
    (List.range' 1 cap).foldl
      (fun best a => if p ^ a ≤ X then max best ((a + 1) * splitDivisorEnvelope cap rest (X / p ^ a)) else best)
      1

private theorem list_le_foldl_max {α : Type*} (score : α → ℕ)
    (valid : α → Prop) [DecidablePred valid] (items : List α) (base : ℕ) :
    base ≤ items.foldl (fun best a => if valid a then max best (score a) else best) base ∧
    ∀ a ∈ items, valid a →
      score a ≤ items.foldl (fun best a => if valid a then max best (score a) else best) base := by
  induction items generalizing base with
  | nil => simp
  | cons a tail ih =>
    let first := if valid a then max base (score a) else base
    have hfirst : base ≤ first := by dsimp [first]; split <;> simp
    refine ⟨hfirst.trans (ih first).1, ?_⟩
    intro b hb hv
    rcases List.mem_cons.mp hb with rfl | hm
    · have h : score b ≤ first := by simp [first, hv]
      exact h.trans (ih first).1
    · exact (ih first).2 b hm hv

private theorem prefixPower_pos (primes exps : List ℕ)
    (h : ∀ p ∈ primes, 0 < p) : 0 < prefixPower primes exps := by
  induction primes generalizing exps with
  | nil => simp [prefixPower]
  | cons p ps ih =>
    cases exps with
    | nil => simp [prefixPower]
    | cons a exps =>
      exact Nat.mul_pos (pow_pos (h p (by simp)) a) (ih exps (fun q hq => h q (by simp [hq])))

private theorem divisorProduct_le_envelope {cap X : ℕ} (primes exponents : List ℕ)
    (hprimes : ∀ p ∈ primes, 0 < p) (hlen : exponents.length ≤ primes.length)
    (hexps : ∀ a ∈ exponents, 0 < a ∧ a ≤ cap)
    (hprod : prefixPower primes exponents ≤ X) :
    divisorProduct exponents ≤ splitDivisorEnvelope cap primes X := by
  induction primes generalizing X exponents with
  | nil =>
    have he : exponents = [] := by simpa using hlen
    simp [he, divisorProduct, splitDivisorEnvelope]
  | cons p primes ih =>
    have hp := hprimes p (by simp)
    let valid := fun a : ℕ => p ^ a ≤ X
    let score := fun a : ℕ => (a + 1) * splitDivisorEnvelope cap primes (X / p ^ a)
    have hfold := list_le_foldl_max score valid (List.range' 1 cap) 1
    cases exponents with
    | nil => exact hfold.1
    | cons a rest =>
      have ha := hexps a (by simp)
      have help : a ∈ List.range' 1 cap := List.mem_range'_1.mpr ⟨ha.1, by omega⟩
      have hrest : prefixPower primes rest ≤ X / p ^ a := by
        apply (Nat.le_div_iff_mul_le (pow_pos hp a)).mpr
        simpa [prefixPower, Nat.mul_comm] using hprod
      have haX : p ^ a ≤ X := by
        have hs := prefixPower_pos primes rest (fun q hq => hprimes q (by simp [hq]))
        exact (Nat.le_mul_of_pos_right _ hs).trans (by simpa [prefixPower] using hprod)
      exact (Nat.mul_le_mul_left _ (ih rest
        (fun q hq => hprimes q (by simp [hq])) (by simpa using hlen)
        (fun b hb => hexps b (by simp [hb])) hrest)).trans (hfold.2 a help haX)


/-- Replace each prime in a sorted split factorization by the
corresponding prime in the initial split-prime prefix. -/
private theorem prefixPower_le_original (split : ℕ → Prop) {upper lower : ℕ} (qs ps : List ℕ)
    (hqs : qs.Pairwise (· < ·)) (hps : ps.Pairwise (· < ·))
    (hupper : ∀ q ∈ qs, lower < q ∧ q ≤ upper)
    (hlower : ∀ p ∈ ps, lower < p) (hps_split : ∀ p ∈ ps, split p)
    (hcomplete : ∀ p, lower < p → p ≤ upper → split p → p ∈ qs)
    (exponent : ℕ → ℕ) :
    prefixPower qs (ps.map exponent) ≤ (ps.map (fun p => p ^ exponent p)).prod := by
  induction qs generalizing ps lower with
  | nil =>
    have hp : 0 < ((ps.map (fun p => p ^ exponent p))).prod := by
      exact prefixPower_pos ps (ps.map exponent) (fun p hp => Nat.zero_lt_of_lt (hlower p hp))
        |> fun h => by simpa [prefixPower, List.zipWith_map_right, List.zipWith_self] using h
    exact hp
  | cons q qs ih =>
    cases ps with
    | nil => simp [prefixPower]
    | cons p ps =>
      have hq : q ≤ p := by
        by_contra! h
        have hmem := hcomplete p (hlower p (by simp)) (h.le.trans (hupper q (by simp)).2) (hps_split p (by simp))
        rcases List.mem_cons.mp hmem with hEq | hmem'
        · omega
        · exact (not_lt_of_ge h.le) ((List.pairwise_cons.mp hqs).1 p hmem')
      have hrec := ih ps (lower := q) (List.pairwise_cons.mp hqs).2 (List.pairwise_cons.mp hps).2
        (fun r hr => ⟨(List.pairwise_cons.mp hqs).1 r hr, (hupper r (by simp [hr])).2⟩)
        (fun r hr => hq.trans_lt ((List.pairwise_cons.mp hps).1 r hr))
        (fun r hr => hps_split r (by simp [hr]))
        (fun r hr hlo hs => by
          have hmem := hcomplete r ((hupper q (by simp)).1.trans hr) hlo hs
          exact (List.mem_cons.mp hmem).resolve_left (ne_of_gt hr))
      exact Nat.mul_le_mul (Nat.pow_le_pow_left hq _) hrec

private theorem primes_prod_le_prefixPower (qs exps : List ℕ)
    (hlen : qs.length ≤ exps.length) (hPrimes : ∀ q ∈ qs, 0 < q) (hExp : ∀ a ∈ exps, 0 < a) :
    qs.prod ≤ prefixPower qs exps := by
  induction qs generalizing exps with
  | nil => simp [prefixPower]
  | cons q qs ih =>
    cases exps with
    | nil => simp at hlen
    | cons a exps =>
      have hp : q ≤ q ^ a := by
        simpa only [pow_one] using Nat.pow_le_pow_right (hPrimes q (by simp))
          (show 1 ≤ a from hExp a (by simp))
      exact Nat.mul_le_mul hp (ih exps (by simpa using hlen)
        (fun p hp => hPrimes p (by simp [hp])) (fun b hb => hExp b (by simp [hb])))


/-- A prefix calculation bounds the divisor product at every integer
in the requested range. The prime bound and the exponent bound are
checked separately to ensure the search includes all factorizations. -/
theorem split_divisors_le_prefix_envelope (split : ℕ → Prop) [DecidablePred split]
    (prefixPrimes : List ℕ) {last base cap X n : ℕ} (hn : 0 < n) (hnX : n ≤ X)
    (hsorted : prefixPrimes.Pairwise (· < ·))
    (hprefix : ∀ q ∈ prefixPrimes, 0 < q ∧ q ≤ last)
    (hcomplete : ∀ p, p ≤ last → split p → p ∈ prefixPrimes)
    (hbase : 1 < base) (hprime : ∀ p, p.Prime → split p → base ≤ p)
    (hcap : X < base ^ (cap + 1)) (hprodPrefix : X < prefixPrimes.prod) :
    (∏ p ∈ n.primeFactors.filter split, (n.factorization p + 1)) ≤
      splitDivisorEnvelope cap prefixPrimes X := by
  let s := n.primeFactors.filter split
  let ps := s.sort
  let e := n.factorization
  let exps := ps.map e
  have hs (p : ℕ) (hp : p ∈ ps) : p ∈ n.primeFactors ∧ split p :=
    Finset.mem_filter.mp ((Finset.mem_sort (· ≤ ·)).mp hp)
  have hps : ps.Pairwise (· < ·) := (Finset.sortedLT_sort s).pairwise
  have hprod₀ : (ps.map (fun p => p ^ e p)).prod ≤ n := by
    rw [← List.prod_toFinset _ (Finset.sort_nodup ..), Finset.sort_toFinset]
    have hall := Nat.prod_factorization_pow_eq_self hn.ne'
    rw [Finsupp.prod, Nat.support_factorization] at hall
    rw [← hall]
    exact Finset.prod_le_prod_of_subset_of_one_le' (Finset.filter_subset _ _)
      (fun p hp _ => Nat.one_le_pow _ _ (Nat.prime_of_mem_primeFactors hp).pos)
  have hprod := (prefixPower_le_original split prefixPrimes ps hsorted hps hprefix
    (fun p hp => (Nat.prime_of_mem_primeFactors (hs p hp).1).pos) (fun p hp => (hs p hp).2)
    (fun p _ hp => hcomplete p hp) e).trans hprod₀
  have hexps (a : ℕ) (ha : a ∈ exps) : 0 < a ∧ a ≤ cap := by
    obtain ⟨p, hp, rfl⟩ := List.mem_map.mp ha
    obtain ⟨hpf, hpsplit⟩ := hs p hp
    have hp' := Nat.prime_of_mem_primeFactors hpf
    have hnp := Nat.le_of_dvd hn (Nat.ordProj_dvd n p)
    have hpow : base ^ e p < base ^ (cap + 1) :=
      ((Nat.pow_le_pow_left (hprime p hp' hpsplit) _).trans (hnp.trans hnX)).trans_lt hcap
    exact ⟨hp'.factorization_pos_of_dvd hn.ne' (Nat.dvd_of_mem_primeFactors hpf),
      Nat.le_of_lt_succ ((Nat.pow_lt_pow_iff_right hbase).mp hpow)⟩
  have hlen : exps.length ≤ prefixPrimes.length := by
    by_contra! hl
    have hsmall := (primes_prod_le_prefixPower prefixPrimes exps hl.le
      (fun q hq => (hprefix q hq).1) (fun a ha => (hexps a ha).1)).trans (hprod.trans hnX)
    exact (not_le_of_gt hprodPrefix) hsmall
  have hvalue : (∏ p ∈ n.primeFactors.filter split, (n.factorization p + 1)) = divisorProduct exps := by
    unfold divisorProduct
    rw [List.map_map, ← List.prod_toFinset _ (Finset.sort_nodup ..), Finset.sort_toFinset]
    rfl
  rw [hvalue]
  exact divisorProduct_le_envelope prefixPrimes exps (fun q hq => (hprefix q hq).1) hlen hexps
    (hprod.trans hnX)

end QuadraticNorms
