import Mathlib.Algebra.Order.Field.Rat
import Mathlib.Data.Finset.Card
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

open scoped BigOperators

namespace QuadraticNorms

/-- A normalized bad-count bound below one leaves a good candidate.
The finite collection can contain promotions of any individual bad set
into the exceptional set, as in the mixed Eisenstein case. -/
theorem exists_unspoiled {α ι : Type*} [DecidableEq α]
    {candidates exceptional : Finset α} {primes : Finset ι} {bad : ι → Finset α}
    (hne : candidates.Nonempty)
    (h : ((exceptional.card : ℚ) + ∑ p ∈ primes, (bad p).card) / candidates.card < 1) :
    ∃ a ∈ candidates, a ∉ exceptional ∧ ∀ p ∈ primes, a ∉ bad p := by
  classical
  have hpos : (0 : ℚ) < candidates.card := by exact_mod_cast Finset.card_pos.mpr hne
  have htotal : (exceptional.card : ℚ) + ∑ p ∈ primes, (bad p).card < candidates.card := by
    simpa using (div_lt_iff₀ hpos).mp h
  have hs : (exceptional ∪ primes.biUnion bad).card ≤
      exceptional.card + ∑ p ∈ primes, (bad p).card :=
    (Finset.card_union_le ..).trans (Nat.add_le_add_left (Finset.card_biUnion_le) _)
  have hcard : (exceptional ∪ primes.biUnion bad).card < candidates.card := by
    exact_mod_cast (show ((exceptional ∪ primes.biUnion bad).card : ℚ) < candidates.card by
      exact lt_of_le_of_lt (by exact_mod_cast hs) htotal)
  by_contra! hnot
  apply (not_le_of_gt hcard)
  apply Finset.card_le_card
  intro a ha
  by_cases he : a ∈ exceptional
  · exact Finset.mem_union_left _ he
  · obtain ⟨p, hp, hbad⟩ := hnot a ha he
    exact Finset.mem_union_right _ (Finset.mem_biUnion.mpr ⟨p, hp, hbad⟩)

/-- The rational margin used above the Eisenstein finite cutoff,
after inserting the three uniform power bounds. -/
theorem eisenstein_spoiled_fraction_bound {u v w : ℚ}
    (hu : u < 2 / 5) (hv : v < 1 / 250) (hw : w ≤ 1 / 10000) :
    29 / 60 + (11 / 10) * u + (111 / 20) * v + (2 / 3) * w < 591 / 625 := by
  linarith

theorem eisenstein_spoiled_fraction_less_than_one {u v w : ℚ}
    (hu : u < 2 / 5) (hv : v < 1 / 250) (hw : w ≤ 1 / 10000) :
    29 / 60 + (11 / 10) * u + (111 / 20) * v + (2 / 3) * w < 1 :=
  (eisenstein_spoiled_fraction_bound hu hv hw).trans (by norm_num)

end QuadraticNorms
