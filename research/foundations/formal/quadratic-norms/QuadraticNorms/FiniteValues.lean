import QuadraticNorms.LocalObstructions
import Mathlib.Data.Nat.Sqrt
import Mathlib.Algebra.Order.Ring.Abs
import Mathlib.Tactic.Linarith

namespace QuadraticNorms

/-- A complete coordinate search for the four monic forms. -/
private def coordinateGrid (bound : ℕ) : List ℤ :=
  let r := (2 * bound).sqrt
  (List.range (2 * r + 1)).map (fun (i : ℕ) => (i : ℤ) - r)

private theorem mem_coordinateGrid {bound : ℕ} {x : ℤ}
    (h : x.natAbs ^ 2 ≤ 2 * bound) : x ∈ coordinateGrid bound := by
  have hle : x.natAbs ≤ (2 * bound).sqrt := Nat.le_sqrt.mpr (by simpa [pow_two] using h)
  have habs : (x.natAbs : ℤ) = |x| := Int.natCast_natAbs x
  have hlow := neg_abs_le x
  have hhigh := le_abs_self x
  change x ∈ (List.range (2 * (2 * bound).sqrt + 1)).map
    (fun (i : ℕ) => (i : ℤ) - (2 * bound).sqrt)
  rw [List.mem_map]
  refine ⟨(x + ((2 * bound).sqrt : ℤ)).toNat, ?_, ?_⟩
  · rw [List.mem_range]
    omega
  · omega

private theorem coordinate_sq_bounds {b c x y : ℤ} {n bound : ℕ}
    (hb : b = -1 ∨ b = 0 ∨ b = 1) (hc : 1 ≤ c)
    (hn : form b c x y = n) (hle : n ≤ bound) :
    x.natAbs ^ 2 ≤ 2 * bound ∧ y.natAbs ^ 2 ≤ 2 * bound := by
  have hc₀ := mul_nonneg (sub_nonneg.mpr hc) (sq_nonneg y)
  have hplus := sq_nonneg (x + y)
  have hminus := sq_nonneg (x - y)
  have hysq := sq_nonneg y
  have hxsq := sq_nonneg x
  have hle' : (n : ℤ) ≤ bound := by exact_mod_cast hle
  have hx : (x.natAbs : ℤ) ^ 2 = x ^ 2 := by rw [Int.natCast_natAbs, sq_abs]
  have hy : (y.natAbs : ℤ) ^ 2 = y ^ 2 := by rw [Int.natCast_natAbs, sq_abs]
  have hbound : x ^ 2 ≤ 2 * (bound : ℤ) ∧ y ^ 2 ≤ 2 * (bound : ℤ) := by
    rcases hb with rfl | rfl | rfl <;> unfold form at hn <;> constructor <;> nlinarith
  constructor
  · exact_mod_cast hx.trans_le hbound.1
  · exact_mod_cast hy.trans_le hbound.2

/-- The list of positive primitive values in the given interval.
Any repeated coordinate witnesses give the same membership test. -/
def finitePrimitiveValues (b c : ℤ) (bound : ℕ) : List ℕ :=
  (coordinateGrid bound).flatMap fun x =>
    (coordinateGrid bound).filterMap fun y =>
      let q := form b c x y
      if Int.gcd x y = 1 ∧ 0 < q ∧ q ≤ bound then some q.toNat else none

theorem mem_finitePrimitiveValues {b c : ℤ} {n bound : ℕ}
    (hb : b = -1 ∨ b = 0 ∨ b = 1) (hc : 1 ≤ c) :
    n ∈ finitePrimitiveValues b c bound ↔ PrimitiveValues b c n ∧ n ≤ bound := by
  simp only [finitePrimitiveValues, List.mem_flatMap, List.mem_filterMap]
  constructor
  · rintro ⟨x, _, y, _, h⟩
    split at h
    · rename_i hq
      have hval := Option.some.inj h
      have hq' : (n : ℤ) = form b c x y := by omega
      exact ⟨⟨by omega, x, y, hq.1, hq'.symm⟩, by omega⟩
    · contradiction
  · rintro ⟨⟨hnpos, x, y, hgcd, hq⟩, hn⟩
    have hbounds := coordinate_sq_bounds hb hc hq hn
    refine ⟨x, mem_coordinateGrid hbounds.1, y, mem_coordinateGrid hbounds.2, ?_⟩
    have htest : Int.gcd x y = 1 ∧ 0 < form b c x y ∧ form b c x y ≤ bound :=
      ⟨hgcd, by omega, by omega⟩
    rw [if_pos htest, hq]
    simp

/-- Test every positive first summand at most the target against
the authenticated list of primitive values. -/
def finitePrimitiveSum (b c : ℤ) (n : ℕ) : Bool :=
  let values := finitePrimitiveValues b c n
  (List.range (n + 1)).any fun a => values.contains a && values.contains (n - a)

theorem finitePrimitiveSum_eq_true {b c : ℤ} {n : ℕ}
    (hb : b = -1 ∨ b = 0 ∨ b = 1) (hc : 1 ≤ c) :
    finitePrimitiveSum b c n = true ↔
      InSumset (PrimitiveValues b c) (PrimitiveValues b c) n := by
  simp only [finitePrimitiveSum, List.any_eq_true, List.mem_range, Bool.and_eq_true,
    List.contains_iff_mem, mem_finitePrimitiveValues hb hc]
  constructor
  · rintro ⟨a, ha, ⟨hP, _⟩, hQ, _⟩
    exact ⟨a, n - a, hP, hQ, by omega⟩
  · rintro ⟨a, d, ha, hd, rfl⟩
    exact ⟨a, by omega, ⟨ha, by omega⟩, by simpa using hd, by omega⟩

theorem P7_small_exceptions :
    ¬ InSumset P7 P7 7 ∧ ¬ InSumset P7 P7 77 ∧ ¬ InSumset P7 P7 203 := by
  have h (n : ℕ) := finitePrimitiveSum_eq_true (n := n) (b := 1) (c := 2)
    (by decide) (by decide)
  refine ⟨fun h₇ => ?_, fun h₇₇ => ?_, fun h₂₀₃ => ?_⟩
  · exact (by decide +kernel : finitePrimitiveSum 1 2 7 ≠ true) ((h 7).mpr h₇)
  · exact (by decide +kernel : finitePrimitiveSum 1 2 77 ≠ true) ((h 77).mpr h₇₇)
  · exact (by decide +kernel : finitePrimitiveSum 1 2 203 ≠ true) ((h 203).mpr h₂₀₃)

theorem P7_sum_necessary {n : ℕ} (h : InSumset P7 P7 n) :
    2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203 := by
  obtain ⟨a, b, ha, hb, hab⟩ := h
  have hsum : InSumset P7 P7 n := ⟨a, b, ha, hb, hab⟩
  have := ha.1
  have := hb.1
  exact ⟨by omega, hab ▸ P7_sum_not_dvd_49 ha hb,
    fun he => P7_small_exceptions.1 (he ▸ hsum),
    fun he => P7_small_exceptions.2.1 (he ▸ hsum),
    fun he => P7_small_exceptions.2.2 (he ▸ hsum)⟩

end QuadraticNorms
