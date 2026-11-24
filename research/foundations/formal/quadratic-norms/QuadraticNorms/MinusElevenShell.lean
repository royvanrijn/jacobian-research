import QuadraticNorms.MinusElevenSieveMargin
import QuadraticNorms.FiniteCoverage
import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.Data.Nat.Sqrt

open scoped BigOperators

namespace QuadraticNorms

/-- A field component of the norm shell, with signed coordinates. -/
abbrev MinusElevenVector := ℤ × ℤ

/-- A point of the quaternary norm shell, ordered by its two field components. -/
abbrev MinusElevenShellPoint := MinusElevenVector × MinusElevenVector

def minusElevenNorm (v : MinusElevenVector) : ℕ :=
  (form 1 3 v.1 v.2).toNat

theorem minusEleven_norm_cast (v : MinusElevenVector) :
    (minusElevenNorm v : ℤ) = form 1 3 v.1 v.2 := by
  apply Int.toNat_of_nonneg
  unfold form
  nlinarith [sq_nonneg (v.1 + v.2), sq_nonneg v.1, sq_nonneg v.2]

def minusElevenCoordinateBox (n : ℕ) : Finset MinusElevenVector :=
  let bound := Nat.sqrt (2 * n)
  let coordinates := (Finset.range (2 * bound + 1)).map
    ⟨fun k : ℕ => (k : ℤ) - (bound : ℤ), fun a b h => by dsimp at h; omega⟩
  coordinates ×ˢ coordinates

theorem mem_minusEleven_coordinate_box {n : ℕ} {v : MinusElevenVector} (hle : minusElevenNorm v ≤ n) :
    v ∈ minusElevenCoordinateBox n := by
  have h := minusEleven_norm_cast v
  have hN : (minusElevenNorm v : ℤ) ≤ n := by exact_mod_cast hle
  have hx : v.1 ^ 2 ≤ 2 * (n : ℤ) := by
    unfold form at h
    nlinarith [sq_nonneg (v.1 + v.2), sq_nonneg v.2]
  have hy : v.2 ^ 2 ≤ 2 * (n : ℤ) := by
    unfold form at h
    nlinarith [sq_nonneg (v.1 + v.2), sq_nonneg v.1]
  have hbound (z : ℤ) (hz : z ^ 2 ≤ 2 * (n : ℤ)) :
      -(Nat.sqrt (2 * n) : ℤ) ≤ z ∧ z ≤ (Nat.sqrt (2 * n) : ℤ) := by
    have hsq : z.natAbs ^ 2 ≤ 2 * n := by
      exact_mod_cast (Int.natAbs_sq z).trans_le hz
    have hle := Nat.le_sqrt'.mpr hsq
    have hcast : (z.natAbs : ℤ) ≤ (Nat.sqrt (2 * n) : ℤ) := by exact_mod_cast hle
    simpa only [Int.natCast_natAbs, abs_le] using hcast
  unfold minusElevenCoordinateBox
  apply Finset.mem_product.mpr
  have hmem (z : ℤ) (hz : -(Nat.sqrt (2 * n) : ℤ) ≤ z ∧ z ≤ (Nat.sqrt (2 * n) : ℤ)) :
      ∃ k ∈ Finset.range (2 * Nat.sqrt (2 * n) + 1), (k : ℤ) - (Nat.sqrt (2 * n) : ℤ) = z := by
    refine ⟨(z + (Nat.sqrt (2 * n) : ℤ)).toNat, Finset.mem_range.mpr ?_, ?_⟩ <;> omega
  exact ⟨Finset.mem_map.mpr (hmem v.1 (hbound v.1 hx)),
    Finset.mem_map.mpr (hmem v.2 (hbound v.2 hy))⟩

/-- The complete odd shell counted by G(n) in the −11 proof. -/
noncomputable def minusElevenOddShell (n : ℕ) : Finset MinusElevenShellPoint :=
  (minusElevenCoordinateBox n ×ˢ minusElevenCoordinateBox n).filter (fun v =>
    minusElevenNorm v.1 % 2 = 1 ∧ minusElevenNorm v.2 % 2 = 1 ∧
      minusElevenNorm v.1 + minusElevenNorm v.2 = n)

theorem mem_minusEleven_odd_shell {n : ℕ} {v : MinusElevenShellPoint} :
    v ∈ minusElevenOddShell n ↔
      minusElevenNorm v.1 % 2 = 1 ∧ minusElevenNorm v.2 % 2 = 1 ∧
        minusElevenNorm v.1 + minusElevenNorm v.2 = n := by
  constructor
  · exact fun h => (Finset.mem_filter.mp h).2
  · intro h
    have h₁ : minusElevenNorm v.1 ≤ n := by omega
    have h₂ : minusElevenNorm v.2 ≤ n := by omega
    exact Finset.mem_filter.mpr ⟨Finset.mem_product.mpr
      ⟨mem_minusEleven_coordinate_box h₁, mem_minusEleven_coordinate_box h₂⟩, h⟩

/-- The full quaternary shell counted by r_F(n). -/
noncomputable def minusElevenFullShell (n : ℕ) : Finset MinusElevenShellPoint :=
  (minusElevenCoordinateBox n ×ˢ minusElevenCoordinateBox n).filter (fun v => minusElevenNorm v.1 + minusElevenNorm v.2 = n)

theorem mem_minusEleven_full_shell {n : ℕ} {v : MinusElevenShellPoint} :
    v ∈ minusElevenFullShell n ↔ minusElevenNorm v.1 + minusElevenNorm v.2 = n := by
  constructor
  · exact fun h => (Finset.mem_filter.mp h).2
  · intro h
    have ha : minusElevenNorm v.1 ≤ n := by omega
    have hb : minusElevenNorm v.2 ≤ n := by omega
    exact Finset.mem_filter.mpr ⟨Finset.mem_product.mpr
      ⟨mem_minusEleven_coordinate_box ha, mem_minusEleven_coordinate_box hb⟩, h⟩

/-- Even norms are exactly the norms of vectors divisible by two. -/
theorem minusEleven_norm_even_iff (v : MinusElevenVector) :
    2 ∣ minusElevenNorm v ↔ (2 : ℤ) ∣ v.1 ∧ (2 : ℤ) ∣ v.2 := by
  have hcase : ∀ x y : ZMod 2, x ^ 2 + x * y + y ^ 2 = 0 ↔ x = 0 ∧ y = 0 := by decide
  have hnorm : (minusElevenNorm v : ZMod 2) =
      (v.1 : ZMod 2) ^ 2 + (v.1 : ZMod 2) * (v.2 : ZMod 2) + (v.2 : ZMod 2) ^ 2 := by
    rw [← Int.cast_natCast, minusEleven_norm_cast]
    have hthree : (3 : ZMod 2) = 1 := by decide
    simp only [form, Int.cast_add, Int.cast_pow, Int.cast_mul,
      Int.cast_ofNat, hthree, one_mul]
  rw [← ZMod.natCast_eq_zero_iff, hnorm, hcase, ZMod.intCast_zmod_eq_zero_iff_dvd,
    ZMod.intCast_zmod_eq_zero_iff_dvd]
  norm_cast

/-- Inert prime divisibility of a field norm forces both coordinates
to be divisible by that prime. -/
theorem minusEleven_inert_dvd_norm_iff {p : ℕ} (hp : MinusElevenInertPrime p)
    (v : MinusElevenVector) :
    p ∣ minusElevenNorm v ↔ (p : ℤ) ∣ v.1 ∧ (p : ℤ) ∣ v.2 := by
  letI : Fact p.Prime := ⟨hp.1⟩
  let x : ZMod p := v.1
  let y : ZMod p := v.2
  have hnorm : (minusElevenNorm v : ZMod p) = x ^ 2 + x * y + 3 * y ^ 2 := by
    rw [← Int.cast_natCast, minusEleven_norm_cast]
    simp [form, x, y]
  have hp11 : p ≠ 11 := by
    rintro rfl
    exact (by decide : ¬ 11 % 11 ∈ [2, 6, 7, 8, 10]) hp.2.2
  have hnosplit : ¬ HasNormRoot 1 3 p := by
    rw [hasNormRoot11_prime_iff hp.1 hp11]
    have h := hp.2.2
    simp only [List.mem_cons, List.not_mem_nil, or_false] at h
    omega
  rw [← ZMod.natCast_eq_zero_iff, hnorm, ← ZMod.intCast_zmod_eq_zero_iff_dvd,
    ← ZMod.intCast_zmod_eq_zero_iff_dvd]
  change x ^ 2 + x * y + 3 * y ^ 2 = 0 ↔ x = 0 ∧ y = 0
  constructor
  · intro h
    by_cases hy : y = 0
    · exact ⟨by simpa [hy] using h, hy⟩
    · have ht : (-(x / y)) ^ 2 - 1 * (-(x / y)) + 3 = 0 := by
        field_simp [hy]
        linear_combination h
      have hroot : HasNormRoot 1 3 p := ⟨-(x / y), by simpa using ht⟩
      exact (hnosplit hroot).elim
  · rintro ⟨hx, hy⟩
    simp [hx, hy]

def minusElevenScaleVector (a : ℕ) (v : MinusElevenVector) : MinusElevenVector :=
  ((a : ℤ) * v.1, (a : ℤ) * v.2)

theorem minusEleven_norm_scale (a : ℕ) (v : MinusElevenVector) :
    minusElevenNorm (minusElevenScaleVector a v) = a ^ 2 * minusElevenNorm v := by
  have heq : (minusElevenNorm (minusElevenScaleVector a v) : ℤ) = (a : ℤ) ^ 2 * (minusElevenNorm v : ℤ) := by
    rw [minusEleven_norm_cast, minusEleven_norm_cast]
    unfold minusElevenScaleVector form
    ring
  exact_mod_cast heq

private def scaleShellPoint (a : ℕ) (ha : 0 < a) : MinusElevenShellPoint ↪ MinusElevenShellPoint where
  toFun v := (minusElevenScaleVector a v.1, minusElevenScaleVector a v.2)
  inj' := by
    rintro ⟨⟨p, q⟩, r, s⟩ ⟨⟨w, x⟩, y, z⟩ h
    simp only [minusElevenScaleVector, Prod.mk.injEq] at h ⊢
    have ha' : (a : ℤ) ≠ 0 := by exact_mod_cast ha.ne'
    exact ⟨⟨mul_left_cancel₀ ha' h.1.1, mul_left_cancel₀ ha' h.1.2⟩,
      ⟨mul_left_cancel₀ ha' h.2.1, mul_left_cancel₀ ha' h.2.2⟩⟩

private abbrev doubleVector := minusElevenScaleVector 2

private theorem norm_doubleVector (v : MinusElevenVector) :
    minusElevenNorm (doubleVector v) = 4 * minusElevenNorm v := minusEleven_norm_scale 2 v

private abbrev doubleShellPoint := scaleShellPoint 2 (by decide)

private theorem even_shell_scales {v : MinusElevenShellPoint}
    (h₁ : 2 ∣ minusElevenNorm v.1) (h₂ : 2 ∣ minusElevenNorm v.2) :
    ∃ w, doubleShellPoint w = v := by
  rcases v with ⟨⟨a, b⟩, c, d⟩
  obtain ⟨⟨x, rfl⟩, ⟨y, rfl⟩⟩ := (minusEleven_norm_even_iff (a, b)).mp h₁
  obtain ⟨⟨z, rfl⟩, ⟨w, rfl⟩⟩ := (minusEleven_norm_even_iff (c, d)).mp h₂
  exact ⟨((x, y), (z, w)), rfl⟩

private theorem full_shell_is_odd_or_doubled {n : ℕ} (he : n % 2 = 0) :
    minusElevenFullShell n = minusElevenOddShell n ∪
      (if 4 ∣ n then (minusElevenFullShell (n / 4)).map doubleShellPoint else ∅) := by
  ext v
  rw [Finset.mem_union, mem_minusEleven_full_shell, mem_minusEleven_odd_shell]
  constructor
  · intro hv
    by_cases h : minusElevenNorm v.1 % 2 = 1
    · exact Or.inl ⟨h, by omega, hv⟩
    · obtain ⟨w, rfl⟩ := even_shell_scales (v := v) (by omega) (by omega)
      have hsum : 4 * (minusElevenNorm w.1 + minusElevenNorm w.2) = n := by
        change minusElevenNorm (doubleVector w.1) + minusElevenNorm (doubleVector w.2) = n at hv
        rw [norm_doubleVector, norm_doubleVector] at hv
        omega
      have hd : 4 ∣ n := ⟨_, hsum.symm⟩
      exact Or.inr (by
        rw [if_pos hd]
        exact Finset.mem_map.mpr ⟨w, mem_minusEleven_full_shell.mpr (by omega), rfl⟩)
  · intro h
    rcases h with ⟨_, _, hsum⟩ | h
    · exact hsum
    · by_cases hd : 4 ∣ n
      · rw [if_pos hd] at h
        obtain ⟨w, hw, rfl⟩ := Finset.mem_map.mp h
        change minusElevenNorm (doubleVector w.1) + minusElevenNorm (doubleVector w.2) = n
        rw [norm_doubleVector, norm_doubleVector, ← Nat.mul_add, mem_minusEleven_full_shell.mp hw]
        exact Nat.mul_div_cancel' hd
      · simp [hd] at h

theorem minusEleven_full_shell_card_add {n : ℕ} (he : n % 2 = 0) :
    (minusElevenFullShell n).card =
      (minusElevenOddShell n).card + (if 4 ∣ n then (minusElevenFullShell (n / 4)).card else 0) := by
  have hdis : Disjoint (minusElevenOddShell n)
      (if 4 ∣ n then (minusElevenFullShell (n / 4)).map doubleShellPoint else ∅) := by
    by_cases hd : 4 ∣ n
    · rw [if_pos hd, Finset.disjoint_left]
      intro v hv heven
      obtain ⟨w, _, rfl⟩ := Finset.mem_map.mp heven
      have hodd := (mem_minusEleven_odd_shell.mp hv).1
      change minusElevenNorm (doubleVector w.1) % 2 = 1 at hodd
      rw [norm_doubleVector] at hodd
      omega
    · simp [hd]
  rw [full_shell_is_odd_or_doubled he, Finset.card_union_of_disjoint hdis]
  by_cases hd : 4 ∣ n <;> simp [hd]

/-- Subtracting the shell at the integral index n/4 gives precisely the
odd shell at an even index. The subtracted count is zero when 4 ∤ n. -/
theorem minusEleven_odd_shell_card {n : ℕ} (he : n % 2 = 0) :
    (minusElevenOddShell n).card =
      (minusElevenFullShell n).card - (if 4 ∣ n then (minusElevenFullShell (n / 4)).card else 0) := by
  rw [minusEleven_full_shell_card_add he, Nat.add_sub_cancel]

noncomputable def minusElevenBadShell (n p : ℕ) : Finset MinusElevenShellPoint :=
  (minusElevenOddShell n).filter (fun v => p ∣ minusElevenNorm v.1 ∨ p ∣ minusElevenNorm v.2)

private theorem odd_shell_scaled_iff {a n : ℕ} (ha : 0 < a) (hodd : a % 2 = 1)
    (v : MinusElevenShellPoint) :
    scaleShellPoint a ha v ∈ minusElevenOddShell n ↔
      minusElevenNorm v.1 % 2 = 1 ∧ minusElevenNorm v.2 % 2 = 1 ∧
        a ^ 2 * (minusElevenNorm v.1 + minusElevenNorm v.2) = n := by
  rw [mem_minusEleven_odd_shell]
  change minusElevenNorm (minusElevenScaleVector a v.1) % 2 = 1 ∧
    minusElevenNorm (minusElevenScaleVector a v.2) % 2 = 1 ∧
      minusElevenNorm (minusElevenScaleVector a v.1) + minusElevenNorm (minusElevenScaleVector a v.2) = n ↔ _
  have hmod (x : ℕ) : (a ^ 2 * x) % 2 = x % 2 := by
    rw [Nat.mul_mod, Nat.pow_mod, hodd]
    simp
  rw [minusEleven_norm_scale, minusEleven_norm_scale, hmod, hmod, ← Nat.mul_add]

/-- When an inert prime divides the target, any spoiled point is scaled
by p in all four coordinates. Thus its count is a smaller odd shell,
or zero when the target is indivisible by p². -/
theorem minusEleven_bad_shell_card_of_dvd {p n : ℕ} (hp : MinusElevenInertPrime p) (hpn : p ∣ n) :
    (minusElevenBadShell n p).card =
      if p ^ 2 ∣ n then (minusElevenOddShell (n / p ^ 2)).card else 0 := by
  let scale := scaleShellPoint p hp.1.pos
  have hodd : p % 2 = 1 :=
    hp.1.mod_two_eq_one_iff_ne_two.mpr (by have := hp.2.1; omega)
  suffices heq : minusElevenBadShell n p =
      if p ^ 2 ∣ n then (minusElevenOddShell (n / p ^ 2)).map scale else ∅ by
    rw [heq]
    by_cases hd : p ^ 2 ∣ n <;> simp [hd]
  ext v
  constructor
  · intro hv
    obtain ⟨hshell, hbad⟩ := Finset.mem_filter.mp hv
    have hsum := (mem_minusEleven_odd_shell.mp hshell).2.2
    have hd : p ∣ minusElevenNorm v.1 ∧ p ∣ minusElevenNorm v.2 := by
      rcases hbad with h | h
      · exact ⟨h, (Nat.dvd_add_iff_right h).mpr (hsum.symm ▸ hpn)⟩
      · exact ⟨(Nat.dvd_add_iff_left h).mpr (hsum.symm ▸ hpn), h⟩
    rcases v with ⟨⟨a, b⟩, c, d⟩
    obtain ⟨⟨u, rfl⟩, ⟨v, rfl⟩⟩ := (minusEleven_inert_dvd_norm_iff hp (a, b)).mp hd.1
    obtain ⟨⟨x, rfl⟩, ⟨y, rfl⟩⟩ := (minusEleven_inert_dvd_norm_iff hp (c, d)).mp hd.2
    let w : MinusElevenShellPoint := ((u, v), (x, y))
    have hw := (odd_shell_scaled_iff hp.1.pos hodd w).mp hshell
    have hd₂ : p ^ 2 ∣ n := ⟨_, hw.2.2.symm⟩
    rw [if_pos hd₂]
    refine Finset.mem_map.mpr ⟨w, mem_minusEleven_odd_shell.mpr ⟨hw.1, hw.2.1, ?_⟩, rfl⟩
    rw [← hw.2.2, Nat.mul_div_cancel_left _ (pow_pos hp.1.pos _)]
  · intro hv
    by_cases hd : p ^ 2 ∣ n
    · rw [if_pos hd] at hv
      obtain ⟨w, hw, rfl⟩ := Finset.mem_map.mp hv
      have h := mem_minusEleven_odd_shell.mp hw
      refine Finset.mem_filter.mpr ⟨(odd_shell_scaled_iff hp.1.pos hodd w).mpr
        ⟨h.1, h.2.1, h.2.2 ▸ Nat.mul_div_cancel' hd⟩, Or.inl ?_⟩
      change p ∣ minusElevenNorm (minusElevenScaleVector p w.1)
      rw [minusEleven_norm_scale]
      exact (dvd_pow_self p (by decide : (2 : ℕ) ≠ 0)).trans (dvd_mul_right _ _)
    · simp [hd] at hv

private theorem split_of_prime {p : ℕ} (hp : p.Prime) :
    MinusElevenSplitPrime p ↔ C11 p := by
  rw [C11, primeCore_prime hp]
  simp [MinusElevenSplitPrime, hp]

private theorem minusEleven_prime_cases {p : ℕ} (hp : p.Prime) :
    p = 2 ∨ p = 11 ∨ MinusElevenSplitPrime p ∨ MinusElevenInertPrime p := by
  by_cases h₂ : p = 2
  · exact Or.inl h₂
  by_cases h₁₁ : p = 11
  · exact Or.inr (Or.inl h₁₁)
  have hmod : p % 11 ≠ 0 := by
    intro h
    rcases (Nat.dvd_prime hp).mp (Nat.dvd_iff_mod_eq_zero.mpr h) with h | h
    · omega
    · exact h₁₁ h.symm
  have hall : ∀ r : Fin 11, r.val ≠ 0 →
      r.val ∈ [1, 3, 4, 5, 9] ∨ r.val ∈ [2, 6, 7, 8, 10] := by decide
  rcases hall ⟨p % 11, Nat.mod_lt _ (by decide)⟩ hmod with h | h
  · exact Or.inr (Or.inr (Or.inl ⟨hp, h⟩))
  · exact Or.inr (Or.inr (Or.inr ⟨hp, lt_of_le_of_ne hp.two_le (Ne.symm h₂), h⟩))

private theorem norm_mem_core {a : ℕ} (ha : a % 2 = 1) (h₁₁ : ¬ 11 ∣ a)
    (hi : ∀ p, MinusElevenInertPrime p → ¬ p ∣ a) : C11 a := by
  refine ⟨by omega, fun p hp hd => ?_⟩
  rcases minusEleven_prime_cases hp with rfl | rfl | hs | hp'
  · omega
  · exact (h₁₁ hd).elim
  · exact (primeCore_prime hp).mp (split_of_prime hp |>.mp hs)
  · exact (hi p hp' hd).elim

/-- The shell lower bound and the three bad-count bounds needed by the
−11 sieve. All counts refer to the actual signed-coordinate shell. -/
structure MinusElevenShellBounds (n : ℕ) (error : ℝ) : Prop where
  shell_lower :
    (18 / 5 : ℝ) * n * (1 - error) ≤ (minusElevenOddShell n).card
  ramified_upper :
    (1 - error) * (minusElevenBadShell n 11).card ≤
      (1 / 3 + 2 * error) * (minusElevenOddShell n).card
  inert_small_upper : ∀ p, MinusElevenInertPrime p →
    (1 - error) * (minusElevenBadShell n p).card ≤
      (minusElevenLocalMass p + (5 / 2) * error) * (minusElevenOddShell n).card
  inert_large_upper : ∀ p, MinusElevenInertPrime p →
    ((minusElevenBadShell n p).card : ℝ) ≤
      (48 / 5) * minusElevenDivisorMaximum n * n / (p : ℝ) ^ 2

private theorem good_shell_point {n : ℕ} {v : MinusElevenShellPoint}
    (hv : v ∈ minusElevenOddShell n) (h₁₁ : v ∉ minusElevenBadShell n 11)
    (hi : ∀ p ≤ n, MinusElevenInertPrime p → v ∉ minusElevenBadShell n p) :
    InSumset C11 C11 n := by
  have ⟨h₁, h₂, hsum⟩ := mem_minusEleven_odd_shell.mp hv
  have havoids (p : ℕ) (h : v ∉ minusElevenBadShell n p) :
      ¬ p ∣ minusElevenNorm v.1 ∧ ¬ p ∣ minusElevenNorm v.2 := by
    rw [← not_or]
    exact fun hd => h (Finset.mem_filter.mpr ⟨hv, hd⟩)
  have hno (p a : ℕ) (hp : MinusElevenInertPrime p) (ha : a ≤ n) (hd₁ : 0 < a)
      (hd₂ : p ∣ a) : v ∉ minusElevenBadShell n p := by
    have hbound : p ≤ n := by
      exact (Nat.le_of_dvd hd₁ hd₂).trans ha
    exact hi p hbound hp
  refine ⟨minusElevenNorm v.1, minusElevenNorm v.2, ?_, ?_, hsum⟩
  · exact norm_mem_core h₁ (havoids 11 h₁₁).1
      (fun p hp hd => (havoids p (hno p _ hp (by omega) (by omega : 0 < minusElevenNorm v.1) hd)).1 hd)
  · exact norm_mem_core h₂ (havoids 11 h₁₁).2
      (fun p hp hd => (havoids p (hno p _ hp (by omega) (by omega : 0 < minusElevenNorm v.2) hd)).2 hd)


/-- The sieve margin turns the shell and bad-count inequalities into
a core representation. The bad sets are finite subsets of the actual shell. -/
theorem minusEleven_core_sum_of_shell_bounds {n P : ℕ} (hn : 0 < n) {E : ℝ} (hE : 0 ≤ E)
    (hB : MinusElevenShellBounds n E) (hmargin : minusElevenSieveExpression n P E < 2 / 3) :
    InSumset C11 C11 n := by
  let s := (Finset.range (n + 1)).filter MinusElevenInertPrime
  let small := s.filter (fun p => p ≤ P)
  let large := s.filter (fun p => ¬ p ≤ P)
  let G : ℝ := (minusElevenOddShell n).card
  let D : ℝ := minusElevenDivisorMaximum n
  let B (p : ℕ) : ℝ := (minusElevenBadShell n p).card
  have hs (p : ℕ) (hp : p ∈ small) : MinusElevenInertPrime p ∧ p ≤ P := by
    obtain ⟨hp₁, hp₂⟩ := Finset.mem_filter.mp hp
    exact ⟨(Finset.mem_filter.mp hp₁).2, hp₂⟩
  have hl (p : ℕ) (hp : p ∈ large) : MinusElevenInertPrime p ∧ P < p := by
    obtain ⟨hp₁, hp₂⟩ := Finset.mem_filter.mp hp
    exact ⟨(Finset.mem_filter.mp hp₁).2, Nat.lt_of_not_ge hp₂⟩
  have h₃E : 3 * E ≤ minusElevenSieveExpression n P E := by
    unfold minusElevenSieveExpression
    have h : E * 3 ≤ E * (3 + (5 / 2) * minusElevenInertCount P) := mul_le_mul_of_nonneg_left
      (le_add_of_nonneg_right (show (0 : ℝ) ≤ (5 / 2) * minusElevenInertCount P by positivity)) hE
    have htail := mul_nonneg (show 0 ≤ (8 / 3 : ℝ) * minusElevenDivisorMaximum n by positivity)
      (minusEleven_inert_tail_nonneg P)
    have hm := minusEleven_inert_mass_nonneg
    linarith only [hm, htail, h]
  have hscale : 0 < 1 - E := by linarith only [h₃E, hmargin]
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hG : 0 < G := lt_of_lt_of_le (by positivity) hB.shell_lower
  have hSmall : (1 - E) * (∑ p ∈ small, B p) ≤
      (minusElevenInertMass + (5 / 2) * E * minusElevenInertCount P) * G := by
    rw [Finset.mul_sum]
    calc
      _ ≤ ∑ p ∈ small, (minusElevenLocalMass p + (5 / 2) * E) * G :=
        Finset.sum_le_sum (fun p hp => hB.inert_small_upper p (hs p hp).1)
      _ = ((∑ p ∈ small, minusElevenLocalMass p) + (5 / 2) * E * small.card) * G := by
        rw [← Finset.sum_mul, Finset.sum_add_distrib, Finset.sum_const, nsmul_eq_mul]
        ring
      _ ≤ _ := by
        apply mul_le_mul_of_nonneg_right _ hG.le
        have hcard : (small.card : ℝ) ≤ minusElevenInertCount P := by
          exact_mod_cast minusEleven_inert_card_le hs
        exact add_le_add (minusEleven_inert_mass_sum_le (fun p hp => (hs p hp).1))
          (mul_le_mul_of_nonneg_left hcard (by positivity))
  have hLarge : (1 - E) * (∑ p ∈ large, B p) ≤
      (8 / 3) * D * minusElevenInertTail P * G := by
    have hpoint (p : ℕ) (hp : p ∈ large) :
        (1 - E) * B p ≤ (8 / 3) * D * (1 / (p : ℝ) ^ 2) * G := by
      calc
        _ ≤ (1 - E) * ((48 / 5) * D * n / (p : ℝ) ^ 2) :=
          mul_le_mul_of_nonneg_left (hB.inert_large_upper p (hl p hp).1) hscale.le
        _ = ((8 / 3) * D * (1 / (p : ℝ) ^ 2)) * ((18 / 5) * n * (1 - E)) := by ring
        _ ≤ _ := mul_le_mul_of_nonneg_left hB.shell_lower (by positivity)
    rw [Finset.mul_sum]
    calc
      _ ≤ ∑ p ∈ large, (8 / 3) * D * (1 / (p : ℝ) ^ 2) * G := Finset.sum_le_sum hpoint
      _ = (8 / 3) * D * (∑ p ∈ large, (1 : ℝ) / (p : ℝ) ^ 2) * G := by
        rw [← Finset.sum_mul, ← Finset.mul_sum]
      _ ≤ _ := mul_le_mul_of_nonneg_right
        (mul_le_mul_of_nonneg_left (minusEleven_inert_tail_sum_le hl) (by positivity)) hG.le
  have hweighted : (1 - E) * (B 11 + ∑ p ∈ s, B p) <
      G * (1 - E) := by
    calc
      _ = (1 - E) * B 11 + ((1 - E) * (∑ p ∈ small, B p) + (1 - E) * (∑ p ∈ large, B p)) := by
        rw [← Finset.sum_filter_add_sum_filter_not s (fun p => p ≤ P) B]
        ring
      _ ≤ (1 / 3 + 2 * E) * G +
          ((minusElevenInertMass + (5 / 2) * E * minusElevenInertCount P) * G +
            (8 / 3) * D * minusElevenInertTail P * G) :=
        add_le_add hB.ramified_upper (add_le_add hSmall hLarge)
      _ = (1 / 3 - E + minusElevenSieveExpression n P E) * G := by
        unfold minusElevenSieveExpression
        change _ = (1 / 3 - E + (minusElevenInertMass + (8 / 3) * D * minusElevenInertTail P +
          E * (3 + (5 / 2) * minusElevenInertCount P))) * G
        ring
      _ < _ := by
        rw [mul_comm G]
        apply mul_lt_mul_of_pos_right _ hG
        linarith only [hmargin]
  rw [mul_comm G] at hweighted
  have hcardR := (mul_lt_mul_iff_of_pos_left hscale).mp hweighted
  dsimp [B, G] at hcardR
  dsimp [G] at hG
  have hcard : (minusElevenBadShell n 11).card + ∑ p ∈ s, (minusElevenBadShell n p).card <
      (minusElevenOddShell n).card := by
    exact_mod_cast hcardR
  have hGQ : (0 : ℚ) < (minusElevenOddShell n).card := by exact_mod_cast hG
  obtain ⟨v, hv, h₁₁, hi⟩ := exists_unspoiled (candidates := minusElevenOddShell n)
    (exceptional := minusElevenBadShell n 11) (primes := s) (bad := minusElevenBadShell n)
    (Finset.card_pos.mp (by exact_mod_cast hG)) ((div_lt_iff₀ hGQ).mpr (by rw [one_mul]; exact_mod_cast hcard))
  exact good_shell_point hv h₁₁ (fun p hpn hp =>
    hi p (Finset.mem_filter.mpr ⟨Finset.mem_range.mpr (Nat.lt_succ_of_le hpn), hp⟩))

/-- The −11 row follows from the analytic shell estimates on the reduced
targets above the scan cutoff. The arithmetic margin and finite range
are discharged by the checked proofs. -/
theorem P11_sum_of_shell_bounds
    (hShell : ∀ t m : ℕ, 0 < t → Odd m → 2 ^ 25 < 2 ^ t * m → ¬ 11 ∣ 2 ^ t * m →
      IsReduced C11 (2 ^ t * m) → MinusElevenShellBounds (2 ^ t * m) (minusElevenFilteredError t m))
    (n : ℕ) :
    InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n := by
  apply P11_sum_classification restoration11
  intro v h₂ he h₁₁
  refine C11_core_coverage ?_ h₂ he h₁₁
  intro u hu₂ hue hu₁₁ hred
  by_cases hle : u ≤ 2 ^ 25
  · exact C11_finite_core_coverage hle hu₂ hue hu₁₁
  obtain ⟨t, m, hodd, rfl⟩ := Nat.exists_eq_two_pow_mul_odd (by omega : u ≠ 0)
  have ht : 0 < t := by
    by_contra! h
    obtain rfl : t = 0 := Nat.eq_zero_of_le_zero h
    obtain ⟨k, hk⟩ := hodd
    simp at hue
    omega
  have hm : 0 < m := by
    obtain ⟨k, hk⟩ := hodd
    omega
  have hI (p : ℕ) (hp : p ∈ m.primeFactors) : MinusElevenInertPrime p := by
    have hpp := Nat.prime_of_mem_primeFactors hp
    have hpd := Nat.dvd_of_mem_primeFactors hp
    have hpdu := hpd.trans (dvd_mul_left m (2 ^ t))
    rcases minusEleven_prime_cases hpp with rfl | rfl | hs | hi
    · exact (hodd.not_two_dvd_nat hpd).elim
    · exact (hu₁₁ hpdu).elim
    · exact (hred p hpp hpdu ((split_of_prime hpp).mp hs)).elim
    · exact hi
  obtain ⟨P, _, hmargin⟩ := minusEleven_sieve_margin_above_scan hm hI (Nat.le_of_lt (lt_of_not_ge hle))
  exact minusEleven_core_sum_of_shell_bounds (Nat.mul_pos (by positivity) hm)
    (minusEleven_filtered_error_nonneg t m) (hShell t m ht hodd (lt_of_not_ge hle) hu₁₁ hred) hmargin

end QuadraticNorms
