import QuadraticNorms.Definitions
import Mathlib.Algebra.QuadraticAlgebra.Basic
import Mathlib.Data.ZMod.Coprime
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination

namespace QuadraticNorms

open QuadraticAlgebra

private abbrev Order (b c : ℤ) := QuadraticAlgebra ℤ (-c) b

private theorem norm_eq_form {b c : ℤ} (z : Order b c) :
    z.norm = form b c z.re z.im := by
  simp [QuadraticAlgebra.norm_def, form]
  ring

private theorem four_norm {b c : ℤ} (z : Order b c) :
    4 * z.norm = (2 * z.re + b * z.im) ^ 2 + (4 * c - b ^ 2) * z.im ^ 2 := by
  rw [norm_eq_form]
  unfold form
  ring

private theorem norm_nonnegative {b c : ℤ} (hd : 0 < 4 * c - b ^ 2) (z : Order b c) :
    0 ≤ z.norm := by
  have h₁ := sq_nonneg (2 * z.re + b * z.im)
  have h₂ := mul_nonneg hd.le (sq_nonneg z.im)
  linarith [four_norm z]

private theorem norm_positive {b c : ℤ} (hd : 0 < 4 * c - b ^ 2) {z : Order b c}
    (hz : z ≠ 0) : 0 < z.norm := by
  have hnon := norm_nonnegative hd z
  by_contra! h
  have hn : z.norm = 0 := le_antisymm h hnon
  have heq := four_norm z
  have h₁ := sq_nonneg (2 * z.re + b * z.im)
  have h₂ := mul_nonneg hd.le (sq_nonneg z.im)
  have him : z.im = 0 := by
    by_contra hi
    have := mul_pos hd (sq_pos_of_ne_zero hi)
    linarith
  have hre : z.re = 0 := by
    simp [him, hn] at heq
    nlinarith [sq_nonneg z.re]
  exact hz (QuadraticAlgebra.ext (by simpa using hre) (by simpa using him))

private theorem exists_nearest (a m : ℤ) (hm : 0 < m) :
    ∃ q : ℤ, -m ≤ 2 * (a - m * q) ∧ 2 * (a - m * q) ≤ m := by
  let q := (a + m / 2) / m
  have heq : (a + m / 2) % m + m * q = a + m / 2 :=
    Int.emod_add_mul_ediv _ _
  have hlo := Int.emod_nonneg (a + m / 2) hm.ne'
  have hhi := Int.emod_lt_of_pos (a + m / 2) hm
  exact ⟨q, by omega, by omega⟩

/-- Integral Euclidean division, by rounding the two coordinates in
the completed-square expression for the norm. -/
private theorem euclidean_step {b c : ℤ}
    (hd : 0 < 4 * c - b ^ 2) (hd' : 4 * c - b ^ 2 ≤ 11)
    (a v : Order b c) (hv : v ≠ 0) :
    ∃ q : Order b c, (a - v * q).norm < v.norm := by
  let N := v.norm
  have hN : 0 < N := norm_positive hd hv
  let A := a * star v
  obtain ⟨j, hjlo, hjhi⟩ := exists_nearest A.im N hN
  let y := A.im - N * j
  have hylo : -N ≤ 2 * y := hjlo
  have hyhi : 2 * y ≤ N := hjhi
  obtain ⟨i, hilo, hihi⟩ := exists_nearest (2 * A.re + b * y) (2 * N) (by omega)
  let q : Order b c := ⟨i, j⟩
  let r : Order b c := A - (algebraMap ℤ (Order b c) N) * q
  let s : Order b c := a - v * q
  have hr : r = s * star v := by
    dsimp only [r, s, A, N]
    rw [algebraMap_norm_eq_mul_star, sub_mul]
    ac_rfl
  have hry : r.im = y := by simp [r, y, q]
  have hrx : r.re = A.re - N * i := by simp [r, q]
  have hrlo : -N ≤ 2 * r.re + b * r.im := by
    rw [hrx, hry]
    linarith
  have hrhi : 2 * r.re + b * r.im ≤ N := by
    rw [hrx, hry]
    linarith
  have hxsq : (2 * r.re + b * r.im) ^ 2 ≤ N ^ 2 :=
    sq_le_sq' hrlo hrhi
  have hysq : (2 * r.im) ^ 2 ≤ N ^ 2 :=
    sq_le_sq' (hry ▸ hylo) (hry ▸ hyhi)
  have hdysq : (4 * c - b ^ 2) * (2 * r.im) ^ 2 ≤ 11 * N ^ 2 :=
    calc
      _ ≤ 11 * (2 * r.im) ^ 2 := mul_le_mul_of_nonneg_right hd' (sq_nonneg _)
      _ ≤ _ := mul_le_mul_of_nonneg_left hysq (by decide)
  have hrs : r.norm = s.norm * N := by rw [hr, map_mul, norm_star]
  have hNsq : 0 < N ^ 2 := sq_pos_of_pos hN
  have hbound : 16 * (s.norm * N) < 16 * (N * N) := by
    rw [← hrs]
    nlinarith [four_norm r]
  exact ⟨q, lt_of_mul_lt_mul_right (by omega : s.norm * N < N * N) hN.le⟩

/-- A common divisor given by a Bézout combination in the quadratic order. -/
private theorem quadratic_bezout {b c : ℤ}
    (hd : 0 < 4 * c - b ^ 2) (hd' : 4 * c - b ^ 2 ≤ 11)
    (a v : Order b c) :
    ∃ g : Order b c, g ∣ a ∧ g ∣ v ∧ ∃ x y, g = a * x + v * y := by
  by_cases hv : v = 0
  · subst v
    exact ⟨a, dvd_rfl, dvd_zero _, 1, 0, by simp⟩
  obtain ⟨q, hq⟩ := euclidean_step hd hd' a v hv
  obtain ⟨g, hgv, hgr, x, y, hg⟩ := quadratic_bezout hd hd' v (a - v * q)
  refine ⟨g, ?_, hgv, y, x - q * y, ?_⟩
  · have hadd := dvd_add hgr (hgv.trans (dvd_mul_right v q))
    simpa using hadd
  · rw [hg]
    ring
termination_by v.norm.toNat
decreasing_by
  exact (Int.toNat_lt_toNat (norm_positive hd hv)).mpr hq

/-- The polynomial defining the quadratic order has a root modulo n. -/
def HasNormRoot (b c : ℤ) (n : ℕ) : Prop :=
  ∃ t : ZMod n, t ^ 2 - (b : ZMod n) * t + (c : ZMod n) = 0

private theorem primitive_of_root {b c : ℤ} {n : ℕ}
    (hd : 0 < 4 * c - b ^ 2) (hd' : 4 * c - b ^ 2 ≤ 11)
    (hn : 0 < n) (hroot : HasNormRoot b c n) : PrimitiveValues b c n := by
  letI : NeZero n := ⟨hn.ne'⟩
  obtain ⟨t, ht⟩ := hroot
  let f : Order b c →ₐ[ℤ] ZMod n :=
    QuadraticAlgebra.lift ⟨t, by
      simp only [zsmul_eq_mul, Int.cast_neg, mul_one]
      linear_combination ht⟩
  let w : Order b c := ⟨-(t.val : ℤ), 1⟩
  have hfN : f (n : Order b c) = 0 := by simp
  have hfw : f w = 0 := by
    simp [f, w, QuadraticAlgebra.lift, zsmul_eq_mul]
  obtain ⟨g, hgN, hgw, u, v, hg⟩ :=
    quadratic_bezout hd hd' (n : Order b c) w
  have hfg : f g = 0 := by rw [hg]; simp [hfN, hfw]
  have hdvd_norm : (n : ℤ) ∣ g.norm := by
    apply (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp
    rw [← map_intCast f]
    change f (algebraMap ℤ (Order b c) g.norm) = 0
    rw [algebraMap_norm_eq_mul_star, map_mul, hfg, zero_mul]
  have hprimitive : Int.gcd g.re g.im = 1 := by
    let d := Int.gcd g.re g.im
    have hgre := Int.gcd_dvd_left g.re g.im
    have hgim := Int.gcd_dvd_right g.re g.im
    obtain ⟨q, hq⟩ := hgw
    have hi : 1 = g.re * q.im + g.im * q.re + b * g.im * q.im := by
      simpa [w] using congrArg QuadraticAlgebra.im hq
    have hdvd : (d : ℤ) ∣ 1 := by
      rw [hi]
      exact dvd_add (dvd_add (hgre.trans (dvd_mul_right ..)) (hgim.trans (dvd_mul_right ..)))
        ((hgim.trans (dvd_mul_left _ _)).trans (dvd_mul_right ..))
    exact Nat.dvd_one.mp (by exact_mod_cast hdvd)
  have hnorm_dvd : g.norm ∣ (n : ℤ) := by
    obtain ⟨q, hq⟩ := hgN
    have hmul :
        (n : Order b c) * star g = (algebraMap ℤ (Order b c) g.norm) * q := by
      rw [algebraMap_norm_eq_mul_star, hq]
      ac_rfl
    have hre : (n : ℤ) * (g.re + b * g.im) = g.norm * q.re := by
      simpa using congrArg QuadraticAlgebra.re hmul
    have him : -(n : ℤ) * g.im = g.norm * q.im := by
      have h := congrArg QuadraticAlgebra.im hmul
      simpa [mul_neg, neg_mul] using h
    have hdi : g.norm ∣ (n : ℤ) * g.im := by
      use -q.im
      linear_combination -him
    have hdr : g.norm ∣ (n : ℤ) * g.re := by
      have heq : (n : ℤ) * g.re = (n : ℤ) * (g.re + b * g.im) - b * ((n : ℤ) * g.im) := by
        ring
      rw [heq]
      exact Int.dvd_sub ⟨_, hre⟩ ((hdi.trans (dvd_mul_left _ _)))
    have hbez := Int.gcd_eq_gcd_ab g.re g.im
    rw [hprimitive, Nat.cast_one] at hbez
    have hnEq : (n : ℤ) =
        ((n : ℤ) * g.re) * Int.gcdA g.re g.im + ((n : ℤ) * g.im) * Int.gcdB g.re g.im := by
      linear_combination (n : ℤ) * hbez
    rw [hnEq]
    exact dvd_add (hdr.trans (dvd_mul_right ..)) (hdi.trans (dvd_mul_right ..))
  have hnorm_nat : (g.norm.toNat : ℤ) = g.norm := Int.toNat_of_nonneg (norm_nonnegative hd g)
  refine ⟨hn, g.re, g.im, hprimitive, ?_⟩
  rw [← norm_eq_form, ← hnorm_nat]
  exact_mod_cast Nat.dvd_antisymm
    (Int.natCast_dvd_natCast.mp (by rwa [hnorm_nat]))
    (Int.natCast_dvd_natCast.mp (by rwa [hnorm_nat]))

/-- The residue obtained by dividing the negative first coordinate by
the second coordinate in a primitive norm representation. -/
def primitiveNormResidue (n : ℕ) (x y : ℤ) : ZMod n :=
  -(x : ZMod n) * (Int.gcdA y (n : ℤ) : ZMod n)

/-- The primitive-coordinate construction supplies both a polynomial
root and the linear congruence of the coordinates at that root. -/
theorem primitiveNormResidue_spec {b c x y : ℤ} {n : ℕ}
    (hxy : Int.gcd x y = 1) (hn : form b c x y = (n : ℤ)) :
    let t := primitiveNormResidue n x y
    t ^ 2 - (b : ZMod n) * t + (c : ZMod n) = 0 ∧ (x : ZMod n) + (y : ZMod n) * t = 0 := by
  have hgcd : Int.gcd y (n : ℤ) = 1 := by
    apply Nat.eq_one_iff_not_exists_prime_dvd.mpr
    intro p hp hd
    obtain ⟨hy, hpn⟩ := Int.dvd_gcd_iff.mp hd
    have hpx : (p : ℤ) ∣ x := by
      apply (Nat.prime_iff_prime_int.mp hp).dvd_of_dvd_pow (n := 2)
      have hxsq : x ^ 2 = (n : ℤ) - (b * x * y + c * y ^ 2) := by
        unfold form at hn
        linear_combination hn
      rw [hxsq]
      exact Int.dvd_sub hpn (dvd_add (hy.trans (dvd_mul_left _ _))
        ((hy.trans (dvd_pow_self _ (by decide : 2 ≠ 0))).trans (dvd_mul_left _ _)))
    have hpd := Int.dvd_gcd hpx hy
    rw [hxy] at hpd
    exact hp.ne_one (Nat.dvd_one.mp hpd)
  have hbez := Int.gcd_eq_gcd_ab y (n : ℤ)
  rw [hgcd, Nat.cast_one] at hbez
  let X : ZMod n := x
  let Y : ZMod n := y
  let U : ZMod n := Int.gcdA y (n : ℤ)
  have hYU : Y * U = 1 := by
    have h := congrArg (fun t : ℤ => (t : ZMod n)) hbez
    simpa [Y, U] using h.symm
  have hn' : X ^ 2 + (b : ZMod n) * X * Y + (c : ZMod n) * Y ^ 2 = 0 := by
    have h := congrArg (fun t : ℤ => (t : ZMod n)) hn
    simpa [X, Y, form] using h
  let t := primitiveNormResidue n x y
  have ht : Y * t = -X := by calc
    _ = -X * (Y * U) := by dsimp [t, primitiveNormResidue, X, Y, U]; ring
    _ = _ := by rw [hYU, mul_one]
  refine ⟨?_, by change X + Y * t = 0; rw [ht, add_neg_cancel]⟩
  change t ^ 2 - (b : ZMod n) * t + (c : ZMod n) = 0
  have heq : Y ^ 2 * (t ^ 2 - (b : ZMod n) * t + (c : ZMod n)) = 0 := by
    calc
      _ = (Y * t) ^ 2 - (b : ZMod n) * (Y * t) * Y + (c : ZMod n) * Y ^ 2 := by ring
      _ = _ := by rw [ht]; linear_combination hn'
  calc
    _ = (Y * U) ^ 2 * (t ^ 2 - (b : ZMod n) * t + (c : ZMod n)) := by rw [hYU]; simp
    _ = U ^ 2 * (Y ^ 2 * (t ^ 2 - (b : ZMod n) * t + (c : ZMod n))) := by ring
    _ = _ := by simp [heq]

private theorem root_of_primitive {b c : ℤ} {n : ℕ} (h : PrimitiveValues b c n) :
    HasNormRoot b c n := by
  obtain ⟨_, x, y, hxy, hn⟩ := h
  exact ⟨primitiveNormResidue n x y, (primitiveNormResidue_spec hxy hn).1⟩

/-- For the four discriminants in the companion paper, a primitive
positive coordinate value is equivalent to a root modulo that value. -/
theorem primitiveValues_iff_normRoot {b c : ℤ} {n : ℕ}
    (hd : 0 < 4 * c - b ^ 2) (hd' : 4 * c - b ^ 2 ≤ 11) :
    PrimitiveValues b c n ↔ 0 < n ∧ HasNormRoot b c n :=
  ⟨fun h => ⟨h.1, root_of_primitive h⟩, fun ⟨hn, h⟩ => primitive_of_root hd hd' hn h⟩

end QuadraticNorms
