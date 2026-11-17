/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeMap
import Mathlib.Algebra.MvPolynomial.Degrees

/-!
# Effective degree bounds for the all-degree quadratic gauge

This module proves the paper's coordinate-degree estimate directly for the
single general `MvPolynomial` map from `GeneralGaugeMap.lean`.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

private theorem totalDegree_add_le_bound
    {p q : GaugePolynomial K} {d : ℕ}
    (hp : p.totalDegree ≤ d) (hq : q.totalDegree ≤ d) :
    (p + q).totalDegree ≤ d :=
  (MvPolynomial.totalDegree_add p q).trans (max_le hp hq)

private theorem totalDegree_sub_le_bound
    {p q : GaugePolynomial K} {d : ℕ}
    (hp : p.totalDegree ≤ d) (hq : q.totalDegree ≤ d) :
    (p - q).totalDegree ≤ d :=
  (MvPolynomial.totalDegree_sub p q).trans (max_le hp hq)

private theorem totalDegree_mul_le_bound
    {p q : GaugePolynomial K} {a b : ℕ}
    (hp : p.totalDegree ≤ a) (hq : q.totalDegree ≤ b) :
    (p * q).totalDegree ≤ a + b :=
  (MvPolynomial.totalDegree_mul p q).trans (Nat.add_le_add hp hq)

private theorem totalDegree_pow_le_bound
    {p : GaugePolynomial K} {a n : ℕ}
    (hp : p.totalDegree ≤ a) :
    (p ^ n).totalDegree ≤ n * a :=
  (MvPolynomial.totalDegree_pow p n).trans (Nat.mul_le_mul_left n hp)

private theorem totalDegree_C_mul_le_bound
    (c : K) {p : GaugePolynomial K} {d : ℕ}
    (hp : p.totalDegree ≤ d) :
    (MvPolynomial.C c * p).totalDegree ≤ d := by
  calc
    (MvPolynomial.C c * p).totalDegree ≤
        (MvPolynomial.C c : GaugePolynomial K).totalDegree + p.totalDegree :=
      MvPolynomial.totalDegree_mul _ _
    _ = p.totalDegree := by simp
    _ ≤ d := hp

/-- The recurrent polynomial `t = 1 + x*y` has total degree at most two. -/
theorem generalGaugeT_totalDegree :
    (generalGaugeT (K := K)).totalDegree ≤ 2 := by
  unfold generalGaugeT
  apply totalDegree_add_le_bound
  · simp
  · simpa using
      (totalDegree_mul_le_bound
        (p := (MvPolynomial.X 0 : GaugePolynomial K))
        (q := (MvPolynomial.X 1 : GaugePolynomial K))
        (a := 1) (b := 1) (by simp) (by simp))

/-- The recurrent polynomial `q` has total degree at most five. -/
theorem generalGaugeQ_totalDegree (G : K[X]) :
    (generalGaugeQ G).totalDegree ≤ 5 := by
  have ht : (generalGaugeT (K := K)).totalDegree ≤ 2 :=
    generalGaugeT_totalDegree
  have ht2 : ((generalGaugeT (K := K)) ^ 2).totalDegree ≤ 4 := by
    have h := totalDegree_pow_le_bound (p := generalGaugeT (K := K))
      (a := 2) (n := 2) ht
    norm_num at h ⊢
    exact h
  have hfirst :
      ((generalGaugeT (K := K)) ^ 2 * MvPolynomial.X 2).totalDegree ≤ 5 := by
    simpa using
      (totalDegree_mul_le_bound ht2
        (q := (MvPolynomial.X 2 : GaugePolynomial K))
        (a := 4) (b := 1) (by simp))
  have hx2 :
      ((MvPolynomial.X 1 : GaugePolynomial K) ^ 2).totalDegree ≤ 2 := by
    have h := totalDegree_pow_le_bound
      (p := (MvPolynomial.X 1 : GaugePolynomial K))
      (a := 1) (n := 2) (by simp)
    norm_num at h ⊢
  have hct :
      (MvPolynomial.C (3 : K) * generalGaugeT).totalDegree ≤ 2 :=
    totalDegree_C_mul_le_bound 3 ht
  have hlast :
      (1 + MvPolynomial.C (3 : K) * generalGaugeT).totalDegree ≤ 2 := by
    apply totalDegree_add_le_bound
    · simp
    · exact hct
  have hcx2 :
      (MvPolynomial.C (G.coeff 1 / G.coeff 3) *
        (MvPolynomial.X 1 : GaugePolynomial K) ^ 2).totalDegree ≤ 2 :=
    totalDegree_C_mul_le_bound _ hx2
  have hsecond :
      (MvPolynomial.C (G.coeff 1 / G.coeff 3) *
          (MvPolynomial.X 1 : GaugePolynomial K) ^ 2 *
          (1 + MvPolynomial.C 3 * generalGaugeT)).totalDegree ≤ 4 := by
    simpa using totalDegree_mul_le_bound hcx2 hlast
  unfold generalGaugeQ
  apply totalDegree_add_le_bound
  · exact hfirst
  · omega

/-- The first coordinate `Π=t*q` has total degree at most seven. -/
theorem generalGaugePi_totalDegree (G : K[X]) :
    (generalGaugePi G).totalDegree ≤ 7 := by
  unfold generalGaugePi
  exact totalDegree_mul_le_bound generalGaugeT_totalDegree
    (generalGaugeQ_totalDegree G)

private theorem generalGaugeB_tail_totalDegree
    (G : K[X]) (_hdeg : 3 ≤ G.natDegree) :
    (∑ k ∈ Finset.Icc 4 G.natDegree,
      MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
        generalGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
          generalGaugeQ G ^ k).totalDegree ≤ 6 * G.natDegree + 2 := by
  apply MvPolynomial.totalDegree_finsetSum_le
  intro k hk
  have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
  have hkN : k ≤ G.natDegree := (Finset.mem_Icc.mp hk).2
  have ht2 : ((generalGaugeT (K := K)) ^ 2).totalDegree ≤ 4 := by
    have h := totalDegree_pow_le_bound (p := generalGaugeT (K := K))
      (a := 2) (n := 2) generalGaugeT_totalDegree
    norm_num at h ⊢
    exact h
  have hx :
      ((MvPolynomial.X 0 : GaugePolynomial K) ^ (k - 2)).totalDegree ≤ k - 2 := by
    simpa using
      (totalDegree_pow_le_bound
        (p := (MvPolynomial.X 0 : GaugePolynomial K))
        (a := 1) (n := k - 2) (by simp))
  have hq : ((generalGaugeQ G) ^ k).totalDegree ≤ 5 * k := by
    have h := totalDegree_pow_le_bound (p := generalGaugeQ G)
      (a := 5) (n := k) (generalGaugeQ_totalDegree G)
    simpa [Nat.mul_comm] using h
  have h0 :
      (MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
        generalGaugeT ^ 2).totalDegree ≤ 4 :=
    totalDegree_C_mul_le_bound _ ht2
  have h1 :
      (MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
          generalGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2)).totalDegree ≤
        4 + (k - 2) :=
    totalDegree_mul_le_bound h0 hx
  have h2 :
      (MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
          generalGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
          generalGaugeQ G ^ k).totalDegree ≤
        (4 + (k - 2)) + 5 * k :=
    totalDegree_mul_le_bound h1 hq
  omega

/-- The second displayed coordinate has the effective degree bound `6N+2`. -/
theorem generalGaugeB_totalDegree (G : K[X]) (hdeg : 3 ≤ G.natDegree) :
    (generalGaugeB G).totalDegree ≤ 6 * G.natDegree + 2 := by
  have hq := generalGaugeQ_totalDegree G
  have htail := generalGaugeB_tail_totalDegree G hdeg
  have hlinear :
      (MvPolynomial.X 1 : GaugePolynomial K).totalDegree ≤
        6 * G.natDegree + 2 := by
    simp
  have hsecond :
      (MvPolynomial.C (3 * (G.coeff 3 / G.coeff 1)) *
          MvPolynomial.X 0 * generalGaugeQ G).totalDegree ≤
        6 * G.natDegree + 2 := by
    have hx :
        (MvPolynomial.C (3 * (G.coeff 3 / G.coeff 1)) *
          (MvPolynomial.X 0 : GaugePolynomial K)).totalDegree ≤ 1 :=
      totalDegree_C_mul_le_bound _ (by simp)
    have h := totalDegree_mul_le_bound hx hq
    omega
  have hthird :
      (MvPolynomial.C (2 * (G.coeff 2 / G.coeff 1)) * generalGaugeT *
          generalGaugeQ G).totalDegree ≤ 6 * G.natDegree + 2 := by
    have ht :
        (MvPolynomial.C (2 * (G.coeff 2 / G.coeff 1)) *
          (generalGaugeT (K := K))).totalDegree ≤ 2 :=
      totalDegree_C_mul_le_bound _ generalGaugeT_totalDegree
    have h := totalDegree_mul_le_bound ht hq
    omega
  unfold generalGaugeB
  apply totalDegree_add_le_bound
  · apply totalDegree_add_le_bound
    · apply totalDegree_add_le_bound
      · exact hlinear
      · exact hsecond
    · exact hthird
  · exact htail

private theorem generalGaugeC_tail_totalDegree
    (G : K[X]) :
    (∑ k ∈ Finset.Icc 4 G.natDegree,
      MvPolynomial.C (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
        (MvPolynomial.X 0 * generalGaugeQ G) ^ k).totalDegree ≤
      6 * G.natDegree := by
  apply MvPolynomial.totalDegree_finsetSum_le
  intro k hk
  have hkN : k ≤ G.natDegree := (Finset.mem_Icc.mp hk).2
  have hbase :
      ((MvPolynomial.X 0 : GaugePolynomial K) * generalGaugeQ G).totalDegree ≤ 6 :=
    totalDegree_mul_le_bound (by simp) (generalGaugeQ_totalDegree G)
  have hpow :
      (((MvPolynomial.X 0 : GaugePolynomial K) * generalGaugeQ G) ^ k).totalDegree ≤
        6 * k := by
    have h := totalDegree_pow_le_bound
      (p := (MvPolynomial.X 0 : GaugePolynomial K) * generalGaugeQ G)
      (a := 6) (n := k) hbase
    simpa [Nat.mul_comm] using h
  have h := totalDegree_C_mul_le_bound
    (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) hpow
  omega

/-- The third displayed coordinate has total degree at most `6N`. -/
theorem generalGaugeC_totalDegree (G : K[X]) (hdeg : 3 ≤ G.natDegree) :
    (generalGaugeC G).totalDegree ≤ 6 * G.natDegree := by
  have hct :
      (MvPolynomial.C (3 : K) * generalGaugeT).totalDegree ≤ 2 :=
    totalDegree_C_mul_le_bound 3 generalGaugeT_totalDegree
  have hparen :
      (MvPolynomial.C (5 : K) - MvPolynomial.C 3 * generalGaugeT).totalDegree ≤ 2 := by
    apply totalDegree_sub_le_bound
    · simp
    · exact hct
  have hfirst :
      (MvPolynomial.X 0 *
        (MvPolynomial.C (5 : K) - MvPolynomial.C 3 * generalGaugeT)).totalDegree ≤
        6 * G.natDegree := by
    have h := totalDegree_mul_le_bound
      (p := (MvPolynomial.X 0 : GaugePolynomial K)) (a := 1) (b := 2)
      (by simp) hparen
    omega
  have hx3 :
      ((MvPolynomial.X 0 : GaugePolynomial K) ^ 3).totalDegree ≤ 3 := by
    have h := totalDegree_pow_le_bound
      (p := (MvPolynomial.X 0 : GaugePolynomial K))
      (a := 1) (n := 3) (by simp)
    norm_num at h ⊢
  have hcx3 :
      (MvPolynomial.C (G.coeff 3 / G.coeff 1) *
        (MvPolynomial.X 0 : GaugePolynomial K) ^ 3).totalDegree ≤ 3 :=
    totalDegree_C_mul_le_bound _ hx3
  have hsecond :
      (MvPolynomial.C (G.coeff 3 / G.coeff 1) *
          (MvPolynomial.X 0 : GaugePolynomial K) ^ 3 *
          MvPolynomial.X 2).totalDegree ≤ 6 * G.natDegree := by
    have h := totalDegree_mul_le_bound hcx3
      (q := (MvPolynomial.X 2 : GaugePolynomial K))
      (a := 3) (b := 1) (by simp)
    omega
  have htail := generalGaugeC_tail_totalDegree G
  unfold generalGaugeC
  apply totalDegree_sub_le_bound
  · apply totalDegree_sub_le_bound
    · exact hfirst
    · exact hsecond
  · exact htail

/-- Every coordinate of the all-degree quadratic-gauge map satisfies the paper's
uniform `6N+2` bound. -/
theorem generalGaugeMap_totalDegree (G : K[X]) (hdeg : 3 ≤ G.natDegree)
    (i : Fin 3) :
    (generalGaugeMap G i).totalDegree ≤ 6 * G.natDegree + 2 := by
  fin_cases i
  · simpa [generalGaugeMap] using
      (show (generalGaugePi G).totalDegree ≤ 6 * G.natDegree + 2 by
        have h := generalGaugePi_totalDegree G
        omega)
  · simpa [generalGaugeMap] using generalGaugeB_totalDegree G hdeg
  · simpa [generalGaugeMap] using
      (show (generalGaugeC G).totalDegree ≤ 6 * G.natDegree + 2 by
        have h := generalGaugeC_totalDegree G hdeg
        omega)

/-- The target-preserving determinant-one normalization has the same effective
coordinate-degree bound. -/
theorem generalGaugeJacobianOneMap_totalDegree
    (G : K[X]) (hdeg : 3 ≤ G.natDegree) (i : Fin 3) :
    (generalGaugeJacobianOneMap G i).totalDegree ≤ 6 * G.natDegree + 2 := by
  fin_cases i
  · simpa [generalGaugeJacobianOneMap, scaleOutput] using
      totalDegree_C_mul_le_bound
        (K := K) (1 : K)
        (generalGaugeMap_totalDegree G hdeg (0 : Fin 3))
  · simpa [generalGaugeJacobianOneMap, scaleOutput] using
      totalDegree_C_mul_le_bound
        (K := K) (-1 / 2 : K)
        (generalGaugeMap_totalDegree G hdeg (1 : Fin 3))
  · simpa [generalGaugeJacobianOneMap, scaleOutput] using
      totalDegree_C_mul_le_bound
        (K := K) (1 : K)
        (generalGaugeMap_totalDegree G hdeg (2 : Fin 3))

#print axioms generalGaugeMap_totalDegree
#print axioms generalGaugeJacobianOneMap_totalDegree

end FiniteEtaleKeller