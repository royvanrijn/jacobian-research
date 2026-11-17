/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GaugeAssembly

/-!
# Assembly of the generic inverse-coordinate identity

This module proves, over an arbitrary commutative ring, the complete finite-sum
identity which rewrites `2*G_π/g₁ - B*S²` into the marked third-coordinate
expression.  It is the remaining coefficientwise bridge between the generic
inverse polynomial and the displayed all-degree map.
-/

noncomputable section

namespace FiniteEtaleKeller

variable {R : Type*} [CommRing R]

private theorem quadraticGauge_inverseCTail
    (S pi : R) (a : ℕ → R) (N : ℕ) :
    2 * (∑ k ∈ Finset.Icc 4 N, a k * pi ^ k * S ^ k) -
        (∑ k ∈ Finset.Icc 4 N, (k : R) * a k * pi ^ k * S ^ (k - 2)) * S ^ 2 =
      -∑ k ∈ Finset.Icc 4 N,
        ((k - 2 : ℕ) : R) * a k * pi ^ k * S ^ k := by
  rw [Finset.mul_sum, Finset.sum_mul, ← Finset.sum_sub_distrib]
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro k hk
  have hk2 : 2 ≤ k := by
    have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
    omega
  have hpow : S ^ (k - 2) * S ^ 2 = S ^ k := by
    rw [← pow_add]
    congr
    omega
  calc
    2 * (a k * pi ^ k * S ^ k) -
        (k : R) * a k * pi ^ k * S ^ (k - 2) * S ^ 2 =
      2 * (a k * pi ^ k * S ^ k) -
        ((k : R) * a k * pi ^ k) * (S ^ (k - 2) * S ^ 2) := by ring
    _ = -(((k - 2 : ℕ) : R) * a k * pi ^ k * S ^ k) := by
      rw [hpow, Nat.cast_sub hk2]
      ring

/-- Complete coefficientwise expansion of the marked third coordinate. -/
theorem quadraticGauge_inverseCExpansion
    (S Q pi c r : R) (a : ℕ → R) (N : ℕ) :
    2 * (S + c * pi * S ^ 2 + r * pi * S ^ 3 +
          ∑ k ∈ Finset.Icc 4 N, a k * pi ^ k * S ^ k) -
        (Q + 2 * c * pi + (3 * r - 1) * pi * S +
          ∑ k ∈ Finset.Icc 4 N,
            (k : R) * a k * pi ^ k * S ^ (k - 2)) * S ^ 2 =
      2 * S - Q * S ^ 2 + (1 - r) * pi * S ^ 3 -
        ∑ k ∈ Finset.Icc 4 N,
          ((k - 2 : ℕ) : R) * a k * pi ^ k * S ^ k := by
  have htail := quadraticGauge_inverseCTail S pi a N
  linear_combination htail

#print axioms quadraticGauge_inverseCExpansion

end FiniteEtaleKeller
