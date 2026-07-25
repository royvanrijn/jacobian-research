/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GaugeAssembly
import Mathlib.Algebra.MvPolynomial.Degrees

/-!
# The all-degree quadratic-gauge polynomial map

This module packages the displayed construction as one genuine three-variable
`MvPolynomial` map for an arbitrary univariate seed `G`.  It also proves the
uniform coordinate-degree estimate from the paper:

`maxᵢ totalDegree(Fᵢ) ≤ 6 * natDegree(G) + 2`.

The proof is termwise and uses the complete finite coefficient sums, rather
than expansion in fixed degrees.  The general Jacobian theorem is intentionally
kept separate from this degree certificate.
-/

noncomputable section

open MvPolynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

abbrev GaugePolynomial (K : Type*) [CommSemiring K] := MvPolynomial (Fin 3) K

private theorem totalDegree_add_le_of_le
    {p q : GaugePolynomial K} {d : ℕ}
    (hp : p.totalDegree ≤ d) (hq : q.totalDegree ≤ d) :
    (p + q).totalDegree ≤ d :=
  (MvPolynomial.totalDegree_add p q).trans (max_le hp hq)

private theorem totalDegree_mul_le_of_le
    {p q : GaugePolynomial K} {dp dq : ℕ}
    (hp : p.totalDegree ≤ dp) (hq : q.totalDegree ≤ dq) :
    (p * q).totalDegree ≤ dp + dq :=
  (MvPolynomial.totalDegree_mul p q).trans (Nat.add_le_add hp hq)

private theorem totalDegree_pow_le_of_le
    {p : GaugePolynomial K} {d n : ℕ} (hp : p.totalDegree ≤ d) :
    (p ^ n).totalDegree ≤ n * d :=
  (MvPolynomial.totalDegree_pow p n).trans (Nat.mul_le_mul_left n hp)

private theorem totalDegree_C_mul_le_of_le
    (c : K) {p : GaugePolynomial K} {d : ℕ} (hp : p.totalDegree ≤ d) :
    (MvPolynomial.C c * p).totalDegree ≤ d := by
  simpa using
    (totalDegree_mul_le_of_le (p := MvPolynomial.C c) (q := p)
      (dp := 0) (dq := d) (by simp) hp)

/-- The source polynomial `t = 1 + xy`. -/
def quadraticGaugeT : GaugePolynomial K :=
  1 + MvPolynomial.X 0 * MvPolynomial.X 1

/-- The recurrent source polynomial of the all-degree gauge. -/
def quadraticGaugeQ (G : Polynomial K) : GaugePolynomial K :=
  quadraticGaugeT ^ 2 * MvPolynomial.X 2
    + MvPolynomial.C (G.coeff 1 / G.coeff 3) * MvPolynomial.X 1 ^ 2 *
      (1 + MvPolynomial.C 3 * quadraticGaugeT)

/-- The first coordinate `Pi = t*q`. -/
def quadraticGaugePi (G : Polynomial K) : GaugePolynomial K :=
  quadraticGaugeT * quadraticGaugeQ G

/-- The second displayed coordinate. -/
def quadraticGaugeB (G : Polynomial K) : GaugePolynomial K :=
  MvPolynomial.X 1
    + MvPolynomial.C (3 * G.coeff 3 / G.coeff 1) *
        (MvPolynomial.X 0 * quadraticGaugeQ G)
    + MvPolynomial.C (2 * G.coeff 2 / G.coeff 1) *
        (quadraticGaugeT * quadraticGaugeQ G)
    + ∑ k ∈ Finset.Icc 4 G.natDegree,
        MvPolynomial.C ((k : K) * G.coeff k / G.coeff 1) *
          (quadraticGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
            quadraticGaugeQ G ^ k)

/-- The third displayed coordinate.  Negative signs are absorbed into the
coefficients so that the degree proof uses only addition and multiplication. -/
def quadraticGaugeC (G : Polynomial K) : GaugePolynomial K :=
  MvPolynomial.X 0 *
      (MvPolynomial.C 5 + MvPolynomial.C (-3 : K) * quadraticGaugeT)
    + MvPolynomial.C (-(G.coeff 3 / G.coeff 1)) *
        (MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2)
    + ∑ k ∈ Finset.Icc 4 G.natDegree,
        MvPolynomial.C (-((k - 2 : ℕ) : K) * G.coeff k / G.coeff 1) *
          ((MvPolynomial.X 0 * quadraticGaugeQ G) ^ k)

/-- The displayed determinant-`-2` quadratic-gauge map before output
normalization. -/
def quadraticGaugePolynomialMap (G : Polynomial K) :
    Fin 3 → GaugePolynomial K :=
  ![quadraticGaugePi G, quadraticGaugeB G, quadraticGaugeC G]

/-- The source polynomial `t` has total degree at most two. -/
theorem quadraticGaugeT_totalDegree_le :
    (quadraticGaugeT : GaugePolynomial K).totalDegree ≤ 2 := by
  unfold quadraticGaugeT
  apply totalDegree_add_le_of_le
  · simp
  · simpa using
      (totalDegree_mul_le_of_le
        (p := MvPolynomial.X 0 : GaugePolynomial K)
        (q := MvPolynomial.X 1 : GaugePolynomial K)
        (dp := 1) (dq := 1) (by simp) (by simp))

/-- The recurrent polynomial `q` has total degree at most five. -/
theorem quadraticGaugeQ_totalDegree_le (G : Polynomial K) :
    (quadraticGaugeQ G).totalDegree ≤ 5 := by
  have ht : (quadraticGaugeT : GaugePolynomial K).totalDegree ≤ 2 :=
    quadraticGaugeT_totalDegree_le
  have ht2 : (quadraticGaugeT ^ 2 : GaugePolynomial K).totalDegree ≤ 4 := by
    simpa using
      (totalDegree_pow_le_of_le (p := quadraticGaugeT : GaugePolynomial K)
        (d := 2) (n := 2) ht)
  have hz : (MvPolynomial.X 2 : GaugePolynomial K).totalDegree ≤ 1 := by simp
  have hfirst :
      (quadraticGaugeT ^ 2 * MvPolynomial.X 2 : GaugePolynomial K).totalDegree ≤ 5 := by
    simpa using totalDegree_mul_le_of_le ht2 hz
  have hy2 :
      ((MvPolynomial.X 1 : GaugePolynomial K) ^ 2).totalDegree ≤ 2 := by
    simpa using
      (totalDegree_pow_le_of_le
        (p := MvPolynomial.X 1 : GaugePolynomial K)
        (d := 1) (n := 2) (by simp))
  have hthreeT :
      (MvPolynomial.C (3 : K) * quadraticGaugeT).totalDegree ≤ 2 :=
    totalDegree_C_mul_le_of_le 3 ht
  have hbracket :
      (1 + MvPolynomial.C (3 : K) * quadraticGaugeT : GaugePolynomial K).totalDegree ≤ 2 := by
    apply totalDegree_add_le_of_le
    · simp
    · exact hthreeT
  have hcoeffY2 :
      (MvPolynomial.C (G.coeff 1 / G.coeff 3) *
          (MvPolynomial.X 1 : GaugePolynomial K) ^ 2).totalDegree ≤ 2 :=
    totalDegree_C_mul_le_of_le _ hy2
  have hsecond :
      (MvPolynomial.C (G.coeff 1 / G.coeff 3) * MvPolynomial.X 1 ^ 2 *
          (1 + MvPolynomial.C 3 * quadraticGaugeT) : GaugePolynomial K).totalDegree ≤ 4 := by
    simpa [mul_assoc] using totalDegree_mul_le_of_le hcoeffY2 hbracket
  unfold quadraticGaugeQ
  exact totalDegree_add_le_of_le hfirst (hsecond.trans (by omega))

/-- The first coordinate has total degree at most seven. -/
theorem quadraticGaugePi_totalDegree_le (G : Polynomial K) :
    (quadraticGaugePi G).totalDegree ≤ 7 := by
  unfold quadraticGaugePi
  simpa using totalDegree_mul_le_of_le
    (quadraticGaugeT_totalDegree_le (K := K))
    (quadraticGaugeQ_totalDegree_le G)

/-- The second coordinate satisfies the paper's uniform `6N+2` bound. -/
theorem quadraticGaugeB_totalDegree_le (G : Polynomial K)
    (hdeg : 3 ≤ G.natDegree) :
    (quadraticGaugeB G).totalDegree ≤ 6 * G.natDegree + 2 := by
  let N := G.natDegree
  have hx : (MvPolynomial.X 0 : GaugePolynomial K).totalDegree ≤ 1 := by simp
  have hy : (MvPolynomial.X 1 : GaugePolynomial K).totalDegree ≤ 1 := by simp
  have ht : (quadraticGaugeT : GaugePolynomial K).totalDegree ≤ 2 :=
    quadraticGaugeT_totalDegree_le
  have hq : (quadraticGaugeQ G).totalDegree ≤ 5 := quadraticGaugeQ_totalDegree_le G
  have hxq :
      (MvPolynomial.X 0 * quadraticGaugeQ G : GaugePolynomial K).totalDegree ≤ 6 := by
    simpa using totalDegree_mul_le_of_le hx hq
  have htq :
      (quadraticGaugeT * quadraticGaugeQ G : GaugePolynomial K).totalDegree ≤ 7 := by
    simpa using totalDegree_mul_le_of_le ht hq
  have htermXQ :
      (MvPolynomial.C (3 * G.coeff 3 / G.coeff 1) *
          (MvPolynomial.X 0 * quadraticGaugeQ G) : GaugePolynomial K).totalDegree ≤ 6 :=
    totalDegree_C_mul_le_of_le _ hxq
  have htermTQ :
      (MvPolynomial.C (2 * G.coeff 2 / G.coeff 1) *
          (quadraticGaugeT * quadraticGaugeQ G) : GaugePolynomial K).totalDegree ≤ 7 :=
    totalDegree_C_mul_le_of_le _ htq
  have hbase12 :
      (MvPolynomial.X 1 +
          MvPolynomial.C (3 * G.coeff 3 / G.coeff 1) *
            (MvPolynomial.X 0 * quadraticGaugeQ G) : GaugePolynomial K).totalDegree ≤ 7 :=
    totalDegree_add_le_of_le (hy.trans (by omega)) (htermXQ.trans (by omega))
  have hbase :
      (MvPolynomial.X 1 +
          MvPolynomial.C (3 * G.coeff 3 / G.coeff 1) *
            (MvPolynomial.X 0 * quadraticGaugeQ G)
          + MvPolynomial.C (2 * G.coeff 2 / G.coeff 1) *
            (quadraticGaugeT * quadraticGaugeQ G) : GaugePolynomial K).totalDegree ≤ 7 :=
    totalDegree_add_le_of_le hbase12 htermTQ
  have ht2 : (quadraticGaugeT ^ 2 : GaugePolynomial K).totalDegree ≤ 4 := by
    simpa using
      (totalDegree_pow_le_of_le (p := quadraticGaugeT : GaugePolynomial K)
        (d := 2) (n := 2) ht)
  have hsum :
      (∑ k ∈ Finset.Icc 4 N,
          MvPolynomial.C ((k : K) * G.coeff k / G.coeff 1) *
            (quadraticGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
              quadraticGaugeQ G ^ k) : GaugePolynomial K).totalDegree ≤ 6 * N + 2 := by
    apply MvPolynomial.totalDegree_finsetSum_le
    intro k hk
    have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
    have hkN : k ≤ N := (Finset.mem_Icc.mp hk).2
    have hxpow :
        ((MvPolynomial.X 0 : GaugePolynomial K) ^ (k - 2)).totalDegree ≤ k - 2 := by
      simpa using
        (totalDegree_pow_le_of_le
          (p := MvPolynomial.X 0 : GaugePolynomial K)
          (d := 1) (n := k - 2) hx)
    have hqpow :
        (quadraticGaugeQ G ^ k).totalDegree ≤ k * 5 :=
      totalDegree_pow_le_of_le (p := quadraticGaugeQ G) (d := 5) (n := k) hq
    have hfirstProduct :
        (quadraticGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) : GaugePolynomial K).totalDegree ≤
          4 + (k - 2) :=
      totalDegree_mul_le_of_le ht2 hxpow
    have hproduct :
        (quadraticGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
            quadraticGaugeQ G ^ k : GaugePolynomial K).totalDegree ≤
          (4 + (k - 2)) + k * 5 :=
      totalDegree_mul_le_of_le hfirstProduct hqpow
    have hcoefficient :=
      totalDegree_C_mul_le_of_le
        ((k : K) * G.coeff k / G.coeff 1) hproduct
    exact hcoefficient.trans (by omega)
  unfold quadraticGaugeB
  apply totalDegree_add_le_of_le
  · exact hbase.trans (by omega)
  · simpa [N] using hsum

/-- The third coordinate satisfies the paper's uniform `6N+2` bound. -/
theorem quadraticGaugeC_totalDegree_le (G : Polynomial K)
    (hdeg : 3 ≤ G.natDegree) :
    (quadraticGaugeC G).totalDegree ≤ 6 * G.natDegree + 2 := by
  let N := G.natDegree
  have hx : (MvPolynomial.X 0 : GaugePolynomial K).totalDegree ≤ 1 := by simp
  have hz : (MvPolynomial.X 2 : GaugePolynomial K).totalDegree ≤ 1 := by simp
  have ht : (quadraticGaugeT : GaugePolynomial K).totalDegree ≤ 2 :=
    quadraticGaugeT_totalDegree_le
  have hq : (quadraticGaugeQ G).totalDegree ≤ 5 := quadraticGaugeQ_totalDegree_le G
  have hminusThreeT :
      (MvPolynomial.C (-3 : K) * quadraticGaugeT).totalDegree ≤ 2 :=
    totalDegree_C_mul_le_of_le (-3) ht
  have hbracket :
      (MvPolynomial.C (5 : K) + MvPolynomial.C (-3 : K) * quadraticGaugeT).totalDegree ≤ 2 := by
    apply totalDegree_add_le_of_le
    · simp
    · exact hminusThreeT
  have hfirst :
      (MvPolynomial.X 0 *
          (MvPolynomial.C 5 + MvPolynomial.C (-3 : K) * quadraticGaugeT) :
        GaugePolynomial K).totalDegree ≤ 3 := by
    simpa using totalDegree_mul_le_of_le hx hbracket
  have hx3 :
      ((MvPolynomial.X 0 : GaugePolynomial K) ^ 3).totalDegree ≤ 3 := by
    simpa using
      (totalDegree_pow_le_of_le
        (p := MvPolynomial.X 0 : GaugePolynomial K)
        (d := 1) (n := 3) hx)
  have hx3z :
      ((MvPolynomial.X 0 : GaugePolynomial K) ^ 3 * MvPolynomial.X 2).totalDegree ≤ 4 := by
    simpa using totalDegree_mul_le_of_le hx3 hz
  have hsecond :
      (MvPolynomial.C (-(G.coeff 3 / G.coeff 1)) *
          (MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2) : GaugePolynomial K).totalDegree ≤ 4 :=
    totalDegree_C_mul_le_of_le _ hx3z
  have hbase :
      (MvPolynomial.X 0 *
          (MvPolynomial.C 5 + MvPolynomial.C (-3 : K) * quadraticGaugeT)
        + MvPolynomial.C (-(G.coeff 3 / G.coeff 1)) *
          (MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2) : GaugePolynomial K).totalDegree ≤ 4 :=
    totalDegree_add_le_of_le (hfirst.trans (by omega)) hsecond
  have hxq :
      (MvPolynomial.X 0 * quadraticGaugeQ G : GaugePolynomial K).totalDegree ≤ 6 := by
    simpa using totalDegree_mul_le_of_le hx hq
  have hsum :
      (∑ k ∈ Finset.Icc 4 N,
          MvPolynomial.C (-((k - 2 : ℕ) : K) * G.coeff k / G.coeff 1) *
            ((MvPolynomial.X 0 * quadraticGaugeQ G) ^ k) : GaugePolynomial K).totalDegree ≤
        6 * N + 2 := by
    apply MvPolynomial.totalDegree_finsetSum_le
    intro k hk
    have hkN : k ≤ N := (Finset.mem_Icc.mp hk).2
    have hpow :
        ((MvPolynomial.X 0 * quadraticGaugeQ G : GaugePolynomial K) ^ k).totalDegree ≤ k * 6 :=
      totalDegree_pow_le_of_le
        (p := MvPolynomial.X 0 * quadraticGaugeQ G)
        (d := 6) (n := k) hxq
    have hcoefficient :=
      totalDegree_C_mul_le_of_le
        (-((k - 2 : ℕ) : K) * G.coeff k / G.coeff 1) hpow
    exact hcoefficient.trans (by omega)
  unfold quadraticGaugeC
  apply totalDegree_add_le_of_le
  · exact hbase.trans (by omega)
  · simpa [N] using hsum

/-- Every coordinate of the all-degree map has total degree at most `6N+2`. -/
theorem quadraticGaugePolynomialMap_totalDegree_le (G : Polynomial K)
    (hdeg : 3 ≤ G.natDegree) (i : Fin 3) :
    (quadraticGaugePolynomialMap G i).totalDegree ≤ 6 * G.natDegree + 2 := by
  fin_cases i
  · simpa [quadraticGaugePolynomialMap] using
      (quadraticGaugePi_totalDegree_le G).trans (by omega)
  · simpa [quadraticGaugePolynomialMap] using quadraticGaugeB_totalDegree_le G hdeg
  · simpa [quadraticGaugePolynomialMap] using quadraticGaugeC_totalDegree_le G hdeg

#print axioms quadraticGaugePolynomialMap_totalDegree_le

end FiniteEtaleKeller
