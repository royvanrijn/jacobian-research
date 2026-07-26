/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeMap

/-!
# Jacobian determinant of the all-degree quadratic gauge

This module proves the constant determinant structurally.  The arbitrary
high-degree part is packaged as a polynomial tail in `u = x*q`; its first two
formal derivatives account for every finite high-degree sum, while a universal
low-degree calculation proves that the determinant is independent of the tail.
-/

noncomputable section

open Matrix Function
open MvPolynomial
open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

private def quadraticGaugeQ (a : K) : GaugePolynomial K :=
  generalGaugeT ^ 2 * MvPolynomial.X 2 +
    MvPolynomial.C a * MvPolynomial.X 1 ^ 2 *
      (1 + MvPolynomial.C 3 * generalGaugeT)

private def quadraticGaugePi (a : K) : GaugePolynomial K :=
  generalGaugeT * quadraticGaugeQ a

private def quadraticGaugeU (a : K) : GaugePolynomial K :=
  MvPolynomial.X 0 * quadraticGaugeQ a

private def quadraticGaugeTail
    (d : ℕ → K) (N : ℕ) (u : GaugePolynomial K) : GaugePolynomial K :=
  ∑ k ∈ Finset.Icc 4 N,
    MvPolynomial.C (d k) * u ^ (k - 2)

private def quadraticGaugeTailDeriv
    (d : ℕ → K) (N : ℕ) (u : GaugePolynomial K) : GaugePolynomial K :=
  ∑ k ∈ Finset.Icc 4 N,
    MvPolynomial.C (((k - 2 : ℕ) : K) * d k) * u ^ (k - 3)

private def quadraticGaugeTailSecond
    (d : ℕ → K) (N : ℕ) (u : GaugePolynomial K) : GaugePolynomial K :=
  ∑ k ∈ Finset.Icc 4 N,
    MvPolynomial.C
      (((k - 3 : ℕ) : K) * (((k - 2 : ℕ) : K) * d k)) *
        u ^ (k - 4)

private theorem pderiv_quadraticGaugeTail
    (d : ℕ → K) (N : ℕ) (u : GaugePolynomial K) (i : Fin 3) :
    pderiv i (quadraticGaugeTail d N u) =
      quadraticGaugeTailDeriv d N u * pderiv i u := by
  rw [quadraticGaugeTail, quadraticGaugeTailDeriv, map_sum, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro k hk
  have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
  have hsub : k - 2 - 1 = k - 3 := by omega
  simp only [pderiv_mul, pderiv_C, zero_mul, zero_add, pderiv_pow]
  rw [hsub, MvPolynomial.C_mul, MvPolynomial.C_eq_coe_nat]
  ring

private theorem pderiv_quadraticGaugeTailDeriv
    (d : ℕ → K) (N : ℕ) (u : GaugePolynomial K) (i : Fin 3) :
    pderiv i (quadraticGaugeTailDeriv d N u) =
      quadraticGaugeTailSecond d N u * pderiv i u := by
  rw [quadraticGaugeTailDeriv, quadraticGaugeTailSecond, map_sum, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro k hk
  have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
  have hsub : k - 3 - 1 = k - 4 := by omega
  simp only [pderiv_mul, pderiv_C, zero_mul, zero_add, pderiv_pow]
  rw [hsub]
  simp only [MvPolynomial.C_mul, MvPolynomial.C_eq_coe_nat]
  ring

private theorem quadraticGaugeTail_B
    (d : ℕ → K) (N : ℕ)
    (t x q : GaugePolynomial K) :
    (t * q) ^ 2 *
        (MvPolynomial.C 2 * quadraticGaugeTail d N (x * q) +
          (x * q) * quadraticGaugeTailDeriv d N (x * q)) =
      ∑ k ∈ Finset.Icc 4 N,
        MvPolynomial.C ((k : K) * d k) * t ^ 2 * x ^ (k - 2) * q ^ k := by
  rw [quadraticGaugeTail, quadraticGaugeTailDeriv]
  simp only [Finset.mul_sum, mul_add]
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro k hk
  have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
  have hpowU :
      (x * q) * (x * q) ^ (k - 3) = (x * q) ^ (k - 2) := by
    calc
      (x * q) * (x * q) ^ (k - 3) =
          (x * q) ^ (1 + (k - 3)) := by rw [pow_add]; simp
      _ = (x * q) ^ (k - 2) := by congr; omega
  have hnat : 2 + (k - 2) = k := by omega
  have hcast : (2 : K) + ((k - 2 : ℕ) : K) = (k : K) := by
    calc
      (2 : K) + ((k - 2 : ℕ) : K) =
          ((2 + (k - 2) : ℕ) : K) := (Nat.cast_add 2 (k - 2)).symm
      _ = (k : K) := congrArg (fun n : ℕ => (n : K)) hnat
  have hcoeff :
      (2 : K) * d k + ((k - 2 : ℕ) : K) * d k = (k : K) * d k := by
    rw [← add_mul, hcast]
  have hcoeffPoly :
      (MvPolynomial.C 2 : GaugePolynomial K) * MvPolynomial.C (d k) +
          MvPolynomial.C (((k - 2 : ℕ) : K) * d k) =
        MvPolynomial.C ((k : K) * d k) := by
    rw [← MvPolynomial.C_mul, ← MvPolynomial.C_add]
    exact congrArg (fun z : K => (MvPolynomial.C z : GaugePolynomial K)) hcoeff
  have hpowUC :
      (x * q) *
          (MvPolynomial.C (((k - 2 : ℕ) : K) * d k) *
            (x * q) ^ (k - 3)) =
        MvPolynomial.C (((k - 2 : ℕ) : K) * d k) *
          (x * q) ^ (k - 2) := by
    calc
      (x * q) *
          (MvPolynomial.C (((k - 2 : ℕ) : K) * d k) *
            (x * q) ^ (k - 3)) =
        MvPolynomial.C (((k - 2 : ℕ) : K) * d k) *
          ((x * q) * (x * q) ^ (k - 3)) := by ring
      _ = MvPolynomial.C (((k - 2 : ℕ) : K) * d k) *
          (x * q) ^ (k - 2) := by rw [hpowU]
  have hPpow :
      (t * q) ^ 2 * (x * q) ^ (k - 2) =
        t ^ 2 * x ^ (k - 2) * q ^ k := by
    calc
      (t * q) ^ 2 * (x * q) ^ (k - 2) =
          (t ^ 2 * q ^ 2) * (x ^ (k - 2) * q ^ (k - 2)) := by
            rw [mul_pow, mul_pow]
      _ = t ^ 2 * x ^ (k - 2) * (q ^ 2 * q ^ (k - 2)) := by ring
      _ = t ^ 2 * x ^ (k - 2) * q ^ (2 + (k - 2)) := by
        rw [pow_add]
      _ = t ^ 2 * x ^ (k - 2) * q ^ k := by rw [hnat]
  calc
    (t * q) ^ 2 *
          (MvPolynomial.C 2 *
              (MvPolynomial.C (d k) * (x * q) ^ (k - 2))) +
        (t * q) ^ 2 *
          ((x * q) *
            (MvPolynomial.C (((k - 2 : ℕ) : K) * d k) *
              (x * q) ^ (k - 3))) =
      (t * q) ^ 2 *
        (((MvPolynomial.C 2 : GaugePolynomial K) * MvPolynomial.C (d k) +
            MvPolynomial.C (((k - 2 : ℕ) : K) * d k)) *
          (x * q) ^ (k - 2)) := by
            rw [hpowUC]
            ring
    _ = (t * q) ^ 2 *
        (MvPolynomial.C ((k : K) * d k) * (x * q) ^ (k - 2)) := by
          rw [hcoeffPoly]
    _ = MvPolynomial.C ((k : K) * d k) *
        t ^ 2 * x ^ (k - 2) * q ^ k := by
          calc
            (t * q) ^ 2 *
                (MvPolynomial.C ((k : K) * d k) * (x * q) ^ (k - 2)) =
              MvPolynomial.C ((k : K) * d k) *
                ((t * q) ^ 2 * (x * q) ^ (k - 2)) := by ring
            _ = MvPolynomial.C ((k : K) * d k) *
                (t ^ 2 * x ^ (k - 2) * q ^ k) := by rw [hPpow]
            _ = MvPolynomial.C ((k : K) * d k) *
                t ^ 2 * x ^ (k - 2) * q ^ k := by ring

private theorem quadraticGaugeTail_C
    (d : ℕ → K) (N : ℕ) (u : GaugePolynomial K) :
    u ^ 3 * quadraticGaugeTailDeriv d N u =
      ∑ k ∈ Finset.Icc 4 N,
        MvPolynomial.C (((k - 2 : ℕ) : K) * d k) * u ^ k := by
  rw [quadraticGaugeTailDeriv, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro k hk
  have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
  have hnat : 3 + (k - 3) = k := by omega
  calc
    u ^ 3 *
        (MvPolynomial.C (((k - 2 : ℕ) : K) * d k) * u ^ (k - 3)) =
      MvPolynomial.C (((k - 2 : ℕ) : K) * d k) *
        u ^ (3 + (k - 3)) := by
          rw [pow_add]
          ring
    _ = MvPolynomial.C (((k - 2 : ℕ) : K) * d k) * u ^ k := by
      rw [hnat]

private def quadraticGaugeBaseB (a c : K) : GaugePolynomial K :=
  MvPolynomial.X 1 +
    MvPolynomial.C (3 * a⁻¹) * quadraticGaugeU a +
    MvPolynomial.C (2 * c) * quadraticGaugePi a

private def quadraticGaugeBaseC (a : K) : GaugePolynomial K :=
  MvPolynomial.X 0 * (MvPolynomial.C 5 - MvPolynomial.C 3 * generalGaugeT) -
    MvPolynomial.C a⁻¹ * MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2

private def quadraticGaugeTailB
    (a : K) (R Rp : GaugePolynomial K) : GaugePolynomial K :=
  quadraticGaugePi a ^ 2 *
    (MvPolynomial.C 2 * R + quadraticGaugeU a * Rp)

private def quadraticGaugeTailC
    (a : K) (Rp : GaugePolynomial K) : GaugePolynomial K :=
  -(quadraticGaugeU a ^ 3 * Rp)

private def quadraticGaugeBaseMap
    (a c : K) : Fin 3 → GaugePolynomial K :=
  ![quadraticGaugePi a, quadraticGaugeBaseB a c, quadraticGaugeBaseC a]

private def quadraticGaugeWithTail
    (a c : K) (R Rp : GaugePolynomial K) : Fin 3 → GaugePolynomial K :=
  ![
    quadraticGaugePi a,
    quadraticGaugeBaseB a c + quadraticGaugeTailB a R Rp,
    quadraticGaugeBaseC a + quadraticGaugeTailC a Rp]

set_option maxHeartbeats 0 in
set_option maxRecDepth 100000 in
private theorem jacobianDet_quadraticGaugeBase
    (a c : K) (ha : a ≠ 0) :
    jacobianDet (quadraticGaugeBaseMap a c) = MvPolynomial.C (-2) := by
  classical
  have hunit :
      (MvPolynomial.C a : GaugePolynomial K) * MvPolynomial.C a⁻¹ = 1 := by
    rw [← MvPolynomial.C_mul]
    simp [ha]
  simp only [jacobianDet, jacobianMatrix, det_fin_three, of_apply,
    quadraticGaugeBaseMap, quadraticGaugeBaseB, quadraticGaugeBaseC,
    quadraticGaugePi, quadraticGaugeU, quadraticGaugeQ, generalGaugeT,
    cons_val_zero, cons_val_one, cons_val_two, head_cons, tail_cons,
    map_add, map_sub, Derivation.map_one_eq_zero,
    pderiv_mul, pderiv_pow, pderiv_C, pderiv_X_self, pderiv_X_of_ne,
    ne_eq, Fin.reduceEq, not_false_eq_true]
  simp only [map_neg, map_ofNat, MvPolynomial.C_mul,
    MvPolynomial.C_eq_coe_nat]
  grobner

set_option maxHeartbeats 0 in
set_option maxRecDepth 100000 in
private theorem quadraticGauge_crossDet
    (a c : K) (ha : a ≠ 0) :
    quadraticGaugePi a ^ 2 *
        jacobianDet ![
          quadraticGaugePi a, quadraticGaugeU a, quadraticGaugeBaseC a] =
      quadraticGaugeU a ^ 2 *
        jacobianDet ![
          quadraticGaugePi a, quadraticGaugeBaseB a c, quadraticGaugeU a] := by
  classical
  have hunit :
      (MvPolynomial.C a : GaugePolynomial K) * MvPolynomial.C a⁻¹ = 1 := by
    rw [← MvPolynomial.C_mul]
    simp [ha]
  simp only [jacobianDet, jacobianMatrix, det_fin_three, of_apply,
    quadraticGaugeBaseB, quadraticGaugeBaseC,
    quadraticGaugePi, quadraticGaugeU, quadraticGaugeQ, generalGaugeT,
    cons_val_zero, cons_val_one, cons_val_two, head_cons, tail_cons,
    map_add, map_sub, Derivation.map_one_eq_zero,
    pderiv_mul, pderiv_pow, pderiv_C, pderiv_X_self, pderiv_X_of_ne,
    ne_eq, Fin.reduceEq, not_false_eq_true]
  simp only [map_neg, map_ofNat, MvPolynomial.C_mul,
    MvPolynomial.C_eq_coe_nat]
  grobner

private theorem pderiv_quadraticGaugeTailB
    (a : K) (R Rp Rpp : GaugePolynomial K)
    (hR : ∀ i, pderiv i R = Rp * pderiv i (quadraticGaugeU a))
    (hRp : ∀ i, pderiv i Rp = Rpp * pderiv i (quadraticGaugeU a))
    (i : Fin 3) :
    pderiv i (quadraticGaugeTailB a R Rp) =
      (MvPolynomial.C 2 * quadraticGaugePi a *
          (MvPolynomial.C 2 * R + quadraticGaugeU a * Rp)) *
          pderiv i (quadraticGaugePi a) +
        quadraticGaugePi a ^ 2 *
          (MvPolynomial.C 3 * Rp + quadraticGaugeU a * Rpp) *
            pderiv i (quadraticGaugeU a) := by
  simp only [quadraticGaugeTailB, map_add, pderiv_mul, pderiv_pow,
    pderiv_C, zero_mul, zero_add, hR i, hRp i]
  simp only [map_ofNat]
  ring

private theorem pderiv_quadraticGaugeTailC
    (a : K) (Rp Rpp : GaugePolynomial K)
    (hRp : ∀ i, pderiv i Rp = Rpp * pderiv i (quadraticGaugeU a))
    (i : Fin 3) :
    pderiv i (quadraticGaugeTailC a Rp) =
      -quadraticGaugeU a ^ 2 *
        (MvPolynomial.C 3 * Rp + quadraticGaugeU a * Rpp) *
          pderiv i (quadraticGaugeU a) := by
  simp only [quadraticGaugeTailC, map_neg, pderiv_mul, pderiv_pow,
    hRp i, map_ofNat]
  ring

set_option maxHeartbeats 0 in
private theorem jacobianDet_quadraticGaugeWithTail
    (a c : K) (ha : a ≠ 0)
    (R Rp Rpp : GaugePolynomial K)
    (hR : ∀ i, pderiv i R = Rp * pderiv i (quadraticGaugeU a))
    (hRp : ∀ i, pderiv i Rp = Rpp * pderiv i (quadraticGaugeU a)) :
    jacobianDet (quadraticGaugeWithTail a c R Rp) = MvPolynomial.C (-2) := by
  classical
  have hB0 := pderiv_quadraticGaugeTailB a R Rp Rpp hR hRp (0 : Fin 3)
  have hB1 := pderiv_quadraticGaugeTailB a R Rp Rpp hR hRp (1 : Fin 3)
  have hB2 := pderiv_quadraticGaugeTailB a R Rp Rpp hR hRp (2 : Fin 3)
  have hC0 := pderiv_quadraticGaugeTailC a Rp Rpp hRp (0 : Fin 3)
  have hC1 := pderiv_quadraticGaugeTailC a Rp Rpp hRp (1 : Fin 3)
  have hC2 := pderiv_quadraticGaugeTailC a Rp Rpp hRp (2 : Fin 3)
  have hbase := jacobianDet_quadraticGaugeBase a c ha
  have hcross := quadraticGauge_crossDet a c ha
  simp only [jacobianDet, jacobianMatrix, det_fin_three, of_apply,
    quadraticGaugeWithTail, quadraticGaugeBaseMap,
    cons_val_zero, cons_val_one, cons_val_two, head_cons, tail_cons,
    map_add] at hbase hcross ⊢
  rw [hB0, hB1, hB2, hC0, hC1, hC2]
  linear_combination hbase +
    (MvPolynomial.C 3 * Rp + quadraticGaugeU a * Rpp) * hcross

set_option maxHeartbeats 0 in
/-- The actual all-degree quadratic-gauge map has constant Jacobian `-2` under
the two nonvanishing hypotheses used by the construction. -/
theorem jacobianDet_generalGaugeMap
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    jacobianDet (generalGaugeMap G) = MvPolynomial.C (-2) := by
  let a : K := G.coeff 1 / G.coeff 3
  let c : K := G.coeff 2 / G.coeff 1
  let d : ℕ → K := fun k => G.coeff k / G.coeff 1
  let u : GaugePolynomial K := quadraticGaugeU a
  let R : GaugePolynomial K := quadraticGaugeTail d G.natDegree u
  let Rp : GaugePolynomial K := quadraticGaugeTailDeriv d G.natDegree u
  let Rpp : GaugePolynomial K := quadraticGaugeTailSecond d G.natDegree u
  have ha : a ≠ 0 := by
    exact div_ne_zero h₁ h₃
  have hBtail :
      quadraticGaugePi a ^ 2 *
          (MvPolynomial.C 2 * R + quadraticGaugeU a * Rp) =
        ∑ k ∈ Finset.Icc 4 G.natDegree,
          MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
            generalGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
              generalGaugeQ G ^ k := by
    simpa [a, d, u, R, Rp, quadraticGaugePi, quadraticGaugeU,
      quadraticGaugeQ, generalGaugeQ] using
      quadraticGaugeTail_B d G.natDegree generalGaugeT
        (MvPolynomial.X 0) (quadraticGaugeQ a)
  have hCtail :
      quadraticGaugeU a ^ 3 * Rp =
        ∑ k ∈ Finset.Icc 4 G.natDegree,
          MvPolynomial.C
              (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
            (MvPolynomial.X 0 * generalGaugeQ G) ^ k := by
    simpa [a, d, u, Rp, quadraticGaugeU, quadraticGaugeQ, generalGaugeQ] using
      quadraticGaugeTail_C d G.natDegree (quadraticGaugeU a)
  have hmap : generalGaugeMap G = quadraticGaugeWithTail a c R Rp := by
    funext i
    fin_cases i
    · simp [generalGaugeMap, quadraticGaugeWithTail, generalGaugePi,
        quadraticGaugePi, generalGaugeQ, quadraticGaugeQ, a]
    · change generalGaugeB G =
        quadraticGaugeBaseB a c + quadraticGaugeTailB a R Rp
      rw [quadraticGaugeBaseB, quadraticGaugeTailB, hBtail]
      simp only [generalGaugeB, quadraticGaugeU, quadraticGaugePi,
        quadraticGaugeQ, generalGaugeQ, generalGaugePi, a, c]
      rw [inv_div]
      ring
    · change generalGaugeC G =
        quadraticGaugeBaseC a + quadraticGaugeTailC a Rp
      rw [quadraticGaugeBaseC, quadraticGaugeTailC, hCtail]
      simp only [generalGaugeC, quadraticGaugeU, quadraticGaugeQ,
        generalGaugeQ, a]
      rw [inv_div]
      ring
  rw [hmap]
  apply jacobianDet_quadraticGaugeWithTail a c ha R Rp Rpp
  · intro i
    simpa [R, Rp, u] using pderiv_quadraticGaugeTail d G.natDegree u i
  · intro i
    simpa [Rp, Rpp, u] using pderiv_quadraticGaugeTailDeriv d G.natDegree u i

/-- The fixed target-preserving output normalization has Jacobian one. -/
theorem jacobianDet_generalGaugeJacobianOneMap [CharZero K]
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    jacobianDet (generalGaugeJacobianOneMap G) = 1 := by
  rw [generalGaugeJacobianOneMap, jacobianDet_scaleOutput,
    jacobianDet_generalGaugeMap G h₁ h₃]
  rw [← MvPolynomial.C_mul]
  norm_num

#print axioms jacobianDet_generalGaugeMap
#print axioms jacobianDet_generalGaugeJacobianOneMap

end FiniteEtaleKeller
