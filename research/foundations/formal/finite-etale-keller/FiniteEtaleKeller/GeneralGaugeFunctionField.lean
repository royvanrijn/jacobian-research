/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeJacobian
import FiniteEtaleKeller.GenericInverseIrreducibility

/-!
# The function-field map of the all-degree quadratic gauge

This module begins the bridge from the generic inverse polynomial to the
function-field definition of geometric degree.  It packages substitution by
the displayed coordinates as an algebra homomorphism and proves that this
homomorphism is injective.  Thus the three target coordinates are
algebraically independent.
-/

noncomputable section

open Matrix Function
open MvPolynomial
open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Substitution by a polynomial self-map. -/
def coordinateSubstitution
    {σ : Type*} (F : σ → MvPolynomial σ K) :
    MvPolynomial σ K →ₐ[K] MvPolynomial σ K :=
  MvPolynomial.aeval F

@[simp]
theorem coordinateSubstitution_X
    {σ : Type*} (F : σ → MvPolynomial σ K) (i : σ) :
    coordinateSubstitution F (MvPolynomial.X i) = F i := by
  simp [coordinateSubstitution]

/-- The coordinate-ring homomorphism induced by the all-degree gauge map. -/
def generalGaugeCoordinateHom (G : K[X]) :
    GaugePolynomial K →ₐ[K] GaugePolynomial K :=
  coordinateSubstitution (generalGaugeMap G)

@[simp]
theorem generalGaugeCoordinateHom_X (G : K[X]) (i : Fin 3) :
    generalGaugeCoordinateHom G (MvPolynomial.X i) =
      generalGaugeMap G i := by
  simp [generalGaugeCoordinateHom]

private theorem pderiv_aeval_chain_rule
    {σ : Type*} [Fintype σ]
    (F : σ → MvPolynomial σ K) (P : MvPolynomial σ K) (i : σ) :
    pderiv i (MvPolynomial.aeval F P) =
      ∑ j, MvPolynomial.aeval F (pderiv j P) * pderiv i (F j) := by
  classical
  induction P using MvPolynomial.induction_on with
  | C a =>
      simp
  | add P Q hP hQ =>
      simp only [map_add, hP, hQ]
      simp_rw [add_mul]
      rw [Finset.sum_add_distrib]
  | mul_X P j hP =>
      classical
      simp only [map_mul, MvPolynomial.aeval_X, pderiv_mul, hP, map_add]
      simp only [pderiv_X, Pi.single_apply]
      simp_rw [add_mul]
      rw [Finset.sum_add_distrib]
      congr 1
      · rw [Finset.sum_mul]
        apply Finset.sum_congr rfl
        intro x hx
        ring
      · simp

private theorem eq_C_of_forall_pderiv_eq_zero
    {σ : Type*} [CharZero K] (P : MvPolynomial σ K)
    (hP : ∀ i, pderiv i P = 0) :
    P = MvPolynomial.C (MvPolynomial.coeff 0 P) := by
  classical
  ext m
  by_cases hm : m = 0
  · subst m
    simp
  · obtain ⟨i, hi⟩ : ∃ i, m i ≠ 0 := by
      simpa [Finsupp.ext_iff] using hm
    let n := m - Finsupp.single i 1
    have hmi : 1 ≤ m i := Nat.one_le_iff_ne_zero.mpr hi
    have hnm : n + Finsupp.single i 1 = m := by
      exact Finsupp.sub_add_single_one_cancel hi
    have hcoeff :
        MvPolynomial.coeff m P * (n i + 1 : ℕ) = 0 := by
      have := congrArg (MvPolynomial.coeff n) (hP i)
      simpa [MvPolynomial.coeff_pderiv, hnm] using this
    have hcast : (n i + 1 : K) ≠ 0 := by
      exact_mod_cast Nat.succ_ne_zero (n i)
    have hcoeff0 : MvPolynomial.coeff m P = 0 := by
      exact (mul_eq_zero.mp hcoeff).resolve_right (by exact_mod_cast hcast)
    simp [Ne.symm hm, hcoeff0]

private theorem totalDegree_pderiv_lt
    {σ : Type*}
    (P : MvPolynomial σ K) (i : σ)
    (hPi : pderiv i P ≠ 0) :
    (pderiv i P).totalDegree < P.totalDegree := by
  classical
  rw [MvPolynomial.totalDegree, Finset.sup_lt_iff]
  · intro m hm
    have hcoeff :
        MvPolynomial.coeff (m + Finsupp.single i 1) P *
            ((m i : K) + 1) ≠ 0 := by
      simpa only [MvPolynomial.coeff_pderiv] using
        (MvPolynomial.mem_support_iff.mp hm)
    have hmem :
        m + Finsupp.single i 1 ∈ P.support := by
      rw [MvPolynomial.mem_support_iff]
      exact fun hz => hcoeff (by rw [hz, zero_mul])
    have hsum :
        (m + Finsupp.single i 1).sum (fun _ e => e) =
          m.sum (fun _ e => e) + 1 := by
      simp [Finsupp.sum_add_index']
    rw [← Nat.add_one_le_iff, ← hsum]
    exact MvPolynomial.le_totalDegree hmem
  · have hpositive : 0 < P.totalDegree := by
      obtain ⟨m, hm⟩ := MvPolynomial.support_nonempty.mpr hPi
      have hcoeff :
          MvPolynomial.coeff (m + Finsupp.single i 1) P *
              ((m i : K) + 1) ≠ 0 := by
        simpa only [MvPolynomial.coeff_pderiv] using
          (MvPolynomial.mem_support_iff.mp hm)
      have hmem :
          m + Finsupp.single i 1 ∈ P.support := by
        rw [MvPolynomial.mem_support_iff]
        exact fun hz => hcoeff (by rw [hz, zero_mul])
      have hle := MvPolynomial.le_totalDegree hmem
      have hsum :
          (m + Finsupp.single i 1).sum (fun _ e => e) =
            m.sum (fun _ e => e) + 1 := by
        simp [Finsupp.sum_add_index']
      rw [hsum] at hle
      omega
    exact hpositive

/-- A polynomial self-map with nonzero Jacobian determinant induces an
injective substitution homomorphism in characteristic zero. -/
theorem coordinateSubstitution_injective_of_jacobianDet_ne_zero
    {σ : Type*} [Fintype σ] [DecidableEq σ] [CharZero K]
    (F : σ → MvPolynomial σ K)
    (hdet : jacobianDet F ≠ 0) :
    Function.Injective (coordinateSubstitution F) := by
  refine (injective_iff_map_eq_zero (coordinateSubstitution F)).mpr ?_
  intro P
  induction hdegree : P.totalDegree using Nat.strong_induction_on generalizing P with
  | h n ih =>
    intro hP
    have hchain (i : σ) :
        ∑ j, MvPolynomial.aeval F (pderiv j P) * pderiv i (F j) = 0 := by
      rw [← pderiv_aeval_chain_rule F P i, show MvPolynomial.aeval F P = 0 from hP]
      simp
    let v : σ → MvPolynomial σ K :=
      fun j => MvPolynomial.aeval F (pderiv j P)
    have hvec : v ᵥ* jacobianMatrix F = 0 := by
      funext i
      simpa [v, Matrix.vecMul, dotProduct, jacobianMatrix] using hchain i
    have hv : v = 0 :=
      Matrix.eq_zero_of_vecMul_eq_zero hdet hvec
    have hpderiv (j : σ) :
        coordinateSubstitution F (pderiv j P) = 0 := by
      exact congrFun hv j
    have hpzero (j : σ) : pderiv j P = 0 := by
      by_contra hp
      exact hp <| ih (pderiv j P).totalDegree
        (by
          rw [← hdegree]
          exact totalDegree_pderiv_lt P j hp)
        (pderiv j P) rfl (hpderiv j)
    have hconst := eq_C_of_forall_pderiv_eq_zero P hpzero
    rw [hconst] at hP ⊢
    simpa using hP

/-- The three coordinates of the all-degree gauge map are algebraically
independent. -/
theorem generalGaugeMap_algebraicIndependent
    [CharZero K] (G : K[X])
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    AlgebraicIndependent K (generalGaugeMap G) := by
  rw [algebraicIndependent_iff_injective_aeval]
  apply coordinateSubstitution_injective_of_jacobianDet_ne_zero
  rw [jacobianDet_generalGaugeMap G h₁ h₃]
  exact MvPolynomial.C_ne_zero.mpr (by norm_num)

/-- The rational function field in the three source (or target) variables. -/
abbrev GaugeFunctionField (K : Type*) [Field K] :=
  FractionRing (GaugePolynomial K)

/-- Coordinate substitution followed by the canonical inclusion into the
source rational function field. -/
def generalGaugeCoordinateHomToFunctionField (G : K[X]) :
    GaugePolynomial K →ₐ[K] GaugeFunctionField K :=
  (IsScalarTower.toAlgHom K (GaugePolynomial K) (GaugeFunctionField K)).comp
    (generalGaugeCoordinateHom G)

theorem generalGaugeCoordinateHomToFunctionField_injective
    [CharZero K] (G : K[X])
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    Function.Injective (generalGaugeCoordinateHomToFunctionField G) :=
  (IsFractionRing.injective (GaugePolynomial K) (GaugeFunctionField K)).comp
    (algebraicIndependent_iff_injective_aeval.mp
      (generalGaugeMap_algebraicIndependent G h₁ h₃))

/-- The injective coordinate-ring map extends canonically to the rational
function fields.  This is the pullback on function fields induced by the
all-degree gauge map. -/
def generalGaugeFunctionFieldHom
    [CharZero K] (G : K[X])
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    GaugeFunctionField K →ₐ[K] GaugeFunctionField K := by
  exact IsFractionRing.liftAlgHom
    (generalGaugeCoordinateHomToFunctionField_injective G h₁ h₃)

@[simp]
theorem generalGaugeFunctionFieldHom_algebraMap
    [CharZero K] (G : K[X])
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (P : GaugePolynomial K) :
    generalGaugeFunctionFieldHom G h₁ h₃
        (algebraMap (GaugePolynomial K) (GaugeFunctionField K) P) =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K)
        (generalGaugeCoordinateHom G P) := by
  rw [generalGaugeFunctionFieldHom, IsFractionRing.liftAlgHom_apply,
    IsFractionRing.lift_algebraMap]
  rfl

/-- The induced function-field pullback is injective. -/
theorem generalGaugeFunctionFieldHom_injective
    [CharZero K] (G : K[X])
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    Function.Injective (generalGaugeFunctionFieldHom G h₁ h₃) :=
  (generalGaugeFunctionFieldHom G h₁ h₃).injective

#print axioms coordinateSubstitution_injective_of_jacobianDet_ne_zero
#print axioms generalGaugeMap_algebraicIndependent
#print axioms generalGaugeFunctionFieldHom_injective

end FiniteEtaleKeller
