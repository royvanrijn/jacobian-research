/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeFullGenericDegree
import FiniteEtaleKeller.GeneralGaugeRawFiber
import FiniteEtaleKeller.GeneralGaugeBaseChange

/-!
# Comparison of the inverse-root and pullback function-field extensions

This module identifies the iterated target presentation `K(Π,B)(C)` with the
three-variable target function field and then identifies the fully generic
inverse-root quotient with the actual source function field.
-/

noncomputable section

open Polynomial MvPolynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K] [CharZero K]

/-- The three-cycle which orders the target variables as `(C, Π, B)` for
`MvPolynomial.finSuccEquiv`. -/
def gaugeTargetCycle : Fin 3 ≃ Fin 3 :=
  (Equiv.swap (1 : Fin 3) 2).trans (Equiv.swap 0 1)

@[simp] theorem gaugeTargetCycle_zero : gaugeTargetCycle 0 = 1 := by
  change (Equiv.swap 0 1) ((Equiv.swap 1 2) 0) = 1
  decide
@[simp] theorem gaugeTargetCycle_one : gaugeTargetCycle 1 = 2 := by
  change (Equiv.swap 0 1) ((Equiv.swap 1 2) 1) = 2
  decide
@[simp] theorem gaugeTargetCycle_two : gaugeTargetCycle 2 = 0 := by
  change (Equiv.swap 0 1) ((Equiv.swap 1 2) 2) = 0
  decide

/-- Polynomial coordinates `K[Π,B,C]` written as `K[Π,B][C]`. -/
def gaugeTargetPolynomialEquiv :
    GaugePolynomial K ≃ₐ[K]
      Polynomial (GaugeTargetParameterPolynomial K) :=
  (MvPolynomial.renameEquiv K gaugeTargetCycle).trans
    (MvPolynomial.finSuccEquiv K 2)

@[simp]
theorem gaugeTargetPolynomialEquiv_X_zero :
    gaugeTargetPolynomialEquiv
        (MvPolynomial.X 0 : GaugePolynomial K) =
      Polynomial.C (MvPolynomial.X 0) := by
  calc
    gaugeTargetPolynomialEquiv
        (MvPolynomial.X 0 : GaugePolynomial K) =
      (MvPolynomial.finSuccEquiv K 2)
        (MvPolynomial.X (gaugeTargetCycle 0)) := by
          simp only [gaugeTargetPolynomialEquiv, AlgEquiv.trans_apply,
            MvPolynomial.renameEquiv_apply, MvPolynomial.rename_X]
    _ = Polynomial.C (MvPolynomial.X 0) := by
      rw [gaugeTargetCycle_zero]
      have h : (1 : Fin 3) = (0 : Fin 2).succ := by decide
      rw [h]
      exact MvPolynomial.finSuccEquiv_X_succ

@[simp]
theorem gaugeTargetPolynomialEquiv_X_one :
    gaugeTargetPolynomialEquiv
        (MvPolynomial.X 1 : GaugePolynomial K) =
      Polynomial.C (MvPolynomial.X 1) := by
  calc
    gaugeTargetPolynomialEquiv
        (MvPolynomial.X 1 : GaugePolynomial K) =
      (MvPolynomial.finSuccEquiv K 2)
        (MvPolynomial.X (gaugeTargetCycle 1)) := by
          simp only [gaugeTargetPolynomialEquiv, AlgEquiv.trans_apply,
            MvPolynomial.renameEquiv_apply, MvPolynomial.rename_X]
    _ = Polynomial.C (MvPolynomial.X 1) := by
      rw [gaugeTargetCycle_one]
      have h : (2 : Fin 3) = (1 : Fin 2).succ := by decide
      rw [h]
      exact MvPolynomial.finSuccEquiv_X_succ

@[simp]
theorem gaugeTargetPolynomialEquiv_X_two :
    gaugeTargetPolynomialEquiv
        (MvPolynomial.X 2 : GaugePolynomial K) =
      Polynomial.X := by
  calc
    gaugeTargetPolynomialEquiv
        (MvPolynomial.X 2 : GaugePolynomial K) =
      (MvPolynomial.finSuccEquiv K 2)
        (MvPolynomial.X (gaugeTargetCycle 2)) := by
          simp only [gaugeTargetPolynomialEquiv, AlgEquiv.trans_apply,
            MvPolynomial.renameEquiv_apply, MvPolynomial.rename_X]
    _ = Polynomial.X := by
      rw [gaugeTargetCycle_two]
      exact MvPolynomial.finSuccEquiv_X_zero

/-- The one-step fraction field of `K[Π,B][C]`. -/
abbrev GaugeTargetOneStepField (K : Type*) [Field K] :=
  FractionRing (Polynomial (GaugeTargetParameterPolynomial K))

/-- The three-variable target function field in the one-step
`Frac(K[Π,B][C])` presentation. -/
def gaugeFunctionFieldOneStepEquiv :
    GaugeFunctionField K ≃ₐ[K] GaugeTargetOneStepField K :=
  IsFractionRing.algEquivOfAlgEquiv gaugeTargetPolynomialEquiv

/-- Embed `K[Π,B][C]` into the iterated fraction field `K(Π,B)(C)`. -/
def gaugeOneStepPolynomialToIterated :
    Polynomial (GaugeTargetParameterPolynomial K) →ₐ[K]
      RatFunc (GaugeTargetParameterField K) :=
  (IsScalarTower.toAlgHom K
      (Polynomial (GaugeTargetParameterField K))
      (RatFunc (GaugeTargetParameterField K))).comp
    (Polynomial.mapAlgHom
      (IsScalarTower.toAlgHom K
        (GaugeTargetParameterPolynomial K)
        (GaugeTargetParameterField K)))

@[simp]
theorem gaugeOneStepPolynomialToIterated_apply
    (p : Polynomial (GaugeTargetParameterPolynomial K)) :
    gaugeOneStepPolynomialToIterated p =
      algebraMap (Polynomial (GaugeTargetParameterField K))
        (RatFunc (GaugeTargetParameterField K))
        (p.map (algebraMap
          (GaugeTargetParameterPolynomial K)
          (GaugeTargetParameterField K))) := rfl

theorem gaugeOneStepPolynomialToIterated_injective :
    Function.Injective (gaugeOneStepPolynomialToIterated (K := K)) :=
  (IsFractionRing.injective
      (Polynomial (GaugeTargetParameterField K))
      (RatFunc (GaugeTargetParameterField K))).comp
    (Polynomial.map_injective _
      (IsFractionRing.injective
        (GaugeTargetParameterPolynomial K)
        (GaugeTargetParameterField K)))

/-- The fraction-field lift from the one-step to the iterated target
presentation. -/
def gaugeOneStepToIterated :
    GaugeTargetOneStepField K →ₐ[K]
      RatFunc (GaugeTargetParameterField K) :=
  IsFractionRing.liftAlgHom gaugeOneStepPolynomialToIterated_injective

/-- Embed `K[Π,B]` into the one-step target fraction field. -/
def gaugeParameterPolynomialToOneStep :
    GaugeTargetParameterPolynomial K →ₐ[K]
      GaugeTargetOneStepField K :=
  (IsScalarTower.toAlgHom K
      (Polynomial (GaugeTargetParameterPolynomial K))
      (GaugeTargetOneStepField K)).comp
    Polynomial.CAlgHom

theorem gaugeParameterPolynomialToOneStep_injective :
    Function.Injective (gaugeParameterPolynomialToOneStep (K := K)) :=
  (IsFractionRing.injective
      (Polynomial (GaugeTargetParameterPolynomial K))
      (GaugeTargetOneStepField K)).comp
    (fun _ _ h => Polynomial.C_injective h)

/-- Embed the coefficient field `K(Π,B)` into the one-step target field. -/
def gaugeParameterFieldToOneStep :
    GaugeTargetParameterField K →ₐ[K] GaugeTargetOneStepField K :=
  IsFractionRing.liftAlgHom gaugeParameterPolynomialToOneStep_injective

/-- Evaluate a polynomial over `K(Π,B)` at the one-step class of `C`. -/
def gaugeIteratedPolynomialToOneStep :
    Polynomial (GaugeTargetParameterField K) →ₐ[K]
      GaugeTargetOneStepField K :=
  Polynomial.eval₂AlgHom gaugeParameterFieldToOneStep
    (algebraMap
      (Polynomial (GaugeTargetParameterPolynomial K))
      (GaugeTargetOneStepField K) Polynomial.X)
    (fun a => (commute_iff_eq _ _).2 (mul_comm _ _))

theorem gaugeIteratedPolynomialToOneStep_leftInverse :
    (gaugeOneStepToIterated (K := K)).comp
        (gaugeIteratedPolynomialToOneStep (K := K)) =
      IsScalarTower.toAlgHom K
        (Polynomial (GaugeTargetParameterField K))
        (RatFunc (GaugeTargetParameterField K)) := by
  apply Polynomial.algHom_ext'
  · apply AlgHom.coe_ringHom_injective
    apply IsFractionRing.ringHom_ext
      (A := GaugeTargetParameterPolynomial K)
    intro p
    simp [gaugeIteratedPolynomialToOneStep, gaugeOneStepToIterated,
      gaugeParameterFieldToOneStep, gaugeParameterPolynomialToOneStep,
      gaugeOneStepPolynomialToIterated]
  · simp [gaugeIteratedPolynomialToOneStep, gaugeOneStepToIterated,
      gaugeParameterFieldToOneStep, gaugeParameterPolynomialToOneStep,
      gaugeOneStepPolynomialToIterated]

theorem gaugeIteratedPolynomialToOneStep_injective :
    Function.Injective (gaugeIteratedPolynomialToOneStep (K := K)) := by
  intro p q h
  apply (IsFractionRing.injective
    (Polynomial (GaugeTargetParameterField K))
    (RatFunc (GaugeTargetParameterField K)))
  have hh := congrArg (gaugeOneStepToIterated (K := K)) h
  have hp := congrArg (fun f => f p)
    (gaugeIteratedPolynomialToOneStep_leftInverse (K := K))
  have hq := congrArg (fun f => f q)
    (gaugeIteratedPolynomialToOneStep_leftInverse (K := K))
  exact hp.symm.trans (hh.trans hq)

/-- The fraction-field lift from the iterated to the one-step target
presentation. -/
def gaugeIteratedToOneStep :
    RatFunc (GaugeTargetParameterField K) →ₐ[K]
      GaugeTargetOneStepField K :=
  RatFunc.liftAlgHom gaugeIteratedPolynomialToOneStep
    (nonZeroDivisors_le_comap_nonZeroDivisors_of_injective _
      gaugeIteratedPolynomialToOneStep_injective)

@[simp]
theorem gaugeIteratedToOneStep_algebraMap
    (p : Polynomial (GaugeTargetParameterField K)) :
    gaugeIteratedToOneStep
        (algebraMap (Polynomial (GaugeTargetParameterField K))
          (RatFunc (GaugeTargetParameterField K)) p) =
      gaugeIteratedPolynomialToOneStep p := by
  exact RatFunc.liftRingHom_algebraMap _ _ p

@[simp]
theorem gaugeOneStepToIterated_algebraMap
    (p : Polynomial (GaugeTargetParameterPolynomial K)) :
    gaugeOneStepToIterated
        (algebraMap (Polynomial (GaugeTargetParameterPolynomial K))
          (GaugeTargetOneStepField K) p) =
      gaugeOneStepPolynomialToIterated p :=
  IsFractionRing.lift_algebraMap
    gaugeOneStepPolynomialToIterated_injective p

theorem gaugeOneStepIterated_left :
    (gaugeOneStepToIterated (K := K)).comp
        (gaugeIteratedToOneStep (K := K)) =
      AlgHom.id K (RatFunc (GaugeTargetParameterField K)) := by
  apply AlgHom.coe_ringHom_injective
  apply IsFractionRing.ringHom_ext
    (A := Polynomial (GaugeTargetParameterField K))
  intro p
  simp only [RingHom.comp_apply, AlgHom.toRingHom_eq_coe,
    AlgHom.coe_toRingHom, AlgHom.comp_apply, AlgHom.id_apply,
    gaugeIteratedToOneStep_algebraMap]
  exact congrArg (fun f => f p)
      (gaugeIteratedPolynomialToOneStep_leftInverse (K := K))

theorem gaugeOneStepIterated_right :
    (gaugeIteratedToOneStep (K := K)).comp
        (gaugeOneStepToIterated (K := K)) =
      AlgHom.id K (GaugeTargetOneStepField K) := by
  apply AlgHom.coe_ringHom_injective
  apply IsFractionRing.ringHom_ext
    (A := Polynomial (GaugeTargetParameterPolynomial K))
  intro p
  apply (gaugeOneStepToIterated (K := K)).injective
  have h := congrArg (fun f => f
      (Polynomial.map
        (algebraMap (GaugeTargetParameterPolynomial K)
          (GaugeTargetParameterField K)) p))
    (gaugeIteratedPolynomialToOneStep_leftInverse (K := K))
  simp only [RingHom.comp_apply, AlgHom.toRingHom_eq_coe,
    AlgHom.coe_toRingHom, AlgHom.comp_apply, AlgHom.id_apply,
    gaugeOneStepToIterated_algebraMap, gaugeIteratedToOneStep_algebraMap]
  simpa [AlgHom.comp_apply, gaugeOneStepPolynomialToIterated_apply] using h

/-- Canonical equivalence between `Frac(K[Π,B][C])` and `K(Π,B)(C)`. -/
def gaugeOneStepIteratedEquiv :
    GaugeTargetOneStepField K ≃ₐ[K]
      RatFunc (GaugeTargetParameterField K) :=
  AlgEquiv.ofAlgHom gaugeOneStepToIterated gaugeIteratedToOneStep
    gaugeOneStepIterated_left gaugeOneStepIterated_right

/-- Canonical target presentation
`K(Π,B,C) ≃ K(Π,B)(C)`. -/
def gaugeTargetPresentationEquiv :
    GaugeFunctionField K ≃ₐ[K]
      RatFunc (GaugeTargetParameterField K) :=
  gaugeFunctionFieldOneStepEquiv.trans gaugeOneStepIteratedEquiv

@[simp]
theorem gaugeTargetPresentationEquiv_X_zero :
    gaugeTargetPresentationEquiv
        (algebraMap (GaugePolynomial K) (GaugeFunctionField K)
          (MvPolynomial.X 0)) =
      algebraMap (GaugeTargetParameterField K)
        (RatFunc (GaugeTargetParameterField K)) (gaugeGenericPi K) := by
  simp [gaugeTargetPresentationEquiv, gaugeFunctionFieldOneStepEquiv,
    gaugeOneStepIteratedEquiv, gaugeOneStepToIterated,
    gaugeOneStepPolynomialToIterated, gaugeGenericPi]

@[simp]
theorem gaugeTargetPresentationEquiv_X_one :
    gaugeTargetPresentationEquiv
        (algebraMap (GaugePolynomial K) (GaugeFunctionField K)
          (MvPolynomial.X 1)) =
      algebraMap (GaugeTargetParameterField K)
        (RatFunc (GaugeTargetParameterField K)) (gaugeGenericB K) := by
  simp [gaugeTargetPresentationEquiv, gaugeFunctionFieldOneStepEquiv,
    gaugeOneStepIteratedEquiv, gaugeOneStepToIterated,
    gaugeOneStepPolynomialToIterated, gaugeGenericB]

@[simp]
theorem gaugeTargetPresentationEquiv_X_two :
    gaugeTargetPresentationEquiv
        (algebraMap (GaugePolynomial K) (GaugeFunctionField K)
          (MvPolynomial.X 2)) =
      (RatFunc.X : RatFunc (GaugeTargetParameterField K)) := by
  simp only [gaugeTargetPresentationEquiv, gaugeFunctionFieldOneStepEquiv,
    AlgEquiv.trans_apply,
    IsFractionRing.algEquivOfAlgEquiv_algebraMap,
    gaugeTargetPolynomialEquiv_X_two,
    gaugeOneStepIteratedEquiv, AlgEquiv.ofAlgHom_apply,
    gaugeOneStepToIterated_algebraMap,
    gaugeOneStepPolynomialToIterated_apply,
    Polynomial.map_X]
  exact RatFunc.algebraMap_X

/-- The fully generic target field acts on the actual source function field
through the displayed map. -/
def generalGaugeFullyGenericTargetHom
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    RatFunc (GaugeTargetParameterField K) →ₐ[K] GaugeFunctionField K :=
  (generalGaugeFunctionFieldHom G h₁ h₃).comp
    gaugeTargetPresentationEquiv.symm.toAlgHom

theorem generalGaugeFullyGenericTargetHom_injective
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    Function.Injective
      (generalGaugeFullyGenericTargetHom G h₁ h₃) :=
  (generalGaugeFunctionFieldHom_injective G h₁ h₃).comp
    gaugeTargetPresentationEquiv.symm.injective

@[simp]
theorem generalGaugeFullyGenericTargetHom_pi
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    generalGaugeFullyGenericTargetHom G h₁ h₃
        (algebraMap (GaugeTargetParameterField K)
          (RatFunc (GaugeTargetParameterField K)) (gaugeGenericPi K)) =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K)
        (generalGaugePi G) := by
  change (generalGaugeFunctionFieldHom G h₁ h₃)
      (gaugeTargetPresentationEquiv.symm
        (algebraMap (GaugeTargetParameterField K)
          (RatFunc (GaugeTargetParameterField K)) (gaugeGenericPi K))) = _
  rw [← gaugeTargetPresentationEquiv_X_zero (K := K),
    gaugeTargetPresentationEquiv.symm_apply_apply,
    generalGaugeFunctionFieldHom_algebraMap]
  simp [generalGaugeCoordinateHom_X, generalGaugeMap]

@[simp]
theorem generalGaugeFullyGenericTargetHom_b
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    generalGaugeFullyGenericTargetHom G h₁ h₃
        (algebraMap (GaugeTargetParameterField K)
          (RatFunc (GaugeTargetParameterField K)) (gaugeGenericB K)) =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K)
        (generalGaugeB G) := by
  change (generalGaugeFunctionFieldHom G h₁ h₃)
      (gaugeTargetPresentationEquiv.symm
        (algebraMap (GaugeTargetParameterField K)
          (RatFunc (GaugeTargetParameterField K)) (gaugeGenericB K))) = _
  rw [← gaugeTargetPresentationEquiv_X_one (K := K),
    gaugeTargetPresentationEquiv.symm_apply_apply,
    generalGaugeFunctionFieldHom_algebraMap]
  simp [generalGaugeCoordinateHom_X, generalGaugeMap]

@[simp]
theorem generalGaugeFullyGenericTargetHom_c
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    generalGaugeFullyGenericTargetHom G h₁ h₃
        (RatFunc.X : RatFunc (GaugeTargetParameterField K)) =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K)
        (generalGaugeC G) := by
  change (generalGaugeFunctionFieldHom G h₁ h₃)
      (gaugeTargetPresentationEquiv.symm
        (RatFunc.X : RatFunc (GaugeTargetParameterField K))) = _
  rw [← gaugeTargetPresentationEquiv_X_two (K := K),
    gaugeTargetPresentationEquiv.symm_apply_apply,
    generalGaugeFunctionFieldHom_algebraMap]
  simp [generalGaugeCoordinateHom_X, generalGaugeMap]

/-- The seed over the full target field, obtained by the two successive
coefficient extensions `K → K(Π,B) → K(Π,B)(C)`. -/
def generalGaugeFullyGenericSeedOverTarget (G : K[X]) :
    (RatFunc (GaugeTargetParameterField K))[X] :=
  (generalGaugeGenericSeed G).map
    (algebraMap (GaugeTargetParameterField K)
      (RatFunc (GaugeTargetParameterField K)))

/-- The fully generic first target coordinate. -/
def gaugeFullyGenericPi : RatFunc (GaugeTargetParameterField K) :=
  algebraMap (GaugeTargetParameterField K)
    (RatFunc (GaugeTargetParameterField K)) (gaugeGenericPi K)

/-- The fully generic second target coordinate. -/
def gaugeFullyGenericB : RatFunc (GaugeTargetParameterField K) :=
  algebraMap (GaugeTargetParameterField K)
    (RatFunc (GaugeTargetParameterField K)) (gaugeGenericB K)

theorem gaugeFullyGenericPi_ne_zero :
    gaugeFullyGenericPi (K := K) ≠ 0 :=
  by
    simpa [gaugeFullyGenericPi] using
      (algebraMap (GaugeTargetParameterField K)
        (RatFunc (GaugeTargetParameterField K))).injective.ne
          (gaugeGenericPi_ne_zero K)

/-- The first generic target coordinate carried as a unit. -/
def gaugeFullyGenericPiUnit : (RatFunc (GaugeTargetParameterField K))ˣ :=
  Units.mk0 gaugeFullyGenericPi gaugeFullyGenericPi_ne_zero

@[simp]
theorem gaugeFullyGenericPiUnit_val :
    (gaugeFullyGenericPiUnit (K := K) :
      RatFunc (GaugeTargetParameterField K)) =
      gaugeFullyGenericPi := rfl

@[simp]
theorem generalGaugeFullyGenericSeedOverTarget_coeff
    (G : K[X]) (n : ℕ) :
    (generalGaugeFullyGenericSeedOverTarget G).coeff n =
      algebraMap K (RatFunc (GaugeTargetParameterField K)) (G.coeff n) := by
  simp [generalGaugeFullyGenericSeedOverTarget, generalGaugeGenericSeed,
    IsScalarTower.algebraMap_apply K (GaugeTargetParameterField K)
      (RatFunc (GaugeTargetParameterField K))]

theorem generalGaugeFullyGenericSeedOverTarget_map
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    (generalGaugeFullyGenericSeedOverTarget G).map
        (generalGaugeFullyGenericTargetHom G h₁ h₃) =
      G.map (algebraMap K (GaugeFunctionField K)) := by
  ext n
  simp [generalGaugeFullyGenericSeedOverTarget_coeff]

theorem generalGaugeSeedPolynomial_map_fieldHom
    {F L : Type*} [Field F] [Field L]
    (f : F →+* L) (G : F[X]) (pi : F) :
    (generalGaugeSeedPolynomial G pi).map f =
      generalGaugeSeedPolynomial (G.map f) (f pi) := by
  ext n
  simp [generalGaugeSeedPolynomial, Polynomial.natDegree_map f,
    Polynomial.coeff_map]
  apply Finset.sum_congr rfl
  intro x hx
  rw [← Polynomial.coeff_map]
  congr 1
  simp

theorem generalGaugeTargetPolynomial_ratFunc_eq_literal
    {F : Type*} [Field F] (G : F[X]) (pi b : F) :
    (generalGaugeTargetPolynomial G pi b).map
        (algebraMap F[X] (RatFunc F)) =
      generalGaugeInversePolynomial
        (G.map (algebraMap F (RatFunc F)))
        (algebraMap F (RatFunc F) pi)
        (algebraMap F (RatFunc F) b)
        RatFunc.X := by
  have hC :
      (algebraMap F[X] (RatFunc F)).comp (Polynomial.C : F →+* F[X]) =
        algebraMap F (RatFunc F) := by
    ext r
    simp
  have hCr (r : F) :
      algebraMap F[X] (RatFunc F) (Polynomial.C r) =
        algebraMap F (RatFunc F) r :=
    DFunLike.congr_fun hC r
  have hhalf :
      algebraMap F (RatFunc F) (G.coeff 1 / 2) =
        algebraMap F (RatFunc F) (G.coeff 1) / 2 := by
    simpa only [map_ofNat] using
      (map_div₀ (algebraMap F (RatFunc F)) (G.coeff 1) (2 : F))
  have hCrhalf :
      algebraMap F[X] (RatFunc F)
          (Polynomial.C (G.coeff 1 / 2)) =
        algebraMap F (RatFunc F) (G.coeff 1) / 2 :=
    (hCr (G.coeff 1 / 2)).trans hhalf
  rw [generalGaugeTargetPolynomial,
    GenericInverse.linearTargetPolynomial,
    generalGaugeTargetFreeInversePolynomial,
    generalGaugeInversePolynomial]
  simp only [Polynomial.map_sub, Polynomial.map_map, Polynomial.map_C,
    Polynomial.map_mul, Polynomial.map_pow, Polynomial.map_X]
  rw [hC]
  rw [generalGaugeSeedPolynomial_map_fieldHom]
  simp only [Polynomial.coeff_map]
  rw [hCr b, hCrhalf]
  simp only [map_mul, hCrhalf, RatFunc.algebraMap_X]
  ring

theorem generalGaugeFullyGenericInversePolynomial_eq_literal
    (G : K[X]) :
    generalGaugeFullyGenericInversePolynomial G =
      generalGaugeInversePolynomial
        (generalGaugeFullyGenericSeedOverTarget G)
        gaugeFullyGenericPi gaugeFullyGenericB RatFunc.X := by
  exact generalGaugeTargetPolynomial_ratFunc_eq_literal
    (generalGaugeGenericSeed G) (gaugeGenericPi K) (gaugeGenericB K)

/-- The three independent source variables inside the source function field. -/
def gaugeSourceVariable (i : Fin 3) : GaugeFunctionField K :=
  algebraMap (GaugePolynomial K) (GaugeFunctionField K) (MvPolynomial.X i)

theorem eval₂_gaugeSourceVariable (P : GaugePolynomial K) :
    MvPolynomial.eval₂ (algebraMap K (GaugeFunctionField K))
        gaugeSourceVariable P =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K) P := by
  let f : GaugePolynomial K →ₐ[K] GaugeFunctionField K :=
    MvPolynomial.aeval gaugeSourceVariable
  have hf :
      f = IsScalarTower.toAlgHom K (GaugePolynomial K)
        (GaugeFunctionField K) := by
    apply MvPolynomial.algHom_ext
    intro i
    simp [f, gaugeSourceVariable]
  exact DFunLike.congr_fun hf P

theorem eval₂_fullyGeneric_generalGaugeMap
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (i : Fin 3) :
    MvPolynomial.eval₂
        (generalGaugeFullyGenericTargetHom G h₁ h₃)
        gaugeSourceVariable
        (generalGaugeMap (generalGaugeFullyGenericSeedOverTarget G) i) =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K)
        (generalGaugeMap G i) := by
  calc
    MvPolynomial.eval₂
        (generalGaugeFullyGenericTargetHom G h₁ h₃)
        gaugeSourceVariable
        (generalGaugeMap (generalGaugeFullyGenericSeedOverTarget G) i) =
      MvPolynomial.eval₂ (RingHom.id (GaugeFunctionField K))
        gaugeSourceVariable
        (MvPolynomial.map
          (generalGaugeFullyGenericTargetHom G h₁ h₃)
          (generalGaugeMap (generalGaugeFullyGenericSeedOverTarget G) i)) := by
            symm
            exact MvPolynomial.eval₂_map
              (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom
              gaugeSourceVariable (RingHom.id (GaugeFunctionField K)) _
    _ = MvPolynomial.eval₂ (RingHom.id (GaugeFunctionField K))
        gaugeSourceVariable
        (generalGaugeMap
          ((generalGaugeFullyGenericSeedOverTarget G).map
            (generalGaugeFullyGenericTargetHom G h₁ h₃)) i) := by
          congr 1
          exact congrFun
            (generalGaugeMap_map
              (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom
              (generalGaugeFullyGenericSeedOverTarget G)) i
    _ = MvPolynomial.eval₂ (RingHom.id (GaugeFunctionField K))
        gaugeSourceVariable
        (generalGaugeMap
          (G.map (algebraMap K (GaugeFunctionField K))) i) := by
          rw [generalGaugeFullyGenericSeedOverTarget_map G h₁ h₃]
    _ = MvPolynomial.eval₂ (algebraMap K (GaugeFunctionField K))
        gaugeSourceVariable (generalGaugeMap G i) := by
          rw [← congrFun
            (generalGaugeMap_map (algebraMap K (GaugeFunctionField K)) G) i]
          exact MvPolynomial.eval₂_map
            (algebraMap K (GaugeFunctionField K))
            gaugeSourceVariable (RingHom.id (GaugeFunctionField K)) _
    _ = _ := eval₂_gaugeSourceVariable (generalGaugeMap G i)

/-- The actual generic source point, viewed as a point of the fully generic
displayed fiber over `K(Π,B)(C)`. -/
def generalGaugeFullyGenericRawPoint
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    letI : Algebra (RatFunc (GaugeTargetParameterField K))
        (GaugeFunctionField K) :=
      (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
    GeneralGaugeRawFiberPoint
      (generalGaugeFullyGenericSeedOverTarget G)
      gaugeFullyGenericPiUnit gaugeFullyGenericB RatFunc.X
      (GaugeFunctionField K) := by
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  refine
    { point := gaugeSourceVariable
      pi_eq := ?_
      b_eq := ?_
      c_eq := ?_ }
  · change MvPolynomial.eval₂
        (generalGaugeFullyGenericTargetHom G h₁ h₃)
        gaugeSourceVariable
        (generalGaugePi (generalGaugeFullyGenericSeedOverTarget G)) =
      generalGaugeFullyGenericTargetHom G h₁ h₃ gaugeFullyGenericPi
    simp only [gaugeFullyGenericPi]
    rw [generalGaugeFullyGenericTargetHom_pi]
    exact eval₂_fullyGeneric_generalGaugeMap G h₁ h₃ 0
  · change MvPolynomial.eval₂
        (generalGaugeFullyGenericTargetHom G h₁ h₃)
        gaugeSourceVariable
        (generalGaugeB (generalGaugeFullyGenericSeedOverTarget G)) =
      generalGaugeFullyGenericTargetHom G h₁ h₃ gaugeFullyGenericB
    simp only [gaugeFullyGenericB]
    rw [generalGaugeFullyGenericTargetHom_b]
    exact eval₂_fullyGeneric_generalGaugeMap G h₁ h₃ 1
  · change MvPolynomial.eval₂
        (generalGaugeFullyGenericTargetHom G h₁ h₃)
        gaugeSourceVariable
        (generalGaugeC (generalGaugeFullyGenericSeedOverTarget G)) =
      generalGaugeFullyGenericTargetHom G h₁ h₃ RatFunc.X
    rw [generalGaugeFullyGenericTargetHom_c]
    exact eval₂_fullyGeneric_generalGaugeMap G h₁ h₃ 2

theorem generalGaugeFullyGenericSeedOverTarget_coeff_one_ne_zero
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) :
    (generalGaugeFullyGenericSeedOverTarget G).coeff 1 ≠ 0 := by
  simpa using
    (algebraMap K (RatFunc (GaugeTargetParameterField K))).injective.ne h₁

theorem generalGaugeFullyGenericSeedOverTarget_coeff_three_ne_zero
    (G : K[X]) (h₃ : G.coeff 3 ≠ 0) :
    (generalGaugeFullyGenericSeedOverTarget G).coeff 3 ≠ 0 := by
  simpa using
    (algebraMap K (RatFunc (GaugeTargetParameterField K))).injective.ne h₃

theorem generalGaugeFullyGenericInversePolynomial_separable
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (hdeg : 3 ≤ G.natDegree) :
    (generalGaugeFullyGenericInversePolynomial G).Separable :=
  PerfectField.separable_of_irreducible
    (generalGaugeFullyGenericInversePolynomial_certificate G h₁ hdeg).1

/-- The fully generic inverse quotient represents the literal generic fiber,
with the pre-existing `K(Π,B)(C)` polynomial exposed in the source type. -/
def generalGaugeFullyGenericRawRepresentingEquiv
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree)
    (A : Type*) [CommRing A]
    [Algebra (RatFunc (GaugeTargetParameterField K)) A] :
    (AdjoinRoot
        (generalGaugeInversePolynomial
          (generalGaugeFullyGenericSeedOverTarget G)
          gaugeFullyGenericPi gaugeFullyGenericB RatFunc.X) →ₐ[
            RatFunc (GaugeTargetParameterField K)] A) ≃
      GeneralGaugeRawFiberPoint
        (generalGaugeFullyGenericSeedOverTarget G)
        gaugeFullyGenericPiUnit gaugeFullyGenericB RatFunc.X A := by
  exact generalGaugeRawRepresentingEquiv
    (generalGaugeFullyGenericSeedOverTarget G)
    gaugeFullyGenericPiUnit gaugeFullyGenericB RatFunc.X
    (generalGaugeFullyGenericSeedOverTarget_coeff_one_ne_zero G h₁)
    (generalGaugeFullyGenericSeedOverTarget_coeff_three_ne_zero G h₃)
    (by
      simp only [gaugeFullyGenericPiUnit_val]
      rw [← generalGaugeFullyGenericInversePolynomial_eq_literal]
      exact generalGaugeFullyGenericInversePolynomial_separable G h₁ hdeg)
    A

/-- Naturality of the fully generic represented-fiber equivalence. -/
theorem generalGaugeFullyGenericRawRepresentingEquiv_natural
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree)
    {A B : Type*} [CommRing A] [CommRing B]
    [Algebra (RatFunc (GaugeTargetParameterField K)) A]
    [Algebra (RatFunc (GaugeTargetParameterField K)) B]
    (f : A →ₐ[RatFunc (GaugeTargetParameterField K)] B)
    (φ : AdjoinRoot
      (generalGaugeInversePolynomial
        (generalGaugeFullyGenericSeedOverTarget G)
        gaugeFullyGenericPi gaugeFullyGenericB RatFunc.X) →ₐ[
          RatFunc (GaugeTargetParameterField K)] A) :
    GeneralGaugeRawFiberPoint.map f
        (generalGaugeFullyGenericRawRepresentingEquiv
          G h₁ h₃ hdeg A φ) =
      generalGaugeFullyGenericRawRepresentingEquiv
        G h₁ h₃ hdeg B (f.comp φ) := by
  exact generalGaugeRawRepresentingEquiv_natural
    (generalGaugeFullyGenericSeedOverTarget G)
    gaugeFullyGenericPiUnit gaugeFullyGenericB RatFunc.X
    (generalGaugeFullyGenericSeedOverTarget_coeff_one_ne_zero G h₁)
    (generalGaugeFullyGenericSeedOverTarget_coeff_three_ne_zero G h₃)
    (by
      simp only [gaugeFullyGenericPiUnit_val]
      rw [← generalGaugeFullyGenericInversePolynomial_eq_literal]
      exact generalGaugeFullyGenericInversePolynomial_separable G h₁ hdeg)
    f φ

/-- The comparison homomorphism in the literal fully generic presentation. -/
def generalGaugeLiteralFunctionFieldComparisonHom
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    letI : Algebra (RatFunc (GaugeTargetParameterField K))
        (GaugeFunctionField K) :=
      (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
    AdjoinRoot
        (generalGaugeInversePolynomial
          (generalGaugeFullyGenericSeedOverTarget G)
          gaugeFullyGenericPi gaugeFullyGenericB RatFunc.X) →ₐ[
            RatFunc (GaugeTargetParameterField K)] GaugeFunctionField K := by
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  exact (generalGaugeFullyGenericRawRepresentingEquiv
    G h₁ h₃ hdeg (GaugeFunctionField K)).symm
      (generalGaugeFullyGenericRawPoint G h₁ h₃)

/-- The comparison homomorphism from the already formalized inverse-root
extension to the actual source function field. -/
def generalGaugeFunctionFieldComparisonHom
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    letI : Algebra (RatFunc (GaugeTargetParameterField K))
        (GaugeFunctionField K) :=
      (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
    AdjoinRoot (generalGaugeFullyGenericInversePolynomial G) →ₐ[
      RatFunc (GaugeTargetParameterField K)] GaugeFunctionField K := by
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  exact (generalGaugeLiteralFunctionFieldComparisonHom
      G h₁ h₃ hdeg).comp
    (AdjoinRoot.algEquivOfEq
      (RatFunc (GaugeTargetParameterField K))
      (generalGaugeFullyGenericInversePolynomial G)
      (generalGaugeInversePolynomial
        (generalGaugeFullyGenericSeedOverTarget G)
        gaugeFullyGenericPi gaugeFullyGenericB RatFunc.X)
      (generalGaugeFullyGenericInversePolynomial_eq_literal G)).toAlgHom

theorem generalGaugeLiteralFunctionFieldComparisonHom_surjective
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    letI : Algebra (RatFunc (GaugeTargetParameterField K))
        (GaugeFunctionField K) :=
      (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
    Function.Surjective
      (generalGaugeLiteralFunctionFieldComparisonHom G h₁ h₃ hdeg) := by
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  let E :=
    generalGaugeInversePolynomial
      (generalGaugeFullyGenericSeedOverTarget G)
      gaugeFullyGenericPi gaugeFullyGenericB
      (RatFunc.X : RatFunc (GaugeTargetParameterField K))
  letI : Fact (Irreducible E) :=
    ⟨by
      dsimp [E]
      rw [← generalGaugeFullyGenericInversePolynomial_eq_literal]
      exact
        (generalGaugeFullyGenericInversePolynomial_certificate
          G h₁ hdeg).1⟩
  let Q := AdjoinRoot E
  let Φ : Q →ₐ[RatFunc (GaugeTargetParameterField K)]
      GaugeFunctionField K :=
    generalGaugeLiteralFunctionFieldComparisonHom G h₁ h₃ hdeg
  let U : GeneralGaugeRawFiberPoint
      (generalGaugeFullyGenericSeedOverTarget G)
      gaugeFullyGenericPiUnit gaugeFullyGenericB RatFunc.X Q :=
    generalGaugeFullyGenericRawRepresentingEquiv
      G h₁ h₃ hdeg Q
        (AlgHom.id (RatFunc (GaugeTargetParameterField K)) Q)
  have hcomp :
      Φ.comp (AlgHom.id (RatFunc (GaugeTargetParameterField K)) Q) = Φ := by
    ext x
    rfl
  have hnat :=
    generalGaugeFullyGenericRawRepresentingEquiv_natural
      G h₁ h₃ hdeg Φ
        (AlgHom.id (RatFunc (GaugeTargetParameterField K)) Q)
  rw [hcomp] at hnat
  have happ :
      generalGaugeFullyGenericRawRepresentingEquiv
          G h₁ h₃ hdeg (GaugeFunctionField K) Φ =
        generalGaugeFullyGenericRawPoint G h₁ h₃ :=
    (generalGaugeFullyGenericRawRepresentingEquiv
      G h₁ h₃ hdeg (GaugeFunctionField K)).apply_symm_apply _
  rw [happ] at hnat
  have hcoord (i : Fin 3) : Φ (U.point i) = gaugeSourceVariable i := by
    exact congrArg (fun p => p.point i) hnat
  intro z
  obtain ⟨p, q, hq, hz⟩ :=
    IsFractionRing.div_surjective (GaugePolynomial K) z
  let p' : Q :=
    MvPolynomial.eval₂ (algebraMap K Q) U.point p
  let q' : Q :=
    MvPolynomial.eval₂ (algebraMap K Q) U.point q
  have hbase :
      Φ.toRingHom.comp (algebraMap K Q) =
        algebraMap K (GaugeFunctionField K) := by
    ext r
    calc
      Φ (algebraMap K Q r) =
          Φ (algebraMap (RatFunc (GaugeTargetParameterField K)) Q
            (algebraMap K (RatFunc (GaugeTargetParameterField K)) r)) := by
              rw [IsScalarTower.algebraMap_apply K
                (RatFunc (GaugeTargetParameterField K)) Q]
      _ = algebraMap (RatFunc (GaugeTargetParameterField K))
          (GaugeFunctionField K)
          (algebraMap K (RatFunc (GaugeTargetParameterField K)) r) :=
            Φ.commutes _
      _ = generalGaugeFullyGenericTargetHom G h₁ h₃
          (algebraMap K (RatFunc (GaugeTargetParameterField K)) r) := rfl
      _ = algebraMap K (GaugeFunctionField K) r :=
        (generalGaugeFullyGenericTargetHom G h₁ h₃).commutes r
  have heval (P : GaugePolynomial K) :
      Φ (MvPolynomial.eval₂ (algebraMap K Q) U.point P) =
        algebraMap (GaugePolynomial K) (GaugeFunctionField K) P := by
    calc
      Φ (MvPolynomial.eval₂ (algebraMap K Q) U.point P) =
          MvPolynomial.eval₂
            (Φ.toRingHom.comp (algebraMap K Q))
            (fun i => Φ (U.point i)) P := by
              simpa using MvPolynomial.hom_eval₂ P
                (algebraMap K Q) Φ.toRingHom U.point
      _ = MvPolynomial.eval₂
            (algebraMap K (GaugeFunctionField K))
            gaugeSourceVariable P := by
              rw [hbase]
              congr
              funext i
              exact hcoord i
      _ = _ := eval₂_gaugeSourceVariable P
  have hq' : q' ≠ 0 := by
    intro hzero
    have himage : algebraMap (GaugePolynomial K)
        (GaugeFunctionField K) q = 0 := by
      calc
        algebraMap (GaugePolynomial K) (GaugeFunctionField K) q =
            Φ q' := (heval q).symm
        _ = Φ 0 := congrArg Φ hzero
        _ = 0 := map_zero Φ
    have qzero : q = 0 := by
      apply IsFractionRing.injective (GaugePolynomial K)
        (GaugeFunctionField K)
      simpa using himage
    exact nonZeroDivisors.ne_zero hq qzero
  refine ⟨p' / q', ?_⟩
  rw [map_div₀, show Φ p' =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K) p from heval p,
    show Φ q' =
      algebraMap (GaugePolynomial K) (GaugeFunctionField K) q from heval q]
  exact hz

/-- Explicit comparison
`K(Π,B)(C)[S]/(E) ≃ K(x,y,z)` over the fully generic target field. -/
def generalGaugeFunctionFieldComparison
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    letI : Algebra (RatFunc (GaugeTargetParameterField K))
        (GaugeFunctionField K) :=
      (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
    AdjoinRoot (generalGaugeFullyGenericInversePolynomial G) ≃ₐ[
      RatFunc (GaugeTargetParameterField K)] GaugeFunctionField K := by
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  let E :=
    generalGaugeInversePolynomial
      (generalGaugeFullyGenericSeedOverTarget G)
      gaugeFullyGenericPi gaugeFullyGenericB
      (RatFunc.X : RatFunc (GaugeTargetParameterField K))
  letI : Fact (Irreducible E) :=
    ⟨by
      dsimp [E]
      rw [← generalGaugeFullyGenericInversePolynomial_eq_literal]
      exact
        (generalGaugeFullyGenericInversePolynomial_certificate
          G h₁ hdeg).1⟩
  let literalEquiv :
      AdjoinRoot E ≃ₐ[RatFunc (GaugeTargetParameterField K)]
        GaugeFunctionField K :=
    AlgEquiv.ofBijective
      (generalGaugeLiteralFunctionFieldComparisonHom G h₁ h₃ hdeg)
      ⟨(generalGaugeLiteralFunctionFieldComparisonHom
          G h₁ h₃ hdeg).injective,
        generalGaugeLiteralFunctionFieldComparisonHom_surjective
          G h₁ h₃ hdeg⟩
  exact (AdjoinRoot.algEquivOfEq
      (RatFunc (GaugeTargetParameterField K))
      (generalGaugeFullyGenericInversePolynomial G) E
      (generalGaugeFullyGenericInversePolynomial_eq_literal G)).trans
    literalEquiv

/-- The same comparison in the geometric direction requested by the
function-field statement:
`K(x,y,z) ≃ K(Π,B)(C)[S]/(E_{Π,B,C})`. -/
def generalGaugeSourceFunctionFieldComparison
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    letI : Algebra (RatFunc (GaugeTargetParameterField K))
        (GaugeFunctionField K) :=
      (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
    GaugeFunctionField K ≃ₐ[RatFunc (GaugeTargetParameterField K)]
      AdjoinRoot (generalGaugeFullyGenericInversePolynomial G) := by
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  exact (generalGaugeFunctionFieldComparison G h₁ h₃ hdeg).symm

/-- The geometric degree of the displayed gauge map, computed using the
canonical target presentation `K(Π,B)(C)`. -/
def generalGaugeGeometricDegree
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) : ℕ :=
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  Module.finrank (RatFunc (GaugeTargetParameterField K))
    (GaugeFunctionField K)

/-- The actual geometric degree of the displayed map is the seed degree. -/
theorem generalGaugeGeometricDegree_eq
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    generalGaugeGeometricDegree G h₁ h₃ = G.natDegree := by
  letI : Algebra (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) :=
    (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
  change Module.finrank (RatFunc (GaugeTargetParameterField K))
      (GaugeFunctionField K) = G.natDegree
  rw [← (generalGaugeFunctionFieldComparison G h₁ h₃ hdeg).toLinearEquiv.finrank_eq]
  exact generalGaugeFullyGenericInverseAdjoinRoot_finrank G hdeg

end FiniteEtaleKeller
