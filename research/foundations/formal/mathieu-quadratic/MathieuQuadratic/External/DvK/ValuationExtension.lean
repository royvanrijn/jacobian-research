/-
Adapted from GMC2/ValuationExtension.lean at MurrellGroup/GMC-2-lean commit
1782de7ff6c97eb1d98e63e7ff34df18b9cd322e. MIT license in this directory.
-/
import Mathlib.RingTheory.DedekindDomain.IntegralClosure
import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.RingTheory.Ideal.GoingUp
import Mathlib.RingTheory.Valuation.AlgebraInstances
import Mathlib.RingTheory.Valuation.Discrete.Basic
import Mathlib.RingTheory.Valuation.Extension

open scoped WithZero

namespace MathieuQuadratic.DvKProof

open IsDedekindDomain MonoidWithZeroHom

section

variable {K L : Type*} [Field K] [Field L] [Algebra K L]

/--
A chosen prime of the integral closure above the maximal ideal of a discrete
valuation ring. This packages the lying-over choice used to extend a discrete
valuation to a finite separable field extension.
-/
noncomputable def primeAboveValuation
    (v : Valuation K ℤᵐ⁰)
    [IsCyclic (valueGroup (.ofClass v))]
    [Nontrivial (valueGroup (.ofClass v))]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] :
    HeightOneSpectrum (integralClosure v.valuationSubring L) := by
  let A := v.valuationSubring
  let B := integralClosure A L
  let p := IsLocalRing.maximalIdeal A
  haveI : IsDiscreteValuationRing A := inferInstance
  haveI : IsDedekindDomain A := inferInstance
  haveI : IsFractionRing A K :=
    (Valuation.valuationSubring.integers v).isFractionRing
  haveI : IsDedekindDomain B :=
    integralClosure.isDedekindDomain A K L
  haveI : Algebra.IsIntegral A B :=
    IsIntegralClosure.isIntegral_algebra A L
  haveI : FaithfulSMul A B :=
    (faithfulSMul_iff_algebraMap_injective ..).2
      (ValuationSubring.integralClosure_algebraMap_injective v L)
  let hPexists :=
    Ideal.exists_maximal_ideal_liesOver_of_isIntegral (S := B) p
  let P := Classical.choose hPexists
  have hPmax : P.IsMaximal := (Classical.choose_spec hPexists).1
  have hPover : P.LiesOver p := (Classical.choose_spec hPexists).2
  exact
    { asIdeal := P
      isPrime := hPmax.isPrime
      ne_bot := by
        intro hP
        have hpbot : p = ⊥ := by
          calc
            p = P.under A := hPover.over
            _ = (⊥ : Ideal B).under A := by rw [hP]
            _ = ⊥ := by
              rw [Ideal.under_def]
              exact Ideal.comap_bot_of_injective _
                (FaithfulSMul.algebraMap_injective A B)
        exact IsDiscreteValuationRing.not_a_field A hpbot }

theorem primeAboveValuation_liesOver
    (v : Valuation K ℤᵐ⁰)
    [IsCyclic (valueGroup (.ofClass v))]
    [Nontrivial (valueGroup (.ofClass v))]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] :
    (primeAboveValuation (L := L) v).asIdeal.LiesOver
      (IsLocalRing.maximalIdeal v.valuationSubring) := by
  let A := v.valuationSubring
  let B := integralClosure A L
  let p := IsLocalRing.maximalIdeal A
  haveI : IsDiscreteValuationRing A := inferInstance
  haveI : IsDedekindDomain A := inferInstance
  haveI : IsFractionRing A K :=
    (Valuation.valuationSubring.integers v).isFractionRing
  haveI : IsDedekindDomain B :=
    integralClosure.isDedekindDomain A K L
  haveI : Algebra.IsIntegral A B :=
    IsIntegralClosure.isIntegral_algebra A L
  haveI : FaithfulSMul A B :=
    (faithfulSMul_iff_algebraMap_injective ..).2
      (ValuationSubring.integralClosure_algebraMap_injective v L)
  let hPexists :=
    Ideal.exists_maximal_ideal_liesOver_of_isIntegral (S := B) p
  simpa only [primeAboveValuation] using
    (Classical.choose_spec hPexists).2

/--
An extension of a rank-one discrete valuation to a finite separable field
extension, obtained from a prime of the integral closure.
-/
noncomputable def extendedDiscreteValuation
    (v : Valuation K ℤᵐ⁰)
    [IsCyclic (valueGroup (.ofClass v))]
    [Nontrivial (valueGroup (.ofClass v))]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] :
    Valuation L ℤᵐ⁰ := by
  let A := v.valuationSubring
  let B := integralClosure A L
  haveI : IsDiscreteValuationRing A := inferInstance
  haveI : IsDedekindDomain A := inferInstance
  haveI : IsFractionRing A K :=
    (Valuation.valuationSubring.integers v).isFractionRing
  haveI : IsDedekindDomain B :=
    integralClosure.isDedekindDomain A K L
  haveI : IsFractionRing B L :=
    IsIntegralClosure.isFractionRing_of_finite_extension A K L B
  exact (primeAboveValuation (L := L) v).valuation L

theorem extendedDiscreteValuation_le_one_iff
    (v : Valuation K ℤᵐ⁰)
    [IsCyclic (valueGroup (.ofClass v))]
    [Nontrivial (valueGroup (.ofClass v))]
    [FiniteDimensional K L] [Algebra.IsSeparable K L]
    (x : K) :
    extendedDiscreteValuation (L := L) v (algebraMap K L x) ≤ 1 ↔
      v x ≤ 1 := by
  let A := v.valuationSubring
  let B := integralClosure A L
  let p := IsLocalRing.maximalIdeal A
  let P := primeAboveValuation (L := L) v
  haveI : IsDiscreteValuationRing A := inferInstance
  haveI : IsDedekindDomain A := inferInstance
  haveI : IsFractionRing A K :=
    (Valuation.valuationSubring.integers v).isFractionRing
  haveI : IsDedekindDomain B :=
    integralClosure.isDedekindDomain A K L
  haveI : IsFractionRing B L :=
    IsIntegralClosure.isFractionRing_of_finite_extension A K L B
  letI : P.asIdeal.LiesOver p :=
    primeAboveValuation_liesOver (L := L) v
  change P.valuation L (algebraMap K L x) ≤ 1 ↔ v x ≤ 1
  constructor
  · intro hxL
    by_contra hx
    have hxne : x ≠ 0 := by
      intro hzero
      apply hx
      simp [hzero]
    have hxgt : 1 < v x := lt_of_not_ge hx
    have hxinvalt : v x⁻¹ < 1 :=
      (v.one_lt_val_iff hxne).mp hxgt
    have hxinvA : x⁻¹ ∈ A := le_of_lt hxinvalt
    let a : A := ⟨x⁻¹, hxinvA⟩
    have haNonunit : ¬IsUnit a := by
      intro ha
      let y : A := ↑(ha.unit⁻¹)
      have hy : (y : K) = x := by
        calc
          (y : K) = ((ha.unit : A) : K)⁻¹ := by
            change
              algebraMap A K (↑ha.unit⁻¹) =
                (algebraMap A K (↑ha.unit))⁻¹
            have hu :
                algebraMap A K (↑ha.unit) ≠ 0 :=
              (map_ne_zero_iff (algebraMap A K)
                (FaithfulSMul.algebraMap_injective A K)).mpr
                (Units.ne_zero ha.unit)
            apply (mul_eq_one_iff_eq_inv₀ hu).mp
            rw [← map_mul]
            simp
          _ = (a : K)⁻¹ := by rw [ha.unit_spec]
          _ = x := by simp [a]
      apply hx
      have hyMem : v (y : K) ≤ 1 := y.property
      simpa [hy] using hyMem
    have haMax : a ∈ p := by
      rw [IsLocalRing.mem_maximalIdeal]
      exact haNonunit
    have haP : algebraMap A B a ∈ P.asIdeal :=
      (Ideal.mem_of_liesOver P.asIdeal p a).mp haMax
    have hinvLt :
        P.valuation L
            (algebraMap B L (algebraMap A B a)) < 1 :=
      (P.valuation_lt_one_iff_mem (K := L)
        (algebraMap A B a)).mpr haP
    have hmapa :
        algebraMap B L (algebraMap A B a) =
          (algebraMap K L x)⁻¹ := by
      calc
        algebraMap B L (algebraMap A B a) =
            algebraMap A L a := by
              rw [IsScalarTower.algebraMap_apply A B L]
        _ = algebraMap K L (algebraMap A K a) := by
              rw [IsScalarTower.algebraMap_apply A K L]
        _ = (algebraMap K L x)⁻¹ := by
              simp [a]
    rw [hmapa] at hinvLt
    have hmapne : algebraMap K L x ≠ 0 :=
      (map_ne_zero (algebraMap K L)).mpr hxne
    have hmapgt :
        1 < P.valuation L (algebraMap K L x) :=
      ((P.valuation L).one_lt_val_iff hmapne).mpr hinvLt
    exact (not_lt_of_ge hxL) hmapgt
  · intro hx
    have hxA : x ∈ A := hx
    let a : A := ⟨x, hxA⟩
    let b : B := algebraMap A B a
    have hb :
        P.valuation L (algebraMap B L b) ≤ 1 :=
      P.valuation_le_one (K := L) b
    have hmapb :
        algebraMap B L b = algebraMap K L x := by
      calc
        algebraMap B L b = algebraMap A L a := by
          dsimp [b]
        _ = algebraMap K L (algebraMap A K a) := by
          rw [IsScalarTower.algebraMap_apply A K L]
        _ = algebraMap K L x := rfl
    rwa [hmapb] at hb

theorem extendedDiscreteValuation_hasExtension
    (v : Valuation K ℤᵐ⁰)
    [IsCyclic (valueGroup (.ofClass v))]
    [Nontrivial (valueGroup (.ofClass v))]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] :
    v.HasExtension (extendedDiscreteValuation (L := L) v) where
  val_isEquiv_comap := by
    rw [Valuation.isEquiv_iff_val_le_one]
    intro x
    exact (extendedDiscreteValuation_le_one_iff (L := L) v x).symm

end

end MathieuQuadratic.DvKProof
