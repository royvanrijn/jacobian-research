/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeFunctionField
import FiniteEtaleKeller.GenericInverseIrreducibility

/-!
# The fully generic inverse extension

The fixed-parameter inverse theorem can be applied over the two-variable
rational function field `K(Π, B)`.  Taking `Π` and `B` to be its independent
coordinate functions and adjoining the remaining rational parameter `C`
produces the inverse equation over the iterated rational function field
`K(Π, B)(C)`.

This module proves that equation irreducible and computes the exact finrank of
its root quotient.  The further comparison between this iterated presentation
of the target function field and the image of `generalGaugeFunctionFieldHom`
is deliberately not built into these statements.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable (K : Type*) [Field K]

/-- A polynomial ring in the independent target parameters `Π` and `B`. -/
abbrev GaugeTargetParameterPolynomial :=
  MvPolynomial (Fin 2) K

/-- The rational function field `K(Π, B)`. -/
abbrev GaugeTargetParameterField :=
  FractionRing (GaugeTargetParameterPolynomial K)

/-- The independent target parameter `Π` in `K(Π, B)`. -/
def gaugeGenericPi : GaugeTargetParameterField K :=
  algebraMap (GaugeTargetParameterPolynomial K)
    (GaugeTargetParameterField K) (MvPolynomial.X 0)

/-- The independent target parameter `B` in `K(Π, B)`. -/
def gaugeGenericB : GaugeTargetParameterField K :=
  algebraMap (GaugeTargetParameterPolynomial K)
    (GaugeTargetParameterField K) (MvPolynomial.X 1)

theorem gaugeGenericPi_ne_zero :
    gaugeGenericPi K ≠ 0 := by
  unfold gaugeGenericPi
  simpa only [map_zero] using
    (IsFractionRing.injective
        (GaugeTargetParameterPolynomial K)
        (GaugeTargetParameterField K)).ne
      (MvPolynomial.X_ne_zero 0)

variable {K}

/-- Extend a seed from `K` to the independent-parameter field `K(Π, B)`. -/
def generalGaugeGenericSeed (G : K[X]) :
    (GaugeTargetParameterField K)[X] :=
  G.map (algebraMap K (GaugeTargetParameterField K))

@[simp]
theorem generalGaugeGenericSeed_coeff (G : K[X]) (n : ℕ) :
    (generalGaugeGenericSeed G).coeff n =
      algebraMap K (GaugeTargetParameterField K) (G.coeff n) := by
  simp [generalGaugeGenericSeed]

theorem generalGaugeGenericSeed_natDegree (G : K[X]) :
    (generalGaugeGenericSeed G).natDegree = G.natDegree := by
  exact natDegree_map_eq_of_injective
    (algebraMap K (GaugeTargetParameterField K)).injective G

/-- The inverse equation over `K(Π, B)(C)`. -/
def generalGaugeFullyGenericInversePolynomial (G : K[X]) :
    (RatFunc (GaugeTargetParameterField K))[X] :=
  generalGaugeGenericInversePolynomial
    (generalGaugeGenericSeed G)
    (gaugeGenericPi K)
    (gaugeGenericB K)

/-- The full three-parameter inverse equation is irreducible and has exactly
the seed degree. -/
theorem generalGaugeFullyGenericInversePolynomial_certificate
    [CharZero K] (G : K[X])
    (h₁ : G.coeff 1 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    Irreducible (generalGaugeFullyGenericInversePolynomial G) ∧
      (generalGaugeFullyGenericInversePolynomial G).natDegree =
        G.natDegree := by
  have h₁' : (generalGaugeGenericSeed G).coeff 1 ≠ 0 := by
    simpa using
      (algebraMap K (GaugeTargetParameterField K)).injective.ne h₁
  have hdeg' : 3 ≤ (generalGaugeGenericSeed G).natDegree := by
    simpa [generalGaugeGenericSeed_natDegree] using hdeg
  simpa [generalGaugeFullyGenericInversePolynomial,
    generalGaugeGenericSeed_natDegree] using
    (generalGaugeGenericInversePolynomial_certificate
      (generalGaugeGenericSeed G)
      (gaugeGenericPi K)
      (gaugeGenericB K)
      h₁' hdeg' (gaugeGenericPi_ne_zero K))

/-- The root quotient of the fully generic inverse equation has extension
dimension exactly the seed degree over `K(Π, B)(C)`. -/
theorem generalGaugeFullyGenericInverseAdjoinRoot_finrank
    (G : K[X]) (hdeg : 3 ≤ G.natDegree) :
    Module.finrank (RatFunc (GaugeTargetParameterField K))
        (AdjoinRoot (generalGaugeFullyGenericInversePolynomial G)) =
      G.natDegree := by
  have hdeg' : 3 ≤ (generalGaugeGenericSeed G).natDegree := by
    simpa [generalGaugeGenericSeed_natDegree] using hdeg
  change
    Module.finrank (RatFunc (GaugeTargetParameterField K))
        (AdjoinRoot
          (generalGaugeGenericInversePolynomial
            (generalGaugeGenericSeed G)
            (gaugeGenericPi K)
            (gaugeGenericB K))) =
      G.natDegree
  calc
    _ = (generalGaugeGenericSeed G).natDegree :=
      generalGaugeGenericInverseAdjoinRoot_finrank
        (generalGaugeGenericSeed G)
        (gaugeGenericPi K)
        (gaugeGenericB K)
        hdeg' (gaugeGenericPi_ne_zero K)
    _ = G.natDegree := generalGaugeGenericSeed_natDegree G

#print axioms generalGaugeFullyGenericInversePolynomial_certificate
#print axioms generalGaugeFullyGenericInverseAdjoinRoot_finrank

end FiniteEtaleKeller
