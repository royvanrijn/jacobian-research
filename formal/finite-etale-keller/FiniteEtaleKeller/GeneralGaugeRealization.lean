/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeNormalization
import FiniteEtaleKeller.AutomaticRealization

/-!
# The actual determinant-one realization map for a prescribed polynomial

This module specializes the literal all-degree map to the rooted translation
`G(S)=P(a+S)-P(a)`.  It proves that the generic inverse polynomial at the
chosen target is exactly `P(a+S)`, transports the quotient back to `K[T]/(P)`,
and obtains a natural equivalence with the literal fiber of the determinant-one
`MvPolynomial` map.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K] [CharZero K]

/-- At first target coordinate one, the coefficientwise generic seed
reconstructs every polynomial with zero constant term. -/
theorem generalGaugeSeedPolynomial_one_eq
    (G : K[X]) (hzero : G.coeff 0 = 0) :
    generalGaugeSeedPolynomial G 1 = G := by
  ext n
  by_cases h0 : n = 0
  · subst n
    simp [generalGaugeSeedPolynomial, hzero]
  by_cases h1 : n = 1
  · subst n
    simp [generalGaugeSeedPolynomial]
  by_cases h2 : n = 2
  · subst n
    simp [generalGaugeSeedPolynomial]
  by_cases h3 : n = 3
  · subst n
    simp [generalGaugeSeedPolynomial]
  have hn4 : 4 ≤ n := by omega
  by_cases hnN : n ≤ G.natDegree
  · have hnmem : n ∈ Finset.Icc 4 G.natDegree := Finset.mem_Icc.mpr ⟨hn4, hnN⟩
    simp [generalGaugeSeedPolynomial, h0, h1, h2, h3, hnmem,
      Finset.sum_eq_single n]
  · have hcoeff : G.coeff n = 0 := by
      by_contra hne
      have hle := Polynomial.le_natDegree_of_ne_zero hne
      omega
    have hnmem : n ∉ Finset.Icc 4 G.natDegree := by
      simp [Finset.mem_Icc, hn4, hnN]
    simp [generalGaugeSeedPolynomial, h0, h1, h2, h3, hnmem, hcoeff]

/-- For the rooted translated seed and the paper's chosen third target, the
literal generic inverse polynomial is exactly `P(a+S)`. -/
theorem generalGaugeInversePolynomial_realization
    (P : K[X]) (a : K) (h₁ : P.derivative.eval a ≠ 0) :
    generalGaugeInversePolynomial
        (rootedTranslate P a) 1 0
        (realizationTargetC P a (P.derivative.eval a)) =
      translatePolynomial P a := by
  rw [generalGaugeInversePolynomial]
  rw [generalGaugeSeedPolynomial_one_eq
    (rootedTranslate P a) (rootedTranslate_coeff_zero P a)]
  simpa using
    (rootedTranslate_inverse_at_target P a (P.derivative.eval a) h₁)

section SuppliedTranslation

variable {A B : Type*}
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- The actual rooted seed used by the displayed map. -/
def realizationSeed (P : K[X]) (a : K) : K[X] :=
  rootedTranslate P a

/-- The literal determinant-one map fiber realizing `K[T]/(P)`, for a supplied
admissible translation parameter. -/
def realizationJacobianOneFiberRepresentingEquiv
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (AdjoinRoot P →ₐ[K] A) ≃
      GeneralGaugeJacobianOneFiberPoint
        (realizationSeed P a) 1
        (realizationTargetC P a (P.derivative.eval a)) A := by
  let G := realizationSeed P a
  let c := realizationTargetC P a (P.derivative.eval a)
  have hG₁ : G.coeff 1 ≠ 0 := by
    simpa [G, realizationSeed] using rootedTranslate_linear_ne_zero P a h₁
  have hG₃ : G.coeff 3 ≠ 0 := by
    simpa [G, realizationSeed] using rootedTranslate_cubic_ne_zero P a h₃
  have hinv : generalGaugeInversePolynomial G 1 0 c = translatePolynomial P a := by
    simpa [G, c, realizationSeed] using
      generalGaugeInversePolynomial_realization P a h₁
  have hsep : (generalGaugeInversePolynomial G 1 0 c).Separable := by
    rw [hinv]
    exact translatePolynomial_separable P a
      ((PerfectField.separable_iff_squarefree).2 hP)
  refine (translatedQuotientHomEquiv (A := A) P a).trans ?_
  rw [← hinv]
  exact generalGaugeJacobianOneRepresentingEquiv G 1 c hG₁ hG₃ hsep A

/-- Naturality of the supplied-translation literal determinant-one
realization. -/
theorem realizationJacobianOneFiberRepresentingEquiv_natural
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0)
    (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A) :
    GeneralGaugeJacobianOneFiberPoint.map f
        (realizationJacobianOneFiberRepresentingEquiv
          (A := A) P a hP h₁ h₃ φ) =
      realizationJacobianOneFiberRepresentingEquiv
        (A := B) P a hP h₁ h₃ (f.comp φ) := by
  let G := realizationSeed P a
  let c := realizationTargetC P a (P.derivative.eval a)
  have hG₁ : G.coeff 1 ≠ 0 := by
    simpa [G, realizationSeed] using rootedTranslate_linear_ne_zero P a h₁
  have hG₃ : G.coeff 3 ≠ 0 := by
    simpa [G, realizationSeed] using rootedTranslate_cubic_ne_zero P a h₃
  have hinv : generalGaugeInversePolynomial G 1 0 c = translatePolynomial P a := by
    simpa [G, c, realizationSeed] using
      generalGaugeInversePolynomial_realization P a h₁
  have hsep : (generalGaugeInversePolynomial G 1 0 c).Separable := by
    rw [hinv]
    exact translatePolynomial_separable P a
      ((PerfectField.separable_iff_squarefree).2 hP)
  change GeneralGaugeJacobianOneFiberPoint.map f
      (generalGaugeJacobianOneRepresentingEquiv G 1 c hG₁ hG₃ hsep A
        (translatedQuotientHomEquiv (A := A) P a φ)) =
    generalGaugeJacobianOneRepresentingEquiv G 1 c hG₁ hG₃ hsep B
      (translatedQuotientHomEquiv (A := B) P a (f.comp φ))
  rw [generalGaugeJacobianOneRepresentingEquiv_natural]
  congr 1
  apply DFunLike.ext _ _
  intro x
  rfl

end SuppliedTranslation

section Automatic

variable {A B : Type*}
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- The actual determinant-one polynomial map chosen automatically from a
squarefree polynomial of degree at least three. -/
def automaticRealizationMap
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    Fin 3 → GaugePolynomial K :=
  generalGaugeJacobianOneMap
    (realizationSeed P (chosenAdmissibleTranslation P hdeg))

/-- The target of the automatically chosen actual realization map. -/
def automaticRealizationTargetC
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) : K :=
  realizationTargetC P (chosenAdmissibleTranslation P hdeg)
    (P.derivative.eval (chosenAdmissibleTranslation P hdeg))

/-- Every squarefree polynomial of degree at least three naturally represents
the literal fiber of an actual determinant-one `MvPolynomial` map.  No
translation parameter or nonvanishing witness remains as an input. -/
def automaticJacobianOneFiberRepresentingEquiv
    (P : K[X]) (hP : Squarefree P) (hdeg : 3 ≤ P.natDegree) :
    (AdjoinRoot P →ₐ[K] A) ≃
      GeneralGaugeJacobianOneFiberPoint
        (realizationSeed P (chosenAdmissibleTranslation P hdeg)) 1
        (automaticRealizationTargetC P hdeg) A :=
  realizationJacobianOneFiberRepresentingEquiv
    (A := A) P (chosenAdmissibleTranslation P hdeg) hP
    (chosenAdmissibleTranslation_linear_ne_zero P hdeg)
    (chosenAdmissibleTranslation_cubic_ne_zero P hdeg)

/-- Naturality of the final automatic literal determinant-one realization. -/
theorem automaticJacobianOneFiberRepresentingEquiv_natural
    (P : K[X]) (hP : Squarefree P) (hdeg : 3 ≤ P.natDegree)
    (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A) :
    GeneralGaugeJacobianOneFiberPoint.map f
        (automaticJacobianOneFiberRepresentingEquiv
          (A := A) P hP hdeg φ) =
      automaticJacobianOneFiberRepresentingEquiv
        (A := B) P hP hdeg (f.comp φ) := by
  exact realizationJacobianOneFiberRepresentingEquiv_natural
    P (chosenAdmissibleTranslation P hdeg) hP
    (chosenAdmissibleTranslation_linear_ne_zero P hdeg)
    (chosenAdmissibleTranslation_cubic_ne_zero P hdeg) f φ

#print axioms generalGaugeInversePolynomial_realization
#print axioms realizationJacobianOneFiberRepresentingEquiv_natural
#print axioms automaticJacobianOneFiberRepresentingEquiv_natural

end Automatic

end FiniteEtaleKeller
