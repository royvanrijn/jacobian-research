/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Translation
import Mathlib.Algebra.Polynomial.Taylor

/-!
# Admissibility coefficients of the rooted translation

For `G(S)=P(a+S)-P(a)`, the linear coefficient is `P'(a)` and the cubic
coefficient is the third Hasse derivative of `P` at `a`.  In characteristic
zero the latter is `P'''(a)/3!`.

For a polynomial of degree at least three over a characteristic-zero field,
the product of these two derivative polynomials is nonzero.  Since the field is
infinite, there is therefore a translation parameter at which both evaluate
nontrivially.  The final definitions choose such a parameter once and for all.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- The rooted translation has zero constant coefficient. -/
@[simp]
theorem rootedTranslate_coeff_zero (P : K[X]) (a : K) :
    (rootedTranslate P a).coeff 0 = 0 := by
  change (Polynomial.taylor a P - C (P.eval a)).coeff 0 = 0
  simp

/-- The linear coefficient of `P(a+S)-P(a)` is `P'(a)`. -/
@[simp]
theorem rootedTranslate_coeff_one (P : K[X]) (a : K) :
    (rootedTranslate P a).coeff 1 = P.derivative.eval a := by
  change (Polynomial.taylor a P - C (P.eval a)).coeff 1 = P.derivative.eval a
  simp

/-- The cubic coefficient of `P(a+S)-P(a)` is the third Hasse derivative at
`a`, i.e. the characteristic-free Taylor coefficient. -/
@[simp]
theorem rootedTranslate_coeff_three (P : K[X]) (a : K) :
    (rootedTranslate P a).coeff 3 = (Polynomial.hasseDeriv 3 P).eval a := by
  change (Polynomial.taylor a P - C (P.eval a)).coeff 3 =
    (Polynomial.hasseDeriv 3 P).eval a
  simp [Polynomial.taylor_coeff]

/-- Nonvanishing of the derivative gives the required nonzero linear gauge
coefficient. -/
theorem rootedTranslate_linear_ne_zero (P : K[X]) (a : K)
    (h : P.derivative.eval a ≠ 0) :
    (rootedTranslate P a).coeff 1 ≠ 0 := by
  simpa using h

/-- Nonvanishing of the third Hasse derivative gives the required nonzero cubic
gauge coefficient. -/
theorem rootedTranslate_cubic_ne_zero (P : K[X]) (a : K)
    (h : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (rootedTranslate P a).coeff 3 ≠ 0 := by
  simpa using h

section AutomaticChoice

variable [CharZero K]

/-- A polynomial of degree at least three has nonzero derivative. -/
theorem derivative_ne_zero_of_three_le_natDegree (P : K[X])
    (hdeg : 3 ≤ P.natDegree) : P.derivative ≠ 0 := by
  apply Polynomial.derivative_ne_zero.mpr
  omega

/-- A polynomial of degree at least three has nonzero third Hasse derivative. -/
theorem hasseDeriv_three_ne_zero_of_three_le_natDegree (P : K[X])
    (hdeg : 3 ≤ P.natDegree) : Polynomial.hasseDeriv 3 P ≠ 0 := by
  have hP : P ≠ 0 := by
    intro hzero
    simpa [hzero] using hdeg
  by_cases htop : P.natDegree = 3
  · have hform : Polynomial.hasseDeriv 3 P = C P.leadingCoeff := by
      simpa [htop] using Polynomial.hasseDeriv_natDegree_eq_C P
    rw [hform]
    exact Polynomial.C_ne_zero.mpr (Polynomial.leadingCoeff_ne_zero.mpr hP)
  · intro hzero
    have hnat := Polynomial.natDegree_hasseDeriv P 3
    rw [hzero] at hnat
    simp only [Polynomial.natDegree_zero] at hnat
    have hlt : 3 < P.natDegree := lt_of_le_of_ne hdeg (Ne.symm htop)
    omega

/-- Two nonzero polynomials over an infinite field are simultaneously nonzero
at some field element. -/
theorem exists_common_nonzero_eval (P Q : K[X]) (hP : P ≠ 0) (hQ : Q ≠ 0) :
    ∃ a : K, P.eval a ≠ 0 ∧ Q.eval a ≠ 0 := by
  have hPQ : P * Q ≠ 0 := mul_ne_zero hP hQ
  by_contra h
  apply hPQ
  apply Polynomial.zero_of_eval_zero
  intro a
  by_cases hPa : P.eval a = 0
  · simp [hPa]
  by_cases hQa : Q.eval a = 0
  · simp [hQa]
  exact (h ⟨a, hPa, hQa⟩).elim

/-- Every polynomial of degree at least three admits a translation at which the
linear and cubic coefficients of `P(a+S)-P(a)` are both nonzero. -/
theorem exists_admissible_translation (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    ∃ a : K,
      P.derivative.eval a ≠ 0 ∧ (Polynomial.hasseDeriv 3 P).eval a ≠ 0 := by
  exact exists_common_nonzero_eval P.derivative (Polynomial.hasseDeriv 3 P)
    (derivative_ne_zero_of_three_le_natDegree P hdeg)
    (hasseDeriv_three_ne_zero_of_three_le_natDegree P hdeg)

/-- A fixed admissible translation parameter, chosen from the preceding
existence theorem. -/
def chosenAdmissibleTranslation (P : K[X]) (hdeg : 3 ≤ P.natDegree) : K :=
  Classical.choose (exists_admissible_translation P hdeg)

/-- The chosen translation has nonzero linear Taylor coefficient. -/
theorem chosenAdmissibleTranslation_linear_ne_zero (P : K[X])
    (hdeg : 3 ≤ P.natDegree) :
    P.derivative.eval (chosenAdmissibleTranslation P hdeg) ≠ 0 :=
  (Classical.choose_spec (exists_admissible_translation P hdeg)).1

/-- The chosen translation has nonzero cubic Taylor coefficient. -/
theorem chosenAdmissibleTranslation_cubic_ne_zero (P : K[X])
    (hdeg : 3 ≤ P.natDegree) :
    (Polynomial.hasseDeriv 3 P).eval (chosenAdmissibleTranslation P hdeg) ≠ 0 :=
  (Classical.choose_spec (exists_admissible_translation P hdeg)).2

#print axioms exists_admissible_translation

end AutomaticChoice

end FiniteEtaleKeller
