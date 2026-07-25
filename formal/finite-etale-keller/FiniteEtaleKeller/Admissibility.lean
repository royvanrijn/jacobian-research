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

end FiniteEtaleKeller
