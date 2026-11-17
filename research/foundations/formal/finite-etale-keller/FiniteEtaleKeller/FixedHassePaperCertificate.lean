/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.FixedHasseArithmetic

/-!
# Certificate accompanying the fixed-map Hasse-failure paper

This is the publication-sized Lean entry point for the algebraic,
local-global, height, and elementary arithmetic layers of
*Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map*.

The Selberg--Delange asymptotic and the Berend--Bilu degree obstruction are
ordinary mathematical inputs and deliberately remain outside this module.
-/

noncomputable section

namespace FiniteEtaleKeller.FixedHasseFamily

/-- The paper-facing endpoint for every admissible parameter.  Its result
contains the exact displayed Jacobian-one map, its geometric degree-five
certificate and inverse affine normalization, the literal rank-five
finite-étale fiber, the Hasse failure over all completions, and the primitive
target of height `32a`. -/
theorem fixedHassePaper_certificate
    (a : ℕ) (ha : 1 < a) (hmod9 : a % 9 = 1)
    (hsupport :
      ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1)
    (hnoncube : ¬∃ r : ℚ, r ^ 3 = (a : ℚ)) :
    PaperParameterCertificate a :=
  paperParameter_certificate a ha hmod9 hsupport hnoncube

/-- Compile-time guard for the prime progression stated in the paper. -/
example (ℓ : ℕ) (hprime : ℓ.Prime) (hmod9 : ℓ % 9 = 1) :
    PaperParameterCertificate ℓ :=
  primeParameter_certificate ℓ hprime hmod9

#print axioms paperMap_normalization_inverse
#print axioms paperMap_geometricDegree
#print axioms jacobianDet_paperMap
#print axioms paperFiberRepresentingEquiv
#print axioms paperFiberPoint_hasse_certificate
#print axioms targetProjectiveContent_eq_one
#print axioms targetProjectiveHeight_eq
#print axioms paperParameter_certificate
#print axioms fixedHassePaper_certificate
#print axioms HasseCoreCondition.mul
#print axioms prime_not_rational_cube
#print axioms primeParameter_certificate

end FiniteEtaleKeller.FixedHasseFamily
