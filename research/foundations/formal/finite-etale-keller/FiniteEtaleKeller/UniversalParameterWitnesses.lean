/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.UniversalParameterCompiler

/-!
# Exact low-rank witnesses for the universal promoted map

These rational witness cards exercise three different rank-four finite étale
algebra types:

* the connected field presentation `T⁴-2`;
* the split presentation `(T-1)(T-2)(T-3)(T-4)`;
* the disconnected product `(T²-2)(T²-3)`.

The theorems check the selected inverse polynomials exactly.  Squarefreeness,
irreducibility, and the represented-fiber equivalence are supplied by the
general modules; these cards deliberately isolate the new promoted target
arithmetic.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller
namespace UniversalParameterWitnesses

/-- Target `(u₄,π,b,c)=(1/4,1,3/2,1/2)` for `T⁴-2`, translated by one. -/
def connectedQuarticTarget : UniversalPromotedTarget ℚ where
  parameter j := if j = 4 then 1 / 4 else 0
  pi := 1
  b := 3 / 2
  c := 1 / 2

/-- Target `(u₄,π,b,c)=(-25/2,1/5,-7/10,24/25)` for the split quartic. -/
def splitQuarticTarget : UniversalPromotedTarget ℚ where
  parameter j := if j = 4 then -25 / 2 else 0
  pi := 1 / 5
  b := -7 / 10
  c := 24 / 25

/-- Target `(u₄,π,b,c)=(-27/32,-2/3,-1/6,2/3)` for
`(T²-2)(T²-3)`, translated by one. -/
def productQuarticTarget : UniversalPromotedTarget ℚ where
  parameter j := if j = 4 then -27 / 32 else 0
  pi := -2 / 3
  b := -1 / 6
  c := 2 / 3

theorem connectedQuartic_inverse :
    universalPromotedInversePolynomial 4 connectedQuarticTarget =
      C (-1 / 4) + X + C (3 / 2) * X ^ 2 + X ^ 3 +
        C (1 / 4) * X ^ 4 := by
  norm_num [universalPromotedInversePolynomial, universalPromotedTail,
    connectedQuarticTarget]
  ring

theorem splitQuartic_inverse :
    universalPromotedInversePolynomial 4 splitQuarticTarget =
      C (-12 / 25) + X + C (-7 / 10) * X ^ 2 +
        C (1 / 5) * X ^ 3 + C (-1 / 50) * X ^ 4 := by
  norm_num [universalPromotedInversePolynomial, universalPromotedTail,
    splitQuarticTarget]
  ring

theorem productQuartic_inverse :
    universalPromotedInversePolynomial 4 productQuarticTarget =
      C (-1 / 3) + X + C (-1 / 6) * X ^ 2 +
        C (-2 / 3) * X ^ 3 + C (-1 / 6) * X ^ 4 := by
  norm_num [universalPromotedInversePolynomial, universalPromotedTail,
    productQuarticTarget]
  ring

#print axioms connectedQuartic_inverse
#print axioms splitQuartic_inverse
#print axioms productQuartic_inverse

end UniversalParameterWitnesses
end FiniteEtaleKeller
