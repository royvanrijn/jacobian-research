/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.Tactic.Ring

/-!
# The polynomial cusp identity

This is the algebraic identity used in Section 8 of the GVC manuscript.  It
is proved over an arbitrary commutative ring, so in particular it does not
depend on analytic sphere coordinates or characteristic zero.
-/

namespace GVC

/-- The defining identity for the homogeneous cusp lift. -/
theorem cusp_identity
    {R : Type*} [CommRing R] (x y t : R) :
    let ρ := t ^ 2 + x * y
    let A := ρ + x ^ 2
    let C := y * ρ ^ 2 - 2 * x * t ^ 2 * ρ - x ^ 3 * t ^ 2
    x * C = ρ ^ 3 - t ^ 2 * A ^ 2 := by
  dsimp
  ring

end GVC
