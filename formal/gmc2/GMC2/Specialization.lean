import Mathlib.RingTheory.FiniteType
import Mathlib.RingTheory.Ideal.Quotient.Operations

/-!
# Good finite-characteristic specializations

The lower-face argument only needs the following reusable interface: a
finitely generated characteristic-zero domain, together with a chosen
nonzero element, has good reductions in arbitrarily large prime
characteristic.  The scheme-theoretic proof is deliberately isolated here.
-/

namespace GMC2

/-- A reduction of `A` to a field of prime characteristic in which `c`
survives. -/
structure GoodReduction (A : Type*) [CommRing A] (c : A) where
  p : ℕ
  hp : p.Prime
  κ : Type
  fieldκ : Field κ
  charκ : CharP κ p
  reduce : A →+* κ
  c_ne_zero : reduce c ≠ 0

attribute [instance] GoodReduction.fieldκ GoodReduction.charκ

/-- Standard specialization theorem for a finite-type `ℤ`-domain.

This is kept as an imported algebraic-geometric interface until the relevant
generic-freeness/Chevalley package is available in mathlib.
-/
axiom exists_goodReduction
    (A : Type*) [CommRing A] [IsDomain A] [CharZero A] [Algebra ℤ A]
    [Algebra.FiniteType ℤ A] (c : A) (hc : c ≠ 0) (N : ℕ) :
    ∃ q : GoodReduction A c, N < q.p

end GMC2
