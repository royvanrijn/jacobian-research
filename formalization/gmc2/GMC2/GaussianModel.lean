import GMC2.ConstantTerm
import Mathlib.Algebra.Polynomial.Degree.TrailingDegree

/-!
# The circular Gaussian model

After the invertible circular change of coordinates, a polynomial in two
Gaussian variables is a Laurent polynomial in the angular variable `T`
whose coefficients are polynomials in the radial variable `U`.
-/

namespace GMC2

open Polynomial LaurentPolynomial

/-- Polynomials in the circular model: Laurent in `T`, polynomial in `U`. -/
abbrev CircularPolynomial (R : Type*) [CommRing R] :=
  LaurentPolynomial R[X]

/-- The radial factorial functional `U^j ↦ j!`. -/
noncomputable def factorialFunctional
    {R : Type*} [CommRing R] (f : R[X]) : R :=
  f.sum fun j a ↦ a * (Nat.factorial j : R)

@[simp] theorem factorialFunctional_monomial
    {R : Type*} [CommRing R] (j : ℕ) (a : R) :
    factorialFunctional (Polynomial.monomial j a) =
      a * (Nat.factorial j : R) := by
  simp [factorialFunctional]

/-- Gaussian expectation in circular coordinates. -/
noncomputable def circularExpectation
    {R : Type*} [CommRing R] (P : CircularPolynomial R) : R :=
  factorialFunctional (constantTerm P)

/-- The `m`-th pure Gaussian moment. -/
noncomputable def pureMoment
    {R : Type*} [CommRing R] (P : CircularPolynomial R) (m : ℕ) : R :=
  circularExpectation (P ^ m)

/-- The mixed moment with a fixed multiplier `Q`. -/
noncomputable def mixedMoment
    {R : Type*} [CommRing R] (Q P : CircularPolynomial R) (m : ℕ) : R :=
  circularExpectation (Q * P ^ m)

def PureMomentsVanish {R : Type*} [CommRing R]
    (P : CircularPolynomial R) : Prop :=
  ∀ m : ℕ, 0 < m → pureMoment P m = 0

def EventuallyMixedMomentsVanish {R : Type*} [CommRing R]
    (P : CircularPolynomial R) : Prop :=
  ∀ Q : CircularPolynomial R,
    ∃ M : ℕ, ∀ m ≥ M, mixedMoment Q P m = 0

/-- The angular support of a circular polynomial. -/
def angularSupport {R : Type*} [CommRing R]
    (P : CircularPolynomial R) : Finset ℤ :=
  P.coeff.support

/-- Every nonzero angular coefficient has positive weight. -/
def PositiveAngularSupport {R : Type*} [CommRing R]
    (P : CircularPolynomial R) : Prop :=
  ∀ k ∈ angularSupport P, 0 < k

/-- Every nonzero angular coefficient has negative weight. -/
def NegativeAngularSupport {R : Type*} [CommRing R]
    (P : CircularPolynomial R) : Prop :=
  ∀ k ∈ angularSupport P, k < 0

def OneSidedAngularSupport {R : Type*} [CommRing R]
    (P : CircularPolynomial R) : Prop :=
  PositiveAngularSupport P ∨ NegativeAngularSupport P

end GMC2
