import GMC2.LowerFaceExtraction
import GMC2.NormalizedMoment
import GMC2.Specialization
import Mathlib.Tactic

/-!
# The finite-type specialization contradiction

This is the arithmetic half of the lower-face argument.  It starts from the
exposed characteristic-zero coefficient, cancels the lowest factorial,
specializes to characteristic `p`, and invokes Frobenius plus factorial
isolation.
-/

namespace GMC2

open Polynomial LaurentPolynomial

/-- Coefficientwise mapping of circular polynomials. -/
noncomputable def circularMap
    {R S : Type*} [CommRing R] [CommRing S] (φ : R →+* S) :
    CircularPolynomial R →+* CircularPolynomial S :=
  AddMonoidAlgebra.mapRingHom ℤ (Polynomial.mapRingHom φ)

@[simp] theorem circularMap_coeff
    {R S : Type*} [CommRing R] [CommRing S]
    (φ : R →+* S) (P : CircularPolynomial R) (k : ℤ) :
    (circularMap φ P).coeff k = (P.coeff k).map φ := by
  exact AddMonoidAlgebra.coeff_mapRingHom _ _ _

@[simp] theorem constantTerm_circularMap
    {R S : Type*} [CommRing R] [CommRing S]
    (φ : R →+* S) (P : CircularPolynomial R) :
    constantTerm (circularMap φ P) = (constantTerm P).map φ := by
  simp [constantTerm_eq_coeff]

theorem factorialFunctional_map
    {R S : Type*} [CommRing R] [CommRing S]
    (φ : R →+* S) (f : R[X]) :
    φ (factorialFunctional f) = factorialFunctional (f.map φ) := by
  induction f using Polynomial.induction_on' with
  | add f g hf hg =>
      simp [hf, hg]
  | monomial j a =>
      simp [factorialFunctional]

theorem circularExpectation_circularMap
    {R S : Type*} [CommRing R] [CommRing S]
    (φ : R →+* S) (P : CircularPolynomial R) :
    φ (circularExpectation P) =
      circularExpectation (circularMap φ P) := by
  rw [circularExpectation, circularExpectation, constantTerm_circularMap,
    ← factorialFunctional_map]

/-- Over a finite-type characteristic-zero domain, vanishing pure moments
exclude a straddling angular support. -/
theorem not_supportStraddlesZero_of_pureMoments_finiteType
    {A : Type*} [CommRing A] [IsDomain A] [CharZero A]
    [Algebra ℤ A] [Algebra.FiniteType ℤ A]
    (P : CircularPolynomial A) (hvan : PureMomentsVanish P) :
    ¬ SupportStraddlesZero P := by
  intro hconv
  let E := exposedLowestTerm_of_straddles P hconv
  let f : A[X] := constantTerm (P ^ E.r)
  let c : A := f.coeff E.n
  have hc : c ≠ 0 := by
    exact E.coeff_ne_zero
  obtain ⟨q, -⟩ := exists_goodReduction A c hc 0
  let φ : A →+* q.κ := q.reduce
  let fbar : q.κ[X] := f.map φ
  have hlowbar : ∀ j < E.n, fbar.coeff j = 0 := by
    intro j hj
    dsimp [fbar, f]
    rw [Polynomial.coeff_map]
    rw [E.coeff_below_eq_zero j hj, map_zero]
  have hcbar : fbar.coeff E.n ≠ 0 := by
    simpa [fbar, c, Polynomial.coeff_map] using q.c_ne_zero
  let H : A[X] := constantTerm (P ^ (E.r * q.p))
  have hHlow : ∀ j < E.n * q.p, H.coeff j = 0 := by
    intro j hj
    exact E.scaled_coeff_below_eq_zero q.p j hj
  have hindex : 0 < E.r * q.p := Nat.mul_pos E.r_pos q.hp.pos
  have hmoment : factorialFunctional H = 0 := by
    have hm := hvan (E.r * q.p) hindex
    simpa [pureMoment, circularExpectation, H] using hm
  have hnormalized :
      normalizedFactorialMoment (E.n * q.p) H = 0 :=
    normalizedFactorialMoment_eq_zero_of_factorialFunctional_eq_zero
      (E.n * q.p) H hHlow hmoment
  have hmapped :
      normalizedFactorialMoment (E.n * q.p) (H.map φ) = 0 := by
    rw [← normalizedFactorialMoment_map φ]
    simp [hnormalized]
  have hFrobenius : H.map φ = fbar ^ q.p := by
    calc
      H.map φ =
          constantTerm (circularMap φ (P ^ (E.r * q.p))) := by
            exact (constantTerm_circularMap φ
              (P ^ (E.r * q.p))).symm
      _ = constantTerm ((circularMap φ (P ^ E.r)) ^ q.p) := by
            congr 1
            rw [← map_pow, ← pow_mul]
      _ = constantTerm (circularMap φ (P ^ E.r)) ^ q.p :=
            constantTerm_frobenius q.hp _
      _ = fbar ^ q.p := by
            rw [constantTerm_circularMap]
  have hprimeMoment :
      primeDilatedMoment q.p E.n fbar = 0 := by
    rw [hFrobenius, normalizedFactorialMoment_frobenius q.hp fbar] at hmapped
    exact hmapped
  exact hcbar
    (lowestCoeff_eq_zero_of_primeDilatedMoment_eq_zero
      q.hp fbar hlowbar hprimeMoment)

end GMC2
