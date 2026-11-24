/-
Adapted from GMC2/AnnulusPartialFractions.lean at MurrellGroup/GMC-2-lean commit
1782de7ff6c97eb1d98e63e7ff34df18b9cd322e. MIT license in this directory.
-/
import MathieuQuadratic.External.DvK.AnnulusSeries
import MathieuQuadratic.External.DvK.SimpleRootPartialFractions

open scoped WithZero

namespace MathieuQuadratic.DvKProof

local instance coeLaurentAnnulusPartialFractions {R : Type*} [Semiring R] : CoeFun (LaurentPolynomial R) (fun _ ↦ ℤ → R) :=
  ⟨fun p ↦ p.coeff⟩


noncomputable section

open Polynomial

namespace AnnulusSeries

variable {K : Type*} [Field K] [DecidableEq K] [Valued K ℤᵐ⁰]
  [NonarchimedeanRing K] [CompleteSpace K]

omit [DecidableEq K] in
/--
An indexed partial-fraction expansion solves the corresponding polynomial
equation in annulus series.
-/
theorem laurentAction_indexedLinearInverseSum_of_identity
    {ι : Type*} (indices : Finset ι) (root c : ι → K)
    (hunit : ∀ i ∈ indices, Valued.v (root i) ≠ 1)
    (U N : Polynomial K) (Q : ι → Polynomial K)
    (hfactor :
      ∀ i ∈ indices,
        U = Q i * (Polynomial.X - Polynomial.C (root i)))
    (hidentity :
      N = ∑ i ∈ indices, Polynomial.C (c i) * Q i) :
    laurentAction (Polynomial.toLaurent U)
        (indexedLinearInverseSum indices root c hunit) =
      laurentAction (Polynomial.toLaurent N) (single 0 1) := by
  classical
  change
    laurentAction (Polynomial.toLaurent U)
        (∑ i ∈ indices.attach,
          scale (c i.1)
            (linearInverse (root i.1) (hunit i.1 i.2))) =
      laurentAction (Polynomial.toLaurent N) (single 0 1)
  change
    (laurentActionHom (Polynomial.toLaurent U))
        (∑ i ∈ indices.attach,
          scale (c i.1)
            (linearInverse (root i.1) (hunit i.1 i.2))) =
      laurentAction (Polynomial.toLaurent N) (single 0 1)
  rw [map_sum]
  change
    ∑ i ∈ indices.attach,
        laurentAction (Polynomial.toLaurent U)
          (scale (c i.1)
            (linearInverse (root i.1) (hunit i.1 i.2))) =
      laurentAction (Polynomial.toLaurent N) (single 0 1)
  calc
    _ = ∑ i ∈ indices,
        laurentAction
          (Polynomial.toLaurent (Polynomial.C (c i) * Q i))
          (single 0 1) := by
      rw [← indices.sum_attach]
      apply Finset.sum_congr rfl
      intro i hi
      have hfactorLaurent :
          Polynomial.toLaurent U =
            Polynomial.toLaurent (Q i.1) *
              (LaurentPolynomial.T 1 -
                LaurentPolynomial.C (root i.1)) := by
        rw [hfactor i.1 i.2, map_mul, map_sub,
          Polynomial.toLaurent_X, Polynomial.toLaurent_C]
      have hscalarLaurent :
          Polynomial.toLaurent
              (Polynomial.C (c i.1) * Q i.1) =
            LaurentPolynomial.C (c i.1) *
              Polynomial.toLaurent (Q i.1) := by
        rw [map_mul, Polynomial.toLaurent_C]
      change
        laurentAction (Polynomial.toLaurent U)
            (c i.1 •
              linearInverse (root i.1) (hunit i.1 i.2)) =
          _
      rw [laurentAction_smul, hfactorLaurent,
        laurentAction_mul,
        laurentAction_T_sub_C_linearInverse]
      rw [hscalarLaurent, laurentAction_mul,
        laurentAction_C]
    _ = laurentAction
          (Polynomial.toLaurent
            (∑ i ∈ indices, Polynomial.C (c i) * Q i))
          (single 0 1) := by
      simp only [map_sum]
      symm
      change
        (laurentActionHom
          (∑ i ∈ indices,
            Polynomial.toLaurent
              (Polynomial.C (c i) * Q i)))
            (single 0 1) =
          ∑ i ∈ indices,
            (laurentActionHom
              (Polynomial.toLaurent
                (Polynomial.C (c i) * Q i)))
              (single 0 1)
      rw [map_sum]
      rw [LinearMap.sum_apply]
    _ = laurentAction (Polynomial.toLaurent N) (single 0 1) := by
      rw [hidentity]

/--
The annulus expansion of the simple-root partial fractions of `X^n / U`
is a solution of `U * u = X^n`.
-/
theorem laurentAction_linearInverseSum
    (U : Polynomial K) (hU : U ≠ 0) (hsplit : U.Splits)
    (hsimple : U.roots.Nodup) (n : ℕ)
    (hdegree : n < U.natDegree)
    (hunit : ∀ a ∈ U.roots.toFinset, Valued.v a ≠ 1) :
    laurentAction (Polynomial.toLaurent U)
        (linearInverseSum U.roots.toFinset
          (fun a ↦ a ^ n / U.derivative.eval a) hunit) =
      single n 1 := by
  classical
  let r := U.roots.toFinset
  let Q : K → Polynomial K := fun a ↦
    Polynomial.C U.leadingCoeff *
      ∏ b ∈ r.erase a, (Polynomial.X - Polynomial.C b)
  have hpoly :
      Polynomial.X ^ n =
        ∑ a ∈ r,
          Polynomial.C (a ^ n / U.derivative.eval a) * Q a := by
    simpa [r, Q] using
      polynomial_eq_sum_over_simple_roots
        U hU hsplit hsimple n hdegree
  have hUprod :
      U = Polynomial.C U.leadingCoeff *
        ∏ a ∈ r, (Polynomial.X - Polynomial.C a) := by
    calc
      U = Polynomial.C U.leadingCoeff *
          (U.roots.map fun a ↦
            Polynomial.X - Polynomial.C a).prod :=
        hsplit.eq_prod_roots
      _ = Polynomial.C U.leadingCoeff *
          ∏ a ∈ r, (Polynomial.X - Polynomial.C a) := by
        congr 1
        have hr :
            @Finset.mk K U.roots hsimple = r :=
          Multiset.toFinset_eq hsimple
        rw [← hr]
        exact (Finset.prod_mk _ _ _).symm
  change
    laurentAction (Polynomial.toLaurent U)
        (∑ a ∈ r.attach,
          scale (a.1 ^ n / U.derivative.eval a.1)
            (linearInverse a.1 (hunit a.1 a.2))) =
      single n 1
  change
    (laurentActionHom (Polynomial.toLaurent U))
        (∑ a ∈ r.attach,
          scale (a.1 ^ n / U.derivative.eval a.1)
            (linearInverse a.1 (hunit a.1 a.2))) =
      single n 1
  rw [map_sum]
  change
    ∑ a ∈ r.attach,
        laurentAction (Polynomial.toLaurent U)
          (scale (a.1 ^ n / U.derivative.eval a.1)
            (linearInverse a.1 (hunit a.1 a.2))) =
      single n 1
  calc
    _ = ∑ a ∈ r,
        laurentAction
          (Polynomial.toLaurent
            (Polynomial.C (a ^ n / U.derivative.eval a) * Q a))
          (single 0 1) := by
      rw [← r.sum_attach]
      apply Finset.sum_congr rfl
      intro a ha
      have hfactor :
          U = Q a.1 *
            (Polynomial.X - Polynomial.C a.1) := by
        calc
          U = Polynomial.C U.leadingCoeff *
              ∏ b ∈ r, (Polynomial.X - Polynomial.C b) :=
            hUprod
          _ = Polynomial.C U.leadingCoeff *
              ((Polynomial.X - Polynomial.C a.1) *
                ∏ b ∈ r.erase a.1,
                  (Polynomial.X - Polynomial.C b)) := by
            rw [Finset.mul_prod_erase r
              (fun b ↦ Polynomial.X - Polynomial.C b) a.2]
          _ = Q a.1 *
              (Polynomial.X - Polynomial.C a.1) := by
            simp only [Q]
            ring
      have hfactorLaurent :
          Polynomial.toLaurent U =
            Polynomial.toLaurent (Q a.1) *
              (LaurentPolynomial.T 1 -
                LaurentPolynomial.C a.1) := by
        rw [hfactor, map_mul, map_sub,
          Polynomial.toLaurent_X, Polynomial.toLaurent_C]
      have hscalarLaurent :
          Polynomial.toLaurent
              (Polynomial.C
                (a.1 ^ n / U.derivative.eval a.1) * Q a.1) =
            LaurentPolynomial.C
                (a.1 ^ n / U.derivative.eval a.1) *
              Polynomial.toLaurent (Q a.1) := by
        rw [map_mul, Polynomial.toLaurent_C]
      change
        laurentAction (Polynomial.toLaurent U)
            ((a.1 ^ n / U.derivative.eval a.1) •
              linearInverse a.1 (hunit a.1 a.2)) =
          _
      rw [laurentAction_smul, hfactorLaurent,
        laurentAction_mul,
        laurentAction_T_sub_C_linearInverse]
      rw [hscalarLaurent, laurentAction_mul,
        laurentAction_C]
    _ = laurentAction
          (Polynomial.toLaurent
            (∑ a ∈ r,
              Polynomial.C (a ^ n / U.derivative.eval a) * Q a))
          (single 0 1) := by
      simp only [map_sum]
      symm
      change
        (laurentActionHom
          (∑ a ∈ r,
            Polynomial.toLaurent
              (Polynomial.C (a ^ n / U.derivative.eval a) * Q a)))
            (single 0 1) =
          ∑ a ∈ r,
            (laurentActionHom
              (Polynomial.toLaurent
                (Polynomial.C (a ^ n / U.derivative.eval a) * Q a)))
              (single 0 1)
      rw [map_sum]
      rw [LinearMap.sum_apply]
    _ = laurentAction
          (Polynomial.toLaurent (Polynomial.X ^ n))
          (single 0 1) := by rw [hpoly]
    _ = single n 1 := by
      rw [Polynomial.toLaurent_X_pow,
        laurentAction_T_single_zero]

end AnnulusSeries

end

end MathieuQuadratic.DvKProof
