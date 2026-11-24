/-
Adapted from GMC2/SimpleRootPartialFractions.lean at MurrellGroup/GMC-2-lean commit
1782de7ff6c97eb1d98e63e7ff34df18b9cd322e. MIT license in this directory.
-/
import Mathlib.FieldTheory.RatFunc.AsPolynomial
import Mathlib.LinearAlgebra.Lagrange

namespace MathieuQuadratic.DvKProof

noncomputable section

open Polynomial

/--
A nonzero polynomial has no repeated roots if its derivative is nonzero at
every root.
-/
theorem nodup_roots_of_derivative_ne_zero_at_roots
    {F : Type*} [Field F] [DecidableEq F]
    (U : Polynomial F) (hU : U ≠ 0)
    (hderiv :
      ∀ a : F, U.IsRoot a → U.derivative.eval a ≠ 0) :
    U.roots.Nodup := by
  rw [Multiset.nodup_iff_count_le_one]
  intro a
  rw [Polynomial.count_roots]
  by_contra hle
  have hmultiple : 1 < U.rootMultiplicity a :=
    Nat.lt_of_not_ge hle
  have hroots :=
    (Polynomial.one_lt_rootMultiplicity_iff_isRoot hU).mp
      hmultiple
  exact (hderiv a hroots.1) hroots.2

/--
Lagrange interpolation over all roots of a nonzero split polynomial with
simple roots, written in the form needed for partial fractions.
-/
theorem polynomial_eq_sum_over_simple_roots
    {F : Type*} [Field F] [DecidableEq F]
    (U : Polynomial F) (hU : U ≠ 0) (hsplit : U.Splits)
    (hsimple : U.roots.Nodup) (n : ℕ) (hdegree : n < U.natDegree) :
    Polynomial.X ^ n =
      ∑ a ∈ U.roots.toFinset,
        Polynomial.C (a ^ n / U.derivative.eval a) *
          (Polynomial.C U.leadingCoeff *
            ∏ b ∈ U.roots.toFinset.erase a,
              (Polynomial.X - Polynomial.C b)) := by
  classical
  let r := U.roots.toFinset
  have hinj : Set.InjOn (id : F → F) r :=
    Function.injective_id.injOn
  have hcard : r.card = U.natDegree := by
    rw [Multiset.toFinset_card_of_nodup hsimple,
      hsplit.natDegree_eq_card_roots]
  have hdegree' :
      (Polynomial.X ^ n : Polynomial F).degree < r.card := by
    rw [Polynomial.degree_X_pow, hcard]
    exact_mod_cast hdegree
  rw [Lagrange.eq_interpolate hinj hdegree',
    Lagrange.interpolate_eq_sum]
  apply Finset.sum_congr rfl
  intro a ha
  have haroots : a ∈ U.roots := by
    simpa [r] using ha
  have herase :
      @Finset.mk F (U.roots.erase a) (hsimple.erase a) =
        r.erase a := by
    ext b
    simp only [Finset.mem_mk, hsimple.mem_erase_iff,
      Finset.mem_erase, r, Multiset.mem_toFinset]
  have hderiv :
      U.derivative.eval a =
        U.leadingCoeff * ∏ b ∈ r.erase a, (a - b) := by
    have hfactor := hsplit.eq_prod_roots
    calc
      U.derivative.eval a =
          (Polynomial.C U.leadingCoeff *
            (U.roots.map fun x =>
              Polynomial.X - Polynomial.C x).prod).derivative.eval a :=
        congrArg (fun Q : Polynomial F => Q.derivative.eval a) hfactor
      _ = U.leadingCoeff *
          ((U.roots.erase a).map fun b => a - b).prod := by
        rw [Polynomial.derivative_mul, Polynomial.derivative_C,
          zero_mul, zero_add, Polynomial.eval_mul,
          Polynomial.eval_C,
          eval_multiset_prod_X_sub_C_derivative haroots]
      _ = U.leadingCoeff * ∏ b ∈ r.erase a, (a - b) := by
        congr 1
        rw [← herase]
        exact (Finset.prod_mk _ _ _).symm
  have hleading : U.leadingCoeff ≠ 0 :=
    Polynomial.leadingCoeff_ne_zero.mpr hU
  have hprod :
      (∏ b ∈ r.erase a, (a - b)) ≠ 0 := by
    apply Finset.prod_ne_zero_iff.mpr
    intro b hb
    exact sub_ne_zero.mpr fun hab =>
      (Finset.mem_erase.mp hb).1 hab.symm
  simp only [Polynomial.eval_pow, Polynomial.eval_X]
  rw [hderiv]
  have hscalar :
      a ^ n / (∏ b ∈ r.erase a, (a - b)) =
        (a ^ n /
          (U.leadingCoeff *
            ∏ b ∈ r.erase a, (a - b))) *
          U.leadingCoeff := by
    field_simp [hleading, hprod]
  change
    Polynomial.C
        (a ^ n / (∏ b ∈ r.erase a, (a - b))) *
        (∏ b ∈ r.erase a,
          (Polynomial.X - Polynomial.C b)) =
      Polynomial.C
          (a ^ n /
            (U.leadingCoeff *
              ∏ b ∈ r.erase a, (a - b))) *
        (Polynomial.C U.leadingCoeff *
          ∏ b ∈ r.erase a,
            (Polynomial.X - Polynomial.C b))
  rw [← mul_assoc, ← Polynomial.C_mul, ← hscalar]

/--
The simple-root interpolation identity with each complementary factor
written as the monic quotient `U / (X-a)`.
-/
theorem polynomial_eq_sum_divByMonic_over_simple_roots
    {F : Type*} [Field F] [DecidableEq F]
    (U : Polynomial F) (hU : U ≠ 0) (hsplit : U.Splits)
    (hsimple : U.roots.Nodup) (n : ℕ) (hdegree : n < U.natDegree) :
    Polynomial.X ^ n =
      ∑ a ∈ U.roots.toFinset,
        Polynomial.C (a ^ n / U.derivative.eval a) *
          (U /ₘ (Polynomial.X - Polynomial.C a)) := by
  classical
  let r := U.roots.toFinset
  let Q : F → Polynomial F := fun a ↦
    Polynomial.C U.leadingCoeff *
      ∏ b ∈ r.erase a, (Polynomial.X - Polynomial.C b)
  have hmain :
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
            @Finset.mk F U.roots hsimple = r :=
          Multiset.toFinset_eq hsimple
        rw [← hr]
        exact (Finset.prod_mk _ _ _).symm
  rw [hmain]
  apply Finset.sum_congr rfl
  intro a ha
  congr 1
  have haRoot : U.IsRoot a := by
    apply Polynomial.isRoot_of_mem_roots
    exact Multiset.mem_toFinset.mp ha
  have hdiv :
      (Polynomial.X - Polynomial.C a) *
          (U /ₘ (Polynomial.X - Polynomial.C a)) =
        U :=
    Polynomial.mul_divByMonic_eq_iff_isRoot.mpr haRoot
  have hfactor :
      (Polynomial.X - Polynomial.C a) * Q a = U := by
    calc
      (Polynomial.X - Polynomial.C a) * Q a =
          Polynomial.C U.leadingCoeff *
            ((Polynomial.X - Polynomial.C a) *
              ∏ b ∈ r.erase a,
                (Polynomial.X - Polynomial.C b)) := by
        simp only [Q]
        ring
      _ = Polynomial.C U.leadingCoeff *
          ∏ b ∈ r, (Polynomial.X - Polynomial.C b) := by
        rw [Finset.mul_prod_erase r
          (fun b ↦ Polynomial.X - Polynomial.C b) ha]
      _ = U := hUprod.symm
  exact mul_left_cancel₀
    (Polynomial.X_sub_C_ne_zero a)
    (hfactor.trans hdiv.symm)

/--
The explicit simple-root partial-fraction expansion of `X^n / U`.
-/
theorem ratfunc_eq_sum_over_simple_roots
    {F : Type*} [Field F] [DecidableEq F]
    (U : Polynomial F) (hU : U ≠ 0) (hsplit : U.Splits)
    (hsimple : U.roots.Nodup) (n : ℕ) (hdegree : n < U.natDegree) :
    algebraMap (Polynomial F) (RatFunc F) (Polynomial.X ^ n) /
        algebraMap (Polynomial F) (RatFunc F) U =
      ∑ a ∈ U.roots.toFinset,
        RatFunc.C (a ^ n / U.derivative.eval a) /
          (RatFunc.X - RatFunc.C a) := by
  classical
  let r := U.roots.toFinset
  have hpoly :=
    polynomial_eq_sum_over_simple_roots
      U hU hsplit hsimple n hdegree
  have hU0 : algebraMap (Polynomial F) (RatFunc F) U ≠ 0 :=
    RatFunc.algebraMap_ne_zero hU
  apply (div_eq_iff hU0).2
  rw [Finset.sum_mul]
  have hUprod :
      U = Polynomial.C U.leadingCoeff *
        ∏ a ∈ r, (Polynomial.X - Polynomial.C a) := by
    calc
      U = Polynomial.C U.leadingCoeff *
          (U.roots.map fun a =>
            Polynomial.X - Polynomial.C a).prod :=
        hsplit.eq_prod_roots
      _ = Polynomial.C U.leadingCoeff *
          ∏ a ∈ r, (Polynomial.X - Polynomial.C a) := by
        congr 1
        have hr :
            @Finset.mk F U.roots hsimple = r :=
          Multiset.toFinset_eq hsimple
        rw [← hr]
        exact (Finset.prod_mk _ _ _).symm
  rw [hpoly]
  simp only [map_sum]
  apply Finset.sum_congr rfl
  intro a ha
  have hpolyfactor :
      U = (Polynomial.X - Polynomial.C a) *
        (Polynomial.C U.leadingCoeff *
          ∏ b ∈ r.erase a,
            (Polynomial.X - Polynomial.C b)) := by
    calc
      U = Polynomial.C U.leadingCoeff *
          ∏ b ∈ r, (Polynomial.X - Polynomial.C b) :=
        hUprod
      _ = Polynomial.C U.leadingCoeff *
          ((Polynomial.X - Polynomial.C a) *
            ∏ b ∈ r.erase a,
              (Polynomial.X - Polynomial.C b)) := by
        rw [Finset.mul_prod_erase r
          (fun b => Polynomial.X - Polynomial.C b) ha]
      _ = (Polynomial.X - Polynomial.C a) *
          (Polynomial.C U.leadingCoeff *
            ∏ b ∈ r.erase a,
              (Polynomial.X - Polynomial.C b)) := by
        ring
  have hUfactor :
      algebraMap (Polynomial F) (RatFunc F) U =
        (RatFunc.X - RatFunc.C a) *
          algebraMap (Polynomial F) (RatFunc F)
            (Polynomial.C U.leadingCoeff *
              ∏ b ∈ r.erase a,
                (Polynomial.X - Polynomial.C b)) := by
    calc
      algebraMap (Polynomial F) (RatFunc F) U =
          algebraMap (Polynomial F) (RatFunc F)
            ((Polynomial.X - Polynomial.C a) *
              (Polynomial.C U.leadingCoeff *
                ∏ b ∈ r.erase a,
                  (Polynomial.X - Polynomial.C b))) :=
        congrArg (algebraMap (Polynomial F) (RatFunc F))
          hpolyfactor
      _ = (RatFunc.X - RatFunc.C a) *
          algebraMap (Polynomial F) (RatFunc F)
            (Polynomial.C U.leadingCoeff *
              ∏ b ∈ r.erase a,
                (Polynomial.X - Polynomial.C b)) := by
        simp
  have hdenom : RatFunc.X - RatFunc.C a ≠ 0 := by
    simpa using
      RatFunc.algebraMap_ne_zero
        (Polynomial.X_sub_C_ne_zero a)
  rw [hUfactor]
  simp only [map_mul, map_prod, map_sub, RatFunc.algebraMap_C,
    RatFunc.algebraMap_X]
  let q : F := a ^ n / U.derivative.eval a
  let d : RatFunc F := RatFunc.X - RatFunc.C a
  let P : RatFunc F :=
    RatFunc.C U.leadingCoeff *
      ∏ b ∈ r.erase a, (RatFunc.X - RatFunc.C b)
  change RatFunc.C q * P = RatFunc.C q / d * (d * P)
  field_simp [d, hdenom]

end

end MathieuQuadratic.DvKProof
