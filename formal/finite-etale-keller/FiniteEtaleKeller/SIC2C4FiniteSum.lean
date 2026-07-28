/-
Copyright (c) 2026 Jacobian Research contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Jacobian Research contributors
-/

import Mathlib.Algebra.Group.ForwardDiff
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Data.Nat.Factorial.DoubleFactorial
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# The residual SIC2C4 finite sum

This module formalizes the discrete part of the direct coefficient-extraction
proof of SIC2C4: finite-difference vanishing, the general remainder principle
and its rank-one endpoint-residue specialization for alternating quotient
sums, the specialized product polynomials, the first-order beta recurrence,
its factorial/double-factorial evaluation, the triangular recurrence for
arbitrary repeated poles, and the two normalized displayed binomial sums.
It also checks the monomial coefficient-extraction and scalar chart identities.
It deliberately does not claim to formalize their linear assembly for the full
witness.
-/

open scoped BigOperators
open Finset

namespace FiniteEtaleKeller.SIC2C4

/-- Algebraic chart identity underlying equation (4.3). -/
theorem chartIdentity {K : Type*} [Field K] [CharZero K]
    (x v : K) (hx : x ≠ 0) :
    (1 + x) * ((1 - v ^ 2) / (2 * x) - (2 + x) * v ^ 2 / 2) =
      (1 + x) / (2 * x) * (1 - v ^ 2 * (1 + x) ^ 2) := by
  field_simp [hx]
  ring

/-- Contraction of a balanced bidegree monomial after matching its two indices. -/
def monomialContraction (n a b : ℕ) : ℚ :=
  if a = b then (a.factorial : ℚ) * ((n - a).factorial : ℚ) else 0

/-- Constant-term and formal-beta value of the same monomial. -/
def monomialChartValue (n a b : ℕ) : ℚ :=
  if a = b then
    ((n + 1).factorial : ℚ) *
      ((a.factorial : ℚ) * ((n - a).factorial : ℚ) /
        ((n + 1).factorial : ℚ))
  else 0

/-- Monomial case of the contraction-to-chart formula (4.2). -/
theorem monomialCoefficientExtraction (n a b : ℕ) :
    monomialChartValue n a b = monomialContraction n a b := by
  unfold monomialChartValue monomialContraction
  split_ifs
  · have hfactorial : ((n + 1).factorial : ℚ) ≠ 0 := by positivity
    field_simp [hfactorial]
  · rfl

/-- The finite-difference cancellation used for the pure SIC2C4 sum. -/
theorem finiteDifference_vanishes (P : Polynomial ℚ) (m : ℕ)
    (hP : P.natDegree < m) :
    (∑ k ∈ range (m + 1),
      ((-1 : ℤ) ^ (m - k) * m.choose k) • P.eval (k : ℚ)) = 0 := by
  have hdiff := congr_fun (Polynomial.fwdDiff_iter_eq_zero_of_degree_lt hP) 0
  rw [Pi.zero_apply, fwdDiff_iter_eq_sum_shift] at hdiff
  simpa using hdiff

lemma negOnePow_sub_eq_mul (m k : ℕ) (hk : k ≤ m) :
    (-1 : ℚ) ^ (m - k) = (-1 : ℚ) ^ m * (-1 : ℚ) ^ k := by
  conv_rhs => lhs; rw [show m = (m - k) + k by omega, pow_add]
  rw [mul_assoc, ← mul_pow]
  norm_num

/-- The sign convention used in the displayed SIC2C4 sums. -/
theorem alternatingPolynomialSum_vanishes (P : Polynomial ℚ) (m : ℕ)
    (hP : P.natDegree < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) * P.eval (k : ℚ)) = 0 := by
  have hdiff := finiteDifference_vanishes P m hP
  have hq :
      (∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ (m - k) * (m.choose k : ℚ) * P.eval (k : ℚ)) = 0 := by
    simpa [zsmul_eq_mul] using hdiff
  calc
    ∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) * P.eval (k : ℚ)
      = (-1 : ℚ) ^ m *
          ∑ k ∈ range (m + 1),
            (-1 : ℚ) ^ (m - k) * (m.choose k : ℚ) * P.eval (k : ℚ) := by
              rw [mul_sum]
              apply sum_congr rfl
              intro k hk
              have hk_le : k ≤ m := by
                simp only [mem_range] at hk
                omega
              rw [negOnePow_sub_eq_mul m k hk_le]
              ring_nf
              rw [show m * 2 = 2 * m by omega, pow_mul]
              norm_num
    _ = 0 := by rw [hq, mul_zero]

/-- A summand in the residual beta sum. -/
def betaTerm (m k : ℕ) : ℚ :=
  (-1 : ℚ) ^ k * (m.choose k : ℚ) / (2 * k + 1)

/-- The residual sum left after polynomial finite-difference cancellation. -/
def betaSum (m : ℕ) : ℚ :=
  ∑ k ∈ range (m + 1), betaTerm m k

/-- Generalized beta summand for a pole of arbitrary positive order. -/
def repeatedBetaTerm (m s k : ℕ) : ℚ :=
  (-1 : ℚ) ^ k * (m.choose k : ℚ) / (2 * k + 1) ^ s

/-- Generalized beta sum governing a repeated pole at `X = -1/2`. -/
def repeatedBetaSum (m s : ℕ) : ℚ :=
  ∑ k ∈ range (m + 1), repeatedBetaTerm m s k

lemma repeatedBetaSum_zero_left (s : ℕ) : repeatedBetaSum 0 s = 1 := by
  simp [repeatedBetaSum, repeatedBetaTerm]

lemma repeatedBetaSum_zero_right (m : ℕ) (hm : 0 < m) :
    repeatedBetaSum m 0 = 0 := by
  have hdegree : (1 : Polynomial ℚ).natDegree < m := by
    simp
    exact hm
  simpa [repeatedBetaSum, repeatedBetaTerm] using
    alternatingPolynomialSum_vanishes (1 : Polynomial ℚ) m hdegree

lemma repeatedBetaSum_one (m : ℕ) :
    repeatedBetaSum m 1 = betaSum m := by
  simp [repeatedBetaSum, repeatedBetaTerm, betaSum, betaTerm]

/-- The mixed finite sum before removing its polynomial part. -/
def mixedTemplateSum (m : ℕ) (A : Polynomial ℚ) : ℚ :=
  ∑ k ∈ range (m + 1),
    (-1 : ℚ) ^ k * (m.choose k : ℚ) *
      (A.eval (k : ℚ) / (2 * k + 1))

/-- Alternating quotient transform with an arbitrary polynomial denominator. -/
def quotientTransform (m : ℕ) (L A : Polynomial ℚ) : ℚ :=
  ∑ k ∈ range (m + 1),
    (-1 : ℚ) ^ k * (m.choose k : ℚ) *
      (A.eval (k : ℚ) / L.eval (k : ℚ))

/-- The quotient transform ignores a denominator multiple whose quotient has
degree less than the finite-difference order. -/
theorem quotientTransform_add_mul_invariant
    (m : ℕ) (L A R D : Polynomial ℚ)
    (hfactor : A = R + L * D)
    (hD : D.natDegree < m)
    (hL : ∀ k ∈ range (m + 1), L.eval (k : ℚ) ≠ 0) :
    quotientTransform m L A = quotientTransform m L R := by
  have hvanish := alternatingPolynomialSum_vanishes D m hD
  rw [quotientTransform, quotientTransform]
  calc
    ∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (A.eval (k : ℚ) / L.eval (k : ℚ))
      = ∑ k ∈ range (m + 1),
          ((-1 : ℚ) ^ k * (m.choose k : ℚ) *
              (R.eval (k : ℚ) / L.eval (k : ℚ)) +
            (-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ)) := by
              apply sum_congr rfl
              intro k hk
              rw [hfactor]
              simp only [Polynomial.eval_add, Polynomial.eval_mul]
              have hden := hL k hk
              field_simp [hden]
              ring
    _ = (∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) *
            (R.eval (k : ℚ) / L.eval (k : ℚ))) +
        ∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ) := by
            rw [sum_add_distrib]
    _ = ∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) *
            (R.eval (k : ℚ) / L.eval (k : ℚ)) := by
              rw [hvanish, add_zero]

/-- The quotient transform factors through Mathlib's monic remainder
operation whenever the corresponding quotient has degree below the
finite-difference order.  For a nonmonic denominator this statement is
tautological because `divByMonic` is zero. -/
theorem quotientTransform_eq_modByMonic_of_quotient_degree
    (m : ℕ) (L A : Polynomial ℚ)
    (hquotient : (A /ₘ L).natDegree < m)
    (hL : ∀ k ∈ range (m + 1), L.eval (k : ℚ) ≠ 0) :
    quotientTransform m L A = quotientTransform m L (A %ₘ L) := by
  apply quotientTransform_add_mul_invariant m L A (A %ₘ L) (A /ₘ L)
  · exact (Polynomial.modByMonic_add_div A L).symm
  · exact hquotient
  · exact hL

/-- A convenient degree criterion for the remainder principle. -/
theorem quotientTransform_eq_modByMonic_of_degree
    (m : ℕ) (L A : Polynomial ℚ)
    (hLmonic : L.Monic)
    (hAdegree : A.natDegree < m + L.natDegree)
    (hL : ∀ k ∈ range (m + 1), L.eval (k : ℚ) ≠ 0) :
    quotientTransform m L A = quotientTransform m L (A %ₘ L) := by
  apply quotientTransform_eq_modByMonic_of_quotient_degree m L A
  · rw [Polynomial.natDegree_divByMonic A hLmonic]
    omega
  · exact hL

theorem mixedTemplate_eq_endpoint_mul_betaSum
    (m : ℕ) (A D : Polynomial ℚ) (c : ℚ)
    (hfactor :
      A = Polynomial.C c +
        (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) * D)
    (hD : D.natDegree < m) :
    mixedTemplateSum m A = c * betaSum m := by
  have hvanish := alternatingPolynomialSum_vanishes D m hD
  rw [mixedTemplateSum, betaSum]
  calc
    ∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (A.eval (k : ℚ) / (2 * k + 1))
      = ∑ k ∈ range (m + 1),
          ((-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ) +
            c * betaTerm m k) := by
              apply sum_congr rfl
              intro k _
              rw [hfactor]
              simp only [Polynomial.eval_add, Polynomial.eval_mul,
                Polynomial.eval_C, Polynomial.eval_X, betaTerm]
              have hden : (2 * (k : ℚ) + 1) ≠ 0 := by positivity
              field_simp [hden]
              ring
    _ = (∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ)) +
        ∑ k ∈ range (m + 1), c * betaTerm m k := by
          rw [sum_add_distrib]
    _ = c * ∑ k ∈ range (m + 1), betaTerm m k := by
      rw [hvanish, zero_add, mul_sum]

theorem mixedTemplate_eq_betaSum (m : ℕ) (A D : Polynomial ℚ)
    (hfactor :
      A = Polynomial.C 1 +
        (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) * D)
    (hD : D.natDegree < m) :
    mixedTemplateSum m A = betaSum m := by
  simpa using mixedTemplate_eq_endpoint_mul_betaSum m A D 1 hfactor hD

theorem mixedTemplate_eq_endpoint_mul_betaSum_of_degree
    (m : ℕ) (A : Polynomial ℚ) (hAdegree : A.natDegree < m) :
    mixedTemplateSum m A =
      A.eval (-1 / 2 : ℚ) * betaSum m := by
  obtain ⟨q, hq⟩ :=
    Polynomial.X_sub_C_dvd_sub_C_eval (p := A) (a := (-1 / 2 : ℚ))
  have hq_factor :
      A - Polynomial.C (A.eval (-1 / 2 : ℚ)) =
        (Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q := by
    exact hq
  let D : Polynomial ℚ := Polynomial.C (1 / 2 : ℚ) * q
  have hfactor :
      A = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
        (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) * D := by
    dsimp [D]
    calc
      A = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
          (A - Polynomial.C (A.eval (-1 / 2 : ℚ))) := by ring
      _ = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
          (Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q := by
            rw [hq_factor]
      _ = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
          (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) *
            (Polynomial.C (1 / 2 : ℚ) * q) := by
              have htwo :
                  Polynomial.C (2 : ℚ) * Polynomial.C (1 / 2 : ℚ) = 1 := by
                norm_num [← Polynomial.C_mul]
              have hminus :
                  Polynomial.C (-1 / 2 : ℚ) =
                    -Polynomial.C (1 / 2 : ℚ) := by norm_num
              simp only [Polynomial.C_1]
              rw [hminus]
              ring_nf
              have htwo' :
                  Polynomial.C (1 / 2 : ℚ) * Polynomial.C (2 : ℚ) = 1 := by
                rw [mul_comm, htwo]
              rw [show
                Polynomial.X * Polynomial.C (1 / 2 : ℚ) * q *
                    Polynomial.C (2 : ℚ) =
                  Polynomial.X * q *
                    (Polynomial.C (1 / 2 : ℚ) * Polynomial.C (2 : ℚ)) by
                      ring,
                htwo', mul_one]
  have hqdegree : q.natDegree ≤ A.natDegree := by
    by_cases hqzero : q = 0
    · simp [hqzero]
    · have hlinear :
          Polynomial.X - Polynomial.C (-1 / 2 : ℚ) ≠ 0 := by
        intro h
        have := congrArg (fun P : Polynomial ℚ ↦ P.coeff 1) h
        norm_num at this
      have hq_le_product :
          q.natDegree ≤
            ((Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q).natDegree := by
        rw [Polynomial.natDegree_mul hlinear hqzero]
        omega
      have hproduct_le :
          ((Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q).natDegree
            ≤ A.natDegree := by
        rw [← hq_factor]
        simpa using Polynomial.natDegree_sub_le A
          (Polynomial.C (A.eval (-1 / 2 : ℚ)))
      exact hq_le_product.trans hproduct_le
  have hDdegree : D.natDegree < m := by
    apply lt_of_le_of_lt _ hAdegree
    exact (Polynomial.natDegree_C_mul_le (1 / 2 : ℚ) q).trans hqdegree
  exact mixedTemplate_eq_endpoint_mul_betaSum m A D
    (A.eval (-1 / 2 : ℚ)) hfactor hDdegree

theorem mixedTemplate_eq_betaSum_of_endpoint (m : ℕ) (A : Polynomial ℚ)
    (hAdegree : A.natDegree < m)
    (hAendpoint : A.eval (-1 / 2 : ℚ) = 1) :
    mixedTemplateSum m A = betaSum m := by
  rw [mixedTemplate_eq_endpoint_mul_betaSum_of_degree m A hAdegree,
    hAendpoint, one_mul]

/-- The normalized product `A_{n+1}` from the SIC2C4 finite sums. -/
noncomputable def coefficientPoly : ℕ → Polynomial ℚ
  | 0 => 1
  | n + 1 =>
      Polynomial.C (1 / (n + 1 : ℚ)) *
        (Polynomial.C 2 * Polynomial.X +
          Polynomial.C ((n + 2 : ℕ) : ℚ)) *
          coefficientPoly n

-- The nested polynomial degree bounds require additional elaboration time.
set_option maxHeartbeats 800000 in
lemma coefficientPoly_degree (n : ℕ) :
    (coefficientPoly n).natDegree ≤ n := by
  induction n with
  | zero => simp [coefficientPoly]
  | succ n inductionHypothesis =>
      rw [coefficientPoly]
      calc
        (Polynomial.C (1 / (n + 1 : ℚ)) *
            (Polynomial.C 2 * Polynomial.X +
              Polynomial.C ((n + 2 : ℕ) : ℚ)) *
              coefficientPoly n).natDegree
          ≤ (Polynomial.C 2 * Polynomial.X +
                Polynomial.C ((n + 2 : ℕ) : ℚ)).natDegree +
              (coefficientPoly n).natDegree := by
                exact Polynomial.natDegree_mul_le.trans
                  (Nat.add_le_add_right
                    (Polynomial.natDegree_C_mul_le
                      (1 / (n + 1 : ℚ))
                      ((Polynomial.C 2 * Polynomial.X +
                        Polynomial.C ((n + 2 : ℕ) : ℚ)) * coefficientPoly n))
                    0)
        _ ≤ 1 + n := by
          apply Nat.add_le_add
          · calc
              (Polynomial.C 2 * Polynomial.X +
                  Polynomial.C ((n + 2 : ℕ) : ℚ)).natDegree
                ≤ max (Polynomial.C 2 * Polynomial.X).natDegree
                    (Polynomial.C ((n + 2 : ℕ) : ℚ)).natDegree :=
                      Polynomial.natDegree_add_le _ _
              _ ≤ 1 := by norm_num
          · exact inductionHypothesis
        _ = n + 1 := by omega

lemma coefficientPoly_endpoint (n : ℕ) :
    (coefficientPoly n).eval (-1 / 2 : ℚ) = 1 := by
  induction n with
  | zero => simp [coefficientPoly]
  | succ n inductionHypothesis =>
      rw [coefficientPoly]
      simp only [Polynomial.eval_mul, Polynomial.eval_add, Polynomial.eval_C,
        Polynomial.eval_X, inductionHypothesis]
      have hnonzero : (n + 1 : ℚ) ≠ 0 := by positivity
      field_simp [hnonzero]
      push_cast
      ring

/-- At a nonnegative integer, the recursive product polynomial is exactly
the binomial coefficient occurring in the chart constant term. -/
theorem coefficientPoly_eval_nat (n k : ℕ) :
    (coefficientPoly n).eval (k : ℚ) =
      ((n + 1 + 2 * k).choose n : ℚ) := by
  induction n with
  | zero => simp [coefficientPoly]
  | succ n inductionHypothesis =>
      rw [coefficientPoly]
      simp only [Polynomial.eval_mul, Polynomial.eval_add, Polynomial.eval_C,
        Polynomial.eval_X, inductionHypothesis]
      have hden : (n + 1 : ℚ) ≠ 0 := by positivity
      have hchooseNat := Nat.succ_mul_choose_eq (n + 1 + 2 * k) n
      have hchoose :
          (n + 2 + 2 * k : ℚ) *
              ((n + 1 + 2 * k).choose n : ℚ) =
            ((n + 2 + 2 * k).choose (n + 1) : ℚ) *
              (n + 1 : ℚ) := by
        exact_mod_cast hchooseNat
      field_simp [hden]
      push_cast
      linear_combination hchoose

/-- The pure finite sum in product-polynomial form vanishes. -/
theorem pureProductSum_vanishes (m : ℕ) (hm : 0 < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        (Polynomial.C (1 / (m : ℚ)) *
          coefficientPoly (m - 1)).eval (k : ℚ)) = 0 := by
  apply alternatingPolynomialSum_vanishes
  have hdegree :
      (Polynomial.C (1 / (m : ℚ)) *
        coefficientPoly (m - 1)).natDegree ≤ m - 1 := by
    exact (Polynomial.natDegree_C_mul_le
      (1 / (m : ℚ)) (coefficientPoly (m - 1))).trans
        (coefficientPoly_degree (m - 1))
  omega

/-- The pure displayed binomial quotient sum vanishes. -/
theorem pureBinomialSum_vanishes (m : ℕ) (hm : 0 < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        (((m + 2 * k).choose m : ℚ) / (2 * k + 1))) = 0 := by
  rw [← pureProductSum_vanishes m hm]
  apply sum_congr rfl
  intro k _
  simp only [Polynomial.eval_mul, Polynomial.eval_C,
    coefficientPoly_eval_nat]
  have htop : m - 1 + 1 + 2 * k = m + 2 * k := by omega
  rw [htop]
  have hchooseNat := Nat.choose_succ_right_eq (m + 2 * k) (m - 1)
  have hpred : m - 1 + 1 = m := by omega
  have hsub : m + 2 * k - (m - 1) = 2 * k + 1 := by omega
  rw [hpred, hsub] at hchooseNat
  have hchoose :
      ((m + 2 * k).choose m : ℚ) * (m : ℚ) =
        ((m + 2 * k).choose (m - 1) : ℚ) * (2 * k + 1 : ℚ) := by
    exact_mod_cast hchooseNat
  have hmQ : (m : ℚ) ≠ 0 := by positivity
  have hodd : (2 * k + 1 : ℚ) ≠ 0 := by positivity
  have hratio :
      ((m + 2 * k).choose m : ℚ) / (2 * k + 1) =
        (1 / (m : ℚ)) * ((m + 2 * k).choose (m - 1) : ℚ) := by
    field_simp [hmQ, hodd]
    exact hchoose
  rw [hratio]

/-- The mixed finite sum in product-polynomial form is the residual beta sum. -/
theorem mixedProductSum_eq_betaSum (m : ℕ) (hm : 0 < m) :
    mixedTemplateSum m (coefficientPoly (m - 1)) = betaSum m := by
  apply mixedTemplate_eq_betaSum_of_endpoint
  · exact lt_of_le_of_lt (coefficientPoly_degree (m - 1)) (by omega)
  · exact coefficientPoly_endpoint (m - 1)

/-- The mixed displayed binomial sum is the residual beta sum. -/
theorem mixedBinomialSum_eq_betaSum (m : ℕ) (hm : 0 < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        (((m + 2 * k).choose (m - 1) : ℚ) / (2 * k + 1))) =
      betaSum m := by
  rw [← mixedProductSum_eq_betaSum m hm, mixedTemplateSum]
  apply sum_congr rfl
  intro k _
  rw [coefficientPoly_eval_nat]
  have htop : m - 1 + 1 + 2 * k = m + 2 * k := by omega
  rw [htop]

lemma betaSum_zero : betaSum 0 = 1 := by
  norm_num [betaSum, betaTerm]

lemma choose_predecessor_mul (n k : ℕ) :
    (n.choose k : ℚ) * (n + 1 : ℕ) =
      ((n + 1).choose k : ℚ) * ((n + 1 - k : ℕ) : ℚ) := by
  exact_mod_cast Nat.choose_mul_succ_eq n k

/-- Triangular creative-telescoping recurrence for every repeated-pole
beta sum; equation (4.8c) in the written proof. -/
theorem repeatedBetaSum_recurrence (n s : ℕ) :
    (2 * (n + 1) + 1 : ℚ) * repeatedBetaSum (n + 1) (s + 1) =
      repeatedBetaSum (n + 1) s +
        (2 * (n + 1) : ℚ) * repeatedBetaSum n (s + 1) := by
  have h_extend :
      repeatedBetaSum n (s + 1) =
        ∑ k ∈ range (n + 2), repeatedBetaTerm n (s + 1) k := by
    calc
      repeatedBetaSum n (s + 1) =
          ∑ k ∈ range (n + 1), repeatedBetaTerm n (s + 1) k := rfl
      _ = (∑ k ∈ range (n + 1), repeatedBetaTerm n (s + 1) k) +
          repeatedBetaTerm n (s + 1) (n + 1) := by
            simp [repeatedBetaTerm]
      _ = ∑ k ∈ range ((n + 1) + 1),
          repeatedBetaTerm n (s + 1) k :=
            (sum_range_succ
              (fun k ↦ repeatedBetaTerm n (s + 1) k) (n + 1)).symm
      _ = ∑ k ∈ range (n + 2), repeatedBetaTerm n (s + 1) k := by
        congr 2
  have h_term (k : ℕ) (hk : k ∈ range (n + 2)) :
      (2 * (n + 1) + 1 : ℚ) * repeatedBetaTerm (n + 1) (s + 1) k =
        repeatedBetaTerm (n + 1) s k +
          (2 * (n + 1) : ℚ) * repeatedBetaTerm n (s + 1) k := by
    have hk_le : k ≤ n + 1 := by
      simp only [mem_range] at hk
      omega
    have hchoose := choose_predecessor_mul n k
    rw [Nat.cast_sub hk_le] at hchoose
    have hinner :
        (2 * (n + 1) + 1 : ℚ) * ((n + 1).choose k : ℚ) -
            (2 * (n + 1) : ℚ) * (n.choose k : ℚ) =
          (2 * k + 1 : ℚ) * ((n + 1).choose k : ℚ) := by
      push_cast at hchoose ⊢
      linear_combination -2 * hchoose
    have hsplit :
        (2 * (n + 1) + 1 : ℚ) * ((n + 1).choose k : ℚ) =
          (2 * k + 1 : ℚ) * ((n + 1).choose k : ℚ) +
            (2 * (n + 1) : ℚ) * (n.choose k : ℚ) := by
      linarith [hinner]
    simp only [repeatedBetaTerm]
    calc
      (2 * (n + 1) + 1 : ℚ) *
          ((-1 : ℚ) ^ k * ((n + 1).choose k : ℚ) /
            (2 * k + 1) ^ (s + 1))
        = ((-1 : ℚ) ^ k / (2 * k + 1) ^ (s + 1)) *
            ((2 * (n + 1) + 1 : ℚ) *
              ((n + 1).choose k : ℚ)) := by ring
      _ = ((-1 : ℚ) ^ k / (2 * k + 1) ^ (s + 1)) *
            ((2 * k + 1 : ℚ) * ((n + 1).choose k : ℚ) +
              (2 * (n + 1) : ℚ) * (n.choose k : ℚ)) := by
                rw [hsplit]
      _ = (-1 : ℚ) ^ k * ((n + 1).choose k : ℚ) /
              (2 * k + 1) ^ s +
            (2 * (n + 1) : ℚ) *
              ((-1 : ℚ) ^ k * (n.choose k : ℚ) /
                (2 * k + 1) ^ (s + 1)) := by
                  have hden : (2 * k + 1 : ℚ) ≠ 0 := by positivity
                  rw [pow_succ]
                  field_simp [hden]
                  ring
  rw [repeatedBetaSum, repeatedBetaSum, h_extend, mul_sum, mul_sum,
    ← sum_add_distrib]
  apply sum_congr rfl
  intro k hk
  exact h_term k hk

/-- Creative-telescoping recurrence for the residual beta sum. -/
theorem betaSum_recurrence (n : ℕ) :
    (2 * (n + 1) + 1 : ℚ) * betaSum (n + 1) =
      (2 * (n + 1) : ℚ) * betaSum n := by
  have h := repeatedBetaSum_recurrence n 0
  rw [repeatedBetaSum_one, repeatedBetaSum_one,
    repeatedBetaSum_zero_right (n + 1) (by omega), zero_add] at h
  simpa using h

/-- The product determined by the beta recurrence. -/
def betaProduct : ℕ → ℚ
  | 0 => 1
  | n + 1 =>
      (2 * (n + 1) : ℚ) / (2 * (n + 1) + 1) * betaProduct n

/-- Closed product evaluation of the residual finite sum. -/
theorem betaSum_eq_betaProduct (m : ℕ) : betaSum m = betaProduct m := by
  induction m with
  | zero => simp [betaSum_zero, betaProduct]
  | succ n inductionHypothesis =>
      rw [betaProduct, ← inductionHypothesis]
      have h := betaSum_recurrence n
      have hden : (2 * (n + 1) + 1 : ℚ) ≠ 0 := by positivity
      rw [div_mul_eq_mul_div]
      apply (eq_div_iff hden).2
      nlinarith [h]

/-- Factorial form of the closed beta product used in the mixed moment. -/
theorem betaProduct_eq_factorial_ratio (m : ℕ) :
    betaProduct m =
      ((2 : ℚ) ^ m * (m.factorial : ℚ)) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  induction m with
  | zero => norm_num [betaProduct, Nat.doubleFactorial]
  | succ n inductionHypothesis =>
      rw [betaProduct, inductionHypothesis, pow_succ, Nat.factorial_succ]
      rw [show 2 * (n + 1) + 1 = (2 * n + 1) + 2 by omega,
        Nat.doubleFactorial_add_two]
      have hodd :
          (Nat.doubleFactorial (2 * n + 1) : ℚ) ≠ 0 := by positivity
      have hnext :
          ((2 * n + 1 + 2 : ℕ) : ℚ) ≠ 0 := by positivity
      field_simp [hodd, hnext]
      push_cast
      ring

/-- Explicit all-order evaluation of the residual beta sum. -/
theorem betaSum_eq_factorial_ratio (m : ℕ) :
    betaSum m =
      ((2 : ℚ) ^ m * (m.factorial : ℚ)) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  rw [betaSum_eq_betaProduct, betaProduct_eq_factorial_ratio]

/-- Equation (4.4), including its normalization, has value zero. -/
theorem pureNormalizedBinomialSum_vanishes (m : ℕ) (hm : 0 < m) :
    (1 / (2 : ℚ) ^ m) *
      (∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (((m + 2 * k).choose m : ℚ) / (2 * k + 1))) = 0 := by
  rw [pureBinomialSum_vanishes m hm, mul_zero]

/-- Equation (4.5), including its normalization and closed evaluation. -/
theorem mixedNormalizedBinomialSum_eq_factorial_ratio
    (m : ℕ) (hm : 0 < m) :
    (1 / (2 : ℚ) ^ m) *
      (∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (((m + 2 * k).choose (m - 1) : ℚ) / (2 * k + 1))) =
      (m.factorial : ℚ) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  rw [mixedBinomialSum_eq_betaSum m hm, betaSum_eq_factorial_ratio]
  have hpow : (2 : ℚ) ^ m ≠ 0 := by positivity
  have hdouble :
      (Nat.doubleFactorial (2 * m + 1) : ℚ) ≠ 0 := by positivity
  field_simp [hpow, hdouble]

end FiniteEtaleKeller.SIC2C4
