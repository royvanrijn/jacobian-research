/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.FieldTheory.Normal.Basic
import Mathlib.FieldTheory.PrimitiveElement
import Mathlib.LinearAlgebra.StdBasis
import Mathlib.NumberTheory.NumberField.DedekindZeta
import Mathlib.NumberTheory.Padics.PadicNumbers
import Mathlib.NumberTheory.SumPrimeReciprocals
import Mathlib.RingTheory.Etale.Field
import Mathlib.RingTheory.KrullDimension.Zero
import Mathlib.RingTheory.Spectrum.Prime.Topology
import Mathlib.RingTheory.TensorProduct.Finite
import Mathlib.RingTheory.TensorProduct.Maps
import Mathlib.RingTheory.TensorProduct.Nontrivial
import Mathlib.RingTheory.TensorProduct.Pi
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

/-!
# The moment-theoretic degree-four barrier

This module formalizes the algebraic and positive-moment core of the
degree-four Hasse barrier.  The analytic input is deliberately isolated:
the Dirichlet prime mean of the local point count is the number of connected
components.  Applying that same input to the tensor square supplies the
second moment.

No Chebotarev theorem, Galois group, or monogenic presentation occurs in this
deduction.
-/

noncomputable section

namespace FiniteEtaleKeller

open scoped TensorProduct

universe u

/-- The number of connected components of a finite étale affine scheme.
For a zero-dimensional reduced algebra, its prime spectrum is finite and
discrete, so its cardinality is exactly its component count. -/
def componentCount (A : Type*) [CommRing A] : ℕ :=
  Nat.card (PrimeSpectrum A)

/-- Component count is invariant under ring equivalence. -/
theorem componentCount_eq_of_ringEquiv
    (A B : Type*) [CommRing A] [CommRing B] (e : A ≃+* B) :
    componentCount A = componentCount B := by
  exact Nat.card_congr (PrimeSpectrum.comapEquiv e).toEquiv

/-- The component count of a finite product is the sum of the component
counts of its factors. -/
theorem componentCount_pi
    (I : Type*) [Fintype I] (A : I → Type*) [∀ i, CommRing (A i)]
    [∀ i, Finite (PrimeSpectrum (A i))] :
    componentCount (∀ i, A i) = ∑ i, componentCount (A i) := by
  classical
  rw [componentCount]
  calc
    Nat.card (PrimeSpectrum (∀ i, A i)) =
        Nat.card (Σ i, PrimeSpectrum (A i)) :=
      Nat.card_congr (PrimeSpectrum.sigmaToPiHomeo A).symm.toEquiv
    _ = ∑ i, componentCount (A i) := by
      rw [Nat.card_sigma]
      rfl

/-- Every nontrivial ring with finite prime spectrum has at least one
connected component. -/
theorem one_le_componentCount
    (A : Type*) [CommRing A] [Nontrivial A] [Finite (PrimeSpectrum A)] :
    1 ≤ componentCount A := by
  exact Finite.card_pos

/-- For an explicit finite product of fields, `componentCount` is the number
of field factors.  This is the concrete finite-étale component adapter used
by the zeta--moment argument. -/
theorem componentCount_eq_natCard_of_algEquiv_pi_fields
    (K A I : Type*) [Field K] [CommRing A] [Algebra K A]
    [Finite I] (Ai : I → Type*) [∀ i, Field (Ai i)]
    [∀ i, Algebra K (Ai i)] (e : A ≃ₐ[K] ∀ i, Ai i) :
    componentCount A = Nat.card I := by
  classical
  letI := Fintype.ofFinite I
  rw [componentCount]
  calc
    Nat.card (PrimeSpectrum A) =
        Nat.card (PrimeSpectrum (∀ i, Ai i)) :=
      Nat.card_congr (PrimeSpectrum.comapEquiv e.toRingEquiv).toEquiv
    _ = Nat.card (Σ i, PrimeSpectrum (Ai i)) :=
      Nat.card_congr (PrimeSpectrum.sigmaToPiHomeo Ai).symm.toEquiv
    _ = Nat.card I := by
      rw [Nat.card_sigma]
      simp

/-- A reduced finite-spectrum algebra admitting a noninjective map to a field
has at least two connected components.  If its spectrum had at most one
point, reducedness would make the algebra itself a field, forcing every map
to a nontrivial ring to be injective. -/
theorem two_le_componentCount_of_noninjective_algHom
    (K B F : Type*) [Field K] [CommRing B] [Nontrivial B] [Field F]
    [Algebra K B] [Algebra K F]
    [IsReduced B] [Finite (PrimeSpectrum B)]
    (f : B →ₐ[K] F) (hf : ¬ Function.Injective f) :
    2 ≤ componentCount B := by
  by_contra h
  have hcard : Nat.card (PrimeSpectrum B) ≤ 1 := by
    rw [show Nat.card (PrimeSpectrum B) = componentCount B by rfl]
    omega
  letI : Subsingleton (PrimeSpectrum B) :=
    Finite.card_le_one_iff_subsingleton.mp hcard
  have hfield : IsField B :=
    PrimeSpectrum.subsingleton_iff_isField_of_isReduced.mp inferInstance
  letI : Field B := hfield.toField
  exact hf f.injective

/-- Multiplication on a field extension, viewed as an algebra map out of its
tensor square. -/
def tensorMultiplication
    (K L : Type*) [Field K] [Field L] [Algebra K L] :
    L ⊗[K] L →ₐ[K] L :=
  Algebra.TensorProduct.productMap (AlgHom.id K L) (AlgHom.id K L)

@[simp]
theorem tensorMultiplication_tmul
    (K L : Type*) [Field K] [Field L] [Algebra K L] (x y : L) :
    tensorMultiplication K L (x ⊗ₜ y) = x * y :=
  rfl

/-- If `L/K` has degree greater than one, multiplication
`L ⊗[K] L → L` is not injective.  The proof is purely dimensional. -/
theorem tensorMultiplication_not_injective
    (K L : Type*) [Field K] [Field L] [Algebra K L]
    [Module.Finite K L] (hdeg : 1 < Module.finrank K L) :
    ¬ Function.Injective (tensorMultiplication K L) := by
  intro hinj
  have hle :=
    (tensorMultiplication K L).toLinearMap.finrank_le_finrank_of_injective hinj
  rw [Module.finrank_tensorProduct] at hle
  nlinarith

/-- The diagonal tensor block of a nontrivial finite separable field
extension has at least two connected components. -/
theorem two_le_componentCount_tensor_self_of_one_lt_finrank
    (K L : Type*) [Field K] [Field L] [Algebra K L]
    [Module.Finite K L] [Algebra.IsSeparable K L]
    (hdeg : 1 < Module.finrank K L) :
    2 ≤ componentCount (L ⊗[K] L) := by
  letI : Algebra.FormallyEtale K L :=
    Algebra.FormallyEtale.of_isSeparable K L
  letI : IsReduced (L ⊗[K] L) :=
    Algebra.FormallyUnramified.isReduced_of_field L (L ⊗[K] L)
  letI : IsArtinianRing (L ⊗[K] L) :=
    IsArtinianRing.of_finite K (L ⊗[K] L)
  exact two_le_componentCount_of_noninjective_algHom K (L ⊗[K] L) L
    (tensorMultiplication K L)
    (tensorMultiplication_not_injective K L hdeg)

/-- Distributing a tensor square over a finite product reduces its component
count to the sum of the component counts of all pairwise tensor blocks. -/
theorem componentCount_tensor_pi_eq_sum
    (K I : Type*) [Field K] [Fintype I]
    (A : I → Type*) [∀ i, CommRing (A i)] [∀ i, Algebra K (A i)]
    [∀ i, Module.Finite K (A i)] :
    componentCount ((∀ i, A i) ⊗[K] (∀ i, A i)) =
      ∑ j, ∑ i, componentCount (A j ⊗[K] A i) := by
  classical
  letI (j : I) : IsArtinianRing ((∀ i, A i) ⊗[K] A j) :=
    IsArtinianRing.of_finite K _
  letI (j : I) : IsArtinianRing (A j ⊗[K] (∀ i, A i)) :=
    IsArtinianRing.of_finite K _
  letI (j i : I) : IsArtinianRing (A j ⊗[K] A i) :=
    IsArtinianRing.of_finite K _
  calc
    componentCount ((∀ i, A i) ⊗[K] (∀ i, A i)) =
        componentCount (∀ j, (∀ i, A i) ⊗[K] A j) :=
      componentCount_eq_of_ringEquiv _ _
        (Algebra.TensorProduct.piRight K K (∀ i, A i) A).toRingEquiv
    _ = ∑ j, componentCount ((∀ i, A i) ⊗[K] A j) :=
      componentCount_pi I _
    _ = ∑ j, componentCount (A j ⊗[K] (∀ i, A i)) := by
      apply Finset.sum_congr rfl
      intro j _
      exact componentCount_eq_of_ringEquiv _ _
        (Algebra.TensorProduct.comm K (∀ i, A i) (A j)).toRingEquiv
    _ = ∑ j, componentCount (∀ i, A j ⊗[K] A i) := by
      apply Finset.sum_congr rfl
      intro j _
      exact componentCount_eq_of_ringEquiv _ _
        (Algebra.TensorProduct.piRight K K (A j) A).toRingEquiv
    _ = ∑ j, ∑ i, componentCount (A j ⊗[K] A i) := by
      apply Finset.sum_congr rfl
      intro j _
      exact componentCount_pi I _

/-- A product of nontrivial finite separable field extensions has the strict
tensor-component surplus used in the degree-four barrier.  Every pairwise
tensor block contributes at least one component, and each diagonal block
contributes at least two. -/
theorem componentCount_tensor_pi_fields_ge_sq_add
    (K I : Type*) [Field K] [Fintype I]
    (L : I → Type*) [∀ i, Field (L i)] [∀ i, Algebra K (L i)]
    [∀ i, Module.Finite K (L i)] [∀ i, Algebra.IsSeparable K (L i)]
    (hdeg : ∀ i, 1 < Module.finrank K (L i)) :
    Fintype.card I ^ 2 + Fintype.card I ≤
      componentCount ((∀ i, L i) ⊗[K] (∀ i, L i)) := by
  classical
  letI (j i : I) : IsArtinianRing (L j ⊗[K] L i) :=
    IsArtinianRing.of_finite K _
  rw [componentCount_tensor_pi_eq_sum K I L]
  have hinner (j : I) :
      (∑ i : I, ((1 : ℕ) + if i = j then 1 else 0)) =
        Fintype.card I + 1 := by
    simp [Finset.sum_add_distrib]
  calc
    Fintype.card I ^ 2 + Fintype.card I =
        ∑ j, ∑ i, ((1 : ℕ) + if i = j then 1 else 0) := by
      symm
      calc
        (∑ j, ∑ i, ((1 : ℕ) + if i = j then 1 else 0)) =
            ∑ _j : I, (Fintype.card I + 1) :=
          Finset.sum_congr rfl (fun j _ ↦ hinner j)
        _ = Fintype.card I ^ 2 + Fintype.card I := by
          simp [pow_two, Nat.mul_add]
    _ ≤ ∑ j, ∑ i, componentCount (L j ⊗[K] L i) := by
      apply Finset.sum_le_sum
      intro j _
      apply Finset.sum_le_sum
      intro i _
      by_cases hij : i = j
      · subst i
        simpa using
          two_le_componentCount_tensor_self_of_one_lt_finrank K (L j) (hdeg j)
      · simp only [hij, if_false, add_zero]
        exact one_le_componentCount (L j ⊗[K] L i)

/-- Tensor products transport algebra equivalences in both variables. -/
def tensorAlgEquiv
    (K A B A' B' : Type*) [CommRing K]
    [CommRing A] [CommRing B] [CommRing A'] [CommRing B']
    [Algebra K A] [Algebra K B] [Algebra K A'] [Algebra K B']
    (eA : A ≃ₐ[K] A') (eB : B ≃ₐ[K] B') :
    A ⊗[K] B ≃ₐ[K] A' ⊗[K] B' :=
  AlgEquiv.ofAlgHom
    (Algebra.TensorProduct.map eA.toAlgHom eB.toAlgHom)
    (Algebra.TensorProduct.map eA.symm.toAlgHom eB.symm.toAlgHom)
    (by
      ext <;> simp)
    (by
      ext <;> simp)

/-- If a product of finite field extensions has no `K`-point, every field
factor has degree greater than one. -/
theorem one_lt_finrank_factors_of_isEmpty_algHom
    (K A I : Type*) [Field K] [CommRing A] [Algebra K A]
    [Finite I] (L : I → Type*) [∀ i, Field (L i)]
    [∀ i, Algebra K (L i)] [∀ i, Module.Finite K (L i)]
    (e : A ≃ₐ[K] ∀ i, L i) (hno : IsEmpty (A →ₐ[K] K)) :
    ∀ i, 1 < Module.finrank K (L i) := by
  classical
  letI := Fintype.ofFinite I
  intro i
  by_contra h
  have hone : Module.finrank K (L i) = 1 := by
    have hpos : 0 < Module.finrank K (L i) := Module.finrank_pos
    omega
  have hb : Function.Bijective (algebraMap K (L i)) :=
    Algebra.finrank_eq_one_iff_bijective_algebraMap.mp hone
  let ei : L i ≃ₐ[K] K :=
    (AlgEquiv.ofBijective (Algebra.ofId K (L i)) hb).symm
  exact hno.false <|
    ei.toAlgHom.comp ((Pi.evalAlgHom K L i).comp e.toAlgHom)

/-- Concrete tensor-component surplus for any supplied decomposition of a
finite étale algebra into nontrivial finite separable field factors. -/
theorem componentCount_tensor_ge_sq_add_of_algEquiv_pi_fields
    (K A I : Type*) [Field K] [CommRing A] [Algebra K A]
    [Finite I] (L : I → Type*) [∀ i, Field (L i)]
    [∀ i, Algebra K (L i)] [∀ i, Module.Finite K (L i)]
    [∀ i, Algebra.IsSeparable K (L i)]
    (e : A ≃ₐ[K] ∀ i, L i)
    (hdeg : ∀ i, 1 < Module.finrank K (L i)) :
    componentCount A ^ 2 + componentCount A ≤
      componentCount (A ⊗[K] A) := by
  classical
  letI := Fintype.ofFinite I
  rw [componentCount_eq_natCard_of_algEquiv_pi_fields K A I L e,
    Nat.card_eq_fintype_card]
  calc
    Fintype.card I ^ 2 + Fintype.card I ≤
        componentCount ((∀ i, L i) ⊗[K] (∀ i, L i)) :=
      componentCount_tensor_pi_fields_ge_sq_add K I L hdeg
    _ = componentCount (A ⊗[K] A) :=
      (componentCount_eq_of_ringEquiv _ _
        (tensorAlgEquiv K A A (∀ i, L i) (∀ i, L i) e e).toRingEquiv).symm

/-- The tensor-component surplus for an arbitrary finite étale algebra with
no rational point.  Mathlib's finite-étale field-product decomposition
supplies all field factors internally. -/
theorem componentCount_tensor_ge_sq_add_of_etale_isEmpty_algHom
    (K A : Type*) [Field K] [CommRing A] [Algebra K A]
    [Algebra.Etale K A] (hno : IsEmpty (A →ₐ[K] K)) :
    componentCount A ^ 2 + componentCount A ≤
      componentCount (A ⊗[K] A) := by
  haveI : Module.Finite K A :=
    Algebra.FormallyUnramified.finite_of_free K A
  obtain ⟨I, hI, L, hField, hAlgebra, e, hseparable⟩ :=
    (Algebra.FormallyEtale.iff_exists_algEquiv_prod K A).mp inferInstance
  letI : Finite I := hI
  letI (i : I) : Field (L i) := hField i
  letI (i : I) : Algebra K (L i) := hAlgebra i
  letI (i : I) : Module.Finite K (L i) :=
    Module.Finite.of_surjective ((LinearMap.proj i).comp e.toLinearMap)
      ((Function.surjective_eval i).comp e.surjective)
  letI (i : I) : Algebra.IsSeparable K (L i) := hseparable i
  exact componentCount_tensor_ge_sq_add_of_algEquiv_pi_fields K A I L e
    (one_lt_finrank_factors_of_isEmpty_algHom K A I L e hno)

/-- The number of algebra maps from `A` to a test field `L`.  For a finite
étale algebra this is the number of local sheets over `L`. -/
def localPointCount
    (K A L : Type*) [Field K] [CommRing A] [Field L]
    [Algebra K A] [Algebra K L] : ℕ :=
  Nat.card (A →ₐ[K] L)

/-- Restrict an algebra map out of a finite product to a factor selected by a
primitive idempotent that maps to one. -/
def factorAlgHomOfPi
    (K I F : Type*) [Field K] [Finite I] [DecidableEq I]
    (L : I → Type*) [∀ i, Field (L i)] [∀ i, Algebra K (L i)]
    [Field F] [Algebra K F]
    (g : (∀ i, L i) →ₐ[K] F) (i : I)
    (hi : g (Pi.single i 1) = 1) :
    L i →ₐ[K] F := by
  classical
  letI := Fintype.ofFinite I
  apply AlgHom.ofLinearMap
    (g.toLinearMap.comp (LinearMap.single K L i))
  · simpa using hi
  · intro x y
    simp only [LinearMap.coe_comp, Function.comp_apply, AlgHom.toLinearMap_apply,
      LinearMap.coe_single]
    rw [← map_mul]
    simp [Pi.single_mul]

/-- A map from a nonempty finite product of fields to a field is supported on
at least one factor.  This is the elementary idempotent argument that replaces
any appeal to Galois theory in the low-degree local count. -/
theorem exists_factorAlgHom_of_piAlgHom
    (K I F : Type*) [Field K] [Finite I]
    (L : I → Type*) [∀ i, Field (L i)] [∀ i, Algebra K (L i)]
    [Field F] [Algebra K F]
    (g : (∀ i, L i) →ₐ[K] F) :
    ∃ i, Nonempty (L i →ₐ[K] F) := by
  classical
  letI := Fintype.ofFinite I
  have hcomplete :=
    (CompleteOrthogonalIdempotents.single L).map g.toRingHom
  have hi : ∃ i, g (Pi.single i 1) = 1 := by
    by_contra h
    have hz : ∀ i, g (Pi.single i 1) = 0 := by
      intro i
      exact (IsIdempotentElem.iff_eq_zero_or_one.mp (hcomplete.idem i)).resolve_right
        (fun hi ↦ h ⟨i, hi⟩)
    have honezero : (1 : F) = 0 := by
      rw [← hcomplete.complete]
      simp [hz]
    exact one_ne_zero honezero
  obtain ⟨i, hi⟩ := hi
  exact ⟨i, ⟨factorAlgHomOfPi K I F L g i hi⟩⟩

/-- Extend an embedding of one field factor to the whole finite étale algebra
by evaluation on that factor. -/
def extendFactorAlgHom
    (K A I F : Type*) [Field K] [CommRing A] [Algebra K A]
    [Finite I] (L : I → Type*) [∀ i, Field (L i)]
    [∀ i, Algebra K (L i)] [Field F] [Algebra K F]
    (e : A ≃ₐ[K] ∀ i, L i) (i : I) :
    (L i →ₐ[K] F) → (A →ₐ[K] F) :=
  fun f ↦ f.comp ((Pi.evalAlgHom K L i).comp e.toAlgHom)

/-- Extending maps from a fixed factor is injective. -/
theorem extendFactorAlgHom_injective
    (K A I F : Type*) [Field K] [CommRing A] [Algebra K A]
    [Finite I] (L : I → Type*) [∀ i, Field (L i)]
    [∀ i, Algebra K (L i)] [Field F] [Algebra K F]
    (e : A ≃ₐ[K] ∀ i, L i) (i : I) :
    Function.Injective (extendFactorAlgHom K A I F L e i) := by
  intro f g hfg
  ext x
  obtain ⟨a, ha⟩ :
      ∃ a, ((Pi.evalAlgHom K L i).comp e.toAlgHom) a = x :=
    ((Function.surjective_eval i).comp e.surjective) x
  rw [← ha]
  exact AlgHom.congr_fun hfg a

/-- Artin independence bounds the number of field-valued algebra maps by the
rank of the source. -/
theorem localPointCount_le_finrank
    (K A F : Type*) [Field K] [CommRing A] [Algebra K A]
    [Module.Free K A] [Module.Finite K A] [Field F] [Algebra K F] :
    localPointCount K A F ≤ Module.finrank K A :=
  card_algHom_le_finrank K A F

/-- A quadratic field factor that embeds in `F` contributes exactly two
`K`-embeddings into `F`. -/
theorem two_le_localPointCount_of_quadratic_factor
    (K A I F : Type*) [Field K] [CommRing A] [Algebra K A]
    [Module.Free K A] [Module.Finite K A]
    [Finite I] (L : I → Type*) [∀ i, Field (L i)]
    [∀ i, Algebra K (L i)] [∀ i, Module.Free K (L i)]
    [∀ i, Module.Finite K (L i)] [∀ i, Algebra.IsSeparable K (L i)]
    [Field F] [Algebra K F]
    (e : A ≃ₐ[K] ∀ i, L i) (i : I)
    (hdeg : Module.finrank K (L i) = 2)
    (f : L i →ₐ[K] F) :
    2 ≤ localPointCount K A F := by
  letI : Algebra.IsQuadraticExtension K (L i) :=
    { finrank_eq_two' := hdeg }
  have hsplits :
      ∀ x : L i, ((minpoly K x).map (algebraMap K F)).Splits := by
    intro x
    have hs :=
      (Normal.splits (inferInstance : Normal K (L i)) x).map f.toRingHom
    rw [Polynomial.map_map] at hs
    have hcomp :
        f.toRingHom.comp (algebraMap K (L i)) = algebraMap K F :=
      RingHom.ext fun r ↦ f.commutes r
    rw [hcomp] at hs
    exact hs
  have hcard : Nat.card (L i →ₐ[K] F) = 2 := by
    rw [AlgHom.natCard_of_splits K (L i) F hsplits, hdeg]
  letI : Finite (A →ₐ[K] F) := Finite.algHom K A F
  rw [← hcard]
  exact Nat.card_le_card_of_injective _
    (extendFactorAlgHom_injective K A I F L e i)

/-- The concrete low-degree local-sheet inequality.  If a finite product of
nontrivial field extensions has total degree at most four and has an
`F`-point, then the number of `F`-points is at least the number of global
components. -/
theorem componentCount_le_localPointCount_of_algEquiv_pi_fields_rank_le_four
    (K A I F : Type*) [Field K] [CommRing A] [Algebra K A]
    [Module.Free K A] [Module.Finite K A]
    [Finite I] (L : I → Type*) [∀ i, Field (L i)]
    [∀ i, Algebra K (L i)] [∀ i, Module.Free K (L i)]
    [∀ i, Module.Finite K (L i)] [∀ i, Algebra.IsSeparable K (L i)]
    [Field F] [Algebra K F]
    (e : A ≃ₐ[K] ∀ i, L i)
    (hdeg : ∀ i, 1 < Module.finrank K (L i))
    (hrank : Module.finrank K A ≤ 4)
    (g : A →ₐ[K] F) :
    componentCount A ≤ localPointCount K A F := by
  classical
  letI := Fintype.ofFinite I
  have hsum : (∑ i, Module.finrank K (L i)) ≤ 4 := by
    rw [← Module.finrank_pi_fintype, ← e.toLinearEquiv.finrank_eq]
    exact hrank
  have hcard_le_two : Fintype.card I ≤ 2 := by
    have htwice :
        2 * Fintype.card I ≤ ∑ i, Module.finrank K (L i) := by
      calc
        2 * Fintype.card I = ∑ _i : I, 2 := by simp [Nat.mul_comm]
        _ ≤ ∑ i, Module.finrank K (L i) :=
          Finset.sum_le_sum fun i _ ↦ Nat.succ_le_of_lt (hdeg i)
    omega
  have hnonempty : Nonempty (A →ₐ[K] F) := ⟨g⟩
  have hpositive : 0 < localPointCount K A F := by
    letI : Finite (A →ₐ[K] F) := Finite.algHom K A F
    exact Finite.card_pos
  rw [componentCount_eq_natCard_of_algEquiv_pi_fields K A I L e,
    Nat.card_eq_fintype_card]
  by_cases hcard : Fintype.card I ≤ 1
  · omega
  have hcard_two : Fintype.card I = 2 := by omega
  let gp : (∀ i, L i) →ₐ[K] F := g.comp e.symm.toAlgHom
  obtain ⟨i, ⟨fi⟩⟩ := exists_factorAlgHom_of_piAlgHom K I F L gp
  obtain ⟨j, hji⟩ :=
    Fintype.exists_ne_of_one_lt_card (hcard_two ▸ by decide) i
  have hpair :
      Module.finrank K (L i) + Module.finrank K (L j) ≤
        ∑ k, Module.finrank K (L k) := by
    calc
      Module.finrank K (L i) + Module.finrank K (L j) =
          ∑ k ∈ ({i, j} : Finset I), Module.finrank K (L k) := by
        simp [hji.symm]
      _ ≤ ∑ k, Module.finrank K (L k) :=
        Finset.sum_le_sum_of_subset (Finset.subset_univ _)
  have hdeg_i : Module.finrank K (L i) = 2 := by
    have hi2 : 2 ≤ Module.finrank K (L i) :=
      Nat.succ_le_of_lt (hdeg i)
    have hj2 : 2 ≤ Module.finrank K (L j) :=
      Nat.succ_le_of_lt (hdeg j)
    omega
  rw [hcard_two]
  exact two_le_localPointCount_of_quadratic_factor K A I F L e i hdeg_i fi

/-- The low-degree local-sheet inequality for an arbitrary finite étale
algebra, with Mathlib's field-product decomposition hidden from the caller. -/
theorem componentCount_le_localPointCount_of_etale_rank_le_four
    (K A F : Type*) [Field K] [CommRing A] [Algebra K A]
    [Algebra.Etale K A] [Field F] [Algebra K F]
    (hno : IsEmpty (A →ₐ[K] K))
    (hrank : Module.finrank K A ≤ 4)
    (g : A →ₐ[K] F) :
    componentCount A ≤ localPointCount K A F := by
  haveI : Module.Finite K A :=
    Algebra.FormallyUnramified.finite_of_free K A
  obtain ⟨I, hI, L, hField, hAlgebra, e, hseparable⟩ :=
    (Algebra.FormallyEtale.iff_exists_algEquiv_prod K A).mp inferInstance
  letI : Finite I := hI
  letI (i : I) : Field (L i) := hField i
  letI (i : I) : Algebra K (L i) := hAlgebra i
  letI (i : I) : Module.Finite K (L i) :=
    Module.Finite.of_surjective ((LinearMap.proj i).comp e.toLinearMap)
      ((Function.surjective_eval i).comp e.surjective)
  letI (i : I) : Algebra.IsSeparable K (L i) := hseparable i
  exact
    componentCount_le_localPointCount_of_algEquiv_pi_fields_rank_le_four
      K A I F L e
      (one_lt_finrank_factors_of_isEmpty_algHom K A I L e hno)
      hrank g

/-- Algebra maps from a tensor product to a commutative algebra are pairs of
algebra maps from its two factors. -/
def tensorAlgHomEquiv
    (K A B C : Type*) [CommRing K] [CommRing A] [CommRing B] [CommRing C]
    [Algebra K A] [Algebra K B] [Algebra K C] :
    ((A ⊗[K] B) →ₐ[K] C) ≃ (A →ₐ[K] C) × (B →ₐ[K] C) where
  toFun f :=
    (f.comp Algebra.TensorProduct.includeLeft,
      f.comp Algebra.TensorProduct.includeRight)
  invFun f :=
    Algebra.TensorProduct.lift f.1 f.2 (fun _ _ ↦ .all _ _)
  left_inv f := by
    ext <;> simp
  right_inv f := by
    ext <;> simp

/-- Local point counts multiply under tensor products. -/
theorem localPointCount_tensor
    (K A B L : Type*) [Field K] [CommRing A] [CommRing B] [Field L]
    [Algebra K A] [Algebra K B] [Algebra K L] :
    localPointCount K (A ⊗[K] B) L =
      localPointCount K A L * localPointCount K B L := by
  rw [localPointCount, localPointCount, localPointCount, ← Nat.card_prod]
  exact Nat.card_congr (tensorAlgHomEquiv K A B L)

/-- In particular, the local point count of the tensor square is the square
of the original local point count. -/
theorem localPointCount_tensor_self
    (K A L : Type*) [Field K] [CommRing A] [Field L]
    [Algebra K A] [Algebra K L] :
    localPointCount K (A ⊗[K] A) L =
      localPointCount K A L ^ 2 := by
  rw [localPointCount_tensor, pow_two]

/-- The normalized Dirichlet sum over rational primes used in the paper's
definition of prime mean. -/
def primeDirichletAverage (a : Nat.Primes → ℝ) (s : ℝ) : ℝ :=
  (∑' p : Nat.Primes, a p * ((p.1 : ℝ) ^ (-s))) /
    (-Real.log (s - 1))

/-- A function on rational primes has Dirichlet prime mean `m` when its
normalized prime Dirichlet sums tend to `m` as `s → 1` from the right. -/
def HasDirichletPrimeMean (a : Nat.Primes → ℝ) (m : ℝ) : Prop :=
  Filter.Tendsto (primeDirichletAverage a)
    (nhdsWithin (1 : ℝ) (Set.Ioi 1)) (nhds m)

/-- A uniform absolute bound on a real-valued function on the rational
primes.  Local sheet counts are bounded by the rank, so this is the natural
domain on which the prime Dirichlet sums are absolutely convergent for
`s > 1`. -/
def IsPrimeBounded (a : Nat.Primes → ℝ) : Prop :=
  ∃ C : ℝ, 0 ≤ C ∧ ∀ p, |a p| ≤ C

namespace IsPrimeBounded

theorem const (c : ℝ) : IsPrimeBounded (fun _ : Nat.Primes ↦ c) :=
  ⟨|c|, abs_nonneg c, fun _ ↦ le_rfl⟩

theorem add {a b : Nat.Primes → ℝ}
    (ha : IsPrimeBounded a) (hb : IsPrimeBounded b) :
    IsPrimeBounded (a + b) := by
  obtain ⟨A, hA, ha⟩ := ha
  obtain ⟨B, hB, hb⟩ := hb
  refine ⟨A + B, add_nonneg hA hB, fun p ↦ ?_⟩
  exact (abs_add_le _ _).trans (add_le_add (ha p) (hb p))

theorem neg {a : Nat.Primes → ℝ} (ha : IsPrimeBounded a) :
    IsPrimeBounded (-a) := by
  obtain ⟨A, hA, ha⟩ := ha
  exact ⟨A, hA, fun p ↦ by simpa using ha p⟩

theorem sub {a b : Nat.Primes → ℝ}
    (ha : IsPrimeBounded a) (hb : IsPrimeBounded b) :
    IsPrimeBounded (a - b) :=
  ha.add hb.neg

theorem smul (c : ℝ) {a : Nat.Primes → ℝ} (ha : IsPrimeBounded a) :
    IsPrimeBounded (c • a) := by
  obtain ⟨A, hA, ha⟩ := ha
  refine ⟨|c| * A, mul_nonneg (abs_nonneg c) hA, fun p ↦ ?_⟩
  simpa [abs_mul] using mul_le_mul_of_nonneg_left (ha p) (abs_nonneg c)

end IsPrimeBounded

/-- Absolute convergence of a bounded prime Dirichlet series to the right of
one. -/
theorem summable_primeDirichlet_of_bounded
    {a : Nat.Primes → ℝ} (ha : IsPrimeBounded a)
    {s : ℝ} (hs : 1 < s) :
    Summable (fun p : Nat.Primes ↦ a p * ((p.1 : ℝ) ^ (-s))) := by
  obtain ⟨C, hC, ha⟩ := ha
  have hw : Summable (fun p : Nat.Primes ↦ (p : ℝ) ^ (-s)) :=
    Nat.Primes.summable_rpow.mpr (by linarith)
  apply (hw.mul_left C).of_norm_bounded
  intro p
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (by positivity :
    0 ≤ (p.1 : ℝ) ^ (-s))]
  exact mul_le_mul_of_nonneg_right (ha p) (by positivity)

/-- Additivity of the normalized prime Dirichlet sum on bounded functions. -/
theorem primeDirichletAverage_add
    {a b : Nat.Primes → ℝ}
    (ha : IsPrimeBounded a) (hb : IsPrimeBounded b)
    {s : ℝ} (hs : 1 < s) :
    primeDirichletAverage (a + b) s =
      primeDirichletAverage a s + primeDirichletAverage b s := by
  have hsa := summable_primeDirichlet_of_bounded ha hs
  have hsb := summable_primeDirichlet_of_bounded hb hs
  unfold primeDirichletAverage
  rw [← add_div, ← hsa.tsum_add hsb]
  congr 1
  apply tsum_congr
  intro p
  simp only [Pi.add_apply]
  ring

/-- Homogeneity of the normalized prime Dirichlet sum on bounded
functions. -/
theorem primeDirichletAverage_smul
    (c : ℝ) {a : Nat.Primes → ℝ}
    (ha : IsPrimeBounded a) {s : ℝ} (hs : 1 < s) :
    primeDirichletAverage (c • a) s =
      c * primeDirichletAverage a s := by
  have hsa := summable_primeDirichlet_of_bounded ha hs
  unfold primeDirichletAverage
  calc
    (∑' p : Nat.Primes, (c • a) p * ((p.1 : ℝ) ^ (-s))) /
          (-Real.log (s - 1)) =
        (∑' p : Nat.Primes, c * (a p * ((p.1 : ℝ) ^ (-s)))) /
          (-Real.log (s - 1)) := by
            congr 1
            apply tsum_congr
            intro p
            simp [mul_assoc]
    _ = (c * ∑' p : Nat.Primes, a p * ((p.1 : ℝ) ^ (-s))) /
          (-Real.log (s - 1)) := by rw [hsa.tsum_mul_left]
    _ = c * ((∑' p : Nat.Primes, a p * ((p.1 : ℝ) ^ (-s))) /
          (-Real.log (s - 1))) := by ring

/-- Dirichlet prime means add on bounded functions. -/
theorem HasDirichletPrimeMean.add
    {a b : Nat.Primes → ℝ} {ma mb : ℝ}
    (ha : IsPrimeBounded a) (hb : IsPrimeBounded b)
    (hma : HasDirichletPrimeMean a ma)
    (hmb : HasDirichletPrimeMean b mb) :
    HasDirichletPrimeMean (a + b) (ma + mb) := by
  refine (Filter.Tendsto.add hma hmb).congr' ?_
  filter_upwards [self_mem_nhdsWithin] with s hs
  exact (primeDirichletAverage_add ha hb hs).symm

/-- Dirichlet prime means commute with real scaling on bounded functions. -/
theorem HasDirichletPrimeMean.smul
    (c : ℝ) {a : Nat.Primes → ℝ} {ma : ℝ}
    (ha : IsPrimeBounded a)
    (hma : HasDirichletPrimeMean a ma) :
    HasDirichletPrimeMean (c • a) (c * ma) := by
  refine (tendsto_const_nhds.mul hma).congr' ?_
  filter_upwards [self_mem_nhdsWithin] with s hs
  exact (primeDirichletAverage_smul c ha hs).symm

/-- Dirichlet prime means subtract on bounded functions. -/
theorem HasDirichletPrimeMean.sub
    {a b : Nat.Primes → ℝ} {ma mb : ℝ}
    (ha : IsPrimeBounded a) (hb : IsPrimeBounded b)
    (hma : HasDirichletPrimeMean a ma)
    (hmb : HasDirichletPrimeMean b mb) :
    HasDirichletPrimeMean (a - b) (ma - mb) := by
  have hneg := HasDirichletPrimeMean.smul (-1) hb hmb
  have hneg' : HasDirichletPrimeMean (-b) (-mb) := by
    simpa using hneg
  have hadd := HasDirichletPrimeMean.add ha hb.neg hma hneg'
  simpa [sub_eq_add_neg] using hadd

/-- Pointwise order is preserved by the normalized prime Dirichlet sum near
`s = 1`, where its logarithmic denominator is positive. -/
theorem primeDirichletAverage_mono
    {a b : Nat.Primes → ℝ}
    (ha : IsPrimeBounded a) (hb : IsPrimeBounded b)
    (hab : ∀ p, a p ≤ b p)
    {s : ℝ} (hs : 1 < s) (hs2 : s < 2) :
    primeDirichletAverage a s ≤ primeDirichletAverage b s := by
  have hsa := summable_primeDirichlet_of_bounded ha hs
  have hsb := summable_primeDirichlet_of_bounded hb hs
  have hsum :
      (∑' p : Nat.Primes, a p * ((p.1 : ℝ) ^ (-s))) ≤
        ∑' p : Nat.Primes, b p * ((p.1 : ℝ) ^ (-s)) := by
    exact hsa.tsum_le_tsum
      (fun p ↦ mul_le_mul_of_nonneg_right (hab p) (by positivity)) hsb
  have hden : 0 < -Real.log (s - 1) := by
    have hlog : Real.log (s - 1) < 0 :=
      Real.log_neg (sub_pos.mpr hs) (by linarith)
    linarith
  exact (div_le_div_iff_of_pos_right hden).2 hsum

/-- The Dirichlet prime mean is positive: pointwise order between bounded
functions passes to their means. -/
theorem HasDirichletPrimeMean.mono
    {a b : Nat.Primes → ℝ} {ma mb : ℝ}
    (ha : IsPrimeBounded a) (hb : IsPrimeBounded b)
    (hab : ∀ p, a p ≤ b p)
    (hma : HasDirichletPrimeMean a ma)
    (hmb : HasDirichletPrimeMean b mb) :
    ma ≤ mb := by
  apply le_of_tendsto_of_tendsto hma hmb
  filter_upwards [self_mem_nhdsWithin,
    (eventually_lt_nhds (by norm_num : (1 : ℝ) < 2)).filter_mono
      nhdsWithin_le_nhds] with s hs hs2
  exact primeDirichletAverage_mono ha hb hab hs hs2

/-- The zero-variance calculation for the actual Dirichlet prime mean.
Unlike `PositiveNormalizedMean.second_moment_eq_sq_of_bounds`, this theorem
works directly with the normalized prime sums and proves their required
linearity and positivity from absolute convergence. -/
theorem second_moment_eq_sq_of_dirichletPrimeMean
    (ν νsq : Nat.Primes → ℝ) (r d second : ℝ)
    (hνbounded : IsPrimeBounded ν)
    (hνsqbounded : IsPrimeBounded νsq)
    (hνsq : ∀ p, νsq p = (ν p) ^ 2)
    (hr : 0 ≤ r)
    (hlower : ∀ p, r ≤ ν p)
    (hupper : ∀ p, ν p ≤ d)
    (hone : HasDirichletPrimeMean (fun _ : Nat.Primes ↦ (1 : ℝ)) 1)
    (hfirst : HasDirichletPrimeMean ν r)
    (hsecond : HasDirichletPrimeMean νsq second) :
    second = r ^ 2 := by
  let one : Nat.Primes → ℝ := fun _ ↦ 1
  let δ : Nat.Primes → ℝ := ν - r • one
  let variance : Nat.Primes → ℝ := νsq - (r ^ 2) • one
  have honebounded : IsPrimeBounded one := by
    simpa [one] using IsPrimeBounded.const (1 : ℝ)
  have hronebounded : IsPrimeBounded (r • one) :=
    honebounded.smul r
  have hrsqonebounded : IsPrimeBounded ((r ^ 2) • one) :=
    honebounded.smul (r ^ 2)
  have hδbounded : IsPrimeBounded δ :=
    hνbounded.sub hronebounded
  have hvariancebounded : IsPrimeBounded variance :=
    hνsqbounded.sub hrsqonebounded
  have honemean : HasDirichletPrimeMean one 1 := by
    simpa [one] using hone
  have hδmean : HasDirichletPrimeMean δ 0 := by
    have hrmean :=
      HasDirichletPrimeMean.smul r honebounded honemean
    have hsub :=
      HasDirichletPrimeMean.sub hνbounded hronebounded hfirst hrmean
    simpa [δ] using hsub
  have hvariancemean :
      HasDirichletPrimeMean variance (second - r ^ 2) := by
    have hrsqmean :=
      HasDirichletPrimeMean.smul (r ^ 2) honebounded honemean
    exact HasDirichletPrimeMean.sub hνsqbounded hrsqonebounded
      hsecond (by simpa using hrsqmean)
  have hzeroBounded :
      IsPrimeBounded (0 : Nat.Primes → ℝ) := by
    convert IsPrimeBounded.const 0 using 1
    ext p
    simp
  have hzeroMean :
      HasDirichletPrimeMean (0 : Nat.Primes → ℝ) 0 := by
    have := HasDirichletPrimeMean.smul 0 honebounded honemean
    convert this using 1 <;> simp
  have hvariance_nonneg :
      ∀ p, (0 : ℝ) ≤ variance p := by
    intro p
    dsimp [variance, one]
    rw [hνsq p]
    nlinarith [hlower p]
  have hvariance_upper :
      ∀ p, variance p ≤ ((d + r) • δ) p := by
    intro p
    dsimp [variance, δ, one]
    rw [hνsq p]
    have h₁ : 0 ≤ ν p - r := sub_nonneg.mpr (hlower p)
    have h₂ : 0 ≤ d - ν p := sub_nonneg.mpr (hupper p)
    nlinarith [mul_nonneg h₁ h₂]
  have hscaledδbounded : IsPrimeBounded ((d + r) • δ) :=
    hδbounded.smul (d + r)
  have hscaledδmean :
      HasDirichletPrimeMean ((d + r) • δ) 0 := by
    have := HasDirichletPrimeMean.smul (d + r) hδbounded hδmean
    simpa using this
  have hnonneg : 0 ≤ second - r ^ 2 :=
    HasDirichletPrimeMean.mono hzeroBounded hvariancebounded
      hvariance_nonneg hzeroMean hvariancemean
  have hnonpos : second - r ^ 2 ≤ 0 :=
    HasDirichletPrimeMean.mono hvariancebounded hscaledδbounded
      hvariance_upper hvariancemean hscaledδmean
  linarith

/-- The local-sheet count at the rational prime represented by `p`. -/
def rationalPrimeLocalPointCount
    (A : Type*) [CommRing A] [Algebra ℚ A] (p : Nat.Primes) : ℝ :=
  letI : Fact p.1.Prime := ⟨p.2⟩
  localPointCount ℚ A ℚ_[p.1]

/-- The rational base field has one connected component. -/
theorem componentCount_rat : componentCount ℚ = 1 := by
  rw [componentCount]
  exact Nat.card_unique

/-- The rational base field has exactly one local sheet at every rational
prime. -/
theorem rationalPrimeLocalPointCount_rat (p : Nat.Primes) :
    rationalPrimeLocalPointCount ℚ p = 1 := by
  letI : Fact p.1.Prime := ⟨p.2⟩
  change (Nat.card (ℚ →ₐ[ℚ] ℚ_[p.1]) : ℝ) = 1
  norm_num

/-- Solubility over every rational `p`-adic field, indexed without an
ambient typeclass assumption on the prime. -/
def RationalPrimeLocallySoluble
    (A : Type*) [CommRing A] [Algebra ℚ A] : Prop :=
  ∀ p : Nat.Primes,
    letI : Fact p.1.Prime := ⟨p.2⟩
    Nonempty (A →ₐ[ℚ] ℚ_[p.1])

/-- Rational local-sheet counts are uniformly bounded by the algebra rank. -/
theorem rationalPrimeLocalPointCount_bounded
    (A : Type*) [CommRing A] [Algebra ℚ A]
    [Module.Free ℚ A] [Module.Finite ℚ A] :
    IsPrimeBounded (rationalPrimeLocalPointCount A) := by
  refine ⟨Module.finrank ℚ A, by positivity, fun p ↦ ?_⟩
  letI : Fact p.1.Prime := ⟨p.2⟩
  change |(localPointCount ℚ A ℚ_[p.1] : ℝ)| ≤ Module.finrank ℚ A
  rw [abs_of_nonneg (by positivity)]
  exact_mod_cast localPointCount_le_finrank ℚ A ℚ_[p.1]

/-- At every rational prime, the local-sheet count of the tensor square is
the square of the original count. -/
theorem rationalPrimeLocalPointCount_tensor_self
    (A : Type*) [CommRing A] [Algebra ℚ A] (p : Nat.Primes) :
    rationalPrimeLocalPointCount (A ⊗[ℚ] A) p =
      (rationalPrimeLocalPointCount A p) ^ 2 := by
  letI : Fact p.1.Prime := ⟨p.2⟩
  simp only [rationalPrimeLocalPointCount, localPointCount_tensor_self]
  norm_cast

/-- The exact analytic theorem still required by the formal development:
the Dirichlet prime mean of the local-sheet count equals the number of
connected components for every finite étale rational algebra.  The explicit
normalization is mathematically its `A = ℚ` case; recording it separately
avoids irrelevant universe-lift bookkeeping. -/
def RationalFiniteEtalePrimeMomentStatement : Prop :=
  HasDirichletPrimeMean (fun _ : Nat.Primes ↦ (1 : ℝ)) 1 ∧
    ∀ (A : Type u) [CommRing A] [Algebra ℚ A] [Algebra.Etale ℚ A],
      HasDirichletPrimeMean
          (rationalPrimeLocalPointCount A) (componentCount A) ∧
        HasDirichletPrimeMean
          (rationalPrimeLocalPointCount (A ⊗[ℚ] A))
          (componentCount (A ⊗[ℚ] A))

/-- Mathlib already supplies the nonzero simple pole of every Dedekind zeta
function.  The remaining analytic work is the Euler-product coefficient
extraction connecting this theorem to
`RationalFiniteEtalePrimeMomentStatement`. -/
theorem dedekindZeta_simplePole_input
    (K : Type*) [Field K] [NumberField K] :
    Filter.Tendsto
      (fun s : ℝ ↦ (s - 1) * NumberField.dedekindZeta K s)
      (nhdsWithin 1 (Set.Ioi 1))
      (nhds (NumberField.dedekindZeta_residue K)) ∧
      NumberField.dedekindZeta_residue K ≠ 0 :=
  ⟨NumberField.tendsto_sub_one_mul_dedekindZeta_nhdsGT K,
    NumberField.dedekindZeta_residue_ne_zero K⟩

/-- A normalized positive mean on a vector space `V` of real-valued
functions.  Taking `V` to be the functions whose Dirichlet prime mean exists
produces exactly this interface; no mean is postulated for functions outside
that domain. -/
structure PositiveNormalizedMean
    (ι : Type*) (V : Submodule ℝ (ι → ℝ)) where
  one_mem : (1 : ι → ℝ) ∈ V
  toLinearMap : V →ₗ[ℝ] ℝ
  map_one : toLinearMap ⟨1, one_mem⟩ = 1
  monotone' : Monotone toLinearMap

namespace PositiveNormalizedMean

variable {ι : Type*} {V : Submodule ℝ (ι → ℝ)}

instance : CoeFun (PositiveNormalizedMean ι V) (fun _ ↦ V → ℝ) :=
  ⟨fun M ↦ M.toLinearMap⟩

/-- The constant function, regarded as an element of the mean's domain. -/
def const (M : PositiveNormalizedMean ι V) (a : ℝ) : V :=
  ⟨fun _ ↦ a, by
    convert V.smul_mem a M.one_mem using 1
    ext i
    simp⟩

@[simp]
theorem apply_one (M : PositiveNormalizedMean ι V) :
    M ⟨1, M.one_mem⟩ = 1 :=
  M.map_one

@[simp]
theorem map_const (M : PositiveNormalizedMean ι V) (a : ℝ) :
    M (M.const a) = a := by
  have heq : M.const a = a • (⟨1, M.one_mem⟩ : V) := by
    ext i
    simp [const]
  rw [heq, M.toLinearMap.map_smul, M.map_one, smul_eq_mul, mul_one]

@[simp]
theorem map_add (M : PositiveNormalizedMean ι V) (f g : V) :
    M (f + g) = M f + M g :=
  M.toLinearMap.map_add f g

@[simp]
theorem map_sub (M : PositiveNormalizedMean ι V) (f g : V) :
    M (f - g) = M f - M g :=
  M.toLinearMap.map_sub f g

@[simp]
theorem map_smul (M : PositiveNormalizedMean ι V) (a : ℝ) (f : V) :
    M (a • f) = a * M f := by
  rw [M.toLinearMap.map_smul, smul_eq_mul]

theorem monotone (M : PositiveNormalizedMean ι V) :
    Monotone M :=
  M.monotone'

/-- The half-page moment calculation used in the degree-four barrier.

If a bounded local sheet count is everywhere at least its first moment, then
its variance under any normalized positive mean is zero. -/
theorem second_moment_eq_sq_of_bounds
    (M : PositiveNormalizedMean ι V) (ν νsq : V)
    (r d second : ℝ)
    (hνsq : ∀ i, νsq.1 i = (ν.1 i) ^ 2)
    (hr : 0 ≤ r)
    (hlower : ∀ i, r ≤ ν.1 i)
    (hupper : ∀ i, ν.1 i ≤ d)
    (hfirst : M ν = r)
    (hsecond : M νsq = second) :
    second = r ^ 2 := by
  let δ : V := ν - M.const r
  let variance : V := νsq - M.const (r ^ 2)
  have hδmean : M δ = 0 := by
    change M (ν - M.const r) = 0
    rw [M.map_sub, hfirst, M.map_const, sub_self]
  have hvariance_nonneg : (0 : V) ≤ variance := by
    intro i
    simp only [variance, Submodule.coe_zero, Pi.zero_apply, Submodule.coe_sub,
      Pi.sub_apply, const]
    rw [hνsq]
    nlinarith [hlower i]
  have hvariance_upper :
      variance ≤ (d + r) • δ := by
    intro i
    simp only [variance, δ, Submodule.coe_sub, Pi.sub_apply, const,
      Submodule.coe_smul_of_tower, Pi.smul_apply, smul_eq_mul]
    rw [hνsq]
    have h₁ : 0 ≤ ν.1 i - r := sub_nonneg.mpr (hlower i)
    have h₂ : 0 ≤ d - ν.1 i := sub_nonneg.mpr (hupper i)
    have hprod : 0 ≤ (ν.1 i - r) * (d - ν.1 i) := mul_nonneg h₁ h₂
    nlinarith
  have hmean_variance_nonneg : 0 ≤ M variance := by
    have := M.monotone hvariance_nonneg
    simpa using this
  have hmean_variance_nonpos : M variance ≤ 0 := by
    have hle := M.monotone hvariance_upper
    calc
      M variance ≤ M ((d + r) • δ) := hle
      _ = (d + r) * M δ := M.map_smul _ _
      _ = 0 := by rw [hδmean, mul_zero]
  have hmean_variance : M variance = 0 :=
    le_antisymm hmean_variance_nonpos hmean_variance_nonneg
  have hvariance_eval : M variance = second - r ^ 2 := by
    change M (νsq - M.const (r ^ 2)) =
      second - r ^ 2
    rw [M.map_sub, hsecond, M.map_const]
  linarith

/-- Moment form of the degree-four contradiction.  A positive component
surplus in the tensor square is incompatible with a local sheet count that is
everywhere at least the number of global components. -/
theorem contradiction_of_component_surplus
    (M : PositiveNormalizedMean ι V) (ν νsq : V)
    (components rank tensorComponents : ℕ)
    (hνsq : ∀ i, νsq.1 i = (ν.1 i) ^ 2)
    (hcomponents : 0 < components)
    (hlower : ∀ i, (components : ℝ) ≤ ν.1 i)
    (hupper : ∀ i, ν.1 i ≤ rank)
    (hfirst : M ν = components)
    (hsecond : M νsq = tensorComponents)
    (hsurplus :
      components ^ 2 + components ≤ tensorComponents) :
    False := by
  have heq :
      (tensorComponents : ℝ) = (components : ℝ) ^ 2 :=
    second_moment_eq_sq_of_bounds M ν νsq components rank tensorComponents
      hνsq (by positivity) hlower hupper (by exact_mod_cast hfirst)
      (by exact_mod_cast hsecond)
  have hstrict :
      (components : ℝ) ^ 2 < tensorComponents := by
    exact_mod_cast
      (lt_of_lt_of_le
        (Nat.lt_add_of_pos_right hcomponents :
          components ^ 2 < components ^ 2 + components)
        hsurplus)
  linarith

end PositiveNormalizedMean

/-- End-to-end algebraic reduction of the degree-four Hasse barrier to the
two prime-moment identities.  All finite-étale decomposition, local-sheet
bounds, tensor-square counting, and strict component surplus are proved
inside Lean; the caller supplies only a positive normalized mean realizing
the first moments of `A` and `A ⊗ A`. -/
theorem no_rank_le_four_hasse_failure_of_moments
    (K A ι : Type*) [Field K] [CommRing A] [Nontrivial A] [Algebra K A]
    [Algebra.Etale K A]
    (F : ι → Type*) [∀ i, Field (F i)] [∀ i, Algebra K (F i)]
    (V : Submodule ℝ (ι → ℝ)) (M : PositiveNormalizedMean ι V)
    (ν νsq : V)
    (hno : IsEmpty (A →ₐ[K] K))
    (hrank : Module.finrank K A ≤ 4)
    (hlocal : ∀ i, Nonempty (A →ₐ[K] F i))
    (hν : ∀ i, ν.1 i = localPointCount K A (F i))
    (hνsq : ∀ i, νsq.1 i = localPointCount K (A ⊗[K] A) (F i))
    (hfirst : M ν = componentCount A)
    (hsecond : M νsq = componentCount (A ⊗[K] A)) :
    False := by
  haveI : Module.Finite K A :=
    Algebra.FormallyUnramified.finite_of_free K A
  letI : IsArtinianRing A := IsArtinianRing.of_finite K A
  apply PositiveNormalizedMean.contradiction_of_component_surplus M ν νsq
    (componentCount A) (Module.finrank K A)
    (componentCount (A ⊗[K] A))
  · intro i
    rw [hνsq i, localPointCount_tensor_self, hν i]
    norm_cast
  · exact one_le_componentCount A
  · intro i
    rw [hν i]
    exact_mod_cast
      componentCount_le_localPointCount_of_etale_rank_le_four
        K A (F i) hno hrank (Classical.choice (hlocal i))
  · intro i
    rw [hν i]
    exact_mod_cast localPointCount_le_finrank K A (F i)
  · exact_mod_cast hfirst
  · exact_mod_cast hsecond
  · exact componentCount_tensor_ge_sq_add_of_etale_isEmpty_algHom K A hno

/-- The isolated zeta--Burnside first-moment statement implies the complete
rational degree-four Hasse barrier directly.  No abstract mean functional is
assumed: linearity, positivity, and the second-moment calculation are all
proved above for the actual normalized Dirichlet prime sums. -/
theorem no_rank_le_four_hasse_failure_of_rationalPrimeMomentStatement
    (A : Type u) [CommRing A] [Nontrivial A] [Algebra ℚ A]
    [Algebra.Etale ℚ A]
    (hmoment :
      HasDirichletPrimeMean (fun _ : Nat.Primes ↦ (1 : ℝ)) 1 ∧
        ∀ (B : Type u) [CommRing B] [Algebra ℚ B] [Algebra.Etale ℚ B],
          HasDirichletPrimeMean
              (rationalPrimeLocalPointCount B) (componentCount B) ∧
            HasDirichletPrimeMean
              (rationalPrimeLocalPointCount (B ⊗[ℚ] B))
              (componentCount (B ⊗[ℚ] B)))
    (hno : IsEmpty (A →ₐ[ℚ] ℚ))
    (hrank : Module.finrank ℚ A ≤ 4)
    (hlocal : RationalPrimeLocallySoluble A) :
    False := by
  haveI : Module.Finite ℚ A :=
    Algebra.FormallyUnramified.finite_of_free ℚ A
  letI : IsArtinianRing A := IsArtinianRing.of_finite ℚ A
  let ν : Nat.Primes → ℝ := rationalPrimeLocalPointCount A
  let νsq : Nat.Primes → ℝ :=
    rationalPrimeLocalPointCount (A ⊗[ℚ] A)
  have hνbounded : IsPrimeBounded ν := by
    exact rationalPrimeLocalPointCount_bounded A
  have hνsqbounded : IsPrimeBounded νsq := by
    exact rationalPrimeLocalPointCount_bounded (A ⊗[ℚ] A)
  obtain ⟨hone, hmoment⟩ := hmoment
  have hfirst :
      HasDirichletPrimeMean ν (componentCount A) := by
    exact (hmoment A).1
  have hsecond :
      HasDirichletPrimeMean νsq (componentCount (A ⊗[ℚ] A)) := by
    exact (hmoment A).2
  have hsq : ∀ p, νsq p = (ν p) ^ 2 := by
    intro p
    exact rationalPrimeLocalPointCount_tensor_self A p
  have hlower : ∀ p, (componentCount A : ℝ) ≤ ν p := by
    intro p
    letI : Fact p.1.Prime := ⟨p.2⟩
    have hp : Nonempty (A →ₐ[ℚ] ℚ_[p.1]) := hlocal p
    change (componentCount A : ℝ) ≤ localPointCount ℚ A ℚ_[p.1]
    exact_mod_cast
      componentCount_le_localPointCount_of_etale_rank_le_four
        ℚ A ℚ_[p.1] hno hrank (Classical.choice hp)
  have hupper : ∀ p, ν p ≤ Module.finrank ℚ A := by
    intro p
    letI : Fact p.1.Prime := ⟨p.2⟩
    change (localPointCount ℚ A ℚ_[p.1] : ℝ) ≤ Module.finrank ℚ A
    exact_mod_cast localPointCount_le_finrank ℚ A ℚ_[p.1]
  have heq :
      (componentCount (A ⊗[ℚ] A) : ℝ) =
        (componentCount A : ℝ) ^ 2 :=
    second_moment_eq_sq_of_dirichletPrimeMean ν νsq
      (componentCount A) (Module.finrank ℚ A)
      (componentCount (A ⊗[ℚ] A))
      hνbounded hνsqbounded hsq (by positivity) hlower hupper
      hone hfirst hsecond
  have hsurplus :=
    componentCount_tensor_ge_sq_add_of_etale_isEmpty_algHom ℚ A hno
  have hpositive : 0 < componentCount A :=
    one_le_componentCount A
  have hstrict :
      (componentCount A : ℝ) ^ 2 <
        componentCount (A ⊗[ℚ] A) := by
    exact_mod_cast
      (lt_of_lt_of_le
        (Nat.lt_add_of_pos_right hpositive :
          componentCount A ^ 2 <
            componentCount A ^ 2 + componentCount A)
        hsurplus)
  linarith

#print axioms tensorAlgHomEquiv
#print axioms componentCount_eq_natCard_of_algEquiv_pi_fields
#print axioms two_le_componentCount_tensor_self_of_one_lt_finrank
#print axioms componentCount_tensor_pi_fields_ge_sq_add
#print axioms one_lt_finrank_factors_of_isEmpty_algHom
#print axioms componentCount_tensor_ge_sq_add_of_algEquiv_pi_fields
#print axioms componentCount_tensor_ge_sq_add_of_etale_isEmpty_algHom
#print axioms exists_factorAlgHom_of_piAlgHom
#print axioms two_le_localPointCount_of_quadratic_factor
#print axioms componentCount_le_localPointCount_of_etale_rank_le_four
#print axioms localPointCount_tensor_self
#print axioms PositiveNormalizedMean.second_moment_eq_sq_of_bounds
#print axioms PositiveNormalizedMean.contradiction_of_component_surplus
#print axioms no_rank_le_four_hasse_failure_of_moments
#print axioms second_moment_eq_sq_of_dirichletPrimeMean
#print axioms no_rank_le_four_hasse_failure_of_rationalPrimeMomentStatement

end FiniteEtaleKeller
