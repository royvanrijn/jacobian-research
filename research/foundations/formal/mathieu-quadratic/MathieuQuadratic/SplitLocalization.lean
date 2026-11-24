import MathieuQuadratic.Laurent
import Mathlib.Algebra.MvPolynomial.Eval
import Mathlib.RingTheory.Localization.Away.Basic

/-!
# Laurent coordinates on a split monomial complement

Inverting a nonzero scalar times a monomial with positive exponent in
each variable inverts each variable, and gives the Laurent polynomial
algebra. This covers the unary and split binary coordinate rings.
-/

open scoped BigOperators

namespace MathieuQuadratic

open MvPolynomial

variable {K J : Type*} [Field K] [Fintype J] [DecidableEq J]

abbrev LaurentAlgebra (K J : Type*) [Field K] :=
  AddMonoidAlgebra K (J → ℤ)

private def basisExponent (i : J) : J → ℤ := Pi.single i 1

private noncomputable def laurentMonomial (v : J → ℤ) : (LaurentAlgebra K J)ˣ where
  val := AddMonoidAlgebra.single v 1
  inv := AddMonoidAlgebra.single (-v) 1
  val_inv := by simp [AddMonoidAlgebra.one_def]
  inv_val := by simp [AddMonoidAlgebra.one_def]

private noncomputable def laurentMonomialHom :
    Multiplicative (J → ℤ) →* (LaurentAlgebra K J)ˣ where
  toFun v := laurentMonomial v.toAdd
  map_one' := by
    apply Units.ext
    rfl
  map_mul' v w := by
    apply Units.ext
    simp [laurentMonomial]

/-- The norm polynomial in Laurent coordinates, before localization. -/
noncomputable def splitNorm (c : Kˣ) (powers : J → ℕ) : MvPolynomial J K :=
  C (c : K) * ∏ i, X i ^ powers i

abbrev SplitComplement (c : Kˣ) (powers : J → ℕ) :=
  Localization.Away (splitNorm c powers)

private noncomputable def laurentEvaluation : MvPolynomial J K →ₐ[K] LaurentAlgebra K J :=
  aeval fun i ↦ (↑(laurentMonomial (K := K) (basisExponent i)) : LaurentAlgebra K J)

private theorem laurentEvaluation_norm_unit (c : Kˣ) (powers : J → ℕ) :
    IsUnit (laurentEvaluation (splitNorm c powers)) := by
  rw [splitNorm, map_mul, map_prod, laurentEvaluation, aeval_C]
  exact (c.isUnit.map (algebraMap K _)).mul <|
    IsUnit.prod_univ_iff.mpr fun i ↦ by rw [map_pow, aeval_X]; exact (laurentMonomial _).isUnit.pow _

private noncomputable def toLaurent (c : Kˣ) (powers : J → ℕ) :
    SplitComplement c powers →ₐ[K] LaurentAlgebra K J :=
  IsLocalization.Away.liftAlgHom (splitNorm c powers) (laurentEvaluation_norm_unit c powers)

private theorem toLaurent_X (c : Kˣ) (powers : J → ℕ) (i : J) :
    toLaurent c powers (algebraMap (MvPolynomial J K) (SplitComplement c powers) (X i)) =
      (↑(laurentMonomial (K := K) (basisExponent i)) : LaurentAlgebra K J) := by
  have h := IsLocalization.Away.lift_eq (S := SplitComplement c powers)
    (splitNorm c powers) (laurentEvaluation_norm_unit c powers) (X i)
  exact h.trans (aeval_X _ _)

omit [DecidableEq J] in
private theorem splitVariable_isUnit (c : Kˣ) (powers : J → ℕ)
    (hp : ∀ i, 0 < powers i) (i : J) :
    IsUnit (algebraMap (MvPolynomial J K) (SplitComplement c powers) (X i)) := by
  apply isUnit_of_dvd_unit _ (IsLocalization.Away.algebraMap_isUnit (splitNorm c powers))
  apply map_dvd
  unfold splitNorm
  exact (dvd_pow_self _ (hp i).ne').trans ((Finset.dvd_prod_of_mem (fun j ↦ X j ^ powers j)
    (Finset.mem_univ i)).trans (dvd_mul_left _ _))

private noncomputable def splitVariable (c : Kˣ) (powers : J → ℕ)
    (hp : ∀ i, 0 < powers i) (i : J) : (SplitComplement c powers)ˣ :=
  (splitVariable_isUnit c powers hp i).unit

omit [DecidableEq J] in
private theorem splitVariable_val (c : Kˣ) (powers : J → ℕ) (hp : ∀ i, 0 < powers i) (i : J) :
    (splitVariable c powers hp i : SplitComplement c powers) =
      algebraMap (MvPolynomial J K) (SplitComplement c powers) (X i) :=
  IsUnit.unit_spec _

private noncomputable def fromLaurentMonomials (c : Kˣ) (powers : J → ℕ) (hp : ∀ i, 0 < powers i) :
    Multiplicative (J → ℤ) →* SplitComplement c powers where
  toFun v := (↑(∏ i, splitVariable c powers hp i ^ v.toAdd i) : SplitComplement c powers)
  map_one' := by simp
  map_mul' v w := by
    change (↑(∏ i, splitVariable c powers hp i ^ (v.toAdd i + w.toAdd i)) :
      SplitComplement c powers) = _
    simp only [zpow_add, Finset.prod_mul_distrib, Units.val_mul]

private noncomputable def fromLaurent (c : Kˣ) (powers : J → ℕ) (hp : ∀ i, 0 < powers i) :
    LaurentAlgebra K J →ₐ[K] SplitComplement c powers :=
  AddMonoidAlgebra.lift K _ (J → ℤ) (fromLaurentMonomials c powers hp)

private theorem fromLaurent_variable (c : Kˣ) (powers : J → ℕ) (hp : ∀ i, 0 < powers i) (i : J) :
    fromLaurent c powers hp (↑(laurentMonomial (K := K) (basisExponent i)) : LaurentAlgebra K J) =
      (splitVariable c powers hp i : SplitComplement c powers) := by
  simp only [fromLaurent, laurentMonomial, AddMonoidAlgebra.lift_single, one_smul]
  change (↑(∏ t, splitVariable c powers hp t ^ basisExponent i t) : SplitComplement c powers) = _
  rw [Finset.prod_eq_single i]
  · simp [basisExponent]
  · intro j _ hj
    simp [basisExponent, Pi.single_eq_of_ne hj]
  · simp

private theorem sum_basisExponent (v : J → ℤ) :
    ∑ i, v i • basisExponent i = v := by
  ext j
  simp [basisExponent, Pi.single_apply]

private theorem toLaurent_fromLaurent (c : Kˣ) (powers : J → ℕ) (hp : ∀ i, 0 < powers i) :
    (toLaurent c powers).comp (fromLaurent c powers hp) = AlgHom.id K (LaurentAlgebra K J) := by
  apply (AddMonoidAlgebra.lift K (LaurentAlgebra K J) (J → ℤ)).symm.injective
  apply MonoidHom.ext
  intro v
  rw [AddMonoidAlgebra.lift_symm_apply, AddMonoidAlgebra.lift_symm_apply]
  change toLaurent c powers (fromLaurent c powers hp (AddMonoidAlgebra.single v.toAdd 1)) = _
  rw [fromLaurent, AddMonoidAlgebra.lift_single, one_smul]
  have hunit (i : J) :
      Units.map (toLaurent c powers).toMonoidHom (splitVariable c powers hp i) =
        laurentMonomial (basisExponent i) := by
    apply Units.ext
    change toLaurent c powers (splitVariable c powers hp i : SplitComplement c powers) = _
    rw [splitVariable_val, toLaurent_X]
  change (↑(Units.map (toLaurent c powers).toMonoidHom
    (∏ i, splitVariable c powers hp i ^ v.toAdd i)) : LaurentAlgebra K J) =
      (↑(laurentMonomial (K := K) v.toAdd) : LaurentAlgebra K J)
  simp_rw [map_prod, map_zpow, hunit]
  congr 1
  change ∏ i, laurentMonomialHom (Multiplicative.ofAdd (basisExponent i)) ^ v.toAdd i =
    laurentMonomialHom v
  simp_rw [← map_zpow, ← ofAdd_zsmul]
  rw [← map_prod, ← ofAdd_sum, sum_basisExponent]
  rfl

private theorem fromLaurent_toLaurent (c : Kˣ) (powers : J → ℕ) (hp : ∀ i, 0 < powers i) :
    (fromLaurent c powers hp).comp (toLaurent c powers) = AlgHom.id K (SplitComplement c powers) := by
  apply IsLocalization.algHom_ext (Submonoid.powers (splitNorm c powers))
  apply MvPolynomial.algHom_ext
  intro i
  change fromLaurent c powers hp (toLaurent c powers (algebraMap _ _ (X i))) = algebraMap _ _ (X i)
  rw [toLaurent_X, fromLaurent_variable, splitVariable_val]

/-- An isomorphism carrying the localized variables to their Laurent monomials. -/
noncomputable def splitComplementEquivLaurent (c : Kˣ) (powers : J → ℕ) (hp : ∀ i, 0 < powers i) :
    SplitComplement c powers ≃ₐ[K] LaurentAlgebra K J :=
  AlgEquiv.ofAlgHom (toLaurent c powers) (fromLaurent c powers hp)
    (toLaurent_fromLaurent c powers hp) (fromLaurent_toLaurent c powers hp)

private noncomputable def scaleVariables (χ : J → Kˣ) :
    MvPolynomial J K →ₐ[K] MvPolynomial J K :=
  aeval fun i ↦ C (χ i : K) * X i

omit [DecidableEq J] in
private theorem scaleVariables_norm (c : Kˣ) (powers : J → ℕ) (χ : J → Kˣ) :
    scaleVariables χ (splitNorm c powers) =
      C (↑(∏ i, χ i ^ powers i) : K) * splitNorm c powers := by
  rw [splitNorm, map_mul, map_prod]
  have hc := map_prod (Units.coeHom K) (fun i ↦ χ i ^ powers i) Finset.univ
  change (↑(∏ i, χ i ^ powers i) : K) = ∏ i, (↑(χ i ^ powers i) : K) at hc
  rw [hc, map_prod]
  simp only [scaleVariables, aeval_C, algebraMap_eq, map_pow, aeval_X, mul_pow,
    Finset.prod_mul_distrib, Units.val_pow_eq_pow_val]
  ring

omit [DecidableEq J] in
private theorem scaleVariables_norm_unit (c : Kˣ) (powers : J → ℕ) (χ : J → Kˣ) :
    IsUnit (algebraMap (MvPolynomial J K) (SplitComplement c powers)
      (scaleVariables χ (splitNorm c powers))) := by
  rw [scaleVariables_norm, map_mul]
  refine IsUnit.mul ?_ (IsLocalization.Away.algebraMap_isUnit _)
  change IsUnit (algebraMap K (SplitComplement c powers) (↑(∏ i, χ i ^ powers i) : K))
  exact (∏ i, χ i ^ powers i).isUnit.map _

/-- The diagonal torus pullback on the split coordinate chart. -/
noncomputable def splitDiagonalPullback (c : Kˣ) (powers : J → ℕ) (χ : J → Kˣ) :
    SplitComplement c powers →ₐ[K] SplitComplement c powers :=
  IsLocalization.Away.liftAlgHom (splitNorm c powers)
    (f := (IsScalarTower.toAlgHom K (MvPolynomial J K) (SplitComplement c powers)).comp
      (scaleVariables χ))
    (scaleVariables_norm_unit c powers χ)

omit [DecidableEq J] in
theorem splitDiagonalPullback_X (c : Kˣ) (powers : J → ℕ) (χ : J → Kˣ) (i : J) :
    splitDiagonalPullback c powers χ (algebraMap (MvPolynomial J K) _ (X i)) =
      algebraMap K _ (χ i : K) * algebraMap (MvPolynomial J K) _ (X i) := by
  have he := IsLocalization.Away.lift_eq (S := SplitComplement c powers) (splitNorm c powers)
    (g := (algebraMap (MvPolynomial J K) (SplitComplement c powers)).comp (scaleVariables χ).toRingHom)
    (scaleVariables_norm_unit c powers χ) (X i)
  change splitDiagonalPullback c powers χ (algebraMap (MvPolynomial J K) _ (X i)) =
    algebraMap (MvPolynomial J K) (SplitComplement c powers) (scaleVariables χ (X i)) at he
  rw [he]
  rw [scaleVariables, aeval_X, map_mul]
  rfl

/-- On Laurent monomials the geometric pullback has the product
character appearing in the manuscript. -/
theorem splitDiagonalPullback_character (c : Kˣ) (powers : J → ℕ)
    (hp : ∀ i, 0 < powers i) (χ : J → Kˣ) (g : J → ℤ) :
    splitDiagonalPullback c powers χ ((splitComplementEquivLaurent c powers hp).symm
      (AddMonoidAlgebra.single g 1)) =
        ((↑(∏ i, χ i ^ g i) : K)) •
          (splitComplementEquivLaurent c powers hp).symm (AddMonoidAlgebra.single g 1) := by
  change splitDiagonalPullback c powers χ (fromLaurent c powers hp (AddMonoidAlgebra.single g 1)) =
    _ • fromLaurent c powers hp (AddMonoidAlgebra.single g 1)
  simp only [fromLaurent, AddMonoidAlgebra.lift_single, one_smul]
  let B := SplitComplement c powers
  let h := splitDiagonalPullback c powers χ
  let u := splitVariable c powers hp
  let φ := Units.map (algebraMap K B).toMonoidHom
  have hu (i : J) : Units.map h.toMonoidHom (u i) = φ (χ i) * u i := by
    apply Units.ext
    change h (splitVariable c powers hp i : B) =
      algebraMap K B (χ i : K) * (splitVariable c powers hp i : B)
    rw [splitVariable_val, splitDiagonalPullback_X]
  change h (↑(∏ i, u i ^ g i) : B) = (↑(∏ i, χ i ^ g i) : K) • (↑(∏ i, u i ^ g i) : B)
  rw [Algebra.smul_def]
  change (↑(Units.map h.toMonoidHom (∏ i, u i ^ g i)) : B) =
    (↑((φ (∏ i, χ i ^ g i)) * (∏ i, u i ^ g i)) : B)
  simp_rw [map_prod, map_zpow, hu, mul_zpow, Finset.prod_mul_distrib]

end MathieuQuadratic
