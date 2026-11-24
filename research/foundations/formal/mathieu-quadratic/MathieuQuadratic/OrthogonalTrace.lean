import MathieuQuadratic.DualTrace
import MathieuQuadratic.Counterexample
import MathieuQuadratic.Infinitesimal
import Mathlib.LinearAlgebra.Matrix.Charpoly.Coeff
import Mathlib.Data.Matrix.Basis

/-!
# Infinitesimal special orthogonal points on the norm complement

The rotations of the draft are realized by determinant-one matrices
over the dual numbers. Invariance of the extended functional on such
points gives the precise trace identity needed by the mixed moments.
-/

open scoped BigOperators DualNumber Matrix
open TrivSqZeroExt

namespace MathieuQuadratic

variable {K : Type*} [Field K] {M : ℕ}

local notation "I" => Fin 2 ⊕ Fin M

private theorem determinant_one_firstOrder {R J : Type*} [CommRing R] [Fintype J]
    [DecidableEq J] (e : R) (N : Matrix J J R) (he : e ^ 2 = 0) (hn : N.trace = 0) :
    (1 + e • N).det = 1 := by
  rw [Matrix.det_one_add_smul, he, hn, zero_mul, mul_zero, add_zero, add_zero]

private theorem matrix_one_add_isometry {R J : Type*} [CommRing R] [Fintype J]
    [DecidableEq J] (e : R) (D N : Matrix J J R) (he : e * e = 0)
    (hn : Nᵀ * D + D * N = 0) :
    (1 + e • N)ᵀ * D * (1 + e • N) = D := by
  calc
    _ = D + e • (Nᵀ * D + D * N) + (e * e) • (Nᵀ * D * N) := by
      simp only [Matrix.transpose_add, Matrix.transpose_one, Matrix.transpose_smul,
        add_mul, mul_add, one_mul, mul_one, Matrix.smul_mul, Matrix.mul_smul, smul_smul,
        smul_add]
      abel
    _ = D := by rw [he, hn, smul_zero, zero_smul, add_zero, add_zero]

private def rotationMatrix (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    Matrix I I K[ε] :=
  Matrix.single (Sum.inl i) (Sum.inr j) (inl (b j)) -
    Matrix.single (Sum.inr j) (Sum.inl i) (inl (a i))

/-- The special orthogonal point realizing a plane-height rotation. -/
def rotationPoint (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    Matrix I I K[ε] :=
  1 + (ε : K[ε]) • rotationMatrix a b i j

theorem rotationPoint_det (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    (rotationPoint a b i j).det = 1 := by
  apply determinant_one_firstOrder _ _ DualNumber.eps_pow_two
  simp [rotationMatrix, Matrix.trace_single_eq_of_ne]

/-- The matrix of the quadratic form in the chosen diagonal coordinates. -/
def rotationDiagonal (a : Fin 2 → K) (b : Fin M → K) : Matrix I I K[ε] :=
  Matrix.diagonal fun k ↦ inl (Sum.elim a b k)

theorem rotationPoint_isometry (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    (rotationPoint a b i j)ᵀ * rotationDiagonal a b * rotationPoint a b i j =
      rotationDiagonal a b := by
  apply matrix_one_add_isometry _ _ _ DualNumber.eps_mul_eps
  simp only [rotationMatrix, rotationDiagonal, Matrix.transpose_sub, Matrix.transpose_single,
    sub_mul, mul_sub]
  apply Matrix.ext
  intro k l
  simp only [Matrix.add_apply, Matrix.sub_apply, Matrix.zero_apply, Matrix.mul_diagonal,
    Matrix.diagonal_mul]
  by_cases hk : k = Sum.inl i <;> by_cases hl : l = Sum.inr j <;>
    by_cases hk' : k = Sum.inr j <;> by_cases hl' : l = Sum.inl i <;>
    simp_all [Matrix.single, eq_comm, mul_comm]

variable (a : Fin 2 → K) (b : Fin M → K)

local notation "R" => QuadPolynomial K M
local notation "A" => NormComplement a b

private noncomputable abbrev mapDualCoefficients : K[ε] →+* A[ε] := dualCoefficientMap

private theorem dual_coeff_inl (x : K) :
    mapDualCoefficients a b (inl x) = inl (algebraMap K A x) := by
  ext <;> simp [mapDualCoefficients]

private theorem dual_coeff_eps : mapDualCoefficients a b (ε : K[ε]) = ε := by
  ext <;> simp [mapDualCoefficients]

private noncomputable def coordinateVector (k : I) : A[ε] :=
  inl (algebraMap R A (MvPolynomial.X k))

/-- A pullback is induced by the displayed dual-number matrix if it
agrees with the matrix action on all coordinate generators. -/
def IsMatrixPullback (T : Matrix I I K[ε]) (f : A →ₐ[K] A[ε]) : Prop :=
  ∀ k, f (algebraMap R A (MvPolynomial.X k)) =
    (T.map (mapDualCoefficients a b) *ᵥ coordinateVector a b) k

theorem firstOrderRotation_isMatrixPullback (i : Fin 2) (j : Fin M) :
    IsMatrixPullback a b (rotationPoint a b i j) (firstOrderRotation a b i j) := by
  intro k
  rw [firstOrderRotation_X, rotation_X, map_sub]
  let d := mapDualCoefficients a b
  change _ = ((rotationPoint a b i j).map d *ᵥ coordinateVector a b) k
  simp only [rotationPoint, rotationMatrix, Matrix.map_add d (map_add d),
    Matrix.map_one d (map_zero d) (map_one d), Matrix.map_smul' d _ _ (map_mul d),
    Matrix.map_sub d (map_sub d), Matrix.map_single _ _ _ d]
  change _ = ((1 + mapDualCoefficients a b (ε : K[ε]) •
    (Matrix.single (Sum.inl i) (Sum.inr j) (mapDualCoefficients a b (inl (b j))) -
      Matrix.single (Sum.inr j) (Sum.inl i) (mapDualCoefficients a b (inl (a i))))) *ᵥ
      coordinateVector a b) k
  simp only [dual_coeff_inl, dual_coeff_eps,
    Matrix.add_mulVec, Matrix.one_mulVec, Matrix.smul_mulVec, Matrix.sub_mulVec,
    Matrix.single_mulVec_eq, coordinateVector, Pi.add_apply, Pi.sub_apply, Pi.smul_apply]
  by_cases hk₁ : k = Sum.inl i <;> by_cases hk₂ : k = Sum.inr j <;>
    apply TrivSqZeroExt.ext <;>
    simp_all [smul_eq_mul, eq_comm, IsScalarTower.algebraMap_apply K R A]

variable [CharZero K]

/-- The part of normalized special orthogonal and scalar invariance used by
the draft. The dual-number condition evaluates the first-order coefficient
of the scalar-extended trace on each special orthogonal pullback. -/
structure GeometricTrace (trace : A →ₗ[K] K) : Prop where
  normalized : trace 1 = 1
  scalarTwo : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f
  orthogonal : ∀ (T : Matrix I I K[ε]) (f : A →ₐ[K] A[ε]),
    T.det = 1 → Tᵀ * rotationDiagonal a b * T = rotationDiagonal a b →
    IsMatrixPullback a b T f → ∀ g, trace (f g).snd = 0

theorem GeometricTrace.complementTrace {trace : A →ₗ[K] K} (h : GeometricTrace a b trace) :
    ComplementTrace a b trace := by
  refine ⟨h.normalized, h.scalarTwo, annihilatesRotations_of_firstOrder_invariance a b trace ?_⟩
  intro i j f
  exact h.orthogonal _ _ (rotationPoint_det a b i j) (rotationPoint_isometry a b i j)
    (firstOrderRotation_isMatrixPullback a b i j) f

def complementDiagonalWeights (ha : ∀ i, a i ≠ 0) (hb : ∀ j, b j ≠ 0) : I → Kˣ :=
  Sum.elim (fun i ↦ Units.mk0 (a i) (ha i)) (fun j ↦ Units.mk0 (b j) (hb j))

omit [CharZero K] in
theorem quadNorm_eq_diagonal (ha : ∀ i, a i ≠ 0) (hb : ∀ j, b j ≠ 0) :
    quadNorm a b = diagonalNormPolynomial (complementDiagonalWeights a b ha hb) := by
  simp only [quadNorm, planeNorm, heightNorm, diagonalNormPolynomial, Fintype.sum_sum_type,
    complementDiagonalWeights, Sum.elim_inl, Sum.elim_inr, Units.val_mk0]

/-- The algebraic trace on the complement in the plane-height coordinates. -/
noncomputable def complementAlgebraicTrace (ha : ∀ i, a i ≠ 0) (hb : ∀ j, b j ≠ 0) : A →ₗ[K] K :=
  algebraicTraceAt (complementDiagonalWeights a b ha hb) (quadNorm a b)
    (quadNorm_eq_diagonal a b ha hb)

/-- The algebraic trace satisfies every geometric trace condition used in
the draft's uniform witness. -/
theorem complementAlgebraicTrace_isGeometric (ha : ∀ i, a i ≠ 0) (hb : ∀ j, b j ≠ 0) :
    GeometricTrace a b (complementAlgebraicTrace a b ha hb) where
  normalized := algebraicTraceAt_one _ _ _
  scalarTwo f := algebraicTraceAt_scalar _ _ _ _
    (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero))
    (scalarPullback_X a b (Units.mk0 (2 : K) two_ne_zero)) f
  orthogonal T ψ hdet hT hψ f := by
    apply algebraicTraceAt_firstOrder_zero _ _ _ _ _ hdet _ _ f
    · have he : (fun j ↦ inl (complementDiagonalWeights a b ha hb j : K)) =
          fun j ↦ (inl (Sum.elim a b j) : K[ε]) := by
        funext j
        cases j <;> rfl
      rw [he]
      exact hT
    · exact hψ

end MathieuQuadratic
