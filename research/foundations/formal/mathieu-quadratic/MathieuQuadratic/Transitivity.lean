import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.Tactic

/-!
# Geometric transitivity on the quadratic complement

Reflections yield determinant-one isometries between any two vectors of
the same nonzero norm in dimension at least two. Over an algebraically
closed field, the scalar action moves between nonzero norm levels.
This is the transitivity step in Section 2 of the draft.
-/

open scoped BigOperators Matrix

namespace MathieuQuadratic

variable {K J : Type*} [Field K] [Fintype J]

/-- The diagonal symmetric pairing, with a unit at every coordinate. -/
noncomputable def diagonalPair (d : J → Kˣ) (x y : J → K) : K :=
  ∑ i, (d i : K) * x i * y i

/-- The value of the diagonal quadratic form. -/
noncomputable def diagonalNorm (d : J → Kˣ) (x : J → K) : K :=
  diagonalPair d x x

theorem diagonalPair_comm (d : J → Kˣ) (x y : J → K) :
    diagonalPair d x y = diagonalPair d y x := by
  unfold diagonalPair
  exact Finset.sum_congr rfl fun i _ ↦ by ring

theorem diagonalPair_add_right (d : J → Kˣ) (x y z : J → K) :
    diagonalPair d x (y + z) = diagonalPair d x y + diagonalPair d x z := by
  simp [diagonalPair, mul_add, Finset.sum_add_distrib]

theorem diagonalPair_sub_right (d : J → Kˣ) (x y z : J → K) :
    diagonalPair d x (y - z) = diagonalPair d x y - diagonalPair d x z := by
  simp [diagonalPair, mul_sub, Finset.sum_sub_distrib]

theorem diagonalPair_smul_right (d : J → Kˣ) (x y : J → K) (c : K) :
    diagonalPair d x (c • y) = c * diagonalPair d x y := by
  simp only [diagonalPair, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ ↦ by ring

theorem diagonalPair_single_right [DecidableEq J] (d : J → Kˣ) (x : J → K) (i : J) :
    diagonalPair d x (Pi.single i 1) = (d i : K) * x i := by
  rw [diagonalPair, Finset.sum_eq_single i]
  · simp
  · intro j _ hj
    simp [Pi.single_eq_of_ne hj]
  · simp

theorem diagonalNorm_smul (d : J → Kˣ) (x : J → K) (c : K) :
    diagonalNorm d (c • x) = c ^ 2 * diagonalNorm d x := by
  unfold diagonalNorm
  rw [diagonalPair_smul_right, diagonalPair_comm d (c • x) x, diagonalPair_smul_right]
  ring

theorem diagonalNorm_add (d : J → Kˣ) (x y : J → K) :
    diagonalNorm d (x + y) =
      diagonalNorm d x + 2 * diagonalPair d x y + diagonalNorm d y := by
  unfold diagonalNorm
  rw [diagonalPair_add_right, diagonalPair_comm d (x + y) x,
    diagonalPair_comm d (x + y) y, diagonalPair_add_right, diagonalPair_add_right,
    diagonalPair_comm d y x]
  ring

theorem diagonalNorm_sub (d : J → Kˣ) (x y : J → K) :
    diagonalNorm d (x - y) =
      diagonalNorm d x - 2 * diagonalPair d x y + diagonalNorm d y := by
  unfold diagonalNorm
  rw [diagonalPair_sub_right, diagonalPair_comm d (x - y) x,
    diagonalPair_comm d (x - y) y, diagonalPair_sub_right, diagonalPair_sub_right,
    diagonalPair_comm d y x]
  ring

private noncomputable def reflection [DecidableEq J] (d : J → Kˣ) (w : J → K) :
    Matrix J J K :=
  1 + Matrix.replicateCol Unit w *
    Matrix.replicateRow Unit (fun i ↦ (-2 / diagonalNorm d w) * ((d i : K) * w i))

private theorem reflection_dot (d : J → Kˣ) (w v : J → K) :
    (fun i ↦ (-2 / diagonalNorm d w) * ((d i : K) * w i)) ⬝ᵥ v =
      (-2 / diagonalNorm d w) * diagonalPair d w v := by
  simp only [dotProduct, diagonalPair, Finset.mul_sum, mul_assoc]

private theorem reflection_action [DecidableEq J] (d : J → Kˣ) (w v : J → K) :
    reflection d w *ᵥ v =
      v + ((-2 / diagonalNorm d w) * diagonalPair d w v) • w := by
  rw [reflection, Matrix.add_mulVec, Matrix.one_mulVec, ← Matrix.mulVec_mulVec,
    Matrix.replicateRow_mulVec_eq_const, reflection_dot]
  ext i
  simp [Matrix.mulVec, dotProduct, mul_comm]

private theorem reflection_det [DecidableEq J] (d : J → Kˣ) (w : J → K)
    (hw : diagonalNorm d w ≠ 0) : (reflection d w).det = -1 := by
  rw [reflection, Matrix.det_one_add_replicateCol_mul_replicateRow, reflection_dot]
  change 1 + -2 / diagonalNorm d w * diagonalNorm d w = -1
  rw [div_mul_cancel₀ _ hw]
  ring

private theorem reflection_preserves [DecidableEq J] (d : J → Kˣ) (w v : J → K)
    (hw : diagonalNorm d w ≠ 0) :
    diagonalNorm d (reflection d w *ᵥ v) = diagonalNorm d v := by
  rw [reflection_action, diagonalNorm_add, diagonalPair_smul_right,
    diagonalNorm_smul, diagonalPair_comm d v w]
  field_simp
  ring

private theorem reflection_self [DecidableEq J] (d : J → Kˣ) (y : J → K)
    (hy : diagonalNorm d y ≠ 0) : reflection d y *ᵥ y = -y := by
  rw [reflection_action]
  change y + ((-2 / diagonalNorm d y) * diagonalNorm d y) • y = -y
  rw [div_mul_cancel₀ _ hy]
  ext i
  simp
  ring

variable [CharZero K]

private theorem exists_nonisotropic_perpendicular [DecidableEq J] (d : J → Kˣ) (y : J → K)
    (hr : 1 < Fintype.card J) (hy : diagonalNorm d y ≠ 0) :
    ∃ w : J → K, diagonalPair d y w = 0 ∧ diagonalNorm d w ≠ 0 := by
  have hn : ∃ i, y i ≠ 0 := by
    by_contra! hw
    exact hy (by simp [diagonalNorm, diagonalPair, hw])
  obtain ⟨i, hi⟩ := hn
  obtain ⟨j, hj⟩ := Fintype.exists_ne_of_one_lt_card hr i
  let u := Pi.single j (1 : K) -
    (((d j : K) * y j) / ((d i : K) * y i)) • Pi.single i (1 : K)
  have hu : diagonalPair d y u = 0 := by
    rw [show u = _ from rfl, diagonalPair_sub_right, diagonalPair_smul_right,
      diagonalPair_single_right, diagonalPair_single_right,
      div_mul_cancel₀ _ (mul_ne_zero (Units.ne_zero _) hi), sub_self]
  by_cases hu0 : diagonalNorm d u = 0
  · let v := Pi.single j (1 : K) -
      ((diagonalPair d y (Pi.single j (1 : K))) / diagonalNorm d y) • y
    have hv : diagonalPair d y v = 0 := by
      rw [show v = _ from rfl, diagonalPair_sub_right, diagonalPair_smul_right]
      change diagonalPair d y (Pi.single j 1) -
        diagonalPair d y (Pi.single j 1) / diagonalNorm d y * diagonalNorm d y = 0
      rw [div_mul_cancel₀ _ hy, sub_self]
    by_cases hv0 : diagonalNorm d v = 0
    · refine ⟨u + v, by rw [diagonalPair_add_right, hu, hv, add_zero], ?_⟩
      rw [diagonalNorm_add, hu0, hv0, zero_add, add_zero]
      have huv : diagonalPair d u v = (d j : K) := by
        rw [show v = _ from rfl, diagonalPair_sub_right, diagonalPair_smul_right,
          diagonalPair_comm d u y, hu, mul_zero, sub_zero, diagonalPair_single_right]
        simp [u, Pi.single_eq_of_ne hj]
      rw [huv]
      exact mul_ne_zero two_ne_zero (Units.ne_zero _)
    · exact ⟨v, hv, hv0⟩
  · exact ⟨u, hu, hu0⟩

/-- The special orthogonal action is transitive on each nonzero norm
level for every diagonal nondegenerate form of dimension at least two. -/
theorem exists_specialOrthogonal_same_norm [DecidableEq J] (d : J → Kˣ) (x y : J → K)
    (hr : 1 < Fintype.card J) (hxy : diagonalNorm d x = diagonalNorm d y)
    (hy : diagonalNorm d y ≠ 0) :
    ∃ T : Matrix J J K, T.det = 1 ∧
      (∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) ∧ T *ᵥ x = y := by
  by_cases hplus : diagonalNorm d (x + y) = 0
  · have hminus : diagonalNorm d (x - y) ≠ 0 := by
      intro hz
      rw [diagonalNorm_add, hxy] at hplus
      rw [diagonalNorm_sub, hxy] at hz
      exact hy (by linear_combination (hplus + hz) / 4)
    have hx' : reflection d (x - y) *ᵥ x = y := by
      rw [reflection_action, diagonalPair_comm d (x - y) x, diagonalPair_sub_right]
      have he : diagonalNorm d (x - y) = 2 * (diagonalNorm d x - diagonalPair d x y) := by
        rw [diagonalNorm_sub, ← hxy]
        ring
      change x + ((-2 / diagonalNorm d (x - y)) * (diagonalNorm d x - diagonalPair d x y)) •
        (x - y) = y
      have hc : (-2 / diagonalNorm d (x - y)) * (diagonalNorm d x - diagonalPair d x y) = -1 := by
        field_simp
        rw [he]
      rw [hc]
      ext k
      simp
    obtain ⟨w, hperp, hw⟩ := exists_nonisotropic_perpendicular d y hr hy
    refine ⟨reflection d w * reflection d (x - y), ?_, ?_, ?_⟩
    · rw [Matrix.det_mul, reflection_det d w hw, reflection_det d _ hminus]
      ring
    · intro v
      rw [← Matrix.mulVec_mulVec, reflection_preserves d w _ hw, reflection_preserves d _ _ hminus]
    · rw [← Matrix.mulVec_mulVec, hx', reflection_action, diagonalPair_comm d w y,
        hperp, mul_zero, zero_smul, add_zero]
  · have hx' : reflection d (x + y) *ᵥ x = -y := by
      rw [reflection_action, diagonalPair_comm d (x + y) x, diagonalPair_add_right]
      have he : diagonalNorm d (x + y) = 2 * (diagonalNorm d x + diagonalPair d x y) := by
        rw [diagonalNorm_add, ← hxy]
        ring
      change x + ((-2 / diagonalNorm d (x + y)) * (diagonalNorm d x + diagonalPair d x y)) •
        (x + y) = -y
      have hc : (-2 / diagonalNorm d (x + y)) * (diagonalNorm d x + diagonalPair d x y) = -1 := by
        field_simp
        rw [he]
      rw [hc]
      ext k
      simp
    refine ⟨reflection d y * reflection d (x + y), ?_, ?_, ?_⟩
    · rw [Matrix.det_mul, reflection_det d y hy, reflection_det d _ hplus]
      ring
    · intro v
      rw [← Matrix.mulVec_mulVec, reflection_preserves d y _ hy, reflection_preserves d _ _ hplus]
    · rw [← Matrix.mulVec_mulVec, hx', Matrix.mulVec_neg, reflection_self d y hy, neg_neg]

/-- Over an algebraically closed field, the special orthogonal and scalar
actions together are transitive on the nonzero-norm locus. -/
theorem exists_isometry_and_scalar [DecidableEq J] [IsAlgClosed K]
    (d : J → Kˣ) (x y : J → K) (hr : 1 < Fintype.card J)
    (hx : diagonalNorm d x ≠ 0) (hy : diagonalNorm d y ≠ 0) :
    ∃ (T : Matrix J J K) (c : Kˣ), T.det = 1 ∧
      (∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) ∧ c • (T *ᵥ x) = y := by
  obtain ⟨r, hr'⟩ := IsAlgClosed.exists_pow_nat_eq (diagonalNorm d y / diagonalNorm d x)
    (n := 2) (by decide)
  have hr0 : r ≠ 0 := by
    intro h
    exact (div_ne_zero hy hx) (by simpa [h] using hr'.symm)
  have heq : diagonalNorm d (r • x) = diagonalNorm d y := by
    rw [diagonalNorm_smul, hr', div_mul_cancel₀ _ hx]
  obtain ⟨T, hd, hf, he⟩ := exists_specialOrthogonal_same_norm d (r • x) y hr heq hy
  exact ⟨T, Units.mk0 r hr0, hd, hf, by rw [Matrix.mulVec_smul] at he; exact he⟩

end MathieuQuadratic
