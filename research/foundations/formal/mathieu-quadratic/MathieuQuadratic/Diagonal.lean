import MathieuQuadratic.Complement
import Mathlib.LinearAlgebra.QuadraticForm.IsometryEquiv

/-!
# The orthogonal splitting used by the uniform construction
-/

open scoped BigOperators

namespace MathieuQuadratic

variable {K V : Type*} [Field K] [CharZero K] [AddCommGroup V] [Module K V]
variable [FiniteDimensional K V]

theorem exists_plane_height_diagonal (q : QuadraticForm K V)
    (hnondegenerate : (QuadraticMap.associated (R := K) q).SeparatingLeft)
    (M : ℕ) (hM : Module.finrank K V = M + 2) :
    ∃ (a : Fin 2 → Kˣ) (b : Fin M → Kˣ),
      q.Equivalent (QuadraticMap.weightedSumSquares K (Sum.elim a b)) := by
  letI : Invertible (2 : K) := invertibleOfNonzero two_ne_zero
  obtain ⟨w, ⟨iso⟩⟩ := q.equivalent_weightedSumSquares_units_of_nondegenerate' hnondegenerate
  let e : Fin 2 ⊕ Fin M ≃ Fin (Module.finrank K V) := Fintype.equivOfCardEq (by simp [hM, add_comm])
  let a : Fin 2 → Kˣ := fun i ↦ w (e (Sum.inl i))
  let b : Fin M → Kˣ := fun i ↦ w (e (Sum.inr i))
  have he (i : Fin 2 ⊕ Fin M) : Sum.elim a b i = w (e i) := by cases i <;> rfl
  let reindex : (QuadraticMap.weightedSumSquares K w).IsometryEquiv
      (QuadraticMap.weightedSumSquares K (Sum.elim a b)) :=
    { __ := LinearEquiv.piCongrLeft' K (fun _ : Fin (Module.finrank K V) ↦ K) e.symm
      map_app' := fun v ↦ by
        simp only [QuadraticMap.weightedSumSquares_apply, he]
        exact e.sum_comp fun i ↦ w i • (v i * v i) }
  exact ⟨a, b, ⟨iso.trans reindex⟩⟩

omit [CharZero K] in
theorem eval_quadNorm {M : ℕ} (a : Fin 2 → Kˣ) (b : Fin M → Kˣ) (v : Fin 2 ⊕ Fin M → K) :
    MvPolynomial.eval v (quadNorm (fun i ↦ (a i : K)) (fun i ↦ (b i : K))) =
      QuadraticMap.weightedSumSquares K (Sum.elim a b) v := by
  simp only [quadNorm, planeNorm, heightNorm, map_add, map_sum, map_mul,
    MvPolynomial.eval_C, MvPolynomial.eval_X, QuadraticMap.weightedSumSquares_apply, Fintype.sum_sum_type,
    Sum.elim_inl, Sum.elim_inr, Units.smul_def, smul_eq_mul, pow_two]

end MathieuQuadratic
