import MathieuQuadratic.Complement

/-!
# The rank-two Lie recurrence

A rotation derivation annihilates the norm. Its sum on the polynomial
numerators used in (8) gives the height recurrence after localization.
-/

open scoped BigOperators

namespace MathieuQuadratic

open MvPolynomial

variable {K : Type*} [Field K] {M : ℕ}

private noncomputable def x (i : Fin 2) : QuadPolynomial K M := X (Sum.inl i)
private noncomputable def z (i : Fin M) : QuadPolynomial K M := X (Sum.inr i)

private noncomputable def planeTerm (a : Fin 2 → K) (i : Fin 2) : QuadPolynomial K M :=
  C (a i) * x i ^ 2

private noncomputable def heightTerm (b : Fin M → K) (i : Fin M) : QuadPolynomial K M :=
  C (b i) * z i ^ 2

/-- The infinitesimal rotation between one plane and one height coordinate. -/
noncomputable def rotation (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    Derivation K (QuadPolynomial K M) (QuadPolynomial K M) :=
  (C (b j) * z j) • pderiv (R := K) (Sum.inl i : Fin 2 ⊕ Fin M) -
    (C (a i) * x (M := M) i) • pderiv (R := K) (Sum.inr j : Fin 2 ⊕ Fin M)

private theorem rotation_apply (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M)
    (f : QuadPolynomial K M) :
    rotation a b i j f =
      (C (b j) * z j) * pderiv (Sum.inl i) f - (C (a i) * x i) * pderiv (Sum.inr j) f := by
  simp [rotation, smul_eq_mul]

private theorem height_pderiv_plane (b : Fin M → K) (i : Fin 2) :
    pderiv (Sum.inl i) (heightNorm b) = 0 := by
  simp [heightNorm]

private theorem plane_pderiv_height (a : Fin 2 → K) (j : Fin M) :
    pderiv (Sum.inr j) (planeNorm (M := M) a) = 0 := by
  simp [planeNorm]

private theorem height_pderiv_height (b : Fin M → K) (j : Fin M) :
    pderiv (Sum.inr j) (heightNorm b) = C (b j) * (2 * z j) := by
  rw [heightNorm, map_sum, Finset.sum_eq_single j]
  · simp [z]
  · intro t _ ht
    simp [pderiv_X_of_ne (show Sum.inr t ≠ Sum.inr j from by simpa)]
  · simp

private theorem plane_pderiv_plane (a : Fin 2 → K) (i : Fin 2) :
    pderiv (Sum.inl i) (planeNorm (M := M) a) = C (a i) * (2 * x i) := by
  rw [planeNorm, map_sum, Finset.sum_eq_single i]
  · simp [x]
  · intro t _ ht
    simp [pderiv_X_of_ne (show Sum.inl t ≠ Sum.inl i from by simpa)]
  · simp

theorem rotation_quadNorm (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    rotation a b i j (quadNorm a b) = 0 := by
  rw [rotation_apply, quadNorm, map_add, map_add, plane_pderiv_plane, height_pderiv_height,
    height_pderiv_plane, plane_pderiv_height]
  ring

private theorem rotation_height (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    rotation a b i j (heightNorm b) = -(2 * (C (a i) * x i) * (C (b j) * z j)) := by
  rw [rotation_apply, height_pderiv_plane, height_pderiv_height]
  ring

private theorem rotation_xz (a : Fin 2 → K) (b : Fin M → K) (i : Fin 2) (j : Fin M) :
    rotation a b i j (x i * z j) = heightTerm b j - planeTerm a i := by
  rw [rotation_apply, pderiv_mul, pderiv_mul]
  simp only [planeTerm, heightTerm, x, z, pderiv_X_self, pderiv_X_of_ne (Sum.inr_ne_inl),
    pderiv_X_of_ne (Sum.inl_ne_inr), one_mul, mul_zero, add_zero, zero_mul, zero_add]
  ring

private theorem rotation_numerator (a : Fin 2 → K) (b : Fin M → K)
    (i : Fin 2) (j : Fin M) (n : ℕ) :
    rotation a b i j (heightNorm b ^ n * (x i * z j)) =
      heightNorm b ^ n * heightTerm b j - heightNorm b ^ n * planeTerm a i -
        (2 * (n : QuadPolynomial K M) * heightNorm b ^ (n - 1)) *
          (planeTerm a i * heightTerm b j) := by
  rw [Derivation.leibniz, Derivation.leibniz_pow, rotation_xz, rotation_height]
  simp only [smul_eq_mul, nsmul_eq_mul, planeTerm, heightTerm]
  ring

theorem sum_rotation_numerators (a : Fin 2 → K) (b : Fin M → K) (n : ℕ) :
    ∑ i : Fin 2, ∑ j : Fin M, rotation a b i j (heightNorm b ^ n * (x i * z j)) =
      2 * heightNorm b ^ (n + 1) -
        ((M : QuadPolynomial K M) + 2 * (n : QuadPolynomial K M)) * planeNorm a * heightNorm b ^ n := by
  have hb : (∑ j, heightTerm b j) = heightNorm b := rfl
  have ha : (∑ i, planeTerm (M := M) a i) = planeNorm a := rfl
  simp_rw [rotation_numerator, Finset.sum_sub_distrib, ← Finset.mul_sum,
    hb, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  simp_rw [← Finset.mul_sum, ← Finset.sum_mul, ha]
  cases n with
  | zero => simp
  | succ n =>
      simp only [Nat.succ_sub_one, pow_succ, Nat.cast_succ]
      ring

section Trace

variable (a : Fin 2 → K) (b : Fin M → K)

/-- An invariant trace annihilates these rotation images with norm-power
denominators. The chosen rotations preserve q by rotation_quadNorm. -/
def AnnihilatesRotations (trace : NormComplement a b →ₗ[K] K) : Prop :=
  ∀ (i : Fin 2) (j : Fin M) (f : QuadPolynomial K M) (n : ℕ),
    trace (algebraMap (QuadPolynomial K M) (NormComplement a b) (rotation a b i j f) *
      (↑((normUnit a b)⁻¹) : NormComplement a b) ^ n) = 0

theorem height_recurrence_of_rotations (trace : NormComplement a b →ₗ[K] K)
    (h : AnnihilatesRotations a b trace) (n : ℕ) :
    (↑(M + 2 * (n + 1)) : K) * trace (height a b ^ (n + 1)) =
      (↑(M + 2 * n) : K) * trace (height a b ^ n) := by
  let A := NormComplement a b
  let f : QuadPolynomial K M →+* A := algebraMap _ _
  let v : A := ↑((normUnit a b)⁻¹)
  let S : A := height a b
  have hS : S = f (heightNorm b) * v := rfl
  have hc : f (planeNorm a) * v = 1 - S := by
    have hu : f (quadNorm a b) * v = 1 := by
      rw [← normUnit_val]
      exact Units.mul_inv _
    change f (planeNorm a + heightNorm b) * v = 1 at hu
    rw [map_add, add_mul, ← hS] at hu
    exact eq_sub_of_add_eq hu
  have hzero : trace (f (∑ i : Fin 2, ∑ j : Fin M,
      rotation a b i j (heightNorm b ^ n * (x i * z j))) * v ^ (n + 1)) = 0 := by
    simp only [map_sum, Finset.sum_mul]
    exact Finset.sum_eq_zero fun i _ ↦ Finset.sum_eq_zero fun j _ ↦ h i j _ _
  rw [sum_rotation_numerators] at hzero
  simp only [map_sub, map_mul, map_pow, map_add, map_natCast, map_ofNat, sub_mul] at hzero
  have halg :
      (2 : A) * f (heightNorm b) ^ (n + 1) * v ^ (n + 1) -
        ((M : A) + 2 * (n : A)) * f (planeNorm a) * f (heightNorm b) ^ n * v ^ (n + 1) =
      (2 : A) * S ^ (n + 1) -
        ((M : A) + 2 * (n : A)) * (1 - S) * S ^ n := by
    rw [← hc, hS, mul_pow, mul_pow, pow_succ v n]
    ring
  rw [← map_sub, halg] at hzero
  have hlin (c : K) (p : A) : trace (algebraMap K A c * p) = c * trace p := by
    rw [← Algebra.smul_def, map_smul, smul_eq_mul]
  have hfor :
      (2 : A) * S ^ (n + 1) - ((M : A) + 2 * (n : A)) * (1 - S) * S ^ n =
        algebraMap K A (M + 2 * ((n : K) + 1)) * S ^ (n + 1) -
          algebraMap K A (M + 2 * (n : K)) * S ^ n := by
    simp only [map_add, map_mul, map_natCast, map_one, map_ofNat, pow_succ]
    ring
  rw [hfor, map_sub, hlin, hlin, sub_eq_zero] at hzero
  simpa using hzero

end Trace

end MathieuQuadratic
