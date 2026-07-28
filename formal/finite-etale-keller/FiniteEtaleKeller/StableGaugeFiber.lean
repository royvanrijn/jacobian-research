/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeJacobian
import FiniteEtaleKeller.GeneralGaugeRealization

/-!
# Fiber-invisible stable gauge deformations

The common power shifts and the cubic lifts change the polynomial map while
leaving its selected fiber over first target coordinate one unchanged.  This
file proves that assertion on literal points and transports the existing
represented-fiber theorem to both stable families.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K] [CharZero K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

private theorem shifted_power_product
    (t q x : A) (r s m : ℕ) (h : t * q = 1) :
    t ^ (m + r) * x * q ^ (s + m) = t ^ r * x * q ^ s := by
  rw [Nat.add_comm m r, pow_add, pow_add]
  calc
    t ^ r * t ^ m * x * (q ^ s * q ^ m) =
        (t * q) ^ m * (t ^ r * x * q ^ s) := by
          rw [mul_pow]
          ring
    _ = t ^ r * x * q ^ s := by rw [h, one_pow, one_mul]

private theorem shifted_power_product_x
    (t q x : A) (k m : ℕ) (h : t * q = 1) :
    t ^ m * x ^ k * q ^ (k + m) = (x * q) ^ k := by
  rw [pow_add]
  calc
    t ^ m * x ^ k * (q ^ k * q ^ m) =
        (t * q) ^ m * (x ^ k * q ^ k) := by
          rw [mul_pow]
          ring
    _ = x ^ k * q ^ k := by rw [h, one_pow, one_mul]
    _ = (x * q) ^ k := (mul_pow x q k).symm

/-- A common power shift is invisible in the second coordinate on
`Π = 1`. -/
theorem eval₂_powerShiftedGaugeB_eq_generalGaugeB
    (G : K[X]) (m : ℕ) (point : Fin 3 → A)
    (hpi : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugePi G) = 1) :
    MvPolynomial.eval₂ (algebraMap K A) point
        (powerShiftedGaugeB G m) =
      MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeB G) := by
  have hterm (k : ℕ) :
      MvPolynomial.eval₂ (algebraMap K A) point
          (MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
            generalGaugeT ^ (m + 2) * MvPolynomial.X 0 ^ (k - 2) *
              generalGaugeQ G ^ (k + m)) =
        MvPolynomial.eval₂ (algebraMap K A) point
          (MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
            generalGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
              generalGaugeQ G ^ k) := by
    simp only [MvPolynomial.eval₂_mul, MvPolynomial.eval₂_pow,
      MvPolynomial.eval₂_X, MvPolynomial.eval₂_C]
    have hs := shifted_power_product
      (MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT)
      (MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeQ G))
      (point 0 ^ (k - 2)) 2 k m (by
        simpa [generalGaugePi] using hpi)
    calc
      (algebraMap K A) ((k : K) * (G.coeff k / G.coeff 1)) *
            MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT ^
              (m + 2) *
            point 0 ^ (k - 2) *
            MvPolynomial.eval₂ (algebraMap K A) point
              (generalGaugeQ G) ^ (k + m) =
          (algebraMap K A) ((k : K) * (G.coeff k / G.coeff 1)) *
            (MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT ^
                (m + 2) * point 0 ^ (k - 2) *
              MvPolynomial.eval₂ (algebraMap K A) point
                (generalGaugeQ G) ^ (k + m)) := by ring
      _ = (algebraMap K A) ((k : K) * (G.coeff k / G.coeff 1)) *
            (MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT ^ 2 *
              point 0 ^ (k - 2) *
              MvPolynomial.eval₂ (algebraMap K A) point
                (generalGaugeQ G) ^ k) := by rw [hs]
      _ = (algebraMap K A) ((k : K) * (G.coeff k / G.coeff 1)) *
            MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT ^ 2 *
            point 0 ^ (k - 2) *
            MvPolynomial.eval₂ (algebraMap K A) point
              (generalGaugeQ G) ^ k := by ring
  have hsum :
      MvPolynomial.eval₂ (algebraMap K A) point
          (∑ k ∈ Finset.Icc 4 G.natDegree,
            MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
              generalGaugeT ^ (m + 2) * MvPolynomial.X 0 ^ (k - 2) *
                generalGaugeQ G ^ (k + m)) =
        MvPolynomial.eval₂ (algebraMap K A) point
          (∑ k ∈ Finset.Icc 4 G.natDegree,
            MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
              generalGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
                generalGaugeQ G ^ k) := by
    simp only [MvPolynomial.eval₂_sum]
    exact Finset.sum_congr rfl (fun k _ => hterm k)
  rw [powerShiftedGaugeB, generalGaugeB]
  simp only [MvPolynomial.eval₂_add]
  rw [hsum]

/-- A common power shift is invisible in the third coordinate on
`Π = 1`. -/
theorem eval₂_powerShiftedGaugeC_eq_generalGaugeC
    (G : K[X]) (m : ℕ) (point : Fin 3 → A)
    (hpi : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugePi G) = 1) :
    MvPolynomial.eval₂ (algebraMap K A) point
        (powerShiftedGaugeC G m) =
      MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeC G) := by
  have hterm (k : ℕ) :
      MvPolynomial.eval₂ (algebraMap K A) point
          (MvPolynomial.C
              (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
            generalGaugeT ^ m * MvPolynomial.X 0 ^ k *
              generalGaugeQ G ^ (k + m)) =
        MvPolynomial.eval₂ (algebraMap K A) point
          (MvPolynomial.C
              (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
            (MvPolynomial.X 0 * generalGaugeQ G) ^ k) := by
    simp only [MvPolynomial.eval₂_mul, MvPolynomial.eval₂_pow,
      MvPolynomial.eval₂_X, MvPolynomial.eval₂_C]
    have hs := shifted_power_product_x
      (MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT)
      (MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeQ G))
      (point 0) k m (by
        simpa [generalGaugePi] using hpi)
    calc
      (algebraMap K A)
            (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
            MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT ^ m *
            point 0 ^ k *
            MvPolynomial.eval₂ (algebraMap K A) point
              (generalGaugeQ G) ^ (k + m) =
          (algebraMap K A)
              (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
            (MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT ^ m *
              point 0 ^ k *
              MvPolynomial.eval₂ (algebraMap K A) point
                (generalGaugeQ G) ^ (k + m)) := by ring
      _ = (algebraMap K A)
              (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
            (point 0 *
              MvPolynomial.eval₂ (algebraMap K A) point
                (generalGaugeQ G)) ^ k := by rw [hs]
  have hsum :
      MvPolynomial.eval₂ (algebraMap K A) point
          (∑ k ∈ Finset.Icc 4 G.natDegree,
            MvPolynomial.C
                (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
              generalGaugeT ^ m * MvPolynomial.X 0 ^ k *
                generalGaugeQ G ^ (k + m)) =
        MvPolynomial.eval₂ (algebraMap K A) point
          (∑ k ∈ Finset.Icc 4 G.natDegree,
            MvPolynomial.C
                (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
              (MvPolynomial.X 0 * generalGaugeQ G) ^ k) := by
    simp only [MvPolynomial.eval₂_sum]
    exact Finset.sum_congr rfl (fun k _ => hterm k)
  rw [powerShiftedGaugeC, generalGaugeC]
  simp only [MvPolynomial.eval₂_sub]
  rw [hsum]

/-- A cubic lift is invisible in the second coordinate on `Π = 1`. -/
theorem eval₂_cubicLiftGaugeB_eq_generalGaugeB
    (G : K[X]) (n : ℕ) (hn : 4 ≤ n) (point : Fin 3 → A)
    (hpi : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugePi G) = 1) :
    MvPolynomial.eval₂ (algebraMap K A) point
        (cubicLiftGaugeB G n) =
      MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeB G) := by
  have hn₁ : n - 1 = (n - 3) + 2 := by omega
  have hn₃ : n = 3 + (n - 3) := by omega
  have hnsub : 3 + (n - 3) - 3 + 2 = (n - 3) + 2 := by omega
  have hs := shifted_power_product
    (MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT)
    (MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeQ G))
    (point 0) 2 3 (n - 3) (by
      simpa [generalGaugePi] using hpi)
  rw [cubicLiftGaugeB]
  simp only [MvPolynomial.eval₂_add, MvPolynomial.eval₂_mul,
    MvPolynomial.eval₂_sub, MvPolynomial.eval₂_pow,
    MvPolynomial.eval₂_X, MvPolynomial.eval₂_C]
  rw [hn₁, hn₃, hnsub, hs]
  ring

/-- A cubic lift is invisible in the third coordinate on `Π = 1`. -/
theorem eval₂_cubicLiftGaugeC_eq_generalGaugeC
    (G : K[X]) (n : ℕ) (hn : 4 ≤ n) (point : Fin 3 → A)
    (hpi : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugePi G) = 1) :
    MvPolynomial.eval₂ (algebraMap K A) point
        (cubicLiftGaugeC G n) =
      MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeC G) := by
  have hn₃ : n = 3 + (n - 3) := by omega
  have hnsub : 3 + (n - 3) - 3 = n - 3 := by omega
  have hs := shifted_power_product
    (MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT)
    (MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeQ G))
    (point 0 ^ 3) 0 3 (n - 3) (by
      simpa [generalGaugePi] using hpi)
  have hs' :
      MvPolynomial.eval₂ (algebraMap K A) point generalGaugeT ^ (n - 3) *
          point 0 ^ 3 *
          MvPolynomial.eval₂ (algebraMap K A) point
            (generalGaugeQ G) ^ (3 + (n - 3)) =
        point 0 ^ 3 *
          MvPolynomial.eval₂ (algebraMap K A) point
            (generalGaugeQ G) ^ 3 := by
    simpa only [Nat.add_zero, pow_zero, one_mul] using hs
  rw [cubicLiftGaugeC]
  simp only [MvPolynomial.eval₂_sub, MvPolynomial.eval₂_mul,
    MvPolynomial.eval₂_pow, MvPolynomial.eval₂_X,
    MvPolynomial.eval₂_C]
  rw [hn₃, hnsub, hs']
  ring

/-- A literal point over target `(1,0,c)` for any three-coordinate
polynomial map. -/
@[ext]
structure StableGaugeFiberPoint
    (F : Fin 3 → GaugePolynomial K) (c : K)
    (A : Type*) [CommRing A] [Algebra K A] where
  point : Fin 3 → A
  pi_eq : MvPolynomial.eval₂ (algebraMap K A) point (F 0) = 1
  b_eq : MvPolynomial.eval₂ (algebraMap K A) point (F 1) = 0
  c_eq : MvPolynomial.eval₂ (algebraMap K A) point (F 2) = algebraMap K A c

/-- Evaluation commutes with a morphism of test algebras. -/
private theorem eval₂_map_algHom
    (f : A →ₐ[K] B) (point : Fin 3 → A) (P : GaugePolynomial K) :
    MvPolynomial.eval₂ (algebraMap K B) (fun i => f (point i)) P =
      f (MvPolynomial.eval₂ (algebraMap K A) point P) := by
  have hcomp : f.toRingHom.comp (algebraMap K A) = algebraMap K B := by
    ext r
    exact f.commutes r
  calc
    MvPolynomial.eval₂ (algebraMap K B) (fun i => f (point i)) P =
        MvPolynomial.eval₂ (f.toRingHom.comp (algebraMap K A))
          (fun i => f (point i)) P := by rw [hcomp]
    _ = f (MvPolynomial.eval₂ (algebraMap K A) point P) :=
      (MvPolynomial.hom_eval₂ P (algebraMap K A) f.toRingHom point).symm

namespace StableGaugeFiberPoint

variable {F : Fin 3 → GaugePolynomial K} {c : K}

/-- Stable literal fibers are functorial in the test algebra. -/
def map (f : A →ₐ[K] B) (p : StableGaugeFiberPoint F c A) :
    StableGaugeFiberPoint F c B where
  point := fun i => f (p.point i)
  pi_eq := by
    rw [eval₂_map_algHom f p.point, p.pi_eq, map_one]
  b_eq := by
    rw [eval₂_map_algHom f p.point, p.b_eq, map_zero]
  c_eq := by
    rw [eval₂_map_algHom f p.point, p.c_eq]
    exact f.commutes c

end StableGaugeFiberPoint

/-- On `Π = 1`, the normalized common power shift and the normalized
undeformed map have identical evaluations. -/
theorem eval₂_powerShiftedGaugeJacobianOneMap_eq
    (G : K[X]) (m : ℕ) (point : Fin 3 → A)
    (hpi : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugePi G) = 1) :
    eval₂Map (powerShiftedGaugeJacobianOneMap G m) point =
      eval₂Map (generalGaugeJacobianOneMap G) point := by
  funext i
  fin_cases i
  · simp [eval₂Map, powerShiftedGaugeJacobianOneMap,
      generalGaugeJacobianOneMap, scaleOutput, powerShiftedGaugeMap,
      generalGaugeMap]
  · change MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (-1 / 2 : K) * powerShiftedGaugeB G m) =
      MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (-1 / 2 : K) * generalGaugeB G)
    simp only [MvPolynomial.eval₂_mul, MvPolynomial.eval₂_C]
    rw [eval₂_powerShiftedGaugeB_eq_generalGaugeB G m point hpi]
  · change MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (1 : K) * powerShiftedGaugeC G m) =
      MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (1 : K) * generalGaugeC G)
    simp only [MvPolynomial.eval₂_mul, MvPolynomial.eval₂_C, map_one,
      one_mul]
    exact eval₂_powerShiftedGaugeC_eq_generalGaugeC G m point hpi

/-- On `Π = 1`, a normalized cubic lift and the normalized undeformed map
have identical evaluations. -/
theorem eval₂_cubicLiftGaugeJacobianOneMap_eq
    (G : K[X]) (n : ℕ) (hn : 4 ≤ n) (point : Fin 3 → A)
    (hpi : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugePi G) = 1) :
    eval₂Map (cubicLiftGaugeJacobianOneMap G n) point =
      eval₂Map (generalGaugeJacobianOneMap G) point := by
  funext i
  fin_cases i
  · simp [eval₂Map, cubicLiftGaugeJacobianOneMap,
      generalGaugeJacobianOneMap, scaleOutput, cubicLiftGaugeMap,
      generalGaugeMap]
  · change MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (-1 / 2 : K) * cubicLiftGaugeB G n) =
      MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (-1 / 2 : K) * generalGaugeB G)
    simp only [MvPolynomial.eval₂_mul, MvPolynomial.eval₂_C]
    rw [eval₂_cubicLiftGaugeB_eq_generalGaugeB G n hn point hpi]
  · change MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (1 : K) * cubicLiftGaugeC G n) =
      MvPolynomial.eval₂ (algebraMap K A) point
        (MvPolynomial.C (1 : K) * generalGaugeC G)
    simp only [MvPolynomial.eval₂_mul, MvPolynomial.eval₂_C, map_one,
      one_mul]
    exact eval₂_cubicLiftGaugeC_eq_generalGaugeC G n hn point hpi

/-- The selected literal fiber of a common power shift is canonically the
selected literal fiber of the undeformed normalized gauge map. -/
def generalGaugeJacobianOneFiberEquivPowerShifted
    (G : K[X]) (m : ℕ) (c : K) :
    GeneralGaugeJacobianOneFiberPoint G 1 c A ≃
      StableGaugeFiberPoint (powerShiftedGaugeJacobianOneMap G m) c A where
  toFun := fun p => by
    have hpi :
        MvPolynomial.eval₂ (algebraMap K A) p.point
          (generalGaugePi G) = 1 := by
      simpa using p.pi_eq
    have heval := eval₂_powerShiftedGaugeJacobianOneMap_eq
      G m p.point hpi
    exact
      { point := p.point
        pi_eq := by
          change eval₂Map (powerShiftedGaugeJacobianOneMap G m) p.point 0 = 1
          rw [heval]
          change MvPolynomial.eval₂ (algebraMap K A) p.point
            (generalGaugeJacobianOneMap G 0) = 1
          simpa using p.pi_eq
        b_eq := by
          change eval₂Map (powerShiftedGaugeJacobianOneMap G m) p.point 1 = 0
          rw [heval]
          exact p.b_eq
        c_eq := by
          change eval₂Map (powerShiftedGaugeJacobianOneMap G m) p.point 2 =
            algebraMap K A c
          rw [heval]
          exact p.c_eq }
  invFun := fun p => by
    have hpi :
        MvPolynomial.eval₂ (algebraMap K A) p.point
          (generalGaugePi G) = 1 := by
      simpa [powerShiftedGaugeJacobianOneMap, scaleOutput,
        powerShiftedGaugeMap] using p.pi_eq
    have heval := eval₂_powerShiftedGaugeJacobianOneMap_eq
      G m p.point hpi
    exact
      { point := p.point
        pi_eq := by
          change eval₂Map (generalGaugeJacobianOneMap G) p.point 0 =
            algebraMap K A ((1 : Kˣ) : K)
          rw [← heval]
          simpa only [eval₂Map, Units.val_one, map_one] using p.pi_eq
        b_eq := by
          change eval₂Map (generalGaugeJacobianOneMap G) p.point 1 = 0
          rw [← heval]
          exact p.b_eq
        c_eq := by
          change eval₂Map (generalGaugeJacobianOneMap G) p.point 2 =
            algebraMap K A c
          rw [← heval]
          exact p.c_eq }
  left_inv := by
    intro p
    apply GeneralGaugeJacobianOneFiberPoint.ext
    rfl
  right_inv := by
    intro p
    apply StableGaugeFiberPoint.ext
    rfl

/-- The selected literal fiber of a cubic lift is canonically the selected
literal fiber of the undeformed normalized gauge map. -/
def generalGaugeJacobianOneFiberEquivCubicLift
    (G : K[X]) (n : ℕ) (hn : 4 ≤ n) (c : K) :
    GeneralGaugeJacobianOneFiberPoint G 1 c A ≃
      StableGaugeFiberPoint (cubicLiftGaugeJacobianOneMap G n) c A where
  toFun := fun p => by
    have hpi :
        MvPolynomial.eval₂ (algebraMap K A) p.point
          (generalGaugePi G) = 1 := by
      simpa using p.pi_eq
    have heval := eval₂_cubicLiftGaugeJacobianOneMap_eq
      G n hn p.point hpi
    exact
      { point := p.point
        pi_eq := by
          change eval₂Map (cubicLiftGaugeJacobianOneMap G n) p.point 0 = 1
          rw [heval]
          change MvPolynomial.eval₂ (algebraMap K A) p.point
            (generalGaugeJacobianOneMap G 0) = 1
          simpa using p.pi_eq
        b_eq := by
          change eval₂Map (cubicLiftGaugeJacobianOneMap G n) p.point 1 = 0
          rw [heval]
          exact p.b_eq
        c_eq := by
          change eval₂Map (cubicLiftGaugeJacobianOneMap G n) p.point 2 =
            algebraMap K A c
          rw [heval]
          exact p.c_eq }
  invFun := fun p => by
    have hpi :
        MvPolynomial.eval₂ (algebraMap K A) p.point
          (generalGaugePi G) = 1 := by
      simpa [cubicLiftGaugeJacobianOneMap, scaleOutput,
        cubicLiftGaugeMap] using p.pi_eq
    have heval := eval₂_cubicLiftGaugeJacobianOneMap_eq
      G n hn p.point hpi
    exact
      { point := p.point
        pi_eq := by
          change eval₂Map (generalGaugeJacobianOneMap G) p.point 0 =
            algebraMap K A ((1 : Kˣ) : K)
          rw [← heval]
          simpa only [eval₂Map, Units.val_one, map_one] using p.pi_eq
        b_eq := by
          change eval₂Map (generalGaugeJacobianOneMap G) p.point 1 = 0
          rw [← heval]
          exact p.b_eq
        c_eq := by
          change eval₂Map (generalGaugeJacobianOneMap G) p.point 2 =
            algebraMap K A c
          rw [← heval]
          exact p.c_eq }
  left_inv := by
    intro p
    apply GeneralGaugeJacobianOneFiberPoint.ext
    rfl
  right_inv := by
    intro p
    apply StableGaugeFiberPoint.ext
    rfl

/-- The power-shift fiber equivalence commutes with change of test algebra. -/
theorem generalGaugeJacobianOneFiberEquivPowerShifted_natural
    (G : K[X]) (m : ℕ) (c : K) (f : A →ₐ[K] B)
    (p : GeneralGaugeJacobianOneFiberPoint G 1 c A) :
    StableGaugeFiberPoint.map f
        (generalGaugeJacobianOneFiberEquivPowerShifted
          (A := A) G m c p) =
      generalGaugeJacobianOneFiberEquivPowerShifted
        (A := B) G m c (GeneralGaugeJacobianOneFiberPoint.map f p) := by
  apply StableGaugeFiberPoint.ext
  rfl

/-- The cubic-lift fiber equivalence commutes with change of test algebra. -/
theorem generalGaugeJacobianOneFiberEquivCubicLift_natural
    (G : K[X]) (n : ℕ) (hn : 4 ≤ n) (c : K) (f : A →ₐ[K] B)
    (p : GeneralGaugeJacobianOneFiberPoint G 1 c A) :
    StableGaugeFiberPoint.map f
        (generalGaugeJacobianOneFiberEquivCubicLift
          (A := A) G n hn c p) =
      generalGaugeJacobianOneFiberEquivCubicLift
        (A := B) G n hn c
          (GeneralGaugeJacobianOneFiberPoint.map f p) := by
  apply StableGaugeFiberPoint.ext
  rfl

/-- The original quotient algebra represents every selected common
power-shifted realization fiber. -/
def powerShiftedGaugeRealizationFiberRepresentingEquiv
    (P : K[X]) (a : K) (m : ℕ) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (AdjoinRoot P →ₐ[K] A) ≃
      StableGaugeFiberPoint
        (powerShiftedGaugeJacobianOneMap (realizationSeed P a) m)
        (realizationTargetC P a (P.derivative.eval a)) A :=
  (realizationJacobianOneFiberRepresentingEquiv
    (A := A) P a hP h₁ h₃).trans
      (generalGaugeJacobianOneFiberEquivPowerShifted
        (A := A) (realizationSeed P a) m
        (realizationTargetC P a (P.derivative.eval a)))

/-- For a cubic polynomial, the original quotient algebra represents every
selected cubic-lift realization fiber. -/
def cubicLiftGaugeRealizationFiberRepresentingEquiv
    (P : K[X]) (a : K) (n : ℕ) (hn : 4 ≤ n)
    (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (AdjoinRoot P →ₐ[K] A) ≃
      StableGaugeFiberPoint
        (cubicLiftGaugeJacobianOneMap (realizationSeed P a) n)
        (realizationTargetC P a (P.derivative.eval a)) A :=
  (realizationJacobianOneFiberRepresentingEquiv
    (A := A) P a hP h₁ h₃).trans
      (generalGaugeJacobianOneFiberEquivCubicLift
        (A := A) (realizationSeed P a) n hn
        (realizationTargetC P a (P.derivative.eval a)))

/-- Naturality of the represented common power-shifted realization fiber. -/
theorem powerShiftedGaugeRealizationFiberRepresentingEquiv_natural
    (P : K[X]) (a : K) (m : ℕ) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0)
    (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A) :
    StableGaugeFiberPoint.map f
        (powerShiftedGaugeRealizationFiberRepresentingEquiv
          (A := A) P a m hP h₁ h₃ φ) =
      powerShiftedGaugeRealizationFiberRepresentingEquiv
        (A := B) P a m hP h₁ h₃ (f.comp φ) := by
  apply StableGaugeFiberPoint.ext
  funext i
  have h := realizationJacobianOneFiberRepresentingEquiv_natural
    P a hP h₁ h₃ f φ
  exact congrArg (fun p => p.point i) h

/-- Naturality of the represented cubic-lift realization fiber. -/
theorem cubicLiftGaugeRealizationFiberRepresentingEquiv_natural
    (P : K[X]) (a : K) (n : ℕ) (hn : 4 ≤ n)
    (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0)
    (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A) :
    StableGaugeFiberPoint.map f
        (cubicLiftGaugeRealizationFiberRepresentingEquiv
          (A := A) P a n hn hP h₁ h₃ φ) =
      cubicLiftGaugeRealizationFiberRepresentingEquiv
        (A := B) P a n hn hP h₁ h₃ (f.comp φ) := by
  apply StableGaugeFiberPoint.ext
  funext i
  have h := realizationJacobianOneFiberRepresentingEquiv_natural
    P a hP h₁ h₃ f φ
  exact congrArg (fun p => p.point i) h

#print axioms eval₂_powerShiftedGaugeJacobianOneMap_eq
#print axioms eval₂_cubicLiftGaugeJacobianOneMap_eq
#print axioms powerShiftedGaugeRealizationFiberRepresentingEquiv
#print axioms cubicLiftGaugeRealizationFiberRepresentingEquiv
#print axioms powerShiftedGaugeRealizationFiberRepresentingEquiv_natural
#print axioms cubicLiftGaugeRealizationFiberRepresentingEquiv_natural

end FiniteEtaleKeller
