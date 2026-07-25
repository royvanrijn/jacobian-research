/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Reconstruction

/-!
# Two-sided source-chart reconstruction

The fiber proof repeatedly passes between source coordinates `(x,y,z)` and the
marked-line coordinates `(S,Q,D)`.  This module packages that passage as an
actual equivalence over an arbitrary commutative ring.  Units are carried as
`Rˣ`, so the statement contains no field division and no hidden localization.
-/

noncomputable section

namespace FiniteEtaleKeller

variable {R : Type*} [CommRing R]

/-- Marked-line chart data.  The unit `d` is the chart factor
`1 - S*Q + pi*S^2`. -/
@[ext]
structure GaugeChart (R : Type*) [CommRing R] (pi : R) where
  S : R
  Q : R
  d : Rˣ
  chart_eq : (d : R) = 1 - S * Q + pi * S ^ 2

/-- Source data on the full fiber over a unit first target coordinate.  The
source polynomial `q` is derived from `(t,y,z)` below, so the structure keeps
only the genuine source coordinates and their two defining relations. -/
@[ext]
structure GaugeSource (R : Type*) [CommRing R] (pi a : R) where
  t : Rˣ
  x : R
  y : R
  z : R
  t_eq : (t : R) = 1 + x * y
  pi_eq : (t : R) *
    ((t : R) ^ 2 * z + a * y ^ 2 * (1 + 3 * (t : R))) = pi

namespace GaugeSource

variable {pi a : R}

/-- The source polynomial `q`. -/
def q (p : GaugeSource R pi a) : R :=
  (p.t : R) ^ 2 * p.z + a * p.y ^ 2 * (1 + 3 * (p.t : R))

/-- The global marked coordinate `S=x/t`, written using the inverse of the
source unit. -/
def S (p : GaugeSource R pi a) : R :=
  p.x * (↑p.t⁻¹ : R)

/-- The global marked coordinate `Q=y+xq`. -/
def Q (p : GaugeSource R pi a) : R :=
  p.y + p.x * p.q

@[simp]
theorem t_mul_q (p : GaugeSource R pi a) : (p.t : R) * p.q = pi :=
  p.pi_eq

end GaugeSource

namespace GaugeChart

variable {pi : R}

/-- Reconstruct the source coordinates from marked-line chart data. -/
def toSource (p : GaugeChart R pi) (a : R) : GaugeSource R pi a where
  t := p.d⁻¹
  x := p.S * (↑p.d⁻¹ : R)
  y := p.Q - pi * p.S
  z := (p.d : R) ^ 2 *
    (pi * (p.d : R)
      - a * (p.Q - pi * p.S) ^ 2 * (1 + 3 * (↑p.d⁻¹ : R)))
  t_eq := by
    have h := unitReconstruction_identities p.S p.Q pi a p.d p.chart_eq
    dsimp only at h
    exact h.1.symm
  pi_eq := by
    have h := unitReconstruction_identities p.S p.Q pi a p.d p.chart_eq
    dsimp only at h
    calc
      (↑p.d⁻¹ : R) *
          ((↑p.d⁻¹ : R) ^ 2 *
              ((p.d : R) ^ 2 *
                (pi * (p.d : R)
                  - a * (p.Q - pi * p.S) ^ 2 *
                    (1 + 3 * (↑p.d⁻¹ : R))))
            + a * (p.Q - pi * p.S) ^ 2 *
                (1 + 3 * (↑p.d⁻¹ : R)))
          = (↑p.d⁻¹ : R) * (pi * (p.d : R)) := by rw [h.2.1]
      _ = pi := h.2.2.2.2.1

@[simp]
theorem toSource_t (p : GaugeChart R pi) (a : R) : (p.toSource a).t = p.d⁻¹ := rfl

@[simp]
theorem toSource_x (p : GaugeChart R pi) (a : R) :
    (p.toSource a).x = p.S * (↑p.d⁻¹ : R) := rfl

@[simp]
theorem toSource_y (p : GaugeChart R pi) (a : R) :
    (p.toSource a).y = p.Q - pi * p.S := rfl

@[simp]
theorem toSource_z (p : GaugeChart R pi) (a : R) :
    (p.toSource a).z = (p.d : R) ^ 2 *
      (pi * (p.d : R)
        - a * (p.Q - pi * p.S) ^ 2 * (1 + 3 * (↑p.d⁻¹ : R))) := rfl

end GaugeChart

namespace GaugeSource

variable {pi a : R}

/-- The reciprocal chart identity derived from the two source equations. -/
theorem toChart_chart_eq (p : GaugeSource R pi a) :
    (↑p.t⁻¹ : R) = 1 - p.S * p.Q + pi * p.S ^ 2 := by
  let d : R := ↑p.t⁻¹
  have hunit : d * (p.t : R) = 1 := by simp [d]
  change d = 1 - (p.x * d) * (p.y + p.x * p.q) + pi * (p.x * d) ^ 2
  linear_combination
    (1 - p.x ^ 2 * d * p.q) * hunit
      - d * p.t_eq
      + d ^ 2 * p.x ^ 2 * p.t_mul_q

/-- Recover marked-line chart data from a source point. -/
def toChart (p : GaugeSource R pi a) : GaugeChart R pi where
  S := p.S
  Q := p.Q
  d := p.t⁻¹
  chart_eq := p.toChart_chart_eq

@[simp]
theorem toChart_S (p : GaugeSource R pi a) : p.toChart.S = p.S := rfl

@[simp]
theorem toChart_Q (p : GaugeSource R pi a) : p.toChart.Q = p.Q := rfl

@[simp]
theorem toChart_d (p : GaugeSource R pi a) : p.toChart.d = p.t⁻¹ := rfl

end GaugeSource

/-- Reconstructing and then returning to marked-line coordinates is the
identity. -/
theorem GaugeChart.toSource_toChart {pi : R} (p : GaugeChart R pi) (a : R) :
    (p.toSource a).toChart = p := by
  have h := unitReconstruction_identities p.S p.Q pi a p.d p.chart_eq
  dsimp only at h
  apply GaugeChart.ext
  · simpa [GaugeSource.S, GaugeChart.toSource] using h.2.2.1
  · simpa [GaugeSource.Q, GaugeSource.q, GaugeChart.toSource] using h.2.2.2.1
  · rw [GaugeSource.toChart_d, GaugeChart.toSource_t]
    simp

/-- Passing to marked-line coordinates and reconstructing recovers every
source coordinate, including `z`. -/
theorem GaugeSource.toChart_toSource {pi a : R} (p : GaugeSource R pi a) :
    p.toChart.toSource a = p := by
  let d : R := ↑p.t⁻¹
  have htd : (p.t : R) * d = 1 := by simp [d]
  have hdt : d * (p.t : R) = 1 := by simp [d]
  have hpi : pi = (p.t : R) * p.q := p.t_mul_q.symm
  have hq : pi * d = p.q := by
    calc
      pi * d = ((p.t : R) * p.q) * d := congrArg (fun r : R => r * d) hpi
      _ = p.q * ((p.t : R) * d) := by ring
      _ = p.q := by rw [htd, mul_one]
  have hy : p.Q - pi * p.S = p.y := by
    simp only [GaugeSource.Q, GaugeSource.S]
    change p.y + p.x * p.q - pi * (p.x * d) = p.y
    calc
      p.y + p.x * p.q - pi * (p.x * d) =
          p.y + p.x * p.q - ((p.t : R) * p.q) * (p.x * d) :=
        congrArg (fun r : R => p.y + p.x * p.q - r * (p.x * d)) hpi
      _ = p.y + p.x * p.q * (1 - (p.t : R) * d) := by ring
      _ = p.y := by rw [htd]; ring
  have hsq : d ^ 2 * (p.t : R) ^ 2 = 1 := by
    rw [← mul_pow, hdt, one_pow]
  have hdiff :
      p.q - a * p.y ^ 2 * (1 + 3 * (p.t : R)) = (p.t : R) ^ 2 * p.z := by
    simp only [GaugeSource.q]
    ring
  apply GaugeSource.ext
  · rw [GaugeChart.toSource_t, GaugeSource.toChart_d]
    simp
  · rw [GaugeChart.toSource_x, GaugeSource.toChart_S, GaugeSource.toChart_d]
    simp [GaugeSource.S, mul_assoc]
  · rw [GaugeChart.toSource_y, GaugeSource.toChart_Q, GaugeSource.toChart_S]
    exact hy
  · rw [GaugeChart.toSource_z, GaugeSource.toChart_d,
      GaugeSource.toChart_Q, GaugeSource.toChart_S]
    change d ^ 2 *
        (pi * d - a * (p.Q - pi * p.S) ^ 2 * (1 + 3 * (p.t : R))) = p.z
    rw [hq, hy, hdiff]
    calc
      d ^ 2 * ((p.t : R) ^ 2 * p.z) =
          (d ^ 2 * (p.t : R) ^ 2) * p.z := by ring
      _ = p.z := by rw [hsq, one_mul]

/-- The source chart and the reconstruction chart are mutually inverse over
any commutative ring. -/
def gaugeChartSourceEquiv (pi a : R) :
    GaugeChart R pi ≃ GaugeSource R pi a where
  toFun := fun p => p.toSource a
  invFun := GaugeSource.toChart
  left_inv := fun p => p.toSource_toChart a
  right_inv := GaugeSource.toChart_toSource

end FiniteEtaleKeller
