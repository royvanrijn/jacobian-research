/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.MarkedLine

/-!
# Coefficientwise assembly of the quadratic gauge

This module formalizes the uniform identities used to pass from the marked
coordinates

`S = x / t`, `Q = y + x*q`, `pi = t*q`

to the displayed all-degree polynomial map.  It isolates the two ingredients
that do not depend on a degree bound:

* every high-degree monomial loses exactly the required powers of `t`;
* the low-degree part of the third coordinate collapses by the cubic identity
  `t + 1 - (t - 1)^2 * (1 + 3*t) = t^2 * (5 - 3*t)`.

The finite polynomial sum and the complete multivariate Jacobian are left to a
later module, but the coefficientwise algebra used by the paper is proved here
over an arbitrary commutative ring.
-/

noncomputable section

namespace FiniteEtaleKeller

variable {R : Type*} [CommRing R]

/-- The cubic cancellation responsible for the low-degree part of the third
quadratic-gauge coordinate. -/
theorem quadraticGauge_cubicCancellation (t : R) :
    t + 1 - (t - 1) ^ 2 * (1 + 3 * t) = t ^ 2 * (5 - 3 * t) := by
  ring

/-- Under `S = x/t` and `pi = t*q`, a degree-`k` monomial in `pi*S`
becomes the polynomial monomial `(x*q)^k`. -/
theorem quadraticGauge_highCMonomial
    (t : Rˣ) (x q : R) (k : ℕ) :
    (((t : R) * q) ^ k) * (x * (↑t⁻¹ : R)) ^ k = (x * q) ^ k := by
  have hbase :
      ((t : R) * q) * (x * (↑t⁻¹ : R)) = x * q := by
    calc
      ((t : R) * q) * (x * (↑t⁻¹ : R)) =
          x * q * ((t : R) * (↑t⁻¹ : R)) := by ring
      _ = x * q := by simp
  calc
    (((t : R) * q) ^ k) * (x * (↑t⁻¹ : R)) ^ k =
        (((t : R) * q) * (x * (↑t⁻¹ : R))) ^ k := by
          rw [mul_pow]
    _ = (x * q) ^ k := by rw [hbase]

/-- Writing the exponent as `m+2`, the high-degree monomial occurring in the
second coordinate becomes `t^2*x^m*q^(m+2)`.  Taking `m = k-2` is the identity
used for every paper coefficient indexed by `k ≥ 4`. -/
theorem quadraticGauge_highBMonomial
    (t : Rˣ) (x q : R) (m : ℕ) :
    (((t : R) * q) ^ (m + 2)) * (x * (↑t⁻¹ : R)) ^ m =
      (t : R) ^ 2 * x ^ m * q ^ (m + 2) := by
  have hbase :
      ((t : R) * q) * (x * (↑t⁻¹ : R)) = x * q := by
    calc
      ((t : R) * q) * (x * (↑t⁻¹ : R)) =
          x * q * ((t : R) * (↑t⁻¹ : R)) := by ring
      _ = x * q := by simp
  calc
    (((t : R) * q) ^ (m + 2)) * (x * (↑t⁻¹ : R)) ^ m =
        ((((t : R) * q) ^ m) * (x * (↑t⁻¹ : R)) ^ m) *
          (((t : R) * q) ^ 2) := by
            rw [pow_add]
            ring
    _ = ((((t : R) * q) * (x * (↑t⁻¹ : R))) ^ m) *
          (((t : R) * q) ^ 2) := by
            rw [mul_pow]
    _ = (x * q) ^ m * (((t : R) * q) ^ 2) := by rw [hbase]
    _ = (t : R) ^ 2 * x ^ m * q ^ (m + 2) := by
      rw [mul_pow, mul_pow, pow_add]
      ring

/-- The linear and cubic contributions to the second coordinate combine with
`Q = y + x*q` to give the displayed coefficient `3*r*x*q`. -/
theorem quadraticGauge_lowBIdentity
    (t : Rˣ) (x y q c r : R) :
    (y + x * q) + 2 * c * ((t : R) * q) +
        (3 * r - 1) * ((t : R) * q) * (x * (↑t⁻¹ : R)) =
      y + 3 * r * x * q + 2 * c * ((t : R) * q) := by
  have hmarked :
      ((t : R) * q) * (x * (↑t⁻¹ : R)) = x * q := by
    calc
      ((t : R) * q) * (x * (↑t⁻¹ : R)) =
          x * q * ((t : R) * (↑t⁻¹ : R)) := by ring
      _ = x * q := by simp
  rw [hmarked]
  ring

/-- The complete low-degree part of the third coordinate.  Here `a` is the
source coefficient `g₁/g₃` and `r` is its inverse `g₃/g₁`. -/
theorem quadraticGauge_lowCIdentity
    (t : Rˣ) (x y z q a r : R)
    (ht : (t : R) = 1 + x * y)
    (hq : q = (t : R) ^ 2 * z + a * y ^ 2 * (1 + 3 * (t : R)))
    (hra : r * a = 1) :
    2 * (x * (↑t⁻¹ : R))
        - (y + x * q) * (x * (↑t⁻¹ : R)) ^ 2
        + (1 - r) * ((t : R) * q) * (x * (↑t⁻¹ : R)) ^ 3 =
      x * (5 - 3 * (t : R)) - r * x ^ 3 * z := by
  let tv : R := t
  let ti : R := ↑t⁻¹
  have hunit : tv * ti = 1 := by
    simp [tv, ti]
  have htv : tv = 1 + x * y := by
    simpa [tv] using ht
  have hxy : x * y = tv - 1 := by
    rw [htv]
    ring
  have hxy0 : tv - 1 - x * y = 0 := by
    rw [htv]
    ring
  have hsq : tv ^ 2 * ti ^ 2 = 1 := by
    calc
      tv ^ 2 * ti ^ 2 = (tv * ti) ^ 2 := by ring
      _ = 1 := by rw [hunit]; ring
  have hq' : q = tv ^ 2 * z + a * y ^ 2 * (1 + 3 * tv) := by
    simpa [tv] using hq
  have hx3y2 : x ^ 3 * y ^ 2 = x * (tv - 1) ^ 2 := by
    calc
      x ^ 3 * y ^ 2 = x * (x * y) ^ 2 := by ring
      _ = x * (tv - 1) ^ 2 := by rw [hxy]
  change
    2 * (x * ti) - (y + x * q) * (x * ti) ^ 2
        + (1 - r) * (tv * q) * (x * ti) ^ 3 =
      x * (5 - 3 * tv) - r * x ^ 3 * z
  calc
    2 * (x * ti) - (y + x * q) * (x * ti) ^ 2
          + (1 - r) * (tv * q) * (x * ti) ^ 3 =
        (2 * x * tv - x ^ 2 * y - r * x ^ 3 * q) * ti ^ 2 := by
          linear_combination
            (-2 * x * ti + (1 - r) * x ^ 3 * q * ti ^ 2) * hunit
    _ = (x * (tv + 1) - r * x ^ 3 * q) * ti ^ 2 := by
      linear_combination (x * ti ^ 2) * hxy0
    _ = (x * (tv + 1) -
          r * x ^ 3 * (tv ^ 2 * z + a * y ^ 2 * (1 + 3 * tv))) * ti ^ 2 := by
      rw [hq']
    _ = x * (tv + 1) * ti ^ 2
          - r * x ^ 3 * z * (tv ^ 2 * ti ^ 2)
          - (r * a) * (x ^ 3 * y ^ 2) * (1 + 3 * tv) * ti ^ 2 := by
      ring
    _ = x * (tv + 1) * ti ^ 2 - r * x ^ 3 * z
          - x * (tv - 1) ^ 2 * (1 + 3 * tv) * ti ^ 2 := by
      rw [hsq, hra, hx3y2]
      ring
    _ = -r * x ^ 3 * z +
          x * (tv + 1 - (tv - 1) ^ 2 * (1 + 3 * tv)) * ti ^ 2 := by
      ring
    _ = -r * x ^ 3 * z + x * (tv ^ 2 * (5 - 3 * tv)) * ti ^ 2 := by
      rw [quadraticGauge_cubicCancellation]
    _ = -r * x ^ 3 * z + x * (5 - 3 * tv) * (tv ^ 2 * ti ^ 2) := by
      ring
    _ = x * (5 - 3 * tv) - r * x ^ 3 * z := by
      rw [hsq]
      ring

#print axioms quadraticGauge_lowCIdentity

end FiniteEtaleKeller
