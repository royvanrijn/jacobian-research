# Lean formalization: finite étale Keller fibers

This project formalizes the scheme-theoretic core of *Every Finite Étale
Algebra Except Rank Two Is a Keller Fiber*. It uses Lean `v4.33.0-rc1` and
Mathlib at the matching release candidate.

## Proof status

| Stage | Scope | Status |
|---|---|---|
| 1 | Explicit quintic map, output scaling, and Bézout inverse | implemented |
| 2 | Universal marked-line identities and Jacobian cancellation core | implemented |
| 3 | Two-sided source/chart reconstruction over arbitrary commutative rings | implemented |
| 4 | Roots versus full source-fiber points, including naturality | implemented |
| 5 | Representation by `K[S]/(E)` and transport to `K[T]/(P)` | implemented |
| 6 | Existence and automatic choice of an admissible translation | implemented |
| 7 | Polynomial-level represented-fiber theorem with no supplied parameter | implemented |
| 8 | Coefficientwise all-degree gauge assembly identities | implemented |
| 9 | One multivariate object for the full all-degree map and its degree bound | paper and exact checker; not yet Lean |
| 10 | Monogenicity and the complete rank classification | paper proof; not yet Lean |
| 11 | Historical degree-two Galois exclusion | external theorem; not yet Lean |

## Central formal theorem

For a field `K`, a separable polynomial `E`, a unit `g₁`, and the derivative
factorization required by the quadratic gauge, the project constructs, for every
commutative test `K`-algebra `A`, an equivalence

```text
(AdjoinRoot E →ₐ[K] A) ≃ GaugeFiberPoint E β pi b a A
```

and proves that it commutes with every algebra homomorphism `A →ₐ[K] B`. Thus
the source-fiber functor is naturally represented by `K[S]/(E)`; this is not
only a bijection on field-valued points.

`RealizationFiber.lean` instantiates this datum for a supplied admissible
translation parameter:

```text
E(S) = P(a + S),
```

and composes the represented-fiber equivalence with the canonical translation
of quotients

```text
K[S]/(P(a+S)) ≃ K[T]/(P(T)).
```

`Admissibility.lean` then proves that every polynomial of degree at least three
over a characteristic-zero field has a parameter where both its first and
third Hasse derivatives evaluate nontrivially. `AutomaticRealization.lean`
chooses such a parameter internally. Its final theorem is

```text
automaticFiberRepresentingEquiv_natural
```

and requires only a squarefree polynomial `P` and `3 ≤ P.natDegree`. It gives,
naturally in every commutative test algebra `A`, a represented-fiber equivalence
from `AdjoinRoot P`; no translation parameter or nonvanishing witness remains
as an external hypothesis.

`GaugeAssembly.lean` formalizes the coefficientwise algebra used to turn the
marked coordinates into the displayed all-degree polynomial map. It proves the
high-degree monomial transport identities, the low-degree second-coordinate
combination, and the complete cubic cancellation in the third coordinate over
an arbitrary commutative ring. This closes the algebraic cancellation that was
previously present only in the manuscript and representative symbolic checks.

The project contains no `sorry` and introduces no project-specific axioms.
`#print axioms` for the final represented-fiber theorem reports only the
standard Lean foundations `propext`, `Classical.choice`, and `Quot.sound`; the
new gauge-assembly certificate is proved by ring algebra and unit identities.

## What the formalization clarifies

The reconstruction works over arbitrary commutative test algebras. Units are
carried explicitly, so there is no hidden localization and no omitted component.
Separability makes the derivative class invertible by Bézout; the unit first
target coordinate makes the source chart global on the fiber.

The abstract represented-fiber argument needs a fixed source coefficient and a
derivative factorization. Nonvanishing of the cubic Taylor coefficient is used
to connect that abstract coefficient to the displayed polynomial gauge
`g₁/g₃`; it is not an additional scheme-reconstruction step. The characteristic-
free cubic coefficient is the third Hasse derivative, equal to `P'''(a)/3!` in
characteristic zero.

The admissible-parameter proof is also logically independent of squarefreeness:
degree at least three makes `P'` and the third Hasse derivative nonzero
polynomials, and a nonzero product cannot vanish at every point of an infinite
field. Squarefreeness enters later, where it supplies separability of the
translated inverse polynomial.

## Scope boundary

The Lean certificate now covers admissible-parameter existence, difficult
scheme structure, reconstruction, naturality, representability, quotient
translation, and the uniform coefficientwise cancellations that assemble the
second and third displayed coordinates. It does not yet package the finite
coefficient sums as one general `MvPolynomial` map, prove that complete object's
Jacobian and `6N+2` coordinate-degree bound in Lean, formalize monogenicity, or
internalize the Campbell--Razar--Wright rank-two exclusion. The explicit optimal
quintic map and its determinant-one normalization are formalized separately.

## Build

```bash
cd formal/finite-etale-keller
lake build
```

Repository CI builds this project independently of the external Lean certificate
for the foundational three-dimensional map.
