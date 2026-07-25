# Lean formalization: finite étale Keller fibers

This project formalizes the scheme-theoretic core of *Finite Étale Algebras as
Keller Fibers*. It uses Lean `v4.33.0-rc1` and Mathlib at the matching release
candidate.

## Proof status

| Stage | Scope | Status |
|---|---|---|
| 1 | Explicit quintic map, output scaling, and Bézout inverse | implemented |
| 2 | Universal marked-line identities and Jacobian cancellation core | implemented |
| 3 | Two-sided source/chart reconstruction over arbitrary commutative rings | implemented |
| 4 | Roots versus full source-fiber points, including naturality | implemented |
| 5 | Representation by `K[S]/(E)` and transport to `K[T]/(P)` | implemented |
| 6 | All-degree displayed map assembly and coordinate-degree bound | paper and exact checker; partial Lean |
| 7 | Monogenicity and the complete rank classification | paper proof; not yet Lean |
| 8 | Historical degree-two Galois exclusion | external theorem; not yet Lean |

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

`RealizationFiber.lean` then instantiates the datum for

```text
E(S) = P(a + S)
```

and composes the represented-fiber equivalence with the canonical translation
of quotients

```text
K[S]/(P(a+S)) ≃ K[T]/(P(T)).
```

The resulting theorem is
`translatedFiberRepresentingEquiv_natural`.

The project contains no `sorry` and introduces no project-specific axioms.
`#print axioms` for the final theorem reports only the standard Lean foundations
`propext`, `Classical.choice`, and `Quot.sound`.

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

## Scope boundary

The Lean certificate now covers the difficult scheme-structure, reconstruction,
naturality, representability, and translation steps. It does not yet package the
entire general displayed polynomial map, its all-degree determinant calculation,
the `6N+2` degree estimate, monogenicity, or the Campbell--Razar--Wright
rank-two exclusion into one Lean theorem. The explicit optimal quintic map and
its determinant-one normalization are formalized separately.

## Build

```bash
cd formal/finite-etale-keller
lake build
```

Repository CI builds this project independently of the external Lean certificate
for the foundational three-dimensional map.
