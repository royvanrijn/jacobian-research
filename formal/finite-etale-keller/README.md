# Lean formalization: finite étale Keller fibers

This project formalizes the polynomial-presentation construction in *Every
Nonzero Finite Étale Algebra Except Rank Two Is a Keller Fiber*. It uses Lean
`v4.33.0-rc1` and Mathlib at the matching release candidate.

## Proof status

| Stage | Scope | Status |
|---|---|---|
| 1 | Explicit quintic map, output scaling, and Bézout inverse | implemented |
| 2 | Universal marked-line identities and Jacobian cancellation core | implemented |
| 3 | Two-sided source/chart reconstruction over arbitrary commutative rings | implemented |
| 4 | Roots versus full abstract source-fiber points, including naturality | implemented |
| 5 | Representation by `K[S]/(E)` and transport to `K[T]/(P)` | implemented |
| 6 | Existence and automatic choice of an admissible translation | implemented |
| 7 | Complete finite-sum all-degree gauge assembly identities | implemented |
| 8 | One arbitrary-degree `MvPolynomial (Fin 3) K` map and exact evaluations | implemented |
| 9 | Generic inverse polynomial, explicit `β`, and derivative factorization | implemented |
| 10 | Actual displayed equations versus represented source equations | implemented |
| 11 | Literal raw map fiber and naturality over arbitrary test algebras | implemented |
| 12 | Determinant-one normalization of the literal fiber | implemented |
| 13 | General arbitrary-degree Jacobian `-2`, normalized Jacobian `1` | implemented |
| 14 | General and final `6N+2` coordinate-degree bounds | implemented |
| 15 | Final automatic actual-map realization from squarefree `P` | implemented |
| 16 | Generic inverse irreducibility and geometric degree | paper proof; not yet Lean |
| 17 | Monogenicity and the passage from arbitrary finite étale algebras | paper proof; not yet Lean |
| 18 | Historical degree-two Galois exclusion | external theorem; not yet Lean |

## Final polynomial-presentation theorem

Let `K` be a characteristic-zero field, let `P : K[X]` be squarefree, and
assume `3 ≤ P.natDegree`. Lean now chooses an admissible translation parameter
internally and defines

```text
automaticRealizationMap P hdeg : Fin 3 → MvPolynomial (Fin 3) K
```

together with its distinguished target. The formal development proves:

```text
jacobianDet (automaticRealizationMap P hdeg) = 1
```

and, for every coordinate `i`,

```text
(automaticRealizationMap P hdeg i).totalDegree ≤ 6 * P.natDegree + 2.
```

For every commutative test `K`-algebra `A`, it constructs an equivalence

```text
(AdjoinRoot P →ₐ[K] A) ≃
  GeneralGaugeJacobianOneFiberPoint ... A
```

where the right side is the literal fiber of the three actual
`MvPolynomial` coordinates at the chosen target. The equivalence commutes with
every algebra homomorphism `A →ₐ[K] B`.

The principal final declarations are:

```text
automaticRealizationMap_certificate
automaticJacobianOneFiberRepresentingEquiv
automaticJacobianOneFiberRepresentingEquiv_natural
```

No translation parameter, coefficient nonvanishing proof, chart unit, or
auxiliary abstract-fiber hypothesis remains as an external input.

## Formal chain

`GaugeAssembly.lean` and `GaugeInverseAssembly.lean` prove the low-degree and
complete finite high-degree coefficient identities over arbitrary commutative
rings. `GeneralGaugeMap.lean` packages them as one actual arbitrary-degree
three-variable polynomial map.

`GeneralGaugeJacobian.lean` expands the complete three-by-three Jacobian,
including all finite coefficient sums, and proves the determinant is exactly
`-2`. The fixed output scaling then has determinant `1`.

`GeneralGaugeInverse.lean` defines `G_π`, `β(π,S)`, and `E_{π,b,c}` and proves
the exact derivative factorization required by the represented-fiber theorem.
`GeneralGaugeDisplayedFiber.lean` proves that evaluating the actual `B` and
`C` coordinates gives precisely the marked equations.

`GeneralGaugeRawFiber.lean` starts from a literal triple satisfying the three
polynomial equations. From `t*q = π` with unit `π`, it constructs the source
chart unit over an arbitrary commutative test ring and proves that the literal
fiber is naturally represented by the generic inverse quotient.
`GeneralGaugeNormalization.lean` transports this theorem to the determinant-one
output normalization at the zero second target coordinate.

`GeneralGaugeRealization.lean` specializes to
`G(S) = P(a+S)-P(a)`, proves that the chosen inverse polynomial is exactly
`P(a+S)`, translates the quotient back to `K[T]/(P)`, and then removes the
supplied parameter through `chosenAdmissibleTranslation`.
`GeneralGaugeRealizationDegree.lean` transports both the Jacobian and degree
certificates to the final automatically chosen map in terms of the original
input degree.

## Foundations and axioms

The project contains no `sorry` and introduces no project-specific axioms.
The final functor-of-points theorem uses only Lean's standard foundations
`propext`, `Classical.choice`, and `Quot.sound`. The determinant and degree
certificates are algebraic theorems with no additional axioms.

The reconstruction works over arbitrary commutative test algebras. Units are
carried explicitly, so there is no hidden localization, omitted component, or
reduction to field-valued points. Separability makes the derivative class
invertible by Bézout; the unit first target coordinate makes the chart global
on the entire fiber functor.

## Remaining formal boundary

The actual map, determinant, effective degree, literal fiber, quotient
translation, and naturality are formalized. The remaining steps required for
a single Lean theorem matching the complete paper statement are:

1. prove the primitive linear-in-the-target-coordinate irreducibility theorem
   for the generic inverse and connect it to the function-field definition of
   geometric degree;
2. formalize monogenicity of arbitrary finite étale algebras over an infinite
   field and compose it with the polynomial-presentation theorem;
3. formalize, or explicitly isolate as a classical theorem interface, the
   Campbell--Razar--Wright degree-two Galois case;
4. formalize the local-number-theoretic and prime-distribution inputs in the
   Hasse-principle applications if those corollaries are to be machine-checked
   end to end.

The current certificate therefore proves the complete constructive and
scheme-theoretic polynomial-presentation layer, while keeping geometric degree,
monogenicity, the classical rank-two obstruction, and arithmetic inputs
explicitly separated.

## Build

```bash
cd formal/finite-etale-keller
lake build
```

Repository CI builds this project independently of the external Lean certificate
for the foundational three-dimensional map.
