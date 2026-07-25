# Every Nonzero Finite Étale Algebra Except Rank Two Is a Keller Fiber

This paper proves a complete classification: every nonzero finite étale
algebra of rank other than two occurs as a full fiber of a polynomial Keller
map. For rank at least three, the realization is explicit in affine
three-space, has Jacobian determinant `1`, and has coordinate degree at most
`6N+2` in rank `N`.

The scheme-theoretic core is a natural equivalence on every commutative test
algebra:

```text
Hom_K-alg(K[T]/(P), A) ≃ literal distinguished map-fiber points over A.
```

Thus the theorem controls the complete represented fiber, not merely its
geometric or rational points. Starting from a squarefree polynomial `P` of
degree at least three, Lean now:

- chooses an admissible translation parameter internally;
- constructs the actual arbitrary-degree `MvPolynomial (Fin 3) K` map;
- proves its general Jacobian is `1` after the fixed output normalization;
- proves every coordinate has total degree at most `6N+2`;
- constructs the literal target fiber of those three polynomial coordinates;
- identifies that fiber naturally with maps from `K[T]/(P)` over every
  commutative test algebra.

The final formal declarations are
`automaticRealizationMap_certificate` and
`automaticJacobianOneFiberRepresentingEquiv_natural`. No translation
parameter, coefficient nonvanishing proof, chart unit, abstract source-fiber
wrapper, or bounded-degree specialization remains as an external input to the
polynomial-presentation theorem.

Three independent exact layers audit the construction:

1. Lean proves the uniform finite sums, actual map, general determinant,
   effective degree bound, literal fiber equivalence, quotient translation,
   and naturality.
2. A structural SymPy checker verifies the source and marked-line Jacobians,
   the generic `k`-th coefficient identities, a six-coefficient bridge, and
   the termwise degree bound.
3. Singular expands the full degree-six map over the rational function field
   in six algebraically independent coefficients, obtains determinant `-2`,
   and checks the `(7,38,36)` degree profile. The concrete checker independently
   expands degrees three, four, and five and reconstructs their quotient fibers
   in both directions.

The arithmetic applications include:

- an explicit degree-five Keller fiber that is everywhere locally soluble
  over `Q` but has no rational point, with degree five proved optimal;
- one fixed Keller map with infinitely many such Hasse-failing fibers;
- exact transfer of connectedness, signatures, splitting fields, local
  factorization data, and intersectivity;
- compatibility with extension of the ground field.

The directory name is retained as a stable repository path from the earlier
draft.

The focused audits accompanying the active draft are:

1. the exact arbitrary-characteristic-zero-field scope of the degree-two
   Galois exclusion, with Campbell--Razar--Wright provenance and a faithfully
   flat descent proof;
2. the primitive linear-in-the-target-coordinate proof of generic inverse
   irreducibility and geometric degree;
3. two-sided reconstruction over arbitrary commutative test algebras, using a
   Bézout inverse of `E'` and proving naturality;
4. existence of an admissible translation and the canonical quotient
   translation `K[S]/(P(a+S)) ≃ K[T]/(P(T))`;
5. a coefficientwise and map-level derivation of the displayed all-degree
   gauge identities, Jacobian, and effective degree bound;
6. the scaling identity `F_displayed = diag(1,19,19) F_normalized` for the
   optimal quintic example;
7. the dated and qualified [literature audit](LITERATURE_AUDIT.md).

The [verification matrix](VERIFICATION.md) records the proof layer supporting
every load-bearing statement and the exact remaining formal boundary.

The Lean project contains no `sorry` and no project-specific axioms. Its final
literal-fiber theorem reports only Lean's standard `propext`,
`Classical.choice`, and `Quot.sound`; the Jacobian and degree certificates use
no additional axioms. The formal scope and remaining nonformalized inputs are
listed in
[`formal/finite-etale-keller/README.md`](../../formal/finite-etale-keller/README.md).

Run the exact independent checkers from the repository root:

```bash
.venv/bin/python scripts/verify_universal_quadratic_gauge.py
Singular -q scripts/verify_universal_quadratic_gauge.sing
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
```

Build the Lean certificate with:

```bash
cd formal/finite-etale-keller
lake build
```

Build the paper with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=tmp/pdfs/common-arithmetic-fibers \
  papers/common-arithmetic-fibers/main.tex
```

The active-paper build copies the PDF to:

```text
output/pdf/common-arithmetic-fibers.pdf
```
