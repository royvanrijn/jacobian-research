# Every Finite Étale Algebra Except Rank Two Is a Keller Fiber

This paper proves a complete classification: every nonzero finite étale
algebra of rank other than two occurs as a full fiber of a polynomial Keller
map. For rank at least three, the realization is explicit in affine
three-space, has Jacobian determinant `1`, and has coordinate degree at most
`6N+2` in rank `N`.

The scheme-theoretic core is a natural equivalence on every commutative test
algebra:

```text
Hom_K-alg(K[T]/(P), A) ≃ distinguished source-fiber points over A.
```

Thus the theorem controls the complete represented fiber, not merely its
geometric or rational points. Lean formalizes the existence and automatic
choice of an admissible translation, the abstract source-equation functor, its
two-sided reconstruction and naturality, representation by the polynomial
quotient, and translation back to `K[T]/(P)`. Its final theorem takes only a
squarefree polynomial of degree at least three; no translation parameter or
nonvanishing witness remains as an external hypothesis.

The displayed all-degree polynomial map is proved coefficientwise in the
paper. Lean checks the low-degree identities and the complete finite sums of
all high-degree terms `4 ≤ k ≤ N`, with an arbitrary coefficient family, over
arbitrary commutative rings. `GeneralGaugeMap.lean` now packages these sums as
one `MvPolynomial (Fin 3) K` map and certifies its three coordinate evaluations
over every commutative test algebra. The remaining map-level steps are its
general Jacobian theorem, `6N+2` `totalDegree` bound, and a direct equivalence
between the raw polynomial-map fiber and the represented-fiber datum.

Two independent exact checkers audit the construction. The structural checker
verifies the source and marked-line Jacobians, the generic `k`-th coefficient,
a six-coefficient bridge, and the termwise degree bound. The concrete checker
fully expands degrees three, four, and five, reconstructs the quotient fiber in
both directions, and audits the arithmetic examples. These are regression and
independence certificates, not replacements for the uniform paper and Lean
proofs.

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
2. two-sided reconstruction over arbitrary commutative test algebras, using a
   Bézout inverse of `E'` and proving naturality;
3. existence of an admissible translation and the canonical quotient
   translation `K[S]/(P(a+S)) ≃ K[T]/(P(T))`;
4. a coefficientwise derivation of the displayed all-degree gauge identities,
   including the complete finite sums and the general `MvPolynomial` object;
5. the scaling identity `F_displayed = diag(1,19,19) F_normalized` for the
   optimal quintic example;
6. the dated and qualified [literature audit](LITERATURE_AUDIT.md).

The [verification matrix](VERIFICATION.md) records the proof layer supporting
every load-bearing statement and the exact remaining formal boundary.

The Lean project contains no `sorry` and no project-specific axioms. Its final
automatic represented-fiber theorem reports only Lean's standard `propext`,
`Classical.choice`, and `Quot.sound`. The formal scope and the remaining
nonformalized parts of the paper are listed in
[`formal/finite-etale-keller/README.md`](../../formal/finite-etale-keller/README.md).

Run the exact checkers from the repository root:

```bash
.venv/bin/python scripts/verify_universal_quadratic_gauge.py
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
