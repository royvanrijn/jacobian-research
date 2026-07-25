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
paper. Lean now checks the generic high-degree monomial transport, the
low-degree second-coordinate identity, and the cubic cancellation producing
the third coordinate over arbitrary commutative rings. The remaining finite-sum
packaging and `6N+2` degree estimate are paper proofs. The exact symbolic
checker independently constructs and audits the complete map in degrees three,
four, and five, together with the explicit arithmetic examples; these are
regression certificates rather than a replacement for the uniform proof.

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
   Galois exclusion, with Campbell--Razar--Wright provenance and an explicit
   descent proof;
2. two-sided reconstruction over arbitrary commutative test algebras, using a
   Bézout inverse of `E'` and proving naturality;
3. existence of an admissible translation and the canonical quotient
   translation `K[S]/(P(a+S)) ≃ K[T]/(P(T))`;
4. a coefficientwise derivation of the displayed all-degree gauge identities,
   with its monomial and cubic cancellations formalized in Lean;
5. the scaling identity `F_displayed = diag(1,19,19) F_normalized` for the
   optimal quintic example;
6. the dated and qualified [literature audit](LITERATURE_AUDIT.md).

The Lean project contains no `sorry` and no project-specific axioms. Its final
automatic represented-fiber theorem reports only Lean's standard `propext`,
`Classical.choice`, and `Quot.sound`. The formal scope and the remaining
nonformalized parts of the paper are listed in
[`formal/finite-etale-keller/README.md`](../../formal/finite-etale-keller/README.md).

Run the exact checker from the repository root:

```bash
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
