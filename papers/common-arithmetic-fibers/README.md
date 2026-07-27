# Prescribed Finite Étale Algebras as Full Fibers of Keller Maps with Symmetric Monodromy

The paper's main theorem realizes every finite étale algebra of rank at least
three as a full fiber of a polynomial Keller map. The realization is explicit
in affine three-space, has Jacobian determinant `1`, and has coordinate degree
at most `6N+2` in rank `N`. As a corollary, the attainable nonzero ranks are
exactly `1,3,4,5,...`; rank two is excluded by the known degree-two theorem.
For every admissible seed, the generic inverse polynomial has geometric and
arithmetic Galois group `S_N`; this is not merely a generic-coefficient
statement. Thus every prescribed algebra is placed on a maximally symmetric
generic cover. Over a Hilbertian field, each fixed constructed map has
infinitely many connected full fibers with splitting-field group `S_N`.

The scheme-theoretic core is organized in three layers.

First, the general localized quadratic-gauge fiber theorem identifies an
arbitrary fiber with `(K[S]/(E))[1/E']`. Its finite-étale (squarefree)
corollary removes the localization. The prescribed-algebra realization theorem
then specializes this quotient to the chosen algebra `K[T]/(P)`.

At the prescribed target this gives a natural equivalence on every commutative
test algebra:

```text
Hom_K-alg(K[T]/(P), A) ≃ literal distinguished map-fiber points over A.
```

Thus the prescribed-algebra realization theorem controls the complete
represented fiber, not merely its geometric or rational points. Starting from
a squarefree polynomial `P` of degree at least three, Lean now:

- chooses an admissible translation parameter internally;
- constructs the actual arbitrary-degree `MvPolynomial (Fin 3) K` map;
- proves its general Jacobian is `1` after the fixed output normalization;
- proves every coordinate has total degree at most `6N+2`;
- constructs the literal target fiber of those three polynomial coordinates;
- identifies that fiber naturally with maps from `K[T]/(P)` over every
  commutative test algebra;
- proves that the representing quotient has dimension `N`;
- proves that the representing quotient is finite étale;
- promotes `Π` and `B` to independent parameters and proves that the fully
  generic inverse equation is irreducible of degree `N` over the iterated
  target field `K(Π,B)(C)`, with root quotient of finrank `N` and
  specialization back to the displayed inverse polynomial;
- proves that the three displayed coordinates are algebraically independent
  and that their coordinate substitution extends to an injective pullback on
  rational function fields;
- identifies the actual source function field over `K(Π,B)(C)` with the
  inverse-root quotient, and transfers its finrank to prove geometric degree
  `N`;
- verifies that the determinant-one normalization is an invertible target
  rescaling and packages determinant, geometric degree, literal fiber,
  naturality, finite étaleness, rank, and the degree bound in one theorem;
- proves that the complete supplied-translation map and distinguished target
  commute coefficientwise with extension of the ground field, and proves the
  corresponding tensor-product base change of the representing quotient; and
- specializes the construction to the paper's exact denominator-free quintic
  map, proving that its literal fiber at `(1,0,-38)` is naturally represented
  by `Q[T]/((T^3-19)(T^2+T+1))`, has rank five, has no rational point, and has
  both a real point and a three-adic point.

The combined final declaration is `automaticRealization_pageOne`, with
`automaticRealizationGeometricDegree_eq` supplying its geometric-degree
field.  The explicit function-field bridge is
`generalGaugeSourceFunctionFieldComparison`, and
`generalGaugeGeometricDegree_eq` transfers the inverse-root degree to the
actual displayed map.  The preceding layers remain separately exposed by
`automaticRealizationMap_certificate`,
`automaticJacobianOneFiberRepresentingEquiv_natural`,
`automaticRepresentingAlgebra_etale`,
`automaticRepresentingAlgebra_finite`,
`generalGaugeFullyGenericInversePolynomial_certificate`,
`generalGaugeFullyGenericInverseAdjoinRoot_finrank`,
`generalGaugeMap_algebraicIndependent`, and
`generalGaugeFunctionFieldHom_injective`.  Supplied-parameter scalar extension
is certified by `realizationMapTarget_map` and
`adjoinRootBaseChangeEquiv`. The concrete quintic fiber
declarations are `integralFiberRepresentingEquiv_natural`,
`p5_quotient_etale`, `p5_quotient_finite`, and `p5_quotient_rank`; its
arithmetic declarations are
`integralFiberPoint_rat_isEmpty` and
`integralFiberPoint_real_nonempty`, while
`integralFiberPoint_threeAdic_nonempty` is the first formal nonarchimedean
local-point certificate. No translation
parameter, coefficient nonvanishing proof, chart unit, abstract source-fiber
wrapper, or bounded-degree specialization remains as an external input to the
polynomial-presentation theorem.

Three independent exact layers audit the construction:

1. Lean proves the uniform finite sums, actual map, general determinant,
   effective degree bound, literal fiber equivalence, finite étaleness,
   quotient translation, naturality, fully independent inverse irreducibility
   and degree over `K(Π,B)(C)`, coordinate algebraic independence, the
   injective function-field pullback, the explicit source/inverse-root
   comparison, actual geometric degree, and the explicit quintic's rational
   obstruction, real point, and three-adic point.
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
- exact transfer of connectedness, signatures, splitting fields, local
  factorization data, and intersectivity;
- compatibility with extension of the ground field.

The active paper is deliberately narrower than the surrounding repository,
but it now includes the exact reduced nonproperness theorem, the complete
`Pi = 0` fiber table, and the boundary-sheet ledger in
`sections/02b-nonproperness.tex`.  These are ordinary mathematical proofs,
not Lean theorems.  The theorem defines the Jelonek locus as a reduced
graph-boundary/non-finite locus, uses Jelonek's complex multiplicity criterion
only after a characteristic-zero descent and base-change argument, and proves
the exact global `Pi`-adic discriminant factor in `K[B,C][Pi]`.

The final addendum retains one direction without theorem status: a candidate
fixed quintic gauge with infinitely many Hasse-failing fibers.  Its working
draft remains in `sections/05-infinite-family.tex` and is not included by
`main.tex`.

The directory name is retained as a stable repository path from the earlier
draft.

The focused audits accompanying the active draft are:

1. the exact arbitrary-characteristic-zero-field scope of the degree-two
   Galois exclusion, with Campbell--Razar--Wright provenance and a faithfully
   flat descent proof;
2. the primitive linear-in-the-target-coordinate proof of generic inverse
   irreducibility, followed by the explicit comparison
   `K(x,y,z) = K(Π,B,C)(S)` proving the actual geometric degree; Lean
   independently formalizes the same comparison over the canonical iterated
   target presentation;
3. the birational quadratic-discriminant parametrization and Morse-polynomial
   proof of full symmetric generic monodromy, with separate specialization and
   regularity arguments;
4. two-sided reconstruction over arbitrary commutative test algebras, using a
   Bézout inverse of `E'` and proving naturality;
5. existence of an admissible translation and the canonical quotient
   translation `K[S]/(P(a+S)) ≃ K[T]/(P(T))`;
6. a coefficientwise and map-level derivation of the displayed all-degree
   gauge identities, Jacobian, and effective degree bound;
7. the scaling identity `F_displayed = diag(1,19,19) F_normalized` for the
   optimal quintic example;
8. the dated and qualified [literature audit](LITERATURE_AUDIT.md).

The [verification matrix](VERIFICATION.md) records the proof layer supporting
every load-bearing statement and the exact remaining formal boundary.

The Lean project contains no `sorry` and no project-specific axioms. Its final
page-one theorem and function-field comparison report only Lean's standard
`propext`, `Classical.choice`, and `Quot.sound`. The formal scope and remaining
nonformalized inputs are listed in
[`formal/finite-etale-keller/README.md`](../../formal/finite-etale-keller/README.md).

Run the exact independent checkers from the repository root:

```bash
.venv/bin/python scripts/verify_universal_quadratic_gauge.py
.venv/bin/python scripts/verify_root_engineered_quadratic_gauge.py
Singular -q scripts/verify_universal_quadratic_gauge.sing
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
```

The independent nonproperness and boundary audit can be run with:

```bash
.venv/bin/python scripts/verify_quadratic_gauge_nonproperness.py
```

Build the Lean certificate with:

```bash
cd formal/finite-etale-keller
lake build
```

Build the paper with:

```bash
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=../../tmp/pdfs/common-arithmetic-fibers \
  papers/common-arithmetic-fibers/main.tex
```

The active-paper build copies the PDF to:

```text
output/pdf/common-arithmetic-fibers.pdf
```
