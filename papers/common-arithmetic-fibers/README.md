# Every Finite Étale Algebra of Rank at Least Three Is a Full Keller Fiber

This is the focused manuscript for the prescribed-fiber theorem.

Given a squarefree polynomial `P` of degree `N >= 3` over a
characteristic-zero field and a supplied admissible translation `a`, the
paper constructs an explicit map of affine three-space with:

- Jacobian determinant `1`;
- geometric degree `N`;
- coordinate degree at most `6N+2`; and
- a distinguished full fiber naturally isomorphic to `Spec K[T]/(P)`.

Monogenicity gives the abstract corollary: every finite étale algebra of rank
at least three occurs as a full Keller fiber. The concrete `(P,a)`
construction commutes with scalar extension. The primitive element and
translation selected from an abstract algebra are noncanonical; Lean makes
these choices noncomputably, and no functorial automatic choice is claimed.

The paper contains only:

1. the polynomial-presentation theorem;
2. the compact inverse-equation design;
3. root-to-source reconstruction and the literal fiber theorem;
4. the function-field comparison and fullness;
5. the finite-étale corollary via monogenicity;
6. base change and the Lean correspondence; and
7. the explicit Berend--Bilu quintic Hasse fiber.

The rank classification is stated only as a corollary. Its exclusion of rank
two uses the classical Campbell--Razar--Wright Galois case and is outside the
Lean certificate.

## Material moved to companion notes

- Symmetric monodromy is in
  [`verified/UNIVERSAL_SYMMETRIC_MONODROMY.md`](../../verified/UNIVERSAL_SYMMETRIC_MONODROMY.md).
- Stable compositional atomicity is in
  [`verified/PRIMITIVE_MONODROMY_ATOMICITY.md`](../../verified/PRIMITIVE_MONODROMY_ATOMICITY.md).
- Exact nonproperness, the complete `Pi = 0` fiber table, and the global
  discriminant order are in
  [`papers/quadratic-gauge-nonproperness`](../quadratic-gauge-nonproperness/).

These results are not used to prove fullness. Fullness follows directly from

```text
rank F^{-1}(y) = N = gdeg(F).
```

## Formal and independent verification

The Lean project in
[`formal/finite-etale-keller`](../../formal/finite-etale-keller/) proves the
polynomial construction end to end, including the literal natural fiber,
finite étaleness, rank, actual function-field comparison, geometric degree,
base change for supplied data, monogenicity, and the abstract finite-étale
corollary. It has no `sorry` and no project-specific axioms.

The public endpoints include:

- `automaticRealization_pageOne`;
- `generalGaugeSourceFunctionFieldComparison`;
- `generalGaugeGeometricDegree_eq`;
- `abstractFiniteEtale_pageOne`;
- `realizationMapTarget_map`; and
- `adjoinRootBaseChangeEquiv`.

The exact symbolic checks are:

```bash
.venv/bin/python scripts/verify_universal_quadratic_gauge.py
.venv/bin/python scripts/verify_root_engineered_quadratic_gauge.py
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
Singular -q scripts/verify_universal_quadratic_gauge.sing
```

Build Lean with:

```bash
cd formal/finite-etale-keller
lake build
```

Build the paper from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The directory name is retained as a stable path from the earlier draft.
