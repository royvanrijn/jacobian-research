# Reproducing the results

The Makefile is the public verification interface.  Run commands from the
repository root after creating the Python environment described in the main
[README](README.md).

## Fast structural check

```bash
make check
```

This compiles the active Python code, checks local Markdown links, and audits
the single status ledger.

## Stable core

```bash
make verify-minimal
make verify-core
make verify-foundations
```

`verify-minimal` uses only the Python standard library for the foundational
map.  `verify-core` adds the cubic marked-root and exact-image implementations.
`verify-foundations` adds the weighted construction and its clean-room checker.

The separately authored Lean certificate is optional because it downloads a
pinned toolchain:

```bash
make verify-lean-foundational
```

GitHub Actions runs this target in the required `formal-lean` job using the
pinned upstream commit and Lean action.  The `papers` job compiles all three
maintained papers, and `macaulay2-independent-check` runs the pinned
Macaulay2 comparison.  Together with the four Python matrix jobs, these are
the complete CI verification pipeline.  The final `verification-complete`
job is the single aggregation check intended for GitHub branch protection.

## Cancellation programme

```bash
make verify-master
```

This runs the construction, parameter arithmetic, boundary, monodromy, and
current-ansatz rigidity regressions.

## External quartic islands

Juntang Zhuang's pinned `F4a`, `F4b`, and `F4c` examples have an independent
compact reconstruction and canonical-boundary audit:

```bash
.venv/bin/python scripts/verify_external_quartic_islands.py
```

This command is also part of `make verify-regressions`.  It requires no
network access and does not copy or execute the upstream checker.

## External consequence identities

Christopher D. Long's direct Gaussian-moment and `(xz)` identities, together
with the exact normalization of the foundational map used in his BCW
discussion, have a dedicated target:

```bash
make verify-external-consequences
```

The Gaussian and `(xz)` scripts use only the Python standard library.  Their
bounded exact regressions are distinguished from Long's written all-exponent
proofs.  A companion symbolic checker proves the `SU(2)=S^3` Haar density in
Hopf coordinates, completing the local integration proof.  The same target
also performs all 18 balanced BCW steps and checks the resulting 79-variable
cubic-homogeneous collision, writes its sparse artifact, and replays it with a
separate standard-library implementation.  A local proof of the
fixed-dimensional DVEZ/Zhao implication, including Gaussian contraction, the
countable-union step, and formal inversion, completes the nonexplicit route to
`not GMC(158)`.  It also verifies the uniform weighted-seed Gaussian bridge:
the exact pencil branch, polynomial determinant correction, and bounded Wick
moments for canonical and split seeds, followed by a separate standard-library
reconstruction.  These checks are part of
`verify-regressions`, not `verify-minimal`.

The same target runs the rank-two Poisson **pre-audit**.  It verifies that the
single displayed output `R=x(2-3xq)` is exactly the foundational third output
after a polynomial source automorphism, and proves that the naive choices
`S=F_1/2`, `T=F_2` have no polynomial `D`-completion.  It does not reconstruct
or certify the unavailable manuscript formulas.

The generated certificate is stored as
[`artifacts/generated-results/long_bcw_79_counterexample.json`](artifacts/generated-results/long_bcw_79_counterexample.json).
It records the sparse cubic map, all reduction-step choices, and the three
exact collision points; regeneration is deterministic.

## Complete active suite

```bash
make verify
```

## Canonical degreewise paper

The canonical statement and proof are in the standalone paper; the five-lemma
verification companion is
[`DEGREEWISE_MULTIPLICITY_AUDIT.md`](DEGREEWISE_MULTIPLICITY_AUDIT.md).  Build
the paper with:

```bash
cd papers/marked-root-multiplicity
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

To compile every maintained paper with the same target used by CI, run:

```bash
make verify-papers
```

To retain an environment record and complete log under `artifacts/`, run:

```bash
make verify-logged
```

Generated outputs, bounded scans, and exploratory search programs are not
part of the public proof navigation.  Existing generated artifacts live under
`artifacts/generated-results/`; historical search tools are preserved under
`archive/tooling/`.
