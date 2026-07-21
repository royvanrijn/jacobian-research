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

## Cancellation programme

```bash
make verify-master
```

This runs the construction, parameter arithmetic, boundary, monodromy, and
current-ansatz rigidity regressions.

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

To retain an environment record and complete log under `artifacts/`, run:

```bash
make verify-logged
```

Generated outputs, bounded scans, and exploratory search programs are not
part of the public proof navigation.  Existing generated artifacts live under
`artifacts/generated-results/`; historical search tools are preserved under
`archive/tooling/`.
