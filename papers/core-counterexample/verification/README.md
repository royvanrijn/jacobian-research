# Verification for the core paper

Run from the repository root:

```bash
make verify-foundations
```

This executes the independent determinant/collision checker, the main exact
counterexample check, both marked-root charts, image and exceptional-fiber
checks, and the weighted schema regressions.

The separately authored pinned Lean 4 verification is optional because it
downloads its own toolchain:

```bash
make verify-foundations-formal
```

Primary scripts:

- `scripts/verify_counterexample_independent.py`
- `scripts/verify_counterexample.py`
- `scripts/verify_marked_root_model.py`
- `scripts/verify_weighted_marked_root_model.py`
- `scripts/verify_weighted_seed_theorem.py`

Exact scripts certify displayed finite identities. Bounded weighted-seed
runs are regressions for the written all-degree proof, not substitutes for it.
