# Sparse Minimality of Gaussian-Moments Counterexamples in Dimension Three

Standalone source for the paper by Roy van Rijn.

Build from this directory with:

```sh
latexmk -pdf main.tex
```

The paper proves by exhaustive exact computation that no Gaussian-Moments
counterexample in a Gaussian Witt frame has both degree at most four and at
most four monomials.  Replay the coefficient-scheme calculation from the
repository root with:

```sh
.venv/bin/python scripts/verify_three_real_gmc_sparse_pareto.py
```

The script first computes the charge/parity contraction Hilbert basis,
quotients by overall coefficient scaling and the Gaussian-preserving
\(W/Z\) torus whenever the support contains distinct charges, and then
computes the remaining rational moment ideals directly in the Laurent
coordinate ring of the coefficient torus.
