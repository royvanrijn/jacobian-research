# The Gaussian Moments Conjecture:
# Dimension Two and Sparse Minimality in Dimension Three

Standalone source for the paper by Roy van Rijn.

Build from this directory with:

```sh
latexmk -pdf main.tex
```

The proof has four dependencies:

1. the circular-coordinate identity
   \(\mathbb E(F)=\mathcal L(\operatorname{CT}_T F)\);
2. an elementary supporting-line argument for the lowest radial orders;
3. the Duistermaat--van der Kallen constant-term theorem;
4. Frobenius and factorial isolation after reduction modulo a prime.

The paper now has two main results:

1. a conceptual proof of GMC(2), requiring no support computation; and
2. an exhaustive exact sparse Pareto theorem in dimension three: no GMC
   failure has both degree at most four and at most four monomials in a
   Gaussian Witt frame.

Replay the computer-assisted second theorem from the repository root with:

```sh
.venv/bin/python scripts/verify_three_real_gmc_sparse_pareto.py
```

The script recomputes all rational saturated moment ideals and takes several
minutes.

A companion Lean 4 package is located at
[`../../formalization/gmc2`](../../formalization/gmc2/README.md).  It checks
the factorial-divisibility, prime-isolation, Frobenius/constant-term, and
eventual one-sided-support modules without `sorry`.  The DvdK theorem,
rational supporting-face extraction, and finite-type specialization remain
explicit named axioms.
