# The Gaussian Moments Conjecture in Two Variables

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
4. the filtered Frobenius isolation lemma after reduction modulo a prime.

A companion Lean 4 package is located at
[`../../formalization/gmc2`](../../formalization/gmc2/README.md).  It checks
the full bivariate theorem, from the circular substitution and Wick formula
through lower-face extraction, finite coefficient-ring descent,
prime-isolation, and eventual one-sided support, without `sorry`.  The
rational supporting-face extraction is proved directly; the DvdK theorem
and finite-type specialization remain explicit named axioms.
