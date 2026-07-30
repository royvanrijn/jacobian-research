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

A contemporaneous manuscript of Christopher D. Long,
[*A Factorially Weighted Constant-Term Theorem on Algebraic
Tori*](https://github.com/octonion/mathematics/blob/main/gmc/gmc2_stronger_arbitrary_torus.tex),
proves a stronger one-radial theorem in arbitrary torus rank.  Its rank-one
specialization overlaps the theorem in this paper.  The paper credits that
overlap explicitly; the repository's local comparison and proof audit are in
[`../../extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md`](../../extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md).

A companion Lean 4 package is located at
[`../../formal/gmc2`](../../formal/gmc2/README.md).  It checks
the full bivariate theorem, from the circular substitution and Wick formula
through lower-face extraction, finite coefficient-ring descent,
prime-isolation, and eventual one-sided support, without `sorry`.  The
rational supporting-face extraction is proved directly; the DvdK theorem
and finite-type specialization remain explicit named axioms.
