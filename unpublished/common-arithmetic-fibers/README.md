# Over Characteristic Zero, Every Finite Etale Algebra of Rank at Least Three Is a Full Keller Fiber

**Unpublished · draft.**

Gives a quadratic-gauge proof and Lean formalization of prescribed full finite etale Keller fibres. The characteristic-zero existence conclusion also follows from Ulam Theorem 5.2; the explicit coefficient specialization is included in the paper and [literature comparison](LITERATURE_AUDIT.md).

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) · [Source bundle](source.zip) · [Current edition checks](verification.json)

The characteristic-zero polynomial-presentation theorem and abstract finite-etale corollary have a Lean development; broader supplied-presentation scope has a separate assurance boundary.

The [formal source supplement](https://github.com/royvanrijn/jacobian-research/tree/c32b9b4d273cf7988e85184a0d5eae9c488727fd/research/foundations/formal/finite-etale-keller)
is pinned to repository commit `c32b9b4d273cf7988e85184a0d5eae9c488727fd`.
[VERIFICATION.md](VERIFICATION.md) states the exact build and axiom-audit
instructions and distinguishes the existing development from a fresh replay.

Run the exact coefficient comparison from this directory with
`python3 verify_comparison.py`; it does not verify the full reconstruction.

Next: Obtain external review of the written/formal correspondence and complete
specialist novelty review of the gauge and supplied-characteristic refinements.

The manuscript source bundle builds with
`latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex`.
