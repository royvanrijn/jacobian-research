# Counterexamples to the Mathieu Property

*From Explicit Quotient Constructions to a General Classification*

**Released preprint · [Zenodo](https://zenodo.org/records/24310836) · [DOI: 10.5281/zenodo.24310836](https://doi.org/10.5281/zenodo.24310836) · Hexagon pending**

Classifies the closed subgroups $H$ of complex linear algebraic groups $G$, with both identity components reductive, for which the kernel of normalized Reynolds integration on $\mathbb C[G/H]$ is a Mathieu–Zhao subspace: exactly those containing $[G^0,G^0]$. When containment fails, the paper constructs one fixed multiplier detecting infinitely many nonzero marked moments.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Exact checker](verify.py) · [Verification](verification.json) · [Source audit](source_audit.json) · [Metadata](metadata.json) · [Deposited edition and checksums](published/)

The main result is Theorem 1.1. The introduction credits Long's earlier SU(2) counterexample, compact connected group classification and spherical-pair theorem, then explains the extra right invariance needed for finite quotients. Sections 2–4 construct the exceptional examples and transfer their coefficient calculation to cubic curves using Klein identities and Cartier compatibility. Section 5 gives an alternative radial pair and its closed moment formula; existence of radial transfer is already due to Long. Sections 6–7 reduce the general classification to these examples and treat disconnected groups and the moment nullcone. Section 8 gives a multiplier chosen from the recomposition degree alone and proves the Cauchy-transform and coefficient-algebra statements. The appendices supply the invariant-form, differential-extraction and polynomial-identity details.

The written proofs cover all orders. The exact checker separately verifies the polynomial identities and invariant substitutions, directly integrates exceptional moments at orders 1 and 2, checks the binary coefficients at 28 primes, and compares radial moments through order 8 using two separately implemented integration methods. It also checks the double-integral evaluation and rejects two deliberately corrupted controls. The checker avoids expanding the full high-degree witness by integrating its weight-zero component. These finite controls do not certify the classification or replace the all-order proofs. Internal mathematical, citation and layout review has been completed; independent replay, formal verification and external review are not claimed. The source audit retains its earlier review snapshot; the verification record pins the final deposited source and PDF and records the internal replay after this move.

From the repository root, replay the checks with the existing Python environment. The recorded run uses Python 3.12.3 and SymPy 1.14.0; the SymPy version is pinned in [the research requirements](../../research/requirements.txt).

```sh
.venv/bin/python results/mathieu-reductive-homogeneous-spaces/verify.py
```

From this paper directory, build the PDF with TeX Live and `latexmk`:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

The command produces the current `paper.pdf`. The exact deposited PDF and `paper.tex` remain in `published/`; both currently match the top-level PDF and `main.tex` byte for byte. The verification record identifies the checked source, checker and PDF by SHA-256 and states the scope of each check.

Next: Obtain an external mathematical review and an independent replay, resolve any findings, and track the pending Hexagon submission. The Zenodo deposit is public; independent review and the scope of verification remain separate from publication.
