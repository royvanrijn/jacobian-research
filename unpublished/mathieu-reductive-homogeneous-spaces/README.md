# Counterexamples to the Mathieu Property

*From Explicit Quotient Constructions to a General Classification*

**Unpublished · ready for review.**

Classifies the closed subgroups $H$ of complex linear algebraic groups $G$, with both identity components reductive, for which the kernel of normalized Reynolds integration on $\mathbb C[G/H]$ is a Mathieu–Zhao subspace: exactly those containing $[G^0,G^0]$. When containment fails, the paper constructs one fixed multiplier detecting infinitely many nonzero marked moments.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Exact checker](verify.py) · [Verification](verification.json) · [Source audit](source_audit.json) · [Metadata](metadata.json)

The main result is Theorem 1.1. The introduction explains the relation to Mathieu's original conjecture and Long's spherical-pair result, then gives the complete reduction strategy. Sections 2–4 prove the exceptional finite-subgroup examples from one binary-form construction. Section 5 starts with a worked first moment and proves the radial identities by phase cancellation and a positive double integral. Sections 6–7 reduce the general classification to these examples and treat disconnected groups and the moment nullcone. Section 8 proves the recomposition, Cauchy-transform and coefficient-algebra statements. The appendices supply the invariant-form, differential-extraction and polynomial-identity details.

The written proofs cover all orders. The exact checker separately verifies the polynomial identities and invariant substitutions, directly integrates exceptional moments at orders 1 and 2, checks the binary coefficients at 28 primes, and compares radial moments through order 8 using two separately implemented integration methods. It also checks the double-integral evaluation and rejects two deliberately corrupted controls. The checker avoids expanding the full high-degree witness by integrating its weight-zero component. These finite controls do not certify the classification or replace the all-order proofs. Internal mathematical, citation and layout review has been completed; independent replay, formal verification and external review are not claimed.

From the repository root, replay the checks with the existing Python environment. The recorded run uses Python 3.12.3 and SymPy 1.14.0; the SymPy version is pinned in [the research requirements](../../research/requirements.txt).

```sh
.venv/bin/python unpublished/mathieu-reductive-homogeneous-spaces/verify.py
```

From this paper directory, build the PDF with TeX Live and `latexmk`:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

The command produces `paper.pdf`; it does not overwrite any published edition. The verification record identifies the checked source, checker and PDF by SHA-256 and states the scope of each check.

Next: Obtain an external mathematical review and an independent replay, resolve any findings, and finalize submission metadata before depositing. This remains an unpublished research draft.
