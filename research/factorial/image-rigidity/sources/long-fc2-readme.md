# Factorial Conjecture in Two Variables

Christopher D. Long

[Paper PDF](factorial_conjecture_two_variables.pdf) · [LaTeX source](paper/factorial_conjecture_two_variables.tex) · [Verification scripts](scripts/) · [Revision record](audits/2026-10-07-exposition.md)

This repository contains an unpublished research manuscript on the factorial functional

$$
\mathcal{L}(x^a y^b)=a!b!.
$$

The **nonflat theorem** proves that a nonzero polynomial $f$ of total degree $D$ has infinitely many nonzero factorial moments whenever its highest homogeneous part $f_D$ is not a scalar multiple of $(x+y)^D$. Equivalently, $A(t)=f_D(t,1-t)$ is nonconstant. This argument uses an exact radial Mittag–Leffler transform, corrected leaf identities with retained continuation histories, polynomial monodromy, simultaneous block traces, and Lindemann–Weierstrass. It is independent of EIT I.

Combining the nonflat theorem with the explicitly stated single-integral nonvanishing theorem of [EIT I, Corollary 1.4 (`cor:single`), at the pinned commit](https://github.com/long-mathematics/exponential-integral-theorem/blob/22d98c4f27e606b5b25fa673c33f9c56c53cc86b/papers/01-exponential-integral-theorem/exponential-integral-theorem.tex) gives the **full two-variable Factorial Conjecture**:

$$
\left(\mathcal{L}(f^n)=0\text{ for every }n\ge1\right)\Longrightarrow f=0.
$$

The remaining flat case has $f_D=c(x+y)^D$. After normalization and algebraic specialization, its proof reduces to the nonvanishing of $\int_0^1\exp(A_1(t)/D)\,dt$. When $\deg A_1\le1$, including total degree $D\le2$, this last step is elementary apart from Hermite–Lindemann and does not require EIT I.

As a formal consequence, the manuscript also settles the two-dimensional product-exponential case of van den Essen's Integral Conjecture. For

$$
M=\lbrace h\in\mathbb{C}[x,y]:
\int_0^\infty\!\int_0^\infty h(x,y)e^{-x-y}\,dx\,dy=0\rbrace
=\ker\mathcal{L},
$$

the paper proves that $M$ is a Mathieu–Zhao subspace; more strongly,

$$
\lbrace f\in\mathbb{C}[x,y]:\mathcal{L}(f^m)=0\text{ for all sufficiently large }m\rbrace=\lbrace0\rbrace.
$$

The FC–Mathieu connection itself is not claimed as new: the [Integral Conjecture](https://arxiv.org/abs/1006.5801v1) and the [Factorial Conjecture](https://arxiv.org/abs/1008.3962v2) arise from the same Image-Conjecture program. The new input to this consequence is the full inhomogeneous FC(2) result proved here.

## Status and dependencies

The manuscript is **unpublished and unrefereed**. Its nonflat argument uses published Lindemann–Weierstrass and Pakovich–Muzychuk results; its general flat conclusion additionally uses a separate research theorem from EIT I. Matching the cited EIT statement is not an independent verification of its proof.

The manuscript proves the endpoint averaging, planar separation, path-module containment, analytic leaf transport, block traces, and arithmetic pairing needed in its own notation. All Pakovich–Muzychuk numbering refers to arXiv:0710.4085v2. Exact source versions, pinpoints, and roles are recorded in [the reference ledger](references/source_references.json).

Supporting computations check finite identities, not FC(2), analytic continuation, polynomial monodromy, or EIT I as universal theorems. Model-assisted audit records are not external refereeing or Lean formalization.

Appendices A and C illustrate continuation conventions and a collision of critical points with incomparable decompositions. Appendix B is retained as an **auxiliary density-certificate appendix**, not a bound for the first nonzero factorial moment. It records Moura's estimate, a formal-power-series proof, and the Bostan–Dumas Wronskian argument. The October 2026 `openai/math` release influenced exploration and supplies an auxiliary proof/formalization reference; none of its results is an input to the main FC(2) proof.

## Build, checks, and PDF snapshot

Install TeX Live with `latexmk` and the standard LaTeX/font packages used by the manuscript, and install the pinned Python dependency:

```sh
python3 -m pip install -r requirements.txt
make
make test
```

`make` builds `output/pdf/factorial_conjecture_two_variables.pdf` and checks the final LaTeX log for unresolved citations/references, duplicate labels, missing glyphs, and fatal diagnostics. The test suite includes the original 17 degree-six checks and 405 symbolic checks, eight supplementary groups (including 70,215 exact radial-coefficient tuples), and 11 release-tool regression tests. The tuple grid is a finite sanity check, not 70,215 independent mathematical proofs.

After changing the manuscript or its build inputs, refresh the browser-facing PDF:

```sh
make snapshot
make check
```

Commit `factorial_conjecture_two_variables.pdf` and `verification/pdf_snapshot.json` together with the source changes. The manifest fingerprints the source, Makefile, snapshot script, and committed PDF. `make snapshot-check` detects stale sources or a modified PDF. It deliberately does not compare a newly built PDF byte-for-byte with the committed one: compiler versions and metadata can change those bytes.

CI runs `make check` with read-only repository permissions and uploads the built PDF, final log, and generated exact-check reports. It does not automatically commit to `main`. Visual inspection remains part of a release; the automated checks do not establish that a proof is correct or a page layout is acceptable.

## Repository layout

- `paper/` and the root PDF contain the manuscript source and tracked snapshot.
- `scripts/` contains the finite mathematical checks and build/snapshot validators.
- `verification/` contains the tracked snapshot manifest and documentation for generated reports.
- `references/` and `audits/` record source provenance and the scope of revision checks.
- `.github/workflows/latex.yml`, `Makefile`, and `requirements.txt` define CI and the local build.

`CITATION.cff` supplies citation metadata. The manuscript and mathematical exposition are licensed under **CC BY 4.0**; code and build tooling are under the **MIT License**. See [LICENSE.md](LICENSE.md).
