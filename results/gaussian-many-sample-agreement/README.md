# Asymptotically Optimal Partitions for Agreement of Many Gaussian Samples

**Hexagon submission compiling.** Reserved `hexagon:2611.01282v1`. Updated November 28, 2026.

Roy van Rijn

Fix \(k\ge3\) and \(0<\rho<1\). When \(m\) standard Gaussian observations
have cross covariance \(\rho I_d\), the supremal probability of unanimous
agreement over balanced measurable \(k\)-cell partitions is \((2+o(1))H_m\)
as \(m\to\infty\), uniformly over \(d\ge1\). Here \(H_m\) is the all-samples
probability for one halfspace of mass \(1/k\), and
\(H_m=m^{-(1-\rho)/\rho+o(1)}\). The two extreme cells of a quantile slab
partition give the leading factor \(2\). The centered regular simplex, for
\(d\ge k-1\), loses by a power of \(m\). The paper also constructs exactly
balanced low-influence vote rules that beat plurality for sufficiently many
fixed correlated ballot copies.

[Review PDF](gaussian-many-sample-agreement.pdf) ·
[LaTeX source](main.tex) · [Source archive](source.zip) ·
[Metadata](metadata.json) · [Verification](verification.json) ·
[Artifact checksums](SHA256SUMS.txt)

The [canonical proof note](../../research/plurality-stability/MULTI_SAMPLE_SIMPLEX_COUNTEREXAMPLE.md)
gives the Gaussian and discrete arguments. The many-observation proposal
appears in [Section 8 of Isaksson–Mossel](https://arxiv.org/pdf/0903.3362v3#page=26).
Their numbered Conjecture 1.4 concerns two observations. The proof uses their
Exchangeable Gaussians Theorem and Prékopa’s three-function inequality. It
establishes asymptotic optimality of the quantile construction. Identification
of the exact finite-sample optimum and classification of all asymptotically
optimal partitions are separate questions.

## Reproduce the PDF

Extract `source.zip` into an empty directory. With TeX Live 2023/Debian and
latexmk 4.83, run from that directory.

```sh
SOURCE_DATE_EPOCH=1795824000 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=gaussian-many-sample-agreement main.tex
```

The source archive contains only the compilation entry point and its embedded
bibliography. The theorem uses no computational evidence.
