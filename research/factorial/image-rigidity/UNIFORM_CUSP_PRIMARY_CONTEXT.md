# Classical context for the uniform cusp proof

The mathematical proof does not import a finite-Gamma Lee--Yang circle law.

The complete cached primary was read: Zakhar Kabluchko, *Lee--Yang zeroes of
the Curie--Weiss ferromagnet, unitary Hermite polynomials, and the backward
heat flow*, arXiv:2203.05533v3, 30 December 2024,
https://arxiv.org/pdf/2203.05533v3 .

Actual inspected locators are equations (1.5), (1.8)--(1.9), Lemma 2.1
(printed p.5), and the critical parameter discussion attached to Theorem
2.18. In his convention

    H_m(z;t)=sum_j (-1)^(m-j) binom(m,j)
                           exp[-t*j*(m-j)/2] z^j.

Therefore the exact Gaussian model used here is

    sum_k binom(m,k) q^(m-k) exp[(2k-m)^2/(2m)]
           =(-1)^m exp(m/2) H_m(-q;4/m).

It is the finite critical Curie--Weiss model at beta=1. Lemma 2.1 gives the
classical circle law for H_m; the present Gamma argument does not transfer
that exact circle law to a finite Gamma polynomial.

The earlier scratch GAMMA_CURIE_WEISS_CRITICAL_BRIDGE.md, SHA256
3702d95a6295da1d0b952eaa6f95c283cf8c80746dd6db836556f58d9ae28698,
already identified the model and proved the genuine O(m/d) obstruction to
relative coefficient transfer. The canonical 402 compact normalized family
and bounded-H theorem already reach the same critical phase. There is no
claim that this correspondence or the Belyi/Curie--Weiss family is new.

The contribution of the new regular-complement argument is uniform transfer
of strict complete-contour geometry using the ACTUAL Gamma critical values,
with no m/s smallness assumption. The primary PDF is retained under
uniform_cusp_sources/ and is context only, not a mathematical premise for the
horizontal maximum calculation or for the actual Gamma theorem.

