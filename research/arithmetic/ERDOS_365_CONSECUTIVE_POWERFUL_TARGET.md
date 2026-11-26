# Erdős #365: a subquarter consecutive-powerful count

Let \(N(X)=\#\{n\le X:n,n+1\text{ are powerful}\}\). Erdős #365 asks, in particular, for a polylogarithmic upper bound. [Reuss, Theorem 4](https://arxiv.org/pdf/1212.3150v2) gives \(N(X)\ll_\varepsilon X^{29/100+\varepsilon}\). [Plymesser, Theorem 1.1 and Remark 9.1](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf) claim \(X^{1/4+\varepsilon}\) and identify the binding dyadic corner: one squarefree cube kernel is bounded, the other has size \(X^{1/4}\), and the square variable on the latter side has size \(X^{1/8}\). The [existing source audit](PLYMESSER_2026_POWERFUL_COUNTS_SOURCE_SCOPE.md) records the literature and formalization boundary. Neither the new descent nor a subquarter exponent is currently proved here.

## Proposed fixed-fiber estimate

For fixed squarefree \(d\ge1\), fixed \(v\ge1\), and either sign, the proposed bound is
\[
 \#\{(u,e)\in\mathbb Z_{>0}^2:d^3u^2\pm1=v^2e^3\}
 \ll_\varepsilon d^{1/2+\varepsilon}v^\varepsilon. \tag{365.1}
\]
The implied constant must be uniform in \(d,v\). This is a **candidate**, not a theorem. The proposed route factors in the appropriate quadratic field, separates ideal classes, and reduces each branch to a binary cubic Thue equation. Audit the integral ideal factorization at primes over \(2dv\), the class-number and unit contributions, the real quadratic case and \(d=1\), and the precise hypotheses and coefficient uniformity of the imported Thue bound. In particular, bounding each cubic branch is insufficient unless the number of branches has the displayed dependence on \(d,v\).

## Global transfer still required

At the corner, summing (365.1) over \(v\ll X^{1/8}\) suggests an \(X^{1/8+\varepsilon}\) contribution for bounded \(d\). A subquarter theorem needs more: establish the estimate uniformly over a **neighborhood** of that corner, treat its mirror, and optimize every remaining dyadic region with the exact hypotheses of Plymesser's elementary, Reuss, Heath-Brown and reduced-component inputs. The source's quarter-exponent proof explicitly uses external literature premises; its internal formal checks are not an audit of those premises. A positive \(\delta\) and a certified global bound \(N(X)\ll_\varepsilon X^{1/4-\delta+\varepsilon}\) remain open. The polylogarithmic target is separate and stronger.

The closest source statement is Plymesser's Remark 9.1: it explicitly calls for a bound uniform in the bounded cube kernel and the opposite square variable on a family of genus-one curves. The fixed-fiber proposal is aimed at that gap, but no equivalence between (365.1) and the required global optimization has been established.
