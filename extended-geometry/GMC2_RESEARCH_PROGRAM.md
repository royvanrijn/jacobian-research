# GMC(2) research program

## Scope

The Gaussian Moments Conjecture is settled negatively in every real Gaussian
dimension \(n\geq 3\): Long's explicit three-variable polynomial

\[
 P=(1+Z)\left(W-\frac{2+Z}{2}T^2\right)
\]

has the five-term expansion

\[
 P=W+WZ-T^2-\frac32ZT^2-\frac12Z^2T^2,
\]

all pure Gaussian moments vanish, and \(Q=Z\) has nonzero mixed moments.
Adding unused independent Gaussian coordinates gives the same conclusion in
every higher dimension.  The exact finite regressions and the all-order
coefficient identity used in Long's proof are transcribed and checked in
finite ranges by
[`verify_long_gaussian_moments.py`](../scripts/verify_long_gaussian_moments.py).

Accordingly, the repository will not start or extend broad witness searches
in dimensions \(n\geq3\).  The sole active Gaussian-dimension problem is
\(\mathrm{GMC}(2)\).

## Active frontier

The cubic Gaussian null-cone theorem proves GMC(2) for every polynomial of
total degree at most three.  Therefore any counterexample must have:

- total degree at least four;
- mixed-sign rotational support; and
- nonlinear dependence on both circular coordinates.

At degree four, the exact three-level family
\(P=WA(U)+C(U)+ZB(U)\) with rotational support \(\{-1,0,1\}\) is also
closed: centering gives a Bessel--factorial functional equation, and moments
two through six generate the unit ideal on all four coefficient charts.  See
the [support-graph analysis](TWO_REAL_GMC_SUPPORT_GRAPH_EXPLORATION.md).

The next positive target is a degree-uniform null-cone theorem that separates
the radial factorial coefficients from zero-sum rotational-weight
combinatorics.  Its sharpest current subproblem is all-degree
three-level Bessel--factorial rigidity: prove that the functional equation in
the support-graph analysis forces both invariant polynomials to vanish.  New
Gaussian search code must either address this two-real frontier directly or
provide a finite, theorem-directed certificate for it.

## Retained high-dimensional material

High-dimensional results remain in the active reference corpus when their
purpose is logical transport rather than witness discovery:

- the fixed-dimensional
  \(\mathrm{GMC}(2r)\Rightarrow\mathrm{SIC}(r)\Rightarrow\) Keller
  invertibility implication;
- the BCW consequence chains to failures of GMC(42) and GMC(158);
- the four-real weighted-seed bridge and its transport of seed geometry into
  moment coordinates; and
- the associated Image/Vanishing consequences.

These are examples of how structure moves between conjectures and
categories.  They are not active attempts to optimize a Gaussian witness
dimension, because the five-term witness already settles every \(n\geq3\).

## Archive

The former three-real weighted-family searches are preserved under
[`archive/high-dimensional-gmc`](../archive/high-dimensional-gmc/README.md).
Their exact finite exclusions remain reproducible historical data, but their
surviving ansatzes are not an open research queue.
