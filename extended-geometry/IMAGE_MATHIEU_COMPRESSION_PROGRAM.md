# Tiny Image-Mathieu counterexamples: compression program

## 1. The absolute problem and the Jacobian-derived problem differ

For

\[
 \mathcal M_r
 =\sum_{i=1}^r(\partial_{z_i}-\zeta_i)
   \mathbb C[\zeta_1,\ldots,\zeta_r,z_1,\ldots,z_r],
\]

the repository now has an eight-term counterexample in three contraction
pairs.  With pairs \((\tau,t),(w,z),(v,y)\), put

\[
 f=\tau^3(t+z)\bigl(wt^3-vy(t+y)^2\bigr),
 \qquad g=z.                                              \tag{1.1}
\]

Then

\[
 f^m\in\mathcal M_3,\qquad gf^m\notin\mathcal M_3
 \quad\text{for every }m\geq1.                            \tag{1.2}
\]

The all-order proof and provenance are in
[`THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
In fact
\([t]\mathcal E_3(gf^m)=(3m+1)!m!\).  Expanded, \(f\) has eight terms and
ordinary degree eight, while \(g\) is one linear monomial.  The
one-pair case is known to satisfy the Image Conjecture, so the remaining
absolute dimension question is exactly whether \(\operatorname{SIC}(2)\)
holds.

The earlier Dvorsky--Long witness remains the smaller-degree separable
example: it uses five pairs, eight degree-six terms, and a linear
multiplier.  Its all-order proof and external provenance are in
[`DVORSKY_GVC5_COUNTEREXAMPLE.md`](DVORSKY_GVC5_COUNTEREXAMPLE.md).

The specifically **Jacobian-derived** frontier is different.  The identity
slice of the stored cubic Keller collision gives a 20-pair witness with a
60-term polynomial of maximum degree four and a linear multiplier.  Its
proof is in
[`IMAGE_VANISHING_COUNTEREXAMPLES.md`](IMAGE_VANISHING_COUNTEREXAMPLES.md).
Linear collision-preserving compression of that stored witness is exhausted;
nonlinear compression remains open.  A small witness such as (1.1) does not
automatically retain Jacobian provenance.

These should be kept as two separate optimization problems:

\[
\begin{array}{c|c|c}
\text{problem}&\text{certified incumbent}&\text{next target}\\ \hline
\text{arbitrary }\mathcal M_r& r=3,\ (8,8;1,1)&r=2\\
\text{derived from the stored Keller collision}&r=20,\ (60,4;1,1)&r<20.
\end{array}                                                \tag{1.3}
\]

Here \((s,d;s_g,d_g)\) records expanded term count and degree for \(f\), then
for \(g\).  It is only a complexity ledger, not a mathematical invariant.

## 2. Exact membership identities

Define the contraction

\[
 \mathcal E_r(\zeta^\alpha q(z))=\partial_z^\alpha q(z).
\]

Zhao's image-kernel identity says

\[
 \mathcal M_r=\ker\mathcal E_r.                            \tag{2.1}
\]

If

\[
 h=\sum_{\alpha,\beta}c_{\alpha,\beta}\zeta^\alpha z^\beta,
\]

then

\[
 \mathcal E_r(h)
 =\sum_{\rho\in\mathbb N^r}
 \left(
  \sum_{\alpha}
  c_{\alpha,\alpha+\rho}\frac{(\alpha+\rho)!}{\rho!}
 \right)z^\rho.                                             \tag{2.2}
\]

Consequently

\[
 \boxed{
 h\in\mathcal M_r
 \Longleftrightarrow
 \sum_\alpha(\alpha+\rho)!\,c_{\alpha,\alpha+\rho}=0
 \quad\text{for every }\rho.
 }                                                          \tag{2.3}
\]

This is the preferred exact finite membership test at any fixed exponent.
It uses integer arithmetic after denominators in \(h\) are cleared.

## 3. The translated Gaussian formulation

Let \(W_i,Z_i\) be independent circular Gaussian pairs normalized by

\[
 \mathbb E(W^\alpha Z^\beta)
 =\delta_{\alpha,\beta}\alpha!.
\]

For \(\tau_a h(\zeta,z)=h(\zeta,z+a)\), the scalar Gaussian functional
\(\mathcal F_r(h)=\mathbb E(h(W,Z))\) satisfies

\[
 \mathcal F_r(\tau_a h)=\mathcal E_r(h)(a).                 \tag{3.1}
\]

Therefore

\[
 \boxed{
 h\in\mathcal M_r
 \Longleftrightarrow
 \mathbb E\!\left[h(W,Z+a)\right]=0
 \quad\text{for every }a\in\mathbb C^r.
 }                                                          \tag{3.2}
\]

This is the precise Gaussian bridge.  Centered Gaussian moment vanishing
alone proves membership in \(\ker\mathcal F_r\), not in \(\mathcal M_r\).
Any direct use of a GMC witness must therefore solve the additional
translation-uniformity problem.  The Dvorsky--Long witness in (1.1) satisfies
this stronger condition; a naive substitution of Long's three-real-variable
Gaussian polynomial does not.

For a candidate pair \(f,g\), the Mathieu failure becomes

\[
\begin{aligned}
 &\mathbb E[f(W,Z+a)^m]=0
   &&\text{for every }a\text{ and }m\geq1,\\
 &\mathbb E[g(W,Z+a)f(W,Z+a)^m]\ne0
   &&\text{as a polynomial in }a
     \text{ for infinitely many }m.                         \tag{3.3}
\end{aligned}
\]

This formulation lets the existing exact Wick, coefficient, and
moment-ideal machinery attack the Image problem without confusing a bounded
moment search with an all-order proof.

## 4. Search tracks

### Track A: the remaining two-pair problem

Start from factorizations

\[
 f=\lambda(\zeta)P(z)
\]

because then

\[
 \mathcal E_r(f^m)
 =\lambda(\partial_z)^mP(z)^m.                              \tag{4.1}
\]

This converts the pure-power condition into generalized-vanishing
identities and exposes recurrences suitable for creative telescoping,
Gröbner elimination, or finite-difference proofs.  The existing bounded
four-variable search in
[`DVORSKY_GVC5_COUNTEREXAMPLE.md`](DVORSKY_GVC5_COUNTEREXAMPLE.md)
only covers one symmetry-preserving lattice slice.

The following four-pair subproblem was closed before the three-pair
construction was found and remains useful structural evidence:
[`FOUR_PAIR_DVORSKY_SLICE_OBSTRUCTION.md`](FOUR_PAIR_DVORSKY_SLICE_OBSTRUCTION.md).
For

\[
 P=(t+a+b+d)(ad+bt)
\]

and the complete six-parameter family
\(\lambda=\zeta_tR(\zeta_a,\zeta_b,\zeta_d)\) with \(R\) quadratic, the
first eight pure identities leave only four rank-one square directions.
Every direction is one-sided and satisfies eventual mixed vanishing for
every fixed multiplier.  Thus this slice contains no counterexample over
any characteristic-zero field; it is no longer merely a bounded negative
search.

For improving the pair dimension, four-pair separable searches are now
obsolete.  The target is \(r=2\), where the three-pair construction suggests:

1. classify one-pair bihomogenizations that can absorb two circular pairs
   without adding a third pair;
2. search bidegree-balanced nonseparable \(f\) directly with (2.3);
3. exploit simultaneous \(\mathrm{GL}_2\) changes on both contraction pairs;
4. test whether the binomial cancellation in the three-pair witness has a
   two-pair polarization or diagonal quotient.

The first goal is not merely many vanishing initial powers.  It is an
all-order recurrence whose zero solution is forced by finitely many initial
conditions, together with a mixed recurrence having infinitely many
certified nonzero terms.

### Track B: minimal translation-uniform lift of a GMC witness

Treat a centered Gaussian counterexample as seed data and add the fewest
auxiliary contraction pairs needed to make every translated pure moment
vanish.  In generating-function language, the determinant or Jacobian
factor must cancel identically as a polynomial in all translation
parameters, not only at the origin.  Long's sparse three-real witness and
the repository's weighted Gaussian families are the first seeds to test.

The benchmark is now the one-pair bihomogenization producing (1.1).  A
two-pair construction would settle the minimum exactly, while a theorem
excluding all two-pair bihomogenizations of fixed-\(W\)-degree circular
witnesses would isolate the genuinely nonseparable remainder.

### Track C: provenance-preserving Keller compression

For a Keller map \(F=z-H(z)\), the canonical contraction candidate is

\[
 f=\zeta\mathbin\cdot H(z),
\]

and a nonpolynomial inverse coordinate supplies \(g=z_i\).  For the stored
identity slice this gives the current 20-pair, 60-term incumbent.

Further linear output or source quotients of that artifact have already
been exhausted.  Resume this track only with one of:

1. a nonlinear collision-preserving quotient;
2. a different degree-lowering/homogenizing circuit before contraction;
3. a new smaller noninvertible Keller presentation;
4. a proof that a contraction witness derived from the inverse recurrence
   can be represented with fewer pairs than the Keller map itself.

Every candidate must retain an explicit map from the Keller collision to
\((f,g)\); otherwise it belongs to Track A or B rather than this track.

## 5. Promotion rules

A promoted counterexample must include:

1. explicit \(f,g\) over a named characteristic-zero field;
2. a written all-order proof of \(f^m\in\mathcal M_r\);
3. a written proof that \(gf^m\notin\mathcal M_r\) infinitely often or
   eventually;
4. an exact sparse replay of a useful finite prefix using (2.3);
5. term counts, degrees, normalization, and provenance.

A bounded search or a long zero prefix remains an experiment.  It must not
change the counterexample scoreboard or `MATH_STATUS.json`.
