# Tiny Image-Mathieu counterexamples: resolved pair dimension and compression

## 1. The absolute problem and the Jacobian-derived problem differ

For

\[
 \mathcal M_r
 =\sum_{i=1}^r(\partial_{z_i}-\zeta_i)
   \mathbb C[\zeta_1,\ldots,\zeta_r,z_1,\ldots,z_r],
\]

the absolute pair-dimension problem is now sharp.  With two pairs, define
\[
 R=\xi _1z_1+\xi _2z_2,\quad Z=\xi _1z_2,\quad
 W=2\xi _2z_1,\quad T=\xi _1z_1-\xi _2z_2
\]
and
\[
 F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right),\qquad Q=Z.
\tag{1.1}
\]
Then
\[
 \mathcal E_2(F^m)=0,\qquad
 \mathcal E_2(QF^m)=\frac{(4m+2)!\,m!}{(2m+1)!!}\ne0
\quad(m\geq1).
\tag{1.2}
\]
The all-order Hopf-coordinate and beta-integral proof is in
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
Since the one-pair case is known, the minimum failing pair dimension is
exactly two.

The previous support-minimal incumbent uses three contraction pairs.  With
pairs \((\tau,t),(w,z),(v,y)\), put

\[
 f=\tau(t-y)(wz+vt),
 \qquad g=y.                                              \tag{1.3}
\]

Then

\[
 f^m\in\mathcal M_3,\qquad gf^m\notin\mathcal M_3
 \quad\text{for every }m\geq1.                            \tag{1.4}
\]

The all-order proof and provenance are in
[`THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
In fact
\([t]\mathcal E_3(gf^m)=(-1)^{m-1}(m+1)!m!\).  Expanded, \(f\) has four
terms and ordinary degree four, while \(g\) is one linear monomial.  The
two-pair witness has sixteen terms, ordinary degree eight, bidegree
\((4,4)\), and full coefficient-matrix rank five.  Thus it closes pair
dimension by using the genuinely nonseparable option rather than linearly
compressing the three-pair form.

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
\text{arbitrary }\mathcal M_r& r=2,\ (16,8;1,2)&
 \text{support/degree minimization}\\
\text{derived from the stored Keller collision}&r=20,\ (60,4;1,1)&r<20.
\end{array}                                                \tag{1.5}
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

### Track A: residual two-pair classification and minimization

The absolute existence problem is closed by (1.1).  The remaining
two-pair questions are sharper: decide the full bidegree-\((3,3)\) stratum,
classify the moment-zero semistable locus in bidegree \((4,4)\), and
minimize support and ordinary degree.

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
The exact reduction of sparse contractions to finite hypergeometric
multisums, the recurrence-certificate protocol, and the separate
applicability gates for Newton-polytope GKZ rank bounds are recorded in
[`HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md`](HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md).

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

For pair dimension, higher-pair and separable searches are now obsolete.
The successful two-pair witness realizes the nonseparable route in the
earlier search list:

1. classify one-pair bihomogenizations that can absorb two circular pairs
   without adding a third pair;
2. search bidegree-balanced nonseparable \(f\) directly with (2.3);
3. exploit simultaneous \(\mathrm{GL}_2\) changes on both contraction pairs;
4. test whether the binomial cancellation in the three-pair witness has a
   two-pair polarization or diagonal quotient.

The first direct frontier is now closed in a complete bidegree.  The
[bidegree-\((2,2)\) frontier theorem](TWO_PAIR_SIC_BIDEGREE22_FRONTIER.md)
proves that the natural four-parameter linear compression of (1.1) is
one-sided whenever its first two pure contractions vanish.  More generally,
exact elimination and radical-power certificates prove that the first six
pure contractions cut out precisely the pair-linear one-sided nullcone for
every bidegree-\((2,2)\) form.  Thus SIC(2) holds on the complete
bidegree-\((2,2)\) stratum, including the nine exact eight-term supports and
the full support.  The \(501\)-chart sparse census and its explicit
six- and seven-term factorizations remain useful regressions.  The new
counterexample at bidegree \((4,4)\) shows that this safe result is sharp as
a low-degree boundary.

The next balanced case is now active.  The
[bidegree-\((3,3)\) frontier](TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md)
eliminates the full seven-dimensional pair-linear nullcone exactly.  On the
highest irreducible \(\operatorname{Sym}^6\) summand, moments of orders
\(2,4,6,10\) have precisely the binary-sextic \(L^4Q\) nullcone radical,
while moments \(2,3\) and \(2\) close the pure
\(\operatorname{Sym}^4\) and \(\operatorname{Sym}^2\) summands.  The
normalized non-null quadratic branch is also closed exactly on
\(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\), where moments through
order six force \(c^6=0\).  On
\(\operatorname{Sym}^6\oplus\operatorname{Sym}^2\), the analogous
\(c^{25}\) certificate is currently verified only over
\(\mathbb F_{32003}\).  The remaining target is the mixed
\(\operatorname{Sym}^6\oplus\operatorname{Sym}^4\oplus
\operatorname{Sym}^2\) representation after the first moment removes the
scalar summand.  Its continuation is now a degree-three classification
problem inside the formerly all-degree
[moment--nullcone program](TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md):
for
\(\operatorname{End}(\operatorname{Sym}^d)
\cong\bigoplus_{r=0}^d\operatorname{Sym}^{2r}\), first force the lowest
nonzero binary-form summand into its nullcone and then force every higher
summand to share its unique high-multiplicity root. This replaces deeper
full \((3,3)\) elimination by a quadratic-anchor lemma and a uniform
common-root synchronization lemma.  The two-pair counterexample falsifies
the uniform statement at \(d=4\), but it does not settle \(d=3\).
For \(d=3\), the first thirteen full moments now have an exact nonzero
\(13\times13\) Jacobian minor, so they are algebraically independent and
attain the invariant-quotient dimension bound. The remaining question is
their zero fiber, not a shortage of independent moment coordinates.

The bidegree-\((4,4)\) witness supplies the required all-order identities.
Bounded zero prefixes remain insufficient for any new candidate or
classification.

For homogeneous two-variable GVC, this search has a smaller
representation-theoretic target.  If \(A\) is the degree-\(d\) binary
symbol of \(\Lambda\) and \(P\) is a degree-\(d\) binary form, then
\(A(\zeta)P(z)\) is a rank-one point of
\(\operatorname{End}(\operatorname{Sym}^d)\), and
\[
 \mathcal E_2\!\left((A(\zeta)P(z))^m\right)
 =\Lambda^m(P^m).
\]
Thus balanced GVC(2) is the Segre restriction of the SIC(2) moments.  The
[balanced cubic theorem](TWO_VARIABLE_CUBIC_GVC_THEOREM.md) now closes the
entire cubic Segre slice: the first four moments force a one-sided form and
give \(\Lambda^m(QP^m)=0\) for \(m>\deg Q\).  The first open separable
balanced degree is four.  The
[low-root continuation](TWO_VARIABLE_LOW_ROOT_GVC_THEOREMS.md) then closes
all degrees for symbols with at most two distinct roots and closes the
quartic \((2,1,1)\) orbit by a first-five-moment radical certificate.
The [split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) goes further:
translated complete polarization and the Laurent constant-term theorem
prove the GVC conclusion for every homogeneous binary operator and every
polynomial, with no degree restriction.  Thus the remaining GVC(2) route
needs a genuinely nonhomogeneous operator; the
[separable escape obstruction](SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
also excludes degree at most the lowest positive order and the
arbitrary-degree factor-unit class.  Its low-order orbit analysis now
closes every binary operator through polynomial degree four, so a GVC
counterexample must begin in degree five.  It further closes every
quadratic-leading pair through degree five, so the first degree-five
counterexample would have lowest positive operator order three or four.
Finite-moment nullcone classification remains a separate SIC-strengthening.
The continuing representation program and the distinction between unequal
bidegree and genuinely nonhomogeneous mixtures are recorded in
[`TWO_VARIABLE_GVC_REPRESENTATION_PROGRAM.md`](TWO_VARIABLE_GVC_REPRESENTATION_PROGRAM.md).

### Track B: minimal translation-uniform lift of a GMC witness

Treat a centered Gaussian counterexample as seed data and add the fewest
auxiliary contraction pairs needed to make every translated pure moment
vanish.  In generating-function language, the determinant or Jacobian
factor must cancel identically as a polynomial in all translation
parameters, not only at the origin.  Long's sparse three-real witness and
the repository's weighted Gaussian families are the first seeds to test.

The one-pair bihomogenization producing (1.3) remains the benchmark for
translation-uniform lifts.  The sharp two-pair witness (1.1) instead has
full coefficient-matrix rank, confirming that genuine nonseparability can
beat this lift architecture.

### Track C: provenance-preserving Keller compression

For a Keller map \(F=z-H(z)\), the canonical contraction candidate is

\[
 f=\zeta\mathbin\cdot H(z),
\]

and a nonpolynomial inverse coordinate supplies \(g=z_i\).  For the stored
identity slice this gives the current 20-pair, 60-term incumbent.

At pair dimension two, the canonical route is now closed under the
normalized Keller constraint.  The
[dual-linear theorem](DUAL_LINEAR_SIC2.md) proves more generally that the
first two pure equations close every
\(p=\zeta\mathbin{\cdot}H\), without a Keller hypothesis.  On the
normalized Keller subclass, its sharper first-moment corollary says that
\(\mathcal E_2(\zeta\mathbin\cdot H)=\operatorname{tr}JH=0\), together with
\(\det J(z-H)=1\), forces
\[
 H=(b,-a)f(ax+by).
\]
It follows that
\(\mathcal E_2(g(\zeta\mathbin\cdot H)^m)=0\) for
\(m>\deg_z g\).  Thus a two-pair SIC counterexample cannot be the canonical
contraction of a normalized plane Keller correction.  This excludes only
Keller provenance; it does not close Track A.

Further linear output or source quotients of that artifact have already
been exhausted.  The first nonlinear identity-slice class is now exhausted
as well: the
[Keller-provenance compression audit](KELLER_PROVENANCE_COMPRESSION.md)
proves
\[
 \{P\in\mathbb Q[X]_{\leq6}:P\circ V=P\}
 =\mathbb Q[X_{20}]_{\leq6}.
\]
Thus no invariant polynomial of degree at most six exposes a second identity
output independent of the known homogenizing coordinate.  The same audit
shows that literal inverse-recurrence dependency pruning still uses all
twenty active coordinates and that no already stored degree-lowering circuit
has a smaller identity slice.  The invariant statement is an exact nonlinear
class obstruction; the recurrence and circuit conclusions have their
narrower scopes stated in that note.

At degree six, both possible correction channels from the lower-degree
near-invariant are also closed. Two exact torus gradings reduce them to
Lie-derivative blocks of sizes 103 and 1604; both blocks are injective and
both required defects lie outside their images. The same gradings classify
all \(220\) sextic sectors. Unique-row peeling certifies \(25\) of the \(28\)
dense sectors outright and reduces the other three to cores of sizes
\(1596,1956,2170\). The complete Lie kernel is
\(\langle s^6,s^4Q,s^2Q^2,Q^3\rangle\); exact pullback defects leave only
\(s^6\) fixed. At degree seven, both lower-degree correction sectors are
exactly \(s\) times the already excluded sextic sectors, so only the pure
homogeneous septic Lie kernel remains open. Reducing that problem modulo
\(s\) gives \(657800\) columns in \(204\) sectors; exact unique-row peeling
removes \(451891\), leaving a \(205909\)-column residual quotient kernel and
its \(s\)-adic lifting problem. The reduced derivation is a constant vertical
direction over the fourteen-variable base. Its first lifting obstruction
lives in the quotient by the six vertical coefficients; it already excludes
the extreme class \(X_9^7\). That coefficient ideal has height two and five
minimal components, including the plane \(X_0=X_1=0\); its degree-eight
quotient has dimension \(158412\). Consequently a regular-sequence/Koszul
shortcut is unavailable. Direct evaluation on the five components excludes
\(71588\) of the \(77520\) support-one base septics. The support-free stacked
component map has rank at least \(61060\), leaving a radical-level subspace
of dimension at most \(16460\); the next attack is embedded-torsion control.

The same audit program now backtraces the near-invariant through the frozen
BCW circuit. In the 24-dimensional rank-compressed homogenization,
\[
Q=c_4s-v_3v_5,
\]
and on the stable source section \(v_3=-xz,\ v_5=-xy\), so \(Q=x^2yz=M\).
Thus \(Q(V)-Q=-Ms^2\) is a shared-factor gate-residual identity. This explains
the low-degree rank-one resemblance to the Hopf construction without
identifying the two cancellation mechanisms.

There is also an exact lower bound for non-invariant semiconjugacies. On the
twenty-variable slice, the first thirteen iterates of both the multiplier
\(X_0\) and the quadratic collision observable \(X_{18}-X_6X_8\) are
algebraically independent. Any rational semiconjugate quotient carrying
either observable therefore has dimension at least thirteen. The modular
rank remains thirteen through twenty-five iterates at the tested point, but
that plateau is experimental and does not prove a thirteen-dimensional
quotient.

Resume this track only with one of:

1. a nonlinear collision-preserving quotient outside the degree-at-most-six
   invariant-slice class;
2. a different degree-lowering/homogenizing circuit before contraction;
3. a new smaller noninvertible Keller presentation;
4. a proof that a contraction witness derived from the inverse recurrence
   can be represented with fewer pairs than the Keller map itself.

Every candidate must retain an explicit map from the Keller collision to
\((f,g)\); otherwise it belongs to Track A or B rather than this track.

### Track D: rank-efficient ordinary-Laplacian polarization

The third-order Dvorsky symbol
\[
 \xi_t(\xi_a\xi_d-\xi_b\xi_c)
\]
suggests a structured quadraticization rather than a generic doubling.
The target is an all-order lift to a nondegenerate quadratic
constant-coefficient operator that preserves both the pure vanishing and
the fixed-multiplier defect.  The auxiliary-variable count should be
minimized by Schur complements and rank stratification, retaining the
determinant factor as one rank-two block.

A constant Schur complement cannot change operator order, so the lift must
couple the operator and polynomial sides or use nonlinear polarization.
Finite moment matching is insufficient.  The first complete ansatz has now
been excluded: after adding one variable \(s\), the nondegenerate operator
\[
 \partial_a\partial_d-\partial_b\partial_c+\partial_t\partial_s
\]
admits no homogeneous cubic lift restricting to the Dvorsky polynomial on
\(s=0\) whose first two pure moments vanish.  The
[exact obstruction](DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md) is the
unrestricted transverse-jet identity \(12t^2-8\rho t\) in the second
moment.  It applies to every polynomial or formal harmonic hyperplane
lift, so degree mixing cannot repair this completion.  Therefore the next
candidate must use a second auxiliary block, a different quadratic
completion, or a nonlinear specialization.
The precise identities, rank-minimization objective, and promotion gates
are in
[`TWO_VARIABLE_GVC_REPRESENTATION_PROGRAM.md`](TWO_VARIABLE_GVC_REPRESENTATION_PROGRAM.md).

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
