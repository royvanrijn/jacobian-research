# Support graphs for the two-real Gaussian moment problem

## 1. Updated starting point

For

\[
Z=(X+iY)/\sqrt2,\qquad W=(X-iY)/\sqrt2,\qquad U=ZW,
\]

one has

\[
\mathbb E(Z^aW^b)=\delta_{a,b}a!.
\]

Give \(Z,W\) rotational weights \(1,-1\).  The repository now proves more
than the preliminary support census: every polynomial of total degree at
most three satisfies \(\operatorname{GMC}(2)\).  Thus a possible
counterexample has degree at least four.

The proposed support-classification direction is still natural, but
“weight-support graph” needs a fixed definition.  Two plausible definitions
give different smallest cyclic supports.

## 2. Canonical circuit-incidence graph

Embed

\[
\mathbb C[Z,W]\hookrightarrow\mathbb C[U][T,T^{-1}],
\qquad Z\mapsto T,\quad W\mapsto UT^{-1}.
\]

Write

\[
P=\sum_{k\in S}T^kB_k(U).
\]

The expectation is

\[
\mathbb E(P^m)=
\mathcal L\!\left(\operatorname{CT}_T P^m\right),
\qquad \mathcal L(U^j)=j!.                         \tag{2.1}
\]

A canonical combinatorial object is the incidence graph between:

- the nonzero weights \(k\in S\setminus\{0\}\); and
- the primitive support-minimal nonnegative relations
  \(\alpha\) satisfying \(\sum_k\alpha_k k=0\).

For a rank-one weight set, every such circuit has one positive and one
negative weight.  If \(p\) positive and \(q\) negative levels occur, the
incidence graph is the subdivision of \(K_{p,q}\).  Consequently it is a
forest exactly when

\[
\min(p,q)\leq1.                                    \tag{2.2}
\]

The zero weight contributes a singleton invariant and does not create a
cycle.  Under this definition, the support \(\{-1,0,1\}\) is a forest.  The
first cycle is \(K_{2,2}\), so it requires at least two positive and two
negative nonzero weight levels.

This gives a precise restricted conjecture:

> **Circuit-forest conjecture.** If the nonzero rotational support of \(P\)
> has at most one positive level or at most one negative level, then \(P\)
> satisfies \(\operatorname{GMC}(2)\).

The known two-weight theorem is only the \(p=q=1\) base case.  The conjecture
also contains genuinely new stars \(K_{1,q}\) and \(K_{p,1}\).

## 3. Expanded coefficient graph

One may instead expand every \(B_k(U)\) into radial monomials and retain
parallel invariant products as separate edges.  This can create cycles even
for weights \(\{-1,0,1\}\).  But contracting the \(U=ZW\) direction merges
those parallel edges and returns the circuit-incidence graph above.

Therefore the two claims

1. “contract invariant \(ZW\)-coefficients,” and
2. “the first cyclic case is \(\{-1,0,1\}\),”

cannot both hold without extra graph structure.  Any theorem or search
should specify whether radial coefficient polynomials are vertices, edges,
labels, or contracted data.

## 4. The exact three-level functional equation

The smallest interaction between two different invariant polynomials is
nevertheless the three-level ansatz

\[
P=W A(U)+C(U)+ZB(U).
\]

Put

\[
D(U)=U A(U)B(U).
\]

Constant-term extraction gives

\[
\begin{aligned}
\mathbb E(e^{tP})
&=\mathcal L\!\left[
 e^{tC(U)}
 \sum_{r\geq0}\frac{t^{2r}D(U)^r}{(r!)^2}
\right]\\
&=\mathcal L\!\left[
 e^{tC(U)}I_0\!\left(2t\sqrt{D(U)}\right)
\right].                                           \tag{4.1}
\end{aligned}
\]

Thus all pure moments vanish exactly when (4.1) equals \(1\).  Equivalently,
for every \(m\geq1\),

\[
\sum_{r=0}^{\lfloor m/2\rfloor}
\frac{m!}{(r!)^2(m-2r)!}
\mathcal L\!\left(C^{m-2r}D^r\right)=0.            \tag{4.2}
\]

This is the desired finite-alphabet functional equation: it depends only on
the two invariant polynomials \(C,D\).  If \(C=0\), the even equations and
the one-variable Factorial Theorem give \(D=0\).  If \(D=0\), the same
theorem gives \(C=0\).  Hence any new phenomenon must be a coupled
cancellation between \(C\) and \(D\).

The ordinary moment series is even simpler:

\[
\boxed{
\sum_{m\geq0}\mathbb E(P^m)t^m
=
\mathcal L\!\left(
\frac1{\sqrt{(1-tC)^2-4t^2D}}
\right).}                                         \tag{4.3}
\]

Indeed, before applying \(\mathcal L\), the right side is the constant term
of

\[
\frac1{1-t(TB+C+T^{-1}UA)}.
\]

This algebraic resolvent is preferable to the Bessel series for
classification.  It also shows that the three-level problem depends only on
\(C\) and \(D\), not on the individual factorization \(D=UAB\).

The mixed GMC moments have an equally useful marked form.  For
\(Q=T^\ell R(U)\), only indices with the appropriate weight imbalance
survive in

\[
\mathcal L\!\left(
\operatorname{CT}_T
T^\ell R(U)
\exp(t(TB+C+T^{-1}UA))
\right).                                           \tag{4.4}
\]

Once a solution of (4.1) is found, (4.4) is a direct finite search for a
multiplier witnessing failure.

## 5. Exact low-degree theorem

In total degree at most \(n\), put

\[
r=\lfloor n/2\rfloor,\qquad s=\lfloor(n-1)/2\rfloor.
\]

After the first moment centers \(C\), the invariant coordinates are

\[
\begin{aligned}
C&=\sum_{j=1}^{r}c_j(U^j-j!),\\
D&=\sum_{j=1}^{2s+1}d_jU^j.                       \tag{5.1}
\end{aligned}
\]

Conversely, every nonzero \(D\) in (5.1) factors over \(\mathbb C\) as
\(UAB\) with \(\deg A,\deg B\leq s\): split the at most \(2s\) roots of
\(D/U\) into two groups.  Thus (5.1) exactly parameterizes this support
after quotienting the relative \(A/B\) scaling.

Coefficient charts cover \(C\ne0,D\ne0\).  Weighted scalar multiplication
normalizes one selected \(c_j\), a selected \(d_j\) is localized, and the
second moment eliminates the highest \(D\)-coefficient.  Exact
characteristic-zero calculations give:

| degree | invariant charts | moment cutoff | result |
|---:|---:|---:|---|
| 4 | 6 | 6 | reduced basis \([1]\) on every chart |
| 5 | 10 | 8 | reduced basis \([1]\) on every chart |
| 6 | 15 | 9 | reduced basis \([1]\) on every chart |

Therefore:

> **Low-degree three-level theorem.** No polynomial of total degree at most six
> with all three weight components \(-1,0,1\) nonzero has all pure Gaussian
> moments zero.

The exact replay is
[`verify_two_real_gmc_three_weight_low_degree.py`](../scripts/verify_two_real_gmc_three_weight_low_degree.py).

This closes the proposed smallest test through degree six.  It does not prove the
circuit-forest conjecture for arbitrary radial degree.

## 6. Next finite searches

The two graph interpretations lead to different next computations.

1. **Test the circuit-forest conjecture.** Fix a star support, beginning with
   \(\{-2,-1,0,1\}\) or its reflection, and increase radial degree.  This
   tests a genuinely new forest beyond a single positive/negative pair.
2. **Test the first circuit cycle.** Use four nonzero levels, beginning with
   \(\{-2,-1,1,2\}\), optionally with weight zero.  The circuit graph is
   \(K_{2,2}\), with first Betti number one.
3. **Test two independent graph cycles.** The smallest canonical graph has
   \((p,q)=(2,3)\) or \((3,2)\), giving
   \[
   b_1(K_{p,q})=pq-p-q+1=2.
   \]
   A minimal weight choice is
   \(\{-3,-2,-1,1,2\}\), up to reflection.
4. **Classify the three-level functional equation.** Seek a theorem that
   (4.1) forces \(C=D=0\).  This would prove the whole
   \(\{-1,0,1\}\) family in all degrees and isolate a reusable
   Bessel-factorial rigidity statement.

The most promising theorem target is item 4.  It is stronger than the
low-degree eliminations, has no graph-definition ambiguity, and turns a
counterexample into two explicit univariate polynomials satisfying (4.2).
Degree seven is the next finite test, but after the invariant reductions its
localized systems still have nine variables.  A sparse coefficient-exponent
generator computes its moments through order ten in seconds, so moment
generation is no longer the obstruction.  On the representative chart
\(c_1d_1\ne0\), however, the nine-equation cutoff-ten system did not finish
modulo \(p=1\,000\,003\) within five minutes.  This timeout is not evidence
of a solution or nonemptiness.  It identifies the next required method:
multihomogeneous elimination, a better variable order, or a theorem based on
the algebraic resolvent (4.3).  Repeating the same dense Gröbner calculation
on the other 20 charts has no evidentiary value.

## 7. Literature connection

The support program sits at the intersection of three existing theorems,
but is not currently implied by any of them.

1. [Derksen--van den Essen--Zhao](https://arxiv.org/abs/1506.05192)
   introduced GMC, proved \(\operatorname{GMC}(1)\), and proved
   \(\operatorname{GMC}(2)\) for homogeneous polynomials.  Their homogeneous
   proof separates radius from angle and applies the torus Mathieu theorem.
   Equation (4.1) identifies why the inhomogeneous case is different: the
   factorial functional couples radial degrees after angular
   constant-term extraction.
2. [Duistermaat--van der Kallen](https://doi.org/10.1016/S0019-3577(98)80020-7)
   classified Laurent polynomials whose positive powers have zero constant
   term; geometrically, vanishing forces the Newton support away from the
   origin.  This is the direct ancestor of the weight-null-cone and
   support-graph viewpoint.  It cannot be applied after (2.1), because
   \(\mathcal L\) can cancel different radial coefficients even when the
   constant term itself is nonzero.
3. [van den Essen--Wright--Zhao](https://arxiv.org/abs/1008.3962)
   introduced the Factorial Conjecture and proved its one-variable case
   (Theorem 4.9).  That result closes the axes \(C=0\) and \(D=0\) in
   (4.1), but does not address their Bessel-weighted coupling.
4. [Long](https://arxiv.org/abs/2607.18186) gives explicit failures in
   every dimension at least three and leaves dimension two as the only
   unresolved dimension.  His examples show that Lagrange/determinant
   cancellation can defeat pure moments, but their auxiliary Gaussian
   source has no analogue in (4.1).  The Bessel-factorial equation is the
   exact two-real-variable replacement problem.

Accordingly, a proof that (4.1) forces \(C=D=0\) would be a new rigidity
theorem interpolating between the Duistermaat--van der Kallen
constant-term theorem and the one-variable Factorial Theorem.  A
nontrivial solution would be comparably significant: it would pass the
pure-moment gate in the only unresolved Gaussian dimension, after which
(4.4) tests whether it is an actual GMC counterexample.

## 8. All-degree proof program

The low-degree theorem changes the role of computation.  Degree seven is a
useful regression target for a new method, but repeating dense localized
Gröbner calculations is not the primary program.  The following three
reformulations isolate more structural attacks.

### 8.1 Laplace--Bessel saddles and the formal-series caveat

On polynomials,

\[
 \mathcal L(f)=\int_0^\infty e^{-u}f(u)\,du.
\]

If

\[
 c=\deg C,\qquad d=\deg D,\qquad q=\max(c,d/2),
\]

the two large-argument Bessel phases are formally

\[
 \Phi_\pm(u,t)=-u+tC(u)\pm2t\sqrt{D(u)}.           \tag{8.1}
\]

For \(q>1\), a generic saddle has

\[
 u\asymp t^{-1/(q-1)},\qquad
 \Phi_\pm(u,t)\asymp t^{-1/(q-1)}.                \tag{8.2}
\]

Thus the leading saddle slopes are separated unless \(d=2c\).  On that
balanced line the first collision occurs precisely when one of the leading
coefficients

\[
 c_c\pm2\sqrt{d_{2c}}
\]

vanishes, equivalently \(c_c^2=4d_{2c}\).  Subsequent collisions are read
from the Newton polygon of

\[
 C(U)^2-4D(U).                                    \tag{8.3}
\]

There is an important analytic qualification.  Equations (4.1)--(4.3) are
identities of formal moment series; in radial degree at least two those
series generally have zero radius of convergence.  All moments vanishing
therefore says that a sectorial Laplace integral has asymptotic expansion
\(1\), not automatically that the analytic integral is identically \(1\).
Nontrivial saddle terms may be exponentially flat at \(t=0\).  A complete
asymptotic proof must consequently establish that every nonzero pair
\((C,D)\) has a nonzero Stokes contribution on some admissible \(t\)-ray.
The Newton classification above gives the finite list of leading collisions
that such a proof must resolve.

### 8.2 Bessel ODE and creative telescoping

Put

\[
 F(t,U)=e^{tC(U)}I_0(2t\sqrt{D(U)}).
\]

The Bessel equation gives the exact differential identity

\[
 \left[
 t^2\partial_t^2+(t-2t^2C)\partial_t
 +t^2(C^2-4D)-tC
 \right]F=0.                                      \tag{8.4}
\]

The factorial functional has the adjoint relation

\[
 \mathcal L(Uf)=\mathcal L((U\partial_U+1)f).      \tag{8.5}
\]

More generally, multiplication by \(U^k\) can be transferred one factor at
a time, producing the rising Euler operator

\[
 (U\partial_U+1)(U\partial_U+2)\cdots
 (U\partial_U+k).                                 \tag{8.6}
\]

For symbolic elimination it is preferable to use the algebraic resolvent

\[
 G(t,U)=Q(t,U)^{-1/2},\qquad
 Q=(1-tC)^2-4t^2D,
\]

which satisfies

\[
 2Q\,\partial_tG+Q_tG=0,\qquad
 2Q\,\partial_UG+Q_UG=0.                          \tag{8.7}
\]

The concrete creative-telescoping target is an operator
\(T(t,\partial_t)\) and an algebraic certificate \(R(t,U)G\) such that

\[
 T(G)=(\partial_U-1)(R\,G)+B(t,U).                \tag{8.8}
\]

After integration against \(e^{-U}dU\), the first term is a boundary term.
Substituting the proposed solution \(\mathcal L(G)=1\) into the resulting
ODE gives algebraic conditions on the coefficients of \(C,D\).  This route
also exposes the irregular singularities where exponentially flat solutions
can occur.  Formal substitution alone does not rule those solutions out;
boundary conditions and Stokes data must still be tracked.  A useful next
computation is to derive the minimal telescoper for generic symbolic degrees
\((c,d)=(2,3)\), then \((3,5)\), and identify the coefficient pattern before
attempting arbitrary degree.

There is already a finite first-order closure before scalar telescoping.
Write

\[
 Q(t,U)=\sum_{j=0}^s q_j(t)U^j,\qquad
 M_k(t)=\mathcal L(U^kG(t,U)).
\]

Integration by parts in the \(U\)-equation of (8.7) gives, for every
\(k\geq0\),

\[
\boxed{
 2\sum_{j=0}^s q_jM_{k+j}
 -\sum_{j=0}^s(2k+j)q_jM_{k+j-1}
 =2\delta_{k,0}Q(t,0)G(t,0),
}                                                  \tag{8.9}
\]

where terms with negative subscripts are omitted.  The coefficient of
\(M_{k+s}\) is \(2q_s\).  Hence, after localizing at the leading coefficient,
every radial moment reduces recursively to

\[
 M_0,\ldots,M_{s-1}.                              \tag{8.10}
\]

Multiplying the \(t\)-equation of (8.7) by \(U^k\) and applying
\(\mathcal L\) gives

\[
 2\sum_{j=0}^s q_jM'_{k+j}
 +\sum_{j=0}^s q'_jM_{k+j}=0.                    \tag{8.11}
\]

Differentiate the reductions (8.9), substitute them into (8.11) for
\(0\leq k<s\), and solve the resulting linear system.  This produces an
explicit meromorphic Pfaffian system

\[
 \partial_t
 \begin{pmatrix}M_0\\ \vdots\\ M_{s-1}\end{pmatrix}
 =
 A(t)
 \begin{pmatrix}M_0\\ \vdots\\ M_{s-1}\end{pmatrix}
 +b(t).                                           \tag{8.12}
\]

The initial data are \(M_k(0)=k!\).  Three-level rigidity is therefore
equivalent to showing that this distinguished formal horizontal section
cannot remain in the affine hyperplane \(M_0=1\) unless \(C=D=0\).
Eliminating the other \(s-1\) coordinates yields the scalar telescoper in
(8.8).  The poles introduced by \(q_s(t)\), which normally vanishes at
\(t=0\), make the expected irregular singularity explicit rather than
hiding it in moment growth.

The exact regression
[`verify_two_real_gmc_resolvent_system.py`](../scripts/verify_two_real_gmc_resolvent_system.py)
constructs the four-dimensional system for a centered degree-\((2,3)\)
pair and checks it against the factorial moment series.

### 8.3 Laguerre operator form

Let \(L_n(U)\) denote the ordinary Laguerre polynomials.  They satisfy

\[
 \int_0^\infty e^{-U}L_n(U)L_m(U)\,dU=\delta_{n,m},
\]

so \(\mathcal L\) selects the \(L_0\)-coefficient.  Multiplication by \(U\)
is the tridiagonal Jacobi operator

\[
 U L_n=(2n+1)L_n-(n+1)L_{n+1}-nL_{n-1}.           \tag{8.13}
\]

If \(J\) denotes this operator and

\[
 H_m(C,D)=
 \sum_{r=0}^{\lfloor m/2\rfloor}
 \frac{m!}{(r!)^2(m-2r)!}C^{m-2r}D^r,
\]

then the moment equations become

\[
 \langle e_0,H_m(C(J),D(J))e_0\rangle=0
 \quad(m\geq1).                                  \tag{8.14}
\]

Here \(C(J)\) and \(D(J)\) are commuting finite-band operators.  This turns
all-degree rigidity into a boundary spectral-moment problem for one fixed
Jacobi matrix.  It also supplies a sparse exact computational basis: the
outermost Laguerre bands and their unique extremal paths can be studied
before any coefficient-chart localization.

Of these approaches, (8.7)--(8.8) is the most concrete theorem engine.
The saddle analysis identifies its irregular singularities and exceptional
Newton polygons; the Laguerre form is the best candidate for exposing a
triangular or positivity-free band argument.

## 9. The first genuine circuit cycle

For discovery beyond the three-level forest, omit weight zero initially and
write the \(K_{2,2}\) support as

\[
 P=T^2B_2+TB_1+T^{-1}UA_1+T^{-2}U^2A_2.
\]

Its four primitive circuit invariants are

\[
\begin{aligned}
 X_{11}&=UA_1B_1,&
 X_{12}&=U^2A_2B_1^2,\\
 X_{21}&=U^2B_2A_1^2,&
 X_{22}&=U^2A_2B_2.
\end{aligned}                                    \tag{9.1}
\]

They obey the toric cycle relation

\[
 \boxed{X_{12}X_{21}=X_{11}^2X_{22}.}             \tag{9.2}
\]

The first constant terms are already sparse:

\[
\begin{aligned}
 \operatorname{CT}_T(P^2)&=2(X_{11}+X_{22}),\\
 \operatorname{CT}_T(P^3)&=3(X_{12}+X_{21}).
\end{aligned}                                    \tag{9.3}
\]

Thus a first-cycle search should use the invariant coordinates (9.1) and
the single toric relation (9.2), rather than the four coefficient
polynomials \(A_1,A_2,B_1,B_2\) with their redundant scalings.  Equation
(9.2) is the algebraic feature absent from the three-level forest and is
the natural location for genuinely new cancellation.  A useful finite
experiment is to enumerate the constant-term polynomials in the
\(X_{ij}\), reduce them modulo (9.2), apply \(\mathcal L\), and only then
introduce bounded radial coefficients.
