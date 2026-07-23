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
