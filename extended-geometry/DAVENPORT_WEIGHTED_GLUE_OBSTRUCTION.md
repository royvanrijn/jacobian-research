# Davenport weighted-glue obstruction

The two proportional tangent charts cover the full quadratic Davenport
base and their marked-root incidence covers glue by constant root
translations.  The smallest affine modification that tries to lift this
gluing to the weighted Keller targets has a precise failure: its slope and
intercept centers are comaximal, so the combined modification makes the
weighted boundary coordinate \(C\) invertible.

Thus the standard weighted charts cannot be glued by filling \(C=0\).
Their minimal common refinement deletes \(C=0\).

Work over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

## 1. The incidence transition

Let \(\overline H_4\) and \(\overline H_2\) be the two reduced primitives
from the
[proportional-section atlas](DAVENPORT_PROPORTIONAL_TANGENT_SECTIONS.md).
For either Davenport polynomial, their root coordinates differ by a
constant translation \(W_2=W_4+\delta\), and

\[
\overline H_4(W)-\overline H_2(W+\delta)
=m(s)W+n(s).                                         \tag{1}
\]

Both \(m\) and \(n\) are nonzero linear polynomials in \(s\).  Exact
resultants give

\[
\operatorname{Res}_s(m_g,n_g)\ne0,\qquad
\operatorname{Res}_s(m_h,n_h)\ne0.                   \tag{2}
\]

Consequently \((m_\bullet,n_\bullet)=(1)\) in \(K[s]\).

For the incidence equation

\[
H(W)-\sigma W+\tau=0,
\]

equation (1) gives the overlap

\[
\sigma_4=\sigma_2+m(s),\qquad
\tau_4=\tau_2-n(s)-\delta\sigma_2.                   \tag{3}
\]

## 2. The promising slope modification

The weighted target monomializes the incidence parameters as

\[
\sigma=BC,\qquad \tau=cAC^2.                         \tag{4}
\]

The slope transition in (3) can be absorbed by adjoining

\[
V=\frac{m(s)}C,
\qquad CV=m(s),\qquad B_4=B_2+V.                     \tag{5}
\]

This modification is as small as possible, and it is affine space.  Since
\(m(s)=m_1s+m_0\) with \(m_1\ne0\), equation (5) solves

\[
s=\frac{CV-m_0}{m_1}.
\]

After retaining the unused weighted coordinates, the hypersurface
\(CV=m(s)\) is a polynomial affine space of the expected dimension.  Thus
the slope pole by itself is not the obstruction.

## 3. The intercept center forces \(C\) to be a unit

The second transition in (3) simultaneously requires a coordinate \(U\)
with

\[
C^2U=n(s)+\delta BC.                                 \tag{6}
\]

Choose \(u(s),v(s)\in K[s]\) satisfying the Bézout identity

\[
u(s)m(s)+v(s)n(s)=1.                                 \tag{7}
\]

In the coordinate ring defined by (5)--(6), substitute

\[
m=CV,\qquad n=C^2U-\delta BC.
\]

Then (7) becomes

\[
\boxed{
1=C\left(uV+v(CU-\delta B)\right).
}                                                     \tag{8}
\]

Therefore \(C\) is a unit on the full common refinement:

\[
\operatorname{Spec}
\frac{K[s,A,B,C,U,V]}
 {(CV-m,\ C^2U-n-\delta BC)}
\subset D(C).                                        \tag{9}
\]

The combined modification contains no point over \(C=0\).  It reproduces
only the Laurent overlap on which the weighted construction was already
regular.

## 4. Splitting the denominators is still non-affine

One can avoid the immediate unit equation (8) by assigning the slope and
intercept corrections to independent boundary variables.  Replace (4) by
the schematic split

\[
\sigma=BC,\qquad \tau=cAD^2
\]

and impose

\[
CV=m(s),\qquad D^2U=n(s)+\delta BC.                  \tag{10}
\]

Because \(m=m_1s+m_0\), eliminate

\[
s=\frac{CV-m_0}{m_1}.
\]

Writing

\[
n(s)=\alpha CV+\beta,\qquad \beta\ne0,
\]

and changing the free coordinate

\[
R=\alpha V+\delta B
\]

turns (10) into

\[
\boxed{D^2U-CR=\beta.}                               \tag{11}
\]

The unused coordinates are free, so the full candidate target is

\[
\mathbb A^2\times
\{D^2U-CR=\beta\}.                                   \tag{12}
\]

The threefold in (11) is smooth.  Its class is computed by separating
\(D\ne0\) and \(D=0\):

\[
\begin{aligned}
[D\ne0]&=(\mathbb L-1)\mathbb L^2,\\
[D=0]&=(\mathbb L-1)\mathbb L,
\end{aligned}
\]

so

\[
\left[\{D^2U-CR=\beta\}\right]
=\mathbb L^3-\mathbb L.                              \tag{13}
\]

Consequently the full split-boundary candidate has class

\[
\boxed{\mathbb L^5-\mathbb L^3\ne\mathbb L^5.}       \tag{14}
\]

It is not affine five-space, and affine stabilization cannot remove this
motivic defect.  Notably, \(\mathbb L^5-\mathbb L^3\) is the same class
that appears for the
[normalized quadratic--cubic factorization slice](QUADRATIC_CUBIC_FACTORIZATION_SLICE.md).
The two problems reach the same \(SL_2\)-type determinant geometry from
different starting points.  Its Euclidean affine-modification charts
therefore provide the relevant existing atlas for any attempt to continue
this split-boundary route.

There is also a direct comparison with the
[oriented cubic affine descent](ORIENTED_CUBIC_AFFINE_DESCENT.md).  The
natural involution of (11) is

\[
D\longmapsto-D.
\]

Its invariant coordinate is \(X=D^2\), and the quotient is

\[
\boxed{XU-CR=\beta.}                                 \tag{15}
\]

After scaling \(\beta\), this is the standard \(SL_2\)-type determinant
threefold.  It again has class \(\mathbb L^3-\mathbb L\), not
\(\mathbb L^3\).  In the cubic Cox chart, forgetting a factor ordering
produces the affine linear--quadratic slice.  Here forgetting the sign
retains the determinant-threefold defect.  The natural symmetric quotient
therefore does not give affine descent.

## 5. Relation to the controlled-boundary programme

This is the weighted analogue of the extra-divisor phenomenon in the
[controlled-boundary suspension](../cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md):
separate reconstruction requirements cannot automatically be assigned to
one shared boundary variable.  Here the failure is even simpler than the
three-term Jacobian divisor in the two-boundary triangular ansatz.  The two
necessary centers are comaximal before any Jacobian calculation, so their
shared modification denominator becomes invertible.

The result also explains why the
[oriented Davenport Cox map](ORIENTED_DAVENPORT_COX_MAPS.md) is genuinely
different.  Its square-root orientation changes the boundary ledger to
\((1,1,1)\) before attempting affine descent; it does not try to glue two
weighted \(C=0\) charts through the monomials \(BC\) and \(AC^2\).

The doubled exponent in (11) does match the \(J^2\) column of the original
Davenport Cox ledger.  The
[derivative-center audit](DAVENPORT_DERIVATIVE_CENTER_MISMATCH.md)
shows that the match is only numerical: the two reduced centers have stable
geometric unit ranks one and two.

## 6. Remaining route

The following approaches are now ruled out:

1. one global proportional tangent chart;
2. polynomial recoordination of the two standard weighted charts;
3. a common refinement obtained by adjoining the forced slope and
   intercept quotients with denominator \(C\); and
4. the naive split assigning those quotients to two independent boundary
   variables, even after its natural sign quotient.

A surviving construction must change at least one structural input:

1. replace the monomial target chart
   \((A,B,C)\mapsto(BC,cAC^2,C)\);
2. split the corrections only together with an additional nonlinear
   modification that removes the determinant-threefold defect;
3. start from the reduced oriented Davenport ledger and find a
   nonsymmetric affine descent; or
4. use a nonlinear higher-dimensional modification whose center changes
   the old boundary rather than localizing it.

The absolute affine-space Sunada pair remains open, but its weighted-gluing
gap is now reduced to these alternatives.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_weighted_glue_obstruction.py
```

The checker reconstructs both overlap cocycles, verifies the nonzero
resultants, proves that the slope-only hypersurface is affine space, and
certifies the explicit inverse (8) for \(C\) in the full modification ring.
It also reduces the independent-boundary alternative to (11), verifies
smoothness, and computes the class (14).
