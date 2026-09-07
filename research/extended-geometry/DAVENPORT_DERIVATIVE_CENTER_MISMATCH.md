# Davenport derivative-center mismatch

The split-boundary weighted gluing problem produces a smooth determinant
threefold

\[
X_\beta=\{D^2U-CR=\beta\}.
\]

Forgetting \(U\) is an affine modification of \(\mathbb A^3_{D,C,R}\)
whose residue Jacobian is \(D^2\).  This doubled exponent exactly matches
the \(J^2\) column in the Davenport branch-pullback ledger.

The match stops at the exponent.  The two modification centers have
different normalized affine geometry and cannot be identified, even after
polynomial stabilization.

Work over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

## 1. The reverse determinant modification

The coordinate ring of \(X_\beta\) is

\[
K[D,C,R]\left[\frac{CR+\beta}{D^2}\right].
\]

Thus

\[
\pi:X_\beta\longrightarrow\mathbb A^3_{D,C,R}
\]

is the affine modification with denominator \(D^2\) and center ideal

\[
(D^2,CR+\beta).
\]

On the hypersurface, the residue form is

\[
\Omega_{X_\beta}
=\frac{dD\wedge dC\wedge dR}{D^2},
\]

so

\[
\pi^*(dD\wedge dC\wedge dR)=D^2\Omega_{X_\beta}.     \tag{1}
\]

The reduced center is

\[
D=0,\qquad CR=-\beta,
\]

hence isomorphic to \(\mathbb G_m\).  Its geometric unit rank is one.

## 2. The Davenport derivative curve

Let

\[
J(T,Y)=g_T'(Y).
\]

Its Newton polygon is the triangle

\[
\operatorname{conv}\{(2,0),(3,0),(0,6)\}.
\]

It has exactly one interior lattice point, \((2,1)\), so its toric
arithmetic genus is one.

The curve has exactly one singularity in the two-dimensional torus:

\[
(T,Y)=\left(-\frac14,-\frac12\right).
\]

After translating \(T=-1/4+t\), \(Y=-1/2+y\), its tangent cone is

\[
-\frac t4\left((-3a-6)t+(2a+8)y\right).
\]

The two factors are distinct, so this is an ordinary node with
\(\delta=1\).  All Newton faces are squarefree.  Therefore the
normalization has genus

\[
1-1=0.
\]

The outer edge joins \((3,0)\) to \((0,6)\) and has lattice length three.
Its face polynomial is a squarefree cubic, so the affine normalization is

\[
\mathbb P^1\setminus\{P_1,P_2,P_3\}
\]

over an algebraic closure.  Consequently

\[
\operatorname{rank}
\left(
\mathcal O(\widetilde{J=0})^*/\overline K^*
\right)=2.                                           \tag{2}
\]

## 3. Why the doubled ledgers cannot be identified directly

The exponent-two agreement is

\[
\begin{array}{c|c}
\text{Davenport branch pullback}&J^2\\
\text{reverse determinant modification}&D^2.
\end{array}
\]

But their reduced centers have geometric unit ranks

\[
2\qquad\text{and}\qquad1.                            \tag{3}
\]

Polynomial extension preserves units, so this mismatch survives adjoining
arbitrarily many affine variables.  In particular, no center-preserving
stable affine identification can pair the Davenport derivative boundary
with the center \(CR=-\beta\).

This explains why the split weighted modification and the Cox ledger meet
at the same multiplicity without actually solving one another.

## 4. Consequence for the remaining construction

The reverse modification remains useful: equation (1) supplies exactly the
doubled Jacobian factor needed by the Cox ledger.  A successful
construction must, however, change its center from a two-puncture rational
curve to the three-puncture normalization of \(J=0\), or merge one of the
three punctures through a nonsymmetric quotient.

The canonical puncture-merging quotient is constructed in the
[normalized-boundary involution audit](DAVENPORT_BOUNDARY_INVOLUTION.md).
It produces \(\mathbb G_m\) on the normalization but fails to preserve the
two-point node conductor, so it does not descend to the actual derivative
curve.

The latter operation is delicate: the global Sunada construction must
retain the nonconjugate point and line stabilizers, so an affine descent
cannot simply forget the boundary label that distinguishes the two covers.

## 5. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_derivative_center_mismatch.py
```

The checker verifies the Newton polygon, unique node, squarefree outer face,
three punctures, the residue Jacobian \(D^2\), and the stable unit-rank
mismatch of the two centers.
