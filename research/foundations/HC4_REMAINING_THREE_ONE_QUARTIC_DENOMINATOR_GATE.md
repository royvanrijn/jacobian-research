# The remaining three-plus-one quartic-denominator gate

## Status

This note proves `HC4NHM6`: the two clean generic-corank-one incidence
packets

\[
P=x^3y,qquad
\det\operatorname{Hess}(h_5)=x^6y^3
\quad\text{or}\quad
x^6y^2z
\]

are empty. Together with `HC4NHM4--5`, this closes every residual-line
incidence in the quartic-denominator partition \(3+1\).

The exclusion occurs at the ternary Hessian boundary, before the Schur
gradient, lower four-variable faces, or collision equations. It does not
treat the partitions \(2+2\), \(2+1+1\), \(1+1+1+1\), positive-defect
packets, lower-Smith components, or unrestricted `HC(4)`.

Replay the exact determinant identities with

~~~bash
.venv/bin/python scripts/verify_hc4_remaining_three_one_quartic_denominator_gate.py
~~~

## 1. Every essential-line kernel is constant

Let

\[
C=\operatorname{Hess}(h_5),\qquad d=\nabla s_3,
\tag{1.1}
\]

and suppose the minimal denominator clearing \(C^{-1}d\) is \(P=x^3y\).
The global cleared vector \(e=PC^{-1}d\) is cubic. On the normalization of
the essential line \(x=0\), its saturated kernel degree \(\kappa\) therefore
satisfies \(0\le\kappa\le3\).

The moving-line and conic arguments of `HC4NHM2` exclude
\(\kappa=1,2\) for any repeated Hessian line. The perfect-square binary
Hessian argument of `HC4NHM3` excludes the defect-free cubic case
\(\kappa=3\). Hence

\[
\kappa=0.
\tag{1.2}
\]

Thus the kernel of \(C|_{x=0}\) is constant.

## 2. The transverse constant kernel is impossible

If the kernel is transverse to \(x=0\), normalize it to \(\partial_x\).
The boundary equations give \(G=H=0\) in the factorial expansion

\[
h_5=F+xG+\frac{x^2}{2}H+\frac{x^3}{6}J
       +\frac{x^4}{24}K+\frac{x^5}{120}c.
\tag{2.1}
\]

Generic boundary rank says the binary Hessian determinant of \(F\) is
nonzero. The first three normal determinant coefficients are successively

\[
J\det\operatorname{Hess}(F),\qquad
\frac K2\det\operatorname{Hess}(F),\qquad
\frac c6\det\operatorname{Hess}(F).
\tag{2.2}
\]

Divisibility by \(x^6\) kills \(J,K,c\), after which the determinant is
zero. Therefore the constant kernel must be tangent to \(x=0\).

## 3. The tangent coefficient ladder

Normalize the tangent kernel to \(\partial_z\). Along \(x=0\), the boundary
equations give

\[
F=\alpha y^5,qquad G=\beta y^4.
\tag{3.1}
\]

Write

\[
\begin{aligned}
H&=h_0y^3+h_1y^2z+h_2yz^2+h_3z^3,\\
J&=j_0y^2+j_1yz+j_2z^2,\\
K&=k_0y+k_1z.
\end{aligned}
\tag{3.2}
\]

### 3.1 The \(\alpha\ne0\) branch

Divisibility by \(x^6\) forces

\[
h_1=h_2=h_3=j_1=j_2=0.
\tag{3.3}
\]

The boundary matrix has generic rank two exactly when

\[
5\alpha h_0-4\beta^2\ne0.
\tag{3.4}
\]

The complete remaining quintic is

\[
h_5=\alpha y^5+\beta xy^4+\frac{h_0}{2}x^2y^3
 +\frac{j_0}{6}x^3y^2+\frac{x^4}{24}(k_0y+k_1z)
 +\frac c{120}x^5,
\tag{3.5}
\]

and its determinant factors exactly as

\[
\boxed{
\det C=-\frac{k_1^2x^6}{108}
\left(60\alpha y^3+36\beta xy^2+9h_0x^2y+j_0x^3\right).
}
\tag{3.6}
\]

For exact \(x\)-multiplicity six, \(k_1\ne0\). Equation (3.6) shows that
the residual cubic on \(x=0\) is necessarily a cube. Hence the root
partition \(2+1\), represented by \(y^2z\), is impossible. This excludes

\[
\det C=x^6y^2z.
\tag{3.7}
\]

For the residual cube \(y^3\), exact equality in (3.6) forces

\[
\beta=h_0=j_0=0.
\tag{3.8}
\]

But then the generic-rank condition (3.4) becomes \(0\ne0\). Thus

\[
\det C=x^6y^3
\tag{3.9}
\]

is impossible as well.

### 3.2 The \(\alpha=0\) branch

Generic boundary rank forces \(\beta\ne0\). The earlier tangent ladder gives

\[
h_1=h_2=h_3=j_2=0,
\tag{3.10}
\]

and the \(x^5\)-coefficient is

\[
-\frac73\beta j_1^2y^4.
\tag{3.11}
\]

Divisibility by \(x^6\) therefore forces \(j_1=0\). The complete determinant
then starts in order at least seven:

\[
\det C=-\frac{k_1^2x^7}{108}
\left(36\beta y^2+9h_0xy+j_0x^2\right).
\tag{3.12}
\]

So this branch cannot have exact multiplicity six.

Combining the constant-kernel alternatives proves:

> **Theorem `HC4NHM6` -- Remaining three-plus-one exclusion.** On the clean
> generic-corank-one Schur branch with minimal denominator \(P=x^3y\), no
> Hessian boundary with determinant \(x^6y^3\) or \(x^6y^2z\) exists.

With `HC4NHM5`, every residual-line incidence for \(P=x^3y\) is now closed:
the \(x^7y^2\) incidence has a nonzero module section but no four-variable
prolongation, while the other two incidences are empty already as Hessian
boundaries.

## 4. Subsequent partitions

The next clean quartic denominator is

\[
P=x^2y^2,
\tag{4.1}
\]

with residual-line incidence types

\[
\det C=x^5y^4,qquad x^4y^5,qquad x^4y^4z.
\tag{4.2}
\]

Exchanging \(x\) and \(y\) identifies the first two, so there are two
geometric cases: a residual line coincident with a repeated component and a
third distinct residual line. Both are excluded in
[`HC4_TWO_TWO_QUARTIC_DENOMINATOR_GATE.md`](HC4_TWO_TWO_QUARTIC_DENOMINATOR_GATE.md).
The \(2+1+1\) partition, including every concurrency stratum, is excluded in
[`HC4_TWO_ONE_ONE_QUARTIC_DENOMINATOR_GATE.md`](HC4_TWO_ONE_ONE_QUARTIC_DENOMINATOR_GATE.md).
The last clean partition is reduced to a finite constant-polar flag problem
in
[`HC4_SQUAREFREE_QUARTIC_DENOMINATOR_FRONTEND.md`](HC4_SQUAREFREE_QUARTIC_DENOMINATOR_FRONTEND.md).
