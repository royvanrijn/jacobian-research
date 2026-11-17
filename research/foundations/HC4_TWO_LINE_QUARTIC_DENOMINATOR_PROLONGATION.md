# The two-line quartic-denominator prolongation obstruction

## Status

This note proves `HC4NHM5`: no genuinely two-component Hessian--Schur pair
from `HC4NHM4` prolongs to a four-variable constant-Hessian potential. The
obstruction occurs in the first face after Schur cancellation and is
independent of every available quartic, cubic, and quadratic repair term.

This closes the clean incidence

\[
P=x^3y,qquad \det\operatorname{Hess}(h_5)=x^7y^2
\]

inside the nonaligned rank-three direct quintic packet. It does not treat the
other residual-line incidences for \(P=x^3y\), other quartic partitions,
positive-defect packets, lower-Smith components, or unrestricted `HC(4)`.

Replay the exact coefficient calculation with

~~~bash
.venv/bin/python scripts/verify_hc4_two_line_quartic_denominator_prolongation.py
~~~

## 1. The complete first prolongation ansatz

Use variables \((x,y,z,t)\). By `HC4NHM4`, the leading pair is

\[
\begin{aligned}
h_5&=Axy^4+\frac{x^4}{24}(By+\Gamma z)+\frac D{120}x^5,\\
s_3&=axy^2+bx^3,
\end{aligned}
\qquad A\Gamma a\ne0.
\tag{1.1}
\]

The rank-three top-kernel equation makes the quartic term affine in \(t\):

\[
h_4=t\,s_3+r_4(x,y,z),
\tag{1.2}
\]

where \(r_4\) is an arbitrary ternary quartic. The Schur face is

\[
\det(C)\,\partial_t^2h_3
-d^{\mathsf T}\operatorname{adj}(C)d=0.
\tag{1.3}
\]

Since the quotient is \(a^2x/(3A)\), the complete cubic term is

\[
h_3=\frac{a^2}{6A}xt^2+tQ_2(x,y,z)+r_3(x,y,z),
\tag{1.4}
\]

with arbitrary ternary quadratic \(Q_2\) and arbitrary ternary cubic
\(r_3\). A \(t^3\)-term is already excluded by (1.3). Finally let \(h_2\)
be an arbitrary homogeneous quadratic in all four variables. Thus

\[
\psi=h_5+h_4+h_3+h_2
\tag{1.5}
\]

retains every lower term allowed in a collision-normalized degree-five
potential.

## 2. The immutable degree-nine coefficient

Introduce a bookkeeping parameter \(\lambda\) and write

\[
M(\lambda)=
\lambda^3\operatorname{Hess}(h_5)
+\lambda^2\operatorname{Hess}(h_4)
+\lambda\operatorname{Hess}(h_3)
+\operatorname{Hess}(h_2).
\tag{2.1}
\]

The coefficient of \(\lambda^{10}\) in \(\det M(\lambda)\) is exactly the
Schur face (1.3), and vanishes after (1.4). The next face is
\([\lambda^9]\det M\). Exact determinant expansion, with every coefficient
of \(r_4,Q_2,r_3,h_2\) retained, gives

\[
\boxed{
[x^8t]\,[\lambda^9]\det M(\lambda)
=-\frac{\Gamma^2a^3}{54A}.
}
\tag{2.2}
\]

The coefficient is independent of \(B,D,b\) as well. The quadratic term
cannot affect it: the only \(3+3+3+0\) determinant contribution is
\(\det\operatorname{Hess}(h_5)\,(h_2)_{tt}\), which is proportional to
\(x^7y^2\), not \(x^8t\). The checker nevertheless retains the complete
quadratic matrix in the expansion.

Because \(A\Gamma a\ne0\), (2.2) is nonzero. Hence the Hessian determinant of
\(\psi\) cannot be constant.

> **Theorem `HC4NHM5` -- Two-line quartic-denominator prolongation
> obstruction.** Every genuinely two-component packet classified by
> `HC4NHM4` has a nonzero immutable \(x^8t\) coefficient in its degree-nine
> Hessian-determinant face. No choice of the arbitrary terms in (1.2),
> (1.4), or \(h_2\) yields a constant-Hessian four-variable potential.

The argument does not use collision equations. The leading module section
remains a valid and useful calibration: the obstruction says precisely that
its self-dual square cancels the first rank-three Schur face but fails at the
next compatibility face.

## 3. Next clean incidences

The partition \(P=x^3y\) still has two residual-line incidences:

\[
\det C=x^6y^3,qquad \det C=x^6y^2z.
\tag{3.1}
\]

They are both excluded at the ternary Hessian boundary in
[`HC4_REMAINING_THREE_ONE_QUARTIC_DENOMINATOR_GATE.md`](HC4_REMAINING_THREE_ONE_QUARTIC_DENOMINATOR_GATE.md).
Consequently the complete \(3+1\) partition is closed, and the next clean
quartic-denominator partition is \(2+2\).
