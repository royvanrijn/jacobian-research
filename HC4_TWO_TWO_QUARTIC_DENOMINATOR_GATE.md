# The two-plus-two quartic-denominator gate

## Status

This note proves `HC4NHM7`: every clean generic-corank-one packet with

\[
P=x^2y^2,\qquad \det\operatorname{Hess}(h_5)=P^2\ell
\]

is empty. Up to exchanging \(x,y\), the residual-line incidences are
\(x^5y^4\) and \(x^4y^4z\). Both fail in the tangent constant-kernel
determinant ladder before Schur or four-variable prolongation equations.

Replay the exact identities with

~~~bash
.venv/bin/python scripts/verify_hc4_two_two_quartic_denominator_gate.py
~~~

## 1. Constant-kernel reduction

The global cleared vector has degree three on either essential repeated line.
As in `HC4NHM6`, `HC4NHM2--3` exclude every positive saturated kernel
degree. The kernel on \(x=0\) is constant. A transverse constant kernel makes
the determinant zero under \(x^4\)-divisibility, so normalize the kernel to
the tangent direction \(\partial_z\).

Use

\[
F=\alpha y^5,\quad G=\beta y^4,\quad
H=h_0y^3+h_1y^2z+h_2yz^2+h_3z^3,
\tag{1.1}
\]

and

\[
J=j_0y^2+j_1yz+j_2z^2,qquad K=k_0y+k_1z.
\tag{1.2}
\]

For \(\alpha\ne0\), exact multiplicity four has residual fifth power
\(-5\alpha j_1^2y^5\), while killing that face raises the first possible
order to six. Thus this branch matches neither exact multiplicity five nor
the \(4+1\) residual of \(x^4y^4z\).

It remains to treat \(\alpha=0\), where generic boundary rank gives
\(\beta\ne0\).

## 2. The coincident residual line: \(x^5y^4\)

Divisibility by \(x^5\) forces

\[
h_1=h_2=h_3=j_2=0.
\tag{2.1}
\]

Exact multiplicity five requires \(j_1\ne0\), and

\[
[x^5]\det C=-\frac73\beta j_1^2y^4.
\tag{2.2}
\]

The complete next coefficients include

\[
[x^7yz]\det C=\frac1{18}j_1^3.
\tag{2.3}
\]

This is nonzero and independent of every remaining coefficient. But
\(x^5y^4\) has no \(x^7yz\)-term. Therefore this incidence is empty.

## 3. The distinct residual line: \(x^4y^4z\)

Exact multiplicity four on the \(\alpha=0\) branch requires \(h_1\ne0\) and

\[
j_2=-\frac{3h_1^2}{4\beta}.
\tag{3.1}
\]

The first coefficient is

\[
[x^4]\det C
=-\frac53h_1y^4
\left((4\beta j_1-3h_0h_1)y-9h_1^2z\right).
\tag{3.2}
\]

The linear factor is automatically distinct from \(y\). A shear of \(z\)
by a linear combination of \(x,y\), preserving both repeated lines and the
kernel direction, normalizes the complete residual line to \(z=0\). In this
chart

\[
4\beta j_1-3h_0h_1=0.
\tag{3.3}
\]

The next determinant coefficient has the immutable term

\[
\boxed{
[x^5y^2z^2]\det C=-\frac{9h_1^4}{2\beta}.
}
\tag{3.4}
\]

It cannot occur in \(x^4y^4z\), and is nonzero because
\(\beta h_1\ne0\). Hence the distinct residual-line incidence is empty.

> **Theorem `HC4NHM7` -- Two-plus-two exclusion.** No clean
> generic-corank-one Hessian--Schur packet with denominator partition
> \(P=x^2y^2\) exists.

The proof uses no gradient, lower-face, or collision equation. The next clean
partition is \(2+1+1\).
