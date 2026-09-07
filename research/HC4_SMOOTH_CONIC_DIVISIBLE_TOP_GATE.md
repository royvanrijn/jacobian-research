# Smooth-conic divisible-top gate

## Status

This note proves `HC4NHM13`, the first exact row beyond split-linear clean
quartic denominators.  It does not close the complete smooth-conic packet.

Let the minimal denominator be

\[
 P=q^2,
\]

where \(q\) is a smooth conic.  The clean normalization theorem `HC4NHM1`
then requires

\[
 \det\operatorname{Hess}(h_5)=q^4\ell.
\tag{0.1}
\]

This note excludes the complete subrow \(h_5\in(q)\).  Consequently every
survivor of (0.1) must restrict nontrivially to the conic.

Replay the exact calculation with

```bash
.venv/bin/python scripts/verify_hc4_smooth_conic_divisible_top_gate.py
```

The checker uses SymPy to construct the coefficient equations and Singular
4.4.1 for exact rational saturation.

## 1. Two residual-line orbits

Over an algebraically closed characteristic-zero field, normalize

\[
 q=xz-y^2.
\tag{1.1}
\]

The projective automorphism group of the conic has two orbits on lines:
tangent lines and secant lines.  Representatives are

\[
 \ell=z
 \qquad\text{and}\qquad
 \ell=y,
\tag{1.2}
\]

respectively.

Every quintic vanishing on the conic has the unique form

\[
 h_5=qG_3,
\tag{1.3}
\]

with \(G_3\) a general ternary cubic.  Write its ten coefficients as
\(g_0,\ldots,g_9\), and introduce a scalar \(k\).  For each representative
in (1.2), let \(I_\ell\) be the coefficient ideal of

\[
 \det\operatorname{Hess}(qG_3)-kq^4\ell.
\tag{1.4}
\]

The exact standard-basis computations give

\[
 \boxed{
 I_z:(k)^\infty=(1),
 \qquad
 I_y:(k)^\infty=(1).
 }
\tag{1.5}
\]

Thus neither orbit has a solution with nonzero Hessian determinant.  The
calculation is an ideal-containment certificate over \(\mathbf Q\), not a
finite-field search or a bounded coefficient sample.

## 2. Deepest support calibration

The most obvious subfamily is \(h_5=q^2L\), with \(L=ax+by+cz\).  Direct
differentiation gives

\[
 \det\operatorname{Hess}(q^2L)
 =16q^3L\,R_2,
\tag{2.1}
\]

where

\[
\begin{aligned}
R_2={}&6a^2x^2+12abxy+(8ac+b^2)xz+(4ac+5b^2)y^2\\
&+12bcyz+6c^2z^2.
\end{aligned}
\tag{2.2}
\]

If \(q\mid R_2\), the \(x^2\) and \(z^2\) coefficients force
\(a=c=0\); then \(R_2=b^2(xz+5y^2)\), which is not a multiple of
\(xz-y^2\) unless \(b=0\).  Hence a nonzero member has conic multiplicity
exactly three, never four.

## 3. Result and next conic row

> **Theorem `HC4NHM13` -- Smooth-conic divisible-top exclusion.**  Let
> \(q\) be a smooth conic and \(h_5\) a ternary quintic divisible by
> \(q\).  Then
> \[
> \det\operatorname{Hess}(h_5)\ne kq^4\ell
> \]
> for every nonzero scalar \(k\) and every line \(\ell\).  Therefore no
> conic-divisible quintic occurs in the clean double-conic minimal-
> denominator packet \(P=q^2\).

The surviving double-conic row has

\[
 h_5|_{q=0}\ne0.
\]

Under the normalization \(q(s^2,st,t^2)=0\), this restriction is a binary
decic.  The next exact target is to stratify that decic by root partition,
solve the fourfold conic-Hessian divisibility on each stratum, and only then
impose the cleared module equations \(Ce=q^2\nabla s_3\) and
\((\nabla s_3)^{\mathsf T}e=q^2a\).  Mixed conic-plus-lines, cubic-plus-line,
and irreducible-quartic denominators remain separate clean packets.
