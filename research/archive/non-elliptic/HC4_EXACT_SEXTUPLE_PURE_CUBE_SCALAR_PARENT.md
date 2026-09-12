# Exact-sextuple pure-cube scalar-parent closure

## Status

This note proves `HC4DIR28`.  It closes the sole degree-five, order-one,
lower-Smith exact-sextuple resonance left by `HC4DIR27`.  It is a
collision exclusion for that scalar-parent family, not an unrestricted
`HC(4)` theorem.

Replay the universal determinant and gradient identities with

```bash
.venv/bin/python scripts/verify_hc4_exact_sextuple_pure_cube_scalar_parent.py
```

The critical-point lemma below is a written linear-algebra proof over the
algebraic closure.  The checker verifies its two scalar equations and the
complete Schur/graph-coordinate identities.

## 1. The remaining resonance

`HC4DIR26--27` reduce the lower-rank exact-sextuple packet to

\[
 f=y^5+cx^3y^2+x^4(uy+vz)+dx^5,
 \qquad cv\ne0,
\]

with first moving direction \(x\partial_z\).  Every four-variable
completion has the scalar-parent form

\[
 \Psi(w,X)=H(X)+wP(X)+\frac{\eta}{2}w^2,
 \qquad X=(x,y,z),
\tag{1.1}
\]

where

\[
 P=P_3+P_{\le2},
 \qquad P_3=\frac{4v}{3}x^3.
\tag{1.2}
\]

Thus it is enough to treat scalar parents whose cubic pivot has a pure-cube
leading form.

## 2. A pure-cube gradient lemma

We use the following degree-three fact.

> **Lemma.** Let \(K\) be algebraically closed of characteristic zero and
> let
> \[
> P(X)=a\ell(X)^3+\frac12A(X,X)+b(X)+c,
> \qquad a\ne0,
> \tag{2.1}
> \]
> where \(A:V\to V^*\) is symmetric.  If \(\nabla P\) has no zero, then
> there is a constant vector \(u\) such that
> \[
> \ell(u)=0,\qquad Au=0,\qquad b(u)\ne0.
> \tag{2.2}
> \]
> In particular \(D_uP=b(u)\in K^\times\).

Put \(N=\ker A\) and \(W=N\cap\ker\ell\).  Suppose conversely that
\(b|_W=0\).

If \(\ell|_N\ne0\), the two linear forms \(b|_N\) and \(\ell|_N\) have
the same kernel containment, so

\[
 b|_N=\beta\ell|_N.
\]

Choose \(t\) with \(\beta+3at^2=0\).  Then
\(b+3at^2\ell\) annihilates \(N\), hence belongs to
\(\operatorname{im}A\).  Choose \(X_0\) satisfying

\[
 AX_0=-b-3at^2\ell.
\]

Adding an element of \(N\) adjusts \(\ell(X_0)\) arbitrarily, so it can be
made equal to \(t\).  The resulting \(X\) is a critical point of \(P\).

If \(\ell|_N=0\), then \(W=N\), and the assumption gives
\(b,\ell\in\operatorname{im}A\).  Choose \(U,V\) with

\[
 AU=-b,
 \qquad AV=-3a\ell.
\]

The vector \(X=U+t^2V\) solves the vector part of the critical equation.
Its remaining condition \(\ell(X)=t\) is

\[
 \ell(U)+t^2\ell(V)-t=0,
\]

a nonconstant polynomial over an algebraically closed field.  It has a
root.  This again produces a critical point.  The contrapositive proves the
lemma.

The argument is invariant under scalar extension.  Therefore if the
gradient ideal of (2.1) is the unit ideal over the original characteristic-
zero field, the vector in (2.2) is already detected after a constant linear
extension, which is enough for collision exclusion.

## 3. The nonzero corner

Assume \(\eta\ne0\) and put

\[
 \psi=H-\frac{P^2}{2\eta},
 \qquad s=w+\frac P\eta.
\]

The Schur identity gives

\[
 \det\operatorname{Hess}\Psi
 =\eta\det\operatorname{Hess}_X(\psi+sP).
\tag{3.1}
\]

Since the left side is a nonzero constant, every member of the ternary
pencil \(\psi+sP\) has nonzero constant Hessian determinant.  Equality of
two parent gradients first fixes \(s\), then becomes equality of the two
ternary pencil gradients.  `HC3` makes those injective.  Hence the nonzero
corner has no collision.

## 4. The zero corner

Now let \(\eta=0\).  With \(g=\nabla P\) and
\(M=\operatorname{Hess}(H+wP)\),

\[
 \operatorname{Hess}\Psi=
 \begin{pmatrix}0&g^{\mathsf T}\\g&M\end{pmatrix}.
\tag{4.1}
\]

If \(g\) vanished at any point, (4.1) would have zero first row and column,
contrary to its nonzero constant determinant.  Thus \(\nabla P\) is
nowhere zero.  The lemma supplies a constant direction \(u\) with
\(D_uP\in K^\times\).

Normalize coordinates to

\[
 X=(U,r_0),\qquad P=r_0+q(U),\qquad U=(u_1,u_2).
\]

Set \(r=P\), write

\[
 C(U,r)=H(U,r-q(U)),
 \qquad \tau=w+C_r(U,r).
\]

The parent gradient becomes exactly

\[
 \left(r,\ \nabla_U(C+\tau q),\ \tau\right),
\tag{4.2}
\]

and

\[
 \det\operatorname{Hess}\Psi
 =-\det\operatorname{Hess}_U(C+\tau q).
\tag{4.3}
\]

For fixed output \((r,\tau)\), the middle map is a binary
constant-Hessian gradient.  `HC2` makes it injective; then \(r_0\) and
\(w\) are recovered from the two coordinate changes.  Hence the zero
corner also has no collision.

## 5. Result

> **Theorem `HC4DIR28` -- Exact-sextuple pure-cube scalar-parent closure.**
> Every four-variable constant-Hessian completion of the degree-five
> order-one exact-sextuple resonance of `HC4DIR26` is collision-free.  The
> nonzero-corner branch reduces to `HC3`; the zero-corner branch has a
> constant unit pivot direction and reduces to `HC2`.

The next rank-three top-Hessian targets are therefore nonlinear repeated
factors, generic lower-Smith components not lying in this scalar-parent
resonance, or interactions among several repeated factors.  Leading-Hessian
rank one and two remain separate synchronization frontiers.
