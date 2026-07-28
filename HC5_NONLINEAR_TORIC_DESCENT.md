# Nonlinear toric descent from the Meng--Yang `HC_5` example

## Status

This note records an exact nonlinear one-variable reduction of the
Meng--Yang five-variable Hessian counterexample.  It succeeds in producing a
polynomial unit pivot and, after a relative toric correction, a
four-variable potential with constant Hessian determinant.  In the same
class, however, constant determinant forces the descended gradient to be a
polynomial automorphism, so the Meng--Yang collision cannot survive.

This is a restricted obstruction.  It does not settle `HC_4` and does not
exclude non-toric symplectic changes, non-coordinate coisotropic embeddings,
or higher-degree critical equations.

## 1. A polynomial unit-pivot coordinate

Put \(u=1+x_1x_2\) and write the Meng--Yang potential as

\[
 \Psi=A^2+11A+2B,
\]

where

\[
 A=u^3y_1+3x_1u^2y_2-x_1^3y_3
\]

and \(B\) is the linear form in the three \(y\)-variables displayed in the
Meng--Yang paper.  Set \(v=x_1x_2\) and

\[
\begin{aligned}
 p(v)&=1-3v+6v^2,\\
 q(x_1,x_2)&=x_2^3(6v^2+15v+10).
\end{aligned}
\]

The Bezout identity

\[
 p(v)u^3-x_1^3q=1
\]

gives the determinant-one matrix

\[
 U=
 \begin{pmatrix}
 u^3&3x_1u^2&-x_1^3\\
 0&1&0\\
 -q&0&p
 \end{pmatrix}.
\]

Define \((t,r,s)^T=U(y_1,y_2,y_3)^T\).  Then \(t=A\), and the transformed
potential is monic quadratic in \(t\).  Consequently the critical equation

\[
 \partial_t\Phi=\sigma
\]

has a polynomial solution for every constant \(\sigma\).

The two Meng--Yang collision points lie on the same critical level,

\[
 \sigma=-\frac{19}{2}.
\]

Thus nonlinear polynomial solvability of the critical equation is not the
obstruction.

## 2. The natural Schur complement

In the above completion, the coefficients of \(r,s,t\) in \(B\) are

\[
\begin{aligned}
 \beta_r&=-x_2P(v),\\
 \beta_s&=x_1Q(v),\\
 \beta_t&=x_2^2(6v^2+15v+4),
\end{aligned}
\]

where

\[
\begin{aligned}
 P(v)&=18v^5+81v^4+120v^3+60v^2-1,\\
 Q(v)&=(v+1)(v+2).
\end{aligned}
\]

After eliminating \(t\), the four-variable potential is linear in \(r,s\).
Its Hessian determinant is

\[
 16J(v)^2,
\]

with

\[
 J(v)=144v^7+945v^6+2394v^5+2910v^4
      +1680v^3+357v^2-6v-2.
\]

This is nonconstant.  Moreover, direct evaluation at the two reduced
Meng--Yang points gives different gradients.  A nonlinear point
transformation preserves the Lagrangian graph but does not automatically
preserve equality of its vertical projections.

## 3. Exact principal-part cancellation

There is nevertheless an explicit relative determinant-one correction.
Let

\[
 C=
 \begin{pmatrix}
 a(v)&x_1^2b(v)\\
 x_2^2c(v)&d(v)
 \end{pmatrix},
\]

where

\[
\begin{aligned}
 a(v)&=(3v^2+6v+2)/2,\\
 b(v)&=2v+3,\\
 c(v)&=3(18v^4+63v^3+69v^2+21v-1)/2,\\
 d(v)&=36v^5+108v^4+87v^3+3v^2-3v+1.
\end{aligned}
\]

Exact expansion gives

\[
 \det C=1,\qquad
 \bigl(-x_2P(v),\,x_1Q(v)\bigr)C=(x_2,\,2x_1).
\]

Replacing \((r,s)^T\) by \(C(r,s)^T\) before eliminating \(t\) cancels the
entire nonconstant principal part.  The descended potential becomes

\[
 \psi_\sigma(x_1,x_2,r,s)
 =f_\sigma(x_1,x_2)+2x_2r+4x_1s
\]

and satisfies

\[
 \det\operatorname{Hess}\psi_\sigma=64.
\]

But its gradient is polynomially invertible.  If its output coordinates are
\((p_1,p_2,p_r,p_s)\), then

\[
\begin{aligned}
 x_1&=p_s/4,&x_2&=p_r/2,\\
 r&=(p_2-\partial_{x_2}f_\sigma)/2,&
 s&=(p_1-\partial_{x_1}f_\sigma)/4.
\end{aligned}
\]

Hence the correction that restores constant determinant necessarily loses
the known collision.

## 4. All-degree toric obstruction

The preceding cancellation is a special case of a simple exact identity.
For arbitrary \(P,Q\in k[v]\), consider

\[
 G(x_1,x_2)=\bigl(-x_2P(v),\,x_1Q(v)\bigr),
 \qquad v=x_1x_2.
\]

Then

\[
 \det DG=\frac{d}{dv}\bigl(vP(v)Q(v)\bigr).
\]

If this determinant is a nonzero constant \(\kappa\), integration in
characteristic zero and evaluation at \(v=0\) give

\[
 vP(v)Q(v)=\kappa v,
\qquad P(v)Q(v)=\kappa.
\]

Since \(k[v]\) is a domain, both \(P\) and \(Q\) are units.  Thus every
constant-determinant member of this toric radial class is linear up to
nonzero scalings, and its doubled four-variable gradient is a polynomial
automorphism.  No collision can survive anywhere in this class.

The missing search space is therefore genuinely non-toric: a successful
descent must use a symplectic change whose reduced two-variable coefficient
map is not radial, or leave the coordinate-graph/linear-in-two-duals model
altogether.  The first bounded non-toric class is treated next.

## 5. Bounded non-toric relative corrections

Keep the natural complementary coefficient row

\[
\beta(x,y)=\bigl(-yP(xy),\,xQ(xy)\bigr)                \tag{5.1}
\]

from Section 2, and let

\[
C(x,y)\in \operatorname{SL}_2(\mathbb Q[x,y]),\qquad
G=\beta C.                                             \tag{5.2}
\]

After the unit-pivot elimination, the corrected four-variable potential has
the form

\[
\psi=f(x,y)+2G_1(x,y)r+2G_2(x,y)s.                    \tag{5.3}
\]

For an arbitrary upper-left Hessian block, the block determinant identity is

\[
\boxed{\det\operatorname{Hess}\psi=16(\det DG)^2.}     \tag{5.4}
\]

Thus constant nonzero Hessian determinant requires \(\det DG\) to be a
nonzero constant.  The known collision bases are

\[
p_+=(1,-3/2),\qquad p_-=(-1,3/2),                     \tag{5.5}
\]

so retaining the collision additionally requires \(G(p_+)=G(p_-)\).

Take all four entries of \(C\) to be arbitrary polynomials of total degree
at most four.  This gives sixty coefficients.  Exact coefficient extraction
produces:

- 45 equations from \(\det C=1\);
- 218 equations from the nonconstant spatial coefficients of \(\det DG\).

The resulting 263-equation ideal over \(\mathbb Q\) has Gröbner basis
\(\{1\}\).  In particular, this class fails before the two collision
equations are imposed:

> **Degree-four non-toric correction obstruction.**  No polynomial
> relative \(SL_2\) correction whose entries have total degree at most four
> can make the descended Hessian determinant a nonzero constant.

There is a complementary perturbative test at the known degree-ten toric
correction \(C_0\), for which \(\beta C_0=(y,2x)\).  Appending an arbitrary
affine \(D(x,y)\in SL_2\) gives twelve parameters.  The determinant-one,
constant-Jacobian, and two collision equations again generate the unit
ideal.  Hence no affine non-toric perturbation of \(C_0\) restores the
collision while retaining constant determinant.

These are bounded obstructions, not an all-degree theorem.  Quadratic and
higher perturbations of \(C_0\), raw corrections of degree at least five,
general nonlinear changes mixing base and dual variables, and non-coordinate
coisotropic embeddings remain open.  The complete homogeneous-quartic mixed
pivot class, including simultaneous pure-base and pure-dual sectors, is
excluded separately by theorem `HC4HQ1` in
[the sparse-quartic audit](HC4_MENG_SPARSE_QUARTIC_AUDIT.md).  The first two
open relative-correction classes already contain genuine
constant-determinant solutions, so their Gröbner systems are substantially
harder than the degree-four raw correction ideal.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_hc5_nonlinear_toric_descent.py
.venv/bin/python scripts/verify_hc4_nontoric_sl2_correction_degree4.py
```
