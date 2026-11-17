# Complete all-degree closure of the rank-two length-three `[3,1]` HC4 branch

## Status and scope

**Current status.**  `HC4RSD63` completes the all-degree rank-two `[3,1]`
closure.  Together with `HC4RSD60` it removes both generic rank-two Jordan
types from the live relative-nilpotent frontier. The remaining rank-three
`[4]` branch in `HC4MR1` is closed by the corrected HC4MRA1 argument and
[HC4MRA2](HC4_NEGATIVE_MOTION_POLYNOMIAL_OBSTRUCTION.md).  None of these statements proves
unrestricted `HC4` or `JC2`.

This note completes the branch reduced in `HC4RSD61--62`.

Let

\[
S=\operatorname{Hess}\psi,\qquad
T=\operatorname{Hess}A,\qquad
N=S^{-1}T,
\]

with

\[
\det S=\delta\in K^*,\qquad
\operatorname{rank}T=2,\qquad
N^3=0,\qquad N^2\ne0
\]

over a characteristic-zero field.

`HC4RSD61` gives constant active coordinates `(x,w)` with

\[
A=A(x,w),\qquad K=\operatorname{Hess}_{x,w}A,
\qquad \det K\ne0,
\]

and passive-kernel coordinates

\[
\psi=\Phi(x,w,h)+a(x,w)z,
\qquad
h=y-q(x,w)z,
\tag{0.1}
\]

where `Phi_hh != 0` on the genuine rank-one passive branch.

`HC4RSD62` further gives

\[
(\nabla q)^T\operatorname{adj}(K)\nabla q=0,
\qquad
da\wedge dq=0.
\tag{0.2}
\]

> **Theorem HC4RSD63 — complete `[3,1]` closure.**
> Every packet above has a constant passive kernel direction.  After a constant
> linear change and affine normalization,
>
> \[
> \boxed{
> \psi=xz+C(x,w,y),
> \qquad
> A=A(x,w),
> \qquad
> \det\operatorname{Hess}_{w,y}C\in K^*.
> }
> \tag{0.3}
>
> Thus the genuine moving `[3,1]` relative-nilpotent branch does not exist.
> The residual packet is an `HC2` suspension and is already closed by the
> two-variable Hessian theorem.

The proof uses two general lemmas which are useful independently of HC4.

## 1. Binary bordered-Hessian lemma

For `f in K[x,w]` define

\[
\mathcal B_2(f)
:=\nabla f^T\operatorname{adj}(\operatorname{Hess}f)\nabla f.
\tag{1.1}
\]

Explicitly,

\[
\mathcal B_2(f)
=f_w^2f_{xx}-2f_xf_wf_{xw}+f_x^2f_{ww}.
\tag{1.2}
\]

> **Lemma.** If `f` is a nonconstant polynomial over a characteristic-zero
> field and `B_2(f)=0`, then
>
> \[
> f=P(\ell)
> \tag{1.3}
> \]
>
> for a univariate polynomial `P` and an affine-linear form `ell`.

### Proof

After extending scalars to an algebraic closure, write

\[
f=P(g)
\]

with `g` non-composite.  The bordered expression scales as

\[
\mathcal B_2(P(g))=P'(g)^3\mathcal B_2(g),
\tag{1.4}
\]

so `B_2(g)=0`.

By Bertini--Krull, a generic fiber `g=lambda` is irreducible.  At a smooth
point its tangent vector can be chosen as

\[
X_g=(g_w,-g_x).
\]

Equation `B_2(g)=0` is exactly

\[
X_g^T(\operatorname{Hess}g)X_g=0,
\]

so the curvature of the generic level curve is zero.  Hence every generic
irreducible fiber is an affine line.

Two distinct fibers of a polynomial are disjoint.  Two nonparallel affine
lines over an algebraically closed field intersect, so all generic fiber lines
are parallel.  Therefore a constant vector `v != 0` is tangent to every generic
fiber.  The polynomial `D_v g` vanishes on infinitely many fibers and hence is
identically zero.  Thus `g` depends on one affine-linear coordinate.  Since `g`
is non-composite, it is affine-linear itself.  This proves (1.3).

The conclusion descends to the original field because the space of constant
annihilating vectors of `f` is defined by linear equations over `K`.

## 2. The full determinant supplies a binary bordered pencil

Return to (0.1).  Put

\[
\lambda=\Phi_h,
\qquad
b_\lambda=a-\lambda q.
\tag{2.1}
\]

A direct Hessian calculation shows that the coefficient of `z` in the full
four-variable determinant is

\[
\boxed{
[z]\det\operatorname{Hess}\psi
=-\Phi_{hh}\,\mathcal B_2(b_\lambda).
}
\tag{2.2}
\]

Here `lambda` is treated as a scalar when the active derivatives in
`B_2(a-lambda q)` are formed.  Since `det Hess psi=delta` is independent of
`z` and `Phi_hh != 0`,

\[
\mathcal B_2(a-\lambda q)=0
\tag{2.3}
\]

for the nonconstant value `lambda=Phi_h`.  As the left side is a cubic
polynomial in the indeterminate `lambda`, and `Phi_h` is transcendental over
the active fraction field unless it is algebraic there (in which case
`Phi_hh=0`), all coefficients vanish.  Equivalently (2.3) holds identically in
`lambda`.

In particular

\[
\mathcal B_2(q)=0,
\qquad
\mathcal B_2(a)=0.
\tag{2.4}
\]

Combined with `da wedge dq=0` from `HC4RSD62`, the binary lemma shows that, if
`q` is nonconstant, there is one constant affine-linear active coordinate
`ell` such that

\[
q=Q(\ell),\qquad a=R(\ell).
\tag{2.5}
\]

If `q` is constant, (2.5) is already true after choosing `ell` from `a` when
needed.  Hence after a constant active change we may assume

\[
q=q(x),\qquad a=a(x).
\tag{2.6}
\]

## 3. Constant determinant kills the last moving ratio

Substitute (2.6) into the complete determinant, not merely its `z` coefficient.
All terms collapse to

\[
\boxed{
\det\operatorname{Hess}\psi
=-\bigl(a'(x)-q'(x)\Phi_h\bigr)^2
\bigl(\Phi_{ww}\Phi_{hh}-\Phi_{wh}^2\bigr).
}
\tag{3.1}
\]

The left side is the nonzero constant `delta`.  Both factors on the right are
polynomials.  Therefore each factor is a unit.

But `Phi_hh != 0`, so `Phi_h` is nonconstant in `h`.  If `q'(x) != 0`, the
first factor in (3.1) is nonconstant.  Hence

\[
q'(x)=0.
\tag{3.2}
\]

Thus the passive kernel line `(q,1)` is constant.  Equation (3.1) then gives

\[
a'(x)\in K^*,
\qquad
\det\operatorname{Hess}_{w,h}\Phi\in K^*.
\tag{3.3}
\]

After affine rescaling and the constant passive change `h=y-qz`,

\[
\psi=xz+C(x,w,y)
\]

with

\[
\det\operatorname{Hess}_{w,y}C\in K^*.
\]

This is exactly the `HC2` suspension in (0.3), proving `HC4RSD63`.

## 4. Consequence for the moving-frame frontier

Together with `HC4RSD60`, both generic rank-two nilpotent types are now closed
in every degree:

\[
\boxed{
\operatorname{rank}T=2
\Longrightarrow
[2,2]\text{ is }JC2/HC2,
\quad
[3,1]\text{ is }HC2.
}
\]

The only relative-nilpotent moving-chain stratum not covered by the rank-one
and rank-two reductions is therefore generic rank three, i.e. the length-four
`[4]` type.

## 5. Verification

The companion checker `scripts/verify_hc4_rank_two_length_three_closure.py`
verifies both universal symbolic identities (2.2) and (3.1), as well as the
composition scaling (1.4).
