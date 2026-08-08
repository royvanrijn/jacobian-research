# The trace–focal dimension-four miracle behind the final HC4 closure

## Status

This note extracts the dimension-dependent linear algebra behind `HC4RSD75`.
It explains why the final Frobenius obstruction disappears automatically in
four variables but would not disappear by the same argument in five or more.

> **Theorem HC4RSD76 — trace–focal dimension-four miracle.**
> Let `N` be a regular nilpotent block, self-adjoint for a nondegenerate
> symmetric form `S`, and put `T=SN`.  Let `B` be a symmetric form on the
> target and `J=BT`.  Assume the last target vector `ell=S ker(N)` is the
> radical direction of the second fundamental form represented by `B`.
>
> In dimension four, after the radical equations are imposed, the following
> three quantities are the same scalar (up to the harmless sign convention in
> `d lambda`):
>
> 1. `tr J`;
> 2. the value of the second fundamental form on the second Jordan/Krylov
>    direction;
> 3. the sole remaining Frobenius coefficient of `ker N^3`.
>
> Therefore nilpotence of `J` forces the complete regular-nilpotent Jordan flag
> to be Frobenius-integrable.

## 1. Canonical pair

In dimension `n`, choose the Jordan chain so that

\[
N e_1=0,\quad N e_{i+1}=e_i,
\]

and take the anti-diagonal form

\[
S_{ij}=\delta_{i+j,n+1}.
\]

Then

\[
T=SN
\]

has its nonzero entries on the shifted anti-diagonal

\[
i+j=n+2.
\]

For symmetric `B=(b_ij)`, one obtains

\[
\operatorname{tr}(BT)
=\sum_{i=2}^{n} b_{i,n+2-i}.
\tag{1.1}
\]

Symmetry of `B` groups the terms in pairs, with one central term when `n` is
even.

## 2. Dimension four

For `n=4`,

\[
\operatorname{tr}(BT)=2b_{24}+b_{33}.
\tag{2.1}
\]

The tangent hyperplane to the gradient image is spanned by target directions
`e2*,e3*,e4*`, and the Gauss line is `ell=e4*`.  Requiring `ell` to be the
radical of the second fundamental form gives

\[
b_{24}=b_{34}=b_{44}=0.
\tag{2.2}
\]

Hence

\[
\operatorname{tr}(BT)=b_{33}.
\tag{2.3}
\]

If

\[
\lambda=Se_1,
\]

then `ker N^3=ker lambda`.  Hessian symmetry gives the three coefficients of
`d lambda` on `span(e1,e2,e3)` as

\[
b_{44},\qquad b_{34},\qquad b_{33}-b_{24}.
\tag{2.4}
\]

After (2.2), exactly one survives:

\[
b_{33}.
\]

Thus

\[
\boxed{
\operatorname{tr}(BT)
=II(m,m)
=\text{Frobenius}_{\ker N^3}
=b_{33}.
}
\tag{2.5}
\]

For the HC4 normal field, `BT` is a nilpotent quasi-translation Jacobian, so
its trace vanishes.  Equation (2.5) is the final closure mechanism.

## 3. Why the same one-line argument stops in dimension five

For `n=5`, before the radical equations,

\[
\operatorname{tr}(BT)=2(b_{25}+b_{34}).
\]

The Gauss-line radical kills the entries involving the last direction, leaving

\[
\operatorname{tr}(BT)=2b_{34}.
\]

But the Frobenius equations for `ker N^4` retain, in addition to `b34`, the
independent coefficients

\[
b_{44},\qquad -b_{24}+b_{33}.
\]

Therefore trace zero removes only one of several obstructions.  Starting in
five variables one needs genuinely new identities or higher focal data.

This is one algebraic explanation for why the four-variable Hessian problem is
exceptionally rigid.

## 4. Verification

The companion checker constructs the canonical pairs for dimensions `3` through
`6` and prints the trace and restricted Frobenius coefficients before and after
the Gauss-line radical equations.  The `n=4` row has exactly one common
surviving scalar; the `n>=5` rows do not.
