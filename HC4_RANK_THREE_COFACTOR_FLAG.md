# Canonical cofactor flag for the final rank-three `[4]` HC4 stratum

## Status and scope

Let

\[
S=\operatorname{Hess}\psi,
\qquad
T=\operatorname{Hess}A,
\qquad
\det S=\delta\in K^*,
\qquad
N=S^{-1}T.
\]

Assume throughout that `N` is nilpotent of Jordan type `[4]`. Equivalently,

\[
N^4=0,
\qquad
\operatorname{rank}N=3.
\]

This is the only moving nilpotent Jordan stratum left after the rank-one and
rank-two closures.

> **Theorem HC4RSD64 — canonical polynomial cofactor flag.**
> Expand
>
> \[
> \operatorname{adj}(S+sT)=C_0+sC_1+s^2C_2+s^3C_3.
> \]
>
> Then
>
> \[
> C_j=\delta(-1)^jN^jS^{-1},\qquad 0\le j\le3,
> \tag{0.1}
> \]
>
> and therefore
>
> \[
> \operatorname{rank}C_j=4-j.
> \tag{0.2}
> \]
>
> Moreover every `C_j` is symmetric and divergence-free, and
>
> \[
> TC_j=-SC_{j+1}\quad(0\le j<3),
> \qquad TC_3=0.
> \tag{0.3}
> \]
>
> In particular
>
> \[
> C_3=\operatorname{adj}T
> \]
>
> is a canonical polynomial rank-one matrix whose image is exactly
> \(\ker T=\ker N\).  Thus the final moving kernel line is encoded without
> choosing a rational Jordan frame.

> **Corollary HC4RSD64a — projectively straight top kernel.**
> On a dense open set write
>
> \[
> C_3=\rho\,k k^{\mathsf T}.
> \]
>
> Then
>
> \[
> Tk=0
> \]
>
> and
>
> \[
> (k\cdot\nabla)k\in K(x)\,k.
> \tag{0.4}
> \]
>
> Hence the top kernel line is invariant along its own characteristic flow.
> After rational rescaling of `k`, the characteristic is a quasi-translation
> direction.

The point is not merely that `ker T` is one-dimensional.  Equations (0.1)--
(0.3) provide the full nested polynomial flag

\[
\operatorname{im}C_3
\subset
\operatorname{im}C_2
\subset
\operatorname{im}C_1
\subset K(x)^4,
\tag{0.5}
\]

of dimensions `1,2,3`.  This is the intrinsic replacement for a moving Jordan
basis.

## 1. Adjugate expansion

Since

\[
S+sT=S(I+sN),
\]

nilpotence gives the finite inverse

\[
(S+sT)^{-1}
=(I+sN)^{-1}S^{-1}
=(I-sN+s^2N^2-s^3N^3)S^{-1}.
\]

The determinant is identically `delta`, so

\[
\operatorname{adj}(S+sT)
=\delta(S+sT)^{-1},
\]

which is exactly (0.1).  A single Jordan block has

\[
\operatorname{rank}N^j=4-j,
\]

proving (0.2).

Because `S+sT` is symmetric, its adjugate and each coefficient `C_j` are
symmetric.

## 2. Piola identities

For any Hessian matrix `Hess f`, every row of its cofactor matrix is
divergence-free:

\[
\sum_j\partial_j\operatorname{adj}(\operatorname{Hess}f)_{ij}=0.
\tag{2.1}
\]

Apply this to

\[
f=\psi+sA.
\]

The identity is polynomial in `s`, so every coefficient separately satisfies

\[
\operatorname{div}C_j=0.
\tag{2.2}
\]

Multiplying (0.1) by `T=SN` gives

\[
TC_j
=\delta(-1)^jSN^{j+1}S^{-1}
=-SC_{j+1},
\]

and `TC_3=0` because `N^4=0`.

## 3. Straightness of the top characteristic

Since `T` has rank three,

\[
C_3=\operatorname{adj}T
\]

has rank one and its image is `ker T`.  Write

\[
C_3=\rho kk^{\mathsf T}.
\]

Taking row divergence gives

\[
0=\operatorname{div}(\rho kk^{\mathsf T})
 =\rho(k\cdot\nabla)k
  +\bigl(\operatorname{div}(\rho k)\bigr)k.
\]

Thus

\[
(k\cdot\nabla)k
=-\frac{\operatorname{div}(\rho k)}{\rho}\,k,
\]

which proves (0.4).

The same conclusion also follows directly from Hessian symmetry: differentiating
`Tk=0` along `k` and using total symmetry of the third derivative tensor gives

\[
T((k\cdot\nabla)k)=0.
\]

As `ker T` is one-dimensional, the derivative is proportional to `k`.

## 4. What remains

The final `[4]` problem is therefore not an arbitrary polynomial matrix
problem.  It consists of a divergence-free nested cofactor flag whose first
line is a quasi-translation characteristic.  Two particularly favorable
subcases now collapse immediately:

1. the top kernel line is constant; this is closed by `HC4RSD65`;
2. `A` is homogeneous; Gordan--Noether makes the top kernel line constant,
   so this is also closed by `HC4RSD65`.

Thus any surviving `[4]` obstruction must be simultaneously

- nonhomogeneous;
- rank three;
- genuinely moving at the top kernel line;
- compatible with the lower cofactor flag `C_2,C_1` and all Piola identities.

This is the narrow frontier for the next attack.
