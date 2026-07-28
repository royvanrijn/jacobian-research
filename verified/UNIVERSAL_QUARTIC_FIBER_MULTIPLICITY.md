# Universal quartic fiber multiplicity

> **Later strengthening.**  The
> [power-shifted gauge theorem](UNIVERSAL_QUARTIC_GAUGE_MULTIPLICITY.md)
> proves the same infinitude statement over every characteristic-zero field.
> The weighted theorem below remains a smaller-degree construction and
> records the exact trace-chord boundary of that mechanism.

This note proves the rank-four case of the universal multiplicity question
for Keller fibers over number fields.  It is stronger than producing three
examples: every quartic finite etale algebra over a number field occurs in
infinitely many stable polynomial left--right classes of Keller maps.

For a characteristic-zero field `K` and a finite etale `K`-algebra `A`, let
`\mathcal R_K(A)` be the set of stable polynomial left--right classes of
Keller maps having a complete fiber isomorphic to `Spec A`, as in
[the common-fiber note](COMMON_ARITHMETIC_FIBERS.md#1-the-invariant).

## The theorem

> **Universal quartic multiplicity theorem.**  
> For every number field `K` and every rank-four finite etale `K`-algebra
> `A`,
> \[
>  \boxed{|\mathcal R_K(A)|=\infty.}
> \]
> The infinitely many classes may all be represented by determinant-one,
> geometric-degree-four, boundary-clean weighted Keller maps
> `A^3_K -> A^3_K`.

Thus the universal lower bound `|\mathcal R_K(A)|>=3` holds in rank four
over every number field, with infinity in place of three.  The theorem
concerns the abstract algebra `A`: the primitive generator and the
surrounding weighted map are allowed to vary.

## 1. The quartic tangent-chord equation

Let `eta in A` be primitive with `Tr(eta)=0`, and write its characteristic
polynomial as

\[
 f_\eta(T)=T^4+pT^2+qT+r.
\]

Newton's identity gives

\[
 p=-\frac12\operatorname{Tr}(\eta^2).                 \tag{1.1}
\]

Choose `u,d in K`, with `d!=0`, and use the affine generator

\[
 T=u+dW.
\]

The weighted presentation condition says that the tangent at `u` meets the
quartic again at `u+d`:

\[
 f_\eta(u+d)-f_\eta(u)=d f_\eta'(u).                  \tag{1.2}
\]

Direct expansion gives

\[
 d^2\left(d^2+4ud+6u^2
                 -\frac12\operatorname{Tr}(\eta^2)\right)=0.
                                                               \tag{1.3}
\]

Put

\[
 e=d+2u.
\]

For `d!=0`, equation (1.3) becomes the trace quadric

\[
 \boxed{
 \operatorname{Tr}(\eta^2)=2e^2+4u^2.
 }                                                        \tag{1.4}
\]

It lives in the five-dimensional `K`-vector space

\[
 A_0\oplus Ke\oplus Ku,
 \qquad A_0=\ker(\operatorname{Tr}_{A/K}).
\]

## 2. Hasse--Minkowski isotropy and the good rational open

The trace pairing on a finite etale algebra is nondegenerate.  Since
`\operatorname{Tr}(1)=4`, its restriction to `A_0=1^\perp` is also
nondegenerate.

At any real place `v` of `K`, write

\[
 A\otimes_{K,v}\mathbb R
 \simeq\mathbb R^{r_1}\times\mathbb C^{r_2},
 \qquad r_1+2r_2=4.
\]

The quadratic form `eta -> Tr(eta^2)` on `A_0` has signature

\[
 (r_1+r_2-1,r_2).
\]

Consequently the five-variable form

\[
 \mathcal Q_A(\eta,e,u)
 =\operatorname{Tr}(\eta^2)-2e^2-4u^2                \tag{2.1}
\]

has signature

\[
 (3,2),\qquad(2,3),\qquad(1,4)
\]

for signatures `(r_1,r_2)=(4,0),(2,1),(0,2)`,
respectively.  Thus it is isotropic at every real completion; it is
automatically isotropic at every complex completion.  At every
nonarchimedean completion `K_v`, every quadratic form in five variables is
isotropic because `u(K_v)=4`.  The Hasse--Minkowski theorem therefore makes
`\mathcal Q_A` isotropic over `K`.

The projective quadric

\[
 Q_A=V(\mathcal Q_A)\subset\mathbb P(A_0\oplus K^2)
                                                               \tag{2.2}
\]

is smooth of dimension three and has a `K`-point.  Hence it is `K`-rational,
and `Q_A(K)` is Zariski dense.

Remove the following proper closed loci:

1. `e=0` or `d=e-2u=0`;
2. nonprimitive `eta`;
3. the finite set of parameter values excluded from the exact-double,
   boundary-clean weighted locus.

The nonprimitive elements form a finite union of proper `K`-linear
subspaces: finite etale subalgebras correspond after Galois closure to
quotients of a finite Galois set.  None contains `Q_A`.  Density therefore
leaves a nonempty `K`-rational good open.

## 3. Reconstruction of the weighted seed and fiber

For a good `K`-point `(eta,e,u)`, put

\[
 d=e-2u
\]

and define

\[
 P(W)=
 -\frac{f_\eta(u+dW)}{2d^3e}.                           \tag{3.1}
\]

Subtract its constant and linear terms:

\[
 H(W)=P(W)-P(0)-P'(0)W.
\]

Equations (1.2)--(1.4) give the exact normal form

\[
 \boxed{
 H(W)=W^2(W-1)(\alpha W-\alpha-1),
 \qquad
 \alpha=-\frac d{2e}=\frac ue-\frac12.
 }                                                        \tag{3.2}
\]

In particular,

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1.
\]

Away from the finite values

\[
 \alpha=0,1,-1
\]

the seed has exact degree four, is weighted-admissible, and has an exact
double zero with otherwise simple primitive roots.  Its Hessian is

\[
 H''(W)=2\bigl(6\alpha W^2-6\alpha W+\alpha-3W+1\bigr),
\]

with discriminant

\[
 \operatorname{Disc}_W(H'')
 =12(4\alpha^2+4\alpha+3).                              \tag{3.3}
\]

The final quadratic has discriminant `-32` as a polynomial in `alpha`, so
it vanishes for at most two parameter values over `K`.  Removing those and
any additional collision with a marked point costs only finitely many
further parameter values.

Let `F_H` be the determinant-one weighted map attached to `H`.  Its inverse
polynomial at the target

\[
 y_{\eta,u,e}
 =\bigl(P(0),-P'(0),1\bigr)                             \tag{3.4}
\]

is exactly `P(W)`.  Since scaling a polynomial and making the affine
generator change `T=u+dW` preserve its quotient algebra,

\[
 K[W]/(P)
 \simeq K[T]/(f_\eta)
 \simeq A.                                              \tag{3.5}
\]

Primitivity makes `f_eta` squarefree of degree four.  The weighted
reconstruction theorem on `C=1` therefore proves that (3.5) is the complete
fiber, not merely a subcollection of inverse points.

## 4. Infinitely many stable classes

On `e!=0`, the normalized seed parameter is the rational function

\[
 \alpha=\frac ue-\frac12:Q_A\dashrightarrow\mathbb A^1.
                                                               \tag{4.1}
\]

It is nonconstant: a constant value would put the irreducible quadric
`Q_A` inside one hyperplane `u=lambda e`, which is impossible for the
nondegenerate form (2.1).  Since `Q_A(K)` is Zariski dense, (4.1) takes
infinitely many `K`-values on the good open.

The
[weighted selected-root Torelli theorem](../extended-geometry/SELECTED_ROOT_TORELLI_AUDIT.md#1-the-precise-pencil-torelli-theorem)
and the
[decorated-normalization stable theorem](../extended-geometry/DECORATED_NORMALIZATION_INVARIANT.md#stable-functoriality-including-multiplicities)
apply on this exact-double, Hessian-clean, boundary-clean locus.  They
recover the normalized seed, including its distinguished affine root sheet,
from the stabilized polynomial map.  Hence two different values of
`alpha` give different stable polynomial left--right classes.

Combining this with Section 3 proves

\[
 |\mathcal R_K(A)|=\infty.
\]

## 5. Scope

The number-field hypothesis is used only to force a rational point on the
five-variable trace-chord quadric by the local--global theorem.  Over an
arbitrary characteristic-zero field, the same construction proves the
result whenever that quadric is isotropic.  The hypothesis is not automatic:
the [low-rank boundary note](LOW_RANK_MULTIPLICITY_BOUNDARIES.md) gives a
connected biquadratic field over `Q((a))((b))` whose trace-chord form is
anisotropic by Springer's theorem.  Ranks at least five need no
quadratic-form input and are treated uniformly in the
[all-rank multiplicity note](UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md).

## 6. Exact regression

Run

```bash
.venv/bin/python scripts/verify_universal_quartic_fiber_multiplicity.py
```

The checker verifies the quartic tangent-chord factorization, the trace
quadric, its three possible real signatures, the normalized weighted seed,
the formula `alpha=u/e-1/2`, and the finite clean-locus exclusions.
Hasse--Minkowski, the local `u`-invariant statement, and rational density
are the written arithmetic proof, not a bounded computation.
