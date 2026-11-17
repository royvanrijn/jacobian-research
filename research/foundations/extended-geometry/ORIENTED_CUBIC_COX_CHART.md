# Oriented cubic Cox chart with two dicritical divisors

The reciprocal three-factor suspension becomes polynomial after the target
remembers the primitive oriented discriminant.  This produces a
constant-Jacobian morphism of smooth affine threefold charts with two
distinct dicritical boundary divisors.

It is still not a polynomial map of affine three-space: the source and target
are Cox/discriminant charts.  In fact, the full normalized source is proved
below not to be affine three-space.

Work over a characteristic-zero field `k`.

## 1. Source and oriented target

For three linear factors, retain the notation

\[
 r_{ij}=u_iv_j-v_iu_j,\qquad
 P=L_1L_2L_3=p_0X^3+p_1X^2Y+p_2XY^2+p_3Y^3,
\]

and set

\[
 m=p_0+p_3.
\]

Unlike the first Cox suspension, do not invert `r_(12)`.  Define

\[
 \overline Y=
 \{r_{13}=r_{23}=m=1\}\subset\mathbb A^6.             \tag{1}
\]

The selected boundary-class minor is unimodular, so (1) is the normalized
projective chart

\[
 (\mathbb P^1)^3\setminus(R_{13}\cup R_{23}\cup E).
\]

In particular it is smooth and affine.  Its boundary-class map is square
unimodular, hence

\[
 \mathcal O(\overline Y)^*=k^*,\qquad
 \operatorname{Pic}(\overline Y)=0.                  \tag{2}
\]

These first invariants do not make it affine space.  Let

\[
 S=\operatorname{Spec}k[\xi]/(\xi^2-\xi+1)
\]

and write `L=[A^1]`.  The exact Grothendieck class is

\[
 \boxed{
 [\overline Y]=\mathbb L^3-2\mathbb L-[S].
 }                                                    \tag{3}
\]

To compute it, use `r_(13)=r_(23)=1` to write

\[
 L_2=L_1+tL_3,\qquad
 (L_1,L_3)=((a,b),(c,d))\in SL_2.
\]

Then

\[
 m=A_0+tA_1,\qquad
 A_0=a^2c+b^2d,\qquad A_1=ac^2+bd^2.
\]

Projection to `SL_2` gives one `t` when `A_1!=0`, no `t` when
`A_1=0,A_0!=1`, and an affine line of `t` values when
`A_1=0,A_0=1`.  Put

\[
 H=(A_1=0)\subset SL_2,\qquad Z=(A_1=0,A_0=1).
\]

Solving the two linear equations `A_1=0`, `ad-bc=1` gives

\[
 H\simeq\{(c,d):c^3+d^3\ne0\},
\]

so

\[
 [H]=\mathbb L^2-\mathbb L-(\mathbb L-1)[S].
\]

On `Z`, set `x=c/d`.  The equation `A_0=1` gives

\[
 d=\frac{x}{x^3+1},
\]

with `x` different from `0`, `-1`, and the degree-two locus `S`.  Hence

\[
 [Z]=\mathbb L-2-[S].
\]

Using `[SL_2]=L^3-L`,

\[
 [\overline Y]=[SL_2]-[H]+\mathbb L[Z]
 =\mathbb L^3-2\mathbb L-[S],
\]

as claimed.  After base change to `C`, the Hodge--Deligne polynomial is

\[
 E_c(\overline Y;u,v)=(uv)^3-2uv-2,
\]

not `(uv)^3`.  Therefore

\[
 \boxed{\overline Y_{\mathbb C}\not\simeq\mathbb A^3_{\mathbb C}.} \tag{4}
\]

For a binary cubic write

\[
 \Delta(P)=
 p_1^2p_2^2-4p_0p_2^3-4p_1^3p_3
 -27p_0^2p_3^2+18p_0p_1p_2p_3.                      \tag{5}
\]

The oriented target is the hypersurface

\[
 \widetilde T=
 \left\{
 p_3=1-p_0,\quad D^2=\Delta(P)
 \right\}\subset\mathbb A^4_{p_0,p_1,p_2,D}.          \tag{6}
\]

Its singular locus is the triple-root locus.  To obtain an affine smooth
chart, put

\[
 h=p_1^2-3p_0p_2
\]

and restrict to

\[
 \widetilde T_h=\widetilde T\cap D(h),\qquad
 \overline Y_h=\{y\in\overline Y:h(P(y))\ne0\}.       \tag{7}
\]

A triple cubic has `h=0`, so `T_tilde_h` is smooth.  The restriction retains
the generic point of the ordinary double-root discriminant.

## 2. Polynomial constant-Jacobian morphism

The discriminant identity

\[
 \Delta(L_1L_2L_3)=(r_{12}r_{13}r_{23})^2
\]

becomes

\[
 \Delta(P)=r_{12}^2
 \quad\hbox{on }\overline Y.                          \tag{8}
\]

Therefore

\[
 \boxed{
 \Psi:\overline Y_h\longrightarrow\widetilde T_h,
 \qquad
 (L_1,L_2,L_3)\longmapsto(P,D=r_{12})
 }                                                    \tag{9}
\]

is polynomial.  No reciprocal boundary function occurs.

### Theorem 2.1

With the complete-intersection residue form on `Y_bar_h` and the hypersurface
residue form on `T_tilde_h`, the Jacobian of (9) is the nonzero constant
`-1/2`, up to orientation.  The morphism is etale and has generic degree
three.

### Proof

The ambient determinant from the three-factor ledger is

\[
 \det D(p_0,p_1,p_2,m,r_{13},r_{23})
 =-r_{12}r_{13}^2r_{23}^2.                           \tag{10}
\]

Let `omega_(Y_bar)` be the residue of the standard six-form along
`m=r_(13)=r_(23)=1`.  Equation (10) gives

\[
 \Psi^*(dp_0\wedge dp_1\wedge dp_2)
 =-r_{12}\omega_{\overline Y}.                       \tag{11}
\]

On the dense chart `D!=0`, a hypersurface residue generator of the canonical
bundle of (6) is

\[
 \Omega_{\widetilde T}
 =\frac{dp_0\wedge dp_1\wedge dp_2}{2D}.              \tag{12}
\]

Using `D=r_(12)` in (11) yields

\[
 \boxed{
 \Psi^*\Omega_{\widetilde T}
 =-\frac12\omega_{\overline Y}.
 }                                                    \tag{13}
\]

Both sides are regular top forms on the smooth charts in (7), so equality on
the dense `D!=0` locus extends across `D=0`.  The determinant is everywhere
nonzero and (9) is etale.

For a squarefree cubic, fixing `D` selects one of the two orientations of
its six ordered factorizations.  Three normalized orderings remain, proving
generic degree three.  QED

This is stronger than the reciprocal suspension in one direction: the
collision `r_(12)=0` is now an ordinary source divisor, and the map remains
polynomial and etale there.

## 3. Two dicritical boundary divisors

At a generic point of `D=0`, the cubic has a double factor `A` and a simple
factor `B`.  Before collision, the three orderings compatible with one
orientation can be represented cyclically as

\[
 (A_1,A_2,B),\qquad
 (B,A_1,A_2),\qquad
 (A_2,B,A_1).                                        \tag{14}
\]

For the first ordering, the colliding pair occupies positions `(1,2)`.
Both normalized resultants `r_(13),r_(23)` stay nonzero, so this branch
extends to the finite divisor `r_(12)=0` inside `Y_bar`.

For the second ordering, the raw `r_(23)` tends to zero.  The unique
normalization

\[
 \lambda_1=\frac{r_{23}}m,\qquad
 \lambda_2=\frac{r_{13}}m,\qquad
 \lambda_3=\frac m{r_{13}r_{23}}                    \tag{15}
\]

has valuations `(1,0,-1)`.  One normalized factor tends to zero and another
to infinity.  This gives a boundary prime over `R_(23)`.

For the third ordering, `r_(13)` tends to zero and (15) has valuations
`(0,1,-1)`, giving a second boundary prime over `R_(13)`.

In both cases the oriented discriminant has order one:

\[
 D=\frac{r_{12}r_{13}r_{23}}{m^2}.
\]

Thus both boundary primes map dominantly and generically with ramification
index one to `D=0`.

### Theorem 3.1

In the normalized inverse graph of (9), exactly two boundary primes dominate
the generic oriented discriminant divisor:

\[
 \mathcal D_{13}\longrightarrow(D=0),\qquad
 \mathcal D_{23}\longrightarrow(D=0).                \tag{16}
\]

They are distinct dicritical divisors.  The third cyclic branch extends
inside the affine source as `r_(12)=0`.

No other boundary prime lies over the generic point of `D=0`, because (14)
exhausts the three generic sheets.  The boundary `E` maps to the omitted
projective hyperplane rather than to the affine target chart.

This gives several genuinely independent dicritical divisors on a Cox
chart, although both have the same target image.

## 4. Fiber and finite-field behavior

The morphism is not finite.  Its geometric fiber size drops from three on
`D!=0` to one at the generic point of `D=0`; the other two sheets are the
dicritical valuations in (16).

Over a finite field of good odd characteristic, oriented squarefree cubics
already prevent permutation behavior:

- a split cubic with square discriminant has three rational affine lifts;
- an irreducible cubic has Frobenius a three-cycle, hence also square
  discriminant, but has no rational ordered-factor lift.

The `F_5` regression gives exactly thirty oriented target points with three
affine lifts and twenty-five discriminant points with one affine lift.
Thus the oriented target removes the weighted `C=0` excess but not the
factorization-cover obstruction.

## 5. What remains

The construction settles the first target-absorption question:

1. the primitive boundary character is placed on the target as `D`;
2. the reciprocal suspension disappears;
3. the map extends etale across one collision divisor; and
4. two other collision branches become distinct dicritical divisors.

The full source is proved not to be affine three-space by (3)--(4).  The
principal source chart in (7) is not separately classified.  The
construction also does not separate the two dicritical divisors by target
image, and its symmetric degree-three inverse cover rules out the desired
finite-field permutation behavior.

The next affine-space problem is therefore narrower: straighten
`Y_bar_h` and `T_tilde_h` to polynomial affine charts, or find an affine
modification which preserves (13) while separating the two valuations.

For the tangent hyperplane `p_1=1`, the natural involution quotient does
produce affine three-space.  The
[affine descent theorem](ORIENTED_CUBIC_AFFINE_DESCENT.md) identifies the
quotient map with the foundational cubic Keller map and proves that the two
dicritical divisors are exchanged and merge downstairs.
The companion
[quotient-rigidity theorem](ORIENTED_CUBIC_QUOTIENT_RIGIDITY.md) proves
that the tangent oriented source is not even stably affine and that no
other factor-permutation quotient can preserve both dicritical valuations
individually.
More generally, the
[linear-hyperplane classification](LINEAR_HYPERPLANE_COX_CLASSIFICATION.md)
shows that no nonzero linear functional of the four cubic coefficients
gives an affine or stably affine normalized source.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_oriented_cubic_cox_chart.py
```

The checker verifies the discriminant cover, the ambient and residue
Jacobians, the unimodular source boundary matrix, the two normalization
valuation vectors, the motivic point-count realization, and the exact `F_5`
fiber profile.
