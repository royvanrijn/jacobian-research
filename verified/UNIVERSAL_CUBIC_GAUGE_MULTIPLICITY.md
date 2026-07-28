# Universal cubic multiplicity by fiber-invisible gauge lifts

The weighted, cancellation, and minimal diagonal quadratic-gauge
presentations of a cubic finite etale algebra all lie in the foundational
stable class.  This note escapes that collapse by changing the lift of the
cubic coefficient away from the selected target while keeping its value at
that target fixed.

Together with the power-shifted quartic theorem and the translated
higher-degree theorem, the result makes rank three the exact universal
multiplicity threshold.

## The theorem

> **Universal cubic gauge-multiplicity theorem.**
> Let `K` be any characteristic-zero field and let `A` be any rank-three
> finite etale `K`-algebra.  Then
> \[
>   \boxed{|\mathcal R_K(A)|=\infty.}
> \]
> All classes may be represented by determinant-one polynomial maps
> `A^3_K -> A^3_K` of geometric degree three.

The maps are indexed by integers `n>=4`.  Over an algebraic closure, the
canonical finite-normalization boundary of the `n`-th map has exactly
`n-1` unramified vertical target components in addition to its irreducible
ramified discriminant component.  This component count is stable under
polynomial left--right equivalence and identity stabilization.

## 1. The fixed cubic algebra

Choose a monic squarefree cubic `f(T)` with

\[
 A\simeq K[T]/(f).
\]

Choose `a in K` with `f'(a)!=0`, and write

\[
 G(S)=f(a+S)-f(a)
     =g_1S+g_2S^2+g_3S^3.                             \tag{1.1}
\]

Then `g_1g_3!=0`.  For every integer `n>=4`, put

\[
\boxed{
 G_{P,n}(S)
 =g_1S+g_2PS^2+
   g_3P\bigl(1+P^{n-1}-P^2\bigr)S^3.
}                                                       \tag{1.2}
\]

At `P=1`, the factor in parentheses equals one, so

\[
 G_{1,n}(S)=G(S)                                       \tag{1.3}
\]

for every `n`.  The deformation is invisible to the selected inverse
quotient.

## 2. The polynomial maps

Put

\[
 t=1+xy,\qquad
 q=t^2z+\frac{g_1}{g_3}y^2(1+3t),\qquad P=tq.          \tag{2.1}
\]

Start with the minimal cubic quadratic gauge and add

\[
\begin{aligned}
\Delta B_n
 &=3\frac{g_3}{g_1}
   \left(t^{n-1}xq^n-t^2xq^3\right),\\
\Delta C_n
 &=-\frac{g_3}{g_1}
   \left(t^{n-3}x^3q^n-x^3q^3\right).                 \tag{2.2}
\end{aligned}
\]

Thus

\[
\begin{aligned}
B_n={}&y+3\frac{g_3}{g_1}xq
          +2\frac{g_2}{g_1}tq+\Delta B_n,\\
C_n={}&x(5-3t)-\frac{g_3}{g_1}x^3z+\Delta C_n,
\end{aligned}                                         \tag{2.3}
\]

and `F_n=(P,B_n,C_n)`.

On `t!=0`, use

\[
 S=\frac{x}{t},\qquad Q=y+xq,\qquad
 D=1-S(Q-PS)=\frac1t.                                  \tag{2.4}
\]

Substitution gives

\[
\boxed{
 E_n(P,B,C;S)
 =G_{P,n}(S)-\frac{g_1}{2}(BS^2+C)=0
}                                                       \tag{2.5}
\]

and, on the source incidence,

\[
\boxed{\partial_SE_n=g_1D.}                            \tag{2.6}
\]

The paired corrections in (2.2) are exactly the slope and intercept
corrections contributed by
`g_3(P^n-P^3)S^3`.  At fixed `P`, the marked-line plane Jacobian is `-2D`;
the reciprocal chart contributes `D^{-1}`.  Consequently

\[
\boxed{\det DF_n=-2.}                                  \tag{2.7}
\]

A fixed target scaling makes the determinant one.

The generic inverse equation (2.5) is an irreducible cubic.  Indeed it is
primitive and linear in `C`, with nonzero constant coefficient of `C`.
Every generic simple root reconstructs through (2.4), so

\[
 \operatorname{gdeg}(F_n)=3.                           \tag{2.8}
\]

## 3. The common complete fiber

At

\[
 y_n=\left(1,0,-\frac{2f(a)}{g_1}\right),              \tag{3.1}
\]

equations (1.3) and (2.5) give

\[
 E_n(y_n;S)=G(S)+f(a)=f(a+S).                          \tag{3.2}
\]

This polynomial is squarefree.  Its derivative at every root is nonzero,
so (2.6) puts every root in the reconstruction chart.  The fiber contains
all three inverse sheets and therefore

\[
\boxed{
 F_n^{-1}(y_n)
 \simeq\operatorname{Spec}K[S]/(f(a+S))
 \simeq\operatorname{Spec}A.
}                                                       \tag{3.3}
\]

Thus every map has the same complete cubic fiber.

## 4. The degree-drop polynomial is squarefree

Let

\[
 h_n(P)=1+P^{n-1}-P^2.                                 \tag{4.1}
\]

It has degree `n-1`, nonzero constant term, and `h_n(1)=1`.  It is
squarefree in characteristic zero.  Indeed, a common root of `h_n` and
`h_n'` would be nonzero and satisfy

\[
 P^{n-3}=\frac2{n-1},\qquad
 P^2=\frac{n-1}{n-3}.                                  \tag{4.2}
\]

The second equation makes `P` real with absolute value greater than one,
while the first makes its absolute value less than one, a contradiction
after embedding the prime field and these algebraic equations into
`\mathbb C`.

Hence `h_n` has exactly `n-1` distinct nonzero geometric roots

\[
 \rho_1,\ldots,\rho_{n-1}.                             \tag{4.3}
\]

None equals one.

## 5. One boundary component over every root of `h_n`

Fix a root `rho` of `h_n` and use `u=P-rho` as a parameter.  At the generic
point of the target divisor `P=rho`, the cubic coefficient in (2.5) has
valuation one, while the quadratic coefficient is a unit.  The lower
Newton polygon is

\[
 (0,0)\longrightarrow(2,0)\longrightarrow(3,1).        \tag{5.1}
\]

The horizontal block gives two finite simple roots.  Since `rho!=0`, they
reconstruct to affine source points.  The final length-one block gives one
root with

\[
 v(S)=-1,\qquad v(D)=-1,\qquad v(q)=v(PD)=-1.          \tag{5.2}
\]

Thus this third branch lies outside affine source space.  Its Newton block
has horizontal length one, so its ramification and residue degrees are

\[
 (e,f)=(1,1).                                          \tag{5.3}
\]

Every root `rho_i` therefore supplies one distinct unramified boundary
prime whose target image is the hyperplane `P=rho_i`.

## 6. No boundary component is hidden over `P=0`

At `P=0`, the inverse equation becomes the generic quadratic

\[
 g_1S-\frac{g_1}{2}(BS^2+C)=0.                         \tag{6.1}
\]

Its two roots are the two affine branches on `q=0`.  The third inverse
branch is also affine: it is the divisor `t=0`.  On that divisor,

\[
 xy=-1,\qquad q=\frac{g_1}{g_3}y^2,\qquad B_n=-2y,
                                                               \tag{6.2}
\]

and

\[
 C_n
 =5x-\frac{g_3}{g_1}x^3z
   +\frac{g_3}{g_1}x^3q^3.                             \tag{6.3}
\]

Since `x` and `y` are nonzero and (6.3) is linear in `z`, this is one
regular residue-degree-one affine branch.  The local degree sum is

\[
 2+1=3,                                                \tag{6.4}
\]

so no boundary prime lies over `P=0`.

## 7. Boundary exhaustion and stable separation

If `P` is nonzero, `h_n(P)` is nonzero, and the inverse root is simple,
(2.4)--(2.6) reconstruct it regularly.  Sections 5 and 6 exhaust every
degree-drop divisor.  The only remaining boundary image is the repeated-root
discriminant.  Its finite-root parametrization by `(P,r)` is irreducible and
its generic boundary prime has ramification index two.

Therefore, over an algebraic closure, the complete list of irreducible
target images of canonical boundary divisors is:

\[
\begin{array}{c|c|c}
\text{target image}&\text{number}&\text{generic boundary label}\\ \hline
\text{repeated-root discriminant}&1&(2,1)\\
P=\rho_i,\ h_n(\rho_i)=0&n-1&(1,1).
\end{array}                                             \tag{7.1}
\]

In particular, `F_n` has exactly `n` geometric boundary target components.

The canonical finite-normalization package and its boundary primes are
functorial under polynomial left--right equivalence.  Identity stabilization
takes every component times affine space and changes neither irreducibility
nor the number of components.  Hence stable equivalence of `F_n` and `F_m`
would force

\[
 n=m.                                                   \tag{7.2}
\]

The maps `F_4,F_5,F_6,\ldots` are pairwise stably inequivalent.  Combining
(7.2) with the common complete fiber (3.3) proves the theorem.

## 8. Consequences and scope

The
[low-rank collapse theorem](LOW_RANK_MULTIPLICITY_BOUNDARIES.md) remains
correct for its three stated minimal mechanisms.  The maps here lie outside
that scope: they deform the cubic gauge lift while fixing its value at the
selected fiber.

Together with
[power-shifted quartic multiplicity](UNIVERSAL_QUARTIC_GAUGE_MULTIPLICITY.md)
and
[higher-rank multiplicity](UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md), this
gives

\[
\boxed{
 |\mathcal R_K(A)|=\infty
 \quad\text{for every characteristic-zero }K
 \text{ and every finite etale }A/K\text{ of rank }N\ge3.
}                                                       \tag{8.1}
\]

Rank two cannot occur as the geometric degree of a noninvertible Keller map,
while rank one is represented by polynomial automorphisms.  Thus rank three
is the exact threshold for universal infinite stable multiplicity among
noninvertible Keller maps.

## 9. Exact regression

Run

```bash
.venv/bin/python scripts/verify_universal_cubic_gauge_multiplicity.py
```

The checker verifies the determinant, inverse and derivative identities,
the fixed connected selected cubic `T^3-T-1`, squarefreeness of `h_n` in
the regression range, the Newton blocks, the affine `P=0` decomposition,
and the boundary counts.  It also calls the public
`compile_polynomial_to_keller_fiber(..., stable_parameter=k)` interface and
checks that the compiled map retains the inverse cubic while returning
exponent `n=k+4` and stable separating value `n`.  Boundary functoriality is
a written theorem.
