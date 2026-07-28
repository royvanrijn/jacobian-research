# Universal power-shifted gauge multiplicity in every rank at least four

The rank-four power shift changes the exponent of its quartic decoration
without changing the selected inverse polynomial.  The same idea works
uniformly in every degree once all coefficients of degree at least four
receive one common extra power of the first target coordinate.

This gives one unconditional construction for all ranks `N>=4`.  It avoids
the quartic trace quadric and the separate translated-coefficient invariants
in higher degrees.

## The theorem

> **Universal power-shifted gauge theorem.**
> Let `K` be a characteristic-zero field, let `N>=4`, and let `A` be a
> rank-`N` finite etale `K`-algebra.  Then `A` occurs as a common complete
> fiber of infinitely many pairwise stably inequivalent determinant-one
> polynomial maps `A^3_K -> A^3_K` of geometric degree `N`.
>
> The maps are indexed by `m>=0`.  On the normalized intrinsic ramified
> stratum, the Newton polygon of the relative Fitting divisor has normalized
> area
> \[
>   \boxed{2N-3+(N-2)m.}
> \]

Together with the fiber-invisible cubic theorem, this proves universal
infinite stable multiplicity in every possible noninvertible degree.

## 1. A fully nonzero translated presentation

Choose a monic squarefree polynomial `f(T)` of degree `N` with

\[
 A\simeq K[T]/(f).
\]

Every derivative `f^{(j)}` for `1<=j<=N` is a nonzero polynomial.  Since
`K` is infinite, choose `a in K` outside the finite union of the zero sets
of

\[
 f',f''',f^{(4)},\ldots,f^{(N)}.                       \tag{1.1}
\]

Write

\[
 G(S)=f(a+S)-f(a)=\sum_{j=1}^Ng_jS^j.                 \tag{1.2}
\]

Then

\[
 g_1g_3g_4\cdots g_N\ne0.                              \tag{1.3}
\]

The coefficient `g_2` may vanish; it is a removable target shear and is not
used by the stable invariant.

## 2. The common power shift

For `m>=0`, define

\[
\boxed{
 G_{P,m}(S)
 =g_1S+P(g_2S^2+g_3S^3)
   +\sum_{j=4}^Ng_jP^{j+m}S^j.
}                                                       \tag{2.1}
\]

At `P=1`,

\[
 G_{1,m}(S)=G(S)                                       \tag{2.2}
\]

for every `m`.  Thus the shift is invisible on the selected fiber.

Put

\[
 t=1+xy,\qquad
 q=t^2z+\frac{g_1}{g_3}y^2(1+3t),\qquad P=tq.          \tag{2.3}
\]

The associated polynomial map `F_m=(P,B_m,C_m)` is

\[
\begin{aligned}
B_m={}&y+3\frac{g_3}{g_1}xq
       +2\frac{g_2}{g_1}tq
       +\sum_{j=4}^N
          j\frac{g_j}{g_1}
          t^{m+2}x^{j-2}q^{j+m},\\
C_m={}&x(5-3t)-\frac{g_3}{g_1}x^3z
       -\sum_{j=4}^N
          (j-2)\frac{g_j}{g_1}
          t^mx^jq^{j+m}.
                                                               \tag{2.4}
\end{aligned}
\]

The inequalities `j+m>=j` are exactly the termwise polynomiality gate.

## 3. Keller identity and geometric degree

On `t!=0`, put

\[
 S=\frac{x}{t},\qquad Q=y+xq,\qquad
 D=1-S(Q-PS)=\frac1t.                                  \tag{3.1}
\]

The inverse equation is

\[
\boxed{
 E_m(P,B,C;S)
 =G_{P,m}(S)-\frac{g_1}{2}(BS^2+C)=0.
}                                                       \tag{3.2}
\]

On the source incidence,

\[
\boxed{\partial_SE_m=g_1D.}                            \tag{3.3}
\]

Each shifted term contributes the marked-line slope/intercept pair

\[
 j\frac{g_j}{g_1}P^{j+m}S^{j-2},
\qquad
 -(j-2)\frac{g_j}{g_1}P^{j+m}S^j.
\]

Consequently the plane Jacobian is `-2D`; the reciprocal chart contributes
`D^{-1}`.  Hence

\[
\boxed{\det DF_m=-2.}                                  \tag{3.4}
\]

A fixed target scaling makes the determinant one.

Equation (3.2) is primitive and linear in `C`, whose coefficient is a
nonzero constant.  It is therefore irreducible over `K(P,B,C)`.  Its degree
in `S` is `N`, and its generic roots are simple and reconstruct by (3.1).
Thus

\[
 \operatorname{gdeg}(F_m)=N.                           \tag{3.5}
\]

## 4. The selected complete fiber

At

\[
 y_m=\left(1,0,-\frac{2f(a)}{g_1}\right),              \tag{4.1}
\]

equations (2.2) and (3.2) give

\[
 E_m(y_m;S)=G(S)+f(a)=f(a+S).                          \tag{4.2}
\]

This polynomial is squarefree.  Equation (3.3) puts every root in the
regular reconstruction chart.  All `N` generic sheets are present, and

\[
\boxed{
 F_m^{-1}(y_m)
 \simeq\operatorname{Spec}K[S]/(f(a+S))
 \simeq\operatorname{Spec}A.
}                                                       \tag{4.3}
\]

## 5. Intrinsic selection of the ramified stratum

Away from `P=0` and the repeated-root discriminant, every inverse root
reconstructs.  At `P=0`, the lower Newton polygon has the three blocks

\[
 (0,0)\longrightarrow(2,0)\longrightarrow(3,1)
 \longrightarrow(N,N+m).                              \tag{5.1}
\]

The first two blocks are the affine `q=0` residue-degree-two branch and the
affine `t=0` residue-degree-one branch.  The last block has horizontal
length `N-3`.  If

\[
 h_m=\gcd(N-3,m+2),
\]

it gives `h_m` boundary primes, each with

\[
 (e,f)=\left(\frac{N-3}{h_m},1\right).                 \tag{5.2}
\]

Over the repeated-root discriminant there is one boundary prime with
`(e,f)=(2,1)` and `N-2` affine residue-degree-one sheets.  This is
intrinsically distinct from the `P=0` ledger, whose affine part is
`(1,2)+(1,1)`.  Hence the two target boundary images are intrinsically
ordered, and the ramified discriminant stratum is selected without using
the displayed coordinate `P`.

After deleting `P=0`, its normalization is

\[
 \widetilde Z_\Delta^\circ
 \simeq\operatorname{Spec}K[P^{\pm1},r^{\pm1}]
 \simeq\mathbb G_m^2.                                  \tag{5.3}
\]

## 6. The Fitting Newton polygon

Put `a_j=g_j/g_1`.  On (5.3), relative differentiation of the tangent-line
parametrization gives

\[
\operatorname{Fitt}_0
\Omega_{\widetilde Z_\Delta^\circ/Z_\Delta^\circ}
=(J_m),
\]

where

\[
\boxed{
 J_m
 =-1+3a_3Pr^2+
   \sum_{j=4}^N j(j-2)a_jP^{j+m}r^{j-1}.
}                                                       \tag{6.1}
\]

By (1.3), its Laurent support is exactly

\[
\mathcal S_{N,m}
=\{(0,0),(1,2)\}
 \cup\{(j+m,j-1):4\le j\le N\}.                       \tag{6.2}
\]

The points in the second set are collinear.  The convex hull has vertices

\[
 (0,0),\quad(4+m,3),\quad(N+m,N-1),\quad(1,2),         \tag{6.3}
\]

with the middle two coinciding when `N=4`.  The shoelace formula gives twice
its Euclidean area, equivalently its normalized lattice area:

\[
\boxed{
 \operatorname{Area}_{\mathbb Z}
 \operatorname{Newt}(J_m)
 =2N-3+(N-2)m.
}                                                       \tag{6.4}
\]

Stable normalization functoriality transports the selected stratum and its
relative Fitting divisor.  Stabilization adds affine variables but no units.
An induced isomorphism acts on the unit lattice by `GL_2(\mathbb Z)`, and
multiplication of `J_m` by a unit translates its support.  Both operations
preserve normalized Newton area.

For fixed `N>=4`, (6.4) is strictly increasing in `m`.  Therefore

\[
 F_m\sim_{\mathrm{stable}}F_{m'}
 \quad\Longrightarrow\quad m=m'.                       \tag{6.5}
\]

Equations (4.3) and (6.5) prove the theorem.

## 7. Relation to the rank-specific proofs

For `N=4`, (6.4) is `5+2m`, the lattice index of the triangular quartic
Fitting support.  Thus the power-shifted quartic theorem is the first case
of the present result.

For `N>=5`, the construction supplies an alternative to translating a
primitive generator through the minimal diagonal coefficient quotient.
The fiber polynomial and marked target can remain fixed while the gauge
lift moves.  No trace moment or nonconstant rational coefficient invariant
is required.

Rank three retains a special low cubic skeleton and is handled by the
[fiber-invisible cubic theorem](UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md).

## 8. Exact regression

Run

```bash
.venv/bin/python scripts/verify_universal_power_shifted_gauge_multiplicity.py
```

The checker verifies the polynomial map, inverse, derivative, and determinant
identities in representative degrees and shifts; the fixed-fiber identity;
the `P=0` Newton ledger; the complete Fitting support; and the normalized
area formula through a broad exact range.  It also exercises the public
`compile_polynomial_to_keller_fiber(..., stable_parameter=m)` interface,
including equality with the minimal quartic gauge at `m=0`, exact support
and area certificates, automatic avoidance of vanishing higher translated
coefficients, and rejection of an inadmissible explicit translation.
