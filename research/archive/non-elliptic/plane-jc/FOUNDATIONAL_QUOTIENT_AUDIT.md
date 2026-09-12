# Plane-quotient audit for the foundational Keller map

This note tests the most direct invariant two-dimensional quotients of the
foundational threefold map.  It does **not** prove the two-dimensional
Jacobian conjecture and it does not classify nonlinear algebraic foliations.

Put

\[
 u=1+xy
\]

and write the foundational map as

\[
 F=(a,b,c)=\left(
 u^3z+y^2u(4+3xy),
 y+3xu^2z+3xy^2(4+3xy),
 2x-3x^2y-x^3z
 \right).
 \tag{1}
\]

Its Jacobian determinant is `-2`.

## 1. The torus quotients are already affine planes

Use the source weights `(1,-1,-2)` and target weights `(-2,-1,1)`.  For the
source action, a weight-zero monomial satisfies

\[
 i-j-2k=0,
\]

so it is `(xy)^j(x^2z)^k`.  Hence

\[
 k[x,y,z]^{\mathbb G_m}=k[U,V],\qquad U=xy,\quad V=x^2z.
 \tag{2}
\]

Likewise a target weight-zero monomial `a^i b^j c^k` satisfies
`k=2i+j`, and therefore

\[
 k[a,b,c]^{\mathbb G_m}=k[P,Q],\qquad P=ac^2,\quad Q=bc.
 \tag{3}
\]

Both categorical quotients are genuine smooth affine planes, not singular
toric surfaces.

More generally, write an equivariant map of these weights on `x != 0` as

\[
 F=\bigl(x^{-2}A(U,V),x^{-1}B(U,V),xC(U,V)\bigr),
 \tag{4}
\]

with the displayed expressions polynomial in `(x,y,z)`.  Its descended map is

\[
 \overline F(U,V)=\bigl(A(U,V)C(U,V)^2,\,B(U,V)C(U,V)\bigr).
 \tag{5}
\]

The invariant determinant identity gives the exact formula

\[
 \boxed{\det D\overline F=C^2\det DF.}
 \tag{6}
\]

Indeed, direct differentiation of `(AC^2,BC)` gives `C^2` times

\[
 \det\begin{pmatrix}
 -2A&A_U&A_V\\
 -B&B_U&B_V\\
 C&C_U&C_V
 \end{pmatrix},
\]

which is `det DF` by the weighted invariant-Jacobian reduction.  Thus every
nonconstant `C` forces ramification in the plane quotient even when the
threefold map is Keller.  Replacing `(U,V)` or `(P,Q)` by other polynomial
generators of the same invariant rings only conjugates by plane polynomial
automorphisms, whose Jacobians are nonzero constants; it cannot remove the
factor `C^2`.

For (1),

\[
\begin{aligned}
 A&=(1+U)^3V+U^2(1+U)(4+3U),\\
 B&=U+3(1+U)^2V+3U^2(4+3U),\\
 C&=2-3U-V.
\end{aligned}
 \tag{7}
\]

Consequently

\[
 \det D\overline F=-2(2-3U-V)^2.
 \tag{8}
\]

The whole line `C=0` is sent to `(P,Q)=(0,0)`.  The original three-point
collision descends to

\[
 (U,V)=(0,0),\quad(-3/2,13/2)
 \longmapsto(0,0),
 \tag{9}
\]

because the two nonzero source points differ by the torus element `-1`.
Thus the torus quotient retains noninjectivity, but only as a visibly
non-Keller plane map with a contracted line.

## 2. The vertical locally nilpotent derivation does not descend

The normalized factorization model supplies the free locally nilpotent
derivation `D=partial_z`, with

\[
 \ker D=k[x,y].
\]

In fact it admits no nonconstant target first integral:

\[
 \boxed{F^*k[a,b,c]\cap k[x,y]=k.}
 \tag{10}
\]

To prove this, work over the generic target field `K=k(a,b,c)` and put

\[
 t=y+1/x=u/x.
\]

The generic inverse equation is

\[
 f(T)=cT^3-2T^2+bT-2a=0.                         \tag{11}
\]

Also

\[
 \partial_zF=x^3(t^3,3t^2,-1).                   \tag{12}
\]

If `h in k[a,b,c]` and `partial_z(h circ F)=0`, then in `K[T]/(f)`

\[
 h_aT^3+3h_bT^2-h_c=0.                            \tag{13}
\]

The generic cubic (11) is the minimal polynomial of `t`, so the left side
of (13) is a scalar multiple `lambda f`.  Comparing the coefficient of `T`
gives `lambda b=0`, hence `lambda=0` in `K`.  Therefore all three partial
derivatives of `h` vanish, and characteristic zero gives `h in k`.

It follows that no dominant target map to `A^2` can make a commuting quotient
square with the source quotient `(x,y)`.  The same conclusion holds for every
nonzero replica `g(x,y)partial_z`, since it has the same kernel.

There is also no quotient square made from rank-two **linear** projections.
If their one-dimensional kernels are spanned by constant vectors `v` and
`w`, commutativity would force

\[
 DF(x,y,z)v\in k[x,y,z]w.                          \tag{14}
\]

The coefficient matrix of `DFv`, viewed as three rows indexed by the target
coordinates, has rank at most one exactly when all its `2 x 2` minors vanish.
For `v=(p,q,r)`, those minors generate

\[
 (p^2,pq,q^2,pr,qr,r^2),                           \tag{15}
\]

whose only point is `v=0`.  This excludes nonzero constant source directions.

## 3. An invariant affine plane exists, but gives an automorphism

The hypersurface `x=0` maps into the target hypersurface `c=0`.  Restriction
of (1) gives

\[
 F|_{x=0}:\mathbb A^2_{y,z}\longrightarrow\mathbb A^2_{a,b},
 \qquad(y,z)\longmapsto(z+4y^2,y).                 \tag{16}
\]

This has determinant `-1` and polynomial inverse

\[
 (a,b)\longmapsto(b,a-4b^2).
\]

So the construction does contain a preserved `A^2` carrying a Keller map,
but the map is triangular and invertible.  It contains only the distinguished
point of the collision fiber.

Moreover this plane is not the full pullback of the target plane:

\[
 F^{-1}(c=0)=V\bigl(x(2-3xy-x^2z)\bigr).
 \tag{17}
\]

The second component has coordinate ring

\[
 k[x,x^{-1},y],
\]

so it is `G_m x A^1`, not `A^2`.  This extra component is where the other
two points of the foundational collision lie.

## 4. What remains open

For the candidates naturally supplied by the construction, the outcomes are
therefore exact:

| candidate | quotient/source | outcome |
|---|---|---|
| residual torus | smooth `A^2` quotient | descended Jacobian `-2C^2`; a line is contracted |
| vertical `G_a` action | `A^2_(x,y)` | no nonconstant target first integral |
| rank-two linear projections | `A^2` | no compatible nonzero constant direction |
| hypersurface `x=0` | `A^2` | triangular plane automorphism |

This does not exclude a nonlinear algebraic foliation, a different LND, or a
rational equivalence relation.  A direct search can be formulated without
choosing vector fields.  Seek algebraically independent pairs
`pi_s=(r,s)` and `pi_t=(p,q)` such that

\[
 p\circ F=G_1(r,s),\qquad q\circ F=G_2(r,s),
 \tag{18}
\]

and impose the differential identity

\[
 d(p\circ F)\wedge d(q\circ F)
 =\lambda\,dr\wedge ds,\qquad\lambda\in k^*.       \tag{19}
\]

Equations (18)--(19), with bounded supports compatible with the torus
filtration but not restricted to weight zero, give a finite coefficient
scheme for the next nonlinear search.  Any successful solution must evade
the square factor (6), the first-integral obstruction (10), and the linear
direction obstruction (15).

The exact symbolic audit is
[`verify_foundational_plane_quotients.py`](../scripts/verify_foundational_plane_quotients.py).
