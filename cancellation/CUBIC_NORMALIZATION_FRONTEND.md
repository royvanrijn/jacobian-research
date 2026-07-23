# Cubic normalization and the affine-linear frontend

This note moves the degree-three classification problem one step upstream.
Instead of assuming a coordinate-preserving suspension, it starts with the
canonical finite normalization of a cubic Keller map.  It proves that the
only possible failure of flatness is zero-dimensional, extracts the
Deligne--Faddeev binary cubic when that defect is absent, and classifies the
case in which the resulting coefficient morphism is affine-linear.

The result does **not** make the proposed minimal-boundary gateway
operational or prove the resulting classification conjecture.  It replaces
one part of the undifferentiated phrase “extract a suspension” by two
concrete intrinsic obligations:

1. remove the zero-dimensional flatness defect of the finite normalization;
2. prove that its binary-cubic coefficient morphism is affine-linear (or
   prove a nonlinear replacement for the hyperplane-slice theorem).

Work over an algebraically closed field `k` of characteristic zero.

## 1. The canonical cubic algebra

Let

\[
 F:U=\mathbb A^3\longrightarrow Y=\mathbb A^3
\]

be a dominant quasi-finite Keller map of geometric degree three.  Put

\[
 A=k[Y],\qquad L=k(U),\qquad
 B=\operatorname{Norm}_A(L),
\tag{1.1}
\]

and write

\[
 \pi:\bar X=\operatorname{Spec}B\longrightarrow Y
\]

for the canonical finite normalization.  The distinguished map
`U -> bar X` is the Zariski--Main open immersion.

The `A`-module `B` is finite, torsion-free, and has generic rank three.
It need not be locally free merely because `B` is normal: normality is an
`S_2` condition, while flatness over the regular threefold `Y` is a
maximal-Cohen--Macaulay condition.

### Proposition 1.1 -- flatness is a codimension-three question

The finite morphism `pi` is flat over every point of `Y` of codimension at
most two.  Consequently its nonflat locus is a finite zero-dimensional
closed subset of `Y`.

Equivalently, since `B` has generic rank three,

\[
 Z_{\mathrm{flat}}
 =
 V\!\left(\operatorname{Fitt}^{A}_{3}(B)\right)
\tag{1.2}
\]

is zero-dimensional, and `pi` is finite flat exactly when
`Z_flat` is empty.

#### Proof

Let `p in Spec A` have height at most two and let `q` be a prime of `B`
above it.  Integrality and equality of fraction-field dimensions give

\[
 \dim B_q=\dim A_p\le2.
\]

The local ring `B_q` is normal.  A normal local ring of dimension at most
two is Cohen--Macaulay.  Since `A_p` is regular and the map is finite,
miracle flatness makes `B_q` flat over `A_p`.  A finite flat module is
locally free, proving the first assertion.

For a finite module of generic rank three, the locally free rank-three
locus is the complement of `V(Fitt_3(B))`; the lower Fitting ideal
`Fitt_2(B)` is zero because it vanishes after tensoring with `Frac(A)` and
`A` is a domain.  The first assertion puts the support of (1.2) in
codimension three.  Since `A` has dimension three, that support is finite.
QED

This is important for the minimal-boundary program.  Height-one valuations,
unit lattices, and divisor ledgers cannot by themselves see
`Z_flat`.  The finite normalization itself does see it through (1.2), so
flatness is an intrinsic finite-stratum condition rather than a new
coordinate marking.

### Proposition 1.2 -- a Cartier boundary removes the point defect

Suppose every point of `bar X\setminus U` has a neighborhood on which the
boundary is an effective Cartier divisor `(t=0)` and `B/(t)` is
Cohen--Macaulay.  Then `Z_flat` is empty.

In particular, it is enough that the canonical boundary be Cartier and
normal as a surface.

#### Proof

The map is already étale, hence flat, at every point of `U`.  Let `q` be a
boundary point.  The Cartier equation `t` is a nonzerodivisor in the
three-dimensional local domain `B_q`.  Since `B_q/tB_q` is
Cohen--Macaulay of dimension two, a regular sequence of length two modulo
`t` lifts to a regular sequence of length three beginning with `t`.
Therefore `B_q` is Cohen--Macaulay.  Proposition 1.1, or miracle flatness
directly, makes it flat over the corresponding regular target local ring.

A normal surface is Cohen--Macaulay, giving the last assertion.  QED

This criterion shows exactly what a strengthened saturation theorem would
need to prove.  Generic DVR saturation is not enough: it must extend the
boundary equation as a Cartier class through the finitely many closed
points and control the scheme-theoretic surface there.

### Proposition 1.3 -- the binary-cubic package

Assume `Z_flat` is empty.  Then:

1. `B` is a finite locally free `A`-algebra of rank three;
2. the trace splitting
   \[
   B=A\oplus M,\qquad
   M=\ker(\operatorname{Tr}_{B/A}),
   \tag{1.3}
   \]
   makes `M` a projective module of rank two;
3. `M` is free by Quillen--Suslin;
4. after choosing a basis of `M`, the cubic algebra is represented by a
   binary cubic
   \[
   f_Y(U,V)
   =aU^3+bU^2V+cUV^2+dV^3,
   \qquad a,b,c,d\in A,
   \tag{1.4}
   \]
   uniquely up to the natural `GL_2(A)` action.

Thus the normalization canonically supplies, up to change of
Tschirnhausen basis, a coefficient morphism

\[
 \kappa_F:Y\longrightarrow
 \operatorname{Sym}^3(k^2)\simeq\mathbb A^4.
\tag{1.5}
\]

#### Proof

Finite flatness gives (1).  Since `3` is invertible, the unit map
`A -> B` is split by `(1/3)Tr`, giving (1.3).  Direct summands of a
finite locally free module are projective, so `M` has rank two.
Quillen--Suslin makes it free over the polynomial ring `A`.

The Deligne--Faddeev correspondence over an arbitrary base identifies
rank-three locally free commutative algebras with twisted binary cubic
forms on the rank-two quotient by the unit.  A basis of `M` untwists the
form and gives (1.4); changing the basis gives exactly the stated
`GL_2(A)` action.  QED

The choice of basis in (1.4) is not intrinsic.  The pair `(M,f_Y)` and its
polynomial `GL_2(A)` orbit, the flatness-defect scheme (1.2), and the
discriminant divisor are intrinsic.  The affine span of the four displayed
coefficients in one basis is **not** intrinsic: a target-dependent
Tschirnhausen shear can make a linear representative look nonlinear.  The
invariant condition used below is the existence of an affine-linear
full-rank representative somewhere in the polynomial `GL_2(A)` orbit.

### Proposition 1.4 -- the cubic critical-divisor degree budget

Let `E` be a boundary prime of the canonical normalization which is ramified
over its target prime divisor `D`.  Then the generic DVR decomposition over
`D` consists of exactly:

1. `E`, with ramification index two and residue degree one;
2. one prime contained in the Keller open `U`, with ramification index and
   residue degree both one.

In particular, no additional boundary prime lies over `D`.

#### Proof

Let `g` be a prime equation of `D` in `k[Y]`.  Dominance makes
`g\circ F` a nonconstant polynomial on `U=A^3`, hence a nonunit.  Any
height-one prime factor of `(g\circ F)` maps dominantly to `D`: its source
divisor has dimension two and the restriction of the quasi-finite map still
has zero-dimensional fibers.  This supplies an affine prime above `D`.
The Keller condition makes it unramified, so its contribution `ef` to the
generic DVR degree sum is at least one.

The critical boundary prime contributes at least two.  Since the function
field degree is three,

\[
 \sum_{\mathfrak q\mid D}e(\mathfrak q/D)f(\mathfrak q/D)=3.
\]

Both lower bounds are therefore equalities.  The critical contribution is
`2*1`, the affine contribution is `1*1`, and there is no room for another
prime.  QED

Consequently, under the hypothesis that there is only one critical boundary
prime, every additional boundary prime is unramified and maps to a target
divisor distinct from `D`.  The extra-simple-boundary obstruction in degree
three is therefore exactly the existence of a second, purely unramified
nonproperness divisor; it cannot hide over the critical discriminant.

### Proposition 1.5 -- point-flatness is cubic fiber-minimality

For every `p in Y`, let

\[
 \lambda(p)=
 \dim_{\kappa(p)}\!\left(B\otimes_A\kappa(p)\right),
\tag{1.6}
\]

the scheme-theoretic length of the finite normalization fiber.  Then

\[
 \pi\text{ is flat at }p
 \quad\Longleftrightarrow\quad
 \lambda(p)=3.
\tag{1.7}
\]

In particular,

\[
 |Z_{\mathrm{flat}}|
 =
 \{p\in Y:\lambda(p)\ge4\}.
\tag{1.8}
\]

Thus the point-flatness obligation can be replaced by the intrinsic
**cubic fiber-minimality** condition `lambda(p)<=3` at the finitely many
closed target points under the boundary.

#### Proof

Localize at `p` and write `R=A_p`, `N=B_p`, and `K=Frac(R)`.  The
torsion-free `R`-module `N` has generic rank three.  Hence every generating
set of `N` has at least three elements, while Nakayama's lemma identifies
the least number of generators with

\[
 \mu_R(N)=\dim_{\kappa(p)}N/\mathfrak m_pN=\lambda(p).
\]

If `lambda(p)=3`, choose a minimal surjection `R^3 -> N`.  After tensoring
with `K` it is a surjection between three-dimensional vector spaces and is
therefore an isomorphism.  Its kernel has generic rank zero, but it is a
submodule of the torsion-free module `R^3`; consequently the kernel is zero.
Thus `N` is free of rank three.  The converse is immediate.

Since `lambda(p)` is always at least three, failure of flatness is
equivalent to `lambda(p)>=4`.  Proposition 1.1 shows that these points are
closed and finite.  QED

Proposition 1.5 is strictly scheme-theoretic.  A bound on the number of
geometric places, or even on the number of reduced points in a fiber, does
not bound `lambda(p)`: multiplicity can create an excess-length fiber
without creating another place.  What suffices is a length-three statement
for the canonical finite fiber itself.  This is weaker than the
Cartier--Cohen--Macaulay atlas in Proposition 1.2 and is the smallest direct
closed-point input that removes the Fitting obstruction.

### Proposition 1.6 -- cubic algebra structure alone does not force flatness

There are integral normal degree-three covers of a smooth factorial
threefold whose trace-free module is reflexive but not locally free.
Consequently neither normality nor the existence of the commutative
rank-three algebra structure eliminates the point defect in Proposition
1.1.

#### Proof

On `X=A^3=Spec k[x,y,z]`, consider the module `M` with minimal presentation

\[
 0\longrightarrow A
 \mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}
 A^3\longrightarrow M\longrightarrow0.
\tag{1.9}
\]

It is the second Koszul syzygy of `(x,y,z)`.  Away from the origin one entry
of the relation is a unit, so `M` is free of rank two.  At the origin the
Auslander--Buchsbaum formula gives depth two; hence `M` is torsion-free and
`S_2`, and therefore reflexive.  Its determinant is trivial.  The relation
vanishes in the origin fiber, so

\[
 \dim_k(M/\mathfrak mM)=3.
\tag{1.10}
\]

To apply Tan's nonflat triple-cover existence theorem without suppressing
its genericity hypothesis, compactify this module.  On
`\mathbb P^3=Proj k[w,x,y,z]`, let `\mathcal M` be the cokernel of

\[
 0\longrightarrow\mathcal O_{\mathbb P^3}(-1)
 \mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}
 \mathcal O_{\mathbb P^3}^{\,3}\longrightarrow\mathcal M\longrightarrow0.
\]

The map vanishes only at `[1:0:0:0]`.  The same depth computation as above
shows that `\mathcal M` is rank-two reflexive, and its restriction to the
chart `w\ne0` is `M`.  Tan's Theorem 7.6 constructs a reduced irreducible
triple cover from a **general** member of the basepoint-free cubic system on
the projectivization of a rank-two reflexive sheaf.  Its projective
corollary says that, after tensoring by a sufficiently negative power of a
very ample line bundle, any rank-two reflexive sheaf on a projective
factorial variety occurs as the trace-free sheaf of such a cover.  Apply
that corollary to `\mathcal M`.  Restriction to `w\ne0` removes the harmless
line-bundle twist and gives an integral normal triple cover of `X` with

\[
 \pi_*\mathcal O_{\bar X}\simeq\mathcal O_X\oplus M.
\tag{1.11}
\]

Its origin fiber has length `1+3=4`, so Proposition 1.5 makes it nonflat
there.  Equivalently, the presentation of the right side of (1.11) has
relation `(0,z,-y,x)`, and its rank-three Fitting ideal is exactly
`(x,y,z)`.  QED

The word **general** is essential here.  Global generation supplies a
basepoint-free linear system; it does not make an arbitrarily prescribed
cubic tensor reduced, irreducible, or normal.  Nor may one normalize an
arbitrary cubic construction and assume without proof that normalization
preserves the displayed trace module.  The existence assertion above uses
Tan's generic reduced-irreducible divisor and his reflexive-extension
argument.  Any search with a prescribed initial cubic must reprove the
corresponding genericity statement inside that constrained subsystem.

This example is not asserted to contain a Keller open isomorphic to `A^3`
or to satisfy the minimal boundary ledger.  It has a narrower purpose: it
rules out any proof of point-flatness that uses only normal finite
degree-three algebra structure.  A successful proof must use the
distinguished Keller open, the scheme intersections of its boundary, or an
equivalent fiber-minimality input.

### Proposition 1.7 -- determinantal classification of every point defect

Let `p` be a point of `Z_flat`, put `R=A_p`, and write the trace splitting

\[
 B_p=R\oplus M.
\tag{1.12}
\]

There is a unique integer `s>=1` and a minimal free resolution

\[
 0\longrightarrow R^s
 \mathop{\longrightarrow}^{\Phi}
 R^{s+2}\longrightarrow M\longrightarrow0
\tag{1.13}
\]

such that

\[
 \lambda(p)=s+3,\qquad
 \operatorname{Fitt}_3^R(B_p)
 =\operatorname{Fitt}_2^R(M)
 =I_s(\Phi).
\tag{1.14}
\]

The maximal-minor ideal `I_s(Phi)` is primary to the maximal ideal of `R`.
In particular, the defect has its smallest possible fiber length four
exactly when

\[
 0\longrightarrow R
 \mathop{\longrightarrow}^{(a,b,c)^{\mathsf T}}
 R^3\longrightarrow M\longrightarrow0,
\tag{1.15}
\]

where `(a,b,c)` is a parameter ideal.  If the Fitting defect is reduced at
`p`, then `(a,b,c)` is the maximal ideal; after completion and a regular
change of parameters, (1.15) is the Koszul model `(x,y,z)^T`.

In fact, **every reduced point defect is automatically this minimal
length-four case**.  For `s>=2`, minimality puts the entries of `Phi` in the
maximal ideal and therefore

\[
 I_s(\Phi)\subseteq\mathfrak m^s\subseteq\mathfrak m^2.
\]

An ideal primary to `mathfrak m` and contained in `mathfrak m^2` cannot
equal its radical `mathfrak m`.  Thus all higher determinantal rungs, and
all non-Koszul `s=1` parameter ideals, are necessarily nonreduced defects.

#### Proof

The finite normal ring `B_p` is torsion-free over `R` and has depth at least
two: at every prime above `p`, normality supplies `S_2`, and depth for a
finite module is the minimum of the depths at those primes.  The trace
summand `M` is therefore a reflexive rank-two `R`-module.

The regular local ring `R` has dimension three and finite global dimension.
Auslander--Buchsbaum gives

\[
 \operatorname{pd}_R(M)+\operatorname{depth}_R(M)=3.
\]

If `M` had projective dimension zero it would be free and `p` would be
flat.  Thus `depth_R(M)=2` and `pd_R(M)=1`.  Its minimal resolution has the
form (1.13), because the difference of the two free ranks is
`rank_R(M)=2`.

Minimality puts every entry of `Phi` in the maximal ideal.  Hence

\[
 \dim_{\kappa(p)}M/\mathfrak mM=s+2
\]

and the unit summand gives `lambda(p)=s+3`.  The direct-sum rule for Fitting
ideals and the presentation (1.13) give (1.14).  Proposition 1.1 says its
vanishing locus is supported only at the closed point, so `I_s(Phi)` is
maximal-ideal-primary.

For `s=1`, this ideal is `(a,b,c)`.  Three generators of a height-three
ideal in the Cohen--Macaulay ring `R` form a parameter regular sequence.
If the ideal is reduced and primary to the maximal ideal, it equals that
maximal ideal.  Its generators are then a regular system of parameters,
giving the completed Koszul normal form.  The preceding containment proves
that a reduced defect cannot have `s>=2`.  QED

This reduces the closed-point part of the conjecture to a concrete
determinantal exclusion problem.  There is no hidden infinite homological
complex: the excess fiber length `lambda(p)-3` is exactly the number of
relations `s`.  The first unresolved local case is therefore the single
Koszul relation (1.15), with compatibility imposed by the distinguished
Keller open and the one-boundary scheme intersections.

### Proposition 1.8 -- every reduced defect has a square-zero fiber

Assume the Fitting defect at `p` is reduced.  Then it automatically has
length four, and its canonical fiber algebra is uniquely

\[
 B\otimes_A\kappa(p)
 \simeq
 \kappa(p)\oplus V,
 \qquad
 \dim_{\kappa(p)}V=3,
 \qquad
 V^2=0.
\tag{1.16}
\]

Consequently the fiber has one geometric point, that point lies in the
boundary, and

\[
 F^{-1}(p)=\varnothing.
\tag{1.17}
\]

#### Proof

By Proposition 1.7, reducedness forces `s=1`; after completion and a regular
change of parameters,

\[
 M=\operatorname{coker}
 \left(R\mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}R^3\right).
\]

Write the product of two trace-free elements as a scalar component and an
`M` component.

First consider an `R`-linear functional `M -> R`.  It lifts to a row
`w=(w_1,w_2,w_3)` satisfying

\[
 w_1z-w_2y+w_3x=0.
\]

Reducing modulo the square of the maximal ideal, the linear independence of
`x,y,z` forces the constant row of `w` to vanish.  Thus every such
functional has image in the maximal ideal, and the scalar component of
`V*V` vanishes in the closed fiber.

Next let `T` be an endomorphism of `M`.  It lifts to a matrix
`C in Mat_3(R)` satisfying

\[
 C(z,-y,x)^{\mathsf T}
 =r(z,-y,x)^{\mathsf T}
\]

for some `r in R`.  Reduction modulo the square of the maximal ideal forces
the constant matrix of `C` to be the scalar matrix `r(0)I`.  Therefore, for
each `u in V`, the `M` component of multiplication by `u` has the form

\[
 v\longmapsto\ell(u)v
\]

for a linear functional `ell:V -> kappa(p)`.  Commutativity gives

\[
 \ell(u)v=\ell(v)u
\quad\text{for all }u,v\in V.
\]

Since `dim V=3`, this identity forces `ell=0`.  Both components of `V^2`
therefore vanish, proving (1.16).

The square-zero algebra in (1.16) is local and has a single geometric
point.  If that point lay in the Keller open, the fiber local algebra would
be reduced of length one because `F` is étale.  Hence it lies in the
boundary, and there is no other fiber point in `U`, proving (1.17).  QED

### Proposition 1.8a -- the first multiplication symbol is a ternary cubic

Retain the completed reduced-defect model

\[
 R=k[[x,y,z]],\qquad
 M=\operatorname{coker}\left(
 R\mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}R^3
 \right).
\tag{1.16a}
\]

A nonflat triple-cover multiplication with trace module `M` is encoded by
the generalized Miranda--Tan tensor

\[
 \phi:\operatorname{Sym}^3(M)\longrightarrow\det(M)\simeq R.
\tag{1.16b}
\]

Choose the displayed generators `e_1,e_2,e_3` of `M` and write

\[
 c_{ijk}=\phi(e_i,e_j,e_k).
\]

Then every `c_ijk` lies in `m^3`.  The degree-three initial symbols form a
ten-dimensional vector space.  More precisely, put

\[
 r=(z,-y,x)\in R_1\otimes k^3.
\]

After choosing a volume form on `k^3`, every degree-three symbol is uniquely
of the form

\[
 \boxed{
 \operatorname{in}_3(\phi)(u,v,w)
 =
 h(r\times u,r\times v,r\times w)
 }
\tag{1.16c}
\]

for a ternary cubic `h`, where the right side means the symmetric
trilinear polarization of `h`.  Thus, up to the harmless determinant twist,
the first multiplication-symbol space is

\[
 \operatorname{Sym}^3(k^3)^\vee.
\tag{1.16d}
\]

#### Proof

The relation in (1.16a) makes well-definedness of `phi` equivalent to

\[
 zc_{1jk}-yc_{2jk}+xc_{3jk}=0
\tag{1.16e}
\]

for every symmetric pair `(j,k)`.  Taking homogeneous degree `d` symbols,
the solution space is the kernel of the natural contraction map

\[
 K_d=
 \ker\left[
 \operatorname{Sym}^d(V^\vee)\otimes\operatorname{Sym}^3(V^\vee)
 \longrightarrow
 \operatorname{Sym}^{d+1}(V^\vee)
 \otimes\operatorname{Sym}^2(V^\vee)
 \right],
\tag{1.16f}
\]

where `V=k^3`.  The Pieri decomposition shows that the displayed map has
zero kernel for `d=0,1,2`.  For `d>=3`, its kernel is the Schur module

\[
 K_d\simeq\mathbb S_{(d,3)}(V^\vee).
\tag{1.16g}
\]

In particular,

\[
 \dim K_3=\dim\mathbb S_{(3,3)}(V^\vee)=10.
\]

For every ternary cubic `h`, formula (1.16c) annihilates `r` in each
argument and therefore satisfies (1.16e).  Polarizing the ten ternary
monomials gives ten linearly independent elements of `K_3`.  Since `K_3`
has dimension ten, these exhaust it and prove (1.16c)--(1.16d).  QED

This removes a large false search space.  At a reduced defect, one should
not enumerate arbitrary bilinear multiplication tables and then impose
associativity: the generalized triple-cover correspondence already packages
associativity, and its first nonzero local datum is just a ternary cubic.
The initial counterexample search therefore begins with the classical
plane-cubic types

\[
 \text{smooth},\quad
 \text{nodal},\quad
 \text{cuspidal},\quad
 \text{reducible}
\tag{1.16h}
\]

of plane cubics, followed by higher-order lifts of (1.16c).

These are geometric types, not literally four `GL_3` orbits.  The smooth
type has its `j`-modulus, line-plus-conic splits into transverse and tangent
types, and three distinct lines split into a triangle and three concurrent
lines.  A useful working list of normal forms is

| type | normal form or family |
|---|---|
| smooth | `X^3+Y^3+Z^3-3 lambda XYZ`, `lambda^3 ne 1` |
| nodal | `Y^2Z-X^2(X+Z)` |
| cuspidal | `Y^2Z-X^3` |
| line plus transverse conic | `Z(XY-Z^2)` |
| line tangent to conic | `Z(YZ-X^2)` |
| three lines forming a triangle | `XYZ` |
| three concurrent lines | `XY(X-Y)` |
| double line plus line | `X^2Y` |
| triple line | `X^3` |

### Proposition 1.8b -- the cubic is the exceptional incidence curve

Blow up the defect point in `Spec R` and let
`E\simeq\mathbb P^2` be the exceptional plane.  After removing pullback
torsion, the transform of `M` on `E` is the universal quotient

\[
 0\longrightarrow\mathcal O_E(-1)
 \longrightarrow\mathcal O_E^{\,3}
 \longrightarrow Q\longrightarrow0.
\tag{1.16i}
\]

For `[r]\in E`, cross product identifies

\[
 Q_{[r]}=k^3/\langle r\rangle
 \mathop{\longrightarrow}^{\sim}r^\perp,
 \qquad [u]\longmapsto r\times u.
\tag{1.16j}
\]

Consequently the exceptional cubic divisor supplied by the order-three
symbol `h` is

\[
 Z_h=
 \left\{
 ([r],[q])\in\mathbb P^2\times\mathbb P^2:
 r\cdot q=0,\ h(q)=0
 \right\}.
\tag{1.16k}
\]

It is the pullback of the plane cubic `C_h=V(h)` along the second projection
of the incidence flag variety.  In particular, scheme-theoretically,
`Z_h\to C_h` is a `\mathbb P^1`-bundle.  Thus reducedness, irreducibility,
and the number of components of this first exceptional divisor are exactly
those of `C_h`.

#### Proof

On the blowup chart the Koszul relation has a common exceptional parameter.
Saturating removes that parameter and leaves the tautological line
`\langle r\rangle\subset k^3`, giving (1.16i).  Cross product has kernel
`\langle r\rangle` and image `r^\perp`, proving (1.16j).  On the diagonal
`u=v=w`, formula (1.16c) becomes `h(r\times u)`.  Under (1.16j), its zero
scheme is therefore (1.16k).  The incidence projection has fiber the line
`r\cdot q=0` over every `[q]`, proving the last assertion.  QED

This gives a first local ranking, but not a classification.  Smooth `h`
produces a ruled exceptional surface over an elliptic curve.  Nodal and
cuspidal `h` are irreducible with rational normalization.  Reducible `h`
produces several exceptional components, and a double or triple line gives
a nonreduced initial exceptional divisor.  None of these facts alone decides
normality of the original finite threefold: an exceptional divisor on a
modification may be nonnormal even when the unmodified cover is normal.

Higher lifts have the form

\[
 \phi=\phi_h+\phi_4+\phi_5+\cdots,\qquad
 \phi_d\in K_d\simeq\mathbb S_{(d,3)}(k^3)^\vee.
\tag{1.16l}
\]

Tan's unconstrained generic-section theorem does not by itself say what
happens after `h` is fixed.  In the present punctured affine model, however,
the constrained genericity statement follows from the same method.

### Proposition 1.8c -- every nonzero cubic type has normal integral lifts

For every nonzero ternary cubic `h`, there is a generalized cubic tensor
`\phi` with initial symbol `h` whose reflexive triple-cover extension is
integral and normal and has trace module `M`.

#### Proof

Work first with the algebraic Koszul module `M_A` over
`A=k[x,y,z]`; completion at the origin recovers (1.16a).  Put

\[
 T_A=
 \operatorname{Hom}_A(\operatorname{Sym}^3M_A,\det M_A)
\]

and let `m=(x,y,z)`.  Formula (1.16c) supplies a polynomial tensor
`\phi_h` with initial symbol `h`.  Every perturbation in `mT_A`
has order at least four and therefore preserves that initial symbol.

On the punctured affine threefold `U=Spec A\setminus\{m\}`, one has
`m\mathcal O_U=\mathcal O_U`.  Hence

\[
 \widetilde{mT_A}|_U=\widetilde{T_A}|_U.
\tag{1.16m}
\]

The module `M_A|_U` is locally free, and `\widetilde{T_A}|_U` is the rank-four
bundle of ordinary binary cubics.  Quasi-compactness permits a
finite-dimensional subspace `W\subset mT_A` which generates this
bundle on `U`.  Enlarging `W` by coordinate multiples makes its induced
linear system on `\mathbb P(M_A|_U)` separate the base directions as well as
restrict to the complete cubic system on every projective-line fiber.

For general `\psi\in W`, the section `\phi_h+\psi` is nowhere the zero
binary cubic on a fiber: vanishing imposes four independent linear
conditions over a three-dimensional base.  Its divisor in
`\mathbb P(M_A|_U)` is therefore finite flat of degree three over `U`.
Characteristic-zero Bertini, applied to this basepoint-free system which is
not composed with a pencil, makes the general divisor smooth and
irreducible.

Tan's codimension-two extension construction now extends this punctured
triple cover to the normalization over `Spec A`.  The extension is integral
and normal.  Its trace-free module and `M_A` agree on `U`; both are
reflexive, so they agree on `Spec A`.  Finally, the perturbation lies in
`mT_A`, so its order-three symbol remains `h`.  Localizing and completing
at the origin gives the asserted cover of the completed model; normality is
preserved because these finite-type local rings are excellent.  QED

Thus initial cubic type alone cannot be excluded by abstract existence of
a normal integral cover.  What remains nontrivial is compatibility of such
a lift with a single global Keller open and its marked boundary.  A
bounded-degree computational ansatz must also verify that it contains
enough of the perturbation space used above; genericity in the full module
does not imply genericity in an arbitrarily small ansatz.

Only after this constrained-lifting gate should a candidate be tested for:

1. preservation of normality inside the chosen global or bounded-degree
   ansatz (`S_2` is already forced by the underlying module `R\oplus M`);
2. a localization isomorphic to a localization of `k[u,v,w]`;
3. vanishing of `Omega_{B/A}` on that affine open;
4. the exact DVR splitting `(e,f)=(2,1)+(1,1)` over the critical divisor;
5. absence of every other boundary prime whose image is a target divisor.

For constructing an example, the nodal and cuspidal types are the first
rational candidates; the smooth family is the cleanest constrained-Bertini
test and may instead expose an elliptic-boundary obstruction.  Reducible
and nonreduced types should follow only after the irreducible cases, because
their extra exceptional components enlarge rather than simplify the
boundary-prime ledger.

There is also a geometric warning.  Ordinary trace and discriminant
identities are built into (1.16b), and the tame different/log-crepancy
ledger is codimension-one automatic.  They cannot eliminate the Koszul
defect.  The genuinely Keller-specific gates begin only after (1.16c):
rationality of the cubic function field, the colored
`(2,1)+(1,1)` boundary decomposition, and an open complement isomorphic to
`A^3`.  If such an `A^3` complement exists and the finite map is étale on
it, its Jacobian is automatically a nonzero constant because every unit of
`k[A^3]` is constant.

### Corollary 1.9 -- reduced point-flatness becomes sheet separation

Let `D` be the target divisor of the critical boundary prime and let
`A_D` be the closure in `bar X` of the unique affine prime over `D` from
Proposition 1.4.  A reduced length-four defect over a closed point `p in D`
forces the critical boundary prime and `A_D` to meet at the unique point of
the fiber over `p`.

Hence any one of the following excludes reduced minimal defects over `D`:

1. `F(A^3)` contains every closed point of `D`;
2. the critical boundary and `A_D` are separated over closed points of `D`;
3. their intrinsic scheme intersection has no fiber whose local algebra is
   the square-zero algebra in (1.16).

#### Proof

Both the critical boundary and `A_D` are finite and dominant over `D`, hence
surjective.  Proposition 1.8 says the fiber over `p` has only one point.
Both closures must contain that point.  Each of the three stated conditions
prevents the square-zero collision.  QED

This is the first direct bridge from the point-flatness defect to the
intrinsic boundary markings.  Generic height-one saturation supplies the
two sheets but cannot say whether they collide at a closed point.  The
scheme-intersection part of the Zariski--Main package can say exactly that.
The correct target is not unconditional sheet separation.  In the flat
foundational model the two sheets meet over a triple-root cubic, but the
fiber there is the curvilinear length-three algebra

\[
 \kappa(p)[\epsilon]/(\epsilon^3),
\tag{1.18}
\]

not the embedding-dimension-three square-zero algebra (1.16).  The remaining
local task is therefore to prove that the intrinsic intersection is
curvilinear/fiber-minimal.  Proposition 1.8 handles every reduced defect;
only nonreduced Fitting defects remain outside the square-zero analysis.

### Proposition 1.10 -- curvilinear fibers force point-flatness

Let `p in Y`.  If the finite algebra

\[
 B\otimes_A\kappa(p)
\tag{1.19}
\]

is generated by one element as a `kappa(p)`-algebra, then `B_p` is free of
rank three over `A_p`.  Consequently, if every closed canonical fiber is
curvilinear, then `Z_flat` is empty.

Here a finite fiber is called curvilinear when each of its local Artin
factors has embedding dimension at most one.  At the closed points relevant
to Proposition 1.1, `kappa(p)=k` is algebraically closed and infinite; a
finite product of such factors is generated by one element, so this
geometric definition implies the hypothesis above.

#### Proof

Put `R=A_p` and choose `bar t` generating `B_p/\mathfrak m_pB_p` as a
residue algebra.  Lift it to `t in B_p` and set `S=R[t]`.  The element `t`
is integral, so `S` is a finite `R`-module.  The generator assumption says

\[
 B_p=S+\mathfrak m_pB_p.
\]

Apply Nakayama's lemma to the finite module `B_p/S`; it gives `B_p=S`.

Let `K=Frac(R)`.  Since `B_p` has generic rank three and is a domain,
the minimal polynomial of `t` over `K` has degree three.  Its coefficients
are integral over `R` and lie in `K`; the regular local ring `R` is
integrally closed, so they lie in `R`.  If this monic cubic is `f(T)`, then
division by `f` shows

\[
 B_p=R[t]\simeq R[T]/(f),
\]

which is free with basis `1,t,t^2`.  The last assertion follows from
Proposition 1.1, since only closed points could be nonflat.  QED

This is the cleanest current point-flatness target.  The stronger
Cartier--Cohen--Macaulay atlas of Proposition 1.2 is sufficient but not
necessary.  It is enough to extract from the intrinsic scheme intersections
that every collision fiber is curvilinear.  The foundational triple-root
fiber (1.18) satisfies this condition, whereas every reduced defect fails it
maximally by Proposition 1.8.  Proposition 1.10 also excludes all
nonreduced and higher determinantal defects without classifying their
individual algebra structures.

### Proposition 1.11 -- the intrinsic cotangent-cyclicity test

For a closed point `p in Y`, put

\[
 C_p=B\otimes_A\kappa(p).
\]

The following conditions are equivalent:

1. the finite fiber `Spec C_p` is curvilinear;
2. `Omega_{C_p/kappa(p)}` is locally generated by at most one element;
3.
   \[
   \operatorname{Fitt}^{C_p}_1
   \left(\Omega_{C_p/\kappa(p)}\right)=C_p;
   \tag{1.20}
   \]
4.
   \[
   \bigwedge\nolimits^2\Omega_{C_p/\kappa(p)}=0.
   \tag{1.21}
   \]

Moreover,

\[
 \Omega_{C_p/\kappa(p)}
 \simeq
 \Omega_{B/A}\otimes_B C_p,
\tag{1.22}
\]

so (1.20)--(1.21) are intrinsic Fitting conditions already contained in the
scheme-theoretic Zariski--Main package.  If they hold at every closed
collision fiber, the canonical normalization is flat.

#### Proof

Decompose `C_p` into its Artin local factors.  For a factor with maximal
ideal `n` and residue field `k`, the standard cotangent identification gives

\[
 \Omega_{C_p/k}\otimes k\simeq n/n^2.
\]

Thus its embedding dimension is at most one exactly when its cotangent
module is cyclic, by Nakayama.  A finite module is locally generated by at
most one element exactly when its first Fitting ideal is the unit ideal,
proving the equivalence of 1--3.  Over each Artin local factor, a module is
cyclic exactly when its second exterior power vanishes, proving equivalence
with 4.  Formula (1.22) is base change for Kähler differentials.  The last
assertion is Proposition 1.10.  QED

The test distinguishes the two collisions without normal forms.  For the
foundational fiber `k[epsilon]/(epsilon^3)`,

\[
 \Omega\simeq k[\epsilon]/(\epsilon^2)\,d\epsilon
\]

is cyclic.  For the reduced defect `k plus V`, `V^2=0` and `dim V=3`, so

\[
 \Omega\otimes k\simeq V
\]

requires three generators.  Intrinsic curvilinearity is therefore exactly
cotangent cyclicity, not an additional chosen root coordinate.

### Proposition 1.12 -- primitive nilradical generation is equivalent

Let `N_p` be the nilradical of the finite collision algebra `C_p`.  Then
`Spec C_p` is curvilinear if and only if `N_p` is locally principal.
Consequently, a primitive conormal element which generates `N_p` on every
closed collision fiber implies point-flatness.

#### Proof

Over the algebraically closed residue field, each Artin local factor of
`C_p` has residue field `k`.  Its maximal ideal is therefore exactly its
nilradical.  The minimal number of generators of that ideal is

\[
 \dim_k N_p/N_p^2,
\]

which is the embedding dimension.  Hence the local factor is curvilinear
exactly when its nilradical is principal.  Proposition 1.10 gives the final
assertion.  QED

For the foundational collision, `N=(epsilon)` and `N^3=0`; for the reduced
defect, `N=V`, `N^2=0`, and `N` needs three generators.  Thus a generic
primitive conormal class is not by itself enough: the closed-point
saturation theorem must say that its specialization generates the whole
nilradical.  Nilradicals, their powers, and their exact nilpotency indices
are already strata data in the formal intrinsic package.

### Proposition 1.13 -- `S_2` extension of the primitive class

Set

\[
 Q=\Omega_{B/A},
\qquad
 T=B/\operatorname{Ann}_B(Q).
\tag{1.23}
\]

Suppose globally on `Spec T` (or on open neighborhoods covering every
closed collision point) that:

1. `T` is pure of dimension two and satisfies Serre's condition `S_2`;
2. `Q` has rank one, full support, and satisfies `S_1` as a `T`-module;
3. the intrinsic primitive conormal class `tau in Q` generates `Q` at every
   point of `Spec T` of codimension at most one.

Then, on those neighborhoods,

\[
 Q=T\tau
\tag{1.24}
\]

through every closed collision point.  If the neighborhoods cover every
collision point, all fiber cotangent modules are cyclic and the canonical
normalization is flat.

#### Proof

Multiplication by `tau` gives a map

\[
 T\longrightarrow Q.
\]

Condition `S_2` implies `S_1`; together with purity this says that `T` has
no embedded associated primes.  The rank-one and codimension-zero part of
assumption 3 therefore makes this map injective.  Its cokernel `K` is
supported in codimension at least two by the rest of assumption 3.

Localize at a point of that support.  The module `T` has depth two, `Q` has
depth at least one, and `K` has finite length.  In

\[
 0\longrightarrow T\longrightarrow Q\longrightarrow K\longrightarrow0,
\]

the first local-cohomology sequence begins

\[
 0=H^0_{\mathfrak m}(Q)\longrightarrow
 K\longrightarrow H^1_{\mathfrak m}(T)=0.
\]

It forces `K=0`.  Thus (1.24) holds.  At points outside `Spec T` the
relative cotangent module is zero.  Base change therefore makes every
fiber cotangent module cyclic, provided the stated neighborhoods cover
every collision point, and Propositions 1.11 and 1.10 give flatness.  QED

The hypotheses are exactly a codimension-two Hartogs package.  The more
concrete conditions

\[
 T\text{ Cohen--Macaulay of pure dimension two},\qquad
 Q\text{ rank-one and without embedded associated points over }T
\tag{1.25}
\]

imply conditions 1--2, but are not necessary.

This theorem is compatible with the foundational triple-root collision.
There `Q` is a cyclic module over its possibly nonreduced hypersurface
support, so both are maximal Cohen--Macaulay; no reduced-support assumption
is imposed.  The remaining extraction problem has now become a precise
depth statement: prove that saturation and boundary monotonicity make the
ramification support `S_2`, make its cotangent module `S_1`, and prevent a
codimension-one zero of the primitive class.  Proposition 1.14 packages the
two closed-point failures as exact `Ext` modules.

### Proposition 1.14 -- the two-`Ext` obstruction ledger

Retain (1.23).  Assume that `T` is pure of dimension two and satisfies
`S_1`, that `Q` has rank one and full support over `T`, and that `tau`
generates `Q` in codimension at most one.  Then multiplication by `tau` is
injective and its cokernel

\[
 K=Q/T\tau
\tag{1.26}
\]

has finite length.  Moreover:

1. the failure of `T` to be `S_2` at closed points is measured by
   \[
   \operatorname{Ext}^2_A(T,A);
   \tag{1.27}
   \]
2. if (1.27) vanishes, then
   \[
   \operatorname{Ext}^3_A(Q,A)
   \simeq
   \operatorname{Ext}^3_A(K,A),
   \tag{1.28}
   \]
   so `Ext^3_A(Q,A)` is the canonical dual of the primitive-generation
   defect `K`;
3. consequently,
   \[
   \boxed{
   \operatorname{Ext}^2_A(T,A)=0,\qquad
   \operatorname{Ext}^3_A(Q,A)=0
   }
   \tag{1.29}
   \]
   imply cotangent cyclicity at every collision and finite flatness of the
   canonical normalization.

#### Proof

The `S_1` hypothesis removes embedded associated primes; as in Proposition
1.13, codimension-zero generation makes `T -> Q` injective, and
codimension-one generation makes `K` zero-dimensional.

Because `A` is a regular ring of dimension three and `T` is a pure
two-dimensional finite `A`-module, local duality identifies its only
possible closed-point `S_2` deficiency with `Ext^2_A(T,A)`.  Purity already
and `S_1` remove `Ext^3_A(T,A)`.  Thus (1.27) vanishes exactly when `T`
is `S_2`.
In that case `T` is Cohen--Macaulay of codimension one over `A`, so

\[
 \operatorname{Ext}^i_A(T,A)=0\qquad(i\ne1).
\tag{1.30}
\]

Apply `Hom_A(-,A)` to (1.26).  A finite-length module over the regular
threefold has `Ext^i_A(K,A)=0` for `i<3`.  Using (1.30), the end of the
long exact sequence gives precisely (1.28).  The right side of (1.28)
vanishes if and only if `K=0`.  Conditions (1.29) therefore give
`Q=T tau`; Propositions 1.11 and 1.10 finish the proof.  QED

This is the smallest current algebraic attack on the collision gap.  It
replaces two qualitative depth assertions by two explicit modules computed
from finite presentations already present in the intrinsic package.

### Proposition 1.15 -- double saturation computes both `Ext` obstructions

Retain the hypotheses and notation of Proposition 1.14.  Let

\[
 C=T^{[2]}
 :=\operatorname{Ext}^1_A(
      \operatorname{Ext}^1_A(T,A),A)
                                                               \tag{1.31}
\]

be the canonical `S_2` hull of `T`, with its evaluation map
`iota:T -> C`.  Then `iota` is injective, is an isomorphism at every point
of `Spec T` of codimension at most one, and `C` is pure of dimension two
and satisfies `S_2`.  Put

\[
 L=C/T,\qquad K=Q/T\tau.                                      \tag{1.32}
\]

Then `L` and `K` have finite length and

\[
 \boxed{
 \operatorname{Ext}^2_A(T,A)
 \simeq \operatorname{Ext}^3_A(L,A).
 }                                                            \tag{1.33}
\]

Consequently,

\[
 \operatorname{Ext}^2_A(T,A)=0
 \quad\Longleftrightarrow\quad
 L=0
 \quad\Longleftrightarrow\quad
 T=C.                                                         \tag{1.34}
\]

After these equivalent conditions hold,

\[
 \boxed{
 \operatorname{Ext}^3_A(Q,A)
 \simeq \operatorname{Ext}^3_A(K,A),
 }                                                            \tag{1.35}
\]

and hence the two-`Ext` certificate of Proposition 1.14 vanishes exactly
when

\[
 \boxed{L=0,\qquad K=0.}                                      \tag{1.36}
\]

In words, the closed-point problem is the conjunction of two literal
saturation statements:

1. **support saturation:** the ramification support already equals its
   `S_2` hull `C`;
2. **conormal saturation:** the primitive class already generates the
   relative cotangent module.

#### Proof

For a pure `S_1` module of codimension one over a Gorenstein ring, the
codimension-one canonical bidual in (1.31) is its `S_2`-ification: the
evaluation map is injective, is an isomorphism in codimension at most one,
and its target is `S_2`.  This can also be read from the biduality spectral
sequence for `RHom_A(RHom_A(T,A),A)`.  Thus `L` is zero-dimensional.
Since `C` is a pure `S_2` module of dimension two over the regular
threefold `A`, it is Cohen--Macaulay of codimension one.  Therefore

\[
 \operatorname{Ext}^i_A(C,A)=0\qquad(i\ne1).            \tag{1.37}
\]

Apply `Hom_A(-,A)` to

\[
 0\longrightarrow T\longrightarrow C\longrightarrow L
 \longrightarrow0.
\]

The segment in cohomological degrees two and three, together with (1.37),
gives (1.33).  A finite-length module over the regular threefold has only
one nonzero `Ext`, in degree three, and its degree-three `Ext` is its
canonical dual.  It vanishes exactly when the module does.  This proves
(1.34).

When `L=0`, the module `T=C` is Cohen--Macaulay.  Applying
`Hom_A(-,A)` to

\[
 0\longrightarrow T\mathop{\longrightarrow}^{\tau}Q
 \longrightarrow K\longrightarrow0
\]

gives (1.35), exactly as in Proposition 1.14.  Canonical duality for the
finite-length module `K` proves (1.36).  QED

Proposition 1.15 is stronger operationally than merely naming the two
`Ext` modules.  Both obstruction modules are now duals of finite quotients
already visible in the formal intrinsic package.  The first quotient is
computed by the canonical bidual of the scheme-theoretic ramification
stratum; geometrically it can be compared with the finite `S_2` model
carried by the normalized critical boundary.  The second is computed by a
single cokernel of the primitive conormal class.

The new geometric target is therefore the following **double-saturation
theorem**:

> Saturated minimal-boundary intersections force `T=C` and `Q=T tau`.

It is enough to prove this at the finite set of closed collisions.  No
global choice of a root coordinate, Tschirnhausen basis, or suspension
chart occurs in the statement.

### Proposition 1.16 -- the two saturation defects are coupled

Retain Proposition 1.15, and let `Z` be the finite closed collision locus.
Put

\[
 P=H^0_Z(Q),
\tag{1.38}
\]

the maximal submodule of the relative cotangent module supported at closed
collisions.  Then there is a canonical exact sequence

\[
 \boxed{
 0\longrightarrow P\longrightarrow K\longrightarrow L
 \longrightarrow H^1_Z(Q)\longrightarrow0.
 }
\tag{1.39}
\]

In particular,

\[
 L=0\quad\Longrightarrow\quad K\simeq P.              \tag{1.40}
\]

Consequently the two-`Ext` certificate is equivalent to

\[
 \boxed{
 L=0,\qquad P=0.
 }
\tag{1.41}
\]

Thus the closed-point theorem can be stated without an independent
generation condition after support saturation:

> The ramification support equals its canonical `S_2` hull, and the
> relative cotangent module has no closed-point torsion.

#### Proof

The exact sequence

\[
 0\longrightarrow T\longrightarrow C\longrightarrow L
 \longrightarrow0
\]

and the `S_1` condition on `T` give `H^0_Z(T)=0`; the `S_2` condition on
`C` gives `H^0_Z(C)=H^1_Z(C)=0`.  Since `L` is supported on `Z`, its
local-cohomology sequence gives a canonical isomorphism

\[
 H^1_Z(T)\simeq L.                                     \tag{1.42}
\]

Apply local cohomology to

\[
 0\longrightarrow T\mathop{\longrightarrow}^{\tau}Q
 \longrightarrow K\longrightarrow0.
\]

Here `H^0_Z(T)=0`, `H^0_Z(Q)=P`, `H^0_Z(K)=K`, and
`H^1_Z(K)=0` because `K` has finite length.  The beginning of the long
exact sequence, followed by (1.42), is exactly (1.39).

If `L=0`, (1.39) gives (1.40).  Proposition 1.15 says that vanishing of the
two `Ext` modules is equivalent to `L=K=0`; by (1.40), this is equivalent
to `L=P=0`.  QED

Proposition 1.16 changes the preferred proof strategy.  It is no longer
necessary to prove closed-point generation of `tau` directly.  Prove:

1. the canonical bidual map `T -> T^[2]` is onto;
2. `Omega_{B/A}` has no associated prime supported at a closed collision.

Then (1.39) forces `Q=T tau`, and the cubic normalization is flat.
Both statements are presentation-theoretic: the first is a canonical
bidual cokernel, while the second is the zeroth local cohomology or,
equivalently, the finite-length torsion of a Jacobian cokernel.

### Proposition 1.17 -- the finite presentation certificate

Let

\[
 F_1\mathop{\longrightarrow}^{\Psi}F_0\longrightarrow Q
 \longrightarrow0
\tag{1.43}
\]

be a finite free `A`-presentation, put `N=im(Psi)`, and let

\[
 I=\operatorname{Fitt}^A_3(B).
\tag{1.44}
\]

The support of `I` is the finite point-flatness locus.  The closed-point
cotangent torsion is

\[
 \boxed{
 H^0_I(Q)\simeq (N:_{F_0}I^\infty)/N.
 }
\tag{1.45}
\]

Therefore Certificate E vanishes if and only if the following two finite
presentation tests pass:

\[
\boxed{
\begin{aligned}
 T&\longrightarrow
 \operatorname{Ext}^1_A(
   \operatorname{Ext}^1_A(T,A),A)
 &&\text{is onto},\\
 N:_{F_0}I^\infty&=N.
\end{aligned}}
\tag{1.46}
\]

#### Proof

An element of `Q=F_0/N` belongs to `H^0_I(Q)` exactly when it is
annihilated by some power of `I`.  A representative `v in F_0` has this
property exactly when

\[
 I^n v\subset N
\]

for some `n`, which is the definition of
`v in N:_{F_0}I^\infty`.  This proves (1.45).  Proposition 1.15 identifies
surjectivity of the first map in (1.46) with `L=0`; Proposition 1.16 then
identifies the second equality with `P=0` and proves the equivalence.  QED

Formula (1.45) is the module form of the general
[support-saturation principle](../verified/SUPPORT_SATURATION_PRINCIPLE.md).
For a finite module \(M\) and an ideal \(\mathfrak a\), that principle
identifies \(H^0_{\mathfrak a}(M)=0\) simultaneously with positive
\(\mathfrak a\)-grade, avoidance of \(V(\mathfrak a)\) by every associated
prime, and saturation of any finite presentation.  Thus the cubic
cotangent obstruction and the degree-forty-two Hessian synchronization
defect are the module and algebra instances of the same theorem.  In both
problems, flatness is a sufficient shortcut, while absence of
support-torsion is the exact weaker target.

This is the current smallest machine-checkable flatness certificate.  It
requires no enumeration of collision algebras: compute one canonical
bidual and one module saturation.  The exact Singular calibration in
`scripts/verify_cubic_double_saturation.sing` separates a pure surface
summand from a single closed-point cotangent summand.

### Proposition 1.18 -- the grade-one perfect shortcut

Retain Proposition 1.14.  Suppose `T` and `Q` admit balanced finite free
presentations

\[
\begin{aligned}
 0&\longrightarrow A^r\mathop{\longrightarrow}^{\Phi}A^r
   \longrightarrow T\longrightarrow0,\\
 0&\longrightarrow A^s\mathop{\longrightarrow}^{\Psi}A^s
   \longrightarrow Q\longrightarrow0,
\end{aligned}
\tag{1.47}
\]

with `det(Phi)` and `det(Psi)` nonzero.  Then Certificate E vanishes and
the canonical cubic normalization is flat.

It is enough, more invariantly, that `T` and `Q` be perfect `A`-modules of
grade one.

#### Proof

The nonzero determinants make both displayed left maps injective over the
domain `A`.  Thus `T` and `Q` have projective dimension one.  Dualizing
(1.47) gives

\[
 \operatorname{Ext}^i_A(T,A)=
 \operatorname{Ext}^i_A(Q,A)=0
 \qquad(i\ge2).
\]

In particular the two modules in (1.29) vanish, and Proposition 1.14 gives
flatness.  The perfect grade-one formulation is equivalent locally to such
balanced resolutions.  QED

Under the one-critical-prime hypothesis, the radicals of the two
determinants in (1.47) are supported on the single branch equation.
Adjugate matrices then turn these presentations into matrix factorizations
of suitable powers of that equation.  This is the strongest convenient
target for the different/conductor attack; Proposition 1.17 remains
available when a perfect presentation cannot be extracted globally.

## 2. From a binary cubic to normalized factorization

### Proposition 2.1 -- affine Hartogs maximality

Let `X` be a normal separated integral variety and let `U=A^n` be a dense
open.  If `X\setminus U` has codimension at least two, then `U=X`.

#### Proof

Normal Hartogs extension gives

\[
 \Gamma(X,\mathcal O_X)=\Gamma(U,\mathcal O_U)
 =k[x_1,\ldots,x_n].
\]

The coordinate functions define a morphism `r:X -> U` whose restriction to
`U` is the identity.  If `j:U -> X` is the open immersion, the morphisms
`jr,id_X:X -> X` agree on dense `U`.  Their equalizer is closed because `X`
is separated, so they agree everywhere.  Thus `j` and `r` are inverse
isomorphisms.  QED

Apply this to the étale locus `V` of the canonical finite normalization.
The Keller open `U=A^3` lies in `V`.  If the intrinsic divisor ledger records
every boundary valuation and every recorded boundary divisor is critical,
then `V\setminus U` has no codimension-one component.  Proposition 2.1
forces `U=V`.

For a finite flat cubic written as a binary cubic, `V` is exactly the full
simple-marked-root locus: the relative differential at the marked root
vanishes precisely when that root is repeated.  Consequently the “no extra
simple boundary” assumption below is automatic from the exhaustive minimal
divisor ledger; it is retained in the theorem so the statement also applies
without that ledger.

### Proposition 2.1a -- the phantom-boundary factor

Let `S_F` be the nonproperness set of `F` and let

\[
 D_F=\pi\!\left(\operatorname{Supp}\Omega_{B/A}\right)
\tag{2.0a}
\]

be the branch set of the canonical finite normalization.  Then

\[
 S_F=\pi(\partial_F),\qquad D_F\subseteq S_F.
\tag{2.0b}
\]

In codimension one, the irreducible components of
`S_F\setminus D_F` are in bijection with the boundary primes at whose
generic points `pi` is unramified.

For a cubic package with one critical boundary prime, `D_F` has one
irreducible divisorial component.  If `delta_F` is its reduced equation and
`j_F` is a reduced equation for `S_F`, then

\[
 j_F=\delta_F\,u_F.
\tag{2.0c}
\]

The following are equivalent:

1. there is no second unramified boundary divisor;
2. the **phantom-boundary factor** `u_F` is a unit;
3. `S_F` and `D_F` agree in codimension one.

In particular, irreducibility of `S_F` closes the unramified-boundary gap.

#### Proof

The restriction `F=pi|_U` fails the valuative criterion for properness
exactly over the finite image of `bar X\setminus U`, giving the first
identity in (2.0b).  The Keller condition makes `pi` étale on `U`, so the
support of its relative cotangent module lies in the boundary and gives the
second inclusion.

Let `E` be a boundary prime mapping to a target divisor.  Its image belongs
to `D_F` exactly when `Omega_{B/A}` is nonzero at the generic point of
`E`, equivalently when the corresponding DVR extension is ramified.  This
proves the codimension-one correspondence.

In geometric degree three, Proposition 1.4 shows that the unique critical
prime has `(e,f)=(2,1)` and exhausts the ramification over its image.
Hence its image is the unique irreducible branch divisor.  Since
`A=k[Y]` is factorial, reduced equations `j_F,delta_F` exist, and
`delta_F` divides `j_F`.  Their quotient records exactly the remaining
divisorial images.  This proves all three equivalences.  QED

Both factors in (2.0c) are intrinsic and computable.  The branch equation is
obtained from the discriminant/different or the Fitting support of
`Omega_{B/A}`; the nonproperness equation is obtained by eliminating the
graph at infinity.  Thus the second-divisor problem is no longer a search
over compactifications: it is the single unit test `u_F in k^*`.

### Proposition 2.2 -- there is no global monogenic shortcut

Assume the canonical normalization is finite flat of degree three and its
étale locus is exactly `U=A^3`.  Then `B` is not generated by one element as
an `A`-algebra.

#### Proof

Suppose `B=A[t]`.  Since `t` has degree three over `Frac(A)`, integrality and
normality of `A` give a monic cubic `f(T) in A[T]` with

\[
 B\simeq A[T]/(f).
\]

The relative differential module is

\[
 \Omega_{B/A}\simeq B/(f'(t))\,dt,
\]

so the étale locus is the principal open `D(f'(t))`.  Its restriction to
`U` is therefore a unit.  But

\[
 \Gamma(U,\mathcal O_U)^*=k^*,
\]

hence `f'(t)=c` for some `c in k^*` in the common function field.  The
nonzero polynomial `f'(T)-c` has degree two and annihilates `t`, contradicting
that its minimal polynomial has degree three.  QED

Thus Proposition 1.10 is deliberately local.  Curvilinear generators prove
flatness point by point, but they cannot be patched into one global affine
root coordinate once the Keller open is the full étale locus.  The
binary-cubic `P^1` of root directions, and therefore its transition
markings, are genuinely necessary.

Call the intrinsic binary-cubic orbit **affine-linear of full rank** if some
Tschirnhausen basis makes `a,b,c,d` affine-linear functions on `Y` and their
affine span has dimension three.  Equivalently, one representative
`kappa_F` is an isomorphism from `Y` onto an affine hyperplane

\[
 H_\ell=\{C\in\operatorname{Sym}^3(k^2):\ell(C)=1\}
\tag{2.1}
\]

for a nonzero `ell in Sym^3(k^2)^*`.

Assume in addition that the distinguished affine open `U` is the full
simple-marked-root locus of the cubic algebra.  This condition is intrinsic:
the omitted divisor is exactly the Fitting support of relative
differentials of the finite cubic cover, with no extra unramified boundary
prime.

On that simple-root locus the cubic factors projectively as

\[
 C=LQ,\qquad \deg L=1,\quad\deg Q=2,
\tag{2.2}
\]

with `Res(L,Q)` nonzero.  The line of marked factors is a line bundle on
`U`.  Since

\[
 \operatorname{Pic}(\mathbb A^3)=0,
\qquad
 k[\mathbb A^3]^*=k^*,
\tag{2.3}
\]

one may choose `L,Q` globally and use their relative scaling to normalize

\[
 \operatorname{Res}(L,Q)=1.
\tag{2.4}
\]

The transformation law

\[
 \operatorname{Res}(\lambda L,\lambda^{-1}Q)
=\lambda\operatorname{Res}(L,Q)
\tag{2.5}
\]

shows both existence and uniqueness up to the harmless normalized scalar.
Thus `U` becomes the normalized factorization slice

\[
 X_\ell=
\{(L,Q):\operatorname{Res}(L,Q)=1,\ \ell(LQ)=1\}.
\tag{2.6}
\]

This extracts the primitive resultant marking from the finite cubic
normalization; it is not supplied as a suspension coordinate.

## 3. Affine-linear cubic classification

### Theorem 3.1 -- unique affine-linear cubic normalization

Let `F:A^3 -> A^3` be a geometric-degree-three Keller map.  Assume:

1. its canonical finite normalization has empty flatness defect (1.2);
2. its Deligne--Faddeev coefficient morphism is affine-linear of full rank;
3. its distinguished affine open is the full simple-marked-root locus.

Then `F` is polynomially left--right equivalent to the foundational map.

#### Proof

Proposition 1.3 gives the binary cubic and (2.1).  Section 2 identifies the
source map with normalized multiplication

\[
 X_\ell\longrightarrow H_\ell,\qquad(L,Q)\longmapsto LQ.
\tag{3.1}
\]

The restriction of `ell` to the twisted cubic has one of the three contact
types

\[
 (1,1,1),\qquad(2,1),\qquad(3).
\]

The cubic hyperplane-orbit theorem computes

\[
\begin{array}{c|c}
\text{contact type}&[X_\ell]\\ \hline
(1,1,1)&\mathbb L^3-\mathbb L,\\
(2,1)&\mathbb L^3,\\
(3)&\mathbb L^3-\mathbb L^2.
\end{array}
\tag{3.2}
\]

Since the distinguished source is `A^3`, only the tangent-nonosculating
type `(2,1)` is possible.  All hyperplanes of this type are in one
`PGL_2` orbit, and normalized factor rescaling restores both equations in
(2.6).  Hence all of their multiplication maps are linearly left--right
equivalent.

For the representative

\[
 L=aT+bS,\qquad Q=cT^2+dTS+eS^2,
\]

the equations are

\[
 a^2e-abd+b^2c=1,\qquad ad+bc=1.
\tag{3.3}
\]

The global coordinates

\[
 b=1+ay,\qquad
 c=1-\frac32ay+a^2z
\tag{3.4}
\]

together with the forced polynomial formulas for `d,e` identify (3.3)
with `A^3`; normalized multiplication is the foundational map after
diagonal source and target changes.  QED

### Corollary 3.2 -- positive labels are consequences in this frontend

Inside the hypotheses of Theorem 3.1, the positive cubic quotient and
conormal labels are automatic:

\[
 y=\frac{b-1}{a},\qquad
 z=\frac{c-1+\frac32ay}{a^2}.
\tag{3.5}
\]

They are the first two saturated coefficients in the unique
tangent-nonosculating factorization slice.  Thus the separate positive-chart
straightening assumptions in the suspension approach are unnecessary in
the affine-linear finite-normalization frontend.

### Corollary 3.3 -- curvilinear normalization gateway

Let `F:A^3 -> A^3` have geometric degree three.  Assume:

1. every closed fiber of its canonical finite normalization is
   curvilinear, equivalently cotangent-cyclic in the sense of Proposition
   1.11; it is enough that the ramification support is `S_2`, the cotangent
   module is `S_1`, and the primitive class generates in codimension one as
   in Proposition 1.13.  Equivalently, the two intrinsic obstruction
   modules in Proposition 1.14 may be required to vanish;
2. the intrinsic binary-cubic orbit, once extracted, has an affine-linear
   full-rank representative;
3. the canonical boundary has no unramified prime over a second target
   divisor.

Then `F` is polynomially left--right equivalent to the foundational map.

#### Proof

Proposition 1.10 gives finite flatness and hence the binary-cubic package of
Proposition 1.3.  The third assumption and affine Hartogs maximality identify
the Keller open with the full simple-marked-root locus.  Theorem 3.1 then
applies.  QED

This is strictly more geometric than assuming an empty Fitting defect:
curvilinearity is read directly from the finite collision algebras in the
scheme-intersection package.  It is also calibrated to keep the
foundational triple-root fiber.

## 4. What remains open

Theorem 3.1 does not construct the finite-normalization witness of
`MINIMAL_BOUNDARY_CLASSIFICATION.md` or prove its eight predicates from an
unmarked boundary-minimal map.  Within this finite-normalization frontend,
the cubic extraction problem is split into exact tests:

1. **point-flatness:** prove `Fitt_3^A(B)=A`, excluding a defect supported at
   finitely many target points.  Equivalently, prove cubic fiber-minimality
   `lambda(p)<=3` there (Proposition 1.5); Proposition 1.2 gives the
   stronger sufficient route through a Cartier--Cohen--Macaulay canonical
   boundary.  For a reduced minimal defect, Propositions 1.7--1.8 reduce
   this further to excluding the square-zero collision of the ramified and
   affine sheets while retaining the foundational curvilinear triple-root
   collision.  Most efficiently, Proposition 1.10 shows that intrinsic
   curvilinearity of all collision fibers eliminates every point defect at
   once;
2. **coefficient linearity:** prove that the intrinsic binary-cubic orbit
   has a full-rank affine hyperplane representative modulo polynomial
   Tschirnhausen gauge;
3. **unramified boundary:** exclude a recorded height-one boundary prime
   inside the étale locus.  By Proposition 1.4 it would have to map to a
   second target nonproperness divisor, distinct from the critical
   discriminant.

Proposition 2.1 shows that there is no additional codimension-two or
codimension-three version of the third obstruction: once unramified boundary
divisors are excluded, the full simple-root locus equals `U`.  If
“divisor-minimal” is strengthened to say that every recorded boundary prime
is critical, the third test is automatic.  If all three tests hold, the
one-place/two-place suspension distinction is no longer needed to prove
cubic uniqueness: Theorem 3.1 goes directly from the canonical normalization
to the foundational map.

Proposition 2.2 rules out a tempting over-simplification after flatness:
the locally monogenic collision charts cannot be replaced by a single
global monic cubic coordinate.  Their `P^1` transition is precisely where
the remaining coefficient-gauge and intrinsic-marking information lives.

The second item is now the main geometric obstruction.  A general finite
flat cubic algebra over `A^3` gives a nonlinear morphism
`A^3 -> A^4`; neither flatness nor the discriminant divisor alone makes its
image a hyperplane.

Here affine-linearity is always understood modulo the polynomial
`GL_2(A)` Tschirnhausen gauge.  The
[cubic gauge-straightening theorem](CUBIC_GAUGE_STRAIGHTENING.md) proves that
every invariant-parameter nonlinear slice

\[
 C_1=q-3C_0h,\qquad q\ne0,\quad h\in\ker D_+,
\]

where `D_+` is the translation locally nilpotent derivation, is
gauge-equivalent by explicit source and target polynomial automorphisms to
the tangent slice `C_1=q`; the opposite shear has the symmetric theorem.
Moreover, a variable-time upper or lower shear is a polynomial automorphism
only when its time is invariant.  Thus these theorems exhaust the
single-unipotent-shear automorphisms, rather than merely giving a broad
family of examples.  The open case is transverse to finite compositions of
these polynomial unipotent gauge orbits or fails to admit such a
Tschirnhausen reduction.

## 5. Reproduction and external inputs

Run

```bash
.venv/bin/python scripts/verify_cubic_normalization_frontend.py
Singular -q scripts/verify_cubic_double_saturation.sing
```

The checker verifies the universal cubic-algebra multiplication table,
trace-zero splitting, trace discriminant, the codimension-three reflexive
module warning with its excess-length-four special fiber, the canonical
`S_2`-hull calibration, the coupling of conormal failure with point torsion,
and the exact tangent-hyperplane quotient coordinates.  The Singular
checker verifies the module-saturation formula (1.45) on a pure surface
module with one closed-point cotangent summand.

External structural inputs:

- the cubic-algebra/binary-cubic equivalence over an arbitrary base:
  Wood, *Parametrizing quartic algebras over an arbitrary base*,
  [Theorem 2.1](https://msp.org/ant/2011/5-8/ant-v5-n8-p05-p.pdf);
- nonflat triple covers with prescribed reflexive trace-free sheaf,
  including the generic reduced-irreducible divisor construction in
  Theorem 7.6 and the generalized cubic tensor used in Proposition 1.8a:
  Tan, *Triple covers on smooth algebraic varieties*,
  [Theorem 7.6, pp. 162--163](https://math.ecnu.edu.cn/~sltan/Triple-Cover.pdf);
- the canonical bidual characterization of the `S_2` condition for modules:
  Hashimoto, *Canonical and n-canonical modules of a Noetherian algebra*,
  [Theorem 8.4](https://www.math.okayama-u.ac.jp/~hashimoto/paper/ncan7.pdf);
- projective modules over a polynomial ring are free:
  Suslin, *Projective modules over a polynomial ring are free*,
  [MathNet](https://www.mathnet.ru/eng/dan40545);
- maximal Cohen--Macaulay modules over regular local rings are free:
  [Stacks Project, Tag 00NT](https://stacks.math.columbia.edu/tag/00NT);
- complements of affine opens and the normal Hartogs argument used in
  Proposition 2.1:
  [Stacks Project, Tag 0BCQ](https://stacks.math.columbia.edu/tag/0BCQ).

The normalized hyperplane-orbit theorem and its motivic exclusions are
internal:
[`../verified/FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md`](../verified/FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md).
