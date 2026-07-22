# The normalized quadratic-cubic factorization slice

This note studies the first consecutive-degree case beyond the foundational
linear-quadratic slice.  The conclusion is negative but informative: the
natural normalized source is a smooth factorial affine fivefold with trivial
canonical class, yet it is not affine five-space.

Work over a field `k` of characteristic zero.  Write

\[
 A=a_0T^2+a_1TS+a_2S^2,\qquad
 B=b_0T^3+b_1T^2S+b_2TS^2+b_3S^3
\]

and use the tangent coefficient

\[
 m=[AB]_{T^4S}=a_0b_1+a_1b_0.
\]

The normalized affine slice is

\[
 X_{2,3}=\{m=1,\ \operatorname{Res}(A,B)=1\}\subset\mathbb A^7. \tag{1}
\]

The unequal-degree theorem makes multiplication etale on (1).  Consecutive
degrees make the normalization unique, so projectivization identifies (1)
with

\[
 U_{2,3}=(\mathbb P^2\times\mathbb P^3)\setminus(R\cup E), \tag{2}
\]

where `R` is the resultant divisor and `E={m=0}`.  A squarefree binary
quintic has `binomial(5,2)=10` quadratic factors, so multiplication is
generically ten-to-one.

All conclusions below concern this natural tangent coefficient.  Other
hyperplanes in the dual binary-quintic space need not have the same
intersection geometry.

## 1. Units and Picard group

The ambient Picard group is

\[
 \operatorname{Pic}(\mathbb P^2\times\mathbb P^3)=\mathbb Z^2.
\]

The binary resultant is irreducible of bidegree `(3,2)`, while the rank-two
bilinear form `m` is irreducible of bidegree `(1,1)`.  Divisor localization
therefore uses the matrix

\[
 \begin{pmatrix}3&1\\2&1\end{pmatrix},\qquad \det=1.  \tag{3}
\]

It is an isomorphism of lattices.  Hence

\[
 \boxed{\mathcal O(U_{2,3})^*/k^*=0,\qquad
 \operatorname{Pic}(U_{2,3})=0.}                      \tag{4}
\]

Thus the first obstructions agree with `A^5`; they do not decide the
question.

## 2. Grothendieck class

Let `L=[A^1]` and `[P^j]=1+L+...+L^j`.  Stratify the resultant divisor by the
degree of `gcd(A,B)`.

- For gcd degree one, choose the common linear form and a coprime residual
  `(1,2)` pair.  The latter has class
  \[
   [P^1][P^2]-[P^1]^2=L^2+L^3.
  \]
- For gcd degree two, choose the quadratic gcd and the remaining linear
  factor, giving `P^2 x P^1`.

Consequently

\[
 [\{\operatorname{Res}\ne0\}]=L^4+L^5.               \tag{5}
\]

It remains to impose `m!=0`.  The coprime part of `E` has class `L^4+L^3`.
Here is a direct stratification.

If `a_0!=0`, make `A` monic.  Polynomial division identifies the hyperplane
of possible `B` with triples `[b_0:r_1:r_0]`, where
`r_1T+r_0` is the remainder modulo `A`.  Coprimality says this nonzero linear
remainder does not vanish at a root of `A`.  Projecting to its root gives a
space of class `L^3`, and the free `b_0` contributes another `L`; this chart
has class `L^4`.

If `a_0=0`, the equation `m=0` and coprimality force `A=S^2` and `b_0!=0`.
The remaining cubic coefficients form `A^3`, contributing `L^3`.  Therefore

\[
 [E\cap\{\operatorname{Res}\ne0\}]=L^4+L^3.          \tag{6}
\]

Subtracting (6) from (5) gives

\[
 \boxed{[U_{2,3}]=L^5-L^3=L^2[\mathrm{SL}_2].}        \tag{7}
\]

The Hodge--Deligne realization sends this to

\[
 E_c(U_{2,3};u,v)=(uv)^5-(uv)^3,\qquad
 E(U_{2,3};u,v)=1-(uv)^2.                              \tag{8}
\]

In particular (7)--(8) prove

\[
 \boxed{U_{2,3}\not\simeq\mathbb A^5}.               \tag{9}
\]

The suggestive equality with `L^2[SL_2]` is motivic only; no product
isomorphism `U_{2,3}=A^2 x SL_2` is asserted.

## 3. Finite-field counts

The same strata are defined over every finite field and give

\[
 \boxed{\#U_{2,3}(\mathbb F_q)=q^5-q^3.}              \tag{10}
\]

This differs from `#A^5(F_q)=q^5` for every prime power.  It is an independent
realization of the motivic obstruction and avoids relying on cancellation in
the Grothendieck ring.

## 4. Low-degree cohomology: what is determined

Over `C`, the standard meridian sequence for the complement of two prime
divisors begins

\[
 0=H^1(\mathbb P^2\mathbin\times\mathbb P^3,\mathbb Q)
 \longrightarrow H^1(U_{2,3},\mathbb Q)
 \longrightarrow \mathbb Q(-1)^2
 \xrightarrow{\delta}H^2(\mathbb P^2\mathbin\times\mathbb P^3,\mathbb Q).
\]

In the two hyperplane bases, `delta` is exactly the divisor-class matrix
(3).  It is invertible.  The open is connected, so this proves

\[
 \boxed{H^0(U_{2,3},\mathbb Q)=\mathbb Q,\qquad
 H^1(U_{2,3},\mathbb Q)=0.}                          \tag{11}
\]

Equation (8) computes the complete **virtual** mixed Hodge polynomial.  Its
first nonconstant term is `-(uv)^2`, as for `A^2 x SL_2`.  It is tempting to
identify this with `H^3=Q(-2)` and all other positive-degree rational
cohomology zero.  That conclusion does not follow from the Grothendieck class:
different cohomological degrees can cancel in the Hodge--Deligne polynomial.
A full normal-crossings or simplicial resolution of `R union E`, including
the degree-two gcd conductor, is still required to separate `H^2` and `H^3`.
Accordingly the proved low-degree statement here is (11), not the stronger
`SL_2`-shaped cohomology guess.

For context, Farb--Wolfson compute the topology and arithmetic of the
`resultant=1` space for **two monic forms of the same degree** in
[`Topology and arithmetic of resultants, II`](https://arxiv.org/abs/1507.01283).
Their purity theorem is useful guidance, but it does not directly cover this
unequal-degree projective complement with the additional tangent divisor.
No cohomology assertion above is imported from that paper.

## 5. Factoriality and canonical class

The slice is smooth because it is the base change of the etale
coefficient-resultant map.  Hence it is regular and locally factorial.  For a
regular affine variety, `Cl=Pic`; (4) therefore proves

\[
 \boxed{k[X_{2,3}]\text{ is a UFD}.}                  \tag{12}
\]

The ambient canonical class in the projective model is `(-3,-4)`.  In the
boundary lattice,

\[
 (-3,-4)=(3,2)-6(1,1).                                \tag{13}
\]

Both boundary classes vanish after restriction to `U`, so

\[
 \boxed{K_{U_{2,3}}=0.}                               \tag{14}
\]

Equivalently, (1) is a smooth codimension-two complete intersection in
affine seven-space, and the Poincaré-residue form trivializes its canonical
bundle.

## 6. Test of the `A^2 x SL_2` guess

The two positive formulations in the product question are not independent if
``SL_2-bundle'' means a principal torsor.  An `SL_2`-torsor on `A^2` is the
oriented-frame torsor of a rank-two vector bundle with a determinant
trivialization.  Quillen--Suslin makes the vector bundle free, and the
determinant trivialization can then be adjusted by a constant change of
basis.  Thus

\[
 H^1_{\mathrm{fppf}}(\mathbb A^2,\mathrm{SL}_2)=1.       \tag{15}
\]

Consequently every Zariski, etale, or fppf locally trivial principal
`SL_2`-bundle over `A^2` is the product.  This does not cover a fiber bundle
whose structure group is the larger group `Aut(SL_2)`.

There is a particularly natural determinant-one morphism in the opposite
direction.  On the affine model (1), the equation `m=1` says

\[
 \det\begin{pmatrix}a_0&a_1\\-b_0&b_1\end{pmatrix}=1.
\]

It therefore defines

\[
 p:X_{2,3}\longrightarrow\mathrm{SL}_2.                \tag{16}
\]

This morphism is **not** an `A^2`-bundle.  Indeed, over the identity matrix we
have `a_0=b_1=1` and `a_1=b_0=0`.  The resultant equation becomes

\[
 b_3^2-2a_2b_3+a_2b_2^2+a_2^2=1.
\]

With `x=b_3-a_2`, `y=b_2`, and `z=a_2`, the fiber is the smooth affine
surface

\[
 S=\{x^2+zy^2=1\}\subset\mathbb A^3.                  \tag{17}
\]

Its `y!=0` part is `A^1 x G_m`, while its `y=0` part is the disjoint union of
two affine lines.  Hence

\[
 [S]=L(L-1)+2L=L^2+L\ne L^2.                          \tag{18}
\]

Thus even the determinant-one coordinates which visibly account for the
`SL_2`-shaped equation do not split off an `SL_2` factor.  This also explains
why Euclidean division does not immediately produce the desired product.
On `a_0!=0` one can write

\[
 B=(q_1T+q_0S)A+S^2(r_1T+r_0S),
\]

and the resultant is controlled by the norm of the linear remainder.  The
chart does not extend regularly across `a_0=0`.  More intrinsically, replacing
`B` by `B+(\ell_0T+\ell_1S)A` changes `m` by

\[
 a_0(2a_1\ell_0+a_0\ell_1).
\]

The space of such additions preserving `m` has rank one when `a_0!=0` and
rank two when `a_0=0`.  The obvious additive quotient therefore jumps rank
on a nonempty divisor and cannot supply a locally trivial `A^2` factor.

Equations (17)--(18) and this rank jump are exact obstructions to the two
natural splitting constructions.  They are not an abstract non-isomorphism
theorem: an isomorphism `X_(2,3) = A^2 x SL_2` could mix all seven coefficient
functions and need not respect (16).  In particular the computations do not
yet separate `H^2` and `H^3`.  The next invariant must therefore go beyond
the natural coefficient maps.

## 7. The complex fundamental group

The fundamental group can in fact be computed without resolving the full
boundary.  Put

\[
 W=\{(A,B)\in V_2\mathbin\times V_3:
             \operatorname{Res}(A,B)m(A,B)\ne0\}.
\]

The complement has two irreducible hypersurface components.  The resultant
is smooth at a pair with exactly one common factor and, at a generic pair with
a squarefree quadratic gcd, has an ordinary two-branch crossing.  Its worse
singular locus has codimension at least three in the coefficient space.  The
rank-two bilinear hypersurface `m=0` has singular locus of codimension four,
and it meets the smooth resultant locus transversely at a general point.
Consequently a general affine two-plane cuts the union
`{Res=0} union {m=0}` in a reduced nodal plane curve with two irreducible
components.

The Zariski--Lefschetz theorem for hypersurface complements identifies the
fundamental group of `W` with that of this general plane section.  The
Fulton--Deligne theorem makes the complement of a nodal plane curve abelian.
Its first homology is freely generated by the meridians of the two irreducible
components, so

\[
 \pi_1(W)=\mathbb Z^2.                                \tag{19}
\]

This is the two-component version of Choudary's computation
`pi_1(V \setminus {Res=0})=Z`; Shimada proves the corresponding generalized
resultant theorem.  See
[`On the resultant hypersurface`](https://doi.org/10.2140/pjm.1990.142.259)
and
[`The fundamental group of the complement of a resultant hypersurface`](https://msp.org/pjm/2003/210-2/pjm-v210-n2-p09-s.pdf).

It remains to remove the normalization torus.  There is an explicit algebraic
product

\[
 W\simeq X_{2,3}\mathbin\times\mathbb G_m^2.          \tag{20}
\]

Indeed, if `r=Res(A,B)` and `s=m(A,B)`, rescale by

\[
 A'=\frac{s^2}{r}A,\qquad B'=\frac{r}{s^3}B.
\]

Then `m(A',B')=Res(A',B')=1`; retaining `(r,s)` gives (20), and the inverse
uses the inverse scalars.  Under (20), the two meridians in (19) map to the
standard generators of the normalization torus.  Hence

\[
 \boxed{\pi_1(U_{2,3}(\mathbb C))=0.}                 \tag{21}
\]

Thus `U_(2,3)` and `A^2 x SL_2` agree also at the level of the complex
fundamental group.  The next subsection extracts the two previously missing
cohomology groups from the same small boundary model.

## 8. The groups `H^2` and `H^3` without a full resolution

Let `Y=P^2 x P^3`, `D=R union E`, and `Z=R intersect E`.  Since `Y` is a
smooth compact complex fivefold, the cohomology sequence with supports and
Poincare--Lefschetz duality give

\[
 H_8(D)\longrightarrow H^2(Y)\longrightarrow H^2(U)
 \longrightarrow H_7(D)\longrightarrow H^3(Y)=0,      \tag{22}
\]

and

\[
 0\longrightarrow H^3(U)\longrightarrow H_6(D)
 \longrightarrow H^4(Y).                              \tag{23}
\]

All groups in this subsection have rational coefficients.

Only the top two homology groups of the boundary are needed.  They can be
computed from the incidence normalization.  The normalization of `R` is

\[
 I=\{(A,B,x):A(x)=B(x)=0\}\longrightarrow R.
\]

It is a `P^1 x P^2`-bundle over `P^1` and has only even cohomology.  The
conductor is the gcd-two locus `P^2 x P^1`; its inverse image marks one of the
two roots of the quadratic gcd and is connected.  The normalization exact
sequence therefore gives

\[
 H_7(R)=0,\qquad \dim H_6(R)=3.                       \tag{24}
\]

Similarly, resolving the rational maps
`[a_0:a_1]` and `[b_1:b_0]` turns `E` into a
`P^1 x P^2`-bundle over `P^1`.  The exceptional locus is too small to affect
degree six or seven, so

\[
 H_7(E)=0,\qquad \dim H_6(E)=3.                       \tag{25}
\]

The intersection `Z` has exactly two irreducible threefold components.  To
see them, write a common linear factor as `L=vT-uS` and put

\[
 A=L(\alpha T+\beta S),\qquad
 B=L(cT^2+dTS+eS^2).
\]

Then

\[
 m=v\bigl(v(\alpha d+\beta c)-2u\alpha c\bigr).       \tag{26}
\]

The factor `v=0` gives

\[
 Z_\infty=\{a_0=b_0=0\}\simeq\mathbb P^1\times\mathbb P^2,
\]

and the second factor in (26) has irreducible closure `Z_0`.  The incidence
normalizations of `Z_0`, its self-conductor, and
`Z_infinity intersect Z_0` again have connected rational components and only
even top cohomology.  Concretely, the normalization of `Z_0` projects to
`P^1 x P^1`; its fiber is `P^1` except at one point, where it is `P^2`, and
the local equation there is the ordinary node `sc+td=0`.  A small resolution
has an iterated projective-bundle paving.  The self-conductor has connected
inverse image.  Moreover `Z_infinity intersect Z_0` is the union of the two
rational surfaces `alpha=0` and `c=0`; their top classes are already
independent in `H_4(Z_infinity)`.  The normalization and Mayer--Vietoris
sequences therefore give

\[
 H_5(Z)=0,\qquad H_6(Z)=\mathbb Q[Z_\infty]
                         \mathbin\oplus\mathbb Q[Z_0]. \tag{27}
\]

These two top classes inject into `H_6(R) direct-sum H_6(E)`.  Indeed, in the
basis `(h_1^2,h_1h_2,h_2^2)` of `H^4(Y)`, their pushforward classes are

\[
 [Z_\infty]=(0,1,0),\qquad [Z_0]=(3,4,2),             \tag{28}
\]

because their sum is
`[R][E]=(3h_1+2h_2)(h_1+h_2)`.  They are independent.
Mayer--Vietoris now yields

\[
 H_7(D)=0,\qquad \dim H_6(D)=3+3-2=4.                \tag{29}
\]

The first map in (22) is the boundary-class matrix (3), hence is an
isomorphism.  This proves `H^2(U)=0`.  The last map in (23) is surjective:
for example the three boundary cycles `E h_1`, `E h_2`, and `R h_1` have
classes

\[
 (1,1,0),\qquad(0,1,1),\qquad(3,2,0),
\]

which span `H^4(Y)`.  Equations (23) and (29) therefore leave a
one-dimensional kernel.  Its generator is the difference of the two edge
classes in the boundary dual graph and has the double-Gysin Tate twist.
Consequently

\[
 \boxed{H^2(U_{2,3},\mathbb Q)=0,\qquad
 H^3(U_{2,3},\mathbb Q)=\mathbb Q(-2).}               \tag{30}
\]

Thus the unresolved groups in the original audit have exactly the
`A^2 x SL_2` values.  This calculation uses only normalizations, conductors,
and top boundary homology; it does not require a full normal-crossings
resolution.  It still does not exclude cancelling classes in degrees four
and five, nor does it produce an algebraic product presentation.

## 9. Consequence

The natural `(2,3)` multiplication map is an etale, generically degree-ten
map

\[
 U_{2,3}\longrightarrow\mathbb A^5,
\]

but its source is not `A^5`.  It therefore does **not** produce a new `JC(5)`
counterexample.  It is still a genuinely new nonproper etale cover model and
its motive, fundamental group, and cohomology through degree three all have
the `A^2 x SL_2` shape.  The remaining questions are whether degrees four and
five vanish and whether the failure of the natural splittings in Section 6
can be upgraded to an abstract product obstruction.

The exact class calculation and direct finite-field enumeration are checked
by
[`verify_quadratic_cubic_factorization_invariants.py`](../scripts/verify_quadratic_cubic_factorization_invariants.py).
The determinant fiber and rank-jump calculations in Section 6 are checked by
[`verify_quadratic_cubic_product_test.py`](../scripts/verify_quadratic_cubic_product_test.py).
The normalization product and a transverse point of the two affine boundary
components used in Section 7 are checked by
[`verify_quadratic_cubic_fundamental_group_inputs.py`](../scripts/verify_quadratic_cubic_fundamental_group_inputs.py).

The determinant-one calculation, normalization uniqueness, Picard
vanishing, and canonical-class identity are consolidated for arbitrary
factor degrees in the
[relative-scaling/boundary-lattice theorem](RELATIVE_SCALING_BOUNDARY_LATTICE.md).
