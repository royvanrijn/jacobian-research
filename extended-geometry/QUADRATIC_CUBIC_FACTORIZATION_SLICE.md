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
The virtual polynomial alone therefore proves only (11), not the stronger
`SL_2`-shaped cohomology guess.  Sections 7--8 below supply the missing
geometric input: a generic-plane argument and the top homology of the
incidence normalization separate `H^2` and `H^3` without constructing a full
normal-crossings resolution.

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

### Integral Chow groups and algebraic `K_0`

The same projective boundary gives stronger integral information.  Put
`h_1=c_1(O(1,0))` and `h_2=c_1(O(0,1))`.  In the basis
`(h_1^2,h_1h_2,h_2^2)` of `CH^2(P^2 x P^3)`, the three cycles supported on
the boundary

\[
 Eh_1=(1,1,0),\qquad Eh_2=(0,1,1),\qquad Rh_1=(3,2,0)
\]

form a unimodular matrix.  The localization sequence therefore gives

\[
 CH^2(U_{2,3})=0                                      \tag{14a}
\]

integrally.  In fact all positive-codimension Chow groups vanish.  In
codimension one this is the boundary-lattice calculation (3).  In
codimensions three through five, multiplication by `E=h_1+h_2` supplies
successive triangular unimodular spanning sets in the standard monomial
bases.  Thus

\[
 \boxed{CH^i(U_{2,3})=0\quad\text{for every }i>0.}    \tag{14b}
\]

This includes `CH_2(U)=CH^3(U)=0` if Chow groups are indexed by dimension.

The localization map `K_0(P^2 x P^3) -> K_0(U)` is surjective because both
varieties are smooth.  The nonvanishing boundary equations trivialize
`O(3,2)` and `O(1,1)` on `U`.  Their classes form a unimodular basis of
`Z^2`, so `O(1,0)` and `O(0,1)` are also trivial on `U`.  The projective
bundle presentation of `K_0(P^2 x P^3)` now shows that its entire image is
generated by `[O_U]`; the rank map prevents any further relation.  Hence

\[
 \boxed{K_0(U_{2,3})=\mathbb Z.}                     \tag{14c}
\]

These integral invariants again agree with `A^2 x SL_2`.

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
resolution.  The next section removes the remaining possibility of
cancelling classes in degrees four and five and determines the integral
cohomology.

## 9. Integral cohomology and the homotopy type

The two elementary coefficient charts give a shorter calculation of the
whole cohomology.  On the open `X^o={a_0!=0}`, make `A` monic and use `m=1`:

\[
 A=T^2+uT+v,\qquad
 B=qT^3+(1-uq)T^2+b_2T+b_3.
\]

Euclidean division has the form

\[
 B=(qT+1-2uq)A+yT+x,
\]

and

\[
 \Delta=\operatorname{Res}(A,B)=x^2-uxy+vy^2.       \tag{31}
\]

The quotient coordinate `q` is free, so

\[
 X^o\simeq\mathbb A^1\times V,\qquad
 V=\operatorname{Spec}
   \mathbb C[u,v,x,y,\Delta^{-1}].                   \tag{32}
\]

Stratify `V` by `y`.  When `y!=0`, replace `v` by `Delta`; when `y=0`,
equation (31) says `x!=0`.  Therefore

\[
 V_{y\ne0}\simeq\mathbb A^2\times\mathbb G_m^2,
 \qquad
 V_{y=0}\simeq\mathbb A^2\times\mathbb G_m.         \tag{33}
\]

The integral Gysin sequence for the smooth divisor `{y=0}` is determined by
two residues.  In logarithmic notation they are

\[
 \operatorname{res}(d\log y)=1,
 \qquad
 \operatorname{res}(d\log y\wedge d\log\Delta)
   =d\log(x^2)=2d\log x.                              \tag{34}
\]

It follows that

\[
 H^i(V,\mathbb Z)=
 \begin{cases}
  \mathbb Z,&i=0,1,\\
  \mathbb Z/2,&i=3,\\
  0,&\text{otherwise}.
 \end{cases}                                         \tag{35}
\]

The order-two class is the integral trace of the quadratic-resultant
monodromy.  It is present on the Euclidean chart but disappears after the
missing divisor is restored.

Indeed, the divisor `D_0={a_0=0}` is smooth.  Normalizing `a_1=b_0=1` gives

\[
 A=T+v,\qquad B=T^3+qT^2+yT+x,
\]

and coprimality is the nonvanishing of

\[
 \rho=-x+vy-qv^2+v^3.
\]

Thus

\[
 D_0\simeq\mathbb A^3\times\mathbb G_m.             \tag{36}
\]

On the overlap, put `r=a_0/a_1`.  Direct substitution gives

\[
 \Delta=\frac{R_1}{r},\qquad R_1|_{r=0}=\rho,
 \qquad
 y_{\rm rem}=\frac{1+r(\cdots)}{r}.                  \tag{37}
\]

The generator `d log Delta` of `H^1(X^o,Z)` consequently has residue `-1`
along `D_0`.  Modulo two, the residue of the torsion class represented on
`y_rem!=0` by
`d log(y_rem) wedge d log(Delta)` is `d log(rho)`, the nonzero generator of
`H^1(D_0,Z/2)`.  The integral Gysin sequence for
`D_0 subset X` therefore contains the nonsplit extension

\[
 0\longrightarrow\mathbb Z
  \longrightarrow H^3(X,\mathbb Z)
  \longrightarrow\mathbb Z/2\longrightarrow0.       \tag{38}
\]

The middle group is `Z`, with the first map multiplication by two.  All
other positive groups vanish:

\[
 \boxed{
 H^i(X_{2,3}(\mathbb C),\mathbb Z)=
 \begin{cases}
  \mathbb Z,&i=0,3,\\
  0,&\text{otherwise}.
 \end{cases}}                                        \tag{39}
\]

In particular `H^4=H^5=0`; there are no hidden cancellations in (8).
Since a complex affine variety has finite CW type, the universal coefficient
theorem gives the same statement for integral homology.  Together with the
simple connectivity proved in Section 7,
Hurewicz supplies a map from `S^3` representing its generator, and the
homological Whitehead theorem gives

\[
 \boxed{X_{2,3}(\mathbb C)\simeq S^3}                \tag{40}
\]

as ordinary topological spaces up to homotopy.  This is again the homotopy
type of `A^2 x SL_2`.

## 10. The Euclidean chart as an affine modification

The chart calculation also makes the modification analogy exact.  Put

\[
 w=xu-vy,\qquad \Delta=x^2-yw.
\]

Then

\[
 k[X^o]=
 k[q,u,x,y,w,v,\Delta^{-1}]/(yv-xu+w).               \tag{41}
\]

Let

\[
 A_0=k[q,u,\alpha,\beta,\gamma,\delta]/
       (\alpha\delta-\beta\gamma-1),
\]

the coordinate ring of `A^2 x SL_2`, and set

\[
 f=\alpha\beta\delta,\qquad
 I=(f,\ \beta\delta^2,\ \alpha^2\beta,\
        \ \alpha^2(u\delta-\gamma)).                 \tag{42}
\]

The birational substitution

\[
 \alpha=x,\qquad \beta=y,\qquad
 \gamma=\frac{w}{\Delta},\qquad
 \delta=\frac{x}{\Delta}
\]

gives the exact affine-modification presentation

\[
 \boxed{k[X^o]=A_0[I/f].}                            \tag{43}
\]

The reduced center has two disjoint components

\[
 C_1=\{\alpha=\delta=0,\ \beta\gamma=-1\},\qquad
 C_2=\{\beta=0,\ \alpha\delta=1,\ \gamma=u\delta\},
\]

each isomorphic to `G_m x A^2`.  The exceptional divisor has the two
components

\[
 E_1=\{x=0\}\simeq\mathbb G_m^2\times\mathbb A^2,
 \qquad
 E_2=\{y=0\}\simeq\mathbb G_m\times\mathbb A^3.     \tag{44}
\]

This `SL_2` is assembled from the Euclidean remainder and is not the natural
determinant-one morphism (16).  Hence the non-`A^2` fiber (17) does not
conflict with (43).

The displayed modification does not extend directly across `a_0=0`.
Equation (37) shows that its `alpha` and `beta` coordinates have poles there,
although `gamma` and `delta` stay finite.  Thus a global affine modification
of `A^2 x SL_2`, if one exists, requires a different birational
identification.  Likewise none of the Chow, `K_0`, or topological
calculations proves or disproves an abstract product isomorphism.  The
Makar--Limanov invariant and flexibility remain the first potentially
separating additive invariants; Section 11 computes the ordinary Derksen
algebra and shows that it does not separate the varieties.

## 11. Additive actions and a global cylinder

For this section abbreviate

\[
 (a,c,z,p,d,y,x)=(a_0,a_1,a_2,b_0,b_1,b_2,b_3).
\]

The Euclidean addition `B -> B+t(aT-2c)A` gives the homogeneous locally
nilpotent derivation

\[
 D_{22}=(0,0,0,a^2,-ac,az-2c^2,-2cz).                \tag{45}
\]

It is not the only global additive action.  There is a second homogeneous
locally nilpotent derivation

\[
 \begin{aligned}
 D_{10}(a)&=D_{10}(p)=0,&D_{10}(c)&=-2a^2,
 &D_{10}(z)&=-ac,\\
 D_{10}(d)&=2ap,&D_{10}(y)&=-2ad+5cp,
 &D_{10}(x)&=-ay+5pz.                                \tag{46}
 \end{aligned}
\]

Both derivations annihilate `m` and the resultant as ambient polynomials,
and they commute.  Their local nilpotence is triangular: `a,p` are fixed by
`D_(10)`, after which its iterates on the remaining five generators vanish
in orders at most four.  In particular (46) moves the quadratic factor and
rules out any argument claiming that every `Ga`-action must preserve `A`.

More importantly, (46) has the global slice

\[
 s={14apx-14cd^2-7cpy-35dpz+9y\over10},
 \qquad D_{10}(s)=1\quad\hbox{in }k[X_{2,3}].          \tag{47}
\]

Indeed the ambient error in (47) is

\[
 D_{10}(s)-1={ (28ad-35cp+10)(m-1)\over10}.
\]

The slice theorem therefore gives a genuine global product

\[
 \boxed{X_{2,3}\simeq \mathbb A^1_s\times Y,
 \qquad k[Y]=\ker D_{10}.}                            \tag{48}
\]

Here `Y={s=0} subset X_(2,3)` is a smooth affine fourfold.  Thus the original
question has already undergone one cancellation: it is whether this
particular cylinder is isomorphic to `A^1 x (A^1 x SL_2)`.

There is also a useful coarse simultaneous quotient.  Set

\[
 \begin{aligned}
 U&=4az-c^2,\\
 K&=4a^2y-4c+8c^2p-4apz,\\
 J&=8a^3x-4a^2cy+4c^2-8c^3p+20acpz-8az.
 \end{aligned}                                      \tag{49}
\]

The four functions `a,U,K,J` are killed by both (45) and (46), and on `X`
they satisfy

\[
 \boxed{J^2+UK^2=64a^3.}                             \tag{50}
\]

Thus the two visible additive actions map to the trinomial threefold

\[
 Q=\{J^2+UK^2=64a^3\}\subset\mathbb A^4.             \tag{51}
\]

This makes the rank jump precise.  After inverting `a`, the two actions have
successive slices `s` and `p/a^2`, and direct reconstruction gives

\[
 k[X]_{a}\simeq k[Q]_{a}[p,s].                       \tag{52}
\]

For example

\[
 c={7Jp+9K-12Up\over20}-2a^2s,                       \tag{53}
\]

and then `m=1` and (49) recover `d,z,y,x` after powers of `a` are inverted.
The failure of (52) to extend naively over `a=0` is concentrated over the
singular line of (51), namely `a=K=J=0`.  This is a cleaner global model for
the same modification phenomenon than the original division chart, but its
base is the singular trinomial (51), not `SL_2`.

The quotient itself does not force `a` to be invariant under additive
actions.  For example

\[
 E(a)=K^2,\qquad E(U)=192a^2,\qquad E(K)=E(J)=0       \tag{54}
\]

is an LND of `k[Q]`.  Its most direct lift through (52), obtained by fixing
`p,s`, does not extend across `a=0`: equations (53) and `d=(1-cp)/a` give

\[
 E(c)=-{576\over5}a^2p-4aK^2s,qquad
 E(d)={576\over5}ap^2+4K^2ps-{dK^2\over a}.          \tag{55}
\]

The pole cannot in fact be removed by a vertical correction.  Suppose a
global derivation `L` induced (54), and allow arbitrary regular functions
`P=L(p)` and `S=L(s)`.  Since

\[
 h={7J-12U\over40},\qquad
 c=2hp+{9K\over20}-2a^2s,
\]

one would have

\[
 L(c)=-{576\over5}a^2p+2hP-4aK^2s-2a^2S.
\]

On `a=0`, where `p=c^(-1)`, `K=4c`, and `h=-2c^2/5`, regularity of
`d=(1-cp)/a` forces

\[
 P=-80cd.
\]

Regularity of `z=(U+c^2)/(4a)` then forces

\[
 64c^2(2c^2d-z)=0.
\]

But `c` is invertible and `d,z` are independent coordinates on the boundary
`{a=0}`.  This is impossible.  Thus the LND (54) of `Q` admits **no** global
derivation lift to `X`, even with arbitrary regular `p,s` corrections.  This
does not rule out an LND moving `a` which fails to descend to the particular
quotient (51).

### Saturation of the visible actions and the Derksen algebra

The apparent rank drop of `D_(10),D_(22)` is not intrinsic.  The effective
coefficient scaling gives the coordinate ring the grading

\[
 \deg(a,c,z,p,d,y,x)=(4,3,2,-3,-4,-5,-6).
\]

The two displayed LNDs have degrees `5` and `11`, while `s` has degree `-5`.
Put

\[
 h=D_{22}(s)={7J-12U\over40}\quad\hbox{in }k[X],
 \qquad
 \Delta=D_{22}-hD_{10}.                              \tag{56}
\]

The function `h` is killed by both derivations, so `Delta` is again an LND
and `Delta(s)=0`.  More is true: `Delta` vanishes identically along `a=0`.
Indeed, there `cp=1`,

\[
 U=-c^2,\qquad K=4c,\qquad J=-4c^2,
 \qquad h=-{2c^2\over5},
\]

and direct substitution in (45)--(46) kills all seven components of
`Delta`.  The divisor `{a=0}` is irreducible, so every component is divisible
by `a` in `k[X]`.  Hence

\[
 \boxed{D_7={D_{22}-hD_{10}\over a}}                 \tag{57}
\]

is a global derivation of degree `7`.  It is locally nilpotent: `D_7(a)=0`
and `Delta=aD_7`, whence `Delta^n=a^nD_7^n`; local nilpotence of `Delta` and
the fact that `a` is a nonzerodivisor force local nilpotence of `D_7`.
This division is primitive.  At
`(a,c,z,p,d,y,x)=(0,1,0,1,0,0,-1)` one has

\[
 D_7(z)=-{2\over5},\qquad D_7(d)=-{1\over5}.
\]

Thus `D_(10),D_7` have rank two there, although `D_(10),D_(22)` have rank
one.  Equivalently,

\[
 D_{22}=aD_7+hD_{10},
\]

so the original pair generated a nonsaturated submodule of the tangent
sheaf along `a=0`.

This also determines the ordinary Derksen algebra.  Since `D_(10)(s)=1`, the
slice theorem says

\[
 k[X]=\ker(D_{10})[s].
\]

But `D_7(s)=0`, so `s` belongs to the kernel of the nonzero LND `D_7`.
Therefore the algebra generated by the kernels of all nonzero LNDs contains
both `ker(D_(10))` and `s`, and hence

\[
 \boxed{\operatorname{HD}(X_{2,3})=k[X_{2,3}].}      \tag{58}
\]

The Derksen algebra therefore agrees with that of
`A^2 x SL_2` and cannot decide the product question.  The full
Makar--Limanov invariant is still only bounded above by the intersections of
the known kernels; all three displayed LNDs fix `a`, but an undiscovered LND
may move `a,U,K`, or `J`.  Flexibility also remains open: saturation repairs
the artificial rank-one drop, but the known orbit distribution still has
rank only two at the displayed boundary point.

An earlier coefficient search found no third LND among ambient-preserving
homogeneous derivations whose generator images have total degree at most two.
The derivation (57) lies outside that ansatz: division takes place in the
quotient ring and its polynomial representatives have higher total degree.
Thus the search remains correct but is not a classification of homogeneous
LNDs.

Finally, a global birational morphism
`X_(2,3) -> A^2 x SL_2`, if constructed, would be an affine modification in
the sense of Kaliman--Zaidenberg; their theorem converts a birational affine
morphism into a modification, but does not supply that morphism.  Equations
(48)--(53) do not produce one.  They instead reduce the open problem to the
smooth fourfold `Y`.  The divisor `a=0` is still the only obstruction to the
generic `Ga^2` product chart (52), but it is no longer a Derksen obstruction:
the saturated action (57) extends nontrivially across it.  See
[`Affine modifications and affine hypersurfaces with a very transitive
automorphism group`](https://arxiv.org/abs/math/9801076),
[`Flexible varieties and automorphism groups`](https://arxiv.org/abs/1011.5375),
and [`Modified Makar--Limanov and Derksen invariants`](https://arxiv.org/abs/2212.05899)
for the general modification, flexibility, and slice-invariant frameworks.

## 12. Consequence

The natural `(2,3)` multiplication map is an etale, generically degree-ten
map

\[
 U_{2,3}\longrightarrow\mathbb A^5,
\]

but its source is not `A^5`.  It therefore does **not** produce a new `JC(5)`
counterexample.  It is still a genuinely new nonproper etale cover model.
Its Grothendieck class, integral Chow groups, algebraic `K_0`, fundamental
group, integral cohomology, and ordinary homotopy type all have the
`A^2 x SL_2` shape.  The Euclidean chart is explicitly an affine modification
of that product, but the modification morphism does not extend through the
second coefficient chart.  Globally, (46)--(48) exhibit `X_(2,3)` as an
`A^1`-cylinder and (49)--(52) identify a generic `Ga^2` quotient by a singular
trinomial threefold.  A different modification of `A^2 x SL_2` has not been
constructed.  The ordinary Derksen algebra is the full coordinate ring by
(58), exactly as for the proposed product; the Makar--Limanov invariant and
flexibility remain open.  What is localized to `a=0` is the failure of the
generic `Ga^2` product chart to extend, not a failure of the additive kernels
to generate the ring.

The exact class calculation and direct finite-field enumeration are checked
by
[`verify_quadratic_cubic_factorization_invariants.py`](../scripts/verify_quadratic_cubic_factorization_invariants.py).
The determinant fiber and rank-jump calculations in Section 6 are checked by
[`verify_quadratic_cubic_product_test.py`](../scripts/verify_quadratic_cubic_product_test.py).
The normalization product and a transverse point of the two affine boundary
components used in Section 7 are checked by
[`verify_quadratic_cubic_fundamental_group_inputs.py`](../scripts/verify_quadratic_cubic_fundamental_group_inputs.py).
The Euclidean modification, complementary-chart transition, and integral
residue coefficients in Sections 9--10 are checked by
[`verify_quadratic_cubic_modification_topology.py`](../scripts/verify_quadratic_cubic_modification_topology.py).
The two commuting LNDs, global slice, trinomial quotient, and the original
nonsaturated rank drop in Section 11 are checked by
[`verify_quadratic_cubic_additive_actions.py`](../scripts/verify_quadratic_cubic_additive_actions.py).
The degree-seven saturated LND, its boundary rank, the full Derksen algebra,
and the no-lift obstruction for (54) are checked by
[`verify_quadratic_cubic_saturated_lnd.py`](../scripts/verify_quadratic_cubic_saturated_lnd.py).

The determinant-one calculation, normalization uniqueness, Picard
vanishing, and canonical-class identity are consolidated for arbitrary
factor degrees in the
[relative-scaling/boundary-lattice theorem](RELATIVE_SCALING_BOUNDARY_LATTICE.md).
