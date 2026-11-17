# The scheme-theoretic decorated-normalization invariant

This note upgrades the reduced Hessian-root support to a scheme-theoretic
stable invariant, computes the off-diagonal self-intersection of the
discriminant normalization, and applies the resulting decoration to the two
previously unresolved split quartic islands.

The core results `D1` and `F2` use only the smaller marked Hessian-divisor
invariant `(P^1; div(H''), 0, infinity)` and, for `F2`, the distinguished
affine root sheet.  Node pairing and conductor are retained here as the
strictly stronger `D2` invariant; they are not load-bearing in the generic
dimension, unramifiedness, exact rerooting-degree, or affine-mark proofs.

Work over an exact field `k` of characteristic zero; geometric point counts
are taken after extension to an algebraic closure.  Irreducible factors are
retained over `k`, so their degrees record residue extensions rather than
being replaced by numerical roots.  Let

\[
 E_{s,t}(W)=H(W)-sW+t,
 \qquad
 \nu_H(r)=\bigl(H'(r),rH'(r)-H(r)\bigr)
\]

and let `D_H` be the reduced discriminant curve.  For an exact-degree
primitive in the generic-discriminant open, `nu_H:A^1 -> D_H` is its affine
normalization; its projective completion is `P^1`, with one smooth point at
infinity.

## 1. The full Fitting divisor

Put

\[
 R_H=k[H'(r),rH'(r)-H(r)]\subset k[r].
\]

The transitivity sequence for differentials gives

\[
 \Omega_{k[r]/R_H}
 =k[r],dr/\bigl(dH'(r),d(rH'(r)-H(r))\bigr).
\]

Since

\[
 dH'=H''(r)dr,
 \qquad
 d(rH'-H)=rH''(r)dr,
\]

there is an equality of modules, not merely an equality of supports,

\[
 \Omega_{k[r]/R_H}\simeq k[r]/(H''(r))\,dr.
\]

Consequently

\[
 \boxed{\operatorname{Fitt}_0
 \Omega_{\widetilde D_H/D_H}=(H''(r)).}                 \tag{1.1}
\]

Thus a root of `H''` of order `m` occurs with coefficient `m` in the Fitting
divisor.  Passing to `sqf(H'')` loses genuine invariant information.

### Weighted intrinsic-boundary reconstruction

Assume that zero is an exact double root of `H`, every nonzero root is
simple, and the weighted reconstruction formulas are taken on their maximal
regular open.  Form, directly from the map,

\[
 \mathcal B(F_H)=
 (\overline X_{F_H}\to Y,
  \mathbb A^3\hookrightarrow\overline X_{F_H},
  \partial_{F_H})
\]

as in the map-intrinsic clause of the
[stable normalization theorem](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md).
The complete target images of its boundary divisors form an intrinsically
ordered pair `(Z_Delta,Z_0)`: `Z_Delta` is the unique member receiving a
ramified boundary prime, and `Z_0` is the other member.

Indeed, on `C!=0`, regular reconstruction is equivalent to `E_W!=0`.  At the
generic point of the reduced discriminant, the normalized double-root branch
is outside the etale source and has `(e,f)=(2,1)`; the other `N-2` roots are
simple affine branches.  The degree sum is therefore `2+(N-2)=N`.  At the
generic point of `C=0`, the exact double root at zero supplies the
residue-degree-two affine cluster and root one supplies a residue-degree-one
affine branch.  For another simple root `rho`, the identity

\[
 y=(W-\gamma)/C,
 \qquad \gamma\longrightarrow-H'(\rho)/c
\]

has a pole unless the numerator cancels; in that exceptional case the
numerator of the reconstruction formula for `z` tends to `rho-1`, so the
branch is still outside the source.  The remaining `N-3` branches are thus
unramified boundary branches, and `2+1+(N-3)=N`.  Off these two target
divisors every root is simple and reconstructs regularly.  The degree sums,
therefore, exclude hidden height-one boundary primes over either divisor or
over any other target divisor.  A boundary component supported only in
higher codimension is also impossible.  After deleting the listed divisors,
choose an affine neighborhood `V` of its hypothetical generic point.  The
intersection of the affine opens `V` and `A^3` is affine because the
normalization is separated, while normal Hartogs extension gives

\[
 \Gamma(V,O)=\Gamma(V\cap A^3,O).
\]

The corresponding inclusion of affine schemes would be an isomorphism, a
contradiction.  Thus the two target vertices and their displayed branches
are complete.

On `C!=0`, the target change

\[
 (A,B,C)\longmapsto(s,t,C)=(BC,cAC^2,C)
\]

therefore identifies the intrinsic finite stratum map with

\[
 \widetilde D_H\times G_m\longrightarrow D_H\times G_m. \tag{1.1a}
\]

This is the family-specific step: the pencil coordinates calculate the
stratum selected from `mathcal B(F_H)`; they do not define the invariant.

### Stable functoriality, including multiplicities

All base-change assertions in this subsection are applications of the
construction-independent
[stable normalization functoriality theorem](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md).
The present note supplies the family-specific identification of the intrinsic
decorated normalization.

For two maps on this weighted clean locus, the preceding proposition
reconstructs the ordered target-boundary pair `(Z_Delta,Z_0)` from the maps.
Removing `Z_0` identifies the discriminant vertex and its normalization with

\[
 D_H\times G_m,
 \qquad
 \widetilde D_H\times G_m.
\]

The map-intrinsic functoriality theorem says that a polynomial left--right
equivalence identifies these finite normalization morphisms.  It therefore
identifies their relative differential modules and their zeroth Fitting ideal
sheaves.  After adjoining `q` identity variables,

\[
 \Omega_{(\widetilde D_H\times A^q)/(D_H\times A^q)}
 \simeq \operatorname{pr}_1^*Omega_{\widetilde D_H/D_H}.
\]

Fitting ideals commute with this flat base change by that theorem.  Hence the effective
Cartier divisor (1.1), including every coefficient
`ord_alpha(H'')`, is preserved under stable polynomial left--right
equivalence.

If the divisor has at least two distinct support points, the unit argument on
`k[r,xi,xi^-1,z_1,...,z_q]` from the degree-five theorem forces the induced
normalization coordinate change to be affine.  The argument applies to the
ideal powers as well as their radicals.  Therefore stable equivalence forces

\[
 \sum_\alpha \operatorname{ord}_\alpha(H'')[\alpha]
 \quad\text{and}\quad
 \sum_\beta \operatorname{ord}_\beta(G'')[\beta]
\]

to be affinely equivalent as effective divisors.  The doubled Hessian root in
the easy degree-five witness is therefore retained with coefficient two.

## 2. The off-diagonal self-intersection scheme

Use two normalization coordinates `r,u`.  The fiber product
`widetilde D_H times_D_H widetilde D_H` is cut out by

\[
 H'(r)-H'(u)=0,
\]

\[
 rH'(r)-H(r)-uH'(u)+H(u)=0.                            \tag{2.1}
\]

Both equations contain the diagonal factor `r-u`.  Define

\[
 D_s(r,u)=\frac{H'(r)-H'(u)}{r-u},
\]

\[
 D_t(r,u)=
 \frac{rH'(r)-H(r)-uH'(u)+H(u)}{r-u}.                 \tag{2.2}
\]

The closure of the genuine off-diagonal scheme is

\[
 \boxed{
 \Gamma_H^{\mathrm{off}}
 =V\bigl((D_s,D_t):(r-u)^\infty\bigr).}               \tag{2.3}
\]

The saturation in (2.3) cannot be omitted.  On the diagonal, the divided
equations specialize to `H''(r),rH''(r)` and retain isolated cusp points.  For
a cubic this would falsely report a double point even though the
off-diagonal scheme is empty.

The involution `(r,u)<->(u,r)` is free on the ordinary off-diagonal locus.
The characteristic-zero finite quotient exists without this freeness; in a
degenerating family its closure can acquire fixed diagonal limits, which must
be retained scheme-theoretically rather than interpreted as ordinary nodes.
In unordered coordinates

\[
 \sigma=r+u,\qquad \pi=ru,
\]

symmetrizing (2.2) and saturating by `sigma^2-4pi` gives the node-pairing
scheme.  For a generic degree-`N` seed it is reduced of length

\[
 \frac{(N-2)(N-3)}2,
\]

while its ordered pullback has length `(N-2)(N-3)`.  Each point records the
two normalization branches belonging to one node, information not present in
the singular support of `D_H` alone.

The fiber product, its diagonal, and hence the saturated off-diagonal closure
are functorial under an isomorphism of normalization morphisms.  Flat
stabilization only takes their product with affine space.  Thus the pairing
scheme is another stable invariant.

## 3. Conductor and the decorated object

The conductor is intrinsically

\[
 \mathfrak c_H=
 \operatorname{Ann}_{\mathcal O_{D_H}}
 \bigl(\nu_{H*}\mathcal O_{\widetilde D_H}/
             \mathcal O_{D_H}\bigr),                 \tag{3.1}
\]

together with its extension to the normalization.  It is functorial under
isomorphism and flat base change.  On the ordinary cusp/node locus, its
pullback to `k[r]` is especially simple.  A cusp `k[[t^2,t^3]]` contributes
`(t^2)` and a branch of a node contributes `(t)`.  If `K_H=H''` is
squarefree and `P_H` is the monic projection polynomial of the ordered node
scheme to either normalization factor, then

\[
 \boxed{\mathfrak c_H k[r]=(K_H(r)^2P_H(r)).}          \tag{3.2}
\]

There is also a uniform formula which does not assume ordinary
singularities.  Let `f_H(S,T)=0` be the primitive implicit equation of the
plane image.  Adjunction for a normalized plane curve gives

\[
 \boxed{
 \bar{\mathfrak c}_H(r)=
 \frac{(\partial f_H/\partial T)
   (H'(r),rH'(r)-H(r))}{H''(r)}.}                     \tag{3.2a}
\]

The quotient is a polynomial, up to a nonzero scalar, whenever the displayed
parametrization is the normalization.  Formula (3.2a) computes the full
conductor for collided cusps, tacnodes, and worse plane-curve singularities;
on the ordinary locus the implementation independently checks that it equals
(3.2).

The invariant is the finite conductor map, not only the polynomial in (3.2).
Writing `A_H=R_H`, `B_H=k[r]`, and
`cbar_H=mathfrak c_H B_H`, retain

\[
 \bar C_H=\operatorname{Spec}(B_H/\bar{\mathfrak c}_H)
 \longrightarrow
 C_H=\operatorname{Spec}
 \bigl(A_H/(A_H\cap\bar{\mathfrak c}_H)\bigr).          \tag{3.3}
\]

The implementation obtains a downstairs presentation by eliminating `r`
from

\[
 \bar{\mathfrak c}_H(r),\qquad
 S-H'(r),\qquad T-rH'(r)+H(r).                          \tag{3.4}
\]

For every ordinary quartic audit, the upstairs algebra in (3.3) has length
six and its downstairs image algebra has length three.  The finite map,
together with the off-diagonal relation, distinguishes the two preimages of a
node from the single thick preimage of a cusp.

On a locus where the ordered boundary pair is intrinsic and every displayed
mark has the intrinsic characterization stated below, the scheme-theoretic
decorated normalization is the tuple

\[
 \mathcal Z(H)=
 \bigl(
  P^1;
  \operatorname{Fitt}_0\Omega;
  \Gamma_H^{\mathrm{off}}/S_2;
  \bar C_H\to C_H;
  \infty;
  \mathbf b_H
 \bigr),                                               \tag{3.5}
\]

where the full center divisor underlying `b_H`, when retained, is the divisor
reconstructed from the finite boundary incidence; in the pencil coordinate
its executable presentation is `gcd(H,H')`, including multiplicities.  The
polynomial gcd alone does not define an invariant.  Point labels are inferred
only when the boundary incidence has a uniquely certified support point; the
extractor never invents a mark at zero for an arbitrary polynomial.  For a
normalized boundary-clean split weighted seed, the canonical zero-cluster
mark is `r=0`.
Further intrinsically labeled boundary branches may be retained as additional
marks.  When the Fitting divisor has at least two distinct support points, the
unit argument above forces every stabilized normalization-coordinate change
to be affine and hence also preserves the unique completion point `infinity`.
Without this two-support condition, infinity must be independently
characterized before it is retained.  Subject to these explicit mark
hypotheses, every entry in (3.5) is constructed from the map-intrinsic
normalization package and its finite stratum maps, so stable left--right
equivalence preserves the decorated isomorphism class.

The assertion that `r=0` is canonical is an exact boundary-chart statement,
not a choice of polynomial coordinates.  Assume that zero is the unique
multiple root of `H`, of exact multiplicity two, and write
`H(r)=h_2r^2+O(r^3)`.  On the normalized discriminant divisor,

\[
 B=\frac{H'(r)}C,
 \qquad
 A=\frac{rH'(r)-H(r)}{cC^2}.
\]

The chart `r=Ck` has specialization

\[
 B=2h_2k,
 \qquad
 A=\frac{h_2}{c}k^2,
 \qquad
 B^2=4ch_2A.                                          \tag{3.6}
\]

Conversely, a finite limit at `C=0` with `r -> rho` forces
`H'(rho)=H(rho)=0`; the boundary-clean hypothesis therefore forces
`rho=0`.  Thus the component of the normalized discriminant meeting the
second boundary has valuation center `r=0`.  An isomorphism of the ordered
boundary pair transports this center, so the mark in (3.5) is intrinsic.
The symbolic checker evaluates (3.6) for both split quartics.

## 4. The split quartic islands

For the two unresolved split examples,

\[
 H_b=W^2(W-1)(W-3),
 \qquad
 H_c=-W^2(W-1)(2W+1).
\]

Their monic Fitting divisors are

\[
 K_b=W^2-2W+\frac12,
 \qquad
 K_c=W^2-\frac14W-\frac1{12}.                         \tag{4.1}
\]

Exact diagonal saturation gives one unordered node pair in each case:

\[
 H_b:\quad \sigma=2,\quad\pi=-\frac12,
\]

\[
 H_c:\quad \sigma=\frac14,\quad\pi=-\frac9{32}.     \tag{4.2}
\]

The cusp pair and node pair by themselves still do not separate the
examples.  In both cases the two pairs have the same midpoint and the squared
node radius is three times the squared cusp radius.  Indeed

\[
 r\longmapsto \frac18+a(r-1),\qquad a^2=\frac{19}{96},
\]

identifies the cusp and node-pair decorations (up to swapping each pair).

The canonical zero-cluster boundary mark breaks this coincidence.  For a
quadratic divisor with roots `m+d,m-d` and a marked point `b`, set

\[
 I(K,b)=\frac{(b-m)^2}{d^2}.                            \tag{4.3}
\]

This is invariant under affine changes of the normalization coordinate.  At
the boundary mark `b=0`, exact calculation gives

\[
 I(K_b,0)=2,
 \qquad
 I(K_c,0)=\frac3{19}.                                  \tag{4.4}
\]

Therefore the decorated normalizations are not isomorphic and

\[
 \boxed{F4b\not\sim_{\mathrm{stable}}F4c.}             \tag{4.5}
\]

This resolves the quartic pair left open by the reduced boundary signature
and by the unmarked two-point Hessian support.

## 5. The minimal moduli map and its generic dimension

Let `A_N` be the normalized admissible seed space

\[
 H(W)=\sum_{j=3}^N h_j(W^j-W^2),
 \qquad H'(1)=-1,
\]

with the degree-drop and weighted exceptional divisors removed.  It has
dimension `N-3`.  Put `d=N-2`, and let `U_d` be the open of squarefree
effective degree-`d` divisors on `P^1` disjoint from zero and infinity.
Scaling the normalization coordinate gives the quotient

\[
 \mathcal Q_d=[U_d/G_m],\qquad r\longmapsto ar,
\]

and the minimal rational moduli map

\[
 \mathfrak h_N:\mathcal A_N\dashrightarrow\mathcal Q_d,
 \qquad H\longmapsto(P^1;\operatorname{div}(H''),0,\infty). \tag{5.1}
\]

On the smaller ordinary cusp-node open, the full decoration (3.5) refines
`mathfrak h_N` by adding node pairing and conductor.  That refinement is D2;
the D1/F2 arguments below concern `mathfrak h_N`.

This further open is nonempty for every `N>=4`; the ordinary and
boundary-clean conditions are not merely formal genericity assumptions.  The
[generic discriminant theorem](GENERIC_DISCRIMINANT_CURVE.md) proves that the
ordinary nodal-cuspidal open meets the normalized admissible slice in every
degree.  For the boundary condition, use the explicit normalized split seed

\[
 H_N^{\mathrm{split}}(W)
 =\frac12W^2(1-W)(1+W^{N-3}).                         \tag{5.1a}
\]

It has exact degree `N`, satisfies

\[
 (H_N^{\mathrm{split}})'(1)=-1,
 \qquad
 (H_N^{\mathrm{split}})''(1)=-(N+1)\ne-2,
 \qquad
 (H_N^{\mathrm{split}})''(0)=1.                       \tag{5.1b}
\]

Writing `P_N=2H_N^split/W^2=(1-W)(1+W^(N-3))`, the polynomial
`P_N` is squarefree and has no zero root.  At any root `alpha` of `P_N`,

\[
 \left.\frac{2(H_N^{\mathrm{split}})'(W)}W
 \right|_{W=\alpha}
 =\alpha P_N'(\alpha)\ne0.                            \tag{5.1c}
\]

Hence

\[
 \gcd(H_N^{\mathrm{split}},(H_N^{\mathrm{split}})')=W,
\]

so zero is the unique multiple primitive root.  Moreover the boundary mark
`r=0` is not a cusp by (5.1b).  If it shared its discriminant image with an
off-diagonal point `u`, then `H'(u)=H(u)=0`, contradicting (5.1c).  Thus it is
not a node branch either.

Boundary-cleanness and avoidance of the cusp/node schemes are open
conditions.  The normalized seed space is irreducible, the ordinary open is
nonempty by the generic discriminant theorem, and (5.1a) proves that the
boundary-clean open is nonempty.  Their intersection is therefore a nonempty
open for every `N>=4`.

The quotient `mathcal Q_d` has dimension

\[
 d-1=N-3.                                             \tag{5.2}
\]

This equals `dim A_N`.  In fact the generic fiber is finite.  Suppose the
Fitting divisors and the marks at zero and infinity of two normalized seeds
`H,G` are isomorphic.  The normalization-line isomorphism has the form
`r -> ar`, and equality of effective Fitting divisors gives

\[
 G''(r)=cH''(ar)
\]

for some nonzero `c`.  Since every seed in `A_N` satisfies
`H(0)=H'(0)=G(0)=G'(0)=0`, integration gives

\[
 G(r)=\frac{c}{a^2}H(ar).                              \tag{5.3}
\]

The endpoint equation `G(1)=0` forces `H(a)=0`.  A generic degree-`N` seed
has only finitely many simple nonzero roots, and for each such `a` the
normalization `G'(1)=-1` determines `c` uniquely.  Hence a generic fiber has
at most `N-2` points.  Consequently

\[
 \boxed{\dim\overline{\operatorname{im}(\mathfrak h_N)}=N-3}
\]

on the generic Hessian-clean locus.

### Generic unramifiedness and the exact rerooting fiber

The preceding finite-fiber argument has a differential strengthening using
only the Fitting divisor and the marks at zero and infinity.  Let
`H+epsilon dot H` be a normalized first-order
seed deformation in its tangent kernel.  Infinitesimal equality of the
projective Fitting divisors, modulo normalization-line scaling, says

\[
 \dot H''(r)=\alpha H''(r)+\beta rH'''(r).             \tag{5.3g}
\]

Here `alpha` is the scalar ambiguity in a defining section and `beta r d/dr`
is the infinitesimal automorphism of `(P^1;0,infinity)`.  Integrating twice
and using `dot H(0)=dot H'(0)=0` gives

\[
 \dot H=\alpha H+\beta(rH'-2H).                        \tag{5.3h}
\]

The remaining normalized endpoint conditions give, successively,

\[
 0=\dot H(1)=-\beta,\qquad
 0=\dot H'(1)=-\alpha.
\]

Thus `dot H=0`.  The map `mathfrak h_N` is unramified; on the
generic locus where the divisor has trivial scaling stabilizer, source and
target are smooth of the same dimension `N-3`, so it is etale after
shrinking.

The generic degree is also exact.  If `a` is a nonzero simple root of `H`, set

\[
 \kappa_a=-\frac1{aH'(a)},\qquad G_a(w)=\kappa_aH(aw). \tag{5.3i}
\]

Then `G_a` is normalized and

\[
 E_{G_a}(w;s,t)=
 \kappa_a E_H\!\left(aw;
       \frac{s}{\kappa_a a},\frac{t}{\kappa_a}\right). \tag{5.3j}
\]

Since `G_a''(w)=kappa_a a^2H''(aw)`, all `G_a` have the same marked Hessian
divisor as `H`.  Generically `H=w^2P` with `P` squarefree of degree
`N-2`, and the `N-2` rerootings are distinct.  Conversely, (5.3) shows that
every seed in the same Hessian-divisor fiber is one of these rerootings.
Therefore

\[
 \boxed{\mathfrak h_N\text{ is generically etale of degree }N-2
        \text{ onto its image}.}                      \tag{5.3k}
\]

Equivalently, over a rerooting-stable generic open the incidence
`{(H,a):H(a)=0, a!=0, H'(a)!=0}` is the finite etale rerooting relation.
Its quotient supplies an ordinary algebraic model for the generic moduli
image, which is therefore reduced.  The affine root sheet separates the
rerootings, as proved next.  In fact node pairing, the complete conductor
morphism, and full Fitting multiplicities cannot replace that sheet: equation
(5.3j) transports the entire decorated finite normalization.  Exact generic
quartic/quintic examples and the first nonordinary Hessian-collision
calculation are recorded in the
[conductor-rerooting counterexample](CONDUCTOR_REROOTING_AMBIGUITY.md).

### Generic affine-mark rigidity

Over `Z_0=V(C)`, the normalized incidence intrinsically distinguishes the
simple root `W=1` from the additional simple roots.  The distinction is a
stratum of the canonical normalized cover: it is the unique simple sheet in
the distinguished affine open outside the double-zero cluster.  No recovery
of a primitive pencil from a pole filtration is needed.

Assume throughout this paragraph that zero has exact multiplicity two, as it
does on the generic open in `A_N`, and write
`H(W)=h_2W^2+O(W^3)` with `h_2!=0`.  Over the generic point
`K=k(A,B)` of `Z_0`, the finite zero-cluster chart `W=CR` has special equation

\[
 h_2R^2-BR+A=0.                                      \tag{5.3a}
\]

Its discriminant `B^2-4h_2A` has odd valuation along its own prime, so (5.3a)
is irreducible over `K`; this affine component has residue degree two.  The
root `W=1` is simple because `H'(1)=-1`, and the completed reconstruction
calculation gives a regular affine component of residue degree one.  Every
additional simple primitive root is a polar boundary component.  Hence:

> **Intrinsic affine-stratum fact.**  On the exact-double-zero,
> boundary-clean locus, the root-one component is the unique
> residue-degree-one component over the generic point of `Z_0` which is
> contained in the distinguished affine open and is not part of the
> zero-cluster component.

This statement is preserved by stable left--right equivalence because the
intrinsic cover retains the open immersion `X -> bar X_F`, not merely its
function field.

It is therefore safe to refine the marked Hessian divisor by adjoining the

\[
 \mathcal A_{1,H}\longrightarrow Z_0,                \tag{5.3aa}
\]

where `mathcal A_(1,H)` denotes this distinguished affine component together
with its embedding in the regular-reconstruction open.  This is the intrinsic
version of the “root-one branch.”  It is a genuine stable invariant and a
stratum of the normalized cover.

The ramification divisor in the incidence is cut out by
`E=partial_W E=0`.  On the root-one affine component,

\[
 C=0,qquad W=1,qquad
 \partial_WE=H'(1)-BC=-1.                            \tag{5.3b}
\]

Thus this component is disjoint from the ramification divisor.  The new proof
does not attempt to recover an arbitrary primitive element from a stabilized
plane automorphism.  Instead it first uses the exact coarse fiber theorem to
reduce every possible identification to one of the finite rerootings below,
and then tests which rerooting preserves this intrinsic affine component.

The ambiguity is visible algebraically.  If `a` is any nonzero simple root of
`H`, define

\[
 \kappa_a=-\frac1{aH'(a)},
 \qquad G_a(w)=\kappa_aH(aw).                         \tag{5.3c}
\]

Then `G_a(1)=0` and `G_a'(1)=-1`, and the bare marked-root pencils satisfy

\[
 E_{G_a}(w;s,t)
 =\kappa_a E_H\!\left(aw;
        \frac{s}{\kappa_a a},\frac{t}{\kappa_a}\right). \tag{5.3d}
\]

Their discriminant normalizations are therefore linearly identified by
`r_H=ar_G`.  Under this identification the distinguished root of `G_a` goes
to the root `a` of `H`, not necessarily to `1`.  When `a!=1`, (5.3d) does not
identify the distinguished regular-reconstruction opens: it sends the affine
root-one component for `G_a` to an extra-root boundary component for `H`.
This shows exactly where the affine mark contains more information than the
marked Hessian divisor.

> **Theorem (generic affine-mark faithfulness).**  On the Hessian-clean
> open, the marked Hessian divisor together with the unique
> unramified root sheet meeting the regular-reconstruction open recovers the
> normalized seed exactly.  After shrinking, the normalized seed open is
> isomorphic to the normalization of the affine-marked Hessian-divisor image, and
> forgetting the affine sheet is finite etale of degree `N-2`.

The proof is now finite.  By the exact fiber theorem, two seeds with the same
marked Hessian divisor differ by one of the rerootings (5.3c).  Under (5.3d), the
root-one sheet of `G_a` maps to the root-`a` sheet of `H`.  The tangent-core
reconstruction identities are

\[
 \gamma=-E_W/c,\qquad x=C/\gamma,
 \qquad y=(W-\gamma)/C.                               \tag{5.3e}
\]

The double-zero sheet is ramified.  For a nonzero simple root `rho`, the
limit of `gamma` is `-H'(rho)/c`.  If this differs from `rho`, the last
expression in (5.3e) has a pole.  If it equals `rho`, the `y`-pole cancels,
but the numerator of the reconstruction formula for `z` tends to `rho-1`.
Thus root one is the unique unramified sheet meeting the affine open.
Preservation of that sheet forces `a=1`, and then `kappa_a=1`, so `G=H`.
Generic unramifiedness plus this one-to-one statement and Zariski's Main
Theorem identify the seed open with the normalization of its marked image.

The former stable-plane descent, low-pole filtration, fixed-curve, and deck-
rigidity proof is no longer an active theorem source.  Its identities remain
available as independent audits, but the generic result uses only the exact
rerooting fiber, tangent-core reconstruction, and stable preservation of the
distinguished affine open.

### Stacky marked extension (complete)

Put `n=N-2`.  On the collision-separating admissible-cover compactification,
the selected root is the tautological marked point `p_1`, and forgetting that
selection is the quotient-stack morphism

\[
 [\overline M_{0,N}/S_{n-1}]
   \longrightarrow[\overline M_{0,N}/S_n].            \tag{5.3q}
\]

This morphism is finite, representable, and etale of degree
`[S_n:S_(n-1)]=n=N-2`.  Thus the selected root already extends across every
collision on the marked stack.

The proper Hurwitz-space package behind this statement—including the LL
branch incidence, collision charts, and conductor boundary—is constructed in
[HURWITZ_LL_COMPACTIFICATION.md](HURWITZ_LL_COMPACTIFICATION.md).  Its formal
comparison theorem identifies the normalized Stein factors with the canonical
root cover and discriminant normalization after contraction.  It also
identifies every completed root-cover chart and the normalization--conductor
square at simultaneous collisions of arbitrarily many distinct root
clusters.  If the multiple roots have `e_i=ord_(rho_i)(H)-1` and
`E=sum_i e_i`, the conductor on the `i`-th completed normalization branch is
exactly

\[
 (u_i^{e_i(E-1)}).                                    \tag{5.3r}
\]

Its boundary pullback splits with coefficient one according as the selected
point lies inside or outside a collision block, with restriction degrees `k`
and `n-k`.  Only after contraction to coefficient space does ramification
appear: a generic pair collision has transposition inertia.  The
`epsilon=delta^m` chart at a total `m`-fold collision is a cyclic
codimension-`m-1` slice, not generic divisorial inertia.

This is the completed stacky marked-extension theorem `H2`.  In particular,
arbitrary collision geometry, normalized-Stein comparison, and conductor
gluing are not part of the remaining open problem.

### Coarse affine-mark descent (complete)

Let `bar S_N` and `bar B_N` be the coarse spaces of the marked and unmarked
collision-separating compactifications.  They are normal: étale-locally they
are finite-group quotients of the smooth space `bar M_(0,N)`, and invariant
subrings of normal domains are normal.  The selected point gives a finite
coarse root cover

\[
 \pi_c:\overline S_N\longrightarrow\overline B_N.   \tag{5.3s}
\]

The contraction morphism from the marked admissible-cover stack to this
algebraic-space cover factors uniquely through its coarse moduli space.  By
the stacky comparison theorem `H2`, its pullback to every completed
repository chart is the canonical finite root cover, and its generic selected
point is exactly the distinguished degree-one component contained in the
regular-reconstruction open.

> **Theorem (coarse affine-mark descent).**  In the relative Hilbert space of
> degree-one subschemes of the contracted finite root cover, the closure `Z`
> of the point selected generically by the regular-reconstruction component is
> finite and separated, and the map `Z -> bar S_N` is an isomorphism.
> Consequently the affine mark has exactly one coarse specialization over
> every DVR limit, independently of the admissible-cover lift and compatibly
> with base change.

There are two formal ingredients.  For a finite separated cover, its relative
Hilbert space of length-one subschemes is canonically the cover itself.  Thus
`Z` is a closed integral subspace of a finite cover and is finite over
`bar S_N`; it is generically the graph of the distinguished mark, so it is
birational.  A finite birational
morphism to a normal integral space is an isomorphism
([Stacks, Lemma 29.55.8](https://stacks.math.columbia.edu/tag/0AB1)).
Second, the coarse moduli map is categorical for morphisms to algebraic
spaces ([Stacks, Section 106.12](https://stacks.math.columbia.edu/tag/0DUF)),
so the contracted stack mark factors uniquely through `bar S_N`.
Separatedness of the finite cover then gives DVR uniqueness.

The collision algebra makes the conclusion explicit.  At a cluster of
`mu` roots, write `R=k[x_1,...,x_mu]`.  The coarse unmarked and marked rings
are

\[
 R^{S_\mu}=k[e_1,\ldots,e_\mu],\qquad
 R^{S_{\mu-1}}
 \simeq k[e_1,\ldots,e_\mu,T]/
 \bigl(T^\mu-e_1T^{\mu-1}+\cdots+(-1)^\mu e_\mu\bigr), \tag{5.3t}
\]

where `T=x_1` is the selected root.  Indeed, if `e'_j` is elementary
symmetric in `x_2,...,x_mu`, the recursion
`e'_j=e_j-T e'_(j-1)` recovers every marked invariant.  Thus (5.3t) is the
universal monic-root incidence, finite flat of degree `mu`.  At a total
collision its fiber is `k[T]/(T^mu)`: it retains the stacky length and
monodromy but has exactly one geometric point.  In particular neither
`mu=2` nor `mu=3` gives a counterexample.  The same calculation handles every
simultaneous cluster independently, and the completed-chart comparison in
`H2` identifies these local covers with the repository contraction.

The normalized-graph construction has this same Hilbert point as its coarse
scheme-theoretic image; the reconstruction open is used to select the generic
point, not to require the limiting point to remain in that open.  For the
nonreduced collision fiber the limit is supported at its sole geometric
point.  This completes the former open problem `OP-MARK`.

Projectivizing the `N-2` nonzero critical values gives a natural matching LL
stratum of the same dimension as the seed space:

\[
 \Lambda_N:\mathcal A_N^\circ\longrightarrow
 \mathbb P(1,2,\ldots,N-2).
\]

The classical polynomial LL count, divided by the generic `mu_N` source
scaling and multiplied by the `N-2` choices of marked unramified point over
zero, gives

\[
 \boxed{\deg\Lambda_N=(N-2)N^{N-3}}.                 \tag{5.3u}
\]

On the normalized graph compactification its weighted hyperplane
intersection is `N^(N-3)/(N-3)!`.  Thus the formerly heuristic connection to
LL degrees is now an exact marked Hurwitz count.

Because the decorated normalization is constant on stable polynomial
left--right classes, generic finiteness has the following stronger
degreewise consequence:

\[
 \boxed{\text{For every }N\ge4,\text{ weighted degree-}N\text{ maps contain
 an }(N-3)\text{-dimensional family of stable classes.}}                 \tag{5.4}
\]

Here “dimension” means the dimension of the image of the generically finite
decorated-invariant map.  The statement does not require the existence of a
coarse moduli space for all Keller maps.

## 6. The full finite-cover boundary invariant

Return now to a dominant quasi-finite Keller map `F:X=A^3 -> Y=A^3`.  The
coordinate-free starting object is

\[
 \mathcal B(F)=
 \bigl(\bar X_F\xrightarrow{\pi_F}Y,
       X\lhook\joinrel\longrightarrow\bar X_F,
       \partial_F\bigr),                                \tag{6.1}
\]

where `bar X_F=Norm_Y(k(X))` and
`partial_F=(bar X_F-X)_red`.  Left--right equivalence transports (6.1), and
normalization commutes with polynomial extension, giving

\[
 \mathcal B(F\times\mathrm{id}_{A^q})
 \simeq\mathcal B(F)\times A^q.                         \tag{6.2}
\]

The map-intrinsic theorem proves the transport assertion here, including the
open immersion; it is not an appeal to uniqueness of the finite cover alone.
The active computable invariant is organized in three functorial layers.
An item in these layers is retained only after its stratum and any label on
it have been characterized uniquely from (6.1).

1. The reduced layer retains the ordered target divisors, every upstairs
   boundary prime and finite edge map, geometric and arithmetic `(e,f)`, and
   every reduced multiple intersection.
2. The scheme layer retains the full upstairs and downstairs intersection
   schemes, their scheme images, and all specialization arrows.
3. The formal layer retains completed stratum maps, the relative different
   `Fitt_0 Omega_(bar X_F/Y)`, finite-stratum conductors, valuation
   filtrations, and derived local monodromy on punctured strict-henselian
   neighborhoods.

There is deliberately no conductor assigned to the open immersion
`X -> bar X_F`: a localization can have zero conductor.  Conductors are used
only for finite schematically dominant stratum maps, such as (3.3).

For a boundary-clean weighted primitive, the executable generic divisorial
profile returns:

- `E_Delta -> Z_Delta` with `(e,f)=(2,1)`, different exponent one,
  completed tame-DVR equation `delta=unit*q^2`, and transposition inertia;
- for an arithmetic extra-root factor of degree `f` and multiplicity `mu`,
  a prime over `Z_0` with `(e,f)=(mu,f)`, different exponent `mu-1`, completed
  equation `C=unit*q^mu`, and `mu`-cycle inertia after geometric splitting;
- for a zero cluster of multiplicity `m>=3`, a prime with
  `(e,f)=(m-1,1)`, different exponent `m-2`, and the analogous
  `(m-1)`-cycle; and
- the affine and boundary contributions whose sum is the full inverse degree
  over each target divisor.

Characteristic zero makes these extensions tame, so `different=e-1`.  The
finite edge data retains irreducible factor polynomials, not only their
degrees, and hence retains residue fields over nonclosed bases.

The generic divisorial layer is fully extracted and tested.  The abstract
scheme/formal layers are defined and proved stable for arbitrary Keller maps.
For cancellation maps the complete prime-intersection diagram is now a
separate exact computation, not something inferred from a target nilpotency
index.

The cancellation extractor supplies the parallel geometric profile for type
`(m,r)`: the discriminant prime has `(e,f)=(r+1,1)`, different exponent `r`,
completed equation `delta=unit*q^(r+1)`, and one `(r+1)`-cycle; the `mr-1`
geometric primes over the second divisor are unramified.  It also verifies the
affine-plus-boundary degree sums.  Arithmetic grouping of those `mr-1` primes
is retained by the separately computed factorization of the cancellation
parameter polynomial.

At the generic point of the thick component `P=Q=0`, the executable
higher-stratum profile also retains the completed target plane `k(R)[[P,Q]]`,
the critical-normalization plane `k(R)[[Y,Q]]`, and

\[
 P=(Q-Y)Y^m.
\]

It records the two branch contributions separately: `Y=0` has divisor
multiplicity `m`, contact order `mr`, and length contribution `m^2r`, while
`Y=Q` has multiplicity one, contact order `mr`, and contribution `mr`.
Their sum is the target nilpotency index `mr(m+1)`.

The complete boundary extractor refines this with

\[
 K_{m,r}(w)=w^{-r-1}\int_0^w v^r(1-v)^{mr}\,dv,
 \qquad
 L_{m,r}(w)=w^{-r-1}\int_0^w
 \{w(1-w)^m-v(1-v)^m\}^r\,dv.
\]

The exact certificate `Res(K,L) != 0` makes `L` a unit in the reduced finite
`K`-branch algebra.  It follows that every geometric `K`-branch meets the
critical prime with completed ring `k(R)[[Q]]/(Q^m)`, while every two distinct
`K`-branches meet in the common reduced stratum `S={P=Q=0}`.  Higher
intersections remain reduced.  This accounts for every boundary prime and
specialization arrow.  Nonvanishing is proved wherever the cancellation
parameter polynomial is known irreducible—in particular for every `mr<=30`,
the full `m=1` ladder, and the uniform arithmetic criteria.  A separate
endpoint-moment, Schur--Cohn, and rational-certificate argument proves all
four complete columns `r=1,2,3,4`, without irreducibility, and direct
resultants are checked for `1<=m,r<=5`.  Outside proven ranges the resultant
is retained as an explicit per-input certificate.

### Completion audit

The following distinction prevents the formal definition from being mistaken
for a completed computation in every family.

| requested layer | theorem/definition | exact extractor and regression |
|---|---|---|
| intrinsic finite normalization cover and stabilization | complete | generic weighted and cancellation divisorial profiles |
| ordered discriminant boundary and normalized pencil | complete under the stated intrinsic-completeness hypothesis | complete for the weighted seeds in the repository |
| full Fitting divisor, including multiplicities | complete | complete |
| saturated off-diagonal scheme and its `S_2` quotient | complete in characteristic zero | complete; quotient pullback and node transversality are checked |
| conductor map | complete as an intrinsic finite-stratum construction, including the admissible-cover contraction | exact implicit-equation formula for arbitrary plane-curve singularities; simultaneous multicluster exponent `e_i(sum_j e_j-1)` and normalized-Stein comparison checked separately |
| infinity and second-boundary marks | complete on the boundary-clean, two-Fitting-support locus; otherwise a mark is not retained without an independent intrinsic characterization | exact for the weighted clean locus; the quartic zero-cluster chart is checked |
| stacky marked extension | complete at arbitrary simultaneous multicluster collisions | selected root retained on the marked admissible-cover stack; quotient map finite representable etale of degree `N-2`; normalized-Stein contractions, completed root-cover charts, and conductor squares identified |
| coarse affine root-one descent | complete after contraction | the finite birational closure is an isomorphism over the normal marked coarse compactification; the local invariant ring is the universal monic-root incidence and its total-collision fiber `k[T]/(T^mu)` has one geometric point |
| upstairs `(e,f)`, different, DVR, and inertia data | complete as invariant data | generic divisorial layer complete |
| higher intersections and completed local extensions | complete as functorial invariant data | full cancellation prime diagram is exact under `Res(K,L) != 0`; this is proved in the five columns `r=1,2,3,4,5`, the six transverse columns `1<=m<=6`, every fixed-`m` eventual tail, and all other parameter-irreducibility ranges; the residual problem has `m>=7`, `r>=6`, with the first open column reduced exactly to a degree-29 endpoint eliminant |

The coarse affine-mark descent theorem is tracked as `H3` and closes the
former `OP-MARK`; the higher cancellation resultant remains `OP-CR` in
[STATUS.md](../STATUS.md).

The stability inputs agree with standard base-change results: normalization
commutes with smooth base change ([Stacks, Lemma 76.25.2](https://stacks.math.columbia.edu/tag/082F)),
and the Kähler different is `Fitt_0 Omega` and commutes with arbitrary base
change ([Stacks, Section 49.7](https://stacks.math.columbia.edu/tag/0BVV)).
The conductor used here is the annihilator of the normalization quotient, in
the standard sense reviewed for curves by the
[Stacks Project](https://stacks.math.columbia.edu/download/curves.pdf).
These references justify the functorial framework; the displayed seed and
boundary formulas remain repository-specific calculations and are therefore
covered by exact regressions rather than outsourced to a general theorem.

## 7. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_decorated_normalization.py
.venv/bin/python scripts/verify_affine_branch_mark_audit.py
.venv/bin/python scripts/verify_coarse_affine_mark_descent.py
.venv/bin/python scripts/verify_stable_generator_rigidity.py
.venv/bin/python scripts/verify_full_boundary_diagram.py
```

The checker retains Fitting multiplicities, performs both ordered and
unordered saturations with Rabinowitsch elimination, verifies the symmetric
quotient by pullback, checks the node-scheme Jacobian, extracts the finite
conductor map upstairs and downstairs, reconstructs all three quartic
decorations, verifies the residual cusp/node coincidence of `F4b,F4c`, proves
their separation by (4.4) and the boundary chart (3.6), and checks generic
boundary-cover `(e,f)`,
different, completed-DVR, inertia, residue-factor, and degree-sum data for
both weighted and cancellation profiles.

The affine-branch audit separately verifies the unramified root-one stratum,
its disjointness from the ramification divisor, and the bare rerooting
identity (5.3d).  The stable-generator audit retains the former low-pole and
target-orientation identities as independent checks.  Its load-bearing
contribution to the new proof is the higher-zero, nonzero-multiple-root, and
final `z`-coordinate analysis showing that every simple extra root is polar.
