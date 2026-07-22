# The scheme-theoretic decorated-normalization invariant

This note upgrades the reduced Hessian-root support to a scheme-theoretic
stable invariant, computes the off-diagonal self-intersection of the
discriminant normalization, and applies the resulting decoration to the two
previously unresolved split quartic islands.

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

### Stable functoriality, including multiplicities

Suppose two split weighted marked-root maps have the intrinsically ordered
target-boundary pair `(Z_Delta,Z_0)` used in the stable boundary theorem.
Removing `Z_0` identifies the discriminant vertex and its normalization with

\[
 D_H\times G_m,
 \qquad
 \widetilde D_H\times G_m.
\]

A polynomial left--right equivalence identifies these finite normalization
morphisms.  It therefore identifies their relative differential modules and
their zeroth Fitting ideal sheaves.  After adjoining `q` identity variables,

\[
 \Omega_{(\widetilde D_H\times A^q)/(D_H\times A^q)}
 \simeq \operatorname{pr}_1^*Omega_{\widetilde D_H/D_H}.
\]

Fitting ideals commute with this flat base change.  Hence the effective
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

The scheme-theoretic decorated normalization is the tuple

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

where the full center divisor underlying `b_H` is `gcd(H,H')`, including its
multiplicities.  Point labels are inferred only when this divisor has a
uniquely certified support point; the extractor never invents a mark at zero
for an arbitrary polynomial.  For a normalized split weighted seed, the
canonical zero-cluster mark is `r=0`.
Further intrinsically labeled boundary branches may be retained as additional
marks.  Every entry in (3.5) is constructed from the canonical normalization,
its finite map, and the ordered boundary pair, so stable left--right
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

## 5. The moduli map and its generic dimension

Let `A_N` be the normalized admissible seed space

\[
 H(W)=\sum_{j=3}^N h_j(W^j-W^2),
 \qquad H'(1)=-1,
\]

with the degree-drop and weighted exceptional divisors removed.  It has
dimension `N-3`.  On the further open where the discriminant has only
ordinary cusps and nodes and the boundary mark avoids them, (3.5) defines
the rational moduli map

\[
 \mathfrak M_N:\mathcal A_N\dashrightarrow
 \mathcal M^{\mathrm{dec}}_{0,
 (N-2)\,\mathrm{cusps},
 \binom{N-2}{2}\,\mathrm{paired\ nodes},
 1\,\mathrm{boundary},1\,\mathrm{infinity}}.           \tag{5.1}
\]

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

The target means decorated rational curves modulo `PGL_2`; the infinity mark
reduces the acting group to `Aff_1`.  The Fitting divisor alone has `N-2`
points modulo this two-dimensional group and therefore has expected image
dimension `N-4`, agreeing with the one-dimensional degree-five result.  The
canonical boundary mark adds one point, so its configuration already has
expected dimension

\[
 (N-2)+1-2=N-3.                                       \tag{5.2}
\]

This equals `dim A_N`; node pairing and the conductor add structure but no
new source parameters.  In fact the generic fiber is finite.  Suppose the
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
 \boxed{\dim\overline{\operatorname{im}(\mathfrak M_N)}=N-3}
\]

on the generic ordinary locus.  This proves the image dimension; it does not
claim that every abstract decorated configuration occurs.

### The affine root-one stratum and the missing cross-stratum lemma

There is a tempting apparent improvement of the preceding argument.  Over
`Z_0=V(C)`, the normalized incidence does intrinsically distinguish the
simple root `W=1` from the additional simple roots.  This distinction is real,
but it does **not** yet supply another point on the discriminant-normalization
line.

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

It is therefore safe to refine the **full-cover** decoration by adjoining

\[
 \mathcal A_{1,H}\longrightarrow Z_0,                \tag{5.3aa}
\]

where `mathcal A_(1,H)` denotes this distinguished affine component together
with its embedding in the regular-reconstruction open.  This is the intrinsic
version of the “root-one branch.”  It is a genuine stable invariant, but it is
a stratum of the normalized cover, not yet a mark on the rational curve in
(3.5).

The missing step is a cross-stratum identification.  The ramification divisor
in the incidence is cut out by `E=partial_W E=0`.  On the root-one affine
component,

\[
 C=0,qquad W=1,qquad
 \partial_WE=H'(1)-BC=-1.                            \tag{5.3b}
\]

Thus this component is disjoint from the ramification divisor.  In contrast,
the mark `r=0` in (3.5) came from an actual valuation center where the
discriminant divisor meets `Z_0`.  There is presently no functorial map from
the disjoint affine component in (5.3b) to a point of
`widetilde D_H`.  Calling its root value `r=1` uses the model's primitive
element `W`; an abstract isomorphism of finite covers need not transport that
primitive element by the same affine formula as its restriction to the
ramification divisor.

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
This shows exactly where the intrinsic open contains more information than
the present decorated-normalization target, but it does not manufacture a
point on that target.

There is nevertheless a clean conditional classification statement.  Suppose
one proves a **cross-stratum generator-rigidity lemma** saying that an
isomorphism of intrinsic covers with their distinguished affine opens carries
the root-one component to a point of the discriminant normalization through
the same affine coordinate change obtained from the Fitting divisor.  That
change already fixes `0` and infinity, so it is `r -> ar`; compatibility with
the root-one component forces `a=1`.  Equation (5.3) then gives `G=cH`, and
`G'(1)=H'(1)=-1` forces `c=1`.  Consequently such a lemma would imply

\[
 F_H\sim_{\mathrm{stable}}F_G\quad\Longrightarrow\quad H=G. \tag{5.3e}
\]

No such cross-stratum lemma is proved here.  Accordingly `r=1` is not added
to (3.5), and the generic-finiteness bound is not upgraded to generic
injectivity.  Establishing or refuting generator rigidity for the normalized
incidence together with its regular-reconstruction open is the precise next
classification problem.

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

The active computable invariant is organized in three functorial layers.

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

The full-cover extractor refines this with

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
the full `m=1` ladder, and the uniform arithmetic criteria—and is also checked
by direct resultants for `1<=m,r<=5`.  Outside proven ranges the resultant is
retained as an explicit per-input certificate.

### Completion audit

The following distinction prevents the formal definition from being mistaken
for a completed computation in every family.

| requested layer | theorem/definition | exact extractor and regression |
|---|---|---|
| intrinsic finite normalization cover and stabilization | complete | generic weighted and cancellation divisorial profiles |
| ordered discriminant boundary and normalized pencil | complete under the stated intrinsic-completeness hypothesis | complete for the weighted seeds in the repository |
| full Fitting divisor, including multiplicities | complete | complete |
| saturated off-diagonal scheme and its `S_2` quotient | complete in characteristic zero | complete; quotient pullback and node transversality are checked |
| conductor map | complete as an intrinsic finite-stratum construction | exact implicit-equation formula for arbitrary plane-curve singularities, checked against the ordinary factorization |
| infinity and second-boundary marks | complete | exact for the weighted seeds; the quartic zero-cluster chart is checked |
| distinguished affine root-one stratum | complete on the exact-double-zero boundary-clean locus | exact cover-stratum and rerooting audit; no cross-stratum point on the discriminant normalization is claimed |
| upstairs `(e,f)`, different, DVR, and inertia data | complete as invariant data | generic divisorial layer complete |
| higher intersections and completed local extensions | complete as functorial invariant data | full cancellation prime diagram is exact under `Res(K,L) != 0`, proved for every `mr<=30`, the `m=1` ladder, and all uniform parameter-irreducibility cases; all-parameter nonvanishing remains open |

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

The affine-branch audit separately verifies the unique unramified root-one
stratum, its disjointness from the ramification divisor, and the exact
rerooting identity (5.3d).  Its scope line records that cross-stratum
generator rigidity remains open.
