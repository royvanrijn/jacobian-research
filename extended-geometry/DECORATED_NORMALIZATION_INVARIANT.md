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
specialization arrow.  Nonvanishing is proved for the full `m=1` ladder and
checked exactly for `1 <= m,r <= 5`; outside proven ranges the resultant is
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
| upstairs `(e,f)`, different, DVR, and inertia data | complete as invariant data | generic divisorial layer complete |
| higher intersections and completed local extensions | complete as functorial invariant data | full cancellation prime diagram is exact under `Res(K,L) != 0`, proved for the `m=1` ladder and certified on `1 <= m,r <= 5`; all-parameter nonvanishing remains open |

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
