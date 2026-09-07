# Global logarithmic second-Chern-character budget

> **Status.**  This note proves the global Chern-character identity for the
> logarithmic differential of a generically etale morphism of smooth
> projective SNC surface pairs, proves its exact transformation law under
> boundary blowups, and separates the signed global class from any effective
> zero-dimensional Fitting defect.  It applies the identity to the currently
> compiled F2 `(75,125)` source and target lower-bound graphs.  On that fixed
> partial model the squarefree and double-row budgets are respectively
> `(7*d-8)/2` and `(7*d-13)/2`; subtracting the exact cyclic root contribution
> `27` gives `-10` at the squarefree degree floor `d=6` and `17/2` at the
> double-row floor `d=12`.  These are **not exclusions**: affine purity forces
> a new component and raises the source floors to `28/49`, but its target
> curve, common-model centers, and Chern module are uncompiled; the final
> geometric degree is also unknown, and a residual becomes a nonnegative
> integer only after an exact filtration realizes it by a finite-length sheaf.

The rational identities, fan blowup types, source intersection calculations,
and F2 specializations are replayed by
[`verify_logarithmic_ch2_budget.py`](../scripts/verify_logarithmic_ch2_budget.py).
The reusable arithmetic is in
[`log_node_profiles.py`](../jcsearch/log_node_profiles.py).

## 1. The perfect boundary complex

Let

\[
 f:(X,D_X)\longrightarrow(Y,D_Y)                 \tag{1.1}
\]

be a proper generically finite morphism of degree `d` between smooth complex
projective surfaces with reduced SNC boundaries.  Suppose the logarithmic
differential is an isomorphism on the selected open surface.  Put

\[
 E=f^*\Omega_Y^1(\log D_Y),\qquad
 F=\Omega_X^1(\log D_X)                           \tag{1.2}
\]

and use the K-theory convention

\[
 \boxed{[\mathcal K_f]=[F]-[E].}                  \tag{1.3}
\]

Thus `K_f` is represented by `E->F`, with `E` in degree `-1`.  It is exact
away from the boundary.  The localized Chern character is therefore defined
with boundary support, and its pushforward to `X` is the ordinary Chern
character of (1.3).  See the Stacks Project,
[*Localized Chern classes*](https://stacks.math.columbia.edu/tag/0FB0).

If the logarithmic differential is generically invertible, it is injective:
locally, multiplication by the adjugate reduces a kernel vector to torsion by
the nonzero determinant in the regular local domain.  Hence

\[
 \mathcal K_f\simeq\mathcal T_f^{\log}
 =\operatorname{coker}(E\longrightarrow F).       \tag{1.4}
\]

The universal purity theorem in
[`UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md`](UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md)
then makes this a pure one-dimensional Cohen--Macaulay sheaf.  This does not
make its codimension-two Chern character positive.

## 2. Global budget theorem

Write

\[
 L_X=K_X+D_X,\qquad L_Y=K_Y+D_Y,                  \tag{2.1}
\]

and put

\[
 e_X=\int_Xc_2\!\left(\Omega_X^1(\log D_X)\right),
 \qquad
 e_Y=\int_Yc_2\!\left(\Omega_Y^1(\log D_Y)\right). \tag{2.2}
\]

### Theorem 2.1 -- logarithmic Chern-character budget

Under the hypotheses of Section 1,

\[
 \boxed{
 \deg\operatorname{ch}_2(\mathcal K_f)
 =\frac12\left(
 L_X^2-dL_Y^2-2e_X+2de_Y
 \right).}                                      \tag{2.3}
\]

If both complements are `A^2`, then `e_X=e_Y=1`, so

\[
 \boxed{
 B_f:=\deg\operatorname{ch}_2(\mathcal K_f)
 =\frac12\left(L_X^2-dL_Y^2+2(d-1)\right).}      \tag{2.4}
\]

#### Proof

For a rank-two bundle `V`,

\[
 \operatorname{ch}_2(V)
 =\frac12\left(c_1(V)^2-2c_2(V)\right).          \tag{2.5}
\]

The first Chern classes in (1.2) are `f^*L_Y` and `L_X`.  Proper
degree-`d` pullback gives

\[
 \int_X(f^*L_Y)^2=dL_Y^2,\qquad
 \int_Xf^*c_2\!\left(\Omega_Y^1(\log D_Y)\right)=de_Y. \tag{2.6}
\]

Subtract the two instances of (2.5).  This proves (2.3).  For an SNC
compactification, the top logarithmic Chern class equals the compactly
supported Euler characteristic of the complement.  Both complements in
(2.4) are `A^2`, whose Euler characteristic is one.  \(\square\)

This theorem concerns `ch_2`.  Calling (2.4) merely a `c_2` formula is
ambiguous: for the virtual rank-zero class (1.3),

\[
 c_2(\mathcal K_f)
 =\frac12c_1(\mathcal K_f)^2
  -\operatorname{ch}_2(\mathcal K_f).            \tag{2.7}
\]

In ramification coordinates this separate virtual Chern class is

\[
 \deg c_2(\mathcal K_f)
 =-f^*L_Y\mathbin{\cdot}R_{\log}+e_X-de_Y,       \tag{2.8}
\]

and for two `A^2` complements it is
`-f^*L_Y dot R_log-(d-1)`.  It should not be conflated with the additive
point-class ledger below.

The point class in a K-theory exact sequence is most directly additive in
`ch_2`, which is why (2.4) is the ledger used here.

## 3. Ramification form

The logarithmic determinant is a section of

\[
 \det(F)\otimes\det(E)^{-1}
 =O_X(L_X-f^*L_Y).                                \tag{3.1}
\]

Let

\[
 R_{\log}=L_X-f^*L_Y.                             \tag{3.2}
\]

Expanding `L_X=f^*L_Y+R_log` in (2.3) gives

\[
 \boxed{
 \deg\operatorname{ch}_2(\mathcal K_f)
 =f^*L_Y\mathbin{\cdot}R_{\log}
  +\frac12R_{\log}^2-e_X+de_Y.}                 \tag{3.3}
\]

For `A^2` pairs this is

\[
 f^*L_Y\mathbin{\cdot}R_{\log}
 +\frac12R_{\log}^2+d-1.                        \tag{3.4}
\]

Thus the global budget is computable from the same intersection matrix,
canonical vector, and target pole vector already used by the intrinsic
boundary compiler.  Formula (2.4) is preferable while the complete
ramification vector is not yet available; (3.4) becomes the independent
cross-check once it is.

## 4. Exact boundary-blowup law

Let `g:X'->X` be the blowup of a boundary point and let `D_X'` be the reduced
boundary with the same open complement.

At an SNC node, the reduced total transform loses one exceptional copy while
`K_X'` gains it, so

\[
 K_{X'}+D_X'=g^*(K_X+D_X).                       \tag{4.1}
\]

The blowup is log crepant and `L_X^2` is unchanged.

At a smooth boundary point,

\[
 K_{X'}+D_X'=g^*(K_X+D_X)+E,
 \qquad E^2=-1,                                  \tag{4.2}
\]

and therefore

\[
 L_{X'}^2=L_X^2-1.                               \tag{4.3}
\]

Consequently, relative to a fixed compiled model, `s_X` additional smooth
source-boundary blowups and `s_Y` additional smooth target-boundary blowups
change (2.4) by

\[
 \boxed{\Delta B_f=\frac{d s_Y-s_X}{2}.}         \tag{4.4}
\]

Any number of source or target boundary-node blowups contributes zero to
(4.4).

This is a covariance law, not an invariance statement.  Local packets must
always be evaluated on the same common source and target model as the global
budget.  In particular, a cyclic `R/(w^m)` presentation computed over a
smooth target-boundary point cannot be carried unchanged into a ledger after
that target point has been blown up and the source fan has been refined: it
may become log-etale or join a larger total-transform determinant packet.

### Example 4.1 -- raw `ch_2` is signed

Blow up one smooth point of the line at infinity in the identity pair
`(P^2,L_infinity)`.  The open map is still the identity on `A^2`, and `d=1`,
but

\[
 L_X^2=3,\qquad L_Y^2=4.
\]

Formula (2.4) gives

\[
 B_f=-\frac12.                                    \tag{4.5}
\]

Therefore no positivity theorem can apply to the raw global
`ch_2(K_f)` or to arbitrary divisorial summands.

## 5. Cyclic divisorial packets

Suppose on an effective Cartier determinant packet `i:D->X` that `Fitt_1`
is the unit ideal and

\[
 \mathcal T_f^{\log}|_D=i_*M
\]

for an invertible `O_D`-module `M`.  Let `K` be the kernel line of the
restricted logarithmic differential.  The kernel-line theorem gives

\[
 M=K\otimes O_D(D),                               \tag{5.1}
\]

and Grothendieck--Riemann--Roch gives

\[
 \boxed{
 \operatorname{ch}_2(i_*M)
 =\deg_D(K)+\frac12D^2.}                         \tag{5.2}
\]

The `D^2/2` term includes both component self-intersections and thickened
node matching.  For

\[
 D=\sum_i m_iC_i
\]

on an SNC curve without triple points,

\[
 \frac12D^2
 =\frac12\sum_i m_i^2C_i^2
  +\sum_{C_i\cap C_j\ne\varnothing}m_im_j.       \tag{5.3}
\]

The last sum is the length of the canonical branchwise quotients

\[
 0\to R/(u^av^b)
 \to R/(u^a)\oplus R/(v^b)
 \to R/(u^a,v^b)\to0.                            \tag{5.4}
\]

Neither the raw node length nor the component term is separately invariant;
their sum (5.3) is.

## 6. Normalization and `Fitt_1` signs

Let `T` be a torsion-free rank-one module on a reduced determinant curve and
let

\[
 \widetilde T=\nu_*(\nu^*T/\text{torsion}).
\]

Then

\[
 0\longrightarrow T\longrightarrow\widetilde T
 \longrightarrow C_T\longrightarrow0            \tag{6.1}
\]

gives

\[
 [T]=[\widetilde T]-[C_T],\qquad
 \operatorname{ch}_2(T)
 =\operatorname{ch}_2(\widetilde T)-\ell(C_T).   \tag{6.2}
\]

For the glued and split node modules,

\[
 0\to R/(uv)
 \to R/(u)\oplus R/(v)
 \to R/(u,v)\to0,                                \tag{6.3}
\]

so the split module exceeds the glued one by one point class.  Their
`Fitt_1` ideals, respectively `R` and `(u,v)`, detect the different local
module structures.  But `Fitt_1` alone supplies only the locus: the exact
presentation or normalization sequence is needed to determine the length,
the comparison module, and hence the sign in the global K-class.

The same warning applies to nonreduced determinant curves.  One must first
record their generic Smith/thickening filtration and only then normalize the
reduced branch data.

## 7. Conditional effective-residual criterion

The useful exclusion statement is deliberately conditional on an exact
K-theory decomposition.

### Proposition 7.1 -- certified point-residual test

Assume a fixed common boundary model admits a proved exact filtration whose
one-dimensional quotients are cyclic packets `i_(alpha)*M_alpha` as in
Section 5 and whose remaining quotient is a finite-length sheaf `Z`.  Then

\[
 [\mathcal T_f^{\log}]
 =\sum_\alpha[i_{\alpha *}M_\alpha]+[Z]           \tag{7.1}
\]

and

\[
 \boxed{
 \ell(Z)=B_f-
 \sum_\alpha\left(
 \deg_{D_\alpha}(K_\alpha)+\frac12D_\alpha^2
 \right).}                                      \tag{7.2}
\]

In particular, the right side of (7.2) must be a nonnegative integer.

#### Proof

Apply additivity of the Chern character to the proved filtration, use (5.2)
for every one-dimensional quotient, and use
`ch_2(Z)=length(Z)` for the finite quotient.  \(\square\)

Without the filtration, the same subtraction produces only a **virtual
residual** in `CH_0(X) tensor Q`.  It need not be effective or integral.
This is the exact missing hypothesis behind a positivity-based kill
criterion.

### Corollary 7.2 -- parity gate

Under Proposition 7.1, assume the log squares and all kernel-line degrees are
integral.  Reducing twice (7.2) modulo two gives the degree-independent
necessary congruence

\[
 \boxed{
 L_X^2-dL_Y^2-\sum_\alpha D_\alpha^2\equiv0\pmod2.} \tag{7.3}
\]

Indeed, the top-log-Chern terms, twice the kernel degrees, and twice the
finite length are all even.  Thus parity can be checked before the kernel
degrees or point lengths are known.  Failure of (7.3) does not by itself
exclude a partial model: it proves that at least one divisorial packet or
model change is still missing from that ledger.

There are two distinct sources of signed point corrections:

1. a contracted cyclic packet has
   `deg(K)=-e_gamma<=0`, where `e_gamma` is its kernel Gauss degree;
2. a noncyclic nodal module can differ from a chosen cyclic/glued reference
   by positive or negative point classes according to the exact sequence
   used.

A global theorem must retain both.  A determinant-only ledger sees neither.

## 8. The current F2 target square

For `(P^2,L_infinity)`, the initial logarithmic square is

\[
 (K_{P^2}+L_\infty)^2=(-2H)^2=4.                 \tag{8.1}
\]

The minimal target fan extracting `(5,2)` has insertion order

\[
 (1,1),(2,1),(3,1),(5,2).                       \tag{8.2}
\]

The first blowup is at a smooth boundary point and the other three are node
blowups.  The carrier ray `(5,36)` is then extracted above the smooth point
`h=125/729` by

\[
 (1,1),(1,2),\ldots,(1,8),(2,15),(3,22),(4,29),(5,36). \tag{8.3}
\]

The first eight are successive smooth-boundary blowups and the final four
are node blowups.  Hence the currently compiled common target has

\[
 \boxed{L_Y^2=4-1-8=-5.}                         \tag{8.4}
\]

This uses both target clusters.  Using only the four-blowup `(5,2)` target
would incorrectly give `L_Y^2=3` and mix it with source packets already
refined against the `(5,36)` fan.

## 9. The current F2 source squares

Solving adjunction on the refined source intersection matrices and squaring
the coefficient vector of `K_X+D_X` gives

\[
 \boxed{
 L_X^2=-6\quad\text{(squarefree)},\qquad
 L_X^2=-11\quad\text{(double row)}.}             \tag{9.1}
\]

Equivalently, the terminal-resolved `19/31` graphs have squares `-4/-10`;
the two squarefree smooth spectator blowups lower the first by two, while the
first smooth blowup of the double-row alpha chain lowers the second by one.
All principal-arm common refinements and the remaining alpha refinements are
boundary-node blowups and are log crepant.

The affine-purity frontier proves that these `27/48` partial-model graphs
must acquire at least one new affine-branch component, giving global floors
`28/49`.  The outgoing terminal tail is unimodular and adds no blowup or
Chern term, but the target equation, proximity chain, and Chern module of the
new component, together with any other genuinely global centers, remain
absent from this fixed model.
The target equation is now restricted to 24 singular normalization charts
`(3k,5k)`, `1<=k<=24`, by the affine target-curve theorem; this does not yet
determine which chart or its source pullback.
For `k=1`, the generic target conductor is exact—four affine nodes and the
delta-2 infinity cusp—and its generic boundary module is now computed by the
[`affine-row Chern theorem`](F2_AFFINE_K1_LOG_CH2.md).  For data
`(e,f,E^2=-n)` and carrier contact count `b`, its cyclic contribution is
`e*f*(b-7)-e^2*n/2`.  The leading residue determines `b=0` off the special
carrier point; on it, the carrier-normalized jet gives the exact rule
`b=min(ord_u(w|_C),8)`, reducing the contact to at most seven coefficient
equations.  It does not determine `e`: the terminal
neighborhood is already resolved and cannot extract the affine divisor.
Both target-contact cases still need their nonterminal source attachment,
self-intersections, and the finite point set where `Fitt_1` is nonunit before
the effective filtration is complete.

The
[`carrier-jet factorization theorem`](F2_AFFINE_K1_CARRIER_JET_FACTORIZATION.md)
now gives the exact fixed-coordinate meaning of those seven equations.  For
fixed `P0,Q0,Gamma`, four rows reconstruct `a,b,c,d` and three weighted rows
are compatibility residuals.  But with the normalization parameters free,
the raw seven-jet map has Jacobian `3*Res(p',q')`; it is étale and dominant
on the immersed locus.  Hence the carrier jet alone supplies no generic
`k=1` obstruction and cannot certify `b`: the next compiler must first fix
or globally relate the target normalization through the affine node fibers,
the completed source pullback, or equivalent global data.  The severe
`E_6+A_1` subfamily is exceptional: its fixed seven-jet locus remains
codimension three even after all target transports, so the actual carrier
center vector gives a finite three-equation cusp test.  Its `E_8` endpoint
is the further hyperplane section `J_1=0` and gives four scale-free raw
carrier equations, hence a codimension-four **maximal-contact** endpoint
test.  Earlier departure only imposes the initial center matches.

<!-- status-consumer: PF2K1JF1 7bc57f390f0531b5 -->

The endpoint monodromy is also rigid.  For meridian type `2+2+1+1`, the
`E_8` torus-knot complement has exactly one transitive degree-six action up
to conjugacy, with image `A_5`.  Its preferred longitude preserves the two
ramified transposition orbits, so its localized-Chern filtration must contain
two distinct `(e,f)=(2,1)` packets; one `(2,2)` packet is impossible.

If `N=n_1+n_2` is the sum of their self-intersection magnitudes, the exact
two-row necessary inequality is `28+4N-8b-s_X>=0`.  Hence maximal contact
`b=8` forces `4N>=36+s_X`, while the minimal `N=2` packet forces `b<=4`.
Eliminating the E8 target transports shows that the first invariant
raw-center equation occurs only at `b>=5`; thus the minimal Chern survivor
lies strictly before the present carrier compatibility gate.

The complete `A_5` coset atlas enlarges this to fixed-sheet degrees
`6,10,15,30`, with `r=2,4,6,14` separate `(2,1)` rows.  Its uniform
squarefree doubled remainder is `7d-62+4N-4r(b-6)-s_X`, and the double-row
version replaces `62` by `67`.  At minimal `N=r`, maximal contact makes the
degree-six and degree-ten squarefree budgets negative, but leaves degrees
15 and 30 positive.  A global bound on negativity or a compatibility
condition with the other boundary packets is therefore indispensable.

<!-- status-consumer: PF2K1E8M1 bbb282c6bcfa62fc -->

The complete simple-inertia orbifold atlas refines `r` into two different
quantities: `q`, the number of source divisors, and `R`, their total residue
degree.  Its order-four central action produces `(2,f)` rows for `f=1,2,4`
and fixed-sheet degrees `6,10,12,15,20,24,30,40,60,120`.  The squarefree
doubled budget is `7d-62+4N-4R(b-6)-s_X`, while the component floor is
`27+q`.  This distinction is precisely the normalization information absent
from a determinant-only decomposition.

<!-- status-consumer: PF2K1E8O1 4251750ed4e43c89 -->

The complete-chain calculation now performs the missing stable subtraction.
For one rational one-puncture `k=1` affine component, the full Cartier square
and conormal degree leave exactly `u-1` point units.  The E8 cusp needs at
least `2R`, and every simple-inertia atlas row has `2R>u-1`.  Consequently
raw self-intersection negativity cannot repair the packet: the only possible
escape is a negative normalization/`Fitt_1` class at the unresolved global
attachment.

<!-- status-consumer: PF2K1CB1 5cc386dba344a867 -->

The isolated-sign ambiguity is now removed.  A generically cyclic `2 x 2`
logarithmic cokernel contains its determinant/kernel-line cyclic module with
an effective finite quotient whose length dominates the `Fitt_1` colength.
Therefore the negative correction required by every simple-inertia E8 row
does not exist; the whole one-component simple-inertia branch is excluded.

<!-- status-consumer: LCSP1 8658eebeb1d65671 -->

The target-side finite correction is nevertheless explicit: the
[`implicit-conductor theorem`](F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md)
gives a twelve-support quintic whose gradient pulls back to the degree-eight
nodal conductor times `(q',-p')`.

<!-- status-consumer: PF2K1I1 a7582c1e36140840 -->

The
[`fixed-coordinate Keller-pullback theorem`](F2_AFFINE_K1_KELLER_PULLBACK.md)
then restores the F2 target chart and closes the affine local correction:
the pullback is squarefree and its affine normalization defect is the
fiberwise sum of ordinary nodes of length one over four explicit target
values.  The
unknown point terms in this budget are therefore boundary-local; the only
remaining affine input is the four finite fiber counts.

<!-- status-consumer: PF2K1PB1 6f837229017243c4 -->

The
[`all-stratum conductor-conservation theorem`](F2_AFFINE_TARGET_K1_CONDUCTOR_CONSERVATION.md)
extends this affine control across the complete `k=1` discriminant.  The
target affine normalization quotient always has length four and its
conductor divisor degree eight; after Keller pullback the respective bounds
are `4(d-1)` and `8(d-1)`.  Only the weighted singular-fiber counts enter the
remaining point budget.

<!-- status-consumer: PF2K1CC1 f152c82ef2d54c32 -->

The
[`tame-node packet theorem`](F2_AFFINE_K1_TAME_NODE_PACKET.md) shows that
none of this ordinary conductor degree is automatically a logarithmic
point correction.  An fs tame Kummer toroidal packet over the resolved node
has zero logarithmic cokernel and localized `ch_2`, including the collided
cyclic packet after its toric resolution.  Any nonzero boundary point term
in the budget must therefore measure a certified failure of that toroidal
chart model.  For a general SNC monomial-with-unit pullback, full-rank
exponent data still force zero cokernel, and a rank-one packet can have
singular determinant support only if two explicit logarithmic-unit first
jets vanish.

<!-- status-consumer: PF2K1TN1 521fb57f7e6abc1f -->

The
[`affine strict-log-étale resolution theorem`](AFFINE_KELLER_STRICT_LOG_ETALE_RESOLUTION.md)
removes the whole affine conductor packet from the relative point ledger.
Embedded resolutions of nodes, cusps, tacnodes, and higher singularities
pull back strictly étale, so their relative logarithmic cokernels vanish.
Ordinary delta-weighted fiber counts remain boundary-escape constraints;
only compactification-boundary matrices enter the filtration (7.1).

<!-- status-consumer: PAER1 60eb24b2232d159e -->

## 10. Provisional F2 budgets

Substituting (8.4) and (9.1) in (2.4) gives

\[
 \boxed{
 B_{\rm sq}(d)=\frac{7d-8}{2},\qquad
 B_{\rm dbl}(d)=\frac{7d-13}{2}.}                \tag{10.1}
\]

The extraction-root packet has

\[
 D_{\rm root}=3E+18L,\qquad
 E^2=-6,\quad L^2=0,\quad E\cdot L=1,            \tag{10.2}
\]

so `D_root^2/2=27`.  The fixed tangential target covector trivializes its
kernel line on the complete nonreduced packet, making its actual contribution

\[
 \boxed{\operatorname{ch}_2(T_{\rm root}^{\log})=27.} \tag{10.3}
\]

The generic multiplicity-three endpoint on this upstream determinant chain
is part of (10.2); the self-intersection term `9E^2/2=-27` and the thickened
root matching term `3*18=54` are already combined in (10.3).  It must not be
booked again as a separate point class.  More generally, every endpoint term
must be transported to the fixed common target model before it is summed.

Subtracting (10.3) gives the current virtual remainders

\[
 \boxed{
 R_{\rm sq}(d)=\frac{7d-62}{2},\qquad
 R_{\rm dbl}(d)=\frac{7d-67}{2}.}                \tag{10.4}
\]

At the known degree floors,

\[
 R_{\rm sq}(6)=-10,
 \qquad
 R_{\rm dbl}(12)=\frac{17}{2}.                  \tag{10.5}
\]

The first number would be an immediate exclusion if Proposition 7.1 were
available on a complete model with the root as the only divisorial quotient.
It is not: purity already forces an additional affine ramification row.  The
half-integral second number itself proves that an unaccounted divisorial
half-intersection term remains before any residual can be interpreted as a
point length.

More generally, on the displayed partial model the root divisor has even
square `54`, so Corollary 7.2 reads

\[
 R_{\rm sq}(d)\in\mathbb Z\Longleftrightarrow d\equiv0\pmod2,
 \qquad
 R_{\rm dbl}(d)\in\mathbb Z\Longleftrightarrow d\equiv1\pmod2. \tag{10.6}
\]

This is a compiler diagnostic, not a restriction on the final degree: a
missing smooth center changes the global parity input, and a missing
divisorial packet contributes its own `D^2` term.

If completion adds `s_X,s_Y` smooth centers, (10.4) becomes

\[
 R_{\rm sq}^{\rm final}(d)
 =\frac{7d-62+d s_Y-s_X}{2},\qquad
 R_{\rm dbl}^{\rm final}(d)
 =\frac{7d-67+d s_Y-s_X}{2},                    \tag{10.7}
\]

before subtracting any new divisorial packet.

## 11. What the next compiler must prove

The next useful computation is now sharply specified.

1. **One fixed common model.**  Supply the complete fixed-coordinate source
   pair, evaluate the four affine node fibers, factor the quintic at the
   unresolved boundary, and complete its source/target smooth-boundary
   centers before comparing local terms.  In particular, fix or globally
   constrain `P0,Q0,Gamma` before interpreting the three normalized carrier
   residuals as obstructions; the free raw seven-jet map is dominant.  The
   seven-center fan ends exactly at raw-parameter saturation, and its first
   normalization-invariant equation lies at jet eight, outside the present
   contact test.
2. **Two global evaluations.**  Compute (2.4) from the adjunction matrices
   and independently compute (3.4) from the complete `R_log` vector.
3. **Generic divisorial filtration.**  For every determinant component,
   record both generic Smith exponents, the component multiplicity, and the
   kernel-line degree or a fixed target covector trivializing it.
4. **Node/cusp exact sequences.**  Record the full `2 x 2` matrix, `Fitt_1`,
   the normalized branch module, and the signed finite quotient at every
   point.  On the generic `k=1` nonimmersion divisor the exact test chart is
   `A_2+3A_1`: its conductor split is `2+6` and its raw carrier Jacobian has
   corank one.  The rank loss is not the localized point length.  Once a
   boundary preimage is located, the source incidence decides the packet.
   At a smooth boundary point the isolated-Fitting lower exponent is
   `2q_p-1`; at an SNC boundary node it is `2q_p`.  The exact ordinary-cusp
   fold realizes point length one, while the minimal two-boundary theorem
   realizes `2q_p`.  For a complete residue-degree-`f` cusp fiber with `h`
   points, `c` of them boundary nodes, the minimal ledger is `2f-h+c`.
5. **Effectivity certificate.**  Construct the actual filtration (7.1).
   Merely subtracting rational numbers is not enough.
6. **Kill test.**  Compare its certified finite length with (7.2).  Negative
   or nonintegral output then gives a proof, not a heuristic.

<!-- status-consumer: PCJDP1 d4c16bb71dfc6b80 -->
<!-- status-consumer: LUAF1 b0279670ffbd3fa5 -->
<!-- status-consumer: LCAD1 7b9c15d3dfae0337 -->

This is degree-independent.  The degree-specific F2 row supplies an unusually
strong regression because one positive cyclic packet is already exact and
nearly every other exposed local fan is log-etale.

## 12. Claim boundary

Theorem 2.1, the blowup law (4.1)--(4.4), Proposition 7.1, and Corollary 7.2
are general exact statements.  The F2 values (8.4)--(10.6) are exact for the
presently compiled lower-bound model.  The generic `k=1` divisorial Chern
module is now known, but the data do not decide that chart, its leading
residue and carrier contact, the nonterminal source attachment, the source
self-intersections, the point-supported `Fitt_1` corrections, the final
geometric degree, or a complete `R_log` vector.  Therefore this note neither
excludes `(75,125)` nor proves `JC(2)`.

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

<!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->

<!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

On the generic cusp face, exact minimal local packets give boundary subledger
`2f-h+c`, where `h` is the number of boundary points over the cusp and `c`
counts those that are boundary nodes.  It ranges from `f` for unramified
smooth folds to `2f` for node-saturated SNC attachments.  The remaining
virtual residual drops by that incidence-sensitive amount; parity is
unchanged.

<!-- status-consumer: PF2K1L1 5221f5659fc19729 -->

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_logarithmic_ch2_budget.py
```

The checker verifies the bundle and ramification formulas over exact
rationals, the complete boundary-blowup transformation law, the signed
smooth-blowup counterexample, both target fan blowup classifications, both
source log-canonical squares, the root charge and kernel trivialization, and
the provisional F2 formulas and parity gate (10.1)--(10.6).

## References

- The Stacks Project,
  [*The Chern character and tensor products*](https://stacks.math.columbia.edu/tag/02UM).
- The Stacks Project,
  [*Localized Chern classes*](https://stacks.math.columbia.edu/tag/0FB0).
- P. Aluffi,
  [*Chern classes of free hypersurface arrangements*](https://arxiv.org/abs/1201.5396),
  including the SNC logarithmic-Chern/CSM-complement identity.
