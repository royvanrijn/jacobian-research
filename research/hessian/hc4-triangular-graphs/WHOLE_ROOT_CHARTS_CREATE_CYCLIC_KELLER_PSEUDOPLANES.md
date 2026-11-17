# Whole root charts create actual cyclic Keller pseudo-planes

18 November 2026. This is an actual realization reduction from the
proved common-coordinate and whole-boundary charts. It is not an
inverse theorem.

If an added unramified boundary component has multiplicity r>1 in its
common-coordinate fibre, its whole plane root chart is itself the source
of a NEW polynomial Keller map. The new generic degree is multiplied by r,
and a cyclic generic deck subgroup of order r has a specified actual
source-image quotient: an affine pseudo-plane with just one multiple
fibre and no other reducible fibre.

The full original finite normalization and every original prime remain
in the construction. The pseudo-plane quotient is a selected etale open,
not the full finite normalization of its intermediate field.

## 1. Actual setup and statement

Let
```
K=(P,Q):X=A2_(x,y)->Y=A2,
Jac_(x,y)(P,Q)=1,
B=C[P,Q], F=C(P,Q), L=C(x,y), d=[L:F].
```
Rescaling a target coordinate achieves the displayed Jacobian without
changing B, F or any normalization.

Let Z be the FULL finite normalization of B in L and U its entire etale
locus. Use the proved
[common source pencil](EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md)
and
[whole boundary charts](COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md).
Fix one added unramified component D. Its whole component open and fibre are
```
M=X union D subset U subset Z,
div_M(y-lambda)=L_lambda+rD,
L_lambda={y=lambda} subset X,   r>=1.
```
M is smooth affine, D and L_lambda are disjoint closed A1 curves, and
Omega=dP wedge dQ=dx wedge dy is a regular exact unit form on M.

Put
```
L'=L(z),   z^r=y-lambda.
```
Let Z' be the FULL finite normalization of B in L'. Over M it is the
normalized root cover Mtilde, with source inverse image A2_(x,z).
Its zero fibre has disjoint components
```
Ltilde, D_1,...,D_r,
```
where each D_j maps isomorphically onto D. Define
```
V=M minus L_lambda,
Vtilde=Mtilde minus Ltilde,
W_j=Vtilde minus union_(k!=j) D_k.
```
Every W_j is the entire plane A2_(t_j,z), with
```
x=z^(1-r)*t_j+h_j(z),   y=lambda+z^r,                    (1)
```
where h_j is finite Laurent. Its nonzero principal part may be chosen
to have no z^-r term and largest pole p>=r+1.

Then:

* K_j:W_j=A2->Y, obtained by restricting Z'->Y, is a polynomial Keller
  map with Jacobian r and generic degree rd. A constant target rescaling
  makes its Jacobian one.
* The extension L'/L is cyclic of degree r, with group mu_r fixing x and
  sending z to zeta*z. It is a generic deck SUBGROUP for K_j; the full
  extension L'/F is not asserted Galois.
* Vtilde->V is a finite etale cyclic cover of degree r, and
  Vtilde is the union of the r translates of W_j. The map W_j->V is
  SURJECTIVE ETALE. Its generic field extension is the stated cyclic
  degree-r extension, although W_j->V is not finite for r>1.
* The morphism c=y-lambda:V->A1 is surjective, with scheme fibre rD
  over0 and exactly one reduced A1 over every nonzero value.
* V is smooth affine, Pic(V)=Z/rZ, pi1_top(V)=Z/rZ and chi(V)=1.
  The form Omega remains exact and nowhere zero on V.

Thus for r>1 any such actual Keller boundary yields an ACTUAL polynomial
plane pseudo-cover of the specified affine pseudo-plane, followed by an
ACTUAL etale morphism V->Y. This is more than a cyclic extension of an
abstract rational field.

For r=1 the same construction gives another plane source W=M minus
L_lambda with the SAME source field L, target ring B, generic degree d
and full normalization Z. It exchanges which one of the two fibre
components is marked as source. It does not manufacture a nontrivial
generic deck subgroup.

## 2. Keep the full normalizations before choosing the opens

The polynomial T^r-(y-lambda) is irreducible over L, by its valuation one
at the original source line L_lambda. Hence [L':L]=r, including r=1.

Let S'=O(Z') be the integral closure of B in L'. It is also the integral
closure of O(Z) in L': integrality is transitive in one direction, and
an element integral over B is integral over any intermediate ring in the
other. Therefore Z'->Z is finite and
```
Z'/mu_r=Z.
```
Over the open M, y is regular and z is integral. Integral closure commutes
with localization, so the inverse image of M in Z' is EXACTLY Mtilde,
the normalized cyclic cover used in the whole-chart theorem.

Globally z need not be integral over O(Z), because y can have poles at
original omitted ramified primes. We do not replace Z' by an unjustified
presentation O(Z)[z]. We normalize the complete finite field extension
first. All original primes of Z and all their primes above them in Z'
are retained. Only the explicitly named opens are then used.

The normalized cover is etale at every D_j. Locally the equation is
y-lambda=u*f^r and normalization adjoins w=z/f with w^r=u a unit-root
extension. Over the source line it is the smooth map
```
(x,z) |-> (x,lambda+z^r).
```
Thus Vtilde is the complete inverse image of V and Vtilde->V is finite
etale. The raw root equation is not substituted for its normalization.

## 3. The new plane maps and their exact generic degrees

Since P,Q are regular on M, their pullbacks are regular on the entire
W_j=A2_(t_j,z), hence are polynomials there. Formula(1) gives
```
dx wedge dy=r*dt_j wedge dz.
```
Consequently
```
Jac_(t_j,z)(P,Q)=r!=0.
```
This proves that K_j is an ACTUAL polynomial Keller map on a WHOLE plane,
not a rational Darboux expression or a generic boundary germ.

The function field of W_j is L'=L(z); the target field is still F.
Therefore
```
[C(t_j,z):F]=[L':L]*[L:F]=r*d.                            (2)
```
The r root rotations fix L and hence P,Q. They supply mu_r as a generic
deck subgroup. In field coordinates on W_j their action is
```
z |-> zeta*z,
t_j |-> zeta^-1*[t_j+z^(r-1)*(h_j(z)-h_j(zeta*z))].      (3)
```
This need not preserve the selected plane W_j: its polar terms can move
that plane to a different whole chart. No polynomial source automorphism
or target Ga action is inferred from(3). Nor does a proper cyclic subgroup
make L'/F a Galois Keller extension.

There is an additional exact necessary restriction on the normalized
principal part h_j(z)=sum a_q*z^-q, q>=r+1:
```
gcd(r, {q : a_q!=0})=1.                                 (3a)
```
If a nontrivial root zeta of order dividing this gcd fixed h_j, (3)
would become the polynomial automorphism
(t_j,z) |-> (zeta^-1*t_j,zeta*z) of the ENTIRE W_j. It fixes the origin,
contradicting etaleness of K_j and the nonidentity deck action. Equivalently
it would stabilize the selected zero line D_j, although the root group
permutes the r distinct copies over D simply transitively. This uses
complete cover geometry. It does NOT assert that the leading pole p
alone is coprime to r: several nonzero exponents can supply primitivity.
Nor is(3a) a sufficient realization criterion.

Let U' be the ENTIRE etale locus of the full Z'->Y. Vtilde is contained
in U', since it is etale over the already etale M. The generic deck action
extends regularly and freely to U' by the proved canonical-completion
theorem. These extensions do not erase any ramified prime of Z'.

## 4. The actual cyclic source image is exactly V

The group mu_r permutes the r copies D_j transitively and preserves their
common complement z!=0. Therefore
```
union_(sigma in mu_r) sigma(W_j)=Vtilde.                 (4)
```
The quotient of Vtilde by this group is the original V. This also
identifies V with the image of W_j in U'/mu_r; it is an open subset of
that quotient, with its original etale map to Y.

Surjectivity of W_j->V includes every closed point. A point of V outside
D has r distinct preimages in the common source overlap z!=0, all of
which belong to W_j. A point of D has one preimage in its selected copy
D_j, again belonging to W_j. Thus
```
# W_j-points over v = r for v in V minus D,
# W_j-points over v = 1 for v in D.                      (5)
```
The map is etale because it is the restriction of Vtilde->V.
For r>1 it cannot be finite: a finite etale morphism over connected V
would have constant fibre cardinality r. This is a surjective etale
plane pseudo-cover, with cyclic generic extension, rather than a
finite Galois cover of V by the selected plane.

The actual image algebra is also exact:
```
O(V)=O(W_j) intersect L inside L'.                       (6)
```
One inclusion is pullback. Conversely a rational function of L regular
on W_j descends to V along the surjective etale map, since its two
pullbacks to W_j times_V W_j agree. Thus no isolated-point completion
or unknown image algebra is introduced.

## 5. Whole fibre geometry, Picard group and fundamental group

On M the complete fibre c=0 is L_lambda+rD. Removing L_lambda leaves
```
c^-1(0)=rD on V.
```
For every c0!=0, M's whole fibre is the source line y=lambda+c0:
there is no other boundary on M. Thus every nonzero fibre on V is
one reduced A1. This proves the precise affine pseudo-plane structure,
not merely a generic A1-fibration.

The curve L_lambda is Cartier on smooth affine M, so its complement V
is affine. It is smooth and etale over Y. By the source localization
argument
```
Pic(M)=Cl(M)=Z[D].
```
The principal divisor of y-lambda gives [L_lambda]=-r[D].
The class-group localization sequence for removing L_lambda therefore
gives
```
Pic(V)=Cl(V)=Z[D]/rZ[D]=Z/rZ.                            (7)
```
This uses the actual full component open M with only D outside its
marked source.

The cover Vtilde is smooth affine and connected. It contains W_j=A2
openly, with complement the r-1 disjoint smooth A1 curves D_k, k!=j.
Moving loops off these real-codimension-two curves makes
pi1(W_j)->pi1(Vtilde) surjective, hence
```
pi1(Vtilde)=1.
```
The finite etale cyclic cover Vtilde->V is consequently the universal
cover, so
```
pi1(V)=mu_r,   chi(Vtilde)=r,   chi(V)=1.                 (8)
```
There is also a free Ga action on Vtilde. Its field is defined by
contraction_Delta(pullback Omega)=r*dz. In each whole W_j it is partial_t_j;
it is regular and locally nilpotent on O(Vtilde), since restriction
embeds that ring into the polynomial ring of W_j. Every fibre is a
union of its entire A1 orbits, and it has no fixed point.

## 6. A specific new omitted ramified prime and the r=1 switch

The original retained source line lifts to the WHOLE Ltilde in Z'.
Relative to the new plane source W_j it is omitted.

For r>1 its normal ramification index over the original source line,
and hence over its target image, is r. In source-cover coordinates
(x,z),
```
dP wedge dQ=r*z^(r-1)*dx wedge dz.
```
Thus Ltilde is an omitted RAMIFIED prime of Z'->Y, and it lies outside
U'. It has
```
v_(Ltilde)(z)=1.
```
If h_j has largest pole p>=r+1, formula(1) gives
```
t_j=z^(r-1)*(x-h_j(z)),
v_(Ltilde)(t_j)=-(p+1-r)<=-2.                            (9)
```
The new source volume dt_j wedge dz has order r-1 there.

For r>1 the new common boundary pencil is z, up to affine change.
Indeed some other D_k is an added unramified component relative to W_j,
and z is a regular source coordinate along it. The intrinsic coordinate
uniqueness from the common-pencil monicity argument makes z affine in
the common coordinate of K_j. Thus the ramified prime in(9) has finite
simple value for that NEW common coordinate. This is a marked ramified
configuration, not an exclusion of it.

For r=1, Z'=Z and Mtilde=M. The curve L_lambda is instead an omitted
UNRAMIFIED prime relative to W=M minus L_lambda, and the canonical
etale locus is the same U. The selected D has entered the new plane
source while the original source line has left it. Formula(1) is the
finite Laurent source-chart switch x=t+h(y-lambda), and the generic
degree remains d. The new map is not asserted equal to the original
source map under a polynomial automorphism.

## 7. The full-normalization limit cannot be dropped

The exact unit form on M restricts to an exact unit form on V.
Thus K_V is trivial even though Pic(V)=Z/rZ and r>1 can occur in the
conditional construction.

This does NOT permit an application of a torsion-canonical theorem for
the FULL FINITE normalization of an intermediate field to V.
Here the field fixed by mu_r is the original L, and its full finite
normalization over B is precisely the original Z. V=M minus L_lambda
is a PROPER selected etale open in that Z. It omits the entire source
line L_lambda and all original boundary components outside M.
Its canonical class and coordinate ring are not those of Z.

In particular neither V nor U'/mu_r is identified with the entire
original U or with a finite normalization order over B.
Their precise selected-open roles in(4)--(6) are retained.
No regularity or codifferent condition at any original ramified omitted
prime is supplied by constructing these plane charts.

The reduction gives the following concrete implication:
```
actual Keller map with some boundary multiplicity r>1
  ==> actual degree-rd Keller map with cyclic deck subgroup mu_r,
      an actual affine-pseudo-plane source image V as above,
      and the marked omitted ramified prime Ltilde.
```
It does not exclude that cyclic realization architecture. If all original
boundary multiplicities are1, only the same-field chart switches are
obtained. Neither case proves a source-finite or integral common
coordinate, injectivity on a source line, properness, or inversion.

## 8. Exact geometric control without a Keller target pair

For r>=2 put p=r+2 and consider the smooth affine surface
```
T_(r,p): u^p*a=b^r-1.
```
The r zero-fibre components u=0,b^r=1 are disjoint affine lines.
Each is contained in a whole plane chart
```
b=epsilon+u^p*s,
a=((epsilon+u^p*s)^r-1)/u^p,   epsilon^r=1.
```
These r charts cover T_(r,p). For epsilon=1 the inverse source chart
uses s=(b-1)/u^p=a/(1+b+...+b^(r-1)). Its omitted components are the
other r-1 lines.

The form
```
omega=db wedge du/u^p
```
is a global unit form. It restricts to ds wedge du on each chart and
has a global regular primitive db/((p-1)*u^(p-1)).
The action
```
(u,a,b) |-> (zeta*u,zeta^(-p)*a,zeta*b)
```
has order r, is free and preserves omega, since p=r+2.
Its quotient is smooth affine with the same one-multiple-fibre
geometry, and T_(r,p) is simply connected by its marked plane chart.
This control shows why exact unit volume, free cyclic action and the
affine pseudo-plane geometry alone are not declared an inverse proof.

No etale target pair to A2 is supplied for this control. The actual
target pair in Sections1--7 comes ONLY from the assumed original
Keller map and the entire proved root chart.

The checker verifies exact field/coordinate/Jacobian identities, cyclic
characters, the marked ramified valuation, fibre and Picard ledgers,
the finite-cover Euler formulas and the whole affine controls.
Written arguments establish the unrestricted geometry. No formal
verification or external human review is asserted.
