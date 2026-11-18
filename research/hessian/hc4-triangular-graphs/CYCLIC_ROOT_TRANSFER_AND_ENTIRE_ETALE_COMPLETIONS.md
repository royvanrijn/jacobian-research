# Cyclic root transfer and the entire canonical etale completion

18 November 2026. This is an exact comparison of the FULL canonical
completions before and after a whole-root-chart transfer of an actual
complex plane Keller map. It does not exclude the transfer architecture,
force a boundary multiplicity to be one, or prove inversion.

The selected pseudo-plane of the
[whole-root transfer](WHOLE_ROOT_CHARTS_CREATE_CYCLIC_KELLER_PSEUDOPLANES.md)
need not exhaust the quotient of the new completion. Here that ENTIRE
quotient and every added component of the new completion are identified.

## 1. Setup and exact statement

Let
```
K=(P,Q): X=A2_(x,y) -> Y=A2,  Jac(P,Q)=1,
B=C[P,Q], F=Frac(B), L=C(x,y), d=[L:F].
```
Let Z be the FULL finite normalization of B in L, and let U be its
ENTIRE etale locus over Y. The
[common-pencil theorem](EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md)
and [whole-chart theorem](COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md)
give disjoint whole smooth affine lines
```
U minus X = D_1 disjoint-union ... disjoint-union D_s,
y|D_i=lambda_i, v_Di(y-lambda_i)=r_i>=1,
div_U(y-lambda)=L_lambda+sum_(lambda_i=lambda) r_i D_i,   (1)
L_lambda={y=lambda} subset X.
```
The displayed source line and all D_i are closed in U. In particular the
supports in(1) are disjoint. The form Omega=dP wedge dQ is a regular exact
unit form on U.

Choose D_1, put lambda=lambda_1, and assume
```
r=r_1>1.
```
Put L'=L(z), z^r=y-lambda. As in the whole-root transfer, let Z' be the
FULL finite normalization of B in L', and choose the entire plane source
W_1=A2_(t,z) supplied by the component chart of D_1. The target map
K':W_1->Y is polynomial Keller, with Jacobian r and generic degree rd.
Let U' be the ENTIRE etale locus of Z'->Y. The finite map q:Z'->Z is
the normalized cyclic root cover, without assuming that y is regular
or integral on all Z.

Partition the original added components into
```
I_good={i : lambda_i!=lambda OR r divides r_i},
I_bad ={i : lambda_i=lambda AND r does not divide r_i},
g=|I_good|>=1.
```
Then the exact quotient and inverse-image identities are
```
W=U minus (L_lambda union union_(i in I_bad) D_i),
U'=q^(-1)(W),             U'/mu_r=W.                     (2)
```
The restricted map U'->W is finite etale of degree r. In particular no
point above Z minus U belongs to U'. Every original omitted ramified
prime remains outside the new etale locus.

For each good original D_i, its inverse image in U' consists of r
disjoint WHOLE affine lines, each mapping isomorphically onto D_i.
Exactly one of those r*g lines belongs to the marked new plane W_1.
Consequently
```
b'=number of components of U' minus W_1 = r*g-1,
chi(U')=r*g,                 chi(W)=g,
pi1_top(W)=mu_r,
Pic(W)=Z^(g-1) direct-sum Z/rZ.                          (3)
```
Here b' counts added unramified divisors; it is not a map degree.

The new common-coordinate multiplicities are also exact. If
lambda_i!=lambda, the r copied lines have distinct constant z-values
alpha with alpha^r=lambda_i-lambda, and each has
```
v(z-alpha)=r_i.
```
If lambda_i=lambda and i is good, all r copies have z-value zero and
```
v(z)=r_i/r.                                             (4)
```
The r copies belonging to D_1 therefore have multiplicity one.

For r=1 the extension and full normalization are unchanged:
Z'=Z, U'=U. The source switch replaces D_1 by the original source line
L_lambda in the added boundary, so b'=s. Formula(2), which removes the
ramified source line, is asserted only for r>1.

## 2. An etale quotient cannot create new original etale points

The root extension has degree r, by the valuation one of y-lambda at
the original source line. Its cyclic group mu_r acts on the full Z'
over Z, and
```
Z'/mu_r=Z.
```
These statements concern full finite normalizations, not a selected
source-image algebra. They were proved in the whole-root transfer.

The group preserves U', since it fixes P,Q and hence the etale property.
Every geometric fibre of the finite quotient q is one group orbit.
Thus q^(-1)(q(U'))=U'. The set W_0=q(U') is open in Z: its complement is
the image of the closed invariant subset Z' minus U' under the finite
map q. Hence U'/mu_r is exactly this open W_0, with its original
scheme structure inside Z.

The generic deck action on the canonical etale completion is free by
[the canonical-topology theorem, Section6](CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md).
For clarity its fixed-point argument applies directly here. The fixed
locus of a nonidentity deck transformation in connected U' is the
inverse image of the diagonal of the separated etale map U'->Y.
That diagonal is open and closed. The fixed locus is therefore open
and closed, and is empty because the generic transformation is
nonidentity.

A free action of the finite constant group mu_r on affine U' makes
U'->W_0 a finite etale torsor. Since U'->Y is etale, its quotient
W_0->Y is etale as well. One may use the precise surjective-etale
descent statement [Stacks02KM](https://stacks.math.columbia.edu/tag/02KM),
or its locally-finite-presentation version
[Stacks02K6](https://stacks.math.columbia.edu/tag/02K6).
The finite-type open W_0 of Z meets that hypothesis. Therefore
```
W_0 subset U,       U' subset q^(-1)(U).                 (5)
```
This is an all-points assertion. It does not merely rule out newly
unramified divisors while overlooking isolated points.

The free quotient statement is also covered by
[Stacks07S7](https://stacks.math.columbia.edu/tag/07S7):
the finite locally free group and affine source hypotheses hold,
and its torsor is finite etale because mu_r is finite etale over C.
The displayed proofs of these primary statements were inspected;
their classical descent and quotient foundations are adopted.

## 3. Normalized local roots classify every point over U

Because y is regular on U, the inverse image q^(-1)(U) is the
normalization of the root equation z^r=y-lambda there. This follows
by localization of the full finite normalization; it makes no assertion
that z is integral on all Z.

Away from the fibre y=lambda, the right side is a unit. The raw equation
already defines a finite etale algebra, since r*z^(r-1) is a unit.
Consequently every point there belongs to U'.

At a point of one of the disjoint smooth curves in(1), choose a local
equation f for that curve. Smoothness of U and the complete fibre
identity give
```
y-lambda=u*f^m
```
with u a local unit. Pass to a strict henselian local neighbourhood,
where u has an rth root. With
```
h=gcd(r,m),     e=r/h,     n=m/h,
```
the normalized equation has h branches. Up to a constant root of unity
and a unit, each has the parameter description
```
f=tau^e,        z=tau^n.                                (6)
```
The transverse map has ramification index e. These descriptions also
apply at every finite point of the whole curve: the remaining curve
parameter is unchanged, and the normalizations in(6) are smooth.
Thus the cover is etale there exactly when r divides m.

For the original source line, m=1 and e=r>1, so none of its inverse
image belongs to U'. For D_i at lambda_i=lambda, m=r_i; its inverse
image is etale precisely for the listed good indices. This proves(2),
including the absence of finite-point exceptions.

More explicitly, the reduced inverse image of a whole such D_i is
a degree-h etale cover of D_i=A1. It therefore splits into h entire
affine lines. This uses the elementary complex A1 covering fact already
used in the whole-chart theorem. Each has transverse index e. At a bad
component these are h=gcd(r,r_i) ramified lines, each of index r/h;
they are retained in full Z', while excluded from U'. At a good
component h=r and e=1, giving the r asserted unramified copies.
The same r-copy assertion off the selected fibre follows by restriction
of its finite etale root cover to D_i.

For a retained copy over a distinct lambda_i, factor
z^r-(lambda_i-lambda) at its nonzero root alpha. Its derivative is
nonzero, so z-alpha and y-lambda_i differ by a unit, proving(4).
At a retained same-value component, (6) has e=1 and n=r_i/r,
which proves the other case of(4).

Every original omitted ramified prime E has target ramification
index e_E>1. For every prime E' above it, valuation ramification
indices multiply:
```
e(E'/Y)=e(E'/E)*e_E>1.
```
Such primes cannot enter U'. Statement(5) also covers every closed
point outside U, without assumptions on smoothness of those original
ramified curves or of the full Z.

## 4. The new source and complete boundary ledger

Over X minus L_lambda the normalized root cover is exactly
```
A1_x times G_m,z.
```
This is also W_1 minus its selected zero line. The entire W_1 is that
open together with the ONE copy of D_1 selected in the whole-root chart.
All other points of U' lie on the r copies of the good D_i. Those
curves are disjoint and wholly covered. Hence the complete new boundary
consists of r*g-1 lines, as in(3).

In particular the entire new common-coordinate fibration, rather than
only its restriction to the selected pseudo-plane cover, has no other
added unramified components. Since r>1, the other r-1 copies of D_1
provide an added component where z is regular. Common-pencil uniqueness
therefore identifies the new common coordinate with z up to affine
change. The multiplicities in(4) are multiplicities for this actual
new common pencil.

If originally s=1, all r-1 new added components have multiplicity one.
The new translation partial_t is then free on the ENTIRE U': it is
translation on the new source, and on every added multiplicity-one
line in its whole chart. This special case still retains the newly
omitted ramified original source line outside U'; it is not an inverse
theorem.

More generally the original source line above lambda is the whole
newly omitted ramified curve described in the root-transfer theorem.
Its index is r and its common-coordinate value is finite:
```
v_(Ltilde)(z)=1.
```
It is not converted into a common-coordinate pole by the completion
comparison.

## 5. Quotient invariants and the saturation of the two Ga actions

The complement defining W is a union of disjoint effective Cartier
divisors in smooth affine U; hence W is affine and smooth. Euler
additivity and chi(U)=1+s give
```
chi(W)=(1+s)-1-|I_bad|=g.
```
The transferred map is an actual Keller map on W_1, so its canonical
etale completion U' is simply connected by the retained topology
theorem. The finite etale cyclic cover U'->W is therefore universal,
giving pi1_top(W)=mu_r and chi(U')=r*g.

For Picard groups use Pic(U)=direct-sum_i Z[D_i].
Removing the bad components kills their basis classes. Removing the
source line kills its class
```
[L_lambda]=-sum_(lambda_i=lambda) r_i[D_i].
```
In the remaining good basis this is one nonzero relation whose
coefficient gcd is exactly r: every same-value good coefficient is
divisible by r, and the chosen coefficient r_1 equals r.
The integer Smith normal form therefore gives the Picard group in(3).
There is no torsion-canonical assertion for the full normalization:
W remains the explicitly identified open subset of the original Z.
When g=1, W is exactly the previously selected pseudo-plane V.
For g>1, W contains V and its complement is the g-1 other good
whole affine lines.

The two Ga fields explain a genuine limit of a tempting descent.
On the new source the divided field is delta_new=partial_t, whereas
the original source translation lifts as
```
delta_old=partial_x=z^(r-1)*partial_t.                    (7)
```
The old field L=C(x,y) is stable under delta_old. It is NOT stable
under delta_new for r>1: from the chart,
```
delta_new(x)=z^(1-r),
```
which is not fixed by mu_r and hence does not belong to L. This
failure holds regardless of the finite Laurent principal part.

The [direct replica index theorem](DIRECT_REPLICA_INDEX_FROM_A_DESCENDED_DIFFERENTIAL.md)
applied to the actual new Keller map and its stable intermediate L
gives exactly
```
[L':L]=r divides deg_z(z^(r-1))+1=r.                      (8)
```
Thus the action comparison saturates its index bound. Dividing the
replica inside the old field, or invoking a minimal-degree map after
a construction that changes d to rd, would require a further argument.
No such descent is asserted.

## 6. Scope and exact controls

The theorem compares full normalizations and ENTIRE etale loci.
It keeps all original ramified primes, all bad copied components,
and all finite points before naming the etale opens. It proves
neither that r=1 nor that all original boundary components disappear.

The fixed field is the original L and its full finite B-normalization
is the original Z. Even after the exact identification U'/mu_r=W,
the quotient is not identified with Z or U. A theorem for the
canonical class of the full finite normalization cannot be applied
to this selected open. No source-finiteness or target integrality of
the common coordinate, properness, injective source line, or inverse
is proved.

The checker records exact local-root normalization exponents and
all-sheet ramification counts, complete divisibility-ledger controls,
Picard relation Smith forms, Euler and covering counts, the r=1
exception, and the action/index saturation. Its finite models do not
assert existence of an actual noninvertible Keller map. The unrestricted
comparison is proved in Sections2--5; no coefficient or companion
degree cutoff is used. No formal verification or external human
review is claimed.
