# Noncontact (6,8): the cubic edge reduces to three necessary fibre profiles

14 November 2026. PF2RULED13 retains a necessary local calculation
for noncontact(6,8). The subsequent
[global cusp-jet theorem PF2RULED14](GLOBAL_CUSP_LEADING_OBSTRUCTION.md)
excludes all three profiles as global finite quartic realizations.
The local curves and this calculation remain exact controls; they are
no longer active realization targets. This source by itself excluded
slope4 and left the three necessary profiles below.

| First slope | Additional low coefficient | High depressed cubic ord(a),ord(b),ord(Disc) | Infinity Q poles |
| --- | --- | --- | --- |
| 6 | no lower cluster | 20,30,61, with leading discriminant cancellation | 3,1 |
| 5 | fifth z-coefficient a unit | >=14,>=21,42 | 1,1,1,1 |
| 5 | fifth z-coefficient vanishes | >=14,20,40 | 1,1,1,1 |

The orders refer to a monic Weierstrass cubic written Z^3+a(w)Z+b(w)
after translating away its quadratic term, on a generic P fibre. They
are exact necessary jet conditions, not numerical evidence of realizability.
The three local curve controls retained below supply no global P,Q,
no proof that the fibre differential has zero residues, and no finite
map on the Wright ring.

## Input and the first-edge alternatives

Use the full compatible geometry and the reductions of
[PF2RULED11](BOUNDARY_LEADING_PERIOD.md) and
[PF2RULED12](POLE_FIVE_EIGHT_EXCLUSION.md). In the noncontact case,
leading symbols of P,Q with poles6,8 are sigma^3,sigma^4. The common
sextic sigma has divisor6p, p!=q, hence sigma=tau^2 with the cubic
tau having divisor3p. Normalize symbols to tau^6,tau^8. At p choose
regular coordinates w,z, with w a base parameter of the ruled surface
and B=(z=0). The volume is a unit times z^-3 dz wedge dw.

For generic c0 write

 F=z^6(P-c0)=w^18*unit+sum_{i=1}^5 A_i(w)z^i+z^6*C(w,z,c0).

Let lambda be its largest first Newton slope from (0,18), as in the
coprime proof. If lambda=b/a>41/12, put c=b-3a>0. The weighted leading
ODE and the actual Q pole bound force

 G^3=constant*H^4,
 H=(A+sum_{i=1}^5 B_i Z^i)/Z^6, A!=0.

Indeed (6+8-2)c>5a is exactly the strict volume inequality. Thus H
must be a rational cube. Its nonconstant numerator, of degree at most
five, must be the cube of a linear polynomial. In particular it has
a nonzero linear term and degree three. Since the z coefficient is
present, lambda is an integer; the nonnegative w-order of the z^3
coefficient bounds it by six. The only such slopes above41/12 are

 lambda=4,5,6, with first face (A+B Z)^3, A,B!=0.          (1)

If3<lambda<=41/12, the complete integer order list is at least
15,12,8,5,1. The possible slopes are17/5,10/3,13/4,16/5.
The normalized fibre differential is

 dQ/dt=unit*z^3*w'(t)/F_z.                               (2)

The general first-edge lower bound from PF2RULED12 gives Q pole bounds
17,11,16,21 respectively at one missing point, all exceeding degree4.
This removes the entire intervening band, retaining repeated face roots.
If lambda<=3, the coefficient cone forces P to be a polynomial in a
V1 lift T of tau, excluded by the polar centralizer. Therefore(1)
is the complete first-edge list for the surviving noncontact pair.

## A global cube completion and strict remainder orders

All possible first slopes are at most six, so the z coefficient of F
has w-order at least12. Lift tau to T in V1 and write P=T^6+R5.
The degree-fifteen symbol of R5 must be divisible by tau^4. Lift its
quotient divided by3 to A1 in V1. Then

 S=T^2+A1 in V2, P=S^3+R4, R4 in V4.                    (3)

This is an exact global statement. Put j=6-lambda in{0,1,2}. The
cubic symbol of A1 has order j at p; in local Laurent form

 S=z^-2*(sigma(w)+beta(w)z+gamma(w,z)z^2),
 ord(sigma)=6, ord(beta)=j.

The perfect-cube first face of P fixes its linear cube root and makes
its first face equal to the first face of S^3. Hence R4 lies strictly
above that edge. If

 z^6*R4=r2(w)z^2+r3(w)z^3+r4(w)z^4+r5(w)z^5+O(z^6),

its necessary orders are

| j | lambda | ord(r2) | ord(r3) | ord(r4) | ord(r5) |
| --- | --- | --- | --- | --- | --- |
| 0 | 6 | >=7 | >=1 | >=0 | >=0 |
| 1 | 5 | >=9 | >=4 | >=0 | >=0 |
| 2 | 4 | >=11 | >=7 | >=3 | >=0 |

The Q first face is also fixed: after the leading target scaling it
is the first face of S^4. Thus Q-S^4 is in V7 and has strictly higher
weighted value than the first face. This does not identify S with a
target coordinate or remove arbitrary odd formal-root terms.

In particular r2 is divisible by sigma. For the surviving j=0,1
after the next exclusion, the restriction sequence gives another exact
necessary global form with auxiliary functions X,Y:

 P=S^3+3S*X+3Y, S,X in V2, Y in V3,                     (4)

where the symbols of X,Y have orders at least2j+1,3j+1 respectively.
All finite-dimensional symbol quotients used here are global sections;
no inversion at p or a choice of a special lift is imposed.

## Slope four is impossible

For j=2, the endpoint of the first edge is (z degree3,w order6).
The three roots of the cubic face form a cluster of total w-projection
degree three. Any branch over it has orders(e,4e), and(2) gives a Q
pole at least e, because the F_z bound is(18-4)e. Summing over the
cluster gives at least three Q poles counted with order.

The next edge starts from (3,6). The strict remainder orders imply
ord(A4)>=3 and ord(A5)>=0. The generic z^6 coefficient is nonzero.
Consequently the largest next slope b/a is between two and three.
There is at least one normalized branch over a nonzero root on that
edge, disjoint from the first cluster. Its F_z order is at least
(6a+2b)e, so its Q pole is at least(5a-b)e>=2.
The total is at least five, impossible for a finite quartic map.
This excludes the entire slope-four branch, with repeated roots allowed.

We henceforth have j=0 or1, so the symbol of A1 is a unit or has a
simple zero at p. Contact does not enter these computations: p!=q.

## Pole accounting by a Weierstrass discriminant

For a local monic Weierstrass factor W of projection degree D, sum
(2) over its normalized branches. If their w,z,F_z orders are a_i,b_i,s_i,
the signed Q pole sum is

 sum(s_i-a_i-3b_i).

For the whole local curve, sum(a_i)=D, sum(b_i) is the w-order of its
constant term, and sum(s_i) is the discriminant order. Units do not
change these values. Finiteness makes every missing point an actual
Q pole, so signed and actual pole sums agree. These formulas remain
valid over all components; no irreducibility is assumed.

### Slope six

For j=0 the entire local Weierstrass degree is three: the coefficient
of z^3 at w=0 is beta(0)^3!=0. All infinity points lie here and the
constant term has w-order18. Thus total degree four requires

 ord(Disc W)=4+3+3*18=61.                               (5)

Translate W to Z^3+a(w)Z+b(w). Its discriminant is -4a^3-27b^2.
If r=ord(a),s=ord(b), unequal3r,2s could not give61, which is divisible
by neither two nor three. Thus3r=2s<61, with r=2t,s=3t, t<=10,
and the leading discriminant cancels. The first cubic cluster makes
t>6. Its leading coalescent cubic has a simple root and a double root;
the simple branch is smooth and has F_z order2t. Its z-order is six,
so its Q pole is2t-19. Finiteness forces positivity, hence t>=10.
Therefore t=10, r=20,s=30 and(5) has leading cancellation at order60.

The simple branch contributes a Q pole of order one. The other two
projection sheets form one cusp branch, contributing Q pole three.
For example the double-factor discriminant has order21: the resultant
with the separated simple factor has order20 and is squared in the
full discriminant61. Its odd discriminant gives ramification two.
Thus the necessary infinity pole partition is(3,1), with w-projection
indices(2,1). This proves the first row of the stated table.

### Slope five: two low-coefficient obstructions

For j=1 the first cubic cluster has z-order five. The actual z^4
coefficient cannot have order zero: the next edge would have slope3
and a simple residual linear root, giving a regular Q at a missing
point. Order one would give slope2 and a nonzero simple pole in dQ,
impossible for a meromorphic derivative. Thus ord(r4)>=2; the S^3
contribution already has order at least two and cannot hide these cases.

Now there are exactly two low-cluster alternatives.

* If r5 is a unit, the next edge has slope3/2. It gives one cusp
  branch with (ord(w),ord(z),ord(F_z))=(2,3,12), hence one Q pole.
  There are no further low branches; the local Weierstrass degree is5.
* If r5 vanishes, the next edge has slope1 and its residual cubic is
  simple for generic c0. Its leading coefficient varies with -c0,
  so its discriminant is not identically zero. It gives three smooth
  branches with orders(1,1,5), each a Q pole of order one. The local
  Weierstrass degree is6.

These statements retain arbitrary regular S lifts and coefficients.
The residual cubic's constant is beta_lead^3!=0, and its discriminant
has a nonzero quadratic leading coefficient in the generic fibre
parameter, so no fixed parameter choice can remove the three branches
for every generic c0.

Factor the local Weierstrass polynomial into its monic high-cluster
cubic H and its low factor. On a high branch, the low factor has
w-order three: the high roots have order five and all low roots have
smaller order, with total low z-order three. Its contribution to the
sum of F_z orders over the high cluster is therefore nine. After
translating H to Z^3+aZ+b, the high Q pole total is

 ord(Disc H)+9-3-3*15=ord(Disc H)-39.                    (6)

In the first low alternative, it must be three, so ord(Disc H)=42.
In the second it must be one, so ord(Disc H)=40.

To sharpen the cubic coefficients, a cancellation below either of
these discriminant orders would require r=2t,s=3t and t<=6. The
simple branch of its coalescent leading cubic would have original
F_z order3+2t and Q pole2t-13, negative. This makes Q regular and violates finiteness. Thus such cancellation is
excluded. A discriminant of order42 then requires r>=14,s>=21,
whose leading cubic after Z=w^7*T has distinct roots, including any
zero root. Its three high branches are smooth, each a Q pole of order one. A discriminant of order40
instead requires r>=14,s=20; its leading binomial has Puiseux type
(3,20). It is one high cusp branch with a Q pole of order one.
Together with the respective low branches, both rows have four
simple Q poles. This proves the last two rows of the table.

## Scope and the next exact gates

Every remaining noncontact(6,8) map must satisfy one of these three
profiles, the global forms(3)-(4), the fixed Q first face, and the
original constant Jacobian and finiteness conditions. The next tests
are global realizability of the stated Weierstrass jets, vanishing of
all required fibre-differential residues, and existence of the actual
quartic Q function. A profile satisfying pole counts is not such a
function. The marked cusp tests of PF2OMIT9/10 also remain necessary.

The [checker](verify_pole68_cube.py) verifies the complete first-edge
and strict-order lists, the slope-four lower-edge budget, the cubic
discriminant arithmetic and explicit local curve controls for all
three profiles. The [receipt](pole68_cube_receipt.json) separates those
controls from the written Weierstrass factorization, global section
arguments, normalization coverage and finiteness proof. No formal
verification, independent full reproof, external review or novelty
priority is claimed.

The first open height-eight pairs remain(6,8),(7,8). The latter has
only its contact72,108 alternative by PF2RULED12. The contact(6,8)
alternatives, both remaining all-height pole-four root types, arbitrary
higher poles and actual ruled compatibility are all still open.

Current handoff: [PF2RULED14](GLOBAL_CUSP_LEADING_OBSTRUCTION.md) globally
excludes the entire noncontact family under the actual marked-model
hypotheses. The local pole/discriminant controls above remain valid.
Only contact72,108 cases remain at the height-eight frontier.
