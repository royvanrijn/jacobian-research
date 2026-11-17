# Ramified fold valuations and poles of the common source coordinate

A ramified omitted prime of a plane Keller normalization is constrained by
the actual open affine-plane source, not just by its canonical label.
This note proves a general determinant obstruction and applies it to
ramification index two. With a nonempty unramified boundary, the minimal
second extraction pattern in Borisov's proof cannot occur at index two.
For the first pattern, a polar common coordinate has transverse degree one,
and only two numerical multiplicity types remain.

These are necessary conditions. No existence of a Keller map with either
remaining type is asserted; no ramified whole-component chart or inverse
criterion is supplied.

## 1. Setting and conventions

Let K=(P,Q):X=A2_C -> Y=A2_C be polynomial with nonzero constant
Jacobian J. Put R=C[x,y], B=C[P,Q], L=Frac R, and let S be the finite
normalization of B in L. The canonical lift X -> Z=Spec S is an open
immersion. Let U be the entire etale locus of Z -> Y. An omitted prime
means a prime divisor E of Z outside X. Its primitive valuation is v_E,
its ramification index over Y is e, and

```
ord_E(dx wedge dy)=e-1,       kappa(E)=e.                   (1)
```

The [canonical boundary proof](CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md)
establishes these facts in the full source field. Normality makes a
neighbourhood of the generic point of E smooth, even when Z has singular
closed points.

Extract all valuations under discussion on one smooth surface obtained
by point blowups outside the marked X from a source-coordinate P2. For
the complete boundary basis, define E^vee by

```
E^vee . F = delta_(E,F),       n_E=-(E^vee)^2.
```

These unnormalized dual pairings are unchanged by further blowups.
Borisov's determinant label d(E) equals (E^vee)^2, so n_E=-d(E).

Actual P and Q are regular at every omitted prime and are algebraically
independent. Resolve [P:Q:1]; the positive-square pullback of target
infinity is supported outside all omitted primes in question. Their
duals are orthogonal to it. The signature (1,N) of a blown-up P2
therefore gives a negative definite simultaneous dual matrix. In
particular n_E>0, and for distinct omitted D,E,

```
n_D*n_E > (D^vee . E^vee)^2.                              (2)
```

This is the same direct positive-pullback argument as
[the common-pencil proof, Section2](EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md).
It uses the actual functions P,Q; no converse valuation-intersection
theorem is being assumed.

## 2. A small dual norm forces every source coordinate to have a pole

For ANY omitted prime E and ANY polynomial source coordinate b,

```
v_E(b)>=0  implies  n_E>=e+1.                             (3)
```

Consequently n_E<=e implies that EVERY polynomial source coordinate is
polar at E.

Choose a polynomial companion a. If both a and b were nonnegative,
the valuation would have a centre on X. Its centre on the separated
Z could then not be the different generic point of E. Hence, whenever
b is nonnegative, v_E(a)=-p<0 for an integer p.

If v_E(b)=0 and its residue were transcendental over C, every nonzero
coefficient in C[b] would have valuation zero. For
f=sum_i f_i(b)*a^i the highest nonzero a-power would then be the unique
term of least valuation. Thus

```
R intersect O_E=C[b],
```

contradicting B subset R intersect O_E and trdeg_C B=2.
Therefore b has constant residue lambda and

```
v_E(a)=-p<0,       v_E(b-lambda)=r>0.                     (4)
```


Translate b by lambda. In the source P2 frame (a,b), the first centre
is the a-infinity point. After that blowup the centre is the finite-b
free point with local parameters

```
u=1/a, b,       v_E(u)=p, v_E(b)=r.
```

Let m_1,...,m_l be the positive point multiplicities in the extraction
over this point. The local canonical calculation and the cancellation
of the first infinity contribution to the dual square give

```
sum_j m_j=2p+e-1,       n_E=sum_j m_j^2.                  (5)
```

Indeed du wedge db=-u^2 da wedge db, so its order is 2p+e-1,
equal to the sum of the coefficients of the successive relative
canonical exceptional divisors. In the total-transform basis,
E^vee=p(H-F_0)-sum m_j F_j, proving the second identity.
Since m_j are positive integers,

```
n_E>=sum m_j=2p+e-1>=e+1,
```

as claimed. No smoothness of the whole curve E or root chart for it
is used.

## 3. The unramified boundary is an anchor for every polar omitted prime

Assume U minus X is nonempty. The proved
[individual-coordinate theorem](INDIVIDUAL_SOURCE_COORDINATES_AND_UNCONDITIONAL_BOUNDARY_CHARTS.md),
[whole-chart theorem](COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md),
and [alignment theorem](EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md)
give one polynomial source coordinate y regular on U. For each
unramified omitted D there is a polynomial companion x_D with

```
v_D(x_D)=-p_D,  v_D(y-lambda_D)=r_D,
p_D>r_D>0,     0<n_D<=p_D^2.                             (6)
```

If ANY other omitted E has v_E(y)=-a<0, then

```
D^vee . E^vee=p_D*a.                                     (7)
```

This does not require an adapted source coordinate at E. In the P2
frame (x_D,y), if E starts at a different infinity point, its infinity
weight is a and the clusters share no centre. If it starts at the
x_D-infinity point, write its infinity weight b>a; its first
multiplicity is b-a and its next centre is the intersection with
strict infinity, whereas D goes to finite y=lambda_D. Thus only the
first point is shared. The dual pairing is respectively p_D*a or
p_D*b-p_D*(b-a), proving (7).

Combining (2), (6), (7) gives

```
n_E > a^2.                                               (8)
```

In particular, if a ramified prime has n_E<=e and e<=a^2, it is
impossible: (3) first forces y to be polar, then (8) contradicts its
norm.

The hypothesis U minus X nonempty is essential to this anchor argument.
Distinct unramified primes may have the SAME lambda_D; their dual
cross pairings need not be zero. No bound obtained by adding one
independent a^2 contribution per component is asserted.

## 4. Minimal second extraction pattern at index two is excluded

Choose a source-coordinate P2 compactification minimizing the number
of ancestors of E. Use precisely Case1/Case2 of Borisov's proof of
Theorem4.3 on printed pp12--15, not an arbitrary blowup presentation.

For minimal Case2, the first positive canonical label is 1, with
determinant c>=0. After m>=0 blowups of the adjacent 0--1 edge, creating
label 2 freely gives

```
d_2=(m+1)^2*c+m^2+m-1.                                   (9)
```

Creating label 2 at a 1--1 edge instead requires m>=1 and gives

```
d_2=(2m+1)^2*c+4m^2-1.                                  (10)
```

Subsequent ancestral canonical labels strictly increase, so an eventual
label 2 must be that first label-2 divisor. Formula (10) is positive.
The negative instances of (9) have m=c=0 and d_2=-1. The actual
negative sign proved in Section1 therefore forces n_E=1.

By (3), every polynomial source coordinate, including the common y,
must be polar at this E. By (8), n_E>v_E(y)^2>=1, a contradiction.
Thus, if U minus X is nonempty, an actual ramified index-two prime
cannot belong to minimal Case2.

The imported part is the same inspected minimality analysis used in
the individual-coordinate proof, now followed through the two explicit
label-2 formulas on printed p15. We adopt its classical contraction
foundations. The primary is Alexander Borisov, *On two invariants of
divisorial valuations at infinity*, author version dated 13 February
2013, published J. Algebraic Combinatorics39(2014),691--710,
DOI10.1007/s10801-013-0462-9:
[retained primary PDF](sources/borisov-divisorial-valuations-2013.pdf),
[author source](https://people.math.binghamton.edu/borisov/documents/papers/divisorialvaluations.pdf).
The formulas and the restriction m>=1 in (10) were read in that primary.

## 5. The remaining first-pattern fold: transverse degree one

Suppose e=2 and E belongs to minimal Case1. Its first two centres are
free, so it has source coordinates (a,b) satisfying (4).
The multiplicities in (5) never increase and each is at most p.
Therefore

```
sum m_j=2p+1,       n_E<=p*(2p+1)<4p^2.                  (11)
```

The last strict inequality holds for every p>=1.

If the common y is polar at E, it is not affine in b. The
coordinate-monicity lemma in the common-pencil proof states that
m=deg_a(y)>=1 and the leading a-coefficient is a nonzero constant.
Consequently

```
v_E(y)=-m*p.
```

Equations (8),(11) force m=1. Thus

```
y=c*a+H(b),       c in C*, H in C[b],                     (12)
```

and (y,b) is itself a polynomial source coordinate pair. The remaining
norm restriction is n_E>p^2.

This statement constrains only primes at which common y is polar.
A first-pattern prime with nonpolar common y is not excluded or
classified here.

## 6. Exact multiplicity types in the polar fold case

Under Section5's polar hypothesis, the complete, unbounded extraction
sequence has just the following possible numerical types:

```
(p,r,n_E)=(2,1,5),
or
r=p+1,       n_E=p^2+p+1,       p>=1.                     (13)
```

Here b has been shifted by its constant residue. These are necessary
types, not existence assertions.

We give the extraction arithmetic and the two actual-source exclusions.
Put g=gcd(p,r). The subtraction Euclidean prefix, including its final
equal-value blowup, has sums

```
sum_prefix m=p+r-g,       sum_prefix m^2=p*r.
```

If the equal-value ratio has transcendental residue the prefix is
terminal; otherwise one continues above its finite nonzero residue.
Every later multiplicity is at most g. By (11), its tail has mass

```
T=p-r+g+1>=0,       n_E<=p*r+g*T.                         (14)
```

For an extraction of at least two point blowups, the last TWO
multiplicities are 1. Indeed the final centre parameters both have
value 1 on the newly created primitive exceptional divisor. Looking
one blowup back, unequal values are related by subtraction, or one
parameter retains the common value at an equal-value blowup. In
either case that preceding centre also has minimum 1.
The sequences here have length at least two since their total
mass is 2p+1>=3.

If p>r, write p=gA,r=gB with A>B>=1. Formula (14) gives

```
n_E-p^2 <= g^2*(1-(A-1)*(A-B))+g.
```

Since n_E>p^2, the positive integers force A=2,B=1.
The prefix is then [g,g] and the tail mass is 2g+1.
Remove its terminal 1; the remaining tail squares are at most
g*(2g). Hence n_E<=4g^2+1. Strict n_E>4g^2 forces equality,
so every remaining tail entry equals g. The penultimate 1 then
forces g=1. This gives (p,r,n_E)=(2,1,5), with five multiplicities 1.

If p<r, nonnegative T means g*(B-A-1)<=1 for
p=gA,r=gB. Either B-A=1 and T=1, or g=1,B-A=2 and T=0.
In the first case the tail is [1]; the preceding last Euclidean
entry is g, so the last-two-1 property forces g=1. This gives
r=p+1 and n_E=p*r+1.
In the second case r=p+2, gcd(p,r)=1, p is odd and the
Euclidean block is terminal. Such a valuation is the MONOMIAL
valuation in the local parameters (u=1/a,b): the last equal-value
ratio is transcendental, and no further finite-residue correction
occurs. Thus for source polynomials it has weights

```
v_E(a)=-p,       v_E(b)=p+2,
```

with the valuation equal to the minimum of the monomial weights.
Every nonconstant a^i term without a factor b has negative value.
Hence every f in R intersect O_E restricts to a constant on b=0.
In particular actual P and Q would both be constant on that entire
source line, forcing their Jacobian to vanish there. This terminal
case is excluded by the actual Keller map.

Finally p=r is excluded by the actual regular primitive.
Scale (a,b) by a constant so da wedge db=J^-1 dP wedge dQ,
and put theta=J^-1 P dQ. On the source polynomial plane,

```
theta=a db+dF,       F in C[a,b],                         (15)
```

by the polynomial Poincare lemma. At the generic point of E,
theta is a regular differential: P,Q are regular on Z there.
After a finite etale extension of the generic DVR take a local
uniformizer z with b=z^r; a separating tangential parameter can be
fixed in the completed smooth neighbourhood. Along the resulting
formal transverse arc, theta has zero dz/z coefficient and so does
the derivative of the Laurent series F(a(z),z^r).
If p=r, the leading a-term is a_-r*z^-r with a_-r nonzero;
a db has dz/z coefficient r*a_-r, contradiction.
The root cover is etale at E because it only extracts an r-th root
of the unit in b=u^r*unit; no ramified whole chart is assumed.
This residue argument in fact excludes p=r whenever the pulled-back
volume order is nonnegative and theta is regular.

## 7. Controls and remaining scope

The two types in (13) pass the elementary local volume/primitive tests.
For r=p+1, the formal model

```
b=z^(p+1),       a=z^-p+z^(1-p)*t
```

has da wedge db=(p+1)*z*dt wedge dz and a db regular.
For (p,r)=(2,1), the model b=z,a=z^-2+z*t has
da wedge db=z*dt wedge dz and

```
a db+d(a*b)=3*z*t*dz+z^2*dt,
```

also regular. These are local formal controls only; they do not give
actual target polynomials with constant Jacobian or global completions.

### A genuine affine geometric survivor

There is a smooth affine surface containing the whole marked source plane,
with no source contraction, that realizes the type (p,r,n_E)=(1,2,3).
Set
```
N=Spec C[a,b,c,d]/(b^2-a*(c+1), b*c-a*d, b*d-c*(c+1)).
```
The map from the source is
```
a=x, b=x*y, c=x*y^2-1, d=x*y^3-y.
```
It identifies X=A2 with D(a) union D(c). Its complement is
E=(a,b,c), an affine line with coordinate d.
On D(c) the inverse source coordinates are x=a,y=d/c.
On D(c+1) put u=b/(c+1),t=d; the equations give
```
c=u*t, b=u*(1+u*t), a=u^2*(1+u*t).
```
Thus these two charts are respectively
C[x,y,(x*y^2-1)^-1] and C[u,t,(1+u*t)^-1].
They cover N and show it is smooth and integral. On their overlap
```
y=u^-1, x=u^2+u^3*t,
dx wedge dy=u*du wedge dt,
x*dy=-(1+u*t)*du.
```
The volume has order one on E and the same primitive is regular
globally. Values at E are v(y)=-1,v(x)=2; the local multiplicities
are [1,1,1] and n_E=3. This realizes the surviving numerical type
in genuine smooth affine marked-source geometry. A finite Keller
map N->A2 extending a constant-Jacobian source pair is NOT provided.
The example therefore illustrates the remaining target condition;
it does not contradict any of the necessary assertions above.

The checker independently controls the local extraction recurrences,
canonical labels, Euclidean sums, terminal multiplicities, the refined
index-two inequalities, the two Borisov determinant formulas and their
negative cases, actual-source contraction in the terminal monomial
case, the logarithmic residue exclusion, and the surviving formal
models. Finite extraction and parameter ranges are regression controls
for the written unbounded arguments, not a proof by exhaustion.

Assurance is an internal written proof with exact symbolic and finite
controls. No formal verification, external review, automatic integrality
of common y, exclusion of all ramified primes, or plane Jacobian
conjecture solution is claimed.
