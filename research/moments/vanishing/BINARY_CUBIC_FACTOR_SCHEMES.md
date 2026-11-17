# Binary cubic and three-factor Gaussian moment schemes

Claim: `G2CUBICSCHEME20261114`. Work over C. The proof combines written
invariant/torus arguments, exact Smith forms, and bounded Singular
calculations. There is no independent replay, external review or formal
verification. The full factor decomposition below is constructive; a
direct expanded45-way CAS intersection is not asserted.

Choose Gaussian null coordinates Z,W with E(ZW)=1 and E(Z^2)=E(W^2)=0.
For example Z=(X+iY)/sqrt(2), W=(X-iY)/sqrt(2) for independent standard
real Gaussians X,Y. This fixed invertible complex change preserves the
coefficient scheme. Write

    P=(aZ+dW)(bZ+eW)(cZ+fW)
     =xW^3+yZW^2+zZ^2W+wZ^3,

    x=def, y=aef+dbf+dec, z=abf+aec+dbc, w=abc.

Let J denote the full pure-moment ideal in C[x,y,z,w], and I its
three-factor pullback to S=C[a,b,c,d,e,f]. Then

    J=(xw+yz, xz^3+y^3w-3y^2z^2, y^3z^3),              (1)
    I=J*S.                                              (2)

These are full all-order ideals, not just finite-moment approximations.
Their sharp generation cutoff6 agrees with the earlier
[dimension-wise coefficient theorem](THREE_FACTOR_COEFFICIENT_PULLBACK.md).

## Results

| Scheme | Minimal primes | All associated primes / primary components | Degree | Nilradical index |
|---|---:|---:|---:|---:|
| Binary-cubic coefficient ideal J, four coefficients |2|5|2|6|
| Three-factor coefficient ideal I, six coefficients |9|45|36|9|

The radical of I is the squarefree monomial ideal

    B=(bcdf,bcde,acef,acde,abef,abdf).

Its nine minimal primes are

    (a,d),(b,e),(c,f),
    (a,b),(a,c),(b,c),(d,e),(d,f),(e,f).                 (3)

The first three are zero-factor planes, with generic local length8.
The other six describe two factors on the same isotropic line, with
generic length2. Thus the reduced degree is9 and the full degree is
3*8+6*2=36.

There are36 embedded primes:18 coordinate primes and18 noncoordinate
primes. The latter are

    (a,d,bf-lambda*ce),                                 (4)

for the six distinct roots lambda of

    D(t)=t^6+4t^5+14t^4+18t^3+14t^2+4t+1,

and the two other choices of the zero factor. All lambda are nonzero.
The intrinsic local torsion at each prime in (4) has length4.

The intrinsic coefficient-origin torsion has Hilbert series
3t^28+6t^29+4t^30, hence length13. The other embedded coordinate types
and their intrinsic local torsion lengths appear in the component table
below. They are not lengths of the chosen primary quotients.

## 1. All-order generation in four binary-cubic coefficients

For even m=2k, E(P^m)=(3k)! times the central coefficient of
(x+yT+zT^2+wT^3)^(2k). Odd pure moments vanish by homogeneity.
The coefficient polynomial is invariant under weights(-3,-1,1,3) and
under reversal (x,y,z,w) to (w,z,y,x).

The weight-zero ring is

    C[s,t,r,q]/(rq-st^3),
    s=xw, t=yz, r=xz^3, q=y^3w.

Indeed, cancel min(exponent(x),exponent(w)) copies of xw in any balanced
monomial. If x remains, balance forces the remaining monomial to be
r^A*t^B; if w remains it is q^A*t^B. The displayed relation is the only
relation, also seen from this normal form.

Reversal exchanges r,q and fixes s,t. Its invariant ring is the
polynomial ring C[s,t,u=r+q]: after q=u-r the relation is the monic
quadratic r^2-u*r+st^3, with basis1,r. Invariance of a+b*r forces b=0
in characteristic zero.

The first three central-coefficient polynomials are

    2(s+t),
    6s^2+24st+6t^2+4u,
    20s^3+180s^2t+180st^2+60su+20t^3+30tu.

The first two give s=-t,u=3t^2; the third becomes -90t^3. Thus their
ideal in the symmetric invariant ring is (s+t,u-3t^2,t^3).
Every subsequent central coefficient has weighted degree k, for weights
wt(s)=wt(t)=1, wt(u)=2, and reduces to a multiple of t^k. All k>=3
therefore already belong to this ideal. Pullback gives exactly (1), and
then (2). No long prefix is extrapolated.

## 2. The five-component four-coefficient scheme

The minimal primary ideals of J are (x,y) and (z,w), both reduced.
There are also the x-axis and w-axis primary ideals and one origin-primary
ideal. Their explicit standard-basis generators are retained in the
[certificate](binary_factor_scheme_certificate.json), in the
`binary_cubic.primary_components` array.

The exact five-way intersection equals J in both directions. The axis
ideals have the asserted coordinate radicals: twelfth powers of the
three transverse variables belong to each, and the axis parameter is a
nonzerodivisor. Because these ideals are homogeneous, this excludes an
embedded origin and proves primaryness. The chosen origin ideal contains
twelfth powers of all four variables and is primary. Deleting any one
component changes the intersection, so the decomposition is irredundant.

The complete monomial tests give B_J^6 contained in J but B_J^5 not
contained in J for B_J=(xz,xw,yz,yw). The intrinsic origin torsion is
one-dimensional in degree8. This small scheme is distinct from its
nonflat factor-coefficient pullback, which has45 components.

## 3. The zero-factor chart and its six exceptional parameters

On the chart b,c,e,f nonzero, rescale factors B,C and apply the Gaussian
torus action Z to rho*Z, W to rho^-1*W. After a finite étale square-root
extension one can set b=c=e=1 and f=t, where

    t=bf/(ce) is invertible.

The two remaining A coefficients are rescaled by units. The normalized
coefficient map is

    x=t*d, y=t*a+(t+1)d, z=(t+1)a+d, w=a.

In the resulting ideal L in C[t,t^-1][a,d], the three generators have
A-degrees2,4,6. Construct every monomial cofactor in each A-degree0,...,7
and take the exact relation-matrix Smith form over Q[t]. After inverting
t, the complete quotient, as a graded base module, is

    C[t,t^-1]^8
      plus C[t,t^-1]/(D), in degree4,
      plus two copies of C[t,t^-1]/(D), in degree5,
      plus C[t,t^-1]/(D), in degree6.                   (5)

The free ranks in degrees0,...,4 are1,2,2,2,1. Degree7 is zero; since
a,d have degree1, every later degree is zero as well. The computed
resultant of the first quadratic and quartic is t^2*D(t)^2, but (5),
not the resultant alone, establishes the torsion and flatness failure.
D is squarefree and D(0)=1.

The D-torsion is killed by D itself. Consequently

    L:D^infinity=L:D,
    L=(L:D^infinity) intersect (L+(D)).                 (6)

For the second identity, if g=i+D*h belongs to the saturation, saturation
also implies h belongs to it, and then D*h belongs to L.
The quotient by the saturation is free of rank8 over C[t,t^-1], with
nilpotent positive-degree ideal (a,d). Multiplication by an element with
nonzero base constant term has nonzero determinant, proving that this
saturation is (a,d)-primary.

Chinese remainders split L+(D) into the six ideals L+(t-lambda). Each is
primary at (a,d,t-lambda), of colength12, while its intrinsic embedded
torsion is only4. Thus this open chart has a full seven-component
decomposition: one minimal component and six embedded point components.
The normalization is a faithfully flat étale monomial change, so the
primary statements descend. Contracting back gives exactly (4).

## 4. Why the associated-prime list is exhaustive

The independent factor scalings and the Gaussian torus form a connected
four-dimensional torus. Its coefficient weights are the columns of

    [1 0 0 1 0 0]
    [0 1 0 0 1 0]
    [0 0 1 0 0 1]
    [1 1 1 -1 -1 -1].

Every associated prime is invariant under this torus, since the
associated-prime set is finite and the torus is connected.

The nine coordinate four-faces complementary to (3) cover the radical.
For a zero-factor face, its open coordinate torus is the chart of (5),
which yields exactly the one minimal prime and its six exceptional
primes. For an isotropic-pair face, its four weights are independent.
A finite étale monomial normalization makes all four free coefficients1;
the transverse quotient is an Artinian local algebra of length2. It has
only its minimal associated prime on that open set.

The remaining boundary faces have at most three free coordinates. Every
set of at most three coefficient weights is linearly independent, so the
torus is transitive on each open coordinate face. An invariant prime
there must be the coordinate prime. There are no hidden noncoordinate
primes left on these boundaries.

Factor permutations and simultaneous Z/W exchange reduce the coordinate
candidates to the following ten exact parameter-field socle tests.
Free variables are put in the coefficient field, and the maximal ideal
m is generated by the others. The complete quotient (I:m)/I detects
whether the candidate is associated. Its dimension need not equal the
number of redundant quotient generators.

| Free-coordinate representative | Orbit size | Socle dimension |
|---|---:|---:|
| b,c,e,f: zero-factor minimal |3|1|
| c,d,e,f: isotropic-pair minimal |6|1|
| a,b,d |12|0|
| a,b,c |2|1|
| a,b,f |6|0|
| a,d |3|5|
| a,b |6|2|
| a,e |6|0|
| a |6|2|
| none |1|4|

Thus there are27 coordinate associated primes, including the nine
minimal primes, and18 noncoordinate associated primes from (4):45 in
total over C. The quotient algorithm is the completeness input for these
socle tests; forward membership checks alone would not prove absence.

## 5. Explicit primary components for all45 primes

For a zero-factor component A, put H=bcef and

    Delta=U^6+4U^5V+14U^4V^2+18U^3V^3+14U^2V^4+4UV^5+V^6,
    U=bf, V=ce.

Its isolated primary ideal is I:(H*Delta)^infinity. Permuting factors
gives three components. For the isotropic-pair prime (a,b), use
I:(cdef)^infinity; factor permutations and Z/W exchange give six.
The length8/2 primaryness arguments are the normalized chart arguments
above followed by contraction.

For every root lambda of D, a primary ideal for (4) is

    (I+(bf-lambda*ce)):(bcef)^infinity.                 (7)

Together with its factor permutations these give18 components, by (6).

For an embedded coordinate prime, let F be its free-coordinate set.
Define Q(F;L)=(I+L):(product of variables in F)^infinity; the empty
product is1. The following table gives all18 remaining components by
symmetry. Weights are on the listed remaining variables after the free
variables are normalized to1.

| F | Copies | Remaining variables / weights | Intrinsic torsion length | Top weight | Pure-power ideal L |
|---|---:|---|---:|---:|---|
| abc |2|def / 1,1,1 |84|13|(d^14,e^14,f^14)|
| ad |3|bcef / 1,1,1,1 |62|18|(b^19,c^19,e^19,f^19)|
| ab |6|cdef / 1,1,1,2 |8|21|(c^22,d^22,e^22,f^11)|
| a |6|bcdef / 1,1,1,2,2 |2|32|(b^33,c^33,d^33,e^17,f^17)|
| empty |1|abcdef / all1 |13|30|(a^31,b^31,c^31,d^31,e^31,f^31)|

Every listed L is generated in weights strictly above the top intrinsic
torsion weight. These are sufficient explicit choices, not optimal
truncation exponents. The resulting ideal is primary at the coordinate
prime complementary to F: on the free-variable chart, torus normalization
makes its quotient a constant Artinian local algebra over a Laurent
base, after a finite étale extension if needed. Primaryness descends and
survives contraction.

The exact local intrinsic torsion Hilbert series used here are:

* F=abc: coefficients1,3,6,9,11,12,12,11,9,6,3,1 in degrees2,...,13.
* F=ad: t^12+8t^13+14t^14+16t^15+14t^16+8t^17+t^18.
* F=ab: t^18+2t^19+3t^20+2t^21 in the displayed weighting.
* F=a: 2t^32 in the displayed weighting.
* F empty: 3t^28+6t^29+4t^30.

They come from complete coordinate-saturation Hilbert differences.
The highest weights independently agree with those of the local socles,
since the top homogeneous part of finite-length torsion is in the socle.

## 6. Proof of the full intersection

Let K be the intersection of the9 isolated,18 noncoordinate, and18
embedded coordinate primary ideals just defined. Each contains I.
At every minimal prime, its isolated component localizes to I. At every
noncoordinate embedded prime, the isolated zero-factor component and
(7) already recover I by (6). Thus (K/I) has zero localization at all
these primes.

Suppose K/I is nonzero and choose an associated prime P of minimal height.
It must be one of the embedded coordinate primes. Localizing at P
leaves no strictly smaller associated prime of K/I, by minimality, so
(K/I)_P has finite length. It is a submodule of the intrinsic local
torsion T_P of (S/I)_P.

On the free-variable chart for P, T_P is the localization of the graded
finite-length torsion tabulated above. If its top weight is B, its
intersection with any ideal generated in weights greater than B is zero:
their homogeneous supports are disjoint. Exact localization preserves
this zero intersection. But K is contained in the selected Q(F;L), whose
image is precisely such a high-weight ideal. Hence (K/I)_P=0, a
contradiction. Therefore K=I.

There is one component for each of the45 distinct known associated
primes, so the decomposition is irredundant. This is the full mathematical
intersection proof; the timed-out general decomposition algorithm is not
used for it.

## 7. Nilpotence and evidence

The nine coordinate primes in (3) are the edges of a triangular-prism
graph. Its six minimal vertex covers give exactly the quartic radical
generators B displayed above. The three moment generators lie in B.
An exhaustive monomial-ideal reduction gives B^9 contained in I, checking
all1540 minimal generators. The element

    (acef)^8

has a nonzero standard remainder modulo I. Thus B^8 is not contained in I
and the exact global nilradical index is9. The full Hilbert numerator
gives degree36, agreeing with the independently computed generic lengths.

The [checker](verify_binary_factor_scheme.py),
[certificate](binary_factor_scheme_certificate.json), and
[receipt](binary_factor_scheme_verification.json) bind the moment
coefficient map, all-order invariant reduction, five-component auxiliary
decomposition and its primary/irredundancy audits, Smith forms, complete
coordinate candidate list, socles, torsion Hilbert series and monomial
power reductions.

```sh
.venv/bin/python research/moments/vanishing/verify_binary_factor_scheme.py
```

Use `--audit-existing` for the retained-output audit and symbolic checks.
The accepted Singular runs each have50-second and768MiB limits and all
complete in under4 seconds individually in the recorded environment.
The full six-variable `primdecGTZ` attempt timed out after50 seconds;
its completed standard-basis/Hilbert checkpoint, inputs and terminal
receipt are retained. No45-component output is attributed to that run.
Benign Singular redefinition comments are removed only when parsing;
original logs and their hashes are preserved, and actual error markers
cause the audit to fail.

This closes the full N=2 coefficient scheme for three centered factors,
and the full four-coefficient homogeneous binary-cubic moment scheme.
The remaining actual-coefficient dimensions in this branch are N=3,4.
Affine factors, changed characteristic, arbitrary coefficient subfamilies
and general ternary cubic GMC retain their separate obligations.
