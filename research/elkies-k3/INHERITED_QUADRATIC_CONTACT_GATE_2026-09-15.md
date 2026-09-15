# Actual quadratic halves from the Q80 basis do not give integral contact carriers

**Bounded constructive attempt, with exact retained points and an intersection
obstruction.** Apply the [quadratic contact constructor](QUADRATIC_CONTACT_INTERPOLATION_2026-09-15.md)
to traces T=P1 and T=P2 in the literal Q80 source basis (zero-based indices).
Both have height8 and an irreducible quadratic pole divisor. Test halves supplied
by O and by each of the34 signed basis sections.

There are exactly two nonpole quadratic contact fields in this signed-basis
bank, one for each trace. In both, exact group arithmetic gives T=2P6. Both
have h(T-2P6)=8 and force reducible members of the proposed contact linear
system. The pole contacts with half O are likewise degenerate. Thus this
specific bank supplies no integral arithmetic-genus-three genus-one bisection.
It does not exclude other inherited words, noninherited halves, other traces,
diagonal contacts or the [full correlated-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md).

## 1. An inherited-half height cutoff

Let T be primitive and pi:S -> Y the quotient by eta(P)=T-P. For a proposed
arithmetic-genus-three bisection B, its quotient section C has C^2=2.
An inherited section R has quotient image C_R, whose full inverse image is
R+(T-R). Translation invariance of height gives

```
h(T-2R)=4+2*R.(T-R),
C_R^2=h(T-2R)/2-4,
C.C_R=(C^2+C_R^2)/2=h(T-2R)/4-1.                 (1)
```

Suppose a quadratic closed contact point supplies both conjugate tangencies,
and its half values are the specializations of this same inherited R. The
image C_R has even contact with the branch curve at both points, because its
inverse image is a pair of sections meeting there. A candidate C also has
contact at least2 at each. Therefore C.C_R>=4 unless C=C_R. The latter has
split inverse image and cannot give an integral bisection. Necessarily

```
h(T-2R)>=20.                                    (2)
```

This argument is for two distinct geometric contacts. For a diagonal divisor
2P, an inherited value at P alone does not supply the higher jet needed to
replace two contacts; that case is not excluded by (2).

For the two retained height8 traces, the positive basis sections all have
distance height below20 except P16, whose distance height is20 for each trace.
The negative basis sections all exceed the cutoff. The complete signed bank
has32 cases below20 and36 at least20. The calculation below keeps both signs
when classifying actual equality points.

## 2. Exact contact search and its finite boundary

For each of the34 unsigned trace/basis pairs form the cancelled rational
function x(T)-x(2P_j). Its nonpole numerator has degree20,28,32 or36. The
checker removes coordinate-pole factors and handles common poles separately.
Only j=i has common poles, exactly the original quadratic poles of T; there
P_j=T=O. Such contacts coincide with the already degenerate half-O contacts.
There are no other common poles.

A fixed first panel of nine primes gave whole-polynomial exclusions for8 of
the34 numerators. The other26 were initially unresolved, not positive matches.
They were factored over Q. The final checker reconstructs each product exactly
and checks its factors at recorded good primes in5..997. For every factor of
degree greater than2 it verifies

```
gcd(f mod p, t^(p^2)-t)=1.
```

The monic factors have p-integral coefficients. The unfactored numerator
witnesses additionally retain a unit leading coefficient. Thus these checks
exclude all characteristic-zero factors of degree1 or2, without treating a
factorization status flag as an irreducibility certificate. Different factors
can use different primes.

Exactly two degree-two factors remain. Both are irreducible by an exact rational
discriminant nonsquare check. Their monic polynomials are

```
T=P1:
t^2
 + (70114129930275017326180596555487071413299835848182172685 /
    14939532517794632950341828404173725687768282987035627928)*t
 + 267127405946075878806222281610025522450828491696804729953 /
   49798441725982109834472761347245752292560943290118759760

T=P2:
t^2
 + (1048492732552684995784555571700541611396025659743719298181 /
    235735765721337163128964754397003895059571754897966013400)*t
 + 442864597984622010770679523559496838832423747875316582469 /
   94294306288534865251585901758801558023828701959186405360.
```

In each field, rational arithmetic in Q[t]/(q) verifies both coordinates of
T=2P6, and verifies T!=-2P6. The fibres are smooth by an exact discriminant gcd.
These are genuine halving points, not finite-field compatibility alone.
Their distance height for R=P6 is8; the opposite sign would have height40 but
does not give the half. Thus neither can pass the necessary condition(2).

Infinity is not a quadratic closed base point. This statement concerns
irreducible quadratic contact divisors, not a diagonal rational contact there.

## 3. Why the jet determinant vanishes without an integral member

The primitive height8 trace has quotient F0. Indeed the image of O is a section
of square0, forcing F0 among the possible F0,F1,F2 quotients. The image of P6
also has square0, since h(T-2P6)=8. Either image is a ruling section.

At either pair of retained quadratic contacts choose the other ruling coordinate
u so that this inherited image is u=0. The smooth branch curve has contact2
there, so u_1=u_2=v_1=v_2=0. The four-jet matrix for

```
A+B*t+C*u+D*t*u=0
```

has rank2 and forces A=B=0. Every member is u*(C+D*t)=0, hence reducible.
This interprets the previously retained rank2 interpolation control on actual
MW17 data. A vanishing norm determinant alone would have produced a false
candidate if the integral-section check were omitted.

The original quadratic pole contacts have the same description using the
image of O. Their failure does not depend on reconstructing quotient equations.

## Evidence and next boundary

The [input packet](../artifacts/generated-results/elkies-k3-inherited-quadratic-contact-v1/input.json)
binds the exact generic equation and17-section basis. The
[certificate](../artifacts/generated-results/elkies-k3-inherited-quadratic-contact-v1/result.json)
retains all68 signed heights and both actual quadratic halves. Initial scripts
and the unresolved-panel receipt are preserved under preflight. The signed-bank and factor producers
were bounded by30 CPU seconds and1GiB and completed in under one second locally.
No old bisection census or missing artifact was reconstructed.

The [replay](scripts/verify_inherited_quadratic_contacts.py) recomputes the
rational functions, factor products, modular witnesses, exact field signs and
height data. The geometric intersection implication remains a written proof,
not formal verification or an independently written whole-proof replay.

```
.venv/bin/python research/elkies-k3/scripts/verify_inherited_quadratic_contacts.py
```

A further inherited-half attempt must escape this fixed signed bank and, for
two distinct contacts, the height20 cutoff. That does not authorize a larger
coefficient search. A high-height word still needs an actual quadratic equality
field and successful integral contact interpolation; its height alone predicts
neither. No new gain or matching cover is supplied here.
