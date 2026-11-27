# The reduced characteristic sequences at conductor twenty

**Proposition ConductorTwentyCharacteristicSequences.** Suppose the
exact global directional coordinate conductor of a complex plane Keller
map in a monic longitudinal source frame has degree twenty. The reduced
generic characteristic sequence belongs to the following list.

| Number of gcd steps | First gcd \(d_2\) | Sequences |
|---|---:|---|
| One | \(1\) | \((6,5),(11,3),(21,2)\) |
| Two | \(2\) | \((8,6,9),(10,4,13),(10,6,5),(14,4,9)\) |
| Two | \(3\) | \((9,6,8),(15,6,5),(21,6,2)\) |
| Two | \(4\) | \((12,8,5)\) |
| Two | \(7\) | \((21,14,2)\) |
| Three | \(4\) | \((12,8,6,9),(12,8,10,5)\) |

## Exhaustive arithmetic

The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives a finite birational plane order with normalization \(K[x]\) for
\(K=\mathbb C(y)\). For the one-place characteristic-sequence conditions
use [Assi–García-Sánchez, Proposition 2](sources/assi-garcia-sanchez-1407.0490v1.pdf)
as in the [degree-fourteen classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md).
Write \(r_0=du,\ r_1=dv\) with \(u>v\ge2\) coprime and \(d=d_2\).
The part of the conductor accounted for by this prefix is
\(c_1=d(u-1)(v-1)\). If the current gcd is \(d_i\), appending \(r_i\)
with \(d_{i+1}=\gcd(d_i,r_i)<d_i\) changes this quantity by
\[
 c_i-c_{i-1}=(d_i/d_{i+1}-1)(r_i-d_{i+1}).
\]
Indeed the sum of the gcd differences telescopes, so the prefix
quantity equals the standard conductor expression
\(\sum_i (d_i/d_{i+1}-1)r_i-r_0+1\) at the final gcd one.
Freeness requires \((d_i/d_{i+1})r_i\in
\langle r_0,\ldots,r_{i-1}\rangle\). It implies \(r_i>d_{i+1}\)
by dividing by \(d_i\), because every previous entry divided by \(d_i\)
is at least two. Thus every increment is positive.

It follows that \(d\le10\) and \(d(u-1)(v-1)\le20\).
Once a prefix of accumulated cost \(c\) has been chosen, each proper
divisor \(d'\) of its current gcd \(d_i\) gives the finite range
\[
 2\le r_i\le d'+\left\lfloor\frac{20-c}{d_i/d'-1}\right\rfloor .
\]
In this range retain the required gcd, the freeness condition and the
strict inequality \(r_{i-1}d_{i-1}>r_i d_i\). Retain a terminating row
only when the cost is exactly twenty. The
[exact enumeration](verify_conductor_twenty_characteristic_sequences.py)
checks this procedure and recovers the degree-eighteen classification
as a control.

## Exclusion of ten embeddings

The [degree-two target-coordinate theorem](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md)
excludes \((21,2)\). The
[short three-to-two Laurent theorem](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md)
excludes \((21,14,2)\) with \(m=7,\ r=2\).

The [filtered parameter exclusions](CONDUCTOR_TWENTY_FILTERED_PARAMETER_EXCLUSIONS.md)
rule out \((14,4,9),(15,6,5),(21,6,2),(8,6,9),(10,6,5)\).
The [common-centroid argument](CONDUCTOR_FOURTEEN_THREE_STAGE_EMBEDDING.md)
excludes \((12,8,6,9)\). The
[degree-five approximate-root argument](CONDUCTOR_TWENTY_TWELVE_EIGHT_FIVE_EXCLUSION.md)
excludes \((12,8,5)\). The
[three-stage Laurent comparison](CONDUCTOR_TWENTY_TWELVE_EIGHT_TEN_FIVE_EXCLUSION.md)
excludes \((12,8,10,5)\). This leaves \((6,5),(11,3),(10,4,13),(9,6,8)\)
for the [four branch proofs completing the degree-twenty theorem](GLOBAL_CONDUCTOR_TWENTY_EXCLUSION.md).
