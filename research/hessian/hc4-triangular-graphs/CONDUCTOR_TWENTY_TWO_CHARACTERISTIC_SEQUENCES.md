# The reduced characteristic sequences at conductor twenty-two

**Proposition ConductorTwentyTwoCharacteristicSequences.** Suppose
the exact global directional coordinate conductor of a complex plane
Keller map in a monic longitudinal source frame has degree twenty-two.
The reduced generic characteristic sequence belongs to the following
list.

| Number of gcd steps | First gcd \(d_2\) | Sequences |
|---|---:|---|
| One | \(1\) | \((23,2)\) |
| Two | \(2\) | \((8,6,11),(10,4,15),(14,4,11)\) |
| Two | \(5\) | \((15,10,4)\) |
| Three | \(4\) | \((12,8,6,11)\) |

For the [one-place sequence conditions in Proposition 2](sources/assi-garcia-sanchez-1407.0490v1.pdf)
use the [exhaustive arithmetic at degree twenty](CONDUCTOR_TWENTY_CHARACTERISTIC_SEQUENCES.md#exhaustive-arithmetic)
with conductor twenty-two. The bound on the initial gcd is now eleven,
and the same positive increments, proper gcd drops, strict inequality
and freeness apply. The [exact enumeration](verify_conductor_twenty_characteristic_sequences.py)
replays all six cases.

## Reduction to four embeddings

The [degree-two target-coordinate theorem](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md)
excludes \((23,2)\). The
[short three-to-two Laurent theorem](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md)
excludes \((15,10,4)\) with \(m=5,\ r=4\). The four remaining
sequences are \((8,6,11),(10,4,15),(14,4,11),(12,8,6,11)\).
