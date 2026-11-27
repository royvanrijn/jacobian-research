# Exclusion of global directional conductor degree twenty

**Theorem KellerDirectionalConductorTwentyExclusion.** Let
\((P,Q)\) be an actual noninvertible complex plane Keller map
in a source frame \((x,y)\) with \(P\) monic of positive degree
in \(x\). Its exact global directional coordinate conductor has
even \(x\)-degree at least twenty-two. The coordinate-order
normalization defect is at least eleven on the generic line and
on every line where the global leading conductor coefficient is
nonzero. On arbitrary specialized lines the lower bounds are
two for conductor degree and one for normalization defect.

## Exclusion of the four remaining embeddings

The [degree-eighteen theorem](GLOBAL_CONDUCTOR_EIGHTEEN_EXCLUSION.md)
gives even global conductor degree at least twenty. The
[degree-twenty classification and reductions](CONDUCTOR_TWENTY_CHARACTERISTIC_SEQUENCES.md)
leave four sequences at degree twenty. Their individual
exclusions use the polynomial parameter identities for the actual
Keller generators.

| Sequence | Branch proof | Limiting condition |
|---|---|---|
| \((6,5)\) | [Quintic–sextic](CONDUCTOR_TWENTY_QUINTIC_SEXTIC_EXCLUSION.md) | The leading invariants yield immersed pole shapes. The homogeneous case gives a fixed singular order. |
| \((11,3)\) | [Cubic–degree eleven](CONDUCTOR_TWENTY_CUBIC_UNDECIC_EXCLUSION.md) | A leading resultant controls finite poles. A nonzero invariant bounds the ramification at infinity and rules out a pole there. |
| \((10,4,13)\) | [Quartic–decic](CONDUCTOR_TWENTY_QUARTIC_DECIC_EXCLUSION.md) | The approximate-root identities control all coefficient poles. The scale has Kummer degree at most seven. |
| \((9,6,8)\) | [Sextic–nonic](CONDUCTOR_TWENTY_SEXTIC_NONIC_EXCLUSION.md) | The square-root coefficient equations isolate the poles and the removable chart. |

Consequently the global degree is at least twenty-two.
The [exact specialization theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the specialized degree with twice the normalization
defect at every parameter. On lines where the leading conductor
coefficient is nonzero, the specialized degree equals the global
degree. The same theorem gives the arbitrary-line bounds.

See the [degree-twenty-two characteristic sequences](CONDUCTOR_TWENTY_TWO_CHARACTERISTIC_SEQUENCES.md)
for the next threshold. Higher conductors and automatic
ordinary trace-dual landing remain open.
