# Exclusion of global directional conductor degree fourteen

**Theorem KellerDirectionalConductorFourteenExclusion.** For an
actual noninvertible complex plane Keller map in a source frame
\((x,y)\) with \(P\) monic of positive degree in \(x\), the exact
global \(x\)-degree of the directional coordinate conductor is
even and at least sixteen. The coordinate-order normalization
defect is at least eight on the generic line and on each line
where the global leading coefficient remains nonzero. On
arbitrary specialized lines the lower bounds are two for
conductor degree and one for normalization defect.

The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the specialized conductor on every line and shows
that its degree is twice the normalization defect. The
[previous global bound](GLOBAL_CONDUCTOR_TWELVE_EXCLUSION.md)
reduces this theorem to excluding degree fourteen.

## Exclusion of all ten embeddings

The [degree-fourteen classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md)
gives the ten necessary characteristic sequences of the actual
generic order over \(\mathbb C(y)\). The reduced generators are
obtained by elementary target changes over that field.
The approximate-root comparisons use the retained
[Assi–García-Sánchez input](sources/assi-garcia-sanchez-1407.0490v1.pdf)
in the one-place presentation of the exact order.

| Sequence | Proof |
|---|---|
| \((15,2),(15,6,2),(15,10,2)\) | The [all-degree exclusion of degree-two elements](DEGREE_TWO_ELEMENTS_AND_KELLER_REFLECTION.md) applies to the quadratic generator or last approximate root. |
| \((8,3)\) | [Cubic–octic parameter equations](CONDUCTOR_FOURTEEN_CUBIC_OCTIC_EXCLUSION.md) and the singular tangent determinant exclude the finite scale degeneration. |
| \((8,6,3),(12,8,3)\) | [Alternate embeddings of the cubic order](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md) contradict the norm identity and the constant parameter Jacobian, respectively. |
| \((12,8,6,3)\) | The [two successive approximate roots](CONDUCTOR_FOURTEEN_THREE_STAGE_EMBEDDING.md) force the constant term of the parameter Jacobian to vanish. |
| \((6,4,11)\) | [Quartic–sextic parameter equations](CONDUCTOR_FOURTEEN_CABLE_EXCLUSION.md) contradict finite source gap conditions. |
| \((9,6,5)\) | The [degree-nine cable proof](CONDUCTOR_FOURTEEN_DEGREE_NINE_CABLE_EXCLUSION.md) uses all pole and zero regimes of the parameter controlling the physical gaps. |
| \((10,4,7)\) | The [quartic–decic proof](CONDUCTOR_FOURTEEN_QUARTIC_DECIC_EXCLUSION.md) uses the exact normalized gap polynomial and the tangent determinant of its remaining fixed order. |

This exhausts degree fourteen. The finite-parameter arguments
use the exact specialized orders and the original monic
generator. They allow the specialized conductor degree to drop.

The three degree-two cases also have separate proofs. For \((15,2)\)
the [singular tangent determinant](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md)
controls the action scale. The other two cases have
[exact coefficient certificates](CONDUCTOR_FOURTEEN_HYPERELLIPTIC_EMBEDDINGS.md).

The checkers cited in the branch proofs replay finite
coefficient, norm, collision and algebraic elimination steps.
The remaining analytic continuation of algebraic branches and
valuation steps are in the written proofs. The
[conductor-sixteen exclusion](GLOBAL_CONDUCTOR_SIXTEEN_EXCLUSION.md)
raises the bound to eighteen. Automatic codifferent landing
and the unrestricted plane Jacobian conjecture remain open.

## The next threshold

**Proposition ConductorSixteenCharacteristicSequences.**
If the exact global directional conductor has degree sixteen
under the same hypotheses, then \((10,4,9)\) is the sole
case left by the approximate-root classification and the
[all-degree exclusion of degree-two elements](DEGREE_TWO_ELEMENTS_AND_KELLER_REFLECTION.md). The
[quartic–decic proof](CONDUCTOR_SIXTEEN_QUARTIC_DECIC_EXCLUSION.md)
excludes that case by its parameter equations.

Use the notation \(r_i,d_i,\mu_i\) of the
[degree-fourteen classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md).
With one gcd step, the equation \((r_0-1)(r_1-1)=16\) and
coprimality give \((17,2)\), excluded because its order contains
a quadratic generator.

With two steps, \(d=d_2\le7\) and
\(16=d\mu_2+(d-1)(r_2-1)\). At \(d=2\) the pairs
\((\mu_2,r_2)\) are \((2,13),(4,9),(6,5)\). In the first
case the prefix is \((6,4)\), so the strict inequality
\(r_1d_1>r_2d_2\) of Assi–García-Sánchez Proposition 2(iii)
would require \(24>26\).
In the last case the possible prefixes are \((14,4),(8,6)\),
and neither generates \(2r_2=10\). The middle case gives
\((10,4,9)\). At \(d=3\) both candidates \(r_2=6,3\) violate
the final gcd one. For \(d=4,5,6,7\), the admissible
even value \(\mu_2=2\) gives no integer \(r_2\).

For three or more steps, \(d_2\le7\) allows precisely three
and \((d_2,d_3)=(4,2),(6,2),(6,3)\). At \(d_2=6\) the
prefix is \((18,12)\). Freeness requires \(e_2r_2\ge12\).
The increment after the prefix is at least four if \(d_3=2\).
If \(d_3=3\), the gcd further gives \(r_2\ge9\) and the
increment is at least six. In both cases the last positive
increment makes the conductor exceed sixteen. At \(d_2=4\)
the prefix is \((12,8)\), and \(r_2+r_3=11\). Freeness gives
\(r_2\ge6\) with \(r_2\equiv2\bmod4\), so the sole possibility
with \(r_3\ge2\) is \((r_2,r_3)=(6,5)\). It would require
\(2r_3=10\in\langle12,8,6\rangle\), which fails.
