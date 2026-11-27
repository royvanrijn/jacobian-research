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
| \((15,2)\) | [Singular tangent determinant](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md) forces the action scale to be constant through every finite specialization. |
| \((15,6,2),(15,10,2)\) | [Degree-two approximate roots](CONDUCTOR_FOURTEEN_HYPERELLIPTIC_EMBEDDINGS.md) give a composition identity whose exact coefficient ideals contain \(1\) in every scaling chart. |
| \((8,3)\) | [Cubic–octic parameter equations](CONDUCTOR_FOURTEEN_CUBIC_OCTIC_EXCLUSION.md) and the singular tangent determinant exclude the finite scale degeneration. |
| \((8,6,3),(12,8,3)\) | [Alternate embeddings of the cubic order](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md) contradict the norm identity and the constant parameter Jacobian, respectively. |
| \((12,8,6,3)\) | The [two successive approximate roots](CONDUCTOR_FOURTEEN_THREE_STAGE_EMBEDDING.md) force the constant term of the parameter Jacobian to vanish. |
| \((6,4,11)\) | [Quartic–sextic parameter equations](CONDUCTOR_FOURTEEN_CABLE_EXCLUSION.md) contradict finite source gap conditions. |
| \((9,6,5)\) | The [degree-nine cable proof](CONDUCTOR_FOURTEEN_DEGREE_NINE_CABLE_EXCLUSION.md) uses all pole and zero regimes of the parameter controlling the physical gaps. |
| \((10,4,7)\) | The [quartic–decic proof](CONDUCTOR_FOURTEEN_QUARTIC_DECIC_EXCLUSION.md) uses the exact normalized gap polynomial and the tangent determinant of its remaining fixed order. |

This exhausts degree fourteen. The finite-parameter arguments
use the exact specialized orders and the original monic
generator. They allow the specialized conductor degree to drop.

The checkers cited in the branch proofs replay finite
coefficient, norm, collision and algebraic elimination steps.
The remaining analytic continuation of algebraic branches and
valuation steps are in the written proofs. Global conductor
degree sixteen and higher, automatic codifferent landing, and
the unrestricted plane Jacobian conjecture remain open.
