# Exclusion of global directional conductor degree twelve

**Theorem KellerDirectionalConductorTwelveExclusion.** For an actual
noninvertible complex plane Keller map in a source frame \((x,y)\) with
\(P\) monic of positive degree in \(x\), the exact global \(x\)-degree
of the directional coordinate conductor is even and at least fourteen.
The generic coordinate-order normalization defect is at least seven.
On arbitrary specialized lines the lower bounds remain two for
conductor degree and one for defect.

The [coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
relates the conductor degree to twice the normalization defect and
identifies the specialized order at every finite parameter. The
[degree-ten exclusion](GLOBAL_CONDUCTOR_TEN_EXCLUSION.md)
establishes the previous global threshold.

## Characteristic sequences of the exact order

**Proposition ConductorTwelveCharacteristicSequences.** If the
exact global conductor has degree twelve, the characteristic
sequence of the reduced generic coordinate order over
\(K=\mathbb C(y)\) is one of
\[
 (13,2),\ (7,3),\ (5,4),\ (6,4,9),\ (9,6,4),\ (10,4,5).
\]

Reduce the two degrees by triangular target changes over \(K\).
The order over an algebraic closure has affine-line normalization
and one place at infinity. The one-place equation meets the
polynomial degree hypotheses in
[Assi–García-Sánchez, Propositions 2 and 4 and §6](sources/assi-garcia-sanchez-1407.0490v1.pdf)
as in the
[cubic–quintic classification](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md#1-the-characteristic-sequences-of-the-exact-order).
The pole semigroup is free and symmetric. Its number of gaps is
the colength of the exact order, by longitudinal degree filtration.
Consequently the semigroup conductor is twelve.

Write the reduced generators of the semigroup as
\(r_0,r_1,\ldots,r_h\), with \(r_0>r_1>\gcd(r_0,r_1)\).
Let \(d_i=\gcd(r_0,\ldots,r_{i-1})\),
\(e_i=d_i/d_{i+1}\), and let \(\mu_i\) denote the conductor of
the degree sequence ending before \(r_i\) and divided by \(d_i\).
The \(d_i\) strictly decrease by divisibility. Freeness gives
\(e_i r_i\in\langle r_0,\ldots,r_{i-1}\rangle\) for \(i\ge2\).
The prefix conductors are even, with
\(\mu_2=(r_0/d_2-1)(r_1/d_2-1)\ge2\), and
\[
 d_{i+1}\mu_{i+1}-d_i\mu_i
  =(e_i-1)r_i-(d_i-d_{i+1})\ge0 .
\]
For \(h=1\) the equation \((r_0-1)(r_1-1)=12\) and coprimality
give \((13,2),(7,3),(5,4)\). For \(h\ge2\) use
\[
 12=d_h\mu_h+(d_h-1)(r_h-1).
\]
Since \(r_h\ge2\) and \(d_h\ge2\), this requires \(d_2\le5\).
If \(h\ge3\), divisibility forces \(d_2=4,d_3=2\) and
the bound forces \(\mu_2=2\), giving \(r_0=12,r_1=8\).
Freeness would require \(2r_2\in\langle12,8\rangle\).
Together with \(\gcd(4,r_2)=2\), this gives \(r_2\ge6\).
Then \(d_3\mu_3=8+r_2-2\ge12\), leaving no room for the final
positive term. Thus \(h\le2\).

For \(h=2\), \(d_2=2\) gives
\((\mu_2,r_2)=(2,9),(4,5)\), and \(d_2=3\) gives \((2,4)\).
For \(d_2=4,5\) the equation has no admissible integer solution.
This proves the list.

## Exclusion in all six sequences

The characteristic sequence is taken for the exact generic
order, so reductions made after a field extension must still
be tested on the original polynomial source. In the arguments
below, any selected generic pair with distinct finite limits
at a finite parameter represents an actual collision of the
specialized order. Monicity of the original \(P\) forces the
mate of a finite conductor root to stay finite. The exact
specialized conductor accounts for every such pair and every
degree drop.

| Reduced sequence | Exclusion |
|---|---|
| \((13,2)\) | The [all-degree exclusion of degree-two elements](DEGREE_TWO_ELEMENTS_AND_KELLER_REFLECTION.md) rules out the quadratic generator. |
| \((7,3)\) | Integrated cubic–septic parameter identities. Finite coefficient poles and nonhomogeneous ramification contradict the finite pairs, and the homogeneous specialization contradicts the fixed-order actions. |
| \((5,4)\) | Integrated quartic–quintic parameter identities and an exact obstruction to the alternate specialized embedding of the fixed order. |
| \((6,4,9)\) | A rational scale and a collision sextic control every finite pole of the quartic–sextic parameter equation. |
| \((9,6,4)\) | Exact square-root tails reduce the parameter identities to a cubic relation. Finite-pair estimates and the homogeneous reflected actions give the final obstruction. |
| \((10,4,5)\) | The polynomial part of the parameter Jacobian and the semigroup degrees are incompatible with a nonzero constant determinant. |

The [hyperelliptic and cubic–septic note](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md)
contains the cubic–septic proof and a separate period proof for \((13,2)\).
The third branch proof is the [quartic–quintic exclusion](CONDUCTOR_TWELVE_QUARTIC_QUINTIC_EXCLUSION.md).
The last three are in the
[cable exclusions](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md).
This exhausts the degree-twelve characteristic sequences.

The exact checkers accompanying these arguments replay the finite
coefficient, collision, and root-scale identities. The reduction to
the six sequences and the passage through finite parameter places
are written arguments. The
[conductor-fourteen exclusion](GLOBAL_CONDUCTOR_FOURTEEN_EXCLUSION.md)
and [conductor-sixteen exclusion](GLOBAL_CONDUCTOR_SIXTEEN_EXCLUSION.md)
raise the bound to eighteen. Automatic codifferent landing and the
unrestricted plane Jacobian conjecture remain open.
