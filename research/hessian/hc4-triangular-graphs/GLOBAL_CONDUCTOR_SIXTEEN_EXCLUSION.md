# Exclusion of global directional conductor degree sixteen

**Theorem KellerDirectionalConductorSixteenExclusion.** For an
actual noninvertible complex plane Keller map in a source frame
\((x,y)\) with \(P\) monic of positive degree in \(x\), the exact
global \(x\)-degree of the directional coordinate conductor is
even and at least eighteen. The coordinate-order normalization
defect is at least nine on the generic line and on each line
where the global leading conductor coefficient remains nonzero.
On arbitrary specialized lines the lower bounds are two for
conductor degree and one for normalization defect.

The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the specialized conductor on every line. Its degree
is twice the normalization defect. The
[degree-fourteen lower bound and the degree-sixteen classification](GLOBAL_CONDUCTOR_FOURTEEN_EXCLUSION.md)
reduce the theorem to the reduced sequence \((10,4,9)\).
The [quartic–decic exclusion](CONDUCTOR_SIXTEEN_QUARTIC_DECIC_EXCLUSION.md)
rules it out. Its parameter invariant makes the scaling rational,
and the two finite degeneration analyses force the remaining
parameter to be a nonzero constant, contrary to the Keller
identity. This proves the new bound.

## Characteristic sequences at degree eighteen

**Proposition ConductorEighteenCharacteristicSequences.** If
the exact global directional conductor has degree eighteen in
the same monic source frame, its reduced characteristic sequence
belongs to the following necessary list:
\[
\begin{gathered}
 (19,2),\ (10,3),\ (7,4),\ (10,4,11),\ (14,4,7),\ (8,6,7),\\
 (10,6,3),\ (9,6,7),\ (15,6,4),\ (15,10,3),\ (12,8,6,7).
\end{gathered}
\]
The [all-degree exclusion of degree-two elements](DEGREE_TWO_ELEMENTS_AND_KELLER_REFLECTION.md)
excludes \((19,2)\). The
[short three-to-two cable theorem](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md)
excludes \((15,10,3)\).

Use the notation of the
[degree-fourteen classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md).
For one gcd step the equation \((r_0-1)(r_1-1)=18\)
and coprimality give \((19,2),(10,3),(7,4)\).
For two steps, \(d=d_2\le8\) and
\(18=d\mu_2+(d-1)(r_2-1)\). The following candidates also
satisfy \(d r_2\in\langle r_0,r_1\rangle\) and the strict
approximate-root inequality \(r_1d_1>r_2d_2\).

| \(d\) | \((\mu_2,r_2)\) | \((r_0,r_1,r_2)\) |
|---:|---|---|
| 2 | \((4,11),(6,7),(8,3)\) | \((10,4,11),(14,4,7),(8,6,7),(10,6,3)\) |
| 3 | \((2,7),(4,4)\) | \((9,6,7),(15,6,4)\) |
| 5 | \((2,3)\) | \((15,10,3)\) |

For \(d=2,\mu_2=2\), the omitted \((6,4,15)\) fails the strict
inequality. At \(\mu_2=8\), the other prefix \((18,4)\) does not
generate \(2r_2=6\). At \(d=4,6,7,8\) there is no integral \(r_2\)
with the conductor and gcd requirements. These arithmetic
conditions are those of
[Assi–García-Sánchez, Proposition 2](sources/assi-garcia-sanchez-1407.0490v1.pdf).

For three or more steps, \(d_2\le8\) leaves \(d_2=4,6,8\).
If \(d_2=8\), then \(d_3=2\) or \(4\) before the final step,
and freeness requires \(e_2r_2\ge r_1\ge16\).
The increment from the prefix of conductor at least sixteen
is at least six or four, respectively. Both are too large.
If \(d_2=6\), the prefix is \((18,12)\). For \(d_3=3\)
freeness and the gcd require \(r_2\ge9\), giving an
increment of at least six before the last positive contribution.
For \(d_3=2\), freeness gives \(r_2\ge4\) and the conductor
bound forces equality. The final entry is then \(r_3=3\),
which fails \(2r_3\in\langle18,12,4\rangle\).

If \(d_2=4\), then \(d_3=2\) and \(h=3\). Freeness and the gcd
force \(r_2\ge6\), so \(\mu_2=4\) would already use at least
twenty before the last positive increment. Thus \(\mu_2=2\),
the prefix is \((12,8)\), and \(r_2+r_3=13\). The only
possibilities with \(r_2\equiv2\bmod4\), \(r_2\ge6\) and
\(r_3\ge2\) are \((6,7)\) and \((10,3)\). The latter fails
\(2r_3\in\langle12,8,10\rangle\). The first gives the listed
three-stage sequence.

## Reduction to four embeddings

The [cubic norm identities](CONDUCTOR_EIGHTEEN_CUBIC_NORM_EXCLUSION.md)
exclude \((10,6,3)\). The
[three-stage parameter equation](CONDUCTOR_FOURTEEN_THREE_STAGE_EMBEDDING.md)
excludes \((12,8,6,7)\), and the
[filtered parameter comparisons](CONDUCTOR_EIGHTEEN_FILTERED_PARAMETER_EXCLUSIONS.md)
exclude \((14,4,7),(15,6,4),(8,6,7)\).
These exclusions reduce degree eighteen to
\((10,3),(7,4),(10,4,11),(9,6,7)\). The
[degree-eighteen theorem](GLOBAL_CONDUCTOR_EIGHTEEN_EXCLUSION.md)
excludes these four embeddings and gives the lower bound twenty. The
[degree-twenty theorem](GLOBAL_CONDUCTOR_TWENTY_EXCLUSION.md) gives the
current global bound twenty-two.
