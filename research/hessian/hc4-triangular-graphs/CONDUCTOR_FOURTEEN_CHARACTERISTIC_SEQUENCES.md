# The characteristic sequences at conductor degree fourteen

Use the original hypotheses of the
[directional conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md).
Thus \(P,Q\) form an actual noninvertible complex polynomial plane
Keller map, \(P\) is monic of positive longitudinal degree, and the
generic order \(K[P,Q]\) in \(K[x]\), for \(K=\mathbb C(y)\), has
exact conductor of degree fourteen. Reduce the two generator
degrees by triangular target operations over \(K\).

**Proposition ConductorFourteenCharacteristicSequences.**
The possible reduced characteristic sequences are among
\[
 (15,2),\ (8,3),\ (6,4,11),\ (10,4,7),\ (8,6,3),\
 (9,6,5),\ (15,6,2),\ (12,8,3),\ (15,10,2),\ (12,8,6,3).
\]
This is a necessary list for the embedding in the original
coordinate order. Different characteristic sequences can
generate the same numerical semigroup.

As in the [degree-twelve classification](GLOBAL_CONDUCTOR_TWELVE_EXCLUSION.md#characteristic-sequences-of-the-exact-order),
the reduced plane relation is monic with the required coefficient
degree bounds and the exact normalization has one place at
infinity. Apply [Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf).
The number of semigroup gaps is the colength of the actual
coordinate order, namely seven. By symmetry the semigroup
conductor is fourteen.

Let \(r_0,\ldots,r_h\) be the characteristic sequence,
\(d_i=\gcd(r_0,\ldots,r_{i-1})\), \(e_i=d_i/d_{i+1}\), and
\(\mu_i\) the numerical semigroup conductor of the prefix
normalized by \(d_i\). For \(h=1\) the product
\((r_0-1)(r_1-1)=14\) gives \((15,2)\) and \((8,3)\).
Otherwise \(\mu_2=(r_0/d_2-1)(r_1/d_2-1)\ge2\) is even.
Every positive-degree entry is at least two since the order
has nonzero colength. The prefix recurrence is
\[
 d_{i+1}\mu_{i+1}=d_i\mu_i+(e_i-1)r_i-(d_i-d_{i+1}).
\]
Its increments are nonnegative, and the final increment gives
\(14=d_h\mu_h+(d_h-1)(r_h-1)\). In particular \(2d_2+1\le14\),
so \(d_2\le6\).

If \(h=2\), use \(14=d_2\mu_2+(d_2-1)(r_2-1)\).
The table gives the candidates satisfying the equation and the
freeness condition \(d_2r_2\in\langle r_0,r_1\rangle\).

| \(d_2\) | \((\mu_2,r_2)\) | \((r_0,r_1,r_2)\) |
|---:|---|---|
| 2 | \((2,11),(4,7),(6,3)\) | \((6,4,11),(10,4,7),(8,6,3)\) |
| 3 | \((2,5),(4,2)\) | \((9,6,5),(15,6,2)\) |
| 4 | \((2,3)\) | \((12,8,3)\) |
| 5 | \((2,2)\) | \((15,10,2)\) |
| 6 | none | |

For the \(d_2=2,\mu_2=6\) row the coprime prefix
\((7,2)\) also solves the prefix conductor equation, but
\(2\cdot3\notin\langle14,4\rangle\).
For \(h\ge3\), strict divisibility and \(d_2\le6\) require
\(h=3\) and \((d_2,d_3)=(4,2),(6,2)\), or \((6,3)\).
At \(d_2=6\) the recurrence bounds \(\mu_2\) by two, and the
prefix is necessarily \((18,12)\). Freeness gives \(e_2r_2\ge12\).
The first increment after the prefix is then at least four
for \(d_3=2\) or at least three for \(d_3=3\), already
beyond fourteen. At \((d_2,d_3)=(4,2)\), the prefix is
\((12,8)\) since \(\mu_2\le3\). Freeness requires
\(2r_2\in\langle12,8\rangle\). Together with
\(\gcd(r_2,4)=2\), this gives \(r_2\ge6\). The formula becomes
\(r_2+r_3=9\). Since \(r_3\ge2\) is odd and
\(r_2\equiv2\bmod 4\), it forces \((r_2,r_3)=(6,3)\).

## Exclusion of the listed embeddings

The [global conductor-fourteen exclusion](GLOBAL_CONDUCTOR_FOURTEEN_EXCLUSION.md)
combines the branch proofs and gives the next lower bound.
