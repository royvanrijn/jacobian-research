# Leading faces of winding-one cubic sphere phases of degree seven

Use the sphere functional and regular presentation
\(f=x^{-1}(C+Bt+Dt^2+Et^3)\) from the
[nonunit-center classification](MIXED_CUBIC_NONUNIT_CENTER_EXCLUSION.md).
Here the field has characteristic zero, \(E(0)\ne0\), and the least sphere
degree of \(f\) is seven. A phase is eventually pure if
\(\Lambda(f^m)=0\) for all sufficiently large \(m\).

**Proposition.** If \(f\) is eventually pure, then exactly one of the
following necessary alternatives holds.

- \(\deg E=4\), with one geometric root of multiplicity four.
- \(\deg E=5\), with at least two geometric roots.

In either alternative the top Newton face at angular exponent zero has
height three and is horizontal. In the degree-five alternative, if
\(\Lambda(x f^l)\ne0\) for some fixed positive \(l\), the geometric root
multiplicities of \(E\) must be \((1,1,1,1,1)\) or \((4,1)\).

## The possible upper faces

The only negative angular exponent in the literal Laurent expansion is
\(-1\), with highest occupied node \((-1,3)\). At nonnegative exponents
the height is at most three. The
[sphere degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md)
gives an occupied nonnegative node of degree \(i+j=7\), and none of larger
degree. If there is no positive angular coefficient, the powers have
constant term equal to the powers of the angular core. The interval
prime argument in the nonunit-center proof either rules out eventual
purity or forces this core to vanish. With zero core and only one negative
angular exponent, the phase is supported only on \(y\) and \(yt\), of
degree at most two, which is impossible here.

Write an upper supporting line at angle zero as \(j+s i=H\). A face with
only one strict angular sign or with exactly two opposite-sign nodes is
excluded by the [prime test for faces](MIXED_CUBIC_NONUNIT_CENTER_EXCLUSION.md#the-prime-test).
For two nodes \((-1,3),(a,j)\), take \(p\equiv1\pmod{2N}\),
\(N=3a+j\), and \(m=(a+1)(p-1)/N\). The same binomial proof applies since
\(N\ge a+1\). Thus a possible pure face has both signs and at least
three occupied nodes. If \(s<0\), the negative node has weight \(3-s>3\)
and all nonnegative nodes have weight at most three. If \(s\ge1\), the
degree-seven node has weight at least seven, exceeding \(3-s\). A
nonhorizontal candidate therefore has \(0<s<1\) and
\[
 7s\le H=3-s,qquad s\le3/8.
\]
In reduced terms write \(s=p_0/q\). On the line, the occupied integer
angles differ from \(-1\) by multiples of \(q\), and each such step
decreases the height by \(p_0\). Two additional occupied nodes force
\(2p_0\le3\), so \(p_0=1\), and then \(q\ge3\). The second possible
step has angle \(2q-1\) and height one. Its degree \(2q\) is at most
seven only when \(q=3\), and any subsequent step has degree at least
eight. Thus the sole nonhorizontal candidate has \(s=1/3,H=8/3\),
with all three angles \(-1,2,5\).

## The nonhorizontal order-eight cover

For this face write its angular Laurent polynomial as \(z^{-1}P(z^3)\)
with a quadratic \(P\) of nonzero constant term. For a prime
\(p=8k+1\) and order \(m=3k\), the highest angular-constant height is
\(p-1\) and its coefficient is
\[
 [z^k]P(z)^{3k}=[z^{p-1}]\bigl(z^7P(z)^3\bigr)^k.
\]
The normalized cyclic cover \(w^8=z^7P(z)^3\) is connected, because the
branch exponent seven is coprime to eight. In the raw degree-sixteen
completion its branch exponents are \((7,3,3,3)\) or \((7,3,6)\),
according as the two quadratic roots are distinct or coincide. All are
less than eight. Thus the character-one summand is \(O(-2)\) and the
displayed coefficient is its
[coherent Frobenius scalar](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#4-the-exact-normalized-cover-frobenius-line).
At an ordinary prime of the primitive factor it is nonzero.

On the primitive factor the multiplicities, in character order
\(1,3,5,7\), are \((1,0,2,1)\) in the distinct-root case and
\((1,0,1,0)\) in the repeated-root case. Indeed for branch exponents
\(a_j\) the entry at \(i\) is \(-1+\sum_j\{i a_j/8\}\).
The deck action preserves the Jacobian polarization, so its Rosati
involution acts as complex conjugation on \(\mathbb Q(\zeta_8)\);
the restricted polarization gives the same action on the primitive
factor. The
[rank-two ordinary-prime proposition](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md#ordinary-primes-for-rank-two-multiplication)
applies in the distinct-root case with one balanced pair. The repeated-root
case has rank one and hence potential CM, to which
[Fité, Remark 12](https://arxiv.org/pdf/2005.10185v3) applies. The passage
from the primitive factor to the selected character uses its rational
cyclotomic idempotent in the deck action, which commutes with Frobenius
at the good primes under consideration.

Adjoin \(\zeta_8\) to a number field of definition of the cover, of its
deck action and, if needed, of its roots. The ordinary inputs give
infinitely many good absolute degree-one primes over rational
\(p\equiv1\bmod 8\), so \(p\Lambda(f^{3k})\) reduces to a nonzero scalar by
the same denominator calculation as in the two-node test. Thus every
nonhorizontal face is excluded.

## The horizontal test at order three

For the horizontal face the angular Laurent polynomial is \(E(x)/x\),
and both angular signs occur, so \(2\le\deg E\le5\). Take
\(p=3k+1\) and moment order \(k\). The top possible height is \(3k=p-1\),
and the selected coefficient is
\[
 A_p=[x^k]E(x)^k.
\]
If \(\deg E\le4\) and all root multiplicities are at most two, the
connected normalized cover \(w^3=x^2E(x)\) has exponents two at zero,
\(4-\deg E\) at infinity, and the indicated finite-root multiplicities.
Its raw completion has degree six, every exponent is below three, and
the character-one summand is \(O(-2)\). Its Frobenius scalar
\([x^{p-1}](x^2E)^k\) equals \(A_p\).
The genus is at most three by the cyclic-cover ramification formula.
The ordinary-prime input in genus one or two follows from
[Sawin, Theorem 3](https://arxiv.org/html/1506.04784v2#S0.Thmstatement3),
using the square of the elliptic curve in genus one. In genus three the
deck multiplication by \(\mathbb Q(\zeta_3)\) meets the hypothesis of
[Fité, Theorem 1](https://arxiv.org/pdf/2005.10185v3). Applying the
theorem over a field containing \(\zeta_3\) supplies a positive-density
set of ordinary primes, meeting the density-one set of absolute
degree-one primes. These primes make \(A_p\ne0\).

When \(\deg E\le4\) and a root \(r\) occurs exactly three times, write
\(E=(x-r)^3R(x)\) with \(\deg R\le1\) and \(R(r)\ne0\).
Over good primes of absolute degree one in a field containing \(r\),
\[
 (x-r)^{p-1}=\sum_{i=0}^{p-1}r^{p-1-i}x^i\pmod p .
\]
Since \(\deg(R^k)\le k\), this gives
\(A_p=r^{p-1-k}R(r)^k\ne0\).
Finally, if \(E\) is a scalar times \((x-r)^5\), the coefficient is a
unit times \(\binom{5k}{k}\). As \(5k=p+2k-1\), its reduction is
\(\binom{2k-1}{k}\ne0\). This gives the two remaining degree and
root-count alternatives in the proposition.

## The second entry when a fixed marker survives

Suppose \(\deg E=5\). For fixed \(l\ge0\) and sufficiently large good
absolute degree-one primes \(p=3k+1\), reduction modulo a prime above
\(p\) gives
\[
 f^{lp+k}\equiv f(x^p,t^p)^l f(x,t)^k .
\]
In multiplying \(\Lambda(f^{lp+k})\) by \(p\), the contributing
heights are \(p-1+pj\), since only their denominators are divisible by
\(p\). For large enough \(p\) relative to \(l\), their scaled denominators
are the units \(j+1\). The \(t\)-height of \(f^k\) is at most
\(3k=p-1\), so only its top \(t\)-part contributes. This part has angular
exponents in \([-k,4k]\), whose only multiples of \(p\) are zero and \(p\).
As \(p-1\) is even and \(p\) is odd, the remaining parity condition is
exactly that used in the interval functional at height \(j\). Consequently
\[
 p\Lambda(f^{lp+k})\equiv
 A_p\Lambda(f^l)+B_p\Lambda(xf^l),
 \qquad B_p=[x^{p+k}]E(x)^k.
\]
Eventual purity makes \(A_p=0\) at all sufficiently large good primes by
\(l=0\), and a nonzero \(\Lambda(xf^l)\) forces \(B_p=0\) there as well.
If \(C(0)\ne0\), the marker hypothesis already holds for \(l=1\), since
\(\Lambda(xf)=C(0)+D(0)/3=2C(0)/3\).

If every multiplicity of \(E\) is at most two, the above order-three cover
has raw completion of degree nine with exponent two at both zero and
infinity. Its character-one summand is \(O(-3)\), with basis
\(w/x,w/x^2\) in cohomology. The first row of Frobenius in that basis
is \((A_p,B_p)\), by taking the coefficient of \(w/x\) after
raising each basis element to the \(p\)th power. The multiplicity
partitions \((2,1,1,1)\) and \((2,2,1)\) give dimensions four and three
with signatures \((2,2)\) and \((2,1)\), respectively, by the same
fractional-part formula. Fité's Theorem 4 in the former case and
Theorem 1 in the latter give infinitely many ordinary absolute
degree-one primes over a field containing \(\zeta_3\). Ordinarity makes
Frobenius invertible on this summand, precluding a zero row.

If a root occurs exactly three times, use \(E=(x-r)^3R(x)\) again.
Now \(\deg R=2\), and, writing \(R^k=\sum_j c_jx^j\), the contributions
to \(A_p\) come from \(0\le j\le k\) and the contributions to
\(B_p\) come from \(k+1\le j\le 2k\). The same expansion of
\((x-r)^{p-1}\) gives
\[
 A_p+r^p B_p=r^{p-1-k}\sum_{j=0}^{2k}c_jr^j
              =r^{p-1-k}R(r)^k\ne0.
\]
This eliminates both partitions containing a triple root under the
nonzero marker premise.

In all sections, one number-field specialization from a finitely
generated coefficient algebra preserves the occupied support, the
nonzero differences of distinct geometric roots and, when used, the
nonzero value of \(\Lambda(x f^l)\). Every putative tail identity survives
specialization because each moment is a rational polynomial in the
coefficients. The prime contradictions therefore apply over every
characteristic-zero field. The [logarithmic boundary reduction](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md)
narrows these necessary alternatives, and the
[degree-five leading-coefficient exclusion](WINDING_ONE_CUBIC_FIVE_DEGREE_CM_EXCLUSION.md)
leaves only the degree-four leading coefficient with one quadruple root.
