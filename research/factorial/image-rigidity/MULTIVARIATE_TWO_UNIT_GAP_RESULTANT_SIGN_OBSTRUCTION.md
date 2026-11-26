# A sign obstruction for the symmetric two-unit-gap resultant

Consider two increasing rows \((a,a+1),(b,b+1)\), allowing
the heights \(a,b\) to be indeterminates. Equal exponent
rows affect only the nonzero moment scalars. With
\(S=a+b,\ P=ab\), the two-row normalized moment polynomial
has monic reversal
\[
 Y_n(w)=\sum_{k=0}^n\binom nk w^{n-k}
              \prod_{j=1}^k(j^2+njS+n^2P).
\]
Put \(R_n(S,P)=\operatorname{Res}_w(Y_n,Y_{n+1})\) with
the usual Sylvester determinant orientation, and
\(H_n=R_n/(1+S+P)\). The factor has multiplicity one
at every start by the
[monic remainder and factor theorem](MULTIVARIATE_UNIT_GAP_MONIC_REMAINDER.md#exact-multiplicity-of-the-forced-symmetric-factor).

**Exact sign result through \(n=7\).** The coefficients
of \(R_n,H_n\) in \(S,P\) are all positive for
\(1\le n\le6\). At \(n=7\), all coefficients remain
positive except at the following four positions in each
polynomial:
\[
 \begin{array}{c|c}
   \text{polynomial}&\text{monomials with negative coefficients}\\ \hline
   R_7 & P^{36},P^{37},P^{38},P^{39}\\
   H_7 & P^{35},P^{36},P^{37},P^{38}.
 \end{array}
\]
For example, the coefficient of \(P^{39}\) in \(R_7\) is
\[
 -9508445284953248328885758152777703463189896134913156560771221402363311847372186386362265625000.
\]
Thus positivity in the elementary symmetric generators
fails even after division by the forced factor.

Both \(R_n(a+b,ab)\) and \(H_n(a+b,ab)\) have strictly
positive coefficients at every \(1\le n\le7\).
The shift \(a=1+x,\ b=1+y\) also gives positive coefficients
in the elementary symmetric generators \(x+y,xy\) at each
of these orders. The shift replaces \((S,P)\) by
\((S'+2,P'+S'+1)\). Consequently the negative coefficients
above leave the positive numerical resultant for
\(a,b\ge0\) intact at \(n=7\). The sign result places no
restriction on either height at these starts. Positivity at
arbitrary starting order in the original heights or the
shifted elementary generators remains open.

## Exact coefficient reconstruction

Put \(D=n(n+1)\). The coefficient of \(w^{j-k}\) in \(Y_j\)
has total degree at most \(k\) in \(S,P\). In any nonzero
term of the Sylvester determinant, the sum of the
coefficient indices is \(D\). Indeed the sum of the selected
column indices is \((2n+1)2n/2\), and the row shifts
sum to \(n^2\). This gives the degree bound
\(\deg_{S,P}R_n\le D\). The highest pure \(P\) coefficient
is \((2n+1)^D\), from retaining the highest \(P\) term in
every rising factor. The resulting monic polynomials
are \((w+n^2P)^n\) and \((w+(n+1)^2P)^{n+1}\).

The [standard-library checker](verify_multivariate_two_unit_gap_resultant_sign.py)
evaluates the full integer Sylvester determinants by
fraction-free Bareiss elimination at every integer
\((S,P)\) with \(S,P\ge0\), \(S+P\le D\).
This triangular grid determines any polynomial of total
degree at most \(D\). Iterated differences at the origin
give its coefficients in the basis
\(\binom Si\binom Pj\), \(i+j\le D\). The checker
converts to monomials with signed falling-factorial
coefficients and a common integer denominator, and checks
exact divisibility. It divides by \(1+S+P\) using that
the divisor is monic in \(S\). It then expands the ordinary
and shifted substitutions by the binomial theorem.
The inspected coefficient counts are
\[
 \begin{array}{c|rrrrrrr}
 n&1&2&3&4&5&6&7\\ \hline
 R_n&6&28&91&231&496&946&1653\\
 H_n&3&21&78&210&465&903&1596 .
 \end{array}
\]
The calculation checks every coefficient of each indicated
polynomial. It is an exact finite sign determination.

    python3 research/factorial/image-rigidity/verify_multivariate_two_unit_gap_resultant_sign.py
