# Central Newton faces for aligned unit-gap factorial tails

Fix an ordinary two-monomial core, consecutive orders
\(n,N=n+1\) with \(n\ge2\), and an integer signed tail gap
\(h\). The tail consists of \(R=S+h\) increasing rows
\((a,a+1)\) and \(S\) decreasing rows \((a+1,a)\).
Exchange the two columns of the core and tail if needed to
take \(h\ge0\).

Choose any prime \(p\) dividing \((n-1)nN\) and an
integer \(E\) with \(p^E>N\). Assume every smaller tail exponent satisfies
\(a+1\equiv0\bmod p^E\). The bases and their multiplicities
can otherwise vary independently, without a common height
window. Write \(C_j(z)=\sum_k c_{j,k}z^k\) for the core
moment at order \(j\). Its coefficients include the binomial
factors. Any fixed factorial scalar may be divided out of
each core moment to give integral coefficients.

**Theorem.** Above a computable threshold in \(S\), the
two tail-augmented Newton polygons have disjoint exterior
slopes. Those slopes are also disjoint from both central
slope sets. The common central slope and residue data are
read from two fixed explicit polynomials, depending only
on the core, \(h\), and \(S\bmod(p-1)\). In particular a
finite residue test can uniformly exclude all common
quadratics, or all common factors, for the stated family.

At \(p=2\), \(n=2^q\) and the empty core, the test excludes
every quadratic for every fixed \(h\). If \(|h|\ne1\) it
excludes all nonconstant common factors. If \(|h|=1\),
every nonconstant common factor has degree at least
\(2^{\lceil q/2\rceil}\); for \(q\le3\) the moments are
coprime. These assertions allow \(S\) arbitrarily large and
the individual bases arbitrarily far apart in their
congruence classes. The first nontrivial cases have \(q\ge2\).

## The fixed central polynomials

Put \(v=v_p\), and for each \(j=n,N\) define, at
\(1\le k<j\),
\[
 u_{j,k}=v(c_{j,k})+
        h\bigl(v((j-1)!)-v((j-k-1)!)\bigr),\qquad
 \mathcal I_j=\left\{k:p\nmid\binom{j-2}{k-1}\right\}.
\]
For \(0\le k<j\), let
\[
 \omega_j(k)=\prod_{t=1}^{k}
         \overline{(t-j)/p^{v(t-j)}}\in\mathbb F_p^\times .
\]
The bar denotes reduction of an integer unit. At \(k\in\mathcal I_j\)
define the residue
\[
 \gamma_{j,k}=
    \overline{c_{j,k}/p^{v(c_{j,k})}}\,
                \omega_j(k)^R\omega_j(j-k)^S .
\]
It is determined by \(S\bmod(p-1)\) since \(R=S+h\).
Use its representative in \(\{1,\ldots,p-1\}\) to form
\[
                      L_j(z)=\sum_{k\in\mathcal I_j}
                                      \gamma_{j,k}p^{u_{j,k}}z^k. \tag{1}
\]
In (1) the same symbol denotes that representative.
The endpoint exponents of (1) are \(1,j-1\).

For a common slope \(a/b\) in lowest terms, compare the
two ordinary stepped residual polynomials of \(L_n,L_N\)
on that slope. The maximum width still available for a
common factor there is \(b\) times the degree of their gcd.
This follows from multiplication of initial forms, as in
the [truncated face comparison](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md).
The sum of these widths bounds the degree of any rational
common factor. For the quadratic test retain also the
requirement that a width on denominator \(b\) be a multiple
of \(b\). This excludes a quadratic even when some higher
degrees remain available. A total bound of one suffices for
physical coprimality since the even factorial moment has
strictly positive values at real arguments.

## Separation of the exterior slopes

Strip the factorial of the smaller endpoint from each tail
row in each order. The factor at offset \(t\) is
\[
                         ja+t=j(a+1)+(t-j).
\]
For \(1\le t<j\), the valuation and residue of the unit
part equal those of \(t-j\), because \(p^E>N\).
At \(t=j\) the valuation is \(v(j)+v(a+1)\).
Only endpoint coefficients sample the latter factor.
For the interior, the increasing row's valuation is
\(v((j-1)!)-v((j-k-1)!)\), and the decreasing row's is
\(v((j-1)!)-v((k-1)!)\). It follows that, after removing
the constant valuation \(S[v((j-1)!)+v(j-1)]\) at order \(j\),
the interior ordinates are exactly
\[
                  S\,v\binom{j-2}{k-1}+u_{j,k}.              \tag{2}
\]
Their coefficient units are given by the formula for
\(\gamma_{j,k}\).

Let \(T_+\) and \(T_-\) be the sums of \(v(a+1)\) over the
two respective tail orientations. Put
\(\mu_j=v(j)-v(j-1)\). After the same scalar removal, the
two endpoint ordinates are
\[
          T_-+S\mu_j+v(c_{j,0}),\qquad
          T_++S\mu_j+h v(j!)+v(c_{j,j}).                       \tag{3}
\]
If \(p\mid n\), then \(\mu_n=v(n)\) and \(\mu_N=-v(n)\).
Otherwise \(\mu_n=-v(n-1)\) and \(\mu_N=v(N)\).
In either case \(\mu_n\ne\mu_N\) and \(E+\mu_j\ge1\).
Since \(T_-\ge ES,\ T_+\ge ER\), we can choose \(K\ge1\) to bound all \(u_{j,k}\) and the fixed
endpoint offsets in (3). Every endpoint is at least \(S\)
after scalar removal. For \(S>4K\), the lower hull between
indices \(1\) and \(j-1\) is that of the entries
\(\mathcal I_j\) alone. In fact it has height at most \(K\)
throughout that interval, whereas each omitted ordinate
is at least \(S\). Its slopes have absolute values at most
\(K\). The endpoint-to-neighbor slopes have opposite signs
and absolute values at least \(S-K\), so adding these two
width-one edges already gives the full lower hull.
The corresponding left and right exterior slopes at the
two orders differ in magnitude by \(S|\mu_n-\mu_N|\) up to
an error of at most \(4K\). When \(p\mid n\) the main
difference is \(2v_p(n)S\). This proves disjointness
and the assertion about (1).

## The dyadic power-of-two pair without a core

Here \(p=2\), \(n=m=2^q\), \(N=m+1\) and
\(c_{j,k}=\binom jk\). For the lower central polynomial
\(\mathcal I_m\) consists of the odd indices: modulo \(2\),
\((1+x)^{m-2}\) is the product of \(1+x^{2^i}\)
over the nonzero binary digits of \(m-2\).
For the upper central polynomial \(\mathcal I_N\) is the
full interval \(1\le k\le m\), since \(m-1\) has all its
lower \(q\) binary digits equal to one.

When \(h=0\), the lower central ordinates are constant.
The upper minimum occurs exactly at \(k=1,m\), since
\[
 v_2\binom{m+1}{k}
       =q-v_2(k)-v_2(k-1)\quad(2\le k\le m).
\]
Thus the horizontal residual polynomials, after stripping
the initial monomial, are \((1+z)^{m-2}\) and
\(1+z^{m-1}\). Their gcd has degree one because \(m-1\)
is odd. The other upper slopes cannot match a horizontal
lower slope. This proves coprimality.

For \(h>0\), subtract the common linear slope \(hk\).
The lower odd ordinates become \(q-hs_2(k)\), where
\(s_2(k)\) is the binary digit sum. Its central hull has
vertices \(2^r-1\), \(1\le r\le q\), with slopes
\(-h/2^r\) and widths \(2^r\), \(1\le r<q\).
Between successive vertices the digit sum is at most
its earlier vertex's value, so every intermediate point
lies strictly above the edge.

The upper shifted ordinate at \(1\) is \(h(q-1)\).
At \(2^r\), \(1\le r\le q\), it is
\(h(q-1-r)+q-r\). At every index strictly between
\(2^r\) and \(2^{r+1}\) it is at least the earlier
vertex's ordinate. This follows from
\(s_2(m-k)=q-s_2(k-1)\), \(s_2(k-1)\le r\), and the
displayed binomial valuation bound \(q-r\).
Consequently all upper edges after its first one have
slopes \(-(h+1)/2^r\) for some \(r\).
The first upper edge has odd width \(2^s-1\), ending
at one of the powers of two.

For \(h\ge2\), an equality between \(-h/2^r\) and a
subsequent upper slope would require \(h+1=2^\ell h\)
for an integer \(\ell\), which is impossible. A common
slope can only be the first upper slope. If it matches,
it has integer value since its reduced denominator
divides an odd width and a power of two. The first
upper face then contains only its endpoints. Otherwise
the same slope would occur among successive powers of two.
The stepped residue polynomials are binomials of odd and
power-of-two degrees over \(\mathbb F_2\), with gcd of
degree one. No other slope remains and parity excludes
a common factor.

For \(h=1\) the first upper slope cannot match a lower
one, again by its odd width and the nonintegral lower
slopes. A later match with lower width \(2^r\) would
have the upper edge starting at \(2^{r+1}\), of slope
\(-1/2^r\). The chord from \(1\) to \(2^{r+1}\) must
then have slope at most \(-1/2^r\). Its slope is
\[
                      \frac{q-2(r+1)}{2^{r+1}-1}.
\]
The condition is equivalent to \(q\le2r\), since
\(q-2r\) is an integer. Thus any matched slope has
reduced denominator \(2^r\ge2^{\lceil q/2\rceil}\).
The stated degree restriction follows from the Newton
polygon product rule. A later upper edge must also end
by \(2^q\), requiring \(r+2\le q\), which rules out all
matches for \(q\le3\). Reversal handles \(h<0\).

## Consequences for diffuse tails and all-degree certificates

For example, use core row \((7,17)\), \(n=16,N=17,p=2\)
and \(h=0\). With the smaller-end factorial stripped
from the core at each order, \(K=173\) suffices. The two
central slope sets are
\[
      \{99/10,21/2,11\},\qquad \{28/3,79/8,10\}.
\]
They are disjoint. Thus any balanced unit-gap tail with
each smaller exponent congruent to \(-1\bmod32\) and
at least \(693\) rows of each orientation has coprime
moments at this pair.

In any positive diffuse limit of unit-gap tails, the
maximum relative gap \(\eta\) tends to zero and each
row contributes at most \(\eta\) to curvature.
Positive limiting curvature and fixed \(h\) force
\(S\to\infty\). Thus every passing fixed central test
excludes accumulating exact ordinary collisions at that
pair over the whole indicated residue class, including
at positive Gaussian exceptional parameters.

For each \(2\le n<N<J\), choose a prime
\(p\mid(n-1)nN\) and exponent \(E\) as above. Requiring all \(a+1\)
to be multiples of the least common multiple of
the finitely many \(p^E\)'s imposes all their conditions
at once. At any specified \(h\), check the finite
central tests, including the allowed residues of \(S\)
modulo the corresponding \(p-1\)'s. A pass at every
pair supplies an effective all-degree threshold in \(S\).
If \(h\) has not been fixed, the
[absolute signed gap bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
restricts it to an explicit finite range at any physical
collision after subtracting the core's signed gap.
The finite central tests can then be run for that range.
The first pair follows by even positivity, and the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
handles \(N\ge J\). The thresholds impose no upper bound
on any ordinary tail height, and use the actual arithmetic
faces rather than replacing the tails by their Gaussian
limits.
