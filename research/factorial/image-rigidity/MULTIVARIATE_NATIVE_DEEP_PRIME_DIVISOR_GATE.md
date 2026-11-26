# Explicit divisor conditions for a native deep prime in ordinary binomials

Consider an ordinary two-monomial profile, with arbitrary shared
and private rows. Select one unequal row and write its smaller
exponent and absolute gap as \(A,d\), \(d>0\). The rest of the
profile may be fixed or vary independently of this row.

Fix adjacent orders \(n,N=n+1\), \(n\ge2\). Choose a marked order
\(m\in\{n,N\}\), write \(m'\) for the other order, and choose
\(1\le h_0\le d\) with \(\gcd(h_0,m)=1\). For every retained row
\(i\) let \(\ell_i,d_i\) be its smaller exponent and absolute gap.
Define the positive integer
\[
 B_j=j!\prod_{i\ {\rm retained}}\prod_{h=1}^{j d_i}(j\ell_i+h)
                         \qquad (j=m,m')
\]
and put
\[
 \begin{aligned}
 C_{\rm mark}&=B_m\prod_{\substack{1\le h\le md\\h\ne h_0}}
                                                   |h-h_0|,\\
 C_{\rm other}&=B_{m'}\prod_{h=1}^{m'd}|mh-m'h_0|,
              \qquad C=C_{\rm mark}C_{\rm other}.
 \end{aligned}                                               \tag{1}
\]
These constants do not involve the marked row's base \(A\).

**Theorem.** A common root at the physical argument forces
\[
                              mA+h_0\mid C^m .               \tag{2}
\]
Thus failure of any one such divisor test certifies consecutive
coprimality. The marked row, marked order, and offset may be
chosen independently at each pair. In particular a possible
common root must satisfy the divisor conditions with \(h_0=1\)
at both \(m=n,N\), whose two left hand sides are coprime.

The proof gives a native prime at which the two Newton polygons
have disjoint slopes whenever (2) fails. Neither the prime nor
the resulting marked slope needs to meet the large-prime or
primitive-denominator hypotheses of earlier angular tests.
All assertions here are for physical coprimality.

## Valuations at a prime dividing a proposed marker

Exchange the two columns if necessary so the selected row
increases. This inverts all nonzero roots and preserves
physical common-root questions. At each degree \(j\) strip
the positive scalar \((j\ell_i)!\) from each row, including
\((jA)!\) from the selected row. The coefficients at \(k\)
take the form
\[
                         b_{j,k}(jA+1)_{dk},
\]
where \(b_{j,k}\mid B_j\). Indeed the retained contribution
of row \(i\) uses a prefix of the displayed full product,
in the appropriate orientation, and the binomial
coefficient divides \(j!\).

Let \(p\) be a prime divisor of \(mA+h_0\), and write
\[
        E=v_p(mA+h_0),\qquad
        U=v_p(C_{\rm mark}),\qquad L=v_p(C_{\rm other}).
\]
The coprimality of \(m,h_0\) gives \(p\nmid m\).
Suppose \(E>m(U+L)\). In the marked degree, for any
\(h\ne h_0\) in the full ramp,
\[
                v_p(mA+h)=v_p(h-h_0),
\]
since the latter value is at most \(U<E\).
The marker itself is in the first block and is sampled
at every nonzero index. Thus the marked coefficient
valuations have the form \(u_0\) at \(k=0\) and
\(E+u_k\) at \(1\le k\le m\), with \(0\le u_k\le U\).

In the other degree the identity
\[
 m'A+h=\frac{m'}m(mA+h_0)+\frac{mh-m'h_0}{m}
\]
gives \(v_p(m'A+h)=v_p(mh-m'h_0)\).
The numerator of the latter term is nonzero:
if \(mh=m'h_0\), then \(m\mid h_0\) by
\(\gcd(m,m')=1\). Its valuation is at most \(L<E\),
and \(p\nmid m\). Consequently every other-degree
coefficient valuation lies in \([0,L]\).

At interior index \(k\) of the marked degree, the scaled
height above the endpoint chord is at least
\[
                      (m-k)E-mU>0 .
\]
Its entire polygon is therefore the endpoint chord,
of slope \((E+u_m-u_0)/m>L\).
Every slope in the other polygon is at most \(L\),
since its ordinates lie in \([0,L]\) and every edge has
positive integer width. The polygons have disjoint
slopes, so the polynomials have disjoint root valuations
at this prime. Thus a common root requires
\(E\le m(U+L)=m v_p(C)\) at every prime dividing
\(mA+h_0\). This proves (2).
Rational coprimality extends to each characteristic-zero
field by a rational Bézout identity.

## Exact local bounds and simultaneous offsets

A sharper divisor, if wanted, follows from the exact
stable ordinates. With \(p\nmid m\), let \(\tau_p\) be
the maximum valuation of any numerator in the offset
products of (1). For \(E>\tau_p\) write the marked
ordinates as \(u_0,E+u_1,\ldots,E+u_m\), as above,
and let \(l_k\) be the fixed ordinate in the other
order. These residual ordinates are computed just
from the retained coefficient factors and the fixed
offsets. The largest other slope is
\[
                     s_{m'}=\max_{0\le k<m'}
                                  \frac{l_{m'}-l_k}{m'-k}.
\]
Indeed the last edge gives equality and its supporting
line lies at or below every other ordinate.
The marked chord is strictly above this slope, with
its interior ordinates strictly above the chord, if
\[
 E>\max\left\{
       0,\tau_p,\ m s_{m'}+u_0-u_m,\
       \max_{1\le k<m}
         \frac{(m-k)u_0+k u_m-m u_k}{m-k}
                     \right\}.                              \tag{4}
\]
Thus one may replace \(C^m\) in (2) by the product of
\(p\) raised to the floor of (4)'s right side, over
\(p\mid C\) with \(p\nmid m\).
For primes \(p\nmid C\), the threshold is zero.

For any collection of offsets on the same row, indexed
by \(i\), write \(Y_i=m_iA+h_i\) and let \(M_i\) be
the associated proposed divisor from (2) or (4).
If the reduced fractions \(h_i/m_i\) are distinct,
then \(\gcd(Y_i,Y_j)\mid |m_i h_j-m_j h_i|\), a
fixed nonzero integer. The simultaneous conditions
give
\[
                 \prod_i Y_i
    \le \operatorname{lcm}_i M_i
                 \prod_{i<j}|m_i h_j-m_j h_i|.              \tag{5}
\]
At each prime, the excess of the sum of the valuations
of the \(Y_i\) above their maximum is at most the sum
of their pairwise minimum valuations. This proves
(5). In particular at \(h_i=1\) in the two orders,
\((nA+1)(NA+1)\mid\operatorname{lcm}(M_n,M_N)\).
The simultaneous test can therefore use a square-root
upper bound on the exceptional bases.

## An all-degree native certification region

Take \(J\ge4\) at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
At each \(2\le n<N<J\), test any available row with
any allowed marked order and offset. If some condition
(2) fails at every such pair, SFC(2) follows at every
order. The first and second moments are coprime by
strict positivity of the quadratic over the reals,
and the absolute theorem handles \(N\ge J\).

In particular, for one designated row with all others
held fixed, set \(m=N,h_0=1\) for every terminal below
\(J\). An explicit finite threshold is
\[
                     A>\max_{3\le N<J} \frac{C_N^N-1}{N},
                                                               \tag{3}
\]
where \(C_N\) is (1) at that terminal.
Then every such pair has a native certifying prime.
This is stronger than certifying one progression:
every integer beyond the threshold is covered without
a complex root or resultant test. The numerical value
of \(J\) can remain unevaluated.

The same argument allows all bases to vary independently.
With the signed gaps fixed and
\(G=\sum_{i\ {\rm retained}}d_i\), the bound \(C^m\)
is a polynomial in the retained bases of total degree
\(m(m+m')G\). If \(B\) bounds all retained smaller
exponents, it is bounded by an explicit constant times
\((1+B)^{m(m+m')G}\). For instance the constant
can be taken from (1) evaluated with all those bases
equal to one. It follows that a single base large
enough relative to this power of all other bases
gives native all-degree coprimality, with no hypothesis
on the retained core. Alternatively the certifying row
can change from one terminal to the next, since (2)
is necessary for each possible choice at that terminal.

There is also a uniform count for failures of this
arithmetic certificate in boxes of independently
varying bases of height \(H\) and fixed gaps.
Condition on all bases except one designated base
and use \(m=N,h_0=1\). The integer in the right side
of (2) is polynomially bounded in \(H\).
The number of positive \(A\le H\) satisfying (2)
is at most the divisor count of that integer.
For each \(\varepsilon>0\) this is \(O(H^\varepsilon)\)
by the elementary estimate
\(\tau(X)=O_\eta(X^\eta)\) for every \(\eta>0\).
For completeness, at each prime above \(2^{1/\eta}\)
use \(e+1\le2^e\le p^{\eta e}\) for its exponent \(e\);
the finite smaller primes contribute a bounded
constant from \(\sup_e (e+1)/p^{\eta e}\).
The union over \(N<J\) gives failure proportion
\(O(H^{-1+\varepsilon})\) for this native physical
certificate. The same count applies to fixed arithmetic
progressions in the independent bases.
This is a count of certificate failures, not of
actual common roots.

## A small certifying prime with an integral slope

At \(n=8,N=9\), take the selected row
\((A,A+1)\), \(A=(2^{36}-1)/9\), and the retained
private row \((1,0)\). The base \(A\) is integral
because \(2^6\equiv1\bmod9\). After stripping the
row scalars, \(b_{m,k}=m!/k!\).
At \(p=2\), the marker \(9A+1\) has valuation \(36\).
In order nine the initial ordinate is \(7\), and the
other ordinates at \(1\le k\le9\) are
\(36+7-v_2(k)\). They give one edge of slope \(4\),
with all interior ordinates strictly above it.
In order eight, the values
\(v_2(9h-8)\) for \(1\le h\le8\) are
\[
                               0,1,0,2,0,1,0,6 .
\]
Adding the retained \(8!/k!\) valuations gives the
array \((7,7,7,7,7,7,7,7,10)\). Its slopes are zero
and three. Equivalently the threshold in (4) is \(27\),
strictly smaller than \(36\).
Thus \(2\) certifies physical coprimality
by disjoint slopes, despite the nonprimitive upper
chord and the nonunit binomial coefficients.
