# A terminal-six \(5\)-adic lift for mixed unit gaps

Let the two exponent columns be distinct nonnegative integer vectors.
Remove equal rows. Suppose every remaining row has gap \(1\) or \(-1\).
Write \(r\) and \(s\) for the numbers of increasing and decreasing rows,
and \(a_i\) for the lower exponent in row \(i\).

**Theorem.** The ordinary moment polynomials of orders five and six
are coprime in each of the following cases.

- Every \(a_i\) is divisible by \(25\), \(r,s\ge2\), and
  \((r,s)\not\equiv(18,18)\pmod{20}\).
- Every \(a_i\) is divisible by \(11\), and at least one of \(r,s\)
  is nonzero modulo \(10\).

In particular, if every lower exponent is divisible by \(275\), the
polynomials are coprime for arbitrary orientation counts. The support
can be fully shared, with distinct lower exponents, and neither the
dimension nor the exponent size is bounded.

## The \(5\)-adic Newton polygons and their unit residues

Factor out the order-dependent nonzero scalar
\(\prod_i(ja_i)!\) from the ordinary order-\(j\) polynomial.
Its coefficient of \(t^k\) is then
\[
 C_{j,k}=\binom jk
 \prod_{\text{increasing }i}(ja_i+1)_k
 \prod_{\text{decreasing }i}(ja_i+1)_{j-k},
 \tag{1}
\]
where \((x)_\ell=x(x+1)\cdots(x+\ell-1)\).
For \(a_i\equiv0\pmod{25}\) and \(r,s\ge2\), the valuation arrays
of \((C_{j,k})_k\) are
\[
 (s,1,1,1,1,r),\qquad (s,s,1,1,1,r,r).
 \tag{2}
\]
Indeed a row factorial ratio has a single factor divisible
by five exactly when its length is at least five, and that factor
has valuation one. The binomial valuations are one on indices
\(1,\ldots,4\) for order five and \(2,\ldots,4\) for order six,
and zero elsewhere. The slopes of the lower convex polygons are
respectively
\[
 -(s-1),\ 0,\ r-1,\qquad
 -(s-1)/2,\ 0,\ (r-1)/2 .
\]
The only shared slope is zero. A common nonzero root over an
algebraic closure of \(\mathbb Q_5\) must therefore have valuation
zero. Otherwise the minimum over coefficient valuations plus
\(k\) times the root valuation would have to be attained at least
twice in each polynomial, forcing another shared slope.

Let \(F_j(t)=\sum_k C_{j,k}t^k/5\).
In the slope-zero reductions, remove the initial monomial
factors \(t\) and \(t^2\) respectively. The resulting cubic and
quadratic have coefficients
\[
 \begin{aligned}
 R_5(t)&=\sum_{k=1}^4
  \frac{\binom5k}{5}(k!)^r((5-k)!)^s t^{k-1},\\
 R_6(t)&=\sum_{k=2}^4
  \frac{\binom6k}{5}(k!)^r((6-k)!)^s t^{k-2}
 \end{aligned}
 \qquad\bmod 5 .
 \tag{3}
\]
Here all factorials pictured in (3) are units modulo five.
The exact Euclidean calculation in the
[checker](verify_multivariate_terminal_six_padic_unit_gaps.py)
uses the sixteen pairs of count classes modulo four. In fifteen
cases \(R_5,R_6\) are coprime. For \(r\equiv s\equiv2\pmod4\)
they are \((1+t)^3\) and \(2(1+t)^2\). In this last case any
common root of the original pair reduces to \(-1\).

## Lifting the exceptional residue

For every coefficient we have the order-two congruence
\[
 F_j(t)\equiv
 \frac15\sum_{k=0}^j\binom jk(k!)^r((j-k)!)^s t^k\pmod{25}.
 \tag{4}
\]
In the interior ranges of (3), the factorial ratios in (1)
agree with their ordinary factorials modulo \(25\), and the
binomial coefficient supplies the factor five. On all other
indices, the product of factorial ratios has valuation at least
two, since at least two rows have length at least five. The
ratio of this product to the corresponding product of ordinary
factorials is \(1\pmod5\), by writing each individual factor
as \(ja_i+\ell=\ell(1+ja_i/\ell)\) with \(1\le\ell\le6\).
After division by five this still gives agreement modulo \(25\).

For the classes \(r\equiv s\equiv2\pmod4\), the checker evaluates
the constant, linear and quadratic coefficients of \(F_6(-1+w)\)
from (4). They are respectively \(15\bmod25\), \(0\bmod5\)
and \(2\bmod5\). If \(w=t+1\) corresponds to a common root, the
valuation of \(w\) is positive. The sextic equation forces
\(v_5(w)=1/2\). When \(v_5(w)<1/2\) the quadratic term has
the unique smallest valuation, and when \(v_5(w)>1/2\) the
constant has the unique smallest valuation. Put \(w=\pi z\)
with \(\pi^2=5\). Division of \(F_6(-1+\pi z)\) by five
and reduction in the residue field gives
\[
 3+2z^2=0,\qquad z^2=1 .
\]
The reduction of \(F_5(-1+w)\) has its first nonzero term
\(4w^3\). If the constant of \(F_5(-1+w)\) has valuation
one, it has the unique smallest valuation after setting
\(w=\pi z\). Otherwise the constant has valuation at least
two. Divide by \(\pi^3\) and reduce to obtain
\[
 \frac{F_5'(-1)}5 z +4z^3=0 .
\]
As \(z^2=1\), a necessary condition for a common root is
\[
 F_5(-1)\equiv0\bmod25,\qquad F_5'(-1)\equiv5\bmod25 .
 \tag{5}
\]

The exact jet calculation uses (4), multiplying coefficients
modulo \(125\) before dividing by five. Its dependence on each
count is determined by the count modulo \(20\), and by whether
the count equals two. For factorial lengths less than five,
the relevant powers are of units modulo \(25\), since the
binomial contributes five in these terms. For the other terms,
each factorial power of exponent at least three vanishes
modulo \(125\). A complete representative list is
\(2,6,10,14,18,22\) for each count. Substitution into (4)
and the binomial expansion of \(t^k=(-1+w)^k\) give (5) only
for \((18,18)\) among these thirty-six pairs. This proves the
first criterion.

## The residue test at eleven

If each lower exponent is \(0\bmod11\), the normalized
polynomial (1) reduces to
\(\sum_k\binom jk(k!)^r((j-k)!)^s t^k\) modulo \(11\),
and all its coefficients, including both endpoints, are units.
The one hundred pairs of count classes modulo ten suffice.
Their reduced polynomials have a nonconstant gcd only for
\((0,0),(2,5),(5,2)\), of respective degrees \(5,1,1\).
This is the third finite Euclidean calculation in the checker.

Any rational common factor preserves its degree after
primitive reduction at eleven, since both polynomial leading
coefficients are units there. The degree of such a factor is
even, because the ordinary order-six moment is strictly
positive at every real coefficient ratio by the Gamma integral.
Thus the degree-one reductions also exclude a common rational
factor. This proves the second criterion.

When every lower exponent is divisible by \(275\), the
prime-eleven criterion leaves only \(r\equiv s\equiv0\pmod{10}\).
If both counts are positive they are at least two, and
the exceptional class \((18,18)\bmod20\) is impossible, so
the prime-five proof applies. If one count is zero, the
[same-orientation post-prime theorem](MULTIVARIATE_UNIT_GAP_PRIME_TERMINALS.md)
with \(p=5\) applies. Rational Bézout identities transfer
the claims to all characteristic-zero fields.

## Preservation under rowwise dilations

The first two criteria also apply after independently dilating the
exponents of each row by a power of the corresponding odd prime, and
then adding sufficiently divisible common offsets to its two exponents.
For the first criterion this uses preservation modulo \(25\), beyond
the first residue transfer for arbitrary digit blocks.

For an odd prime \(p\) and any integer \(q\ge0\),
\[
 \prod_{\ell=1}^{p-1}(pq+\ell)\equiv(p-1)!\pmod{p^2},
\]
because the linear correction involves
\(\sum_{\ell=1}^{p-1}\ell^{-1}=0\) in \(\mathbb F_p\).
With \(E_e=(p^e-1)/(p-1)\) it follows by iterating the removal of
multiples of \(p\) from a factorial that
\[
 \frac{(p^e x)!}{p^{E_e x}x!}
 \equiv ((p-1)!)^{E_e x}\pmod{p^2}\qquad(x\ge0).
 \tag{6}
\]
The quotient in (6) is a \(p\)-adic unit. For a row at index \(k\),
\(x=(j-k)a+kb\) is affine in \(k\), with the same slope \(b-a\)
at both orders. Thus after independent powers \(p^{e_i}\),
using the common argument substitution
\[
 t\longmapsto [p(p-1)!]^{-\sum_i E_{e_i}(b_i-a_i)}t
\]
and separate scalars for each order, every coefficient retains its
valuation and its unit modulo \(p^2\). In particular integral
polynomials normalized to minimum valuation zero retain their
coefficient reductions modulo \(p^2\).

For a common shift \(M\) after the dilation, write the row's smaller
exponent as \(a\) and the absolute gap as \(d\). After dividing
out the constant \((ja)!\) before the shift and \((j(a+M))!\)
after it, each row-factor ratio changes by a product of
\(1+jM/(ja+\ell)\) with \(1\le\ell\le jd\). Choose \(M\) divisible
by a power of \(p\) at least two larger than the maximal
\(v_p(ja+\ell)\) in this finite set for \(j=5,6\).
Then every such product is \(1\bmod p^2\).
This gives the same precision for the shifted profile and preserves
the corresponding prime criterion. For example, the prime-five
condition gives profiles with independently chosen signed
power-of-five gaps, and the prime-eleven condition gives the
analogous power-of-eleven gaps.
