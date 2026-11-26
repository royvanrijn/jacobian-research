# When a fixed pure factorial term can dominate in every order

Fix a prime \(p\), a finite support of distinct nonnegative
integer exponent columns, and a distinguished column \(\sigma\).
For a count vector of replacements \(m=(m_\nu)_{\nu\ne\sigma}\)
write \(h=\sum m_\nu\) and
\(D_i=\sum_{\nu\ne\sigma}m_\nu(\nu_i-\sigma_i)\).
The support may include the zero column and may be fully shared.
A row can decrease if \(D_i<0\) for some count vector.

**Theorem.** The following conditions are equivalent.

1. For any nonzero coefficient ratios, the valuation of each
   fixed replacement term divided by the pure \(\sigma\) term
   is bounded below as the integer order \(n\ge h\) varies.
2. Every replacement count vector decreases at most one row,
   and the distinguished exponent in each row that can decrease
   is a power of \(p\), including \(p^0=1\).
3. Some choice of nonzero coefficient ratios makes the pure
   \(\sigma\) term uniquely least in valuation at every positive
   order.

When these conditions hold, increasing all the nonpivot
coefficient valuations by a sufficiently large common amount
gives the conclusion in item 3, with an effective threshold.
If the effective slope of the pivot is strictly smaller than
every other effective slope, its strict dominance at every
order can be decided by finitely many count-vector tests, each
over a finite residue ring. If item 2 fails, every coefficient
choice has arbitrarily late orders at which a fixed competitor
has strictly lower valuation than the pivot.

## The denominator poles

Write \(r_\nu=c_\nu/c_\sigma\) for the coefficient ratios.
The exact factor multiplying the pure term for each fixed
count vector is the rational function
\[
 R_m(z)=
   \frac{z^{\underline h}\prod_\nu r_\nu^{m_\nu}}
        {\prod_\nu m_\nu!}
   \frac{\prod_{D_i>0}\prod_{u=1}^{D_i}(\sigma_i z+u)}
        {\prod_{D_i<0}\prod_{u=0}^{-D_i-1}(\sigma_i z-u)} .
                                                               \tag{1}
\]
This is the [endpoint count formula](FACTORIAL_ENDPOINT_NORMALIZATION_AND_INTEGER_ANCHORS.md).
Its specialization at integers \(n\ge h\) is well defined
since \(-D_i\le h\sigma_i\). Such integers are dense in
\(\mathbb Z_p\).

If a count vector decreases two or more rows, every one of
those denominator products has a simple zero at \(z=0\).
The numerator has only the one zero from \(z^{\underline h}\).
This gives a pole at zero and valuation tending to minus
infinity along \(n=p^k\) sufficiently large.

If a row can decrease and \(\sigma_i=p^b u\), \(u>1\) prime
to \(p\), take sufficiently many copies of a column that
decreases it to make \(-D_i>p^b\). The denominator now has
a zero at \(z=1/u\in\mathbb Z_p\). A nonnegative rational
root of the numerator must be a nonnegative integer from
\(z^{\underline h}\). Thus \(1/u\) is an uncanceled pole.
Integers \(n\ge h\) approaching \(1/u\) have valuations
tending to minus infinity. In both cases the coefficient
ratios are fixed nonzero scalars and cannot remove the pole.

Conversely, assume item 2. For the sole decreasing row in
a count vector let \(q=\sigma_i=p^b\), \(L=-D_i>0\), and
\(r=\lceil L/q\rceil\le h\). Its factors with \(q\mid u\)
are exactly \(q^r z^{\underline r}\), which cancel from
the numerator. The remaining factors are scalar multiples
of \(1-qz/u\) with \(v_p(q/u)\ge1\), which have constant
unit valuation on \(\mathbb Z_p\). Hence there is no pole
on \(\mathbb Z_p\). This proves item 1, and the preceding
two obstructions show item 3 requires item 2 as well.

The support condition is a finite test. Check the pivot
height at each row with a lower entry in some other column.
For each pair of rows, solve for nonnegative count weights
that strictly decrease both. If a real solution exists,
a rational one exists by rational linear feasibility,
and multiplication by a common denominator gives an
integer count vector.

## The uniform comparison at strict slope

Let \(K_i=D_i\) on each increasing row of the count. After the
cancellation above the denominator valuation apart from
the count factorials is
\[
          M=\sum_{\ell=1}^b\lceil L/p^\ell\rceil
                \le \frac{L(1-1/q)}{p-1}+b .
\]
Indeed an index \(0\le u<L\) contributes one at each level
\(p^\ell\mid u\) for \(1\le\ell\le b\), with \(u=0\)
counted at every level. Its canceled scalar contributes
the same \(b\) for each index divisible by \(q\).
Use \(L=r=b=M=0\) when no row decreases.

The remaining falling product at \(n\ge h\) is divisible
by \((h-r)!\). Each positive row interval is divisible by
\(K_i!\). Use \(v(\prod m_\nu!)\le h/(p-1)\), Legendre's
digit sum formula and \(r-L/q\le1\) to get
\[
 v(R_m(n))\ \ge\
   \sum_\nu m_\nu\left(v(r_\nu)+\frac{|\nu|-|\sigma|}{p-1}\right)
              -C(1+\log_p(1+h)).                               \tag{2}
\]
Here \(C\) is computable from the support and \(p\):
all \(K_i\) are bounded by a fixed constant times \(h\),
\(h-r\le h\), and \(b\) has a fixed maximum over rows.
This bound applies uniformly for every \(n\ge h\).

A common increase in the \(v(r_\nu)\) increases the linear
term in (2) by that amount times \(h\), and
\((1+\log_p(1+h))/h\) is bounded for \(h\ge1\).
This proves item 3 and the effective coefficient threshold.
If every slope gap
\(v(r_\nu)+(|\nu|-|\sigma|)/(p-1)\) is positive, (2) is
strictly positive for \(h\) above a computable threshold.
There are finitely many count vectors at or below it.
For one such vector the cancellation expresses the
valuation as \(\kappa+v_p(P(n))\) with a rational
constant \(\kappa\) and an integer polynomial \(P\):
take the remaining falling and rising products, and
put the coefficient valuations, count factorials and
the scalar denominator valuations into \(\kappa\).
If \(\kappa>0\) the comparison always passes. Otherwise
check \(P(n)\equiv0\bmod p^{\lfloor-\kappa\rfloor+1}\)
in every residue class. Each class contains integers
\(n\ge h\), so this is an exact all-order decision.
The [digit-window test](FACTORIAL_ALL_ORDER_PURE_TERM_VALUATION_TEST.md)
also allows slope ties and a different pure term at
different orders.

## A uniform dyadic unit family with arbitrarily many terms

**Proposition.** At a place over two, suppose every pivot
exponent is even, each decreasing row has pivot exponent two,
and every replacement count vector decreases at most one row.
Suppose also that \(|\nu|-|\sigma|\ge d+1\) for every other
column and that the coefficient ratios are integral. Then the
moment normalized by the pure \(\sigma\) term is a unit at every
positive order.

A replacement count of size \(h>0\) increases the total
exponent by at least \((d+1)h\). With no decrease, some rising
product has length at least two, so it supplies a factor of
two and the multinomial is integral. If one row decreases by
\(L>0\), the positive changes \(K_i\) sum to at least
\((d+1)h+L\). Each rising product starts at an odd value,
since its pivot exponent is even. It supplies at least
\(\lfloor K_i/2\rfloor\) factors of two. There are at most
\(d-1\) positive rows, the decreasing row loses at most
\(\lceil L/2\rceil\), and the count factorials lose at most
\(v_2(h!)\le h-1\). The resulting integral lower bound is at least
\[
     \frac{(d+1)h+L-(d-1)}2-\frac{L+1}2-(h-1)
                 =\frac{(d-1)(h-1)+1}{2}>0 .
\]
The remaining falling product is integral. Thus every
nonpure term vanishes in the residue of the normalized sum.

For example, \(x^2y^2+cxy^6+d x^6y\) works with integral
\(c,d\). More generally, one can take a pivot \((2,\ldots,2)\)
in any dimension at least two and arbitrarily many distinct
columns that decrease one row by one and increase every
other row by at least four. Indeed two decreases would require
more than half of the replacements to decrease each of two
different rows, which is impossible. Such columns are all
fully shared and their degree gaps are at least \(d+1\).
