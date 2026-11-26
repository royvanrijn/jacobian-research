# Bilinear three-windows through the first dyadic cross threshold

Let \(f=BXY+dX+eY\) with \(Bde\ne0\) over a characteristic-zero
field. Write \(P_j=\Lambda_2(f^j)/j!\), with
\(\Lambda_2(X^iY^j)=i!j!\). At a dyadic place,
assume \(0<\min(v(d),v(e))-v(B)\le1/2\), where \(v(2)=1\).
Then no three consecutive positive moments are zero.
The proof works over an arbitrary dyadic residue field,
including ramified extensions.

The [exact bilinear recurrence](ZERO_CONSTANT_BINARY_QUADRATIC_MATCHING_AND_BOUNDARY_LIFTS.md#a-valuation-one-lift-for-the-bilinear-three-term-phase)
and the forcing from three zeros ending at \(k\) give
\(d^k+e^k=0\). In particular \(v(d)=v(e)\).
Scale by \(d\), and put \(q=e/d,\ b=B/d,\ s=-v(b)\),
\(\sigma=1+q\) and \(H(T)=(1-T)(1-qT)\).
Thus \(q^k=-1\) and \(0<s\le1/2\).
Use the exact series
\[
         \sum_{n\ge0}P_n T^n
                  =\sum_{j\ge0}j!b^j T^j H(T)^{-j-1}.      \tag{1}
\]
We use \(h_m(1,q)=1+q+\cdots+q^m\), with a negative
coefficient index interpreted as zero.

## Before the threshold

If \(0<s<1/2\), the unique lowest cross weight in (1)
is at \(j=1\), by the
[bitwise cross weight formula](BINARY_ZERO_CONSTANT_QUADRATIC_WEIGHTED_MATCHING_LOCAL.md#exact-cross-weights-below-the-unit-valuation).
After dividing by \(b\), the leading residue at the odd
indices \(2m+1\) is \(h_m(1,\bar q)^2\).
It cannot vanish at two successive values of \(m\),
so three zeros cannot end at odd \(k\).

At an even index \(n\), every odd \(j\) has an extra
factor \(2\) in \([T^{n-j}]H^{-j-1}\), since \(j+1\)
is even. Those weights become strictly positive,
while among even \(j\) the unique smallest weight
is zero at \(j=0\). Thus the residue of \(P_k\)
at an even \(k\) is \(h_k(1,\bar q)=1\):
if \(\bar q\ne1\) use \(\bar q^k=1\), and otherwise
use \(k+1\) odd. This excludes even \(k\) too.

## The first tied layer

Suppose \(s=1/2\), and set \(u=2b^2\), a unit.
The minimizing cross groups are \(j=1,3\).
For odd \(k=2m+1\), use \(k\) if \(m\) is even,
and \(k-2\) if \(m\) is odd. At the selected index
\(2l+1\), \(l\) is even. The \(j=3\) contribution
in the leading residue is absent, since
\(2l-2\) is not a multiple of four in \(\overline H^{-4}\).
The \(j=1\) contribution is \(h_l(1,\bar q)^2\ne0\).
If \(\bar q=1\), this follows from \(l\) even;
otherwise \(\bar q^{2(l+1)}\) is \(\bar q\) or
\(\bar q^{-1}\), using \(\bar q^k=1\).

Now \(k=2m\) is even. If \(\bar q=1\), the leading
residue at odd indices, after dividing by \(b\), is
\[
                 \frac{T(1+\bar u T^2+T^4)}{1+T^8}.
\]
Its zero at \(k-1\) requires \(k\equiv0\bmod8\).
At the even \(k\), the lowest possible terms are
\(j=0,2\), of weight zero. The \(j=2\) term
has factor \([T^{k-2}](1+T^2)^{-3}
=\binom{m+1}{2}=0\) in the residue field, while
the \(j=0\) term is one.
All odd terms are higher at the even index because
their coefficients have a factor \(2\).
Thus this residue is excluded.

Let \(\bar q\ne1\). At even \(k\), the zero condition
on the two weight-zero terms is
\(1+\bar u[T^{k-2}]\overline H^{-3}=0\).
Use \(\overline H^{-3}=\overline H/\overline H^4\).
If \(k\equiv0\bmod4\), the indicated coefficient is
\(\bar q\,h_{k/4-1}(1,\bar q^4)=0\), by \(\bar q^k=1\),
a contradiction. If \(k\equiv2\bmod4\), it is
\(h_{(k-2)/4}(1,\bar q^4)=1/\bar\sigma^2\).
These evaluations follow by summing the geometric
series and using \(\bar q\ne1\).
The remaining case has
\[
            m\ \text{odd},\quad m\ge3,\quad
        \bar q^m=1,\quad \bar u=\bar\sigma^2.                \tag{2}
\]
The inequality \(m\ge3\) follows because a
three-zero block has \(k\ge3\), and \(\bar q^m=1\)
follows by Frobenius from \(\bar q^{2m}=1\).

## The next tied precision

Evaluate \(P_{k-1}/b\) at valuation one.
For brevity let \(c^{(a)}_r=[T^r]H(T)^{-a}\).
The \(j=0\) term has valuation at least \(3/2\):
\(h_{k-1}(1,q)=2/(1-q)\), and \(v(b^{-1})=1/2\).
For \(j=2\), the odd coefficient
\(c^{(3)}_{k-3}\) is zero modulo \(2\mathcal O\):
the odd part of \(\overline H/\overline H^4\)
has indices \(1\bmod4\), whereas \(k-3\equiv3\bmod4\).
Thus this term also has valuation at least \(3/2\)
in \(P_{k-1}/b\).
For every even \(j\ge4\) the weight in \(P_{k-1}\)
is at least one and increases by \(1/2\) after division
by \(b\). For odd \(j\ge9\) the weight increases to
at least three, by the same bitwise formula.
It remains to evaluate \(j=1,3,5,7\).

For the \(j=1\) term we have the full congruence
\[
                      c^{(2)}_{k-2}
                               \equiv \frac{2}{1-q^2}
                                      \pmod{4\mathcal O}.    \tag{3}
\]
Indeed let \(S=T^2\), \(K=1+\sigma^2 T^2+q^2T^4\),
and \(J=-\sigma T+qT^2-\sigma qT^3\).
Then \(H^2=K+2J\), so \(H^{-2}\equiv K^{-1}-2JK^{-2}\)
modulo \(4\mathcal O\).
The correction at \(T^{k-2}\) vanishes modulo
\(4\mathcal O\) because the even part of \(J\) is
\(qT^2\) and \(K^{-2}\) has indices divisible by
four modulo \(2\mathcal O\).
Set \(K_0(S)=(1-S)(1-q^2 S)\).
Then \(K-K_0=2(1+q+q^2)S\).
The difference of \(K^{-1}\) and \(K_0^{-1}\)
at \(S^{m-1}\) also vanishes modulo \(4\mathcal O\),
because its correction uses the odd coefficient
at \(S^{m-2}\) of \(K_0^{-2}\).
The resulting coefficient is
\(h_{m-1}(1,q^2)=2/(1-q^2)\) by \(q^{2m}=-1\).

For the \(j=3\) term, \(H^{-4}\equiv K_0(S)^{-2}\)
modulo \(4\mathcal O\), by squaring the previous
denominator comparisons. The odd coefficient at
\(S^{m-2}\), divided by two, reduces to
\[
 \bar\sigma^2 [S^{m-3}]
       \frac{1+\bar q^2 S^2}
            {(1-S^2)^2(1-\bar q^4 S^2)^2}
                            =\frac1{\bar\sigma^4}.         \tag{4}
\]
The expression uses the odd numerator terms of
\((1+S)^2(1+q^2 S)^2\) over the displayed denominator.
To check the last value, if \(m=4a+1\) the coefficient
after removing \(\bar\sigma^2\) equals
\(\bar q^2h_{a-1}(1,\bar q^8)\), using
\(\bar q^{8a}=\bar q^{-2}\).
If \(m=4a+3\), it equals \(h_a(1,\bar q^8)\),
using \(\bar q^{8(a+1)}=\bar q^2\).
Both yield (4) by the geometric formula.
Thus after division of \(P_{k-1}/b\) by two, the
\(j=1,3\) residues are \(1/\bar\sigma^2\) and
\(\bar u/\bar\sigma^4=1/\bar\sigma^2\), which cancel.

The \(j=5\) multiplier after the two divisions is
\(60 b^4=15u^2\), with residue \(\bar u^2\).
Its \(\overline H^{-6}\) coefficient at \(k-6\)
is the square of \([T^{m-3}]\overline H^{-3}\).
This last coefficient is \(1/\bar\sigma^3\):
for \(m=4a+1\) it is
\(\bar q h_{a-1}(1,\bar q^4)\), and for \(m=4a+3\)
it is \(h_a(1,\bar q^4)\).
Here \(\bar q^{4a}=\bar q^{-1}\) and
\(\bar q^{4(a+1)}=\bar q\), respectively.
Thus the \(j=5\) residue is \(1/\bar\sigma^2\).
For \(j=7\), the multiplier after the divisions is
\(2520b^6=315u^3\), but the coefficient of
\(\overline H^{-8}\) is zero, since \(k-8\)
is not divisible by eight. All four contributions
are now accounted for. The nonzero \(j=5\) residue
gives \(v(P_{k-1}/b)=1\), a contradiction.

This covers the bilinear three-term family whenever
the dyadic cross dominance lies in \(0<s\le1/2\).
The opposite local ordering is covered by the
[valuation-one lift](ZERO_CONSTANT_BINARY_QUADRATIC_MATCHING_AND_BOUNDARY_LIFTS.md#a-valuation-one-lift-for-the-bilinear-three-term-phase).
The stronger cross dominance \(s>1/2\) remains open
in this argument.
