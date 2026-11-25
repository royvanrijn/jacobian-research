# Mixed unit-gap Newton widths at prime-plus-one terminals

For distinct exponent columns \(a,b\) with nonnegative integer entries,
write the ordinary factorial moment as
\[
 M_j(t)=\sum_{k=0}^j\binom jk
             \prod_i ((j-k)a_i+kb_i)!\,t^k .
\]
Equal rows contribute only a scalar and can be removed. At an odd prime
\(p\), the results below concern the pair \(M_p,M_{p+1}\). They give a
valuation-only criterion when unit gaps occur in both orientations.
The criterion is preserved under adding certain complementary pairs
of arbitrary-width exponent rows.

## A width criterion for mixed unit gaps

Assume each retained row has gap \(1\) or \(-1\). Denote its lower
exponent by \(c_i\) and its orientation by \(+\) or \(-\), respectively.
Both orientations are present. Put
\[
 R=\sum_{i\in+}(1+v_p(c_i+1)),\qquad
 S=\sum_{i\in-}(1+v_p(c_i+1)).
\]
For a nonnegative integer \(\ell\le p+1\), let
\[
 H_i(\ell)=\sum_{h=1}^{\ell}v_p((p+1)c_i+h),\qquad
 V(k)=\mathbf 1_{\{2\le k\le p-1\}}
              +\sum_{i\in+}H_i(k)+\sum_{i\in-}H_i(p+1-k).
\]
For an integer \(\lambda\), let \(W(\lambda)\) be the largest minus
the smallest minimizing index of \(V(k)-\lambda k\) over
\(0\le k\le p+1\). There are at most three widths to evaluate. If
\[
 \min\{W(0),p-\mathbf 1_{\{S>1\}}-\mathbf 1_{\{R>1\}}\}
 +\mathbf 1_{\{S>1,\ W(1-S)>0\}}
 +\mathbf 1_{\{R>1,\ W(R-1)>0\}}\ \le 1,
 \tag{1}
\]
then \(M_p\) and \(M_{p+1}\) are coprime. This criterion allows every
residue class of the lower exponents.

Here is a simpler sufficient test. Suppose no \(c_i\equiv-1\pmod p\).
Let \(h_i\in\{2,\ldots,p\}\) be the unique index with
\(p\mid(p+1)c_i+h_i\), and put \(w_i=v_p((p+1)c_i+h_i)\).
If there are \(r\) increasing and \(s\) decreasing rows, require
\[
 \sum_{i\in+} w_i\le2r-2\quad\text{when }r\ge2,
 \qquad
 \sum_{i\in-} w_i\le2s-2\quad\text{when }s\ge2.
 \tag{2}
\]
Then only \(W(0)\) is needed. It is sufficient that the minimum
of \(V\) is attained at a single index or at two adjacent indices.

In particular, for orders five and six, take \(r,s\ge2\), with
every \(c_i\equiv0\) or \(1\pmod5\) and at least one
\(c_i\equiv1\pmod5\). The two moments are coprime whenever (2)
holds. The bounds in (2) hold if all hits are simple. That condition
is \(c_i\not\equiv20\pmod{25}\) for the bases divisible by five,
and \(c_i\not\equiv16\pmod{25}\) for the bases congruent to one.
For example, one can take increasing bases \(1,6\) and decreasing
bases \(5,10\). The families allow arbitrarily many rows and arbitrarily
large exponents.

### Proof

Divide \(M_j\) by \(\prod_i (j c_i)!\). Each increasing row then
contributes the product of the first \(k\) integers above \(j c_i\),
and a decreasing row contributes the analogous product of length
\(j-k\). At \(j=p\), the only possible hit in a row's length-\(p\)
product is at position \(p\), with valuation \(1+v_p(c_i+1)\).
The interior binomial coefficients have valuation one. The ordinates
of the lower-order polynomial are consequently
\[
 (S,1,\ldots,1,R).
\]
Its lower Newton polygon has a width-one slope \(1-S\) if \(S>1\),
a width-one slope \(R-1\) if \(R>1\), and a horizontal edge of width
\(p-\mathbf 1_{\{S>1\}}-\mathbf 1_{\{R>1\}}\). At \(j=p+1\)
the binomial valuations are one at \(2\le k\le p-1\) and zero
elsewhere, giving exactly the array \(V\).

For a polynomial \(P(t)\) over \(\mathbb Q\) and an integer \(\lambda\),
its width at \(\lambda\) is the width of the valuation-minimum
indices for \(P(p^{-\lambda}t)\). After division by a coefficient
of minimum valuation, the two endpoints are the first and last
nonzero coefficients of its reduction modulo \(p\). These indices
add under multiplication of nonzero reduced polynomials. Hence the
width at each integer slope adds under multiplication. Widths at
distinct slopes occupy disjoint index intervals apart from shared
endpoints, since minimizing indices move to the right as \(\lambda\)
increases. In a factorization of \(M_p\), the widths at the displayed
integer slopes already sum to \(p\). They therefore sum to the degree
in each factor. This gives the left side of (1) as an upper bound on
the degree of any rational common factor. This is the Newton polygon
factorization principle in the constant-to-leading orientation
used here. [Filaseta and Rocks, §2](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf)

For every real \(t\), the Gamma integral gives
\[
 M_{p+1}(t)=\int_{(0,\infty)^d}(x^a+t x^b)^{p+1}
               e^{-\sum_i x_i}\,dx>0.
\]
The integrand is nonnegative and the monomial combination, whose
two exponent columns are distinct, is nonzero on an open set.
Thus every rational divisor of \(M_{p+1}\) has even degree. An upper
bound of one on the gcd degree proves coprimality over \(\mathbb Q\)
and a rational Bézout identity transfers this to every
characteristic-zero field.

Under the hypotheses of (2), \(R=r\), \(S=s\), and each upper-order
row has exactly one p-divisible factor in the full length-\((p+1)\)
product. Its index is neither \(1\) nor \(p+1\).
The first two ordinates of \(V\) are both \(\sum_{i\in-} w_i\)
and the last two both \(\sum_{i\in+} w_i\). Every ordinate is at
least one. Thus every lower-polygon slope is at least
\(-(\sum_{i\in-}w_i-1)/2\) and at most
\((\sum_{i\in+}w_i-1)/2\). For the lower bound, use the least slope
from the left endpoint, with numerator zero at index one and with
any negative numerator at a distance at least two. The upper bound
follows from the analogous slopes into the right endpoint. If
\(s\ge2\), (2) puts the lower bound strictly above \(1-s\).
For \(r\ge2\) the upper bound is strictly below \(r-1\).
Consequently (1) follows from \(W(0)\le1\).

For the claimed terminal-six specialization, let \(R_4,S_4\)
be the sums of hit weights at index four in the two orientations,
and \(R_*,S_*\) their total hit weights. The hits are at index
four when \(c_i\equiv1\pmod5\), and five otherwise. This gives
\[
 V=(S_*,S_*,1+S_4,1,1+R_4,R_*,R_*).
\]
Both endpoint weights are at least two. Index three is a minimum
and at least one of indices two and four has greater valuation,
so \(W(0)\le1\). The residues for simple hits follow by testing
\(6c_i+h_i\) modulo \(25\).

## Preservation under changing or adjoining rows

Independently multiplying the two exponents in a row by \(p^e\)
adds an affine function of \(k\) to its valuation contribution
at each order. Indeed, for \(x\ge0\) Legendre's formula gives
\[
 v_p((p^e x)!)=v_p(x!)+\frac{p^e-1}{p-1}\,x .
\]
At \(x=(j-k)a_i+kb_i\), the slope of the added affine function
is independent of \(j\). It translates the slopes of the two
Newton polygons by the same amount. Therefore any criterion
established above is preserved by arbitrary rowwise dilations
by powers of \(p\).

For a second operation, take positive \(A,B\) with \(A\ne B\).
Choose \(P=p^m>(p+1)\max(A,B)\) and adjoin the complementary rows
\[
 (A,B),\qquad(P-A,P-B).
\]
At index \(k\) of order \(j\in\{p,p+1\}\), set
\(x=(j-k)A+kB\), so \(0<x<P\). The two factorial arguments are
\(x,jP-x\). The identity
\[
 \binom{jP}{x}=\frac{jP}{x}\prod_{u=1}^{x-1}\frac{jP-u}{u}
\]
shows \(v_p\binom{jP}{x}=v_p(jP)-v_p(x)\), since
\(v_p(u)<m\le v_p(jP)\) for \(0<u<P\). Thus the sum of the two
factorial valuations is a constant in \(k\) plus \(v_p(x)\).
In particular, if
\[
 v_p(B-A)>v_p(A)+1,
\]
then \(v_p(x)=v_p(jA)\) at every index of both moment orders.
The two rows only add a scalar baseline to each valuation array
and preserve the coprimality criterion. One can adjoin arbitrarily
many such pairs. At terminal six they include arbitrary positive
unequal \(A,B\) with \(A\equiv B\not\equiv0\pmod5\) and
\(A\equiv B\pmod{25}\), so the resulting profiles can contain
arbitrarily large gaps that are not powers of five.

The width test uses neither residual polynomial calculations nor a
finite search. A bound of two or larger in (1) leaves open whether
the available Newton edges support a common factor.
