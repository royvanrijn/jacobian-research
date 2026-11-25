# Shared terminal-eight separation with outer mean at least .95

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are all shared. Let \(u\) be the order-seven outer geometric mean
after normalizing both endpoint coefficients to one and removing
binomial coefficients. If \(u\ge19/20\), the order-seven and
order-eight factorial moment polynomials have disjoint root
arguments under independent positive rescalings. The dimension and
exponents are arbitrary. Rational Bézout gives coprimality over every
characteristic-zero field.

## Scaling and matched reciprocals

Put \(p=1-u=\lambda^2\), \(q=1-w=rp\), with \(w\) the order-eight
outer mean. Use \(0\le\lambda\le9/40\),
\(1019/1000\le r\le128/125\), and \(q<.055\). All decimals in
this proof denote exact rationals. The
[linked outer-mean comparison](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md)
gives
\[
 w\le u(1+(1-u^2)/48)^{-1/2},\qquad
 r\ge1+\frac{u(1+u)}{96(1+(1-u^2)/48)^{3/2}} .
\]
For a shared row of lower exponent \(\ell\), gap \(d\), and Gauss
index \(1\le h\le d\), let \(x=d/(7\ell+h)\) and \(x'=d/(8\ell+h)\).
Then \(x\le(u^{-2}-1)/6<.0185\), and
\[
       \frac{7x'}{6x}\le\frac{49}{6(8-.0185)} .
\]
The power comparison on the outer product, followed by the tangent
bound \(u^\gamma\ge1-\gamma(1-u)\) for \(\gamma\ge1\), gives
the upper bound on \(r\). The checker verifies both endpoints.

Use \(t=(\lambda X-1)/(\lambda X+1)\) and multiply the normalized
order-\(j\) polynomial by \((1+\lambda X)^j/(2^j\lambda^j)\).
For the lower order use the septic coefficients of the
[terminal-seven near-flat proof](MULTIVARIATE_TERMINAL_SEVEN_NEAR_FLAT.md)
with ratio argument equal to one. That enclosure's within-order
deficit bounds apply, since the septic deficit is now \(p\).

## Octic linked means and skews

Write the normalized multipliers as
\[
 (1,me^\alpha,ye^\beta,ze^\gamma,c,ze^{-\gamma},
                                ye^{-\beta},me^{-\alpha},1),
 \qquad m=w.
\]
For a single increasing factor put \(g_h=\log(1+h x')\) and
\(s=g_7/2\). The normalized log multiplier at \(k\) is
\(\sum_{h=0}^{k-1}g_h-(k/8)\sum_{h=0}^{7}g_h\).
Subtracting reflections gives the single-factor links
\[
 L=5\alpha-7\gamma
 =\Delta^4g_0+\tfrac52\Delta^4g_1+\tfrac52\Delta^4g_2+\Delta^4g_3,
 \qquad
 K=3\alpha-7\beta+7\gamma=(\Delta^6g_0+\Delta^6g_1)/4 .
\]
Thus locally \(|L|\le42(x')^4\) and \(|K|\le60(x')^6\), by
integrating the respective derivatives of \(\log(1+x't)\) over
unit cubes. The local \(\gamma\) is the sum of
\(g_h+g_{7-h}-g_3-g_4\) for \(h=0,1,2\), divided by eight.
Pairing each product shows its magnitude is at most
\((12+6+2)(x')^2/8\). Reflections only reverse the signs.
The sum of all order-eight reciprocals is at most
\((w^{-2}-1)/7\le.312q\). Summing the bounds therefore gives
\[
 |\gamma|\le.244q^2,\quad |L|\le.399q^4,\quad |K|\le.056q^6.
\]
Using \(\alpha=(7\gamma+L)/5\) and
\(\beta=(8/5)\gamma+(3/35)L-K/7\), the checker also obtains
\(|\alpha|\le.343q^2\) and \(|\beta|\le.391q^2\).

With \(G_a(s)=as-\operatorname{artanh}(a\tanh s)\), averaging the
reflected log coefficients gives
\[
\begin{aligned}
 H&=\log(y/w^{12/7})=\sum G_{5/7}(s),\\
 J_3&=\log(z/w^{15/7})-2H=\sum (G_{3/7}-G_{5/7})(s),\\
 J_4&=\log(c/w^{16/7})-\tfrac{12}{5}H
      =J_3+\sum(G_{1/7}-\tfrac25G_{5/7})(s).
\end{aligned}
\]
The derivative is
\(G'_a(s)=a(1-a^2)\tanh^2(s)/(1-a^2\tanh^2(s))\).
The cubic leading terms cancel in each difference. In both
differences the first derivative is nonpositive, since the
numerators after normalizing the factors are equal and the
second term has the smaller denominator. To bound their
absolute values, use \(\tanh s\le s\) and put
\(s_* = q/(1-q)\ge\sum s\). The derivative magnitude for
\(G_{3/7}-G_{5/7}\) is bounded above by
\[
 \frac{(120/343)(16/49)}{(1-(25/49)s_*^2)(1-(9/49)s_*^2)}\,s^4 .
\]
For \(G_{1/7}-(2/5)G_{5/7}\) use the same expression with
\((120,16,9)\) replaced by \((48,24,1)\).
Integration and \(\sum s^5\le s_*^5\) give the verified caps
\[
      0\le H\le.14q^3,\qquad -.031q^5\le J_3\le0,\qquad
      -.05q^5\le J_4\le0 .
\]
For the first cap bound \(G'_{5/7}\) by
\((120/343)s^2/(1-(25/49)s_*^2)\).
All scalar caps are checked by rational arithmetic.

## Stable Cayley families

Let \(e=(0,1,12/7,15/7,16/7)\) and let the pure five means be
\(w^{e_i}\), \(0\le i\le4\). The actual log corrections to those
means are
\[
       D=(0,0,H,2H+J_3,\tfrac{12}5H+J_4).
\]
Put \(a=(0,7/5,8/5,1,0)\),
\(b=(0,1/5,3/35,0,0)\),
and \(d=(0,0,-1/7,0,0)\).
The reflected-pair skew at index \(i\) is
\(a_i\gamma+b_iL+d_iK\). Each even Cayley numerator uses the
means times \(\cosh\) of their skew; each odd numerator uses
the means times \(\sinh\) of their skew.

In the even part expand a corrected multiplier as
\[
 w^{e_i}\big[1+D_i+\tfrac12a_i^2\gamma^2\big]
        +q^6\,\epsilon_i .
\]
In the odd part use
\[
 w^{e_i}(a_i\gamma+b_iL+d_iK)+q^5\,\delta_i .
\]
The absolute bounds for \(\epsilon_i\), in order of \(i=0,\ldots,4\),
are \((0,.028,.027,.042,.058)\), and those for \(\delta_i\)
are \((0,.001,.057,.070,0)\).
The scalar check uses the second-order Taylor remainder for \(e^D\),
the difference between the actual skew square and \(a_i^2\gamma^2\),
and the product of the shape defect and the hyperbolic correction.
On \(|s|\le.01\) the series give
\(\cosh s-1\le.501s^2\),
\(|\cosh s-1-s^2/2|\le s^4/20\), and
\(|\sinh s-s|\le |s|^3/5\). The replay checks that all skew
arguments satisfy this bound and verifies each resulting error cap.

The code transforms the following seven shared parameters
separately, with their displayed leading powers of \(q\).
\[
\begin{array}{c|c|c}
\text{parameter}&\text{weight vector across }i&\text{factor}\\ \hline
 H/q^3&(0,0,1,2,12/5)&q^3\\
 J_3/q^5&(0,0,0,1,0)&q^5\\
 J_4/q^5&(0,0,0,0,1)&q^5\\
 \gamma^2/q^4&(0,49/50,32/25,1/2,0)&q^4\\
 \gamma/q^2&a&q^2\\
 L/q^4&b&q^4\\
 K/q^6&d&q^6
\end{array}
\]
The directions of all shared contributions are kept intact
when evaluating each projected polynomial.

A Cayley coefficient is obtained by the binomial expansion
of \((Z-1)^i(Z+1)^{8-i}\) with \(Z=\lambda X\), then pairing
\(i\) and \(8-i\). Within any of the seven contributions it
is thus a linear combination of the powers \((1-q)^{e_i}\).
The checker sums their first nine Taylor coefficients exactly,
removes the leading zero coefficients, and evaluates the
remaining polynomial with interval arithmetic. If the first
nonzero coefficient is in degree \(d\), the tail after dividing
by \(q^d\) is bounded in absolute value by
\[
 \frac{q^{9-d}}{1-q}
       \sum_i |c_i\tbinom{e_i}{9}|,
\]
where \(c_i\) are the Cayley weights of that combination.
Successive absolute binomial coefficients decrease after
degree nine since \(0\le e_i\le16/7\). The associated powers of
\(q\) are replaced by \(r\lambda^2\) before division by the
Cayley scaling, and every surviving exponent of \(\lambda\)
is checked to be nonnegative. This certifies the central
cancellations without interval subtraction of nearly equal
unscaled terms.

## Ray checks and anchor

Leading dominance on \(|X|=2\) bounds all the scaled roots
strictly inside that disk, putting the physical roots in the
left half-plane. A ray with \(\theta=\pi-2\arctan(\lambda b)\)
has the equation
\((Y-b)(1+\lambda^2bY)=-b\lambda^2x^2\) for \(X=x+iY\).
The checker encloses \(Y\) with three updates and excludes roots
by strictly signed projections
\(\Im(P_j(X)\overline X^{\,h})\), \(0\le h\le j\).
As every opposite-parity parameter and remainder can be
reflected in sign, it suffices to check \(x\ge0\).
The six ray parameters interpolate affinely in \(p\).
The replay used \(2845,1049,4055,1351,837,455\) boxes in order
seven and \(26945,3551,28767,2345,3463,705\) in order eight.

To reach the flat anchor multiply all matched Gauss reciprocals
by a common parameter decreasing to zero. The factor ratios
stay fixed, so the certified bounds on \(r\) continue to hold.
The scaled limits are the imaginary Hermite polynomials
\[
\begin{aligned}
 &X^7+\tfrac74X^5+\tfrac{35}{48}X^3+\tfrac{35}{576}X,\\
 &X^8+2rX^6+\tfrac{15}{14}r^2X^4
                   +\tfrac{15}{98}r^3X^2+\tfrac{15}{5488}r^4 .
\end{aligned}
\]
These also follow by taking the zero-deficit values of the
stable coefficient formulas. The checker verifies exact brackets
for the positive squared imaginary coordinates, dividing those
of order eight by \(r\), and compares them with the squared
ray heights at the anchor. The successive upper sectors contain
orders \(8,7,8,7,8,7,8\); the additional septic root is real.
The verified exclusions preserve these counts, counting the
last sector together with its conjugate.

The sole septic root in that joined sector is real by conjugacy.
The actual octic is strictly positive on the real axis by its
eighth-power factorial moment integral for distinct monomials.
Thus its two roots in the sector are nonreal and have arguments
distinct from the real septic root.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_eight_near_flat.py
~~~
