# A second axis pairing at reciprocal dyadic endpoints

Use the notation and hypotheses of [final-even nonvanishing](BINARY_EVEN_QUADRATIC_CLOSE_AXIS_FINAL_EVEN.md). In particular \(d=2^e\ge4\), both scaled axis residues are one, \(2d\mid N\), and
\[
       r=1/d,\qquad \tau=v_2(N/(2d)),\qquad h=v((c-a)/a).
\]
We retain the abbreviations \(D=2d,J=N-1\), \(\gamma=(c-a)/a\) and \(\epsilon=(b-a-c)/a\).

**Theorem.** A four-zero block is excluded in each of the following cases.

- \(\tau\ge e-4\).
- \(2/d<h<4/d\), except at \(h=(e-\tau-3+r)/(d-1)\). That strength is in the displayed interval only for \(e-\tau=5,6\).
- \(e-\tau\ge5,\ h\ge2/d\), and
\[
  h>\max\left\{\frac{3/2+r}{2^{\tau+5}-1},
                    \frac{e-\tau-3+r}{d-1}\right\}.
\]
For \(h\ge3/d\) one can replace \(3/2+r\) by \(1+r\) in the first numerator.

Together with the [small-\(h\) exclusion](BINARY_EVEN_QUADRATIC_CLOSE_AXIS_ENDPOINT.md), the first case closes the full branch at the listed \(\tau\). If \(e-\tau=5\) or \(6\), only \(h=2/d\) and \(h=(e-\tau-3+r)/(d-1)\) can survive.

## The boundary cancels the next two odd pairs

Assume the four moments vanish and \(3/d\le h\le1/2\). This section needs \(d\ge8\). We show that in the radial expression for \(M_N\), the odd degree contributions with each fixed power \(m<d/2\) and each fixed number of quadratic picks gain at least two beyond \(\eta_{\mathrm{odd}}+\beta(m)\), where \(\eta_{\mathrm{odd}}=2-r+\tau\) and \(\beta(m)=mh-v_2(m+1)\) as in the final-even proof.

Let \(a_0=d-1\), \(b_0=d/2-1\), and continue to use \(G(Z)=(1+u_* Z^d)P(Z)+E(Z)\) from that proof. For \(0\le i<d\), let \(u_i\) be the paired coefficient ratio for the indices \(i,i+d\) in \(G\), so that \(u_*=u_{a_0}\). They are odd rational units since the scale satisfies \(t^d=2\). The coefficient identity
\[
 \frac{\binom J{i}(2i)!}{\binom J{i-1}(2i-2)!}
                =2(J-i+1)(2i-1)
\]
shows that \(u_i/u_{i-1}\) is one modulo four except at \(i=d/2\), where it is three. Indeed its first factor is \(1-d/(N-i)\) and its second is \(1+2d/(2i-1)\), with \(v_2(N-i)=v_2(i)\) for \(0<i<d\). Thus \(u_i\equiv u_*\pmod4\) above or at \(d/2\), and \(u_i\equiv-u_*\pmod4\) below \(d/2\).

The coefficient \(\lambda_B\) of \(Z^{b_0}\) in \(P\) has valuation \(1/2\), since it omits exactly the bit of value \(d/2\) from \(a_0\). All lower indices below \(d/2\) other than \(b_0\) have gap at least \(1/2+3/4\) in the first block. This follows by omitting that bit and one other profitable bit, whose cost is at least \(3/4\). In the high blocks the gap is at least one. As every paired error gains one from subtracting odd rational units, it follows that
\[
                     E(Z)\equiv
             \lambda_B(u_{b_0}-u_*)Z^{b_0+d}                  \tag{1}
\]
with coefficient errors of valuation at least two.

Write \(P[X,Y]=(P(X)-P(Y))/(X-Y)\) and use the same notation for other divided differences. The boundary equation gives \(v(G[X,Y])>2\), since \(G(Y)\) has valuation greater than one by the preceding proof, and \(g_N(tZ/4)-1\) has coefficients of valuation at least one. The divided difference of \(Z^d\) has valuation at least
\[
            \sum_{i=0}^{e-1}\min(1,2^ih)\ge 1+(d/2-1)3/d>2 .
\]
The first term on the right of (1) has coefficient of valuation \(3/2\). In its divided difference, the ratio of the power divided differences at exponents \(b_0+d\) and \(a_0\) is \(X^{d/2}\) times a unit congruent to one to valuation at least \(1/2\). In fact \(Y/X=1+\gamma\), both exponents are odd, and the two numerator powers differ by \(d/2\). The relative error therefore has valuation at least \(v((1+\gamma)^{d/2}-1)-h\ge1-h\ge1/2\). Replacing the leading divided difference at degree \(a_0\) by \(P[X,Y]\) also gains at least \(1/2\). It follows that
\[
                 1+u_*X^d\equiv
                 -\lambda_B(u_{b_0}-u_*)X^{d/2}              \tag{2}
\]
to valuation at least two.

In the radial formula at \(N\), let \(C_k=\binom Nk(2k+1)!a^k\) and \(Q_k=C_{k+d}/C_k\). The factorial shift in the preceding proof gives \(Q_k\equiv u_{k-1}X^d\pmod4\) in the first block. At \(k=a_0,b_0\), moving from \(u_{k-1}\) to \(u_k\) preserves the congruence since \(k\) is odd. Also
\[
 \frac{C_{b_0}}{C_{a_0}}=
 \lambda_B X^{-d/2}\,
 \frac{2b_0+1}{2a_0+1}\frac{N-a_0}{N-b_0},
\]
and the final two ratios are units congruent to one modulo four. For each \(m<d/2\), the two multinomial coefficients of a \(+d\) pair agree modulo four, and the coefficients at the lower indices \(a_0,b_0\) agree modulo two. These assertions follow by expanding \((1+U+V)^d\) and \((1+U+V)^{d/2}\) in total degree less than \(d/2\).

Consider just the \(a_0,b_0\) pairs, divided by \(C_{a_0}\) and their common powers of the linear and quadratic coefficients. The second pair contributes, to precision two, the common parity of the multinomial coefficient times \(\lambda_B(u_{b_0}-u_*)X^{d/2}\), since its extra factor \(C_{b_0}/C_{a_0}\) has valuation \(1/2\) and \(1+Q_{b_0}\equiv (u_{b_0}-u_*) X^d\) to precision \(3/2\). The first pair contributes the negative of this by (2).

In every other odd \(+d\) pair the same precision holds without mixing pairs. In the first block above \(d/2\), its degree costs at least an extra \(1/2\) and \(u_{k-1}\equiv u_*\pmod4\) gives pair gain \(3/2\). Below or at \(d/2\), any odd degree other than \(b_0\) costs at least \(1/2+3/4\), and the pair gains one. High blocks cost at least one and also gain one on pairing. This proves the claimed gain of two.

For even degrees and odd \(m<d/2\) the corresponding lower bound, excluding \(\beta(m)\), is strictly greater than \(4-r+\tau\). If \(v_2(k)=1\), remove the factor \(k/i\) from the multinomial coefficient, where \(i\) is its odd linear count. After division by two the paired coefficients agree modulo four, since \(d/2\) is divisible by four and their remaining multinomial count is below \(d/2\). A first block with \(k-1\ge d/2\) therefore gains \(3/2\) on pairing, in addition to the factor two and \(\eta_N(k)\ge2-2r+\tau\). A first block with \(k-1<d/2\) instead has an extra \(1/2\) axis cost, from omitting the bit of value \(d/2\) as well as the bit of value two, and gains one on pairing. High blocks have an extra cost of at least one. When \(2\le v_2(k)<e\), the earlier paired bound \(2-2^{v_2(k)}r+v_2(k)+1+\tau\) suffices. For \(d\mid k\), the raw bound \(1+e+\tau\) suffices.

## The strength bounds

For \(h\ge2/d\), the global minimum of \(\beta\) is at least \(2-2r-e\). On \(m\ge d\) it is at least \(3-2r-e\): if \(v_2(m+1)<e\) use \(m\ge d\), if \(v_2(m+1)=e\) use \(m\ge3d-1\), and for higher denominator valuations compare the dyadic endpoints starting at \(2d-1\). The radial bounds of the final-even proof give minimum at least \(5-3r-e+\tau\) for the odd-power groups. At even powers the nonconstant terms are positive directly. For odd \(m<d\) this uses paired gain one, and for odd \(m\ge d\) the larger bound on \(\beta\). The even terms use their multinomial gain \(v_2(k)\), with paired gain one when \(m<d,d\nmid k\). At \(d\mid k\), that multinomial gain is at least \(e\). Thus \(\tau\ge e-4\) gives a unit \(M_N\).

Now take \(2/d<h<4/d\); only \(e-\tau\ge5\) needs consideration. At \(m=d-1\) the leading odd pair \(k=d-1,2d-1\) gains exactly one for the pure linear term, because its radial ratio is within valuation \(>1\) of \(-1\) and its binomial ratio is three modulo four. Its valuation is
\[
                            3-r+\tau+(d-1)h-e.               \tag{3}
\]
Every other term with this \(m\) is higher, by the positive curvature penalty or the nonminimal degree bound. For \(m<d/2\), the minimum \(\beta\) is at \(m=d/2-1\). If \(h<3/d\), the original pair gains put that range at least \(3/2-dh/2>0\) above (3). If \(h\ge3/d\), the new pairing makes the margin \(2-dh/2>0\). The even degrees lie higher by the corresponding even paired bounds. For \(d/2\le m<d\) except \(d-1\), the denominator valuation is at most \(e-2\) and the least \(\beta\) is bounded below by \(\beta(3d/4-1)\), since \(m+1\ge d/2+2^{v_2(m+1)}\) and \(h<4/d\). Its margin over \(\beta(d-1)\) is \(2-dh/4>0\).

For \(m\ge d\), the candidate minima of \(\beta\) are bounded below by the values at \(3d/2-1\), \(3d-1\), and \(2d-1\), for denominator valuations below, equal to, or above \(e\), respectively. Relative to the paired value (3), their raw odd bounds have margins \(dh/2,\ 2dh-1,\ dh-2\), all positive. The even terms have a higher base valuation. All other groups either lie strictly above (3) or have positive valuation, so the only possible equality with the constant is \(h=(e-\tau-3+r)/(d-1)\). Comparing it with \(2/d\) and \(4/d\) gives exactly \(e-\tau=5,6\).

For the last stated sufficient condition set
\[
               \beta_*(h)=\min_{s\ge0}\bigl[(2^s-1)h-s\bigr].
\]
The small-\(m\) groups are positive if \(\beta_*(h)>-(c-r+\tau)\), with \(c=7/2\), or \(c=4\) when \(3/d\le h\le1/2\) by the second pairing. The exact strength threshold for this strict inequality is \((5-c+r)/(2^{\tau+5}-1)\). Indeed the successive ratios \((s-c+r-\tau)/(2^s-1)\) increase through \(s=\tau+5\) and decrease thereafter wherever their numerators are nonnegative. This follows by cross multiplication of consecutive ratios. When \(h>1/2\), all \(\beta(m)\ge-1/2\), so even the raw bound \(\eta_N(k)\ge1\) suffices.

For \(d/2\le m<d\) besides \(d-1\), the bound \(\beta(m)\ge (d/2)h+\beta_*(h)\) automatically suffices. At \(m=d-1\), positivity of (3) is exactly the second strength threshold, and the other degrees lie higher. For \(m\ge d\) with denominator valuation below \(e\), use \(\beta(m)\ge dh+\beta_*(h)\). At \(e\) or above it is enough to test the bound at \(2d-1\), since the case of equality with \(e\) has \(m\ge3d-1\) and higher \(\beta\). This last threshold is \((e-\tau-1+r)/(2d-1)\), strictly smaller than the second stated threshold when \(e-\tau\ge5\). This proves the sufficient conditions. In particular, for \(e-\tau=5,6\) they cover \(h\ge4/d\).

## Exact checks

The [finite replay](verify_binary_even_quadratic_close_axis_second_pairing.py) checks 744 coefficient crossings, 1,488 axis gaps, 45 radial ratios, 2,695 paired counts, 127 strength thresholds and 24,363 frontier comparisons.
