# Shared terminal-eleven separation with outer mean at least .95

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be the normalized order-ten outer geometric mean after
removing the binomial factors. If \(u\ge19/20\), the tenth and eleventh
factorial moment polynomials have disjoint root arguments, including
after independent positive rescalings. Rational Bézout gives coprimality
over every characteristic-zero field.

Write \(p=1-u\), and let \(q=1-w\) be the order-eleven deficit. For
\(p>0\) put \(r=q/p\). The match from an order-ten reciprocal
\(x=d/(10\ell+h)\) to order eleven is \(x'=10x/(11-hx/d)\). Since
\(u^{-2}\) is the product of the factors \(1+9x\), we have
\(x\le c=(u^{-2}-1)/9\). Thus
\[
                 \frac{100}{99}\le \frac{10x'}{9x}
                          \le \frac{100}{9(11-c)}=b_0.
\]
Writing \(A_i=1+9x_i\), the lower match bound gives
\(1+10x_i'\ge A_i(1+(1-A_i^{-1})/99)\).
Multiplying and using \(\sum_i(1-A_i^{-1})\ge1-u^2\) gives
\(w\le u(1+(1-u^2)/99)^{-1/2}\).
The upper bound and \((1+s)^{b_0}\ge1+b_0 s\) for \(s\ge0\)
give \(w\ge u^{b_0}\). Integrating these two powers yields
\[
 1+\frac{u(1+u)}{198(1+(1-u^2)/99)^{3/2}}
              \le r\le\frac{100}{9(11-c)} .
\]
The checker uses \(r\in[10093/10000,5058/5000]\),
\(\sqrt p\le 2236068/10^7\), and \(q\le q_0=51/1000\).
The lower bound increases and the upper bound decreases as \(u\) increases.
These bounds also hold under common positive contraction of the
reciprocals at the two orders, which keeps the individual ratios
\(10x'/9x\) fixed.

## Grouped coefficient enclosure and exact sectors

At either order \(j\), with its own deficit \(h=p\) or \(q\), use the
[stable Cayley transform and grouped Taylor enclosure](MULTIVARIATE_TERMINAL_TEN_NEAR_FLAT.md#a-grouped-taylor-enclosure-through-deficit-weight-seven).
The same weight-seven expansion keeps the signed even and unsigned
odd log-moment directions and every exponential monomial of weight
at most seven. The log moment and raw error bounds use \(C=1/(1-q_0)\).
The checker uses the exact rational Cayley transform of each retained
direction at orders ten and eleven. It verifies cancellation before
division by \(\lambda_j^j\), with \(\lambda_j=\sqrt h\).
The independent raw remainder starts at \(h^8\), so its transformed
powers \(\lambda_j^{16+a-j}\) are positive at every degree
\(0\le a\le j\).

Leading coefficient dominance gives \(|Z|<2\) in the Cayley scale,
which maps into the physical left half-plane since \(4q_0<1\).
Nine separating rays use levels \(b\) that vary piecewise affinely
in \(p\). The level in the eleventh-order scale is \(b/\sqrt r\).
The checker uses the height equation and centered interval
projections from the same stable Cayley proof. It verifies strict
point exclusion throughout every ray within the root disk.

At the zero-deficit anchor exact signs between consecutive rays give
upper root orders \(11,10,11,10,11,10,11,10,11,10\).
The common contraction of the reciprocals preserves these sector
counts away from the degenerate endpoint by continuity of the
Cayley polynomials. The odd-order root in the joined sector through
the negative real axis stays real by conjugacy. The even-order
factorial moment is strictly positive on the real axis by its
Gamma integral, so its two roots there are nonreal. Thus the
arguments are disjoint in that sector as well.

The [exact replay](verify_multivariate_terminal_eleven_stable_upper.py)
certifies 143,580 projection cells in 180 order, ray, and deficit
subband checks. It compares each count with the frozen receipt.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_stable_upper.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_stable_upper.py 11 -1 25 50
~~~
