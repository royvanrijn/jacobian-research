# Close-axis exclusion at reciprocal dyadic endpoints

Normalize the constant coefficient to one in \(f=ax^2+bxy+cy^2+1\). Put \(M_j=\Lambda_2(f^j)\), and let \(v=v_2\) with \(v(2)=1\). Suppose \(a\ne c\), \(d=2^e\ge4\), and \(v(4a)=v(4c)=1/d\). Choose \(t^d=2\) in a valued extension, and put
\[
 X=4a/t,\quad Y=4c/t,\quad
 \gamma=(c-a)/a,\quad \epsilon=(b-a-c)/a,\quad h=v(\gamma).
\]

**Theorem.** Assume \(X,Y\) both reduce to one, and \(2d\mid N\). If \(M_{N-3},\ldots,M_N\) all vanish at positive indices, then \(h\ge 2/d\). At \(h=2/d\), the residue of \(\gamma/t^2\) must be one. The conclusions hold at every \(d=2^e\ge4\).

## Boundary equations and the curvature valuation

Write \(g_j(u)=\Lambda_1((1+ux^2)^j)\), \(J=N-1\), and \(D=2d\). The exact weighted boundary equations give
\[
 g_J(a)g_N(c)=g_J(c)g_N(a),\qquad
 b-a-c=(a-c)\frac{g_N(a)-g_N(c)}{g_N(a)+g_N(c)} .                  \tag{1}
\]
For the determinant the two weights do not both vanish since \(a\ne c\). For the second identity use the boundary equation at \(N\). In the present hypotheses the denominator has valuation one.

Indeed the [even-wrap axis calculation](BINARY_EVEN_QUADRATIC_UNIFORM_ENDPOINTS.md#the-case-nd-even) shows \(g_N(tZ/4)\) has unit constant, every nonconstant term has valuation at least one, and only the even indices \(d,D\) can have valuation one. Since \(X,Y\) have the same residue and all nonconstant coefficients are divisible by two, their normalized sum satisfies \((g_N(a)+g_N(c))/2\equiv1\bmod\mathfrak m\), where \(\mathfrak m\) is the maximal ideal. Also \(v(X-Y)=h>0\). At even positive indices \(i\), the factorization of \(X^i-Y^i\) as a difference of squares gives \(v(X^i-Y^i)>h\). Thus \(v(g_N(a)-g_N(c))>1+h\), and
\[
                                 v(\epsilon)>2h.                  \tag{2}
\]

## Pairing the two low axis degrees

Use \(s\) for binary digit sum. The valuation of the \(k\)-th coefficient in \(g_J(tZ/4)\) is
\[
                         \eta(k)=k/d+s(J-k)-s(J).
\]
Write \(J=DH+D-1\), \(k=Dq+w\), \(0\le w<D\). Then
\[
                  \eta(k)=2q+s(H-q)-s(H)+w/d-s(w).                \tag{3}
\]
The high part is at least \(2q-s(q)\), positive if \(q>0\). The low part minimizes exactly at \(w=d-1,D-1\), since all bits of size less than \(d\) have negative cost and the bit of size \(d\) has zero cost. Denote the overall minimum by \(\eta_*\). Within each block, \(k\) and \(k+d\), for \(w<d\), have equal valuation.

Divide \(g_J(tZ/4)\) by its coefficient at \(Z^{d-1}\), calling the resulting polynomial \(G(Z)\). Let \(u_*\) be the ratio of its coefficients at \(Z^{D-1}\) and \(Z^{d-1}\). For every low pair the ratio of the coefficient at \(Z^{k+d}\) to that at \(Z^k\) is an odd rational unit, since \(t^d=2\) and the two coefficients have the same valuation. Hence there is an exact decomposition
\[
                G(Z)=(1+u_* Z^d)P(Z)+E(Z),                        \tag{4}
\]
where \(P\) takes the lower term of each pair and reduces to \(Z^{d-1}\). Every coefficient of \(E\) has valuation strictly greater than one. At the leading pair there is no error, while at every other pair there is a positive valuation gap and the two odd rational ratios agree modulo two.

Take divided differences between \(X,Y\) in (1). Since \(G(Y)\) has positive valuation and the coefficients of \(g_N(tZ/4)-1\) have valuation at least one, the divided difference of \(G\) has valuation strictly greater than one. The divided difference of \(P\) is a unit because its reduction is \(Z^{d-1}\) and both residues are one. From (4),
\[
 v(1+u_*X^d)\ \ge\
          \min\left(\rho,\ v\!\left(\frac{X^d-Y^d}{X-Y}\right)\right)
 \quad\text{for some }\rho>1 .                                   \tag{5}
\]
For \(h<2/d\), factoring the difference into successive square differences gives \(v((X^d-Y^d)/(X-Y))=(d-1)h\). At each \(0\le i<e\), \(v(X^{2^i}-Y^{2^i})=2^ih<1\), so the corresponding sum has this same valuation. At \(h=2/d\) only the last sum might have greater valuation, giving the lower bound \(2-2/d>1\) in (5).

## Original radial moments

The change of variables \(s=x+y,\ z=y/(x+y)\) in the two exponential integrals, applied coefficientwise, gives
\[
 M_J=\sum_{k=0}^J C_k\int_0^1[1+(\gamma+\epsilon)z-\epsilon z^2]^k\,dz,
 \qquad C_k=\binom Jk(2k+1)!a^k.                                  \tag{6}
\]
This is a polynomial identity over characteristic zero. The valuation of \(C_k\) is \(\eta(k)\). A term contributing the power \(z^m\) from the integrand has valuation, including integration,
\[
                    \ge \eta(k)+\beta(m),\qquad
                    \beta(m)=mh-v_2(m+1),\quad 0\le m\le2k.       \tag{7}
\]
The bound is strict when it uses a quadratic factor, by (2). For the pure linear term, equality holds when the binomial coefficient \(\binom km\) is odd.

For \(m<d\), pair the contributions from \(k,k+d\) for every block in (3). For each pair and each fixed number of quadratic picks, the two multinomial coefficients have the same parity, with a missing term counted as zero. Indeed the generating polynomials for the two coefficients differ modulo two by multiplication by \(1+u^d+v^d\) in the two count variables, and their relevant total count is at most \(m<d\). The ratio \(C_{k+d}/C_k\) is \(X^d\) times an odd rational unit, so summing each pair gains at least
\[
                        \min(1,v(1+u_*X^d))                      \tag{8}
\]
above (7). This also works if both multinomial coefficients are even. For integer \(k\) in a high block \(q\ge1\), the increase in \(\eta\) is at least \(2q-s(q)\ge1+\lfloor\log_2 q\rfloor\). A term with \(m\ge D-1\) can gain at most \(1+\lfloor\log_2 q\rfloor\) in its denominator relative to \(m=D-1\): indeed \(m+1\le2k+1<D(2q+2)\). If \(m>D-1\) it still costs strictly more than \(\eta_*+\beta(D-1)\) by the larger linear exponent, and if \(m=D-1\) the upper block cost is strictly positive.

For a fixed \(v_2(m+1)=s_0\) the smallest \(\beta(m)\) occurs at \(m=2^{s_0}-1\). Along these indices the successive increments are \(2^{s_0}h-1\). When \(h<1/d\), among \(m\le2D-2\) the unique minimum is at \(m=D-1\). All larger-degree blocks are excluded by the high-block bounds just proved, using this strict inequality also for \(m<D-1\). Only \(k=m=D-1\) attains \(\eta_*+\beta(D-1)\), since a pure linear term needs \(m\le k\), and (2) discards quadratic picks there.

For \(1/d\le h<2/d\), \(\beta(d-1)\) is the smallest valuation for \(m<d\); the target \(\beta(D-1)\) lies \(dh-1\) higher. By (5), the gain (8) is strictly larger than \(dh-1\). Thus all \(m<d\) are suppressed in the paired sums. For \(d\le m<D-1\) one has \(v_2(m+1)\le e-1\), and therefore
\[
                    \beta(m)-\beta(D-1)
                        \ge 2-(D-1-m)h>0 .
\]
The preceding high-block bounds cover \(m\ge D-1\), and in the lowest block \(m>D-1\) has larger cost. Again only \(k=m=D-1\) remains. It has nonzero residue after normalization, contradicting \(M_J=0\).

Finally let \(h=2/d\). Now the minimum for \(m<d\) is reached exactly at \(m=d/2-1,d-1\), and the target \(\beta(D-1)\) is one higher. By (5), \(v(1+u_*X^d)>1\), so the gain in (8) is at least one. All pairs except the leading one \(k=d-1,D-1\) are strictly higher by (3). Curvature terms are also strictly higher. For the leading pair, put \(R=C_{D-1}/C_{d-1}\). This is \(u_*X^d\) times \((4d-1)/(2d-1)\), which is one modulo four. Thus \(v(1+R)>1\). The binomial ratio
\[
 \frac{\binom{2d-1}{m}}{\binom{d-1}{m}}
    =\prod_{i=1}^m(1+d/(d-i))
\]
is one modulo four for \(m=d/2-1\) and three modulo four for \(m=d-1\), because only \(i=d/2\) contributes a nontrivial factor modulo four. The pair at \(d/2-1\) therefore lies above the target, and the pair at \(d-1\) gains exactly one. The other least term is \(k=m=D-1\). Its residue relative to the surviving pair is the residue of \(R(\gamma+\epsilon)^d/4\), namely \(\overline{\gamma/t^2}^{\,d}\). Cancellation forces this residue to be one.

## Exact checks

The [finite replay](verify_binary_even_quadratic_close_axis_endpoint.py) checks 243 radial identities, 3,720 block bounds, 380,800 frontier terms, and 20,088 parity and binomial-ratio identities. The proof supplies the uniform exclusions.
