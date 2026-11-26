# Final-even nonvanishing for close axes at dyadic endpoints

Normalize the constant to one in \(f=ax^2+bxy+cy^2+1\), and put \(M_j=\Lambda_2(f^j)\). Let \(v=v_2\) with \(v(2)=1\). Assume \(a\ne c\), \(d=2^e\ge4\), and \(v(4a)=v(4c)=1/d\). Take \(t^d=2\) and assume both \(X=4a/t\) and \(Y=4c/t\) have residue one. Write
\[
 D=2d,\quad T=N/D,\quad h=v((c-a)/a),\quad J=N-1 .
\]
We consider \(D\mid N\), with \(N-3\ge1\).

**Theorem.** Under these assumptions, the four moments ending at \(N\) cannot all vanish in any of the following ranges.

- For \(d=4,8,16\), there is no restriction on \(h\).
- For \(d=32\), all \(h\) are covered if \(T\) is even. If \(T\) is odd, only the two strengths \(h=2/d\) and \(h=(2+1/d)/31\) remain.
- For \(d\ge64\), a four-zero block requires \(2/d\le h\le (3/2+1/d)/31\).

For \(h<2/d\), these statements use the [earlier close-axis exclusion](BINARY_EVEN_QUADRATIC_CLOSE_AXIS_ENDPOINT.md). Below we handle \(h\ge2/d\) using the original final even moment \(M_N\).

## The boundary quotient at this strength

Put \(\gamma=(c-a)/a\) and \(\epsilon=(b-a-c)/a\). Assume the four moments vanish. The [weighted boundary and curvature calculation](BINARY_EVEN_QUADRATIC_CLOSE_AXIS_ENDPOINT.md#boundary-equations-and-the-curvature-valuation) gives \(v(\epsilon)>2h\) and \(g_J(a)g_N(c)=g_J(c)g_N(a)\), where \(g_j(v)=\Lambda_1((1+vx^2)^j)\).

Let \(G(Z)\) be \(g_J(tZ/4)\) divided by its coefficient of \(Z^{d-1}\), and let \(u_*\) be the ratio of its coefficients at \(2d-1\) and \(d-1\). The [axis pairing](BINARY_EVEN_QUADRATIC_CLOSE_AXIS_ENDPOINT.md#pairing-the-two-low-axis-degrees) writes
\[
                  G(Z)=(1+u_*Z^d)P(Z)+E(Z).
\]
Here \(u_*\) is an odd rational unit and \(P\) reduces to \(Z^{d-1}\). The coefficients of \(E\) have valuation at least \(3/2\). In fact every nonleading pair has valuation gap at least \(1/2\) by (3) in the cited note, since omitting any profitable low bit costs at least \(1/2\), and the error gains one from subtracting odd rational units.

For \(h\ge2/d\), the squaring factors give
\[
       v\!\left(\frac{X^d-Y^d}{X-Y}\right)
          \ge \sum_{i=0}^{e-1}\min(1,2^ih)\ \ge\ 2-2/d .
\]
At a unit argument \(g_N(tZ/4)=1+\) terms of valuation at least one. By the boundary equality the divided difference of \(G\) between \(X,Y\) has valuation at least \(1+v(G(Y))>1\). Since \(P[X,Y]\) is a unit, the displayed pairing first gives \(v(1+u_*X^d)>1\) and the same inequality at \(Y\). Thus \(v(G(Y))>1\). Using this back in the same divided difference and the coefficient bound on \(E\) yields
\[
                      v(1+u_*X^d)\ge3/2.                     \tag{1}
\]
In particular \(X^d\equiv1\bmod 2\), where this congruence means valuation of the difference at least one.

## Radial coefficients and their paired valuations

The radial formula at \(N\) is
\[
 M_N=\sum_{k=0}^N C_k \int_0^1[1+(\gamma+\epsilon)z-\epsilon z^2]^k\,dz,
 \qquad C_k=\binom Nk(2k+1)!a^k.                               \tag{2}
\]
The \(k=0\) term is one. Use \(s\) for binary digit sum and \(\eta_j(k)=k/d+s(j-k)-s(j)\). Thus \(v(C_k)=\eta_N(k)\), and for \(k\ge1\)
\[
   \eta_N(k)=v_2(N)+1/d-1+\eta_J(k-1),\qquad
   \eta_J(Dq+w)=2q+s(T-1-q)-s(T-1)+w/d-s(w).                 \tag{3}
\]
The first equality follows either from the digit sums or \(\binom Nk=(N/k)\binom{J}{k-1}\). The high-block cost in the second is at least \(2q-s(q)\). The low-bit costs are \(b/d-1\) for a present bit of value \(b\). Thus the minima of \(\eta_N\) are at \(k=d,D\), with valuation \(1+v_2(T)\). The odd \(k\) have minimum
\[
                            A=2-1/d+v_2(T)
\]
at \(k=d-1,D-1\), because \(k-1\) omits its bit of value one. Every other odd \(k\) costs at least \(A+1/2\), as it omits another profitable bit or uses a high block. If \(1\le\nu=v_2(k)<e\), the bit of value \(2^\nu\) in \(k-1\) is zero, giving
\[
                        \eta_N(k)\ge 2-2^\nu/d+v_2(T).        \tag{4}
\]

To produce \(z^m\) in (2), take \(i\) linear and \(\lambda\) quadratic factors, \(m=i+2\lambda\). Set
\[
                  \beta(m)=mh-v_2(m+1).
\]
The term valuation is at least \(\eta_N(k)+\beta(m)\) plus the valuation of its multinomial coefficient, strictly higher when \(\lambda>0\). If \(m\) is even then \(\beta(m)\ge0\), so all such terms at \(k\ge1\) have positive valuation. For odd \(m\), \(i\) is odd and the multinomial coefficient contributes at least \(v_2(k)\) by taking out the factor \(k/i\).

For \(m<d\), pair degrees \(k,k+d\) using the full blocks of \(k-1\) in (3). They have equal \(\eta_N\). If \(d\nmid k\), then \(\nu=v_2(k)=v_2(k+d)\); after division by \(2^\nu\) the two multinomial coefficients have the same parity. Indeed after taking out \(k/i\) they are coefficients of \((1+U+V)^{k-1}\) and \((1+U+V)^{k+d-1}\) in total degree less than \(d\), and \((1+U+V)^d\equiv1+U^d+V^d\bmod2\). Also \(C_{k+d}/C_k\) is \(X^d\) times an odd rational unit, so (1) gives a gain of at least one in the sum of each such pair, relative to \(\eta_N(k)+\nu+\beta(m)\).

For the minimal odd pair \(k=d-1,D-1\) and \(m<d/2\), the gain is at least \(3/2\). The multinomial coefficients now agree modulo four, since the below-\(d/2\) coefficients of \((1+U+V)^d-1\) are divisible by four. To check the ratio of the signed radial coefficients, write \(B_j(k)=\binom jk(2k)!\). The exact identity \(C_k=2Na(2k+1)(2k-1)B_J(k-1)a^{k-1}\) gives
\[
                  C_{D-1}/C_{d-1}\equiv u_* X^d\bmod4 .       \tag{5}
\]
Indeed \(B_J(k+1)/B_J(k)=2(J-k)(2k+1)\), and at \(k=d-2,D-2\) the two quotients of factors \(J-k\) and \(2k+1\) are one modulo four. The two radial factors \(2k+1,2k-1\) are also unchanged modulo four by increasing \(k\) by \(d\). Apply (1) to obtain the claimed pair gain.

## Exclusion by the constant term

For a fixed \(j=v_2(m+1)\), the smallest value of \(\beta\) occurs at \(m=2^j-1\). Along these indices the successive increments are \(2^j h-1\).

For \(d\le16\) and \(h\ge2/d\), the minimum over all \(m\) is at least \(2-2/d-e\ge-17/8\). All odd \(m<d\) give positive terms or paired sums by (3)–(4) and the gain of one. At degrees divisible by \(d\), the multinomial gain of at least \(e\) suffices without pairing. For \(m\ge d\), \(\beta(m)>-(2-1/d)\). In fact if \(j\le e-1\), then \(\beta(m)\ge3-e\); if \(j=e\), then \(m\ge3d-1\); and if \(j\ge e+1\) the dyadic endpoint bound is smallest at \(j=e+1\), of value \(3-2/d-e\). For odd \(m\ge d\), \(\eta_N(k)\ge2-1/d\) at odd \(k\); at even \(k\) the multinomial adds at least one to \(\eta_N(k)\ge1\). All terms are therefore positive, and \(M_N\) is a unit.

For \(d\ge64\), suppose \(h>H_d=(3/2+1/d)/31\). At \(H_d\) the smallest dyadic value of \(\beta\) is at \(m=31\), so for the given \(h\) all positive \(m\) satisfy \(\beta(m)>-(7/2-1/d)\). If \(m\ge d/2\ge32\) then \(\beta(m)>-(3-1/d)\): for \(j\le4\) use \(m\ge32\), for \(j=5\) use \(m\ge95\), and for \(j\ge6\) use the least such dyadic endpoint \(m=63\). If \(m\ge d\ge64\), the stronger bound \(\beta(m)>-(2-1/d)\) follows using \(m\ge64\) for \(j\le5\), \(m\ge191\) for \(j=6\), and the dyadic endpoint \(m=127\) for \(j\ge7\). Each comparison already holds at \(H_d\).

For \(m\ge d\) the last bound gives positivity as in the smaller \(d\) case. For odd \(m<d\), degrees divisible by \(d\) gain at least \(e\), and the other even pairs have lower bound \(2-2^\nu/d+\nu+1+\beta(m)>0\) by (4). The nonminimal odd pairs have base valuation at least \(A+1/2\) and gain at least one. The minimal odd pair gains \(3/2\) for \(m<d/2\), or gains at least one and uses \(\beta(m)>-(3-1/d)\) otherwise. Consequently \(M_N\) is a unit.

## The remaining middle size

Take \(d=32\), \(h\ge1/16\). The general dyadic bound is \(\beta(m)\ge-49/16\). Thus even degree terms or pairs with \(m<d\) are positive by (4). Nonminimal odd pairs gain at least one over \(A+1/2\), and the minimal odd pair at \(m<16\) gains \(3/2\). Terms with \(16\le m<32\) have \(\beta(m)\ge-2\) unless \(m=31\), because in \(16\le m\le30\) one has \(v_2(m+1)\le4\) and the value four does not occur. This leaves the \(m=31\) term in the minimal odd pair.

For \(m\ge32\), the only possible value \(\beta(m)\le- (2-1/d)\) is \(m=63\). For \(j\le3\) use \(m\ge32\); for \(j=4,5\) use \(m\ge47,95\), respectively. For \(j=6\) the endpoint is \(63\) and the next value is at least \(191\); for \(j\ge7\) the dyadic minimum is at \(127\). The \(m=63\) terms at even \(k\) are positive by (4), or by \(\eta_N(k)+v_2(k)\ge1+e\) at \(d\mid k\). At odd \(k\), the pure linear term at \(k=63\) attains the bound \(A+\beta(63)\) from (3). The other terms are strictly higher since \(k=31\) cannot supply \(m=63\) without curvature and curvature has a positive gain.

In the \(m=31\) term of the paired minimal odd degrees, \(\binom{63}{31}\equiv3\bmod4\). By (1) and (5) one has \(v(1+C_{63}/C_{31})>1\). Thus the pure linear pair gains exactly one, with valuation
\[
                         A+\beta(31)+1=v_2(T)+31h-2-1/d .
\]
At this \(m\), curvature pairs gain at least one and also the positive quadratic penalty, so lie strictly higher. The value of the pure \(m=63\) term exceeds the displayed value by \(32h-2\). If \(T\) is even both values are positive. Otherwise, for \(h>1/16\) the \(m=31\) pair is uniquely first among the possible nonpositive nonconstant groups, and \(M_N\) is nonzero unless that value is zero. Its zero is \(h=(2+1/d)/31\). The possible tie at \(h=1/16=2/d\) is the other stated exception.

## Exact checks

The [finite replay](verify_binary_even_quadratic_close_axis_final_even.py) checks 14,112 factorial shifts, 42 coefficient ratios, 38,099 paired counts and 8,152 dyadic inequalities. The proof applies to all orders.
