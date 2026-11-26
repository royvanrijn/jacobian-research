# Universal factorial marks from the coefficient envelope

Consider an original ordinary two-monomial profile. Remove equal rows and write \(a_i\) for the smaller exponent of row \(i\), \(d_i>0\) for its gap, and \(\sigma_i=+1\) or \(-1\) for its orientation. Write
\[
 A=\max_i a_i,\quad B=A+1,\quad
 K_{j,i}(k)=\begin{cases}k,&\sigma_i=+1,\\ j-k,&\sigma_i=-1,\end{cases}
 \qquad E_j(k)=\max_i(ja_i+d_iK_{j,i}(k)).
\]
This is the upper envelope of the factorial arguments at order \(j\). The results below apply to the unchanged integer profile, with no bounds on the heights, gaps or multiplicities.

## An exact mark window in every profile

For \(j\ge3\) put
\[
 u_j=\min_{0\le k\le j}E_j(k),\qquad
 T_j=\min_{0\le k\le j-2}\max(E_j(k),E_j(k+2)),\qquad
 G_j=T_j-u_j .
\]
For any integer mark \(x\in(u_j,T_j]\), let \(I(x)=\{k:E_j(k)<x\}\). This is a nonempty interval of integer indices of width at most one. At any index in \(I(x)\) the literal factor \(x\) occurs in no row, while at every index outside it at least one row contributes \(x\). Indeed \(u_j\ge jA\), so \(x\) is above every smaller-end argument \(ja_i\). After stripping \((ja_i)!\) from the row, the literal factor is present exactly when \(x\le ja_i+d_iK_{j,i}(k)\). Convexity of \(E_j\) makes \(I(x)\) an interval. It contains a minimizing index and cannot contain \(k,k+2\) by the definition of \(T_j\). Conversely, any integer \(x>T_j\) puts three consecutive indices in \(I(x)\), since the middle ordinate is bounded by the maximum of its neighbors.

One always has
\[
                  \min_i d_i\ \le\ G_j\ \le\ 2\max_i d_i.       \tag{1}
\]
In every triple, choose a row attaining the middle envelope ordinate. At the endpoint in the direction of its slope it is higher by at least \(\min_i d_i\). This proves the lower bound. For the upper bound, choose a triple containing a minimizing index. The envelope is Lipschitz with constant \(\max_i d_i\), so both endpoints of that triple have ordinate at most \(u_j+2\max_i d_i\).

We require \(j\nmid x\) when comparing orders \(j,j'\) so that no literal factor at order \(j'\) is a proportional copy of \(x\). The window fails to contain a permissible mark exactly when \(G_j=1\) and \(u_j+1\) is divisible by \(j\).

**Corner classification.** An order \(j\ge3\) lacks a permissible mark if and only if, after a possible common reversal of the columns, the following three conditions hold.

- There is an increasing unit-gap row of smaller exponent \(A\).
- All increasing rows have larger exponent at most \(B\).
- The decreasing rows satisfy
  \[
             \max_{\sigma_i=-1}\bigl(2d_i-j(B-a_i)\bigr)
                             \in\{-1,0\}.                    \tag{2}
  \]
  In particular a decreasing row exists.

To prove necessity, \(x=u_j+2\) gives three consecutive indices where the envelope is at most \(u_j+1\). At the middle index a maximizing row must have gap one and ordinate \(u_j\), since its increasing endpoint inside the triple would otherwise exceed \(u_j+1\). Because \(u_j\equiv-1\bmod j\), this middle index is \(j-1\) if the row increases, or \(1\) if it decreases. Reverse in the latter case. The unit row has smaller exponent \((u_j+1)/j-1\), which equals \(A\) because \(jA\le u_j\). The endpoint \(j\) has \(E_j(j)\le jB\), so every increasing larger exponent is at most \(B\). At index \(j-2\) each increasing row has ordinate at most \(u_j-1\), by comparing with \(j-1\). Thus a decreasing row attains there the envelope, which lies in \(\{u_j,u_j+1\}\), proving (2).

Conversely, the unit row gives ordinate \(jB-1\) at \(j-1\). Every increasing row has ordinate at most \(jB-1\) there by the bound on its larger exponent. Equation (2) bounds the envelope of the decreasing rows at \(j-2\) between \(jB-1\) and \(jB\), and at \(j-1\) it is at most \(jB-1\). For earlier indices the decreasing envelope is at least \(jB-1\). At index \(j\) the increasing unit row gives \(jB\), and no row exceeds that value at \(j-2,j-1,j\). Hence \(u_j=jB-1,\ G_j=1\), and the unique mark is divisible by \(j\).

**Consecutive classification.** For consecutive orders \(n,N=n+1\), \(n\ge3\), at least one order has a permissible mark unless the following holds after a common reversal of the columns:
\[
 \begin{gathered}
      n\ \text{is even},\qquad
      \text{there is an increasing row }(A,A+1),\\
      a_i+d_i\le B\quad(\sigma_i=+1),\qquad
      2d_i\le n(B-a_i)\quad(\sigma_i=-1),\\
      \text{there is a decreasing row }(A+n/2,A).
 \end{gathered}                                               \tag{3}
\]
For simultaneous corners the distinguished directions must agree. If they differed, the upper-exponent bound would apply to all rows. Then every expression in (2) at order \(N\ge4\) would be at most \((2-N)(B-a_i)\le-2\), a contradiction. In a common direction the maximum in (2) decreases by at least one per order, since each \(B-a_i\ge1\). To be in \(\{-1,0\}\) at both orders it must be \(0,-1\), respectively. A row attaining \(-1\) at \(N\) must have \(B-a_i=1\) and \(2d_i=n\), or its expression at \(n\) would be positive. The upper bound \(2d_i\le n(B-a_i)\) comes from the maximum at \(n\). Conversely (3) gives maxima \(0,-1\), so both windows fail by the first classification. In particular every even terminal \(N\ge4\) has a permissible mark in at least one parent order.

## An arithmetic certificate from the actual ramps

Fix a permissible mark at order \(m\) of a pair, and write \(n'\) for the other order. Let \(I=I(x)\), \(r_k\) be the number of copies of the literal \(x\) contributed at index \(k\), \(r_*=\min_{k\notin I}r_k\), and \(w=\max_{0\le k\le m}\operatorname{dist}(k,I)\). The set of outside indices is nonempty since \(m\ge3\), and \(r_*\ge1\).

For each row \(i\), let \(L_i(x)\) be the least common multiple of the factors of its two full stripped ramps, deleting the literal \(x\) from the marked-order ramp only. Take the least common multiple of an empty set to be one. A physical common root forces
\[
               x^{r_*}\ \mid\ m!(n'!)^{w}
                                         \prod_i L_i(x)^{\,w+2}.   \tag{4}
\]
There is a more selective prime test. At any prime \(p\mid x\) put \(e=v_p(x)\), \(H_i=v_p(L_i(x))\), \(R=\sum_i H_i\), and \(b_j=v_p(j!)\). If for every \(k\notin I\),
\[
 \sum_{i:\ x\text{ occurs at }k}(e-H_i)_+
           > \bigl(1+\operatorname{dist}(k,I)\bigr)R+
                         b_m+\operatorname{dist}(k,I)b_{n'},      \tag{5}
\]
there is no common root.

In fact, truncate the valuation layers of row \(i\) at \(H_i\) in both orders. The factorial floor formula shows that subtracting the same slope
\(\sigma_i d_i\sum_{t=1}^{H_i}p^{-t}\) times the coefficient index leaves an oscillation of at most \(H_i\) in each order. Each layer contributes, up to a constant, a fractional part with range at most one. Every truncated factor other than the marked literal was included to its full valuation, by definition of \(H_i\). Thus the only missing layers are the contributions on the left of (5). In the unmarked order, after common slope removal and with the binomial coefficient included, the ordinates range over an interval of width at most \(R+b_{n'}\). The lower polygon stays in that interval, so its slopes have absolute value at most \(R+b_{n'}\).

Fix one of those slopes and subtract it times the index at order \(m\). The marked ordinate at any \(k\notin I\) exceeds the ordinate at its nearest \(q\in I\) if (5) holds, by the marked truncated oscillation bound \(R+b_m\). A marked edge of shared slope must therefore have both endpoints in \(I\). The total shared width is at most one. The Newton polygon product rule bounds the degree of a common factor by this width. The even moment is strictly positive on the real line, since it is the expectation of an even power of a nonzero real binomial under independent exponentials. A positive-degree rational common factor therefore has degree at least two. This proves the prime test. For (4), the left of (5) is at least \(r_*e-R\). If the divisibility (4) fails, its failing prime satisfies \(r_* e>(w+2)R+b_m+w b_{n'}\), which suffices. [Filaseta and Rocks, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

The caps in (5) only charge primes of the literal mark that also occur off the mark. Where desired, the caps can be bounded by cross-differences. For any \(y\) in the marked order distinct from \(x\), use \(y-x\). For any \(y\) in the other order use \(my-n'x\), which is nonzero since \(m\) and \(n'\) are coprime and \(m\nmid x\). If \(e\) exceeds the valuation of every cross-difference for a row, then all the actual off-mark valuations for that row are bounded by those cross-difference valuations. The use of actual caps in (4) and (5) requires no depth assumption.

There is also a factorization-free test using the exact coefficient ratios. Let \(c_{j,k}\) be the stripped integer coefficient, including the binomial factor. For \(k\notin I\), choose its nearest \(q\in I\) and set \(d=|k-q|\). For each \(1\le\ell\le n'\) form the rational number
\[
 \left(\frac{c_{m,k}}{c_{m,q}}\right)^\ell
 \begin{cases}
 \left(c_{n',n'-\ell}/c_{n',n'}\right)^d,&k>q,\\
 \left(c_{n',\ell}/c_{n',0}\right)^d,&k<q .
 \end{cases}                                                  \tag{6}
\]
If the positive reduced numerators of all these numbers have a common prime divisor, there is no common root. At that prime, positivity of all the corresponding valuations says exactly that, after subtracting any unmarked slope times the index, the ordinate at \(k\) exceeds the one at \(q\). The largest unmarked slope is the maximum of all suffix secants from its right endpoint, and the smallest is the minimum of all prefix secants from its left endpoint. These statements follow by comparing each secant with the lower polygon and using the endpoints of its last or first edge for equality. This proves (6) by the same width argument. Unlike (5), the test takes account of cancellation at every coefficient index.

These criteria extend interior mark tests that require preselected rows and first-block tests anchored at an endpoint. The common-slope step is the same factorial floor comparison used there. [Deep blocks](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md) [First-block depths](MULTIVARIATE_FIRST_BLOCK_DEPTH_AND_MOVING_RESIDUES.md) [Oriented secants](MULTIVARIATE_ORIENTED_NATIVE_TWO_SLOPE_BUDGET.md) [Common-slope removal](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md)

## A compulsory large-prime load when the window grows

Fix \(J\ge4\) and a pair with \(n,N<J\). Assume total curvature is below twelve and that the window width \(G=G_m\) of one order tends to infinity along a family. Take all marks \(x\equiv1\bmod m\) in its window, and write \(P\) for their product and \(Q\) for their count. The marks have \(Q\sim G/m\), and
\[
                  \log P_{>G}\ge (1/2-o(1))Q\log A,             \tag{7}
\]
where the subscript denotes the prime part above \(G\). Indeed curvature below twelve gives \(d_i^2<12(a_i+d_i)\), hence \(\max_i d_i=O(\sqrt{A+1})\). By (1), \(G=O(\sqrt{A+1})\); all the marks have size between \(A\) and \(O_J(A+1)\). At each prime \(p\le G\) dividing the product, the step \(m\) is a unit, since every mark is \(1\bmod m\). Remove one mark having the largest valuation at that prime. At every prime power \(p^v\) the remaining indispensable positions are spaced by \(p^v\), and their count is at most \(\lfloor(Q-1)/p^v\rfloor\). Thus the remaining small-prime part divides \((Q-1)!\). The removed factors cost at most \(\pi(G)\log(\max x)=o(Q\log A)\) by the elementary bound \(\pi(G)=O(G/\log G)\). Finally \(\log((Q-1)!)\le(1/2+o(1))Q\log A\). [Elementary prime count](MULTIVARIATE_DENSE_CLUSTERED_UNIT_GAP_CLOSURE.md#elementary-prime-estimates)

A prime \(p>G\) dividing \(P\) divides a unique mark \(x_p\) because the mark span is less than \(G\). Let \(r_{\min}\) be the smallest \(r_*(x)\) among the chosen marks. For large \(G\), the prime \(p>G\) is greater than both orders. The necessary inequality (4) at \(p\), summed over these primes, gives
\[
 \sum_{\substack{p>G\\p\mid P}}\sum_i v_p(L_i(x_p))\log p
                 \ \ge\ \frac{r_{\min}}{m+2}\log P_{>G}.        \tag{8}
\]
This is an exact budget on the original rows at the primes of the marks, with no upper cluster assumption. The growing-window inequality gives a uniform lower bound on that budget. A matching upper bound is still needed for arbitrary profiles.

## The odd-terminal limitation is achieved by ordinary rows

Let the lower order be \(m=2q\) with \(q\ge2\). At one smaller exponent \(a\), take equally many increasing gap-one rows and decreasing gap-\(q\) rows. At the smaller exponent \(a+q+1\) take equally many copies of the reversed pair. The two groups' full ramps are disjoint at either order, because the upper endpoint of the first group is below the start of the second.

At either order \(j=m,m+1\) the literal count in one group, divided by its common multiplicity, has the form
\[
                  {\bf 1}_{k\ge h}+{\bf 1}_{k\le j-\lceil h/q\rceil}
\]
or its reversal, where \(h\) is the offset from the smaller end \(ja\) of that group. If \(1\le h\le j-1\), the two indicated index sets cover the whole range and both endpoints are minimizing. If \(h>j\), the unit row is absent and the gap-\(q\) row leaves at least three consecutive minimizing indices. The offset \(h=j\) is divisible by \(j\) and fails the nonmultiple requirement. Thus at this consecutive pair **every** nonmultiple literal has minimum width at least two, even if the mark need not come from the envelope window.

The signed degree difference of the full profile is zero. With the same multiplicity \(\lfloor\tau a/(2(1+q^2))\rfloor\) for each type, the largest relative gap tends to zero and the curvature tends to any prescribed \(\tau>0\) as \(a\to\infty\). The example is an obstruction to a single isolated nonmultiple-literal argument. It makes no claim about a common root of the ordinary moments.

## Exact checks

For the two mirrored rows \((1,2),(2,1)\) at orders five and six, the order-five mark is \(x=9\) and \(I=\{2,3\}\). At \(p=3\) the two valuation arrays are
\[
         (3,4,2,2,4,3),\qquad (3,3,3,4,3,3,3).
\]
The first polygon's slope zero has width one, while the second is horizontal. The numerator gcd in (6) is \(3\). Both omitted-literal ramp least common multiples equal \(27720\), so the divisibility test (4) passes and gives no exclusion. This verifies that the exact secant test can certify a pair that satisfies the coarser necessary divisibility.

The [finite checker](verify_multivariate_universal_envelope_marks.py) verifies 20,468 individual windows and 17,544 consecutive classifications, including 204 simultaneous corners. It replays 998 secant certificates, 83 layer and divisibility certificates, 415 progression-prime bounds, and 3,374 literal failures in the sharp odd-terminal family. These are finite arithmetic checks. The unrestricted window classification and necessary prime-load inequalities are proved above.
