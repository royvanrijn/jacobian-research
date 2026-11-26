# Anchored derivative constraints for deep mixed even quadratics

Normalize the constant coefficient to one and write
\[
 f=ax^2+bxy+cy^2+1,\qquad M_j=\Lambda_2(f^j),\qquad v=v_2,\quad v(2)=1 .
\]
Suppose \(v(4a)\le-3\), \(v(c)>-1\), and \(N\ge4\). The assertions also hold with the two axes exchanged.

**Theorem.** If \(M_{N-3},\ldots,M_N\) all vanish, let
\[
 \tau=v(N),\quad L=2^\tau,\quad w=\tau-3,\quad W=wL,\quad
 Y=N/L-1 .
\]
Then \(\tau\ge4\) and \(Y\) is a positive even integer. Write \(B=\lfloor\log_2 Y\rfloor+1\) for its binary length, and set
\[
 V(q)=Wq+2s_2(q)-\operatorname{carr}_2(q,Y-q),\qquad 0\le q\le B.
\]
Here \(s_2\) is the binary digit sum, and \(\operatorname{carr}_2\) counts binary carries. For every set \(\mathcal A\subseteq\{0,\dots,B\}\) of at most two distinct anchors, the minimum
\[
              \min_{0\le q\le B}
                    \left[V(q)+\sum_{u\in\mathcal A}v_2(q-u)\right]       \tag{1}
\]
is attained at least twice. Use \(v_2(0)=+\infty\). This strengthens the [unanchored borrow constraint](BINARY_EVEN_QUADRATIC_DEEP_MIXED_BORROW_GATE.md), and uses three original moments together with the forced boundary identities. For instance, \(Y=2^{W+2}\) passes the unanchored test but fails (1) with anchor zero.

## The exact adjacent derivative identities

Write \(Z=4a\) and \(J=N-1\). By the [coupled two-adic calculation](BINARY_EVEN_QUADRATIC_COUPLED_TWO_ADIC_SECTORS.md#a-dominant-axis-and-approximation-by-its-single-face-plane), \(v(Z)=1-\tau\). Write \(M_j^0\) for the moment on the plane \(b=2a\) with the same \(a,c\), and normalize it by
\[
 P_j=\frac{M_j^0}{(j!)^2 Z^j},\qquad
 H_m=\sum_{i=0}^m \binom{2i}{i}(c/Z)^i .
\]
For \(j=J,J-1,J-2\), the perturbation bound for \(M_j-M_j^0\) and the assumed vanishing of \(M_j\) give
\[
 v(P_j)\ge (j-N)v(Z)+s_2(N)-s_2(j)+wj+2s_2(j)\ge w(J-2).             \tag{2}
\]
Define \(C(d)=(J-d)!/[J!\,d!\,Z^d]\) for \(0\le d\le J\). With \((t)_{\underline s}=t(t-1)\cdots(t-s+1)\) and \((t)_{\underline 0}=1\), the exact plane formula gives
\[
 T_s:=\frac{P_{J-s}}{(J)_{\underline s}Z^s}
    =\sum_{d=0}^{J}(d)_{\underline s} C(d)H_{J-d},\qquad s=0,1,2.       \tag{3}
\]
Indeed the term at reverse degree \(d-s\) in \(P_{J-s}\) has the same \(H_{J-d}\) as degree \(d\) in \(P_J\), and its factorial ratio is multiplied by \((J)_{\underline s}(d)_{\underline s} Z^s\). The denominator on the left has valuation at most zero for \(s=0,1,2\), since \(J\equiv -1\bmod L\), \(v(J-1)=1\) and \(v(Z)=1-\tau\). Thus every \(T_s\) has valuation at least \(w(J-2)\).

## Truncating the exact sums

For \(d=Lq+b\) with \(0\le b<L\), the [borrow calculation](BINARY_EVEN_QUADRATIC_DEEP_MIXED_BORROW_GATE.md#the-two-low-blocks-have-the-same-minimum) gives
\[
                 v(C(Lq+b))=V(q)+wb+2s_2(b).                         \tag{4}
\]
Here \(V\) can be evaluated for all \(q\le Y\), and \(H_m\equiv1\bmod 2\). The unanchored result provides some positive \(q\) with \(V(q)\le0\). This requires at least \(W+2\) carries, whereas a \(B\)-bit \(Y\) permits at most \(B-1\) carries. Therefore \(B\ge W+3\). In particular \(Y\ge 2B\) and \(W\ge16\).

For comparison use the integer precision bound
\[
                           R=2W+2+2\tau+B.                           \tag{5}
\]
Every \(q>B\) has \(V(q)\ge (W-1)B+W>R\), since the carry count is at most \(B\), and \(\tau\le W/4\). Such blocks can be removed from (3) to precision strictly above \(R\). On \(q\le B\), the least \(H\) index in the block is at least \(L(Y-B)\). Replacing every such \(H\) by the unit
\[
                            H_\infty=\sum_{i\ge0}\binom{2i}{i}(c/Z)^i
\]
has valuation error strictly above \(2L(Y-B)-B>R\) in (3), using \(v(c/Z)>\tau-2\ge2\) and \(V(q)\ge-B\). The sum can be taken in the completed valued field. The bound in (2) satisfies \(w(J-2)\ge WY>R\). These three strict comparisons also follow from \(B\ge W+3\), \(Y\ge2B\), and \(L\ge16\). Thus, after division by \(H_\infty\), the truncated versions of the three sums on the right of (3) all have valuation strictly above \(R\).

## Block weights and divided differences

Put \(x=Lq\), and divide the terms in its block by \(C(Lq)\). Their relative weights are
\[
 D_b(x)=\frac{Z^{-b}}{(J-x)_{\underline b}(x+b)_{\underline b}},
                  \qquad b=0,\ldots,L-1.
\]
Define
\[
 F=\sum_b D_b,\qquad
 g=\frac{\sum_b bD_b}{F},\qquad h=\frac{\sum_b b^2 D_b}{F},\qquad
 p(x)=x+g(x),\qquad k(x)=h(x)-g(x)^2.
\]
The functions \(F-1,g,h,k\), and each of their divided differences of order one and two at distinct points of \(L\mathbb Z_2\), have valuation at least two.

To prove this, take \(b>0\). At \(x=0\) the denominators consist of \(i\) and \(J+1-i\) for \(1\le i\le b\). Since \(J+1\) is divisible by \(L\), their valuations are at most \(u=\lfloor\log_2 b\rfloor\le\tau-1\). In the power series in \(x\) at zero, the coefficient of degree \(m\) in \(D_b\) has valuation at least \(wb+2s_2(b)-mu\), by expanding each inverse factor. The divided difference of order \(s\le2\) of the monomial \(x^m\) vanishes if \(m<s\); otherwise it is the sum of monomials of total degree \(m-s\) in the \(s+1\) nodes. All nodes have valuation at least \(\tau\), so each nonzero contribution has valuation at least
\[
            wb+2s_2(b)-su+(m-s)(\tau-u)\ge2 .
\]
The last inequality uses \(2\lfloor\log_2 b\rfloor\le b\), \(w\ge1\), and \(s_2(b)\ge1\). The coefficients have increasing valuation after evaluation, so the divided differences converge. The sums \(F-1,\sum bD_b,\sum b^2D_b\) inherit the bound. The product rule for divided differences and \(F(1/F)=1\) give the same bounds for the positive-order divided differences of \(1/F\), and then for \(g,h,k\).

Consequently \(p\) satisfies \(v(p(x)-p(y))=v(x-y)\). Moreover the order-two divided difference of \(k\) with respect to \(p\) still has valuation at least two. Indeed its order-one divided differences in the \(p\) variable are \(k[x,y]/p[x,y]\), with brackets denoting divided differences in \(x\). The denominators are units. The difference between \(k[x,y]/p[x,y]\) and \(k[x,z]/p[x,z]\) has valuation at least \(2+v(y-z)\), since \(k[x,y]-k[x,z]=(y-z)k[x,y,z]\) and the same formula applies to \(p\). Divide by \(p(y)-p(z)\) to obtain the claimed second-order bound.

After replacing \(H\) by its common value, the three block coefficients in (3), divided by the common factor \(C(Lq)F(x)\), are exactly
\[
                          1,\qquad p(x),\qquad p(x)^2-p(x)+k(x).       \tag{6}
\]
For one anchor \(u\), subtracting \(p(Lu)\) times the zeroth sum from the first sum makes the corresponding coefficient \(p(x)-p(Lu)\), of valuation \(v(x-Lu)\) away from that anchor. For two anchors, subtract from the third sum the affine interpolation of its coefficient in the \(p\) variable through the two anchor nodes, using the zeroth and first sums. The two affine coefficients are integral: the slope of \(p^2-p\) is integral at integral nodes, and the slope of \(k\) in the \(p\) variable has valuation at least two. The resulting coefficient vanishes at both anchors. Away from them its valuation is exactly the sum of the two distances \(v(x-Lu)\). Indeed the monic quadratic part leaves the product of the two differences in the \(p\) variable; the corresponding interpolation residual of \(k\) is that product times its second divided difference, which has strictly positive valuation.

For every allowed anchor set this gives a sum of valuation strictly above \(R\). Each nonvanishing summand has valuation exactly \(V(q)+|\mathcal A|\tau+\sum_{u\in\mathcal A}v_2(q-u)\). One of \(q=0,1,2\) is outside the anchors, at which \(V(q)\le2W+2\) and the total distance is at most \(2\tau+2\lfloor\log_2 B\rfloor\). Thus the smallest summand valuation is at most \(R\). A unique minimum would give the sum that valuation, a contradiction.

For the displayed example \(Y=2^{W+2}\), the unanchored minimum is zero at \(q=0,1\). For \(q>1\) the carry count is at most \(W+2\), which gives \(V(q)\ge W>0\). The anchor zero therefore leaves \(q=1\) as the unique minimum in (1).

## A limitation of the anchored minima

For every \(W\ge16\), there is an even \(Y\) for which all minima in (1) tie. This concerns the integer costs only.

Choose three zero runs \([a_i,b_i)\) of the binary digits, starting at \(a_1=0\), with
\[
 b_i=a_i+W2^{a_i}+2,\qquad
 a_2=b_1+W2^{b_1}+3,\qquad
 a_3=b_2+W2^{b_2}+a_2+3 .
\]
Fill the gaps \([b_1,a_2),[b_2,a_3)\) with ones, and put a one at the highest bit \(b_3\). This defines an even \(Y\) with \(B=b_3+1\).

In the carry automaton, if the carry state stays fixed and a first-addend bit is optional, choosing zero is strictly cheaper by \(W2^i+2\). Hence a least-cost path for a prescribed sequence of carry states uses a first-addend bit precisely to begin a carry at a zero of \(Y\) or to extend it through a one of \(Y\). Beginning at bit \(a\) costs \(W2^a+1\), propagating through a zero costs \(-1\), and propagating through a one at \(i\) costs \(W2^i+1\). Stopping at a one costs zero.

A carry that begins in one zero run and stops at its first following one has nonnegative cost, with equality when it begins at the first bit of that run. Beginning later increases the starting cost while shortening the carry. Continuing through extra ones and stopping inside the gap only adds positive terms. If the carry bridges an entire all-one gap \([b,a)\) and enters the next zero run, break it instead at \(b\) and restart at \(a\). The bridged path costs
\[
                     (a-b)-W2^b-2
\]
more than the broken path. This follows by subtracting the new starting cost \(W2^a+1\) from \(\sum_{i=b}^{a-1}(W2^i+1)-1\), including the original propagation cost at the first zero \(a\). The two gaps have positive bridge penalties, the second equal to \(a_2+1\). Breaking all bridges shows \(V(q)\ge0\) for every \(q\). Its zero-cost addends are exactly the eight subset sums of \(2^{a_1},2^{a_2},2^{a_3}\).

For any choice of at most two anchors the weighted minimum is at most \(a_2\). With one anchor, or two anchors of the same parity, use a zero-cost addend of opposite parity. With two anchors of opposite parity, choose one parity class. Within that class choose the bit at \(a_2\) opposite the corresponding anchor. The opposite-parity distance has valuation zero, and the other distance has valuation at most \(a_2\).

Every weighted minimizer therefore has \(V(q)\le a_2\). Its carry path cannot bridge the final gap, since breaking that bridge alone would save at least \(a_2+1\) and leave a valid path of nonnegative cost. The path enters the final zero run with carry zero. Holding the bits below \(a_3\) fixed, the two cheapest continuations either avoid a carry or carry through that entire run. They give two addends differing by \(2^{a_3}\), both at most \(2^{a_3+1}-1<B\), and preserve all anchor distances of valuation at most \(a_2<a_3\). Each minimizes the weighted cost, as replacing the original continuation cannot increase its cost. Thus every permitted anchor choice still has a tie, despite the exact three-moment interpolation constraints.

## Exact checks

The [finite replay](verify_binary_even_quadratic_anchored_borrow.py) checks 684 exact derivative ratios and identities, 4,718 block valuations and divided differences, 378 interpolation identities, and 12 borrow fixtures. The inequalities above apply to all orders.
