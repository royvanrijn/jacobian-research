# Exact-five carriers in sums of two 3-full integers

Let \(\mathcal F_3\) be the positive integers whose nonzero prime valuations are at least three, including \(1\), and let \(r_3(N)\) count ordered positive pairs from this set with sum \(N\). Put
\[
 \mathcal P_{3,4}=\{a^3b^4:a,b\ge1\},\qquad
 v(m)=\prod_{p:\,v_p(m)=5}p.
\]
The integer \(v(m)\) is the exact-five carrier of \(m\). The current [full-\(k\)-full theorem](POWERFUL_SUM_GRAPHS.md) gives \(r_3(N)\ll_\epsilon N^{1/4+\epsilon}\). The following result improves that exponent in restricted carrier sectors; it makes no unrestricted improvement claim.

**Theorem 1.** For every \(N\ge2\), every \(H\ge1\), and every \(\epsilon>0\),
\[
 \#\{(m_1,m_2)\in\mathcal F_3^2:m_1+m_2=N,
                         \ v(m_1)v(m_2)\le H\}
 \ll_\epsilon N^\epsilon
       \{N^{1/5}H+N^{3/14}H^{13/28}\}.                       \tag{1}
\]
For \(H=N^\rho\), its exponent is at most
\[
                  \max\{1/5+\rho,\ 3/14+13\rho/28\},          \tag{2}
\]
which is strictly below \(1/4\) for \(0\le\rho<1/20\). Moreover, uniformly for fixed actual carriers \(v,w\), the representation fiber is
\[
 \ll_\epsilon N^\epsilon
       \max\{N^{1/5},\ N^{3/14}(vw)^{-15/28}\}.              \tag{3}
\]
There is no coprimality condition on the original summands or restriction on the shape of \(N\).

## 1. The carrier bijection

**Lemma 2.** A positive integer is \(3\)-full if and only if, uniquely,
\[
 m=v^5u,\qquad v\text{ squarefree},\quad (v,u)=1,
                                     \quad u\in\mathcal P_{3,4}. \tag{4}
\]
Here \(v=v(m)\). Indeed, \(\langle3,4\rangle=\{0,3,4,6,7,8,\ldots\}\): the integers \(6,7,8\) belong to the semigroup, and adding three covers all larger integers. Among the positive \(3\)-full valuations, only five is missing. Removing exactly the corresponding primes proves (4), and conversely (4) gives precisely that carrier. This differs from the usual canonical exponent-five coordinate, which can contain primes of valuations \(8,11,\ldots\).

Every representation therefore yields a unique ordered pair of squarefree carriers and
\[
                    v^5u+w^5z=N,\qquad u,z\in\mathcal P_{3,4},
                  \quad(v,u)=(w,z)=1.                        \tag{5}
\]
In particular \(v,w<N^{1/5}\). Extra restrictions may be discarded for an upper bound.

## 2. A numerical coefficient-sensitive envelope

Use nonnegative logarithmic scales \(\alpha,\beta,\gamma,\delta\), put
\(X=\alpha+\gamma\), \(Y=\beta+\delta\), and define
\[
\begin{split}
 K&=\max\{\beta,\delta,
       \min(\beta+\delta/2,\delta+\beta/2)+(X-1/2)^+\},\\
 D&=\max\{\alpha,\gamma,
       \min(\alpha+\gamma/2,\gamma+\alpha/2)+(Y-1/2)^+\},\\
 M_1&=\beta+\gamma+\max(\alpha,\delta)/4,\\
 M_2&=\alpha+\delta+\max(\beta,\gamma)/4,\\
 F&=\min\{X,Y,K,D,M_1,M_2\}.
\end{split}                                                   \tag{6}
\]
Here \(t^+=\max(t,0)\).

**Lemma 3.** If \(c_1,c_2\in[0,1]\), \(M=\max(c_1,c_2)\), and
\[
             3\alpha+4\beta\le c_1,\qquad
             3\gamma+4\delta\le c_2,
\]
then
\[
                  F\le\max\{M/5,\ 3(c_1+c_2)/28\}.            \tag{7}
\]
This is a statement about the numerical menu alone.

If \(M=0\), all coordinates vanish. Otherwise put \(m=\min(c_1,c_2)\) and suppose \(F\) exceeds the right side of (7). Both height losses vanish. Indeed \(Y\le(c_1+c_2)/4\le1/2\). If \(X>1/2\), then
\[
 F\le Y<(c_1+c_2-3/2)/4
       \le(2M-3/2)/4\le M/5,
\]
where the last inequality holds for \(M\le1\).

We exclude every individual-coordinate branch of \(K,D\). If \(\beta,\gamma\ge F\) and \(\beta\ge\gamma\), then
\[
 M_2\le(c_1-4\beta)/3+(c_2-3\gamma)/4+\beta/4
           \le(4c_1+3c_2)/12-11F/6.
\]
If \(\gamma\ge\beta\), the same bound follows with coefficients \(-4/3,-1/2\) on \(\beta,\gamma\). Thus \(F\le M_2\) implies \(F\le(4c_1+3c_2)/34\). Interchanging the sides treats the opposite pair \(\alpha,\delta\). Both are at most \((4M+3m)/34\), and
\[
 \frac{4M+3m}{34}\le
 \begin{cases}
 M/5,&m/M\le14/15,\\
 3(M+m)/28,&m/M\ge5/9.
 \end{cases}
\]
The two ranges cover \([0,1]\), so every opposite-individual case contradicts the supposition. Individuals on one side force \(7F\le c_i\le M\).

If \(\beta\ge F\) and the root graph term is mixed, either \(\gamma\ge\beta\), already an opposite-individual case, or
\(\alpha+\beta/2\ge\alpha+\gamma/2\ge F\). Hence
\(c_1\ge3\alpha+4\beta\ge3F+5\beta/2\ge11F/2\), giving \(F\le2M/11<M/5\).
If \(\alpha\ge F\) and the kernel graph term is mixed, either \(\delta\ge\alpha\), already opposite individuals, or
\(\beta+\alpha/2\ge\beta+\delta/2\ge F\). Then
\(c_1\ge3\alpha+4\beta\ge\alpha+4F\ge5F\), giving \(F\le M/5\).
Interchanging sides exhausts the cases.

Consequently both graph terms must be mixed. Each mixed term is at most three quarters of its coordinate sum, so
\[
 X,Y\ge4F/3,\qquad
              28F/3\le3X+4Y\le c_1+c_2.
\]
This is the last contradiction. The strict-case argument proves (7), including equality and all height boundaries.

## 3. From primitive boxes to complete carrier fibers

The [mixed-product proof](ERDOS_943_DUAL_GRAPHS_MIXED_POWER_PRODUCTS.md), Sections 1 and 4--5, supplies canonical \(3\)-powerfree-kernel coordinates and the actual two-graph, fixed-pair and absolutely irreducible mixed-curve estimates. Applied to
\[
 \kappa_1a^3b^4+\kappa_2c^3d^4=n,
 \qquad(\kappa_1a^3b^4,\kappa_2c^3d^4)=1,
\]
the box exponent is bounded by (6). The resources are
\(3\alpha+4\beta\le1-\log_n\kappa_1\) and its right counterpart. Nonempty equations have \(n\ge2\) and \(\kappa_i<n\), so Lemma 3 applies. Summing the logarithmically many boxes gives the primitive value count
\[
 \ll_\epsilon n^\epsilon
       \max\{n^{1/5}\kappa_{\min}^{-1/5},\ n^{3/14}(\kappa_1\kappa_2)^{-3/28}\}. \tag{8}
\]
Here \(\kappa_{\min}=\min(\kappa_1,\kappa_2)\), since \(M=1-\log_n\kappa_{\min}\). Constants are uniform in both positive coefficients. This uses the written actual graph and curve bounds, not a new external input. Dyadic endpoint constants can be treated by the same constant-width boxes, and the logarithmic factors are absorbed into \(n^\epsilon\).

To transfer (8) to (5), let \(h=(v^5u,w^5z)\) and \(n=N/h\). Use the full residue normalization in Lemma 7 and Section 8 of that proof. Outside \(h\), the new coefficients \(\widetilde\kappa_i\) retain the original coefficient prime parts. At a prime dividing \(h\), use the total quotient valuation modulo three for the coefficient, and put the multiple of three into the root. The residual values remain in \(\mathcal P_{3,4}\); the two new full summands are coprime. There are at most
\(\tau(N)5^{\omega(N)}\ll_\epsilon N^\epsilon\) normalization choices, and each normalized value pair recovers the original values uniquely by multiplication by \(h\) and division by the fixed original coefficients.

The carrier conditions retain a coefficient-size gain. If a prime of \(h\) divides exactly one carrier, that summand has valuation five and the other has a positive valuation in \(\langle3,4\rangle\), at least three. Thus \(v_p(h)\ge3\). If it divides both carriers, their residuals are units, so \(v_p(h)=5\). At other primes of \(h\) there is no carrier contribution. Outside \(h\) the original coefficient parts are retained. Hence, prime by prime,
\[
                  v^5w^5\mid h^2\widetilde\kappa_1
                                      \widetilde\kappa_2.     \tag{9}
\]
In particular
\[
 (N/h)^{3/14}(\widetilde\kappa_1\widetilde\kappa_2)^{-3/28}
                     \le N^{3/14}(vw)^{-15/28}.
\]
The first term of (8) is at most \(N^{1/5}\). Summing the subpower number of normalization choices proves (3). No assertion that the extracted coefficients divide \(h\) is needed.

## 4. Small products and restricted prime supports

There are \(\ll H\log(2H)\) positive ordered pairs with \(vw\le H\). Partial summation gives
\[
          \sum_{vw\le H}(vw)^{-15/28}
                         \ll H^{13/28}\log(2H).               \tag{10}
\]
Sum (3), replacing a maximum by a sum. This proves (1). Actual carriers satisfy \(vw<N^{2/5}\), so the cutoff may first be truncated to that range and logarithms absorbed into \(N^\epsilon\). The inequalities in (2) are below \(1/4\) when \(\rho<1/20\), since the second permits the larger range \(\rho<1/13\).

For a separate useful sector, fix \(C>0\) and an integer \(1\le Q\le N^C\). Carriers restricted by \(v\mid Q,w\mid Q\) give
\[
                         \ll_{C,\epsilon}N^{3/14+\epsilon}.   \tag{11}
\]
There are at most \(\tau(Q)^2\ll N^\epsilon\) carrier pairs, and the coefficient-uniform full value count in Theorem 1 of the mixed-product proof gives \(N^{3/14+\epsilon}\) per pair, independently of Lemma 3. This includes the sector in which every exact-five prime of either summand divides \(N\). It restricts the prime support of the whole sector; no such common polynomial-height support is asserted for all representations.

## 5. Large carriers and a conditional global saving

**Proposition 4.** For \(V\ge1\), the representations with \(\max(v,w)\ge V\) number
\[
                              \ll N^{1/3}V^{-2/3}.            \tag{12}
\]
Indeed
\[
 \#(\mathcal P_{3,4}\cap[1,X])
 \le\sum_{b\le X^{1/4}}(X/b^4)^{1/3}
 \le\zeta(4/3)X^{1/3}.
\]
Using (4), the number of possible summands with carrier at least \(V\) is at most
\(O(N^{1/3}\sum_{v\ge V}v^{-5/3})=O(N^{1/3}V^{-2/3})\). For \(V>N^{1/5}\) it is zero. Fixing one summand fixes its complement, and a union bound over the two positions proves (12). This is independent of the mixed-product counting inputs.

Let \(\mathcal C_N(H,V)\) be the ordered carrier pairs occurring in actual representations with \(vw>H\) and \(\max(v,w)<V\), each pair counted once. Define the support capacity
\[
 \mathfrak C_N(H,V)=
  \sum_{(v,w)\in\mathcal C_N(H,V)}
       \max\{N^{1/5},\ N^{3/14}(vw)^{-15/28}\}.              \tag{13}
\]
The weight is a fiber majorant, not the actual fiber multiplicity.

**Proposition 5.** For all \(H,V\ge1\),
\[
 r_3(N)\ll_\epsilon N^\epsilon
 \{N^{1/5}H+N^{3/14}H^{13/28}
              +N^{1/3}V^{-2/3}+\mathfrak C_N(H,V)\}.          \tag{14}
\]
This follows by partitioning into the small-product sector, the large-carrier sector, and the remaining support, and applying (1), (12), and (3).

For \(0<\eta\le1/28\), put
\[
 \rho_\eta=\min\{1/20-\eta,\ (1-28\eta)/13\},\qquad
 H=N^{\rho_\eta},\quad V=N^{1/8+3\eta/2}.
\]
Both entries defining \(\rho_\eta\) are nonnegative. An additional estimate
\[
                \mathfrak C_N(H,V)
                        \ll_\epsilon N^{1/4-\eta+\epsilon}   \tag{15}
\]
would therefore imply \(r_3(N)\ll_\epsilon N^{1/4-\eta+\epsilon}\). Independent epsilon losses are divided in advance. Estimate (15), for any fixed positive \(\eta\) uniformly over all \(N\), remains open. It is only a sufficient condition: a saving for \(r_3(N)\) alone does not imply it, since the weights can substantially exceed their actual fiber sizes.

## 6. The balanced limitation of independent carrier freezing

**Proposition 6.** Independent application of the mixed-product box menu followed by unrestricted carrier-pair summation does not improve the current full-\(3\)-full bound in its balanced numerical box.

In the usual canonical presentation \(m=a_0^3a_1^4a_2^5\), take all three coordinate scales to be \(N^{1/12}\) on both sides. In the subcase \((a_0,a_2)=1\), the exact-five carrier is \(a_2\). After fixing both carriers, the residual scales in (6) are
\(\alpha=\beta=\gamma=\delta=1/12\). Thus
\[
 X=Y=1/6,\qquad K=D=1/8,\qquad M_1=M_2=3/16,
                              \qquad F=1/8.                 \tag{16}
\]
There are at most \(O(N^{1/6})\) carrier pairs, so independent summation gives exponent \(1/8+1/6=7/24>1/4\). This is a calculation about a numerical majorant and feasible scales, not an assertion that any fixed \(N\) realizes all those carrier pairs, or that every way of using their correlations fails.

For primitive actual summands in this residual subbox, an actual support bound of \(N^{1/8-\eta+\epsilon}\) carrier pairs would combine with the fiber estimate \(N^{1/8+\epsilon}\) to give exponent \(1/4-\eta\). That frequency estimate is not proved here. The existing grouped full-\(3\)-full graphs retain information discarded by this independent freezing step. No fixed-difference, global squarefull-triple, or new ordinary-prime theorem is invoked.
