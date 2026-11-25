# Cube-divisor bounds for original noncoprime square-full sums

Let r(N) count ordered positive square-full integers u,v with u+v=N. For fixed positive integers kappa_1,kappa_2, let R_kappa(n) count positive representations with coprime summands
\[
 \kappa_1a^2b^3+\kappa_2c^2d^3=n,
 \qquad b,d\text{ squarefree}.
\]
The coefficient-uniform fixed-target estimate of [Theorem 3](ERDOS_943_FIXED_TARGET_LINES_AND_SEVENTH_POWER_SUMS.md) is
\[
 R_\kappa(n)\ll_\epsilon n^\epsilon
 \left(q+\frac{n^{1/2}}
 {(\kappa_1\kappa_2)^{1/4}q^{3/4}}\right)
 \qquad(q\ge1,\ q^3\mid n).                            \tag{1}
\]
Its q=1 branch follows by counting the first square-full summand; no old three-tenths theorem is used.

**Theorem 1 (Original square-full sums).** For EVERY N>=2 and EVERY integer q>=1 satisfying q^3|N,
\[
 r(N)\ll_\epsilon N^\epsilon
                 \left(q+\frac{N^{1/2}}{q^{3/4}}\right). \tag{2}
\]
There is no extra squarefree or radical-divisibility hypothesis. The estimate does not assert that a suitably sized cube divisor exists for every target.

For an actual original pair set h=gcd(u,v), n=N/h. Each prime of h has valuation at least two, since both original summands are square-full. Thus h is square-full: for every prime p,
\[
                  t=v_p(h)\text{ is }0\text{ or }t\ge2. \tag{3}
\]
The quotients u/h,v/h are coprime. Remove from each its primes of exponent exactly one; their products are kappa_1,kappa_2. These products are squarefree and coprime and divide rad(h). What remains is square-full, with its unique square-times-squarefree-cube expression. This is an injection into R_kappa(n), for at most tau(N)^3 normalization choices. Only the injection, not a converse, is used.

Write E=v_p(N), l=v_p(q), and t=v_p(h). For p|q, E>=3l. Define
\[
 j_p=\min\left(l,\left\lfloor(E-t)/3\right\rfloor\right),
 \qquad q_h=\prod_{p\mid q}p^{j_p}.                     \tag{4}
\]
Then q_h^3|N/h and q_h<=q. The crucial inequality is
\[
                         q_h\ge qh^{-1/2}.             \tag{5}
\]
For p|q it suffices to prove j_p>=l-t/2. If t=0, the minimum in (4) is l. If t=2 or3, the floor is at least l-1, which is at least l-t/2. If t>=4, then
\[
 \left\lfloor(E-t)/3\right\rfloor
 \ge (E-t-2)/3
 \ge l-(t+2)/3
 \ge l-t/2.
\]
The other entry l in the minimum has the same lower bound. These cases exhaust (3). Primes of h not dividing q only make the right side of (5) smaller. This proves (5).

The forbidden case t=1 is EXACTLY where the naive floor argument could fail when E=3l. Its exclusion is supplied by the actual square-full gcd; it is not true for arbitrary integer gcd normalizations. Consequently requiring q^3 rad(q)|N would be sufficient but unnecessarily strong here.

Apply (1) with q_h, including its q_h=1 case. Since kappa_i>=1, the coefficient factor can be discarded in an upper bound. Equations (4)--(5) give
\[
\begin{split}
 R_\kappa(N/h)
 &\ll_\epsilon N^\epsilon
  \left(q+\frac{(N/h)^{1/2}}{(qh^{-1/2})^{3/4}}\right)\\
 &=
 N^\epsilon\left(q+\frac{N^{1/2}}{q^{3/4}}h^{-1/8}\right)
 \le
 N^\epsilon\left(q+\frac{N^{1/2}}{q^{3/4}}\right).
\end{split}
\]
If n=1 there is no positive pair. Summation over divisor-many h,kappa_1,kappa_2 and absorption of the divisor factors proves (2).

For N=M^7 choose q=M^2, recovering r(M^7)<<M^(2+epsilon). More generally, if M^7|N, the same q gives
\[
              r(N)\ll_\epsilon
                 N^\epsilon(M^2+N^{1/2}M^{-3/2}).
\]
Theorem (2) is stronger than this particular divisor-family corollary.

Balancing the two terms in (2) gives q=N^(2/7). For any eta>=0, a cube divisor in
\[
 N^{2/7-\eta}\le q\le N^{2/7+\eta}
\]
gives r(N)<<N^(2/7+eta+epsilon). More generally, for any 0<=delta<=1/70, a cube divisor in
\[
 N^{4/15+4delta/3}\le q\le N^{3/10-delta}
\]
gives r(N)<<N^(3/10-delta+epsilon). These are sufficient target conditions, not assertions that every target has such a divisor. Cube-free targets have only q=1 and receive only the elementary square-root scale from this route.

The only counting input is (1), with the actual gcd normalization and valuation arguments above. No old arbitrary-target three-tenths bound or cycle-frequency hypothesis is used. No new all-target exponent is asserted.
