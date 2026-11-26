# Same-orientation block gcds for common-factor denominators

Consider the ordinary moment polynomials \(P_n,P_N\) of any two
distinct exponent columns, with \(N=n+1\), \(n\ge2\). The exposition
uses the original common argument. Rows with equal exponents can be
removed. Write each remaining row in its native orientation, with
smaller exponent \(a_i\), absolute gap \(d_i\), and orientation
\(+\) when its exponent increases from the first to the second column,
or \(-\) when it decreases.

For an irreducible rational common factor \(H\) of degree \(d\), let
\(L_H\) be its [projective denominator](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#exact-projective-denominators)
and let \(D_H\) be the [denominator of its root-quotient polynomial](MULTIVARIATE_DEEP_BINOMIAL_PROJECTIVE_QUOTIENTS.md#exact-quotient-denominator).
The even moment has no real zeros, so \(d\) is even. The following
gcds bound both denominators using the exact factorial blocks.

For \(j=n,N\) and \(1\le k\le j\), put
\[
\begin{aligned}
 U_{j,k}^+&=(j-k+1)\prod_{i:\,+}
                  \prod_{h=(k-1)d_i+1}^{kd_i}(ja_i+h),\\
 U_{j,k}^-&=k\prod_{i:\,-}
                  \prod_{h=(j-k)d_i+1}^{(j-k+1)d_i}(ja_i+h),\\
 B_{j,r}^\pm&=\operatorname{lcm}_{|I|=r}
                              \prod_{k\in I}U_{j,k}^\pm,\\
 G_r^\pm&=\gcd(B_{n,r}^\pm,B_{N,r}^\pm),\qquad 0\le r\le n .
\end{aligned}
\]
The index set \(I\) ranges through subsets of \(\{1,\ldots,j\}\).
In particular \(B_{j,0}^\pm=G_0^\pm=1\).

**Theorem.** With \(q=d/2\) one has the exact divisibility bounds
\[
 \begin{split}
 L_H&\ \mid\
 \operatorname{lcm}_{0\le k\le d}
         (G_k^+)^{d-k}(G_{d-k}^-)^k,\\
 D_H&\ \mid\
      G_q^+G_q^-\prod_{k=1}^{q-1}(G_k^+G_k^-)^2.                \tag{1}
 \end{split}
\]
Consequently a prime \(p>N\) can occur in either denominator only if
one orientation has a factorial factor divisible by \(p\) in
**both** parent orders. This support condition still permits the
primes, heights, orientations and row counts to vary with the profile.

## Means of consecutive blocks of valuations

Write \(c_{j,k}\) for the coefficient at index \(k\) of \(P_j\).
The exact ratio identity
\[
                  \frac{c_{j,k}}{c_{j,k-1}}
                            =\frac{U_{j,k}^+}{U_{j,k}^-}        \tag{2}
\]
follows from the ratio of the binomial coefficients and a single
factorial block in each row. In the decreasing orientation the
removed block has the reversed index shown in \(U_{j,k}^-\).

Fix a prime \(p\). A segment of the lower coefficient Newton polygon
on the indices \(l,\ldots,u\) has slope the mean of the consecutive
first differences of the ordinates between \(l\) and \(u\).
List that slope \(u-l\) times. The resulting \(j\)-entry list is
obtained from
\((v_p(c_{j,k}/c_{j,k-1}))_{k=1}^j\) by replacing disjoint intervals
by their means. Thus a sum of any \(r\) entries of the list is at
most the sum of the \(r\) largest original first differences.
Indeed a sum using \(t\) entries from an interval of width \(b\)
is the original sum on that interval weighted by \(t/b\).
The weights are in \([0,1]\) and sum to \(r\), so the maximum
is achieved by giving weight one to the \(r\) largest entries.
The same comparison applies to the negatives of the slopes.

The Newton polygon root theorem in its ordinary case makes the
negatives of these slopes precisely the root valuations, with
multiplicity. For a common factor use its valuations
\(s_1\le\cdots\le s_d\) in a single \(p\)-adic algebraic closure.
They occur in each parent's list. Equation (2) and the same averaging
comparison now give
\[
 \begin{split}
  -\sum_{i=1}^r s_i &\le
             \min_{j=n,N}\max_{|I|=r}\sum_{k\in I}v_p(U_{j,k}^+)
                 =v_p(G_r^+),\\
  \sum_{i=d-r+1}^d s_i &\le
             \min_{j=n,N}\max_{|I|=r}\sum_{k\in I}v_p(U_{j,k}^-)
                 =v_p(G_r^-).                                  \tag{3}
 \end{split}
\]
The bounds take the minimum *separately* for each orientation
and number \(r\).

For the first line of (1), use
\[
 v_p(L_H)=\max_{0\le k\le d}
       \left((d-k)\Big(-\sum_{i=1}^k s_i\Big)
                                +k\sum_{i=k+1}^d s_i\right)
\]
and insert (3). For the second line, use
\[
\begin{aligned}
 v_p(D_H)
   &=\sum_{i<j}|s_i-s_j|
     =\sum_{i=1}^q(d+1-2i)(s_{d+1-i}-s_i)\\
   &=S_q+2\sum_{r=1}^{q-1}S_r,\qquad
            S_r=\sum_{i=d-r+1}^d s_i-\sum_{i=1}^r s_i .
\end{aligned}
\]
Again (3) gives the result, prime by prime. If \(p>N\) it sees none
of the binomial ratio factors \(j-k+1\) and \(k\). If in each
orientation one of the parent factorial strips is a \(p\)-unit,
all \(G_r^\pm\) have valuation zero.

One can cancel within each ratio in (2) before computing the block
LCMs. Replace \(U_{j,k}^+,U_{j,k}^-\) by the reduced numerator and
denominator of that ratio. The resulting \(G_r^\pm\) can only
decrease by divisibility, and (3) and (1) still hold. The native
uncancelled version displays the separate orientation budgets and
can be computed without forming the moment polynomial.

## Opposed types at independently varying heights

Suppose all unequal rows consist of \(R\) increasing copies of
\((a,a+d_+)\) and \(S\) decreasing copies of \((b+d_-,b)\).
There is no relation between \(a,b,d_+,d_-,R,S\).
At a prime \(p>nN\max(d_+,d_-)\), a factor \(na+h\) in the
increasing strip of order \(n\) and a factor \(Na+h'\) in its
order-\(N\) strip can both be divisible by \(p\) only if
\[
                        Nh=nh'.
\]
Indeed \(Nh-nh'\) is divisible by \(p\) and its absolute value is
less than \(nN d_+\). Thus \(h=nv,\ h'=Nv\) with \(1\le v\le d_+\).
The primed and unprimed multiplier factors \(n,N\) are \(p\)-units.
Each strip contains at most one such nonunit because its length
is less than \(p\). Consequently, for every \(k\ge1\) the \(p\)-part
of \(G_k^+\) is
\[
       p^{R v_p(\prod_{v=1}^{d_+}(a+v))},
\]
and the corresponding \(G_k^-\) is given by \(S,b,d_-\).
By (1), above \(nN\max(d_+,d_-)\) the full quotient denominator
divides the corresponding prime part of
\[
              \left(\prod_{v=1}^{d_+}(a+v)^R
                      \prod_{v=1}^{d_-}(b+v)^S\right)^{d-1}.
                                                               \tag{4}
\]
This restriction is independent of any congruences linking the
two base heights.

## Consequences at a positive diffuse contact

Fix a core and a compact positive strength interval in the
[physical diffuse comparison](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md).
Let \(B_H\) be the degree-\(d\) factor's projective polynomial, and
suppose its limiting \(B_{H_*}\) at an algebraic reference branch
has an irrational coefficient of degree \(e\). Use the
[projective sieve's matching exponent](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#consequences-at-finite-physical-diffuse-contacts)
\(\gamma\). Then (1) and the height estimate of that sieve imply
\[
  \log\operatorname{lcm}_{0\le k\le d}
               (G_k^+)^{d-k}(G_{d-k}^-)^k
                 \ \ge\ \frac{\gamma}{e}\log(1/\eta)-O(1).     \tag{5}
\]
Thus any uniformly smaller exponent in the aligned-orientation
budget excludes that branch, even if the large primes in the
ramps change and the gap and row counts grow. If every candidate
reference selection of every possible degree has such an irrational
coordinate, a uniform exponent strictly below the corresponding
finite minimum gives one positive exclusion tolerance.

The second line of (1) also bounds the denominator on the right
of the [depth-weighted binomial-face inequality](MULTIVARIATE_DEEP_BINOMIAL_PROJECTIVE_QUOTIENTS.md#improvement-near-a-finite-positive-diffuse-contact).
For instance, with
\(
 V=G_q^+G_q^-\prod_{k=1}^{q-1}(G_k^+G_k^-)^2
\)
and that result's notation, a common factor must satisfy
\[
 \sum_p U_{p,\ell}\log p+
       2\gamma(t_{\ell,*}-t_\ell)\log(1/\eta)
                    \ \le\ \ell\log V+O(1).                   \tag{6}
\]
The denominator budget on the right charges primes above \(N\)
only for same-orientation factorial overlap across the two orders.
The bound can be applied
alongside the direct shared-slope budget of the depth-weighted
test. It gives profile-dependent restrictions; the aligned
factorial block gcds are not yet controlled for every profile.

## Exact comparison with the ratio-range budget

Take orders three and four and rows \((65,68)\) and \((77,74)\).
At \(p=101\), the increasing strip at order three contains 202
at offset seven, and the decreasing strip at order four contains
303 at offset seven. No other strip factor is divisible by \(p\).
The two lists of root valuations are \((-1,0,0)\) and
\((0,0,1/2,1/2)\), so a degree-two common candidate is locally
permitted, with both valuations zero. Both parent lists of first
differences have range one, so the gcd obtained from the lcm of
numerators of all cross ratios of neighboring coefficients has
valuation one. In (1) every \(G_r^\pm\) is a \(101\)-unit because
the two strip hits have opposite orientations. Thus the new
bound excludes this prime from both projective denominators
without relying on its size or on a presumed global factor.
