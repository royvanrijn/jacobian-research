# The prime-only entropy test at polynomial dilation scales

Write
\[
 C_2(X)=\sum_{\substack{X<r\le2X\\r,r+2\ {\rm prime}}}\frac1r,\qquad
 \mathfrak S=2\prod_{q>2,\ q\ {\rm prime}}\left(1-\frac1{(q-1)^2}\right).
\]
Let \(H\) be a positive even integer tending to infinity, let \(\mathcal P_H\) consist of the primes in \((H/2,H]\), and put
\[
 L_H=\sum_{p\in\mathcal P_H}\frac1p,\qquad
 R_H=HL_H,\qquad b_{X,H}=\frac{R_H}{(\log X)^2},\qquad
 \rho_H=\frac{\sum_{p\in\mathcal P_H}p^{-2}}{L_H}.
\]
Choose \(N\) in \([\lceil XH/3\rceil,\lfloor3XH\rfloor]\) with law \(1/(ZN)\), where \(Z\to\log 9\). For \(1\le j\le H\) and \(p\in\mathcal P_H\), put \(m=N+j\). Define the prime-only column \(A_p=(A_{j,p})_{j=1}^H\) by restricting \(m\) to \(pX<m\le2pX\) and setting its entry to one in either of these cases:

- \(p\nmid m\), and \(m,m+2p\) are prime.
- \(p\mid m\), and \(m/p,m/p+2\) are prime.

The entry is zero otherwise. For \(H=o(X)\) and for \(H>4X+4\) this agrees with the [permissive prime-only column](ROUGH_TWIN_PARITY_ENTROPY.md), since the primes produced by division in the second case are distinct from \(p\). The definition here uses a single division and gives exact identities even when \(p\asymp X\).

## The information cost of the detection test

Put \(Y_p=N\bmod p\). Let \(\widetilde Y_p\) be independent of \(A_p\) with the same marginal law as \(Y_p\). The same test under both probability laws is
\[
 F_p(A_p,v)=\sum_{j=1}^H A_{j,p}\mathbf1_{v=-j\bmod p}.
\]
Each nonzero entry has \(N+j\) odd for large \(X,H\). Two distinct indices in the same residue class would differ by \(p\), because \(p>H/2\), and \(p\) is odd. Thus \(F_p\in\{0,1\}\) on every column in the support and every residue. No randomization is needed. Write
\[
 a_p=\Pr(F_p(A_p,Y_p)=1),\quad
 b_p=\Pr(F_p(A_p,\widetilde Y_p)=1),\quad
 A=\sum_p a_p,\quad B=\sum_p b_p,\quad
 \mathcal I(X,H)=\sum_p I(A_p;Y_p).
\]

**Theorem.** For \(B>0\), with \(\phi(u)=1-u+u\log u\) and \(\phi(0)=1\),
\[
 \mathcal I(X,H)\ \ge\
 A\log\frac AB+(|\mathcal P_H|-A)\log\frac{|\mathcal P_H|-A}{|\mathcal P_H|-B}
 \ \ge\ B\phi(A/B).                                                \tag{1}
\]
The formula extends to zero probabilities by the usual limits. In particular, if \(A=o(B)\), the cost is at least \((1-o(1))B\).

Indeed, the test is a function of the column and the tested residue. Data processing on the joint law against the product of its marginals gives the Bernoulli relative entropy with parameters \(a_p,b_p\). Sum and apply the log-sum inequality to both outcomes. The second inequality uses \(s\log(s/t)\ge s-t\).

Write the resampled prime-shift mean as
\[
 M(X,H)=\frac1{L_H}\sum_{p\in\mathcal P_H}\frac1p
             \sum_{\substack{pX<r\le2pX\\r,r+2p\ {\rm prime}}}\frac1r.
\]
The detection masses are
\[
 A=\frac{R_H}{Z}C_2(X)(1+O(1/X)),\qquad
 B=\frac{R_H}{Z}\big(M(X,H)+\rho_H C_2(X)\big)(1+O(1/X)).              \tag{2}
\]
For the actual residue, \(N=pn-j\) and \(1/N=(1+O(1/X))/(pn)\) uniformly. For an independent uniform residue, every entry has detection probability \(1/p\). Reindex \(N=m-j\), replace its weight by \((1+O(1/X))/m\), and use \(m=pn\) on the divisible part. Every index from the reindexing lies in the sampling interval for large \(X\), since \(p>H/2\). Finally \(\Pr(Y_p=v)=1/p+O(1/(XH))\) by harmonic summation on a progression. Replacing uniform residues by the true marginals changes the independent mass by a factor \(1+O(1/X)\). All replacements in (2) apply to nonnegative summands and are independent of the size of \(H\).

## A computable resampled mean

**Theorem.** For fixed \(\theta>8/25\) and \(H=2\lfloor X^\theta/2\rfloor\), the prime-restricted mean satisfies
\[
 M(X,H)\sim\frac{\mathfrak S\log 2}{(\log(XH))^2}.                    \tag{3}
\]

Put \(Y=XH\). The averaged prime-pair theorem of Matomäki, Radziwiłł and Tao evaluates \(\sum_{U<r\le2U}\Lambda(r)\Lambda(r+h)\) as \(\mathfrak S(h)U+O_D(U(\log U)^{-D})\) for all but \(O_D(K(\log U)^{-D})\) shifts \(|h|\le K\), whenever \(U^{8/33+\varepsilon}\le K\le U^{1-\varepsilon}\) with fixed \(\varepsilon>0\). On our range \(\theta/(1+\theta)>8/33\), so we may take \(K=2H\) for every \(U\) between \(Y/(\log Y)^5\) and \(2Y\), using one fixed sufficiently small \(\varepsilon\). Here \(\mathfrak S(2p)=\mathfrak S(p-1)/(p-2)=\mathfrak S+O(1/H)\) for every \(p\in\mathcal P_H\). [Matomäki, Radziwiłł and Tao, Theorem 1.3(i) and the singular series in (7)](https://arxiv.org/pdf/1707.01315v3)

Take grid points \(t\) on \([Y/2,2Y]\) with spacing \(O(Y/\log Y)\), including the endpoints. At each grid point, sum the dyadic estimates on \((t/2^{j+1},t/2^j]\) down to the smallest lower endpoint still at least \(Y/(\log Y)^5\). The remaining prefix has length at most \(2Y/(\log Y)^5\), so its trivial contribution costs \(O(Y/(\log Y)^3)\). With \(D=10\), the union of exceptional shifts across all these grids and dyads has size \(O(H/(\log Y)^8)\). The errors on the retained dyadic intervals sum to \(o(Y)\) at each grid point. Monotonicity of the weighted prefix and boundedness of \(\mathfrak S(2p)\) give
\[
 \sum_{r\le t}\Lambda(r)\Lambda(r+2p)=\mathfrak S(2p)t+o(Y)
\]
uniformly for \(Y/2\le t\le 2Y\) and every prime shift outside the union. This step allows the endpoints \(pX,2pX\) to depend on \(p\).

Removing proper prime powers costs \(O(\sqrt{Y}(\log Y)^3)\) uniformly in the prefix. Partial summation on \((pX,2pX]\) with weight \(1/r\) gives the weighted main term \(\mathfrak S(2p)\log 2+o(1)\). On this interval, both prime logarithms are \(\log Y+O(1)\). Dividing by \((\log Y)^2\) gives the unweighted-prime harmonic sum in (3), uniformly outside the exceptional shifts. For the exceptional primes that harmonic sum is \(O(1)\). Their total weight in the average is \(O((\log Y)^{-7})\) by the bound on the exceptional set and \(L_H\sim(\log2)/\log H\). This proves (3).

Merikoski's Corollary 2 gives an unconditional lower bound over all gaps \(2k\) with \(k\le E\) when \(E\ge C\log x\), \(C>1/2\). Its all-shift lower bound does not control the prime subsample. In (3) the power of the logarithm in the exceptional-set bound is what allows the restriction to prime shifts. [Merikoski, Averaged Form of the Hardy–Littlewood Conjecture](https://arxiv.org/pdf/1605.04757v1)

## Two scales force an information loss

**Theorem.** Put \(c_\theta=\mathfrak S\log 2/(1+\theta)^2\). For any distinct fixed exponents \(\theta_1,\theta_2>8/25\), set \(H_i=2\lfloor X^{\theta_i}/2\rfloor\). Then
\[
 \frac{\mathcal I(X,H_1)}{b_{X,H_1}}+\frac{\mathcal I(X,H_2)}{b_{X,H_2}}
 \ \ge\ \frac{(\sqrt{c_{\theta_1}}-\sqrt{c_{\theta_2}})^2}{\log9}+o(1).       \tag{4}
\]
In particular \((\theta_1,\theta_2)=(1/3,1/2)\) gives a lower bound \(\mathfrak S\log 2/(144\log 9)+o(1)\).

In fact \(\rho_H\le2/H\) and the twin upper sieve gives \(t_X=(\log X)^2 C_2(X)=O(1)\). By (2) and (3),
\[
 \frac{A}{b_{X,H}}=\frac{t_X}{\log9}+o(1),\qquad
 \frac{B}{b_{X,H}}=\frac{c_\theta}{\log9}+o(1).
\]
Thus (1) gives \(\mathcal I(X,H)/b_{X,H}\ge c_\theta\phi(t_X/c_\theta)/\log9+o(1)\), with uniform continuity also at \(t_X=0\). For \(c_1,c_2>0\), the function \(c_1\phi(t/c_1)+c_2\phi(t/c_2)\) on \(t\ge0\) is minimized at \(t=\sqrt{c_1c_2}\) and has minimum \((\sqrt{c_1}-\sqrt{c_2})^2\). This proves (4).

For any fixed collection of \(m\) banks, use the normalization \(\beta_i=b_{X,H_i}/Z_i\), where \(Z_i\) is the harmonic sampling normalizer for that bank. Then \(A_i/\beta_i=t_X+o(1)\) and \(B_i/\beta_i=c_{\theta_i}+o(1)\). The minimizing value of \(t\) for \(\sum_i c_{\theta_i}\phi(t/c_{\theta_i})\) is \((\prod_i c_{\theta_i})^{1/m}\), so
\[
 \sum_{i=1}^m\frac{\mathcal I(X,H_i)}{\beta_i}
 \ \ge\ \sum_{i=1}^m c_{\theta_i}
             - m\left(\prod_{i=1}^m c_{\theta_i}\right)^{1/m}+o(1).
\]
In this normalization the two-bank example has constant \(\mathfrak S\log 2/144\). With \(b_{X,H_i}\) in the denominators as in (4), the constant has the factor \(1/\log9\).

A single fixed polynomial bank with \(\mathcal I=o(b_{X,H})\) would force \(t_X\to c_\theta\), whereas the Hardy–Littlewood prediction at height \(X\) is \(t_X\to\mathfrak S\log 2\). The factor \((1+\theta)^{-2}\) in the resampled density comes from evaluating both primes at height \(XH\). The same column-wise proof as in the smaller bank can still be used here for the transfer error. The total uniform-residue variance is at most the sum of the detection probabilities in (2), which is \(O(b_{X,H})\). The per-column marginal divergences are \(O(X^{-2})\), with sum \(o(b_{X,H})\). Formula (4) identifies a limitation of the unadjusted independence target even when the resampled average can be computed.

## A one-bank obstruction from the full detection pattern

The Bernoulli lower bound (1) makes \(\mathcal I\le(1-\varepsilon)B\) a sufficient condition for a quantitative twin-prime lower bound whenever \(B\asymp b_{X,H}\). Indeed \(\phi(0)=1\), so the condition forces \(A/B\gg_\varepsilon1\). Under the Hardy–Littlewood prediction the ratio \(A/B\) tends to \((1+\theta)^2\) at polynomial scale, and \(\phi((1+\theta)^2)<1\) for \(0<\theta<\sqrt e-1\). This calculation concerns only the lower bound from the binary output. The full column contains substantially more information.

Let
\[
 T_p(N)=\{-j\bmod p:1\le j\le H,\ A_{j,p}(N)=1\}
\]
be the set of residues that would fire the detector. The residues in this set occur without multiplicity by the parity argument before (1). Write \(K_p(N)=|T_p(N)|=\sum_{j=1}^HA_{j,p}(N)\).

**Theorem.** For every fixed \(\theta>8/25\) and \(H=2\lfloor X^\theta/2\rfloor\),
\[
 \sum_{p\in\mathcal P_H} I(T_p(N);Y_p)\asymp_\theta H,\qquad
 \mathcal I(X,H)\asymp_\theta H,\qquad
 \frac{\mathcal I(X,H)}{B}\asymp_\theta(\log X)^3.                  \tag{5}
\]
Thus the one-bank bound \(\mathcal I\le(1-\varepsilon)B\) fails on every fixed polynomial scale in this range. The same lower bound on information holds for any statistic of the column from which the exact outcome at every tested residue can be read.

For the sieve estimates, take \(z=H^\kappa<X\) with a sufficiently small fixed \(\kappa>0\), also small enough that \(z<H/2\). Every active entry forces both \(N+j\) and \(N+j+2p\) to avoid prime factors \(q\le z\) with \(q\ge7\). This holds on the divisible branch as well, since its endpoints are \(pn,p(n+2)\) with \(n,n+2>X>z\) prime and \(p>z\). The two roots are distinct modulo each sieve prime, since their difference is \(2p\). The resulting dimension-two upper sieve on an interval of length \(H\) gives
\[
 K_p(N)\ll_\theta\frac H{(\log X)^2}                              \tag{6}
\]
uniformly over \(N,p\).

Suppose the same residue set is observed at \(N\) and \(N'\), with \(N\not\equiv N'\pmod p\). Each active \(j\) at \(N\) matches a unique active \(j'\) at \(N'\) with \(j'-j=sp\) for some \(s\in\{-1,0,1\}\). Fix \(s\) and put \(d=N'-N+sp\). In this case \(p\nmid d\), so \(0,2p,d,d+2p\) are distinct. The four integers at those offsets from \(N+j\) all avoid the sieve primes. Put \(D=\prod_{a<b}|a-b|\) over these four offsets. Then \(1\le D\ll(XH)^6\) and for every prime \(7\le q\le z\) with \(q\nmid D\) the four forbidden classes are distinct. The same upper sieve gives, for each \(s\),
\[
 \#\{j:\text{both entries are active and }j'-j=sp\}
 \ll_\theta H\left(\frac{\log\log X}{\log X}\right)^4.              \tag{7}
\]
For comparison, if the entire columns at two distinct starts agree, take \(j'=j\). The offsets then have at least three distinct values even when \(p\mid N'-N\), and the corresponding upper bound is \(O_\theta(H(\log\log X/\log X)^3)\).

For completeness, the upper sieve in (6)–(7) uses \(k=2,3,4\) distinct classes at each retained prime. The Selberg argument with squarefree weight cutoff \(R=H^{1/8}\) gives a main bound \(O(H\prod_q(1-k/q))\) and error \(O_k(R^2(1+\log R)^{2k-2})\). To see the uniformity, interval counting by the Chinese remainder theorem gives an error at most \(k^{\omega([d,e])}\) for weight indices \(d,e\). The sum of these errors is at most \((\sum_{d\le R}\tau_k(d))^2\). For the main term, the Euler-product probability law of the Selberg denominator chooses each retained prime independently with probability \(k/q\). Its expected logarithm is \(O_k(\log z)\), at most \((\log R)/2\) for a small enough \(\kappa\), so truncation preserves at least half its mass. For the collision bound, applying Mertens to all primes in \([7,z]\) gives \(O_k((\log z)^{-k})\) before removing the divisors of \(D\). Their inverse factors cost \(O_\theta((\log\log X)^k)\). The primes at most \(\log X\) contribute at most that much by Mertens. Above \(\log X\), at most \(O_\theta(\log X/\log\log X)\) primes divide \(D\), with reciprocal sum \(O_\theta(1/\log\log X)\). The sieve error is smaller than the main upper-bound scale. [Selberg weights and Mertens input](SQUAREFREE_TWIN_CORRELATION_NEAR_ZERO.md)

By (7), there is a threshold
\[
 u_X=O_\theta\big(H(\log\log X/\log X)^4\big)=o\big(H/(\log X)^2\big)
\]
such that any residue set of size \(>u_X\) determines \(Y_p=N\bmod p\) exactly. The resampled detection probability satisfies \(b_p=(1+O(1/X))\mathbb E K_p/p\). By (2)–(3), \(B\asymp_\theta|\mathcal P_H|/(\log X)^2\), and so \(\sum_p\mathbb E K_p\gg_\theta H|\mathcal P_H|/(\log X)^2\). Using (6) to bound \(K_p\) on the event \(K_p>u_X\) gives
\[
 \sum_{p\in\mathcal P_H}\Pr(K_p>u_X)\gg_\theta |\mathcal P_H|.
\]
On such a set the conditional law of \(Y_p\) given \(T_p\) is a point mass. Its relative entropy against the marginal law is \(\log p+O(1/X)\), by the progression estimate used in (2). The other conditional relative entropies are nonnegative. Thus the sum of mutual informations of the residue sets is \(\gg_\theta |\mathcal P_H|\log H\asymp H\). Data processing from \(A_p\) to \(T_p\) gives the same lower bound for the full columns, and \(\sum_p\log p\ll H\) bounds both sums above. Dividing by \(B\) proves (5).

## Applicability of the weighted-dilation graph

The averaged prime-pair theorem applies before the selection on the residue \(p\mid m\). On the actual selected residue the corresponding endpoints are \(pn,p(n+2)\), each with at most two distinct prime factors. On the off-residue part they each have just one. In the weighted-dilation Theorem 5.1 the factor \(W_\ell\) vanishes unless the endpoint has at least \(\ell\) distinct primes in each of \(\asymp\log\log x\) disjoint groups. It vanishes identically on both supports. Definition 3.1 also places the group primes below \(\exp((\log x)^d)\) for fixed \(d<0.47\) and requires the harmonic mass of every group to be bounded below by a fixed positive constant. At \(x\asymp XH\), the polynomial bank violates the size restriction, and its harmonic mass \(L_H\sim(\log2)/\log H\) tends to zero. A graph estimate that could exploit (3) would need to work on these low-factor endpoints and allow for the change in prime density between the two heights. [Weighted-dilation preprint, Definition 3.1 and Theorem 5.1](sources/weighted-dilation-graphs-sep24-2026.pdf)
