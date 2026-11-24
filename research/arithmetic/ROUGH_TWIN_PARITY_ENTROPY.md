# Entropy decrement for power-rough twin pairs

Fix \(0<\eta<1/2\) and write \(y=X^\eta\). Let \(P_y\) be the product of odd primes at most \(y\), and set
\[
 f_y(n)=\mu(n)\mathbf1_{n\ {\rm odd}}\mathbf1_{(n,P_y)=1}.
\]
Thus \(f_y\) is supported on odd squarefree integers with no prime factor at most \(y\). For a real weight \(|G(a,b)|\le 1\), put
\[
 C(G)=\sum_{X<n\le2X}\frac{f_y(n)f_y(n+2)G(n,n+2)}n,
 \qquad
 U(G)=\sum_{X<n\le2X}\frac{|f_y(n)f_y(n+2)|G(n,n+2)}n.
\]
Write \(G_3(a,b)=\mathbf1_{\omega(a)+\omega(b)\le3}\) and
\(C_2(X)=\sum_{X<r\le2X,\ r,r+2\ {\rm prime}}1/r\).
On the support of \(G_3\) in these sums, the only configurations are two primes or a prime and a squarefree semiprime. Their parity products are \(+1\) and \(-1\), respectively. Hence, for large \(X\),
\[
 C(G_3)+U(G_3)=2C_2(X).                                           \tag{1}
\]
A bound \(C(G_3)\ge -U(G_3)+c/(\log X)^2\) with a fixed \(c>0\) would give \(\gg X/(\log X)^2\) twin primes in \((X,2X]\). The factor-count restriction matters. Two squarefree semiprimes also have parity product \(+1\).

## Entropy comparison at the rough-pair scale

Let \(H\) be an integer tending to infinity with \(H\le\sqrt{\log X}\), and let
\[
 \mathcal P=\{p\ {\rm prime}:H/2<p\le H\},\qquad
 R=H\sum_{p\in\mathcal P}\frac1p\asymp\frac H{\log H},\qquad
 b=\frac{R}{(\log X)^2}.
\]
Choose \(N\) in \([\lceil XH/3\rceil,\lfloor3XH\rfloor]\) with probability \(1/(ZN)\), where \(Z\asymp1\) normalizes the mass. Let \(Y_p=N\bmod p\) and \(V_p\) denote the uniform law on residues modulo \(p\). For \(m=N+j\), \(1\le j\le H\), let \(m^{(p)}=m/p^{v_p(m)}\) and define
\[
 \begin{split}
 f_{y,p}(m)&=\mu(m)\mathbf1_{m\ {\rm odd}}\mathbf1_{(m,P_y/p)=1},\\
 A_{j,p}&=\mathbf1_{pX<m\le2pX} f_{y,p}(m)f_{y,p}(m+2p)
                    G\big(m^{(p)},(m+2p)^{(p)}\big).
 \end{split}
\]
The array allows \(p\) at either endpoint and has no explicit divisibility indicator for \(p\). Write \(A_p=(A_{j,p})_{j=1}^H\) and set, with logarithms in nats,
\[
 d=\sum_{p\in\mathcal P}
   \left(I(A_p;Y_p)+D_{\rm KL}(\mathcal L(Y_p)\Vert V_p)\right).
\]

**Theorem.** Uniformly for these weights,
\[
 C(G)=\frac{Z}{R}\,\mathbb E\sum_{p\in\mathcal P}\sum_{j=1}^H\frac{A_{j,p}}p
   +O_\eta\left(\frac{\sqrt{db}+d}{R}+\frac1X\right).             \tag{2}
\]
In particular \(d=o(b)\) suffices for an error \(o((\log X)^{-2})\). Each marginal relative entropy in \(d\) is \(O((p/(XH))^2)\).

For the exact change of variables, note that \(f_{y,p}(pn)=-f_y(n)\). If \(p\mid N+j\) and both factors in \(A_{j,p}\) are nonzero, stripping \(p\) gives \(n,n+2\). With \(N=pn-j\), the interval for \(N\) contains all the terms and \(1/(pn-j)=(1+O(X^{-1}))/(pn)\). Consequently
\[
 \mathbb E\sum_{p,j}A_{j,p}\mathbf1_{Y_p=-j\bmod p}
       =(R/Z)C(G)+O(R/X).                                       \tag{3}
\]
An identical identity with absolute values applies to \(U(G)\) for \(G\ge0\).

For each fixed column, take \(v\) uniform modulo \(p\) and put
\(S_p(v)=\sum_j A_{j,p}(\mathbf1_{v=-j\bmod p}-1/p)\).
Its mean is zero, \(|S_p(v)|\le4\), and
\(\mathbb E_v S_p(v)^2\le(2/p)\sum_j |A_{j,p}|^2\). At most two indices \(j\) occupy any residue because \(p>H/2\). Taylor's inequality therefore bounds
\(\log \mathbb E_v e^{sS_p(v)}\le C s^2\sum_j|A_{j,p}|^2/p\)
for \(|s|\le1/4\), with an absolute \(C\).

The unsigned sieve bound gives
\[
 \mathbb E\sum_{p,j}|A_{j,p}|^2/p\ll_\eta b.                     \tag{4}
\]
Indeed, \(m,m+2p\) avoid the two distinct classes \(0,-2p\bmod q\) at each odd \(q\le y\) other than \(p\), and \(N\) has mass \(O(1/(XH))\) at each point. The dimension-two Selberg bound on an interval of length \(O(XH)\) is \(O_\eta(XH/(\log X)^2)\), uniformly in \(p\). One may run the [Selberg argument for fixed forbidden classes](SQUAREFREE_TWIN_CORRELATION_NEAR_ZERO.md) up to a sufficiently small fixed power of \(X\) below \(y\), omitting \(p\). Omitting one prime changes the product by a bounded factor.

Apply the relative-entropy inequality to the conditional law of \(Y_p\) given \(A_p\) against \(V_p\) in both signs. Averaging the conditional logarithmic moment bound and summing over \(p\) gives \(s|\mathbb E\sum_p S_p(Y_p)|\le d+C_\eta s^2b\) for \(0<s\le1/4\). Choosing \(s\) of order \(\min(1,\sqrt{d/b})\), with the case \(d=0\) obtained by taking \(s\) to zero, proves (2) using (3). Summation by parts with \(1/N\) shows that each residue modulo \(p\) has mass \(1/p+O(1/(XH))\). A Taylor bound for relative entropy gives the marginal estimate. The argument uses separate conditional laws for each prime and requires no independence between the true \(Y_p\).

For \(G(a,b)=\mathbf1_{a,b\ {\rm prime}}\), one has \(C(G)=C_2(X)\) and an explicit expression for the reference average. The nonzero entries have either \(p\nmid m\), in which case \(m,m+2p\) are prime, or \(p\mid m\), in which case the stripped pair consists of twin primes. Both cases give \(A_{j,p}=1\). Averaging over \(N\), and replacing \(1/(m-j)\) by \(1/m\) as above, gives
\[
 \begin{split}
 (1-\rho_H)C_2(X)
  &=\frac{\displaystyle\sum_{p\in\mathcal P}p^{-1}
          \sum_{\substack{pX<r\le2pX\\ r,r+2p\ {\rm prime}}}r^{-1}}
         {\displaystyle\sum_{p\in\mathcal P}p^{-1}}
    +O_\eta\left(\frac{\sqrt{db}+d}{R}+\frac1X\right),\\
 \rho_H&=\frac{\sum_{p\in\mathcal P}p^{-2}}{\sum_{p\in\mathcal P}p^{-1}}
           \le\frac2H .
 \end{split}                                                     \tag{5}
\]
Thus one analytic target for the independent-residue average is a lower bound of order \((\log X)^{-2}\) for prime pairs with shifts \(2p\) averaged over the small primes \(p\) in \(\mathcal P\).

## Information revealed by divisibility

Write \(Y=(Y_p)_{p\in\mathcal P}\) and \(V=\prod_p V_p\). Define the array recording actual divisibility by the permitted prime by
\[
 B_{j,p}=\mathbf1_{A_{j,p}\ne0}\mathbf1_{Y_p=-j\bmod p},\qquad
 B_p=(B_{j,p})_{j=1}^H .
\]

**Theorem.** For \(G=1\) and fixed sufficiently small \(\eta>0\), the marked columns satisfy
\[
 \sum_{p\in\mathcal P} I(B_p;Y_p)\gg_\eta b\log H.               \tag{6}
\]
For every fixed \(\eta\) in the original range and any weight for which
\(\mathbb E\sum_{p,j}|A_{j,p}|^2/p\gg b\), the full permissive array satisfies
\[
 I(A;Y)\gg b\,\frac H{\log H}.                                 \tag{7}
\]
The hypothesis of (7) holds for \(G=1\) and sufficiently small fixed \(\eta\). It also holds for the prime-pair projection in (5) if that reference average has the proposed lower bound.

For (6), the dimension-two lower-bound sieve gives \(U(1)\gg_\eta(\log X)^{-2}\). The [lower bound for rough pairs](SQUAREFREE_TWIN_CORRELATION_NEAR_ZERO.md) applies on a fixed-ratio interval with this cutoff. Squares of primes \(q>y\) remove \(O(X/y+\sqrt X)\) pairs. By (3) for absolute values and the upper sieve, \(M=\mathbb E\sum_{p,j}B_{j,p}\asymp_\eta b\). Conditional on a fixed marked column, \(\sum_j B_{j,p}\mathbf1_{v=-j\bmod p}\) for uniform \(v\) is at most two and is nonzero with probability at most \(\sum_j B_{j,p}/p\). Under its true residue this sum is \(\sum_j B_{j,p}\). The relative-entropy inequality for the column therefore implies, after summation,
\[
 \sum_p\big(I(B_p;Y_p)+D_{\rm KL}(\mathcal L(Y_p)\Vert V_p)\big)
     \ \ge\ sM-(e^{2s}-1)(2/H)M .
\]
Take \(s=\tfrac14\log H\) and use the negligible marginal terms.

For (7), let \(K=\sum_{p,j}\mathbf1_{A_{j,p}\ne0}\) and \(L=|\mathcal P|\). Its mean is \(\gg H L/(\log X)^2\) by the hypothesis and \(p\asymp H\). For \(G=1\) and small fixed \(\eta\), the uniform lower bound at the shifts \(2p\) follows by applying the same dimension-two lower sieve on each interval \(pX<m\le2pX\) of length \(pX\). Each odd prime other than \(p\) removes two classes. Both conditions \(p^2\mid m\) and \(p^2\mid m+2p\) imply \(p\mid m\). The upper-bound sieve for the other odd primes on this progression gives \(O_\eta(XH/(p(\log X)^2))\). This is \(o(XH/(\log X)^2)\). Squares at \(q>y\) cost \(O(XH/y+\sqrt{XH})\). For the prime-pair projection the assertion follows directly from the average in (5) before the entropy comparison.

For two distinct entries of \(A\), their endpoint offsets are \(\{j,j+2p\}\) and \(\{i,i+2q\}\). There are \(O(HL^2)\) pairs sharing an endpoint and \(O(H^2L^2)\) pairs with four distinct endpoints. All offsets lie in \([1,3H]\), so any sieve prime greater than \(3H\) gives distinct roots for distinct endpoints. For \(k=3,4\) offsets, dropping primes up to \(3H\) and using the dimension-\(k\) upper sieve bounds the joint probability by \(O_\eta((\log H/\log X)^k)\). For example the same Selberg proof as for (4) works with \(k/q\) and the product is bounded using Mertens' estimate. Hence
\[
 \mathbb E K(K-1)\ll_\eta HL^2\left(\frac{\log H}{\log X}\right)^3+
                    H^2L^2\left(\frac{\log H}{\log X}\right)^4
             =o\left(\frac{HL}{(\log X)^2}\right).
\]
Here \(L\asymp H/\log H\) and \(H\le\sqrt{\log X}\). The first two inclusion–exclusion terms give \(\Pr(K>0)\gg HL/(\log X)^2\).

Given any array with \(K>0\), choose its first nonzero \((j,p)\) in a fixed ordering. Its true conditional law for \(Y\) excludes two residues \(-j,-j-2p\) at every \(q\in\mathcal P\setminus\{p\}\). The product-uniform probability of those simultaneous exclusions is \(\prod_{q\ne p}(1-2/q)\le \exp(-c/\log H)\). A probability law supported on an event of reference mass \(u\) has relative entropy at least \(\log(1/u)\). Average over \(A\) to obtain a lower bound \(\gg HL/((\log X)^2\log H)\) for \(I(A;Y)+D_{\rm KL}(\mathcal L(Y)\Vert V)\). By the Chinese remainder theorem, the marginal of \(Y\) corresponds to a residue modulo \(Q=\prod_{p\in\mathcal P}p=\exp(O(H))\). Its relative entropy against \(V\) is \(O((Q/(XH))^2)\) by the same residue mass estimate. Since \(b\asymp L/(\log X)^2\), this proves (7).

## Relation to Chowla and weighted dilation

For \(p\le y\), the original rough weight vanishes at \(pn\). Allowing \(p\) restores an exact multiplicative change of variables. Recording its occurrence defeats even the column-wise entropy budget by (6). Leaving \(p\) optional still makes the full-array budget fail under (7), because nonzero entries reveal restrictions at the other primes. Neither lower bound controls the column-wise information \(d\) of the permissive array.

Tao's logarithmically averaged two-point Chowla proof uses entropy decrement with genuinely multiplicative inputs and a dense error scale. Formula (2) gives the smaller information budget needed here. Lichtman's November 2026 preprint reports a binary Möbius bound \(O_\epsilon(X/(\log X)^{1-\epsilon})\) for a fixed nonzero shift. Its proof obtains much stronger savings on integers with factors in several subpower prime bands, and treats the contribution with both endpoints outside that part by an unsigned sieve. Power-rough pairs miss all those bands for every fixed \(\eta>0\). [Tao, The logarithmically averaged Chowla and Elliott conjectures for two-point correlations](https://arxiv.org/abs/1509.05422v4), [Lichtman, Improved bounds for the binary Chowla conjecture, Theorem 1.1 and §2](https://arxiv.org/pdf/2611.28392v1).

The [weighted-dilation audit](TWIN_CORRELATION_DOUBLING_AND_DILATION.md) identifies the actual marked prime factors required by the weighted-dilation theorem. If the marking primes are below \(y\), power-rough pairs are outside those marked sums. The variance estimate (4) uses only an unsigned sieve. The signed independent-residue average in (2) still needs an analytic bound.

**Open question.** Obtain the one-sided parity margin in (1) at scale \((\log X)^{-2}\). The column-wise route requires a bound \(d=o(b)\) for the permissive array and a lower bound of order \((\log X)^{-2}\) for the reference average in (5), or corresponding bounds directly for the signed and unsigned degree-three projection in (1).
