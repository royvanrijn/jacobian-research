# Semigroup-scale binary variance and the thin-window energy problem

This note proves a fixed-width variance theorem and exact first moments in every
fixed power window. It derives both from published results, without claiming
priority for the corollaries. The remaining power-window variance problem is
reduced to one sharp four-form energy estimate. The [general theory note](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md)
contains the every-target conjecture, the positive-proportion theorem and the
prime-based binary and ternary baselines.

Throughout, counts are unweighted and ordered. A statement about density one
means all but \(o(X)\) compatible targets in \((X,2X]\); it supplies neither
an effective cutoff nor every-target coverage.

## Definitions and the proved fixed-width theorem

Fix \(q\ge3\), and put
\[
L=\operatorname{lcm}(2,q),\quad
\alpha=1/\varphi(q)=1/\varphi(L),\quad
f=\mathbf1_{C_q},\quad \kappa=2-2\alpha.
\]
Let \(r_{q,H}(N)\) count ordered positive pairs \(a,b\in C_q\) with
\(a+b=N\) and \(|a-b|\le H\). Define the positive constant \(c_q\) by
\[
\sum_{n\le x}f(n)\sim c_qx(\log x)^{\alpha-1},\qquad
\mu_X=c_q(\log X)^{\alpha-1}.
\]
The mean follows from [Loughran--Matthiesen (LM), arXiv:1904.12845v2,
Lemma 2.8](https://arxiv.org/html/1904.12845v2#S2). An equivalent formula,
by Abelian partial summation of that mean, is
\[
c_q=\frac1{\Gamma(\alpha)}
\lim_{s\downarrow1}(s-1)^\alpha
\prod_{p\equiv1\ (q)}(1-p^{-s})^{-1}.
\]
Defining \(c_q\) by the mean avoids a normalization ambiguity in the source's
displayed Euler product. The Abelian formula includes \(1/\Gamma(\alpha)\).

Write \(\mathcal B=\{p:p\nmid L,\ p\not\equiv1\pmod q\}\). These primes
are odd. Put
\[
\sigma_p(N)=
\begin{cases}1+1/(p-1),&p\mid N,\\1-1/(p-1)^2,&p\nmid N,\end{cases}
\qquad
S(N)=L\mathbf1_{N\equiv2\ (L)}\prod_{p\in\mathcal B}\sigma_p(N).
\tag{V1}
\]
For fixed \(N\), the product converges: only finitely many factors are of
the first type, and the deviations of the others from one are summable.

**Fixed-width variance theorem Vq.** For every fixed \(0<\lambda\le1\),
with \(H=\lambda X\),
\[
\sum_{X<N\le2X}|r_{q,H}(N)-H\mu_X^2S(N)|^2
=o(XH^2\mu_X^4).\tag{V2}
\]
Incompatible targets contribute zero. This is a derived consequence of LM,
with the character contribution evaluated below; no novelty claim is made.
It gives a semigroup-scale lower bound for all but \(o(X)\) compatible
targets. It does not give a fixed-power sublinear window or an explicit rate.

## Local laws and the Loughran--Matthiesen constant

For \(z\) larger than every prime factor of \(L\), put
\[
W_z=L\prod_{p\in\mathcal B,p\le z}p,\qquad
w_z(t)=L\mathbf1_{t\equiv1\ (L)}
\prod_{p\in\mathcal B,p\le z}\frac{\mathbf1_{p\nmid t}}{1-1/p}.
\]
The mean of \(w_z\) over its period is one.
For a primitive Dirichlet character \(\chi\), average it over the subgroup of
units modulo \(\operatorname{lcm}(L,\operatorname{cond}\chi)\) reducing to
one modulo \(L\). The allowed primes occupy these classes equally, outside
finitely many primes. Character orthogonality gives
\[
m(f\chi)=\begin{cases}\alpha,&\chi\text{ factors through }L,\\0,&\text{otherwise}.
\end{cases}
\tag{V3}
\]
Thus LM's exceptional primitive characters are precisely those inducing the
\(\varphi(L)\) characters modulo \(L\). Every exceptional twist equals
\(f\) on positive integers. For fixed \(z\), the nonexceptional twists
have sums \(O_z(x/\log x)=o(x\mu_x)\), by LM Lemma 2.8. Since all prime
divisors of \(W_z\) are forbidden, orthogonality yields
\[
\sum_{\substack{n\le x\\n\equiv t\ (W_z)}}f(n)
=\left\{\frac{w_z(t)}{W_z}+o_z(1)\right\}x\mu_x.\tag{V4}
\]
Nonunit classes have zero mass. Subtraction proves the same law on fixed
macroscopic intervals, and rectangle approximation proves product-measure
convergence on bounded Jordan-measurable regions.

Here is the constant comparison needed for four-form correlations. Apply
[LM Theorem 1.9, Proposition 3.11 and Remark 3.12](https://arxiv.org/html/1904.12845v2#S3),
including equations (3.19)--(3.20), to \(r\) fixed linear polynomials
\(\psi_i\) with pairwise independent nonconstant parts, on a fixed convex
body of positive volume where all arguments are positive. Compare this
with \(r\) independent coordinate forms on a unit cube. For a fixed
sufficiently large \(B_0=z\), both applications have the common outer factor
\[
\left\{\frac1{\log Y}\prod_{z<p\le Y}(1+f(p)/p)\right\}^{r}.
\]
The finite coefficient in Remark 3.12 is \(A_zP_z(\boldsymbol\psi)\), where
\[
A_z=\left\{\frac{\varphi(L)e^{-\gamma\alpha}}{\Gamma(\alpha)}\right\}^{r}
\prod_{p\le z}(1-1/p)^{-r},
\]
and \(P_z\) is the Haar probability that every \(\psi_i\) is one modulo
\(L\), and no forbidden \(p\le z\) outside \(L\) divides any \(\psi_i\).
For the independent coordinates the probability is \(t_z^r\), with
\[
t_z=L^{-1}\prod_{p\in\mathcal B,p\le z}(1-1/p).
\]

To verify this specialization rather than assume an Euler product, use the
character sum in Remark 3.12. On the units modulo its \(Q_0\), it equals
\(\varphi(L)\mathbf1_{A_i\equiv1\ (L)}\). The modulus \(Q_0\) is divisible
by \(L\): each prime divisor already occurs, and a power \(p^k\Vert L\),
\(k\ge2\), admits a primitive character of that conductor, also for
\(p=2,k\ge2\). This covers the case where the least common multiple of
primitive conductors is smaller than \(L\).
The variable \(u_i\) in the source enumerates the complete \(p\le z\)-smooth
part of \(\psi_i\). A point away from the zero hyperplanes determines a
unique \(u_i\), then a unique unit \(A_i\) modulo \(Q_0\). The weight
\(f(u_i)\) forces all these small prime factors to be allowed; hence
\(u_i\equiv1\pmod L\), and \(A_i\equiv1\pmod L\) is exactly
\(\psi_i\equiv1\pmod L\). Summing the source's normalized congruence
probabilities gives \(P_z\). Zero hyperplanes have Haar measure zero, and
nonnegative summands justify the countable interchange.

Take \(Y\to\infty\) for each fixed \(z\), then let \(z\to\infty\).
The source's \(O(z^{-1/2})\) errors lie inside the brackets multiplying the
common outer factor. The reference coefficient \(A_zt_z^r\) is a fixed
positive factor times \(\prod_{p\le z,\ f(p)=1}(1-1/p)^{-r}\), so is bounded
below. Provided the local ratio has a finite limit, comparison with this
independent-coordinate reference proves that the arithmetic constant relative
to \(\mu_X^r\) is
\[
\lim_{z\to\infty}\frac{P_z(\boldsymbol\psi)}{t_z^r}
=\lim_{z\to\infty}\mathbf E_{\boldsymbol v}
 \prod_iw_z(\psi_i(\boldsymbol v)).\tag{V5}
\]
The independent reference has constant \(c_q^r\) by the long mean. A dilation
\(Y=2X\), or any sufficiently large fixed multiple, places each body within
the source's cube; \(\mu_Y/\mu_X\to1\). This argument resolves the character
issue for this specific semigroup; it makes no blanket local-to-global claim
for arbitrary Frobenian functions.

Take the four forms \((a,N-a,c,N-c)\) in \((N,a,c)\). Their coefficient
vectors are pairwise nonproportional. The body
\[
X<N\le2X,\qquad |2a-N|\le H,\quad |2c-N|\le H
\]
has volume \(XH^2\) and fixed shape when \(H=\lambda X\).
Conditionally on \(N\), the residues of \(a,c\) are independent. Define
\[
S_z(N)=\mathbf E_a w_z(a)w_z(N-a)
=L\mathbf1_{N\equiv2\ (L)}\prod_{p\in\mathcal B,p\le z}\sigma_p(N).
\]
Indeed the number of forbidden residues for \(a\) is one when \(p\mid N\)
and two otherwise. With
\(c_p(N)=p-1\) for \(p\mid N\), and \(-1\) otherwise,
\(\sigma_p=1+c_p/(p-1)^2\),
\(\mathbf E_pc_p=0\) and \(\mathbf E_pc_p^2=p-1\). Thus
\[
K_z:=\mathbf E_NS_z(N)^2
=L\prod_{p\in\mathcal B,p\le z}(1+(p-1)^{-3})
\longrightarrow K_q<\infty.\tag{V6}
\]
Equations (V5)--(V6) give
\[
\sum_{X<N\le2X}r_{q,H}(N)^2
=\{K_q+o(1)\}XH^2\mu_X^4.\tag{V7}
\]
For fixed \(z\), apply (V4) to the two independent variables \((a,b)\)
in the polygon \(X<a+b\le2X,\ |a-b|\le H\). Its area is \(XH\).
The mixed arithmetic mean is \(\mathbf E_{a,b}w_z(a)w_z(b)S_z(a+b)=K_z\).
Together with periodicity of \(S_z^2\), this proves
\[
\lim_{X\to\infty}\frac1X\sum_{X<N\le2X}
\left|\frac{r_{q,H}(N)}{H\mu_X^2}-S_z(N)\right|^2=K_q-K_z.\tag{V8}
\]

## Arithmetic mean-square convergence and the fixed-width conclusion

For completeness the passage to the infinite product is along actual integers.
Let \(A=\prod_{p\in\mathcal B}(1-(p-1)^{-2})>0\), and define a nonnegative
squarefree multiplicative function \(h\) by \(h(p)=1/(p-2)\) on
\(\mathcal B\), zero at other primes. Then
\[
S(N)=LA\mathbf1_{N\equiv2\ (L)}\sum_{d\mid N}h(d).
\]
Its truncated version uses \(A_z\) and \(h\) restricted to primes at most
\(z\). The square expansion is dominated by
\[
\sum_{d,e}\frac{h(d)h(e)}{[d,e]}
=\prod_{p\in\mathcal B}\left\{1+\frac{2h(p)+h(p)^2}{p}\right\}<\infty.
\]
The fraction, relative to \(X\), of \(N\in(X,2X]\) divisible by \([d,e]\)
is at most \(2/[d,e]\). Dominated convergence and the Chinese remainder
theorem prove
\[
\lim_{z\to\infty}\limsup_{X\to\infty}
\frac1X\sum_{X<N\le2X}|S(N)-S_z(N)|^2=0.\tag{V9}
\]
A uniform \(O_q(1/z)\) estimate is proved again below when needed for the
power-window reduction. Applying the normalized triangle inequality to
(V8)--(V9) proves Vq.

A diagonal consequence gives a sublinear window without an explicit rate:
choose \(X_j\) so that (V2) at \(H=X/j\) has normalized variance at most
\(1/j\) for all \(X\ge X_j\), and \(\log j/\log X_j\to0\). The piecewise
window \(H=X/j\) for \(X_j\le X<X_{j+1}\) is \(X^{1-o(1)}=o(X)\).
This specifies no fixed power and no logarithmic rate.

## Exact first moments in every power window

Let \(f(n)\) indicate positive integers all of whose prime factors belong to a fixed Frobenian set \(\mathcal P\) of prime density \(0<\alpha\le1\). Write
\[
m_X=\frac2X\sum_{X/2<n\le X} f(n).
\]
**Power-window first-moment theorem Aq.** For fixed \(0<\theta<1\) and \(H=X^\theta\), some \(\eta>0\) satisfies
\[
\sum_{X<N\le 2X} r_{\mathcal P,H}(N)
=XH m_X^2(1+O_{\mathcal P,\theta}((\log X)^{-\eta})).
\tag{1}
\]
In particular the constant is the square of the one-variable leading coefficient in \(m_X\sim c_{\mathcal P}(\log X)^{\alpha-1}\).

**One-variable input.** Apply [Mangerel, Research in the Mathematical Sciences 10 (2023), Theorem 1.7, equations (1.6)–(1.11)](https://link.springer.com/article/10.1007/s40687-023-00376-0). Use \(A=\alpha,B=C=1\): \(|f(p)|\le 1,\ |f(n)|\le d_1(n)^1=1\). Taking endpoint differences in [Loughran–Matthiesen, Lemma 2.4(2)](https://arxiv.org/html/1904.12845v2#S2) for \(f\) and the constant function gives \(\sum_{z<p\le w}f(p)/p=\alpha\sum_{z<p\le w}1/p+O_{\mathcal P}(1/\log z)\) uniformly for \(2\le z\le w\). This proves (1.6) in Mangerel. Its inclusion after (1.9) gives \(\gamma=1\) and (1.8) for any fixed \(0<\sigma<\frac{\alpha}{3}\{1-\mathrm{sinc}(\pi\alpha/2)\}\). We may take \(t_0=0\) because \(f\) is nonnegative. Mangerel's products satisfy
\[
D_X=\prod_{p\notin\mathcal P,p\le X}(1+1/p)\asymp(\log X)^{1-\alpha},
\qquad
P_X=\prod_{p\notin\mathcal P,p\le X}(1-1/p)\asymp(\log X)^{\alpha-1}.
\]
Both estimates follow by logarithms from the preceding prime sums. The dyadic one-variable asymptotic in LM Lemma 2.8 gives \(P_X\asymp m_X\). Here \(c_{\mathcal P}\) is defined by the mean. If its Euler product is used, include the Selberg–Delange factor \(1/\Gamma(\alpha)\).

Put \(h_0=h/D_X\). Uniformly for \(H/2\le h\le 2H\), \(10\le h_0\le X/(10D_X)\) and \(\log h_0\asymp_\theta\log X\) once \(X\) is large. The normalized mean-square error in Theorem 1.7 is bounded by
\[
P_X^2\left\{\left(\frac{\log\log h_0}{\log h_0}\right)^\alpha+
\left(\frac{\log\log X}{(\log X)^{\kappa}}\right)^\alpha\right\},
\qquad \kappa=\min(1,\sigma)/29.
\]
Thus, for some \(d>0\), with \(U_h(x)=h^{-1}\sum_{x-h<n\le x}f(n)\),
\[
\frac2X\int_{X/2}^X |U_h(x)-m_X|^2\,dx\ll m_X^2(\log X)^{-d}.
\tag{2}
\]
For example \(d=\alpha\kappa/2\) works for large \(X\). All constants are permitted to depend on the fixed data.

**From a local mean to the sharp binary first moment.** Let \(T(h)=\sum_{X<a+b\le2X}f(a)f(b)(h-|a-b|)_+\). Integrating \((hU_h(x))^2\) counts each ordered pair with weight
\[
|[X/2,X]\cap[\max(a,b),\min(a,b)+h)|.
\]
For \(u=|a-b|<h\) and \(m=(a+b)/2\), the second interval is \([m+u/2,m+h-u/2)\subset[m,m+h]\), and its length is \(h-u\). The weight differs from \(\mathbf1_{X/2<m\le X}(h-u)\) only for \(m\) within \(h\) of one of the two endpoints. Each discrepant pair has both coordinates within \(2h\) of that endpoint. There are \(O((h+1)^2)\) such pairs, and each changes the total by at most \(h\), since \(0\le f\le1\). Hence
\[
T(h)=\int_{X/2}^X(hU_h(x))^2\,dx +O(h(h+1)^2).
\]
By (2) and Cauchy–Schwarz the integral equals \(\frac X2 h^2m_X^2(1+O((\log X)^{-d/2}))\). The relative geometric error is \(O(H/(Xm_X^2))\), smaller than any fixed inverse log power. This estimate is uniform for all three widths used next. For \(R(H)=\sum_{X<N\le2X}r_{\mathcal P,H}(N)\) and \(0<\epsilon<1/2\),
\[
\frac{T(H)-T((1-\epsilon)H)}{\epsilon H}\le R(H)\le
\frac{T((1+\epsilon)H)-T(H)}{\epsilon H}.
\]
Indeed each summand is the nonnegative weight \(f(a)f(b)\) times a positive part, whose right slope at \(H=u\) counts the equality \(u=H\). The main slopes are \(XH m_X^2(1\pm\epsilon/2)\), and the relative error is \(O((\log X)^{-d/2}/\epsilon)\). Taking \(\epsilon=(\log X)^{-d/4}\) proves (1).

## Fixed arithmetic projections in every power window

**Projection theorem Pq.** Now \(f=\mathbf1_{C_q}\), \(L=\mathrm{lcm}(2,q)\), \(\alpha=1/\varphi(q)\). Let \(W\) be a fixed squarefree product of primes \(p\nmid L\) with \(p\not\equiv1\bmod q\), allowing \(W=1\). Every supported integer is a unit modulo \(W\). For \(u,v\in(\mathbf Z/W)^*\), one has, for some \(e_W>0\),
\[
\sum_{\substack{X<a+b\le2X,\ |a-b|\le H\\ a\equiv u,\ b\equiv v\ (W)}} f(a)f(b)
=\frac{XH m_X^2}{\varphi(W)^2}\{1+O_{q,W,\theta}((\log X)^{-e_W})\}.
\tag{3}
\]

To verify (3), put \(g_\chi=f\chi\) for characters modulo \(W\). Since \(|g_\chi|=f\), it satisfies Mangerel's hypotheses with the same \(A,B,C,\gamma,\sigma\), and the scale and product are still \(D_X,P_X\). For \(\chi\) principal, (2) applies. For \(\chi\ne\chi_0\) of order \(k\), we claim that for \(Y=X/2,X\) and \(|s|\le 2X\),
\[
E_\chi(s;Y):=\sum_{\substack{p\le Y\\p\equiv1\ (q)}}\frac{1-\Re(\chi(p)p^{-is})}{p}
\ge c_W\log\log X-O_{q,W}(1).
\tag{4}
\]
For \(|s|\le1\), fixed-modulus Chebotarev (LM Lemma 2.4(1) for each reduced residue class mod \(qW\)) and CRT give \(\sum_{p\le y,\,p=1\ (q)}\chi(p)=O_{q,W}(y e^{-c_{q,W}\sqrt{\log y}})\), since the main coefficients sum to zero by character orthogonality. Partial summation against \(y^{-1-is}\), whose derivative is at most \(2y^{-2}\), makes the twisted prime sum with weight \(p^{-1-is}\) bounded absolutely in size, uniformly for \(|s|\le1\). The untwisted prime mass is \(\alpha\log\log X+O(1)\). For \(1\le |s|\le2X\), \(|1-z^k|\le k|1-z|\) on the unit circle yields \(E_\chi(s;Y)\ge k^{-2}E_{\chi_0}(ks;Y)\). Use (1.8) of Mangerel for \(f\) on the prime range up to \(4kX\), for which the minimizer is zero. The result is at least \(k^{-2}\sigma\log\log X-O_{q,W}(1)\): \(|ks|\le 2kX\), and removing the primes \(Y<p\le 4kX\) costs at most \(2\sum_{Y<p\le4kX}1/p=O_{k}(1)\). This proves (4).

Let \(t_0\in[-X,X]\) be any minimizer for \(g_\chi\) in Mangerel's theorem. Apply the [Halász–Montgomery–Tenenbaum theorem, Granville–Soundararajan, p. 1192](https://doi.org/10.4153/CJM-2003-047-0) to the 1-bounded function \(g_\chi(n)n^{-it_0}\), on each of \(Y=X/2,X\), with \(T=(\log X)^4\). Its full distance, at any \(|v|\le2T\), is
\[
\sum_{p\le Y}\frac{1-\Re(g_\chi(p)p^{-i(t_0+v)})}{p}
=(1-\alpha)\log\log X+E_\chi(t_0+v;Y)+O(1).
\]
For large \(X\), \(|t_0+v|\le2X\). Thus the minimum \(M\) of the full distance satisfies \(M\ge(1-\alpha+c_W)\log\log X-O(1)\). The bound \((1+M)e^{-M}+T^{-1/2}\) for each cumulative average is \(O_{q,W}(m_X(\log X)^{-c'_W})\) for some \(c'_W>0\); we can reduce \(c_W\) so that \(c_W\le\alpha\). Taking the two endpoints gives the same bound on Mangerel's dyadic twisted mean. Its factor \(h^{-1}\int_{x-h}^x t^{it_0}dt\) has modulus at most 1. Theorem 1.7 therefore bounds the normalized local mean of every nonprincipal \(g_\chi\) in \(L^2(dx/(X/2))\) by \(O_{q,W}(m_X(\log X)^{-d_W})\) for some \(d_W>0\), uniformly for \(H/2\le h\le2H\). Character orthogonality gives the same error for the local mean of \(f(n)\mathbf1_{n=u\ (W)}\) around \(m_X/\varphi(W)\). Testing the product of the two local means for \(u,v\) by Cauchy–Schwarz, using exactly the Fubini and nonnegative finite-difference argument proving (1), proves (3).

Consequently for any residue \(t\bmod W\) (including nonunits), write
\[
s_W(t)=\frac{W}{\varphi(W)^2}\#\{(u,v)\in ((\mathbf Z/W)^*)^2:u+v=t\}.
\]
Then
\[
\sum_{\substack{X<N\le2X\\ N\equiv t\ (W)}}r_{q,H}(N)=\frac{XH m_X^2}{W}s_W(t)(1+O_{q,W,\theta}((\log X)^{-e_W})).
\tag{5}
\]
The constant for each \(W\) is fixed before taking \(X\to\infty\). All representations occurring already have \(N\equiv2\bmod L\). On the latter progression the classwise mean of \(r/(H m_X^2)\) tends to \(L s_W(t)\). Indeed there are \(X/(LW)+O(1)\) targets in each compatible class.

Here \(s_W(t)=\prod_{p\mid W}(1+c_p(t)/(p-1)^2)\), where \(c_p(t)=p-1\) for \(p\mid t\) and \(-1\) otherwise. This follows from CRT and the counts \(p-1,p-2\) of nonzero ordered summands modulo \(p\). The published input gives no uniform savings in \(W\). All fixed projections are sufficient to identify the mixed moment with the infinite singular factor in conjunction with the bounded second moment F8 and a mean-square tail. In particular, for each fixed \(W\), conditional variance gives
\[
\liminf_{X\to\infty}\frac1X\sum_{X<N\le2X}\left(\frac{r_{q,H}(N)}{H m_X^2}\right)^2
\ \ge L\prod_{p\mid W}\left(1+\frac1{(p-1)^3}\right).
\]
The equality of the product to the normalized second moment of the projection uses CRT, \(\mathbf E c_p=0\), and \(\mathbf E c_p^2=p-1\) modulo \(p\). Increasing \(W\) over the forbidden coprime primes gives the infinite product lower bound. Equality in the limit requires control of the four-form additive energy in the thin region.


## The exact power-window variance reduction

**Reduction theorem Eq.** Let \(\mathcal B=\{p:p\nmid L,\ p\not\equiv1\bmod q\}\), let \(W_z=\prod_{p\in\mathcal B,p\le z}p\), and write \(S_z(N)=L\mathbf1_{N\equiv2\ (L)}s_{W_z}(N)\). Define \(S\) by taking the infinite product in the same formula. Every \(p\in\mathcal B\) is odd. Since
\[
s_{W_z}(N)=A_z B_z(N),\quad
 A_z=\prod_{\substack{p\in\mathcal B\\p\le z}}\left(1-\frac1{(p-1)^2}\right),
\quad B_z(N)=\prod_{\substack{p\in\mathcal B,\ p\le z\\p\mid N}}\left(1+\frac1{p-2}\right),
\]
the pointwise limiting product is \(A B(N)\), where \(A=\lim A_z>0\) and the divisor product \(B(N)\) has finitely many factors. For \(z\ge3\) we have \(|A_z-A|\ll 1/z\) by summing \((p-1)^{-2}\) over all integers beyond \(z\). On squarefree integers supported on \(\mathcal B\) define multiplicative \(b_2\) by \(b_2(p)=2/(p-2)+1/(p-2)^2\), \(b_2(1)=1\), with value zero off that support. The product \(\sum_d b_2(d)/d=\prod_{p\in\mathcal B}(1+b_2(p)/p)\) is finite since \(b_2(p)\ll 1/p\). For \(\mathbf E_X=\frac1X\sum_{X<N\le2X}\) we have uniformly for \(X\ge1\)
\[
\mathbf E_X B^2\le 2\sum_d b_2(d)/d,\qquad
\mathbf E_X(B-B_z)^2\le \mathbf E_X(B^2-B_z^2)
\le 2\sum_{P^+(d)>z}\frac{b_2(d)}d\ll \sum_{p>z}p^{-2}\ll z^{-1}.
\]
The first two upper bounds use \(\#\{N\le2X:d\mid N\}\le 2X/d\); the tail estimate follows by selecting a prime \(p>z\) of \(d\) and factoring out \(b_2(p)/p\). In particular
\[
\sup_{X\ge1}\mathbf E_X|S-S_z|^2\ll_q z^{-1}.
\tag{6}
\]
This proves actual arithmetic \(L^2\) approximation along the dyadic ranges.

Put \(R_X(N)=r_{q,H}(N)/(H m_X^2)\) and \(K_z=L\prod_{p\mid W_z}(1+(p-1)^{-3})\). For fixed \(z\), (5) applied to each \(t\bmod W_z\) gives
\[
\lim_{X\to\infty}\mathbf E_X R_X S_z
=\frac L{W_z}\sum_{t\bmod W_z}s_{W_z}(t)^2=K_z.
\]
Periodicity modulo \(LW_z\) gives \(\lim_{X\to\infty}\mathbf E_X S_z^2=K_z\) with precisely the same normalizing \(L\). (F8) in the general theory note and the one-variable mean imply \(\sup_{X\gg1}\mathbf E_X R_X^2<\infty\); the divisor estimate proves the corresponding uniform second-moment bound on \(S\) and \(S_z\). Cauchy–Schwarz in (6) now permits passage \(X\to\infty\) then \(z\to\infty\) in both mixed and model moments. Thus for \(K=L\prod_{p\in\mathcal B}(1+(p-1)^{-3})\),
\[
\mathbf E_X R_X S=K+o(1),\qquad
\mathbf E_X S^2=K+o(1),\qquad
\mathbf E_X|R_X-S|^2=\mathbf E_X R_X^2-K+o(1).
\tag{7}
\]
Equivalently \(\sum|r_{q,H}-H m_X^2S|^2=\sum r_{q,H}^2-(K+o(1))XH^2m_X^4\). All sums may be taken over the compatible progression since both counts and models vanish outside it. Positivity gives \(\liminf\mathbf E_X R_X^2\ge K\); the complete variance follows if and only if \(\limsup\mathbf E_X R_X^2\le K\). Since \(m_X/(c_q(\log X)^{\alpha-1})=1+O_q(1/\log X)\) and \(\mathbf E_X S^2=O(1)\), this criterion is identical upon replacing \(m_X\) in the main term by \(c_q(\log X)^{\alpha-1}\).


## What the variance results establish

For the fixed-width theorem, \(\mathfrak S_q(N)\) denotes the function called
\(S(N)\) above. Its target-independent lower bound on compatible targets is
\[
\mathfrak S_q(N)\ge
L\prod_{p\in\mathcal B}\left(1-\frac1{(p-1)^2}\right)>0.
\]
Consequently Vq gives \(r_{q,\lambda X}(N)\gg_q X/(\log X)^{\kappa_q}\)
for all but \(o(X)\) compatible targets. More precisely, some
\(\varepsilon(X)\to0\) gives
\[
r_{q,\lambda X}(N)
=c_q^2\lambda X(\log X)^{-\kappa_q}\mathfrak S_q(N)
\{1+O(\varepsilon(X))\}
\]
uniformly outside \(o(X)\) targets. Chebyshev's inequality and the displayed
positive lower bound prove this assertion; choosing a decreasing error threshold
requires no rate from Loughran--Matthiesen. At \(q=12\), these are
\(\gg X/(\log X)^{3/2}\) semigroup pairs, rather than the
\(\gg X/(\log X)^2\) prime-pair baseline. The choice \(\lambda=1/3\)
places both summands in \([N/3,2N/3]\), since \(X/3\le N/3\).

For every fixed \(0<\theta<1\), Aq and Pq determine the first moment and every
fixed forbidden-prime projection at the semigroup scale. Eq proves that the
requested variance theorem is equivalent to the sharp upper energy bound
\[
\begin{split}
&\sum_{X<N\le2X}
\sum_{\substack{|2a-N|\le H\\|2c-N|\le H}}
 f(a)f(N-a)f(c)f(N-c)\\
&\qquad\le\{\mathcal K_q+o(1)\}XH^2m_X^4,
\quad H=X^\theta,\qquad
\mathcal K_q=L\prod_{p\in\mathcal B}\left(1+\frac1{(p-1)^3}\right).
\end{split}\tag{E}
\]
The inner sums range over positive summands; extending \(f\) by zero to
nonpositive integers gives the same formula. The matching lower bound is already
proved by Pq. The current upper bound (F8) supplies only an unspecified constant,
which cannot replace the exact constant \(\mathcal K_q\) in (E).

The fixed-width proof uses a fixed positive-volume convex body. At
\(H=X^\theta\), its volume after scaling by \(X\) tends to zero. The stated
Loughran--Matthiesen theorem gives no error uniform in that shrinking body.
Thus Vq cannot simply be applied with \(\lambda=X^{\theta-1}\).
The Mangerel application supplies a logarithmic saving in Aq and Pq for each
fixed projection; it supplies no thin four-form energy asymptotic. No
power-window density-one theorem, power-saving exceptional set or every-target
binary theorem is concluded here.

A power-saving centered variance, at the same main-term scale, would bound the
number of zero or very small representation counts by Chebyshev. Eliminating
all exceptions still requires a separate uniform argument. The every-target
conjecture FBS-U and the effective completion problem remain open.

## Proof and review scope

The mathematical inputs used here were inspected in the identified primary
sources: Loughran--Matthiesen arXiv:1904.12845v2, Mangerel's published 2023
Theorem 1.7, and the Halász--Montgomery--Tenenbaum estimate cited in the
projection proof. The conductor calculation, local constants, integer
mean-square tails, interval geometry and short-interval hypotheses have
received internal AI review with separate reconstructions of the two main
arguments. Small exact residue enumerations at odd primes through 19 and
modulus 12 agreed with the displayed two- and four-form factors. Those checks
are diagnostics; the identities are proved algebraically above.

There is no external human review, formal verification, new finite-range
coverage certificate or originality claim for these corollaries. The full
Iwaniec--Pomykała priority audit remains unresolved in the case study.
