# An explicit small-ratio region for ordinary binomial Strong Factorial

For integers \(a\geq1,d\geq2\), write \(D=a+d\) and \(J=d^2/D\). The actual reversed moment polynomials are
\[
 A_m(x)=\sum_{k=0}^m\binom mk(am+1)_{dk}x^{m-k},\qquad m\geq1.
\]

**Theorem EffectiveSmallRatioOrdinaryBinomialSFC.** If \(J\geq128\), then \(A_{N-1},A_N\) are coprime over \(\mathbb Q\) for every \(N\geq2\). The associated ordinary binomial Strong Factorial statement holds over every characteristic-zero field and at every positive starting index.

**Corollary.** The same Strong Factorial conclusion holds if \(a/(d-1)^2\leq1/512\). For \(d<3a\) this follows from \(J>d^2/(4a)>(d-1)^2/(4a)\). For \(d\geq3a\), use the [relative-gap theorem](RELATIVE_GAP_THREE_STRONG_FACTORIAL.md) when \(a\geq2\) and the [lower-one theorem](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md) when \(a=1\).

## 1. Exact normalization

Theorem \(\mathrm{thm:sparse}\) in the [initial and boundary argument](lower_exponent_initial_and_boundary.tex) excludes a common root at terminal \(N\) if \(J\geq4\log(8N)\). This covers \(N\leq10^8\), since \(8N<2^{30}\) and \(4\log(8N)<120\). Seek a common root \(x\) with \(N>10^8\) and \(J<4\log(8N)\). Write \(R=D/a\). If \(R\geq2\), then \(d\geq D/2\) and \(R\leq D\leq4J\). The bound \(R\leq4J\) also holds for \(R<2\). At each of \(m=N-1,N\) we have
\[
 m\geq10^8,\qquad R\leq4J<16\log(8N)<128N\leq256m.                \tag{1}
\]

For \(\ell\geq1\) set \(f_\ell(S)=(S/\ell)^\ell e^{\ell-S}\). Put
\[
 W=\exp(D\log D-a\log a-d),\qquad R_0=D^D e^{-D},\qquad
 S_0=W^{1/d},\qquad q_m=x/(m^d W).
\]
The factorial integral followed by the change of variables \(X=mS\) gives
\[
 \frac{(am)!}{m^{Dm+1}}A_m(m^dWq)
 =\int_0^\infty[S^a(Wq+S^d)e^{-S}]^m\,dS
 =R_0^m\int_0^\infty(q f_a+f_D)^m\,dS.                              \tag{2}
\]
Both arguments \(q_m\) are nonzero by the constant terms of \(A_m\).

The chord and midpoint inequalities for \(\log\) on \([a,D]\), and its average on \([0,D]\), give
\[
 \log S_0=\frac1d\int_a^D\log u\,du,\qquad
 \sqrt{aD}\leq S_0\leq(a+D)/2,\qquad S_0\geq D/e.
\]
The ratio \(f_D/f_a\) equals \((S/S_0)^d\). The elementary inequalities
\[
 \log(1+t)-t\leq-t^2/2\quad(-1<t\leq0),\qquad
 \log(1+t)-t\leq-t^2/[2(1+t)]\quad(t\geq0)
\]
give
\[
 f_D(S_0)\leq e^{-J/8},\qquad f_a(S_0)\leq e^{-J/18}.               \tag{3}
\]
For the first, \(D-S_0\geq d/2\). For the second, view the identric mean as \(I(t)=\exp(d^{-1}\int_0^d\log(t+u)\,du)\) for \(t>0\). Then \(I'(t)=I(t)d^{-1}\int_0^d(t+u)^{-1}\,du\geq1\), since the arithmetic mean of the reciprocal is at least its geometric mean by Jensen's inequality. As \(I(0)=d/e\) by the right limit, \(S_0-a\geq d/e>d/3\), giving the second bound in (3). Thus \(f_D\) on the left half of the split at \(S_0\) and \(f_a\) on the right half are each at most \(1/900\). Indeed \(J/18>7\) and \(e^7>(8/3)^7>900\).

Also \(\int_0^\infty f_\ell(S)\,dS\leq6\sqrt\ell\). On \(0<S\leq2\ell\), use \(v=(S-\ell)/\sqrt\ell\) and the bound \(e^{-v^2/4}\). The remaining integral is at most \(2\), since \((\log f_\ell)'(S)\leq-1/2\) there. The windows \(|S-\ell|\leq\sqrt\ell/8\) at \(\ell=a,D\) lie in their respective halves, since
\[
 (S_0-a)/\sqrt a\geq\sqrt D-\sqrt a\geq\sqrt J/2,\qquad
 (D-S_0)/\sqrt D\geq\sqrt J/2.
\]

## 2. Two local stationary points

Near \(S=\ell+\sqrt\ell\,v\), use the normalized function \(g(S)=f_\ell(S)+c f_{\ell'}(S)\), where \(\{\ell,\ell'\}=\{a,D\}\) and \(|c|\leq3\). On \(|v|\leq1/4\), both cross ratios \(h(v)=f_{\ell'}(S)/f_\ell(S)\) satisfy \(|h(v)|\leq e^{-J/4}\). To verify this, set \(K_L=\log(W/a^d)=D\log(D/a)-d\) and \(K_R=\log(D^d/W)=d-a\log(D/a)\). Both are at least \(J/2\), by integrating \((D-u)/u\) and \((u-a)/u\) over \([a,D]\). The logarithm series gives \(|\log(1+v/\sqrt\ell)|\leq2|v|/\sqrt\ell\). The right cross ratio is at most \(e^{-K_R+\sqrt J/2}\). For \(R\leq4\), the left ratio is at most \(e^{-K_L+\sqrt{JR}/2}\leq e^{-J/2+\sqrt J}\). For \(R\geq4\), the function \(u\log u/(u-1)-1\) increases for \(u>1\), since \(u-1>\log u\). At \(u=4\) its value is at least \(7/9\), using \(\log 4\geq4/3\) from convexity of \(1/u\) on \([1,2]\). Hence \(K_L\geq3d/4\) in this case and the left ratio is at most \(e^{-d/4}\leq e^{-J/4}\).

Put \(\epsilon=10^{-12}\), \(\eta=10^{-8}\) and \(b=1/\sqrt\ell\). A local logarithm of \(g\) is
\[
 F(v)=H_b(v)+L(v),\quad H_b(v)=\frac{\log(1+bv)-bv}{b^2},\quad
 L(v)=\log(1+c h(v)).
\]
The power-series logarithms are defined on \(|v|\leq1/4\) and \(|L|<6e^{-J/4}<\epsilon\), using \(e^{J/4}\geq e^{32}>(8/3)^{32}>50^8\). On \(|v|\leq1/6\), the Cauchy bounds for \(L',L'',L'''\) are \(12\epsilon,288\epsilon,10368\epsilon\), respectively. Since \(H_b'(v)=-v/(1+bv)\), Rouché's theorem on \(|v|=\eta\) gives exactly one simple zero \(v_*\) of \(F'\) inside. The analytic implicit function theorem gives local holomorphic dependence on \(c\). Moreover
\[
 |F(v_*)|<2\epsilon,\qquad
 \alpha=-F''(v_*)/2,\quad |\alpha-1/2|<10^{-6},\qquad
 |F'''(v)|<6\ \ (|v|\leq1/6).                                      \tag{4}
\]
These use the series for \(H_b\) and
\(H_b''=-(1+bv)^{-2}\), \(H_b'''=2b(1+bv)^{-3}\).
Write \(\lambda=\exp(F(v_*))\).

## 3. Explicit estimates for the two integrals

Integrate \(g^m\) on either the half containing \(\ell\), for \(|c|\leq3\), or on the full positive ray, for \(|c|\leq1/2\). Subject to (1), the integral in either case has the form
\[
 \sqrt\ell\,\lambda^m\sqrt{2\pi/m}\,(1+E),\qquad |E|<1/100.           \tag{5}
\]
To prove this, put \(r=1/8\). On the selected half outside \(|v|\leq r\), the selected \(f_\ell\) is at most \(e^{-1/144}\) by the logarithm bounds and monotonicity. The other term contributes at most \(1/300\) by (3). In the full-ray case the bound on the opposite half is \(1/2+1/900\). Both bounds are at most \(|\lambda|e^{-1/500}\), because \(e^{-1/144}\leq1-1/145\) and \(|\lambda|e^{-1/500}\geq1-1/500-2\epsilon\). The full \(L^1\) norm of \(g\) is at most \(24\sqrt D\). Removing one power from \(|g|^m\), and using (1) and \(1/|\lambda|<2\), bounds the relative contribution from outside the window by
\[
 48\sqrt{mR}\,e^{-(m-1)/500}\leq768m e^{-(m-1)/500}<10^{-6}.
\]
For the last bound, use \(m-1\geq m/2\) and \(e^t\geq t^4/24\) with \(t=(m-1)/500\).

Deform the segment \([-r,r]\) in the \(v\)-plane to the translated segment \(v_*+[-r,r]\) and the two straight connectors. This preserves the integral of the entire function \(g^m\). The connectors have length at most \(\eta\) each, stay within \(\eta\) of \(\pm r\), and satisfy \(\Re(F(v)-F(v_*))\leq-1/500\). Indeed \(|F'|<1\) on \(|v|\leq1/6\) and \(\Re F(\pm r)\leq-1/144+\epsilon\). Their total relative contribution is at most \(2\eta\sqrt m e^{-m/500}<10^{-6}\).

On the translated segment write \(F(v_*+t)-F(v_*)=-\alpha t^2+\mathcal R(t)\). Taylor's remainder and (4) give \(|\mathcal R(t)|\leq |t|^3\). After factoring \(\lambda^m\), the absolute error from replacing the exponential by the quadratic exponential is at most
\[
 m\int_{\mathbb R}|t|^3e^{-mt^2/3}\,dt=9/m
\]
since \(\Re\alpha-r>1/3\). The two Gaussian tails in extending \([-r,r]\) to \(\mathbb R\) have absolute integral at most \(3/(mr)=24/m\), by \(\Re\alpha>1/3\) and integration of \(t e^{-mt^2/3}\) on \(t\geq r\). The full Gaussian is \(\sqrt{2\pi/m}(2\alpha)^{-1/2}\), from the positive-real formula by holomorphic continuation in \(\Re\alpha>0\). Its normalized factor differs from \(1\) by less than \(10^{-5}\), from (4) by integrating the derivative on the straight segment from \(1/2\) to \(\alpha\). Summing the relative bounds gives \(2\cdot10^{-6}+33/\sqrt m+10^{-5}<1/100\), proving (5).

## 4. Consecutive orders

For \(|q|\leq1/2\), the full-ray case of (5) at \(\ell=D\) applies to \(qf_a+f_D\). For \(|q|\geq2\), divide by \(q\) and apply (5) to \(f_a+q^{-1}f_D\) at \(\ell=a\). Thus both arguments of a putative common root in (2) satisfy \(1/2<|q_{N-1}|,|q_N|<2\).

For \(1/2<|q|<2\), split the integral at \(S_0\) and apply (5) on each half, including the factor \(q^m\) on the left. Let \(\lambda_L=q\exp(F_L(v_{*,L}))\) and \(\lambda_R=\exp(F_R(v_{*,R}))\) be the two critical values of \(qf_a+f_D\). At either zero integral the two partial integrals have equal modulus, so
\[
 \left|m\log|\lambda_R/\lambda_L|+\tfrac12\log R\right|
 \leq\log(101/99)<1/40.                                             \tag{6}
\]
Set \(m=N-\theta\), \(q_\theta=(N/m)^d q_N\) for \(0\leq\theta\leq1\), and
\[
 E(\theta)=m\log|\lambda_R(q_\theta)/\lambda_L(q_\theta)|
                  +\tfrac12\log R .
\]
The interpolation points are on the same radial line and their moduli are between the endpoint moduli. This interpolation concerns critical values. The integrals were used only at the two integer orders.

For \(w=Wq\), the phase \(\Phi(S,w)=a\log S+\log(w+S^d)-S\) obeys \(dw\Phi_w=D-S\) at a stationary point. Consequently, for the stationary points \(S_L,S_R\) and local logarithm branches,
\[
 d q\,\frac{d}{dq}\log(\lambda_R/\lambda_L)=-(S_R-S_L).
\]
By (4), \(\left|\log|\lambda_R/\lambda_L|\right|\leq\log2+4\epsilon<1\). Also \(S_L=a+\sqrt a\,v_{*,L}\) and \(S_R=D+\sqrt D\,v_{*,R}\) give \(\Re(S_R-S_L)\geq d-2\eta\sqrt D>0.99d\). Therefore
\[
 E'(\theta)=-\log|\lambda_R/\lambda_L|-\Re(S_R-S_L)
 <1-0.99d\leq-0.98.
\]
This contradicts (6) at the two endpoints.

For a nonzero binomial \(X^a(A+B X^d)\) with \(B\ne0\), its factorial moment of order \(m\) is \(B^m(am)!\,A_m(A/B)\). The rational Bezout identity persists over every characteristic-zero coefficient field. A zero coefficient leaves a monomial with nonzero positive moments. This proves the stated ordinary Strong Factorial assertion in the explicit region.

The proof quantifies the two-well argument of the [qualitative small-ratio theorem](UNRESTRICTED_SMALL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md). The other qualitative sectors of the support-uniform cutoff require separate effective estimates.
