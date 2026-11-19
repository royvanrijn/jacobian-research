# Imported analytic inputs for pointwise primitive norm correlations

19 November 2026. Source-validation record for the two recent inputs used in the binary-correlation argument. **Source role: acceptable attributed imports after internal written proof validation.** The status recorded here is internal mathematical validation: the applicable proof chains below were read and reconstructed, and no concrete mathematical gap was found in those chains. This is not external human peer review, formal verification, or validation of every theorem in either manuscript. Finite computations supplement the analytic reconstruction; they do not prove the infinite estimates.

The applicable source claims and the deductions from them are kept separate. In particular, the graph input does not itself assert the binary correlation asymptotic, and the stronger seven-eighths result is unnecessary here.

## 1. Exact imported statements and primary editions

### ZF: a fixed all-conductor zero-free half-plane

The input needed by the correlation argument is the existence of one fixed \(\sigma_0<1\) such that every Dirichlet \(L\)-function is nonzero on \(\Re s>\sigma_0\), apart from a possible principal pole at one. Both conductor and height are unrestricted. An effective constant or a uniform threshold for the correlation target is not part of this assertion.

The inspected source is OpenAI, *The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane \(\Re(s)>7/8\)*, 30 September 2026. Its **Theorem 3.1**, printed p. 10, states that every finite-order Hecke \(L\)-function over \(F=\mathbb Q(\sqrt{-3})\), and every Dirichlet \(L\)-function including \(\zeta\), has no zero in \(\Re s>11/12\); the principal pole is allowed. This Part I result is sufficient for ZF. Theorem 1.1 and Part II's sharper \(7/8\) assertion are not used or validated by this record.

Primary [PDF](https://raw.githubusercontent.com/openai/math/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/paper.pdf), SHA-256 `8fe93046f8cf5ef1ba5969c89addc02d76311adc4ee907509ff9cd96f7ec99e7`.
Primary [LaTeX](https://raw.githubusercontent.com/openai/math/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex), SHA-256 `42a5ee0febca59fd1def55cfd6c6808c322ef4237d7726303ea52f711deac6a3`.
The LaTeX was inspected where text extraction loses conjugation bars.

### G: weighted-dilation root-moment transference

The graph input is OpenAI, *Weighted dilation graphs, smooth shifted primes and totient fibers*, 24 September 2026, **Definitions 3.1–3.3, Lemmas 3.4 and 3.6, and Theorem 3.5**, with the endpoint consequence Corollary 3.11 when used. Its Section 3 proof is independent of ZF and of the desired binary asymptotic.

Primary [PDF](https://raw.githubusercontent.com/openai/math/main/preprints/Weighted-Dilation-Graphs-Smooth-Shifted-Primes-and-Totient-Fibers-September-24-2026/paper.pdf), SHA-256 `351692885d112277e969901e75fe56eeee0d263d5d5847d0b980fe2c8c7744c0`.

For clarity, the exact scope is as follows. Put \(L=\log x\), \(T=\log L\). There are \(\asymp T\) disjoint groups of each of two sizes, with primes respectively in \([\exp(L^a),\exp(L^b)]\) and \([\exp(L^c),\exp(L^d)]\), where
\[
0<a<b<c<d<0.47,\qquad v_-\le V_g=\sum_{p\in P_g}p^{-1}\le v_+.
\]
The laws are \(\mu_g(p)=(V_gp)^{-1}\). The finite sets need not be intervals in this theorem. There are \(M=J+\ell\) ordered physical marks per group, distinct within each group, all dividing the physical integer. The signed comparison coefficients may depend only on the specified big labels, with at most \(m\) freed first-\(J\) slots per big group and \(\sum_\nu\sup|K_\nu|\le L^{C_1}\). Both \(m,C_1\) are fixed independently of \(J\). Source and target permutation projections are retained.

If the corresponding independent-label ideal operator satisfies \(\sup_\Theta\|T(\Theta)\|\le L^{-E_{\rm id}}\), then for every requested fixed \(E>0\), a sufficiently strong \(E_{\rm id}\) and then sufficiently large fixed \(J\) give, with \(R=\lceil L^{0.52}\rceil\), \(N=2R\),
\[
\frac1{|I\cap\mathbb Z|}\sum_{n\in I\cap\mathbb Z}
 \langle u_n,(AA^*)^Ru_n\rangle\le L^{-EN},
\qquad I=[X,2X),\quad x^\gamma\le X\le x^{1/\gamma}.
\]
Here \(u_n\) is one on **all** mark lists at integer position \(n\). Thus only the integer position returns; initial and final mark lists need not agree. The source permits \(0<|k|\le L^{C_2}\); the imported application needs only \(k=1\). A target-sized shift, a fixed-progression ideal construction, and a quotient character frame are additional application arguments, not assertions silently added to G.

## 2. Part I proof dependencies and their inspected scope

The substantial published inputs were checked in their applicable primary statements.

- Dunn–Radziwill, *Bias in cubic Gauss sums: Patterson's conjecture*, [arXiv:2109.07463v3](https://arxiv.org/abs/2109.07463v3), PDF SHA-256 `df888f8527b3647a3428d6677a3e3070aab7439f79f737f3e894218fb64cfef0`: the theta automorphy, three cusp expansions and Appendix A formulas. These formulas are unconditional. Their GRH-conditional prime theorem is not imported.
- Blomer–Goldmakher–Louvel, *L-functions with n-th order twists*, [arXiv:1112.1650v1](https://arxiv.org/abs/1112.1650v1), PDF SHA-256 `1ab29e362c592b9f4078eb85e495ecce6a4665972e6427e7cc2ab0fb214e303b`: Theorem 1.3 and Sections 2–3, including the \(K+D+(KD)^{2/3}\) squarefree higher-order sieve and its recursion. The Eisenstein sextic case is in the stated field/order range.
- Goldmakher–Louvel, *A quadratic large sieve inequality over number fields*, [arXiv:1112.1642v2](https://arxiv.org/abs/1112.1642v2), PDF SHA-256 `5cc5f03388d3996b4333f07be51f989f60779fe172f8bd24c7aa6ddad64dd304`: Definition 1, Theorem 1.1/Corollary 1.2, the fixed-family conditions and the relevant monotonicity steps.

This does not replace the source's new reflection or primary-symbol calculation by a citation. Those interfaces require, and received, the case checks below.

### 2.1 Elementary analytic controls and the detector

The strip bound in Lemma 4.8 follows from the primitive Hecke functional equation, Euler convergence, Stirling, and a Gaussian-damped rectangle argument. It gives \(|L(\sigma+it,\psi)|\ll Q^{3/5}(3+|t|)^2\) on the stated strip. Fixed-conductor qualitative growth is used only before taking the rectangle height to infinity; it does not enter the resulting uniform constant.

Lemma 4.9 starts with a **known zero-free disk**. Borel–Carathéodory bounds its holomorphic logarithm by \(O(\log\mathcal C)\) on an interior disk. Three-circles interpolation from an Euler disk gives \(O((\log\mathcal C)^\theta)\), with fixed \(\theta<1\), on the retained disk. Exponentiation yields both \(L\) and \(L^{-1}\ll_\epsilon\mathcal C^\epsilon\). This does not assume the sought global zero-free conclusion.

Section 8's buffered bins, Gamma detector and saturation were checked in full at the interface to Section 9. The complete Möbius divisor sum kills every coefficient with \(1<n\le D_*\); the terminal cutoff therefore removes precisely the unit term. The product witness has length \(\le U^{21+o(1)}\) before saturation. The two factors share the same zero and twist frequency. Comparing their product spike to the two pointwise bounds forces the plain length \(m\le1/2+O(\epsilon)\) and inverse length \(r\ge1/2-O(\epsilon)\). Only finitely many smooth parameters vary rowwise; their restoration costs a fixed height power.

### 2.2 Completed reflection and low energy

Proposition 5.1's three rational-cusp cases were reconstructed: \(v_\lambda(c)\ge2\), \(v_\lambda(c)=1\), and \(v_\lambda(c)=0\), including a zero or nonunit fixed Fourier frequency. The same three unconditional cusp expansions suffice. The derivative in the completed series removes every cusp constant; entire continuation does not pretend that cubic theta is cuspidal.

The conjugation of cusp coefficients, \(\lambda\)-mask, primal/dual Gauss orientation, imprimitive finite transforms, and all six local exponents were retained. In particular a principal local exponent is zero-extended, so its nonzero Fourier frequency is \(-1/Q\), not zero. The exceptional exponent-four branch retains its collision indicator. At the dual prime above two only a cubic character is used; no nonexistent sextic character of its residue field is introduced.

The Mellin reflection constant is \(-i/81\), with scale \((2\pi)^4/27\), in the stated self-dual normalization. A common joint kernel is separated before the positive row norm is enlarged. The powerful-row contribution is counted once, giving the source's \(M=N=1\) row energy \(\ll Z^{1+\epsilon}\). The planar additive norm has coefficient mass \(Y^{-1}\) and spacing \(Y^{-1}\). Consequently Proposition 6.3 gives
\[
I_\eta(Z^{1/2},Z^{1/2},Z)\ll_{\eta,\epsilon}Z^{1/4+\epsilon}.
\]
The exponent, unlike its fixed-target constant, is common to the targets.

### 2.3 Scalar Euler calibration and a uniform nonzero correction

Section 7's finite scalar was reconstructed for all prime-power cases, including noncoprime indices. Poisson normalization leaves \(\sqrt X/(\sqrt{q_B}q_A)\), with no residual \(q_s\). Primitivity at the fixed bad conductor removes the zero Fourier frequency. The bad Gauss factors cancel with the stated character normalizations. Extracting a primary sixth power is unique and introduces no factor six.

The global cross-prime calculation leaves precisely the calibrated local phase
\[
(-1)^e\gamma_1(p)^{-e}\eta(p)^e
 (a_p\overline{\gamma_3(p)})^l\rho^{k-t}
 \omega_p^{tk-el-l(l-1)/2}.
\]
The bar on \(\gamma_3\), the diagonal reciprocity, the bar on \(\chi_p(4)\) in \(b_p\), and the unbarred target \(\eta\) were checked against the primary LaTeX. All ten valuation-family branches give the source's exact local rational expression.

Here is a quantitative check of the analytic correction rather than an appeal to a qualitative Euler-product constant. With \(V=Q^{-6z}\), \(W=\chi_p(u)Q^{-w}\), \(D=\eta(p)\overline{\chi_p(u)}Q^{-x}\), and \(E_p=P_p^*+D\), exact algebra gives
\[
H_p-1=\frac{D(V+W-VW)-VW+(1-V)(1-W)E_p}{1-D}.
\]
Thus no \(1-W\) denominator remains when \(w\) is moved through zero. In the principal-row region \(\Re x\ge7/8\), \(\Re z\ge33/200\), \(\Re w\ge19/20\), direct termwise estimation gives
\[
|H_p-1|\le144Q^{-363/200}.
\]
The elementary Eisenstein ideal count \(A(T)\le2T\) implies that deleting all prime ideals of norm \(\le65536\) gives
\[
\sum_{p\notin S}|H_p-1|<\frac{3267}{41728},\qquad
 |H_\eta(s)-1|<\frac{3267}{38461}<\frac12
 \quad(\Re s>7/8).
\]
Additional conductor deletions only decrease this majorant. The correction is holomorphic and nowhere zero uniformly in targets and heights. For nonprincipal physical rows the source's wider-region bounds retain the finite \(q_u^\epsilon\) factor.

For \(u=1\) the exact Euler identity is
\[
F_{\eta,1}(x,w,z)=
 \frac{\zeta_F^S(6z)\zeta_F^S(w)}{L_F^S(x,\eta)}H_{\eta,1}(x,w,z).
\]
There is exactly one principal numerator row. The double residue has a positive coefficient
\(c_S=\widehat W_1(1)M(1/6)(\operatorname{Res}_{s=1}\zeta_F^S(s))^2/6\).
The factor \(1/6\) is the Jacobian at \(\zeta_F(6z)\), counted once. Positivity of \(M(1/6)\) follows by Fourier duality against the positive Riesz kernel. Hence the principal signal contains the reciprocal of the actual target, with a coefficient that cannot cancel a zero.

### 2.4 Sextic recursion and physical row envelope

A raw sextic symbol can be nontrivial on units, so it is not simply identified with a published Hecke family. For squarefree denominators in the same fixed primary class, their paired character is unit-trivial, has the exact primitive conductor, and has a separated Gauss factor. The source's paired Poisson recursion was checked against the full published recursion, including conjugates, transpose/inverse orientation, monotonicity and the closure \(\{1,2,4\}\) under doubling modulo six. No principal kernel is accidentally introduced. Iterating the exponent improvement finitely many times for a prescribed loss suffices; no uniform limiting constant is required.

For Proposition 9.2 the powerful row part is frozen only inside the current positive norm sum. Its character remains in one common coefficient sequence; its primes are not added to the fixed bad set. The independent powerful-row count and the sieve moment give exactly
\[
R(\delta)=\min\{1,\max(1-\delta/2,4/3-\delta)\}.
\]
The smooth row parameters have a fixed-order Sobolev restoration before the final small height power is chosen. Constants may depend on the fixed target; the power losses do not.

## 3. Explicit closure of Part I

This closes the infinite implication from the checked estimates; it is separate from finite normalization tests.

Let \(\beta_*\) be the supremum of real parts of zeros over all primitive finite-order Hecke targets, with \(1/2\) adjoined, and suppose \(\Delta=\beta_*-11/12>0\). Put \(C(s)=s-2/3\). The normalized physical probe \(J_\eta\), its principal signal
\[
f_\eta(Z)=\frac1{2\pi i}\int_{\Re s=2}
 Z^{C(s)}e^{(s-5/6)^2}\frac{H_\eta(s)}{L_F^S(s,\eta)}\,ds,
\]
and the correction above are independent of the later analysis height. The low estimate gives \(J_\eta(Z)\ll_\eta Z^{1/4+\Delta/2}\).

The central, principal and outer-row estimates in Sections 10–11 have respective target-independent error-free margins
\[
1021/25000,\quad1/40,\quad1/1200,\quad43/300,\quad>1
\]
relative to \(C(\beta_*)\). In the central range the exponent is
\[
-\frac34+\frac{\delta+R(\delta)}2
 +(d-\tfrac12)(R(\delta)+\delta/2-17/50),
\]
with the actual ceiling \(\delta\le2\beta_*-1\); this gives the first displayed margin, including the whole retained \(d\)-range. The real losses are fixed before the target: \(46\cdot10^{-6}<1/4800\). Detector ranges up to length exponent 22 are fixed before saturation.

For each fixed target all retained terms therefore satisfy, for every external integer \(N\),
\[
|J_\eta-f_\eta|\ll_{\eta,N}
 Z^{C(\beta_*)-1/2400}(1+T_1)^{A_\eta}
 +Z^{B_\eta}T_1^{-N},
\]
with \(A_\eta,B_\eta\) finite and fixed before \(N\). Choose \(T_1=Z^{\tau_\eta}\) only now, with \(\tau_\eta>0\) small enough to satisfy both the detector ceiling and \(\tau_\eta(A_\eta+1)\le1/9600\). Next choose external \(N\) large enough. Increasing this external order changes only test seminorms, not \(A_\eta,B_\eta\). The result is
\[
|J_\eta-f_\eta|\ll_\eta Z^{C(\beta_*)-1/4800}.
\]
The saving \(1/4800\) is common to all targets, although \(\tau_\eta\), constants and lower thresholds may depend on the target.

For completeness, this really yields a contradiction without moving a contour through a zero. The triangle inequality bounds \(f_\eta(Z)\ll Z^{C(\beta_*)-\epsilon_*}\), where \(\epsilon_*=\min(\Delta/2,1/4800)>0\). At zero, moving its defining line right gives arbitrary power decay. Thus \(\int_0^\infty f_\eta(Z)Z^{-C(s)}dZ/Z\) is holomorphic on \(\Re s>\beta_*-\epsilon_*\). Fourier inversion on the original line identifies it there by continuation with \(e^{(s-5/6)^2}H_\eta(s)/L_F^S(s,\eta)\). Since \(H_\eta\ne0\), this continues every target reciprocal holomorphically to that half-plane. The common \(\epsilon_*\) contradicts the definition of the supremum; attainment is unnecessary. The principal pole becomes a harmless zero of the reciprocal.

The Dirichlet transfer follows from \(L_F^S(s,\chi\circ N)=L^S(s,\chi)L^S(s,\chi\chi_{-3})\). There is no pole-zero cancellation in the open critical strip. The only possible cancellation at one reduces to \(L(1,\chi_{-3})=\pi/(3\sqrt3)>0\). Deleted Euler factors are nonzero for \(\Re s>0\). This gives the exact Part I statement, with no claim on its boundary.

## 4. Explicit closure of G's root-moment proof

The full Section 3 proof was reconstructed, including the signed memory argument. No zero-free input or desired binary asymptotic enters it.

First the integer-root average is replaced by independent prime residues. For a path of length \(N\), active moduli and a Bonferroni truncation of the remaining damping have total logarithm \(O(NT L^d)=o(L)\), since \(0.52+d<1\). Active damping is factored exactly before comparison. The independent avoidance probability is at least \(\exp(-O(NT))\); the interval CRT error, divided by this probability, is still negligible. Absolute path row/column bounds permit summation of this relative error. Fourier integration imposes only integer return.

For a never-active prime the exact one-prime expansion is
\[
b_p+\frac1p\sum_{j<i}\eta_j\eta_i1_{h_j=h_i\ (p)}
 \prod_{j<t<i}(1-\eta_t1_{h_t=h_j\ (p)}),
\qquad b_p=1-p^{-1}\sum_i\eta_i.
\]
It supplies exactly one absent-or-orphan lifespan. An active prime with \(r\) maximal active runs costs \(V_g^{-r}/(pb_p)\), not \(p^{-r}\). The memory Hilbert space uses reciprocal-prime mass times **counting** residue measure and factorial measure on pending particles. With \(\rho^2>\sqrt q\), its birth, hit and termination factors reproduce every original damping factor, including the two boundary factors \(\sqrt q\). Distinct births reproduce each original path uniquely.

Memory truncation and ghost creation have uniform weighted row/column bounds. Transfers between pending memory and active labels are rare: the product of omitted small labels is less than every big prime, so each pending residue permits at most one of \(\binom{J+\ell}{\ell}^{|S|}\) omission choices. The coefficient's independence from small labels is essential. Clean edges Fourier-diagonalize to the assumed ideal operator; changes of prime measure and auxiliary collisions cost less than every fixed logarithmic power.

Distinct births are restored by exact equality-partition Möbius inversion. Rank \(r\) uses at most \(2r\) addresses and total coefficient \(\le L^{C'r}\). At large rank, every dependent fresh birth costs a maximum atom \(\ll\exp(-L^c)\); backward absolute row bounds multiply these savings. At small rank, equality is encoded by phases at the birth addresses only. At most \(2r\) edges lose signed contraction, while all other edges retain the complete signed operator. With \(r_0/N=O(T/L^c)+O(1/N)=o(1)\), the rank restoration preserves contraction on \(N-o(N)\) edges.

The order of choices closes: ghost exponent \(C_v\) and small-state exponent \(C_s\) are independent of the eventual \(J,m,C_1\); choose \(G>E+C_v+3\), then \(E_{\rm id}>G+C_s+2\). After \(m,C_1\) are fixed, choose \(J\) large enough for the rare-transfer bound. Rank cutoffs and residue precision are chosen last. This proves the stated root moment with its quantifiers. Corollary 3.11 follows from positive-semidefinite trace domination and spectral Hölder for a first endpoint constant on the mark fiber; no equality of fiber norms or final mark lists is presumed.

The binary application must still verify its own big-label coefficient restriction, \(J\)-independence, group masses, endpoint norms, fixed-progression character frame and large-target rescaling. G does not waive these hypotheses.

## 5. Scope of validation

The applicable Part I and G proof chains have written internal
reconstructions. Their exact statements are imported on the basis of those
analytic arguments and the hypothesis checks above. The source editions
and preprint dependency remain explicit. Part II and unrelated
applications are outside this inspection. No external human peer review,
external consensus or local formal-verification run is claimed.

The binary application verifies its separate interfaces in
[the supporting arguments](PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md).
The complete assembly and resulting leading terms are in
[the main theorem](PRIMITIVE_NORM_BINARY_CORRELATION.md).
