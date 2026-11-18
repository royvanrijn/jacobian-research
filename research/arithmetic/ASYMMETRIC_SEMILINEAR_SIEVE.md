# A half-dimensional sieve criterion for two modulus-12 summands

Let
\[
 C=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod{12}\},\qquad
 A=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod3\},
\]
and put \(U=\tau^2\mathbf1_C\), \(V=\tau\mathbf1_A\).
This note proves an exact reduction from a sufficient distribution estimate
to representations by \(C+C\). The distribution estimate remains open.
The results below do not establish an unconditional two-sided theorem.
The separate [critical switching bound](ASYMMETRIC_CRITICAL_SWITCHING.md)
now supplies an alternative sufficient route using only square-root-level
distribution up to logarithmic losses. Its two-large-prime upper bound is
unconditional; that distribution estimate remains open. The stronger
beyond-half criterion below remains valid.

Write
\[
 X_N=N\mathfrak S_{\rm asym}(N)
\]
for the positive main term in Theorem AsymPointwise. Its asymptotic is
available for every sufficiently large \(N\equiv2\pmod6\).
Throughout this note, \(N\equiv2\pmod{12}\), and
\[
 w_N(b)=U(N-b)V(b),\quad
 A_d(N)=\sum_{\substack{1\le b<N\\d\mid b}}w_N(b).
 \tag{1}
\]
All sieve divisors are squarefree and supported on
\(\mathcal P_7=\{p:p\equiv7\pmod{12}\}\).

## Exact local density and the target factors

If \((d,N)>1\), then \(A_d(N)=0\): any common prime would divide
\(N-b\in C\). Define the multiplicative squarefree density by
\[
 g_N(p)=
 \begin{cases}
 \dfrac{2p-1}{p^2-p+1},&p\in\mathcal P_7,\ p\nmid N,\\
 0,&\text{otherwise},
 \end{cases}
 \qquad r_d(N)=A_d(N)-X_Ng_N(d).
 \tag{2}
\]
Thus \(r_d(N)=0\) whenever \((d,N)>1\).

To verify (2), take \(p\in\mathcal P_7\), put \(x=1/p\), \(q=1-x\),
and use the local coefficients
\(U(p^0)=1\), \(U(p^j)=0\) for \(j\ge1\), and \(V(p^j)=j+1\).
For a unit target, its total local convolution mass and the mass with
\(p\mid b\) are
\[
 I_{p,0}=q+q^{-1}-1=\frac{1-x+x^2}{1-x},\qquad
 J_p=q\sum_{j\ge1}(j+1)x^j=\frac{2x-x^2}{1-x}.
 \tag{3}
\]
Their ratio is the nonzero value in (2).
For every positive target valuation \(r\ge1\),
\[
 I_{p,r}=q.
\]
Consequently, with \(g(p)=(2p-1)/(p^2-p+1)\),
\[
 \frac{I_{p,r}/I_{p,0}}{1-g(p)}
 =\frac{p-1}{p-2}>1.
 \tag{4}
\]
This identity retains every target valuation.

Put
\[
 P_7(z)=\prod_{\substack{p<z\\p\in\mathcal P_7}}p,\qquad
 \mathcal V_N(z)=\prod_{p<z}(1-g_N(p)).
\]
Since \(g(p)=2/p+O(p^{-2})\) and the prime class \(7\bmod12\) has density
\(1/4\),
\[
 \mathcal V_0(z):=\prod_{\substack{p<z\\p\in\mathcal P_7}}(1-g(p))
 \asymp(\log z)^{-1/2},\qquad
 \mathcal V_N(z)\ge\mathcal V_0(z).
 \tag{5}
\]

There is also the uniform main-mass lower bound
\[
 \boxed{X_N\mathcal V_N(\sqrt N)\gg
               \frac{N}{\sqrt{\log N}}.}
 \tag{6}
\]
Here is a direct verification of the target factors in this assertion.
The unit-target singular product
\[
 \mathfrak S_0=6\rho_U\rho_V
    \prod_{p\ge5}\frac{I_{p,0}}{M_{U,p}M_{V,p}}
\]
is a positive absolute constant. At \(p\equiv1\pmod{12}\), both local
coefficient sequences are increasing, and \(I_{p,r}\ge I_{p,0}\).
Indeed, writing
\(\alpha_i=U(p^i)-U(p^{i-1})\),
\(\beta_j=V(p^j)-V(p^{j-1})\), gives
\[
 I_{p,r}-I_{p,0}
 =\sum_{\substack{i,j\ge1\\\min(i,j)\le r}}
       \alpha_i\beta_jp^{-\max(i,j)}\ge0.
 \tag{7}
\]
At \(p\equiv5,11\pmod{12}\), the target ratio is
\((p-1)/(p-2)>1\). At a target prime \(p\in\mathcal P_7\) below \(z\),
its combined singular-series and sieve ratio is exactly (4).
There is at most one target prime \(p\in\mathcal P_7\) with
\(p\ge z=\sqrt N\), and its uncancelled factor
\[
 \frac{I_{p,r}}{I_{p,0}}
   =\frac{(p-1)^2}{p^2-p+1}
   =1+O(p^{-1})
\]
is bounded below by an absolute positive constant.
Multiplying these observations with (5) proves (6).

## A sufficient distribution criterion

Fix \(0<\sigma<1/2\), and set
\[
 z=\sqrt N,\qquad D=N^{1/2+\sigma},\qquad
 s=\frac{\log D}{\log z}=1+2\sigma.
 \tag{8}
\]
The open distribution condition is
\[
 \sum_{\substack{d<D\\d\mid P_7(z)}}|r_d(N)|
 \ll X_N(\log N)^{-A}
 \quad\text{for every sufficiently large }N\equiv2\pmod{12},
 \tag{9}
\]
where \(A\) is a sufficiently large fixed constant.

**Conditional representation theorem.** If (9) holds for some fixed
\(\sigma>0\), then
\[
 \sum_{\substack{a+b=N\\a,b\in C}}U(a)V(b)
 \ge
 \bigl(f_{1/2}(1+2\sigma)+o(1)\bigr)
      X_N\mathcal V_N(\sqrt N)>0
 \tag{10}
\]
for every sufficiently large \(N\equiv2\pmod{12}\), where
\[
 f_{1/2}(s)=
 \sqrt{\frac{e^{\gamma_E}}{\pi s}}
       \int_1^s\frac{dt}{\sqrt{t(t-1)}}>0
 \qquad(1<s\le2).
 \tag{11}
\]
In particular this would give \(\gg_\varepsilon N^{1-\varepsilon}\)
distinct ordered value pairs for every fixed \(\varepsilon>0\).

To identify the sifted set exactly, observe that every nonzero \(w_N(b)\)
has \(b\equiv1\pmod{12}\), while every prime factor of \(b\in A\) is
\(1\) or \(7\bmod12\). Hence
\[
 \sum_{p\in\mathcal P_7}v_p(b)\equiv0\pmod2.
 \tag{12}
\]
If such a \(b<N\) has no \(\mathcal P_7\)-prime below \(\sqrt N\), any
positive sum in (12) is at least two and gives a prime-factor product
greater than \(N\), including multiplicities. Thus it must be zero.
Since \(N\equiv2\pmod{12}\) is not a square, no prime endpoint is lost.
Therefore sifting to \(z=\sqrt N\) selects precisely \(b\in C\).

For completeness, the uniform sieve hypotheses can be checked in the
original half-dimensional formulation. Set \(\omega_N(p)=p g_N(p)\).
Then \(0\le\omega_N(p)<p\), and for all \(z\ge w>1\),
\[
 -L_N+\tfrac12\log(z/w)
 \le \sum_{w\le p<z}g_N(p)\log p
 \le \sum_{w\le p<z}\frac{g_N(p)}{1-g_N(p)}\log p
 \le K+\tfrac12\log(z/w).
 \tag{13}
\]
The upper constant \(K\) is absolute and independent of \(N\):
deleting target primes only decreases that sum. The fixed-progression
prime number theorem gives the required \(O(1)\) error for the undeleted
prime sums; their \(O(p^{-2})\) corrections are summable with \(\log p\).
A finite admissible lower constant is, for example, an absolute constant
plus \(\sum_{p\mid N,\ p\in\mathcal P_7}g(p)\log p\).

Iwaniec's Theorem 1 has lower error
\[
 O\!\left(e^{2K-s}(\log z)^{-\delta_I}\right),\qquad
 \delta_I=\frac12-\frac1{(2-\log2)e}>0,
 \tag{14}
\]
and its implied constant is absolute. This error contains \(K\), not
\(L_N\). Thus the lower estimate is uniform for this target-dependent
density. The weights \(w_N(b)\) are nonnegative integers, so the sieve
applies with multiplicity. Formally, for the finite product \(P_7(z)\),
each labelled copy can be replaced by a distinct positive integer in the
same residue class modulo \(P_7(z)\). All divisibility counts occurring
in the theorem are thereby preserved.

Equations (5), (9), and (11)--(14) give (10), with a remainder smaller than
\(X_N\mathcal V_N(z)\). Divisor bounds for \(U(a)V(b)\), together with
(6), give the stated number of distinct pairs.

## The signed lower-sieve condition

The lower half-dimensional weights may be chosen with
\[
 |\lambda^-_d|\le1,\qquad
 \operatorname{supp}(\lambda^-)\subset
     \{d<D:d\mid P_7(z)\},
\]
so that
\[
 \mathbf1_{(b,P_7(z))=1}\ge\sum_{d\mid b}\lambda^-_d,
 \qquad
 \sum_d\lambda^-_d g_N(d)
 \ge\bigl(f_{1/2}(s)+o(1)\bigr)\mathcal V_N(z).
 \tag{15}
\]
These are the weights \(\mu(d)\mathbf1_{d\in\mathcal D_D^-}\)
in Iwaniec's equations (3.5)--(3.6) and Lemma 17.
Consequently it suffices to establish the weaker open condition
\[
 \sum_{d<D}\lambda^-_d r_d(N)
 \ge-o\!\left(X_N\mathcal V_N(\sqrt N)\right).
 \tag{16}
\]
An absolute \(o\)-bound for this signed sum also suffices.
Condition (9) is a convenient stronger criterion.

## Exact finite-kernel compression

For every squarefree \(d\) supported on \(\mathcal P_7\) and every
positive integer \(t\),
\[
 V(dt)=2^{\omega(d)}
       \sum_{k\mid(d,t)}\frac{\mu(k)}{2^{\omega(k)}}V(t/k).
 \tag{17}
\]
At a prime \(p\mid d\), its local identity for \(e=v_p(t)\ge1\) is
\(e+2=2(e+1)-e\); at \(e=0\) both sides equal two. If \(t\) has
a forbidden prime factor for \(A\), every term vanishes. Multiplication
of these local identities proves (17), with arbitrary valuations in \(t\).

Define
\[
 C_c(N)=\sum_{\substack{a+ct=N\\a,t\ge1}}U(a)V(t).
\]
Substitution of (17) into (1) gives the exact identity
\[
 A_d(N)=2^{\omega(d)}
       \sum_{k\mid d}\frac{\mu(k)}{2^{\omega(k)}}C_{dk}(N).
 \tag{18}
\]
For a coefficient \(c\) supported on \(\mathcal P_7\), define its local
main model by zero when \((c,N)>1\), and otherwise by
\[
 H(c)=\prod_{p\mid c}(1-p^{-1}+p^{-2})^{-1},\qquad
 \mathcal M_c(N)=\frac{X_N}{c}H(c).
 \tag{19}
\]
This definition is an exact local model; no growing-coefficient
asymptotic is asserted. At \(p\mid c\), a unit target forces the first
coordinate to be a unit, so its mass is
\(M_{V,p}=(1-p^{-1})^{-1}\), independently of \(v_p(c)>0\).
Dividing this by (3) gives the factor in \(H(c)\).
Therefore
\[
 2^{\omega(d)}
       \sum_{k\mid d}\frac{\mu(k)}{2^{\omega(k)}}\mathcal M_{dk}(N)
 =X_Ng_N(d).
 \tag{20}
\]
For \((d,N)=1\), this follows also by factoring the sum:
\[
 \frac{2}{p}\left(1-\frac1{2p}\right)
       \frac1{1-p^{-1}+p^{-2}}
 =\frac{2p-1}{p^2-p+1}.
\]
For \((d,N)>1\), both sides of (20) are zero.

Let \(|\lambda_d|\le1\), supported on squarefree \(\mathcal P_7\)-divisors
\(d<D\), and let \(2\le K\le D\le N\). Then
\[
\begin{split}
 \sum_{d<D}\lambda_d r_d(N)
 ={}&
 \sum_{\substack{k\le K\\k\ {\rm squarefree}\\
                  p\mid k\Rightarrow p\in\mathcal P_7\\(k,N)=1}}\mu(k)
 \sum_{\substack{\ell<D/k\\\ell\ {\rm squarefree}\\
                  p\mid\ell\Rightarrow p\in\mathcal P_7\\(\ell,kN)=1}}
 2^{\omega(\ell)}\lambda_{k\ell}
 \left(C_{k^2\ell}(N)
       -\frac{X_NH(k\ell)}{k^2\ell}\right)\\
 &+O_\varepsilon(N^{1+\varepsilon}/K).
\end{split}
 \tag{21}
\]
The coefficients may in particular be the lower-sieve weights, extended
by zero off \(d\mid P_7(z)\).

Without the restriction \(k\le K\), equation (21) follows exactly from
(18)--(20) by writing \(d=k\ell\). To bound its omitted actual terms,
use positivity and
\[
 C_{dk}(N)\ll_\varepsilon\frac{N^{1+\varepsilon}}{dk},
\]
which follows by counting \(t\le N/(dk)\) and applying divisor bounds.
It remains valid when \(dk>N\), since the correlation is then zero.
Thus their sum is
\[
 \ll_\varepsilon N^{1+\varepsilon}
       \sum_{k>K}\frac1{k^2}
       \sum_{\ell<D/k}\frac{2^{\omega(\ell)}}{\ell}
 \ll_\varepsilon\frac{N^{1+\varepsilon}}K
                   (\log(2D))^2.
\]
The model terms obey the same bound using
\(X_N\ll_\varepsilon N^{1+\varepsilon}\) and
\(H(d)\ll_\varepsilon d^\varepsilon\).
Reducing the auxiliary divisor-bound exponent absorbs the logarithms
and proves (21). This bound is uniform in \(N,D,K\) and in the
coefficients \(\lambda_d\).

## The remaining distribution estimate

Take \(K=N^\kappa\), with a sufficiently small fixed \(\kappa>0\).
By (6), the tail in (21) is negligible relative to the positive sieve
main term. The remaining open task is therefore a lower bound for
\[
 \sum_{\substack{k\le N^\kappa\\(k,N)=1}}\mu(k)
 \sum_{\substack{\ell<D/k\\(\ell,kN)=1}}
 2^{\omega(\ell)}\lambda^-_{k\ell}
 \left(C_{k^2\ell}(N)-\frac{X_NH(k\ell)}{k^2\ell}\right),
 \tag{22}
\]
with the squarefree and prime-support conditions of (21), sufficient
to imply (16). It is a specific signed average. No individual
large-coefficient asymptotic or \(L^1\) error over all coefficients
is being assumed.

Expanding the two-character divisor series in \(C_c\) gives a progression
modulo \(cq\) with
\[
 q\le\sqrt{N/c},\qquad cq\le\sqrt{Nc}.
\]
For \(c\) near \(N^{1/2+\sigma}\), this reaches
\(N^{3/4+\sigma/2}\); the small \(k\) factors change that exponent
by \(O(\kappa)\). The asymmetric pointwise proof permits a small
auxiliary coefficient, and does not cover this range. Its BGLM input
also does not admit arbitrary weights on the combined modulus \(cq\).
The periodic three-factor bound at product scale \(N\) does not reach
a modulus near \(N^{3/4}\).

A small-coefficient part can be isolated if the uniform asymmetric
family estimate is available. More precisely, suppose
\[
 C_c(N)-\mathcal M_c(N)=O(N^{1-\eta})
 \quad(c\le N^\theta,\ (c,N)=1)
 \tag{23}
\]
uniformly for coefficients supported on \(\mathcal P_7\).
The portion of (22) with \(k^2\ell\le T\le N^\theta\) is then
\[
 O\!\left(N^{1-\eta}T\log(2T)\right),
 \tag{24}
\]
because
\[
 \sum_{k\le\sqrt T}\sum_{\ell\le T/k^2}2^{\omega(\ell)}
 \ll T\log(2T)\sum_{k\ge1}k^{-2}.
\]
For \(T\le N^{\eta/4}\), this saves a fixed power.
This controls only the small-coefficient portion; (23) supplies
no estimate for the remaining large-\(\ell\) region of (22).

The stronger [full auxiliary average](ASYMMETRIC_AUXILIARY_AVERAGE.md)
now controls the weighted absolute coefficient errors for \(c\le N^\xi\),
with \(\xi\) chosen before the analytic error saving. Its application to
the exact compression gives an \(L^1\) level \(D=N^{\xi-\delta}\)
for every fixed \(0<\delta<\xi\). The critical square-root range in
the switching criterion remains open.

## Primary sieve input

H. Iwaniec, *The half dimensional sieve*, Acta Arithmetica **29**
(1976), 69--95, [DOI](https://doi.org/10.4064/aa-29-1-69-95),
[primary scan](https://matwbn.icm.edu.pl/ksiazki/aa/aa29/aa2914.pdf).

Theorem 1 and conditions \((\Omega_1),(\Omega_2)\) are on printed p.71,
in the right half of PDF page2. The statement that all implied constants
are absolute is on printed p.70 in the same spread.
Equations (3.5)--(3.6) on printed p.75 give the signed lower weights.
Lemma 17 on printed p.88 gives their lower main term with an error
depending only on the upper constant \(K\).
The primary PDF has 14 spread-pages and SHA-256
\[
 \texttt{ca6e84aafb32ee62dcd5ce76679c65652bf677ccd4e36d200576024fac1eca67}.
\]
No beyond-half distribution estimate is imported from this sieve theorem.

