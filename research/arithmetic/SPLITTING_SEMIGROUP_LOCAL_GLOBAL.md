# A Goldbach theory for abelian prime-splitting semigroups

The principal result is a pointwise correlation theorem for two arbitrary
fixed abelian subgroup prime prescriptions, with a positive logarithmic
remainder and an exact congruence-admissibility criterion.
Its defining specialization is the full \(C_{12}+C_{12}\) asymptotic.
Together with classical higher correlations, this gives strong and weak
Goldbach laws, explicit singular series and exact asymptotic basis orders.
The higher-fold counts are classical consequences; the binary
pointwise correlation is the part requiring the new analytic adaptation.

The analytic proof uses the precisely scoped weighted-dilation and fixed
Dirichlet zero-free theorems recorded in the
[input audit](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md).
The density model, progression frame, growing-target comparison and
target-grid estimates are proved in the
[binary-correlation support](PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md).
Module J there supplies the additional quantitative cutoff argument.
The source applications and new arguments have internal written review.
No external human review, formal verification, numerical onset or
historical-priority claim is asserted.

## 1. The defining case

Let
\[
 C_{12}=\{n\ge1:p\mid n\Longrightarrow p\equiv1\pmod{12}\},\qquad
 {\cal E}_{12}=\{p:p\nmid6,\ p\not\equiv1\pmod{12}\},
\]
and count ordered pairs \(R_{12}(N)=\#\{(a,b)\in C_{12}^2:a+b=N\}\).
Define the absolutely convergent positive product and target factor
\[
 A_{12}=\frac{4\pi}{\Gamma(1/4)^2}
       \sqrt{\frac{\log(2+\sqrt3)}6}
       \prod_{p\in{\cal E}_{12}}
                    \left(1-\frac2{p(p-1)}\right),\qquad
 J_{12}(N)=\prod_{\substack{p\mid N\\p\in{\cal E}_{12}}}
                                      \frac{p-1}{p-2}.          \tag{S.1}
\]

**Theorem C12QuantitativeCorrelation.** Uniformly over sufficiently large
\(N\equiv2\pmod{12}\),
\[
 R_{12}(N)=A_{12}J_{12}(N)\frac N{(\log N)^{3/2}}
                    \left(1+O((\log N)^{-1/32768})\right).       \tag{S.2}
\]
Every other residue class has zero count, identically.
Moreover \(N\equiv2\pmod{12}\) is exactly the condition for a
representation in every finite quotient. It follows that every
sufficiently large congruence-admissible target has an actual positive
representation.

This is an unweighted count of numerical pairs. It is neither a
coordinate representation count nor a weighted divisor convolution.
The threshold and the implicit constant are not given numerically.
The theorem gives no exact classification of the small targets.

### Exact constant

The four literal characters modulo twelve give the Euler identity
\[
 F_{12}(s)^4=
 \zeta(s)L(s,\chi_{-4})L(s,\chi_{-3})L(s,\chi_{12})
 (1-2^{-2s})(1-3^{-2s})
 \prod_{p\in{\cal E}_{12}}(1-p^{-2s})^2,                         \tag{S.3}
\]
where \(F_{12}(s)=\sum_{n\in C_{12}}n^{-s}\) and the three displayed
nonprincipal characters are primitive. Thus the marginal constant is
\[
 c_{12}=\frac1{\Gamma(1/4)}
 \left[\frac{\pi^2\log(2+\sqrt3)}{54}
                  \prod_{p\in{\cal E}_{12}}(1-p^{-2})^2
 \right]^{1/4}.                                                \tag{S.4}
\]
Here \(L(1,\chi_{-4})=\pi/4\),
\(L(1,\chi_{-3})=\pi/(3\sqrt3)\), and
\[
 L(1,\chi_{12})
 =\int_0^1\frac{1-x^4}{1+x^6}\,dx
 =\left[\frac1{2\sqrt3}
   \log\frac{x^2+\sqrt3x+1}{x^2-\sqrt3x+1}\right]_0^1
 =\frac{\log(2+\sqrt3)}{\sqrt3}.
\]
The bad Euler factor is \((1-1/4)(1-1/9)=2/3\).
The finite additive density is twelve. Therefore the general theorem
below gives \(12c_{12}^2 C_{{\cal E}_{12}}\), where
\(C_{\cal E}=\prod_{\cal E}(1-(p-1)^{-2})\).
Multiplying each factor \(1-p^{-2}\) into this product gives exactly
(S.1). This also fixes the Gamma normalization; there is no additional
beta-integral factor.

## 2. General quantitative correlation

Fix an even integer \(L_0\ge2\), put \(G=(\mathbb Z/L_0)^\times\),
and let \(H_1,H_2\le G\) be arbitrary subgroups. Define
\[
 {\cal P}_i=\{p:p\nmid L_0,\ p\bmod L_0\in H_i\},\quad
 C_i=\langle{\cal P}_i\rangle,\quad
 \rho_i=\frac{|H_i|}{|G|},\quad \kappa=2-\rho_1-\rho_2.
\]
The two subgroups need not be nested.
Fix \(U\subset H_1,V\subset H_2\), and coprime positive
\(L_0\)-supported coefficients \(A,B\). Put
\[
 R_{U,V}^{A,B}(N)
 =\#\{Ax+By=N:x\in C_1,\ y\in C_2,\
                       x\bmod L_0\in U,\ y\bmod L_0\in V\}.
\]
Let
\[
 \begin{split}
 c_i&=\frac1{\Gamma(\rho_i)}
       \lim_{s\downarrow1}(s-1)^{\rho_i}
                             \prod_{p\in{\cal P}_i}(1-p^{-s})^{-1},\\
 \Gamma_N&=\frac{L_0}{|H_1||H_2|}
        \#\{(u,v)\in U\times V:Au+Bv=N\bmod L_0\},\\
 {\cal E}&=\{p:p\nmid L_0,\ p\bmod L_0\notin H_1,\
                                  p\bmod L_0\notin H_2\},\\
 C_{\cal E}&=\prod_{p\in{\cal E}}
                         \left(1-\frac1{(p-1)^2}\right),\qquad
 J_{\cal E}(N)=\prod_{\substack{p\mid N\\p\in{\cal E}}}
                                      \frac{p-1}{p-2}.
 \end{split}                                                   \tag{S.5}
\]
All main-term factors are explicit. For example, if \(r_i(p)>1\) is the order
of \(pH_i\) in \(G/H_i\), then
\[
 c_i=\frac{(\varphi(L_0)/L_0)^{\rho_i}}{\Gamma(\rho_i)}
 \prod_{\substack{\eta\in H_i^\perp\\\eta\ne1}}L(1,\eta)^{\rho_i}
 \prod_{\substack{p\nmid L_0\\p\notin{\cal P}_i}}
                         (1-p^{-r_i(p)})^{1/r_i(p)}.            \tag{S.6}
\]
The characters in (S.6) are literal modulo \(L_0\). Positive real
continuation from \(s>1\) specifies the product. The last product is
absolutely convergent, and \(c_i>0\).

**Theorem QuantitativeSplittingCorrelation.** With all displayed data
fixed, put \(\delta=1/(8192\varphi(L_0))\). Then
\[
 R_{U,V}^{A,B}(N)
 =c_1c_2\Gamma_N C_{\cal E}
       \frac{NJ_{\cal E}(N)}{AB(\log N)^\kappa}
 +O\!\left(\frac{NJ_{\cal E}(N)}
                  {AB(\log N)^{\kappa+\delta}}\right).           \tag{S.7}
\]
The error is uniform in the target and all its prime valuations.
If \(\Gamma_N=0\), the count is exactly zero.
For \(\Gamma_N>0\), the finitely many possible positive densities have
a positive minimum, so (S.7) is a relative asymptotic.

For every interval \(I\subset[0,1]\), restricting \(Ax/N\in I\)
multiplies the main term by \(|I|\), with the same absolute error,
uniformly in both endpoints of \(I\). Hence a sequence of intervals satisfying
\(|I|(\log N)^\delta\longrightarrow\infty\) has the corresponding
relative asymptotic. This does not cover power-short intervals.

**Proof.** Modules D/I/G/P/E/F establish the binary correlation with
the correct marginal constants, finite root factor and target-grid
normalization. Module J upgrades every error to the explicit rate.
Here are the decisive points, including the extra argument needed for
that upgrade.

The rough pole order is \(\rho_i\); every character in \(H_i^\perp\)
has the same principal branch. The projected rough model retains all
these branches, with fixed Hankel depth chosen after the required
precision. Identity-class prime banks preserve the original indicator,
fixed residue filter and model exactly.

Choose \(\lfloor\log\log N/(32\varphi(L_0))\rfloor\) mass-one banks of
each size type. The exact two-root record envelope gives success
probability approaching \(e^{-2}>1/8\); Chernoff discards bad bank
sets with relative error \(O((\log N)^{-2\delta})\).
The total exact reconstruction coefficient cost is at most \(\log N\).

For the moving bulk width \(\epsilon=(\log N)^{-\delta}\), the Mellin
norm costs at most \(\epsilon^{-1}\). On each contributing reflected
edge the target/source ratio is \((1-t)/t\), so a TWO-SIDED target cutoff
confines the target to dyads between \(c\epsilon X\) and
\(C\epsilon^{-1}X\). This cutoff is inserted before Mellin inversion.
The target norm and frequency-window changes are explicitly bounded
by fixed powers of \(\epsilon^{-1}\); Module J gives their complete
budget and noncircular order of choices.

The graph comparison retains two free big-prime products. On target
grids its cofactor length/modulus ratio exceeds every fixed logarithmic
power, even in the smallest target dyad and in every gcd profile.
The fixed zero-free theorem supplies prime saving \(P^{-1/24}\).
Exceptional fixed characters are canceled by the projected model;
other characters are bounded after extracting one unique small-bank
prime. The profile sum costs \(\sum_{d\mid q}1/d\), not \(\tau(q)\).
The quotient progression frame and residue rescaling verify the graph
hypotheses for the growing target \(N\).

For model times model, CRT and the fundamental sieve lemma give the
finite root probability. A prime forbidden on only one side has
relative factor one. A prime forbidden on both sides contributes
\[
                \frac{1-(2-\mathbf1_{p\mid N})/p}{(1-1/p)^2}.
\]
Their product is \(C_{\cal E}J_{\cal E}(N)\).
The real line has parameter length \(N/(AB)\) and uniform
archimedean density. Its integral is one.
Leading-density replacement costs \(O((\log N)^{-1/2})\).
Partial summation costs bounded total variation, and direct endpoint
sieves cost \(O(\epsilon)\) of the main scale plus a power-smaller
CRT error. The inequalities in Module J then give exactly (S.7).
The same construction with two additional transition points proves
the uniform interval assertion. \(\square\)

For noncoprime fixed \(A,B\), divide by their fixed gcd when it divides
\(N\); otherwise the count is zero. This reduces to the stated theorem.

### Balanced fixed-modulus families

**Corollary FixedModulusBalancedExistence.** For every fixed \(q\ge3\), let
\(L=\operatorname{lcm}(2,q)\), \(H=\{1\}\subset(\mathbb Z/L)^\times\), and
\(C_q=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod q\}\).
Let \(c_q,{\cal E}_q,C_{{\cal E}_q},J_{{\cal E}_q}\) be (S.5) for this
subgroup. Uniformly over sufficiently large \(N\equiv2\pmod L\),
\[
 \#\{a+b=N:a,b\in C_q,\ N/3\le a,b\le2N/3\}
 =\frac{Lc_q^2C_{{\cal E}_q}}3 J_{{\cal E}_q}(N)
       \frac N{(\log N)^{2-2/\varphi(q)}}
       \left(1+O((\log N)^{-1/(8192\varphi(q))})\right).
\]
In particular the count is positive. The threshold may depend on \(q\).

**Proof.** Every allowed prime is \(1\bmod L\): for odd \(q\), it is
odd and \(1\bmod q\); the even case is immediate. Conversely a prime
\(1\bmod L\) is \(1\bmod q\), so these are exactly the subgroup cores.
Also \(\varphi(L)=\varphi(q)\), \(\rho_1=\rho_2=1/\varphi(q)\), and
\(\Gamma_N=L\) on the stated targets. Apply (S.7) and its uniform interval
version with \(A=B=1\), \(U=V=\{1\}\), and \(I=[1/3,2/3]\).
This interval has length \(1/3\), and its restriction is equivalent to
the displayed two-sided balance. All leading factors are positive.
\(\square\)

This proves the macroscopic family FBS-G. At \(q=12\) the leading
constant is \(A_{12}/3\); it proves balanced S12, and dropping balance
proves S12U. The power-short FBS-U requirement \(|a-b|\le N^\theta\)
for each \(0<\theta<1\) does not follow from this corollary. Nor does it
prove the separate coordinate-weight hypotheses H1/H2 or the natural
raw norm-count lower bound and marked-prime marginals.

## 3. The exact meaning of local admissibility

Separate \(p\)-adic closures need not retain the coupling between
different prime components of a residue subgroup.
The correct general condition is solvability in every finite quotient.

For one subgroup \(H\le G\), put
\(C_H(U)=\{n\in C_H:n\bmod L_0\in U\}\).

**Lemma SplittingFiniteImages.** A residue \(r\bmod M\) lies in the image
of \(C_H(U)\) if and only if:

1. \(r\bmod\gcd(M,L_0)\) is the reduction of some \(h\in U\);
2. every \(p\mid M\) with \(p\nmid L_0\) forbidden by \(H\) satisfies
   \(p\nmid r\).

**Proof.** Necessity is immediate. Choose \(h\in U\) as in condition 1. For each allowed \(p^k\Vert M\)
choose \(e_p=\min(v_p(r),k)\), taking \(e_p=k\) if \(r=0\bmod p^k\).
Put \(D=\prod_{p\in{\cal P}_H,\ p\mid M}p^{e_p}\in C_H\).
Choose a prime \(q\) in the class \(hD^{-1}\bmod L_0\).
Where \(e_p<k\), impose the unit congruence
\[
 q\equiv(r/p^{e_p})(D/p^{e_p})^{-1}\pmod{p^{k-e_p}}.
\]
At primes not dividing \(D\), use \(q\equiv rD^{-1}\pmod{p^k}\);
conditions 1--2 make it a unit. Where \(e_p=k\), impose \(q\equiv1\bmod p\).
At primes dividing \(L_0\), compatibility is exactly condition 1.
CRT therefore produces a reduced residue class.
Dirichlet's theorem supplies infinitely many such primes \(q\).
Choose one outside the finite exceptional set.
Then \(Dq\in C_H(U)\) has the required residue. \(\square\)

**Theorem SplittingCongruenceLocalGlobal.** For fixed data in (S.7),
\(\Gamma_N>0\) is equivalent to the following statement:
for every \(M\ge1\), there exist actual \(x\in C_1(U)\), \(y\in C_2(V)\)
with \(Ax+By\equiv N\pmod M\).
Consequently, for all sufficiently large \(N\), congruence admissibility
is equivalent to actual positive additive representability.

**Proof.** Taking \(M=L_0\) proves necessity. Conversely fix a root
\((u,v)\) counted by \(\Gamma_N\), and put \(K=\operatorname{lcm}(M,L_0)\).
At \(p\mid L_0\), retain \(x=u,y=v\) modulo the \(L_0\)-part.
At least one of \(A,B\) is a unit; adjusting its variable lifts the
equation to \(p^{v_p(K)}\) while preserving those unit states.
At \(p\mid K\), \(p\nmid L_0\), both coefficients are units and \(p\ge3\).
The equation modulo \(p\) has \(p\) pairs. Requiring a unit on each
forbidden side excludes at most two, so a pair survives and lifts.
CRT gives simultaneous residue pairs meeting the finite-image lemma.
Realize each coordinate by an actual semigroup value using that lemma.
This proves every finite congruence. Finally (S.7) is positive on all
sufficiently large targets with \(\Gamma_N>0\). \(\square\)

For \(C_{12}\), take \(L_0=12,H_1=H_2=\{1\},A=B=1\).
Then \(\rho_1=\rho_2=1/4,\Gamma_N=12\) exactly on \(N=2\bmod12\).
Equations (S.3)--(S.7) prove the defining theorem (S.2).

**Proposition CoupledCongruenceObstruction.** Separate admissibility at
every \(p\)-adic closure is insufficient even for abelian subgroup
splitting semigroups.

Take \(L_0=12,H=\{1,11\}=\ker\chi_{12}\).
Every value of \(C_H\) is \(1\) or \(11\bmod12\), so
\[
                        H+H=\{0,2,10\}\pmod{12}.
\]
Thus every target \(N=4\bmod12\) has zero count.
Yet \(H\) projects onto the full unit groups modulo four and three.
The finite-image lemma gives closures \(\mathbb Z_2^\times\) and
\(\mathbb Z_3^\times\). At a good allowed prime the closure is
\(\mathbb Z_p\); at a good forbidden prime it is \(\mathbb Z_p^\times\).
Two units sum to every element of \(\mathbb Z_p\) for odd \(p\), and
two 2-adic units sum to every even element.
Every such \(N\) therefore passes every separate \(p\)-adic test.
The missing condition is the joint residue state modulo twelve.
These are exactly the good primes splitting in \(\mathbb Q(\sqrt3)\).

## 4. Splitting fields and finite conjunctions

Let \(F_j/\mathbb Q\) be finitely many finite abelian extensions and
\(K_j\le\operatorname{Gal}(F_j/\mathbb Q)\).
In their compositum \(F\), the simultaneous prescription
\(\operatorname{Frob}_p(F_j)\in K_j\) is exactly the subgroup condition
\(\operatorname{Frob}_p(F)\in\bigcap_j\pi_j^{-1}(K_j)\).
Kronecker--Weber places \(F\) in a cyclotomic field.
Taking its Artin inverse image gives precisely \(H_i\le G\) above.
The exact source is Milne, *Class Field Theory*, Chapter I, Theorem 4.16,
printed p. 49 in the retained edition.

Thus (S.7) covers arbitrary finite conjunctions of abelian subgroup
splitting conditions, independently on both sides and without nesting.
In an abelian extension, such a subgroup prescription is complete
splitting in its fixed field.
The theorem discards the fixed finite bad set, including two.
Adding a freely generated bad prime requires separate analysis below.
The arithmetic reduction is classical; the target-uniform binary
indicator correlation is the substantive analytic application.

## 5. Controlled bad-prime profiles

Let \(S=\{p:p\mid L_0\}\).
For finite sets \({\cal R}_1,{\cal R}_2\) of \(S\)-supported integers
and states \(U_r\subset H_1,V_s\subset H_2\), define disjoint exact-prefix
sets
\[
 P_1=\coprod_{r\in{\cal R}_1}rC_1(U_r),\qquad
 P_2=\coprod_{s\in{\cal R}_2}sC_2(V_s).
\]
An empty profile set gives identically zero counts; assume both sets
are nonempty below.
For \(g=(r,s)\), \(g\mid N\), use the ORIGINAL states in
\(\Gamma_{r/g,s/g;U_r,V_s}(N/g)\), and set
\[
 {\cal B}(N)=
 \sum_{\substack{r\in{\cal R}_1,\ s\in{\cal R}_2\\(r,s)\mid N}}
       \frac{(r,s)}{rs}
       \Gamma_{r/(r,s),s/(r,s);U_r,V_s}\!\left(\frac N{(r,s)}\right).
                                                                  \tag{S.8}
\]

**Theorem FiniteBadProfileCorrelation.** With \(\kappa,\delta\) as above,
\[
 \#\{a+b=N:a\in P_1,b\in P_2\}
 =c_1c_2C_{\cal E}{\cal B}(N)\frac{NJ_{\cal E}(N)}{(\log N)^\kappa}
 +O\!\left(\frac{NJ_{\cal E}(N)}
                       {(\log N)^{\kappa+\delta}}\right).        \tag{S.9}
\]
The rational-valued function \({\cal B}(N)\) is periodic, with period
dividing \(L_0\operatorname{lcm}\{\gcd(r,s):r\in{\cal R}_1,s\in{\cal R}_2\}\).
Its positivity is equivalent to solvability in every finite quotient.
Its zero classes have zero count identically; positive classes have
the eventual local-global property.

**Proof.** Each stratum is empty unless \(g\mid N\). Divide by \(g\)
and apply (S.7) with coefficients \(r/g,s/g\) and target \(N/g\).
Since \(g\) is bad-supported, \(J_{\cal E}(N/g)=J_{\cal E}(N)\);
the real factor is \(Ng/(rs)\).
Replacing \(\log(N/g)\) by \(\log N\) costs relative \(O(1/\log N)\),
smaller than the displayed remainder. There are finitely many strata.

Each term in (S.8) has period dividing \(gL_0\), proving periodicity.
Positivity supplies all congruence solutions by the core theorem.
Conversely take modulus
\(L_0\operatorname{lcm}\{\gcd(r,s):r\in{\cal R}_1,s\in{\cal R}_2\}\).
Any finite-profile solution modulo it has \(g\mid N\) and, after division
by \(g\), a root counted by the appropriate \(\Gamma\).
The finite positive values of \({\cal B}\) have a positive minimum.
This proves all assertions. \(\square\)

### Why unrestricted bad generators require a residual scale

Let \(C^*=\langle2,\ p:p\equiv1\pmod4\rangle\).
For every \(t\ge1\),
\[
                         R_{C^*}(2^t)=1.                       \tag{S.10}
\]
Indeed the two summands must have equal 2-adic valuation; otherwise the
smaller valuation would equal \(t\), impossible for positive summands.
If that common valuation were below \(t-1\), their odd cores would
sum to a multiple of four although each is \(1\bmod4\).
Thus both summands equal \(2^{t-1}\).

This is not an inadmissibility example: the target is represented.
It shows that an unrestricted bad generator can force a bounded
residual core even while \(N\) grows.
For completeness, its 2-adic value closure is
\(V=\{0\}\cup\bigcup_{e\ge0}2^e(1+4\mathbb Z_2)\).
With Haar measure one, the local convolution has
\[
 d_2(N)=\tfrac14\ (N\ {\rm odd}),\qquad
 d_2(N)=\tfrac38\ (v_2(N)=1),\qquad
 d_2(N)=\tfrac12d_2(N/2)\ (v_2(N)\ge2).
\]
These follow by separating odd-odd and even-even pairs in
\(V=(1+4\mathbb Z_2)\cup2V\).
In particular \(d_2(2^t)=3/2^{t+2}\).
A naive fixed-\(\log N\) relative main term proportional to
\(Nd_2(N)/\log N\) would tend to zero, contradicting (S.10).
For infinite profiles, a uniform positive admissible-density bound
together with a proved uniform tail bound is one sufficient extension
route. A different theory can retain the residual-core parameter.
Neither extension follows from summability of formal constants alone.

## 6. Proven scope and source comparison

The binary proof covers fixed abelian SUBGROUP prime conditions, their finite
conjunctions and finite bad-prime profiles.
For binary correlations, arbitrary unions of Frobenius classes are not a verbatim specialization:
their fixed-character twists have additional fractional singular
exponents, and may lack a common identity-class bank.
General nonabelian binary prescriptions require further twisted prime
distribution input and a new model/frame verification.

Loughran--Matthiesen, *Frobenian multiplicative functions and rational
points in fibrations*, JEMS 26 (2024), Theorems 1.9 and 3.10, prove
correlations for systems whose nonconstant linear parts are pairwise
linearly independent. The pair \(n,N-n\), at a fixed target, does not
meet that hypothesis. Averaging the target is a different conclusion.
Their theorem is therefore a close general comparison, not the present
binary pointwise input.

Walsh's inspected divisor-correlation theorem is for fixed shifts and
integer divisor orders; norm-indicator correlations and larger shifts
are listed as extensions. The present proof supplies the projected
fractional model, progression quotient frame, target-sized transfer,
relative target-grid energy, and now the quantitative cutoff transport.
The older asymmetric subgroup theorem by itself does not establish
(S.2). None of these comparisons certifies historical priority.

All groups, fields, moduli and coefficient/profile sets are fixed.
The binary remainder is logarithmic, not a power saving in \(N\).
The current proof retains its exact imported 2026 analytic theorems and
does not provide a numerical effective cutoff or an exact small-target
exception list.

## 7. Strong and weak Goldbach existence laws

For fixed \(k\ge2\), fix subgroups \(H_i\le(\mathbb Z/L_0)^\times\),
their cores \(C_i\), and residue filters \(U_i\subset H_i\). Write
\[
 U_1+\cdots+U_k
 =\{u_1+\cdots+u_k\bmod L_0:u_i\in U_i\}.
\]
All sums of numerical values below use positive summands and permit
repetition.

**Theorem SplittingGoldbachExistence.** For every such fixed list and
all sufficiently large \(N\),
\[
 N\in C_1(U_1)+\cdots+C_k(U_k)
 \quad\Longleftrightarrow\quad
 N\bmod L_0\in U_1+\cdots+U_k.                         \tag{G.1}
\]
The residue criterion is also equivalent, for every integer \(N\), to
solvability with actual semigroup values modulo every finite modulus.
The onset may depend on the whole fixed list.

**Proof.** Every actual or finite-quotient representation modulo \(L_0\)
gives the necessary residue condition. For sufficiency, choose one
tuple \(u_i\in U_i\) with sum \(N\bmod L_0\). For \(i\ge3\), choose an
actual fixed \(z_i\in C_i(U_i)\) with \(z_i=u_i\bmod L_0\), using the
finite-image lemma. Set \(T=N-\sum_{i=3}^kz_i\).
The pair \((u_1,u_2)\) gives a positive binary density for \(T\).
The binary finite-congruence theorem therefore supplies all finite
congruences for the first two coordinates, proving the stated
all-moduli criterion. Its congruence argument uses only residue
arithmetic and also applies when \(T\) is nonpositive.
As \(N\) grows, \(T\) grows. The binary correlation theorem then gives
an actual positive pair with sum \(T\).
There are only finitely many residue tuples; fixing representatives
once and taking the maximum of their thresholds makes the conclusion
uniform in \(N\). Empty filters give empty residue sets and zero counts.
\(\square\)

For one subgroup \(H\), (G.1) states that the eventual \(k\)-fold sumset
has exactly the residue classes \(kH\). In particular,
\[
 kC_{12}\text{ eventually consists exactly of }N=k\bmod12
 \qquad(k\ge2\text{ fixed}).                            \tag{G.2}
\]
The qualifiers “eventually” and “fixed” are essential: this is not an
exact assertion about the initial exceptions or about \(k\) growing
with \(N\).

The binary instance is the strong Goldbach law for these semigroups.
The ternary instance is its weak Goldbach counterpart, retaining the
whole splitting congruence rather than parity alone.
Here the summands may be composite, and the unit value \(1\) belongs
to every core.

The proof just given is an existence corollary of the binary theorem.
It does not by itself give a natural-scale \(k\)-fold count: the
\(k-2\) anchors were fixed. Higher leading constants require the
separate correlation argument below.


## 8. Higher Goldbach counts: the classical finite-complexity range

For fixed \(k\ge3\) retain the data of (G.1). Let \(R_k(N)\) be the
number of ordered positive representations and put
\[
 \kappa_k=k-\sum_{i=1}^k\rho_i,\qquad
 \Gamma_k(N)=\frac{L_0}{\prod_i|H_i|}
       \#\{(u_1,\ldots,u_k)\in\prod_iU_i:\ \sum_i u_i=N\bmod L_0\}.
\]
Let \({\cal E}_k\) be the good primes forbidden on every side, and set
\[
 c_p(N)=\begin{cases}p-1,&p\mid N,\\-1,&p\nmid N,\end{cases}
 \qquad
 {\mathfrak S}_k(N)=
 \prod_{p\in{\cal E}_k}
       \left(1+\frac{(-1)^kc_p(N)}{(p-1)^k}\right).       \tag{G.3}
\]

**Theorem SplittingHigherGoldbachCount.** For every fixed \(k\ge3\),
uniformly in \(N\) and its target valuations,
\[
 R_k(N)=
 \frac{(\prod_i c_i)\Gamma_k(N){\mathfrak S}_k(N)}{(k-1)!}
       \frac{N^{k-1}}{(\log N)^{\kappa_k}}
 +o\!\left(\frac{N^{k-1}}{(\log N)^{\kappa_k}}\right).    \tag{G.4}
\]
A zero \(\Gamma_k\) gives zero count identically; otherwise the formula
is a relative asymptotic with a uniformly positive leading factor.
No logarithmic remainder exponent is asserted in this higher-fold
statement.

This is a classical finite-complexity consequence. Its proof uses
Matthiesen, *Linear correlations of multiplicative functions*,
[arXiv:1606.04482v4](https://arxiv.org/pdf/1606.04482v4),
Theorem 2.1, printed pp. 6–7, and Remark 2.2, together with
Loughran–Matthiesen, *Frobenian multiplicative functions and rational
points in fibrations*, [published primary edition](https://ems.press/content/serial-article-files/48228),
Proposition 3.3, Lemmas 3.4–3.5 and the fixed-character proof on
pp. 4799–4800.
It does not use the binary graph or 2026 zero-free inputs.
Matthiesen's Theorems 2.4 and 3.5 have a sublinear-intercept restriction
and are not the theorem applied here.

**Proof.** Write \(f_i=\mathbf1_{C_i}\).
These bounded nonnegative Frobenian multiplicative functions have
positive prime means \(\rho_i\). Proposition 3.3 puts them in
Matthiesen's class \({\cal F}\). The growth and nontrivial-prime-mass
conditions also follow directly from \(f_i(p^a)\in\{0,1\}\) and the
prime number theorem in fixed residue classes.

Residue filters need a separate check because they need not be
multiplicative. Expand each fixed filter on the unit group in Dirichlet
characters. Every term \(f_i\chi\), with fixed \(\chi\bmod L_0\), is in
\({\cal F}\). Indeed \(|f_i\chi|=f_i\). For the stable-progression-mean
condition, split a reduced class modulo \(q\) into classes modulo
\(\operatorname{lcm}(q,L_0)\). The ratio of these moduli is bounded.
Once \(q\) contains the source primorial they have the same prime
support; the normalizing mean \(E_{f_i}\) is unchanged. The character
is constant on each lifted class, so stability for \(f_i\) implies
stability for \(f_i\chi\). This verifies every fixed character term.

Apply Theorem 2.1 with \(r=k,s=k-1,B_0=0,T=N\), and forms
\[
 \varphi_i(x)=x_i\quad(i<k),\qquad
 \varphi_k(x)=N-\sum_{i<k}x_i.
\]
Their nonconstant parts are pairwise linearly independent when
\(k\ge3\). The moving intercept \(N\) meets that theorem's explicit
\(O(N)\) allowance at \(B_0=0\), and its conclusion is uniform in the
intercepts. The fixed open simplex
\[
 K=\{t_i>0,\ \sum_{i<k}t_i<1\}
\]
has volume \(1/(k-1)!\). On \(NK\) all form values are in \((0,N]\);
its integer points give exactly the positive ordered representations.
The source allows a convex positive-volume set and does not require
it to be closed. Thus no endpoint approximation is being assumed.

The printed theorem first gives an \(\varepsilon+o(1)\) error and a
finite \(W\)-model. Removing \(\varepsilon\) requires the additional
fixed-character argument. For \(f_i\), the characters of maximal prime
mean are precisely \(H_i^\perp\), and they satisfy
\(f_i(n)\eta(n)=f_i(n)\) identically. All their primitive conductors
divide fixed \(L_0\); a fixed twist just translates this family.
Choose the source function \(w(N)=\log\log N\), which meets the
progression hypotheses of Lemma 3.4. The argument on pp. 4799–4800
makes the \(W\)-trick independent of nilsequence complexity in exactly
this fixed-character situation, refining it to include \(L_0\).
It applies to Theorem 2.1 itself and preserves its moving-intercept
uniformity. We evaluate its finite model next; an arbitrary
existential \(W\) is not being assumed automatically to have the
needed refinement.

The principal-branch mean in a reduced class \(a\bmod q\), for source
moduli \(q\) divisible by \(L_0\), is
\[
 S_{f_i}(N;q,a)
 =c_i(\log N)^{\rho_i-1}\frac q{\varphi(q)}
   \left(\sum_{\eta\in H_i^\perp}\eta(a)\right)
   \prod_{\substack{p\mid q\\p\text{ allowed on side }i}}
                    (1-1/p)+o(E_{f_i}(N;q)).             \tag{G.5}
\]
Here \(S_f(N;q,a)=(q/N)\sum_{n\le N,\ n=a\bmod q}f(n)\).
Lemmas 3.4–3.5 give the uniform source-modulus error and the
coprimality product. Alternatively the principal coefficient follows
by twisting the Euler series (S.6). Character orthogonality makes
the displayed sum \(|G|/|H_i|\) on \(H_i\), zero elsewhere.
The errors are absolute at the \(E_{f_i}\) scale, including when this
branch vanishes. The product comparison uses only convergent
\(1+O(p^{-2})\) corrections, so these errors can be summed in the
finite model as in the cited fixed-character proof.

At the fixed modulus, a compatible root tuple has lattice density
\(L_0^{-(k-1)}\), and each coordinate state has marginal density
\(L_0/|H_i|\) relative to its unrestricted mean.
Their product gives exactly \(\Gamma_k\), retaining the joint
congruence states. At a good prime, an allowed coordinate has uniform
Haar distribution on \(\mathbb Z_p\); a forbidden coordinate has
density \(p/(p-1)\) on \(\mathbb Z_p^\times\).
If at least one coordinate is allowed, integrating that free uniform
coordinate gives local factor one.

If all coordinates forbid \(p\), character orthogonality gives
\[
 \#\{(a_1,\ldots,a_k)\in(\mathbb F_p^\times)^k:
                            \sum a_i=N\}
       =\frac{(p-1)^k+(-1)^kc_p(N)}p.
\]
Divide this count by \(p^{k-1}(1-1/p)^k\) to obtain (G.3).
There is no higher-valuation condition here: each coordinate need only
be a unit, and a unit solution lifts with uniform density. This also
explains uniformity in arbitrary target valuations.

Two tail checks justify passing from the source model to these full
local factors. First, for every \(p\ge3\),
\[
 \left|\frac{(-1)^kc_p(N)}{(p-1)^k}\right|
                  \le (p-1)^{1-k}.
\]
All factors are positive, with a uniform positive lower product.
Consequently the Euler product and its inverse are uniformly bounded,
and the omitted prime tail tends to zero uniformly in \(N\).

Second, the source model truncates small-prime prefixes at
\((\log N)^{B_2}\); this truncation can be removed uniformly.
In its normalized nonnegative local root model, an allowed coordinate
has
\(\mathbb E\,v_p(x_i)\ll_k1/(p-1)\).
If another coordinate is allowed this follows from Haar integration.
If it is the only allowed coordinate, the other \(k-1\ge2\) unit
coordinates give a mod-\(p\) probability \(O_k(1/p)\) for \(p\mid x_i\);
the higher digits are geometric. Forbidden coordinates have valuation
zero. Bad-prime prefixes vanish since every core is coprime to \(L_0\).
For each prefix \(w_i\), therefore,
\[
 \mathbb E\log w_i
   \ll_k\sum_{p\mid W}\frac{\log p}{p-1}
   \ll_{k,B_1}\log\log\log N.
\]
The last bound follows by splitting at \(z=\log\log N\):
the small-prime sum is \(O(\log z)\), and the rest is
\(O(\log W/z)=O(B_1)\), since \(W\le(\log N)^{B_1}\).
Markov's inequality bounds the lost normalized mass by
\[
 O_{k,B_1,B_2}\!\left(
                \frac{\log\log\log N}{\log\log N}\right)=o(1).
\]
For the finite Fourier expansion of the filters use the unfiltered
nonnegative root model to bound absolute tails. Thus the prefix
removal remains valid for every character term.

For each fixed \(\varepsilon>0\), the finite model is now the closed
expression in (G.4), independent of the auxiliary \(W,B_2\), up to
uniform \(o(1)\). The error scale in Theorem 2.1 is comparable to
\((\log N)^{-\kappa_k}\) by the fixed prime means.
Remark 2.2 consequently gives normalized limsup error \(O(\varepsilon)\).
Let \(\varepsilon\) tend to zero. The simplex volume supplies
\(1/(k-1)!\), proving (G.4).
Finally, \(\Gamma_k=0\) is an exact residue obstruction; its positive
values have a fixed positive minimum. The uniformly positive Euler
product makes the absolute asymptotic relative on those classes.
\(\square\)

For example, on \(N=3\bmod12\),
\[
 R_{12,3}(N)\sim
   6c_{12}^{\,3}\frac{N^2}{(\log N)^{9/4}}
   \prod_{\substack{p\in{\cal E}_{12}\\p\mid N}}
              \left(1-\frac1{(p-1)^2}\right)
   \prod_{\substack{p\in{\cal E}_{12}\\p\nmid N}}
              \left(1+\frac1{(p-1)^3}\right).             \tag{G.6}
\]
For every fixed \(k\ge3\), the corresponding \(C_{12}\) coefficient
is \(12c_{12}^k/(k-1)!\), with logarithmic exponent \(3k/4\).

The distinction between binary and higher correlation is structural:
for \(k=2\) the homogeneous parts \(x\) and \(-x\) are dependent,
and the source requires at least two ambient variables.
The classical higher-fold theorem therefore supplies neither the
binary leading asymptotic nor its target-uniform error.

## 9. Exact asymptotic basis orders

There are two different basis questions. A fixed number of positive
core values has a forced affine congruence. Adding zero allows the
number of positive summands to vary and gives an ordinary asymptotic
basis of all sufficiently large nonnegative integers. We compute
both orders from a finite cyclic group.

Fix \(H\le(\mathbb Z/L_0)^\times\), \(C=C_H\), and put
\[
 d=\gcd\bigl(L_0,\{h-1:h\in H\}\bigr),\qquad
 m=L_0/d,\qquad S=(H-1)/d\subset\mathbb Z/m.
\]
The gcd includes \(L_0\), so \(d=L_0\) when \(H=\{1\}\).
The definition of \(S\) is independent of representatives.
It contains zero, has size \(|H|\), and generates \(\mathbb Z/m\)
additively. Every core value is \(1\bmod d\). Since \(L_0\) is even,
\(d\ge2\).

**Theorem SplittingRelativeBasisOrder.** Let
\[
 t(H)=\min\{k\ge2:kS=\mathbb Z/m\}.
\]
This is exactly the smallest order \(k\ge2\) for which every sufficiently
large \(N=k\bmod d\) lies in \(kC\). Moreover
\[
 t(H)\le \max(2,m-|H|+1).                              \tag{G.7}
\]
For all fixed \(k\ge t(H)\), every sufficiently large \(N=k\bmod d\)
has a \(k\)-term representation.

**Proof.** Exactly
\[
                       kH=k+d(kS)\pmod{L_0}.
\]
Thus \(kH\) is the whole compatible affine class precisely when
\(kS=\mathbb Z/m\); (G.1) proves the equivalence with eventual
representability.

For the bound, let \(A_j=jS\). As \(0\in S\), these sets increase.
If \(A_j=A_{j+1}\), every translation by an element of \(S\) preserves
\(A_j\): it maps \(A_j\) into itself and preserves cardinality.
The additive group generated by \(S\) is the whole \(\mathbb Z/m\);
a nonempty set invariant under all those translations is the whole
group. Consequently each step before saturation adds at least one
element. Starting from \(|A_1|=|H|\), saturation occurs by
\(j=m-|H|+1\). Once full it stays full. Taking the minimum order at
least two proves (G.7). \(\square\)

This bound is sharp in an unbounded family. For even \(L_0\ge8\) and
\(H=\{1,-1\}\),
\[
 d=2,\quad m=L_0/2,\quad
 kH=\{k-2j:0\le j\le k\}\pmod{L_0},\quad t(H)=m-1.
\]
There is therefore no universal relative Goldbach order over all fixed
abelian subgroup prescriptions, even when parity is the only common
affine congruence. At smaller orders the additional obstructions remain
visible modulo \(L_0\).

For an ordinary basis, set \(C^0=C\cup\{0\}\). Its asymptotic basis
order is the least \(s\) such that every sufficiently large integer is
a sum of exactly \(s\) elements of \(C^0\), equivalently at most \(s\)
positive elements of \(C\).

**Theorem SplittingOrdinaryBasisOrder.** If \(H\) is proper, the exact
order is
\[
 b(H)=\min\left\{s\ge2:
                   \bigcup_{j=2}^s jH=\mathbb Z/L_0\right\},
 \qquad d+1\le b(H)\le L_0+1.                         \tag{G.8}
\]
If \(H=(\mathbb Z/L_0)^\times\), the order is two when \(L_0\) is a
power of two, and three otherwise.

**Proof.** Suppose first that \(H\) is proper. Choose a good prime
\(q\) in a residue class outside \(H\), by Dirichlet's theorem.
For every residue \(r\bmod L_0\), CRT gives infinitely many positive
integers \(N=r\bmod L_0\) divisible by \(q\). None belongs to \(C\).
Thus a single summand cannot supply eventual coverage of any whole
residue class.

If a residue \(r\) is absent from \(\bigcup_{j=2}^s jH\), infinitely
many such \(N\) cannot be represented by any positive number of
summands at most \(s\). This proves necessity in (G.8). Conversely,
if the union is full, choose a \(j\in[2,s]\) for each residue and
apply (G.1). There are finitely many residues and orders, so one
threshold suffices; pad with zeros to reach order \(s\).

All values in \(jH\) are \(j\bmod d\). For the class \(1\bmod d\),
no order \(2\le j\le d\) is possible, and the one-term obstruction
above proves \(b(H)\ge d+1\).
For the upper bound, \(1\in H\), so \(j\bmod L_0\in jH\).
The integers \(j=2,\ldots,L_0+1\) cover every residue.

For full \(H\), \(C\) is exactly the set of positive integers coprime
to \(L_0\). If \(L_0\) is a power of two these are all odd integers.
Large odd targets need one summand and large even targets need two;
the basis order is two. If \(L_0\) has an odd prime divisor, infinitely
many odd integers fail the one-term condition. Two odd summands cannot
represent them, so the order is at least three.
Modulo an odd prime power every residue is a sum of two units (avoid
at most two forbidden residues modulo the prime), while modulo a power
of two every even residue is a sum of two units. CRT gives
\(2H=\{\text{even residues}\}\) and hence
\(3H=\{\text{odd residues}\}\). Apply (G.1) to obtain order three.
\(\square\)

These are finite, exact formulae for the basis orders; the onsets in
the numerical semigroup sumsets remain unspecified.
They give a terminating finite algorithm: start from \(H\), compute
\(jH=(j-1)H+H\bmod L_0\), and track the cumulative union for
\(2\le j\le L_0+1\). For relative order, stop when \(|jH|=m\);
for ordinary order with proper \(H\), stop when the cumulative union
has size \(L_0\). The proofs above certify both stopping criteria.

| Core | Binary target classes | Relative order \(t(H)\) | Ordinary order of \(C\cup\{0\}\) |
|---|---|---:|---:|
| \(C_4\), \(H=\{1\}\bmod4\) | \(2\bmod4\) | 2 | 5 |
| \(C_3\), equivalently \(H=\{1\}\bmod6\) | \(2\bmod6\) | 2 | 7 |
| \(C_{12}\), \(H=\{1\}\bmod12\) | \(2\bmod12\) | 2 | 13 |
| \(H=\{1,11\}\bmod12\) | \(0,2,10\bmod12\) | 5 | 6 |
| All good primes, \(L_0=6\) | All even classes | 2 | 3 |

For the singleton rows, \(jH=\{j\bmod L_0\}\), and the least union
from \(j=2\) covering every residue ends at \(L_0+1\).
For \(H=\{1,11\}\), the formula \(jH=\{j-2r:0\le r\le j\}\)
shows that \(5H\) fills the odd classes, \(6H\) fills the even
classes, and \(6\bmod12\) is absent from the union through \(j=5\).
These direct residue calculations prove the table without an integer
search. For \(C_{12}\), minimality can also be seen explicitly:
\(25+60r\), \(r\ge0\), is \(1\bmod12\) and has forbidden prime factor
five. At most twelve positive summands would force just one summand
on this residue class. Conversely at most eleven copies of \(1\)
bring a large target to \(2\bmod12\), where two core values suffice.

In particular, for every fixed \(q\ge3\),
\[
 \boxed{\ \{0\}\cup C_q\text{ has exact asymptotic basis order }
                         \operatorname{lcm}(2,q)+1.\ }          \tag{G.9}
\]
This is a consequence of the correlation theory and residue arithmetic,
rather than a further analytic breakthrough.

## 10. What the Goldbach theory establishes

| Question | Result and proof source |
|---|---|
| Strong binary law | Explicit target-uniform leading count and logarithmic remainder, Sections 1–5 and Modules D–J |
| Weak and fixed higher laws | Explicit \(k\ge3\) asymptotic by the classical moving-target finite-complexity theorem, Section 8 |
| Local obstruction | Joint root states modulo \(L_0\), equivalent to all finite congruences; independent primewise projections can fail |
| Minimal summand count | Exact finite cyclic sumset formulae for relative and ordinary basis orders, Section 9 |
| Exact exceptions | Finitely many unspecified exceptions in admissible classes; numerical completeness needs an effective onset |
| General power-short Goldbach law | Not implied by these theorems; specialized proved ranges remain as recorded in the [fixed-modulus programme](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md#the-central-binary-conjecture) |
| Ramified or conductor generators | Finite profiles are covered; unbounded valuations can force a bounded residual target, as (S.10) shows |

All the existence and basis statements inherit the written-proof
assurance of their stated analytic route. The classical higher-fold
route does not require the imported 2026 binary inputs.
No historical novelty is claimed for these higher-fold consequences,
the finite-group reductions, or the Goldbach terminology.
Fixed prime-factor multiplicities, prime-only summands, growing
moduli and arbitrary nonabelian binary prescriptions are separate
problems with different hypotheses and counting scales.

The [coupled-system extension](COUPLED_SPLITTING_LOCAL_GLOBAL.md) treats
several linear equations with shared semigroup variables, including the
finite coefficient-prime tests, image-lattice normalization and an exact
45-class eventual criterion for shared primitive Gaussian norm sums.

The [multiplicative Hardy–Littlewood theory](MULTIPLICATIVE_HARDY_LITTLEWOOD.md)
extends the arithmetic model to arbitrary prime residue prescriptions and
Euler exponent rules, with a separate elementary theorem for dense sets.

The [actual-prime decompositions](PRIME_SPLITTING_DECOMPOSITIONS.md)
use separate prime-distribution inputs: full prime-star constants,
natural-order short-window counts, and a class-number-three application.
