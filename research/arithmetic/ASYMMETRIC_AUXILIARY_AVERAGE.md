# Averaging the auxiliary coefficient in the asymmetric convolution

The full auxiliary-modulus average in the dispersion estimate gives a
larger coefficient range than summing individual error terms would allow.
It controls a weighted sum of absolute correlation errors. Combined with
the exact finite-kernel identity, this also gives a positive level of
distribution for the asymmetric sequence. The critical square-root level
needed for the two-sided modulus-12 theorem remains open.

Include 1 in
\[
 C=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod{12}\},\qquad
 A=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod3\},
\]
and put \(U=\tau^2\mathbf1_C\), \(V=\tau\mathbf1_A\).
For \(N\equiv2\pmod{12}\), let
\(X_N=N\mathfrak S_{\rm asym}(N)\), using the exact positive main
term in [AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md).
For a positive integer \(c\) whose prime factors are all 7 modulo 12,
define
\[
 \begin{split}
 C_c(N)&=\sum_{\substack{a+ct=N\\a,t\ge1}}U(a)V(t),\\
 H(c)&=\prod_{p\mid c}(1-p^{-1}+p^{-2})^{-1},\\
 \mathcal M_c(N)&=
 \begin{cases}X_NH(c)/c,&(c,N)=1,\\0,&(c,N)>1.\end{cases}
 \end{split}
 \tag{1}
\]
The actual count is also zero when \((c,N)>1\): a common prime would
divide \(a\in C\). Arbitrary prime-power exponents in \(c\) and
\(N\) are allowed.

## Theorem AsymAuxiliaryAverage

For every fixed \(A_0\ge0\), there are constants \(\xi,\eta>0\)
such that, for every sufficiently large admissible target,
\[
 \boxed{\displaystyle
 \sum_{\substack{c\le N^\xi\\p\mid c\Rightarrow p\equiv7\ (12)}}
   \tau(c)^{A_0}|C_c(N)-\mathcal M_c(N)|
       \ll_{A_0}N^{1-\eta}.}
 \tag{2}
\]
The coefficient range \(\xi\) is chosen inside the auxiliary-modulus
and three-factor ranges before the smaller error saving. It is not
required to be smaller than \(\eta\); the proof below chooses
\(\eta\ll\xi\). In particular (2) controls any signed linear form
with a fixed divisor bound on its coefficients, uniformly even when
they depend on \(N\).

## Corollary AsymDistributionLevel

Fix \(\xi,\eta\) from (2) with \(A_0=1\). Use the exact density
and remainder in [AsymSemilinear](ASYMMETRIC_SEMILINEAR_SIEVE.md):
\[
 \begin{gathered}
 A_d(N)=\sum_{\substack{1\le b<N\\d\mid b}}U(N-b)V(b),\qquad
 r_d(N)=A_d(N)-X_Ng_N(d),\\
 g_N(p)=\begin{cases}
 (2p-1)/(p^2-p+1),&p\equiv7\ (12),\ p\nmid N,\\
 0,&\text{otherwise}.
 \end{cases}
 \end{gathered}
 \tag{3}
\]
Extend \(g_N\) multiplicatively on squarefree integers.
For every fixed \(0<\delta<\xi\), put
\[
 D=N^{\xi-\delta},\qquad
 \eta'=\min(\eta,\delta/2).
\]
Then
\[
 \sum_{\substack{d<D,\ d\ {\rm squarefree}\\
                  p\mid d\Rightarrow p\equiv7\ (12)}}
           |r_d(N)|\ll_\delta N^{1-\eta'}.
 \tag{4}
\]
The same estimate holds with any additional restriction
\(d\mid P_7(z)\), for arbitrary \(z\). Thus every fixed level exponent
strictly below \(\xi\) is available, rather than only \(\xi/2\).
This positive level remains below the square-root hypothesis in
[AsymCriticalCriterion](ASYMMETRIC_CRITICAL_SWITCHING.md).

## 1. The full dispersion input

The primary input is BGLM Proposition 5.1, printed pp.1192–1193 in
V. Blomer, L. Grimmelt, J. Li and S. Rydin Myerson,
*Additive problems with almost prime squares*, GAFA **33** (2023),
1173–1242: [DOI](https://doi.org/10.1007/s00039-023-00635-w),
[cached arXiv v5](sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf),
PDF p.15. Its auxiliary-modulus average, including the location of
its absolute value, is essential.

In its notation, for unit residues \(a\bmod q\), \(w\bmod d\),
\[
 u_R(n;q,a,d,w)=
 \mathbf1_{\substack{n\equiv a\ (q)\\n\equiv w\ (d)}}
 -\frac1{\varphi(q)\varphi(d)}
 \sum_{\substack{\chi\ (q),\ \psi\ (d)\\
                  \operatorname{cond}(\chi\psi)\le R}}
       \chi(n\bar a)\psi(n\bar w).
\]
For a sufficiently small fixed \(\eta_0>0\), the proposition supplies
\(\delta_0>0\). At product scale \(x=M_1M_2\), assume the short
factor scale satisfies
\[
 x^{\eta_0}\le M_2\le x^{1/4+\eta_0},\quad
 Q\le x^{1/2+\delta_0},\quad R\le x^{\delta_0},\quad
 |a_1|\le x^{1+\delta_0},\quad |a_2|\le x^{\delta_0}.
\]
For divisor-bounded \(\alpha_m,\beta_n\), and nonnegative
divisor-bounded \(\Lambda_d\), its specialization to the fixed
source modulus 12 is
\[
 \begin{split}
 &\sum_{d\le x^{\delta_0}}\Lambda_d
  \sup_{\substack{w\bmod d\ {\rm unit}\\a_3\mid d}}
  \left|
   \sum_{\substack{Q\le q\le2Q\\(q,a_1a_2d)=1\\
                   q\equiv b_0\ (12)}}
   \sum_{\substack{m,n\\(n,a_2)=1}}
    \alpha_m\beta_n
    u_R(mn;q,a_1\overline{a_2a_3},d,w)
  \right|\\
 &\hspace{35mm}\ll x(\log x)^{O(1)}R^{-1/3}.
 \end{split}
 \tag{5}
\]
Here \((b_0,12)=1\); the source's factor 12 is absorbed in the
constant. The condition \((n,a_2)=1\) is part of the theorem.
There is an absolute value for each auxiliary \(d\), outside its
\(q\)-sum. There is no sum of absolute values over \(q\), and
no arbitrary coefficient on the inner varying modulus.

The other analytic input is
[FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md), with
\(\sigma_0=1/10000\). Its error for one to three bounded
12-periodic interval factors of product scale \(Y\) is
\(\ll Y^{1-\sigma_0}/L\) for a unit residue, \((L,12)=1\), and
\(L\le Y^{1/2+\sigma_0}\).
The remaining character means and absolute kernel bounds are proved
in AsymPointwise; their precise uses are restated below.

## 2. Parameters and coefficient-independent families

Take \(\eta_0<\sigma_0/1000\) permitted by (5), and its associated
\(\delta_0\). Choose
\[
 \begin{gathered}
 0<\xi<\tfrac1{100}\min(\delta_0,\sigma_0,1/100),\\
 0<r<\tfrac1{10000}\min(\eta_0,\delta_0,\sigma_0,\xi),
 \qquad 0<\kappa<r/10000,\qquad 0<\nu<\kappa/10^6,\\
 R=N^r,\quad K=N^\kappa,\quad Z=N^\xi,\quad
 Q_{\mathrm{aux}}=N^{10\kappa},\quad
 \Delta=N^{-\zeta},\quad\zeta=\kappa/100.
 \end{gathered}
 \tag{6}
\]
The free divisor-bound exponent \(\nu\) can be reduced further to
absorb fixed logarithmic powers. Keep \(\xi\) fixed while decreasing
\(\kappa\) to ensure
\[
 \begin{gathered}
 \xi+11\kappa<(1-13\kappa)\delta_0,\\
 \frac{1+\xi+\kappa}{2}
 <(1-3\kappa)(1-8\eta_0)(1/2+\sigma_0),\\
 1<(1-13\kappa)(1+\delta_0),\quad
 \tfrac12<(1-13\kappa)(1/2+\delta_0),\quad
 r<(1-13\kappa)\delta_0,\\
 \eta_0(1-3\kappa)\ge10\kappa(1-\eta_0).
 \end{gathered}
 \tag{7}
\]
These inequalities are strict after the choices in (6). A final
\(\eta<\kappa/10000\) suffices. The choices are existential and
give no numerical cutoff.

For an allowed divisor \(g\) of the target define
\[
 F_g(n)=U(gn)/U(g),\qquad G_g(n)=V(gn)/V(g).
\]
The first sequence, with the target deletion \((n,M)=1\), is a
convolution of the four literal characters modulo 12 with a kernel
\(e_{g,M}\). The second is
\[
 G_g=e_g*g_{E,6},\qquad
 g_{E,6}(n)=(1*\chi_{-3})(n)\mathbf1_{(n,6)=1}.
 \tag{8}
\]
Both kernels have fixed divisor bounds and absolute \(3/5\)-norm
\(\ll_\nu N^\nu\). In the second kernel the only prime factors
are inert primes modulo 3, at square powers, and allowed primes of
\(g\). In particular no prime \(7\pmod{12}\) occurs there.
These families depend on \(g,M\), not on the coefficient \(c\).

## 3. Aggregate truncations and the exact hyperbola

For any \(b<N\),
\[
 \sum_{c\mid b}\tau(c)^{A_0}\ll_\nu N^\nu.
 \tag{9}
\]
Thus original summands \(a\le N^{1-\kappa}\) contribute
\(O(N^{1-\kappa+O(\nu)})\) to the weighted absolute family sum.
This counts actual divisors \(c\mid N-a\), rather than multiplying
an individual tail by \(Z\).

Stratify \(g=(a,N)\). Since \(a\in C\) and \((c,N)=1\),
\(g\in C\) divides \(N\) and \((c,g)=1\).
Discard \(g>K\) with the same bound using (9) and
\(\sum_{g\mid N,g>K}N/g\). For \(g\le K\), write
\[
 M=N/g,\qquad a=gm,\qquad t=gs,\qquad
 m+cs=M,\qquad(m,M)=1.
 \tag{10}
\]
The normalized weights are \(F_g(m)G_g(s)\), with external factor
\(U(g)V(g)\), independent of \(c\).

The first kernel tail can be bounded after recombining the \(c\)
divisors of each original value. For fixed \(g,M\), the corresponding
absolute expression is bounded by
\[
 \begin{split}
 &\sum_{m<M}\sum_{j\mid m,\ j>K}|e_{g,M}(j)|\tau_4(m/j)
     \sum_{c\mid M-m}\tau(c)^{O(1)}\tau((M-m)/c)^{O(1)}\\
 &\hspace{12mm}\ll N^{1+O(\nu)}
           \sum_{j>K}|e_{g,M}(j)|/j
   \ll N^{1-2\kappa/5+O(\nu)}.
 \end{split}
 \tag{11}
\]
The harmless external divisor-stratum factors are absorbed in
\(N^{O(\nu)}\). There is no assumption that an absolute truncated
convolution retains the original prime support.

Expand (8), truncate its index \(v\le K\), and discard all
\((v,M)>1\) terms before any source or model is applied. Such a term
has \(m=M-cvu\), hence a prime common to \(v,M\) would divide
\(m\), contradicting \((m,M)=1\). The second-kernel tail retains
\(N/(cv)\); summing \(c\) uses only
\(\sum_{c\le Z}\tau(c)^{A_0}/c\ll(\log N)^{O(1)}\).
It has the saving in (11).

The remaining residual argument of \(g_{E,6}\) is 1 modulo 6:
\(c\equiv1\pmod6\), \(M\equiv2\pmod6\),
\(m\equiv1\pmod6\), and every actual second-kernel index is 1
modulo 6. The strict divisor-pairing formula is
\[
 g_{E,6}(u)=2\sum_{q\mid u,\ q<\sqrt u}\chi_3(q)
              +\mathbf1_{u=\square}\chi_3(\sqrt u).
\]
Its aggregate square-boundary error is
\[
 \ll N^{1/2+O(\nu)}
       \sum_{c\le Z}\tau(c)^{O(1)}c^{-1/2}
       \sum_{v\le K}|e_g(v)|v^{-1/2}
 \ll N^{(1+\xi)/2+\kappa/10+O(\nu)}.
 \tag{12}
\]
The resulting hyperbola expression, at fixed \(g,v\), has the inner
sum for each \(c\)
\[
 2e_g(v)\mathbf1_{(v,M)=1}
 \sum_{\substack{q\le\sqrt{M/(cv)}\\(q,6M)=1}}\chi_3(q)
 \sum_{\substack{m\le M-cvq^2\\m\equiv M\ (cvq)}}
       F_g(m)\mathbf1_{(m,M)=1}.
 \tag{13}
\]
For (2), the absolute value is taken after this \(q\)-sum for each
\(c\). Triangle inequalities over the following fixed decompositions
will preserve that order.

## 4. Using the outer average without changing the two sequences

The first kernel in (13) depends only on \(g,M\). At index
\(j\le K\), a four-character product has scale
\(X\ge N^{1-3\kappa}\), since small original \(a\) was removed.
If a character variable has scale in
\[
 [X^{2\eta_0},X^{1/4+\eta_0}],
 \tag{14}
\]
use it as the short factor in (5). Combine the other three into the
arbitrary longer coefficient. If no such variable exists, use the
three-factor argument in Section 5.

The inverse multiplier is \(j\), not \(c\). Extract the entire part
\(c_0\) of the selected short variable supported on the primes of
\(j\), leaving a variable coprime to \(j\). Discard
\(c_0>Q_{\rm aux}\) by the smooth Rankin bound, and put the inverse
multiplier equal to \(jc_0\). The product scale is now
\(x\ge N^{1-13\kappa}\). The lower buffer in (14) and (7) keeps
the short factor in the source's original range. Its upper relative
exponent decreases on division by \(c_0\).
The two source sequences depend on \(g,j,c_0\) and the selected
character boxes, but are **independent of \(c\)**.

Write \(q=q_{\rm aux}q'\), where \(q_{\rm aux}\) is the entire part
of \(q\) supported on \(cv\), and truncate
\(q_{\rm aux}\le Q_{\rm aux}\). Fix \(g,v,j,c_0,q_{\rm aux}\)
and the fine boxes. Set
\[
 a_1=M,\qquad a_2=jc_0,\qquad a_3=1,\qquad
 d=cvq_{\rm aux},\qquad\text{fixed source modulus}=12.
 \tag{15}
\]
For these fixed scalars, \(c\mapsto d\) is injective, with inverse
\(c=d/(vq_{\rm aux})\). The restriction \((c,N)=1\), the required
support of \(q_{\rm aux}\), and the selected \(c\)-interval can all
be imposed in \(\Lambda_d\). They depend only on \(d\).
The weight \(\tau(c)^{A_0}\) is bounded by a fixed power of
\(\tau(d)\). The unit residue \(w\bmod d\) may depend on \(c\),
which is exactly what the supremum in (5) permits.

If \((j,cvq)>1\), the unit-residue congruence has no solution;
its unit model is also zero. Omit those terms. In every remaining
term \((q',a_1a_2d)=1\), and the reduced short factor is coprime to
\(a_2\). The latter follows because \(c_0\) is supported on
\(j\). The progression modulo \(cvq=dq'\) splits by CRT into
the two source congruences, with unit residue
\(M\overline{jc_0}\bmod d\). In particular the condition
\((v,M)=1\), together with \((c,N)=1\), ensures that this residue
is a unit.

The size conditions follow from
\[
 d\le N^{\xi+11\kappa}\le x^{\delta_0},\quad
 |a_2|\le N^{11\kappa},\quad M\le N,\quad
 q'\le N^{1/2},\quad x\ge N^{1-13\kappa},
\]
and (7). No inverse of \(c\) is placed in \(a_2\): it belongs to
the auxiliary modulus \(d\).
On each reduced \(q'\)-class modulo 12, \(\chi_3(q')\) is constant,
while \(\chi_3(q_{\rm aux})\) is fixed. Thus the inner modulus sum
in (5) has no unsupported phase.

The full outer average now bounds the weighted \(L^1\) over all
\(c\) in that block by a single error \(xR^{-1/3}(\log x)^{O(1)}\).
It is not multiplied by \(Z\). Summing the fixed scalar choices and
fine boxes costs a fixed \(N^{O(\kappa)}\); the deliberately generous
bound is
\[
 O\bigl(N^{1-r/3+100\kappa+O(\nu)}\bigr).
 \tag{16}
\]
A subinterval of the inner modulus range is handled by differences of
prefixes, each decomposed into \(O(\log N)\) intervals \([T,2T]\).
This costs only logarithms, not arbitrary \(q'\)-weights.

The discarded smooth tails retain their reciprocal weights. For example,
\[
 \sum_{\substack{h>Q_{\rm aux}\\h\mid(cv)^\infty}}\frac1h
 \ll_\nu N^{-5\kappa}(cv)^\nu.
\]
Its actual main contribution has \(1/c\), and hence aggregate bound
\(N^{1-5\kappa+O(\nu)}\). Its integer endpoints sum to at most
\(N^{(1+\xi)/2+O(\kappa+\nu)}\), using the \(c^{-1/2}\) sum.
For the \(c_0\)-tail, recombine divisors \(c\mid N-a\) first and
then count actual multiples of \(c_0\). Its Rankin saving is again
\(5\kappa\), with no number-of-coefficients factor. Nonprincipal
model tails use the four-character twist estimate and smooth Euler
sums at exponent \(3/4\), rather than \(R^2\) times these actual
error estimates.

## 5. Three-factor boxes, cumulative cuts, and model errors

If (14) contains no variable scale, at most three variables are long.
Fix the short variables, whose product is at most \(X^{8\eta_0}\).
The remaining product scale is at least
\(N^{(1-3\kappa)(1-8\eta_0)}\). The full modulus is
\[
 cvq\le\sqrt{Mcv}\le N^{(1+\xi+\kappa)/2}.
\]
Condition (7) puts it inside FixedPeriodicBoxes. Its individual error
has denominator \(cvq\), so the \(c\)-sum is harmonic. This gives
\[
 O\bigl(N^{1-\sigma_0+8\eta_0+O(\kappa+\nu)}\bigr).
 \tag{17}
\]
A fixed scalar nonunit at the full modulus gives zero actual and unit
model terms. No \(Z\) factor is introduced.

Use relative-width \(\Delta\) boxes in \(c,q'\), and the four
character variables. There are at most \(N^{6\zeta+\nu}\) boxes.
A \(c\)-box is selected only by \(\Lambda_d\) in (15), so the
other two source sequences stay fixed over that auxiliary average.
To count boundary boxes at \(m+cvq^2=M\), restore divisibility and
write \(M-m=cvqt\). Their near-boundary condition forces
\(\lvert cvq(t-q)\rvert\ll\Delta N\). Direct integer counting gives
\[
 O\left(\frac{\Delta N\log N}{cv}+\sqrt{\frac N{cv}}\right).
\]
At each original \(m\), recombine its first-kernel divisors with
\(\sum_{j\mid m}|e_{g,M}(j)|\tau_4(m/j)\ll_\nu N^\nu\).
Summing \(c\) and the second kernel gives aggregate boundary error
\[
 O\left(N^{O(\nu)}
  \left(\Delta N(\log N)^{O(1)}
       +N^{(1+\xi)/2}K^{1/10}\right)\right).
 \tag{18}
\]
The six-coordinate boxes form one partition; this boundary estimate
is not multiplied by their number. Principal model volumes retain
\(1/c\) and \(1/j\), and obey the same estimate.

For the nonprincipal truncated-character model, a primitive conductor
\(h\le R\) gives a literal four-character sum
\(\ll hN^{3/4+O(\nu)}\), with unit-removal divisor losses.
There are \(O(R^2)\) primitive characters, independent of \(c\).
The characters on \(q'\) and \(d\) have coprime moduli; their
primitive product determines the inflated pair uniquely. Their model
coefficients retain \(1/(\varphi(q')\varphi(d))\). Therefore
summing \(c\) costs only
\[
 \sum_{c\le Z}\frac{\tau(c)^{O(1)}}{\varphi(c)}
              \ll(\log N)^{O(1)}.
\]
The generous total budget is
\[
 O\bigl(R^6N^{3/4+100\kappa+6\zeta+O(\nu)}\bigr).
 \tag{19}
\]
All moduli are prime to 6, so a nonprincipal primitive conductor in
this model cannot divide 12. The elementary four-character twist
bound therefore applies. Unit removal and kernel convolution retain
\(j^{-3/4}\); smooth scalar sums retain their \(3/4\) weights.
No \(R^2\) multiple of an actual boundary length, or extra \(Z\),
is used.

## 6. Completion and exact normalization

Let \(H_{g,p}(z)=\sum_{i\ge0}F_g(p^i)z^i\), and set
\[
 \rho_{g,M}=\rho_{F_g}\prod_{p\mid M}H_{g,p}(1/p)^{-1},
 \qquad E_g(L)=\prod_{p\mid L}H_{g,p}(1/p)^{-1}.
\]
After the preceding estimates, the principal term at fixed \(g,v,c\)
is
\[
 2\rho_{g,M}e_g(v)\mathbf1_{(v,M)=1}
 \sum_{\substack{q\le\sqrt{M/(cv)}\\(q,6M)=1}}
  \chi_3(q)\frac{E_g(cvq)(M-cvq^2)}{\varphi(cvq)}.
 \tag{20}
\]
For \(W_p=[(1-1/p)H_{g,p}(1/p)]^{-1}=1+O(1/p)\), the
\(q\)-coefficients are \(\chi_3(q)\) times the corresponding
divisor weights. Their Dirichlet series is a fixed
\(L(s,\chi_3^{(6)})\) times a correction absolutely convergent for
\(\Re s>0\), with finite modifications at \(cvM\).
Its norms grow only by \(N^{O(\nu)}\). Convolution with bounded
partial sums of the fixed character gives a partial-sum bound
\(O((NcvQ)^{O(\nu)})\), also for growing \(c\).

Partial summation completes \(q\) to infinity and removes the
\(q^2\) term with error
\(O(N^{O(\nu)}\sqrt{M/(cv)})\) for each \(c,v\).
Its aggregate is (12). Completing the second kernel retains \(1/c\)
and has error \(N^{1-2\kappa/5+O(\nu)}\). The completed positive
main tail \(g>K\) retains the same harmonic \(c\)-sum and has
error \(N^{1-\kappa+O(\nu)}\).

At primes outside \(c\), the Euler computation is that of
AsymPointwise. If \(p^j\Vert c\), then \(p\equiv7\pmod{12}\)
and \(p\nmid M\). The first local sequence is unit-only, with
\(M_{U,p}=1-1/p\). The second kernel has no \(p\)-factor, by (8).
Since \(\chi_3(p)=1\), the completed local \(q\)-series is
\[
 \sum_{e\ge0}\frac1{\varphi(p^{j+e})}
                  =\frac{p^{-j}}{(1-p^{-1})^2}.
\]
After the first mean, this is \(p^{-j}M_{V,p}\), where
\(M_{V,p}=(1-p^{-1})^{-1}\). The base unit-target local factor is
\(I_{p,0}=(1-p^{-1}+p^{-2})/(1-p^{-1})\). Their ratio is
\[
 p^{-j}\frac{M_{V,p}}{I_{p,0}}
           =p^{-j}(1-p^{-1}+p^{-2})^{-1}.
 \tag{21}
\]
Multiplication over \(p\mid c\) gives exactly \(c^{-1}H(c)\),
for every exponent \(j\ge1\). At all other target primes the
completed \(g\)-strata restore the full all-valuation local factor.
The arithmetic factor 6 at 2 and 3 is unchanged. Thus the completed
main is precisely \(\mathcal M_c(N)\) in (1).

Every error above has a strict margin under (6)–(7). Absolute kernels
save \(2\kappa/5\), small summands and large gcd strata save
\(\kappa\), scalar Rankin tails save \(5\kappa\), and boundaries
save \(\zeta\). The dispersion saving is \(r/3-100\kappa\);
the three-factor saving is \(\sigma_0-8\eta_0-O(\kappa)\).
The model saving is \(1/4-6r-O(\kappa)\), while endpoints have
exponent \((1+\xi)/2+O(\kappa)<1\).
Choose a sufficiently small fixed \(\eta<\kappa/10000\).
The argument proves the weighted absolute family estimate (2), with
\(\xi\) fixed before that saving.

## 7. Extracting a positive distribution level

The exact compression identity (21) of AsymSemilinear says that, for
\(2\le L\le D\le N\) and coefficients \(|\lambda_d|\le1\)
supported on squarefree \(7\pmod{12}\)-prime divisors \(d<D\),
\[
 \begin{split}
 \sum_{d<D}\lambda_dr_d(N)
 ={}&\sum_{\substack{k\le L\\(k,N)=1}}\mu(k)
    \sum_{\substack{\ell<D/k\\(\ell,kN)=1}}
      2^{\omega(\ell)}\lambda_{k\ell}
       \bigl(C_{k^2\ell}(N)-\mathcal M_{k^2\ell}(N)\bigr)\\
 &+O_\varepsilon(N^{1+\varepsilon}/L).
 \end{split}
 \tag{22}
\]
Both \(k,\ell\) are squarefree with all prime factors 7 modulo 12.
The coefficients may be restricted further by \(d\mid P_7(z)\).
The cutoff \(L\) in this section is independent of the internal
analytic cutoff \(K=N^\kappa\) used to prove (2).

Take \(L=N^\delta\) and \(D=N^{\xi-\delta}\).
Every retained coefficient satisfies
\[
 c=k^2\ell=(k\ell)k<DL=N^\xi.
\]
Moreover the map \((k,\ell)\mapsto c\) is injective: the primes
of exponent 2 in \(c\) belong exactly to \(k\), and those of
exponent 1 belong exactly to \(\ell\). Thus the main sum in (22)
is one linear form in the errors from (2), with coefficient
\[
 a_c=\mu(k)2^{\omega(\ell)}\lambda_{k\ell},\qquad
 |a_c|\le2^{\omega(\ell)}\le\tau(c).
 \tag{23}
\]
Apply (2) with \(A_0=1\). The main sum is
\(O(N^{1-\eta})\).

If \(L\le D\), choose the free exponent
\(0<\varepsilon<\delta/2\) in (22); its tail is
\(O(N^{1-\delta/2})\). If \(L>D\), use the exact full identity
before truncation: every divisor \(k\mid d<D\) already satisfies
\(k<L\), so the tail is identically zero. This handles all
\(0<\delta<\xi\) without invoking (22) outside its stated
cutoff hypothesis. Terms with \((c,N)>1\) are zero throughout.

Finally choose \(\lambda_d=\operatorname{sgn}(r_d(N))\), with
\(\operatorname{sgn}(0)=0\), on the desired squarefree support.
The left side of (22) becomes its sum of absolute remainders.
The two bounds give exactly (4), with
\(\eta'=\min(\eta,\delta/2)\). This proves the corollary.

## Dependencies and scope

The proof uses the four-character means and kernels in AsymPointwise,
BGLM Proposition 5.1 with its full outer auxiliary-modulus average,
FixedPeriodicBoxes, and the exact finite-kernel identity in
AsymSemilinear. The input dependency is
\[
 \begin{gathered}
 \text{AsymPointwise character/kernel lemmas}
 +\text{BGLM full auxiliary average}
 +\text{FixedPeriodicBoxes}\\
 \Downarrow\quad\text{AsymAuxiliaryAverage},\\
 \text{AsymAuxiliaryAverage}
 +\text{AsymSemilinear exact kernel identity}
 \quad\Longrightarrow\quad\text{AsymDistributionLevel}.
 \end{gathered}
\]
The proof does not sum a previously known pointwise error over all
coefficients, and it introduces no arbitrary phase into BGLM's inner
modulus sum. Its range is limited by the auxiliary-modulus cap and
the full-modulus condition for the three-factor boxes.
The square-root remainder in AsymCriticalCriterion remains open.
The separate [subgroup correlation theorem](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md)
proves unweighted \(C+C\) counts for every sufficiently large admissible
target, including fixed macroscopic balance, under the retained ZF and G
imports. A numerical onset, exact finite coverage and historical priority
remain separate questions.
