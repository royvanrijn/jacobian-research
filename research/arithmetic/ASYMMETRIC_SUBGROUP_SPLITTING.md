# Asymmetric sums from abelian splitting classes

A subgroup of the reduced residue classes prescribes complete splitting
in an abelian extension. Pairing its prime-factor semigroup with the
splitting semigroup of a quadratic character gives a weighted asymptotic
for every sufficiently large admissible target. A finite residue density
determines admissibility; the two allowed prime sets need not be nested.
Removing the weights gives a count at the corresponding logarithmic scale.

The analytic inputs are
[BGLM Proposition5.1](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md) and
[Lemma FixedPeriodBoxesGeneral](ASYMMETRIC_FIXED_MODULUS_THEORY.md).
The fixed periodic filter, its character expansion and its interaction
with the principal main term are proved below. All data \(L,H,\chi\)
are fixed independently of the target.


Fix an even \(L\), let \(G=(\mathbb Z/L)^*\), and let \(H<G\) be a
proper subgroup of index \(\kappa=[G:H]\ge2\). Choose a fixed nonprincipal
quadratic literal character \(\chi\bmod L\), and write \(K=\ker\chi\).
No inclusion between \(H\) and \(K\) is assumed.
Define
\[
 \mathcal C_H=\{n\ge1:p\mid n\Rightarrow\bar p\in H\},\qquad
 \mathcal A_\chi=\{n\ge1:p\mid n\Rightarrow\chi(p)=1\}.
\]
All primes dividing \(L\) are excluded from both sets. Let
\[
 U_H(p^e)=
 \begin{cases}\binom{e+\kappa-1}{\kappa-1},&\bar p\in H,\\
 \mathbf1_{e=0},&\text{otherwise},
 \end{cases}
 \qquad V_\chi(n)=\tau(n)\mathbf1_{\mathcal A_\chi}(n).
\]

For a residue \(t\bmod L\), put
\[
 \Gamma(t)=\frac{L}{|H||K|}
       \#\{(a,b)\in H\times K:a+b=t\pmod L\}.
 \tag{H1}
\]
Call an integer admissible when \(\Gamma(N)>0\).
All admissible targets are even, but evenness alone need not suffice.

## 1. Theorems

Put \(F_B(m)=U_H(Bm)/U_H(B)\) and
\(G_C(l)=V_\chi(Cl)/V_\chi(C)\).
**Theorem SubgroupAsymUniform.** There exist fixed \(\eta,\theta>0\),
with \(\theta\) chosen after and
much smaller than \(\eta\), such that
\[
 \sum_{\substack{um+vl=N\\m,l\ge1}}F_B(m)G_C(l)
 =\frac{N}{uv}\mathfrak S_{B,C;u,v}(N)+O_{L,H,\chi}(N^{1-\eta})
 \tag{H2}
\]
for all sufficiently large positive targets, uniformly in
\[
 B,u\in\mathcal C_H,\qquad C,v\in\mathcal A_\chi,\qquad
 \max(B,C,u,v)\le N^\theta,
\]
all coefficient gcds and all target valuations.
When \(\Gamma(N)=0\), the count is exactly zero.
The main is
\[
 \mathfrak S_{B,C;u,v}(N)=
 \Gamma(N)\rho_{F_B}\rho_{G_C}
       \prod_{p\nmid L}\frac{D_p}{M_{A,p}M_{B,p}}.
 \tag{H3}
\]
Here \(A_e=F_B(p^e)\), \(B_e=G_C(p^e)\),
\(M_A=(1-p^{-1})\sum_eA_ep^{-e}\), similarly \(M_B\), and
\[
 D_p=\lim_{h\to\infty}p^{-h}
 \sum_{\substack{x,y\bmod p^h\\ux+vy=N\bmod p^h}}
 A_{\min(v_p(x),h)}B_{\min(v_p(y),h)}.
 \tag{H4}
\]
The usual truncated valuation at the zero residue is used.
Here \(\rho_F=\operatorname*{Res}_{s=1}\sum_{n\ge1}F(n)n^{-s}\).
The product converges absolutely; all local obstruction factors and
coefficient-gcd gains are included.

**Theorem SubgroupAsymWeighted.** For every sufficiently large positive
integer \(N\),
\[
 \sum_{a+b=N}U_H(a)V_\chi(b)
 =N\mathfrak S_{H,\chi}(N)+O_{L,H,\chi}(N^{1-\eta}),
 \qquad
 \mathfrak S_{H,\chi}=\mathfrak S_{1,1;1,1}.
 \tag{H4a}
\]
The main is positive when \(\Gamma(N)>0\). If \(\Gamma(N)=0\), both
the count and main are exactly zero. This is the specialization
\(B=C=u=v=1\) of SubgroupAsymUniform.

**Theorem SubgroupSemigroupScale.** One has

\[
 \#\{a+b=N:a\in\mathcal C_H,\ b\in\mathcal A_\chi\}
 \asymp_{L,H,\chi}
 \frac{N J_{H,\chi}(N)}{(\log N)^{3/2-1/\kappa}}
 \quad\text{for every sufficiently large admissible }N,
 \tag{H5}
\]
where
\[
 J_{H,\chi}(N)=
 \prod_{\substack{p\mid N,\ p\nmid L\\\bar p\notin H,\ \chi(p)=-1}}
                         \frac{p-1}{p-2}.
 \tag{H6}
\]
Thus this product is over primes outside both allowed sets.
No target prime merely outside \(K\) may automatically be put in it.

**Corollary SubgroupDivisibility.** When \(b\in\mathcal C_H\),
\(c\in\mathcal A_\chi\), the exact substitutions \(m=bx,l=cy\) give
\[
\begin{split}
 \sum_{\substack{um+vl=N\\b\mid m,\ c\mid l}}F_B(m)G_C(l)
 ={}&F_B(b)G_C(c)\frac{N}{ubvc}
          \mathfrak S_{Bb,Cc;ub,vc}(N)+O(N^{1-\eta}).
\end{split}
 \tag{H6a}
\]
After decreasing \(\theta\), this holds when each of
\(B,C,u,v,b,c\) is at most \(N^{\theta/10}\).
A forbidden divisor gives zero. Powers and exact valuations follow by
finite differences provided every tested product stays in range.
These are small-power tests, not a large auxiliary-modulus estimate.


## 2. Character factorization and residue measures

Let \(H^\perp\) be the \(\kappa\)-element group of characters trivial on
\(H\). For \(p\nmid L\), let \(f_p\) be the order of its coset in \(G/H\).
Then
\[
 D_{U_H}(s)=E_H(s)\prod_{\eta\in H^\perp}L(s,\eta),
 \qquad
 E_{H,p}(s)=
 \begin{cases}
 1,&f_p=1,\\
 (1-p^{-f_ps})^{\kappa/f_p},&f_p>1 .
 \end{cases}
 \tag{H7}
\]
This follows because the \(\eta(p)\) comprise all \(f_p\)-th roots
\(\kappa/f_p\) times. Since \(f_p\ge2\) at every corrected prime,
the correction is absolutely convergent for \(\Re s>1/2\).
All first kernel indices lie in \(H\) modulo \(L\), rather than
necessarily being 1. A term \(p^{f_p}\) is in \(H\).
The shifted and Euler-deleted corrections have the same bounded
polynomial degree and coefficient properties as in the identity-class
proof, now with degree depending on \(\kappa\).

The second factorization is exactly
\[
 D_{V_\chi}(s)=E_\chi(s)L(s,\chi_0^L)L(s,\chi),\qquad
 E_\chi(s)=\prod_{\chi(p)=-1}(1-p^{-2s}).
\]
Every second kernel index is an \(L\)-unit with \(\chi(k)=1\).
For all shifted/deleted kernels,
\[
 \sum_j|e(j)|j^{-3/5}\ll_{L,H,\chi,\epsilon}(BCT)^\epsilon,
 \qquad |e(j)|\le\tau(j)^{O_{L,H}(1)}.
 \tag{H8}
\]
First shifts are at primes in \(H\); second shifts at primes in \(K\).
The local coefficients are bounded uniformly in every shift valuation.

Put \(\phi=\varphi(L)\), and
\[
 \beta_\kappa=\max(3/5,1-1/\kappa),\qquad
 \delta_\kappa=1-\beta_\kappa>0.
\]
For any fixed character \(\psi\bmod L\), twisting the first series gives
the \(\kappa\) character factors \(\psi H^\perp\), with kernel
\(e(j)\psi(j)\), whose absolute norm is unchanged.
If \(\psi\in H^\perp\), then \(F_B(n)\psi(n)=F_B(n)\) pointwise.
If \(\psi\notin H^\perp\), all these factors are nonprincipal and
\[
 \sum_{n\le Y}F_B(n)\mathbf1_{(n,T)=1}\psi(n)
       \ll_{L,H,\epsilon}(BTY)^\epsilon Y^{\beta_\kappa}.
 \tag{H9}
\]
This is the largest-variable argument followed by (H8).
If an additional primitive character has conductor \(h\nmid L\),
none of the factors can become principal, for any fixed \(\psi\);
the same bound has an additional factor \(O_L(h)\).
In the BGLM model, its nonprincipal conductors are prime to \(L\),
so all such conductors satisfy this condition.

The untwisted first mean has error \(Y^{\beta_\kappa+\epsilon}\),
from hyperbola splitting the principal character from the other
\(\kappa-1\) characters.
This proof also covers ODD \(\kappa\), including \(\kappa=3\).
Conjugate characters occur in pairs, and all real nonprincipal
\(L\)-values at 1 are positive. Hence both residue means are positive.

Orthogonality in (H9) shows that the normalized residue measures are
uniform on \(H\) and \(K\), respectively.
The principal-factor residues are
\[
 \rho_{F_B}=\frac{\phi}{L}\prod_{p\nmid L}M_{A,p},\qquad
 \rho_{G_C}=\frac{\phi}{L}\prod_{p\nmid L}M_{B,p}.
 \tag{H10}
\]
The prefactor is \(\phi/L\), rather than \(\kappa/L\), unless \(H=\{1\}\).

## 3. The necessary fixed periodic filter

A nonzero second weight requires \((l,L)=1\) and \(\chi(l)=1\).
Before expanding its kernel and pairing divisors, insert
\[
 w_{M,u}(x)=
 \mathbf1_{\,(M-ux,L)=1,\ \chi(M-ux)=1}.
 \tag{H11}
\]
This does not change the convolution. It is independent of the second
kernel \(k\) and of \(v\), because both are \(L\)-units and
\(\chi(k)=\chi(v)=1\).
For every surviving term,
\[
 t=(M-ux)/(vk)\quad\text{has}\quad (t,L)=1,\ \chi(t)=1.
\]
Thus the exact symmetric divisor pairing applies. Without (H11), this
would fail when \(H\not\subset K\).

On the first unit variable expand the periodic filter, not the
second divisor weight:
\[
 w_{M,u}(x)=\sum_{\psi\bmod L}\widehat w_{M,u}(\psi)\psi(x),
 \quad
 \widehat w(\psi)=\phi^{-1}\sum_{r\in G}w(r)\overline{\psi(r)}.
 \tag{H12}
\]
There are only \(\phi\) terms, with coefficients of absolute value at
most one. Each term retains \(\kappa\) first character factors,
now \(\psi H^\perp\), and an absolutely summable twisted kernel.
After first kernel \(j\) and the coprimality scalar \(c_0\), the
additional factor is just the fixed scalar \(\psi(jc_0)\).
It is not a varying-modulus phase.

All fixed twists \(\psi\notin H^\perp\) have zero principal mean by
(H9); their principal and nonprincipal BGLM model terms have power-size
bounds. For \(\psi\in H^\perp\), the first sequence is unchanged.
The coefficient of its principal main term is
\[
 \sum_{\psi\in H^\perp}\widehat w_{M,u}(\psi)
 =\frac1{|H|}\sum_{a\in H}w_{M,u}(a)
 =:\mu(M).
 \tag{H13}
\]
Multiplication by the allowed coefficient \(u\in H\) permutes \(H\),
so this fraction is independent of \(u\). Also
\[
 \mu(M)=\frac{\#\{a\in H:M-a\in K\}}{|H|},\qquad
 \Gamma(M)=\frac{2L}{\phi}\mu(M).
 \tag{H14}
\]
This proves the fixed-place density, rather than assuming it survives a
new filter. Allowed coefficients \(u\in H\), \(v\in K\) permute
the two residue sets, so \(\Gamma(M)\) is independent of them.

No assertion is made for arbitrary coefficient cosets outside \(H,K\).
Those would have finite density for \(uH+vK\), rather than (H1).

## 4. Exact coefficient and gcd reductions

Write \(\mathcal C_I=\mathcal C_H\cap\mathcal A_\chi\), where
\(I=H\cap K\). This notation specifies prime support in the intersection,
not merely the residue of the whole integer.

For \(r=(u,v)\), every prime of \(r\) is in \(I\).
If \(r\nmid N\), there are no solutions. Otherwise divide coefficients
and target by \(r\).
For coprime coefficients let \(d=(u,N)\), \(e=(v,N)\).
There are two zero tests:
\[
 e\notin\mathcal C_H\Longrightarrow0,\qquad
 d\notin\mathcal A_\chi\Longrightarrow0.
 \tag{H15}
\]
The second is new relative to the nested case and is not automatic.
Indeed the equation forces \(e\mid m\), \(d\mid l\).
Since already \(e\in\mathcal A_\chi\), \(d\in\mathcal C_H\),
every surviving \(d,e\) has prime support in \(I\).

Substitute \(m=em'\), \(l=dl'\), with exact prefactors
\(F_B(e)G_C(d)\), new shifts \(Be,Cd\), coefficients \(u/d,v/e\),
and target \(N/(de)\).
All normalized weights remain in their allowed families and
\((u'v',N')=(u',v')=1\).
The prefactors are at most \(\tau_\kappa(e)\tau(d)\).

In the target-unit core, \(g=(m,N)\) divides \(l\) as well.
A nonzero stratum therefore has \(g\in\mathcal C_I\), not merely
\(g\in\mathcal C_H\).
Writing \(m=gx,l=gy,M=N/g\) gives the exact prefactor
\(F_B(g)G_C(g)\), shifts \(Bg,Cg\), and \((x,M)=1\).

For every \(s\in\mathcal C_I\),
\[
 \Gamma(N/s)=\Gamma(N)\quad\text{if }s\mid N.
 \tag{H16}
\]
Its residue belongs to \(H\cap K\), so multiplication by it permutes
both \(H\) and \(K\); this bijects the finite solution sets.
Thus all \(r,d,e,g\) preserve the finite density, although they need
not preserve the target's individual residue modulo \(L\).
Consequently \(\mu(M)\) in (H13) is constant throughout the gcd sum.

If \(\Gamma(N)=0\), reducing the original equation modulo \(L\)
gives an immediate exact contradiction. For \(\Gamma(N)>0\), every
divided target remains admissible by (H16).

## 5. Source parameters, all valuations and the principal main

Use the canonical fixed-period lemma at period \(L\) and BGLM with
progression modulus \(c=L\). These inputs do not require \(L\) to be
divisible by3, nor \(\kappa\) to be even.
Finite expansion (H12) only multiplies their constants by \(\phi\).

Choose \(\eta_0<\min\{\sigma_0/(1000\kappa),1/(4\kappa)\}\),
and its BGLM output \(\delta_0\). Choose
\[
 r<\min(\eta_0,\delta_0,\sigma_0,\delta_\kappa)/10000,\quad
 \lambda<r/10^6,\quad
 \zeta=\lambda/[100(\kappa+1)],
\]
and a core-family exponent \(\vartheta\) last and sufficiently small.
Put \(K_0=N^\lambda,H_0=N^{10\lambda},R=N^r,\Delta=N^{-\zeta}\).

Discard large gcd strata and small original summands, and truncate both
kernels at \(K_0\), exactly as in the fixed-modulus proof.
All weights and all new prefactors remain divisor-bounded.
The first character-product scale is at least
\(X=N^{1-3\lambda-\vartheta}\).
In every Fourier term, a medium variable in
\([X^{2\eta_0},X^{1/4+\eta_0}]\) is treated by BGLM after extracting
its entire \(uj\)-smooth part \(c_0\le H_0\).
The complement has Rankin saving \(N^{-5\lambda}\) up to divisor losses.
The resulting product scale is at least \(N^{1-13\lambda-\vartheta}\).

With no medium variable, at most three variables are long; their product
is at least \(X^{1-2\kappa\eta_0}\).
The fixed-period box lemma applies to the full modulus \(vka\) if
\[
 \tfrac12+\lambda+\vartheta
 <(1-3\lambda-\vartheta)(1-2\kappa\eta_0)(\tfrac12+\sigma_0).
 \tag{H17}
\]
Here \(a\) is the paired divisor variable, and \(k\) the second kernel.
It is an \(L\)-unit by (H11), as are \(k,v\).

Extract \(a=a_{\rm aux}a'\), where \(a_{\rm aux}\) is supported on
\(vk\), and truncate \(a_{\rm aux}\le H_0\) by Rankin.
The exact BGLM parameters remain
\[
 a_1=M,\quad a_2=ujc_0,\quad a_3=1,\quad
 d=vka_{\rm aux},\quad c=L.
\]
The auxiliary and varying moduli and inverse multipliers are \(L\)-units;
the residue \(M\) need NOT be an \(L\)-unit.
All source coprimalities follow after zero-solution terms are removed.
Selecting \(a'\) in a reduced class modulo \(L\) makes \(\chi(a')\)
constant. The filter coefficients and scalar \(\psi(jc_0)\) are fixed
in that application and do not compromise the signed modulus sum.

The size inequalities are exactly
\[
\begin{gathered}
 11\lambda+\vartheta<(1-13\lambda-\vartheta)\delta_0,\quad
 1<(1-13\lambda-\vartheta)(1+\delta_0),\\
 \tfrac12<(1-13\lambda-\vartheta)(\tfrac12+\delta_0),\quad
 r<(1-13\lambda-\vartheta)\delta_0,\\
 \eta_0(1-3\lambda-\vartheta)>10\lambda(1-\eta_0),\quad
 L\le N^{(1-13\lambda-\vartheta)\delta_0}.
\end{gathered}
 \tag{H18}
\]
They hold with the stated parameter order and eventually fixed \(L\).

The \(\kappa+1\) fine-box count is \(N^{\lambda/100+\epsilon}\).
Each Fourier term has the same actual boundary
\(ux+vka^2=M\), bounded by
\[
 N^\epsilon(\Delta N\log N+\sqrt N K_0^{1/10}).
\]
The actual multiplicity is
\(\sum_{j\mid x}|e(j)|\tau_\kappa(x/j)\ll N^\epsilon\).
The aggregate dispersion is
\(N^{1-r/3+100\lambda+100\vartheta+\epsilon}\).
All varying-conductor models and the fixed zero-mean twists from (H9)
cost at most
\[
 N^{\beta_\kappa+6r+100\lambda+100\vartheta+\epsilon}.
\]
Both save a power; the principal model retains \(j^{-1}\), the
zero-mean models \(j^{-\beta_\kappa}\). No conductor count multiplies
a boundary length or an actual Rankin tail.

Only the \(\psi\in H^\perp\) main terms survive, with factor \(\mu(M)\).
Their local radial profiles are exactly the original \(A_e,B_e\), because
\(F_B\psi=F_B\) pointwise. The standard divisor-series completion gives
the fixed-g raw coefficient
\[
 2(\phi/L)\mu(M)
\]
times the product of the same radial local factors as in AsymUniform:
\(1-p^{-1}\) at \(p\mid M\), \(M_A\) at \(p\mid u\), \(M_B\) at
\(p\mid v\), and \(I_{p,0}\) elsewhere.
The real prefactor is \(M F_B(g)G_C(g)/(uv)\).

For two unit coefficients the exact local integral is
\[
\begin{split}
 I_{p,n}={}&(1-x)\sum_{i<n}A_iB_ix^i+
 x^n(1-2x)A_nB_n\\
 &+(1-x)\left(B_n\sum_{i>n}A_ix^i+
                    A_n\sum_{j>n}B_jx^j\right),\quad x=1/p .
\end{split}
\]
The gcd choices only have primes in \(I\). This gives exactly the same
Euler reassembly, because \(A_iB_i=0\) for \(i\ge1\) at primes
outside \(I\).
At coefficient primes, the usual \(p^t\) common-coefficient gain and
shifted one-unit-coefficient formulas give (H4). Thus all target
valuations are included.

By (H10) and (H14),
\[
 \Gamma(N)\rho_{F_B}\rho_{G_C}
   =2(\phi/L)\mu(N)\prod_{p\nmid L}M_A M_B,
\]
so the preceding coefficient is exactly (H3).
Take a core saving below \(\zeta/1000\), then an external family power
smaller than both it and \(\vartheta/100\).
Coefficient reduction replaces the target by at least
\(N^{1-3\theta}\) and parameters by at most \((N')^{3\theta}\);
its divisor prefactors preserve a smaller saving.
This proves (H2) and the small-divisibility interface.

Primes dividing \(L\) need no assumption on their target valuations.
The finite density depends only on \(N\bmod L\), and all other
moduli are prime to \(L\). Weights at fixed places are periodic modulo
the corresponding powers in \(L\); higher target digits introduce
no additional line-density factor.
For example \(N\equiv0\pmod{14}\) is permitted in the test below even
with arbitrarily large \(v_7(N)\).

## 6. Unweighting for nonnested allowed sets

Put
\[
 a=|H|/|G|=1/\kappa,\qquad b=|K|/|G|=1/2,\qquad
 c=|H\cap K|/|G|.
\]
At a unit target, rough bad densities in the four prime classes are
\[
 \begin{array}{c|c}
 I&(\kappa+2)/p+O(p^{-2})\\
 H\setminus K&\kappa/p+O(p^{-2})\\
 K\setminus H&2/p+O(p^{-2})\\
 G\setminus(H\cup K)&0 .
 \end{array}
 \tag{H19}
\]
Their averaged dimension is
\[
 c(\kappa+2)+(a-c)\kappa+(b-c)2=\kappa a+2b=2.
\]
At positive target valuations, the bad density is
\(2\kappa/p+O(p^{-2})\) on \(I\) and zero in each other class.
The exact fully rough mass is always \(1-2/p\) or \(1-1/p\).
The finite small primes, including3 if it is an \(L\)-unit, are
treated exactly; their unit-target mass is at least1/3.

Use common cofactors \(d\in\mathcal C_I\), first directives
\(u\in\mathcal C_H\), second directives \(v\in\mathcal A_\chi\),
all squarefree and pairwise coprime, with \(d\mid N,(uv,N)=1\).
Normalize exact first valuation-one weights by
\(\kappa^{\omega(u)}\), second weights by \(2^{\omega(v)}\).
The prescribed-directive ratio relative to the fully rough unit model
is universally
\[
 a_p=(p-1)/(p(p-2)),\qquad E_p=(p-1)/(p-2).
\]
The exact local-predicate main includes the same finite factor
\(\Gamma(N)\): multiplication by \(u\in H\), \(v\in K\) permutes
the residue sets, while division by \(d\in I\) preserves it by (H16).

The finite-block lower sieve works with the densities in (H19).
Its support and error conditions are the same
\(2\alpha+\xi<(1-\alpha)\theta,\ \alpha+\xi<\eta/4\).
The exact fully rough main is the analogue of Q24 in the fixed-modulus
theorem, with positive uniformly bounded tail and with \(\Gamma(N)\)
in its fixed constant. Over the finitely many admissible classes,
\(\Gamma(N)\) has an absolute positive lower bound depending on the data.

At a unit target, restoration gives two \(a_p\) options on \(I\), one
first-side option on \(H\setminus K\), and one second-side option on
\(K\setminus H\). The total exponent restored is
\[
 2c+(a-c)+(b-c)=a+b=1/\kappa+1/2.
\]
At target primes in \(I\), the common-cofactor compensation is
\[
 E_p+p^{-1}E_{p,e-1}\ge1+2a_p.
\]
At target primes in either one-sided class,
\[
 E_p\ge1+a_p,\qquad E_p-(1+a_p)=1/[p(p-2)].
\]
At primes outside both classes, \(E_p\) remains as the target factor
(H6). This proves the target normalization in (H5), including first-only
allowed primes whose quadratic character is -1.

The actual positive cofactor model has inclusion probability \(O(1/p)\).
Markov on its logarithmic size retains a fixed fraction in the tiny
cofactor budget. The lift uses common \(d\) on VALUES only, with no
factor \(U_H(d)V_\chi(d)\).
Its recoverable patterns and normalized residual weight bound remain
\((1,1),(1,0),(0,1)\) and \((2\kappa)^{1/\gamma}\).
Thus the logarithmic lower exponent is
\(2-(a+b)=3/2-1/\kappa\).

The upper sieve excludes residue0 at primes outside \(H\), and residue
\(N\) at primes outside \(K\). Its dimension is
\[
 (1-a)+(1-b)=3/2-1/\kappa.
\]
Only primes outside BOTH sets merge two forbidden residues when they
divide the target. Their exact enhancement is \(E_p\), giving (H6).
The fixed-power upper block sieve and CRT endpoint estimate give the
matching bound, proving (H5).

## 7. Corrected odd-index example and negative controls

For \(L=14\),
\[
 G=\{1,3,5,9,11,13\},\quad H=\{1,13\},\quad
 K=\{1,9,11\},\quad\kappa=3
\]
with the literal \(\chi_{-7}\).
The input \(\{1,6\}\) is not a subgroup of \(G\), because6
is a nonunit; it is the corresponding subgroup modulo7.
This is a useful rejection control.

The six sums \(H+K\), with multiplicities, are
\[
 2,10,12,\ 0,8,10\pmod{14}.
\]
Hence
\[
 \Gamma(t)=
 \begin{cases}
 7/3,&t=0,2,8,12,\\
 14/3,&t=10,\\
 0,&\text{otherwise}.
 \end{cases}
\]
In particular \(N=4,6\bmod14\) gives exact zero, while \(N=0\bmod14\)
is admissible even though the target is not an \(L\)-unit.
At \(N=2\bmod14\), \(\mu(N)=1/2\), giving \(\Gamma=7/3\).
An allowed coefficient \(u=13\) swaps the contributing first residue,
but does not change this density.

The densities of \(I,H\setminus K,K\setminus H,\) and neither are
\(1/6,1/6,1/3,1/3\). Roughness has dimension2, restoration has
exponent5/6, and (H5) has exponent7/6.
The target product (H6) uses classes3 and5 modulo14, not class13.
A class13 target prime is first-only allowed and gets the one-option
compensation, not an inert two-exclusion target factor.

For an additional first-only prime3 control, take
\(L=16,H=\langle3\rangle=\{1,3,9,11\}\), \(\chi=\chi_{-4}\),
\(K=\{1,5,9,13\}\).
Then \(3\in H\setminus K\). Its target factor is compensated by
\(E_3\ge1+a_3\), and it must not appear in (H6), despite \(\chi(3)=-1\).

### Route comparison: \(C_3+C_4\) versus \(C_{12}+C_{12}\)

Write \(C_r=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod r\}\).
Take
\[
 L=12,\qquad H=\{1,5\},\qquad K=\{1,7\},\qquad\kappa=2.
\]
The first semigroup is \(C_4\), the second is \(C_3\).
Each sum of a residue in \(H\) and one in \(K\) occurs once, giving
\[
 \Gamma(N)=
 \begin{cases}3,&N\equiv0,2,6,8\pmod{12},\\0,&\text{otherwise}.
 \end{cases}
\]
Thus SubgroupSemigroupScale specializes to
\[
 \#\{a+b=N:a\in C_4,\ b\in C_3\}
 \asymp
 \frac{N}{\log N}
 \prod_{\substack{p\mid N\\p\equiv11\ (12)}}\frac{p-1}{p-2}
 \tag{R1}
\]
for every sufficiently large target in these four classes.
It does not by itself put both components in \(C_{12}\).

For \(N\equiv2\pmod{12}\), each component has residue1 modulo12,
but that is insufficient for the prime-factor restriction:
\[
 25+49=74,\qquad25=5^2\in C_4,\quad49=7^2\in C_3.
\]
Neither component belongs to \(C_{12}\).
The additional tasks are to exclude primes \(5\bmod12\) from the
first component and primes \(7\bmod12\) from the second.

At either type of prime, one side is always a unit and the permitted
side has local weight \(\tau(p^j)=j+1\). For \(x=1/p\) and \(p\nmid N\),
the total mass and the mass of divisibility on the permitted side are
\[
 I_0=\frac{1-x+x^2}{1-x},\qquad
 J=\frac{2x-x^2}{1-x}.
\]
The exact predicted event density is therefore
\[
 g_N(p)=
 \begin{cases}
 (2p-1)/(p^2-p+1),&p\nmid N,\\0,&p\mid N.
 \end{cases}
 \tag{R2}
\]
At a target prime the other, forbidden side is a unit and forces this
side to be a unit as well, proving the zero event.
Each family has density \(2/p+O(p^{-2})\) on one quarter of the primes,
hence sieve dimension \(1/2\); together they have dimension1.

Starting instead with \(C_{12}+C_3\) enforces the first of these two
restrictions before the final sieve, leaving only the second
half-dimensional family. This explains the distinction between the two
routes. It does not rule out a method using \(C_3+C_4\) and additional
joint distribution or switching estimates.


## 8. Scope of any almost-splitting consequence

The identity-class bounded-excluded-prime corollary does not retain
the exponent \(2-2/\kappa\) when \(H\not\subset K\).
Restricting second cofactors to \(I\) restores density \(a+c\),
not \(2a\). If \(\chi|_H\) is nonprincipal, then \(c=a/2\).
The unrestricted asymmetric result (H5) still restores \(a+b\).

No subgroup bounded-excluded-prime theorem is asserted here.
Statements (H2)--(H6) hold for all fixed admissible target
classes and all valuations, with odd or even subgroup index.
The abelian interpretation and precise dependencies are recorded below.
The full modulus-12 critical remainder is not resolved.


## 9. Abelian interpretation and primary scope comparison

Identify \(G\) with the Galois group of \(\mathbb Q(\zeta_L)/\mathbb Q\).
For \(F=\mathbb Q(\zeta_L)^H\), the primes in \(\mathcal C_H\) are
precisely its completely split primes away from primes dividing \(L\).
The quadratic fixed field of \(K\) gives the second splitting condition,
with the same finite exclusions. In particular the interpretation does
not silently include an extra prime dividing an auxiliary conductor.
Every nontrivial fixed abelian extension can be represented this way
after choosing a suitable common conductor and retaining these exclusions.
There is no nonabelian extension of the analytic argument here.

Drappeau--Topacogullari,
*Combinatorial identities and Titchmarsh's divisor problem for
multiplicative functions*, ANT13(2019),2383--2425,
[DOI](https://doi.org/10.2140/ant.2019.13.2383),
[primary article](https://msp.org/ant/2019/13-10/ant-v13-n10-p05-s.pdf),
Theorem1.2 on printedp2386, treats prime-periodic multiplicative functions
against \(\tau(an-h)\), for positive \(a\) and power-small shifts.
The same page explicitly separates the dual additive sum
\(\sum_{n<N}f(n)\tau(N-n)\) as outside its method.
It therefore does not directly supply the all-target family (H2).

Topacogullari,
*The shifted convolution of generalized divisor functions*, IMRN2018,
7681--7724, [DOI](https://doi.org/10.1093/imrn/rnx111),
[primary preprint](https://arxiv.org/pdf/1605.02364v1),
Theorems1.1--1.2 on preprintp2, gives power savings for each fixed
number \(k\ge4\) of divisor factors in a small-shift problem.
Those results concern ordinary divisor functions, not the additive
quadratic-character family and subgroup filter used here.
These comparisons distinguish the exact interfaces; they do not establish
a historical priority claim.

## 10. Dependencies and proof evidence

The exact BGLM statement, including its signed modulus sum,
large-residue range and short-variable coprimality condition, is recorded
with its [primary v5 PDF](sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf)
in [AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md).
[FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md) records
Heath-Brown's all-modulus complete correlation and three box estimates
with the [primary scan](sources/heath-brown-d3-progressions-1986.pdf).
Their arbitrary fixed-period transfer is proved in
[FixedModulusTheory](ASYMMETRIC_FIXED_MODULUS_THEORY.md).
The generic radial coefficient integrals and exact substitutions are
developed in [AsymUniform](ASYMMETRIC_UNIFORM_COEFFICIENTS.md);
the finite block sieve, cofactor cap and removal of weights are proved in
[AsymSemigroupScale](ASYMMETRIC_SEMIGROUP_SCALE.md).

The new subgroup obligations are the fixed filter (H11)--(H14), both zero
tests (H15), finite-density invariance (H16), and the four local classes
(H19). Their proofs are given here. Separate internal mathematical
reviews reconstructed these steps, including the zero-mean fixed-character
models and odd subgroup index. Finite exact controls test selected local
integrals, finite residue densities and invalid inputs. Such controls
do not verify the global analytic estimate. External review, formal
verification, publication and novelty are not asserted.

The [exact diagnostic checker](eisenstein-restricted-prime-sums/verify_subgroup_asymmetric.py)
and its [current receipt](eisenstein-restricted-prime-sums/subgroup_asymmetric_checks.json)
cover 13 nested and nonnested configurations, including odd subgroup index,
finite residue densities, coefficient and common-divisor obstructions,
144 local line integrals at 3, 5 and 7, and recoverable directive fibres.
They also check the mixed quadratic specialization (R1)--(R2).
The finite targets stop at 1800, with directive components at most 300;
these tests do not supply an eventual cutoff or prove the analytic estimates.

All constants concern fixed \(L,H,\chi\). There is no uniformity in a
growing conductor, effective cutoff, identified unweighted leading
constant, or fixed-ratio interval theorem. The full two-sided
modulus-12 critical remainder remains open.

