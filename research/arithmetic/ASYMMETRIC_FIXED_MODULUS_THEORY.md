# Asymmetric sums with a fixed congruence restriction and quadratic splitting

For any fixed congruence modulus, one summand may be restricted to primes
in its identity class while the other uses primes split by a chosen
quadratic character. A weighted pointwise asymptotic and a uniform
small-coefficient family give the unweighted count at its logarithmic
scale. Restricting the restoration further gives a count at the full
two-sided scale with only boundedly many excluded prime factors on the
second side.

Fix \(q\ge3\), put \(L=\operatorname{lcm}(2,q)\), and set
\(\kappa=\varphi(q)=\varphi(L)\).
Choose a fixed nonprincipal quadratic literal character \(\chi\bmod L\):
on units it has values \(+1,-1\), and it is zero at every nonunit.
Define
\[
 \mathcal C_q=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod q\}
              =\mathcal C_L,\qquad
 \mathcal A_\chi=\{n\ge1:p\mid n\Rightarrow\chi(p)=1\}.
\]
Include 1 in both semigroups. In particular all primes dividing \(L\)
are excluded from \(\mathcal A_\chi\), and
\(\mathcal C_q\subset\mathcal A_\chi\).
Let \(U_q,V_\chi\) be multiplicative, with
\[
 U_q(p^e)=
 \begin{cases}\binom{e+\kappa-1}{\kappa-1},&p\equiv1\pmod L,\\
 \mathbf1_{e=0},&\text{otherwise},\end{cases}
 \qquad V_\chi(n)=\tau(n)\mathbf1_{\mathcal A_\chi}(n).
 \tag{Q1}
\]

## 1. Theorems and exact scope

**Corollary FixedModulusAsymWeighted.** For every sufficiently large \(N\equiv2\pmod L\),
\[
 \sum_{a+b=N}U_q(a)V_\chi(b)
       =N\mathfrak S_{q,\chi}(N)+O_{q,\chi}(N^{1-\eta})
 \tag{Q2}
\]
for some fixed \(\eta>0\), with an everywhere positive local main term.
Here \(\mathfrak S_{q,\chi}=\mathfrak S_{1,1;1,1}\) from (Q8).

More generally, put
\[
 F_B(m)=U_q(Bm)/U_q(B),\qquad
 G_C(l)=V_\chi(Cl)/V_\chi(C).
\]
**Theorem FixedModulusAsymUniform.** There are fixed \(\eta,\theta>0\), with \(\theta\) chosen after and much
smaller than \(\eta\), such that
\[
 \sum_{\substack{um+vl=N\\m,l\ge1}}F_B(m)G_C(l)
 =\frac{N}{uv}\mathfrak S_{B,C;u,v}(N)
       +O_{q,\chi}(N^{1-\eta})
 \tag{Q3}
\]
uniformly for every sufficiently large \(N\equiv2\pmod L\) and
\(B,u\in\mathcal C_q\), \(C,v\in\mathcal A_\chi\),
each at most \(N^\theta\), and all target valuations.
The main is defined in Section2. Coefficient gcds are not excluded.
All admissible shifts and coefficients are automatically \(L\)-units.
A shift \(C\) with a prime dividing \(L\) is not admissible:
\(V_\chi(C)=0\), so its normalized ratio would be undefined.

Small divisibility tests are included:
\[
\begin{split}
 &\sum_{\substack{um+vl=N\\b\mid m,\ c\mid l}}F_B(m)G_C(l)\\
 &\qquad=F_B(b)G_C(c)\frac{N}{ubvc}
             \mathfrak S_{Bb,Cc;ub,vc}(N)+O(N^{1-\eta}),
\end{split}
 \tag{Q4}
\]
when \(b\in\mathcal C_q,c\in\mathcal A_\chi\), with every involved
parameter at a sufficiently smaller fixed power. A forbidden divisor
gives count zero. Powers and exact valuations are covered by finite
differences. The error is uniform enough for polynomially many tiny-power
tests, with the test exponent chosen smaller than \(\eta\).

**Theorem FixedModulusSemigroupScale.** For every sufficiently large
\(N\equiv2\pmod L\), one has
\[
 \boxed{\displaystyle
 \#\{a+b=N:a\in\mathcal C_q,\ b\in\mathcal A_\chi\}
 \asymp_{q,\chi}
 \frac{N J_\chi(N)}{(\log N)^{3/2-1/\kappa}},}
 \qquad N\equiv2\pmod L,
 \tag{Q5}
\]
where
\[
 J_\chi(N)=
 \prod_{\substack{p\mid N\\p\nmid L,\ \chi(p)=-1}}
                       \frac{p-1}{p-2}.
 \tag{Q6}
\]
Products in this note run over all \(p\nmid L\), potentially including 3.
The prime 2 is absent because \(L\) is even.
There is no imposed sign condition on \(\chi(-1)\).

## 2. The local main, including coefficient gcds

For \(p\nmid L\), set
\[
 A_e=F_B(p^e),\quad B_e=G_C(p^e),\qquad
 M_A=(1-p^{-1})\sum_eA_ep^{-e},\quad
 M_B=(1-p^{-1})\sum_eB_ep^{-e}.
\]
Define
\[
 D_p=\lim_{h\to\infty}p^{-h}
   \sum_{\substack{x,y\bmod p^h\\ux+vy\equiv N\ (p^h)}}
      A_{\min(v_p(x),h)}B_{\min(v_p(y),h)},
 \tag{Q7}
\]
with truncated valuation \(h\) at the zero residue. Then
\[
 \mathfrak S_{B,C;u,v}(N)=
 \frac{2L}{\kappa}\rho_{F_B}\rho_{G_C}
       \prod_{p\nmid L}\frac{D_p}{M_{A,p}M_{B,p}}.
 \tag{Q8}
\]
Here \(\rho_F=\operatorname*{Res}_{s=1}\sum_{n\ge1}F(n)n^{-s}\).
Outside the finitely many primes dividing \(BCuvN\), the factor is
\(1+O_{q,\chi}(p^{-2})\); hence the product converges absolutely.
The local limits exist by polynomial-geometric tail bounds.

For two unit coefficients and \(n=v_p(N)\), let \(x=p^{-1}\), \(q_p=1-x\).
The exact radial line integral is
\[
\begin{split}
 I_{p,n}(A,B)={}&q_p\sum_{i<n}A_iB_ix^i
       +x^n(1-2x)A_nB_n\\
 &+q_p\left(B_n\sum_{i>n}A_ix^i+
                  A_n\sum_{j>n}B_jx^j\right).
\end{split}
 \tag{Q9}
\]
This formula is valid at \(p=3\), where \(1-2/p=1/3>0\).
For general coefficient valuations \(a,b\), remove
\(t=\min(a,b)\). If \(v_p(N)<t\), then \(D_p=0\).
Otherwise multiply the residual answer by \(p^t\).
At most one residual coefficient has positive valuation. If
\(a>0,b=0\), the residual factor is
\[
 \begin{cases}
 B_n M_A,&n<a,\\
 B_a I_{p,n-a}(A,B^{[a]}),&n\ge a,
 \end{cases}
 \qquad B^{[a]}_j=B_{a+j}/B_a.
 \tag{Q10}
\]
There is the symmetric formula for \(a=0,b>0\).
A zero normalizing coefficient makes the factor zero.
These identities follow by solving for the unit-coefficient variable,
and by the normalized \(p^t\) gain on dividing the congruence.
Thus all target valuations and common coefficients occur in (Q8).

In the target-unit coefficient case,
\[
 (u,v)=(uv,N)=1,
\]
the factors simplify to
\[
 p\mid u:D_p=M_A,\qquad
 p\mid v:D_p=M_B,\qquad
 p\nmid uv:D_p=I_{p,v_p(N)}(A,B).
 \tag{Q11}
\]
In the coefficient-target-unit core, the unit-unit stratum has mass at
least \(1-2/p\) at a unit target and \(1-1/p\) at a divisible target.
These bounds are positive at every prime in the product. General
nonobstructed positivity follows by the exact coefficient reductions.

## 3. Dirichlet series and uniform shifted kernels

For \(p\nmid L\), let \(f_p\) be its order in \((\mathbb Z/L)^*\).
Character orthogonality gives
\[
 D_{U_q}(s)=E_q(s)\prod_{\psi\bmod L}L(s,\psi),\qquad
 E_{q,p}=
 \begin{cases}
 1,&f_p=1,\\
 (1-p^{-f_ps})^{\kappa/f_p},&f_p>1.
 \end{cases}
 \tag{Q12}
\]
Indeed the character values at \(p\) comprise every \(f_p\)-th root
\(\kappa/f_p\) times. All corrections have degree at least two in
\(p^{-s}\), so \(E_q\) is absolutely convergent for \(\Re s>1/2\).

For the second series,
\[
 D_{V_\chi}(s)=E_\chi(s)L(s,\chi_0^L)L(s,\chi),\qquad
 E_\chi(s)=\prod_{\substack{p\nmid L\\\chi(p)=-1}}(1-p^{-2s}).
 \tag{Q13}
\]
At \(\chi(p)=1\) the two-character local factor is \((1-z)^{-2}\);
at \(\chi(p)=-1\) it is \((1-z^2)^{-1}\), cancelled by the correction.
No local factor occurs at \(p\mid L\).

At a shifted first prime \(p^b\Vert B\), replace the first correction by
\[
 (1-z)^\kappa\sum_{e\ge0}
 \frac{\binom{b+e+\kappa-1}{\kappa-1}}
      {\binom{b+\kappa-1}{\kappa-1}}z^e.
\]
Its degree is at most \(\kappa-1\) and its coefficients are bounded
independently of \(b\). The ratio is a polynomial in \(e\) with uniformly
bounded coefficients, since every denominator \(b+j\ge j\).
At \(p^c\Vert C\), the second correction is
\[
 1-\frac c{c+1}z.
\]
First-side Euler deletion at \(T\) substitutes the polynomial cancelling
all first character factors. Consequently every relevant kernel satisfies
\[
 |e(j)|\le\tau(j)^{O_q(1)},\qquad
 \sum_j|e(j)|j^{-3/5}\ll_{q,\chi,\epsilon}(BCT)^\epsilon.
 \tag{Q14}
\]
The first kernel has indices \(1\bmod L\).
The second has indices coprime to \(L\) and satisfying \(\chi(j)=1\):
an inert local index is a square, while a shifted local prime has
\(\chi(p)=1\). This property replaces the special \(1\bmod6\) support;
no mod-3 support assertion is made.

Set
\[
 \beta_\kappa=\max(3/5,1-1/\kappa),\qquad
 \delta_\kappa=1-\beta_\kappa>0.
\]
For primitive conductor \(h\nmid L\), every first character twisted by
that primitive character is nonprincipal. Summing the largest of the
\(\kappa\) variables last, then using (Q14), gives
\[
 \sum_{n\le Y}F_B(n)\mathbf1_{(n,T)=1}\xi(n)
 \ll_{q,\chi,\epsilon}
     h(BTY)^\epsilon Y^{\beta_\kappa}.
 \tag{Q15}
\]
The exponent \(3/5\) is retained when \(\kappa=2\), so the square-root
box error is not incorrectly convolved against a .6 norm.
The untwisted mean has the same error, by hyperbola splitting the
principal factor from the other \(\kappa-1\) fixed characters.
In natural prime order,
\[
 \rho_{F_B}=\frac{\kappa}{L}\prod_{p\nmid L}M_A,\qquad
 \rho_{G_C}=\frac{\kappa}{L}\prod_{p\nmid L}M_B.
 \tag{Q16}
\]

These residues are positive for BOTH even and odd \(\chi\).
The real nonprincipal Dirichlet values at 1 are positive; complex first
characters occur in conjugate pairs. If \(\chi\) is imprimitive, its
value at 1 is its positive primitive value multiplied by positive real
factors \(1-\chi^*(p)/p\). Every real kernel correction at \(s=1\)
is positive. No use is made of \(\chi(-1)\).

## 4. The fixed-place factor

The first sequence is supported exactly on \(1\bmod L\).
For the second sequence, twist by an arbitrary character \(\psi\bmod L\).
A pole occurs only when \(\psi=\chi_0\) or \(\psi=\chi\):
these twists are identical to \(G_C\), since \(\chi(n)=1\) on its support.
For every other \(\psi\), both \(\psi\) and \(\psi\chi\) are
nonprincipal, and the two-character box argument plus (Q14) gives a
power error. Orthogonality therefore gives the normalized residue masses
\[
 \mu_F(a)=\mathbf1_{a=1\bmod L},\qquad
 \mu_G(b)=\begin{cases}2/\kappa,&(b,L)=1,\ \chi(b)=1,\\0,&\text{otherwise}.
             \end{cases}
\]
The kernel of the nonprincipal quadratic character has size \(\kappa/2\).
For \(u\equiv1\pmod L\), \(v\in\mathcal A_\chi\), and \(N\equiv2\pmod L\),
the first variable has residue1 and the second must have residue \(v^{-1}\),
which has \(\chi(v^{-1})=1\). Thus
\[
 L\sum_{ua+vb=N\bmod L}\mu_F(a)\mu_G(b)=2L/\kappa.
 \tag{Q17}
\]
This proves the fixed-place factor in (Q8), including for even \(\chi\).

## 5. Divisor boxes with any fixed period

**Lemma FixedPeriodBoxesGeneral.** Let \(c\) be any fixed positive integer.
Let \(\omega_1,\omega_2,\omega_3\) have period \(c\), with
\(|\omega_i|\le1\), and let their interval
factors have scales \(N_i\), product \(X\), as in
[FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md).
For \(\sigma_0=1/10000\), their convolution satisfies
\[
 \left|\sum_{n\equiv a\ (m)}f(n)
  -\frac1{\varphi(m)}\sum_{(n,m)=1}f(n)\right|
 \ll_c X^{1-\sigma_0}/m
 \tag{P1}
\]
when \((m,c)=(a,m)=1\), \(m\le X^{1/2+\sigma_0}\).
Singleton padding and bounded-variation weights are permitted.
In this article it is used with \(c=L\).

The following proves the fixed-period transfer, including the modulus
range. Its deep input is exactly Heath-Brown's all-modulus complete
correlation Lemma4, recorded with its statement and primary scan in
[FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md).

Put \(Q=cm\) and
\[
 W_i(h)=\sum_{n\in I_i}e_Q(-hn),\qquad
 G_i(t)=\sum_{r\bmod c}\omega_i(r)e_c(tr).
\]
Then \(|G_i|\le c\) and
\(|W_i(h)|\ll\min(\#I_i,Q/|h|)\).
CRT requires only \((c,m)=1\), and gives
\[
 S(a)=Q^{-3}\sum_{h_1,h_2,h_3\bmod Q}
  \prod_iW_i(h_i)G_i(m^{-1}h_i)
  K_2(c^{-1}h_1,c^{-1}h_2,ac^{-1}h_3;m).
 \tag{P2}
\]
Unit changes of variables replace the kernel by
\(K_2(h_1,h_2,a'h_3;m)\), \(a'=ac^{-3}\).
The centered frequency range is \(|h_i|\le cm/2\).
If a frequency is zero modulo \(m\), its Ramanujan-sum factor is
independent of \(a\); remove such terms into an \(a\)-independent main.

In the principal estimate, bound \(G_1,G_2\) before combining their
frequencies. Retain \(G_3\) in coefficients \(b_t\) with
\(|b_t|\ll_c\min(K,m/T)\).
Write \(r=\rho r'\), \(s=\sigma s'\), with \(r's'\) units and
\(\rho,\sigma\) supported on prime divisors of \(m\).
The identity
\[
 K_2(\rho r',\sigma s',a't;m)
     =K_2(\rho,\sigma,a'r's't;m)
\]
allows the product frequency \(h=r's'\).
Cauchy and triangle inequality then leave the unweighted correlation
\[
 \sum_{\substack{H<h\le4H\\(h,m)=1}}
 K_2(\rho,\sigma,a'ht_1;m)
 \overline{K_2(\rho,\sigma,a'ht_2;m)}.
 \tag{P3}
\]
Its completion is the original complete correlation, with unrestricted
additive parameter \(a'^{-1}j\). No fixed-period coefficient remains
inside it.

The extended ranges require exactly three controls. Smooth parts extend
only to \(cm\), costing \(O_{c,\epsilon}(m^\epsilon)\). Congruent
diagonal frequencies give \(O_c(1)\) partners per frequency. For other
pairs a divisor class in an interval of length \(O_c(T)\) has
\(O_c(1+T/d)\) representatives. Completion of (P3) is valid for every
\(H\), including \(H>m\); its zero Fourier mode has the existing
factor \(H/m\). Thus the divisor and diagonal calculations in the
linked proof retain their exponents.

For the small-variable estimate, bound the two completed Gauss factors
absolutely and apply the ordinary Weil--Estermann Kloosterman estimate.
For the long-variable estimate, bound \(G_3\) in positive collision
weights before Hölder. The weighted fourth moment is bounded by the
unweighted reciprocal energy, since every expanded coefficient has
absolute value at most one. The original three estimates therefore hold
with constants depending on \(c\):
\[
\begin{split}
 E_{\rm pr}&\ll_{c,\epsilon}X^\epsilon
     \{m^{5/6}+m^{1/4}(IJ)^{1/2}+m^{1/2}K^{1/2}\},\\
 E_{\rm small}&\ll_{c,\epsilon}X^\epsilon m^{1/2}I,\\
 E_{\rm long}&\ll_{c,\epsilon}X^\epsilon
 I^{3/4}J^{1/2}(K+U)^{1/4}
 (1+V^2/m+V^{3/2}/m^{1/2})^{1/4}.
\end{split}
 \tag{P4}
\]
Here \(U,V,W\) are interval scales, \(I,J,K\) their numbers of integers,
and variables may be permuted. These are the transferred versions of
Heath-Brown's Lemmas5,6,7, printed pp44,47.

For the exponent selection, order scales \(S\le M\le T\), \(SMT=X\).
Take \(\nu=1/100\), \(\delta=1/10000\),
\(\epsilon=1/1000\), and \(m\le X^{1/2+\delta}\).
If \(S\le X^{1/4-\nu}\), the small-variable estimate gives
\[
 mE_{\rm small}\ll X^{1-\nu+3\delta/2+\epsilon}
                    \ll X^{1-\delta}.
\]
If \(T\le X^{1/2-\nu}\), attach \(K\) to the largest interval in
the principal estimate. Its three terms, after multiplication by \(m\),
have exponents at most
\[
 \frac{11}6(\tfrac12+\delta),\qquad
 \frac54(\tfrac12+\delta)+\frac13,\qquad
 \frac32(\tfrac12+\delta)+\frac12(\tfrac12-\nu).
\]
They are \(0.91685\), less than \(0.95846\), and \(0.99515\).
Adding \(\epsilon=0.001\) leaves them below \(1-\delta=0.9999\).

Otherwise
\[
 S>X^{0.24},\quad T>X^{0.49},\quad
 SM<X^{0.51},\quad M<X^{0.27},\quad T<X^{0.52}.
\]
Use the long-variable estimate with \(U=S,V=M,W=T\).
Its three expanded terms are at most
\[
 X^{1/4}S^{1/2}M^{1/4},\quad
 X^{1/4}S^{1/2}M^{3/4}m^{-1/4},\quad
 X^{1/4}S^{1/2}M^{5/8}m^{-1/8}.
\]
The scale bounds give respective non-modulus exponents
\(0.44125,0.5725,0.53875\). After multiplication by \(m\),
the total exponents are at most
\[
 0.94135,\qquad0.947575,\qquad0.9763375.
\]
Each remains below \(0.9999\) after \(\epsilon\).
This proves the saving in (P1). Different main terms are harmless:
averaging an \(a\)-independent main over unit \(a\) changes the error
by at most a factor two. Bounded \(X\) is absorbed in the \(c\)-constant,
and \(m=1\) has zero discrepancy.
Singleton padding and tensor-product partial summation complete the proof.

In particular, if \(3\nmid L\), (P1) permits arbitrarily large powers
of 3 in its variable modulus. This uses the actual all-modulus
Heath-Brown input. No artificial coprimality to 3 or 6 remains.
BGLM likewise imposes no universal coprimality to 3; its fixed
progression modulus below is \(c=L\).


## 6. Exact arithmetic reduction and divisor pairing

For (Q3), first remove \(r=(u,v)\). If \(r\nmid N\), the count is zero;
otherwise divide both coefficients and the target by \(r\).
Here \(r\in\mathcal C_q\), so \(r\equiv1\pmod L\).
For coprime coefficients put \(d=(u,N)\), \(e=(v,N)\).
The equation forces \(e\mid m\), \(d\mid l\). If
\(e\notin\mathcal C_q\), the count is zero. Otherwise substitute
\(m=em'\), \(l=dl'\), giving exact prefactors
\[
 F_B(e)G_C(d)
\]
and new shifts \(Be,Cd\), coefficients \(u/d,v/e\), target \(N/(de)\).
They have \((u'v',N')=(u',v')=1\).
All divided target factors lie in \(\mathcal C_q\), so \(N'\equiv2\pmod L\).
The prefactors are bounded by \(\tau_\kappa(e)\tau(d)\).
The local transformation gives the same coefficient factors in (Q7).
For a common coefficient \(r\), the main expressed with the original
coefficients has scale \(rN/(uv)\).

In this target-unit case, partition \(g=(m,N)\in\mathcal C_q\).
Writing \(m=gx,l=gy,M=N/g\), the exact prefactor is
\(F_B(g)G_C(g)\), with shifts \(Bg,Cg\) and condition \((x,M)=1\).
After second kernel \(k\), define
\[
 t=(M-ux)/(vk).
\]
For every contributing first integer \(x\), one has \(ux\equiv1\pmod L\),
\(M\equiv2\pmod L\). Thus the numerator is \(1\bmod L\).
Since \(v,k\) are \(L\)-units with \(\chi(v)=\chi(k)=1\),
\[
 (t,L)=1,\qquad\chi(t)=1.
 \tag{Q18}
\]
This is the entire character-sign condition needed for divisor pairing.
It holds when \(\chi(-1)=1\), as well as when \(\chi(-1)=-1\).
Therefore
\[
 (\chi_0^L*\chi)(t)
 =2\sum_{\substack{a\mid t\\a<\sqrt t}}\chi(a)
       +\mathbf1_{t=\square}\chi(\sqrt t).
\]
Impossible common factors are removed exactly as in AsymUniform.
The progression modulus is \(vka\), its first inverse multiplier is \(uj\),
and the residue/deletion target is \(M\).
The restrictions are \((k,uM)=1\), \((a,LuM)=1\), with NO extra
factor3 in this last restriction.

## 7. Every source condition and power budget

Choose a permitted BGLM exponent \(\eta_0\) with
\[
 \eta_0<\min\{\sigma_0/(1000\kappa),1/(4\kappa)\},
\]
and let \(\delta_0\) be its source output. Choose
\[
 r<\min(\eta_0,\delta_0,\sigma_0,\delta_\kappa)/10000,
\]
then sufficiently small \(\lambda>0\), and finally a core-family power
\(\vartheta>0\). Put
\[
 K=N^\lambda,\quad H=N^{10\lambda},\quad R=N^r,\quad
 \zeta=\lambda/[100(\kappa+1)],\quad \Delta=N^{-\zeta}.
\]
Large gcd strata and small original summands cost \(N^{1-\lambda+\epsilon}\).
Both kernels truncate at cost \(N^{1-2\lambda/5+\epsilon}\).
After discarding \(um\le N^{1-\lambda}\), the first character product
has scale \(X\ge N^{1-3\lambda-\vartheta}\).

If a variable lies in \([X^{2\eta_0},X^{1/4+\eta_0}]\), combine the
other \(\kappa-1\) factors into an arbitrary divisor-bounded coefficient.
Extract its full \(uj\)-smooth part \(c_0\); the part \(c_0>H\) costs
\(N^{1-5\lambda+\epsilon}\) by Rankin's inequality. For \(c_0\le H\),
the new short variable is coprime to \(ujc_0\), and the product scale is
at least \(N^{1-13\lambda-\vartheta}\).

If no such scale exists, there are at most three long variables and the
product of all short scales is at most \(X^{2\kappa\eta_0}\).
Apply P(L) to the one to three long factors, at full modulus \(vka\).
Its sufficient inequality is
\[
 \tfrac12+\lambda+\vartheta
 <(1-3\lambda-\vartheta)(1-2\kappa\eta_0)
                          (\tfrac12+\sigma_0).
 \tag{Q19}
\]
Even the crude short-value sum has a power saving.

Write \(a=a_{\rm aux}a'\), with \(a_{\rm aux}\) the entire part supported
on \(vk\). Its part \(a_{\rm aux}>H\) has the same Rankin saving.
The exact BGLM parameters are
\[
 a_1=M,\quad a_2=ujc_0,\quad a_3=1,\quad
 d=vka_{\rm aux},\quad c=L,\quad\text{varying modulus }a'.
\]
All source coprimalities hold after zero-solution terms are removed:
\((a',a_1a_2d)=1\), the auxiliary residue is a unit, and the new short
variable is coprime to \(a_2\).
All varying and auxiliary moduli and inverse multipliers are prime to
\(L\), but may be divisible by 3. The residue \(a_1=M\) is excluded
from this assertion; in particular it is even.
The fixed value of \(\chi(a')\) in a reduced class modulo \(L\)
supplies the required signed modulus coefficient.

The source size conditions are
\[
\begin{gathered}
 11\lambda+\vartheta<(1-13\lambda-\vartheta)\delta_0,\\
 1<(1-13\lambda-\vartheta)(1+\delta_0),\\
 \tfrac12<(1-13\lambda-\vartheta)(\tfrac12+\delta_0),\\
 r<(1-13\lambda-\vartheta)\delta_0,\\
 \eta_0(1-3\lambda-\vartheta)>10\lambda(1-\eta_0),\\
 L\le N^{(1-13\lambda-\vartheta)\delta_0}.
\end{gathered}
 \tag{Q20}
\]
All hold for sufficiently small \(\lambda,\vartheta\), with the last
eventually true because \(L\) is fixed.

There are \(\kappa+1\) fine-box coordinates.
The box count is \(N^{(\kappa+1)\zeta+\epsilon}=N^{\lambda/100+\epsilon}\).
The two kernel masses cost \(N^{6\lambda/5+\epsilon}\), and
\(c_0,a_{\rm aux}\) each have at most \(N^{10\lambda}\) choices.
Thus even \(100\lambda+100\vartheta\) safely bounds the total polynomial
loss, after absorbing the fixed-\(\kappa\) logarithms in \(\epsilon\).
The dispersion and nonprincipal models are bounded respectively by
\[
 N^{1-r/3+100\lambda+100\vartheta+\epsilon},
 \qquad
 N^{\beta_\kappa+6r+100\lambda+100\vartheta+\epsilon}.
 \tag{Q21}
\]
For the model, every nonprincipal conductor is prime to \(L\), so cannot
divide \(L\); (Q15) applies. Choose \(\lambda<r/10^6\) and then
\(\vartheta\) small enough. Both exponents are below 1.

The moving boundary \(ux+vka^2=M\) is handled by the same actual integer
count after writing \(M-ux=vkat\):
\[
 |vka(t-a)|\ll\Delta N.
\]
Recombining first kernels at \(j\mid x\), with
\(\sum|e(j)|\tau_\kappa(x/j)\ll N^\epsilon\), gives error
\[
 N^\epsilon(\Delta N\log N+\sqrt N K^{1/10}).
 \tag{Q22}
\]
The principal model retains \(j^{-1}\); nonprincipal models retain
\(j^{-\beta_\kappa}\), allowed by (Q14). No \(R^2\) multiple of a
boundary length or scalar tail is used.
Square endpoints and the completed divisor-character series have the
same square-root-scale errors as in AsymUniform.

Choose a core saving smaller than \(\zeta/1000\), then a final family
power \(\theta\) smaller than both this saving and \(\vartheta/100\).
Coefficient reduction replaces \(N\) by at least \(N^{1-3\theta}\),
and new parameters by at most \((N')^{3\theta}\).
Its divisor prefactors preserve a smaller uniform saving.
This verifies the analytic estimates behind (Q3)--(Q4).

## 8. Reassembly of the main term

In a fixed gcd stratum, the first principal density is its deleted mean
times \(E_A(vka)/\varphi(vka)\), with real length \((M-vka^2)/u\).
The signed divisor series uses the fixed \(L(s,\chi)\), with finite
Euler deletions and an absolutely convergent correction.
Its partial sums have only \(N^\epsilon\) size, so completing it and
the second kernel has the same saved errors as before.

The local double sums are \(I_{p,0}/M_A\) outside \(uvM\), 1 at \(p\mid u\),
\(H_B/H_A\) at \(p\mid v\), and first-mean deletion at \(p\mid M\).
Therefore the fixed-g coefficient is \(2\kappa/L\) times the product of
\[
 p\mid M:1-p^{-1};\quad
 p\mid u:M_A;\quad p\mid v:M_B;\quad
 p\nmid uvM:I_{p,0}.
\]
Its real prefactor is \(M F_B(g)G_C(g)/(uv)\).
Completing \(g\) gives (Q9) at every target valuation: smaller valuations
give \((1-p^{-1})A_iB_ip^{-i}\), and the largest one gives the shifted
unit-target integral. This is exactly the line integral (Q7).
Finally, (Q16) implies
\[
 \frac{2L}{\kappa}\rho_{F_B}\rho_{G_C}
       =\frac{2\kappa}{L}\prod_{p\nmid L}M_A M_B.
\]
This proves (Q8), hence (Q2)--(Q4), with every fixed prime and coefficient
prefactor retained.

## 9. Unweighting: every local and logarithmic exponent

Every prime in this section satisfies \(p\nmid L\), so \(p\ge3\).
The exact fully rough mass, where both components are units, is
\[
 c_e(p)=
 \begin{cases}1-2/p,&e=0,\\1-1/p,&e\ge1.\end{cases}
\]
For the original \(U_q,V_\chi\), the rough bad density
\(1-c_e/I_{p,e}\) is
\[
 \begin{array}{c|cc}
 &e=0&e\ge1\\\hline
 p\equiv1\pmod L&(\kappa+2)/p+O_q(p^{-2})&
                         2\kappa/p+O_q(p^{-2})\\
 \chi(p)=1,\ p\not\equiv1\pmod L&2/p+O(p^{-2})&0\\
 \chi(p)=-1&0&0 .
 \end{array}
 \tag{Q23}
\]
All bounds are uniform in \(e\), with \(0\le g_e(p)<1\).
At finitely many small primes, including 3 when present, treat events
exactly. Their unit mass is at least \(1/3\), so the model stays positive.

The class densities are \(1/\kappa\), \(1/2-1/\kappa\), \(1/2\).
Hence the total rough dimension is
\[
 (\kappa+2)/\kappa+2(1/2-1/\kappa)=2.
\]
For
\[
 \mathcal V(z)=\prod_{\substack{p\le z\\p\nmid L}}
                 \frac{1-2/p}{M_A M_B},
\]
fixed-progression Mertens gives \(\mathcal V(z)\asymp(\log z)^{-2}\).
The tail beyond \(z=N^\gamma\), including the density
\(2L/\kappa\) and positive residues, is bounded above and below uniformly
in \(N\): unit factors are \(1+O_q(p^{-2})\), and at most \(1/\gamma\)
target primes exceed \(z\), each changing its factor by \(1+O_q(1/z)\).

Use pairwise coprime squarefree cofactors \(d,u,v\), with
\(d\mid N\), \(d,u\in\mathcal C_q\), \(v\in\mathcal A_\chi\),
\((uv,N)=1\), primes at most \(w<z\), and \(duv\le N^\alpha\).
At target \(M=N/d\), prescribe valuation one on first cofactor primes
and on second cofactor primes, then sieve all remaining primes through
\(z\). Normalize these directive weights by
\(\kappa^{\omega(u)}2^{\omega(v)}\).
The exact main is
\[
 C_N(z)\frac Nd\,\mathcal V(z)
  \prod_{\substack{p\mid N/d\\p\le z,\ p\nmid L}}E_p
  \prod_{p\mid uv}a_p,
 \qquad
 E_p=\frac{p-1}{p-2},\quad
 a_p=\frac{p-1}{p(p-2)}.
 \tag{Q24}
\]
Indeed valuation one has mass \((1-p^{-1})U_q(p)/p\) or
\((1-p^{-1})V_\chi(p)/p\); after division by \(\kappa\) or 2,
its ratio to the unit-target rough mass is exactly \(a_p\).

A finite-block lower sieve proves a fixed fraction of (Q24).
Split primes into \((z^{2^{-j-1}},z^{2^{-j}}]\).
Their local bad-density sums are bounded by a fixed constant.
Use even Bonferroni \(U_j\), odd Bonferroni \(L_j\), and the lower polynomial
\[
 \prod_jU_j-\sum_j(U_j-L_j)\prod_{i\ne j}U_i .
\]
Its model is at least half the independent unit probability by factorial
bounds on the block errors. Degrees \(2(s+j),2(s+j)+1\) give divisor
support \(C_0z^{4s+6}\), and coefficients \(O(\log\log N)\).
Choose \(\gamma\) so this lies below \(N^\xi\).
All valuation and roughness predicates expand into (Q4) with parameters
at most \(N^{2\alpha+\xi}\). It suffices that
\[
 2\alpha+\xi<(1-\alpha)\theta,\qquad
 \alpha+\xi<\eta/4 .
\]
The total analytic error over the cofactor and sieve predicates is
\(N^{1-\eta+\alpha+\xi+o(1)}\), hence negligible.

At a unit target, restoration options give \(1+2a_p\) in
\(\mathcal C_q\) and \(1+a_p\) in
\(\mathcal A_\chi\setminus\mathcal C_q\).
The restored exponent is
\[
 2/\kappa+(1/2-1/\kappa)=1/2+1/\kappa.
\]
At every positive target valuation in \(\mathcal C_q\), the common
divisor option gives
\[
 E_p+p^{-1}E_{p,e-1}\ge E_p+1/p=1+2a_p.
\]
At an allowed target prime outside \(\mathcal C_q\),
\[
 E_p-(1+a_p)=1/[p(p-2)]>0.
\]
At inert target primes the uncompensated enhancement is exactly \(E_p\);
their product is (Q6). These identities are valid at \(p=3\).

Thus the lower model has order
\[
 N J_\chi(N)(\log N)^{-2}(\log w)^{1/2+1/\kappa}
 \asymp \frac{N J_\chi(N)}{(\log N)^{3/2-1/\kappa}}.
\]
Under the actual positive cofactor model each prime is included with
probability \(O_q(1/p)\). Taking \(w\) to a sufficiently smaller fixed
power and applying Markov to \(\log(duv)\) retains a fixed fraction
inside \(duv\le N^\alpha\).

Lift residual values by \(d\), without inserting a weight
\(U_q(d)V_\chi(d)\). The valuation patterns
\((1,1),(1,0),(0,1)\) recover \(d,u,v\) uniquely.
The remaining rough weights are bounded, because
\[
 \binom{e+\kappa-1}{\kappa-1}\le\kappa^e,\qquad e+1\le2^e,
\]
and rough integers have at most \(1/\gamma\) factors with multiplicity.
The resulting bound is \((2\kappa)^{1/\gamma}\).
This proves the lower half of (Q5).

For the upper half, inert primes exclude both residues \(0,N\);
allowed primes outside \(\mathcal C_q\) exclude only the first residue0.
The dimension is \(2(1/2)+(1/2-1/\kappa)=3/2-1/\kappa\).
Only inert target primes merge two residues, giving exactly \(E_p\).
The nonnegative even-block upper sieve, at divisor level at most
\(N^{1/2}\), has CRT endpoint error \(O(N^{1/2}\log^C N)\).
It yields the upper bound in (Q5).

## 10. Bounded excluded prime factors at the full two-sided scale

Define
\[
 J_q(N)=\prod_{\substack{p\mid N\\p\nmid L\\p\not\equiv1\ (q)}}
                         \frac{p-1}{p-2}.
 \tag{B1}
\]
This includes \(p=3\) whenever it is outside \(L\) and divides the target.
The prime 2 is excluded because \(L\) is even.

**Theorem BoundedExcludedPrimeBridge.** For the fixed \(q,\chi\), there
exists some sufficiently small fixed \(\gamma>0\) such that
\[
\begin{split}
 R_\gamma(N)=\#\{a+b=N:\;&a\in\mathcal C_q,\ b\in\mathcal A_\chi,\\
             &p\mid b,\ p\not\equiv1\pmod q\Rightarrow p>N^\gamma\}
\end{split}
\]
satisfies
\[
 \boxed{\displaystyle
 R_\gamma(N)\asymp_{q,\chi}
       \frac{N J_q(N)}{(\log N)^{2-2/\kappa}}}
 \qquad(N\equiv2\pmod L,\ N\ \text{sufficiently large}).
 \tag{B2}
\]
In every counted pair,
\[
 \Omega_{\rm excluded}(b):=
       \sum_{\substack{p\mid b\\p\not\equiv1\ (q)}}v_p(b)<1/\gamma.
 \tag{B3}
\]
Thus the assertion counts multiplicities in its bounded number of
excluded prime factors. It does not assert that these factors disappear.

For the lower bound, use the same rough construction as in Section9, but
allow only cofactors
\[
 d,u,v\in\mathcal C_q,\quad d\mid N,\quad (uv,N)=1,
 \quad duv\le N^\alpha.
\]
In particular omit every second cofactor prime outside \(\mathcal C_q\).
The first normalization remains \(\kappa^{\omega(u)}\), the second
\(2^{\omega(v)}\). The rough sieve and its distribution errors are
unchanged: it still makes every residual prime at most \(z=N^\gamma\)
absent, except the exact valuation-one directives.

At a unit target prime in \(\mathcal C_q\), the available options are
\(1,a_p,a_p\). At a target prime in that class the common divisor option
still gives
\[
 E_p+p^{-1}E_{p,e-1}\ge1+2a_p.
\]
There are no cofactor options at any other prime. All nonidentity-class
target enhancements \(E_p\) therefore remain in the main, whether
\(\chi(p)=1\) or \(\chi(p)=-1\).
Their product below \(z\) is \(J_q(N)\) up to a uniform \(1+o(1)\):
at most \(1/\gamma\) target primes exceed \(z\), and each missing
factor is \(1+O(1/z)\).

The restoring product is now only
\[
 \prod_{\substack{p\le w\\p\equiv1\ (q)}}(1+2a_p)
       \asymp_q(\log w)^{2/\kappa}.
\]
Together with the rough dimension two this gives lower model mass
\[
 \gg_{q,\chi}N J_q(N)(\log N)^{-2}(\log w)^{2/\kappa}
 \asymp_{q,\chi}
       \frac{N J_q(N)}{(\log N)^{2-2/\kappa}}.
\]
The actual positive model still has inclusion probabilities \(O_q(1/p)\).
Choose the fixed power \(w=N^\omega\) small after \(\gamma,\alpha\);
Markov retains a fixed fraction on \(duv\le N^\alpha\).
The same inequalities
\(2\alpha+\xi<(1-\alpha)\theta\),
\(\alpha+\xi<\eta/4\)
make the full summed error negligible.

The lift by \(d\) and the recoverable patterns
\((1,1),(1,0),(0,1)\) are unchanged.
All \(d,u,v\) primes are in \(\mathcal C_q\), so any excluded prime
in the resulting \(b\) lies in its rough residual and exceeds
\(N^\gamma\). Normalized fibres remain at most
\((2\kappa)^{1/\gamma}\). This proves the lower bound in (B2).
Taking the product of the excluded prime powers and using \(b<N\)
proves (B3).

For the upper bound, every counted pair has both components coprime
to each nonidentity-class prime at most \(N^\gamma\).
At such a prime, sieve the first variable out of residues \(0,N\);
the two coincide exactly when \(p\mid N\).
The dimension is \(2(1-1/\kappa)=2-2/\kappa\), and target coincidences
supply exactly \(E_p\), hence the product (B1).
Use the nonnegative block upper sieve at a sufficiently small fixed
power cutoff, no larger than \(N^\gamma\), with divisor level at most
\(N^{1/2}\). CRT endpoint errors are
\(O(N^{1/2}\log^C N)\), giving the matching upper bound.

Only the existence of some sufficiently small fixed \(\gamma\) is
proved. There is no claim for an arbitrary prescribed \(\gamma\),
or for \(\gamma\ge1/2\), and no fixed-ratio restriction on the pair.
For \(\kappa>2\), this is an every-target almost-splitting result at
the full two-sided logarithmic count scale, rather than a proof of
\(\mathcal C_q+\mathcal C_q\).


## 11. Specializations and scope

For \(q=5,L=10\), take the even quadratic \(\chi_5\).
The allowed unit classes are 1 and9 modulo10, the fixed factor is5,
and 3 is inert. A target divisible by3 has \(J_\chi\)-factor
\((3-1)/(3-2)=2\).

For \(q=8,L=8\), take \(\chi_{-4}\) inflated to8.
Its allowed unit classes are1 and5, its fixed factor is4, and3 is inert.
The period-8 theorem must therefore allow moduli divisible by powers of3,
as Section5 proves. As an additional split-prime3 control,
\(\chi_{-8}\bmod8\) allows classes1 and3; then a target prime3 has the
one-sided allowed-prime compensation \(E_3-(1+a_3)=1/3\), rather than
an inert \(J_\chi\)-factor.

At \(\kappa=2\), the kernel of \(\chi\) consists only of1 modulo \(L\),
so \(\mathcal A_\chi=\mathcal C_q\). The theorem includes those exact two-sided unrestricted
baseline cases (\(q=3,4,6\)). For \(\kappa>2\), no general
\(\mathcal C_q+\mathcal C_q\) conclusion follows; in particular the
critical two-sided modulus-12 remainder is still open.

A nonprincipal quadratic character exists for every \(q\ge3\):
if an odd prime divides \(q\), inflate its Legendre character to \(L\);
if \(q\) is a power of2, use \(\chi_{-4}\).
The theorems concern each chosen fixed character and modulus.
There is no uniformity for growing \(q\), effective cutoff, identified
unweighted leading constant, external review, or novelty assertion.


## 12. Relation to shifted divisor correlations

Drappeau and Topacogullari,
*Combinatorial identities and Titchmarsh's divisor problem for
multiplicative functions*, ANT13(2019),2383--2425,
[DOI](https://doi.org/10.2140/ant.2019.13.2383),
[primary published article](https://msp.org/ant/2019/13-10/ant-v13-n10-p05-s.pdf),
Theorem1.2 on printedp2386, treat prime-periodic multiplicative \(f\)
against \(\tau(an-h)\), for positive \(a\) and power-small \(|h|\).
Their unconditional final error is logarithmic to arbitrary fixed order.
The same page explicitly distinguishes the dual additive sum
\(\sum_{n<N}f(n)\tau(N-n)\) as outside their method.
Thus this theorem does not directly supply (Q3), even before replacing
ordinary \(\tau\) by \(\chi_0*\chi\) and its square correction.

Topacogullari, *The shifted convolution of generalized divisor functions*,
IMRN2018,7681--7724,
[DOI](https://doi.org/10.1093/imrn/rnx111),
[primary preprint](https://arxiv.org/pdf/1605.02364v1),
Theorems1.1--1.2 on preprintp2, give power savings for each fixed
number of divisor factors \(k\ge4\), with shift size \(x^{15/19}\).
The discussion and equation(2.2) permit weaker sublinear-shift ranges.
These are shifted products with ordinary divisor function, not the
all-target additive quadratic-character family stated here.
Arbitrary fixed numbers of divisor factors therefore have substantial
precedent; no priority claim follows from these scope comparisons.

## 13. Analytic dependencies

The exact BGLM Proposition5.1 statement, including the signed modulus
sum and coprimality requirements, is retained in
[AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md).
The original periodic-box derivation, all-modulus complete correlation,
three source estimates and primary Heath-Brown scan are retained in
[FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md).
Section5 supplies the full fixed-period transfer needed here, rather than
assuming a period-12 theorem for other periods.
The exact coefficient and local-density method is also developed in
[AsymUniform](ASYMMETRIC_UNIFORM_COEFFICIENTS.md);
the finite block sieve and cofactor mechanism are recorded in
[AsymSemigroupScale](ASYMMETRIC_SEMIGROUP_SCALE.md).
The new quadratic signs, prime3 cases, fixed-place density and logarithmic
exponents have been derived above.

All statements concern each fixed \(q,\chi\). They give no uniformity
for moduli growing with the target, effective cutoff, identified
unweighted leading constant, or fixed-ratio interval theorem.
The critical and beyond-half weighted/signed remainder obligations of
this sieve route remain unresolved. The full unrestricted
\(\mathcal C_{12}+\mathcal C_{12}\) value-count asymptotic is now
proved by the separate [quantitative subgroup theorem](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md),
under its stated analytic inputs. That result does not settle these
weighted remainder questions or power-short windows.

