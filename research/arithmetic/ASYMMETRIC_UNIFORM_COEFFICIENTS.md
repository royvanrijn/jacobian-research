# Uniform small coefficients in the asymmetric convolution

The asymmetric pointwise convolution admits a uniform extension to
sufficiently small coefficients, shifted weights and divisibility tests.
The external family exponent is chosen after the analytic error saving.
This result uses the same analytic inputs as
[Theorem AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md):
BGLM Proposition 5.1 and
[Theorem FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md).
It gives no distribution estimate beyond the square-root range and
does not prove the two-sided modulus-12 theorem.

## Theorem AsymUniform

Let
\[
 \mathcal C=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod{12}\},\qquad
 \mathcal A=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod3\},
\]
including 1 in both sets. Put
\[
 U(n)=\tau(n)^2\mathbf1_{\mathcal C}(n),\quad
 V(n)=\tau(n)\mathbf1_{\mathcal A}(n),\quad
 F_B(m)=\frac{U(Bm)}{U(B)},\quad
 G_C(l)=\frac{V(Cl)}{V(C)}.
 \tag{1}
\]
There are absolute constants \(\eta,\theta>0\) such that, uniformly for
every sufficiently large \(N\equiv2\pmod6\) and
\[
 B,u\in\mathcal C,\qquad C,v\in\mathcal A,\qquad
 \max(B,C,u,v)\le N^\theta,
\]
one has
\[
 \boxed{\displaystyle
 \sum_{\substack{um+vl=N\\m,l\ge1}}F_B(m)G_C(l)
 =\frac{N}{uv}\mathfrak S_{B,C;u,v}(N)+O(N^{1-\eta}).}
 \tag{2}
\]
The constant is uniform in all parameters and target valuations.
No coprimality assumption is imposed. The exponent \(\theta\) may be
decreased arbitrarily after fixing \(\eta\).

Here is the local definition of the main term. For \(p\ge5\), put
\[
 A_e=F_B(p^e),\quad B_e=G_C(p^e),\qquad
 M_{A,p}=(1-p^{-1})\sum_{e\ge0}A_ep^{-e},\quad
 M_{B,p}=(1-p^{-1})\sum_{e\ge0}B_ep^{-e}.
\]
For residues modulo \(p^h\), set
\(\nu_h(x)=\min(v_p(x),h)\), with \(\nu_h(0)=h\), and define
\[
 D_p=\lim_{h\to\infty}p^{-h}
       \sum_{\substack{x,y\bmod p^h\\ux+vy\equiv N\ (p^h)}}
             A_{\nu_h(x)}B_{\nu_h(y)}.
 \tag{3}
\]
Then
\[
 \mathfrak S_{B,C;u,v}(N)
 =6\rho_{F_B}\rho_{G_C}
       \prod_{p\ge5}\frac{D_p}{M_{A,p}M_{B,p}},
 \qquad
 \rho_F=\operatorname*{Res}_{s=1}\sum_{n\ge1}F(n)n^{-s}.
 \tag{4}
\]
The residue means are positive. The local limit exists: after removing
the fixed coefficient valuations, the remaining zero-residue tails are
bounded by a polynomial in \(h\) times \(p^{-h}\).
Outside the primes dividing \(BCuvN\), the product factor is
\(1+O(p^{-2})\), with an absolute constant. Thus (4) converges absolutely.
It vanishes exactly in the coefficient/target local obstruction cases
described below; all remaining factors are positive.

The factor 6 is the fixed small-modulus density. Modulo 12, the normalized
residue measures of \(F_B\) and \(G_C\) are respectively point mass at 1
and equal masses \(1/2\) at 1 and 7. Therefore
\[
 12\sum_{\substack{a,b\bmod12\\ua+vb\equiv N\ (12)}}
              \mu_F(a)\mu_G(b)=6
 \tag{5}
\]
for \(u\equiv1\pmod{12}\), \(v\equiv1\pmod6\), \(N\equiv2\pmod6\).
Shifts of \(C\) at primes \(7\bmod12\) do not alter this factor.
The equal second-side means follow also from the literal mod-6
two-character factorization: its nonprincipal mod-4 twists have no
pole, and its finite shifted kernel preserves that equality.

### Divisibility events

If \(b\notin\mathcal C\) or \(c\notin\mathcal A\), the count restricted
by \(b\mid m,c\mid l\) is zero. Otherwise
\[
\begin{split}
 \sum_{\substack{um+vl=N\\b\mid m,\ c\mid l}}F_B(m)G_C(l)
 ={}&F_B(b)G_C(c)\frac{N}{ubvc}
          \mathfrak S_{Bb,Cc;ub,vc}(N)
       +O(N^{1-\eta}).
\end{split}
 \tag{6}
\]
After decreasing \(\theta\), this is uniform when each of
\(B,C,u,v,b,c\) is at most \(N^{\theta/10}\).
A sufficiently small fixed bound on their total product is an equivalent
usable restriction. The prefactor in (6) is exact.
The divisors may have powers. Exact prescribed valuations follow by
finite differences between divisibility by \(p^e\) and \(p^{e+1}\),
provided each tested product remains in range.
Summing at most \(N^{\varepsilon_{\rm family}}\) events adds that
exponent to the error; \(\eta\) can accordingly be reduced.
This is a small positive-power interface, not a large divisor level.

## 1. Local factors at arbitrary valuations

For each \(p\), write \(a=v_p(u)\), \(b=v_p(v)\), \(n=v_p(N)\),
\(t=\min(a,b)\). If \(n<t\), then \(D_p=0\). Otherwise subtract \(t\)
from \(a,b,n\), and multiply the answer below by \(p^t\).
At least one of the remaining coefficient valuations is zero.

For \(x=p^{-1}\), \(q=1-x\), define
\[
\begin{split}
 I_{p,n}(A,B)={}&q\sum_{i<n}A_iB_ix^i
        +x^n(1-2x)A_nB_n\\
 &+q\left(B_n\sum_{i>n}A_ix^i+
                A_n\sum_{j>n}B_jx^j\right).
\end{split}
 \tag{7}
\]
If \(a=b=0\), then \(D_p=I_{p,n}(A,B)\).
If \(a>0,b=0\), then
\[
 D_p=
 \begin{cases}
 B_nM_{A,p},&n<a,\\
 B_a I_{p,n-a}(A,B^{[a]}),&n\ge a,
 \end{cases}
 \qquad B^{[a]}_j=B_{a+j}/B_a.
 \tag{8}
\]
If \(b>0,a=0\), then
\[
 D_p=
 \begin{cases}
 A_nM_{B,p},&n<b,\\
 A_b I_{p,n-b}(A^{[b]},B),&n\ge b,
 \end{cases}
 \qquad A^{[b]}_j=A_{b+j}/A_b.
 \tag{9}
\]
A zero normalizing coefficient in these expressions gives \(D_p=0\).
In these families all subsequent coefficients at that prime are then
also zero.

To prove these identities, remove the common coefficient \(p^t\)
from the congruence. The \(p^{2t}\) lifts, with the normalization
\(p^{-h}\), give a net factor \(p^t\).
For \(a>0,b=0\), write the second variable as
\(v^{-1}(N-ux)\). Its valuation is exactly \(n\) when \(n<a\).
When \(n\ge a\), remove \(p^a\) from that variable, producing its shifted
weight in (8). The other case is symmetric. Unit changes of variables
remove coefficient unit parts, so only valuations occur.
The valuation partition for two unit coefficients gives (7).

In the coefficient-target-unit case
\[
 (u,v)=(uv,N)=1,
\]
these formulas reduce to
\[
 p\mid u:\ D_p=M_{A,p},\qquad
 p\mid v:\ D_p=M_{B,p},\qquad
 p\nmid uv:\ D_p=I_{p,v_p(N)}(A,B).
 \tag{10}
\]
For \(r=(u,v)\), the common-coefficient factors multiply to \(r\)
when \(r\mid N\). Thus the scale using the original coefficients is
\(rN/(uv)\).

## 2. Exact coefficient reduction

First put \(r=(u,v)\). If \(r\nmid N\), there are no solutions.
Otherwise divide \(u,v,N\) by \(r\), keeping variables and shifts
unchanged. Here \(r\in\mathcal C\), so \(r\equiv1\pmod6\).

Now assume \((u,v)=1\), and put \(d=(u,N)\), \(e=(v,N)\).
If \(e\) contains a prime \(7\bmod12\), every solution forces that
prime to divide \(m\), so the count is zero.
Otherwise \(d,e\in\mathcal C\), \((d,e)=1\), and
\[
 m=em',\quad l=dl',\qquad
 u'=u/d,\quad v'=v/e,\quad N'=N/(de)
 \tag{11}
\]
gives exactly
\[
 u'm'+v'l'=N',\qquad
 F_B(e)G_C(d)F_{Be}(m')G_{Cd}(l').
\]
Primewise minimum valuations give
\((u'v',N')=(u',v')=1\).
The target remains \(2\bmod6\).
The prefactors satisfy
\[
 F_B(e)\le\tau(e)^2,\qquad G_C(d)\le\tau(d)
 \tag{12}
\]
uniformly in the shift valuations.

Starting with parameters at most \(N^\theta\), the complete reduction
has \(N'\ge N^{1-3\theta}\), and transformed parameters at most
\((N')^{3\theta}\). A theorem in the unit case with family power
\(\vartheta\) therefore proves (2) with \(\theta=\vartheta/100\)
and a slightly smaller error saving.

The local model is preserved exactly. Removing a common coefficient
gives \(p^{v_p(r)}\). The substitutions in (11) transform \(D_p\)
into
\[
 F_B(p^{v_p(e)})G_C(p^{v_p(d)})D'_p.
\]
The variable Jacobian cancels the density factor from dividing the
equation. After coprime coefficients have been reached,
\(N'/(u'v')=N/(uv)\). The residue mean changes occur only in the
corresponding finite local \(M\)-factors. Thus (4) retains every
coefficient and gcd prefactor.

## 3. Uniform kernels and character means

At a prime \(p\in\mathcal C\), the first shifted local coefficients are
\[
 A_e=\left(\frac{v_p(B)+e+1}{v_p(B)+1}\right)^2.
\]
They vanish for \(e\ge1\) at other primes.
At \(p\in\mathcal A\), the second coefficients are
\[
 B_e=\frac{v_p(C)+e+1}{v_p(C)+1},
\]
and they vanish for \(e\ge1\) outside \(\mathcal A\).

The first Dirichlet series retains the four literal mod-12 character
factors, with shifted correction of degree at most three and uniformly
bounded coefficients. The second retains the two literal mod-6 factors.
At every shifted prime \(p\in\mathcal A\), including \(p\equiv7\pmod{12}\),
its replacement correction is
\[
 (1-z)^2\sum_{e\ge0}\left(1+\frac e{v_p(C)+1}\right)z^e
       =1-\frac{v_p(C)}{v_p(C)+1}z.
\]
Deleting first-side Euler factors at \(T\) replaces the whole local
correction by the polynomial canceling the four \(L\)-factors.
Every kernel \(e\) therefore satisfies
\[
 |e(j)|\le\tau(j)^{O(1)},\qquad
 \sum_{j\ge1}|e(j)|j^{-3/5}
        \ll_\varepsilon(BCT)^\varepsilon.
 \tag{13}
\]
The constants and exponent are independent of all valuations.
First kernels remain \(1\bmod12\); undeleted second kernels remain
\(1\bmod6\). This latter support is sufficient for divisor pairing.

The four-character means and twists from the pointwise proof hold
uniformly with \((BCT)^\varepsilon\) losses. With
\(H_{B,p}(z)=\sum_{i\ge0}F_B(p^i)z^i\), the principal density after
first-side deletion at \(T\), in a reduced class modulo \(L\), is
\[
 \frac{\rho_{F_B}}{\varphi(L)}
 \prod_{p\mid T}H_{B,p}(1/p)^{-1}
 \prod_{p\mid L}H_{B,p}(1/p)^{-1},
 \qquad (L,T)=1,\ (L,6)=1.
 \tag{14}
\]
For primitive conductor \(h\nmid12\), the cumulative nonprincipal
twist error is
\[
 \ll_\varepsilon h(BCTY)^\varepsilon Y^{3/4}.
 \tag{15}
\]
Kernel summation retains \(j^{-3/4}\).
The smooth scalar decomposition also has finite reciprocal sum at
\(3/4\), up to \((uj)^\varepsilon\). No conductor count is multiplied
by an actual boundary length or Rankin tail.

## 4. Unit-case gcd strata and hyperbola

Assume \((u,v)=(uv,N)=1\), and put \(g=(m,N)\).
Then \(g\in\mathcal C\), \(g\mid N\), \((g,uv)=1\).
Writing \(m=gx,l=gy,M=N/g\) gives
\[
 ux+vy=M,\qquad (x,M)=1,
\]
and the exact weight
\[
 F_B(g)G_C(g)F_{Bg}(x)G_{Cg}(y).
 \tag{16}
\]
For \(K=N^\lambda\), strata \(g>K\) cost
\(O(N^{1-\lambda+\varepsilon})\): all weights and prefactors are
divisor-bounded, and the number of pairs in a stratum is at most
\(N/g+1\).

For \(g\le K\), truncate the second kernel at \(k\le K\), with error
\(O(N^{1-2\lambda/5+\varepsilon})\).
The original conditions force
\[
 (k,uM)=1,\qquad(q,6uM)=1.
\]
Every nonzero residual
\(t=(M-ux)/(vk)\) is \(1\bmod6\).
Hence \(\chi_3(t)=1\), and strict divisor pairing is exact.
Its square endpoint costs
\(O(N^{1/2+\lambda/10+\varepsilon})\).
The remaining normalized fixed-\(g\) count is
\[
\begin{split}
 2\sum_{\substack{k\le K\\(k,uM)=1}}e_{Cg}(k)
 \sum_{\substack{q\le\sqrt{M/(vk)}\\(q,6uM)=1}}\chi_3(q)
 \sum_{\substack{x\le(M-vkq^2)/u\\ux\equiv M\ (vkq)}}
               F_{Bg}(x)\mathbf1_{(x,M)=1}.
\end{split}
 \tag{17}
\]
The residue and first-side deletion modulus are both \(M\).

## 5. Source conditions with small coefficients

Choose \(\eta_0,\delta_0,r\) as in the pointwise proof, then an internal
\(\lambda>0\). Put
\[
 K=N^\lambda,\quad H=N^{10\lambda},\quad R=N^r,\qquad
 \zeta=\lambda/100,\quad\Delta=N^{-\zeta}.
\]
Choose the external unit-family exponent \(\vartheta\) last.
A sufficient initial hierarchy is
\[
 \begin{gathered}
 r<\min(\eta_0,\delta_0,\sigma_0,1/100)/10000,\\
 \lambda<r/10^6,\qquad \vartheta<\lambda/10^{12}.
 \end{gathered}
 \tag{18}
\]
Decrease \(\lambda\) first and \(\vartheta\) last for the strict
inequalities below. In this section all four parameters are at most
\(N^\vartheta\).

Discard \(um\le N^{1-\lambda}\), at cost
\(O(N^{1-\lambda+\varepsilon})\).
After \(g,j\le K\), a first four-character product has scale
\[
 X\ge N^{1-3\lambda-\vartheta}.
\]
If a variable has scale in
\([X^{2\eta_0},X^{1/4+\eta_0}]\), write it uniquely as \(c_0n'\),
where \(c_0\) is its entire part supported on primes of \(uj\).
The tail satisfies
\[
 \sum_{\substack{c_0>H\\c_0\mid(uj)^\infty}}\frac1{c_0}
 \le H^{-1/2}\prod_{p\mid uj}(1-p^{-1/2})^{-1}
 \ll_\varepsilon N^{-5\lambda}(uj)^\varepsilon.
 \tag{19}
\]
Counting multiples gives an actual-count error
\(O(N^{1-5\lambda+\varepsilon})\), without a number-of-scalars
endpoint loss. Nonprincipal model tails use (15) and the smooth Euler
product at \(3/4\).
For \(c_0\le H\), the new short variable is coprime to the full inverse
multiplier \(ujc_0\), and
\[
 x_{\rm BGLM}\ge N^{1-13\lambda-\vartheta}.
\]

Write \(q=q_{\rm aux}q'\), with \(q_{\rm aux}\) the entire part supported
on \(vk\). Discard \(q_{\rm aux}>H\) by the same Rankin bound.
The actual progression count is at most
\(N^\varepsilon(N/(uvkq)+1)\).
The reciprocal term has the Rankin saving, while the endpoint is
\(O(N^{1/2+\lambda/10+\varepsilon})\) after (13).
Model tails again use their \(3/4\) twist and the smooth Euler product.

The exact BGLM parameters are
\[
 a_1=M,\quad a_2=ujc_0,\quad a_3=1,\qquad
 d=vkq_{\rm aux},\quad c=12,\quad\text{varying modulus }q'.
 \tag{20}
\]
Terms with \((uj,vkq)>1\) have no solution, since \(M\) is a unit
there; their principal unit model is also zero.
In every remaining term,
\((q',a_1a_2d)=1\), the auxiliary residue
\(M\overline{ujc_0}\bmod d\) is a unit, and \((n',a_2)=1\).
All moduli are prime to 6.
The signed \(\chi_3(q')\) coefficient is constant in a fixed reduced
class modulo 12.

The sizes are
\[
 |a_2|,d\le N^{11\lambda+\vartheta},\quad |a_1|\le N,\quad
 Q'\le N^{1/2},\quad R=N^r,\quad
 x_{\rm BGLM}\ge N^{1-13\lambda-\vartheta}.
\]
A sufficient list of source inequalities is
\[
\begin{gathered}
 11\lambda+\vartheta<(1-13\lambda-\vartheta)\delta_0,\\
 1<(1-13\lambda-\vartheta)(1+\delta_0),\\
 \tfrac12<(1-13\lambda-\vartheta)(\tfrac12+\delta_0),\\
 r<(1-13\lambda-\vartheta)\delta_0,\\
 \eta_0(1-3\lambda-\vartheta)\ge10\lambda(1-\eta_0).
\end{gathered}
 \tag{21}
\]
The upper relative short exponent decreases on division by \(c_0\);
the last inequality supplies its lower bound.
Thus \(\eta_0\) is unchanged after choosing \(\delta_0\).

If no medium variable exists, fix the short variables, whose product is
at most \(X^{8\eta_0}\). At most three long variables remain, with product
\(Y\ge X^{1-8\eta_0}\). Apply FixedPeriodicBoxes to the full modulus
\(vkq\). Indeed
\[
 vkq\le\sqrt{vkM}\le N^{1/2+(\lambda+\vartheta)/2},
\]
and it suffices that
\[
 \tfrac12+\lambda+\vartheta
 <(1-3\lambda-\vartheta)(1-8\eta_0)
       (\tfrac12+\sigma_0).
 \tag{22}
\]
This follows from \(\eta_0<\sigma_0/1000\) and (18).
A nonunit fixed scalar gives zero in both the progression and its
principal unit model. The total periodic-box error is, crudely,
\[
 O(N^{1-\sigma_0+8\eta_0+100\lambda+100\vartheta+\varepsilon}).
\]

A fixed \(d\) in (20) is selected by an indicator in the auxiliary
average of BGLM Proposition 5.1. No arbitrary \(q'\)-dependent
coefficient is introduced.
The two kernel masses cost at most \(N^{6\lambda/5+\varepsilon}\);
\(c_0,q_{\rm aux}\) each have at most \(H\) choices.
Fine boxes have at most \(N^{5\zeta+\varepsilon}\) choices.
A safe aggregate BGLM error is
\[
 O(N^{1-r/3+100\lambda+100\vartheta+\varepsilon}).
 \tag{23}
\]
Nonprincipal models, including the box count, are bounded by
\[
 O(N^{3/4+6r+100\lambda+100\vartheta+\varepsilon}).
 \tag{24}
\]
Every nonprincipal primitive conductor here fails to divide 12,
because the moduli are prime to 6. Thus (15) applies.
Both (23) and (24) save a fixed power under (18).

## 6. Moving boundaries with coefficients

Use relative-width \(\Delta\) boxes in the four character variables and
\(q\). Crossing boxes satisfy
\[
 |ux+vkq^2-M|\ll\Delta N.
\]
Restore divisibility first and write \(M-ux=vkqt\). Then
\[
 |vkq(t-q)|\ll\Delta N.
\]
Dropping only the integrality requirement on \(x\), the number of
integer pairs is
\[
 O\!\left(\frac{\Delta N\log N}{vk}+\sqrt{\frac N{vk}}\right).
\]
At each original \(x\), the first-kernel indices \(j\mid x\) and all
decomposition multiplicities have total weight \(N^\varepsilon\).
The scalar decompositions are unique; disjoint boxes add no actual
multiplicity. Summing the second kernel by (13) gives
\[
 O_\varepsilon\!\left(
 N^\varepsilon(\Delta N\log N+\sqrt N K^{1/10})\right).
 \tag{25}
\]
The principal-model interval at kernel \(j\) has length
\(O(\Delta N/(uj))\), and its reciprocal modulus retains \(1/(vk)\).
Its volume therefore uses the convergent sums
\(\sum|e(j)|/j\), \(\sum|e(k)|/k\).
Integer endpoints are far smaller than \(N^{1-\zeta}\).
Nonprincipal models use (15) box by box, retaining \(j^{-3/4}\).
This checks the moving boundary with both coefficients present.

## 7. Extraction of the local main term

For fixed \(g\), write \(A=F_{Bg}\), \(B=G_{Cg}\).
The principal part of (17) is
\[
 \frac2u\rho_{A,M}
 \sum_{\substack{k\\(k,uM)=1}}e_B(k)
 \sum_{\substack{q\le\sqrt{M/(vk)}\\(q,6uM)=1}}
 \chi_3(q)\frac{E_A(vkq)}{\varphi(vkq)}(M-vkq^2).
 \tag{26}
\]
Put \(W_p=[(1-p^{-1})H_A(1/p)]^{-1}\). Then
\[
 \frac{E_A(vkq)}{\varphi(vkq)}
 =\frac1{vkq}\prod_{p\mid vk}W_p
                \prod_{\substack{p\mid q\\p\nmid vk}}W_p.
\]
The \(q\)-series is the fixed \(L_6(s,\chi_3)\) times an absolutely
convergent correction for \(\Re s>0\), with finite deletions at \(uM\).
Its partial sums are
\(O_\varepsilon((Bg\,Cg\,uvMkQ)^\varepsilon)\).
Completing \(q\) and removing its \(q^2\) term costs
\[
 O\!\left(N^\varepsilon\frac{\sqrt M}{u\sqrt{vk}}\right).
\]
Summation over \(k\) costs \(K^{1/10}\).
Completing the \(k\)-sum costs
\(O(N^{1-2\lambda/5+\varepsilon})\).

At \(p\nmid uvM\), the local double sum equals
\[
 1+W_p\left[
   \frac{\sum_i e_B(p^i)p^{-i}}{1-\chi_3(p)/p}-1\right]
 =\frac{M_{A,p}+M_{B,p}-1}{M_{A,p}}.
\]
At \(p\mid u\), both \(k\) and \(q\) are excluded, so that sum is 1.
At \(p\mid v\), \(W_p\) is present even at zero local \(k,q\) valuations,
and the local sum equals
\[
 W_p\frac{\sum_i e_B(p^i)p^{-i}}{1-\chi_3(p)/p}
       =\frac{H_B(1/p)}{H_A(1/p)}.
\]
At \(p\mid M\), both indices are excluded and the first mean deletes
\(H_A\). Since
\(\rho_A=\frac13\prod_{p\ge5}M_{A,p}\) in natural prime order,
the fixed-\(g\) main is \(M/(uv)\,F_B(g)G_C(g)\) times
\(2/3\) times the product of
\[
 \begin{cases}
 1-p^{-1},&p\mid M,\\
 M_{A,p},&p\mid u,\\
 M_{B,p},&p\mid v,\\
 M_{A,p}+M_{B,p}-1,&p\nmid uvM.
 \end{cases}
 \tag{27}
\]

Complete the positive \(g\)-divisor sum. Its main tail is
\(O(N^{1-\lambda+\varepsilon})\), since shifted means and prefactors
are \(N^\varepsilon\).
At \(p^n\Vert N\), the choices \(v_p(g)=i<n\) contribute
\[
 (1-p^{-1})A_iB_ip^{-i}.
\]
The choice \(i=n\) contributes
\[
 p^{-n}A_nB_n
       (M_{A^{[n]},p}+M_{B^{[n]},p}-1).
\]
Expanding these shifted means gives exactly (7), with every higher
valuation retained. If \(A_n=0\), that \(g\)-choice is absent and the
remaining expression still equals (7).
At coefficient primes \(g\) has valuation zero, and (27) gives (10).
Thus the reassembled local factors are exactly (3).

Finally, in natural prime order,
\[
 \rho_{F_B}\rho_{G_C}
       =\frac19\prod_{p\ge5}M_{A,p}M_{B,p}.
\]
Hence the coefficient \(\frac23\prod D_p\) is exactly (4).
This proves (2) in the unit case. The exact reductions in Section 2
prove the stated general case.

## 8. Error saving and the external family

Kernels save \(2\lambda/5\); large gcd and small original-summand tails
save \(\lambda\); smooth scalar tails save \(5\lambda\).
Actual and principal boundaries save
\(\zeta=\lambda/100\).
The BGLM, periodic-box, nonprincipal model, square, \(q\)-completion
and gcd-completion errors save larger powers under (18), (21), and (22).

Fix
\[
 \eta_{\rm core}=\lambda/100000,\qquad
 \vartheta<\lambda/10^{12},\qquad
 \theta=\vartheta/100,\qquad \eta=\eta_{\rm core}/4,
 \tag{28}
\]
leaving a strict margin in each estimate, and take every divisor-loss
exponent \(\varepsilon<\eta_{\rm core}/10^6\).
The exact coefficient reduction has
\(N'\ge N^{1-3\theta}\) and transformed parameters
at most \((N')^{3\theta}\), inside the proved \(\vartheta\)-family.
Its divisor prefactors are absorbed by the reduction from
\(\eta_{\rm core}\) to \(\eta\).
The same reserve proves (6) and its finite-difference tests.
The external \(\theta\) is much smaller than the internal cutoff
\(\lambda\) and the final error saving. This makes the absolute error
useful at the family main scale and allows polynomially many small tests.

For locally admissible parameters, the main is positive and of size
\[
 \frac{rN}{uv}N^{o(1)},\qquad r=(u,v),
\]
uniformly. Normalized allowed-prime coefficients lie between 1 and
their unshifted divisor bounds. The means vary by products
\(1+O(p^{-1})\); target unit-unit strata give positive lower bounds;
and coefficient reductions introduce positive shift ratios.
The resulting finite products are bounded above and below by
\(N^{\pm\varepsilon}\) for every fixed \(\varepsilon>0\).
A fixed positive singular-series lower bound is not asserted.

## Scope and sources

The proof supplies only sufficiently small positive-power coefficient
and divisibility families. It does not supply moduli near
\(N^{1/2}\) as free auxiliary coefficients, a beyond-half sieve
remainder, the two-sided \(C+C\) theorem, or an unweighted
\(N/(\log N)^{5/4}\) estimate by itself.

The precise source statements and their primary PDFs are recorded in
[AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md) and
[FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md).
BGLM Proposition 5.1 is on printed pp.1192--1193 of
*Additive problems with almost prime squares*, GAFA **33** (2023),
1173--1242, [DOI](https://doi.org/10.1007/s00039-023-00635-w);
it is PDF p.15 of the
[cached arXiv v5](sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf).
Its absolute value is outside the varying-modulus sum, as used above.
The periodic-box result supplies \(\sigma_0=1/10000\) for moduli
prime to 12.

