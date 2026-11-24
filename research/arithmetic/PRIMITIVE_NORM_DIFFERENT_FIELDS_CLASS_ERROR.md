# Primitive norm class errors across different quadratic fields

Fix two primitive positive definite integral binary forms \(Q_1,Q_2\).
Their discriminants \(\Delta_i=f_i^2D_i<0\) may belong to different
imaginary quadratic fields. The class error is smaller than the natural
binary scale \(NJ_{12}(N)/\log N\), uniformly in every target valuation.
The only new factor relative to the same-field proof comes from a prime
dividing the target that is split in the marked field and inert in the
other field. This factor is explicit and costs at most a quarter power
of \(\log\log N\). No independence of ring-class primes across fields
is assumed, and no joint ring-class prime theorem is used.

## 1. Statement and fixed parameters

Put
\[
 G_i=\operatorname{Pic}(\mathcal O_i),\quad B_i=2G_i,\quad
 \chi_i=\chi_{D_i},\quad S=\{p:p\mid2\Delta_1\Delta_2\},\quad
 L_N=\log(3N),\quad \ell_N=1+\log L_N.
\]
Let \(P_i\) be the positive primitive numerical values of \(Q_i\),
and \(\Lambda_i\) their primitive local shadows. Counts
\(R_{A,B}(N)=\#\{a+b=N:a\in A,b\in B\}\) are ordered counts of
values, not representations or ideals. Define
\[
 J_{12}(N)=\prod_{\substack{p\mid N\\p\notin S\\
                    \chi_1(p)=\chi_2(p)=-1}}\frac{p-1}{p-2},
 \qquad
 W_i(N)=\prod_{\substack{p\mid N\\p\notin S\\
                    \chi_i(p)=1,\ \chi_j(p)=-1}}
                      (1-1/p)^{-1},\quad j\ne i.
 \tag{D.1}
\]
Only common inert primes appear in \(J_{12}\). A prime excluded on
just one side does not merge two baseline roots when it divides \(N\).
Local admissibility means that, for every prime \(p\), the target is a
sum of one primitive \(\mathbb Z_p\)-value of each form, together with
the positive real condition.

For a nontrivial finite abelian group \(B\), use **all proper subgroups**
and define
\[
 r_H=[B:H],\qquad d_H=\frac{r_H-1}{2r_H},\qquad K_H=r_H-2,
 \qquad \Phi_B(N)=\sum_{H<B}\frac{\ell_N^{K_H}}{L_N^{d_H}}.
 \tag{D.2}
\]
Set \(\Phi_0=0\). The [endpoint note](PRIMITIVE_NORM_ENDPOINT_CLASS_ERROR.md),
(E.6)--(E.7), proves
\[
 \Phi_B(N)\asymp_B
 \frac{\ell_N^{q-2}}{L_N^{(q-1)/(2q)}}
 \quad(B\ne0),
 \tag{D.3}
\]
where \(q\) is the smallest prime divisor of \(|B|\).

**Theorem MixedFieldClassError.** For every fixed pair of forms and all
integers \(N\ge3\),
\[
 0\le R_{\Lambda_1,\Lambda_2}(N)-R_{P_1,P_2}(N)
 \ll_{Q_1,Q_2}\frac{NJ_{12}(N)}{L_N}
       \bigl(W_1(N)\Phi_{B_1}(N)+W_2(N)\Phi_{B_2}(N)\bigr).
 \tag{D.4}
\]
If the fields agree, both \(W_i=1\), recovering the same-field
endpoint bound. If they differ, then
\[
 W_i(N)\ll_{D_1,D_2,S}\ell_N^{1/4}.
 \tag{D.5}
\]
Consequently the class error is
\[
 \ll_{Q_1,Q_2}\frac{NJ_{12}(N)}{L_N}\ell_N^{1/4}
                     (\Phi_{B_1}(N)+\Phi_{B_2}(N))
       =o_{Q_1,Q_2}(NJ_{12}(N)/L_N).
 \tag{D.6}
\]
For any fixed \(0<\eta<\min_{i:B_i\ne0}\delta_*(B_i)\), it is also
\(O_{Q_1,Q_2,\eta}(NJ_{12}(N)/L_N^{1+\eta})\). If both \(B_i=0\),
the error is exactly zero for every target. All orders and forms are
fixed; no numerical onset, optimal error or varying-field uniformity
is claimed.

## 2. The ideal-class cover remains one-sided

[FormIdealShadow](PRIMITIVE_NORM_ALL_ORDERS.md), (F.2), identifies the
local-shadow condition with the existence of a primitive integral
invertible ideal in \(c_{Q_i}+2G_i\). This remains valid at conductor-
divisible norms. [StabilizerDeficiency](PRIMITIVE_NORM_ENDPOINT_CLASS_ERROR.md),
(E.8)--(E.8a), then gives, for each exceptional value in slot \(i\),
some proper \(H<B_i\) with
\[
 k_{i,H}(a)=\#\{p\in\Omega_{i,H}:v_p(a)=1\}\le K_H.
 \tag{D.7}
\]
Here \(\Omega_{i,H}\) consists of good split primes in field \(i\)
whose inverse-pair class \(g,-g\) has \(2g\notin H\). All primes in
the common \(S\) are omitted. To justify the extra deletion, freeze
those factors of the chosen ideal before forming its class-switch sum.
Conductor-local ideals remain unchanged. The numerical norm need not
determine a unique genus; the chosen compatible ideal suffices.

The fixed-order class-prime Mertens estimate is still
\[
 \sum_{\substack{p\le y\\p\in\Omega_{i,H}}}\frac1p
       =d_H\log\log y+O_{Q_i,S}(1).
 \tag{D.8}
\]
Its precise primary input is Bourgain--Fuchs, arXiv:1105.3698v2,
Section 2 part (5), Lemma 2, equation (5.2), page 15, as retained
and checked in the endpoint source ledger. No subdivision of (D.8)
by the other field's character is made.

## 3. Exact local overlap and a uniform mixed sieve

For a good prime put \(I_i(p)=\mathbf1_{\chi_i(p)=-1}\).
On the integer line \(a+b=N\), the baseline excludes \(a=0\bmod p\)
if \(I_1=1\), and \(a=N\bmod p\) if \(I_2=1\). Its number of roots is
\[
 \rho_N(p)=I_1(p)+I_2(p)-I_1(p)I_2(p)\mathbf1_{p\mid N}.
 \tag{D.9}
\]
Mark only slot 1 for now. A marking prime lies in \(\Omega_{1,H}\),
so \(I_1=0\). In addition to the baseline, a marked prime excludes
\(v_p(a)=1\). Relative to the baseline surviving density, its exact
incremental exclusion probability is
\[
 b_{p,N}=\begin{cases}
 (p-1)/p^2,&\chi_2(p)=1,\\
 1/p,&\chi_2(p)=-1,\ p\nmid N,\\
 0,&\chi_2(p)=-1,\ p\mid N.
 \end{cases}
 \tag{D.10}
\]
Indeed, in the middle case the baseline forbids \(p\) classes modulo
\(p^2\), and the exact-one condition forbids \(p-1\) further classes,
disjoint because \(N\) is a unit. The surviving density is
\((p-1)^2/p^2\), whose ratio to \(1-1/p\) is \(1-1/p\).
In the last case the baseline already forbids every multiple of \(p\),
including all exact-one classes. An additional marking has no effect.
The slot-2 calculation follows by \(a\mapsto N-a\).

**Lemma MixedFieldMarkedSieve.** Let \(\mathcal M\) be any finite
subset of \(\Omega_{1,H}\cap[1,y]\). Write \(S_{\mathcal M}(X;y)\)
for the number of integers in an interval of length \(X\) surviving
the baseline restrictions at good primes up to \(y\), and the exact-one
exclusions at \(\mathcal M\). For \(y=z^{1/16}\),
\[
 S_{\mathcal M}(X;y)\ll_{Q_1,Q_2}
 \frac{XJ_{12}(N)}{\log y}\prod_{p\in\mathcal M}(1-b_{p,N})
            +O_\epsilon(z^{4+\epsilon}).
 \tag{D.11}
\]
The constants are independent of the mark set, target valuations and
root positions.

**Proof.** Adapt the finite quadratic majorant in MixedPrimeSquareSieve,
(M.14). Use modulus \(p\) at an unmarked active prime, and modulus
\(p^2\) at a marked one. Let \(\nu(p)\) be the number of forbidden
classes in this modulus \(m_p\), and \(f(p)=\nu(p)/m_p\).
Always \(0<f(p)\le2/p<1\). At a marked prime the larger count
\(\nu(p)=2p-1\) is possible: the old bound \(\nu(d)\le d\)
must not be reused. For squarefree prime indices,
\[
 \nu(d)\le2^{\omega(d)}d\ll_\epsilon d^{1+\epsilon}.
 \tag{D.12}
\]
CRT gives an interval count \(Xf(d)+O(\nu(d))\).

Set \(w(d)=\prod_{p\mid d}f(p)/(1-f(p))\) and
\(G(z)=\sum_{d\le z}\mu^2(d)w(d)\), with indices supported on the
active primes up to \(y\). The explicit Selberg weights from (M.14)
give quadratic main \(1/G(z)\) and
\( |\lambda_d|\le\prod_{p\mid d}(1-f(p))^{-1}
\ll_\epsilon d^\epsilon\). Thus, choosing the small exponents so
their total is the prescribed \(\epsilon\),
\[
 \sum_{d,e\le z}|\lambda_d\lambda_e|\nu([d,e])
 \le\sum_{d,e\le z}|\lambda_d\lambda_e|\nu(d)\nu(e)
 \ll_\epsilon z^{4+\epsilon}.
\]
This explicitly pays for the enlarged prime-square root counts.

The positive Euler expansion has independent inclusion probabilities
\(f(p)\le2/p\), hence mean logarithmic index at most
\(2\sum_{p\le y}(\log p)/p\le4\log y\), as proved in (M.14).
Markov therefore retains at least three quarters of its mass on
\(d\le z\). The full inverse Euler product is
\[
 E_{\mathcal M,N}(y)=
 \prod_{\substack{p\le y\\p\notin S}}
          (1-\rho_N(p)/p)^{-1}
 \prod_{p\in\mathcal M}(1-b_{p,N})^{-1}.
\]
Its baseline factor equals
\[
 \frac{\displaystyle\prod_{p\le y,p\notin S}
                  (1-(I_1(p)+I_2(p))/p)^{-1}}
      {\displaystyle\prod_{p\le y,p\mid N,p\notin S,
                    I_1(p)=I_2(p)=1}(p-1)/(p-2)}
 \gg_{Q_1,Q_2}\frac{\log y}{J_{12}(N)}.
 \tag{D.13}
\]
Each quadratic inert set has harmonic coefficient one half, so the
numerator is comparable to \(\log y\); its quadratic Euler-log error
is summable. This uses the two separate quadratic Mertens estimates,
not independence. Combine \(G(z)\ge3E_{\mathcal M,N}(y)/4\) with the
quadratic majorant to prove (D.11). \(\square\)

## 4. Random marking, target divisors and the full remainder

Let \(f_p=(p-1)/p^2\). If \(b_{p,N}\ne0\), (D.10) gives
\(b_{p,N}\ge f_p\). If it vanishes, the prime divides \(N\) and is
split on side 1 but inert on side 2. Thus, uniformly for \(0\le t\le1\),
\[
 \prod_{p\le y,p\in\Omega_{1,H}}(1-tb_{p,N})
 \le \prod_{p\le y,p\in\Omega_{1,H}}(1-tf_p)
       \prod_{\substack{p\le y,p\in\Omega_{1,H}\\
                      p\mid N,\ \chi_2(p)=-1}}(1-tf_p)^{-1}
 \ll_{Q_1,Q_2}(\log y)^{-td_H}W_1(N).
 \tag{D.14}
\]
For the last inequality use UniformMarkProduct (E.10)--(E.12) and
\((1-tf_p)^{-1}\le(1-1/p)^{-1}\). One can keep the smaller exact
target factor with \(p\in\Omega_{1,H}\); the displayed \(W_1\)
avoids any joint class-prime distribution assertion.

Set \(z=N^{1/10}\), \(y=N^{1/160}\). Independently mark eligible
primes with probability \(t<1\). A value deficient as in (D.7) survives
with probability at least \((1-t)^{K_H}\). Finite expectation in
(D.11), with \(\epsilon=1/2\), therefore yields
\[
 D_{1,H}(N)\ll(1-t)^{-K_H}
 \left(\frac{NJ_{12}(N)W_1(N)}{(\log y)^{1+td_H}}
                         +N^{9/20}\right).
 \tag{D.15}
\]
The survival factor multiplies both terms. For a positive cap, take
\(t=1-1/\log\log y\) once \(y\ge e^e\). For a zero cap use
deterministic full marking, with \(t=1\) and no division by zero.
Exactly as in (E.15)--(E.16),
\[
 D_{1,H}(N)\ll
 \frac{NJ_{12}(N)W_1(N)(\log\log y)^{K_H}}
                  {(\log y)^{1+d_H}}
          +N^{9/20}(\log\log y)^{K_H}.
 \tag{D.16}
\]
Because \(J_{12},W_1\ge1\), the multiplied remainder divided by the
first term is \(O(N^{-11/20}(\log y)^{1+d_H})\), tending to zero.
It is fully absorbed. Repeat for side 2 and sum the finite proper-
subgroup covers. Replace the logarithms by \(L_N,\ell_N\) and enlarge
the constant over the finite remaining \(N\ge3\). If both doubled
groups are trivial, use exact individual shadow equality instead.
This proves (D.4).

## 5. The quarter-loglog target factor

Suppose \(D_1\ne D_2\). The two quadratic characters and their product
are nonprincipal on a common fixed modulus. The exclusive splitting set
\(\mathcal E_i=\{p\notin S:\chi_i(p)=1,\chi_j(p)=-1\}\) has indicator
\((1+\chi_i)(1-\chi_j)/4\) at good primes. Fixed-modulus prime Mertens
therefore gives
\[
 \sum_{p\le x,p\in\mathcal E_i}\frac1p
       =\tfrac14\log\log x+O_{D_1,D_2,S}(1).
 \tag{D.17}
\]
Only its upper bound is needed. Alternatively this upper bound follows
directly from the fixed Dirichlet Euler products: at real
\(s=1+1/\log x\), the exclusive prime sum is
\(\tfrac14\log(1/(s-1))+O(1)\), since the three fixed nonprincipal
\(L(1,\chi)\)'s are nonzero and prime powers give a bounded error.
For \(p\le x\), replacing \(p^{-s}\) by \(p^{-1}\) costs at most
\((\log x)^{-1}\sum_{p\le x}(\log p)/p=O(1)\).
The tail in the exclusive prime sum is nonnegative and can be discarded.
Thus no uniform growing-modulus theorem enters this estimate.

For large \(N\), split its prime divisors at \(x=\log N\). Small
divisors contribute at most \(\tfrac14\log\log\log N+O(1)\) by
(D.17). There are at most \((\log N)/\log x\) larger prime divisors,
so their reciprocal sum is at most \((\log N)/(x\log x)=O(1)\).
Finally \(-\log(1-1/p)=1/p+O(p^{-2})\), uniformly for these odd primes.
Exponentiating proves (D.5), with small targets covered by the constant.
For equal fields the exclusive sets are empty and \(W_i=1\) exactly.

## 6. Classical lower counts with the actual finite filters

The lower bound uses the classical
[nonnested subgroup theorem](ASYMMETRIC_SUBGROUP_SPLITTING.md),
SubgroupAsymUniform (H2), its fixed-twist estimate (H9), and the positive
unweighting proof of SubgroupSemigroupScale (H5), Section 6. Set
\(H_1=\ker\chi_1\), \(H_2=\ker\chi_2\), both of index two in the
unit group of a common even modulus. Their overlap can be arbitrary.
The additional finite-filter and bad-coefficient interface is proved
here; the underlying weighted estimates retain their existing classical
BGLM and fixed-period box premises.

### 6.1. Filtered split cores and fixed coefficients

Take \(L=8\prod_{p\in S,\ p\text{ odd}}p\). Let \(\mathcal C_i\)
be generated by the good primes split in field \(i\). Fix subsets
\(U_i\subset H_i\), and fixed coprime positive \(S\)-supported
coefficients \(A,B\). Put
\[
 \Gamma_{A,B;U_1,U_2}(M)=
 \frac{L}{|H_1||H_2|}
 \#\{(a,b)\in U_1\times U_2:Aa+Bb=M\pmod L\}.
\]

**Lemma FilteredQuadraticCoreCounts.** On sufficiently large targets
with positive displayed finite density,
\[
 \#\{Ax+By=M:x\in\mathcal C_1,\ y\in\mathcal C_2,
                      \ x\bmod L\in U_1,\ y\bmod L\in U_2\}
 \asymp_{L,U_1,U_2,A,B}\frac{M J_{12}(M)}{\log M}.
 \tag{D.18}
\]
Zero finite density gives exactly zero count. All fixed data are held
constant, and the estimate is uniform in target valuations.

**Proof.** Use weights \(V_i=\tau\mathbf1_{\mathcal C_i}\) and their
normalized shifted families, as in (H2). To retain residual divisors
\(u\in\mathcal C_1\), \(v\in\mathcal C_2\) in the full positive
unweighting family, consider \(Au x+Bv y=M\), with its **current**
masks \(V_1\subset H_1,V_2\subset H_2\). Insert on the original
first unit variable the exact filter
\[
 w(x)=\mathbf1_{x\bmod L\in V_1}
       \mathbf1_{B\mid M-Au x}
       \mathbf1_{z\bmod L\in vV_2},\qquad z=(M-Au x)/B.
\]
This has fixed period \(BL\). It forces \(\chi_2(z)=1\), so the
second paired-divisor quotient also has \(\chi_2=1\), because its
removed kernel indices and \(v\) lie in \(H_2\). This is the precise
nonnested quotient condition (H11), now with the actual two masks.
Positive finite density also forces \((M,AB)=1\), by reducing the
equation at every prime dividing one of the coprime bad coefficients.

Expand \(w\) in characters modulo \(BL\). Only the two first-side
twists \(1,\chi_1\) retain a principal mean; for every other fixed
twist, (H9) has zero principal mean and a power-saving remainder.
Keep every kernel phase and every state, including nontrivial genus
states. This is the proof of (M.25)--(M.26) in
[the maximal-order note](PRIMITIVE_NORM_MAXIMAL_ORDERS.md), with
\(\chi_1\) on the first side and the independent \(\chi_2\) test in
the second quotient. Its weighted main is
\(M\Gamma_{Au,Bv;V_1,V_2}(M)\rho_{T_1}\rho_{T_2}/(ABuv)\)
times the same normalized good local factors. This current finite
density need not equal \(\Gamma_{A,B;U_1,U_2}(M)\) when coefficients
or masks are changed independently.

The period normalization has no extra factor \(B\). Since the prime
divisors of \(B\) already divide \(L\), both lifted kernels modulo
\(BL\) have size \(B\varphi(L)/2=B|H_i|\). Every solution pair
modulo \(L\) has exactly \(B\) lifts satisfying the equation modulo
\(BL\), since \((A,B)=1\) and \(u,v\) are good units. Thus
\(BL\) divided by the two lifted kernel sizes, multiplied by the lifted
root count, is precisely the current \(\Gamma_{Au,Bv;V_1,V_2}(M)\).

For fixed bad coefficients, the BGLM parameters are exactly those
checked after (M.26): \(a_1=M\), \(a_2=Au j c_0\), \(a_3=1\),
\(d=v k a_{\rm aux}\), \(c=BL\). The actual v5 Proposition 5.1,
printed page 15, does not require \((a_2,c)=1\). Its varying modulus
and short variables are \(BL\)-units; the added fixed factor \(A\)
preserves its required coprimalities after the existing gcd extractions.
The fixed-period box estimates allow this fixed enlarged period. The
small-power parameter can be decreased after the weighted error saving
has been chosen, just as in (M.25).

During every gcd or cofactor removal, transport the masks by inverse
multiplication. For restored values \(x=du x'\), \(y=dv y'\), the
residual masks are \((du)^{-1}U_1\), \((dv)^{-1}U_2\); all these
factors are units modulo \(L\). The resulting map of finite solution
pairs is a bijection. When \(d\mid M\), the exact identity is
\[
 \Gamma_{Au,Bv;(du)^{-1}U_1,(dv)^{-1}U_2}(M/d)
       =\Gamma_{A,B;U_1,U_2}(M).
 \tag{D.18a}
\]
Multiply the residual coordinates by \(du,dv\) and its equation by
\(d\) to obtain this identity. Keeping the original masks or the base
density fixed before this simultaneous transport would be invalid.
For example, modulo 8 take \(A=1,B=2\), \(U_1=\{1\}\subset\ker\chi_{-4}\),
\(U_2=\{3\}\subset\ker\chi_{-8}\), and \(u=5,v=1\). At \(M=3\),
the shifted equation with unchanged masks has a solution, while the
base equation has none. Common cofactors belong to
\(\mathcal C_1\cap\mathcal C_2\), and side-specific ones remain in
their respective cores.

Now apply the positive finite-block unweighting in (H5), Section 6, to
this full filtered shifted family. Its rough dimension is two. If
\(c=|H_1\cap H_2|/|G_0|\), restoration contributes
\(2c+(1/2-c)+(1/2-c)=1\), irrespective of overlap. Its exact target
compensation leaves \((p-1)/(p-2)\) only at primes forbidden by both
sides, hence precisely \(J_{12}\). The finite masks preserve the
positive normalization and only reduce the upper bounds for residual
weights; their positive densities have a positive minimum. All cofactor
supports and summed errors remain within the same small-power budget.
This proves the lower bound in (D.18). The unmarked version of (D.11)
gives the matching upper bound. \(\square\)

### 6.2. A finite cover of local targets

**Theorem MixedFieldValueCounts.** For every fixed pair of forms there
are constants \(c,C>0\) and a threshold such that every sufficiently
large locally admissible target satisfies
\[
 c\frac{NJ_{12}(N)}{\log N}
 \le R_{P_1,P_2}(N)\le R_{\Lambda_1,\Lambda_2}(N)
 \le C\frac{NJ_{12}(N)}{\log N}.
 \tag{D.19}
\]
Every nonadmissible target has zero count.

**Proof.** The form-specific local-coset argument in
[the all-orders note](PRIMITIVE_NORM_ALL_ORDERS.md), (O.9) and Section
4.1, gives full local value cosets modulo \(p^{e_p}\), where
\(e_p=1+2\max_i v_p(\Delta_i)\) suffices for both forms. Given a
local decomposition, perturb its two values by \(\pm p^{e_p}t\).
One can make both valuations at most \(e_p+1\): modulo \(p^2\), at
most two values of \(t\) would create valuation at least \(e_p+2\),
and \(p^2>2\). Values of valuation below \(e_p\) are unchanged in
valuation. The full-coset property preserves primitive local solubility.

Therefore a finite collection of \(S\)-supported prefixes \(r,s\),
whose prime exponents are at most \(e_p+1\), covers all admissible
targets. FormGenusCoreFilters (F.4) for each form and the common \(S\)
gives the actual masks \(U_{r,1}\subset H_1\),
\(U_{s,2}\subset H_2\). Divide by \(d=(r,s)\): put \(A=r/d\),
\(B=s/d\), \(M=N/d\). Local admissibility supplies a profile with
\(d\mid N\) and positive finite density in (D.18). The map
\((x,y)\mapsto(rx,sy)\) injects these numerical core pairs into the
shadow count. Since \(d\) is \(S\)-supported,
\(J_{12}(N/d)=J_{12}(N)\). The prefixes are from a fixed finite
collection, so their positive constants and thresholds can be made
uniform. This proves the shadow lower bound. No infinite prefix tail
is truncated in this argument.

Dropping all bad-prime conditions and using the unmarked mixed inert
sieve proves the shadow upper bound. Finally (D.6) is smaller than that
lower scale, so subtracting the class error gives the actual lower
bound. Local necessity gives zero outside admissibility. \(\square\)

**Corollary MixedFieldGenusInvariance.** If \(Q_i'\) lies in the same
genus as \(Q_i\) within its own discriminant, then on their common
sufficiently large admissible targets,
\[
 R_{P_1,P_2}(N)=R_{P_1',P_2'}(N)
 \left(1+O_{\mathbf Q}(W_1\Phi_{B_1}+W_2\Phi_{B_2})\right).
 \tag{D.20}
\]
The same relative error compares actual and shadow counts. Subtract
(D.4) for both pairs, whose shadows agree slot by slot, and divide by
(D.19). If both doubled groups are trivial, equality holds at every
target without division. No genus is identified across discriminants.

The separate pointwise leading-constant theorem (B.4) in
[the binary-correlation note](PRIMITIVE_NORM_BINARY_CORRELATION.md)
allows nonnested subgroups, but uses its explicitly recorded ZF and graph
G imports. It is unnecessary for (D.4). Any use of it for a leading
asymptotic must retain those premises and the bad-prefix tail argument;
no new leading coefficient or quantitative leading remainder is claimed
here.

## 7. The discriminant minus 23 and Gaussian pair

Take \(Q_1=x^2+xy+6y^2\) and \(Q_2=x^2+y^2\). The retained
[class-number-three proof](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md)
gives \(G_1=C_3\); the Gaussian class group is trivial. Hence \(B_1=C_3\),
\(B_2=0\), \(q_1=3\), and (D.6) reads
\[
 R_{\Lambda_1,\Lambda_2}(N)-R_{P_1,P_2}(N)
 \ll \frac{NJ_{12}(N)\ell_N^{5/4}}{L_N^{4/3}},
 \quad
 J_{12}(N)=\prod_{\substack{p\mid N,\ p\ne2,23\\
                  \chi_{-23}(p)=\chi_{-4}(p)=-1}}\frac{p-1}{p-2}.
 \tag{D.21}
\]

Every positive integer is an additive locally admissible target for
this pair. At 2 the first field is split because \(-23\equiv1\bmod8\),
so its primitive norm represents every 2-adic value. At 23, Gaussian
primitive norms represent every unit; choose the first form's value 1
or 4 with residue different from the target. Both are primitive local
values at 23. At any other odd prime, a split side represents all local
values, and two inert sides require a sum of two units, always possible
for an odd prime. Positive real solubility is immediate. This is local
admissibility, not an effective global cutoff. The classical theorem
(D.19) now supplies \(NJ_{12}(N)/\log N\)-order numerical-value counts
for every sufficiently large integer, and (D.20) gives relative class
error \(O(\ell_N^{5/4}/L_N^{1/3})\). Eventual separately primitive
binary-pair existence already has a classical quaternary route in this
repository; the assertions here concern numerical-value counts and
quantitative class/genus comparison.

## 8. Evidence and limits

The [exact verifier](verify_primitive_norm_different_fields.py) enumerates local prime-square residues,
checks the three overlap probabilities, baseline root collisions,
finite marking expectations and target-factor inequalities, and tests a
small CRT example against direct enumeration. Deliberate controls reject
the old \(\nu(d)\le d\) estimate, a spurious marker gain at a target-
dividing opposite-inert prime, and an incorrect common-inert product.
These finite checks do not prove distribution statements or asymptotics;
those are the written argument with the explicit retained premises.

Run from the repository root:

    python3 research/arithmetic/verify_primitive_norm_different_fields.py

The [receipt](primitive_norm_different_fields_checks.json) records the
actual finite replay. The [source ledger](primitive_norm_different_fields_source_ledger.json)
pins the proof, checker and retained mathematical premises, distinguishing
classical count comparisons from the separate leading-asymptotic inputs.

No optimality, varying-order uniformity, numerical cutoff, independent
external review, formal verification or historical priority is asserted.
