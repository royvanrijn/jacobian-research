# Proof audit: pair-primitive binary sums

Reviewed revision: 841f6e73699bee3d5fafe657b4ff0aec23a8bbf8
(18 November 2026). This is an internal proof reconstruction and
primary-source application audit. It is not external human review or
formal verification. Deep classical inputs remain imported theorems.
Pre-existing worktree changes, including the separate -11 specialization,
were preserved; that specialization is outside this audit. A separate
-8 completion appeared during this work; its proof and full finite
certificate were subsequently audited in the addendum below.

The [canonical proof](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md) contains the
proved improvements. Equation references retain its original numbering.

| Investigation | Finding and status |
|---|---|
| General asymptotic | Supported by the exact all-index source estimate; no failed hypothesis identified in (13)--(20). |
| Exponent | The original Rouse route gives every eta<1/18; the subsequent Waibel application now gives every eta<1/8, including1/9; see the novelty comparison. |
| Normalization | \(4\pi^2/\sqrt{\det A}\) is correct for \(Q=x^{\mathsf T}Ax/2\); the displayed Rouse Section 3 prefactor is inconsistent with that convention. |
| Finite local test | Valid at arbitrary target valuations; the maximum determinant valuation can even be replaced by the minimum. |
| Effective -8 specialization | The class number, normalization, dyadic factorization and normalizer average reconstruct correctly under the identified classical spectral imports. |
| Exact -8 classification | Open in the reviewed revision; the contemporaneous exact proof now passes the proof/completeness audit and full finite replay described below. |
| Stronger consequences | Same-field linear and distinct-field square-root log-log coordinate bounds; subsequent corrected-Henriot application strengthens value-pair counts to N/log^(log4+epsilon)N and gives the NJ12/logN upper bound. |
| Novelty | See the thorough follow-up comparison: classical shell counts and methods, a modest arbitrary-pair deduction, and exact -8 classification as the main unverified-priority candidate. |

## Exact analytic input and hypothesis audit

The source is Rouse,
[*Integers represented by positive-definite quadratic forms and Petersson inner products*,
arXiv:1802.03437v1](https://arxiv.org/pdf/1802.03437v1).
The input is not its representation conclusion, Theorem 1.

Theorem 3, PDF p.10, applies to every primitive positive-definite
integer-valued quaternary \(T=x^{\mathsf T}Ax/2\), with integral
even-diagonal Gram matrix, determinant \(\Delta\), and level \(L\).
It states
\[
 \langle C_T,C_T\rangle\ll_\epsilon
 \max\{(L^2\Delta)^{1/4+\epsilon},L^{1+\epsilon}\}.
\]
No squarefree level, fundamental discriminant, class number, shape
condition or target-coprimality hypothesis is imposed.

Section 5, PDF pp.15--16, uses the general coefficient bound of
Schulze-Pillot--Yenirce. In their retained
[arXiv:1602.01803v6](https://arxiv.org/pdf/1602.01803v6)
it is Theorem 12, PDF p.12, rather than Theorem 10 in Rouse's
earlier citation. It applies to any integral-weight
\(C\in S_k(\Gamma_0(L),\chi)\), including oldforms, arbitrary level
and arbitrary character, at every positive coefficient index.
For weight two it gives
\[
 |c(n)|\ll_\epsilon
 \langle C,C\rangle^{1/2}L^{1+\epsilon}\tau(n)\sqrt n,
\]
using \(\dim S_2(\Gamma_0(L),\chi)\ll_\epsilon L^{1+\epsilon}\).
The Petersson normalizations differ by the fixed factor \(3/\pi\).
Combining the two bounds gives exactly (8), after renaming arbitrarily
small epsilon parameters. Neither local solubility nor
\((n,L)=1\) is a hypothesis of this coefficient bound. The sharper
coprime-index estimate and its ineffective newform norm lower bound
are not used.

The v1 Theorem 1 strong-local-solubility bullet has a different first
level exponent from its proof on p.16. The present argument does not
import that representation threshold; its precise input is Theorem 3
plus the general all-index coefficient bound.

For \(g=(d,e)\), \(u=d/g\), \(v=e/g\), first require \(g^2\mid N\).
The divided form \(u^2Q_1\oplus v^2Q_2\) is primitive, because its
coefficient content is \(\gcd(u^2,v^2)=1\). Its determinant is
\(\det A_1\det A_2(uv)^4\), its level divides \(L_0u^2v^2\),
and its target is \(N/g^2\). These satisfy all source hypotheses.
The bound is uniform in \(d,e,N\); its constant depends only on the
fixed forms and selected epsilon parameters. A source parameter \(h\)
gives at most \((uv)^{3+6h}\sqrt N/g\), before the divisor loss.

The Möbius tail keeps both blocks nonzero. For a fixed \(d>D\),
there are \(O(N/d^2)\) positive first-block values; each binary
coordinate multiplicity and the divisor count of the other block gcd
cost only an arbitrarily small power of \(N\). Summing \(d^{-2}\)
proves (12). Zero blocks would invalidate this tail and are correctly
handled only after truncation, at cost \(O(D^2N^\epsilon)\).

The Jacobian in (15) is correct: scaling a two-coordinate block by \(p\)
changes Haar measure by \(p^{-2}\). In particular
\[
 \lambda_p(1,1;N)
 ={\bf1}_{p^2\mid N}p^{-2}\beta_p(F;N/p^2).
\]
Primitive binary densities are at most 2 at good primes; (18) is
therefore uniform in target valuations. Nonnegative Euler factorization
and the absolutely convergent weighted sum justify (20); no
unrestricted density is divided out.

## Real normalization

The four-dimensional ellipsoid \(x^{\mathsf T}Ax/2\le n\) has volume
\(2\pi^2n^2/\sqrt{\det A}\), and derivative
\(4\pi^2n/\sqrt{\det A}\). This pairs with the stated finite densities.

For a separate exact challenge, four squares at \(n=1\) have
\(\det A=16\), \(\beta_2=1\), and odd-prime product \(8/\pi^2\).
The corrected leading product is 8, the actual representation number.
The prefactor displayed in Rouse's Section 3 would instead give 2.
It cannot be copied under these definitions. This discrepancy does
not affect a big-O coefficient estimate.

## 1. Weighted Euler product: alpha < 1 is valid

Let

\[
 H_p(\alpha,N)=\lambda_p(0,0;N)+p^\alpha(\lambda_p(1,0;N)+\lambda_p(0,1;N))+p^{2\alpha}\lambda_p(1,1;N).
\]

At every good prime, (17)--(18) hold with absolute constants uniformly in the valuation of N. If p does not divide N, then lambda_p(1,1;N)=0 exactly: both divisible blocks force p^2 | F. Therefore

\[
 H_p(\alpha,N)=1+O(p^{-2+\alpha}) \quad(p\nmid N).
\]

This is absolutely summable over primes whenever alpha < 1. At target primes,

\[
 H_p(\alpha,N)\le 1+C(p^{-1}+p^{-2+\alpha}+p^{-2+2\alpha})
 \le 1+C_\alpha p^{-\delta},\qquad \delta=\min(1,2-2\alpha)>0.
\]

For each fixed epsilon > 0, choose P so large that log(1+C_alpha p^{-delta}) <= epsilon log p for p>P. The product of these factors over p|N, p>P, is at most rad(N)^epsilon <= N^epsilon. The finitely many p<=P contribute a constant depending only on alpha and epsilon. Bad primes are a fixed finite set, bounded by the unrestricted-density bound. This proves

\[
 \sum_{d,e\ge1}\mu(d)^2\mu(e)^2(de)^\alpha E_{d,e}(N)
 \ll_{Q_1,Q_2,\alpha,\epsilon}N^{1+\epsilon}
\]

for every fixed 0<alpha<1, uniformly in N and all target valuations. One must change the existing sentence describing the target-prime error as O(1/p): for alpha>1/2 the double-divisibility error is generally only O(p^{-2+2alpha}), but its product over target primes remains N^{o(1)}.

Absolute convergence and nonnegativity justify the Euler factorization and completion. On d>D or e>D, (de)^alpha >= D^alpha, giving tail O(N^{1+epsilon}D^{-alpha}). The endpoint alpha=1 does not follow: at nontarget primes p^alpha lambda_p(1,0) is of order 1/p, and the positive absolute Euler product diverges.

## 2. Exponent optimization

For g=(d,e), a=d/g, b=e/g, the primitive content reduction has level at most L_0(ab)^2 and determinant D(F)(ab)^4. Applying (8) with a small source parameter u gives coefficient factor at most (ab)^{3+6u}; the second max term has exponent 3+2u. With tau(n) <<_v n^v and n=N/g^2 this gives, more sharply,

\[
 |c_{d,e}(N)|\ll (de)^{3+6u}g^{-7-12u-2v}N^{1/2+v}.
\]

Dropping the favorable g power and summing d,e<=D yields O(N^{1/2+v}D^{8+12u}). Thus the claimed candidate budget is rigorous using the imported, hypothesis-checked (8):

\[
 N^{1+\epsilon}D^{-1}+N^{1/2+\epsilon}D^{8+O(\epsilon)}+D^2N^\epsilon+N^{1+\epsilon}D^{-\alpha}.
\]

Put D=N^d. Neglecting the arbitrarily small epsilon losses, the relevant saving is min(alpha*d, 1/2-8d). For fixed alpha it is optimized by d=1/(16+2alpha), with saving alpha/(16+2alpha). Its supremum as alpha approaches 1 is 1/18. Hence every eta<1/18 is available from this budget; the endpoint is not established.

A concrete eta=1/20 choice is alpha=15/16, D=N^{1/18}. The four pre-epsilon savings are 1/18, 1/18, 8/9, and 5/96, respectively. Their minimum 5/96 is strictly greater than 1/20, so fixed sufficiently small u,v,epsilon absorb all losses. This uses no improvement of the imported theta estimate.

## 3. Bad-prime criterion and a smaller admissibility modulus

The stated max t_p is safe. In fact, because BOTH blocks are primitive, the derivative valuation is at most

\[
 u_p=\min\{v_p(\det A_1),v_p(\det A_2)\}.
\]

Each block separately has a derivative component with valuation bounded by its own determinant valuation. Consequently congruences modulo p^{2u_p+1} suffice, and the Hensel density lower bound improves to p^{-5u_p-3}. No target valuation enters.

To check the exponent in the lower bound, fix the other three coordinates modulo p^{2s+1}, where one selected derivative has exact valuation s. For h>=2s+1 there are p^{3(h-2s-1)} choices of those three coordinates, each with p^s roots of the remaining congruence near the chosen exact root. Dividing by p^{3h} gives p^{-5s-3}. Primitive residue classes are preserved. The recursion (7) supplies a uniform unrestricted upper bound from a bounded primitive density and the geometric factor sum p^{-2j}. No failure of the local criterion or its uniformity was found.

## 4. Exact mixed target ratios

Let a_p be the good-prime factor for p not dividing N. Direct factorization of (5) gives the exact target ratio b_p/a_p:

\[
 T_p=\begin{cases}
 (p+3)/(p+2),&\varepsilon_1=\varepsilon_2=1,\\
 (p-1)/(p-2),&\varepsilon_1=\varepsilon_2=-1,\\
 1-1/p,&\varepsilon_1=-\varepsilon_2.
 \end{cases}
\]

For the mixed-sign case, a_p=1-p^{-2} and b_p=(1-p^{-1})(1-p^{-2}). Put C_good=product_{p notin S} a_p>0. Then exactly

\[
 S^*(N)=C_{\rm good}\prod_{p\in S}\beta_p^*(N)
                         \prod_{p\mid N,\ p\notin S}T_p.
\]

If det A_1 det A_2 is a square, all good-prime ratios exceed one. This covers arbitrary forms in the same imaginary quadratic field, including different ideal classes and orders, not only self-sums. For every admissible N,

\[
 S^*(N)\ge C_{\rm good}\prod_{p\in S}p^{-5u_p-3}=:c_{\rm loc}>0.
\]

Under the asymptotic, sufficiently large admissible N therefore satisfy

\[
 R^*(N)\ge (Jc_{\rm loc}/2)N.
\]

This is an explicit positive lower coefficient; its numerical onset still needs an explicit constant in the analytic error.

## 5. Sharp uniform lower order for genuinely different fields

When det A_1 det A_2 is nonsquare, let chi be the associated nonprincipal real quadratic character. Exact ratios and positivity of the other target ratios give

\[
 S^*(N)\gg_{Q_1,Q_2}\prod_{p\mid N,\ \chi(p)=-1}(1-1/p).
\]

The Mertens theorem for arithmetic progressions implies

\[
 \prod_{p\le x,\ \chi(p)=-1}(1-1/p)\asymp_\chi (\log x)^{-1/2}.
\]

Primary source inspected: Kenneth S. Williams, *Mertens' Theorem for Arithmetic Progressions*, Journal of Number Theory 6 (1974), 353--359, Theorem 1, printed p.356, author-hosted PDF https://people.math.carleton.ca/~williams/papers/pdf/057.pdf . Its fixed-modulus theorem applies to each reduced residue class with chi=-1; there are exactly phi(q)/2 such classes.

Take L=log N. The product over target primes p<=L is at least the full product over negative-character primes <=L. For p>L,

\[
 \sum_{p\mid N,p>L}p^{-1}\le \frac{\omega_{>L}(N)}L
 \le \frac{\log N}{L\log L}=1/\log L.
\]

Using log(1-1/p)>=-2/p shows that this latter product is bounded below by a positive absolute constant for large N. Therefore

\[
 S^*(N)\gg_{Q_1,Q_2}(\log\log(3N))^{-1/2},\qquad
 R^*(N)\gg_{Q_1,Q_2}\frac{N}{\sqrt{\log\log(3N)}}
\]

for sufficiently large admissible N, using the main asymptotic.

This uniform order cannot generally be improved to cN. Fix an admissible residue a modulo the finite admissibility modulus M. Let P_x be the product of all good negative-character primes <=x. Choose 1<=c_x<=M with c_x P_x=a mod M, and N_x=c_x P_x. All N_x are admissible; c_x contributes only boundedly many possible local factors. Exact ratios and the Mertens product give S^*(N_x) asymp (log x)^(-1/2). The prime number theorem in arithmetic progressions gives log N_x asymp x, hence S^*(N_x) asymp (loglog N_x)^(-1/2). In fact, to refute the linear bound it already suffices that the negative-character Mertens product tends to zero; the PNT comparison is only used to identify the sharp envelope.

## 6. Stronger elementary conversion to distinct value pairs

Let f_K(n) count integral ideals of norm n in an imaginary quadratic field K. Its Dirichlet series is zeta(s)L(s,chi_K). The norm embedding already used for (11) gives

\[
 r_Q(n)\le w_K f_K(4an),\quad Q=ax^2+bxy+cy^2,
\]

where w_K is the finite unit count. The map from coordinate solutions to generators is injective.

At a split prime, f_K(p^j)=j+1, whereas the coefficient of p^{-js} in (zeta(s)L(s,chi_K))^4 is binom(j+7,7). The inequality

\[
 (j+1)^3\le\binom{j+7}{7}
\]

follows by induction, since (j+8)(j+1)^2-(j+2)^3=4j^2+5j>=0. At inert primes the local factor has only even exponents and coefficients one; taking its fourth power only increases those coefficients. At ramified primes the same domination is immediate. By multiplicativity,

\[
 f_K(n)^3\le f_K^{*4}(n).
\]

The classical ideal-count bound A_K(X)=sum_{n<=X}f_K(n)<<_K X is elementary here: f_K=1*chi_K, and Dirichlet's hyperbola method plus bounded character partial sums gives A_K(X)=L(1,chi_K)X+O_K(sqrt X). If A_K(X)<=C_K X, then partial summation gives sum_{n<=X}f_K(n)/n<=C_K(1+log X). Therefore

\[
 \sum_{n\le X} f_K^{*4}(n)
 \le C_KX\left(\sum_{n\le X}\frac{f_K(n)}n\right)^3
 \le C_K^4X(1+\log X)^3.
\]

Replacing n by 4an and bounding the resulting subsequence by the full sum proves

\[
 \sum_{n\le X} r_Q(n)^3\ll_Q X(1+\log X)^3.
\]

Let V(N) count distinct ordered pairs (a,b) with a+b=N and each value primitively represented by its respective form, and let w(a)=r_{Q_1}^*(a)r_{Q_2}^*(N-a). Hölder gives

\[
 R^*(N)\le V(N)^{1/3}\left(\sum_a w(a)^{3/2}\right)^{2/3}.
\]

Cauchy and the third moments show sum_a w(a)^{3/2} << N(1+log N)^3. Hence

\[
 V(N)\gg_{Q_1,Q_2}\frac{R^*(N)^3}{N^2(1+\log N)^6}
 \gg_{Q_1,Q_2}\frac{N S^*(N)^3}{(1+\log N)^6}
\]

on sufficiently large admissible targets. Consequently:

* Same-field pairs, including self-sums: V(N) >> N/(log N)^6.
* Distinct-field pairs: V(N) >> N/((log N)^6 (loglog N)^(3/2)).

These are stronger than N^(1-epsilon), while still far from the expected/logarithmic-scale value-pair asymptotic. They require no theorem on shifted divisor correlations. The same argument applies to a fixed positive-length balance interval once the manuscript's angular corollary is accepted; the third-moment upper bound applies to any subset of summand indices. The lower constant then also depends on the interval.


## Independent -8 spectral reconstruction

### Verdict

No false step was identified in (32)--(38). The class number statement, normalization, dyadic factorization, and normalizer average can be proved directly. The only substantial imports for the spectral bound are strong approximation for the norm-one group, global Jacquet--Langlands with the archimedean weight identification, and weight-two Ramanujan at unramified primes. Their hypotheses are satisfied even though the dyadic order is nonmaximal. Nonmaximality affects the transferred conductor at 2, whereas all operators to which Ramanujan is applied are at odd primes different from the reduction prime p.

The reviewed revision's assertion that all other constituents transfer was compressed; the canonical proof now supplies the bridge directly. Removing the two finite-group characters by itself does not prove this assertion: one must use strong approximation to rule out an automorphic constituent having a one-dimensional factor at any split place without being a global norm character. The earlier Eisenstein lemma contains that argument; it is reconstructed for the present order below. This was a compressed imported proof dependency, now expanded, rather than a failed hypothesis.

### 1. The order and its right ideal classes

Put iota=I+J and j=K in the Hamilton algebra B. As (I+J)K=I-J, the order generated by 1,iota,j,iota*j is exactly
O={a+bI+cJ+dK in Z^4: b=c mod 2}.
Its index in the Lipschitz order is 2, hence its index in the Hurwitz maximal order H is 4. In coordinates u+v*j, with u=a+b*iota and v=c+d*iota, its reduced norm is a^2+2b^2+c^2+2d^2.

Here is a proof of h(O)=1 for locally principal right ideals. H is norm Euclidean: round the four coordinates to integers with errors r_i of magnitude at most 1/2. If s=sum r_i^2, the squared distance to a suitable half-integer vector is s-sum |r_i|+1. Since sum |r_i|>=2s, the minimum of the integer and half-integer distances is at most 1/2. Thus every right H-ideal is principal.

At 2, additive volume normalized by vol(H_2)=1 gives vol(H_2^x)=3/4. The half-integer coset always has odd norm, and exactly half of the integer coset does. Also vol(O_2)=1/4 and half of O_2 has odd norm, so vol(O_2^x)=1/8. Multiplication by these local units preserves additive Haar measure. Consequently [H_2^x:O_2^x]=6. H has 24 global units, of which exactly +/-1,+/-K lie in O_2^x. Its orbit on H_2^x/O_2^x therefore has size 6 and is transitive.

For a locally principal right O-ideal I, principalize IH=alpha*H and replace I by alpha^-1 I. Away from 2 the resulting ideal equals H_l=O_l. At 2 it is x O_2 for x in H_2^x. Write x=h u with h in H^x and u in O_2^x, using transitivity. Multiplying on the left by h^-1 makes every localization equal to O_l. Hence I was principal. This proves the exact class-number statement needed, without a mass formula and without assuming O maximal.

### 2. The finite operator and its normalization

For inert p>=5, B splits at p and O/pO is M_2(F_p). The projective image Gamma of the four O-units has order 2; only +/-1 are scalar. Thus Y=Gamma\G has p(p^2-1)/2 elements.

The sum
T_n f(Gamma g)=1/4 sum_{alpha in O,nrd alpha=n} f(Gamma rho(alpha)g)
is well-defined: replacement of g by rho(u)g, u in O^x, permutes the summation by alpha -> alpha u. It commutes with the right G action. Quaternion conjugation permutes the shell and is projective inversion, so T_n is self-adjoint for counting measure. This measure also agrees up to scale with the invariant measure on the adelic quotient, because the right G action is transitive.

For odd n coprime to p, all integral O-ideals of norm n are locally principal and are principal by Section 1. Their generators have norm n and are unique modulo the four right units. Principalization acts by alpha^-1 at p; conjugation replaces this by the concrete alpha sum. This proves the operator identity including both direction of action and 1/4 normalization.

The number of such ideals is sigma(n), because O_l=M_2(Z_l) at every l dividing n and there are 1+l+...+l^a integral right ideals of norm l^a. Accordingly T_n 1=sigma(n)1. The determinant character chi(Gamma g)=(det(g)/p) is well-defined and satisfies T_n chi=(n/p)sigma(n)chi.

The operator includes imprimitive ideals. It is not just the primitive sphere correspondence. In particular
T_(l^(a+2))=T_l T_(l^(a+1))-l T_(l^a)
on this projective rational quotient. The central scalar ideal acts trivially.

### 3. Exact analytic imports and hypothesis check

Let G_B=PB^x, let V_l be the image of O_l^x, and replace V_p by the kernel of reduction modulo p to obtain U(p). Class number one identifies
G_B(Q)\G_B(A_f)/U(p) = Gamma\G.
The stabilizer is exactly the projective image of O^x: if a rational quaternion lies projectively in every V_l, local scalar corrections have even norm valuations and can be combined into one rational scalar, after which the quaternion lies in O_l^x at every finite l.

The strong approximation input is Voight, Quaternion Algebras, Main Theorem 28.5.3 (the author PDF retained in the retained arithmetic source directory). It applies to B^1 away from {infinity,l} for any odd split place l: B is split there, so the S-indefinite/noncompact hypothesis holds.

An irreducible finite-dimensional smooth representation of GL_2(Q_l) factors through determinant: its kernel is open, and normality plus diagonal conjugation forces the kernel to contain the upper and lower unipotents. If an automorphic constituent has a one-dimensional local factor at a split place l, strong approximation therefore makes it invariant under B^1(A); it is a global norm character. The restriction of such a character to the transitive finite quotient is a character of PGL_2(F_p). There are exactly two: the trivial and the quadratic determinant characters, since elementary unipotents generate SL_2(F_p). Their realizations on Y are 1 and chi. Thus the complement of their span contains no automorphic constituent with a finite-dimensional local factor at a split place.

The applicable primary source is Jacquet--Langlands, Automorphic Forms on GL(2), Theorem 14.4, together with Section 5, the n=0 archimedean case around Proposition 5.20, for the infinite local correspondence. The theorem asks for infinite-dimensional local factors at split places; the preceding paragraph supplies this hypothesis. The trivial definite infinite type transfers to weight two. No maximal-order assumption at 2 occurs in this transfer statement. See the IAS source, theorem on printed p.248 (PDF p.254 in the currently accessible 2026 typesetting):
https://publications.ias.edu/sites/default/files/automorphic-forms-on-gl2_rpl_9.pdf

At every odd l dividing n, l != p, the transferred factor is unramified because V_l is maximal. Its two unnormalized Hecke roots have absolute value sqrt(l). This is the classical weight-two Ramanujan bound; the pinned references in the earlier proof are Deligne, Formes modulaires et representations l-adiques, Theorems 4.9 and 5.1, combined with La conjecture de Weil I, Theorem 1.6:
https://www.numdam.org/item/SB_1968-1969__11__139_0.pdf
https://www.numdam.org/item/PMIHES_1974__43__273_0.pdf

The all-ideal recurrence then bounds the l^a eigenvalue by (a+1)l^(a/2); multiplicativity gives ||T_n|| <= tau(n)sqrt(n) on {1,chi}^perp. The implied constant is exactly 1 and is independent of the reduction prime p and the odd index n. Voight Propositions 41.3.1 and 41.3.6 explicitly give the coprime factorization and all-ideal recurrence; their local maximality hypothesis is only needed at the primes of the index.

This verifies the imported statements as applications of classical theorems. It does not reconstruct the proofs of global Jacquet--Langlands or Ramanujan.

### 4. Dyadic factorization

Let N=2^e m with m odd and p not dividing N. The operator A_e=T_(2^e) is still the exact quaternion-shell sum. Its row sum is c_e=r(2^e)/4, equal to 1,2,6,6 for e=0,1,2,3; positivity and self-adjointness give ||A_e||<=c_e. Both 1 and chi are eigenfunctions, so their orthogonal complement is A_e-invariant.

For q in O of norm N, form the right O-ideal J whose localizations are O_2 at 2 and qO_l at odd l. It contains qO, is integral and locally principal, and has norm m. It is unique with these local specifications. Principalize J=beta O with nrd beta=m. Then q=beta alpha for alpha in O of norm 2^e. Exactly four pairs arise, by replacing beta with beta u and alpha with u^-1 alpha for u in O^x. Constructing the other superideal (qO_2 at 2, O_l at odd l) gives the reversed factorization. Hence the convolution sums satisfy
T_N=A_e T_m=T_m A_e.
This does not assume unique factorization of elements at the nonmaximal dyadic place. The complementary operator norm is bounded by c_e tau(m)sqrt(m).

### 5. Normalizer average

Write C=(F_(p^2)^x)/(F_p^x). Its image is a nonsplit torus of order p+1. The projective image of j acts by the nontrivial Galois involution, so K=C union C rho(j) is its normalizer of order 2(p+1). For q=u+v j invertible modulo p, its image lies in K exactly when u=0 or v=0. This is precisely the event that one binary coordinate block vanishes. Gamma lies in K.

For f=1_(Gamma\K), its mean on Y is mu=2/[p(p-1)]. The norm map from F_(p^2)^x to F_p^x is onto, so f is orthogonal to chi. Let h=f-mu*1. Then ||h||^2=|Gamma\K|(1-mu), and h is in the complementary subspace. Right K invariance and right G equivariance give
T_N f(Gamma)=<f,T_N f>/|Gamma\K|.
Also b_p(N)=4 T_N f(Gamma) and r(N)=4c_e sigma(m). Therefore
b_p(N)/r(N) <= mu + [tau(m)sqrt(m)/sigma(m)](1-mu).
This proves (35), with no missing factor of 2 or 4.

When p divides N, anisotropy of x^2+2y^2 modulo p forces both blocks to vanish if either does. Thus there are no spoiled points when p^2 does not divide N. Otherwise their exact ratio is sigma(m/p^2)/sigma(m)<p^-2<=mu.

### 6. Sieve and numerical cutoff

On the dyadically primitive shell both blocks are nonzero. Counting the p-divisible block first gives at most 4N/p^2 choices after division by p; each other positive norm has at most 2M(N) coordinate witnesses. With both orientations this is 16NM(N)/p^2. Since r(N)=4c_e sigma(m)>=4c_e m, (36) follows with L_e=4*2^e/c_e<=16/3.

The elementary progression estimates and the rational bound for B in (37) have the claimed direction. Combining the small-prime estimate with the dyadic shell lower bound 1/2-E(m)/2 gives (38): the coefficient P/4+3/2 consists of P/4+1 from the prime count and 1/2 from the dyadic error.

The final arithmetic in the global cutoff also checks. The bounds M(N)<=9m^(1/6), tau(m)<=42m^(1/6), P=m^(1/4) give
B+(45/2)m^(-1/12)+159m^(-1/3).
For m>=2^75, the manuscript's rational upper bound is below 1/2 with margin
2620327/234881024 > 1/100.
N>=2^78 and 16 not dividing N imply m>=2^75. Thus the announced onset follows from the stated divisor-factor certificates and the spectral proof above. It is not a proof that the residual family below the onset has been exhausted.

The distinct-value bound follows as written from r(N)>=3N and multiplicity <=144N^(1/3), giving N^(2/3)/4800.


## Exact classification and novelty comparison

### A proved elementary refinement of the residual family

Write C for the odd positive integers supported on primes congruent to 1 or 3 modulo 8. Unique factorization in Z[sqrt(-2)] gives P=C union 2C. This identity is valid for values, not for arbitrary products of fixed primitive coordinate witnesses. Its proof uses norm Euclideanity, the split/inert criterion, and the ramified exponent bound v2(P)<=1. The rounding proof in the existing exact-classification note supplies Euclideanity without an additional source.

A least admissible counterexample N to P+P={N>=2:16 does not divide N} cannot satisfy v2(N)=2. Indeed N=4m with odd m. The smaller eligible target 2m is represented by leastness. A primitive norm is either odd or 2 modulo 4, so a pair summing to 2 modulo 4 must consist of two odd values a,b in C. Then N=2a+2b, and 2a,2b belong to 2C subset P. This contradicts nonrepresentation. This argument is elementary and does not use the spectral input.

Consequently, conditional on the previously claimed all-prime theorem and finite core cutoff, the exact classification is equivalent to coverage of the smaller residual family

    N=2^e m>=2, e in {0,1,3}, m<2^60 odd,
    every prime divisor of m congruent to 5 or 7 modulo 8.

Equivalence here is by least counterexample, not a direct representation-preserving division for every target. The e=2 family can be removed from a least-counterexample search; it cannot be called exhausted by computation.

### The precise remaining parity subproblems

Every member of C is 1 or 3 modulo 8. Thus C+C is contained in the even residue classes 2,4,6 modulo 8, and C+2C consists of odd integers. The exact identities proposed by the full conjecture reduce to the following obligations:

* Odd N: show N belongs to C+2C.
* v2(N)=1: show N belongs to C+C.
* v2(N)=2: the existence problem follows by doubling a C+C representation of N/2, whose valuation is 1.
* v2(N)=3: both primitive summands must be even, since two elements of C cannot sum to 0 modulo 8. Therefore N belongs to P+P if and only if N/2 belongs to C+C. The latter target is 4 modulo 8.

In particular, the e=3 residual problem is a C+C problem at targets 4m; ordinary doubling of a representation at valuation 2 does not settle it. There was no proof of exact C+C coverage in the reviewed revision. The contemporaneous completion is audited separately below.

### Why the natural multiplicative descent stops

Multiplication by an odd split prime preserves P as a set of values because C and 2C are preserved. Multiplication by an inert prime does not preserve P: already 1 belongs to P while 5 does not. Multiplication by 2 is only safe from C into 2C; its next iteration produces a value divisible by 4, outside P. Thus the split-prime descent in the manuscript is justified, but it supplies no inert-prime descent or general dyadic descent.

Composition identities on Q or its quaternionic sum do not automatically preserve separate primitivity. For example multiplying a primitive binary value by an inert prime fails at the level of values, so no identity claiming such unconditional closure can work. Quaternionic rotations can create new inert factors in individual summand values even when they preserve the total norm. Any proposed stronger descent must show how it avoids those factors in both summands.

The existing unit-coordinate ternary reduction is correctly limited. The exact counterexample 29=x^2+u^2+2v^2 only has v=0 and {|x|,|u|}={2,5}, so no representation has gcd(u,v)=1. Ordinary ternary representability therefore cannot justify the required primitive-block statement. More strongly, the existing elementary equivalence

    t^2 represented by x^2+u^2+2v^2 with gcd(u,v)=1
    iff 2t belongs to C+C

shows that exact primitive ternary regularity on odd squares already contains the valuation-1 part of the original binary conjecture. This is a precise circularity obstruction for that proposed route, not a proof that every ternary approach must fail.

### Novelty comparison against primary sources

Jorg Brudern, Sums of Squares Coprime in Pairs, Bull. Polish Acad. Sci. Math. 62 (2014), 215-231, DOI 10.4064/ba62-3-3, Theorem 2 (printed p.219), proves a power-saving asymptotic for four positive squares with all six pairwise gcd conditions. Theorem 1 (p.215) gives a uniform linear lower bound on the six admissible classes modulo 24. Lemma 2 gives the large-common-factor tail; equations (3)-(8) supply truncated Mobius inversion; Lemma 4 imports a uniform weighted-quaternary estimate; Lemma 5 completes an absolutely convergent weighted series; Section 5 recovers its local-density product.

Primary PDF inspected from the repository and independently located online:
https://www.impan.pl/shop/en/publication/transaction/download/product/85773
Primary journal landing page:
https://www.impan.pl/en/publishing-house/journals-and-series/bulletin-polish-acad-sci-math/all/62/3/85773/sums-of-squares-coprime-in-pairs

This is a very close methodological antecedent. Its coprimality graph is stronger than the two disjoint block conditions here, but its quadratic form is the particular four-square form. It does not state the arbitrary fixed positive-definite binary-pair theorem, arbitrary ideal-class scope, or nonmaximal-order specialization. Conversely its theorem cannot imply the proposed arbitrary-form theorem merely because it imposes more gcd conditions: the quadratic forms and admissible sets differ. Even for four squares its local congruence classes are a proper subset of those for only two disjoint gcd conditions.

Jeremy Rouse, Integers represented by positive-definite quadratic forms and Petersson inner products, arXiv:1802.03437v1, Theorem 1, states representation criteria for a single primitive quaternary form under successively weaker local conditions. Its primitive local condition concerns the entire four-vector. It does not impose separate primitive conditions on two fixed binary blocks. Therefore the manuscript must import the uniform all-coefficient cusp estimate and supply the two-block sieve; it cannot cite Rouse Theorem 1 alone as the claimed theorem.

Primary source:
https://arxiv.org/pdf/1802.03437v1

No exact antecedent for the arbitrary-binary-pair, separately primitive theorem was identified in the targeted searches. This is not a priority determination. The defensible scope statement is that the proof architecture and its analytic ingredients are classical; the exact arbitrary-pair deduction may be useful, with novelty unestablished. The explicit -8 spectral constants, onset, and prime-family proof are separate specialized assertions requiring their own spectral audit. The reviewed revision left the full -8 identity open; the completion below resolves it using additional input. Neither unweighted distinct-value asymptotics nor simultaneous membership in two norm semigroups follows from the coordinate asymptotic alone.



## Verification record and remaining obligations

The existing standard-library quaternary checker was replayed at integer
bound 256. It passed 3282 good-prime residue checks, 68 lifting-residue
checks, 1792 binary Möbius identities, 1020 double-Möbius targets, the
four-square normalization, and the dyadic Sturm/divisor certificates.
The all-prime certificate covered all 295947 primes through \(2^{22}\),
with the same retained membership and sumset hashes and no fallback
targets. This checks the finite prime family; it does not exhaust
the composite inert-core residual family.

A separate finite matrix implementation challenged reduction primes
5 and 7, odd indices 1,3,9,11, and dyadic factors 1,2,4,8.
The quotient sizes were 60 and 168. Tested matrices have integral
entries exactly representable in binary64. Symmetry, both coprime
factorizations and normalizer fractions held exactly in these finite
calculations; complementary spectra were checked numerically with
NumPy 1.26.4 and the stated tolerance. These are normalization
challenges, not a uniform theorem. The script and receipt remain local
support under .scratch/primitive-norm-audit; no theorem depends on
this exploratory calculation.

Rouse PDF p.16 and Schulze-Pillot--Yenirce PDF p.12 were rendered and
visually inspected for exponents and normalizations. Exact primary
editions and hashes are retained in
[the source ledger](primitive_norm_quaternary_source_ledger.json).
The proofs of Rouse's Petersson norm estimate, global
Jacquet--Langlands and Ramanujan remain classical imports.

The exact -8 classification was open in the reviewed revision and is
resolved by the addendum below. A representation-preserving inert-prime
descent, unweighted distinct-value leading asymptotics and a complete
priority determination remain open. The \(C_{12}+C_{12}\) problem
requires two norm conditions per summand and is outside the theorem.
No external-review, formal-verification, or independent-replay flag
is upgraded by this internal audit.


## Addendum: audit of the contemporaneous exact -8 completion

The separate [exact proof](PRIMITIVE_NORM_MINUS8_EXACT.md), its
[checker](verify_primitive_norm_minus8.py) and
[receipt](primitive_norm_minus8_checks.json) appeared in the shared
working tree after the audit of revision 841f6e73699b had begun.
They were not produced by the audit's classification subtask.
Their new claims were assessed separately before accepting the
current classification.

The completion proves the stronger identities
\[
 C+C=\{N\ge2:N\text{ even},8\nmid N\},\qquad
 C+2C=\{N\ge3:N\text{ odd}\}.
\]
Together with \(P_{-8}=C\cup2C\), these imply exactly
\[
 \boxed{P_{-8}+P_{-8}=\{N\ge2:16\nmid N\}}.
\]

The proof is not an exhaustion of all inert cores below \(2^{60}\).
It uses three selected shells, at \(m,2m,4m\), with respective
surviving proportions \(1/2-E(m)/2,1/2,2/3\). The latter two
proportions are exact odd-odd shell counts
\(r(2m)-r(m)\) and \(r(4m)-r(2m)\). Multiplication or division
by \(\sqrt{-2}\) justifies the dyadic maps on coordinates and preserves
all required odd-prime divisibility conditions.

The sharper lattice bounds
\(S_Q(X)\le(8/3)X\), \(S_o(X)\le2X\), the inert-core estimate
\(E(m)<(18/5)m^{-1/3}\), and the exact split-divisor frontier
give 107 rational interval inequalities. The checked interval
indices are 27--59 for odd cores, 26--59 for \(2m\), and 20--59
for \(4m\). Every inequality has margin greater than \(1/256\).
The independently audited original large-core theorem supplies all
larger cores. Its use at dyadic exponents 0,1,3 gives precisely these
three strong-core conclusions by parity. No exact-classification
conclusion is fed back as a premise.

For the remaining small range, the complete finite certificate covers
both strong identities through \(2^{27}\): 50,331,648 even strong
targets and 67,108,863 mixed odd targets. It sieves every forbidden
prime, constructs exact core membership flags, and verifies complete
periodic target masks using shifts by actual positive short summands.
The corruption checks detect deleted eligible targets and inserted
forbidden positions. A short-summand assumption is unnecessary because
the computed union actually equals each complete required mask.

This audit inspected the full checker, separately replayed the sharp
bounds and all interval certificates, and replayed the **entire**
finite computation from a byte-identical checker snapshot with its
output directed outside the retained receipt. The resulting JSON
equals the retained receipt exactly, including all coverage hashes.
The checker SHA-256 is
0a7cdfbd1642e54d12c6882c119b1de6168bf438a70b9b15467f343fdfecab96.

Finally, a least missing eligible strong-core target with a split
prime factor descends to a smaller eligible target by division by
that prime. The odd mixed problem has the additional prime boundary,
already supplied by the independently audited all-prime theorem.
The finite certificate, interval families and large-core theorem
therefore cover every remaining least counterexample. This closes
both strong identities and the full primitive identity.

Verdict: the exact classification is now supported by a complete
internal mathematical/computational proof under the explicitly
identified classical spectral imports. There is no remaining
unexhausted finite obligation for this identity. This addendum is an
internal audit and replay, not external human review, formal
verification, or an established novelty claim.


The subsequent [novelty comparison](PRIMITIVE_NORM_NOVELTY.md) applies newer classical inputs and supersedes the earlier quantitative choices with eta<1/8 and the general distinct-value lower bound N/(log N)^(log4+epsilon). The Rouse-based derivation above remains a valid independent weaker route.
