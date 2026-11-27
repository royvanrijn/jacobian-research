# A critical rectangular six-cycle criterion for the five-seventeenths bound

For positive integers kappa_1,kappa_2 let R_kappa(n) count the ordered coprime positive representations
\[
 \kappa_1a^2b^3+\kappa_2c^2d^3=n,\qquad b,d\ \hbox{squarefree}.
                                                               \tag{1}
\]
Coprimality means gcd(kappa_1*a^2*b^3,kappa_2*c^2*d^3)=1. The coefficients are fixed during each count; write kappa_max=max(kappa_1,kappa_2) and kappa_min=min(kappa_1,kappa_2). The prescribed-sum count r(N) allows noncoprime positive square-full summands.

Put
\[
 t_0=\frac5{17},\qquad
 A_n=2n^{11/51},\quad B_n=2n^{32/153},\quad
 \Gamma=\frac{26}{51}.                                        \tag{2}
\]
Use the coefficient-shaped upper cutoffs
\[
\begin{aligned}
 A_{1,n}&=A_n/(\kappa_1^2\kappa_2)^{1/3},&
 A_{2,n}&=A_n/(\kappa_1\kappa_2^2)^{1/3},\\
 B_{1,n}&=B_n/(\kappa_1^4\kappa_2^2)^{1/9},&
 B_{2,n}&=B_n/(\kappa_1^2\kappa_2^4)^{1/9}.
\end{aligned}                                                 \tag{2s}
\]
The rectangular graph G_kappa(n) has kernel vertices b<=B_(1,n), d<=B_(2,n), and an edge when (1) has a representation with a<=A_(1,n), c<=A_(2,n). Keep one selected actual representation per edge. All cycles below are simple six-cycles. Statements uniform in the coefficients are for each fixed coefficient pair, not an aggregate over varying pairs. Since kappa_i>=1, the shaped cutoffs are bounded by A_n,B_n. A cutoff below one makes the graph empty.

Consider the following still OPEN hypothesis:
\[
\begin{gathered}
 \#C_6(G_\kappa(n))\ll_\varepsilon n^{\Gamma+\varepsilon}
 \quad\hbox{for every }n\ge2,\\
 \kappa_1,\kappa_2\ \hbox{positive, squarefree and coprime},\qquad
                         \kappa_{\max}^{134}\kappa_{\min}^{117}\le n.
\end{gathered}                                                 \tag{H}
\]
A fixed multiplicative change in the two cutoff constants in (2) is allowed if used consistently.

**Theorem 1 (CriticalRectangularCriterion, conditional).** Hypothesis (H) implies
\[
                         r(N)\ll_\varepsilon N^{5/17+\varepsilon}
                         \quad\hbox{for every }N\ge2.           \tag{3}
\]
No unconditional improvement on the current three-tenths exponent is asserted. The hypothesis is weaker than a full T^(2+epsilon) cycle estimate in the larger maximum box T=n^(11/50): the cycle exponent allowed in (H) is 26/51, rather than 22/50.

We prove both the reduction and exact counted sectors that can be removed from (H).

## 1. The five existing box estimates and their critical locus

For a nonempty dyadic box write alpha,beta,gamma,zeta for the logarithms to base n of the LOWER endpoints of a,b,c,d. Set X=alpha+gamma and Y=beta+zeta. The five coefficient-uniform estimates in [the dual-graph argument](ERDOS_943_DUAL_KERNEL_GRAPHS.md) give
\[
 R_{\kappa,\mathrm{box}}(n)\ll_\varepsilon n^{F+\varepsilon},
\]
where
\[
\begin{split}
 F=\min\{&
 Y,\ X,\\
 &\max(\beta,\zeta,
       \min(\beta+\zeta/2,\zeta+\beta/2)+(X-1/2)_+),\\
 &\max(\alpha,\gamma,
       \min(\alpha+\gamma/2,\gamma+\alpha/2)+(Y-1/2)_+),\\
 &\min(\alpha+\zeta+\max(\beta,\gamma)/3,\,
       \beta+\gamma+\max(\alpha,\zeta)/3)\}.
\end{split}                                                    \tag{4}
\]
The base polytope has 2alpha+3beta<=1 and 2gamma+3zeta<=1. More precisely, for k_i=log_n(kappa_i), nonemptiness gives
\[
 L:=2\alpha+3\beta\le1-k_1,\qquad
 R:=2\gamma+3\zeta\le1-k_2.                                   \tag{5}
\]
The letter R in (5) is a scalar resource, not a counting function.

**Lemma 2 (CriticalLocalizationAndCoefficientGate).** If F>=t>5/17, then
\[
\begin{gathered}
 \alpha\le1-\frac83t-\frac{2k_1+k_2}{3},\qquad
 \gamma\le1-\frac83t-\frac{k_1+2k_2}{3},\\
 \beta\le\frac23-\frac{14}{9}t-\frac{4k_1+2k_2}{9},\qquad
 \zeta\le\frac23-\frac{14}{9}t-\frac{2k_1+4k_2}{9},\\
 F\le\frac3{10}-\frac{2k_{\max}+k_{\min}}{10}.
\end{gathered}                                                 \tag{6}
\]
Consequently F>5/17 forces
\[
 \alpha,\gamma<\frac{11}{51},\qquad
 \beta,\zeta<\frac{32}{153},\qquad
 \kappa_{\max}^2\kappa_{\min}<n^{1/17}.                         \tag{7}
\]
In particular kappa_i<n^(1/34).

Here is a complete case argument. Fixed-pair bounds give X,Y>=t, so 2X+3Y<=2 implies
\[
 X\le1-\frac32t,\qquad Y\le(2-2t)/3<1/2.                       \tag{8}
\]
Thus the dual graph has no height loss. Suppose beta>=t. Then alpha<=(1-3t)/2. If gamma<t, its dual-graph square-root term is at most alpha+gamma/2<1/2-t<t, while alpha,gamma<t; this contradicts F>=t. Hence gamma>=t, and zeta<=(1-2t)/3. When gamma>=beta, the first mixed bound is at most
\[
 \alpha+(1-\gamma)/3\le5/6-11t/6.
\]
When beta>=gamma, it is at most
\[
 (1-3\beta)/2+\beta/3+(1-2t)/3\le5/6-11t/6.
\]
Both contradict t>5/17. By symmetry beta,zeta<t.

If X>1/2, the kernel square-root term is at most
\[
 3Y/4+X-1/2\le X/2\le1/2-3t/4<t,
\]
since t>2/7. Together with beta,zeta<t this contradicts F>=t. Thus X<=1/2. If alpha>=t, then beta<=(1-2t)/3, and the kernel square-root term is at most beta+zeta/2<1/3-t/6<t, another contradiction. Similarly gamma<t.

All four coordinates are below t and both graph height losses are zero. Both square-root terms must therefore be at least t:
\[
 2\alpha+\gamma,\ \alpha+2\gamma,\ 
 2\beta+\zeta,\ \beta+2\zeta\ge2t.                             \tag{9}
\]
Using the actual resource upper bounds (5) gives
\[
\begin{aligned}
 2\alpha+\gamma&\le3/2-3t-(2k_1+k_2)/2,\\
 \alpha+2\gamma&\le3/2-3t-(k_1+2k_2)/2,\\
 2\beta+\zeta&\le1-4t/3-(2k_1+k_2)/3,\\
 \beta+2\zeta&\le1-4t/3-(k_1+2k_2)/3.
\end{aligned}                                                 \tag{9s}
\]
Twice one upper bound minus the opposite lower bound in (9), divided by three, gives the coordinate bounds in (6). The lower inequalities also give X,Y>=4t/3 and hence t<=3/10.

Finally, weighting the first and third inequalities of (9) by two and three gives
\[
             2L+R\ge10t,\qquad L+2R\ge10t.
\]
Equation (5) yields 10t<=3-2k_1-k_2 and 10t<=3-k_1-2k_2. Taking t=F proves the coefficient assertion. Substituting t>5/17 gives (7).

**Corollary 3 (CoefficientSensitiveRepresentationBound).** Uniformly in positive coefficients,
\[
 R_\kappa(n)\ll_\varepsilon n^\varepsilon
   \max\left\{n^{5/17},
       \frac{n^{3/10}}{(\kappa_{\max}^2\kappa_{\min})^{1/10}}
       \right\}.                                             \tag{10}
\]
For a nonempty graph kappa_i<=n. Apply Lemma 2 when F>5/17, and the first term otherwise, then sum the logarithmically many boxes. The constants remain independent of coefficient sizes. This is an unconditional coefficient refinement, while the coefficient-one case still has exponent 3/10.

## 2. From the rectangular cycle budget to representations

For fixed b,d, equation (1) has O_epsilon(n^epsilon) positive root pairs uniformly in the coefficients. Indeed put U=kappa_1*b^3, V=kappa_2*d^3. In the imaginary quadratic field Q(sqrt(-UV)), the algebraic integer Ua+c sqrt(-UV) has norm U*n. Prime-ideal factorization bounds the number of ideals of that norm by tau(Un), and a principal ideal has at most six generators of the specified norm. A nonempty pair has U<=n, so the divisor count is subpower in n. This is the fixed-kernel argument in the same dual-graph source.

Morris--Saxton's simple-cycle supersaturation [Theorem 1.5(a), ell=3](https://arxiv.org/pdf/1309.2927v3) gives
\[
                    e(G)\ll v^{4/3}+v\,\#C_6(G)^{1/6}         \tag{11}
\]
for a graph on v vertices. The first term accounts for edge counts below its supersaturation threshold.

In G_kappa(n), v<=2B_n=4n^(32/153). Under (H),
\[
 e(G_\kappa(n))\ll_\varepsilon
       n^{128/459}+n^{32/153+(26/51)/6+\varepsilon}
       \ll_\varepsilon n^{5/17+\varepsilon},                 \tag{12}
\]
because 128/459<135/459=5/17 and (32+13)/153=5/17.

For a coefficient pair allowed by (H), every box with F>5/17 is contained in the coefficient-shaped rectangle by (6): its UPPER root endpoints are at most A_(1,n),A_(2,n), and its UPPER kernel endpoints are at most B_(1,n),B_(2,n). The factors two in (2s) are explicit. There is no silent identification of lower and upper dyadic endpoints. The fixed-kernel multiplicity and (12) bound all these boxes by n^(5/17+epsilon). Other boxes already have that bound by (4). Summing boxes proves
\[
 R_\kappa(n)\ll_\varepsilon n^{5/17+\varepsilon}
 \quad\hbox{when kappa_i are squarefree, coprime and }
                     \kappa_{\max}^{134}\kappa_{\min}^{117}\le n,             \tag{13}
\]
under (H).

## 3. Square-full gcd normalization and the exact coefficient gate

Take an original ordered positive square-full pair m_1+m_2=N, and let h=gcd(m_1,m_2), n=N/h. The gcd h is itself square-full: every prime dividing both original summands has valuation at least two in each, so its minimum valuation is at least two. In particular
\[
                         \operatorname{rad}(h)^2\mid h.        \tag{14}
\]
The positive quotients m_i/h are coprime. Let kappa_i be the product of primes occurring in m_i/h with exponent exactly one. Then kappa_1,kappa_2 are squarefree and coprime; after removing them each quotient is square-full. Thus the original pair injects into a representation (1) for its uniquely specified h,kappa_i,a,b,c,d.

Each exponent-one prime divides h, since m_i itself is square-full. Put K=kappa_1*kappa_2. Then
\[
                  K\mid\operatorname{rad}(h),\qquad K^2\mid h. \tag{15}
\]
Not every arbitrary normalized tuple must yield two square-full original summands; only this injection is used for an upper bound.

Write Q=kappa_max^2*kappa_min and W=K^100*Q^17=kappa_max^134*kappa_min^117. If W>=n, then h^50*Q^17>=n by (15). Consequently
\[
 \left(\frac{n^{3/10}}{Q^{1/10}}\right)^{170}
       =\frac{n^{51}}{Q^{17}}
       \le n^{50}h^{50}=N^{50}.
\]
Corollary 3 therefore bounds every such normalization by
\[
 R_\kappa(n)\ll_\varepsilon n^\varepsilon
     \max(n^{5/17},N^{5/17})\le N^{5/17+\varepsilon}.            \tag{16}
\]
If W<n, equation (13) applies under (H), giving the same bound because n<=N. There are at most tau(N)^3 coefficient/gcd normalization choices; the actual square-full and divisor restrictions only reduce this number. Divisor and dyadic factors are absorbed by choosing smaller intermediate epsilons. This proves Theorem 1 for all prescribed targets, including noncoprime summands.

For perspective, the coefficient family in (H) is small:
\[
  W=K^{251/2}(\kappa_{\max}/\kappa_{\min})^{17/2}\le n
  \quad\Longrightarrow\quad
     K\le n^{2/251},\qquad \kappa_{\max}\le n^{1/134}.            \tag{17}
\]
The exact weighted gate is stronger than the product-only condition. This restriction comes from an unconditional coefficient-sensitive estimate and the square-full gcd, not from discarding normalizations.

Theorem 1 is a sufficient cycle criterion. No converse from the representation bound to (H) is asserted. The next sections give an actual equivalence between (H) and a more restricted cycle count.

## 4. Actual rectangular codegrees and opposite-root-gcd frequency

All counts below are uniform in positive coefficients, whether or not they satisfy the coefficient condition in (H). Restricting them to (H)'s coefficient family is legitimate. For large n, 16A_n^4<n. The finitely many smaller targets can be absorbed uniformly, since a nonempty equation has kappa_i<=n.

**Lemma 4 (RectangularCodegreeRecovery).** In the rectangle (2), either pair of distinct same-side kernels has O_epsilon(n^epsilon) actual common neighbors.

For distinct b_1,b_2 sharing d, elimination gives
\[
 b_1^3(a_1c_2)^2\equiv b_2^3(a_2c_1)^2\pmod n.
\]
All factors and coefficients are units modulo n by summand coprimality. The reduced positive fraction u/v=(a_1c_2)/(a_2c_1) has u,v<=A_n^2 and lies in one of at most 2^(omega(n)+1) unit square-root classes. Two distinct such fractions in one class have cross difference a nonzero multiple of n, but its absolute value is at most 2A_n^4<n. Hence each class has at most one fraction.

Write a_i=f alpha_i, c_i=g gamma_i, removing the gcd of each pair. Pairwise and actual within-edge coprimalities make u=alpha_1*gamma_2 and v=alpha_2*gamma_1 already reduced. There are divisor-many choices for these four factors. The exact identity is
\[
 \kappa_1 f^2(b_1^3u^2-b_2^3v^2)=n(\gamma_2^2-\gamma_1^2).
\]
Its coefficient is nonzero: equality would force b_1=b_2 by squarefree prime-valuation parity. It determines f at most once. The original edge then determines
\[
 g^2d^3=(n-\kappa_1 b_1^3f^2\alpha_1^2)/(\kappa_2\gamma_1^2).
\]
A positive integer has at most one squarefree-cube-times-square decomposition, so this determines g,d. Nonintegral candidates are discarded. All factorized integers are polynomially bounded in n. Swapping summands gives the other side.

**Theorem 5 (RectangularLargeOppositeRootGcdSector).** Cycles containing any shared-neighbor path whose two opposite roots have gcd at least n^(43/306) number O_epsilon(n^(26/51+epsilon)).

At a shared right neighbor let a_i=alpha*r_i with alpha=gcd(a_1,a_2)>=n^(43/306). Then r_i<=2n^(23/306). Actual subtraction and gcd(alpha,kappa_2*d)=1 imply
\[
 \ell=(c_1^2-c_2^2)/\alpha^2\in\mathbb Z\setminus\{0\},
 \qquad |\ell|\le4n^{46/306}.
\]
The difference cannot be zero, since that would force the distinct squarefree b_i to coincide. The primitive kernel triple satisfies
\[
       \kappa_1 r_1^2b_1^3-\kappa_1 r_2^2b_2^3
                                      +\kappa_2\ell d^3=0.      \tag{18}
\]
There are O(n^(46/153)) triples (r_1,r_2,ell). Heath-Brown's [coefficient-uniform primitive cubic theorem, Theorem 1.3](https://arxiv.org/abs/2601.07817v1), bounds each triple's kernel completions by n^epsilon. The third kernel is coprime to each of the first two, ensuring primitivity. A fixed kernel triple has at most one selected actual path.

Every specified two-edge path extends to at most O_epsilon(B_n*n^epsilon) simple cycles: choose the third kernel on that side, then use the two actual codegrees. The exponent is 46/153+32/153=26/51. Sum over the six path positions, and interchange summands for the other orientation. This counts actual frequency, not just completions for a fixed set of labels.

## 5. Low kernel lcm and small coherent gcd

**Lemma 6 (RectangularLowLcmSector).** Cycles with either kernel-side lcm at most n^(26/51) number O_epsilon(n^(26/51+epsilon)).
For each lcm l, all three kernels are divisors of l, so at most tau(l)^3 triples occur. Each distinct triple has subpower many lifts by Lemma 4. Sum over l<=n^(26/51), and use the same argument on the other side.

A cycle is singular when the products of its selected a/c ratios along its two alternating three-edge matchings agree. Equivalently it admits positive rational vertex labels p_i,q_j with p_i/q_j=a_ij/c_ij on every edge, obtained by propagation around the cycle. Normalize these labels to positive integers of global gcd one. Connectedness makes them unique; all are units modulo n. A five-edge spanning tree produces integral labels of height at most A_n^5 by multiplying its five root numerators/denominators, followed by primitive normalization. Let
\[
                       g_*=\min_{\mathrm{edge}}\gcd(p_i,q_j).
\]

**Lemma 7 (RectangularSmallCoherentGcdSector).** Singular cycles with
\[
                              g_*\le n^{107/510}               \tag{19}
\]
number O_epsilon(n^(26/51+epsilon)).

For completeness, fix one selected actual edge and prescribe its intrinsic edge gcd s. It fixes its two coherent labels p_1=s*a_11,q_1=s*c_11. Put L_i=kappa_1*b_i^3*p_i^2 and M_j=kappa_2*d_j^3*q_j^2. Each edge says L_i+M_j=n*g_ij^2. Choose successively
\[
             g_{21}\mid q_1,\quad g_{22}\mid p_2,\quad
             g_{32}\mid q_2,\quad g_{33}\mid p_3.
\]
Each choice forces the next positive potential. Divide by its known coefficient and apply the unique squarefree-cube-times-square decomposition to recover the next kernel and coherent label. Reject nonpositive or nonintegral values and prune every label to A_n^5 BEFORE the next divisor choice. The final edge and primitive normalization are checked. All dividends are polynomially bounded in n, so this gives n^epsilon branches for a fixed edge and s. The same actual reconstruction is used in [the balanced divisor argument](ERDOS_943_SIX_CYCLES.md); the explicit height bound here supplies the rectangular extension without a lower root bound.

Orient a cycle at an edge attaining g_*. Summing over s<=Y gives O_epsilon(n^epsilon*e(G)*Y). The known coefficient-uniform three-tenths theorem bounds e(G)<=R_kappa(n)<<n^(3/10+epsilon). Taking Y=n^(107/510) works because 3/10+107/510=26/51.

## 6. The coefficient-uniform low weighted-ratio sector

For a singular cycle put P=prod_i p_i, Q=prod_j q_j, and
\[
                           \rho=\kappa_1^3P^2/(\kappa_2^3Q^2).
\]
Its height is that of its reduced positive numerator/denominator.

**Corollary 8 (RectangularSmallWeightedRatioSector).** Singular cycles with
\[
                               H(\rho)\le n^{11/68}            \tag{20}
\]
number O_epsilon(n^(26/51+epsilon)), uniformly in the coefficients.

Use [the weighted product-ratio theorem, Theorem 8](ERDOS_943_COHERENT_PRODUCT_RATIO.md). It needs only a kernel cutoff K and gives
\[
 n^\varepsilon\left(K^3X^{4/3}/n^{1/3}+K^6X/n\right).
\]
At K=B_n the two base exponents are 5/17 and 13/51. With X=n^(11/68), the first exponent is 26/51 and the second is 5/12<26/51. The factors two in B_n affect only constants. This supplies an actual cycle-frequency bound with no lost coefficient factor.

### An unweighted ratio with no coefficient loss

**Corollary 9 (ShapedUnweightedProductRatioSector).** In the coefficient-shaped graph, singular cycles with
\[
                         H(P/Q)\le n^{11/136}                  \tag{20s}
\]
number O_epsilon(n^(26/51+epsilon)), uniformly in the coefficients.

Here is the anisotropic counting argument, rather than a removal of a coefficient factor from an unshaped estimate. Set K_1=B_(1,n), K_2=B_(2,n). If a cutoff is below one there are no cycles. Write the reduced positive ratio P/Q=r/s with r,s<=Y, and put B=prod b_i, D=prod d_j. Multiplication along an actual matching gives
\[
                     n\mid\kappa_1^3r^2B^3+\kappa_2^3s^2D^3.
\]
All these factors are units modulo n. Let g=gcd(B,D), u=B/g, v=D/g. Then u,v are positive coprime integers with u<=K_1^3/g,v<=K_2^3/g. Positivity gives
\[
 g\le\frac{Y^{2/3}}{n^{1/3}}
        (\kappa_1^3K_1^9+\kappa_2^3K_2^9)^{1/3}
 \le\frac{Y^{2/3}}{n^{1/3}}
                         (\kappa_1K_1^3+\kappa_2K_2^3).
\]
The ratio u/v is in a subpower number of unit cube-root classes modulo n. In one such class, positive primitive fractions in a rectangular range u<=H_1,v<=H_2 number O(1+H_1H_2/n): rescale the vectors to (u/H_1,v/H_2), whose determinants are nonzero multiples of n/(H_1H_2), then use the same angular packing proof as for a square. Apply this at H_i=K_i^3/g. Summing g up to its positive cutoff gives, for each r,s,
\[
 \ll_\varepsilon n^\varepsilon\left[
 \frac{(\kappa_1K_1^3+\kappa_2K_2^3)Y^{2/3}}{n^{1/3}}
                        +\frac{(K_1K_2)^3}{n}\right].
\]
There are at most Y^2 centers. Each product pair B,D has divisor-many kernel factorizations because actual edges give kappa_1*B<n,kappa_2*D<n. A fixed ordered kernel sextuple has at most one cycle in the selected-edge graph. Thus the cycle frequency is
\[
 \ll_\varepsilon n^\varepsilon\left[
 \frac{(\kappa_1K_1^3+\kappa_2K_2^3)Y^{8/3}}{n^{1/3}}
                        +\frac{(K_1K_2)^3Y^2}{n}\right].       \tag{20u}
\]
This general count uses only the two kernel caps, not any root cap or codegree estimate.

The shaped caps give
\[
 \kappa_1K_1^3\le8n^{32/51},\quad
 \kappa_2K_2^3\le8n^{32/51},\quad
 (K_1K_2)^3=\frac{64n^{64/51}}{(\kappa_1\kappa_2)^2}
                                                        \le64n^{64/51}.
\]
The first inequality, for example, retains the factor
kappa_1^(-1/3)*kappa_2^(-2/3)<=1. Consequently (20u) is at most
\[
          n^\varepsilon
           (n^{5/17}Y^{8/3}+n^{13/51}Y^2)
\]
up to an absolute constant. Taking Y=n^(11/136) makes the first exponent 26/51 and the second 5/12, proving the corollary. The cancellation of the coefficient loss comes from the coefficient-shaped kernel caps. It does not give the same claim in an arbitrary unshaped maximum box.

## 7. Both directions of the remaining rectangular cycle problem

In the graph (2), let E_kappa(n) consist of cycles satisfying all the following:
each kernel-side lcm exceeds n^(26/51); every opposite-root gcd is less than n^(43/306); and, if the cycle is singular, g_*>n^(107/510), H(rho)>n^(11/68), and H(P/Q)>n^(11/136).

**Proposition 10 (CriticalRectangularResidualEquivalence).** Hypothesis (H) is equivalent to
\[
             \#E_\kappa(n)\ll_\varepsilon n^{26/51+\varepsilon}
                                                               \tag{21}
\]
for exactly the same targets and squarefree coefficient family in (H).

The forward direction is subset inclusion. Conversely the cycles outside E_kappa(n) lie in the union of the low-lcm sector, the opposite-root-gcd sector, or the three singular sectors (19), (20), (20s). Lemma 6, Theorem 5, Lemma 7 and Corollaries 8--9 bound each by n^(26/51+epsilon); the union bound proves (H). Both cycle estimates remain OPEN, and hence (3) remains conditional on them.

This equivalence concerns the specified rectangular graph and exponent 26/51. It is not an equivalence with the older balanced T^(2+epsilon) obligation, or with the representation bound (3). No prime-only sign choice, proper-power removal, or coefficient-dependent estimate is used in these two directions.

## 8. Classical inputs and the remaining gap

The five box estimates, uniform three-tenths theorem, fixed-kernel multiplicity and gcd normalization are in the dual-graph source linked above. The cycle conversion uses Morris--Saxton's precise simple six-cycle supersaturation statement. The opposite-root-gcd frequency uses Heath-Brown's primitive diagonal-cubic theorem uniformly in nonzero signed integer coefficients. Primitivity follows from actual coprime summands; coefficients are fixed during each graph count.

The new localization, coefficient gate, large-gcd split and rectangular-sector arithmetic are written in full here. The two cycle estimates in Section 7 remain mathematical obligations. No finite experiment, count for a fixed label tuple, or imported average over targets is substituted for their pointwise frequency bound.

The [rectangular continuation argument](ERDOS_943_CRITICAL_RECTANGULAR_CONTINUATION.md) proves further sectors at the same pointwise budget and gives a narrower equivalent cycle count within this coefficient family.
