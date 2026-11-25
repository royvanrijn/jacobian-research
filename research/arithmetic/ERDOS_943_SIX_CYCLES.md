# Six-cycles in the square-full kernel graph

Fix a coprime representation box \(a,b,c,d\asymp T\) for
\[
                          a^2b^3+c^2d^3=N,\qquad N\asymp T^5 ,
\]
with \(b,d\) squarefree. The simple graph \(G\) has kernel vertices \(b\) on the left and \(d\) on the right. Choose one representation for each edge. The fixed-kernel and common-neighbor bounds for actual representations are \(O_\varepsilon(N^\varepsilon)\) on this scale. [Kernel graph and recovery](ERDOS_943_KERNEL_GRAPH.md)

**Proposition.** A pointwise \(O_\varepsilon(T^{2+\varepsilon})\) bound for the number of simple six-cycles would imply \(O_\varepsilon(N^{4/15+\varepsilon})\) representations in the balanced box. On average in \(N\) on a dyadic interval, with \(T=X^{1/5}\) fixed, the six-cycle count is \(O_\varepsilon(T^{1/2+\varepsilon})\).

For the first assertion, take the graph to have \(O(T)\) vertices. Above a sufficiently large constant times the \(T^{4/3}\) edge threshold, even-cycle supersaturation gives \(\gg (e(G)/T)^6\) simple six-cycles. This is Morris–Saxton's Theorem 1.5(a) with \(\ell=3\), applied with any needed isolated vertices added so that the vertex parameter is comparable to \(T\). Thus a \(T^{2+\varepsilon}\) bound gives \(e(G)\ll_\varepsilon T^{4/3+\varepsilon}\) after changing \(\varepsilon\). Each edge has at most subpower multiplicity and \(T\asymp N^{1/5}\). [Morris–Saxton, Theorem 1.5](https://arxiv.org/pdf/1309.2927v3)

For the mean, every simple six-cycle projects to a triangle on its three left kernels. Each pair of this triangle must be compatible in the sense of \(K_N(T)\), with the appropriate fixed constants in its kernel intervals and fraction height. The three choices of shared neighbors lift a compatibility triangle to at most \(O_\varepsilon(N^\varepsilon)\) cycles by common-neighbor recovery. A graph with \(K\) ordered adjacencies has \(O(K^{3/2})\) triangles, for instance by bounding the trace of the cube of its adjacency matrix by the \(3/2\)-power of the trace of its square. Therefore
\[
                       \#C_6(G)\ll_\varepsilon N^\varepsilon K_N(T)^{3/2}.
\]
Since \(K_N(T)\ll T^2\) pointwise and has mean \(\ll_\varepsilon T^{1+\varepsilon}\), this already gives \(O_\varepsilon(T^{2+\varepsilon})\) on average. A uniform \(K_N(T)\ll_\varepsilon T^{4/3+\varepsilon}\) would suffice for the desired pointwise cycle bound by the same inequality. [Kernel-pair mean](POWERFUL_SUM_GRAPHS.md)

For the stronger mean, let \(L\asymp T\) bound both side sizes, and let \(K\ge1\) bound the two codegrees of \(G\). We show
\[
                           \#C_6(G)\ll K^{5/2} L^{3/2}e(G).
\]
Call a right vertex low if its degree is at most \(D=4\sqrt{K L}\).
To count cycles through a low right vertex \(d\), choose its two left neighbors \(b,b'\). The number of paths of length two from \(b\) to another left vertex is at most \(K L\), by the left codegree bound. Once that path ends at \(b''\), the final right vertex joining \(b''\) and \(b'\) has at most \(K\) choices. The total contribution is at most \(K^2 L\sum_{d\ \mathrm{low}}\deg(d)^2\le K^2 LD\,e(G)\).

If there are \(h\) high right vertices and \(e_h\) edges incident to them, the right codegree bound gives \(e_h^2\le L(e_h+K h^2)\) by Cauchy–Schwarz on the left degrees. Thus \(hD\le e_h\le L+h\sqrt{K L}\), so \(h\ll\sqrt{L/K}\). To count cycles with all three right vertices high, choose one edge and the other two right vertices. The other two left vertices are common neighbors of specified pairs, so there are at most \(O(e(G)h^2 K^2)\) such cycles. These two estimates prove the displayed inequality.

Across all \(N\) in the dyadic interval, each quadruple of labels \(a,b,c,d\) defines one target \(N\). There are \(O(T^4)\) such quadruples, hence \(\sum_N e(G)\ll T^4\). Since \(K\ll_\varepsilon T^\varepsilon\), the graph inequality gives total \(O_\varepsilon(T^{11/2+\varepsilon})\) six-cycles. Dividing by \(X=T^5\) proves the stronger mean.

## Exact equations for cycle labels

**Proposition.** For a cyclically ordered six-cycle let its edges be \((b_i,d_i)\), \((b_{i+1},d_i)\) for \(i\in\mathbb Z/3\mathbb Z\), with labels \((a_i,c_i)\), \((\tilde a_i,\tilde c_i)\), respectively. Set
\[
 A_i=(a_i\tilde c_i)^2,\quad \widetilde A_i=(\tilde a_i c_i)^2,\quad
 E_i=\tilde c_i^2-c_i^2,\qquad
 U=\prod_i a_i\tilde c_i,\quad V=\prod_i\tilde a_i c_i .
\]
Then \(N\mid U^2-V^2\). If \(U\ne V\), the labels determine all six kernel cubes uniquely, and the left kernel cubes satisfy
\[
 b_i^3=
 \frac{E_i A_{i+1}A_{i+2}+\widetilde A_i E_{i+1}A_{i+2}
                      +\widetilde A_i\widetilde A_{i+1}E_{i+2}}
      {(U^2-V^2)/N}.                                                       \tag{1}
\]
If \(N\) is an odd prime, the denominator in (1) is \(m(U+\sigma V)\) for the unique \(\sigma\in\{1,-1\}\) with \(N\mid U-\sigma V\), where \(m=(U-\sigma V)/N\) satisfies \(|m|\ll T\). When \(U=V\), any fixed set of labels has \(O_\varepsilon(N^\varepsilon)\) lifts to kernel cycles.

To prove this, multiply the two equations at \(d_i\) by \(\tilde c_i^2\) and \(c_i^2\) and subtract. With \(B_i=b_i^3\) this yields
\[
                         A_i B_i-\widetilde A_i B_{i+1}=N E_i .
\]
Successive substitution proves (1) before dividing, with \((U^2-V^2)B_i\) equal to \(N\) times its numerator. Coprimality gives \((B_i,N)=1\), so \(N\mid U^2-V^2\). If the determinant is nonzero this yields \(B_i\) uniquely and each original equation then yields \(d_i^3\). The condition that these quotients are cubes of positive squarefree integers is an exact factor-matching constraint on the labels. For prime \(N\), both \(U,V\) are units, so exactly one sign works and factorization gives the asserted denominator. The bound on \(m\) follows from \(U,V\ll T^6\). For arbitrary \(N\), at each odd prime power dividing \(N\) at most one sign can work. Composite moduli can require different signs on different prime-power factors, so a single congruence \(N\mid U-\sigma V\) cannot be assumed.

If \(U=V\), the first reduced equation still imposes the diagonal cubic
\(A_1 b_1^3-\widetilde A_1 b_2^3=N E_1\). Here \(E_1\ne0\): equality of the two \(c\) labels would give \(a_1^2 b_1^3=\tilde a_1^2 b_2^3\) and hence \(b_1=b_2\), contrary to a simple cycle. Heath-Brown's coefficient-uniform diagonal-cubic bound with third coordinate 1 gives \(O_\varepsilon(N^\varepsilon)\) choices for \(b_1,b_2\). The other cubes then follow uniquely from the linear equations. [Heath-Brown, Theorem 1.3](https://arxiv.org/abs/2601.07817v1)

In the singular case there is also an exact factorization of the edge ratios. The rational relation \(U=V\) lets us write them as \(p_i/q_j\) with positive integers attached to left and right vertices after a common clearing of denominators. This follows by assigning a value at one vertex and propagating along the six edges; consistency at the return is exactly \(U=V\). Put \(g_{ij}=(p_i,q_j)\) on each edge. Coprimality of the original roots gives the labels \(p_i/g_{ij}\), \(q_j/g_{ij}\). The edge equations become \(b_i^3p_i^2+d_j^3q_j^2=N g_{ij}^2\), so the alternating sums of the three \(g_{ij}^2\) on either parity agree. This exact condition remains in addition to the rational product relation.

**Open question.** Bound the number of simple six-cycles in each balanced kernel graph by \(O_\varepsilon(T^{2+\varepsilon})\). The uniform bound from the current codegree theorem is \(O_\varepsilon(T^{3+\varepsilon})\), by choosing three left kernels and their shared neighbors. In the nonsingular case the congruence and \(|m|\ll T\) leave the simultaneous perfect-cube tests in (1) and the right-kernel tests. In the singular case the exact ratio and square-sum identities remain. A bound of the requested size still needs to control the number of label tuples passing these tests.

## Restricted sectors with pointwise counts

For the literal constants below, take all four variables in the common interval \( [T,2T) \), with \(T\ge2\) and \(c_0T^5\le N\le C_0T^5\) for fixed positive balance constants. The same exponents hold for fixed-comparable boxes after replacing the numerical constants. Keep one selected actual representation per edge throughout. Counts are of simple cycles; changing the cyclic starting point or orientation costs only an absolute factor. All implied constants may depend on the fixed balance constants in \(N\asymp T^5\). A bound for arbitrary congruence witnesses does not replace an actual-edge count.

### Prime quotient and proper rational powers

**Lemma PrimeQuotientFilter.** If \(N\) is an odd prime and \(U\ne V\), the unique sign and nonzero integer quotient in the preceding proposition satisfy
\[
 U-\sigma V=mN,\qquad 0<|m|\ll T,\qquad \gcd(U,V)\mid m.                    \tag{2}
\]
For \(\sigma=-1\), one has \(m\asymp T\). In particular, a sublinear cutoff \(|m|=o(T)\) has only the positive-sign branch for sufficiently large \(T\).

Indeed, all selected roots, and therefore \(U,V\), are units modulo \(N\). The common divisor \(h=\gcd(U,V)\) is coprime to \(N\), so \(h\mid mN\) implies \(h\mid m\). Also \(T^6\le U,V<64T^6\); hence \((U+V)/N\asymp T\). Primality is essential to the global sign assertion: composite targets can have different signs at their prime-power factors.

**Lemma UniformPowerDifferenceCount.** For a nonzero integer \(D\), the number of positive pairs satisfying \(R^2-S^2=D\) is \(O(\tau(|D|))\). The same bound holds for \(R^2+S^2=D>0\). For every odd integer \(j\ge3\), uniformly in \(j\), both \(R^j-S^j=D\ne0\) and \(R^j+S^j=D>0\) have \(O(\tau(|D|))\) positive solutions.

The difference of squares factors into two integer divisors of \(D\). For the sum of squares, unique factorization in \(\mathbb Z[i]\) gives at most four times the product of \(e+1\) over split prime powers \(p^e\Vert D\): inert primes have forced even exponents and the ramified prime adds only a constant factor. This is \(O(\tau(D))\), including imprimitive pairs. For odd \(j\), first take a positive difference by interchanging \(R,S\) if necessary. Then \(d=R-S\mid D\), and \((S+d)^j-S^j\) is strictly increasing for \(S>0\). Each divisor gives at most one pair. In the sum case \(d=R+S\mid D\), and \(R^j+(d-R)^j\), for \(0<R<d\), is strictly convex and symmetric about \(d/2\). Each divisor gives at most two pairs. These constants are independent of \(j\).

**Theorem ProperPowerCycleSector.** For an odd prime target \(N\), the number of nonsingular actual cycles for which \(U/V\) is a positive rational proper power and \(0<|m|\le M\) is
\[
             O_\varepsilon(MN^\varepsilon),\qquad 1\le M\ll T.            \tag{3}
\]
Consequently this entire sector has \(O_\varepsilon(T^{1+\varepsilon})\) cycles. The exponent of the rational power may vary without any fixed bound.

Write \(U=hR^j,V=hS^j\), where \(h=\gcd(U,V)\), \((R,S)=1\), and \(j\ge2\). Equation (2) gives
\[
                    R^j-\sigma S^j=(m/h)N,\qquad h\mid m.                \tag{4}
\]
Every even exponent is already covered by \(j=2\). For an odd exponent, \(U\ne V\) ensures that at least one of \(R,S\) is at least two, so \(j\le\log_2(64T^6)=O(\log N)\). For each fixed nonzero \(m\), there are only \(\tau(|m|)\) choices of \(h\). The preceding lemma, applied to the nonzero integer on the right of (4), gives \(O_\varepsilon(N^\varepsilon)\) product pairs \((U,V)\), after absorbing the logarithmic exponent count and reallocating \(\varepsilon\). Negative differences are handled by interchanging the variables; a negative sum has no positive solution.

For fixed \(U,V\), the twelve ordered positive root labels have at most \(d_6(U)d_6(V)\ll_\varepsilon N^\varepsilon\) factorizations. Here \(d_6\) counts ordered factorizations into six positive factors, and \(U,V\ll T^6\). The nonzero-determinant equations (1) give at most one kernel-cycle lift for each label tuple. All interval, coprimality and squarefree-cube tests only discard tuples. Thus each fixed \(m\) contributes subpower many actual cycles, and summing over the at most \(2M\) nonzero quotients proves (3). This argument does not count the complementary non-power ratios.

### Kernel overlap

**Proposition LowKernelLcmSector.** For every target \(N\), the number of actual simple cycles with
\[
                    \operatorname{lcm}(b_1,b_2,b_3)\le L
\]
is \(O_\varepsilon(N^\varepsilon L)\), for \(L\ge1\). The same bound holds on the right. In particular, the sector where either side has lcm at most \(8T^2\) has \(O_\varepsilon(T^{2+\varepsilon})\) cycles, irrespective of the determinant.

For a fixed integer \(h\le\min(L,8T^3)\), the ordered left triples of exact lcm \(h\) are among the at most \(\tau(h)^3\) triples of divisors of \(h\). Summing gives \(O_\varepsilon(N^\varepsilon L)\) triples. Each fixed distinct triple lifts to subpower many cycles by the three actual common-neighbor bounds. Repeated right vertices are discarded. The argument on the right is identical. Moreover, if \(g=\gcd(b_1,b_2,b_3)\), then the lcm is at most \(b_1b_2b_3/g^2<8T^3/g^2\). Thus a common kernel gcd at least \(T^{1/2}\) is contained in the displayed sector.

**Proposition SharedNeighborOverlapFilter.** Suppose distinct left kernels \(b_1,b_2\) share an actual right neighbor \(d\), with selected roots \((a_1,c_1),(a_2,c_2)\), and put \(g=\gcd(b_1,b_2)\). Then
\[
 (g,c_1c_2d)=1,\qquad c_1\ne c_2,\qquad g^3\mid c_2^2-c_1^2.             \tag{5}
\]
For fixed \(g\), the number of possible ordered off-diagonal root pairs \((c_1,c_2)\in[T,2T)^2\) in (5) is
\[
                         O_\varepsilon(N^\varepsilon T^2/g^3).          \tag{6}
\]
Every odd prime \(p\mid g\) also satisfies \(p^3<4T\).

Coprimality of each pair of summands proves the first assertion. Subtraction of the edge equations gives (5). Equality of the two \(c\) labels would force \(a_1^2b_1^3=a_2^2b_2^3\); prime valuations modulo two would give \(b_1=b_2\), a contradiction. Let \(g_0\) be the odd part of the squarefree integer \(g\). Each odd \(p\mid g_0\) divides exactly one of \(c_2-c_1,c_2+c_1\). Splitting these primes gives coprime \(r,s\) with \(rs=g_0\), \(r^3\mid c_2-c_1\), \(s^3\mid c_2+c_1\). These are nonzero multiples of sizes less than \(T\) and \(4T\); their number is \(O(T^2/(rs)^3)\). An interval shorter than its first nonzero multiple contributes zero. Sum over at most \(2^{\omega(g_0)}\) splits and use \(g^3/g_0^3\le8\) to obtain (6). The same nonzero-multiple observation proves the prime-size restriction. No full cycle-frequency saving follows from this local constraint alone.

### Singular coherent labels and divisor reconstruction

**Lemma PrimitiveCoherentHeight.** Every singular actual cycle has unique positive integral coherent vertex labels \(p_i,q_j\) with global gcd one and
\[
 \frac{a_{ij}}{c_{ij}}=\frac{p_i}{q_j},\qquad
 g_{ij}=\gcd(p_i,q_j),\qquad p_i=g_{ij}a_{ij},\quad q_j=g_{ij}c_{ij}.       \tag{7}
\]
Every coherent label and edge gcd is a unit modulo \(N\), for arbitrary \(N\). They satisfy
\[
 p_i,q_j<(2T)^5\ll N,\qquad g_{ij}<32T^4.                                 \tag{8}
\]
The intrinsic integer
\[
                             g_*(C)=\min_{(i,j)\text{ edge}}g_{ij}        \tag{9}
\]
is independent of scale and starting point, and every edge satisfies \(g_*\le g_{ij}<8g_*\).

The singular ratio identity supplies rational labels by propagation on the cycle. Connectedness makes them unique up to a common positive rational scalar. Clearing denominators and dividing the global gcd gives the unique primitive integral labeling. Original root coprimality yields (7). For a prime \(\ell\mid N\), all original roots are \(\ell\)-units, so (7) equates the valuations at the endpoints of every edge. Connectedness equates all six valuations; primitive normalization makes their common value zero.

For (8), orient a five-edge spanning tree away from one root vertex. The ratio child/parent on an edge is a fraction with numerator and denominator below \(2T\). Give the root the product of the five denominators. Each resulting vertex value is the product of the numerators on its path and the denominators off that path: exactly five factors, all below \(2T\). This is an integral coherent labeling, and primitive normalization can only decrease its entries. Dividing by the original roots, which are at least \(T\), bounds the edge gcds. At an edge attaining \(g_*\), both endpoints are below \(2Tg_*\). Every vertex is within two edges of one of these endpoints, and every edge ratio is between \(1/2\) and two. All labels are consequently below \(8Tg_*\), proving the final comparison.

**Lemma SingularDivisorReconstruction.** Fix a selected actual edge and a positive integer \(s\) prescribed to be its intrinsic coherent gcd. The number of ordered actual singular cycles extending it is \(O_\varepsilon(N^\varepsilon)\), uniformly in \(s\).

Set \(L_i=b_i^3p_i^2\), \(R_j=d_j^3q_j^2\). The actual edge equations become
\[
                             L_i+R_j=Ng_{ij}^2.                           \tag{10}
\]
Every positive integer has at most one decomposition \(n=b^3p^2\) with \(b\) squarefree: a prime exponent \(e\) forces its exponent in \(b\) to be \(e\bmod2\), and its exponent in \(p\) is then \((e-3(e\bmod2))/2\). If this is negative, there is no such decomposition.

Start from \((b_1,d_1)\) and its selected roots \(a_{11},c_{11}\). The prescribed gcd gives \(p_1=sa_{11}\), \(q_1=sc_{11}\), hence \(L_1,R_1\). Choose successively
\[
 g_{21}\mid q_1,\qquad g_{22}\mid p_2,\qquad
 g_{32}\mid q_2,\qquad g_{33}\mid p_3.                                    \tag{11}
\]
The first choice forces \(L_2=Ng_{21}^2-R_1\), whose unique possible decomposition supplies \(b_2,p_2\). The second forces \(R_2=Ng_{22}^2-L_2\), supplying \(d_2,q_2\); the next two similarly supply \(b_3,p_3,d_3,q_3\). Reject nonpositive potentials, nonexistent decompositions, labels violating (8), and candidates violating the balanced intervals, selected-edge convention, simplicity or primitive normalization. Prune a recovered label to (8) before making the next divisor choice. No actual cycle is lost.

Each dividend in (11) is \(O(N)\), so the four divisor choices give subpower many branches. The final edge gcd is already determined as \(\gcd(p_1,q_3)\); check \(L_1+R_3=N\gcd(p_1,q_3)^2\). This proves the count. In particular, the six gcd-square identity alone is not used as a substitute for the exact edge equations.

**Theorem SingularWeightedAndSmallGcdSector.** For all targets \(N\) and all \(Y\ge1\),
\[
 \#\{C\text{ singular}:g_*(C)\le Y\}
       \ll_\varepsilon N^\varepsilon e(G)Y,\qquad
 \sum_{C\text{ singular}}\frac1{g_*(C)}
       \ll_\varepsilon N^\varepsilon e(G).                                \tag{12}
\]
Thus the singular sector \(g_*\le T^{1/2}\) has \(O_\varepsilon(T^{2+\varepsilon})\) cycles, and \(g_*=1\) has \(O_\varepsilon(T^{3/2+\varepsilon})\) cycles.

Orient a cycle from an edge attaining its minimum gcd. For each selected edge and each \(s=g_*\le Y\), the reconstruction lemma gives subpower many cycles; this proves the first bound. For the weighted bound sum these counts with weight \(1/s\) up to \(32T^4\), using (8). The harmonic factor is absorbed into \(N^\varepsilon\). Actual codegree bounds and the elementary bipartite edge count give \(e(G)\ll_\varepsilon T^{3/2+\varepsilon}\), proving the displayed sector consequences after reallocating \(\varepsilon\). These are counted divisor reconstructions for actual edges; no metric gcd-graph theorem is a premise.

### Exact remaining domain and both implications

**Proposition ResidualCycleEquivalence.** For arbitrary targets let \(\mathcal R_N(T)\) consist of the actual cycles for which both kernel-side lcms exceed \(8T^2\) and either

- \(U\ne V\); or
- \(U=V\) and \(g_*>T^{1/2}\).

Then a uniform bound \(\#C_6(G)\ll_\varepsilon T^{2+\varepsilon}\) for all balanced targets is equivalent to the same uniform bound for \(\#\mathcal R_N(T)\). Both statements remain open. For odd prime targets the analogous equivalence holds after further requiring \(U/V\) to be a non-power rational number on the nonsingular part of \(\mathcal R_N(T)\).

The forward directions hold because these are subsets of the actual cycles. Conversely, all cycles outside \(\mathcal R_N(T)\) lie either in a low-lcm sector or in the small-coherent-gcd singular sector. Each has just been bounded by \(O_\varepsilon(T^{2+\varepsilon})\); adjoining a bound for the residual sector gives the full bound. On odd prime targets the removed proper-power sector costs only \(O_\varepsilon(T^{1+\varepsilon})\) by (3), so the prime specialization has the same converse. These implications do not establish the remaining estimate or an improved representation exponent.

Large coherent \(g_*\) has not been shown to force small kernel lcm. The weighted estimate in (12) cannot be turned into an unweighted bound for large \(g_*\) without another input. A cubic completion theorem for fixed labels likewise does not count how often the nonsingular cube tests in (1) hold. The remaining arithmetic problem retains these exact constraints.
