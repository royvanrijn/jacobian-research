# Outward continuation and matching products in critical six-cycles

Two consecutive selected edges can be extended in opposite directions, using the same seed for both continuation congruences. For singular cycles, the primitive coherent labels give a second restriction. A large continuation divisor uses root height on one alternating matching. The two matching products therefore give a lower cutoff for the coherent gcd than the combined root height.

Work in the coefficient-shaped selected graph and coefficient family of the [critical rectangular criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md). Write \(A=2n^{11/51}\), \(T=A^{1/33}\), and \(R_1,R_2\) for the two shaped root caps \(A_{1,n},A_{2,n}\). Let the vertices be \(v_i\), the edges \(e_i=v_iv_{i+1}\) with indices modulo six, and set
\[
 d_i=\min\bigl(|R-R'|/\gcd(R,R'),\,n^{11/51}\bigr)
\]
at a vertex with selected roots \(R,R'\). The gcd of the two opposite roots on this path is \(\alpha_i\). The quantities \(F_i=F_{v_i}\), \(K_i\), and \(g_*\) on singular cycles are those of the [rectangular continuation count](ERDOS_943_CRITICAL_RECTANGULAR_CONTINUATION.md).

## Extending both ends of a two-edge path

On singular cycles set \(m_i^+=F_iK_{i-2}\), \(m_i^-=F_iK_{i-1}\). The signs refer to forward and reverse order respectively. Use the step costs of the [joint continuation count](ERDOS_943_JOINT_RECTANGULAR_CONTINUATION.md):
\[
 H_i^\pm=\operatorname{lcm}(\alpha_i,m_i^\pm),\qquad
 L_i^\pm=\max(1,d_i/(H_i^\pm)^2),\qquad
 S_i=\max(1,d_i/\alpha_i^2).
\]
On nonsingular cycles put \(H_i^\pm=\alpha_i\).

**Proposition (Outward continuation count).** For every \(Q\ge1\), the cycles satisfying \(L_i^-L_{i+2}^+\le Q\) for some \(i\) number \(O_\varepsilon(n^\varepsilon n^{64/153}Q)\). This applies to both determinant branches and includes the critical budget at \(Q=T^{14}\).

Choose the two entering edges \(e_i,e_{i+1}\) at cost \(O_\varepsilon(n^{64/153+\varepsilon})\) by choosing their two outer kernels and an actual common neighbor. Read these edges in reverse to extend at \(v_i\), and in forward order to extend at \(v_{i+2}\). Each extension is from two known edges. Its divisor branches are recovered from those edges, and its reduced-difference residue list gives the stated \(L_i^-\) or \(L_{i+2}^+\) cost up to a subpower factor. Fix dyadic ranges of the two reduced differences and the two moduli before counting, as in the joint continuation proof. In a box admitting the product condition, the product of the two branch counts is \(O_\varepsilon(n^\varepsilon Q)\). The four-edge path \(e_{i-1},e_i,e_{i+1},e_{i+2}\) closes by the actual codegree bound on its endpoints.

The inequality is at two same-side vertices with opposite step orientations. For instance, in the clipped size inequalities used below, two values \(f_0=f_2=17\) with all other \(f_i,k_i=0\) satisfy the consecutive non-counting inequalities and violate the outward non-counting inequality. This is an example for the size variables.

## A bound involving both matchings

Write \(g_i=\gcd(h_i,h_{i+1})\) for the edge gcds of the primitive coherent labels, and put \(F=\prod_i F_i\). Let \(\mathcal A,\mathcal B\) be the even and odd edge matchings. Write \(K_{\mathcal A}=\prod_{i\ {\rm even}}K_i\) and \(K_{\mathcal B}=\prod_{i\ {\rm odd}}K_i\). Let \(R_{\mathcal A},R_{\mathcal B}\) be the products of the six selected roots on each matching.

**Proposition (Matching product inequality).** For a singular cycle the two ratios
\[
 P_{\mathcal A}=\frac{R_{\mathcal A}}{F K_{\mathcal A}^{2}},\qquad
 P_{\mathcal B}=\frac{R_{\mathcal B}}{F K_{\mathcal B}^{2}}
\]
are positive integers and
\[
 \prod_i g_i\le
       F^4(K_{\mathcal A}K_{\mathcal B})^3
                    \sqrt{P_{\mathcal A}P_{\mathcal B}}\min(P_{\mathcal A},P_{\mathcal B}).
 \tag{1}
\]
In particular
\[
                      g_*\le R_1 R_2\,\bigl(F\min(K_{\mathcal A},K_{\mathcal B})\bigr)^{1/3}.
 \tag{2}
\]

At each prime-power level mark the vertices whose coherent labels are divisible by that prime power. The marked set is proper by primitivity. Let \(e\) be the number of edges with both endpoints marked, and let \(c_{\mathcal A},c_{\mathcal B}\) be the numbers of edges in each matching with exactly one marked endpoint. This level contributes \(e\) to the exponent of \(\prod_i g_i\), and \(c_{\mathcal A},c_{\mathcal B}\) to the exponents of \(R_{\mathcal A},R_{\mathcal B}\). The latter facts follow from the exact selected roots \(h_i/g_i,h_{i+1}/g_i\).

Five marked vertices give \(e=4\) and one cut edge in each matching. These levels are counted by \(F\), where the missing vertex determines \(F_i\). Four consecutive marked vertices give \(e=3\) and two cut edges in one matching. These levels are counted by \(K_i\) when the two missing vertices are \(v_{i+1},v_{i+2}\), and the cut edges \(e_i,e_{i+2}\) are in the matching of \(i\). Indeed the four complementary labels satisfy the divisibility defining \(J_i\). If \(J_i\) is divisible by the prime, at least one of the two other labels is a unit. Subtracting the valuation of \(\gcd(J_i,h_{i+1}h_{i+2})\) then removes exactly the levels at which at least one additional vertex is marked.

Remove these two kinds of level. This proves the asserted integrality of \(P_{\mathcal A},P_{\mathcal B}\). In every remaining level \(e\le2\). If both cut counts are positive, then \(e\le(3c_{\mathcal A}+c_{\mathcal B})/2\) and \(e\le(c_{\mathcal A}+3c_{\mathcal B})/2\). If one cut count vanishes, the marked vertices are a union of pairs of that matching. The only nonempty remaining possibility is one pair: two pairs leave an adjacent pair unmarked, making the marked four consecutive. Hence \(e=1\) and the other cut count is two. Both inequalities still hold. Sum the two inequalities separately over all primes and levels. Restore the contributions of \(F\) and \(K_{\mathcal A}K_{\mathcal B}\), and take the smaller bound to get (1).

Each root product is at most \((R_1R_2)^3\). Substituting this for both \(R_{\mathcal A}\) and \(R_{\mathcal B}\) in (1) gives
\[
             g_*^6\le (R_1R_2)^6 F^2\min(K_{\mathcal A},K_{\mathcal B})^2.
\]
This proves (2).

## The coefficient-sensitive singular cutoff

**Proposition (Improved critical gcd cutoff).** The singular cycles with
\[
                         g_*\ge \frac{8n^{2749/4590}}{\kappa_1\kappa_2}
 \tag{3}
\]
number \(O_\varepsilon(n^{26/51+\varepsilon})\).

Delete the counted path sectors with some \(\alpha_i\ge n^{43/306}\), with a consecutive two-step product at most \(T^{14}\) or a consecutive three-step product at most \(T^{321/10}\) in either orientation, and with an outward product at most \(T^{14}\). The consecutive sectors follow from the [joint continuation proof](ERDOS_943_JOINT_RECTANGULAR_CONTINUATION.md) with thresholds multiplied by fixed constants. The first sector is counted by [the large opposite-root gcd bound](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#4-actual-rectangular-codegrees-and-opposite-root-gcd-frequency).

Set \(f_i=\log_T F_i\), \(k_i=\log_T K_i\). We have \(0\le f_i\le43/2\) since \(F_i\mid\alpha_i\), and \(0\le f_i+k_{i-2},f_i+k_{i-1}\le33\) by the [two orientations of the root divisibilities](ERDOS_943_CRITICAL_RECTANGULAR_CONTINUATION.md#coefficient-uniform-sectors). Write \(f_*=\sum_i f_i\), \(k_{\mathcal A}=\sum_{i\ {\rm even}}k_i\), \(k_{\mathcal B}=\sum_{i\ {\rm odd}}k_i\). The root products give \(f_*+2k_{\mathcal A},f_*+2k_{\mathcal B}\le198\). Put
\[
 u_i=\min(f_i,33/2),\quad
 z_i^+=\min(f_i+k_{i-2},33/2),\quad z_i^-=\min(f_i+k_{i-1},33/2).
\]
Since \(d_i\le A\), the actual step costs are bounded by \(T^{33-2u_i}\) and \(T^{33-2z_i^\pm}\). The failures of the five continuation tests thus imply, for every \(i\),
\[
\begin{array}{rclcrcl}
 z_i^++z_{i+1}^+&\le&26, && z_i^-+z_{i-1}^-&\le&26,\
 z_i^-+z_{i+2}^+&\le&26,\
 u_i+z_{i+1}^++z_{i+2}^+&\le&669/20, &&
 u_i+z_{i-1}^-+z_{i-2}^-&\le&669/20 .
\end{array}                                                         \tag{4}
\]
Weak inequalities suffice for this upper bound.

These rational inequalities give
\[
                        2(f_*+\min(k_{\mathcal A},k_{\mathcal B}))\le769/5.
 \tag{5}
\]
For the finite linear optimization, mark the positions reaching \(33/2\) in each of \(z^+,z^-,u\). The two \(z\) sets are independent in the six-cycle by their adjacent tests. A reverse high position at \(i\) forbids a forward high at \(i+2\), and the \(u\) high set is a subset of their intersection. Enumerate these triples of sets. In each branch the clipped variable is the constant \(33/2\) at a marked position and its own underlying linear form elsewhere. Impose the respective lower or upper clipping inequality to cover the boundary. Maximize \(2(f_*+b)\) with \(0\le b\le k_{\mathcal A},k_{\mathcal B}\) in each resulting rational polytope. The [exact checker](verify_erdos943_outward_continuation.py) uses SageMath 10.9 with the rational PPL backend. Replay from the repository root:

    sage -python research/arithmetic/verify_erdos943_outward_continuation.py

The checker covers 250 branches. Grouped by the number of high \(u\) positions \(0,1,2\), the respective maxima are \(764/5,719/5,769/5\).

The largest maximum occurs with two high \(u\) positions, which must be opposite by the adjacent and outward inequalities. Take them at 0 and 3. Every other vertex is adjacent to one of these two and its two step exponents are at most \(9.5\), so none of those exponents is clipped. The four triple inequalities starting from these high positions in their two orientations give
\[
 f_1+f_2+f_4+f_5+k_0+k_3+\max(k_1,k_5)+\max(k_2,k_4)\le 339/10 .
\]
The terms involving \(k\) on the left sum to at least \(\tfrac12\sum_i k_i\), and \(f_0+f_3\le43\). This gives the bound \(769/5\) in this branch. The vector with \(k_i=0\) and \(f=(21.5,7.45,9.5,21.5,7.45,9.5)\) attains it in the weak size inequalities.

By (2) and (5), every singular cycle surviving the counted tests has
\[
                    g_*\le R_1R_2 T^{769/30}
                         =\frac{2^{2749/990} n^{2749/4590}}{\kappa_1\kappa_2}
                         <\frac{8 n^{2749/4590}}{\kappa_1\kappa_2}.
\]
The removed sectors all meet the stated cycle budget.

## Remaining pointwise estimate

**Corollary (Outward-continuation residual equivalence).** In the [matching-ratio residual](ERDOS_943_MATCHING_RATIO_RECTANGULAR.md#removing-these-sets-from-the-cycle-problem), one can delete the outward sector at \(Q=T^{14}\) and the singular sector (3). The \(O_\varepsilon(n^{26/51+\varepsilon})\) count on the surviving set is equivalent to (H) for the same shaped graph and coefficient family.

One direction follows by restriction. The converse adds the counted sectors and uses the matching-ratio equivalence. The surviving nonsingular cycles and the singular cycles below the new cutoff still require a pointwise count.
