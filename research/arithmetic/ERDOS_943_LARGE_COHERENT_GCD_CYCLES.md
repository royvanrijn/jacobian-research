# Large coherent gcds in singular square-full six-cycles

Fix the balanced graph of [actual selected coprime square-full representations](ERDOS_943_SIX_CYCLES.md), with \(a,b,c,d\in[T,2T)\) and \(N\asymp T^5\). Put \(B=2T\). The known [coherent labeling](ERDOS_943_SIX_CYCLES.md#singular-coherent-labels-and-divisor-reconstruction) on a singular simple six-cycle assigns positive integers \(h_v\) to its six vertices, with global gcd one. On edge \(uv\), set \(g_{uv}=\gcd(h_u,h_v)\). Its original selected roots at \(u,v\) are \(h_u/g_{uv},h_v/g_{uv}\). The intrinsic coherent gcd \(g_*\) is the smallest edge gcd.

**Theorem.** Uniformly in the target and every \(G\ge1\),
\[
 \#\{C:\ C\text{ singular},\ g_*(C)\ge G\}
 \ll_\varepsilon N^\varepsilon
       \min\{B^3,\ B^{3/2}+B^{51/4}/G^3\}.                           \tag{1}
\]
In particular the count is \(O_\varepsilon(T^{2+\varepsilon})\) for \(G\ge B^{43/12}\), and \(O_\varepsilon(T^{3/2+\varepsilon})\) for \(G\ge B^{15/4}\). The constant depends on the balance constants. The argument uses the existing actual codegree and coupled path-gcd bounds.

## Common divisors on four edges

For a vertex \(v\), let \(F_v\) be the gcd of the four edge gcds off \(v\). These edges form a path covering the other five vertices. If \(\ell^k\mid F_v\), the five corresponding labels are all divisible by \(\ell^k\), so \(\ell\nmid h_v\) by primitive normalization. It follows that for each edge \(uv\) the selected root at \(u\) is divisible by \(F_v\). In particular the two roots at the outer vertices of the path through \(v\) have a common divisor at least \(F_v\). Each \(F_v<B\), and the six \(F_v\) are pairwise coprime.

For each prime \(\ell\) and level \(k\ge1\), mark the edges whose gcd is divisible by \(\ell^k\). Let \(m\) be the number of marked edges and \(j\) the number of changes of status between consecutive edges around the cycle. A nonempty marking is proper because the labels have global gcd one. If it consists of one run, it has \(j=2\) and length at most four. A run of five would force all six vertices to be divisible by \(\ell^k\). If there are two or more runs, \(j\ge4\) and \(m\le4\). Consequently
\[
 m\le \tfrac32 j+\boldsymbol{1}_{\{\text{one run of length four}\}} .
\]
In the exceptional case the run is the complement of the edges at a unique \(v\), and its indicator contributes exactly one to \(\nu_\ell(F_v)\). Let \(R_v,R'_v\) be the selected roots at \(v\) on its two incident edges, and \(\gamma_v=\gcd(R_v,R'_v)\). The sum of changes of status at \(v\) over levels is
\[
 |\nu_\ell(g_{e})-\nu_\ell(g_{e'})|
       =\nu_\ell(R_vR'_v/\gamma_v^2).
\]
Summing over levels and then primes proves the exact product inequality
\[
             \prod_{e}g_e
   \le \left(\prod_v R_v R'_v/\gamma_v^2\right)^{3/2}\prod_v F_v
   \le B^{18}\prod_v F_v .                                      \tag{2}
\]
In particular high coherent content forces a large product of actual root divisors.

The exponent three in the height bound \(g_*\le B^3\) when all \(F_v=1\) is sharp for primitive gcd label data. Choose six pairwise coprime integers \(q_i\) in \([T,2T)\) and take \(h_i=q_iq_{i+1}q_{i+2}q_{i+3}\), with indices modulo six. Such choices exist at arbitrarily large scales by translating six fixed pairwise coprime integers by a common multiple of their pairwise differences. The edge \(i,i+1\) has gcd \(q_{i+1}q_{i+2}q_{i+3}\) and selected-root analogues \(q_i,q_{i+4}\). Every prime support is a run of three edges, \(F_v=1\) for every vertex, and equality holds in the first inequality in (2). This construction represents primitive coherent label data. The common-target addition equations and square-free kernels are extra arithmetic conditions.

## Reconstructing from consecutive large divisors

Put \(L_v=\max(1,B/F_v^2)\). The cycles for which some three consecutive vertices \(v_1,v_2,v_3\) satisfy, for a given \(Y\ge1\),
\[
                             L_{v_1}L_{v_2}L_{v_3}\le Y     \tag{3}
\]
number \(O_\varepsilon(N^\varepsilon B^{3/2}Y)\). To see this, start from a selected edge entering \(v_1\) and choose \(F_{v_1}\) as a divisor of its root at the preceding vertex. This costs \(O_\varepsilon(N^\varepsilon)\) choices per edge. At a current vertex \(v\) with known kernel \(b_v\) and divisor \(F_v=f\), the root \(R\) at \(v\) on the next edge satisfies
\[
                              b_v^3 R^2\equiv N\bmod f^2
\]
because the other root is a multiple of \(f\). For any actual edge both \(b_v\) and \(N\) are units modulo \(f\). There are \(O_\varepsilon(N^\varepsilon)\) square-root classes modulo \(f^2\), including the class for \(f=1\). Each allowed \(R<B\) gives at most one actual next edge by unique decomposition of the complementary square-full summand. Once the edge is known, choose \(F_{\rm next}\) from the divisors of the root at \(v\). Thus a step costs \(O_\varepsilon(N^\varepsilon L_v)\), including the choice of the next divisor. Fix the dyadic ranges of the three divisors first. The product of step bounds in a dyadic box admitting (3) is at most a constant times \(Y\), with a subpower allowance. There are logarithmically many boxes.

These three steps give four consecutive edges. The last two are a common-neighbor problem and have \(O_\varepsilon(N^\varepsilon)\) choices by the [actual codegree theorem](ERDOS_943_KERNEL_GRAPH.md#common-neighbors-and-reconstruction). The condition of that theorem holds for all sufficiently large \(T\) at the fixed balance constants. The two codegree bounds also give \(e(G_N)\ll_\varepsilon N^\varepsilon B^{3/2}\) for the whole graph \(G_N\) by counting unordered pairs of neighbors on one side. This proves the count on (3).

## The coherent gcd tail

Outside (3), put \(f_v=\log_B F_v\), \(d_v=(1-2f_v)_+\) and \(y=\log_B Y\). Every three consecutive \(d_v\) sum to more than \(y\). Adding the six inequalities gives \(\sum_v d_v>2y\), and at most four of the \(d_v\) vanish. If \(F_{\max}=\max_v F_v\le B^{1/2}\), then \(\sum_v f_v=3-\frac12\sum_v d_v<3-y\), and (2) implies \(g_*^6<B^{21}/Y\). If \(F_{\max}>B^{1/2}\), the at most four vertices with \(d_v=0\) give
\[
         \prod_v F_v\le (B^3/Y)(F_{\max}/B^{1/2})^4
                          = B F_{\max}^4/Y.                   \tag{4}
\]
So whenever \(G^6\ge B^{21}/Y\), a cycle with \(g_*\ge G\) outside (3) has \(F_{\max}^4\ge YG^6/B^{19}\).

At the vertex realizing \(F_{\max}\), use the [coupled path weight](ERDOS_943_COUPLED_PATH_ROOT_GCD.md#9-weighted-path-and-cycle-sums). For the two roots \(c_1,c_2\) at that vertex and the two roots \(a_1,a_2\) at the outer vertices, after swapping the two sides if necessary, let \(\alpha=\gcd(a_1,a_2)\), \(\gamma=\gcd(c_1,c_2)\), and \(e=|c_1^2-c_2^2|/(\alpha^2\gamma^2)\). Here the selected summands are primitive and \(e\) is a positive integer. Its weight is \(W=\alpha^3/\sqrt e\). Because \(\alpha\ge F_{\max}\) and \(\sqrt e\le B/(\alpha\gamma)\), this gives \(W\ge YG^6/B^{20}\).

The cited path weight theorem has a tail \(O_\varepsilon(N^\varepsilon B^4/Z)\) for a shared-vertex weight at least \(Z>0\) in this box. Its actual graph condition \(2B^4<N\) holds for sufficiently large \(T\). Combining it with the count on (3) gives
\[
 \#\{g_*\ge G\}\ll_\varepsilon N^\varepsilon
             \left(B^{3/2}Y+\frac{B^{24}}{YG^6}\right)
                         \quad \text{when } G^6Y\ge B^{21}.      \tag{5}
\]
For \(G\ge B^{13/4}\), take \(Y=\max(1,B^{45/4}/G^3)\), which meets this condition. This gives \(B^{51/4}/G^3\) when the second argument of the maximum is at least one, and \(O(B^{3/2})\) otherwise. For \(G<B^{13/4}\) the trivial codegree bound \(O_\varepsilon(N^\varepsilon B^3)\) suffices by choosing three same-side vertices. This also supplies the upper cap in (1).

## Congruences from three complementary edges

Index vertices \(v_i\) and edges \(e_i=v_i v_{i+1}\) modulo six, and write \(h_i=h_{v_i}\). For the path \(e_i,e_{i+1},e_{i+2}\), put
\[
 J_i=\gcd(g_{e_{i+3}},g_{e_{i+4}},g_{e_{i+5}}),\qquad
 K_i=\frac{J_i}{\gcd(J_i,h_{i+1}h_{i+2})}.                         \tag{6}
\]
Both outer roots, at \(v_i\) and \(v_{i+3}\), are divisible by \(K_i\). Indeed at a prime with \(a=
u_\ell(J_i)>0\), all four labels on the complementary path have valuation at least \(a\). One of \(t_1=
u_\ell(h_{i+1}),t_2=
u_\ell(h_{i+2})\) is zero by primitivity. The exponent of \(K_i\) is \(k=\max(0,a-t_1-t_2)\). If \(k>0\), the gcd on each outer edge of the original path has valuation \(t_1\) or \(t_2\), respectively, so the outer roots have valuations at least \(a-t_1,a-t_2\). In particular \(K_i<B\).

**Proposition (Three-edge continuation).** Given the first two actual selected edges of the path, and prescribed values \(f=F_{v_{i+2}}\), \(K=K_i\), the possible third edges in singular cycles lie in \(O_\varepsilon(N^\varepsilon)\) residue classes for their root at \(v_{i+2}\) modulo \((fK)^2\). Consequently the number of such edges is \(O_\varepsilon(N^\varepsilon(1+B/(fK)^2))\). Summing over the possible \(f,K\) incurs only divisor-many branches. The reverse statement uses \(F_{v_{i+1}}\) and the same \(K_i\).

Write the ordered roots on these edges, from left to right in path order, as \((a_0,b_0),(c_0,d_0),(e_0,f_0)\). Their coherent gcds are in the ratio
\[
                   g_{e_i}:g_{e_{i+1}}:g_{e_{i+2}}
                  =c_0e_0:b_0e_0:b_0d_0.
\]
This follows by equating the two coherent labels at each shared vertex. The alternating sum of all six \(g_e^2\) is zero, since each vertex potential \(h_v^2 b_v^3\) occurs with opposite coefficients in the sum of the six edge equations after multiplying the equation on edge \(uv\) by \(g_{uv}^2\). The three complementary gcds are divisible by \(J_i\). The gcd of the three gcds on the original path is coprime to \(J_i\), because the six edge gcds have gcd one. Divide the ratio triple by its gcd, reduce the alternating square sum modulo \(J_i^2\) and multiply by the square of the gcd of the ratio triple to obtain
\[
                   (c_0^2-b_0^2)e_0^2+b_0^2d_0^2\equiv 0
                                                     \pmod{J_i^2}.          \tag{7}
\]

At a prime \(\ell\mid K\) use \(a=t_1+t_2+k\) as above. Since \(h_i,h_{i+3}\) have valuations at least \(a\), the roots \(b_0,e_0\) are units. The middle-edge gcd is a unit, so \(c_0,d_0\) have valuations \(t_1,t_2\), visible from the known edge. If \(t_1=0\), (7) forces \(
u_\ell(c_0^2-b_0^2)=2t_2\), and division gives a unit square congruence for \(e_0\) modulo \(\ell^{2k}\). If \(t_1>0\), then \(t_2=0\), the leading coefficient is a unit, and (7) is a unit square congruence modulo \(\ell^{2a}\). In this latter case \(
u_\ell(f)=t_1\), since the four labels besides the interior ones have valuations at least \(a>t_1\). In the former case \(\ell
mid f\) because \(h_{i+1}\) is a unit. Thus the square congruence at every \(\ell\mid K\) works modulo \(\ell^{2
u_\ell(fK)}\). At primes of \(f\) outside \(K\), the ordinary edge equation gives the unit square congruence modulo \(\ell^{2
u_\ell(f)}\) as in the proof of (3). There are at most four root classes per prime power, including at 2. CRT gives the assertion, and each selected root \(e_0\) gives at most one third edge by exact complementary square-full decomposition. Finally \(K\mid a_0\) and \(f\mid c_0\) enumerate both unknown divisors from the given edges.

**Corollary (Two successive continuations).** Singular cycles for which for some \(i\)
\[
              F_{v_{i+2}}K_i\ge B^{1/2},\qquad
              F_{v_{i+3}}K_{i+1}\ge B^{1/2}                     \tag{8}
\]
number \(O_\varepsilon(N^\varepsilon B^2)\). There are \(O_\varepsilon(N^\varepsilon B^2)\) choices for the first two edges by choosing the outer pair and applying the actual codegree bound. Use the proposition twice, then close the four-edge path by the codegree bound. In particular two adjacent \(F_v\ge B^{1/2}\) give this estimate.

Each level in (2) obeys the stronger profile bound
\[
 m\le j+\boldsymbol{1}_{\{\text{one run of length three}\}}
             +2\boldsymbol{1}_{\{\text{one run of length four}\}}.
\]
A run of length three can only be the sole run. A missing edge with both neighbors marked is forced to be marked by its endpoint valuations, so a length-three run plus a second run cannot fit in six edges. The length-three indicator counts the corresponding level of \(K_i\) by (6). Thus
\[
                        \prod_e g_e\le B^{12}\prod_i K_i\prod_v F_v^2.        \tag{9}
\]
If \(\prod_i K_i\ge B^{9/2}\), two adjacent \(K_i\) are at least \(B^{1/2}\): otherwise at most three reach the threshold, each less than \(B\), and the others are less than \(B^{1/2}\). This is a sector of (8). Consequently for any fixed \(b\ge0\), the singular cycles with \(\prod_v F_v\le B^b\) and \(g_*\ge B^{11/4+b/3}\) have the \(O_\varepsilon(N^\varepsilon B^2)\) bound. The case \(F_v=1\) everywhere reaches \(B^{11/4}\).

## Refined coherent gcd tail

**Corollary (Refined coherent gcd tail).** For \(G\ge B^{10/3}\), the number of singular cycles at \(g_*\ge G\) is
\[
 O_\varepsilon\!\left(N^\varepsilon
       \left(B^2+\frac{B^{45/2}}{G^6}+\frac{B^{43}}{G^{12}}\right)\right).
                                                                    \tag{10}
\]
In particular the \(T^2\) budget holds for \(G\ge B^{41/12}\). Bound (1) is also available, including its \(T^{3/2}\) estimate at the larger cutoff.

To prove (10), first remove (8) in both orientations. Also remove the sector \(\max_v F_v\ge B^{3/4}\), where the weight \(W\ge F_{\max}^4/B\) from the proof of (5) gives the \(T^2\) budget. Write \(f_i=\log_B F_{v_i}\), \(k_i=\log_B K_i\), and let \(A=\{i:f_i\ge1/2\}\). The set \(A\) has no adjacent vertices by (8). If \(v_i\) has a neighbor in \(A\), then in the remainder
\[
                    k_{i-2}+f_i<\tfrac12,\qquad
                    k_{i-1}+f_i<\tfrac12.                        \tag{11}
\]
Indeed in forward order the continuation through \(v_i\) uses \(F_{v_i}K_{i-2}\), and in reverse order it uses \(F_{v_i}K_{i-1}\). In either order the continuation through the high neighbor has multiplier at least \(B^{1/2}\), so the one through \(v_i\) must fall below that threshold.

Each relevant prime-power level is counted by exactly one of the \(F_v\) or \(K_i\) and has two changes of status in (2). Hence
\[
 (\textstyle\prod F_v\prod K_i)^2\ \mid\ \textstyle\prod_v (R_v R'_v/\gamma_v^2),
 \qquad\qquad \textstyle\sum_i (k_i+f_i)\le 6 .                   \tag{12}
\]
If \(A\) is empty, outside (3) we have \(\sum f_i<3-\log_B Y\). Then (9) and (12) give \(g_*^6<B^{21}/Y\). Take \(Y=\max(1,B^{21}/G^6)\) and use the bound on (3). This contributes at most \(O_\varepsilon(N^\varepsilon(B^{3/2}+B^{45/2}/G^6))\).

If \(A=\{0\}\) after rotation, (11) at vertices 1 and 5 restricts four distinct \(k_i\). The other two satisfy \(k_i\le1\). Therefore
\[
                \sum(k_i+2f_i)<4+2(f_0+f_2+f_3+f_4)\le7+2f_0.
\]
By (9), \(G^6\le B^{19}F_{v_0}^2\). At the vertex \(v_0\) the coupled path weight used for (5) satisfies \(W\ge F_{v_0}^4/B\ge G^{12}/B^{39}\). Its tail contributes \(O_\varepsilon(N^\varepsilon B^{43}/G^{12})\).

For \(A=\{0,2\}\) after rotation and reversal, the three low vertices 1, 3 and 5 give disjoint pairs of \(k_i\) in (11). Consequently \(\sum(k_i+2f_i)<3+2(f_0+f_2+f_4)<8\). For \(A=\{0,3\}\), the four low vertices 1, 2, 4 and 5 give restrictions on all six \(k_i\). Assign \(k_5\) to vertex 1, \(k_1\) to 2, \(k_2\) to 4, \(k_4\) to 5, and \(k_0,k_3\) to 1 and 4 respectively. Summing their six bounds from (11) gives
\[
                    \sum(k_i+2f_i)<3+2(f_0+f_3)+f_2+f_5<8.
\]
If \(A\) has three vertices it is an alternating triple. The three low vertices again give disjoint pairs of \(k_i\) in (11), so \(\sum(k_i+2f_i)<3+2\sum_{i\in A}f_i<8\). In each case with \(|A|\ge2\), (9) gives \(g_*^6<B^{20}\), excluding \(g_*\ge B^{10/3}\) in the remainder. This proves (10).

The original low-gcd sector and both kernel-lcm restrictions remain available in the [residual-cycle criterion](ERDOS_943_SIX_CYCLES.md#exact-remaining-domain-and-both-implications). The complementary singular sector retains the coherent gcd interval below the cutoff in (10), the failures of (8) in both directions, and the original exact common-target conditions. The nonsingular residual is a separate count.
