# Continuation in the critical rectangular six-cycle graph

Use the coefficient-shaped graph (2s) and target and coefficient family of [the critical rectangular criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md). Write
\[
 A=2n^{11/51},\qquad B=2n^{32/153},\qquad
 M=n^{13/153},\qquad Y=n^{107/510},\qquad F_0=n^{43/306}.
\]
The selected roots are at most \(A\) and the kernels at most \(B\). All counts below are per target and coefficient pair, with constants uniform in both. We use the actual codegree bound \(O_\varepsilon(n^\varepsilon)\) in Lemma 4 of the criterion and the coefficient-uniform three-tenths edge bound used in its Lemma 7.

For a simple six-cycle write \(e_i=v_iv_{i+1}\) modulo six. At each shared vertex let \(D=|R-R'|/\gcd(R,R')\) for its two selected roots. Denote by \(D_*\) the minimum product of three such numbers across four consecutive edges. On a singular cycle use the positive coherent labels \(h_i\) of global gcd one from Lemma 7 of the criterion. Write \(g_e=\gcd(h_i,h_j)\) for \(e=v_iv_j\) and \(g_*=\min_e g_e\). Let \(F_v\) be the gcd of the four edge gcds off \(v\), and use the three-edge divisors from the [balanced argument](ERDOS_943_LARGE_COHERENT_GCD_CYCLES.md#congruences-from-three-complementary-edges):
\[
 J_i=\gcd(g_{e_{i+3}},g_{e_{i+4}},g_{e_{i+5}}),\qquad
 K_i=J_i/\gcd(J_i,h_{i+1}h_{i+2}).
\]

## Coefficient-uniform sectors

**Theorem (Rectangular continuation sectors).** Each of the following sectors has \(O_\varepsilon(n^{26/51+\varepsilon})\) cycles.

1. All cycles with \(D_*\le Y\), including nonsingular cycles.
2. Singular cycles with \(F_{v_{i+2}}K_i\ge M\) and \(F_{v_{i+3}}K_{i+1}\ge M\) for some \(i\) in either orientation.
3. Singular cycles with \(g_*\ge 16n^{6059/9180}\).

Furthermore, every singular cycle outside sector 2 satisfies
\[
                g_*< 16n^{89/153}\left(\prod_v F_v\right)^{1/6}.   \tag{1}
\]
In particular the smaller cutoff \(16n^{89/153}\) suffices when all \(F_v=1\).

For sector 1, each \(D\) is positive because equal roots at a shared kernel fix the complementary summand and hence its squarefree kernel. The [reduced-difference argument](ERDOS_943_REDUCED_ROOT_DIFFERENCES.md) applies with fixed coefficients: for a known root \(R\) and prescribed \(D\), write the next root as \(R'=rR/s\) with \((r,s)=1\). Then \(s\mid R\) and \(r=s\pm D\). After fixing \(R'\), subtract its summand from \(n\), divide by the complementary coefficient and recover the unique positive square-times-squarefree-cube decomposition if it exists. There are divisor-many choices per step. There are \(O(Y(1+\log Y)^2)\) ordered triples of positive integers with product at most \(Y\). Seed with one edge, take three steps and use the actual codegree to close between the two endpoints of the four-edge path. The edge bound gives the exponent \(3/10+107/510=26/51\).

For sector 2, the proof of the [three-edge continuation congruence](ERDOS_943_LARGE_COHERENT_GCD_CYCLES.md#congruences-from-three-complementary-edges) extends with root cap \(A\). The coherent relation between the three path edge gcds uses only the selected roots. The alternating sum of the squares of the six edge gcds is still zero: after multiplying each edge equation by its gcd squared, the vertex potential is \(\kappa(v)b_v^3 h_v^2\), where \(b_v\) is its kernel and \(\kappa(v)\) the fixed coefficient on its side. Each potential cancels exactly in the alternating sum. Thus the congruence modulo \(J_i^2\) and its prime-power reductions on \(K_i\) are identical. At primes of \(f=F_{v_{i+2}}\) outside \(K_i\), the condition on the joining root \(R\) is \(\kappa(v_{i+2})b_{v_{i+2}}^3R^2\equiv n\bmod f^2\) locally. Both the coefficient times the kernel cube and \(n\) are units at these primes by coprimality of the selected summands, since the other root is divisible by \(f\). The first two edges contain roots divisible by \(K_i\) and \(f\), so these divisors also have subpower many branches. The resulting continuation count, including exact recovery of the complementary summand, is \(O_\varepsilon(n^\varepsilon(1+A/(fK_i)^2))\) for prescribed divisors. There is no use of a lower root cutoff.

Two edges with a common vertex cost \(O_\varepsilon(n^\varepsilon B^2)\) choices by choosing the outer kernels and then their actual common neighbor. Each of the next two steps costs \(O_\varepsilon(n^{7/153+\varepsilon})\), since \(A/M^2=2n^{7/153}\). Close by the actual codegree. The resulting exponent is \(64/153+14/153=26/51\).

Inequality (9) of the balanced argument depends only on the primitive coherent labels and their exact gcds. Here it gives
\[
 \prod_e g_e\le A^{12}\prod_i K_i\prod_v F_v^2 .                    \tag{2}
\]
There is a joint height restriction on the two factors of each step:
\[
             F_{v_{i+2}}K_i\mid R_{v_{i+3},e_{i+2}},\qquad
             F_{v_{i+1}}K_i\mid R_{v_i,e_i},                       \tag{3}
\]
where \(R_{v,e}\) is the selected root at \(v\) on \(e\). To see the first divisibility at a prime \(\ell\mid K_i\), put \(a=\nu_\ell(J_i)\), \(t_1=\nu_\ell(h_{i+1})\), \(t_2=\nu_\ell(h_{i+2})\), and \(k=\nu_\ell(K_i)>0\). The four complementary labels have valuation at least \(a=t_1+t_2+k\). The four-edge divisor \(F_{v_{i+2}}\) has valuation \(t_1\), since its supporting path covers the four complementary vertices and \(v_{i+1}\). The outer root at \(v_{i+3}\) has valuation at least \(a-t_2=t_1+k\). At primes outside \(K_i\) use the divisibility by \(F_{v_{i+2}}\) alone. Reversing the path proves the second assertion.

Outside sector 2 at most three of the six forward multipliers \(F_{v_{i+2}}K_i\) reach \(M\). All are bounded by \(A\) by (3), and \(M\le A\). Hence \((\prod_i K_i)(\prod_v F_v)\le A^3M^3\), and (2) gives \(g_*^6\le A^{15}M^3\prod_v F_v\). Taking sixth roots proves (1), using \(2^{15/6}<16\).

## The large coherent gcd bound

Remove sector 2. Theorem 5 of the [critical rectangular criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#4-actual-rectangular-codegrees-and-opposite-root-gcd-frequency) also allows deletion of all cycles with any \(F_v\ge F_0\), since \(F_v\) divides both opposite roots of the path through \(v\). Put \(L_v=\max(1,A/F_v^2)\). Cycles with \(L_{v_i}L_{v_{i+1}}L_{v_{i+2}}\le Y\) for any three consecutive vertices already meet the budget. Indeed, seed with one edge and choose the first divisor \(F_v\) from the preceding root. At each step an actual joining root satisfies the same unit square congruence modulo \(F_v^2\) used above. The number of residue classes is at most \(4^{\omega(F_v)}\), and a root fixes the complementary summand. Choose the next divisor from the newly found edge, then close after three steps by the codegree bound. In each dyadic box of the three divisors the product of the step bounds is \(O(n^\varepsilon Y)\) whenever the box contains a cycle in this sector. Summing the logarithmically many boxes and using the three-tenths edge bound gives \(O_\varepsilon(n^{26/51+\varepsilon})\).

In the remaining cycles let \(S=\{i:F_{v_i}\ge M\}\). It has no consecutive indices by sector 2. Multiplication of the six inequalities \(L_{v_i}L_{v_{i+1}}L_{v_{i+2}}>Y\) gives
\[
               \prod_v\min(F_v,\sqrt A)<A^3/Y.
\]
When \(|S|\le1\) at most one \(F_v\) can exceed \(\sqrt A\), since \(M<\sqrt A\). With \(F_v<F_0\) and \(\max(1,F_0/\sqrt A)\le n^{5/153}\), this gives \(\prod_vF_v< (A^3/Y)n^{5/153}\). Combining with the product bound from (3) gives
\[
                    g_*^6< A^{18}M^3 n^{5/153}/Y
                         = 2^{18} n^{6059/1530}.
\]

If \(F_{v_i}<M\) and \(v_i\) has a neighbor indexed by \(S\), the failure of sector 2 in both orientations gives
\[
                      K_{i-2}<M/F_{v_i},\qquad K_{i-1}<M/F_{v_i}.   \tag{4}
\]
In each orientation the step at the high neighbor already has multiplier at least \(M\). The forward step at \(v_i\) uses \(K_{i-2}\), and the reverse step uses \(K_{i-1}\). Use (4) in the second bound of (2).

- For \(S=\{0,2\}\), the pairs at 1, 3 and 5 cover all six \(K_i\). This gives \((\prod K_i)(\prod F_v^2)\le M^8F_0^4\). For \(S=\{0,3\}\), use (4) at 1 for \(K_5,K_0\), at 2 for \(K_1\), at 4 for \(K_2,K_3\), and at 5 for \(K_4\). This leaves one factor \(F_{v_2}F_{v_5}\) bounded by \(M^2\), giving the same estimate.
- For an alternating triple \(S\), take the two bounds at each of the three low vertices. This gives \((\prod K_i)(\prod F_v^2)\le M^6F_0^6\).

These estimates give powers of \(n\) in \(g_*^6\) of \(586/153\) and \(603/153\), respectively, both less than \(6059/1530\), with constants at most \(2^{18}\). Since \(16^6>2^{18}\), this proves sector 3.

## The remaining critical rectangular count

**Corollary (Continuation residual equivalence).** In the [previous residual](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#7-both-directions-of-the-remaining-rectangular-cycle-problem) \(E_\kappa(n)\), keep only cycles with \(D_*>Y\). For singular cycles also require \(g_*<16n^{6059/9180}\) and failure of the two successive multiplier thresholds in both orientations. Call this subset \(E'_\kappa(n)\). Hypothesis (H) of the criterion is equivalent to the \(O_\varepsilon(n^{26/51+\varepsilon})\) count on \(E'_\kappa(n)\) in the same coefficient family. In the singular branch the extra inequality (1) also holds.

The forward implication is subset inclusion. For the converse, the theorem counts every discarded sector pointwise in the target and coefficients. Proposition 10 of the criterion then yields (H). Both the surviving pointwise cycle count and the resulting \(5/17\) representation bound remain open.

The [joint continuation bounds](ERDOS_943_JOINT_RECTANGULAR_CONTINUATION.md) count reduced differences inside the continuation residue classes and give a narrower equivalent residual.
