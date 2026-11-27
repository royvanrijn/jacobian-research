# Joint continuation bounds in the critical six-cycle rectangle

Use the target, coefficient family and selected actual graph of the [critical rectangular criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md). Write \(A_*=n^{11/51}\), so selected roots are at most \(2A_*\) and the kernel cap is \(2n^{32/153}\). Put \(Y=n^{107/510}\) and \(Q=n^{14/153}\). The results use the coefficient-uniform three-tenths edge bound and the actual codegree bound in that criterion.

In an oriented simple six-cycle, index vertices \(v_i\) and edges \(e_i=v_iv_{i+1}\) modulo six. Let \(R,R'\) be the roots at \(v_i\) on \(e_{i-1},e_i\), and let \(D_i=|R-R'|/\gcd(R,R')>0\). Write \(\alpha_i\) for the gcd of the two opposite roots on this path with shared vertex \(v_i\), and \(d_i=\min(D_i,A_*)\). On a singular cycle use the divisors \(F_{v_i},K_i\) and \(g_*\) of the [rectangular continuation proof](ERDOS_943_CRITICAL_RECTANGULAR_CONTINUATION.md). In the chosen orientation put \(m_i=F_{v_i}K_{i-2}\), and define
\[
 H_i=\begin{cases}
       \operatorname{lcm}(\alpha_i,m_i),&\text{singular},\
       \alpha_i,&\text{nonsingular}.
      \end{cases}
 \qquad
 L_i=\max(1,d_i/H_i^2),\qquad S_i=\max(1,d_i/\alpha_i^2).
\]
Reversing the orientation replaces the continuation divisor at \(v_i\) by \(F_{v_i}K_{i-1}\), in the original indexing. The other quantities are defined in the same way in either orientation.

## Counting with the reduced difference and the modulus together

**Theorem (Joint continuation sectors).** Each of the following sets has \(O_\varepsilon(n^{26/51+\varepsilon})\) cycles, uniformly for the target and coefficients.

- Cycles with \(L_iL_{i+1}\le Q\) for some two consecutive steps in either orientation.
- Cycles with \(S_iL_{i+1}L_{i+2}\le Y\) for some three consecutive steps in either orientation.
- Singular cycles with \(g_*\ge16n^{2539/4080}\).

For the first two assertions, suppose first that only the entering edge at \(v_i\) is known. There are divisor-many possibilities for \(\alpha_i\) from the preceding opposite root. For a specified \(\alpha_i\), the next root at \(v_i\) satisfies a unit square congruence modulo \(\alpha_i^2\) because the next opposite root is divisible by \(\alpha_i\). In particular \(\alpha_i\) is coprime to the next root. If two entering edges are known on a singular cycle, the continuation congruence in the previous proof supplies subpower many classes modulo \(m_i^2\), for divisor-many possible \(F_{v_i},K_{i-2}\) recovered from the entering edges. The divisibility in equation (3) there places \(m_i\) in the next opposite root. Combining the residue lists prime by prime gives subpower many classes modulo \(H_i^2\), and \(H_i\) divides the next opposite root. Thus the next joining root is a unit modulo \(H_i\). On a nonsingular cycle use the single congruence modulo \(\alpha_i^2\).

To count from a residue list modulo \(h^2\) for either \(h=\alpha_i\) or \(h=H_i\), fix \(\gamma=\gcd(R,R')\) by a divisor choice from the entering root \(R\). Put \(s=R/\gamma\), so \(R'=\gamma(s\pm D_i)\). Since \((R',h)=1\), \((\gamma,h)=1\). Consequently the list for \(R'\) gives a subpower list of classes for the positive integer \(D_i\) modulo \(h^2\) for either sign. If \(d_i\) lies in a specified dyadic interval of size \(D_0\), then \(D_i\le4D_0\) since \(D_i\le2A_*\). The number of candidates is bounded by a subpower multiple of \(1+D_0/h^2\). Each actual root determines the complementary summand, its squarefree kernel and its root. This gives the step costs \(n^\varepsilon S_i\) with one entering edge and \(n^\varepsilon L_i\) with two, up to absolute and \(\varepsilon\)-dependent factors. Fixing all dyadic intervals before counting keeps the product of the step bounds within a constant factor of the stated product, and the number of intervals contributes a subpower factor.

For the two-step assertion, choose two entering edges in \(O_\varepsilon(n^{64/153+\varepsilon})\) ways by picking their two outer kernels and an actual common neighbor. For the three-step assertion choose a single entering edge, using the \(O_\varepsilon(n^{3/10+\varepsilon})\) edge bound. The resulting four-edge path in either case closes by the actual codegree bound on the endpoints. The exponents \(64/153+14/153\) and \(3/10+107/510\) are both \(26/51\).

This reconstruction uses coprimality to pass from a root residue to a residue of its reduced difference. It applies in the selected graph before graph supersaturation is invoked.

## The strengthened coherent gcd cutoff

Remove the first two counted sets and the already counted sector with any \(\alpha_i\ge n^{43/306}\), using Theorem 5 of the critical rectangular criterion. In a singular cycle of the remaining set one then has \(F_{v_i}\le\alpha_i<n^{43/306}\). Write \(T=n^{1/153}\), \(f_i=\log_T F_{v_i}\), \(w_i=\log_T m_i\) in a fixed orientation, and
\[
                      u_i=\min(f_i,16.5),\qquad z_i=\min(w_i,16.5),
\]
To avoid dependence on \(D_i\) in the cutoff calculation, upper bound the actual step costs by
\(\max(1,A_*/F_{v_i}^2)=T^{33-2u_i}\) and \(\max(1,A_*/m_i^2)=T^{33-2z_i}\), respectively. The failures of the two product tests give
\[
        z_i+z_{i+1}<26,\qquad
        u_i+z_{i+1}+z_{i+2}<33.45,\qquad 0\le u_i\le z_i.           \tag{1}
\]
Also \(m_i\le2T^{33}\) by its divisibility into one opposite root, and \(f_i<21.5\). Denote by \(k\) the number of \(w_i\ge16.5\), and by \(s\) the number of \(f_i\ge16.5\). The \(k\) high positions are independent in the cycle by the pair inequalities, and \(s\le k\). The product of the six edge gcds is bounded by
\[
            \prod_e g_e\le (2T^{33})^{12}\prod_i m_i F_{v_i}.        \tag{2}
\]
We show that the last product is at most \(2^3 T^{175.275}\).

If \(k\le2\) and \(s\le1\), let \(V_0=\sum_i z_i\), \(U_0=\sum_i u_i\). Summing (1) gives \(U_0\le\min(V_0,200.7-2V_0)\), hence \(U_0+V_0\le133.8\). Undoing the clips increases the exponent of the product by at most \(2(16.5)+5=38\), apart from the factor \(2^k\) from the height bounds for \(m_i\). This gives exponent at most \(171.8\).

If \(k=3\), the high positions alternate. The three low \(z_i\) are at most \(9.5\). For a high position \(i\), (1) gives \(u_i\le16.95-z_{i+1}\). For a low position \(i\), it gives \(u_i\le\min(z_i,16.95-z_{i+2})\). Summing the high bounds cancels the sum of the low \(z_i\) in \(U_0+V_0\). Undoing the clips shows that the product in (2), after removal of the factor \(2^k\), has exponent at most
\[
                           149.85+5s+\sum_{i\ {\rm low}}u_i.       \tag{3}
\]
Each high \(f_i\) forces the next low \(z_{i+1}\le0.45\). For \(s=0\), sum the three inequalities \(\min(z_i,16.95-z_{i+2})\le(z_i+16.95-z_{i+2})/2\) to obtain a low \(u_i\) sum at most \(25.425\). For \(s=1\), one low \(z_i\) is at most \(0.45\), and the two other \(u_i\) sum to at most \(16.95\) by their consecutive low-pair bounds. For \(s=2\) or \(3\), use \(u_i\le z_i\) to bound the sum by \(10.4\) or \(1.35\). The four bounds from (3) are \(175.275,172.25,170.25,166.2\), all at most \(175.275\).

It remains to have \(k=s=2\). If \(v_i\) neighbors a vertex with high \(F_{v_j}\), the two-step test in both orientations gives
\[
                        F_{v_i}K_{i-2},\ F_{v_i}K_{i-1}
                                              < T^{9.5}.          \tag{4}
\]
In fact the cost upper bound at its high neighbor is \(1\), so the upper bound at \(i\) must exceed \(Q=T^{14}\). If the two high vertices are \(0,2\), the three low vertices \(1,3,5\) give disjoint pairs of \(K_i\) in (4). Use also the three-step test starting at \(2\) forward and the one starting at \(0\) in reverse. Since \(u_2=u_0=16.5\) and the following two steps on each path include \(4\), the clipped exponents for the forward step at \(3\) and reverse step at \(5\) are at most \(16.95-f_4\). If this number is less than \(9.5\) it is below the clipping threshold and gives the same bound on the unclipped exponents. Thus for \(K_1\) at vertex \(3\) and \(K_4\) at vertex \(5\), we can replace \(9.5\) in (4) by \(\min(9.5,16.95-f_4)\). The exponent of \((\prod_i K_i)(\prod_i F_{v_i}^2)\) is at most
\[
             4\cdot9.5+4\cdot21.5+
                        2 f_4+2\min(9.5,16.95-f_4)
                                                    \le157.9.
\]
If the highs are \(0,3\), assign \(K_5,K_0\) to vertex 1, \(K_1\) to 2, \(K_2,K_3\) to 4 and \(K_4\) to 5 in (4). Vertices 2 and 5 each retain one factor \(F_{v_i}\le T^{9.5}\). This gives exponent \(8\cdot9.5+4\cdot21.5=162\). These are the two cases up to rotation and reversal.

Equation (2) now gives \(g_*^6\le 2^{15} T^{571.275}\) outside the counted sectors. Since \(571.275/(6\cdot153)=2539/4080\) and \(16^6>2^{15}\), the third assertion follows.

## The remaining cycle criterion

**Corollary (Joint continuation residual equivalence).** In the [previous residual](ERDOS_943_CRITICAL_RECTANGULAR_CONTINUATION.md#the-remaining-critical-rectangular-count) \(E'_\kappa(n)\), delete the three sets in the theorem. The bound \(O_\varepsilon(n^{26/51+\varepsilon})\) for the resulting subset is equivalent to (H), for the same shaped graph, targets and coefficient gate. The remaining cycles fail both oriented product tests, and the singular ones have \(g_*<16n^{2539/4080}\). The conditional \(5/17\) representation conclusion still requires this pointwise residual count.

One direction is subset inclusion, and the other follows by adding the three uniform sector counts and applying the previous equivalence.

The [matching-ratio criterion](ERDOS_943_MATCHING_RATIO_RECTANGULAR.md) further restricts this residual by testing either matching in singular and nonsingular cycles.
