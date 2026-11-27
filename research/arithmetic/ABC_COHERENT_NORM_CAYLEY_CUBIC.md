# Prime denominators on the Cayley cubic of a coherent norm triangle

Take three distinct odd unramified split primes \(p_1,p_2,p_3\) and [coherent witnesses](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients). Write \(P_i=p_i^3\), \(W=P_1P_2P_3\), and \((v_{12},v_{23},v_{13})=\delta(x,y,z)\) with \(\delta>0\) their common positive gcd. Thus \(x,y,z\) are pairwise coprime, \(E=d\delta^2>0\) is coprime to \(p_1p_2p_3\), and
\[
 A=P_1y^2,\qquad B=P_2z^2,\qquad C=P_3x^2,\qquad t=E/W.
\]
The three normalized trace coordinates, indexed by the vertex opposite their edge, are
\[
 \xi_1=2-tA=\frac{u_{23}^2}{P_2P_3}-2,\quad
 \xi_2=2-tB=\frac{u_{13}^2}{P_1P_3}-2,\quad
 \xi_3=2-tC=\frac{u_{12}^2}{P_1P_2}-2.
\]
On edge \(ij\) this is the field trace of \(\alpha_{ij}/\bar\alpha_{ij}\). The symmetric numerator of the reconstruction formula gives
\[
 2(AB+BC+CA)-(A^2+B^2+C^2)=tABC .
\]
Expanding after the change of coordinates gives the affine Cayley cubic
\[
                   \xi_1^2+\xi_2^2+\xi_3^2-\xi_1\xi_2\xi_3=4 .       \tag{1}
\]

## A Vieta orbit contains at most one coherent prime triangle in trace coordinates

**Proposition.** Starting at any such trace point, every nonempty reduced word in the three Vieta involutions leaves the set of trace points of coherent triangles on distinct odd primes, even with the primes and field free to change. Distinct reduced words give distinct points from the same start.

Here \(s_i\) fixes \(\xi_j,\xi_k\) and replaces \(\xi_i\) by \(\xi_j\xi_k-\xi_i\), with \(\{i,j,k\}=\{1,2,3\}\). These maps exchange the two roots of (1) as a quadratic in the changed coordinate. A word is reduced when consecutive maps differ.

At \(p_i\), the two incident-edge coordinates have valuation exactly \(-3\) and the opposite coordinate \(\xi_i\) is integral. This follows from \(\xi_j=2-tB\) or the corresponding formula: the incident coefficient is a unit at \(p_i\), and \(E\) is a unit there. The opposite coordinate has denominator dividing the product of the other two prime cubes. Start a reduced word with \(s_i\) and keep this one prime \(p_i\) throughout. After \(s_i\) the three valuations are \(-6,-3,-3\), with the strict minimum at the changed coordinate. At every later step the product used in the next involution contains the previous minimum and another coordinate of negative valuation. Its valuation is strictly below all three previous valuations, so it is also the valuation of the new difference. It is the new strict minimum. Consequently all three resulting traces have negative valuation at the same odd prime, which cannot happen for any coherent triangle on distinct primes. If two reduced words gave the same point, composing one with the inverse of the other and cancelling adjacent repetitions would give a nonempty reduced word returning to the start.

On the positive side in \((A,B,C)\), the map changing \(C\) replaces it by \((A-B)^2/C\). It decreases \(C\) when \(C\ge\max(A,B)\), but that move already takes the point out of the prime-denominator locus. The valuation argument applies to the rational points with the prescribed denominators. [Ghosh and Sarnak, Sections 2 and 4](https://arxiv.org/pdf/1706.06712) use the Vieta maps for integral Markoff-type points and describe the integral orbits at level \(4\).

## The odd content of a node congruence

**Proposition.** Define
\[
\begin{array}{ll}
 Q_0=W(2+\xi_1+\xi_2+\xi_3)=u_{12}u_{23}u_{13},&
 Q_1=W(2-\xi_1-\xi_2+\xi_3)=E u_{12}yz,\
 Q_2=W(2+\xi_1-\xi_2-\xi_3)=E u_{23}xz,&
 Q_3=W(2-\xi_1+\xi_2-\xi_3)=-E u_{13}xy.
\end{array}                                                        \tag{2}
\]
All four entries are nonzero, and
\[
    \sum_{i=0}^3 Q_i=8W,\qquad \sum_{i=0}^3\frac1{Q_i}=0,\qquad
    \gcd(Q_1,Q_2,Q_3)_{\rm odd}=E_{\rm odd},\qquad
    \gcd(Q_0,E)_{\rm odd}=1.                                      \tag{3}
\]
The positive odd part of an integer is the quotient of its absolute value by its maximal power of \(2\). In particular for each odd prime power \(\ell^k\mid E\), the primitive projective point reduces to the same coordinate node \([1:0:0:0]\) modulo \(\ell^k\). The gcd of all four entries divides \(8\), and their primitive projective height is between \(W/4\) and \(8W\).

For the last three expressions in (2), substitute \(2-tA,2-tB,2-tC\) and use the reconstruction trace formulas. For the first, multiply \(u_{12}u_{23}=2P_2u_{13}+Exy\) by \(u_{13}\) and use the norm and the formula for \(u_{13}\). Each edge trace \(u_{ij}\) is a unit at either of its odd incident primes, by the norm equation and the units of \(E\) and its normalized coefficient there. In particular no entry of (2) vanishes. The four linear forms in (2) sum to \(8W\). Their third elementary symmetric function is \(8W^3(4-\xi_1^2-\xi_2^2-\xi_3^2+\xi_1\xi_2\xi_3)\), so (1) gives the reciprocal identity in (3).

At an odd prime \(\ell\mid E\), the norms make all three \(u_{ij}\) units. Pairwise coprimality of \(x,y,z\) makes at least one of \(yz,xz,xy\) a unit. Thus the common valuation of the last three entries is exactly \(\nu_\ell(E)\), and \(Q_0\) is a unit. Now suppose \(\ell\nmid E\). If \(\ell\) divides one normalized coefficient, say \(x\), then \(y,z\) are units and \(\ell\nmid p_1p_2\) by the incident prime condition. The edge \(12\) norm makes \(u_{12}\) a unit, so \(Q_1\) is a unit. For \(y\) use \(Q_2\) and for \(z\) use \(Q_3\). If all three normalized coefficients are units, a common prime divisor of \(Q_1,Q_2,Q_3\) would divide all three \(u_{ij}\). This contradicts \(u_{12}u_{23}-Exy=2P_2u_{13}\). This proves the two odd-part statements. The overall content is therefore a power of \(2\) dividing \(8W\) with \(W\) odd. Since \(-2<\xi_i<2\) by the norms and nonvanishing, \(|Q_i|\le8W\) and the sum gives \(\max_i|Q_i|\ge2W\).

There is a further restriction at the prime denominators. At each \(p_i\) the incident coefficients and incident traces are units. The opposite coefficient and opposite trace cannot both be divisible by \(p_i\), by the opposite-edge norm. If either is divisible by \(p_i\), the entries split into two units and two entries divisible by precisely that power of \(p_i\) by (2). If both are units, all four entries are units. More explicitly at \(p_1\) the paired possible common factors are \(y\) on \((Q_1,Q_3)\) and \(u_{23}\) on \((Q_0,Q_2)\), while the two other entries in either case are units.

## Three square variables in the Cayley torsor

Let \(q\) be the primitive integer vector obtained by dividing \(Q\) by its positive content. In the [Cayley torsor of Heath-Brown, Section 2](https://arxiv.org/pdf/math/0210333), the absolute values of the four variables \(y_i\) are the gcds of the three coordinates other than \(q_i\). Thus \(|y_0|_{\rm odd}=E_{\rm odd}\). The other three are squares at every odd prime:
\[
\begin{aligned}
 |y_1|_{\rm odd}&=\gcd(x,u_{23},u_{13})_{\rm odd}^{\,2},\
 |y_2|_{\rm odd}&=\gcd(y,u_{12},u_{13})_{\rm odd}^{\,2},\
 |y_3|_{\rm odd}&=\gcd(z,u_{12},u_{23})_{\rm odd}^{\,2}.
\end{aligned}                                                       \tag{4}
\]
To prove the first identity, it suffices to work at an odd prime \(\ell\) dividing all of \(Q_0,Q_2,Q_3\). Such a prime does not divide \(E\). If it divides \(y\), the edge \(23\) norm makes \(u_{23}\) a unit, so \(Q_2\) is a unit. If it divides \(z\), the edge \(13\) norm makes \(u_{13}\) a unit, so \(Q_3\) is a unit. Thus \(y,z\) are units. If \(x\) were a unit as well, \(Q_2,Q_3\) would force both \(u_{23},u_{13}\) to be nonunits, which contradicts \(u_{12}u_{23}-2P_2u_{13}=Exy\). Therefore \(x\) is a nonunit. Its incident prime condition and its edge norm make \(P_2,u_{12}\) units. The three valuations of \(x,u_{23},u_{13}\) now have their minimum at least twice by the same coherence relation. Equation (2) shows that the common valuation of \(Q_0,Q_2,Q_3\) is the smallest pairwise sum of these three valuations, twice their minimum. Conversely every odd prime dividing \(x,u_{23},u_{13}\) divides all three entries. Division by the content of \(Q\) has no effect at odd primes. Relabel the vertices to obtain the other two identities. In particular each of these three odd square roots divides the corresponding normalized coefficient, so their heights are bounded by the coefficient box rather than by the projective height.

The point also retains the prime direction and the absolute normalized coefficients. Recover the trace coordinates from the projective sums, for example \(\xi_3=4(q_0+q_1)/\sum_iq_i-2\). The reduced denominator of each trace is the product of the two incident prime cubes, by their valuations. Taking pairwise gcds recovers the three ordered primes. Next the quantities \(P_1P_2(2-\xi_3),P_2P_3(2-\xi_1),P_1P_3(2-\xi_2)\) are \(Ex^2,Ey^2,Ez^2\). Their gcd is \(E\), and division by \(E\) recovers the three squares. Thus a count of the projective points with these constraints bounds the number of such normalized signed triangles up to a factor of at most eight.

## Exact dyadic content and a lifted congruence

**Proposition (Full square torsor coordinates).** Put \(g=\gcd(Q_0,Q_1,Q_2,Q_3)>0\) and \(m_i=\gcd(q_j:j\ne i)\) for \(q=Q/g\). An even \(E\) is divisible by \(4\). The content and the first node variable are given by the following cases.

| Condition | \(g\) | \(m_0\) |
| --- | ---: | ---: |
| \(E\) odd and \(x,y,z\) all odd | \(1\) | \(E\) |
| \(E\) odd and one of \(x,y,z\) even | \(2\) | \(E\) |
| \(E\) even | \(8\) | \(E/4\) |

Let \(T_{ij}=u_{ij}\) if \(E\) is odd and \(T_{ij}=u_{ij}/2\) otherwise. Then
\[
 m_1=\gcd(x,T_{23},T_{13})^2,\qquad
 m_2=\gcd(y,T_{12},T_{13})^2,\qquad
 m_3=\gcd(z,T_{12},T_{23})^2.                                      \tag{6}
\]
In particular \(m_0\) has squarefree part \(d\) exactly and the other three \(m_i\) are integer squares.

For the dyadic assertions, pairwise coprimality permits at most one even normalized coefficient. When \(E\) is odd, each edge trace is odd exactly when its normalized coefficient is odd, by (7) of the reconstruction. If all three coefficients are odd, all entries of \(Q\) are odd. Otherwise relabel so that \(x\) is even and let \(v=\nu_2(x)\). In (2) the valuations of \(Q_2,Q_3\) are both \(v\); those of \(Q_0,Q_1\) are both \(\nu_2(u_{12})\). If \(v\ge2\), then \(\nu_2(u_{12})=1\) since \(u_{12}^2=4P_1P_2-Ex^2\). If \(v=1\), both terms on the right divided by \(4\) are odd, so \(u_{12}/2\) is even. In either case \(g=2\) and two entries of \(q\) are odd. Each three-entry gcd is therefore odd.

When \(E\) is even an odd coefficient in (7) of the reconstruction forces its trace to be even and \(4\mid E\). The other traces are also even. Write \(E'=E/4\). Then
\[
 P_1P_2=T_{12}^2+E'x^2,\quad P_2P_3=T_{23}^2+E'y^2,\quad
 P_1P_3=T_{13}^2+E'z^2,\quad T_{12}T_{23}-E'xy=P_2T_{13}.
\]
If \(E'\) is even, all \(T_{ij}\) are odd. The first entry of \(Q/8\) is odd, and the other three have a common valuation exactly \(\nu_2(E')\), since at least one of \(yz,xz,xy\) is odd. This gives \(g=8\) and the claimed gcds at \(2\). If \(E'\) is odd and \(x,y,z\) were all odd, all \(T_{ij}\) would be even, contradicting the last identity modulo \(2\). Relabel so that \(x\) is even. Then \(T_{12}\) is odd, \(T_{23},T_{13}\) are even, and \(Q_1/8=E'T_{12}yz\) is odd. Thus \(g=8\) and \(m_0,m_2,m_3\) are odd. The valuations of \(x,T_{23},T_{13}\) have their minimum at least twice by the last identity and the oddness of \(P_2,T_{12},E',y\). The valuations of \(q_0,q_2,q_3\) are the three pairwise sums of these valuations. Their minimum is twice the smallest, proving (6) at \(2\) in the last case. All odd valuations follow from (3)–(4).

The square roots also give explicit lifts of the original congruences. For example \(b_1=\sqrt{m_1}\) divides \(x\) and \(T_{23}\). The trace formula gives \(A-B=x(P_3x-zu_{23})\), so \(|x|b_1\mid A-B\). At odd primes this lift is exact in the following sense:
\[
              (b_1)_{\rm odd}=\gcd\bigl(x,(A-B)/x\bigr)_{\rm odd}.
\]
Indeed \(z\) is a unit at every prime dividing \(x\), and the two trace formulas give
\(\gcd(x,u_{23})=\gcd(x,u_{13})=\gcd(x,(A-B)/x)\). The other two lifts follow by relabeling. Thus a square node variable records the amount of extra congruence beyond the modulus of the corresponding coefficient, up to its stated dyadic adjustment.

## The six pair factors and their marked prime cubes

Index the three nonzero nodes by edges \(12,23,13\), in that order. Use \(c_i\) for the normalized coefficient \(x,y,z\) of that edge, and \(t_i\) for its trace \(u_{ij}\). Define adjusted coefficients \(C_i\) and traces \(U_i\) as follows.

- If \(E\) is even, take \(C_i=c_i\) and \(U_i=t_i/2\) for all three indices.
- If \(E\) is odd and all \(c_i\) are odd, take \(C_i=c_i\) and \(U_i=t_i\).
- If \(E\) is odd and \(c_i\) is even, divide this \(c_i\) and its own \(t_i\) by \(2\), leaving the other two pairs at their original values.

Put \(b_i=\sqrt{m_i}\) at the three nonzero nodes. For \(\{i,j,k\}=\{1,2,3\}\), define
\[
              a_i=|C_i|/b_i,\qquad e_i=|U_i|/(b_jb_k).             \tag{7}
\]
These are positive integers. The Heath-Brown pair variables, defined as the gcds of the two corresponding coordinates after removing their triple-node factors, are exactly
\[
                              z_{0i}=e_i,\qquad z_{jk}=a_i.         \tag{8}
\]
In particular their pairwise coprimality holds with the specified powers of \(2\). This identifies the remaining six variables in terms of the three normalized coefficients and their coherent traces.

Indeed, the adjusted coordinates obey
\[
                       q_0=U_1U_2U_3,\qquad
                       |q_i|=m_0\,|U_i C_j C_k|.
\]
The divisions in (7) follow from (6) and the pairwise coprimality of the \(b_i\). In the only case where a \(C_i\) is divided by \(2\), all the \(b_i\) are odd, so these divisions preserve integrality. Hence \(|q_0|/(m_1m_2m_3)=e_1e_2e_3\) and \(|q_i|/(m_0m_jm_k)=e_i a_j a_k\).

To verify the pair gcds, it suffices to see that a prime divides at most one of the six residual factors in (7). At an odd prime of \(m_0\) all traces are units, and at most one \(c_i\) is a nonunit. At an odd prime of neither \(m_0\) nor the coefficients, the coherence \(u_{12}u_{23}-2P_2u_{13}=Exy\) excludes two nonunit traces. If a prime outside \(m_0\) divides a coefficient, relabel to \(x\); then \(u_{12}\) is a unit and the valuations of \(x,u_{23},u_{13}\) have their minimum at least twice. Removing \(b_1\) from those three entries leaves at most one nonunit. The same argument at \(2\) uses the adjusted trace relation in the even \(E\) case. When \(E\) is odd with one even coefficient, the dyadic analysis in (5) shows that exactly one of its adjusted coefficient and adjusted trace is odd, while all other entries are odd. This proves the claim and (8).

For each \(i\), let \(p\) be the prime at the vertex opposite edge \(i\). There is an explicit \(\kappa_i\) such that both complementary pair sums have the same prescribed cube factor:
\[
       q_0+q_i=\kappa_i\,p^3\,m_jm_k\,z_{0i}^2,\qquad
       q_j+q_k=\kappa_i\,p^3\,m_0m_i\,z_{jk}^2 .                   \tag{9}
\]
Here \(\kappa_i=1\) for even \(E\), and \(\kappa_i=2\) for odd \(E\) with all coefficients odd. For odd \(E\) with one even coefficient, its edge has \(\kappa_i=4\) and the other two have \(\kappa_i=1\). For example the two sums with edge \(12\) are \(2P_3 u_{12}^2/g\) and \(2E P_3 x^2/g\) by (2). Substitute (5), (7) and (8) to get (9). Thus each prime cube is an explicit normalized pair sum inside the torsor, with the same cube prescribed on both sides of its partition of the four nodes.

## Relation to the rank-one count

The reciprocal model is the projective Cayley surface studied by [Heath-Brown, Theorem](https://arxiv.org/pdf/math/0210333). These points avoid its nine lines. Six require two zero entries. The other three require two pairs of opposite entries, hence \(Q_0+Q_i=0\) for some \(i\). Each such sum is \(2W(2+\xi_j)>0\). The unrestricted count up to height \(B\) on this open subset has order \(B(\log B)^6\), and [Zhao, Theorem 1.1](https://arxiv.org/pdf/2611.14226) gives an asymptotic with Peyre's constant for the same unrestricted count. Here the height is already comparable to \(W\asymp Z^9\) for primes in \([Z,2Z)\), so the unrestricted count does not give a saving on the \(V^3=Z^{15/4}\) record count at the moment crossover. The representation offers an explicit restricted problem on this surface. The sum of the four integral coordinates is prescribed by three prime cubes and, for a field window \(d\in[X,2X)\), the first node variable has squarefree part \(d\) and the remaining three node variables are squares. The [large-lift count](ABC_COHERENT_NORM_LARGE_LIFT_COUNT.md) gives a uniform saving for sufficiently large lifts near the rank-one crossover. The aggregate estimate still needs to control configurations with three small lifts.


## Exact replay on the retained examples

Run the [coherent norm checker](verify_coherent_norm_pencils.py) from the repository root. For every generated triangle it checks (1)–(9) and the lifted congruences with exact rational and integer arithmetic, including the three square torsor variables. On the first triangle of each tested field it checks the common prime denominator after every reduced Vieta word through length six. The proofs above give the uniform claims.
