# Squareclasses of paired path products and their norms

A three-edge path with internal product \(J=p/q\) in lowest terms has the positive endpoint norm \(p^2+q^2\). The squareclass of \(pq(p^2+q^2)\) controls both factors jointly. Prescribing this one class leaves a subpower number of possible path products at arbitrary polynomial height. This gives a pointwise cycle sector for Erdős #943, including some cycles with large relative product height.

## A uniform bound in one squareclass

For a positive rational \(J\), put \(\Theta(J)=J(1+J^2)\). Write \([z]\) for the squareclass of a positive rational \(z\), and \(H(p/q)=\max(p,q)\) for coprime positive integers \(p,q\). Since
\[
                  \Theta(p/q)=\frac{pq(p^2+q^2)}{q^4},              \tag{1}
\]
the positive squarefree representative of \([\Theta(p/q)]\) divides \(pq(p^2+q^2)\).

**Proposition (Uniform path-product squareclass bound).** For fixed \(A>0\) and \(\varepsilon>0\), uniformly in the positive squarefree integer \(d\) and \(n\ge2\), there are \(O_{A,\varepsilon}(n^\varepsilon)\) positive rationals \(J\) with \(H(J)\le n^A\) and \([\Theta(J)]=[d]\).

If the set is nonempty then \(d\le2n^{4A}\). Each such \(J\) determines the unique point with \(x=dJ\) and \(y>0\) on
\[
                           E_d:\quad y^2=x^3+d^2x.                 \tag{2}
\]
Indeed \(d^3\Theta(J)\) is a square, and a point with \(x,y>0\) on (2) recovers \(J=x/d\). If \(pq(p^2+q^2)=d k^2\), the positive ordinate is \(y=d^2 k/q^2\). In particular \(h(x)\le\log d+A\log n\), where \(h\) is logarithmic rational height.

For \(y^2=x^3+a x^2+b x\), the two-isogeny descent uses the isogenous curve with coefficients \(-2a,a^2-4b\). The images of the two descent maps are represented by signed squarefree divisors of \(b\) and \(a^2-4b\), respectively, and their cardinalities have product \(2^{r+2}\) for rank \(r\). Here \(a=0,b=d^2\), so the two prime supports lie in the primes dividing \(2d\), and
\[
                  r(E_d)\ll 1+\omega(d)\ll
                           \frac{\Delta}{\log(2+\Delta)},\qquad
                  \Delta=\log(2d).                               \tag{3}
\]
The last inequality follows, for instance, by bounding the product of \(\omega(d)\) distinct primes below by \((\omega(d)+1)!\). The descent input is [Cremona, §3.6, (3.6.2) and the squarefree-divisor test](sources/cremona-modular-elliptic-curves-chapter3.pdf).

The coefficient \(d^2\) is fourth-power-free. [Voutier–Yabuta, Theorems 1.2 and 1.5 and §8](sources/voutier-yabuta-1104.4645v4.pdf) give \(\hat h(P)\ll h(x(P))+\Delta\), \(\hat h(P)\gg\Delta\) for nontorsion points, and at most four rational torsion points. For the uniform lower bound use \(\frac{1}{16}\log(d^2)-\frac18\log 2\) for \(d\ge3\). The two fixed curves \(d=1,2\) have a positive minimum on the nonzero free lattice if their rank is positive.

The canonical height gives a positive definite quadratic form on the real free Mordell–Weil space by [Milne, Chapter IV, Proposition 4.5 and Theorems 4.7, 6.1](sources/milne-elliptic-curves-2006.pdf). In its Euclidean norm, all image points have distance \(O_A(\sqrt{\log n})\) from zero and distinct free classes have mutual distance \(\gg\sqrt{\Delta}\). Disjoint balls of half the latter radius, and the bounded torsion multiplicity, bound the number of points by
\[
                O\!\left((1+C_A\sqrt{\log n/\Delta})^{r(E_d)}\right).
\]
Write \(L=\log n\). If \(\Delta\le\sqrt L\), the logarithm of this bound is \(O_A(\sqrt L\log(2+L))\) by (3). If \(\Delta>\sqrt L\), it is \(O_A(L/\log(2+L))\) by (3), since \(u\log(1+C_A/\sqrt u)\) is bounded for \(0<u=\Delta/L\le O_A(1)\). This proves the subpower bound uniformly in \(d\). The argument only uses upper height and spacing on the rational points. It counts actual positive rational products whether or not they come from a graph.

## Counts for pairs of actual paths

Use the [selected graph and path notation](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-endpoint-equation), with left and right root caps \(R_1,R_2\) satisfying \(2R_1R_2<n\). All counts fix the target \(n\) and coefficient pair. For a path oriented from left to right, its selected middle edge has initial root \(x_1\le R_2\) and terminal root \(y_1\le R_1\). Every internal product of an actual path has height at most \(R_1R_2\).

For a fixed middle edge and \(J=p/q\), the number of paths is
\[
                   O_\varepsilon(n^\varepsilon\, W(J)),
 \qquad W(p/q)=\frac{\min(p,R_2)\min(q,R_1)}{\min(p,q)}.            \tag{4}
\]
Indeed the [factorization for the endpoint equation](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-endpoint-equation) writes \(p=uz,\ q=vw\) with \(u\mid x_1,\ v\mid y_1\). Choosing \(u,v\) costs a divisor factor and determines \(w,z\). The positive integer \(h\) in that factorization is less than \(\sqrt{p^2+q^2}/(wz)=\sqrt{p^2+q^2}\,uv/(pq)\). Use \(u\le\min(p,R_2),v\le\min(q,R_1)\) and reconstruct both outer edges from \(h\) and the known middle edge. This also proves \(W(J)\le \min(H(J),R_1R_2/\min(p,q))\).

**Theorem (Paired path norm squareclasses).** For \(U\ge1\), there are \(O_\varepsilon(n^\varepsilon e(G) U)\) actual simple six-cycles admitting an ordered pair of paths between opposite endpoints with \(H(J_1)\le U\) and a prescribed class of either \(\Theta(J_2)\) or \(\Theta(J_2)/\Theta(J_1)\). The same bound holds with the condition \(W(J_1)\le U\) in place of \(H(J_1)\le U\) if \([\Theta(J_1)]\) is also prescribed. No additional height bound is imposed on \(J_2\) or the relative product.

In the height version, first choose the first path through a selected middle edge at a cost \(O_\varepsilon(n^\varepsilon e(G)U)\) by [fixed middle-edge recovery](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#recovery-and-pointwise-counts). In the version with prescribed \([\Theta(J_1)]\), the uniform squareclass bound gives only a subpower number of \(J_1\), and (4) gives the same cost. For each first path both endpoints and \(J_1\) are known, and either class test fixes \([\Theta(J_2)]\). Now there are subpower many \(J_2\) by the same squareclass bound, and subpower many second paths for each \(J_2\) by fixed-endpoint recovery. Retain pairs that form a simple cycle.

One can sum over a prescribed family of classes, paying linearly for its size. In particular the height version still costs \(O_{B,\varepsilon}(n^\varepsilon e(G)U)\) if the positive squarefree representative of the relative class must divide a positive integer \(M\le n^B\) determined by the first path. There are at most \(\tau(M)\) allowed classes once that path is known. If instead the representative is at most \(Z\), the cost is \(O_\varepsilon(n^\varepsilon e(G)UZ)\).

For the [critical coefficient-shaped graph](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#5-low-kernel-lcm-and-small-coherent-gcd), \(e(G)\ll_\varepsilon n^{3/10+\varepsilon}\), so \(U\le n^{107/510}\) fits the cycle budget \(n^{26/51+\varepsilon}\). In the version with both classes prescribed this includes large reduced products. For example \(R_1R_2=4n^{22/51}/(\kappa_1\kappa_2)\), so the condition on \(W\) holds when \(\min(p_1,q_1)\ge4n^{113/510}/(\kappa_1\kappa_2)\).

## Singular cycles and opposite middle gcds

In a [singular cycle](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#5-low-kernel-lcm-and-small-coherent-gcd), use the primitive positive coherent vertex labels \(h_i\) in cyclic order, and let \(g_i=\gcd(h_i,h_{i+1})\) for edge \(i\), with indices modulo six. The oriented roots on edge \(i\) are \(h_i/g_i,h_{i+1}/g_i\). Choose the paths from vertex \(0\) to vertex \(3\) so that \(0\) is left. The first has ratios \(s_1=g_0/g_1,\ t_1=g_1/g_2\), and the second \(s_2=g_5/g_4,\ t_2=g_4/g_3\).

**Corollary (Prescribed middle gcd ratio).** Singular cycles with an ordered pair of such paths satisfying \(H(J_1)\le U\), a prescribed positive ratio \(g_4/g_1\), and a prescribed ordinary relative squareclass \([J_2/J_1]\) number \(O_\varepsilon(n^\varepsilon e(G)U)\). In particular, equal opposite middle gcds allow an ordinary relative squareclass test at arbitrary relative height.

For either path set \(D_i=s_i^2+t_i^{-2}-1>0\). Positivity follows by dividing the right-hand factor of the endpoint equation by \(t_i^2\). Alternating the six norm equations for the coherent labels cancels the six vertex terms, giving
\[
       g_1^2D_1=g_0^2+g_2^2-g_1^2
               =g_5^2+g_3^2-g_4^2=g_4^2D_2.
\]
Also \(1+J_i^2=t_i^2(1+D_i)\). Therefore
\[
        [\Theta(J_2)/\Theta(J_1)]
           =[J_2/J_1]\left[\frac{1+D_2}{1+D_1}\right].             \tag{5}
\]
Given the first path and the middle gcd ratio, \(D_2\) is fixed. The ordinary relative class thus fixes the combined relative class for the count of second paths. If \(g_1=g_4\), then \(D_1=D_2\) and the ordinary and combined relative classes agree.

## Exact examples and the residual criterion

The [path at \(n=157\)](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#exact-examples-and-remaining-count) has \(J=11/6\), \(pq(p^2+q^2)=10362=2\cdot3\cdot11\cdot157\), and point \((18997,1727^2)\) on \(E_{10362}\). The two paths at \(n=26497\) in the same source have products \(14/45,45/14\). Both give \(pq(p^2+q^2)=9\cdot155470\) with \(155470=2\cdot5\cdot7\cdot2221\) squarefree. In general \(\Theta(1/J)=\Theta(J)/J^4\). The example gives the same combined class on actual nonsingular paths.

**Corollary (Residual after equal norm squareclasses).** From the [paired-product residual](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-remaining-cycle-criterion), delete cycles with an ordering of the paths at any opposite endpoints such that \(H(J_1)\le n^{107/510}\) and \([\Theta(J_2)/\Theta(J_1)]=[1]\). The pointwise \(O_\varepsilon(n^{26/51+\varepsilon})\) bound for this smaller residual is equivalent to (H) for the same coefficient family.

Apply the paired path count at \(U=n^{107/510}\) to the deleted cycles. Conversely, a bound for (H) bounds any subset. The same counting argument allows any specified subpower family of relative classes, and allows first paths of large height when both classes are specified and (4) meets the budget. The unrestricted joint distribution of the classes on actual pairs of paths is open.
