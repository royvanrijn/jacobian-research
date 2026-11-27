# Products of internal root ratios in two three-edge paths

A six-cycle consists of two three-edge paths with the same opposite endpoints. For each path, the product of the two internal root ratios gives a positive binary norm equation at those endpoints. With the product fixed, a square congruence fixes the endpoint root ratio up to a subpower number of choices. The same norm equation bounds the number of small products through a selected middle edge.

Work in a selected actual graph for \(\kappa_1a^2b^3+\kappa_2c^2d^3=n\), with positive fixed coefficients, squarefree kernels, and coprime summands. The selected left and right roots are bounded by \(R_1,R_2\), respectively, with \(2R_1R_2<n\). Write \(e(G)\) for the number of selected edges. All counts are per target and coefficient pair. The [critical shaped graph](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md) meets the root-product restriction for sufficiently large \(n\).

## The endpoint equation

Take a simple oriented three-edge path from a left kernel to a right kernel. Write the weighted kernel values of its four vertices as \(B_0,B_1,B_2,B_3\), incorporating the coefficient and the cube of the kernel. On edge \(i\), in path order, the roots are \(x_i,y_i\), so
\[
                      B_i x_i^2+B_{i+1} y_i^2=n,\qquad 0\le i\le2.
\]
Define the internal ratios \(s=x_1/y_0\), \(t=x_2/y_1\), their reduced positive product \(J=st=p/q\), and the endpoint roots \(X=x_0\), \(Y=y_2\). Eliminating \(B_1,B_2\) gives
\[
                B_0 J^2 X^2+B_3 Y^2=n(1+J^2-t^2).                 \tag{1}
\]
In fact \(B_1x_1^2=s^2(n-B_0X^2)\). Substitute this on the middle edge and multiply the resulting expression for \(B_2y_1^2\) by \(t^2\) before applying the last edge equation.

Put \(\alpha=\gcd(x_1,y_0)\), \(\beta=\gcd(x_2,y_1)\) and write
\[
       x_1=\alpha u,\qquad y_0=\alpha h w,\qquad
       y_1=\beta v,\qquad x_2=\beta h z,\qquad
       h=\gcd(y_0/\alpha,x_2/\beta).
\]
The reduced numerators and denominators at the two internal vertices are coprime, and \(u,v\) are coprime by the coprimality of the middle edge. After extracting \(h\), the products \(uz,vw\) are therefore coprime. Thus \(p=uz,\ q=vw\) and \(qt=hwz\) is an integer. Positivity in (1) gives the strict bound
\[
                          hwz<\sqrt{p^2+q^2}.                       \tag{2}
\]

Set \(X/Y=\ell/m\) in lowest terms and \(X=\lambda\ell,\ Y=\lambda m\). Reducing (1) after multiplication by \(q^2\) shows that
\[
 C=\frac{B_0 p^2\ell^2+B_3q^2m^2}{n}
\]
is a positive integer. Indeed \(\lambda\) is a unit modulo \(n\) since it divides selected roots, and \(qt\) is an integer. The endpoint equation is exactly
\[
                             (qt)^2+C\lambda^2=p^2+q^2.             \tag{3}
\]
The ratio of the products of initial and terminal roots on the path is \(J\ell/m\). A cycle is singular precisely when this last ratio agrees for the two paths oriented from the same left endpoint to the same right endpoint.

## Recovery and pointwise counts

**Proposition (Fixed endpoint and product recovery).** For two prescribed opposite endpoint kernels and a prescribed \(J\), there are \(O_\varepsilon(n^\varepsilon)\) selected simple three-edge paths, uniformly in the product \(J\). For a prescribed selected middle edge and \(U\ge1\), the number of such paths with \(H(J)\le U\) is \(O_\varepsilon(n^\varepsilon U)\), where \(H(p/q)=\max(p,q)\) in lowest terms.

For the first assertion, \(B_0,B_3,p,q\) are units modulo \(n\) for an actual path. Thus \(n\mid B_0p^2\ell^2+B_3q^2m^2\) fixes \(\ell/m\) in a subpower number of unit square-root classes modulo \(n\). In one class there is at most one reduced fraction with \(\ell\le R_1,\ m\le R_2\), since the cross determinant of two distinct fractions is a nonzero multiple of \(n\) of absolute value at most \(R_1R_2\).

For each such fraction, (3) is a positive binary norm equation. For any positive integer \(C\), solutions \(a,b>0\) to \(a^2+Cb^2=L\) give distinct algebraic integers \(a+b\sqrt{-C}\) of norm \(L\) in the maximal imaginary quadratic order. There are at most \(\tau(L)\) integral ideals of norm \(L\), by the split, inert and ramified prime factorizations, and at most six unit multiples per generator. This bounds the number of pairs in (3) uniformly, since \(p,q\le R_1R_2\le n\) for any path. This is the same norm-counting method used for fixed kernels in [Proposition 3.1 of the squarefull-sum note](https://www.ulam.ai/research/erdos943.pdf). From \(qt,\lambda\) one knows \(X,Y,t\) and hence \(s=J/t\). The root \(X\) and the prescribed starting kernel recover the first edge by subtracting its summand from \(n\) and factoring the complementary summand into a square times a squarefree cube after removing its fixed coefficient. Multiplying \(y_0\) by \(s\) gives \(x_1\), so the same procedure recovers the middle edge. The ratio \(t\) and the middle terminal root recover the last edge. The selected-edge and simple-path tests can only discard candidates.

For the middle-edge assertion, fix \(u\mid x_1,\ v\mid y_1\) at subpower cost. Then \(\alpha,\beta\) are fixed. The integers \(1\le z\le U/u\), \(1\le w\le U/v\) give the two other factors of \(p,q\). By (2) there are at most \(\sqrt2 U/(wz)\) possible positive integers \(h\). Summing over \(w,z\) costs \(O(U(1+\log U)^2)\). Each choice fixes \(y_0,x_2\) and reconstructs the two outer edges from their respective known internal kernels and the target. One may cut \(U\) off at the actual product cap before absorbing logarithms.

**Theorem (Paired internal-product sectors).** For \(U,V\ge1\), the actual selected simple cycles for which one pair of opposite endpoints has path products \(J_1,J_2\) with \(H(J_1)\le U\) and either \(H(J_2)\le V\) or \(H(J_2/J_1)\le V\) number
\[
                       O_\varepsilon(n^\varepsilon e(G)U V^2).     \tag{4}
\]
If the bounded-height quantity on the second path belongs to one prescribed rational squareclass, or is a rational proper power, the bound is \(O_\varepsilon(n^\varepsilon e(G)U V)\). For a fixed relative product \(J_2/J_1\) the bound is \(O_\varepsilon(n^\varepsilon e(G)U)\), uniformly in that product.

Choose the first path through its middle edge using the \(O_\varepsilon(n^\varepsilon U)\) bound. Now both endpoints and \(J_1\) are known. There are \(O(V^2)\) reduced positive fractions of height at most \(V\), and either specification gives \(J_2\). Each \(J_2\) leaves a subpower number of second paths by the fixed-endpoint assertion. A fixed squareclass gives \(O_\varepsilon(n^\varepsilon V)\) fractions by allocating the squarefree part of the numerator times denominator to the two sides and counting their square roots. A proper power gives the same bound by summing \(O(V^{2/j})\) for powers \(j\ge2\), with \(1\) counted once. All heights can be restricted to the actual polynomial caps in \(n\). Retain disjoint internal vertices for a simple cycle.

In the critical shaped graph \(e(G)=O_\varepsilon(n^{3/10+\varepsilon})\), by the coefficient-uniform edge bound used in the [critical criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#5-low-kernel-lcm-and-small-coherent-gcd). Thus (4) meets the cycle budget whenever \(UV^2\le n^{107/510}\), including when both path products have height at most \(n^{107/1530}\). The squareclass or proper-power sector meets it at \(UV\le n^{107/510}\). For equal products it suffices that their common height be at most \(n^{107/510}\). These tests use the product inside each path and the shared endpoint congruence.

## Equal products with prime norm are forbidden

**Proposition (Prime sum of two squares).** Suppose \(R_1R_2<\sqrt n\). If two simple three-edge paths at the same opposite endpoints have the same internal product \(p/q\) and \(p^2+q^2\) is prime, the paths coincide. More generally, for a pair at those endpoints with \(H(J_2/J_1) R_1R_2<\sqrt n\), their initial-to-terminal root product ratios agree, so a simple cycle formed by the pair is singular.

To see the more general statement, each effective ratio \(J_i\ell_i/m_i\) has square congruent to \(-B_3/B_0\) modulo \(n\). Write \(J_2/J_1=r/s\) in lowest terms. The squares of \(s\ell_1m_2\) and \(r\ell_2m_1\) are congruent modulo \(n\) and both are less than \(n\). Equality gives the rational agreement. In particular equal internal products under the stated root cap force the same \(\ell/m\), hence the same \(C\) in (3).

Write \(P=p^2+q^2\) prime. Two positive representations \(a_i^2+C b_i^2=P\) have \(b_i\) units mod \(P\), so \(P\) divides \(a_1b_2-a_2b_1\) or \(a_1b_2+a_2b_1\). Both absolute values are at most \(P/\sqrt C\) by Cauchy--Schwarz on the vectors of length \(\sqrt P\). If \(C>1\), the sum is positive and too small, and the difference must vanish. The norm then fixes the pair. If \(C=1\), the same argument allows only equality or coordinate swap, the latter by equality in Cauchy--Schwarz. As \((p,q)\) already represents \(P\), any solution has \(qt=p\) or \(qt=q\). These imply \(s=J/t=1\) or \(t=1\). In either case two consecutive edges have the same root at their common kernel and hence coincide. Thus the \(C=1\) case supplies no simple path.

## Exact examples and remaining count

At \(n=157\), use the path with kernels \((3,1,1,3)\), alternating left and right, and oriented roots \((2,7),(11,6),(7,2)\). Here \(J=11/6\), \(\ell/m=1\), \(C=27\), and (3) is \(7^2+27\cdot2^2=157\). For this endpoint ratio the prime norm has a single positive pair, \((7,2)\).

At \(n=26497\), two paths from left kernel \(1\) to right kernel \(1\) have internal kernels \((13,3)\) and \((3,13)\) in path order. Their oriented root pairs are respectively
\[
 (82,3),(1,30),(28,73)
 \qquad\text{and}\qquad
 (73,28),(30,1),(3,82).
\]
They give reciprocal products \(14/45,45/14\) and endpoint root ratios \(82/73,73/82\). In both cases the norm is \(42^2+457=2221=p^2+q^2\), while the initial-to-terminal ratios are distinct. These are six actual coprime edges on distinct kernel pairs. This example is for the exact identities in an unshaped graph.

## The remaining cycle criterion

**Corollary (Paired-product residual equivalence).** In the [outward-continuation residual](ERDOS_943_OUTWARD_CONTINUATION.md#remaining-pointwise-estimate), set \(Z_*=n^{107/510}\). For each pair of opposite endpoints and either ordering of the two products \(J_1,J_2\), one may delete cycles with
\[
      H(J_1) H(W)^2\le Z_*,\qquad W\in\{J_2,J_2/J_1\}.
\]
If \(W\) is a rational proper power, one may also delete cycles with \(H(J_1)H(W)\le Z_*\). The latter threshold also works for a fixed squareclass of \(W\), with a subpower allowance for the number of prescribed classes. The pointwise \(O_\varepsilon(n^{26/51+\varepsilon})\) bound for the remaining cycles is equivalent to (H) for the same coefficient family.

In each case dyadically partition the two heights and apply the corresponding count in (4). Each contributing box has the required product bounded by a fixed multiple of \(Z_*\). There are only logarithmically many choices of each height, since the selected root caps bound them by powers of \(n\). The sum of the removed counts meets the cycle budget. Subset inclusion gives the reverse implication. Equal internal products with prime \(p^2+q^2\) are also excluded in the critical root range by the prime obstruction. The open count concerns the surviving actual cycles pointwise in the common target.
