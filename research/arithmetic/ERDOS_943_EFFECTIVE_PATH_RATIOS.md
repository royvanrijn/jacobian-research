# Effective path ratios and opposite kernels in six-cycles

A three-edge path has an effective ratio, the product of its initial roots divided by the product of its terminal roots. Its residue class links the kernels at the two endpoints. Specifying this ratio at one opposite pair in a six-cycle recovers one kernel from the other. At two opposite pairs it leaves only actual codegrees after an edge is selected. The same path has a primitive integer label vector whose endpoint gcd gives a separate way to count paired paths.

Use the [selected coprime graph](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-endpoint-equation) at target \(n\) with fixed positive coefficients \(\kappa_1,\kappa_2\). The left and right kernel caps are \(K_1,K_2\), and the left and right root caps are \(R_1,R_2\). All caps are at most \(n\). Where codegrees are used, assume the actual graph has subpower same-side codegrees and \(K_1K_2<n\). These hold in the [critical coefficient rectangle](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#4-actual-rectangular-codegrees-and-opposite-root-gcd-frequency) for sufficiently large \(n\). All counts fix the target and coefficient pair. The underlying fixed-kernel norm count of [Chojecki, Proposition 3.1](https://www.ulam.ai/research/erdos943.pdf) bounds representations along one kernel pair. Here the effective slope relates the nonadjacent kernels of an actual three-edge path.

## Primitive labels on a three-edge path

Write the path from left to right with weighted kernel cubes \(B_0,B_1,B_2,B_3\) and selected roots \((x_i,y_i)\) on edge \(i\), \(0\le i\le2\), in that order. Thus \(B_i x_i^2+B_{i+1}y_i^2=n\). Put \(X=x_0\), \(Y=y_2\), \(s=x_1/y_0\), \(t=x_2/y_1\) and \(J=st=p/q\) in lowest terms. Use the [path factorization](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-endpoint-equation)
\[
  x_1=\alpha u,\quad y_0=\alpha h w,\quad
  y_1=\beta v,\quad x_2=\beta h z,\qquad
  p=uz,\quad q=vw,\quad A=qt=h w z .
\]
The effective ratio from left to right is \(\rho=J X/Y\). It agrees for the two paths at opposite endpoints exactly when the cycle is singular. Define the endpoint gcd \(G=\gcd(pX,qY)\), and write \(U=pX/G\), \(V=qY/G\), so \(\rho=U/V\).

**Proposition (Primitive path labels and endpoint reconstruction).** The unique primitive positive integer label vector with successive ratios \(x_i/y_i\) is
\[
        (pX,\ A x_1,\ A y_1,\ qY).
\]
In particular \(\gcd(A,G)=1\). If \(C=(B_0U^2+B_3V^2)/n\), then \(C\) is a positive integer, and
\[
         A^2+C G^2=p^2+q^2,\qquad
         U G\le R_1^2R_2,\quad V G\le R_1R_2^2.                    \tag{1}
\]
For fixed endpoint kernels, \(G\) and \(\rho\) leave \(O_\varepsilon(n^\varepsilon)\) paths. For fixed endpoint kernels and \(G\) alone the count is
\[
              O_\varepsilon\left(n^\varepsilon
                  \left(1+\frac{(R_1R_2)^3}{n G^2}\right)\right).  \tag{2}
\]

The successive ratios follow from \(p y_0=A x_1\) and \(q x_2=A y_1\). A prime dividing \(A\) cannot divide both endpoint labels. If it divides \(w\), it divides \(q\) and is prime to \(p\) and \(X\). If it divides \(z\), it divides \(p\) and is prime to \(q\) and \(Y\). If it divides \(h\), it is prime to both \(X,Y\) and at least one of \(p,q\). Here the outer-edge coprimalities give exactly the units at \(X\) and \(Y\). This proves \(\gcd(A,G)=1\), and the middle-edge coprimality gives \(\gcd(Ax_1,Ay_1)=A\), so the vector is primitive.

The endpoint equation for the internal product is
\[
              B_0(pX)^2+B_3(qY)^2=n(p^2+q^2-A^2).
\]
All path roots, their divisors and endpoint kernels are units modulo \(n\). Dividing the congruence by \(G^2\) gives
\[
                         B_0U^2+B_3V^2\equiv0\pmod n ,            \tag{3}
\]
so \(C\) is a positive integer and (1) follows. For given \(U,V,G\), choose \(p\mid U G,\ q\mid V G\) at subpower cost. This fixes \(X,Y\), so the outer edges are recovered by subtracting the chosen summand and factoring the complementary summand. Equation (1) then fixes the positive \(A\), and the two formulas for \(p y_0\) and \(q x_2\) recover both middle roots. Invalid divisibility, edge or simplicity conditions discard the candidate.

If only \(G\) is fixed, (3) gives a subpower number of unit square-root classes for \(U/V\). In each such class, distinct primitive pairs in the rectangle of (1) have cross determinants nonzero and divisible by \(n\). Scale by the two rectangle widths. Consecutive directions then have angular separation at least a constant times the reciprocal of the rectangle area divided by \(n\). Equivalently the number of directions is \(O(1+\text{area}/n)\), by the [rectangular fraction-packing argument](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#an-unweighted-ratio-with-no-coefficient-loss). The area bound is \((R_1R_2)^3/G^2\). Applying the preceding recovery for each direction proves (2).

## Starting a path at its outer edge

**Proposition (Outer-edge recovery at arbitrary internal product height).** Given one selected outer edge of a three-edge path and its internal product \(J\), there are \(O_\varepsilon(n^\varepsilon)\) such paths, with no further height restriction on \(J\).

Start with the first edge and its known terminal root \(y_0\). For \(p/q=J\), the choices of \(u\mid p,\ v\mid q\) fix \(z=p/u,\ w=q/v\). If \(w\mid y_0\), choose \(h\mid y_0/w\), which fixes \(\alpha=y_0/(hw)\). The root \(x_1=\alpha u\) at the known middle-edge kernel recovers the middle edge by the exact target equation. Now \(\beta=y_1/v\) is forced, and \(x_2=\beta h z\) recovers the last edge. Each choice has polynomial-size divisor inputs since \(p,q\le R_1R_2\). Reverse the path for an anchor at its last edge, replacing \(J\) by \(1/J\).

Write \(\Theta(z)=z(1+z^2)\) and \([z]\) for the positive rational squareclass. By the [uniform path-product squareclass bound](ERDOS_943_PAIRED_PATH_NORM_CLASSES.md#a-uniform-bound-in-one-squareclass), a prescribed \([\Theta(J)]\) has only a subpower number of possible \(J\), since every actual product has polynomial height. Hence one prescribed absolute class for \(J_1\), and one prescribed absolute class for \(J_2\) or relative class \([\Theta(J_2)/\Theta(J_1)]\), give \(O_\varepsilon(n^\varepsilon e(G))\) cycles. Seed the first path at an outer edge and use fixed endpoints and \(J_2\) to recover the second path. The two prescribed classes require no cutoff on either product or on the weight in the earlier middle-edge estimate.

## Cycle counts from one or two effective slopes

For fixed \(\rho=U/V\), equation (3) gives at most subpower many possible opposite kernels once one endpoint kernel is known. In fact \(B_0=\kappa_1 b^3,\ B_3=\kappa_2 d^3\), and the equation fixes \(b/d\) to a subpower number of unit cubic-root classes. Each root class fixes one residue for the missing \(b\) or \(d\), whose cap is at most \(n\). If neither kernel is known but their gcd \(g\) is given, use \((b/g,d/g)=1\). A cubic-root class now contains at most one such reduced pair, since the cross difference is bounded in absolute value by \(K_1K_2/g^2<n\).

**Theorem (One or two effective slopes).** Let \(K=\min(K_1,K_2)\). Cycles with one prescribed effective slope on either path at an opposite pair number \(O_\varepsilon(n^\varepsilon e(G)K)\). If that opposite pair has kernel gcd at most \(L\), the bound improves to
\[
                          O_\varepsilon(n^\varepsilon
                                        \min(e(G)K,LK^2)).       \tag{4}
\]
For two prescribed effective slopes at distinct opposite pairs the bound is \(O_\varepsilon(n^\varepsilon e(G))\). If the two associated kernel gcds are at most \(L_1,L_2\), it is
\[
                          O_\varepsilon(n^\varepsilon
                                        \min(e(G),L_1L_2)).      \tag{5}
\]
These statements hold with prescribed classes \([\Theta(\rho)]\) in place of prescribed slopes. With two classes, one absolute class and the relative class suffice. All estimates allow either determinant branch and arbitrary effective slope height.

For (4), label the relevant path \(v_0v_1v_2v_3\) with \(v_0\) left. Seed \(v_0v_1\) by a selected edge. The slope gives \(v_3\), and \(v_2\) is an actual common neighbor of \(v_1,v_3\). To find the complementary path, choose the left neighbor \(v_4\) of \(v_3\) at cost at most \(K_1\). Then \(v_5\) is an actual common neighbor of \(v_4,v_0\). Seeding with the last edge of the first path interchanges the two kernel caps in this count. If the kernel gcd of \(v_0,v_3\) is at most \(L\), choose it and recover both endpoints from the cube congruence. At fixed endpoints each path is recovered by choosing its first neighbor on the side of smaller kernel cap and then using one codegree. This costs at most \(K n^\varepsilon\) per path and gives the second term in (4).

For (5), there is an edge between the two prescribed opposite pairs, and its opposite edge joins the same two pairs. Seed with one such edge. Both other endpoints of the pairs follow from (3). Each of the two remaining vertices is a common neighbor of two known distinct vertices. If both kernel gcds are bounded, choose them first at a cost \(L_1L_2\), then recover the four kernels using the two cubic-root tests and use the same two codegrees. The root choices already describe actual selected edges. Close only simple cycles. Finally \(H(\rho)\le\max(R_1^2R_2,R_1R_2^2)\) by (1), so replacing prescribed slopes by prescribed classes costs a subpower factor by the same uniform squareclass bound.

In the critical rectangle \(K\le 2n^{32/153}\) and \(e(G)\ll_\varepsilon n^{3/10+\varepsilon}\), so one prescribed class costs \(O_\varepsilon(n^{779/1530+\varepsilon})\), below the budget \(n^{26/51+\varepsilon}\). With a kernel gcd cap \(L\) one can sum over a family of \(Z\) classes whenever \(LZ\le n^{14/153}\) by (4). At two opposite pairs one can sum over a prescribed collection of class pairs of size at most \(n^{107/510}\) by (5), and larger collections when the kernel gcd caps make \(L_1L_2\) smaller than the edge bound.

## A second-path count without prescribing its internal product

**Corollary (Prescribed endpoint gcd on the second path).** Suppose the first internal product \(J_1\) belongs to one prescribed class \([\Theta(J_1)]\). For each first path let a finite list of at most \(L\) possible second endpoint gcds \(G_2\) be determined by that path. If every \(G_2\) on the list is at least \(T\), the count in either determinant branch is
\[
           O_\varepsilon\left(n^\varepsilon e(G)L
                      \left(1+\frac{(R_1R_2)^3}{nT^2}\right)\right). \tag{6}
\]
For singular cycles the bound is \(O_\varepsilon(n^\varepsilon e(G)L)\) without the lower cutoff.

The outer-edge recovery counts the first paths. Each such path fixes the opposite endpoints and hence, by (2), the upper bound at each prescribed \(G_2\). On a singular cycle the same effective slope must hold on both paths, so fixed-endpoint recovery with both \(\rho\) and \(G_2\) gives the sharper bound. For example, in the critical rectangle a subpower list with \(T\ge n^{43/1020}\) suffices for (6) to meet the budget. On the singular branch one can take the entire list from (1). This unconditionally covers that branch within a prescribed first internal-product class when the common effective slope has height at least \(8n^{223/510}\), since there are then at most \(n^{107/510}\) possible values of \(G_2\). There is no second internal-product class condition here.

## Examples and a residual criterion

For the [three-edge path at \(n=157\)](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#exact-examples-and-remaining-count), the four primitive labels are \((22,77,42,12)\), \(A=7\), and \(G=2\). The effective slope is \(11/6\) and (1) reads \(7^2+27\cdot2^2=157\). For the actual nonsingular pair at \(n=26497\) in the same source, the two primitive vectors are \((1148,42,1260,3285)\) and \((3285,1260,42,1148)\). Both have \(G=1\), and \(C=457\) in (1). Their effective slopes \(1148/3285,3285/1148\) are distinct and reciprocal, so their effective norm squareclasses agree. At another opposite pair, left kernel \(13\) and right kernel \(13\), the two directed paths also give those same reciprocal effective slopes. Thus two distinct opposite pairs can carry the same prescribed effective class on a genuine nonsingular cycle.

**Corollary (Residual after effective slope classes).** In the [paired norm squareclass residual](ERDOS_943_PAIRED_PATH_NORM_CLASSES.md#exact-examples-and-the-residual-criterion), one can delete cycles admitting a path at opposite endpoints for which the positive squarefree representative \(z_\rho\) of \([\Theta(\rho)]\) satisfies
\[
              z_\rho\le n^{1/1530}
             \qquad\text{or}\qquad
                  z_\rho\gcd(b,d)\le n^{14/153}.
\]
The pointwise critical-budget count for the smaller residual is equivalent to (H) in the same coefficient family.

The first test sums the one-class count over at most \(n^{1/1530}\) representatives. For the second, partition the representative and the kernel gcd dyadically. A box satisfying the test has a count at most \(O_\varepsilon(n^\varepsilon K^2 ZL)\) by (4), with \(ZL\ll n^{14/153}\). The logarithmic number of boxes costs a subpower factor. Both deleted sectors meet the critical budget. Subset inclusion gives the converse. Bounding the frequency of arbitrary effective slope classes on common-endpoint paths remains open.
