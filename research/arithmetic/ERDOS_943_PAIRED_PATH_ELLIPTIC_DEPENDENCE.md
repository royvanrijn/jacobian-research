# Elliptic dependence of paired path products

When two paths share the same class of \(\Theta(J)=J(1+J^2)\), their internal products give points on one elliptic curve \(E_d:y^2=x^3+d^2x\). Apart from equal and reciprocal products, dependence of those points forces one path product into the image of an elliptic multiplication map. Divisibility by the actual outer root of that path gives a pointwise count for this sector.

Work at a fixed target and coefficient pair in the [selected coprime graph](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-endpoint-equation), with root caps \(R_1,R_2\) satisfying \(2R_1R_2<n\), and set \(R_{\max}=\max(R_1,R_2)\). All congruences and reconstructions use selected actual paths. Use \([c]\) for the positive rational squareclass and \(H(a/b)=\max(|a|,b)\) for reduced fractions with \(b>0\).

## Paths in the image of elliptic multiplication

For an integer \(m\ge2\) let \(f_m\) be the rational map defined by \(x([m]P)=f_m(x(P))\) on \(E_1:y^2=x^3+x\). The same map describes multiplication in the scaled coordinate \(J=x/d\) on \(E_d\), by the isomorphism \((x,y)\mapsto(x/d,y/(d\sqrt d))\) over \(\mathbb Q(\sqrt d)\).

**Lemma (Multiplication maps and primitive residue lattices).** For each fixed \(m\ge2\), \(f_m\) has degree \(m^2\). Take either homogeneous coordinate \(F\) of \(f_m\) in a coprime integral homogeneous pair \((F_m,G_m)\). For any modulus \(w\ge1\), the primitive integer vectors \((a,b)\) with \(w\mid F(a,b)\) lie in \(O_{m,\varepsilon}(w^\varepsilon)\) homogeneous congruence lattices of index at least \(c_m\sqrt w\), where \(c_m>0\). Moreover \(H(f_m(z))\asymp_m H(z)^{m^2}\) for positive finite values.

The map descends through the degree-two \(x\)-map, the quotient by \(P\mapsto-P\). Over \(\mathbb C\), multiplication is unramified of degree \(m^2\), since it is induced by \(z\mapsto mz\) on a torus \(\mathbb C/\Lambda\) with kernel \((m^{-1}\Lambda)/\Lambda\). This uses [Milne, Chapter III, Theorem 3.10 and Corollary 3.9](sources/milne-elliptic-curves-2006.pdf). Comparing degrees and local ramification in \(f_m\circ x=x\circ[m]\) gives degree \(m^2\) for \(f_m\) and ramification index at most two at every point. Thus all projective roots of each of \(F_m,G_m\) have multiplicity at most two.

For the lattice assertion, write \(F=cF_1F_2^2\) with coprime squarefree homogeneous \(F_1,F_2\). They can be taken primitive integral by Gauss's lemma. Put \(S=F_1F_2\). At a prime \(\ell^e\Vert w\), at least one coordinate of a primitive pair is a unit. Choose that chart. The divisibility implies \(\nu_\ell(S(a,b))\ge\lceil(e-\nu_\ell(c))/2\rceil\). In the chosen chart this is a congruence for one of \(S(z,1)\) and \(S(1,z)\), which are fixed squarefree polynomials. Outside the primes dividing \(c\) and their derivative resultants, every root modulo \(\ell\) has at most one lift modulo each higher power. Where a chart polynomial is constant, there is no condition unless its constant is divisible by \(\ell\), and the possible valuations are bounded.

At any remaining fixed prime let \(v\) be the valuation of the nonzero resultant of a nonconstant chart polynomial \(s_0\) and its derivative. Bézout gives \(\nu_\ell(s_0'(z))\le v\) whenever \(\nu_\ell(s_0(z))>v\). If two such solutions to \(s_0(z)\equiv0\bmod\ell^j\) agree modulo \(\ell^{v+1}\), polynomial subtraction gives
\[
 \nu_\ell(s_0(z)-s_0(z'))=\nu_\ell(z-z')+\nu_\ell(s_0'(z)).
\]
They therefore agree modulo \(\ell^{j-v}\). This leaves \(O_m(1)\) classes at this prime after a bounded loss in the modulus exponent. Small \(j\) are absorbed by the same bounded loss. Combining charts and residues gives linear congruences \(a\equiv r_\ell b\) or \(b\equiv r_\ell a\) modulo powers of the distinct primes. Each such lattice has index the product of these powers, at least \(c_m\sqrt w\), and the number of lattice choices is subpower.

Finally \(F_m,G_m\) have no common projective zero. For a primitive pair their evaluated gcd divides a fixed integer, by applying Bézout to the two homogeneous coordinates in both affine charts and choosing a unit chart at each prime. On the real compact set \(\max(|a|,|b|)=1\), the maximum of the two absolute values of the forms has a positive minimum and a finite maximum. Homogeneity and the bounded evaluated gcd give the claimed height comparison.

**Proposition (Actual paths with a multiplication image as product).** For any fixed \(m\ge2\), there are
\[
               O_{m,\varepsilon}\bigl(n^\varepsilon e(G)
                                 (1+R_{\max}^{2/m^2})\bigr)       \tag{1}
\]
actual three-edge paths whose internal product \(J\), or its reciprocal, is a finite positive value of \(f_m\) on \(\mathbb P^1(\mathbb Q)\). The same bound counts cycles if the relative norm squareclass \([\Theta(J_2)/\Theta(J_1)]\) is prescribed and one of their two paths has such an internal product.

For the path count, orient the path so that the denominator \(q\) of \(J=p/q\) is the larger of \(p,q\). This interchanges the roles of the two root caps if necessary. Choose its outer edge first, with actual terminal root \(y_0\). In the [internal product factorization](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-endpoint-equation), \(q=v w\), where \(w\mid y_0\) and \(v\) divides the actual middle terminal root. Choose \(w\) among the divisors of \(y_0\). Then \(q\le R_{\max}w\). When \(J=f_m(a/b)\) or \(J=1/f_m(a/b)\), the height comparison gives
\[
              \max(|a|,b)\ll_m (R_{\max}w)^{1/m^2}.
\]
Also \(w\) divides the corresponding evaluated denominator before reduction. Apply the lemma to this denominator form. A lattice of index \(\gg_m\sqrt w\) contains \(O_m(1+B^2/\sqrt w)\) primitive pairs with \(|a|,b\le B\). In one line there are at most two primitive vectors. In rank two this follows by scaling the lattice and using the two-dimensional lattice-point bound, or by angular packing of the directions with nonzero determinants. Take \(B\ll_m(R_{\max}w)^{1/m^2}\). Since \(2/m^2\le1/2\), the count is \(O_m(1+R_{\max}^{2/m^2})\) per lattice. For each resulting \(J\), [outer-edge recovery](ERDOS_943_EFFECTIVE_PATH_RATIOS.md#starting-a-path-at-its-outer-edge) gives a subpower number of paths. Divisor, lattice and orientation choices give (1).

For the cycle count choose the path in the multiplication image first, in either orientation. Reversing a path preserves its combined class, since \(\Theta(1/J)=\Theta(J)/J^4\). The first path fixes both endpoints and the second path's norm squareclass by the relative class condition. The [uniform squareclass count](ERDOS_943_PAIRED_PATH_NORM_CLASSES.md#a-uniform-bound-in-one-squareclass) and fixed-endpoint recovery leave only subpower many second paths.

## Dependent points on the common twist

**Theorem (Dependent points in the equal norm squareclass sector).** In the critical rectangle, consider cycles with a pair of paths between opposite endpoints having equal \([\Theta(J_i)]\). Let \(d\) be the positive squarefree representative, and use the rational points \(P_i\) on \(E_d\) with \(x(P_i)=dJ_i\) and \(y(P_i)>0\). The cycles for which the two points are dependent modulo torsion and \(J_2\notin\{J_1,1/J_1\}\) number
\[
                               O_\varepsilon(n^{104/255+\varepsilon}). \tag{2}
\]
This applies at the full allowed product heights, to both determinant branches and uniformly per target and coefficient pair.

For \(d\le n^{1/10}\), sum the height-free count for two prescribed absolute classes from [outer-edge recovery](ERDOS_943_EFFECTIVE_PATH_RATIOS.md#starting-a-path-at-its-outer-edge). This costs \(O_\varepsilon(n^\varepsilon e(G)n^{1/10})\).

For \(d>n^{1/10}\), take \(n\) sufficiently large, so \(d>2\). The torsion is generated by \(T=(0,0)\), by [Voutier–Yabuta, §8](sources/voutier-yabuta-1104.4645v4.pdf). The two \(P_i\) are nontorsion as their abscissas are positive. If they are dependent in the free group, the subgroup their images generate is cyclic. Lift a generator to a rational point \(P_0\) and write \(P_i=[m_i]P_0+T_i\) with \(m_i\ne0\) and \(T_i\in\{O,T\}\). The canonical heights satisfy
\[
               \hat h(P_i)=m_i^2\hat h(P_0),\qquad
               \hat h(P_i)\ll\log n,\qquad
               \hat h(P_0)\gg\log(2d).
\]
The upper bound uses the polynomial heights of the path products and \(d\). The lower bound is [Voutier–Yabuta, Theorem 1.2](sources/voutier-yabuta-1104.4645v4.pdf), and the comparison for the upper bound is their Theorem 1.5. Thus \(|m_i|\) ranges over a fixed finite set. Negating the point preserves its abscissa. Adding \(T\) sends the scaled abscissa \(J\) to \(1/J\), since the addition formula gives \(x(P+T)=d^2/x(P)\). If both \(|m_i|=1\), the two products are equal or reciprocal. Otherwise one product or its reciprocal is \(f_m(x(P_0)/d)\) for an \(m\ge2\) in that finite set.

Apply (1) with the relative class trivial. The largest exponent in the path count is at \(m=2\), giving \(O_\varepsilon(n^\varepsilon e(G)(1+\sqrt{R_{\max}}))\). Since \(R_{\max}\le2n^{11/51}\) and \(e(G)\ll_\varepsilon n^{3/10+\varepsilon}\), the bound becomes (2), which also dominates the small-\(d\) part.

## Example and residual

At \(n=113569\), take kernels \((2,1,1,2)\) with oriented roots \((84,239),(288,175),(239,84)\). All three edge sums are \(113569\) and their summands are coprime. The path has \(J=288/175=f_2(7)\), since
\[
               f_2(z)=\frac{(z^2-1)^2}{4z(1+z^2)}.
\]
The squarefree representative is \(d=14\). The point with scaled abscissa \(7\) is \((98,980)\) on \(E_{14}\), and its double has scaled abscissa \(288/175\). This example illustrates the actual path condition in an unshaped graph. The [pair at \(n=26497\)](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#exact-examples-and-remaining-count) illustrates the separate reciprocal case.

**Corollary (Critical residual after the effective-slope and dependence tests).** From the paired norm squareclass residual, remove the cycles covered by the [effective slope class tests](ERDOS_943_EFFECTIVE_PATH_RATIOS.md#examples-and-a-residual-criterion), and also cycles having paired internal products in the same combined class with nontrivially dependent points as in (2). The \(O_\varepsilon(n^{26/51+\varepsilon})\) estimate for the remaining actual cycles is equivalent to (H) with the same coefficient family.

The bounds in the two sources give the critical budget for the deleted cycles. The reverse direction is subset inclusion. In the equal-class part of this residual, one must count equal products, reciprocal products or points independent in \(E_d(\mathbb Q)\) modulo torsion. The unrestricted frequency of the distinct relative classes is also open.
