# Joint power classes in coherent norm pencils

Several power-class comparisons at the same new vertex can be counted on one curve. The compared coefficients are linear functions of the same two pencil coordinates. Each additional anchor contributes a disjoint set of branch points, and the equation for the completion prime contributes a further pair of branch sets.

## Uniform completions with several prescribed ratios

Fix an odd prime \(\ell\). Use coherent witnesses \(\alpha_{ij}=(u_{ij}+v_{ij}\sqrt{D_0})/2\) on distinct odd unramified split primes in \([Z,2Z)\), with \(D_0\) a signed squarefree nonsquare integer and \(|v_{ij}|\le V\). Each witness has the ideal \((\overline{\mathfrak p}_i\mathfrak p_j)^\ell\), and \(\alpha_{ij}\alpha_{jk}=p_j^\ell\alpha_{ik}\). For the real sign one can take witnesses from a common unit-logarithm interval as in the [signed prime-path construction](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#coherent-witnesses-with-a-bounded-unit-logarithm). The argument uses the coefficients of the chosen witnesses and works with any additional path restrictions.

Fix \(k\ge2\) anchors \(1,\ldots,k\) and their witnesses. Let \(p\) be a proposed completion prime, and write \(x_i=v_{ip}\). Prescribe classes of \(x_i/x_j\) in \(\mathbb Q^\times/\mathbb Q^{\times\ell}\) for a set of anchor pairs incident to every anchor. These pairs may be disjoint. Put \(B=\max(Z,V^{2/\ell})\).

**Proposition (Joint power-class completion bound).** The number of completion primes for the given data is
\[
                 O_{\ell,k,\varepsilon}\big((ZV)^\varepsilon B^{\,1/\ell^{k-1}}\big).
                 \tag{1}
\]
The constant is uniform in the field, the anchor data and the prescribed ratio classes. In particular, for \(\ell=3\) and \(V^{2/3}\le Z\), two comparisons involving three anchors give \(O_\varepsilon(Z^{1/9+\varepsilon})\). Two comparisons on disjoint pairs of four anchors give \(O_\varepsilon(Z^{1/27+\varepsilon})\).

For a prescribed pair \(i,j\), the [pencil identity](ABC_COHERENT_NORM_PENCILS.md#a-pencil-through-two-split-primes) gives \(\gcd(x_i,x_j)\mid v_{ij}\). This uses primitivity of the anchored generator and holds for either sign and every \(\ell\) here. Divide out that gcd. The reduced numerator and denominator are coprime, so their \(\ell\)-free parts have support contained in the primes occurring in the prescribed ratio class. At each such prime there are at most two choices of which side it divides, and then its exponent modulo \(\ell\) is specified. The product of primes in the support is at most \(V^2\) if the class is realized. The divisor bound gives a subpower number of choices of the gcd and the two \(\ell\)-free parts for each comparison. Since every anchor appears, at cost \(O_{\ell,k,\varepsilon}((ZV)^\varepsilon)\) we can fix all these gcds and write
\[
                         x_i=A_i r_i^\ell,\qquad A_i>0,\quad
                         0<|r_i|\le V^{1/\ell},
\]
with \(A_i\) fixed. Absorb signs in \(r_i\), and take \(1,2\) to be a pair on which a comparison is prescribed.

Let \(v=v_{12}\), \(u=u_{12}\), \(U=r_1/r_2\), \(t_i=r_i/r_2\) for \(i\ge3\), and \(h=p/r_2^2\). Plücker and the norm from the two-anchor pencil give
\[
\begin{aligned}
 t_i^\ell&=F_i(U)=\frac{v_{1i}A_2-v_{2i}A_1U^\ell}{v A_i}\quad (3\le i\le k),\
 h^\ell&=G(U)=
 \frac{p_2^\ell A_1^2U^{2\ell}-u A_1 A_2 U^\ell+p_1^\ell A_2^2}{v^2}.
\end{aligned}                                                        \tag{2}
\]
For the norm equality, expand \(N(x_1\alpha_{12}-x_2p_1^\ell)=(v^2)(p_1p)^\ell\) and divide by \(p_1^\ell r_2^{2\ell}\). Each \(F_i\) has \(\ell\) simple zeros above its nonzero rational root in \(U^\ell\). These roots are pairwise distinct because
\[
             v_{1i}v_{2j}-v_{1j}v_{2i}=v_{12}v_{ij}\ne0.
\]
The two roots of \(G\) in \(U^\ell\) are distinct, nonzero and nonrational. Indeed the discriminant before division by \(v^2\) is \(A_1^2 A_2^2D_0 v^2\). Thus \(G\) has \(2\ell\) simple zeros, disjoint from all the zeros of the \(F_i\). The disjointness holds in real fields as well.

Over \(\mathbb C(U)\), the extension obtained by adjoining the roots in (2) has degree \(d_*=\ell^{k-1}\). In successive adjunctions, the previously adjoined roots are unramified at a simple zero of the new radicand, since their radicands are units there. The new radicand still has valuation one and its \(\ell\)-th root has ramification index \(\ell\). The full extension is Galois and the roots can be independently multiplied by \(\ell\)-th roots of unity. In particular \(Y=h+\sum_{i=3}^k t_i\) generates the extension. If an automorphism fixes \(Y\), the resulting linear relation among the radicals can be subtracted from its image under an automorphism changing only one radical. The coefficient of that radical must be zero. Doing this for each radical proves that the stabilizer is trivial.

On the smooth model there are \(d_*\) points over \(U=\infty\), because each radicand has pole order divisible by \(\ell\), with unit part admitting \(\ell\)-th roots over the completed field with residue field \(\mathbb C\). At each of these points \(h\) has pole order two and all the \(t_i\) have pole order one. Consequently \(Y\) has pole order two there and no poles at finite \(U\). The map \([U:Y:1]\) is birational to its plane image. The pullback of a general line has degree \(2d_*\), by the poles of its affine equation, so the image is a geometrically irreducible plane curve of degree \(2d_*\) defined over \(\mathbb Q\). Its equation can also be formed starting with \(Y^\ell-G(U)\) and successively taking the resultant of \(T^\ell-F_i(U)\) with the previous polynomial at \(Y-T\).

A completion gives the rational point
\[
                  [\,r_1r_2:p+r_2\sum_{i=3}^k r_i:r_2^2\,]
\]
of height \(O_k(B)\). For a fixed choice of parameters, \(U\) fixes \(x_1/x_2\), and the fixed gcd of \(x_1,x_2\) fixes the pair up to common sign. The homogeneous norm equation fixes \(p\). Thus counting the points bounds the number of completion primes without a scaling multiplicity. Apply [Walsh, Theorem 1.2](https://arxiv.org/pdf/1308.0574v2) to a primitive integral equation of this plane curve. Its coefficient factor is bounded in fixed degree, giving \(O_{\ell,k}(B^{2/(2d_*)})\) uniformly. Summing over the fixed parameters proves (1).

An equality of two opposite products in power class on \(\{1,2,3,p\}\) can prescribe \([x_1/x_2]=[v_{13}/v_{23}]\). An equality on \(\{1,3,4,p\}\) can simultaneously prescribe \([x_3/x_4]=[v_{13}/v_{14}]\). These involve disjoint pairs and so use all four fixed anchors in (1). This is a count of common completions to both comparisons on the same witnesses.

## Density of an adjusted class

Fix one reference vertex \(0\), and on \(s\) other vertices use the adjusted class of \(v_{ij}/(v_{0i}v_{0j})\). Each class is unchanged by swapping the two endpoints, since \(\ell\) is odd.

**Corollary (Sparsity of an adjusted class on large cosets).** For \(m\ge2\), the graph of any one adjusted class has at most
\[
                  O_m\big(s+s^{2-1/m}L^{1/m}\big)
\]
edges, with \(L=O_{\ell,m,\varepsilon}((ZV)^\varepsilon B^{\,1/\ell^m})\). In particular if \(V^{2/\ell}\le Z\) and \(s\ge Z^\delta\) for fixed \(\delta>0\), the density of every class is \(O_{\ell,\delta}(Z^{-\sigma})\) for some \(\sigma=\sigma(\ell,\delta)>0\).

For \(m\) fixed vertices \(j\), a neighbor \(p\) in the specified class \(\chi\) satisfies \([v_{jp}/v_{0p}]=\chi[v_{0j}]\) for every \(j\). Apply (1) with the \(m+1\) anchors to bound their common neighborhood by \(L\). Counting the choices of \(m\) neighbors of a vertex gives \(\sum_p\binom{\deg(p)}m\le L\binom{s}m\). If the mean degree is at least \(2m\), convexity bounds its \(m\)-th power by \(O_m(Ls^{m-1})\). Below \(2m\) the edge count is \(O_m(s)\). For the power saving choose \(m\) so that \(\ell^{-m}<\delta\), then choose \(\varepsilon\) sufficiently small in terms of \(\ell,\delta,m\).

The density bound gives a power saving inside cosets of polynomial size. The existing lower bound for the prime-path argument near the second-moment crossover does not ensure cosets of that size. An aggregate application to the rank-one count would need a way to exploit joint comparisons across the moving fields and the smaller cosets.

## Exact configurations

In \(\mathbb Q(\sqrt{-23})\), use primes \((101,127,131,139,151)\) with oriented prime cube generators given by the doubled coefficient pairs
\[
    (1602,260),\quad (2608,246),\quad (2796,226),\quad (412,678),\quad (3616,174).
\]
With the first four as anchors, the last vertex has \((x_i)=(-330706,-217872,-165356,-1189980)\). The positive cube-free parts are \((330706,27234,165356,1189980)\), so \(r=(-1,-2,-1,-1)\). The \(1,2\) anchor edge has \((u,v)=(2824548,-141994)\); the curve (2) maps this completion to \([2:155:4]\) on the degree-54 plane image. Omitting the fourth anchor gives \([2:153:4]\) on the degree-18 image. For these two inputs, elimination by the stated resultants gives polynomials of the indicated total degrees and vanishes at the displayed points.

For a real fifth-power example in \(\mathbb Q(\sqrt2)\), take the fifth powers of \(3+\sqrt2,5+2\sqrt2,5+\sqrt2,7+2\sqrt2\), of base norms \(7,17,23,41\). With the first three as anchors the new coefficients have positive fifth-power-free parts \((846778,23132072,28745454)\) and roots \((-1,-1,1)\). The coprime radicands in (2) give a degree-50 resultant vanishing at \([1:40:1]\).
