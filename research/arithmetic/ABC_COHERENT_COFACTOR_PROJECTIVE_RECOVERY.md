# Recovering coherent prime-path moduli from a projective cofactor

The three terms of the cofactor identity determine the coefficient moduli up to divisor choices when their ratios are specified. The recovery uses both the plane containing the prime direction and the coherent norm congruences. The cubic quotient equation alone allows unbounded families with the same ratios.

## Recovery from the projective terms

Use the [exact-modulus cofactor plane](ABC_COHERENT_NORM_DUAL_PRODUCT_LATTICE.md#products-and-the-three-larger-moduli). For a fixed odd prime \(\ell\), let \(r_i\) be the three opposite primes of a coherent prime triangle and \(q_i\) the absolute normalized coefficients. They satisfy
\[
 q_i\mid r_j^\ell q_j^2-r_k^\ell q_k^2,\qquad
 \gcd(q_i,r_jr_kq_jq_k)=1\qquad(\{i,j,k\}=\{1,2,3\}).
 \tag{1}
\]
The cofactor vector after dividing its coordinates by \(q_i\) satisfies
\[
              \sum_i r_iq_i a_i=0,\qquad A_i=a_i^\ell q_i^{\ell-2}.
 \tag{2}
\]
Work with nonzero \(a_i\) and nonzero pair sums \(A_j+A_k\). The signs of the original coefficients are unrestricted.

**Proposition (Projective cofactor recovery).** Suppose \(q_i\ll V\), \(r_i\in[Z,2Z)\), \(Z,V\ge2\), and the chosen plane has \(|a_i|\le (ZV)^C\) for a fixed \(C\). Fix \([A_1:A_2:A_3]\). At most \(O_{\ell,C,\varepsilon}((ZV)^\varepsilon)\) triples of exact moduli can occur. In the full-modulus rank-one count the same bound holds for occupied records. The result applies to either quadratic field sign.

Write \(q_i=c_i u_i^\ell\), with \(c_i>0\) \(\ell\)-power-free and \(u_i>0\). The \(c_i\) are already determined by the projective triple. Indeed at any prime at least two of the valuations of the \(A_i\) are zero modulo \(\ell\), by the coprimality of the \(q_i\). If the representative of the projective triple is multiplied by an unknown common rational scalar, the most frequent residue among the three valuations determines its valuation modulo \(\ell\) by requiring that this residue become zero. The resulting residues are \((\ell-2)v_p(q_i)\) modulo \(\ell\). Invert \(\ell-2\) modulo \(\ell\) and take representatives between \(0\) and \(\ell-1\). This recovers all \(c_i\).

Now \(A_i=c_i^{\ell-2}(a_i u_i^{\ell-2})^\ell\). The ratios of the \(A_i/c_i^{\ell-2}\) determine the primitive signed integer triple \(b=(b_i)\) proportional to \((a_i u_i^{\ell-2})_i\), up to common sign. The real odd roots of the ratios are unique, and the recovery only keeps data for which these roots are rational. Put \(D_i=c_i^{\ell-2}b_i^\ell\), so \([D_1:D_2:D_3]=[A_1:A_2:A_3]\). Define \(B_i=r_i c_i b_i u_i^2\). After cancellation of the common factor in \(a_i u_i^{\ell-2}\), (2) gives
\[
                    \sum_i B_i=0,\qquad B_i^\ell=D_i r_i^\ell q_i^2.
 \tag{3}
\]
Fix \(i\), and reduce modulo \(c_i u_i^2\), which divides \(q_i\). In (3), \(B_i\) is zero and \(B_j\equiv-B_k\). Raising to the odd power and applying (1) gives
\[
 0\equiv B_j^\ell+B_k^\ell
   \equiv (D_j+D_k)\,r_k^\ell q_k^2\pmod {c_i u_i^2}.
\]
The last factor is a unit, so the full divisibility is
\[
                              c_i u_i^2\mid D_j+D_k.               \tag{4}
\]
The three right-hand sides are fixed, nonzero, and of polynomial size in \(ZV\). For the height assertion, use any actual record in the indicated box. Its \(a_i u_i^{\ell-2}\) are polynomially bounded, and dividing out their gcd only decreases the bounds. Thus there are divisor-many \(u_i\), at subpower cost. The moduli \(q_i=c_i u_i^\ell\) follow. For fixed moduli the unit-root congruences in (1) and CRT give subpower many lattice records, allowing the bounded refinements at two in the cubic construction. The rank-one restriction allows at most one primitive positive prime direction on each. This yields the record count, with any further field-window or witness restriction.

## Two cube classes and their combined support

For nonzero factors, put \(\chi_{jk}=[A_j A_k(A_j+A_k)]_3\), where \([\,\cdot\,]_3\) denotes the rational cube class. A positive cube-free integer represents each class, since signs are cubes. The radical of the product of two such representatives measures the union of their prime supports.

**Corollary (Two cofactor product classes).** In the range of projective cofactor recovery, two prescribed classes \(\chi_{ij},\chi_{ik}\) for distinct pairs leave \(O_{\ell,C,\varepsilon}((ZV)^\varepsilon)\) occupied rank-one records. If the radical of the product of their representatives is at most \(R\ge1\), the bound is \(O_{\ell,C,\varepsilon}((ZV)^\varepsilon R^{1+\varepsilon})\). In particular, at both the third- and fifth-torsion crossovers this gives a power saving on the nondegenerate sector where two such classes have combined radical at most \(V^{3-\sigma}\), for any fixed \(0<\sigma<3\).

For one pair, put \(J=A_j/A_k\). Its class \([J(1+J)]_3\) equals \(\chi_{jk}\). The [uniform cube-class slope bound](ABC_COHERENT_CUBE_PRODUCTS_AND_COFACTOR_COUNTS.md#the-class-of-a-product-across-a-sum) gives a subpower number of exact \(J\)'s, because their heights are bounded by fixed powers of \(ZV\). That bound applies [Heath-Brown, Theorem 1.3](sources/heath-brown-squarefull-2601.07817v1.pdf) to the primitive diagonal cubic left after the two fraction terms and their sum have prescribed cube-free parts. Applying it to both pairs fixes the projective triple at subpower cost. For combined radical \(m\), each prime dividing \(m\) has at most eight choices of its pair of exponents, namely the pairs in \(\{0,1,2\}^2\) other than \((0,0)\). Summing \(8^{\omega(m)}\) over \(m\le R\) proves the second bound. In the crossover applications \(Z\) is a fixed power of \(V\); the subpower allowances can be absorbed into \(V^\varepsilon\). The shortest cofactor plane satisfies the requisite polynomial height by \(P^2\ll Z q_1q_2q_3\), as in the [plane-index estimate](ABC_COHERENT_NORM_COFACTOR_INDEX.md#counting-with-the-index-and-a-pair-of-coordinates).

The count uses the support of the two classes together. It covers some ranges where neither individual radical falls under the [one-class bound](ABC_COHERENT_CUBE_PRODUCTS_AND_COFACTOR_COUNTS.md#application-to-occupied-rank-one-records). It leaves the frequency of large combined support open. A further route to the full rank-one saving would be to count the occupied projective triples themselves, using (4) to discard candidate directions before testing the prime box and field window.

## Exact norm witness and relaxed cofactor records

In \(\mathbb Q(\sqrt{-23})\), take the primes \((p_1,p_2,p_3)=(41,47,73)\) and the principal prime-cube generators with doubled coordinates \((426,64),(48,134),(1226,48)\). Taking conjugate products gives traces \((108848,103392,296466)\) and imaginary coefficients \((27006,-80990,-29008)\), in edge order \(12,23,13\). The common scale is \(14\), so \(r=(73,41,47)\) and \(q=(1929,5785,2072)\). The primitive vector \(a=(88,128,-439)\) satisfies (2). Here \(c=(1929,5785,259)\), \(u=(1,1,2)\), and \(b=(44,64,-439)\). The third modulus in (4) is \(259\cdot 2^2=1036\), which divides \(D_1+D_2=1680822976\). This test detects the cube factor \(2^3\) of the third \(q_i\) using the same projective direction as the other two tests. The example uses an arbitrary primitive cofactor vector and needs no successive-minimum assumption.

For a relaxed family at \(\ell=3\), choose any integer \(u>1\) and put \(q=(u^3,1,1)\), \(a=(1,u,-2u)\). The moduli are pairwise coprime, \(a\) is primitive, and \(A=u^3(1,1,-8)\) satisfies \(\sum A_i=-6q_1q_2q_3\). The projective triple is independent of \(u\). Yet (4) would require \(u^2\mid 1-8\), which is impossible. Thus the cofactor quotient equation and the elementary primitivity data permit unbounded families that cannot occur with a prime direction satisfying both (1) and (2).
