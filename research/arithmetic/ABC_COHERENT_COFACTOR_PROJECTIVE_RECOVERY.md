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

## Counting a relative cube class at the third power

For \(\ell=3\), write \(Q=q_1q_2q_3\), \(A_i=a_i^3q_i\) and \(m=(A_1+A_2+A_3)/Q\in\mathbb Z\). Define the normalized pair sums
\[
             L_i=\frac{A_j+A_k}{q_i},\qquad
             L_i+a_i^3=m q_jq_k\qquad(\{i,j,k\}=\{1,2,3\}).        \tag{5}
\]
These are integers by the [cofactor congruence](ABC_COHERENT_NORM_COFACTOR_INDEX.md#exact-index-and-divisibility). Since \(\chi_{jk}=[Q L_i]_3\), the relative class of two cofactor pair products is exactly the cube class of \(L_i/L_j\) for the two opposite indices. The shared \(m\) in (5) gives an aggregate count even when the two absolute cube classes have large prime support.

**Theorem (Relative cube-class count for occupied cofactor profiles).** Consider exact-modulus occupied rank-one records for coherent cubic prime triangles, with the primes in \([Z,2Z)\), normalized moduli \(q_i\ll V\), and \(Z,V\ge2\). Use the cofactor of the plane of the first two successive-minimum directions. Fix distinct indices \(i,j\) with \(a_i a_j L_i L_j m\ne0\). If \([L_i/L_j]_3\) is prescribed and the radical of its positive cube-free representative is \(S\), the number of records is
\[
   O_\varepsilon\left((ZV)^\varepsilon
         \left(Z^2 V^{4/3}S^{-1/3}+ZV^2S^{-1/6}\right)\right).      \tag{6}
\]
In particular, the bound with this radical at most \(R\ge1\) is
\[
   O_\varepsilon\left((ZV)^\varepsilon
         \left(Z^2 V^{4/3}R^{2/3}+ZV^2R^{5/6}\right)\right).        \tag{7}
\]
Both bounds are uniform under further field-window or witness conditions. At \(V=Z^{5/4}\), two equal pair-product cube classes cost \(O_\varepsilon(V^{44/15+\varepsilon})\). More generally, \(R\le V^{1/10-\eta}\) for fixed \(\eta>0\) gives a fixed power saving over \(V^3\) at this crossover.

Write \(L_i=c d_i u_i^3,\ L_j=c d_j u_j^3\) with \(u_i,u_j\) nonzero signed integers and \(c,d_i,d_j\) positive. Here \(c\) is the gcd of the two positive cube-free parts, and \(d_i,d_j\) are their coprime quotients. The prescribed class \([d_i/d_j]_3\) fixes \(d_i,d_j\) up to subpower many choices. Each prime of its support can lie on at most one side, with its exponent one or two then specified. In particular \(\operatorname{rad}(d_id_j)=S\), so \(d_id_j\ge S\). No coprimality between \(c\) and \(d_i d_j\) is assumed.

Take dyadic sizes \(Q_i\asymp q_i\) and \(P\asymp\|n\|\) for the moduli and the cofactor before dividing its coordinates by the moduli. The [successive-minimum bound](ABC_COHERENT_NORM_COFACTOR_INDEX.md#counting-with-the-index-and-a-pair-of-coordinates) gives \(P^2\ll ZQ_1Q_2Q_3\). Use the caps
\[
         T_i\asymp P/Q_i,\qquad
         M_i\asymp P^3/(Q_i\min(Q_j,Q_k)^2)
\]
with fixed constants so that \(|a_i|\le T_i,\ |L_i|\le M_i\). All integers subject to a divisor bound below are polynomially bounded in \(ZV\). For any recovered triple of moduli, cube-root congruences and CRT leave subpower many full records, allowing signs and the bounded refinements at two. Occupation in rank one leaves at most one positive primitive prime direction on each such record. Each record supports only \(O(1)\) successive-minimum planes, so the same estimates bound the number of occupied projective profiles of these planes. In fact, project the exact lattice perpendicular to its shortest vector \(r\). Any nonzero projection lifts after subtracting a multiple of \(r\) to a vector with longitudinal part at most \(\lambda_1/2\). Since its length is at least \(\lambda_2\), the projected minimum is at least \((\sqrt{3}/2)\lambda_2\). Only \(O(1)\) projected vectors can then have length at most \(\lambda_2\), by lattice packing.

### Large common cube-free part

Take \(c>Q_k\). Reduce \(u_i,u_j\) by their positive common gcd to \(v_i,v_j\). The same \(m q_k\) in the two instances of (5) gives
\[
 m q_k(d_j v_j^3q_j-d_i v_i^3q_i)
                       =d_j(v_j a_i)^3-d_i(v_i a_j)^3.             \tag{8}
\]
Choose \(a_i,a_j,v_i,v_j\). If the right-hand side is zero then \(d_i=d_j=1\) by coprimality and cube-freeness. Because \(m\ne0\), this forces \(q_i/q_j=(v_j/v_i)^3\), so both moduli are perfect cubes. This gives \(O_\varepsilon((ZV)^\varepsilon V^{5/3})\) records with \(S=1\), absorbed by (6).

Otherwise factor the right-hand side of (8) to recover \(m,q_k\) and \(d_jv_j^3 q_j-d_iv_i^3q_i\) at subpower cost. Put \(D=\gcd(d_i v_i^3,d_j v_j^3)\). The last linear equation in \(q_i,q_j\) has at most
\[
                O\left(1+\frac{V D}{\max(d_i |v_i|^3,d_j|v_j|^3)}\right)
\]
solutions in the box. The unit contribution, using \(|v_l|\le (M_l/(Q_k d_l))^{1/3}\) for \(l=i,j\), costs
\[
                   \frac{T_iT_j(M_iM_j)^{1/3}}
                        {Q_k^{2/3}(d_id_j)^{1/3}}                  \tag{9}
\]
up to the subpower factor.

For the progression contribution, \(D=D_iD_j\) with \(D_i=\gcd(d_i,v_j^3)\) and \(D_j=\gcd(d_j,v_i^3)\), since \(d_i,d_j\) and \(v_i,v_j\) are each coprime pairs. Fix \(D_i,D_j\) at divisor cost. Then \(v_i\) is a multiple of \(\operatorname{rad}(D_j)\) and \(v_j\) of \(\operatorname{rad}(D_i)\). Bound the maximum below by the geometric mean of the two terms. Summing the resulting reciprocal powers \(3/2\) over nonzero multiples gives
\[
 \sum_{v_i,v_j}\frac{D}{\max(d_i |v_i|^3,d_j|v_j|^3)}
    \ll_\varepsilon (ZV)^\varepsilon (d_i d_j)^{-1/4}.
\]
Indeed \(D_l/\operatorname{rad}(D_l)^{3/2}\le D_l^{1/4}\le d_l^{1/4}\) since \(D_l\) is cube-free. Thus this contribution is at most \(O_\varepsilon((ZV)^\varepsilon T_iT_j V/(d_i d_j)^{1/4})\).

### Small common cube-free part

Fix \(c\le Q_k\) and write \(U_i=(M_i/(cd_i))^{1/3}\), \(U_j=(M_j/(cd_j))^{1/3}\). The common index factor at the third modulus is \(t=\gcd(q_k,a_i,a_j)\). For \(p\mid t\), put \(e=v_p(q_k)\) and \(s=v_p(t)\). At \(e>3s\), the congruence \(q_k\mid a_i^3q_i+a_j^3q_j\) and the definition of \(t\) force \(v_p(a_i)=v_p(a_j)=s\). Since \(A_k\) is divisible by \(p^e\) and \(q_i,q_j\) are units, \(v_p(L_i)=v_p(L_j)=3s\). The same conclusion with \(s=0\) holds at a prime dividing \(q_k\) outside \(t\). At \(e\le3s\) both \(L_i,L_j\) are divisible by \(p^e\). If \(e<3s\), their valuations are exactly \(e\), because primitivity of the scaled cofactor implies that \(a_k\) is a unit at \(p\mid t\).

Choose \(t\) and the depths \(\min(e,3s)\) for the primes of \(t\). For a fixed \(t\), the depth choices cost at most \(\prod_{p^s\Vert t}(2s+1)=\tau(t^2)\). These depths fix \(b=\prod_{p^s\Vert t} p^{\lfloor\min(e,3s)/3\rfloor}\). We have \(t\mid a_i,a_j\) and \(b\mid u_i,u_j\), since the positive coefficients \(cd_i,cd_j\) are cube-free. Choose \((a_i,u_i)\) at cost \(O(T_i U_i/(tb))\). By (5), factoring their nonzero cubic sum recovers the ordered moduli \(q_j,q_k\) and \(m\) at subpower cost.

For the other pair use \(q'_k=q_k/\gcd(q_k,t^3)\) and
\[
            (t/b)^3(a_j/t)^3+cd_j(u_j/b)^3\equiv0\bmod q'_k .       \tag{10}
\]
Indeed \(b^3\mid q_k\), and division by \(b^3\) in the second identity of (5) gives this congruence. All bases and coefficients in (10) are units at \(q'_k\) by the case \(e>3s\). Hence there are subpower many unit cube-root ratio classes for \(a_j/t,u_j/b\), each giving a lattice of index \(q'_k\). If the actual candidate pairs in the scaled rectangle span rank two, the planar successive-minimum bound gives \(O(T_j U_j/(tbq'_k))\) pairs per lattice. If they span one line, the unscaled pairs \((a_j,u_j)\) are integer multiples of a primitive integer vector. Each scale \(f>0\) then satisfies \(f^3\mid m q_i q_k\). With \(N_0=|m|q_k\) already fixed, \(f/\gcd(f,N_0)\) has cube dividing \(q_i\) by prime valuations. Hence \(f\) has \(O_\varepsilon((ZV)^\varepsilon V^{1/3})\) choices, by choosing a divisor of \(N_0\) and an integer of size \(O(V^{1/3})\). In both cases the pair recovers \(q_i\) by (5).

The rank-two term for fixed depths contains the multiplier \(\gcd(q_k,t^3)/(t^2 b^2)\) after extracting \(T_iU_iT_jU_j/q_k\). At \(p^s\Vert t\) this is \(p^{e'-2s-2\lfloor e'/3\rfloor}\), with \(e'=\min(e,3s)\). The exponent is at most \(-s\), except at \(e'=3s-1\) where it is \(1-s\). In the exceptional case \(v_p(L_i)=e\), so \(p^2\mid cd_i\). Thus the multiplier is bounded by \(\gcd(t,\operatorname{rad}(cd_i))/t\). Its sum over \(t\ll V\) costs at most \(\tau(\operatorname{rad}(cd_i))(1+\log O(V))\), by first choosing a divisor of \(\operatorname{rad}(cd_i)\). The line term sums \(1/(tb)\) over \(t\), which also costs only a logarithm. This argument retains the entire index even when \(t\) has large prime factors.

For fixed \(c,d_i,d_j\), the resulting count up to subpower factors is
\[
               \frac{T_iU_iT_jU_j}{Q_k}+T_iU_i V^{1/3}.
\]
We may start with either \(i\) or \(j\). Sum over \(c\le Q_k\). The volume term gives (9) again. Taking the better starting index and bounding the minimum by its geometric mean, the line term is bounded by
\[
          \frac{\sqrt{T_iT_j(M_iM_j)^{1/3}}\,V^{1/3}Q_k^{2/3}}
               {(d_i d_j)^{1/6}} .
\]

### Summing the profile bounds

In (9) the height factor satisfies
\[
 \frac{T_i T_j (M_iM_j)^{1/3}}{Q_k^{2/3}}
      \ll Z^2\bigl(\max(Q_i,Q_k)\max(Q_j,Q_k)\bigr)^{2/3}
      \ll Z^2 V^{4/3}.
\]
This follows from \(P^2\ll ZQ_1Q_2Q_3\). The same bound implies that the small-\(c\) line term is at most \(O(ZV^2/(d_id_j)^{1/6})\). For the large-\(c\) progression term, \(T_iT_j\ll ZQ_k\), giving \(O(ZV^2/(d_id_j)^{1/4})\). These estimates prove (6) after the dyadic sums. There are \(2^{\omega(S)}\) possible relative classes with radical \(S\); summing their weights \(S^{-1/3},S^{-1/6}\) gives (7). The cutoff \(R\) may first be capped at a fixed polynomial in \(ZV\), since the \(L_i\) have polynomial height.

The count includes the equal-class case at any combined radical. It also handles relative classes with small support even if their common factor \(Q\) in the two absolute classes has large support. It leaves the frequency of profiles for which all relative radicals are large. The separately counted \(m=0\) and zero-coordinate sectors are covered by the [plane-index bounds](ABC_COHERENT_NORM_COFACTOR_INDEX.md#counting-with-the-index-and-a-pair-of-coordinates).

## Exact norm witness and relaxed cofactor records

In \(\mathbb Q(\sqrt{-23})\), take the primes \((p_1,p_2,p_3)=(41,47,73)\) and the principal prime-cube generators with doubled coordinates \((426,64),(48,134),(1226,48)\). Taking conjugate products gives traces \((108848,103392,296466)\) and imaginary coefficients \((27006,-80990,-29008)\), in edge order \(12,23,13\). The common scale is \(14\), so \(r=(73,41,47)\) and \(q=(1929,5785,2072)\). The primitive vector \(a=(88,128,-439)\) satisfies (2). Here \(c=(1929,5785,259)\), \(u=(1,1,2)\), and \(b=(44,64,-439)\). The third modulus in (4) is \(259\cdot 2^2=1036\), which divides \(D_1+D_2=1680822976\). This test detects the cube factor \(2^3\) of the third \(q_i\) using the same projective direction as the other two tests. The example uses an arbitrary primitive cofactor vector and needs no successive-minimum assumption.

In the same example with \(r=(73,41,47)\), formula (5) gives \(m=-7\) and \(L=(-84587112,-30075368,6489664)\). Here \(q_3=2072\) and \(\gcd(q_3,a_1,a_2)=8\), so both \(L_1,L_2\) have 2-adic valuation three, the case \(e<3s\). An equal-class example uses \(r=(11,13,17)\), \(q=(1,1,1)\), \(a=(17,17,-24)\). The signed reconstruction for primes in order \((13,17,11)\) has positive field defect \(9778603\) coprime to their product. It gives \(m=-3998\) and \(L=(-8911,-8911,9826)\). These two checks illustrate the normalized sums on genuine witnesses with arbitrary cofactor vectors, independently of the successive-minimum requirement in the count.

For a relaxed family at \(\ell=3\), choose any integer \(u>1\) and put \(q=(u^3,1,1)\), \(a=(1,u,-2u)\). The moduli are pairwise coprime, \(a\) is primitive, and \(A=u^3(1,1,-8)\) satisfies \(\sum A_i=-6q_1q_2q_3\). The projective triple is independent of \(u\). Yet (4) would require \(u^2\mid 1-8\), which is impossible. Thus the cofactor quotient equation and the elementary primitivity data permit unbounded families that cannot occur with a prime direction satisfying both (1) and (2).
