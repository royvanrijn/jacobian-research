# A product lattice for coherent prime paths

The exact prime-path congruence lattice has a dual lattice whose coordinatewise products recover the original root ratios. On the plane perpendicular to the prime vector after scaling, products satisfy three stronger congruences, one for each opposite prime. The same plane equation also measures the part of a scaled-normal pair gcd that the original moduli cannot account for.

## Products and the three larger moduli

Use the [exact coefficient moduli](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#exact-cubic-and-quadratic-forms). Reorder the prime vector as \(r=(r_1,r_2,r_3)=(p_3,p_1,p_2)\), put \(q=(|x|,|y|,|z|)\), and write \(D_q=\operatorname{diag}(q_i)\), \(Q=q_1q_2q_3\), and \(N=r_1r_2r_3\). The integers \(r_iq_i\) are pairwise coprime by the incident prime units. For \(\{i,j,k\}=\{1,2,3\}\), the full lattice \(\Lambda\) in this ordering has relation
\[
                    w_j\equiv (r_j/r_k)w_k\bmod q_i .
\]
The coherent norm congruence at \(q_i\) is \(r_j^3q_j^2\equiv r_k^3q_k^2\bmod q_i\). Fractions in congruences denote unit quotients.

**Proposition (Exact dual and product lattices).** Let \(\Lambda^*\) denote the Euclidean dual and put \(\Lambda^\sharp=QD_q^{-1}\Lambda^*\). This is an integer lattice of index \(Q\) with relations
\[
                  a_j\equiv-\frac{r_k q_k}{r_j q_j}a_k\bmod q_i.    \tag{1}
\]
In particular the coordinatewise product of any two elements of \(\Lambda^\sharp\) lies in \(\Lambda\). On the plane
\[
                           \sum_i r_iq_i a_i=0                     \tag{2}
\]
the products generate exactly the lattice
\[
        \Lambda^+=
        \{w\in\mathbb Z^3:\ w_j\equiv(r_kq_k/(r_jq_j))^2w_k
                              \bmod r_iq_i\text{ for each }i\}.      \tag{3}
\]
Thus \([\Lambda:\Lambda^+]=N\), and the prime vector \(r\) generates the quotient \(\Lambda/\Lambda^+\).

For the dual, \((Q/q_i)e_i\in\Lambda\) for the \(i\)-th coordinate vector \(e_i\), so membership in \(QD_q^{-1}\Lambda^*\) forces integral coordinates. With these coordinates, membership says that \(Q\mid\sum_i q_i a_i w_i\) for all \(w\in\Lambda\). Modulo \(q_i\) the \(i\)-th term is zero, and the local pairs \((w_j,w_k)\) range over scalar multiples of \((r_j,r_k)\), independently by CRT. This gives (1) and its index. Multiplying two instances of (1) gives the corresponding relation for \(\Lambda\) because \((r_j/r_k)^3=(q_k/q_j)^2\bmod q_i\).

The integer lattice \(L\) of (2) consists of the scaled cofactors \(D_q^{-1}(r\times w)\) with \(w\in\Lambda\). Indeed \(r\times\mathbb Z^3=r^\perp\cap\mathbb Z^3\), by primitivity of \(r\) and the vector triple product identity. Requiring \(q_i\mid(r\times w)_i\) for each \(i\) is precisely the full congruence condition on \(w\), since \(r_j,r_k\) are units at \(q_i\). For an integral basis \(a,b\) of \(L\), the cross product \(a\times b\) is \(\pm D_qr\), the primitive plane normal. Modulo \(r_iq_i\), equation (2) gives \(a_j=-(r_kq_k/(r_jq_j))a_k\), and the same for \(b\). Hence the three coordinatewise products \(a^2,ab,b^2\) lie in \(\Lambda^+\). Their determinant is, up to sign, the product of the three minors \(a_jb_k-a_kb_j\), by the quadratic Vandermonde identity. Its absolute value is \(QN\). Relation (3) has index \(QN\) by CRT, proving equality.

It remains to identify the class of \(r\). Its defect at \(r_iq_i\) in (3), after multiplying by the unit \((r_jq_j)^2\), is \(r_j^3q_j^2-r_k^3q_k^2\). The \(r_i\)-adic valuation of this difference is exactly \(v_{r_i}(q_i)\). For instance with \(i=1\), in the original signed coordinates it is
\[
     p_1^3 y^2-p_2^3 z^2=p_3^3 x^2-xz u_{23}.
\]
Both \(z\) and \(u_{23}\) are \(p_3\)-adic units, the latter by the edge norm and \(p_3\nmid E y\). The other two identities follow by relabeling. Hence \(r\) satisfies the congruence at \(q_i\), and its defect has order \(r_i\) at the one larger modulus \(r_iq_i\). This remains true when \(r_i\mid q_i\), because the valuation is exactly one short of that modulus. The quotient class of \(r\) has order \(N\), as required.

## The extra prime in a pair gcd

For a rational plane through \(r\), use the [exact plane index](ABC_COHERENT_NORM_COFACTOR_INDEX.md#exact-index-and-divisibility). Write its primitive normal as \(\nu\), its cofactor as \(n=g\nu\), and put \(d_i=\gcd(q_i,\nu_i)\), \(t_i=q_i/d_i\), \(a_i=n_i/q_i\). Then \(g=t_1t_2t_3\) and \(a\) is primitive. In particular
\[
                    \frac{\gcd(a_i,a_j)}{t_k}\in\{1,r_k\}
                       \quad(\{i,j,k\}=\{1,2,3\}).                  \tag{4}
\]
Indeed \(\gcd(a_i,a_j)\) divides \(r_k q_k\) by (2) and primitivity of \(a\). Its gcd with \(q_k\) is \(t_k\) by the index formula. In the second case of (4), if \(r_k\mid q_k\), then \(r_k\mid t_k\).

**Proposition (Power saving for extra prime factors).** In an occupied rank-one record take the cofactor plane through the first two successive-minimum directions of the exact lattice. At \(V=Z^{5/4}\), the number of records with a pair satisfying \(\gcd(a_i,a_j)=t_k r_k\) is \(O_\varepsilon(V^{14/5+\varepsilon})\). The saving holds uniformly under additional field-window and witness restrictions. More generally the bound is
\[
                      O_\varepsilon((ZV)^\varepsilon
                               (V^3/Z+ZV^2+Z^2)).                  \tag{5}
\]

The zero-coordinate and pair-cancellation cases are already bounded by \(O_\varepsilon((ZV)^\varepsilon(V^2/Z+V^{5/3}))\) in the plane-index count. Assume all coordinates and the relevant pair sum are nonzero. Fix dyadic sizes \(Q_\ell\asymp q_\ell\), \(Q_*=Q_1Q_2Q_3\), \(g\asymp G\), and \(P\asymp\|n\|\). Rank one and successive minima give \(P^2\ll ZQ_*\). Fix \(g\) and its split \(t_1t_2t_3\). In the exceptional case of (4), write
\[
    a_i=(g/t_i)r_k b_i,\qquad a_j=(g/t_j)r_k b_j,\qquad
    a_k=(g/t_k)b_k.
\]
The additional prime factor can be removed in the first two entries because it occurs beyond the depth contributed by \(t_k\). Also \(\gcd(b_i,b_j)=1\), by the equality for the pair gcd. The choices of \((q_\ell,b_\ell)\) for \(\ell=i,j\) cost \(O(P/(gZ))\) each, since \(t_\ell\mid q_\ell\) and \(|b_\ell|\ll Pt_\ell/(Q_\ell g Z)\). These choices leave the nonzero integer
\[
            S'=(g/t_i)^3 b_i^3 q_i+(g/t_j)^3 b_j^3 q_j .
\]
The cofactor congruence implies \(q_k\mid r_k^3 S'\).

If \(r_k\mid q_k\), (4) gives \(r_k\mid t_k\), so there are only subpower many choices of \(r_k\) after fixing the split of \(g\). Then \(q_k\) has subpower many choices as a divisor of \(r_k^3 S'\). The moduli and the signed coefficients give subpower many root classes and at most one positive primitive prime direction per occupied rank-one record. Summing \(g\asymp G\) and its splits costs at most a subpower factor times \(P^2/(G Z^2)\ll Q_* /(G Z)\).

If \(r_k\nmid q_k\), choose \(q_k\mid S'\) without choosing \(r_k\). Substituting into (2) and cancelling \(g r_k\) now gives
\[
                          d_i b_i r_i+d_j b_j r_j=-d_k b_k .
\]
Thus the pair \((r_i,r_j)\) belongs to a lattice of index \(d_k\). Its defining functional is surjective modulo \(d_k\), because \(d_i,d_j\) are units there and \(b_i,b_j\) are coprime. There are \(O(1+Z^2/d_k)\) primitive prime pairs. Indeed, if the points lie on a line, there are at most two primitive vectors, and if they span rank two the successive-minima count in the box is \(O(Z^2/d_k)\). For each pair, the original unit cube congruences modulo \(q_iq_j\) give a subpower multiple of \(1+Z/(q_iq_j)\) choices for \(r_k\). This uses the same cube-root and lattice bounds as [Heath-Brown–Pierce, Sections 4.2–4.3](https://arxiv.org/pdf/1409.3177).

Summing \(g\asymp G\) and its factorizations therefore gives, up to the subpower factor,
\[
           \frac{P^2}{G Z^2}
                \left(1+\frac{Z^2G}{Q_k}\right)
                \left(1+\frac{Z}{Q_iQ_j}\right),
\]
since \(d_k=q_k/t_k\gg Q_k/G\). Using \(P^2\ll ZQ_*\) and \(Q_i\ll V\) gives (5) after the dyadic sums. The same formula dominates the bound when \(r_k\mid q_k\). At the crossover each term in (5) is at most \(V^{14/5}\).

For an exact example use the [coherent triple](ABC_COHERENT_NORM_COFACTOR_INDEX.md#exact-index-and-divisibility) with \(r=(29,19,23)\), \(q=(2,3,5)\) and \(d=860107\). The dual plane has basis \((1,1,-1),(57,-58,0)\). Its three products have determinant \(380190=30\cdot19\cdot23\cdot29\). The plane normal with scaled coordinates \((23,253,-137)\) has \(g=1\) and pair gcd \(\gcd(23,253)=23=r_3\). This illustrates the extra prime in (4) on a plane without a successive-minimum restriction.

The case where an opposite prime already divides its modulus occurs with \((p_1,p_2,p_3)=(3,5,7)\), \((x,y,z)=(7,17,4)\), and \(E=d=59\). The signed traces are \((-103,393,-190)\). In reordered coordinates the dual plane is \(49a_1+51a_2+20a_3=0\). The difference \(p_1^3y^2-p_2^3z^2=5803\) has \(7\)-adic valuation one, exactly the valuation in \(q_1=7\). The scaled normal \((-71,49,49)\) on this dual plane has \(g=t_1=7\) and \(\gcd(49,49)=49=t_1r_1\). It realizes the other case of the divisor argument.

## Transfer to higher odd powers

**Corollary (Higher odd powers).** For any fixed odd prime \(\ell\) in the [signed prime-path construction](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#exact-full-moduli-and-the-lattice-count), the product lattice (3) is contained in the congruence lattice with incident ratios \((r_j/r_k)^{\ell-2}\) at \(q_i\). Its quotient has order \(N\) and is generated by the vector with coordinates \(r_i^{\ell-2}\). The extra-prime bound (5) holds for either sign and gives \(O_\varepsilon(V^{18/7+\varepsilon})\) at the fifth-power crossover.

The proof of the exact dual lattice only uses the unit ratios of the original congruences. Their \(\ell\)-th powers now satisfy \((r_j/r_k)^\ell=(q_k/q_j)^2\bmod q_i\), so the squared dual ratios reduce to \((r_j/r_k)^{\ell-2}\). Equation (2), the product determinant, and the order of the quotient are independent of \(\ell\). The valuation calculation for the quotient generator uses \(r_j^\ell q_j^2-r_k^\ell q_k^2\) and gives the same exact valuation \(v_{r_i}(q_i)\): the term \(r_i^3x^2\) in the displayed trace identity becomes \(r_i^\ell x^2\), with the other term still a unit times \(x\).

The index formula and (4) also use only the occupied lattice and the prime direction. In the count at \(\ell\), the pair sum is \(a_i^\ell q_i^{\ell-2}+a_j^\ell q_j^{\ell-2}\). Its vanished case forces \(q_i,q_j\) to be perfect \(\ell\)-th powers, because \(\gcd(\ell,\ell-2)=1\) and the moduli are coprime. This contributes at most a subpower factor times \(V^{1+2/\ell}\). In the nonzero case remove \(r_k^\ell\) from the pair sum and use the same divisor argument and linear congruence on \((r_i,r_j)\). The remaining prime now satisfies unit \(\ell\)-th power congruences with the same subpower bound on root classes. This proves (5), and setting \(Z=V^{4/7}\) gives the asserted fifth-power bound. The full rank-one estimate still requires the records with \(\gcd(a_i,a_j)=t_k\) for every pair.
