# Cube classes of products in coherent prime paths

The product of the six determinant coefficients of a coherent four-point configuration has a cube class determined by its cross-ratio. Each such class contains only a subpower number of cross-ratios of bounded height. A related product on a cofactor plane gives pointwise power savings for some occupied rank-one records at the quadratic third- and fifth-torsion crossovers.

## The class of a product across a sum

Write \([w]_3\) for the class of \(w\ne0\) in \(\mathbb Q^\times/\mathbb Q^{\times3}\). Its representative is the unique positive cube-free integer in that class. Write \(H(a/b)=\max(|a|,b)\) in lowest terms with \(b>0\).

**Lemma (Subpower multiplicity of the class of \(J(1+J)\)).** Uniformly in a rational cube class, there are \(O_\varepsilon(M^\varepsilon)\) rationals \(J\ne0,-1\) of height at most \(M\ge2\) such that \([J(1+J)]_3\) is that class. The trivial class never occurs.

Write \(J=a/b\) in lowest terms. The class is \([ab(a+b)]_3\), and \(a,b,a+b\) are pairwise coprime. Every prime in the class belongs to exactly one of the three integers, and its exponent modulo three is prescribed. Its prime support has product at most \(2M^3\) if there is a solution. Thus by the divisor bound there are \(O_\varepsilon(M^\varepsilon)\) assignments of the prime support, each fixing positive cube-free integers \(A,B,C\) with
\[
            a=A u^3,\qquad b=B v^3,\qquad a+b=C w^3 .
\]
The roots are primitive as a triple, nonzero, signed, and bounded in absolute value by \((2M)^{1/3}\). Apply [Heath-Brown, Theorem 1.3](sources/heath-brown-squarefull-2601.07817v1.pdf) to \(Au^3+Bv^3-Cw^3=0\), with bound uniform in these coefficients. Since the assignment and roots recover \(a/b\), the lemma follows. In the trivial class \(A=B=C=1\), and the signed Fermat equation for cubes rules out every candidate.

If the positive representative is only required to have radical at most \(R\ge1\), the class count costs \(O_\varepsilon(R^{1+\varepsilon})\). For each possible squarefree radical \(m\le R\) there are \(2^{\omega(m)}\) choices of the exponents one and two.

## Four and five coherent vertices

Use the [coherent norm pencils](ABC_COHERENT_NORM_PENCILS.md#a-pencil-through-two-split-primes), with pairwise distinct oriented prime ideals over distinct rational primes. All coefficients \(v_{ij}\) are nonzero, and within a selected coset they are represented by \(2\)-by-\(2\) determinants of rational vectors up to a common scalar. For a quadruple define
\[
 A=v_{12}v_{34},\quad B=v_{14}v_{23},\quad C=v_{13}v_{24},\quad
 J=A/B,\qquad \chi=[ABC]_3 .
\]
By Plücker, \(A+B=C\), so \(\chi=[J(1+J)]_3\). Changing the sign of any edge, hence any ordering of the four vertices, preserves the cube class of the product of the six edges. The product class also depends only on the projective points: rescaling one of the four determinant vectors rescales the product by a cube, and a common change of basis rescales it by a sixth power.

**Theorem (Four-point cube products and diversity).** Suppose \(|v_{ij}|\le V\) with \(V\ge2\). The cube class of the product of the six edges is nontrivial. For any fixed three vertices and prescribed class, there are \(O_\varepsilon(V^\varepsilon)\) possible fourth vertices from a coherent set. A coherent set of \(s\ge4\) vertices contains \(\gg_\varepsilon s^{24/11}V^{-\varepsilon}\) distinct such product classes on its four-element subsets.

The height of \(J\) is at most \(V^2\). The class gives only \(O_\varepsilon(V^\varepsilon)\) possible \(J\) by the lemma. Fixing the three anchors, \(J\) is a fractional linear coordinate of the fourth projective direction, so each value gives at most one vertex. In fact the determinants between the anchors are nonzero and the other projective directions are distinct, as witnessed by the nonzero \(v_{ij}\). This gives the completion count even when allowed primes or the field window impose further restrictions on the vertices of the coherent set.

Over all ordered quadruples, \(A/C=J/(1+J)\) gives the cross-ratios of these projective directions. Move all the finitely many directions into one affine real chart. [Rudnev, Theorem 1](sources/rudnev-distinct-cross-ratios-1705.01830v1.pdf) gives \(\gg s^{24/11}\) distinct real cross-ratios, and each product class accounts for only subpower many of them. This proves the lower bound on classes. In particular, in a coset of size at least \(V^\delta\) for fixed \(\delta>0\), each class has power-saving density among its quadruples, since it occurs on only \(O_\varepsilon(s^3 V^\varepsilon)\) of them.

## Application to occupied rank-one records

Take an occupied full-modulus record for the third or fifth power, with primes in \([Z,2Z)\), \(Z\ge2\), and exact positive normalized moduli \(q_i\ll V\), \(V\ge2\). Use the [plane through the first two lattice directions](ABC_COHERENT_NORM_COFACTOR_INDEX.md#counting-with-the-index-and-a-pair-of-coordinates) and the [same construction for fifth powers](ABC_FIFTH_COFACTOR_CUBIC_COUNTS.md#small-cofactor-quotients-in-rank-one). Reorder the cofactor coordinates so that the \(i\)-th coordinate is opposite modulus \(q_i\). Write them \(n_i^*=q_i a_i\), put \(P=\|n\|\), and write \(Q=q_1q_2q_3\). The exact lattice refines the original lattice by a bounded index at two. The actual prime vector is primitive in its first direction and is the unique such direction in the rank-one ball. Therefore \(\lambda_1\asymp Z\), \(P\asymp\lambda_1\lambda_2\), and \(P^2\ll ZQ\).

For \(\ell=3,5\), put
\[
                 A_i=a_i^\ell q_i^{\ell-2},\qquad
                 \chi_{ij}=[A_iA_j(A_i+A_j)]_3 .
\]
Use \(\chi_{ij}\) only when the three displayed factors are nonzero. The cofactor congruence gives \(\sum_i A_i\equiv0\bmod Q\), so \(q_k\mid A_i+A_j\) when \(\{i,j,k\}=\{1,2,3\}\). The index \(g\) of the plane in the integer plane lattice has a factorization \(g=t_1t_2t_3\) with \(t_i\mid q_i\), \(g/t_i\mid a_i\). Also \(\gcd(a_i,a_j)\) equals \(t_k\) or \(t_k r_k\), where \(r_k\) is the opposite prime. These statements follow from the [plane index](ABC_COHERENT_NORM_COFACTOR_INDEX.md#exact-index-and-divisibility) and the [pair gcd formula and its odd-power transfer](ABC_COHERENT_NORM_DUAL_PRODUCT_LATTICE.md#the-extra-prime-in-a-pair-gcd). The index formula uses only the congruence ratios of the actual prime direction and so applies unchanged at the fifth power.

**Proposition (A small radical on a cofactor pair).** Fix a pair \(i,j\) with \(A_iA_j(A_i+A_j)\ne0\). For any prescribed cube class \(\chi_{ij}\), the rank-one record count, up to a factor \(O_\varepsilon((ZV)^\varepsilon)\), is at most \(V^{5/3}\) at \(\ell=3\) and \(Z^{1/3}V^{4/3}\) at \(\ell=5\). If \(g\le G\), \(G\ge1\), and \(\gcd(a_i,a_j)=t_k\), one can replace the factor \(V\) in these bounds by \(\min(V,G)\). In particular, the sector where the radical of the positive representative of \(\chi_{ij}\) is at most \(V^{3/2}\) at \(V=Z^{5/4}\) contains \(O_\varepsilon(V^{179/60+\varepsilon})\) records. At the fifth-power crossover \(V=Z^{7/4}\), the threshold \(V^{9/5}\) gives \(O_\varepsilon(V^{619/210+\varepsilon})\) records. These estimates permit arbitrary further field or path restrictions.

To prove the count, fix dyadic scales \(Q_i\asymp q_i\) and \(P\), and choose the ordered pair \(i,j\). Since \(|a_i|\ll P/Q_i\), both \(A_i,A_j\) are bounded by fixed powers of \(ZV\). The lemma gives at most a subpower number of possible *exact* ratios \(J=A_i/A_j\) within the prescribed class.

At the third power, \([q_i/q_j]_3=[J]_3\). As \(q_i,q_j\) are coprime, allocating the support of the class to these integers fixes their cube-free parts at subpower cost. The remaining cube roots give \(O_\varepsilon((ZV)^\varepsilon(Q_iQ_j)^{1/3})\) pairs. Each pair fixes the primitive signed ratio \(a_i/a_j\), since its cube is \(Jq_j/q_i\). At the fifth power, \([a_i/a_j]_3=[J]_3^2\). Count primitive ratios \(a_i/a_j=u/v\), \(v>0\), with \(|u|\ll P/Q_i\), \(v\ll P/Q_j\) and this class, by the same allocation. There are \(O_\varepsilon((ZV)^\varepsilon(P^2/(Q_iQ_j))^{1/3})\) possibilities. The resulting \(J/(u/v)^5\) determines \((q_i/q_j)^3\), so the coprime pair \(q_i,q_j\) is unique if it exists. The factor \((P^2/(Q_iQ_j))^{1/3}\) is \(O((ZV)^{1/3})\).

In either case write \(a_i/a_j=u/v\) in lowest terms, \(v>0\), and \(h=\gcd(a_i,a_j)\). Then \(q_k\) divides
\[
                   h^\ell (u^\ell q_i^{\ell-2}+v^\ell q_j^{\ell-2}),
\]
a nonzero integer. The common sign of \(a_i,a_j\) has no effect on divisibility. We may charge \(O(V)\) possible \(q_k\) at once. If \(h=t_k\) and \(g\le G\), instead choose \(h\le G\), then choose \(q_k\) by the divisor bound. For every resulting modulus triple the coherent cube or fifth-power congruences allow only a subpower number of root-class choices, and each occupied rank-one record has at most one primitive positive prime direction. Summing the dyadic scales proves the per-class bounds.

For the aggregate sector, the complementary large-index count, restricted to the fixed pair with nonzero sum, is
\[
                  O_\varepsilon((ZV)^\varepsilon ZV^3/G)
                  \qquad (g>G).
\]
Indeed, fix \(g,t_1,t_2,t_3\) in each dyadic box \(g\asymp G'\), and choose \((q_i,a_i)\), \((q_j,a_j)\). The coordinate divisibilities give \(O((Q_i/t_i)(Pt_i/(Q_i g)))=O(P/g)\) choices for the first pair, and the same for the second pair. Choose \(q_k\) among divisors of the nonzero \(A_i+A_j\). The choices of \(g\) and the split cost \(O_\varepsilon(G'^{1+\varepsilon})\), so the box costs \(O_\varepsilon((ZV)^\varepsilon P^2/G')\), and \(P^2\ll ZQ\) gives the bound. This proof works at both powers because \(\ell-2\ge1\). The exceptional \(h=t_k r_k\) sector costs \(O_\varepsilon(V^{14/5+\varepsilon})\) at the third power or \(O_\varepsilon(V^{18/7+\varepsilon})\) at the fifth by the [extra-prime count](ABC_COHERENT_NORM_DUAL_PRODUCT_LATTICE.md#the-extra-prime-in-a-pair-gcd).

For \(g\le G\) and \(h=t_k\), sum over the \(O_\varepsilon(R^{1+\varepsilon})\) classes of radical at most \(R\). Here \(R\) can be cut off at a fixed power of \(ZV\). With \(Z=V^{4/5}\), use \(G=V^{49/60}\) and \(R=V^{3/2}\). The small-index contribution is \(O_\varepsilon(V^{2/3+49/60+3/2+\varepsilon})\), and the large-index one is \(O_\varepsilon(V^{4/5+3-49/60+\varepsilon})\). Both exponents are \(179/60\). At \(Z=V^{4/7}\), use \(G=V^{131/210}\) and \(R=V^{9/5}\). The exponents \(11/21+131/210+9/5\) and \(4/7+3-131/210\) both equal \(619/210\). More generally, balancing \(R V^{2/3}G\) against \(ZV^3/G\) at the third power gives \(O_\varepsilon(V^{67/30+\varepsilon}R^{1/2}+V^{14/5+\varepsilon})\) at its crossover. At the fifth power balance \(R(ZV)^{1/3}G\) against \(ZV^3/G\) to obtain \(O_\varepsilon(V^{43/21+\varepsilon}R^{1/2}+V^{18/7+\varepsilon})\). The balancing choice has \(G\ge1\) whenever these estimates improve on \(V^3\). Thus the small-radical bounds save a fixed power up to \(R\le V^{23/15-\eta}\) and \(R\le V^{40/21-\eta}\), respectively, for any fixed \(\eta>0\). A union over the three pairs costs only a constant. The frequency of the classes of large radical remains open.

## Exact examples

In \(\mathbb Q(\sqrt{-23})\) take the oriented prime-cube generators with doubled coordinates
\[
\begin{array}{c|rrrrr}
p & 3 & 13 & 29 & 31 & 41\ \hline
u & 4 & 74 & 282 & 344 & 426\
v & 2 & 12 & 28 & 6 & 64
\end{array}
\]
They satisfy \(u^2+23v^2=4p^3\), with \(p\nmid v\). The \((3,13,29,31)\) and \((3,13,29,41)\) quadruples share three anchors. Their triples \((A,B,C)\), divided by the common positive gcd, are \((49625,54448,104073)\) and \((-19125,24436,5311)\), respectively. The product classes have positive representatives \(281203372686\) and \(19856278188\).

For the [fifth-power example](ABC_FIFTH_COFACTOR_CUBIC_COUNTS.md#counting-a-box-of-cofactor-normals) on primes \(19,23,29\) with moduli \(2,3,5\), the exact cofactor plane has \(a=(1,1,-1)\). The cubic cofactor terms for that example are \((8,27,-125)\); \(J=8/27\) has trivial cube class, but the cube class of \(J(1+J)\) has representative \(35\). For the third-power example with the same primes, moduli and plane, the cofactor terms are \((2,3,-5)\) and the corresponding representative is \(30\). These plane examples illustrate the invariant independently of the rank-one restriction.
