# Coefficient counts and recovery for the cofactor forms

The scaled cofactor plane gives binary forms with integral divided coefficients at odd prime powers. For cubics, two bounded evaluations give a small dictionary of occupied ideal classes. From degree five onward, the catalecticant recovers the three marked linear forms, and their weights and minors recover the normalized moduli and primes.

## Integral coefficients in a cofactor basis

Use the [dual cofactor plane](ABC_COHERENT_NORM_DUAL_PRODUCT_LATTICE.md#products-and-the-three-larger-moduli) with weights \(r_iq_i\), where the \(r_i\) are the distinct opposite primes and the \(q_i\) the absolute values of the normalized signed coefficients. Write \(Q=\prod_i q_i\), \(N=\prod_i r_i\), and choose an integral basis \(u,v\) of the plane lattice. Put \(L_i(S,T)=u_iS+v_iT\). For every odd prime \(\ell\), the coherent \(\ell\)-power congruences give the form
\[
        F_\ell(S,T)=\frac1Q\sum_i q_i^{\ell-2}L_i(S,T)^\ell
                =\sum_{h=0}^\ell \binom{\ell}{h} f_h S^{\ell-h}T^h
        \quad\text{with all } f_h\in\mathbb Z.                      \tag{1}
\]
Indeed modulo \(q_i\), the plane makes \(L_j\equiv-(r_kq_k/(r_jq_j))L_k\) coefficientwise. The coherent norm relation \(r_j^\ell q_j^2\equiv r_k^\ell q_k^2\bmod q_i\) gives cancellation between the two other terms, also before multiplication by each binomial coefficient. CRT gives (1). At \(\ell=3\) this is the [cofactor cubic](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#exact-cubic-and-quadratic-forms) up to an overall sign from the signed coefficient product.

Each \(L_i\) is primitive, since the coordinate projection of the plane lattice onto the \(i\)-th coordinate is surjective by \(\gcd(r_jq_j,r_kq_k)=1\). If \(\Delta_{jk}=u_jv_k-u_kv_j\), the primitive plane normal gives \(|\Delta_{jk}|=r_iq_i\).

## A two-value count for quadratic three-torsion

Assume \(V=Z^{5/4}\) and \(q_i\asymp V\) at the original crossover, with the original discriminant window and allowances. Use the [irreducible path selection](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#irreducible-path-selection), which keeps the original lower-bound scale. Choose a basis of the projected lattice \(\Lambda/\mathbb Zr\) in \(r^\perp\) with product of Euclidean lengths \(O(Q/\|r\|)\). Such a basis exists by two-dimensional lattice reduction. Taking \(r\times w\) and then dividing the coordinates by \(q_i\) gives a basis \(u,v\) of the scaled cofactor plane. In particular
\[
       U=\|u\|_\infty,\qquad W=\|v\|_\infty,\qquad UW\ll ZV.          \tag{2}
\]
Indeed the projected lattice has covolume \(Q/\|r\|\). The cross product with \(r\) scales lengths by \(\|r\|\), and each \(q_i\asymp V\).

**Proposition (Small products of two cubic evaluations).** For any two fixed independent integer vectors \(e,e'\in\mathbb Z^2\), the occupied rank-one records in this balanced irreducible range for which
\[
                          |F_3(e)F_3(e')|\le M
\]
in a cofactor basis satisfying (2) number \(O_\varepsilon(V^{12/5+\varepsilon}M)\) for \(M\ge1\). The implied constant may depend on the two fixed vectors and the constants in (2). In particular the product bound \(M\le V^{1/2}\) gives a saving of \(V^{1/10}\).

The mixed coefficients satisfy
\[
             |f_1|\ll U^2W/V^2,\qquad |f_2|\ll UW^2/V^2.
\]
Irreducibility implies \(F_3(1,0)\) and \(F_3(0,1)\) are nonzero integers, so \(U,W\gg V^{2/3}\). Both mixed-coefficient scales are therefore at least a constant, and the number of \((f_1,f_2)\) choices is \(O((UW)^3/V^4)=O(Z^3/V)=O(V^{7/5})\). Dyadically fixing \(U,W\) only costs logarithms.

The pair \((F_3(e),F_3(e'))\) consists of nonzero integers, so there are \(O(M(1+\log M))\) choices. For fixed \(f_1,f_2\), the two evaluations determine \(f_0,f_3\) by a nonsingular linear system. Its determinant is \((e_1e'_2)^3-(e_2e'_1)^3\), nonzero by independence and injectivity of the real cube map. Thus there are \(O_\varepsilon(V^{7/5+\varepsilon}M)\) possible cubic coefficient tuples in the range where the bound is nontrivial.

A cubic tuple determines the cofactor order and the ideal class up to orientation through its primitive Hessian. Equivalently apply [Bhargava–Varma, Theorem 9 and Corollary 11](https://arxiv.org/pdf/1401.5875) to (1). The [fixed-class mark count](ABC_COHERENT_NORM_MARKED_REPRESENTATIONS.md#a-count-with-one-or-two-prescribed-primes) then gives \(O(V)\) paths per tuple. Multiplying proves the result. The argument counts the condition on any basis satisfying (2) and does not require a short second lattice direction. With the same proof one can instead require \((1+|f_1|)(1+|f_2|)\le M\), since the product of the pure-coefficient height bounds is also \(O((UW)^3/V^4)\).

## Exact fifth-power recovery from the catalecticant

**Proposition (The quintic determines the primes and normalized moduli).** For a coherent fifth-power path of either sign, form \(F_5\) in (1) in any integral scaled-cofactor basis. Its coefficient tuple determines the three pairs \((r_i,q_i)\) up to permutation. Put \(C=(f_{i+j})_{0\le i\le2,\ 0\le j\le3}\). The kernel polynomial of \(C\) is, up to scalar, the split cubic
\[
                           B(X,Y)=\prod_i(v_i X-u_iY).
\]
It is primitive, and the gcd of the absolute maximal minors of \(C\) is exactly \((QN)^2=\operatorname{disc}(B)\). In fact the three nonzero Smith invariants of \(C\) are
\[
                  1,\qquad \gcd(Q,N)^2,\qquad
                            \operatorname{lcm}(Q,N)^2.             \tag{3}
\] Consequently a dictionary of integral quintic tuples in arbitrary cofactor bases gives a bound for the number of records of the same order of magnitude as the dictionary size.

To prove this, put \(\rho_i=q_i^3/Q>0\). The matrix \(C\) factors as the product of the \(3\)-by-\(3\) evaluation matrix with columns \((u_i^2,u_iv_i,v_i^2)\), the diagonal matrix of \(\rho_i\), and the \(3\)-by-\(4\) evaluation matrix with rows \((u_i^3,u_i^2v_i,u_iv_i^2,v_i^3)\). The first determinant is, up to sign, \(\prod_{j<k}\Delta_{jk}\), hence has magnitude \(QN\). The second has determinant one, since \(\prod_i q_i^3=Q^3\). The last matrix has rank three because its three primitive projective points are distinct. This gives the stated kernel. It is the catalecticant calculation behind Sylvester's decomposition of a binary quintic into three fifth powers. See [du Plessis–Wall, Lemma 1.1](https://d-nb.info/1149632828/34) for that classical decomposition.

The signed maximal minors of the last matrix are the coefficients of \(B\) times \(\prod_{j<k}\Delta_{jk}\), up to a common sign. For instance the minor in the first three columns is the quadratic evaluation determinant times \(\prod_i u_i\), which also gives the coefficient of \(Y^3\) in \(B\) up to sign. The coefficient-vector identity follows first when every \(u_i\ne0\), by the one-dimensional kernel, and then throughout by the polynomial identity. Gauss's lemma gives primitivity of the product of the three primitive factors. This proves the asserted gcd of minors of \(C\).

For (3), work at a prime \(p\) of \(QN\), belonging to a unique \(r_iq_i\). Put \(e=v_p(q_i)\) and \(\eta=v_p(r_i)\in\{0,1\}\). Use an integral \(p\)-adic change of variables so that \(L_k=S\), \(L_j=cS+hT\), with \(c\) a unit and \(v_p(h)=e+\eta\). The \(T\)-coefficient of \(L_i\) is a unit because its minor with \(L_k\) is \(r_jq_j\) up to sign. These basis changes act through invertible matrices on the symmetric powers in the catalecticant, so they preserve the Smith invariants. In the divided coefficients \(f_j\), the weights of \(L_i\) and the two incident forms have valuations \(2e,-e,-e\), respectively.

If \(\eta=0\), then \(e\ge1\), \(f_1\) is a unit, and all \(f_j\) with \(j\ge2\) are divisible by \(p\). The catalecticant therefore has rank two modulo \(p\). If \(\eta=1\), the combined pure \(S^5\) coefficient of the incident pair is a unit: \(c\equiv-r_kq_k/(r_jq_j)\bmod p^{e+1}\) and the coherent defect \(r_j^5q_j^2-r_k^5q_k^2\) has valuation exactly \(e\), by the [opposite-prime lift](ABC_COHERENT_NORM_DUAL_PRODUCT_LATTICE.md#transfer-to-higher-odd-powers). For \(e=0\), \(F_5\bmod p\) is now a sum of two nonzero fifth powers of independent directions, again giving rank two.

For \(\eta=1,e\ge1\), the term from \(L_i\) vanishes modulo \(p\), and \(f_0\) is a unit, \(v_p(f_1)=1\), and \(v_p(f_j)\ge2\) for \(j\ge2\). Every \(2\)-by-\(2\) minor is divisible by \(p^2\), since a diagonal product involving \(f_0\) must use a coefficient with index at least two. When \(e\ge2\), \(f_2\) is divisible by \(p^4\), so \(f_0f_2-f_1^2\) has valuation two. When \(e=1\), \(v_p(f_5)=2\) by the unit \(T\)-coefficient in \(L_i\), and \(f_0f_5-f_2f_3\) has valuation two. Thus the gcd of the entries is one, and the gcd of the \(2\)-by-\(2\) minors has valuation two exactly at the primes of \(\gcd(Q,N)\). Together with the maximal minors this proves (3).

Factor \(B\) to recover the three primitive \(L_i\) up to signs and order. Since the fifth powers of three distinct linear directions are linearly independent, matching coefficients then gives \(\rho_i\) up to the corresponding signs. The absolute coefficient ratios \(|\rho_i/\rho_j|=(q_i/q_j)^3\) recover the positive primitive triple \(q\). The three \(r_i\) follow from \(r_i=|\Delta_{jk}|/q_i\). The signs of the original normalized coefficients have only boundedly many choices, and the prime tuple and signed coefficients determine the field by the [exact trace recovery](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#exact-full-moduli-and-the-lattice-count). This proves the dictionary bound, including records whose quotient on a shortest cofactor is large.

The integrality in (1) also gives a consistency check on an inverse candidate. Recover the signs of \(L_i\) by making all \(\rho_i\) positive. If the signed cross product lies in the positive prime direction up to an overall sign, and the recovered \(r_iq_i\) are pairwise coprime, then modulo \(q_i\) one has \(L_j\equiv-r_kq_k L_k/(r_jq_j)\). For each prime power in \(q_i\), one coefficient of \(L_k\) is a unit by primitivity. Using the corresponding pure coefficient of the integral \(F_5\) gives \(r_j^5q_j^2\equiv r_k^5q_k^2\bmod q_i\). Thus the original congruences are forced; the exact norm lift and the field and rank conditions can then be tested on the recovered tuple.

For the [fifth-power example](ABC_FIFTH_COFACTOR_CUBIC_COUNTS.md#counting-a-box-of-cofactor-normals) with \(r=(29,19,23)\) and \(q=(2,3,5)\), take \(u=(1,1,-1)\), \(v=(57,-58,0)\). The quintic is
\[
             \frac{8(S+57T)^5+27(S-58T)^5-125S^5}{30}.
\]
Its kernel is generated by \(Y(57X-Y)(-58X-Y)\). Its maximal-minor gcd is \(380190^2\), recovering \(Q N=30\cdot29\cdot19\cdot23\), and the weight ratios give \(q=(2,3,5)\).

The overlapping-prime case of (3) occurs over \(\mathbb Q(\sqrt{-2})\) on primes \((p_1,p_2,p_3)=(3,11,17)\). Use prime generators \(1+\sqrt{-2},3+\sqrt{-2},3+2\sqrt{-2}\) and take the fifth powers of their coherent products. The traces are \((-10550,-6214,20174)\), the normalized coefficients are \((-58,8247,-269)\), and \(\delta=82\). Here \(r_2=3\) divides \(q_2=8247\), and the gcd of the \(2\)-by-\(2\) catalecticant minors is \(9\).

## Higher odd prime powers

**Corollary (Determinantal recovery from degree seven onward).** For every odd prime \(\ell\ge7\), the tuple of \(F_\ell\) determines the three pairs \((r_i,q_i)\) up to permutation for either sign. For \(C_\ell=(f_{i+j})_{0\le i\le2,\,0\le j\le\ell-2}\), the three nonzero Smith invariants are
\[
                  1,\quad \gcd(Q,N)^2,\quad
                            Q^{\ell-5}\operatorname{lcm}(Q,N)^2.
\]
The maximal-minor gcd alone determines \(Q\) and \(N\).

Put \(m=\ell-2\). The factorization for \(C\) now uses the same quadratic evaluation matrix, the weights \(q_i^{\ell-2}/Q\), and the degree-\(m\) evaluation matrix \(M_m\) with rows \((u_i^m,u_i^{m-1}v_i,\ldots,v_i^m)\). Its kernel consists of the homogeneous forms of degree \(m\) that vanish at the three points \([u_i:v_i]\), so it is the space of multiples of \(B\) in that degree. The gcd of those forms recovers \(B\) up to scalar. Matching the \(\ell\)-th powers again recovers the weights, since the three degree-\(\ell\) evaluation rows are independent. The absolute weight ratios are \((q_i/q_j)^{\ell-2}\), giving \(q\) and then \(r\) by the same primitive-coordinate and minor argument. The trace recovery gives bounded record multiplicity for each fixed \(\ell\).

The maximal-minor gcd of \(M_m\) is \(QN\). For \(p\nmid QN\), the three points remain distinct modulo \(p\), and homogeneous evaluation at them has rank three in every degree at least two, by interpolation using the products of linear factors vanishing at the other two points. For \(p\mid r_iq_i\), use the same \(p\)-adic coordinates as in the quintic proof and write \(L_i=A S+B_0 T\), with \(B_0\) a unit. The rows for \(L_j\) and \(L_k\) coincide modulo \(h\) after scaling the latter by \(c^m\). Thus all maximal minors are divisible by \(h\) over \(\mathbb Z_p\). In the columns of indices \(0,1,m\), the minor is, up to sign,
\[
                        c^{m-1}h B_0^m-h^m A^{m-1}B_0.
\]
Since \(m\ge2\), its valuation is \(v_p(h)\), exactly the valuation of \(QN\). The determinant of the diagonal weight matrix is \(Q^{\ell-5}\), so the maximal-minor gcd of \(C_\ell\) is \(Q^{\ell-5}(QN)^2\).

For the first two determinantal divisors use the same local coordinates as before. The weight valuations are \((\ell-3)e,-e,-e\). If \(\eta=0\), the same \(f_1\) unit gives rank two modulo \(p\). If \(e=0,\eta=1\), the exact opposite-prime defect again gives two independent powers with unit weights modulo \(p\) in divided coefficients. For \(e,\eta\ge1\) the exact defect gives \(f_0\) a unit, \(v_p(f_1)=1\), and \(v_p(f_j)\ge2\) for all \(j\ge2\), just as at degree five. Here \(f_2\) has valuation at least \(\min((\ell-3)e,e+2)\ge3\), so \(f_0f_2-f_1^2\) has valuation two. This gives the first two divisors as at degree five. Finally the valuation of the maximal-minor gcd at each prime is \((\ell-3)e+2\eta\). The residue modulo \(\ell-3\) determines \(\eta\), and subtraction and division give \(e\). This separates the two products, including their shared primes.

For the same prime generators over \(\mathbb Q(\sqrt{-2})\) at \(\ell=7\), the signed edge coefficients are \((130148,101629662,-1319162)\). Thus \(\delta=2\) and \(q=(65074,50814831,659581)\). The \(2\)-by-\(2\) and \(3\)-by-\(3\) determinantal divisors are \(9\) and \(Q^4(3\cdot11\cdot17)^2\), respectively.
