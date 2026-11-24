# Target-independent collisions modulo a cube

Let \(E(B)\) be the number of ordered collisions
\[
 b^3\mid c_1^2 d_1^3-c_2^2d_2^3,\qquad
 b,d_1,d_2\in[B,2B),\quad c_1,c_2\in[C,2C),\quad C\asymp B^{3/4},
\]
with \(b,d_1,d_2\) squarefree, \((c_1c_2d_1d_2,b)=1\), and distinct coordinate pairs. The count is independent of the target in Erdős #943. It is useful because for any collection of \(R_b\) distinct inputs in the same target fiber, \(R_b(R_b-1)\le E_b\), where \(E_b\) is the collision count at this modulus.

**Theorem.** The off-diagonal count satisfies
\[
                            E(B)\ll_\varepsilon B^{71/28+\varepsilon}.        \tag{1}
\]
The proof combines reduced rational heights with the bilinear lattice covering of Heath-Brown. It uses his uniform bound for primitive diagonal cubics, Theorem 1.3 of the cited version. [Heath-Brown, *Counting Square-full Solutions to \(x+y=z\)*, Theorem 1.3 and §4](https://arxiv.org/pdf/2601.07817v1)

## Reduced ratios and elementary averages

Write \(c_i=g_c x_i,\ d_i=g_d y_i\), where \((x_1,x_2)=(y_1,y_2)=1\). Let \(U,V\) be dyadic heights for the two pairs, respectively. Positivity and the initial intervals give comparability within each pair. For a given primitive quadruple there are \(O(CB/(UV))\) choices of \(g_c,g_d\). The congruence after cancelling these unit factors is \(b^3\mid F\), where \(F=x_1^2 y_1^3-x_2^2 y_2^3\). When \(F=0\), the squarefree ratio \(y_2/y_1\) is a rational square, so \(y_1=y_2=1\) and \(x_1=x_2=1\). Those inputs give exactly the diagonal. Throughout the off-diagonal count \(F\ne0\).

The exact unit parametrization is
\[
 r=x_1/x_2,\quad s=y_2/y_1,\quad t=r/s,\qquad r=t^3,\quad s=t^2\pmod{b^3}. \tag{2}
\]
Two primitive fractions of height \(O(V)\) that agree modulo \(b^3\) coincide, since their cross difference is \(O(V^2)=o(b^3)\). The same holds at height \(O(U)\). Specifying either coordinate in (2) leaves at most \(O_\varepsilon(B^\varepsilon)\) choices for the other, by the unit square- and cube-root bounds modulo \(b^3\). For squarefree \(b\) the latter bounds follow from CRT and the at most two or three roots over each prime and unique lifting at primes exceeding 3, with an absolute factor at 2 and 3.

Let \(Q(U,V)\) count the primitive quadruples together with \(b\). The preceding argument gives \(Q\ll_\varepsilon B^\varepsilon B\min(U^2,V^2)\). For each primitive quadruple \(b^3\mid F\ne0\) gives at most \(O_\varepsilon(B^\varepsilon)\) choices of \(b\), so \(Q\ll_\varepsilon B^\varepsilon U^2V^2\). Finally fix \(x_1,x_2\) and the nonzero integer \(F/b^3\), whose absolute value is \(O(U^2V^3/B^3)\). The triple \((y_1,y_2,b)\) is primitive and solves a diagonal cubic with all three coefficients nonzero. Theorem 1.3 gives \(O_\varepsilon(B^\varepsilon)\) choices uniformly in the coefficients. Thus
\[
 Q(U,V)\ll_\varepsilon B^\varepsilon
 \min\left(BU^2,\ BV^2,\ U^2V^2,\ \frac{U^4V^3}{B^3}\right).               \tag{3}
\]
In particular bounded \(g_d\) gives \(O_\varepsilon(B^{5/2+\varepsilon})\) after the \(CB/(UV)\) weight. Other ranges already at the desired \(B^{2+\varepsilon}\) include \(U\ll B^{5/12}\) or \(V\ll B^{1/2}\), by the last term. The largest exponent allowed by (3) overall is \(8/3\), at \(U\asymp B^{3/4}\) and \(V\asymp B^{5/6}\).

For a subfamily of squarefree moduli each possessing a divisor in \([H,2H)\), one can reduce modulo the cube of that divisor first. A specified unit fraction residue modulo \(h^3\) has at most \(O(1+V^2/h^3)\) distinct primitive representatives at height \(O(V)\). Indeed their real values, in a bounded interval here, are spaced by \(\gg h^3/V^2\). Fix the \(x\)-pair and use the cube-root count, then recover \(b\) among the divisors of \(F\). This gives the further bound \(Q\ll_\varepsilon B^\varepsilon U^2(H+V^2/H^2)\) in that subfamily. For example \(H\asymp B^{2/3}\), \(U\asymp C\) and \(V\asymp B\) give \(B^{13/6+\varepsilon}\) including scale factors. This argument requires the stated divisor and leaves prime \(b\) untreated.

## Bilinear covering and averaging short generators

Fix \(b\) and the unit \(\kappa=t\bmod b\). We use the lifted bilinear congruences of Heath-Brown's Lemma 4.1, with coefficient signs \(+1,-1\) and his choice of \(\kappa\) taken with the opposite sign. Their inputs here are the unit collision and the choice of \(t\bmod b\), without any requirement on the quotient \(F/b^3\). In the integer substitutions \(X_1=x_1,\ X_2=\kappa^3x_2,\ Y_1=\kappa^2 y_1,\ Y_2=y_2\), these give
\[
 X_1\equiv X_2,\quad Y_1\equiv Y_2\bmod b,
\]
\[
 X_2(Y_1-Y_2)\equiv0\bmod b,\quad
 2X_1Y_1+X_2Y_1-3X_2Y_2\equiv0\bmod b^2,\quad
 5X_1Y_1-X_1Y_2+X_2Y_1-5X_2Y_2\equiv0\bmod b^3.                         \tag{4}
\]
Here \(\kappa\) is lifted to any integer prime to \(b\). The modulus-\(b^3\) form can also be obtained by setting \(X_1/X_2=1+\alpha,\ Y_2/Y_1=1+\beta\). Both increments are multiples of \(b\), and \((1+\alpha)^2=(1+\beta)^3\bmod b^3\) gives \(4\alpha=6\beta+\alpha\beta\bmod b^3\). The converse part of Lemma 4.1 ensures that **every** integer vector satisfying (4), including nonunit vectors, obeys
\[
                        b^3\bigl(\gcd(y_1,y_2,b)\bigr)^2\mid 18 F.            \tag{5}
\]
The factor 18 keeps the statement valid at primes dividing 6.

For clarity, the projective slicing from the proof of Heath-Brown's Lemma 4.2 applies to this box as follows. The four products \(x_i y_j\) are \(O(UV)\). The three bilinear constraints in (4) are inner-product conditions with rational vectors of denominators \(b,b^2,b^3\). In the coordinate order \((x_1 y_1,x_1 y_2,x_2 y_1,x_2 y_2)\), use the rows
\[
 (0,0,1,0),\quad b^{-1}(0,0,\kappa^2,-1),\quad
 b^{-2}(2,0,\kappa^3,-3\kappa),\quad
 b^{-3}(5\kappa^2,-1,\kappa^5,-5\kappa^3).
\]
The first is an integer coordinate and the next two normalize their bilinear forms by units. The fourth is the last form in (4) itself. The row lattice has determinant \(2b^{-6}\) in absolute value by the four triangular pivots. Scaling the four coordinates by \(UV\) gives determinant \(O((UV)^4 b^{-6})\). Applying the projective slicing lemma (Lemma 2.8) to this lattice, with constants absorbing the box endpoints, gives
\[
                         L=1+(UV)^{4/3} B^{-2}
\]
up to a constant number factor of bilinear relations \(x^{\mathsf T}T y=0\) per \((b,\kappa)\). The integer coefficient matrices \(T\ne0\) are obtained by clearing at most \(b^3\), and have entries of size \(B^{O(1)}\).

A rank-one \(T\) fixes one primitive fraction and gives \(O_\varepsilon(B^\varepsilon)\) points by (2). If \(T\) is invertible, orthogonality to \(T^{\mathsf T}x\) and primitivity of \(y\) give \(hy=J T^{\mathsf T}x\), where \(J\) is rotation by ninety degrees and \(h\) is an integer divisor of \(\det T\), with either sign. The divisibility follows by multiplying by the adjugate and using primitivity of \(x\). This gives \(O_\varepsilon(B^\varepsilon)\) rational planes per relation. On the saturated integral lattice of a plane, the two linear congruences of (4) cut out a sublattice. The remaining forms are binary quadratic forms in any basis of that sublattice. Their zero congruences can be imposed by \(O_\varepsilon(B^\varepsilon)\) sublattices on each of which all the congruences hold identically. At each prime digit a binary quadratic either vanishes identically, vanishes only at the origin, or vanishes on at most two lines modulo the prime. One iterates on the corresponding sublattices, dividing the form by the power already imposed. Divisibility of its values implies divisibility of its coefficients by evaluation at the two basis vectors and their sum. At most six rounds per prime are needed for (4). This is the refinement in Lemma 2.6 and the proof of Lemma 4.2. Always refine the same saturated lattice only once for each given plane and \((b,\kappa)\). All counted planes have injective projection onto both coordinate pairs, and planes on which \(F\) vanishes identically can be discarded. [Heath-Brown, Lemmas 2.6, 2.8, 4.1 and the proof of 4.2](https://arxiv.org/pdf/2601.07817v1)

For a rank-two lattice use the norm \(\max(|x_i|/U,|y_i|/V)\) and its successive minima \(R\le S\). If the box contains no two independent points, it contains at most one positive primitive point since its \(x\)-pair is primitive. Otherwise \(S\ll1\) and the lattice-point bound is \(O((RS)^{-1})\). These are the standard two-dimensional successive-minima bounds, applied after scaling the lattice along the coordinates. The estimates below are for \(U\ll B^{3/4},V\ll B\).

For a short vector of norm \(\le R\ll1\) on such a lattice, \(x_1,x_2\) are both nonzero. If one were zero, the congruence modulo \(b\) would force \(b\) to divide the other, whose size is \(o(b)\), contradicting injectivity of the projection. Put \(r_0=(y_1,y_2,b)\) and suppose \(F\ne0\). By (5), \(f=18F/(b^3 r_0^2)\) is a nonzero integer of size \(O(R^5 U^2 V^3/(B^3r_0^2))\). For fixed \(x_1,x_2,r_0,f\) the triple \((y_1/r_0,y_2/r_0,b/r_0)\) is primitive and again satisfies a diagonal cubic with all coefficients nonzero. Since \(b\) is squarefree, given \(y_1,y_2,b\) the second linear constraint of (4) gives \(O_\varepsilon(B^\varepsilon r_0)\) choices of \(\kappa\). Dividing it by \(r_0\), the remaining \(y\)-coordinates are units modulo \(b/r_0\) whenever the constraint has a solution, so this follows from the unit square-root bound. Summing over \(r_0\) proves that the total number of nonzero-form short vectors, counted with their \((b,\kappa)\), is
\[
                   \ll_\varepsilon B^\varepsilon Q_0 R^7,\qquad
                   Q_0=U^4 V^3/B^3.                                        \tag{6}
\]
Indeed there are \(O((RU)^2)\) choices of the nonzero \(x\)'s and the sum over \(r_0^{-1}\) costs a logarithm.

For a zero-form short vector, \(x_1=j a^3,\ x_2=j c^3\) with coprime signed \(a,c\), since the reduced squared \(x\)-ratio is a rational cube. At dyadic height \(H=\max(|a|,|c|)\), let \(l=(b,j)\). Then \(l\ll RU/H^3\), and for each \(a,c,l\) there are \(O(B/l)\) bases and at most \(O_\varepsilon(B^\varepsilon l)\) possible \(\kappa\) for each base, by \(x_1\equiv \kappa^3 x_2\bmod b\). Summing the \(O(H^2)\) pairs and the choices of \(l\), then the dyadic \(H\), shows that at most \(O_\varepsilon(B^{1+\varepsilon} RU)\) pairs \((b,\kappa)\) admit such a vector. Each class has \(O_\varepsilon(B^\varepsilon(1+U^2/B))\) positive primitive collisions, by spacing the primitive \(x\)-ratios in the prescribed residue class modulo \(b\) and recovering \(y\) by (2).

Split the lattices at \(RS=K\), where \(0<K\le1\). Rank-one relations and lattices with at most one point, together with lattices with \(RS\ge K\), give \(O_\varepsilon(B^{2+\varepsilon}L/K)\). The classes with a zero-form vector of norm \(O(\sqrt K)\) contribute within the preceding zero-form estimate. For the remaining lattices with two independent box points, choose a nonzero-form shortest vector and a second independent vector of size \(O(S)\) with nonzero form. If the vector for the second minimum has zero form, add between one and six copies of the first vector. The resulting polynomial in that integer multiplier has degree at most five and cannot vanish identically, as the form on the plane is nonzero. For fixed dyadic \(R,S\), a pair of these vectors determines the plane, while the number of sublattices for that plane and class is \(O_\varepsilon(B^\varepsilon)\). In each class, the number of resulting lattices is bounded by \(O_\varepsilon(B^\varepsilon)\) times the smaller of \(L\) and the product of the two short-vector counts. The inequality \(\min(L,mn)\le\sqrt{Lmn}\), Cauchy–Schwarz over \((b,\kappa)\), and (6) give a lattice count \(O_\varepsilon(B^\varepsilon L^{1/2} Q_0(RS)^{7/2})\). Multiplying by \(O((RS)^{-1})\) and summing over dyadic \(RS\ll K\) gives the third term of
\[
 Q(U,V)\ll_\varepsilon B^\varepsilon
 \left(\frac{B^2 L}{K}+B U(1+U^2/B)\sqrt K+L^{1/2}Q_0 K^{5/2}\right).          \tag{7}
\]
There are \(O((\log B)^2)\) dyadic pairs of minima, since a nonzero integer vector in the original coordinates has scaled norm \(\gg B^{-1}\). This only costs a subpower factor. The difference from the corresponding averaging in Heath-Brown's square-full sum problem is that all points on the lattices must be counted here, with weight \((RS)^{-1}\).

To optimize (3) and (7), write \(U=B^u,\ V=B^v\) up to constants, with \(u\le3/4,\ v\le1\). The multiplier \(CB/(UV)\) has exponent \(7/4-u-v\). If \(u+v\le3/2\), the last term of (3) gives exponent \(3u+2v-5/4\le5/2\). If \(v\ge27/28\), the \(BU^2\) term gives \(11/4+u-v\le71/28\). Otherwise \(u+v>3/2\), \(3/4<v<27/28\), and \(u>1/2\). Put \(K=B^{-k}\) with \(k=2v/3-3/7>0\). The first term of (7) with the scale multiplier has exponent \(7/4+(u+v)/3+k\le11/7+v\). The third has at most that exponent because
\[
 \tfrac72 k\ \ge\ (10u+7v)/3-4.
\]
The second term has exponent \(7/4+2u-v-k/2\le5/2\). This proves (1) after the dyadic height sum.

## Character orthogonality and structured families

**Proposition.** In the balanced box, comparable \(c\)-pairs with a rational-cube ratio contribute \(O_\varepsilon(B^{2+\varepsilon})\). The off-diagonal with equal \(d\)'s is empty for sufficiently large \(B\), and that with equal \(c\)'s contributes \(O_\varepsilon(B^{3/4+\varepsilon})\). The latter is empty for prime \(b>3\) and sufficiently large \(B\).

If \(c_1/c_2=(a/c)^3\) in lowest terms, of height \(H\ll C^{1/3}\), the reduced \(d\)-ratio gives \(b^3\mid(c^2y_2)^3-(a^2y_1)^3\). The integer on the right can vanish only on the diagonal by squarefreeness. Consequently \(VH^2\gg B\). For each base and root ratio, squarefree \(b\) gives only \(O_\varepsilon(B^\varepsilon)\) primitive \(d\)-ratios by cube roots and uniqueness at height \(V\). The scale factors are \(O(C/H^3)\) and \(O(B/V)\ll H^2\). The \(O(H^2)\) root ratios and \(O(B)\) bases therefore give \(O_\varepsilon(B^{1+\varepsilon} C H)\), which is at most \(O_\varepsilon(B^{2+\varepsilon})\). If \(d_1=d_2\), divisibility of \(c_1^2-c_2^2\) by \(b^3\) forces equality. If \(c_1=c_2\), the nonzero quotient \((d_1^3-d_2^3)/b^3\) is bounded, and \((b,d_1,d_2)\) is primitive by the unit condition. Theorem 1.3 bounds the possible triples by \(O_\varepsilon(B^\varepsilon)\). For prime \(b>3\) neither of \(d_1-d_2,\ d_1^2+d_1d_2+d_2^2\) is divisible by \(b^3\) in these ranges, and their gcd is prime to \(b\) under the unit condition whenever \(b\) divides the first. This proves the prime-modulus assertion.

Any lift of \(t\) in (2) to a primitive rational number \(a/c\) of height \(H\) with \(H^2 B\) less than a sufficiently small constant times \(b^3\) must have \(y_2/y_1=(a/c)^2\) exactly. Squarefreeness forces \(y_1=y_2\) and then \(c_1=c_2\) by the equal-\(d\) argument. Thus off-diagonal parameters in this box have rational lift height \(\gg B\). There is a large exact-equality family if the squarefree condition on \(d_i\) is dropped. For example \((c_1,c_2)=(125m,64m)\), \((d_1,d_2)=(16n,25n)\) fits the two dyadic boxes for \(\gg C\) choices of \(m\) and \(\gg B\) choices of \(n\), and gives exact equality. Every prime \(b\in[B,2B)\) for large \(B\) is a unit modulus for them. The resulting lower bound is \(\gg B^{11/4}/\log B\) for that relaxed problem.

For the character formula, write \(q=b^3\), let \(C_b,D_b\) be the cardinalities of the unit sets of \(c,d\) in the boxes with the prescribed squarefreeness, and set
\[
 A_\chi=|\textstyle\sum_c\chi^2(c)|^2-C_b,\qquad
 T_\chi=|\textstyle\sum_d\chi^3(d)|^2-D_b.
\]
Let \(J_c,J_d\) count the off-diagonal pairs giving equality modulo \(q\) of the powers on the respective sides alone. Orthogonality gives the exact off-diagonal identity
\[
 E_b=\frac1{\varphi(q)}\sum_{\chi\bmod q} A_\chi T_\chi+D_bJ_c+C_bJ_d.           \tag{8}
\]
This uses the full characters modulo \(q\), including imprimitive ones. In the present box \(J_c=0\) and the summed \(C_bJ_d\) term has the bound just proved. The summed contribution of the principal character is \(O_\varepsilon(B^{3/2+\varepsilon})\).

The marginal second moments in (8) are bounded by \(O_\varepsilon(B^\varepsilon C^2)\) and \(O_\varepsilon(B^\varepsilon B^2)\). For the first, \(A_\chi\) sums the characters of powered nonidentity ratios, with their scale weights. Orthogonality counts the squared weights aggregated by the powered residue. Each aggregated residue has at most \(O_\varepsilon(B^\varepsilon)\) primitive fractions by the root bound and height uniqueness. At dyadic height \(H\) the sum of the squared weights is \(O(H^2(C/H)^2)\). Summing in \(H\) proves the bound. The argument for \(T_\chi\) is the same with height scale \(B\). The respective moments are at least \(C_b(C_b-1)\) and \(D_b(D_b-1)\), by pairing each ordered off-diagonal ratio term with itself in the squared-character sum. Thus even after the diagonal is removed, separate marginal bounds and Cauchy–Schwarz give only the \(B^{11/4+\varepsilon}\) scale after summing moduli. To reach the target through (8) requires a mixed estimate for the signed products of the two centered quantities, with a power saving over the estimate using separate absolute moments.

The same fixed unit coefficient on both monomials cancels. If the two coefficients can differ, the exact-equality locus must be subtracted separately unless each term is in its canonical form with the new squarefree coefficient supported on primes of exponent exactly one. For example \(5(49m)^2(5n)^3=7(25m)^2(7n)^3\) where \(5n,7n\) can both be squarefree, with both ratios lying inside the dyadic width. Under the canonical condition, an exact equality forces the coefficients to coincide by unique factorization.

A sufficient asymmetric estimate for a full \(1/3\) theorem would be
\[
 E(B;C,D)\ll_\varepsilon M^\varepsilon\bigl(B^2+(CD)^2/B^2\bigr),             \tag{9}
\]
where the input scales for \(c,d\) are \(C,D\), the target is \(M\), and the estimate is only needed in the remaining ranges \(B\le D,\ M^{1/3}<BD<M^{4/9}\) from the fixed-pair bounds. The previous reduction also permits imposing \(B>M^{4/21}\). Indeed \(\sum_b R_b\ll B+\sqrt{B E(B;C,D)}\), and (9) bounds this by \(M^\varepsilon(B^{3/2}+CD/\sqrt B)\). Both terms are \(O(M^{1/3})\) before the subpower factor, using \(B^2\le BD\), \(C^2D^3\le M\) and the two bounds on \(BD\). The fixed coefficient multiplying \(c^2d^3\) is a unit in the coprime reduction and cancels in the collision relation. [Sixth-power carry and the active pointwise ranges](ERDOS_943_SIXTH_POWER_CARRY.md)

**Open question.** Prove \(E(B)\ll_\varepsilon B^{2+\varepsilon}\) in the original unit box. A sufficient average after (3) is the improvement of (7) in the large-height range, including prime bases with \(U\) near \(B^{3/4}\) and \(V\) near \(B^{27/28}\). Equivalently, one can bound the summed nonprincipal mixed term of (8) from above by \(O_\varepsilon(B^{2+\varepsilon})\). The proposed \(E(B)\) estimate would give \(O_\varepsilon(B^{3/2+\varepsilon})\) representations in the balanced box by Cauchy–Schwarz. A full \(1/3\) pointwise theorem also needs corresponding asymmetric box estimates. The \(71/28\) bound alone gives a weaker pointwise bound in the critical box than the sixth-power carry note's cube-projection bound.
