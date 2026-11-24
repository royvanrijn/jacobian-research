# Six-cycles in the square-full kernel graph

Fix a coprime representation box \(a,b,c,d\asymp T\) for
\[
                          a^2b^3+c^2d^3=N,\qquad N\asymp T^5 ,
\]
with \(b,d\) squarefree. The simple graph \(G\) has kernel vertices \(b\) on the left and \(d\) on the right. Choose one representation for each edge. The fixed-kernel and common-neighbor bounds for actual representations are \(O_\varepsilon(N^\varepsilon)\) on this scale. [Kernel graph and recovery](ERDOS_943_KERNEL_GRAPH.md)

**Proposition.** A pointwise \(O_\varepsilon(T^{2+\varepsilon})\) bound for the number of simple six-cycles would imply \(O_\varepsilon(N^{4/15+\varepsilon})\) representations in the balanced box. On average in \(N\) on a dyadic interval, with \(T=X^{1/5}\) fixed, the six-cycle count is \(O_\varepsilon(T^{1/2+\varepsilon})\).

For the first assertion, take the graph to have \(O(T)\) vertices. Above a sufficiently large constant times the \(T^{4/3}\) edge threshold, even-cycle supersaturation gives \(\gg (e(G)/T)^6\) simple six-cycles. This is Morris–Saxton's Theorem 1.5(a) with \(\ell=3\), applied with any needed isolated vertices added so that the vertex parameter is comparable to \(T\). Thus a \(T^{2+\varepsilon}\) bound gives \(e(G)\ll_\varepsilon T^{4/3+\varepsilon}\) after changing \(\varepsilon\). Each edge has at most subpower multiplicity and \(T\asymp N^{1/5}\). [Morris–Saxton, Theorem 1.5](https://arxiv.org/pdf/1309.2927v3)

For the mean, every simple six-cycle projects to a triangle on its three left kernels. Each pair of this triangle must be compatible in the sense of \(K_N(T)\), with the appropriate fixed constants in its kernel intervals and fraction height. The three choices of shared neighbors lift a compatibility triangle to at most \(O_\varepsilon(N^\varepsilon)\) cycles by common-neighbor recovery. A graph with \(K\) ordered adjacencies has \(O(K^{3/2})\) triangles, for instance by bounding the trace of the cube of its adjacency matrix by the \(3/2\)-power of the trace of its square. Therefore
\[
                       \#C_6(G)\ll_\varepsilon N^\varepsilon K_N(T)^{3/2}.
\]
Since \(K_N(T)\ll T^2\) pointwise and has mean \(\ll_\varepsilon T^{1+\varepsilon}\), this already gives \(O_\varepsilon(T^{2+\varepsilon})\) on average. A uniform \(K_N(T)\ll_\varepsilon T^{4/3+\varepsilon}\) would suffice for the desired pointwise cycle bound by the same inequality. [Kernel-pair mean](POWERFUL_SUM_GRAPHS.md)

For the stronger mean, let \(L\asymp T\) bound both side sizes, and let \(K\ge1\) bound the two codegrees of \(G\). We show
\[
                           \#C_6(G)\ll K^{5/2} L^{3/2}e(G).
\]
Call a right vertex low if its degree is at most \(D=4\sqrt{K L}\).
To count cycles through a low right vertex \(d\), choose its two left neighbors \(b,b'\). The number of paths of length two from \(b\) to another left vertex is at most \(K L\), by the left codegree bound. Once that path ends at \(b''\), the final right vertex joining \(b''\) and \(b'\) has at most \(K\) choices. The total contribution is at most \(K^2 L\sum_{d\ \mathrm{low}}\deg(d)^2\le K^2 LD\,e(G)\).

If there are \(h\) high right vertices and \(e_h\) edges incident to them, the right codegree bound gives \(e_h^2\le L(e_h+K h^2)\) by Cauchy–Schwarz on the left degrees. Thus \(hD\le e_h\le L+h\sqrt{K L}\), so \(h\ll\sqrt{L/K}\). To count cycles with all three right vertices high, choose one edge and the other two right vertices. The other two left vertices are common neighbors of specified pairs, so there are at most \(O(e(G)h^2 K^2)\) such cycles. These two estimates prove the displayed inequality.

Across all \(N\) in the dyadic interval, each quadruple of labels \(a,b,c,d\) defines one target \(N\). There are \(O(T^4)\) such quadruples, hence \(\sum_N e(G)\ll T^4\). Since \(K\ll_\varepsilon T^\varepsilon\), the graph inequality gives total \(O_\varepsilon(T^{11/2+\varepsilon})\) six-cycles. Dividing by \(X=T^5\) proves the stronger mean.

## Exact equations for cycle labels

**Proposition.** For a cyclically ordered six-cycle let its edges be \((b_i,d_i)\), \((b_{i+1},d_i)\) for \(i\in\mathbb Z/3\mathbb Z\), with labels \((a_i,c_i)\), \((\tilde a_i,\tilde c_i)\), respectively. Set
\[
 A_i=(a_i\tilde c_i)^2,\quad \widetilde A_i=(\tilde a_i c_i)^2,\quad
 E_i=\tilde c_i^2-c_i^2,\qquad
 U=\prod_i a_i\tilde c_i,\quad V=\prod_i\tilde a_i c_i .
\]
Then \(N\mid U^2-V^2\). If \(U\ne V\), the labels determine all six kernel cubes uniquely, and the left kernel cubes satisfy
\[
 b_i^3=
 \frac{E_i A_{i+1}A_{i+2}+\widetilde A_i E_{i+1}A_{i+2}
                      +\widetilde A_i\widetilde A_{i+1}E_{i+2}}
      {(U^2-V^2)/N}.                                                       \tag{1}
\]
If \(N\) is an odd prime, the denominator in (1) is \(m(U+\sigma V)\) for the unique \(\sigma\in\{1,-1\}\) with \(N\mid U-\sigma V\), where \(m=(U-\sigma V)/N\) satisfies \(|m|\ll T\). When \(U=V\), any fixed set of labels has \(O_\varepsilon(N^\varepsilon)\) lifts to kernel cycles.

To prove this, multiply the two equations at \(d_i\) by \(\tilde c_i^2\) and \(c_i^2\) and subtract. With \(B_i=b_i^3\) this yields
\[
                         A_i B_i-\widetilde A_i B_{i+1}=N E_i .
\]
Successive substitution proves (1) before dividing, with \((U^2-V^2)B_i\) equal to \(N\) times its numerator. Coprimality gives \((B_i,N)=1\), so \(N\mid U^2-V^2\). If the determinant is nonzero this yields \(B_i\) uniquely and each original equation then yields \(d_i^3\). The condition that these quotients are cubes of positive squarefree integers is an exact factor-matching constraint on the labels. For prime \(N\), both \(U,V\) are units, so exactly one sign works and factorization gives the asserted denominator. The bound on \(m\) follows from \(U,V\ll T^6\). For arbitrary \(N\), at each odd prime power dividing \(N\) at most one sign can work. Composite moduli can require different signs on different prime-power factors, so a single congruence \(N\mid U-\sigma V\) cannot be assumed.

If \(U=V\), the first reduced equation still imposes the diagonal cubic
\(A_1 b_1^3-\widetilde A_1 b_2^3=N E_1\). Here \(E_1\ne0\): equality of the two \(c\) labels would give \(a_1^2 b_1^3=\tilde a_1^2 b_2^3\) and hence \(b_1=b_2\), contrary to a simple cycle. Heath-Brown's coefficient-uniform diagonal-cubic bound with third coordinate 1 gives \(O_\varepsilon(N^\varepsilon)\) choices for \(b_1,b_2\). The other cubes then follow uniquely from the linear equations. [Heath-Brown, Theorem 1.3](https://arxiv.org/abs/2601.07817v1)

In the singular case there is also an exact factorization of the edge ratios. The rational relation \(U=V\) lets us write them as \(p_i/q_j\) with positive integers attached to left and right vertices after a common clearing of denominators. This follows by assigning a value at one vertex and propagating along the six edges; consistency at the return is exactly \(U=V\). Put \(g_{ij}=(p_i,q_j)\) on each edge. Coprimality of the original roots gives the labels \(p_i/g_{ij}\), \(q_j/g_{ij}\). The edge equations become \(b_i^3p_i^2+d_j^3q_j^2=N g_{ij}^2\), so the alternating sums of the three \(g_{ij}^2\) on either parity agree. This exact condition remains in addition to the rational product relation.

**Open question.** Bound the number of simple six-cycles in each balanced kernel graph by \(O_\varepsilon(T^{2+\varepsilon})\). The uniform bound from the current codegree theorem is \(O_\varepsilon(T^{3+\varepsilon})\), by choosing three left kernels and their shared neighbors. In the nonsingular case the congruence and \(|m|\ll T\) leave the simultaneous perfect-cube tests in (1) and the right-kernel tests. In the singular case the exact ratio and square-sum identities remain. A bound of the requested size still needs to control the number of label tuples passing these tests.
