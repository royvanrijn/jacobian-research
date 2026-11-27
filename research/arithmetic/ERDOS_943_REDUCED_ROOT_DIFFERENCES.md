# Six-cycles with small reduced root differences

Fix \(N\) and the [balanced coprime kernel graph](ERDOS_943_SIX_CYCLES.md) on left kernels \(b\) and right kernels \(d\). Select one actual representation \(a^2b^3+c^2d^3=N\) per edge, with all four coordinates in \([T,2T)\), \(N\asymp T^5\). The root at a left vertex on an incident edge is \(a\), and the root at a right vertex is \(c\). For consecutive edges \(e,e'\) on a simple cycle, let \(R,R'\) be their roots at the shared vertex, and set
\[
                  D(e,e')=\frac{|R-R'|}{\gcd(R,R')}.
\]
This is a positive integer: equal roots at a shared vertex would give equal complementary summands, and the uniqueness of the squarefree kernel would give the same edge. For a simple six-cycle \(C\), let \(D_*(C)\) be the smallest product of the three \(D(e,e')\) across any run of four consecutive edges.

**Proposition.** For every target \(N\) in the balanced range and every \(Y\ge1\),
\[
 \#\{C:D_*(C)\le Y\}\ll_\varepsilon N^\varepsilon e(G)Y,\qquad
                   \sum_C D_*(C)^{-1}\ll_\varepsilon N^\varepsilon e(G).
\]
In particular, \(D_*(C)\le T^{1/2}\) is a pointwise \(O_\varepsilon(T^{2+\varepsilon})\) sector, including singular and nonsingular cycles and composite targets.

Fix a selected edge and which endpoint is shared with the next edge. Given the root \(R\) at that vertex and an exact positive difference \(D=D(e,e')\), write \(R'/R=r/s\) in lowest terms. The denominator \(s\) divides \(R\), and \(r=s\pm D\). Thus there are at most \(2\tau(R)\) possibilities for \(R'\). For each one, subtract its actual summand from \(N\). If the complementary summand is positive and has the required square-times-squarefree-cube decomposition, it determines the other endpoint and its root uniquely. Check coprimality, the box constraints and the selected-edge convention before continuing. This recovers only actual edges.

Orient a cycle so that its first four edges realize \(D_*(C)\). Starting at \(e_1=(b_1,d_1)\), the reconstructed edges can be written \((b_2,d_1),(b_2,d_2),(b_3,d_2)\), or with left and right reversed. The last two edges are determined by a common neighbor of \(b_3,b_1\). The coefficient-uniform actual codegree bound in the [kernel codegree theorem](ERDOS_943_KERNEL_GRAPH.md#common-neighbors-and-reconstruction) gives \(O_\varepsilon(N^\varepsilon)\) choices for this neighbor. Its hypothesis \(16T^4<N\) holds for all sufficiently large \(T\) at the fixed balance constants. Bounded \(T\) is absorbed in the implied constant. Reject repeated vertices. For each fixed triple of positive differences the three reconstruction steps and this closing step cost subpower many possibilities per first edge. The number of ordered triples with product at most \(Y\) is \(O(Y(1+\log Y)^2)\) by summing the two reciprocal variables. Every difference is \(O(T)\), so the logarithms may be absorbed into \(N^\varepsilon\). The weighted assertion also follows by summing the reciprocal of the ordered product over three differences up to \(O(T)\). Finally the codegree bound gives \(e(G)\ll_\varepsilon T^{3/2+\varepsilon}\) by counting pairs of neighbors.

Reduced differences can be small even when the numerator and denominator of the ratio are both large. The reconstruction here uses the reduced denominator as a divisor of an actual root and then enforces the exact addition equation, in place of a congruence-only count. This criterion is local along four edges and makes no assumption on the global determinant \(U^2-V^2\). By partitioning the original residual cycles at the threshold, a bound of order \(T^{2+\varepsilon}\) for the original residual is equivalent to the same bound on its subset with \(D_*>T^{1/2}\) in the [remaining high-lcm domain](ERDOS_943_SIX_CYCLES.md#exact-remaining-domain-and-both-implications). The subset direction is immediate, and this proposition bounds the complementary sector.
