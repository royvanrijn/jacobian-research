# Path congruences and deep target valuations in original additive graphs

Target divisors constrain endpoints along whole paths in the
original additive graph. At a fixed prime the relevant precision
is the least precision on a path. Different paths can supply
independent primes to the same endpoint comparison. This gives a
count for fixed overlapping divisor records and transfers deep
target valuations to same-side differences.

## Exact projection of overlapping edge records

Let \(G=(V,E)\) be a finite graph with vertex values \(z_v\in\mathbb Z\).
For each edge \(e\) impose \(z_v-z_u\equiv r_e\pmod{M_e}\) in a
prescribed direction \(u\to v\), with \(M_e\ge1\).
For a prime \(p\) and depth \(j\ge1\) let \(G_{p,j}\) keep precisely
the edges with \(v_p(M_e)\ge j\), retaining all vertices.
The system is consistent if and only if for every \(p,j\)
the oriented residue sum around each cycle of \(G_{p,j}\)
is zero modulo \(p^j\).

For consistency at one prime, let \(e_*\) be the largest
depth on an edge. Construct the vertex labels from low
digits to high digits. At step \(j\) the labels chosen
modulo \(p^{j-1}\) already obey all requirements to that
depth. The equations on their extra digits modulo \(p\)
on \(G_{p,j}\) again have zero sum around every cycle,
by the assumed consistency modulo \(p^j\).
Assign one extra digit freely on each component and
propagate along a tree. After step \(e_*\), combine primes
by CRT. Necessity follows by adding the equations around
the cycle.

For two different vertices \(u,v\) define
\[
 Q(u,v)=\prod_p p^{\max\{j:u,v\text{ lie in the same component of }G_{p,j}\}},
\]
with the maximum zero when there is no such depth.
Every solution satisfies a fixed congruence on \(z_v-z_u\)
modulo \(Q(u,v)\), obtained by adding along a qualifying
path at each prime. It is independent of the paths by
cycle consistency. If values are given on a subset
\(A\subseteq V\), they extend to an integer solution if
and only if they satisfy all these path congruences on
pairs in \(A\). To see sufficiency subtract one full
solution, leaving homogeneous edge congruences. At each
digit step the already prescribed vertices in a component
agree modulo that depth, so they are compatible with the
same digit. The unprescribed digits can be filled in
component by component. This also proves that for any
order \(v_1,\ldots,v_\ell\) of the vertices, after the
first \(i-1\) coordinates are fixed compatibly, the
allowable \(i\)-th coordinates form one residue class modulo
\[
                   Q_i=\operatorname{lcm}_{k<i} Q(v_k,v_i).
                                                               \tag{GL1}
\]
The lattice of homogeneous solutions has index
\[
 [\mathbb Z^V:\Lambda]=
       \prod_i Q_i
       =\prod_p p^{\sum_{j\ge1}(|V|-\#\pi_0(G_{p,j}))}.        \tag{GL2}
\]
For example, at each depth one vertex per component
is first and pays no exponent at that depth. Every
other vertex contributes one, giving the second
identity. This exponent is also the maximum forest
weight with edge weights \(v_p(M_e)\), by taking a spanning
forest in descending weight order.

The homogeneous path and projection formulas are those
of generalized integer splines. The affine version above
records the consistency condition and makes the extension
of prescribed coordinates explicit.
[Rose–Suzuki, generalized integer splines on arbitrary graphs](https://arxiv.org/abs/2111.02471)
· [Kruskal, descending-weight spanning forests](https://www.cs.yale.edu/homes/lans/readings/routing/kruskal-1956.pdf)

## Counting capped additive configurations with shared locks

Fix a bipartite graph of summand labels
\(A_i=F x_i,\ B_j=G y_j\), with target
\(C_{ij}=A_i+B_j>0\) on every designated edge.
The positive original letters are primitive along edges.
The signs may instead be fixed after permuting the
original three roles. Suppose a divisor \(M_{ij}\)
of each \(C_{ij}\) is recorded. For positive \(A_i,B_j\)
as written, set \(z_{A_i}=A_i,\ z_{B_j}=-B_j\).
The homogeneous path congruences apply and the moduli
are prime to \(FG\) by primitivity. Thus for each
\(x_i\) or \(y_j\), the next allowable cofactor is in
one class modulo its \(Q_v\) from (GL1). With cofactor
dictionaries of sizes \(N_v\) and cofactor intervals of
width \(B_v\), the number of label configurations
compatible with a fixed record is at most
\[
                        \prod_v\min\{N_v,1+B_v/Q_v\}.
                                                               \tag{GL3}
\]
Each edge target is then fixed by its two original labels.
The bound is uniform in the positive \(F,G\) and can be
minimized over the vertex order, or applied after
specifying some anchor vertices.

In a fixed finite pattern at height \(X\), use
\(n=X^L\) and assume the union radical of the
recorded divisors is at most \(n^\mu\). There are at
most \(O(n^\mu X^\varepsilon)\) records. First
enumerate the union radical, then allocate its primes
to their occupied edges and record their exponents.
At polynomial height and for a fixed number of edges
the latter costs are subpower. In particular, suppose
the dictionaries have sizes \(N_v\le n^{d_v}\), and
each qualifying configuration admits a record for
which an ordering reconstructs vertices of combined
cost at least \(b\) by \(Q_v>B_v\). Then the total
number of configurations is
\[
                    O(n^{\,\sum_v d_v-b+\mu}X^\varepsilon).
                                                               \tag{GL4}
\]
The record may depend on the configuration. Thus
a near-full dictionary pattern with an exponent
deficit less than \(b-\mu\) cannot satisfy that
joint-divisor condition.

Paths can combine divisors even where two endpoints
have no common neighbor. For example, take two
disjoint length-four paths between two distinct
\(A\)-vertices in an eight-cycle. If the full divisor
\(P\) is present on all four target edges along the
first path and \(Q\) on all four along the second,
then the two \(A\)-cofactors agree modulo
\(\operatorname{lcm}(P,Q)\). When this modulus
exceeds their containing interval's width, distinct
endpoints rule out the configuration. Neither path
has to supply such a modulus by itself, and an
individual common-neighbor comparison between the
endpoints is unavailable. This use of (GL1) only
needs exact target divisibility for the involved edges.

## A difference-moment inequality for the actual edge layers

Now take any actual bipartite additive graph with
distinct positive cofactor labels on each side, and
write \(N_A,N_B\) for its part sizes.
Use the edge set of the graph and put
\(M_{ij}=C_{ij}\) at all primes outside \(FG\).
For each \(p,j\) denote the part sizes on a
component \(D\) of \(G_{p,j}\) by \(a_D,b_D\).
Let \(e_D\) be the number of edges at that depth in \(D\).
For any set \(\mathcal J\) of prime-depth pairs
with \(p\nmid FG\), write
\[
 \mathcal W_{\mathcal J}(t)
     =\sum_{(p,j)\in\mathcal J,\ p^j\mid t}\log p
                                      \qquad(t\ne0).
\]
Then
\[
 \begin{split}
 \sum_{(p,j)\in\mathcal J}\log p\!
      \sum_D\!\left\{\binom{a_D}{2}+\binom{b_D}{2}\right\}
 \ \le{}&
 \sum_{i<i'}\mathcal W_{\mathcal J}(x_i-x_{i'})
       +\sum_{j<j'}\mathcal W_{\mathcal J}(y_j-y_{j'}).
 \end{split}                                                    \tag{GL5}
\]
The sets of distinct pairs on a side are disjoint
across components at each prime and depth. Every
pair counted on the left is congruent modulo \(p^j\)
by the path rule and the fact that \(p\nmid FG\).
This proves the inequality. The right side is
at most \(\binom{N_A}{2}\log(2B_A)+
\binom{N_B}{2}\log(2B_B)\) when the two cofactor
intervals have widths at most \(B_A,B_B\ge1\).

For a component of total size \(a_D+b_D\ge h>2\)
one has
\[
 \binom{a_D}{2}+\binom{b_D}{2}
        \ge (1-2/h)a_Db_D \ge (1-2/h)e_D .
\]
Indeed the first inequality follows from
\((a-b)^2(1-1/(a+b))\ge0\). Smaller components
have at most \(h(a_D+b_D)/4\) edges. Assume the
complete targets are at most \(X^B\) for fixed \(B\)
and all tested primes are at most \(n^\beta\).
Across the tested layers, the weighted number of
edges in those smaller components is at most
\(O_B(h(N_A+N_B)n^\beta\log X)\), since the sum of
\(\log p\) over the tested possible depths at one
prime is at most \(B\log X\). Consequently
\[
 \begin{split}
 &(1-2/h) \sum_{\text{edges }(i,j)}
                         \mathcal W_{\mathcal J}(C_{ij})\\
 &\quad\le
 \sum_{i<i'}\mathcal W_{\mathcal J}(x_i-x_{i'})
 +\sum_{j<j'}\mathcal W_{\mathcal J}(y_j-y_{j'})
 +O_B\!\left(h(N_A+N_B)n^\beta\log X\right).
 \end{split}                                                    \tag{GL6}
\]
Select, for example, all depths at least
\(\lceil\kappa/L\rceil\) at primes up to \(n^\beta\).
This gives a direct transfer from deep target divisor
weight to deep same-side difference weight, retaining
the original role of each factor in the equation.

There is also a power bound for imbalanced components.
In a component with edges and
\(\max(a_D/b_D,b_D/a_D)\ge n^\delta\ge2\),
the larger side supplies at least \(n^\delta e_D/4\)
same-side pairs. Thus the tested prime-depth weight
on edges in these components is at most
\[
      4n^{-\delta}\!
       \left(\sum_{i<i'}\mathcal W_{\mathcal J}(x_i-x_{i'})
        +\sum_{j<j'}\mathcal W_{\mathcal J}(y_j-y_{j'})\right).
                                                               \tag{GL7}
\]
For the near-minimum complete dictionaries of sizes
\(n^{1/6+o(1)}\) with \(n^{1/3-o(1)}\) edges, take
\(0<\beta<1/6\) and \(h=n^\theta\) with
\(0<\theta<1/6-\beta\). The error in (GL6) is
power smaller than the edge count times \(\log X\).
The imbalanced part is also power smaller, for
each fixed \(\delta>0\), by (GL7).
It follows that any tested deep target mass of full
order must appear in large components balanced at the exponent scale
and also pass to the summed difference moment. These are
constraints within the original arithmetic graph.
They allow components that are sparse in the
ambient dictionaries, and (GL6) by itself gives
no power bound on the resulting difference moments.
[Common full factors and the complete-label reduction](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass)
· [Direct LCM amplification on aggregate stars](ABC_AGGREGATE_SHIFT_DICTIONARY_STARS.md)
