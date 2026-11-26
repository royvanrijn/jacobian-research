# Arithmetic derivative minors on a shared-value box

When the selected primes are private to complete labels of a
box, its arithmetic derivative matrix is a cycle matrix with
weighted columns. Its maximal minors have an exact
spanning-tree formula. Large prime powers on an individual
label do not contribute to their gcd. For a prime to contribute
beyond the gcd of the local valuation exponents, its high-order
occurrences must contain a cycle of the difference graph.

## The cycle matrix and exact minor formulas

Use an actual shared-value box of size \(m\ge2\), with all
\(m^3\) original positive primitive rows in \(X/2<C\le X\).
Its \(h=3m\) potentials have complete tripartite difference
graph \(G\), with \(E=3m^2\) edges and positive labels \(L_e\).
Orient each edge from the smaller to the larger potential.
Every original additive equation is the corresponding
signed triangle circulation. Select a set \(S_e\) of at least
one prime on each edge label, with each selected prime
private to its edge among all labels of the box. Set
\(S=\bigcup_e S_e\). The private-prime bounds give such
selections on generic capped boxes.
[Actual shared-value arrays (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem)
· [Private radical of each label (PV2)](ABC_SHARED_BOX_PRIVATE_PRIME_SUPPORT.md#simultaneous-privacy-and-the-outside-cap)

The triangle circulations span the integral cycle lattice
of \(G\). A simple cycle of length at least five has a
chord: for a chosen vertex, if the vertices in positions
two and three along the cycle both missed it, they would
be in its own axis and could not be adjacent. Split along
a chord and descend. A chordless four-cycle uses two axes,
since its two pairs of nonadjacent vertices belong to the
respective same axes. Its oriented edge vector is the sum
of four triangles obtained by choosing an apex on the
third axis. Finally, every integer circulation is a sum
of oriented simple cycles, by decomposing its balanced
positive directed edge multiset.

Let \(r=E-h+1\) and take a matrix \(B\) of \(r\)
independent triangle circulations. Write \(I_B\) for the
index of their integral span in the cycle lattice.
For \(p\in S_e\), put
\[
 a_{e,p}=(L_e/p)v_p(L_e),\qquad
 g_e=\gcd_{p\in S_e} a_{e,p},\qquad
 w_e=\sum_{p\in S_e}a_{e,p}^2.
\]
These quantities are positive. Use the same triangle rows
to form the arithmetic derivative matrix \(A\) on \(S\)
as in (AD1). The \(p\)-column of \(A\) is \(a_{e,p}\)
times the \(e\)-column of \(B\), because \(p\) divides
exactly that complete label. The selected prime sets on
each edge are nonempty, so \(A\) has rank \(r\).
[Family arithmetic derivative inequality (AD1)](ABC_ARITHMETIC_DERIVATIVE_FAMILY_LATTICES.md#the-family-inequality)

**DM1.** With \(T\) ranging over the spanning trees of \(G\),
\[
 \Delta(A)=I_B\gcd_T \prod_{e\notin T}g_e,\qquad
 \det(AA^{\mathsf T})=I_B^2\sum_T\prod_{e\notin T} w_e .
 \tag{DM1}
\]
In fact an \(r\)-element edge set \(J\) indexes a nonzero
minor of \(B\) exactly when its complement is a spanning
tree. The restriction of the cycle space to coordinates
on \(J\) is injective precisely when its complement has
no cycle. That complement has \(h-1\) edges, so acyclicity
is equivalent to being a spanning tree. For such a
complement, the fundamental cycles of its tree are an
integral basis whose restriction to \(J\) is unimodular.
Thus the corresponding minor of \(B\) has absolute value
\(I_B\). An \(r\)-minor of \(A\) is zero if two selected
prime columns come from the same edge. Otherwise its
absolute value is the corresponding \(B\)-minor times
the product of the selected coefficients \(a_{e,p}\).
Taking the gcd over the prime choices on any edge set
yields the product of the \(g_e\). Summing the squares
instead and applying Cauchy–Binet gives the second formula.

## The prime contribution is a cycle count

Fix a rational prime \(p\), not necessarily in \(S\).
Let \(G_{p,j}\) be the spanning subgraph on the edges with
\(v_p(g_e)\ge j\), and let \(\beta_1\) denote cycle rank,
so \(\beta_1(H)=|E(H)|-|V(H)|+\#\operatorname{comp}(H)\).

**DM2.** The exact valuation is
\[
       v_p(\Delta(A)/I_B)=\sum_{j\ge1}\beta_1(G_{p,j}).
 \tag{DM2}
\]
Indeed (DM1) minimizes the sum of the orders on edges
outside a spanning tree, or equivalently maximizes the
sum on its edges. Greedily adding edges by descending
order, avoiding cycles, uses a maximal forest in each
level subgraph after processing that level. At each \(j\)
the number of remaining edges in the level subgraph is
then its cycle rank. Summing over the levels gives the
minimum.

There is a direct bound in terms of the original label
valuations. Put \(u_e=\gcd_{q\in S_e}v_q(L_e)\), so that
\(1\le u_e\le\log_2 X\). For every prime \(p\),
\[
 v_p(g_e)\le v_p(L_e)+v_p(u_e),
\]
with equality if \(p\notin S_e\). This follows by taking
the minimum of
\(v_p(L_e)+v_p(v_q(L_e))-\mathbf1_{p=q}\)
over \(q\in S_e\). If the edges divisible by \(p\) in the
original labels form a forest, extend that forest to a
spanning tree in (DM1). It follows that
\(v_p(\Delta(A)/I_B)\le\sum_e v_p(u_e)\). Consequently the
part of \(\Delta(A)/I_B\) on all such primes divides
\(\prod_e u_e\le(\log_2 X)^E\). This includes every
original private prime simultaneously, as well as primes
shared along a forest. The bound concerns the gcd of all
maximal minors, independently of the sizes of the private
prime powers.

Suppose now the box is in a fixed three-packet tuple with
exact full divisors \(F,G,M\) and product \(D=FGM\).
For a prime of one of the packets, the original label
valuations equal its packet order on the entire bipartite
subgraph of that role and are zero on the other two roles.
This follows from exact fullness and original primitivity.
That subgraph has cycle rank \((m-1)^2\). If the product
of the full prime powers on shared primes outside \(D\)
in each complete label is at most \(X^\rho\), then
\[
         \Delta(A)/I_B
             \le D^{(m-1)^2}X^{E\rho}(\log_2 X)^E .
 \tag{DM3}
\]
To check this, for each packet prime use a spanning tree
containing a maximal forest of its role subgraph.
The contribution from original valuations on the
remaining edges is exactly its packet order times
\((m-1)^2\). For other shared primes bound the original
valuations by their sum over the labels. Private primes
were handled above. Each prime also contributes at most
\(\sum_e v_p(u_e)\) from the valuation exponents,
which accounts for the logarithmic factor. The residue
block and privacy bounds control \(\rho\), for example
when all shared nonpacket primes have scaled order bounded
by \(S_0\). In that case (PV4) gives
\(\rho\le S_0(2d+\kappa+\sigma)\).
[Shared support left after removing the packet and private primes (PV4)](ABC_SHARED_BOX_PRIVATE_PRIME_SUPPORT.md#simultaneous-privacy-and-the-outside-cap)

## Slack in the derivative inequality on tall private boxes

The private-coordinate model also computes the quotient
dimension \(u\) from (AD1) exactly:
\[
                            u=h-2=3m-2 .
 \tag{DM4}
\]
The map from prime coordinates to the logarithmic
derivatives of the \(E\) labels is onto, because the
edge support sets are disjoint and nonempty. On its
additive kernel the actual derivatives are in the cut
space of \(G\), of dimension \(h-1\), since they pair
trivially with every triangle circulation. The subspace
where the logarithmic derivatives agree within each
original row maps to the single dimension of the common
value. To see the agreement is global, all edges in a
triangle agree, and the graph of triangles joined by a
common edge connects all edges. Explicitly one can change
one index at a time in the complete array; consecutive
triangles then share an edge. Subtracting the two image
dimensions gives (DM4).

For a numerical comparison, suppose all original labels
are at least \(X^{1-\gamma}\), and their radicals are
at most \(X^{\ell L}\) for fixed \(\ell\).
Each \(w_e^{1/2}\) is at least
\(X^{1-\gamma-\ell L}\) by using any selected prime
\(p\le X^{\ell L}\) in \(a_{e,p}\).
It follows from (DM1) and (DM3) that
\[
 \log_X\frac{\sqrt{\det(AA^{\mathsf T})}}{\Delta(A)}
 \ \ge\ r(1-\gamma-\ell L)
      -(m-1)^2\log_X D-E\rho-o_X(1).
\]
The logarithm of the left side of (AD1) is at most
\(u+E\ell L\), because its height threshold is at most
\(X\) and the selected primes are private to their labels.
For fixed \(m\) and \(D\le X\) the resulting slack is
at least
\[
     2(m-1)^2-r\gamma-(E+r)\ell L-E\rho-o_X(1).
 \tag{DM5}
\]
Thus, in the tall sector with small shared nonpacket
numerical weight and sufficiently small fixed \(L\),
the arithmetic derivative inequality using these private
coordinates holds with fixed positive slack even when
all eligible coordinate flags are used. A surviving
packet cell with more than one row has \(D\le X\)
because the difference in at least one complete letter
is a nonzero multiple of \(D\).

This describes the private-coordinate specialization
of the derivative inequality. A saving from a different
selection of prime coordinates or extra arithmetic
must use information beyond those private supports.


## Geometric mean heights and general private weight

Write
\[
 {\cal L}=\prod_e L_e,\qquad
 {\cal L}_{\rm priv}
   =\prod_e\prod_{\substack{p\mid L_e\\p\nmid L_f\ (f\ne e)}}
              p^{v_p(L_e)},\qquad
 a_m=(1-1/m)^2 .
\]
The private part here includes every private prime power in
the original box, regardless of whether its prime was
selected as a derivative coordinate.

**DM6.** For the private-coordinate matrix of (DM1),
\[
           \Delta(A)/I_B\le
   (\log_2 X)^E ({\cal L}/{\cal L}_{\rm priv})^{a_m}.
 \tag{DM6}
\]
This uses neither a common packet tuple nor a bound on
the scaled valuation orders.

For the proof fix \(p\) and apply the maximal-forest
formula to the edge weights \(v_p(L_e)\).
At level \(j\ge1\), the edges of depth at least \(j\)
join potentials with equal residue modulo \(p^j\).
One such residue class cannot meet all three axes:
choosing one potential from each axis would give an
original primitive triangle whose three labels are
divisible by \(p\). A residue class meeting two axes
therefore gives one complete bipartite component
\(K_{k,l}\), with \(1\le k,l\le m\). A class meeting one
axis contributes only isolated vertices. The cycle
rank in a component with \(kl\) edges is
\[
                   (k-1)(l-1)\le a_m kl .
\]
If \(p\) is private, every positive level has at most
a single edge and its cycle rank is zero. Sum over
the levels of a shared prime to obtain the upper bound
\(a_m\sum_e v_p(L_e)\) for the minimum cotree weight.
Passing to the \(v_p(g_e)\) costs at most
\(\sum_e v_p(u_e)\): use a tree minimizing the original
cotree weight and the pointwise inequality
\(v_p(g_e)\le v_p(L_e)+v_p(u_e)\).
Now multiply over primes and use
\(\prod_e u_e\le(\log_2 X)^E\). This proves (DM6).

Choose one selected private coordinate \(p_e\in S_e\)
on each edge. Define
\[
 \bar h=\frac1E\log_X{\cal L},\quad
 \bar p=\frac1E\sum_e\log_X p_e,\quad
 V_{\rm priv}=\log_X{\cal L}_{\rm priv},\quad
 B_{\rm raw}=\log_X\operatorname{rad}({\cal L}).
\]
Let \(\lambda_{\rm row}=\log_X R\), where \(R\) is the
maximum complete row radical in (AD1). Then the logarithm
of its quotient of right to left sides is at least
\[
  a_m V_{\rm priv}
      -u(1-\bar h-\lambda_{\rm row})
      -r\bar p-B_{\rm raw}-o_X(1).
 \tag{DM7}
\]
The same lower bound holds for every derivative coordinate
set containing the chosen private coordinate on every edge,
even if it also includes nonprivate primes, and for every
eligible set of prime flags in (AD1).

Indeed \(w_e^{1/2}\ge L_e/p_e\).
A uniform spanning tree of the equal-size complete
tripartite graph has an edge in its complement with
probability \(r/E\), the same for all edges by symmetry.
Average the logarithm of the corresponding cotree
product. The maximum cotree product is at least its
geometric mean
\[
        \max_T\prod_{e\notin T} (L_e/p_e)
              \ \ge\ \left(\prod_e L_e/p_e\right)^{r/E}.
\]
Use this in the Gram formula (DM1) and divide by (DM6).
Since \(r-a_m E=3m-2=u\), the logarithm of the
right side of (AD1) is at least
\(u\bar h+a_m V_{\rm priv}-r\bar p-o_X(1)\).
Its left side has logarithm at most
\(u(1-\lambda_{\rm row})+B_{\rm raw}+o_X(1)\):
the first term uses the exact \(H\) of (AD1), and
every possible prime flag is supported on
\(\operatorname{rad}({\cal L})\).

For the assertion with additional primes, use the same
independent triangle basis. Adding matrix columns
adds a positive semidefinite matrix to \(AA^{\mathsf T}\)
and can only decrease the gcd of the maximal minors.
The rank remains \(r\) because the chosen private
coordinates already make the logarithmic-derivative
map onto the edge space. The quotient dimension
\(u=h-2\) is unchanged for the same reason. The
flag bound just used already allows every prime in
the original box. It suffices in this argument to form
(DM1)–(DM6) with just the chosen private coordinate
on each edge. A coordinate set lacking such a
private coordinate on some edge is outside this
comparison, even if its own matrix happens to have
full row rank.

For instance, suppose each label has radical at most
\(X^{CL}\), and its private full-power part is at least
\(X^\omega\) for a fixed \(\omega>0\). Then
\(V_{\rm priv}\ge E\omega\), \(\bar p\le CL\), and
\(B_{\rm raw}\le ECL\). Even without a lower bound on
\(\bar h\) beyond nonnegativity, (DM7) gives slack
\[
             a_m E\omega-u-(r+E)CL-o_X(1).
\]
First choose fixed \(m\) with \(3(m-1)^2\omega>3m-2\),
and then take \(L\) sufficiently small. In this range,
every matrix carrying a private coordinate per label
satisfies (AD1) with a fixed positive numerical margin.


## Private coordinate sets with inactive labels

The exact formulas also apply when the derivative coordinates
are all private but some complete labels have no selected
prime. Let \(J\) be the set of edges with at least one selected
private coordinate. Call the other edges inactive. Contract
each connected component of the inactive spanning subgraph of
\(G\), including its isolated vertices. This gives a connected
multigraph \(G_J\) with \(c_J\) vertices and edge set \(J\).
Parallel edges and loops are retained. For \(e\in J\) keep
\(a_{e,p},g_e,w_e\) as in (DM1).

**DM8.** The rank of the derivative matrix on these coordinates
is \(r_J=|J|-c_J+1\). The formulas (DM1) and (DM2) hold with
\(G_J\) in place of \(G\), and \(I_B\) now the index of the
selected projected triangle basis in the integral cycle lattice
of \(G_J\). In (AD1) the quotient dimension is
\[
  u_J=
  \begin{cases}
     c_J-2,& J=E(G),\\
     c_J-1,& J\ne E(G).
  \end{cases}
\]
A spanning tree of \(G_J\) excludes every loop. Thus each loop
\(e\) contributes a factor \(\sqrt{w_e}/g_e\) to the quotient
\(\sqrt{\det(AA^{\mathsf T})}/\Delta(A)\). If there are inactive
labels, this factor alone bounds the product of all prime flags
in (AD1) selected from \(S_e\).

The restrictions of the original integral circulations to \(J\)
form exactly the integral circulation lattice of \(G_J\).
Contraction gives the inclusion in one direction. In the other
direction, lift an integer circulation on \(G_J\) to its
weights on the original active edges. The induced divergence
has sum zero within each inactive component, so a spanning
tree of that component can route the compensating flows
integrally. This gives an original integral circulation.
As the original triangles generate the original cycle lattice,
their projections generate the contracted lattice. The
derivative column formula therefore gives rank \(r_J\) and the
same fundamental-cycle proof as in (DM1). The level forest
formula for (DM2) works for multigraphs, with loops contributing
one to the cycle rank at every level on which they occur.

The private prime coordinate map onto the logarithmic
derivatives of the active labels is surjective. On its
additive kernel, after multiplying by label lengths, its
image is the cut space of \(G_J\), of dimension \(c_J-1\).
Equivalently an additive derivative potential is constant
across each inactive edge. In the equal-logarithmic-derivative
subspace \(Z\), all labels have the same logarithmic derivative,
as in (DM4). If any label is inactive, that common value is
zero. This proves the quotient dimension formula.

A loop is always in the complement of a spanning tree, giving
the stated factor. For an inactive-label selection, the image
of \(Z\) on the active logarithmic derivatives is zero,
and its prime-coordinate description is the direct sum of
the single-edge kernels. Thus the independent flag primes
\(I_e\) on edge \(e\) leave at least one coordinate
\(q\in S_e\setminus I_e\). For every \(p\in I_e\),
eligibility gives \(p\nmid v_p(L_e)\), since \(p\) occurs
on no other complete label. It follows that
\[
       v_p(g_e)\le v_p(a_{e,p})=v_p(L_e)-1,\qquad
       v_p(a_{e,q})\ge v_p(L_e).
\]
Hence \(\prod_{p\in I_e}p\) divides \(a_{e,q}/g_e\) and
is at most \(\sqrt{w_e}/g_e\). In particular, if the inactive
subgraph is connected, the contracted graph has only loops.
Here \(u_J=0\) and the loop factors account for every
possible flag in (AD1).

A private prime supported on a nonloop active edge can still
be omitted from a cotree product by taking a tree containing
that edge. On a loop it can instead contribute a large
valuation to the gcd. This is the extra possibility that
arises when the derivative coordinate set omits original
labels. The contraction argument uses private coordinates;
additional shared-prime columns need a separate analysis.
