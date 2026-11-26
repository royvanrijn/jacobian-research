# Shared-value boxes and power-class cycles in primitive abc graphs

Three simultaneous full-divisor records give a dense system
of actual difference arrays, with one complete letter at each
pair of indices and a primitive additive row at every triple
of indices. The box count attains the exponent allowed by a
spanning tree. On most of the boxes, the rational products
defined by all cycles simultaneously have large power-class
defect, even when the order and the cycle weights can vary.

## A finite difference-array theorem

Let \(F,G,M\) be pairwise coprime positive integers and let
\(\mathcal E\) be a set of distinct positive primitive rows
\(A+B=C\le X\), with \(F\mid A,\ G\mid B,\ M\mid C\).
Use the absolute residues
\[
           i=B\bmod F,\qquad j=A\bmod G,\qquad k=A\bmod M.
\]
Suppose these coordinates take values in alphabets of sizes
at most \(R_F,R_G,R_M\), respectively. For each role and
each pair of its opposite residue indices suppose there are
at most \(K\) complete values in the whole family. Thus the
\(A\)-values for given \((j,k)\) number at most \(K\), the
\(B\)-values for given \((i,k)\) number at most \(K\), and
the \(C\)-values for given \((i,j)\) number at most \(K\).
Take integer bounds \(R_F,R_G,R_M,K\ge1\), and use the
actual normalized density
\[
                 \delta=\frac{|\mathcal E|}{R_F R_G R_M K^2},
                 \qquad 0<\delta\le1 .
\]
Write \(R=R_F R_G R_M\). An ordered box of size \(m\)
consists of arrays \(A_{jk},B_{ik},C_{ij}\), indexed by
\(1\le i,j,k\le m\), such that every
\((A_{jk},B_{ik},C_{ij})\) is in \(\mathcal E\).
Its three lists of residues are required to be distinct
within each axis. Its labels in one role are therefore
distinct between different pairs of indices in that role.

**LB1.** For every integer \(m\ge2\) there are \(c_m,C_m>0\)
such that, if \(\min(R_F,R_G,R_M)\ge C_m\delta^{-C_m}\),
the number of ordered boxes is between
\[
         c_m\delta^{m(m+1)^2}R^m K^{3m-1}
              \quad\hbox{and}\quad \delta^m R^m K^{3m-1}.             \tag{LB1}
\]
In the lower bound it is the original rows that occupy the
boxes. The same complete value is used on every cell sharing its two indices.

### Common-target tuples and their difference vectors

For each fixed \(C\), write \(\deg(C)\) for the number of its
rows. At most \(K\) of these rows can have the same \(k\)
residue, since \(C\) fixes \(j\) and each \(A\) then
determines the row. Write \(E=|\mathcal E|\).
Discard centers of degree below \(E/(2R_F R_G K)\).
There are at most \(R_F R_G K\) centers, so at least
\(E/2\) rows remain. A retained center has degree at least
\(\delta R_M K/2\). If \(\delta R_M\ge4m\), there are at
least \((\deg(C)/2)^m\) ordered \(m\)-tuples of its rows
having distinct \(k\) residues: at each step at most
\((m-1)K\) choices are forbidden. The total number of
these ordered \(m\)-tuples is at least
\[
  2^{-2m} E^m/(R_F R_G K)^{m-1}
   \ \ge\ \rho R_F R_G R_M^m K^{m+1},\qquad
                 \rho=2^{-2m}\delta^m .                     \tag{LB2}
\]

Make a simple bipartite graph whose left node is the
ordered list of the \(m\) complete \(A\)-values in such
a tuple, and whose right node lists its \(m\) complete
\(B\)-values. An edge specifies the whole tuple of rows,
since \(C\) is determined. For either side, the residues
on the \(j\) axis, respectively on the \(i\) axis, are
constant over the list, because \(C\) is fixed. There
are at most \(R_G(R_M K)^m\) possible left nodes and at
most \(R_F(R_M K)^m\) right nodes. Remove nodes iteratively
if their degrees fall below \(\rho R_F K/4\) on the left
or \(\rho R_G K/4\) on the right. Their total deletions
are at most half the lower bound in (LB2).

An edge connects only lists in the same group, given
by the two pieces of data
\[
     (0,h_2,\ldots,h_m),\quad r=A_1\bmod M,\qquad
                h_a=A_a-A_1=B_1-B_a.
\]
The right node reads \(r=-B_1\bmod M\). At a fixed group a
left node is determined by \(A_1\), a right node by \(B_1\).
Color the left node by \(j=A_1\bmod G\) and the right
node by \(i=B_1\bmod F\). Each color class has at most
\(K\) nodes, because \(r\) is fixed in its group.
Thus the two part sizes are at most \(R_G K,R_F K\).
In a surviving nonempty group the minimum degrees imply
part sizes at least \(\rho R_G K/4,\rho R_F K/4\).
They also imply edge density at least \(\rho/4\).
The edge bound \(R_F R_G K^2\) in each group and (LB2)
give at least \(\rho R_M^m K^{m-1}/2\) surviving groups.

### From common differences to shared values

In a bipartite graph of density \(\theta\), choose \(m\)
independent uniform vertices in each part, with replacement.
The chance that all \(m^2\) edges occur is at least
\(\theta^{m^2}\). First sample one left vertex and \(m\)
right vertices, and apply Jensen to the conditional
edge probability given the left vertex. Then condition
on the \(m\) right vertices and apply Jensen in the
left sample. In a surviving group the chance that
some two samples on one side have the same color is at
most
\[
          \binom m2\frac4\rho(1/R_F+1/R_G).
\]
With the stated lower bound on the \(R\)'s this is at most
\(\frac12(\rho/4)^{m^2}\). Thus each group supplies
\(\gg_m \rho^{m^2+2m}(R_F R_G)^m K^{2m}\) ordered
complete bipartite graphs with distinct colors on both sides.

For every such choice put
\[
        A_{jk}=a_j+h_k,\qquad
        B_{ik}=b_i-h_k,\qquad C_{ij}=a_j+b_i ,
        \qquad h_1=0,
\]
where \(a_j,b_i\) are the first entries of the selected
lists. The group was built using actual rows, so this
supplies all \(m^3\) original edges. The \(k\)-residues
are distinct by the tuple selection, and the \(i,j\)
residues are distinct by the color test. The output
arrays recover the group and the two ordered lists,
so this map is injective. Multiplying by the number
of groups gives a lower bound of order
\(\rho^{(m+1)^2}R^m K^{3m-1}\), as required.

For the upper bound use three groups of potentials
\[
               u_i=b_i,\qquad v_j=-a_j,\qquad w_k=-h_k.
\]
Then \(A_{jk}=w_k-v_j,\ B_{ik}=u_i-w_k,\ C_{ij}=u_i-v_j\).
These potentials are unique up to a common translation
because the tripartite graph of differences is connected.
Fix the ordered list of actual rows in the \(m\)
diagonal cells \((a,a,a)\). This sets all three residue
lists and fixes the lengths of two edges in each of the
\(m\) vertex-disjoint diagonal triangles. Those \(2m\)
edges form a forest. Extend it to a spanning tree of the
potential graph, filling the remaining \(m-1\) edges from
their pair classes, with at most \(K\) choices each.
The resulting lengths recover every potential difference.
There are at most \(E^m\) choices of the diagonal rows,
giving the upper bound \(E^m K^{m-1}=\delta^m R^m K^{3m-1}\).

For any prescribed ordered residue lists, the unrestricted
spanning-tree bound is \(K^{3m-1}\), by enumerating one
label per edge. Thus at least \(c_m\delta^{m(m+1)^2}R^m\)
ordered residue boxes in (LB1) admit shared-value lifts.
At fixed \(m\) and positive price margins, the two bounds
in (LB1) make a power saving in the box count equivalent
to a power saving in the original row count within the
marked regime. In the sparse direction the upper bound
uses the actual density; in the dense direction the
hypothesis on the alphabets is satisfied when the tested
density loss is sufficiently small.
[The eight-edge common-difference argument](ABC_BALANCED_LATIN_TRADE_OBSTRUCTION.md#eight-edge-supersaturation-for-actual-additive-values)
· [Residue box counts](ABC_SIMULTANEOUS_RESIDUE_PRODUCTS.md#nondegenerate-residue-boxes)
· [Bipartite complete-graph density inequality](https://doi.org/10.1017/S000497272500019X)

## Cycle power classes on actual boxes

Orient each edge of the complete tripartite graph on the
potentials from its smaller endpoint to its larger
endpoint. Thus the three edge types are \(v\to w,\ w\to u\)
and \(v\to u\), with positive edge label \(L_e\) equal to
the potential difference. By an integer circulation
\(\mathbf b=(b_e)\) mean an integer vector with sum of its
incoming weights equal to sum of its outgoing weights at
each vertex. For a prime pool \(P\) and \(q\ge2\), use
\[
    \mathcal D_{q,P}(z)
        =\prod_{\substack{\ell\notin P\\ q\nmid v_\ell(z)}}\ell
                  \qquad(z\in\mathbb Q_{>0}) .
\]
For example, an alternating ratio of four labels on an
oriented rectangle is a circulation product.

**LB3.** Fix \(m\ge2\), \(\varepsilon>0\), and a pool \(P\)
whose prime product is at most \(X^B\) for fixed \(B\).
Uniformly in \(Z\ge1\), the number of ordered boxes,
in the setup of (LB1), with
\[
        \mathcal D_{q,P}\biggl(\prod_e L_e^{b_e}\biggr)\le Z
\]
for some \(q\ge2\) and some circulation with an entry
nonzero modulo \(q\) is at most
\[
 O_{m,B,\varepsilon}\left(
    Z X^\varepsilon R^m K^{3m-2}
              (R_F^{-1}+R_G^{-1}+R_M^{-1})\right).            \tag{LB3}
\]
This estimate needs only the alphabet and class upper
bounds. The same estimate applies to arbitrary lists of
integer exponents whose reduction after dividing by a
common prime power of the order is a unit on at least
two edges incident to some vertex, as made explicit below.

### Recovering one vertex in a weighted power relation

Choose a prime \(p\mid q\) with some \(b_e\) nonzero modulo
\(p^{v_p(q)}\), and put \(h=\min_e v_p(b_e)<v_p(q)\),
using \(+\infty\) for an entry zero. Divide all weights by
\(p^h\) and consider the resulting product modulo \(p\)-th
powers. Its defect support outside \(P\) is contained in
the original defect support. Indeed at each prime outside
that support the sum of valuations with the old weights
was divisible by \(p^{v_p(q)}\). At least one new weight
is a unit modulo \(p\). Since the new weights are still a
circulation, at either endpoint of an edge of unit weight
another incident weight must also be a unit modulo \(p\).
Choose that vertex and those two edges. There are only
\(O_m(1)\) choices of this combinatorial data.

Fix all other potentials up to common translation. The
edges incident to the remaining potential \(x\) give
positive affine forms of slope \(1\) or \(-1\).
Their roots are the fixed neighboring potentials.
They are distinct. Within one other group they are
distinct by the residue tests, and any neighbors in
different groups are separated by a positive edge in
the box. Enlarge the prescribed pool by the support of
all fixed edge labels and all nonzero differences
between neighboring potentials. Its prime product is
bounded by \(X^{O_m(1)+B}\). Indeed all potential
differences are \(O_m(X)\) by connectedness.
Outside this pool each prime divides at most one of
the forms. Consequently every prime of either selected
form outside the enlarged pool and the defect support
occurs to an exponent divisible by \(p\). This assertion
uses only the two unit weights at the chosen vertex and
holds uniformly in the values of all weights.

Enumerate the outside defect support by a squarefree
integer at most \(Z\). The two canonical \(p\)-free
coefficients can be supported only on its union with the
enlarged pool. Their values are at most \(X\) and cost
\(X^{o_X(1)}\) by uniform supported-value counting. For
\(p\le\log_2 X\) there are at most \(\log_2 X\) orders
to enumerate. At fixed coefficients the two affine
forms become \(A y^p\) and \(B_1 z^p\).
Eliminating \(x\) gives a binomial of fixed nonzero
value, because the two affine roots differ. For orders
bounded in terms of \(\varepsilon\) the uniform
fixed-value binary norm count gives \(O_\varepsilon(X^\varepsilon)\)
choices, with a decreased allowance if necessary.
For all larger orders, one of the positive roots has
at most \(X^{1/p}\) choices and determines \(x\).
For \(p>\log_2 X\), both forms themselves have their
whole support on the specified union. Enumerate their
bounded values directly, without any cost for the
orders. In total the possible values of \(x\), for all
orders and weights, number \(O_{m,B,\varepsilon}(Z X^\varepsilon)\).
[Uniform fixed-value binary counts](ABC_RECTANGLE_POWER_CLASS_DEFECT.md#the-fixed-value-bound-with-moving-coefficients)
· [Uniform supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For an omitted vertex in axis \(i\), choose the ordered
residues of the other vertices at cost at most \(R^m/R_i\).
The tripartite difference graph after deleting one vertex
is connected, since \(m\ge2\). A spanning tree on its
remaining vertices determines their potentials up to
translation, at cost at most \(K^{3m-2}\).
Apply the preceding bound to each anchor and sum over
the choices of vertex and incident edges. This proves
(LB3), without counting the possibly unbounded weights.
The same proof permits a pool chosen separately from
each such anchor, as long as the pool used on the
corresponding omitted vertex has the stated height bound.

## Consequences for simultaneous packet strata

In a surviving fixed tuple of the three-role case of
(SB2), use \(n=X^L\), prices \(s_F,s_G,s_M>0\), and put
\(d=s_F+s_G+s_M\), \(t=1/6-d>0\). At fixed positive
margins, the SB2 proof gives
\[
        R_i=n^{s_i+O(\kappa+\eta)},\qquad
        K=n^{t+O(\eta)},\qquad
        E\ge n^{d+2t-\eta}.
\]
Here \(R_i\) use the one-role residue cap and \(K\)
uses the three-role ambient class pruning. They can
be chosen as upper bounds even when a further power-small
fraction of rows or their cells have been removed.
For any fixed \(m\), take the tolerances small in terms
of \(m,\min_i s_i,t\) before taking the height limit.
Then (LB1) forces
\[
                 n^{md+(3m-1)t-O_m(\kappa+\eta)}
                                                               \tag{LB4}
\]
ordered actual shared-value arrays, up to a constant.
The same lower bound holds with an additional
\(O_m(\xi)\) in the loss on every subfamily with at least
\(n^{d+2t-\xi}\) edges, for sufficiently small fixed \(\xi\).
Every row of an array has all three original radical caps
and full original packets. Its normalized difference edges
therefore retain their three outside radical caps.

By (LB3) almost all the boxes at the exponent in (LB4)
have the following stronger property. For any fixed
\(0<\sigma<t+\min_i s_i\), every circulation
nonzero modulo its order \(q\ge2\) has defect outside
a prescribed polynomial-height pool greater than
\(n^{t+\min_i s_i-\sigma}\).
Choose the fixed tolerances and the normalized inner
allowance \(\varepsilon/L\) small in terms of \(\sigma,m\).
The bad boxes are power fewer than the lower bound (LB4).
In particular, for every prime order \(p\) the map
\[
  \mathbf b\ \longmapsto\
       \left(\sum_e b_e v_\ell(L_e)\bmod p\right)_{\ell\notin P}
\]
on the circulation space over \(\mathbb F_p\) is injective
on those boxes. Every circulation over \(\mathbb F_p\)
lifts to an integer one by using the fundamental cycles
of a spanning tree. Its image thus has dimension
\(3m^2-3m+1\) simultaneously at every prime order.
This includes independent power-class information
around all cycles, on the same actual residue box.

The reduction now gives a test using genuine integer
arrays of radical-constrained differences. A power upper
bound below the scale in (LB4) for those arrays, or a
deduction of cheap cycle defects on enough of them, would
rule out their packet strata. The existence of the arrays
and the large cycle defects give structural requirements
for saturation. They do not supply the upper bound.
[Simultaneous full packets and the shortened abc lattice](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md#three-marked-roles)
