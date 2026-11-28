# Orthogonal product groups cannot occur in Gaussian noise maximizers

Fix \(0<\rho<1\) and positive prescribed Gaussian cell masses.
Let \(A=(A_1,\ldots,A_k)\) be a maximizing partition for
\[
 S_\rho(A)=\sum_i\langle 1_{A_i},T_\rho 1_{A_i}\rangle_\gamma .
\]
Suppose \(\mathbb R^d=W\oplus Y\) orthogonally and a group \(I\)
of at least two labels has
\[
 A_i=B_i\times C\quad(i\in I),                                \tag{OP1}
\]
where the \(B_i\) partition \(W\) into positive-mass cells
and \(C\subset Y\) has mass strictly between zero and one.
Equalities of cells are understood modulo null sets.
This is impossible if \(d\ge k-1\) or \(\dim W\ge2\).
The cells with labels outside \(I\) can depend on both
coordinate blocks.

In particular, arbitrary positive-mass maximizing partitions
in dimension \(d\ge k-1\) cannot be Cartesian products of two
partitions with at least two labels each. The same conclusion
holds if the proposed product allows different fine partitions
over different coarse cells. This last case follows from the
recognition argument below.

## Ties on common essential supports

For a measurable set \(D\), let \(E(D)\) be its closed
essential support, consisting of points where every open
ball meets \(D\) in positive measure. For a finite union
the essential support is the union of these supports.
For products it satisfies \(E(D\times F)=E(D)\times E(F)\).
Gaussian and Lebesgue essential supports agree.

A maximizer has adjusted winning scores
\(U_i-\lambda_i\), where \(U_i=T_\rho 1_{A_i}\) and the
\(\lambda_i\) are real prices. Indeed, positivity of
\(T_\rho\) as a quadratic form implies a maximizing
partition maximizes the linear functional given by
its own scores among assignments of the same masses.
For the prices, minimize
\(\int\max_i(U_i-\lambda_i)\,d\gamma+\sum_i\gamma(A_i)\lambda_i\).
This function is invariant under adding a common constant
to the prices. With \(\min_i\lambda_i=0\), its sublevel
sets are bounded because all cell masses are positive
and the integral is nonnegative. A minimum exists.
By analyticity and Gaussian injectivity, pair ties have
zero measure for all prices. Differentiating the dual
function at its minimum shows that the winning cells
have the prescribed masses. The original partition
therefore attains equality in the pointwise assignment
bound and agrees almost everywhere with the winners.
Continuity implies equality of the two adjusted scores
on \(E(A_i)\cap E(A_j)\). This use of the tie
equation needs neither smoothness nor a regular pair
interface.
[Adjusted-score assignment](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)
[Gaussian convolution regularity, Heilman–Tarter, §2](https://arxiv.org/html/2011.05583#S2)

The scores are real analytic. A contrast that is constant
on a positive-measure set in a coordinate block is constant
throughout that block. We will also use this injectivity
consequence: for nontrivial measurable sets \(D,F\)
and positive coefficients \(a,b\),
\[
 aT_\rho 1_D-bT_\rho 1_F=\text{constant}
 \quad\Longrightarrow\quad D=F\text{ a.e. and }a=b.           \tag{OP2}
\]
Indeed, \(T_\rho\) is injective on Gaussian \(L^2\) by its
Hermite expansion, so \(a1_D-b1_F\) is constant. If
\(D\setminus F\) has positive measure the constant must be
\(a\), which would force \(D\) to have full measure. Exchanging
the two sets gives equality, and their nontriviality then
gives \(a=b\).

## Rotating one group preserves stability exactly

For any \(R\in O(W)\), replace \(A_i\) by
\(A_i^R=(RB_i)\times C\) for \(i\in I\), leaving
all other labels fixed. Every cell keeps its mass.
Each summand of \(S_\rho\) with index in \(I\) is
the product of the individual stabilities of its two
factors. Rotating its \(W\)-factor preserves both of
these numbers. The rotated partition is therefore
another maximizer with prices \(\lambda_i^R\).

There is \(y_0\in E(C)\cap E(C^c)\) since these
two nonempty closed sets cover the connected space \(Y\).
For any \(i,h\in I\) with \((RB_i)\cap B_h\)
of positive measure, consider the points \((x,y_0)\)
where \(x\in E(RB_i)\cap E(B_h)\). This set of \(x\)
has positive measure. Each \((x,y_0)\) belongs to the essential
support of some outside label, since the outside labels fill
\(W\times C^c\). There are finitely many such labels,
so one label \(j\notin I\) accounts for a positive-measure
subset of these \(x\). On this subset the original
and rotated tie equations are
\[
 \begin{aligned}
 T_\rho1_C(y_0)T_\rho1_{B_h}(x)-U_j(x,y_0)
     &=\lambda_h-\lambda_j,\\
 T_\rho1_C(y_0)T_\rho1_{RB_i}(x)-U_j(x,y_0)
     &=\lambda_i^R-\lambda_j^R .
 \end{aligned}
\]
Subtracting eliminates the unchanged outside score.
The resulting inner-score difference is constant on a
positive-measure set because \(T_\rho1_C(y_0)>0\).
Analyticity and (OP2) give \(RB_i=B_h\) almost everywhere.

Thus every orthogonal map of \(W\) permutes the fine sets
modulo null sets. The action of a continuous path of rotations
on any indicator is continuous in Gaussian \(L^1\), by
approximation by bounded continuous functions.
Since the orbit lies in a finite set, every element of
the connected group \(SO(W)\) fixes every \(B_i\).
For \(\dim W\ge2\), a vector fixed by all these rotations
is zero. In particular all Gaussian centroids of the
fine cells vanish.

On the other hand, the fine partition itself maximizes
noise stability at its prescribed fine masses in \(W\).
An improvement in \(W\) would improve \(A\) by the
strictly positive individual stability of \(C\).
The translation argument of Heilman–Tarter rules out
a nontrivial maximizing partition whose cell centroids
all vanish. Specifically, all translations would have
zero linear volume flux, and their Eq. (38) forces the
wall normals to be orthogonal to all translation vectors.
The reduced boundaries would then have zero perimeter,
which forces every fine cell to be null or full.
This proves the exclusion in (OP1).
[Heilman–Tarter, Lemma 5.5 and Eq. (38)](https://arxiv.org/html/2011.05583#S6)

If \(d\ge k-1\) and \(\dim W=1\), extend the original
partition to one more independent Gaussian coordinate.
This extension is a maximizer by dimension reduction:
an improvement in dimension \(d+1\) would have a
maximizer depending on at most \(k-1\) coordinates
and could therefore be realized in dimension \(d\).
The new coordinate may be included in \(W\), giving
a two-dimensional fine block with the same product
group and leaving the outside labels unchanged.
[Heilman–Tarter, Theorem 1.9](https://arxiv.org/html/2011.05583#S1)

## Recognizing a common factor across coarse branches

Suppose each cell instead has the form
\(B_i^{(j)}\times C_j\), with pairwise disjoint label
groups over the coarse positive-mass partition
\((C_j)_{j=1}^q\) of \(Y\), \(q\ge2\). The fine
partitions of \(W\) may initially differ across groups.
Join two coarse indices when their essential supports
intersect. The graph is connected, since the finitely
many nonempty closed essential supports cover \(Y\)
and a disconnected graph would give a disjoint closed
cover by two nonempty sets.

For adjacent \(j,l\), choose
\(y_0\in E(C_j)\cap E(C_l)\). If
\(B_i^{(j)}\) and \(B_h^{(l)}\) overlap in
positive measure, continuity of the full adjusted
scores gives a constant contrast
\[
 T_\rho1_{C_j}(y_0)T_\rho1_{B_i^{(j)}}(x)
   -T_\rho1_{C_l}(y_0)T_\rho1_{B_h^{(l)}}(x)
      =\lambda_{ij}-\lambda_{hl}
\]
on the positive-measure set
\(E(B_i^{(j)})\cap E(B_h^{(l)})\).
If one fine set is full and the other proper,
analyticity and injectivity make this impossible.
If both are proper, (OP2) makes them equal almost
everywhere. Thus complete fine partitions match
across adjacent groups up to relabeling. Connectivity
makes all groups share a common fine partition.
In particular, if some branch has a nontrivial split,
all branches do and \(k=pq\), where \(p\ge2\) is
their common number of fine labels.

For a Cartesian product with both factors nontrivial,
any independent unused coordinate may be included in the
fine block. If \(d\ge3\), either one of the factor blocks
already has dimension at least two or there is such an
unused coordinate. In the former case use that block
for the fine partition, exchanging the two factors if
needed. In either case (OP1) then applies. In particular
it applies in the standard simplex dimension regime
\(d\ge k-1\) since \(k=pq\ge4\).
The argument leaves the product of two one-dimensional
factors in dimension two open. For a single product group
with arbitrary outside cells in a one-dimensional fine
block below dimension \(k-1\), the rotation argument
only forces reflection to permute its fine cells.

## An explicit mixed rotation instability for polyhedral products

The exact rotation argument above uses global maximality.
For a finite polyhedral Cartesian product, an interface
calculation also rules out cylindrical second-variation
stability when the two factors have nonzero weighted
self-centroids. This condition holds for regular simplex
factors.

Let \(B=(B_i)\) and \(C=(C_j)\) be finite polyhedral
partitions of the two blocks. Assume the product partition
is an ordinary adjusted-score winner and give it one
independent Gaussian coordinate. On the regular walls
of the product without this last coordinate, use Gaussian
surface measure \(d\nu\) and the oriented incidence vector
\(d_e=e_a-e_b\) for a wall oriented from cell \(a\) to
cell \(b\). With \(a_e\) the magnitude of the pair score
gradient on wall \(e\), the polarized first-coordinate
form is
\[
 Q^{(1)}(f,g)
 =\sum_e\int_e a_e f_e g_e\,d\nu
 -\rho\sum_{e,h}(d_e\cdot d_h)
       \int_{e\times h}M_\rho(x,x') f_e(x)g_h(x')\,
          d\nu(x)\,d\nu(x').
                                                               \tag{OP3}
\]
Stability implies this form is positive semidefinite.
In the finite polyhedral setting one may first cut
off on compact subsets of regular facets and
multiply by compactly supported odd functions of the
extra coordinate. The resulting linear Gaussian
volume flux is zero. The connected facet-adjacency
graph supplies spanning-tree bumps for exact volume
correction. The bounded score gradients and the Gaussian
decay of the joint density restricted to finitely
many supporting hyperplanes allow the facet cutoffs
to be removed. The odd functions converge in Gaussian
\(L^2\) to the linear coordinate and its Mehler eigenvalue
is \(\rho\), yielding (OP3).
[First extra-coordinate form, G1](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Fix \(v\in W\). Write \(s_j\) for the normal speed
of translation along \(v\) on \(B\)-walls inside
the coarse cell \(C_j\), with zero speed on all
other walls. When \(j\ne l\), the two wall supports
are disjoint and the incidence vectors of the two
families of walls lie in disjoint blocks of label
coordinates. Therefore \(Q^{(1)}(s_j,s_l)=0\).
The full translation \(\sum_j s_j\) satisfies the
translation identity \(a s=\rho K s\), where \(K\)
is the interaction operator in (OP3) without the
factor \(\rho\). Thus
\[
 Q^{(1)}(s_j,s_j)=Q^{(1)}(s_j,\sum_l s_l)=0.
\]
Likewise, for \(w\in Y\) the normal translation speed
\(t_i\) on \(C\)-walls inside \(B_i\) satisfies
\(Q^{(1)}(t_i,t_i)=0\).
[Normal translation identity, G2](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Heilman–Tarter, Lemma 5.4](https://arxiv.org/html/2011.05583#S5)

The mixed term pairs a \(B\)-wall with a \(C\)-wall.
The local term vanishes. A \(B\)-wall oriented
from \((a,j)\) to \((b,j)\) and a \(C\)-wall
oriented from \((i,c)\) to \((i,d)\) have incidence
pairing \((\delta_{ai}-\delta_{bi})(\delta_{jc}-\delta_{jd})\).
The two signs turn the respective normals into the
outward normals of \(B_i,C_j\). The Mehler kernel
factors in the two coordinate blocks. Integrating in
the interior variables therefore gives
\[
 Q^{(1)}(s_j,t_i)
 =-\rho
    \left\langle v,\int_{\partial^*B_i}
             N_{B_i}T_\rho1_{B_i}\,d\nu_W\right\rangle
    \left\langle w,\int_{\partial^*C_j}
             N_{C_j}T_\rho1_{C_j}\,d\nu_Y\right\rangle .
                                                               \tag{OP4}
\]
For either factor cell \(D\), Gaussian divergence
rewrites the wall vector as
\[
 \int_{\partial^*D} N_D T_\rho1_D\,d\nu
 =-\frac{1}{1+\rho}\int_D x T_\rho1_D(x)\,d\gamma(x).
                                                               \tag{OP5}
\]
Indeed divergence gives the integral over \(D\) of
\(\nabla T_\rho1_D-xT_\rho1_D\).
Differentiate the conditional density in its base
variable to obtain the score derivative with multiplier
\(\rho(y-\rho x)/(1-\rho^2)\). Integrating over
\(D\times D\) and using exchangeability makes the
gradient term \(\rho/(1+\rho)\) times the weighted
centroid in (OP5).

If each factor has a nonzero vector in (OP5),
choose the relevant \(i,j,v,w\) so that (OP4) is
nonzero. One of \(s_j+t_i\) and \(s_j-t_i\) has
negative form, proving instability. A regular simplex
cone lies in a closed origin halfspace, with its
interior in the open halfspace. Its weighted centroid
has a strictly positive projection on the inward
halfspace normal because its score is strictly positive.
This verifies the claimed condition for regular simplex
products at every \(0<\rho<1\).
