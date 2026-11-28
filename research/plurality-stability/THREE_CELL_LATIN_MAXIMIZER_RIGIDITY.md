# Three-way selectors and pointwise-distinct Gaussian maximizers

Fix \(0<\rho<1\) and let \(s_*\) denote the balanced three-cell
Gaussian optimum. Suppose a centered regular simplex partition
\(T\) attains \(s_*\). A three-way selector gives a necessary
comparison for any collection of three maximizers. The geometric
equality case in which the three label assignments are pointwise
distinct has a rigidity property when one of them is \(T\).

## The selector inequality

Write \(Q(f)=\sum_i\langle f_i,T_\rho f_i\rangle\) for a vector
assignment. Let \(F^1,F^2,F^3\) be balanced maximizers in a common
Gaussian space. At each point \(x\), permute these assignments by
any measurable permutation to obtain \(G_1(x),G_2(x),G_3(x)\).
The pointwise mean \(U=(F^1+F^2+F^3)/3\) is unchanged. In an
independent Gaussian plane use the simplex label \(J\) as the
selector, so the constructed partition is \(G_{J(z)}(x)\).
It is ordinary and balanced since its conditional mean given
\(x\) is \(U(x)\).

The simplex transition matrix has equal diagonal entries and
equal off-diagonal entries. Put \(\lambda=(3s_*-1)/2>0\).
Its conditional channel is \(\lambda I+(1-\lambda)\mathbf1
\mathbf1^{\mathsf T}/3\). Conditioning on the two auxiliary
labels and using independence of the Gaussian blocks gives
\[
      (1-\lambda)Q(U)+\frac{\lambda}{3}\sum_{j=1}^3 Q(G_j)
                \le s_* .                                 \tag{LR1}
\]
Equality forces the selector partition to depend only on the
base or only on the auxiliary block, by the
[strict finite recombination comparison](THREE_CELL_PRODUCT_RECOMBINATION.md#finite-recombinations-using-two-blocks-lose-strictly).
The base input for that comparison consists of the finite
tuple of all three labels and its permutation. Base-only
dependence occurs when all three original assignments coincide.
Auxiliary-only dependence occurs when each selected row is a
constant label. Balance then makes the constants distinct, so
the three original assignments are pointwise distinct up to
null sets.

For a concrete necessary comparison, take a maximizer \(A\)
and the two cyclic relabelings of \(T\). With \(d=1_A-1_T\),
the three available labels at \(x\) are all labels once, except
that the label of \(T\) is replaced by the label of \(A\).
Thus the row with default label \(i\) can be chosen as
\(G_i=e_i+1_{T_i}d\). Put \(m=\gamma(A\ne T)\). Since \(d\)
has zero mean, \(Q(U)=1/3+Q(d)/9\), and expanding the row
energies gives \(\sum_iQ(G_i)=3-2m+\sum_iQ(1_{T_i}d)\).
If \(A\) and \(T\) differ on a positive-measure set, then
\[
     \frac{1-\lambda}{3\lambda}Q(d)+\sum_i Q(1_{T_i}d)<2m .
                                                               \tag{LR2}
\]
The strict inequality follows because the selector always uses
at least two outputs on every fiber, and each row has constant
label \(i\) outside \(T_i\), so base independence would give
\(A=T\). Taking \(T\) and the two cyclic relabelings of \(A\)
instead gives the same bound with \(1_{A_i}d\) in place of
\(1_{T_i}d\).

## A pointwise-distinct triple lies in one active plane

Suppose \(A,B,T\) are ordinary balanced maximizers and their
three labels are distinct at almost every point. For each label
\(i\), this says \(1_{A_i}+1_{B_i}+1_{T_i}=1\). Let their
two-dimensional active planes be \(V_A,V_B,V_T\), as supplied
by the [dimension reduction and full rank results](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible).
For \(v\perp V_T\) and \(w\perp V_B\), the pointwise sum gives
\(D_vD_w1_{A_i}=0\) in distributions. Restrict to \(V_A\), where
these are derivatives in the projected directions.

If the two projections are parallel and nonzero, the second
derivative along their common direction vanishes. Convolving
with a smooth compactly supported mollifier gives bounded
functions affine on each line in that direction. They are
constant along such lines, and the mollification limit is
cylindrical. If the projections are independent, use them as
coordinate axes \(s,t\). The mixed derivative equation gives
\(1_{A_i}(s,t)=p(s)+q(t)\) almost everywhere. Indeed the
mollifications satisfy the zero rectangular difference identity,
which passes to the limit in local \(L^1\) in the four rectangle
coordinates. Fubini then gives the representation by fixing two
generic base coordinates. If \(p\) takes two different values
at generic \(s\), the two binary values \(p(s)+q(t)\) differ by
a fixed nonzero number for almost every \(t\), so \(q\) is
constant. Otherwise \(p\) is constant. The function is again
cylindrical, contrary to the fact that
[each maximizing cell uses both active coordinates](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions).

Thus for every such \(v,w\) at least one projection vanishes.
If \(V_A\ne V_T\), choose \(v\) with nonzero projection and
conclude \(V_A\subseteq V_B\). Their dimensions give equality.
The pointwise sum would then also put \(T\) in this plane,
contradicting \(V_A\ne V_T\). Consequently \(V_A=V_T\), and the
pointwise sum gives \(V_B=V_T\) as well.

## Angular rigidity at the origin

Work in the common active plane. Write \(V_i^A,V_i^B,V_i^T\) for
the analytic scores with their stationary prices subtracted.
No label can be a strict winner of \(A\) or \(B\) at the origin,
since every label of \(T\) occurs arbitrarily close by.
We use the [local junction analysis and reciprocal-gradient bound (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions).
At a genuine triple the gradients are pairwise distinct. With
noncollinear gradients their linearization assigns each label
an angular interval of positive length less than \(\pi\), by
intersection of the two independent winning half-planes. With
collinear gradients only the two extreme gradients win on
angular intervals of positive length. All angular statements
refer to open intervals away from the finitely many tangent
directions, where strict first-order inequalities persist on
neighborhoods of sufficiently small points.

Suppose that label \(k\) of \(A\) does not approach the origin.
In a small disk only labels \(i,j\) of \(A\) occur, where
\(\{i,j,k\}=\{1,2,3\}\). The \(T_i,T_j\) sectors force \(A\) to
choose \(j,i\), respectively. This gives an actual \(A_i,A_j\)
interface on their common ray up to the origin. Its score
gradient is nonzero at the origin by (G20), and its zero set
there continues on the same line by analyticity. The two
labels of \(A\) therefore win on opposite sides of this line
near the origin. The line enters \(T_k\) along its bisector,
so \(B\) has both of its labels \(i,j\) on positive angular
intervals in \(T_k\). On the two other sectors \(B\) must use
\(k\), a union with angular length \(4\pi/3\). Thus \(B\) has a
genuine triple with positive angles for all three labels, and
one angle exceeding \(\pi\), which neither gradient case allows.

The same contradiction applies if \(A\) has a collinear
genuine triple. In this case its two extreme gradients \(i,j\)
win in opposite open half-planes of the linearized scores.
The middle label \(k\) has no positive angular interval.
The interiors of \(T_i,T_j\) must fall in the respective
half-planes for \(j,i\). Therefore their adjacent ray is on
the dividing line, which bisects \(T_k\) in the other direction.
This again forces the forbidden angular pattern for \(B\).
Exchanging \(A,B\) proves that both have noncollinear genuine
triples at the origin.

Let \(I_i\) be the circular interval of directions for \(T_i\).
For each \(i\), the corresponding angular intervals of \(A_i\)
and \(B_i\) partition the complementary interval of length
\(4\pi/3\), apart from endpoint directions. Each therefore
touches exactly one endpoint of \(I_i\). The three endpoint
choices for \(A\) are distinct. Otherwise rotate and relabel
the circular coordinate so the shared endpoint is \(2\pi/3\)
between \(I_1=(0,2\pi/3)\) and \(I_2=(2\pi/3,4\pi/3)\).
The interval of \(A_1\) starts just beyond \(2\pi/3\), and
that of \(A_2\) ends there. If \(A_3\) chooses the endpoint
\(4\pi/3\) of \(I_3\), it ends there from below, between \(A_1\)
and the start of \(A_2\). Then \(A_2\) occupies the interval
from \(4\pi/3\) to \(2\pi/3\) through zero, of length \(4\pi/3\).
If \(A_3\) chooses endpoint zero instead, it starts there and
ends between zero and the end of \(A_2\). Then \(A_1\) occupies
the interval from \(2\pi/3\) to zero of the same length.
Either choice contradicts the angle bound. Thus the three
endpoints of \(A\) are exactly the three rays of \(T\), and its
labels are a cyclic relabeling. The pointwise-distinct property
gives the other cyclic relabeling for \(B\) in the tangent sectors.

## Continuation inside the simplex sectors

Normalize each row of prices to have sum zero, taking the
simplex prices equal to zero. All three adjusted scores of
each partition now tie at the origin, because the origin is a
genuine triple. The ordinary pointwise sum gives
\(T_\rho1_{A_i}+T_\rho1_{B_i}+T_\rho1_{T_i}=1\), and the
simplex scores at zero are all \(1/3\). Hence the two price
vectors of \(A,B\) sum to zero, and \(V_i^A+V_i^B+V_i^T=1\).

Inside \(T_k\) the only possible labels of \(A,B\) are the other
two, \(i,j\). Any actual \(i/j\) jump of \(A\) there is the same
jump with reversed labels for \(B\), by the pointwise-distinct
property. Adding the two score equalities shows that
\(V_i^T=V_j^T\) on the jump. The simplex contrast has a strict
sign on each side of its mirror line. This follows by folding
the positive Gaussian kernel across the reflection exchanging
labels \(i,j\). Inside \(T_k\) its zero set is the bisector ray
of that sector. Hence all jump points inside \(T_k\) lie on
that ray.

The jump closure is relatively closed on the interior ray.
It is also relatively open. Local finiteness of the analytic
stratification gives an incident regular jump arc at each
closure point. By (G20) the \(i/j\) score difference has nonzero
gradient there, so the implicit function theorem continues its
zero curve through the point. On both sides of this curve the
sign of that difference determines the label almost everywhere,
since label \(k\) is absent from \(T_k\). Thus the whole local
curve is a sign jump, including at any dormant third-score tie.
The ray is connected,
so a nonempty jump closure would be the whole ray. But on its
initial segment the cyclic tangent sectors give a strict winner.
Therefore \(A\) has no jump inside \(T_k\). Its label is constant
almost everywhere there, since the connected region has no
distributional boundary between the two allowed labels. This
holds in all three sectors and completes the classification of
\(A,B\) as the two cyclic relabelings of \(T\).

## Uniform defect and the scope of the third-maximizer condition

Fix a compact subinterval of \((0,1)\). There is a positive
lower bound for the Gaussian measure of the event that the
three labels are not distinct, uniformly over \(T\) attaining
the optimum, nonsimplex balanced maximizers \(A\), and balanced
maximizers \(B\) on this interval and in any ambient dimension.
Otherwise the three active planes span at most six coordinates,
so all partitions may be represented in one fixed \(\mathbb R^6\).
Rotate to fix the labelled simplex \(T\) and pass to convergent
correlations. Strong maximizer compactness gives a limit of each
of \(A,B\), and the vanishing defect passes to the limit. The
preceding result identifies the limit \(A\) with a relabeling
of \(T\). In particular its moment vectors span its active plane,
and convergence of moment vectors shows that the active planes
of the sequence converge to this plane. Align them by orthogonal
maps tending to the identity. The
[uniform measurable local comparison](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)
with the all-correlation three-cell gap (S37) excludes nonsimplex
maximizers approaching that orbit, a contradiction.

If \(A,T\) are pointwise disjoint maximizing assignments, their
ordinary complementary partition has indicator vector
\(1_B=\mathbf1-1_A-1_T\) and is balanced. Expansion of its energy
gives \(S_\rho(B)=4s_*-1-Q(1_A-1_T)\). Therefore
\[
                 Q(1_A-1_T)\ge 3s_*-1,                     \tag{LR3}
\]
with equality exactly when \(A\) is a cyclic relabeling of the
labelled simplex \(T\). Indeed equality would make \(B\) an
actual maximizer and invoke the pointwise-distinct rigidity.
For nonsimplex maximizers in a compact positive band the strict
gap in (LR3) is uniform. If it tended to zero, after embedding
the two active planes in four dimensions the complementary
partitions would approach the optimum as well. Strong
compactness and local isolation would give the same
contradiction as for the uniform defect.

A simplex \(T\) and its central
inversion \(-T\) choose the largest and smallest of the three
linear simplex scores, respectively. Their unique complementary
partition chooses the middle score. It is balanced by cyclic
rotation symmetry and invariant under central inversion, so all
its first moments vanish. Full active rank excludes it from the
maximizers. This gives disjoint optimizing assignments with strict
inequality in (LR3) even when both are simplices.
