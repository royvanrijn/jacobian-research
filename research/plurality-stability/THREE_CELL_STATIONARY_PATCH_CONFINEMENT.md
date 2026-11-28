# Finite relabeling supported near regular Gaussian ties

Fix \(0<r<1\). Let \(A=(A_i)_{i=1}^3\) be a planar Gaussian
partition with prescribed positive masses. Assume it satisfies the
winner equations for adjusted scores \(V_i=T_r1_{A_i}-\beta_i\).
Thus for almost every \(x\in A_i\), \(V_i(x)=\max_j V_j(x)\).
Let \(K\) be compact and suppose \(\nabla(V_i-V_j)(x)\ne0\)
whenever \(x\in K\) and \(V_i(x)=V_j(x)=\max_l V_l(x)\).
This condition holds at ordinary two-label interfaces with strict
spectator gap and at genuine triple junctions of an actual maximizer.
At a genuine triple the three pair gradients are nonzero even in
the collinear case.
[Pair-gradient nonvanishing and junction description](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

There are \(C,r_0>0\) such that any finite cover by patches
\(B_b=B(z_b,s_b)\) with \(z_b\in K\) and \(0<s_b\le r_0\)
has the following comparison. Define
\[
  M_{bc}=\sup_{\overline{B_b}\times\overline{B_c}}M_r,\qquad
  G_{bc}=2C\sqrt{s_b s_c}\,M_{bc},
\]
where \(M_r\) is the Gaussian-relative Mehler kernel. If the largest
eigenvalue of \(G\) is less than one, every change confined to
\(\bigcup_b B_b\) that preserves each cell mass strictly decreases
noise stability unless the change is null. The comparison holds for
ordinary partitions and for simplex-valued label assignments.

The score-regret mechanism also gives a small-support comparison at
the Euclidean \(Y\)-cone. Here it applies at arbitrary Gaussian
winner-stationary geometry with the stated regularity, and the matrix
keeps the pair-dependent interactions of separated patches.
[Euclidean \(Y\)-cone comparison](THREE_CELL_NEAR_ONE.md#a-compact-euclidean-junction-test)
[Gaussian score-regret expansion at the simplex](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#global-score-margins-and-facet-normal-derivatives)

## Score sublevels in a ball

Write \(a(x)\) for the label of \(A\) and
\(\eta(x)=\min_{j\ne a(x)}[V_{a(x)}(x)-V_j(x)]\ge0\) almost
everywhere. There is a uniform bound
\[
             \gamma_2(B(z,s)\cap\{\eta\le t\})\le Cst
 \qquad(z\in K,\ 0<s\le r_0,\ t>0).                        \tag{PCF1}
\]
For each pair \(i,j\), put
\(g_{ij}=\max_l V_l-\min(V_i,V_j)\).
The set \(\{\eta\le t\}\) lies in the union of \(\{g_{ij}\le t\}\).
On \(K\), at every point and for each pair one has either
\(g_{ij}\ge 2b\) or \(|\nabla(V_i-V_j)|\ge 2b\), for a
common \(b>0\). This follows from compactness and the top-tie
hypothesis. Choose \(r_0\) by uniform continuity so that in the first
case \(g_{ij}\ge b\) on the ball, and in the second case the
directional derivative in the direction of the contrast gradient
at its center stays at least \(b\) throughout the ball.
In the first case a nonempty sublevel requires \(t\ge b\), and the
ball's Gaussian area is at most a constant times
\(s^2\le (r_0/b)st\). In the second case \(g_{ij}\le t\) implies
\(|V_i-V_j|\le t\). Each line parallel to the chosen direction has
a slice of the latter sublevel in the ball of length at most
\(2t/b\) by one-dimensional monotonicity. Integrate in the
orthogonal direction over length \(2s\), bound the Gaussian density,
and sum over the three pairs. This proves (PCF1).

## Finite-amplitude comparison

Let \(P\) be a simplex-valued assignment with the prescribed masses
that agrees with \(A\) outside the cover. Put \(e=P-1_A\) and
\(\lambda(x)=1-P_{a(x)}(x)\) so that \(0\le\lambda\le1\).
Choose disjoint measurable pieces \(D_b\subset B_b\) covering the
union of the balls, and put \(m_b=\int_{D_b}\lambda\,d\gamma_2\).
The exact quadratic expansion and the cancellation of prices give
\[
 S_r(P)-S_r(A)=
 \sum_i\langle e_i,T_re_i\rangle
       -2\int\sum_{j\ne a(x)}P_j(x)[V_{a(x)}(x)-V_j(x)]\,d\gamma_2(x).
\]
The integrand of the loss is at least \(\lambda\eta\). Using
(PCF1) and integrating the sublevel bound in the threshold gives
\[
 \int_{D_b}\lambda\eta\,d\gamma_2
     \ge\int_0^\infty(m_b-Cs_b t)_+\,dt
     =\frac{m_b^2}{2Cs_b}.
\]
The standard inner product on label vectors and the pointwise bound
\(\|e(x)\|_2\le\sqrt2\,\lambda(x)\) give
\(\sum_i\langle e_i,T_re_i\rangle\le2\sum_{b,c}M_{bc}m_bm_c\).
With \(w_b=m_b/\sqrt{s_b}\), the exact comparison yields
\[
                   S_r(P)-S_r(A)\le -C^{-1}w^{\mathsf T}(I-G)w. \tag{PCF2}
\]
A nonnull change gives \(w\ne0\), which proves strictness under
the stated eigenvalue condition.

## Two junctions and thin interface patches

Let \(K_{r_0}=\{x:\operatorname{dist}(x,K)\le r_0\}\) and
\(\overline M=\sup_{K_{r_0}\times K_{r_0}}M_r<\infty\).
The largest eigenvalue of \(G\) is at most
\(2C\overline M\sum_b s_b\), by Cauchy–Schwarz. For a finite
set of junctions there is thus a strictly stable neighborhood
consisting of sufficiently small balls at those points. This applies
to two transverse junctions of opposite posterior orientation,
including a pair where the unadjusted scores agree exactly. A
transfer supported entirely in those neighborhoods decreases the
stability, regardless of how much mass it moves there.
[Posterior orientations and common score vectors](THREE_CELL_POSTERIOR_DEGREE_AND_TRIPLE_BISECTOR.md)

More generally, if
\(\mathcal H^1(K)<(2C\overline M)^{-1}\),
then there is an open neighborhood of \(K\) supporting no
improving mass-preserving change. Cover \(K\) by sets of diameter
less than \(r_0/2\) whose total diameter is arbitrarily close
from above to \(\mathcal H^1(K)\). Enlarge each to an open ball
centered at a point of \(K\) with radius at most its diameter plus
a summably small error. Compactness gives a finite subcover
obeying the needed bound on the sum of radii. Its union
contains a thin tube around \(K\). In particular the comparison
also applies to thin perturbation patches along a collection of
regular interface subarcs with sufficiently small total length
relative to the local gradient and kernel bounds.
