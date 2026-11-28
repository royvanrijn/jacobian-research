# Dimension-free compactness of Gaussian near-maximizers

Fix \(k\ge2\), \(0<r<1\), and positive Gaussian cell masses \(m_i\).
Every sequence of ordinary partitions of these masses whose noise
stabilities tend to the unrestricted maximum has a subsequence
converging in Gaussian label distance to cylindrical copies of one
maximizer in at most \(k-1\) coordinates. The ambient dimensions may
tend to infinity. For three cells the embeddings can be chosen in
the planes spanned by the cell centroids. The result also holds
along convergent sequences of correlations and positive masses, with
the limit optimizer taken at the limiting parameters.

The proof uses Chris Jones's graded Gaussian regularity theorem,
specifically the isometric ultrapower statement in Theorem 3.5.
Jones's Section 4 uses the same compactification and promotion of
higher grades for the two-cell Borell inequality. The additional
steps here are dimension reduction for countable products and
recovering cylindrical copies through the grade-one embedding.
[Chris Jones, *Regularity lemma for Gaussian polynomials*, arXiv:2610.16180v2](https://arxiv.org/html/2610.16180v2)

## Countable Gaussian products

Write \(V_r(m)\) for the unrestricted finite-dimensional optimum.
It is attained in \(k-1\) coordinates by ordinary dimension
reduction. The supremum over simplex-valued tuples of the same
masses on a countable Gaussian product is also \(V_r(m)\).
To see the upper bound, condition the tuple on the first \(N\)
coordinates and use \(L^2\) martingale convergence. The resulting
finite-dimensional fractional tuple of masses \(m\) has stability at
most \(V_r(m)\). An optimizing tuple on the countable product must
be ordinary. Otherwise splitting two positive fractional coordinates
in opposite directions on two subsets of equal measure gives two
feasible tuples averaging to the optimizer. Strict positive
definiteness of \(T_r\) makes one of their values larger.
[Heilman–Tarter, Theorem 1.9](https://doi.org/10.1017/fms.2021.56)
[Finite-dimensional relaxed maximum](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#optimization-and-volume-repair)

The following extension of the bulk translation argument reduces
any such countable-product maximizer \(f\) to its centroid span.
Take mass prices \(\beta_i\) as in finite dimensions, applying a
supporting supergradient to the value of the fixed linear
assignment at the prescribed interior mass vector. The tuple
chooses winners of \(s_i=T_rf_i-\beta_i\) almost everywhere. Give the
label graph a directed edge \(i\to j\) when
\[
                  \mathop{\rm ess\,inf}_{\{f_i=1\}}(s_i-s_j)=0 .
                                                               \tag{GC1}
\]
This graph is strongly connected. Indeed if a nonempty proper label
group \(I\) has no outgoing edge, finiteness of the labels gives
an \(\eta>0\) such that
\[
 b=\max_{i\in I}s_i-\max_{j\notin I}s_j
 \quad\text{takes values in }(-\infty,0]\cup[\eta,\infty)
 \quad\text{almost surely}.
\]
The label group is \(\{b\ge\eta\}\) up to null sets. Each score
has a continuous version along almost every line in any fixed
coordinate, since the product noise operator factors into smoothing
in that coordinate and the bounded noise operator in the others.
Fubini and continuity forbid a line from meeting both sides of
the omitted interval on sets of positive one-dimensional measure.
Thus the group indicator depends on none of the coordinates. For
example all its nonconstant Hermite coefficients vanish, by testing
in a coordinate of positive degree. Its Gaussian mass is then zero
or one, contradicting positivity of all masses. This proves (GC1)'s
connectivity claim.

Translate \(f\) by a small multiple \(t v\) of a Cameron–Martin
direction \(v\) orthogonal to all centroids. Gaussian translation in
the coordinate \(v\cdot X\) gives mass errors \(O(t^2)\), since their
first derivatives are \(-\mathbb E[f_i(v\cdot X)]=0\). For every
edge (GC1) and every \(\epsilon>0\) there is a positive-measure
source patch where the score loss is at most \(\epsilon\). The
translated tuple converges to \(f\) in \(L^2\), by boundedness and
one-coordinate Gaussian translation continuity. The translated
scores converge in \(L^2\) as well. Each source patch therefore
still has a positive reservoir of the right source label with
translated score loss at most \(2\epsilon\). Route corrections of
total mass \(O(t^2)\) along fixed directed paths, using disjoint
subsets of these reservoirs. The noise interaction of the
correction with itself is nonnegative. Consequently the priced
second derivative along the uncorrected translation is nonpositive,
by the same estimate as (SD1) with \(\epsilon\downarrow0\).

The calculation (SD2)–(SD3) applies on the countable product.
One may choose \(v\cdot X\) as the first standard Gaussian
coordinate, after normalizing \(v\) and completing it to an
orthonormal basis of the first chaos. Along that coordinate the
maximum of the adjusted scores is locally Lipschitz, and the
active scores have its derivative almost everywhere. This gives
the distributional stationarity identity along almost every
line, then (SD2) by Gaussian integration in that coordinate.
The same smoothing estimates bound all needed derivatives
uniformly, including over the remaining coordinates. Hermite
lowering and its exponential summability give the priced second
derivative divided by two as
\[
        (1-r)\sum_{i,j\ge0}r^j\|P_jD_vf_i\|_2^2 .           \tag{GC2}
\]
Each term must vanish. Use a countable orthonormal basis of the
centroid complement to conclude that no Hermite coefficient
uses those directions. Thus \(f\) is a cylindrical copy of a
maximizer in its centroid span, which has dimension at most
\(k-1\) because the centroids sum to zero.
[Bulk translation identity](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#the-bulk-translation-identity)

## Graded limits and original linear coordinates

Let \(f^{(n)}\) be a sequence of ordinary near-maximizers. Take a
nonprincipal ultraproduct of their finite-dimensional Gaussian chaos
spaces. Jones's proof of Theorem 3.5 applies with the dimension of
each tensor space equal to that of its corresponding factor. Apply
the theorem simultaneously to the degree components
\([P_j f^{(n)}_i]\) of all labels. It identifies each component
isometrically with a graded Hermite component \(p_{ij}\).
The sum \(g_i=\sum_{j\ge0}p_{ij}\) exists in \(L^2\), since
for every finite truncation the squared norms are the ultralimits
of the corresponding ordinary squared norms. The means of \(g_i\)
are \(m_i\).

For \(0<s<1\), the vector of graded-smoothed \(g_i\) is the joint
distributional ultralimit of the vector of ordinary-smoothed
\(f^{(n)}_i\). Indeed the degree truncation errors are uniformly
bounded in \(L^2\) by a constant times \(s^{D+1}\) beyond degree
\(D\), and the finite truncations converge jointly by Theorem 3.5.
The latter ordinary vectors take values in the closed probability
simplex. Therefore so do the graded-smoothed limits, and letting
\(s\uparrow1\) gives the simplex property for \(g\). Likewise the
graded stability of \(g\) is the ultralimit of the stabilities
of the sequence, by summing the squared graded norms with
multipliers \(r^j\). This is the labelwise extension of Jones's
Lemma 4.2.

Treating every graded variable as an ordinary grade-one variable
does not change the marginal masses or the probability simplex
constraint. It strictly increases stability if any graded Hermite
coefficient uses a higher-grade variable, since that coefficient
then receives a strictly greater multiplier. Countably many
ordinary coordinates cannot exceed \(V_r(m)\) by the argument
above, and the original graded stability of \(g\) is already
\(V_r(m)\). Thus \(g\) uses only grade-one variables. It is an
ordinary optimizer, cylindrical over its centroid span by (GC2).

Write its active dimension as \(a\) and its profile as \(h\).
An orthonormal basis of its active span consists of \(a\) elements
of the grade-one ultrapower in Jones's representation. Choose
linear representatives \(U_{n,1},\ldots,U_{n,a}\) in the original
Gaussian spaces. They are asymptotically orthonormal along the
ultrafilter by the grade-one isometry. On an ultrafilter-large set,
Gram–Schmidt gives
exactly orthonormal representatives with the same ultrapower
classes.

For a polynomial \(q\) in \(a\) coordinates, the graded
representative of each homogeneous Hermite degree of \(q(U_n)\)
is that degree of \(q\) in the active coordinates of \(h\). This
follows from the Wick product of the orthonormal grade-one
variables in Theorem 3.5. The isometry gives
\[
       \lim_{\mathcal U}\mathbb E[f^{(n)}_i q(U_n)]
                     =\mathbb E[h_i q].                    \tag{GC3}
\]
Approximate the indicator \(h_i\) in \(L^2(\gamma_a)\) by such
polynomials. The approximation error has the same norm on
each \(U_n\) because \(U_n\) is a standard Gaussian vector.
It follows that \(\lim_{\mathcal U}\mathbb E[f^{(n)}_i h_i(U_n)]=m_i\).
After summing in \(i\), the Gaussian label discrepancy tends to
zero along the ultrafilter. Take a subsequence with ordinary
convergence. No weak compactness of the full noise operator
in the countable product is needed.

For \(k=3\) the centroids of \(h\) span its active plane. Their
components in the recovered plane converge by (GC3) with
linear \(q\), and their full squared norms converge by the
degree-one isometry. The centroids of \(f^{(n)}\) consequently
have vanishing components perpendicular to the recovered plane.
The limiting rank is two, so their actual spanning planes have
orthogonal projections converging to those of the recovered
planes in operator norm. Replacing the recovered orthonormal
coordinates by ones in those centroid planes preserves label
convergence. This follows for polynomials by Gaussian moment
continuity, then for indicators by \(L^2\) approximation.
[Full active moment rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

For varying \(r_n\to r\in(0,1)\) and \(m^{(n)}\to m\) with
positive limit masses, the same proof applies. The ordinary
operator norms satisfy
\(\sup_{j\ge0}|r_n^j-r^j|\to0\), and the optimum is continuous
in the masses by direct transfers. Thus the stability ultralimit
is \(V_r(m)\), and the degree-zero components give the limiting
masses. In particular compactness gives a uniform qualitative
near-optimality implication on any compact set of such
parameters, measured against cylindrical copies of limiting
optimizers along convergent parameter sequences.

## Uniform full moment rank for three cells

Fix a compact set of positive correlations and positive three-cell
mass triples. There exist \(\varepsilon,\delta>0\) such that in
every ambient dimension admitting the unrestricted value, any
partition with stability at least \(V_r(m)-\varepsilon\) has second
centroid singular value at least \(\delta\). Otherwise a sequence
of near-maximizers with vanishing second singular value would
converge in label distance, after identification of its recovered
plane, to a planar maximizer. The full centroid Gram matrices
converge by the degree-one isometry and (GC3). Their limit has
rank two by (G8), which contradicts the vanishing singular value.
In particular every sufficiently high-dimensional conical
maximizer has full centroid rank uniformly on the compact set,
since (AX3) makes its stability deficit at most \(O(d^{-1})\).
[Uniform conical approximation](GAUSSIAN_CONE_AXIAL_INSTABILITY.md#a-uniform-expansion-for-radial-spectator-blocks)

## A dimension-free gap away from a unique maximizing orbit

Let \(K\) be a compact set of correlations in \((0,1)\) at which
the centered simplex is the unique balanced three-cell global
maximizer up to rotation and relabeling. For each \(\eta>0\) there
is \(\delta(K,\eta)>0\) such that every balanced partition at
Gaussian label distance at least \(\eta\) from the simplex orbit
has deficit at least \(\delta(K,\eta)\) at every \(r\in K\),
uniformly in the ambient dimension. Interpret both the distance
and the simplex by cylindrical extension if the original
dimension is one.

Otherwise a sequence of these partitions has deficits tending
to zero. Extract a subsequence of correlations converging in
\(K\) and apply the compactness theorem. The limit maximizer is
a simplex by uniqueness at the limiting correlation. Convergence
to its cylindrical copies contradicts the lower bound \(\eta\).
At a candidate transition correlation a sequence staying away
from the simplex orbit can therefore approach optimal value only
by converging, after orthogonal identification and passage to a
subsequence, to a genuine nonsimplex planar maximizer.
