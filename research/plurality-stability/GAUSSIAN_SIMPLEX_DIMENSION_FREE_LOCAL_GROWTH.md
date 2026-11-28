# Dimension-free measurable local growth at the Gaussian simplex

For a fixed number of cells, suppose the regular simplex interface
has a positive gap away from rotations, uniformly in ambient
dimension, on a compact correlation band \(J\subset(0,1)\).
Then the measurable local comparison (B11) holds uniformly in
dimension and for all correlations in \(J\). In particular, for
three cells the hypothesis holds on every such band by (S37).
There are positive \(c,c',\eta\), independent of the ambient
dimension, such that every balanced three-cell partition \(A\)
with Gaussian label distance less than \(\eta\) from the simplex
orbit satisfies
\[
 S_r(T)-S_r(A)\ge c L_r(A;T_A)\ge c' d_3(A)^2,
                  \qquad r\in J .                         \tag{DL1}
\]
Here \(T_A\) is the nearest first-moment alignment, \(L_r\) is
the linear score regret, and \(d_3\) is label distance to the
simplex orbit.
[The active-space proof and its quantitative truncation](POSITIVE_CORRELATION_CONTINUATION.md#an-effective-local-radius-from-truncated-normal-flux)
[The ambient-dimension interface gap (S37)](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

The proof combines the positive-correlation score margins with
the conditional mass estimates and spectator smoothing used for
opposite simplex pairs. The saturated positive margin permits
every compact interior positive correlation band.
[Dimension-free transfer estimates for opposite simplex pairs](NEGATIVE_SIMPLEX_INTERFACE_LIFT.md#nearby-measurable-pairs-in-every-ambient-dimension)

## Normal flux with a fixed active cutoff

Align the reference simplex and write a point as \((x,z)\) in
the active \(p=k-1\) dimensions and the remaining spectator
dimensions. Write \(M_r^p\) for the active Mehler kernel and
\(T_r^\perp\) for the spectator operator. Let \(h_{ij}\) be the
transfer indicator from reference cell \(i\) to destination \(j\)
and put \(s^2=L_r/r\). Score margins depend only on \(x\), so
the saturated lower bound on each margin controls
\(r\min(t,1)\), where \(t=x\cdot n_{ij}\ge0\) is the active
normal coordinate. The normalization is used when \(s>0\).
Zero regret makes all transfers null by the margin bound.

For each spectator value \(z\), regard the signed transfer
tuple divided by \(s\) as active measures \(m_z\). For an active
measure tuple \(m_z\) let \(F_m(z)\) be its total weighted
variation summed over labels, using \(w_\alpha(x)\) from (B2).
Choose the upper correlation bound \(b<1\) for \(J\) in that
weight. Applied in the active variables for each \(z\), the
Gaussian one-sided mass bound (B5) and its global form give
\[
 \|F_m\|_{L^2(\gamma_\perp)}\le C,
 \qquad
 \|F_{m,\mathrm{out}}\|_2\le C\sqrt{\beta_{p,\alpha}(R)} ,   \tag{DL2}
\]
where the exterior cutoff is \(|x|>R\) and the same bounds
hold for the sum of the individual transfer variations.
Integration in \(z\) uses the full regret divided by \(rs^2= L_r\).

Apply the active truncation from (B15) to every slice, keeping
only \(|x|\le R\), \(0\le t\le h<1\) and a normal projection
onto the corresponding active facet. The normalized total
mass with \(|x|\le R\) and \(t>h\) is \(O(s/h)\) by the score
margin. Its weighted mass function in \(z\) has \(L^1\) norm
\(O_R(s/h)\) and bounded \(L^2\) norm by (DL2). In each slice
the total normal mass density on the whole pair plane has
squared \(L^2\) norm controlled by the slice regret divided by
\(rs^2\), by the one-sided normal density inequality. Thus its
squared norm integrated in \(z\) is bounded. The discard for
\(t\le h\) with projection outside the facet lies in the fixed
active pair-plane strips of Gaussian measure \(O_k(h)\) in
(B15). Cauchy–Schwarz in that plane bounds the discarded
weighted mass function in spectator \(L^2\) by \(O_R(\sqrt h)\).

For \(a\le r\le b\), the discarded weighted mass functions
therefore satisfy
\[
 \|T_{\sqrt r}^\perp F_{\mathrm{discard}}\|_2
  \le C\sqrt{\beta_{p,\alpha}(R)}
       +C_R\bigl(\sqrt h+(s/h)^{(1-b)/(1+b)}\bigr).        \tag{DL3}
\]
The first and third cutoff types use contraction on \(L^2\).
For the large-normal-distance part, Gaussian hypercontractivity
bounds the smoothed norm by \(\|F\|_{1+r}\). Interpolation
between \(L^1\) and \(L^2\) gives the stated exponent for
\(s/h\le1\). The bounds are independent of spectator dimension.
[Active exterior and normal density bounds (B3)–(B5)](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)
[Nelson hypercontractivity as stated in Gross, Theorem 1.1](https://doi.org/10.1007/BF02392573)

## Interaction and moment errors

Use the uncentered product Mehler kernel and divide the interaction
by \(r\). On \(J\) the active kernel divided by \(r\) is bounded
by a constant times \(w_\alpha(x)w_\alpha(y)\). The absolute
bilinear interaction of active measure families \(m_z,n_{z'}\)
is therefore at most
\[
 C\langle F_m,T_r^\perp F_n\rangle
 \le C\|T_{\sqrt r}^\perp F_m\|_2\|T_{\sqrt r}^\perp F_n\|_2 .
                                                               \tag{DL4}
\]
Interactions involving a discard cost at most the right side
of (DL3) times a constant. Move each retained mass to its active
normal projection with its spectator coordinate fixed. Both
active points lie in a fixed ball and move by at most \(h\).
The change in the active kernel is bounded by \(C_Rh\) by
compact derivative bounds. The interaction error is then
\(O_R(h)\) by the same spectator estimate.

Let \(f_{ij}\) be the signed projected flux on the product
facet. The normal regret comparison used in (B16) gives
\[
                  A_r(f)\le 2/(1-C_R h).                 \tag{DL5}
\]
Its \(L^2\) norm is uniformly bounded by the lower bound on
the normal score coefficient. The balance and active first
moment defects are controlled by (DL3) plus \(O_R(h)\).
Indeed \(w_\alpha\) dominates the active mass and moment weights,
its integral over a discard is controlled by the smoothed norm
by Cauchy–Schwarz and preservation of the mean, and projection
moves active coordinates by at most \(h\).

The retained projection preserves spectator first moments.
For each label the discarded spectator moment vector consists
of the degree-one Hermite coefficients of its discarded signed
mass function. Its Euclidean norm is at most \(r^{-1/2}\) times
the smoothed \(L^2\) norm of that function, since \(T_{\sqrt r}^\perp\)
multiplies the first level by \(\sqrt r\). Positivity and the
weighted variation bound control that smoothed norm by (DL3).
Summing over the fixed label count controls the full moment
defect in Frobenius norm.

## Uniform repair and comparison

The mass correction in (B15) uses a fixed spanning tree of
facets and constant speeds independent of the spectators.
Its cost and induced moment change are bounded uniformly.
For rotations within the active block the moment correction
belongs to a fixed finite-dimensional space. For a rotation
mixing the active and spectator blocks, the facet speed has
the form \((Vz)\cdot n_{ij}\) for a linear map from spectator
to active coordinates. Its squared norm is
\(\sum_{i<j}\nu_{ij}(F_{ij})|V^T n_{ij}|^2\), equivalent to
\(\|V\|_F^2\) with uniform constants because the active normals
span. The moment tangent has norm comparable to the same
Frobenius norm because the regular active moment matrix has
full rank and fixed singular values. The two rotation blocks
are orthogonal in the speed space and the moment matrix space.
The moment map on speeds is itself uniformly bounded, by
Bessel in the spectator first chaos and the fixed active
moment weights. Thus the gauge coercivity in (B12) and both
repair maps are uniform in the spectator dimension.

Put \(E=C\sqrt{\beta_{p,\alpha}(R)}+
C_R(\sqrt h+(s/h)^{(1-b)/(1+b)}+h)\). The full unprojected
transfers satisfy the volume and moment constraints by balance
and nearest-moment alignment. The repairs give a balanced
speed \(g\) with moments perpendicular to the rotation tangent
and with \(\|g-f\|_2\le CE\). The two interface terms are
uniformly bounded on the speed space because their active
operators tensor with identity or \(T_r^\perp\). The centered
and uncentered quadratic terms coincide for \(g\) and for the
original transfers by balance. Before repair use the uncentered
interaction comparison (DL4). The uniform gauge-space estimate
\(B_r(g)\le(1-\theta)A_r(g)\) yields
\[
 \frac{Q_r(e)}{L_r}\le
        \frac{2(1-\theta)}{1-C_Rh}+CE .                   \tag{DL6}
\]
Choose \(R\) large, then \(h\) small, then \(s\) small so the
right side is at most \(2-\theta\). The deficit identity
\(S_r(T)-S_r(A)=2L_r-Q_r\) gives the first inequality in (DL1).

The active slice mass estimate (DL2) bounds the squared aligned
label distance by \(Cs^2\) by Cauchy–Schwarz in spectator
variables. To reach small \(s\) from small distance to the
orbit, use the aligned linear regret estimate after (B16).
Bessel for the first Gaussian chaos bounds the moment distance
from any reference simplex by a constant times the square
root of their label distance. Nearest-moment alignment decreases
this moment distance. The scalar product with the fixed-norm
simplex moment matrix bounds the aligned linear regret, which
bounds \(s^2\) from above by the upper reflection bound. All of
these bounds are independent of dimension. This proves (DL1).

## Exact conical delay at a possible first planar tie

Suppose \(r_0\in(0,1)\) is the endpoint of the initial interval
of balanced three-cell simplex optimality. At \(r_0\) the simplex
\(T\) is a global maximizer by continuity. Let \(\operatorname{Opt}(r_0)\)
be the balanced planar maximizers and define
\[
 \tau_* = \inf\left\{
 \frac{C_{r_0}(A)}{S'_{r_0}(A)-S'_{r_0}(T)}:
 A\in\operatorname{Opt}(r_0),\ S'_{r_0}(A)>S'_{r_0}(T)\right\},  \tag{DL7}
\]
where the infimum of the empty set is \(+\infty\), and the
derivatives are taken with the partitions fixed. This number is
positive. Indeed the nonsimplex maximizers are separated from
the simplex orbit by local growth. They form a compact set in
the planar label metric, on which the dilation loss has positive
minimum, and the correlation derivatives are bounded.
[The dilation coefficient and its continuity in the active law expansion](GAUSSIAN_CONE_AXIAL_INSTABILITY.md#exact-leading-loss-for-the-entire-radial-profile-class)

For every \(0\le t<\tau_*\), all sufficiently large \(n\) have
the simplex as the unique conical maximizer up to its orbit
throughout \([r_0,r_0+t/n]\). Otherwise choose nonsimplex conical
competitors \(F_n\) with at least simplex stability at some
\(r_n=r_0+t_n/n\), \(0\le t_n\le t\). They are unrestricted
near-maximizers at \(r_0\). By dimension-free compactness they
converge after a subsequence and rotation to a cylinder over a
planar maximizer \(A\), and \(t_n\to u\in[0,t]\). If \(A\) is a
simplex, (DL1) at \(r_n\) contradicts that the competitors are
nonsimplex. Otherwise the unconditional lower bound (DC2) and
the same Taylor estimate as in (DC7) give
\[
 \limsup_n n\bigl(S_{r_n}(F_n)-S_{r_n}(T)\bigr)
   \le u\bigl(S'_{r_0}(A)-S'_{r_0}(T)\bigr)-C_{r_0}(A)<0.
\]
The strict inequality follows from \(u\le t<\tau_*\) and
positive dilation loss for every nonsimplex \(A\), including
when its derivative is at most the simplex derivative.

If \(t>\tau_*\) is finite, a fixed planar maximizer in (DL7)
with ratio less than \(t\) has balanced radial lifts beating
the simplex at \(r_0+t/n\) for all sufficiently large \(n\),
by the sharp radial coefficient and the same Taylor estimate.
Consequently the endpoint \(r_{2+n}^{\mathrm{cone}}\) of the
initial conical simplex-optimal interval satisfies
\[
         \lim_{n\to\infty} n(r_{2+n}^{\mathrm{cone}}-r_0)=\tau_* .
                                                               \tag{DL8}
\]
If \(\tau_*=+\infty\), the limit is understood as divergence
to infinity. This determines the scale of the conical delay at
a transverse planar transition and remains compatible with such
a transition at \(r_0\).
[Dimension-free compactness of near-maximizers](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)
[Full conical comparison and correlation expansion (DC2), (DC7)](GAUSSIAN_CONE_DECONING_ASYMPTOTICS.md)
