# Local rigidity of Gaussian score winners

For \(k\ge3\) work first in dimension \(n=k-1\). An adjusted-score
winner at correlation \(\rho\) is a balanced partition that chooses
a largest score \(T_\rho1_{A_i}-\beta_i\) almost everywhere for
some constant prices \(\beta_i\). Suppose at \(r\in(0,1)\) the
simplex interface form has a positive gap in the moment gauge
of (B12). Then any adjusted-score winner sufficiently close
to the simplex orbit at \(\rho\) sufficiently close to \(r\) is
the simplex. This also holds at \(r=0\) for small positive
\(\rho\), with the endpoint interface form. In particular it
holds locally at every positive correlation for three cells.

The argument uses the facet limit for arbitrary balanced
partitions, with an additional equality from the winner equation.
It assumes neither second variation nor an energy comparison
for the nearby competitor.
[Uniform weighted exterior bounds and facet compactness, B3–B14](POSITIVE_CORRELATION_CONTINUATION.md)
[Endpoint moment-gauge gap, L2](SMALL_POSITIVE_CORRELATION_INTERVAL.md#an-interface-estimate-in-the-moment-gauge)
[All-correlation three-cell gap, S37](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

## The winner equation identifies the limiting facet speed

Suppose otherwise, with balanced nonsimplex winners \(A^{(j)}\)
approaching the simplex orbit and \(\rho_j\to r\).
Align all references to \(T\) by the exact moment gauge of
(B12). Put \(L_j=L_{\rho_j}(A^{(j)};T)>0\),
\(s_j^2=L_j/\rho_j\), \(e_{j,i}=1_{A_i^{(j)}}-1_{T_i}\)
and \(\mu_{j,i}=e_{j,i}\gamma_n/s_j\). Then \(s_j\to0\).
For \(\rho_j\) bounded away from zero this follows from bounded
scores and the vanishing changed measure. Near zero use the
normalized kernel bound (B4) and Gaussian integrability of the
weight in the fixed active plane. The constant contribution
to the normalized simplex score differences is zero by balance.

The compactness proof of (B13)–(B14) uses the transfer estimates,
balance and the moment gauge up to its final contradiction.
It gives, along a subsequence, \(\mu_{j,i}\rightharpoonup
\lambda_i^f\) for a facet speed \(f\) in the balanced
moment-gauge space \(\mathcal G\), and
\[
 Q_{\rho_j}(e_j)/L_j\longrightarrow B_r(f).
\]
At \(r=0\) the same limit uses the continuous extension of
\((M_\rho-1)/\rho\) to \(K_0(x,y)=x\cdot y\) in (B4).
Here \(Q_\rho(e)=\sum_i\langle e_i,T_\rho e_i\rangle\);
the means of each \(e_i\) vanish. The same product-kernel
argument, with one argument in a compact set, gives locally
uniform convergence
\[
 w_{j,i}(x)=\frac{T_{\rho_j}e_{j,i}(x)}{\rho_j s_j}
            \longrightarrow
 w_i(x)=\int K_r(x,y)\,d\lambda_i^f(y),\qquad
 K_r=(M_r-1)/r \text{ for }r>0 .
\]
To obtain uniformity, truncate the \(y\) tails by (B3)–(B4)
and use uniform continuity and weak convergence on each compact
product. The same bound makes all \(w_{j,i}\) uniformly bounded
on the compact observation set.

Write \(\lambda_{j,i}=\beta_{j,i}/\rho_j\) in the zero-sum
gauge. The normalized price differences are bounded: every label
wins somewhere in one fixed compact ball, by positive masses
and convergence in label measure, and differences of normalized
scores are bounded there. Any limit prices leave \(T\) an
adjusted-score winner for its limiting normalized scores.
Indeed the scores converge uniformly on compacts and the labels
converge in measure. Continuity approaching a relative interior
of any reference facet forces the limiting price difference
there to vanish, since the reference pair scores tie on the
facet. Thus \(\lambda_j\to0\).

More precisely \(|\lambda_{j,i}-\lambda_{j,l}|=O(s_j)\).
On a fixed compact interior patch of the \(i,l\) facet, take
signed normal distance \(t>0\) into cell \(T_i\). The normalized
reference gap there is bounded above by \(C t\), and all
other reference labels have a fixed positive gap on a small
patch neighborhood. The perturbation of normalized pair scores
is bounded by \(C's_j\). If \(d=\lambda_{j,i}-\lambda_{j,l}\)
is positive and exceeds \(4C's_j\), then \(A_i^{(j)}\) loses
throughout \(0<t<c d\) on that patch for a sufficiently small
fixed \(c>0\). Every such relabeling costs at least \(c't\)
in normalized reference regret, by the positive reference
normal derivative and the gaps to other labels. The strip
therefore contributes at least \(c''d^2\) to \(L_j/\rho_j=s_j^2\).
For \(d<0\) use the other side. This proves the bound for
each pair, and after a subsequence \(\lambda_{j,i}/s_j\to b_i\).

On a compact interior pair patch use \(t=s_j\tau\).
Only the two labels can win in a fixed sufficiently small
tube around that patch. Dividing their adjusted contrast by
\(\rho_j s_j\) gives, uniformly for bounded \(\tau\), the limit
\[
 \alpha_r(z)\tau+w_i(z)-w_l(z)-b_i+b_l ,
 \qquad \alpha_r=a_{r,il}/r ,
\]
with \(\alpha_0=|m_i(T)-m_l(T)|\). The reference normal
derivative is uniformly positive on the patch. The uniform
\(O(s_j)\) score and price perturbations also confine all
relabeling in the tube to \(|\tau|\le C\). The limiting winner
indicator is thus given almost everywhere by an affine threshold
in \(\tau\). Gaussian integration in the normal coordinate
identifies the signed flux from \(l\) to \(i\) as
\[
 \alpha_r(z)f_{il}(z)=w_i(z)-w_l(z)-b_i+b_l.                 \tag{LW1}
\]
The plane weak limit gives this same flux: normalized transfers
away from the fixed tube vanish by the reference regret bound.
Exhausting facet interiors gives (LW1) surface-almost everywhere.
There is no junction mass in the limit by the plane \(L^2\)
bound in (B13).

Multiply (LW1) by \(f_{il}\) and integrate over facets.
The price terms cancel by balance. The remaining equality is
\(A_r(f)=B_r(f)\), so the gap on \(\mathcal G\) forces \(f=0\).
On the other hand \(A^{(j)}\) wins the linear assignment with
its own scores against the balanced reference \(T\). Therefore
\[
 0\le\sum_i\langle e_{j,i},T_{\rho_j}1_{A_i^{(j)}}\rangle
       =-L_j+Q_{\rho_j}(e_j),
\]
and \(B_r(f)\ge1\), a contradiction. The sequential proof also
gives uniformity on compact ranges of correlations where the
gap holds, and near the endpoint \(r=0\) by the endpoint gap.

## A lower bound on the balanced reassignment gain

Under the same gap hypothesis and in the same active dimension,
there is \(c>0\) such that every balanced partition \(A\) in a
sufficiently small neighborhood of the simplex orbit satisfies
\[
 R_\rho(A)=\max_{B\ {\rm balanced}}
     \sum_i\langle 1_{B_i}-1_{A_i},\,T_\rho1_{A_i}\rangle
       \ \ge\ c L_\rho(A;T_A).                              \tag{LW2}
\]
Here \(T_A\) is the closest-moment orientation. The bound is
uniform as \(\rho\to r\), including the endpoint \(r=0\) as
above. In particular the inequality applies to approximate score
winners without any interface assumptions on their cells.

Suppose a sequence violates the conclusion with
\(R_{\rho_j}(A^{(j)})/L_j\to0\). Use the orientation,
\(s_j^2=L_j/\rho_j\), and compactness of the previous section
for the \(A^{(j)}\). Let \(B^{(j)}\) be a balanced linear maximizer
for its scores. A minimizer of the finite-dimensional price dual
exists by the coercive lower bound in the gauge of minimum price
zero. The score differences are real analytic and nonconstant
by injectivity, so differentiating the objective in prices gives
the prescribed volumes. Thus the ordinary winner for the adjusted
scores realizes the maximum.

Write \(e^B=1_{B^{(j)}}-1_T\), \(e^A=1_{A^{(j)}}-1_T\) and
\(L_j^B=L_{\rho_j}(B^{(j)};T)\), with the same reference \(T\).
Comparison of the linear objective at \(B^{(j)}\) and \(T\)
gives \(L_j^B\le \sum_i\langle e^B_i,T_{\rho_j}e^A_i\rangle\).
The global weighted mass bound of (B5) and the normalized
kernel bound (B4) bound this last expression by
\(C\sqrt{L_j^B L_j}\). Hence \(L_j^B=O(L_j)\).
The saturated margin and (B5) give convergence of \(B^{(j)}\)
to \(T\) in label distance, and bound its weighted discrepancy
by \(O(s_j)\). The facet compactness for \(B^{(j)}\) now
applies with the same normalization \(s_j\). Its limiting
flux \(f^B\) is balanced. Its orientation need not satisfy the
moment gauge for \(B^{(j)}\).

The price-strip argument above applied to \(B^{(j)}\) and driven
by the score perturbation \(w_{j,i}\) of \(A^{(j)}\) makes its
normalized prices \(O(s_j)\), because \(L_j^B/\rho_j=O(s_j^2)\).
It gives the equation (LW1) for \(f^B\) with the field \(w_i\)
from the \(A^{(j)}\) flux \(f^A\). To see \(f^A=f^B\), note
that \(R_{\rho_j}(A^{(j)})\) is exactly the integral of the
nonnegative pointwise adjusted gap between the label of
\(B^{(j)}\) and the label of \(A^{(j)}\). Their masses agree,
so the price terms do not affect the integral. In the rescaled
normal coordinate \(t=s_j\tau\) on a compact interior facet
patch, outside the set where the absolute limiting two-label
contrast is smaller than \(\varepsilon\), any disagreement
costs at least \(c_\varepsilon \rho_j s_j\) for \(|\tau|\le M\).
This follows from the compact score convergence and the fixed
separation to the other labels. Thus the normalized mass of
such disagreements tends to zero. The limiting small-contrast
slab has normalized mass \(O(\varepsilon)\) since
\(\alpha_r\) has a positive lower bound there. The mass of
disagreements in the fixed tube with \(|\tau|>M\), after
normalization, is \(O(1/M)\): either partition must change
the reference label, and its normalized reference loss there
is at least a constant times \(|t|\). Use \(L_j^B=O(L_j)\)
and \(L_j=\rho_j s_j^2\). Sending \(\varepsilon\downarrow0\)
and \(M\to\infty\) proves equality of the two local normal
fluxes. Exhausting compact facet interiors proves equality
in the flux space.

The equation for \(f^B=f^A\) now implies \(A_r(f^A)=B_r(f^A)\),
so \(f^A=0\) by the gap and the gauge on \(A^{(j)}\). On the
other hand using \(T\) as a competitor in the reassignment
definition gives \(R_{\rho_j}(A^{(j)})\ge L_j-Q_{\rho_j}(e^A)\).
Hence \(B_r(f^A)\ge1\) in the limit, a contradiction.

For fixed positive \(\rho\) with the gap, balanced score
reassignment therefore converges geometrically in label
distance to a simplex if initialized sufficiently near its
orbit. Indeed the local comparison gives
\(c_0 L_\rho(A;T_A)\le S_\rho(T)-S_\rho(A)\le2L_\rho(A;T_A)\).
With \(B\) a best response the exact stability increment is
\(2R_\rho(A)+Q_\rho(1_B-1_A)\). The estimate
\(L_\rho(B;T_A)\le C L_\rho(A;T_A)\) and the saturated
margin bound both the distance from \(B\) to \(T_A\) and
the successive label displacement by
\(C'\sqrt{L_\rho(A;T_A)}\) at this fixed correlation.
Choose the initial loss so that \(C'\sqrt{\Delta/c_0}\)
is less than the local comparison radius whenever
\(0\le\Delta\le S_\rho(T)-S_\rho(A_0)\). Then \(B\) lies
in the neighborhood by its distance to \(T_A\). The
local comparison and the increment bound show that its
loss decreases by a fixed factor less than one. This
keeps each iterate within the comparison neighborhood.
The square roots of their losses are summable, and the
iterates converge in label measure. Their distance to the
simplex orbit tends to zero by local coercivity. These
bounds are uniform over compact positive correlation bands
where the simplex gap holds.

## Stable planar winners at sufficiently small correlation

There is \(\varepsilon>0\) such that a balanced planar
three-cell adjusted-score winner at \(0<\rho<\varepsilon\)
is the simplex if it satisfies the ordinary constrained second
variation on its compact regular interface speeds and the
regular-interface hypotheses of the global interface bounds.
The full active-plane case uses (IB2), which is available for
stationary stable partitions. In the one-dimensional case,
the unused planar coordinate gives the nonnegative first
extra-coordinate form (G1) directly by testing odd tensor
speeds, with compact cutoffs. The one-transition argument (G8)
then excludes one active dimension. It uses only the analytic
winner rule, positivity of the third mass, and that form.

The Gaussian multi-bubble theorem gives \(P_\gamma(A)\ge
P_\gamma(T)\), with \(P_\gamma(T)^2=F_3=W_1(T)\).
Consequently (IB2) gives \(S'_\rho(A)\ge F_3\). The Hermite
expansion gives uniformly \(S'_\rho(A)=W_1(A)+O(\rho)\), since
for \(0<\rho\le1/2\) the multipliers of degree \(j\ge2\)
are bounded by \(C\rho\). The sharp first-moment bound and
its quantitative alignment therefore force convergence
to the simplex orbit in Gaussian label distance along any
such sequence with \(\rho\downarrow0\). The endpoint local
winner theorem just proved excludes any sequence of
nonsimplex winners with these hypotheses.
[Interface and perimeter bound for stable stationary partitions, IB2](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#the-same-bounds-for-ordinary-stable-stationary-partitions)
[Quantitative first-moment alignment](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md)

This also classifies every measurable balanced three-cell
local maximum in Gaussian label distance in any finite
dimension \(d\ge2\) for the same small correlations.
Such a local maximum wins its own linear assignment. Indeed
if a balanced competitor gave strict positive linear gain,
select the same small fraction of each directed transfer
from the upper quantiles of its linear gain. This preserves
masses and has strictly positive linear gain at arbitrarily
small distance. The quadratic Mehler term is nonnegative.
The price minimization for the bounded scores attains its
minimum in the gauge \(\min_i\beta_i=0\), since its objective
is bounded below by \(\max_i\beta_i/3\). The analytic pair
contrasts are nonconstant by Mehler injectivity. Thus the
price derivatives give the required masses and equality in
the pointwise score bound gives the winner rule almost
everywhere.

The bulk dimension reduction proof (SD1)–(SD3) then applies
to this local maximum: its exact-mass repaired translations
converge in label distance to the original partition, so
local maximality suffices in (SD1). The partition is
cylindrical off its centroid span. Project to a plane
containing that span. Planar competitors lift isometrically
in label distance, so the projected partition is a planar
local maximum. The analytic adjusted scores give locally
finite perimeter and regular two-label pieces covering the
surface measure. For this consequence use the stratification
of the finitely many analytic score inequalities as in the
proof of (G20). The codimension-one strata have two sides
with fixed winners, and the lower-dimensional strata have
zero surface measure. The polynomial length estimate in the
interface bounds gives local finiteness of the perimeter.
Exact mass-preserving variations on the regular pieces
give the ordinary stability hypothesis. The planar
classification now applies.
[Bulk repair and translation identity](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md)
[Analytic interface stratification](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

For positive correlations bounded away from zero, local winner
rigidity rules out stationary branches approaching \(T\)
through a positive spectral gap. It places no upper bound
on the stability of a distant stationary branch.
