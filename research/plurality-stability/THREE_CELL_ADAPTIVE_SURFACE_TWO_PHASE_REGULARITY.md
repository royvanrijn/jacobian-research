# Regularity of two-phase regions of adaptive surface minimizers

Fix \(0<r<1\). Let \(A\) be any balanced planar minimizer
of the original adaptive surface functional \(\mathcal F_r\),
using two direction coordinates. Suppose there is an open
neighborhood of \(p\) in which only labels \(i,j\) occur up
to null sets, and the tagged interface measure satisfies
\[
 \liminf_{t\downarrow0}\frac{\nu_{ij}(B(p,t))}{t}>0 .
\]
Then the interface in a neighborhood of \(p\) is real
analytic, for the representative given by the support of
the perimeter measure. In particular this applies at a
reduced-boundary point with a two-phase neighborhood.
The proof first gives a local almost-minimizing inequality
for ordinary perimeter. The stationary metric perimeter
equation then makes the resulting graph analytic.
[The adaptive functional and dual test representation](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md#a-dual-formula-with-metric-surface-weights)
[Attaining tests, positive interface weights and the pressure equation](THREE_CELL_ADAPTIVE_SURFACE_STATIONARITY.md)

## Nearby dual tests and an exterior mass correction

Write \(s=\sqrt r\) and \(Z_A(u)=T_s m_A^u\) for the tensor in
the adaptive functional. The set of \(Z_A(u)\) for admissible
directions is norm compact. Each \(T_s\nu_e\) is in \(L^2\),
as follows by taking a direction supported on the pair \(e\).
Weak-* convergence of \(u\) gives pointwise convergence after
smoothing. The sum of the smoothed unsigned pair measures
dominates these tensors in \(L^2\), giving norm convergence.
The set \(\mathcal M(A)\) of normalized tensors
\(Z_A(u)/\mathcal F_r(A)\) at attaining directions is
therefore compact. The minimum is finite and positive by
the simplex competitor and the perimeter lower bound.
[Existence and the perimeter bound](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md#existence-and-a-strict-perimeter-bound)

Choose a small closed disk \(D\) about \(p\) contained in the
two-phase neighborhood, with each label still having positive
mass in its open complement. There are two smooth compactly
supported vector fields in this complement whose Gaussian
mass first variations span the two-dimensional mass plane.
Indeed a linear combination of indicators that annihilates
these variations has zero Gaussian distributional gradient
on the open complement of \(D\). It is constant there since
this region is connected. Each label occurs there with
positive mass, so the three coefficients are equal.
Two local flow parameters and the implicit function theorem
give a mass correction outside \(D\) for every sufficiently
small error \(\delta\) in the masses. The parameter norm is
\(O(|\delta|)\). A competitor \(B\) which agrees with \(A\)
outside \(D\) has the same mass response under those exterior
flows, since their supports are disjoint from \(D\). The
corrected competitor will be denoted by \(B^*\).

Let \(U\subset D\) be a sufficiently small ball or closed
disk where the competitor is changed. Denote by \(P(B;U)\)
the sum of the unweighted tagged perimeters in \(U\), and
similarly for \(A\). Assume the two partitions are identical
on the open complement of the closure of the change
region, with all changed tagged measures counted in \(U\).
For any admissible directions \(u\) of \(B^*\), copy their
values to \(A\) where its interfaces are unchanged, pull
them back under the mass-fixing flow on the exterior
interfaces, and set them equal to zero in the change region.
This gives directions \(v\) on \(A\) with
\[
 \|Z_{B^*}(u)-Z_A(v)\|_2
       \le C\big(P(A;U)+P(B;U)+|\delta|\big).               \tag{AR1}
\]
For the changed pieces this follows from the uniform
\(L^2\) bound of the point-source Mehler kernel on compact
source sets and the triangle inequality. On the compact
flow support, that kernel and the Gaussian tangential
Jacobian vary at most linearly in the flow parameters
in \(L^2\) and uniformly, respectively. This bounds the
exterior contribution by \(C|\delta|\), since the
unchanged reference has finite Gaussian perimeter there.
All these compact perturbations have finite adaptive energy
by the same compact-source estimate and the unchanged tails.

If the right side of (AR1) tends to zero along balanced
corrected competitors, every attaining test of \(B^*\)
is uniformly close in \(L^2\) to \(\mathcal M(A)\).
To see this, choose an attaining \(u\) with norm
\(\|Z_{B^*}(u)\|_2=\mathcal F_r(B^*)\). Minimality of \(A\)
and (AR1) put this norm between \(\mathcal F_r(A)\)
and \(\mathcal F_r(A)+o(1)\). Norm compactness of the
reference tensors gives, on any sequence, a subsequence
of \(Z_A(v)\) converging to a tensor of norm
\(\mathcal F_r(A)\). Dividing by the norms proves the
assertion by contradiction.

For a unit-ball test \(G\) put
\(w^G_{ij}=|T_sG_i-T_sG_j|\) and
\(f^G_{ij}=\phi_2 w^G_{ij}\), where \(\phi_2\) is the
planar Gaussian density. Mehler smoothing maps the Hilbert
unit ball to uniformly bounded locally Lipschitz fields on
real compact subsets. Evaluation there is continuous in
\(G\) in the \(L^2\) norm. The reciprocal weight bound for
attaining tests implies \(w^G_{ij}(p)>0\) for every
\(G\in\mathcal M(A)\). Compactness makes the infimum of
these values strictly positive. Shrink to a neighborhood
\(D_0\) of \(p\) with closure inside \(D\). There are
uniform positive lower bounds of \(f^G_{ij}\) in \(D_0\)
for tests \(G\) close enough to \(\mathcal M(A)\). All
these weights have a common local Lipschitz bound.

## Ordinary perimeter almost-minimality

For a competitor \(B\) changed in \(D_0\), use an attaining
test \(G\) of \(B^*\). The extremal inequalities give
\(W_G(A)\le\mathcal F_r(A)\le\mathcal F_r(B^*)=W_G(B^*)\),
where \(W_G\) is the dual weighted perimeter. The exterior
mass correction changes \(W_G\) by at most \(C|\delta|\),
uniformly for all unit-ball tests, by the local Lipschitz
bounds and Gaussian tangential Jacobian on its support.
For a two-phase competitor in \(U\) this proves
\[
 \int_{U\cap\partial^*A_i} f^G_{ij}\,d\mathcal H^1
 \ \le\
 \int_{U\cap\partial^*B_i} f^G_{ij}\,d\mathcal H^1
       +C|\delta|.                                          \tag{AR2}
\]
Both tagged interfaces in \(U\) are the respective reduced
boundaries of the \(i\)-cell, up to surface-null sets.
All constants may depend on \(A,r,p\).

First let \(U=\overline{B(x,t)}\subset D_0\) and fill
the ball with label \(i\), leaving the other cells unchanged
outside. The new local perimeter is at most \(2\pi t\);
for almost every \(t\) this follows by the usual slicing
and gluing formula for finite-perimeter sets on the circle.
The mass error is \(O(t^2)\). The old local perimeter
\(P(A;\overline{B(x,t)})\) tends to zero uniformly for \(x\)
in a compact neighborhood of \(p\). Indeed the finite
local perimeter measure has no atoms. If uniform smallness
failed, a sequence of centers would converge to a point
carrying an atom. Thus (AR1) guarantees the positive
lower bound for the weight of the attaining competitor
test for all small \(t\), uniformly in those centers.
In (AR2), the lower bound on the old interface and the
upper bound on the new one give
\(P(A;B(x,t))\le C t\) for almost every small \(t\).
Monotonicity extends the estimate to every small \(t\).

For any two-phase competitor with
\(A_i\mathbin{\triangle} B_i\Subset B(x,t)\subset D_0\),
the case \(P(B;B(x,t))\ge P(A;B(x,t))\) already satisfies
the desired comparison. Otherwise the two local perimeters
are \(O(t)\) by the preceding bound. Again (AR1) gives
an attaining test of \(B^*\) with a common positive lower
bound on its weight near \(p\). In (AR2) replace the weight
on both sides by \(f^G_{ij}(x)\). The resulting error is at
most the product of \(Ct\) and the sum of the two local
perimeters. The mass error is \(O(t^2)\). Dividing by the
uniform lower bound proves
\[
 P(A;B(x,t))\le P(B;B(x,t))+C t^2                         \tag{AR3}
\]
for all sufficiently small balls centered in a smaller
neighborhood of \(p\). Every ordinary competitor for the
set \(A_i\) in this neighborhood is admissible in (AR3)
by assigning its local complement to label \(j\). This
is the planar almost-minimality estimate for Euclidean
perimeter. Apply Simmons's Theorem 1.1 with the identity
matrix and, for example, Hölder parameter \(1/2\).
The error \(Ct^2\) is bounded by \(Ct^{3/2}\) on unit
scales. There is no singular set in dimension two, so
\(\partial^* A_i\) there agrees with the support of its
perimeter and is locally a \(C^{1,1/8}\) graph.
[Simmons, Regularity of Almost-Minimizers of Hölder-Coefficient Surface Energies, Theorem 1.1](https://arxiv.org/pdf/2008.06150)

## The analytic graph equation

Use the stationary mixture of attaining tests and its
pressure vector. The averaged Euclidean weight is
\(f=\phi_2\int w^G_{ij}\,d\pi(G)\). It is positive and
real analytic near \(p\). In detail, on a fixed complex
neighborhood of \(p\), the vectors \(T_sG_i-T_sG_j\)
have holomorphic extensions with locally uniform bounds
and derivative bounds for all unit-ball \(G\). This
follows by Cauchy–Schwarz against the complexified
Mehler kernel, whose source Gaussian makes its \(L^2\)
norm locally bounded. On the real neighborhood the
norms for attaining tests are bounded below. On a smaller
common complex neighborhood, the sum of the component
squares stays in the right half-plane by the uniform
derivative bound. Its positive square root extends the
real norm. Averaging the holomorphic extensions against
\(\pi\) gives an analytic extension of \(f\).

In coordinates where the regular interface is the graph
\(y=g(x)\), the first variation of stationary weighted
perimeter against graph bumps gives in distributions
\[
 \frac{d}{dx}
       \frac{f(x,g(x))g'(x)}{\sqrt{1+g'(x)^2}}
   = f_y(x,g(x))\sqrt{1+g'(x)^2}
       -b\,\phi_2(x,g(x)),                                  \tag{AR4}
\]
where \(b\) is the oriented difference of the two mass
pressures. The graph is \(C^1\), so the right side is
continuous. The momentum on the left is therefore
\(C^1\), and it can be solved locally for \(g'\) since
its derivative with respect to slope is
\(f(x,g)/(1+g'^2)^{3/2}>0\). This makes the ordinary
differential equation for \((g,g')\) classical, with
analytic coefficients. Its unique local solution is
analytic. The argument uses a two-phase open neighborhood
for the comparisons and leaves points with arbitrarily
close third-phase incursions open.
