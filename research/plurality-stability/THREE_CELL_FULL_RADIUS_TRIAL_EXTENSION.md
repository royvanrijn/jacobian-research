# An unrestricted three-cell comparison through correlation 2611/4000

Every ordinary balanced Gaussian three-cell partition in dimension at
least two has noise stability at most that of the simplex for
$0<\rho\le2611/4000$. On the added interval $[2589/4000,2611/4000]$,
equality requires a simplex up to rotation, relabeling and null sets.
Existence and dimension reduction apply as in the preceding result.
[Projected-trial comparison](THREE_CELL_PROJECTED_TRIAL_EXTENSION.md)

The [directed checker](check_three_cell_full_radius_trial_extension.py)
reuses the spectral block and the degree-one variance credit. It uses the
full radial trial through a bounded core to strengthen the clipping
cost. It also retains the forced trial component of the projected
remainder and uses a one-sided central-profile reserve. The result has
internal verification only.

## A full-radius table for the clipping cost

Use the notation of (PT1), with $t=d r^2$, $d=1-\rho^2$, and
$x_0=\sqrt{d_u}r$ at slab upper endpoint $u$.
The physical trial is $\phi_\rho(t)=\sqrt{d_u/d}\,\widehat\phi(x_0^2)$.
Write $F(x)=x\widehat\phi(x^2)$, $\beta_u=1/d_u$ and let $Q(x)$ be the
positive mixed row polynomial at the anchor. Substituting the anchor
variables in (PT1) gives its two integrands

\[
 \frac{Q(x)e^{-\beta_u x^2}}{2d_u}
    (|F(x)|-v x)_+^2,
 \qquad
 \frac{Q(x)e^{-\beta_u x^2}}{2d_u}
    |F(x)|(|F(x)|-v x)_+,
 \qquad v=\frac{\sqrt{2/3}}{A}\sqrt{d/d_u}.                  \tag{FR1}
\]

The checker bounds $v$ above by using $d_l=1-l^2$. For each interval
in a partition of $[0,4]$, it bounds $|F|$ below by $f_i$, $x$ above
by $x_i$ and the integral of the positive weight in (FR1) below by
$w_i$. Intervals with $f_i=0$ can be omitted. On the subset certified
to have $f_i/x_i\ge v$, each cell contributes at least
$w_i(f_i-vx_i)^2$ and $w_i f_i(f_i-vx_i)$ to the two costs.
Sorting the ratios $f_i/x_i$ makes this subset a suffix. Its three
stored sums are $\sum w_i f_i^2$, $\sum w_i f_ix_i$ and
$\sum w_i x_i^2$. They give both lower costs by evaluating a quadratic
and a linear polynomial in $v$. The sort and threshold use outward
rounding, so every selected interval meets the stated inequality.
The complement contributes a nonnegative cost.

The checker takes the maximum of each table lower bound with the
corresponding origin-window bound (PT2). It uses their combination
$F_1+F_2^2$ to bisect for an upper amplitude cap. Each bound tested
in bisection is sufficient by monotonicity of the underlying exact
costs. The same lower bounds at an amplitude-bin left endpoint apply
throughout that bin.

## Retaining the forced component of the remainder

Let $u_*=\alpha/A$ be the trial label direction. Decompose the
projected remainder in the mixed weighted space as

\[
             g-h=A\theta u_*\phi+k,
            \qquad k\perp u_*\phi .                          \tag{FR2}
\]

The projection inequalities in (PT1) give $\theta\ge F_2(A)\ge0$.
For lower clip costs $f_1,f_2$ on an amplitude bin and upper ratio
$r_+$, they also give

\[
  \theta\le\theta_+=(r_+^2-f_1)^{1/2},\qquad
  \|k\|_J/A\le r_\perp=(r_+^2-f_1-f_2^2)^{1/2}.             \tag{FR3}
\]

A bin with negative radicand is infeasible. When using the smaller
orthogonal budget $r_\perp$, the bounds keep the axial part
$A\theta u_*\phi$ explicitly in the reference. In particular the
upper raw second moment takes the minimum with

\[
 \bigg[\sqrt{\mathbb E_{\mu_1}\min\{A_+^2\phi^2,2/3\}}
   +A_+\theta_+\sqrt{M_\phi}+A_+\sqrt{B_1}\,r_\perp\bigg]^2,  \tag{FR4}
\]

where $M_\phi$ bounds the trial's raw tilted second moment. This is
the triangle inequality for the clipped reference, the axial part
and the orthogonal remainder. The preceding bounds (PT4) and (SR19)
remain available.

For the median alternative, the axial-adjusted reference
$\pi_\Delta(\bar p)+A\theta u_*\phi$ has the same coordinate order
as $\bar p$, since both summands have that order. Its median coordinate
is at most $1/3+A(1+\theta_+)|\phi|/\sqrt6$. For a positive axial
component use that projection only subtracts a threshold, and for a
negative component the median is at most $1/3$. Thus (PT6) applies
with its guard amplitude increased to $A_+(1+\theta_+)$ and with
the orthogonal error $r_\perp$. Inflating the amplitude in the
quadratic cost and in the cubic purity estimate retains upper median
cost and lower purity bounds. The variance bounds in (PT5) still
use the original residual ratio. The checker compares this alternative
with the one in the projected-trial checker.

## A separate central-profile reserve

Let $G_0$ be the radial set where the actual masses lie in
$[1/4,1/2]^3$. By (WV11), the positive Taylor slack can contribute
an extra

\[
 \frac3{100}\frac{\rho\pi(1-\rho)}{2(1+\rho)^2}
       \mu_0(G_0)\operatorname{Var}_{\mu_0}(p)                \tag{FR5}
\]

along with either preceding remainder bound. At a pair with one
member in $G_0$, every same-dominant-label product used in those
bounds vanishes. Thus (FR5) uses positive Taylor slack left by
either alternative. [One-sided central reserve](THREE_CELL_WEIGHTED_ANGULAR_VARIANCE.md)

The checker bounds the trial's $\mu_0$ variance from the anchor core
$[0,4]$. The density in anchor $x$ is proportional to
$x^2e^{-(1+\rho)x^2/d_u}$ and its normalizing factor is
$4[(1+\rho)/d_u]^{3/2}/\sqrt\pi$.
For an anchor cell with left endpoint $a$, its density relative
to the auxiliary weighted measure is at most

\[
 \frac{8(1+u)^{3/2}}{\sqrt{\pi d_u}\,Q(a)}e^{-la^2/d_u}.      \tag{FR6}
\]

Here $Q$ is increasing. The same expression at zero without the
exponential bounds the density on the whole half-line. With it,
(SR15) yields a lower variance, still using the original ratio $r_+$.

Set $b_0=1/(12\sqrt{2/3})$. For a chosen $0<q<1$, select the anchor
cells where $A_+(1+\theta_+)|\widehat\phi|\le q b_0$ is certified.
The axial-adjusted reference lies at distance at most $q b_0$ from
$p_*$ on these cells, by nonexpansiveness of the projection. Any
centered mass vector at distance at most $b_0$ from $p_*$ has each
coordinate in $[1/4,5/12]$. A selected cell can fail to be central
only where $|k|>(1-q)b_0$. Consequently

\[
 \mu_0(G_0)\ge m_S-B_S\frac{A_+^2r_\perp^2}{(1-q)^2 b_0^2},
\]

where $m_S$ is a lower mass of the selected cells and $B_S$ is the
maximum of the upper density bounds (FR6) on those cells. The checker
takes the best lower bound over $q=1/10,\ldots,9/10$ and zero.

## Directed cover and discrete transfer

The checker uses the slab spectral block of the excess-purity proof
with 120-bit Arb, 168 kernel degrees, 54 row terms, and 4096 radial
panels through radius 14. On each of 128 residual intervals it takes
the maximum of the quartic and median alternatives after subtracting
the same fixed variance credit. When required, the 16 amplitude bins
also use (FR1)–(FR5). The central reserve is added only after taking
the maximum. All 22 rational slabs of width $1/4000$ satisfy the
strict comparison.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_full_radius_trial_extension.py

For every fixed $0<\rho\le2611/4000$ and $\varepsilon>0$, there exist
$D$ and $\tau>0$ such that any balanced $f:[3]^m\to\Delta_3$ with
$\operatorname{Inf}_i^{\le D}(f_j)\le\tau$ for all $i,j$ satisfies
$S_\rho(f)\le\lim_{k\to\infty}S_\rho(\operatorname{PLUR}_{k,3})+\varepsilon$.
The noise independently retains each coordinate with probability
$\rho$ and otherwise resamples it uniformly. Balance is
$\mathbb E f=\mathbf1/3$. The invariance theorem preserves this balance
in a Gaussian realization, and the discrete plurality limit has the
simplex stability.
[Isaksson–Mossel, Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362)
