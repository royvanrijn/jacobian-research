# An unrestricted three-cell comparison through correlation 2589/4000

For every ordinary balanced Gaussian three-cell partition in dimension at
least two, simplex stability is maximal for $0<\rho\le2589/4000$.
On the added interval $[513/800,2589/4000]$, equality requires a simplex
up to rotation, relabeling and null sets. The previous comparison covers
correlations through $513/800$ and supplies the reduction to the plane.
[Excess-purity comparison](THREE_CELL_EXCESS_PURITY_EXTENSION.md)

The [directed checker](check_three_cell_projected_trial_extension.py)
covers the added interval with internal verification. It reuses the
anchored finite operator and the unconditional degree-one variance credit
of the preceding comparison. The new steps are Euclidean projection of
the adverse trial and an alternative median-label bound.

## The two expenses forced by clipping

For a fixed profile write $z=p-p_*$ with $p_*=\mathbf1/3$ and
$z=\alpha\phi+g$ in the mixed weighted space of (EP2). Use its norm
$\|\cdot\|_J$, and put $A=|\alpha|$, $G=\|g\|_J$, $r=G/A$
when $A>0$. Let $\bar p=p_*+\alpha\phi$, $\pi_\Delta$ be the Euclidean
projection onto the probability simplex and

\[
                h=\pi_\Delta(\bar p)-\bar p .
\]

Patchwise projection gives $\|g\|_J^2\ge\|h\|_J^2+\|g-h\|_J^2$.
There is a further lower bound on $\|g-h\|_J$. For $u=\alpha/A$,
the variational inequality of the projection against $p_*$ gives
$\phi\,u\cdot h\le0$. When $A|\phi|>r_*=\sqrt{2/3}$, the circumradius
bound sharpens it to $\phi\,u\cdot h\le-|\phi|(A|\phi|-r_*)$.
Since $g$ is weighted-orthogonal to the trial, pairing $g-h$ with the
weighted unit vector $u\phi$ yields

\[
 \frac{G^2}{A^2}\ge F_1(A)+F_2(A)^2,
 \quad
 F_1(A)=\int J(|\phi|-r_*/A)_+^2d\nu,
 \quad
 F_2(A)=\int J|\phi|(|\phi|-r_*/A)_+d\nu.                 \tag{PT1}
\]

Use the same origin window and constants $X,c,C,j$ as in the preceding
checker. On $0\le x\le X$, put $b=C+r_*/A$ and
$y=\min\{X,c/b\}$. The checker uses the lower bounds

\[
 F_1(A)\ge 2j\int_0^y(c-bx)^2dx,
 \qquad F_2(A)\ge 2j\int_0^y(c-Cx)(c-bx)dx .                \tag{PT2}
\]

The two right sides increase in $A$ and have elementary polynomial
antiderivatives. At $y=c/b$ the second is
$jc^3(3b-C)/(3b^2)$. At $y=0$ use the continuous value zero.
The bisection for the amplitude cap uses the sum in (PT1).

For a certified amplitude bin $[A_-,A_+]$ and an upper residual ratio
$r_+$, let $f_1(A_-)$ be the first lower bound in (PT2). Then

\[
 \|g-h\|_J^2/A^2\le r_+^2-f_1(A_-).                      \tag{PT3}
\]

A bin with negative right side is infeasible. The forced component
$F_2^2$ is used in the amplitude cap, while (PT3) charges only the
projection cost $F_1$ to the available residual norm.

## The clipped raw second moment

Set $r_{\rm eff}=\sqrt{r_+^2-f_1(A_-)}$ on a feasible bin. The projection
contracts distance to $p_*$, so
$|\pi_\Delta(\bar p)-p_*|\le\min(A|\phi|,r_*)$.
If $B_1$ bounds the tilted density relative to the mixed weight, then

\[
 (\mathbb E_{\mu_1}|z|^2)^{1/2}
 \le \left(\mathbb E_{\mu_1}\min\{A_+^2\phi^2,r_*^2\}\right)^{1/2}
               +\sqrt{B_1}\,A_+ r_{\rm eff}.                 \tag{PT4}
\]

This follows from (PT3) and the triangle inequality on the tilted
space. The checker takes the minimum of this upper bound squared and
the unprojected estimate (SR19). It uses the smaller raw bound in
(EP4) and the usual purity coefficient. The lower variance still uses
(SR15) with the original ratio $r_+$.

To enclose the clipped expectation, use the anchor coordinates of
the preceding checker. On $0\le x_0=\sqrt{t_0}\le4$ each subinterval
has an upper tilted mass and an upper bound on the absolute anchor
trial. The physical absolute trial is at most its anchor value.
Sort the subintervals by these upper values. A prefix can use the
square of that bound times $A_+^2$ and the complement the constant
$r_*^2$. Both bounds apply at every radius, so the choice of the
prefix is allowed to be approximate. The tail after four uses
$r_*^2$ times its exact gamma survival probability.

## The median-label alternative

For $n=1,3$, let $\mu_n$ be the normalized tilt
$t^{n/2}e^{-\kappa t}d\nu$. Put

\[
 c_1=1/\pi^2,\quad c_3=7/1000,\qquad
 \omega_n=c_n(2\rho)^n d\Gamma(1+n/2)^2,\qquad
 q_3'=\left(\frac7{54\pi^2}-c_3\right)(2\rho)^3d\Gamma(5/2)^2.
\]

Choose a median coordinate $j$ of $\alpha$ and write
$m_{jn}=\mathbb E_{\mu_n}(p_j-1/2)_+$. The $L=3$ weight in (WV3)–(WV5)
and the leftover cubic quartic term give the Taylor slack bound

\[
 \frac12(\omega_1 V_1+\omega_3 V_3)
   -5(\omega_1m_{j1}^2+\omega_3m_{j3}^2)
       +q_3'\,\mathcal P_3 V_3,                           \tag{PT5}
\]

where $V_n=\operatorname{Var}_{\mu_n}(p)$ and
$\mathcal P_3=\mathbb E_{\mu_3}(1-|p|^2)$.
The absorbed scalar credit is $(12/25)\omega_1 V_1$, so (PT5) leaves
$(1/50)\omega_1 V_1$, the displayed cubic terms and the median error.
[Weighted angular variance](THREE_CELL_WEIGHTED_ANGULAR_VARIANCE.md)

Projection preserves coordinate order. Its coordinate formula is
$(\bar p_i-\lambda)_+$ with $\lambda\ge0$, since $\sum_i\bar p_i=1$.
For either sign of $\phi$, the predicted label $j$ is a median both
before and after projection. Thus $(p_j-1/2)_+\le|g-h|/\sqrt2$:
if $p_j>1/2$ there is a label ranked at least as high in the reference
whose actual mass is at most $1-p_j$. The difference of those two
residual coordinates is at least $2(p_j-1/2)$.

Also $|\alpha_j|\le A/\sqrt6$ by the centered label constraint. Given
$0<\eta<1/6$, on the set $|A\phi|\le\sqrt6(1/6-\eta)$ the predicted
median stays below $1/2-\eta$, including after projection. On that set
use $(p_j-1/2)_+\le(r_*|g-h|-\eta)_+\le|g-h|^2/(6\eta)$.
If $B_n$ bounds the tilted density relative to the mixed weight and
$b_n$ is an upper mass of the exceptional set, this proves

\[
 \frac{m_{jn}}{A}\le
       \frac{B_n A_+ r_{\rm eff}^2}{6\eta}
             +\sqrt{B_n b_n/2}\,r_{\rm eff}.               \tag{PT6}
\]

There is also the bound $\sqrt{B_n/2}\,r_{\rm eff}$ from the unconditional
median inequality. The checker uses the best of these upper bounds from
a fixed finite list of $\eta$, with the same tilt mass partition as in
(PT4). The exceptional mass includes every interval whose upper trial
bound exceeds the threshold and the gamma tail after four.

For $B_3$, let $Q(x)=\sum_{i\ge1}b_i^{\rm mix}x^{i-1}$ be the mixed row
polynomial at the anchor. In the anchor coordinate the density ratio for
$\mu_3$ is at most $4x^2/[d_u^{3/2}\Gamma(5/2)Q(x)]$. The lower bound

\[
 Q(x)/x^2\ge b_3^{\rm mix}+2\sqrt{b_1^{\rm mix}b_5^{\rm mix}}
                            +2\sqrt{b_2^{\rm mix}b_4^{\rm mix}}
\]

gives the value in the checker. The physical trial variance under $\mu_3$
is bounded below by its anchor variance times $d_u/d_l$, using the core
through four for the anchor variance. For $V_3/A^2$ use (SR15) with
$B_3$ and the original ratio $r_+$, and for $\mathcal P_3$ use
$2/3-B_3 A_+^2(1+r_+^2)$ or zero, whichever is larger.

## Directed cover

On each residual bin, the checker compares the scalar upper bound from
(EP2) to the largest certified remaining Taylor cost. It can use the
excess-quartic bound from (EP3)–(EP4), with (PT4) on an amplitude bin, or
the median alternative (PT5)–(PT6). These are lower bounds on the same
remainder after the degree-one credit. The checker takes their maximum.
It uses 128 residual bins and 16 amplitude bins whenever the full-amplitude
bound requires refinement. The actual first angular deficit is available
as a separate nonnegative term.

The cover consists of 24 correlation slabs of width $1/4000$. Spectral
inputs use the 120-bit Arb enclosure, 168 kernel degrees, 54 row terms,
and 4096 radial panels through radius 14 from the preceding checker.
The replay command is

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_projected_trial_extension.py

If $A=0$, the strictly negative complement bound excludes every nonuniform
profile. At a uniform profile equality reduces to the equal-circle
rigidity statement of (DG6).

## Discrete plurality

For fixed $0<\rho\le2589/4000$ and $\varepsilon>0$ there exist $D$ and
$\tau>0$ such that any balanced $f:[3]^m\to\Delta_3$ with
$\operatorname{Inf}_i^{\le D}(f_j)\le\tau$ for all $i,j$ satisfies
$S_\rho(f)\le\lim_{k\to\infty}S_\rho(\operatorname{PLUR}_{k,3})+\varepsilon$.
Balance is $\mathbb E f=\mathbf1/3$. The noise independently retains a
coordinate with probability $\rho$ and otherwise resamples it uniformly.
The invariance theorem gives a balanced Gaussian realization and its
plurality limit is the simplex stability, as in the preceding transfers.
[Isaksson–Mossel, Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362)
