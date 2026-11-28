# Signed shifted clipping in the variance bootstrap

The signed projection costs sharpen the radial variance bootstrap
when used both for the clipped second moment and for the amplitude
cap. The amplitude cap uses shifted clipping in the adverse trial
direction. The resulting stages use the full quartic moment bound
while retaining the actual first angular deficit. The directed
cover proves the balanced three-cell comparison through $671/1000$.

## Feasible amplitudes

Use the notation of the [variance bootstrap](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md).
At each stage let $H=L^2(J'_m\,d\nu)$, let $u$ denote the unit adverse
radial trial and write $z=p-\mathbf1/3=Avu+g$, where $A=|\alpha|$,
$v=\alpha/A$ for $A>0$ and $g\perp vu$. A profile requiring a higher
gain has $\|g\|_H\le Ar$ by (SV7) with the radial reserve. Put
$C=\Delta_3-\mathbf1/3$ and $h=\Pi_C(Avu)-Avu$. For fixed $v$ define
\[
 I_2(A)=A^{-2}\|h\|_H^2,\qquad
 I_1(A)=-A^{-1}\langle vu,h\rangle_H .
\]
Both are nonnegative and nondecreasing in $A$. The signed
formulas (SA1) give this directly after dividing the active terms
by $A$. Projection and orthogonality imply
\[
 I_2(A)+I_1(A)^2\le r^2.
\]
For every $\eta\ge0$, (SV15) also gives
\[
 (1+\eta)^2 I_2((1+\eta)A)-\eta^2\le r^2.                \tag{SS1}
\]
Indeed $z\in C$ almost everywhere and $\|z-(1+\eta)Avu\|_H^2
=\|g\|_H^2+\eta^2A^2$. Minimizing this squared distance pointwise
over $C$ gives the inequality. The upper amplitude bound is found
by rejecting amplitudes that violate any of these inequalities.
The checker uses the five shifts $1/25,2/25,3/25,4/25,1/5$.
Each test is nondecreasing in $A$. The shifts are only used to
restrict amplitude. They do not enter the tilted residual budget.
[Signed projection costs](THREE_CELL_SIGNED_AXIS_EXTENSION.md)

## Directed radial and direction bounds

Relabel $v$ to the sorted arc (SA3) and split it into 16 equal
angular intervals. On each interval take the lower endpoints of
the two increasing and two decreasing functions in (SA3), with a
nonnegative lower bound for each parameter. These bounds give
pointwise lower bounds on each of the positive-part terms in
(SA1) and on the subtracted terms in (SA2). Use the parameters
of $v$ when $u(x^2)>0$ and of $-v$ when $u(x^2)<0$.

Let $Q(x)$ be the positive modified row polynomial of (SV8),
$d=1-\rho^2$, $\kappa=\rho^2/d$, $\beta=1+\kappa$ and
$\mathrm{tilt}=1/Z_1$. The trial on the actual correlation slab is
\[
 f(x)=x u(x^2)=\frac{4d[P(x)-c e^{\kappa x^2}]}{\sqrt D\,Q(x)},
\]
with $P$ the polynomial of the chosen finite feature vector,
$c$ its balance correction and $D$ its squared norm. The weights in
$x=\sqrt t$ are
\[
 J'_m(x^2)d\nu(x^2)=\frac{x^2}{2d}Q(x)e^{-\beta x^2}dx,
 \qquad d\mu_1(x^2)=2\,\mathrm{tilt}\,x^2e^{-\beta x^2}dx.
\]
On each of 1200 panels in $[0,4]$ whose trial sign is determined,
take a lower bound $f_i$ on $|f(x)|$, an upper bound $x_i$ on $x$,
a lower bound $w_i$ on the panel integral of $Q(x)e^{-\beta x^2}/(2d)$
and a lower bound $m_i$ on its $\mu_1$ mass. Define these bounds
separately for the two signs. For a direction coefficient
$k\in\{3a,b,3a',b'\}$, with its associated sign, the contributions
to $I_2(A)$ and $I_1(A)$ are bounded below, respectively, by the
corresponding (SA1) multipliers times
\[
 \sum_i w_i(k f_i-x_i/A)_+^2,\qquad
 \sum_i w_i f_i(k f_i-x_i/A)_+ .                           \tag{SS2}
\]
In particular the $x^2$ in the first weight cancels the origin
pole in $u$. The numbers in a table are fixed lower or upper
bounds for the entire correlation slab. For each sum use only
indices with a lower bound on $f_i/x_i$ above an upper bound on
$1/(kA)$. Sorting by that lower bound reduces (SS2) to three
suffix sums. Omitting a panel of uncertain sign, an uncertain
threshold, or the region $x>4$ is safe for these lower bounds.

On an amplitude bin $[A_-,A_+]$, evaluate the lower costs of
(SS2) at $A_-$ and call them $l_2,l_1$. The clipped profile
$z_*=\Pi_C(Avu)$ satisfies $\|g-h\|_H^2/A^2\le r^2-l_2$ and
$\langle g-h,vu\rangle_H/A\ge l_1$. Use (SV13) on the remaining
tilted residual. In particular it is infeasible when
$r^2-l_2<l_1^2$. Equation (SA2) makes $|z_*|^2$ nondecreasing along
each ray, so its $\mu_1$-integral at $A_+$ bounds the whole bin.
Starting from $A_+^2 M$, where $M$ bounds $\int u^2d\mu_1$,
subtract the (SA2) multipliers times
\[
 \sum_i m_i\big[(kA_+ f_i/x_i)^2-1\big]_+ .                \tag{SS3}
\]
This uses the same signed ordering and the suffix sums of $m_i$
and $m_i(f_i/x_i)^2$. The result is an upper bound for the clipped
second moment on the bin. Independently (SV12) supplies another
upper bound for the same moment using the smaller of the panelwise
raw trial moment and the simplex-radius-squared times the panel
mass, with the explicit chi-square tail after $x=4095/1024$. The minimum
of the two upper bounds is valid.

For the untruncated raw moment $M$ the checker integrates the
directed interval enclosures on 4096 panels in $[0,12]$. On the
remaining half-line let $R=12$, $n=\deg P$ and
$S=\sum_i|P_i|R^i$. Then $|P(x)|\le S(x/R)^n$. Separate the two
terms in $P(x)-c e^{\kappa x^2}$ using $(a-b)^2\le2a^2+2b^2$,
and use $Q(x)\ge Q_jx^j$ from any positive coefficient. For
$p=0,2n$ and the respective Gaussian decay parameter
$\zeta=1-\kappa,\beta$, the resulting monomial integral obeys
\[
 \int_R^\infty\!x^{p-2j}e^{-\zeta x^2}dx
 \le \frac{R^{p-2j-1}e^{-\zeta R^2}}
 {2[\zeta-\max(p-2j-1,0)/(2R^2)]} .                       \tag{SS4}
\]
To see this, put $q=x^2$ and use
$\log(q/R^2)\le(q-R^2)/R^2$ for the power of $q$ in the
integrand. The checker verifies both denominators are positive.

Finally combine the clipped second moment with (SV13) and take
the minimum with (SR19), which uses the same untruncated $M$.
Apply both tilted mean bounds in (SV11). The gain bound (SV9)
follows on each of 32 amplitude bins and each direction interval.
The minimum over feasible bins is valid for a subsequent stage.
Every stage uses the entire new gain in the same variance penalty.
[Scalar remainders](THREE_CELL_SHARP_SCALAR_REMAINDERS.md)

## Verified correlation cover

The [checker](check_three_cell_signed_shift_bootstrap.py) and
[complete receipt](three_cell_signed_shift_bootstrap_verification.json)
cover $[329/500,671/1000]$ in 290 adjacent rational slabs. The
largest full operator norm upper bound is below $1.998836$.
The checker starts each slab from gain $12/25$, uses a preliminary
radial reserve of $10^{-5}$ and requires a final norm strictly
less than $1999/1000$. A rejected slab is subdivided from the initial
width $1/16000$. The spectral stage is the unchanged directed stage
of the [tilted-orthogonality bootstrap](check_q_directmean_bootstrap_spent.py),
with 110-bit Arb, 192 kernel modes, 96 row terms and 512 cubic
Hermite radial panels. Only the adverse trial step uses (SS1)–(SS4).
The signed projection identities enter through freshly computed
actual-correlation tables. Neither the earlier signed-axis amplitude
checker nor the graded or localized checker is a computation source.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_signed_shift_bootstrap.py --workers 8

The adjoining verified lower interval and the existing fixed-mass
dimension reduction give simplex optimality and equality rigidity for
balanced three-cell Gaussian partitions in every dimension $d\ge2$
and every $0<\rho\le671/1000$. For the newly covered interval the
retained radial reserve and first angular deficit give the planar
transfer-cost margin by (FA27), using
$w(t)=\sqrt t\,e^{-t-4\sqrt t}$ because $\kappa<1$ and
$\gamma\sqrt2<4$. Equation (CT5) gives a positive unevaluated right
neighborhood of the endpoint.
[Planar transfer cost](THREE_CELL_FINITE_AMPLITUDE_COMPARISON.md)
[Correlation continuation](THREE_CELL_CORRELATION_CONTACT.md)

For each fixed $0<\rho\le671/1000$ and $\varepsilon>0$, there are
$D\in\mathbb N$ and $\tau>0$ such that the following holds uniformly
in $n$. A simplex-valued rule $F:[3]^n\to\Delta_3$ with uniform mean
$\mathbf1/3$ and $\operatorname{Inf}_i^{\le D}(F_j)\le\tau$ for all
$i,j$ has stability under retain-and-resample noise of correlation
$\rho$ at most the infinite-voter plurality limit plus
$\varepsilon$. The transfer uses the balanced ordinary Gaussian
partition in Theorem 7.1 and the plurality limit in Lemma 7.3.
[Isaksson–Mossel](https://arxiv.org/pdf/0903.3362)

