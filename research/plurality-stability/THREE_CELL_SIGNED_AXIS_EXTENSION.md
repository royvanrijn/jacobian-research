# Signed-direction clipping for the adverse radial trial

The exact simplex projection gives different clipping costs on the
two signs of a centered radial trial. The formulas below bound the
two costs and the resulting clipped second moment.

The attempted cover of $[2611/4000,2641/4000]$ had an invalid amplitude
bound. The amplitude-loop index overwrote the origin weight used on
later residual bins. The [corrected checker](check_three_cell_signed_axis_extension.py)
rejects the slab $[1307/2000,523/800]$. The signed-direction estimates
remain available as analytic inputs. The [full-radius result](THREE_CELL_FULL_RADIUS_TRIAL_EXTENSION.md)
uses its own unchanged amplitude bound.

## Exact clipping along a label direction

As in (FR2), for $A>0$ let $v=\alpha/A$ be the unit centered label
direction of the adverse profile. Relabel so $v_1\ge v_2\ge v_3$ and
write $a=-v_3$, $b=v_1-v_2$. Thus $3a^2+b^2=2$ and $0\le b\le3a$.
For $s\ge0$ put
$D_v(s)=\operatorname{dist}(p_*+sv,\Delta)^2$ and
$h_s=\pi_\Delta(p_*+sv)-(p_*+sv)$. The two simplex faces in (AP1)
give

\[
 D_v(s)=\frac16(3as-1)_+^2+\frac12(bs-1)_+^2,\qquad
 -v\cdot h_s=\frac{D_v'(s)}2
                 =\frac a2(3as-1)_++\frac b2(bs-1)_+.        \tag{SA1}
\]

The derivative also follows from the gradient of squared distance
to a closed convex set. For $-v$, the corresponding parameters are
$a'=(a+b)/2$ and $b'=(3a-b)/2$.
[Exact simplex projection](THREE_CELL_FINITE_AMPLITUDE_PROJECTION.md)

In particular the clipped reference's squared centered norm is

\[
 \begin{split}
 |\pi_\Delta(p_*+sv)-p_*|^2
 &=s^2-sD_v'(s)+D_v(s)\\
 &=s^2-\frac16((3as)^2-1)_+-\frac12((bs)^2-1)_+ .          \tag{SA2}
 \end{split}
\]

As a function of $s$, this expression is nondecreasing. Its coefficient
of $s^2$ on any interval of activity is nonnegative since
$\frac16(3a)^2+\frac12 b^2=1$. These formulas allow the positive and
negative parts of the radial trial to charge different clipping costs.

## Directed direction and radius tables

The sorted arc of centered unit vectors can be written

\[
 v(\theta)=\cos\theta\,\frac{(2,-1,-1)}{\sqrt6}
         +\sin\theta\,\frac{(0,1,-1)}{\sqrt2},
                        \qquad 0\le\theta\le\pi/3 .          \tag{SA3}
\]

On this arc, $a=\cos\theta/\sqrt6+\sin\theta/\sqrt2$ and
$b'=\sqrt2\sin\theta$ are increasing, while
$b=3\cos\theta/\sqrt6-\sin\theta/\sqrt2$ and
$a'=2\cos\theta/\sqrt6$ are decreasing. Divide the arc into 24
intervals of equal angular length and use lower endpoints of the
enclosures of all four parameters on each bin. Both lower costs in
(SA1) increase in the indicated nonnegative parameters. The removed
second-moment terms in (SA2) increase as well.

Let $F(x)=x\widehat\phi(x^2)$ in the anchor coordinates of (FR1).
On each of 2400 intervals in $[0,4]$ the checker encloses its sign and
bounds $|F|$ below by $f_i$, the right endpoint above by $x_i$, and
the integral of $Q(x)e^{-\beta_u x^2}/(2d_u)$ below by $w_i$.
An interval of uncertain sign is omitted from the signed lower bounds.
When the sign is positive use $(a,b)$, and when it is negative use
$(a',b')$. Here $w_i$ is unrelated to the ancillary density weights
used for the median bound.

For each coefficient $k=3a,b,3a',b'$ the resulting clip terms have
the polynomial form

\[
   \sum_{i} w_i(kf_i-v_0x_i)_+^2,
   \qquad \sum_i w_i f_i(kf_i-v_0x_i)_+,
   \qquad v_0=\sqrt{d_l/d_u}/A.                            \tag{SA4}
\]

Indeed, original mass amplitude $s=A|\phi_\rho|$ becomes
$A\sqrt{d_u/d_\rho}|\widehat\phi|$ at the anchor, and the scalar
weight contributes the compensating factor $d_\rho/d_u$.
Using $d_l$ gives an upper threshold, so (SA4) is a lower bound.
The coefficients $1/6,1/2$ from (SA1) multiply the first terms and
$a/2,b/2$ or $a'/2,b'/2$ multiply the second terms. The resulting
bounds are the two expenses $\|h\|_J^2/A^2$ and
$-\langle u_*\phi,h\rangle_J/A$ from (PT1). They can be combined
by taking the maximum of each expense with its direction-free bound
in (FR1). Sorting $f_i/x_i$ separately for the two signs gives the
same three suffix sums as in the full-radius trial checker.

For the raw clipped second moment, use (SA2) at the upper amplitude
endpoint $A_+$ and the anchor magnitude $|\widehat\phi|$. This gives
an upper bound for the smaller physical amplitude and magnitude by
monotonicity. From its raw upper bound $A_+^2 M_\phi$, subtract lower
terms from (SA2) with the lower direction parameters. In a cell of
known sign, use $f_i/x_i$ as a lower trial magnitude and a lower
bound on its $\mu_1$ mass. Sorting $f_i/x_i$ again gives the two
suffix sums of mass and mass times squared lower trial magnitude.
They evaluate the sum of the active terms $((kA_+ f_i/x_i)^2-1)_+$.
This procedure gives a certified upper clipped moment without a
tail assumption about the signed distribution, since the omitted
tail terms only reduce the actual clipped norm.

## Finite feasibility test

On each residual bin, the checker first tests the nonprojected
excess-quartic bound for all directions. For an unresolved bin it
uses the directional costs (SA4) in each of the 24 angular bins
to bisect the amplitude cap. It then divides the possible amplitudes
into 12 bins. Lower clip costs at the left amplitude endpoint bound
the remaining axial and orthogonal budgets by (FR2)–(FR3). A bin
whose clip costs exceed the residual budget is infeasible. Otherwise
the upper projected raw moment from (SA2) and the two residual
budgets give (FR4), followed by the excess-quartic estimate (EP4).

For an unresolved bin the checker also evaluates the median
alternative of (PT5)–(PT6) on the axial-adjusted reference. This
uses the universal median coordinate bound $|v_2|\le1/\sqrt6$.
The median and quartic estimates are alternatives, with the
unconditional variance credit already subtracted. The original
one-sided central reserve (FR5) can then be added to their maximum.
The computed scalar upper bound must be strictly smaller than this
certified remainder on every feasible bin. The case $A=0$ follows
from the negative scalar complement. Uniform circular masses reduce
equality to the equal-circle rigidity statement.

The attempted cover uses 30 rational correlation slabs of width $1/4000$,
128 residual bins and the direction and amplitude subdivisions stated
above. Spectral data use the same 120-bit Arb, 168 kernel modes,
54 row terms and radial moment quadrature as the preceding result.
The replay command is

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_signed_axis_extension.py

## Projection along a simplex ray and exact axial expense

Use (FR2) from the [full-radius comparison](THREE_CELL_FULL_RADIUS_TRIAL_EXTENSION.md) and (SA1). The original trial in the centered
label space is \(A\phi v\) with \(|v|=1\) and \(\|\phi\|_J=1\).
The original remainder is orthogonal to \(v\phi\). Therefore the
coefficient of the projected remainder is exactly
\[
  \theta=-\langle h,v\phi\rangle_J/A,\qquad
  h=\Pi_{\Delta}(p_*+A\phi v)-(p_*+A\phi v).
 \tag{AR1}
\]
The one-sided expense in (FR3) is a lower bound for this scalar. It also
has a useful upper bound that depends on the trial alone.

### A projection identity

If \(v_1\ge v_2\ge v_3\), the simplex projection of \(p_*+tv\) has
coordinates \((1/3+tv_i-\lambda(t))_+\) with
\[
 \lambda(t)=\max_{k=1,2,3}
  \left[ \frac13-\frac1k+t\frac{\sum_{i\le k}v_i}{k}\right].
\]
For \(0<t<s\), convexity and the active value at \(t\) give
\(\lambda(s)-\lambda(t)\ge(s-t)\bar v_k\), where \(k\) is the number
of positive projected coordinates at \(t\) and \(\bar v_k\) is the
average of the \(k\) greatest coordinates. For every previously zero
coordinate \(i>k\), \(v_i\le\bar v_k\), so the threshold
\(\lambda(s)-\lambda(t)\) keeps that coordinate zero after applying
\((s-t)v\) to the projected point. On the other coordinates this
threshold yields precisely \((1/3+sv_i-\lambda(s))_+\). These
coordinates sum to one. The same conclusion at \(t=0\) follows from
\(\Pi_\Delta(p_*)=p_*\). Thus
\[
  \Pi_\Delta\!\bigl(\Pi_\Delta(p_*+tv)+u v\bigr)
                       =\Pi_\Delta(p_*+(t+u)v)
                       \qquad(t,u\ge0).
 \tag{AR2}
\]

Now write the actual profile as
\(p=\Pi_\Delta(p_*+A\phi v)+ A\theta\phi v+k\), as in (FR2).
Since \(p\) is feasible, projecting \(p\) is the identity. Apply
projection nonexpansiveness and (AR2) with \(v\) or \(-v\) according
to the sign of \(\phi\), to obtain the pointwise bound
\[
 \left|p-\Pi_\Delta(p_*+A(1+\theta)\phi v)\right|\le |k|.
 \tag{AR3}
\]
Hence the raw second moment in the \(\mu_1\) feature obeys
\[
 \mathbb E_{\mu_1}|p-p_*|^2
 \le
 \left(\sqrt{\mathbb E_{\mu_1}
           |\Pi_\Delta(p_*+A_+(1+\theta_+)\phi v)-p_*|^2}
       + A_+\sqrt{B_1}\,r_\perp\right)^2 .
 \tag{AR4}
\]
Here \(0\le\theta\le\theta_+\), \(A\le A_+\), and
\(\|k\|_J/A\le r_\perp\), with the constants of (FR3). Monotonicity
of the projected norm in the nonnegative trial amplitude follows
from (SA2). In the median and central comparisons one can also
use the projected enlarged trial as the reference with the same
pointwise error bound. Its median coordinate is at most
\(1/3+A_+(1+\theta_+)|\phi|/\sqrt 6\), as in (FR4).

### An upper bound for the axial expense

For a sorted direction \(v\) put \(a=-v_3\), \(b=v_1-v_2\), and
\(s=A|\phi|\), using the sorted direction of \(-v\) when \(\phi<0\).
For \(s>0\), differentiating the squared clipping distance as in
(SA1) and using \(3a^2+b^2=2\) give
\[
 \frac{v\cdot(\Pi_\Delta(p_*+s v)-p_*)}{s}
 =
 \frac{a}{2}\min(3a,1/s)+\frac{b}{2}\min(b,1/s)
 =:D_v(s).
 \tag{AR5}
\]
This expression is nonnegative, decreases in \(s\), and increases
separately in the two nonnegative parameters on the right.
For \(\phi=0\) its value is unused. Equation (AR1) becomes
\[
 \theta = 1-\int \phi^2
      D_{\operatorname{sgn}(\phi)v}(A|\phi|)\,dJ .
 \tag{AR6}
\]
In particular an upper bound for \(\theta\) needs only a lower
integral on any selected known-sign panels. This bound is
available jointly with (FR3).

In the anchor variable of (FR1), the integral \(\phi^2dJ\) is
\(Q(x)e^{-\beta_u x^2}|F(x)|^2\,dx/(2d_u)\), and the physical
\(|\phi|\) is at most \(|F(x)|/x\). On each known-sign panel
with positive lower endpoint take a lower bound \(w_j\) for
its \(\phi^2dJ\) mass and an upper bound \(m_j\) for \(|F(x)|/x\).
Then the retained integral is bounded below by the sum of
the panel masses \(w_j\) times the right side of (AR5) at
\(A_+m_j\) and the certified lower directional parameters.
The contribution of any other panel is nonnegative. Sorting
the \(m_j\)'s within each sign gives a prefix table of \(w_j\)
and a suffix table of \(w_j/m_j\), so each term of type
\(\sum_j w_j\min(k,1/(A_+m_j))\) uses two cumulative evaluations
with outward threshold rounding. Taking the minimum of the
resulting upper bound with the earlier FR3 upper bound gives
the tighter \(\theta_+\) used in (AR4).

### Retaining the shifted clipping charge

A separate bound is obtained by selecting any \(\xi\ge0\) and recentering
the pointwise trial at \(y_\xi=p_*+A(1+\xi)\phi v\). Orthogonality of
the original trial residual to \(\phi v\) gives
\[
 \|p-y_\xi\|_J^2/A^2\le r^2+\xi^2,\qquad r\ge\|g\|_J/A.
\]
Let \(C(A)\) be the exact normalized squared clipping distance
of the unshifted trial at amplitude \(A\). The projection
inequality yields
\[
 \frac{\|p-\Pi_\Delta(y_\xi)\|_J^2}{A^2}
   \le r^2+\xi^2-(1+\xi)^2C((1+\xi)A).                      \tag{AR7}
\]
Thus any finite list of fixed shifts improves the feasibility test:
a bin is empty when any right side is negative. The same right
side also bounds the remainder norm in the analogue of (AR4)
with the \((1+\xi)A\) clipped reference instead of the
\((1+\theta)A\) reference. The density ratio \(B_1\) converts the
\(J\)-norm to the \(\mu_1\)-norm. On an amplitude bin the normalized
distance \(C\) is nondecreasing, as seen from the signed hinge
formula (SA1), so its lower endpoint supplies the lower distance
in (AR7). The upper amplitude endpoint supplies the clipped
raw moment upper bound via (SA2).

The first-tangent estimate in FR3 is a special case after convex
linearization. The normalized distance as a function of \(\xi\)
has derivative \(2\theta\) at zero, so evaluating (AR7) at
\(\xi=\theta\) after keeping only the tangent incurs the old
charge \(C(A)+\theta^2\). The exact shifted distance keeps the
nonnegative convex remainder. Finite shifts chosen using any
fixed lower cost data on a bin are valid without solving the
dual optimization problem.

### A scalar dual for the tilted raw moment

Let \(C=\Delta-p_*\), let \(q(t)\in C\) be any measurable reference,
and assume \(d\mu=m\,dJ\) with nonnegative integrable \(m\). For
\(b\ge0\), the raw tilted moment under a residual budget has the upper bound
\[
 \sup_{\substack{z(t)\in C\\ \|z-q\|_J\le b}}\int m|z|^2\,dJ
 \le
 \inf_{\lambda\ge0}\left\{\lambda b^2+
   \int \sup_{z\in C}\bigl[m(t)|z|^2-\lambda|z-q(t)|^2\bigr]\,dJ(t)\right\}.
 \tag{AR8}
\]
The bound follows by adding \(\lambda(b^2-\|z-q\|_J^2)\ge0\)
to the objective at each feasible \(z\), then taking the
pointwise supremum and the infimum over \(\lambda\).

For \(\lambda>m(t)\) the inner optimizer is the projection onto
\(C\) of \(\lambda q/(\lambda-m(t))\), by completing the
square. For \(\lambda\le m(t)\) it is a simplex vertex with
the largest \(q_i\), since the objective is convex and the
vertices have the same norm. These prescriptions also give
measurable selections, resolving ties by the first index.

If \(q=\Pi_C(s v)\) for a fixed direction and \(s\ge0\), the
inner value in (AR8) is nondecreasing in \(s\), so an
upper amplitude endpoint gives an upper integrand on each bin.
For the projection case, first use
\[
 \Pi_C\bigl(\kappa\Pi_C(sv)\bigr)=\Pi_C(\kappa s v),
                  \qquad \kappa\ge1.                         \tag{AR9}
\]
Indeed, on the original positive support the projection
threshold for the rescaled clipped vector is
\(\lambda(\kappa s)-\kappa\lambda(s)\), where \(\lambda(s)\)
is the threshold from (AR2). This difference is nonnegative
by convexity and \(\lambda(0)=0\). Previously inactive
components of the rescaled clipped simplex vector are
\((1-\kappa)/3\le0\), so the same threshold leaves them zero.
Thus the optimizer is a further projected point on the same
ray. Along the ray each top-\(k\) cumulative coordinate sum is
nondecreasing: whenever the active size is \(j\ge k\), its
derivative is \(\sum_{i\le k}(v_i-\bar v_j)\ge0\), and
whenever \(j<k\) it is constant. Hence the optimizer majorizes
\(q\) and stays in its support. For the vertex case a top
vertex likewise majorizes \(q\) and stays in its support.
On an interval of fixed active support, the derivative of the
optimal pointwise value in \(s\) is \(2\lambda(z^*-q)\cdot q'\).
The derivative \(q'\) has coordinates \(v_i-\bar v_j\) on the
support and is decreasing there, so this dot product is
nonnegative by summing over the successive top-coordinate
partial sums of \(z^*-q\). Continuity across the support
changes proves the claim.

The application takes \(q=p_{\mathrm{ref}}-p_*\) from (AR3)
or (AR7), and \(b\) from its residual bound. Thus the scalar
dual keeps the full density \(m\) and the exact simplex
constraint instead of replacing either by a uniform bound.
