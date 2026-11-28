# Localized shape responses at a tangential simplex tie

Let \(A\) be an ordinary balanced planar three-cell maximizer at
\(0<r<1\). Use the Gaussian interface measure \(\nu\), outward normal
speeds, positive score-gradient weight \(a\), mass flux \(C\), and
interface form
\[
 Q_A(f)=\|f\|_a^2-\langle f,Kf\rangle
\]
from the [maximizer analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds).
On each oriented pair interface write \(g\) for the fixed-shape
correlation derivative of the corresponding score contrast, and
write \(G(f)=\int gf\,d\nu\). Compact regular speeds are dense in
the Gaussian weighted interface space. The form and the flux extend
continuously to that space, and the form is nonnegative on the
balanced space \(Cf=0\).
[Volume-flux continuity](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux).
[Mixed pairing continuity](THREE_CELL_CORRELATION_CONTACT.md#anisotropic-dilation-and-mass-tangent-directions).

Choose finitely many disjoint compact regular interface subarcs
\(I_\ell\), each strictly inside an open analytic two-cell arc
\(\Gamma_\ell\). Require an open buffer on the same analytic arc
that is disjoint from all chosen subarcs. The third adjusted score
is strictly lower on these arcs. Write \(I=\bigcup_\ell I_\ell\),
and let \(\mathcal F_I\) be the balanced weighted speeds supported
on \(I\). Then there is a strict spectral gap on \(\mathcal F_I\),
even when the full balanced Hessian of \(A\) has zero modes.

Suppose further that \(A\) ties the regular tripod \(T\) at \(r\),
that other simplex-optimal correlations accumulate at \(r\), and
that the fixed-shape deficit \(D(s)=S_s(T)-S_s(A)\) satisfies
\(D'(r)=0\) and \(D''(r)>0\). For sufficiently thin fixed normal
tubes about \(I\), every ordinary balanced partition agreeing
with \(A\) outside those tubes loses to \(T\) at all nearby
scalar correlations other than \(r\). The comparison includes
arbitrary measurable changes inside the tubes.

## A strict gap on buffered interface pieces

The unrestricted stability inequality at the maximizer is applied
only to balanced speeds. If \(f\in\mathcal F_I\) has \(Q_A(f)=0\),
polarization gives \(Q_A(f,w)=0\) for all balanced \(w\), including
speeds supported elsewhere. Because the global flux \(C\) is onto
the centered label space, the functional \(Q_A(f,\cdot)\) then
factors through \(C\). Equivalently there are label constants
\(\lambda_i\) with
\[
              a_e f_e-(Kf)_e=\lambda_i-\lambda_j
                    \quad(e=(i,j))                           \tag{LR1}
\]
almost everywhere on the interfaces. One can obtain the
factorization using two fixed regular mass-correcting bumps on
a spanning tree.

For a compactly supported speed \(f\), the function \((Kf)_e\)
extends analytically in the observation point across the full
plane, since it is a Gaussian smoothing of a compactly supported
signed surface measure. On the buffer of each \(\Gamma_\ell\)
we have \(f=0\), so (LR1) makes \((Kf)_e\) a constant there.
Analytic continuation along \(\Gamma_\ell\) gives that same
constant on \(I_\ell\). Thus \(a_e f_e=0\) on \(I_\ell\) for
every \(\ell\). The positive \(a_e\) gives \(f=0\).

The same continuation observation applies to noncompact zero
speeds. Indeed the bound \(\int a^{-1}d\nu<\infty\) from (IB1)
and weighted Cauchy–Schwarz make the signed surface measure of
any weighted speed finite. On a bounded complex observation
set the Gaussian kernel and its derivatives are bounded in the
source point by Gaussian damping. Thus \(Kf\) remains analytic.
Equation (LR1) and the positive analytic weight \(a\) give an
analytic representative of \(f\) on every connected regular
arc. Such a speed either vanishes on the whole arc or has
only isolated zeros there.

On \(I\), the weight \(a\) is bounded above and below by positive
constants. The interaction on that compact support is compact
in the weighted norm, since its surface kernel is bounded.
If a sequence of unit speeds in \(\mathcal F_I\) had forms tending
to zero, a weakly convergent subsequence would have interaction
converging to one. Its limit would have norm at most one and
nonnegative form, hence norm one and zero form. This contradicts
(LR1). Consequently
\[
 Q_A(f)\ge c_I\|f\|_a^2\quad(f\in\mathcal F_I),\qquad
 \Gamma_I=\max_{f\in\mathcal F_I}\{2G(f)-Q_A(f)\}
                                                               \tag{LR2}
\]
for some \(c_I>0\), and the maximum is attained.

## Strict loss of localized response at a tie

The [tangential contact inequality (G10)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)
and weighted approximation with fixed mass bumps give
\[
     2G(f)-Q_A(f)\ \le\ \tfrac12 D''(r)
                   \quad\text{for every balanced speed } f .
                                                               \tag{LR3}
\]
If \(\Gamma_I=D''(r)/2\), its maximizing compactly supported
speed \(u\) would maximize this functional on the full balanced
space. Testing its first derivative makes \(Q_A(u,\cdot)-G\)
a mass-flux functional. The analogue of (LR1) is
\[
              a_e u_e-(Ku)_e-g_e=\lambda_i-\lambda_j .
                                                               \tag{LR4}
\]
Both \((Ku)_e\) and \(g_e\) are analytic in the observation point
along the carrying arc. For \(g_e\) this follows by differentiating
the Gaussian kernel at the fixed \(r\), with a locally uniform
integrable bound. On every buffer, (LR4) makes their sum constant.
Analytic continuation through \(I_\ell\) therefore gives \(u=0\).
This contradicts the positive value \(D''(r)/2\). Hence
\[
                         \Gamma_I<\tfrac12D''(r).
                                                               \tag{LR5}
\]
In particular no differentiable moving competitor with a nonzero
speed of this support can saturate the tangential shape bound.
For a general saturated weighted response, (LR4) yields the
same analytic representative on each regular arc.

## Measurable competitors in thin tubes

In arc-length and signed normal-distance coordinates \((s,u)\)
near \(I_\ell\), let its tube be the image of
\(I_\ell\times[-\delta,\delta]\). Such tubes can be chosen
disjoint with regular normal coordinates up to their boundaries.
Write \(K_\delta\) for their union and let \(\mathcal C_\delta\)
consist of all ordinary balanced partitions \(P\) that agree
with \(A\) outside \(K_\delta\). The local response has the
exact value
\[
 \lim_{t\to0}
 \frac{\sup_{P\in\mathcal C_\delta} S_{r+t}(P)-S_{r+t}(A)}{t^2}
                       =\Gamma_I
                                                               \tag{LR6}
\]
for every sufficiently small fixed \(\delta>0\). The argument
works at any balanced global maximizer \(A\); a tie is needed
only for (LR5) and the comparison with \(T\).

To prove (LR6), put \(e_i=1_{P_i}-1_{A_i}\), and let \(L\) be
the integral of the adjusted-score loss from changing \(A\)'s
label to \(P\)'s at the base correlation \(r\). The score prices
cancel in this integral by exact balance. We have
\[
 S_r(A)-S_r(P)=2L-\sum_i\langle e_i,T_re_i\rangle .
                                                               \tag{LR7}
\]
Transfers between the two labels of a tube have score loss
at least \(c|u|\). Transfers to the third label have score loss
at least \(c>0\) in a sufficiently thin tube. These bounds follow
from the positive normal derivative and the strict third-score
separation on the compact subarcs. The Gaussian Jacobian in normal
coordinates is a smooth positive factor times \(d\nu(s)\,du\),
with the factor equal to one at \(u=0\).

Here is the localized normal-flux limit used below. Suppose
\(L_j=O(\eta_j^2)\) for \(\eta_j\downarrow0\). For each one-sided
binary transfer, integrating its density in \(u\) and dividing
by \(\eta_j\) gives bounded \(L^2(d\nu)\) densities on the
carrying arc. This follows from the one-dimensional bound
\(\int_0^\delta u\,p(u)\,du\ge(\int_0^\delta p(u)\,du)^2/(2\sup p)\)
for nonnegative subdensities \(p\). In particular
\(\sum_i\|e_{j,i}\|_1=O(\eta_j)\). The normalized mass of
third-label transfers tends to zero. For every fixed positive
normal distance from the interface, the normalized binary mass
past that distance also tends to zero by the score margin.
Passing to weak limits of the one-sided densities, the signed
measure limit of \(e_{j,i}\gamma_2/\eta_j\) is
\(\lambda_i^f=\sum_{e=(k,l)}(e_k-e_l)_i f_e\nu_e\) with
\(f\in\mathcal F_I\). The sign is outward speed for label \(k\).
Exact balance passes to the limit because the supports are compact.

Using the normal score expansion \(a(s)|u|+O(u^2)\) and the
Jacobian factor \(1+O(|u|)\), first restrict the mass-cost
estimate to \(|u|\le\delta'\) and then let \(\delta'\) decrease
to zero. Weak lower semicontinuity of the square norms gives
\[
 \liminf_j 2L_j/\eta_j^2\ \ge\ \|f\|_a^2.
\]
Indeed the two one-sided squares bound the square of their
antisymmetric difference. The truncated weak limits are the
same as the full ones since the normalized mass omitted for
each fixed \(\delta'>0\) tends to zero. On the compact product
of the tubes, continuous kernels can be approximated uniformly
by finite sums of separated continuous functions. Bounded
total variation and weak convergence therefore give convergence
of the interaction in (LR7) divided by \(\eta_j^2\) to
\(\langle f,Kf\rangle\).
[Planar normal-flux method (B13)–(B14)](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison).

These estimates imply \(S_r(A)-S_r(P)\ge c L\) for all \(P\)
in sufficiently thin tubes. In fact a sequence contradicting
this assertion with \(L_j\downarrow0\) could be normalized by
\(\eta_j=\sqrt{L_j}\). Its limiting speed satisfies
\(\|f\|_a^2\le2\), and (LR2) bounds its limiting interaction
by \(2(1-c_I)\), after reducing \(c_I\) to lie in \((0,1)\).
Equation (LR7) gives the contradiction for \(c<c_I\).
To choose a single tube width, prove this assertion for all
sufficiently small \(L\) in an initial fixed tube. Then shrink
the tube so every such perturbation has small \(L\), since
the adjusted scores are bounded and the tube mass tends to zero.

At the perturbed correlation the exact quadratic expansion
in \(e\) also gives
\[
\begin{split}
 S_{r+t}(P)-S_{r+t}(A)
   ={}&-2L+\sum_i\langle e_i,T_r e_i\rangle
                      +2t\sum_i\int(\partial_r U_i)e_i\,d\gamma_2\\
      &+O\!\left(|t|\,\|e\|_1^2+t^2\|e\|_1\right),
\end{split}                                                   \tag{LR8}
\]
where \(\|e\|_1=\sum_i\|e_i\|_{L^1(\gamma_2)}\). The kernel
derivatives for the interaction are bounded on the compact
product. The derivatives of each smoothed score are bounded
on the compact observation set, with an integrable Gaussian
bound on the source kernels on a fixed small correlation band.
Thus the constants in the remainder are uniform over ordinary
partitions in the tubes.

For partitions within \(o(t^2)\) of the supremum on the left
of (LR6), the strict base loss and \(\|e\|_1\le C\sqrt L\)
in (LR8) give \(L=O(t^2)\). Normalize now by \(|t|\) and
apply the flux limit. The linear term divided by \(|t|\)
converges to \(G(f)\), so the upper limit in (LR6) is bounded
by \(\sup_{f\in\mathcal F_I}\{2\operatorname{sgn}(t)G(f)-Q_A(f)\}
=\Gamma_I\). The sign is immaterial by replacing \(f\) by \(-f\).

For the lower limit, approximate a maximizing speed in
\(\mathcal F_I\) by smooth speeds supported strictly inside
the selected subarcs. Use compact bumps to keep exact first
flux zero within the chosen tubes. The incidence vectors of
the carrying pairs span exactly the possible fluxes, and one
bump on each edge of a spanning forest of these pairs repairs
all such fluxes. Move the interfaces by the signed normal
strips of amplitude \(t\) given by the approximate speed.
Only these binary transfers are used. The same bumps correct
the remaining \(O(t^2)\) mass errors exactly, since each
connected component of the carrying-pair graph conserves
its total volume. Formula (LR8) and the one-sided normal strip
integrals give the quadratic coefficient \(2G(f)-Q_A(f)\)
for a fixed approximate speed. Passing to the weighted limit
gives (LR6).

Combining (LR5)–(LR6) with the fixed-shape Taylor expansion
gives a strict deficit
\[
 S_{r+t}(T)-\sup_{P\in\mathcal C_\delta} S_{r+t}(P)
          = \bigl(\tfrac12D''(r)-\Gamma_I\bigr)t^2+o(t^2).
                                                               \tag{LR9}
\]

## A quantitative escape estimate

The outside changed-label mass of any competitor beating the
tripod at \(r+t\) is at least \(c|t|^{1+r}\) for small nonzero
\(t\). Borell's half-space bound strengthens this to
\[
 \gamma_2\big(\{P\ne A\}\setminus K_\delta\big)
       \ \ge\
     c'\,|t|^{1+r}\big(\log(e/|t|)\big)^{r/2}.               \tag{LR10}
\]
Both estimates use the same tubes as (LR9), reduced in width
if necessary, and require no condition on which labels
occur in the tubes.

First extend the compact growth estimate to transfers in
the tubes that may change the masses by \(h\). If \(L_{\rm in}\)
is their base score regret and \(e_{\rm in}\) their indicator
change, then
\[
 \sum_i\langle (e_{\rm in})_i,T_r(e_{\rm in})_i\rangle-2L_{\rm in}
                    \le -c_0 L_{\rm in}+C|h|^2.              \tag{LR11}
\]
The normal mass-cost bound gives
\(\|e_{\rm in}\|_1\le C\sqrt{L_{\rm in}}\), and the local
Mehler interaction is at most \(C L_{\rm in}\). If (LR11)
failed for every \(C\) and arbitrarily small regrets, one
could take a sequence with regrets decreasing to zero and
\(|h|/\sqrt{L_{\rm in}}\to0\). The normalized flux limit is
then balanced. Its interaction divided by \(L_{\rm in}\) has
upper limit at most \(2(1-c_I)\) by (LR2), as in the proof of
the balanced growth bound. This contradicts failure of
(LR11) for \(0<c_0<2c_I\). Shrinking the tubes ensures all
interior regrets are in the range of (LR11).

For a balanced competitor \(P\), split its change into
\(e=e_{\rm in}+e_{\rm out}\) and put
\(\delta_{\rm out}=\gamma_2(\{P\ne A\}\setminus K_\delta)\).
The interior flux \(h=\int e_{\rm in}\,d\gamma_2\) is the
negative of the outside flux, so \(|h|\le C\delta_{\rm out}\).
The price terms cancel before the split. In the exact
stability expansion of \(P\) about \(A\), the outside base
score regret is nonnegative and can be dropped in an upper
bound. Uniformly for \(\rho=r+t\) near \(r\), the cross
interaction of \(e_{\rm in},e_{\rm out}\) is at most
\(C\|e_{\rm in}\|_1\delta_{\rm out}\) in absolute value.
Indeed at a point \(x\) in the compact tubes, completing
the square in the source \(y\) bounds
\(M_\rho(x,y)\) by \((1-\rho^2)^{-1}e^{|x|^2/2}\).
The covariance derivative of the scores of the fixed
partition \(A\) satisfies
\(\left|\partial_\rho T_\rho1_{A_i}(x)\right|\le C(1+|x|)\)
on a compact positive correlation band, by differentiating
the conditional Gaussian density. Therefore the outside
linear term from changing the correlation is
\(O(|t|\delta_{\rm out}\sqrt{\log(e/\delta_{\rm out})})\).
To see this, split the integral over the outside changed
set at radius \(\sqrt{2\log(1/\delta_{\rm out})}\) and bound
the radial Gaussian tail directly. Inside the tubes the
same Taylor estimates as (LR8) apply, with
\(O(|t|\|e_{\rm in}\|_1^2+t^2\|e_{\rm in}\|_1)\) remainder.

For the outside self-interaction, Gaussian hypercontractivity
gives
\[
 \sum_i\|T_{\sqrt\rho}(e_{\rm out})_i\|_2^2
                            \le C \delta_{\rm out}^{2/(1+\rho)}.
                                                               \tag{LR12}
\]
Here \(T_{\sqrt\rho}\) contracts \(L^{1+\rho}(\gamma_2)\) to
\(L^2(\gamma_2)\), and the functions are bounded by one and
supported on the outside changed set.
[Nelson's theorem as stated by Gross, Theorem 1.1](https://doi.org/10.1007/BF02392573).

For the sharper bound, dominate the absolute values by the
indicator of that set and apply Borell's two-set theorem with
equal masses. At \(\delta=\bar\Phi(a)\) the relevant half-space
probability is at most \(C a^{-2}\exp(-a^2/(1+\rho))\).
In fact writing the two coordinates as \(a+u,a+v\), \(u,v>0\),
the Gaussian exponent is the sum of
\(a^2/(1+\rho)\), \(a(u+v)/(1+\rho)\), and a nonnegative
quadratic form in \(u,v\). Dropping the last term and
integrating proves the bound uniformly on compact positive
correlation bands. Using the one-dimensional Gaussian tail
bounds at \(a\) shows the right side of (LR12) can be replaced
by
\[
 C\delta_{\rm out}^{2/(1+\rho)}
                  \big(\log(e/\delta_{\rm out})\big)^{-\rho/(1+\rho)}.
                                                               \tag{LR13}
\]
[Borell's comparison as proved in Mossel–Neeman, Theorems 1.1–1.2](https://ems.press/content/serial-article-files/32017).

If (LR10) failed, there would be a sequence \(t_j\to0\)
of improving competitors with
\(\delta_{\rm out}=o(|t_j|^{1+r}(\log(e/|t_j|))^{r/2})\).
In particular \(\delta_{\rm out}=o(|t_j|)\). The bound (LR13)
is then \(o(t_j^2)\): raising the asserted scale to
\(2/(1+\rho)\), with \(\rho=r+t_j\), gives
\(|t_j|^{2+O(t_j)}(\log(e/|t_j|))^{r/(1+\rho)}\).
The logarithmic factor in (LR13) cancels this factor up to
a bounded ratio, since \(\log(e/\delta_{\rm out})\) is bounded
below by a constant times \(\log(e/|t_j|)\).
Likewise (LR12) alone gives \(o(t_j^2)\) at
\(\delta_{\rm out}=o(|t_j|^{1+r})\).

The exterior linear error is \(o(t_j^2)\) in either case,
since \(u\sqrt{\log(e/u)}\) is increasing for small \(u>0\)
and either asserted scale times its square-root logarithmic
factor is \(o(|t_j|)\).
Improvement gives \(S_{r+t_j}(P)-S_{r+t_j}(A)\ge
D(r+t_j)\), and \(D(r+t_j)=O(t_j^2)\).
Using (LR11) and \(\|e_{\rm in}\|_1\le C\sqrt{L_{\rm in}}\)
in the split expansion yields \(L_{\rm in}=O(t_j^2)\).
In fact the upper bound has the form
\(-c_0L_{\rm in}+C(|t_j|+\delta_{\rm out})\sqrt{L_{\rm in}}
+C|t_j|L_{\rm in}+O(t_j^2)\). The cross interaction is now
\(o(t_j^2)\) and the compact weighted remainders remain
\(o(t_j^2)\). Normalize the interior changes by \(|t_j|\).
Their limiting flux is balanced since
\(|h|/|t_j|\to0\). The same compact flux limit as in (LR6)
then bounds
\(\limsup_j [S_{r+t_j}(P)-S_{r+t_j}(A)]/t_j^2\) by
\(\Gamma_I\). This contradicts improvement and (LR5),
proving (LR10) and the hypercontractive bound.

If all changed points are confined to one fixed compact
set, then the exterior self-interaction is at most
\(C\delta_{\rm out}^2\) by the bounded kernel on that compact
product. In this case the same argument works under the
weaker smallness condition \(\delta_{\rm out}=o(|t|)\).
The exterior covariance derivative is bounded on the compact
set, so its linear contribution is \(O(|t|\delta_{\rm out})\).
The interior regret is again \(O(t^2)\) and the normalized
interior flux is balanced. Thus a compactly supported improving
family must change at least \(c|t|\) Gaussian mass outside
the selected tubes, with \(c\) allowed to depend on the fixed
compact set.
