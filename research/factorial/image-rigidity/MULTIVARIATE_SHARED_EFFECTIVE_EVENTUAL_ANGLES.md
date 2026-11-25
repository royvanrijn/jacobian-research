# An effective eventual angular cutoff for fully shared moments

Let two distinct nonnegative integer exponent columns have only shared
unequal rows. There is an absolute computable integer \(J_{\rm ang}\)
such that the order-\(N-1\) and order-\(N\) factorial moments have
disjoint root arguments for every \(N\ge J_{\rm ang}\), including after
independent positive rescalings. The bound is uniform in the rows and
dimension, and is unevaluated.

The proof uses the scale hierarchy of the
[effective absolute eventual comparison](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).
The additional comparisons below allow the two fields to have
arbitrary real difference at equal imaginary height.

## Reduction to the near-critical charts

Use \(W,\rho,\delta,w\) from the shared curvature proofs,
with \(\delta=\max\rho-1\).
For \(W\ge128\) one has \(S\ge W/2\ge64\).
Each exclusion in the
[absolute high-curvature proof](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md)
is angular: the small-order and geometric-half cases subtract phase
equations at equal arguments, and the remaining case first compares
the two normalized radii by the modulus equations. That latter
comparison only uses the two actual radii and their common angle.
For a positive computable small \(W_0\), the
[weak angular theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
covers \(W\le W_0\).
On \(W_0\le W\le128\) use the
[fixed-band angular theorem](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md)
outside any prescribed positive critical margin.

The [effective critical entry](MULTIVARIATE_EFFECTIVE_SHARED_CRITICAL_ENTRY.md)
now applies to equal-angle pairs as well. In its exactly critical
limiting family each pair has either a one-term exclusion at some
field or a common regular chart away from the cusp.
The comparison in a common chart is the modulus and phase
comparison from the fixed-band angular theorem, since both real
slope and phase margins are strict there. At the upper boundary
the field outside the prescribed cusp neighborhood has a
one-maximum witness. Thus the same compact cover procedure,
with those angular local witnesses, gives joint entry into
any prescribed physical cusp disk for sufficiently small
\(|\delta|\). Use the real \(q_*,w_*\) and put
\(\beta_j=w_*-w_j\). The imaginary parts of \(\beta_n,\beta_N\)
agree, since the original roots have the same argument.

## Regular comparisons in the scaled quartic fields

Use the two order packets and their common-contour comparison
from the [escaping quartic proof](MULTIVARIATE_EFFECTIVE_SHARED_ESCAPING_QUARTIC_FIELDS.md#linked-comparison-at-a-common-root),
with physical order increment multiplier now
\[
                  e^{f_{\beta_n}(h)+\omega h},
                  \qquad \omega=N(\beta_N-\beta_n)\in\mathbb R .
\]
All stationary and shifted contour estimates there apply with
\(d=|\omega|\epsilon\), where \(\epsilon=|\beta_n|^{1/3}\).
Indeed confinement of both fields in the same cone used their
equal imaginary parts. The first real-height estimate (6)
and real saddle separation (7) bound \(d\) by a fixed constant.
This is the only bound on the drift used in the multiplier errors
(9), including the exterior loss: on the physical connectors
\(|\omega|\le C/\epsilon\), and \(M=n\ell_0\epsilon^4\ge1\)
gives \(\epsilon^{-1}=O(n^{1/4})\).
Here \(M\) and \(\ell_0\) are the linked scales in that proof.
The stationary leading multipliers are consequently bounded above
and away from zero.

Let \(m_R,m_L\) be the upper-to-lower ratios for the two actual
packets on the common lower-order contours. The two exact zero
equations imply \(|m_R-m_L|=O(e^{-c n})\), after adjusting the
positive exponential constant. Taking the real log modulus of
their quotient, and using (6), (7) and (9), gives
\[
 d\le C\bigl(\epsilon^4/M+d/\sqrt M+n^{-2}+e^{-c'n}\bigr).
\]
In fact the leading log difference is
\(D(\beta_n)+\omega(h_R-h_L)\), whose real drift term has magnitude
at least a fixed multiple of \(d\); the real phase difference is
\(O(1/n)=O(\epsilon^4/M)\).
The relative ratio errors in (9) are small so the real logarithms
cost the same errors. Absorb \(d/\sqrt M\) for large \(M\).
Since \(\epsilon^4/M=1/(n\ell_0)\), this gives
\(d=O(\epsilon^4/M)\).
The imaginary leading log difference now has magnitude at least
\(c_0\epsilon^4-O(\epsilon^4/M)\), by the imaginary margin (7).
The logarithms of the relative errors cost \(O(\epsilon^4/M)\)
as well. The whole log difference is small on the chosen physical
field cap, so exponentiation preserves this lower bound up to a
fixed factor. It contradicts the exponentially small exact ratio
difference. This proves the escaping quartic exclusion for
equal-angle pairs with the same adjustable caps.

The same real log comparison applies to the
[subcritical regular thermal proof](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_REGULAR.md#effective-comparison-on-a-common-two-saddle-chart).
There the notation is \(D(b)\), \(M=nc\),
\(\omega=M_N(b_N-b_n)\in\mathbb R\).
The one-saddle versus regular versus fold cover was built on pairs
of fields of equal imaginary part. In a common convex regular
chart, \(|\Re D(b_j)|=O(1/M_j)\), \(\Re D'\) has definite
magnitude, and \(|\Im D|\) has a positive lower bound. Formula (6)
bounds \(|\omega|\), and its proof only needs equal imaginary
parts. Formula (7) holds for this drift on the shared contours,
with leading log difference \(cD(b_n)+\omega D'(b_n)\).
The real log comparison gives
\[
        |\omega|\le O(c/M+|\omega|/\sqrt M+n^{-2}),
\]
with the exponential physical error included in \(O(c/M)\).
Hence \(\omega=O(c/M)\). The imaginary leading difference then
has size at least a fixed multiple of \(c\), larger than the
remaining \(O(c/M)\) ratio error. This proves the regular
exclusions and the localization of any equal-angle pair to
the same prescribed fold disk.
The [subcritical thermal fold proof](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_FOLD.md#airy-localization-and-the-linked-order-increment)
already uses only real drift in its final comparison.
Equal imaginary fields imply a common Airy label and contradict
the nonzero imaginary displacement of that branch over one
mass increment.

## The bounded supercritical thermal fields

Use the common seam and exact upper-to-lower packet ratios from
the [supercritical thermal proof](MULTIVARIATE_EFFECTIVE_SHARED_BOUNDED_THERMAL_CUSP.md#adjacent-multiplier-and-the-two-kinds-of-drift).
With \(A=M_n/n\), \(b=b_n\) and \(T=M_N(b_N-b_n)\in\mathbb R\),
their leading log difference is
\(A d(b)+T d'(b)\).
The first real comparison gives bounded \(|T|\) using only
\(|\Re d(b_j)|=O(1/M_j)\) and equality of imaginary fields.
The ratio and exact cancellation estimates (12) and (13) thus
apply. Take the real log of the quotient of the two actual
ratios. Since \(\Re d'\) is bounded below, these estimates give
\[
                 |T|\le C(A/M_n+|T|/\sqrt{M_n}+n^{-2}).
\]
The true integer-order remainder is absorbed here into \(O(A/M_n)\).
Increase the mass threshold to get \(|T|=O(1/n)\).

For each fixed real \(T\) in that narrowed range, regard the
packet ratios as analytic functions of \(b\) on a smaller strip with
fixed slack. The upper integrand uses the field \(b+T/M_N\) and
the exact lower order seam \(z_c(b)\). The seam is analytic in \(b\),
so these ratios are analytic and nonzero by the packet bounds.
The relative errors in (12) are \(O(1/n)\) on the strip, with
constants independent of \(T\) in the narrowed range.
On the real axis the lower and upper half integrals are positive,
since the seam is real and \(T\) is real. The exact log corrections
are therefore real there. Cauchy's estimate and vertical integration
give an imaginary error \(O(|\Im b|/n)\).
Also \(d'\) is real on the real axis, so
\(|\Im (T d'(b))|=O(|\Im b|/n)\).
The leading term has
\(|\Im(A d(b))|\ge c_d A|\Im b|\).
For large mass this is strictly larger than both errors.
On the small physical cap the imaginary log difference has
modulus less than one, so the two exact multipliers differ by
at least a fixed multiple of \(A|\Im b|\).
The [zero-free interval of imaginary heights](MULTIVARIATE_EFFECTIVE_SHARED_BOUNDED_THERMAL_CUSP.md#adjacent-multiplier-and-the-two-kinds-of-drift)
gives \(|\Im b|\ge c_0/M_n\) at a true zero.
Thus the multiplier difference is at least a positive multiple
of \(1/n\), contradicting (13).

## The bounded quartic analytic roots

Let \(T_j=-2a_2\ell^2\sqrt j\), \(H_j=\ell j^{3/4}\beta_j\)
and \(\lambda_j=j^{-1/4}\) be the variables in the
[bounded quartic expansion](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md#compact-residuals-and-termination).
For prescribed bounded \(T_j,H_j\), the common polynomial truncation
in the small coefficient parameter is
\[
        J(\lambda,T,H)=Z(T,H)+\sum_{p=1}^4 \lambda^p C_p(T,H).
\]
Its model \(Z\) has simple purely imaginary roots at real \(T\).
The coefficients \(C_p\) are common to the two orders for the
fixed profile and have computable analytic bounds on compact
complex neighborhoods of \(T,H\).
At an actual normalized zero one has
\(|J(\lambda_j,T_j,H_j)|\le C\lambda_j^5\).

The two fields of a putative equal-angle pair can be put in
a common model-root chart for all large orders. Indeed \(T_N-T_n\)
and \(\Im(H_N-H_n)\) are \(O(1/n)\), the latter because
\(\Im H_N=(N/n)^{3/4}\Im H_n\). A limiting pair of model roots
therefore has the same real \(T\) and imaginary part.
Both fields are on the imaginary axis, so they coincide.
An effective finite cover of these limiting pairs uses simple
Rouché disks at the model roots and strict nonzero bounds off
those disks. The isolation tests use the bounded quartic integrals
and their Cauchy or differentiated-integral bounds.
At a root choose a convex field disk where \(Z_H\) differs
little from its value at the center. A small complex \(T\)-disk
continues the model root, with a strict Rouché boundary margin.
The same margin continues to \(J\) for small complex
\(\lambda\) by the uniform \(C_p\) bounds, giving an analytic
root \(h(\lambda,T)\). These open witnesses cover the compact
limiting pair set. On a smaller disk, differentiating by Cauchy
gives
\[
 h-h(0,T)=O(\lambda),\qquad
 \partial_T h-\partial_T h(0,T)=O(\lambda),\qquad
                         \partial_\lambda h=O(1).
\]
The uniform derivative lower bound follows by bounding the
variation of \(J_H\) from the nonzero center derivative
on the whole convex field disk.
Integration along the segment from \(H_j\) to \(h\) then gives
\(|H_j-h(\lambda_j,T_j)|=O(\lambda_j^5)\).
All estimates keep the profile-dependent \(C_p\) fixed while
varying \(\lambda,T\), and the chart and threshold choices
are uniform by the computable derivative bounds.

Return to the physical field on this common chart with
\(\Xi(j)=\ell^{-1}j^{-3/4} h(\lambda_j,T_j)\), for real
\(n\le j\le N\). Its derivative is
\[
 \Xi'(j)=\ell^{-1}j^{-7/4}
     \left(\tfrac12 T_j h_T-\tfrac34h-\tfrac14\lambda_j h_\lambda\right).
\]
For an upper imaginary model root \(h(0,T)=iy(T)\),
the [quartic root-speed inequality](MULTIVARIATE_DIFFUSE_MOVING_CRITICAL_WINDOW.md#the-bounded-moving-critical-field)
gives \(3y-2Ty'>0\). Its positive margin on each small chart
is computable from the simple root and analytic derivative.
The lower model root has the opposite sign by conjugation,
and the real model field has no roots.
Thus the imaginary derivative of \(\Xi\) has fixed sign
and modulus at least a positive multiple of \(n^{-7/4}\).
At each integer the difference of the actual physical field
from this surrogate is \(O(j^{-3/4}\lambda_j^5)=O(j^{-2})\).
It cannot cancel the imaginary displacement over one order.
This gives angular exclusion throughout the bounded quartic window.

## Uniform cutoff

Use exactly the order of choices in the
[effective absolute comparison](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md#a-finite-hierarchy-of-cutoffs)
on \(W_0\le W\le128\).
The angular escaping estimates first give a physical field cap
and outer scaled thresholds. Then choose a large bounded thermal
field cap, the two signed thermal estimates and subcritical fold
parameters, and a thermal mass threshold.
The complementary bounded quartic window has fixed bounds on
temperature and field, so its angular estimate above gives a further
profile-independent order. Finally choose the critical entry margin
for that physical field cap, and apply the fixed-band angular theorem
outside the resulting critical neighborhood.
All finite contours, residual bounds, and compact witnesses in
these choices have the same computable tests as their cited
moment comparisons. The added inequalities are strict reductions
of those caps or increases of those thresholds.
Taking the maximum order completes the asserted uniform result
for fully shared profiles.
