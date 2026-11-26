# Angular separation of blended factorial symbols through the asymmetric cusp

Let \(a,b\) be distinct fully shared ordinary exponent columns.
Omit equal rows and set \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\),
\(W=\sum_i d_i^2/m_i\), and
\[
 G(q)=\sum_i[(a_i+d_iq)\log(a_i+d_iq)]_{\rm c},\qquad
 F(q)=G(q)-q\log q-(1-q)\log(1-q).
\]
Here \([\,\cdot\,]_{\rm c}\) subtracts the endpoint chord.
Put \(\rho=q(1-q)G''\) and \(\delta=\max\rho-1\).
Use the normalized polynomial at \(m=n,N=n+1\) with
coefficients \(\binom mk e^{m\widehat L_m(k/m)}\), where
\(\widehat L_m=(1-\lambda_m)L_n+\lambda_m L_N\) and the
two blend parameters belong independently to \([0,1]\).

**Theorem.** On each band \(0<W_-\le W\le W_+<\infty\)
there are \(\delta_0,\gamma>0\) and a uniform terminal
threshold such that, if \(|\delta|\le\delta_0\), the two
root argument sets have circular distance at least
\(\gamma/n^2\) above that threshold, for every choice
of the blends. Consequently there is one absolute
eventual terminal threshold for angular separation
of all ordinary fully shared profiles at unrestricted
curvature, for all such blends. The constants are
unevaluated.

## Physical lattice integral and regular complement

Use the centers \(q_*, w_*=F'(q_*)\) with \(F'''(q_*)=0\)
from the
[shared moving quartic argument](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md#uniform-centers-and-finite-integral-reduction).
Put \(a_2=F''(q_*)/2\), \(g_4=F''''(q_*)/24\),
\(\ell_0=(-4g_4)^{-1/3}\), and
\(\ell_4=(-4g_4)^{-1/4}\).
Here \(-g_4\) is bounded above and away from zero,
\(|a_2|\asymp|\delta|\) has its sign, and
\(q_*,1-q_*\) are bounded away from zero.
With \(w=-\log z\) and \(\beta=w_*-w\), the centered
physical phase is
\[
        f_\beta(h)=\beta h+a_2h^2+g_4h^4+h^5E(h),
                       \qquad h=q-q_*,
\]
where the analytic bounds on a fixed disk are uniform.

On a fixed small complex \(\beta\)-disk the cotangent
panels and physical tail estimate of that argument
give the normalized integer polynomial as
\[
       I(m,\beta)+\mathcal E_m(\beta),\qquad
       I(m,\beta)=\int_{-h_0}^{h_0}
               A_m(q_*+h)e^{m f_\beta(h)}\,dh,\qquad
                \mathcal E_m=O(e^{-c n}).                    \tag{1}
\]
The nonvanishing normalization can be chosen real for
real \(w\). Use the affine blend parameter in complex
order joining the two prescribed integer values.
The
[analytic blend estimate](MULTIVARIATE_FIXED_SHARED_BAND_COMMON_SYMBOL_ANGLE_SEPARATION.md#the-same-stationary-geometry-for-both-blends)
gives uniformly bounded nonzero analytic \(A_m\) on
the fixed physical disk and a fixed relative order
disk \(|m-n|\le\eta n\).
The leading phase and fixed real tail loss are unchanged,
so the same cotangent error in (1) applies, analytically
on the physical field disk. At integer order it is real
for real \(\beta\), by reality of the polynomial and the
real-interval integral. The full integral \(I\) is the
complex-order interpolation, including any exterior
pieces in its scaled contours. Only \(\mathcal E_m\)
is discarded at integer order.

Consider root angle pairs at distance at most
\(\gamma/n^2\). By conjugation use both angles in
\([0,\pi]\), so their imaginary \(\beta\)'s are
nonnegative with difference at most \(\gamma/n^2\).
For sufficiently small \(|\delta|\) and large \(n\),
any such pair is in the physical cusp disk of (1).
Indeed the
[blend annulus](MULTIVARIATE_FIXED_SHARED_BAND_COMMON_SYMBOL_ANGLE_SEPARATION.md#a-uniform-annulus-for-blended-factorial-symbols)
confines the fields to a compact strip.
Any counterexample sequence with \(\delta\to0\)
then has a critical retained-and-diffuse limit.
The
[critical regular complement](MULTIVARIATE_CRITICAL_REGULAR_COMPLEMENT.md#the-critical-limiting-graph)
gives a divider, one active balance per interior
height, and the strict regular phase margin away
from the cusp. Its dominant-term and regular
witnesses apply to the blends because the leading
phase is unchanged. The same modulus subtraction
as in the fixed-band blend proof excludes such
a regular pair even with the allowed height
difference. The real upper boundary away from the
cusp has one nondegenerate maximum and a one-term
exclusion. The lower boundary is regular.
This proves the localization by compactness.

## The bounded quartic window

Suppose \(n^{1/2}|a_2|\) and \(n^{3/4}|\beta_m|\)
are bounded for the two orders. Use
\(T_m=-2a_2\ell_4^2\sqrt m\) and \(H=\ell_4m^{3/4}\beta\).
Divide \(I\) by \(\ell_4m^{-1/4}A_m(q_*)\) and put
\(h=\ell_4m^{-1/4}u\).
The resulting integral differs from
\[
                 Z(T_m,H)=\int_{\mathbb R}
                            e^{-u^4/4-T_m u^2/2+Hu}\,du
\]
by \(O(n^{-1/4})\), uniformly on compact \(H\)-disks
and the relative complex-order disk.
This is the
[bounded quartic dominated estimate](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md#bounded-quartic-fields-and-the-adjacent-imaginary-drift).
For complex \(m\), the scaled segment rotates by
\(\arg(m)/4\), preserving uniform quartic decay
for small \(\eta\). The normalized amplitude
difference and quintic remainder have the same
\(O(n^{-1/4})\) quartic-moment bound there.
For the leading entire quartic integrand the
contour rotates back to the real line with
exponentially small tails. These estimates
use the amplitude at the center for normalization.

The relevant model zeros are the simple
imaginary labels \(H=i y_k(T_m)\) and their
conjugates. On the compact real temperature range,
only finitely many labels enter a fixed larger
field disk. Simplicity and compactness give
analytic continuations of those labels in
neighborhoods of the real temperature range.
Decrease \(\eta\) so that each \(T_m\) stays
in those neighborhoods. Rouché gives an analytic
continuum label
\[
                 H_k(m)=i y_k(T_m)+e_k(m),
             \qquad |e_k|=O(n^{-1/4}),\quad
                    |e_k'|=O(n^{-5/4})
\]
on the relative disk and its smaller interior,
respectively. The derivative bound is Cauchy's
estimate in \(m\). The actual integer roots
differ from the continuum label exponentially
little by (1), and compact nonzero model bounds
exclude other labels in the bounded field.
The imaginary physical field on that branch has
real-order derivative
\[
 \ell_4^{-1} m^{-7/4}
      \bigl(-3y_k(T_m)/4+T_m y_k'(T_m)/2+O(n^{-1/4})\bigr).
\]
The
[quartic root-speed inequality](MULTIVARIATE_DIFFUSE_MOVING_CRITICAL_WINDOW.md#the-bounded-moving-critical-field)
gives a uniform strict sign and a positive
lower bound after dividing by \(n^{-7/4}\).
Therefore same-label adjacent imaginary fields
differ by at least \(c_1 n^{-7/4}\).
Distinct labels have gap of order \(n^{-3/4}\).
Either gap excludes the candidate pair.

## Comparison in a regular two-saddle chart

Write \(h=\ell_0 r v\) for a common scale \(r>0\),
\(A=\ell_0 r^4\), \(b=\beta/r^3\), and \(M_m=mA\).
Consider any of the scaled two-saddle charts
below, with both candidate fields in one convex
chart and uniform strict seam, exterior and
quadratic losses. The right-minus-left scaled
stationary action \(D(b)\) satisfies
\(\Re D'(b)\ge c_2>0\), since \(D'\) is the
difference of the continued saddle coordinates.
Choose the analytic split of the continuum
integral there, evaluating the two halves on
their own horizontal saddle paths. Their
exact normalized amplitudes after removal of
\(r M_m^{-1/2}e^{M_m P_s(b)}\) are analytic,
bounded, and bounded away from zero on a fixed
relative complex-order disk when \(M_n\) is large.
Indeed the established quadratic and exterior
losses persist under small relative complex
rotation: the modulus of the phase difference is
bounded by a multiple of its real loss near the
stationary point, on bounded strict-loss
pieces, and on the far quartic tails.
The same Gaussian estimate then applies
with the amplitude in (1).

Let \(U_m\) be the analytic log of the ratio of
the left normalized half to the right normalized
half, with one branch on the common polydisk. Its field derivative is
uniformly bounded on smaller charts, and
\(\partial_m U_m=O(1/n)\) by Cauchy on the
relative disk. The integer remainder in (1)
is exponentially small relative to either
packet: the packet factors cost powers of
\(n\) and \(\exp(O(M_n))\), absorbed by reducing
the physical field and thermal scale caps.
Thus the exact root equations have the form
\[
                e^{M_mD(b_m)}
                         =-\exp(U_m(b_m)+e_m(b_m)),
                    \qquad e_m=O(e^{-c' n}).                \tag{2}
\]
The correction \(e_m\) can be obtained by adding
the physical remainder to one exact packet
before taking the logarithm.

At a root \(|\Re D(b_m)|\le C/M_n\). Put
\(\Delta b=b_N-b_n\). The real slope first gives
\(|\Re\Delta b|\le C/M_n+C|\Im\Delta b|\).
Subtracting the signed modulus equations
improves it to
\[
                   |\Re\Delta b|\le C/(nM_n)
                                      +C|\Im\Delta b| .      \tag{3}
\]
The extra phase term \(A\Re D(b_n)\) and the
order change of \(U_m\) are \(O(1/n)\);
its field change is \(O(|\Delta b|)\), whose
horizontal part is absorbed for large \(M_n\).
If the chart keeps \(|\Im D|\) bounded below,
division in (2) and (3) yield
\[
     |e^{A D(b_n)}-1|\le C/n+C M_n|\Im\Delta b|
                                   \le (C+C'\gamma r)/n .
\]
Here \(|\Im\Delta b|\le\gamma/(n^2r^3)\).
But \(|e^{A D(b_n)}-1|\ge c_3 A\) for small
enough \(A\). Thus large \(M_n=nA\) excludes
the candidate pair.

In the scaled supercritical chart through
the real coexistence line use the real
analytic seam of the thermal proof cited
below. For real \(b\) the exact continuum
halves are positive and their actions and
normalized amplitudes are real.
Take \(U_m\) real there at real orders;
the integer error logs are real there as
well. The complex field bound
\(U_N-U_n=O(1/n)\) then gives imaginary
difference \(O(|\Im b|/n)\) by Cauchy.
The same imaginary bound holds for the
individual integer corrections. Put
\(s=\Im b_n\). In the divided equations
unwrap the phase using small \(A\) and (3).
Reality bounds the horizontal imaginary
changes of \(D\) and of the log amplitudes
by \(C |s|\,|\Re\Delta b|\); the vertical
changes are bounded by \(C|\Im\Delta b|\).
The leading contribution \(A|\Im D(b_n)|\)
is at least \(c_4 A|s|\) by the real
vertical derivative. Hence
\[
                    c_4 A|s|\le C|s|/n+C M_n|\Im\Delta b| .
\]
Absorb the first term for large \(M_n\).
Positivity of the coefficients requires
\(r^3|s|\ge\pi/n\): at smaller arguments
the summands lie in an open half-plane.
The left side is therefore at least a
fixed multiple of \(r/n\), while the last
term on the right is at most \(C'\gamma r/n\).
Choose \(\gamma\) small enough to exclude
the pair.

## Escaping and thermal witnesses

In the escaping window take \(r=|\beta_j|^{1/3}\)
at one of the fields. The scaled phase is
\[
 -v^4/4+bv+(a_2\ell_0/r^2)v^2+r\ell_0^4v^5E(\ell_0rv).
\]
If \(r^2\) is sufficiently large compared to
\(|a_2|\) and \(n r^4\) is large, the
[escaping-quartic witnesses](MULTIVARIATE_EFFECTIVE_SHARED_ESCAPING_QUARTIC_FIELDS.md#effective-finite-contour-near-the-center)
give one-saddle exclusion outside a narrow cone
about the imaginary axis. Both the Gaussian losses
and this exclusion depend only on the indicated
phase and nonzero analytic amplitude bounds.
Thus a root there has imaginary field comparable
to \(r^3\). A second angle within \(\gamma/n^2\)
has comparable imaginary field and meets the
same one-saddle tests if the initial lower
thresholds are enlarged by fixed factors.
Both fields therefore lie in the same cone
using the first field's scale, in a common convex
chart with the strict seam gap and nonzero
imaginary action of the escaping-quartic witness.
The regular scaled comparison applies.

For the bounded thermal window take
\(r^2=2|a_2|\ell_0\).
The scaled phase is a perturbation of
\(-v^4/4\pm v^2/2+bv\) with the sign of \(a_2\).
When \(n r^4\) is large and both \(b\)'s bounded,
the one-field and paired-field charts are
those of the
[supercritical thermal comparison](MULTIVARIATE_EFFECTIVE_SHARED_BOUNDED_THERMAL_CUSP.md#effective-horizontal-contours-and-the-analytic-seam)
or the
[subcritical thermal comparison](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_REGULAR.md#computing-a-cover-of-the-nonfold-field-pairs).
These strict phase and contour tests apply
to (1), including the complex-order extension
of its amplitudes. The same finite cover of
equal-height pairs allows uniformly close
imaginary fields, since its neighborhoods
are open. Indeed \(|\Im\Delta b|\le\gamma/(n^2r^3)\)
tends to zero uniformly when \(n r^4\) is
bounded below. At nonreal regular two-saddle
points the nonzero scaled imaginary action
is a strict witness test. For the supercritical
model it follows from integrating the positive
real saddle separation along the imaginary
axis, starting from zero at \(b=0\).
Cut out a fixed small real-coexistence chart
to apply that test uniformly on the rest.
In the real-coexistence chart use the analytic
seam of the supercritical proof, which is
real for real \(b\), as required above.
In the subcritical regular charts, the nonzero
action is the strict above-fold margin of
the cited witness. This leaves only the
thermal subcritical fold disks.

In those disks use the
[relative complex-mass Airy construction](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_FOLD.md#constructive-airy-chart-and-relative-remainder).
Its scaled finite continuum contour has the
same elementary Airy class and amplitude
bounds, by (1) and the affine complex-order
blend estimate. The entire continuum
exterior remains in the interpolating
contour at complex mass \(M=mA\).
Only the physical integer error is
exponentially small in \(n\), after the
same small-scale shrinkage.
Write \(z_m=M_m^{2/3}\Lambda(b_m)\)
at an integer root, with \(\Lambda\)
the uniform Airy coordinate. Then
\(z_m=a_k+O(M_n^{-1/3})\) near a
negative Airy zero, and
\(\Im\Lambda(b_m)=O(1/M_n)\).
The imaginary derivative of \(\Lambda\)
on horizontal segments is bounded away
from zero. Since
\(|\Im\Delta b|\le C\gamma r/(n M_n)\),
this gives \(|\Delta b|=O(1/M_n)\)
and \(|z_N-z_n|=O(M_n^{-1/3})\).
Shrinking the fixed Airy field disk
forces the same label by the Airy
spacing bound, as in the cited proof.
Cauchy in the relative complex mass
disk on that label gives the same
imaginary branch derivative, of
fixed sign and magnitude at least
\(c_5 |a_k| M_n^{-5/3}\) in \(b\).
The amplitude is only required to be
analytic with a nonzero central value
and uniform bounds in that argument.
Thus the actual adjacent imaginary
displacement in the physical field is
\[
       \ge c_6 r^3 A M_n^{-5/3}
             \asymp r^{1/3}n^{-5/3}\ge c_7 n^{-7/4}.
\]
Here \(A\asymp r^4\) and \(nA\) is
bounded below. The exponential
integer errors cannot remove it.

## Exhaustion and the unrestricted-curvature result

Let \(r_e^3=\max(|\beta_n|,|\beta_N|)\).
Choose fixed large cutoffs for the escaping
test \(n r_e^4\ge B_e\),
\(r_e^2\ge C_e|a_2|\).
If it fails and \(n a_2^2\) is bounded,
then \(n r_e^4\) is bounded as well,
which is the bounded quartic window.
Otherwise take \(n a_2^2\) above a fixed
large cutoff. Failure of the escaping
test then implies \(r_e^2\le C_e|a_2|\),
also when \(n r_e^4<B_e\) by enlarging
the thermal cutoff. With \(r^2=2|a_2|\ell_0\)
this bounds both \(b\)'s and gives the
bounded thermal window. Thus the
constant choices proceed from the
escaping witnesses to the fixed
thermal field bound and its large
mass threshold, and then to the
bounded quartic label ranges.
This exhausts the local field disk.
With its regular complement it
proves the near-cusp theorem.

Combine that theorem with
[fixed shared-band blend separation](MULTIVARIATE_FIXED_SHARED_BAND_COMMON_SYMBOL_ANGLE_SEPARATION.md)
away from \(\max\rho=1\), the
[all-degree weak cap](MULTIVARIATE_WEAK_ASYMMETRIC_COMMON_SYMBOL_ANGLE_SEPARATION.md)
for small \(W\), and the
[absolute-curvature blend theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27_COMMON_SYMBOL.md)
when \(W\ge54\), since \(S\ge W/2\).
This gives one eventual threshold for
all fully shared ordinary profiles
with unrestricted curvature and both
independent blend parameters.
