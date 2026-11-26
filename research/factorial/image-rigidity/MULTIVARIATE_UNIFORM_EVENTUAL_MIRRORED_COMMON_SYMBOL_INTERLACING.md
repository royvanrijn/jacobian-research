# Uniform eventual interlacing of fully shared mirrored profiles

Let the two distinct exponent columns be fully shared and
mirrored, omitting equal rows. Write \(\ell_i\ge1\) for the
smaller exponents and \(d_i\) for the positive integer gaps.
Use \(W=\sum_i d_i^2/(\ell_i+d_i/2)\) over both orientations,
and let \(L_j\) be the endpoint-normalized log-factorial symbol
per degree at order \(j\).

**Theorem.** For each fixed \(W_+<\infty\) there is one
terminal threshold, uniform over profiles with \(W\le W_+\)
and real smaller exponents, for strict circle-cut angular
interlacing of
\[
 \sum_{k=0}^n\binom nk
      e^{n[(1-\lambda)L_n(k/n)+\lambda L_{n+1}(k/n)]} z^k,
 \qquad
 \sum_{k=0}^{n+1}\binom{n+1}k e^{(n+1)L_{n+1}(k/(n+1))}z^k,
 \qquad 0\le\lambda\le1.
\]
The roots are simple and on the unit circle. For ordinary
integer exponents there is one absolute terminal threshold for
the same assertion, including every \(0\le\lambda\le1\), over
all fully shared mirrored profiles with no bound on \(W\). The thresholds
are unevaluated.

The fixed [subcritical](MULTIVARIATE_SUBCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
and [supercritical](MULTIVARIATE_SUPERCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
bands and the [moving critical window](MULTIVARIATE_CRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
leave a range in which \(W\) approaches four but
\(\sqrt n\,|W-4|\) grows. The missing comparison is resolved
at effective order \(n|W-4|^2\).

## The rescaling in the intermediate thermal range

Put \(\epsilon=|W-4|>0\), \(\sigma=\operatorname{sgn}(W-4)\),
\(\nu=m\epsilon^2\), where \(m\) is an order near \(n\), and let
\[
             \theta=\epsilon^{3/2}s,\qquad
             x=\sqrt{\epsilon}\,u,\qquad
 F_\epsilon(u;s)=
       \epsilon^{-2}\bigl(H(\sqrt\epsilon\,u)-H(0)\bigr)+isu .
\]
Here \(H\) is the centered entropy phase from the moving
critical window. We take \(\epsilon\le\epsilon_0\) small and
\(\nu_n=n\epsilon^2\ge B\) large, with fixed constants.
On bounded \(u\)-disks, \(F_\epsilon\) and its derivatives
are an \(O(\epsilon)\) perturbation of
\[
                   F^0(u;s)=\sigma u^2/2-\kappa u^4+isu,
          \qquad 11/12\le\kappa\le4/3 .
\]
For \(\theta\) on a fixed small complex neighborhood of
zero the exact lattice-to-integral estimate of
[the critical proof](MULTIVARIATE_CRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md#small-angles-outside-the-bounded-critical-field)
still holds with error \(O(e^{-c n})\) and the same bound
for its angular derivative. In that proof the horizontal
phase bound uses \(\max(H(x)-H(0))=O(\epsilon^2)\) in
place of \(O(K^2/n)\). Choosing \(\epsilon_0\) small
preserves the strict horizontal-kernel loss. The real
tail and vertical-panel losses also persist. As a result,
the normalized integer polynomial is the value of
\[
       I(m,s)=\sqrt{\epsilon}
        \int_{-\delta/\sqrt\epsilon}^{\delta/\sqrt\epsilon}
          a_m(\sqrt\epsilon\,u)e^{\nu F_\epsilon(u;s)}\,du       \tag{1}
\]
plus an analytic-in-angle \(O(e^{-c n})\) error. The
amplitude \(a_m\) is the analytic complex-order Gamma
amplitude of that proof, with the same \(\lambda\) and
upper order fixed for the two integer evaluations.
The entire finite interval integral defines the
interpolation for complex \(m\) in \(|m-n|\le\eta n\).
Thus its exterior pieces are included in the
complex-order interpolation. Only the lattice error
is discarded at the integer orders.

Angles separated by a fixed positive amount from
zero are covered by the regular charts at the critical
family, by the
[regular-angle argument](MULTIVARIATE_CRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md#regular-angles-and-combination-of-labels).
The strict contour, amplitude, and phase tests persist
for every sufficiently small \(\epsilon\), uniformly.
Otherwise one could take a sequence with \(\epsilon\to0\)
and a limiting critical profile and angle in the compact
cover, where all tests are strict.

For \(s\ge S_0\) large and \(\theta\le\theta_0\) small,
the outer quartic proof applies directly: its two
parameters are
\[
            |W-4|/\theta^{2/3}=s^{-2/3},\qquad
            n\theta^{4/3}=\nu_n s^{4/3}.
\]
These meet its small-perturbation and simple-saddle
conditions by choosing \(S_0,B\) large. It remains to
exclude common angles with \(s\in[0,S_0]\).

## Regular saddles on the intermediate scale

For a reflected pair of stationary points \(a+ib,-a+ib\)
with \(a>0\), the quartic model gives
\[
            a^2=3b^2+\sigma/(4\kappa),\qquad
            s=32\kappa b^3+2\sigma b .
\]
On the horizontal through them the real quartic has
its two maxima exactly at \(\pm a\). Its quadratic
coefficient in the horizontal variable is
\(\sigma/2+6\kappa b^2\), and its quartic coefficient
is \(-\kappa\).

For \(\sigma=1\) these saddles exist from \(s=0\) onward
with \(a>0\) and \(b\ge0\). For \(\sigma=-1\), the first
pair occurs at the fold
\[
       b_*=(12\kappa)^{-1/2},\qquad s_*=\tfrac23 b_* .
\]
The pair with \(a>0\) exists for \(s>s_*\). Below this
fold, the relevant saddle is purely imaginary \(ib\)
on the smaller branch \(s=b-4\kappa b^3\), \(0\le b<b_*\).
On its horizontal the real difference from the saddle
is \((-1/2+6\kappa b^2)u_{\rm re}^2-\kappa u_{\rm re}^4\),
strictly negative off zero. The perturbed fold
\(s_\epsilon=s_*+O(\epsilon)\) remains real, by solving
the stationary and curvature equations on the
imaginary axis and using their nondegenerate
Jacobians there.

Away from the fold the appropriate saddle and its
strict real-height loss persist in (1) on compact
scaled ranges. On the bounded portions of the
contour this is by the quartic model and analytic
implicit function theorem. At large \(|u|\), bounded
imaginary height has the loss \(-c|u|^4\) after
shrinking \(\delta\): the higher Taylor terms are
\(O(\epsilon |u|^6)\), so their relative size is
\(O(\delta^2)\) up to the endpoint strip.
The fixed vertical endpoint connectors have strict
unscaled height loss. These losses are also valid
for small relative complex-order variations, by the
bound on the ratio of the modulus of the phase
difference to the real loss used in the critical
quartic proof.

With one saddle below the fold, the single Gaussian
term is nonzero and dominates the contour and lattice
errors for \(\nu_n\ge B\). In the two-saddle ranges split
the deformed horizontal into two reflected halves at
the imaginary axis, splitting the real lattice error
equally between the halves. The right packet before
that exponentially small integer correction is
\[
      \sqrt{\epsilon}\,\nu^{-1/2}C_m(s)
          e^{\nu F_\epsilon(u_+;s)}(1+O(1/\nu_n)).
\]
This is the horizontal Gaussian estimate from the
critical proof, with \(r\) replaced by \(\sqrt\epsilon\),
and the error is analytic on a relative complex
\(m\)-disk. \(C_m,C_m^{-1}\) are uniformly bounded
there. The exponent at the right saddle has imaginary
part \(d_\epsilon(s)\) with \(d_\epsilon'(s)=\Re u_+>0\).
In the model it begins at zero at \(s=0\) on the
supercritical side and at \(s=s_*\) on the above-fold
side. Hence \(d_\epsilon\) is bounded below by a
positive constant on compact two-saddle intervals
away from those endpoints. Cauchy's estimate in
\(m\) gives the adjacent packet phase difference
\[
                        \epsilon^2 d_\epsilon(s)+O(1/n).
                                                               \tag{2}
\]
For \(B\) large this is strictly between zero and
\(\pi\) on these intervals.

On the supercritical side near \(s=0\), use instead
the exact right half of the original real continuum
interval as the analytic packet. The other half is
its conjugate for real \(s,m\). At zero the right
saddle is strictly interior and separated from the
half endpoint by a fixed positive loss in the scaled
phase. The Gaussian chart for the right half thus
holds on a fixed complex disk about \(s=0\).
The contribution from the half endpoint is an
exterior part of this same analytic integral.
The reciprocal packet normalization costs at most
a power of \(n\) times \(\exp(O(\nu_n))\) there.
Shrinking \(\epsilon_0\) absorbs this into the
\(O(e^{-cn})\) lattice error. Its Gaussian factor
and relative error are real at \(s=0\) for real \(m\),
as is the exponentially small integer correction.
Taking a complex \(s\)-derivative bound for the
order differences by Cauchy's estimate improves the
error in (2) to \(O(s/n)\). Here both the coefficient
and remainder differences already have size
\(O(1/n)\) on the fixed complex \(s\)-disk by
the complex-\(m\) bound. Since
\(d_\epsilon(0)=0\) and \(d_\epsilon'(0)\) is uniformly
positive, (2) again lies strictly between zero and
\(\pi\) for positive \(s\) near zero, with large \(B\).
No positive ray is a zero by positivity. In these
two-saddle ranges a common argument would make both
right packets purely imaginary, which the phase
differences exclude.

## A real Airy chart at the intermediate subcritical fold

Consider \(\sigma=-1\) and \(s\) near \(s_\epsilon\).
The exact interpolating integral (1) has the local
Airy form on a fixed small complex \(s\)-disk about
\(s_\epsilon\). To check the contour and uniform
normal form, at the quartic fold move the horizontal
contour to \(ib_*\). With displacement \(v\) about
this point the phase difference is
\(-i c_* v^3-\kappa v^4\), \(c_*=4\kappa b_*>0\).
The right ray may be tilted to argument \(-\pi/12\)
and the left to \(-11\pi/12\) on a fixed small
radius, with joining arcs back to the horizontal.
The cubic real part is nonpositive on the swept
wedges and strictly negative on the selected rays.
The quartic makes the loss strict even at the
horizontal ends of those arcs. This persists under
the uniform perturbation at the actual fold:
the phase has vanishing first and second derivatives
there and a negative imaginary cubic coefficient.
The \(-\kappa v^4\) coefficient perturbs by
\(O(\epsilon)\) and is real. Higher Taylor terms
are uniformly smaller on a small fixed circle.
Beyond this local contour the horizontal and
endpoint loss estimates used above hold at the
fold, and persist in a small complex field disk.
In the far scaled horizontal tails the quartic
loss bounds the absolute integral uniformly.
Small relative complex variation in \(m\) preserves
the exterior losses.

The analytic cubic normal form is
\[
           F_\epsilon(u;s)=C(s)+w^3/3-\chi(s)w .
\]
At the fold choose \(w=i(3c_\epsilon)^{1/3}v+O(v^2)\)
with \(c_\epsilon=c_*+O(\epsilon)>0\). This sends
the indicated left and right rays into the two
elementary Airy decay sectors in their standard
orientation. The normal form nearby is given by
the
[analytic simple-fold construction](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_FOLD.md#testing-a-fold-chart),
a local proof of the
[Chester–Friedman–Ursell cubic transformation](https://doi.org/10.1017/S0305004100032655).
All disk sizes and derivative bounds are uniform,
since the quartic data range over a compact set and
the perturbations are uniformly small. Reflection
\(u\mapsto-\bar u\) becomes conjugation of \(w\).
In particular \(\chi(s_\epsilon)=0\), \(\chi\) is
real for real \(s\), and
\(\chi'(s_\epsilon)=-(3c_\epsilon)^{-1/3}\).
This derivative follows by comparing the terms
linear in displacement and in \(s-s_\epsilon\)
in the normal form, and persists with a uniform
negative bound on the smaller real chart.

With \(z=\nu^{2/3}\chi(s)\), divide the integral by
its nonzero factor \(\sqrt\epsilon\,e^{\nu C(s)}
\nu^{-1/3}\) and the elementary Airy scalar. The
result is
\[
 P(\nu,s)\operatorname{Ai}(z)
   +\nu^{-1/3}Q(\nu,s)\operatorname{Ai}'(z)+E(\nu,s),
 \qquad
 |E|\le C\nu_n^{-1}
 \left(|\operatorname{Ai}(z)|
          +(1+|z|)^{-1/2}|\operatorname{Ai}'(z)|\right).         \tag{3}
\]
This is the
[constant and linear Airy amplitude division](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md#the-subcritical-fold-and-complex-order)
on the fixed rescaled contour, with \(P,P^{-1},Q\) uniformly
bounded on the complex-order disk. The Stirling amplitudes
have the needed bounds also in the effective order variable,
because \(m=\nu/\epsilon^2\) and the relative disks are
unchanged. The amplitude at the central point is nonzero.
Here (3) describes the full interpolation (1).
Its exterior terms are analytic in the same relative
complex order and, after shrinking \(|\chi|\), smaller
than \(e^{-c'\nu_n}\) times the
[joint Airy norm](AIRY_EFFECTIVE_ELEMENTARY_JOINT_NORM.md).
Indeed the strict-loss tails have rate \(c'\nu_n\)
before division by that norm, whose inverse is bounded
by \(C\exp(|z|\sqrt{1+|z|})\). The subtraction in the
amplitude division uses the same fixed cubic rays
and the same strict relative loss at the finite
radial exits. The lattice error, which is only
available at integer orders, is exponentially small
in \(n\) relative to the normalized joint norm:
the factors from \(\sqrt\epsilon,\nu_n\) cost at most
powers of \(n\), since \(\nu_n\ge B\), and
\(\exp(O(\nu_n))\le \exp(O(\epsilon_0^2 n))\) is
absorbed by decreasing \(\epsilon_0\).
This distinction retains the full effective-order
error in the analytic branch comparison.

At a true integer-order zero in this fold chart the
Airy equation and (3) give
\(|z-a_k|\le C\nu_n^{-1/3}\) for a negative Airy
zero \(a_k\), as in the
[fold label estimate](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md#the-subcritical-fold-and-complex-order).
At a proposed common angle, the same \(s\) is used
at both orders, and the change of \(\nu^{2/3}\chi(s)\)
has size \(O(\epsilon^2\nu_n^{-1/3}|\chi(s)|)\).
Thus the Airy spacing bound \(c(1+|a_k|)^{-1/2}\)
forces the same label after a fixed shrinkage of
the \(\chi\)-disk. The \(C\nu_n^{-1/3}\) bound is
uniform before this shrinkage, and
\((1+|a_k|)\nu_n^{-2/3}\) is then uniformly small
for large \(B\).

For this label apply Rouché in the scaled \(z\)-disk of
radius a suitable fixed multiple of \(\nu_n^{-1/3}\)
about \(a_k\) to the continuum interpolation (3),
substituting \(s=\chi^{-1}(\nu^{-2/3}z)\).
The affine Airy approximation and remainder estimates
of the same label proof are uniform on a fixed
smaller relative complex-order disk. The unique
zero is \(a_k+b_k(\nu)\) with
\[
            b_k=O(\nu_n^{-1/3}),\qquad
            \partial_\nu b_k=O(\nu_n^{-4/3})
\]
by Cauchy's estimate. For real \(\nu\) this zero
is real by reflection of (1) and uniqueness.
Its scaled physical angle \(s(\nu)\) has derivative
\[
   \frac{d s}{d\nu}
        =(\chi^{-1})'\left(
                   -\tfrac23 a_k\nu^{-5/3}
                           +O(\nu_n^{-2})\right).
\]
The derivative is negative and has magnitude at least
\(c|a_k|\nu_n^{-5/3}\) for large \(B\). The two
integer orders change \(\nu\) by \(\epsilon^2\).
The resulting displacement in \(s\) is
at least \(c\epsilon^2|a_k|\nu_n^{-5/3}\), whose
size is bounded below by a multiple of \(n^{-5/3}\)
because \(\epsilon^2=\nu_n/n\), \(\nu_n\le n\),
and \(|a_k|\) is bounded away from zero.
The true integer zeros differ from these continuum
roots by exponentially small errors in \(s\):
on the boundary of the isolating disk the affine
model bounds \(I\), after division by the Airy
normalization and \(\operatorname{Ai}'(a_k)\), below
by a multiple of the disk radius. Divide by the
linear factor of its sole zero. The reciprocal of
this analytic quotient is bounded inside by the
maximum principle. The normalized lattice error
is \(O(e^{-c n})\) there, since the Airy joint norm
divided by \(|\operatorname{Ai}'(a_k)|\) is bounded
on the disk by Airy's equation. Thus the continuum displacement
excludes a common angle in the Airy chart.

## Uniform combination and order of the labels

The regular intermediate ranges and the subcritical
Airy chart give disjoint arguments whenever
\(0<|W-4|\le\epsilon_0\), \(n|W-4|^2\ge B\), for
all sufficiently large \(n\).
The remaining part of this neighborhood is the
moving critical theorem with \(K=\sqrt B\).
Outside it the fixed subcritical and supercritical
bands apply. In the compact subcritical band the
same-grid freezing bound multiplied by \(\lambda\)
in the jump measure makes all the intermediate
lower-grid movements smaller than the robust actual
angular gap. At small \(W\) the
[uniform weak theorem](MULTIVARIATE_UNIFORM_WEAK_COMMON_SYMBOL_INTERLACING.md)
gives the corresponding statement. For bounded \(W\)
the same-grid theorem keeps the lower-grid roots
simple on the circle throughout. These facts give
the required uniform disjointness over each bounded
upper range of \(W\).

To order the labels at any such profile, increase
all smaller exponents by a common real factor while
keeping the integer gaps and the mirror pairing
fixed. This is a continuous path on which \(W\)
decreases to zero. Stop in the weak band, where
the actual and frozen pairs have the known
interlacing order and the same-grid bound preserves
it at every intermediate \(\lambda\). The uniform
angle exclusion applies throughout the path and
the roots never meet the positive ray. Thus the
interlacing order holds at the original profile.

## Unrestricted curvature of ordinary integer profiles

For ordinary integer exponents,
\(S=\sum_i d_i^2/(\ell_i+d_i)\ge W/2\). When
\(W\ge54\), the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
applies at \(S\ge27\). Write
\(\mathcal L_\lambda=(1-\lambda)L_n+\lambda L_N\).
The factorial symbol satisfies, for \(r\ge2\), with
\(\delta_i=b_i-a_i\) and \(x_i(q)=a_i+q\delta_i\),
\[
  (-1)^r L_\nu^{(r)}(q)
       =\sum_i\delta_i^r\int_0^\infty
         \frac{y^{r-1}e^{-x_i(q)y}}{\nu(e^{y/\nu}-1)}\,dy .
\]
The positive kernel is at most \(1/y\) and its
\(\nu\)-derivative lies in \([0,1/(2\nu^2)]\).
To see the derivative bound, put \(g(t)=t/(e^t-1)\).
The derivative in \(\nu\) is \(-g'(y/\nu)/\nu^2\).
The numerator of \(g'\) is nonpositive, since it is
zero at zero and has derivative \(-t e^t\).
Also \(2g'(t)+1\ge0\), since after multiplying by
\((e^t-1)^2\) it equals \(2e^t(\sinh t-t)\).
Thus \(L_N''\ge L_n''\). The two symbols vanish at
both endpoints, so \(L_N\le L_n\) on \([0,1]\).
The interpolated lower coefficients are positive
and are bounded coefficientwise by the actual ones.

Let \(H_m\) be the left half with half the middle
coefficient when \(m\) is even, using the interpolated
lower coefficients at \(m=n\) and the actual coefficients
at \(m=N\). Reflection gives the exact decomposition
\(P_m(z)=H_m(z)+z^mH_m(1/z)\).

In the coefficient regions of the absolute curvature proof the bounds
\(|H_m-1|<37/200\), \(|zH_m'|<1/5\) on the circle
hold also for the symmetric half, since the estimate
uses the absolute coefficient sums. Thus
\(|\partial_\theta\log H_m|<40/163<1/4\).
In the
[large slope region](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate),
the geometric bound \(c_{m,k}\le8^{-k}\) for the
left half gives
\(|H_m-1|\le1/7\), \(|zH_m'|\le8/49\), and
\(|\partial_\theta\log H_m|\le4/21<1/4\).

In the remaining cluster region the
[paired radius proof](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#paired-radius-and-phase-comparison)
gives \(\|\log H_m^{\rm left}\|\le3(D+1)\) and
\(\|\log H_N^{\rm left}-\log H_n^{\rm left}\|\le1/55\)
on \(|z|\le6/5\) for the actual original left halves,
where \(D=N e^{-S/2}\). The interpolated lower half
also satisfies the individual log bound and differs
from the actual lower half in log norm by at most
\(1/55\). Here are the checks for applying the same
endpoint-log estimates to this change of symbol.

Use \(K,\alpha\) of the cluster proof, write
\(V=\sum_i d_i^2/\ell_i\) and \(\mathcal W=1+V\)
for that proof's parameters, and let
\(L=\max(1,\max_i|\delta_i|/a_i)\) so \(\alpha=2L/n\).
For \(x\in[0,K]\) the row argument \(x_i(x/n)\)
is at least \(a_i/2\), and \(n-x\ge n/2\).
With \(B_\lambda=c_{n,1}^{\lambda}\), the interpolating
function for \(\log(c_{n,k}^{\lambda}k!/B_\lambda^k)\)
at \(x=k\) is
\[
 f_\lambda(x)=\log\Gamma(n+1)-\log\Gamma(n+1-x)
        -x\log n+n\bigl(\mathcal L_\lambda(x/n)
                            -x\mathcal L_\lambda(1/n)\bigr).
\]
For \(2\le r\le K\), the kernel \(1/y\) bound and
\(x_i(x/n)\ge a_i/2\) give
\[
 |\partial_x^r f_\lambda|\le
       H_r:=(r-2)!\alpha^{r-2}4\mathcal W/n .
\]
The binomial term here is bounded by
\((r-2)!(n-x)^{1-r}\) from the same integral
kernel. For the \(\lambda\)-derivative, integrate
the \(\nu\)-derivative bound from \(n\) to \(N\).
Since \(a_i\ge1\),
\(\sum_i|\delta_i|^r/a_i^r\le V L^{r-2}\), giving
\[
 \begin{split}
 |\partial_\lambda\partial_x^r f_\lambda|
 &\le n^{1-r}(r-1)!2^{r-1} V L^{r-2}/n^2\\
 &\le H_r(r-1)/(2n^2)\le (2r/n)H_r .
 \end{split}
\]
Consequently the forward differences and their
exponential edge weights satisfy exactly the
majorants \(M_r=2H_r\) and \((2r/n)M_r\) used for the
ordinary order derivative in the
[cluster log proof](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum).
Also \(0<B_\lambda\le B_0\le101D/100\).
On the full unit interval in \(q\), the curvature
difference is bounded by \(V/(2n^2)\), from
\(x_i(q)\ge\ell_i\ge1\). Its endpoint chord
and evaluation at \(q=1/n\) give
\(|\partial_\lambda\log B_\lambda|\le V/(4n^2)<1/n\)
by the cluster bounds. Thus the identical
hypergraph and logarithmic series estimates give
\(\|A_{n,\lambda}-A_{n,0}\|\le13D/n\) for the
truncated logs on \(|z|\le6/5\).
The truncation and recovery errors are each
less than \(1/120\), since the upper bounds for
\(B_\lambda\) and the polynomial tail only improve
under coefficient domination. This proves
\(\|\log H_n^{\lambda,{\rm left}}-
         \log H_n^{0,{\rm left}}\|\le13D/n+1/60<1/55\)
and the individual bound \(3(D+1)\).

Making each half symmetric changes it by at
most \((3/20)^{m/2}\) on the disk, since
the middle coefficient is bounded by \(8^{-m/2}\).
The relative change is at most
\((3/20)^{m/2}\exp(3(D+1))<e^{-m/2}\), using
\(D/n<1/720000\) and \(n>59000\) in this region.
The log change is at most \(2e^{-m/2}\). Hence the
symmetric halves satisfy
\[
     \|\log H_m\|\le4(D+1),\qquad
     \|\log H_N-\log H_n\|<1/27.
\]
Cauchy from \(6/5\) to the circle gives angular
derivative bounds \(20(D+1)\) for each log and
\(5/27\) for their difference.

In all these regions the circle phase
\[
                    p_m(\theta)=m\theta-2\Im\log H_m(e^{i\theta})
\]
is strictly increasing, starts at zero, and satisfies
\(0<p_N(\theta)-p_n(\theta)<2\theta\) for
\(0<\theta\le\pi\). Indeed the uniform angular
derivative bound on \(p_N-p_n-\theta\) is \(<1\),
using the separate \(1/4\) bounds in the coefficient
regions or the \(5/27\) difference bound in the
cluster region. Reflection gives
\(0<p_N-p_n<2\pi\) also on \((\pi,2\pi)\).
The logarithms are normalized at \(z=0\), so their
circle phases come back to \(2m\pi\) at \(2\pi\).
The zeros occur where \(p_m\) is an odd multiple
of \(\pi\). Therefore the two degree sequences
strictly interlace for \(N\ge3\).
Combining this with the bounded-\(W\) theorem at
\(W_+=54\) gives one eventual threshold for all
ordinary fully shared mirrored profiles, uniformly
in the interpolation parameter.
