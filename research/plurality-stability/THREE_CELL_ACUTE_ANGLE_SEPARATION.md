# Separation from the acute angle threshold

Fix \(0<r<1\). Under the exact-three bounded wedge-end
hypotheses in the [acute-end projection analysis](THREE_CELL_ACUTE_END_PROJECTION.md),
suppose that the three directions of an actual maximizer lie
in an open semicircle. Write \(\theta_0\) for the smaller
successive gap, and set \(x_*=(\sqrt{1+4r^2}-1)/2\).
There is an \(\epsilon(r)>0\) such that
\[
                  \tan^2\theta_0\le x_*-\epsilon(r).
\]
The proof below excludes a whole interval immediately below
the limiting bound, without uniform geometric cutoffs or
uniform constants in the estimates for individual partitions.
The range at or above \(x_*\) is excluded by
[the angle bound](THREE_CELL_ACUTE_END_PROJECTION.md#an-angle-bound-in-the-open-semicircle-geometry)
and [its equality case](THREE_CELL_ACUTE_ANGLE_BOUNDARY.md).

## A uniform gap in the exponents

Use the notation \(R=r^2\), \(\sigma^2=1-R\),
\(\Lambda=R/(2\sigma^2)\), and
\(a(\lambda)=R\lambda/(1+2\sigma^2\lambda)\).
Let \(\lambda_*=x_* /(2\sigma^2)\) and
\(\kappa_*=a(\lambda_*)\), so
\(\lambda_*+\kappa_*=\Lambda\) and
\(2\lambda_*>\Lambda\).
Write \(\lambda_0=\tan^2\theta_0/(2\sigma^2)\) for the
current network. Let \(\theta_*=\arctan\sqrt{x_*}\).
We take \(\lambda_0\) sufficiently close to \(\lambda_*\)
from below.

Fix a common baseline exponent \(\lambda<\lambda_*\)
sufficiently close to \(\lambda_*\), and take the
current \(\lambda_0>\lambda\).
The previous bootstrap applies with this common
height rate on all arms. For small enough loss
\(\varepsilon>0\), the residual rate in (AU3),
\[
   d_\lambda=\min\{\mu(\lambda),
        \lambda+a(\lambda_0)-\varepsilon,
        \lambda_0+\lambda-\varepsilon\},
\]
is as close below \(\Lambda\) as desired by choosing
the baseline and then the neighborhood. The first
and third quantities are strictly above \(\Lambda\)
at the limiting parameters, so only the second
quantity determines this limit.
The shifted-contour update for \(H=V-J\) in (AU4)
tends to \(\Lambda\) as \(d_\lambda\uparrow\Lambda\).
Indeed \(B=d_\lambda/[R(1+2d_\lambda)]\) tends to
\(1/2\), and the scalar \(c\) in (TF8) tends to
\(1/2\) as well. Hence for a fixed
\(\rho_0\in(\lambda_*,\Lambda)\) the baseline and
neighborhood can be chosen to give rate at least
\(\rho_0\) on the right for every \(H\).
By making \(\rho_0\) smaller inside that interval
if needed, the direct replacement error in (AU2)
also has rate strictly greater than \(\rho_0\),
by \(2\lambda_*>\Lambda\).
[The simultaneous bootstrap](THREE_CELL_ACUTE_END_PROJECTION.md#the-simultaneous-bootstrap)
[The contour update, TF8](THREE_CELL_TAIL_FLATNESS.md#improving-the-gaussian-rate-by-shifting-the-inversion-contour)

Every term of \(J\) also has rate strictly above
\(\rho_0\) on a sufficiently small neighborhood.
Indeed the initial affine CDF rates are at least
\(\lambda_0\); their first inverse rates, when
an inverse is needed, are at least
\(\lambda_0/(R-2\sigma^2\lambda_0)\), which
tends to \(1/(2\sigma^2 x_*)>\Lambda\).
From \(Q=V-D=H+J-D\), we obtain for every arm
\[
    u(t)=-C_{\rm aff}(t)/\psi'(d)+O(e^{-\rho_0t^2}).
                                                               \tag{NG1}
\]
Here the quadratic error when recovering \(u\) from
\(Q\) is faster, since the common baseline can be chosen
with \(2\lambda>\rho_0\).
The exponent \(\rho_0\) is fixed as \(\theta_0\) varies
in the neighborhood. Constants in the big-O estimates,
and the points beyond which they hold, may vary with
the partition.

If the larger successive angle exceeds \(\theta_*\),
consider its associated extreme arm. Its direct rate
is then strictly greater than \(\lambda_*\), and
(NG1) gives its height rate strictly greater than
\(\lambda_*\) as well. The local direct derivative
rate after convolution is greater than \(\kappa_*\).
The squares in (AU3) have rates greater than
\(\Lambda\). The mixed observation term also has rate
greater than \(\Lambda\), by
\(\lambda_*+\kappa_*=\Lambda\). For the direct
replacement in (AU2), use the local greater
angle rate plus the common baseline as close to
\(\lambda_*\) as needed. This gives all parts
of \(F_{\rm err}\) faster than \(\Lambda\) on this
arm. The open-semicircle sign contradiction in the
[weak angle bound](THREE_CELL_ACUTE_END_PROJECTION.md#an-angle-bound-in-the-open-semicircle-geometry)
therefore applies to this arm directly.

## Exact analytic reconstruction just below the threshold

We may now suppose that both gaps are at most
\(\theta_*\) and the smaller one is strictly below.
Their sum is less than \(\pi/2\), since
\(\theta_*<\pi/4\). Both other arms therefore
project forward relative to either extreme arm.
In particular the sufficiently far-left half-plane
of an extreme arm consists entirely of its incident
spectator. Its exact projected contrast \(V(s)\) is
the nonzero constant \(\eta-\beta\) for all sufficiently
negative \(s\), up to null sets.

Choose the extreme arm adjacent to the smaller gap.
Its nearest angle is \(\theta_0\) and its other
relative angle is the sum of the two gaps.
As \(\theta_0\) approaches \(\theta_*\) with the
larger gap in \([\theta_0,\theta_*]\), the other
angle's direct Gaussian rate stays strictly
above \(\Lambda\), by the strict \(2\theta_*\) bound
in the equality analysis.
Write \(d,\eta,\tau,c_1=\psi'(d)\) for this arm's
pressure coordinate, step sign, direct shift, and
binary derivative, as in (AB1). Set
\[
   m=\tan\theta_0,\quad \alpha=r\sin\theta_0,\quad
   c=\tau\cos\theta_0,\qquad
   X(t)=\Phi(-(mt+\tau)/\sigma),\quad
   Y(t)=\Phi(-(\alpha t+c)/\sigma).
\]
The same expansion as (AB1) gives
\[
 F_{\rm err}(t)=
     \frac{r\psi''(d)}{c_1^2}X(t)Y(t)
    +\frac{\eta r\cos\theta_0}{\sigma c_1}
                          X(t)\phi((\alpha t+c)/\sigma)
    +O(e^{-(\Lambda+\delta)t^2}) .                         \tag{NG2}
\]
To see the strict error bound, (NG1) makes the error
in replacing \(u\) by \(-\eta X/c_1\) have rate at
least \(\rho_0>\lambda_*\). Its smoothed error rate
can be taken arbitrarily close to \(a(\rho_0)>\kappa_*\),
since the left leakage has rate \(\Lambda\).
Thus the mixed errors in the two products have
rates above \(\Lambda\) for angles sufficiently
close to critical, by the inequalities
\(\rho_0+\kappa_*>\Lambda\) and
\(\lambda_*+a(\rho_0)>\Lambda\).
The affine step derivative replacement itself gains a
strict positive rate from the remote strip using
its positive common height rate. The remaining
square, mixed squared-source, and direct replacement
errors are strictly faster at the limiting parameters
by \(2\lambda_*>\Lambda\) and
\(\lambda_*+a(2\lambda_*)>\Lambda\). The far arm's
smoothed derivative likewise has a strict rate gain.
All inequalities persist on a small enough
angle neighborhood.
[The exact Taylor expansion, AB1](THREE_CELL_ACUTE_ANGLE_BOUNDARY.md#critical-reconstruction-for-equal-gaps)

Each product in (NG2) has an exact bounded,
real-analytic Mehler inverse with right Gaussian rate
strictly above \(\Lambda\) on a sufficiently small
left neighborhood. For this, put \(q=\sqrt{m^2+\alpha^2}\).
The two-coordinate rotations in the equality proof
express each product in the form
\(P_\sigma B(qt)\), with \(B\) bounded and
supported on a left half-line. In the present case
\(q<r\), so its Mehler inverse is
\[
     s\longmapsto
          P_{\sigma\sqrt{1-q^2/r^2}}B((q/r)s).
\]
Indeed integrating this function under \(N(rt,\sigma^2)\)
adds independent variances summing to \(\sigma^2\).
The inverse is real analytic, and its right Gaussian
rate is at least \(q^2/[2\sigma^2(R-q^2)]\) up to a
linear shift, which tends to infinity as \(q\uparrow r\).
Denote the weighted sum of these two inverses,
with the coefficients in (NG2), by \(B_{\rm err}\).

In the model (AU4) the nearest affine CDF has just
one inverse term and that inverse is supercritical.
The far affine CDF is already supercritical.
Thus \(\mathcal E\) has rate strictly above \(\Lambda\).
For \(H_1=V-J-B_{\rm err}\), (AU3), (AU4), and (NG2)
give a residual \(T_rH_1-rH_1\) with rate strictly
above \(\Lambda\) on the right.
The first Hardy step gives left support of \(H_1\).
Then \(T_rH_1\) itself has the strict supercritical
rate, and the second Hardy step gives \(H_1=0\)
globally. This argument only uses bounded
extensions on the left.
[The two Hardy steps for a supercritical remainder](THREE_CELL_ACUTE_END_PROJECTION.md#an-angle-bound-in-the-open-semicircle-geometry)

Consequently \(V=J+B_{\rm err}\) is real analytic
on the whole line and tends to zero at \(+\infty\).
It agrees on a left interval, up to null sets,
with the nonzero constant \(\eta-\beta\). Continuity
and analytic continuation make it constant globally,
a contradiction. All neighborhood choices above
depend on \(r\) through strict inequalities of
exponents. This gives the claimed positive separation
from the limiting angle bound.
