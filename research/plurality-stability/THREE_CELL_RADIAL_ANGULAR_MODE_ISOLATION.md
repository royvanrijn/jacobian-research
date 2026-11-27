# Isolating radial obstructions and retaining their angular loss

For every compact correlation interval in \((0,1)\), the adverse scalar
radial directions lie in a uniformly finite-dimensional space. A global
angular upper bound retaining polynomial moments improves the first-mode
bound. At quadratic order about balanced circular masses, some finite
number of these moments outweighs **every** adverse radial direction,
uniformly on the interval. The required order is qualitative here.

At \(\rho=99/100\), two angular modes already give a directed negative
matrix on the explicit radial space spanned by the two functions used
in DG5 and DG9. The older bound has one positive direction in that space.
This provides a concrete mode-isolation and angular-domination calculation.
The finite-amplitude domination needed for the unrestricted global
inequality remains open. All results below have internal checking only.

Use \(d\nu(t)=e^{-t}dt\), \(r=\sqrt{2t}\), and the scalar
\(K_+,J_+\) of [DG15](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md).
The standard-interface kernel in the sources has a different definition.

## RMA1. The scalar adverse space is finite at every positive correlation

Fix \(I=[a,b]\subset(0,1)\). There is an explicit \(j_I>0\) such that
\[
 J_\rho(t)\ge j_I\min(\sqrt t,1),\qquad
 j_I=\frac{e^{-17/8}}{64\sqrt{2\pi}}
             \min\left(\frac{\sqrt2a}{\sqrt{1-a^2}},1\right).
 \tag{RMA1}
\]
In particular, the exponential lower bound used in the earlier scalar
certificate is unnecessary for qualitative compactness at larger \(\rho\).

Here is a direct proof with constants. For the angular weight \(w\) of
(DG15), integration by parts gives
\[
 w(u)=\frac{u\pi}{2I_0(u)}\int_0^{2/3}
       (v-3v^2/2)\sin(\pi v)e^{u\cos(\pi v)}dv
                 \ge\frac{\min(u,1)}{64e}.
\]
For \(u\le1\), restrict to \([1/6,1/3]\), where the polynomial
factor is at least \(1/8\), the sine at least \(1/2\), the exponential
at least one, and \(I_0(u)\le e\). The result is
\(u\pi/(192e)\ge u/(64e)\). For \(u\ge1\), restrict to
\([1/(6\sqrt u),1/(3\sqrt u)]\), use
\(v-3v^2/2\ge v/2\), \(\sin(\pi v)\ge2v\), and
\(\cos(\pi v)\ge1-\pi^2v^2/2\). Also
\[
 I_0(u)\le \frac{e^u\sqrt\pi}{2\sqrt{2u}},
\]
by \(1-\cos\theta\ge2\theta^2/\pi^2\) in its angular integral.
This gives \(7\sqrt{2\pi}e^{-\pi^2/18}/648>1/(64e)\).
The elementary bounds \(3<\pi<4\), \(\pi^2<18\) suffice for the last
inequality.

Conditional on \(|X|=\sqrt{2t}\), write
\(Y=\rho X+\sqrt{1-\rho^2}Z\), with \(X\) in the first coordinate
direction. On \(1\le|Z_2|\le3/2\), the radial variable
\(s=|Y|^2/2\) is at least \((1-\rho^2)/2\). This event has probability
at least \(e^{-9/8}/\sqrt{2\pi}\). The angular parameter is consequently
at least \(\sqrt2\rho\sqrt t/\sqrt{1-\rho^2}\). Apply the weight
bound, drop the nonnegative \(\theta\) term in \(J_+\), and use
\(\min(cx,1)\ge\min(c,1)\min(x,1)\) to obtain (RMA1).

On \(H=L^2(\nu)\), define the **uncentered** operator
\[
 A_\rho(t,s)=\frac{K_+(t,s)}{\sqrt{J_\rho(t)J_\rho(s)}},\quad
 w_\rho=J_\rho^{-1/2},\quad
 B_\rho=P_\rho A_\rho P_\rho,
 \tag{RMA2}
\]
where \(P_\rho\) is orthogonal projection onto \(w_\rho^\perp\).
This \(A_\rho\) has a different additive center from (DG41).
Both operators in (RMA2) are compact and self-adjoint for every
\(0<\rho<1\), and are Hilbert--Schmidt continuous on \(I\).
Indeed, put \(c_I=(1-b)/(1+b)\) and \(C_0=2+2j_1/e\), with
\(j_1\) as in (DG15). The elementary kernel envelope gives
\[
 \|A_\rho\|_{HS}^2\le
 \left[\frac{C_0}{(1-b^2)j_I}
            \left(c_I^{-1}+\sqrt\pi c_I^{-1/2}\right)\right]^2.
 \tag{RMA3}
\]
To check the tail for all correlations, after squaring the kernel and
including both exponential measures its exponent is
\[
 -\frac{(1+\rho^2)(t+s)-4\rho\sqrt{ts}}{1-\rho^2}
                         \le-c_I(t+s).
\]
The reciprocal weights add at most
\(j_I^{-2}(1+t^{-1/2})(1+s^{-1/2})\), which is integrable.
This proves (RMA3), and dominated convergence proves continuity.
It also gives \(w_\rho\in H\) and its continuity in \(H\), hence
continuity of the projections. Bounded mean-zero scalar functions are
dense in the weighted hyperplane, by truncation and subtraction of their
ordinary means, as in (ST3).

For any \(0<\eta<2\), let \(E_\rho\) be the eigenspace of
\(B_\rho\) with eigenvalues greater than \(2-\eta\), inside the
mean-zero hyperplane. It has uniformly finite rank, at most the right
side of (RMA3) divided by \((2-\eta)^2\). For a feasible mass profile
\(z_i=p_i-1/3\), put \(y_i=\sqrt J z_i=a_i+g_i\), its orthogonal
decomposition into \(E_\rho\) and its complement. Then
\[
 \sum_i q_\rho(z_i)\le
 \sum_i\langle a_i,(B_\rho-2)a_i\rangle
                               -\eta\sum_i\|g_i\|_H^2.
 \tag{RMA4}
\]
Thus all possibly positive scalar contribution is isolated in finitely
many radial coefficients. Their summed squared norm is at most
\(2j_{\max}/3\), where \(j_{\max}=1/2+j_1/(2e)\), since
\(\sum_i(p_i-1/3)^2\le2/3\). No numerical count of the full adverse
space at \(99/100\) is asserted by the finite-block checker below.

The nonzero compressed eigenmodes also belong to the natural Hessian
space \(L^2(\sqrt t\,d\nu)\) used below. Their equation is
\(K_+z=\lambda Jz+c\). The map \(z\mapsto J^{-1}K_+z\), from
\(L^2(Jd\nu)\) to that natural space, is Hilbert--Schmidt: its squared
normalized kernel adds \(t^{1/2}/(J(t)^2J(s))\) to the Gaussian
envelope, bounded by a constant times
\(\max(t^{-1/2},t^{1/2})\max(s^{-1/2},1)e^{-c_I(t+s)}\).
Also \(J^{-1}\in L^2(\sqrt t\,d\nu)\). The eigen-equation with
\(\lambda\ne0\) gives the assertion. Thus no unproved domain
extension is needed when applying the angular quadratic gap to these
adverse eigenmodes.

## RMA2. Global angular loss from polynomial moments

For angular parameter \(u>0\), masses \(p,q\in\Delta_3\), set
\(z_0=\min_i\cos(\pi(p_i+q_i))\) and \(D=1+z_0\). Define the
Hermite interpolation polynomial
\[
 H_n(z)=e^{-u}\sum_{k=0}^{n-1}\frac{u^k(1+z)^k}{k!}
              +\frac{e^{uz_0}-e^{-u}\sum_{k=0}^{n-1}(uD)^k/k!}{D^n}
                           (1+z)^n,\qquad n\ge1.
 \tag{RMA5}
\]
At \(D=0\) use its continuous limit, the order-\(n\) Taylor polynomial
of \(e^{uz}\) at \(-1\). Every coefficient in powers of \(1+z\)
is nonnegative. The residual \(e^{uz}-H_n(z)\) is nonpositive for
\(-1\le z\le z_0\) and nonnegative and nondecreasing for \(z\ge z_0\).
For example, with \(v=1+z\), its series is
\[
 e^{-u}\sum_{k=n+1}^\infty\frac{u^k}{k!}
                 v^n(v^{k-n}-D^{k-n}).
\]
This proves all the signs and derivative assertions, including the limit
at \(D=0\), without an interpolation-error hypothesis.

Let \(B_H(p_i,q_i)\) denote the aligned two-arc pairing for the circle
kernel \(H_n(\cos\theta)\), and let
\[
 \mathcal D_{n,u}(p,q)=\sum_i B_H(p_i,q_i)
       -\max_{F,G}\sum_i\langle1_{F_i},H_n(\cos(\cdot))1_{G_i}\rangle,
 \tag{RMA6}
\]
where the maximum ranges over ordinary circular partitions of masses
\(p,q\). This is a joint moment deficit, not three independently
optimized moment vectors. It is nonnegative by circular rearrangement,
since \(H_n\) is nondecreasing in the cosine. Replacing the residual
by zero below \(z_0\) gives a nondecreasing majorant agreeing with it
on every aligned pair's angular differences. Apply
[Friedberg--Luttinger, Theorem I, printed pp. 1--2, (2)--(5)](https://digital.library.unt.edu/ark:/67531/metadc1025197/m2/1/high_res_d/4446751.pdf)
to that nonnegative symmetric decreasing circular majorant, exactly as
in DG9. For sums greater than
one, complement both sets first; the effective sum \(2-p_i-q_i\)
has the same cosine. The resulting global bound is
\[
 S_\rho(A)-S_\rho(T)\le\mathcal E_\rho^{[n]}(p):=
 \iint R_0\left[\sum_i B_u(p_i,q_i)-3B_u(1/3,1/3)
                      -\mathcal D_{n,u}(p,q)/I_0(u)\right]d\nu d\nu.
 \tag{RMA7}
\]
Here \(p=p(\sqrt{2t})\), \(q=p(\sqrt{2s})\) and
\(u=2\rho\sqrt{ts}/(1-\rho^2)\). Measurable attainment is not needed
to apply the pointwise maximum. The maximum value is continuous in the
mass vectors and polynomial coefficients, by compactness and continuous
mass transfer. The total loss lies between zero and the full aligned
arc sum: test a maximizing pair for (RMA6) in the proved remainder bound
and use nonnegativity of the original exponential pairing. This also
proves integrability of (RMA7).

This angular maximization is finite. One maximizing pair can be chosen
with at most \(6(n+1)\) switching points in each circle. To prove this,
add \(\epsilon\sum_i i\int\cos((n+1)\theta)1_{F_i}d\sigma\) and the
same linear term for \(G\). A soft maximizing pair exists, and an ordinary
one follows by separately taking extreme points of its linear maximizing
faces. With the other slice fixed, its winning scores minus mass prices
are trigonometric polynomials of degree \(n+1\). Every active pair
contrast has a nonzero highest harmonic from the added term, hence at
most \(2(n+1)\) zeros. There are at most three pairs. Zero-mass labels
are omitted, and a single active label needs no switching points.
Let \(\epsilon\downarrow0\): the finite family of label words and the
compact ordered-endpoint parameters have a convergent subsequence.
The limiting partition preserves its masses and attains the unperturbed
maximum. Thus (RMA6) is a finite endpoint optimization allowing multiple
arcs per label. It is not the false global consecutive-arc reduction.

For \(n=2\), write \(H_2(\cos\theta)=c_0+c_1\cos\theta+c_2\cos2\theta\).
Then
\[
 c_2=\frac{e^{uz_0}-e^{-u}(1+uD)}{2D^2},\quad
 c_1=ue^{-u}+4c_2,\quad 0\le\frac{c_2}{c_1+c_2}\le1/5.
 \tag{RMA8}
\]
The limit is \(c_2=u^2e^{-u}/4\) at \(D=0\).
The secant of this quadratic between \(-1\) and any
\(z_i=\cos(\pi(p_i+q_i))\ge z_0\) is
\(ue^{-u}+2c_2(1+z_i)\), at least the exponential secant
\((e^{uz_0}-e^{-u})/D\). The convex-majorant argument of DG9, now
applied to this quadratic, shows that (RMA6) for \(n=2\) is at least
that exponential secant times \(\Gamma_\sharp(p,q)\). Consequently
\[
                    \mathcal E_\rho^{[2]}(p)\le\mathcal E_\rho(p)
 \tag{RMA9}
\]
for every feasible profile, where the right side is (DG39).
This is a global improvement, although its remaining finite angular
optimization has not been evaluated for every profile.

## RMA3. Finite angular order dominates all radial directions quadratically

For each fixed \(n\), the kernels
\((H_n(\cos\theta)-H_n(-1))/(H_n(1)-H_n(-1/2))\)
lie in the compact convex hull of
\((1+\cos\theta)^k/(2^k-2^{-k})\), \(1\le k\le n\).
They have nonnegative Fourier coefficients. At equal masses, their
secant between \(-1\) and \(-1/2\) is at least
\(2/(4^n-1)>0\). The secant rearrangement therefore forces an extremal
cell to saturate its first angular moment, hence to be an arc.
Angular Cauchy--Schwarz forces a maximizing balanced pair to be the
same equal arcs. A fractional self-maximizer is excluded by a feasible
transfer with nonzero first harmonic. Compactness, bounded assignment
prices and the simple endpoint ties now give the same uniform local
arc-envelope argument as [JA2--JA4](THREE_CELL_JOINT_ANGULAR_ENVELOPE.md).
The kernels form a compact polynomial family; their balanced endpoint
derivative is controlled by the normalized difference \(K(0)-K(2\pi/3)=1\).
Their balanced arc scores are strictly decreasing in distance from
their center, directly by endpoint differentiation. Thus no extra
switching point survives in the limiting argument.

Subtracting the joint Hessian from the independent aligned-arc Hessian
gives, uniformly for fixed \(n\),
\[
 \mathcal D_{n,u}(p,q)=
 \frac{H_n(1)-H_n(-1/2)}{12}\|p-q\|^2
       +O_n\big((H_n(1)-H_n(-1/2))(|p-\mathbf1/3|+|q-\mathbf1/3|)^3\big)
\]
near the balanced pair. The formula uses the interpolating node of the
actual masses; changing it from \(z_0=-1/2+O(|p-\mathbf1/3|+|q-\mathbf1/3|)\)
affects only the cubic remainder. For bounded radial directions the
integrated remainder is uniform on compact correlation intervals,
because the coefficient derivatives are bounded by a polynomial in
\(u\) times \(e^u\), integrable against the radial Gaussian factor.

Hereafter set the interpolating node to \(-1/2\) and write
\(S_n(u)=H_n(1)-e^{-u/2}\), \(Q=d^{-1}e^{-\rho^2(t+s)/d}\).
The quadratic coefficient of (RMA7), per squared label-direction norm,
is
\[
 q_\rho^{[n]}(h)=\frac12(N_1+N_{-1/2}-D_1+D_{-1/2})
       -\frac1{12}\iint Q S_n(u)(h(t)-h(s))^2d\nu(t)d\nu(s).
 \tag{RMA10}
\]
The notation \(N_z,D_z\) is that of (DG40), extended by integration
to any scalar \(h\). The first term keeps the exact aligned-arc
interaction. One has \(S_n\uparrow e^u-e^{-u/2}\), so these quadratic
forms decrease to the exact joint Hessian (JA4). By (JA6) and the
standard-sector gap (S37), its limit satisfies, for mean-zero \(h\),
\[
 q_\rho^{[\infty]}(h)\le-k_I\int\sqrt t\,h(t)^2d\nu(t),\qquad k_I>0.
 \tag{RMA11}
\]
The standard sector contains no active-plane rotation. The speed
normalization is \(f(r)=r b(r)\), as in (JA6); its squared edge norm
is a positive fixed multiple of \(\int r|b(r)|^2d\mu(r)\).

There is a **finite** \(n_I\) and \(c_I'>0\) for which
\[
 q_\rho^{[n_I]}(h)\le-c_I'\int\sqrt t\,h(t)^2d\nu(t)
               \quad(\rho\in I,\ \int h\,d\nu=0).
 \tag{RMA12}
\]
The proof needs compactness rather than pointwise convergence of forms.
Write (RMA10) as \(\langle h,K_nh\rangle-\int D_nh^2d\nu\).
The off-diagonal kernel is
\(K_n=Q(e^u+e^{-u/2})/2+Q S_n/6\), and
\(D_n(t)=\int[Q(e^u-e^{-u/2})/2+QS_n/6]d\nu(s)\).
Uniformly on \(I\),
\(d_I\sqrt t\le D_n(t)\le C_I\sqrt t\), with \(d_I>0\).
For the lower bound on \(t\le1\), restrict \(\sqrt s\) to \([1,2]\)
and use \(e^u-e^{-u/2}\ge3u/2\). For \(x=\sqrt t\ge1\), complete
the square in \(Qe^u\) and restrict \(\sqrt s\) to
\([\rho x,\rho x+\sqrt{1-\rho^2}]\). Its row integral is at least
a constant times \(x\), and \(1-e^{-3u/2}\) is uniformly positive
there. The upper bound near zero follows by differentiating the
Gaussian row integral, and for \(x\ge1\) by its elementary completed
square bound \(C_I(1+x)\).

Under \(v=t^{1/4}h\), the normalized off-diagonal kernels converge
in Hilbert--Schmidt norm to \(K_\infty/(ts)^{1/4}\), uniformly on
\(I\), by the same all-correlation Gaussian domination as in (RMA3).
The diagonal multipliers \(D_n/\sqrt t\) increase pointwise and have
the common positive lower bound \(d_I\). If (RMA12) failed, there would
be \(n_j\to\infty\), \(\rho_j\to\rho_*\in I\) and unit vectors
\(v_j\rightharpoonup v\) with forms tending to zero or positive values.
The mean-zero condition passes to the limit since \(t^{-1/4}\in L^2(\nu)\).
Compactness gives convergence of the off-diagonal terms. For any fixed
order \(m\) and compact radial interval, weak lower semicontinuity
applied to \(D_m/\sqrt t-d_I\ge0\) gives, after exhausting the axis
and letting \(m\to\infty\),
\[
 \limsup q_{\rho_j}^{[n_j]}(h_j)
 \le q_{\rho_*}^{[\infty]}(t^{-1/4}v)-d_I(1-\|v\|^2)
 \le-k_I\|v\|^2-d_I(1-\|v\|^2)<0.
\]
This contradiction proves the claimed uniform finite order and gap.
For the varying correlations in this step, each fixed \(D_m\) converges
uniformly on compact positive radial intervals by dominated Gaussian
integration. This justifies the same lower-semicontinuity argument
before the two exhaustions.
Thus finitely many angular moments outweigh all positive scalar radial
directions **at quadratic order** on every compact positive interval.
No uniform order as \(\rho\uparrow1\), or finite-amplitude extension
of (RMA12), is asserted.

## RMA4. A directed isolated block at correlation 99/100

Let \(f\) be (DG21), and \(g\) the six-exponential function following
(DG40). Both are exactly mean zero and bounded. Displayed matrix entries
below are rounded; the assertions use directed enclosures. At \(\rho=99/100\),
the older coefficient form (DG40) on \(h=c_f f+c_g g\) has matrix
\[
 \begin{pmatrix}
 -0.03925812235691134&0.000008697386761891848\\
 0.000008697386761891848&0.000001202983953209141
 \end{pmatrix}
\]
with determinant in
\([-4.730253576497131,-4.730253576497129]\,10^{-8}\).
It has exactly one positive eigenvalue. In coefficient Euclidean norm
that eigenvalue is in \([1.2049107446686638,1.2049107446686642]10^{-6}\),
with \(c_f/c_g\) near \(0.000221536826206832\).

For order two at the balanced interpolating node,
\(S_2=(3/2)[10e^{-u/2}-(10+4u)e^{-u}]\). The resulting coefficient is
\[
 q_\rho^{[2]}(h)=\tfrac12(N_1+N_{-1/2}-D_1+D_{-1/2})
 -\tfrac14[10(D_{-1/2}-N_{-1/2})-10(D_{-1}-N_{-1})-4(D_a-N_a)],
 \tag{RMA13}
\]
where \(N_a=\iint Q u e^{-u}h(t)h(s)d\nu d\nu\) and \(D_a\) uses
\(h(t)^2\). These are the derivative terms in (DG24).
The directed matrix on the same radial space is
\[
 \begin{pmatrix}
 -0.039692869758179995&-0.000092691631595778\\
 -0.000092691631595778&-0.00003654343184932419
 \end{pmatrix}.
\]
The checker proves that this matrix plus \(I/30000\) is negative
definite: its upper-left entry is negative and its determinant lies in
\([1.1871928045439813,1.1871928045439816]10^{-7}\). Hence
\[
                q_{99/100}^{[2]}(c_f f+c_g g)
                        \le-(c_f^2+c_g^2)/30000.
 \tag{RMA14}
\]
This is a bound on the entire radial block, including its previously
positive direction, rather than a test of \(g\) alone. The written
uniform Taylor remainder gives one \(\epsilon_0>0\) for every unit
\((c_f,c_g)\): all ordinary partitions with circular masses
\(\mathbf1/3+\epsilon(1,-1,0)(c_f f+c_g g)\) strictly lose to the
simplex for \(0<|\epsilon|<\epsilon_0\), regardless of angular shapes.
The checker uses 160-bit Arb closed quadrant and derivative formulas,
with no radial cutoff or numerical quadrature. Its
[receipt](three_cell_radial_angular_modes_verification.json) records all
inputs, enclosures and limits. Replay from the root:

```sh
.venv/bin/python research/plurality-stability/check_three_cell_radial_angular_modes.py
```

## The remaining global comparison

Let \(U_\rho(p)=\sum_iq_\rho(p_i-1/3)\) and define
\(\mathcal A_\rho^{[2]}(p)=U_\rho(p)-\mathcal E_\rho^{[2]}(p)\ge0\),
using (RMA9), DG39 and the scalar arc envelope. Combining (RMA4) and
(RMA7) bounds the actual excess by
\[
 \sum_i\langle a_i,(B_\rho-2)a_i\rangle
              -\eta\sum_i\|g_i\|^2-\mathcal A_\rho^{[2]}(p).
\]
The requested global step is to prove that the angular loss and the
coercive remainder outweigh the finite adverse radial contribution for
**every feasible amplitude**, allowing a higher polynomial order where
needed. RMA12 proves the directional quadratic assertion uniformly;
RMA14 makes one adverse block fully quantitative. Neither assertion
proves this finite-amplitude inequality. The all-correlation target
\(\widetilde Q_\rho(e)\le2L_\rho(A;T)\) remains open.
