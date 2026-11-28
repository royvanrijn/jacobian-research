# A direct three-cell comparison and the angular information it requires

The unrestricted inequality
\(\widetilde Q_\rho(e)\leq2L_\rho(A;T)\) remains open here.
The exact expansion (B1) in the [continuation note](POSITIVE_CORRELATION_CONTINUATION.md)
makes it equivalent to \(S_\rho(A)\leq S_\rho(T)\), including for
partitions far from the simplex. The
[central shifted-clipping comparison](THREE_CELL_CENTRAL_SHIFT_BOOTSTRAP.md)
proves the unrestricted inequality through \(\rho=339/500\)
and on an unevaluated right neighborhood.
The independent direct comparison adopted in DG21 holds for every balanced three-cell partition in every
dimension \(d\geq2\) throughout \(0<\rho\leq61/100\).
Its [same-label variance proof](THREE_CELL_SAME_LABEL_VARIANCE_PENALTY.md)
uses the sharp arc envelope, a global angular variance reserve, and
positive radial features with directed tails. DG20 records the earlier
all-mass scalar comparison through \(53/100\), and DG17 gives a
finite-amplitude comparison through
\(221/500\): it isolates one adverse radial direction, limits its amplitude
by feasibility, and pays for its positive contribution with unavoidable
angular loss. Both arguments are global in partition distance.
The [variance and bootstrap note](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md)
records the unconditional variance absorption and its own certified covers.
The earlier [scalar-threshold argument](THREE_CELL_BEYOND_SCALAR_THRESHOLD.md),
adopted in DG15 below, crosses the first obstruction to its separated
scalar criterion qualitatively. It gives no further explicit numerical
endpoint beyond the current all-mass scalar comparison.
A separate [near-one result](THREE_CELL_NEAR_ONE_EXACT.md) records a qualitative
interval adjacent to one. The full all-positive-correlation statement remains
open; no overlap closing the intermediate range is established here. The argument also proves the comparison
at every positive correlation when the three labels occupy equal portions
of almost every circle. For general balanced partitions it isolates an
exact missing inequality and proves that discarding the angular
Cauchy--Schwarz remainder cannot establish it near correlation one.
The counterexample is to an upper-bound strategy, not to simplex optimality.
Propositions DG6, DG8 and DG9 retain progressively more angular information,
including a globally valid measure-dependent rearrangement bound with the
exact arc interaction. The certified restricted-profile comparisons and
limitations of those sufficient criteria do not extend the unrestricted cutoff.
The joint angular comparison in DG11 excludes all sufficiently small
squared circular-mass deviations on a fixed annulus, with arbitrary angular
shapes. Together with the near-one theorem, it confines any unresolved global
failure to a fixed intermediate correlation range and a definite annular
mass discrepancy.

## Primary results and scope

The closest inspected global theorem is Heilman,
[arXiv:2608.00333v1](https://arxiv.org/html/2608.00333v1), Theorem A.1,
which proves the balanced comparison for \(0<\rho\leq2/5\).
Corollary A.3 is an all-parameter angular inequality; the endpoint
\(2/5\) enters through the directed radial certificate in Lemma A.4.
Proposition DG4 adapts that certificate to \([2/5,21/50]\), with its
analytic tail domain checked explicitly. This is an application of that
method, not a new rearrangement mechanism or a historical-priority claim.

The earlier inspected text is Heilman,
[arXiv:2306.03312v3](https://arxiv.org/html/2306.03312v3), Theorem 1.5,
Sections 1.4, 5, 6 and 8. Its current theorem statement gives the
unrestricted positive comparison for \(0<\rho\leq0.1082\), not every
positive correlation. Section 1.4 identifies the competing radial-mean
and angular terms. The proof of Lemma 6.1 uses the angular
Cauchy--Schwarz step studied below.

The circular rearrangement input needed for the restricted comparison is
only this fixed-measure statement: if \(E\subset S^1\) has normalized
length \(p\), and \(I_p\) is an arc of length \(p\), then

\[
 \iint 1_E(\theta)1_E(\varphi)e^{\kappa\cos(\theta-\varphi)}
             d\sigma(\theta)d\sigma(\varphi)
 \leq
 \iint 1_{I_p}(\theta)1_{I_p}(\varphi)e^{\kappa\cos(\theta-\varphi)}
             d\sigma(\theta)d\sigma(\varphi),\qquad \kappa>0.
 \tag{DG0}
\]

The statement, hypotheses and its application were inspected in the
rearrangement step of Heilman's proof of Lemma 6.1. That paper attributes
the input to Baernstein--Taylor, Theorem 2; the 1976 original was not
inspected. We use (DG0) as the standard imported circular rearrangement
inequality. No numerical constant or variable-mass optimization from
Heilman's Corollary 4.2 is needed below.

For general ambient dimension, existence and reduction of an actual
maximizer to a plane follow from Heilman--Tarter, Lemma 2.3 and Theorem 1.9,
[Forum of Mathematics, Sigma 9 (2021), e65](https://doi.org/10.1017/fms.2021.56).
The theorem concerns maximizers, not a dimension reduction of an arbitrary
competitor. The polar identities below are planar; the equal-circle
comparison also has a direct cylindrical extension. These comparisons
are not a literature-exhaustiveness or priority claim.

## Exact separation of radial masses and angular changes

Fix \(0<\rho<1\), put \(s_\rho=\sqrt{1-\rho^2}\), and use
\(d\sigma=d\theta/(2\pi)\) and \(d\mu(r)=r e^{-r^2/2}dr\).
For an ordinary measurable three-cell partition of the plane write

\[
 a_i(r,\theta)=1_{A_i}(r\cos\theta,r\sin\theta),\quad
 p_i(r)=\int a_i(r,\theta)d\sigma(\theta),\quad
 b_i=p_i-\tfrac13,\quad g_i(r,\theta)=a_i(r,\theta)-p_i(r).
\]

Balance says \(\int b_i\,d\mu=0\); it does not say \(b_i=0\).
Define

\[
\begin{split}
 \kappa&=\frac{\rho r u}{s_\rho^2},\qquad
 I_0(\kappa)=\int e^{\kappa\cos\theta}d\sigma(\theta),\qquad
 k_\kappa(\theta)=\frac{e^{\kappa\cos\theta}}{I_0(\kappa)},\\
 (U_\kappa f)(\theta)&=\int k_\kappa(\theta-\varphi)f(\varphi)d\sigma(\varphi),\\
 P_\rho(r,u)&=s_\rho^{-2}
  \exp\!\left[-\frac{\rho^2(r^2+u^2)}{2s_\rho^2}\right]I_0(\kappa).
\end{split}
 \tag{DG1}
\]

Then \(P_\rho\) is the reversible Markov kernel for the radius under
the planar OU operator, relative to \(\mu\). Write
\(dW_\rho=P_\rho(r,u)d\mu(r)d\mu(u)\), a probability measure.
Let \(g^*\) denote the centered indicator vector of a fixed partition
of the circle into three equal arcs, and set

\[
\begin{split}
 H_\kappa(r,u)&=\sum_i\langle g_i(r),U_\kappa g_i(u)\rangle_{\sigma},
 &H_\kappa^*&=\sum_i\langle g_i^*,U_\kappa g_i^*\rangle_{\sigma},\\
 \mathcal R_\rho(A)&=\sum_i\langle b_i,P_\rho b_i\rangle_{\mu},\\
 \mathcal D_\rho(A)&=\int\left[H_\kappa^*
                -\frac{H_\kappa(r,r)+H_\kappa(u,u)}2\right]dW_\rho,\\
 \mathcal C_\rho(A)&=\frac12\int\sum_i
       \langle g_i(r)-g_i(u),U_\kappa[g_i(r)-g_i(u)]\rangle_{\sigma}
                                                            dW_\rho.
\end{split}
 \tag{DG2}
\]

In \(H_\kappa(r,r)\), the angular operator still has parameter
\(\kappa=\rho r u/s_\rho^2\), not \(\rho r^2/s_\rho^2\).
The diagonal notation refers only to its two angular arguments.

### Proposition DG1 (exact polar deficit)

For every ordinary balanced
three-cell partition of \(\mathbb R^2\),

\[
 S_\rho(T)-S_\rho(A)
    =\mathcal D_\rho(A)+\mathcal C_\rho(A)-\mathcal R_\rho(A).
 \tag{DG3}
\]

Both \(\mathcal C_\rho\) and \(\mathcal R_\rho\) are nonnegative.
Consequently the unrestricted direct comparison is exactly the remaining
inequality

\[
 \mathcal R_\rho(A)\leq\mathcal D_\rho(A)+\mathcal C_\rho(A).
 \tag{DG4}
\]

To prove the identity, write the planar Mehler density in polar
coordinates. The angular conditional density is \(k_\kappa\), so

\[
 S_\rho(A)=\int\left[\sum_i p_i(r)p_i(u)+H_\kappa(r,u)\right]dW_\rho
          =\tfrac13+\mathcal R_\rho(A)+\int H_\kappa(r,u)dW_\rho.
\]

The second equality uses balance and the Markov property. For the
simplex the corresponding expression is
\(1/3+\int H_\kappa^*dW_\rho\). Polarizing the self-adjoint
quadratic form of \(U_\kappa\) proves (DG3).
Its Fourier eigenvalues are \(I_m(\kappa)/I_0(\kappa)>0\), as follows
directly by expanding
\(e^{(\kappa/2)(z+z^{-1})}\):

\[
 I_m(\kappa)=\sum_{j\geq0}
   \frac{(\kappa/2)^{2j+|m|}}{j!(j+|m|)!}>0.
\]

This proves \(\mathcal C_\rho\geq0\); positivity of the OU quadratic
form proves \(\mathcal R_\rho\geq0\). All angular functions are
bounded and all operators are Markov contractions, so all terms are
finite. Fubini and the polarization identity apply without a limiting
regularity assumption. This proof needs no rearrangement theorem.

The qualitative self-comparison in Heilman's Lemma 6.1 also gives
\(\mathcal D_\rho\geq0\). That additional imported fact alone does
not compare its size with \(\mathcal R_\rho\).

## The comparison under equal circular masses

### Proposition DG2 (equal-circle comparison)

Suppose that
\(p_i(r)=1/3\) for every label and almost every radius. Then
\(S_\rho(A)\leq S_\rho(T)\), equivalently
\(\widetilde Q_\rho(e)\leq2L_\rho(A;T)\), for every \(0<\rho<1\).
The same holds in \(\mathbb R^d=V\oplus V^\perp\), \(\dim V=2\),
if the three labels occupy one third of almost every circle
\(\{(r\cos\theta,r\sin\theta,z):\theta\in S^1\}\).

Indeed (DG0), separately for each label, bounds
\(H_\kappa(r,r)\) by \(H_\kappa^*\), since all three arc lengths
are \(1/3\). Hence \(\mathcal D_\rho\geq0\), while
\(\mathcal R_\rho=0\); now apply (DG3).
For the higher-dimensional assertion, apply the same angular argument
to the two slices \((r,z),(u,w)\). The remaining Mehler factor in
\((z,w)\) is positive and integrates to one. No reduction of arbitrary
partitions or stationarity assumption is used.

In particular this settles partitions equivariant under a planar
rotation by \(2\pi/3\) that cyclically permutes the three labels and
fixes \(V^\perp\). It permits arbitrary measurable angular changes
and does not require invariance under reflection. It also settles every
equal-angle sector partition with a measurable radius-dependent rotation.
These are consequences of circular rearrangement, not a solution for
all balanced partitions. The three-seed cycles inside a four-cell
tetrahedron in the [swap note](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md)
have a different slice geometry and are not covered by this assertion.

## Why dropping the angular remainder fails

The angular Cauchy--Schwarz bound replaces
\(H_\kappa(r,u)\) by
\([H_\kappa(r,r)+H_\kappa(u,u)]/2\). Thus it gives the valid upper bound

\[
 S_\rho(A)\leq\mathcal B_\rho(A)
 :=\tfrac13+\mathcal R_\rho(A)
       +\frac12\int[H_\kappa(r,r)+H_\kappa(u,u)]dW_\rho,
 \qquad \mathcal B_\rho(A)-S_\rho(A)=\mathcal C_\rho(A).
 \tag{DG5}
\]

### Proposition DG3 (failure of the angular Cauchy relaxation)

There
is an exactly balanced family of ordinary planar partitions, with each
label a single arc on every circle, for which
\(\mathcal B_\rho(A)>S_\rho(T)\) for every sufficiently large fixed
\(\rho<1\) and all sufficiently small nonzero family parameters.
Thus the stronger assertion
\(\mathcal R_\rho\leq\mathcal D_\rho\) is false, even in this class.

Put

\[
 b(r)=e^{-r^2/2}-\tfrac12,\qquad c=(1,-1,0),\qquad
 p_i(r)=\tfrac13+t c_i b(r).
 \tag{DG6}
\]

Take consecutive arcs of these lengths, with endpoints
\(0,2\pi p_1(r),4\pi/3,2\pi\). For \(|t|<2/3\) they define an
ordinary partition. Since \(\int b\,d\mu=0\), every cell has mass
exactly \(1/3\). At \(t=0\) the partition is \(T\).

For an arc of normalized length \(p\), its centered self-energy is

\[
 F_\kappa(p)=2\int_0^p(p-v)k_\kappa(2\pi v)dv-p^2,
 \qquad F_\kappa''(p)=2[k_\kappa(2\pi p)-1].
 \tag{DG7}
\]

Therefore, at each fixed \(\rho<1\),

\[
 \mathcal B_\rho(A(t))-S_\rho(T)
   =t^2\|c\|^2 E_\rho+o(t^2),\qquad
 E_\rho=\langle b,P_\rho b\rangle
                  -\int[1-m_\rho(r)]b(r)^2d\mu(r),
 \tag{DG8}
\]

where

\[
 m_\rho(r)=\int P_\rho(r,u)k_\kappa(2\pi/3)d\mu(u).
\]

The first derivatives cancel because \(\sum_i c_i=0\).
Differentiation under the integral is valid: near \(p=1/3\),
\(k_\kappa(2\pi p)\) and its needed angular derivatives are bounded
uniformly in \(\kappa\geq0\), and \(b\) is bounded.

Gaussian integration gives

\[
 \|b\|_\mu^2=\tfrac1{12},\qquad
 \langle b,P_\rho b\rangle=\frac1{4-\rho^2}-\tfrac14,\qquad
 \langle b,(I-P_\rho)b\rangle=\frac{s_\rho^2}{3(4-\rho^2)}.
 \tag{DG9}
\]

For example the joint Laplace transform of the two correlated planar
Gaussians gives
\(\mathbb E e^{-(|X|^2+|Y|^2)/2}=1/(4-\rho^2)\).
With \(a=\rho r/s_\rho\), the remaining multiplier is

\[
 m_\rho(r)=M(a),\qquad
 M(a)=e^{-a^2/2}\int_0^\infty v e^{-v^2/2-av/2}dv.
 \tag{DG10}
\]

In particular \(0<M(a)\leq e^{-a^2/2}\). Substitution
\(r=s_\rho a/\rho\) and dominated convergence show that

\[
 \lim_{\rho\uparrow1}\frac1{s_\rho^2}\int m_\rho b^2d\mu
   =\tfrac14 C,\qquad
 C=\iint_{(0,\infty)^2}av\,e^{-(a^2+v^2+av)/2}da\,dv
   =\frac43-\frac{4\pi}{9\sqrt3}.
 \tag{DG11}
\]

The last integral can be evaluated by differentiating the elementary
quadrant Gaussian integral in its cross coefficient, or by polar
coordinates. For completeness the general formula is

\[
 J(A,c):=\iint uv\,e^{-(Au^2+v^2+2cuv)/2}du\,dv
 =\frac1{A-c^2}-\frac{c\arccos(c/\sqrt A)}{(A-c^2)^{3/2}},
 \quad A>c^2,\ c\geq0.
\]

Here both integrals are on the positive quadrant and \(C=J(1,1/2)\).
Combining (DG8)--(DG11) proves

\[
 \boxed{\displaystyle
 \lim_{\rho\uparrow1}\frac{E_\rho}{1-\rho}
    =\frac{4-2\pi/\sqrt3}{9}>0.}
 \tag{DG12}
\]

Positivity is the elementary inequality \(\pi<2\sqrt3\).
First fix any sufficiently large \(\rho<1\), and then take sufficiently
small \(t\ne0\) in (DG8). This proves Proposition DG3 with the stated
order of quantifiers. It asserts nothing adverse about the actual
stability of \(A(t)\), since (DG5) has discarded \(\mathcal C_\rho\).

An exact finite-correlation expression, useful for checking the
normalization independently of the limit, is

\[
 E_\rho=s_\rho^2\left[
    J(1+2s_\rho^2,\rho/2)-J(1+s_\rho^2,\rho/2)
      +\tfrac14J(1,\rho/2)-\frac1{3(4-\rho^2)}\right].
 \tag{DG13}
\]

It follows by expanding \(b^2=e^{-r^2}-e^{-r^2/2}+1/4\) in the
quadrant integral for \(\int m_\rho b^2d\mu\).
Unlike the independent-arc relaxation in the separate
[global obstruction note](THREE_CELL_GLOBAL_OBSTRUCTIONS.md), the
partition here remains disjoint on every circle. The loss occurs in
replacing interactions between different circles by angular self-energies.

## A direct global interval through twenty-one fiftieths

### Proposition DG4 (global comparison through rho equal to 21/50)

For every \(d\geq2\), every ordinary measurable Gaussian three-cell
partition with cell masses \(1/3\), and every \(0<\rho\leq21/50\),

\[
 S_\rho(A)\leq S_\rho(T),\qquad
 \widetilde Q_\rho(1_A-1_T)\leq2L_\rho(A;T).
 \tag{DG14}
\]

This is a direct global comparison on the stated interval. It does not
use the continuation theorem or a hypothesis that \(A\) is close to
the simplex. The proof imports Heilman's sharpened aligned-arc inequality
and maximizer dimension reduction, and uses a new directed interval run.

Here are the exact analytic inputs and their range. Put
\(d_\rho=1-\rho^2\), \(\kappa_\rho=\rho^2/d_\rho\), and use
\(d\nu(t)=e^{-t}dt\). The planar radius is \(r=\sqrt{2t}\).
For \(a\geq0\) define

\[
\begin{split}
 j_1&=\frac{9-\pi\sqrt3}{4\pi^2}>0,\qquad
 \theta(a)=\frac{j_1 a e^{-a}}{2I_0(a)},\\
 \tau_+(a)&=\frac1{2I_0(a)}\int_0^{2/3}(2+3v)e^{a\cos(\pi v)}dv,\\
 w(a)&=\frac1{2I_0(a)}\int_0^{2/3}(1-3v)e^{a\cos(\pi v)}dv,\\
 a(t,s)&=\frac{2\rho\sqrt{ts}}{d_\rho},\qquad
 R_0(t,s)=\frac{e^{-\kappa_\rho(t+s)}}{d_\rho}I_0(a(t,s)),\\
 K_+(t,s)&=R_0(t,s)[\tau_+(a(t,s))-4\theta(a(t,s))],\\
 J_+(t)&=\int R_0(t,s)[w(a(t,s))+\theta(a(t,s))]d\nu(s).
\end{split}
 \tag{DG15}
\]

Write \(\mathcal K_+\) for this kernel operator, and put
\(z_i(t)=p_i(\sqrt{2t})-1/3\). The fixed masses give
\(\int z_i\,d\nu=0\). Circular rearrangement and Corollary A.3,
equation (142), in the inspected primary text imply

\[
 S_\rho(A)-S_\rho(T)
 \leq\sum_i\left[\langle z_i,\mathcal K_+z_i\rangle_\nu
                    -2\int J_+z_i^2d\nu\right].
 \tag{DG16}
\]

Corollary A.3 is stated for every angular parameter \(a\geq0\) and
every pair of three-component probability vectors. Consequently this
reduction is valid for every fixed \(0<\rho<1\), whenever its integrals
are evaluated, rather than only the interval of Theorem A.1. The proof
of Lemma A.2, the first-mode decomposition (140)--(142), and the
assembly of (146) were inspected. We import that inequality with partial
proof-audit coverage; no independent review of the preprint is asserted.

Let \(\ell_\rho\) be the exact binary64 polynomial specified by
`CENTER` in the retained upstream positive certificate. No approximation
property of this choice is assumed. Set
\(\widetilde J_+(x)=J_+(x^2)/x\), using its positive continuous limit
at zero. Since additive centering vanishes on mean-zero functions,
the Hilbert--Schmidt estimate gives

\[
\begin{split}
 \langle z,\mathcal K_+z\rangle_\nu
   &\leq h_+(\rho)\int J_+z^2d\nu,\\
 h_+(\rho)^2
   &=4\iint_{(0,\infty)^2}e^{-x^2-y^2}
       \frac{[K_+(x^2,y^2)-\ell_\rho(x)-\ell_\rho(y)]^2}
            {\widetilde J_+(x)\widetilde J_+(y)}dx\,dy.
\end{split}
 \tag{DG17}
\]

This is the same calculation as equations (151)--(152) of the source.
Its use of the mean-zero constraint is essential. Boundedness of the
\(z_i\), positivity of \(J_+\), and the finite bound below justify
the integrals and the weighted Cauchy--Schwarz step.

The full directed run, on the entire interval
\([2/5,21/50]\), gives the conservative bounds

\[
                  h_+(\rho)^2<3.999335108,
             \qquad h_+(\rho)<1.999834<2.
 \tag{DG18}
\]

It accepted all 85 covering correlation subintervals after 166 attempts.
The run used 3200 radial panels per coordinate on \([0,7]^2\), 64
second-order angular midpoint panels, and the same Arb precision and
auxiliary table sizes as the original certificate. All transcendental
tables and polynomial-center coefficients were reevaluated over the new
correlation interval. The displayed source output was rounded to twelve
significant digits; (DG18) rounds it further upward.

The change of range also needs an analytic check. The source's lower
bound and envelope, valid here, are

\[
\begin{split}
 \widetilde J_+(x)&\geq c_+\frac{\rho}{\sqrt{d_\rho}}
                         e^{-\kappa_\rho x^2},
 &c_+&=\frac{6-\pi}{4\pi^{3/2}}>0,\\
 |K_+(x^2,y^2)|&\leq\frac{1+2j_1/e}{d_\rho}
       e^{-\kappa_\rho(x^2+y^2)+(2\rho/d_\rho)xy}.
\end{split}
 \tag{DG19}
\]

After division by the weights, the squared-kernel tail is bounded by
a constant times
\(\exp[-(x^2+y^2-4\rho xy)/d_\rho]\).
Completing the square leaves a positive Gaussian rate because
\(\rho\leq21/50<1/2\); at the upper endpoint that rate is
\(736/2059\). For the polynomial-center tails the rate is
\(1-\kappa_\rho\geq1618/2059>0\). Thus both tail arguments
remain valid beyond the original advertised endpoint. At the slab with
the largest certified total, the reported squared-norm tail contributions were
\(8.5647947864\times10^{-5}\) and
\(2.44677710152\times10^{-11}\), respectively. Endpoint expansion by
one binary64 step encloses both exact rational endpoints. Failed slabs
were subdivided on both sides; no unaccepted slab or sampled-only gap
remained when the process exited successfully.

Combining (DG16)--(DG18) proves the planar comparison on
\([2/5,21/50]\), with the stronger radial-mass estimate

\[
 S_\rho(T)-S_\rho(A)
       \geq(2-h_+(\rho))\sum_i\int J_+(t)z_i(t)^2d\nu(t).
 \tag{DG20}
\]

Theorem A.1 of the source supplies the interval below \(2/5\).
For each ambient dimension, take an actual maximizing partition using
Heilman--Tarter Lemma 2.3, and apply their Theorem 1.9 to reduce that
maximizer to a plane. The planar bound then bounds every competitor in
that ambient dimension. This is a reduction of a maximizer, not of an
arbitrary competitor. Finally apply the exact quadratic expansion (B1)
to get the second inequality in (DG14).

The [verification record](positive_correlation_extension_verification.json)
pins the complete command, environment, exact source hashes, interval
inputs, output and assurance limits. The two original upstream modules
are retained unchanged in
[the source directory](sources/cut_bounds_6b720fe57bfe/plur_cert_positive.py),
from [the source repository](https://github.com/sheilman77/cut_bounds)
at commit `6b720fe57bfe2a3c8fffc51521d00b31011aa1a3`. The short
[extension wrapper](sources/cut_bounds_6b720fe57bfe/check_positive_extension.py)
changes only the correlation-range constants and mesh arguments. The
executed ignored copies and retained copies have identical bytes.
Replay from the repository root with the existing environment:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/sources/cut_bounds_6b720fe57bfe/check_positive_extension.py
```

This is an internally executed computer-assisted proof under the stated
analytic inputs. It is not formal verification, an independent replay,
external mathematical review, or an unrestricted positive-correlation
result. In particular the tail estimate's validity below \(1/2\)
does not prove the norm bound on that larger interval.

## The scalar radial bound fails at correlation one half

### Proposition DG5 (a bounded radial obstruction at rho equal to 1/2)

For the exact operators \(K_+,J_+\) in (DG15), at \(\rho=1/2\) put

\[
 f(t)=e^{-t}+2e^{-16t}-\frac{21}{34},\qquad
 N=\langle f,\mathcal K_+f\rangle_\nu,\qquad
 D=\int J_+f^2d\nu.
 \tag{DG21}
\]

Then \(f\) is bounded, \(\int f\,d\nu=1/2+2/17-21/34=0\), and

\[
           D>0,\qquad N-2D>0.0015,\qquad N/D>2.07.
 \tag{DG22}
\]

Thus no choice of additive center, or sharper estimate of the same
scalar operator norm, can extend the sufficient condition
\(\langle f,\mathcal K_+f\rangle\leq2\int J_+f^2d\nu\) to
\(\rho=1/2\). This is not a partition beating the simplex. It is a
counterexample to a sufficient radial estimate used after relaxation.

Here is an exact reduction of (DG22) to a one-dimensional interval
calculation. For \(P>c^2\) define

\[
 \Xi(P,c)=\frac1{P-c^2}
       -\frac{c\arccos(c/\sqrt P)}{(P-c^2)^{3/2}}.
 \tag{DG23}
\]

For \(A,B>0\), \(P=AB\), this is the positive-quadrant integral
\(\iint uv\exp[-(Au^2+Bv^2+2cuv)/2]du\,dv\).
The formula is valid for either sign of \(c\). It follows by scaling
the quadrant Gaussian integral and differentiating in \(c\).
For nonnegative \(\lambda,\mu\), put
\(P_{\lambda\mu}=(1+\lambda d_\rho)(1+\mu d_\rho)\). Directly
substituting \(t=d_\rho u^2/2\), \(s=d_\rho v^2/2\) gives

\[
 \iint R_0(t,s)\frac{e^{a(t,s)z}}{I_0(a(t,s))}
           e^{-\lambda t-\mu s}d\nu(t)d\nu(s)
       =d_\rho\Xi(P_{\lambda\mu},-\rho z).
\]

Write \(F_{\lambda\mu}(v)=d_\rho
\Xi(P_{\lambda\mu},-\rho\cos(\pi v))\), and set

\[
\begin{split}
 B_{\lambda\mu}
  &=\frac{j_1\rho d_\rho}{2}
    \left[\frac{(P_{\lambda\mu}+2\rho^2)
              \arccos(\rho/\sqrt{P_{\lambda\mu}})}
                    {(P_{\lambda\mu}-\rho^2)^{5/2}}
            -\frac{3\rho}{(P_{\lambda\mu}-\rho^2)^2}\right],\\
 N_{\lambda\mu}
  &=\frac12\int_0^{2/3}(2+3v)F_{\lambda\mu}(v)dv-4B_{\lambda\mu},\\
 D_{\lambda\mu}
  &=\frac12\int_0^{2/3}(1-3v)F_{\lambda\mu}(v)dv+B_{\lambda\mu}.
\end{split}
 \tag{DG24}
\]

The first expression is the pairing of the kernel
\(R_0\theta\) against the two exponentials: use
\(\theta(a)=(j_1/2)ae^{-a}/I_0(a)\) and differentiate the preceding
identity in \(z\) at \(z=-1\). For
\(\lambda=(0,1,16)\), \(c=(-21/34,1,2)\), the required values are

\[
 N=\sum_{i,j}c_i c_jN_{\lambda_i\lambda_j},\qquad
 D=\sum_{i,j}c_i c_jD_{\lambda_i+\lambda_j,0}.
 \tag{DG25}
\]

All radial integrals have been evaluated exactly; no radial cutoff or
unbounded numerical tail is used here. The
[coupled-radial checker](check_three_cell_coupled_radial.py) evaluates
the remaining angular integrals by Arb interval Riemann enclosures over
all 8192 subintervals of \([0,2/3]\), at 100-bit precision. The
[retained verification](three_cell_coupled_radial_verification.json)
gives outward enclosures, in particular
\(0.00156437027<N-2D<0.00191546337\). This proves (DG22).

## An explicit correction for angular moment closure

The [first-moment compatibility argument](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#circular-balance-and-compatibility-of-first-angular-moments)
defines a positive deficit by a compact vector optimization. The next
bound replaces that optimization by an elementary expression. It applies
to every pair of actual circular partitions, at every positive angular
parameter.

### Proposition DG6 (an explicit angular compatibility penalty)

Let \(p,q\in\Delta_3\). Initially assume all their coordinates are
positive. Put

\[
 a_i=\frac{\sin(\pi p_i)}\pi,\qquad
 b_i=\frac{\sin(\pi q_i)}\pi,\qquad r_i=b_i/a_i.
\]

Let \(G(a)\) be the symmetric matrix with
\(G(a)_{ii}=a_i^2\) and
\(G(a)_{ij}=(a_k^2-a_i^2-a_j^2)/2\) for distinct \(i,j,k\).
Define

\[
\begin{split}
 M(r)&=\max\left(0,\max_{\{i,j,k\}=\{1,2,3\}}
                         \frac{(r_i-r_j)(r_i-r_k)}{r_i}\right),\\
 \Gamma_0(a,b)&=\frac{r^TG(a)r}{2(\sum_i r_i+M(r))},\\
 \Gamma(p,q)&=\max\{\Gamma_0(a,b),\Gamma_0(b,a)\}.
\end{split}
 \tag{DG26}
\]

Set \(\Gamma(p,q)=0\) if either profile has a zero coordinate.
For two circular partitions with masses \(p,q\), let \(u_i,v_i\)
be their first angular vector moments, with probability measure
\(d\theta/(2\pi)\). Then

\[
 \sum_i a_i b_i-\sum_i u_i\cdot v_i\ \geq\ \Gamma(p,q)\geq0.
 \tag{DG27}
\]

For interior profiles, \(\Gamma(p,q)>0\) exactly when \(p\ne q\).
For fixed vectors \(x,y\) of coordinate sum zero, as \(\epsilon\to0\),

\[
 \Gamma(\tfrac13\mathbf1+\epsilon x,
        \tfrac13\mathbf1+\epsilon y)
       =\frac{\epsilon^2}{16}\|x-y\|^2+O(\epsilon^3).
 \tag{DG28}
\]

**Proof.** Projection in the direction of the angular moment gives
\(|u_i|\leq a_i\), \(|v_i|\leq b_i\), by selecting the largest
cosine values. Since the slices are partitions,
\(\sum_i u_i=\sum_i v_i=0\). Write
\(\Delta=\sum_i a_i b_i-\sum_i u_i\cdot v_i\),
\(d_i=a_i^2-|u_i|^2\), and \(e_i=b_i^2-|v_i|^2\).
The exact square identity is

\[
 2\Delta=\sum_i\left|\sqrt{r_i}u_i-v_i/\sqrt{r_i}\right|^2
                      +\sum_i r_i d_i+\sum_i e_i/r_i.
\]

In particular all three terms on the right are nonnegative. Closure
and weighted Cauchy--Schwarz imply
\( |\sum_i r_i u_i|^2\leq2(\sum_i r_i)\Delta\).
For a closing triple, its Gram matrix is determined by its three
squared lengths. With
\(\eta_i=(r_i-r_j)(r_i-r_k)\), this gives

\[
 r^TG(a)r=\left|\sum_i r_i u_i\right|^2+\sum_i\eta_i d_i
 \leq2(\sum_i r_i)\Delta+M(r)\sum_i r_i d_i
 \leq2(\sum_i r_i+M(r))\Delta.
\]

This proves the first one-sided bound. Interchanging the triples gives
the other. The matrix \(G(a)\) is positive semidefinite: take the
angular moments of consecutive arcs of widths \(p_i\). They close
and have lengths \(a_i\), so their Gram matrix is exactly \(G(a)\).
For interior profiles they form a nondegenerate triangle, with kernel
\(\operatorname{span}(\mathbf1)\). Thus its quadratic form vanishes
only if all \(r_i\) agree. The side lengths then define similar
triangles with angles \(\pi p_i\), \(\pi q_i\); uniqueness of their
ordered angle triples gives \(p=q\). The boundary definition uses
only the trivial nonnegative bound from the square identity or individual
Cauchy--Schwarz, and makes no strictness claim there.

Finally at the balanced profile, with \(a_*=\sqrt3/(2\pi)\),
\(G(a)=\tfrac32a_*^2P+O(\epsilon)\), where
\(P=I-\mathbf1\mathbf1^T/3\), and
\(r_i=1+\epsilon(y_i-x_i)/(2a_*)+O(\epsilon^2)\).
The numerator of \(\Gamma_0\) is therefore
\(3\epsilon^2\|x-y\|^2/8+O(\epsilon^3)\), while its denominator
is \(6+O(\epsilon^2)\). The same leading term holds for the reverse
bound, proving (DG28). The expansion is uniform for bounded \(x,y\).

To apply this to the noise kernel, write

\[
 k_a(\vartheta)=c_a\cos\vartheta+R_a(\vartheta),\qquad
 c_a=\frac{ae^{-a}}{I_0(a)}.
\]

As a function of \(z=\cos\vartheta\), the remainder is
\([e^{az}-ae^{-a}z]/I_0(a)\), nonnegative and nondecreasing on
\([-1,1]\). Its derivative is \(a(e^{az}-e^{-a})/I_0(a)\).
Theorem I of
[Friedberg--Luttinger](https://digital.library.unt.edu/ark:/67531/metadc1025197/m2/1/high_res_d/4446751.pdf),
printed pp. 1--2, equations (2)--(5), therefore bounds each same-label
pairing for the remainder by aligned arcs. The theorem and its
nonnegativity/monotonicity hypotheses were inspected in the primary
report; use two functions and set the other kernel factor to one.
Finite unions of arcs approximate arbitrary measurable sets in \(L^1\),
and the bounded kernel makes this extension continuous. The separated
cosine pairing is exactly \(u_i\cdot v_i\), while aligned arcs give
\(a_i b_i\). Hence the actual angular pairing is at most its
aligned-arc envelope minus \(c_a\Gamma(p,q)\).

Define the positive radial kernel
\(V_\rho(t,s)=R_0(t,s)a(t,s)e^{-a(t,s)}/I_0(a(t,s))\).
Applying Corollary A.3 only to the remaining envelope gives the explicit
global upper bound

\[
\begin{split}
 S_\rho(A)-S_\rho(T)\leq{}&
 \sum_i[\langle z_i,\mathcal K_+z_i\rangle-2\int J_+z_i^2d\nu]\\
 &-\iint V_\rho(t,s)\Gamma(p(\sqrt{2t}),p(\sqrt{2s}))d\nu(t)d\nu(s).
\end{split}
 \tag{DG29}
\]

Every term is finite at fixed \(0<\rho<1\). In particular
\(c_a\leq1/e\), and the penalty is bounded by the actual first-moment
deficit, which is at most \(2\sum_i a_i b_i\leq6/\pi^2\).
Thus its integral is controlled by the radial Markov measure. No
stationarity or closeness of the partition is assumed.

### Corollary DG7 (the coupled bound handles the explicit profile at rho equal to 1/2)

At \(\rho=1/2\), there exists \(\epsilon_0>0\) such that every
ordinary planar partition with circular mass profile

\[
 p_i(\sqrt{2t})=\tfrac13+\epsilon(1,-1,0)_i f(t),
                  \qquad 0<|\epsilon|<\epsilon_0,
 \tag{DG30}
\]

strictly loses to the simplex. This permits arbitrary angular shapes,
not only consecutive arcs. The profiles are balanced by (DG21) and
are feasible, for example by consecutive arcs, for sufficiently small
\(\epsilon\), since \(|f|\leq81/34\).

To check the sign, define
\(N_B=\sum c_i c_jB_{\lambda_i\lambda_j}\) and
\(D_B=\sum c_i c_jB_{\lambda_i+\lambda_j,0}\), with the same
coefficients and exponents as in (DG25). Since
\(R_0\theta=(j_1/2)V_\rho\), (DG28) shows that the coefficient of
\(2\epsilon^2\) in the right side of (DG29) is

\[
       N-2D-\frac{D_B-N_B}{4j_1}<-0.0038.
 \tag{DG31}
\]

The directed checker certifies this inequality as well as (DG22).
Its enclosure for the subtracted correction is
\([0.005810135707537683,0.005810135707537685]\), and its enclosure
for the corrected coefficient lies inside
\([-0.004246,-0.003894]\). The remainder is \(O(\epsilon^3)\)
uniformly over partitions with profile (DG30): (DG28) is uniform for
bounded profiles near the balanced one, \(f\) is bounded, and
\(V_\rho\) has finite integral. Therefore the strictly negative
quadratic coefficient proves the assertion.

The checker also performs exact rational controls of the Gram identity,
both one-sided bounds, the coincident-triple case and a rejected
nonclosing triple. These finite controls are supplementary to the
unrestricted algebraic proof of Proposition DG6. The all-profile
claim that the right side of (DG29) is always nonpositive is false.
At \(\rho=3/5\), the same checker and the same mean-zero function give

\[
       N-2D-\frac{D_B-N_B}{4j_1}>0.0038.
\]

The enclosure is
\([0.0038888594281296417,0.0043384122595988344]\).
Consequently the upper bound in (DG29) is positive for small nonzero
profiles (DG30) at that correlation. This does not establish the sign of
any actual stability deficit: (DG29) remains a valid upper bound.
It shows that its particular elementary penalty is not sufficient for
the full goal. Corollary DG7 does not extend the unrestricted cutoff in
Proposition DG4.

There is a concrete algebraic refinement to pursue. On
\(H=\mathbf1^\perp\), let \(X\) be positive definite and extend it
and \(X^{-1}\) by zero on \(\mathbf1\). Put
\(\lambda_i(X)=2X_{ii}+X_{jk}\) for distinct \(i,j,k\).
If all \(\lambda_i(X)\) and \(\lambda_i(X^{-1})\) are nonnegative,
the same length constraints give the additional bound

\[
 \Delta\geq\sum_i a_i b_i
       -\tfrac12[\operatorname{tr}_H(XG(a))
                         +\operatorname{tr}_H(X^{-1}G(b))].
 \tag{DG32}
\]

Here \(G(u)\) denotes the actual Gram matrix \((u_i\cdot u_j)\).
Matrix Cauchy--Schwarz gives
\(\sum_i u_i\cdot v_i\leq
\tfrac12[\operatorname{tr}(XG(u))+\operatorname{tr}(X^{-1}G(v))]\).
For a closing triple,
\(\operatorname{tr}(XG(u))=\sum_i\lambda_i(X)|u_i|^2\);
the nonnegative coefficients allow replacement by \(a_i^2\), and
similarly for \(v\). Taking a maximum with zero and with (DG26) is
therefore valid. The following refinement makes a profile-dependent
choice explicit, handles possible negative coefficients, and computes
its exact local strength. A uniform radial comparison is still absent.

### Proposition DG8 (sharp first-moment penalty near balanced profiles)

For interior profiles, put \(A=G(a)|_H\), \(B=G(b)|_H\),
\(h=\sqrt{\det_H A\det_H B}\), and

\[
 F=\sqrt{\operatorname{tr}_H(AB)+2h},\qquad
 X=\frac{B+hA^{-1}}{F}.
 \tag{DG33}
\]

These matrices are positive definite on \(H\). Extend \(X,X^{-1}\)
by zero on \(\mathbf1\), and use the coefficients from (DG32).
Writing \(c_+=\max(c,0)\), define

\[
 \Gamma_\sharp(p,q)=\max\!\left\{\Gamma(p,q),\,
 \sum_i a_i b_i-\frac12\sum_i
  [\lambda_i(X)_+a_i^2+\lambda_i(X^{-1})_+b_i^2]\right\}.
 \tag{DG34}
\]

Set \(\Gamma_\sharp=0\) on the boundary. This is a nonnegative,
globally valid first-moment deficit bound, so (DG29) remains valid when
\(\Gamma\) is replaced by \(\Gamma_\sharp\). At the balanced profile,
uniformly for bounded zero-sum \(x,y\),

\[
 \Gamma_\sharp(\mathbf1/3+\epsilon x,
                    \mathbf1/3+\epsilon y)
       =\frac{\epsilon^2}{8}\|x-y\|^2+O(\epsilon^3).
 \tag{DG35}
\]

This coefficient is optimal if the only constraints used on the first
moments are closure and their individual norm bounds. In fact, for all
\(p,q\) in a sufficiently small fixed neighborhood of \(\mathbf1/3\),
the exact minimum deficit under those constraints is
\(\sum_i a_i b_i-F\).

The fixed-Gram step is the two-matrix orthogonal Procrustes problem.
[Ten Berge (1977)](https://three-mode.leidenuniv.nl/pdf/t/tenberge1977pmet.pdf),
Theorem 2, printed p. 269, was inspected: the matrices are real and the
rotation ranges over all orthogonal matrices; the optimum makes the
cross product symmetric positive semidefinite. That result does not
supply the individual norm inequalities needed here. The calculation
and the additional sign issue are proved explicitly below; no novelty
is claimed for the matrix optimization.

For positive \(2\times2\) matrices, Cayley--Hamilton gives
\(XAX=B\) and
\(\operatorname{tr}(XA)=\operatorname{tr}(X^{-1}B)=F\).
For example, apply Cayley--Hamilton to
\(A^{1/2}BA^{1/2}\); its positive square root is
\((A^{1/2}BA^{1/2}+hI)/F\). Matrix Cauchy as in (DG32) still
applies when some length coefficients are negative. In that case use
\(\lambda_i(X)|u_i|^2\leq\lambda_i(X)_+a_i^2\), and likewise
for \(v\). This proves (DG34) without a sign assumption. Taking its
maximum with the already valid \(\Gamma\) preserves the inequality.
The bound is measurable, and it is bounded by the same actual deficit
as \(\Gamma\), so the integrability argument for (DG29) still applies.

At \(p=q=\mathbf1/3\), \(X=I_H\) and all six length coefficients
are one. By continuity they remain strictly positive near this pair.
In that neighborhood (DG34)'s second entry is \(\sum_i a_i b_i-F\).
Choose closing two-dimensional vector triples with Gram matrices
\(G(a),G(b)\), which exist because these are Gram matrices of the
three arc moments. A singular value decomposition of the cross matrix
shows that an orthogonal transformation of the second triple attains
pairing \(F\): the maximum is the sum of its two singular values,
and the square of this sum is \(\operatorname{tr}(AB)+2h\). Hence the
norm-and-closure relaxation attains the bound, proving its stated
optimality. This also shows that this entry is at least \(\Gamma\)
in that neighborhood. The optimizing moment triples themselves are
realizable by consecutive arc partitions, after a rotation or reversal
of one circle; the remainder kernel need not attain its rearrangement
envelope at this relative orientation.

Here are details of the local expansion. Write
\(a_*=\sqrt3/(2\pi)\), \(s=3a_*^2/2\), and let \(M(x)\)
be the symmetric matrix with diagonal \(a_*x_i\) and off-diagonal
entry \(M(x)_{ij}=a_*x_k\), for distinct indices. Then
\(A=sI_H+\epsilon M(x)|_H+O(\epsilon^2)\),
\(\operatorname{tr}M(x)=0\), and
\(\|M(x)-M(y)\|_F^2=3a_*^2\|x-y\|^2\).
Expanding (DG33), or its two-dimensional square-root formula, gives

\[
 F=\tfrac12\operatorname{tr}(A+B)
             -\frac{\|A-B\|_F^2}{8s}+O(\epsilon^3).
\]

Also
\(\operatorname{tr}G(a)=3a_*^2-\epsilon^2\|x\|^2/2+O(\epsilon^3)\)
and
\(\sum_i a_i b_i=3a_*^2+\epsilon^2[x\cdot y/4
-3(\|x\|^2+\|y\|^2)/8]+O(\epsilon^3)\).
Substitution proves (DG35). The expansions are uniform on bounded
tangent sets because the matrices stay in a compact positive definite
neighborhood and the coefficient signs stay positive.

For the profile (DG30), replacing \(\Gamma\) by
\(\Gamma_\sharp\) therefore replaces (DG31)'s coefficient by

\[
             N-2D-\frac{D_B-N_B}{2j_1}.
 \tag{DG36}
\]

At \(\rho=3/5\), its directed enclosure is
\([-0.0019184952902936525,-0.001468942458824461]\).
Consequently every partition with (DG30)'s circular mass profile
strictly loses to the simplex for sufficiently small nonzero
\(\epsilon\) at that correlation, with arbitrary angular shapes.
The same uniform Taylor argument as in DG7 proves this assertion.

Even this sharper separated-first-moment criterion is not sufficient
at every correlation. At \(\rho=7/10\), (DG36) lies in
\([0.012105848280513554,0.012685028889351772]\).
Thus the right side of (DG29), using \(\Gamma_\sharp\), is positive
for the same balanced profiles and sufficiently small nonzero
\(\epsilon\). Since its first-moment penalty is exact in that
neighborhood, improving a bound that uses only those norm and closure
constraints cannot fix this coefficient. The remainder rearrangement
envelope or additional angular information must also be strengthened.
This is a limitation of the specified upper-bound criterion, not a
counterexample to the original inequality. This angular refinement does
not itself extend the unrestricted cutoff; DG10 supplies a separate
spectral extension to \(11/25\).

### Proposition DG9 (measure-dependent rearrangement with the exact arc interaction)

Let \(B_a(p,q)=\langle1_{I_p},U_a1_{I_q}\rangle_\sigma\), where
the two arcs have the same center. For \(a>0\) and \(0\le s\le2\), put

\[
 c_a(s)=\frac{e^{a\cos(\pi s)}-e^{-a}}
                  {(1+\cos(\pi s))I_0(a)},\qquad
 c_a(1)=\frac{ae^{-a}}{I_0(a)}.
 \tag{DG37}
\]

The value at one is the continuous limit. For arbitrary measurable
\(E,F\subset S^1\) of masses \(p,q\), with first moments \(u,v\),

\[
 \langle1_E,U_a1_F\rangle
 \le B_a(p,q)-c_a(p+q)
       \left[\frac{\sin(\pi p)\sin(\pi q)}{\pi^2}-u\cdot v\right].
 \tag{DG38}
\]

This improves the mass-independent coefficient used in DG6:
\(c_a(s)\ge ae^{-a}/I_0(a)\) for every \(s\).
For three-component profiles define
\(c_a^*(p,q)=\min_i c_a(p_i+q_i)\). For every ordinary balanced
planar three-cell partition, with \(p_t=p(\sqrt{2t})\), one therefore has

\[
\begin{split}
 S_\rho(A)-S_\rho(T)\le\mathcal E_\rho(p):=
 \iint R_0(t,s)\bigg[&\sum_i B_{a(t,s)}(p_{t,i},p_{s,i})
                         -3B_{a(t,s)}(1/3,1/3)\\
 &-c_{a(t,s)}^*(p_t,p_s)\Gamma_\sharp(p_t,p_s)\bigg]d\nu(t)d\nu(s).
\end{split}
 \tag{DG39}
\]

This keeps the exact arc interaction instead of its scalar quadratic
envelope in (DG29), as well as a stronger moment penalty. The bound is
valid for all profiles, including those far from the balanced one.
It does not assume that the actual circular cells are arcs.

Here is the proof of (DG38). First suppose \(p+q\le1\), and set
\(z_0=\cos(\pi(p+q))\). Temporarily omit the positive factor
\(I_0(a)^{-1}\), and let \(c\) be the secant slope of \(e^{az}\)
between \(-1\) and \(z_0\). With
\(R(z)=e^{az}-cz\), convexity gives
\(R(-1)=R(z_0)\), \(R(z)\le R(z_0)\) for \(z\in[-1,z_0]\),
and \(R\) is nondecreasing on \([z_0,1]\). Consequently the function
which equals \(R(z_0)\) below \(z_0\) and \(R(z)\) above it is a
nondecreasing majorant \(\widehat R\) of \(R\). Add a constant if
necessary to make this kernel nonnegative; the added constant cancels
from the comparison because both pairings have masses \(p,q\).
Apply the circular rearrangement theorem used in DG6 to
\(\widehat R(\cos\vartheta)\). On the aligned pair of arcs,
\(|\vartheta|\le\pi(p+q)\), so \(\widehat R=R\) there.
Restoring the separated cosine term proves (DG38). The case
\(p+q=1\) follows either directly or by the continuous secant limit.

If \(p+q>1\), apply the same argument to \(E^c,F^c\). Their first
moments are \(-u,-v\); their sine caps are unchanged. Since \(U_a\)
is Markov,
\(\langle1_E,U_a1_F\rangle=p+q-1+
\langle1_{E^c},U_a1_{F^c}\rangle\), and the identical relation
holds for the aligned arcs. The complement masses sum to \(2-p-q\),
which has the same cosine in (DG37). This proves all cases, including
empty or full sets. Convexity also gives the asserted lower bound on
the secant slope by the derivative at \(-1\).

Each individual cosine deficit in (DG38) is nonnegative by its norm
caps. Their sum is at least \(\Gamma_\sharp\) by DG8. Summation
therefore proves the angular version of (DG39), and polar integration
proves the displayed bound. All its terms are integrable. Indeed
\(0\le\sum_iB_a(p_i,q_i)\le\sum_i\min(p_i,q_i)\le1\), and (DG38)
shows that its total subtracted penalty, including
\(c_a^*\Gamma_\sharp\), lies between zero and this arc sum.
The radial Markov measure is a probability measure.

The arc interaction is itself an explicit one-dimensional integral.
Define the even function
\(\Psi_a(u)=\int_0^{|u|}(|u|-v)k_a(\pi v)dv\). Direct endpoint
integration gives
\(B_a(p,q)=[\Psi_a(p+q)-\Psi_a(p-q)]/2\). In particular, at
\(p=q=1/3\), its second derivatives are

\[
 (B_a)_{pp}=\tfrac12[k_a(2\pi/3)-k_a(0)],\qquad
 (B_a)_{pq}=\tfrac12[k_a(2\pi/3)+k_a(0)].
\]

These identities permit a closed, nonquadrature check of the local
coefficient of (DG39). For any bounded finite exponential sum
\(f(t)=\sum_j d_j e^{-\lambda_jt}\) with \(\lambda_j\ge0\) and
\(\sum_jd_j/(1+\lambda_j)=0\), use the profiles (DG30) with this
\(f\). For \(z\in\{1,-1/2,-1\}\), define

\[
\begin{split}
 H_z(\lambda,\mu)&=d\,\Xi(P_{\lambda\mu},-\rho z),\\
 N_z&=\sum_{j,k}d_jd_k H_z(\lambda_j,\lambda_k),\qquad
 D_z=\sum_{j,k}d_jd_k H_z(\lambda_j+\lambda_k,0).
\end{split}
\]

Here \(d=1-\rho^2\), and \(\Xi,P\) are the exact quadrant formulas
in DG5. Then

\[
 \mathcal E_\rho(p)=2\epsilon^2 E_\rho(f)+O(\epsilon^3),\qquad
 E_\rho(f)=\tfrac12[N_1+2N_{-1/2}-N_{-1}-D_1+D_{-1}].
 \tag{DG40}
\]

For completeness, the exact arc part contributes
\([N_1+N_{-1/2}-D_1+D_{-1/2}]/2\). At the balanced profile
\(R_0c_a(2/3)=2R_0(e^{-a/2}-e^{-a})/I_0(a)\).
Combining this identity with the coefficient \(1/8\) from (DG35)
subtracts
\([D_{-1/2}-D_{-1}-N_{-1/2}+N_{-1}]/2\). The sum is (DG40).
The expansion is uniform over angular shapes. For each fixed
\(\rho<1\), derivatives of the arc kernels and the local secant
coefficient are bounded by a polynomial in \(a\) times
\(e^a/I_0(a)\). The product with \(R_0\) is integrable, since its
radial exponent is a positive definite Gaussian quadratic form.
Thus the bounded profile and the uniform expansion (DG35) justify
the integrated \(O(\epsilon^3)\), even though the minimum in
\(c_a^*\) need not be differentiable.

For the explicit \(f\) of (DG21), the directed checker gives, at
\(\rho=9/10\),
\(E_\rho(f)\in[-0.06164098388637688,-0.061640983886376864]\).
It follows that all partitions with this sufficiently small circular
mass profile strictly lose to the simplex at \(\rho=9/10\), with
arbitrary angular shapes. This is again a restricted-profile result.

The sufficient assertion \(\mathcal E_\rho(p)\le0\) for every
balanced profile and every correlation is false. Set

\[
 g(t)=\frac1{1000}\sum_{j=1}^6 n_j
        \left(e^{-\lambda_jt}-\frac1{1+\lambda_j}\right),\qquad
 (\lambda_j)=(1,4,16,64,256,1024),\quad
 (n_j)=(27,37,74,181,145,1000).
\]

This is exactly mean zero and \(|g|\le183/125\), so the profiles
\(1/3+\epsilon(1,-1,0)g(t)\) are balanced and feasible for small
\(\epsilon\). At \(\rho=99/100\), the same closed formula gives

\[
 E_{99/100}(g)\in
 [1.202983953209141\mathbin{\cdot}10^{-6},
  1.2029839532091414\mathbin{\cdot}10^{-6}].
\]

Both new enclosures use 160-bit Arb evaluation of closed quadrant
expressions; there is no angular or radial numerical quadrature.
The strictly positive value makes (DG39)'s right side positive for
small nonzero \(\epsilon\). It does not determine the actual noise
stability of those partitions. The bound (DG39) is valid but still
loses angular information needed for the full theorem.

The closest inspected quantitative rearrangement input is
[Christ--Iliopoulou](https://preprint.press.jhu.edu/ajm/sites/default/files/AJM-christ-iliop.pdf),
Theorem 1.3, printed p. 4, with Definitions 1.6--1.7.
For fixed \(\eta>0\), the three measures must be \(\eta\)-strictly
admissible and \(\eta\)-bounded, including a lower measure bound and
a total measure at most \(2-\eta\). A sufficiently small convolution
deficit then gives distance at most \(C_\eta\sqrt\delta\) to
compatibly centered parallel rank-one Bohr triples. Its actual statement
and hypotheses were inspected. Converting such stability to a useful
explicit moment-and-remainder penalty, while controlling admissibility
across the radial parameters, is a further obligation. This theorem is
context here and is not used to prove (DG37)--(DG40).

Run the directed witnesses, closed-form coefficients and rational controls from the
repository root with the existing environment:

```sh
.venv/bin/python research/plurality-stability/check_three_cell_coupled_radial.py
```

The verification record pins this checker, its exact inputs and the
current proof source. These new arguments have internal checking only;
independent replay, formal verification and external review remain
unperformed.


## A spectral refinement of the unrestricted direct interval

### Proposition DG10 (global comparison through rho equal to 11/25)

For every $d\ge2$, every ordinary balanced Gaussian three-cell partition,
and every $0<\rho\le11/25$, the inequalities in (DG14) hold.
The new interval $[21/50,11/25]$ is global in partition distance.
It uses (DG16) and a directed spectral estimate rather than continuation.
The all-positive-correlation target remains open. A separate
[near-one analysis](THREE_CELL_NEAR_ONE_EXACT.md) records a qualitative
interval adjacent to one; the present certificate addresses the lower interval.

The distinction from DG4 is that a Hilbert--Schmidt norm exceeding two
need not mean that the largest eigenvalue exceeds two. Define the compact
self-adjoint operator on $L^2(\nu)$ with kernel

$$
 A_\rho(t,s)=\frac{K_+(t,s)-\ell_\rho(\sqrt t)-\ell_\rho(\sqrt s)}
                    {\sqrt{J_+(t)J_+(s)}}.
 \tag{DG41}
$$

Its squared Hilbert--Schmidt norm is exactly (DG17). For any bounded
mean-zero $z$, additive centering disappears, so
$\langle z,\mathcal K_+z\rangle\le\lambda_{\max}(A_\rho)\int J_+z^2d\nu$.
It therefore suffices to bound the largest eigenvalue, using the rest of
the spectrum to control the error in a Rayleigh quotient.

Here is the spectral estimate, with its isolation hypothesis explicit.
Suppose $\|A\|_{HS}^2\le H^2$ and a unit vector $\phi$ satisfies
$a_-\le a=\langle\phi,A\phi\rangle\le a_+$ and
$\|(A-2)\phi\|^2\le r^2$. Put $c=\sqrt{H^2-a_-^2}$.
If $a_->c\ge0$, then

$$
 \lambda_{\max}(A)\le a_++\frac{r^2}{a_--c}.
 \tag{DG42}
$$

Indeed the largest eigenvalue $\lambda_1$ is at least $a_->0$;
removing its squared contribution from the Hilbert--Schmidt sum bounds
every other eigenvalue in absolute value by $c$. In particular it is
simple. For the spectral probability measure of $\phi$,
$(\lambda_1-\lambda)(\lambda-c)\le0$ at every spectral point.
Taking expectations gives
$(\lambda_1-a)(a-c)\le\|(A-a)\phi\|^2$.
The right side is $\|(A-2)\phi\|^2-(a-2)^2\le r^2$, proving (DG42).
This is the usual one-sided Temple argument, derived here.
[Kuroda, Theorem A](https://www.kurims.kyoto-u.ac.jp/~prims/pdf/43-2/43-2-24.pdf),
printed p. 506, was inspected: its isolated-eigenvalue interval and
residual-gap hypotheses are essential. No spectral isolation is inferred
from a floating eigenvalue calculation; it is supplied by $H^2$ and $a_-$.

The exact auxiliary trial is

$$
 f(t)=t^{-1/2}\sum_{j=1}^5 b_j e^{-\lambda_jt}-M,\qquad
 M=\sqrt\pi\sum_j\frac{b_j}{\sqrt{1+\lambda_j}},
$$

where

$$
 (\lambda_j)=(0,1/4,1,4,16),\qquad
 (b_j)=(151/250,-1,28/125,-33/1250,121/2500).
 \tag{DG43}
$$

Thus $\int f\,d\nu=0$ exactly. Although $f$ is unbounded at zero,
$\sqrt{J_+}f\in L^2(\nu)$: $J_+(x^2)/x$ has a finite positive limit
at zero, and $f(x^2)=O(1/x)$. At infinity $f$ is bounded and $J_+$ is
bounded. This trial is used only as a Hilbert-space vector; the partition
mass functions to which (DG16) is applied remain bounded. Set
$D=\int J_+f^2d\nu$ and $\phi=\sqrt{J_+}f/\sqrt D$.

All trial integrals are one dimensional after an exact Gaussian radial
integration. The following formulas specify the reduction checked by the
new program. Write $d=1-\rho^2$, $\kappa=\rho^2/d$,
$q=\rho x/\sqrt d$, $w_j=1+\lambda_jd$, and
$h(x)=\sum_jb_je^{-\lambda_jx^2}-Mx$. For real $q\ge0$ put

$$
\begin{split}
 L(q)&=\tfrac12\int_0^{2/3}(2+3v)e^{q^2\cos^2(\pi v)}
                          \operatorname{erfc}(-q\cos(\pi v))dv,\\
 P(q)&=\tfrac12\int_0^{2/3}(2+3v)
    [1+\sqrt\pi q\cos(\pi v)e^{q^2\cos^2(\pi v)}
                          \operatorname{erfc}(-q\cos(\pi v))]dv,\\
 R(q)&=1-\sqrt\pi q e^{q^2}\operatorname{erfc}(q),\\
 S(q)&=\sqrt\pi(1+2q^2)e^{q^2}\operatorname{erfc}(q)-2q.
\end{split}
$$

The displayed expression for $S$ is also used at negative arguments.
Let $H_+(q)$ be the positive angular integral `positive_q_integrand`
in the retained upstream source. Explicitly it is
$\frac\pi2\int_0^{2/3}v(1-3v/2)\sin(\pi v)S(-q\cos(\pi v))dv$.
Then

$$
\begin{split}
 \bar J(x)&=e^{\kappa x^2}J_+(x^2)/x
      =\frac\rho{\sqrt d}[H_+(q)+(j_1/2)S(q)],\\
 G(x)&=e^{\kappa x^2}(\mathcal K_+f)(x^2)\\
 &=\sum_jb_j\left[\frac{\sqrt\pi}{\sqrt{dw_j}}L(q/\sqrt{w_j})
           -\frac{4j_1q}{\sqrt d\,w_j}R(q/\sqrt{w_j})\right]
       -M[P(q)-2j_1qS(q)].
\end{split}
 \tag{DG44}
$$

For example, substitute $y=\sqrt{d/w_j}\,u$ in the radial integral
with $s=y^2$. The elementary integrals of $1,u,u^2$ against
$e^{-u^2+2bu}$ give $L,P,R,S$ and the displayed factors.
The additive center is integrated exactly by gamma moments:

$$
 c_0=\int\ell_\rho(\sqrt t)f(t)d\nu(t)
 =\sum_{n=0}^4e_n(\rho)\left[
   \sum_j\frac{b_j\Gamma((n+1)/2)}{(1+\lambda_j)^{(n+1)/2}}
                         -M\Gamma(1+n/2)\right],
$$

where $\ell_\rho(x)=\sum_ne_n(\rho)x^n$ is the original binary64
polynomial, treated exactly as in DG4. With
$U(x)=G(x)-c_0e^{\kappa x^2}-2\bar J(x)h(x)$, the quantities used are

$$
\begin{split}
 D&=\int_0^\infty2e^{-(1+\kappa)x^2}\bar J(x)h(x)^2dx,\\
 I&=\int_0^\infty2e^{-(1+\kappa)x^2}h(x)U(x)dx,\\
 B&=\int_0^\infty\frac{2e^{-(1+\kappa)x^2}}{\bar J(x)}U(x)^2dx,\\
 a&=2+I/D,\qquad \|(A_\rho-2)\phi\|^2=B/D.
\end{split}
 \tag{DG45}
$$

The middle identity uses the exact mean zero of $f$. These integrands
extend continuously to $x=0$, so no singular point is omitted.

The finite domain is $[0,7]$. The functions $L,P,H_+,S,R$ are enclosed
at exact binary64 grid nodes by 100-bit Arb arithmetic; the angular
integrals use validated `acb.integral` on entire functions, so ignoring
its analytic flag introduces no branch assumption. Between nodes the
program uses the linear interpolant with the Lagrange error
$\frac12\sup|F''|\,|(q-q_0)(q-q_1)|$. This remains valid when a query
interval crosses a cell boundary, by taking the curvature bound on its
whole convex hull with the two interpolation nodes. All curvature
majorants are explicit increasing functions of $q\ge0$. They follow
from differentiating $e^{z^2}\operatorname{erfc}(-z)$ and, for example,

$$
 S''(z)=\sqrt\pi(6+24z^2+8z^4)e^{z^2}\operatorname{erfc}(z)
                         -20z-8z^3.
$$

The bound for $H_+''$ uses
$\frac\pi2\int_0^{2/3}v(1-3v/2)dv=\pi/27$.
The remaining one-dimensional integrations are complete directed interval
Riemann sums, not point sampling. Positivity used to clamp lower bounds
follows from the positive Gaussian integral representations of $L,P,R,S$
and $H_+$. The lower bound $\bar J\ge c_+\rho/\sqrt d$ is also kept.

Here are explicit tail bounds, including their uniform parameter domain.
Let $X=7$, $C=\sum_j|b_j|$, $F=C/X+|M|$,
$j_{\max}=1/2+j_1/(2e)$ and $k_{\max}=2+2j_1/e$.
The angular definitions give $0<J_+\le j_{\max}$ and
$|K_+|\le k_{\max}R_0$. Moreover

$$
 \sup_t P_\rho(t^{-1/2})\le\sqrt{\pi/d}.
$$

To verify this bound, integrate the identity
$t^{-1/2}=\pi^{-1/2}\int_0^\infty u^{-1/2}e^{-ut}du$ against the
conditional Gaussian radial distribution. Its Laplace transform is
$(1+du)^{-1}\exp[-\rho^2tu/(1+du)]\le(1+du)^{-1}$; the last integral
is $\sqrt{\pi/d}$. Therefore let
$K_*=k_{\max}(C\sqrt{\pi/d}+|M|)$ and
$U_*=K_*+|c_0|+2j_{\max}F$. The tails of $D,I,B$ are bounded by

$$
 0\le D_{\rm tail}\le j_{\max}F^2e^{-X^2},\qquad
 |I_{\rm tail}|\le FU_*e^{-X^2},\qquad
 0\le B_{\rm tail}\le
 \frac{U_*^2e^{-(1-\kappa)X^2}}
      {c_+(\rho/\sqrt d)(1-\kappa)X}.
 \tag{DG46}
$$

All rates are positive on $[21/50,11/25]$. The coarse two-dimensional
Hilbert--Schmidt cover uses the original (DG19) tail argument, still
valid because $11/25<1/2$. At the exact upper endpoint the squared-kernel
completion rate is $47/168$, and $1-\kappa=383/504$; the outward-expanded
inputs retain positive rates. No table from a smaller parameter range is
reused without enclosing the new domain.

The complete run accepted 397 covering correlation subintervals after
774 attempts, using 20 coarse Hilbert--Schmidt parent bounds. Its largest
outward-rounded spectral upper bound was `1.999981001453614`, so uniformly

$$
                 \lambda_{\max}(A_\rho)<1.999981002<2.
                 \tag{DG47}
$$

The run used 262144 radial intervals for the one-dimensional trial integrals,
4096 interpolation cells, 800 intervals per coordinate for each coarse
Hilbert--Schmidt square, 1024 upstream q bins, 4096 angular-parameter bins,
131072 exponential bins and 64 upstream angular midpoint intervals.
It completed in 524.30 seconds. Every accepted child lies in its
certified parent, and an exact Fraction audit of the retained binary64 bounds
checks the complete rational endpoint cover, all spectral caps and all
397 assemblies of (DG42). Positive and rejected controls are recorded;
these finite controls supplement the written unrestricted spectral argument.


Substitute the accepted bounds in (DG42), then (DG16), to prove the
planar comparison on $[21/50,11/25]$. DG4 supplies the lower correlations;
Heilman--Tarter maximizer existence and dimension reduction, with exactly
the hypotheses checked in DG4, extend this bound to every $d\ge2$.
The exact quadratic expansion gives the requested inequality.

The [spectral checker](sources/cut_bounds_6b720fe57bfe/check_positive_spectral_extension.py)
imports both retained upstream modules unchanged. Its
[verification record](positive_spectral_extension_verification.json)
pins the program, sources, exact inputs, environment and full covering
output. Replay from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/sources/cut_bounds_6b720fe57bfe/check_positive_spectral_extension.py
```

This is an internal computer-assisted result. It is not independent replay,
formal verification, external review or a proof for all positive correlations.

## Adopting the joint angular and near-one comparisons

### Corollary DG11 (annular localization of any unresolved global failure)

For every compact correlation interval \(J\subset(0,1)\), there are
\(0<a_J<R_J<\infty\), \(\delta_J>0\) and \(c_J>0\) such that
every ordinary balanced planar partition satisfying

\[
 E_J(A):=\int_{a_J}^{R_J}
       \sum_i|p_i(r)-1/3|^2d\mu(r)\le\delta_J
\]

obeys

\[
 S_\rho(T)-S_\rho(A)\ge c_J\left[
   \int_0^{R_J}r\sum_i|p_i(r)-1/3|^2d\mu(r)
   +\int_{R_J}^\infty\sum_i|p_i(r)-1/3|^2d\mu(r)\right],
 \qquad \rho\in J.
 \tag{DG48}
\]

Equality in the simplex comparison forces a common rotated, labeled
simplex up to null sets. No angular closeness, profile derivative bound,
or bound on the circular masses inside or outside the annulus is required.
This is the annular mass theorem (JA1), including its final squared-mean
extension, in the [joint angular note](THREE_CELL_JOINT_ANGULAR_ENVELOPE.md).
It is adopted here as an analytic input; its proof remains in that single
source.

The useful change of mechanism is its exact joint pair Hessian. Write
\(b=p-\mathbf1/3\), \(c=q-\mathbf1/3\),
\(k_0=k_\kappa(0)\) and \(k_\ell=k_\kappa(2\pi/3)\).
For bounded angular parameter and masses near \(\mathbf1/3\), the
maximizing two-slice pair consists of compatible consecutive arcs, and
its quadratic change is

\[
 \tfrac13[-(k_0-k_\ell)(|b|^2+|c|^2)
                       +(2k_0+k_\ell)b\cdot c].
\]

After polar integration this is (JA6), with standard facet speed
\(f(r)=r b(r)\) and edge measure
\(d\eta(r)=(2\pi)^{-1}e^{-r^2/2}dr\). Its interface kernel is
\(\mathsf K_\rho^{\rm st}=2M_1+M_{-1/2}\) of (S14); the radial
kernel \(K_+\) keeps its definition in (DG15). The interface gap (S37)
controls this joint Hessian on every compact positive correlation interval.
The radial proof of (JA1) handles arbitrary inner and outer masses by
small-parameter angular coercivity and Gaussian kernel tails. Its
squared-mean extension handles the exceptional annular radii by the
all-mass centered angular penalty in Heilman's Lemma 6.1. Thus it retains
the label compatibility lost in the sufficient criteria of DG5--DG9.
The constants and the admissibility neighborhood are allowed to depend
on \(J\); no uniform limit as \(\rho\uparrow1\) is asserted in (DG48).

There is also one annulus and one positive threshold sufficient to
localize every possible failure of the original all-correlation inequality.
More precisely, there exist \(u\in(11/25,1)\),
\(0<a<R<\infty\) and \(\delta>0\) with the following property:
if the inequality fails at any correlation in any \(d\ge2\), then

\[
 11/25<\rho<u,
 \quad\text{and a planar balanced global maximizer }B\text{ exists with}
 \quad\int_a^R\sum_i|p_{B,i}(r)-1/3|^2d\mu(r)>\delta.
 \tag{DG49}
\]

Every planar strict improver at that correlation has the same annular
lower bound. Consequently its Gaussian symmetric-difference distance
from the centered simplex orbit, summed over labels, is greater than
\(\delta\).

To prove this assertion, let \(\rho_*<1\) be the qualitative threshold
in the [exact near-one theorem](THREE_CELL_NEAR_ONE_EXACT.md). Set
\(u=(1+\max(11/25,\rho_*))/2\), and apply (DG48) on
\(J=[11/25,u]\). DG10 excludes every strict improvement below or at
\(11/25\), and the near-one theorem excludes it above \(\rho_*\).
Thus every strict failure lies in this fixed \(J\). The contrapositive
of (DG48) gives the annular lower bound for every planar strict improver.
For an improvement in higher dimension, maximizer existence and
Heilman--Tarter dimension reduction supply a planar global maximizer with
at least that value. These are exactly the existence and dimension
hypotheses used in DG4 and DG10. Finally, for any centered simplex \(U\),
conditional Jensen gives

\[
 \int_a^R\sum_i|p_{B,i}(r)-1/3|^2d\mu(r)
 \le\sum_i\gamma_2(B_i\mathbin\triangle U_i).
\]

Take the infimum over its rotations and labelings to get the distance
assertion. The same annular exclusion applies to a nonsimplex tied
global maximizer whenever its correlation belongs to \(J\), using the
equality statement in (DG48).

In particular, fix any bounded radial vector \(v\) with
\(\sum_i v_i=0\) pointwise and \(\int v_i d\mu=0\). Uniformly for
\(\rho\in J\), all sufficiently small feasible profiles
\(p=\mathbf1/3+\epsilon v\) satisfy the simplex comparison, with
arbitrary angular shapes. Indeed
\(E_J(A)=\epsilon^2\int_{a_J}^{R_J}|v|^2d\mu\), so (DG48) applies.
This includes both the profile (DG21) and the six-exponential profile in
DG9 at every fixed positive correlation. Their positive separated
upper-bound coefficients remain valid limitations of those particular
criteria; they are not actual stability counterexamples.

The near-one core and exterior estimates enter (DG49) through their
proved exact comparison theorem. Their small-noise remainders and
qualitative thresholds do not provide numerical constants extending that
theorem down to \(11/25\). The outstanding direct task is the comparison
for annular mass profiles separated from the balanced one as in (DG49).
No exclusion of that entire class, or certified overlap of the correlation
intervals, follows from this localization.

## Exact comparisons within reflection orbits

### Proposition DG12 (reflection cuts and the limit of type-constant sorting)

Fix \(0<\rho<1\), a reflection \(R\) through a hyperplane containing
the origin, and one of its open half-spaces \(H\). A reflection-pair
swap exchanges the labels of \(x\) and \(Rx\), for an arbitrary
measurable set of pairs. Such swaps preserve every cell's circular mass
at every radius where that mass is defined, and hence preserve Gaussian
balance. They have an exact quadratic comparison and an exact cut
criterion, given below. If at most two unordered label-pair types occur,
a common orientation increases stability, strictly unless the partition
already has that orientation or its reverse almost everywhere on the
active pairs.

With all three pair types present this sorting statement is false.
For every fixed \(\rho\in(0,1)\), there is an ordinary balanced planar
partition and a reflection for which every choice of a single orientation
per unordered pair type has strictly smaller stability. This is a
counterexample to that sorting operation, not to the simplex comparison.

For context, scalar polarization places the maximum of two reflected
values on one side and the minimum on the other; see
[Bianchi--Gardner--Gronchi--Kiderlen, equation (1), printed p. 2](https://arxiv.org/pdf/1912.07498).
That primary definition concerns a scalar function. No three-label
monotonicity theorem is imported here; all the comparisons below follow
directly from the Gaussian kernel.

Write \(F=(1_{A_1},1_{A_2},1_{A_3})\), and let the Gaussian-relative
Mehler kernel in dimension \(d\) be

\[
 M_\rho(x,y)=(1-\rho^2)^{-d/2}
 \exp\left[-\frac{\rho^2(|x|^2+|y|^2)-2\rho x\cdot y}
                  {2(1-\rho^2)}\right].
\]

For \(x,y\in H\), put
\(W(x,y)=M_\rho(x,y)-M_\rho(x,Ry)\). This is strictly positive:
if \(u,v>0\) are the normal coordinates, then
\(M_\rho(x,Ry)=M_\rho(x,y)e^{-2\rho uv/(1-\rho^2)}\).
On \(H\), define

\[
 m(x)=\tfrac12(F(x)+F(Rx)),\qquad
 E=\{x:F(x)\ne F(Rx)\}.
\]

The unordered pair at \(x\in E\) is denoted \(e(x)\in\{12,23,31\}\).
Choose incidence vectors
\(d_{12}=(1,-1,0)\), \(d_{23}=(0,1,-1)\),
\(d_{31}=(-1,0,1)\). Every partition with these same reflection-pair
multisets is specified by a measurable sign \(\sigma:E\to\{-1,1\}\),
through \(F(x)-F(Rx)=\sigma(x)d_{e(x)}\). Its stability is exactly

\[
 S_\rho(A_\sigma)=C_R(m)+\frac12
 \int_{E\times E}W(x,y)
   (d_{e(x)}\cdot d_{e(y)})\sigma(x)\sigma(y)
                  d\gamma_d(x)d\gamma_d(y),
 \tag{DG50}
\]

where
\(C_R(m)=2\int_{H\times H}(M_\rho(x,y)+M_\rho(x,Ry))
m(x)\cdot m(y)d\gamma_d(x)d\gamma_d(y)\).
Indeed, putting \(\delta_x=F(x)-F(Rx)\), the two sums of label
inner products are
\[
 \begin{split}
 F(x)\cdot F(y)+F(Rx)\cdot F(Ry)
     &=2m(x)\cdot m(y)+\tfrac12\delta_x\cdot\delta_y,\\
 F(x)\cdot F(Ry)+F(Rx)\cdot F(y)
     &=2m(x)\cdot m(y)-\tfrac12\delta_x\cdot\delta_y.
 \end{split}
\]
Splitting the full Gaussian integral into its four half-space pairs
proves (DG50). All integrals are finite since \(0<W\le M_\rho\)
on \(H\times H\) and the label vectors are bounded.

If \(\sigma'\) reverses \(\sigma\) on \(B\subset E\), then

\[
 S_\rho(A_{\sigma'})-S_\rho(A_\sigma)
 =-2\int_{B\times(E\setminus B)}W(x,y)
   (d_{e(x)}\cdot d_{e(y)})\sigma(x)\sigma(y)
                         d\gamma_d(x)d\gamma_d(y).
 \tag{DG51}
\]

Thus nonnegativity of every integral on the right before the factor
\(-2\) is necessary and sufficient for maximality within this fixed
reflection-orbit class: every other sign function differs on one such
set \(B\). Every balanced global maximizer satisfies these inequalities
for every reflection. This is a finite-variation condition, requiring no
interface regularity or closeness assumption. Preservation of circular
masses follows by pairing points under the radius-preserving reflection
on each sphere; it does not require an integration in the radius.

If only two pair types occur, orient both incidence vectors away from
their common label. Their inner products are then \(2\) for equal
types and \(1\) for different types. The kernel multiplying the signs
in (DG50) is strictly positive. Consequently constant signs maximize it;
the gain over any other sign choice is
\(2\int_{\{\sigma=1\}\times\{\sigma=-1\}}W(d_e\cdot d_f)
d\gamma_d d\gamma_d>0\) when both sets have positive measure.
The same conclusion holds with one type. With all three cyclic incidence
vectors no choice of three fixed sign gauges makes every inner product
positive: their three off-diagonal sign requirements have product \(-1\).
The following construction proves that this sign obstruction is realized
by ordinary Gaussian partitions, even after exact balance is imposed.

Work in \(\mathbb R^2\), with \(H=\{(t,v):v>0\}\) and
\(R(t,v)=(t,-v)\). Put \(h=1-\rho^2\) and
\[
 a=\sqrt{\frac{h\log2}{2\rho}},\qquad
 x_j=((2j-5)a,1),\quad 1\le j\le4,\qquad
 w_j=W(x_j,x_j)^{-1/2}.
\]
The exact normalized interaction is
\[
 w_iw_jW(x_i,x_j)
 =\exp\left[-\frac{\rho((2i-5)a-(2j-5)a)^2}{2h}\right]
 =2^{-(i-j)^2}.
 \tag{DG52}
\]
Assign the four patches types \((12,23,31,12)\). The unordered atomic
orientation energy is
\[
 \mathcal E(s)=\sum_{i<j}2^{-(i-j)^2}
                 (d_{e_i}\cdot d_{e_j})s_i s_j.
\]
Exhausting its sixteen signs gives
\[
 \mathcal E(1,-1,1,-1)=\frac{351}{256},\qquad
 \max_{s_1=s_4}\mathcal E(s)=\frac{161}{256},\qquad
 \sum_{i<j}2^{-(i-j)^2}|d_{e_i}\cdot d_{e_j}|=\frac{417}{256}.
 \tag{DG53}
\]
The restriction \(s_1=s_4\) is exactly a single orientation for type
\(12\); the other two types each occur once and their signs are free.
The strict gap is therefore \(95/128\).

Here are explicit choices transferring this gap to positive-area patches.
Let
\[
 R_0=\sqrt{1+9a^2},\qquad
 L=2\left[\frac{\rho(R_0+1)}{1-\rho}+2\right],\qquad
 b=\min\{1/4,a/4,1/(100L)\}.
\]
For points in the \(b\)-balls about any two centers, the normal
coordinates exceed \(1/2\) and the norms are at most \(R_0+1\).
Writing \(\beta=2\rho/h\), direct differentiation gives
\[
 \nabla_x\log W(x,y)
 =\frac{-\rho^2x+\rho y}{h}
       +\frac{\beta v}{e^{\beta uv}-1}(0,1).
\]
The second term has norm at most \(1/u\le2\), and the first at most
\(\rho(R_0+1)/(1-\rho)\). The same bounds hold for \(\nabla_y\).
Integrating along the two segments from the centers gives
\[
 |\log W(x,y)-\log W(x_i,x_j)|\le Lb\le1/100.
\]
Hence the relative kernel error is at most \(\eta=1/50\), since
\(e^{1/100}\le100/99<1+\eta\) and
\(e^{-1/100}\ge99/100>1-\eta\).

Set \(C=\sum_jw_j\) and choose the positive number
\[
 \epsilon=\frac12\min\left\{\frac1{12C},
       \min_j\frac{\gamma_2(B(x_j,b))}{w_j}\right\}.
\]
There is a unique radius less than \(b\) giving a ball \(B_j\)
centered at \(x_j\) of Gaussian mass \(m_j=\epsilon w_j\).
The four balls are disjoint and lie in \(H\). Their total mass is
at most \(1/24\). On \(B_j\cup RB_j\), assign the two labels of
type \(e_j\), one to each reflected ball, with orientation \(s_j\).
On the remaining half-space region \(G=H\setminus\bigcup_jB_j\),
choose three measurable sets \(C_k\) of masses
\[
 \gamma_2(C_k)=\frac12\left(\frac13-
                         \sum_{j:k\in e_j}m_j\right).
\]
These masses are positive and sum to \(\gamma_2(G)\). For a concrete
choice take successive quantiles of the tangent coordinate \(t\) under
Gaussian measure restricted to \(G\). This measure has no atoms on
coordinate levels, so the specified masses are attained. Label both
\(C_k\) and \(RC_k\) by \(k\). Each full cell then has mass exactly
\(1/3\), independently of the four signs. Null boundaries may be
assigned arbitrarily.

All sixteen resulting partitions have the same reflection-pair multisets
and the same circular masses. The normalized cross-patch couplings
\(\epsilon^{-2}\int_{B_i\times B_j}W\,d\gamma_2d\gamma_2\)
differ from (DG52) by at most the relative error \(\eta\).
The within-patch terms and \(C_R(m)\) are independent of the signs.
It follows from (DG50)--(DG53) that the split signs \((1,-1,1,-1)\)
beat every type-constant choice by at least
\[
 \epsilon^2\left[\frac{95}{128}
                 -2\frac1{50}\frac{417}{256}\right]
       =\frac{4333}{6400}\epsilon^2>0.
 \tag{DG54}
\]
This proves the stated ordinary balanced counterexample at every fixed
positive correlation. It excludes only a universal type-constant sorting
step; the exact all-cut condition (DG51) remains valid with all three
types. It neither supplies nor excludes a global simplex improver.

The [exact finite checker](check_three_cell_reflection_orientation.py)
exhausts the 81 four-label identities, sixteen orientation patterns,
eight type-constant choices and 256 sign/cut identities, and checks the
rational robustness margin. The Gaussian normalization, explicit patch
construction and unrestricted statements use the written proof above.
Its [verification record](three_cell_reflection_orientation_verification.json)
records the internal replay and its limits. Run from the repository root:

```sh
python3 research/plurality-stability/check_three_cell_reflection_orientation.py
```

## Reflection order and a necessary three-type witness

### Theorem DG13 (reflection-triangle rigidity for actual maximizers)

Fix \(0<\rho<1\). An ordinary balanced planar global maximizer is a
centered regular simplex if and only if, for every reflection through the
origin, at least one of the three unordered reflected label-pair types
has zero Gaussian measure. Thus every nonsimplex global maximizer has
some reflection with all three types present on positive-measure sets.
The four-patch partition in DG12 is not asserted to be a global maximizer.

First there is a quantitative three-type defect valid for any partition.
Let \(P_H A_i\) be scalar polarization of the individual cell toward
\(H\); the three resulting sets need not form a partition. Define
\(\mathcal P_R(A)=\sum_i S_\rho(P_H A_i)-S_\rho(A)\), where the
stability of a single set means its indicator quadratic form. With
\(\delta_i(x)=1_{A_i}(x)-1_{A_i}(Rx)\) on \(H\),

\[
 \begin{split}
 \mathcal P_R(A)&=\frac12\int_{H^2}W(x,y)
  [\sum_i|\delta_i(x)||\delta_i(y)|-\delta(x)\cdot\delta(y)]
                                     d\gamma_2(x)d\gamma_2(y)\\
 &=2\sum_i\int_{\{\delta_i=1\}\times\{\delta_i=-1\}}
                                     W(x,y)d\gamma_2(x)d\gamma_2(y)\ge0.
 \end{split}
 \tag{DG55}
\]

In three labels the absolute product sum is
\(|d_{e(x)}\cdot d_{e(y)}|\) on active pairs. Thus this is the
defect between (DG50)'s absolute-kernel envelope and the actual energy.
If \(K\subset H\), \(W\ge w>0\) on \(K^2\), and
\(m_e=\gamma_2(K\cap E_e)\), then

\[
 \mathcal P_R(A)\ge2w\min\{m_{12}m_{23},m_{23}m_{31},m_{31}m_{12}\}.
 \tag{DG56}
\]

For \(K=\{|x|\le B,x\cdot n\ge s\}\), \(s>0\), an explicit
choice is
\(w=(1-\rho^2)^{-1}e^{-\rho B^2/(1-\rho)}
(1-e^{-2\rho s^2/(1-\rho^2)})>0\), directly from the Mehler formula.
This bound does not compare the polarized-set sum with the simplex value.

**The integrated defect.** Replacing scalar \(\delta_i\) by
\(|\delta_i|\) in (DG50) proves (DG55). Since \(W>0\), its
vanishing is equivalent to each \(\delta_i\) having one sign almost
everywhere on \(H\). For the cyclic incidences of DG12, the cross-type
defect on \(K\) is at least
\(w\sum_{e<f}\int_{(E_e\cap K)\times(E_f\cap K)}
(1+\sigma(x)\sigma(y))d\gamma_2d\gamma_2\); same-type terms are
nonnegative. For every triple of signs,
\((1+s_1s_2)+(1+s_2s_3)+(1+s_3s_1)\ge2\).
Integrating over the three type sets gives a weighted sum at least
\(2m_{12}m_{23}m_{31}\). Its coefficients are at most
\(M=\max_e m_e\), so the unweighted cross-type sum is at least
\(2m_{12}m_{23}m_{31}/M\). This is (DG56). If a mass is zero,
the asserted lower bound is zero and no division is required.

**Reflection order gives fixed-axis caps.** Suppose now that \(A\)
is a balanced global maximizer and every reflection has at most two
active types. DG12's strict two-type comparison gives
\(\mathcal P_R(A)=0\) for every \(R\). Thus each cell is ordered
with its reflection on either half-space. Put
\(b_i=\int x1_{A_i}d\gamma_2\). The
[full active-rank result](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)
gives \(\operatorname{span}\{b_i\}=\mathbb R^2\). Since
\(\sum_i b_i=0\), all three vectors are nonzero and no pair is
parallel. Set \(c_i=b_i/|b_i|\). For every unit normal \(n\),

\[
 b_i\cdot n=\int_{H_n}(x\cdot n)
                  [1_{A_i}(x)-1_{A_i}(R_nx)]d\gamma_2(x).
 \tag{DG57}
\]

The bracket has a fixed sign, so its order is determined by
\(c_i\cdot n\); it vanishes if that dot product is zero. Hence,
for almost every radius, the section of \(A_i\) is an arc centered
at \(c_i\), possibly empty or full. To check the measurable step,
reflection with normal orthogonal to \(c_i\) gives axis symmetry.
Rotate this axis to angle zero. For \(0<\theta_1<\theta_2<\pi\),
reflection with axis \(\alpha=(\theta_1+\theta_2)/2\) sends the
first point to the second. Its normal \((\sin\alpha,-\cos\alpha)\)
has positive dot product with both \(c_i\) and the first point.
Fubini in \((r,\alpha,\theta_1)\), followed by this nonsingular
angular change of variables, gives the decreasing binary order for
almost every \((r,\theta_1,\theta_2)\). A binary function with
this pairwise order is an initial interval up to null sets; axis
symmetry gives the claimed arc.

This is the familiar foliated Schwarz order characterization. The
measurable Lebesgue \(L^p\) statement in
[Brock--Chiacchio--Croce--Mercaldo, Theorem 2.3](https://www.degruyterbrill.com/document/doi/10.1515/anona-2025-0103/html)
was inspected. A Gaussian cell need not be Lebesgue \(L^p\) on the
whole plane, so that theorem is context; (DG57) and the Fubini argument
prove the weighted measurable version needed here.

**Maximizing scores exclude radial changes.** Three such disjoint arcs
partitioning a circle either have one full and two empty members or
have all three nonempty and proper. Exactly two nonempty proper arcs
would be complementary and have antipodal centers, contrary to the
nonparallel \(c_i\). In the three-arc case their order is the cyclic
order of the fixed centers. Their widths satisfy
\(\omega_i+\omega_{i+1}=2g_i\), where \(g_i\) are the consecutive
center gaps and indices are cyclic. This system has a unique solution,
independent of radius. Write \(I_i\) for these fixed open arcs.
Three-arc radii have positive measure, since otherwise every cell is
radial and all \(b_i=0\).

The first circular moment of \(1_{A_i}\) is \(a_i(r)c_i\), with
\(a_i(r)\ge0\) and strict positivity on the three-arc radii.
For \(U_i=T_\rho1_{A_i}\), positive radial Gaussian factors and
\(I_1(\kappa)>0\) give

\[
 \int_{S^1}\theta U_i(r\theta)d\sigma(\theta)
                  =q_i(r)c_i,\qquad q_i(r)>0\quad(r>0).
 \tag{DG58}
\]

Thus \(U_i-U_k\) cannot be constant on an entire positive-radius
circle, since its first circular moment is
\(q_i(r)c_i-q_k(r)c_k\ne0\).

The actual maximizing-score rule assigns \(A_i\) by the largest
\(U_i-\lambda_i\), up to null sets, with constant prices and
real-analytic Gaussian scores, as in the
[conical rigidity argument](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-cells-and-conical-partitions).
For completeness, the exact quadratic expansion at a global maximum
shows that \(A\) maximizes the linear assignment with scores \(U_i\).
Minimize \(\frac13\sum_i\lambda_i+
\int\max_i(U_i-\lambda_i)d\gamma_2\) in the gauge
\(\sum_i\lambda_i=0\). This convex functional is coercive, since
\(\int U_i=1/3\) bounds it below by \(1/3-\min_i\lambda_i\).
Every adjusted pair difference is nonconstant analytic: its first
Gaussian moment is \(\rho(b_i-b_k)\ne0\). Ties are null, so
differentiation of the price functional makes every winning mass
\(1/3\). The winning assignment is unique almost everywhere, and
\(A\), having the same optimal linear value, equals it.
On each \(I_i\), the label depends only on radius: it is \(i\) at
three-arc radii and the full-circle label otherwise. For almost every
radius its winning inequalities hold at every angle in \(I_i\), by
continuity from the almost-everywhere version. If two radial labels had
positive measure, their essential supports, together with those of the
remaining label, form a finite closed cover of \((0,\infty)\).
Connectedness forces two nonempty supports to meet at a positive radius.
Approach it through good radii of both labels. Their adjusted scores
then agree on all of \(I_i\). Angular analyticity extends equality
to the whole circle, contradicting (DG58). The radial label on each
\(I_i\) is therefore constant, and is \(i\) because three-arc radii
have positive measure. Thus the partition consists of three fixed
sectors. Balance forces their widths to be \(2\pi/3\).

Conversely each simplex cell is the fixed cap
\(c_i\cdot x/|x|>1/2\), so its reflection difference has one sign
on every half-space. Formula (DG55) vanishes. If three types were
present, exhaustion by bounded sets separated from the reflecting line
would capture positive mass of each, contradicting (DG56). This proves
the stated equivalence.

**Uniform witnesses for the remaining failure class.** Together with
DG11, this gives finitely many fixed planar reflections
\(R_1,\ldots,R_N\) and \(\eta>0\) such that every hypothetical
strict failure has a planar balanced maximizing representative \(B\)
at \(11/25<\rho<u\) for which some \(j\) satisfies

\[
 \gamma_2(B_i\cap R_jB_k)\ge\eta\qquad(1\le i<k\le3).
 \tag{DG59}
\]

Here is the compactness proof, including its relaxation step. The
relaxed balanced set \(0\le f_i\le1\), \(\sum_i f_i=1\) is
weakly compact in \(L^2\). At positive correlation the OU operator
is compact and its form is strictly positive on nonzero functions.
A relaxed maximizer is ordinary: where two coordinates stay a fixed
distance from zero and one, two disjoint equal-mass pieces admit a
nonzero mean-zero transfer \(h\), with both \(f+h\) and \(f-h\)
feasible. Their average energy exceeds that of \(f\) by the positive
form of \(h\). If no such set exists, the vector is ordinary almost
everywhere. This also proves equality of relaxed and ordinary maxima.
On \([11/25,u]\), the operators are norm continuous, since
\(\sup_{m\ge0}|\rho^m-s^m|\to0\) uniformly away from one.
Energies on the relaxed set change by at most this operator-norm
difference, so the maximum value is continuous as well.
A weak limit of maximizers is therefore a relaxed maximizer at the
limiting correlation, hence ordinary. Its summed squared norm is one,
as for the sequence, so convergence is strong in \(L^2\) and \(L^1\).

The planar maximizing family in this compact correlation interval is
therefore \(L^1\)-compact. Impose DG11's closed condition
\(E_J(B)\ge\delta\), continuous in \(L^1\) by conditional
expectation and bounded circular masses. This leaves a compact family
containing every strict-failure representative and no centered simplex.
On it
\(\max_R\min_{i<k}\gamma_2(B_i\cap RB_k)\) is strictly positive
by the rigidity above and continuous: overlaps obey a uniform
symmetric-difference bound in \(B\), and rotations are strongly
continuous on Gaussian \(L^1\). A positive minimum and a finite open
subcover give fixed witness reflections and (DG59). If this compact
family is empty the assertion is vacuous.

Each overlap in (DG59) equals the mass of its unordered type in either
half-space. Removing a thin strip and a distant exterior loses at most
\(\eta/2\) of each type, uniformly in \(B\) and the reflection.
Then (DG56), uniformly over the compact correlation interval, gives a
fixed positive scalar polarization defect at the selected reflection.
The constants are qualitative. No estimate excluding this defect or
the entire remaining maximizing class is proved; the original global
inequality remains open.

The reflection checker exhausts the eight triangle signs and checks
the absolute-incidence identity in all 81 label cases. These support
(DG55)--(DG56). The measurable order, analytic scores and compactness
arguments are written proofs, not conclusions from finite sampling.

## A global obstruction to reducing every circular slice to consecutive arcs

### Theorem DG14 (balanced global consecutive-arc reduction is false)

For every fixed \(0<\rho<1\), there is an ordinary balanced planar
partition \(A\) such that
\[
 S_\rho(A)>\sup\{S_\rho(C):p_C(r)=p_A(r)\text{ a.e., and each slice
 of }C\text{ consists of consecutive arcs}\}.
\]
The competitors may choose arbitrary measurable rotations and either
cyclic label order at each radius. Empty arcs and a full-circle arc are
allowed. There is also such an example with every circular mass bounded
below by one fixed positive constant. Thus reducing a single global
partition to consecutive arcs while preserving its circular masses is
not a universally increasing operation. No comparison of this example
with the simplex is asserted.

The [joint angular note](THREE_CELL_JOINT_ANGULAR_ENVELOPE.md) already
excludes the analogous reduction for a bilinear pair of slices away from
balanced masses. The issue here is different: radial self-interactions
could in principle repair a failed pair reduction. The proof controls
all those self-interactions and arbitrary outer rotations explicitly.
The inspected primary Corollary A.3 in
[Heilman, arXiv:2608.00333v1](https://arxiv.org/html/2608.00333v1)
bounds the sum of individually aligned same-label arc pairings; it does
not assert their joint attainment by a three-cell partition on every
circle. Its upper bound and (DG16) remain valid.

**An angular gap.** Write \(k_\kappa=e^{\kappa\cos\theta}/I_0(\kappa)\).
Let \(h_\kappa\) be its convolution with the half-circle centered at
zero, and let \(g_\kappa(t)\) be the pairing of that half-circle
with the arc of length \(1/3\) centered at angle \(t\), using
probability measure on both circles. Direct endpoint differentiation gives
\[
 h_\kappa'(t)=-\frac{\sinh(\kappa\sin t)}{\pi I_0(\kappa)},\qquad
 g_\kappa'(t)=\frac{h_\kappa(t+\pi/3)-h_\kappa(t-\pi/3)}{2\pi}.
\]
Both functions are even and strictly decreasing in circular distance
from zero on \((0,\pi)\). For \(0\le t\le\pi/6\), the interval
\([\pi/3-t,\pi/3+t]\) lies in \([\pi/6,\pi/2]\), where
\(\sin v\ge1/2\). Integrating the displayed derivatives therefore gives
\[
 g_\kappa(0)-g_\kappa(\pi/6)
                  \ge\frac{\sinh(\kappa/2)}{72 I_0(\kappa)}>0.
 \tag{DG60}
\]
For the monotonicity of \(g\), the circular distance of
\(t+\pi/3\) from zero is larger than that of \(t-\pi/3\) for
\(0<t<\pi\); use the monotonicity of \(h\).

Give the outer slice masses \((1/2,1/2,0)\), with its first two
half-circles centered at zero and \(\pi\). The inner slice of masses
\((1/3,1/3,1/3)\) can put label 1's arc at zero and label 2's arc at
\(\pi\), leaving label 3 the two complementary arcs. Its cross
pairing is \(2g_\kappa(0)\). In any consecutive inner slice the
first two centers differ by \(\pm2\pi/3\), so their offsets from the
two outer centers have circular separation \(\pi/3\). At least one
offset is at distance at least \(\pi/6\). Its total pairing is at
most \(g_\kappa(0)+g_\kappa(\pi/6)\). Formula (DG60) is thus a
uniform strict gap over both cyclic orders and all inner rotations.

**The radial construction and constants.** Put \(d=1-\rho^2\),
\(d\mu(r)=r e^{-r^2/2}dr\), and
\[
 Q(r,s)=d^{-1}e^{-\rho^2(r^2+s^2)/(2d)},\qquad
 \kappa(r,s)=\rho rs/d.
\]
The polar Mehler kernel is \(Qe^{\kappa\cos\theta}\).
Take an outer annulus \(V=[2,3]\) of radial mass
\(\beta=e^{-2}-e^{-9/2}\), and an inner annulus
\(U=[1,r_\alpha]\subset[1,5/4]\) of mass \(\alpha\).
On \(U\), use the nonconsecutive inner slice just described; on
\(V\), use the two fixed half-circles. At each remaining radius use
a single label for the full circle. Partition those remaining radii
into three sets of respective \(\mu\)-masses
\[
 1/3-\alpha/3-\beta/2,\quad
 1/3-\alpha/3-\beta/2,\quad 1/3-\alpha/3.
\]
They are positive if \(\alpha\le1/12\), since \(\beta<1/3\), and
sum to \(1-\alpha-\beta\). Successive radial quantiles on the
complement of \(U\cup V\) specify these sets. The full partition
is ordinary and exactly balanced. All interactions with those pure
radii depend only on the circular masses and hence cancel in every
comparison under consideration.

For arbitrary rotating outer half-circles with center \(t(s)\), put
\(z=\beta^{-1}\int_V e^{it(s)}d\mu(s)\) and \(v=1-|z|^2\).
Their same-annulus energy is at most the aligned value minus
\[
                        c\beta^2v,\qquad
 c=\frac{8\rho}{\pi^2d^2}e^{-9\rho^2/d}>0.
 \tag{DG61}
\]
Indeed the two-label angular pairing at relative rotation \(t\) is
\(1/2+(4/\pi^2)\sum_{m\ge1\text{ odd}}I_m(\kappa)
\cos(mt)/(m^2I_0(\kappa))\). All Fourier coefficients are positive.
Keep its first mode, use \(I_1(\kappa)\ge\kappa/2\), and use
\(r,s\in[2,3]\) to obtain (DG61).

The following constants are positive and finite for each fixed \(\rho\):
\[
 \begin{split}
 \delta&=\frac{e^{-169\rho^2/(32d)}}{72d}\sinh(\rho/d),\\
 L&=\frac{2}{3\pi d}\sinh(15\rho/(4d)),\qquad
 C=\sqrt2\pi L,\\
 M&=d^{-1}e^{25\rho/(16(1+\rho))},\qquad
 K=M+C^2/(4c).
 \end{split}
 \tag{DG62}
\]
For \(r\in U,s\in V\), the factor \(Q\) is at least
\(d^{-1}e^{-169\rho^2/(32d)}\), and \(\kappa\ge2\rho/d\).
Thus (DG60), after multiplication by \(QI_0\), gives cross gap at
least \(\delta\) against any consecutive inner slice and a fixed
outer phase. The cross pairing changes by at most \(L|t-t_0|\)
under rotation of the outer half-circles: its derivative is bounded by
\((2/3)Q\sinh\kappa/\pi\), since the inner first two masses sum
to \(2/3\), and \(\kappa\le15\rho/(4d)\).
The inner Mehler kernel is at most \(M\), using
\(|x|,|y|\le5/4\) and
\(M_\rho(x,y)\le d^{-1}
e^{\rho(|x|^2+|y|^2)/(2(1+\rho))}\).

Choose \(t_0=\arg z\), arbitrarily if \(z=0\). For circular
distance \(|t-t_0|\le\pi\),
\[
 \beta^{-1}\int_V|t(s)-t_0|d\mu(s)
 \le\frac\pi2\sqrt{2(1-|z|)}\le\frac\pi{\sqrt2}\sqrt v.
\]
This follows from \(|t-t_0|\le(\pi/2)|e^{it}-e^{it_0}|\) and
Cauchy--Schwarz. Rotate the constructed partition to phase \(t_0\);
its stability is unchanged. Combining the uniform cross gap, the
cross rotation error, (DG61), and the worst possible inner self cost
gives, for every consecutive-arc competitor \(C_*\) with this profile,
\[
 \begin{split}
 S_\rho(A)-S_\rho(C_*)&\ge
 2\delta\alpha\beta+c\beta^2v-C\alpha\beta\sqrt v-M\alpha^2\\
 &\ge2\delta\alpha\beta-K\alpha^2.
 \end{split}
 \tag{DG63}
\]
The factor two in the cross terms counts both ordered radial pairs.
The final step completes the square in \(\beta\sqrt v\).
This estimate controls all measurable outer rotations and all inner
rotations and cyclic orders; it is not a comparison only with a fixed
pair of reference angular shapes.

For an explicit admissible choice set
\[
 \alpha=\tfrac12\min\{1/12,e^{-1/2}-e^{-25/32},\delta\beta/K\},
 \qquad r_\alpha=\sqrt{-2\log(e^{-1/2}-\alpha)}.
 \tag{DG64}
\]
Then \(1<r_\alpha<5/4\), its annulus has exactly the stated mass,
and (DG63) is at least \(\delta\alpha\beta>0\).
This proves the all-parameter assertion with the order of quantifiers
stated: the construction and its positive gap may depend on \(\rho\).

**Strictly positive profiles.** Let \(G>0\) be the lower gap in
(DG63). Replace the profile by
\(p^\tau=(1-\tau)p+\tau(1/3,1/3,1/3)\). It is balanced and every
coordinate is at least \(\tau/3\). Reallocate angular mass within
each slice to realize it by an ordinary \(A^\tau\) satisfying
\(\sum_i\gamma_2(A_i^\tau\triangle A_i)\le2\tau\).
For example remove the required mass from surplus labels in label
order, starting clockwise at angle zero, and assign the pieces to
deficit labels in label order. The profile and angular shapes here
have finitely many radial classes, so these choices are measurable.

Conversely every consecutive-arc partition of profile \(p^\tau\)
can be changed to one of profile \(p\) while keeping its starting
endpoint and cyclic order, with summed symmetric difference at most
\(8\tau\). In fact shifting the two cumulative endpoints bounds
the sum by \(4\int\|p^\tau-p\|_1d\mu\le8\tau\).
For any two partitions \(B,D\), Markov contraction gives
\(|S_\rho(B)-S_\rho(D)|\le2\sum_i\gamma_2(B_i\triangle D_i)\).
Thus the new strict gap is at least \(G-20\tau>0\) when
\(0<\tau<G/40\). The obstruction persists with uniformly interior
circular masses.

The [directed checker](check_three_cell_global_arc_reduction.py) certifies
one explicit instance at \(\rho=1/2\), \(\alpha=10^{-8}\), with
\[
 5.633936166484404\,10^{-12}
 <2\delta\alpha\beta-K\alpha^2
 <5.633936166484406\,10^{-12}.
 \tag{DG65}
\]
It also checks the exact six labeled consecutive orders, radial
admissibility, the remaining positive masses, and a rejected oversized
annulus control. The all-parameter and all-competitor proof is
(DG60)--(DG64); the selected directed run supplies a concrete numerical
instance only. The [verification record](three_cell_global_arc_reduction_verification.json)
pins its environment, inputs and internal scope. Replay from the root:

```sh
.venv/bin/python research/plurality-stability/check_three_cell_global_arc_reduction.py
```

The example obstructs an operation on arbitrary balanced partitions.
It does not disprove a consecutive-arc classification restricted to
actual global maximizers. In particular its original unregularized
form has first moments on one line and is excluded as a global
maximizer by the full active-rank theorem. The obstruction also exists
at correlations already covered by DG10, where the simplex comparison
is known to hold.

## Adopting the direct comparison across the scalar obstruction

### Corollary DG15 (a direct interval beyond the first scalar threshold)

Let \(\rho_{\mathrm s}\) be the first separated scalar obstruction
defined by (ST1) in the
[threshold proof](THREE_CELL_BEYOND_SCALAR_THRESHOLD.md). Then
\(11/25<\rho_{\mathrm s}<1/2\), and there is \(\varepsilon>0\)
such that every ordinary balanced Gaussian three-cell partition in
every \(d\ge2\) satisfies
\[
 \widetilde Q_\rho(1_A-1_T)\le2L_\rho(A;T)
                  \qquad(0<\rho<\rho_{\mathrm s}+\varepsilon).
 \tag{DG66}
\]
The threshold and extension are implicit. On a neighborhood of
\(\rho_{\mathrm s}\), the planar comparison has the stronger weighted
radial-mass deficit (ST2), and equality forces the simplex. This does
not assert that the separated scalar criterion itself remains valid
above its first obstruction.

The proof source is (ST1)--(ST6) and their subsequent compactness
argument; its complete proof remains there. Its essential new point
is strict angular loss on every nonzero scalar null profile. At the
threshold the scalar form is nonpositive. A bounded null profile
satisfies the Euler equation (ST4), hence has analytic representatives
on the positive radial axis. Its mass coordinates are interior except
at isolated zeros and are not all constant. The penalty \(\Gamma\)
of DG6 is strictly positive on a positive-measure set of unequal
interior profile pairs, so (DG29) is strictly negative on that null
profile.

Uniformity across the threshold uses the compact normalized radial
operator for correlations strictly below \(1/2\). A sequence nearly
violating (ST2) either has weighted mass norm tending to zero, in which
case the adopted annular comparison (DG48) gives a uniform contradiction,
or has a nonzero strongly convergent scalar null limit. In the latter
case positivity and continuity of \(\Gamma\) at almost every limiting
pair, followed by Fatou, give a strictly positive angular loss, again a
contradiction. The source also proves density of bounded mean-zero
tests in the weighted hyperplane, compactness and continuity of its
uncentered operator, and the relaxation-to-null step. These hypotheses
are needed; pointwise strictness alone would not supply a neighborhood.

This adopts the direct threshold comparison, rather than the
fixed-correlation interface isolation theorem. Maximizer existence and
dimension reduction extend its planar bound to all competitors in
\(d\ge2\), and (B1) gives (DG66). The internal audit checked the
weighted origin and Gaussian tail bounds, feasible bounded limiting
profiles, analytic interiority, and the two norm cases in the uniformity
argument. No new numerical spectral cutoff, independent replay, formal
verification or external review is claimed.

In particular the lower bound \(11/25\) in the necessary failure
localizations (DG49) and (DG59) may be strengthened to one fixed
\(v>\rho_{\mathrm s}\), by taking
\(v=\rho_{\mathrm s}+\varepsilon/2\) after shrinking
\(\varepsilon\) so that \(v<1/2\). Their qualitative upper bound
and their annular and reflection witnesses remain valid. The global
comparison beyond this direct interval and below the near-one interval
remains open.

## Isolating radial directions and proving their angular cost

### Theorem DG16 (finite adverse modes and directional angular domination)

The [radial and angular mode proof](THREE_CELL_RADIAL_ANGULAR_MODE_ISOLATION.md)
establishes the following three statements, with their proofs and
normalizations in RMA1--RMA14.

First, the compressed weighted scalar operator is compact at every
\(0<\rho<1\), not only below \(1/2\). The new lower bound
\(J_\rho(t)\ge j_I\min(\sqrt t,1)\) on each compact correlation
interval restores an integrable Hilbert--Schmidt tail at all such
correlations. All positive or nearly positive scalar modes therefore
lie in a uniformly finite-dimensional space. On its complement the
scalar form has a fixed negative gap.

Second, Hermite interpolation of the exponential angular kernel gives
global moment-loss bounds (RMA5)--(RMA7), with exact angular optimizations
over finitely many endpoints and allowing disconnected cells. At
quadratic order about balanced circular masses, some finite angular
order outweighs every positive scalar radial direction uniformly on
each compact correlation interval, by (RMA12). The order is qualitative;
it is not uniform as \(\rho\uparrow1\). This directional assertion
does not bound every finite-amplitude mass profile.

Third, at \(\rho=99/100\) the two-dimensional radial space spanned by
the explicit \(f\) in (DG21) and \(g\) following (DG40) has exactly
one positive direction for the older exact-arc/secant quadratic bound.
The second angular mode makes the entire block negative:
\[
 q_{99/100}^{[2]}(c_f f+c_g g)\le-(c_f^2+c_g^2)/30000.
 \tag{DG67}
\]
The directed closed-form checker proves the matrix inequality, and the
uniform Taylor argument gives a common small-profile comparison over
all unit directions in that block, with arbitrary angular shapes.
This repairs an isolated adverse block rather than only one test
function. It does not claim that this block contains every adverse
mode of the unrestricted scalar operator.

The remaining direct step is quantitative domination at every feasible
amplitude: the unavoidable angular loss plus the negative radial
remainder must exceed the positive finite-mode contribution in (RMA4).
The source proves the complete quadratic directional statement and
one explicit block certificate; this global nonlinear amplitude
inequality remains open.

## Remaining obligation and reproducible controls

Completing the all-positive-correlation statement requires a comparison
across the remaining intermediate range above the qualitative direct
interval in DG15 and below the separate near-one interval. DG10 still
supplies the explicit numerical cutoff \(11/25\). In the polar
route, the unresolved all-partition comparison is (DG4), retaining enough of
\(\mathcal C_\rho\) or using a different mechanism. Global balance
alone does not remove \(\mathcal R_\rho\); equal circular masses do.
Proposition DG5 rules out extending the particular scalar operator
estimate in (DG17) through \(\rho=1/2\), even with optimal additive
centering. Proposition DG6 supplies an explicit nonlinear correction
retaining moment closure, and Corollary DG7 shows it handles that
obstructing profile. Its elementary version is itself insufficient at
\(\rho=3/5\). Proposition DG8 obtains the optimal first-moment
coefficient near balanced profiles, repairs that particular profile at
\(\rho=3/5\), and proves that even this sharpened separated criterion
fails at \(\rho=7/10\). Proposition DG9 retains the exact arc
interaction and uses a stronger measure-dependent secant penalty. It
handles the explicit profile at \(\rho=9/10\), but its own sufficient
all-profile criterion fails at \(\rho=99/100\) for a different bounded
mean-zero profile. A uniform proof therefore still needs a sharper
joint angular estimate, additional modes or a different mechanism.
Proposition DG3 rules out proving the stronger mass-only inequality after
the angular Cauchy step throughout \((0,1)\). It does not rule out
stronger coupled angular estimates, a different direct proof, or the
classification of actual maximizers.
Theorem DG14 rules out a universal global reduction to consecutive arcs
at fixed circular masses, even for uniformly positive profiles. A proof
using such a reduction specifically for actual global maximizers would
still need their maximizing-score information; it cannot follow from a
monotonic rearrangement available to every balanced partition.

The supplementary [checker](check_three_cell_direct_global.py) checks the
quadrant formula against independent numerical integration, exact radial
moments, and finite-correlation values of \(E_\rho\). It also checks the
near-one limit with the factor \(s_\rho^2/(1-\rho)\to2\).
Run from the repository root:

```sh
python3 research/plurality-stability/check_three_cell_direct_global.py
```

The retained [verification record](three_cell_direct_global_verification.json)
states the environment, inputs and actual internal outcome. Numerical
quadrature is a normalization diagnostic, not an interval certificate or
a proof of (DG12). The unrestricted identities and obstruction use the
written arguments above. The circular comparison imports (DG0) with the
source-inspection limits stated above. No independent replay, formal
verification or external mathematical review is asserted.


## Finite-amplitude radial and angular comparison

### Proposition DG17 (global comparison through rho equal to 221/500)

For every ordinary balanced Gaussian three-cell partition in every
\(d\ge2\), the direct inequalities (DG14) hold for
\(0<\rho\le221/500=0.442\). This improves the explicit cutoff in DG10.
It uses no continuation or closeness hypothesis.

The proof is [FA1--FA7](THREE_CELL_FINITE_AMPLITUDE_COMPARISON.md).
On the new interval \([11/25,221/500]\), decompose the weighted circular
mass profile into one explicit radial trial and its orthogonal remainder,
with coefficient norm \(A\) and remainder norm \(G\). The directed
certificate and the written Hilbert--Schmidt block argument give

\[
 \sum_i q(z_i)\le\frac{31}{10000}A^2
          +\frac{2}{\sqrt{50000}}AG-\frac75G^2.
\]

Any profile with nonnegative scalar contribution must have
\(A<27/1000\), since larger amplitude of the pole at the origin would
force a correction violating this quadratic bound. Uniform angular
coercivity on a mass cube extends to arbitrary actual profiles by charging
the error outside the cube to the negative remainder. On an explicit
annulus it yields

\[
 \mathcal L(p)\ge\frac3{800}A^2-\frac{20}{27}G^2.
\]

Consequently the full stability excess is bounded by

\[
 -\frac{13}{20000}A^2+\frac{2}{\sqrt{50000}}AG
                         -\frac{89}{135}G^2\le0.
\]

The exact determinant is positive. FA7 and the
[verification record](three_cell_finite_amplitude_verification.json)
cover the whole new correlation interval and all integration domains,
including analytic infinite tails. The internal full run passed in
243.72 seconds. Together with DG10 and the same maximizer dimension
reduction used there, this proves the stated all-partition cutoff.
This is a finite-amplitude result on an explicit interval, not the
all-correlation endpoint of the mode-isolation program in DG16.
The near-one interval and the qualitative scalar-threshold extension
remain separately scoped; no overlap closing the remaining range is
established.


## Angular compatibility shared across radii

### Proposition DG18 (global moment loss and an all-amplitude angular hierarchy)

For every fixed positive correlation and every ordinary balanced planar
partition, (GAC7) of the
[global angular compatibility note](THREE_CELL_GLOBAL_ANGULAR_COMPATIBILITY.md)
is a valid upper bound for the stability excess, at least as strong as
(DG29). It retains the actual common first-moment field through a global
rank-one minorant of the positive radial loss kernel. Its exact optimization
reduces to one compact two-dimensional matrix disk and pointwise convex
three-point geometric-median supports. No small-amplitude hypothesis is used.

The shared-radius constraint can matter even when every pair envelope is
exact. For the three half-mass types, the sharp common-slice deficit is
(GAC8), strictly positive whenever all three weighted types occur. Every
separately optimized pair penalty, including the polynomial ones, vanishes
on this family. The loss is stable under integrated mass errors by (GAC9).

There is also a global finite-angular-order envelope (GAC12), using the same
two radial fields in every pairing. For all mass profiles and all
\(0<\rho\le b<1\), its excess over the true prescribed-profile optimum is
at most
\(4(1-b)^{-2}[(4b)/(1+3b)]^{n+1}\).
Its supremum can be approached with at most\(6(n+1)\) switches per slice.
This uniform error is global in amplitude; it does not prove the envelope
is nonpositive. The coupled radial optimization and the all-correlation
comparison remain open. DG18 supplies no new correlation cutoff beyond DG17.


## Sharp retention of the second angular mode

### Proposition DG19 (global two-mode support and an explicit adverse-profile comparison)

For a circular set of mass \(p\), and nonnegative weights \(r,s\), the
[sharp two-mode theorem](THREE_CELL_SHARP_TWO_MODE_ENERGY.md) proves

\[
 \max_E\{r|\widehat{1_E}(1)|^2+s|\widehat{1_E}(2)|^2\}
 =\frac{\sin^2(\pi p)}{\pi^2}
                 \max\{r+s\cos^2(\pi p),s\}.
\]

The maximizer is an arc or two equal antipodal arcs. The optimal
second-mode excess coefficient relative to the first-mode squared cap
deficit is \(\sin^2(\pi p)\), strengthening AM7's sufficient coefficient.
TME3--TME7 give global pair and shared-radius matrix bounds and a further
nonnegative angular loss in the polar comparison at every amplitude.
They retain the two harmonics together; treating the second arc moment
as an independent norm cap would be false.

A full directed four-dimensional mass cover certifies both penalties on
\(|p_i-1/3|\le1/20\). In particular, at \(\rho=3/5\), for the radial
function \(f(t)=e^{-t}+2e^{-16t}-21/34\), every centered unit label vector
\(v\) and every \(|\epsilon|\le1/40\), arbitrary angular shapes with
masses \(1/3+\epsilon v_i f(t)\) have deficit at least
\(\epsilon^2/500\). The scalar contribution is positive; the two
certified angular losses outweigh it. This is an explicit restricted-profile
comparison beyond DG7's qualitative amplitude statement, not a global
optimality theorem at \(3/5\).

TME6 also proves a limit of the first-mode framework: for the bounded
six-exponential witness at \(99/100\) in a two-equal-mass label direction,
one common arc field attains every exact first-mode pair optimum.
Replacing that penalty by its fully shared version therefore leaves
DG39's positive quadratic coefficient unchanged. The assertion concerns
that remainder envelope and does not falsify actual simplex optimality.
Higher angular information or a stronger remainder comparison is needed.
The all-positive-correlation inequality remains open.


## Adopting the all-mass scalar comparison through 0.53

### Corollary DG20 (unrestricted scalar comparison through 53/100)

For every \(d\ge2\), every ordinary measurable balanced Gaussian
three-cell partition \(A\), and every \(0<\rho\le53/100\),

\[
 S_\rho(A)\le S_\rho(T),\qquad
 \widetilde Q_\rho(1_A-1_T)\le2L_\rho(A;T).
\]

The first assertion is the current theorem of the
[all-mass scalar note](THREE_CELL_GLOBAL_SCALAR_THROUGH_HALF.md),
GH1--GH10, adopting its written proof and internal directed certificates.
It applies to all circular mass profiles and arbitrary angular shapes.
The second assertion follows by the exact expansion (B1), as in DG4.
No continuation of maximizing partitions or closeness hypothesis is used.

The new angular envelope (GH2) is global in both probability vectors.
Its positive radial feature expansion and mean-zero centering have a
certified norm below two on the full interval. The high-mode tail above
one half uses the exponential lower bound for the full row weight in
(GH9), rather than extending the polynomial trace estimate beyond its
range. The retained through-half cover reports an upper norm below
1.92573 on 240 covering intervals; the adjoining cover through53/100
reports an upper norm below1.98823 on160 covering intervals, including
the analytic tail. Both saved runs pass and the registered checker hash
matches the current source. This adoption is not an independent replay,
formal verification or external review.

The explicit lower interval is therefore stronger than DG17. The
shared-radius and two-mode bounds DG18--DG19 remain useful potential
corrections on the remaining range; they are not needed for (GH2)'s
current certified interval. No overlap with the qualitative near-one
interval closes all positive correlations here.


## Paying for the radial obstruction with a stronger global variance loss

### Corollary DG21 (unrestricted direct comparison through 61/100)

For every \(d\ge2\), every ordinary measurable balanced Gaussian
three-cell partition \(A\), and every \(0<\rho\le61/100\),

\[
 S_\rho(A)\le S_\rho(T),\qquad
 \widetilde Q_\rho(1_A-1_T)\le2L_\rho(A;T).
\]

The [same-label comparison](THREE_CELL_SAME_LABEL_VARIANCE_PENALTY.md),
SLV1--SLV7, strengthens the sharp scalar bound by
\(17\rho\pi(1-\rho^2)\operatorname{Var}_{\mu_\rho}(p)/1000\),
where \(d\mu_\rho\) is proportional to
\(\sqrt t e^{-\rho^2t/(1-\rho^2)}d\nu(t)\).
This loss holds for every feasible amplitude and every \(0<\rho<1\).
It combines the positive slack when two slices favor the same label
with the separate quartic difference slack. The key cumulative-density
comparison SLV2 prevents the same-label contribution from being
weakened by an exponentially small factor.
The convex Laplace-transform calculation SLV3 also brackets the optimal
pure same-label linear constant between 0.1952578742 and 0.1952578830.
Combined with the dominant-label bound SV3, it gives the larger
all-correlation reserve \(781\rho\pi(1-\rho^2)/32000\) in SLV8.
The independent radial cover below uses only the smaller reserve.

The resulting corrected radial operator has a directed norm upper
bound below \(1999/1000\) throughout \([29/50,61/100]\).
The full 40-slab cover, including the mean-zero transport error and
infinite feature and radial tails, reports an upper bound below
1.998842. The preceding sharp scalar theorem supplies correlations
through \(29/50\). Existence and planar reduction of maximizers,
equality for uniform circular masses, and the exact expansion B1
then give the displayed unrestricted comparison. Equality on the
new interval determines the simplex, up to the usual symmetries.

This is a global repair of the scalar criterion by angular loss;
it uses neither continuation nor a small-distance hypothesis.
The proof and certificate have internal verification only. The
all-positive-correlation inequality remains open, since no sign proof
for the corrected operator on the whole remaining interval is supplied.
