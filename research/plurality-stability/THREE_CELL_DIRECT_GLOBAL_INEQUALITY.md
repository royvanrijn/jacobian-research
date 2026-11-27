# A direct three-cell comparison and the angular information it requires

The unrestricted inequality
\(\widetilde Q_\rho(e)\leq2L_\rho(A;T)\) remains open here.
The exact expansion (B1) in the [continuation note](POSITIVE_CORRELATION_CONTINUATION.md)
makes it equivalent to \(S_\rho(A)\leq S_\rho(T)\), including for
partitions far from the simplex. Proposition DG4 below now gives a direct
computer-assisted proof for every balanced three-cell partition in every
dimension \(d\geq2\) throughout \(0<\rho\leq21/50\). The full range
\(21/50<\rho<1\) remains open here. The argument also proves the comparison
at every positive correlation when the three labels occupy equal portions
of almost every circle. For general balanced partitions it isolates an
exact missing inequality and proves that discarding the angular
Cauchy--Schwarz remainder cannot establish it near correlation one.
The counterexample is to an upper-bound strategy, not to simplex optimality.

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
dominance of the right side of (DG29) remains open; Corollary DG7 does
not extend the unrestricted cutoff in Proposition DG4.

Run both the directed witness and the exact rational controls from the
repository root with the existing environment:

```sh
.venv/bin/python research/plurality-stability/check_three_cell_coupled_radial.py
```

The verification record pins this checker, its exact inputs and the
current proof source. These new arguments have internal checking only;
independent replay, formal verification and external review remain
unperformed.

## Remaining obligation and reproducible controls

Above the newly certified cutoff \(21/50\), an unrestricted direct proof
must establish (DG4), retaining enough of
\(\mathcal C_\rho\) or using a different mechanism. Global balance
alone does not remove \(\mathcal R_\rho\); equal circular masses do.
Proposition DG5 rules out extending the particular scalar operator
estimate in (DG17) through \(\rho=1/2\), even with optimal additive
centering. Proposition DG6 supplies an explicit nonlinear correction
retaining moment closure, and Corollary DG7 shows it handles that
obstructing profile. A bound uniform over all profiles is still needed.
Proposition DG3 rules out proving the stronger mass-only inequality after
the angular Cauchy step throughout \((0,1)\). It does not rule out
stronger coupled angular estimates, a different direct proof, or the
classification of actual maximizers.

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
