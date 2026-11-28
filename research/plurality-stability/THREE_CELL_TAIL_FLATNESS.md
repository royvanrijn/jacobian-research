# Straightness of bounded isolated ends

Fix a balanced planar three-cell Gaussian maximizer at \(0<r<1\).
Suppose that in an orthonormal coordinate system about the Gaussian
origin almost every vertical section of a half-plane \(x>S\)
has one switch between the same two labels at \(y=h(x)\), with
\(\operatorname*{ess\,sup}_{x>R}|h(x)|/x\to0\). Then the interface
is exactly the straight ray \(h=b\) on a tail, where \(b\)
is the offset determined by the pair pressure. The sublinear
graph theorem first gives an analytic representative with
bounded height. The straightening proof uses the constant
pair-score equation on this graph.
[Boundedness of sublinear isolated graphs](THREE_CELL_ISOLATED_END_RIGIDITY.md#every-sublinear-isolated-graph-is-bounded)

The pressure equation yields a quadratic residual in the transformed
height at the observation point, with a further term that is linear
there and quadratic in the source height. The Mehler telescoping
identity improves the decay of a Gaussian convolution of the
transformed height. Strip interpolation gives an initial Gaussian
decay exponent. Sector interpolation and a real shift of the
Gaussian inversion contour then increase that exponent to the
threshold of the one-sided Hardy theorem. The shift recovers the
Gaussian decay lost in passing from the wider convolution back to
the narrower one in the pressure equation.

## The stationary residual

Let \(\sigma=\sqrt{1-r^2}\), let \(p_\sigma\) be the centered
normal density of standard deviation \(\sigma\), and write
\(\psi(u)=2\Phi(u/\sigma)-1\). The pair contrast equals its
pressure \(\beta\) on \(h\). The bounded-end theorem gives
\(h(x)\to b\) and \(\psi((1-r)b)=\beta\). Put
\(q=h-b\), \(d=(1-r)b\), and
\[
   v(x)=\psi(d+q(x))-\psi(d),\qquad
   M(X)=\sup_{x\ge X}|v(x)|.
\]
The functions \(|v|\) and \(|q|\) are comparable on the tail because
\(\psi'\) is positive on every compact interval. In particular
\(M(X)\to0\).

Choose \(0<a<b_0<c<r\) and take \(X\) large enough that \(aX>S\).
Define \(q_X=q\,1_{\{x\ge aX\}}\) and
\(Q=v\,1_{\{x\ge aX\}}=\psi(d+q_X)-\psi(d)\). Write
\(T_r w(x)=\int w(s)p_\sigma(s-rx)\,ds\). For \(x\ge X b_0/r\),
the conditional pair equation gives
\[
 \left|\int [\psi(d+q_X(s)-r q(x))-\psi(d)]
                       p_\sigma(s-rx)\,ds\right|
      \le 2\Pr[N(rx,\sigma^2)<aX].
\]
Indeed on \(s\ge aX\) the two-label formula is exact, and replacing
the conditional pair contrast elsewhere changes a bounded integrand
only on \(s<aX\). Also \(rx-aX\ge (1-a/b_0)rx\). Expanding the
observation-side displacement and using the pressure equation to
eliminate its mixed quadratic term gives, for \(f_X=T_rQ-rQ\),
\[
 |f_X(x)|\le
    C\bigl(|Q(x)|^2+|Q(x)| T_r(Q^2)(x)\bigr)+C e^{-\kappa x^2}
       \qquad (x\ge X b_0/r).                              \tag{TF1}
\]
Here and below constants may depend on the fixed parameters and the
end. In fact,
\(\psi'(d+q_X)=\psi'(d)+kQ+O(Q^2)\) with
\(k=\psi''(d)/\psi'(d)\), by Taylor expansion and comparability of
\(|q_X|\) with \(|Q|\). Expanding the pressure equation to first
order in \(r q(x)\) gives
\[
  (1-krq(x))\,T_rQ(x)=r q(x)\psi'(d)
       +O\!\left(q(x)^2+|q(x)|T_r(Q^2)(x)+e^{-\kappa x^2}\right).
\]
The prefactor is bounded away from zero for large \(X\).
Dividing and using \(rQ(x)=r q(x)\psi'(d)+O(q(x)^2)\) gives
(TF1). On the entire line \(\|f_X\|_\infty\le (1+r)M(aX)\).

## The improved smoothed bound

Write \(P_u\) for convolution with the centered normal density of
standard deviation \(u\). The Mehler telescoping identity for
bounded \(Q\) and \(f_X=T_rQ-rQ\) is
\[
   G(t):=P_1Q(t)=\sum_{k\ge1} r^{k-1}P_1 f_X(t/r^k).       \tag{TF2}
\]
[Derivation of the identity](GAUSSIAN_EXTERIOR_LINEAR_RIGIDITY.md#iteration-at-expanding-observation-points)

Set \(m=M(aX)\), \(m_1=M((b_0/r)X)\), and \(D_0=m_1^2+m_1 m^2\).
On \(x\ge (b_0/r)X\), (TF1) gives \(C D_0+C e^{-\kappa_0 X^2}\)
as a uniform bound on \(f_X\), because \(T_r(Q^2)\le m^2\).
For \(t\ge cX\) and \(k\ge1\), the conditional center \(t/r^k\)
in (TF2) is at least \(cX/r\). Its variance-one Gaussian
therefore sees the complementary half-line with probability at
most \(C e^{-\kappa_1 X^2}\). Summing the geometric coefficients
gives, with some \(\kappa_2>0\),
\[
       |G(t)|\le C(D_0+ e^{-\kappa_2 X^2})\qquad(t\ge cX).
                                                               \tag{TF3}
\]

## Interpolation on an imaginary strip

Fix \(0<\nu<1-r^2\). Choose a small constant \(\eta>0\), with
\(2\eta\le\kappa_2\), and put
\[
   A=\max(m,e^{-\eta X^2}),\qquad
   D=\max(D_0,e^{-2\eta X^2}),\qquad B=D/A .
\]
For large \(X\), \(0<A,B<1\) since \(D_0\le A(m_1+m_1m)\).
The entire Gaussian convolution
\(G\) satisfies \(|G(u+iy)|\le A e^{y^2/2}\) for every real
\(u,y\), because \(\|Q\|_\infty\le A\). By (TF3),
\(|G(u)|\le C A B\) when \(u\ge cX\). At any \(t\ge rX\)
let
\[
   R=(r-c)X,\qquad L=\log(1/B),\qquad
   \tau=L/R^2,\qquad H=\sqrt{L/(\tau+1/2)} .
\]
The entire function \(G(z)e^{-\tau(z-t)^2}/A\) has modulus
at most \(C' B\) on the real line. On \(|u-t|\le R\) this
comes from (TF3), and outside that interval it follows from
the real bound \(|G(u)|\le A\) and the multiplier.
On the line of imaginary part \(H\) it has modulus at most
\(e^{(\tau+1/2)H^2}=e^L\).
Interpolating the two logarithmic bounds gives
\[
 |G(t+iy)|\le C' A
             \exp\!\big(-L+2\sqrt{L(\tau+1/2)}\,|y|\big)
                \qquad (|y|\le H).                         \tag{TF4}
\]
For the interpolation, multiply the function by an exponential
linear in \(z\) so the two boundary bounds become one. Apply
the maximum modulus principle on finite rectangles in the strip
and let their vertical sides tend to infinity. The factor
\(e^{-\tau(z-t)^2}\) makes the vertical bounds tend to zero
uniformly. At \(t+iy\) its modulus is at least one, so removing it
preserves (TF4). The lower strip follows by conjugation.

For standard real normal \(Z\), direct Gaussian integration gives
\[
                  P_\sigma Q(t)=\mathbb E G(t+i rZ).
                                                               \tag{TF5}
\]
Indeed insert the Gaussian integral for \(G\). Fubini is justified
by \(\mathbb E e^{r^2 Z^2/2}<\infty\). Averaging the resulting
complex exponential multiplies \(p_1(t-s)\) by
\(\sigma^{-1}\exp[-r^2(t-s)^2/(2\sigma^2)]\), giving \(p_\sigma\).
On \(|rZ|\le H\), (TF4) and the normal exponential moment bound
the integral by \(C A B^{1-r^2-2r^2\tau}\). On the complement,
the bound \(|G(t+iy)|\le A e^{y^2/2}\) gives
\[
  C A \exp[-(1-r^2)H^2/(2r^2)]
     = C A B^{(1-r^2)/(r^2(1+2\tau))}.
\]
Since \(A\le1\), \(L\le 2\eta X^2\) and
\(\tau\le 2\eta/(r-c)^2\). Choose \(\eta\) sufficiently small
that both displayed exponents are at least \(\nu\).
Consequently \(|P_\sigma Q(t)|\le C A B^\nu\) for \(t\ge rX\).

## Decay by iteration

For \(x\ge X\), (TF1) and \(M(aX)\to0\) allow its factor of
\(|Q(x)|\) to be absorbed into \(r|Q(x)|\). Using
\(T_rQ(x)=P_\sigma Q(rx)\) gives
\[
 M(X)\le C M(aX)^{1-\nu}
       \left[M((b_0/r)X)^2+M((b_0/r)X)M(aX)^2\right]^\nu
       + C e^{-\kappa_3 X^2}
                                                               \tag{TF6}
\]
for some \(\kappa_3>0\). Here the terms from the exponential floors
in \(A,B\) are bounded by \(C e^{-\kappa_3 X^2}\), using boundedness
of \(M\) and positivity of \(\nu\).

Choose \(a,b_0,c\) close enough to \(r\) and \(\nu\) close
enough to \(1-r^2\) that
\[
  \theta=\min\{(1-\nu)a^2+2\nu(b_0/r)^2,\,
                  (1+\nu)a^2+\nu(b_0/r)^2\}>1 .
\]
The two expressions tend to
\(1+(1-r^2)^2\) and \(1+r^2(1-r^2)\), respectively.
Choose \(K>0\) with \(2C e^{-(\theta-1)K}\le 1/2\).
For sufficiently large \(X_0\),
\(M(x)\le e^{-K}\) on \(x\ge aX_0\). Therefore
\(M(x)\le e^{-K(x/X_0)^2}\) on \([aX_0,X_0]\).
Enlarge \(X_0\) so also
\(C e^{-\kappa_3 X^2}\le \frac12 e^{-K(X/X_0)^2}\)
for all \(X\ge X_0\). It suffices to take
\(K/X_0^2\le \kappa_3/2\) and \(C e^{-\kappa_3 X_0^2/2}\le 1/2\).
If the displayed bound on \(M\) holds at \(aX\) and
\((b_0/r)X\), (TF6) implies it at \(X\).
Induct over intervals with successive endpoints
\(X_0(b_0/r)^{-n}\) to obtain it for every \(X\ge X_0\).
Comparing \(|v|\) to \(|q|\) proves the claimed height bound.

## Exact straightness at a sufficiently fast Gaussian rate

Suppose in addition that \(h(x)-b=O(e^{-\lambda x^2})\) with
\(\lambda>0\) such that
\[
  \mu=\min\!\left\{2\lambda,\,
               \lambda+\frac{2\lambda r^2}{1+4\lambda(1-r^2)}\right\}
     \ \ge\ \frac{r^2}{2(1-r^2)}.                          \tag{TF7}
\]
Then the end is exactly straight on a half-line.

Use a fixed positive cutoff \(S_0>S\) for the transformed height
\(Q=v\,1_{\{x\ge S_0\}}\). Its global absolute bound is
\(C e^{-\lambda x^2}\) after enlarging \(C\) to cover the finite
initial interval. For sufficiently large \(x\), the derivation of
(TF1) with this cutoff bounds the error from the left by
\(C T_r 1_{(-\infty,S_0]}(x)\). Also
\[
 T_r(Q^2)(x)\le C\exp\!\left(
               -\frac{2\lambda r^2 x^2}{1+4\lambda(1-r^2)}\right)
\]
by Gaussian integration. Since the residual \(f=T_rQ-rQ\) is
bounded on the whole line, choosing \(S_1\ge S_0\) sufficiently
large gives
\[
  |f|\le C\left(1_{(-\infty,S_1]}+
           T_r 1_{(-\infty,S_0]}+ e^{-\mu x^2}\right)
\]
everywhere. In (TF2), the first term on the right contributes
at most \(C\exp[-(t/r-S_1)^2/2]\) for large positive \(t\).
For the second term use \(P_1T_r w(u)=P_1w(ru)\), since
the two variances sum to \(r^2+(1-r^2)=1\). Its contribution
is bounded by \(C\exp[-(t-S_0)^2/2]\).
The last contribution is bounded by
\(C\exp[-\mu t^2/(r^2(1+2\mu))]\), which by (TF7) is at most
\(C e^{-t^2/2}\). Thus \(P_1Q(t+S_0)=O(e^{-t^2/2})\) on the
positive half-line. The critical one-sided Gaussian deconvolution
theorem implies \(\operatorname{supp}Q\subseteq(-\infty,S_0]\).
It applies also at equality in (TF7), since the translation has
absorbed the shift in the Gaussian bound and the Fourier-side
hypothesis is the Gaussian multiplier of a bounded function.
Therefore \(v=0\) for \(x>S_0\), and continuity gives \(h=b\) there.
[Critical one-sided Gaussian deconvolution](GAUSSIAN_EXTERIOR_LINEAR_RIGIDITY.md#one-sided-gaussian-deconvolution)

## Improving the Gaussian rate by shifting the inversion contour

Fix \(R=r^2\) and \(q_0=1-R\). Start with any estimate
\(|h(x)-b|\le C e^{-\lambda x^2}\) on the tail, with \(\lambda>0\).
Such an estimate was proved from boundedness in (TF6). Put
\(\Lambda=R/(2q_0)\), and let \(\mu=\mu(\lambda)\) be the minimum
in (TF7). If \(\mu\ge\Lambda\), the preceding argument applies.
Assume \(\mu<\Lambda\) and use the same fixed cutoff \(Q\) as in
that argument. We will improve \(\lambda\) to
\[
 \lambda_{\rm new}
    =\min\!\left\{\mu,\frac{R c}{1-2R c}\right\},
 \quad
 a=\frac{\lambda}{1+2\lambda},\quad
 B=\frac{\mu}{R(1+2\mu)},\quad
 c=B-\frac{(B-a)(1/2-B)}{1/(2R)-B}.                         \tag{TF8}
\]
The right endpoint \(B<1/2\) follows from \(\mu<\Lambda\), and
\(0<a<B\) follows from \(\mu\ge\lambda\) and \(R<1\).

The same bounds on \(f=T_rQ-rQ\) and the same evaluation of
(TF2) as in the proof of (TF7) give
\[
  |G(t)|\le C e^{-B t^2}\qquad(t\ge0),\qquad G=P_1Q .
                                                               \tag{TF9}
\]
Indeed, the other two contributions to (TF2) are bounded by
shifted Gaussians of rate at least \(1/2\), so they are bounded
by the right side since \(B<1/2\). Enlarging \(C\) handles
the bounded initial interval. On the entire complex plane,
completing the square in the Gaussian integral with
\(|Q(x)|\le C e^{-\lambda x^2}\) gives
\[
            |G(s+iy)|\le C e^{-a s^2+y^2/2}.                \tag{TF10}
\]

Write \(d=2\sqrt{(B-a)(1/2-B)}\) and choose
\(\theta\in(0,\pi/2)\) by \(\tan\theta=\sqrt{(B-a)/(1/2-B)}\).
The entire function \(G(z)e^{(B+id/2)z^2}\) is bounded on
the positive real axis by (TF9), and on the ray of angle
\(\theta\) by (TF10). In fact on that ray
\((B-a)s^2+(1/2-B)y^2-dsy=0\).
The function has growth at most \(C e^{C|z|^2}\).
For the sector maximum principle, multiply it by
\(\exp[-\varepsilon(e^{-i\theta/2}z)^p]\) with
\(2<p<\pi/\theta\). Choose the power continuous on the
closed sector and positive on its central ray after the
rotation. Its real part is at least a positive constant
times \(|z|^p\) throughout the sector, so the maximum
principle on truncated sectors applies. Letting the outer
radius tend to infinity and then \(\varepsilon\downarrow0\)
makes the original function bounded in the sector.
Conjugating gives the lower sector bound as well. We obtain
\[
 |G(s+iy)|\le C\exp(-B s^2+B y^2+d s|y|)
       \qquad(s\ge0,\ |y|\le s\tan\theta).                 \tag{TF11}
\]

For real \(t,s\), (TF5) and Cauchy's theorem give the identity
\[
 P_\sigma Q(t)
   =\frac{1}{\sqrt{2\pi R}}
      \int_{\mathbb R}
          G(s+iy)\exp\!\left(\frac{(s+iy-t)^2}{2R}\right)dy .
                                                               \tag{TF12}
\]
The line with \(s=t\) is (TF5). To shift to any fixed \(s\),
integrate the entire integrand on the rectangle between the
two vertical lines. Its horizontal sides vanish as
\(|y|\to\infty\) by (TF10), since \(1/(2R)-1/2>0\).

For \(s\ge0\), take absolute values in (TF12). Factor out
\(\exp[(s-t)^2/(2R)]\). On \(|y|\le s\tan\theta\), (TF11)
and the Gaussian integral after completing the square give
the remaining bound \(C e^{-c s^2}\), because
\[
 -B s^2-\bigl(1/(2R)-B\bigr)y^2+d s|y|
     \ \le\ -c s^2 .
\]
On the complementary interval, (TF10) and the Gaussian
tail bound give
\[
 C\exp[-a s^2-(1/(2R)-1/2)s^2\tan^2\theta].
\]
At \(|y|=s\tan\theta\) the two pointwise exponents agree,
so this bound is at most \(C e^{-c s^2}\) as well. Hence
\[
 |P_\sigma Q(t)|\le C\exp[(s-t)^2/(2R)-c s^2]
        \qquad(s\ge0).
\]
Here \(0<c\le B<1/2\), as also follows from
\(c=[q_0B+R a(1-2B)]/(1-2RB)\). For \(t\ge0\) choose
\(s=t/(1-2R c)\). This yields
\[
 |P_\sigma Q(t)|\le C e^{-c t^2/(1-2R c)}.                \tag{TF13}
\]
At \(t=rx\), \(P_\sigma Q(t)=T_rQ(x)\). Since
\(rQ=T_rQ-f\), the residual bound used for (TF9) now
gives \(Q(x)=O(e^{-\lambda_{\rm new} x^2})\). The
left-tail terms have rate \(\Lambda>\mu\), so they do
not change this conclusion. This proves (TF8).

## Reaching the straightness threshold

The inequality \(R c/(1-2R c)>\lambda\) is equivalent
to \(R c>a\). Using the displayed expression for \(c\)
and the definitions of \(a,B\), this is equivalent to
\[
              \mu(1+2q_0\lambda)>(1+R)\lambda .
                                                               \tag{TF14}
\]
For the branch \(\mu=2\lambda\) it follows from \(R<1\).
For the other branch set \(t_0=2q_0\lambda\). After
division by \(\lambda\), the difference of the two
sides of (TF14) is \(t_0+R/(1+2t_0)>0\).
Also \(\mu>\lambda\), so every application of (TF8)
below the threshold strictly improves the exponent.

The function \(\mu(\lambda)\) is continuous and strictly
increasing from zero to infinity, and \(\mu(\lambda)>\lambda\)
for \(\lambda>0\). The threshold value
\(\lambda_*=\mu^{-1}(\Lambda)\) therefore satisfies
\(0<\lambda_*<\Lambda\). The update in (TF8) is continuous
on \((0,\lambda_*)\) and extends to \(\lambda_*\) with
value \(\Lambda\): as \(B\uparrow1/2\), \(c\to1/2\).
An increasing sequence of iterates staying below
\(\lambda_*\) would have a positive limit in
\((0,\lambda_*]\) fixed by the extended update, which
strictly increases every point in that interval. Thus
after finitely many updates (TF7) applies, and \(h=b\)
on a tail.

In particular, for a maximizer with a bounded isolated
graph for each of the three pairs, all three ends are
straight. The geometric argument of (PG10) gives
agreement with a centered simplex outside a compact
set. A nonsimplex discrepancy requires a bounded
strict-winner component in each label. If any one
strict-winner set is connected, or if any label has
no bounded strict-winner component, the entire
partition is therefore the simplex. First contacts
without these hypotheses remain in scope.
[Exterior geometry from straight ends](THREE_CELL_PRESSURE_GRADIENT_TAIL.md#exterior-geometry-when-all-three-ends-straighten)
[Bounded receiving components in the compact discrepancy](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)
