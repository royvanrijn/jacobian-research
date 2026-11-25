# A global derivative sign for the half-line cubic

For real \(T\), let
\[
              F(T,H)=\int_0^\infty e^{-u^3/3-Tu^2/2+Hu}\,du.
\]
**Theorem.** At every zero with \(\Im H>0\), one has
\[
                              \Im F_H(T,H)<0.             \tag{1}
\]
In particular, for every real \(T\), all zeros are simple
and distinct zeros have distinct imaginary heights.
On an upper simple zero branch, \(\Im H_T>0\).
For \(T\le0\) every zero also satisfies
\(2H-T H_T\ne0\).

The proof uses a positive stationary density on the real
projective line. The Riccati and Fourier framework is
classical in random Schrödinger theory. The construction
below gives the required positivity directly for the
combination of Gaussian and Cauchy translations.
[Comtet and Tourigny, sections 4.4 and 6.2](https://arxiv.org/abs/1601.01822)

## A positive density for the normalized Airy solution

Fix real \(X\), \(Y>0\), and write \(Z=X+iY\).
Let \(v_Z\) be the solution on \(r\ge0\) of
\[
                       v_Z''+(Z+i r)v_Z=0
\]
that is bounded at positive infinity and has \(v_Z(0)=1\).
We show that its Hermitian extension to the real line
is the Fourier transform of a strictly positive density.
The convention used below is
\(\widehat\mu(r)=\int e^{-i r x}\,d\mu(x)\).

Use the real projective line with angle coordinate
\(\theta\bmod\pi\), where \(x=\tan\theta\) on the finite
chart. The drift \(x^2-X\) and the translation vector
field of \(x\) extend to smooth vector fields
\[
 V_0=(\sin^2\theta-X\cos^2\theta)\partial_\theta,\qquad
 V=\cos^2\theta\,\partial_\theta .
\]
Write \(\tau_s\) for translation \(x\mapsto x+s\),
extended to fix the point at infinity. For \(h>0\),
compose three Markov operators on this compact circle:
the time-\(h\) deterministic \(V_0\) flow, translation
by a centered Gaussian of variance \(2h\), and
translation by a centered Cauchy variable of scale
\(Yh\). The order of this composition can be fixed
arbitrarily. Each operator maps continuous functions
to continuous functions by dominated convergence, and
the composition has an invariant probability \(\mu_h\).
For completeness, use any initial probability, average
its first \(m\) iterates, and take a weak limit as
\(m\to\infty\). The difference between the transformed
and original averages tends to zero against every
continuous function.

For a smooth function \(f\), the increments of the
three operators, divided by \(h\), converge uniformly
to \(V_0 f\), \(V^2 f\), and
\[
 {\cal C}f
 =\frac{Y}{\pi}\int_0^\infty
          \frac{f\circ\tau_s+f\circ\tau_{-s}-2f}{s^2}\,ds
                                                               \tag{2}
\]
respectively. For the Gaussian term use Taylor expansion
of \(f(\tau_s(\theta))\), whose translation derivatives
are uniformly bounded on the compact circle. For the
Cauchy term its exact difference quotient has the same
numerator in (2) with denominator \(s^2+Y^2h^2\).
The numerator is bounded by \(s^2\|V^2f\|_\infty\)
at small \(s\), and by \(4\|f\|_\infty\) at large \(s\).
This proves uniform convergence. The three operators
also converge uniformly to the identity on each
continuous function. Expanding the difference of their
composition shows that its generator limit is the
sum of the three generators.

Pass to a weak limit \(\mu\) along \(h\downarrow0\).
Invariance of each \(\mu_h\) and the uniform generator
limit give
\[
             \int (V_0 f+V^2 f+{\cal C}f)\,d\mu=0        \tag{3}
\]
for every smooth \(f\) on the projective line. The measure
has no atom at infinity. To see this, choose a smooth
compactly supported function \(\chi\) on \((-1,1)\)
with \(\chi'(0)=1\), and use the tests
\[
                  f_\delta(\theta)=
                    \delta\chi((\theta-\pi/2)/\delta)
\]
on the local chart, extended by zero.
Their sup norms are \(O(\delta)\). On their supports
the coefficient of \(V\) is \(O(\delta^2)\) and its
derivative is \(O(\delta)\), giving
\(\|V^2 f_\delta\|_\infty=O(\delta^3)\).
Splitting (2) at \(s=1\) gives
\(\|{\cal C} f_\delta\|_\infty=O(\delta)\).
The drift term stays bounded and tends pointwise to
the indicator of the point at infinity, because
its coefficient there is one. Dominated convergence
in (3) proves the assertion.

Thus \(\mu\) is a probability on the finite real chart.
Every Schwartz test function on this chart extends
smoothly by zero at infinity. In that chart,
\(V_0=(x^2-X)\partial_x\) and \(V^2=\partial_x^2\).
Equation (2) has Fourier multiplier \(-Y|r|\), using
\(\int_0^\infty(1-\cos s)s^{-2}\,ds=\pi/2\).
For example the latter integral follows by integrating
by parts and taking the Abel limit of
\(\int_0^\infty e^{-\varepsilon s}\sin(s)s^{-1}\,ds
=\arctan(1/\varepsilon)\).
The symmetric second difference in (2) justifies the
Fourier interchange on Schwartz tests.
Writing \(\psi=\widehat\mu\), the adjoint of (3)
in Fourier coordinates is the distributional identity
\[
     i r\,(\psi''+X\psi)-(r^2+Y|r|)\psi=0 .
\]
The product \(x^2\mu\) causes no moment assumption,
as it is a tempered distribution. On \(r>0\) this
equation gives the ordinary differential equation
\[
                       \psi''+(Z+i r)\psi=0 .
\]
Indeed its right-hand expression for \(\psi''\) is
continuous there since \(\psi\) is a characteristic
function.

Let \(\omega=e^{2\pi i/3}\), \(\lambda=\omega^2 Z\),
and \(\alpha=e^{-i\pi/6}\).
Two solutions are \(\operatorname{Ai}(\lambda+\alpha r)\)
and \(\operatorname{Bi}(\lambda+\alpha r)\), since
\(\alpha^2\lambda=-Z\) and \(\alpha^3=-i\).
In the sector near angle \(-\pi/6\) the former decays
exponentially in \(r^{3/2}\) and the latter grows.
The bounded solution with value one at zero is therefore
\[
             v_Z(r)=
                \frac{\operatorname{Ai}(\lambda+\alpha r)}
                     {\operatorname{Ai}(\lambda)}
                          \quad (r\ge0).
\]
The Airy denominator is nonzero. The zeros of Ai are
negative real, and \(\omega^2 Z<0\) real would place
\(Z\) in the lower half-plane.
[Airy zero law](https://dlmf.nist.gov/9.9)
and [sectorial Airy asymptotics](https://dlmf.nist.gov/9.7)

Since \(\psi(0)=1\) and \(|\psi|\le1\), its positive
restriction is \(v_Z\). Its negative restriction is
the conjugate reflection. Thus \(\psi\) and all its
polynomially weighted functions are integrable. Fourier
inversion gives a smooth real density \(k_Z\ge0\),
integrating to one, with
\[
               \pi k_Z(x)=\Re \int_0^\infty
                                         e^{i x r}v_Z(r)\,dr .
                                                               \tag{4}
\]
This density is strictly positive everywhere. Indeed its
stationary equation holds pointwise as
\[
              k_Z''-\partial_x((x^2-X)k_Z)+{\cal C}k_Z=0 .
\]
The integral for the last term is continuous and equals
the Fourier multiplier just used, by the integrability
of \(|r|\psi(r)\). At a hypothetical zero of this
nonnegative smooth density, the derivative is zero
and the second derivative is nonnegative.
The drift divergence vanishes there, whereas
\({\cal C}k_Z>0\) there by (2), since \(Y>0\) and the
density is positive somewhere. This contradicts the
stationary equation and proves (4) with strict positivity.

## The Wronskian at a cubic zero

Let \(H\) be an upper zero of \(F(T,\cdot)\) and put
\(a=T/2\), \(Z=H+a^2\), \(X=\Re Z\), \(Y=\Im Z\).
The field equation
\(F_{HH}+T F_H-HF=1\) gives
\[
 \bigl(e^{a h}F(T,h)\bigr)''
       -(h+a^2)e^{a h}F(T,h)=e^{a h}.
\]
Set \(G(r)=e^{a(H+i r)}F(T,H+i r)\) for \(r\ge0\).
This solves
\[
              G''+(Z+i r)G=-e^{aH}e^{i a r},
                        \qquad G(0)=0 .
\]
The Wronskian satisfies
\[
          (v_Z G'-v_Z' G)'=-e^{aH}e^{i a r}v_Z(r).
\]
It tends to zero at infinity. The Airy solution and
its derivative decay there, and \(G,G'\) are bounded
by the absolutely convergent Fourier integrals on the
fixed vertical line. Since \(G'(0)=i e^{aH}F_H(T,H)\),
integration gives
\[
             i F_H(T,H)=\int_0^\infty e^{i a r}v_Z(r)\,dr .
\]
Equation (4) now gives the exact sign identity
\(-\Im F_H(T,H)=\pi k_Z(a)>0\).

Real fields have positive \(F\), and the lower zeros
follow by conjugation, so this proves simplicity.
Any upper same-height pair would have a left member
with \(\Im F_H>0\), by the
[finite-interval resolvent comparison](MULTIVARIATE_HALF_LINE_CUBIC_HORIZONTAL_PAIR_SIGNS.md).
The two signs exclude the pair at every real temperature.
The heat equation and the field equation at a zero give
\(H_T=(1/F_H-T)/2\), whose imaginary part is now strictly
positive in the upper half-plane. If \(T\le0\) then
\(\Im(2H-T H_T)>0\), proving the last assertion.

These results remove the distinct-branch height
obstruction from the
[bounded moving cubic window](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
A possible simple zero with \(2H-T H_T=0\) at positive
temperature needs a further estimate.
