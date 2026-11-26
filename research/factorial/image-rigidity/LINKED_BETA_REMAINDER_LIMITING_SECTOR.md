# The limiting sector and its fold correction for linked beta remainders

Fix positive integers \(a,b\), assume \(a\le b\), and use the
[fixed-degree linked-beta remainder](LINKED_BETA_FIXED_DEGREE_REMAINDER.md)
\(E_{n;a,b}(s)\). Put \(T=a+b\), \(d=T+1\), and use the native
coefficient variable \(z=s/(s-1)\). Define
\[
 R=\frac{a+1}{b+1},\qquad c=\frac{ab}{(a+1)(b+1)},\qquad
 \theta_*=2\arcsin\sqrt c.
\]
Every order-\(n\) linked moment root in this variable satisfies
\[
              a/b\le |z|\le R,\qquad
                 \sin^2(\arg z/2)>c
\]
by the [linked quadratic pencil](LINKED_BETA_VISHNYAKOVA_PENCIL.md).
In the balanced case the radius equals one.
The limiting remainder as \(n\to\infty\) has only its two fold
roots in the closed sector
\[
                  a/b\le |z|\le R,\qquad
                          \theta_*\le|\arg z|\le\pi .
                                                               \tag{1}
\]
These roots lie exactly at \(|z|=R\) and
\(\sin^2(\arg z/2)=c\). The finite-\(n\) angular correction
at either fold is strictly smaller than the first correction
in the exact moment-root sector bound. In particular this gives
a separate, fixed-\((a,b)\) sufficiently-large-order
coprimality proof. It still leaves the finite-order sector
comparison for general \(a,b\) open.

## Limiting divided difference and its sector

Write \(p_s(t)=(s-t)t^a(1-t)^b\). The characteristic roots of
the limiting recurrence for \(C_k\) in the fixed-degree
elimination are \(t_\pm\), the two solutions of
\[
                    dt^2-(Ts+a+1)t+as=0 .
\]
Thus \(C_k(\infty,s)=(t_+^k-t_-^k)/(t_+-t_-)\), with the
derivative interpretation at coincidence, and
\[
 E_{\infty;a,b}(s)=
                  \frac{p_s(t_+)-p_s(t_-)}{t_+-t_-}.
                                                               \tag{2}
\]
Equivalently, \(t_\pm\) are the nonendpoint stationary points
of \(p_s\). In the \(z\)-variable the discriminant, times
\((1-z)^2\), is
\[
  (b+1)^2 z^2-2(a+b+1-ab)z+(a+1)^2.
\]
Its roots are precisely the two boundary points
\(R e^{\pm i\theta_*}\).

Consider \(0<\arg z<\pi\) with \(|z|<R\). Then
\(\Im s<0\), and both \(t_\pm\) are in the lower
half-plane. Indeed their stationary parametrization is
\[
 s(t)=\frac{d}{T}t+\frac{q-1}{T}
                +\frac{q(1-q)}{T(q-t)},\qquad q=\frac aT .
\]
It maps each open half-plane to the same half-plane and
takes real values at real points where it is defined.
The identity
\[
                  s-t=\frac{t(1-t)}{a-Tt}
\]
allows the logarithms of the two stationary values to
be defined using principal logarithms of \(t,1-t,a-Tt\).
Label \(t_-\) by its limit zero as \(z\to0\).
Let \(D(z)\) be the logarithm of \(p_s(t_-)\) minus the
logarithm of \(p_s(t_+)\) with these choices. Since no
logarithm is singular at the fold, \(D\) continuously
takes the value zero there. Direct differentiation of
the stationary values gives
\[
 zD'(z)=(a+1)\frac{
          \sqrt{(1-e^{i\theta_*}Z)(1-e^{-i\theta_*}Z)}
                           }{1-RZ},\qquad Z=z/R,
                                                               \tag{3}
\]
where the square root is the product of the two principal
square roots, and is one at \(Z=0\). In fact
\((s-t_-)(s-t_+)=s(s-1)/d\), so the logarithmic
derivative in \(s\) is minus the square root of the
stationary discriminant divided by \(s(s-1)\).
The conversion \(ds/dz=-1/(z-1)^2\) gives (3).

The real part of the quotient in (3) is positive for
\(|Z|<1\). For \(R<1\), verify this on the boundary by
writing \(Z=e^{i\phi}\), \(0<\phi<\pi\).
The numerator has argument \(\phi/2\) when
\(\phi<\theta_*\) and \((\phi-\pi)/2\) when
\(\phi>\theta_*\).
The denominator has argument between
\((\phi-\pi)/2\) and zero. Thus the quotient has
nonnegative real part at the boundary. Conjugation covers
the other half, and the harmonic minimum principle gives
strict positivity inside. For \(R=1\) take the limit
\(R\uparrow1\) in this quotient and again use its positive
value at zero.

For \(R<1\), while \(\theta_*<\phi<\pi\) along
\(|z|=R\), the quotient in (3) has negative imaginary
part: its numerator argument is \((\phi-\pi)/2\) and
the denominator argument is strictly larger.
It follows that \(\Re D\) increases strictly from zero
along this outer arc. Along each fixed angle it increases
strictly with \(\log |z|\) by (3), starting at minus
infinity as \(z\to0\). Thus the locus \(\Re D=0\)
is one radial graph on this angular interval. It starts
at the fold and ends on the negative real ray.
On the graph \(\Im D\) increases strictly with the angle:
if \(zD'=U+iV\), then \(U>0\) and differentiation of
\(\Re D=0\) gives imaginary derivative \((U^2+V^2)/U\).
On the negative real ray \(0<s<1\), the two stationary
points are real, one in \((0,q)\) and one in \((q,1)\).
The chosen logarithms give \(\Im D=\pi\) there by
the sign of \(a-Tt\). Thus on the graph \(0<\Im D<\pi\),
and the two stationary values cannot agree.
For \(\phi=\theta_*\), \(\Re D<0\) radially inside.
When \(R=1\), the outer arc has \(\Re D=0\) throughout,
the interior has \(\Re D<0\), and the same differentiation
of the imaginary part on the outer arc gives the
strict increase from zero to \(\pi\).
This proves (1) using (2), with conjugation for the
other half-plane.

## Exact first correction at a fold

Let \(s_*\) be the upper-half-plane fold in the \(s\)-variable
and let \(u\) be its repeated saddle. Write
\(\kappa=a-Tu\), \(w(t)=t^a(1-t)^b\) and
\(f(t)=w(t)/(t(1-t))\). The stationary discriminant has
derivative \(-4d\kappa\) at the fold, and
\[
                p_{s_*}'(t)=d f(t)(t-u)^2 .
\]
Expanding (2) about its midpoint gives
\[
 E_\infty'(s_*)=w'(u)+
       p_{s_*}'''(u)\frac{-4d\kappa}{24d^2}
                 =\frac23 w'(u)\ne0.
\]
Here \(p_{s_*}'''(u)=2d f(u)\) and \(w'(u)=\kappa f(u)\).
The fold remainder root is therefore simple.

For the first order in \(\lambda=1/n\), write the
recurrence at this fold as
\[
 C_k(\lambda,s_*)=k u^{k-1}
                      +\lambda u^{k-1}\eta_k+O(\lambda^2).
\]
Put \(A=(s_*+1)/u\), \(B=s_*/u^2\).
The second forward difference of \(\eta_k\), multiplied
by \(d\), is
\[
     (k+2)(k+1)A-(k+1)kB-(k+3)(k+2).
\]
The initial values \(\eta_0=\eta_1=0\) give
\[
 \eta_k=\frac{A-B-1}{12d}k^{\underline4}
        +\frac{4A-2B-6}{6d}k^{\underline3}
        +\frac{A-3}{d}k^{\underline2}.
\]
Inserting this in the linear combination defining the
remainder replaces \(k^{\underline j}u^{k-1}\) by
\(u^{j-1}p_{s_*}^{(j)}(u)\). Use
\(p_{s_*}''(u)=0\),
\(p_{s_*}''''(u)=6d f'(u)\), and
\(w'(u)/w(u)=1/(s_*-u)\). The first remainder correction
is
\[
       \partial_\lambda E(0,s_*)
                =f(u)\frac{2u-3u^2-s_*+2u s_*}{6}.
\]
The derivative of the simple root with respect to
\(\lambda\) is therefore
\[
                s_1=-\frac{2u-3u^2-s_*+2u s_*}{4\kappa}.
\]
For the native fold root on the lower half-plane, the
logarithmic derivative with respect to \(\lambda\) is
\(s_1/[s_*(1-s_*)]\). Substituting
\(s_*=d u^2/a\) and \(u=(a+i\sqrt{ab/d})/T\) gives
\[
 \frac{s_1}{s_*(1-s_*)}
   =\frac{a-b}{2(a+1)(b+1)}
     -i\frac{ab(T+4)-Td}{4ab(a+1)(b+1)}
                                       \sqrt{\frac{ab}{d}} .
\]
Consequently the first-order increase of
\(\sin^2(\arg z/2)\) at the fold remainder root is
\[
                  c_E=\frac{ab(T+4)-Td}{4(a+1)^2(b+1)^2}.
\]
The exact root sector lower bound of the linked pencil
is at least
\[
 \frac{b(na+1)}{(a+1)(n(b+1)-1)}
       =c+\frac{1}{n}\frac{bd}{(a+1)(b+1)^2}+O(n^{-2}).
\]
Its first correction is strictly larger than \(c_E\),
since \(4 b d(a+1)>ab(T+4)-Td\).
The finite-\(n\) remainder roots near the two folds
are therefore outside the allowable root sector for all
sufficiently large \(n\) at the fixed exponent pair.
Equation (1) and compact convergence of the fixed-degree
remainder exclude all other possible roots there.
