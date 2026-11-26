# Shifted midpoint laws in the angular equations of both-private factorial moments

Reserve one factorial on each side of an ordinary binomial
with private variables in both monomials. The
[positive Gamma factorization](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)
gives a geometric moment polynomial with residual positive
bases \(U,V\). The following reduction applies to arbitrary
orientations of their ordinary shared rows.

At a proposed common root, the angular equation in each
degree is exactly the equation for an even reflection of
its shifted midpoint law. An independent odd equation
fixes a unique candidate radius for each angle in a
right-half angle box. The exact matched factorial kernels
also give a centered convex-order comparison across
degrees. A stronger translated likelihood-ratio comparison
would rule out a common root. It can fail
at an arbitrary physical radius even for a single ordinary
shared Gamma factor.

## The two exact equations

Let \(\rho_j\) be the law of \(H=\log(V/U)\) under the
normalized tilt by \((UV)^{j/2}\). At \(z=r e^{i\theta}\)
put \(S=H+\log r\), \(D_\theta(s)=\cosh s-\cos\theta\).
For \(0<\theta<\pi\), define
\[
 \begin{aligned}
 A_e&=\mathbb E_{\rho_j}\frac{\cosh((j/2+1)S)}{2D_\theta(S)},&
 B_e&=\mathbb E_{\rho_j}\frac{\cosh(jS/2)}{2D_\theta(S)},\\
 A_o&=\mathbb E_{\rho_j}\frac{\sinh((j/2+1)S)}{2D_\theta(S)},&
 B_o&=\mathbb E_{\rho_j}\frac{\sinh(jS/2)}{2D_\theta(S)} .
 \end{aligned}
\]
All four integrals are finite since their exponential
growth has rate at most \(j/2\).
Set \(K_e=A_e-B_e e^{-i\theta}\) and
\(K_o=A_o-B_o e^{-i\theta}\). After removal of a positive
scalar, the moment polynomial at \(z\) equals
\[
                 (K_e-K_o)+e^{ij\theta}\overline{(K_e+K_o)}.
\]
Indeed, split \(\sum_{k=0}^j w^k\) at \(w=e^{S+i\theta}\)
as \(1/(1-w)+w^j/(1-1/w)\). The midpoint tilt inserts
the factor \(e^{-jS/2}\) up to a positive scalar.
Multiplying by \(e^{-ij\theta/2}\) therefore turns the
two real equations at a root into
\[
  \Re(e^{-ij\theta/2}K_e)=0,\qquad
  \Im(e^{-ij\theta/2}K_o)=0.                                  \tag{1}
\]

In particular, the first equation is the symmetric
geometric resolvent equation for the law of \(|S|\).
At a root \(K_o/K_e\) is purely imaginary. Equivalently,
the product of the endpoint packets with the same
phase convention is
\(K_e^2-K_o^2=K_e^2(1-(K_o/K_e)^2)\), a strictly positive
multiple of \(K_e^2\). Thus their phase sum is exactly
\(2\arg K_e\). Their individual phases belong to
\((0,\pi-\theta)\), and the branch is fixed by the
[angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).

Suppose a common root were to occur at consecutive
orders \(n,N=n+1\). By the boxes its upper-half-plane
angle is in the same numbered box of both degrees and
lies below \(\pi/2\). The two first equations in (1)
use the same physical shift \(\log r\). If
\(\sigma_j(r)\) is the law of \(|H+\log r|\) under the
order-\(j\) midpoint tilt, then the comparison
\[
                   (n/N)\sigma_n(r)
                       \preceq_{\rm LR} \sigma_N(r)            \tag{2}
\]
at this shared radius would suffice to rule out the root.
Indeed the
[linked phase proof](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)
uses only this lower likelihood-ratio comparison of
the folded laws, endpoint log moment integrability,
and the same angle labels. The integrability holds
for the Gamma variables here, since the positive
endpoint shapes give strict exponential margins.

## A unique radius for the odd equation

Fix \(j\ge3\) and an angle \(\theta\le\pi/2\) in an
open positive angle box. Put
\[
     \lambda=\cos\theta+\sin\theta\cot(j\theta/2)<1,\qquad
     h(s)=\frac{\sinh((j/2+1)s)-\lambda\sinh(js/2)}{D_\theta(s)}.
\]
The second equation of (1) is exactly
\(\mathbb E_{\rho_j}h(H+\log r)=0\). In this equation,
\(\lambda\) depends on the angle and degree alone.

The function \(h\) is odd and strictly increasing.
For \(s>0\) write \(g(s)\) for its numerator and
\(a=j/2\). Then \(g>0\), and
\[
     g'(s)-a g(s)=\cosh((a+1)s)+a e^{-as}(e^{-s}-\lambda)
                       \ge 1-a s e^{-as}>0 .
\]
Here \(e^{-s}-\lambda\ge-s\).
Also \(D_\theta'(s)/D_\theta(s)\le1/\sin\theta\), by
\[
       D_\theta(s)^2-\sin^2\theta\,\sinh^2s
                         =(\cos\theta\cosh s-1)^2 .
\]
The angle box gives \(\theta>\pi/j\).
On the right half-circle, sine concavity on
\([0,\pi/3]\) yields \(j\sin\theta\ge3\sqrt3/2>5/2\).
Therefore \(h'(s)\ge(j/10)h(s)>0\) for \(s>0\).
Oddness gives strict increase also on the other half-line.

The function \(\mathbb E_{\rho_j}h(H+x)\) is continuous,
strictly increasing in \(x\), and tends to opposite
infinities at the two ends of the real axis. Continuity
uses domination on compact \(x\)-sets by a constant
times \(1+e^{j|H|/2}\), and the limits follow by monotone
convergence after subtracting a fixed integrable value.
Thus the odd equation specifies a unique radius
\(r_j(\theta)\). An angular contact of the shifted-even
equations can give a common ordinary root only when
\(r_n(\theta)=r_N(\theta)\).

## The exact centered transport and a limitation

For each Gamma variable in the positive representation,
let its shape before the tilt be \(\gamma>0\), and let
its two exponents have midpoint \(m\) and nonzero
signed gap \(d\). Under the order-\(j\) midpoint tilt
its shape is \(\gamma+jm\). The contribution of this
variable to the centered \(\log(V/U)\) has log-Gamma
Lévy density, in jump magnitude \(u>0\), given by
\[
            w_j(u)=\frac{e^{-(\gamma+jm)u/|d|}}
                              {u(1-e^{-u/|d|})},
\]
with jump sign opposite to \(d\). Indeed integrate
the trigamma series twice to get the centered
log-Gamma cumulant in the form
\[
 \log\frac{\Gamma(\alpha+z)}{\Gamma(\alpha)}
      -z\psi(\alpha)
   =\int_0^\infty (e^{-zu}-1+zu)
                    \frac{e^{-\alpha u}}{u(1-e^{-u})}\,du .
\]
One can split the integral into positive measures,
construct centered compound Poisson laws by
truncation, and pass to the limit using the integrable
second jump moment.

For \(t=n/N\), the densities satisfy
\(w_n(u)\ge w_N(u)\ge t^{-1}w_n(u/t)\).
The second comparison uses
\((\gamma+nm)/t=Nm+\gamma/t\) and the increase in
\(1-e^{-u/(t|d|)}\).
Sum separately for each jump sign. Consequently for
the centered midpoint laws \(H_j^\circ=H-\mathbb E_{\rho_j}H\)
there are independent centered noises in the identities
\(H_n^\circ\stackrel{\rm law}=H_N^\circ+Z\) and
\(H_N^\circ\stackrel{\rm law}=tH_n^\circ+Y\).
In particular
\[
              tH_n^\circ\preceq_{\rm cx}H_N^\circ
                               \preceq_{\rm cx}H_n^\circ,     \tag{3}
\]
where convex order means comparison of the expectations
of convex functions. The matched tanh variables recover
the exact midpoint scale \(m/|d|\), since for every
matched Gauss factor of an ordinary row
\[
              \frac{n}{T'}-\frac{n-1}{T}=\frac{2m}{|d|}.
\]
The Gauss offset \(h/|d|\) is the corresponding positive
shape offset in its Gamma representation.

Condition (3) is a centered comparison. In general it
cannot be upgraded to (2) at every common physical shift.
Take the ordinary shared row \((1,2)\) and the private
unit pair, and let \(n=2,N=3\), \(r=1/10\).
The midpoint shapes of \(X\) are \(4\) and \(11/2\),
and \(S=\log(rX)\). For shape \(a\), let \(f_a\) denote
the density of \(|S|\). The density before folding is
proportional to \(\exp(a s-10e^s)\). Thus
\[
           (\log f_a)''(0)=(a-10)^2-10 .
\]
The first derivatives of the folded densities vanish
at zero. In the ratio of the upper folded density to
the lower one dilated by \(2/3\), the second logarithmic
derivative at zero is
\[
                    \frac{41}{4}-\frac94\,26
                                          =-\frac{193}{4}<0 .
\]
The ratio decreases immediately to the right of zero.
The exact same-profile translated fold therefore fails
the monotonicity in (2) at this radius. This example
concerns an arbitrary radius. It does not assert a
failure of (2) at two coincident radial candidates,
which is the narrower condition relevant to common roots.
