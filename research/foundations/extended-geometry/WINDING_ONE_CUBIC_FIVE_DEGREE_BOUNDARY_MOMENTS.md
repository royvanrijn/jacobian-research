# Fixed-deficit moment conditions for the degree-five cubic boundary

Let \(K\) have characteristic zero. Use the Laurent presentation
\[
 f=x^{-1}(C+Bt+Dt^2+Et^3),\qquad C,B,D,E\in K[x],
\]
of a regular winding-one phase on \(t^2+2xy=1\). The functional
\(\Lambda\) takes the angular constant term and applies normalized
integration over \([-1,1]\) in \(t\). Suppose that
\(E=e(x-r)^4(x-s)\), where \(e,r,s\) are nonzero and \(r\ne s\).
We allow scalar extension of the coefficient field to name the roots.

By the [logarithmic boundary theorem](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md),
eventual purity in this branch forces \(v^3\mid D\), \(v^2\mid B\),
and \(v\mid C\) for \(v=x-r\). Write
\[
 C+Bt+Dt^2+Et^3=vR(x,vt),\qquad
 R(x,T)=U_0+U_1T+U_2T^2+UT^3,
\]
where \(U=e(x-s)\), \(U_2=D/v^3\), \(U_1=B/v^2\), and
\(U_0=C/v\). In the degree-seven case
\(\deg U_j\le7-2j\) for \(0\le j\le2\).

## Proposition (Fixed-deficit rational conditions)

Define on \(K(x)\) the linear operator
\[
 \mathcal L(h)=h'-\frac h3\left(\frac2x+\frac1{x-r}+\frac1{x-s}\right).
\]
For each integer \(d\ge0\), put
\[
 \begin{split}
 G_d(x)&=[z^{3d}]\left(
    1+\frac{U_2}{U}z+\frac{U_1}{U}z^2+\frac{U_0}{U}z^3
                     \right)^{d-1/3},\\
 \phi_d(x)&=\frac{v^{d-1}U^d}{x^d}\,G_d(x).
 \end{split}                                                       \tag{1}
\]
The coefficient is defined by the formal binomial series. For each \(d\)
there is a unique reduction
\[
                \phi_d=\frac{a_d}{x-r}+b_d+\mathcal L(h_d),\qquad
                       h_d\in K(x).                                \tag{2}
\]
If \(\Lambda(f^m)=0\) for all sufficiently large \(m\), then
\(b_d=0\) for every \(d\ge0\). In particular the first condition on the
lower coefficients is the vanishing of the constant part in (2) for
\[
 \phi_1=\frac{2}{3x}
       \left(U_0-\frac{U_1U_2}{3U}+\frac{2U_2^3}{27U^2}\right).      \tag{3}
\]

### Reduction by rational derivatives

The finite poles of \(\phi_d\) are among \(0,r,s\). At any of these
points \(a\), the coefficient of \((x-a)^{-1}\) in the multiplier
\((2/x+1/(x-r)+1/(x-s))/3\) is \(\alpha=2/3\) or \(1/3\).
The leading term of \(\mathcal L(c(x-a)^{-n})\) for \(n\ge1\) is
\(-c(n+\alpha)(x-a)^{-n-1}\). Successively subtract these images to
remove all finite poles of order at least two. This introduces at most
simple poles at the other finite points. The leading term at infinity
of \(\mathcal L(x^n)\) is \((n-4/3)x^{n-1}\), so polynomial images
remove the polynomial part of positive degree. Finally use the
image of an affine polynomial to cancel the simple poles at \(0\)
and \(s\), prescribing its values at those two distinct points.
The resulting representative is \(a_d/v+b_d\).

For uniqueness, if \(a/v+b=\mathcal L(h)\), a finite pole of \(h\)
would give a pole of order at least two of \(\mathcal L(h)\).
This also holds at a finite point outside \(\{0,r,s\}\), where the
derivative supplies the leading pole. Thus \(h\) is a polynomial.
The infinity calculation gives \(\deg h\le1\), and its values at
\(0\) and \(s\) must vanish since no poles occur there. Hence \(h=0\).
The reduction only uses field operations with specified nonzero
denominators. Geometrically, with \(F_0=x^2vU\) and \(W^3=F_0\),
the operator identity is \(d(h/W)=\mathcal L(h)\,dx/W\).

### The selected prime coefficient

First work over a number field. Fix \(d\), and take sufficiently large
good primes \(p=3k+1\), with \(m=k+d\). The angular-constant polynomial
of \(f^m\) has \(t\)-height at most \(3m=p-1+3d\). For \(p>3d\) its
only height with an interval denominator divisible by \(p\) is
\(p-1\), which is even. Thus \(p\Lambda(f^m)\) reduces to
\[
       [x^{k+d}]v^{4k+d}[T^{3k}]R^{k+d}.
\]
In characteristic \(p\), \(v^{4k+d}=(x^p-r^p)v^{k+d-1}\).
Since \(k+d<p\) and the rest of the product is polynomial for
sufficiently large \(p\), only the constant \(-r^p\) contributes
from \(x^p-r^p\).

The deficit from the top height of \(R^{k+d}\) is \(3d\).
Its coefficient is \(U^{k+d}\) times the coefficient of \(z^{3d}\)
in the series in (1) with exponent \(k+d\). Each contributing
factorial denominator is prime to \(p\) when \(p>3d\), so this
coefficient reduces to \(G_d\) because \(k+d\equiv d-1/3\pmod p\).
The angular selection therefore gives, up to \(-r^p\),
\[
 [x^{k+d}]v^{k+d-1}U^{k+d}G_d
    =[x^k](vU)^k\phi_d
    =[x^{p-1}]F_0^k\phi_d.                                       \tag{4}
\]
All coefficient operations on rational functions use their Laurent
expansions at zero. The first line has \(G_0=1\) and \(\phi_0=1/v\).
For \(d=1\) the three contributions to \(G_1\) have respective
coefficients \(2/3,-2/9,4/81\) on
\(U_0/U,U_1U_2/U^2,U_2^3/U^3\), giving (3).

The retained coefficient of a rational derivative disappears in
characteristic \(p\). Indeed \(k\equiv-1/3\pmod p\) gives
\(F_0^k\mathcal L(h)=(F_0^k h)'\), whose coefficient of \(x^{p-1}\)
is zero. So (4) reduces to
\[
 a_d [x^k](vU)^k/v + b_d [x^k](vU)^k.                            \tag{5}
\]

### The nonzero coefficient and the ratio restriction

Put \(\mu=r/s\), so \(\mu\ne0,1\). In the scaled coordinate
\(x=rX\) define
\[
 H_p(\mu)=[X^k](1-X)^k(1-\mu X)^k
         =(-1)^k\sum_{j=0}^k\binom{k}{j}^{\!2}\mu^j .
\]
The first coefficient in (5) is a nonzero multiple of
\[
 [X^k](1-X)^{k-1}(1-\mu X)^k
                         =\frac{\mu}{k}H_p'(\mu).
\]
Thus \(d=0\), with \(a_0=1\) and \(b_0=0\), forces
\(H'_p(\mu)=0\). The coefficient recurrence of \(H_p\) gives
\[
 \mu(1-\mu)H_p''+\bigl(1+(2k-1)\mu\bigr)H_p'-k^2 H_p=0.
\]
At any \(\mu\ne0,1\) in the residue field, simultaneous vanishing of
\(H_p,H_p'\) would force all derivatives through order \(k\) to vanish
by this equation and its derivatives. Since \(k<p\), the polynomial
would be zero, contrary to its constant coefficient. Hence the second
coefficient in (5) is nonzero at every sufficiently large good prime.
Vanishing of (5) forces \(b_d=0\) modulo all such primes. This proves
the claimed rational identity. Over a general characteristic-zero field,
if \(b_d\) were nonzero, specialize a finitely generated coefficient
algebra to a number field retaining \(b_d\ne0\) and all required
nonvanishings. The moment identities persist, giving the same
contradiction. Absolute residue degree one is not needed in (4).

The [quadruple-root marker and ratio analysis](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_MARKERS.md)
gives the independent fixed-marker consequence and the exact known
restrictions on \(\mu\).
