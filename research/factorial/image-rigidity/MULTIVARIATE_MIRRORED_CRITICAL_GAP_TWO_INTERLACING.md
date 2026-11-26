# All-degree interlacing for the shared mirrored pair \((1,3),(3,1)\)

For every integer \(n\ge1\), let
\[
 F_n(z)=\sum_{k=0}^n\binom nk(n+2k)!(3n-2k)!\,z^k .
\]
The roots of \(F_n\) and \(F_{n+1}\) strictly interlace on the
unit circle in the angular ordering cut at \(1\). In particular,
these consecutive polynomials are coprime at every order. This is
the ordinary fully shared mirrored binomial with exponent rows
\((1,3),(3,1)\), with arbitrary constant rows allowed.

The proof combines the exact factorial recurrence between orders with
the hypergeometric differential equation at the lower roots. The
required signs then follow from a mesh bound on two Gegenbauer factors.

## Exact signs in terms of the lower roots

Set \(N=n+1\), \(h=2n+1\), \(A=nh(h+1)/2\),
\(c=\cos(\theta/2)\), \(s=\sin(\theta/2)\), and
\(t=s/c\) where \(c\ne0\). Write
\[
 f(\theta)=e^{-in\theta/2}F_n(e^{i\theta}),\qquad
 g(\theta)=e^{-iN\theta/2}F_N(e^{i\theta}),\qquad
 \partial=\frac{d}{d\theta}.
\]
Both functions are real by reciprocity. With \(E=z\,d/dz\), the
successive coefficient ratio inside \(F_n\) gives
\[
 E(3n+2-2E)(3n+1-2E)F_n
   =z(n-E)(n+2+2E)(n+1+2E)F_n .
\]
Substitute \(E=n/2-i\partial\) in the centered expression. The
two differential operators are complex conjugates. The real and
imaginary parts of the first are respectively \(A+6N\partial^2\)
and \(4\partial^3-(3n+2)\partial\). The resulting angular equation is
\[
 c\big(4f'''-(3n+2)f'\big)=s\big(6Nf''+Af\big).                 \tag{1}
\]

There is a second identity between the two degrees. The coefficient
ratio from \(F_n\) to \((1-E/N)F_N\), at the same index \(k\le n\), is
\[
 (n+2k+1)(3n-2k+1)(3n-2k+2)(3n-2k+3).
\]
The other summand \(EF_N/N\) is its reversal in degree \(N\).
If \(q=(h^2+4\partial^2)f\), centering this decomposition gives
\[
 g=2c\big((h+1)(h+2)q-4q''\big)+4(2h+3)s q'.                \tag{2}
\]
Indeed the centered first summand is
\(e^{-i\theta/2}(h+1+2i\partial)(h+2+2i\partial)q\), and the
centered reversal is its complex conjugate.

At a zero of \(f\) with \(0<\theta<\pi\),
(1) gives \(4 f'''=(3n+2)f'+6Nt f''\).
Its derivative at that zero is
\[
 4 f''''=(3n+2)f''+6Nt f'''+3N(1+t^2) f''+At f'.
\]
Inserting both formulas into (2) yields the exact identity at the zero
\[
                  g=12c N^2 t\big((4n+3)f'+2t f''\big).       \tag{3}
\]
We will show that the factor in parentheses has the sign of \(f'\).

## A mesh bound rules out a wrong sign

By the
[factorial circle-root theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order),
the lower roots are simple and have arguments \(\theta_j\in(0,2\pi)\).
At a lower zero \(\theta_i\), set
\[
 r=\frac{f''(\theta_i)}{f'(\theta_i)},\qquad
 S=\sum_{j\ne i}\csc^2\!\left(\frac{\theta_i-\theta_j}{2}\right).
\]
The real centered factorization \(f(\theta)=K\prod_j
\sin((\theta-\theta_j)/2)\) gives
\[
 \frac{f'''(\theta_i)}{f'(\theta_i)}=\frac{3(r^2-S)-1}{4}.
\]
To see this, remove the vanishing sine from the product.
The first two logarithmic derivatives of the remaining factor at
\(\theta_i\) are \(r/2\) and \(-S/4\). The third derivative of the
removed sine divided by its first derivative is \(-1/4\).
Using (1) at the zero now gives
\[
                         S=r^2-2Nt r-N .                       \tag{4}
\]

Grouping the rising factorials in \(F_n\) by parity writes it, up to
a positive factor, as the degree-\(n\) Schur–Szegő composition of
two Gegenbauer circle polynomials of parameters \(a=(n+1)/2\) and
\(a+1/2\). The
[composition preserves the angular mesh of each factor](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#actual-mirrored-reference-profiles).
For the first factor, its centered polynomial is proportional to
\(C_n^a(\cos(\theta/2))\). Multiplication by \(\sin(\theta/2)^a\)
takes its equation to \(y''+V y=0\), where
\[
 V(\theta)=\frac{(n+a)^2}{4}
               -\frac{a(a-1)}{4\sin^2(\theta/2)}
       \ \le\ \frac L4,\qquad L=n^2+2an+a .
\]
Every interval between consecutive zeros of \(y\), including the
two intervals from an endpoint \(0\) or \(2\pi\) to the first or last
interior zero, has length at least \(2\pi/\sqrt{L}\). Integration by
parts gives \(\int |y'|^2\le (L/4)\int |y|^2\) on such an interval,
and the Dirichlet Poincaré inequality gives the length bound. At the
endpoints the integration applies because \(a\ge1\),
\(y=O(\theta^a)\), \(y'=O(\theta^{a-1})\) near zero, and the same
bounds hold with \(2\pi-\theta\) near \(2\pi\). The circular gap through
the cut is the sum of both endpoint intervals. Thus the angular mesh
of \(F_n\) is at least \(\mu=2\pi/\sqrt L\).

Count neighbors of a lower root in each direction up to angular
distance \(\pi\), assigning an antipodal neighbor to just one direction.
The \(j\)-th neighbor in either list has distance
\(\delta\ge j\mu\), so \(\sin(\delta/2)\ge\delta/\pi\) implies
\[
 S\le\frac{2\pi^2}{\mu^2}\sum_{j\ge1}j^{-2}
   =\frac{\pi^2 L}{12}
   <4n^2 .
\]
Here \(L=2n^2+\tfrac32 n+\tfrac12\le4n^2\) for \(n\ge1\).
If the factor in parentheses in (3) had zero or opposite sign
relative to \(f'\), then \(r\le -(4n+3)/(2t)\) since \(t>0\).
Equation (4) would give \(S\ge N(4n+2)>4n^2\), a contradiction.
This proves the sign agreement for \(0<\theta<\pi\), and the
reciprocity symmetries
\(g(2\pi-\theta)=(-1)^Ng(\theta)\) and
\(f'(2\pi-\theta)=(-1)^Nf'(\theta)\) give it for \(\pi<\theta<2\pi\).

## The central root and the cut

When \(f(\pi)=0\), differentiating (1) at \(\pi\) gives
\[
 (6N+2)f'''(\pi)=\big((3n+2)/2-A\big) f'(\pi).
\]
At the same point (2) gives \(g(\pi)=4(2h+3)(h^2 f'(\pi)+4 f'''(\pi))\).
The ratio \((6N+2)(h^2f'+4f''')/f'\) here is strictly positive, since
\[
 h^2(6N+2)+2(3n+2)-4A
   =h\big((4n+8)h-2n\big)+2(3n+2)>0 .
\]
Thus \(g\) has the sign of \(f'\) at every lower zero. Since \(f(0)>0\),
the signs of \(f'\) at its roots in ascending angular order are
\((-1)^i\), \(1\le i\le n\). The upper endpoint values satisfy
\(g(0)>0\) and \(g(2\pi)=(-1)^Ng(0)\).
Hence \(g\) has one zero in each of the \(n+1\) open intervals
between the lower zeros and the cut endpoints, as claimed.
