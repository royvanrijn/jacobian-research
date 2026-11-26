# All-degree interlacing for a private unit pair with mirrored core curvature four

Take an ordinary private unit pair \((1,0),(0,1)\) and a finite
collection of ordinary shared mirrored pairs
\((\ell_i+d_i,\ell_i),(\ell_i,\ell_i+d_i)\) with \(\ell_i,d_i\ge1\).
Set
\[
                  S=2\sum_i\frac{d_i^2}{\ell_i+d_i}.
\]
At each degree \(m\), one may additionally multiply the normalized
coefficients by \(\exp[-\tau k(m-k)/(2m)]\) for a common \(\tau\ge0\).
If \(S+\tau\ge4\), every consecutive pair strictly interlaces in
symmetric angular order. At \(\tau=0\) this gives SFC(2) for the
indicated ordinary profile over any characteristic-zero field.
Equal exponent rows affect only scalar factors.

The proof uses the exact private-unit resolvent with a finite
geometric expansion. This improves the previous curvature
threshold \(21/4\) and uses no individual-degree root checks.
[Previous angle boxes and curvature theorem](MULTIVARIATE_ONE_PRIVATE_MIRRORED_ANGLE_BOXES.md)

## The angular comparison at a single benchmark

Integration over the private unit pair cancels the binomial
coefficient up to the scalar \(m!\). Normalize the remaining
polynomial \(P_m(z)=\sum_{k=0}^m a_{m,k}z^k\) by
\(a_{m,0}=a_{m,m}=1\). Without the Gaussian factor,
\(a_{m,k}=\mathbb E_m T^k\), where \(T=V/U\) for the shared core
and \(\mathbb E_m\) means expectation under the endpoint tilt
\(U^m/\mathbb E[U^m]\). This law satisfies
\[
                       \mathbb E_m[T^m f(T)]=\mathbb E_m[f(1/T)] .
\]
For the Gaussian factor multiply \(T\) by an independent
log-normal variable with log mean \(-\tau/2\) and log variance
\(\tau/m\), or by one when \(\tau=0\). Its moments give the
specified coefficient multiplier. The same exchange identity
holds, since tilting the log-normal variable by its \(m\)-th
power reflects its log, with normalization one. Thus it holds
for the product.

For \(z=e^{i\theta}\) with \(0<\theta<\pi\), put
\[
 K_m(z)=\mathbb E_m(1-zT)^{-1},\quad
 c=\cos\theta,\quad D_\theta(T)=1-2cT+T^2,\quad
 A=\mathbb E_m D_\theta(T)^{-1},\quad B=\mathbb E_m T/D_\theta(T).
\]
The geometric sum and exchange identity give
\[
 P_m(z)=K_m(z)+z^m\overline{K_m(z)},\qquad K_m(z)=A-Be^{-i\theta}.
\]
The angle boxes from this resolvent are
\[
 \frac{(2i-1)\pi}{m}<\theta_{m,i}<\frac{2i\pi}{m+1},
      \qquad 1\le i\le\lfloor m/2\rfloor .
 \tag{1}
\]
For completeness, \(A-B\) equals
\(\frac12\mathbb E_m[(1-T)(1-T^{m+1})/D_\theta(T)]>0\)
by exchange. The strict inequality holds since \(S+\tau\ge4\)
makes the law nonconstant. Thus \(\phi_m=\arg K_m\)
belongs to \((0,(\pi-\theta)/2)\). The two endpoints in
each interval (1) give opposite signs to
\(2|K_m|\cos(m\theta/2-\phi_m)\). Reflection and the structural
root \(-1\) in odd degree account for the full degree. This
also assigns each positive root its phase label
\(m\theta_{m,i}/2-\phi_m(\theta_{m,i})=(i-1/2)\pi\).

For consecutive degrees \(m-1,m\), (1) already gives
\(\theta_{m-1,i}<\theta_{m,i+1}\) whenever both positive roots
exist. It also gives \(\theta_{m,i}<\theta_{m-1,i}\) when
\(\theta=(2i-1)\pi/(m-1)\ge\pi/2\), since then \(4i\ge m+1\).
For the remaining indices it suffices to show
\(\phi_m(\theta)<\theta/2\) at that benchmark. By continuity the
upper phase then reaches its \(i\)-th label below the benchmark,
which lies below the lower root. In terms of the resolvent,
\[
             \phi_m(\theta)<\theta/2
             \quad\Longleftrightarrow\quad (1+2c)B<A.
 \tag{2}
\]
This follows by taking the imaginary part of
\(e^{-i\theta/2}(A-Be^{-i\theta})\) and using
\(\sin(3\theta/2)/\sin(\theta/2)=1+2c\).

## Bounding the resolvent with endpoint moments

The normalized coefficient attenuation is
\[
             a_{m,k}\le \exp\bigl(-2k(1-k/m)\bigr).
 \tag{3}
\]
The factorial curvature bound gives its exponent with \(S/2\)
in place of \(2\). The Gaussian factor adds \(\tau/2\).
[Factorial endpoint curvature estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)

For any integer \(1\le L\le m\), the exact finite expansion is
\[
 K_m(z)=\sum_{k=0}^{L-1}a_{m,k}z^k+
                    z^L\mathbb E_m[T^L/(1-zT)] .
\]
On \(0<\theta\le\pi/2\) one has \(|1-zT|\ge\sin\theta\).
Taking the imaginary part and dividing by \(\sin\theta\) gives
\[
 B\le\sum_{k=1}^{L-1}
           \frac{\sin(k\theta)}{\sin\theta}\,a_{m,k}
                     +\frac{a_{m,L}}{\sin^2\theta}.
 \tag{4}
\]
Also \(D_\theta^{-1}\ge 2-D_\theta\) for positive \(D_\theta\),
so \(A\ge1+2c a_{m,1}-a_{m,2}\).

We use the elementary inequality
\[
                (1+2\cos\theta)/\sin^2\theta\le3/\theta^2
                        \qquad (0<\theta\le\pi/2).                \tag{5}
\]
Put \(x=\tan(\theta/2)\) and compare \(\arctan x\) with
\(x/\sqrt{(1+x^2)(1-x^2/3)}\). The latter has derivative
\((1+x^4/3)/[(1+x^2)(1-x^2/3)]^{3/2}\), at least \(1/(1+x^2)\)
because \((1+x^2)(1-x^2/3)^3\le1\) by the arithmetic-geometric
mean inequality. Both functions vanish at zero, and
substitution gives (5).

For \(5\le m\le8\), take \(L=2\) in (4). At the benchmark,
\(\theta\ge\pi/(m-1)\). Using \(3/\pi^2<61/200\) in (5),
\[
 (1+2c)B-A
 \le a_{m,1}+\bigl(1+(61/200)(m-1)^2\bigr)a_{m,2}-1<0 .
\]
In fact \(a_{m,1}\le e^{-8/5}<101/500\). For
\(b=61/200\), the function
\((1+b(m-1)^2)e^{-4+8/m}\) increases for \(m\ge5\).
The sign of its derivative is the sign of
\(2b(m-1)(m-2)^2-8\). Hence the second term is less than
\((1+49b)/20\), using \(e^{-3}<1/20\) at \(m=8\).
These bounds sum to \(3997/4000<1\).

For \(9\le m\le13\), take \(L=3\) in (4) and use
\(\sin(2\theta)/\sin\theta=2c\). Equation (5) gives
\[
 (1+2c)B-A
 \le a_{m,1}+7a_{m,2}+(61/200)(m-1)^2 a_{m,3}-1<0 .
\]
Here \(a_{m,1}\le e^{-16/9}<17/100\) and
\(a_{m,2}\le e^{-28/9}<9/200\). The function
\((m-1)^2e^{-6+18/m}\) increases for \(m\ge9\), as follows
from its logarithmic derivative. At \(m=13\) its value is
less than \(144/100\), since \(e^{-60/13}<1/100\).
The three positive terms are bounded by \(17/100\),
\(63/200\), and \(44/100\), respectively.

For \(m\ge14\), take \(L=\lfloor m/2\rfloor\) and use
\(|\sin(k\theta)|\le k\sin\theta\) in (4). Since
\(\theta\ge\pi/m\) and \(\theta\le\pi/2\), one has
\(\sin\theta>3/m\), by the cubic Taylor lower bound for sine.
The first six summands \(k a_{m,k}\), evaluated at the worst
curvature envelope \(m=14\), are bounded above by
\(160,66,27,14,9,7\) divided by \(1000\).
For \(7\le k<L\), use \(a_{m,k}\le e^{-k}\).
The sum \(\sum_{k\ge7}k e^{-k}\) is less than \(12/1000\),
by summing the geometric derivative at \(e^{-1}<37/100\).
The remainder satisfies
\[
 \frac{m^2}{9}a_{m,L}
    \le \frac{(2L+1)^2}{9}e^{-L}\le25e^{-7}<24/1000.
\]
The middle inequality follows from decrease in \(L\ge7\).
Consequently \(B<319/1000\), whereas
\(A\ge1-a_{m,2}>967/1000\), which proves \(3B<A\).
The exponential bounds in all three cases follow from positive
Taylor sums for \(e^x\) through degree twenty. This completes (2)
at every required benchmark with \(m\ge5\).

For \(m\le4\) use the coefficient identities directly.
The quadratic \(1+a z+z^2\) has positive-half cosine
\(-a/2<0\), and the cubic \(1+bz+bz^2+z^3\) has positive-half
cosine \((1-b)/2\ge0\). The structural roots give interlacing
through degree three. For the cubic \(1+a z+a z^2+z^3\) and
quartic \(1+bz+c z^2+bz^3+z^4\), one has \(c\le b\) by
coefficient log-convexity and symmetry, and \(b\le a\) by the
adjacent first-coefficient identity. Indeed, a row pair
contributes
\[
                    \prod_{h=1}^{d_i}
       \frac{m\ell_i+h}{m(\ell_i+d_i)-d_i+h}
\]
to \(a_{m,1}\), and every displayed factor decreases with \(m\).
The additional Gaussian factor at \(k=1\) also decreases.
At the cubic positive-half cosine \(x=(1-a)/2\), the quartic
cosine quadratic \(4x^2+2bx+c-2\) is at most
\((1-a)^2+a(1-a)+a-2=-1\). This places the cubic root
strictly between the two quartic positive-half angles. The
structural root and reflection finish the small degrees.

At \(\tau=0\) the normalized polynomials are rational.
Disjoint root arguments give coprimality over \(\mathbb Q\),
and a rational Bézout identity gives coprimality after scalar
extension. A nonzero coefficient ratio therefore cannot
annihilate two consecutive factorial moments. The cases of a
single nonzero coefficient have nonzero factorial moments.
