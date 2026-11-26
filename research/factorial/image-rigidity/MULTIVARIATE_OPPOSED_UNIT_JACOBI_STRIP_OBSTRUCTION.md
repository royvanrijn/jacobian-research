# A linked remainder with zeros inside the Jacobi strip

The [Jacobi strip criterion](MULTIVARIATE_OPPOSED_UNIT_JACOBI_STRIP.md)
cannot exclude all common roots for arbitrary opposed unit rows
using the real-part strip alone. Take rows \((22,23),(2,1)\) at
orders \(n=23\), \(N=24\). Then the linked elimination polynomial
\(H_{23;22,1}\) has a conjugate pair strictly inside the strip
\[
                -\frac{37926}{38503}<\Re x<-\frac{21}{25}.
\]
The two roots supplied by this certificate are outside the known
moment-root annulus in the native coefficient variable.

## Exact disk certificate

In the notation of the strip criterion, \(T=554\), \(D=504\),
\(S=602\) at the upper order. These give the two displayed strip
endpoints. Let \(H=H_{23;22,1}\), and use the center and radius
\[
       c=\frac{-4209+3898i}{5000},\qquad r=\frac1{100000}.
\]
The closed disk lies strictly inside the strip. If
\(H(c+w)=\sum_j h_j w^j\), the exact check below gives
\[
 \bigl(|\Re h_0|+|\Im h_0|\bigr)
 +\sum_{j\ge2}\bigl(|\Re h_j|+|\Im h_j|\bigr)r^j
       <\frac r2\max\{|\Re h_1|,|\Im h_1|\}.
\]
Hence \(H(c+w)\) and \(h_1w\) have the same zero count in
\(|w|<r\) by Rouché's theorem. Conjugation gives the second zero in
the disjoint reflected disk.

The [paired beta pencil](LINKED_BETA_VISHNYAKOVA_PENCIL.md#2-root-radius-and-a-pairwise-obstruction-certificate)
gives \(|z|<2/23\) for the roots of the actual order-24 polynomial
with these rows. In the Cayley coordinate
\(x=(z+1)/(z-1)\), the disk \(|z|\le R<1\) maps to a disk of
radius \(2R/(1-R^2)\) centered on the real axis. With \(R=2/23\)
that radius is less than \(3/4\). Both certified disks have
\(|\Im x|>3/4\), so these zeros of the linked remainder cannot be
common moment roots.

~~~python
import sympy as sp

a, b, n = 22, 1, 23
N = n + 1
x = sp.Symbol("x")
T, D = N*(a+b)+2, N*(a-b)
S = T+2*N
low = -sp.Rational(D*S, (T+2)*T)
high = -sp.Rational(D, S-2)
assert (low, high) == (-sp.Rational(37926, 38503), -sp.Rational(21, 25))

p = sp.Poly((x-1)**a*(x+1)**b, x)
ratios = [sp.Poly(1, x), sp.Poly(-((T+2*n)*x+D)/sp.Integer(T+n-1), x)]
for k in range(2, min(n, a+b)+1):
    j = n-k+1
    first = (x*x-1)*(j+1)*ratios[-2]
    second = ((T+2*j)*x+D)*ratios[-1]
    ratios.append(sp.Poly(-(first+second)/(T+j-1), x))
H = sum(p.diff((0, k))*ratio/sp.factorial(k)
        for k, ratio in enumerate(ratios))
H = sp.Poly(H, x).clear_denoms()[1].primitive()[1]
assert H.degree() == 23

center = (-4209+3898*sp.I)/5000
radius = sp.Rational(1, 100000)
assert low < sp.re(center)-radius < sp.re(center)+radius < high
assert sp.im(center)-radius > sp.Rational(3, 4)
R = sp.Rational(2, 23)
assert 2*R/(1-R*R) < sp.Rational(3, 4)

shifted = [sp.S.Zero]*(H.degree()+1)
for coefficient in H.all_coeffs():
    for j in range(len(shifted)-1, 0, -1):
        shifted[j] = (center*shifted[j]+shifted[j-1]).expand()
    shifted[0] = (center*shifted[0]+coefficient).expand()

def taxicab(z):
    real, imag = z.as_real_imag()
    return abs(real)+abs(imag)

linear_reserve = radius*max(map(abs, shifted[1].as_real_imag()))
remainder = taxicab(shifted[0]) + sum(
    taxicab(value)*radius**j for j, value in enumerate(shifted[2:], 2))
assert remainder < linear_reserve/2
print("Two nonreal remainder zeros in the strip and outside the physical annulus")
~~~
