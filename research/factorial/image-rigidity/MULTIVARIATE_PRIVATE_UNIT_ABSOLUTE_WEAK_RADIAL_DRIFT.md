# A cutoff-free weak-curvature bound for two private unit rows

Consider an ordinary two-monomial profile with one private unit
row on each side and arbitrary shared rows. Put
\[
 V=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\min(a_i,b_i)},
 \qquad \Delta=\sum_{\rm shared}(b_i-a_i).
\]

**Theorem.** If \(V\le10^{-11}\), the ordinary consecutive
moments are coprime in every degree. There is no restriction on
the number, orientations or signed total gap of the shared rows.
No eventual terminal cutoff enters the bound.

The proof shows that, at degree \(m\ge4\), endpoint-normalized
roots lie at log radii of order \(V/m\), uniformly even when
\(mV\) is unbounded. The coefficient estimates used when
\(mV\) is bounded are replaced by exponential attenuation and
a lower bound for the exact symmetric resolvent. The physical
integer drift then reduces the question to the
[balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md).

## Symmetrization and coefficient bounds

Cancel the private unit factorials against the binomial
coefficient and normalize the constant and leading coefficients
to one. Write \(P_m(z)=\sum_{k=0}^m c_k z^k\), and let
\(f(x)\) be the continuous endpoint-normalized log coefficient
on \(0\le x\le m\). In this proof \(0<V\le1/80\).
The trigamma series and its upper and lower integral bounds give
\[
       \frac{2V}{3m}\le f''(x)\le\frac Vm,\qquad
       -V/2\le f'(x)\le V/2,\qquad f'(0)\le-V/3 .
 \tag{1}
\]
For the lower curvature, on each row use
\(d^2/(m(\ell+d)+1)\), where \(\ell\) is its smaller exponent
and \(d\) its absolute gap. Both \(d/\ell\) and \(1/\ell\)
are at most \(V\) for an unequal ordinary row. The derivative
bounds follow by integrating between zero endpoint values.
The lower bound in (1) gives
\(f(x)\le-V\min(x,m-x)/6\). All the same bounds hold for
\(g(x)=f(m-x)\).

Symmetrize the coefficients to obtain
\(C_m(z)=\sum_{k=0}^m(c_k+c_{m-k})z^k/2\).
Both arrays are geometric moments of positive variables after
private cancellation. More explicitly, if \(T>0\) has
\(\mathbb E T^x=e^{f(x)}\), then tilting its law by \(T^m\)
and applying inversion gives log moments \(g(x)\); the tilt
is a probability law since \(c_m=1\). Let \(\mu\) be the
equal mixture of these two laws. Its moment array is the
reference \(C_m\), and it is invariant under normalized
\(T^m\) tilt followed by inversion. The
[geometric argument boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)
give one simple root per argument box. Reciprocal conjugation
fixes the argument, so the reference roots lie on the unit
circle. Their circular angular spacing, and their distance in
argument from the positive ray, are at least \(\pi/m\).
[Spacing at a reciprocal reference](MULTIVARIATE_PRIVATE_UNIT_WEAK_RADIAL_INTEGER_DRIFT.md#a-quadratic-bound-on-the-radial-displacement)

Put \(t(x)=\min(x,m-x)\) and \(w(x)=(e^{f(x)}-e^{g(x)})/2\).
The [trigamma skew estimate](MULTIVARIATE_PRIVATE_UNIT_WEAK_RADIAL_INTEGER_DRIFT.md#normalization-and-a-bound-for-the-coefficient-skew)
gives \(|f'-g'|\le V^2/2\); its proof uses only the shared
curvature and has no bound on \(mV\). Thus
\[
 |w(x)|\le e^{-Vt(x)/6}\frac{V^2t(x)}4,\qquad
 |w'(x)|\le e^{-Vt(x)/6}
                   \left(\frac{V^2}{4}+\frac{V^3t(x)}8\right).
 \tag{2}
\]
The second estimate follows by differentiating the difference,
using (1) and the bound for \(w\).
Since \(\sum_{k\ge1} k e^{-sk}\le s^{-2}\) for \(s>0\),
summing the first bound gives \(\sum_{k=0}^m |w(k)|\le18\).
Integrating the second bound on the two halves gives
\(\sum_{k=1}^m|w(k)-w(k-1)|\le12V\).

## A reference slope at every root

Assume \(mV\ge1\) and take a reference root
\(\zeta=e^{i\theta}\), \(0<\theta\le\pi\).
Set \(Y=|\log T|\) under the symmetric moment law \(\mu\).
Differentiation of its two components' log Mellin transforms
at zero gives means at most \(-V/3\) with absolute values at
most \(V/2\), and variances at most \(V/m\), by (1).
Consequently
\[
             \mathbb E_\mu Y\ge V/3,\qquad
             \mathbb E_\mu Y^2\le 5V^2/4 .
\]
Also \(f,g\) extend to \(x\in[-1,1]\), where the upper
trigamma bound gives \(f'',g''\le(5/4)V/m\).
Indeed the smallest factorial argument before adding one is
at least \(m\ell-d\) on this extension, and \(d/\ell\le V\).
Using \(V/m\le V^2\) and the slope bound at zero, Taylor's
identity gives
\[
 \mathbb E_\mu\cosh Y
       \le e^{5V^2/8}\cosh(V/2)
       \le e^{3V^2/4}\le 1+3V^2/2 .
\]

Let
\[
 E(\theta)=\mathbb E_\mu(1-e^{i\theta}T)^{-1},\qquad
 D_\theta(y)=\cosh y-\cos\theta .
\]
The exact finite geometric sum and inversion identity give
\[
 e^{-im\theta/2}C_m(e^{i\theta})
           =2\operatorname{Re}(e^{-im\theta/2}E(\theta)).
 \tag{3}
\]
Conditional on \(Y=y>0\), the two signs of \(\log T\) have
mass ratio \(e^{-my}\) for the positive sign relative to the
negative sign. This follows from the invariance of \(\mu\).
Pairing their real and imaginary parts yields
\[
 \operatorname{Re}E=\frac12+
    \mathbb E_\mu\frac{\tanh(mY/2)\sinh Y}{2D_\theta(Y)},\qquad
 \operatorname{Im}E=\mathbb E_\mu\frac{\sin\theta}{2D_\theta(Y)} .
 \tag{4}
\]
In particular \(\mathbb E_\mu D_\theta(Y)\le3V^2/2+\theta^2/2\).
The inequality \(\tanh(mY/2)\ge Y/(2V+Y)\) follows from
\(e^u\ge1+u\) and \(mV\ge1\). Hölder's inequality gives
\[
 \mathbb E_\mu Y
 \le \left(\mathbb E_\mu\frac{Y}{\sqrt{2V+Y}}\right)^{2/3}
                     \bigl(\mathbb E_\mu Y(2V+Y)\bigr)^{1/3}.
\]
Now use \(\sinh Y\ge Y\), this bound, and Cauchy–Schwarz
in the real part of (4), obtaining
\[
 \operatorname{Re}E
   \ge \frac{(\mathbb E_\mu Y)^3}
              {2\,\mathbb E_\mu Y(2V+Y)\,\mathbb E_\mu D_\theta(Y)}
   \ge \frac{2V}{207\,\mathbb E_\mu D_\theta(Y)} .
\]
For the last step \(r=\mathbb E_\mu Y/V\ge1/3\) and
\(r^3/(2r+5/4)\) is increasing in \(r>0\).
If \(\theta\le V\), this gives \(\operatorname{Re}E\ge1/(207V)\).
If \(V\le\theta\le\pi/2\), Jensen's inequality in the imaginary
part gives \(\operatorname{Im}E\ge1/(8\theta)\), using
\(\sin\theta\ge\theta/2\). At \(\theta\ge\pi/2\) use
\(\operatorname{Re}E\ge1/2\). Thus
\[
                        |E(\theta)|\ge\frac1{207(V+\theta)}.
 \tag{5}
\]

Write \(\phi=\arg E\). Equations (4) express this argument
as the argument of \(R-\cos\theta+i\sin\theta\), where \(R\)
is the expectation of
\(\cosh Y+\tanh(mY/2)\sinh Y\)
under the law reweighted by \(D_\theta(Y)^{-1}\).
This integrand is increasing in \(Y\). Differentiating that
expectation gives its covariance with
\(-\sin\theta/D_\theta(Y)\), another increasing function of
\(Y\) for \(0<\theta<\pi\). Hence \(R'\ge0\).
With \(u=(R-\cos\theta)/\sin\theta\ge0\), differentiation
of the argument shows
\[
 \phi'=\frac{R\cos\theta-1-R'\sin\theta}{(R-\cos\theta)^2+\sin^2\theta}
     \le\frac{u\cot\theta-1}{1+u^2}
     \le \tfrac12\max(0,\cot\theta).
\]
At the reference roots \(\theta\ge\pi/m\), so \(\phi'\le m/6\).
At \(\pi\) the same upper bound follows by continuity in
the first expression, whose denominator stays positive.
Differentiate (3) at a root and use (5) to obtain
\[
                        |C_m'(\zeta)|\ge\frac{m}{400(V+\theta)} .
 \tag{6}
\]

## Uniform radial disks and the physical integer drift

For \(|z-\zeta|\le1/m\) one has \(|z|^m<3\),
\(|z|^{m+1}<4\) and \(|1-z|\ge\theta/4\).
The last inequality uses
\(|1-\zeta|\ge2\theta/\pi\) and \(\theta\ge\pi/m\).
Absolute summation of the attenuated coefficients bounds
\(|C_m(z)|\) by \(37/V\). The symmetric coefficients are
log-convex, as positive moments, and are one at the endpoints.
Their total variation is at most two. Thus multiplying \(C_m\)
by \(1-z\) bounds the modulus by \(1+4+6=11\).
Combining the two estimates gives
\[
     |C_m(z)|\le\frac{100}{V+\theta},\qquad
     |P_m(z)-C_m(z)|\le\frac{200V}{V+\theta}.
\]
For the latter use respectively the sum and the total
variation from (2), obtaining upper bounds \(54\) and
\(144V/\theta\). The minimum of \(a/V\) and \(b/\theta\)
is at most \((a+b)/(V+\theta)\), which gives both combined
estimates in the stated form.

Cauchy's estimate for the second derivative on the half-size
disk bounds it by \(800m^2/(V+\theta)\).
On the circle of radius \(\eta=160000V/m\) about
\(\zeta\), use (6) and Taylor's integral remainder. If
\(V\le10^{-11}\),
\[
 |C_m(z)|\ge\frac{m\eta/400-400m^2\eta^2}{V+\theta}
           \ge\frac{297V}{V+\theta}
           >|P_m(z)-C_m(z)|.
\]
The circle is inside the half-size disk, and the reference
spacing makes all these circles disjoint. Conjugation covers
the lower half-plane. Rouché's theorem gives one actual root
in each reference disk. It follows that every endpoint-normalized
root in this regime satisfies
\[
                     |\log|z||\le320000V/m .
 \tag{7}
\]
For \(mV\le1\) the
[weak radial bound](MULTIVARIATE_PRIVATE_UNIT_WEAK_RADIAL_INTEGER_DRIFT.md#a-quadratic-bound-on-the-radial-displacement)
gives \(2V^2\le2V/m\), so (7) holds there as well.

Let \(g_m\) be the positive endpoint scale of the physical
degree-\(m\) moments. The
[exact center comparison](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
states
\[
          \left|\log(g_N/g_n)-\Delta\log(N/n)\right|
                         \le V/(nN).
\]
At a common physical root with \(n\ge4\), (7) would imply
\[
 |\Delta|\log(N/n)\le
             320000V(1/n+1/N)+V/(nN).
\]
Since \(\log(N/n)>1/N\) and \(N/n\le5/4\),
this gives \(|\Delta|<720001V<1\).
The signed degree difference is an integer, so it must be
zero. The balanced Euclidean theorem then excludes the common
root. Terminals at most four follow from the
[low-order separation theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md).
At \(V=0\) there is no shared gap and the moments are geometric.
All coefficients are rational in the ordinary case, so the
complex coprimality transfers to every characteristic-zero
field.
