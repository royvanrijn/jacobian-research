# All-degree weak curvature separation for two private units

Consider an ordinary two-monomial profile with exactly one private
unit row on each side. Shared rows may have either orientation.
Set
\[
                V=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\min(a_i,b_i)}.
\]
Let \(J\) be any integer at least \(10^8\) and at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Theorem.** If \(V\le17/J\), all ordinary consecutive physical
moments are coprime over every characteristic-zero field.
There is no restriction on the signed total degree difference
or the number and orientations of the shared rows.
The numerical value of the underlying absolute angular cutoff
is unevaluated.

The proof expands the weak-curvature
[integer drift argument](MULTIVARIATE_PRIVATE_UNIT_WEAK_RADIAL_INTEGER_DRIFT.md)
by bounding the resolvent derivative of the symmetrized coefficient
array. Once signed degree difference is zero, the exact
[positive-radius test](MULTIVARIATE_PRIVATE_UNIT_LINKED_POSITIVE_RADIUS_GATE.md)
uses curvature attenuation throughout the coefficient sum.

## Symmetrization and skew at bounded \(mV\)

After private cancellation and endpoint normalization write
\(P_m(z)=\sum_{k=0}^m c_k z^k\) with \(c_0=c_m=1\).
Let \(f(x)\) be the continuous endpoint-normalized log
coefficient on the index variable \(0\le x\le m\).
Throughout suppose \(2\le m<J\), \(0<V\le17/J\), so
\(mV\le17\) and \(V<1/4\).
The case \(V=0\) has consecutive geometric sums.

Writing \(\ell=\min(a,b)\) in each shared row, the trigamma
series gives \(0\le f''(x)\le V/m\); the derivative
is between \(-V/2\) and \(V/2\) by the zero endpoint values.
For the log skew \(h(x)=f(x)-f(m-x)\) one has
\[
                         |h''(x)|\le V^2/m .
\]
Indeed the arguments of each row's two trigammas differ
by at most \(m|d|\) for gap \(d=b-a\), and the second
polygamma series, in absolute value, is bounded by
\((m\ell)^{-2}\) on the interval between the arguments.
Then \(\sum d^2 |d|/\ell^2\le V^2\) since \(|d|/\ell\le V\).
The chord Green kernel and zero endpoints give
\(\lvert h\rvert\le mV^2/8\) and \(\lvert h'\rvert\le V^2/2\).

Put \(s_k=(c_k+c_{m-k})/2\) and \(P_m^s(z)=\sum s_k z^k\).
The antisymmetric coefficient difference extends to
\(w(x)=(e^{f(x)}-e^{f(m-x)})/2\). As \(f\le0\), its
derivative is bounded by
\[
                       |w'(x)|\le V^2/4+mV^3/32
                                            \le(25/32)V^2 .
\]
Differentiate \(2w=e^{f(m-x)}(e^h-1)\), using instead
\(2w'=e^{f(x)}h'+(e^{f(x)}-e^{f(m-x)})\,\partial_x f(m-x)\)
to bound the two terms directly. For the second term
use \(1\)-Lipschitzness of the exponential on the
nonpositive real axis and \(|f'|\le V/2\).
Thus the absolute coefficient mass of
\((1-z)(P_m-P_m^s)\) is at most \(25mV^2/32\).

The coefficients \(c_k\) are the moments of the
positive ratio \(T\) under the lower endpoint tilt,
with \(T\) scaled so \(\mathbb E[T^m]=1\).
The \(s_k\) are moments of the equal mixture of this
law and its reciprocal after tilting by \(T^m\).
This symmetrized law has the exact exchange identity.
It follows by the geometric resolvent and the
[both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
that the reference \(P_m^s\) has simple circular roots.
Their minimum angular spacing is at least \(\pi/m\),
also after adding the point \(1\), as follows immediately
from those boxes and reflection.

## A linear derivative bound on the circle

For either orientation the log moment of the normalized
\(T\) extends to \(x=-1\). For \(-1\le x\le1\), every
Gamma argument before adding one is at least
\(m\ell-|d|\ge(1-V/m)m\ell\). Thus
\(f''(x)\le V/[m(1-V/m)]\le2V/m\).
Taylor at zero gives
\[
 \mathbb E[\cosh(\log T)]\le e^{V/m}\cosh(V/2)
                          \le e^{V/m+V^2/8}.
\]
The last inequality follows from \(\tanh u\le u\).
The exponent is below \(1/2\) for \(V<1/4\) and
\(m\ge2\). Since \(mV\le17\), the right side is at most
\(1+107/m^2\), using \(e^x-1\le2x\) on \(0\le x\le1/2\).
Both orientations obey this bound, so it also holds
for the symmetrized law.

At \(z=e^{i\theta}\), \(0<\theta<\pi\), use its
exchange resolvent \(K=A-B e^{-i\theta}\), with
\(A=\mathbb E D^{-1}\), \(B=\mathbb E T/D\), and
\(D=1-2T\cos\theta+T^2\).
The centered reference is \(2|K|\cos\psi_m(\theta)\).
As in the
[uniform angular reserve](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_UNIFORM_ANGULAR_RESERVE.md#the-resolvent-and-its-phase-slopes),
the ratio \(R=A/B\) is at least one and increasing
in \(\theta\), by the even midpoint resolvent.
Consequently at a reference upper root
\(\psi_m'\ge m/4\): the positive contribution to
the derivative of \(\arg(R-e^{-i\theta})\) is at most
\(1/(2\sin\theta)\) for \(\theta\le\pi/2\) and zero
for \(\theta\ge\pi/2\). The first angle box gives
\(\sin\theta\ge2/m\) in the first range.
The same reasoning applies at the structural root
\(-1\) by continuity.

For \(\pi/m\le\theta\le\pi/2\), the identity
\(2B=\mathbb E[1/(\cosh(\log T)-\cos\theta)]\) and
Jensen's inequality give, with \(\delta=1-\cos\theta\),
\[
 |(1-z)K|
    \ \ge\ |1-z|B\sin\theta
    \ \ge\ \frac{\delta\cos(\theta/2)}{\delta+107/m^2}
    \ \ge\ \frac1{109}.
\]
Here \(\delta\ge2/m^2\) and \(\cos(\theta/2)\ge1/2\).
For \(\theta\ge\pi/2\), use \(|K|\ge A\ge1/4\)
and \(|1-z|\ge1\); the bound on \(A\) follows by
\(D\le(1+T)^2\) and \(\mathbb E T\le1\) as before.
It follows that at a root \(\zeta\) of \(P_m^s\),
\[
                         |[(1-z)P_m^s]'(\zeta)|\ge m/218 .     \tag{1}
\]

## Radial displacement and the signed degree lock

Draw disks of radius \(r=1200 V^2\) around the
reference roots of \(P_m^s\).
For \(G=(1-z)P_m^s\), the circular mesh of its roots
implies
\(\sum_{\xi\ne\zeta}|\xi-\zeta|^{-1}\le m(1+\log m)\)
at any fixed root \(\zeta\). Indeed each of the
two directions has chord distance at least \(2j/m\)
for cyclic index distance \(j\) in that direction
when it is the shorter circular direction.
Our bound on \(J\) gives
\[
     r m(1+\log m)
                 \le 1200\cdot17^2(1+\log J)/J < 1/2 .
\]
The disks are disjoint and exclude the root \(1\).
On a boundary the product factorization of \(G\)
and (1) give \(|G|\ge r m/436\), since the product
of \(1-r/|\xi-\zeta|\) is at least one minus the
sum of their losses. The perturbation
\((1-z)(P_m-P_m^s)\) is at most
\(3\cdot25mV^2/32<r m/436\) on the boundary,
using \((1+r)^m<3\).
Thus Rouché pairs all the ordinary roots with the
symmetric reference. In particular every normalized
root has absolute log radius at most \(2r=2400 V^2\).

Let a common physical root occur at \(2\le n<N<J\),
and write \(\Delta=\sum_{\rm shared}(b_i-a_i)\).
The
[exact endpoint scale drift bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
and the two normalized radius bounds give
\[
       |\Delta| <
             4800N V^2+V/n
             \le4800\cdot17^2/J+17/(2J)<1.
\]
This uses \(\log(N/n)>1/N\).
The integer \(\Delta\) is therefore zero.

## Attenuation in the positive-radius test

With \(\Delta=0\), the raw physical arrays are already
ordered in both orientations.
Use the
[scalar comparison](MULTIVARIATE_PRIVATE_UNIT_LINKED_POSITIVE_RADIUS_GATE.md#a-scalar-comparison-using-both-orientations)
with lower and upper leading coefficients \(a,b\),
positive test radius \(r_*=a/b\), and \(L=a r_*^n\),
reversing if necessary so \(L\ge1\).
Its coefficient difference \(d_k\) at index \(k\)
satisfies the sharper use of the attenuation bound
\[
 \frac{d_k r_*^k}{L}
         \le \frac{V k(k+1)}{2nN}\,
                     e^{-S k(1-k/n)/2},
 \qquad
 S=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\max(a_i,b_i)} .
\]
Indeed the physical derivative bound multiplies the
lower coefficient by \(V k(k+1)/(2nN)\).
Compared with its endpoint chord, that lower coefficient
has the displayed attenuation by the
[ordinary discrete factorial curvature bound](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate).
The endpoint chord at \(r_*\) is \(L^{k/n}\le L\).
Here \(S\ge V/(1+V)\ge(999/1000)V\), since
\(|d|/\ell\le V\) rowwise.

For \(n<J/4\) the earlier \(V\le6/(n+2)\) test applies
directly. For \(n\ge J/4\), put \(\lambda=nV\le17\),
\(q=999/1000\), and
\(F_\lambda(x)=x^2 e^{-q\lambda x(1-x)/2}\).
Then the sum of the normalized losses is at most
\[
 \frac{\lambda}{2}
      \left[\int_0^1 F_\lambda(x)\,dx + \frac{12}{n}\right].
                                                               \tag{2}
\]
For the Riemann sum use \(|F_\lambda'|\le2+17/2\)
and \(k(k+1)/(nN)\le(k/n)^2+1/n\).

The integral term multiplied by \(\lambda/2\) is less
than \(199/200\). To verify this, bound
\(e^{-q\lambda x(1-x)/2}\) by the degree-eight alternating
Taylor polynomial, whose remainder is nonpositive
for every nonnegative argument. Integrating yields the
polynomial
\[
 \frac{\lambda}{2}\sum_{i=0}^8
       \left(-\frac{q\lambda}{2}\right)^i
                        \frac{(i+2)!}{(2i+3)!}.
\]
Every Bernstein coefficient on \(\lambda\in[0,17]\)
is less than \(199/200\), as verified by the rational
calculation below. The Riemann error in (2) is at most
\(102/n<1/200\). The total normalized loss is thus below
one, giving the positive-radius criterion and excluding
the common physical root.

The pair at \(N=2\) follows by real positivity, and
\(N\ge J\) follows from the absolute angular cutoff.
The ordinary arrays have rational coefficients, hence
the resulting rational Bézout identities give field
descent as asserted.

~~~python
from fractions import Fraction as F
from math import factorial, comb

maximum = F(17)
q = F(999, 1000)
powers = [F(0)] + [
    maximum/2 * (-q*maximum/2)**i * F(factorial(i+2), factorial(2*i+3))
    for i in range(9)
]
for k in range(10):
    bernstein = sum(powers[i] * F(comb(k,i), comb(9,i)) for i in range(k+1))
    assert bernstein < F(199, 200)
assert F(1200 * 17**2 * 20, 10**8) < F(1, 2)
assert F(3*25, 32) < F(1200, 436)
assert F(4800*17**2, 10**8) + F(17, 2*10**8) < 1
assert F(102, 10**8//4) < F(1, 200)
~~~
