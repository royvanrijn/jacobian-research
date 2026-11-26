# All degree physical separation through shared curvature twenty three sixtieths

Take an ordinary two-monomial profile with exactly one private
unit row on each side and arbitrary shared rows.
For the shared signed gaps \(d_i=b_i-a_i\) and lower exponents
\(\ell_i=\min(a_i,b_i)\), put \(V=\sum_i d_i^2/\ell_i\) and
\(\Delta=\sum_i d_i\).

**Theorem.** If \(V\le23/60\), all ordinary consecutive physical
factorial moments are coprime. This holds for every signed integer
degree difference. The proof has no terminal cutoff and uses no
finite terminal verifier. It sharpens the
[three-tenths criterion](MULTIVARIATE_PRIVATE_UNIT_WEAK_THREE_TENTHS.md)
by using an evenized second exponential moment of the exact
adjacent Beta noise.

## Evenized exponential moment

The [linked endpoint law](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#exact-signed-log-noise)
in either orientation has independent additive log noise \(Z\).
Write \(M(t)=\mathbb E e^{tZ}\), \(\Lambda=\log M\),
\(\mu=\mathbb E Z\), and \(L_0=\Lambda''(0)\).
The same proof of the positive signed mean gives
\[
       0\le\mu\le V^2/(4nN),\qquad L_0\le V/(nN).
\]
In the mean bound use
\(\sum_i d_i^2/\ell_i^2\le V^2\) since
each \(|d_i|\ge1\).
For the variance there is an exact finite expression in either
orientation, because the second Beta shape \(a_i\) is integral:
\[
                \Lambda''(t)=
                \sum_i\sum_{h=1}^{a_i}
                     \frac{d_i^2}{(n a_i+h+d_i t)^2}.
\]
The bound on \(L_0\) follows by integral comparison.
All exponential moments differentiated here are finite on
\(|t|<n/V\), since \(|d_i|/a_i\le V\).

At an angle \(\theta\) of a candidate common root, write
\(b=n\sin\theta\), \(\alpha=1/\sin\theta\), and
\(\eta=V/b<1\). The displayed finite expression gives
\[
 \frac{\Lambda''(\alpha)+\Lambda''(-\alpha)}2
      \le L_0 \frac{1+\eta^2}{(1-\eta^2)^2},\qquad
 |\Lambda'(\pm\alpha)-\mu|\le \frac{\alpha L_0}{1-\eta}.
\]
Indeed the even average of each reciprocal square is bounded by
its zero-shift value times \((1+\eta^2)/(1-\eta^2)^2\). For
the derivative integrate the bound
\(L_0(1-V|t|/n)^{-2}\).
Integrating the bound a second time gives
\[
     \Lambda(\pm\alpha)\le K=\alpha\mu+
                             \frac{\alpha^2 L_0}{2(1-\eta)}.
\]
Here \(\int_0^1(1-u)/(1-\eta u)^2\,du\le1/[2(1-\eta)]\),
by the power series with coefficients \(1/(j+2)\).

Set \(C=(M''(\alpha)+M''(-\alpha))/2
     =\mathbb E[Z^2\cosh(\alpha Z)]\).
The estimates above yield \(C\le c\,\mathbb E Z^2\) with
\[
 c=e^K\left[\frac{1+\eta^2}{(1-\eta^2)^2}
         +\frac{2\mu\alpha}{1-\eta}
         +\frac{\alpha^2 L_0}{(1-\eta)^2}\right].
                                                               \tag{1}
\]
Use \(M''=M(\Lambda''+(\Lambda')^2)\) and
\(\mathbb E Z^2=L_0+\mu^2\).
By Cauchy–Schwarz and
\(\sinh^2 u/\cosh u\le2(\cosh u-1)\),
\[
  \mathbb E[Z^2 e^{\alpha|Z|}]
       \le C+\sqrt{2C(C-\mathbb E Z^2)}
       \le (c+\sqrt{2c(c-1)})\,\mathbb E Z^2 .                  \tag{2}
\]
This improves the use of the sum of the two tilted moments
in the relative resolvent error.

## Integer drift with the sharper packet estimate

The two packet differences at the same log radius are bounded
as in the [linked resolvent comparison](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#resolvent-variation-at-any-horizontal-position),
using relative error
\[
         \varepsilon=\mu/\sin\theta+
                       \mathbb E[Z^2 e^{\alpha|Z|}]/\sin^2\theta.
\]
Both orientations satisfy the same upper bound.
The normalized log radii \(x_n,x_N\) at two roots of the
same angle satisfy the
[linked radius inequality](MULTIVARIATE_PRIVATE_UNIT_WEAK_THREE_TENTHS.md#modulus-estimate-with-the-weak-row-cap)
\[
 N|x_N-x_n|\le
 \frac{b(1+1/n)}{b(1+1/n)-2}
       \left(\frac{V(n-1)}{2n}+2[-\log(1-\varepsilon)]\right).
                                                               \tag{3}
\]
The first term inside uses the strengthened Eneström bound.

At a physical common root, \(|x_N-x_n|=|\log(G_N/G_n)|\).
Here the endpoint drift estimate sharpens to
\[
               |\log(G_N/G_n)-\Delta\log(N/n)|\le V/(2nN).
                                                               \tag{4}
\]
Indeed \(1-u\psi'(1+u)\le1/(2u)\) follows from
\(1-h(v)\le v/2\) in the trigamma integral
\(\psi'(1+u)=\int_0^\infty e^{-uv}h(v)\,dv\).
Integrate \(x\psi'(1+mx)\) over the two endpoint exponents
of each row. The error relative to its signed gap divided
by \(m\) is at most \(\log((\ell_i+|d_i|)/\ell_i)/(2m^2)
\le d_i^2/(2m^2\ell_i)\). The private units cancel.
Now integrate \(m\) from \(n\) to \(N\).
Since \(\log(N/n)>1/N\), (3) plus \(V/(2n)\)
strictly bounds any nonzero \(|\Delta|\).

Put \(v=23/60\) for the numerical caps below.
For \(n\ge8\), \(b\ge61/20\) by
\(n\sin(\pi/n)\ge\pi-\pi^3/(6n^2)\).
Using that bound, \(N\ge9\) and \(n/N\le1\) gives
\(c<8/7\) in (1), as checked below. Thus (2) is at most
\((7/4)V/(nN)\). With \(x=1/n\) bound the error by
\[
 E(x)=
 \frac{v}{1+x}\left[\frac{7/4}{(61/20)^2}+
                                 \frac{v x}{4(61/20)}\right].
\]
As in the
[convex endpoint bound](MULTIVARIATE_PRIVATE_UNIT_WEAK_THREE_TENTHS.md#uniform-bounds-for-the-integer-drift),
the resulting drift majorant is convex in \(x\in[0,1/8]\):
\[
 \frac{(61/20)(1+x)}{(61/20)(1+x)-2}
        \left(\frac v2+2[-\log(1-E(x))]\right)+\frac{v x}{2}.
\]
The two factors of its product are positive, decreasing and convex.
Both endpoint values are below one, using
\(-\log(1-u)\le u+u^2/[2(1-u)]\).

For \(4\le n\le7\) use \(b\ge14/5\), \(N\ge5\) and
\(n/N\le7/8\).
This gives \(c<7/6\) and
\(\mathbb E[Z^2 e^{\alpha|Z|}]\le(11/6)V/(nN)\).
The contraction factor in (3) is at most \(8/3\).
Combining with \(V(n-1)/(2n)\le3v/7\),
the drift majorant is again below one.

For \(n=2,3\), the strengthened Eneström bound and (4)
give \(|\Delta|\log(N/n)\le (1-1/N)V\).
At \(n=2\) use \(\log(3/2)>2/5\).
At \(n=3\) use \(\log(4/3)>2/7+2/(3\cdot7^3)\),
by the positive series for \(2\operatorname{artanh}(1/7)\).
Both show \(|\Delta|<1\).
Terminal two follows from the real first root and strict
positivity of the second moment on the real axis.
Thus all possible common roots would require \(\Delta=0\),
and the [balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
excludes them.

## Exact scalar checks

~~~python
from fractions import Fraction as F

v = F(23, 60)

def cap(b, N0, p):
    eta = v/b
    K = v*v/(4*b*N0) + v*p/(2*b*b*(1-eta))
    G = ((1+eta*eta)/(1-eta*eta)**2
              + v*v/(2*b*N0*(1-eta))
              + v*p/(b*b*(1-eta)**2))
    return G/(1-K)

b = F(61, 20)
assert F(157, 50)-F(22, 7)**3/(6*8**2) > b
assert cap(b, 9, F(1)) < F(8, 7)
assert F(12, 7)*(1+v/(16*8*9)) < F(7, 4)
assert F(7,4)/b**2 > v/(4*b)
for x in [F(0), F(1,8)]:
    e = v*(F(7,4)/b**2+v*x/(4*b))/(1+x)
    assert b*(1+x)/(b*(1+x)-2)*(v/2+2*e+e*e/(1-e))+v*x/2 < 1

b = F(14, 5)
assert cap(b, 5, F(7,8)) < F(7,6)
assert 2*F(7,6)*(F(7,6)-1) < F(19,30)**2
assert F(9,5)*(1+v/(16*4*5)) < F(11,6)
e = v*v/(4*b*5)+F(11,6)*v*F(7,8)/b**2
assert F(8,3)*(3*v/7+2*e+e*e/(1-e))+v/8 < 1

assert F(3,4)*v < F(2,7)+F(2,3*7**3)
assert F(2,3)*v < F(2,5)
~~~
