# Angle boxes and all-order interlacing with one mirrored private unit pair

Take ordinary independent variables \(X,Y\) and any
number of independent shared mirrored rows with
smaller exponents \(\ell_i\ge1\) and gaps \(d_i\ge1\).
Write the two monomials as \(X U\) and \(Y V\), where
\(U,V\) are the shared mirrored products.
Additional equal rows contribute a nonzero scalar.
Put
\[
                      S_{\rm sh}=2\sum_i\frac{d_i^2}{\ell_i+d_i}.
\]

**Theorem.** If \(S_{\rm sh}\ge21/4\), the consecutive
ordinary moment polynomials for \(XU+zYV\) strictly
interlace on the circle in every degree. Thus SFC(2)
holds over every characteristic-zero field for the
indicated two monomials and arbitrary coefficients.
The total absolute curvature condition is \(S\ge29/4\)
in this family.

The private unit pair removes the binomial coefficient
from the normalized polynomial. It therefore gives
a single resolvent representation in the angle
variable. The resulting angle boxes and one-sided
alternation apply more generally to exchangeable
positive \(U,V\) with finite moments, independently
of the factorial form of the shared core.

## The exact resolvent and angle boxes

For the independent ordinary private variables,
\(\binom mk k!(m-k)!=m!\), so after positive
normalization the polynomial is
\[
 P_m(z)=\sum_{k=0}^m a_{m,k} z^k,\qquad
         a_{m,k}=\frac{\mathbb E[U^{m-k}V^k]}{\mathbb E[U^m]} .
\]
The expectation is with independent unit exponentials
for the shared variables. The following representation
only uses exchangeability of \(U,V\).
Under probability tilted by \(U^m/\mathbb E[U^m]\)
write \(\mathbb E_m\) for expectation and \(T=V/U\).
Then \(a_{m,k}=\mathbb E_m T^k=a_{m,m-k}\),
and the exchange identity is
\[
                   \mathbb E_m[T^m f(T)]
                             =\mathbb E_m[f(1/T)].
\]
For \(z=e^{i\theta}\) with \(0<\theta<\pi\), the
finite geometric sum and this identity give
\[
 P_m(z)=K_m(z)+z^m\overline{K_m(z)},\qquad
 K_m(z)=\mathbb E_m(1-zT)^{-1}.
\]
Indeed the second term of the geometric sum is
\(-z^{m+1}\mathbb E_m[T^{m+1}/(1-zT)]\),
which becomes \(z^m\mathbb E_m[(1-z^{-1}T)^{-1}]\).

Set \(D_\theta(T)=1-2T\cos\theta+T^2\),
\(A_m=\mathbb E_m[1/D_\theta(T)]\),
\(B_m=\mathbb E_m[T/D_\theta(T)]\), and
\(\mu_m=B_m/A_m\). Then
\(K_m=A_m-e^{-i\theta}B_m\).
Applying the exchange identity to \(A_m\) and
\(B_m\) gives
\[
 A_m-B_m=\tfrac12\mathbb E_m
                 \left[\frac{(1-T)(1-T^{m+1})}{D_\theta(T)}\right]>0
\]
whenever \(T\) is nonconstant. Thus \(0<\mu_m<1\)
and the phase \(\phi_m=\arg K_m\) satisfies
\[
       0<\phi_m(\theta)
         =\arg(1-\mu_m e^{-i\theta})
         <(\pi-\theta)/2 .
\]
The upper endpoint is the phase for \(\mu_m=1\).
If \(T=1\) almost surely, the moments give the
geometric-sum polynomials, which already strictly
interlace.

The centered polynomial on the circle is
\(2|K_m|\cos(m\theta/2-\phi_m)\).
For each \(1\le j\le\lfloor m/2\rfloor\)
its signs are opposite at
\((2j-1)\pi/m\) and \(2j\pi/(m+1)\),
by the strict phase bounds. These intervals
are disjoint. Conjugation, and the structural
root at \(-1\) in odd degree, account for the
full degree, with one simple root per indicated
interval. In increasing positive-half angle
order the boxes are therefore
\[
                  \frac{(2j-1)\pi}{m}
                        <\theta_{m,j}<\frac{2j\pi}{m+1}.      \tag{1}
\]
In particular
\(\theta_{n,j}<\theta_{n+1,j+1}\) whenever both
roots are in the open positive half.
The other inequality
\(\theta_{n+1,j}<\theta_{n,j}\)
also follows from (1) whenever
\(4j\ge n+2\). It remains to control the
early upper roots, whose phase depends on the
matched endpoint laws.

## Bounding the early upper phases

The
[ordinary discrete curvature estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
applied to the shared core gives, with \(S_{\rm sh}\)
as above,
\[
             a_{m,k}\le
                  \exp(-S_{\rm sh} k(m-k)/(2m)).
                                                               \tag{2}
\]
For \(m\ge5\) and \(\pi/m\le\theta\le\pi/2\), we
prove \(\mu_m(\theta)<1/3\) when \(S_{\rm sh}\ge21/4\).
For this range \(A_m\ge 1-a_{m,2}\) since
\(1/D_\theta(T)\ge1/(1+T^2)\ge1-T^2\).
Writing \(c=\cos\theta\), define \(C_\theta=1/[2(1-c)]\).
The two squares
\[
      D_\theta(T)-2(1-c)T=(T-1)^2,\qquad
      D_\theta(T)-2(1-c)(2c-T)=(T-(2c-1))^2
\]
give the bounds \(T/D_\theta(T)\le C_\theta\) and
\(T/D_\theta(T)\le T+C_\theta T^2\), respectively.
In this angular range \(C_\theta\le m^2/9\),
because
\(\sin(\pi/(2m))\ge3/(2m)\) for \(m\ge5\).
The last inequality follows from
\(\sin u\ge u-u^3/6\), \(\pi>31/10\) and
\(\pi^2<10\).

For \(5\le m\le12\) use the second-moment bound
\[
                   B_m\le a_{m,1}+(m^2/9)a_{m,2}.
\]
When \(5\le m\le8\), (2) gives \(a_{m,1}<1/8\),
\(a_{m,2}<1/20\), and \((m^2/9)a_{m,2}<1/7\).
For the last bound the function
\(m^2 e^{-(21/4)(1-2/m)}\) has no interior
maximum in the interval, as is immediate
from its logarithmic derivative.
The endpoint comparisons are
\((25/9)e^{-63/20}<1/7\) and
\((64/9)e^{-63/16}<1/7\).
For \(9\le m\le12\), the bounds are
\(a_{m,1}<1/10\), \(a_{m,2}<1/50\), and
\((m^2/9)a_{m,2}<21/100\). Here the last
bound uses \(16e^{-35/8}<21/100\).
In both intervals the resulting upper bound
for \(B_m/A_m\) is below \(1/3\).

For \(m\ge13\), on \(T\le2/3\) use
\(T/D_\theta(T)\le T/(1-T)^2\le T+12 T^2\).
For the other part use \(C_\theta\) times the
probability bound
\(\mathbb P_m(T>2/3)\le(3/2)^k a_{m,k}\)
with \(k=\lfloor m/2\rfloor\). Thus
\[
                  B_m\le a_{m,1}+12 a_{m,2}
                             +(m^2/9)(3/2)^k a_{m,k}.
\]
Now \(a_{m,1}\le e^{-63/26}<9/100\),
\(a_{m,2}\le e^{-231/52}<3/250\), and
\((3/2)e^{-21/16}<81/200\).
The tail term is bounded by
\((2k+1)^2(81/200)^k/9\), a decreasing
function on integers \(k\ge6\).
Its value at six is less than \(17/200\).
Therefore
\[
                 \frac{B_m}{A_m}
                   <\frac{9/100+12(3/250)+17/200}{1-3/250}
                   <\frac13 .
\]
All exponential comparisons in these bounds
follow from the positive Taylor sum for \(e^x\)
through degree twelve. The arguments bound
whole degree ranges and do not require root
verification at a finite list of degrees.

## Consecutive order comparison

For terminal \(m=n+1\ge5\), every upper
positive-half root with angle \(\theta\le\pi/2\)
satisfies \(\theta\ge\pi/m\) by (1).
At such a root \(\mu_m<1/3\) gives
\[
      \phi_m=\arctan\frac{\mu_m\sin\theta}{1-\mu_m\cos\theta}
                                               <\theta/2 .
\]
At the root in box \(j\), the phase is \((j-1/2)\pi\).
The phase root equation then puts the upper
\(j\)-th angle strictly below
\((2j-1)\pi/(m-1)\), which is in turn below
the lower \(j\)-th angle whenever it exists
in the open half-circle.
For an upper angle above \(\pi/2\), its
upper box in (1) implies \(4j>m+1\), so
the same order follows directly from the
boxes. Together with the other side of
alternation in (1) and the odd structural
root this gives strict interlacing.

For \(n=1,2\), compare the elementary
polynomials \(1+z\), \(1+a z+z^2\), and
\(1+bz+bz^2+z^3\) with \(0<a,b\le1\).
The positive-half cosine at degree two
is \(-a/2<0\), while at degree three it is
\((1-b)/2\ge0\), and the odd root is
structural. At \(n=3\) write the lower
polynomial \(1+a z+a z^2+z^3\) and the
upper \(1+b z+c z^2+b z^3+z^4\).
One has \(c\le b\le a\):
\(c\le b\) follows by log convexity and
symmetry, and \(b\le a\) follows from
the product of row factors for \(a_{m,1}\),
\[
                     a_{m,1}=
          \prod_i\prod_{h=1}^{d_i}
             \frac{m\ell_i+h}{m(\ell_i+d_i)-d_i+h},
\]
each decreasing in \(m\).
The upper cosine quadratic
\(4x^2+2b x+c-2\) evaluated at the
lower positive-half cosine
\(x=(1-a)/2\) is at most
\((1-a)^2+a(1-a)+a-2=-1\).
Thus the lower root lies strictly between
the two upper roots. This proves all
remaining orders by coefficient identities.

The resulting rational consecutive polynomials
are coprime, and their Bézout identity transfers
to every characteristic-zero field.
