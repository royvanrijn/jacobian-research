# Limits of the first-sector density comparison

The first-sector result for
[comparable ordinary profiles](MULTIVARIATE_COMPARABLE_FIRST_SECTOR_HCM_RADIAL_COMPARISON.md)
uses an upper hyperbolically completely monotone density and an
independent contraction coupling to its lower density. Those
properties allow a common physical root beyond the first sector
for general positive bases, even with finite moments of every
order. A separate fixed positive law that is a sum of two
independent gamma variables also gives a common root.

## An upper gamma law with a completely monotone density ratio

Let \(X\) be gamma of shape \(10\) and rate one. Write
\(L_m(z)=\mathbb E(1+zX)^m\), so
\[
 \begin{split}
  L_4(z)&=1+40z+660z^2+5280z^3+17160z^4,\\
  L_5(z)&=1+50z+1100z^2+13200z^3+85800z^4+240240z^5 .
 \end{split}
\]
Let \(z_*\) be the unique zero of \(L_5\) in the disk
of radius \(10^{-6}\) centered at
\((-54556+57802i)/10^6\). The disk and the following bounds
are certified by the rational Taylor calculations below.
For \(t=(1/20,4/5,19/20,1)\), put \(v_j=L_4(t_j z_*)\).
Then the rectangular enclosures are
\[
\begin{array}{c|cc}
j & \Re v_j & \Im v_j\\ \hline
1 &(0.88,0.90)&(0.10,0.11)\\
2 &(-0.15,-0.13)&(0.09,0.10)\\
3 &(-0.11,-0.10)&(-0.04,-0.03)\\
4 &(-0.08,-0.07)&(-0.07,-0.06)
\end{array}
\]
The notation \([v,w]=\operatorname{Im}(\bar v w)\)
is the determinant of the two real coordinate vectors.
The three boxes give lower bounds on the cyclic determinants
of \(v_1,v_2,v_3\), namely
\(922/10000,\ 129/10000,\ 143/10000\).
The point \(-v_4/99\) is strictly inside that triangle.
Its coordinate magnitudes are below \(1/1000\), and translating
each vertex by the negative of this point changes each pair
determinant by less than \(3/1000\). The three translated
determinants remain positive, which gives the stated interior
condition.

Take a positive random variable \(B\) independent of \(X\).
Give it mass \(1/100\) at \(1\). Divide its remaining mass among
\(t_1,t_2,t_3\) using the barycentric coordinates of \(-v_4/99\)
in their triangle. This specifies positive algebraic probabilities
exactly, and yields \(\mathbb E L_4(Bz_*)=0\).

Let \(f_0,f_1\) be the densities of \(BX,X\), respectively.
The upper gamma density is hyperbolically completely monotone.
Their ratio is
\[
 h(t)=f_0(t)/f_1(t)
       =\mathbb E\big[B^{-10} e^{-(B^{-1}-1)t}\big].
\]
It is strictly decreasing and completely monotone. It is bounded
below by \(1/100\) and above by the finite value \(h(0)\).
To realize both laws as tilts of one fixed positive base,
draw \(T\) with density proportional to \(h(t)^4 f_0(t)\),
and set \(U=1/h(T)\), \(V=T/h(T)\). The tilt of this base
by \(U^4\) gives \(T\) density \(f_0\), and the tilt by \(U^5\)
gives density \(h^{-1}f_0=f_1\). All mixed moments are finite,
since \(U\) is bounded and the law of \(T\) has exponentially
decaying tails. It follows that
\[
               \mathbb E(U+z_*V)^4=\mathbb E(U+z_*V)^5=0 .
\]

Ordinary comparable rows have a further restriction on the
contraction. In the ordinary beta-gamma coupling, a row of
lower exponent \(a>0\) and gap \(d\) contributes
\(B_0^d\) for \(B_0\) beta of parameters \(na+1,a\).
For \(s\ge0\) its \(s\)-moment is
\[
  \mathbb E B_0^{ds}
    =\prod_{h=0}^{a-1}\frac{na+1+h}{na+1+h+ds}.
\]
Thus the negative logarithm of the contraction is a sum
of independent exponential variables with the indicated rates
\((na+1+h)/d\). Rows with lower exponent zero contribute
no factor. The four-point law in this example fails that
restriction. Here hyperbolic complete monotonicity is asserted
for the upper gamma density, and ordinary complete monotonicity
for the density ratio.

## A fixed sum of two gamma variables

Let \(s_1<s_2\) be the two roots of
\[
                     x^2-\frac{713}{47}x+\frac{1580}{141}=0.
\]
The roots are positive and bracket \(21/19\), by evaluating the
polynomial there and using its positive constant and positive
sum of roots. Define \(\rho_i>0\) by
\(\rho_1+\rho_2=19\) and \(s_1\rho_1+s_2\rho_2=21\).
The recurrence from the quadratic gives the first four values
of \(\sum_i\rho_i s_i^k\) for \(k=0,1,2,3\) as
\[
                           19,\quad21,\quad317/3,\quad4103/3 .
\]
Let \(T\) be the sum of independent gamma variables with
scales \(s_i\) and shapes \(\rho_i/s_i^2\). Its cumulants
at orders \(2,3,4,5\) are the four displayed values
multiplied respectively by \(1!,2!,3!,4!\), hence
\(19,42,634,32824\). Its mean is \(1884/79\), by using
the quadratic to calculate \(\sum_i\rho_i/s_i\).
Consequently, with \(Y=T-\mathbb ET\),
\[
     \mathbb E(x+Y)^5=(x+4)(x^2-2x+101)^2,\qquad
     5\mathbb E(x+Y)^4=\frac{d}{dx}\mathbb E(x+Y)^5 .
\]
The two polynomials share the nonreal roots \(1\pm10i\).
Taking \(U=1\), \(V=T\) gives common nonzero roots of the
two physical moment polynomials by
\(z=1/(x-\mathbb ET)\). This fixed law is a generalized
gamma convolution. Both examples have roots in the open left
half-plane, beyond the respective first sectors.

## Exact arithmetic replay

~~~python
from fractions import Fraction as Q
from math import comb, prod

def taylor(m, center):
    re, im = center
    powers = [(Q(1), Q(0))]
    for k in range(m):
        x, y = powers[-1]
        powers.append((x*re-y*im, x*im+y*re))
    coefficients = [comb(m,k)*prod(range(10,10+k)) for k in range(m+1)]
    return [
        (
            sum(comb(k,j)*coefficients[k]*powers[k-j][0]
                for k in range(j,m+1)),
            sum(comb(k,j)*coefficients[k]*powers[k-j][1]
                for k in range(j,m+1)),
        )
        for j in range(m+1)
    ]

def norm1(value):
    return sum(map(abs,value))

c = (Q(-54556,10**6), Q(57802,10**6))
r = Q(1,10**6)
a = taylor(5,c)
error = norm1(a[0]) + sum(norm1(a[k])*r**k for k in range(2,6))
assert error < Q(1,400000) < r*max(map(abs,a[1]))

boxes = [
    (Q(1,20),  88,  90, 10, 11),
    (Q(4,5),  -15, -13,  9, 10),
    (Q(19,20),-11, -10, -4, -3),
    (Q(1),     -8,  -7, -7, -6),
]
for t, lx, hx, ly, hy in boxes:
    a = taylor(4,(t*c[0],t*c[1]))
    error = sum(norm1(a[k])*(t*r)**k for k in range(1,5))
    assert Q(lx,100) < a[0][0]-error <= a[0][0]+error < Q(hx,100)
    assert Q(ly,100) < a[0][1]-error <= a[0][1]+error < Q(hy,100)
assert Q(8,9900) < Q(1,1000)
assert Q(202,100)*Q(1,1000) < Q(3,1000) < Q(129,10000)

p, q = Q(713,47), Q(1580,141)
assert p > 0 and q > 0 and Q(21,19)**2 - p*Q(21,19) + q < 0
m0, m1 = 19, 21
m2 = p*m1 - q*m0
m3 = p*m2 - q*m1
assert (m2,m3) == (Q(317,3),Q(4103,3))
assert (p*m0-m1)/q == Q(1884,79)
assert 6*m2+3*m0**2 == 1717
assert 24*m3+10*m0*(2*m1) == 40804
~~~
