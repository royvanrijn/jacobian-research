# Signed imbalance two beyond order 122 for private unit pairs

Consider an ordinary factorial binomial with one private unit
row on each side and shared rows of arbitrary orientations.
If the absolute signed degree difference is two, the physical
moments of adjacent degrees \(n,n+1\) are coprime for
\(n\ge123\). The proof uses explicit radius and coefficient
bounds, with no eventual terminal cutoff.

For \(n\ge123\), the
[higher-label theorem](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md)
excludes all labels above one. Reverse the columns if
necessary so that the signed degree difference is \(\Delta=2\).
The
[centered first-angle comparison](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_THREE.md#reduction-to-the-middle-coefficient-quotient)
can use any integer \(j\) with \(j\theta<\pi\) and
\((N-j)\theta<\pi\), where \(N=n+1\) and \(\theta\) is the
proposed first upper root angle. Its coefficient-ratio test
is \(v_n\ge v_j\). The linked factorial loss is bounded by
the cost in (2) of the
[signed loss comparison](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_THREE.md#the-signed-loss-over-the-upper-half)
with \(j/n\) as the integration lower limit.

Write \(S=\sum_i d_i^2/(\ell_i+d_i)\) for the shared
absolute curvature, with \(d_i=|b_i-a_i|\) and
\(\ell_i=\min(a_i,b_i)\) on unequal shared rows.
The [curvature fifteen-eighths theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
already covers \(S\ge15/8\).
If \(S<4/3\), every row has \(d_i/\ell_i\le1\),
by the integer packing in the imbalance-three proof.
Use \(j=\lceil n/2\rceil\) and its cost at most
\((26/19)S<2\). This rules out the first box by (3)
of that proof.
It remains to treat
\[
                            4/3\le S<15/8 .
\]

## A sharper first angle at this curvature

Take \(m\ge123\). Normalize both endpoint coefficients
of the degree-\(m\) moment to one. The
[absolute normalized radius bound](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md)
gives \(|\log r|\le11/m\) at any normalized root
\(z=r e^{i\theta}\).
In each orientation, the normalized coefficient \(c_k\)
has the bound
\[
      c_k\le\exp\left[-\frac{\kappa}{2}k(1-k/m)\right],
       \qquad \kappa=\sum_i\frac{d_i^2}{\ell_i+d_i+1/m}
                        \ge\frac{33}{25}.
 \tag{1}
\]
In fact the row log interpolation has second derivative
at least \(d_i^2/[m(\ell_i+d_i)+1]\) by the trigamma
integral bound. Comparing its zero-endpoint chord with
that of this constant curvature gives (1).
The coefficient of \(S\) in this bound is at least
\(2m/(2m+1)>99/100\), since \(\ell_i+d_i\ge2\).

Split the polynomial into the left half and the reverse
half times \(z^m\), giving half weight to the middle
term when \(m\) is even. Use the angles \(+\theta\)
for both half arguments, conjugating the evaluated reverse
half. Denote their arguments by \(\phi_L,\phi_R\).
The half coefficients at the root radius, before their
angle factors, satisfy the following bounds in either
orientation:
\[
    c_k r^{\pm k}\le
        \begin{cases}
        (13/20)^k,& k\le m/5,\\
        (4/5)^k,&m/5<k\le m/2.
        \end{cases}                                           \tag{2}
\]
Indeed the first exponent bound from (1) is
\((2/5)(33/25)-11/123\), and the second is
\((1/4)(33/25)-11/123\). Their exponentials exceed
\(20/13\) and \(5/4\), respectively, by the positive
Taylor series.

At a first root \(0<\theta<2\pi/(m+1)\), the early
terms \(k\le m/5\) have positive real parts.
The late terms \(k>m/5\) start at \(k\ge25\).
For each half, (2) therefore gives real part at least
\(1-5(4/5)^{25}\). Its imaginary part is nonnegative,
and is at most
\(\theta[260/49+145(4/5)^{25}]\), by
\(\sin(k\theta)\le k\theta\) and the geometric series
for the two bounds. Consequently
\[
               0\le\phi_L,\phi_R<6\theta .
\]
The root equation now gives \(m\theta=\pi+\phi_L+\phi_R\):
the phase \(m\theta-\phi_L-\phi_R\) is an odd
multiple of \(\pi\), positive because \(m>12\),
and less than \(2\pi\) by the first box.
In particular
\[
                              (m-12)\theta<\pi .
 \tag{3}
\]

## The cost at index twelve

For a proposed common first angle, apply (3) at the upper
order \(m=N\). This makes \(j=12\) allowable in the
centered imaginary projection, with \(12\theta<\pi\)
by the first box and \((N-12)\theta<\pi\) by (3).
We have \(j/n\le12/123<49/500\).
The increasing-row cost is at most its \(S_i\), and
for a decreasing row with \(r_i=d_i/\ell_i\) the upper
cost is the average of \(d_i r_i y/[1+r_i(1-y)]\)
over \([\eta,1]\), where \(\eta=49/500\).
This upper bound increases with the left integration limit,
so it suffices here.

For \(r_i\le1\), the chord bound for this convex integrand
gives the cost-to-\(S_i\) factor
\[
              \frac{1+r_i}{2}
                   \left(1+\frac{\eta}{1+(1-\eta)r_i}\right)
                        \le\frac{1000}{951}.
\]
This also bounds the increasing-row factor.
If there is no exceptional decreasing row with \(r_i>1\),
the total cost is below \((1000/951)(15/8)<2\).
The only possible exceptions are
\((\ell_i,d_i)=(1,2)\) and \((2,3)\), at most one,
by the earlier integer packing. Integrating this row
exactly gives respectively
\[
               \frac{1500}{451}\log(701/250)-2,\qquad
               \frac{2500}{451}\log(2353/1000)-3 .
\]
The remaining \(S\) costs less than \(13/24\) and \(3/40\),
respectively. With
\(\log(701/250)<10311/10000\) and
\(\log(2353/1000)<107/125\), the two combined bounds
are strictly below two. Formula (3) of the
[signed loss comparison](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_THREE.md#the-signed-loss-over-the-upper-half)
now gives \(v_n>v_{12}\), so the centered imaginary
projection excludes the first root. This proves the physical
coprimality assertion. The rational Bézout identity transfers
it to every characteristic-zero field.

## Exact scalar replay

```python
from fractions import Fraction as Q
from math import factorial

def log_upper(a, b):
    x = Q(a-b, a+b)
    return 2*(sum(x**(2*j+1)/(2*j+1) for j in range(8))
              + x**17/(17*(1-x*x)))

def exp_lower(t):
    return sum(t**j/factorial(j) for j in range(4))

assert Q(246,247)*Q(4,3) > Q(33,25)
assert exp_lower(Q(2,5)*Q(33,25)-Q(11,123)) > Q(20,13)
assert exp_lower(Q(1,4)*Q(33,25)-Q(11,123)) > Q(5,4)
tail = Q(4,5)**25
assert Q(260,49)+145*tail < 6*(1-5*tail)
eta = Q(49,500)
assert Q(12,123) < eta
soft = 1+eta/(2-eta)
assert soft == Q(1000,951)
assert soft*Q(15,8) < 2
assert log_upper(701,250) < Q(10311,10000)
assert log_upper(2353,1000) < Q(107,125)
assert Q(1500,451)*Q(10311,10000)-2+soft*Q(13,24) < 2
assert Q(2500,451)*Q(107,125)-3+soft*Q(3,40) < 2
```
