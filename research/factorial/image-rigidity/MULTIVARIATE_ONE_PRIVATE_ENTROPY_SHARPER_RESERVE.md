# A stronger entropy bound for one private side

Suppose an ordinary two-monomial profile has private variables
on exactly one side, of total private degree \(2\le P\le5\).
Write
\[
          S_{\rm sh}=\sum_{\rm shared}
                   \frac{(b_i-a_i)^2}{\max(a_i,b_i)} .
\]
If \(S_{\rm sh}+(35/18)P\ge12\), the consecutive factorial
moment polynomials have disjoint root arguments in every
degree. The resulting thresholds for \(P=2,3,4,5\) are
\(73/9,37/6,38/9,41/18\), respectively.
Physical SFC(2) follows over every characteristic-zero field.

This is a direct strengthening of the
[one-private entropy proof](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md).
The change retains more of the private factorial bound for
long halves and uses three short ranges.

## Long halves

If \(P+S_{\rm sh}\ge12\), use the
[absolute-curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md).
Otherwise the actual derivative budget of the entropy proof
satisfies \(V<168\), by its bound \(V\le S(S+2)\) for
\(S=P+S_{\rm sh}<12\). At real orders \(m\ge82\), its
\(\delta_m=\log(m)/(2(m-1))\) is strictly less than \(1/36\),
since it decreases and \(\log82<9/2\).
Consequently the private factor in its coefficient envelope
allows the chord parameter \(s_0=(35/18)P\).
This gives \(S'=s_0+S_{\rm sh}\ge12\) and directed budget
\(V_d'\le V+s_0-P<173\).
The actual Gamma derivatives use \(V<168\), \(P\le5\), and
the same relative slope bounds as in the entropy proof.
The quadratic derivative budget is covered by \(V_d'\),
since \(s_0\ge P\).

For the lighter orientation, the reserve bound of the entropy
proof gives \(c_{m,k}^{1/k}\) at most the shared envelope times
\(e^{3/2}e^{-P(1-w)}\), with \(w=k/m\).
This is at most \(e^{3/2}Q(w)\) for the new chord parameter,
since \(s_0\le2P\). Therefore the
[long-phase comparisons](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#matching-the-long-phase-hypotheses)
apply without changing their endpoint constants, derivative
bounds, or thresholds. They cover \(N=n+1\ge83\).

## Short halves

At integer \(13\le m\le82\), the coefficient cap in either
orientation is \(\binom mk e^{-\sigma k(1-k/m)}\), where
\[
 \sigma=P(1-\delta_m)+S_{\rm sh}/2
                    \ge 6+P(1/36-\delta_m).
\]
Use \(P\le5\) and upper bounds \(13/120,7/124,21/650\) on
\(\delta_m\) over \([13,31],[32,65],[66,82]\), respectively.
Each upper bound exceeds \(1/36\), so a common bound for the
group is \(\sigma\ge6+5(1/36-\delta_{\rm upper})\).
The same
[short-half coefficient and phase comparison](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#short-phase-intervals)
then applies. For \(k\le5\) take the two order endpoints
of the factorial envelope, and for \(k\ge6\) use the
entropy endpoints with geometric ratios \(1/4,11/50,203/1000\)
for the three groups. The replay tests the radii and within-group
and adjacent-group phase margins. This covers \(14\le N\le82\).
The [earlier-terminal arguments](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#matching-the-long-phase-hypotheses)
give the same conclusion for \(N\le13\). Positive rescaling
does not change root arguments.

## Exact scalar replay

~~~python
from fractions import Fraction as Q
from math import factorial


def lo(x):
    return sum(Q(x)**j / factorial(j) for j in range(65))


def hi(x):
    x = Q(x)
    return lo(x) + x**65 / factorial(65) / (1 - x/66)


for logarithm, order in [(Q(13, 5), 13), (Q(7, 2), 32),
                         (Q(21, 5), 66), (Q(9, 2), 82)]:
    assert lo(logarithm) > order

groups = [(13, 31, Q(13, 120), Q(1, 4)),
          (32, 65, Q(7, 124), Q(11, 50)),
          (66, 82, Q(21, 650), Q(203, 1000))]
phases = []
for lower, upper, delta, ratio in groups:
    sigma = 6 + 5*(Q(1, 36) - delta)

    def envelope(order, index):
        return Q(order) / lo(sigma*(1 - Q(index, order)))

    assert 4/lo(sigma/2) < ratio
    assert hi(1)*envelope(upper, 6)/6 < ratio
    radius = hi(Q(13, 20*lower))
    x = radius*ratio
    head = [(radius*max(envelope(lower, k), envelope(upper, k)))**k
            / factorial(k) for k in range(1, 6)]
    mass = sum(head) + x**6/(1 - x)
    tail_euler = sum(k*head[k-1] for k in range(2, 6))
    tail_euler += x**6*(6/(1 - x) + x/(1 - x)**2)
    assert mass < 1
    assert (1 + mass)/(1 - mass) < lo(Q(13, 20))
    assert 2*(head[0] + 2*tail_euler)/(1 - mass) < 1
    phases.append(((head[0] + tail_euler)/(1 - mass),
                   tail_euler/(1 - mass)))

for (up, down), (next_up, next_down) in zip(phases, phases[1:]):
    assert 2*max(up + next_down, next_up + down) < 1
~~~
