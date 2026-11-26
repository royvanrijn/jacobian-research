# Angular separation from the joint matched private-unit packet equations

Take an ordinary factorial binomial with exactly one private unit row in
each orientation and arbitrary shared rows. For each unequal shared row
write \(d_i=|b_i-a_i|\) and \(\ell_i=\min(a_i,b_i)\), and set
\[
 S=\sum_i\frac{d_i^2}{\ell_i+d_i},\qquad V=\sum_i\frac{d_i^2}{\ell_i}.
\]
In matched tanh coordinates, \(S=\sum_{\rm slots}2v/(1+v)\) and
\(V=\sum_{\rm slots}2v/(1-v)\), where
\(v=(n/T_{n+1}-(n-1)/T_n)^{-1}=d_i/(2\ell_i+d_i)\).
[Exact matched coordinates](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

**Theorem.** At adjacent orders \(n,n+1\), the moments have disjoint
root arguments in each of the following ranges. Positive rescalings may
be chosen independently.

- For \(n\ge13\), there is no shared argument with
  \(b=n|\sin\theta|\ge8\) and \(b^3\ge62n\).
- For every \(n\ge3\), \(V\le 1/n\) excludes every shared argument. For
  \(n\ge13\) it suffices to have \(V\le 13/(5n)\).

The first range includes every upper angle box \(k\) with
\(4k-2\ge\max\{8,(62n)^{1/3}\}\). The pairs at \(n=1,2\) separate
in argument by the angle boxes and strict positivity of the even
moment on the real axis. In particular, if \(J\ge14\) is at least
the [absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
then \(V\le\min\{1/12,13/[5(J-2)]\}\) is a sufficient single
matched-moment condition for all-degree angular separation.
The numerical value of \(J\) is unevaluated. The two local
comparisons below have no eventual cutoff.

## Eliminate the radii in the complex packet equation

Write \(N=n+1\). The [angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
place a shared argument in the same label after conjugation, with
\(0<\theta<\pi/2\) and \(\theta>(2k-1)\pi/n\).
At \(S=0\) the normalized moments are geometric sums. At \(S\ge15/8\)
the [curvature theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
gives disjoint arguments. Take \(0<S<15/8\). Ordinary
[row packing](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md)
gives \(V<V_*:=61/12\) and the normalized log root radius satisfies
\(|x_m|\le11/m\) for \(m\ge4\).
The orientation budgets \(M_+=\sum_i d_i^2/a_i\) and
\(M_-=\sum_i d_i^2/b_i\) satisfy \(M_\pm\le V\) and
\(M_++M_-=V+S\).

Let \(E_{m,\pm}(x)=\mathbb E(1-e^{x+L_{m,\pm}+i\theta})^{-1}\)
be the geometric endpoint packets, where \(L_{m,\pm}\) is the
log of the corresponding endpoint variable. Each packet has
argument in \((0,\pi-\theta)\). Set
\[
 H_m(x)=\log E_{m,+}(x)-\overline{\log E_{m,-}(-x)}.
\]
The packet equations at the root with log radius \(x_m\) give
\(m x_m=\Re H_m(x_m)\) and
\(m\theta=(2k-1)\pi+\Im H_m(x_m)\). The
[resolvent derivative bound](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#resolvent-variation-at-any-horizontal-position)
gives \(|H_m'|\le2/\sin\theta\) for every real \(x\).
[Exact endpoint packet equations](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)

Suppose the roots share an argument. Write
\(\delta=x_N-x_n\),
\(D=H_N(x_N)-H_n(x_N)\), and
\(A=\int_0^1 H_n'(x_n+t\delta)\,dt\). Then
\[
       N\delta+x_n+i\theta=D+A\delta .
\]
If \(|D|\le E\) and \(\beta=2/(N\sin\theta)<1\), set
\(y=N\delta+x_n\). The joint modulus bound is
\(\sqrt{y^2+\theta^2}\le E+\beta|x_n|+\beta|y|\).
By Cauchy–Schwarz,
\(\sqrt{y^2+\theta^2}\ge\beta|y|+\theta\sqrt{1-\beta^2}\).
Hence every shared ray satisfies
\[
               \theta\sqrt{1-\beta^2}\le E+\beta|x_n|.        \tag{1}
\]
The radius difference is unrestricted in this bound.

## Matched log increment

Use the independent adjacent Beta log noise \(Z_\pm\) in the
two orientations. The [linked endpoint comparison](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md#the-two-orientation-budgets)
gives nonnegative means \(\mu_\pm\) satisfying
\(\mu_++\mu_-\le S/(2nN)\) and \(\mu_\pm\le M_\pm/(2nN)\).
Use \(\alpha=1/\sin\theta\).
Suppose the sum of the second exponential moments
\(\mathbb E[Z_\pm^2 e^{\alpha Z_\pm}]+\mathbb E[Z_\pm^2 e^{-\alpha Z_\pm}]\)
is at most \(P M_\pm/(nN)\) in each orientation. The
[linked resolvent Taylor estimate](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#resolvent-variation-at-any-horizontal-position)
then bounds the relative change of each packet at any fixed
real shift by
\[
     \epsilon_\pm=\frac{\mu_\pm}{\sin\theta}
                +\frac{P M_\pm}{nN\sin^2\theta}.
\]
If \(\epsilon_\pm\le e_*<1\), the difference between packet
logarithms is at most \(-\log(1-\epsilon_\pm)\), by integrating
the reciprocal along the segment joining the two packet values.
The segment stays in the packet sector, so the bound uses their
original logarithm branches. Since \(n/N<1\) and \(\sin\theta\le1\),
one may use
\[
     E=\frac{S/2+P(V+S)}{(1-e_*)b^2},\qquad b=n\sin\theta,    \tag{2}
\]
in (1).

For \(n\ge13,\ b\ge8\), use \(P=18/5\).
Equation (6) of the
[endpoint exponential bound](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md#explicit-exclusion-at-a-common-angle)
gives this constant with \(n_0=13,\ b_0=8\), using the two row
slope caps \(2/3,2\).
The replay verifies the required rational inequality.
Put \(S_*=15/8\) and take \(e_*=(S_*/2+PV_*)/64\), using
\(\mu_\pm\le S/(2nN)\) separately for the maximum.
From (1)–(2), \(\beta\le2/b\), \(|x_n|\le11/n\), and \(b\le n\),
\[
 b^2\theta
 \le \frac{4}{\sqrt{15}}
       \left(\frac{S_*/2+P(V_*+S_*)}{1-e_*}+22\right)<62 .
\]
The strict inequality follows by squaring rational quantities
as below. Since \(\theta>\sin\theta=b/n\), a shared argument
forces \(b^3<62n\). The box inequality and concavity of the sine
on \([0,\pi/2]\) give \(b>2(2k-1)\), proving the label claim.

## Weak curvature in every angle box

For the weak bounds, \(|d_i|/\ell_i\le V\) by integrality of the
nonzero \(d_i\). Thus both slope caps in the linked noise
bound can be replaced by \(V\). For lower bounds \(n_0,b_0\)
on \(n,b\), assume \(V\le v_0\), put \(q=v_0/b_0\), and set
\[
 \Lambda_*=\frac{v_0}{2(n_0+1)b_0}+
                      \frac{v_0}{2b_0^2(1-q)}.
\]
If \(q,\Lambda_*<1\), the exponential bound gives
\[
 P\ \ge\
 \frac{2}{1-\Lambda_*}\left((1-q)^{-2}
       +v_0\left(\frac1{2n_0}+\frac1{b_0(1-q)}\right)^2\right).
\]
Indeed \(M_\pm\le V\), and (4)–(6) of the endpoint comparison
apply in both signs with this same \(q\); the positive mean
bound in \(\Lambda_*\) may be used even at the negative argument.
Take \(e_*=(P+1/2)v_0/b_0^2\) since
\(\mu_\pm\le M_\pm/(2nN)\). The
[weak chord estimate](MULTIVARIATE_PRIVATE_UNIT_WEAK_THREE_TENTHS.md#modulus-estimate-with-the-weak-row-cap)
and the Eneström–Kakeya bound give \(|x_n|\le V/2\).
With \(S\le V\) we obtain from (1)–(2)
\[
 \theta\le\frac{V}{\sqrt{1-(2/b_0)^2}}
        \left(\frac{1/2+2P}{(1-e_*)b_0^2}+\frac1{b_0}\right).
\]
When \(V\le c/n\) take \(v_0=c/n_0\). The following two
choices make the right side strictly smaller than \(\pi/n\).

| \(n_0\) | \(b_0\) | \(c\) | \(P\) |
| ---: | ---: | ---: | ---: |
| \(3\) | \(5/2\) | \(1\) | \(16/5\) |
| \(13\) | \(31/10\) | \(13/5\) | \(12/5\) |

In both cases \(b\ge b_0\) follows from
\(n\sin(\pi/n)\ge\pi-\pi^3/(6n^2)\), using the angle boxes.
The strict upper bound contradicts \(\theta>\pi/n\).
This excludes every candidate label in the stated weak ranges.

## Exact scalar replay

~~~python
from fractions import Fraction as Q

n0, b0 = 13, Q(8)
V, S, P = Q(61, 12), Q(15, 8), Q(18, 5)
a = S / (2 * (n0 + 1) * b0)
bound = 0
for slope, positive in ((Q(2, 3), True), (Q(2), False)):
    q = slope / b0
    ell = (a if positive else 0) + V / (2*b0**2*(1-q))
    assert 0 < q < 1 and 0 < ell < 1
    square = V * (Q(1, 2*n0) + 1/(b0*(1-q)))**2
    bound += ((1-q)**-2 + square) / (1-ell)
assert bound < P
cap = (S/2 + P*V) / b0**2
assert 0 < cap < 1
H = (S/2 + P*(V+S)) / (1-cap)
assert 16*(H+22)**2 < 15*62**2

for n0, b0, c, P in (
    (3, Q(5, 2), Q(1), Q(16, 5)),
    (13, Q(31, 10), Q(13, 5), Q(12, 5)),
):
    assert b0 < Q(157, 50) - Q(22, 7)**3 / (6*n0**2)
    v0 = c/n0
    q = v0/b0
    ell = v0/(2*(n0+1)*b0) + v0/(2*b0**2*(1-q))
    assert 0 < q < 1 and 0 < ell < 1
    square = v0 * (Q(1, 2*n0) + 1/(b0*(1-q)))**2
    assert 2*((1-q)**-2 + square)/(1-ell) < P
    cap = (P+Q(1, 2))*v0/b0**2
    assert 0 < cap < 1
    H = (Q(1, 2)+2*P)/((1-cap)*b0**2) + 1/b0
    assert (c*H)**2 < Q(157, 50)**2 * (1-(2/b0)**2)
~~~
