# An all-degree binomial entropy criterion for arbitrary retained shared rows

Use the row curvature \(S\), Gamma curvature \(V\), and \(W=1+V\)
of the [oriented endpoint theorem](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md).
An optional common Gaussian strength \(\tau\ge0\) is included
in both \(S,V\). In each orientation retain any subset of the
shared rows. Let \(s_0\) be the curvature of unretained rows
plus \(\tau\), and let \(Q_\pm(w)\) be the oriented chord
envelopes of that theorem. Let \(q=\max Q_\pm(0)\) and
\(L_0=\max(1,\max_{\rm shared}|d|/\ell)\), with \(d=b-a\)
in an orientation and \(\ell=\min(a,b)\).

For each orientation put
\[
 L_\pm=\frac{s_0}{2}+\sum_{\rm retained}
       \left(\frac{d^2}{2\min(a,a+d/2)}
                              +\frac{d/a-\log(b/a)}2\right).
\]
Let \(L=\max L_\pm\).
One may replace \(L\) throughout by any nonnegative bound
for the two increments
\(\log Q_\pm(w)-\log Q_\pm(0)\le Lw\) on \(0\le w\le1/2\),
provided \(L\le3V/4\).

Define, for \(\delta>-1\),
\[
 E(\delta)=8\log(1+\delta/2)-\frac{8\delta}{2+\delta}
              +\frac{4\delta}{(2+\delta)^2}-\log(1+\delta).
\]

**Theorem.** The consecutive moment polynomials have disjoint
root arguments in all degrees if in both orientations
\[
 \begin{gathered}
 s_0+\sum_{\rm retained}E(d/a)\ \ge\ 8(1-\log2),\\
 q\le1/7776,\quad Wq\le51/7776,\quad L_0q\le5/7776,
                                      \quad Lq\le25/7776,\\
 Q_\pm(1/2)\le1/32,\quad 256Q_\pm(1/256)\le1/16,\quad
                         128eQ_\pm(1/128)\le1/8 .
 \end{gathered}                                                   \tag{1}
\]
A profile satisfying these tests remains certified under addition
of any ordinary rows or any increase of the Gaussian strength.
For \(\tau=0\) it satisfies SFC(2) over every characteristic-zero
field. The inequalities allow arbitrary numbers and orientations
of retained rows and do not require a mirrored pair.

For example, take \((1,5),(1,5),(5,1)\) as retained shared rows
and \((2,3)\) as an unretained row. This asymmetric profile passes
(1) and has \(S=149/15\). None of its relative gaps reaches five.
The interval checks below certify the inequalities without
evaluating any of its terminal polynomials.

## Combining the binomial entropy and the row chord

The stated \(L_\pm\) bounds follow from
\[
 (A(w)/w)'=d^2\int_0^1 u/(a+dwu)\,du
                     \le d^2/(2\min(a,a+d/2))
\]
and the tangent bound
\(G(w)/2\le (d/a-\log(b/a))w/2\) for the concave log chord.
Also \(L_\pm\le3V/4\), since the first term per retained row
is at most \(d^2/(2\ell)\) and the second at most
\(d^2/(4a\ell)\).

Write \(B(w)=\exp(H(w)/w)\),
\(H(w)=-w\log w-(1-w)\log(1-w)\).
For integer \(k=mw\), the binomial theorem gives
\(\binom mk^{1/k}\le B(w)\), by the single term in the
binomial sum with success probability \(w\).
In particular \(B(1/2)=4\) and \(B(w)\le e/w\).
The
[Gamma chord curvature calculation](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#a-gamma-chord-bound-in-both-orientations)
bounds the derivative of \(w(\log Q_\pm)'(w)\) below by
\(\tfrac12(s_0+\sum E(d/a))\) for \(w\le1/2\).
On the other hand
\[
 \frac{d}{dw}\bigl[w(\log B)'(w)\bigr]
     =-\frac{w/(1-w)+\log(1-w)}{w^2}
     \ge -4(1-\log2),
\]
by the positive power series of the numerator.
Thus \(w(\log (BQ_\pm))'(w)\) is nondecreasing.
The maxima of \(B(w)Q_\pm(w)\) on closed intervals in
\((0,1/2]\) are attained at the endpoints. The same holds
for \(B(2/m)Q_\pm(2/m)\) as a function of \(m\ge4\),
since its logarithmic derivative multiplied by \(m\) is
the negative of that increasing quantity at \(w=2/m\).
The oriented chord proof gives the corresponding
endpoint maximum for \(mQ_\pm(1/m)\) on \(m\ge2\).

## Separation in all degrees

We can assume \(S<16\) by the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md),
so \(V\le S(S+2)<288\).
Write \(N=n+1\) and, for \(n\ge2\), \(m=n,N\).
For \(N\le256\), the endpoint tests in (1) and the
preceding extrema give \(c_{m,1}\le1/16\) and
\(c_{m,k}\le8^{-k}\) up to the midpoint.
Indeed the first test also uses \(2Q_\pm(1/2)\le1/16\)
at \(m=2\). For \(k\ge2\), the test uses
\(B(1/2)Q_\pm(1/2)\le1/8\), and the remaining endpoint
is bounded by \(128e Q_\pm(1/128)\) at \(m=256\).
The [short phase argument](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-short-ranges)
excludes a common angle.

For \(N>256\), put \(D=Nq\).
If \(D\le4/25\), use \(L/n\le27/32\) and \(LD/n<323/100000\).
For \(k\le4\), the maximum of
\(D\exp(k\min(27/32,(323/100000)/D))\) is at \(D=4/25\).
In fact the only other candidate is the switch
\(D_0=(323/100000)/(27/32)\), since the derivative on the
upper segment can only change from negative to positive;
the comparison \(D_0e^{4(27/32-(323/100000)/(4/25))}<4/25\)
checks every \(k\le4\).
The first three half coefficients are thus bounded by
\(\big((4/25)e^{(323/100000)k/(4/25)}\big)^k/k!\).
At \(k=4\), use \(B(w)\le e/w\) and
\(e(4/25)e^{4(323/100000)/(4/25)}/4<1/8\).
Together with the midpoint this bounds the remaining terms
by \(8^{-k}\). On the disk \(501/500\) these are precisely
the short phase bounds checked in the
[opposed-gap proof](MULTIVARIATE_OPPOSED_SHARED_GAPS_ENTROPY_INTERLACING.md#separation-at-every-order).

If \(D>4/25\), then \(n\ge1244\).
Put \(K=\lceil32D\rceil+14\le126D\).
All four bounds on \(D/n, WD/n, L_0D/n, LD/n\) in the
[long part of the opposed-gap proof](MULTIVARIATE_OPPOSED_SHARED_GAPS_ENTROPY_INTERLACING.md#separation-at-every-order)
follow from the budgets in (1).
Its bound \(c_{m,k}\le6^{-k}\) for \(k\ge K\) follows by
the coefficient checks at \(K\) and the midpoint, with
exactly the same inequalities.
For a positive starting endpoint the first-factor derivative
is bounded there by \(5/(2n)\) in modulus also here.
Indeed \(5V+S(3+\log N)<1440+16(3+\log N)\);
its ratio to \(n=e^{\log N}-1\) decreases for \(\log N\ge7\)
and is below \(3/2\) at \(7\).
This gives the same bound \(B_m<(103/100)D\) on the
interpolated first coefficient. The differentiated Gamma
and connected logarithm comparison then applies with
exactly the parameters of the opposed-gap proof.
For a private starting endpoint, the ratio bounds there
follow from \(q e^{LK/n}\) and the geometric tail.
Thus both orientations have the individual log bound \(2D\)
and the matched difference bound \(1/100\) on radius \(6/5\).
The same paired radius and phase contradiction applies.
Order one versus two follows from positivity of the even
quadratic on the real axis, preserved by Gaussian averaging.

## Preservation when adding rows

Keep the previously retained rows and add the new curvature
\(\sigma\) to \(s_0\). The entropy inequality only improves,
and each \(Q_\pm(w)\) is multiplied by \(e^{-\sigma(1-w)/2}\).
The bounds on \(Wq\) are preserved, by the
[excess-curvature product estimate](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#preservation-under-arbitrary-added-rows)
and \(W\ge5\). The latter follows already from \(q\le1/7776\)
and \(-\log Q_\pm(0)\le V/2\).
The new largest relative gap satisfies
\(L_{0,\rm new}\le\max(L_0,1+\sigma)\).
Since \((1+\sigma)e^{-\sigma/2}<2\), its product budget
also persists, whether \(L_0\ge2\) or \(L_0<2\).
The new increment bound is \(L+\sigma/2\).
Its product with \(e^{-\sigma/2}\) is nonincreasing if
\(L\ge1\), and otherwise is at most one.
This preserves the fourth budget. Increasing \(\tau\)
obeys the same bounds. Therefore every test persists.

## Exact scalar checks and the example

The new scalar checks use rational Taylor bounds.

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    return sum(F(x)**j/factorial(j) for j in range(81))

def hi(x):
    x=F(x)
    return lo(x)+x**81/(factorial(81)*(1-x/82))

assert hi(7)<1245 and F(1600)/(lo(7)-1)<F(3,2)
cap, g, d=F(27,32),F(323,100000),F(4,25)
assert g/cap*hi(4*(cap-g/d))<d
~~~

The example uses directed Arb intervals.

~~~python
from flint import arb, ctx, fmpq
ctx.prec=180
rows=[(1,5),(1,5),(5,1)]
V=arb(fmpq(97,2))
W=1+V
s0=arb(fmpq(1,3))
def row_log(a,b,w):
    d=b-a
    if w == 0:
        return d-b*(b/a).log()
    x=a+d*w
    A=x*x.log()-(1-w)*a*a.log()-w*b*b.log()
    G=x.log()-(1-w)*a.log()-w*b.log()
    return A/w+G/2
def E(d):
    return 8*(1+d/2).log()-8*d/(2+d)+4*d/(2+d)**2-(1+d).log()
for retained in [rows,[(b,a) for a,b in rows]]:
    retained=[(arb(a),arb(b)) for a,b in retained]
    L=s0/2+sum((b-a)**2/(2*(a if b>a else (a+b)/2))
               +((b-a)/a-(b/a).log())/2 for a,b in retained)
    assert L<22
    assert s0+sum(E((b-a)/a) for a,b in retained)>8*(1-arb(2).log())
    def Q(w):
        w=arb(w)
        return (-s0*(1-w)/2+sum(row_log(a,b,w) for a,b in retained)).exp()
    q=Q(0)
    assert q<arb(fmpq(1,7776)) and W*q<arb(fmpq(51,7776))
    assert 4*q<arb(fmpq(5,7776)) and 22*q<arb(fmpq(25,7776))
    assert Q(fmpq(1,2))<arb(fmpq(1,32))
    assert 256*Q(fmpq(1,256))<arb(fmpq(1,16))
    assert 128*arb(1).exp()*Q(fmpq(1,128))<arb(fmpq(1,8))
~~~
