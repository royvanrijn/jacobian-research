# All higher angle labels for two private unit rows

Consider an ordinary factorial binomial with exactly one
private unit row on each side and arbitrary shared rows.
At adjacent degrees \(n,N=n+1\), a common physical nonreal
root has the same positive angle-box label \(k\) in the two
degrees after conjugation if necessary.

**Theorem.** A common physical root cannot have label \(k\ge2\)
in any degree. This holds for every orientation and signed
degree difference of the shared rows. The analytic comparison
covers every \(n\ge13\). Orders through terminal thirteen use
the [existing ordinary terminal theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md).
Any remaining common rational factor is a single irreducible
quadratic, with its upper root in the first angle box.

The proof strengthens the
[linked endpoint noise comparison](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md).
The two noise laws share one curvature budget. Only the
lower-order log radius enters the linked modulus comparison,
so angle-dependent radial walls supply the remaining margin
at small and medium degrees.

## Geometry and the selected radial walls

Use \(d_i=b_i-a_i\), \(\ell_i=\min(a_i,b_i)\) over unequal
shared rows and put
\[
 S=\sum_i\frac{d_i^2}{\ell_i+|d_i|},\qquad
 M_+=\sum_i\frac{d_i^2}{a_i},\qquad
 M_-=\sum_i\frac{d_i^2}{b_i},\qquad
 V=\sum_i\frac{d_i^2}{\ell_i}.
\]
The [curvature fifteen-eighths theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
handles \(S\ge15/8\). At \(S=0\) the normalized moments are
geometric sums and consecutive argument sets are disjoint.
For \(0<S<15/8\) the
[integer packing bound](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md#shared-curvature-and-the-endpoint-laws)
gives
\[
 \frac{|d_i|}{\ell_i}\le2,\qquad
 V<V_*=\frac{61}{12},\qquad S<S_*=\frac{15}{8}.
 \tag{1}
\]
Each \(M_\pm\le V\), and the sum satisfies \(M_++M_-=V+S\).

After canceling the private unit factorials and scaling the
endpoints to one, the coefficients are moments of a positive
law \(\mu_m\) with \(\mathbb E T^m=1\); the reverse law
\(\mu_m^\vee\) is obtained by order-\(m\) tilt followed by
inversion. The two
[geometric endpoint packets](MULTIVARIATE_PRIVATE_UNIT_ABSOLUTE_SIGNED_IMBALANCE.md#endpoint-packets-on-the-two-radial-boundaries)
are expectations of \((1-rTe^{i\theta})^{-1}\) and its
reverse version at reciprocal radius and conjugate angle.
Both laws have \(\mathbb E T\le1\).

At \(0<\theta\le\pi\), the kernel's argument belongs to
\([0,\pi-\theta]\): its numerator after rationalization is
\(1+rT e^{i(\pi-\theta)}\), and its denominator is positive.
Rotating to the middle of this sector and applying Jensen
therefore gives the lower packet bound
\[
 \left|\mathbb E(1-rTe^{i\theta})^{-1}\right|
     \ge\sin(\theta/2)\,\mathbb E(1+rT)^{-1}
     \ge\frac{\sin(\theta/2)}{1+r}.
 \tag{2}
\]
The upper bound is \(1/\sin\theta\) at \(\theta\le\pi/2\)
and one at \(\theta\ge\pi/2\). These bounds remain valid
under mixtures of the endpoint laws.

At any normalized root of degree \(m\) and label \(k\ge2\),
write \(x_m\) for the log radius. The
[angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#consecutive-labels-and-the-remaining-comparison)
give \(\theta>3\pi/m\). We claim \(|x_m|\le L/m\) with
\(L=3\) for \(13\le m\le18\), and \(L=4\) for
\(18\le m\le30\). To prove this, use the two radial walls
\(r=e^{\pm L/m}\).
For \(13\le m\le18\), \(\theta\ge\pi/6\), so (2) and
the upper packet bound give ratio at most
\(8(1+e^{3/13})\le8(1+13/10)<e^3\).
For \(18\le m\le30\), \(\theta\ge\pi/10\), and the
corresponding ratio is at most
\[
          (200/9)(1+e^{4/18})\le 3200/63<e^4 .
\]
Here \(\sin(\pi/10)>3/10\) and
\(\sin(\pi/20)>3/20\) follow from the cubic sine bound;
the exponential upper estimates use \(e^u\le(1-u)^{-1}\).
At higher angles the upper packet bound one only improves
these comparisons. Thus the two endpoint packets in the
finite geometric sum cannot cancel on either wall, since
the factor on the reverse packet has modulus \(e^{\pm L}\).
Along the linear interpolation from \(\mu_m\) to the
reciprocal mixture \((\mu_m+\mu_m^\vee)/2\), the law is
nonconstant with normalized endpoints. The single simple
root in each open angle box varies continuously and keeps
its label. At the reciprocal reference it is on the unit
circle. The wall exclusion therefore traps that root in the
indicated annulus. In larger degrees use the
[absolute normalized radius bound](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md),
which gives \(|x_m|\le11/m\).

## The two orientation budgets

Let \(G_m\) be the physical argument normalization and
\(\Delta=\sum_i d_i\). A common physical root has normalized
radii satisfying \(x_N-x_n=\log(G_N/G_n)\) and a shared
argument \(0<\theta<\pi/2\). Put
\(J(y)=\psi(Ny+1)-\psi(ny+1)\), where \(\psi\)
is the logarithmic derivative of the Gamma function.
The [exact adjacent Beta links](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#exact-signed-log-noise)
give independent log noises \(Z_+,Z_-\) for the two
orientations. In either orientation,
\(\log T_n\) has the law of \(\log T_N+Z_\pm\).
Write \(\mu_\pm=\mathbb E Z_\pm\). The signed mean formula
and \(0\le J'(y)\le1/(2nNy^2)\) yield
\[
 \begin{aligned}
 \mu_+&=\sum_i d_i^2\int_0^1(1-v)J'(a_i+d_i v)\,dv,\\
 \mu_-&=\sum_i d_i^2\int_0^1 vJ'(a_i+d_i v)\,dv .
 \end{aligned}
\]
The weights exchange on reversing the row. Thus
\[
 0\le\mu_\pm,\qquad
 \mu_++\mu_-\le\frac1{2nN}\sum_i\frac{d_i^2}{a_i b_i}
                   \le\frac{S}{2nN},\qquad
 \mu_\pm\le\frac{M_\pm}{2nN}.                                  \tag{3}
\]
The row integral of \((a_i+d_i v)^{-2}\) is \(1/(a_i b_i)\).

For either orientation write \(M_0\) for its \(M_+\) or
\(M_-\), \(c_0=M_0/(nN)\), and let
\(\Lambda(t)=\log\mathbb E e^{tZ_\pm}\).
It is well defined for \(|t|<n/2\).
The beta integral, followed by integrating the difference
of trigammas in the order, gives the sharper one-sided bound
\[
 0\le\Lambda''(t)\le\frac{c_0}{(1-q)^2},\qquad
 q=\begin{cases}(2/3)|t|/n,&t\ge0,\\ 2|t|/n,&t\le0.\end{cases}
 \tag{4}
\]
For the positive side only decreasing exponent rows can
decrease the shifted denominator, and \(|d_i|/a_i\le2/3\)
for such a row by (1). The same reasoning applies after
reversal. On the negative side use \(|d_i|/a_i\le2\).
The shifted shapes stay positive in the indicated open
domain. The derivative estimates follow by integration:
\[
 |\Lambda'(t)|\le c_0\left(\frac12+\frac{|t|}{1-q}\right),
 \qquad
 \Lambda(t)\le t\mu_\pm+\frac{c_0t^2}{2(1-q)} .
 \tag{5}
\]
For the second bound the integral of the second derivative
majorant costs \(-q-\log(1-q)\le q^2/[2(1-q)]\) after
the change of variable. The linear term is nonpositive
when \(t\le0\).

## Explicit exclusion at a common angle

Put \(b=n\sin\theta\) and \(\alpha=1/\sin\theta\).
On one interval of \(n\) take \(n_0\le n\), \(b_0\le b\)
and the corresponding radius constant \(L\).
For the two signs set \(q_+=(2/3)/b_0\) and \(q_-=2/b_0\).
By (4)–(5) the sum
\[
                \mathbb E Z_\pm^2 e^{\alpha Z_\pm}
                   +\mathbb E Z_\pm^2 e^{-\alpha Z_\pm}
\]
is at most \(P M_\pm/(nN)\) whenever
\[
 \begin{split}
       P &\ge \sum_{\sigma=\pm}
         \frac{(1-q_\sigma)^{-2}
          +V_*\bigl(1/(2n_0)+1/[b_0(1-q_\sigma)]\bigr)^2}{1-L_\sigma},
       \\
       L_\sigma&=\mathbf 1_{\sigma=+}\frac{S_*}{2(n_0+1)b_0}
                            +\frac{V_*}{2b_0^2(1-q_\sigma)} .
 \end{split}                                                    \tag{6}
\]
In fact the second derivative of the exponential moment is
\(e^{\Lambda(t)}(\Lambda''(t)+\Lambda'(t)^2)\).
Use \(e^u\le(1-u)^{-1}\) at the positive upper bounds
\(L_\sigma<1\).

For the packet difference at a common real log radius,
the [linked resolvent Taylor inequality](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#resolvent-variation-at-any-horizontal-position)
bounds the relative error in orientation \(\pm\) by
\(\mu_\pm/\sin\theta+P M_\pm/(nN\sin^2\theta)\).
In particular it is at most \(a+cM_\pm\) with
\[
                 a=\frac{S_*}{2(n_0+1)b_0},\qquad
                 c=\frac{P}{b_0^2}.
\]
In each tier \(a+cV_*<1\).
The log error is bounded by the increasing convex function
\(v\mapsto-\log(1-a-cv)\). Since
\(M_++M_-\le V_*+S_*\) and each \(M_\pm\le V_*\),
its two-orientation sum is bounded by
\[
             U(a+cV_*)+U(a+cS_*),\qquad
             U(y)=y+y^2/2+y^3/[3(1-y)].
 \tag{7}
\]
Indeed monotonicity first permits filling the total budget,
and convexity then permits moving all available mass to
one orientation until it reaches \(V_*\).
The formula for \(U\) bounds the logarithm by its power
series.

For the physical drift the elementary inequality
\(h(u)=u/(e^u-1)\ge1-u/2\) in the trigamma integral gives
\[
                  0\le 1/m-y\psi'(my+1)\le1/(2m^2y).
\]
Integrate from \(n\) to \(N\) and along each oriented row.
The private unit contributions to the endpoint scale cancel.
Since \(\sum_i|\log(b_i/a_i)|\le V\), this gives the
half-sized center error
\[
             \left|\log(G_N/G_n)-\Delta\log(N/n)\right|
                               \le V/(2nN).
 \tag{8}
\]
Now subtract the two modulus equations as in the
[linked packet comparison](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#subtracting-the-two-modulus-equations).
The log-radius derivative is at most \(2/\sin\theta\),
and \(N-2/\sin\theta\ge (1-2/b_0)N\).
Since \(\log(N/n)>1/N\), (7)–(8) give for nonzero
\(\Delta\) the strict bound
\[
 |\Delta|<
  \frac{b_0}{b_0-2}
       \bigl(L/n_0+U(a+cV_*)+U(a+cS_*)\bigr)+\frac{V_*}{2n_0}.
 \tag{9}
\]
Use the following rational data. The lower bound on \(b\)
follows from \(\theta>3\pi/n\), \(\theta<\pi/2\), and
\(\sin u\ge u-u^3/6\).

| Lower order \(n\) | \(n_0\) | \(b_0\) | \(L\) | \(P\) | Strict upper bound in (9) |
| --- | ---: | ---: | ---: | ---: | ---: |
| \(13\)–\(17\) | \(13\) | \(859/100\) | \(3\) | \(17/5\) | \(99/100\) |
| \(18\)–\(29\) | \(18\) | \(449/50\) | \(4\) | \(13/4\) | \(17/20\) |
| \(n\ge30\) | \(30\) | \(46/5\) | \(11\) | \(16/5\) | \(47/50\) |

Thus (9) forces the integer \(\Delta\) to vanish.
The [balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
excludes a physical common root when \(\Delta=0\).
The only possible common zeros left are the first upper
root and its conjugate, since the even polynomial has no
real zeros. The angle boxes give simplicity, so a nontrivial
rational gcd must be quadratic and irreducible.

## Exact scalar replay

```python
from fractions import Fraction as Q

V, S = Q(61, 12), Q(15, 8)

def log_upper(x):
    assert 0 <= x < 1
    return x + x*x/2 + x**3/(3*(1-x))

tiers = (
    (13, Q(859, 100), 3, Q(17, 5), Q(99, 100)),
    (18, Q(449, 50), 4, Q(13, 4), Q(17, 20)),
    (30, Q(46, 5), 11, Q(16, 5), Q(47, 50)),
)
for n0, b, L, P, bound in tiers:
    assert 3*Q(157, 50)-(3*Q(22, 7))**3/(6*n0*n0) > b
    a = S/(2*(n0+1)*b)
    second = 0
    for slope, positive in ((Q(2, 3), True), (Q(2), False)):
        q = slope/b
        ell = (a if positive else 0)+V/(2*b*b*(1-q))
        assert 0 <= q < 1 and 0 <= ell < 1
        prime_sq = V*(Q(1, 2*n0)+1/(b*(1-q)))**2
        second += ((1-q)**-2+prime_sq)/(1-ell)
    assert second < P
    c = P/(b*b)
    drift = b/(b-2)*(Q(L,n0)+log_upper(a+c*V)+log_upper(a+c*S))
    drift += V/(2*n0)
    assert drift < bound < 1

assert 8*(1+Q(13,10)) < Q(27,10)**3
assert Q(3200,63) < Q(27,10)**4
for d in (10,20):
    assert Q(157,50)/d-(Q(22,7)/d)**3/6 > Q(3,d)
```
