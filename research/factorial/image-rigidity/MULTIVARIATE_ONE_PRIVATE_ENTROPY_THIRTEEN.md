# All-degree angular separation with a private side of large degree

Suppose ordinary nonnegative integer exponent columns \(a,b\) have
private rows on exactly one side. Write \(P\) for the sum of the
exponents on these rows and
\[
 S_{\rm sh}=\sum_{i:a_i b_i>0}\frac{(b_i-a_i)^2}{\max(a_i,b_i)} .
\]
**Theorem.** If \(2P+S_{\rm sh}\ge13\), the consecutive factorial
moment polynomials have disjoint root arguments in every degree.
This implies ordinary SFC(2) over every characteristic-zero field.
In particular, arbitrary shared rows are allowed when \(P\ge7\).

The [absolute-curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
already applies when \(P+S_{\rm sh}\ge12\). In the rest of the proof
\(P+S_{\rm sh}<12\), so \(2\le P\le11\). The private factors yield
a stronger coefficient bound on both halves. We use it with the
same linked phase estimates as in the absolute-curvature theorem.

## The private factorial estimates

Use positive endpoint normalization, so the normalized full
coefficient at \(k\) in degree \(m\) is \(c_{m,k}\) with
\(c_{m,0}=c_{m,m}=1\). For real \(m\ge13\), define
\[
  f_{m,j}=\frac{\Gamma(1+j)}{\Gamma(1+m)^{j/m}},\qquad
  C_m=\frac{1+\frac12\log m}{m},\qquad
  \delta_m=\frac{\log m}{2(m-1)}.
\]
A private row of exponent \(e\) contributes a factor at most
\(f_{m,j}^e\) to \(c_{m,k}\), with \(j=k\) or \(j=m-k\).
Indeed \(x\mapsto\log\Gamma(1+ex)-e\log\Gamma(1+x)\) is
convex and vanishes at zero. Its second derivative is nonnegative
by grouping the terms \((x+r/e)^{-2}\) of the trigamma series
in blocks of \(e\), each block ending at \(r/e\in\mathbb Z\).
Subtracting the endpoint chord proves the claim.

In each half write \(k\) for the index from its chosen endpoint,
\(1\le k\le m/2\), and put \(w=k/m\). With \(k\) integral,
\[
 \frac{1}{k}\log f_{m,m-k}
 \le -1+C_m+\frac1k\sum_{j=0}^{k-1}-\log(1-j/m)
 \le -1+C_m+\frac23(w-1/m).
\tag{1}
\]
The first inequality follows from the Gamma shift by \(k\) and
\[
 m\log m-m\le\log\Gamma(1+m)\le m\log m-m+1+\tfrac12\log m.
\tag{2}
\]
For real \(m\ge13\), (2) follows from the first signed Stirling
remainder bounded by \(1/(12m)\), since \(2\pi<44/7<e^{19/10}\).
For positive integer \(m\), (2) also follows from the integral and
trapezoidal comparisons for \(\log x\).
[DLMF §5.11(ii)](https://dlmf.nist.gov/5.11.ii)

For the second inequality in (1), use
\(-\log(1-t)\le t+t^2/[2(1-w)]\) for \(0\le t\le w\).
The averaged linear term is \((w-1/m)/2\) and the quadratic term
is at most \((w-1/m)w/[6(1-w)]\le(w-1/m)/6\).
Now \(\delta_m\) decreases for \(m>1\) and \(\delta_m<1/3\) for
\(m\ge13\). The two linear functions of \(w\) in the following
comparison agree at \(1/m\), and the left one has smaller slope:
\[
 -1+C_m+\tfrac23(w-1/m)\le -(1-\delta_m)(1-w).
\tag{3}
\]
We have \(f_{m,k}\le f_{m,m-k}\) also for real \(m\).
In fact \(t\mapsto \psi(1+t)+\psi(1+m-t)\) is symmetric about
\(m/2\) and increases up to that point by the decreasing
trigamma function. Its average over \([k,m-k]\) is at least
its whole-interval average. Integrating this inequality gives
\(\log\Gamma(1+m-k)-\log\Gamma(1+k)
\ge (1-2k/m)\log\Gamma(1+m)\), as needed. At \(k=m/2\)
the comparison is equality.

## Short-phase intervals

For integer \(m\ge13\), the private bound (3) in both
orientations and the [shared discrete curvature bound](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
give
\[
 c_{m,k}\le \binom m k e^{-\sigma k(1-k/m)},\qquad
 \sigma=P(1-\delta_m)+S_{\rm sh}/2
\tag{4}
\]
as an upper envelope in either half. For \(13\le m\le31\),
\(\delta_m\le 13/120\), using \(\log13<13/5\). If \(P\le6\),
the assumption yields \(\sigma\ge 13/2-6(13/120)=117/20\).
If \(P\ge7\), the private contribution already exceeds \(6\).
For \(32\le m\le82\), \(\delta_m<7/124\) using \(\log32<7/2\),
so the same two cases give \(\sigma\ge6\).

Use the short-half argument of the [curvature twelve proof](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#early-degrees-and-the-sign-of-the-first-term)
on the groups \([13,31]\) and \([32,82]\), with exponents
\(117/20\) and \(6\), respectively. The resulting scalar caps
are checked below. For \(k\le5\) bound the coefficient by
\([m e^{-\sigma(1-k/m)}]^k/k!\). The maximum at fixed \(k\)
over each group is at an endpoint by logarithmic differentiation.
For \(k\ge6\), use the binomial entropy bound
\((B(w)e^{-\sigma(1-w)})^k\), where
\(B(w)=\exp([-w\log w-(1-w)\log(1-w)]/w)\).
The derivative of \(w(\log B(w)-\sigma(1-w))'\) is at least
\(\sigma-4(1-\log2)>0\) for \(0<w\le1/2\), by differentiation
of \(\log(1-w)/w=-\int_0^1(1-wt)^{-1}\,dt\).
Thus the maximum is at an endpoint of \([6/m_1,1/2]\)
for the upper group endpoint \(m_1\). Use \(B(w)\le e/w\)
at \(6/m_1\) and \(B(1/2)=4\), giving tail ratios \(11/50\)
and \(1/5\) in the two groups.

Split the polynomial into halves as \(H_m^+(z)+z^m H_m^-(1/z)\),
giving each half half the middle coefficient when \(m\) is even.
At radius \(R_m=e^{13/(20m)}\) let \(A\) bound the absolute
coefficient norm of \(H_m^\pm-1\), \(F\) bound the first
coefficient times \(R_m\), and \(T\) bound the absolute Euler
norm of the higher terms. The check \((1+A)/(1-A)<R_m^m\)
confines the root moduli to \((1/R_m,R_m)\) by Rouché in both
orientations. The two half arguments at a root of argument
\(0<\theta<\pi\) have sum \(\phi_m\), with each half argument
between \(-T\theta/(1-A)\) and \((F+T)\theta/(1-A)\).
Here one evaluates both halves on the upper ray at reciprocal
radii and uses complex conjugation on the reciprocal evaluation
in the root equation. The bounds follow because the first sine
is positive, \(|\sin(k\theta)|\le k\theta\) and the real part is
at least \(1-A\). The two within-group checks and the adjacent
group check give \(|\phi_N-\phi_n|<\theta\). The root equations
would force \(\theta-(\phi_N-\phi_n)\in2\pi\mathbb Z\),
which is impossible. This covers \(14\le N\le82\).


## Matching the long-phase hypotheses

For \(N\ge83\) use the oriented function \(Q(w)\) of the
[curvature twelve proof](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#setup-and-a-geometric-coefficient-prefactor),
with all shared rows retained and the private contribution in
its definition set to \(s_0=2P-1\). This is a parameter in
the coefficient envelope. The actual Gamma factors in the
polynomial remain the original ones. Put
\[
 S'=S_{\rm sh}+2P-1,\qquad
 V_d'=2P-1+\sum_{\rm sh}\frac{d_i^2}{a_i},\qquad
 W_d'=1+V_d',\qquad q=Q(0).
\]
Here \(a_i\) is the starting exponent of the chosen orientation.
Now \(S'\ge12\), \(V_d'\ge S'\), and \(V_d'\le
V+P-1<178<195\), where \(V=P+\sum_{\rm sh}d_i^2/\min(a_i,b_i)\)
is the actual bound. For \(V\) use \(V\le S(S+2)\) with
\(S=P+S_{\rm sh}<12\). The rowwise bound follows from
\(|d_i|/\ell_i\le S_i+1\) for \(\ell_i=\min(a_i,b_i)\ge1\)
and \(S_i=d_i^2/(\ell_i+|d_i|)\), then summing the resulting
\(d_i^2/\ell_i\le S_i(S_i+2)\). The actual private degree
satisfies \(P<13\).

The inequalities for \(Q(w)\), including its directed increment,
the endpoint maximum of \(eQ(w)/w\), the midpoint bound and the
cubic endpoint budgets, hold with \(S',V_d',W_d'\).
These are inequalities for a real parameter \(s_0\ge0\) in the
definition of \(Q\). The shared-row calculations are identical,
and \(s_0\) enters those inequalities only through its explicit
chord term. The applicable budgets are
\[
 \begin{split}
 q&\le e^{-6},\qquad W_d' q\le 13e^{-6},\qquad
 (W_d')^3q\le 13^3e^{-6},\\
 L_d q&\le (9/8)e^{-6},\qquad Q(1/2)\le e^{-3},\qquad
 L_d\le\sqrt{V_d'} ,
 \end{split}
\tag{5}
\]
using \(L_d=\max(1,\max_{\rm sh}|d_i|/a_i)\) and the
decrease of the budgets in \(S'\ge12\).
[Directed increment and endpoint budgets](MULTIVARIATE_ABSOLUTE_CURVATURE_13.md#directed-endpoint-estimates)

For real \(m\ge82\), \(P\delta_m<1/2\) since
\(P\le11\), \(\log82<5\) and \(\delta_m\) decreases.
Consequently (3) gives the factor
\(\exp(-s_0 k(1-w)/2)\) in either orientation. Together with
the shared-row estimates this supplies the
\(c_{m,k}\le[mQ(w)]^k/k!\) hypothesis, including the
real-order coefficients used in differentiating the heads.
Use \(k\) integral here, as the head differences involve only
coefficients at integral indices.

For the extra bound at the private starting endpoint, reserve one
factor \(f_{m,k}\) and use
\[
   \bigl[\binom m k f_{m,k}\bigr]^{1/k}
     =\frac{[(m)_{\underline{k}}]^{1/k}}{\Gamma(1+m)^{1/m}}
     \le e .
\]
The remaining \(P-1\) factors obey
\[
  \frac1k\log f_{m,k}\le -(1-w)-\frac{1/2-w}{P-1}.
\tag{6}
\]
Indeed, (2) gives \(\log f_{m,k}/k\le\log w+C_k\), with
\(C_k=(1+\frac12\log k)/k\) and the integer upper estimate
also valid at \(k=1\). To verify (6) for \(k\le6\), use
\(C_k\le1,\ w\le6/82,\ P-1\ge1\), giving an upper
bound \(\log(6/82)+5/2<0\) after moving its right
side to the left. For \(7\le k\le11\), use
\(C_k<2/7\) and \(w\le11/82\) to bound the same difference
by \(\log(11/82)+25/14<0\). For \(k\ge12\), \(C_k<3/16\).
The function \(\log w+1-w+(1/2-w)/(P-1)\) increases on
\(0<w\le1/2\) and takes its maximum \(1/2-\log2\)
there, which is less than \(-3/16\).
Thus the reserved factor and (6) yield
\[
 \bigl[\binom m k f_{m,k}^P\bigr]^{1/k}
      \le e^{3/2} e^{-P(1-w)} .
\]
This supplies \(c_{m,k}^{1/k}\le e^{3/2}Q(w)\) with the same
\(s_0=2P-1\) as in the other orientation.

The rest of the [direct](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#direct-estimates-in-one-orientation),
[connected](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#connected-estimates-in-one-orientation),
and [joining](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#joining-direct-and-connected-estimates)
proofs now applies. The range cut \(D=Nq=2/5\), the increment
bound, the coefficient tails and the connected bounds use (5)
and \(V_d'<195\). The exact positive-start logarithmic derivative
of the first coefficient uses the actual \(V<195\) and actual
private degree \(P<13\). It gives the same \(2/n\) and \(1/n\)
bounds in the direct and connected ranges, respectively.
The polygamma derivative and finite-difference bounds use the
slopes of the actual Gamma arguments. Their quadratic budgets
are bounded above by \(W_d'\), since \(2P-1\ge P\), and their
relative slope bound is the same \(L_d\). A private ending has
relative slope one. The extra factor \(e^{13k/328}\) in the
direct real-order head bound is harmless since the coefficient
bound above already holds on that real-order interval.
Thus the unchanged scalar checks in the known proof give the
same half-log and paired-phase estimates. This covers \(N\ge83\).

For terminals five through thirteen use the
[angular terminal results](MULTIVARIATE_ALL_DEGREE_REDUCTION.md#the-fixed-terminal-obligation).
At terminal three the quadratic has no real zero.
A shared nonreal argument gives divisibility of the
quadratic into the cubic after independent positive rescaling,
which the strictly decreasing scale-invariant \(E_m\) forbids.
The terminal-four contradiction likewise uses only
scale-invariant coefficient bounds, quartic real positivity
and the common unit conjugate factor obtained by these
rescalings. At terminal two the linear root is negative real
and the even moment has no real root.
[Small-terminal coefficient proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)

On the real rays, positive coefficients exclude positive roots
and the even-order Gamma moment is strictly positive on the
real line. Conjugation covers the lower half plane. Over a
characteristic-zero field the ordinary integral coefficients
give a rational Bézout identity for each pair.

## Exact scalar replay for the new estimates

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    x = F(x)
    return sum(x**j/factorial(j) for j in range(65))

def hi(x):
    x = F(x)
    return lo(x) + x**65/factorial(65)/(1-x/66)

assert lo(F(13, 5)) > 13
assert lo(F(7, 2)) > 32
assert lo(5) > 82
assert lo(F(19, 10)) > F(44, 7)
assert 11 * F(5, 2*81) < F(1, 2)
assert hi(F(5, 2)) < F(82, 6)
assert lo(2) > 7
assert hi(F(25, 14)) < F(82, 11)
assert lo(F(5, 2)) > 12
assert F(56, 81) > F(11, 16)  # two positive terms of log 2

short = []
for m0, m1, sigma, r in [(13, 31, F(117, 20), F(11, 50)),
                         (32, 82, F(6), F(1, 5))]:
    def envelope(m, k):
        return F(m)/lo(sigma*(1-F(k,m)))

    assert 4/lo(sigma/2) < r
    assert hi(1)*envelope(m1,6)/6 < r
    R = hi(F(13,20*m0))
    x = R*r
    A = first = T = F(0)
    for k in range(1,6):
        start = max(2*k,m0)
        c = (R*max(envelope(start,k),envelope(m1,k)))**k/factorial(k)
        A += c
        if k == 1:
            first = c
        else:
            T += k*c
    A += x**6/(1-x)
    T += x**6*(6/(1-x)+x/(1-x)**2)
    assert A < 1 and (1+A)/(1-A) < lo(F(13,20))
    assert 2*(first+2*T)/(1-A) < 1
    short.append(((first+T)/(1-A), T/(1-A)))

(up, down), (next_up, next_down) = short
assert 2*max(up+next_down, next_up+down) < 1
print("Private endpoint and short-phase scalar checks passed.")
~~~
