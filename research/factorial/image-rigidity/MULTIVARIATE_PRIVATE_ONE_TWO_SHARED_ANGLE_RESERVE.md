# Angular separation with private totals one and two and a shared reserve

Let the ordinary two-monomial profile have private exponent totals
one and two, after exchanging the two columns if needed. Put
\[
 S_{\rm sh}=\tau+\sum_{\rm shared}
                    \frac{(b_i-a_i)^2}{\max(a_i,b_i)},
\]
including any common normalized Gaussian strength \(\tau\ge0\).

**Theorem.** Consecutive moment polynomials have disjoint root arguments
at every terminal \(N\ge14\) under either of the following conditions.

- If the total two comes from one private row of exponent two, assume
  \(S_{\rm sh}\ge13/20\).
- If the total two comes from two private unit rows, assume
  \(S_{\rm sh}\ge83/80\).

For ordinary profiles the
[complete terminal-thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
gives all-degree angular separation and SFC(2) under these conditions.
The shared orientations and signed total gap are arbitrary. The new
proof tests scalar envelopes only, with a uniform tail bound.

The exact matched tanh slots give the normal form directly. A private
slot has \((n/T'-(n-1)/T)^{-1}=1\). A shared slot with
\(u=(n/T'-(n-1)/T)^{-1}\) contributes \(2u/(1+u)\) to
\(S_{\rm sh}-\tau\), since its row of gap \(d\) contains \(d\) such
slots and \(u=d/(2\min(a_i,b_i)+d)\).
[Exact matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

## The exact private residual and its two envelopes

Orient the side with total two toward index \(k\). The opposite
private unit provides \((m-k)!\). The two possibilities on the
chosen side supply \((2k)!=4^k k!(1/2)_k\) or \((k!)^2\).
One \(k!\) cancels the binomial together with \((m-k)!\).
The remaining normalized private factor is
\[
                    f_{m,k}^{(a)}=\frac{(a)_k}{(a)_m^{k/m}},
                    \qquad a=\tfrac12\ \text{or}\ 1 .
\]
Use \(s=13/20\) or \(83/80\), respectively, as the lower bound on
\(S_{\rm sh}\). A shared row of gap \(d>0\) and larger endpoint \(b\)
has second log-factorial difference at least \(d^2/(mb)\).
In the increasing orientation, with \(x\) the factorial argument
at an interior index, that difference is
\[
          \sum_{j=1}^{d}\log\frac{x+j}{x-d+j}
                     \ge\sum_{j=1}^{d}\frac{d}{mb}.
\]
Reversal gives the same bound for the decreasing orientation.
Comparing with the quadratic with this second difference and with
the same endpoint chord proves the normalized coefficient bound
\[
                  c_{m,k}\le f_{m,k}^{(a)}
                              e^{-(s/2)k(1-k/m)}.
\]
The Gaussian contribution has exactly this quadratic form.

Take \(R_m=e^{7/(10m)}\) and form the endpoint halves
\(P_m=L_m+z^m U_m(1/z)\). Define
\[
 B_{L,k}=R_m^k f_{m,k}^{(a)}e^{-(s/2)k(1-k/m)},\qquad
 B_{U,k}=R_m^k f_{m,m-k}^{(a)}e^{-(s/2)k(1-k/m)}
\]
for \(1\le k\le m/2\), with half weight at the even middle index.
The descending cap dominates: for
\(\phi(x)=\log\Gamma(a+x)-\log\Gamma(a)-x\log (a)_m/m\),
the function \(\phi(m-x)-\phi(x)\) is concave on \([0,m/2]\)
and vanishes at both endpoints, since the trigamma decreases.

For fixed \(k\), both unhalved caps decrease with \(m\). This
uses the increasing geometric mean \(((a)_m)^{1/m}\) in the
ascending direction. In the descending direction the logarithmic
derivative of the private factor equals
\[
 \frac{k}{m^2}\int_0^m x\psi'(a+x)\,dx
                    -\sum_{h=1}^k\frac1{m+a-h}<0.
\]
Indeed \(a\ge1/2\) and the midpoint integral bound for the
convex function \((x+t)^{-2}\) gives \(\psi'(a+x)<1/x\).
Since \(a\le1\) the sum is at least \(k/m\).
The radial and quadratic factors also decrease at fixed \(k\).

For \(m\ge20\), retain the first nine unhalved caps at \(m=20\)
and use \(2^{-k}\) for \(k\ge10\). At an even midpoint \(m=2h\)
the right cap before halving is
\[
              e^{7/20-sh/4}\frac{(a)_h}{\sqrt{(a)_{2h}}}.
\]
It is below \(2^{-h}\) at \(h=10\), and its consecutive ratio
is at most \(e^{-s/4}(h+a)/(2h+a)<1/2\) thereafter.
The log cap is convex in \(k\) before halving and vanishes at
zero, so the midpoint cap implies the geometric cap for the
whole half. An odd degree uses the previous even unhalved cap
by fixed-index monotonicity.

## Signed phases on the root annulus

Take the following coefficients \(\alpha_H,\beta_H\) for
\(H=L,U\). All entries are exact decimals.
\[
\begin{array}{c|cc|cc}
 a & \alpha_L&\alpha_U&\beta_L&\beta_U\\ \hline
 \tfrac12&.203&.559&.038&.194\\
 1       &.254&.518&.051&.166
\end{array}
\]
Let \(D(\beta)=2/[9\beta(1-(22\beta/7)^2/6)]\).
The cap estimates and exact replay below give
\[
 \sum_kB_{U,k}<\tfrac12,\qquad
 \sum_k(k-\alpha_H)B_{H,k}<\alpha_H,\qquad
 B_{H,1}+D(\beta_H)\sum_{k\ge2}(k+\beta_H)B_{H,k}<1              \tag{1}
\]
for \(m\ge13\).
The lower and upper argument pieces give the corrections
\(\phi_m\) at the actual root radii. Here are the signs and
radii needed for that comparison.

On \(|z|=R_m\), compare \(L_m(z)\) with \(U_m(1/\bar z)\)
on the same ray. The coefficient difference at each index
is bounded by the corresponding right cap. Thus the
difference of moduli is bounded by \(\sum B_{U,k}\), and
both half moduli are at least \(1-\sum B_{U,k}\).
The ratio is less than \(2<R_m^m\). Both halves are
nonzero on the disk by the same bounds. Rouché in the
two orientations gives \(R_m^{-1}<|z|<R_m\) for every root.

For \(0<\theta\le\pi\), \(0\le r\le R_m\), the two phase
bounds in (1) imply
\[
                 -\beta_H\theta<\arg H(r e^{i\theta})
                                           <\alpha_H\theta.
\]
For the upper bound use decreasing sinc on \([0,\pi]\)
and \(|\sin x|/x\le1/\pi\) past \(\pi\). Since \(\alpha_H\le3/5\),
\(\sin(\alpha_H\theta)/(\alpha_H\theta)>1/\pi\), by the
cubic lower bound for sine and \(25/8<\pi<22/7\).
The imaginary part after rotation by \(-\alpha_H\theta\)
is therefore strictly negative by the second inequality of (1).

For the lower bound rotate by \(\beta_H\theta\). The first term
can contribute at most \(\sin(\beta_H\theta)\) times its
coefficient in the negative direction. For \(k\ge2\), use
\(-\sin x\le2x/9\) when \(x\ge\pi\). This is immediate for
\(x\ge9/2\). For \(\pi\le x\le9/2\) put \(y=3\pi/2-x\);
then \(-\sin x=\cos y\le1-3y^2/8\), and
\[
 \frac{2x}{9}-1+\frac38y^2
         =\frac38(y-8/27)^2+\frac{\pi}3-1-\frac8{243}>0.
\]
The cubic lower bound on \(\sin(\beta_H\theta)\) gives the
factor \(D(\beta_H)(k+\beta_H)\) after division. The sum
in (1) now makes the rotated imaginary part strictly positive.
The rotation signs use the principal arguments because the
unrotated halves satisfy \(|H-1|<1/2\) and both rotated
arguments lie in \((-\pi,\pi)\).

At a root of angle \(\theta\) the correction to \(m\theta\)
is the sum of the arguments of the halves on that ray at
reciprocal radii. For each row of the table the width of its
allowed interval is strictly below \(\theta\), because the
sum of the four constants is less than one. Subtraction at
adjacent roots of common argument would give a number strictly
between zero and \(2\theta\) equal to a multiple of \(2\pi\).
Positive coefficients exclude the positive real ray, and the
lower half-plane follows by conjugation. This proves
the analytic theorem for \(n\ge13\) and \(N=n+1\).

## Exact scalar replay

~~~python
from fractions import Fraction as F
from math import factorial, prod

def rising(a, k):
    return prod(a+j for j in range(k))

def exp_lo(x):
    return sum(x**j/factorial(j) for j in range(21))

def f_upper(a, m, j):
    den = 10**8
    target = rising(a, j)**m * den**m
    weight = rising(a, m)**j
    lo, hi = -1, den
    while hi-lo > 1:
        mid = (hi+lo)//2
        if mid**m*weight >= target:
            hi = mid
        else:
            lo = mid
    assert hi**m*weight >= target
    return F(hi, den)

def cap(a,s,m,k,j):
    loss = s*k*(1-F(k,m))/2-F(7*k,10*m)
    assert loss > 0
    return f_upper(a,m,j)/exp_lo(loss)

def check(head, a, beta, t0=F(0), t1=F(0)):
    A = sum(head)+t0
    B = sum(k*v for k,v in enumerate(head,1))+t1
    D = F(2,9)/(beta*(1-(F(22,7)*beta)**2/6))
    assert A < F(1,2)
    assert B-a*A < a
    assert head[0]+D*(B-head[0]+beta*(A-head[0])) < 1

assert 1+F(7,10)+F(7,10)**2/2+F(7,10)**3/6 > 2
assert 1-(F(22,7)*F(3,5))**2/6 > F(8,25)
assert F(25,8)/3-1-F(8,243) > 0
for a,s,angles,betas in [
    (F(1,2),F(13,20),(F(203,1000),F(559,1000)),(F(38,1000),F(194,1000))),
    (F(1),F(83,80),(F(254,1000),F(518,1000)),(F(51,1000),F(166,1000))),
]:
    assert sum(angles)+sum(betas) < 1
    assert F(7,13)/exp_lo(s/4) < F(1,2)
    assert cap(a,s,20,10,10) < F(1,2)**10
    for m in range(13,20):
        for side,(alpha,beta) in enumerate(zip(angles,betas)):
            head = [
                cap(a,s,m,k,m-k if side else k)/(2 if 2*k == m else 1)
                for k in range(1,m//2+1)
            ]
            check(head,alpha,beta)
    q = F(1,2)
    t0 = q**10/(1-q)
    t1 = t0*(10+q/(1-q))
    for side,(alpha,beta) in enumerate(zip(angles,betas)):
        head = [cap(a,s,20,k,20-k if side else k) for k in range(1,10)]
        check(head,alpha,beta,t0,t1)
print("All private totals one/two signed half-phase bounds passed.")
~~~
