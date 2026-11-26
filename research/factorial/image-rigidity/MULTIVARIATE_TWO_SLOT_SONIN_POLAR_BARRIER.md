# A root-relative Sonin barrier for the two-slot polar comparison

For a circular convolution of two Gegenbauer factors
of shapes \(\lambda_1,\lambda_2\) in degree \(n\), assume
\(\lambda_1\lambda_2\ge n^2/4\).
The exact two-slot differential equation gives quantitative
bounds for the force at each positive-half root.
The bounds apply to the native polar comparison and can
be written entirely in the recovered matched tanh slots.
For the ordinary doubled shared unit pair
\[
                  (1,2),(2,1),(1,2),(2,1),
\]
they prove the strict lower-root sign of the actual upper
polynomial at every lower root with
\(\tan(\theta/2)\ge1/5\) in \(0<\theta<\pi\).
The small-angle comparison for this family remains open.

## Two upper bounds on the force

Let \(f\) be the centered reciprocal polynomial,
\[
 K=n/2,\quad A_j=K+\lambda_j,\quad B=K+\lambda_1+\lambda_2,\quad
 C=K^2-\lambda_1\lambda_2,\quad L=K A_1 A_2 .
\]
Set \(\Omega^2=L/B\), \(D=\Omega^2-C>0\), and
\(t=\tan(\theta/2)\). The
[two-slot equation and force sign](MULTIVARIATE_TWO_SLOT_MIRRORED_FORCE_COMPARISON.md#negative-force-for-a-two-slot-convolution)
give
\[
         f'''+C f'=B t(f''+\Omega^2 f),\qquad
                     r=f''(\theta_i)/f'(\theta_i)<0
\]
at each positive-half root \(\theta_i\).
Write \(w=f''+\Omega^2 f\) and \(g=f'\), so
\(w'=B t w+Dg\). In particular \(w(\pi)=0\) by the
original equation before division by \(\cos(\theta/2)\).
The Sonin functional
\[
          {\cal E}=g^2+\Omega^2 f^2-w^2/D
\]
satisfies \({\cal E}'=-2B t w^2/D\le0\).
At any interior maximum of \(|g|\) of positive value,
\(f''=0\), so \(w=\Omega^2f\) and
\({\cal E}=g^2-(\Omega^2 C/D)f^2\ge g^2\).
Also \({\cal E}(\pi)\ge g(\pi)^2\).
Since \({\cal E}(\theta_i)=g(\theta_i)^2(1-r^2/D)\),
these observations show
\[
               |g(\varphi)|\le |g(\theta_i)|
                             \qquad(\theta_i\le\varphi\le\pi).
\]

With \(p(\theta)=\cos^{2B}(\theta/2)\),
\((pw)'=Dpg\). Integrate to \(\pi\) and use the last
bound. The substitution \(x=\cos(\varphi/2)\) gives
\[
 \frac{1}{p(\theta)}\int_\theta^\pi p(\varphi)\,d\varphi
   <\frac1{(B+1/2)\tan(\theta/2)} .
\]
Indeed \((1-x^2)^{-1/2}\) under the integral is at most
its endpoint value at \(x=\cos(\theta/2)\), with strict
inequality away from that endpoint. At a root, with
\(\rho=-r>0\), we obtain
\[
                  \rho<\frac{D}{(B+1/2)t}.                      \tag{1}
\]

A separate comparison gives the quadratic bound
\[
                            \rho^2+B t\rho<D.                   \tag{2}
\]
Normalize by a real constant so \(g(\theta_i)=1\).
Then \(w(\theta_i)=-\rho\), and \(f>0\) until its
next positive-half root or \(\pi\), whichever comes first.
Put \(a=B\tan(\theta_i/2)\).
While \(w\le0\) on this interval,
\(g'\le w\) and \(w'\le Dg+a w\), since \(t\) increases.
The constant comparison matrix on \((g,w)\) has
nonnegative off-diagonals, so its exponential has
nonnegative entries and respects componentwise
differential inequalities. If
\(\rho\ge\rho_*=(\sqrt{a^2+4D}-a)/2\), the two
functions are bounded above by the exact solution
\[
                  (e^{-\rho_*s},-\rho_* e^{-\rho_*s}),
                         \qquad s=\theta-\theta_i .
\]
The second component is strictly negative. By continuity
the comparison holds through the entire interval and
at its endpoint (by a limit if the endpoint is \(\pi\)).
At a next positive-half root \(g<0\) and \(r<0\), so
\(w=r g>0\). At \(\pi\) instead \(w=0\).
Either contradicts the upper bound on \(w\), proving (2).

## Exact matched tanh input and differential reduction

In one matched Gauss slot the recovered shape is
\(\lambda_n=(n-1)(T^{-1}-1)/2\) and
\(\lambda_N=n((T')^{-1}-1)/2\), by the
[exact tanh coordinates](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid).
In particular
\[
       A_j=\tfrac12+\frac{n-1}{2T_j},\qquad
       n/2+\lambda_{N,j}=\frac{n}{2T'_j},\qquad
       D=\frac{\lambda_1\lambda_2(n+\lambda_1+\lambda_2)}{B}.
\]
Thus (1) and (2) use the actual linked shape data.
The integer Gauss grid also gives all symmetric updates
in the [native Euler factorization](MULTIVARIATE_MIRRORED_POLAR_EULER_FACTORS.md#factorization-and-the-attached-shapes).
A symmetric update for a gap \(d\), at the relevant shape
\((n\ell+h)/d\), acts in centered variables, up to a
positive scalar, as
\(((n\ell+h)/d+n/2)^2+\partial_\theta^2\).
The residual coefficient \(a_j/\beta_j\) in that
factorization is \(n/(2T'_j)\).
This uses the recovered ordinary slots themselves rather
than only inequalities between their aggregate moments.

For a constant-coefficient differential operator, reduce
its action on \(f\) to \((a,b,c)\), denoting
\(a(t)f+b(t)f'+c(t)f''\). One more derivative applies
\[
 (a,b,c)\longmapsto
       ({\cal T}a+Lt c,\ a+{\cal T}b-C c,\
                      b+{\cal T}c+B t c),
              \qquad {\cal T}=\tfrac12(1+t^2)\partial_t .
                                                               \tag{3}
\]
This is exact away from \(\theta=\pi\).
At a lower root only \(b+c r\) remains, which can be
bounded by (1) or (2) when \(c>0\).
The structural central root has distinct parity at
consecutive orders.

## An angular range for two equal shared unit pairs

For the displayed ordinary profile the two lower shapes
are \(n+1\); put \(A=3n/2+1\).
The two symmetric updates act, after positive
normalization, as \((A^2+\partial_\theta^2)^2\).
The upper residuals then give the centered evaluation
at a lower root, divided by \(f'\), with the sign of
\[
     \frac{((A+1)^2-\partial_\theta^2
                     +2(A+1)t\partial_\theta)
                    (A^2+\partial_\theta^2)^2 f}{f'} .
\]
Indeed the centered residual factors are proportional to
\(A+1+i\partial_\theta\), and the upper lift takes
\(\operatorname{Re}(e^{-i\theta/2}\cdot)\), with positive
\(\cos(\theta/2)\) on this half-circle.
The reduction (3) gives the value \(P(t)+Q(t)r\), where
\[
 \begin{split}
 P(t)&=\frac{t(n+1)^3(3n+2)}4
                 \left((5n+6)t^2+21n+22\right),\\
 Q(t)&=\frac{(n+1)^2}{16}
       \left(5(5n+4)(5n+6)t^4
         +2(117n^2+234n+112)t^2+5(3n+2)(3n+4)\right).
 \end{split}
\]
Both are positive when \(t>0\).
For (1), the exact difference after clearing denominators is
\[
 B(B+\tfrac12)t P-(L-BC)Q
  =\frac{(n+1)^4(3n+2)}{16}
            \left((291n^2+502n+216)t^2-5(3n+2)(3n+4)\right).
\]
It is positive for \(t\ge1/2\).

For \(1/5\le t\le1/2\) and \(n\ge3\), use (2).
It suffices that \(P^2+B t P Q-D Q^2>0\).
Writing \(x=t^2\), this quantity times
\(256 B/((n+1)^5(3n+2))\) is a degree-four polynomial
in \(x\), with coefficients \(q_0,\ldots,q_4\).
Here \(q_1,q_2>0\), \(q_4<0\), and \(q_3+q_4/4>0\)
for \(n\ge3\), so on \(1/25\le x\le1/4\) it is
bounded below by \(q_0+q_1/25+q_2/625\).
That last expression is positive because its coefficients
in \(n-3\) are positive.
The exact symbolic check below verifies these identities
and inequalities.
For \(n=2\), the positive lower root has \(\theta>\pi/2\)
from the quadratic with positive coefficients, so the
first bound applies. At \(n=1\) only the structural
central lower root occurs. Thus all lower roots in the
stated angular range satisfy the strict upper sign.

## Exact symbolic check

~~~python
import sympy as s

t, n, x, u = s.symbols("t n x u")
K = n/2
A = K+n+1
B = K+2*(n+1)
C = K**2-(n+1)**2
L = K*A**2
D = L/B-C

def diff(v):
    a,b,c = v
    def T(p):
        return (1+t*t)*s.diff(p,t)/2
    return s.Matrix([T(a)+L*t*c,a+T(b)-C*c,b+T(c)+B*t*c]).applyfunc(s.expand)

f = s.Matrix([1,0,0])
f2 = diff(diff(f))
advanced = A**4*f+2*A**2*f2+diff(diff(f2))
actual = (A+1)**2*advanced-diff(diff(advanced))+2*(A+1)*t*diff(advanced)
Pn=(n+1)**3*(3*n+2)*((5*n+6)*x+21*n+22)/4
Qn=(n+1)**2*(5*(5*n+4)*(5*n+6)*x*x
              +2*(117*n*n+234*n+112)*x+5*(3*n+2)*(3*n+4))/16
P=t*Pn.subs(x,t*t)
Q=Qn.subs(x,t*t)
assert s.expand(actual[1]-P)==s.expand(actual[2]-Q)==0
assert s.expand(B*(B+s.Rational(1,2))*t*P-(L-B*C)*Q
                -(n+1)**4*(3*n+2)*((291*n*n+502*n+216)*t*t
                                   -5*(3*n+2)*(3*n+4))/16)==0
assert all(c>0 for c in s.Poly(291*n*n+502*n+216
                              -20*(3*n+2)*(3*n+4),n).all_coeffs())

q=s.Poly(s.cancel(256*B*(x*Pn**2+B*x*Pn*Qn-D*Qn**2)
                         /((n+1)**5*(3*n+2))),x)
for item in [q.nth(1),q.nth(2),-q.nth(4),
             q.nth(3)+q.nth(4)/4,
             q.nth(0)+q.nth(1)/25+q.nth(2)/625]:
    assert all(c>0 for c in s.Poly(item.subs(n,u+3),u).all_coeffs())
~~~
