# Quantitative endpoint attenuation and shared Gauss-factor links

There are exact consecutive-order controls beyond monotonicity
of the normalized multipliers. First, at any fixed distance from
an endpoint their decrement has a uniform curvature lower bound
in every order. Second, at terminal six the two outer means,
the signed outer skews, and both middle shape defects satisfy
linked estimates from the same Gauss factors of an ordinary row.

## An all-order endpoint inequality

Let \(F_{n,k}\) be the endpoint-normalized factorial coefficient
after removing the binomial coefficient. Write
\[
 S=\sum_{a_i\ne b_i}\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 M=\min_{a_i\ne b_i}\max(a_i,b_i).
\]
For integer \(n\ge k\ge1\) the following holds, both from the
left endpoint and from the right endpoint:
\[
 F_{n+1,k}\le F_{n,k}\left(
   \frac{n(n+1+1/M)}{(n+1)(n+1/M)}
 \right)^{Mk^2 S/2}.                                             \tag{1}
\]
In particular, at orders five and six for fully shared ordinary
profiles, \(F_{6,k}\le F_{5,k}(65/66)^{k^2 S}\).

Use the [endpoint integral](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md#bounds-at-a-fixed-distance-from-an-endpoint).
For a row \((a,b)\), write \(f(t)=\psi'(1+t)\). The derivative
of its contribution to \(G_{n,k}=-\log F_{n,k}\) is
\[
 k\int_a^b \big[x f(nx)-a f(y)\big]\,dx,\qquad
 y=(n-k)a+kx.
\]
If \(x\ge a\), then \(nx\ge y\). Since \(t f(t)\) increases,
\(x f(nx)\ge y f(y)/n\), and the integrand is at least
\(k(x-a) f(y)/n\). If \(x\le a\), both inequalities reverse,
with the orientation of the integral also reversed. In
either case \(y\le n\max(a,b)\), and integral comparison gives
\(f(y)\ge1/(1+y)\). Hence the derivative is at least
\[
 \frac{k^2(b-a)^2}{2n(1+n\max(a,b))}.
\]
Sum over the rows and integrate
\(k^2 S/(2n(n+1/M))\) across the consecutive orders to get (1).
For fully shared unequal integer rows one may use \(M=2\)
as a lower bound on the actual minimum. Reversing both
columns gives the estimate at the other endpoint.

## Four linked quantities from each ordinary factor

Write \(x=d/(5a+h)\), \(y=d/(6a+h)\), \(h=1,\ldots,d\),
for a shared ordinary row of lower exponent \(a\) and gap \(d\).
These satisfy the exact relation
\[
               y=\frac{5x}{6-\gamma x},\qquad \gamma=h/d\in(0,1].
 \tag{2}
\]
It is the Gauss shift \(\gamma\) that couples the two orders;
an independent choice of positive Gamma shapes need not
satisfy (2).

Use the normalized shapes and signed skews from the
[terminal-six coefficient formulas](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
Put \(U=-\log u,\ W=-\log w\) and
\[
 V=\log(v/u^{3/2}),\qquad
 Y=\log(y_{\mathrm{mean}}/w^{8/5}),\qquad
 Z=\log(z_{\mathrm{mean}}/(y_{\mathrm{mean}}w^{1/5})).
\]
Here \(y_{\mathrm{mean}},z_{\mathrm{mean}}\) are the order-six
geometric multiplier means, distinct from the reciprocal \(y\)
in (2). Let \(M_j=\sum x^j\) across all Gauss factors. Whenever
every \(x\le2\), the following global inequalities hold:
\[
 \begin{aligned}
 \big|W-\tfrac{25}{24}U\big|&\le\tfrac{25}{144}M_2,\\
 \big|\beta-\tfrac{125}{108}\alpha\big|&\le\tfrac{125}{648}M_3,\\
 \big|Y-\tfrac{125}{108}V\big|&\le\tfrac{125}{432}M_4,\\
 \big|Z-\tfrac{125}{216}V\big|&\le\tfrac{125}{864}M_4 .
 \end{aligned}                                                    \tag{3}
\]
In particular the cap \(x\le2\) includes all fully shared
ordinary profiles with \(S<6\), by the elementary row bounds
in the [hybrid curvature proof](MULTIVARIATE_TERMINAL_SIX_HYBRID_CURVATURE.md#bounds-forced-by-the-integer-rows).

For a single increasing factor, write
\[
\begin{aligned}
 A(x)&=\frac{2\sum_{i=1}^3\log(1+ix)-3\log(1+4x)}{10},\\
 B(y)&=\frac{\sum_{i=1}^4\log(1+iy)-2\log(1+5y)}6 .
\end{aligned}
\]
These are the unsigned outer skews. Their signs in the
normalized arrays agree for the two orders and reverse
together when the row is reversed. Put
\(U(x)=\tfrac12\log(1+4x)\), \(W(y)=\tfrac12\log(1+5y)\).
The local defect derivatives are
\[
 V'(x)=\frac{3x^2}{(1+x)(1+3x)(1+4x)},\qquad
 Y'(y)=\frac{6y^2}{(1+y)(1+4y)(1+5y)},\qquad
 Z'(y)=\frac{3y^2}{(1+2y)(1+3y)(1+5y)} .
\]
They follow by differentiating the single-factor geometric
mean formulas. The other derivatives used in the proof are
\[
 \begin{aligned}
 A'(x)&=\frac{2x(18x^2+20x+5)}
           {5(1+x)(1+2x)(1+3x)(1+4x)},\\
 B'(y)&=\frac{y(5y+2)(24y^2+25y+5)}
           {3(1+y)(1+2y)(1+3y)(1+4y)(1+5y)} .
 \end{aligned}
\]

To check a single-factor bound in (3), denote its difference by
\(E(x,\gamma)=H(y)-c L(x)\). Here \(H\) is respectively
\(W,B,Y,Z\), \(L\) is \(U,A,V,V\), and \(c\) is the
corresponding leading ratio. The function \(H(y)\) increases
strictly, so at fixed \(x\) the two extremes of \(E\) occur at
\(\gamma=0,1\) in (2). For the respective power \(m\) and
constant \(b\) on the right of (3), the
[exact derivative check](verify_multivariate_terminal_six_joint_factor_links.py)
proves
\[
 \frac{d}{dx}(b x^m+E(x,0))\ge0,\qquad
 \frac{d}{dx}(b x^m-E(x,1))\ge0
 \quad(0\le x\le2).
\]
For each rational derivative the check substitutes
\(x=2z/(1+z)\), \(z\ge0\), and clears the denominator.
The resulting numerator and denominator have nonnegative
coefficients, with positive denominator constant.
This uses exact rational polynomial arithmetic.
Both functions vanish at zero, giving the desired local
bounds. Sum them across rows, taking absolute values of
signed skew contributions where orientations differ.
This proves (3).

For example, if \(x\le t\) throughout a profile, the
defect error may be bounded without retaining \(M_4\) as
an independent parameter. The formula for \(V'\) gives
\(V(x)\ge x^3/[(1+t)(1+3t)(1+4t)]\) for \(x\le t\).
Thus
\[
 M_4\le t(1+t)(1+3t)(1+4t)V .
\]
The estimates link the signed and symmetric coefficient
data on both sides of the terminal, leaving the remaining
angle or relative-center check to a root comparison.
