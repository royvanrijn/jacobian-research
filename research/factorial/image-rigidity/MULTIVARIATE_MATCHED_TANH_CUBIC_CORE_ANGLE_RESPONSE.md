# A third-order matched tanh test for simple core rays

Fix an ordinary two-monomial core and write \(Q_m\) for its order-\(m\)
moment polynomial. Append a nonempty collection of unequal shared
ordinary rows. For each appended row let \(d\) be its absolute gap,
\(\epsilon\) its orientation, and \(v=d/(a+b)\) its rate. The rate
is recovered from its exact matched Gauss pair by
\(v=(n/T_{n+1}-(n-1)/T_n)^{-1}\).
The [quadratic angular response](MULTIVARIATE_MATCHED_TANH_SIMPLE_CORE_ANGULAR_RESPONSE.md)
gives an all-degree weak-tail criterion when its leading signed
test is strict. The following test also handles some contacts
where the quadratic response reaches zero.

The tail statistics below are sums over rows:
\[
 W=2\sum d v,\quad P=\sum\epsilon d v^2,\quad E=4\sum v^2,
 \qquad U_3=\sum d v^3,\quad V_3=\sum\epsilon v^3,\quad R_3=\sum v^3/d.
\]
Only \(W\) is first-order. The next two statistics are
\(O(W^2)\), and the last three are \(O(W^3)\), uniformly
in the number and size of the rows.

## The cubic response

At degree \(m\) put \(\xi=k/m-1/2\), \(q=k(m-k)\) and define
the diagonal coefficient multipliers
\[
\begin{aligned}
 A(k)&=-q/(2m),&
 B(k)&=2q\xi/(3m),&
 C(k)&=q/(4m^2),\\
 U(k)&=-2q(\xi^2+1/4)/(3m),&
 V(k)&=-4q\xi/(3m^2),&
 R(k)&=-2q/(3m^3).
\end{aligned}
\]
After positive tail normalization, the multiplier of the \(k\)-th
core coefficient has logarithm
\[
 W A(k)+P B(k)+E C(k)+U_3 U(k)+V_3 V(k)+R_3 R(k)+O_m(W^4).
                                                               \tag{1}
\]
To obtain (1), let \(s=k/m\), \(M=(a+b)/2\), and
\(x=M(1+2\epsilon v(s-1/2))\) on a row. Its contribution to the
second derivative of the per-order chord is
\(d^2 m\psi_1(1+mx)\). The trigamma integral gives
\[
          m\psi_1(1+mx)=
        x^{-1}-(2m x^2)^{-1}+(6m^2x^3)^{-1}+O_m(x^{-5}).
\]
For example this follows from
\(0\le 1+z^2/3-z\coth z\le z^4/45\) for \(z>0\).
Multiply the bounds by \(\sinh z\) and compare power series.
For the upper bound on the difference, the coefficient of
\(z^{j+1}\) in the difference of the two sides, multiplied by
\((j+1)!\), is \(j(j-2)(j-4)(j+4)/45\) for even \(j\ge0\).
The lower bound has coefficients \(j(j-2)/3\) in the same
normalization. Substituting \(z=y/2\) bounds the remainder
in the Gamma kernel by \(y^4/720\).

Since \(1/M=2v/d\le2v\), expand the first reciprocal through
degree two in \(2\epsilon v(s-1/2)\), the second through
degree one, and the third at degree zero. The error on each row
is \(O_m((d^2/M)v^3)\), uniformly on \(0\le s\le1\).
It sums to \(O_m(W^4)\) since \(v\le W/2\).
The retained second derivative is
\[
 W-4P(s-1/2)-E/(2m)+8U_3(s-1/2)^2
              +(8/m)V_3(s-1/2)+4R_3/(3m^2).
\]
Integrating with zero endpoint values and multiplying by \(m\)
gives (1).

Let \(\zeta\) be a simple core root, write \(Q_m(z)=\sum a_k z^k\),
and use the evaluation
\[
       \langle f\rangle=
       \frac{\sum_k a_k\zeta^k f(k)}{\sum_k k a_k\zeta^k}.
\]
Put \(\alpha=-\langle A\rangle\), \(p=-\langle B\rangle\),
\(q_0=-\langle C\rangle\), and \(F(k)=k\alpha+A(k)\).
The next coefficients of the root response are
\[
\begin{aligned}
 \beta&=-\langle F^2/2\rangle,&
 \gamma&=-\langle F\,k\beta+F^3/6\rangle,\\
 r_1&=-\langle F(kp+B)\rangle,&
 r_2&=-\langle F(kq_0+C)\rangle,\\
 g_1&=-\langle U\rangle,&g_2&=-\langle V\rangle,&
 g_3&=-\langle R\rangle .
\end{aligned}
\]
The logarithmic displacement of the nearby perturbed root is
\[
 W\alpha+W^2\beta+Pp+Eq_0+
     W^3\gamma+WP r_1+WE r_2+U_3g_1+V_3g_2+R_3g_3+O(W^4).
                                                               \tag{2}
\]
Indeed, evaluate the perturbed polynomial at \(\zeta e^w\)
using (1) and compare each power of the assigned orders of the
statistics. The coefficient of \(w\) at the original parameters is
the nonzero denominator of \(\langle\cdot\rangle\). A fixed
Rouché disk gives the unique nearby perturbed root and bounds
its displacement uniformly. Substitution of the displayed
response leaves residual \(O(W^4)\); integration of the
derivative on the short displacement segment converts residual
to root error, because that derivative stays near its nonzero
original value. Positive scaling does not affect the imaginary
part of (2).

## A sufficient all-degree criterion

Take an integer \(J\) at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For every core ray coincidence below \(J\) assume both roots are
simple. Write \(\Delta\) for upper minus lower in (2).
If \(\Im\Delta\alpha\ne0\), the first-order test applies.
Otherwise put
\[
\begin{aligned}
 H(u)&=\Im\left(\Delta\beta+\tfrac u4\Delta p+u^2\Delta q_0\right),\\
 L(u)&=\Im\left(\tfrac u4\Delta r_1+u^2\Delta r_2\right),\\
 M(u)&=\Im\left(\tfrac18(u^2\Delta g_1+
                                  u^3\Delta g_2+u^4\Delta g_3)\right),\\
 c_*&=\Im\Delta\gamma,\qquad
 \mathcal U=\{0\}\cup\{\pm1/d:d\in\mathbb Z_{>0}\}.
\end{aligned}
\]
One may reverse the sign of all these quantities in an individual
comparison. Assume \(H\ge0\) on \(\mathcal U\). The following
conditions each suffice in the specified case:

- If \(H(0)>0\), require \(c_*+L(u)+M(u)>0\) at every zero of
  \(H\) on \(\mathcal U\).
- If \(H(0)=0\) and \(H\) is not identically zero, require \(c_*>0\).
- If \(H\) is identically zero, require
  \(c_*+t L(u)+t^2 M(u)>0\) for \(0\le t\le1\), \(u\in\mathcal U\).

If each comparison passes, every sufficiently weak appended
ordinary tail gives disjoint root arguments in every degree.
The tolerance is computable by algebraic root isolation and
the uniform coefficient and root errors.

For the proof, put \(\lambda_i=2d_iv_i/W\) and \(u_i=\epsilon_i/d_i\)
on each appended row. The \(\lambda_i\) sum to one. The exact
normalized quadratic coefficient is
\[
       \sum_i\lambda_i^2 H(u_i)+(1-\sum_i\lambda_i^2)H(0)\ge0,
\]
and the normalized cubic coefficient is
\(c_*+\sum_i[\lambda_i^2L(u_i)+\lambda_i^3M(u_i)]\).
If \(H(0)>0\), quadratic coefficient tending to zero forces
\(\max_i\lambda_i\to1\); that row's atom must eventually lie
among the finitely many zeros of \(H\), since \(H(0)>0\).
The cubic coefficient then has a positive lower bound by the
first test. Quantitatively a quadratic coefficient at most
\(\eta\) forces \(\max_i\lambda_i\ge1-\eta/H(0)\).
If \(H(0)=0\), nonnegativity at both signs as \(u\to0\) forces
the linear term of \(H\) to vanish. In the nonzero case
\(H(u)=h u^2\) with \(h>0\). If the quadratic coefficient is small,
the sum of \(\lambda_i^2\) on atoms bounded away from zero
is small, so the cubic coefficient is close to \(c_*\).
In the identically zero case the cubic coefficient is the
\(\lambda_i\)-weighted average of
\(c_*+\lambda_i L(u_i)+\lambda_i^2M(u_i)\), uniformly positive
by compactness. These bounds and (2) prove the strict sign in all
three cases. Distinct core rays stay separate by continuity,
and the absolute cutoff handles the remaining degrees.

The first two boundary tests are sharp for a uniform lower
bound at order \(W^3\) in their positive quadratic orientation.
For a zero atom with \(H(0)>0\) and cubic value zero, a pure
ordinary row of that type gives \(o(W^3)\) splitting.
If its cubic value \(K<0\), append also a unit-gap row whose
weight \(2v\) is asymptotic to
\((-K/(2H(0))) W_{\rm pure}^2\), by rounding its smaller exponent
to an integer. Its normalized share is asymptotic to
\((-K/(2H(0)))W\). The resulting positive quadratic coefficient
cancels \(K\) at order \(W^3\).
For \(H(u)=h u^2\) and \(c_*<0\), take one ordinary row of increasing
gap \(d\) and smaller exponent asymptotic to
\((-c_*/h)d^4\). Then \(H(1/d)/W\to-c_*\).
If \(c_*=0\), use smaller exponent asymptotic to \(d^3\).
These constructions concern angular margins and do not imply
an exact coincidence at nonzero tail strength.

## Exact algebraic comparison that needs the cubic test

For an illustration independent of the factorial core
restriction, put \(x=(z+z^{-1})/2\) and use the polynomials
\[
 \begin{aligned}
 Q_{12}(z)&=z^6(x+12x^2+x^5+2000x^6),\\
 Q_{13}(z)&=z^6(1+z)(x+(25/2)x^2+(595/24)x^4+x^5+x^6).
 \end{aligned}
\]
The powers \(x^5,x^6\) supply strictly positive coefficients
in every degree after expansion. The point \(z=i\) is a
simple root of each, since \(x=0\), \(x'(i)=1\), and the
factor in parentheses has linear coefficient one.

At this pair of roots the common first angular coefficient
is \(-1\). The unit-gap second angular coefficient is \(12\)
in both degrees. For reciprocal cores the response \(p\) is
real at a simple unit-circle root, since \(BQ_m\) is
anti-reciprocal. Also \(q_0=-\alpha/(2m)\).
It follows that \(H(u)=(1-u^2)/312\), so only the
atoms \(\pm1\) require the cubic test. The cubic difference
at the positive unit atom is \(17158441/2530944>0\).
Reversing the tail orientation takes each perturbed root
to its reciprocal conjugate because the core is reciprocal.
It leaves its angle unchanged, so the negative unit atom
gives the same value. Thus the comparison at this ray
uses the cubic test where the strict quadratic test does
not apply.

For an exact replay take an increasing unit-gap tail. Then
\(\ell=1/W-1/2\), and the normalized logarithmic coefficient
multiplier is exactly
\[
 \sum_{h=1}^k\log(1+W(h/m-1/2))
   -\frac{k}{m}\sum_{h=1}^m\log(1+W(h/m-1/2)).
\]
The following symbolic check expands these logarithms and
uses the evaluation formula above.

~~~python
import sympy as s

z = s.symbols("z")
x = (z+1/z)/2
R = s.Rational
data = {12: z**6*(x+12*x**2+x**5+2000*x**6),
        13: z**6*(1+z)*(x+R(25, 2)*x**2+R(595, 24)*x**4+x**5+x**6)}
responses = {}
for m, expression in data.items():
    polynomial = s.Poly(s.expand(expression), z)
    assert all(polynomial.nth(k) > 0 for k in range(m+1))
    assert polynomial.eval(s.I) == 0 and polynomial.diff().eval(s.I) != 0
    wt = [polynomial.nth(k)*s.I**k for k in range(m+1)]
    denominator = sum(k*w for k, w in enumerate(wt))

    def ev(values):
        return s.cancel(sum(w*v for w, v in zip(wt, values))/denominator)

    def power(k, j):
        return sum((R(h, m)-R(1, 2))**j for h in range(1, k+1))

    def log_coefficient(k, j):
        return R((-1)**(j+1), j)*(power(k, j)-R(k, m)*power(m, j))

    l1, l2, l3 = [[log_coefficient(k, j) for k in range(m+1)]
                  for j in (1, 2, 3)]
    alpha = -ev(l1)
    f = [k*alpha+l1[k] for k in range(m+1)]
    beta = -ev([l2[k] + f[k]**2/2 for k in range(m+1)])
    gamma = -ev([l3[k] + f[k]*(k*beta+l2[k]) + f[k]**3/6
                 for k in range(m+1)])
    responses[m] = [s.cancel(s.im(s.expand_complex(t)))
                    for t in (alpha, beta, gamma)]

assert responses[12][:2] == responses[13][:2] == [-1, 12]
assert responses[13][2] - responses[12][2] == R(17158441, 2530944)
~~~
