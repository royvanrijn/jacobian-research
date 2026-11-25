# Linked coefficient corrections at a double cubic limiting zero

Use the ordinary profiles and the actual interpolation of the
[moving one-private cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
Here the only private row is \((0,1)\). Within the adjacent pair
the profile is fixed. Assume bounded \(S\), \(\sigma\to1\), and
\(n\to\infty\), with bounded scaled temperatures \(T_n\) and
bounded fields at both orders. Write \(\epsilon=j^{-1/3}\).

**Proposition.** In this setting, a common-root sequence cannot
have both scaled fields tending to the same double zero
\(H_0\) of \(F(T_0,\cdot)\) at their real limiting temperature
\(T_0\). The degree difference \(\Delta\) is arbitrary.
In particular the double-zero obstruction is excluded whenever
the two limiting scaled fields agree, including when
\(\Delta=o(n^{1/3})\).

The proof uses the first two corrections of the actual linked
coefficient integral. They put the local quadratic discriminant
strictly off the real axis on the scale needed to distinguish
its two branches.

## The actual first two corrections

Use \(A_0,\ell,E\) and the normalized interpolation
\(J(T,H,\epsilon)\) from the moving cubic proof. Its leading
function is \(F(T,H)\). Define the real profile quantities
\[
 \alpha=\ell A_0'(0),\quad
 \beta=\ell^4 E(0)
       =\frac{\ell^4}{12}\left(\sum_i d_i^4/a_i^3-1\right),
 \quad
 U=\ell^2 A_0''(0)/2,\quad V=\ell^5 E'(0).
\]
Use \(F_k=\partial_H^k F\). On bounded \(T,H\) sets the actual
interpolation has the uniform expansion
\[
 \begin{split}
 J={}&F+\epsilon(\alpha F_1+\beta F_4)\
    &+\epsilon^2\left(
           U F_2+(V+\alpha\beta)F_5+\tfrac12\beta^2 F_8
                  +\frac1{2\ell}\right)+O(\epsilon^3).
 \end{split}                                                 \tag{1}
\]
The parameters \(T,H\) here are independent in the local
interpolation. The estimate also holds with their derivatives
on smaller complex neighborhoods and with the remainder's
\(\epsilon\)-derivative of order \(O(\epsilon^2)\).
To verify (1), put \(q=\ell\epsilon y\) in the finite integral.
The amplitude from (3) of the moving cubic proof is
\(1+\alpha\epsilon y+U\epsilon^2 y^2+O(\epsilon^3(1+|y|^3))\).
Its complex-order Stirling derivative bounds control the
\(\epsilon\)-derivative of this error. The extra exponent is
\(\epsilon\beta y^4+\epsilon^2 V y^5+O(\epsilon^3|y|^6)\).
Taylor expansion of its exponential gives (1) after multiplying
the amplitude. The differentiated terms are integrable against
\(e^{-c y^3+C y^2+C|y|}\) on the retained interval, by reducing
the fixed physical cutoff. The moving upper endpoint has
strict loss. Finally the cotangent correction divided by
\(\ell\epsilon\) is \(\epsilon^2/(2\ell)+O(\epsilon^4)\),
since \(\xi=H\epsilon^2/\ell\). Its derivative follows
by the same sector estimate. This establishes (1) for the
actual interpolating integral, with the lattice error at
integer orders still exponentially small. The formula for
\(\beta\) follows by integrating the cubic term of
\(\log(1-q)+\sum_i d_i\log(1+d_iq/a_i)\).

Set \(\sigma_i=d_i^2/a_i\) and \(r_i=d_i/a_i\). Then
\(|r_i|\le\sigma_i\) because the nonzero gaps are integers, and
\[
           \sum_i d_i^4/a_i^3=\sum_i\sigma_i r_i^2
                         \le\sum_i\sigma_i^3\le\sigma^3 .
\]
Pass to a subsequence on which \(\beta\) has a limit.
Either its limit is strictly negative, or
\(\sum_i\sigma_i r_i^2\to1\).
In the latter case the maximum row mass tends to one by
\(\sum_i\sigma_i^3\le(\max_i\sigma_i)^2\sigma\).
The other masses tend to zero in sum and the largest row
has \(|r_i|\to1\). Thus \(|d_i|\to1\) and \(a_i\to1\),
so the row is eventually precisely \((1,2)\) by integrality
and positive shared exponents. This is the sole vanishing
case of the first correction.

## The quadratic discriminant

Let \((T_0,H_0)\) be the common limiting double zero and,
by conjugation, take \(\Im H_0>0\). The
[central interval certificate](MULTIVARIATE_HALF_LINE_CUBIC_BOUNDED_TEMPERATURES.md)
in particular ensures \(T_0\ne0\). Let \(t=T-T_0\).
The
[double-zero expansion of the model](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#isolated-exceptional-temperatures)
gives, for the square of the half-difference of its two
nearby zeros,
\[
                  D_0(t)=t+(H_0/6+T_0^2/18)t^2+O(t^3).
 \tag{2}
\]
In each sufficiently small parameter neighborhood, factor
the two local zeros by their sums and products computed
from a fixed field-circle integral of \(J_H/J\).
This gives the exact representation
\(H=H_0+A(T,\epsilon)\pm\sqrt{D(T,\epsilon)}\)
for the two interpolation zeros. All coefficient
expansions and derivative bounds in (1) pass to
\(A,D\) by integration on that fixed circle.
The complementary analytic factor is nonzero and uniformly
bounded away from zero on a smaller neighborhood.
In particular the first variation of \(D\) at \(t=0\) is
minus twice the corresponding variation of the value at
\(H_0\), since the model local factor is
\((H-H_0)^2/2+O((H-H_0)^3)\).

When the limit of \(\beta\) is negative, use
\(F_1=0\) and \(F_4=H_0+T_0^2\) at the double zero.
Equations (1) and (2) give
\[
 D=D_0(t)-2\epsilon\beta(H_0+T_0^2)
                        +O(\epsilon|t|+\epsilon^2).
\]
The imaginary part is bounded below by
\(c(t^2+\epsilon)\) for real \(t\) near zero.

In the remaining case the actual profile contains
the row \((1,2)\). Its other row masses have sum
\(\sigma-1=O(\epsilon_n)\) at the lower order because
the scaled temperature is bounded. Hence along the
whole adjacent order interval
\(\beta=O(\epsilon_n^3)=O(\epsilon^3)\), using
\(\sum_{\rm other}\sigma_i r_i^2\le(\sigma-1)^3\).
Also \(T_0<0\), since \(\sigma\ge1\) and \(T_0\ne0\).
Translate the field by the real amount \(\alpha\epsilon\),
using \(K=H+\alpha\epsilon\).
For the retained row, \(\ell\to1\), \(\alpha\to1\),
\(U\to1/2\), and \(V\to-1/10\).
The latter follows from
\(\log(1+q)+\log(1-q)=\log(1-q^2)\), whose integral
has fifth-order coefficient \(-1/10\).
The other rows contribute vanishingly to these coefficients
because their total mass tends to zero.

Under the translation, the first term \(\epsilon\alpha F_1\)
cancels as a function of \(K\), and the coefficient of
\(\epsilon^2 F_2\) becomes \(U-\alpha^2/2\).
The first term proportional to \(\beta\) is negligible
on the adjacent interval. At the double zero
\(F_5=3-2T_0H_0-T_0^3\), by differentiating the model
differential equation. Thus the second correction in
the translated field, evaluated at the double zero,
tends to
\[
                  \frac15+\frac{T_0 H_0}{5}+\frac{T_0^3}{10}.
]
The imaginary part of minus twice this correction is
strictly positive because \(T_0<0\).
The quadratic factor receives only vanishing contributions
to its second-order coefficient from the first-order term:
that term after translation is proportional to \(\beta\).
Thus the factor at the double zero has the limiting second
variation obtained by dividing the value correction by its
unit \(1/2\). In the translated discriminant the first correction is
\(O(\epsilon\beta)\) and the second has the stated nonzero
limiting imaginary coefficient, with its variation in
\(t\) of order \(O(|t|)\). Therefore \(\Im D\ge c(t^2+\epsilon^2)\)
on the adjacent interval. The real translation does not
change imaginary differences of the two fields.

In both cases, with \(m=1\) in the first and \(m=2\) in
the second, the discriminant and its derivatives satisfy
\[
 |D|\asymp|t|+\epsilon^m,\quad
 \Im D\ge c(t^2+\epsilon^m),\quad
 D_T=1+o(1),\quad D_\epsilon=O(1).
 \tag{3}
\]
The modulus bounds follow from the imaginary lower bound
when \(|t|=O(\epsilon^m)\) and the real linear term in
(2) when \(|t|\) dominates it.
The derivatives of the mean \(A\) are also bounded.

## Adjacent orders stay on one branch

Choose the square root in the first quadrant in (3).
The imaginary separation of the two zeros is
\(2\Im\sqrt D\), bounded below by a positive multiple of
\[
              \frac{t^2+\epsilon^m}{\sqrt{|t|+\epsilon^m}}
                       \ \,\gtrsim\ \epsilon^{3m/4}.
 \tag{4}
\]
The first bound uses \(2\Re\sqrt D\,\Im\sqrt D=\Im D\);
the second follows by splitting at \(|t|=\epsilon^{m/2}\).
The derivative of each branch in \(T\) or \(\epsilon\) is
\(O(\epsilon^{-m/2})\). The fixed profile has
\(T(\epsilon)=\ell^2(1-\sigma)/\epsilon\), so the scaled
branch motion between adjacent orders is
\(O(n^{-1}\epsilon^{-m/2})\). This is smaller than
(4) for both \(m=1,2\).
The lattice error moves each integer-order zero by
an exponentially small distance, since the local
two roots have separation at least a power of \(\epsilon\)
and the analytic unit is nonvanishing.
At a common original root the fields' imaginary parts
after scaling differ only by a factor \((1+1/n)^{2/3}\).
Thus the two actual zeros must use the same local branch.

For that branch the physical field is
\(\xi=\epsilon^2 H(T(\epsilon),\epsilon)/\ell\), where
\(H=K-\alpha\epsilon\) in the translated case.
Its derivative in \(\epsilon\) is
\[
        \frac\epsilon\ell
              (2H-T H_T+\epsilon H_\epsilon).
]
It has magnitude \(O(\epsilon^{1-m/2})\).
The consecutive displacement is hence
\(O(n^{-1}\epsilon^{2-m/2})=o(n^{-1})\), incompatible
with the exact nonzero integer degree drift
\(\Delta\log(1+1/n)\) when \(\Delta\ne0\).
When \(\Delta=0\), the bracket is dominated by its
term \(\mp T_0/(2\sqrt D)\), because \(T_0\ne0\) and
(3) holds. This leading vector, after a fixed sign
and rotation, has a strictly positive real projection
comparable to \(|D|^{-1/2}\) throughout the adjacent
interval, since \(1/\sqrt D\) stays in a quadrant.
The relatively vanishing error cannot cancel this
projection. The resulting consecutive physical
displacement has modulus at least a constant times
\(\epsilon^2/n\), again exceeding the lattice errors.
This excludes the zero degree drift as well and proves
the proposition. \(\square\)
