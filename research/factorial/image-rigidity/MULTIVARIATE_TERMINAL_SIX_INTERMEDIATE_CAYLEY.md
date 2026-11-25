# Shared terminal-six separation on the intermediate outer-mean interval

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are shared. Let \(u\) be the order-five outer geometric mean after
endpoint normalization and removal of binomial coefficients. If
\[
                        \frac{18}{25}\le u\le\frac9{10},
\]
the order-five and order-six polynomials have disjoint root arguments
under independent positive rescalings. The dimension and the exponents
are arbitrary. The estimate joins the upper endpoint of the
[small-reciprocal bands](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md)
to the lower endpoint of the
[near-flat extension](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT_EXTENSION.md).

The [exact checker](verify_multivariate_terminal_six_intermediate_cayley.py)
runs with Python 3 from the repository root. Its interval endpoints
are integers on the lattice \(2^{-48}\mathbb Z\), rounded outwards
at every arithmetic operation. The checker uses neither numerical
roots nor optimization in its replay.

## Bounds from the ordinary rows

In a shared row put \(\ell=\min(a_i,b_i)>0\), \(d=|a_i-b_i|\), and
\(x_h=d/(5\ell+h)\), \(x'_h=d/(6\ell+h)\) for \(1\le h\le d\).
The row contributes at most
\[
 \left(\frac{5\ell+d}{5\ell+5d}\right)^{d/2}
\]
to \(u\), since each squared outer-mean factor is
\((5\ell+h)/(5\ell+4d+h)\). When \(d\ge2\) this gives
\(d/(5\ell)\le(1-u)/(5u-1)\). For \(d=1\) use the exact row mean
\(\sqrt{(5\ell+1)/(5\ell+5)}\). In the absence of a row
\((\ell,d)=(1,1)\), these comparisons give

| Range of \(u\) | Upper cap \(t\) for the order-five reciprocals |
| --- | ---: |
| \([.857,.9]\) | \(1/16\) |
| \([.75,.857]\) | \(1/11\) |
| \([.72,.75]\) | \(27/250\) |

All decimals denote exact rationals. In the first band, \(d=1\)
requires \(\ell\ge3\), because \(.857^2>11/15\).

For any cap \(x\le t\), we have \(x'\le5t/(6-t)\). More precisely,
\[
 \frac{25}{24}\le\frac{5x'}{4x}
 = \frac{25}{4(6-xh/d)}\le\frac{25}{4(6-x)}.
\]
Let \(U=-\log u\), \(W=-\log w\) and \(\rho=W/U\), where \(w\) is
the order-six outer mean. Since the upper endpoint in the preceding
display gives \(1+5x'\le(1+4x)/(1-x/6)\), taking logarithms gives
\[
 \frac{\log(1+(25/6)t)}{\log(1+4t)}
 \ \le\ \rho\ \le\
 1-\frac{\log(1-t/6)}{\log(1+4t)}.                       \tag{1}
\]
For the lower bound the quotient of the concave function
\(z\mapsto\log(1+(25/24)(e^z-1))\) by \(z\) decreases.
For the upper bound
\(-\log(1-x/6)/\log(1+4x)\) increases because its numerator is convex,
its denominator concave, and both vanish at zero. Both comparisons
sum across factors. The interval (1) contains \(25/24\).

We use the
[linked arrays](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
At either order let \(x\) denote that order's reciprocal,
\(t_j\) its cap, \(m=j-1\), and \(s_t=\tfrac12\log(1+mt_j)\).
The sum \(T=\sum x\) is at most each of
\((e^{2U_j}-1)/m\) and \(t_jU_j/s_t\), where \(U_5=U\), \(U_6=W\).
The first uses the product for the outer mean. The second uses
the decreasing quotient \(\log(1+mx)/x\). The linked skew and
defect bounds are
\[
\begin{array}{c|ccc}
  &\text{outer skew}&\text{linked skew}&\text{log defect}\\ \hline
 j=5&|\alpha|\le t_5T&|\alpha-2\chi|\le (3/5)t_5^3 T
   &0\le V=\log(v/u^{3/2})\le s_t^2 U/6\\
 j=6&|\beta|\le (5/3)t_6T&|4\beta-5\eta|\le 6t_6^3 T
   &0\le H=\log(y/w^{8/5})\le s_t^2 W/5 .
\end{array}
\]
Also \(0\le J=-\log z+\tfrac32\log y-\tfrac35\log w
\le(6/5)t_6^4 T\). These follow from the single-factor power bounds
in the [near-flat proof](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT.md),
summing \(x^k\le t_j^{k-1}x\); for the mean defects use
\(G_{1/2}(s)\le s^3/6\), \(G_{3/5}(s)\le s^3/5\), and
\(s\le s_t\).

## Joint Cayley projections

Put \(p=1-u\), \(q=1-w\), \(\lambda=\sqrt p\), and use the
[exact Cayley expansion](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT.md#stable-cayley-coefficients)
under \(t=(\lambda X-1)/(\lambda X+1)\). The symmetric terms are
centered at \(v_0=u^{3/2}\), \(y_0=w^{8/5}\), \(z_0=w^{9/5}\).
To bound the even parts without losing their cancellations, use
\[
\begin{aligned}
 a_5&=u+u(\cosh\alpha-1),&
 b_5&=v_0+(v-v_0)+v(\cosh\chi-1),\\
 a_6&=w+w(\cosh\beta-1),&
 b_6&=y_0+D+y(\cosh\eta-1),\\
 z&=z_0+\tfrac32w^{1/5}D+E,& D&=y-y_0 .
\end{aligned}
\]
Here \(D\le y_0H/(1-H)\),
\(-w^{-3/5}J\le E\le (3/8)w^{-7/5}D^2\), and the order-five
analogue bounds \(v-v_0\). For \(E\), expand
\(w^{-3/5}y^{3/2}\) at \(y_0\), then apply \(e^{-J}\).
The lower error bound uses \(y\le1\). The checker confirms
\(w>7/10\) on the enclosures.

For the opposite-parity terms retain one shared variable in each
order. In the notation \(c,d\) of the Cayley expansion, use
\[
\begin{array}{c|cc}
   &c&d\\ \hline
5&2u\sinh\chi+\epsilon_5&v\sinh\chi\\
6&(5/4)w\sinh\eta+\epsilon_6&y\sinh\eta .
\end{array}
\]
On \(|s|\le1/25\), the hyperbolic series give
\(\cosh s\le1.001\), \(|\sinh s|\le1.001|s|\),
\(\cosh s-1\le .501s^2\), and
\(|\sinh s-s|\le |s|^3/5\). The checker confirms the small-argument
bounds before using these estimates. They give
\[
\begin{aligned}
 |\epsilon_5|&\le
   u(1.001|\alpha-2\chi|+1.01|\chi|^3),\\
 |\epsilon_6|&\le
   w(1.001|4\beta-5\eta|/4+|\eta|^3).
\end{aligned}
\]
For the first cubic error use
\(\sinh(2s)-2\sinh s=2\sinh s(\cosh s-1)\).
For the second, applying the cubic error bound to both terms gives
coefficient \(((5/4)^3+5/4)/5<1\).

The checker substitutes these expressions in the Cayley coefficients
and groups every perturbation into a linear contribution. The same
\(D\) multiplies its contributions from both \(y\) and \(z\), and the
same skew multiplies its contributions at all opposite-parity
degrees. The nominal terms use the binomial series for
\((1-p)^a\), \(1<a<2\). Starting with the positive coefficient
\(a(a-1)(2-a)/6\) at degree three, successive tail coefficients
multiply by \((k-a)/(k+1)\). The checker bounds the remainder by
a geometric series and applies the same procedure to \(w=u^\rho\)
and \(w^{1/5}\). Bounds on logs use the series for
\(2\operatorname{artanh}((x-1)/(x+1))\) with an explicit geometric
remainder.

## The marked row \((\ell,d)=(1,1)\)

At most one such row occurs in the stated band, since its outer-mean
factor is \(\sqrt{3/5}\) and two would give \(u\le3/5\). Orient it
from smaller to larger exponent by interchanging the columns if
necessary. For the residual array define \(u_R=u/\sqrt{3/5}\);
then \(u_R>.929\), so \(p_R=1-u_R<.071\).
A residual row with \(d=1\) has \(\ell\ge5\), because the row
mean would be at most \(\sqrt{21/25}\) otherwise. For \(d\ge2\)
the preceding row comparison applies with \(u_R\) in place of \(u\).
Thus every residual order-five reciprocal is at most \(1/26\).
Apply (1) and the linked skew and defect bounds to the residual
array.

Write the normalized marked-row coefficients as
\[
 A_{j,k}=\binom jk\frac{(j+k)!}
       {(j!)^{\,1-k/j}((2j)!)^{k/j}}.
\]
The checker obtains each by enclosing the \(j\)-th root of a rational
number. Multiply \(A_{j,k}\) by the appropriate symmetric Gaussian
power of \(u_R\) or \(u_R^\rho\) and denote the resulting
coefficient by \(B_{j,k}\). The actual coefficient is
\(B_{j,k}\exp(\sigma_{j,k})\), where the log correction arrays are
\[
\begin{aligned}
 \sigma_5&=(0,\alpha_R,V_R+\chi_R,V_R-\chi_R,-\alpha_R,0),\\
 \sigma_6&=(0,\beta_R,H_R+\eta_R,\tfrac32 H_R-J_R,
                          H_R-\eta_R,-\beta_R,0).
\end{aligned}
\]
Expand \(\exp\sigma=1+\sigma+e\), with
\(0\le e\le \sigma^2/[2(1-|\sigma|)]\). In the linear term use
\(\alpha_R=2\chi_R+(\alpha_R-2\chi_R)\) and
\(\beta_R=(5/4)\eta_R+(4\beta_R-5\eta_R)/4\).
The checker transforms each shared linear contribution before
multiplying by its parameter interval. It bounds the Gaussian
part by the binomial series grouped by powers of \(p_R\).
This removes the large cancellations that arise if each of the
seven sextic coefficients receives an independent skew interval.

## Ray exclusion and sector counts

The symmetric anchors, actual profiles, and the straight segments
between them have normalized interior multipliers at most \(19/20\).
For the first branch this follows from log convexity against the
two outer coefficients and \(u e^{|\alpha|}, w e^{|\beta|}
\le .9/(1-.04)<19/20\). For the marked branch the marked row's
interior multipliers are below \(4/5\), which the checker verifies.
The constant and leading coefficients are one. The inequality
\((19/20)((9/8)^6-1)+(1/8)^5<1\) then proves the root annulus
\(1/8<|t|<8\) by dominance at both ends.

In the first branch, interpolate the four scaled ray parameters
\(b\) in \(p\) using the rational knots in the checker. In the
marked branch, interpolate their original cosines in \(p_R\);
use the same Cayley transform with \(\lambda=1\). The ray parameters
are strictly ordered. The original cosine of each ray is at most
\(-2/5\), so for \(R=|t|\) in the verified annulus,
\[
 |\lambda\Re X|=
 \frac{|1-R^2|}{1-2R\cos\theta+R^2}\le\frac{15}{17}.
\]
It suffices to use \(R\le1\) by inversion; on that interval the
maximum of the bound using \(\cos\theta=-2/5\) is at \(R=1/8\).
Writing \(X=x+iY\), \(\xi=\lambda x\), the ray equation is
\[
 (Y-b)(1+\lambda^2bY)=-b\xi^2,\qquad 0<Y\le b.
\]
Starting with \(b(1-\xi^2)\le Y\le b\), the checker applies this
equation three times to enclose \(Y\). For the first branch only
\(\xi\ge0\) is needed, since reflection \(X\mapsto-\overline X\)
reverses the opposite-parity perturbations and their intervals
are symmetric. The marked branch checks both signs of \(\xi\).

On each interval rectangle or box, the checker excludes a root
by giving some projection \(\Im(P_j(X)\overline X^{\,h})\),
\(0\le h\le j\), a strictly positive or strictly negative enclosure.
It bisects boxes until the projections prove every one. The ordinary
replay has eleven knot intervals. The marked replay has five.

For the first branch, at a given \(p,\rho\), deform actual
coefficients by a straight segment to the nominal symmetric arrays
at \(p,\rho\). All perturbation intervals contain zero and the
certified projections continue to apply. Vary \(\rho\) to \(25/24\)
inside (1), and then \(p\) to \(1/10\) along the symmetric arrays.
At the anchor the checker evaluates \(P_5(iY)/(iY)\) and
\(P_6(iY)\) at the four ray parameters. The respective signs are
\((+,-,-,+)\) and \((+,+,-,-)\), in order of increasing physical
argument. The sextic value at zero is positive. With the leading
signs, these prove that the successive upper sectors contain
sextic, quintic, sextic, quintic and sextic roots. At the marked
anchor, first deform to the marked array times the symmetric
residual by the same straight segment, then vary \(p_R\) to zero
at fixed \(\rho\). The checker proves the same sector assignment
for the marked array by disjoint rational disks. At each disk
center \(c\) and radius \(r\) it verifies
\[
 |P'_j(c)|r>|P_j(c)|+\sum_{k=2}^j|P_j^{(k)}(c)|r^k/k!.
\]
Rouché's theorem gives one root per disk, and the oriented distance
from each center to each ray places its entire disk in the indicated sector.

In both branches the ray exclusions preserve the sector counts.
Count the first and last upper sectors together with their lower
reflections across the real axis. In the last joined sector the
single quintic root is real by conjugacy. The actual sextic is
strictly positive on the real line: its factorial coefficients
are the moments of the two distinct monomials against the
product exponential measure, and the sixth power of their
real linear combination is nonnegative and positive on a set
of positive measure. Its two roots in that joined sector are
therefore nonreal. The root arguments are disjoint.
The original factorial polynomials have rational coefficients,
so their rational Bézout identity gives coprimality over every
characteristic-zero field.
