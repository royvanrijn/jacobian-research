# Factorial mark products for independent opposed heights

Fix an ordinary two-monomial core \(C\), with arbitrary orientations.
Append \(r\ge1\) copies of the increasing row \((L,L+d)\), where
\(L\ge0\) and \(d\ge1\) are integers. Each other moving row is decreasing.
Write these rows as \((Y_i+e_i,Y_i)\), of multiplicity \(s_i\ge1\).
The row list, its heights \(Y_i\ge0\), its gaps \(e_i\ge1\), and its
multiplicities are independent. Equal fixed core rows may be dropped.
Reversing both columns gives the other choice of orientation.

For adjacent orders \(m,m'\ge2\), let
\[
 \begin{split}
 B_C&=m!m'!\!\prod_{(a,b)\in C}\prod_{j\in\{m,m'\}}
             \prod_{u=1}^{j|b-a|}\big(j\min(a,b)+u\big),\\
 B_m^-&=\prod_i\prod_{u=1}^{m e_i}(mY_i+u)^{s_i},\qquad
 K=mm'd,\qquad D_0=\operatorname{lcm}(1,\ldots,K).
 \end{split}
\]
All products over the core count row repetitions. For any integer
\(q\) with \(1\le q\le\min(m,m')\), form
\[
 T_q=\{t:1\le t\le qd,\ m\nmid t\},
 \qquad x_t=mL+t,\qquad a_t=\frac{x_t}{\gcd(x_t,D_0)} .
\]

**Theorem.** If the two original moment polynomials have a rational
common factor of degree at least \(q\), then
\[
                  \left(\prod_{t\in T_q}a_t\right)^r
                 \ \bigm|\ B_m^- D_0^{(m+1)r} B_C^{m+1}.       \tag{1}
\]
In particular \(q=2\) gives an exact coprimality certificate because the
even-order moment is strictly positive for real arguments. There is no
restriction on the prime factorization of the marks or on the relative
sizes of the opposing rows. Only the opposing factorial strips at the
marked order \(m\) occur on the right.

The [paired first-mark argument](MULTIVARIATE_PAIRED_PRIME_MARK_AFFINE_DESCENT.md)
uses the direction of one opposing row to constrain a chosen prime
power. Here every prime power of every mark above the fixed cutoff is
charged once, and the first \(q\) blocks are usable for a degree-\(q\)
test. The [common-slope estimate](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md)
controls arbitrary retained orientations by pricing both orders.
Monotonicity of all the moving opposing rows gives the one-order
charge in (1).

## Proof of the product condition

Strip the smaller-end factorial from each row at each order.
For each prime \(p\) put \(b=v_p(D_0)\). Suppose \(v_p(a_t)=E>0\)
for some \(t\). The prime power \(p^{b+1}\) exceeds \(K\).
At order \(m\), distinct selected-row ramp factors differ by at
most \(md-1\), so the only factor in that row divisible by
\(p^{b+1}\) is \(x_t\). At order \(m'\), for every selected-row
factor \(m'L+u\),
\[
                     m(m'L+u)-m'x_t=mu-m't .
\]
This is nonzero because \(m,m'\) are coprime and \(m\nmid t\).
Its absolute value is at most \(mm'd\). Thus the selected row at
order \(m'\) has no factor divisible by \(p^{b+1}\).

The valuation layers beyond \(b\) therefore contribute exactly
\(rE\) at marked indices \(k\ge\lceil t/d\rceil\) and zero at every
other marked index. They contribute nothing at order \(m'\).
For layers \(1,\ldots,b\) of the selected rows, remove the common
affine slope \(rd\sum_{u=1}^b p^{-u}\) from both ordinate arrays.
Indeed the valuation of an initial factorial strip at layer \(u\)
is a difference of floors with linear part \(dk/p^u\). Subtracting
this part leaves a function of \(k\) with oscillation at most one.
The binomial coefficients and retained core strips add ordinates
in \([0,v_p(B_C)]\). Together these two contributions, after the
common slope is removed, have oscillation at most
\[
                             \delta_p=rb+v_p(B_C)
\]
in each order.

At order \(m'\), every opposing-row ordinate is nonincreasing
because the sampled factorial strip shrinks as \(k\) increases.
The last edge of the shifted lower Newton polygon ends at its
last coefficient index, and its left endpoint is an index of the
ordinate array. Its slope is at most \(\delta_p\), since the
opposing contribution cannot increase and the residual oscillation
is at most \(\delta_p\). All preceding slopes are at most the
last slope. This reasoning uses the entire opposing-row valuations
in the companion order and places no size condition on them.

For any candidate common shifted slope \(\xi\le\delta_p\), compare
the marked ordinate at \(k\ge\lceil t/d\rceil\) with the one at
\(k=0\), after subtraction of \(\xi k\). The difference is at least
\[
                      rE-v_p(B_m^-)-(m+1)\delta_p .
\]
The bound uses the full possible drop \(v_p(B_m^-)\) of the opposing
rows, the residual oscillation, and \(k\xi\le m\delta_p\).
If the bound is positive, every endpoint of a face of slope \(\xi\)
has index less than \(\lceil t/d\rceil\le q\). All edges with
shared slopes then have total width at most \(q-1\).
The Newton polygon product rule makes the degree of any rational
common factor at most this width. This holds after subtracting the
same affine slope from both polygons.
[Filaseta and Rocks, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

Thus a common factor of degree at least \(q\) forces
\[
                 rE\le v_p(B_m^-)+(m+1)(rb+v_p(B_C)).
\]
A given prime is charged for at most one \(t\), since two distinct
marks differ by less than \(K\) and cannot both be divisible
by \(p^{b+1}\). Summing the prime exponents proves (1). If (1)
fails at \(q=2\), a common factor could have degree at most one.
The even polynomial is the integral of an even power of the
nonzero real binomial under independent exponential measures,
so it has no real root. This excludes a common linear factor.
A rational Bézout identity transfers coprimality to every
characteristic-zero field.

## Consequences for arbitrarily separated scales

Put \(M=\max(1,\max_i(Y_i+e_i))\), taking \(M=1\) when the
opposing row list is empty. Let \(D_-=\sum_i s_i e_i\), and
\(h_q=qd-\lfloor qd/m\rfloor\). For \(L\ge1\), taking logarithms
in (1) gives the necessary condition
\[
 h_q\log(mL)
 \ \le\ \frac{mD_-}{r}\log(mM)
      +(h_q+m+1)\log D_0+\frac{m+1}{r}\log B_C .                \tag{2}
\]
Here \(a_t\ge mL/D_0\), and \(B_m^-\le (mM)^{mD_-}\).

Fix \(d\) and the core, and let \(J\ge16\) be an integer at least
the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The [absolute signed-degree bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
implies \(D_-\le rd+\delta\) at any collision, with
\(\delta=146J-73+|\Delta_C|\) and \(\Delta_C\) the core signed gap.
Consequently there are explicit computable \(c>0,\ \beta>0\)
depending only on \(C,d,J\) such that any physical common root in
any degree requires
\[
                                  M\ge c L^\beta .
\]
For example \(q=2\) in (2) allows
\[
              \beta=\min_{2\le m<J}
                         \frac{2d-\lfloor 2d/m\rfloor}{m(d+\delta)}.
\]
To obtain \(c\), solve (2) with \(D_-/r\le d+\delta\),
\(1/r\le1\), and take the minimum of the resulting positive
constants over \(2\le m<J\). The first pair is excluded by even
positivity. At a failure the adjacent orders \(m,m'\ge2\) are
below \(J\), and either order may be marked.

Along any sequence of possible failures with a fixed marked order
and \(r\to\infty\), (2) yields the sharper asymptotic restriction
\[
              \liminf_{L\to\infty}\frac{\log M}{\log L}
                     \ \ge\ \frac{2d-\lfloor 2d/m\rfloor}{md}.
\]
This follows by dividing by \(\log L\) and using
\(D_-/r\le d+\delta/r\). For precisely two opposed types with fixed
gaps, exchange columns and apply the same argument to the other
type. Both heights are then bounded by explicit powers of the
other at any failure, uniformly in both counts. These conclusions
also apply when one type's curvature contribution vanishes while
the total tail curvature stays positive. The necessary condition
allows distinct power scales and leaves the comparable-scale case
unresolved.

