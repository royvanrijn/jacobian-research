# Row cancellation in matched tanh angular transport

Let \(n\ge13\), \(N=n+1\), and consider the exact moment polynomials
of an ordinary two-monomial profile. On an unequal row \(i\) of
smaller exponent \(\ell_i\) and gap \(d_i\), put
\[
 \rho_i=\frac{d_i}{2\ell_i+d_i},\qquad
 b_{i,h}=\frac{2h-d_i}{d_i},\qquad
 r_{i,h}(j)=\frac{\rho_i}{1+b_{i,h}\rho_i/j}\quad(1\le h\le d_i).
\]
At integer \(j\) these are the normalized Gauss factors
\(jT_{i,h}/(j-1)\). Write \(A=\sum_{i,h}r_{i,h}(n)\),
\(V_p=\sum_{i,h}r_{i,h}(n)^p\),
\(P=\sum_{i,h}\epsilon_i r_{i,h}(n)^2\) with the row orientation
sign \(\epsilon_i\), and \(B_2=\sum_i\rho_i^2\).

**Theorem.** The two moment polynomials have disjoint root arguments
if either of the following conditions holds.
\[
\begin{array}{c|c}
 2A\ge13 & 22|P|+324V_3+4\min(B_2,V_2)/n\le10^{-8}\\[2mm]
 16/3\le2A\le6,\quad n\ge10^7 &
 262154|P|+46000000000V_3+4\min(B_2,V_2)/n\le10^{-18}.
\end{array}
\]
These tests use the same analytic angular lifts as the earlier
[matched tanh criterion](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md)
and
[fixed curvature strip](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md#using-the-midpoint-hierarchy-in-a-paired-phase-lift).
The change is in the order variation of the exact coefficient
kernel, where opposite slots in a single row can be summed
before taking absolute values.

## The row pairing estimate

Let \(L_j(s)\) be the endpoint-normalized log coefficient per
degree, with the binomial coefficient stripped. Write \(g=2s-1\)
and \(\zeta=3/4\). On the disks \(|s|\le5\) and \(|s|\le65538\)
in the two criteria, respectively, \(|g|\le11\) and
\(|g|\le131077\). The bound on \(V_3\) in each criterion gives
\(r_{i,h}(n)<10^{-3}\) and \(r_{i,h}(n)<10^{-9}\),
respectively. Since \(\rho_i\le1\), the terminal slot
\(r_{i,d_i}(n)=\rho_i/(1+\rho_i/n)\) yields
\(\rho_i\le2r_{i,d_i}(n)\). The whole interpolated slot list is
bounded by \(2\rho_i\). Thus \(\Re(1+\epsilon_i rg)\ge\zeta\)
on every such parameter interval and on the disk in use.

Freeze the order in the exact curvature kernel at \(N\). For
each row orientation set
\[
 F(r)=\int_0^\infty e^{-xu}\frac{u/N}{2\sinh(u/(2N))}\,du,
 \qquad x=\frac{1+\epsilon_i rg}{2r}.
\]
This is the [Gauss chord curvature](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#exact-curvature-of-the-assembled-symbol).
Let \(U(r)=r^2 F'(r)\), where the derivative is in \(r\).
Differentiating under the integral gives
\[
 |U(r)|\le\frac{2r^2}{\zeta^2},\qquad
 |U'(r)|\le\frac{4r}{\zeta^3}.
\]
For the second bound use the exact cancellation
\(2F'+rF''=(4r^3)^{-1}\int_0^\infty u^2 e^{-xu}
[u/N\,/\,(2\sinh(u/(2N)))]\,du\).
All bounds are uniform on the stated disk.

The derivative of the row sum of \(F(r_{i,h}(j))\) in real \(j\)
is \(j^{-2}\sum_h b_{i,h}U(r_{i,h}(j))\). Pair \(h\) with
\(d_i-h\), leaving the terminal slot \(h=d_i\) and a possible
zero shift. Its terminal contribution has modulus at most
\(2\rho_i^2/(\zeta^2 j^2)\).
Suppress the row subscript and put \(q=\rho/n\).
For a paired shift \(b>0\), the absolute difference of its two
parameters is \(2b\rho^2/[j(1-b^2\rho^2/j^2)]\).
The sum of \(b^2\) over these pairs is
\((d-1)(d-2)/(6d)\). Bounding \(|U'|\) on the interval of
parameters by \(4\rho/[\zeta^3(1-q)]\) shows that the pairs have
total derivative bounded by
\[
 \frac{4 d\rho^3}{3\zeta^3 j^3(1-q)(1-q^2)} .
\]
Since \(d\rho^3\le(1+q)^3\sum_h r_h(n)^3\) and \(q<1/44\),
the factor \((1+q)^3/[(1-q)(1-q^2)]\) is less than two.
Integrate from \(n\) to \(N\), using
\(\int_n^N j^{-2}dj=1/(nN)\) and
\(\int_n^N j^{-3}dj\le1/(n^2N)\).
The changing-parameter part of \(N(L_N''-L_n'')\) has norm at most
\[
                  \frac{2B_2}{\zeta^2 n}
                         +\frac{8V_3}{3\zeta^3 n^2}.        \tag{1}
\]

## Angular conditions and a strict enlargement

Use the same fixed-parameter grid drift and lower-order bounds as
in the [earlier comparison](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#size-and-adjacent-change-of-the-log-coefficient-correction).
With \(G=\sup|g|\) and \(E_j=j[L_j+A s(1-s)]\), these give,
after adding (1),
\[
 \max\{\|E_n''\|/n,\|E_N''\|/N,\|E_N''-E_n''\|\}
 \le 2G|P|+\frac{4B_2}{n}
  +\left[\frac{2G^2}{\zeta}
       +\frac{2}{3\zeta^3}\left(\frac14+\frac5{12}\right)
       +\frac{8}{3\zeta^3\cdot13^2}\right]V_3 .
\]
Here the lower-order grid error uses \(n^{-2}\le1/4\), the
adjacent grid error uses \(N(n^{-2}-N^{-2})\le5/12\), and
\(2/\zeta^2<4\). The bracket is below \(324\) for \(G=11\)
and below \(46000000000\) for \(G=131077\).
The two cited lifts apply to the resulting exact sampled
coefficient multipliers. Combining this bound with the earlier
\(V_2\) bound gives the minimum in the statement.

For a strict enlargement in the lower strip, take \(n=10^7\),
\(D=497000000000000\), and the mirrored rows
\((\ell,\ell+D),(\ell+D,\ell)\) with \(\ell=9D^2/25\).
Put \(a=D/(2\ell+2D)\), \(b=D/(2\ell)\), and
\(\rho=D/(2\ell+D)\). The \(2D\) slots satisfy \(a\le r\le b\)
and \(P=0\). The following exact rational checks imply that
the new criterion holds and the corresponding \(4V_2/n\)
criterion fails.
~~~python
from fractions import Fraction as Q

n, D = 10**7, 497000000000000
ell = 9*D*D//25
a, b = Q(D, 2*ell+2*D), Q(D, 2*ell)
rho = Q(D, 2*ell+D)
assert 9*D*D == 25*ell
assert 2*D*D < 12*(ell+D)
assert Q(16, 3) <= 4*D*a <= 4*D*b <= 6
assert 46000000000*2*D*b**3 + 8*rho**2/n <= Q(1, 10**18)
assert 46000000000*2*D*a**3 + 8*D*a**2/n > Q(1, 10**18)
zeta = Q(3, 4)
for G, cap in ((11, 324), (131077, 46000000000)):
    coeff = (2*G*G/zeta + 2/(3*zeta**3)*(Q(1, 4)+Q(5, 12))
             + 8/(3*zeta**3*13**2))
    assert coeff < cap
q = Q(1, 44)
assert (1+q)**3 < 2*(1-q)*(1-q*q)
~~~

## The scalar clock and its limitation

Write \(A_j=\sum_{i,h}r_{i,h}(j)\) and
\(S=\sum_i d_i^2/(\ell_i+d_i)\).
For \(S<12\), at every integer \(n\ge13\) there is a strict
linked correction
\[
            \frac{B_2}{5nN}\le A_N-A_n\le\frac{B_2}{nN}.
\]
Indeed on each row, in the derivative of the raw sum take
\(F(r)=r\) and write \(t=2\ell+d,\ q=\rho/j\) for this
calculation. The terminal contribution is
\(\rho^2/[j^2(1+q)^2]\). The paired contributions are
nonpositive. Using the exact difference of the two squared
parameters and the sum of \(b^2\) gives the lower derivative
bound
\[
 \frac{\rho^2}{j^2(1+q)^2}
   \left(1-\frac{2(d-1)(d-2)}{3jt(1-q)^2}\right).
\]
For \(j\ge13\), \(q\le1/j\) and \(d^2/t\le S<12\), so the
subtracted fraction is at most \(8j/(j-1)^2\le13/18\).
Also \((1+q)^{-2}\ge(13/14)^2\). Their product in the lower
bound exceeds \(1/5\), and the upper derivative is bounded by
the terminal term. Integration proves the claim.
Moreover \(r_{i,h}(j)/j=\rho_i/(j+b_{i,h}\rho_i)\) strictly
decreases, so these reference clocks have ratio in \((n/N,1)\).

At \(n=2\) the raw drift can reverse for actual ordinary rows.
The row \((1,7)\) gives \(A_2=20417/4620\) and
\(A_3=17631/4004\), with \(A_3-A_2=-239/15015\).
Adding the reflected row doubles this strict decrease and makes all
signed moments vanish. Its mirrored curvature is \(72/7<12\).
The exact [row concentration comparisons](MULTIVARIATE_MATCHED_TANH_ROW_CONCENTRATION.md)
give further scalar constraints. General angular separation
still requires positive root gaps and a transfer to the
exponentiated factorial coefficients beyond the two stated regimes.

## An explicit all-degree strip from one matched pair

For each matched slot set \(v=(n/T_N-(n-1)/T_n)^{-1}\), which
equals the invariant row rate \(\rho_i\). Define the following sums
over slots, with orientations \(\epsilon_i\):
\[
 W=2\sum v,\qquad B_\infty=-4\sum\epsilon_i v^2,\qquad
 D=8\sum v^3.
\]
They can be evaluated from any one adjacent pair with lower order
at least two. The following bounds imply disjoint consecutive
root arguments in every degree:
\[
 \frac{16}{3}+10^{-12}\le W\le 6-10^{-12},\qquad
                    |B_\infty|\le 10^{-24},\qquad D\le10^{-28}.       \tag{2}
\]
The proof uses the row-cancelled criterion at \(n\ge10^7\)
and the finite-degree quadratic skew estimate at \(n<10^7\).

For the high-order comparison let \(V_{2,\infty}=\sum v^2\).
Cauchy–Schwarz gives \(V_{2,\infty}^2\le WD/16<10^{-28}\).
Each \(v<10^{-9}\) follows from \(v^3\le D/8\).
In \(r=v/(1+bv/n)\), \(|b|\le1\). For every \(n\ge2\),
the identity for \(r-v\) gives \(|2A-W|<2V_{2,\infty}\).
The mean-value bound for \((1+x)^{-2}\) on \(|x|\le10^{-9}/2\)
gives
\[
 V_3\le D/7,\qquad
 |P|\le |B_\infty|/4+3D/(8n).
\]
Here \((1-10^{-9}/2)^3>7/8\) also gives the bound on \(V_3\).
In particular \(2A\) lies in \([16/3,6]\), and
\(B_2\le V_{2,\infty}<10^{-14}\).
At \(n\ge10^7\) the second criterion in the theorem follows by
substituting these bounds.

For \(2\le n<10^7\), apply the
[quadratic skew estimate](MULTIVARIATE_QUADRATIC_SKEW_TANH_ANGLE_TRANSFER.md)
to both \(j=n,n+1\). Its exact clocks satisfy \(t_j\ge W/4>4/3\).
Since \(\pi<22/7\) implies \(\pi^3<128/3\), the chart parameters
are \(r_j=r=1/128\), \(K_j=3\). Every row rate
\(1/(2v)>6\), so \(S_j\le32r^{-4}D\) in that estimate.
With \(H=|B_\infty|+D<2\cdot10^{-24}\) and \(M=10^7\), we have
\[
 E_j\le2r^{-3}H+32r^{-4}D<r/M,\qquad
 Q_j\le2Mr^{-5}H^2+33r^{-4}D<1/(126M^2).
\]
The two angular errors sum to less than \(1/(63M^2)\), which
is below the available comparison gap
\(t_n/[36nN(1+t_n)]>1/(63nN)\).
For the pair of orders one and two, the quadratic is strictly
positive on the real axis by its representation as the expected
square of a nonzero real binomial. This separates it from
the linear root.

All scalar comparisons in the proof can be checked directly.
~~~python
from fractions import Fraction as Q

M, D, B = 10**7, Q(1, 10**28), Q(1, 10**24)
assert 6*D/16 < Q(1, 10**28)
assert D/8 < Q(1, 10**27)
assert (1-Q(1, 2*10**9))**3 > Q(7, 8)
assert Q(22, 7)**3 < Q(128, 3)
assert 262154*(B/4 + 3*D/(8*M)) + 46000000000*D/7 + Q(4, 10**14*M) < Q(1, 10**18)
r, H = Q(1, 128), 2*B
assert B + D < H
assert 2*r**-3*H + 32*r**-4*D < r/M
assert 2*M*r**-5*H**2 + 33*r**-4*D < Q(1, 126*M*M)
~~~

The strip contains asymmetric ordinary profiles. Take \(m=10^{14}\),
\(L=40m+1\), and increasing unit-gap rows with smaller exponents
\((3L-1)/2,(7L-1)/2\) and multiplicities \(27m,343m\).
Use \(250m\) decreasing unit-gap rows with smaller exponent
\((5L-1)/2\). They have \(W=216m/L\), \(B_\infty=0\), and
\(D=32m/L^3\), satisfying (2) by rational comparison.
