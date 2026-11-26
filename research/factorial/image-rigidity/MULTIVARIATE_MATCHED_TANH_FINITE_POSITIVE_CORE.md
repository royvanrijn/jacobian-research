# Recovering the large-gap ordinary rows from finitely many positive matched tanh moments

Fix an integer lower order \(n\ge2\), a positive rational curvature cap \(K\),
and a rational \(0<\eta\le1\). On each unequal ordinary row let \(d\) be
its positive gap and \(\ell\) its smaller exponent, and set
\[
 v=\frac{d}{2\ell+d},\qquad
 W=\sum_{\rm rows}\frac{d^2}{\ell+d/2}=\sum_{\rm rows}2dv .
\]
Include private rows by allowing \(\ell=0\). At each of the \(d\) Gauss
factors let \(t=T_n\), \(s=T_{n+1}\), with the orientation of the row.

**Theorem.** On all ordinary profiles with \(W\le K\), a computable finite
set of positive mixed moments
\[
                    \sum_{\text{one orientation}} t^p s^q,
                  \qquad p\ge1,\quad q\ge0 ,
\]
determines both orientation counts of every row with \(v\ge\eta\).
There is a computable positive tolerance on the input moments for exact
recovery of all these integer counts. The moments can equivalently be
supplied as the unsigned and signed moments over both orientations.

## Finite support above a positive tanh cutoff

The [joint Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md)
gives
\[
 t=\frac{(n-1)v}{n+(2h/d-1)v},\qquad
 s=\frac{nv}{n+1+(2h/d-1)v}.
\]
In particular \(v/3\le t\le v\le1\), and \(\sum t\le W/2\).
Take \(b=\eta/6\). If an ordinary atom has \(t\ge b\), then its
row has \(v\ge b\) and \(2dv\le K\), so
\(d\le K/(2b)\) and \(\ell\le d/(2b)\).
Enumerate the distinct rational points \(A\) with \(t\ge b\) from
these finitely many rows, keeping the row cap \(2dv\le K\).
Duplicate points contribute to the same atom count.

For \(x=(t_x,s_x)\in A\) with \(t_x\ge2b\), put
\[
 L_x(t,s)=\prod_{y\in A\setminus\{x\}}
              \frac{(t-t_y)^2+(s-s_y)^2}{(t_x-t_y)^2+(s_x-s_y)^2},
 \qquad F_x(t,s)=(t/t_x)^m L_x(t,s).
\]
The polynomial \(L_x\) is nonnegative everywhere on the real plane.
It is one at \(x\) and zero at the other listed points. On the unit
square it is bounded by the computable rational
\[
               C_x=\prod_{y\in A\setminus\{x\}}
                     \frac{2}{(t_x-t_y)^2+(s_x-s_y)^2}.
\]
Since the remaining atoms have \(t<b\), in either orientation
\[
 0\le \sum F_x(t,s)-\#\{(t,s)=x\}
      \le \frac{C_x K}{2t_x}\left(\frac{b}{t_x}\right)^{m-1}.
\]
Use \(\sum t\le K/2\) in the last step. Choose \(m\ge1\) making the
last bound less than \(1/4\). The coefficients of \(F_x\) are rational,
and each monomial in it has positive \(t\)-degree. The linear
combination of the input moments that gives \(\sum F_x\) then recovers
the integer atom count by rounding to the nearest integer. For
perturbed input it suffices to bound the absolute error of every
input moment by \(1/(4\max_x\|F_x\|_1)\), where \(\|F_x\|_1\)
is the sum of the absolute coefficients. If there are no target
points, the desired row counts are all zero.

Every atom of every row with \(v\ge\eta\) has \(t\ge2b\). At each such
\(v\), recover \(h/d\) from the inverse joint formula and apply the
finite divisor inversion for the row multiplicities. Any row sharing
these grid points has the same \(v\). The recovered atoms therefore
give the complete data used by that inversion, regardless of the
sizes or multiplicities of the rows with \(v<\eta\).

## A consequence for weighted test profiles

The same argument separates a fixed finitely supported weighted
profile from ordinary profiles with \(W\le K\) whenever it has an
atom of positive \(t\) with noninteger mass in either orientation.
Use a smaller \(b\) such that the target has \(t_x\ge2b\), and add
the weighted profile's points with \(t\ge b\) to \(A\). If the target
is absent from the ordinary list, its ordinary count is zero.
For the fixed weighted profile the same tail estimate uses its own
finite first moment in place of \(K/2\). By taking \(m\) large,
both errors can be made smaller than one quarter of the distance
from the fixed mass to the nonnegative integers. The single
linear combination of positive moments then separates the values
with a positive margin. For a fixed mass that is positive but
located at an inadmissible grid point, the same argument uses its
distance from zero. This applies to every fixed profile in the
[finite inverse-moment collision construction](MULTIVARIATE_MATCHED_TANH_FINITE_INVERSE_MOMENT_OBSTRUCTION.md),
whose atom at \(v=15/16,\ h/d=1\) has positive fractional mass in
each orientation. The argument depends on the curvature cap and
the observed hard atoms and gives no angular estimate for the
remaining diffuse profile.
