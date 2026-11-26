# Angular separation from matched tanh data near a hard critical core

Consider a fully shared ordinary two-monomial profile with unequal
rows retained. Fix any adjacent order \(n\ge2\), \(N=n+1\).
For each Gauss factor let \((T_i,T_i')=(T_{n,i},T_{N,i})\)
be the exact matched tanh pair, let \(\epsilon_i\) be its
orientation sign, and put
\[
 v_i=\left(\frac{n}{T_i'}-\frac{n-1}{T_i}\right)^{-1},\qquad
 W=2\sum_i v_i,\qquad B_\infty=-4\sum_i\epsilon_i v_i^2.
\]
These three quantities, including \(\max_i v_i\), are invariant
under the choice of adjacent orders. Let \(J\ge3\) be an integer
at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Theorem.** If
\[
   \max_i v_i>\frac13,\qquad |B_\infty|<\frac12,\qquad
                          W<4+10^{-24}J^{-4},                  \tag{1}
\]
then the profile has disjoint consecutive root arguments in
every degree. Hence it satisfies SFC(2) over every
characteristic-zero field. The residual rows may be asymmetric.

## Integer packing forces the exact critical pair

The
[matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
gives \(v=d/(2\ell+d)\) for each of the \(d\) factors of an
ordinary row of smaller exponent \(\ell\ge1\) and gap \(d\ge1\).
The row therefore contributes
\[
         W_{\rm row}=2dv,\qquad
         |B_{\infty,{\rm row}}|=4dv^2=2vW_{\rm row}.
\]
Call a row hard here when \(v>1/3\). A row with \(d=1\) has
\(v\le1/3\). Hard \(d=2\) forces \(\ell=1\), giving \(W_{\rm row}=2\).
Hard \(d=3\) forces \(\ell\le2\), giving \(W_{\rm row}\ge18/7\).
Hard \(d\ge4\) gives \(W_{\rm row}>8/3\).
In particular every hard row costs at least two units of \(W\).

Under the weaker bound \(W<33/8\), a hard \(d\ge3\) row
leaves less than \(33/8-18/7<2\) units of \(W\).
All remaining rows then have \(v\le1/3\). Thus their total
opposing skew has absolute value at most
\((2/3)(W-W_{\rm row})\), and
\[
   |B_\infty|>\frac23(2W_{\rm row}-W)
       >\frac23\left(\frac{36}{7}-\frac{33}{8}\right)
       =\frac{19}{28}.
\]
Consequently (1) excludes every hard row with \(d\ge3\).
Its required hard row has \(d=2,\ell=1\). If there is no
opposite-orientation hard row, all opposing skew comes from
rows with \(v\le1/3\), and
\[
   |B_\infty|\ge 2-\frac23(W-2)>\frac7{12}.
\]
The hypotheses therefore force both orientations \((1,3),(3,1)\).
They consume four units of \(W\), and the remaining rows have
\(W_{\rm tail}=W-4\ge0\).

## Angular transport from the exact core

The exact consecutive
[critical core reserve](MULTIVARIATE_MIRRORED_GAP_TWO_FORCE_RESERVE.md#an-angular-bound-at-the-ordinary-shape)
is greater than \(10^{-10}N^{-2}\). For the residual rows set
\(S_{\rm tail}=\sum d_i^2/(\ell_i+d_i)\), with the sum now over
rows. Since each midpoint is at most the larger exponent,
\[
          S_{\rm tail}\le W_{\rm tail}<10^{-24}J^{-4}.
\]
If the tail is empty, the exact core proof applies directly.

Otherwise, the
[asymmetric tail angular estimate](MULTIVARIATE_MIRRORED_UNIT_CORE_ANGULAR_TRANSPORT.md#asymmetric-ordinary-tails)
pairs each full root argument with some core root argument at
distance at most \(24\sqrt{S_{\rm tail}}\), at each degree
\(2\le m<J\), provided \(S_{\rm tail}<10^{-16}/(4J^2)\).
This part of the estimate uses the matched Gauss disks for the
tail and the Grace–Walsh–Szegő polarization of the core.
It applies to any circular core, since the choice of core enters
only through its angles and the final comparison with its gap.
The displayed bound on \(S_{\rm tail}\) implies the disk
hypothesis. For \(2\le n<N<J\) the two angular errors total
at most
\[
            48\sqrt{S_{\rm tail}}
             <48\cdot10^{-12}J^{-2}<10^{-10}N^{-2}.
\]
They therefore cannot exhaust the exact core reserve.
At orders one and two, the first root is negative real,
whereas the quadratic moment is strictly positive at every
real argument because the two distinct ordinary monomials
are linearly independent under the positive factorial integral.
For terminals \(N\ge J\) the absolute angular theorem applies.

This criterion uses the exact linked tanh rates and signed skew
together with integer row packing. The
[weighted matched collision](MULTIVARIATE_WEIGHTED_TANH_ANGLE_OBSTRUCTION.md)
satisfies the numerical inequalities in (1) and is excluded by
the row-packing step.

## A certificate from finitely many positive moments

One can also certify this ordinary neighborhood using only positive
mixed moments at any fixed lower order \(n\ge2\). Let \(M_1=\sum_i T_i\)
be the unsigned first moment and set
\[
 M_{\rm core}
       =2(n-1)\left(\frac1{2n}+\frac1{2n+1}\right),\qquad
 \varepsilon=\frac{10^{-24}}{6J^4}.
\]
For a row of gap two and smaller exponent one, \(v=1/2\) at its
two Gauss shifts \(h/d=1/2,1\). The exact lower tanh values are
\((n-1)/(2n)\) and \((n-1)/(2n+1)\), respectively. There is one
copy of each shift in each orientation in the critical core.

Require \(M_1<M_{\rm core}+\varepsilon\), and use the
[finite positive-moment extraction theorem](MULTIVARIATE_MATCHED_TANH_FINITE_POSITIVE_CORE.md)
with \(K=7\) and \(\eta=1/2\). Require that it recovers at
least one row of smaller exponent one and gap two in each
orientation. The first-moment condition itself verifies
the curvature cap used in the extraction. Indeed the joint
grid gives \(v_i\le (n+1)T_i/(n-1)\), so
\[
 W\le\frac{2(n+1)}{n-1}M_1
   <4(n+1)\left(\frac1{2n}+\frac1{2n+1}\right)+6\varepsilon
   <7 .
\]
After removing the two recovered rows, the nonnegative
remaining first moment is \(M_1-M_{\rm core}\).
Since \(v_i\le3T_i\), the residual curvature satisfies
\(W_{\rm tail}\le6(M_1-M_{\rm core})<10^{-24}J^{-4}\).
The core transport above proves angular separation at all orders.
The extraction permits private rows in the starting profile.
The residual bound excludes them, since a private row contributes
\(2d\ge2\) to \(W\).

This certificate needs only the extraction moment set and the
unsigned first moment. It also applies with bounded input error
when the extraction tolerance holds and the error bound on \(M_1\)
is smaller than the margin in its displayed upper test.
