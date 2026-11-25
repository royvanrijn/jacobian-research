# Terminal fourteen on the fully shared upper bridge

Let the two monomials be distinct and have equal support. At order
thirteen set \(x_i=|a_i-b_i|/(13\min(a_i,b_i)+1)\) for each unequal
shared row, and let \(u\) be the product of the row endpoint mean
factors. Then the order-thirteen and order-fourteen moment polynomials
have disjoint zero arguments whenever
\[
                 \frac{19}{25}\le u\le\frac{171}{200}.
\]
As usual, each degree uses its own positive endpoint scaling. Equal
rows only multiply a moment by a positive scalar.

The [high-curvature theorem and retained profiles](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
cover profiles with some \(x_i>1/8\). The
[direct upper-bridge estimate](MULTIVARIATE_TERMINAL_FOURTEEN_HIGH_U.md)
covers profiles with every \(x_i\le1/60\). The residual complement is
proved below. Together with the adjacent
[middle retained](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_MIDDLE.md),
[middle direct](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_DISKS.md)
and [upper retained](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_UPPER.md)
theorems, this covers the fully shared range \(16/25\le u<1\).

## Retained profiles and ray replay

Retain every row with \(x_i>1/60\). A row of base \(a\) and gap \(d\)
has squared mean factor
\[
           B(a,d)=\prod_{h=1}^d\frac{13a+h}{13a+12d+h}.
\]
The exact enumeration in the
[terminal-thirteen retained-band argument](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md)
applies at lower order thirteen with \(1/60<x_i\le1/8\) and
\(\prod B(a,d)\ge(19/25)^2\). It lists thirteen nonempty oriented
multisets modulo simultaneous column reversal. If the
squared product of the listed rows is \(B\), the residual mean is
\(v=u/\sqrt B\) intersected with \(v\le1\) and enclosed outwards.
The exact row factorial products of the retained list are kept in
each coefficient model.

Twelve lists use the
[grouped ray certificate](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md)
with the terminal-fourteen order constants. A cell of residual
width at most \(1/1000\) with lower endpoint at least \(17/20\) can
instead use the [signed disk certificate](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_UPPER.md).
Intersecting with the loss bound gives a smaller residual reciprocal
cap wherever applicable. The ordinary loss-ratio model explicitly
enumerates residual reciprocals above \(1/40\) for the first six
lists, and above \(1/150\) for the final six lists. Below those floors
it uses the continuous row bound.

## Centered signed disks for the remaining list

The retained singleton has \((a,d)=(4,1)\). Its residual enclosure is
\([526/625,9469/10000]\). For each order let \(T_i\) be the residual
tanh variables, bounded by \(T_*\), and let \(W\ge\sum_i
\operatorname{artanh}T_i\). Write
\[
 M_r=\sum_i\varepsilon_i T_i^r,\qquad V_r=\sum_i T_i^r,\qquad P=-M_2/2.
\]
The sign of \(P\) is the convention for the logarithmic directions
used in the checker. Transport \(P\) between orders using the interval
in the signed disk certificate above.

For \(l=2,3,4\) choose \(\gamma_l=3/4,2/3,5/8\), respectively.
Then
\[
 \left|-\frac{M_{2l}}{2l}-\frac{\gamma_l}{l}T_*^{\,2l-2}P\right|
       \le \frac{1-\gamma_l}{2l}T_*^{\,2l-1}W.
\]
Indeed, with \(q=2l-2\), the maximum of
\(y^{q+1}-\gamma_l y\) on \([0,1]\) is \(1-\gamma_l\).
The magnitude of its interior minimum is at most \(1-\gamma_l\)
by the exact inequality
\[
 q^q\gamma_l^{q+1}\le(q+1)^{q+1}(1-\gamma_l)^q.
\]
Apply the resulting absolute bound to \(y=T_i/T_*\), multiply by
\(T_iT_*^{2l-1}\), and sum using \(\sum_iT_i\le W\).
When \(T_*=0\) the moment bound holds directly. The centered directions
are included in the common parameter \(P\), preserving their
cancellations when shifting the polynomial for the disk test.

The unsigned moment intervals also depend on \(P\). If
\(S=\sum_iT_i>0\), Jensen's inequality with weights \(T_i/S\) and
\(\sum_iT_i^2\ge2|P|\) gives
\[
             \frac{(2|P|)^{2l}}{W^{\,2l-1}}
        \le V_{2l+1}\le T_*^{\,2l}W,\qquad 1\le l\le4.
\]
The zero-loss endpoint is included by setting the moments to zero.
For each signed cell the lower bound uses the least possible
\(|P|\). The upper bound and the tails of the logarithmic series
use the outward loss and tanh bounds for that order. The signed
interval is bisected as needed to depth six. The same quartic
grouped expansion and certified Rouché and angular inequalities as
in the upper retained disk proof finish the singleton.

## Exact replay

The [checker](verify_multivariate_terminal_fourteen_retained_bridge.py)
verifies thirteen lists, 449 outer cells, 68,912 radial cells and 211
outer disk covers. The latter includes 160 centered singleton cells,
split into counts \(32,64,64\) at \(v=87/100\) and \(v=91/100\).
The checker uses the audited 15-slot native projector for the ray
cells and outward ball arithmetic for the disk cells.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_retained_bridge.py
~~~

To run selected retained lists, append their zero-based indices.
