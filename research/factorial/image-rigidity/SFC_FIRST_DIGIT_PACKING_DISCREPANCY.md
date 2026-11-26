# Exact first residue discrepancy under factorial row packing

Fix a prime \(p\) and distinct nonzero ordinary exponent
columns \(\nu\), with no constant term. Use the natural
coefficient variables \(u_\nu\), so the phase coefficient
of \(X^\nu\) is \(u_\nu/\nu!\). Its normalized factorial
moment series is
\[
                 E(T)=\sum_{n\ge0}
                       \frac{\Lambda(f^n)}{n!}T^n .
\]
Choose a high row \(A\) and a low row \(B\) and pack them
into \(p^e A+B\), with \(p^e>\max_\nu B_\nu\). Use the same
\(u_\nu\) after changing the exponent factorials. The
coefficient ring in this comparison is the universal
polynomial ring \(\mathbb F_p[(u_\nu)_\nu]\).

**Theorem.** If \(A\) is zero or the low row contains no
column of the form \(p^b e_B\), the original and packed
series agree modulo \(p\) at all orders. Otherwise put
\[
 b_*=\max_{\nu=p^b e_B} b,\qquad
 c_*=\min_{A_\mu>0}v_p(A_\mu),\qquad
                         n_*=p^{e+c_*-b_*}+1 .
\]
The series agree modulo \(p\) at all orders below \(n_*\)
and disagree at order \(n_*\). This is exact for the
universal coefficients. A specialization can remove
the exhibited discrepancy.

In particular the absence of low singleton columns in the
[all-order packing theorem](SFC_ALL_ORDER_PRIME_POWER_ROW_PACKING.md)
is necessary for universal residue congruence if the
high row is nonzero. The formula gives the exact
modulo-\(p\) cutoff for the individual packing operation
even when that condition fails. The argument does not
require any finite terminal search.

## Comparing count vectors by base-\(p\) digits

For a vector of monomial counts \(m_\nu\), let
\(w_\nu=\sum_i s_p(\nu_i)\) be the digit weight of the
column. Its native coefficient of \(\prod u_\nu^{m_\nu}\)
in order \(\sum m_\nu\) is
\[
   C(m)=\frac{\prod_i(\sum_\nu m_\nu \nu_i)!}
                  {\prod_\nu m_\nu!\prod_{\nu,i}(\nu_i!)^{m_\nu}} .
\]
This is an integer, as in the
[native moment recurrence](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#the-integral-partition-recurrence).
By the
[digit formula at residue precision](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#degree-bounds-in-characteristic-p),
it can be nonzero modulo \(p\) only if the following
two conditions hold.

1. Every count on a column of weight greater than
   one is less than \(p\).
2. Adding the expanded digits in each output row
   produces no carry. A count of a pure column
   \(p^b e_i\) contributes its own digits, shifted
   by \(b\). Any other counted column contributes
   its exponent digits multiplied by its one-digit
   count.

When those conditions hold, the residue is the
product of the factorials of the output digit loads,
divided by the digit factorials from the count and
column factors. For completeness, Legendre's formula
gives
\[
 (p-1)v_p(C(m))=
 \sum_\nu (w_\nu-1)m_\nu+\sum_\nu s_p(m_\nu)
                   -\sum_i s_p\!\left(\sum_\nu m_\nu\nu_i\right).
\]
The final digit sum is at most
\(\sum_\nu w_\nu s_p(m_\nu)\), with equality exactly
in the absence of carries in the expanded sums.
Thus the valuation is zero exactly under the two
conditions. The claimed residue formula follows
by Wilson's theorem applied to the \(p\)-free
factorial blocks. The signs cancel because the
total valuation is zero.

Packing moves the high row digits by \(e\) places
and superimposes them on the low row. The digits
of any individual column stay disjoint by the
hypothesis \(p^e>\max B_\nu\). Therefore every
column retains its digit weight and the same
digit factorials. The first condition is unchanged.
An original carry cannot disappear: on the packed
row each raw digit load is the sum of at most two
nonnegative original digit loads after moving the
high row. For a count vector with no carries in the
original profile, the residue is unchanged whenever
the two occupied sets of digit slots stay disjoint.

## The first overlap

In an original carry-free count vector the low
row contributions from all columns of weight greater
than one have positions less than \(e\). Their
counts are single digits and their low exponent
entries are less than \(p^e\). No new carry is made
from those positions. The only way to reach a
high position in the low row is therefore by the
shifted digits of the counts of columns \(p^b e_B\).
This proves all-order equality when such columns
are absent. If \(A\) is zero the packing changes
no nonzero row at all.

When the two kinds of columns required in the
theorem exist, a moved high digit lies at position
at least \(e+c_*\). Indeed each high exponent
is divisible by \(p^{c_*}\), and expanding counts
can only move digits upwards. If an original
carry-free count vector overlaps a high digit,
some low pure count must be at least
\(p^{e+c_*-b_*}\). There must also be a nonzero
count on a different column carrying a high
digit. This proves absence of all discrepancies
below \(n_*\).

At \(n_*\), take the count \(p^{e+c_*-b_*}\) on
a low column with position \(b_*\) and count one
on a column \(\mu\) with \(v_p(A_\mu)=c_*\).
All other counts are zero. In the original
profile the low pure count contributes exactly
one digit at position \(e+c_*\). It is separate
from the low digits of \(\mu\), so the original
count vector is carry-free and its residue is one.
In the packed profile that digit instead meets
the nonzero lowest digit \(a\) of \(A_\mu\).
The residue becomes \(a+1\) modulo \(p\): if
\(a=p-1\) a carry makes it zero, and otherwise
the change in the slot factorial is
\((a+1)!/(a!\,1!)\). This differs from one and
occurs at the displayed order. Distinctness of
the packed columns follows by reading the low
digits and then the high digits, so no identification
of universal coefficient variables is needed.

## Minimum dimension under separated digit packing

Discard zero rows and allow a partition of the other
rows into groups. In each group sum the rows after
placing their exponent entries in disjoint digit
intervals. Write \(q\) for the number of original
rows containing at least one pure \(p\)-power column.
The minimum possible number of packed rows preserving
the universal series at every order modulo \(p\) is
\(\max(1,q)\).

For the lower bound, if two pure-column rows end
in the same packed row, retain just one pure column
from each. Their packed exponents are distinct
powers \(p^a,p^b\), since their digit intervals
are disjoint. For \(a<b\), the count \(p^{b-a}\)
on the lower column and count one on the other
produce two loads in the same final slot. The
residue changes from one to two, read modulo \(p\),
as in the first-discrepancy proof. All remaining
column counts are zero, so no other row operation
affects this coefficient. This proves necessity
even if intermediate series were allowed to change.

For attainment, keep the \(q\) pure-column rows
separate and pack each other row below one of
them without overlap. If \(q=0\), pack every row
into one. A singleton column remains a singleton
under these operations, and a column of greater
digit weight retains that weight. Thus each packing
has a low row without a pure column, and the theorem
preserves the universal residue series at every
step.
