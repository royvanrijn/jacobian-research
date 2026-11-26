# Simultaneous factorial records at finitely many primes and orders

Fix a support size \(s\ge2\), a finite nonempty set of positive
moment orders \(\mathcal J\), and precisions \(\ell_p\ge1\)
at a finite nonempty set \(\mathcal P\) of primes.
For a list of ordinary exponent rows \(a_i\in\mathbb Z_{\ge0}^s\)
put
\[
 R_\nu=\prod_i\frac{(\nu\cdot a_i)!}{(j a_{i,1})!},
 \qquad j=|\nu|\in\mathcal J,\qquad
 \Delta_c=\sum_i(a_{i,c}-a_{i,1}).
\]
These are the constant-normalized coefficients after removing
the fixed multinomial coefficients. We omit \(\nu=j e_1\)
since then \(R_\nu=1\). A record consists of the valuations
\(V_{p,\nu}=v_p(R_\nu)\) and unit parts
\(U_{p,\nu}=p^{-V_{p,\nu}}R_\nu\bmod p^{\ell_p}\).

For each \(p\) group the coordinates by the pair
\[
 \left(\nu/j,\quad j/p^{v_p(j)}\right).
\]
Let \(C_p\) be the number of these classes, and put
\[
 Q=\operatorname{lcm}_{p\in\mathcal P}(2p^{\ell_p-1}),\qquad
 M=\operatorname{lcm}\left(Q,\ \{\varphi(p^{\ell_p}):p\in\mathcal P\}\right).
\]

**Theorem.** Every ordinary record at these primes and orders
can be reproduced in the original coefficient variables and
with the same signed column differences by a shared ordinary
profile with at most
\[
                     (Q^s-Q)(M-1)+1+2\sum_{p\in\mathcal P}C_p          \tag{1}
\]
rows. The row count is independent of the original exponent
sizes and the sizes of the recorded valuations.

The same bound applies to the number of rows that must be
appended to an arbitrary fixed core to impose any attained
target record with any specified target signed column
differences. When the core and target have the same signed
differences, the appendage can have arbitrarily small curvature
in the binary case, and its constant-normalized real
coefficient multipliers can be arbitrarily close to one
through any prescribed finite range of orders. In the binary
case two additional rows suffice to ensure that the resulting
profile has disjoint consecutive root arguments in every
degree.

The exact compatibility conditions on a prospective record are
the scaling identities in (2) and (3) below, with one common
finite table of row residue counts. Thus simultaneous records
at arbitrary finite collections of moment orders have a finite
description, even when the original profiles have unbounded
dimension. In contrast to the
[one-prime packing theorem](SFC_PADIC_BOUNDED_ROW_PACKING.md),
the bound (1) depends on the chosen orders. Here the same
profile works at all the given primes with no variable rescale.

## Compatibility along the prime-power chains

Two coordinates are in the same class at \(p\) precisely when
one is \(p^e\) times the other, for an integer \(e\ge0\).
For \(j=|\nu|\) write \(D_e=1+\cdots+p^{e-1}\) if \(e\ge1\).
Iterating the factorization of factorials into \(p\)-divisible
and \(p\)-unit factors gives
\[
        V_{p,p^e\nu}-V_{p,\nu}=D_e\sum_c\nu_c\Delta_c.              \tag{2}
\]
In the unit group \(G_p=(\mathbb Z/p^{\ell_p})^\times\), let
\[
 F_p(t)=\prod_{\substack{1\le u\le t\\p\nmid u}}u,\qquad
 D_{p,e}^{\,\mathrm{unit}}(x)=\prod_{h=1}^e F_p(p^h x)\bmod p^{\ell_p}.
\]
The same factorial factorization gives
\[
 \frac{U_{p,p^e\nu}}{U_{p,\nu}}
   =\prod_i\frac{D_{p,e}^{\,\mathrm{unit}}(\nu\cdot a_i)}
                  {D_{p,e}^{\,\mathrm{unit}}(j a_{i,1})}.       \tag{3}
\]
The product of the units in a complete block of length
\(p^{\ell_p}\) has square one. Increasing \(x\) by
\(2p^{\ell_p-1}\) in \(F_p(p^h x)\) adds \(2p^{h-1}\)
complete blocks, for every \(h\ge1\). Consequently each
factor on the right of (3) depends only on the row modulo \(Q\).

In particular the compatibility test for a prospective record
and signed difference vector \(\Delta\) is (2) together with
the existence of counts \(m_u\in\{0,\ldots,M-1\}\) for
\(u\in\{0,\ldots,Q-1\}^s\) such that (3) holds with
these residue rows and multiplicities and
\[
                \sum_u m_u(u_c-u_1)\equiv\Delta_c\pmod Q .
                                                                    \tag{4}
\]
Necessity follows by reducing the actual row counts modulo \(M\).\nResidues with every entry equal contribute trivially to both (3)\nand (4), so their counts can be set to zero.
It is enough to compare any one representative with the
remaining members in each class.

To extend a fixed ordinary core to a prescribed compatible
target, work on the quotient of their local \(R_\nu\) data,
subtract the signed column differences, and subtract the two
tables modulo \(M\). The resulting data have the same
compatibility conditions. Thus it suffices to construct an
appendage for a compatible record.

Take \(m_u\) rows of residue \(u\), all shared. Next take at most
one row with all entries divisible by \(Q\) to correct the
signed differences exactly. This is possible by (4), using a
large first exponent and the required differences for the
other entries. This row contributes unit quotient one in (3).
The rows so far realize the signed differences and all the
prescribed between-order identities. After dividing out their
record, the residual valuations and units are constant
on each class at the respective prime. We will correct these
residual coordinates while adding zero signed difference.

## Bounded interpolation on the remaining classes

Choose a vector \(d=(d_1,\ldots,d_s)\) of distinct nonnegative
multiples of \(Q\), with \(d_1=0\), such that \(r\mapsto r\cdot d\)
is injective on the finite set of vectors \(r=\nu/j\),
including \(e_1\). Such a choice is obtained by taking
sufficiently separated powers of an integer for the last
\(s-1\) entries, and multiplying by \(Q\). Put
\(d_{\max}=\max_c d_c\) and \(t(r)=r\cdot d/d_{\max}\).
All tested coordinates have \(0<t(r)\le1\). By multiplying
\(d\) by a common positive integer, ensure that for each
tested \(r\) and each prime-free index \(b\) under consideration
there is a rational \(\alpha\) of reduced denominator \(b\)
with
\[
  0<\alpha<d_{\max},\qquad
       t(r_{\mathrm{prev}})<\alpha/d_{\max}<t(r),
\]
where \(r_{\mathrm{prev}}\) is the preceding point in the order
of \(t\), with \(t(e_1)=0\) at the start.
In fact the admissible values \(\alpha=m+1/b\) have spacing
one, so it suffices to make each original open interval for
\(\alpha\) longer than one.

At prime \(p\) reserve one pair of rows for each class \((r,b)\)
with \(b=j/p^{v_p(j)}\), using such a pole \(\alpha\).
The two rows have components
\[
              a_++d_c,\qquad a_-+d_{\max}-d_c .
\]
For the first row the ratio that defines \(R_\nu\) samples
factors \(ja_++h\) for \(1\le h\le \nu\cdot d\).
For the second row it uses the inverse of the final
\(\nu\cdot d\) factors of \(ja_-+1,\ldots,ja_-+jd_{\max}\).
At \(p\) choose large exact depths \(H_+,H_-\) of
\(a_++\alpha,a_-+d_{\max}-\alpha\), with independently
chosen unit parts \(c_+,c_-\) modulo \(p^{\ell_p}\).
The denominator \(b\) is prime to \(p\), so both prescribed
classes are ordinary congruences. For the first row,
\[
                     ja_++h=j(a_++\alpha)+(h-j\alpha).
\]
Unmarked offsets have fixed valuations and units once
\(H_+\) exceeds \(\ell_p\) plus the maximum of their valuations
in the finite collection of ramps. Marked offsets occur
at \(h=j\alpha\) when \(b\mid j\). The second row has the same
description with pole \(d_{\max}-\alpha\). The two marked
factors occur in the respective products at exactly the
same coordinates, since no pole lies at a tested threshold.
At every other monitored prime \(q\), impose \(a_\pm\equiv0\)
modulo a large enough power of \(q\) that those valuations
and units are exactly the ones at base zero. The Chinese
remainder theorem combines all these conditions.

A given pair produces the mask
\[
       w_{r,b}(r',b')={\bf1}_{\,t(r')\ge t(r)}{\bf1}_{\,b\mid b'} .
\]
Relative to reference depths \(H_+=H_-\) sufficiently large
and reference units \(c_+=c_-=1\), the change in valuation
is \((H_+-H_-)w_{r,b}\), and the unit multiplier is
\((c_+/c_-)^{w_{r,b}}\). Contributions at the other primes
have not changed. Each pair has zero signed column
differences, and all its differences within a row are
multiples of \(Q\). Its aggregate valuation drift in (2)
is therefore zero and each of its unit quotients in (3)
is one, for any chosen depths and units. This also holds
for the reference rows.

Order the classes at a prime first by \(t(r)\), then by
increasing \(b\). The matrix of the masks on these classes
is lower triangular with diagonal one. Within a block of
equal \(r\), the off-diagonal entries are the divisibility
indicators \(b\mid b'\). Thus the masks form an integral
basis. They also give a bijective parametrization of
multiplicative unit corrections, by the same triangular
elimination in the finite abelian unit group. Every
prescribed integer difference of depths is realized
with both depths above their thresholds. No repeated
rows are required to represent large integers. This
corrects all coordinates, proves sufficiency of the
compatibility test, and gives the row bound (1).

When no two distinct tested orders are related by a power of
any monitored prime, there are no between-order identities.
The residue calibration is unnecessary. Use at most one row
for the gap repair and the same marker construction, with the
bound \(1+2\sum_p C_p\) independent of precision. This includes
every adjacent pair starting at order at least two. The argument
then refines the [adjacent approximation theorem](MULTIVARIATE_REAL_LOCAL_COEFFICIENT_APPROXIMATION.md)
by bounding the number of rows independently of the prescribed
local valuations or precisions.

For prescribed data with unknown signed differences, (2) gives
linear integer equations in \(\Delta\). Enumerate the finite
residue tables allowed by (3) and solve these equations together
with (4). Thus realizability is decidable in this case as well.

## Real approximation and all-degree regularity

Translate both entries of each binary row, or all entries
of any row, by a common large positive integer \(T\). At
each monitored order every factorial quotient changes by
a product of factors \(1+jT/h\) or their reciprocals
over its original shorter factorial interval. If \(T\) is
divisible by sufficiently high powers of the monitored
primes, all these factors are units congruent to one at
the chosen precisions. We can take \(T\) arbitrarily large.

In the real constant-normalized coefficients, a fixed
new row with exponent vector \(a\) then contributes
\((jT)^{\nu\cdot a-j a_1}(1+O(T^{-1}))\) at each fixed order
and coordinate. For zero net signed differences these powers
cancel over the whole appendage. For binary rows the total
added curvature tends to zero even if their signed gap is
nonzero. Applied to the entire compressed binary profile,
this also gives all-degree angular separation by the
[uniform weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md).

For binary appendages preserving an arbitrary fixed core,
use the [universal opposed-row extension](MULTIVARIATE_UNIVERSAL_OPPOSED_PADIC_ANGLE_MARKING.md)
at one auxiliary prime to guarantee all-degree angular
separation. Choose its equal gap \(d_*\) divisible by \(Q\)
and at least the absolute angular cutoff. Fix a reference
shared mirrored pair of this gap and divide its local
multipliers out of the desired record before interpolation.
Its signed difference is zero and its local unit
quotients in (3) are one. The compatibility conditions
are preserved. After building the interpolation rows,
choose the pair's common baseline congruent to the
reference baseline modulo large powers of the monitored
primes, and sufficiently close to \(-1\) at the auxiliary
prime as in the extension theorem. The Chinese remainder
theorem again permits arbitrary real height while keeping
all local multipliers of the reference pair. This gives
the stated two-row overhead and retains real approximation
in the zero-gap case.

The finite compatibility test describes constant-normalized
local factorial data at all tested orders together. Any
such ordinary binary record can coexist, at those same
primes and orders, with a profile that satisfies SFC(2)
in every degree. Exact coprimality of an arbitrary original
profile still requires information beyond these finitely
many valuations and residues.

In particular, at signed gap zero one can prescribe \(R_\nu\)
to have valuation zero and unit part one at every tested
coordinate. Use the zero residue table. This realizes the
binomial polynomials locally at all tested orders at once,
with all-degree angular separation over characteristic zero.
At an arbitrary signed gap one may instead prescribe arbitrary
local data on one adjacent pair and freely extend the targets
to any other tested orders. No chain has two coordinates
from different orders of that pair, since their ratio is
never an integer prime power for adjacent orders starting
at two. Choose the residue table of any ordinary profile of
the desired gap and use (2) and (3) to fill each chain from
its one prescribed coordinate, or any representative if
none is prescribed. In both cases the appendage may retain
any fixed core with that same gap and approximate its real
coefficients.

## Two rows suffice to preserve an existing binary record

For an existing binary profile, one can keep its complete
tested record and signed gap while ensuring all-degree angular
separation by adding only two shared rows. Fix an auxiliary
prime \(q\notin\mathcal P\), and choose an even gap \(d\) at least
the absolute angular cutoff. Require \(d\) to be divisible by
\(p^{A_p}\) for all monitored primes, where
\[
 A_p\ge \ell_p+\max_{1\le u\le \max\mathcal J}v_p(u).
\]
Take rows \((a,a+d)\), \((H-1-a,H-1-a-d)\) with \(H>a+d+1\).
Choose \(a\) divisible by the same \(p^{A_p}\) and congruent
to \(-1\) modulo the power \(q^E\) given by the
universal opposed-row extension for this gap and core.
Then take \(H\equiv a+d\pmod{q^E}\). This makes the second
smaller exponent \(H-1-a-d\) congruent to \(-1\) as well.

Write \(x=ja+kd\) and \(x_0=ja\). At a monitored prime
choose \(v_p(H)\) large enough that \(jH/h\) has valuation
at least \(\ell_p\) for every integer \(1\le h\le x+j-1\)
in the tested orders. These divisibilities are simultaneous
with the auxiliary congruence. The shifted version of the
[complementary factorial identity](MULTIVARIATE_MIRRORED_COMPLEMENTARY_ENDPOINT_CONGRUENCE.md#exact-paired-complementary-factors)
is
\[
 \frac{x!(jH-x-j)!}{(jH-1)!}
     =\frac{(-1)^{x+j-1}}
        {\prod_{u=1}^{j-1}(x+u)
         \prod_{h=1}^{x+j-1}(1-jH/h)} .
\]
The second product consists of units congruent to one. In
the first product each factor has the same valuation and
unit part as \(u\), since \(x\) is divisible by \(p^{A_p}\).
Also \(x-x_0=kd\) is even. Hence the ratio between the
displayed expression at \(x\) and at \(x_0\) has valuation
zero and unit part one modulo \(p^{\ell_p}\). This proves
exact preservation of the entire monitored record.

Finally translate both smaller exponents by an arbitrarily
large common multiple of \(q^E\) and of the moduli needed
to preserve the monitored data. This fixes the auxiliary
congruences. The real coefficient multipliers tend to one
at every fixed order and the added curvature tends to zero.
The opposed-row theorem gives all-degree angular separation
for the new ordinary profile.
