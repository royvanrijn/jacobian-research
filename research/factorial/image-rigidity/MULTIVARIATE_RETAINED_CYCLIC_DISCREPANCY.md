# Cyclic discrepancy in retained factorial valuations

Fix adjacent orders \(n,N=n+1\) of an ordinary two-monomial
profile. Partition its rows into a selected group and retained
rows. The groups may have arbitrary independent heights.
Use the [common-slope rule](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md)
to compare the two full Newton polygons at a prime \(p\).
The remaining oscillation can be measured by cyclic interval
discrepancy, accounting for cancellation among retained rows.

## Exact retained oscillation

Choose \(H\ge0\) at least the largest \(p\)-valuation of any
integer in any full smaller-end-stripped retained ramp at either
order. Write \(d_i=b_i-a_i\) for the signed gap. For \(1\le e\le H\)
let \(q=p^e\). For each \(\delta\in\mathbb Z/q\mathbb Z\) form the
multiset of residues
\[
 B_{j,q,\delta}=\{j a_i\bmod q : i\text{ retained},
                       \ d_i\equiv\delta\bmod q\},\qquad j=n,N .
\]
For a multiset \(B\) of residues modulo \(q\), define
\[
  \operatorname{disc}_q(B)
    =\max_I\left|\#(B\cap I)-\frac{|I|}{q}\#B\right|,
\]
where \(I\) is a cyclic interval of residues and counts
include multiplicities. Put
\[
 D_{j,p}=\sum_{e=1}^H\sum_{\delta\bmod p^e}
                       \operatorname{disc}_{p^e}(B_{j,p^e,\delta}).
\]

**Theorem.** After subtracting the common slope
\[
                    \left(\sum_{i\ {\rm retained}} d_i\right)
                                    \sum_{e=1}^H p^{-e},
\]
times the coefficient index, the retained ordinate oscillation
at order \(j\) is at most \(D_{j,p}\). Thus a common factor degree
is bounded by the size of a maximum matching between the selected
slope lists with tolerance \(D_{n,p}+D_{N,p}\).

Indeed the unstripped argument of row \(i\) at index \(k\) is
\(j a_i+k d_i\). Removing its smaller-end factorial changes no
ordinate differences. By the floor formula for factorial
valuation, the difference of the retained contributions at
indices \(k,k'\), after removal of the displayed slope, equals
\[
 -\sum_{e=1}^H \sum_{i\ {\rm retained}}
 \left(\left\{\frac{j a_i+k d_i}{p^e}\right\}
                  -\left\{\frac{j a_i+k' d_i}{p^e}\right\}\right).
                                                               \tag{1}
\]
To justify truncation, for \(e>H\) the two arguments of any
one row lie in \([j\min(a_i,b_i),j\max(a_i,b_i)]\).
The floor changes only on passing a multiple of \(p^e\),
none of which occur in the full stripped ramp at that order.
Each omitted level therefore contributes zero.

For one residue class \(\delta\) at level \(q\), set
\(t=(k-k')\delta\bmod q\), represented in \([0,q)\).
Translate the multiset \(B=B_{j,q,\delta}\) by \(k'\delta\).
For each resulting residue, its fractional part increases by
\(t/q\), with an extra subtraction of one precisely on
\([q-t,q)\). Consequently its group's contribution to (1)
is the count in that length-\(t\) interval minus \(t\#B/q\).
In the original \(B\) this is a cyclic interval discrepancy,
of absolute value at most \(\operatorname{disc}_q(B)\).
Summing proves the oscillation bound.
The matching conclusion follows by centering the residual
ordinate interval at each order and applying the
common-slope rule. A maximum matching of size at most one
excludes a nonconstant rational common factor by even
positivity, as in that rule.

There is also a version using different truncations: if row
\(i\) is truncated at any height \(H_i\) at least its largest
full-ramp prime depth at the two orders, use at level \(p^e\)
only rows with \(H_i\ge e\), and use the signed sum of their
gaps in the linear slope for that level. The same proof
applies. In particular taking each \(H_i\) equal to the
largest prime depth recovers the setting of the
[per-row least common multiple bound](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md#arbitrary-retained-rows-at-one-prime).
Both bounds may be computed and the better one used at each
prime. In the applications there, \(R_p\) can be replaced by
\(\max(D_{n,p},D_{N,p})\) when this improves the charge.
The charged primes still come from the selected marks.

## Large exact additions preserving a fixed prime certificate

Fix a finite set \(\mathcal P\) of primes greater than an
absolute angular cutoff \(J\). Let \(d\ne0\) be any signed
integer gap. Choose \(H_p\) so that \(q_p=p^{H_p}\ge J|d|\)
and put \(Q=\prod_{p\in\mathcal P}q_p\).
One can append \(Q\) rows \((a,a+d)\), with bases exhausting
the residues modulo \(Q\), whose Newton slope effect at
every \(p\in\mathcal P\) is a common translation in all
orders \(j<J\). The rows' heights can increase independently
on arithmetic progressions.

To construct the rows, start with each base residue modulo
\(Q\). At \(p\) lift that base residue \(a_0\bmod q_p\) to
\(a_0+t q_p\bmod p q_p\). For any one order \(j<J\), the
full stripped ramp has length \(j|d|\le q_p\). Thus at most
one block offset meets zero modulo \(q_p\). At that offset,
exactly one last digit \(t\bmod p\) would give zero modulo
\(p q_p\), since \(p\nmid j\). There are fewer than \(p\)
orders, so a last digit avoiding all the forbidden values
exists. The Chinese remainder theorem chooses the last
digits simultaneously over all primes; each resulting
base can then be taken arbitrarily large in its residue
class modulo \(\prod_{p\in\mathcal P}p q_p\). This permits
positive endpoints for either sign of \(d\).
No full stripped ramp has a prime power beyond \(p^{H_p}\)
at these orders.

At every level \(p^e\) with \(e\le H_p\), the bases are
uniformly distributed modulo \(p^e\). Multiplication by
\(j\) permutes these residues, so every relevant multiset
\(B_{j,p^e,\delta}\) is uniform or empty. All discrepancies
are zero. The theorem gives the claimed exact common slope
translation. Arbitrarily many such packets with different
signed gaps and independently chosen heights can be added.
For opposite packets of equal absolute gap, take the same
\(H_p\) and counts for each orientation; their common
translations cancel.

In particular, suppose a selected profile has a prime
\(p\in\mathcal P\) at each low adjacent pair for which the
selected unit slope lists have shared width at most one.
Any such additions preserve this exact coprimality
certificate at the corresponding pairs. The first pair and
the pairs above the absolute angular cutoff are covered by
the [ordinary low order](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)
and [absolute angular](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
results. For example, take balanced opposite packets with
fixed gap size \(d_0>0\). For a common large baseline \(M\),
take all lifted bases in \([M,M+T]\), where \(T\) is the
fixed common CRT modulus. Taking
\(\lfloor\tau M/(2Qd_0^2)\rfloor\) copies of each opposite
packet adds curvature tending to any prescribed \(\tau>0\).
The exact invariance is independent of this choice of
baseline or multiplicity.
