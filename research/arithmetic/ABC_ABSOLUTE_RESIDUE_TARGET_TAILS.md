# Absolute residue counts for popular complete targets

When a full divisor on one summand is shared by many
targets, the opposite summands give a dictionary of target
residues. The target's own radical cap charges that
dictionary at the square root of its size. This gives a
cross-target estimate with no requirement that the divisor
reach the target's full numerical height.

## One full divisor and the multiplicity threshold

Let \(n=X^L\), and consider original positive primitive
rows \(A+B=C\), \(X/2<C\le X\). Keep any prescribed
eligibility conditions, including a least-repair cutoff.
Fix nonnegative \(u,\alpha,r_C,b\) and positive \(w,\kappa\).
Use a bin in which every complete target has between
\(n^{u-\tau}\) and \(n^{u+\tau}\) eligible distinct rows.
Suppose the complete \(B\) dictionary has cardinality at
most \(n^{\alpha+\tau}\), and the target radical is at
most \(n^{r_C+\tau}\).

Give each row a record \(D\mid A\). The record dictionary
has size at most \(n^{b+\tau}\), and the record is a full
divisor in the sense that \(D\ge X^w\) and it is
\(\lfloor\kappa/L\rfloor\)-full. A prime from \(D\) need
not account for its entire power on \(A\). All quantities
are fixed before the height limit, except the integer
records and the members of their dictionaries. The
roles of \(A,B\) may be exchanged.

**AR1.** For fixed \(\eta>0\), at sufficiently small
fixed bin width \(\tau>0\) and all sufficiently small
fixed \(L>0\), write
\[
        \beta=\max\{0,\alpha-(u-b)_+\}.
\]
The count of complete targets in the bin is
\[
                   O\big(n^{\,b+(r_C+\beta)/2+\eta}
                                    X^\varepsilon\big).        \tag{AR1}
\]
In every nonempty sequence of these bins one has
\(\alpha\ge u-o(1)\) because the opposite summands in
a target fibre are distinct. The maximum with zero
only handles allowances or empty parameter ranges.

For every target choose a most frequent marked \(D\)
among its own witnesses. It occurs on at least
\(n^{(u-b)_+-O(\tau)}\) rows, interpreting a lower bound
below one by one. Partition the targets by the
selected record. At a fixed \(D\) all those witnesses
have \(B\equiv C\bmod D\), and primitivity makes this
a unit class. A single target requires at least
\(n^{(u-b)_+-O(\tau)}\) distinct members of the
complete \(B\) dictionary in its class. Thus at most
\(n^{\beta+O(\tau)}\) unit classes can contain such a
target. Apply the degree-one case of (AT1) to the
target integers, with modulus \(D\), no excluded
prime support, and that unit residue dictionary. It
gives exponent \((r_C+\beta)/2\) before allowances.
Sum over \(D\). The cutoff in (AT1) depends on the
fixed numerical lower bound \(w\), which is selected
before taking \(L\) small.
[Absolute residue pricing, degree one of (AT1)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md#algebraic-integers-with-a-ray-residue-dictionary)

For clarity, the degree-one estimate tests the full
integers rather than the separate residues of their roots.
In the usual exact-layer boxing, expand a fixed
\(2h\)-moment for the character sum of one degree-\(e\)
root position. The two collapsed products each have
size \(X^{O(L)}\) and are units modulo \(D\). Agreement
of their \(e\)-th powers puts their ratio among the
\(e\)-th roots of unity modulo \(D\). The size of that
group is height-subpower since \(e\) and the moment
orders are fixed at the given budget. Each product is
smaller than \(D\) for sufficiently small fixed \(L\),
so each root of unity fixes at most one other product.
This gives the diagonal moment bound. Cauchy–Schwarz
on the sum of the allowed unit classes then contributes
the square root of the number of classes. The moment
orders for the layer sums allow Hölder on the remaining
character weight. This is why the bound need not price
individual target residues independently.

If \(A\) has positive nonatomic mass in a compact
positive scaled-degree band, an arbitrarily small fixed
\(b>0\) is available. On an anatomic box select original
degree layers in a continuity window with positive
nonatomic mass and total radical mass smaller than \(b\).
Their full powers have numerical weight bounded below
by a positive constant after choosing the window, and
the radical cap bounds the record dictionary by
\(n^{b+o(1)}X^\varepsilon\). The window can have fixed
positive width. No bound on the number of root positions
in the marked packet, or on the size of the packet's
largest prime, is needed for this argument. The full
rank of a CRT list is also unnecessary.
[Selection of affordable full divisors from diffuse mass](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md#completing-all-balanced-profiles-after-small-modulus-support-recovery)

In a balanced near-saturating sequence with \(r_C=1/3\),
letting \(b\) tend to zero after taking the height limit
gives the necessary dictionary inequality
\[
                  \alpha\ \ge\ \frac13-u.                    \tag{AR2}
\]
Indeed the row exponent of the bin is at most
\((r_C+\alpha+u)/2+o(1)\), by multiplying the target
count by its maximal fibre. A fixed deficit in
\(\alpha+u\) below \(1/3\) gives a strict coefficient
saving for this bin, even if all three numerical
heights are one and all observed scaled degrees lie
in a bounded diffuse band. At the occurrence endpoint
\(u=0\) any fixed compression of this opposite
dictionary below \(1/3\) suffices. At \(u=1/6\),
(AR2) agrees with the lower bound \(1/6\) of the
fixed-letter method, so that endpoint still needs
another input.
[Fixed-letter fibre bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## A second charge from the target's nonatomic mass

Suppose in addition that the boxed exact valuation
measures of the targets converge in the outer limit
to a fixed measure \(\nu_C\). Write \(d_C\) for its
nonatomic mass. For a fixed full record \(D\) and
at most \(n^\beta\) prescribed unit classes, the target
count has exponent at most \(d_C+\beta+\eta\), with
an arbitrarily small fixed allowance \(\eta\).

To see this, select finitely many positive atoms
containing all but a small mass of the positive atomic
part. Enclose them in slowly shrinking windows and
include one shrinking window at zero. Using weak
convergence and the uniform numerical moment bound,
fix all the exact factors outside these windows at
cost \(n^{d_C+\eta_0+o(1)}X^\varepsilon\), where the
positive atom list omits mass at most \(\eta_0\).
Choose that list after prescribing the desired fixed
allowance and before shrinking its windows. A sufficiently
large cutoff at each fixed budget makes the radical
of the remaining high tail arbitrarily small, so it
can be included in the same fixing. Each prescribed
residue is a unit, so the fixed factor can be divided
out modulo \(D\). Group the roots in each positive
window at a representative exponent of scale \(1/L\).
The discrepancy and the zero window together have
height \(X^{o(1)}\). By amplified (GW) the resulting
ratio, in this case the remaining positive integer,
has subpower normalized enumeration cost for each
residue. There are finitely many narrow coordinates
and their prime supports are disjoint by the original
valuation decomposition. If no positive atom is
retained, direct reconstruction of the small
correction under \(D\) suffices.
[Grouped short roots and wide cofactors (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

**AR3.** Under the profile hypothesis, one can replace
\((r_C+\beta)/2\) in (AR1) by
\[
                  \min\{(r_C+\beta)/2,\ d_C+\beta\}.
\]
Consequently if a positive nonatomic mass on the
marked summand permits arbitrarily cheap records,
saturation at row exponent \(1/3\) forces
\(\alpha+d_C\ge1/3\), in addition to (AR2).
This condition matters at high multiplicity when
the target has a sufficiently small nonatomic mass.
For a purely atomic target it forces the opposite
complete dictionary to reach the entire normalized
row exponent.

## Using a joint residue dictionary given by a packet list

The same estimate can be used on a fixed list
extracted by (TC3). Suppose the extraction loses
\(n^\gamma\) in the target exponent, so the common
target set on a fixed list has size at least its
original count divided by \(n^{\gamma+o(1)}\).
If the compatible residue vectors of those targets
have joint dictionary size at most \(n^\beta\), the
complete target count is at most
\[
 n^{\,\gamma+\min\{(r_C+\beta)/2,\ d_C+\beta\}+\eta}
     X^\varepsilon ,
\]
where the second term of the minimum requires the
specified target profile. Indeed the Chinese remainder
theorem identifies each compatible vector with one
unit class modulo the least common multiple. The first count uses the least common multiple, which
has polynomial height because the list length is fixed.
For the second count, project the vectors to any one
positive-weight full component and use (GW) there. The
projected list still has cardinality at most the joint
dictionary size.
This version needs no numerical surplus from the list.
If a sufficiently short real interval and large least
common multiple are available, the direct CRT bound
can be used as a further alternative. In this
comparison \(\beta\) counts compatible target residues,
not independently all combinations of residue and
height-bin labels.
[Conditional list extraction (TC3)](ABC_TARGET_PACKET_CRT_MOMENTS.md#many-targets-with-one-packet-list)
· [Short progression cells (TC4)](ABC_TARGET_PACKET_CRT_MOMENTS.md#unit-progression-cells-with-a-large-small-radical-step)

Since each fixed target radical has height-subpower
many possible complete representatives, these bounds
also apply to a bin containing the exact primitive
quotients under least repair. If the target exponent
obtained is \(v\), its \(k\)-fold scalar image with
\(\lceil n^u\rceil\) retained rows per radical has
exponent at most \(v+ku\), with the stated allowances.
For instance in the all-tall profile sector with
\(u=1/8,\ \alpha=9/50,\ r_C=1/3\), and arbitrarily
cheap packets, the row exponent from (AR1) tends
to at most
\[
                       \frac{1/3+9/50+1/8}{2}
                              =\frac{383}{1200}<\frac13.
\]
The numerical calculation is a sufficient criterion,
without any claim of arithmetic realization of the
profile. The remaining unrestricted sector includes
fully diffuse tall profiles with neither summand
dictionary compressed below \(1/3-u\). There the
single absolute-residue estimate has no strict margin,
and a gain would require additional arithmetic
information about the compatible joint residues.
[Scalar image of the true target multiplicity (SR3)](ABC_SCALAR_MASK_TARGET_MOMENTS.md#a-scalar-image-for-weighted-target-occurrence)
· [Exact least-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
