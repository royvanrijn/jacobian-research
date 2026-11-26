# Empirical prices of packets supported on the nonatomic part

At the minimum nonatomic mass, a packet taken from diffuse valuation
coordinates has its full radical size as an empirical price. The
fixed-letter nonatomic bound gives this directly on original fibres
and retains the radical caps on all three roles. In particular, the
one-letter packet price obstruction persists for these packets in
actual near-full additive graphs, if such graphs exist.

## A price bound from one original fibre

Use \(c_0=1/6\), \(n=X^L\), and the height-first small-budget
limit. Consider a sequence of actual primitive additive graphs
\(a+b=c\), \(X/2<c\le X\), with all three complete dictionary
exponents at most \(c_0+o(1)\) and all three individual radical
caps \(n^{2c_0+o(1)}\). Consider a cohort of edge exponent
\(2c_0-o(1)\) with a designated packet \(M\) of full prime powers
of one original role \(T\). The packet can vary in the cohort.
Suppose its dictionary of distinct moduli has exponent at most
\(s+o(1)\). Assume that the complete scaled valuation measures
on \(T\) on the cohort converge to the same profile \(\nu_T\),
as in anatomic height-box refinements. Write
\(c_T=(\nu_T)_{\mathrm{na}}([0,\infty))\).

Pin one complete value of a different role \(U\), and pin one
packet \(M\), retaining a fixed-\(U\) original fibre of exponent
at least \(c_0-s-o(1)\). Such a choice follows by averaging the
cohort over at most \(n^{c_0+s+o(1)}\) pairs of records. If
\(s<c_0\), pass to a subsequence on which the packet measures
\[
          \mu_M=\sum_{p^e\parallel M}\log_n p\,\delta_{Le}
\]
have a weak limit \(\mu\), on this pinned choice. Their first
moments are bounded by \(1\), since \(M\le X\), so the measures
are tight. Their masses are bounded by the original target
cap. Domination by the exact valuation measures on \(T\)
gives \(\mu\le\nu_T\). Put
\(b=\mu_{\mathrm{na}}([0,\infty))\).

**RP1.** For every such pinned subsequential limit,
\[
                           s\ge c_0-c_T+b .                  \tag{RP1}
\]
For the proof use the prime support of \(M\) as the prescribed
pool on \(T\) in the original fixed-\(U\) fibre. Its product
is at most \(X\). The packet contains the full powers on its
support, so after stripping the pool the residual \(T\)-measure
is exactly the difference of the two measures before taking
limits. The residual limit is \(\nu_T-\mu\). The nonatomic
part has mass \(c_T-b\), by domination and decomposing the two
measures into their atomic and nonatomic parts. Apply the
fixed-fibre bound (FT1) to this residual \(T\)-measure, which
gives \(c_0-s\le c_T-b\). If \(s\ge c_0\) the conclusion
follows from \(b\le c_T\) without pinning a positive-exponent
fibre.
[Prescribed-pool nonatomic mass on a fixed fibre (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

For a near-full empirical bin in (AP2), the packet dictionary
has exponent at most the price \(s\) plus the prescribed cap
tolerance. Thus (RP1) applies to this price when the tolerances
tend to zero. Together with the separate fixed-letter radical
cap argument there, a bin of packet radical mass approaching
\(\beta\) satisfies
\[
                  s\ge\max\{\beta/2,\ c_0-c_T+b\}.
                                                               \tag{RP2}
\]
The same proof applies to the price of each marked role in a
simultaneously binned near-full family. It does not require
an independence assumption on the marked roles, since one
can pin just the relevant packet and one different complete
letter on the same retained edge cohort.
[Individual packet price from fixed-fibre radical caps (AP2)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#composite-packet-inverse)
· [Simultaneous packet bins (SB2)](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md#simultaneous-marking-theorem)

The same averaging gives the exponent-dependent form
\[
                       s\ge D-\alpha-c_T+b ,
\]
where the edge cohort has limiting exponent \(D\), the chosen
opposite complete dictionary has exponent at most \(\alpha\),
and all three original radical caps are bounded in base \(n\).
Indeed the pinned fibre has exponent at least \(D-\alpha-s\).
Whenever that number is positive, (FT1) applies directly.
Otherwise the inequality follows from \(b\le c_T\).

Several packet records on the same target may be treated as one
by taking their least common multiple, still using full prime
powers on that target. Thus (RP1) also lower-bounds the joint
price by the nonatomic mass of their union at the minimum
endpoint. The records may come from different bands or may
overlap. The natural cost of a shared prime is counted once.


## Consequences at the minimum mass

At total nonatomic mass \(1/2\) one has \(c_T=c_0\) on every
role. A purely nonatomic packet of radical mass \(\beta\)
then has \(b=\beta\), so \(s\ge\beta\). The usual support
enumeration gives the opposite limiting inequality
\(s\le\beta\) when all packet supports have sizes
\(n^{\beta+o(1)}\). Thus the empirical packet price of such a
bin equals its natural prime cost.

A precise way to guarantee the hypothesis is to take the
packet from a fixed positive compact continuity region of
the scaled degrees after excluding neighborhoods of the atoms.
If the closure of the retained region has \(\nu_T\)-atomic
mass at most \(\epsilon\), any resulting packet of mass
\(\beta\) in it has \(b\ge\beta-\epsilon\), by domination.
One can let \(\epsilon\) tend to zero in a diagonal sequence
of such refinements. For every fixed forbidden price margin
this uses fixed profile tolerances first, sufficiently small
fixed budgets next, and arbitrarily large heights at each
such budget. Nothing requires packets to be large individual
prime powers.

In the same minimum-mass sector, the nonatomic support lies
in \([0,6]\). For purely nonatomic packets in compact degree
ranges, their numerical weight \(w\) in the limit therefore
satisfies \(w\le6\beta=6s\). The conditions
\(kw>1\) and \(ks<1/6\) in the direct-height version of
(AP4) cannot both hold. Thus discounting a purely diffuse
packet alone cannot make that height test eliminate the
minimum-mass endpoint. More generally \(s\) near \(\beta/2\)
at \(c_T=c_0\) requires \(b\le \beta/2+o(1)\). At least
about half of such a packet's radical mass must come from
the atomic part. Away from the minimum nonatomic mass the
bound loses exactly the term \(c_T-c_0\) in (RP1).
[Upper support of the minimum nonatomic degrees (MS11)](ABC_DIFFUSE_FACTOR_GRAPH.md#consequence-for-all-three-nonatomic-supports)
· [Height exclusion with empirical packet price (AP4)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#height-exclusion-and-power-classes)
