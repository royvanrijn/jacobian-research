# R17 higher-$2$-power Selmer boundary

This is a conditional higher-descent protocol for record fibres 356 and 385,
not a completed computation or a claim in
[`MATH_STATUS.json`](../../MATH_STATUS.json).  Their complete $2$-Selmer
groups and every higher-filtration dimension remain `UNKNOWN`.  The current
global prerequisite is the [relative Selmer filter](ELKIES_R17_RELATIVE_2SELMER_PIPELINE.md).

Let $H$ be the certified generic rank-17 Kummer image and let $W$ be the
known twelve-dimensional residual rational-point block.  Once complete descents
exist, use the image filtration

\[
 F_j=\operatorname{im}(\operatorname{Sel}_{2^j}(E/\mathbb Q)
 \longrightarrow\operatorname{Sel}_2(E/\mathbb Q))/H,\qquad j=1,2,3,
\]

with $F_1\supseteq F_2\supseteq F_3\supseteq W$.  This is not a literal
inclusion of groups named $S_2,S_4,S_8$.  Rational point classes already lift
at every level; the computation concerns only the complement to $W$.  The first
Cassels--Tate drop has even codimension.  If a certified stage equals $W$, stop:
the rank upper bound meets the known rank-29 lower bound.  An image larger than
$W$ leaves room but does not construct a rational point.

After a complete global $F_1$ and all-place matrix are certified, freeze their
basis before opening the exceptional-point labels.  Compute pairing and cover
data only on a complementary quotient of $F_1/W$, retain its full linear
structure, and compare the two curves through basis-independent dimensions and
local contribution profiles.  A failed basis cover, incomplete BNF, capped
higher-descent call, raw MW17 local fingerprint, or bounded point-search miss
is `UNKNOWN`, not a Selmer or rank bound.

The complete protocol, precise decisive profiles, literature boundary and
historical command surface are preserved byte-for-byte in the
[archive](../../archive/elliptic-curves/notes/R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md.txt).
Do not launch it merely because the protocol exists; a future run needs a
certified global envelope, complete local conditions, immutable inputs and a
separately scoped certificate plan.
