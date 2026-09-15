# The minimum-distance inherited contact layer for trace P3

**Bounded exact construction attempt; the correlated-gain objective remains
open.** For the direct11952 trace T=P3 (zero-based), none of the 230 generic
sections R satisfying h(R)=4 and h(T-2R)=14 supplies a rational contact
2R(t)=T(t). This includes the trace's coordinate pole and infinity. It does
not exclude other traces, larger distance heights, higher-height inherited
sections, or halves outside the generic section image.

The input is selected by geometry, not by a specialization target. The
[first singular-bisection constructor](SINGULAR_BISECTION_TRACE_GEOMETRY_2026-09-15.md)
needs a rational point of the genus-nine halving curve. The
[inherited-contact intersection bound](Q80_RATIONAL_CONTACT_P3_ATTEMPT_2026-09-15.md)
is h(T-2R)>=14. Thus the present bank is the minimum inherited height on this
rootless K3, at the minimum eligible distance. It extends the earlier signed
basis attempt to the entire retained height-four shell at this distance.
It is not a bound on the height of an unknown half.

## Exact calculation

The retained [height-four roster](../artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/norm4-sections.json)
contains 1,313 sections up to sign, with their generic lattice words. Its
coordinates are only over F131. The new calculation uses those words and
the actual rational generic basis; it does not lift finite-field coordinates
by guessing rational coefficients. The shell's completeness is inherited
from the [retained complete replay](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md).

For each orientation the exact height matrix tests h(T-2R)=14. There are
230 eligible words. At seven smooth rational fibres, away from all basis
poles, exact rational group addition evaluates each word. Polynomial
interpolation gives x of degree at most four and y of degree at most six.
Literal substitution verifies y²=x³+Ax+B over Q(t).

These points are attached to the words, rather than merely to the curve:
two distinct height-four sections on this rootless K3 intersect at most
six times. Indeed their pairing is at least -4 by positive definiteness,
and their intersection number is 2 minus that pairing. Agreement at seven
distinct fibres therefore proves equality. This uses the retained geometric
height formula and the full marked generic basis.

For every reconstructed R, cancel the poles from x(T)-x(2R). Its numerator
is a product of a quintic and a degree-thirteen polynomial, each with
multiplicity one. The checker verifies the product identity and supplies,
for each factor, a prime at most 997 where its coefficients are integral,
its leading coefficient is a unit, and it has no residue root. Consequently
neither factor has a rational root. This argument does not depend on an
irreducibility label returned by the factorization routine.

No ordinate of R vanishes at the trace pole, so 2R cannot equal T=O there.
The weighted leading abscissas at infinity are unequal in all 230 cases.
Thus the affine cancellation loses no rational contact. Checking abscissas
also covers the larger-distance opposite sign; no claim of completeness for
that other distance layer is needed here.

## Reproduction and next boundary

The [checker](scripts/verify_q80_minimal_inherited_contact.py) uses python-flint,
with a 60 CPU-second and 1 GiB limit; production and replay finish in a few
seconds. The compact [result](../artifacts/generated-results/elkies-k3-q80-minimal-inherited-contact-v1/result.json)
retains the source and roster hashes, selected words, interpolation sites,
coordinate and factor hashes, and all 460 modular root-obstruction primes.
Every polynomial is recomputed from the retained rational input during replay.
The geometry and inherited shell completeness are not independently proved
by this script; no formal verification or external review is claimed.

```
.venv/bin/python research/elkies-k3/scripts/verify_q80_minimal_inherited_contact.py
```

No tangent or residual quartic can be constructed from this bank. A new
attempt needs an actual half from outside this completed layer before doing
quotient interpolation. Neither this result nor the previous basis attempt
provides an all-height obstruction. Matching a literal low-genus cover with
a second independent section and proving infinitely many rational base
points remain requirements of the
[original objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md).
