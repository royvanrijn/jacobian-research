# The small-prime splitting-curve obstruction does not close the quartic strata

For each of the four retained MW17 parents, its cubic splitting curve has
a Q_p point at every prime in the fixed panel

```
2,3,5,7,11,13,17,19,23,29,31,37,41,43.
```

All56 witnesses lie above smooth parent fibres and have **three distinct
Q_p-rational cubic roots**. The certificate independently checks168 disjoint
Hensel balls, including weighted infinity charts. Thus none of these primes
can exclude the required splitting-curve divisors merely by having an empty
local splitting curve. This is stronger than the previous discovery of
single local cubic roots.

**No global splitting-curve point or low-degree rational divisor is supplied.**
The surviving coefficient equations and the actual
[MW17 plus two-gain construction](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remain open. The test does not establish local solubility of those entire
coefficient systems, nor solubility of the splitting curve at untested places.

## 1. Why complete splitting was a separate gate

The [branch-stratum theorem](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md#5-what-remains-and-why-it-is-a-different-arithmetic-problem)
requires the cubic to split completely at every branch value where the two
abscissas select different roots. The subsequent
[reciprocity gate](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md) eliminates
the degree-one disagreement cases on Q80 by using squareclass conditions
from the actual section identities. Subsequent
[coefficient collision](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md) and
[nodal orientation](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md) proofs
close the remaining genus-zero and the genus-one k2 cases. The later
[cubic reciprocity proof](Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md)
also closes k1. The five remaining genus-one allocations require
disagreement divisors of degree4; their
global existence is unknown. This does not change any local split fibre
certified here.

A cubic root by itself does not imply complete splitting. A square
discriminant by itself does not imply a cubic root: a cyclic irreducible
cubic still has square discriminant. Earlier apparent small-prime root
obstructions also lost exceptional base and abscissa charts. The present
test first searched for local discriminant squares, then imposed actual
cubic roots, resolving the remaining exceptional cells. It never inferred
an exclusion from an unresolved cell or from a bounded failure to find one.

## 2. Exact local certificate

Each frozen choice specifies a parent, prime and integer coordinate in
either the t chart or the reciprocal chart `v=1/t`. In the latter use

```
X=v^4*x, Y=v^6*y,
A_v(v)=v^8*A(1/v), B_v(v)=v^12*B(1/v).
```

At v=0 these are the leading coefficients A_8,B_12. A further rational
abscissa scaling `X=p^s*z` gives a monic p-integral cubic
`f(z)=z^3+a*z+b`. Its discriminant is nonzero.

For each of three rational integral centres c_i, the exact checker verifies

```
u_i=v_p(f(c_i)) > 2*v_p(f'(c_i))=2*d_i,
v_p(c_i-c_j) < min(u_i-d_i,u_j-d_j)  for i!=j.       (1)
```

An exact root is allowed, with infinite u_i and nonzero derivative.
[Hensel's lemma](https://kconrad.math.uconn.edu/blurbs/gradnumthy/hensel.pdf)
gives a root within valuation radius `u_i-d_i` of each centre. The second
inequality makes these roots distinct. This proves complete splitting over
Q_p without treating a p-adic approximation as an exact root. The criterion
applies at2 and3 as well as the odd primes and does not assume good reduction
of the entire K3 surface at p. Only the selected characteristic-zero parent
fibre must be smooth.

The final input keeps the first valid witnesses from the retained local
previews. The following selected values illustrate the exceptional charts;
the four checkpoints contain the complete56-row table.

| Parent | Prime | Local base value |
|---|---:|---|
| Published R17 |2|`v=2`, so `t=1/2`|
| Direct11952 alternate Q80 |5|`v=80`, so `t=1/80`|
| Direct11952 alternate Q80 |19|`t=130359`|
| Recovered Curve302 |17|`t=597`|
| X1092 class1 |29|`t=0`|

These are local splitting witnesses at rational parameter values. Their
cubic roots are asserted in the specified completions, not over Q.

## 3. Finitely many local splitting requirements can hold simultaneously

At each witness the three roots are distinct, so the same strict Hensel
inequalities persist in an open p-adic neighbourhood of its base value.
An infinity witness can first be replaced by a finite rational value in
that neighbourhood. For any one parent, choose a common denominator D for
these finitely many rational centres. The conditions on `t=a/D` become
congruences on the integer a modulo powers of distinct primes. The Chinese
remainder theorem supplies infinitely many a satisfying all of them.

Thus, for each parent, infinitely many rational parameter values have
cubics splitting over **every tested completion simultaneously**. This
does not make the cubic split over Q or produce a global point on the
splitting curve. It explains why imposing these finite local conditions
alone cannot certify the required global branch incidence.

## 4. Retained computation and next boundary

The [input](../artifacts/generated-results/elkies-k3-common-quartic-splitting-local-gate-v1/input.json)
binds the literal four generic parents, all choices, both implementations,
and the four preview scripts and receipts. No exceptional specialization
data were used. The local previews had30 CPU-second caps; their base-cell
resolution was bounded at depth7 and2000 nodes per stage. All requested
local witnesses were found. Discovery preceded the frozen proof packet.

The [Sage producer](scripts/certify_common_quartic_local_splitting.sage)
uses40-digit p-adic roots only to select rational centres, then verifies
their exact inequalities. It used0.053 CPU seconds under a40-second/4GiB
cap. The [independent Python checker](scripts/verify_common_quartic_local_splitting.py)
uses integer/Fraction arithmetic and the original reciprocal-chart identity;
it uses no p-adic library. Its
[receipt](../artifacts/generated-results/elkies-k3-common-quartic-splitting-local-gate-v1/independent-replay.json)
records all168 root balls in0.012 seconds. The timings describe these
components, not total research time. The written Hensel and approximation
arguments are not proof-assistant verification.

```
python3 research/elkies-k3/scripts/verify_common_quartic_local_splitting.py
python3 -m unittest discover -s research/tests -p 'test_common_quartic_local_splitting.py' -q
```

The six tests reject duplicated roots, overlapping approximations to one
root, equality in the strict Hensel bound, nonintegral cubics and singular
fibres, and accept three exact roots. Two exploratory coefficient-display
errors are retained in the execution record and were not mathematical
witnesses.

The next fixed-parent gate needs global branch-divisor arithmetic or the
full coefficient incidence. Repeating the empty-local-set test on this
prime panel cannot supply it. Other primes and stronger local conditions
on the full common-quartic equations remain separate, unproved possibilities.
