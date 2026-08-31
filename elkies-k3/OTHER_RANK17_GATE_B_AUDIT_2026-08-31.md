# Other rank-17 fibration: Gate B audit — 2026-08-31

## Result

The target-specific short-bridge probe found no usable crossover from the
equation-explicit H3 corridor to either the alternate `A1/MW16` parent or its
rootless `MW17` child.

The cheapest literal graph-ray intersections are already too large:

| target | cheapest equation node | old-fibre degree | presentation `q` |
| --- | --- | ---: | ---: |
| alternate `A1/MW16` | current rootless | 391 | 5,660,116 |
| alternate rootless `MW17` | current rootless | 11,511 | 4,905,228,474 |

These are exact integral-lattice costs in the fixed pinned marking.  They are
not substitutes for resolved Riemann--Roch dimensions.  In particular, the
alternate q4 equation divisor is presently pinned geometrically by its
degree-47 CM24 specialization, while the final q6 has its separate certified
two-reflection reduction in the alternate A1 chamber.

## Meet-in-the-middle probe

The bounded two-sided search covered all fourteen equation-explicit nodes at
and after the equation-D13 common marking, including the physical q323-free
branch.  It searched both alternate endpoints with

```text
old-fibre degrees 1,2,3,4 on both sides,
q in {4,6,8,12} on both sides,
41 target-ray scales from 20% through 180%,
64 exact closest vectors per scale.
```

There were 955,136 closest-vector visits and 6,998 distinct exact fixed-norm
presentations.  None was even a two-sided lattice match, before the exact
component, closest-section, and augmented-lattice horizontal-nef gates became
necessary.  Thus the run found no physically nef bridge.

This is a target-directed bounded negative computation, not an exhaustive
nonexistence theorem for every low-degree ray.  A naive complete short-vector
shell was also attempted and rejected as the route to a certificate: the
first low-root-rank node exceeded a 1 GiB PARI stack, and the complete
Mordell--Weil-quotient shell remained too large.  No capped enumeration is
presented as exhaustive.

The reproducible artifacts are

```text
artifacts/generated-results/elkies-k3-other-r17-gate-b-direct-costs.json
artifacts/generated-results/elkies-k3-other-r17-gate-b-mitm.json
```

and are produced by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_other_r17_gate_b_direct_costs.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_other_r17_gate_b_mitm.sage
```

## Q80 fallback boundary

Gate B therefore selects the Q80 reconstruction route.  The repository does
contain an exact characteristic-zero **CM24 specialization** through the
terminal q6, with final child `4A2+A3+A5/MW2`; it must not be mistaken for the
determinant-948 generic rootless `MW17` equation required here.

The generic alternate suffix is currently reproducible only as an exact
lattice route plus a `GF(73)` equation shadow.  Its finite-field equation
scripts hard-code both the prime and prime-specific section coordinates.  The
intermediate generated inputs beginning with

```text
q80-fourth-q12-cm24-moving-cubic-gf73.json
```

are also absent from this checkout.  Consequently the present code cannot be
rerun at a second prime by changing a command-line option, and one residue
cannot support CRT/LLL reconstruction.

The next construction step is therefore earlier than CRT: recover and persist
the generic fourth-q12 equation/marking in a prime-independent form, then
parameterize the compensated pair23 and final-q6 compilers by the prime.  Only
after compatible outputs exist at several split good primes is canonical
marking alignment and rational reconstruction justified.  Until that generic
equation is obtained, no alternate-frame `j`-map exists and no recognition
claim for the rank-29 curve or ICARM 398--400, 273, or 302 is made.
