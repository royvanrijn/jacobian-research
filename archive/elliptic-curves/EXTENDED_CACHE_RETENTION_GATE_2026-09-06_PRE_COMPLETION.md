# Extended-score retention and the11952 trace cache

A retrospective control comparison exposes substantial reordering before the
stronger score is available. Against the **same967 saved outer11952 candidates**,
the known29 control at89074/31895 has hypothetical insertion position435 under
S1 through4093 and2 under extended S1 through32749. The exact score call,
four independent character sums and both orderings replay. Its denominator
residue150 modulo1024 differs from the sampled positive residue754.

The [comparison certificate](../../artifacts/generated-results/elliptic-curves/outer_known29_retention_comparison_v1.json)
is conditional on the saved population. It proves neither a full-population
quantile nor a discard probability or prospective success rate. The completed
[outer48 experiment](OUTER131072_TRIAL_2026-09-06.md) remains unchanged:2160 boxes,
2223 retained point witnesses, all48 certified bounds still17 modulo2,3,5.

## Fixed cache cost gate

Calling PARI independently for each retained rational parameter repeats many
finite-field fibres. For a selection prime p, a projective table has onlyp+1
entries. Homogeneity of degrees8and12 makes the raw fibre at a primitive(n,d)
isomorphic over F_p to the table entry n/d when d is a unit, and to the
leading-coefficient entry when d=0 modulo p. This preserves the existing
raw-good-reduction score policy, including its singular markers.

An eight-prime prototype on11952 spans4099through32749. All projective residue,
discriminant and Hasse frames pass, together with40 independent full character
sums. Its complete run takes2.830643 seconds; replay takes0.431832 seconds.
The largest measured wall-time-per-residue projects to1619.190 seconds for
all selection primes on one worker, passing the declared1800-second gate.
This extrapolation includes startup overhead and is **not a runtime bound**.

The [prototype certificate](../../artifacts/generated-results/elliptic-curves/extended_projective_trace_cache_benchmark_v2.json)
retains all eight cases. Version one returned no table because a GP function
definition consumed the following execution loop on the same line. That failed
call is preserved. Version two separates the statements; a small signed-trace
syntax regression passes, and the eight primes and budgets remain unchanged.

## One complete cache and one score benchmark

The active cache protocol fixes **2948 primes and52989620 projective residues**
for11952 alone. It reuses the eight benchmark tables and allows one20-second GP
call for each other prime, two workers, a3600-second build and1200-second exact
replay. Every table receives complete frame checks and five independent
character sums. Checkpoints close each sixteen-prime block; any failure stops
further blocks, with no retries or broader prime range.

Only after that full replay passes does the separate cached-score controller
encode a little-endian integer cache, verify every encoded byte against the
source tables, compile the retained-list scorer and compare its result on all967
previously scalar-scored outer11952 candidates. Quantization preserves the exact
Python operation order of the existing extended-score implementation. Validation
primes never enter this cache. Signed-residue and infinity arithmetic and
fail-closed truncated-cache handling pass separate compiled regressions.

This controller is a tool and cost check. It launches no parameter enumeration,
new finalist selection or point search. A larger survivor population still needs
its own frozen selection, budget and point-proof protocol. The source scripts
are `benchmark_extended_projective_trace_cache_v2.py`,
`build_extended_projective_trace_cache_11952.py`,
`benchmark_retained_extended_cache.py` and `finish_extended_cache_benchmark.py`
under `../cas/`. Their immutable protocols and logs are in the corresponding
local elliptic artifact directories.
