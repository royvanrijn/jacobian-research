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

The complete cache covers **2948 primes and52989620 projective residues**
for11952 alone. Build and replay pass in556.278 and93.302 seconds under the
declared3600/1200-second limits. All14740 independent character sums pass.
The eight benchmark tables are reused; the remaining tables each have one
bounded GP call. All raw tables, frames and supervisor records are retained.

The cached-score controller also passes: all967 previously scalar-scored
outer11952 candidates have exactly matching sums and good-prime counts. The
476930184-byte integer cache is checked against every source entry. Scoring
takes0.483 seconds. A separate short-band cache exactly matches all967 original
short scores. Quantization preserves each existing implementation's arithmetic
order; validation primes never enter the selection cache. Signed-residue and
infinity arithmetic and truncated-cache rejection pass compiled regressions.

A full-size fixture repeats these967 inputs cyclically to1048576 rows. Every
output is checked: extended scoring takes11.408 seconds and short scoring2.230
seconds, both passing their120-second gates. This repeated fixture measures
cost and agreement; it supplies no new candidates or rank information.

The separately frozen [full11952 trial](FULL11952_RETENTION_TRIAL_2026-09-06.md)
now completes the wider short population and all1048576 retained extended
scores. Its64 finalists have fresh scalar checks; point exposure has its own
fixed protocol. This cache gate alone is a tool and cost check. The source scripts
are `benchmark_extended_projective_trace_cache_v2.py`,
`build_extended_projective_trace_cache_11952.py`,
`benchmark_retained_extended_cache.py` and `finish_extended_cache_benchmark.py`
under `../cas/`. Their immutable protocols and logs are in the corresponding
local elliptic artifact directories.

The [pre-completion note](../../archive/elliptic-curves/EXTENDED_CACHE_RETENTION_GATE_2026-09-06_PRE_COMPLETION.md)
is historical. No isolated portable replay of the complete cache is claimed.
