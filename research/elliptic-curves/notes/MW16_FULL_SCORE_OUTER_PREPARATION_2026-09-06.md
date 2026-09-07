# Preparing longer-score selection in new MW16 territory

**The balanced fifteen-table benchmark and independent arithmetic replay pass.
All 14,740 five-family tables are complete; full replay, binary encoding and
every saved-score comparison pass. No new parameter or point search runs under this cache protocol.**

The fresh outer trial supplies three new high-rank fibres beyond its former
signed short-score prefix of 512, including the 26-point curve at position
528. The separate 11952 benchmark finds 366 of 512 longer-score leaders
outside an earlier 4,096-candidate short prefix. This motivates scoring before
truncation on untouched higher MW16 parameter bands. It does not estimate
rank density or prove a better rank predictor.

A projective trace table contains every residue t in P1(Fp) of the fixed
homogeneous model. It can be used at arbitrary rational parameter heights;
there is no need to repeat the same finite-field calculation for each new
rational parameter. Singular raw-model residues remain marked bad.

The same three predetermined primes, 4099, 17749 and 32749, are tested on all
five existing MW16 atlas models. Every residue, infinity value, singularity
marker and Hasse bound is checked. All 75 independent character sums pass.
The descriptive serial-time projections are 1,613–1,659 seconds per family,
with all five below the frozen 1,800-second gate. Actual finite time limits
are enforced separately. No record equations, parameters, points, ranks,
j-invariants or jump labels enter these calculations.

The cache protocol fixes all 2,948 primes from 4099 through 32749 for each
of the five families: 14,740 tables, of which 15 reuse the frozen benchmark.
It permits at most 14,725 new PARI calls, five workers, 20 seconds per call,
80-case checkpoints and 7,200 seconds for the build. Each table has five
independent character-sum checks. The subsequent whole-cache replay has a
3,600-second cap and runs no fresh PARI calls. Failure or censoring stops
without retry. Later parameter slices and point exposure require separate
frozen protocols after the cache replay passes.

Sources: `../cas/benchmark_mw16_extended_projective_tables.py` and
`../cas/build_mw16_extended_projective_caches.py`. The benchmark certificate is
`mw16_extended_projective_benchmark_v1.json`; protocols and raw tables are
under `mw16-extended-projective-benchmark-v1` and
`mw16-extended-projective-caches-v1` in local artifacts. This computation
supplies reusable finite-field values, not new rational points or rank bounds.

The [next higher-annulus protocol](MW16_JOINT_HIGHER_ANNULI_2026-09-06.md)
now freezes twenty untouched slices and their subsequent scalar/point budgets.
Its controller waits for every cache and saved-score proof above.

The complete build covers 264,948,100 projective residues and passes 73,700
independent character sums, retaining all raw calls and the 15 reused tables.
The whole-cache replay and binary-cache comparisons both pass. All
2,432,852,660 encoded bytes replay, and both score components on all 40,960
saved target-free candidates match their original canonical/scalar proofs.
