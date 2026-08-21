# Two split-infinity Mestre rank-19 frontiers

## Exact status

Two bounded specialization searches in the split-infinity six-root Mestre
families produced exact rank-at-least-19 curves below the operational
`log(N)<182.72` cutoff.

| family roots | `u` | `T` | minimal model | `log(N)` |
| --- | ---: | ---: | --- | ---: |
| `(0,25,95,143,168,205)` | `483` | `-8441/42` | `[1,0,0,-639860186233795396843608642357,185944471483117795919620263190737781477557921]` | `157.759935999987...` |
| `(0,43,128,197,231,289)` | `660` | `-12655/44` | `[1,-1,1,-14853112884170608407469777362587,21717789706067374061109415673917936898882673611]` | `164.053646218834...` |

For each curve, exact rational quartic points are mapped to its short
Jacobian.  Stacked good-reduction images in `E(F_p)/3E(F_p)` have full column
rank 19, and a separate good reduction excludes rational 3-torsion.  Infinite
descent therefore proves the nineteen selected rational points independent.
PARI/GP directly replays each minimal model, conductor, and root number.  Both
root numbers are `-1`.

These are rank lower bounds, not unconditional exact-rank statements.  They
are two points short of the rank-21 target and do not change the tracked list
of rank-at-least-20 scores.

## Conditional analytic diagnostic

The audited sinc-squared explicit formula at `Delta=11/5` gives conservative
values

```text
family2_u483  20.3010039962309377161...
family3_u660  20.4139420468395791917...
```

Under GRH these bound analytic rank.  Since both functional-equation signs
are odd, each analytic rank is then at most 19.  Together with BSD this would
identify algebraic rank exactly 19.  Neither GRH nor BSD is used in the
unconditional rank lower bounds.

## Replay

```sh
PYTHONPATH=elliptic-curves/cas .venv/bin/python \
  elliptic-curves/cas/certify_mestre_dsquare_rank19_frontiers.py --check
```

The pinned artifact is
[`elliptic_mestre_dsquare_rank19_frontiers.json`](../../artifacts/generated-results/elliptic_mestre_dsquare_rank19_frontiers.json),
with file SHA-256
`e78613cc35ad523242a6d3af529a4b59bece136a8f2b7880bead3ef7094144be`.
The script stores the exact points, finite-reduction matrices, direct global
data, prime sums, and the distinction between unconditional and conditional
claims.

## Search scope

Discovery used finite `ratpoints` boxes and exact post-search independence
tests.  Alternate-cover passes on the fixed rank-18/19 fibers returned only
points in the certified subgroups.  The negative searches are bounded
experiments and are not rank upper bounds; the GRH-conditional explicit
formula, rather than those misses, is the reason to deprioritize these fixed
fibers.
