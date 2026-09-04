# Residual-Selmer fingerprints and the monotone search sieve

## Result and boundary

The complete groups

\[
  \operatorname{Sel}_2(E_t/\mathbf Q)/
  \operatorname{im}(MW17/2MW17)
\]

for ICARM curves 356 and 385 remain `UNKNOWN`.  The exact cubic BNF and full
local-solubility calculation have not completed, so neither a Selmer upper
bound nor an exact-rank-29 claim is made.

What is now complete is the richer control layer needed around that descent:

- the two rank-29 fibres are exact inputs to the checkpointed complete
  2-descent;
- the twelve certified point classes modulo generic MW17 are embedded as the
  displayed residual lower-bound space;
- the two rigid-character directions are recorded inside that twelve-space;
- every exact control retains its known residual local Kummer subspaces,
  cumulative intersections, delete-one-place ranks, component data, available
  Hilbert information, and sampled `t mod p^k` strata;
- the complete descent worker now emits the analogous delete-one-place ranks
  on the *full* global norm-square ambient space when the BNF is available;
- an incomplete monotone residual sieve may authorize only an explicitly
  bounded point search.  It cannot authorize a Selmer or rank claim.

The machine-readable record is
[`../artifacts/generated-results/elkies-k3-r17-residual-selmer-fingerprints-v1.json`](../artifacts/generated-results/elkies-k3-r17-residual-selmer-fingerprints-v1.json).

## Exact known residual and rigid quotients

For both record fibres, the first seventeen displayed points are the exact
specialized generic MW17 subgroup and `P18,...,P29` are independent modulo it.
Thus the certified realized residual subgroup is `F_2^12`.

For curve 356 the rigid plane has rows

```text
000000010000
011010010100
```

in `P18,...,P29` coordinates.  Its leftmost-pivot complement is

```text
P18,P20,P21,P22,P23,P24,P26,P27,P28,P29.
```

For curve 385 the rigid rows are

```text
000100000000
000100110000
```

and the complement is

```text
P18,P19,P20,P22,P23,P25,P26,P27,P28,P29.
```

Consequently the known subgroup modulo the rigid plane is exactly
ten-dimensional.  This is a lower-bound statement inside the residual Selmer
quotient.  Once a complete descent returns total 2-Selmer dimension `d`, the
three decisive dimensions are

```text
residual modulo MW17              d - 17
after quotienting the rigid plane d - 19
additional beyond all 29 points   d - 29.
```

If `d=29`, the known twelve classes exhaust the residual Selmer group and the
rigid quotient is exactly the displayed ten-space.  Together with the already
certified trivial rational 2-torsion and saturation/independence gates, this
would supply the desired rank upper bound 29.  No such value of `d` is asserted
here.

## Control fingerprint: `+12` versus `+5`

The exact control dimensions are

```text
curve     351  356  376  377  385   12
gain        8   12    5    6   12   12
```

At the eleven places audited for both record fibres and the `+5` control 376,
the following scalar features separate *both* `+12` samples from 376:

| place | separating exact features |
|---:|---|
| 2 | Kodaira symbol, minimal-discriminant valuation, Tamagawa 2-part, known-residual localization-kernel dimension |
| 13 | known-residual localization-kernel dimension |
| 37 | localization-kernel dimension, nontrivial component-Hilbert pair count |
| 53, 67, 71, 83, 113 | localization-kernel dimension |
| 79, 97, 101 | ambient local Kummer dimension, known image dimension, kernel dimension, selected-block image dimension |

This confirms that 2 is structurally discriminating in the current 356
comparison, but not that it alone creates the jump.  The odd-place pattern is
substantial, including for 385.  With two high samples and one low sample this
is a theorem about the controls, not a sufficient local criterion or a
statistical classifier.

The stacked known-point localization matrix has full source rank on all six
controls.  Deleting any one audited place leaves that rank unchanged.  This
means every individual place is redundant *for separating the certified known
classes in the full stored place set*.  It does not answer which place cuts
down the unknown global Selmer ambient space.  That requested matrix is now
computed by `run_elkies_2026_relative_2selmer_checkpointed.py` after a
certified BNF; until then its ranks remain `UNKNOWN`.

## Pairing and decomposition boundary

The available componentwise Hilbert forms do not descend through the rigid
plane: the exact obstruction primes are `13,23,37,139` for curve 356 and
`5,29,37,41,73,109,127` for curve 385.  Their corestricted local Tate controls
are zero.  Therefore there is no canonical pairing on the ten-dimensional
quotient from these data and no invariant meaning to an “indecomposable ten”
or a split into smaller blocks.  The certificate records this as `NOT_DEFINED`
instead of diagonalizing an arbitrary coordinate complement.

A complete Cassels pairing on the complete residual Selmer group would make
the decomposition question meaningful and is retained as a separate field in
the fingerprint schema when available.

## Parameter strata and CRT prototypes

For every common audited prime the artifact records the exact affine or
infinity-chart residue of each control modulo `p`, `p^2`, and `p^3`, together
with the complete local-feature hash.  These are finite sampled strata; the
data do not prove that a fingerprint is constant on an entire residue
cylinder.

Two exact CRT classes preserve five discriminating `p^3` residues:

```text
356 prototype:
t = 14503794234702288112 + 47438163879590960216*n
places 2,13,37,53,71

385 prototype:
t = 13329277794157146704 + 39863665779809550328*n
places 2,13,37,53,67
```

The congruences are exact; preservation of the local Kummer/Selmer fingerprint
is a search hypothesis.  Every manufactured parameter must be re-audited
before it inherits any arithmetic label.

## Monotone residual sieve

`elliptic-curves/cas/elkies_residual_selmer_gate.py` now separates two kinds of
authorization.

1. A complete unconditional descent may issue the existing exact residual
   gate and support theorem-directed cover/search work.
2. An incomplete sieve stores a sequence of proved residual upper bounds.
   Bounds must be nonincreasing and carry evidence provenance.  A missing BNF
   is stored as “no finite upper bound yet”, not as a numerical estimate.  The
   fibre is rejected as soon as the proved bound is below 15.  Otherwise only
   a point search with explicit height/time/resource limits is authorized.

The point-search entrypoints report their requested limits back to the gate;
an open authorization is rejected if a required limit is absent or the run
would exceed its recorded allowance.

An open sieve always stores `theorem_claim_authorized=false` and
`expensive_search_authorized=false`.  A bounded search can improve the
Mordell--Weil lower bound; it cannot turn missing global descent data into an
upper bound.

## Replay

```bash
sage -python \
  elkies-k3/scripts/build_r17_residual_selmer_fingerprints.sage --check

python3 -m unittest -v \
  elliptic-curves/tests/test_elkies_residual_selmer_gate.py \
  elliptic-curves/tests/test_elkies_relative_2selmer_checkpointed.py
```

The record-pair complete-descent input suite can be regenerated with

```bash
python3 elliptic-curves/cas/build_elkies_2026_relative_2selmer_suite.py \
  --record-pair-only \
  --output-dir \
  artifacts/local/elliptic-curves/r17-074d9-record-residual-2selmer-v1/programs \
  --manifest \
  artifacts/local/elliptic-curves/r17-074d9-record-residual-2selmer-v1/manifest/input.json \
  --overwrite
```

If the PARI 2.19 relation collector produces a `bnfcertify`-validated binary
checkpoint, continue in that same GP build with
`elliptic-curves/cas/run_elkies_2026_pari219_selmer_from_bnf.py`. This avoids
binary-checkpoint compatibility assumptions and emits the requested full
leave-one-place-out matrix before any known-point alignment.

Long BNF/descent outputs remain under `artifacts/local/elliptic-curves/`; only
a completed certified result should be promoted to the generated-results
tree.
