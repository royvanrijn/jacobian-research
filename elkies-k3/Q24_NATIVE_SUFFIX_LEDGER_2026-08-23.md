# H3 q24-native D12 -> A1 suffix ledger

Date: 2026-08-23

Status: **exact lattice/root-Weyl replay from the q24 equation-side D12 frame**.

This note records the result introduced by commit
`bfb96bf526b0b289e12d531e7aff697cae5f39c6` (`lifting q24 path`) and keeps its proof
boundary separate from characteristic-zero equation execution.

## Source frame

The suffix starts from `D24eq`, exported by

```text
scripts/export_h92_q24_native_d12_frame.sage
```

from the selected q24/D13 equation-side geometry. The exporter deliberately does not
identify this frame with q32 or with a later historical D12 frame solely from common
`D12/MW5` root data.

This distinction matters: an ADE/MW label does not determine the embedded elliptic `U`,
zero section, or the next neighbour marking.

## Replay

Primary driver:

```text
scripts/run_h92_q24_native_suffix_to_a1.py
```

Supporting search/replay scripts added in the same commit include

```text
scripts/search_h92_q24_native_q4_candidates.py
scripts/search_h92_q24_native_selected_q.sage
scripts/search_h92_q24_native_historical_weyl.sage
scripts/search_h92_q24_native_historical_weyl_exact.py
scripts/test_h92_q24_native_candidate.sage
```

The selected q-sequence downstream of D12 is

```text
6, 8, 4, 4, 4, 4, 4, 4, 4
```

and the replayed children are

| step | q | child | root rank | MW rank |
|---:|---:|---|---:|---:|
| 1 | 6 | `A11` | 11 | 6 |
| 2 | 8 | `2A5` | 10 | 7 |
| 3 | 4 | `3A3` | 9 | 8 |
| 4 | 4 | `A3+2A2` | 7 | 10 |
| 5 | 4 | `5A1` | 5 | 12 |
| 6 | 4 | `4A1` | 4 | 13 |
| 7 | 4 | `3A1` | 3 | 14 |
| 8 | 4 | `2A1` | 2 | 15 |
| 9 | 4 | `A1` | 1 | 16 |

The reported final step is

```text
ROOTWEYL_HIST|q=4|root_data=1,2,2|MW=16|orbits=1
Q24NATIVE_SUFFIX|step=9|q=4|orbit=17593|ADE=A1|root_data=1,2,2|MW=16|status=PASS
Q24NATIVE_SUFFIX_RESULT|steps=9|final=A1/MW16|status=PASS_Q24_NATIVE_D12_TO_A1
```

Local output artifact:

```text
artifacts/local/elkies-k3/q24-native-suffix/q24-native-d12-to-a1.json
```

## What this certifies

The result certifies that the selected q24 equation-side D12 marking has a concrete
native lattice/root-Weyl continuation through the same desired ADE/MW corridor to
`A1/MW16`.

It therefore removes the need to use either

- a q32-derived D12 model, or
- a historical D12 frame matched only by lattice invariants

as the preferred bridge from the q24 front end to the known high-rank suffix.

The already-certified final q6 neighbour from `A1/MW16` to the rootless `MW17` frame can
then be appended at lattice level.

## What this does not certify

This replay is downstream lattice/root-Weyl work. It does **not** by itself provide
characteristic-zero Weierstrass equations, transported zero sections, or RR pencils for
the nine arrows from D12 to A1.

The immediate equation boundary remains:

```text
exact D13/MW4 equation
  -- selected q24 orbit 85 -->
native D12/MW5 equation + marking
```

Once that native D12 equation certificate is complete, the correct equation programme is
to execute the stored q24-native suffix in order while preserving the marking at every
stage.

## Preferred-route consequence

The preferred H3 path is therefore the q24-native path documented in
[`H3_PREFERRED_PATH.md`](H3_PREFERRED_PATH.md). q32 remains useful as a compiler and
marking-recovery experiment, but is no longer the preferred continuation while the native
q24 route remains viable.

## Reproduction

From a checkout containing the prerequisites, the intended top-level replay is

```bash
python3 elkies-k3/scripts/run_h92_q24_native_suffix_to_a1.py
```

Expected terminal status:

```text
PASS_Q24_NATIVE_D12_TO_A1
```
