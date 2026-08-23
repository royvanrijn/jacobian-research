# H3 preferred construction path

Status: **preferred route and proof-boundary note**.

This file records the route that should be treated as the default H3 construction path after commit `bfb96bf526b0b289e12d531e7aff697cae5f39c6` (`lifting q24 path`). It is intentionally narrower than `CONSTRUCTION_ROUTES.md`: that file retains comparison routes and discovery provenance, while this file answers only **which path should we compile next?**

## Preferred route

```text
H3 E7+E8/MW2
 --q6 --> E8+E6/MW3
 --q8 --> D13/MW4
 --q24 orbit 85 --> D12/MW5          [native equation-side frame D24eq]
 --q6 --> A11/MW6
 --q8 --> 2A5/MW7
 --q4 --> 3A3/MW8
 --q4 --> A3+2A2/MW10
 --q4 --> 5A1/MW12
 --q4 --> 4A1/MW13
 --q4 --> 3A1/MW14
 --q4 --> 2A1/MW15
 --q4 --> A1/MW16
 --q6 --> rootless/MW17 = R17
```

The preferred D12 suffix has q-sequence

```text
6, 8, 4, 4, 4, 4, 4, 4, 4
```

through `A1/MW16`. The terminal `A1/MW16 -> rootless/MW17` q6 arrow remains the already-certified final H3 lattice neighbor.

## What changed in `bfb96bf`

The important change is not a new abstract ADE sequence; that sequence was already known. The commit adds a **native replay from the q24 equation-side D12 frame**:

- `scripts/export_h92_q24_native_d12_frame.sage` exports `D24eq`, the `D12/MW5` frame obtained directly from the selected q24/D13 equation-side geometry;
- `scripts/run_h92_q24_native_suffix_to_a1.py` searches/replays the nine-step suffix from that frame;
- the replay reaches `A1/MW16` with terminal output `PASS_Q24_NATIVE_D12_TO_A1`;
- the final q4 step has root data `(1,2,2)`, MW rank 16, and the reported selected orbit is `17593`.

This removes a route-selection ambiguity: downstream work should be compiled from the q24-native `D24eq` frame, not from q32 and not from a later historical D12 frame that merely shares the same ADE/MW label.

## Proof boundary

The route has three different levels of completion and they must not be conflated.

1. **H3 source through D13:** characteristic-zero equations and markings are exact.
2. **D13 -> D12:** q24 orbit 85 is the selected equation frontier. The equation-side D12 lattice frame `D24eq` is pinned, but the full characteristic-zero D12 Weierstrass/Jacobian + marked-section certificate is still the immediate equation task unless/until separately certified.
3. **D12 -> A1:** the new q24-native suffix is an exact lattice/root-Weyl replay from `D24eq`. It certifies which neighbors to execute, not their characteristic-zero equations.
4. **A1 -> R17:** the final q6 lattice arrow is already certified and the endpoint is pinned to the recovered rank-17 lattice.

Therefore `PASS_Q24_NATIVE_D12_TO_A1` means **the preferred suffix is now pinned to the correct q24-native D12 marking**. It does not mean the nine downstream characteristic-zero equation transforms have already been executed.

## Route priority

Use this order unless new evidence changes the cost materially:

1. finish/certify the characteristic-zero q24 `D13 -> D12` equation in the native `D24eq` marking;
2. execute the stored native `D12 -> A11` q6 neighbor and preserve the marking needed for the next step;
3. continue the native q-sequence `8,4,4,4,4,4,4,4` to `A1/MW16`, checking equation/fibre/marking data at every stage;
4. execute the certified final q6 to rootless `MW17` and compare with pinned R17;
5. only then specialize to the curve-273 endpoint.

The q32 D12 route is now **alternate/regression work**, not the preferred H3 continuation. Its spinor-quartic marking experiments remain useful compiler techniques, especially for recovering a zero section on degree-two genus-one quartics, but q32 should not displace the q24-native path while the latter remains viable.

## Canonical reproduction entry points

```text
scripts/export_h92_q24_native_d12_frame.sage
scripts/run_h92_q24_native_suffix_to_a1.py
```

The runner writes its local replay artifact to

```text
artifacts/local/elkies-k3/q24-native-suffix/q24-native-d12-to-a1.json
```

with terminal status

```text
PASS_Q24_NATIVE_D12_TO_A1
```

For route provenance and superseded approaches, see `SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`. For the larger source/comparison map, see `CONSTRUCTION_ROUTES.md`.