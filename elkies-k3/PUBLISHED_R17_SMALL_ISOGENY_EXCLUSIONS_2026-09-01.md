# Published-R17 small-prime isogeny exclusions

<!-- status-consumer: EC-K3-ELKIES-2026-R17-SMALL-ISOGENY-EXCLUSIONS fc2c4caaa79fb36c -->

## Outcome

The 2024 Elkies--Klagsbrun rank-29 curve and ICARM curves 398--400, 273,
and 302 are not cyclic `ell`-isogenous images of a rational fibre of the
published R17 fibration for any

```text
ell in {3,5,7,11}.
```

For each target and degree, the exact checker specializes the classical
modular polynomial at the target `j`-invariant and forms

```text
Phi_ell(j_R17(t),j(E)).
```

Writing the reduced published map as `j_R17=N/D`, it clears
`D^(ell+1)` and factors the resulting projective polynomial modulo clean
primes.  The witness primes are:

| target | `ell=3` | `ell=5` | `ell=7` | `ell=11` |
| --- | ---: | ---: | ---: | ---: |
| Elkies--Klagsbrun 2024 / ICARM 12 | 131 | 137 | 137 | 131 |
| ICARM 398 | 131 | 137 | 137 | 211 |
| ICARM 399 | 137 | 157 | 157 | 151 |
| ICARM 400 | 151 | 131 | 151 | 137 |
| ICARM 273 | 151 | 151 | 131 | 157 |
| ICARM 302 | 137 | 227 | 157 | 137 |

Every displayed reduction retains the full projective degree
`24*(ell+1)`, has no linear factor over `F_p`, and has no root at infinity.
Thus it has no point in `P1(F_p)`.  A rational characteristic-zero parameter
would reduce to such a point at every clean prime, so one witness suffices for
each exclusion.

The machine certificate records the complete factor-degree multiset,
coefficient hash, projective-infinity gate, target invariants, input hashes,
and software version for all 24 tests:

[`elkies_2026_published_r17_isogeny_exclusions_v1.json`](../artifacts/generated-results/elliptic-curves/elkies_2026_published_r17_isogeny_exclusions_v1.json).

## Boundary

This is a statement only about rational parameters on the published R17
fibration.  It excludes cyclic isogeny degrees `3,5,7,11`; it does not exclude
composite degrees, algebraic fibre parameters, another elliptic fibration on
the K3, or another K3 family.

In particular, it supplies no evidence for or against the equation-open
alternate Q80-derived rootless frame.  The available `GF(73)` terminal Q80
equation is the CM24 specialization with roots `4A2+A3+A5`, not the generic
rootless endpoint, and cannot be used for the proposed alternate-fibration
recognition gate.

## Replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_elkies_2026_r17_isogeny_exclusions.sage --check
```
