# Exact residual-descent inputs for four low-conductor near misses

[MATH_STATUS.json](../../MATH_STATUS.json) is authoritative. This note is the
canonical source for `EC-CF-NEARMISS-DESCENT-INPUTS`. The complete dated
execution diary, including resource terminals and historical commands, is
[archived](../../archive/elliptic-curves/notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.before-2026-09-13.md.txt).

## Certified inputs

The retained
[`conductor_first_near_miss_descent_targets_v1.json`](../../artifacts/generated-results/elliptic-curves/conductor_first_near_miss_descent_targets_v1.json)
pins exact descent models, 2-division cubics, trivial rational 2-torsion
witnesses, source hashes, point bases, and finite-quotient certificates for
these four fixed curves:

| Target | Certified known rank | Certified image in `E(Q)/2E(Q)` |
|---|---:|---:|
| ICARM 245 | 20 | 20 |
| Fermigier--Mestre `u=28917/20` | 20 | 20 |
| split-infinity family 2, `u=483` | 19 | 19 |
| split-infinity family 3, `u=660` | 19 | 19 |

For each target, exact finite products of `E(F_p)/2E(F_p)` certify the stated
binary rank of the displayed point basis. The lower bounds and Kummer-image
dimensions are therefore exact inputs to a relative descent; they do not
assert global Mordell--Weil saturation.

## Kummer-image preflight

The older mod-3-selected bases for the two split-infinity fibres had only
one-dimensional images in the bounded mod-2 quotient product. A bounded PARI
`ellsaturation(E,P,3)` pass supplied shorter candidate point lists. Those
lists were admitted only after exact point membership and independent
finite-quotient rank-19 certificates. The bounded saturation routine is a
point-basis discovery aid, never a global saturation theorem.

This preflight prevents a relative 2-descent from spending its output on a
known subgroup index. It does not compute a Selmer group, residual cover,
rank upper bound, or new rational point.

## Replay boundary

[`build_conductor_first_near_miss_targets.py`](../cas/build_conductor_first_near_miss_targets.py)
is the retained reproducer. Its `--check` mode reconstructs the manifest from
the pinned source artifacts, redoes the exact finite-quotient certificates,
and invokes bounded PARI candidate saturation for the two Mestre targets
before byte-comparing the result. It is a mathematical replay, not a cheap
metadata check. Its recorded source digest is
`5095ac3c28a2f2f57d02bcefc3923bb25e3684748e936f39535fa6d66b563035`; the
retained output digest is
`cfd2ac0d0ec995c102df534fc216fbb849a21970ed68ce6785ce000f12cba0ff`.

[`build_conductor_first_near_miss_magma.py`](../cas/build_conductor_first_near_miss_magma.py)
can generate separately scoped Magma jobs using the pinned subgroup as
`RemoveGens`. No complete 2- or 3-Selmer computation, Cassels--Tate matrix,
residual cover, or new point has been produced by this repository. The
archived resource limits and collector experiments are historical evidence;
they do not authorize a replay or new campaign.
