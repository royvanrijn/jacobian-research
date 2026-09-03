# Alternate-Q80 rational V4 bases and product-twist laboratory

<!-- status-consumer: EC-K3-R17-NORM12-11952-V4-RATIONAL-PAIR-SHORTLIST-64 e14368b602eebedb -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-V4-BASE-RANK-SCREEN-64 b12acf0b90056c18 -->

## Exact result

The cheapest-1,024 native alternate-Q80 bisection prefix contains many
arithmetically usable genus-one pair bases despite having no shared branch
fibres.

For a norm-ten bisection

```text
C_w = (2,2,w),
```

the direct-frame intersection is

```text
C_w.C_v = 8 - <w,v>.
```

Exact evaluation of all `523,776` unordered pairs finds `10,362` with
intersection number one.  Two rational curves with intersection number one
have a unique geometric intersection point.  It is fixed by Galois and gives
a rational point on the corresponding `V4` fibre product.

The production certificate selects the first 64 such pairs in the declared
priority/complexity order.  For every pair it also recovers the point from the
two explicit lifted-section formulas.  Substitution verifies

```text
s^2 = q_i(u),
t^2 = q_j(u),
(st)^2 = q_i(u)q_j(u),
```

and verifies that both lifts give the same point of the alternate elliptic
surface.  No bounded rational-point search is used in this proof.

All 64 bases are connected genus-one curves: their two irreducible quadratic
branch divisors are disjoint.  The two existing character sections therefore
give generic rank at least 19 over each base.  The shortlist certificate alone
does not determine the Jacobian ranks.

The exact artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json),
SHA-256
`4739e0b24b00e276269b228d7c85a01743caa1dcbc62ef58ea59ca3152d116a0`.

## Exact base-rank screen

Every base Jacobian was isolated in its own Sage process and given an exact
PARI `ellrank` call with effort zero and a hard ten-second limit.  The bounded
screen gives:

| interval | bases |
|---|---:|
| `[1,1]` | 17 |
| `[0,2]` | 27 |
| `[0,4]` | 8 |
| `[1,3]` | 10 |
| timeout / `UNKNOWN` | 2 |

Thus seventeen shortlisted bases have certified Mordell--Weil rank exactly
one.  No ambiguous interval is treated as an exact rank.  The rank artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json),
SHA-256
`ed41526fca5724f529aaf73d3f695fe2a781eaf4ea6b6fb792a2af5676759f69`.

## Product-twist compiler and pilot

For each shortlisted pair the missing third character is the quartic twist

```text
E_alt^(q_i q_j):
Y^2 = X^3 + A(u)(q_iq_j)^2 X + B(u)(q_iq_j)^3.
```

This has `chi=4`.  A section disjoint from the zero section lies in the
polynomial box

```text
deg X <= 8,   deg Y <= 12.
```

The modular exporter now accepts direct alternate product keys.  It requires
squarefree good reduction away from the `24I1` discriminant, chooses a smooth
chart fibre without rational 2-torsion, fixes every affine leading point, and
exports the recursively reduced systems.  The batch driver records all
limits and delegates each system to the existing checkpointed `msolve`
runner.

For the first certified rank-one pair

```text
alternate-orbit-1463f : alternate-orbit-19bad
```

the exact `p=131` export has 140 leading-point blocks and 70 distinct systems,
each with twelve equations in eight variables.  A single declared 15-second
pilot group timed out after covering its two sign-equivalent blocks; its
sampled process-tree high-water mark in the promoted replay was 999,956 KiB.
This runtime measurement can vary across replays and is a complexity
measurement only.  It is not a negative section result.

Accordingly, launching the full shortlist-by-prime collection of unsliced
systems is not justified by this pilot.  The next computational gate is a dimension-reducing ansatz or a
cheaper elimination/specialization sieve applied to these exact exports.

The checkpoint manifests are

- [`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-quartic-po0-campaign-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-quartic-po0-campaign-v1.json);
- [`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-po0-direct-product-alternate-orbit-1463f--alternate-orbit-19bad-p131-msolve.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-po0-direct-product-alternate-orbit-1463f--alternate-orbit-19bad-p131-msolve.json).

## Replay

```bash
sage -python \
  elkies-k3/scripts/select_r17_norm12_11952_v4_pair_shortlist.sage

sage -python \
  elkies-k3/scripts/select_r17_norm12_11952_v4_pair_shortlist.sage --check

.venv/bin/python \
  elkies-k3/scripts/screen_r17_norm12_11952_v4_base_jacobian_ranks.py \
  --limit 64 --timeout 10 --jobs 2

.venv/bin/python \
  elkies-k3/scripts/screen_r17_norm12_11952_v4_base_jacobian_ranks.py \
  --check

.venv/bin/python \
  elkies-k3/scripts/run_r17_norm12_11952_v4_quartic_po0_campaign.py \
  --pair-limit 1 --primes-per-pair 1 --max-groups 1 --timeout 15

.venv/bin/python \
  elkies-k3/scripts/run_r17_norm12_11952_v4_quartic_po0_campaign.py --check
```

## Claim boundary

The 64 rational pair-base points, their genus-one geometry, and the seventeen
exact rank-one base Jacobians are established.  The modular campaign remains
an incomplete bounded search.  No section on any quartic product twist, no
third nontrivial `V4` character carrying a section, and no generic
rank-at-least-20 `V4` base change is claimed.
