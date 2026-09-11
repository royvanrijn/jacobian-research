# Accessibility gains: exact events, all subsets, and basis sensitivity

The original-basis contrast survives removal of acquisition order, but the large contrast in relative gains is sensitive to modest generic-basis changes. The old zero on 11952 was partly order-dependent. These are sparse effects in a fixed coordinate atlas, not evidence of an intrinsic avalanche or of a rank upper bound.

All labels below are one-based. M1,...,M17 are the original specialized generic basis, and E1,... are the original displayed exceptional directions. Full integer words in both the public D basis and the (M,E) basis, exact points, and N/F values are retained in [events.json](events.json). Curve302 E1,...,E14 correspond to public points 1,3,4,6,9,13,14,17,19,22,26,29,30,31; 11952 E1,...,E8 are D18,...,D25.

## Four successful events

Every contribution is before minus after, in bits; positive values help the drop. The residual column is one quarter of the rational x-height on the frozen short Weierstrass model. The identity used is k = h_x/4 + a + c, with c = log2(gcd(N,F))/4 and a = k - log2(max(|N|,|F|))/4. All N/F pairs are fixed primitive integral homogeneous quartics before target evaluation.

| Hidden target | Added direction | Q before | Q after | Total drop | Residual x-height term | Archimedean term a | Finite term c |
|---|---|---|---|---:|---:|---:|---:|
| E10 | E1 | -M3+M14 | -E1 | 17.805600 | -1.591692 | -4.630262 | +24.027553 |
| E14 | E1 | -M12+M17 | -E1 | 13.767574 | +12.385440 | -3.983428 | +5.365562 |
| E7 | E2 | M10 | E2 | 3.113518 | +4.953308 | -1.193063 | -0.646727 |
| E12 | E6 | M10-M13 | M17+E6 | 6.772843 | +0.341911 | +22.900864 | -16.469932 |

- E10 after E1: the residual x-height and archimedean terms worsen. A 24.028-bit fall in c more than offsets them. The exact parameter-height ratio is about 229,097; this is not a runtime measurement.
- E14 after E1: the main favorable term is the 12.385-bit decrease in residual x-height/4. The finite term contributes another 5.366 bits, while the archimedean term offsets 3.983 bits.
- E7 after E2: the 4.953-bit residual x-height improvement is partly offset by both remaining terms.
- E12 after E6: residual x-height/4 changes by only 0.342 bits. The archimedean term improves by 22.901 bits, while c increases and offsets 16.470 bits. The new chart is Q=M17+E6.

The signs matter: a larger c is not automatically favorable. These accounts are exact coordinate-height identities, not fitted explanations or canonical-height causal claims. The verifier checks the multiplicative identity H_before^4/H_after^4 = residual_ratio * archimedean_ratio * finite_ratio and checks x(2P-Q)=N/F using a separate rational group law.

## Complete subset landscape in the original basis

| Metric | Curve302 | 11952 |
|---|---:|---:|
| Known exceptional directions | 14 | 8 |
| Held-out target/subset cells | 114,688 | 1,024 |
| Directed subset edges | 745,472 | 3,584 |
| Targets helped somewhere | 4/14 | 1/8 |
| Largest single-addition drop (bits) | 17.805600 | 1.654608 |
| Positive-G cells | 18,432 (16.071%) | 64 (6.250%) |
| Median G over cells (bits) | 0.000000 | 0.000000 |
| Mean positive part of G (bits/cell) | 1.549438 | 0.103413 |
| Positive edges | 18,432/745,472 | 64/3,584 |
| Positive-edge AUC / all edges (bits) | 0.227800 | 0.029547 |

Each target has every subset of the other exceptional directions represented. At subset size k the control has exactly 2(17+k)^2 generic-only charts. Full discrete survival curves, their exact rational height ratios, and all winner indices are in the compressed per-policy files. Cell and edge means have explicitly different denominators; neither is an independent-sample statistic.

The original-basis functions reduce to the following minimal offers. For any target, subtract the largest listed improvement whose required subset is supplied from its initial kappa; use zero when no offer is available. Every unlisted target stays at its initial kappa for every eligible subset. The independent verifier checks this compact description against every stored cell.

| Fibre | Target | Required exceptional direction(s) | Winning anchor | Improvement over M17 (bits) |
|---|---|---|---|---:|
| 302 | E7 | E2 | E2 | 3.113518 |
| 302 | E10 | E1 | -E1 | 17.805600 |
| 302 | E10 | E2 | E2 | 3.849481 |
| 302 | E12 | E6 | M17+E6 | 6.772843 |
| 302 | E14 | E1 | -E1 | 13.767574 |
| 11952 | E2 | E4 | E4 | 1.654608 |

11952 E2 gains 1.654608 bits from Q=E4. The old ordering supplied E2 before E4 and could not expose this held-out improvement. On 302, E2 also offers E10 a 3.849481-bit improvement, which the old path masked by supplying the better anchor -E1 first. All six minimal offers across both fibres have singleton exceptional support. Thus no pair of exceptional directions supplies an additional gain beyond these offers in the original full signed-unit/pair atlas.

## Frozen basis sensitivity

The five perturbations use disjoint consecutive pairs (a,b) -> (a+s*b,b), with s=+1 or -1. They are unimodular, preserve M17, and never mix generic and exceptional blocks. Generic transformations are shared across targets. Exceptional transformations pair the other exceptional indices separately for each held-out target, keeping that original target fixed and excluding its coordinates from every supplied generator. These are leave-one-target-out basis sensitivity probes, not a single simultaneous rebased acquisition path. All policies were pinned before their measurements.

| Basis policy | 302 targets helped | 11952 targets helped | 302 largest drop / max G | 11952 largest drop / max G | 302 mean positive G | 11952 mean positive G |
|---|---:|---:|---:|---:|---:|---:|
| original | 4/14 | 1/8 | 17.806 / 17.806 | 1.655 / 1.655 | 1.549 | 0.103 |
| generic_plus | 7/14 | 4/8 | 18.293 / 18.293 | 16.380 / 16.380 | 2.396 | 3.316 |
| generic_minus | 6/14 | 4/8 | 18.293 / 18.293 | 12.094 / 12.094 | 2.154 | 2.551 |
| exceptional_plus | 4/14 | 1/8 | 17.806 / 17.806 | 1.655 / 1.655 | 0.986 | 0.052 |
| exceptional_minus | 4/14 | 1/8 | 17.806 / 17.806 | 1.655 / 1.655 | 0.986 | 0.052 |
| both_plus | 7/14 | 4/8 | 18.293 / 18.293 | 16.380 / 16.380 | 1.737 | 2.338 |

Relative drops are sensitive to a worse initial bank. For 11952 E2 under generic_plus, initial kappa increases from 72.929942 to 87.655814 bits, while the best held-out chart remains Q=E4 at 71.275334 bits. Its new 16.380480-bit drop is the original 1.654608-bit effect plus a 14.725872-bit degradation of the starting bank. The endpoint has not improved.

More generally, the four helped 11952 targets under generic_plus end at equal or worse absolute costs than under the original full held-out bank. The rebase changes which good generic charts are available; exceptional anchors can restore accessibility relative to that altered baseline. The original-basis advantage therefore does not establish intrinsic exceptional subgroup geometry. Exceptional-block shears alone preserve the strongest original gains but alter how many subsets make the required anchor available.

The supported conclusion is a small collection of useful anchor relations with presentation-dependent accessibility gains. Acquisition order explains the old 11952 zero; generic-basis choice explains much of the dramatic cross-fibre contrast in relative drops. The fixed raw-coordinate metric and one prescribed generic control policy remain limits. These results give no information about missing rank on 11952.

## Reproduction and verification

Run from the repository root:

```bash
sage -python research/elliptic-curves/cas/rank_accessibility_subsets.py prepare
sage -python research/elliptic-curves/cas/rank_accessibility_subsets.py original
sage -python research/elliptic-curves/cas/rank_accessibility_subsets.py all
python3 research/elliptic-curves/cas/verify_rank_accessibility_subsets.py --check
python3 research/elliptic-curves/cas/report_rank_accessibility_subsets.py --check
```

The producer resumes immutable hash-bound checkpoints. The separate verifier uses Python Fraction arithmetic for group operations and checks every distinct rational slope, every subset minimum via direct sparse support predicates, generic chart-prefix minima, and complete edge/G survival counts. It checks the four original exact event identities separately. Logs are descriptive float displays; all cost comparisons use exact integers.
Across six policies on two fibres, verification replays 149,916 distinct target–anchor slopes, 694,272 held-out cells, and 4,494,336 directed edges. Repeated states are not treated as independent observations.
