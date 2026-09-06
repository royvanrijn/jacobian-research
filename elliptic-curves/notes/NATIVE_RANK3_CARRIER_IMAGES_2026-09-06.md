# Twelve fresh parameter images of the marked rank-three carrier

**Twelve distinct fibres now have certified rank at least19. All exceed
the fixed400-bit normalized integral model gate; further point exposure
remains open.**

A separate finite proof certifies all17 inherited points plus both native
points on every fibre. All twelve19-point proofs and their replay pass.
The earlier height gate governed further point searching; it does not limit
certification of existing witnesses. The
[point certificates](../../artifacts/generated-results/elliptic-curves/native_rank3_carrier_subgroups_v1.json)
and [executable twelve-curve Sage export](../../artifacts/generated-results/elliptic-curves/new_native_carrier_rank19_curves.sage)
are available, and the export executes successfully. No exact-rank or
record claim follows.

The new [marked genus-one carrier](../rank-jump/MINIMAL_CARRIER_AND_RATIONAL_SPLITTING_OF_A_TWO_DIRECTION_BLOCK.md)
for native directions `0b2d0` and `19e45` has a degree-four map to the
published08234 parameter line. Its auxiliary Jacobian has exact rank three;
the existing descent returns two explicit points. This audit uses only those
two recorded points, without claiming they generate the full auxiliary group.

The fixed words are the four signed single points, four signed pair sums,
and four signed doubles. There is no auxiliary point search or word
extension. Exact pointed-quartic inversion gives twelve finite rational
parameters. Both defining quadratic forms are squares at every image, and
the native point maps satisfy the original elliptic equation exactly.
Published-to-compact08234 transport is also exact. All twelve compact
parameter heights exceed4096; their bit lengths range from33 to314.
Their twelve j-invariants are different.

The smallest displayed integral equation comes from word[1,0], at

```
compact08234 parameter = -7119612289/2394065174
```

Its largest coefficient has543 bits. The invariant inequality proves that
any normalized integral model in its rational isomorphism class has a
coefficient of at least519 bits. Every other sampled image also exceeds400
bits after arbitrary rational scaling and normalization. This comparison was
frozen against the current outer-MW16 models, whose displayed sizes reach396
bits; it is a finite search-budget gate and has no rank-upper-bound meaning.
No image entered the original height-gated finite-independence branch.
The subsequent existing-witness proof above certifies rank at least19
without any further point search.

For integral source invariants c4,c6 put G=gcd(abs(c4),abs(c6)). For any
rational scale u giving integral target invariants,

```
abs(c6_target)^2 >= c6^2 / G^3.
```

A normalized integral equation with a1,a3 in{0,1}, a2 in{-1,0,1}, and
abs(a4),abs(a6)<=M has abs(c6_target)<=1224M+521. The strict integer
inequality therefore excludes all normalized400-bit integral models of
each image, including normalized global minimal and integral short models.
It gives no height bound for other points on the rank-three carrier.

The [image certificate](../../artifacts/generated-results/elliptic-curves/native_rank3_carrier_images_v1.json)
retains every parameter, root, equation and native point. Its Sage build and
replay complete under180 seconds per stage and2GiB. A separate
[ordinary-Python replay](../../artifacts/generated-results/elliptic-curves/native_rank3_carrier_images_replay_v1.json)
independently reconstructs the homogeneous models, square roots, point
identities, rational scales and all twelve strict height inequalities.
Both independent replay commands pass. A separate
[post-construction comparison](../../artifacts/generated-results/elliptic-curves/native_rank3_carrier_images_comparison_v1.json)
finds no rational-isomorphism match among593 pinned catalogue equations and
917 prior equations, with its replay also passing. This is relative novelty;
the auxiliary descent is not rerun and universal novelty is not claimed.

Sources: `../cas/audit_native_rank3_carrier_images.sage` and
`../cas/verify_native_rank3_carrier_images.py`. The frozen twelve-word
protocol and execution ledger are under
`artifacts/local/elliptic-curves/native-rank3-carrier-images-v1`.
The current point-search priority remains the
[sixty outer-band MW16 fibres](MW16_OUTER_PARAMETER_BANDS_2026-09-06.md).

The [explicit third auxiliary direction](NATIVE_CARRIER_THIRD_DIRECTION_2026-09-06.md)
now closes the missing-point gap in the rank-three carrier subgroup and
supports a separate fixed125-word parameter-image audit.
