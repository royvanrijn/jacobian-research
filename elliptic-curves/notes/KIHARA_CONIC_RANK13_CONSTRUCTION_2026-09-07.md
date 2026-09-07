# A primitive rank-thirteen conic family over the new Kihara parent

The first new Kihara K3 parent now has an explicit rational conic base
change carrying a **primitive rank-thirteen subgroup over Q(s)**. Its
surface has **chi(O)=4**, with forty I1 fibres and two I4 fibres. This
extends the parent's full rational rank-twelve group by one independent
direction. The full generic rank after base change remains **UNKNOWN**.

A single deterministic control on this family has **thirteen certified
independent rational points** on a 998-bit model. It has no match in the
pinned 620-equation catalogue or 201-equation inventory. No point search
ran, no near-record curve was obtained, and the control is not added to the
high-rank inventory. The tested K3-parent count remains unchanged: this is
a new base change of an existing parent, not an additional K3 parent.

Authority: `EC-KIHARA-CONIC-PRIMITIVE-RANK13-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json). The
[construction](../../artifacts/generated-results/elliptic-curves/kihara_split_infinity_v1.json),
[control](../../artifacts/generated-results/elliptic-curves/kihara_conic_control_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/kihara_conic_replay_v1.json)
retain equations, points, transformations and exact witnesses.

## Why this construction was selected

The [first-parent proof](KIHARA_FIRST_PARENT_PICARD_AND_RANK_2026-09-07.md)
shows that changing fibrations on this K3 cannot supply arithmetic generic
MW16 or MW17: its rational NS rank is seventeen. That bound does not apply
to a higher-chi base change. The existing Kihara rank-fourteen path varies
both parent and fibre and has chi36. A low-degree base change of a fixed
new parent is a different, explicit way to add a direction.

This does not establish that chi4 has better coefficient sizes or discovery
yield than the rank-fourteen path. It begins with one fewer certified
direction than that path. The purpose of this calculation is to construct
and validate the extra direction before proposing population exposure.

Splitting the infinity points of a quartic by a conic is established prior
art, already used by the repository's
[Mestre D-square family](MESTRE_DSQUARE_FOUR_SCREEN.md). The result here is
its exact application, group calculation and primitive-basis proof on the
new fixed Kihara parent. No new general construction method is claimed.

## Explicit conic and section

Use the normalized six-root quartic of the parent attached to the Kihara
path parameter 3/2. Its coefficient of x^4 is

```
a(T) = lambda^2 (T^2+c),
c = 30982012044572382281613 / 26507921388877509235232.
```

The nonzero rational lambda is given exactly in the certificate. Put

```
T = (s^2-c)/(2s),
w = lambda (s^2+c)/(2s).
```

Then w²=a(T), so the two quartic infinity points are rational over Q(s).
This is a degree-two extension of Q(T), with deck involution
`gamma(s)=-c/s`.

Translate the old quartic origin to x=0, writing its equation as

```
y^2 = a x^4+b x^3+c0 x^2+d x+q0^2.
```

In the short Jacobian model used by the parent certificate, the positive
infinity point maps to

```
x_Q = 18 q0 w + 3 c0,
y_Q = 27 (q0 b + d w).
```

The exact equation is verified. The polynomial model clears the base-map
denominator with scale `(2s)^2`, multiplying x and y by that scale squared
and cubed. Its A,B,Delta degrees are 16,24,44. Delta has order four at
s=0 and forty remaining simple zeros; c4 is coprime to Delta. Infinity
has I4. Thus chi=4 and the complete configuration is 40I1+2I4.

All twelve old sections and all four infinity images under

```
s, -s, c/s, -c/s
```

were fixed before the height calculation. The four maps preserve T² and
form a V4 group, with the required Weierstrass scaling retained. The
sixteen displayed sections span rank thirteen, not fourteen or sixteen.
Every dependency is checked by exact function-field group law.

## Independence and saturation from the deck action

Let K be the first element of the parent's full rational section basis,
the hyperelliptic involution of the old quartic origin. On the cover,

```
gamma(Q)=K-Q,
R=2Q-K,
gamma(R)=-R.
```

All old sections are gamma-invariant. The exact coordinate identities show
R is nonzero. On 40I1+2I4, Shioda's height formula gives every nonzero
geometric section height at least

```
2 chi - maximum component correction = 8-2 = 6.
```

In particular geometric torsion is trivial. Applying 1-gamma to an integer
relation between Q and the old twelve sections shows that Q is independent
of the old group. Hence the constructed subgroup has rank thirteen.

The independently computed height of R is **10**. It is obtained by
multiplying R by four, so both I4 component corrections vanish, and using

```
height(R) = (8 + 2(4R.O))/16,
2(4R.O) = max(deg denominator(x(4R)), deg numerator(x(4R))-8).
```

If R were divisible by n>=2 in the geometric group, the divisor's height
would be at most 10/4, contradicting the lower bound six. Thus R is
primitive in its rational span.

Now suppose nP belongs to the subgroup M generated by the old full group
and Q, say nP=mQ+B with B old. Applying 1-gamma gives
`n(P-gamma(P))=mR`. Primitivity of R forces n to divide m; this can also be
seen by taking a Bezout combination if m,n share a proper common divisor.
Then `P-(m/n)Q` is fixed by gamma and lies in the old group, which is already
proved full over Q(T). Therefore P belongs to M.

**M is saturated in its rational span inside E(Q(s)).** This does not claim
that M is the whole generic MW group; additional independent directions
remain possible. The original constructor left saturation unknown. This
subsequent independent deck-action argument resolves that specific gap
without changing its frozen artifact.

The height Gram matrix also replays independently. Heights on the old
basis double under the degree-two base change. R is orthogonal to every
old section by the deck action, so

```
<Q,P> = <K,P>/2  for old P,
height(Q)=(height(K)+height(R))/4=13/2.
```

This reconstructs the complete thirteen-dimensional Gram matrix from the
prior parent theorem and one new height. Its determinant is **1935360**.
It matches the constructor's complete sixteen-point height calculation,
without repeating that calculation's 240 sum/difference evaluations.

## One balanced rational control

Write c=d*r² with squarefree positive integer d. Here

```
d = 599323759674,
r = 321543 / 230251694408.
```

The protocol fixes `v=floor(sqrt(d))+1=774161`, then s=r*v. The resulting
old parent parameter is

```
s = 248926050423 / 230251694408,
T = 480464663121 / 356503763989183376.
```

This choice balances the conic factors near sqrt(d); it uses no arithmetic
score, known high-rank fibre, point-search outcome or validation prime.
Large v does not mean large real T. The displayed rational T has large
arithmetic height, and the exact fibre model has 998 coefficient bits.
No claim of optimal normalization or matched height is made.

The intake undoes the polynomial model's artificial scale, clears only
necessary denominator valuations, and removes common powers at those
same known primes and at 2,3. It does not factor the discriminant or claim
a global minimal model. All thirteen points are transported exactly.

The fixed prime bound 997 supplies a full mod-two finite-column rank of
thirteen and a rational-two-torsion exclusion witness. The independent
replay enumerates complete finite elliptic groups and their quotients by
doubles. It checks every point and the model transport. This proves a
rational rank lower bound thirteen, not exact rank.

The catalogue comparison is post-construction and limited to the pinned
620 catalogue equations and 201 inventory equations. No matches were
found; literature-wide novelty is not established. This one control is
not a population experiment or a point-search rank gain.

## Reproducibility and decision boundary

Copy the [bundle](../../artifacts/generated-results/elliptic-curves/kihara_conic_replay_bundle_v1.json)
and [independent verifier](../cas/verify_kihara_conic.sage) to an empty
directory and run:

```sh
sage -python verify_kihara_conic.sage --input kihara_conic_replay_bundle_v1.json
```

The parent rank and full-basis theorem is an explicit dependency, with its
prior independent replay bound into this bundle. The conic verifier
imports no repository code. It verifies the maps, all infinity images,
deck identities, new height, primitive subgroup argument's arithmetic
witnesses, and the independent thirteen-point control certificate.

The constructor completed in113.997334907 seconds under its120-second cap.
The control took1.307479128 seconds. The independent proof took0.900870137
seconds. All stages used one worker and at most2GiB; the total recorded
supervised time is **116.205684172 seconds**.

```sh
python3 elliptic-curves/cas/record_kihara_conic_replay.py --check
```

This supplies a usable new family and a certified seed. It does not yet
justify a larger scan or establish an advantage over the existing
rank-fourteen path. Further exposure needs its own mathematical question,
fixed selection and completed-exposure accounting. No point search or
following campaign was launched here.

The elliptic-surface and height ingredients are standard; see
[Schütt–Shioda, *Elliptic Surfaces*](https://arxiv.org/abs/0907.0298).
