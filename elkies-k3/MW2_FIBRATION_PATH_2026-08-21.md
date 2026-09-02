# An exact MW2 fibration path (2026-08-21)

## Status

The target `MW <= 2` has been reached.  The path and terminal lattice data in
this note are **exact lattice computations**.  The beam that discovered the
path was bounded and order-dependent, so this is neither an exhaustive
classification nor a proof that the path or terminal frame is optimal.

A subsequent widened search reaches rank `MW=1` on the older rank-19,
discriminant-948 lattice from the canonical Kumar entrance.  See
[`MW1_OPTIMAL_FIBRATION_PATH_2026-08-21.md`](MW1_OPTIMAL_FIBRATION_PATH_2026-08-21.md).
At this 2026-08-21 checkpoint the explicit neighbour bridge from the rank-17
chain in this note to the Kumar frame was still open, so the
`q=25,4,4,4` chain was then the preferred transported reconstruction path.
That routing statement is superseded: the selected H3 equation corridor and
its rootless q12/orbit5867 endpoint are now complete.  The chain below remains
an exact low-MW lattice ancestry and historical comparison, not the current
construction or arithmetic frontier.  See [`README.md`](README.md).

The exact rational model displayed below has since acquired a third
independent section and is actually a Picard-rank-20, discriminant-43 K3.
Thus its fibration has MW rank three, not two; the rank-18 frame in
[`data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt`](data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt)
is the corrected active search target.  The rank-19 path remains an exact
lattice ancestry for the two-section sublattice, but its old rank labels do
not describe the full NS of this specialization.

The preferred path is

```text
MW17
  -- q=25, 5*5 --> MW7,  A3 + A1^7
  -- q=4,  2*2 --> MW4,  D4 + A3 + 2 A2 + 2 A1
  -- q=4,  2*2 --> MW3,  A5 + D4 + 2 A2 + A1
  -- q=4,  2*2 --> MW2,  E6 + D4 + 2 A2 + A1.
```

Thus the widened search improves the earlier terminal MW3 path by one more
small `q=4` neighbor.  The last step promotes the `A5` factor to `E6` and
leaves the other root factors unchanged.

## Exact neighbor witnesses

Coordinates in each row refer to the positive frame Gram in the preceding
row.  For `NS = U + (-F)`, the isotropic class is `f=(a,b,v)` and satisfies
`a*b=q(v)=v*F*v/2`.  Every row below passes
[`scripts/verify_fibration_neighbor.sage`](scripts/verify_fibration_neighbor.sage).

| step | `q=a*b` | `v` | child root data `(rank,count,det)` | child ADE |
|---|---:|---|---|---|
| MW17 -> MW7 | `25=5*5` | `(-1,0,-4,3,0,0,0,0,0,-1,1,0,0,0,-3,0,0)` | `(10,26,512)` | `A3 + A1^7` |
| MW7 -> MW4 | `4=2*2` | `(-1,-2,1,0,1,1,2,-3,0,-2,0,1,0,0,-1,0,0)` | `(13,52,576)` | `D4 + A3 + 2 A2 + 2 A1` |
| MW4 -> MW3 | `4=2*2` | `(-1,0,0,2,0,2,-1,1,0,0,1,0,0,0,0,0,0)` | `(14,68,432)` | `A5 + D4 + 2 A2 + A1` |
| MW3 -> MW2 | `4=2*2` | `(0,-2,-2,0,1,0,3,2,4,0,3,0,-1,-3,-1,-4,0)` | `(15,110,216)` | `E6 + D4 + 2 A2 + A1` |

The canonical frames are

- [`data/fibrations/q25_mw7_frame.txt`](data/fibrations/q25_mw7_frame.txt),
- [`data/fibrations/q25_mw4_frame.txt`](data/fibrations/q25_mw4_frame.txt),
- [`data/fibrations/mw3_a5_d4_a2a2_a1_frame.txt`](data/fibrations/mw3_a5_d4_a2a2_a1_frame.txt), and
- [`data/fibrations/mw2_e6_d4_a2a2_a1_frame.txt`](data/fibrations/mw2_e6_d4_a2a2_a1_frame.txt).

The third and fourth frame SHA-256 hashes are respectively
`0ea8f0b7df3acc76a58536bce2d9730f74ab7c2889b2e0a0a8549cd328a45dd9`
and
`d7499c13010015eaf09f62bc9c53a8e5354917f2cd512b803ec63eeabe0384a0`.

## Exact MW data and section profiles

The intermediate MW3 height lattice has reduced Gram

```text
(1/6) * [ 3  0 -1]
        [ 0  9  1]
        [-1  1 18],       det = 79/36.
```

Its expected fiber presentation is

```text
I6 + I0* + 2 I3 + I2 + 4 I1.
```

The terminal saturated MW2 height lattice has reduced Gram

```text
(1/6) * [9  2]
        [2 18],           det = 79/18.
```

The exact frame glue gives the following convenient profiles, ordered as
`(E6,D4,A2,A2,A1; P.O)`:

```text
P1 = (1, 0,  0, 1, 1; 0)
P2 = (2, d1, 1, 0, 0; 1),
P1.P2 = 2.
```

Here `1,2` are the two inverse nonzero `E6` discriminant classes and `d1` is
one of the three triality-equivalent nonzero `D4` classes.  The two `A2`
factors may be swapped, and simultaneous component inversion gives an
equivalent presentation.  The exact glue recovery and Shioda-height replay
are implemented in
[`scripts/recover_mw2_component_glue.sage`](scripts/recover_mw2_component_glue.sage).

The corresponding expected fibers are

```text
IV* + I0* + 2 I3 + I2 + 2 I1.
```

Putting `IV*` at infinity gives short-Weierstrass bounds
`deg(A)<=5`, `deg(B)<=8`.  As a reconstruction target this is materially
smaller than the `D6+A5+A3` MW3 route: it has two sections rather than three,
one polynomial and one simple-pole section rather than two polynomial and one
simple-pole section, one pair gate rather than three, and lower coefficient
degree bounds (`5,8` instead of `6,9`).  Its one pair intersection is `2`
rather than the minimal value `1`, so this comparison is about total system
size, not a claim that every individual incidence condition is simpler.  No
explicit Weierstrass model was assumed in selecting the frame; one has now
been reconstructed over `QQ`.

## Explicit rational reconstruction

The complete normalized charts at `p=13,17,23` each have exactly one target.
CRT reconstructs their common rational fiber and pole locations as

```text
lambda = 9/25,   mu = 49/25,   pole(P2) = 16/25.
```

Hensel lifting with `lambda=9/25` fixed through `23^24` reconstructs every
remaining coordinate, and literal substitution makes all 14 integral
equations vanish over `QQ`.  The resulting model is

```text
A = -(32447500/583443)t^2 -(906250/194481)t^3
    +(31250000/194481)t^4 -(19531250/194481)t^5,

B = (300827000000/2315685267)t^3 +(340001171875/1029193452)t^4
    -(498857421875/257298363)t^5 +(29541015625/10501974)t^6
    -(152587890625/85766121)t^7 +(152587890625/343064484)t^8.
```

Its sections are

```text
P1.X = 1 -(800/1323)t +(625/147)t^2,
P1.Y = 1 -(400/441)t -(394375/18522)t^2
       +(484375/9261)t^3 -(390625/18522)t^4,

P2.X = N/(t-16/25)^2,
P2.Y = M/(t-16/25)^3,

N = (77824/33075)t +(12400/1323)t^2 -(30500/1323)t^3 +(5000/441)t^4,
M = (4096/343)t^2 +(281408/9261)t^3 -(1517000/9261)t^4
    +(1296875/6174)t^5 -(1015625/9261)t^6 +(390625/18522)t^7.
```

The discriminant factorization is

```text
constant * t^6 (t-1)^3 (t-9/25)^3 (t-49/25)^2
         * (t^2 -(512/675)t -4096/16875).
```

The quadratic is squarefree, the section identities hold exactly, the
profiles are the frame profiles, and

```text
gcd(P1.X*q^2-N, P1.Y*q^3-M) = t^2 -(67/25)t +1008/625,
```

so `P1.P2=2`.  Shioda's formula gives the exact height Gram
`(1/6)[9,2;2,18]`, determinant `79/18`.  This is independently checked by
`scripts/verify_mw2_e6d4a2a2a1_qq.sage`.

### Complete `GF(23)` audit

The rational surface was first found as the unique exact anchor in the
complete normalized `GF(23)` chart.  Put the fibers `I0*`, `I3`, `I3`, and
`I2` at `t=0,1,16,13`, respectively, and `IV*` at infinity.  The surface is

```text
A = 13 t^2 + 17 t^3 + 11 t^4 + 19 t^5,
B = 20 t^3 + 7 t^4 + 3 t^6 + 19 t^7 + t^8.
```

The two sections are

```text
P1.X = 1 + 10 t + 3 t^2,
P1.Y = 1 + 15 t + 14 t^2 + 15 t^3 + t^4,

P2.X = (15 t + 6 t^2 + 19 t^3 + 8 t^4) / (t-8)^2,
P2.Y = (22 t^2 + 14 t^3 + 13 t^4 + 2 t^5 + 19 t^6 + 22 t^7)
       / (t-8)^3.
```

The discriminant has exact valuations `6,3,3,2,8` at the four finite
reducible fibers and infinity, and its remaining degree-two factor is
squarefree.  The component and pole gates recover the frame profiles, while
the polynomial gcd of the two section-incidence equations is
`t^2+t+22`, so `P1.P2=2`.

The normalized `P1` equations reduce to three equations in five variables.
All 21 possible nonzero normalized `I3` locations were exhausted: 9,366
algebraic chart hits reduce to 248 exact surfaces after the Kodaira and
squarefree-residual gates.  The fused `P2` search finds 160 square-root
orientations.  On the inverse-`IV*` orientation, the final intersection
degrees `2,3,4` occur `1,4,7` times; hence the displayed frame target is the
unique degree-two hit.

The full two-section deformation system has 13 variables, 14 nonzero
residual equations, and Jacobian rank 12 at the displayed point.  Thus it is
a smooth one-dimensional local branch over `GF(23)`.  Hensel lifting verifies
the branch through `23^24`, and a structured formal expansion verifies all 18
raw residuals through order 60.  Fixing the normalized `I3` coordinate at its
residue does not rationally reconstruct, no coordinate has a bidegree at most
`8 x 8` relation with that parameter, and none of the 78 coordinate pairs has
a bidegree at most `6 x 6` relation.  These are bounded negative tests for
that incorrect integer transverse slice; the three-prime CRT comparison is
what exposes the rational value `9/25`.

The preceding `D6+A5+A3` MW3 route was also subjected to a deterministic
bounded audit: 20,000 normalized slices gave 14,213 algebraic hits, 5,343
exact `P1` surfaces, 183 exact `P1+P2` surfaces, and no target `P3`.  This is
why the MW2 endpoint is now the preferred reconstruction frontier.

## Widened beam

The beam contained 24 distinct MW3 input frames, representing 22 root-data
types from the lower-`q` searches together with the existing
`E6+A3^2+A1^2` and `A6+A4+A1^4` frames.  For each parent it sampled proper
factor neighbors with `4 <= q <= 16`, one `a<=b` factor ordering, a baseline
cap of 200, four restarted caps of 200, and at most two retained frames per
root-data triple.  It retained 1,920 frames and found four MW2 root-data
types:

| parent type | `q` | terminal ADE | root `(count,det)` | reduced MW2 Gram |
|---|---:|---|---|---|
| `(90,192)` | `8` | `A9+A5+A1` | `(122,120)` | `(1/30)[117,24;24,248]` |
| `(68,432)` | `4` | `E6+D4+2A2+A1` | `(110,216)` | `(1/6)[9,2;2,18]` |
| `(68,432)` | `12` | `A6+A4+A3+A2` | `(80,420)` | `(1/105)[81,-39;-39,326]` |
| `(102,96)` | `14` | `A5+A4+2A3` | `(74,480)` | `(1/60)[67,-27;-27,117]` |

The `q=4` endpoint was selected because it lies on the clean
`q=25,4,4,4` ancestry and has the smallest coefficient degree bounds among
these audited endpoints.  The complete retained hit table is
[`../artifacts/generated-results/elkies-k3-mw2-beam-round1/beam-hits.tsv`](../artifacts/generated-results/elkies-k3-mw2-beam-round1/beam-hits.tsv),
with SHA-256
`3644e8a2766c2992cb63668403ec2957188b98e7f48be96ab1301347c254442e`.

## Reproduction

Run the bounded beam from the recorded input manifest with

```bash
.venv/bin/python elkies-k3/scripts/run_fibration_beam.py \
  --input artifacts/generated-results/elkies-k3-mw2-beam-round1-input.tsv \
  --out-dir artifacts/generated-results/elkies-k3-mw2-beam-round1 \
  --qmin 4 --qmax 16 \
  --enum-baseline-cap 200 --enum-restarts 4 --enum-cap 200 \
  --enum-seed 20260821 --per-root-data-cap 2 --report 80 --workers 4
```

The output directory must not already exist.  Its `commands.txt` records the
exact per-parent commands.  The search uses exact PARI `qfminim(2)` root
enumeration for each constructed frame, but its norm-vector traversal is the
bounded discovery step.

Replay the terminal invariants and component glue with

```bash
sage elkies-k3/scripts/analyze_mw3_branch.sage \
  --frame elkies-k3/data/fibrations/mw2_e6_d4_a2a2_a1_frame.txt \
  --name mw2-e6-d4-a2a2-a1

sage elkies-k3/scripts/recover_mw2_component_glue.sage
```

Build and exhaust the normalized finite-field reconstruction chart with

```bash
sage elkies-k3/scripts/build_mw2_e6d4a2a2a1_p1.sage \
  --export=artifacts/local/elkies-k3/mw2-e6d4a2a2a1-p1/p23-p1.ms

sage elkies-k3/scripts/search_mw3_a10_p1_multislice.sage \
  --seed-start=1 --seed-end=21 \
  --input=artifacts/local/elkies-k3/mw2-e6d4a2a2a1-p1/p23-p1.ms \
  --open-input=artifacts/local/elkies-k3/mw2-e6d4a2a2a1-p1/p23-p1.open.ms \
  --dir=artifacts/local/elkies-k3/mw2-e6d4a2a2a1-p1 \
  --fixed-names=lam --nonzero-keep=mu,s --prefix=p1 \
  --deterministic-fixed --max-hits=1000

sage elkies-k3/scripts/search_mw2_e6d4a2a2a1_chain_from_logs.sage

sage elkies-k3/scripts/analyze_mw2_e6d4a2a2a1_target_lift.sage \
  --hensel-digits=24 --fixed-lam=9/25

sage elkies-k3/scripts/verify_mw2_e6d4a2a2a1_qq.sage
```

The four neighbor witnesses in the table are independently replayable with
[`scripts/verify_fibration_neighbor.sage`](scripts/verify_fibration_neighbor.sage).
Full exact analyses are preserved under
`artifacts/generated-results/elkies-k3-mw2-beam-round1/`.
