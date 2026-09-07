# R17 rational-quadratic MW20 handoff (2026-09-04)

## Stop state

This programme is **paused with the construction still `UNKNOWN`**.  No
rational-base quadratic character of the native R17 fibration has a certified
twist rank at least three.  Consequently there is no certified MW20 surface
over `QQ(t)`, and the requested tail-survival comparison has not been run.
No entry is promoted into `MATH_STATUS.json`; `STATUS.md` is therefore left to
the repository renderer and was not edited for this paused exploratory result.

The detailed result ledger is
[`R17_RATIONAL_QUADRATIC_MW20_SEARCH_2026-09-04.md`](R17_RATIONAL_QUADRATIC_MW20_SEARCH_2026-09-04.md).
This handoff records the coordinate choice, completed gates, exact proof
boundaries, and the first unrun continuation step.

## Correct target and coordinate

For a quadratic character `q`, the quadratic-base rank split is

```text
rank E(QQ(s)) = rank E(QQ(u)) + rank E^(q)(QQ(u)).
```

Thus twist rank at least three is the required lower-bound certificate.  A
Nagao score, a repeated modular cover, or three isolated finite-field sections
does not certify it; the endpoint needs three characteristic-zero sections and
an exact rank-three height pairing.

The five requested controls are fibres of the exact
`norm12-orbit-074d9` lineage.  They are not rational fibres of the compact
published Weierstrass coordinate.  Search and control transport must therefore
be performed on the native `074d9` fibration or one of its eight certified
rational-`PGL2` coordinates.  This distinction is essential.

## Completed search ledger

### Bounded heuristic screens

Six disjoint blocks of eight primes from 211 through 491 were used only as a
ranking sieve.

| family | characters scored, with chart multiplicity | best weakest-block score | conclusion |
|---|---:|---:|---|
| compact published coordinate, `q=t^2+b*t+c`, `H=100` | 40,380 | 0.251 | no rank claim |
| all eight `074d9` coordinates, same box | 323,040 | 0.256 | no rank claim |
| five control anchors, all eight coordinates, `H=50` | 407,440 | 0.436 | no rank claim |

The top anchored diagnostic is in chart `norm12-orbit-08aaa`, anchored at
curve 356:

```text
q(z) = 1 + 46*(z-z_356) + 40*(z-z_356)^2.
```

It has no known new section and is not an MW20 candidate.  The three stored
ledgers are, respectively:

- `elkies-k3-r17-rational-quadratic-twist-nagao-h100-v1.json`, SHA-256
  `8db623ca4f66102a4c470f89c6240520809fb8d2f1e9c737e15103fdc32a22a2`;
- `elkies-k3-r17-rational-quadratic-twist-pgl8-nagao-h100-v1.json`, SHA-256
  `9a372e5840c43415c02524a15ab334442f4f8c9feb4f66849f54046148485dfd`;
- `elkies-k3-r17-rational-quadratic-twist-control-anchored-pgl8-nagao-h50-v1.json`,
  SHA-256
  `3d2fee8c139926f498cc265cffd670ee579fd786b149da8af6fef96cfb8e4b69`.

All lie under `artifacts/generated-results/`.  Their `proof_boundary` fields
explicitly deny a Mordell--Weil rank conclusion.

### Exact low-genus collision gate on alternate chart `11952`

The complete stored smooth layer and the processed affine arithmetic-genus
two and three layers have no scalar-sensitive quadratic cover occurring in
three directions under the stated simultaneous-good-reduction hypotheses.
For same-trace-mask triples, the surviving mask intersection over
`17,23,29,31,37` is `65,4,0,0,0`.  The remaining four distinct-mask
genus-three pairs are eliminated at `43,47,53`.

This is exact only for the processed `11952` charts.  It is neither a result on
`074d9` nor a global twist-rank upper bound.  The compact certificate is
`elkies-k3-r17-norm12-11952-low-genus-rank3-cover-collision-v1.json`, SHA-256
`16ee728a17b0257472a8a0b5ca490a0fbc226c173dce26c9d9651214d7f36f76`.

### Exact native `074d9` polynomial-x subchart

At the good prime 19, the trace-zero polynomial-x census checks all
`19^5=2,476,099` coefficient vectors.  Three vectors survive, but their three
scalar-sensitive quadratic cover keys are distinct.  Hence this subchart has
no pair or triple collision modulo 19.

The omitted general trace-zero form is `x=X/q`; this calculation covers only
polynomial `x`.  It is not a rank upper bound.  The certificate is
`elkies-k3-r17-074d9-tracezero-genus5-polynomial-x-p19-v1.json`, SHA-256
`c22068ae09b96b58b5298f943e60f6c00d5a79d04b48d64eb079a3696c10bda2`.

## Control transport and the unperformed tail test

For a future character `q(z)=a*z^2+b*z+c`, a control at `z_i` has a rational
lift precisely when `q(z_i)` is a rational square.  The exact preflight is
implemented by `scripts/check_r17_mw20_control_transport.py`.

For the top anchored diagnostic above, only curve 356 has rational lifts:

| curve | old tail beyond 17 | rational lift on this diagnostic | reference budget if a future exact-rank-20 cover splits |
|---:|---:|---|---:|
| 351 | 8 | no | 5 |
| 356 | 12 | yes, two lifts | 9 |
| 376 | 5 | no | 2 |
| 377 | 6 | no | 3 |
| 385 | 12 | no | 9 |

The last column is only the numerical budget obtained by subtracting 20 from
the certified fibre-rank lower bound.  It becomes a tail statement only after
an exact generic rank-20 subgroup specializes injectively.  If the generic
rank is larger than 20, the displayed budget must be recomputed.

The preflight certificate is
`elkies-k3-r17-mw20-control-transport-anchored-heuristic-top-v1.json`, SHA-256
`3f3541566044d6193b02dbd9b0e96ad44b6face93e388f3adaa96c7d9aa7f05a`.
It proves only rational splitting.  The actual survival test remains:

1. certify three independent twist sections over `QQ(u)`;
2. specialize them at every rational control lift;
3. compare their span exactly with the certified exceptional quotient basis
   `P18,...` at that fibre.

No step in this three-part comparison has been possible because step 1 is
open.

The alternate-Q80 controls—12/395 in class `11952` and 363/364/378 in class
`08f72`—are different rational-`PGL2` `j`-map classes.  Transport through a
`074d9` quadratic base change is undefined until an explicit common-K3
birational transport is supplied.  They were therefore not mixed into the
five-control test.

## Exact replay commands

The compact exact gates and diagnostic control preflight replay with:

```bash
python3 elkies-k3/scripts/analyze_r17_norm12_low_genus_rank3_collisions.py --check

sage -python elkies-k3/scripts/search_r17_tracezero_genus5_normalizations_modp.sage \
  --model artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json \
  --prime 19 \
  --output artifacts/generated-results/elkies-k3-r17-074d9-tracezero-genus5-polynomial-x-p19-v1.json \
  --check

python3 elkies-k3/scripts/check_r17_mw20_control_transport.py \
  --candidate-artifact artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-control-anchored-pgl8-nagao-h50-v1.json \
  --finalist-rank 1 --generic-rank 20 \
  --output artifacts/generated-results/elkies-k3-r17-mw20-control-transport-anchored-heuristic-top-v1.json \
  --check
```

Each heuristic artifact stores its full `reproducing_command`; append
`--check` for a bitwise payload replay modulo the intentionally ignored runtime
field.  These are longer searches and were not rerun merely to prepare this
handoff.

## Resume point

The first unrun theorem-directed gate is the full native monic-quadratic
trace-zero `P.O=0` chart:

```text
q=u^2+q1*u+q0,
Y^2=X^3+q^2*A*X+q^3*B,
deg(X)<=6, deg(Y)<=9.
```

`scripts/export_r17_rational_quadratic_tracezero_po0_msolve.sage` exports the
complete finite coefficient schemes modulo a chosen good prime, including the
degree-drop block `y9=0`, and saturates away `disc(q)=0`.  Its final formulation
has passed syntax/help loading but has not been exported or solved into a
pinned artifact.  Old local recursive-elimination exports predate this
formulation and must not be treated as evidence.

A clean restart begins with:

```bash
sage -python elkies-k3/scripts/export_r17_rational_quadratic_tracezero_po0_msolve.sage \
  --prime 19 \
  --output-dir artifacts/local/elkies-k3/r17-074d9-monic-q-po0/p19/systems-full-y \
  --summary artifacts/local/elkies-k3/r17-074d9-monic-q-po0/p19/export-full-y.json
```

Then solve the exported blocks, group finite-field solutions by the full
scalar-sensitive `(q0,q1)` character, intersect surviving collisions across
good primes, lift any collision to characteristic zero, and run the height
pairing.  Even a complete miss in this `P.O=0` chart would leave higher
intersection-height sections and anchored nonmonic cover charts open.
