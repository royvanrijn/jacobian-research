# Historical H3 q=8 frontier and repair

> **Historical snapshot (2026-08-22).** This note records the q8 repair and is
> still authoritative for that edge, but it is not the current programme
> frontier. The equation route now passes q4/orbit164, q8/orbit376, and
> q12/orbit5867 to the certified rootless `24I1/MW17` endpoint.  The current
> programme frontier is residual 2-descent in the compact published `t`
> chart; start from [`README.md`](README.md) and
> [`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md).

Status at snapshot: 2026-08-22, after the binary-quartic 2-cover and
q-normalizer repairs.

The authoritative audit/repair note is [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md).

## Exact route now certified

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4.
```

Both displayed hops are now exact at equation level over characteristic zero.

The q8 certificate is

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --representative component-nef \
  --output artifacts/local/elkies-k3/q8-target-component-nef-audit.json
cmp artifacts/local/elkies-k3/q8-target-component-nef-audit.json \
  elkies-k3/data/fibrations/h3_q8_component_nef_physical_root_target.json
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_corrected2cover_qq.sage
```

The first command rebuilds the pinned 17 KB component-nef target and `cmp`
checks it byte-for-byte. The q8 checker consumes the tracked fibration-data
copy rather than an ignored `artifacts/local/` target prerequisite.  However,
the final q8 command still requires the ignored q6 child-Jacobian, transported-
zero, and E7-infinity JSON intermediates, none of which is retained in the
current clean checkout.  Rebuild the `EC-K3-H3-Q6` artifact chain first.  This
historical block therefore records the intended full-chain replay, not a
standalone clean-checkout command.

with endpoint

```text
Q8QQRR|ambient=13|rows=11|rank=11|kernel=2|...
Q8QQQUARTIC|degree=4
Q8QQCHILD|finite=[(1,[2,3,15],'I9*'),(9,[0,0,1],'I1')]|infinity=((0,0,0),'smooth')|root_rank=13|root_euler=24|root_det=4|MW_rank=4|status=PASS_EXACT_CORRECTED_Q8_D13_CHILD
```

Thus the second child is exactly `D13/MW4`.

## Two repaired bugs

### 1. Binary-quartic covariant 2-cover

The old marking treated differences of covariant images as primitive MW differences and then doubled them again. In fact the covariant map is the degree-two covering map to the Jacobian.

With

```text
Pmap = phi(E7_7)-phi(old_O),
Qmap = phi(E7_7)-phi(affine_E7),
```

the actual heights are

```text
height(Pmap)=32/3,
height(Qmap)=32/3,
<Pmap,Qmap>=4/3.
```

Therefore these are doubled primitive directions. The corrected q8 section is

```text
S=Pmap+Qmap,
MW=(-2,-2,0),
height=24,
S.O=10.
```

The withdrawn point was exactly `2*S`, with height `96` and collision degree `46`.

Canonical marking regression:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_marking.sage
```

which routes to `derive_h92_q6_child_q8_marking_2cover.sage`.

### 2. Missing `Dx` in q normalization

For

```text
x(S)=Nx/Dx,
y(S)=Ny/Dy,
p=-y(S)/x(S),
q=(m-p)/h,
```

the correct `Nx` pole cancellation is

```text
R*h*Dy == Ny*Dx mod Nx.
```

The old formula omitted `Dx`, leaving a degree-24 vertical pole. That produced generic modular branch degree about `100`; the corrected normalization gives branch degree `4` at every nonsingular level modulo both `43` and `59`.

## Corrected q8 geometry

```text
collision divisor degree h = 10
Nx degree = 24
Ny degree = 36
component-nef vertical fibre coefficient = -2
II* finite ideal = (u^2,X,Y)
IV* finite ideal = (u^2,X,Y)
```

After the corrected globally regular q frame, exact infinity orders imply

```text
deg(s) <= 6
deg(t) <= 5
ambient = 13.
```

Over `QQ` the complete matrix has rank `11`, hence `h0=2`.

## Source representative audit remains valid

The old source-side degree-18 and degree-16 classes are different Weyl representatives:

```text
dominant D13 hit --122 finite-root reflections--> degree 18
classifier-nef   --120 finite-root reflections--> degree 16.
```

The degree-16 terminal vector remains

```text
(22,16,-14,-20,-27,-40,-33,-26,-18,-4,-5,-7,-10,-8,-6,-4,-2,8,0).
```

This lattice result is exact, but the source-side `true1600` and hand-built `corrected1278` compilers are no longer the canonical equation route.

## Superseded q8 experiments

Treat these as historical diagnostics only:

- degree-46 marking/collision chain;
- old q-regular normalization without `Dx`;
- source `true1600 -> 18` pipeline as the final q8 pencil;
- experimental `corrected1278 -> 14 -> 7` q6^8 pipeline and identity-component tests.

## Next exact gate

Continue from the exact `D13/MW4` child toward the rootless/high-rank target. Preserve two regression rules in every later compiler:

1. binary-quartic covariant point differences carry the 2-cover multiplier;
2. derive base-pole residues from the fully cleared rational expression before CRT normalization—never drop section-coordinate denominators such as `Dx`.
