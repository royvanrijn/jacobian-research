# Genuine quadratic neighbours of the `E6+A1` Picard-rank-19 K3 — 2026-09-02

<!-- status-consumer: EC-K3-E6A1-RHO19-GENUINE-Q2-MW3 cd4314040bb028f7 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-EQUATION 8cfa9387612ac443 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-ARITHMETIC-RANK2 2130bc147519ac6b -->

## Result

The first genuine quadratic-neighbour shell of the new determinant-36 K3 is
complete.  In the split marking

```text
NS(X) = U + M(-1),     F=e,     O=f-e,
```

the correct zero-neutral classes are

```text
D = 2e+2f-w,     w.M.w=8.
```

There are `119` Weyl-dominant norm-eight classes.  Three are divisible by two,
so they are degree-one pencils written twice.  Of the remaining `116`
primitive classes, physical chamber reduction leaves `90` at old-fibre degree
two, sends `14` to degree one, and sends `12` to degree zero.

Exactly `18` of the genuine degree-two classes have child root rank `14`.
Every one passes a complete nefness proof and therefore defines a Jacobian
elliptic fibration of Mordell--Weil rank three on the same generic
Picard-rank-19 K3.  The four resulting frame types are

| child roots | number | torsion | representative free MW height lattice | regulator |
|---|---:|---:|---|---:|
| `4A1+2D5` | 4 | `Z/2` | `[[1/2,-1/4,-1/4],[-1/4,5/4,1/2],[-1/4,1/2,5/4]]` | `9/16` |
| `2A2+2D5` | 4 | trivial | `[[1/6,-1/12,-1/12],[-1/12,17/12,2/3],[-1/12,2/3,17/12]]` | `1/4` |
| `2A1+D6+E6` | 4 | trivial | `[[2/3,-1/6,-1/6],[-1/6,7/6,-1/3],[-1/6,-1/3,7/6]]` | `3/4` |
| `A7+D7` | 6 | trivial | `diag(3/8,1,3)` | `9/8` |

The exact replay is
[`scripts/certify_e6a1_rho19_genuine_q2_neighbors.sage`](scripts/certify_e6a1_rho19_genuine_q2_neighbors.sage),
with generated certificate
[`../artifacts/generated-results/elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json`](../artifacts/generated-results/elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json).

## Why norm four and norm six were false starts

For a split class

```text
D=a*e+b*f+w
```

one has `D.O=a-b`.  If `a<b`, the old zero is a fixed component.  Removing
`(b-a)O` exchanges the two hyperbolic coefficients:

```text
D-(b-a)O = b*e+a*f+w.
```

Thus `e+2f+w`, `w^2=4`, reduces to degree one, as does
`e+3f+w`, `w^2=6`.  They are not quadratic or cubic neighbours.  The first
presentation with nonnegative zero intersection at old degree `q` has
`a=q`, hence norm `2q^2`; this gives norm eight for `q=2` and norm eighteen
for `q=3`.

This correction is now checked in the dissection certificate and recorded as
a reusable chamber lemma in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).

## Complete orbit enumeration

The norm-eight shell is not enumerated as a large signed vector list.  The
current root system is `2E6+A3`, of rank `15` in the rank-17 frame.  A
Weyl-dominant vector is specified by nonnegative Dynkin labels and a vector in
the rank-two orthogonal kernel.

The exact component label counts through dual norm eight are

```text
A3: 20,     E6: 17,     E6: 17.
```

Combining them leaves `553` labels of total projected norm at most eight.
The simple-root pairing matrix has Smith invariants

```text
1^14, 6,
```

and its integral kernel has Gram matrix

```text
[ 12 -18 ]
[-18  30 ].
```

An exact closest-vector enumeration in this rank-two lattice produces the
`119` dominant norm-eight vectors.  This is a complete Weyl-orbit census, not
a bounded coordinate-box search.

The three nonprimitive classes have the form

```text
2(e+f-r),     r^2=2,
```

and are discarded before splitting the child hyperbolic plane.

## Physical chamber and child frames

Each primitive class is tested against the old zero, all fifteen nonidentity
components, and the three affine identity components of the two `IV*` fibres
and the `I4` fibre.  Deterministic fixed-component removal gives

```text
old degree 2: 90,
old degree 1: 14,
old degree 0: 12.
```

For every primitive class, a unimodular hyperbolic split constructs the child
rank-17 frame.  Exact norm-two enumeration gives the complete child root
system.  The full distribution over all `116` primitive classes is serialized
in the certificate.  The eighteen root-rank-14 classes all remain at degree
two with an empty physical reduction sequence.

Root saturation is part of the computation.  The `4A1+2D5` root lattice has
index two in its primitive closure, giving child torsion `Z/2`; the other
three types are primitive and have trivial torsion.  The displayed regulators
satisfy exactly

```text
disc(NS) = disc(R) * Reg(MW) / |MW_tors|^2 = 36.
```

## Complete nefness proof

After the physical component gate, a negative irreducible `(-2)`-curve can
have old-fibre degree only one or two.

Degree one curves are old sections.  Write their Mordell--Weil coordinates as
`n=(a,b)` in the old height lattice

```text
H = diag(1/3,3).
```

For each candidate the exact MW projection gives an integral centre `delta`.
Proposition C2 implies that a negative section must lie in the strict finite
ellipsoid

```text
(n-delta/2)^T H (n-delta/2) < 2.
```

There are only five, six, or eight such sections per candidate.  Their
integral divisor classes are reconstructed with the exact `IV*` and `I4`
Shioda corrections, and every intersection is nonnegative.

For degree two, suppose an irreducible root `C` had `D.C<0`.  It would be a
fixed component, so `V=D-C` would be effective and vertical.  Since the
vertical intersection form is negative semidefinite,

```text
V^2 = -2-2*D.C >= 0
```

forces `D.C=-1` and `V^2=0`.  Hence `V=lF`.  Primitivity of `F` makes `l` an
integer, while

```text
C^2=(D-lF)^2=-4l=-2
```

forces `l=1/2`, a contradiction.  This excludes every remaining horizontal
root and proves global nefness.  The unimodular `U` split then supplies an
effective degree-one `(-2)` class, so each pencil is Jacobian.

## Compiler targets

The sparse Pareto winner is orbit `103`:

```text
D = P0 + P1 + A3_2.
```

Its old-basis coefficient maximum and `L1` norm are `(1,3)`, its child root
system is `2A2+2D5`, and its free MW regulator is `1/4`.  This is the preferred
resolved Riemann--Roch target because the pencil is built from two already
explicit sections and one old fibre component.

This target has now been compiled exactly.  Its resolved basis is `<1,z>` and
its clean Weierstrass model has fibres `2I1*+2I3+4I1`; see
[`E6A1_RHO19_ORBIT103_WEIERSTRASS_2026-09-02.md`](E6A1_RHO19_ORBIT103_WEIERSTRASS_2026-09-02.md).
Its arithmetic descent is also complete: the two displayed `QQ(k)(r)` points
are independent, while an explicit third geometric direction is
anti-invariant over `QQ(sqrt(-3))`.  Thus its generic arithmetic MW rank is
exactly two, not three.

If the number of child reducible fibres matters more than source-divisor
sparsity, orbit `96` is the secondary target.  It has roots `A7+D7`, MW height
`diag(3/8,1,3)`, and old-basis complexity `(2,12)`.

The abstract `E6` root basis used by this census has not yet been matched to a
complete resolved physical `E6` marking on the equation.  The obvious
equation-side `P0` tangent trace has been compiled as a fail-closed test: its
fibres are `2IV*+I4+4I1`, with roots `2E6+A3`, so it is not orbit `96` and
cannot be used as a low-coefficient `A7+D7` competitor.  The declared Pareto
comparison is currently exact only for source-divisor complexity and fibre
count.  Weierstrass coefficient optimality requires the missing physical
marking and a genuine `A7+D7` compilation.

## Reproduction and boundary

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_genuine_q2_neighbors.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_genuine_q2_neighbors.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage --check
```

Pinned SHA-256 values are

```text
checker  e7f8eb8a336fba9ee457245816a6ab11c58de17e1939ee011234c96b0a7c1428
artifact 07a6cdb181e9c2febdfcf085954608ef40b07a654fc652624bae195e4128d179
descent/audit checker  23988e2a5e608301ebf393fb157c0c03a02ae353760aadf1fbe5ed5a689304c1
descent/audit artifact 1daae6cc371f04b541ac1127748db75924098bc959c30ffbe067bbd3027a41d9
```

The orbit census, physical reductions, child root systems, torsion, free MW
height lattices, and nefness of the eighteen rank-three frames are exact.  The
orbit-103 equation and arithmetic rank-two descent are exact in the dependent
certificates.  The naive orbit-96 tangent is exactly rejected, but a genuine
physically marked `A7+D7` equation and coefficient comparison remain open.
