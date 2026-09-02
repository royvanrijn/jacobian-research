# Golay-octad design of a rank-17 frame — 2026-09-01

<!-- status-consumer: EC-K3-GOLAY-DET720-NS-SATURATION 306320cad1fd8e6e -->
<!-- status-consumer: EC-K3-GOLAY-DET720-PHYSICAL-CORRIDOR 0868df67fe8c37ad -->

## Outcome

Seven octads of the extended binary Golay code give a primitive rank-seven
auxiliary lattice in the `24A1` Niemeier lattice whose orthogonal complement
has the designed profile

```text
rank                         17
determinant                 720
minimum squared norm          4
signed norm-four vectors   3064
unoriented norm-four pairs 1532
```

This is an exact lattice construction, not merely a bounded search hit.  The
bounded septuple proposal that found it is not exhaustive, so no determinant
or short-shell optimality is claimed.

The frame discriminant group has Smith invariants `(2,6,60)`.  Among the
twenty even ternary genera of signature `(2,1)` and determinant `-720`, exactly
one has the required discriminant form.  One representative is

```text
[ 8 -2  2]
[-2 -4 10]
[ 2 10 -4].
```

Thus `NS = U + frame(-1)` passes the exact K3-lattice realizability gate.  A
general complex K3 with this Neron--Severi lattice has a rootless Jacobian
fibration of Mordell--Weil rank 17 and height lattice equal to the constructed
frame.  This does not supply a rational model, arithmetic marking, explicit
sections over `QQ(t)`, or any exceptional-specialization rank claim.

## Saturation and exact discriminant form

The two integral saturation statements both have index one.  The seven octad
vectors are primitive in `N(24A1)`, and their rank-17 orthogonal complement is
primitive as well.  The latter is also formal: it is the kernel of the
integral pairing map from the Niemeier lattice to the dual of the auxiliary,
so its quotient is torsion-free.  The independent coordinate-basis saturation
calculation returns index one.

For the K3 intersection lattice, distinguish the literal determinant from the
signed convention used in the elliptic-source ledgers:

```text
NS = U + frame(-1),       signature (1,18)
det Gram(NS) = +720,
disc_signed(NS) = -|det Gram(NS)| = -720.
```

Thus the requested `-720` is a signed discriminant, not the determinant of a
signature-`(1,18)` Gram matrix.  The discriminant group is

```text
A_NS = Z/2 + Z/6 + Z/60
     = (Z/2)^2 + Z/4 + (Z/3)^2 + Z/5.
```

In the canonical primary normal basis, the exact quadratic form
`q_NS : A_NS -> Q/2Z` is the negative of

```text
[1   1/2]  +  [1/4]  +  [2/3]  +  [4/3]  +  [4/5].
[1/2 1  ]
```

Equivalently, the displayed matrix is `q_frame=q_T`.  The rank-three lattice

```text
[ 8 -2  2]
[-2 -4 10]
[ 2 10 -4]
```

has signature `(2,1)`, determinant `-720`, and exactly that finite quadratic
module.  Hence `q_T=-q_NS`; gluing along this anti-isometry gives an even
unimodular lattice of signature `(3,19)`, so the embedding of `NS` in the K3
lattice is primitive.  Since the discriminant length is three and the
indefinite rank is nineteen, the usual rank-at-least-length-plus-two criterion
also makes this `NS` isometry class unique in its genus.  A generic period in
the displayed rank-three complement avoids every extra rational hyperplane,
and therefore has **exact** Neron--Severi lattice `NS` and Picard rank 19.

The exact replay is

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_ns_saturation.sage --check
```

with certificate
[`../artifacts/generated-results/elkies-k3-golay-det720-ns-saturation-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-ns-saturation-v1.json).

This theorem is geometric, but it does not make the K3 rational over `QQ`.
In fact the first exact rational `3I6` specialization is a sharp negative
control: it has rational 3-torsion and a rational half of the displayed
height-four section, hence index six and signed discriminant `-20`.  It is a
rational Picard-19 K3, but not the determinant-720 Golay K3.

## Physical corridor to the rootless frame

There is now an exact six-edge marking-level corridor inside the saturated
determinant-720 Neron--Severi lattice:

```text
3A5 -> 4A2+A5 -> 3A1+2A2+A3 -> 4A1+A2
    -> 3A1 -> 2A1 -> rootless MW17.
```

All six neighbours have old-fibre degree two and `q=4`.  Every fibre is
primitive; every finite/affine component, all-section, and finite horizontal
wall gate passes; and every edge carries an exact determinant-one transport.
The first five edges have horizontal `P.O. <= 2`.  The final `2A1`-to-rootless
edge is physical after two component-Weyl repairs, but this audited
presentation has `P.O.=4`.  A bounded rank-first search reaches `A1/MW16`
while keeping `P.O. <= 2`, and the tested order-eight automorphism group of
the displayed `2A1` bridge does not lower its final pole; neither computation
is an optimality proof.  Thus a physical corridor is proved, while the
stronger end-to-end `P.O. <= 2` corridor remains open.

Replay the selected witness with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_physical_corridor.sage --check
```

The certificate is
[`../artifacts/generated-results/elkies-k3-golay-det720-3a5-to-mw17-physical-corridor-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-3a5-to-mw17-physical-corridor-v1.json).
This is a corridor from the prescribed `G720-S0128` lattice marking, not from
the stored rational `s6=10` equation: determinant is invariant under these
unimodular transports, whereas that equation has saturated determinant 20.

## The seven-octad design

Use roots `e_i` with `e_i^2=2`.  For an octad `O`, put

```text
g_O = (1/2) * sum(e_i, i in O).
```

The selected supports are

```text
{1,13,15,17,18,19,23,24}
{3,4,16,18,20,21,22,24}
{4,5,6,10,14,20,21,22}
{1,3,4,8,12,14,16,24}
{5,6,13,15,17,19,20,22}
{2,4,8,10,11,12,14,21}
{3,5,6,7,9,16,20,22}.
```

Every support is a weight-eight word of `G_24`, and the seven supports cover
all 24 coordinates.  Their Gram matrix follows directly from

```text
(g_O,g_O') = |O intersect O'|/2
```

and is

```text
[4 1 0 1 2 0 0]
[1 4 2 2 1 1 2]
[0 2 4 1 2 2 2]
[1 2 1 4 0 2 1]
[2 1 2 0 4 0 2]
[0 1 2 2 0 4 0]
[0 2 2 1 2 0 4].
```

Its determinant is `720`.  Exact saturation in `N(24A1)` has index one.
Because the seven supports cover every coordinate, no coordinate root is
orthogonal to the auxiliary.  The full short-vector enumeration independently
confirms that the complement has no norm-two vectors and has exactly 3,064
norm-four vectors.  A second count using only the Construction-A description
splits those vectors as

```text
24   signed vectors of coordinate-pair type
3040 signed vectors supported with signs on a Golay octad
3064 total.
```

This also explains the proper role of the coding layer.  Octad intersections
are an excellent proposal language, but the acceptance gate must still use the
full integral Niemeier lattice: primitive closure, complement determinant,
local discriminant form, and exact short-shell enumeration are not determined
by a visual support pattern alone.

## Relation to the determinant-948 H3 problem

This construction varies the Neron--Severi discriminant.  It does not produce
a third frame on the fixed H3 surface.  The completed determinant-948 `J2`
classification proves that both rootless H3 frames occur in
`N(2A7+2D5)` and that no rootless determinant-948 target-genus complement
occurs in `N(24A1)`.

The Golay construction is therefore a foundry for new Picard-19 lattice
classes, not a shortcut around the fixed-H3 classification.  Compared with
the earlier `24A1` determinant-952 near-miss (1,313 norm-four pairs) and the
broader determinant-950 foundry target (1,322 pairs), the present octad design
has both smaller determinant and a substantially larger minimal shell: 1,532
pairs.

## Replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_octad_rank17_design.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_octad_rank17_design.sage --check
```

The deterministic certificate is
[`../artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json`](../artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json).

## Immediate source-first gate

The next required test was run on the same determinant-720 frame genus.  A
deterministic eighteen-generation Kneser-neighbour beam visited 34,101 exact
reduced classes.  It found no root-rank-15--17 source in that declared window,
so it did not reach the preferred MW0--2 band.  Its best primitive-root source
is nevertheless

```text
root type      2A1 + A2 + 3A3
root rank      13
MW rank         4
MW regulator   15/16
```

The exact MW height Gram is

```text
[ 5/6  5/12 -1/3   0]
[5/12   5/6 1/12   0]
[-1/3  1/12  5/6   0]
[   0     0    0   3].
```

The root lattice is primitive, and
`det(root lattice) * regulator = 768 * 15/16 = 720` as required.  The all-`A`
semistable fibre profile would be

```text
3I4 + I3 + 2I2 + 5I1.
```

It has six reducible-fibre supports.  This is a viable MW4 source, but it is
less equation-friendly than the foundry's preferred two/three-support MW2
targets.

Because the discriminant module is noncyclic with invariants `(2,6,60)`, the
cyclic-unit glue shortcut is invalid.  Exact Smith-generator enumeration gives
96 anti-isometries.  The selected three-generator graph glue embeds the same
rank-seven auxiliary primitively into

```text
N(4A5 + D4),
```

and its saturated orthogonal complement is integrally isometric to the
displayed `2A1+A2+3A3/MW4` source.  Thus the companion is a certified
same-auxiliary Kneser--Nishiyama source, not merely a genus mate.

Replay the bounded source search and the noncyclic Niemeier gate with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/hunt_golay_octad_rootful_source.sage \
  --allow-below-target

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_octad_noncyclic_niemeier_source.sage
```

The exact outputs are
[`../artifacts/generated-results/elkies-k3-golay-octad-det720-source-hunt.json`](../artifacts/generated-results/elkies-k3-golay-octad-det720-source-hunt.json)
and
[`../artifacts/generated-results/elkies-k3-golay-octad-det720-source-niemeier.json`](../artifacts/generated-results/elkies-k3-golay-octad-det720-source-niemeier.json).

The correct next search, if this class is retained, is a direct full-Niemeier
prescribed-support enumeration for a rank-15 two/three-support companion.  A
longer undirected Kneser beam is lower-value: the present beam already shows
that root rank grows readily to 13 but not to the preferred source band.

## Direct all-Niemeier source closeout — 2026-09-02

The prescribed-support search has now been run directly for the fixed
determinant-720 auxiliary.  The auxiliary is rootless, so the older `D5`-anchor
enumerator does not apply.  Instead, its ordered norm-four octad basis is
embedded one vector at a time.  At each prefix the next vector is placed in a
dominant chamber for the residual Weyl group; at the final vector, zero Dynkin
labels prescribe root rank 15 or 16 and one to three root supports before the
remaining exact shifted ellipsoid is solved.

The unrestricted run covers all 23 rooted Niemeier lattices and has no
candidate or label limit.  It retains 4,823 deterministic reduced-Gram rows:

```text
MW rank 1       32 rows
MW rank 2     4791 rows
total         4823 rows
```

These are stored rows, not an assertion of 4,823 integral-isometry classes.
Equal deterministic reduced Grams are merged, but diagram/umbral
automorphisms are not quotiented and a general rank-17 isometry classification
is not used as a discovery gate.

The separate pole audit resolves primitive and nonprimitive root lattices
uniformly.  It computes the exact Smith quotient by the roots, enumerates
every torsion/free MW class whose exact height lower bound can permit frame
norm at most eight, and repeats each affine closest-vector decision with
double-double and MPFR-256 GSO arithmetic.  Returned norms, root pairings,
Smith classes, and unimodular MW bases are checked exactly.  The strict result
is

```text
complete MW basis with every P.O. <= 2       1587 rows
  MW1                                           16 rows
  MW2                                         1571 rows
```

The useful Pareto leaders are:

| representative | root type | MW | supports | all-`A` | complete basis P.O. | MW height Gram |
|---|---|---:|---:|---:|---:|---|
| `G720-S0046` | `A4+A8+D4` | 1 | 3 | no | `[1]` | `[[4]]` |
| `G720-S0224` | `A1+A11+D4` | 1 | 3 | no | `[2]` | `[[15/2]]` |
| `G720-S0260` | `A11+A4` | 2 | 2 | yes | `[0,0]` | `[[4,-2],[-2,4]]` |
| `G720-S0780` | `A2+D13` | 2 | 2 | no | `[2,2]` | `[[8,-6],[-6,12]]` |
| `G720-S0422` | `A1+A4+D10` | 2 | 3 | no | `[1,1]` | `[[9/2,-3/2],[-3/2,9/2]]` |

All five displayed rows have primitive roots and trivial torsion.  The
equation-shape leader is `G720-S0260`: geometrically its semistable fibre
profile is

```text
I12 + I5 + 7I1,
```

its two displayed pole-zero sections give a unimodular basis of the free MW
quotient, and `det(A11)*det(A4)*det(MW) = 12*5*12 = 720`.
Four stored `A11+A4/[0,0]` rows occur in `N(A15+D9)`; the table names the
first deterministic representative.  The MW1 leader instead has only one
section condition to impose, but trades the second support for a `D4` fibre.

This clears the requested lattice gate and is strictly more attractive than
the older NS0024 `A3+A4+A6/MW4` three-support baseline on MW rank, support
count, and complete-basis pole cost.  It is not yet a construction-level
replacement: NS0024 has a certified marked degree-two corridor, while the
Golay-720 rows still have no rational source marking, Weierstrass equation, or
marked route to the rootless target.

Several negative and near-miss slices remain useful routing information:

- `N(D24)` gives 23 `D15/MW2` rows but no non-torsion section through pole
  two.
- `N(A24)` gives 90 one-support `A15/MW2` rows, none with a complete basis
  through pole two.
- `N(4A6)` and `N(4E6)` give no retained source in the declared window.
- `N(4A5+D4)` gives 59 three-support `3A5/MW2` rows with 3-torsion.  Each has
  a pole-zero non-torsion direction, but the cheap sections span only one
  free-MW line, so this known glue home is a near miss rather than the winner.
- The older 34,101-class Kneser beam therefore failed because it did not enter
  the required root-rank band; direct prescribed-root enumeration was the
  decisive change.

### First equation gate for `G720-S0260`

The normalized semistable fibre ansatz has also passed its first finite-field
gate.  Over `GF(5)`, fix `I12` at zero, `I5` at infinity, and normalize
`A(0)=-3` in

```text
y^2 = x^3 + A(t)x + B(t),  deg(A)<=8, deg(B)<=12.
```

An exhaustive scan of all `5^8 = 390625` normalized `A` polynomials solves
the rank-13 Hermite system for the 17 prescribed discriminant jets.  Of 916
compatible signed local branches, 208 have the exact orders `(12,5)` and 164
have squarefree residual degree-seven discriminant.  Thus the bare
`I12+I5+7I1` equation stratum is nonempty and abundant modulo 5.

This is encouraging but is not yet the determinant-720 source equation.  The
two pole-zero sections, the full lattice marking, a characteristic-zero lift,
and a rational parameter remain open.  Replay the exact modular census with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_golay_det720_a11_a4_source_ansatz_modp.sage \
  --prime 5 --check
```

The artifact is
[`../artifacts/generated-results/elkies-k3-golay-det720-a11-a4-source-ansatz-mod5-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-a11-a4-source-ansatz-mod5-v1.json),
SHA-256
`52fc4cf278bb0d8631ec9de32c82603af3bd9e7f4ffab7b5a12fc473fe77e17e`.

Replay the source enumeration and its independent pole audit with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_golay_det720_prescribed_root_sources.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_golay_det720_source_poles.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_golay_det720_source_poles.sage --check
```

The generated artifacts are
[`../artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json`](../artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json)
with SHA-256
`f203027cae98df3f1cf69286dd149a73ef93f35ee9408130d89e39a856dbb7af`,
and
[`../artifacts/generated-results/elkies-k3-golay-octad-det720-source-poles-v1.json`](../artifacts/generated-results/elkies-k3-golay-octad-det720-source-poles-v1.json)
with SHA-256
`ddeda8fd0059ec7b1e45f09d87be1a0dd6815152ee12bcd2bb3fb313e6e46bf8`.
The 1587-row success count is a complete-basis statement; the source
artifact's larger minimum-section count is the weaker metric and should not be
substituted for it.

## Literature boundary

The ambient construction is the standard Golay-code description of
`N(24A1)`, while the passage from a primitive Niemeier embedding to a K3 frame
is the Kneser--Nishiyama method.  The repository contribution here is the
specific septuple and its exact determinant, saturation, short-shell, and
ternary-realizability certificate.

Primary context is Nishiyama's original K3/Niemeier construction
([J. Math. Kyoto Univ. 22 (1982), 293--304](https://www.jstage.jst.go.jp/article/math1924/22/2/22_2_293/_article)),
the lattice-polarized K3 review of Braun--Kimura--Watari
([arXiv:1312.4421](https://arxiv.org/abs/1312.4421)), and Shimada's exact
elliptic-K3 ADE/Mordell--Weil classification framework
([arXiv:math/0505140](https://arxiv.org/abs/math/0505140)).
