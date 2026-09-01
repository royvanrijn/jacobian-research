# Golay-octad design of a rank-17 frame — 2026-09-01

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

## Literature boundary

The ambient construction is the standard Golay-code description of
`N(24A1)`, while the passage from a primitive Niemeier embedding to a K3 frame
is the Kneser--Nishiyama method.  The repository contribution here is the
specific septuple and its exact determinant, saturation, short-shell, and
ternary-realizability certificate.
