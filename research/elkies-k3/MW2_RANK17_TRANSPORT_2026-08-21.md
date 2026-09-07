# Exact transport from the reconstructed MW2 frame to rank 17 (2026-08-21)

## Status

This is an **exact integral Neron--Severi lattice certificate**.  It composes
the four pinned neighbors

```text
rank17 --q=25--> MW7 --q=4--> MW4 --q=4--> MW3 --q=4--> MW2
```

and expresses the result in a terminal basis made from the fiber, zero
section, reducible-fiber components, and the two reconstructed section
classes.  It supplies the divisor classes required to begin the inverse
geometric neighbor construction.  It does not by itself construct the
corresponding pencils or Weierstrass coordinate changes.

Run

```bash
sage elkies-k3/scripts/verify_mw2_rank17_transport.sage
```

to reconstruct all four child frames, all four unimodular transition
matrices, the terminal section lifts, and the composite certificate from the
pinned witnesses.

## Conventions and terminal explicit basis

For a positive frame Gram `M`, write

```text
N(M) = U + (-M),       U = [0 1; 1 0].
```

In the standard terminal basis the fiber and isotropic mate are `e0,e1`, and
the zero section is `O=e1-e0`.  The pinned explicit basis `B` has row order

```text
F, O,
E6[1..6], D4[1..4], A2_P2[1..2], A2_P1[1..2], A1[1],
P1, P2.
```

It is stored in
[`data/fibrations/mw2_e6_d4_a2a2_a1_explicit_basis.txt`](data/fibrations/mw2_e6_d4_a2a2_a1_explicit_basis.txt).
Its rows give the explicit basis in the standard
`U+(-mw2_e6_d4_a2a2_a1_frame)` basis, and `det(B)=-1`.

The roots are deterministic simple-root systems.  Their positive Grams in
the displayed row order are

```text
E6 = [ 2 -1  0  0  0  0]    D4 = [ 2 -1  0  0]
     [-1  2  0 -1  0  0]         [-1  2 -1 -1]
     [ 0  0  2  0 -1  0]         [ 0 -1  2  0]
     [ 0 -1  0  2 -1 -1]         [ 0 -1  0  2]
     [ 0  0 -1 -1  2  0]
     [ 0  0  0 -1  0  2].
```

Both `A2` blocks have Gram `[2,-1;-1,2]`, and the `A1` block has Gram
`[2]`.  With the negative-definite signs in `NS`, the section/root
intersections are

```text
P1: E6[1], A2_P1[2], A1[1]
P2: E6[3], D4[3], A2_P2[2].
```

All other section/root intersections vanish.  Moreover

```text
P1.O=0, P2.O=1, P1.P2=2,
height Gram = (1/6) [9 2; 2 18].
```

This identifies the basis with the reconstructed profiles
`P1=(1,0,0,1,1)` and `P2=(2,d1,1,0,0)`, up to the independent Dynkin
diagram orientations and simultaneous section negation.  The opposite
section orientation is induced on the Weierstrass model by `(x,y)->(x,-y)`.

## Composite transport

For each neighbor the verifier constructs

```text
T_i = [f_i; g_i; K_i],
```

where `f_i` is the pinned primitive isotropic class, `g_i` is the
deterministically constructed isotropic mate, and `K_i` is an integral basis
of their common orthogonal complement.  It checks `det(T_i)=1` and the raw
child frame entry by entry.

Let

```text
C = T4 T3 T2 T1.
```

Then `det(C)=1` and

```text
C N(rank17) C^T = N(MW2).
```

The requested explicit composite is

```text
T = B C.
```

Its rows express the displayed terminal trivial-plus-sections basis in the
original `U+(-rank17_gram)` basis.  It is integral, has determinant `-1`, and
is stored in
[`data/fibrations/mw2_e6_d4_a2a2_a1_ns_transport_to_rank17.txt`](data/fibrations/mw2_e6_d4_a2a2_a1_ns_transport_to_rank17.txt).

## Old rank-17 classes in the terminal explicit basis

Inverting `T` over the integers gives the original rank-17 fiber

```text
F17 =
(322560,21927,63147,79421,9900,104071,55990,48836,7502,
 15891,1148,23070,-427,833,8390,15505,15422,24189,-30612).
```

The original isotropic `U`-mate is

```text
G17 =
(322604,21930,63156,79432,9901,104085,55997,48843,7503,
 15893,1148,23073,-427,833,8391,15507,15424,24192,-30616).
```

Consequently the original zero section `O17=G17-F17` is the much smaller
class

```text
O17 = (44,3,9,11,1,14,7,7,1,2,0,3,0,0,1,2,2,3,-4).
```

Exact pairing checks give

```text
F17^2=G17^2=0, F17.G17=1,
O17^2=-2, F17.O17=1.
```

## Stagewise inverse-neighbor pencils

For direct execution of the inverse `q=4,4,4,25` path, the verifier also
inverts every partial composite.  All coordinates below use the same
terminal explicit basis.

```text
previous MW3 fiber =
(42,3,8,10,1,13,7,6,1,2,0,3,0,0,1,2,2,3,-4)
previous MW3 mate  =
(42,3,8,10,2,13,7,6,1,2,0,3,0,0,1,2,2,3,-4)
previous MW3 zero  =
(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

previous MW4 fiber =
(811,55,159,200,25,262,141,123,19,40,3,58,-1,2,21,39,39,61,-77)
previous MW4 mate  =
(811,55,159,200,25,262,141,123,19,40,3,58,-1,2,22,40,39,61,-77)
previous MW4 zero  =
(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)

previous MW7 fiber =
(15375,1045,3010,3786,472,4961,2669,2328,358,758,55,1100,
 -20,40,400,739,735,1153,-1459)
previous MW7 mate  =
(15373,1045,3010,3786,472,4961,2669,2328,358,758,55,1100,
 -21,39,400,739,735,1153,-1459)
previous MW7 zero  =
(-2,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0).
```

The final stage is `(F17,G17,O17)` above.  A machine-readable copy of all
twelve vectors is
[`data/fibrations/mw2_e6_d4_a2a2_a1_inverse_neighbor_classes.tsv`](data/fibrations/mw2_e6_d4_a2a2_a1_inverse_neighbor_classes.tsv).

For every stage the verifier checks `F^2=G^2=0`, `F.G=1`,
`O=G-F`, `O^2=-2`, and `F.O=1` against the terminal explicit Gram.

## Pinned hashes

```text
explicit basis:
4e2bb7f755c8769978e7a96e65e2bfd638134d780b1f8f103b1643229e0059a8

explicit terminal-to-rank17 transport:
c712135c3695d045584b1c27f0f5a5bec360cc42b6bbf22a88991b58c388f463

inverse-neighbor class table:
3e6cdb4307bb64504d1d449128fca82a0a37b3545d665cba3ea012f1b9e91717
```
