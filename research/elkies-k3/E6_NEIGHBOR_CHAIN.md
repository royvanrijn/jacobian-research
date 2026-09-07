# Exact neighbor chain from the rank-17 frame to E6/MW3

## Status

This note is an **exact lattice computation**.  It recovers the actual
neighbor path by which the stored `E6/MW3` frame was found; it is not yet a
birational or Weierstrass reconstruction of the corresponding K3 fibration.

Run:

```text
sage elkies-k3/scripts/verify_e6_neighbor_chain.sage
```

The verifier starts from the committed matrices, reconstructs every child
frame from its pinned primitive isotropic vector, and checks the composite
integral isometry against
[`data/fibrations/e6_ns_transport_from_rank17.txt`](data/fibrations/e6_ns_transport_from_rank17.txt).

## Recovered path

Write `N(M) = U + (-M)`, where `M` is a positive-definite rank-17 frame.
The exact path is

```text
rank17_gram
  -- q=90, (a,b)=(9,10) --> q90_mw7_frame
  -- q=4,  (a,b)=(2,2)  --> q90_mw4_frame
  -- q=4,  (a,b)=(2,2)  --> mw3_e6_a3a3_a1a1_frame.
```

For an input frame `M`, each arrow uses a vector `v` with
`v M v^T = 2q` and the primitive isotropic class

```text
f = (a,b,v) in U + (-M),       ab = q.
```

The three pinned vectors are:

```text
q=90: (0,0,0,0,0,0,-2,-1,0,0,6,-5,1,0,0,0,0)
q=4:  (-1,-1,-1,-3,0,1,1,-1,0,2,1,-1,-2,0,1,0,0)
q=4:  (0,0,-2,0,-2,3,-2,2,1,-1,-1,0,0,-1,-1,0,0)
```

For each `f`, the verifier deterministically constructs an isotropic mate
`g0` with `f.g0=1`, takes an integral basis `K` of their common orthogonal
complement, and forms

```text
T = [f; g0; K].
```

It checks `det(T)=+/-1` and

```text
T N(M_parent) T^T = N(M_child)
```

with the raw child matrix equal, entry for entry, to the committed frame.  If
the three matrices are `T1,T2,T3`, the committed composite is

```text
C = T3 T2 T1,
det(C) = 1,
C N(rank17_gram) C^T = N(E6frame).
```

The rows of `C` express the E6 Neron--Severi basis in the original rank-17
basis.  In particular, its first row is the E6 fiber class in the original
coordinates:

```text
(1458,1619,-10,1,-12,9,3,-1,-318,-181,10,1,978,-822,158,-9,1,2,3).
```

## What this corrects

The previous E6 equation search began from the abstract Kodaira configuration

```text
IV* + I4 + I4 + I2 + I2 + 4 I1
```

and imposed a convenient split normalization.  That chart was not obtained
by transporting the original K3 fibration through the neighbor path above.
Consequently, the exhaustive small-field result proves emptiness only for the
declared rational split chart with the imposed section ansatz; it does not
reject the exact E6 neighbor.

The explicit integral transport closes the lattice-provenance gap, but an
integral Neron--Severi isometry is not an explicit map of elliptic surfaces.
It does not by itself determine:

- the rational function giving the new elliptic parameter;
- the Weierstrass coefficients after each neighbor step;
- the fields of definition and Galois permutation of fiber locations and
  components;
- the rational representatives of the three final Mordell--Weil generators.

## Correct geometric backtrack

The next reconstruction must execute the recovered path geometrically, in
one direction or the other, from a genuine explicit model.  At each arrow it
must construct the pencil associated to `f`, compute its Jacobian/Weierstrass
model, and transport the zero section, reducible components, and section
classes.  Only the resulting final model may be normalized and fed to the
P1/P2/P3 equations, with the full height condition `(P1+P2).O=1` imposed.

Until an explicit starting model for the rank-17 K3 (or an independently
identified model at an intermediate node) is available, the lattice
certificate cannot manufacture those rational equations.  The current clean
anchor is the discriminant-3 CM model, but its former all-IV deformation is
parity-obstructed; see
[`E8_A2_KODAIRA_CORRECTION.md`](E8_A2_KODAIRA_CORRECTION.md).  Recovering the
correct non-isotrivial `I3/IV` deformation is therefore the precise missing
input.  Enlarging the old finite-field scan cannot supply it.
