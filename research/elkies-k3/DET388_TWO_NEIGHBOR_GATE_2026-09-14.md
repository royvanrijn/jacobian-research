# Determinant 388: one complete rootless 2-neighbour obstruction

For the pinned five-A1 frame at index106 of the
[bounded probe](../artifacts/generated-results/elkies-k3-det388-rootless-probe-v1/checkpoint.json),
**every even integral 2-neighbour contains a root**. This excludes a
one-step rootless route from this particular frame, not a rootless frame
elsewhere in its genus and not an MW17 fibration on the abstract
[fully rational determinant388 K3 source](DET388_RATIONAL_MARKING_SOURCE_2026-09-14.md).

The exact finite cover has 1792 eligible parity classes. Each has an
explicit norm-eight witness. The certificate therefore proves a negative
one-step statement without relying on the completeness of a neighbour walk.

## Why the parity certificate covers every even 2-neighbour

Write F for the pinned positive frame. An integral 2-neighbour N has
K=F intersect N of index two in both. Choose x in N outside F and put
v=2x in F. The parity of v is nonzero. If (v,F) is not contained in 2Z,
integrality forces

    K = {y in F : (v,y) is even}.

All five displayed roots must be outside K for N to be rootless, so
(v,r_i) must be odd for every i. Evenness of N requires v² divisible by
eight. In particular q(v)=v²/2 is even; this last condition depends only
on v modulo2F. Exhausting all 2^17 parity vectors with these six conditions
gives exactly1792 candidates.

The bad-prime exception is checked, not assumed away. The radical of
F/2F is one-dimensional. Its nonzero binary representative has norm4
modulo8. Every integral lift of this radical vector has the same norm
modulo8, because its pairings are even and F is even. Consequently v/2
would have odd norm, which is impossible in N. Thus no even 2-neighbour
can arise from the radical, and the preceding parametrization is exhaustive.

For each surviving parity a, the certificate gives w in F with
w congruent to a modulo2F and w²=8. Let v be the even-neighbour lift of a.
Write w=v+2z. Since w² and v² are divisible by eight,

    w²-v² = 4(v,z)+4z²

implies (v,z) even. Therefore z is in K and w/2=v/2+z lies in N.
Its norm is two: it is a root. One such witness for each parity excludes
all1792 candidates. Neighbours not satisfying the five odd-pairing
conditions retain a displayed old root.

## Evidence and boundary

The [certificate](../artifacts/generated-results/elkies-k3-det388-two-neighbor-gate-v1/certificate.json)
contains the full frame Gram, five orthogonal roots, the complete parity
list and1792 norm-eight witnesses. The checker also verifies the retained
rational neighbour transports connecting this frame to the admitted NS
source, verifies the initial integral change of basis, and checks equality
of the full local genera. The stabilized lattices U+(-F) have the same
signature and discriminant form; indefinite lattice uniqueness at rank19
and discriminant length one identifies the NS lattice. It checks the
mod2 radical and every witness exactly. This replay needs no short-vector
search and grants no global genus classification.

The producer enumerated806586 signed nonzero vectors of norm at most eight
in about six seconds. That enumeration found the witnesses; its completeness
is not required for this exclusion. A supplementary shell census in
`selected-child.json` has minimum five root pairs among the1792 root-removing
neighbours. It provides a same-genus child and remains supporting search data;
the registered theorem is the witness-based one-step exclusion above.

The previous 128-neighbour miss is retained separately. No new equation,
rootless frame, global maximum MW rank, or arithmetic exclusion follows.
The subsequent [global root-energy theorem](DET388_GLOBAL_ROOT_OBSTRUCTION_2026-09-14.md)
now excludes rootless existence in this NS. That stronger conclusion does
not follow from the one-frame gate alone. Revisiting
this same one-step search without a changed frame cannot produce a hit;
a further route needs another frame or a genuinely global lattice gate.

```sh
sage -python research/elkies-k3/scripts/check_det388_two_neighbor_gate.py
```
