# Determinant1020: arithmetic MW17 existence

**Subsequent construction:** the [explicit source certificate](DET1020_EXPLICIT_RATIONAL_SOURCE_2026-09-15.md)
now supplies the full split source equation and Picard19. The unresolved
source-equation statements below describe the earlier scope of this note;
physical MW17 transport and section compilation remain open.

The [fully rationally marked determinant1020 source](DET1020_RATIONAL_MARKING_SOURCE_2026-09-15.md)
admits an elliptic fibration over Q with Mordell–Weil group Z^17 over Q(t).
Its height lattice has determinant1020 and minimum4. This proves existence
on a different geometric NS from the determinant948 and1092 surfaces.
An explicit equation, a nef fiber in the source basis, and the seventeen
section coordinates remain uncomputed. This is not yet an experimental
surface for specialization or point search.

## Exact positive witness

In the positive even unimodular rank24 lattice with root system4E6,
embed one E6 component and a vector v perpendicular to it, with v²=340.
On the other three E6 components the Dynkin pairings of v are

```text
1 1 1 1 1 1
1 1 1 1 3 1
1 1 2 2 1 1
```

The retained [witness](../artifacts/generated-results/elkies-k3-det1020-root-gate-v1/4e6-probe.json)
fixes the actual ambient Gram, component bases, seven embedding rows and
seventeen complement rows. Thus component numbering and glue are explicit.
All seven Smith factors of the embedding are1. Its Gram is E6 + <340>.
The integral orthogonal complement F has determinant1020 and no norm2
vectors. Exact enumeration gives2524 norm4 vectors, so its minimum is4.
Positivity of the Dynkin labels also explains the construction: no root
in any of the remaining E6 components is perpendicular to v.

The final checker independently checks positive even unimodularity,
embedding primitivity, orthogonality, the saturated complement, its
Gram and the norm2 enumeration. It does not require completeness of an
auxiliary genus search or the discovery script's root-system labels.

## Identity with the admitted NS

Put T=[[-2,1,0],[1,2,0],[0,0,204]]. The discriminant groups of F and T
are cyclic of order1020. The certificate supplies explicit generators
with quadratic values

    q_F = 10877/1020,       q_T = -403/1020  (mod 2Z).

Multiplication by59 is a discriminant isometry: gcd(59,1020)=1 and
(q_F-59²q_T)/2 is integral. Hence U + F(-1) and the admitted
S=(-T)+E8(-1)^2 have the same signature(1,18) and discriminant form.
Both are even and indefinite. The discriminant length is1 and rank19;
Nikulin's uniqueness theorem applies, proving

    S ≅ U + F(-1).

The theorem is stated, with its length hypothesis, in
[Hassett–Tschinkel, Proposition2.2](https://cims.nyu.edu/~tschinke/papers/yuri/22autoeq/autoeq.pdf)
(citing Nikulin1.14.2). No explicit integral isometry to the source basis
has been computed. Such an isometry is needed for equation construction,
but not for this existence conclusion.

## Arithmetic fibration and saturation

The source proof gives a projective K3 X/Q with geometric NS exactly S
and every integral NS class represented by a line bundle over Q. Transfer
the displayed U abstractly to NS(X). Choose the effective sign of a
primitive isotropic generator and reflect into the nef chamber. The
result is a primitive nef class f of square0, still of divisibility1.
These operations preserve the integral marking; all resulting classes
are actual rational line bundles.

The complete pencil of f is defined over Q and has geometric base P1.
Since the line bundle itself descends, its two-dimensional space of
sections gives a map to P1_Q. The other U generator gives a class of
intersection1 with f; after adding a multiple of f it is effective.
Its horizontal part has degree1, giving a section over Q. Alternatively,
the usual U-to-elliptic-fibration theorem followed by uniqueness and full
rational marking gives the same rational section.

The frame f-perp/Zf is F(-1), with no roots, so every geometric fiber is
irreducible. Shioda–Tate gives geometric MW rank19-2=17. With a section
and no reducible fibers, the trivial lattice is the primitive unimodular
U; NS/U is free, so MW torsion is zero and its full saturated height
lattice is F. Every geometric section has a rational divisor class.
Its unique curve representative descends, and the degree-one map to
P1 makes it a Q(t)-point. Thus the full group Z^17 is rational over Q(t).
This is stronger than rationality of a finite-index subgroup.

## Replay and retained exploration

```sh
sage -python research/elkies-k3/scripts/certify_det1020_rootless_frame.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det1020-root-gate-v1/certificate.json)
contains the complete frame and finite-form witnesses. The arithmetic
admission checker is separate. Lattice uniqueness, nef reduction and
elliptic K3 theory are named written inputs; independent implementation,
formal verification, external review and literature novelty are unclaimed.

Discovery first tried a bounded p7 auxiliary genus exploration. PARI
stack failures in large-vector automorphism/isometry calculations were
retained. The root-splitting workaround and checkpoints are in the same
packet. The final bounded run stopped at118 classes and100000 neighbor
evaluations, with mass18329/7680 versus target123721/51840. It is not a
complete genus classification. The separate exact E6+<340> genus test
suggested the successful4E6 witness; no unfinished mass claim enters the
proof. No additional genus enumeration is needed for this positive gate.

The next gate is an explicit source equation at the intrinsic non-CM
period, followed by an integral transport, nef chamber calculation,
equation compilation and seventeen saturated sections. The overall
explicit foundry objective remains OPEN.
