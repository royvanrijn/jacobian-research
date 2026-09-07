# Embed the generic discriminant-948 NS lattice into the discriminant-3 CM endpoint

The CM surface from Utsumi No.1 has

    NS_CM = U + E8(-1) + E8(-1) + A2(-1),

with determinant -3.

A primitive vector w of square -316 and divisibility 1 has orthogonal complement
of determinant

    3 * 316 = 948.

The script chooses an explicit norm-316 primitive vector in one E8 summand,
constructs S=w^perp, and compares its cyclic discriminant form with the recovered
generic target

    U + (-M_rank17).

If the discriminant forms agree, indefinite lattice uniqueness strongly supports
that S is the desired generic NS lattice.

Crucially, this embedding leaves the second E8 and the A2 untouched, but destroys
most/all roots in the first E8. Thus the useful deformation of Utsumi No.1 should
preserve only one II* fiber plus the IV/A2 structure, not both II* fibers.

Run:

    sage elkies-k3/scripts/embed_generic_ns_in_disc3_cm.sage

Key lines:

    CMEMBED|stage=generic_ns|...|det=...
    CMEMBED|stage=discform|...|unit_hits=...
    CMEMBED|stage=first_E8_roots|...
    CMEMBED|stage=fibration_prediction|...

The predicted MW rank is for the inherited elliptic fibration after removing the
extra CM class; it gives us the right complexity target for an explicit deformation.
