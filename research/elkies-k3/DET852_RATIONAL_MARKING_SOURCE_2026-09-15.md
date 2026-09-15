# Determinant 852: full saturated rational NS marking

The subsequent [global root obstruction](DET852_GLOBAL_ROOT_OBSTRUCTION_2026-09-15.md)
now excludes every MW17 fibration on this NS. Its full rational marking
and non-CM period remain valid; any rootless UNKNOWN below belongs to
this earlier gate alone.

There exists a projective K3 surface X/Q of geometric Picard rank19 whose
full saturated geometric Neron–Severi lattice is represented by rational
divisor classes and has determinant852. Its transcendental lattice is

```text
-2 0   1
 0 4   0
 1 0 106
```

This admits K3-4ff75fec54d01662 to the geometric rootless-frame gate.
It does not yet give a K3 equation or an MW17 fibration. The period is
specified intrinsically; its elliptic-model coordinates remain uncomputed.

## Primitive integral marking

Let F be the positive rank17 Gram of retained source frame F001 for this
surface, and let S=U+(-F). The certificate stores the entire matrix; its
signature is(1,18) and determinant852. Both S and the displayed T have
cyclic discriminant group of order852. Their exact quadratic generators
admit a checked anti-isometry. The checker constructs the graph overlattice
of S+T of index852 and verifies that its Gram is even, integral, unimodular
and of signature(3,19), hence the K3 lattice.

The graph projections are injective on both discriminant groups, so S
and T remain primitive and are actual orthogonal complements. This is not
an unsaturated divisor subgroup or a determinant-only identification.

## Non-CM period and effective surface descent

The [full stable marked curve](DET852_FULL_MARKED_CURVE_2026-09-15.md)
is C=X_0^6(71)/<w426> over Q. Its complete arithmetic level, including all
local units and the determinant-minus-one stable reflection, is certified.
The [intrinsic non-CM construction](DET852_INTRINSIC_NONCM_POINT_2026-09-15.md)
provides the following particular rational period. Choose a rational
CM(-67) point O. Take the unique rational CM(-163) point R such that R-O
is divisible by two in Jac(C)(Q). The point O+2(R-O) is rational non-CM.
Every allowed choice of O works; numerical coordinates are unnecessary
for this existence/descent statement.

Choose an ample chamber in S. The canonical arithmetic period map and
Torelli identify its lattice-polarized moduli component with the appropriate
open part of C, with the same explicit theorem inputs as the
[determinant388 descent proof](DET388_RATIONAL_MARKING_SOURCE_2026-09-14.md#3-the-actual-primitive-ns-and-descent-to-q).
For a ternary transcendental lattice, extra rational algebraic classes
force a CM period. Thus the specified non-CM point lies in the moduli image
and its geometric NS is exactly S.

An automorphism fixing S fixes an ample class and has finite order. On
rank-three T, any finite Hodge isometry other than +I or -I leaves a
nonzero rational eigenspace orthogonal to the period plane, giving an extra
algebraic class. This contradicts NS=S. The alternative -I is incompatible
with the identity on S under primitive gluing: it acts nontrivially on the
cyclic discriminant group of exponent852. Therefore a marked automorphism
acts identically on all H² and is the identity by Torelli. The non-CM
marked object has trivial inertia; its unique marked isomorphisms give an
effective descent datum at the rational moduli point. This yields X/Q
with every class of S Galois invariant.

## Invariant classes descend to actual rational divisors

Write e,f for the U basis, e.f=1, and w for the fifth frame coordinate in
-F. It has square-6. The two integral classes

    r=e-f,        s=e+2f+w

satisfy r²=s²=-2 and r.s=1, checked against the full NS Gram.
K3 Riemann–Roch gives chi(r)=chi(s)=1. For an invariant geometric line
bundle, the index of its Brauer descent obstruction divides its Euler
characteristic, since the cohomology spaces are twisted descent spaces.
Both r and s therefore descend. Their intersection is a zero-cycle class
of degree one. Restriction and corestriction imply that Br(Q)->Br(X) is
injective. The Picard-Brauer exact sequence now makes every invariant
geometric line bundle descend. Thus all saturated NS classes are actual
rational divisor classes, not merely Galois-invariant classes.

## Verified boundary and next step

```sh
sage -python research/elkies-k3/scripts/certify_det852_rational_marking.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det852-marking-v1/certificate.json)
contains the primitive glue, full NS/T Grams and descent intersection pair.
It imports the already-certified full curve and non-CM point. Moduli,
Torelli and Brauer descent are written proof inputs; independent
implementation, formal verification, external review and novelty are
unclaimed.

The next question is whether this actual NS admits a rootless primitive U.
Arithmetic admission alone does not imply MW17, as the globally rootful
388 and622 sources demonstrate. No maximum MW rank or MW17 equation is
asserted here. This admits one of the601 remaining historical rows,
leaving600 other rows unadmitted.
