# Six inequivalent MW15 triangle fibrations exclude 302

Authority: `EC-K3-CURVE302-SIX-MW15-TRIANGLES`.

Six explicit rational triangle pencils on the determinant-948 K3 have
exact generic arithmetic Mordell–Weil rank 15, trivial torsion and full
abstract MW determinant 316. Their branch divisors prove them pairwise
inequivalent, even over the algebraic closure. **None specializes to 302
at a rational parameter.** These are six fibrations on one surface.

This extends the [completed MW14 dictionary](CURVE302_SHORTWORD_TRIANGLES_2026-09-07.md)
beyond its two-term section words. It is a finite construction and inverse
result, not a classification of triangles or an exclusion of all parents.
The requested family, full coordinate basis and specialization to 302
remain unresolved for generic ranks 14–16.

## Construction and exact rank

Use the certified seventeen-section basis of the
[11952 presentation](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md),
with height matrix G. In coordinates `(O,F,v1,...,v17)`, the full NS matrix is

`N = diag([[-2,1],[1,0]], -G)`.

The class of the section belonging to a word w is
`S(w) = (1, w*G*w/2, w)`. Fix `P=P6` and use the six explicit Q words in the
[frozen protocol](../artifacts/generated-results/elkies-k3-curve302-six-mw15-triangles-protocol-v1.json).
For example, its row 0 is

`Q = P6 - P9 + P12 - P16 + P17`.

Every pair satisfies `height(P)=height(Q)=6` and `<P,Q>=3`.
Consequently `O,P,Q` intersect pairwise once. Their sum D is an effective
nef isotropic divisor. In rows 0–4, the old section `P13` has `D.P13=1`;
in row 5 the section `P6-P13` has intersection one. Thus D is primitive
and defines a Jacobian genus-one fibration over Q with a rational zero.

For each row the exact orthogonal frame has determinant 948. PARI and
independent rational LDL enumeration give precisely six roots, spanning
a primitive A2 root lattice. The root quotient has no torsion. Shioda–Tate
therefore gives `19-2-2=15`, and the full abstract MW determinant is
`948/3=316`. A Smith decomposition supplies the entire abstract MW Gram
matrix, recorded in the protocol. This is not a list of fifteen explicit
Weierstrass-coordinate sections.

The rational pencil compiler verifies finite distinct intersection
parameters for P, Q and P-Q, cancels the vertical poles, and checks that
the new base function restricts to a degree-one map on the chosen zero.
All six rational pencil functions, implicit cubic equations, zero maps,
frame matrices and root witnesses are stored before the target tests.

The preceding local height-shell pilot motivated these six words. The
theorem checks the fixed words directly and does not depend on claiming
that pilot to classify all possible triangles or all higher-rank frames.

## Inverse tests and inequivalence

The screen uses exact generic cubic conversion at the frozen primes
1013, 1021 and 1009. It preserves the preliminary result of five exclusions
and one UNKNOWN in a [separate certificate](../artifacts/generated-results/elkies-k3-curve302-six-mw15-triangles-v1.json).

| Row | Final exclusion witness |
|---:|---|
| 0 | Degree-preserving modular j-map at 1021 |
| 1 | Degree-preserving modular j-map at 1021 |
| 2 | Degree-preserving modular j-map at 1021 |
| 3 | Degree-preserving modular j-map at 1021 |
| 4 | Exact rational inverse polynomial, no projective root modulo 97 |
| 5 | Degree-preserving modular j-map at 1013 |

Row 4 survives the three initial modular screens. Its exact generic
conversion over Q(s) supplies a plane cubic, both coordinate directions,
explicit Jacobian Weierstrass coefficients and a j-map of degrees 24/21.
The companion-matrix identities and discriminant change are verified by
the [generic conversion method](CURVE302_SHORTWORD_TRIANGLES_2026-09-07.md#exact-generic-conversion-closes-row17).
The primitive integral comparison with j302 has degree 24 and no root
in the projective line over F97. This excludes every rational parameter;
no height bound or heuristic search is involved in that conclusion.

At 1013 all six j-maps retain degree 24. Their normalized binary branch
discriminants have degree 44 and are pairwise different. Because the
branch divisor on the fixed target j-line is invariant under a change of
base coordinate, this proves inequivalence of the six fibrations. It
does not assert that all MW15 fibrations have now been enumerated.

The [final certificate](../artifacts/generated-results/elkies-k3-curve302-six-mw15-triangles-complete-v1.json)
retains the exact rational coordinate-change witness, inverse polynomial,
branch signatures and eighteen independently normalized fibre values.
Those three values per pencil check agreement with the separate fibre
normalization method; they alone are not an interpolation proof of a j-map.

## Replay

```sh
sage -python elkies-k3/scripts/verify_curve302_six_mw15_triangles.sage
```

The replay rebuilds all six frames and rational pencils, recomputes the
modular j-maps, checks the supplied rational birational conversion and
exact inverse obstruction, and reconstructs the six branch signatures and
eighteen independent fibre checks. Each invocation is limited to 120
seconds with one worker; the software is SageMath 10.9. No timeout or
failed conversion is treated as an exclusion.

The construction-recovery objective remains open. MW14, MW15 and MW16
are all eligible; rank 17 is not a prerequisite.
