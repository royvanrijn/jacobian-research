# Full rational and geometric NS lattices for the six Mestre parents

The six Q-distinct parents `u=11,13,17,19,23,29` have **identical rational
NS intersection matrices** in the supplied divisor bases. Their geometric
NS extensions also have **identical matrices**, determined by a unique
index-two glue. Thus these six different surfaces share one additional
geometric lattice type, of determinant468. They do not supply six different
NS lattice types.

This closes the Gram-matrix gap in the
[Picard, saturation and determinant proof](MESTRE_DETERMINANT468_PARENTS_2026-09-07.md).
It does not construct another elliptic fibration, determine a marked
transcendental complement, or identify the quadratic field defining the
missing geometric section. No frame, neighbour or coefficient search runs.
The rational generic-rank ceiling16 remains; no new high-rank elliptic curve
is claimed.

Follow-up: [nine explicit chord presentations](MESTRE_EXPLICIT_CHORD_FIBRATIONS_2026-09-07.md)
now compile on the u11 parent. Their exact generic Q ranks are7–10, so this
bounded change-of-fibration trial gives no improvement over the original11.
The maps and complete fibre geometry replay; no point search follows.

The subsequent [complete degree-two rank bound](MESTRE_DEGREE_TWO_RANK_GATE_2026-09-07.md)
uses these matrices to prove that every Q-Jacobian fibration of old degree
two has arithmetic rank at most13. This goes beyond the chord dictionary;
higher-degree fibrations remain open.

The [degree-three component gate](MESTRE_DEGREE3_COMPONENT_GATE_2026-09-07.md)
now bounds all old-section triangles by MW13 and all degree-three
Q-Jacobian fibrations by MW15. The split I4 fibre forces an untouched
rational component; useful MW14–15 constructions must be nontriangle.

## Rational divisor basis and intersection matrix

Use the basis

```
F, O, R_minus, R_plus, R4_1, R4_2, R4_3,
C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C12.
```

The two finite I2 fibres are ordered by their rational base values.
`R_minus,R_plus` are their nonidentity component classes. The infinity I4
components are labelled by the rational tangent used in the earlier height
audit. Each `Ci` is a rational section: `C0=P0`, `2Ci=Pi-P0` for `i>0`,
where the `Pi` are the coherently labelled quartic covariant images.

The [matrix certificate](../../artifacts/generated-results/elliptic-curves/mestre_rational_ns_gram_v2.json)
contains the full18-by-18 matrix, all eleven generic section coordinates,
their intersection and component profiles, and an integral splitting of U.
The six matrices agree entry for entry.

The generator computes the intersections using

```
F^2=0, O^2=Ci^2=-2, F.O=F.Ci=1,
Ci.Cj = 2 + Ci.O + Cj.O - <Ci,Cj> - contr(Ci,Cj),  i != j.
```

The five fibre components have intersection form `-A1 -A1 -A3`.
A section meets the component specified by its exact local reduction.
For In components i,j, the local correction is
`min(i,j)-i*j/n`. Intersections with O are computed from poles of the
minimal weighted x-coordinate. Every resulting entry is integral and
the intersections of distinct sections are nonnegative.

The full matrix G is even, has determinant **−468**, and has signature
**(1,17)**. The signature is verified by splitting off U using
`F,O+F`, then checking that the remaining sixteen-dimensional form is
negative definite. This is an intersection-form check, not a rootless-frame
test. Its discriminant group is

```
A_G = Z/6 + Z/78.
```

The previous exact generic-rank and saturation certificate proves this
basis generates the full rational fixed NS lattice. The present matrix
calculation uses that result; determinant equality alone would not prove
saturation.

## The geometric extension has only one possible glue

The preceding proof supplies a primitive anti-invariant geometric divisor
D of square **−4**. It is orthogonal to the rational fixed lattice. The
full geometric lattice has index2 over `G + <-4>`; both the existence of
this glue and its index were proved using the geometrically rational
elliptic quotient and its narrow height-two anti-invariant section.

Because the fixed and anti-invariant lattices are primitive, an index-two
glue generator must have the form `(v+D)/2`, with

```
G*v = 0 mod 2,       v != 0 mod 2,       v^t*G*v = 4 mod 8.
```

The first condition makes its pairings integral and the last makes its
square even. The mod-two kernel has dimension2, so exactly four classes,
including zero, must be checked. The three nonzero representatives have
`v^2/4 = 2,0,9`; only the last gives the required odd value.
Consequently the geometric extension is unique relative to the displayed
rational marking and choice of sign of D.

Its numerator is

```
v = F + O + R_minus + R_plus
    + C1 + C4 + C5 + C6 + C7 + C8 + C9.
```

Here `v^2=36`. Replace F in the nineteen-element basis of
`G + <-4>` by `H=(v+D)/2`, keeping the other seventeen rational basis
elements and D. The resulting integral even Gram matrix has

```
signature (1,18), determinant +468,
discriminant group Z/3 + Z/156.
```

The certificate gives all nineteen rows and the rational change-of-basis
matrix implicitly through v and its replacement index. The nontrivial
element of the Galois **image** fixes the rational classes, sends
`D -> -D`, and sends `H -> H-D`. Its integral matrix preserves the Gram
form, squares to identity and has fixed rank18.

This determines the nontrivial image generator, not the quadratic character
of the absolute Galois group: the field over which D becomes rational is
still **UNKNOWN**. An explicit equation for the geometric anti-invariant
section is also **UNKNOWN**. The abstract divisor class and the existence
of its index-two extension suffice for the lattice conclusion.

## Independent replay and preserved failures

[`certify_mestre_rational_ns_gram_v2.sage`](../cas/certify_mestre_rational_ns_gram_v2.sage)
builds the matrices for precisely the six existing parents. It completes
in1.606167292 seconds under a180-second,2-GiB,one-worker protocol.

[`verify_mestre_ns_gram_v2.sage`](../cas/verify_mestre_ns_gram_v2.sage)
runs in a fresh directory with only Sage and the two JSON inputs. It
imports none of the producer's code. It checks:

- all supplied rational-function section and doubling identities;
- local component intersections using reciprocal rational functions at
  infinity and direct finite-node tests;
- **396 independent generic heights**, using quadrupling into the narrow
  subgroup, and the height form as the negative Schur complement of the
  seven-dimensional trivial lattice;
- the exact U splitting, integral discriminant groups, all four possible
  two-primary classes, the unique even extension, and its Galois-image
  generator;
- equality of all six rational and geometric matrices.

The independent run passes in **7.125597383 seconds**, within its
120-second and2-GiB cap. The
[replay record](../../artifacts/generated-results/elliptic-curves/mestre_ns_gram_independent_v2.json)
binds the copied sources and inputs to the completed transcript.

Two failed initial runs are preserved under the local evidence directories.
The generator V1 completed the six arithmetic calculations but could not
serialize a Sage integer. The independent verifier V1 computed Smith
invariants over QQ instead of ZZ and rejected its comparison. V2 explicitly
converts JSON integers and the already-integral Gram respectively. Neither
failure changes the models, section data, geometric formulas or glue roster.
No failed run is counted as a passing certificate.

To replay, place the verifier and these two inputs in one directory:

```sh
sage -python verify_mestre_ns_gram_v2.sage \
  --bundle mestre_468_replay_bundle_v1.json \
  --certificate mestre_rational_ns_gram_v2.json
```

The independent rank/Picard/saturation replay in the preceding note remains
a prerequisite for identifying the lattice as the full NS fixed lattice.
This supplemental check supplies its intersection matrix and the unique
geometric extension; it does not replace that earlier proof.

## Consequence for the parent search

The portfolio has six additional Q-distinct surfaces, but this cohort
occupies one further geometric NS lattice type. A parameter expansion
inside this component should not be advertised as more lattice diversity.
This finite result does not classify every parameter on the component.

The new matrices supply concrete divisor data for subsequent construction
work. A useful rank16 fibration is still unconstructed; the marked
transcendental complement and arithmetic-MW17 foundry gates are not closed
by these matrices. The completed point and score experiments remain frozen.
