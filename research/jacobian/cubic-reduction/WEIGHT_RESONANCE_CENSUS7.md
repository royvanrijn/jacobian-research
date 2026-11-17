# Complete degree-seven weight census and fourteen remaining improvement types (KCUB36)

Subsequent [KCUB40](EQUAL_NEGATIVE_FULL_DEGREE7.md) closes full111 through
degree seven by a changed curve/Euler-Wronskian proof. The current
degree<=7 diagonal-weight raw<=12 frontier is empty:only112 can collide,
with exactly13 terms. All larger frontier and full111-unresolved descriptions
below are historical stages;original timeout records are unchanged.

Subsequent [KCUB39](ZERO_WEIGHT_CLASSIFICATION.md) closes every one-zero-weight
three-variable class in all degrees. The current degree<=7 diagonal-weight
raw<=12 frontier contains only full111. All fourteen/twelve/ten-type
frontiers and unresolved zero-weight descriptions below are historical
stages, not current next steps; their solver records remain unchanged.

Subsequent [KCUB37](WEIGHTED_124_125_EXCLUSION.md) closes both twelve-slot
weights(1,-2,-4),(1,-2,-5) by a changed boundary proof. That step gave twelve types; subsequent KCUB38 closes312/412, leaving
ten current diagonal-weight improvement types. The fourteen-type stage and all four
original native timeout records below are preserved as historical evidence.

Let K have characteristic zero. Consider all polynomial Keller maps in three
variables of ordinary degree at most seven, equivariant for a nontrivial
diagonal G_m action in the displayed source and target coordinates. Up to
permutation and simultaneous reversal of weights, every algebraic-closure
collision lies either in the known weight(1,-1,-2) class or in these fourteen
unresolved weight types:

| types | weights | number |
|---|---|---:|
| mixed, distinct nonzero weights | (1,-2,-4),(1,-2,-5),(3,-1,-2),(4,-1,-2) | 4 |
| repeated negative weight | (1,-1,-1) | 1 |
| one zero weight | (0,p,-q), gcd(p,q)=1, 1<=p<=q, p+q<=7 | 9 |

The nine zero-weight pairs are (1,1),(1,2),(1,3),(1,4),(1,5),(1,6),
(2,3),(2,5),(3,4). KCUB15 proves that every collision in the separate
(1,-1,-2) class has exactly thirteen nonlinear occurrences. Therefore:

> Any degree<=7 diagonally equivariant three-variable Keller map with a
> collision and at most twelve nonlinear occurrences must have one of the
> fourteen weight types above.

This is a necessary finite frontier, not existence at any retained weight
or a global support lower bound. The two full mixed models(1,-2,-4) and
(1,-2,-5) have twelve slots each, so they are direct raw13 improvement
targets. The other two mixed models have fifteen and fourteen slots.
No smaller positive witness is constructed here.

A new universal lemma also settles every rank-one active monomial-relation
span, in every dimension and degree: after identity-linear normalization,
such a Keller map is either the identity or a single elementary monomial
shear. Its proof is given below, independently of the degree-seven census.

## Normalization and the raw-count boundary

Translate F(0) in the target fixed subspace and identify the source and
target weight multisets through the invertible differential at zero.
Normalize the linear part by its inverse on the output. This intertwines
the same weights, preserves degree and collisions, and gives F_i=x_i+H_i
with no constant or linear terms in H and det JF=1. No source translation
or nonlinear linearization is needed.

Within repeated-weight blocks this output normalization can change the
number of nonlinear occurrences. Every repeated-weight class discarded
here has a full no-collision theorem, independent of sparsity. The sole
count-sensitive import is KCUB15 for distinct weights(1,-1,-2), where
normalization is diagonal and preserves occurrence counts. The unresolved
repeated-weight111 class is retained in full. Thus the fourteen-weight
necessary condition is valid for the original raw count too.

This count statement requires the action to be diagonal in the very
coordinates where occurrences are counted. An arbitrary source linear
change used to diagonalize a nondiagonal action can change occurrence
counts; in particular a different source frame of the known112 map is
not excluded from having fewer terms by this result. General source-frame
sparsification, nonlinear conjugacies and maps without the displayed
diagonal symmetry remain outside the raw-count frontier.

## All-degree rank-one active-relation lemma

For each actual nonzero coefficient of x^alpha in H_i put v=alpha-e_i.
If the span is zero the map is identity. If it has rank one, choose the
primitive integer generator v of positive coordinate sum. Every actual
relation is k*v for a positive integer k: integrality follows from
primitivity, and positivity from |alpha|-1>0.

If v has a negative coordinate, any native relation has at most one
negative entry, namely -1 in its output position. Nonnegativity of
alpha=e_i+k*v consequently forces exactly one negative entry v_j=-1,
output i=j and k=1. All other coordinates are unchanged, and the added
monomial is independent of x_j. This is a single elementary shear.

If v has no negative entries, put t=x^v. Then F_i=x_i*f_i(t), f_i(0)=1.
Conjugating its Jacobian on the coordinate torus by diag(x_i) gives a
diagonal matrix plus the rank-one update t*(f_i')*(v_j). Multilinearity
of the determinant gives the polynomial identity

    det JF=product_i f_i+t*sum_i v_i*f_i'*product_(j!=i) f_j.

If d_i=deg f_i, its leading coefficient has scalar
1+sum_i v_i*d_i times product_i lc(f_i). This integer is positive.
Thus det JF=1 forces sum_i d_i=0, and every f_i=1. The identity on the
torus extends as a polynomial identity; no coordinate-zero boundary is
lost. This proves the lemma in all dimensions/degrees over every
characteristic-zero field. It is a field argument, not a statement about
nilpotent coefficient rings.

The checker verifies the three-dimensional formal rank-one determinant
identity, all 343 degree triples through degree six that could occur in
the present census, and all 99 native negative-relation generic shear
inverses. The all-dimension/all-degree assertion rests on the written
rank-one determinant and positive-leading-degree proof, not those tests.

## Complete finite weight-ray coverage

Every native nonlinear monomial of degree two through seven gives
(alpha-e_i).w=0. There are 348 labelled native atoms, 182 distinct relation
vectors and 154 primitive relation planes. These are generated by degree
compositions, not an arbitrary bound on entries of w.

If the actual relations have rank two, two independent ones determine
ker(w) uniquely: their perpendicular line gives a primitive weight ray.
If they have rank at most one, the preceding lemma already gives an
automorphism. Rank three is impossible for nonzero w. Hence it is complete
to enumerate kernels of every pair of native relations and then take the
full allowed source basis at each ray, up to coordinate permutation and
overall sign. There are exactly 896 such weight rays.

The primary census uses integer cross products. The independent checker
uses rational row reduction for every relation pair, compares the entire
ray set, and reconstructs every full allowed basis separately. It also
reconstructs all generic rank-one-plane supports. No bounded weight scan
or hand-selected primitive ratio is a coverage premise.

The [archived degree-four note](../../foundations/verified/DEGREE_FOUR_EQUIVARIANT_CLASSIFICATION.md)
and its script supplied the arrangement approach. They were read during
preflight. No old radical computation or unregistered degree-four theorem
is imported into this proof; the new census and rank-one argument are
reconstructed here.

## First structural and imported gates

The exact partition of all 896 full rank-two bases is:

| gate | rays |
|---|---:|
| full cross-dependency graph acyclic | 193 |
| degree>3 output coefficient matrix has structural rank below3 | 199 |
| prior full weighted exclusions | 40 |
| known112 sharp13 class | 1 |
| initially unresolved | 463 |

For the acyclic gate, order coordinates topologically in the graph with
edge j->i whenever a nonlinear term of H_i uses x_j, j!=i. The determinant
is the product of diagonal factors 1+partial_i H_i. All are polynomial
units, hence constant, and their values at zero are one. In characteristic
zero H_i is independent of x_i. The map is triangular and invertible.
Self-dependence is deliberately not a graph edge; it is removed by this
unit argument.

For the structural-rank gate, form the coefficient matrix of all source
monomials of degree>3, with three output rows. Its allowed nonzero pattern
has maximum matching below three. Every 3-by-3 minor is consequently zero,
so every parameter value has a nontrivial output combination of degree
at most three. [GLDC10](../quartic-three-audit/DEGREE7_COMPONENT_EXCLUSION.md)
excludes collisions in this degree-seven domain. A possible matching of
three is not a Keller completion and is only a necessary frontend test.

The prior weighted gates, applied only within their complete scopes, are
[KCUB32](EQUAL_NEGATIVE_ALL_RATIOS.md) (10 rays),
[KCUB33](OPPOSITE_WEIGHT_CLASSIFICATION.md) (14),
[KCUB34](EXCEPTIONAL_UNIT_WEIGHT_CLASSIFICATION.md) (4),
[KCUB35](MIXED_WEIGHT2_CLASSIFICATION.md) (9), and one each for
[KCUB16](WEIGHTED_123_DEGREE7_EXCLUSION.md),
[KCUB17](WEIGHTED_213_DEGREE7_EXCLUSION.md),
[KCUB18](WEIGHTED_113_DEGREE7_EXCLUSION.md).
These total forty full exclusions after the earlier gates.
[KCUB15](WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md) contributes the separate
one-ray sharp13 class, excluded only for the fewer-than13 objective.
The broader111 degree-seven sector is not closed by the imported strict
KCUB23 result and remains on the frontier.

## 447 new complete native coefficient closures

Of the 463 initially unresolved bases, 451 have distinct nonzero mixed
weights other than111. Orient/permutate each as (p,-q,-r), q<=r.
For each entire native basis, form every coefficient equation of det JF-1,
over Q, retaining all native coefficients and every possible repair.
For each nonlinear coefficient c_j in the positive output, an exact identity
has the form

    c_j^k=sum_i M_i*I_i.

The bounded primary run obtained complete identities for all positive-output
coefficients in 447 full models. There are 594 targets, 153,878 rational
multiplier monomials, maximum multiplier degree15 and maximum target power16.
Every identity is directly multiplied by the independent checker.
Consequently all those coefficients vanish at every characteristic-zero
Keller field point. The first output is linear, so GLDC10 excludes collisions.
No generic coefficient, nonzero-stratum assumption, finite coefficient grid
or primary ideal-membership answer is trusted.

The producer forms source determinants using the existing sparse polynomial
helper and uses Singular for bounded discovery/lifts. The independent
checker reconstructs every derivative atom from the monomial exponents and
expands all six determinant permutations using integer exponent bookkeeping.
It then verifies the generator lists and multiplies every rational identity.
Thus the accepted units are checked against a separate determinant
construction, not just copied source equations or a primary status label.

Four full queries timed out at five seconds each:

| weights | full nonlinear slots | retained model |
|---|---:|---|
| (1,-2,-4) | 12 | [input](weight-resonance7/mixed-cases/p1q2r4.json) |
| (1,-2,-5) | 12 | [input](weight-resonance7/mixed-cases/p1q2r5.json) |
| (3,-1,-2) | 15 | [input](weight-resonance7/mixed-cases/p3q1r2.json) |
| (4,-1,-2) | 14 | [input](weight-resonance7/mixed-cases/p4q1r2.json) |

Their scripts, logs and terminal timeout records remain unchanged. A timeout
is no exclusion. There is no live process to await. The exact identities
exclude only the other447 models; the four untouched models are part of the
fourteen-weight frontier. Before another query on them, change the method
by deriving structural boundary/coefficient restrictions rather than merely
raising the unchanged native Groebner time cap.

## Two elementary zero-weight reductions

The remaining twelve bases consist of111 and eleven zero-weight types.
Two of the latter close without a coefficient search.

For weights(0,0,1), equivariance gives
F=(A(x,y),B(x,y),z*C(x,y)) and det JF=C*Jac(A,B)=1.
Polynomial units force C=1, so the third output is linear and GLDC10 applies.
For weights(0,1,1), the map is (A(x),M(x)*(y,z)), with a 2-by-2 polynomial
matrix. The determinant factors as A'*det M=1. Normalization gives A=x
and det M=1. Its inverse is polynomial by the adjugate formula. Both generic
block determinant identities are instances of the independently checked
block-triangular determinant formula. Other same-sign zero-weight bases
already passed through the acyclic gate. Exactly nine opposite-sign
zero-weight types remain, as listed in the opening statement.

Thus 463-447-2=14 unresolved weight types. Retaining the known112 collision
class gives fifteen possible collision types in total; removing it only for
the raw<=12 objective yields the asserted necessary fourteen-type frontier.

## Reproduction, exact boundaries and remaining goal

Run from this directory:

    python verify_weight_resonance7.py

The [receipt](weight-resonance7/validation.json) pins the complete census,
all451 native models, all447 certificate files, both terminal checkpoints,
the source producers, arithmetic helper, imported theorem scopes and this
note. The [artifact record](weight-resonance7/README.md) distinguishes the
four timeouts and the first replay's corrected variable-shadowing error.
The checker completes the full rational replay in a few seconds locally.

The universal rank-one proof, normalization, arrangement coverage,
triangular/unit arguments and rank-to-output-combination implication are
written mathematics, not a formal proof. Imported GLDC10 and prior weighted
results retain their own written/external theorem boundaries; they are not
reproved by this replay. No formal verification, external independent full
review or novelty priority is claimed.

The fourteen types are unresolved, not constructed. Arbitrary maps without
a displayed linear torus symmetry are outside this census. Whole unrestricted
strict-profile counts stay62 excluded/12 remaining, and all positive bounds
are unchanged: raw3D13, normalized3D[10,48], identity-linear-part cubic12/36,
free-linear-part cubic12/31+17 and cubic-linear89. Global degree, support and
special-form dimension minima remain open.
