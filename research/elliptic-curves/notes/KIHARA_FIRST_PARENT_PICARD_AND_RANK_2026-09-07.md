# First Kihara parent: Picard ranks 17/18 and a full rational basis

For the new fixed Kihara parent attached to path parameter **3/2**, the
rational/geometric Néron–Severi ranks are exactly **17/18**, and the current
fibration has generic Mordell–Weil ranks exactly **12/13** over
Q(T)/Qbar(T). Its full rational MW basis has height determinant **189**;
the full rational divisor lattice has determinant **+756**.

This is a geometric lattice type different from both production X948 and
the six determinant468 parents, which have geometric Picard rank nineteen.
It is not a claim that the surface or its lattice is new to the literature.
Only this first fresh Kihara parent is covered; the other three retain
unknown exact generic ranks and lattice types.

**Every Jacobian elliptic fibration over Q on this parent has generic
rank at most fifteen.** A change of fibration on it cannot provide an MW16
or MW17 construction. This does not bound the ranks of its rational
fibres, and does not establish inferior discovery yield.

Authority: `EC-KIHARA-FIRST-PARENT-PICARD-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json). The
[rank certificate](../../artifacts/generated-results/elliptic-curves/kihara_first_parent_rank_v1.json)
and [standalone replay](../../artifacts/generated-results/elliptic-curves/kihara_first_parent_rank_replay_v1.json)
contain the exact witnesses.

A [subsequent conic base change](KIHARA_CONIC_RANK13_CONSTRUCTION_2026-09-07.md)
now supplies a primitive generic rank13 subgroup on a chi4 surface and
one certified rank13 fibre. This is a base change, so the K3 fibration
upper bound below does not apply to its generic rank.

## Symmetry closure of the twelve sections

Start with the first new parent from the
[completed four-parent pilot](KIHARA_PARENT_EXPANSION_AND_SEED_INDEX_2026-09-07.md).
Its even short Weierstrass coefficients have degrees eight and twelve,
and its discriminant has degree twenty, is squarefree and is coprime to
c4. Infinity is split I4. Thus the trivial lattice is U plus A3, of rank
five and determinant of absolute value four.

The experiment froze all twelve supplied sections and all twelve images
under the base involution sigma:T->-T before computing any heights. The
[height certificate](../../artifacts/generated-results/elliptic-curves/kihara_section_involution_v1.json)
checks every section equation, local component, doubling height,
parallelogram identity and exact relation. The 24 sections still span
rank twelve: the reflected sections add no rational direction.

The original twelve-section Gram matrix is positive definite, with

```
determinant = 6804,
sigma^2 = 1,
dim ker(sigma-1) = 6,
dim ker(sigma+1) = 6.
```

The action is verified in the function-field elliptic group, not just as
an isometry of a numerical matrix. The independent verifier recomputes
78 heights from multiples landing on identity components:

```
height(P) = (4 + 2 (4P.O)) / 16,
2 (4P.O) = max(deg denominator(x(4P)), deg numerator(x(4P))-4).
```

All reducible component groups have exponent four. Consequently the local
correction vanishes on 4P. This replay avoids the producer's blowup and
tangent-direction component implementation.

## The geometric thirteenth section and the missing Galois sign

The quotient s=T^2 has coefficient degrees four and six and discriminant
degree ten, with ten I1 fibres and split I2 at infinity. It is a rational
elliptic surface. Its geometric divisor rank is ten and its geometric MW
rank is seven. The six invariant rational sections of the K3 descend as
six independent rational sections on this quotient.

Pullback therefore gives seven geometric invariant MW directions on the
K3. The six verified anti-invariant sections are independent of them.
The K3 has geometric MW rank at least thirteen and geometric NS rank at
least eighteen. This initial proof did not construct the additional geometric section
or its field. The [subsequent exact section proof](KIHARA_QUADRATIC_SECTION_AND_PARENT_COVERAGE_2026-09-07.md)
now gives the full geometric NS lattice of determinant-756 and its exact
constant field Q(sqrt(-3)).

The quotient's known rational divisor rank is nine: U, its I2 component,
and six rational sections. Its remaining geometric divisor direction has
Galois character of order at most two. At a good prime, its eigenvalue is
therefore eps*p, eps in {+1,-1}; the quotient point count determines eps.
The same character occurs on the additional geometric section pulled back
to the K3. This determines the otherwise missing sign in its Frobenius
calculation without counting a third extension field.

## Complete finite-field witnesses

The protocol chooses the first two primes in 31..251 passing the fixed
20I1+splitI4 good-reduction gate and smooth-T=0 condition: **53 and 83**.
There is no prime extension or refill. At both primes:

| p | Quotient #Y(Fp) | K3 #X(Fp) | K3 #X(Fp²) | eps |
|---:|---:|---:|---:|---:|
| 53 | 3234 | 3576 | 7936320 | -1 |
| 83 | 7554 | 8124 | 47587176 | -1 |

Every finite fibre is counted, with the nodal correction at infinity
included. The producer uses norm-character sums in Fp² and checks smooth
fibre counts with PARI. The independent replay reconstructs the same
finite fields, recomputes all smooth elliptic cardinalities in Sage/PARI,
and counts nodal cubics by their split type. It also reselects the primes
and checks the quotient counts.

The K3's known seventeen rational divisor directions contribute 17p to
the first trace. The extra quotient direction contributes -p. After
removing these eighteen directions, the residual degree-four Frobenius
polynomials are

```
p=53: z^4 + 82 z^3 + 5724 z^2 + 230338 z + 7890481
p=83: z^4 + 94 z^3 + 1992 z^2 + 647566 z + 47458321.
```

Their middle coefficients are nonzero, so orthogonal reciprocity forces
determinant sign +1; the first two traces determine every coefficient.
The verifier checks the weight-two Weil bounds and checks that the
normalized residual polynomial has no cyclotomic factor of degree at most
four. Orders 1,2,3,4,5,6,8,10,12 exhaust the possibilities.

Good-reduction specialization injects the characteristic-zero geometric
NS group into cohomology. Divisor eigenvalues are p times roots of unity.
The only such eigenspaces here are the seventeen rational ones and the
one negative quotient direction. Thus each prime separately bounds the
geometric NS rank by eighteen and the rational NS rank by seventeen.
The lower bounds match. No Tate-conjecture converse or Artin–Tate
squareclass argument is needed.

Shioda–Tate with U+A3 now gives exactly MW12 over Q(T) and MW13 over
Qbar(T). This Frobenius calculation identifies dimensions. The subsequent section
proof linked above separately closes the full geometric lattice and its
defining quadratic character.

## Integral basis, not merely a rational span

Let P1 be the quartic origin among the twelve product-root points, and let
S_i=P_(i+1)-P1, i=1,...,12, be the original section basis. If K is the
hyperelliptic involution of that origin viewed on the Jacobian, then

```
6K = S1 + ... + S11.
```

The divisor identity was explained in the preceding Kihara intake note;
here it is verified exactly over Q(T). Replacing S1 with K has index six,
so the corrected Gram determinant is 6804/36=189 and the associated
rank-seventeen divisor determinant is 4*189=756.

Every nonzero section has height at least 4-1=3 on 20I1+I4. Hence geometric
torsion is trivial. If the corrected divisor lattice had index n in the
full rational NS lattice, integrality would require n² to divide 756.
The only possible indices are 1,2,3,6.

The corrected twelve sections specialize to explicit integer combinations
of the independently certified fourteen-point fibre seed. The stored
12-by-14 matrix is checked by exact rational group law. Independent
complete finite groups show that these twelve specialized sections inject
modulo both two and three. Any generic relation witnessing divisibility
by either prime would survive specialization and contradict that finite
injection. This excludes all nontrivial possible indices. Thus the new
basis is the full rational MW basis, and **+756 is the full rational NS
determinant**, not an unsaturated sublattice determinant.

This repairs the generic seed subgroup. It does not increase the rank of
the already searched specialization or add a discovered direction.

## Consequence for searches

For any other Jacobian elliptic fibration over Q on this same surface,
rational Galois invariants in Shioda–Tate give

```
generic MW rank <= rational NS rank - 2 = 15.
```

That is a complete upper bound across fibrations, unlike a bounded divisor
or neighbour enumeration. It rules out this parent as a generic MW16–17
source. A useful MW15 presentation remains unconstructed, and its existence
is not proved by the upper bound. Higher-rank specializations remain open;
the preceding rank14 seed has two directions beyond this generic rank12,
already built into the retained Kihara path. The pilot discovered zero
further directions. These separate counts must not be combined into a new
rank-gain claim.

No rootless-frame enumeration, different-NS foundry admission, point
search or new parameter population was launched in this proof experiment.

## Replay and costs

Copy the [bundle](../../artifacts/generated-results/elliptic-curves/kihara_first_parent_rank_bundle_v1.json)
and [standalone verifier](../cas/verify_kihara_parent_rank.sage) into an
empty directory and run:

```sh
sage -python verify_kihara_parent_rank.sage --input kihara_first_parent_rank_bundle_v1.json
```

The completed independent replay takes **27.019968436 seconds** within a
120-second, 2-GiB, one-worker cap. It imports no repository modules. The
height calculation, two-prime count, saturation build and independent
replay total **40.339151416 supervised seconds**. Protocols, all finite
counts and exact function-field witnesses are retained. Check the bindings:

```sh
python3 elliptic-curves/cas/record_kihara_parent_rank_replay.py --check
```

The geometric inputs are standard elliptic-surface, height and
Shioda–Tate results, reviewed in
[Schütt–Shioda, *Elliptic Surfaces*](https://arxiv.org/abs/0907.0298).
The independent computation establishes their hypotheses on this exact
parent; the conclusions are restricted to it.
