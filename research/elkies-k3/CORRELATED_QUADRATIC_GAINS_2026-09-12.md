# Two independent gains on one quadratic cover

The construction target remains **OPEN**. The
[Mestre construction](R17_MESTRE_CORRELATED_SECTIONS_AND_GENUS_GATE_2026-09-13.md)
now supplies two independent new sections on one quadratic cover of published
R17, with height matrix `diag(24,24)` and subgroup rank at least19. Its covering
base has genus21. Every rational auxiliary function in that identity still
forces genus at least9, so it cannot meet the infinite-rational-base requirement.

The first prospective calculation excludes all **300 pairs of 25 complete
genus-one pencils** on published R17. The exclusion covers all rational pencil
parameters, including infinity, for this fixed bank and supplies no new section.

The subsequent [one-node construction gate](R17_ONE_NODE_CORRELATED_COVER_GATE_2026-09-12.md)
also excludes all 1,675 pairs of 67 height-six regular chord nets with these
25 smooth pencils, including every rational node position and all projective
parameter boundaries. Its retained construction attempts supply no gain.

A different [shared-ordinate tangent construction](R17_SHARED_ORDINATE_TANGENT_GENUS_GATE_2026-09-13.md)
from every pair of the 17 generic basis sections also fails the base condition:
each of its 136 covers has genus at least86. This tests one tangent per basis
pair; arbitrary generic words and further iterations remain outside the gate.

On the alternate-Q80 parent, the
[complete smooth genus-one comparison](Q80_ALL_SMOOTH_GENUS_ONE_BISECTIONS_2026-09-13.md)
now covers all 63,966 pencils modulo inherited section translation, including
the 49 norm-twelve pencils with a moving intersection with zero. All
2,045,792,595 pairs are excluded over rational parameters and every individual
pencil is injective. A solution on this parent with a genus-one quadratic
base must therefore use a singular bisection image of arithmetic genus at
least two in at least one independent direction. This does not exclude that
possibility or a construction on another parent.

## Required positive endpoint

Construct a nonconstant quadratic extension of the original parameter field
of a high-rank parent, with two new sections independent modulo the inherited
subgroup. Construct the cover and sections without fitting the desired
exceptional points. The smooth projective covering base must have an actual
source of infinitely many rational points: a rational parametrization or a
genus-one model with a certified rational nontorsion point. Verify the section
identities, the rank increment and the specialization hypotheses separately.

For an MW17 parent this targets a subgroup of rank at least19 over one
quadratic cover. Neither exact total rank nor an explanation of all fourteen
Curve302 directions is required. The motivating
[ancestry calculation](../elliptic-curves/notes/PRINCIPAL28_EXCEPTIONAL_ANCESTRY_2026-09-12.md)
only separates its selected covers; it does not exclude this endpoint.
An alternative source from the [different-NS foundry](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md)
must still pass the rational marking gate and demonstrate an arithmetic
advantage beyond another rootless frame or MW17 equation.

The retained [biquadratic constructions](BISECTION_PAIR_COVER_GEOMETRY_2026-08-31.md)
already give two directions over a degree-four cover of the original line.
Their diagonal quotient is a useful degree-two map between genus-one curves,
but it keeps degree four over that line. Reinterpreting that tower does not
close the present stronger target. The original published construction has
the same degree-four character structure; see
[Elkies, abstract](https://arxiv.org/abs/2608.25406).

## Why test complete smooth genus-one pencils

The completed smooth rational-bisection maps on published R17 and alternate
Q80 are injective. The
[rational-normalization boundary](R17_NORM12_RATIONAL_NORMALIZATION_BOUNDARY_2026-09-04.md)
also excludes the recorded singular arithmetic-genus-one route to a rational
normalization. Those calculations do not exclude common **smooth genus-one**
quadratic covers. This is the reason for the present new test.

The input uses only the published equation, generic basis and its certified
height Gram. Among words `P_i + epsilon*P_j`, with `i<j` and
`epsilon in {1,-1}`, retain those of height8. Deduplicate modulo2, preferring
`epsilon=1`, and sort their 17-coordinate words lexicographically. There are
exactly25; the predeclared cap was64. No exceptional point, specialization
parameter or quotient label enters the worker packet. This is field-level
allowlisting, not process isolation or a claim that the parent was selected
without knowledge of its published significance.

For each trace `T=(Nx/h^2,Ny/h^3)`, `deg h=2`, put

```
M = M0 + lambda*h^2,
M0*Nx + Ny = 0 mod h^2,          deg M0 < 4,
q_lambda = (M^4-6*M^2*Nx-8*M*Ny-3*Nx^2-4*A*h^4)/h^6.
```

The [regular chord theorem](RANK_MUTATION_AND_LIFT_THEOREMS.md#proposition-f11-the-height-eight-genus-one-bisection-pencil)
identifies `s^2=q_lambda(t)` and its residual section. The exact coefficients
give a 5-by-5 matrix `B_T` with

```
coeff_t(q_lambda) = B_T*(1,lambda,lambda^2,lambda^3,lambda^4)^t.
```

Each matrix is invertible, and its member at `lambda=0` is a squarefree quartic,
so the pencil has a nonempty smooth genus-one open set.
Homogenization therefore embeds the whole pencil
as a rational normal quartic in the projective space of binary quartics in
the **fixed t-coordinate**. A smooth genus-one double cover has four simple
branch points, so equal quadratic extensions require proportional binary
quartics. This condition is weaker than equality of rational squareclasses,
because it forgets the constant twist; excluding it also excludes equality
of extensions. It must not be used to accept a cover match.

Within one pencil the embedding is injective. Opposite residual points have
the inherited trace and add at most one rational direction. Thus the frozen
bank needs a match between two distinct pencils before independence testing.

## Exact exclusion and independent replay

At an odd prime where a rational branch matrix is integral with invertible
reduction, it defines an embedding of `P1` over `Z_p`. Every rational pencil
parameter extends to a projective `Z_p` point. A rational equality of branch
images therefore reduces to equality of images of two points in `P1(F_p)`.
Disjoint finite images give an exact exclusion, even if a branch quartic
becomes singular at that prime. Nonintegral or singular matrices are deferred.

The certificate assigns every pair to one good-prime witness:

| Prime | Newly excluded pairs | Pairs left |
|---:|---:|---:|
|101|294|6|
|103|4|2|
|107|2|0|

The constructor retains every rational coefficient and pair witness. Its
independent checker uses only Python fractions and integer arithmetic. It
reconstructs the17 generic polynomial points, all25 trace sums by direct
chord formulas, the RR congruences and all125 polynomial branch coefficient
identities and checks the25 smooth members at `lambda=0`.
It then recomputes all good-reduction rosters and all projective
finite images, with an explicit infinity point, and checks exact pair coverage.
It imports neither Sage nor the producer. The original generic-rank17 theorem
is an inherited dependency, not independently reproved here.

The seven negative/control tests reject altered trace ordinates, altered branch
coefficients, missing exclusions and duplicated pencils, defer singular matrix
reductions, preserve constant-scaling matches and include projective infinity.

Artifacts: [frozen input](../artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/input.json),
[25 exact pencils](../artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/pencils.json),
[pair certificate](../artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/result.json),
[independent replay with smoothness checks](../artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/replay-with-smoothness.json).

```
python3 research/elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py
python3 -m unittest discover -s research/tests -p 'test_r17_correlated_pencil_replay.py' -q
```

Discovery uses `construct_r17_correlated_genus_one_pencils.sage freeze` and then
`run`, with `--output` naming a new directory. Set `OPENBLAS_NUM_THREADS=1`
for the run. The frozen limits are120 CPU seconds and4GiB virtual address
space. An initial startup attempt failed before constructing a pencil because
the default numerical-library thread pool reserved about5.7GiB virtual memory
while using about214MiB resident memory. The corrected invocation uses the
same input and arithmetic, with one BLAS thread; no mathematical candidate or
budget was replaced. The failure is retained in the
[execution receipt](../artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/execution.json).
Construction plus exclusion took0.174 component CPU seconds; independent
replay took0.123 wall seconds. These exclude startup and authoring and are not
an end-to-end performance comparison.

## Remaining construction gate

The result concerns these25 entire pencils. It does not exhaust all norm-eight
trace classes, singular bisections of higher arithmetic genus, arbitrary
quadratic twists, the determinant1092 parent, changes of exceptional basis,
or other arithmetic sources. No rank bound for a specialized fibre follows.

The proposed height-six one-node route has now been derived and tested
against these 25 pencils; its full coefficient correspondence is closed for
the frozen 67 nets. The [new note](R17_ONE_NODE_CORRELATED_COVER_GATE_2026-09-12.md)
retains the constructive halving formulas, the height-16/20 independence
criterion, all exact exclusions and their boundaries. The subsequent Mestre
identity supplies the required two directions, but its irreducible degree8
and degree12 coefficient divisors remain in the branch locus for every rational
auxiliary function. The [valuation proof](R17_MESTRE_CORRELATED_SECTIONS_AND_GENUS_GATE_2026-09-13.md)
closes that identity for the required base condition without a degree or height
bound. A further attempt needs a different carrier identity or arithmetic
parent, with a proved low-genus branch condition before a parameter campaign.
No enlargement or heavy computation is scheduled by this note.

The positive endpoint remains unverified: no shared quadratic cover with both
two independent new sections and infinitely many rational base points has
been constructed here.
