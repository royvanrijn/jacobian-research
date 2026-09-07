# Adaptive half-lattice visibility: definitions and proved statements

The frozen selector, followed by a separately frozen budget-only continuation,
autonomously gives **17 → 19 → 22 → 23 → 26 → 28 → 28** on 302. All 540 charts
complete. Ten score-stratified determinant-1092 controls each complete their
first 50 charts and stop at 17. Thus all 1,040 charts complete without timeout.
Every stage passes separate exact mod-2, mod-3 and mod-5 rank checks. The
selector itself is rederived for every chart under an artifact-read guard;
a deliberate attempt to open the retrospective target artifact is rejected.
The original four-wave budget reached 26; the continuation changes only the
wave limit to twelve, and stops after two further waves at the first no gain.
This is eleven discovered independent directions above the generic 17. Full
autonomous recovery to 31 is **not established**.

Results: [independent replay and certificates](../../artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v1.json),
[final audit](../../artifacts/generated-results/elliptic-curves/visibility_cascade_summary_v1.json),
[complete 308-cell table](../../artifacts/generated-results/elliptic-curves/visibility_cascade_table_v1.tsv),
[earliest-support matrix](../../artifacts/generated-results/elliptic-curves/curve302_visibility_cascade_support_matrix_v1.json),
and [scientific figure](../../artifacts/generated-results/elliptic-curves/visibility_cascade_v1.svg).

The seven original residual directions have the following first cheap
representatives in the **declared retrospective finite atlas**, after removing
the confound of centres being introduced later than their nonzero support:

| Direction | M17 vetted height digits | M24 vetted height | First cheap stage | Coordinate | Generic base shell |
| --- | ---: | ---: | --- | --- | ---: |
| strict-01 | 176 | 2,108,115,427 | M27 | 373/4 | 10 |
| strict-02 | 67 | 729,881,574 | M26 | -2411/2010 | 10 |
| strict-03 | 152 | 4,294,281 | M27 | 22/37 | 8 |
| strict-04 | 147 | 1,215 | M24 | 1215/1103 | 10 |
| strict-05 | 165 | 262,960,073 | M28 | -368/219 | 10 |
| strict-06 | 4 | 2,945 | M17 | 2945/581 | 10 |
| strict-07 | 88 | 16,705,612 | M28 | -286/249 | 8 |

In particular, strict-04 did not require M25 for this witness, strict-03 did
not require M29, and strict-07 did not require M30. The exact historical order
contains scheduling delays as well as subgroup effects. All seven winning
parities at M17 and M24 were absent from the corresponding old selected
49-chart lists. Later old-policy entries are blank because that policy was
not executed at those stages.

The two matrices have 154 cells each (14 directions at 11 historical subgroup
stages). Independent exact rational group and quartic checks cover 5,418 and
5,544 chart/translation trials respectively. A second pass recomputes every
retained minimum and cheap-witness count. Once an exact containment endpoint
is found, subsequent cheap coordinates are not counted as evidence of a new
direction; the figure greys those cells. Historical M28 and autonomous M28
are different subgroups and must not be identified just by their rank.

The controls were selected by parameter ID within six fixed score/height
strata: six strong, two moderate and two lower-fixed. Quartic coefficient
medians are 113–130 bits on controls and 118 bits across the 302 cascade, so
gross coefficient size does not by itself explain the contrast. This small,
purposive comparison establishes the recorded difference in recovery, not
that these fibres have exact rank17 or that no ordinary fibre can cascade.

The remaining policy bottleneck is explicit: only 32 generic anchors feed
the extension lane, and it considers at most 64 binary extension masks. At
rank28 the full extension space over M17 has 2^11=2,048 masks per base parity.
The terminal result therefore leaves substantial unscheduled geometry. The
experiment did not retune that shortlist after inspecting missing labels.

The experiment is implemented in
[`adaptive_visibility_cascade.sage`](../cas/adaptive_visibility_cascade.sage).
Its selector is calibrated on the completed 302 history but receives only
generic section formulas, the generic degree-two orbit table, the curve
parameter, and its own discoveries. Exceptional point coordinates, residual
labels, and privileged orbit IDs are absent from its input and ranking.
The separate retrospective matrix is implemented in
[`measure_visibility_cascade.sage`](../cas/measure_visibility_cascade.sage).

The following propositions are elementary deductions, with proofs included.
They describe what subgroup enlargement guarantees and where empirical
coordinate visibility requires additional information. They do not give a
rank upper bound, a point-existence theorem for a Selmer class, or a guarantee
that the bounded selector recovers a specified rank.

## 1. Degree-two orbits and parity

Let E be an elliptic curve over a characteristic-zero field K with origin O.
The sum map identifies Pic^2(E)(K) with E(K): the class O+q has label q.
Pushing a divisor forward by translation by a sends its label to q+2a.
Consequently the translation orbits under a subgroup M of E(K) on Pic^2(E)(K)
are E(K)/2M. The subset of divisor classes with labels in M has orbit set
M/2M. If M is free of rank r, this subset has 2^r orbits.

Proof. Under Pic^0(E)(K)=E(K), a degree-two divisor D corresponds to D-2O.
Translation adds a to each of its two points, hence adds 2a to its sum.
This calculation extends from geometric point divisors to their divisor
classes and is defined over K. Orbit equivalence is exactly difference in
2M. Restrict labels to M and quotient to obtain M/2M. ∎

Thus 131,072 labels in the generic rank-17 lattice are a parity universe.
Effectivity, genus, irreducibility, rational parametrization and splitting
of a specialized fibre are additional conditions. A norm-8 or norm-10 label
of the original generic lattice describes a base orbit. An extended centre
using newly found specialized points need not be a generic section of that
norm; the experiments record the generic base label and extension separately.

## 2. Inclusion and the precise doubling statement

For torsion-free finitely generated groups M subset N, the inclusion-induced
map M/2M -> N/2N has kernel (M intersect 2N)/2M. More precisely there is an
exact sequence

\[
0\longrightarrow (N/M)[2]\longrightarrow M/2M
\longrightarrow N/2N\longrightarrow (N/M)/2(N/M)\longrightarrow0.
\]

Proof. Send a class n+M of order dividing two to 2n+2M. This is well-defined;
injectivity follows because N has no 2-torsion. Its image is precisely the
displayed kernel. The remaining exactness follows by taking classes in N/M. ∎

If N=M direct-sum ZQ_1 direct-sum ... direct-sum ZQ_k, old parity labels embed
and every chosen base parity has exactly 2^k extensions in these coordinates.
For a nonprimitive inclusion, old labels can merge. For example 2Z subset Z
maps both old parity classes to zero. The doubling claim therefore needs the
direct-summand condition (or an equivalent primitive inclusion).

For the stored historical 24-to-31 chain each later independent basis has the
earlier basis as its literal prefix. Together with exact independence, this
proves the direct-summand condition for these *displayed subgroups*. It does
not assert saturation in the whole Mordell--Weil group.

The translation quotient N/M is Z^k in this case: it has infinitely many
cosets. It is the parity extension space that has 2^k possibilities, not the
whole translation space. CVP still optimizes over infinitely many integral
translations within each parity.

## 3. Intrinsic half-lattice distance and covering

Work modulo torsion in V=E(K) tensor R, with a positive definite height pairing
on the finite-dimensional span under consideration (for a number field, the
canonical height). Write the squared norm as h. For a lattice M define

\[
d_M(P)^2=\inf_{m,q\in M}h(P+m-q/2)
        =\operatorname{dist}(P,\tfrac12M)^2.
\]

Proof. The set of m-q/2 is exactly (1/2)M: one inclusion follows from
2m-q in M, and the other by setting m=0. Lattices are discrete, so the
infimum is attained. ∎

For M subset N, d_N(P)<=d_M(P). If N=M+sum ZQ_j, then

\[
d_N(P)^2=\inf_{k\in\mathbf Z^s}
d_M(P-\tfrac12\sum_j k_jQ_j)^2.
\]

Both statements follow by inclusion or decomposition of the half-lattice.
If W=span_R(M), P=P_parallel+P_perp orthogonally, and rho(M) is M's covering
radius in W, then

\[
\|P_\perp\|^2\le d_M(P)^2
 =\|P_\perp\|^2+\operatorname{dist}(P_\parallel,\tfrac12M)^2
 \le\|P_\perp\|^2+\rho(M)^2/4.
\]

Proof. Orthogonality gives the equality; the definition of covering radius
and scaling by one-half give the upper bound. ∎

Rank growth has no positive universal improvement factor at the lattice
level. Take M=Ze_1, N=Ze_1+Z(2Te_2), P=Te_2+epsilon e_3 in Euclidean space.
Then d_M(P)^2=T^2+epsilon^2 and d_N(P)^2=epsilon^2. The ratio can approach zero.
There need not be any drop: take P=e_3 with M=Ze_1 and N=Ze_1+Ze_2. These are
lattice counterexamples to a guarantee based only on rank; they do not assert
realization of arbitrary Gram matrices by elliptic-curve heights.

## 4. Coordinate complexity and monotonicity

A chart is an exact degree-two function t_q:E->P^1 together with its retained
horizontal PGL2(Q) coordinate change A. For a rational point R use projective
height H(A^{-1}t_q(R))=max(|a|,|b|) for primitive [a:b]. Known pointed endpoints
are recorded separately from discoveries. A finite atlas consists of exact
chart/translation pairs (q,A,m); define

\[
C_\mathcal A(P)=\min_{(q,A,m)\in\mathcal A}H(A^{-1}t_q(P+m)).
\]

For nested atlases A_M subset A_N, C_A_N(P)<=C_A_M(P), and the number of
distinct retained witnesses of height at most B cannot decrease.

Proof. A minimum over a superset cannot increase, and every witness counted
before remains counted afterward. This also holds for an infinite union
using infima. ∎

This proposition applies to the matrix's persistent witnesses. It does not
imply monotonicity of the *newly computed* single CVP representative, of a
top-k scheduling list, or of a newly reduced horizontal coordinate system.
If a CVP solver or basis selection drops old data, inclusion of search spaces
must be checked again. Equal-rank saturation can change chart choices too.

There is no coordinate-height conclusion from lattice distance alone with
unrestricted horizontal charts. The same geometric point and centre can have
coordinate 0 or T under an affine coordinate translation, while their height
pairing is unchanged. A fixed map A gives finite distortion bounds (for an
integer matrix, H(Ax)<=2 max|A_ij| H(x)); uniform conclusions require uniform
control of those matrices and their inverses. This is why the matrix records
exact coordinate transports and quartic coefficient bit sizes in addition to
rounded height scores.

## 5. Bounded experiment and what it measures

The frozen pilot enumerates the complete exported generic shells 8 and 10,
ranks canonical representatives by specialized rounded height, and takes the
next 25 unsearched representatives per shell. At later stages it also reduces
the first 64 extension masks of the first 16 anchors per shell and takes 48
deep representatives. Every wave retains its exact witnesses. The budget is
height 125,000, ten seconds per chart, at most 98 charts per wave and four
waves. The first complete no-gain wave is a policy stop, not an arithmetic
absence result. A timeout is censored. In particular this finite prefix does
not exhaust large extension spaces or generic shells.

Every returned point is replayed through its exact map. Each completed cloud
is independently recertified modulo 2, 3 and 5; the mod-2 independent cloud
seeds the next wave. Complete-cloud certification can find directions that
online parity admission missed. The frozen policy uses no target rank to
choose its centres; 32 is only a declared terminal budget endpoint.

The retrospective matrix uses all fourteen original diagnostic directions
at M17, M19, M22, M24, and every M25 through M31. Its atlas is the original
32 vetted charts per direction plus all historical recovery centres as their
seed subgroups become available. It retains exact earlier witnesses and tries
a new target-relative CVP translation in each chart. It reports finite-atlas
minima, cheap-witness counts, coordinate matrices and coefficient sizes. This
design measures certified visibility drops, but does not establish global
minimum coordinate heights across M/2M.

## 6. Conditional eventual completeness and the finite-budget boundary

Fix any nonzero rational centre q and its exact degree-two chart. Every
rational point other than the pointed endpoints has a rational projective
chart coordinate of finite height and a rational ordinate. Therefore an
exact point-search procedure that exhausts all projective heights on that
single chart, with complete ordinate and endpoint handling, eventually
encounters every rational point through the degree-two map.

Proof. Apply the rational inverse-chart function to the chosen rational
point. Its reduced numerator and denominator are finite integers; some
finite height box contains them. Exact square testing and the map recover
the point. The finitely described pointed endpoints are handled separately. ∎

The same conclusion holds for a fair adaptive atlas that keeps this chart
and gives it unbounded height. It supplies no finite stopping criterion or
useful time bound. Extra parity charts change bounded computational visibility;
they are not needed for this weak eventual-completeness statement. The present
policy is explicitly finite and does not meet the unbounded-height hypothesis.
