# Rank-jump mechanisms: solubility, dimension, and search gates

This is the canonical theorem layer distilled from the September 5–6 searches
and cover experiments. Its foundations are the
[specialization quotient theorems](../../elkies-k3/SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md)
and [residual Selmer theorems](../../elkies-k3/RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md).
Those remain authoritative for saturation, heights and Cassels–Tate radicals.
Here the question is: **what must a construction prove to explain or produce
several independent directions on a chosen fibre?**

The results below are deductions from standard algebra and cover geometry,
with exact finite replays of the retained examples. They are not new general
existence theorems for large jumps. Full ranks, the global squareclass images
of the higher-genus carriers, and a predictor for the largest jumps remain
UNKNOWN. Computational regression checks do not replace the written proofs.

## J1. The exact contribution of a specialized mechanism

Let E be an elliptic curve over a characteristic-zero function field K(B).
Let C→B be a finite cover, c a K-rational point above a smooth fibre b, and
assume the specified sections specialize at c. Tensor all groups in this
paragraph with Q. Let V be a finite-dimensional marked section space over
K(C), containing the pullback U of a marked subgroup over K(B). Set

\[
G=\operatorname{sp}_c(U),\quad k=\dim(V/U),\quad
\bar s:V/U\longrightarrow (E_b(K)\otimes\mathbf Q)/G.
\]

**Theorem J1 (mechanism dimension and specialization loss).** If
κ=dim ker(\bar s) and q=dim im(\bar s), then

\[
\boxed{q=k-\kappa,\qquad
\dim\operatorname{sp}_c(V)=\dim G+k-\kappa.} \tag{J1}
\]

If G lies in a certified independent witness space W and J=dim W−dim G,

\[
\dim W/(W\cap\operatorname{sp}_c(V))\ge\max(0,J-q).
\]

Equality holds when sp_c(V)⊂W. If only k≤K is proved, at least
max(0,J−K) witness directions are outside this mechanism.

**Proof.** The quotient map is well-defined because sp_c(U)=G.
Rank–nullity gives J1. The intersection of W/G with im(\bar s) has dimension
at most q; subtract from J. The capacity version uses q≤k≤K. ∎

No injectivity on U is needed for J1. Calling a quotient gain a jump from
the *full original generic rank* additionally requires U to span that group
and specialization to preserve its rank. Also κ is a kernel **modulo G**;
it is not automatically the kernel on all of V. In the retained quartets
U does specialize injectively, k=4 and q=3. Thus κ=1, the marked generic
rank is 21, and its specialized image has rank 20. Witness ranks 24 and 25
leave respectively four and five directions unexplained.

**Search use.** Record U, the base change, k, q or its bounds, and the
independence certificate. A mechanism with certified capacity below the
requested contribution can be retired *for that contribution*. It does not
exclude the fibre: other sections and mechanisms may contribute.

## J2. Relations bound an entire block

Let P₁,…,Pₙ∈E_b(K), G a fixed generic rational span, and R a matrix over Q
whose rows are verified relations ∑rᵢPᵢ∈G after clearing denominators.

**Theorem J2 (relation capacity).**

\[
\boxed{\dim\langle[P_1],\ldots,[P_n]\rangle
       \le n-\operatorname{rank}_{\mathbf Q}R.} \tag{J2}
\]

For relations aPᵢ+bPⱼ∈G with a,b∈{±1}, give edge i—j sign −ab.
Let b(Γ) count connected components in which every signed cycle has product
+1; isolated vertices count. Then n−rank_Q R=b(Γ).

**Proof.** R lies in the kernel of Qⁿ→(E_b(K)⊗Q)/G, proving the bound.
In the nullspace of the pair-relation matrix each edge says xⱼ=−ab xᵢ.
A spanning tree determines every coordinate from the root coordinate.
A consistent cycle leaves that coordinate free; an inconsistent cycle
forces x=−x, hence x=0 over Q. Each balanced component contributes one
to the nullity and each unbalanced component contributes zero. ∎

This is an upper bound until independent points attain it. Raw edge count
is wrong when rows are dependent. F₂ rank is wrong for an inconsistent
signed cycle: characteristic two loses the sign obstruction. Torsion and
saturation are discarded only for rational rank; residual Kummer dimension
requires the correction in the foundational theorem R1.

The 165-address panel has 18 edges in total. At −70/61 six edges have rank
five. At −4112/1937 the forced pair supplies one quotient line despite
two independent characters before specialization. The +7 quartet's triple
relation shows why pair graphs must be extended to general relation matrices.

**Search use.** Maintain the rank of verified relation rows. If
n−rank R<d, the current block cannot supply d directions. Select a basis
of the certified image for subsequent visibility work; opposite branches
and generic translations of a known quotient line are not new dimensions.

## J3. The carrier required by marked quadratic characters

Let f₁,…,fₘ∈K[t] be squarefree quadratics with disjoint geometric root sets,
m≥2. Let C be the smooth projective normalization of uᵢ²=fᵢ(t), and H the
smooth product curve y²=∏fᵢ(t). Suppose the marked directions have nonzero
independent singleton character components in E(K(C))⊗Q. For example,
a native point Pᵢ with conjugate wᵢ−Pᵢ gives anti-invariant vector
Pᵢ−wᵢ/2, provided it is non-torsion.

**Theorem J3 (marked carrier and product quotient).**

\[
[K(C):K(t)]=2^m,\quad
g(C)=1+2^{m-1}(m-2),\quad g(H)=m-1.
\]

The map C→H has degree 2^(m−1) and is unramified. Its deck group is the
even-sign subgroup, the unique largest native subgroup acting freely.
For m≥3, H has minimal genus among unramified native deck quotients.
For m=2 every unramified quotient has genus one.

Making the whole specified singleton block rational in an algebraic
extension of K(t) requires that extension to contain K(C). None of these
singleton vectors descends to K(H).

**Proof.** A root unique to fᵢ has odd valuation in fᵢ and even valuation
in all other factors, proving independence of the squareclasses even
geometrically. The deck group is (Z/2)^m. There are 2m branch points,
each with inertia generated by one coordinate sign change, and no geometric
ramification at infinity. Riemann–Hurwitz gives
2g(C)−2=−2·2^m+2m·2^(m−1). The product double cover has 2m simple branch
points, giving g(H)=m−1. The even-sign subgroup avoids every inertia
generator and acts freely. A proper subgroup has order at most 2^(m−1);
an index-two kernel avoids every coordinate vector precisely when all
coefficients of its defining linear form are one. This proves uniqueness
and the genus assertion, with the genus-one exception.
The intersection of the kernels of all singleton characters is trivial,
so Galois invariance of their independent nonzero vectors forces the whole
multiquadratic field to be fixed. Each singleton is nontrivial on the
even-sign subgroup when m≥2, proving the last assertion. ∎

The curve calculation uses standard tame
[Riemann–Hurwitz](https://stacks.math.columbia.edu/tag/0C1B).
Minimality is for the **marked block** and native quotients, not unrelated
auxiliary maps or alternative sections. A genus-one carrier can have no
rational point. For m≥3 no genus-zero or genus-one curve can dominate C,
again by Riemann–Hurwitz. This does not exclude isolated rational components
of an intersection on the original surface.

**Search use.** Four independent native characters demand a genus-17
carrier; its genus-three product curve leaves a degree-eight lift.
Removing equations while forgetting that lift does not remove the problem.

## J4. A finite global squareclass obstruction

More generally let Qᵢ(T,Z) be integral binary forms of positive even degree
2dᵢ with nonzero pairwise homogeneous resultants. Retain constant factors:
only rational-square scalings preserve native squareclasses. Put
R=∏ᵢ<ⱼRes(Qᵢ,Qⱼ)≠0 and S={p:p divides R}. Take a primitive integer
pair (T,Z) with all Qᵢ(T,Z) nonzero and ∏Qᵢ(T,Z) an integer square.

**Theorem J4 (finite collision target).** Every valuation of each Qᵢ(T,Z)
outside S is even. At p∈S its odd-valuation mask is an even subset of

\[
I_p([T:Z])=\{i:Q_i(T,Z)=0\pmod p\}.
\]

All native roots are rational if and only if all Qᵢ(T,Z)>0 and the first
m−1 values have even valuation at every p∈S. The last value follows from
the product. Signs and the **separate** masks at every prime determine
the entire global squareclass tuple.

**Proof.** Outside S no two forms vanish at the same projective point
modulo p, by the resultant criterion. At a primitive pair at most one
valuation is positive; its parity must be even because their sum is even.
At any prime odd valuations are supported on the vanishing forms and
have even sum. A nonzero rational number is a square exactly when it is
positive and every prime valuation is even. Homogeneous evaluation differs
from affine evaluation by Z^(2dᵢ), a square. At infinity use the homogeneous
chart directly. ∎

If Aₚ contains the allowed masks at p and A∞ the allowed sign vectors,
the number of possible global squareclass tuples is at most

\[
|A_\infty|\prod_{p\in S}|A_p|.
\]

This is a necessary ambient set, not a count of realized rational points
or Selmer classes. Defects at distinct primes cannot cancel in Q*/Q*².
For the positive quartets |A∞|=1; the bounds are 262,144 and 18,874,368,
despite the three-coordinate isogeny target.

**Factorization-free corollary.** Repeatedly divide |Qᵢ| by its gcd with
|R|, recomputing the gcd on the remaining quotient until it is one.
Write |Qᵢ|=aᵢbᵢ, with aᵢ supported on S and bᵢ coprime to R.
Then bᵢ is an integer square, and the lift criterion is positivity and
aᵢ square. This supplies portable obstruction representatives without
factoring R or the values. It need not be faster than direct square tests.

**Boundary.** A branch value requires a separate normalization chart.
A rational point on H in the prescribed class is still required.
Even valuations of a merely Qₚ-valued product point do not remove nonsquare
units; J4 does not assert a local-to-global principle for H.

## J5. Which square tests are necessary?

At each p∈S list every cluster of factors vanishing at the same
projective Fₚ-point. Join every pair in a cluster to form a graph on
the m labels. Suppose positivity has separately been established.

**Theorem J5 (collision graph compression).** If I is a vertex cover,
checking that Qᵢ(T,Z) is a rational square for i∈I proves all values square
at a rational product point. Among fixed subsets of tested indices,
vertex covers are exactly those that rule out all nonzero masks permitted
by these cluster bounds.

**Proof.** An undetected odd mask is a positive even subset of untested
indices in one cluster. It contains an edge with both endpoints untested,
which a vertex cover forbids. Conversely an uncovered edge itself is an
allowed even mask invisible on I. This converse concerns the mask model,
not global realizability. Apply J4 at every prime and retain positivity. ∎

For sparse graphs this can improve the usual m−1 tests. For a complete
graph the minimum is m−1. All three retained quartets have complete graphs,
and each pair mask has an independently verified local witness. Thus a
smaller fixed subset of individual tests cannot suffice using this local
information alone. A smaller *global* image on H(Q) is still possible.

The local witness principle is also general. At p≥5 suppose exactly two
forms have a simple affine collision r, their first-order lifts
aᵢ+bᵢs have distinct roots, and all other values are units. Put
L=bᵢbⱼ∏ₖ≠ᵢ,ⱼQₖ(r,1). Exactly

\[
(p-2-\chi(L))/2>0
\]

values of s modulo p give valuation two and square leading unit for the
product at t=r+ps. Hensel supplies a Qₚ product point with odd mask {i,j},
which cannot lift natively. Indeed the quadratic character sum for two
distinct linear factors is −χ(L); remove its two zeros and solve for
the number of +1 values. The identity ∑χ(u(u−1))=−1 follows by counting
the p−1 solutions of x²−y²=1. Zero mask requires a separate witness.

## J6. Short trace cosets force pair solubility, with a dimension cost

Let M be a positive definite lattice with minimum μ and quadratic
height h. In a coset modulo 2M, vectors v,w not equal up to sign satisfy

\[
\boxed{h(v)+h(w)\ge4\mu.} \tag{J6}
\]

**Proof.** Both v+w and v−w are nonzero in 2M and have height at least
4μ. The parallelogram identity gives 2h(v)+2h(w)≥8μ. ∎

A height λ<2μ vector is therefore unique up to sign among vectors of that
height in its coset. A known height α excludes a distinct height β when
α+β<4μ. The strict inequality matters: (1,1) and (1,−1) in the square
lattice attain equality in the same parity coset.

In the marked R17 geometry μ=4. For a proper intersection of distinct
Q-defined native bisections the established class formula is

\[
B_{w_i}\cdot(S-B_{w_j})=\tfrac12h(2S-w_i-w_j)-2.
\]

A norm-six R in wᵢ+wⱼ+2M gives S=(R+wᵢ+wⱼ)/2 and intersection length
one. A proper zero-dimensional Q-scheme of length one is Spec Q,
forcing a rational intersection point. If it lies on a smooth nonbranch
fibre, its native branches satisfy Pᵢ+Pⱼ=S; J2 leaves at most one
quotient direction. Singular and branch intersections need separate handling.
Norm-six representatives are unique up to sign; norm four excludes norm six.

This positive mechanism does not locate its point at a preassigned
parameter. The uniform panel has 131 globally eligible pairs among 519
co-split pairs, but only 18 incidences at the supplied parameters.
At −2/377 all six globally rational pair components miss the specified
fibre. The mechanism could instead construct candidate parameters from
generic intersections, followed by nonsingularity and rank checks.

For triples satisfying the existing disjointness and pair-sum image-class
hypotheses, degree=h(R)+2 with R in the triple trace-sum coset. A nonzero
coset gives degree≥6 here. This is the **total proper intersection length**,
not a lower bound on each component degree. The degree-12 example with
factorization 1+11 refutes that stronger inference. The geometric class
formulas are proved in the
[native carrier calculation](MINIMAL_CARRIER_AND_RATIONAL_SPLITTING_OF_A_TWO_DIRECTION_BLOCK.md)
and checked across the panel in the
[degree barrier](DEGREE_ONE_RELATIONS_DO_NOT_EXPLAIN_THE_LARGEST_JUMPS.md).

## Search actions and evidence requirements

These are gates and small certified constructions for future campaigns;
this theorem work launches none. The machine-readable companion is
[SEARCH_THEOREM_GATES.json](SEARCH_THEOREM_GATES.json).

| Priority | Concrete next task | Success certificate | What a miss permits |
|---|---|---|---|
| 1 | Extend the frozen triple dictionary to residual norms four and six, hence total degrees six and eight | Exact identities, complete declared coset enumeration, rational-component or branch-relation verification, and quotient relation rank | Exclude only the finite dictionary; other degrees and covers remain open |
| 2 | Compute a product curve's image in its **native affine** squareclass target | Rational product point in the zero class and exact roots, or a complete labelled obstruction | Complete obstruction excludes this carrier; incomplete Jacobian descent gives UNKNOWN |
| 3 | Use generic norm-six pair intersections to construct parameters | Exact length-one intersection, smooth rational parameter, native maps, independent quotient witness | Failed independence leaves a rational construction without the promised gain |
| 4 | Certify each block's q before increasing visibility work | Exact relations and independence, with generic subgroup and base change named | Insufficient capacity retires this mechanism for the target, not the curve |
| 5 | Freeze selection, then use separately masked known directions to measure exposure at the actual coefficient sizes | Exact inverse-chart heights and completed-coverage records | No recovery leaves exceptional-point existence UNKNOWN |

For item 1 the 2,853 existing triples and 165 existing addresses are the
reproducible retrospective population, with no additional parameters.
A future run must declare node/time limits and checkpoint each coset;
a timeout is UNKNOWN. Larger eliminations, descents and parameter sweeps
still require explicit campaign scope. Prospective performance needs a
separately frozen holdout with rank labels and exceptional coordinates
excluded until selection is fixed.

A further gate comes from the
[elliptic multiplicity theorem](ONE_AUXILIARY_POINT_HAS_A_MULTIPLICITY_BOUND.md):
one point of a positive-dimensional auxiliary curve, evaluated through
fixed maps to E, contributes at most rank Hom(J_C,E) modulo the full
generic group. Do not apply this to a finite parameter cover whose generic
fibre is zero-dimensional. Nor can a large Selmer space, a restricted
Cassels–Tate radical, or a Jacobian rational point replace the required
point on the embedded product curve in its native affine class.

## Verification and exact scope

[mechanism_theory.py](mechanism_theory.py) implements rational relation
ranks, signed graph bounds, exact resultants, a factorization-free lift
certificate at one supplied parameter, and collision test compression.
Adversarial tests cover inconsistent cycles, overlapping branches,
omitted signs, infinity, nonsquare constants and unproved product points.

[verify_mechanism_theory.py](verify_mechanism_theory.py) replays the two
quartet independence/point-relation certificates and three local defect
certificates through their existing independent verifiers. It recomputes
signed bounds on all 165 rows, quartet quotient losses, retained product
lifts, pair resultants, squareclass upper counts, and the new minimum-test
result. The existing Sage verifier separately checks generic identities,
the lattice minimum and all panel point relations.

Replay from the repository root:

~~~sh
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_mechanism_theory.py
python3 elliptic-curves/rank-jump/verify_mechanism_theory.py check
sage -python elliptic-curves/rank-jump/verify_degree_one_relation_panel.py check
sage -python elliptic-curves/rank-jump/triple_degree_barrier.py check
python3 scripts/render_status.py --check
~~~

The compact certificate is
[rank_jump_mechanism_theory_v1.json](../../artifacts/generated-results/elliptic-curves/rank_jump_mechanism_theory_v1.json).
It binds this note, gates, code, tests and retained inputs by hash.
This is written proof plus exact arithmetic replay, not Lean verification
or external review. No bounded enumeration proves J1–J6 in general;
the proofs above provide the quantified statements.

