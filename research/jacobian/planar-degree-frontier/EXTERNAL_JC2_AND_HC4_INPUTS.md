# External JC(2) results and the HC(4) equivalence

## Common partial degree at most eleven

Over an algebraically closed characteristic-zero field K, take P,Q in K[x,y]
with nonzero constant Jacobian and with both partial y-degrees at most 11 in
the same source coordinates. The conclusion K[P,Q] = K[x,y] gives a polynomial
inverse by expressing both x and y in P,Q. The independent coordinate x has
no degree bound. The pinned Lean declaration
`Max11DegreeRoutes.Max11PlaneKellerGenerationWithClassicalInterfaces`
has exactly two classical theorem parameters in addition to its field and
pair hypotheses. It supplies all five low-scale leaves, (4,6), (4,10),
(6,8), (6,10), (8,10), each for normalized scales 0 and 2.
The accompanying Comparator configuration names the headline declaration.
This is a source inspection of the pinned declaration and its assembly.
[Lean declaration](https://github.com/dcposch/jc2-lean/blob/6ee763cd623649a5af0517dd12d7202ce2e2ea3c/max11-partial-y/Max11UnconditionalAssemblyScratch.lean) ·
[Subdirectory scope](https://github.com/dcposch/jc2-lean/blob/6ee763cd623649a5af0517dd12d7202ce2e2ea3c/max11-partial-y/README.md) ·
[Comparator configuration](https://github.com/dcposch/jc2-lean/blob/6ee763cd623649a5af0517dd12d7202ce2e2ea3c/max11-partial-y/comparator.json)

### Nagata prime total-degree gcd interface

`PlaneKellerPrimeTotalDegreeGCDRoute` asserts generation for a Keller pair
with prime gcd of ordinary total degrees. Nagata’s revised Kyoto account,
p. 154 and Theorem 7.3(1), pp. 169–170, expressly takes K algebraically
closed of characteristic zero and proves this for determinant one.
Multiplying Q by the inverse of the nonzero determinant makes it one and
preserves both total degrees and K[P,Q]. This identifies the classical
mathematical input over the headline’s field. Lean still passes it as a
hypothesis. [Primary text](https://repository.kulib.kyoto-u.ac.jp/server/api/core/bitstreams/9ef8e868-5526-4830-b19f-543c0af09e7c/content) ·
[Lean interface](https://github.com/dcposch/jc2-lean/blob/6ee763cd623649a5af0517dd12d7202ce2e2ea3c/max11-partial-y/TotalDegreeRoutes.lean)

### Guccione–Guccione–Valqui standard-endpoint interface

`PlaneKellerStandardEndpointGCDObstruction` assumes coprime a,b > 1,
partial y-degrees (au,bu), x-degrees (av,bv) of their top y-coefficients,
total degrees (aS,bS), u < v, and a Keller pair. It concludes gcd(u,v) > 2.
After swapping x and y, the right endpoint of the first polynomial is
(au,av). Guccione–Guccione–Valqui Corollary 7.9 gives this gcd conclusion
for a standard (a,b)-pair. Definition 4.3 also imposes the negative
(1,-1)-weight condition on the *start* of the right face. The displayed
interface directly controls the end. The remaining source audit is to
derive the start condition, or produce an admissible normalization that
keeps the needed endpoint gcd, at least on the actual large-shear inputs
used in the max-11 proof. The classical implication remains conditional
on this interface in this DAG. The same corollary’s assertion about B
concerns the globally minimal counterexample degree gcd.
[Primary text, Definition 4.3 and Corollary 7.9](https://arxiv.org/pdf/1401.1784v3) ·
[Lean interface and large-shear adapter](https://github.com/dcposch/jc2-lean/blob/6ee763cd623649a5af0517dd12d7202ce2e2ea3c/max11-partial-y/EndpointGCDRoute.lean)

### Transfer to the repository’s plane families

To use the theorem on a family, exhibit an actual polynomial source frame
and polynomial target coordinates whose full components have both partial
degrees at most eleven in the same source variable. The
[exceptional immersed six-sheet fold](../monodromy-theory/KELLER_IMMERSED_FOLD_DEGREE_PAIR.md)
has restriction degrees (8,10) in its equality case, but full longitudinal
degrees (48,60) after the stated affine target change. The established
bound is at least 60 for the maximum longitudinal degree in the
normalized source frame. The transverse variable or a new source frame
would need a separate argument.

In the shared-resultant planar (12,42) setting, the existing
[momentum-graded algebra](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md)
has generators of X-degrees 12, 34, 42 and 99. Its degree-zero part is k.
Every nonconstant element of k[P,Q] therefore has X-degree at least 12.
Any polynomial target change has both components in this algebra, so the
max-11 bound cannot hold in the fixed momentum direction. A source
automorphism preserving the position pencil k[x] sends x to ax+b and X
to cX+r(x), with a,c nonzero. It preserves X-degree on each polynomial by
its leading-coefficient action, so the obstruction also holds throughout
that source class. Other source directions remain to be checked.

## Cusp-filtration and minimal-resolution comparisons

The companion cusp theorem uses B=C[x,y], w=x^p-y^q for coprime p,q ≥ 2,
and the semidegree on the unique normal form with weights
(q,p,pq-p-q) on (x,y,w). It assumes that pullback by a Keller map scales
this filtration by one integer e ≥ 1 on every nonzero polynomial. It
concludes e=1 and a diagonal map with multipliers satisfying α^p=β^q.
The producer and agent reviewer label it provisional. A transfer to the
repository’s cusp and Poisson work needs that full pullback property or
a suitable analogue of the filtration-preserving Euler derivation.
Boundary exponents alone give no such property.
[Reviewed source](https://github.com/dcposch/jc2/blob/15dcad7b1e65e8cb3d39e9aaac1d4c44063ec4fc/xmodel/uniform-cusp-filtration-rigidity-swarmHQ-root-sol56-20261107T1825Z.md)

The companion positivity theorem concerns a generically finite morphism
f:Y→P² of degree N>1, étale on the whole affine-plane source open, which
maps into the complement of a target line L. Here Y is smooth projective
and f contracts no smooth rational (-1)-curve. Put H=f*L, K=K_Y,
R=R_f for the full determinant divisor including exceptional
multiplicities, and Q=2K+3H=K+R in divisor class. The provisional report
gives R and Q nef and big, with each null locus consisting exactly of
the contracted smooth rational (-2)-curves, p_a(R)≥2, and h⁰(Y,Q)≥2.
It uses its preceding arithmetic-genus theorem and named surface
foundations. The proposed Keller application takes Y to be the actual
minimal resolution of the finite projective normalization. The
[local codifferent criterion](../../hessian/hc4-triangular-graphs/CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md)
asks for a velocity to lie in the finite-order codifferent. Connecting
that membership or the whole-source boundary constraints to these
intersection inequalities requires a comparison on this same minimal
resolution. Extra boundary blowups change the ramification divisor
and its numerical constraints.
[Reviewed source](https://github.com/dcposch/jc2/blob/15dcad7b1e65e8cb3d39e9aaac1d4c44063ec4fc/xmodel/minimal-resolution-ramification-positivity-swarmHQ-root-sol56-20261107T122700Z.md)

## Six-sheet genus-one watchlist

A user-supplied report dated 25 November describes a genus-one exclusion
in a specified six-sheet cusp setting. The exact statement and its
unavailable imports need to be matched to the
[retained directional cusp-slice hypotheses](../monodromy-theory/KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md)
before using it. This is a source-audit entry with no proof dependencies
on the reported exclusion.

## Formal HC(4) and JC(2) equivalence

Over an algebraically closed field of characteristic zero,
`HC4.StatementAudit.jc2_iff_hc4_explicit` equates polynomial inversion for
all two-variable unit-Jacobian maps with polynomial inversion for all
gradients of four-variable potentials of unit Hessian determinant. Both
conclusions explicitly require a two-sided polynomial inverse. The
26 November 2026 GitHub Actions clean reproduction succeeded at commit
`1377de427f2801e2f9f6b49ce1a07371e60dabb8`. The rebuild from source, axiom audit in the workflow, and
negative control belong to that run. The repository’s latter Palomar
submission snapshot `6b9526a446b4f2998307fb8d211436ccfbf723e3` exposes the
same statement as `HC4.Palomar.jc2_iff_hc4`. The snapshot and the dated
workflow run are distinct evidence. No local Lean replay was performed.
[Expanded Lean theorem](https://github.com/JamieLittle16/HessianConjecture4/blob/1377de427f2801e2f9f6b49ce1a07371e60dabb8/HC4/FinalEquivalence.lean) ·
[Workflow definition](https://github.com/JamieLittle16/HessianConjecture4/blob/1377de427f2801e2f9f6b49ce1a07371e60dabb8/.github/workflows/clean-reproduction.yml) ·
[Successful run](https://github.com/JamieLittle16/HessianConjecture4/actions/runs/45855370863) ·
[Palomar-facing Lean source](https://github.com/JamieLittle16/HessianConjecture4/blob/6b9526a446b4f2998307fb8d211436ccfbf723e3/Solution.lean)

Taking K=C supplies both directions of the existing DAG equivalence. The
[retained Tang structural audit](../../hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md)
and its exact degree-bounded transfers have their own scope and
verification obligations. Each conjecture remains open.
