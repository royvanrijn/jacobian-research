# Dixmier conjecture bounds and research connections

Reviewed through 14 November 2026. The generalized Dixmier conjecture fails
for every Weyl rank at least three. The remaining ranks are one and two.
Our retained results explain several failures of proposed rank reduction
and quantization routes; they leave both questions unresolved.

This requested follow-up does not reactivate the archived search queue.
The [status registry](../../MATH_STATUS.json) remains authoritative for local
claims. Literature statements below retain their stated hypotheses and
preprint status. The exact calculation here reproduces classical formulas
and tests their uncorrected normal ordering; it makes no novelty claim.

[Concrete experiments and connections with the newer research](NEXT_EXPERIMENTS.md)
adds an ordering pilot and corrects the earlier next-step list using QRF2.

## Rank and degree bounds

Over a characteristic-zero field, write

\[
A_n=K\langle x_1,\ldots,x_n,\partial_1,\ldots,\partial_n\rangle,
\qquad [\partial_i,x_j]=\delta_{ij}.
\]

The assertion \(DC_n\) says that every unital endomorphism of \(A_n\) is
surjective. Injectivity follows already from simplicity. Generator images
must satisfy every Weyl relation; those identities alone do not prove that
the images generate the ambient algebra.

| Question | Current boundary | Consequence for our work |
|---|---|---|
| \(DC_1\) | Open in the reviewed literature; Zheglov v5 claims a proof | Keep that claim separate from established inputs |
| \(DC_2\) | Open in the reviewed literature | Rank-two classical maps still need a valid quantum construction |
| \(DC_n\), \(n\ge3\) | False by the inverse-Jacobian lift of the three-variable Keller collision | Reuse the completed lift rather than search for its existence |
| First false Weyl rank | At most three; one or two is not excluded | Assuming \(DC_1\) would narrow this to two or three |
| Possible \(A_1\) counterexample \((P,Q)\) | \(\gcd(\deg P,\deg Q)>15\) in the Guccione–Guccione–Valqui theorem | This is a gcd bound for total degrees, not an exhaustive finite search bound |

The last row is the precise result of
[Guccione, Guccione and Valqui, 1111.6100v3](https://arxiv.org/abs/1111.6100v3).
It implies that both total degrees exceed 15 for a counterexample, but
provides neither a finite upper bound nor a termination criterion for search.
This review does not assert that it is the strongest bound in every
normalization or degree convention.

[Han, Pan and Chen, July 2024](https://arxiv.org/abs/2407.11291v1) exclude
Weyl mates for the stated Joseph normal forms when either normal-form
integer is prime, extending the known coprime case. Their integers are
normal-form parameters, and must not be substituted for arbitrary total
degrees in the preceding gcd statement.

The dimension-sensitive implications remain

\[
JC_{2n}\Longrightarrow DC_n\Longrightarrow JC_n.
\]

[Belov-Kanel and Kontsevich](https://arxiv.org/abs/math/0512171v2), following
Tsuchimoto's independent 2005 result, establish the first implication.
Thus \(JC_2\Rightarrow DC_1\), while \(DC_2\Rightarrow JC_2\).
The failure of \(JC_4\) does not decide \(DC_2\): a false sufficient
condition says nothing about the conclusion. Stable equivalence does not
give same-rank equivalence or a downward transfer of counterexamples.

## Recent developments

Alpöge's July 2026 three-variable example, and Gao's account revised
4 November, give actual Keller collisions. Gao also constructs examples
of arbitrarily large geometric degree in each dimension greater than two.
See [Gao, 2608.00222v2](https://arxiv.org/abs/2608.00222v2) and our
[retained symplectic and Weyl lift](../../archive/non-elliptic/extended-geometry/SYMPLECTIC_WEYL_LIFT.md).
The old all-dimensional conjecture is therefore false, while its low-rank
questions survive.

[Long's rank-two Poisson construction](https://arxiv.org/html/2608.23777v1)
has four polynomial outputs \((R,T,D,S)\), preserving two canonical pairs,
with \(R=x(2-3xq)\). The reducible first output prevents polynomial
invertibility. Its differential-operator lift acts on four position
variables and consequently produces an \(A_4\) endomorphism. This does not
produce an \(A_2\) endomorphism.

[Kanel-Belov and Fedorov, 24 October 2026](https://arxiv.org/html/2610.33919v1)
prove \(JC_2\Rightarrow DC'_1\), with \(DC'_1\) equivalent to invertibility
of endomorphisms of \(\operatorname{Mat}_2(A_1)\). Their converse requires
existence of specified characteristic-zero graph-module lifts. Theorem 1.8
characterizes lifting along a faithful arithmetic family by bounded
Berest–Wilson invariants; Conjecture 1.9 proposes a degree-dependent uniform
bound. Theorem 4.2 excludes lifts of Long's prescribed Frobenius center map
at every matrix size. This is an arithmetic lifting obstruction, distinct
from a finite-order symbol obstruction. These are preprint results.

[Zheglov v5, revised 19 January 2026](https://arxiv.org/abs/2410.06959v5)
continues to claim ordinary \(DC_1\). Kanel-Belov and Fedorov explicitly
report difficulty understanding its proof, and retain \(DC_1\) among the
remaining questions. A proof claim, ordinary \(DC_1\), and generalized
\(DC'_1\) must therefore remain separate.

## Exact classical and quantum comparison

The new [Fraction arithmetic checker](verify_long_pbw.py) independently
transcribes Long's formulas. It checks all six Poisson identities and three
distinct rational points mapping to \((0,1/8,0,0)\). It then interprets each
expanded polynomial in normal order \(x^a q^b p^c z^d\), with
\([p,x]=[z,q]=\hbar\), and computes the complete PBW commutators.

| Pair | Powers in the nonzero commutator defect | Nonzero terms after setting \(\hbar=1\) |
|---|---|---:|
| R,T | 2,3 | 15 |
| R,D | 2,3 | 90 |
| R,S | 2 | 3 |
| T,D | 2,3,4,5 | 578 |
| T,S | 2,3 | 88 |
| D,S | 2,3,4,5 | 359 |

Here the defect subtracts \(\hbar\) times the required Poisson constant.
For example, although \(\{R,S\}=0\), the normal-ordered images satisfy

\[
[R,S]=-6\hbar^2x^6(3qx-2)(9qx-4)\ne0.
\]

This short identity has a second derivation: because R depends only on
positions and S has momentum degree two, normal multiplication gives
\([R,S]/\hbar^2=3qS_{pp}+6xS_{pz}\); direct differentiation yields the
displayed factorization. A separate SymPy calculation reproduced it.
Thus literal normal ordering fails. Corrections or a different construction
require further arguments; this calculation does not exclude them.
The subsequent symmetric-order pilot removes this R,S defect exactly,
while five other relation defects survive; see the linked experiments.

For comparison, the existing rank-three checker passes: if \(J=DF\) and
\(B=J^{-1}\), then

\[
\Phi_F(x_i)=F_i(x),\qquad
\Phi_F(\partial_i)=\sum_r B_{ri}(x)\partial_r
\]

satisfies all Weyl relations. Differential order is preserved exactly,
because the associated cotangent pullback is injective. Hence a hypothetical
preimage of each multiplication operator \(x_i\) would itself have order
zero and supply a polynomial inverse to F, contradicting the collision.
The [full written proof](../../archive/non-elliptic/extended-geometry/SYMPLECTIC_WEYL_LIFT.md)
also proves persistence on adjoining identity Weyl pairs.

## Connections with our results

| Retained result | What it establishes | Remaining boundary |
|---|---|---|
| Q1 and Q3 | A quartic classical completion and exact obstructions for specified corrections | Other supports, orderings and filtrations remain possible |
| QO1, QO2 | Degree-five fixed-symbol and low-support obstructions | They do not cover every possible Weyl endomorphism |
| QO10 | A nonempty open degree-six family is obstructed at order five in its declared filtration | Odd corrections, wider filtrations and other polarizations remain outside scope |
| QRF1 and QRF2 | Relative residues and an exact order-seven obstruction on the full reduced cubic exceptional component | Boundary charts, other filtrations and correction classes remain separate |
| DC2HL3 | The R21 Hamiltonian locally nilpotent locus is exactly K[R] | Formal Darboux triviality does not supply a polynomial chart |
| DC2QB1 | R21 has no strict filtered PBW Weyl frame; its graph centralizer is exactly the transported A2 | An abstract identification with collapsing filtration remains outside this theorem |

Sources are the [quartic note](../../archive/non-elliptic/extended-geometry/QUARTIC_WEIGHTED_GEOMETRY.md),
[filtered obstruction note](../../archive/non-elliptic/extended-geometry/RANK_TWO_FILTERED_QUANTIZATION_OBSTRUCTION.md),
[degree-six family](../../archive/non-elliptic/extended-geometry/DEGREE_SIX_RELATIVE_QUANTIZATION_OBSTRUCTION.md),
[quantum residue note](../../archive/non-elliptic/extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md),
and [R21 sections 19–21](../../archive/non-elliptic/extended-geometry/DC2_HIGHER_NILPOTENCE_R21_FRONTIER.md).
The receipt pins ten current claim scopes and their canonical source hashes.
Their retained proofs were read for scope, without rerunning the large archived campaigns.

The connection is a shared failure of passage between categories: formal
solutions, finite Laurent solutions, polynomial solutions, strict filtered
quantizations and arithmetic lifts impose different requirements. Our residue
classes cannot be identified with a Berest–Wilson invariant without a new
comparison theorem. R21 and Long's map are also distinct objects; the new
nonlifting theorem cannot be applied to R21 solely because both have
rank-two classical geometry.

[Lauritzen and Thomsen, 2308.09384v2](https://arxiv.org/abs/2308.09384v2)
prove that finite generation of A over an endomorphism image on either side
forces surjectivity. Combining that theorem with the verified rank-three
counterexample shows that its image cannot make A3 a finite left or right
module. Finite generic degree of a Keller map concerns function fields;
it does not give this module finiteness. The paper's automatic bimodule
finiteness is a different property.

## Concrete next questions

1. For a possible DC2 counterexample, require one explicit four-generator
   assignment satisfying every Weyl relation, followed by a proof of
   nonsurjectivity. Neither a classical collision nor a bounded failure to
   recover generators supplies both requirements.
2. For the plane lifting route, construct graph splitting modules for
   specified finite-field reductions and measure their Berest–Wilson
   invariants with certified presentations. A finite prime sample cannot
   prove a uniform bound; an automorphism control and field/degree
   hypotheses must accompany such an experiment.
3. Use QRF2's completed characteristic-zero obstruction rather than
   reopening OP-QO's superseded cubic exception. The full replay artifact
   is absent locally; its scalar Bezout identity passes a fresh check,
   while full replay remains unavailable. Degree nine has a separate
   unresolved characteristic-zero section-zero reconstruction.
4. For R21, any new bridge must address the precise strictness or ambient
   identification hypothesis. Repeating formal Darboux corrections or
   restricting to its exhausted graph centralizer cannot cross that gate.

## Reproduction and verification limits

From the repository root:

```sh
.venv/bin/python research/jacobian-followup/dixmier-audit/verify_long_pbw.py
.venv/bin/python research/archive/non-elliptic/scripts/verify_symplectic_weyl_lift.py
```

[receipt.json](receipt.json) records the new exact coefficient checks;
[rank_three_replay.json](rank_three_replay.json) retains the small archived
checker's output and source hash. These checks have no search component.
They do not reprove the full October preprint, verify Zheglov's proof, or
independently reprove the archived R21 and residue theorems.

The repository-wide research command currently fails before this audit on
the missing historical source
`artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json`
in METHOD-EC-CENTRE-POLICY. Local links and certificates are checked separately;
the absent unrelated evidence is not reconstructed.
