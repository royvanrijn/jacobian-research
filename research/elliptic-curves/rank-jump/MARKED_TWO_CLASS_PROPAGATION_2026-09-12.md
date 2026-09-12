# Marked two-class propagation: theorem and bounded carrier experiment

The simultaneous-splitting implication is valid with the hypotheses below.
The prescribed first experiment has **not constructed a marked rank18
subfamily**: the retained prefix contains182 exact smooth geometrically
rational bisections, all nonsplit above `t0=3/17`. There are73 reducible
residuals and one unresolved interrupted compilation. Both lattice shells
remain incomplete at the frozen two-million-node limit. This is a finite
visibility result, not an exclusion of the proposed construction.

The carrier branch is now **closed at that stop**. The
[final incidence diagnostic and arithmetic redirection](CARRIER_CLOSURE_AND_DEPENDENCY_CONTINUATION_2026-09-12.md)
find 182 distinct extensions and record the translation/pairing obstruction.
The next compatibility check concerns a complete principal-atom dependency;
no shell, candidate or Riemann–Roch computation is resumed.

The [earlier blind recovery](BLIND_CONSTRUCTED_CLASS_RECOVERY_2026-09-12.md)
already proves rational solubility of the two fixed constructed covers.
Solubility at this one fibre does not produce the parameter curves required
here. Rational-versus-Sha controls remain applicable.

## 1. Verified simultaneous-splitting theorem

Let `E/Q(t)` be non-isotrivial, and let `G1,...,Gr` be rational sections.
Let `t0` be a smooth rational fibre with no rational2-torsion. Suppose
`delta(G1(t0)),...,delta(Gr(t0)),beta6,beta7` are independent in
`E_t0(Q)/2E_t0(Q)`, where the last two classes are rationally soluble.

Let `Ci` be geometrically integral degree-two multisections defined over Q,
whose smooth projective normalizations have genus zero. Supply exact
presentations and point maps

```
Ci: si^2=di(t),     Pi(t,si) in E(Q(t,si)),
```

with `di` a nonconstant squarefree polynomial of degree1 or2, including
its rational constant factor. Require `di(t0)` to be nonzero rational
squares. At chosen lifts `ci,0`, require exact cubic-field identities

```
delta(Pi(ci,0)) / (beta6^ai beta7^bi product_j delta(Gj(t0))^eij) = xi_i^2,
```

where the two rows `(ai,bi)` have rank2 over F2. In the present application,
the cubic is irreducible, so the usual etale cubic algebra is a field.
These identities require actual square roots, not finite-character agreement.

Let `B` be the smooth projective normalization of the component of
`C1 x_(P1_t) C2` through the chosen lifts, and require it to be geometrically
integral and dominate the t-line. Write `b0` for its rational point over
the chosen pair and `tau:B->P1` for the finite nonconstant map. Then:

1. The displayed `r+2` sections are Z-independent in `E(Q(B))`.
2. If the quadratic extensions over the **fixed t-line** are equal, the
   selected component is a genus-zero graph component, rational over Q.
3. Otherwise the extensions are geometrically independent. If `b` is the
   number of geometric branch points in their union, `g(B)=b-3`. Thus
   `b=3` gives genus zero and `b=4` gives genus one.
4. If `B` has genus zero, or genus one with a chosen rational origin and
   a certified rational nontorsion point, then infinitely many rational
   parameters in `tau(B(Q))` have rank at least `r+2`. Only finitely many
   rational points of B, and finitely many values in its image, can fail
   the specialization conclusion.

**Proof.** The two point maps become sections after base change. At b0 the
invertible label matrix and exact generic corrections give `r+2` independent
Kummer classes. To see that these imply Z-independence, divide the coefficients
of a hypothetical relation by their maximal common power of2. The resulting
point sum is2-primary torsion, hence zero because `E_t0(Q)[2]=0`. Reduction
modulo2 then contradicts independence. Any relation between the sections
would specialize to a relation at b0, so the sections are independent.

If the extensions are equal, their tensor product has two graph components;
the chosen unramified rational lifts select one. If they are distinct but
become equal over Qbar(t), their ratio is `c*r(t)^2` with `c in Q*`.
Evaluation at t0 makes c a rational square. Thus this apparent third case
is impossible under the simultaneous-splitting hypothesis. This also explains
why dropping a rational constant twist would be an error.

In the independent case the compositum is a geometrically connected degree4
V4 cover. Each geometric branch point has inertia order2, including a point
ramified in both quadratic subcovers, and contributes2 to ramification.
Riemann--Hurwitz gives `2g(B)-2=-8+2b`. Each double cover has two branch
points; independent covers share at most one. Infinity is a branch point
precisely for the degree1 presentations.

A genus-zero curve with a rational point is Q-rational. A genus-one curve
with a rational nontorsion point has infinitely many rational points.
Non-isotriviality survives the finite nonconstant base change. Silverman's
specialization theorem gives injectivity outside bounded height; Northcott
finiteness turns this into finitely many rational exceptions. Finite fibres
of tau then give infinitely many distinct rational t. The application is
standard specialization theory; see [Schuett--Shioda, §14.7][SS] and the
original [Silverman reference][S].

The original Silverman publisher page was not retrievable in this review;
the specialization statement was cross-checked in the accessible primary
survey. The proof above is not an independent reproving of that theorem.

For MW16-05, r=16. Once actual labelled carriers are supplied, the two square
conditions describe the image of the selected common base, with finite
exceptional images excluded. No such pair is supplied by the present run.

If the original fibration has exact generic rank16, two additional independent
sections over the unchanged Q(t) are impossible. Any rational linear relation
with the sixteen sections can be cleared of denominators and specialized at
the smooth control, contradicting the prescribed Kummer block. A nontrivial
base cover or isolated specializations are therefore essential. Opposite
branches of one bisection have a common generic trace and supply one quotient
direction; they cannot be used as the two independently labelled carriers.

## 2. An effective strengthening: congruences that certify every output

There is a distinction between an explicit pair of square conditions with
an **unspecified finite exceptional set**, and an explicit condition which
certifies each output. Silverman's qualitative theorem alone gives the former.
The following sufficient refinement can give the latter once a pair exists.

Suppose the18 specialized points at b0 have a finite-reduction certificate:
their images in a finite product `product_p E_t0(Fp)/2E_t0(Fp)` are independent.
Also choose a good odd prime at which the two-division cubic has no root.
Spread the family and its18 sections over suitable local integral models.
For every certificate prime choose an explicit p-adic neighbourhood `Up`
of b0 on B such that the reduced equation and all18 reduced points are the
same as at b0. Such neighbourhoods exist by continuity; rational-function
formulas with denominators prime to p give explicit residue conditions.
At a chart pole, use a regular projective chart and sufficient p-adic precision.

Every rational `b in intersection_p Up` then has the same independence
certificate and no rational2-torsion. Consequently

```
b in B(Q) and b in every certified Up  ==>  rank E_tau(b)(Q) >=18.
```

This condition has no uncomputed rank-exception set. It does not assert that
the old finite primes already work on an as-yet-unconstructed common base.
The models, neighbourhoods and reduction identities still require certificates.

The restriction retains infinitely many points. On a rational B, move b0
to parameter0 and choose an integer multiple M small enough p-adically that
`u=M*z`, `z in Z`, lies in all Up. On an elliptic B take b0 as origin. Each
Up contains an open subgroup of the compact group `B(Qp)`, hence a subgroup
of finite index. A common multiple N puts a fixed nontorsion point Q into
all of them; every `b=[Nk]Q` lies there. N can be certified using local
reduction/formal-group calculations, including bad reduction of B. It is
not necessary to assume good reduction of B at every certificate prime.

The inherited mod2 control certificate can be used only after the chosen
point representatives and their exact generic corrections have been
transported. Independently choosing signs or generic corrections at
different primes does not supply a single global point or correction.

## 3. Full integral lattice and equation compiler

Use native11952 coordinates `F_old,O_old+F_old,M17`, with Gram
`U + (-M17)` and determinant948. The sanitized MW16-05 trace word w has
norm8, and its actual fibre is `(2,2,w)`. The stored equation-side zero
is an old section. The fibre at old `lambda=infinity` is exactly
`O_old+P_w`; intersection with the new zero identifies its nonidentity
I2 component. Taking the integral orthogonal complement and reducing it
gives a rank17 positive frame of determinant948, with two roots ±R.
Both coordinate transports are integral and unimodular. Projecting the
sixteen marked sections away from U and R recovers the stored height Gram
of determinant474. The root and integral glue have not been replaced by
that rank16 height matrix.

For `D.F=2`, `D.O=h`, write `D=2O+(h+4)F+v`. Direct intersection gives
`D^2=4h+8+v^2`, hence frame norm `4h+10`. The two prescribed shells are
therefore10 and14, with both signs of every vector retained. Translation
deduplication was not used.

The equation step adapts the existing
[residual-chord compiler](../../elkies-k3/scripts/construct_elkies_2026_bisections.sage)
and [trace RR interpolation](../cas/construct_curve302_parent_cheapest_lattice_bisection.sage).
There is a useful explicit A1 extension. Put `m=D.R`, `epsilon=m mod2`.
For the generic trace T, its height and zero intersection satisfy

```
height(T)=norm(v)-m^2/2,
c=T.O=(height(T)-4+epsilon/2)/2.
```

The marked class identity is

```
D + S_(-T) = 3O + (h+c+6)F - ((m+epsilon)/2)R.
```

For surviving vertical-wall candidates, `m in {0,1,2}`. Write `k=h+c+6`.
The ambient polynomial channels `f0+f1*x+f2*y` have degree bounds
`(k,k-4,k-6)`. Vanishing on `-T` is a linear rational-function identity.
When R must be subtracted, impose first-order vanishing at the physical
I2 node, using the infinity chart when required. Eliminating y and removing
the known trace gives an exact quadratic in x. Its discriminant is retained
as a rational square times a squarefree branch polynomial, with its constant.

An irreducible quadratic with degree1/2 branch gives an actual geometrically
rational horizontal curve C. The line construction gives `D=C+V` for an
effective vertical V. The wall checks give `D.V>=0`, while `C.V>=0`.
Adjunction gives `C^2>=-2`; hence `-2=D^2=C^2+C.V+D.V` forces equality.
Then `V^2=0`, so V is a sum of whole fibres, which would have positive
intersection with C unless V=0. Thus C has class D and `C^2=-2`, and is
smooth by adjunction. A lattice vector by itself was never called a carrier.

## 4. Frozen finite experiment and its boundary

The protocol uses one worker and8 GiB RSS, a four-CPU-hour aggregate ceiling,
two million enumeration nodes or60 CPU-minutes,128 compilation candidates
per shell and30 CPU-seconds per candidate. No shell enlargement occurred.

The exact integer LDL traversal uses the existing
[integer lattice arithmetic](../cas/visibility_lattice_fast.py). It stops
after2,000,000 nodes in the norm-at-most14 ellipsoid, retaining16,356 norm10
and562,451 norm14 vectors. Both shells remain incomplete. Known effective
walls reject4,775 and190,150 respectively. These are individual exact
intersection rejections, not a complete chamber test.

The128 candidates per shell are the first128 **within that retained prefix**
under the frozen generic RR-cost proxy: native fibre degree, native zero
intersection, absolute-coordinate sum, then native coordinates. They are
not claimed to be the first128 of either completed shell. No control
splitting, exceptional coordinates, recovered point representatives, V3
clouds or class labels enter this ordering.

The first export failed on a Sage integer in the JSON kernel dimension.
Its171 failed receipts and source are preserved. The replacement reuses the
same arithmetic, converts Sage integers for JSON, and subtracts prior CPU
from each original30-second allowance. The interrupted172nd child has no
complete usage receipt, so it is charged its full30 seconds and left UNKNOWN.
It is not repeated. No fresh30-second budget or new candidate was introduced.

The separate verifier reconstructs trace sums by its own rational group-law
code, rechecks the marked divisor identities and all coefficient maps,
and checks exact rational nonsplitting. The resulting endpoints are:

| Stage | Outcome |
|---|---|
| Full integral frame | rank17, determinant948, exact U/A1 marking |
| Retained candidates |256, from two incomplete shells |
| Exact smooth geometrically rational bisections |182 |
| Split/reducible residual quadratics over Q(t) |73 |
| Unresolved interrupted compilation |1 |
| Completed carriers splitting at3/17 |0 |
| Nonzero marked labels |0; no candidate reaches square-identity matching |
| Eligible pairs, infinite common bases, fresh fibres |Not reached |
| Strict and ordinary ideal-class transfer |UNKNOWN; no fresh fibres |

Phase B therefore stops at its exact splitting prerequisite. No prime-character
match was promoted to a class label, and no independent points on unrelated
fibres were substituted. Phases C and D, the local/ideal audit, and V3 comparison
are not run without their prerequisite pair.

For each of the182 carriers the replay also retains a small-prime nonsquare
witness with integral branch coefficients at that prime. If M is the product
of those distinct primes, every carrier is nonsplit throughout the explicit
progression `t=3/17+M*z`, `z in Z`. The exact M and witnesses are in the replay
certificate. This strengthens the finite visibility exclusion only. It says
nothing about rank drops, other bisections, or the existence of the marked
block by some other mechanism on that progression.

Here `M=1217269416735962936106761`, the product of
`7,11,13,19,23,29,31,37,41,43,47,53,59,61,67,79`.

The prescribed stops are substantive: this run does not complete the norm10
or norm14 shells, does not exclude the one unresolved candidate, and does
not address higher zero-intersection curves or singular rational curves
of positive arithmetic genus.

## 5. Strictness and ordinary ideal classes remain separate

For a fresh fibre with independent18-dimensional Kummer span H and generic
subspace G of dimension16, let loc be the single combined map to the complete
local squareclass spaces at all bad primes,2 and the real place. Then

```
k_strict = dim ker(loc|H)-dim ker(loc|G)
         = 2-(rank loc(H)-rank loc(G)).
```

Indeed, the strict directions visible in `H/G` are the image of `ker(loc|H)`;
the kernel of that projection is exactly `ker(loc|G)`. Since `loc(G)` is a
subspace of `loc(H)`, `k_strict=2` precisely when their images are equal.
The correcting generic word must be the same at every place.

At good odd primes, multiplication by2 on the abelian scheme is finite
etale. Pulling it back along a rational point gives an unramified Kummer
torsor. This justifies good-prime unramifiedness in a compatible integral
model; every additional prime introduced by the actual model/transport must
be accounted for. New bad primes cannot be replaced by the control's list.

The strict half-ideal map has a unit-squareclass kernel. Two independent
strict Kummer classes therefore still require the direct half-ideal/Artin
verification to establish two additional ordinary ideal-class directions.
The [existing block verifier](CONSTRUCTED_CLASS_BLOCK_AND_RATIONAL_LIFTS.md)
is the relevant positive reference. A successful finite panel would not
prove uniform strictness or ideal independence on an infinite set.

## 6. Provenance, replay and next gate

The family and control block remain retrospectively selected. Geometry uses
an allowlisted equation/marking packet with source snapshots and hashes;
there is no claim of operating-system isolation. Frozen compaction columns6
and7 retain their full1,676/1,572 atom dependencies and square-equivalence
circuits as provenance. The label prerequisite failed before any new field
square-root computation was needed.

The [compact result](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/result.json),
[182 explicit carriers](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/carriers.json),
[independent replay](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/independent-replay.json)
and [frozen block provenance](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/frozen-block.json)
retain the exact mathematics. The [12.3 MB evidence archive](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/evidence.zip)
and [492-member manifest](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/evidence-manifest.json)
include completed and failed work; every member hash was checked.
The archive SHA256 is
`d9a8244d51aa1c66590ba8acce959cbfd7a66609c592a828124548d93970de55`.

Recorded arithmetic plus the conservative interruption charge totals269.613
CPU seconds: marking0.061, enumeration25.902, compilation children169.202,
export-repair parent1.948 and independent replay72.500. Process startup,
the interrupted first parent's overhead, source retrieval, authoring,
packaging and navigation are not separately component-timed. This is not an
end-to-end benchmark. All eighteen small regression tests pass, including
finite ranks3,3,2,1 at the four prescribed anchors, unequal-ordinate and
duplication controls, constant twists, shared branch points at infinity,
and the need for one global local-correction word.

Run the narrow independent replay in an empty restored replay checkout:

```sh
sage -python research/elliptic-curves/rank-jump/verify_marked_two_class.sage
sage -python research/elliptic-curves/rank-jump/test_marked_two_class.sage
```

The replay's default output is created exclusively; preserve any prior output
before an explicitly requested replay. It does not enumerate another shell,
repeat equation discovery, compute a class group, or search for points.

The low-degree existence proposition remains open, but this carrier branch
is closed. The [next arithmetic gate](CARRIER_CLOSURE_AND_DEPENDENCY_CONTINUATION_2026-09-12.md)
returns to the successful principal-relation construction: define a complete
dependency with a family-level continuation and specialization maps, then
test its compatibility. The soluble classes do not imply the low-shell
geometric hypothesis. No larger shell or interrupted-candidate repair is
scheduled.

The theorem is a certificate-oriented application of established descent,
elliptic-surface and specialization theory. [Salgado's related rank-gain
constructions][Sa] have rational-surface hypotheses which are not assumptions
proved here for the production K3. [Cremona--Fisher--Stoll][CFS] supply the
relevant genus-one minimization and reduction algorithms when an actual
genus-one common base is reached; small coefficients alone never certify
a rational point or positive rank.

[SS]: https://arxiv.org/html/0907.0298v3#S14.SS7
[S]: https://doi.org/10.1515/crll.1983.342.197
[Sa]: https://arxiv.org/pdf/1307.5912
[CFS]: https://arxiv.org/pdf/0908.1741
