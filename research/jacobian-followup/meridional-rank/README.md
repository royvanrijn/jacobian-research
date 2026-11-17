# Meridional rank conjecture

The [known-family audit](ORIGINAL_SIXTEEN_KNOWN_GRAPH_COVER_FAMILY.tex)
now establishes that **all sixteen original targets lie in the existing
Boileau--Jang--Weidmann graph-cover MRC family**. Closed-manifold rigidity
identifies each actual double cover with its prime infinite-group model,
so the published2018 theorem applies. These are not cases outside known
MRC families. The exact cover-rank3 and bridge4 certificates stand;
priority of those named computations is still unaudited.
[Scope proof PDF](output/pdf/ORIGINAL_SIXTEEN_KNOWN_GRAPH_COVER_FAMILY.pdf),
[finite bindings](graph_cover_family_scope_verification.json), and
[classical source scopes](graph_cover_family_source_audit.json).
The four earlier sharp simple3cover seeds require their own family audit.
Full MRC remains open; current catalogue46, no formal/external review.

The [final-five proof](FINAL_FIVE_DOUBLE_COVERS_SHARP_FOUR.tex) completes
the **original sixteen-case panel**: every target has double-cover group
rank3 and **mr = bridge = 4**. The five added direct source models are
`13a_234,242,243,244,265`; all actual marked central fillings have
nonabelian quotients, excluding the only applicable rank2 JSJ type.
[Proof PDF](output/pdf/FINAL_FIVE_DOUBLE_COVERS_SHARP_FOUR.pdf),
[independent inputs](final_five_cover_verification.json), and
[source scope](final_five_source_audit.json).
Eight full/32fragment reductions, three positive controls, five individual
uppers and six corrupted-input checks pass. No cross-knot transfer is used.
Twenty of293panel names are certified here;273have no certificate here.
These working proofs use explicit classical imports and remain
unformalized and externally unreviewed. Full MRC remains open;
named computational priority remains unestablished. Its unknown-family
checkpoint is superseded by the actual-cover result above.

The [three-piece filling proof](THREE_PIECE_DOUBLE_COVERS_SHARP_FOUR.tex)
adds **mr = bridge = 4** for `13a_233,239,253`, with double-cover rank3.
Its actual central filling has a nonabelian quotient, excluding the only
three-piece type in the classical rank2 JSJ list. Three pieces alone
would not suffice. [Proof PDF](output/pdf/THREE_PIECE_DOUBLE_COVERS_SHARP_FOUR.pdf),
[independent inputs](three_piece_cover_verification.json), and
[source scope](three_piece_source_audit.json).
This historical checkpoint gave eleven original targets; its five
remaining cases are now closed by their direct models above.
All five proposed prefix-based transfers failed and have different
determinants; no rank conclusion transfers. Full MRC remains open.
These are unformalized, externally unreviewed working proofs using
explicit classical imports; worldwide novelty is unestablished.

The [four-piece proof](FOUR_PIECE_DOUBLE_COVERS_SHARP_FOUR.tex) now gives
**mr = bridge = 4** for eight original targets:
`13a_129,135,136,137,194,196,197,198`.
Their two shared double-cover groups have rank3. Exact peripheral
models realize each group by a graph manifold with four canonical
Seifert pieces; Boileau--Weidmann's rank2 necessary classification
allows at most three. Independent normal-form/port/slope checks and
all eight individual upper4 traces pass; five further corruptions reject.
This is a written working proof using classical topology, without
formal verification, external specialist review or a worldwide novelty claim.
[Proof PDF](output/pdf/FOUR_PIECE_DOUBLE_COVERS_SHARP_FOUR.pdf),
[independent inputs](four_piece_cover_verification.json), and
[classical source scopes](four_piece_source_audit.json).
This eight-case checkpoint is retained; the three further cases above
supersede its remaining-target count.

The checkpoints below retain their original experiments; their rank
and open-target statements are superseded by these results for all16knots.

The [compact presentation audit](COMPACT_DOUBLE_COVER_PRESENTATIONS.tex)
gives exact three-generator presentations for the two shared cover groups,
with relator lengths `(9,13,15)` and `(12,15,11)`.
Both also have an exact assembly from the same two formal groups with
two generators and one relator each, through four retained port words.
Independent replay checks all moves, full arc gluing, abelianizations
`C391/C401` and nonabelian quotient order5616; nine corruptions reject.
The outer generators stall on five inner arcs, so no closed rank2 or
rank3 certificate follows. Peripheral injectivity and geometric piece
classification remain the next obligations; numerical recognition gave
no target classification. All original16 equalities remain open.
[Proof PDF](output/pdf/COMPACT_DOUBLE_COVER_PRESENTATIONS.pdf),
[compact replay](compact_cover_verification.json), and
[four-port replay](port_assembly_verification.json).

The [mutation audit](MUTATION_COMMON_DOUBLE_COVERS.tex) now proves
common double branched covers within two four-knot blocks:
`{13a_129,136,194,197}` and `{13a_135,137,196,198}`.
Six exact Tait-graph certificates pass independent corner, spanning-tree
and full cycle-space replay, with the classical graph-to-mutation and
cover-homeomorphism inputs explicitly sourced. These **eight cover-rank
questions reduce to two**. Each rank remains between2 and3; rank3 for
either shared group would settle its four knots using their individual
four-bridge certificates. No additional target is settled.
[Proof PDF](output/pdf/MUTATION_COMMON_DOUBLE_COVERS.pdf),
[independent replay](mutation_cover_verification.json), and
[primary source scopes](mutation_source_audit.json).

The [thirteen-point audit](PROJECTIVE_THIRTEEN_COVER_HOMOLOGY.tex)
now proves the Betti-gain rank bound **in every characteristic**:
`rank G >= 1+ceil((b1(H)-b1(G))/(r-1))` for an index-r subgroup.
The splitting identity can fail when the characteristic divides r;
the rank implication remains valid. A degree13/characteristic13
positive control has gain23 and detects rank3.
The projective-cover test from the new ternary quotients is complete.
All32 source actions give16 inequivalent covers on8 knots. Their
homology is `Z^2 + C2 + C882878` or `Z^2 + C2 + C3128602`.
The largest Betti gain is4, below13; direct characteristic13
augmentation cocycles also fail. This excludes these augmentation
modules **over every field**, without determining actual cover rank.
The retained maps give one quotient kernel per knot among these maps.
[Proof PDF](output/pdf/PROJECTIVE_THIRTEEN_COVER_HOMOLOGY.pdf),
[independent replay](projective13_cover_verification.json), and
[scoped structural source cues](projective13_source_audit.json).
The two recurring homology/polynomial blocks are now certified mutation
blocks by the subsequent audit above. All original16 targets and full
MRC remain open.

The [ternary affine audit](TERNARY_AFFINE_AND_COVER_CERTIFICATES.tex)
exhausts both rank-one `GL3(F3)` meridian classes on the original16.
Every admissible affine lift is principal; no sharp case results.
It produces **16 full GL3(F3) representation orbits on8 targets**,
inducing full `SL3(F3)` double-cover quotients. Natural/dual cover
modules all have `dim Z1=4`, below the required7. Known triple-trefoil
controls give explicit transitive27point affine actions of orders
151632 and648 with sharp rank4. The verified13point projective control
passes its homology-gain gate; target covers are now tested above.
[Proof PDF](output/pdf/TERNARY_AFFINE_AND_COVER_CERTIFICATES.pdf),
[complete base replay](ternary_base_verification.json), and
[ordinary cover replay](ternary_cover_module_verification.json).
All sixteen targets and full MRC remain open.

The [permutation-index gate](PERMUTATION_INDEX_GATE.tex) gives an
explicit full `A8` quotient for the already-known torus family
`T(4,60*j+15)`, with double-transposition meridians and sharp bound4.
Its meridian class differs from binary4 transvections by class size.
On the **original sixteen targets**, degree8 and degree9 actions of
this type are now completely excluded: **18,450 and6,561 transitive
four-seed orbits**, independently checked by full raw-tuple marking
and geometric relator replay. These are the only transitive degrees
where this graph bound could reach4 under the retained upper bound.
[Proof PDF](output/pdf/PERMUTATION_INDEX_GATE.pdf) and
[independent replay](double_transposition_verification.json).
This adds no target closure; enhanced cover-rank bounds, other
meridian cycle types and full MRC remain open.

The [simple-cover additivity proof](SIMPLE_COVER_RANK_ADDITIVITY.tex)
establishes `C(K#J)=C(K)+C(J)-1`, where `C` is the strongest ordinary-rank
lower bound over all connected simple branched covers. The fixed-field
homology version is additive too. Any sum of `t` members of the four
periodic families and `s` trefoils has **mr = bridge = 3*t+s+1**.
The method has a strict limit: **C(T(4,7)) = 3 < mr(T(4,7)) = 4**;
its deficit grows under repeated sums. General meridional additivity
and full MRC remain open.
[Proof PDF](output/pdf/SIMPLE_COVER_RANK_ADDITIVITY.pdf) and
[independent input replay](cover_sum_verification.json).

The [periodic-twist proof](PERIODIC_TWIST_FAMILIES.tex) now gives
**four explicit infinite families with mr = bridge = 4**. In each retained
seed PD, replace the first crossing by `1+6*N` alternating crossings.
A 75-point finite block action preserves the lower bound; the same four
outside seeds give the upper bound for every `N>=0`. Determinants grow
as `543+1776*N` or `537+1680*N`, separating members within each family.
[Proof PDF](output/pdf/PERIODIC_TWIST_FAMILIES.pdf),
[literal inputs and sample diagrams](periodic_twist_family_results.json),
[independent replay](periodic_twist_family_verification.json).
Periodic quotient labellings have prior precedent; worldwide novelty and
known-family exclusion remain unestablished. Full MRC remains open.

The [trefoil-amalgam proof](TREFOIL_AMALGAMS_SHARP_FOUR.tex) gives
**mr=bridge=4 for13a_3145,13a_3147,13a_3148,13a_3149**. Their simple
three-sheeted cover groups are nonabelian rank-two peripheral amalgams
of trefoil groups. Exact group isomorphisms, independent named PDs,
branch words and four-seed upper certificates replay successfully.
[Proof PDF](output/pdf/TREFOIL_AMALGAMS_SHARP_FOUR.pdf) and
[exact replay](three_sheet_splice_verification.json).
There are also explicit **25-point finite commutator certificates for
all four cover groups**: [literal relators and permutations](finite_cover_certificates.json),
[standalone Python checker](check_finite_cover_certificates.py), and
[independent full-cover replay](finite_cover_quotient_verification.json).
Run `python3 check_finite_cover_certificates.py` from this directory.
Each pair satisfies all six checkpoint relators; the full replay checks
all39 original branched relators and rejects four corruptions.
These strengthen the same four case proofs; they add no further knot.
These four lie elsewhere in the293-name panel; the original sixteen
remain unchanged. All26 tested three-sheeted cover-group ranks are now
determined:22 cyclic and4 rank two. Worldwide priority/known-family
exclusion remains unestablished; full MRC is still open.

The [rank-two equivalence proof](RANK_TWO_AFFINE_EQUIVALENCE.tex)
now characterizes that part of the original criterion exactly:
for a fixed S3 colouring and simple three-sheeted cover Y, a transitive
dimension-k certificate exists iff `b1(Y;F2) >= k-1`. Two distinct S4
quotient kernels above the same S3 kernel already force lower4.
The complete293-name symmetric audit has no sharp case: 21 names admit
only dimension2/bound3, and272 admit no original certificate in any
dimension. Every binary three-dimensional affine base is now excluded
as a sharp4 route, including proper GL3 images.
[Proof PDF](output/pdf/RANK_TWO_AFFINE_EQUIVALENCE.pdf) and
[independent replay](rank_two_equivalence_verification.json).
The [earlier nonlinear checkpoint](three_sheet_nonabelian_verification.json)
proved19 cyclic and retained seven unresolved presentations. Their ranks
are now resolved by the trefoil-amalgam continuation above.

The [singular affine certificate proof](GL3_QUOTIENTS_AND_SINGULAR_AFFINE_CERTIFICATES.tex)
gives explicit sharp binary certificates for the known torus family
`T(4,7*(2r+1))`. Their dependent meridian image directions explain a
precise coverage restriction in the external four-dimensional search.
All293 pinned open names were then exhausted for GL3(F2) quotient
kernels and specified affine lifts: 191 maps on154 knots, maximum7
kernels versus the required58, and no sharp root lift over F2 or F7.
[Independent replay](gl3f2_verification.json) and
[proof PDF](output/pdf/GL3_QUOTIENTS_AND_SINGULAR_AFFINE_CERTIFICATES.pdf).
Every characteristic-two coefficient module through these epimorphisms
also fails the ordinary rank-three gate. Proper images and other fields
are outside that exclusion. No new open case or full MRC proof results.

The [simple-module proof](SIMPLE_MODULE_BOUNDS_AND_A5_EXHAUSTION.tex)
now exhausts ordinary linear cocycle bounds through **every single A5
map** on the sixteen-target panel, over all fields and coefficient
dimensions, including nonsplit extensions: `dim Z1 <= 3*dim V/2`.
Every A5-induced cover degree also fails the Betti rank gates; all have
zero rational b1. This does not prove actual group rank two.
[Exact replay](a5_all_modules_verification.json) and
[proof PDF](output/pdf/SIMPLE_MODULE_BOUNDS_AND_A5_EXHAUSTION.pdf).
The four two-kernel knots have separate natural tensor tests over F4,
all below threshold; other product modules and quotients remain open.

The [meridian-deformation proof](MERIDIAN_DEFORMATIONS_AND_STABILIZATION.tex)
now computes all sixteen retained H3 points on twelve target knots.
Each has one eigenvalue-varying tangent modulo conjugation; keeping
meridians involutive gives formal rigidity. A uniform splitting theorem
also rules out irreducible higher-rank families through these points
plus any finite trivial summand, even with varying meridian eigenvalue.
Both restricted cohomology hypotheses are checked on all sixteen bases.
[Exact tangent/cohomology replay](meridian_tangent_verification.json).
All32 degree-three residue jets pass, with
[independent polynomial replay](meridian_jet_verification.json), but
actual complex integrability is unresolved. A transverse slice reduces
it to one-parameter residual identities or a finite obstruction.
[Proof PDF](output/pdf/MERIDIAN_DEFORMATIONS_AND_STABILIZATION.pdf).
No retained knot is settled; full MRC remains open.

The [actual peripheral-curve proof](TORUS_WHITEHEAD_PERIPHERAL_CURVES.tex)
gives a positive source-derived application: all iterated untwisted
Whitehead doubles of positive torus knots `T(p,q)`, `2<=p<q`, have
`augmentation rank = meridional rank = bridge = p*2^N`. Explicit
companion matrices realize Cornwell's known maximal-rank curve;
its ideal-point orders meet Gao-Gu-Li-Wang's iteration criterion.
The November2026 preprint and pattern conventions are explicit
dependencies, with no novelty claim. More generally a rank-n
augmentation on an n-strand closure with nonzero writhe suffices.
An actual `T(3,5)` curve passes through an H3 point and escapes
orthogonality by varying the meridian eigenvalue. This does not
construct curves or settle any of the sixteen retained table targets.
[Exact controls](torus_whitehead_verification.json),
[source scope](torus_whitehead_source_audit.json),
[proof PDF](output/pdf/TORUS_WHITEHEAD_PERIPHERAL_CURVES.pdf).

The [six-point and formal-orthogonality proof](SIX_POINT_COVERS_AND_FORMAL_ORTHOGONALITY.tex)
uses the symmetric-square character of the icosahedral representation.
All sixteen A5-induced index-six covers have cyclic homology `4D` or
`28D`, with no sharp rank gain. This also proves that every formal complex
deformation through these rank-three reflection points which keeps
meridians involutive preserves a varying symmetric form and keeps
`lambda^2=1`. The new deformation note above strengthens the fixed-
involution conclusion to formal rigidity and computes the varying-
eigenvalue tangent; actual integrability remains open.
[Independent replay](six_point_cover_verification.json).
[Proof PDF](output/pdf/SIX_POINT_COVERS_AND_FORMAL_ORTHOGONALITY.pdf).

The [augmentation homology proof](AUGMENTATION_COVER_HOMOLOGY.tex)
proves the user's stronger index-five Betti-gain test without knot-root
geometry. All sixteen A5 kernel classes on the current sixteen targets
give gain one at three and zero otherwise. Their abelianizations force
every connected index-five cover to occur in this panel, so the proposed
certificate fails throughout degree five. [Direct module and Fox replay](augmentation_cover_verification.json).
The note also proves that orthogonal rank-one bases have longitude
`lambda^2=1`, excluded by the Whitehead doubling step. The
[original-source dependency audit](whitehead_dependency_audit.json)
checks Dutra's p7 subgroup contradiction and its p12 use; the affected
2022 theorem is not an established-family dependency.

The [joint-kernel proof](JOINT_SIMPLE_QUOTIENT_KERNELS.tex) proves the
user's twenty-A5-kernel gate. All 293 candidates have at most five each;
the known pretzel `P(2,3,3,5)` has 23, giving sharp four where raw Hom
gives only three. [Separate matrix replay](joint_a5_verification.json).

The [quadratic-coefficient proof](ICOSAHEDRAL_LIFTS_AND_QUADRATIC_COEFFICIENTS.tex)
closes the fixed twelve H3 bases: every affine lift factors through the
120-element group in every field. It compares restricted/full cocycle
bounds on the same representation and recovers the original
nondegenerate conditional bound. [Replay](icosahedral_lift_verification.json).

The [branched-cover proof](BRANCHED_COVER_LIFT_BOUNDS.tex) identifies the
root-lift bound as `r-1 + dim H1(Y; F)` for the simple degree-r branched
cover associated to a transposition quotient. It also proves
`meridional rank >= r-1 + rank pi1(Y)` and gives a nonabelian wreath-lift
count. The degree-two version is the classical Kelly-Wada core-group
method, with historical credit retained. These are working proofs;
formal verification, external review and general priority audits remain
open. [Proof PDF](output/pdf/BRANCHED_COVER_LIFT_BOUNDS.pdf).

The [all-characteristic lift proof](ALL_CHARACTERISTIC_LIFTS.tex) gives
`meridional rank >= dim Z` over every field. Integer Smith certificates
find five-dimensional lift spaces over F5 for all 16 earlier five-seed
table diagrams, proving `meridional rank = bridge = 5` on those inputs.
The [separate replay](integer_lift_verification.json) checks all 512 lift
matrices and supplies [explicit sharp certificates](sharp_five_certificates.json).
The [named-diagram and family audit](diagram_audit_verification.json)
identifies all 16 as five-tangle Montesinos knots, already covered by
Boileau-Zieschang (1985). A separate rational-tangle constructor and exact
planar-map bijections verify these identifications. They are controls,
not newly resolved cases.
[Proof PDF](output/pdf/ALL_CHARACTERISTIC_LIFTS.pdf).

The [lift-dimension proof](LIFT_DIMENSION_AND_RADICAL_TOWERS.tex) now proves
`meridional rank >= dim Z` for the entire space of signed lifts of a fixed
symmetric transposition colouring. For an S4 colouring of dimension `ell`,
it constructs a transitive affine action on `2^(ell-1)` points preserving
a rank-two alternating form. The first new target has order 1536 on 16
points; the original nondegenerate dimension-four criterion is already
covered by sharp A5 or D5 quotients. [New proof PDF](output/pdf/LIFT_DIMENSION_AND_RADICAL_TOWERS.pdf).

A [frozen panel](lift_panel.json) of 16 table cases with five-seed upper
certificates and no reported sharp finite Coxeter quotient gives binary
bounds four; the F5 extension above reaches five. Both the source extraction and the complete
S3/S4 colouring and lift calculations have independent replay in
[the new verification receipt](lift_dimension_verification.json).

The full conjecture remains open. The [working proof](PROOF.tex) constructs
explicit sharp affine representations for a known pretzel family, extends
the quadratic bound to nonzero degenerate alternating forms, and gives an
exact linear test for signed lifts of symmetric transposition quotients.
It also proves that the entire affine criterion cannot apply to
`T(p,q)` when `p >= 3` and `q` is prime with `q > 2^(p-1)-1`.

The [literature map](LITERATURE.tex) and [45-source catalogue](sources.json)
include the November 2026 Whitehead-double preprint, its notice of a gap
in earlier work, and separate unrefereed computations. Audit depth and
unfinished source checks are explicit. No new unresolved knot is claimed
settled, and no novelty claim is made for these reformulations.
All 243 October matrix certificates now pass a
[separate algebraic and upper-bound replay](current13_certificate_verification.json).
Their [named diagrams match the independent census](current13_name_audit.json),
and [both archived KnotInfo snapshots](knotinfo_snapshot_audit.json) confirm
all 243 PDs and the author's complete 536-case uncertainty partition.
A [frozen 13-crossing panel](current13_panel.json) now has an
[all-degree, all-field exclusion](branched_cover_verification.json):
four seeds force symmetric degree at most five; degrees three through
five have no base colouring, and degree two has maximum lift dimension
two in every characteristic. The original binary common-form affine
criterion is impossible in every dimension on these targets, since
`k+1 <= bridge <= 4` forces a rank-two form and the absent S3 quotient.
The earlier suggestion of higher-rank binary forms on this panel is
withdrawn. A [separate A5 core-count replay](nonabelian_core_verification.json)
also gives no sharp four-bound, using the same sixteen targets.
The [index-five cover audit](core_subgroup_homology_verification.json)
uses the twelve retained first nonconstant A5 witnesses. All those
covers have cyclic homology over the integers, so none reaches the
Schreier threshold of seven homology classes for a sharp four-bound.
The later augmentation audit above extends this to all sixteen kernel
classes and every connected index-five cover of the same base groups;
other degrees and nonabelian subgroup information remain untested.
The [odd-characteristic SL2 root-lift panel](parabolic_panel_verification.json)
also gives no sharp four-dimensional lift space over F3, F5 or F7. The
general restricted-cocycle bound and these finite negative scopes are
written in the all-characteristic proof; other fields remain untested.
The [cross-programme scope review](CROSS_PROGRAMME_SCOPE.json) checks the
radical-tower centralizer lemma and retains the extra geometric hypotheses
of the proposed Keller exclusion; the knot representations remain valid.

- [Proof PDF](output/pdf/PROOF.pdf)
- [Literature PDF](output/pdf/LITERATURE.pdf)
- [Finite verification receipt](verification.json)
- [Pinned primary inputs](sources/manifest.json)

From the repository root, using Python 3.10 or later:

```sh
python3 research/jacobian-followup/meridional-rank/verify_signed_lifts.py
python3 research/jacobian-followup/meridional-rank/verify_integer_lifts.py
python3 research/jacobian-followup/meridional-rank/verify_current13.py
python3 research/jacobian-followup/meridional-rank/verify_diagram_audit.py
python3 research/jacobian-followup/meridional-rank/verify_parabolic_panel.py
python3 research/jacobian-followup/meridional-rank/verify_branched_cover_bounds.py
python3 research/jacobian-followup/meridional-rank/verify_nonabelian_core.py
python3 research/jacobian-followup/meridional-rank/verify_core_subgroup_homology.py
python3 research/jacobian-followup/meridional-rank/verify_icosahedral_lifts.py
python3 research/jacobian-followup/meridional-rank/verify_six_point_covers.py
python3 research/jacobian-followup/meridional-rank/verify_core_generator_probe.py
python3 research/jacobian-followup/meridional-rank/verify_joint_a5.py
python3 research/jacobian-followup/meridional-rank/verify_augmentation_covers.py
python3 research/jacobian-followup/meridional-rank/audit_whitehead_dependencies.py
python3 research/jacobian-followup/meridional-rank/signed_lifts.py research/jacobian-followup/meridional-rank/example_pretzel_5.json
```

The example presents a published Coxeter quotient of a pretzel group;
the knot conclusion uses the published meridian-preserving epimorphism.
For a general input, the code checks the proposed symmetric quotient and
every relator, but the knot identity, meridian labels and bridge upper
certificate require their own provenance. Finite checks support the
written proofs and are not formal or external verification.

The sixteen sharp cases are now identified controls in a known family.
The next gate is a certified nonabelian cover-rank bound or another
representation beyond the excluded routes, with independent knot and
prior-family checks. Derive a rank-drop or group certificate before
enlarging a coefficient search. The retained negative panels should not
be repeated as discovery. Exact many-twig weights/tree conditions and
the Lustig-Moriah numerical plat criterion are now extracted; original
geometric and recent full-proof audits remain unfinished.
The mathematical status is registered as `MRC-AFFINE-SIGNED-LIFTS-20261116`,
`MRC-LIFT-DIMENSION-RADICAL-20261116`, `MRC-ALL-CHARACTERISTIC-LIFTS-20261116`,
`MRC-BRANCHED-COVER-LIFTS-20261116`, and
`MRC-QUADRATIC-COEFFICIENT-LIFTS-20261116`,
`MRC-JOINT-SIMPLE-QUOTIENT-KERNELS-20261116`, and
`MRC-AUGMENTATION-COVER-HOMOLOGY-20261116`, and
`MRC-FORMAL-ORTHOGONALITY-20261116`,
`MRC-TORUS-WHITEHEAD-CURVES-20261116`,
`MRC-MERIDIAN-DEFORMATIONS-STABILIZATION-20261116`,
`MRC-SIMPLE-MODULE-A5-EXHAUSTION-20261116`,
`MRC-GL3-SINGULAR-AFFINE-CERTIFICATES-20261116`,
`MRC-RANK-TWO-AFFINE-EQUIVALENCE-20261117`,
`MRC-TREFOIL-AMALGAM-SHARP-FOUR-20261117`,
with `OP-MRC` still open, in
[MATH_STATUS.json](../../MATH_STATUS.json).
