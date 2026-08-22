# Exact elliptic-neighbour compiler

[`scripts/elliptic_neighbor_compiler.sage`](scripts/elliptic_neighbor_compiler.sage)
is the strict equation-level core for a fibration hop. Given an ambient
Riemann--Roch basis and maps to finite quotients of *resolved* fibre charts,
it stacks exact vertical conditions, computes their rational rank and (when
requested) materializes their rational kernel, and permits an `h0(D)=2`
conclusion only when the chart cover is explicitly marked complete. It never
infers a chart from a Kodaira symbol or component label. Dimension-only mode
is deliberate: it certifies the nullity without forcing Sage to expand a
huge rational nullspace before a later simplification has chosen a pencil.
The same core now supplies the marked-chord discriminant, exact squarefree
binary-quartic extraction (including the scalar-square check that rejects an
accidental quadratic twist), classical Jacobian invariants, and marked-section
pencil evaluation for degree-two hops.  Its resolved marked-chord adapter
also compiles a caller-supplied local expression `a+b*m` in an explicit
finite quotient ideal; this is used by the q6-child II*/IV* jet blocks and
works for both multivariate chart quotients and univariate jet rings.
`compile_resolved_chart_cover` collects a declared set of those quotient
charts and explicitly derived overlap rows into the one matrix from which an
`h0(D)=2` conclusion may be made; declaring a complete cover remains a
separate geometric certificate.

The first H3 lattice search calls its shortest orbit a `q6` shell. After the
recorded Weyl reduction, however, the divisor is `D=O+(-P1)-F` and `D.F=2`.
Thus `q6` is a search label, not the degree used in the Riemann--Roch bound.
The generic degree-two space requires the marked chord
`(y-y(P1))/(x-x(P1))`; substituting the five monomials of `L(5O)` would silently
construct the wrong problem.

Run the core regression and H3 preflight with:

```sh
sage -python elkies-k3/scripts/verify_elliptic_neighbor_compiler.sage
sage -python elkies-k3/scripts/compile_h3_first_q6_preflight.sage
```

The preflight writes
[`../artifacts/generated-results/elkies-k3-h3-q6-compiler-preflight.json`](../artifacts/generated-results/elkies-k3-h3-q6-compiler-preflight.json).
It verifies the exact divisor data and marked `P1` input, and records the
22 recorded Weyl reflections plus the E7/E8 affine and simple-wall pairings.
It also records the actual E7 quotient block and marked-point module. Its
status is `PASS_EXACT_Q6_ACTUAL_E7_LOCAL_INPUTS`. The actual eight-blow-up E8 chart
tree is checked by
[`scripts/derive_h92_q6_e8_resolution.sage`](scripts/derive_h92_q6_e8_resolution.sage),
while [`scripts/derive_h92_q6_e8_p1_branch_module.sage`](scripts/derive_h92_q6_e8_p1_branch_module.sage)
derives the complete local E8 module `u<1,Q>` in the common `F_infinity`
representative. The four smooth P1.O collision modules are certified by
[`scripts/derive_h92_q6_smooth_po_module.sage`](scripts/derive_h92_q6_smooth_po_module.sage).
The preflight deliberately stops at these local inputs; the complete resolved
RR cover and the first-hop certificate below supply the vertical-matrix,
child, and transport claims.

Important correction: [`scripts/audit_h92_q6_actual_e7_marked_chord_order.sage`](scripts/audit_h92_q6_actual_e7_marked_chord_order.sage)
substitutes the raw chord into the actual `E7_2--E7_5` chart and proves
`ord_Z(m/t)=-1` at the generic point of `E7_5`. The corrected calculation
uses `Z*m/t=unit/W`, and the new all-edge certificate proves the raw chord
has nonnegative exceptional order on every actual E7 component. Thus `<1,m>`
is the resolved q6 module cover: generator `1` away from `-P1`, generator
`m` at that marked point. The rerun cover has rank eight on the ten-term
ambient, certifying `h0(D)=2`; its eliminated child and Weyl transport now
certify `E8+E6/MW3` and the stated rank-three Gram.

The historical global saturation is reproduced by

```sh
sage -python elkies-k3/scripts/assemble_h92_q6_global_rr.sage
sage -python elkies-k3/scripts/derive_h92_q6_p1_actual_e7_quotient_block.sage
sage -python elkies-k3/scripts/certify_h92_q6_actual_resolved_rr_cover.sage
sage -python elkies-k3/scripts/eliminate_h92_q6_global_pencil.sage
sage -python elkies-k3/scripts/certify_h92_q6_child_jacobian.sage
sage -python elkies-k3/scripts/derive_h92_q6_e7_resolution.sage
sage -python elkies-k3/scripts/derive_h92_q6_e7_valuation_atlas.sage
sage -python elkies-k3/scripts/derive_h92_q6_third_e7_local_target.sage
sage -python elkies-k3/scripts/assemble_h92_q6_third_generic_rr_ambient.sage
sage -python elkies-k3/scripts/derive_h92_q6_third_e7_cartier_charts.sage
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_resolution.sage
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_resolution_full.sage
sage -python elkies-k3/scripts/derive_h92_q6_third_actual_e7_cartier_charts.sage
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_chart_pullbacks.sage
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_valuation_atlas.sage
sage -python elkies-k3/scripts/trace_h92_q6_p1_actual_e7.sage
sage -python elkies-k3/scripts/derive_h92_q6_p1_actual_e7_marked_module.sage
sage -python elkies-k3/scripts/reject_h92_q6_third_old_monomial_e7_ideal.sage
sage -python elkies-k3/scripts/derive_h92_q6_third_actual_e7_quotient_block.sage
sage -python elkies-k3/scripts/evaluate_h92_q6_third_marked_chord_actual_e7_quotient.sage
sage -python elkies-k3/scripts/evaluate_h92_q6_third_generic_ambient_actual_e7_quotient.sage
sage -python elkies-k3/scripts/evaluate_h92_q6_third_e7_point_series.sage
sage -python elkies-k3/scripts/certify_h92_q6_third_e7_chord_units.sage
sage -python elkies-k3/scripts/derive_h92_q6_child_zero_section.sage \
  --component-output artifacts/generated-results/elkies-k3-h92-q6-child-e7-infinity-sections.json
sage -python elkies-k3/scripts/certify_h3_q6_component_section_lattice.sage
sage -python elkies-k3/scripts/certify_h3_q6_weyl_section_transport.sage
sage -python elkies-k3/scripts/certify_h3_q6_actual_neighbor_hop.sage
sage -python elkies-k3/scripts/lift_h92_p2_hensel.sage \
  --precision 1024 \
  --output artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json
sage -python elkies-k3/scripts/verify_h92_p2_coordinates.sage \
  --input artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json
sage -python elkies-k3/scripts/transport_h92_q6_third_modp.sage
```

It produces a 10-dimensional ambient, an exact collision matrix of rank and
codimension 8, and a two-dimensional kernel. Its ratio has squarefree
old-base degree 4. The Jacobian is the displayed minimal short Weierstrass
model in the generated child artifact; its reducible fibres are `II*+IV*`,
so its root lattice is `E8+E6` and its Mordell--Weil rank is three when
combined with the rank-19 source Néron--Severi certificate. Transported
old-zero coordinates are also certified: `D.O=1`, its new-base map is
`u=T`, and the binary-quartic covariant map yields an exact child section in
the displayed Jacobian. The rank-three lattice transport and its exact Shioda
Gram are certified; minimized equation-level coordinates for every member of
that basis remain open.

The resolved-cover certificate replays its finite smooth collision block
through `compile_resolved_conditions` with a complete cover and exact rank
calculation. It uses dimension-only mode because the already pinned global-RR
artifact is the canonical materialized two-function basis; this avoids
duplicating its enormous rational nullspace in a second backend.

The final hop certificate now normalizes that resolved-RR evidence together
with the child equation and the marked-section transport through
`certify_exact_neighbor_hop` in
[`scripts/elliptic_neighbor_compiler.sage`](scripts/elliptic_neighbor_compiler.sage).
The shared gate checks rank-nullity and complete-cover status, exact child
root data, section identities, all three new-fibre degrees, and the full
Shioda Gram.  It records `nef_on_declared_walls`, not an unproved assertion
that those declared walls exhaust the effective cone; this is the reusable
equation-level handoff for later neighbours.

`compile_elliptic_neighbor_rr_pencil` is the reusable lattice-to-pencil
orchestrator.  Starting with the raw isotropic class, it replays the ordered
Weyl record, checks the declared nef walls of the resulting primitive class,
stacks caller-supplied resolved-chart blocks, and binds a displayed pencil to
the same matrix.  It stops there: a raw lattice vector cannot supply the
vertical correction, chart trivialization, or a child equation.  The H3 q6
preflight now replays its 22 reflections through this interface before the
surface-specific chord conversion.

Before choosing that ambient, `certify_generic_fibre_divisor_decomposition`
now binds the caller-supplied marked section and vertical correction to the
actual class identity `D=(q-1)O+P+V+kF`.  It checks every declared vertical
component has old-fibre degree zero and refuses a reconstruction mismatch.
For the first q6 hop this records the exact identity
`D=O+(-P1)-F`; it makes the marked chord input explicit, but does not claim
to discover `P` or the vertical correction from a lattice vector alone.
Given that certified input, `marked_section_generic_fibre_basis` records
`L((q-1)O+P)` as the ordinary basis of `L((q-1)O)` plus a separately labelled
marked chord.  Thus q6 has the generic-fibre basis `(1,m)` rather than a
monomial surrogate; the routine deliberately leaves all base powers and
resolved vertical conditions to the surface-specific chart calculation.

For a natural representative with repeated marked sections,
`certify_generic_fibre_horizontal_support` instead records the literal
effective horizontal support before any group-law rewrite.  The first q8
source uses this to certify `9O+9(-P1)+V-11F` in the Néron--Severi frame;
this preserves the required 18-element chord-power basis and prevents an
unrecorded replacement by a one-point divisor from changing the ambient.
For a certified quadratic marked-chord extension,
`balanced_marked_chord_power_basis(n)` supplies the checked `2n` functions
for `L(nO+nP)`; the q8 source uses its `n=9` specialization.

`exact_neighbor_certificate_handoff` also accepts the versioned output of
the repository's degree-independent Néron--Severi engine. It retains the
engine's precise supplied-wall boundary, fixed-component record, old-fibre
degree, root data, and saturated MW height before an equation-level ambient
is attached. In particular, it rejects a chain/search label such as `q24`
when the certificate's actual old-fibre degree is `2`.

`certify_component_pairing_transport` is the companion final-stage guard for
declared component-pairing data.  It accepts only caller-supplied source-curve
identities, degree-one intersections with the new fibre, a named ordered
component basis, and exact pairing rows.  It pins those entries and checks
the displayed horizontal-plus-vertical degree balance.  In particular it
does not turn an E6/E7/E8 label into a component correspondence.  The first
q6 replay records the Néron--Severi old-E7 pairing rows of the third correction
and its `4812-4811=1` balance.  This adapter does not assert that a lattice
pairing row is a new resolved-chart trace.

For every degree-two marked-chord hop, the preceding equation conversion is
now one core operation: `compile_degree_two_chord_hop` validates the marked
point on the old Weierstrass fibre, solves the two-generator pencil for the
chord, retains the exact square-factor check while extracting the binary
quartic, computes the Jacobian invariants, and evaluates explicitly supplied
marked chords on the new parameter.  The H3 elimination script calls this
operation, while the pinned child-Jacobian calculation is its regression
model; minimization and local Kodaira classification remain surface-specific
certificate stages.

`compile_resolved_degree_two_child_jacobian` is the stricter equation hand-off:
it first binds the chord conversion to the exact two-row resolved RR kernel,
then performs finite-place minimization and Kodaira aggregation over the new
base.  Its infinity data remain a boundary diagnostic until a caller supplies
the corresponding line-bundle trivialization and resolved additive-fibre
charts.  The q6 child-Jacobian replay uses this route directly.

More specifically, `compile_resolved_degree_two_chord_hop` accepts the
complete resolved condition matrix, its independently displayed two-row
kernel, and one `(a,b)` chord expansion for every ambient coordinate.  It
first verifies that the displayed pencil is annihilated by that very matrix,
then derives both `a_i+b_i*m` generators before invoking the conversion.
The first H3 q6 elimination now replays its pinned 8-by-10 collision matrix
through this bridge and reproduces its pinned genus-one artifact hash.

For higher old-fibre degree, `eliminate_cleared_weierstrass_pencil` is the
degree-independent first elimination stage. Given the old Weierstrass
relation and two common-denominator pencil representatives in one declared
affine polynomial ring, it returns the full ideal obtained after eliminating
only the named old variables. An optional saturation requires the explicit
list of resolved-chart denominator factors; the unsaturated ideal, product,
and saturated ideal are retained in the result. It does not choose one
relation from a multi-generator ideal, compute a Jacobian, or minimize a
model: each of those actions needs its own resolved-chart provenance. This
provides a safe hand-off for the q8, q24, and later genus-one adapters.

`compile_resolved_genus_one_elimination` makes that hand-off strict.  It
combines the two rows of a certified complete RR kernel with caller-supplied
common-denominator representatives in one affine ring, then invokes the
elimination stage.  Consequently a higher-degree relation cannot be taken
from an unrelated displayed pencil.  It still retains the complete
elimination ideal: choosing a plane cubic, a rational origin, a Jacobian, or
a minimal model requires the separate exact adapters below.

When a retained elimination relation is a pointed plane cubic rather than a
binary quartic, `pointed_plane_cubic_to_weierstrass` implements the exact
non-flex normalization used elsewhere in this repository.  The caller must
provide the cubic point and its two successive third tangent intersections;
the core verifies all three points and the resulting Weierstrass equation,
then returns the rational cubic-to-Weierstrass map.  Limits at a map base
point must still be taken in the caller's resolved chart, and minimalization
or identification with a target child remains an explicit later certificate.

`finite_minimize_short_weierstrass` performs the finite-place part of the
next stage for a short child model.  It records every irreducible base factor,
the raw and minimized `(ord(A),ord(B),ord(Delta))`, and the exact
fourth/sixth-power scaling unit.  It separately reports the infinity orders
instead of falsely declaring a global minimal model from an affine
calculation; the appropriate projective line-bundle scaling and fibre
classification must still be certified for each surface.

`kodaira_data_from_short_orders` converts a minimal characteristic-zero
short-Weierstrass valuation triple into its Kodaira symbol, root rank, Euler
number, and root determinant.  It is intentionally not a component-label
oracle: child components and section intersections are accepted only from the
transported resolved geometry.

`classify_finite_short_weierstrass_fibres` combines those two finite-place
steps into one auditable root-signature record: irreducible factor, degree,
raw and minimized valuations, Kodaira symbol, and aggregate finite root rank,
Euler number, and determinant.  It reports the normalized infinity valuation
only as a boundary record.  Thus it cannot silently turn an affine
minimization into a global fibre or a section-component assertion.

`resolved_chart_overlap_condition` is the corresponding local-gluing block.
It accepts two finite chart-quotient evaluators and two explicit maps to one
common overlap quotient, then imposes equality there.  Thus the forthcoming
q8 node residual can use the actual completed-chart transition data without
substituting a rectangular corner jet or a fibre-type heuristic.

`transport_binary_quartic_point_to_jacobian` supplies the next hand-off for
an already transported resolved divisor: it checks `w^2=f` at the supplied
binary-quartic point, evaluates the classical `(H,G)` covariants, and checks
the fourth/sixth-power minimalizing change on the child coordinates.  It does
not look for a section on the child equation.  The existing q6 zero and
E7-infinity transports are its surface-specific regression model.

The last command pins the two low-height component sections that must be
transported next: old `E7_7` and old affine `E7` each meet `D` once. Their
Mordell--Weil difference gives the exact leading height block
`[[8/3,1/3],[1/3,8/3]]`; it is a component-label and intersection certificate,
not yet their child-equation coordinate transport. Its optional
`--component-output` now transports both points at old-base infinity through
the binary-quartic covariant map to exact child-Jacobian coordinates. The
E7 chart replay now resolves the sign: the affine E7 component is `plus`,
and `E7_7` is `minus`. Its seven ordinary blow-ups and four final nodes are
recorded explicitly, so this is not an inference from the III* label.

The Weyl-aware transport applies the recorded reflections to the three
canonical child-frame lifts before reading old Mordell--Weil coordinates. It
therefore gives actual degree-one divisor classes in the nef model, not raw
complement coordinates. Their projections to the old marked MW quotient are
`4*(-P1)`, `(-P1)`, and `22*(-P1)-P2`, and their exact Shioda Gram is
`[[8/3,1/3,-1],[1/3,8/3,1],[-1,1,46]]`. The p-adic Newton lift recovers all
139 normalized P2 coefficients exactly and the verifier proves its H92
Weierstrass identity, with pole profile `(46,42)` for `x` and `(69,63)` for
`y`. Thus the remaining coordinate task is a prescribed divisor-aware
transport, not a reconstruction or search for another child section.

The final command is a deliberately negative regression guard. Evaluating the
third old MW-*projection* alone modulo `100003` creates old coordinates of
degrees `4996,7494` and gives q=6-pencil degree `4769`, not one. The H92
coordinate has the frame sign `-P2`; this calibration agrees with the
Néron--Severi pairing. The Weyl certificate now identifies the full third
divisor as `43*O + (22*(-P1)-P2) + V`, with
`V=(-22,-33,-44,-66,-55,-44,-33;0,...,0;-2389)` in the old
`E7,E8,F` vertical basis. Its vertical contribution is `-4811`, reducing
the horizontal degree `4812` to one. This guards exactly against the invalid
“evaluate a large MW multiple” shortcut.

[`scripts/derive_h92_q6_e7_valuation_atlas.sage`](scripts/derive_h92_q6_e7_valuation_atlas.sage)
supplies a formal-model regression atlas for this correction. It
derives the old-coordinate orders `(2,2,4,3,1,2,3)` of `Z` on
`E7_1,...,E7_7`, with the corresponding orders of `U,Y`, directly from the
seven *formal* blow-up charts. In particular it accounts for the nonreduced
first three intermediate exceptional curves. It is not a transported H92
atlas and cannot support an H92 Riemann--Roch claim by itself.

The follow-on
[`scripts/derive_h92_q6_third_e7_local_target.sage`](scripts/derive_h92_q6_third_e7_local_target.sage)
performs that divisor-class bookkeeping in the formal resolved numbering. It proves
that the E7 exceptional cycle is
`(-22,-44,-66,-44,-33,-33,-55)=22*c_q6`, with component-degree vector
`(0,0,0,0,22,0,0)`. This is an exceptional-cycle identity only: the integral
vertical divisor still requires a transported H92 chart trivialization and
P2's smooth-chart condition remains separate.

[`scripts/assemble_h92_q6_third_generic_rr_ambient.sage`](scripts/assemble_h92_q6_third_generic_rr_ambient.sage)
then writes the exact 44-dimensional generic ambient as 43 monomials of
`L(43O)` plus the marked chord for `P=22*(-P1)-P2`. The high-height point is
recorded as a group-law expression DAG, evaluated in each resolved chart by
the core's `evaluate_marked_point_dag` helper. This avoids a global
`QQ(t)` GCD normalization which is neither needed nor suitable for preserving
chart trivializations. The positive source height `5260` certifies that the
marked chord is nondegenerate.

[`scripts/derive_h92_q6_third_e7_cartier_charts.sage`](scripts/derive_h92_q6_third_e7_cartier_charts.sage)
turns the integral E7 cycle into six formal-model node-chart conditions. For
example, the E1--E4 chart is `Y^2-Z-U*Z^2=0`, with local factor `U^44Y^22`;
all six E7 edges are recorded with their own smooth chart equation and
Cartier factor. They are useful regression data, but are not H92
trivializations until an explicit coordinate transport is exhibited.

[`scripts/derive_h92_q6_actual_e7_resolution.sage`](scripts/derive_h92_q6_actual_e7_resolution.sage)
is the required safety gate. It replays the initial blow-ups on the exact H92
germ and proves that the purported second-Z-chart node is smooth (its
`Z`-derivative is `-B1`, which is nonzero). The direct replacement is now
[`scripts/derive_h92_q6_actual_e7_resolution_full.sage`](scripts/derive_h92_q6_actual_e7_resolution_full.sage): its fourth node is in the second-U
chart at `Z=-A1/B1`, and together with the three third-stage nodes it gives
the seven actual H92 exceptional curves and all six E7 edge charts. Pulling
the degree-44 divisor into those charts is now done by
[`scripts/derive_h92_q6_third_actual_e7_cartier_charts.sage`](scripts/derive_h92_q6_third_actual_e7_cartier_charts.sage),
which records all six actual Cartier factors. Evaluating the marked chord in
these charts requires the exact global-to-local maps now certified by
[`scripts/derive_h92_q6_actual_e7_chart_pullbacks.sage`](scripts/derive_h92_q6_actual_e7_chart_pullbacks.sage): each map pulls the H92 equation back to the edge-chart
equation times its recorded total-transform factor. The corresponding marked
point is evaluated as an exact 96-jet Laurent-series point by
[`scripts/evaluate_h92_q6_third_e7_point_series.sage`](scripts/evaluate_h92_q6_third_e7_point_series.sage).
[`scripts/certify_h92_q6_third_e7_chord_units.sage`](scripts/certify_h92_q6_third_e7_chord_units.sage)
then certifies that its chord is a unit at every actual edge-chart origin.
Selecting finite quotient jets remains the next matrix step.

[`scripts/derive_h92_q6_actual_e7_valuation_atlas.sage`](scripts/derive_h92_q6_actual_e7_valuation_atlas.sage)
reads the valuations of the original `(t,x,y)` directly from those pullbacks.
It obtains the fibre-multiplicity vector `(2,2,4,3,1,2,3)` and checks the
actual E7 Cartan relation, furnishing the input for complete-ideal quotient
construction without invoking a normal-form atlas.

The atlas also rules out a common but invalid shortcut for the degree-44
third marked divisor. The actual anti-nef E7 cycle has complete-ideal
colength 363, whereas the ideal generated solely by old-coordinate monomials
with sufficiently large valuation has colength 529; see
[`scripts/reject_h92_q6_third_old_monomial_e7_ideal.sage`](scripts/reject_h92_q6_third_old_monomial_e7_ideal.sage).
The correct actual translated coordinate is instead
`U=x-c2*t^2-c3*t^3`: the third cycle is exactly `11*v(U)`, and
[`scripts/derive_h92_q6_third_actual_e7_quotient_block.sage`](scripts/derive_h92_q6_third_actual_e7_quotient_block.sage)
constructs its 23-generator, length-363 H92 quotient. Thus raw coordinate
orders are insufficient, while the transported local coordinate gives the
finite quotient required for the high-degree condition block.
[`scripts/evaluate_h92_q6_third_marked_chord_actual_e7_quotient.sage`](scripts/evaluate_h92_q6_third_marked_chord_actual_e7_quotient.sage)
now evaluates the exact DAG chord in that quotient and verifies
`(x-x(P))*m_-P=y+y(P)` there, including its 33-step nilpotent denominator
inverse. The 44 generic degree-44 functions are then evaluated by
[`scripts/evaluate_h92_q6_third_generic_ambient_actual_e7_quotient.sage`](scripts/evaluate_h92_q6_third_generic_ambient_actual_e7_quotient.sage):
the exact E7 map has rank 33 and local kernel dimension 11.  Its entire
local kernel is now materialized without a dense rational nullspace call:
it is the eleven individually vanishing monomials
`x^17,...,x^21,y*x^15,...,y*x^20`, as certified by exact rank-nullity.
The block is replayed through `compile_resolved_conditions`, whose sparse
zero-column path produces that coordinate kernel. This is only the E7 block,
not a global high-degree pencil.

The first q=6 marked branch has also been transported to this same actual
resolution by
[`scripts/trace_h92_q6_p1_actual_e7.sage`](scripts/trace_h92_q6_p1_actual_e7.sage):
both `P1` and `-P1` pass through the nonzero second-U node and meet `E7_5` at
distinct smooth points away from the `E7_2--E7_5` edge. This replaces a
component-label inference with an exact blow-up trace. The old marked-frame
script stated that `m/t` was regular in the `E7_5` chart, but the exact audit
rejects that statement. The corrected derivation
[`scripts/derive_h92_q6_p1_actual_e7_marked_module_corrected.sage`](scripts/derive_h92_q6_p1_actual_e7_marked_module_corrected.sage)
proves `Z*m/t=unit/W`; because `t/Z` is a unit, the raw q6 frame `<1,m>` is
the marked frame. The all-edge certificate
[`scripts/derive_h92_q6_actual_e7_all_edge_module.sage`](scripts/derive_h92_q6_actual_e7_all_edge_module.sage)
then proves nonnegative exceptional chord orders on every actual E7 component.
The resolved RR cover, child Jacobian, and transport certificate therefore
complete the first equation-level q6 hop.

For the next q=8 hop, the endpoint data now has an exact smooth-collision
algebra rather than a Kodaira-form guess.  In
[`scripts/derive_h92_q8_smooth_collision_frame.sage`](scripts/derive_h92_q8_smooth_collision_frame.sage),
put `h=Z4`, `q=(m-y(P1)/x(P1))/h`, and `X=h^2*x`.  After clearing the chord
equation and removing `X-h^2*x(P1)`, the residual equation is monic quadratic
in `X` with coefficients regular at `h=0`.  The resulting 18-element
`q/X` frame is only a local algebra frame: the q=8 divisor lattice and its
smooth quotient matrix have intentionally not been inferred from it.

The all-component generic E7 layer is now also an executable compiler block:
the 22 forced singleton coefficients and 983 actual-chart residue-cancellation
rows stack to a 1005-by-54 exact matrix of rank 54. Hence the least q8
endpoint ambient has no generic-E7-compatible direction, independently of
the smooth obstruction. This rejects that finite ambient only; a genuine q8
compiler must enlarge the Riemann--Roch space before applying its node and
overlap conditions.

<!-- status-consumer: EC-K3-H3-Q8-GENERIC-COMPONENT-CONDITION-BLOCK 33693f196eb13091 -->

The first smooth-surviving enlargement (`r=7`, 558 columns) is rejected too:
the smooth block has mod-43 rank 556, and the complete actual generic E7 row
set has rank two on its two-dimensional kernel. The resulting full-column
good reduction proves that this enlarged smooth-plus-generic block has zero
characteristic-zero kernel. It is an ambient rejection, not a statement about
larger ambient spaces or the missing node/overlap compiler.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA7-SMOOTH-GENERIC-REJECTION e0e67861f23d4f24 -->

The next smooth-surviving enlargement (`r=8`, 630 columns) is rejected by the
same actual-chart interface. Its smooth block has mod-43 rank 624, while the
complete generic E7 layer has rank six on its six-dimensional smooth kernel.
The full-column good reduction is an exact QQ rank certificate for this one
ambient, not a replacement for the missing node/overlap compiler or a q8
pencil.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA8-SMOOTH-GENERIC-REJECTION 0d66d55dd153a589 -->

At `r=9` (702 columns), the smooth block has mod-43 rank 690 and the complete
generic E7 layer has rank 12 on its 12-dimensional kernel. Thus this third
bounded ambient has an exact full-rank smooth-plus-generic certificate, still
without claiming a q8 pencil or the missing node/overlap implementation.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA9-SMOOTH-GENERIC-REJECTION 8035f0b465038995 -->

With SageMath 10.9 and the default commands above, the current generated
artifacts have SHA-256
`1465de949187ad22bb8fca47daad7f66f3023dfb498d27ea8130225f5ae1ca8f`
(E8 resolution),
`97dbadd5a9e00b95106266aab04e3d37911b9800f8f219b9028a34b32530fb85`
(E8 marked-chord module),
`91259edd627dea02445c272416ad4a76fae1e20f133f25c0254d65b3438b47d4`
(formal E7 resolution),
`0a408e706820a62a4c2e3290e615f5883808803565912b5b3e3e0867761f5f58`
(formal E7 valuation atlas),
`a0699e4ec75930cc93a9706ddf96f4ffc744954809e53a24029bd8c6668843f7`
(third marked-divisor E7 lattice target),
`598480530e0c845fcecb1bdcdef4756dec61f668f0aed07c3c71ccc1049b6878`
(third marked-divisor symbolic degree-44 generic ambient),
`731a5b57cd0a11dd72da85d1da8dfce19fedfdc692014a6dcfd3264bd78ff5a4`
(third marked-divisor formal E7 Cartier chart atlas),
`a898cdd4147ef915c45b2eaeb63ff1542a7c89076d368be8a4736532b2b8619a`
(actual-H92-versus-formal-E7 rejection diagnostic),
`14378f4718d3fbe781d5b351ba4943a962a430d9827c0ce285cd1125a9e8c500`
(full actual-H92 E7 resolution tree),
`9b2c3e5356065babe6f3fe336f8ac8400b9e0225e8f226a4e59885c519fd03eb`
(actual-H92 third E7 Cartier chart factors),
`a8f19a1e205ed2250c83c312a4bf722462867d289ddd598987dd5112f1c33177`
(actual-H92 E7 edge-chart pullbacks),
`ae7eb1e79a2fb41ab05d0092cb8c04663307bd8aa1e0cb562a8d3a014e94f451`
(actual-H92 E7 old-coordinate valuation atlas),
`3796ee20121a94ce6d3a707c0cd119983b64fce79336fa99a5a894729174900c`
(first q=6 actual-H92 E7 quotient block),
`928955731227eead4006392cbe13fc5b95a1e2e951acbc2834eea6bdc1bbaacc`
(first q=6 complete actual resolved RR cover),
`7848c3a506b2255fc1e42cab9ced9b72f8852e26aa703f9b23e2b2417474d2ed`
(third marked divisor's actual-H92 E7 quotient block),
`bbfdd659a38dbdbe6f0fd272970df251049ceb4f9344068de7f2f4abe9a6ed9f`
(third marked chord evaluated in that actual-H92 quotient),
`c1cb5457c781af53d08522144b15d95ec2452a90ec0c6ee65377c3f51f835a10`
(third degree-44 generic ambient evaluated in that quotient),
`744d194f7a2799bed4e65aa0369e3a36ad99a9374518531a5aff1f8b9adbcc5b`
(first H3 q=6 actual neighbour-hop certificate),
`e27f10c238d0700e3fa236dbe616aeec2575e4f17dd19bd2f934404ece6188d2`
(rejection of the insufficient old-monomial E7 complete ideal),
`a73ccb1c729814219f172df4c6feb49c05859125db1cff7591eeb8544fb664e1`
(P1 and -P1 actual-H92 E7 blow-up trace),
`7e3e8ad89b31a76ccaaee078f7f2830b5d5f09f875ceda734b302a760a9982a8`
(P1 actual-H92 E7 marked smooth-point module),
`d2d2af0e172d5d3efc06c79e4e1d30b3c325921cc5c2ab1de0d0efd071cf4b49`
(third marked-point exact E7 Laurent-series evaluation),
`e860672b80d9fc933f04c4d7b825b186702226986cd17160cc788996ce44a022`
(third marked-chord actual E7 unit certificate),
`a950defe18be876d96b215d287d02aebfe4c0a375b788b1c3bfbf0ff864b839d`
(smooth `P1.O` module),
`fd4738a5a2a3ae8ed2971c9b48f2bc38702c3930391d26b107a6413f0749c5aa`
(q=8 smooth collision algebra frame),
`88466369aea8f838dde85051b174f6f2ea4c6df3edb29828b70ee8e79b9975cb`
(global Riemann--Roch kernel),
`9217d0492046954d8522680f9ddf2a0f5c0177921db761b584c281527b41145d`
(pencil elimination),
`5eb43d9a0d04195e7a6e38ebd337b0e10a3b1a2eb9246a3b02cce4331bcd36ac`
(child Jacobian), and
`3d8d79c24bccf74ce0f2878ac4e02c3edbaa05bcd1537bf9c99608b5b95f8221`
(transported child zero section), and
`156821384b45fd5e731dce130686b549030d64d450ece78bfb9f9083bbaf3005`
(the two transported E7-infinity sections), and
`335a9cb6c1060ac170c063f99bb02d4c4357fa2426d37b4dc3efd447ac2b62ad`
(component-section lattice identities), and
`c4b7e38f0ea9fc3f748200ca9923ea3ffe5c0028c979e5f81be6507954d7c822`
(Weyl-aware transported divisor classes and their rank-three Gram), and
`e02e2803387d3a7f53907f548b275bb592d366f653f630f6ba8c9ef2611f3e37`
(the exact 1024-adic P2 reconstruction), and
`6cf559999a256cc8c323786d833c099ac38b3238e468ae885c981b10fdae9d9b`
(the missing-vertical-correction regression), and
`baefad65dbda1043005631a42af3bcd9e94e674603a0377b0d74132c683d09f4`
(H3 preflight).

The first equation-level H3 q=6 hop is complete.  The corrected actual E7
trivialization is replayed in
[`scripts/certify_h92_q6_actual_resolved_rr_cover.sage`](scripts/certify_h92_q6_actual_resolved_rr_cover.sage): its complete resolved cover has the
eight-by-ten collision matrix, exact rank eight, and a two-dimensional
kernel.  The neighbour-hop certificate then records the E8+E6 child,
rank-three Mordell--Weil target, predicted transported Gram matrix, and
torsion-free Shioda--Tate discriminant `948` (root determinant `3`, height
determinant `316`).  It also replays the exact Néron--Severi old-E7 pairing
rows for the third correction and the degree balance `4812-4811=1`.

For q=8, the next unresolved operation is deliberately narrower.  The
actual E7_4--E7_3 principal frame is cleared in
[`scripts/derive_h92_q8_e7_4_3_principal_node_clearing.sage`](scripts/derive_h92_q8_e7_4_3_principal_node_clearing.sage): each bounded ambient
combination is regular at that node precisely when its cleared numerator is
in `(t^17)`.  In the completed regular parameters this is product
divisibility by `Z^51*Y^68`.  The accompanying
[`scripts/certify_h92_q8_e7_4_3_node_divisibility_geometry.sage`](scripts/certify_h92_q8_e7_4_3_node_divisibility_geometry.sage)
certifies that `R/(t^17)` has Krull dimension one, whereas the tempting
finite rectangle `R/(Z^51,Y^68)` has length 3468 and is a different quotient.
It must therefore not be used as a node condition.  The modular local-normal-
form probe
[`scripts/probe_h92_q8_e7_4_3_principal_node_local_normal_form_modp.sage`](scripts/probe_h92_q8_e7_4_3_principal_node_local_normal_form_modp.sage)
uses Singular's local standard basis for the actual ideal `(surface,t^17)`.
It records the finite image of a chosen ambient in that infinite quotient;
on the 54-column seed modulo `43` this image has rank 54 in 22,126 normal-
monomial coordinates.  Thus it is a chart-faithful modular regression, not a
finite quotient substitution or a characteristic-zero node block.  The q=8
implementation still lacks the remaining all-chart, overlap-compatible
residual maps.
In particular, it does not claim q=8 `h0=2`, a pencil, a child equation, or
transported child sections.

For that pinned seed, the local normal-form map has full column rank modulo
`43`, while the source denominators and all five common-clearing factors are
units at the chart origin after reduction.  The companion
[`scripts/certify_h92_q8_e7_4_3_principal_node_good_reduction.sage`](scripts/certify_h92_q8_e7_4_3_principal_node_good_reduction.sage)
therefore proves injectivity of this one actual principal-node condition over
QQ by primitive-vector reduction.  This is a bounded-ambient rejection, not
a characteristic-zero coordinate presentation or a q8 pencil.

The same unit-clearing construction is now derived on all six actual E7 edge
charts in
[`scripts/derive_h92_q8_e7_node_principal_clearings.sage`](scripts/derive_h92_q8_e7_node_principal_clearings.sage).
It includes the cancellation-sensitive `E7_2--E7_5` chart and records each
local chord numerator and its principal `(t^17)` condition.  This is an exact
input atlas for the local-condition evaluator, not yet six quotient matrices
or their common kernel.

[`scripts/probe_h92_q8_e7_node_principal_local_normal_form_modp.sage`](scripts/probe_h92_q8_e7_node_principal_local_normal_form_modp.sage)
now consumes that atlas one chart at a time.  It agrees with the specialized
E7_4--E7_3 test (rank 54 modulo 43) and has full rank on the same 54-column
seed at `E7_1--E7_4`, `E7_2--E7_5`, `E7_3--E7_7`, and `E7_7--E7_2` as well.
For the expensive `E7_3--E7_6` local-degree standard basis, the compiler also
has a rigorously one-way Artinian-corner obstruction: the actual local ideal
maps to `(surface,Z^34,U^34)`, whose finite image likewise has rank 54.  It
is used only in the forward direction—an actual local solution vanishes in
that quotient—and never as a replacement for the principal node condition.
The associated good-reduction checker proves the resulting fixed-ambient
characteristic-zero injectivity statement, while retaining that one-way
boundary.
[`scripts/certify_h92_q8_all_e7_node_modular_regression.sage`](scripts/certify_h92_q8_all_e7_node_modular_regression.sage)
aggregates the five actual local images and this one-way sixth-chart witness.
These are modular finite-ambient obstructions; they are not an all-chart
characteristic-zero matrix or a completed Čech calculation.

The first actual overlap maps are now available in
[`scripts/derive_h92_q8_e7_sibling_chart_transitions.sage`](scripts/derive_h92_q8_e7_sibling_chart_transitions.sage).
For the two sibling U/Z chart pairs, it derives
`(Z,U,Y) -> (U*Z,1/Z,Y/Z)` from the blow-up coordinates, verifies the old
H92 `(t,x,y)` maps exactly, and records the q8 Cartier-factor ratios
`Y^4/Z^2` and `1/Y`.  These are the required transition functions for future
overlap evaluators; they are not assertions that the ratios are units, nor a
complete E7 Čech cover.
<!-- status-consumer: EC-K3-H3-Q6 177cd6e614c8b8e0 -->
