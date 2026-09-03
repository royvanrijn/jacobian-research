# Reproducing the results

<!-- status-consumer: EC-K3-H3-Q12O5867-POINT-FACTORY 9399c93ee42ee2a4 -->

The Makefile is the public verification interface.  Run commands from the
repository root after creating the Python environment described in the main
[README](README.md).

## Separate elliptic-curve programme

The elliptic-curve rank/conductor work has its own
[reproduction catalogue](elliptic-curves/REPRODUCE.md).  Its combined
dependency-light gate is:

```bash
make verify-elliptic-curves PYTHON=python3
```

The active Elkies rank-32 programme is the compact-`t` positive-control and
residual 2-Selmer section of that catalogue. It requires an actual completed
global/local descent before any expensive point search; a score or timeout is
never an authorization. The direct q12 point-search entry points enforce this
same-parameter, same-minimal-model gate; the obsolete pre-descent x-ansatz
search is parked. The rank-28 bad-place ledger and factor-supplied PARI backend
are reproducible exact descent inputs, but remain fail-closed until a complete
Selmer dimension is returned. The same catalogue now includes a stage-aware
factor-supplied `S`-class worker and an exact BNF-free rank-28 pilot. The first
stops in class-group relation generation before certification; the second is
far below a valid factor-base generation bound. An exact depressed-cubic
variant lowers the defining-order index by 27 but reaches the same PARI
relation plateau. None is a Selmer bound or permission to search. The exact
rank-28 local calibration additionally shows
that eleven certified global quotient directions add zero bad-place signature
rank. Four odd places and infinity have certified full known-point coverage;
a bounded resumable norm-one-cover pilot records 60 selected local witnesses
and 24 inconclusive place tests, but no everywhere-local or Selmer class. The
eleven public complement directions now provide genuine cover controls with
exact rational witnesses and certify residual Selmer dimension at least 11;
the complete upper bound and threshold 15 remain open.

The exact q12/orbit5867 arbitrary-point map and its backward calibration on
all 42 public-complement points are replayed by:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/q12o5867_genus_one_point_factory.sage \
  --mode controls \
  --output artifacts/generated-results/elkies-k3-q12o5867-genus-one-point-factory-controls.json
```

This certifies the birational round trips and records numerical parent heights;
it does not yet provide coordinates in a specialized equation-level MW13 basis.

## Level-474 branch in the Kumar H92 chart

This exact Sage replay downloads the pinned Humbert-21 equation and extracts
the pinned Elkies--Kumar H92 chart, then identifies the correct component
modulo 11.  It normally completes in under one minute.

```bash
mkdir -p artifacts/local/humbert-inputs
curl -L https://www.maths.usyd.edu.au/u/davidg/ThesisData/SatakeHumbert/level1_21.txt \
  -o artifacts/local/humbert-inputs/level1_21.txt
curl -L https://export.arxiv.org/e-print/1209.3527 \
  -o artifacts/local/humbert-inputs/elkies-kumar-source.tar
tar -xf artifacts/local/humbert-inputs/elkies-kumar-source.tar \
  -C artifacts/local/humbert-inputs 21/21.txt 92/92.txt 92/igusa92.txt
sage -python elkies-k3/scripts/verify_h21_h92_level474_branch.sage \
  --h21 artifacts/local/humbert-inputs/level1_21.txt \
  --h92 artifacts/local/humbert-inputs/92/igusa92.txt \
  --output artifacts/generated-results/elkies-k3-h21-h92-level474-branch-mod11.json
```

Expected terminal status:

```text
H21H92|stage=complete|status=PASS_COMPUTATIONAL_BRANCH_IDENTIFICATION
```

The result is an exact one-prime computational identification.  The separate
characteristic-zero reconstruction and normalization below promote this to an
exact birational identification; see
[`elkies-k3/KUMAR_E7E8_BACKTRACK.md`](elkies-k3/KUMAR_E7E8_BACKTRACK.md).

The corresponding characteristic-zero component equation is reconstructed
and proved with the following modular workers.  Each worker normally finishes
in under one minute; the jobs are independent and two may be run concurrently
on a typical workstation.

```bash
mkdir -p artifacts/local/humbert-level474-modp
for p in 17 19 23 29 31 41 43 47; do
  sage -python elkies-k3/scripts/factor_h21_h92_level474_modp.sage \
    --prime "$p" \
    --h21 artifacts/local/humbert-inputs/level1_21.txt \
    --h92 artifacts/local/humbert-inputs/92/igusa92.txt \
    --output "artifacts/local/humbert-level474-modp/target-p${p}.json"
done
```

Combine those images and run the exact 961-specialization degree-bound
certificate:

```bash
sage -python elkies-k3/scripts/reconstruct_h21_h92_level474_qq.sage \
  --h21 artifacts/local/humbert-inputs/level1_21.txt \
  --h92 artifacts/local/humbert-inputs/92/igusa92.txt \
  --modular-factor artifacts/local/humbert-level474-modp/target-p17.json \
  --modular-factor artifacts/local/humbert-level474-modp/target-p19.json \
  --modular-factor artifacts/local/humbert-level474-modp/target-p23.json \
  --modular-factor artifacts/local/humbert-level474-modp/target-p29.json \
  --modular-factor artifacts/local/humbert-level474-modp/target-p31.json \
  --modular-factor artifacts/local/humbert-level474-modp/target-p41.json \
  --modular-factor artifacts/local/humbert-level474-modp/target-p43.json \
  --modular-factor artifacts/local/humbert-level474-modp/target-p47.json \
  --output artifacts/generated-results/elkies-k3-h21-h92-level474-factor-qq.json
```

Expected terminal status:

```text
H21H92QQ|stage=complete|status=PASS_CHARACTERISTIC_ZERO_FACTOR
```

This proves the characteristic-zero degree-21 component equation.  Normalize
it and verify the exact map to the published level-474 sextic with:

```bash
sage -python elkies-k3/scripts/normalize_h21_h92_level474_qq.sage \
  --factor artifacts/generated-results/elkies-k3-h21-h92-level474-factor-qq.json \
  --adjoint-cache artifacts/local/humbert-level474-adjoint.json \
  --output artifacts/generated-results/elkies-k3-h21-h92-level474-normalization.json
```

The first run computes and caches a degree-nine adjoint basis and normally
takes under one minute.  Cached runs take a few seconds.  The script takes the
exact involution quotient, applies a pinned Cremona reduction from degree 13
to degree 11, saturates and LLL-reduces the adjoint lattice, recovers a
degree-one parameter from its order-eight/order-nine osculating kernels, and
checks the Padé inverse by exact substitution.  It then verifies an exact
Möbius-plus-square identity with the published sextic.

Expected terminal status:

```text
H21H92NORM|stage=complete|status=PASS_LEVEL474_NORMALIZATION
```

The published non-CM point is also an exact rational source anchor for the
H3 construction.  The following checker maps `(13/7,12048/343)` to the H92
chart, recovers a rational H21 presentation, proves that the short
Weierstrass models are isomorphic over `QQ` (the twist parameter itself is a
rational square), and replays the bidegree-`(3,3)` H21 entrance cubic from
the pinned `21/21.txt` source.  The displayed rational point is a nonflex.
It also records that the H21 and H92 oriented Hilbert-cover coordinates both
have square class `-52203427`, so their common orientation field is quadratic
rather than biquadratic.

```bash
sage -python elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage \
  --h92 artifacts/local/humbert-inputs/92/igusa92.txt \
  --h92-entrance artifacts/local/humbert-inputs/92/92.txt \
  --h21-entrance artifacts/local/humbert-inputs/21/21.txt \
  --output artifacts/generated-results/elkies-k3-h3-noncm-q6-source-anchor.json
```

Expected terminal status:

```text
H3NONCMQ6|stage=complete|status=PASS_H3_NONCM_Q6_SOURCE_ANCHOR
```

The authoritative artifact SHA-256 is
`0560b1921c87ad2d8db6c293ce070cb30aa75626315801e3e4a71cad59573ea5`.
This certificate stops at the `A2+A6+E8 -> E7+E8` H21 entrance.  It proves
the unmarked H21/H92 surfaces are isomorphic over `QQ`, but by itself does not
identify the entrance cubic with the desired `E7+E8 -> E8+E6` q=6 pencil.
The separate exact descent checker below now proves that the height-`21/2`
section and both signed q=6 divisor classes are individually defined over
`QQ(r,s)`.  Constructing an explicit basis of `H0(O+(-P1)-F)` is the remaining
equation-level step.  The common nonsquare orientation field is compatibility
data for the H21/H92 intersection, not a K3-section descent obstruction.

```bash
sage -python elkies-k3/scripts/verify_h21_q6_section_descent.sage
```

The pinned ancillary model starts with split `A2+A6+E8` fibers and a
3-neighbor over `QQ(r,s)`.  Its trivial lattice has squarefree determinant
`3*7=21`, so these rational curves already generate the generic H21
Neron--Severi lattice integrally.  The transported `E7+E8` MW generator is
therefore individually rational and has height `21/2`; `P1`, `-P1`, and
`D=O+(-P1)-F` are all Galois fixed.  The generated artifact
`artifacts/generated-results/elkies-k3-h21-q6-section-descent.json` has
SHA-256
`9ccdfc7b7a1ca79d549c161e9922051e9f90d7b89ddc81057e17188eedc2a4d2`.

The complementary H92 direction in the rank-two H3 source is certified by:

```bash
sage -python elkies-k3/scripts/verify_h92_section_descent.sage
```

The pinned H92 ancillary construction begins with split `D6+A8+A1` fibers
and an explicit section of height `4-1/2-20/9=23/18`.  Adjoining that section
gives determinant `92` and Smith factors `(2,46)`.  The checker exhausts all
92 discriminant classes and finds no nonzero isotropic class, proving that
there is no proper even overlattice.  It also replays the split-square and
section identities at the exact non-CM target.  Since the ancillary two- and
three-neighbor parameters are rational over `QQ(r,s)`, the height-`46`
generator on the final `E7+E8` fibration is individually rational.  Thus both
directions in the H3 height Gram `[[21/2,3],[3,46]]` are rational, rather than
only their span.  Expected terminal status:

```text
H92DESCENT|status=PASS_H92_SECTION_Q_DEFINED
```

The generated artifact
`artifacts/generated-results/elkies-k3-h92-section-descent.json` has SHA-256
`fe525f75fa87c31afb34755fe63fc778349d2843010eb5c9b17ce6d8b8712e40`.
This is an exact lattice and field-of-definition certificate; it does not
recover explicit height-`46` section coordinates on the final short H92
model.  Those coordinates are not needed for the first q=6 chord, which uses
the height-`21/2` section below, but may be needed for later equation-level
section transport.

Before attempting coordinate reconstruction of that second direction, replay
the exact marked-basis pole reduction:

```bash
sage elkies-k3/scripts/analyze_h3_p2_pole_profile.sage
```

It proves that `P2+m*P1` has its unique minimal zero-section intersection at
`m=0`: `P2.O=21`, while the nearest translate `(P2-P1).O=24`.  Hence the
pole-reduced target is `P2` itself; no small lattice basis change lowers the
denominator degree.

The exact marked intermediate section used in the remaining degree-29 divisor
transport is reconstructed from modular records by:

```bash
H92P2_PRIME=2305843009213693967 \
H92P2_MODULAR_OUTPUT=artifacts/generated-results/h92-p2-modular/intermediate-2305843009213693967.json \
sage elkies-k3/scripts/sample_h92_p2_intermediate_modp.sage
sage elkies-k3/scripts/crt_h92_p2_intermediate.sage
H92P2_PROBE=0 sage elkies-k3/scripts/probe_h92_p2_final_divisor.sage

# Check that the residual-origin hyperplane class is not a small P1 multiple.
sage elkies-k3/scripts/diagnose_h92_p2_normalization_modp.sage

# Extract the pole-reduced modular half from the doubled canonical class.
H92P2_CANDIDATE_INPUT=artifacts/generated-results/h92-p2-candidate-mod-100003-500.json \
H92P2_DOUBLE_X_DEGREE=184 \
H92P2_HALF_OUTPUT=artifacts/generated-results/h92-p2-half-mod-100003-v2.json \
sage elkies-k3/scripts/extract_h92_p2_half_modp.sage

# Hensel lift that fixed-degree half.  The first two outputs are reusable
# p-adic checkpoints; only the final record is retained as a generated result.
sage elkies-k3/scripts/lift_h92_p2_hensel.sage --precision 128 \
  --output /tmp/h92-p2-hensel-p128.json
sage elkies-k3/scripts/lift_h92_p2_hensel.sage --precision 512 \
  --seed /tmp/h92-p2-hensel-p128.json --output /tmp/h92-p2-hensel-p512.json
sage elkies-k3/scripts/lift_h92_p2_hensel.sage --precision 1024 \
  --seed /tmp/h92-p2-hensel-p512.json \
  --output artifacts/generated-results/elkies-k3-h92-p2-hensel.json
sage elkies-k3/scripts/verify_h92_p2_hensel_lift.sage
```

The worker constrains each record to the marked degree-three divisor and the
fixed coordinate degrees `(22,18)` and `(33,27)`.  The final probe verifies
the exact intermediate Weierstrass identity and constructs the degree-29
divisor.  Its canonical class `3D-29H` is the doubled `P2`; the fixed-profile
Hensel lift and final verifier recover the exact height-46 H92 section.

The complete defining equation of the H3 source family is exported by:

```bash
sage -python elkies-k3/scripts/export_h3_level474_source_family.sage
```

Over the published genus-two curve

```text
Y^2=-27*X^6+198*X^4-171*X^2+576,
```

the normalization artifact gives rational functions `t(x)`, `a(x)` and a
published-`Y` multiplier.  The checker inverts the linear-fractional
`X`-map, sets `Y0=Y/m(x)`, and proves the exact H92 chart formulas

```text
r=(a+Y0)/2,   s=2/(Y0-a).
```

It then proves the degree-21 H21/H92 component equation after composition and
exports the short elliptic K3 family

```text
v^2=u^3+(A1*tau^3+A*tau^4)*u+(B1*tau^5+B*tau^6+B2*tau^7),
```

where the five coefficient formulas are the pinned H92 functions evaluated
at this `(r,s)`.  The fibers are `E7` at `tau=0` and `E8` at infinity; the
two individually rational MW directions have height Gram
`[[21/2,3],[3,46]]`.  Exact specialization of `(X,Y)=(13/7,12048/343)` gives
`(r,s)=(-3621005/690947,158286/143585)` and the same five coefficients as the
source-anchor certificate.  Expected status:

```text
H3SOURCE|status=PASS_EXACT_H3_SOURCE_FAMILY
```

The generated artifact
`artifacts/generated-results/elkies-k3-h3-level474-source-family.json` has
SHA-256
`8f5afd11e1d8979d57cb1a569833309f9664c19cd47194af0581a5cbbf8f1d59`.
This certifies the genus-two H3 `E7+E8/MW2` source family, not the downstream
rootless MW17 equation or its specialization to curve 273.
<!-- status-consumer: EC-K3-H3-SOURCE a4bb40c9c9d0ff09 -->

A global exact determination of the rational points uses Magma's two-cover
descent and elliptic Chabauty implementation:

```bash
magma elkies-k3/scripts/prove_h3_level474_rational_points.m \
  | tee artifacts/generated-results/elkies-k3-h3-level474-rational-points.txt
```

The certificate uses the Q-isomorphic model

```text
y^2=-3*x^6+22*x^4-19*x^2+64.
```

Its two-cover descent has three locally soluble classes, whose quartic-factor
elliptic quotients reduce to two covers over a cubic field.  At the good prime
`41`, elliptic Chabauty proves the complete rational `x`-image sets
`{0}` and `{-13/7,-1,1,13/7}`; exact substitution then gives all affine
fibres, and the nonsquare leading coefficient excludes rational infinity.
This requires Magma `2.29-9` or a compatible later release, and uses no GRH or
BSD assumption.  Expected terminal status:

```text
H3GLOBAL|status=PASS_GLOBAL_H3_LEVEL474_RATIONAL_POINTS
```

The pinned output is
[`artifacts/generated-results/elkies-k3-h3-level474-rational-points.txt`](artifacts/generated-results/elkies-k3-h3-level474-rational-points.txt).
The earlier quotient-sieve command remains available as a bounded independent
cross-check:

```bash
sage -python elkies-k3/scripts/sieve_h3_level474_rational_points.sage
```
<!-- status-consumer: EC-K3-H3-PTS 8f0a27c947843b4a -->

The marked section itself is now recovered exactly on the smaller-coefficient
H92 short model.  First generate nine disjoint modular windows.  Exceptional
primes at which interpolation degenerates are skipped deterministically; the
commands below retain 204 good primes in total.

```bash
sage -python elkies-k3/scripts/recover_h21_p1_from_entrance_cubic.sage \
  --prime-start 100 --prime-count 12 --target h92 \
  --output artifacts/generated-results/elkies-k3-h92-p1-mod-window-100.json
for START in 199 383 557 733 887 1069 1300 1600; do
  sage -python elkies-k3/scripts/recover_h21_p1_from_entrance_cubic.sage \
    --prime-start "$START" --prime-count 24 --target h92 \
    --output "artifacts/generated-results/elkies-k3-h92-p1-mod-window-${START}.json"
done
```

CRT-lift the windows and perform the characteristic-zero check with:

```bash
sage -python elkies-k3/scripts/lift_h21_p1_modular.sage \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-100.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-199.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-383.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-557.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-733.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-887.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-1069.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-1300.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-1600.json \
  --output artifacts/generated-results/elkies-k3-h92-p1-lift.json
```

The structured simultaneous-LLL reconstruction uses the exact pole identity
`D(u)=d4*u^4*Z4(u)^2`, with `Z4(0)=1`.  The 204-prime CRT modulus has 1,945
bits.  The resulting H92 coordinates have degrees `(10,12)` for `x` and
`(15,18)` for `y`; direct substitution proves their Weierstrass identity over
`QQ`, and an exact marked-fiber incidence selects the square-root sign.  The
expected terminal status is

```text
H21P1LIFT|primes=204|modulus_bits=1945|exact_square=1|status=PASS_EXACT_H92_P1
```

The pinned artifact
`artifacts/generated-results/elkies-k3-h92-p1-lift.json` has SHA-256
`c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397`.
It supersedes historical SHA-256
`0602c3b199629c6f460c9b7c728e048822418ecf85bf54807852be3d97b66616`,
which lacked only the exact orientation-incidence block and used the older
H21 status label.  This certifies the exact height-`21/2` section on the
rational H92 model. The associated q=6 lattice divisor and chord ambient now
have an all-edge actual E7 cover; the certified ten-dimensional matrix has a
two-dimensional kernel and yields the first exact fibration hop below.

### Resolved-chart q=6 compiler preflight

The strict equation-level compiler is deliberately separate from the
lattice-only neighbour engine.  It accepts vertical conditions only as maps
from an explicit Riemann--Roch ambient basis to finite quotients of actual
resolved blow-up charts; a Kodaira label or Smith saturation cannot stand in
for such a map.  Run its core regression and the first H3 gate with:

```bash
sage -python elkies-k3/scripts/verify_elliptic_neighbor_compiler.sage
sage -python elkies-k3/scripts/compile_h3_first_q6_preflight.sage \
  --output artifacts/generated-results/elkies-k3-h3-q6-compiler-preflight.json
```

The H3 search's `q=6` label is not its old-fibre degree: after the recorded
reflections, `D=O+(-P1)-F` has `D.F=2`.  The preflight therefore records a
degree-two generic-fibre basis `(1,m)` with marked chord
`(y-y(P1))/(x-x(P1))`, rather than treating the five monomials of `L(5O)` as
the relevant space.  It verifies the exact P1 input and the E7 resolved
module, replays all 22 Weyl reflections, and checks the E7/E8 affine and
simple fibre-wall pairings. The exact eight-step E8 chart tree itself is
reproduced by:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_e8_resolution.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-e8-resolution.json
sage -python elkies-k3/scripts/derive_h92_q6_e8_p1_branch_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-e8-p1-branch-module.json
sage -python elkies-k3/scripts/derive_h92_q6_smooth_po_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-smooth-po-module.json
```

The middle command derives the complete E8 local module, including the affine
II* component; the last certifies the four smooth P.O local changes of basis.
The preflight is an exact local-input gate, but its former marked E7
trivialization has been rejected by the following actual-chart audit:

```bash
sage -python elkies-k3/scripts/audit_h92_q6_actual_e7_marked_chord_order.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-e7-marked-chord-order-audit.json
```

It proves `ord_Z(m/t)=-1` at the generic point of actual `E7_5`; the
corrected all-edge cover below supplies the missing line-bundle
trivialization rather than treating this audit as a failure of the q6 hop.

The corrected marked-chart calculation is reproducible with:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_p1_actual_e7_marked_module_corrected.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module-corrected.json
```

It proves `Z*m/t=unit/W` at `-P1`, so `<1,m>` is locally valid on that chart
because `t/Z` is a unit. The all-edge continuation is:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_all_edge_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-e7-all-edge-module.json
```

It proves `m` has exceptional orders at least `(1,1,3,2,0,2,2)` and no
horizontal pole except the marked `-P1`, so it completes the actual E7 module
cover without adding a q6 matrix row.

The exact global assembly and child-Jacobian gates are:

```bash
sage -python elkies-k3/scripts/assemble_h92_q6_global_rr.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-global-rr.json
sage -python elkies-k3/scripts/eliminate_h92_q6_global_pencil.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-global-pencil-elimination.json
sage -python elkies-k3/scripts/audit_h92_q6_pencil_marked_section_degrees.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-marked-section-degrees.json
sage -python elkies-k3/scripts/certify_h92_q6_child_jacobian.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json
sage -python elkies-k3/scripts/derive_h92_q6_e7_resolution.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-e7-resolution.json
sage -python elkies-k3/scripts/derive_h92_q6_e7_valuation_atlas.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-e7-valuation-atlas.json
sage -python elkies-k3/scripts/derive_h92_q6_third_e7_local_target.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-e7-local-target.json
sage -python elkies-k3/scripts/assemble_h92_q6_third_generic_rr_ambient.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-generic-rr-ambient.json
sage -python elkies-k3/scripts/derive_h92_q6_third_e7_cartier_charts.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-e7-cartier-charts.json
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_resolution.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-diagnostic.json
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_resolution_full.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json
sage -python elkies-k3/scripts/derive_h92_q6_third_actual_e7_cartier_charts.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-actual-e7-cartier-charts.json
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_chart_pullbacks.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_valuation_atlas.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-e7-valuation-atlas.json
sage -python elkies-k3/scripts/trace_h92_q6_p1_actual_e7.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-trace.json
sage -python elkies-k3/scripts/derive_h92_q6_p1_actual_e7_marked_module_corrected.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module-corrected.json
sage -python elkies-k3/scripts/derive_h92_q6_actual_e7_all_edge_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-e7-all-edge-module.json
sage -python elkies-k3/scripts/certify_h92_q6_actual_resolved_rr_cover.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-actual-resolved-rr-cover.json
sage -python elkies-k3/scripts/reject_h92_q6_third_old_monomial_e7_ideal.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-rejected-old-monomial-e7-ideal.json
sage -python elkies-k3/scripts/derive_h92_q6_third_actual_e7_quotient_block.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-actual-e7-quotient-block.json
sage -python elkies-k3/scripts/evaluate_h92_q6_third_marked_chord_actual_e7_quotient.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-marked-chord-actual-e7-quotient.json
sage -python elkies-k3/scripts/evaluate_h92_q6_third_generic_ambient_actual_e7_quotient.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-generic-ambient-actual-e7-quotient.json
sage -python elkies-k3/scripts/evaluate_h92_q6_third_e7_point_series.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-e7-point-series.json
sage -python elkies-k3/scripts/certify_h92_q6_third_e7_chord_units.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-third-e7-chord-units.json
sage -python elkies-k3/scripts/derive_h92_q6_child_zero_section.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-zero-section.json \
  --component-output artifacts/generated-results/elkies-k3-h92-q6-child-e7-infinity-sections.json
sage -python elkies-k3/scripts/certify_h3_q6_component_section_lattice.sage \
  --output artifacts/generated-results/elkies-k3-h3-q6-component-sections.json
sage -python elkies-k3/scripts/certify_h3_q6_weyl_section_transport.sage \
  --output artifacts/generated-results/elkies-k3-h3-q6-weyl-section-transport.json
sage -python elkies-k3/scripts/certify_h3_q6_actual_neighbor_hop.sage \
  --output artifacts/generated-results/elkies-k3-h3-q6-actual-neighbor-hop.json
sage -python elkies-k3/scripts/lift_h92_p2_hensel.sage \
  --precision 1024 \
  --output artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json
sage -python elkies-k3/scripts/verify_h92_p2_coordinates.sage \
  --input artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json
sage -python elkies-k3/scripts/transport_h92_q6_third_modp.sage
```

These commands certify the 10-by-8 vertical matrix, `h0(D)=2`, a degree-four
genus-one model, and the Jacobian fibre signature `II*+IV*+6I1`, hence
`E8+E6` and MW rank three in the rank-19 source Néron--Severi lattice. The
final first-hop certificate replays the rank-three Néron--Severi height Gram,
although minimized equation-level coordinates for every rank-three section
remain a separate transport task. The zero-section command has `D.O=1`,
`u=T` on the new base, and exact coordinates on the minimized child Jacobian.
The component
certificate pins the two old E7 components which form the
`[[8/3,1/3],[1/3,8/3]]` leading block of the target height matrix. The optional
component output transports both binary-quartic infinity points exactly. The
preceding E7 chart replay assigns them without a fibre-type guess: affine E7
is `plus`, and E7_7 is `minus`. The final Weyl transport identifies the
remaining marked child direction as `4*(-P1)`, `(-P1)`, and `22*(-P1)-P2`;
the same certificate replays the full predicted rank-three
Gram. These projections are not old-model section coordinates. The p-adic
lift reconstructs all 139 normalized P2 coefficients and independently
verifies the exact H92 equation; only divisor-aware conversion of those
classes to minimized child coordinates remains.
The final command is a negative guard: after calibrating the reconstructed
coordinate as `-P2` in the H3 frame, naive generic-fibre evaluation of the
third MW projection has q=6-pencil degree 4769 modulo 100003. The Weyl
transport certificate supplies the exact E7-plus-fibre correction which turns
the corresponding degree-4812 horizontal divisor into the child section.
The valuation-atlas command records the `(Z,U,Y)` orders in the *formal* E7
normal-form resolution. The third-E7 target records the matching formal
exceptional cycle `22*c_q6`; neither is a transported H92 chart.
The following symbolic ambient has exact degree 44 and 44 basis elements:
43 old Weierstrass monomials plus the chord at the exact group-law expression
`22*(-P1)-P2`. It is evaluated inside each resolved chart, avoiding an
irrelevant globally reduced rational-function coordinate.

The child-Jacobian default pin is the SageMath 10.9 replay hash
`5eb43d9a0d04195e7a6e38ebd337b0e10a3b1a2eb9246a3b02cce4331bcd36ac`.
It is generated through the reusable exact finite-place minimization and
Kodaira aggregation route, bound directly to the same certified q6 resolved
RR kernel that produces the pencil; the q6 invariants and transported Gram
are unchanged.
The corresponding regenerated final-hop artifact has SHA-256
`744d194f7a2799bed4e65aa0369e3a36ad99a9374518531a5aff1f8b9adbcc5b`.
Its certified Shioda--Tate data are root determinant `3`, height determinant
`316`, trivial torsion/glue index, and absolute Néron--Severi discriminant
`948`.  It additionally pins the three source curves meeting the new fibre
once, the exact Néron--Severi old-E7 pairing rows of the third correction, and its
horizontal-plus-vertical degree balance `4812-4811=1`. The first two
low-height sections are also pinned as exact child-Jacobian points by the
actual affine-E7 blow-up chart and its two binary-quartic infinity signs.
The Cartier-chart atlas writes the six formal E7 node-chart factors. The
following actual-H92 diagnostic proves that its second standard node is smooth
in the H92 germ. The full actual-H92 tree then resolves the displaced
second-U node at `Z=-A1/B1` plus three third-stage nodes; these exact charts,
rather than the formal atlas, are the input for a finite quotient matrix. The
following actual Cartier atlas supplies the six transported integral-vertical
factors. The following pullback atlas then maps the old H92 coordinates into
each edge chart, before the remaining marked-chord jet conditions are imposed.
<!-- status-consumer: EC-K3-H3-Q6 177cd6e614c8b8e0 -->

In the blocked state the preflight makes no
claim about vertical codimension, `h0(D)`, a child equation, or transported
sections.  See [`elkies-k3/ELLIPTIC_NEIGHBOR_COMPILER.md`](elkies-k3/ELLIPTIC_NEIGHBOR_COMPILER.md).
The commands use SageMath 10.9; the corresponding E8-resolution and H3-
preflight artifact hashes are recorded in that compiler note.

The exact signed lattice gate is reproduced by:

```bash
sage -python elkies-k3/scripts/verify_h3_q6_signed_descent_gate.sage
```

It constructs the two inverse section classes with their correct, different
E7 corrections.  The associated primitive isotropic degree-two divisors
`D-` and `D+` satisfy `D-.D+=21`; their sum is primitive with square `42`,
old-fiber degree four, `h0=23`, and arithmetic genus `22`.  Consequently, if
equation-level Galois conjugation exchanges the two signs, the natural
tensor/trace-norm descent has ranks four/three and gives a degree-21 surface
map, not a rank-two elliptic pencil.  The Galois sign rule remains unproved,
and the preceding descent certificate proves that it is not the actual H21
marking: both signs are fixed.  The calculation remains a useful rejection of
that hypothetical shortcut and does not exclude accidental elliptic factors
in the genus-22 Jacobian.  The generated artifact
`artifacts/generated-results/elkies-k3-h3-q6-signed-descent-gate.json` has
SHA-256
`3be2de6e2f7c722bc04dde0ad5eba81924b130b93d6850009cd266398b4b60d7`.

## Preferred q=6 transport from H3

The exact labeled `H3` frame has a much smaller first neighbor than the q80
route.  This command exhausts its q=6 sign-pair shell with marked coordinate
15 equal to one and retains four representative output frames:

```bash
mkdir -p artifacts/local/elkies-k3/h3-q6-mw3-frames
sage -python elkies-k3/scripts/search_alternate_fibrations.sage \
  --frame elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt \
  --min-qnorm 6 \
  --max-qnorm 6 \
  --proper-factors-only \
  --one-factor-order \
  --fixed-coordinate 15:1 \
  --per-root-data-cap 4 \
  --quiet-candidates \
  --report 4 \
  --out artifacts/local/elkies-k3/h3-q6-mw3-search.txt \
  --frames-dir artifacts/local/elkies-k3/h3-q6-mw3-frames
```

The decisive summary is:

```text
FIBSEARCH|stage=restart_done|q=6|restart=0|new_pairs=56|union_pairs=56|nodes=441|exhaustive=true
FIBSEARCH|stage=summary|tested=56|unique_frames=4|sampled_vectors=56
```

All 56 tested vectors have root data `(14,312,3)`, identifying `E8+E6`, and
MW rank three.  To inspect all 56 rows rather than four retained
representatives, omit `--per-root-data-cap 4` and use `--report 56`.  The
four-frame count is a retention cap, not a proved count of Weyl or isometry
orbits.  The clean first witness has `(a,b)=(2,3)` and

```text
(0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,1,0).
```

Check the degree-two chamber interpretation, saturated MW height, distinction
from the old H2/q60 child, and the nef q=8 continuation with:

```bash
sage -python elkies-k3/scripts/analyze_h3_first_q6_chamber.sage
```

The two terminal status lines are:

```text
H3Q6|...|status=PASS
H3Q8|...|child=D13/MW4|root_data=13,312,4|status=PASS
```

The raw q=8 shell is too large to enumerate directly.  This exact Weyl
quotient classifies it in seconds and writes all coordinate changes, neighbor
bases, source-H3 divisor coordinates, and both D13/MW4 child frames:

```bash
sage -python elkies-k3/scripts/classify_h3_q6_child_q8_orbits.sage \
  --output artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json
```

Expected terminal status:

```text
H3Q6Q8|raw_norm16=219758670|horizontal_norm16=139006800|mw_projections=10|dominant_orbits=63|primitive_neighbors=61|mw4_hits=2|status=PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION
```

The checker also verifies the pinned continuation frame
`elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt`.  Its first
13 coordinates are D13 simple roots and its last four quotient directions
have reduced MW height Gram
`[[3/4,1/4,-1/4,0],[1/4,11/4,1/4,1],[-1/4,1/4,11/4,-1],[0,1,-1,46]]`.

Before deriving the global q=8 pencil, reconstruct its exact generic-fibre
ambient directly on the H92 source. The selected class is first source-chamber
reduced by 122 fixed-component reflections; its resulting source-nef class
has generic restriction `9(O)+9(-P1)`, and the 18-element chord basis uses
only the already explicit `P1`.  The artifact first reconstructs the literal
source class as `9O+9(-P1)+V-11F`, with every nonzero simple-component
coefficient of `V` recorded and checked before choosing this basis:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_generic_rr_ambient.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json
```

It reports `H92Q8AMBIENT|degree=18|basis=18|...|status=PASS_EXACT_Q8_GENERIC_RR_AMBIENT`.
This is the generic-fibre ambient, not an equation-level q=8 pencil: the
finite vertical and resolved E7/E8 conditions for the source-nef vertical
divisor remain to be derived.
<!-- status-consumer: EC-K3-H3-Q8-AMBIENT 2e14dd27b9a3dd79 -->

The actual all-edge q6 E7 module also yields the first reusable source-q8
generic-component condition layer. It evaluates the cleared ninth-power
comparison on all seven actual E7 components. On the corrected 54-term seed,
139 negative-order groups occur; 41 singleton groups give a rank-22 exact
coordinate block, while the non-singletons are retained for later
chart-function-field residue calculations.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions.json
```

<!-- status-consumer: EC-K3-H3-Q8-ALL-COMPONENT-GENERIC 6f60b40f3b99f693 -->

The non-singleton groups now have pinned actual generic-component chart
frames: each of the seven components is assigned a blow-up chart, reduced
component equation, and normal `(Z,U,Y)` weight that reproduces its certified
orders of `t,x,y`. In particular, the three `Y=0` components retain the
actual quadratic normal branch rather than a Kodaira-label proxy.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_generic_component_chart_frames.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-generic-component-chart-frames.json
```

<!-- status-consumer: EC-K3-H3-Q8-GENERIC-COMPONENT-CHART-FRAMES 48fc1dbd9486f9a2 -->

The selected q=8 class also has an exact degree-two generic marking directly
on the q=6 child.  The binary-quartic covariant is a 2-covering map, so the
primitive geometric section is `S=Pmap+Qmap`, not the previously used
`2Pmap+2Qmap`.  Its relative Mordell--Weil coordinate is `(-2,-2,0)`, its
height is 24, and it meets the zero section at the degree-ten divisor `h`:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_marking.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json
```

Expected geometry:

```text
MW coordinate = (-2,-2,0)
height = 24
O-intersection = 10
smooth collision degree = 10
II*, IV* = identity component
```

Compile the corrected characteristic-zero pencil and exact child with:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --representative component-nef \
  --output artifacts/local/elkies-k3/q8-target-component-nef-audit.json
cmp artifacts/local/elkies-k3/q8-target-component-nef-audit.json \
  elkies-k3/data/fibrations/h3_q8_component_nef_physical_root_target.json
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_corrected2cover_qq.sage
```

Expected terminal records:

```text
Q8QQRR|ambient=13|rows=11|rank=11|kernel=2|...
Q8QQQUARTIC|degree=4
Q8QQCHILD|finite=[(1,[2,3,15],'I9*'),(9,[0,0,1],'I1')]|infinity=((0,0,0),'smooth')|root_rank=13|root_euler=24|root_det=4|MW_rank=4|status=PASS_EXACT_CORRECTED_Q8_D13_CHILD
```

This proves the exact second H3 neighbour `E8+E6/MW3 -> D13/MW4` over
`QQ`.  The older degree-46 marking and the endpoint-envelope obstruction
experiments below remain historical diagnostics only; they do not describe
the repaired pencil.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-MARKING 745bf011cb47e7f3 -->

Rewrite the same exact divisor in finite component roots relative to the
explicit transported old zero:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-physical-root-target.json
```

Expected terminal status:

```text
H92Q6CHILDQ8ROOTS|E6_cycle=3,5,6,4,2,3|E6_degrees=-1,-1,0,0,0,0|E8_cycle=4,5,7,10,8,6,4,2|E8_degrees=-1,0,0,0,0,0,0,0|status=PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET
```

This pins only the lattice-component target. The simple-root labels have not
yet been mapped to resolved II*/IV* charts, and no finite quotient module or
global q8 pencil is asserted.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-PHYSICAL-ROOT-TARGET 064318c2afe537fd -->

For the Weyl-nef q8 fibre (not the dominant D13 image), use the invariant IV*
ideal `(u^2,X,Y)`. Its finite q-regular module is
`<(1,lift(R/Nx)),(0,f_II^2*f_IV^2)>`; the infinity lattice remains open.

<!-- status-consumer: EC-K3-H3-Q8-CHILD-NEF-LOCAL-MODULE e2887bd2bd4f6c27 -->

For the explicit q6 child, make the abstract Weyl-nef q8 class nonnegative on
the actual `E6+E8` components, then certify all remaining fixed-component
walls.  The component reduction takes 102 reflections.  The certificate
enumerates every possible negative old-fibre section wall in its exact
rank-17 short coset; a separate parity identity excludes bisection walls.
Together with the finite and affine component degrees, this proves the
resulting primitive isotropic class is nef and defines a genus-one pencil of
old-fibre degree two.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --representative component-nef \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-physical-root-target.json
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_section_walls.sage \
  --target artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-physical-root-target.json \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-section-walls.json
```

Expected terminal status:

```text
H92Q6CHILDQ8SECTIONS|representative=component-nef|short=274731|coset_short=0|negative=0|vertical_nef=1|status=PASS_PRIMITIVE_NEF_DEGREE_TWO_CLASS
```

This is a lattice bisection pencil only.  Its standard-Weierstrass equation
also needs the exact NS transport of the translation from the transported old
zero to the Weierstrass infinity section; the existing marked chord belongs
to the translated divisor. Its branch divisor, quadratic-extension hash,
collision analysis, and rank remain open.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-BISECTION-PENCIL 84142196f27e5e2d -->

The exact generic chord for that physical divisor is now available. Let `P0`
be the transported old zero as a finite standard-Weierstrass point, let `S`
be the existing marked point, and set `Q=P0+S`. The divisor is `P0+Q` and
translation by `-P0` carries it to `O_standard+S`; pull the standard chord
back through that translation:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_component_nef_chord.sage
```

Expected terminal status:

```text
H92Q6CHILDCOMPNEFCHORD|translation=tau_-P0|generic_basis=2|status=PASS_EXACT_COMPONENT_NEF_OLD_ZERO_CHORD
```

The translated resolved local modules and infinity module are still required
for an equation or branch divisor.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-CHORD-TRANSPORT c3896078bcdd432d -->

Construct and screen level zero of the generic chord over the old base:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_component_nef_bisection_branch.sage
```

Expected terminal status:

```text
H92Q6CHILDCOMPNEFBRANCH|level=0|quadratic=1|canonicalized=True|branch_degree=192|status=PASS_EXACT_COMPONENT_NEF_GENERIC_LEVEL_BRANCH
```

The result is a quadratic in the translated coordinate `x'`; its physical
equation is its pullback by `tau_-P0`. Its exact squareclass has degree-192
finite branch divisor, so this level is rejected as a rational-bisection
collision candidate.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-GENERIC-BISECTION-BRANCH 6d0dbb4b90b14710 -->

At both additive cusps, the exact chord denominator is a unit; translation by
`-P0` has the prescribed additive-order cusp image and a unipotent
determinant-one tangent action.  In particular the translation centre is
smooth at both II* and IV*, so translation by `-P0` extends on their Néron
smooth loci:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_component_nef_translation_additive_jets.sage
```

Expected terminal status:

```text
H92Q6CHILDCOMPNEFJETS|II*=P0_smooth|IV*=P0_smooth|status=PASS_EXACT_COMPONENT_NEF_TRANSLATION_ADDITIVE_JETS
```

This is a singular-germ and tangent prerequisite only; resolved blow-up-chart
pullbacks, the infinity module, pencil equation, and branch divisor remain
open.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-TRANSLATION-ADDITIVE-JETS f76374817b493c9a -->

The II* part of that target has a unit-normalized complete ideal:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_iistar_vertical_ideal.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-iistar-vertical-ideal.json
```

It reports `H92Q6CHILDQ8II|ideal=(u2,X,Y)|colength=2|...`. This fixes only
the vertical II* ideal; the generic-chord trivialization and IV* module still
have to be derived before forming a q8 global condition matrix.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-IISTAR-VERTICAL-IDEAL ea029da7275e86da -->

The IV* vertical condition is determined up to the unresolved E6 arm
orientation.  Its two conjugate colength-four ideals are
`(Y-c*u^2,u*X,X^2,u^3)` and `(Y+c*u^2,u*X,X^2,u^3)`, where `c^2=b(0)` in the
unit-normalized IV* germ.  Do not select either one as the q8 module until the
physical E6-root-to-chart attachment is derived.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_ivstar_vertical_ideal.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-vertical-ideal.json
```

It reports `H92Q6CHILDQ8IV|ideals=2|colength=4|orientation=unresolved_E6_arm|...`.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-VERTICAL-IDEAL-PAIR 4300b94b65dabc64 -->

The transported old `E7_7` section orients this pair: it meets physical E6
root five and the IV* `Y/u^2=c` branch.  The selected vertical condition is
`(Y+c*u^2,u*X,X^2,u^3)`, with colength four.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_ivstar_orientation.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-orientation.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-ORIENTATION 4e931cdcb76201c3 -->

The actual generic chord `m=(y+y(S))/(x-x(S))` now has exact finite
coefficient blocks in both selected additive quotients.  They have ranks two
at II* and four at IV*; at IV*, reducing `Y=-c*u^2` contributes the essential
`u^2` correction to the residue of `m`.  Both blocks are compiled through the
resolved marked-chord quotient interface with their actual ideals and
normal-form bases. This is local jet data, not a global base-function bound or
a q8 pencil.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_additive_chord_blocks.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-additive-chord-blocks.json
```

Expected terminal status:

```text
H92Q6CHILDQ8ADDITIVE|II_rows=2|IV_rows=4|II_rank=2|IV_rank=4|status=PASS_EXACT_Q6_CHILD_Q8_ADDITIVE_CHORD_BLOCKS
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-ADDITIVE-CHORD-BLOCKS bf6fc5b74f51fc0f -->

As bounded modular reconnaissance only, the first saturated coefficient window
already has the expected dimension: at each of `p=43,53,59`, take
`B in <1,T,...,T^7>` and no `h^2*C` correction.  The smooth congruence plus
the six additive rows has rank six in dimension eight, hence kernel dimension
two.  The neighbouring `B` degrees six and eight give dimensions one and
three.  These degree bounds have not been derived in characteristic zero, so
this does not certify a q8 pencil.

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_ansatz.sage \
  --prime 43 --max-b-degree 7 --max-c-degree -1 \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-saturated-ansatz-probe.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-SATURATED-ANSATZ-PROBE cdcfe4ca7989839d -->

Screen every finite level of the resulting deterministic two-dimensional
ratio for the degree-four branch divisor required of a genus-one double
cover.  This is a three-prime bounded modular obstruction, not a global q8
pencil result:

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_pencil_modp.sage --prime 43
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_pencil_modp.sage --prime 53
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_pencil_modp.sage --prime 59
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-SATURATED-PENCIL-OBSTRUCTION 5a4beade3e2caff1 -->

Compute the entire finite coefficient module, rather than a bounded
polynomial window.  The expected profile at each listed good prime is six
independent finite rows and Smith degrees `(1,5)`:

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_finite_module_modp.sage --prime 43
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_finite_module_modp.sage --prime 53
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_finite_module_modp.sage --prime 59
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-FINITE-MODULE-MODULAR 78fc7f298da9eaf6 -->

Derive the complete smooth O.S collision module for that marking:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_smooth_collision_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-smooth-collision-module.json
```

Expected terminal status:

```text
H92Q6CHILDQ8SMOOTH|degree_h=46|squarefree=1|base_regular=h_divides_b|saturated_quotient=92|status=PASS_EXACT_Q6_CHILD_Q8_SMOOTH_COLLISION_MODULE
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-SMOOTH-MODULE 559fade7fa0cb618 -->

The former depth-two IV* entrance and branch-orientation scripts are
historical only: the corrected marking is smooth at IV*, so neither is a
valid q8 local-module input.
<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-ENTRANCE 86bf6cc2487fd0b4 -->
<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-COMPONENT ee2d661e3dc73ef0 -->

The first resolved q=8 local target is also exact.  It identifies the source
and chart-resolved E7 component orders and proves that the q=8 E7 class is the
ninth tensor power of `Z*J_-P1^dual` followed by the integral exceptional
twist `(2,5,6,4,6,3,5)`.  The correction has exceptional intersection degrees
`(0,1,0,0,-7,0,1)`, so it is not anti-nef: a later quotient compiler must use
resolved-chart gluing rather than substitute one complete ideal downstairs:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_local_target.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-local-target.json
```

Expected terminal status:

```text
H92Q8E7TARGET|degrees=0,1,0,0,2,0,1|twist=2,5,6,4,6,3,5|status=PASS_EXACT_Q8_E7_LOCAL_TARGET
```

<!-- status-consumer: EC-K3-H3-Q8-E7-TARGET a6dd94428dfa14e4 -->

That non-anti-nef correction is now attached to the six actual H92 E7 edge
charts.  The generated factors are line-bundle transitions, not a substitute
complete ideal: if `g` is the displayed factor then the q8 chart
representative `f` must satisfy `g*f` in the ninth q6 marked module.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_actual_e7_gluing.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json
```

Expected terminal status:

```text
H92Q8ACTUALE7GLUING|charts=6|twist=2,5,6,4,6,3,5|status=PASS_EXACT_Q8_ACTUAL_E7_GLUING
```

For example, on the actual `E7_2--E7_5` chart the factor is `Z^6*Y^5`.
This fixes the resolved-chart orientation required before a finite q8 E7
condition matrix can be formed.

The fractional ninth power can now be supplied to that compiler without
collapsing its marked branch. With `J_P1=(x-xP1,y-yP1)`, clear the q6
denominator chartwise by testing `(x-xP1)^9*g*f` against the ten generators
of `t^9*J_P1^9`. The following compact manifest records those generators,
all six actual pullbacks, and the cleared template for every generic q8 frame
term; it is not yet a membership matrix.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_actual_e7_power_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-actual-e7-power-module.json
```

Expected terminal status:

```text
H92Q8E7POWERMODULE|charts=6|power_generators=10|generic_terms=18|status=PASS_EXACT_Q8_E7_CLEARED_POWER_MODULE
```

It preserves the non-anti-nef gluing problem for the five non-marked edges;
in particular it does not replace the target by a scalar complete ideal.

With the current 54-term source envelope and the deeper `h^-15` smooth
frame, the declared `r=7` enlargement has a 558-column smooth block with a
two-dimensional kernel modulo each of `43,53,59,89`. The generic `E7_3`
condition of the actual
`E7_4--E7_3` chart eliminates both directions at all four primes. This is a
one-chart modular obstruction, not a characteristic-zero pencil result.

```bash
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --extra-h-power 7 --include-kernel \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra7.json
sage -python elkies-k3/scripts/probe_h92_q8_e3_generic_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra7.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e3-generic-module-mod-43-extra7.json
```

Expected terminal status from the second command:

```text
H92Q8E3GENERIC|prime=43|smooth_kernel=2|constraints=2|survivor=0|status=EXPERIMENTAL_MODULAR_E3_GENERIC_MODULE_OBSTRUCTION
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA7-E3-GENERIC-OBSTRUCTION 41e4df0f0b84986e -->

The all-component singleton layer independently gives a mod-43 obstruction
to this corrected `r=7` smooth kernel: its six exact singleton coordinate
rows restrict with rank two to the two-dimensional smooth kernel. Thus no
direction remains even before resolving non-singleton residues or edge nodes.

```bash
sage -python elkies-k3/scripts/assemble_h92_q8_endpoint_rr_ambient.sage \
  --extra-h-power 7 \
  --output artifacts/generated-results/elkies-k3-h92-q8-extra7-endpoint-rr-ambient.json
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage \
  --ambient artifacts/generated-results/elkies-k3-h92-q8-extra7-endpoint-rr-ambient.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra7.json
sage -python elkies-k3/scripts/probe_h92_q8_all_component_singleton_modp.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-singleton-mod-43-extra7.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA7-ALL-COMPONENT-SINGLETON-OBSTRUCTION 2b07df0cd2f4e19f -->

The broader declared `r=10` envelope has 774 columns and a 16-dimensional
smooth kernel modulo each of `43,53,59,89`. Necessary generic conditions on
all six unmarked E7 components leave four directions, but the actual
`E7_4--E7_3` node has four successive unique Pareto-leading negative terms
and eliminates all four. This is again a bounded modular leading-term
obstruction, not a complete resolved-chart or characteristic-zero result.

```bash
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --extra-h-power 10 --include-kernel \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_unmarked_e7_generic_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_e7_4_3_node_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json \
  --generic artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-module-mod-43-extra10.json
```

Expected terminal status from the final command:

```text
H92Q8E743NODE|prime=43|generic_candidate=4|constraints=4|survivor=0|status=EXPERIMENTAL_MODULAR_E7_4_3_NODE_MODULE_OBSTRUCTION
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA10-NODE-OBSTRUCTION 3152bae4a804d3d8 -->

The declared `r=13,16,19` envelopes are also excluded modulo `43`. Their
smooth-kernel/generic-survivor dimensions are respectively `38/26`, `50/38`,
and `50/38`; the `E7_4--E7_3` node eliminates every survivor.

```bash
for r in 13 16 19; do
  sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
    --prime 43 --extra-h-power "$r" --include-kernel \
    --output "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r.json"
  sage -python elkies-k3/scripts/probe_h92_q8_unmarked_e7_generic_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r.json"
  sage -python elkies-k3/scripts/probe_h92_q8_e7_4_3_node_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r.json" \
    --generic "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-module-mod-43-extra$r.json"
done
```

Expected final status:

```text
H92Q8E743NODE|prime=43|generic_candidate=26|constraints=26|survivor=0|status=EXPERIMENTAL_MODULAR_E7_4_3_NODE_MODULE_OBSTRUCTION
```

This is bounded and modular only.
<!-- status-consumer: EC-K3-H3-Q8-EXTRA13-NODE-OBSTRUCTION 0cd4414b09443b94 -->

Allowing one individual E7-pole slack unit does not open the initial windows:
`r=4` has zero smooth kernel, `r=7` has `9 -> 0` after unmarked generic E7
rows, and `r=10` has `25 -> 13 -> 0` after generic then node rows.

```bash
for r in 4 7 10; do
  sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
    --prime 43 --extra-h-power "$r" --extra-e7-pole 1 --include-kernel \
    --output "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r-e7slack1.json"
  sage -python elkies-k3/scripts/probe_h92_q8_unmarked_e7_generic_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r-e7slack1.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r-e7slack1.json"
  sage -python elkies-k3/scripts/probe_h92_q8_e7_4_3_node_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r-e7slack1.json" \
    --generic "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r-e7slack1.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-module-mod-43-extra$r-e7slack1.json"
done
```

These are bounded one-prime necessary-condition screens only.
<!-- status-consumer: EC-K3-H3-Q8-E7-SLACK1-OBSTRUCTION 419a01f2b59d42ea -->

Before using that node screen as anything more than a leading-term
obstruction, the local module itself is now transported from the actual H92
chart.  On `E7_4--E7_3`, the surface equation is `Y^2-U*H(Z,U)=0` with
`H(0)=1`; the exact P1 entrance orders show
`x-x(P1)=x*unit` and `y-y(P1)=y*unit`.  Thus the q6 module is `t*R`, and the
q8 node condition has the actual frame `(Z^4*Y^6)*f/t^9 in R`.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_4_3_node_frame.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-frame.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-3-NODE-FRAME c3d6b194bb72dd59 -->

The same chart-level calculation now covers all five unmarked E7 edge nodes.
Their q6 module is uniformly `t*R`; the actual q8 Cartier factors are
`U^4Y^2`, `Z^4Y^6`, `U^5Y^6`, `Z^5Y^5`, and `Z^3Y^6` in edge-chart order.
The sole exception remains the marked `E7_2--E7_5` node, whose leading P1
cancellation is deliberately not forced into this principal frame.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_unmarked_e7_node_frames.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-node-frames.json
```

<!-- status-consumer: EC-K3-H3-Q8-UNMARKED-E7-NODE-FRAMES 59977f99242c92c6 -->

The remaining `E7_2--E7_5` node has the required leading cancellation, so it
is checked separately rather than treated by strict orders.  The actual
second-U chart gives `x-x(P1)=Z^3U^2*unit` and
`y-y(P1)=Z^3U^2Y*unit`; consequently its node module is also `t*R` and its
q8 condition is `Z^6Y^5*f/t^9 in R`.  This is distinct from the marked
smooth point `-P1` on E7₅, whose frame remains separate.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_2_5_node_frame.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-2-5-node-frame.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-2-5-NODE-FRAME eaff4f299d528774 -->

All six frames now feed one exact two-parameter node template.  For each
term `u^i*x^a*m^b/h^k`, it records the actual leading bidegree
`ord(g)+(4k-i-9)ord(t)+a ord(x)+b ord(m)`.  The seed has 260 negative groups
(196 initially singleton), but Pareto-minimality leaves one independent exact
initial coordinate row: five charts see the same coefficient. Equal or
non-minimal groups are not promoted to rows until the finite local quotient is
evaluated.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_all_e7_node_principal_bidegrees.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-e7-node-principal-bidegrees.json
```

<!-- status-consumer: EC-K3-H3-Q8-ALL-E7-NODE-PRINCIPAL-BIDEGREE-TEMPLATE f9faa97d398f8f02 -->

That safe Pareto row is now materialized through the shared resolved-condition
compiler.  Its normalized `1 x 54` matrix has rank one and kernel dimension
53.  It is intentionally only the initial node block, not the full node
quotient.

```bash
sage -python elkies-k3/scripts/compile_h92_q8_initial_node_conditions.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-initial-node-conditions.json
```

<!-- status-consumer: EC-K3-H3-Q8-INITIAL-NODE-CONDITION-BLOCK a2294ddca25d7344 -->

For the remaining `E7_4--E7_3` node quotient, the actual principal frame now
clears every unit denominator before any ideal computation.  On the 54-column
seed, `g*f/t^9` is regular exactly when the common-cleared numerator lies in
`t^17 R`; this is a single actual-chart divisibility problem, not the raw
ten-generator singular-model power ideal.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_3_principal_node_clearing.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-clearing.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-3-PRINCIPAL-NODE-CLEARING dc4e3312a8eb2cef -->

This does not make `R/(t^17)` finite.  In the actual completed chart,
`t^17=Z^51Y^68*unit`, so that quotient has Krull dimension one.  The finite
corner jet `(Z^51,Y^68)` has length 3468 but is not product divisibility;
the following guard records that distinction before a residual quotient is
introduced.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_4_3_node_divisibility_geometry.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-divisibility-geometry.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-3-NODE-DIVISIBILITY-GEOMETRY e287ab81214ee330 -->

There is now a chart-faithful modular prototype for this principal condition.
It clears the certified units, discards only terms already in the actual
ideal `(t^17)`, and reduces each ambient numerator by a local (`ds` order)
standard basis of `(surface,t^17)`.  Its finite coordinate space is the
image of the supplied ambient—not a finite quotient of `R/(t^17)` and not
the rectangular corner jet.  On the 54-column seed modulo `43`, the image
has 17,612 displayed normal monomials and rank 54 (zero kernel); this is
consistent with, but weaker in scope than, the existing generic-E7 rejection.

```bash
sage -python elkies-k3/scripts/probe_h92_q8_e7_4_3_principal_node_local_normal_form_modp.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-local-normal-form-mod-43.json
```

The script accepts an enlarged endpoint ambient as well.  It remains a
single-prime local regression until the corresponding characteristic-zero
and all-chart overlap maps are supplied.

For the pinned 54-column seed, the full-rank local image now has the standard
good-reduction consequence: after primitive normalization, any nonzero
characteristic-zero kernel vector would reduce to a nonzero vector in the
mod-43 local kernel.  The source denominators and every common-clearing
factor are units at the node modulo `43`, so this yields an exact bounded
injectivity certificate for the actual principal node condition:

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_4_3_principal_node_good_reduction.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-good-reduction.json
```

```text
H92Q8E743LOCALGOODREDUCTION|prime=43|ambient=54|rank=54|kernel=0|status=PASS_EXACT_Q8_E7_4_3_PRINCIPAL_NODE_INJECTIVITY
```

This rejects that seed at this one genuine resolved node even without its
already-known smooth and generic-E7 obstructions. It does not give a
characteristic-zero coordinate matrix, an enlarged-ambient result, or the
remaining node/overlap maps.

The same normalized chord calculation is now replayed on all six actual E7
edge charts, including the cancellation-sensitive `E7_2--E7_5` node.  It
records the local formula for `m`, the unit multiplier, and the common
principal target `(t^17)` for each chart; it does not replace the remaining
local images or their overlap compatibility with component data.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_node_principal_clearings.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-node-principal-clearings.json
```

```text
H92Q8E7NODECLEARINGS|nodes=6|ambient=54|T=17|K=6|status=PASS_EXACT_Q8_E7_NODE_PRINCIPAL_CLEARINGS
```

The shared modular evaluator consumes this atlas directly.  These two runs
check the ordinary and cancellation-sensitive node charts without replacing
their actual local quotients by finite jets.  Both have full 54-column rank
modulo `43` (the normal-form coordinate counts are chart-dependent).

```bash
sage -python elkies-k3/scripts/probe_h92_q8_e7_node_principal_local_normal_form_modp.sage \
  --chart E7_4--E7_3 \
  --output artifacts/local/elkies-k3-h92-q8-e7-4-3-generic-node-local-mod-43.json
sage -python elkies-k3/scripts/probe_h92_q8_e7_node_principal_local_normal_form_modp.sage \
  --chart E7_2--E7_5 \
  --output artifacts/local/elkies-k3-h92-q8-e7-2-5-generic-node-local-mod-43.json
```

```text
H92Q8E7NODELOCALNF|chart=E7_4--E7_3|prime=43|ambient=54|rows=15636|rank=54|kernel=0|status=EXPERIMENTAL_MODULAR_Q8_E7_NODE_LOCAL_NORMAL_FORM_BLOCK
H92Q8E7NODELOCALNF|chart=E7_2--E7_5|prime=43|ambient=54|rows=23553|rank=54|kernel=0|status=EXPERIMENTAL_MODULAR_Q8_E7_NODE_LOCAL_NORMAL_FORM_BLOCK
```

The same default local-degree computation has full rank for
`E7_1--E7_4`, `E7_3--E7_7`, and `E7_7--E7_2`; their normal-form coordinate
counts are respectively `11841`, `19669`, and `16757`.  The remaining
`E7_3--E7_6` local standard basis is deliberately expensive.  The command
therefore provides a distinct, one-way Artinian-corner obstruction mode:
the finite quotient `(surface,Z^34,U^34)` contains the actual `(t^17)`
ideal and is supported at the chart origin.  Thus every true local solution
maps to zero, while the converse is intentionally not used.

```bash
sage -python elkies-k3/scripts/probe_h92_q8_e7_node_principal_local_normal_form_modp.sage \
  --chart E7_3--E7_6 --mode finite-corner-obstruction \
  --output artifacts/local/elkies-k3-h92-q8-e7-3-6-corner-obstruction-mod-43.json
```

```text
H92Q8E7NODECORNER|chart=E7_3--E7_6|prime=43|ambient=54|rows=1806|rank=54|kernel=0|status=EXPERIMENTAL_MODULAR_Q8_E7_NODE_FINITE_CORNER_OBSTRUCTION
```

The good-reduction checker verifies the corner containment and support from
the actual resolved chart, then promotes this full-rank finite image to a
characteristic-zero injectivity statement for true local relations in the
fixed 54-column ambient.  Its conclusion remains one-way and does not call
the corner the q8 node quotient.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_3_6_corner_obstruction_good_reduction.sage \
  --corner artifacts/local/elkies-k3-h92-q8-e7-3-6-corner-obstruction-mod-43.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-3-6-corner-obstruction-good-reduction.json
```

```text
H92Q8E736CORNERGOODREDUCTION|prime=43|ambient=54|rank=54|kernel=0|status=PASS_EXACT_Q8_E7_3_6_CORNER_OBSTRUCTION_INJECTIVITY
```

The all-node regression combines the five actual local normal-form images
with that sixth-chart one-way obstruction.  It records the distinction rather
than calling the finite corner a node quotient.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_all_e7_node_modular_regression.sage \
  --node-image artifacts/local/elkies-k3-h92-q8-e7-1-4-generic-node-local-mod-43.json \
  --node-image artifacts/local/elkies-k3-h92-q8-e7-2-5-generic-node-local-mod-43.json \
  --node-image artifacts/local/elkies-k3-h92-q8-e7-3-6-corner-obstruction-mod-43.json \
  --node-image artifacts/local/elkies-k3-h92-q8-e7-3-7-generic-node-local-mod-43.json \
  --node-image artifacts/local/elkies-k3-h92-q8-e7-4-3-generic-node-local-mod-43.json \
  --node-image artifacts/local/elkies-k3-h92-q8-e7-7-2-generic-node-local-mod-43.json \
  --output artifacts/local/elkies-k3-h92-q8-all-e7-node-resolved-obstruction-mod-43.json
```

```text
H92Q8ALLE7NODERESOLVED|prime=43|nodes=6|ambient=54|status=EXPERIMENTAL_MODULAR_Q8_ALL_E7_NODE_RESOLVED_OBSTRUCTION
```

These are finite-ambient modular regressions only. The characteristic-zero
residual matrices and the overlap-compatible common kernel still need to be
compiled.

Two actual sibling-chart overlap maps are now pinned as well.  The U- and
Z-chart pairs `E7_1--E7_4`/`E7_4--E7_3` and
`E7_3--E7_7`/`E7_7--E7_2` are related on `Z != 0` by
`(Z,U,Y) -> (U*Z,1/Z,Y/Z)`.  The checker verifies the old H92 `(t,x,y)`
pullbacks identically across each overlap and records the q8 Cartier-frame
ratios `Y^4/Z^2` and `1/Y`.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_sibling_chart_transitions.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-sibling-chart-transitions.json
```

These are transition functions, not units or an already-complete Čech
matrix. The other component overlaps and finite residual maps remain needed.

At the recorded prime `43`, the actual all-component generic screen is
stronger still: it includes marked `E7_5` with its audited exact order
`ord(m)=0`.  The successive unique-live leading-coefficient cuts have rank
16 on the 16-dimensional smooth kernel, so the `r=10` envelope has no
survivor before a node condition is used. This remains a one-prime bounded
generic-point obstruction, not a characteristic-zero or complete-chart
result.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage \
  --ambient artifacts/generated-results/elkies-k3-h92-q8-extra10-endpoint-rr-ambient.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_all_component_generic_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json \
  --conditions artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-module-mod-43-extra10.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA10-ALL-COMPONENT-GENERIC-OBSTRUCTION 8fa617eca25dc71e -->

The first non-singleton residue compiler goes beyond those valuation rows on
the actual `E7_4--E7_3` and `E7_3--E7_7` charts. It normalizes the components
using their resolved surface equations, substitutes the audited chord, and
splits each normal-order group by its component-parameter power.  On the
`r=10` ambient this yields 124 exact rows on `E7_4` and 131 on `E7_7`; after
reduction mod `43`, their restriction to the 16-dimensional smooth kernel
has rank 16 and no survivor, without using a node condition.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_7_generic_residue_rows.sage \
  --conditions artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residue-rows-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_e7_4_7_generic_residues_modp.sage \
  --residues artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residue-rows-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residues-mod-43-extra10.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-7-GENERIC-RESIDUES f4ea27a1597d3df7 -->
<!-- status-consumer: EC-K3-H3-Q8-EXTRA10-E7-4-7-GENERIC-RESIDUE-OBSTRUCTION 146750c0a39e15f4 -->

The conic components are also now evaluated in their actual coordinate rings,
not through an inferred component parameter. Setting `Z=0` on the
`E7_2--E7_5` and `E7_3--E7_6` charts gives exact rings
`QQ(U,Y)/(F(0,U,Y))`; the compiler clears a stated common denominator and
reduces there. The 54-column seed gives 228 E7₅ rows and 57 E7₆ rows. At E7₅
the calculation retains the audited `x-x(P1)` and `y-y(P1)` cancellation.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_5_6_generic_residue_rows.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-5-6-generic-residue-rows.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-5-6-GENERIC-RESIDUES de263de56ed5c2dc -->

The complete generic-component residue cover now uses actual chart equations
on all seven E7 components of the base 54-column seed.  The Y-branch solver
covers E7₁--E7₃ (including the E7₂ entrance cancellation), the conic rings
cover E7₅--E7₆, and the existing normalized charts cover E7₄ and E7₇.  Their
exact residue rows total 983, distributed as
`(42,189,391,42,228,57,34)` in component order.  This certifies generic
component conditions only; boundary nodes, the marked branch, overlaps, and
global kernel computation remain outstanding.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_1_3_generic_residue_rows.sage
sage -python elkies-k3/scripts/derive_h92_q8_e7_5_6_generic_residue_rows.sage
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_7_generic_residue_rows.sage
sage -python elkies-k3/scripts/certify_h92_q8_all_generic_e7_residue_cover.sage
```

Expected terminal status from the final command:

```text
H92Q8ALLGENERICRESIDUES|components=7|rows=983|status=PASS_EXACT_Q8_ALL_GENERIC_E7_RESIDUE_COVER
```

<!-- status-consumer: EC-K3-H3-Q8-ALL-GENERIC-E7-RESIDUE-COVER e3ea804e3f4dd675 -->

The compiler now stacks those 983 resolved-chart residue rows with the 22
singleton generic-component rows. On the least 54-column endpoint ambient,
the resulting 1005-by-54 characteristic-zero matrix has rank 54 and zero
kernel. Thus this bounded ambient cannot contain the q8 pencil even before
node, marked-branch, overlap, E8, or smooth conditions are imposed; this does
not rule out the required enlarged ambient.

```bash
sage -python elkies-k3/scripts/compile_h92_q8_generic_component_conditions.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-generic-component-condition-block.json
```

```text
H92Q8GENERICCOMPONENTBLOCK|ambient=54|singleton_rows=22|residue_rows=983|rank=54|kernel=0|status=PASS_EXACT_Q8_GENERIC_COMPONENT_CONDITION_BLOCK
```

<!-- status-consumer: EC-K3-H3-Q8-GENERIC-COMPONENT-CONDITION-BLOCK 33693f196eb13091 -->

The first enlargement with a nonzero smooth modular kernel is also rejected
exactly. At `extra_h_power=7`, the 558-column smooth block has rank 556 mod
43; all 2,487 generic E7 rows restrict with rank two to that kernel. The
stacked good reduction is therefore full column rank, proving that the
characteristic-zero smooth-plus-generic block has zero kernel for this
particular enlargement.

```bash
mkdir -p artifacts/local/elkies-k3/q8-extra7
sage -python elkies-k3/scripts/assemble_h92_q8_endpoint_rr_ambient.sage --extra-h-power 7 --output artifacts/local/elkies-k3/q8-extra7/endpoint.json
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage --ambient artifacts/local/elkies-k3/q8-extra7/endpoint.json --output artifacts/local/elkies-k3/q8-extra7/generic-template.json
sage -python elkies-k3/scripts/derive_h92_q8_e7_1_3_generic_residue_rows.sage --conditions artifacts/local/elkies-k3/q8-extra7/generic-template.json --output artifacts/local/elkies-k3/q8-extra7/e7-1-3-rows.json
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_7_generic_residue_rows.sage --conditions artifacts/local/elkies-k3/q8-extra7/generic-template.json --output artifacts/local/elkies-k3/q8-extra7/e7-4-7-rows.json
sage -python elkies-k3/scripts/derive_h92_q8_e7_5_6_generic_residue_rows.sage --conditions artifacts/local/elkies-k3/q8-extra7/generic-template.json --output artifacts/local/elkies-k3/q8-extra7/e7-5-6-rows.json
sage -python elkies-k3/scripts/certify_h92_q8_all_generic_e7_residue_cover.sage --y-branch artifacts/local/elkies-k3/q8-extra7/e7-1-3-rows.json --simple artifacts/local/elkies-k3/q8-extra7/e7-4-7-rows.json --conics artifacts/local/elkies-k3/q8-extra7/e7-5-6-rows.json --allow-enlarged --output artifacts/local/elkies-k3/q8-extra7/generic-cover.json
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage --extra-h-power 7 --include-kernel --output artifacts/local/elkies-k3/q8-extra7/smooth-mod43.json
sage -python elkies-k3/scripts/certify_h92_q8_smooth_generic_good_reduction.sage --ambient artifacts/local/elkies-k3/q8-extra7/endpoint.json --template artifacts/local/elkies-k3/q8-extra7/generic-template.json --cover artifacts/local/elkies-k3/q8-extra7/generic-cover.json --smooth-probe artifacts/local/elkies-k3/q8-extra7/smooth-mod43.json --output artifacts/generated-results/elkies-k3-h92-q8-extra7-smooth-generic-good-reduction.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA7-SMOOTH-GENERIC-REJECTION e0e67861f23d4f24 -->

The next endpoint-compatible enlargement is rejected by the same exact
good-reduction calculation. At `extra_h_power=8`, the 630-column smooth block
has rank 624 mod 43, and the complete generic E7 cover has rank six on its
six-dimensional kernel. Consequently the stacked smooth-plus-generic block is
full rank modulo 43 and therefore over QQ. This rejects this bounded ambient
only; node, marked-branch, overlap, and E8 conditions are still not a full q8
compiler.

```bash
mkdir -p artifacts/local/elkies-k3/q8-extra8
sage -python elkies-k3/scripts/assemble_h92_q8_endpoint_rr_ambient.sage --extra-h-power 8 --output artifacts/local/elkies-k3/q8-extra8/endpoint.json
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage --ambient artifacts/local/elkies-k3/q8-extra8/endpoint.json --output artifacts/local/elkies-k3/q8-extra8/generic-template.json
sage -python elkies-k3/scripts/derive_h92_q8_e7_1_3_generic_residue_rows.sage --conditions artifacts/local/elkies-k3/q8-extra8/generic-template.json --output artifacts/local/elkies-k3/q8-extra8/e7-1-3-rows.json
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_7_generic_residue_rows.sage --conditions artifacts/local/elkies-k3/q8-extra8/generic-template.json --output artifacts/local/elkies-k3/q8-extra8/e7-4-7-rows.json
sage -python elkies-k3/scripts/derive_h92_q8_e7_5_6_generic_residue_rows.sage --conditions artifacts/local/elkies-k3/q8-extra8/generic-template.json --output artifacts/local/elkies-k3/q8-extra8/e7-5-6-rows.json
sage -python elkies-k3/scripts/certify_h92_q8_all_generic_e7_residue_cover.sage --y-branch artifacts/local/elkies-k3/q8-extra8/e7-1-3-rows.json --simple artifacts/local/elkies-k3/q8-extra8/e7-4-7-rows.json --conics artifacts/local/elkies-k3/q8-extra8/e7-5-6-rows.json --allow-enlarged --output artifacts/local/elkies-k3/q8-extra8/generic-cover.json
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage --extra-h-power 8 --include-kernel --output artifacts/local/elkies-k3/q8-extra8/smooth-mod43.json
sage -python elkies-k3/scripts/certify_h92_q8_smooth_generic_good_reduction.sage --ambient artifacts/local/elkies-k3/q8-extra8/endpoint.json --template artifacts/local/elkies-k3/q8-extra8/generic-template.json --cover artifacts/local/elkies-k3/q8-extra8/generic-cover.json --smooth-probe artifacts/local/elkies-k3/q8-extra8/smooth-mod43.json --output artifacts/generated-results/elkies-k3-h92-q8-extra8-smooth-generic-good-reduction.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA8-SMOOTH-GENERIC-REJECTION 0d66d55dd153a589 -->

Applying the same eight-command sequence with every `q8-extra8` path and
`extra_h_power 8` replaced by `q8-extra9` and `extra_h_power 9` produces the
702-column endpoint certificate
`artifacts/generated-results/elkies-k3-h92-q8-extra9-smooth-generic-good-reduction.json`.
Here the smooth rank is 690 mod 43 and the generic restriction has rank 12 on
the 12-dimensional smooth kernel, again giving full column rank. This is a
separate bounded-ambient rejection, not a q8 pencil or a replacement for
node/overlap conditions.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA9-SMOOTH-GENERIC-REJECTION 8035f0b465038995 -->

The same full resolved E7 cover rejects the next `r=11` envelope: its
846-column smooth block has rank 822 modulo `43`, and the generic restriction
has rank 24 on the 24-dimensional smooth kernel.  Thus the stacked matrix is
full rank modulo `43` and over `QQ`.  Use the preceding eight-command sequence
with every `q8-extra9` path and `extra_h_power 9` replaced by `q8-extra11` and
`extra_h_power 11`; the final certificate is
`artifacts/generated-results/elkies-k3-h92-q8-extra11-smooth-generic-good-reduction.json`.
This remains a bounded ambient rejection before the node, marked-branch,
overlap, and E8 layers.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA11-SMOOTH-GENERIC-REJECTION c755bb2f34d68f41 -->

At the marked smooth point `-P1` on that same chart, the q8 generic basis is
normalized using the corrected q6 generator `Z*m/t=unit/W`.  The `m^b`
generators have local forms `m^b/t^6`, and the `x*m^b` generators have forms
`x*m^b/t^8`.  After multiplying by `Z^6*Y^5`, each is a unit times
`(Z*m/t)^b` in the ninth q6 marked module:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_actual_e7_marked_frame.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-actual-e7-marked-frame.json
```

Expected terminal status:

```text
H92Q8ACTUALE7MARKED|basis=18|m_denominator=6|xm_denominator=8|status=PASS_EXACT_Q8_ACTUAL_E7_MARKED_FRAME
```

This is an actual marked-chart frame and records pole orders through nine;
the other E7 edges and global compatibility are still required.

Combining those marked-E7 bounds with the actual E8 floors gives the first
finite source-q8 coefficient envelope.  For each generic `x^a*m^b`, it uses
all `u^i/h(u)^k` with the least `k` satisfying
`e8_floor <= i <= 4*k + e7_pole`; the resulting endpoint-compatible ambient
has dimension 54:

```bash
sage -python elkies-k3/scripts/assemble_h92_q8_endpoint_rr_ambient.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json
```

Expected terminal status:

```text
H92Q8ENDPOINTAMBIENT|families=18|basis=54|h_degree=4|status=PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT
```

This is a bounded seed for resolved-condition compilation, not a completed
global Riemann--Roch space: it still needs the smooth P1.O collision block
and the remaining actual E7 gluing conditions.

The actual smooth block rules out this endpoint seed by itself.  Its
1080-by-54 principal-part map has full column rank modulo `43` (which already
certifies characteristic-zero injectivity), so no nonzero element of the
54-term seed is smooth at all four collisions.  This is a bounded obstruction
to the seed, not to the q8 pencil: the endpoint construction may require a
larger global envelope.

```bash
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43.json
```

<!-- status-consumer: EC-K3-H3-Q8-ENDPOINT-SMOOTH-OBSTRUCTION d629df5d215d009c -->

For the same nested enlargement, replace every family denominator power `k`
by `k+r` and include all `e8_floor <= i <= 4*(k+r)+e7_pole`.  With the
corrected marked frame, `r=4` has 342 columns and smooth rank 342 modulo
`43`, hence no modular kernel.  The earlier 335-column, three-dimensional
kernel used the invalid `m/t` normalization and is withdrawn.  This is a
bounded modular obstruction, not a characteristic-zero pencil.

```bash
sage -python elkies-k3/scripts/assemble_h92_q8_endpoint_rr_ambient.sage \
  --extra-h-power 4 \
  --output artifacts/generated-results/elkies-k3-h92-q8-extra4-endpoint-rr-ambient.json
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --extra-h-power 4 --include-kernel \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra4.json
```

<!-- status-consumer: EC-K3-H3-Q8-ENLARGED-ENDPOINT-SMOOTH-KERNEL 3a89050b29ce8eac -->

For this corrected `r=4`, zero-E7-slack ambient, every one of the 342
generators satisfies the exact actual marked-E7 inequality.  Hence the
marked chart adds no row; only the other five E7 edges remain to be compiled.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_enlarged_endpoint_marked_e7.sage \
  --extra-h-power 4 --extra-e7-pole 0 \
  --output artifacts/generated-results/elkies-k3-h92-q8-extra4-marked-e7-cover.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA4-MARKED-E7-COVER 5deb19aa922fd23b -->

The same enlargement automatically retains the E8 cover: `h(0)=1` and the
certified lower `u`-floor for every generic monomial is unchanged.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_enlarged_endpoint_e8_cover.sage \
  --extra-h-power 4 --output artifacts/generated-results/elkies-k3-h92-q8-extra4-e8-cover.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA4-E8-COVER 135e8e4da32d56b6 -->

The smooth collision contribution is now bounded exactly in the q6 saturated
coordinate `q=(m-y(P1)/x(P1))/h`.  Substituting `m=r/h+h*q` into all 54 seed
generators produces no principal part beyond `h^-15`, so this records the
finite jet range for the future smooth condition matrix:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_smooth_collision_principal_parts.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-collision-principal-parts.json
```

Expected terminal status:

```text
H92Q8SMOOTHPROFILE|ambient=54|max_h_pole=15|q_frame=actual_q6_saturated|status=PASS_EXACT_Q8_SMOOTH_COLLISION_PRINCIPAL_PARTS
```

The resulting smooth condition block is now compiled exactly through the
shared finite-ambient-image interface, rather than only screened modulo a
prime. It imposes every negative `h`-principal-part
coefficient in the actual `(q,X)` frame.  The `1080 x 54` template has full
column rank modulo `43`; because reduction cannot increase rank, this proves
that its characteristic-zero rank is 54 and that the least endpoint seed has
no smooth-compatible direction.  This obstructs this seed locally only; a
larger endpoint ambient and the remaining resolved E7 conditions are still
needed for the q8 hop.

```bash
sage -python elkies-k3/scripts/compile_h92_q8_smooth_principal_parts_exact.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-exact.json
```

Expected terminal status:

```text
H92Q8SMOOTHEXACT|extra_h=0|rows=1080|columns=54|rank=54|kernel=0|status=PASS_EXACT_Q8_SMOOTH_PRINCIPAL_PART_CONDITION_BLOCK
```

<!-- status-consumer: EC-K3-H3-Q8-SMOOTH-PRINCIPAL-PART-CONDITION-BLOCK 6888f78205d642be -->

On its own this profile does not prescribe the allowed q8 submodule or assert
a condition rank; the following line-bundle calculation supplies that local
identification. No q8 pencil is asserted.

The exact collision algebra itself can now be regularized before specifying
that submodule.  With `h=Z4`, set `q=(m-y(P1)/x(P1))/h` and `X=h^2*x`.
Clearing the chord equation and removing its marked `x=x(P1)` factor gives a
monic quadratic in `X` whose coefficients are regular at every smooth
`h=0` collision:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_smooth_collision_frame.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-collision-frame.json
```

Expected terminal status:

```text
H92Q8SMOOTHFRAME|quadratic_in_X=1|frame=18|h_regular=1|status=PASS_EXACT_Q8_SMOOTH_COLLISION_FRAME
```

The source-nef divisor has no vertical modification at these smooth fibres:
its E7/E8 vertical terms are supported at the additive fibres, while its
`-11F` term can be represented away from a chosen collision.  Thus the
regular `q,X` frame is the actual q8 line-bundle lattice, and all negative
`h` principal parts in it must vanish.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_smooth_line_bundle_lattice.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-line-bundle-lattice.json
```

<!-- status-consumer: EC-K3-H3-Q8-SMOOTH-LINE-BUNDLE-LATTICE 340ac8bee0e38750 -->

The matching source-E8 target is:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e8_local_target.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e8-local-target.json
```

Expected terminal status:

```text
H92Q8E8TARGET|degrees=1,0,0,0,0,0,0,0|twist=-4,-5,-7,-10,-8,-6,-4,-2|status=PASS_EXACT_Q8_E8_SOURCE_TARGET
```

It identifies the cycle in both source and actual blow-up-chart component
orders; its finite chart-level quotient map remains to be derived.
<!-- status-consumer: EC-K3-H3-Q8-E8-TARGET 7f4bac8ed72db930 -->

The complete E8 local module is then derived by:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e8_complete_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e8-complete-module.json
```

Expected terminal status:

```text
H92Q8E8MODULE|ideal=(u2,X,Y)|colength=2|status=PASS_EXACT_Q8_E8_COMPLETE_MODULE
```

<!-- status-consumer: EC-K3-H3-Q8-E8-MODULE 74327bc7489c8ca6 -->

The generic q8 basis now has exact, per-generator E8 coefficient floors in
the same resolved chart.  For `m^b` the floor is `u^(11+2b)`; for `x*m^b` it
is `u^(13+2b)`.  These follow from `m=u^-2*Q` with `Q` a unit and the actual
module `u^9*(u^2,X,Y)`:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e8_ambient_weights.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e8-ambient-weights.json
```

Expected terminal status:

```text
H92Q8E8WEIGHTS|basis=18|m_floors=11..29|xm_floors=13..27|status=PASS_EXACT_Q8_E8_AMBIENT_WEIGHTS
```

These are only the E8 end of a q8 coefficient ambient; the non-anti-nef E7
gluing and the remaining finite conditions still have to be imposed before a
global kernel or a child equation is claimed.

For a bounded modular reconnaissance on the exact `E8+E6` q6 child, the
following first checks that the selected prime retains the additive-fibre
valuations, then exhausts the polynomial `x` ansatz on an explicitly chosen
`IV*` branch:

```bash
sage -python elkies-k3/scripts/search_h92_q6_child_polynomial_sections_modp.sage \
  --prime 43 --max-x-degree 4 --require-iv-star-singular \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-mod-43-iv-singular.json
```

Expected terminal status:

```text
H92Q6CHILDPOLYMOD|prime=43|x_degree=4|iv_singular=1|x_space=79507|sections=6|status=EXPERIMENTAL_EXHAUSTIVE_MODULAR_ANSATZ
```

This is a finite-field, bounded-ansatz experiment; it does not prove a
characteristic-zero nonexistence statement or produce the q8 pencil.

The six mod-43 residues have coefficient-Jacobian rank 11, so they must not
be treated as unique p-adic lifts.  The diagnostic is:

```bash
sage -python elkies-k3/scripts/lift_h92_q6_child_polynomial_sections.sage \
  --input artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-mod-43-iv-singular.json \
  --precision 64 \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-hensel-43.json
```

It reports `records=6|isolated=0` and is experimental p-adic evidence only.

The reusable exact neighbor engine packages the supplied-wall reduction,
primitive `U`-split, root/MW minimization, transported component classes, and
versioned certificate separately from the global chamber proof.  Rebuild its
pinned q80 and H3 certificates, then replay the engine regressions with:

```bash
sage elkies-k3/scripts/build_exact_neighbor_engine_certificates.sage
sage elkies-k3/scripts/verify_exact_neighbor_engine.sage
```

The certificate builder refuses to replace changed output.  Its q80 and H3
artifacts have payload hashes
`b6ea4c8b421cf782bf57416935b20bb3424118c0531a236c4e66548bc07895c3` and
`093a0e1b7fe8a1ef93cfffaa758762ae7e7ff83278ee74d439e1ff4ea052c01c`.
The full `D13/MW4` to rootless/MW17 replay below uses this same engine for its
eleven primitive splits while retaining the pinned root-adapted chain bases.

The deterministic lateral q=4 presentation from this D13 frame is certified
with:

```bash
sage -python elkies-k3/scripts/analyze_h3_d13_q4_chamber.sage
```

It proves, using an exact D13/MW closest-vector calculation rather than a
bounded section scan, that the raw `(a,b)=(2,2)` class is already nef and has
old-fiber degree two.  Its child is `A12+A1/MW4`, with root data
`(13,158,26)`.  This is a genuine but lateral presentation, not a rank gain.
The same checker also certifies the preferred q=24 rank-growing class below.

The D13 proper-presentation no-growth barrier now extends through q=20.  The
q=20 endpoint is reproduced in both factor orders by:

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt \
  --root-rank 13 --q 20 --degree 2 --rank-growth-only \
  --output artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q20-degree2.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt \
  --root-rank 13 --q 20 --degree 4 --rank-growth-only \
  --output artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q20-degree4.json
```

Each run has 1,567 dominant Weyl orbits.  Degree two has 1,533 primitive
neighbors and degree four has 1,567; neither has MW rank above four.  Together
with the previously completed proper-presentation quotient shells through
q=18, this closes all proper D13 presentations through q=20.  This is an
exact bounded barrier, not a global obstruction; continuation should use a
Weyl quotient rather than a larger blind shell.

The intervening q=21 degree-three and q=22 degree-two proper presentations
also have maximum MW rank four; q=23 has no proper factor presentation.  The
first rank growth is q=24, already with old-fiber degree two:

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt \
  --root-rank 13 --q 24 --degree 2 --rank-growth-only --adapt-mw-at-least 5 \
  --frames-dir artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2-frames \
  --output artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2.json
```

The exact quotient has 2,709 dominant orbits and 2,653 primitive neighbors.
Exactly three have root data `(12,264,4)`, identifying `D12/MW5`.  The
preferred orbit 85 is component-nef without reflections, and the chamber
checker proves full nefness from an empty shifted MW ball of radius squared
two plus the degree-two bisection parity identity.  The q=24 artifact SHA-256
is `66d5a7ff6ec26f8aa8344cdbd779a6c96707b041ba4f89d7dbfe460c95485a93`.
Thus the exact lattice/chamber prefix is now
`MW2 --q6--> MW3 --q8--> MW4 --q24--> MW5`, with every geometric step of
old-fiber degree two.

The equation-side q24 horizontal section is also now exact.  It is recovered
from the modular degree-46 q24 section by the structured 24-by-24 Newton lift,
not by the failed rational trace interpolation path:

```bash
sage -python elkies-k3/scripts/lift_q24_structured_newton.sage --precision 8192
```

Expected terminal lines include:

```text
Q24STRUCT_RECON|resolved=156/156|complete=1|modulus_bits=136067|status=PASS
Q24STRUCT_RESULT|identity=PASS|modp=PASS|...|status=PASS_EXACT_Q24
```

The p-adic checkpoint is
`artifacts/local/elkies-k3/q24-structured-hensel-p8192.json`, with status
`PASS_Q24_STRUCTURED_HENSEL`.  The exact characteristic-zero section is
`artifacts/local/elkies-k3/q8-q24-horizontal-section-qq.json`, with status
`PASS_EXACT_Q24_HORIZONTAL_SECTION`; it has `deg(Z,X,Y)=(24,52,78)` and
rational affine denominators of degrees `(48,72)`.

The D12 equation-level child is recovered from that exact section by the
resolved component-valuation RR compiler:

```bash
sage -python elkies-k3/scripts/lift_q24_d13_to_d12_resolved_rr_qq.sage
```

Expected terminal lines include:

```text
Q24DIVVALQQ|stage=CHILD|...|root_rank=12|root_det=4|euler=24|MW=5|status=PASS_D12
Q24DIVVALQQ|stage=MODULAR_SIGNATURE|...|plane=1|quartic=1|jacobian=1|status=PASS
Q24DIVVALQQ_RESULT|ambient=56|collision=48|post=8|resolved=6|kernel=2|quartic=4|root_rank=12|root_det=4|euler=24|MW=5|status=PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR
```

The output is
`artifacts/local/elkies-k3/q24-d13-to-d12-component-valuation-qq.json`.  This
closes the selected `D13/MW4 --q24 orbit85--> D12/MW5` equation hop.  The next
equation gate is the fixed `D12/MW5 --q6 orbit42--> A11/MW6` hop, with
correction `3`, `P.O=3`, fibre twist `0`, and no zero-pole/dx4 search.

The exact orbit42 prerequisites are reproduced by:

```bash
sage -python elkies-k3/scripts/preflight_h92_q24_orbit42_component_valuation_qq.sage
sage -python elkies-k3/scripts/map_h92_q24_orbit42_i8star_physical_components_qq.sage
```

Expected terminal statuses are `PASS_Q42_DIVVAL_PREFLIGHT` and
`PASS_Q42_EXACT_I8STAR_PHYSICAL_MARKING`.  The physical marking has twelve
components and exactly two surviving spinor-arm orientations; the marked
section meets `C11` in one and `C10` in the other.

The exact resolved-RR equation edge and its equation-side child marking are
reproduced by:

```bash
sage -python elkies-k3/scripts/lift_h92_q24_orbit42_resolved_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q24_orbit42_a11_equation_marking.sage
```

Expected terminal statuses are
`PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR` and
`PASS_Q42_A11_EQUATION_MARKING_ORBIT64_MOD100003`.  The exact RR dimension
chain is `9 -> 3 -> 2`, using the weighted C01 valuation `(2,2,3)`.  The
degree-four quartic has a minimized Jacobian with `I12 + 12 I1`, root lattice
A11, Euler number 24, and MW rank 6 under the rank-19 marking.  The ordered
identity-shell new-fibre degrees select orbit64/mapping7 in the chosen C10
orientation; orbit65/mapping6 is its spinor conjugate.  This closes the
selected `D12/MW5 --q6 orbit42--> A11/MW6` equation edge.  The next equation
gate is `A11/MW6 --q8--> 2A5/MW7` in the orbit64 child frame.

<!-- status-consumer: EC-K3-H3-Q24-O42-QQ-A11 ffa4308117c55056 -->

The construction-compatible equation-side q8 target and the reduced modular
section chart are reproduced by:

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/local/elkies-k3/q24-downstream-lift/d12-c10a-zero-q6-frames/q6-o0064-r11-n132-d12-ad4a027cb197.txt \
  --root-rank 11 --q 8 --degree 2 --rank-growth-only --adapt-mw-at-least 7 \
  --frames-dir artifacts/local/elkies-k3/q24-a11-orbit64-q8-frames \
  --output artifacts/local/elkies-k3/q24-a11-orbit64-q8-all.json
sage -python elkies-k3/scripts/certify_h92_q24_a11_q8_construction_fingerprint.sage
sage -python elkies-k3/scripts/build_h92_q24_a11_q8_component6_section_modp.sage \
  --prime 43 --branch 0 --eliminate-r-count 8
sage -python elkies-k3/scripts/certify_h92_q24_a11_target_coset_bridge.sage
sage -python elkies-k3/scripts/build_h92_q24_a11_q8_component6_section_modp.sage \
  --prime 100003 --free-infinity --pole-order 5 --component-depth 3 \
  --eliminate-r-count 0
```

The exact construction certificate reports two root-chain orientations, an
MW automorphism group of order 16, eight integral-glue isometries, and four
nef transports in exactly equation-side orbits `12` and `2162`.  Orbit `12`
is the declared minimum-MW-L1 target and retains the historical formula
`O+P-2F` with zero vertical-root correction.  The modular builder reports a
26-variable, 28-equation chart with 16 infinity branches.  These are a partial
target/marking result and a discovery system; they do not prove the section or
q8 equation edge.

<!-- status-consumer: EC-K3-H3-A11-Q8-CONSTRUCTION-TARGET c892eec88af45f08 -->

The target-coset certificate additionally proves that the eighteen exact
identity-shell points generate an index-five sublattice, that the old
pole-order-four missing-direction vector is in the wrong coset, and that the
minimum-pole bridge is `M=(1,0,0,0,0,1)` with `P.O=5`, height `47/4` and
correction `9/4`.  It verifies the exact equation-marked word
`P12=M+S6-2*S2-2*S8`.  The free-infinity command emits a 36-variable,
37-equation system covering all leading smooth-fibre points at once.  The
recorded 600-second msolve benchmark used `-t 4 -l 42 -v 2`, completed two
degree-eight reductions and stopped in a third `257692 x 2857438` matrix; no
output section or non-existence conclusion was obtained.

<!-- status-consumer: EC-K3-H3-A11-Q8-TARGET-COSET-BRIDGE 8d17ab150a7e3567 -->

The later exact A11 q8 closeout and the first three physical-suffix equations
are replayed in construction order by:

```bash
sage -python elkies-k3/scripts/lift_h92_q24_a11_q8_residual_resolved_hensel.sage
sage -python elkies-k3/scripts/derive_h92_q24_a11_q8_difference_qq.sage
sage -python elkies-k3/scripts/lift_h92_q24_a11_q8_resolved_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q24_a11_q8_equation_marking_qq.sage

sage -python elkies-k3/scripts/certify_h92_a5a5_physical_q4o208_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_a5a5_physical_q4o208_equation_marking_qq.sage

sage -python elkies-k3/scripts/certify_h92_q4o208_physical_q4o1584_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o1584_equation_marking_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o1584_physical_q4o164_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o164_c8_equation_marking_qq.sage

# Fixed-ADE replacement audit for the invalid q4/orbit323 identification:
sage -python elkies-k3/scripts/certify_h92_q4o323_horizontal_marking_qq.sage
sage -python elkies-k3/scripts/lift_h92_q4o323_horizontal_via_summands_qq.sage --precision 120
sage -python elkies-k3/scripts/compile_h92_q4o208_q4o1599_a3_2a2_qq.sage

# Corrected fixed-corridor q4/orbit323 recovery and exact equation compiler:
sage elkies-k3/scripts/construct_h92_q4o323_horizontal_by_halving_qq.sage
sage elkies-k3/scripts/compile_h92_q4o208_q4o323_a3_2a2_qq.sage

# Carry the physical q323 wall through the fixed suffix and certify the
# wall-corrected q207 continuation (q12, old degree two, child 5A1/MW12):
sage elkies-k3/scripts/build_h92_q4o323_reflected_fixed_suffix_marking.sage
sage -python elkies-k3/scripts/certify_h92_q4o323_component2_pointing_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_degree2_candidate.sage \
  --source-marking artifacts/local/elkies-k3/q4o323-reflected-fixed-suffix-component2-marking.json \
  --source-frame artifacts/generated-results/elkies-k3-h3-q4o208-corrected-a3-2a2-old_a11_component_2-frame.txt \
  --fibre 6,2,-2,3,3,4,3,9,0,3,0,4,2,-4,2,-1,1,1,-1 \
  --candidate-label physical_q12_reflected_q4o207_after_q4o323 \
  --target current_4A1 \
  --frame-output artifacts/local/elkies-k3/q4o323-physical-q12-reflected-q4o207-5a1-frame.txt \
  --output artifacts/local/elkies-k3/q4o323-physical-q12-reflected-q4o207-5a1-certificate.json

# Complete polynomial P.O=0 shell on the component-2-pointed q323 child.
# The script rejects p=31 because both I3 fibres specialize to I4; p=61
# preserves the marked fibre multiplicities.
sage -python elkies-k3/scripts/construct_h92_q4o323_p0_shell_modp.sage --prime 61
```

The component-2 pointing replay selects `W=+L0(u)` in the stored positive
square-root normalization, verifies `81*A_pointed=A_child` and
`729*B_pointed=B_child`, and constructs the opposite branch exactly over
`QQ(u)`.  The good-prime shell has 602 signed polynomial sections; 120 have
ordinary coefficient-Jacobian rank 12.  This is a modular construction
frontier, not yet the q12 horizontal or a characteristic-zero lift.

The exact q8/orbit376 horizontal is reconstructed without a large Groebner
calculation.  The following good-prime set gives the pinned 566-bit CRT
certificate; every trace uses 91 training fibres and nine holdouts:

```bash
for p in 131 137 151 157 167 173 181 \
  1000003 1000033 1000037 \
  1000000007 1000000009 1000000021 \
  2000000011 2000000033 2000000063 \
  2000000087 2000000089 2000000099 \
  1000000000000000003 1000000000000000009 1000000000000000031; do
  sage -python elkies-k3/scripts/probe_h92_q4o164_inherited_p1_abel_trace_modp.sage \
    --prime "$p" --interpolate --good-fibre-limit 100
  sage -python elkies-k3/scripts/identify_h92_q4o164_q8_horizontal_mod131.sage \
    --prime "$p"
done
sage -python elkies-k3/scripts/reconstruct_h92_q4o164_q8_horizontal_crt_qq.sage
```

The final command reports
`PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT`, exact compact degrees
`(12,8)/(18,12)`, literal `QQ(t)` substitution, and replay at all 22 primes.

The first four commands close `A11/MW6 --q8/orbit12--> 2A5/MW7` with an
exact `14 -> 2` resolved plane and a `2I6+12I1` Jacobian.  The q4/orbit208
pair closes the physical `3A3/MW8` equation and its effective `C5` marking.
The q4/orbit1584 pair gives `D4+A3+3A1/MW7`, with finite
`I4+3I2+8I1`, `I0*` at infinity, and an exact second-I6-affine zero.  The
q4/orbit164 pair gives `2A3+2A1/MW9`, with finite `I4+2I2+12I1`, `I4` at
infinity, and exact `C8` pointing.  Each fibre list has Euler number 24.
The replacement audit proves that branch 33 inverse is not q4/orbit323, then
uses lifted branch 16 instead.  Its q4/orbit1599 `5 -> 3 -> 2` RR plane has
`I4+2I3+14I1`, hence `A3+2A2/MW10`.  Both complete signed graph solutions
give branch 16 the same sign and NS class.  The raw-to-physical transport then
checks `D=O+P+V` with zero MW tail in `V`, so this equation edge is aligned
unambiguously with its own marked q4/orbit1599 lattice edge; only an unused
shell branch remains twofold.  That marked child is not the stored canonical
q4/orbit323 frame which feeds q4/orbit207.

The corrected fixed-corridor replay uses a different construction.  In the
marked Mordell--Weil quotient it verifies
`2*T=P8+2*P18+P33-2*C7`; the two inherited global orientations give rational
halves with `P.O=3` and `P.O=1`, so the latter is selected exactly.  The
duplication quartic factors as degrees `1+3`, and the marked half has compact
degrees `x=(6,2)`, `y=(8,3)`.  The second command verifies
`D=O+T+C2+C3` on the second old `I4`, computes the resolved `5 -> 3 -> 2`
plane, and returns a minimal `I4+2I3+14I1` Jacobian.  Both commands avoid a
Groebner basis; the first also avoids a new finite-field shell search.
<!-- status-consumer: EC-K3-H3-Q4O208-Q4O323-QQ-A3-2A2 a903147a9023d49f -->

The resulting exact marked lattice suffix continues as

```text
3A3/MW8 --q4/orbit1584--> D4+A3+3A1/MW7
          --q4/orbit164--> 2A3+2A1/MW9
          --q8/orbit376--> 4A1/MW13
          --q12/orbit5867--> rootless/MW17.
```

The first three arrows have complete characteristic-zero equations. For q8,
the inherited-`P1` unsplit degree-seven Abel strategy gives the exact
horizontal, the resolved `12 -> 4 -> 2` pencil, quartic, and `4A1/MW13`
Jacobian. The finite `T=0` I4 chain is oriented exactly, and a rational arc on
the first finite-I2 exceptional conic points the child directly at `P1229`,
with quartic sign `-1` and exact `81/729` invariant identities. The
q12/orbit5867 edge, endpoint isometry, and characteristic-zero rootless
equation are exact. q12/orbit5867 is the preferred optional final edge. Its
nominal lattice compiler word has four physical `P.O=0` parent branches of
degrees `(3,2,1,2)` and parent `a-b` values `(2,2,1,1)`, but complete shells
at four good primes did not realize its Q1 branch. The executed exact word has
degrees `(4,2,1,5)` and `a-b=(4,2,1,4)`; q12/orbit4484 remains the certified
fallback. The direct 17-section R17 and endpoint certificates are now exact.
The full proof and endpoint requirements are in
[`elkies-k3/PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md`](elkies-k3/PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md).

The historical parent-degree-one q12/orbit5867 compiler branch is exact over QQ. Its
marked word is `2*C8opp+B1+2*B2-2*B3+B4-2*B5+B6+B7`; resolved-q8 restriction
and P1229-pointed quartic transport give the mod-131 `(4,6)` section, and a
regular rank-12 Hensel lift reconstructs 800-bit exact coefficients:

```bash
sage -python elkies-k3/scripts/construct_h92_q12o5867_degree1_branch_mod131.sage
sage -python elkies-k3/scripts/lift_h92_q12o5867_degree1_branch_qq.sage
```

The equation-effective word and terminal rootless equation replay with:

```bash
sage -python elkies-k3/scripts/classify_h92_q12o5867_p0_shell_against_lattice_mod89.sage
sage -python elkies-k3/scripts/lift_h92_q12o5867_replacement_word_seeds_qq.sage
sage -python elkies-k3/scripts/construct_h92_q12o5867_target_horizontal_qq.sage
sage -python elkies-k3/scripts/compile_h92_q12o5867_smooth_rr_qq.sage
```

The last command proves the exact `22 -> 2` RR plane and a minimal
degree-`(8,12,24)` Jacobian with geometrically `24I1`. No new parent RR
calculation, elimination, or Groebner basis is used.

The direct rank-17 endpoint basis replays with:

```bash
sage -python elkies-k3/scripts/construct_h92_q12o5867_rootless_p0_shell_mod131.sage
python3 elkies-k3/scripts/select_h92_q12o5867_rootless_mod131_basis.py
sage -python elkies-k3/scripts/lift_h92_q12o5867_rootless_selected_basis_qq.sage
sage -python elkies-k3/scripts/certify_h92_q12o5867_rootless_height_basis_qq.sage
```

This enumerates 2,622 signed regular modular sections, lifts only a selected
17-section basis, and proves its exact determinant-948 QQ height Gram and
determinant-minus-one isometry to pinned R17. Thus rank at least 17 and trivial
torsion are unconditional. The final endpoint replay is:

```bash
sage -python elkies-k3/scripts/certify_h92_q12o5867_endpoint_qq.sage
```

It verifies an exact rational point on the q12 binary quartic and the
`81/729` pointed-invariant identities, counts the K3 over `F_p` and
`F_{p^2}` for `p=131,137`, applies the two-prime discriminant-square-class
test to prove geometric Picard rank 19, and excludes the unique possible
index-two even overlattice. Its terminal status is
`PASS_EXACT_Q12O5867_SOURCE_IDENTITY_RHO19_FULL_MW_R17`; hence the full
geometric Mordell--Weil group is saturated R17 of exact rank 17, trivial
torsion, and determinant 948. No Groebner basis or surface elimination is
used.

<!-- status-consumer: EC-K3-H3-Q4O208-Q4O1584-QQ 8463d62d0e9f2b83 -->
<!-- status-consumer: EC-K3-H3-Q4O208-Q4O1599-QQ-A3-2A2 2018f08fd6b8e2a9 -->
<!-- status-consumer: EC-K3-H3-Q4O1584-Q4O164-QQ bafe854a24d6762b -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-HORIZONTAL 688f0a5f6d989e9c -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-SMOOTH-RR c0bc752848961743 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-DEGREE1-SECTION-QQ 28056772646e9fc7 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-QQ-ROOTLESS 6a3dbee5942ddd0a -->
<!-- status-consumer: EC-K3-H3-Q12O5867-QQ-R17-BASIS a2097150acf00645 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-ENDPOINT-QQ a83b08acd921c32b -->
<!-- status-consumer: EC-K3-H3-Q4O208-R17-CURRENT-MARKED-ROUTE 432b34c44c78bcb9 -->

The retained zero-pole boundary audit is reproduced by:

```bash
sage -python elkies-k3/scripts/recover_h92_q24_orbit42_zero_pole_smallprime.sage --prime 43
sage -python elkies-k3/scripts/scan_h92_q24_orbit42_zero_pole_model_modp.sage --prime 53
sage -python elkies-k3/scripts/lift_h92_q24_orbit42_zero_pole_sections_qq.sage --precision 65536
sage -python elkies-k3/scripts/lift_h92_q24_orbit42_spinor_zero_pole_sections_qq.sage
sage -python elkies-k3/scripts/construct_h92_q24_orbit42_exact_section_candidates_qq.sage
```

The two final terminal statuses are
`PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ`, with nine signed
identity-class pairs, and `PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ`, with
the remaining opposite spinor-class pair.  All twenty zero-pole sections now
have exact characteristic-zero Weierstrass identities.  This supplied a
construction aid for the subsequently completed A11 child; it was not itself
the resolved-RR proof.  The two additional shortcut audits are reproduced by:

```bash
python3 elkies-k3/scripts/run_h92_q24_orbit42_fast_parallel.py
sage -python elkies-k3/scripts/analyze_h92_q24_orbit42_identity_halving.sage
sage -python elkies-k3/scripts/recover_h92_q24_orbit42_by_identity_halving_qq.sage
```

The first command terminates with `STOPPED_INVALID_FAST_Q6_ROUTE` after an
exact passing audit: the converted q6 degrees are `435` and `703`, so no q6
rational-point transport is run.  The identity-halving lattice gate is exact;
the following equation-shell matching and chord census is modulo `100003`.
Its terminal status `Q42_IDENTITY_HALVING_HAS_NO_A11_CHORD` records four
rational degree-three candidates whose branch polynomials remain squarefree
of degree `18`.  That modular result rejects this shortcut only and is not a
characteristic-zero non-existence theorem.

The working route and exact artifact hashes can be checked without rerunning
the long Sage calculations via:

```bash
python3 elkies-k3/scripts/success-path/verify_ledger.py
```

The same root-adapted quotient search continues through three further
rank-growing degree-two steps.  The selected path is reproduced by:

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2-frames/q24-o0085-r12-n264-d4-add367fba084.txt \
  --root-rank 12 --q 6 --degree 2 --rank-growth-only --adapt-mw-at-least 6 \
  --frames-dir artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2-frames \
  --output artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2-frames/q6-o0042-r11-n132-d12-e7e61e5dd4c2.txt \
  --root-rank 11 --q 8 --degree 2 --rank-growth-only --adapt-mw-at-least 7 \
  --frames-dir artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2-frames \
  --output artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2-frames/q8-o0922-r10-n60-d36-c9cd5a498117.txt \
  --root-rank 10 --q 4 --degree 2 --rank-growth-only --adapt-mw-at-least 8 \
  --frames-dir artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2-frames \
  --output artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2.json
sage -python elkies-k3/scripts/analyze_h3_rank_growing_degree2_chain.sage
```

The chamber checker proves full nefness for the selected q6, q8, and q4
classes.  Their children are respectively `A11/MW6`, `A5+A5/MW7`, and
`3A3/MW8`, with root data `(11,132,12)`, `(10,60,36)`, and `(9,36,64)`.
The authoritative search-artifact hashes are, in the same order,
`1b8d7f37794bcf49c48949cec3bffe7baba69cb55e369f7af9f5002908a75b7f`,
`d336b0d32a07907ec61464c7d5ace4c76f257ddcf16b904a96a3d9064f408323`,
and `98fdd553768b27d5800f247b41a6e2a28f0ee2787ad5959d96ef420a2eb09185`.
The selected `3A3/MW8` continuation frame has SHA-256
`e535e5abc8c70c79be9b088c1217a307c7c4940de1acd4b2cad4cb2b9fda22bb`.

The exact low-degree lattice chain is therefore

```text
H3 MW2 -q6-> MW3 -q8-> MW4 -q24-> D12/MW5
       -q6-> A11/MW6 -q8-> A5+A5/MW7 -q4-> 3A3/MW8,
```

and every arrow has old-fiber degree two.

One further selected q=4 shell from the `3A3/MW8` endpoint is reproduced by:

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2-frames/q4-o0472-r9-n36-d64-4841c34fa442.txt \
  --root-rank 9 --q 4 --degree 2 --rank-growth-only --adapt-mw-at-least 10 \
  --frames-dir artifacts/generated-results/elkies-k3-h3-a3x3-q4-degree2-frames \
  --output artifacts/generated-results/elkies-k3-h3-a3x3-q4-degree2.json
sage -python elkies-k3/scripts/analyze_h3_rank_growing_degree2_chain.sage
```

The quotient has 2,481 dominant primitive orbits.  Selected orbit 323 has
child root data `(7,24,36)`, identifying `A3+2A2/MW10`.  The search artifact
SHA-256 is
`8f1d5105831cc3356bc4598932380295bcc6e91629e0d05a6a9d64c0d840d29d`.
The final q=4 presentation is also nef of old-fiber degree two.  Its source
component pairings are `(1,0,1,0,0,0,1,1,0)`, its affine pairings are
`(0,1,1)`, and exact shifted root/MW CVP gives section distances
`2,2,3` repeated sixteen times, then `4,4`.  Bisection parity closes the last
case.  Thus the geometric chamber certificate currently reaches
`A3+2A2/MW10` with every arrow a degree-two pencil.

The exact Weyl-quotient lattice continuation reaches a rootless frame:

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-a3x3-q4-degree2-frames/q4-o0323-r7-n24-d36-87b284dff2bc.txt \
  --root-rank 7 --q 4 --degree 2 --rank-growth-only --adapt-mw-at-least 11 \
  --frames-dir artifacts/generated-results/elkies-k3-h3-mw10-a3a2a2-q4-degree2-frames \
  --output artifacts/generated-results/elkies-k3-h3-mw10-a3a2a2-q4-degree2.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-mw10-a3a2a2-q4-degree2-frames/q4-o0207-r5-n10-d32-a462a553a1e9.txt \
  --root-rank 5 --q 4 --degree 2 --rank-growth-only --adapt-mw-at-least 13 \
  --stop-after-first-growth \
  --frames-dir artifacts/generated-results/elkies-k3-h3-mw12-5a1-q4-degree2-first-hit-frames \
  --output artifacts/generated-results/elkies-k3-h3-mw12-5a1-q4-degree2-first-hit.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-mw12-5a1-q4-degree2-first-hit-frames/q4-o0052-r4-n8-d16-066f47d7fff3.txt \
  --root-rank 4 --q 4 --degree 2 --rank-growth-only --adapt-mw-at-least 14 \
  --stop-after-first-growth \
  --frames-dir artifacts/generated-results/elkies-k3-h3-mw13-4a1-q4-degree2-first-hit-frames \
  --output artifacts/generated-results/elkies-k3-h3-mw13-4a1-q4-degree2-first-hit.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-mw13-4a1-q4-degree2-first-hit-frames/q4-o0114-r3-n6-d8-018b225c409b.txt \
  --root-rank 3 --q 4 --degree 2 --rank-growth-only --adapt-mw-at-least 15 \
  --stop-after-first-growth \
  --frames-dir artifacts/generated-results/elkies-k3-h3-mw14-3a1-q4-degree2-first-hit-frames \
  --output artifacts/generated-results/elkies-k3-h3-mw14-3a1-q4-degree2-first-hit.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-mw14-3a1-q4-degree2-first-hit-frames/q4-o0498-r2-n4-d4-e86cf7c2d2f9.txt \
  --root-rank 2 --q 4 --degree 2 --rank-growth-only --adapt-mw-at-least 16 \
  --stop-after-first-growth \
  --frames-dir artifacts/generated-results/elkies-k3-h3-mw15-2a1-q4-degree2-first-hit-frames \
  --output artifacts/generated-results/elkies-k3-h3-mw15-2a1-q4-degree2-first-hit.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-mw15-2a1-q4-degree2-first-hit-frames/q4-o0981-r1-n2-d2-4f02793cfc09.txt \
  --root-rank 1 --q 6 --degree 2 --rank-growth-only --adapt-mw-at-least 17 \
  --stop-after-first-growth --stream-first-growth --stream-skip 1000 \
  --stream-limit 2000 --stream-progress-every 500 --mw-vector-cap 10000 \
  --mw-vectors-cache artifacts/generated-results/elkies-k3-h3-mw16-a1-q6-mw-vectors-cap10000.json \
  --frames-dir artifacts/generated-results/elkies-k3-h3-mw16-a1-q6-degree2-stream-first-hit-frames \
  --output artifacts/generated-results/elkies-k3-h3-mw16-a1-q6-degree2-cap10000-stream-chunk001.json
sage -python elkies-k3/scripts/verify_h3_d13_to_mw17_path.sage
sage -python elkies-k3/scripts/analyze_h3_mw10_to_rootless_chambers.sage
sage -python elkies-k3/scripts/verify_rank17_to_h3_reverse_transport.sage
```

The selected suffix is

```text
A3+2A2/MW10 -q4-> 5A1/MW12 -q4-> 4A1/MW13
             -q4-> 3A1/MW14 -q4-> 2A1/MW15
             -q4-> A1/MW16 -q6-> rootless/MW17.
```

All six factor presentations have old-fiber degree two.  The q=6 hit is
streamed witness 2,247 in an explicitly capped 10,000-vector MW quotient
sample; PARI reports 7,187,438 vectors in the full shell, so the search is
bounded but the displayed hit and child are exact.  Before that pivot, q=4
was tested with an exact 9,000-orbit prefix plus 1,000-orbit strata starting
at 16,000, 32,000, 48,000, 64,000, 80,000, 112,000, and 144,000; none was
rootless.  This is a bounded obstruction for the selected MW16 marking, not
an exhaustive q=4 theorem.

The eleven-step replay has determinant-one composite transport with SHA-256
`6542f74b2780b4143999e346d519bb72690fac2eeeb99293a97192f305d24c40`.
Its terminal status is:

```text
H3D13MW17|steps=11|final=rootless|MW=17|...|status=PASS_H3_D13_TO_MW17_LATTICE_PATH
```

The generated replay SHA-256 is
`f6eac2339c86de84b79a0ddfec3229df9b9c1617110bdd9c474443e7e39fd484`.
The second checker proves that all six raw suffix divisors are already in the
chosen chambers, require no reflections, and are nef by exact full-frame CVP
plus bisection parity.  Its terminal status is
`PASS_H3_MW10_TO_ROOTLESS_NEF`.  Thus the entire displayed path is an exact
nef degree-two geometric chain.  Characteristic-zero equation execution
remains a separate gate.

The third checker identifies the unnamed rootless endpoint with the pinned
recovered `rank17_gram.txt` by a determinant-one positive-frame isometry and
then inverts all thirteen H3 corridor transports.  It exports a lossless
fourteen-stage ledger in both H3 and pinned-R17 coordinates, including the
exact NS bridge between the distinct dominant and component-nef q8/D13
markings.  Expected status and pinned artifact hash are

```text
PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT
db9518ee9ba5ffb520898242cbff06894900ea1fea2908476e40433a212af4d2
```

See
[`elkies-k3/RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md`](elkies-k3/RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md)
for the exact proof boundary and retained fields.
<!-- status-consumer: EC-K3-H3-D13-MW17-LATTICE-CHAIN 2c6a2a36699933ab -->

The H3 q8-child finite additive gate also has a characteristic-zero q-frame
calculation; it proves the six finite rows have rank six, leaving only the
global infinity condition:

```bash
sage elkies-k3/scripts/derive_h92_q6_child_q8_finite_q_module_qq.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-FINITE-Q-MODULE-QQ 83452e0f702d6d9a -->

The same q frame has a degree-96 generic vertical pole divisor. A global
coefficient pair must either make `B` vanish there or use matching base
principal parts in `C`; this is an exact global gluing condition, not part of
the finite II*/IV* module.

```bash
sage elkies-k3/scripts/derive_h92_q6_child_q8_q_pole_profile.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-FRAME-POLE-PROFILE 3685b882ed19702a -->

The modular principal-part normalization has degree `95` at primes `43` and
`59` (normalized infinity order one); prime `53` has a leading-coefficient
drop and is excluded for that reconstruction step.

```bash
sage elkies-k3/scripts/probe_h92_q6_child_q8_q_pole_normalization_modp.sage --prime 43
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-FRAME-NORMALIZATION-MODULAR 6bd8f54f7e8887a9 -->

The exact degree-95 correction can be reconstructed without a raw rational
extended Euclidean computation by deterministic modular CRT.  The following
is a resumable *work checkpoint*, not a certificate: it saves the 96 CRT
residues locally and deliberately suppresses expensive interim rational
reconstruction.  A later run may add `--resume-from` with a larger total
`--maximum-primes`; only a run ending `PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION`
has performed the withheld-prime and exact `QQ[T]` checks.

```bash
sage elkies-k3/scripts/reconstruct_h92_q6_child_q8_q_pole_normalization_crt.sage \
  --prime-bits 31 --maximum-primes 4000 --minimum-primes 5000 \
  --reconstruct-every 5000 \
  --checkpoint artifacts/local/h92-q8-q-normalizer-4000.crt.json \
  --output artifacts/local/h92-q8-q-normalizer-4000.json
```

The exact continuation is pinned at 4,600 good 31-bit primes: it reconstructs
the degree-95 correction, passes five incorporated and three withheld modular
checks, then proves its `QQ[T]` congruence exactly.  Replay the certificate:

```bash
sage -python elkies-k3/scripts/reconstruct_h92_q6_child_q8_q_pole_normalization_crt.sage \
  --prime-bits 31 --maximum-primes 4600 --minimum-primes 4001 \
  --reconstruct-every 100 --accepted-validation-primes 5
```

This degree-95 normalizer belongs to the withdrawn degree-46 section and is
retained only as a historical diagnostic. It omits the `Dx` factor required
for the primitive degree-ten q8 section and must not be used in the canonical
q8 construction; the corrected exact checker above proves
`R*h*Dy == Ny*Dx mod Nx` directly.

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-FRAME-NORMALIZATION-CRT 6f6a4e097d4bddd5 -->

The resulting `q_regular=q-R/Nx` has no generic base vertical pole and has
infinity order one.  This is a coordinate-frame certificate, not yet the
transformed q8 local-module assembly:

```bash
sage -python elkies-k3/scripts/certify_h92_q6_child_q8_q_regular_frame.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-FRAME 1cead360048f47f6 -->

Transport the complete finite II*/IV* q8 module to this normalized frame:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_q_regular_finite_module_qq.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-FINITE-MODULE 794f629d93420ed4 -->

The same `q_regular` coordinate also trivializes the complete smooth
degree-46 collision module: the transition identities between
`a=A/h^2,b=B/h` and `C+B*q_regular` are exact in `QQ[T]_(h)`, so there is no
additional `h`-supported quotient row in this frame.

```bash
sage -python elkies-k3/scripts/certify_h92_q6_child_q8_q_regular_smooth_frame.sage
```

Expected terminal status:

```text
H92Q6CHILDQREGSMOOTH|h_degree=46|smooth_quotient=0|status=PASS_EXACT_Q_REGULAR_SMOOTH_FRAME
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-SMOOTH-FRAME b394a23706a914ce -->

The canonical normalized finite-generator ratio, and the bounded monomial
family with `0 <= d,e <= 4`, fail the necessary genus-one branch-degree
screen at every constant level in `GF(43)` and `GF(59)`:

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_q_regular_generator_modp.sage \
  --prime 43 --max-a-monomial-degree 4 --max-b-monomial-degree 4
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_q_regular_generator_modp.sage \
  --prime 59 --max-a-monomial-degree 4 --max-b-monomial-degree 4
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-GENERATOR-OBSTRUCTION 487f98d8a254b683 -->

The unnormalized diagonal finite-module ratio is ruled out by its modular
branch degree (484 rather than 4):

```bash
sage elkies-k3/scripts/probe_h92_q6_child_q8_diagonal_candidate_modp.sage --prime 43 --v 1
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-DIAGONAL-PENCIL-OBSTRUCTION ec76cc0097b76a4c -->

The independent rational Q80 route now has an exact unmarked first `q=4`
equation step. It proves that the reconstructed coefficient curve lies on the
first-child collision divisor, yielding `I5*+I5+8I1` over `QQ(u)`:

```bash
sage elkies-k3/scripts/verify_q80_unmarked_first_q4_collision_qq.sage
```

This does not provide its later pencils, a rootless equation, or bisections.

<!-- status-consumer: EC-K3-Q80-UNMARKED-FIRST-Q4-COLLISION d18185784da1e93d -->

## H3 D12 to retained-Q80 crossover audit

The exact marked transport through the initial Q80 frame and pinned R17 frame
scores all eleven retained Q80 fibre classes in the current equation-side H3
`D12/MW5` coordinates.  It computes exact old-fibre intersections, shortest
D12 section profiles, vertical support, pole bounds, and a connected-layer RR
ambient estimate.  The best degree is still
`16328023738263177`; candidate1 has degree
`370213639961146392704841338`, so the audit rejects the proposed direct
crossover as an equation route.

This machine has the conda-forge SageMath 10.9 environment at the path used
below.  The environment is outside the repository and does not change system
packages:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/score_h3_d12_q80_crossovers.sage
```

The full exact table, CM24 comparison, and claim boundary are in
[`elkies-k3/H3_D12_Q80_CROSSOVER_AUDIT_2026-08-24.md`](elkies-k3/H3_D12_Q80_CROSSOVER_AUDIT_2026-08-24.md).

<!-- status-consumer: EC-K3-H3-D12-Q80-CROSSOVER-AUDIT 34f8f8038e591f00 -->

## Rootless MW17 bisection-orbit enumeration

<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: EC-K3-ELKIES-2026-HIGH-RANK-CALIBRATIONS 345b9fb977057133 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R18-COVER 6b4ee5bbc1afc01e -->
<!-- status-consumer: EC-K3-ELKIES-2026-R19-PAIRED f1e135d2ba803e80 -->
<!-- status-consumer: EC-K3-ELKIES-2026-NAGAO-POSITIVE-CONTROL f99c98cdb6b8cd7d -->
<!-- status-consumer: EC-K3-ELKIES-2026-R28-BAD-PLACE-KUMMER 611e63935d2340bc -->
<!-- status-consumer: EC-K3-ELKIES-2026-R28-S-CLASS-PILOT 8c88abe96881b79d -->
<!-- status-consumer: EC-K3-ELKIES-2026-R28-LOCAL-COVERAGE c078c1aa8e97df47 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R28-PUBLIC-SELMER-CONTROLS 56509673b9eb1940 -->
<!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE 855128c3da8d2b41 -->

The rootless `U + (-M)` lattice has a finite degree-two quotient under section
translation. The following exact lattice calculation enumerates its
section-nonnegative `(-2)` bisection-class orbits and writes their canonical
norm-10 representatives in the pinned `rank17_gram` coordinates:

```bash
sage -python elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage \
  --output artifacts/generated-results/elkies-k3-rootless-bisection-orbits.json \
  --orbits-output artifacts/generated-results/elkies-k3-rootless-bisection-orbits.tsv
```

It is a completed lattice quotient, not a bounded height search. In any
rootless K3 realization, the section-nonnegative classes are consequently
irreducible smooth rational bisections.  The equation stage for the complete
39,120-class quotient is replayed below; this lattice command itself does not
construct branch divisors, quadratic extensions, collision hashes, or a
rank-19 family.  See
[`elkies-k3/BISECTION_COLLISION_SEARCH.md`](elkies-k3/BISECTION_COLLISION_SEARCH.md).
<!-- status-consumer: EC-K3-BISECT-ORBIT 81da2fd80c3623b6 -->

Compute the exact low-intersection priority frontier among the exported
rootless bisection orbits:

```bash
sage -python elkies-k3/scripts/analyze_rootless_bisection_disjoint_frontier.sage \
  --output artifacts/generated-results/elkies-k3-rootless-bisection-disjoint-frontier.json
```

Expected terminal status:

```text
R17BISECTDISJOINT|orbits=39120|norm4_masks=1311|active_masks=1311|pairs=8895801|status=PASS_EXACT_ROOTLESS_BISECTION_DISJOINT_FRONTIER
```

The graph itself ranks equation-level work only.  The complete equation and
squareclass replay below now proves that none of these pairs has an equal
quadratic extension.
<!-- status-consumer: EC-K3-BISECT-DISJOINT-FRONTIER c7ad7497253ac0b3 -->

Rank every surviving pinned orbit by the exact published-section group-law
and chord-input score and compute the complete disjoint-pair graph:

```bash
sage -python elkies-k3/scripts/rank_elkies_2026_bisection_orbits.sage \
  --pool-size 39120 \
  --output artifacts/generated-results/elkies-2026-bisection-equation-priority-full.json \
  --table-output artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv \
  --pairs-output artifacts/generated-results/elkies-2026-bisection-equation-priority-disjoint-pairs-full.tsv
```

Expected terminal status:

```text
ELKIES2026BISECTIONPRIORITY|orbits=39120|pool=39120|disjoint_pairs=8895801|status=PASS_EXACT_R17_BISECTION_EQUATION_PRIORITY
```

Construct all 39,120 exact residual-chord bisections and normalize their
quadratic extensions:

```bash
sage -python elkies-k3/scripts/construct_elkies_2026_bisections.sage \
  --priority-table artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv \
  --limit 39120 \
  --output artifacts/generated-results/elkies-2026-equation-bisections-full.json \
  --orbits-output artifacts/generated-results/elkies-2026-equation-bisections-orbits-full.tsv

.venv/bin/python elkies-k3/scripts/hash_bisection_extensions.py \
  --compact \
  --input artifacts/generated-results/elkies-2026-equation-bisections-full.json \
  --output artifacts/generated-results/elkies-2026-equation-bisection-collisions-full-compact.json
```

The direct compiler solves `M*Nx+Ny=0 mod h^2` for each height-ten trace,
then verifies the exact quadratic relation and the lifted section over its
double cover.  It uses the reciprocal base chart for the unique trace whose
finite denominator loses a pole at infinity.  All 39,120 records pass and all
branch at smooth fibres.  Exact normalization finds 39,120 distinct
squareclasses and no collision.  This exhausts the complete pinned survivor
set, proving injectivity of its bisection-to-extension map and excluding every
one of the 8,895,801 disjoint pairs as a common-cover rank-19 route.  Each
individual smooth cover has one height-12 anti-invariant direction, so the
batch gives 39,120 explicit generic-rank-at-least-18 base changes.

<!-- status-consumer: EC-K3-BISECT-EQUATION-BATCH a0570a5a4ea8e02b -->

Evaluate the complete equation-level bisection atlas at the four rank-25--28
positive controls and at ICARM curve 394 (`t=3/8`):

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-SPECIALIZATION-CONTROLS 04f49e48e1c1dd88 -->

```bash
.venv/bin/python \
  elliptic-curves/scripts/evaluate_elkies_2026_bisections_at_controls.py \
  --check
```

The 195,600 exact square tests find respectively `6,3,2,1,25` split
bisections. Their finite-quotient classes span `5,3,2,1,4` known directions
beyond the generic seventeen and produce no finite-quotient escape. In
particular, the rank-28 control sees only one of its eleven known exceptional
directions, while the `t=3/8` splits span all four directions certifying rank
at least 21. Full-rank relation blocks, verified by exact group addition, show
that adjoining all split points leaves the displayed subgroup ranks at
`25,26,27,28,21`. These are not upper bounds for the full curves. See
[`elliptic-curves/notes/ELKIES_BISECTION_SPECIALIZATION_CONTROLS.md`](elliptic-curves/notes/ELKIES_BISECTION_SPECIALIZATION_CONTROLS.md).

Build the quotient-first rank-jump fingerprints for the four R17 controls:

```bash
PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_elkies_2026_rank_jump_fingerprints.py

PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_elkies_2026_rank_jump_fingerprints.py --check
```

Expected terminal status:

```text
ELKIES2026RANKJUMPFINGERPRINTS|quotient_ranks=8,9,10,11|degree2_visible=5,3,2,1|status=PASS_CERTIFIED_SUBGROUP_QUOTIENT_FINGERPRINTS|output=artifacts/generated-results/elliptic-curves/elkies_2026_rank_jump_fingerprints_v1.json
```

The Smith and cover-code fields are exact for the displayed certified
subgroups. The projected canonical-height Grams and successive minima use 80
decimal digits. Degree-three and degree-four cover visibility are recorded as
missing, not zero. See
[`elliptic-curves/notes/ELKIES_RANK_JUMP_FINGERPRINTS.md`](elliptic-curves/notes/ELKIES_RANK_JUMP_FINGERPRINTS.md).

Resolve the rank-28 visibility quotient, prove that translated trace shells
cannot enlarge it, and run exact `j`-recognition for the 2024 rank-29 curve
and ICARM 273, 302, and 398--400 against the published `R17` fibration:

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-VISIBILITY-RECORD-CURVES 1c39220ee5fedc77 -->

```bash
.venv/bin/python \
  elliptic-curves/scripts/analyze_elkies_bisection_visibility_and_record_curves.py \
  --check
```

The exact row reduction gives a ten-dimensional canonical complement to the
sole visible rank-28 class.  All six recognition equations are primitive
irreducible degree-24 polynomials.  The rank-28 control instead has the exact
factor `5471*t+9529`, recovering `t=-9529/5471`.  Thus none of the six target
curves is a rational fibre of the published fibration, including after
quadratic twisting.  This does not exclude another fibration, another family,
or an isogeny construction.  See
[`elliptic-curves/notes/ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md`](elliptic-curves/notes/ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md).

Exclude cyclic degree `3`, `5`, `7`, and `11` isogeny images of rational
published-R17 fibres for the same six target curves:

<!-- status-consumer: EC-K3-ELKIES-2026-R17-SMALL-ISOGENY-EXCLUSIONS fc2c4caaa79fb36c -->

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_elkies_2026_r17_isogeny_exclusions.sage --check
```

Each of the 24 tests has a clean prime at which the full projective
recognition polynomial has no finite linear factor and no root at infinity.
This excludes the displayed small prime degrees on the published fibration;
it does not test the equation-open alternate Q80-derived fibration, composite
isogeny degrees, or algebraic fibre parameters.  See
[`elkies-k3/PUBLISHED_R17_SMALL_ISOGENY_EXCLUSIONS_2026-09-01.md`](elkies-k3/PUBLISHED_R17_SMALL_ISOGENY_EXCLUSIONS_2026-09-01.md).

Canonicalize the known alternate rootless rank-17 frame in pinned H3/R17 NS
coordinates:

<!-- status-consumer: EC-K3-H3-OTHER-R17-J2-CANDIDATE f1884d1f6168a934 -->

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python -c \
  'import sys; from sage.all import *; sys.argv=["verify_q80_alternate_fifth_q6_rootless.sage","--write-artifact"]; globals()["__file__"]="/home/royvanrijn/src/jacobian-research/elkies-k3/scripts/verify_q80_alternate_fifth_q6_rootless.sage"; load(__file__)'

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python -c \
  'from sage.all import *; globals()["__file__"]="/home/royvanrijn/src/jacobian-research/elkies-k3/scripts/canonicalize_other_rank17_candidate.sage"; load(__file__)'
```

This proves a second rootless `J2` frame class of rank 17 and determinant 948,
not a complete fibration classification or a rank-29 fibre identification.
See
[`elkies-k3/OTHER_RANK17_FIBRATION_RECOVERY_2026-08-31.md`](elkies-k3/OTHER_RANK17_FIBRATION_RECOVERY_2026-08-31.md).

Replay the exact equation-compilation handoff for the final alternate q6
divisor, including its physical zero, complete NS transport, rootless
MW17 frame, and the explicit equation-open boundary:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_q80_alternate_final_divisor_handoff.sage --check
```

Audit the denominator-integral normalization of the exact quadratic descent
field at the current third-q12 equation frontier:

```bash
python3 \
  elkies-k3/scripts/audit_q80_third_q12_descent_field_normalization.py --check
```

The latter proves `QQ(sqrt(q1*q2))=QQ(sqrt(D))` for the displayed reduced
denominator `D`, without factoring it.  It does not claim an exact third-q12
Jacobian or a coefficient-height improvement.  See
[`elkies-k3/Q80_ALTERNATE_ROOTLESS_EQUATION_HANDOFF_2026-09-03.md`](elkies-k3/Q80_ALTERNATE_ROOTLESS_EQUATION_HANDOFF_2026-09-03.md).

Run the corresponding 63-term raw coefficient-height gate:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_q80_third_q12_pencil_basis_heights.py --check
```

The denominator-integral `delta` basis is worse than the existing `omega`
basis on 54 of 63 terms and raises the raw maximum coordinate height by 5,929
bits.  Rational projective normalization reverses this to a small 10,888-bit
(about 0.7 percent) improvement, while exposing only 7--12 bits of integer
content.  This is not the large compression required for reconstruction and
does not test quadratic-field, base-`PGL2`, or integral-ideal normalization.

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-LINEAR-CONDUCTOR 957479f39bedd57b -->

Certify the first exact generic conductor factor of the third-q12
discriminant and replay the retained quartic-denominator candidate:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --base-value 0 --certify-generic-linear --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_q80_third_q12_quartic_denominator_candidate.sage \
  --check
```

The first command proves over the full characteristic-zero base that the
reconstructed linear factor occurs with exact multiplicity three.  The
second verifies directly from the exact pencil at untouched inert primes
163, 191, and 199 that all four nonleading coefficients of the exponent-two
quartic have the predicted common linear denominator.  That denominator is
still candidate data: exact characteristic-zero recovery and division of the
quartic square, the Jacobian, and the maps remain open.  See
[`elkies-k3/Q80_ALTERNATE_ROOTLESS_EQUATION_HANDOFF_2026-09-03.md`](elkies-k3/Q80_ALTERNATE_ROOTLESS_EQUATION_HANDOFF_2026-09-03.md).

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-SPECIALIZED-QUARTICS 725664f9e36ae8a7 -->

Continue the exact quartic-square recovery with a resumable Brown
subresultant sequence:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --attempt-subresultant-prs \
  --output artifacts/generated-results/elkies-k3-q80-third-q12-exact-discriminant-specialization-v1.json
```

This long command writes its latest exact state atomically under
`artifacts/local/elkies-k3/` and resumes it by default.  Use
`--restart-subresultant-prs` only when intentionally discarding compatibility
with the retained state.  The completed command proves at `V=0` that the
linear-stripped monic discriminant is `Q^2*D` for exact monic quartics `Q,D`.
Replay the completed checkpoint, exact division, and 18 MB artifact by adding
`--check`.  Generic recovery over the full `V`-line, the Jacobian, and its maps
remain open.

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-GENERIC-QUARTICS aa704dc4685e4c9b -->

Recover and certify the full generic exponent-two quartic by exact first-order
lifting from `V=0` followed by fraction-free two-variable division:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --attempt-generic-quartic-division \
  --output artifacts/generated-results/elkies-k3-q80-third-q12-exact-generic-quartic-factorization-v1.json
```

This reconstructs all sixteen cleared numerator coordinates, replays the
`19^12288` lift, and proves exactly that the full characteristic-zero
linear-stripped discriminant is divisible by `Q(V,W)^2`.  The complementary
factor is polynomial and quartic in `W`.  The Jacobian, minimization, and maps
remain subsequent gates.

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-COMPLETE c6f054948b04b507 -->
<!-- status-consumer: EC-K3-H3-ROOTLESS-J1-UNIFORM-BOUND b71330a75ad2c9ad -->

Audit the two-control rootless `J2` corpus, its exact local genus, the mass
obstruction to an unfiltered full-genus traversal, and the complete
Niemeier-first classification:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_rootless_j2_completeness_track.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rootless_j2_niemeier_controls.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rooted_niemeier_catalog.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_niemeier_d5_anchor_orbits.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_niemeier_auxiliary_sixth_dominant.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_rootless_j2_niemeier_first.sage --check

python3 elkies-k3/scripts/certify_rootless_j1_uniform_bound.py --check
```

The target genus mass is
`77731517730627488307787/925557271717478400`, which forces at least 167,967
full-genus isometry classes because every automorphism group contains
`+/-identity`.  The audit retains the two mandatory non-isometric rootless
controls and rejects all 65 old `target-genus-rootless-pneighbor` files as one
different 2-adic and 79-adic genus (nineteen isometry classes there).  This is
an exact corpus correction and route-selection certificate, not a complete
Kneser--Nishiyama classification.  The second replay pins the unique required
rank-seven genus, certifies all cyclic discriminant gluings for both positive
controls, and places every control embedding in `N(2A7+2D5)`.  It finds two
glue double-coset orbits for the published frame and one for the alternate,
while explicitly retaining primitive embeddings and saturated complements.
Its exact 329,206,692-vector-pair cost probe shows why the complete
enumeration had to be orbit based.  The final four replays certify all 23 rooted
Niemeier Gram models, reduce the thirteen `D5`-admissible classes to sixteen
anchor orbits, enumerate 3,220 primitive residual-Weyl sixth vectors, and
reduce 167 positive-label seventh-vector cases to twelve primitive rootless
embeddings.  Exact integral-isometry deduplication gives exactly two frame
classes, the published R17 control and the alternate Q80 control, both with
`N(2A7+2D5)` provenance.  This is a complete rootless `J2` frame-isometry
classification; the retained embedding cover is not deduplicated to full
automorphism embedding-orbit or `J1` counts.  See
[`elkies-k3/ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md`](elkies-k3/ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md).
The final dependency-free replay uses rank-three Hodge rigidity and the exact
eight-element discriminant-form isometry group to give a uniform `J1`
multiplicity bound of four per frame.  Thus the complete rootless `J1` count
is rigorously in `[2,8]`; exact surface-automorphism representatives remain
open.

Build or byte-check the surface-first rank-seven auxiliary catalogue:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_prefix_orbits.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_prefix_orbits.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/run_24a1_octad_rank7_completion_frontier.py --jobs 3

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/run_24a1_octad_rank7_completion_frontier.py --jobs 3 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 0 --prefix-stop 250

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 0 --prefix-stop 250 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 250 --prefix-stop 500

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 250 --prefix-stop 500 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 500 --prefix-stop 750

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 500 --prefix-stop 750 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 750 --prefix-stop 1000

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 750 --prefix-stop 1000 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1000 --prefix-stop 1250

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1000 --prefix-stop 1250 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1250 --prefix-stop 1500

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1250 --prefix-stop 1500 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1500 --prefix-stop 1750

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1500 --prefix-stop 1750 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1750 --prefix-stop 2000

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1750 --prefix-stop 2000 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_24a1_octad_completion_manifest.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_24a1_octad_completion_manifest.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_24a1_weyl_m24_shard.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_24a1_weyl_m24_shard.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_4a_fixed_rank7.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_4a_fixed_rank7.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_4d6_swap_fixed_high_mw_seed.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_4d6_swap_fixed_high_mw_seed.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_6a4_double_swap_fixed_high_mw_seed.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_6a4_double_swap_fixed_high_mw_seed.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_4a5_d4_order4_fixed_high_mw_seed.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_4a5_d4_order4_fixed_high_mw_seed.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_4a6_4e6_residual_sections.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_4a6_4e6_residual_sections.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_4a6_4e6_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_4a6_4e6_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_4a6_4e6_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_4a6_4e6_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_8a3_glue_code_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_8a3_glue_code_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_8a3_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_8a3_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_8a3_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_8a3_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_6d4_hexacode_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_6d4_hexacode_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_6d4_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_6d4_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_6d4_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_6d4_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3e8_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3e8_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3e8_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3e8_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3d8_glue_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3d8_glue_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3d8_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3d8_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_3d8_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_3d8_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2d12_glue_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2d12_glue_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2d12_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2d12_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_d10_2e7_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_d10_2e7_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_d10_2e7_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_d10_2e7_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a12_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a12_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a12_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a12_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a9_d6_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a9_d6_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a9_d6_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_2a9_d6_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_2a9_d6_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_2a9_d6_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3a8_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3a8_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3a8_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_3a8_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_3a8_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_3a8_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_12a2_ternary_golay_residual_group.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_12a2_ternary_golay_residual_group.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_12a2_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_12a2_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_12a2_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_12a2_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_eta_only_niemeier_residual_groups.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_eta_only_niemeier_residual_groups.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_eta_only_niemeier_fixed_coordinate_shells.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_eta_only_niemeier_fixed_coordinate_shells.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_cross_niemeier_mod2_priority.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_cross_niemeier_mod2_priority.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_co0_backend.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_co0_backend.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_minimal_line_action.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_minimal_line_action.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_leech_minimal_basis_coordinate_shell.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_leech_minimal_basis_coordinate_shell.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_t_arithmetic.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_t_arithmetic.sage --check

python3 \
  elkies-k3/scripts/compare_leech_coordinate_targets_to_rooted_catalogue.py

python3 \
  elkies-k3/scripts/compare_leech_coordinate_targets_to_rooted_catalogue.py --check

python3 elkies-k3/scripts/build_rank7_surface_pareto.py

python3 elkies-k3/scripts/build_rank7_surface_pareto.py --check

# Ordered all-Niemeier factory; stops before equation solvers.
python3 elkies-k3/scripts/build_rank7_all_niemeier_factory.py

python3 elkies-k3/scripts/build_rank7_all_niemeier_factory.py --check
```

The first two commands enumerate exact residual-`M24` orbits of unordered
Golay-octad subsets through size five; the orbit counts are
`1,3,16,206,10547`, and every orbit-stabilizer mass equals
`binomial(759,k)`. The next two commands run or byte-check all 10,547
five-prefix orbits in 43 contiguous shards. The following sixteen commands
are independent first-eight-shard examples and may be omitted after the full
runner. After exact saturation and local residual-`M24` deduplication the full
frontier retains 3,051 records, with MW-rank distribution `387,2423,241` in
ranks 12, 13, and 14; 3,015 pass the ternary genus gate. The next two commands
build and byte-check the gap-free, hash-pinned shard manifest. The following
two commands apply the full `2^24 semidirect M24` quotient across all 43
shards: the 3,051 local records become 24 full embedding orbits, with MW-rank
distribution `5,13,6` in ranks 12, 13, and 14, and 18 pass the ternary-genus
gate. These yield two `(T,NS)` keys and add one global surface, two auxiliaries,
and eighteen frames after catalogue deduplication. The next two commands close the primitive corank-one
family inside the common rank-eight `4A` fixed lattice of `N(2A7+2D5)`.
They enumerate 336 embeddings, with MW-rank distribution `16,16,304` in
ranks 13, 15, and 17. Every embedding has nontrivial `2B,2C,4A` complement
action modulo two, but all five frame classes fail the ternary gate because
their discriminant groups have length seven. The following two commands test
all 11,440 coordinate direct summands of a pinned `2C` fixed-lattice basis.
The 97 high-MW survivors give 97 exact `Dih_4`-section orbits and 73 admissible
`(T,NS)` surfaces; post-deduplication they contribute 86 frames with MW-rank
distribution `16,66,1,3` in ranks 12 through 15. The next two commands test
the analogous `4D6` component-transposition coordinate shell. Of 11,440
summands, 183 fail discriminant length and all remaining 11,257 have MW rank
below 12, an exact negative result only for that shell. The following two
commands derive all 240 chamber-preserving residual automorphisms of `N(6A4)`
and canonicalize a double-swap fixed-lattice shell under the whole group. Its
161 MW12--13 orbits yield 42 local `(T,NS)` surfaces, 55 auxiliary classes,
and 128 frames. All selected involutions act nontrivially modulo two. The
next two commands derive the 48 chamber-preserving residual
automorphisms of `N(4A5+D4)` and close an order-four fixed-lattice coordinate
shell.  Its 39 high-MW seeds give 39 residual orbits and nine local
`(T,NS)` surfaces, auxiliary classes, and rootless frames, with distribution
`MW13:7, MW17:2`; all nine surfaces are globally new.  The order-four action
moves eleven complement dimensions modulo two.  The next six commands recover
the exact `N(4A6)` and `N(4E6)` residual groups, exhaust every fixed-coordinate
shell attached to a nonidentity matrix conjugacy class of fixed rank at least
seven, and apply the residual, ternary, and `(T,NS)` gates.  The groups have
orders 24 and 48 with component images `A4` and `S4`.  Of 26,064 coordinate
summands, only order-three classes give high-MW seeds: 86 for `4A6` and 48 for
`4E6`.  Residual canonicalization and the ternary gate yield nine `4A6`
surfaces with ten frames (`MW12:2, MW13:8`) and one `4E6` MW12 surface/frame;
all ten are globally new. The next six commands recover the order-256
`N(8A3)/A3^8` glue code, exhaust all `2^8 8!` signed component maps, and
certify the complete order-2,688 residual group with order-1,344 component
image. Seven eligible fixed-rank classes contribute 24,600 coordinate
summands and 1,166 high-MW/mod-two seeds. Full residual canonicalization gives
1,162 orbits, of which 1,160 pass the ternary gate; `(T,NS)`-first
deduplication gives 435 local surfaces, 523 auxiliary classes, and 574 frames.
Twenty-four surfaces overlap prior backends, so the global catalogue gains
411. The next six commands recover the order-64 `N(6D4)/D4^6` hexacode,
exhaust all `6^6 6!` triality-permutation maps, and certify the complete
order-2,160 residual group with kernel three and full `S6` component image.
Eleven eligible classes contribute 25,416 coordinate summands and 472
high-MW/mod-two seeds. Full quotienting gives 466 orbits, 456 K3-compatible
orbits, and 218 local surfaces with 255 auxiliaries and 289 frames. Fifty
surfaces overlap prior backends, so `6D4` adds 168 global surfaces, 251
auxiliaries, and 285 frames. The next four commands intrinsically recover the
three `E8` components of `N(3E8)`, certify the complete order-six residual
`S3`, and scan both nonidentity fixed-coordinate classes. Their fixed ranks
16 and 8 contribute 11,448 rank-seven coordinate summands; none survives the
determinant/length/MW12--17/nontrivial-mod-two gates, so canonicalization and
catalogue merging are vacuous for this declared shell. The next six commands
recover the order-eight `N(3D8)/D8^3` glue code, reject 42 of the 48 possible
component/diagram maps, and certify the natural order-six residual `S3`.
The two nonidentity classes again contribute 11,448 coordinate summands. The
transposition class gives 40 qualified seeds (`MW12:28, MW13:12`), while the
three-cycle class is negative. Full residual quotienting gives 40 orbits; 25
pass the ternary gate and `(T,NS)`-first deduplication gives seven local
surfaces, auxiliaries, and frames. One surface overlaps the existing
`2A7+2D5`/`6A4` class, but its auxiliary and frame are new; the other six
surfaces are new. The next four commands recover the order-four
`N(2D12)/D12^2` glue code, reject six of eight component/diagram maps, and
certify the natural order-two component-swap residual group. Its rank-twelve
fixed lattice contributes 792 coordinate summands, none of which survives
all determinant/length/MW12--17/nontrivial-mod-two gates. Canonicalization and
catalogue merging are therefore vacuous for this declared shell. The next
four commands recover the mixed `D10+2E7` root decomposition and test all four
diagram/permutation chamber maps. Glue rejects the `D10` flip and `E7` swap
separately and retains only their simultaneous product. Its rank-sixteen
fixed lattice contributes 11,440 coordinate summands, again with no seed
surviving all pre-quotient gates. The next four commands recover the
index-thirteen `N(2A12)/A12^2` quotient and certify the cyclic order-four
residual group. Its order-four classes have fixed rank six; the sole eligible
order-two class has a rank-twelve fixed lattice whose 792 coordinate summands
also yield no seed. The next six commands recover the index-twenty
`N(2A9+D6)/(A9^2+D6)` quotient, test all sixteen component/diagram maps, and
certify the complete cyclic order-four residual group. The three nonidentity
fixed lattices have ranks `16,10,10`; their 11,680 coordinate summands yield
no order-two seed and 32 qualified seeds for each inverse order-four class.
Exact residual quotienting gives 32 orbits, 13 K3-compatible orbits, and five
local surfaces, auxiliaries, and frames (`MW13:4, MW17:1`). One surface
overlaps `8A3`; four surfaces but all five auxiliaries and frames are new.
The next six commands recover the index-27 `N(3A8)/A8^3` quotient, test all
48 component/diagram maps, and certify the exact order-twelve residual group
`{+/-1} x S3`. Four nonidentity classes have fixed rank at least seven; their
13,032 coordinate summands yield 189 qualified transposition-fixed seeds
(`MW12:135, MW13:54`) and no seed in the central-involution,
signed-transposition, or three-cycle shells. Exact residual quotienting gives
189 size-six orbits, all K3-compatible, and 25 local surfaces with 30
auxiliaries and 64 frames. Twenty surfaces overlap earlier backends, so five
surfaces, twelve auxiliaries, and all 64 frames are new. The next six commands
recover the intrinsic ternary Golay `[12,6,6]` glue of `N(12A2)`, certify its
full residual group `2.M12` of order 190,080 and 26 conjugacy classes, and
scan all 13,968 eligible fixed-coordinate summands. The 237 qualified seeds
become 214 exact residual orbits, 210 K3-compatible orbits, and 99 local
surfaces with 108 auxiliaries and 151 frames. Forty-seven surfaces overlap
earlier backends, so 52 surfaces, 86 auxiliaries, and 143 frames are new.
The next four commands intrinsically recover all six remaining rooted
residual groups. Exhausting their component diagram maps proves trivial
groups for `D24` and `D16+E8` and order-two eta groups for `A24`, `A17+E7`,
`A15+D9`, and `A11+D7+E6`. The four nontrivial fixed-coordinate languages
contain 35,112 subsets but no MW12--17 seed: `A24` fails discriminant length,
and every length-admissible mixed-system frame has MW rank below 12. The
following two commands
build the cross-Niemeier mod-2 scheduler, prioritizing `2B`, `2C`, `4A` and
analogous component permutations but requiring the exact gate
`rank_GF2(g_M-I)>0`.  These fixed-lattice shells
remain bounded coordinate languages; they are not full
backend or Weyl-orbit censuses. Likewise, the `24A1` prefix frontier and
non-positive-octad generator languages remain open.

The following two commands pin the separate Leech ambient directly from the
AtlasRep `2.Co1=Co0` action: the invariant form is one-dimensional and its
primitive positive integral generator is the even unimodular rank-24 lattice
of minimum four with 196,560 minimal vectors. The next four commands certify
the exact `Co1` action on 98,280 antipodal minimal lines and exhaust all
346,104 rank-seven coordinate summands of one norm-four determinant-one
ambient basis. The finite language gives 221 signed-basis types, 194
ternary-compatible types, and 150 preliminary `(T,NS)` keys; it is still
pre-`Co1` and is not an all-primitive Leech census. The later crosswalk finds
43 exact rooted-catalogue matches, including 25 with a catalogued MW12 rooted
frame, and leaves 107 keys absent from the current catalogue. No rank-seven
Leech embedding orbit is claimed before the Conway quotient.

The catalogue imports exact primitive embeddings from the current
`N(2A7+2D5)` mutation and `2C` fixed-lattice seed shells, the `N(6A4)`
double-swap, `N(4A5+D4)` order-four, and
`N(2A9+D6)`/`N(3A8)`/`N(3D8)`/`N(4A6)`/`N(4E6)`/`N(6D4)`/`N(8A3)`/`N(12A2)` all-class
coordinate shells, the exact negative `N(2A12)`, `N(2D12)`, `N(D10+2E7)`,
`N(3E8)`, `N(A24)`, `N(A17+E7)`, `N(A15+D9)`, and `N(A11+D7+E6)` coordinate
scans, and the
determinant-720 `N(24A1)` Golay design,
deduplicates first by `(T,NS)` and only then by auxiliary/frame isometry, and
retains legacy `NS....` identifiers as aliases. Three of the 42 local `6A4`
surfaces overlap prior `2A7+2D5` classes, while all nine local `4A5+D4`
and all ten `4A6`/`4E6` surfaces are new. Of the 435 local `8A3` surfaces, 24
overlap earlier backends and 411 are new. Of 218 local `6D4` surfaces, 50
overlap and 168 are new. Of seven local `3D8` surfaces, one overlaps and six
are new, while all seven auxiliaries and frames are new. Of five local
`2A9+D6` surfaces, one overlaps `8A3` and four are new, while all five
auxiliaries and frames are new. Of 25 local `3A8` surfaces, 20 overlap earlier
backends and five are new; twelve auxiliaries and all 64 frames are new. The
99 local `12A2` surfaces include 47 overlaps, so 52 are new; 86 auxiliaries
and 143 frames are new. The completed positive seven-octad `24A1` subfamily
adds one new determinant-480 surface and thirteen frames, while its
determinant-500 surface overlaps `12A2` and gains five frames. The global
catalogue therefore has 827 exact surface classes, 1,074 auxiliaries, and
1,840 MW12--17
frames. The 23 rooted
backends and separate Leech
backend are split into 96 determinant-band shards; all remain open because the
imported searches have bounded, narrower completeness statements. See
[`elkies-k3/RANK7_AUXILIARY_CATALOGUE_2026-09-01.md`](elkies-k3/RANK7_AUXILIARY_CATALOGUE_2026-09-01.md).

The final two commands build and byte-check the typed surface-wide discovery
ordering. Its exact four-metric core frontier has fourteen of 827 surfaces.
Coverage-restricted enriched ledgers contain 39 surfaces with exact
minimum-pole evidence and 787 with nontrivial stabilizer evidence; no certified
physical neighbour route is currently available. Missing equation, field,
route, conductor, and moduli data are retained as typed unknowns rather than
imputed.

Build or byte-check the first determinant-varying Picard-19 lattice-foundry
shell:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_lattice_foundry.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_lattice_foundry.sage --check
```

The JSON specification declares exactly 768 one-root mutations of the two
stored `N(2A7+2D5)` control embeddings.  Saturation, complements, root data,
integral-isometry deduplication, and ternary discriminant-form realization are
exact inside that shell.  The generator is not a complete determinant-5,000
classification.  See
[`elkies-k3/LATTICE_FOUNDRY_REPORT_2026-09-01.md`](elkies-k3/LATTICE_FOUNDRY_REPORT_2026-09-01.md).

The default low-MW source-discovery workflow is the direct prescribed-root
Niemeier enumerator
`elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage`.
It fixes each selected NS class's rank-seven auxiliary and searches its
primitive Niemeier embeddings with the desired rank-15--17 complement roots
built into the enumeration.  Every retained source is attached to all
catalogued MW15--17 target frames in the same NS class.  Unlike the older
cyclic-gluing Kneser scout, this workflow also covers the foundry classes with
noncyclic discriminant group.

The calculation is exact inside the root types, support counts, ambient
classes, and auxiliary-embedding shell declared by the command.  It certifies
the primitive auxiliary embedding, saturated complement, complete complement
root system, and resulting geometric MW rank.  It does not by itself prove a
rational marking, an equation over a number field, or a physical
elliptic-neighbour route, and a miss is not a classification beyond that
declared slice.  Exact repetitions with the same deterministic reduced Gram
are merged; distinct reduced Grams are not claimed to be distinct integral-
isometry or `J2` classes.

Run or byte-check the narrow determinant-948 `NS0001` positive control over
all rooted Niemeier ambients with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ambients \
  --source-root-rank-min 15 --source-root-rank-max 15 \
  --source-support-min 2 --source-support-max 2 \
  --require-hit --require-h3-control \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-ns0001-all-ambients-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ambients \
  --source-root-rank-min 15 --source-root-rank-max 15 \
  --source-support-min 2 --source-support-max 2 \
  --require-hit --require-h3-control \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-ns0001-all-ambients-v1.json \
  --check
```

Here `--require-h3-control` requires an exact `E7+E8/MW2` root profile and an
integrally isometric binary Mordell--Weil height form.  It is deliberately not
a full rank-17 frame-isometry assertion.

Run or byte-check the declared production slice over every foundry NS class,
restricted to two- or three-support all-`A` sources in `N(3E8)`, with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ns --ambient-label 3E8 \
  --source-support-min 2 --source-support-max 3 --all-a-only \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ns --ambient-label 3E8 \
  --source-support-min 2 --source-support-max 3 --all-a-only \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json \
  --check
```

This exact declared slice finds 97 reduced-Gram MW2 source representatives in
23 of the 48 NS classes: 64 with root type `A2+A6+A7` and 33 with root type
`A1+2A7`.  It finds no MW0 or MW1 source in this slice.  Distinct reduced
Grams may still duplicate an integral-isometry or `J2` class, and the 25
slice misses are not non-existence theorems.

The full-support rank-16/17 (MW0--1 at Picard rank 19) census is split into
four disjoint ambient shards so it can run in parallel:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ns --ambient-label 2A7_2D5 --ambient-label 2A9_D6 \
  --ambient-label 2D12 --ambient-label 3D8 \
  --source-root-rank-min 16 --source-root-rank-max 17 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-a-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ns --ambient-label 3E8 --ambient-label 4D6 \
  --ambient-label 4E6 --ambient-label A11_D7_E6 \
  --source-root-rank-min 16 --source-root-rank-max 17 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-b-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ns --ambient-label A15_D9 --ambient-label A17_E7 \
  --ambient-label D10_2E7 \
  --source-root-rank-min 16 --source-root-rank-max 17 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-c-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ns --ambient-label D16_E8 --ambient-label D24 \
  --source-root-rank-min 16 --source-root-rank-max 17 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-d-v1.json
```

Append `--check` to any shard command for an exact byte replay.  Audit the
four shards and regenerate the compact aggregate accounting with

```bash
python3 elkies-k3/scripts/summarize_lattice_foundry_prescribed_root_shards.py \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-a-v1.json \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-b-v1.json \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-c-v1.json \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-d-v1.json \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-all-ambients-summary-v1.json
```

The summary distinguishes shard-local occurrences from exact repeated
`(NS, reduced-Gram digest)` identities across ambients.  It does not identify
unequal reduced Grams up to integral isometry or `J2`.

The audited census has 2,134 reduced-Gram source representatives, all MW1,
covering all 48 NS classes.  It finds no MW0 representative.  There are no
repeated exact `(NS, reduced-Gram digest)` identities between shards.  Of the
MW1 inventory, 245 representatives in 33 NS classes have at most two supports
and only `A`-type components.  These are exact results inside the declared
thirteen-ambient, sixteen-D5-anchor, sixth-norm-at-most-24 embedding cover;
they are not pairwise integral-isometry counts or a global MW0 non-existence
theorem.  The aggregate summary has SHA-256
`c5e610ac5baf12e01f86d506a6b42b6593a48f8949311eee095dfc27b55f9ad6`.

Audit the minimum section pole on every MW1 row with primitive root lattice:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_lattice_foundry_rank1_section_poles.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_lattice_foundry_rank1_section_poles.sage --check
```

The audit certifies exact affine-CVP norms, section poles, and independent
256-bit MPFR agreement for 1,342 primitive-root rows.  It leaves 792
nonprimitive-root rows open rather than guessing their torsion/glue section
lattice.  The result is
`artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json`.

Audit the cheapest complete MW basis on all 97 primitive rank-two rows:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_lattice_foundry_rank2_section_basis_poles.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_lattice_foundry_rank2_section_basis_poles.sage --check
```

The exact height bound makes the tail enumeration finite.  Each affine CVP is
repeated with double-double and MPFR-256 arithmetic, every returned frame norm
is recomputed integrally, and all determinant-one tail pairs are compared.
The result is
`artifacts/generated-results/elkies-k3-lattice-foundry-rank2-section-basis-poles-v1.json`
with SHA-256
`387a95156a06acb342fa4233f0aa69fb09e4c53fd69020a19e686da7cc4bcf38`.

Replay the NS0011 `A2+A6+A8/MW1` equation gate over `GF(5)` with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0011_source_ansatz_modp.sage --examples 7

python3 - <<'PY'
import json, subprocess
from pathlib import Path
source = Path('artifacts/generated-results/elkies-k3-lattice-foundry-ns0011-source-ansatz-mod5.json')
data = json.loads(source.read_text())
exe = '/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
script = 'elkies-k3/scripts/build_lattice_foundry_ns0011_pole2_section_modp.sage'
for example, model in enumerate(data['examples']):
    a8 = model['A_coefficients_low_to_high'][8]
    b12 = model['B_coefficients_low_to_high'][12]
    nonzero = sum((c**3 + a8*c + b12) % 5 in (1, 4) for c in range(5))
    for branch in range(nonzero):
        subprocess.run([exe, script, '--example', str(example), '--branch',
                        str(branch), '--eliminate-r-count', '9'], check=True)
    subprocess.run([exe, script, '--example', str(example),
                    '--zero-y-branch', '0'], check=True)
PY

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0011_pole2_sections_modp.sage
```

The first command exhausts all `5^8` normalized fibre coefficients and stores
all seven exact squarefree `I9+I7+I3+5I1` models.  The builder encodes the
exact I9/I7 component depths as a sparse pole-two square-root system.  The
final tensor-product scan exhausts 19 infinity charts and 1,484,375 affine
tuples.  Its empty result is a finite-field obstruction for this displayed
split normalized chart only; it is not a characteristic-zero nonexistence
claim.  The bounded `GF(7)` pilot and its exact section-chart census are stored
as the corresponding `mod7-pilot-v1.json` artifacts; its fibre scan covers
only 500,000 of `7^8` coefficient tuples.

Replay the still cheaper NS0007 `A1+A3+2A6/MW1` fibre gate with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage \
  --examples 30

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage \
  --examples 30 --check
```

This source has exact height `11/4`, minimum pole zero, component corrections
`1/2,3/4,0,0`, and rootless same-NS MW17 endpoints `NS0007-F003` and
`NS0007-F014`.  The complete split `GF(5)` scan exhausts all three
cross-ratios and `3*5^8` normalized `A` polynomials.  It finds 966
Hermite-compatible branches but no squarefree `I2+I4+2I7+4I1` model.  This is
an exact obstruction for the displayed characteristic-five fibre chart, not
a characteristic-zero nonexistence theorem.  The corresponding
`mod7-pilot-v1.json` artifact covers only 100,000 coefficient/cross-ratio
cases and is a bounded negative result.

Replay the bounded arithmetic-support samples, in which the two `I7` fibres
form one irreducible quadratic Frobenius orbit, with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage \
  --prime 7 --support-configuration conjugate-i7 \
  --max-a-samples-per-lambda 100000 --sample-stride 104729 \
  --sample-offset 1 --examples 20 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0007-source-ansatz-conjugate-i7-mod7-sample-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage \
  --prime 11 --support-configuration conjugate-i7 \
  --max-a-samples-per-lambda 100000 --sample-stride 104729 \
  --sample-offset 17 --examples 20 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0007-source-ansatz-conjugate-i7-mod11-sample-v1.json
```

The prime-7 run checks four support orbits and 400,000 normalized `A`
polynomials; 133 signed branches pass the descended Hermite equations and none
has the exact prescribed fibre orders.  The prime-11 run checks six support
orbits and 600,000 polynomials; its corresponding counts are 24 and zero.
Neither sample produces a squarefree residual quartic.  These are bounded
negative arithmetic-chart experiments, not exhaustive prime obstructions.
The `A(0)=-3` normalization fixes only one local `I2` twist at each prime,
while the `I4` and conjugate `I7` tangent characters are unrestricted.

Build and audit the globally reduced fixed-`lambda=2` `GF(7)` pole-zero
system with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_lattice_foundry_ns0007_pole0_reduced_modp.sage \
  --prime 7 --lambda-value 2

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_lattice_foundry_ns0007_pole0_reduced_modp.sage \
  --prime 7 --lambda-value 2 --compact-factored-msolve

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_lattice_foundry_ns0007_compact_msolve_encoding.sage \
  --check
```

The audit proves equality of all 19 equations after Sage parses and expands
them.  It does not license the compact syntax for msolve: that syntax produced
12 false nonunit cases in the first 10,000 assignments, while the expanded
input proves that all 10,000 are unit ideals.  The census runner refuses
factored input by default.  Replay the pinned expanded prefix with:

```bash
python3 \
  elkies-k3/scripts/run_lattice_foundry_ns0007_p7_fixed_case_census.py \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0007-pole0-fixed-case-leading-prefix-10000-mod7-v1.json \
  --workers 16 --solver-mode leading-ideal --timeout-seconds 10 \
  --start-index 0 --limit 10000 --progress-every 5000
```

On the pinned run, seven cases timed out under transient host pressure.  Each
was rerun as a singleton with a 60-second bound and resolved to the unit
ideal.  Replay the checked merge with:

```bash
python3 \
  elkies-k3/scripts/repair_lattice_foundry_ns0007_p7_fixed_case_census.py \
  --base artifacts/generated-results/elkies-k3-lattice-foundry-ns0007-pole0-fixed-case-leading-prefix-10000-mod7-v1.json \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0007-pole0-fixed-case-leading-prefix-10000-resolved-mod7-v1.json \
  artifacts/local/elkies-k3/ns0007-pole0-leading-replays/p7-lambda2-009065.json \
  artifacts/local/elkies-k3/ns0007-pole0-leading-replays/p7-lambda2-009066.json \
  artifacts/local/elkies-k3/ns0007-pole0-leading-replays/p7-lambda2-009067.json \
  artifacts/local/elkies-k3/ns0007-pole0-leading-replays/p7-lambda2-009068.json \
  artifacts/local/elkies-k3/ns0007-pole0-leading-replays/p7-lambda2-009069.json \
  artifacts/local/elkies-k3/ns0007-pole0-leading-replays/p7-lambda2-009071.json \
  artifacts/local/elkies-k3/ns0007-pole0-leading-replays/p7-lambda2-009072.json \
  --check
```

The resolved artifact has SHA-256
`7d30c0302547240427f8eecf8e5d38ffe36f0f11546f1104e5b2c2a44b18331d`.
It is an exact bounded prefix of the `7^6=117649` base-field assignments in
the displayed fixed-lambda chart, not an exhaustive chart obstruction and not
a characteristic-zero nonexistence result.  Long continuations should use
nonoverlapping `--start-index/--limit` shards and then
`combine_lattice_foundry_ns0007_p7_fixed_case_census.py`; the combiner accepts
only contiguous, fully expanded, terminal shards.

Build the reduced NS0034 nodal-Hermite system, replay the complete fixed
`GF(7)` fibre slice, and exhaust its pole-zero section chart with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_lattice_foundry_ns0034_fibre_hermite_reduced_modp.sage \
  --prime 7 --lambda-value 2 --A8 1 --hi0 3 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage \
  --candidate ns0034 --prime 7 --fixed-lambda-value 2 --fixed-A8 1 \
  --examples 100 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0034-source-ansatz-mod7-lambda2-A8-1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0034_pole0_sections_modp.sage \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-ns0034-source-ansatz-mod7-lambda2-A8-1.json \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0034-pole0-sections-mod7-lambda2-A8-1.json \
  --check
```

The fixed fibre slice exhausts `7^7=823543` normalized `A` polynomials.  Of
1,210,104 locally eligible signed branches, 505 satisfy the Hermite equations
and exactly one has the prescribed `I4+I8+I3+I5+4I1` orders with squarefree
residual quartic.  Its complete `7^8=5764801` section chart contains eight
polynomial sections but no `NS0034-S008` marked section.  This is an exact
obstruction only for `lambda=2,A8=1` in the displayed characteristic-seven
normalization.  The reduced 20-equation Hermite ideal is positive-dimensional
before residual-order saturation and is not itself a source-family
certificate.

Replay both local twist classes of the equation-friendlier NS0043
`A2+A6+A8/MW1` pole-zero marking on the complete `GF(5)` fibre census with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0043_pole0_sections_modp.sage \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0043_pole0_sections_modp.sage \
  --quadratic-twist 2 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0043-pole0-sections-mod5-nonsquare-twist.json \
  --check
```

The source has height four, zero component corrections, determinant 756, and
four same-NS MW15 targets.  The square twist contains no polynomial section;
the nonsquare twist contains 54, but every one meets a nonidentity component
at one or more of the `I9,I7,I3` supports.  Together these are an exact
two-twist obstruction for the displayed split `GF(5)` normalization, not a
characteristic-zero nonexistence result.  The corresponding two `mod7` pilot
artifacts exhaust the section charts only on the two stored bounded-sample
fibre models.

Replay the completed NS0030 pole-zero source census attached to an MW16 target:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0030_source_ansatz_modp.sage \
  --max-a-samples-per-support-pair 100000 --sample-stride 17 \
  --examples 100 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0030-source-ansatz-mod5-pilot100k-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0030_source_ansatz_modp.sage \
  --max-a-samples-per-support-pair 290625 --sample-stride 17 \
  --sample-offset 137500 --examples 100 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0030-source-ansatz-mod5-suffix290625-v1.json \
  --check

python3 elkies-k3/scripts/combine_lattice_foundry_ns0030_source_ansatz_modp.py \
  --check
```

The two adjacent segments of the same coprime-stride permutation cover all
`5^8` normalized `A` polynomials in each of the six ordered `GF(5)` support
pairs for the profile `2I2+I3+2I7+3I1`.  Among 2,343,750 rows, 1,536 signed
branches satisfy Hermite compatibility and none has the exact prescribed
orders.  This is a complete obstruction for the displayed normalized
characteristic-five chart, not a characteristic-zero nonexistence theorem.

Replay the current marked MW1--MW16 equation lead `NS0048-S030` with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0048_source_ansatz_modp.sage \
  --prime 5 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod5.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0048_source_ansatz_modp.sage \
  --prime 7 --examples 100 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod7.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0048_pole0_sections_xonly_modp.sage \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0048_pole0_sections_xonly_modp.sage \
  --quadratic-twist 3 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-pole0-sections-xonly-mod7-nonsquare-twist.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_lattice_foundry_ns0048_marked_family_modp.sage \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_lattice_foundry_ns0048_marked_family_modp.sage \
  --fibres artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod11-suffix600k-v1.json \
  --sections artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-pole0-sections-xonly-mod11-suffix600k-v1.json \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-marked-family-jacobian-mod11-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/lift_lattice_foundry_ns0048_marked_family_padic.sage \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_lattice_foundry_ns0048_tate_family_modp.sage \
  --prime 11 --fixed-lambda 10 --check
```

The source has root type `A1+A4+A6+D5`, MW height `37/14`, pole zero,
determinant 740, and same-NS MW16 endpoint `NS0048-F001`.  The complete
`GF(7)` fibre census has six models; the X-only section scan finds one marked
sign pair in the square twist and none in the nonsquare twist.  The marked
family has Jacobian rank 18 in 19 variables at both the characteristic-seven
point and a point in the bounded characteristic-eleven suffix.  The latter
lifts through `11^80` with one parameter fixed but does not rationally
reconstruct.  These are modular and p-adic source-locus certificates, not a
rational K3 equation or a neighbour corridor.

Replay the equation-level audit of the multisection-leading MW2--MW17
candidate NS0028 with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --prime 5 --examples 0 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0028_pole0_section_pairs_modp.sage \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0028_pole0_section_pairs_modp.sage \
  --quadratic-twist 2 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod5-nonsquare-twist.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --prime 7 --examples 0 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-source-ansatz-mod7.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0028_pole0_section_pairs_modp.sage \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-source-ansatz-mod7.json \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod7.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0028_pole0_section_pairs_modp.sage \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-source-ansatz-mod7.json \
  --quadratic-twist 3 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod7-nonsquare-twist.json \
  --check
```

The exact lattice source has profile `A2+A6+A7/MW2`; both generators have
pole zero and require only depth-one node contacts.  The complete fibre scans
find 25 models at 5 and 112 at 7.  Neither twist at 5 contains either marked
section.  At 7 the square twist has ten copies of one generator and the
nonsquare twist has ten copies of the other plus two of the first, but never
both on one model.  This is an exact two-prime obstruction for the displayed
normalized charts only.

## Section-first MW1/MW2 normal-form compiler

Generate and byte-check the generalized marked-section frontend with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_section_first_normal_forms.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_section_first_normal_forms.sage --check
```

The MW1 chart is
`y^2+a1*x*y+a3*y=x^3+a2*x^2` with the marked point `(0,0)`.
The MW2 chart additionally compiles
`Q=(h*r,h^2*s)` by a Bezout relation, so both point equations and their
affine intersection divisor `h` are identities before fibre tuning.  The
control artifact has SHA-256
`6edfa5f3487f020a05772bd2b1a6b5d74586c126d852b87418285ef46f843c34`.
It replays the rational Golay `3I6+6I1` model and the marked NS0031
`GF(7)` model 157, including their degree-two smooth intersections and exact
semistable fibre orders.  It does not change the existing saturation and
characteristic-zero proof boundaries.  See
[`elkies-k3/SECTION_FIRST_NORMAL_FORM_COMPILER_2026-09-02.md`](elkies-k3/SECTION_FIRST_NORMAL_FORM_COMPILER_2026-09-02.md).

For bounded-search provenance, replay or byte-check the exact `NS0024`
rootful source found by the older Kneser scout and its certified degree-two
route to a new rootless MW17 frame with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/hunt_lattice_foundry_rootful_source.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/hunt_lattice_foundry_rootful_source.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_lattice_foundry_route.sage \
  --manifest elkies-k3/data/lattice-foundry/ns0024-nef-route-v1.json
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_lattice_foundry_route.sage \
  --manifest elkies-k3/data/lattice-foundry/ns0024-nef-route-v1.json --check
```

The source is `5A1+A2+A5/MW5` in `N(A11+D7+E6)`. The route has eleven
primitive-nef old-degree-two edges, uses only `q=4,6`, requires no physical
Weyl repair, and carries full integral NS markings to a rootless MW17 frame.
Its equation-side resolved-RR dimensions remain planning estimates; no
Weierstrass equation is asserted.

Replay the exact `A3+A4+A6/MW4` source route and its modular fibre-ansatz
gate with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_lattice_foundry_route.sage \
  --manifest elkies-k3/data/lattice-foundry/ns0024-r13-nef-route-v1.json --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0024_source_ansatz_modp.sage --check
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0024_source_ansatz_modp.sage \
  --prime 13 --max-samples 200000 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-ansatz-mod13.json \
  --check
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0024_source_ansatz_modp.sage \
  --prime 17 --max-samples 300000 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-ansatz-mod17.json \
  --check
```

The MW4 route has thirteen primitive-nef degree-two edges, uses only
`q=4,6`, has zero physical Weyl repairs, and lands on rootless catalogue
frame `NS0024-F005`. The modular checks prove only the exact
`I7+I5+I4+8I1` fibre profile; four MW-section conditions, the `NS0024`
marking, and characteristic-zero lifting remain open.

Rank all stored same-NS source fibrations against the MW0--2-first equation
objective, and replay the low-degree multisection spectrum on the selected
rootless foundry targets, with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/score_lattice_foundry_sources.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/sample_lattice_foundry_multisection_spectrum.sage \
  --sample-count 256 --height-slack 4 \
  --frame-id NS0001-F001 --frame-id NS0002-F007 \
  --frame-id NS0005-F008 --frame-id NS0011-F002 \
  --frame-id NS0022-F011 --frame-id NS0024-F005 \
  --frame-id NS0028-F005 --frame-id NS0032-F011 \
  --frame-id NS0033-F026 --check
```

The source ranking attaches every catalogued MW15--17 frame in the same
Neron--Severi class and never infers a rational marking or route from the
lattice. It combines 75 individual source certificates, the 97 exact MW2 rows
of the declared `3E8`/all-`A` slice, and all 2,134 MW1 rows of the full-support
census. Primitive-root MW1 rows consume the exact pole audit; nonprimitive
rows retain an open pole status.  The 97 primitive MW2 rows consume the exact
complete-basis audit, which ranks the minimum possible maximum basis pole
before the cheapest single-section pole.  Thus NS0028 has profile `[0,0]`,
NS0005 `[0,1]`, and the equation-first NS0011 source `[1,1]`; a different
NS0011 row has a cheap individual section but only profile `[0,2]`. No row
receives inferred rational-marking or route data. The
multisection replay is complete for degree-two low-height translation orbits;
degree-three/four results in that pilot artifact are exact only for the
declared 256 sampled cosets per frame.

Recover the umbral stabilizer images and orbit-resolved short-vector and
multisection-coset data for R17, Q80, and four selected foundry frames with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_lattice_foundry_umbral_orbits.sage \
  --d3-orbit-seeds 512 --pari-stack-gb 4
```

The stabilizer, norm-four, and degree-two results are exact; degree three is
exact only inside the deterministic invariant sample. The interpretation and
comparison with lambency-eight umbral trace data are recorded in
[`elkies-k3/UMBRAL_COMPLEMENT_ORBIT_PILOT_2026-09-02.md`](elkies-k3/UMBRAL_COMPLEMENT_ORBIT_PILOT_2026-09-02.md).

Exhaust all `3^17` degree-three translation cosets on the pinned five-surface
batch selected by the earlier cheapest-single-section MW2 ranking with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --frame-id NS0028-F005 --frame-id NS0011-F002 \
  --frame-id NS0022-F011 --frame-id NS0005-F008 \
  --frame-id NS0001-F001 \
  --workers 8 --chunk-size 1000000 --float-type dd \
  --audit-precision 256 --audit-stride 4096 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --frame-id NS0028-F005 --frame-id NS0011-F002 \
  --frame-id NS0022-F011 --frame-id NS0005-F008 \
  --frame-id NS0001-F001 \
  --workers 8 --chunk-size 1000000 --float-type dd \
  --audit-precision 256 --audit-stride 4096 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json \
  --check
```

An initial run may add
`--reuse-checkpoint artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-top5-v1.json.partial`
to reuse the three already-complete frame blocks. The final artifact and its
adjacent `.partial` checkpoint are independent of that acceleration: the
check requires complete task coverage, exact histogram totals, pinned Gram
hashes, and the deterministic cross-precision audit accounting. The stable
whole-file SHA-256 of the current-top-five artifact is
`8be0e881f5c170366dada6319aed9a09fed689eacc032fcaf5ee70878d735fd0`.
See
[`elkies-k3/LATTICE_FOUNDRY_SOURCE_FIRST_OBJECTIVE_2026-09-01.md`](elkies-k3/LATTICE_FOUNDRY_SOURCE_FIRST_OBJECTIVE_2026-09-01.md).

Compute the first structure-sensitive R17 multisection profile, including the
complete degree-two quotient metric and graph, the complete degree-three
one-vertex spectrum, deterministic degree-three/four graph samples, exact
degree-two-to-four overlap, equation-complexity weights, and Gauss-local
squareclass signatures, with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_r17_multisection_diversity.py
```

The output is
[`artifacts/generated-results/elkies-k3-r17-multisection-diversity-v1.json`](artifacts/generated-results/elkies-k3-r17-multisection-diversity-v1.json).
Degree two is complete, including all `2^17` cosets, 8,895,801
zero-intersection edges, and 157,553,175 graph triangles.  The degree-three
graph and degree-four data outside the embedded two-torsion subset are sampled
and are labelled accordingly.  See
[`elkies-k3/R17_MULTISECTION_DIVERSITY_2026-09-02.md`](elkies-k3/R17_MULTISECTION_DIVERSITY_2026-09-02.md).

Run the exact degree-two diversity comparison on R17, `NS0032-F011`, and
`NS0028-F005` with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_r17_multisection_diversity.py \
  --comparison-only
```

The output is
[`artifacts/generated-results/elkies-k3-r17-foundry-degree2-diversity-comparison-v1.json`](artifacts/generated-results/elkies-k3-r17-foundry-degree2-diversity-comparison-v1.json).
All `2^17` cosets per lattice, automorphism orbits, quotient-pair
distributions, rational zero-intersection graph components and triangles, and
the inherited degree-two-to-four mass are complete.  The stable artifact
SHA-256 is
`408947ebb3e67048767005a005c3f283cc2bd2b4971e12006716809134a7146c`.

Calibrate the R17 graph on the exact bisections that split at the rank-21 and
rank-25--28 controls with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_r17_multisection_diversity.py \
  --control-calibration-only
```

This requires the preceding global degree-two comparison and the complete
bisection-specialization control artifact.  It writes
[`artifacts/generated-results/elkies-k3-r17-bisection-control-diversity-calibration-v1.json`](artifacts/generated-results/elkies-k3-r17-bisection-control-diversity-calibration-v1.json),
with SHA-256
`57945b3431cec5227ff245ae0fa238576529f23ec95bedc9f2667b7d51bc742e`.
Every induced quotient-distance distribution, graph component, clique count,
lattice span, equation weight, and displayed exceptional quotient class is
exact.  The five-control correlations are descriptive and are not a fitted
rank predictor.

The representative historical random high-rank-frame source scout is

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/hunt_lattice_foundry_rootful_source.sage \
  --ns-id NS0005 --target-frame-id NS0005-F001 \
  --generations 12 --beam 12 --samples-per-parent 60 \
  --primes 3,7,11,13,17,23 --seed 20262906 \
  --target-root-rank 15 --allow-below-target \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-mw2-source-from-high-mw-scout-v1.json \
  --root-adapted-frame-output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-mw2-source-from-high-mw-scout-root-adapted.txt
```

It starts at MW15 and emits the exact best bounded source even though the
MW0--2 target is missed; it finds `A1+2A3+A6/MW4`. This is a Kneser discovery
walk, not a certified physical elliptic-neighbour route, and is retained for
bounded provenance/replay rather than as the default foundry workflow.

Classify all distinct biquadratic pair bases, build the complete exact
5,566-row immediate-point arithmetic catalogue, replay the completed bounded
rank-lower-bound ledger, and verify the simplest rank-at-least-nine base:

```bash
/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/analyze_elkies_2026_bisection_pair_covers.sage

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/catalogue_elkies_2026_immediate_point_pairs.sage

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/screen_elkies_2026_immediate_pair_ranks.sage \
  --start 1 --limit 5566 --backend pari-only --pari-effort 1 \
  --prime-bound 200 --checkpoint-every 100

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/verify_elkies_2026_rank19_rank9_base.sage

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/search_elkies_2026_rank9_paired_base.sage

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/analyze_elkies_2026_high_rank_control_pair_bases.sage

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/sieve_elkies_2026_rank9_paired_base.sage \
  --height-bound 60 --certify all

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/sieve_elkies_2026_rank9_paired_base.sage \
  --height-floor 60 --height-bound 150 --gram-scale 10000000 \
  --certify split \
  --output artifacts/generated-results/elkies-2026-rank9-paired-base-sieve-height150.json

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/catalogue_elkies_2026_control_pair_bases.sage

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/screen_elkies_2026_immediate_pair_ranks.sage \
  --input artifacts/generated-results/elkies-2026-control-pair-base-catalogue.json \
  --output artifacts/generated-results/elkies-2026-control-pair-base-rank-ledger.json \
  --start 1 --limit 300 --backend pari-only --pari-effort 2 \
  --search-timeout 30 --checkpoint-every 10

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/search_elkies_2026_control_pair_base_points.sage \
  --certify all

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/search_elkies_2026_control_pair_base_points.sage \
  --coefficient-radius 2 --certify none \
  --output artifacts/generated-results/elkies-2026-control-pair-base-point-search-radius2.json

.venv/bin/python \
  elkies-k3/scripts/select_elkies_2026_control_spanning_bisections.py

.venv/bin/python \
  elkies-k3/scripts/certify_elkies_2026_control_pair_conductors.py \
  --limit 1 --tate-timeout 300 --checkpoint-every 1
```

The complete geometry proves all 765,167,640 distinct pair bases genus one
with generic surface rank at least 19. The bounded point ledger is complete
but does not turn empty searches into rank-zero claims. Its record lower bound
is 9, attained by two bases. The standalone verifier replays the simpler one,
constructs the paired base's 2-isogenous minimal model, and certifies an exact
birational map from nine independent points to nine rational `(t,u,v)` values.
The promoted-base search then performs bounded height-60 enumeration with a
separate `t(P)`-complexity ordering, specializes all nineteen visible points,
tests the remaining 39,118 covers at each of 100 retained parameters, and
applies Nagao only after the explicit split count. The control-pair replay
constructs all 19 pairs selected by the four high-rank fibres and computes
their finite-quotient incidence with the public exceptional complements.
The modular promoted-base sieve exhausts all 1,640 height-at-most-60
parameters with a rank-at-least-19 certificate at every fibre, then tests the
99,200 new parameters in the shell `60 < h <= 150`; neither shell contains an
extra split.  The `t=3/8` catalogue contains all 300 control-selected pairs.
Its complete bounded rank ledger certifies positive rank for 251 bases and
rank at least 6 for two.  Exact degree-two isogeny and pointed-quartic inverses
map the radius-one LLL boxes to 6,676 distinct parameters.  Every one has
exactly two split bisections and a certified surface rank lower bound 19.
Nagao scoring is performed only after those exact split counts and
certificates.
The radius-two discovery box expands this to 75,504 exact parameters and
1,144,616 post-sieve integer square tests; it also contains no third split.
Because this larger artifact uses discovery policy, it does not assert new
specialized-rank lower bounds at the no-extra-split fibres.
The dependency-free spanning selector exhausts the 619 minimum-cardinality
four-cover bases of the exceptional quotient and selects masks
`19735,22912,30787,66034` under its declared equation/priority-rank
complexity order.  It also writes exact rational parametrizations of the four
conics and retains their lifted-section formulas.  The conductor gate consumes
pair-base points in increasing projective height of the resulting `t`, proves
the generic seventeen plus the two defining cover points independent, and
only then calls exact global minimalization and local Tate reduction.  On the
declared one-candidate prefix the 81-bit-height parameter has certified rank
at least 19, while PARI times out after 300 seconds before producing a complete
conductor.  The artifact is therefore fail-closed and authorizes no wider
point search; the timeout is not a conductor estimate.

<!-- status-consumer: EC-K3-BISECT-BIQUADRATIC-R19 707bffd8b85f8f3e -->

Replay the strict target-only held-out gate for product `27431:92937`.  This
uses twelve pairwise-disjoint eight-prime blocks above 491, extends through
1151, and does not score any other product:

<!-- status-consumer: EC-K3-BISECT-MULTIQUADRATIC-CHARACTERS dc58103d8d2494cf -->

```bash
.venv/bin/python \
  elkies-k3/scripts/screen_elkies_2026_quadratic_twist_ranks.py \
  --products-only --product-key 27431:92937 \
  --prime-block 499,503,509,521,523,541,547,557 \
  --prime-block 563,569,571,577,587,593,599,601 \
  --prime-block 607,613,617,619,631,641,643,647 \
  --prime-block 653,659,661,673,677,683,691,701 \
  --prime-block 709,719,727,733,739,743,751,757 \
  --prime-block 761,769,773,787,797,809,811,821 \
  --prime-block 823,827,829,839,853,857,859,863 \
  --prime-block 877,881,883,887,907,911,919,929 \
  --prime-block 937,941,947,953,967,971,977,983 \
  --prime-block 991,997,1009,1013,1019,1021,1031,1033 \
  --prime-block 1039,1049,1051,1061,1063,1069,1087,1091 \
  --prime-block 1093,1097,1103,1109,1117,1123,1129,1151 \
  --output artifacts/generated-results/\
elkies-2026-quadratic-twist-product-27431-92937-holdout-p499-1151.json
```

The replay gives weakest/mean scores `-0.494/0.017` on `499--821` and
`-0.077/0.150` on `823--1151`.  This fails the first gate, so no `chi=4`
section solve or denominator layer is authorized.  It is a bounded heuristic
rejection, not a product-twist rank-zero theorem.

Separately, replay the broader discovery census, descend the known singleton
twist section for mask 18075, and solve its complete reduced `P.O=0`
polynomial-section scheme over `F_37`:

```bash
.venv/bin/python \
  elkies-k3/scripts/screen_elkies_2026_quadratic_twist_ranks.py

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/derive_elkies_2026_singleton_twist_section.sage \
  --mask 18075 --prime 37

/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage \
  --singleton-mask 18075 --prime 37

.venv/bin/python \
  elkies-k3/scripts/run_elkies_2026_twist_polynomial_sections_msolve.py \
  --export artifacts/local/elkies-k3/twist-polynomial-sections/\
singleton-18075/p37/export.json \
  --threads 4 --jobs 2 --timeout 180
```

The default census and the separate broader holdout artifacts are heuristic;
they are not part of the strict target-only gate above.  Their
`reproducing_command` fields contain argv without the Python interpreter, so
use `.venv/bin/python` before those stored strings.  The modular singleton
solve is exact in the displayed finite field and polynomial degree box: 23
distinct systems are empty and the sole degree-one system is the known
`Q/-Q` pair.  It is not a characteristic-zero rank upper bound.  See
[`elkies-k3/QUADRATIC_TWIST_RANK_CENSUS_2026-08-31.md`](elkies-k3/QUADRATIC_TWIST_RANK_CENSUS_2026-08-31.md).

The separate norm-twelve `0x103b2` audit descends its known anti-invariant
section and exhausts the `P.O=0` polynomial shell without Magma. The `p=19`
run takes about one minute; the meet-in-the-middle `p=29` run takes several
minutes.

```bash
sage -python elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage \
  --genus-one-label norm12-orbit-103b2 --prime 37
sage -python elkies-k3/scripts/derive_r17_genus_one_bisection_twist_section.sage \
  --prime 37 \
  --export artifacts/local/elkies-k3/twist-polynomial-sections/\
genus-one-norm12-orbit-103b2/p37/export.json \
  --output artifacts/generated-results/\
elkies-k3-norm12-orbit-103b2-twist-section-v1.json

sage -python elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage \
  --genus-one-label norm12-orbit-103b2 --prime 19
sage -python elkies-k3/scripts/derive_r17_genus_one_bisection_twist_section.sage \
  --prime 19 \
  --export artifacts/local/elkies-k3/twist-polynomial-sections/\
genus-one-norm12-orbit-103b2/p19/export.json \
  --output artifacts/generated-results/\
elkies-k3-norm12-orbit-103b2-p19-twist-section-v1.json
python3 elkies-k3/scripts/run_twist_polynomial_sections_bruteforce.py \
  --export artifacts/local/elkies-k3/twist-polynomial-sections/\
genus-one-norm12-orbit-103b2/p19/export.json \
  --known-section artifacts/generated-results/\
elkies-k3-norm12-orbit-103b2-p19-twist-section-v1.json \
  --output artifacts/generated-results/\
elkies-k3-norm12-orbit-103b2-p19-polynomial-section-bruteforce-v1.json

sage -python elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage \
  --genus-one-label norm12-orbit-103b2 --prime 29
sage -python elkies-k3/scripts/derive_r17_genus_one_bisection_twist_section.sage \
  --prime 29 \
  --export artifacts/local/elkies-k3/twist-polynomial-sections/\
genus-one-norm12-orbit-103b2/p29/export.json \
  --output artifacts/local/elkies-k3/twist-polynomial-sections/\
genus-one-norm12-orbit-103b2/p29/known.json
python3 elkies-k3/scripts/run_twist_polynomial_sections_bruteforce.py \
  --export artifacts/local/elkies-k3/twist-polynomial-sections/\
genus-one-norm12-orbit-103b2/p29/export.json \
  --known-section artifacts/local/elkies-k3/twist-polynomial-sections/\
genus-one-norm12-orbit-103b2/p29/known.json \
  --output artifacts/generated-results/\
elkies-k3-norm12-orbit-103b2-p29-polynomial-section-bruteforce-v1.json
sage -python elkies-k3/scripts/hensel_lift_r17_103b2_po0_sections.sage \
  --precision 800
```

The `p=19` shell has only the known signed pair. The `p=29` shell has 276
unsigned solutions and seven full-tangent-rank branches; only the known
branch converges and reconstructs exactly in the reduced-branch Hensel audit.
The 269 singular branches and all sections with positive `P.O` remain open,
so these computations do not determine the full twist rank. See
[`elkies-k3/R17_103B2_ANTI_INVARIANT_RANK_AUDIT_2026-09-03.md`](elkies-k3/R17_103B2_ANTI_INVARIANT_RANK_AUDIT_2026-09-03.md).

The equation-friendlier alternate q80 q6 endpoint has a distinct rootless
rank-17 lattice.  Its short shell is too large for PARI's materialized vector
list, so the same enumerator uses an LLL-reduced streaming traversal with
exact leaf norms and a PARI count-only cross-check.  It finds 39,147
section-nonnegative bisection orbits and 805,466 unoriented minimal
representatives.  This second quotient remains lattice-only: its available
equation is finite-field and no characteristic-zero branch cover is claimed.

```bash
sage -python elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage \
  --frame-artifact artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json
```

<!-- status-consumer: EC-K3-ALT-BISECT-ORBIT eca5fc0bfee5038d -->

The exact squareclass and collision-height gate has a synthetic regression:

```bash
.venv/bin/python elkies-k3/scripts/hash_bisection_extensions.py --self-test \
  --require-collision-heights --require-rank-at-least 19
```

It rejects split or non-rational-bisection covers, incomplete orbit tables,
duplicate declared translation-orbit records even with inconsistent branch
data, missing anti-invariant height data, and indefinite declared height
lattices.
Its optional complete-coverage mode also accepts the alternate q80 table when
the input declares its rootless-frame artifact and uses
`alternate_rank17_w`.  The self-test exercises that alternate-frame schema
with a temporary exact frame, while the actual q80 frame enumeration remains
the separate Sage replay above.
For a collision with smooth branch fibres and rootless double pullback, it can
instead compute the anti-invariant height form from exact lifted-section
intersections using `2*(P_i.tau(P_j)-P_i.P_j)`, checking the diagonal values
`P_i^2=-4` and `P_i.tau(P_i)=2`.  The synthetic test itself supplies no new
rootless cover; the production 39,120-cover input is the batch above.

<!-- status-consumer: EC-K3-BISECT-EXTENSION-PROTOCOL 90dc72ea57ae22dc -->

## Fast structural check

```bash
make check
```

This compiles the active Python code, checks local Markdown links, and audits
the single status ledger.

## Characteristic-two plane Keller counterexample

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_characteristic_two.py
```

The first block independently replays the Huq--Kuruvilla threefold Jacobian,
collision, inverse cubic, projective normalization charts, wild-radicial
boundary, reconstruction, and zero-pole determinant ledger.  The second block
internally replays Mondello's external theorem
[arXiv:2608.02634v1](https://arxiv.org/abs/2608.02634): the source and target
coordinate changes, skew product and preserved plane fibre, plane Jacobian,
three-point collision, hidden cubic, rational recovery identities,
irreducibility coprimality certificate, and separability witness.  The
degree-one-in-the-actual-target-parameter irreducibility argument in the
canonical note is not inferred from a bounded search.  The same argument is
written separately for arbitrary characteristic-two base fields as
`HKM2-ALLFIELDS`; no perfectness is assumed.  The command also checks that
neither displayed integer formula is a characteristic-zero Keller map.  This
internal replay is distinct from the Lean kernel-check and computationally
separate Harmonic Aristotle replay recorded by Mondello; none is independent
human peer review.

## Characteristic-two plane normalization and wild boundary

Requires Singular:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_boundary.py
```

The dependency-light block verifies the hidden-cubic discriminant, the
normalized source presentation, three compatible reconstruction charts, the
global normalized formula for the second source coordinate, the retained and
missing primes over `Q=0`, their three reduced intersections, and the local
equation giving generic different exponent one.  The Singular block certifies
the integral closure, the primitive-order conductor `(P,T)`, the two upstairs
conductor branches, and the exact reconstruction-boundary ideal.  This is a
finite exact calculation, not a bounded search.  The checker assumes a
working Singular installation in addition to the repository Python
environment.  The SymPy identities and Singular certificate do not constitute
a second independent implementation of the normalization algorithm; the
separate audit requested in the canonical note remains open.

## Characteristic-two plane modulo-four lift obstruction

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_w2_obstruction.py
```

This exact symbolic check computes the full cokernel of the first Jacobian
variation as `H^2_dR(F_2[x,y]) = xy F_2[x^2,y^2]` and reduces the integral
error to `xy(1+x^3(1+xy))^2`.  It checks the exterior-form identity, the
monomial quotient, the nonzero Cartier class, a dense
representative-independence regression, and exact regressions for affine and
triangular source generators.  The
written Jung--van der Kulk argument proves invariance under every polynomial
plane left--right equivalence; this is not a bounded correction or
equivalence search.  The command also verifies the stabilized primitive,
the explicit determinant-one lift after adjoining one identity coordinate,
the sharp degree-18 lower bound at `W_2` using the all-degree functional
`[x^13 y^4]+[x^14 y^5]`, the improved degree-25 construction over `W_3`, and
the degree-34 leading-error certificate proving that 25 is sharp in the
canonical first-digit gauge, and the geometric-series identity producing a
compatible tower over all finite Witt levels.

## Unrestricted stabilized `W_3` degree-18 exclusion

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py
```

This exact finite check expands every coefficient of an arbitrary degree-18
first Witt correction that can affect the high-degree `z^0` or `z^1`
determinant layers.  It constructs 1,639 necessary Boolean equations in
1,083 variables.  The pinned Z3 4.15.3 solver returns `unsat`, proving
`d_3 >= 19`.  This is not a sparse-support search.

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 300000
```

The second command solves all determinant equations in a degree-19 ansatz
whose first correction is linear in `z`.  Its 2,685-variable,
4,513-constraint system is satisfiable.  The checker then reconstructs the
three polynomials separately from those equations and directly recomputes
their full Jacobian as an odd constant modulo 8.  Scaling one target
coordinate by that constant gives determinant one without changing degree or
reduction, since every odd residue squares to one modulo 8.
Thus `d_3=19`; add `--show-model` to print every support.  The ansatz is used
only for the existence witness.  Neither calculation has an independent
implementation or external human review.

The preferred 440-term determinant-one witness replays without SAT:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --replay-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280_second_160.json
```

Its correction support is `280+160`: first-layer counts
`p:41,q:54,c:52,a:27,b:61,r:45`, and second-layer counts
`u0:41,v0:12,w1:39,u1:28,v1:26,u2:2,v2:5,w3:7`.  The certificate has
SHA-256
`a79984550854ce01d903783156baa6f7d4720f56ec2824bd5823cd13088a5d7f`.
For its fixed first digit, the exact second-support minimum is 160, both for
an arbitrary odd constant Jacobian and with determinant one.  No global
support minimum over all first digits is claimed.

The exact affine completion-space audit is:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280_second_160.json \
  --audit-second-linear
```

It fixes the 280-term first Witt digit and reduces the remaining equations
over `F_2`, obtaining 1,485 variables, 709 nontrivial equations, rank 681,
and nullity 804.  It also verifies that no equation mixes the `z0`, `z1`, and
`z2` blocks.  Their equation/rank/nullity triples are `(275,275,335)`,
`(205,177,203)`, and `(229,229,266)`.  This is an exact rank computation, not
by itself a support minimum.  Their incidence graphs have 59, 14, and 50
nontrivial components; the largest `z0` component has only 90 variables and
33 equations.

The generated support-reduction chain was produced with:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 60000 \
  --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness.json \
  --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_reduced.json

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 90000 \
  --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_reduced.json \
  --second-support-bound 280 --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_280.json

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 120000 \
  --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_280.json \
  --second-support-bound 180 --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_180.json

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 120000 \
  --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_180.json \
  --second-support-layer z2 --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_177.json

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 120000 \
  --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_177.json \
  --second-support-layer z1 --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_172.json

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 60000 \
  --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_172.json \
  --component-minimize-second --component-determinant one \
  --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_second_165.json

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --first-support-bound 280 --timeout-ms 120000 \
  --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280.json

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --timeout-ms 60000 \
  --minimize-second \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280.json \
  --component-minimize-second --component-determinant one \
  --write-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280_second_160.json
```

The last command proves the componentwise block minima
`(s_(2,0),s_(2,1),s_(2,2))=(92,54,14)` and hence `s_2=160`.  Replacing
`--component-determinant one` by `any` independently proves the same minimum
when any odd constant Jacobian is allowed.  Direct sparse-polynomial replay
checks that the pinned minimum-support representative already has determinant
one modulo eight.  The full nonlinear first-support search is only an upper
search: support 280 is attained, while bound 275 times out under the stated
cap, so 440 is not claimed globally minimal.

Because checked-in outputs already exist, use fresh output filenames when
rerunning the generation commands; the checker deliberately refuses to
overwrite a certificate.

## Fixed-representative stabilized `W_4` boundary

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w4_extension_obstruction.py
```

This exact sparse-polynomial checker reconstructs two `W_3` representatives
modulo 16.  For the preferred degree-19 certificate it computes the next
determinant digit, with support 1,250 and degree 54, and projects it onto
`H^3_dR(F_2[x,y,z])`.  The 48-term odd--odd--odd class is nonzero (it contains
`xyz`) and has canonical support SHA-256
`e219eeab3de8badeaf76c9cb393a1b1f0d8a791ae794ef66e97fa4b70c77b9fb`.
Thus this fixed degree-19 representative has no extension to `W_4` in any
degree.

For the block-triangular degree-25 representative, its next error `L` is
`z`-independent, has support 35 and degree 51, and has support SHA-256
`c9e4b2139db532cac4af47e976f0b2373da24996c66a3557ad00d77f61e6655d`.
The reusable affine decoder compiles the complete necessary degree-51 `z0`
system (4,082 variables and 1,677 coefficient equations) and extracts the
two-target dual certificate

```text
[(39, 12), (40, 13)]
```

with SHA-256
`6b5bb61d2fc7e79fce3d1623eef931c4a398dfa70376caf67afdc9ab50acd280`.
It is the all-degree identity
`Lambda_4(D_F(R,S))=[x^41 y^13]R`; together with `Lambda_4(L)=1`, it excludes
every constant-Jacobian extension of degree at most 51.  The checker then
directly verifies that adding `8*z*L` to the third coordinate has degree 52
and determinant one modulo 16.  Hence 52 is exact for extensions of this
fixed degree-25 representative, while the unrestricted conclusion is only
`19 <= d_4 <= 52`.

The next-class compiler and the attempted unrestricted UNSAT search are:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --audit-w4-class-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280_second_160.json \
  --timeout-ms 120000

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-class-zero \
  --w4-fix-first-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_witness_first_280_second_160.json \
  --timeout-ms 120000 --random-seed 0

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-class-zero \
  --timeout-ms 300000 --random-seed 0
```

The compiler collects 38,760 symbolic terms into 241 odd--odd--odd
coefficient equations and reproduces the preferred 48-term class exactly.
With the preferred first digit fixed, the 5,954-constraint augmented ansatz
is `unsat`.  Unrestricted within the degree-19 existence ansatz, however, the
4,754-constraint system is `sat`, so the proposed cohomological UNSAT theorem
is false.  The pinned model has support `327+491=818`, constant determinant 5
modulo 8, and certificate SHA-256
`9ad15068593af7cca87169c25eed2ff53068cc466d183ce320c2fc7d0e2c1aaa`.

The exact degree-19 joint master--subproblem experiment is:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-degree 19 \
  --w4-seed-cuts 32 --w4-cut-limit 64 \
  --timeout-ms 600000 --random-seed 0

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-degree 19 \
  --w4-seed-cuts 599 --w4-cut-limit 599 \
  --timeout-ms 600000 --random-seed 0
```

For every master model, the checker directly reconstructs the next
determinant digit, eliminates the full degree-19 affine correction system,
extracts a replayable dual certificate when it is inconsistent, and compiles
that coefficient parity back into the master.  In the 64-cut run every dual
used is a singleton zero row of the correction operator; the master remains
SAT, but its returned model again has an UNSAT completion subproblem.  The
model-independent structural audit enumerates 5,396 possible determinant
targets and 4,340 singleton holes, with layer counts
`z0:1223,z1:1228,z2:1165,z3:724` and structural-hole SHA-256
`de442207ad627a8202168496c37fcd2b9af7bb8cf03cbeb96bf90a662097ab99`.
The 599 seeded singleton targets have SHA-256
`8faad19cd5212c598c270233f1d4407791cc75f00a689f768f690c4267c2bcb1`
and compile from 6,095,343 determinant terms.  The second command returns
`unknown` at the stated solver bound.  It is deliberately an experiment and
does not prove UNSAT or change `19 <= d_4 <= 52`.

Use the shared-minor quotient compiler for the scalable version:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-degree 19 \
  --w4-seed-cuts 599 --w4-cut-limit 599 \
  --w4-quotient-compiler --timeout-ms 600000

.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-degree 19 \
  --w4-seed-source structural \
  --w4-seed-cuts 4340 --w4-cut-limit 4340 \
  --w4-quotient-compiler --timeout-ms 600000
```

The first command represents the 599 targets using 13,804 shared minor
coefficients, 899,303 minor products, and 623,490 final products.  The second
imposes the complete structural singleton quotient using 15,564 shared
minors, 926,438 minor products, and 4,200,036 final products.  Both master
decisions are `unknown` at the stated bound.  The factored coefficients are
replayed against independent direct determinant expansion on both pinned
`W_3` representatives before solving.

The structural quotient can also be split by output `z`-degree.  Use
`--w4-structural-z-layer K` and set both cut counts to the layer size
`1223,1228,1165,724` for `K=0,1,2,3`, respectively.  At 600 seconds the
recorded results are `z0:unknown`, `z1:sat`, `z2:unknown`, and `z3:sat`.
The SAT models in the odd layers still have inconsistent full degree-19
completion systems.  Adding `--solve-equations-first` to the `z0` command
also returns `unknown`; it is a solver experiment, not an elimination proof.

To export a smaller seeded system for a separate SAT solver, use an
untracked temporary path:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --degree 19 --seek-lift --require-w4-degree 19 \
  --w4-seed-cuts 8 --w4-cut-limit 8 --timeout-ms 600000 \
  --write-w4-dimacs /tmp/hkm2-w4-degree19-seed8.cnf
```

This export has 465,654 variables and 2,431,415 clauses before the external
solver's own preprocessing.  MiniSat 2.2 remained indeterminate after
409.63 CPU seconds and 3,364,085 conflicts in the recorded trial.

Replay it without SAT using:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w3_degree19_lower.py \
  --replay-certificate \
  artifacts/generated-results/huq_kuruvilla_w3_degree19_w4_class_zero.json
```

Its exact fixed-representative extension degree is computed and replayed by:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_w4_extension_obstruction.py \
  --search-degree19-extension

.venv/bin/python scripts/verify_huq_kuruvilla_w4_extension_obstruction.py \
  --replay-extension-certificate \
  artifacts/generated-results/huq_kuruvilla_w4_degree52_from_degree19.json
```

The full next error has support 1,027 and degree 52.  The unrestricted affine
correction system for this fixed map is `unsat` at degree 51 (10,255
variables and 5,791 equations) and `sat` at degree 52 (10,663 variables and
5,972 equations).  The pinned correction has supports `R:314,S:98,T:674`
and directly replays to constant determinant 13 modulo 16.  Its certificate
SHA-256 is
`438a189da33fdd081f61f9410186ca7d1b22c454dfb2cae5fbc02060f1b838ae`.
This proves exact degree 52 only after fixing the class-zero degree-19
representative; unrestricted `d_4` remains between 19 and 52.

## Plane wild-boundary atlas

```bash
.venv/bin/python scripts/verify_wild_boundary_atlas.py
.venv/bin/python scripts/verify_wild_boundary_atlas.py --singular
.venv/bin/python scripts/verify_wild_boundary_atlas.py --balanced-singular
.venv/bin/python scripts/verify_wild_boundary_atlas.py --thickened-singular
.venv/bin/python scripts/compile_plane_wild_boundary_survivors.py
.venv/bin/python scripts/verify_plane_wild_boundary_p3_degree7.py
.venv/bin/python scripts/verify_conductor_jet_truncation.py
```

The first command verifies the characteristic-divisible hidden-order
identities, the split `(f_sep,f_insep)` degree/different rows,
Artin--Schreier--Witt different formula, and the exact
numerical-semigroup conductor of the normalized missing boundary.  The second
also normalizes the seven rows
`(p,N)=(2,2),(2,4),(2,6),(3,3),(3,6),(5,5),(7,7)`, checks
the generator `P^N/T`, computes primitive-order conductor `(P,T)`, decomposes
the pullback of `P=0`, and exactly saturates the relative Kähler different
away from the fierce boundary.  In the requested rows `p=3,5,7`, the residual
ideals are `(P,W,T)`, `(P,W,T^3)`, and `(P,W,T^5)`, proving that a companion
tame branch remains with different lengths `1,3,5`.  The all-prime
normalization theorem and the uniform local row `(N-1,1,1,N-2,N-1)` are
proved separately by the `S_2+R_1` and DVR arguments in
[`PLANE_WILD_BOUNDARY_ATLAS.md`](extended-geometry/PLANE_WILD_BOUNDARY_ATLAS.md);
the bounded Singular rows are exact global regressions.  This constructs
finite plane covers, not new polynomial Keller maps: deleting only the fierce
boundary leaves the companion ramification for every `N>2`.

The third command replaces `PQT` by the balanced gluing `P^(N-1)*Q*T` and
normalizes `(p,N)=(2,2),(2,6),(3,3),(3,6),(5,5),(7,7)`.  It computes
conductors `(P,T)^(N-1)`, verifies the
base-field counts `#C=p^2+p`, `#E=p`, and `#(C-E)=p^2`, and performs the full
relative-Fitting saturation in the `p=3` row.  The written Newton-polygon and
purity argument proves uniformly that the normalized complement is étale over
the target.  The dependency-light check also proves that the natural
birational affine-plane chart, using the polynomial quotient
`(r*u^N-1)/x^(N-1)`, has Jacobian `-u^(2*N-4)`, so it is Keller only for
`N=2`.  Finally, the divisor-class localization on the Laurent chart
`D(x)=Spec(k[x^(+/-1),u^(+/-1)])` gives `Cl(C-E)=Z/(N-1)`.  Hence the balanced
complements are not affine planes for `N>2`; the matching point counts simply
show that this obstruction is invisible to the coarsest arithmetic test.
The exact named-class calculation gives `ord([L1])=N-1`, including order
five in both `N=6` controls.
The same arithmetic audit checks the full monomial band `P^a*Q*T`,
`1<=a<=N-1`: its companion row has tame index
`(N-1)/gcd(a,N-1)` and zero different only at `a=N-1`, exactly where the
class-group obstruction applies.
The written extension closes every remaining nonnegative exponent: `a=0`
has affine-UFD-core row `(N+1,N+1)`, hence one free unit and class group
`Z/(N+1)`, while `a>=N` has a wild index-`N` branch over `P=0`.  For a
general coefficient `C(P)`, the identity
`(T-1)*H_T-H=-C(P)*Q` proves that every factor away from `P=0` contributes an
additional different divisor.  At a normalized noncollision prime its exact
coefficient is `h*ord_D(R)`, not automatically the factor multiplicity `h`.
Thus the audit exhausts all one-variable gluing coefficients with no extra
target support.  Replacing `C(P)` by arbitrary `C(P,Q)` reduces the same
support-constrained search to the two-parameter monomials `cP^aQ^b`, and the
base-change theorem closes the whole quadrant.  If `b+1=p^s*d`, `d>1` is
excluded by the compactly supported Euler characteristic
`(N+1)*d-N`; for pure `p`-powers and `N>2`, finite push--pull preserves the
exact order `N-1` of `[L1]`.  In the remaining `p=N=2` tower, the corrected
core is `D(x*u)`, not `D(x)`, and its three double fibres have exact vertical
relation lattice with Smith diagonal `(1,2,2)`.  The explicit geometric
generic-fibre parameterization `R=s^2+s, y=s^(2^h)` proves its class group is
zero, so localization upgrades the vertical calculation to
`Cl(D(x*u))=(Z/2)^2` for every Frobenius exponent.  The complement of this
core in the source has two reduced primes, and the core units `x/u,u` have
unimodular valuation matrix `[[1,0],[-1,1]]`; hence restriction is an
isomorphism and `Cl(U_(2,c))=Cl(D(x*u))=(Z/2)^2`.  The dependency-light
checker replays both the generic parameterization and this source-fill
matrix.  The fourth
command supplies exact Singular normalization, conductor, and point-count
controls for `(p,N,b+1)=(2,2,2),(2,2,4),(2,2,8),(3,3,2),(3,3,3),(5,5,2)`.
Together with the unchanged generic `P=0` rows, this makes
`N=2,a=1,b=0` the unique affine-plane Keller row under the stated boundary
support hypothesis.

The fifth command checks the pinned
[`plane_wild_boundary_survivor_atlas.json`](artifacts/generated-results/plane_wild_boundary_survivor_atlas.json).
It compiles `46` proved monomial controls, `23` prescribed-degree covers,
`23` balanced prescribed-degree rows, and `25` additive/Artin--Schreier/Witt/
Kummer comparisons without promoting local or reconstruction-free rows to
Keller candidates.  The identity
`A*H_T-A'*H=P^(N-1)*Q*(A-T*A')` forces
`A=a0+T*B(T^p)` in the balanced retained-polynomial search.  The normalization
chart then proves
`#(C_A-E_A)(F_q)=q^2+(n_q(A)-1)q` and
`chi_c(C_A-E_A)=deg(A)` after splitting `A`.  Consequently all five former
support-only rows through degree `15`,
`(3,7),(3,10),(3,13),(5,11),(7,15)`, are geometrically obstructed and the
balanced reconstruction queue is empty.  The bounded packet scan rejects
`186/240`, `103/142`, and `43/77` packets in characteristics `3,5,7`; its
remaining packets are exactly pairwise coprime.  Regenerate the artifact
only after an intentional theorem or bound change with
`.venv/bin/python scripts/compile_plane_wild_boundary_survivors.py --write`.

The sixth command exhausts the six monic support-admissible retained
polynomials `A=T^4+bT+a0` over `F_3`.  Singular proves every normalization
smooth with conductor `(P,T)^2`, fierce boundary `A^1`, and relative
different supported only on that boundary.  Exact prime-field counts leave
precisely `(a0,b)=(1,1),(1,2)` with nine points on the open.  The command
then normalizes those two rows over `F_9` and `F_27`: their open counts are
`81` and `810`, respectively.  Since `810!=27^2`, neither is geometrically
an affine plane.  Regenerate its pinned artifact only with
`.venv/bin/python scripts/verify_plane_wild_boundary_p3_degree7.py --write`.

The program-wide normal-core lattice gate is replayed separately by:

```bash
.venv/bin/python plane-jc/cas/boundary_lattice_prefilter.py
python3 scripts/verify_boundary_package_compiler.py
.venv/bin/python scripts/verify_toroidal_boundary_feasibility.py
.venv/bin/python scripts/verify_f20_colored_cox_packets.py \
  --output artifacts/generated-results/f20_colored_cox_packets.json
.venv/bin/python scripts/verify_f20_global_multi_rees_cox_algebra.py \
  --output artifacts/generated-results/f20_global_multi_rees_cox_algebra.json
.venv/bin/python scripts/verify_f20_normalized_cox_conductor.py \
  --output artifacts/generated-results/f20_normalized_cox_conductor.json
.venv/bin/python scripts/verify_f20_exceptional_cox_atlas.py \
  --output artifacts/generated-results/f20_exceptional_cox_atlas.json
.venv/bin/python scripts/verify_f20_exceptional_cox_corners.py \
  --output artifacts/generated-results/f20_exceptional_cox_corners.json
.venv/bin/python scripts/verify_f20_strict_boundary_attachments.py \
  --output artifacts/generated-results/f20_strict_boundary_attachments.json
```

Besides the class-trivial-core and balanced-wild rows, these commands check
the presented-core block `[[V,A],[0,R]]`.  With `V=R=(2)`, changing only the
lift correction from `A=(0)` to `A=(1)` changes the Smith group from
`Z/2 + Z/2` to the nonsplit `Z/4`, and the compiler verifies exact order four
for the lifted named class.  The prefilter also replays the unequal
multiple-fibre formula: `(4,6,9)` forces `Z/6`, whereas the pairwise-coprime
packet `(4,9,25)` forces no vertical torsion.
The package compiler also replays the shared retained-root Euler gate:
certified degree four is obstructed by `chi_c=4`, degree one passes this
gate, and incomplete proof data remains `uncertified` rather than being
rejected.
The toroidal command compiles the \(A_4\), \(D_5\), Davenport, and corrected
Lecacheux \(F_{20}\) boundary data as colored fans and valuation matrices.
It checks smooth-cone Smith diagonals, both exceptional \(D_5\) rows, the
primitive bounded \(D_5\) model, and the Davenport unimodular block.  For
\(F_{20}\) it independently factors the discriminant as
\(d^3q^2r^2/256\), distinguishes the unramified \(q\)-crossing from the
ramified \(d\)- and \(r\)-colors, and verifies all ten generic geometric
rows.  Exact elimination then reduces the finite base-incidence locus to one
node, one ramphoid \(A_4\) cusp, two conjugate triple tangencies, one
transverse intersection, and three further tangencies.  Blowing up the node
adds four unramified derivative-order-one slope colors and one simple color,
all checked from the exact residual quadratic; the compiler therefore has
fifteen \(F_{20}\) rows.  Four Newton residuals over the ramphoid-cusp
resolution then prove the profiles \((5)\), \((5)\), \((1,2,2)\), and
\((1,1,1,1,1)\), adding ten more rows and recovering derivative sums
\((4,8,10,20)\).  The two conjugate triple centers add fourteen exact rows.
For the last cubic orbit, one residue-field computation proves the center
profile \((3,2)\), the first-exceptional profile \((2,1,2)\), and five
unramified colors on the second exceptional; its eight-row template repeats
at three centers and adds twenty-four rows.  The final \(F_{20}\) fixture has
six primitive rays, six smooth maximal cones, and sixty-three color rows.
The same checker then normalizes the \(q\)-curve by
\(t=(y^2-9)/8,\ s=4(y+2)/((y+1)(y+3))\), proves that its two crossing
slopes form the connected rational cover \(w^2=(y-5)/(y+3)\), and computes
the rank-three-to-rank-four conductor-unit pullback.  The selector \(w-1\)
completes that lattice unimodularly.  Four proposed mask columns then define
the exhaustive 1,458-assignment architecture
\(d^a q^b r^c(w-1)^e\).  It has no model, so its inverse-adjugate and
affine-recognition report sets are exactly empty.  The core fixture remains
`feasible`; the separate mask fixture is `obstructed` only within this
certificate-scoped Laurent-monomial architecture.  A genuinely colored Cox
divisor remains open.
The same run also checks the general colored divisor-span theorem.  For any
declared generator matrix \(A\) and colored target \(\tau\), it computes the
class of \(\tau\) in \(\mathbf Z^m/A\mathbf Z^n\) and emits one violated
proportional-row relation per primitive row class.  In the \(F_{20}\)
fixture, the generator rank is three, the augmented rank is four, and six
representative witnesses force color separation on the \(d\), \(q\), \(r\),
triple-\(E_1\), triple-\(E_2\), and \(q\)-\(r\)-tangent classes.
The colored-packet command constructs the first packet-supported
continuation.  Six
primitive orbit packets break all six witnesses but leave ranks `9 -> 10`.
The complete positive derivative support has sixteen packets; their
indicator columns give a saturated rank-nineteen lattice and the unique
nonnegative derivative model.  They compress to three Cartier-compatible
different-factor columns satisfying
`3*D_d+D_q+D_r=div(P_X)`, with unique model `(3,1,1)`.  The checker verifies
the determinant-minus-one \(q\)-conductor unit completion and proves that the natural `4*X-1` and
linear \(q\)-collision selectors acquire extra interior norm factors.  The
generated certificate deliberately stops at the missing conductor residue
cocycle and global Cox algebra, so it does not run entrywise adjugate or
affine-space recognition.  Regenerate the pinned JSON only with the displayed
command.
The multi-Rees command requires Singular.  It proves the general saturated
presentation for a multi-Rees algebra of two-generated ideals and applies it
to the natural (d,q,r) incidence ideals.  At a triple-(E_1) color it
certifies the local orders (4,4,2) of the base/root parameters and the
value-one gap, so the natural (d)-ideal has order two rather than the
required compact-column order one.  Exact quotient, elimination, and radical
calculations then show that (P_X\notin I_d^3I_qI_r), that the cyclic residue
has length 57 and base length 33, that (dqr) is its unique squarefree
boundary annihilator, and that its reduced base support is exactly the eight
known collision centers.  This is a normalization-first obstruction to the
ordinary incidence algebra, not an obstruction to the normalized divisorial
Cox algebra.
The normalization command also requires Singular, including `normal.lib`
and `primdec.lib`.  It computes the exact global normalization module and
conductor, decomposes the root/base residues into five collision packets,
factors the connected conductor cover, and proves that the total derivative
slope residue is anti-invariant.  Both nontrivial normalization generators
have triple-`E1` order two, and the natural normalized degree-`(3,1,1)`
product still excludes `P_X`.  The checker therefore stops at the required
exceptional Cox algebra; it does not run inverse-adjugate or affine-space
recognition.
The exceptional-atlas command uses exact SymPy arithmetic over the rational,
quadratic, and cubic chart fields.  It verifies thirteen strict-transform
types covering all `48` positive exceptional colors, constructs one
primitive value-one Rees variable on each chart, and proves
`P_X=tau^(3*a+b+c)*dF/dY` for the compact local orders `(a,b,c)`.  It also
checks a parity-compatible three-frame factorization of the total derivative
residue on the punctured conductor.  The generated certificate stops at the
missing Čech overlap maps; it does not promote chartwise cancellation to a
global Cox ring, full inverse-adjugate polynomiality, or affine-space
recognition.
The exceptional-corner command verifies the two-parameter continuation.  It
constructs four exact corners in the cusp graph and three in the
`q-r`-tangency graph, proves bivariate strict-transform and derivative
identities, and checks the compact `(3,1,1)` monomial on every corner.  It
also proves the complementary `q`-node transition law and certifies that the
positive triple-`E1` and triple-`E2` colors are root-center separated despite
base-ray adjacency.  Strict-boundary attachments and the global Čech class
remain outside the certificate.
The strict-boundary command requires Singular for one saturated incidence
calculation and otherwise uses exact SymPy arithmetic.  It proves the
weighted Taylor--Cox attachment theorem and instantiates it on the cusp
`E4`-to-`r`, triple-`E2`-to-`d`, and `q-r` `A`-packet-to-`q` families.  It
also proves that the `q`-node scaled incidence loses its root-coordinate
coefficient, so that edge remains a conductor-normalization problem.  It
locates Cartier-compatible double-root fibres for the remaining two triple
and six `q-r` strict-`r` candidate colours, without claiming their full
family saturation.  Those saturations, conductor transitions, global Čech class,
inverse-adjugate polynomiality, and affine-space recognition are not
claimed.
The boundary-package compiler command also verifies the sharp scalar
conductor/contact-loss bound
`n_i >= c_i+d_i+ell_i+epsilon_i`, its stronger per-input/output dependency
bound `n_(i,j) >= c_i+lambda_(i,alpha,j)`, node and cusp quotients, arbitrary
numerical-semigroup gap bases, the four audit states, strict scalar and
detailed JSON parsing, asymmetric `P/Q` losses, exact deficit reporting, the
certified omitted-support-frontier normal-valuation adapter, and invariance of a finite
matching-map cokernel and distinguished class.  It writes no artifact.

<!-- status-consumer: BL1 e86cdcd66993bccc -->
<!-- status-consumer: PWB1 4ce9a0bf6d277321 -->
<!-- status-consumer: PWB2 2346d64d0f1eaa07 -->
<!-- status-consumer: PWB3 a6bcf405759ddd5d -->
<!-- status-consumer: PWB4 ebddf245e65b62a7 -->
<!-- status-consumer: PWB5 142b02344181fed3 -->
<!-- status-consumer: PWB6 35636805d73e0bec -->
<!-- status-consumer: PWB7 19f4f4ffc96227a3 -->
<!-- status-consumer: CJT1 afb70f90ff10f3d7 -->

## Six-variable quartic HN Waring rigidity

Requires SymPy in the repository Python environment; the loop-closure
command also requires Singular:

```bash
.venv/bin/python scripts/verify_quartic_hn_waring_rigidity.py
.venv/bin/python scripts/verify_quartic_hn_rank9_one_zero.py
.venv/bin/python scripts/verify_quartic_hn_rank9_top_determinant.py
.venv/bin/python scripts/verify_quartic_hn_rank10_parallel_obstructions.py
.venv/bin/python scripts/verify_quartic_hn_rank10_matroid_survivors.py
.venv/bin/python scripts/verify_quartic_hn_rank10_loop_survivors.py
.venv/bin/python scripts/verify_quartic_hn_rank10_loop_closure.py
.venv/bin/python scripts/verify_quartic_hn_rank10_simple_survivors.py
```

These commands exactly replay the finite codimension-two Gale calculation,
the rank-nine one-zero trace coefficients and complementary-minor gate, the
rank-ten parallel-class obstructions, and the rank-ten nonsimple-survivor
audit.  The last check gives literal and loopless characteristic-zero
counterexamples to the proposed cyclic-complement lemma, constructs the
normalized realization ideals for the frozen catalogue slice, and excludes
all 35 characteristic-zero types from the full Gram branch using the first
two HN traces.  The
catalogue completeness is relative to `matroid-database==0.3`.  The
penultimate command replays the frozen complete Gale-loop census: 115
abstract coloured types, 111 characteristic-zero types, 103 universally
closed by splitting, the self-square support obstruction, or exact SymPy
saturation, with eight special types left.  The final Singular-backed command
closes those eight by four exact characteristic-zero saturations, one
disconnected Waring splitting, and three loop--triple Witt obstructions.
The final simple-survivor command freezes the exact 23-extension/five-type
census, constructs a rational realization of every type, and closes all five
rank-six Gram branches using the universal six-point-flat obstruction.  Thus
the rank-ten branch is empty and every essential six-variable quartic HN
counterexample has Waring rank at least eleven.  See
[`QUARTIC_HN_RANK10_MATROID_SURVIVORS.md`](extended-geometry/QUARTIC_HN_RANK10_MATROID_SURVIVORS.md).
Pass `--details` to either survivor verifier to print normalized matrices,
ideal generators, basis-minor data, and rational witnesses as exact JSON.

To regenerate the external catalogue extraction rather than replay its
pinned artifact, run

```bash
.venv/bin/python scripts/enumerate_quartic_hn_rank10_loop_survivors.py \
  --database-root /path/to/matroid_database/_all
.venv/bin/python scripts/enumerate_quartic_hn_rank10_simple_survivors.py \
  --database-root /path/to/matroid_database/_all
```

The source must be the `matroid-database==0.3` rank-four files.  The
loop enumerator uses the files through nine elements; the simple enumerator
uses the complete nine-element file.  The canonical note records the
source-wheel, catalogue-file, and generated-artifact hashes.

## LND-image Mathieu finite-fiber replay

```bash
.venv/bin/python scripts/verify_lnd_radical_slice_fibers.py
```

For the slice LND `D=d/ds`, this constructs a radical complete intersection
of six points in three vertical fibers.  It checks reducedness, verifies on
a generic degree window that primitive membership is exactly the three
vertical interval conditions, and replays a safe seed and multiplier through
exponent twelve.  It also replays two nonreduced length-two residual schemes
for the carrier `q=s^2`, one on and one off the carrier.  The powers are
bounded regression checks; the all-order arbitrary finite-residual and
monic-carrier finite-residual theorems are proved in
[`LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md`](extended-geometry/LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md).

## LND-image nonmonic degree-drop search

```bash
.venv/bin/python scripts/search_lnd_nonmonic_degree_drop.py
```

This uses exact rational arithmetic for the carrier `p=x*s-1`, powers
`x^c*p^d` with `(c,d)=(0,1),(0,2),(1,1)`, three primary residual schemes
at the degree-drop point, 256 sparse seeds, and six fixed multipliers.  The
primitive-carrier assertions replay the all-order support-weight exclusion
theorem.  All reported finite prefixes are bounded regressions or candidate
generation only.

## LND-image plinth-divisor search

```bash
.venv/bin/python scripts/search_lnd_plinth_ideal_images.py
```

This searches the linear LND `D=x*d/dy+y*d/dz` on `Q[x,y,z]`, whose local
slice `y/x` has plinth element `x`.  In each required homogeneous degree it
constructs `D(I_n)` exactly by rational linear algebra.  Five homogeneous
ideals, 45 sparse seeds, pure powers through six, and four multipliers are
tested.  The output is candidate generation only.

## LND-image reducible-plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_reducible_plinth.py
```

This searches
`D_r=x*(x-1)*d/dy+y^r*d/dz` for `r=1,2` on `Q[x,y,z]`.  With weights
`wt(x)=0`, `wt(y)=1`, and `wt(z)=r+1`, image membership in each weight is
an exact finite module calculation over `Q[x]`; Singular computes the
module standard bases.  Five ideals couple the fibers `x=0,1`.  The
`r=1` profile tests 205 seeds through weight three, and the `r=2` profile
tests 210 through weight four; both use pure powers through six and five
multipliers.  A third profile tests
`D=x*(x-1)*d/dy+(y^2+x)*d/dz` on eight ideals and 350 mixed-weight seeds.
It uses bounded normalized primitive lifts and finite quotient/kernel
images, so individual membership decisions remain exact even though the
grading is broken.  Every exponent range is bounded.

## LND-image crossing-plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_crossing_plinth.py
```

This searches the four-variable LND
`D=u*v*d/dy+(y^2+u)*d/dz`, whose plinth divisor has intersecting
components `u=0` and `v=0`.  The compiler normalizes primitives modulo
`Q[u,v,3*u*v*z-y^3-3*u*y]`, computes an exact bounded lift by a Singular
module standard basis over `Q[u,v]`, and decides the remaining kernel
correction in each finite quotient exactly.  Five zero-dimensional
crossing ideals and 956 sparse seeds are tested through six pure powers;
mixed powers four through six use the multipliers `1,u,v,y,z`.  Individual
membership decisions are exact, but both exponent windows are bounded.

## LND-image nonprincipal-plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_nonprincipal_plinth.py
```

This searches `D=u*d/dx+v*d/dy` on `Q[u,v,x,y]`, with invariant
`w=u*y-v*x` and nonprincipal plinth ideal `(u,v)` in the kernel.
Homogeneous primitive membership is an exact Singular module calculation
over `Q[u,v]`; kernel corrections are decided exactly in five finite
quotients by closure under `u,v,w`.  The search tests 1,055 sparse seeds,
pure powers through six, and mixed powers four through six for
`1,u,v,x,y,w`.  Individual membership decisions are exact; the exponent
windows are bounded.

## LND-image positive-dimensional plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_positive_dimensional_plinth.py
```

This retains a free `y`-direction for `D=u*d/dx+v*d/dy`.  The five ideals
begin with `(u,v,x)` and include four nilpotent or tilted plinth jets.
Although their quotients are positive-dimensional, the images of
`u,v,w=u*y-v*x` are nilpotent, so the kernel-image span and every
membership decision are exact without a `y`-degree cutoff.  The checker
also verifies through total degree eight the coefficient-functional
identity used in the free-line corollary of the all-order plinth-power
saturation theorem, and the filtration identity through total degree
seven.  Pure powers through six and mixed powers four through six are
bounded regressions; the written theorem independently proves all five
displayed ideals safe.

## LND-image principal-conductor search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_principal_conductor.py
```

For `D=u*d/dx+v*d/dy` and `I=(x)`, the image of the invariant ring modulo
`I` is exactly `Q[u,v,u*y]`.  Consequently a primitive residue monomial
`u^a*v^b*y^c` is correctable exactly when `a>=c`.  The script combines
this exact valuation-face test with Singular primitive lifts for 1,055
seeds, pure powers through six, and mixed powers four through six.  The
membership decisions are exact and untruncated; the exponent windows are
bounded.  The support census also checks 43 `y`-free survivors, thirteen
invariant survivors, four forms touching the slope-zero face, and exact
membership of all 48 survivors in `u*ker(D)[x]`.  It additionally verifies
the exact algebraic square-gate failure: a homogenized shifted-Legendre
face gives `f,f^2 in D((x))` but `f^3 not in D((x))`.  This is not an
LNED counterexample.  It also checks the local-slice identity
`T(D(x*y*a))=(u*y)*a` underlying the exact two-branch criterion, together
with `D(u*x*a)=u*D(x*a)` and `D(u*y*a)=u*D(y*a)` for the aligned and
crossed invariant-content carriers.  For the rational-root carrier
`q_1=u*x+w`, it normalizes `A/(q_1)` by
`x=u*t, y=-(u-v)*t`; the invariant image is
`Q[u,v,u^2*t]`, so the exact residue test is `a>=2*c`.  Seventeen seeds
survive the genuinely eventual window `m=4,5,6`; all seventeen lie in
`u*ker(D)[q_1]`, and no bounded mixed-tail obstruction occurs.  The
checker also verifies the normalization identities through the ladder
`q_n=u^n*x+w`, `1<=n<=4`.  For the first tied carrier
`q=u*v*x+w`, its exact quotient cone is `a>=2*c, b>=c`.  Eight seeds
survive powers `4,5,6`; all eight lie in `u*v*ker(D)[q]`, with no bounded
mixed-tail obstruction.  The checker verifies the two-prime
normalization identities for `1<=r,s<=3`.  Finally, for a sample
invariant-affine coordinate `h=b_0+b_1*x+b_2*y`, it verifies
`D(h)=b_1*u+b_2*v`, `D^2(h)=0`, and both inverse-chart identities
expressing `x,y` in `ker(D)_(D(h))[h]`.

The all-order proof in the canonical note uses the
full eventual-power hypothesis and the one-variable polynomial moment
lemma to prove that `D((x))`, and hence `D((ell(x,y)))` for every nonzero
linear form `ell`, is Mathieu--Zhao.  Divisibility bootstrapping closes
all powers `(ell^d)`, while primitive evaluation at two generic roots
proves zero Mathieu radical for carriers such as `(x*y)` and
`(x*(x-1))`.  A fixed-denominator extension of the local-slice proof also
closes `(a*ell)` for `a` in `{u,v}` and `ell` in `{x,y}`.  A lowest-face
moment argument closes every rational-root ladder carrier
`q_n=u^n*x+w`, `n>=1`, and the paired lowest-face argument closes the
two-prime grid `q_(r,s)=u^r*v^s*x+w`, `r,s>=1`.  None of these all-order
results is inferred from the bounded census.  The same lowest-face proof
closes `u^r*v^s*x+b` for every invariant intercept `b in ker(D)`, with
`r>=1` and `s>=0`; this arbitrary-intercept extension is a written
theorem rather than a separate bounded search.  Prime-by-prime lowest
faces close every coprime `q=a*x+b`, `a,b in ker(D)`.  On the aligned
condition `v*(b mod u)=(a mod u)*w`, one has `q in u*A`, so pure
membership itself forces the missing `u`-content.  A nontrivial common
invariant factor is removed by the Mathieu scaling lemma: if `M` is
Mathieu--Zhao and `c*M` is contained in `M`, then `c*M` is
Mathieu--Zhao.  Consequently every `q=a*x+b` with nonzero
`a in ker(D)` is closed, without a coprimality assumption.  The intrinsic
condition for generic orbit degree one is `D^2(h)=0`, `D(h)!=0`.
Prime-local inverse charts and the same moment-face argument close every
irreducible such `h`; divisibility bootstrapping closes all powers, and
the scaling lemma restores invariant content.  Factoring an arbitrary
principal carrier over the generic orbit now gives either at least two
distinct roots (zero eventual-power radical) or one invariant-affine
irreducible factor with multiplicity.  Thus the canonical note proves
`D(q*A)` Mathieu--Zhao for every nonzero `q in A` for this model
derivation.  The checker only replays the coordinate identities and
bounded searches; the all-order conclusion is deductive.

## Checked-in Lean projects

All five local Lean packages use the pinned Lean/Mathlib `v4.32.1` release.
Build their default targets and audit their source policies with:

```bash
make verify-lean-local
```

This rejects `sorry` and `admit` throughout `formal/`, rejects unexpected
explicit axioms, checks the finite-étale publication-certificate import
boundary, and builds `discriminant-pencils`, `finite-etale-keller`, `gmc2`,
`gvc`, and `support-saturation`.  The GMC(2) package deliberately exposes
exactly two mathematical
inputs as axioms; their names and roles are documented in
[`formal/gmc2/README.md`](formal/gmc2/README.md).  The GVC package contains no
explicit axioms, but is currently a partial audit.  It constructs and proves
both the minimal concrete and arbitrary cusp-profile quadric phase bridges,
giving the full profile failure family over every characteristic-zero field
and an unconditional counterexample in every finite dimension at least
three.  Its sole remaining bridge structure isolates the unformalized binary
envelope obligations: global envelope closure and the `delta = 0` equal-face
ordering from shifted-ray separation.  The complete finite-support
common-threshold cutoff is now proved in Lean: finite maximization constructs
the integral coordinate cut, which is refined to a strict positive weight
whose unit gap is amplified under powers.  Empty equality faces are handled
directly.  The finite core of Lemma 3.1 is also checked: Mathlib's Hall
theorem extracts a deficient set, two-dimensional linear algebra localizes it
in one direction class, and exact counting yields the sharp `d - e + 1`
annihilator bound.  Lean also proves the coordinate-free power divisibilities
for both displayed normal forms.  The translated
Duistermaat--van der Kallen/polarization
step that supplies the absence of a matching remains explicit.  Its
checked core now includes the algebraic beta and full endpoint-profile
coefficient calculations, the literal multivariate profile family's
degree/order formula, the rational `p`-adic factorial-valuation lemma,
the coefficientwise Reynolds/Laurent phase identity and concrete endpoint
extraction, the apolar contraction and operator-composition laws,
characteristic-zero base change, unused-variable padding, and the final
negative-slope envelope crossing.  The full-profile proof enforces the
manuscript's declared-degree condition `S.natDegree <= e`, proves arbitrary
even-phase extraction and the shifted primitive identity, and constructs
the previously missing bridge.  Its exact coverage is
documented in [`formal/gvc/README.md`](formal/gvc/README.md).

Build only the support-saturation theorem with:

```bash
make verify-lean-support-saturation
```

This kernel-checks the associated-prime equivalence, the regular-element
criterion, the no-embedded-primes support theorem, the explicit torsion
counter-witness, and the quotient/presentation saturation equivalence.  The
precise boundary of the formalization is recorded in
[`formal/support-saturation/README.md`](formal/support-saturation/README.md).

<!-- status-consumer: SST1F 838e558b5fcb9d81 -->

Build only the GVC audit with:

```bash
make verify-gvc-lean
```

## Stable core

```bash
make verify-minimal
make verify-core
make verify-foundations
```

`verify-minimal` uses only the Python standard library for the foundational
map.  `verify-core` adds the cubic marked-root and exact-image implementations.
Its normalized-factorization certificate checks both polynomial compositions
across `a=0`, residual-torus equivariance, determinant `-1` for normalized
multiplication, and the two explicit linear changes recovering the announced
map.  The same target runs the scoped ordinary-degree-six boundary audit;
it can be replayed separately with:

```bash
.venv/bin/python scripts/verify_ordinary_degree_six_boundary_audit.py
```

This verifies the exact asymmetric `(1,2)` determinant and degree floor,
both balanced `2+2` resolution charts and their transition in addition to
the cone/Cox identities, the affine-linear Wronskian reduction, the
impossibility of the `(1,2)` cubic profile, the removable-jet residue
obstruction for `(0,3)`, the rational-map-degree obstruction for a
nonconstant boundary jet, and the double-pole residue obstruction for a
nonzero constant boundary jet.  It also checks the pure-`C` weighted `D^3`
divisor ledger, the degree-eleven floor for every `z`-linear standard
reciprocal clearing, and the nodal conductor character rank.  It does not
enumerate arbitrary affine modifications, nonmonomial Wronskian profiles,
or all maps with two boundary relations.

It also checks the normalized `(2,3)` factorization slice: the unimodular
boundary lattice, class `L^5-L^3`, direct counts `q^5-q^3` for four small
prime fields, and generic degree ten.  A separate two-chart certificate
checks the Euclidean quadratic norm, its affine-modification presentation
over `A^2 x SL_2`, the complementary-chart transition, and the integral
residue coefficients `1`, `2`, and the nonzero mod-two boundary used to prove
that integral cohomology is `Z` in degrees zero and three only.  The scripts
check the algebraic, arithmetic, and Gysin inputs; the written audit supplies
the localization sequences and homotopy argument.
The same slice now has an exact invariant-kernel certificate.  It constructs
the primitive saturated LND `D7`, verifies
`ker(D10) intersect ker(D7)=k[K,H,V]`,
`ker(D7)=k[K,H,V,s]`, and `ker(D22)=k[K,H,V,W]`, and checks the boundary
identities used by the minimal-pole proof.  The finite/non-finite control
target also replays Maubach's cusp-base ladder:

```bash
make verify-hilbert14-invariants
```

The bounded ladder replay is regression evidence.  The uniform
modulo-`T^4` degree argument in
[`HILBERT14_INVARIANT_KERNEL_PROGRAM.md`](extended-geometry/HILBERT14_INVARIANT_KERNEL_PROGRAM.md)
is the non-finite-generation proof.  The same target verifies the next
normalized `(2,4)` experiment: three triangular gauge LNDs, the exact
generic quotient `M^2+4*U*N^2=256*a^4`, the regular boundary classes
`C,Q,S`, and the boundary-linear relation that terminates the saturation
ladder.  It also checks the induced third LND and the exact finitely
generated triple intersection `k[a,U,N,M,C]`, whose boundary is a cusp.
The written minimal-pole proofs are in
[`QUADRATIC_QUARTIC_HILBERT14_SLICE.md`](extended-geometry/QUADRATIC_QUARTIC_HILBERT14_SLICE.md).
Finally, the target runs the genuine multiboundary control.  It checks the
two commuting cusp LNDs, the invariant grid
`s^2*t^2*(X+sY)^m*(U+tV)^n`, and the conductor-square replay modulo
`(s^4,t^4)`.  The arbitrary-bidegree rectangle escape in
[`MULTIBOUNDARY_HILBERT14_CONTROL.md`](extended-geometry/MULTIBOUNDARY_HILBERT14_CONTROL.md)
is the non-finite-generation proof.  The same note computes the exact
finite-generation ideal as the conductor
`s^2*t^2*k[s,t,P,Q]` to the normalized-ambient invariant algebra and gives
its four infinite return ladders and infinite monomial SAGBI basis.  The
checker replays the monomial conductor criterion in a configurable box; the
written localization/specialization argument proves the arbitrary-degree
statement.  The same argument is proved for every
`tensor_i(k+t_i^2*k[t_i,P_i])`: its finite-generation ideal is the product
conductor and has `2^r` infinite return ladders.  The checker replays this
formula through `r=4` by default.  The note also proves that the two leading
divisors in a tangent-normalized factorization slice are disjoint and
explains why the coupled three-boundary Cox fill is a different branch.
The same target verifies the general weight-`(1,-1,-k)` invariant-coordinate
Jacobian reduction for `k=1,2,3,4`, including the foundational
`(-2,-1,1)` output weights.  It then reconstructs the complete
sixteen-monomial coefficient ideal, proves the gauge-fixed dual-number
presentation, extracts the infinitesimal deformation and its quadratic
obstruction, and separates it from the affine left--right orbit.  The same
target independently rewrites that normalized ideal as three univariate
weighted-Wronskian layers and checks its exact Poisson-square and tangent-pencil
identities.  The leading layer exposes the quadratic obstruction directly,
and two further unit-ideal checks eliminate both one-sided nonconstant-`C`
boundary charts.

## Collision-axis unimodular frontend

Verify the normalized one-variable reduction, the sharp three-occurrence
pure-axis bound, the exact minimum elementary length two, and the induced
pruning of the pinned degree-seven support ledger with:

```bash
.venv/bin/python scripts/verify_collision_axis_unimodular_frontend.py
```

The audit shows that the elementary orbit problem itself is trivial and that
the zero-moment constraint is not preserved by the relative elementary
action.  The useful gate is instead
`gcd(h_1',h_2',h_3')=1` for `h(t)=F(t*e_1)`.  It eliminates all balanced
supports of sizes four and five and leaves `900` labelled (`450` orbit)
supports in size six.  It does not strengthen the already proved global lower
bound seven, because completion over `k[t]` and first-jet integrability impose
no further obstruction.  See
[`COLLISION_AXIS_UNIMODULAR_FRONTEND.md`](extended-geometry/COLLISION_AXIS_UNIMODULAR_FRONTEND.md).

## Global low-degree support census below `(7,6,4)`

Generate the eight support-first stage ledgers and their manifest with:

```bash
.venv/bin/python scripts/compile_global_low_degree_census.py
```

This enumerates the 74 invariant degree flags, the complete raw-degree-seven
exact-support strata through six nonlinear monomial occurrences, every
determinant bucket, and every integer infinity weight modulo exposed Newton
faces and coordinate strata.  It then runs the sign SMT gate, exact Singular
coefficient-torus algebra over `F_11`, `F_13`, `F_17`, and `QQ`, plus an
independent SymPy rational Gröbner replay.  The pinned result has `30`, `85`,
and `1694` determinant-balanced labelled supports in sizes four, five, and
six; their `913` residual-symmetry representatives all have unit exact ideals.
The dense quadratic collision ideal is also `(1)`.

Replay every pinned JSON decision with:

```bash
make verify-global-low-degree-census
```

The result is complete only through nonlinear support six and for the dense
degree-at-most-two row.  It proves a support lower bound of seven below
`(7,6,4)`, without asserting attainment at seven;
it does not claim the cardinality-unbounded census is complete.

`verify-foundations` adds the weighted construction and its clean-room checker.
It also runs the all-degree rational-fiber checker, whose symbolic odd/even
identities prove uniform admissibility and whose exact degrees `3,...,100`
remain as a regression:

```bash
.venv/bin/python scripts/verify_padic_inverse_branches.py
.venv/bin/python scripts/verify_foundational_arithmetic_dynamics.py
.venv/bin/python scripts/verify_composite_degree_twelve.py
.venv/bin/python scripts/verify_degree_twelve_wreath_elimination.py
.venv/bin/python scripts/verify_all_degree_rational_fibers.py
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
.venv/bin/python scripts/verify_common_arithmetic_fibers.py
.venv/bin/python scripts/verify_locally_prescribed_common_fibers.py
.venv/bin/python scripts/search_cross_family_collision.py
.venv/bin/python scripts/verify_universal_quartic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_quartic_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_cubic_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_power_shifted_gauge_multiplicity.py
.venv/bin/python scripts/verify_whole_plane_stable_multiplicity.py
.venv/bin/python scripts/verify_universal_quintic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_higher_degree_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_multiplicity_witness_cards.py
.venv/bin/python scripts/verify_universal_relative_keller_map.py
.venv/bin/python scripts/verify_generic_tschirnhaus_non_descent.py
.venv/bin/python scripts/verify_rank_five_tschirnhaus_transition_locus.py
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --module-resolution
.venv/bin/python scripts/verify_keller_tschirnhaus_descent_567.py
.venv/bin/python scripts/verify_rank_three_collision_descent.py
.venv/bin/python scripts/verify_rank_four_collision_cross_ratio.py
.venv/bin/python scripts/verify_low_rank_multiplicity_boundaries.py
.venv/bin/python scripts/verify_real_fiber_spectrum.py
.venv/bin/python scripts/verify_adelic_fiber_engineering.py
.venv/bin/python scripts/verify_local_global_keller_fibers.py
.venv/bin/python scripts/verify_a5_grunwald_keller_fiber.py
.venv/bin/python scripts/verify_hasse_keller_fiber.py
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
.venv/bin/python scripts/verify_multiplicative_hasse_artifact.py
```

The finite-étale Keller-fiber checker includes the exact ordered collision
fiber. For the degree-three, degree-four, and degree-five presentations it
verifies the diagonal/off-diagonal Chinese-remainder decomposition, its
explicit separability idempotent, and ranks `N^2`, `N`, and `N*(N-1)`. It
also verifies three cubic `S_3` normal-closure sheets and the optimal Hasse
fiber decomposition `A5 tensor A5 = A5 times (N6^3 times L2)`. The
presentation-independent collision algebra, diagonal kernel, ordered-pair
universal property, and obstruction rank are checked in Lean.

The universal-relative checker applies this interface to the Osada
`T^N-T-1` root covers in ranks three through eight.  It verifies the
divided-difference idempotent and the exact collision, diagonal, and
off-diagonal standard-monomial ranks `N^2`, `N`, and `N*(N-1)`.  It also
enumerates every `S_N` orbit of ordered distinct `m`-tuples for
`1<=m<=N<=8`, checking the rank `N!/(N-m)!` and stabilizer `(N-m)!`.
The rank-three descent checker then identifies the cubic ordered-pair sheet
with the full `S_3` frame torsor and verifies the exact projective
interpolation cocycle, quadratic Tschirnhaus boundary ledger,
target-localized factorization transport, saturated global stabilizer, and
fixed-map scaling equivariance.  It does not classify nonlinear polynomial
self-equivalences outside the canonical factorization transport.
The rank-four continuation checks that ordered triples, rather than ordered
pairs, give the full `S_4` frame.  It factors the fourth-root projective
interpolation residual and labeled cross-ratio difference by the same exact
defect, separates that defect from the primitive-element boundary, and
clears it into the universal quartic Keller target coordinates.  It does not
assume that every Keller-incidence equivalence is projective on the root
line.

The all-rank continuation is:

```bash
.venv/bin/python scripts/verify_all_rank_collision_projective_descent.py
```

It completes `Conf_(N-1)` to the full `S_N` frame, verifies the intrinsic
rank-at-most-three criterion for the columns `1,r,u,r*u`, constructs its
normalized polynomial coefficient matrix, recovers the automatic cubic and
quartic cross-ratio cases, and checks the `N-3` independent framed residuals.
It also supplies exact projective and primitive-nonprojective witnesses in
every tested rank.  The bounded replay supports the written all-rank
linear-algebra proof; it does not claim that every Keller equivalence acts
projectively on the root line.

The generic stable non-descent continuation is:

```bash
.venv/bin/python scripts/verify_generic_tschirnhaus_non_descent.py
```

It verifies that the split change `r -> r+r^2` has a nonzero projective
minor in every rank, checks the quintic `I_5` base case, and proves
symbolically for every `N>=6` that

```text
J_N(P_(r+r^2))-J_N(P_r)
=-(N-1)(7N+11)/(30N(N+1)(N+2)).
```

Direct root-polynomial calculations through rank twenty are regressions.
Dominance of the presentation-to-boundary map and the resulting generic
codimension statement are written geometric proofs, not bounded searches.

The rank-five transition-locus continuation is:

```bash
.venv/bin/python scripts/verify_rank_five_tschirnhaus_transition_locus.py
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py
```

It computes the ambient stable-equivalence hypersurface and the two
labelled projective residuals, verifies local dimensions `4`, `3`, and `2`
for the ambient, projective, and intersection loci, and checks the explicit
coefficient-torus equivalence.  It also proves that this canonical
equivalence carries the selected complete fibre exactly on the root-scaling
locus.  The canonical note then determines the fixed map's standard marked
stable target orbit completely.  It does not classify vertical
automorphisms of the added identity factors.

The second command supplies the exact fixed-map calculations.  It
factors the prime quintic ramified discriminant, computes the exact
logarithmic-vector-field spaces through quotient degree twelve using FLINT
integer nullspaces, and applies a three-point characteristic-zero Jacobian
Groebner test.  An exact triple-root point forces the boundary multiplier to
be one in every degree.  Exact recursive Newton-face pruning proves that the
stable marked-target orbit is a point through total target degree twenty-eight;
every unstabilized target self-equivalence in every degree is the identity.
The all-degree unstabilized conclusion uses the coordinate-polynomial
intruder theorem at `P^2*B^5*C`.  Kuroda's stable-invariant theorem, applied
to conjugates of every stable translation by the target automorphism and its
inverse, makes the standard marked orbit a point for arbitrary stabilization.
All bounded stable
branches expose one of `P^12*C^4` and `P^2*B^5*C`.
The checker also proves in all degrees that these are the only positive
upper Newton vertices.  They tie on
`10*w_P-5*w_B+3*w_C=0`, and it verifies an explicit logarithmic field whose
two leading contributions cancel there.  More sharply, its `P`-zero Koszul
ladder first ties at target degree fifty, where the UFD cube condition
fails, and first admits leading cancellation at degree fifty-five.  Thus
unrestricted monomial avoidance is false; exact boundary preservation on
that one binomial wall remains open.  This is a comparatively expensive
exact regression.  The optional third command requires Macaulay2.  It
proves that the homogenized logarithmic module has two generators in
quotient degree seven, thirteen in degree eight, eighteen first relations
in degree nine, and six second relations in degree ten.  Equivalently, its
filtered Hilbert numerator is `2*t^7+13*t^8-18*t^9+6*t^10`.

The optional Newton-topology certificate requires Singular:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --generic-fibre-newton
```

It checks all forty-six nontrivial Newton faces of `H-h` and its coordinate
restrictions, obtains normalized-volume contributions
`8,2,0;-38,-52,-2;328`, and certifies `chi(H=h)=246`.  Thus the
vanishing-`H^2` stable-rigidity shortcut does not apply; this calculation is
independent of the Kuroda descent proof.

The optional all-degree wall research calculations are:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-homology
make verify-rank-five-singular-support
```

The singular-support target is the short form of these three commands:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-triple-root-prime
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-two-double-root-prime
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-singular-boundary
```

These use exact characteristic zero and require Macaulay2.  Add
`--research-characteristic=1000003` for the much faster good-prime
discovery pass.  The homology command proves that the non-Koszul quotient
has dimension two, degree 296, and the Betti table recorded in the canonical
note.  The next two commands construct the prime triple-root and
two-double-root curves by contraction; their projective degrees are
seventeen and nineteen.  The final command proves that the affine
`P=0` chart is empty and that the radical at infinity is `(Z,P*C)`.
Together with the root-partition argument in the canonical note, these
three targeted commands give the four minimal supports without asking
Macaulay2 for a blind primary decomposition.  The older
`--research-singular-primes` mode remains available as an expensive
independent comparison, but is not part of the proof chain.

The first exact continuation of the cancellable target-degree-55 Koszul
wall is:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-hensel --research-depth=4
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-hensel --research-depth=9 \
  --research-zero-l=l_3_0 --research-continue-constraints
```

The first command projects to `B`-degree zero and forces `l_3_0^4=0`.
After imposing `l_3_0=0`, the second uniquely lifts all divisible residuals
with lower homogeneous pieces of `G`, records the depth-eight plane
constraint, and proves that it and the five depth-nine equations generate
the unit ideal.  This excludes the normalized two-generator family
`V_B=L*H_C`, `V_C=-L*H_B-eta*G*H`.  It does not include the independent
`H`-multiple in `V_B` or any non-Koszul singular-support class.

The first nonprojective rank-four continuation is:

```bash
.venv/bin/python scripts/verify_rank_four_nonprojective_keller_lift.py
```

It verifies the exact fifth-power ground-field orbit class, shows why the
`(1,2,3,4)` quadratic witness carries a separate rational Kummer twist, and
constructs an arithmetic-neutral witness with seed ratio `4^5`.  It checks
the explicit source--target scaling to the single fixed map
`F_(-124416)`, the two endpoint fibers, and their residual target
translation.  The straight target line has inverse polynomial
`-(S-12)(S+12)(S^2+24S+108*lambda)/3456`; its discriminant and wrong framed
sheet partition are checked exactly.  A label-preserving rational path has
a divergence-free polynomial first-order source lift of degrees
`(55,53,55)`, while the fixed-map target translation has first-order degrees
`(31,29,31)`.  All-finite-order liftability is supplied by the separate
formal-orbit theorem.  At `lambda=-4`, the checker reconstructs the exact
two-point affine fiber and uses integer-orbit fiber invariance to rule out a
polynomial lift of the straight target translation.  It does not claim a
high-degree endpoint self-equivalence: the prime ramified discriminant has
ordinary degree thirteen, so divisibility and the exact `mu_5` orbit test
exclude every target candidate through degree twelve.

The next exact frontier requires Singular:

```bash
.venv/bin/python scripts/verify_rank_four_degree_eighteen_target_obstruction.py
```

It computes the logarithmic-derivation nullities
`(0,0,0,0,1,7)` through multiplier degree five.  For target degree at most
eighteen, the endpoint condition leaves four parameters.  Ten exact
constant-Jacobian evaluations in those parameters generate the unit ideal
over `QQ`; Singular returns the reduced basis `[1]`.  This excludes all
endpoint target symmetries through degree eighteen and leaves degree
nineteen as the first unresolved case.  It is an exact Gröbner
inconsistency, not a bounded coefficient search.

Build the Lean interface with:

```bash
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.CollisionFiber
lake build FiniteEtaleKeller.PaperCertificate
```

The relative whole-plane statement and the exact stable-separation
certificates are checked independently in Lean:

```bash
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.WholePlaneStableMultiplicity
lake build FiniteEtaleKeller.StableSeparationCertificates
```

The first build treats all `(B,C)` simultaneously and restricts naturally to
the universal inverse-discriminant open.  The second proves the Fitting
shoelace formula, Laurent translation and unimodular invariance, strict
power-shift separation, and cubic boundary-count separation.  The geometric
normalization/Fitting and boundary-exhaustion inputs remain explicit
interfaces.

The first reconstruction step beyond the common `P=1` plane is checked by:

```bash
.venv/bin/python scripts/verify_two_marked_fiber_gauge_reconstruction.py
```

For a root-marked fiber over `P=p`, the checker extracts the linear
coefficient `g_1/(g_N*p^(N+m))` of the monic root annihilator.  It verifies
that the planes `P=1,2` recover the normalized seed, the gauge exponent, and
therefore the stable Fitting area.  It also checks the sharp periodic
counterexamples on finite collections of torsion planes and recovers `m`
from the pole order `N+m` on a transverse affine line.  The result retains
the inverse-root generator and the base character `P`; it is not an
unmarked-cover Torelli theorem.

The unrestricted finite-sampling question is answered negatively by:

```bash
.venv/bin/python scripts/verify_finite_marked_plane_nonreconstruction.py
```

The checker replaces the common monomial power shift by an arbitrary
polynomial multiplier `R(P)`.  Exact interpolation makes `R=1` on every
prescribed sample plane, while simple roots of `R` each create a Newton
block `(0,0)->(3,0)->(N,1)` and one boundary prime with
`(e,f)=(N-3,1)`.  Squarefree interpolants of increasing degree therefore
give maps agreeing on any finite collection of complete marked inverse
planes but having strictly increasing stable boundary counts
`deg(R)+2`.  The direct quartic calculation independently checks the full
determinant, inverse, and reconstruction identities.

The final full-boundary reconstruction layer is checked by:

```bash
.venv/bin/python scripts/verify_polynomial_gauge_decorated_torelli.py
```

On the clean polynomial-multiplier locus, the boundary ledger selects
`P=0` intrinsically and hence recovers `P` up to scalar.  The checker then
verifies that the unmarked ramified-stratum Fitting divisor excludes root
inversion, kills every punctured-base unit twist, reconstructs the seed and
`R(P)` coefficientwise, and yields exactly the ordinary source--target
scaling action.  The written theorem also records the general converse:
the full finite-normalization morphism carrying its reconstruction boundary
restricts to a left--right equivalence, so that top layer is a complete
stable invariant without a chosen inverse root.

To intentionally refresh the pinned count artifact after changing its
generator, run:

```bash
.venv/bin/python scripts/count_multiplicative_hasse_parameters.py \
  --bound 1000000 \
  --output artifacts/generated-results/multiplicative_hasse_parameters_1000000.json
```

Then record the changed file hash in the canonical note and the Hasse paper
before committing. Continue the broader verifier catalogue with:

```bash
.venv/bin/python scripts/verify_fixed_quintic_arithmetic_zoo.py
.venv/bin/python scripts/verify_stratified_adelic_engineering.py
```

The first command generates the explicit decomposable degree-twelve map
`F_4 o F_3`, checks both determinant-one factors, and records the expanded
coordinate fingerprints and the `4*3` intermediate-field tower.  Pass
`--print-map` to print all three expanded coordinates.  The second command
reduces the pulled-back cubic discriminant modulo the quartic inverse
equation, factors the saturated resultant as `C^8 Q` with `Q` irreducible
of exponent one, separates the other boundary image, and certifies
`Mon(F_4 o F_3)=S_3 wr S_4`.

The cubic and power-shifted multiplicity checkers also exercise the public
`compile_polynomial_to_keller_fiber(..., stable_parameter=k)` path.  They
compare its maps with the symbolic constructions, preserve the selected
inverse polynomial, and audit the returned boundary-count or Newton-area
record.  Stable functoriality is supplied by the corresponding written
theorems rather than inferred from these finite regressions.

The local-to-global checker audits the ramified quintic coefficient CRT:
the prescribed algebras at `2` and `3`, signature `(1,2)`, cycle types `(5)`
at `5` and `(2,2,1)` at `7`, and the determinant-one quadratic-gauge
compilation with complete target `(1,0,-98/809)`.  It first reconstructs the
polynomial through the generic prime-power coefficient synthesizer in
`jcsearch.local_global`, including its common denominator `1261`, and then
feeds it through the shared end-to-end compiler in `jcsearch.keller_fiber`.
It also derives the universal radii `2^5` and `3^3` from the two local
discriminants and reconstructs the fully automatic witness with common
denominator `30241`.

The fixed-map Hasse checker verifies the determinant, target-line
factorization, modulo-`9` Hensel reduction, and the first prime and composite
parameters.  The multiplicative enumerator then lists both the full sufficient
family and the clean prime-support subfamily through `a=10^6`, checks
primitive coordinates and height `32*a`, counts one- and two-prime members,
and records stable SHA-256 digests in the generated JSON certificate.

The locally prescribed common-fiber checker keeps both maps fixed.  It
derives parameter radii `2^9`, `3^3`, and `5`, constructs
`u=95231/69121`, verifies the ramified completions at `2` and `3`, proves
inertness at `5` and signature `(2,2)`, and checks both transported common
targets.

The common-fiber checker synthesizes the arithmetic transfer and stable
boundary results.  It verifies the fixed all-degree pair over `Q`, the
fixed quartic triple over `Q(sqrt(-2))`, the small rational quartic, and the
mod-`17` irreducibility certificates for the connected triple fibers.  The
following search command enumerates the declared rational tangent-chord,
scale, and constant-term boxes, checks the weighted and quadratic
presentation gates, filters at a split prime of `Q(sqrt(-2))`, and recovers
the coefficient-minimal shared polynomial
`9W^4-19W^3+10W^2-8W-4`.

The universal-quartic checker verifies the trace-zero quartic
tangent-chord factorization, its diagonal rank-five trace quadric, the three
possible indefinite real signatures, the normalized weighted parameter
`alpha=u/e-1/2`, and the finite clean-locus exclusions.  Over a number
field, the uniform existence and infinitude proof uses local isotropy,
Hasse--Minkowski, and rational-point density on the resulting smooth quadric,
as recorded in
[`verified/UNIVERSAL_QUARTIC_FIBER_MULTIPLICITY.md`](verified/UNIVERSAL_QUARTIC_FIBER_MULTIPLICITY.md);
it is not a bounded-search conclusion.

The universal quartic gauge checker verifies the second, unconditional
rank-four mechanism.  It raises the quartic lift from `P^4*S^4` to
`P^(m+4)*S^4`, checks the denominator-free determinant and inverse identities,
keeps the selected inverse polynomial fixed at `P=1`, and computes the
ramified-stratum Fitting-support indices `2*m+5`.  The written stable
normalization argument proves infinite multiplicity for every rank-four
finite etale algebra over every characteristic-zero field, including the
anisotropic trace-chord example.

The universal cubic gauge checker verifies the final low-rank step.  It adds
`g_3(P^n-P^3)S^3` to the cubic lift, checks the paired polynomial corrections
and determinant, and confirms that the selected inverse cubic at `P=1` is
unchanged.  The degree-drop polynomial `1+P^(n-1)-P^2` has `n-1` simple
nonzero geometric roots; the written Newton and boundary-exhaustion argument
turns them into `n-1` intrinsic unramified boundary target components.  Their
count separates the maps stably and proves universal infinite multiplicity
in rank three.

The all-degree power-shift checker verifies the uniform replacement
`g_j*P^j*S^j -> g_j*P^(j+m)*S^j` for every `j>=4`.  Representative
three-variable expansions check the determinant, inverse, and reconstruction
identities.  Exact convex-hull calculations verify that the normalized
ramified Fitting Newton polygon has area
`2*N-3+(N-2)*m`; the written stable-normalization argument makes this a
strict stable invariant for every `N>=4`.  This unifies the quartic and
higher-rank multiplicity mechanisms without trace or translation input.

The universal-quintic checker verifies the translated quintic derivative
jets, the primitive relation `(-1,-6,5)` among the three quadratic-gauge
stable-moduli weights, the invariant
`a_5^5/(a_3*a_4^6)=g_5^5*g_1^2/(g_3*g_4^6)`, and its forced pole after
choosing a trace-zero primitive generator with nonzero second trace moment.
The written argument then gives infinitely many stable classes for every
rank-five finite etale algebra over every characteristic-zero field.  The
higher-degree checker verifies the universal top-weight relation and

```text
J_N=a_(N-2)*a_N/a_(N-1)^2
   =(N-1)/(2N)+c_(N-2)/(N^2*s^2),
```

which proves the same conclusion in every rank `N>=6`.  Together with the
power-shifted quartic argument, this gives infinite universal multiplicity
over every characteristic-zero field in every rank at least three, as recorded in
[`verified/UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md`](verified/UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md).
The next checker supplies connected degree-four, degree-five, and degree-six
three-map witness cards, including modular irreducibility, exact targets,
complete inverse identities, and distinct stable invariant values.  The
universal-relative checker verifies the block-triangular determinant through
the compact reciprocal chart, the unchanged-coordinate promotion to one
absolute map of `A^N`, the exact `U_N=V_N x A1` normalization, and the sharp
`N-3`-parameter inverse-polynomial specialization through degree twelve.  The
written theorem proves all-rank finite-etale universality and imports `S_N`
monodromy and primitive-monodromy atomicity; those last two steps are not
inferred from a bounded symbolic computation.  It also distinguishes
presentation dominance from stack descent and essential dimension.  The
adversarial extension gives a closed-form all-rank target for the Osada
`S_N` family `T^N-T-1`, checks it through rank twelve, pins additional
connected, split, and disconnected targets in ranks three through six, and
verifies the genuine degree-drop, bad-translation, and repeated-root
boundaries.

The formalized promoted map, coefficient compiler, and witness cards are:

```bash
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.UniversalPromotedBlock
lake build FiniteEtaleKeller.UniversalParameterCompiler
lake build FiniteEtaleKeller.UniversalParameterQuotient
lake build FiniteEtaleKeller.UniversalParameterWitnesses
lake build FiniteEtaleKeller.UniversalPromotedMap
lake build FiniteEtaleKeller.UniversalPromotedGauge
```

These modules prove the abstract unchanged-coordinate block determinant, the
promoted inverse-polynomial identity, its selected degree, nonvanishing of
compiled top parameters, automatic admissible translation, invariance of the
quotient algebra under the compiler's nonzero normalization and translation,
three exact quartic targets, and the literal promoted map on an `N`-element
coordinate type with its actual full Jacobian block and determinant-one
identity.  They do not formalize the literal promoted full-fiber/compiler
bridge, its geometric degree, `S_N` monodromy, or stable atomicity.  The
low-rank checker verifies the collapse of all three present cubic mechanisms
and the exact biquadratic trace form used in the written two-step Springer
anisotropy proof.

The Hasse-fiber command expands an explicit degree-eight weighted map, checks
its determinant `-38`, proves that its complete target fiber has no rational
point, and audits roots over `R` and every `Q_p` through the elementary
quadratic-residue covering and the two exceptional Hensel lifts.

The normal-covering front end and its first two exact certificates are
replayed by

```bash
.venv/bin/python scripts/verify_normal_covering_certificates.py
python3 scripts/verify_banks_degree_5_10_candidates.py
.venv/bin/python scripts/verify_degree_six_normal_cover_keller.py
```

The first command independently enumerates the groups and all subgroups in
the `S_3` quintic and `C_2^2` sextic actions, proves conjugate coverage and
trivial common core, computes `gamma(S_3)=2` and `gamma(C_2^2)=3`, and checks
the exact ramified-prime Hensel witnesses.  The second validates the pinned
necessary-candidate transcription of Banks' Table C.1; it does not assert
that every candidate row is arithmetically realized.  The third compiles
`(T^2-2)(T^2-17)(T^2-34)` with the shared quadratic-gauge compiler, verifies
the target `(1,0,528/577)`, and expands its determinant-one Keller map.
For larger finite groups, the GAP front end is loaded with
`Read("scripts/normal_covering_certificate.g");`; the checked-in small
certificates deliberately have a dependency-free Python replay.

The fixed-quintic commands check one determinant-`-2` map and its finite
certificate ledger: all three real quintic signatures, all five transitive
groups `C_5`, `D_5`, `F_{20}`, `A_5`, and `S_5`, split and
quadratic-times-cubic fibers, and all seven unramified partitions modulo
`7`.  The group certificates use witness-prime factor patterns, an explicit
order-five automorphism, a pair-sum resolvent, and Cayley's sextic resolvent.
The clean Hasse row has normalized polynomial
`(T^2-8T+47)(T^3+8T^2+12T+8)` and common quadratic resolvent
`Q(sqrt(-31))`; only `2` and `31` need special local witnesses.  The
original `Q(sqrt(-3))` row remains an independent regression and supplies
the `Q_5` trace obstruction to the standard pure-cubic infinitude route.
Infinitude inside this particular split-seed pencil remains open.

The local `Q_2` action-certificate branch has a separate three-layer replay:

```bash
python3 scripts/verify_gq2_permutation_action.py \
  arithmetic/certificates/gq2_s3_x3_minus_2.json --json
python3 scripts/verify_gq2_permutation_action.py \
  arithmetic/certificates/gq2_s4_mixed_action.json --json
python3 scripts/verify_gq2_permutation_action.py \
  arithmetic/certificates/gq2_common_quintic_stable_pair.json --json
.venv/bin/python scripts/verify_gq2_action_first_keller.py
.venv/bin/python scripts/verify_gq2_s4_quartic_keller.py
.venv/bin/python scripts/verify_marked_q2_stable_separation.py
gp -q scripts/verify_gq2_s4_local_models.gp
gp -q scripts/verify_gq2_local_decompositions.gp
```

The first command is dependency-free and evaluates the exact Roe--Turturean
word ledger, including the finite `omega_2` powers and the normal 2-core
condition.  Its named comparison proves that the marked `S_3` action is the
splitting action of the tame Eisenstein cubic `T^3-2` over `Q_2`.  The second
command translates that polynomial into the degree-three quadratic-gauge
formula, verifies determinant one after output scaling, and checks the complete
inverse polynomial `(S+1)^3-2`.  The mixed `S_4` checker enumerates the three
marked `x_0` orbits over the fixed tame frame, evaluates the candidate
quadratic obstruction as `(1,0,0)`, and compiles
`T^4+4T^2-4T+2` into a determinant-one complete quartic fiber.  The first
PARI/GP command proves that the three classified local quartics have closure
group `S_4`, inertia `A_4`, wild inertia `V_4`, the displayed exact
ramification groups, and normalized relative Stiefel--Whitney bits `(0,0,1)`.
It also verifies the three resolvent Kummer square classes and their sole
product relation over the tame `S_3` closure.  The unique nonzero obstruction
therefore matches the worked `x_0=(12)(34)` orbit to
`T^4+4T^2-4T+2`; the other two orbits remain unordered.
The common-quintic checker verifies the exact unramified marking
`sigma=(1234)(5), tau=x_0=x_1=1`, global irreducibility modulo `17`, both
determinant-one inverse equations, and the stable unit-rank separation
`1 != 2`.  The final PARI/GP
command recomputes exact
ramification-index/residue-degree decompositions at `2` for the selected
quartics, all ten fixed-quintic zoo rows, and the separate quintic witness
card.  It uses maximal-order prime-ideal decomposition rather than bounded
`2`-adic precision.  The checked-in table was last recomputed with PARI/GP
2.17.4 on arm64 Darwin (GMP 6.3.0).  The combined target is:

```bash
make verify-gq2-local-fibers
```

The height-`21` five-row witness card and its separate bounded discovery
audit are reproduced by

```bash
.venv/bin/python scripts/verify_universal_quintic_calculator.py
.venv/bin/python scripts/search_universal_quintic_calculator.py --bound 21
```

The first command uses only exact rational arithmetic and finite-field
factorization.  The second requires PARI/GP, enumerates primitive
projective targets through height `21` modulo the sign involution, and uses
`polgalois` only after exact discriminant and Frobenius-pattern prefilters.
It is bounded computational minimality evidence, separate from the
oracle-free certificates for the five displayed rows.

The mechanically generated finite ledger is checked by

```bash
.venv/bin/python scripts/verify_fixed_quintic_certificate_ledger.py
```

It recomputes all ten rows, their real-root counts and witness-prime
patterns, the seven modulo-`7` partitions, and the `-48*Pi^8` coefficient
Jacobian.  It also runs the three canonical exact checkers and compares both
the Markdown table and
`artifacts/generated-results/fixed_quintic_certificate_ledger.json` with
the generated data.  Pass `--write` only when intentionally refreshing both
generated forms.

Its bounded height search requires PARI/GP:

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_targets.py
.venv/bin/python scripts/search_fixed_quintic_hasse_curves.py
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_seven.py
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_thirty_one.py
.venv/bin/python scripts/analyze_fixed_quintic_hasse_minus_thirty_one.py
.venv/bin/python scripts/search_fixed_quintic_hasse_rational_curves.py
.venv/bin/python scripts/search_fixed_quintic_hasse_elliptic_slice.py
.venv/bin/python scripts/search_fixed_quintic_hasse_rank_one_slice.py
```

The first command's default box is stated in its output.  It reports the two
sign-related Hasse targets of projective height `257280` and no other target
below the previous height `458080` in that box.  This is search evidence,
not a global height-minimality claim.  The second command verifies an exact
rational parametrization of the common-quadratic-resolvent incidence, checks
rank-two and rank-one elliptic slices, and searches a bounded proportional
family for irreducible candidates having cubic roots over `Q_2`, `Q_3`, and
`Q_5`.  Its only small-prime survivors are four presentations of the known
Hasse target; it does not test every completion and is not an infinitude
proof.
The rank-two elliptic-slice command closes the
`kappa/A=-1, R=1` route exactly:
the Mordell equation forces `v_2(Pi)=-(2m+1)` and
`v_2(A)=-(3m+1)`, while the cubic factor has a single Newton-polygon
slope `3m+7/3`.  Neither it nor the discriminant-`-3` quadratic has a
`Q_2` root.  Its default 624-point Mordell--Weil enumeration is only a
regression for that all-points proof.
The final command closes the rank-one `kappa/A=5/4, R=4` route at `Q_5`:
the elliptic equation forces `v_5(Pi)=-2m` and `v_5(A)=-3m`; after
translating the cubic by `T=2A+Y`, its unique Newton slope is
`3m-1/3`.  The quadratic has discriminant `-48`, so neither factor has a
`Q_5` root.  Its 24-point multiple regression is again secondary to the
uniform valuation proof.
The third command varies squarefree shared quadratic resolvents.  In its
default integral box it finds a new `Q(sqrt(-7))` Hasse target
`(-7,387/14,400/2401)` of projective height `132741`.  The fourth command
independently checks its factorization, irreducibility, common resolvent, and
exact local witnesses at `2`, `5`, `7`, and `79`.
The wider command

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py \
  --d-bound 60 --r-bound 20 --a-bound 50 --pi-bound 40
```

also finds a `Q(sqrt(-31))` target
`(5,-144/5,-188/3125)` of projective height `90000`.  The fifth command
above audits it independently; only the exceptional primes `2` and `31`
require local witnesses.
The sixth command verifies two exact continuation reductions for the
`Q(sqrt(-31))` row: a genus-two cube curve on the fixed normalized-factor
slice, and a rational trace quadric with a quartic cube condition for affine
variation of the two field generators.  It also enumerates rational
coordinates of height at most `600` on the genus-two slice; only the known
coordinate `Pi=5` occurs.  This last statement is bounded search evidence.
The seventh command proves an exact obstruction to every base line through
`(A,R,Pi)=(-8,2,5)` on the fixed-`-31` common-resolvent double cover: after
recursive square-root reconstruction, the residual ideals have Groebner
basis `[1]` on all three projective direction charts.  It also excludes
every degree-at-most-two curve on each coordinate-fixed slice
`A=-8`, `R=2`, and `Pi=5`: all twelve weighted-projective charts have empty
fiber modulo the good prime `32003`, hence empty characteristic-zero
generic fiber by properness.  Finally, it tests 15024 genuine general
quadratic parametrizations with six integral coefficients in `[-2,2]` and
finds no square pullback.  The line and coordinate-slice results are exact;
the general quadratic result is only bounded evidence.  New primes ramifying
in cubic specializations remain an additional all-prime local obstacle.

The larger fixed-discriminant integral search is

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py \
  --d-value -31 --r-bound 100 --a-bound 200 --pi-bound 200
```

It finds only the certified point and its sign mate.  The low-denominator
rational search is

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py \
  --d-value -31 \
  --r-bound 20 --r-denominator 4 \
  --a-bound 30 --a-denominator 4 \
  --pi-bound 30 --pi-denominator 4 \
  --show-failures
```

It finds two further common-resolvent presentations, both failing at `17`.
Both commands are bounded search evidence, not finiteness theorems.

The last command audits the constructive CRT/weak-approximation lift and an
explicit nonsurjective type-`(3,2)` quintic seed with trivial Hessian
stabilizer and complete fibers of all three quintic signatures, each with
cycle types `(5)` at `7` and `(2,2,1)` at `11`. The preceding adelic command
audits an explicit totally imaginary quartic complete fiber that is inert at
`7` and has unramified splitting type `(2,1,1)` at `11`.

The linear-torus-free quadratic-gauge specialization has a separate exact
certificate:

```bash
make verify-linear-torus-free
```

It checks determinant one, a four-point rational collision, all 734
coefficient equations in `B F = JF A x`, and a displayed `18 x 18` primitive
integer minor of determinant `-5`.  Thus every infinitesimal linear
source-target symmetry vanishes; conjugation makes the result invariant
under independent linear coordinate changes.  A dependency-free clean-room
replay rebuilds the sparse rational polynomial calculation and verifies the
matrix by Bareiss elimination.  A separate parameter calculation proves
within the SymPy checker that the same rational minor has determinant
`(10935/4) g_4^6/g_1^6`, so every admissible quartic quadratic-gauge map has
the same linear-symmetry exclusion.  For the displayed small-coefficient
map, both implementations also verify that the complete `785 x 24` system
allowing constant terms in both vector fields has full column rank.  Hence
the example remains free of affine-linear torus equivariance after
independent affine coordinate changes.

The intrinsic algebraic-torus strengthening is checked by:

```bash
make verify-algebraic-torus-free
```

On the canonical ramified normalization stratum it reconstructs
\(J(P,r)=-1-3Pr^2+4P^4r^3\), computes the scheme-theoretic stabilizer
\(\beta=\alpha^{-2}\), \(\alpha^5=1\), and verifies that its tangent matrix
has determinant \(5\).  The Newton-support pass checks all six permutations:
only the identity and the involution
\(\left(\begin{smallmatrix}-2&-1\\3&2\end{smallmatrix}\right)\) are integral
unimodular, and the ordered second-boundary image rejects the involution
because it does not preserve the intrinsic base character \(P\).  The
checker also verifies the resulting explicit
\(\mu _5\) source--target symmetries of the displayed map.  The canonical
boundary argument then upgrades the conclusion from affine-linear actions:
a connected torus acts trivially on the decorated stratum, hence fixes the
prime nonnormal discriminant hypersurface pointwise; the weight-space lemma
in the canonical note forces the target action to be trivial, and
\(S_4\) deck rigidity forces the source action to be trivial.  Thus no
polynomial left--right representative is algebraic-torus-equivariant.  This
does not classify all discrete or unipotent polynomial self-equivalences,
and literal symmetry-freeness after identity stabilization is impossible:
the added identity coordinates carry tautological torus actions.  The exact
stable conclusion is that every connected action on the pulled-back
decoration is vertical over its intrinsic two-torus; no splitting of such a
vertical action is claimed.

The same checker expands the rational-root sparse representative from
`G(S)=S(S-1)(S-2)(3S+2)`, verifies determinant one and its displayed
four-point rational fiber, and obtains component support counts `(7,51,38)`
and ordinary degrees `(7,26,24)`.  A symbolic coefficient audit shows that
exactly seven generic support coefficients contain `g_2/g_1`, while every
other coefficient is a nonzero Laurent monomial in the admissible
`g_3/g_1,g_4/g_1`.  Hence `g_2=0` is support-minimal in this fixed normal
form.  No absolute sparsity claim under arbitrary polynomial left--right
changes is made.

The bounded exact polynomial left--right sparsity search is:

```bash
.venv/bin/python scripts/search_quartic_lr_sparsity.py \
  --source-degree 2 \
  --target-degree 2 \
  --tadic-max-exponent 12 \
  --scaling-bound 16 \
  --output artifacts/generated-results/quartic_lr_sparsity_search.json
```

It tests every rational exceptional parameter in 15 one-monomial source
shears and 15 one-monomial target shears, all 76 two-term representatives of
the essential source jet through exponent 12, and 25281 rational diagonal
scalings.  No searched nonidentity shear improves support `(7,51,38)`.
This is bounded computational evidence, not an absolute minimum.  The
scaling search finds the balanced exact height improvement
`(alpha,beta)=(1/4,12/5)`.  The generated record has SHA-256
`3fea20be042106fb5fe452ebe241dc5c3316eed6a893b00bf7ca2bcc0bef1b70`.

The continued two-move circuit search is:

```bash
.venv/bin/python scripts/search_quartic_lr_two_move_circuits.py \
  --source-degree 2 \
  --source-parameter-bound 4 \
  --jet-max-exponent 12 \
  --workers 4 \
  --output artifacts/generated-results/quartic_lr_two_move_circuits.json
```

It checks 330 rational source-shear/optimal-linear-target circuits and all
286 three-term essential jets through exponent 12, exactly at every rational
exceptional jet parameter.  The best nonidentity jet is
`z -> z+(16/7)y^2`, with support `(7,51,39)`.  Target monomial cleanup
through degree three does not improve it; among all exceptional elementary
second source shears through degree two, only the literal inverse returns to
96.  Adding a second structured monomial through exponent 12 also fails to
improve the 97-term near miss.  This remains bounded evidence.  The generated
record has SHA-256
`28c9e1c2ed9c765fef7c51d7e8ace3262c6fac337c2b11d723f3dccdb3781826`.

The remaining constant-`C` boundary has a separate exact Singular
certificate:

```bash
make verify-weighted-boundary
```

It computes exactly two primary components, checks their declared radicals,
and verifies that the reduced affine-three-space components meet in an
affine plane.

The reduced global attachment of the open torus orbit is checked by

```bash
Singular -q scripts/verify_foundational_reduced_gluing.sing
```

This verifies the degree-ten toric closure and its two boundary lines.

The heavier regression target also checks the explicit degree-five family and
its rank-two symplectic descent:

```bash
.venv/bin/python scripts/verify_degree_five_rank_two_descent.py
```

This exact calculation constructs the relative Hamiltonian over
`Q(lambda)`, extracts all four negative-`X` residue coefficients, proves the
unique parameter-dependent shear cancels them, and verifies the normalized
base brackets and polynomial source automorphism.  It normally takes roughly
half a minute in the pinned symbolic environment.

The smaller classical degree-drop viability test is:

```bash
.venv/bin/python scripts/verify_quartic_weighted_map.py
```

Besides the quartic inverse and collision, it identifies the seed with
\((\kappa,\tau)=(-5,0)\), specializes the rank-two completion, checks all six
Poisson brackets and the canonical coordinate change, and transports the
generic degree-four cover and an explicit two-point collision.  It reports
fiber orders \((4,3)\).  This is an exact classical certificate, not an
\(A_2\) quantization.

The rebuilt restricted quantization test for those exact \((4,3)\) symbols
is:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_quartic_degree_drop_quantization.py \
  --certificate artifacts/generated-results/quartic_degree_drop_quantization.json
```

It derives the specialized Bernstein bounds, solves the complete
parity-preserving \(\hbar^3\) affine equation, and writes an exact six-term
dual cocycle proving the rank jump \(143\to144\) at \(\hbar^5\).  It also
rebuilds the unrestricted first-order kernel, removes the complete
target-Hamiltonian gauge, projects the next Maurer--Cartan quadrics, and
tests all five surviving coordinate axes through the coupled
\(\hbar^2/\hbar^3\) equations.  Its bounded low-support pass classifies one
coordinate \(\mathbf P^4\), nine isolated rational directions, and no
algebraic support-two directions.  All nine isolated directions fail, while
a uniform third-order relaxation reduces the \(\mathbf P^4\) to an explicit
residual \(\mathbf P^2\).  Its genuine compatibility obstruction factors as
\((21a+28b+64c)^3/21^3\), so one rational \(\mathbf P^1\) reaches
\(\hbar^3\).  Parameterizing that line by
\((4,-3,0)+t(0,16,-7)\), the complete 38-dimensional lower-lift calculation
over \(\mathbb Q(t)\) gives the exact fourth-order rank jump
\(143\to144\).  Its six-term cocycle has sole denominator factor \(t\);
exact audits at \(t=0\) and projective infinity give the same jump.  Thus
the complete projective resonance line is eliminated at \(\hbar^4\).
The discrete base ranks are computed over \(\mathbb Q\) and independently
repeated over \(\mathbf F_{32003}\).  This certifies obstructions only for
the displayed normal ordering and inherited filtration; it is not an
\(A_2\) nonexistence theorem.  The recorded JSON certificate has SHA-256
`04646808c526697e7538a268605da6df1b5e3a66c51a6a5e3c1d68c80ab20ab9`.

The degree-five filtered contact problem has a separate two-invariant audit:

```bash
.venv/bin/python scripts/verify_degree_five_torus_module.py
```

It verifies the torus-gauge root recurrence over `Q[u,gamma]`, proves the
all-order profile `24m+1`, and checks survival of the candidate class in the
invariant-ring-saturated equivariant target quotient.

The minimal opposite-weight quadratic Rees witness is reproduced by

```bash
.venv/bin/python scripts/search_rees_torsion_witnesses.py --max-target-degree 0
```

It finds the unique constant opposite-weight pair
`(partial_B,partial_C)` of weights `(1,-1)`, computes its exact second
fundamental form by two independent formulas, and returns the nonzero leading
normal symbol `-146880u^5/7` in the third saturated summand `R/(gamma)`.

The finite LR Rees/SAGBI module calculation is reproduced by

```bash
.venv/bin/python scripts/compute_lr_rees_sagbi_modules.py
# or, including the dependency-free separator replay:
make verify-lr-rees-sagbi
```

It constructs the three-generator target-invariant SAGBI basis, the target
modules and initial lifts in weights `p=+-1,+-2`, and the saturated normal
quotient.  It certifies a linear weight-one degree drop `39 -> 34`, computes
the further subduction to a new degree-`29` initial-module generator, computes
the complete `3 x 24` quadratic matrices for `p=1,2`, proves the structural
cutoff `|p|>=3`, and performs exact Singular module membership.  The sole
new `p=2` column modulo the full `p=1` image is
`II_(F,2,-2)(partial_A,A^2 partial_A)`, with remainder
`-987/395*e_C`.  The generated JSON certificate is
`artifacts/generated-results/lr_rees_sagbi_module_computation.json`.

The decisive independence statement has a dependency-free replay:

```bash
python3 scripts/audit_lr_rees_sagbi_module_certificate.py
```

The all-order constant-direction rooted-tree normal classes are reproduced by

```bash
.venv/bin/python scripts/compile_lr_rooted_tree_classes.py --max-order 12
python3 scripts/audit_lr_rooted_tree_normal_classes.py
```

The compiler works in exact torus semi-invariant coordinates, reproduces the
known `II_F(partial_B,partial_C)` residue at order two, and constructs the
weight-zero ladders `tau_2=B(C)`, `tau_3=A(C(C))`,
`tau_(n+2)=B(C(tau_n))`.  A fixed `3 x 3` transfer matrix at
`(u,gamma)=(1/6,0)` and a positive-coefficient Cayley--Hamilton recurrence
prove that the third saturated normal residue, hence its associated-graded
symbol, is nonzero for every order `n>=2`.  This is an all-order theorem for
the individual tree classes, not a proof that the same class survives the
sum and lower-jet variation in a mixed BCH/LR forcing coefficient.

The balanced linear-in-`X` mixed BCH sector is reproduced by

```bash
.venv/bin/python scripts/compile_lr_mixed_bch_classes.py --max-k 3
python3 scripts/audit_lr_mixed_bch_classes.py
```

Here `X=N*(x,0,-3z)`, `D_B=ell_F(partial_B)`, and
`D_C=ell_F(partial_C)`.  The checker proves `[D_B,D_C]=0`, collapses the
balanced order-`2k+1` BCH sum to
`binomial(2k,k)*(ad(D_B)*ad(D_C))^k*X`, and derives the exact third-normal
recurrence

```text
c_(k+1) = -73440*(k+3)*(2k+7)*c_k,
c_1 = 14438891520/2401.
```

Thus the actual multihomogeneous BCH coefficient is nonzero in every odd
order.  It survives the saturated linear target quotient, but with target
amplitudes `s,t` it is multiplied by `s^k*t^k`; consequently this sector
alone is not universal over the full lower-jet scheme.

At `(u,gamma)=(1/6,0)`, the covector `(0,-144/79,1)` descends through the
saturated normal relations, kills all 24 `p=1` columns, and takes the value
`-987/395` on the new `p=2` column.  This matches the Singular normal
remainder without using SymPy or Singular in the replay.  The main checker
also computes the exact annihilator `(gamma,6*u-1)`, proving that the
`p=2` image modulo `p=1` is one reduced residue-field copy of `Q`.

The full normalized degree-five seed surface is checked by

```bash
.venv/bin/python scripts/verify_degree_five_flux_surface.py
```

This exact two-parameter calculation works over `Q(a,tau,s_2)`, verifies the
uniform adapted coordinate and quotient brackets, extracts the complete four
term Laurent obstruction, and proves that its unique quadratic shear makes
the Hamiltonian polynomial.  It takes several minutes in the pinned symbolic
environment.

The exceptional `kappa=-1` chart and its pole-filtered monomial shear
responses through cubic degree are replayed by

```bash
for degree in 0 1 2 3; do
  .venv/bin/python scripts/explore_kappa_minus_one_flux.py --shear-degree "$degree"
done
```

Each run verifies the replacement determinant and quotient brackets, all
three Hamiltonian components, and every negative-`X` residue coefficient.
The full exceptional-divisor completion is checked by

```bash
.venv/bin/python scripts/explore_kappa_minus_one_flux.py \
  --x-degree 1 --shear-degree 1
```

It proves that the complete principal part is canceled by
`2(2*tau^2-15*tau-18)*X*Q/105`.

The full degree-six generic chart, exceptional divisor, and fixed-`gamma`
specialization are checked by

```bash
.venv/bin/python scripts/verify_degree_six_flux_surface.py
.venv/bin/python scripts/verify_degree_six_kappa_minus_one_descent.py
.venv/bin/python scripts/verify_degree_six_fixed_gamma_descent.py
```

These verify the three-parameter generic seed chart and the full exceptional
divisor, componentwise Hamiltonian identities, complete residues, and unique
completing shears.  The generic symbolic replay is a heavy calculation.

The all-degree Laurent recurrence and exact fixed-`kappa=-9` probes in degrees
seven and eight are checked by

```bash
.venv/bin/python scripts/verify_four_residue_recurrence.py
.venv/bin/python scripts/explore_all_degree_fixed_gamma.py 7
.venv/bin/python scripts/explore_all_degree_fixed_gamma.py 8
```

The direct second-Weyl-algebra parity test is replayed by

```bash
.venv/bin/python scripts/explore_degree_five_a2_subprincipal.py
```

It solves the `hbar^3` equation exactly, retains its full 42-dimensional
solution space, and proves that the `hbar^5` cokernel contains `1=0`.  This is
an obstruction only to the declared parity-preserving filtered ansatz.

The parameter-uniform third-order lift, four bounded fifth-order periods,
their common cubic locus, and the genuine nonlinear fifth-order equations on
that locus are replayed by the commands in
[`extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md`](extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md#10-reproduction).
The small relative-family package, including the two exact coprime Fitting
charts, the interior Kuranishi shadow, and the root-at-infinity valuation
filtration, is checked independently by

```bash
make verify-degree-five-relative-quantization-family
```

It verifies the valuation weights
`(X,Q,W,R,gamma)=(1,-1,-2,1,0)` and the induced pure correction weights
`(S_2,T_2)=(4,5)` and `(S_4,T_4)=(10,11)`.  It does not promote the modular
length-218 Fitting computation to characteristic zero.
The decisive exact cubic-field check is

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --print-radical-basis --seventh-line
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --seventh-component-elimination
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --seventh-component-elimination \
  --seventh-component-program-output \
    artifacts/generated-results/degree_five_cubic_h7_unit_certificate.sing
.venv/bin/python \
  scripts/verify_degree_five_cubic_h7_unit_certificate.py
```

It verifies all 680 projected quadratic equations at an explicit
two-coordinate lower lift and solves the unreduced fifth-order correction
equation with a 14-term particular solution.  The exact radical is one
affine 27-space with a six-linear-form nonlinear core, and the entire
explicit one-parameter line is obstructed at order seven.  The second
command proves over `GF(32003)` that the full 20-column order-seven matrix
has constant rank six and that its 401-polynomial consistency ideal, in
only ten effective parameters, is the unit ideal.  Repeat it with
`--prime 31991 --a 109 --tau 28672` for the second good-prime certificate.
The third command performs the characteristic-zero lift over the cubic
field.  Batched constant-field elimination replaces the former
27-variable function-field solve.  The resulting 401 equations contain 27
nonzero constants; the selected \(X^{18}\) residual has an explicit Bézout
inverse, and Singular verifies both a direct one-generator identity and a
one-term degree-zero lift of \(1\).  The final command replays the pinned
1.3 MB Singular certificate in under a second.  Hence the complete reduced
affine 27-space of fifth-order lifts is obstructed at order seven.

The bounded audit of the standard-support parameter Fitting scheme is

```bash
.venv/bin/python scripts/compute_degree_five_qper_fitting.py \
  --timeout 120
```

It finishes at the three fixed good primes, checks the common
21-generator leading-monomial staircase, saturation exponent 12, dimension
zero, and length 218.  This is a stable modular certificate, not a
characteristic-zero proof.  The opt-in rational reconstruction experiment is
checkpointed by

```bash
.venv/bin/python scripts/compute_degree_five_qper_fitting.py \
  --prime 0 --method modular-rebuild --timeout 900 \
  --basis-output \
    artifacts/generated-results/degree_five_qper_fitting_basis_Q.sing

.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check shape
.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check groebner --jobs 8 --timeout 600
.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check boundary-unit --timeout 600
```

The rebuilt 20,840,615-byte rational candidate has SHA-256
`25788668021f563e17373b55703a08ef5693576077ebdbe53c4c3f2c659d98e6`.
The 20 adjacent staircase \(S\)-pairs give an exact Gröbner certificate, and
the exact boundary-unit check proves that the candidate itself is saturated.
These checks do not prove that it equals the saturated maximal-minor ideal.
That last identification requires fraction-free quotient identities in both
containment directions.

The bounded modular audit of the 16 input-containment quotients is built
incrementally by:

```bash
.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --no-resume --skip-diagnostic

.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --primes 70001 70003 70009 70019 70039 70051 \
    70061 70067 70079 70099 70111 70117 \
  --jobs 8 --checkpoint-every 64 --skip-diagnostic --timeout 600

.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --prime-start 1000000000 --prime-count 100 \
  --jobs 8 --checkpoint-every 64 --skip-diagnostic --timeout 600

.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --prime-start 1000002043 --prime-count 500 \
  --jobs 8 --checkpoint-every 64 --timeout 600
```

It resumes
`artifacts/generated-results/degree_five_qper_input_quotients_modular.json`,
uses deterministic prime order even with parallel Singular workers, and
keeps one good image out of the CRT pool.  The recorded run has 613
support-stable good primes, two support-unlucky primes, 18,116 CRT bits, and
30 of 11,701 balanced reconstructions confirmed at the held-out prime.
These data are modular evidence and a coefficient-height diagnostic, not an
exact containment proof.  Use `--skip-diagnostic` when only extending the
checkpoint; the held-out reconstruction pass uses GMP-backed FLINT
arithmetic.  The recorded compact checkpoint has SHA-256
`70a690fd53b4b3a15d4eebf5116acf57b7d0079a8f96a1aadfb2826da86d0481`.

The low-support unrestricted odd audit is replayed by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_mixed_quantization.py
```

It reconstructs the 38-dimensional gauge quotient and its 41 quadratic
obstruction equations and classifies every exact-support-two branch,
including the nine quadratic closed points.  The only three mixed
support-two directions reaching the simultaneous third-order equation retain
63 lower-lift parameters; after adding all 2079 enlarged obstruction
coefficients and every bounded next correction, the constant raises the span
rank from 626 to 627 over each of the good primes 31991, 32003, and 65521.

The generic residual-line and exact support-three audits are replayed by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_mixed_function_field.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_support_three.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_support_three_points.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_support_three_curves.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_support_three_curves.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_residual_five_space.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_residual_support_three.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_residual_fourth_identity.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_cones.py L1 --exact
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_cones.py L2 --exact
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_branch_lifts.py L1 --order 8
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_branch_lifts.py L2 --order 8
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_l1_high_support.py
```

The first constructs a fixed three-monomial residue over `QQ(r)`, finds its
sole lower-basis pole at `r=-3/4`, and eliminates that exceptional
specialization exactly; hence the complete residual projective line is
closed.  The next two commands classify all 8436 support-three coordinate
charts and eliminate all 66 isolated closed-point classes.  The curve
commands compress the 149 positive-dimensional line/conic charts to 23
closed points, force all 19 quadratic classes to zero scale, and eliminate
the four rational survivors by exact `646->647` rank jumps.  Thus every
exact-support-three branch is closed.

The last three commands treat the residual projective four-space inside `L_2`.
The first proves exactly that its nonzero-scale locus is the union of the hyperplane
`z2+2*z3-9838*z4/105=0` and one explicit primitive quadric.  Uniform
fourth-order obstruction on those two threefolds is supplied by the last
command: the quadric is a binary form in the two coupling coordinates and
splits into conjugate hyperplanes over `QQ(sqrt(-2))`; a fixed three-term
residue handles every nonzero-coupling chart, and an exact 16-term residue
handles their rank-zero intersection plane.  The middle command
intersects that locus with every residual exact-support-three chart: the
already-closed residual line and twelve closed points result, and all twelve
points have exact `626->627` next-order rank.  Thus exact support three is
closed even inside the coordinate planes contained in `L_1 union L_2`, and
the final command eliminates every nonzero-scale branch in the residual
projective four-space, in every support.  The final five commands attack
higher support: they compute the exact `L1` and `L2` normal cones; show that
a generic `L1` normal branch is an exact 26-support classical solution while
the analogous `L2` branch is obstructed at its next Kuranishi equation; and
isolate a projective high-support `P6` that reaches nonzero filtered scale.
Uniform fourth-order continuation on that `P6` is the next open calculation.
These are statements about the displayed classical symbol only; they do not
prove `(DC_2)`.

The separately authored foundational Lean certificate is optional because it
downloads an external pinned checkout:

```bash
make verify-lean-foundational
```

GitHub Actions runs this target in the required `formal-lean` job using the
pinned upstream commit and Lean action.  The `formal-local-lean` job
audits both publication certificates and runs `make verify-lean-local`, so
every checked-in Lean package is built.  The `papers` job compiles the
finalized and active manuscripts listed in `papers/README.md`, while parked
manuscripts remain available for direct local builds.  The
`macaulay2-independent-check` runs the pinned Macaulay2 comparison.  Together
with the four Python matrix jobs, these are the complete CI verification
pipeline.  The final `verification-complete` job is the single aggregation
check intended for GitHub branch protection.

## Cancellation programme

```bash
make verify-master
```

The direct three-puncture reciprocal ledger and its all-degree
polynomiality obstruction are replayed separately by:

```bash
.venv/bin/python scripts/verify_puncture_rank_frontier.py
```

The checker derives the universal determinant ledger, enumerates all
two-character coefficient matrices in `[-2,2]` into 129 rank failures,
392 nonsaturated class-lattice failures, and 104 saturated bases, and lists
the 44 primitive positive `(r,a,b)` ledgers with coefficients at most four.
It then verifies the boundary-moment eigenvalue recurrence used by the
written all-degree proof and retains the eleven degree-four through
degree-seven coefficient factorizations as regressions.  No candidate
reaches polynomiality, so the calculation does not claim a Keller map or a
complete collision.  Singular is not required.

The surviving two-reconstruction-variable `A^6` core and its nonlinear
screens are checked by:

```bash
.venv/bin/python scripts/verify_three_puncture_nonlinear_frontier.py
```

In addition to the two dimension-free rank-drop gates and the 80 fixed
degree-at-most-three skeletons, this computes the zero-modification slice of
the proposed coupled ansatz.  For arbitrary affine `P,Q` over `Q(c,v)`,
arbitrary polynomial `H,S`, and two arbitrary affine transverse outputs, it
forms eleven determinant coefficient equations and the Plücker quadric.
Their exact Gröbner basis is `(1)`.  This proves that the next search must
give at least one transverse output degree at least two; it does not exclude
the remaining degree-two through degree-four systems.  The checker then
keeps `P,Q` arbitrary affine, makes the fourth output a completely general
degree-at-most-two polynomial, and proves unit coefficient ideals for the
eight transverse skeletons
`u,z,w,u+z,u+w,z+w,D0+z,D1+w`.  The two-general-quadratic-output system
remains open.  Uniformly for every nonconstant affine direction `C`, four
projective pivot trees avoid a slow monolithic Gröbner calculation.  The
nonzero-`r` chart `C=r+g*u+a*z+b*w` follows the exceptional divisors
`p1,q1,p2,q2,g,p0-a`; the `r`-free charts
`C=u+a*z+b*w,z+b*w,w` have chains `(p1,q1)`, `(q1-b*p1)`, and `(p1)`.
The coefficient and augmented ranks differ on every open and terminal
branch.  Thus no affine `C` can be paired with a general
degree-at-most-two fourth output.  Arbitrary quadratic `C`, two
simultaneously general quadratic outputs, and degree-at-least-three fourth
outputs remain open.

The checker also closes the exposed-`r` simultaneous-quadratic boundary
`P=Q=0`: for two general degree-at-most-two outputs, an eight-step
coefficient pivot tree ends with coefficient/augmented ranks `6/7`.
Degree three is the first zero-slice escape.  With `q=1-c*v` and `C=w`, it
verifies
`D3=(u*(-q^2-v*r*q+2*v^2*r^2)-6*v^2*z)/q^3` has slice determinant one.
The polynomial numerator instead has full determinant
`q^3+6*r*v^2*w`; its cofactor derivation has an explicit common zero.
Thus this exact rational survivor neither satisfies polynomiality nor lifts
to a full Keller map, and a next quadratic/cubic search must use nonzero
`P,Q` or `H`.

This target includes the exact quadratic-gauge/cancellation intersection
regression.  To run its symbolic `N=4,5,6,7` discriminant and all-factorization
checks directly:

```bash
.venv/bin/python scripts/verify_quadratic_cancellation_intersection.py
```

The all-rank clean quadratic-gauge stable-moduli and marked-stabilizer
certificate is:

```bash
.venv/bin/python scripts/verify_quadratic_gauge_stable_moduli.py
```

Besides the two-torus quotient and its saturated compiler-slice invariants,
it verifies the weight-one global receiver slice `lambda=u_5/u_4`, the
exact finite-etale descent identity, the universal discriminant inequalities
through rank 128, and exact discriminant supports in ranks four through
eight.  These checks certify the all-rank written proof that
`D_N=(2,N,1)`, corresponding to `P^2*B^N*C`, is uniquely exposed by
`(1,N+1,N)`.  Kuroda's and
Derksen--Hadas--Makar-Limanov's theorems are external mathematical inputs,
not re-proved by the script.

The minimal-boundary gateway and classification program has a separate fast
cubic certificate:

```bash
make verify-minimal-boundary
```

The eight-predicate invariant pipeline on finite canonical-normalization
exports has its own dependency-light target:

```bash
make verify-minimal-boundary-pipeline
```

It checks weighted degrees `3,...,8`, six cancellation parameter pairs,
quadratic-gauge degrees `3,...,8`, five single-defect perturbation/spectator
records, and a relabeling/order blindness regression.  Regenerate
`artifacts/generated-results/minimal_boundary_pipeline.json` with

```bash
.venv/bin/python scripts/verify_minimal_boundary_pipeline.py \
  --write-artifact
```

This target starts from exact finite exports; it does not compute a canonical
normalization or its intrinsic marking from a bare polynomial map.

It proves that the weighted geometric-degree-three seed has no modulus, that
the cancellation degree equation forces `(m,r)=(1,1)` and `h=3+9A`, and that
both maps are carried to the foundational polynomial by explicit diagonal
source and target automorphisms.  It also verifies the cubic two-place toric
defect atlas and the diagonal reciprocal-lift obstruction.  The accompanying
proof uses Abhyankar--Moh to make the one-place plane-core marking automatic.
It also checks the positive quotient tower and the target-polynomiality jet
that forces `gamma=1-3xy/2 mod x^2`; the written LND/Stein argument supplies
the slice under explicit intrinsic saturation labels.  On the reciprocal
side it checks the coefficient valuations `(n-1,2n-1)` and the extraction of
`Y=Q-Ps` from a primitive quadratic conormal coefficient.  The eight
minimal-boundary predicates are formalized in the accompanying note, but
this checker does not construct their finite-normalization witness, verify
`PC`, `NC`, or `CS` for an arbitrary boundary-minimal map, or extract a
suspension from the unmarked canonical normalization.

The same target also checks the finite-normalization frontend: the
Deligne--Faddeev cubic-algebra table and discriminant, the codimension-three
reflexive-module warning and its minimal excess-length-four special fiber,
whose exact module-theoretic defect is `Fitt_3=(x,y,z)`, the unique
critical-divisor DVR budget `(2,1)+(1,1)`, and the
tangent-hyperplane quotient coordinates.  The written local argument proves
that cubic point-flatness is equivalent to every canonical scheme fiber
having length three.  The cited nonflat triple-cover correspondence also
shows that normal cubic algebra structure alone cannot remove this defect.
The local structure theorem further identifies every defect with an
`(s+2)`-by-`s` determinantal presentation, where the excess fiber length is
exactly `s`; the checker includes the origin-primary `s=2`, length-five
rung in addition to the minimal Koszul rung.
For a reduced minimal defect it also verifies the linear-algebra inputs
forcing the unique square-zero fiber `k plus k^3`.  The written corollary
then identifies such a defect with a closed-point collision of the
ramified boundary sheet and the affine sheet over the critical divisor.
It distinguishes this from the allowed foundational collision, whose
triple-root fiber is curvilinear of length three.
The maximal-minor order argument proves that every reduced defect is
automatically this minimal Koszul rung; only nonreduced Fitting defects
remain outside the square-zero classification.
The local monogenicity theorem then closes all of those cases
simultaneously under intrinsic curvilinearity of the collision fibers:
Nakayama lifts a fiber generator and the resulting monic cubic algebra is
free.
Before either saturation test, the written tame-local proposition removes
every simple-normal-crossing point of the critical discriminant: after
strict henselization the cubic normalization is the finite-free sum
`R[s]/(s^2-t_1...t_r) plus R`.  Consequently the point-defect computation
first reduces to closed non-SNC points of the discriminant.  The written
ordinary-cusp proposition then classifies the two possible three-sheet
braid representations.  Equal meridian transpositions give the finite-free
`2+1` Kummer algebra, while distinct transpositions give the finite-free
monic cubic root cover.  Thus only worse-than-ordinary-cusp points remain.
For a reduced Koszul defect, the next written proposition identifies the
projectivized branch tangent cone with the discriminant of line sections of
the ternary cubic `h`.  The frontend checker verifies the complete
degree-six factor table for smooth, nodal, cuspidal, conic-plus-line,
triangle, and concurrent-line symbols, together with vanishing for double
and triple components.  Thus every reduced defect forces branch
multiplicity six, or at least seven in the non-squarefree case.
The frontend checker verifies that the foundational discriminant's singular
ideal is the expected triple-root locus.  The local cusp models and all nine
three-letter braid pairs are checked by:

```bash
.venv/bin/python plane-jc/cas/test_cubic_cusp_local_model.py
```

The remaining nonzero ternary-cubic symbol strata and both canonical
saturation modules are audited by:

```bash
.venv/bin/python scripts/verify_cubic_symbol_double_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_deformation_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_quartic_tangent_saturation.py
.venv/bin/python scripts/verify_smooth_cubic_quartic_plane_saturation.py
.venv/bin/python scripts/verify_singular_cubic_quartic_plane_saturation.py
.venv/bin/python scripts/verify_smooth_cubic_quartic_three_space_saturation.py
.venv/bin/python scripts/research_universal_cubic_quartic_kernel_saturation.py
.venv/bin/python scripts/verify_universal_cubic_filtered_syzygy_frontier.py
.venv/bin/python scripts/verify_cubic_quartic_ext_tail_absorption.py
.venv/bin/python scripts/verify_universal_cubic_quartic_different_complex.py
.venv/bin/python scripts/verify_universal_cubic_kahler_annihilator.py
.venv/bin/python scripts/verify_cubic_symbol_dense_quartic_plane_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_affine_dense_quartic_plane_saturation.py
.venv/bin/python scripts/verify_universal_cubic_cotangent_saturation.py
.venv/bin/python scripts/verify_cubic_formal_gauge_cokernel_atlas.py
.venv/bin/python scripts/verify_cubic_double_saturation_stratification.py
.venv/bin/python scripts/verify_nodal_cubic_formal_slice.py
```

The last checker also quotients the degree-five curvature by the
five-dimensional kernel of the quartic gauge lift.  It verifies a
rank-four action on the six slice--gauge coefficients, extracts two
intrinsic cross forms and three intrinsic pure-gauge quadrics, and checks
that the pure-curvature zero scheme is two reduced rational planes with
one embedded quadratic socle class.  Relative to the stored quartic lift,
it continues both reduced planes through degree six, constructs the exact
quadratic corrections, and obtains the two cubic Veronese classes.  It then
changes the plus-branch lift by `p` times the first quartic stabilizer and
proves that the resulting degree-six class changes literally while retaining
the same origin-only zero locus.  The full stabilizer-orbit classification
and embedded-socle continuation remain open.

For the homogeneous tensor, all seven squarefree strata have saturated
cotangent presentation and a length-six `Ext_A^2(T,A)` support defect with
Hilbert function `3+3t` (three-dimensional top and zero `m^2` action);
double and triple lines instead have a one-dimensional support defect.
The zero homogeneous tensor passes both module tests but is nowhere
generically étale.  One explicit order-four kernel tensor makes the support
defect finite of length six in all ten strata, while cotangent saturation
still passes.  This is an exact leading-model computation.  It neither
proves lift-independence nor constructs a normal lift with a Keller open.
The support presentation is computed directly as the module preimage
`modulo(H1,H2)`, namely the kernel of the four action columns in the
threefold direct sum of the cotangent quotient.  This is exactly equivalent
to extracting four coordinates from the combined syzygy module, but avoids
the unnecessarily large 97-column elimination.
The second command works over `Q[t,x,y,z]`.  On each of the seven
squarefree lines `phi_h+t*psi_4`, it verifies uniform cotangent saturation,
no parameter torsion in relative `Ext_A^2(T,A)`, radical support equal to
the collision axis, multiplicity six, and equality of the relative
presentation with the scalar extension of its central specialization.
This proves constancy on those lines, not on the full 24-parameter
order-four space.
The filtered-syzygy frontier command resolves the smooth central
unit-pruned cotangent presentation with cokernel ranks `7 -> 13 -> 6`,
checks that `x+y+z` is regular there, and then applies the 24 unchanged
central input syzygies to the universal 6-by-25 matrix.  Twelve exact
remainders survive modulo the central image.  This certifies that
entrywise collision-order growth and two-jet agreement do not themselves
give a coefficient-independent Rees-strict resolution.  It does not prove
boundary torsion or disprove the smooth universal saturation theorem,
which the formal-gauge command below now proves.
The quartic-tangent command is the longer tangent-direction audit.  It
tests every
one of the 24 exact nullspace-basis axes for every squarefree symbol.  All
168 families have uniform cotangent saturation, no parameter torsion,
collision-axis radical support, and relative multiplicity six.  Literal
presentation equality changes in four rows, without changing any of those
invariants.  The basis axes span the order-four kernel, but this computation
does not test all their linear combinations.
The fourth command is the four-worker mixed-direction audit for the smooth
symbol.  Over `Q[u,v,x,y,z]` it tests all 276 full coordinate planes.  On
every plane the cotangent presentation is saturated and the relative
length-six `Ext^2` presentation is pulled back from `u=v=0`; this includes
every specialization on the plane.  Directions supported on three or more
basis tensors are not tested.
The fifth command is the four-worker coordinate-plane audit for the six
singular squarefree symbols.  It verifies that 1,652 pruned presentations
are pulled back from the origin.  Four ambient presentations jump, so the
checker forms their exact finite `Q[p0,p1]` presentations using the
verified zero `m^2` action and proves `Fitt_6=(1), Fitt_5=(0)`.  Thus all
1,656 singular-squarefree planes are flat of relative rank six; the four
exceptional modules are in fact free by Quillen--Suslin.
The sixth command is the longer four-worker three-space audit.  It tests
all 2,024 smooth coordinate three-spaces over
`Q[p0,p1,p2,x,y,z]`.  After pruning contractible free summands, every
relative rank-three `Ext^2` presentation is pulled back from the parameter
origin with multiplicity six.  Directions supported on four or more basis
tensors are not tested.
The seventh command is an exact frontier calculation, not a universal
theorem.  It checks four full-support dense lines for every squarefree
symbol and the full first-ten-coordinate subspace for the smooth symbol.
All 28 lines have uniform cotangent saturation, no parameter torsion,
collision-axis support, multiplicity six, and central Ext presentation.
The smooth parameter ten-space has the central pruned rank-three
presentation.  It also constructs the full universal cotangent matrix,
checks that its parameter-dependent terms have bidegrees `(1,3)`, `(1,5)`,
and `(2,6)` in parameters/collision variables, and removes six
parameter-independent unit pivots to obtain a cokernel-equivalent
6-by-25 presentation.  This command is an exact input reduction, not a
saturation or Ext calculation over all 24 parameters; the final command
below supplies the universal theorem by a different method.
The eighth command computes the last nonzero differential in the minimal
support resolution on the seven full-support squarefree planes.  Six rows
are parameter-independent and linear; their transposes present a
length-six quotient killed by `(x,y,z)^2`.  The seventh row has a central
quadratic part and a parameter-linear cubic part, so it reduces to zero
modulo those six rows.  The resulting 12-generator parameter module has
six independent constant relations and therefore
`Fitt_6=(1), Fitt_5=(0)`.  It also verifies that the seven canonical
different generators equal the complete annihilator on each plane.
The ninth command constructs the canonical different matrix
`[(0,z,-y,x),(s_ij,2*mu_ij)]` over all 24 quartic-kernel parameters and
all seven squarefree symbols.  Its explicit universal syzygy matrix
satisfies the Buchsbaum--Eisenbud grade conditions.  The resulting
canonical-different support has constant length-six `Ext^2` and
`Fitt_6=(1), Fitt_5=(0)` over the full parameter ring.  Identifying these
seven generators with the complete annihilator `Ann(Omega)` is conditional
at this stage; the eighth command verifies that equality on the seven
full-support planes, and the final command closes it universally.
The same ninth command now checks the intrinsic form of the fixed tail.
On the exceptional plane it presents `Sym^2(Q)`, whose second Chern number
is six; its transpose vertex quotient has Hilbert layers `3+3` and is
killed by `(x,y,z)^2`.  This is the localized-Chern explanation of the
length-six block, independent of the squarefree cubic orbit.  The adjacent
binary-cubic discriminant class is `6H`, so branch multiplicity and support
length are the first- and second-Chern forms of the same universal quotient.

<!-- status-consumer: KLC6 88e899a9645c4a70 -->

The tenth command computes the universal Deligne--Faddeev locally free
cubic algebra and proves that the Kähler different `Fitt_0(Omega)` equals
the full annihilator `Ann(Omega)`.  On the punctured Koszul base this
identifies the canonical different with the actual support ideal.  Together
with the depth of the ninth command's exact complex, relative cotangent
saturation extends the equality across the collision axis and closes the
actual Ext Fittings.  The command does not prove that remaining universal
cotangent saturation; the final command does.
The eleventh command tests one low-height full-support plane for all seven
squarefree symbols.  If `psi_plus` is the sum of the 24 fixed kernel-basis
tensors and `psi_minus` their alternating sum, it computes the complete
family `phi_h+u*psi_plus+v*psi_minus` over `Q[u,v,x,y,z]`.  In every row
the cotangent presentation is saturated and the pruned rank-three relative
Ext presentation is pulled back from the origin with multiplicity six.
On `u^2-v^2!=0`, all 24 basis coordinates are nonzero.  This proves a
two-parameter full-support result, not a Zariski-open theorem in the
24-dimensional kernel.  Low-height dense parameter-three-space and
higher-height parameter-four-space runs reached their declared 600-second
timeouts and provide no mathematical evidence.
The twelfth command translates the same sum/alternating-sum plane by the
deterministic generic quartic lift and tests all nine nonzero cubic-symbol
orbits plus the zero symbol.  Over `Q[u,v,x,y,z]`, all ten rows have
saturated cotangent presentation, support exactly equal to the parameter
plane, relative Ext multiplicity six, and pruned rank-three presentation
pulled back from `u=v=0`.  The checker also reconstructs the generic lift
in the fixed primitive 24-element basis and verifies that every coordinate
is nonzero.  This proves uniform purity restoration for the double-line,
triple-line, and zero symbols on one affine plane; it does not prove
normality or Keller-open compatibility.

The final command proves the smooth-symbol 24-parameter theorem without
computing the universal saturation.  For the graded module `K` of all
compatible tensor corrections and the exact `10`-by-`9` simultaneous
coordinate/coefficient gauge matrix `G`, it verifies
`K=im(G)+A*eta` and `(x,y,z)*eta subset im(G)`.  An explicit matrix `L`
satisfies `G*L=[x*eta,y*eta,z*eta]`, so every compatible tensor term of
collision degree at least four is gauge.  It independently derives all
nine columns of `G` by expanding the determinant-twisted finite action over
the dual numbers.  It also stores an explicit linear-polynomial
`9`-by-`24` matrix `Q` with `G*Q=[psi_1,...,psi_24]`; both the quartic
compatible space and this gauge image have rank `24`, with gauge kernel
dimension three.  Successive homogeneous changes formally identify the
universal quartic family with its saturated central fiber.  Since
completion detects `(x,y,z)`-power torsion, this proves
`H^0_(x,y,z)(Omega)=0`; the canonical-different argument then gives the
universal annihilator equality and actual-support
`Fitt_6=(1), Fitt_5=(0)`.
The atlas command derives the determinant-twisted gauge differential over
the dual numbers for all ten ternary-cubic symbols and computes the exact
graded modules `ker(C)/im(G_h)`.  Their Hilbert series prove that smooth is
the unique symbol formally rigid above collision degree three.  The exact
quartic nongauge dimensions are `0`; `2,4,4,6,6,8` on the six singular
squarefree symbols; and `11,16,24` on the double-line, triple-line, and
zero symbols.  It also proves the exact singular-squarefree annihilator
sequence `(x),(x^2),(yz),(y^3),(xyz),(x^3)`; the three non-squarefree
quotients have zero annihilator and generic ranks `1,2,4`.  These data
delimit the formal-triviality method; they do not assert failure of
cotangent saturation.
The final nodal command proves the cyclic refinement
`ker(C)/im(G_nodal)=Q[y,z](-3)` with generator given by the tensor of
`Z^3`.  In quartic degree the 24-dimensional compatible space splits as a
22-dimensional gauge image plus the slice generated by `y*eta,z*eta`.
Only the first two fixed quartic basis directions survive in the quotient.
The sum/alternating-sum plane is a second transverse slice, with
change-of-slice determinant two.  Both slices replay the saturated
cotangent presentation and constant length-six Ext block.  This is a
first-stage slice theorem.  For the stored row-reduced quartic gauge lift,
the command also computes the complete degree-five normal curvature in the
basis `y^2*eta,y*z*eta,z^2*eta`: its components have `14,16,13` quadratic
terms and 30 nonzero cross-parameter pairs.  It vanishes on the coordinate
slice and has the explicitly recorded nonzero restriction on the dense
slice.  The command then quotients changes in the five-dimensional
gauge-lift kernel: its rank-four action on the six slice--gauge
coefficients leaves two intrinsic cross forms, while the three
pure-gauge quadrics are already invariant.  Their reduced zero scheme is
two rational planes, and their unreduced ideal has one embedded quadratic
socle class.  On both reduced planes the checker constructs an exact
quadratic correction of the degree-five term.  Its 15-dimensional
ambiguity acts trivially on the degree-six quotient, and the resulting
classes are `27/8*(q*y+p*z)^3*eta` and
`27/8*(q*y-p*z)^3*eta`.  Independence from the earlier quartic-lift
ambiguity and continuation of the embedded socle remain open.

The equivalent coordinate-free test is that each collision cotangent module
has unit first Fitting ideal (or vanishing second exterior power); the
checker separates the cyclic triple-root cotangent from the three-generator
square-zero cotangent.
The equivalent nilradical test has one generator and nilpotency index three
for the foundational collision, versus three generators and index two for
the reduced defect.
The written Hartogs extension theorem proves that a primitive cotangent
generator in codimension one extends through closed collisions whenever the
pure two-dimensional ramification support is `S_2` and its rank-one
cotangent module is `S_1`.  The companion two-`Ext` theorem identifies the
only closed-point obstruction modules as `Ext_A^2(T,A)` and
`Ext_A^3(Omega_{B/A},A)`.  Its double-saturation refinement forms the
canonical `S_2` hull `C=Ext_A^1(Ext_A^1(T,A),A)` and identifies those
obstructions successively with the canonical duals of `C/T` and
`Omega_{B/A}/T tau`.  The coupled local-cohomology sequence shows that
after `C=T` the latter is exactly the closed-point torsion of
`Omega_{B/A}`.  If `N` is the image of a free presentation and
`I=Fitt_3(B)`, the exact test is `N:I^infinity=N`; the Singular regression
checks this module saturation directly.
The phantom-boundary theorem identifies the quotient between the reduced
nonproperness and branch equations as the exact extra-divisor detector.
The checker calibrates it on the foundational map: boundary elimination and
the cubic discriminant give the same irreducible equation, so the quotient
is one.  The written boundary-minimality corollary then closes this
certificate for every boundary-minimal cubic: the foundational competitor
gives upper bound one for the number of target boundary components, while
nonproperness gives the matching lower bound.  Thus no second unramified
target divisor remains in the minimality problem.
This is not a global arbitrary-cubic closure.  Proposition 1.4 shows that a
second boundary sheet cannot lie over the critical divisor, because the
ramified `(2,1)` and affine `(1,1)` sheets exhaust degree three.  An
arbitrary cubic can still have a distinct unramified nonproperness divisor;
excluding that factor, or reducing it to the minimal stratum while
preserving genuine ungradedness, is the separate `OP-UG3` obligation.

The universal flat ungraded coefficient cell is checked separately by:

```bash
.venv/bin/python scripts/verify_universal_cubic_ungraded_testbed.py
```

It verifies the seven-parameter degree-at-most-four cell, and the written
argument extends the same identities to arbitrary
`A(P),gamma(P) in k[P]`.  The checker proves the determinant-minus-two
identity, inverse cubic and derivative reconstruction, reciprocal chart,
finite free Deligne--Faddeev multiplication table, discriminant, smooth
Laurent ramification parametrization, universal polynomial `GL_2`
discriminant transformation, and the exact equivalence
`G automorphic <=> phantom factor unit`.  It also separates this flat cell
from the 24-dimensional Koszul order-four tensor kernel by their intrinsic
third Fitting ideals.  It does not construct a Keller open for an arbitrary
Koszul-kernel combination.

The diagonal one-parameter orbit boundary of the displayed foundational
map is checked by:

```bash
.venv/bin/python scripts/verify_foundational_toric_degenerations.py
```

This exact weight-cone calculation leaves four faces: the foundational map
itself and three triangular automorphisms.  Its scope is diagonal source and
target weights in the displayed coordinates; it does not cover constant-
conjugate or nonlinear polynomial degenerations.

<!-- status-consumer: FTD3 f4b5bf44c04dba69 -->

The split one-coordinate foundational base-change calculation is checked
by:

```bash
.venv/bin/python scripts/verify_foundational_split_base_changes.py
```

For a coefficient map `(w,u,v) -> (f(w),u,v)`, it verifies the global
factor-space pullback equation, its smoothness and Laurent-UFD chart, and
the divisor-valuation presentation
`Cl(U_f)=Z^r/Z(m_1,...,m_r)`.  The written localization proof applies in
every degree and shows that the source is affine three-space only when `f`
is affine-linear.

<!-- status-consumer: SBC3 f5cbae00b4e87623 -->

The written no-global-monogenicity proposition then shows why these local
generators cannot be patched into one root coordinate: the derivative would
be a constant unit on `A^3` and would contradict cubic degree.
The written theorem proves uniqueness without a supplied suspension when
the intrinsic flatness defect is empty, the binary-cubic coefficient map is
affine-linear of full rank, and no extra simple boundary is omitted.
It also checks the nonlinear gauge-straightening theorem: every slice
`C_1=q-3C_0h` with `q!=0` and translation-invariant `h` is carried to
`C_1=q` by explicit polynomial source and target automorphisms.  It checks
the symmetric lower-unipotent family, the discriminant invariant, and the
variable-time Jacobian formula `1+D(h)`, which makes invariance necessary
for a single shear to be an automorphism.  The stress-test family
`C_1+tC_0^2=1` is verified directly to have source `A^3` and Jacobian `-1`
before being reduced to the foundational class.
The Borel corollary now exhausts every polynomial upper- or
lower-triangular `GL_2` gauge as well: its diagonal entries must be
constants, leaving exactly one classified invariant shear.
On the invariant coefficient hyperplane `C_0=0`, the same checker restricts
the time `h=4*C_1*C_3-C_2^2` to an explicit `A^3` automorphism, verifies its
inverse, determinant, and multidegree `(1,3,5)`, and checks the exact linear
conjugacy to the Nagata automorphism with parameter `-4`.  Wildness uses the
external Shestakov--Umirbaev theorem.  The first three swapped iterates are
expanded exactly and have multidegrees `(4s-3,4s-1,4s+1)`; this is a
known-family calibration, not a resolution of the open `(7,8,12)` case.
The alternating regression verifies the exact two-shear rank-two Jacobian
formula.  When the first time is invariant, conjugation reduces the second
factor to the single-shear theorem and gives an if-and-only-if transported
kernel criterion.  A Gröbner coefficient audit excludes every normalized
linear-time cancellation between two individually noninvertible factors.
The all-degree support theorem excludes every pair of nonzero monomial
times; the checker exhausts all 1,156 pairs through degree three.  Exact
graded ranks through degree eight give cokernel dimensions
`0,0,0,1,0,0,0,1`, confirming that the general recursive gauge equation has
only one discriminant obstruction in every fourth degree.  A second exact
checker parametrizes the ten-dimensional quadratic cancellation kernel and
proves that its degree-four discriminant projection vanishes identically;
all coupled basis directions admit recursive corrections through degree
eight.  The written `sl_2` divergence identity proves the bilinear
vanishing in every degree.  The ranked next attacks are recorded in
[`cancellation/CUBIC_CLOSURE_ATTACKS.md`](cancellation/CUBIC_CLOSURE_ATTACKS.md).

This runs the construction, parameter arithmetic, boundary, monodromy, and
current-ansatz rigidity regressions.  It includes the endpoint-moment
reduction of the cancellation contact resultant: the general triangular
identity is checked exactly on a bounded grid, while the complete
`r=1,2,3,4` columns are proved uniformly in `m`.  It also checks the
irreducibility transfer proving every `1<=m<=1000` column uniformly in `r` and
an explicit effective `r`-tail for each fixed `m`:

```bash
.venv/bin/python scripts/verify_parameter_irreducibility.py
.venv/bin/python scripts/verify_parameter_irreducibility_dusart_frontier.py
.venv/bin/python scripts/verify_parameter_irreducibility_sharp_dusart_frontier.py
.venv/bin/python scripts/verify_parameter_irreducibility_adaptive_dusart_frontier.py
.venv/bin/python scripts/verify_contact_resultant_irreducible_ranges.py
```

The second command is a slow four-process exact replay of the 2192 residual
pairs in `301<=m<=499`; it is kept out of the ordinary `verify-master`
target.  The third is a slower six-process exact replay of the 2899 residual
pairs in `500<=m<=741`; the fourth replays the 3335 adaptive residual pairs
in `742<=m<=1000`.  Both are likewise kept out of `verify-master`.

The `r=3` certificate checks
coefficientwise positivity of all six principal minors of the reciprocal
eliminant's Schur--Cohn matrix.  The heavier `r=4` certificate computes the
degree-eleven eliminant's `(9,2)` Schur--Cohn inertia, runs a 228-cell rational
Rouche localization, and proves the remaining argument separation by exact
angle and Bernstein-sign certificates.

The complete `r=5` column is a separate, substantially heavier exact replay.
It requires Singular for its boundary resultants:

```bash
.venv/bin/python scripts/verify_contact_resultant_r5.py
```

The formerly first open fixed-`r` column has an exact bounded-degree
reduction: the following Singular-backed checker constructs the
quintic--sextic endpoint equations and verifies that their residual eliminant
has degree 29 in `y` and degree 90 in `m`.

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_reduction.py
```

The branch-at-infinity replay then proves eventual nonvanishing in that
column.  It checks the complete Newton edge after `y=1+c/m`, the squarefree
degree-29 edge polynomial, and the linear reconstruction of the limiting
`z`.  Lindemann--Weierstrass separates algebraic `z` from `exp(c)`.  This
intermediate argument does not by itself provide an explicit threshold in
`m`.

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_asymptotic.py
```

For an effective certificate on the limiting edge itself, the following
atlas computes the first two `y` terms and first two `z` terms for all 29
branches (compressed to 15 complex-conjugation classes).  It uses 29 disjoint rational
Rouche disks and rational exponential enclosures to prove the strict modulus
gap `|z_0|!=|exp(c)|` branch by branch.  This limiting-edge atlas does not by
itself extract a positive tail threshold in `t=1/m`.

```bash
.venv/bin/python scripts/explore_contact_resultant_r6_branch_atlas.py
```

The effective continuation requires the pinned `python-flint` dependency in
`requirements.txt`.  It certifies 29 disjoint Rouche tubes on each of 256
rational cells covering `0<=t<=1/41`, separates the sixth-power identity by
modulus or phase on every tube, and checks the finite range `1<=m<=40`
modulo `1,000,003`:

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_effective.py
```

A separate bounded structural audit constructs the limiting endpoint systems
for `r=5,6,7,8`.  At `r=7,8` it verifies squarefree branch polynomials of
degrees 42 and 55, excludes `c=0` and `z=infinity`, and reconstructs a unique
finite `z`.  The first command is the cross-column limiting-system audit.
The next two commands construct the full bidegree-`(42,126)` and
`(55,200)` endpoint eliminants, identify their complete top Newton edges,
and prove eventual nonvanishing.  These three commands do not claim an
effective threshold or a continuation uniform in `r`:

```bash
.venv/bin/python scripts/verify_contact_resultant_fixed_r_branch_schema.py
.venv/bin/python scripts/verify_contact_resultant_r7_asymptotic.py
.venv/bin/python scripts/verify_contact_resultant_r8_asymptotic.py
```

The effective fixed-column template is replayed for `r=7,8` by

```bash
.venv/bin/python scripts/verify_contact_resultant_r7_effective.py
.venv/bin/python scripts/verify_contact_resultant_r8_effective.py
```

For `r=8`, the second command certifies 55 disjoint roots on each of 1,024
rational cells covering `0<=t<=1/1001`.  Its 56,320 Arb tubes split into
52,224 modulus and 4,096 phase separations.  It then uses FLINT over
`GF(1,000,003)` for 1,000 degree-preserving endpoint gcd certificates,
closing every integer `m>=1`.  The shared checker is configuration-driven
for `r=6,7,8`; subsequent columns reuse the endpoint-chart construction,
Rouche tubes, logarithmic separation, and finite-field completion after
supplying their exact degree/edge data and a certified partition.

The additional finite `5<=r<=12` endpoint grid is quick to replay.  It checks
203 monic gcd certificates modulo `1,000,003`, including denominator and
leading-coefficient unit conditions:

```bash
.venv/bin/python scripts/verify_contact_resultant_modular_grid.py
```

It also runs the log-geometric bridge regression, including the reciprocal
determinant, canonical Jacobian-LND exponent, the degree-two plinth/Stein
countermodel, spectral squarefreeness, and Laurent-tail descent.  The reusable
classifier additionally checks exact prime valuations, both localized chart
compositions, boundary elimination, the displayed residue degree, the full
Stein field via local-slice invariantization, hidden covers, and the spectral
gcd obstruction.  It also checks the unsliced divided-difference Hensel
multiplier which upgrades the boundary value to the complete cancellation
jet and global slice.  Its built-in examples can be inspected directly:

```bash
.venv/bin/python scripts/classify_reciprocal_link.py cancellation
.venv/bin/python scripts/classify_reciprocal_link.py masuda
.venv/bin/python scripts/classify_reciprocal_link.py masuda-hidden
```

The arithmetic portion also checks the fixed-row Newton-ramification
extraction:

```bash
.venv/bin/python scripts/verify_fixed_r_newton_ramification.py
```

It verifies the reciprocal numerator and prime-power congruence on a bounded
grid and exact cyclotomic-cluster Newton edges for derivative orders one
through eight.  The analytic density estimate is the cited external theorem
input, not a finite computation.

## External quartic islands

Juntang Zhuang's pinned `F4a`, `F4b`, and `F4c` examples have an independent
compact reconstruction and canonical-boundary audit:

```bash
.venv/bin/python scripts/verify_external_quartic_islands.py
```

This command is also part of `make verify-regressions`.  It requires no
network access and does not copy or execute the upstream checker.

## Decorated normalization, affine-mark faithfulness, and Hurwitz--LL calculations

The LL critical-value incidence, low-pole filtration, contravariant
triangular target convention, affine pencil transport, higher-zero Newton
polygons, nonzero multiple-root collisions, and normalized rerooting
identities are checked exactly by

```bash
.venv/bin/python scripts/verify_stable_generator_rigidity.py
.venv/bin/python scripts/verify_generic_affine_mark_faithfulness.py
.venv/bin/python scripts/verify_weighted_stable_moduli_chart.py
.venv/bin/python scripts/verify_quartic_biprojective_graph.py
.venv/bin/python scripts/verify_quartic_rees_stratification.py
.venv/bin/python scripts/verify_quartic_rees_corner_atlas.py
.venv/bin/python scripts/verify_quartic_rees_side_corner_matching.py
.venv/bin/python scripts/verify_intrinsic_selector_attack.py
.venv/bin/python scripts/verify_hasse_typical_seed_recovery.py
python3 scripts/verify_positive_characteristic_deformation_landscape.py
.venv/bin/python scripts/verify_multicluster_ll_comparison.py
.venv/bin/python scripts/verify_labelled_node_saturation.py
.venv/bin/python scripts/verify_branch_wonderful_pullback.py
.venv/bin/python scripts/verify_source_vertex_rigidity.py
.venv/bin/python scripts/verify_general_radial_source_atlas.py
.venv/bin/python scripts/verify_polynomial_monodromy_forests.py
.venv/bin/python scripts/verify_monodromy_inertia_characters.py
.venv/bin/python scripts/verify_recursive_resonance_atlas.py
.venv/bin/python scripts/verify_h1_h2_comparison_obstruction.py
.venv/bin/python scripts/verify_branch_scale_fan.py
.venv/bin/python scripts/verify_degree_six_branch_target_graph.py
.venv/bin/python scripts/verify_degree_six_admissible_equal_scale.py
.venv/bin/python scripts/verify_degree_six_admissible_radial_atlas.py
.venv/bin/python scripts/verify_degree_six_admissible_maxwell_atlas.py
.venv/bin/python scripts/verify_degree_six_central_hurwitz_selection.py
.venv/bin/python scripts/verify_degree_six_stack_inertia.py
.venv/bin/python scripts/verify_degree_six_stacky_fan_descent.py
.venv/bin/python scripts/verify_rerooting_groupoid_boundary.py
.venv/bin/python scripts/verify_coarse_affine_mark_descent.py
.venv/bin/python scripts/verify_restricted_ll_degree.py
.venv/bin/python scripts/verify_caustic_maxwell_boundary.py
```

The quartic biprojective-graph command requires Singular `4.4.1`.  It first
checks the source base triangle and generic line orders `(6,4,2)`, checks
two radical length-`28` sections exactly over `Q`, then separately computes
the saturated graph dimension, minimal-generator count, and multidegrees
over `F_32003`; the latter is a modular calibration, not a
characteristic-zero normalization certificate.

The Rees-stratification command is an exact SymPy calculation over rational
function fields.  It computes the complete generic-side base-point clusters,
the integral-closure colengths `(43,19,17)` and defects `(9,6,1)`, the six
compact vertex Newton facets, and the quadratic-transform principalization
of all six generic face colors.

The Rees-corner command completes the finite vertex atlas.  It verifies the
eight maximal Newton cones and their unit witnesses, the two `P_18` boundary
branch packets and conductor colengths `(12,15)`, the residual rational and
quadratic cusp conductors of colengths `(1,2)`, their explicit quadratic
normalizations, and the smooth colored `V_23` overlaps.  It requires
Singular `4.4.1` for four exact conductor calculations.  It does not glue the
closures of the three generic side clusters, compute the global bigraded
Hilbert polynomial, or prove simultaneous normalization.

The side-corner matching command compactifies all rational and algebraic
infinitely-near center curves over the three side parameter lines.  It finds
the three nontrivial packet bidegrees `(4,4)`, `(6,2)`, `(6,2)`, verifies
endpoint conductor colengths `(6,3,3)` at both ends and the corresponding
packet Hilbert corrections, and checks the exact `V_12`, `V_13`, and `V_23`
transitions.  It requires Singular `4.4.1` for six curve-conductor
calculations.  It does not assemble the rational exceptional surfaces or
compute the global normalized-Rees Hilbert polynomial.

For at least eight labelled finite branch values, the target receiver in the
wonderful-pullback construction is `Mbar_0,b+2` with `b+2>=10`.  Its Cox ring
is not finitely generated by the externally cited non-polyhedral
pseudo-effective-cone theorem.  This corollary has no local checker and does
not assert the same statement for the normalized pullback graph.

These checks support the generic affine-mark faithfulness theorem: the coarse
fiber is the exact rerooting orbit and every nontrivial rerooting moves the
unique unramified affine sheet into the reconstruction boundary.  The
Hasse-typical checker separately proves the sharp
`floor(log_p(N))+1`-channel coefficient repair in positive characteristic
and replays the degree-eight `F_5` collision.  It also verifies a clean
five-member degree-twelve `F_5` family whose distinct seeds have identical ordinary
derivatives and therefore define the identical weighted polynomial map,
proving that the channels cannot be made map-intrinsic without enriching
the construction.  The same checker then constructs the dimension-preserving
correction `A -> A-K(W)/(cC^2)`, verifies polynomiality, determinant one, and
all five intended inverse pencils, and certifies the marked transverse node
in the `c=2` member.  Five distinct reduced equal-image Groebner bases then
prove that the enriched maps are pairwise stably left--right inequivalent.
The written theorem upgrades this example to every odd-characteristic tame
clean degree by recovering the complete primitive-root factor from the
intrinsic second-boundary edge data.  It also records the characteristic-two
parity reconstruction, the identically singular old suspension parameter,
and the scalar-ansatz no-go theorem.  The checker then verifies the
weight-redistributed replacement
`u=1+x^2y`, `gamma=1+xz`: its coordinates are polynomial, its Jacobian is
one, and its inverse pencil is the prescribed normalized seed.  The final
characteristic-two block verifies the radicial discriminant factor
`W^2-T` and explicit squarefree, compressed-birational wild-clean witnesses
in every degree from five through sixteen; the formulas prove the resulting
stable-faithfulness theorem uniformly for all `N>=5`.  It also verifies the
complete symbolic quartic slice
`(1+lambda)W^2+W^3+lambda W^4`; the two affine marks on the radicial edge
remove the former low-support ambiguity, while the normalized cubic is
unique.  The same weight redistribution is checked in characteristics
three, five, and seven on the universal singular-parameter quartic
`2W^2-3W^3+W^4`, confirming that it complements the original chart exactly
on `2+H''(1)=0`.  The full-edge theorem then removes the old Hessian
degree/support restriction on every boundary-clean generically birational
locus.  Finally `d(WH'-H)=W dH'` proves that critical birationality is
automatic for every odd-characteristic exact-double seed; the
characteristic-two checker exhaustively regresses the corresponding
clean-implies-compressed-birational lemma through degree twelve.
Repeated-root examples in characteristics two, three, and five then verify
that the normalized second-boundary prime retains the complete primitive
root divisor with multiplicities, even when its critical image collides
with the zero cluster.  This is the executable collision regression for the
full theorem: the smaller marked-edge quotient
`(A^1_W;(W),(W-1),div(H/(W^2(W-1))))` reconstructs the normalized seed
exactly on every declared stratum.
The
selected root extends on the marked corrected graph, and the
normalized-Stein, completed-chart, and conductor comparisons are complete at
arbitrary simultaneous collisions.  Coarse affine-mark descent is also
complete over that graph: the marked invariant ring is the universal
monic-root incidence,
and the total-collision fiber `k[T]/(T^mu)` has one geometric point.  The
specialized restricted-LL
degree and caustic/Maxwell boundary-class calculations have no recorded
external review.

The companion affine-stratum audit verifies that the root-one component is
regular and that a nontrivial rerooting sends it to an extra-root boundary
component.  The multicluster audit checks distinct tangent lines, all pairwise
intersection numbers, the conductor exponent
`e_i(sum_j e_j-1)`, and regularity of the full marked-root incidence at
collisions.  The H1/H2 obstruction checker recovers the degree-five
`(x^3,y^2)` normalized blowup.  The branch-scale checker then computes the
degree-six `(2,2,2)` moving critical values, all six weighted braid-fan
chambers, and a triple-resonance cross-ratio proving that the radial fan is
only the first layer of the full logarithmic graph.  The wonderful-pullback
checker enumerates the complete `Mbar_0,n` boundary building set and maximal
nested sets for four through seven target marks, verifies permutation
equivariance, and recovers both the degree-five weighted blowup and the
degree-six six-line/four-center target from that one construction.  The
source-vertex checker exhausts 2,024 zero/pole divisor profiles in degrees
one through seven and proves computationally that two fibers reconstruct a
rational component map up to scale while one third-flag point fixes the
scale.  The general radial-source
checker then verifies the connector/local-polynomial-tail/identity-strand
rule for 780 multiplicity profiles and 48,580 ordered scale types, including
all component degrees, Riemann--Hurwitz identities, node partitions, lcm
saturations, label permutations, and independent dynamic verification of
the full-chain inertia formula.  It finds 42,158 nontrivial
unequal-multiplicity types in this range; equal multiplicities remain
trivial.  The monodromy-forest checker then
exhausts all 1,441 reduced polynomial transposition factorizations through
degree six and proves that every nested resonance source tree and node
partition is the corresponding edge subforest; it recovers pairwise
Maxwell, triple Maxwell, and caustic nodes from one rule.  The recursive
resonance-atlas checker then verifies framed residue coordinates on all 534
nested families with two through five branch labels, all 534 affine gauge
changes, 1,453 one-step and
2,926 two-step contractions, normalized flag equations through degree seven,
84 source/target frame transitions, 63 nonfactorized smoothing families,
automatic tame character extraction, 76 bounded full-centralizer radial
charts, all 89 degree-six
interval-nested families, and the order-four pair--triple inertia.  This
closes the former
explicit-stack gap.  The finite-normalization theorem
uses finiteness of the fully marked admissible-cover branch morphism to prove
that the normal wonderful graph is already the complete coarse source graph;
no additional source-side coarse blowup is possible, and corrected H2/H3
are unconditional.  The monodromy-centralizer checker computes all
polynomial tree deck groups through degree six, all cyclic connector groups
through degree eight, and anchored/unanchored inertia on every collision
node.  The recursive checker corrects the full-chain radial calculation:
equal multiplicities have trivial inertia, while an ordered partition
\(B_0|\cdots|B_k\) of arbitrary multiplicities has order
\(\prod_jL_j/M_j\); it checks 76 bounded equal and unequal charts.  The
centralizer checker gives one generic formula covering Maxwell and caustic
resonance.  The
complete-target
checker identifies the radial target with the three-coordinate-point blowup
of `P^2`, the stable target with the additional diagonal-point blowup
`Mbar_0,5`, and its source pullback with four reduced triple-Maxwell
branches.  The equal-scale admissible checker constructs the central
degree-six component and three quadratic tails, verifies all
Riemann--Hurwitz counts, and proves that the three index-two source nodes
normalize into exactly the same four Kummer branches.  The radial-atlas
checker then enumerates all thirteen ordered scale types, verifies degree six
and Riemann--Hurwitz on every target-bubble preimage, and checks every
node-index partition and Kummer saturation count.  The Maxwell-atlas checker
handles all three pairwise collision divisors and the triple collision,
matches their two- and four-branch source-node normalizations, and proves
that their residual radial intersections are transverse while their
coordinate intersections are already radial equality faces.  The central
Hurwitz-selection checker finds two ambient degree-six cover classes with
the required profiles, then proves by an exact square-cubic branch invariant
that the labelled source-root cross-ratio selects the polynomial class as a
reduced local branch.  The stack-inertia checker separates normalization
branches from genuine label-preserving cover inertia: every radial lift in
the equal-multiplicity degree-six chart has trivial inertia, while pairwise
and triple Maxwell lifts each retain one diagonal `mu_2`.  The stacky-fan
checker constructs the four-divisor Maxwell
root complex, proves all pair--triple face inclusions and `S_3` equivariance,
computes the four radial quotient orbit types, verifies the pair--triple and
radial--Maxwell inertia ranks needed for smooth tame-stack reconstruction,
and keeps the local
`(S_2)^3 semidirect S_3` pair-root stabilizer separate.  The
general labelled-node checker exhausts 1,554 index profiles, proves the
phase-quotient and label-preserving inertia formulas, checks permutation
equivariance, and verifies that the corrected marked/unmarked quotient over
any labelled normalized graph has degree `N-2`.  This makes label gluing and
the finite H2 factor formal over the `H1-COARSE` graph, independently of the
substantially stronger `H1-STACK` theorem.  The
rerooting-groupoid
audit separately checks the quotient degree `N-2`,
the selected-in/selected-out boundary pullbacks, generic transposition
ramification after coefficient contraction, and the distinction between a
cyclic total-collision slice and generic divisor inertia.  These three audits
and the companion affine-stratum audit are part of `make verify-regressions`.
The restricted-LL audit checks the Cayley/marking count and independently
computes degrees `8` and `75` from the quartic and quintic critical-value
eliminants.  The caustic--Maxwell audit checks the unique invariant Keel
relation, every collision and infinity valuation, both boundary
presentations, and the exact factorization `LL-discriminant=C^3 M^2` in
degrees four and five.  All displayed commands are part of
`make verify-regressions`.

## External consequence identities

Christopher D. Long's direct Gaussian-moment, `(xz)`, `SU(2)`, and `SO(3)`
identities, together with the exact normalization of the foundational map
used in his BCW discussion, have a dedicated target:

```bash
make verify-external-consequences
```

The Gaussian, `(xz)`, spherical `SO(3)`, and algebraic Haar scripts use only
the Python standard library.  Their bounded exact regressions are
distinguished from the all-exponent proofs in the canonical notes.  The
`SO(3)` replay checks the displayed moments through order fifteen and the
endpoint-jet identity through order one hundred.  The two algebraic checkers
verify the unique normalized functional on `UV+T^2=1`, its three
infinitesimal `so3` identities, the factorial functional on
`k[SL2]`, all six left/right `sl2` identities, and the explicit
`SL2/T` pullback.  The proof and the quotient/transfer theorem are in
[`ALGEBRAIC_HAAR_QUADRIC_AND_SL2.md`](extended-geometry/ALGEBRAIC_HAAR_QUADRIC_AND_SL2.md).
A separate symbolic checker proves the `SU(2)=S^3` Haar density in Hopf
coordinates, retaining an independent compact integration proof.  The same target
also performs all 18 balanced BCW steps and checks the resulting 79-variable
cubic-homogeneous collision, writes its sparse artifact, and replays it with a
separate standard-library implementation.  It then runs the shared-factor
optimization, which introduces 13 variables, reaches degree three in
dimension 16, and writes and replays a 33-variable baseline artifact.  It then
computes the exact rational rank 7 of the cubic component vector, constructs
the rank-compressed 24-variable cubic collision, and independently replays
the factorization, sparse map, and collision using only the standard library.
It then removes the two-dimensional constant Jacobian kernel, constructs the
22-variable quotient, and independently replays `BK=0`, `BC=I`, `H=HCB`,
cubic homogeneity, the descended collision, and the triangular determinant
factorization using only the standard library.
Finally, the essential-dimension search freezes a different 17-dimensional
trace of cubic-output rank six, homogenizes it in 24 variables, removes its
three-dimensional constant kernel, and independently replays the resulting
21-variable collision from the original map using only the standard library.
The backward-cubic continuation keeps the nonhomogeneous and homogeneous
dimension objectives separate.  It audits MacFarlane's displayed `F13` and
`G20`, restricts the sole fixed covector `tau` at the collision level, and
verifies the exact stable factorization
`M19=A_B o (F13 x I_6) o S_gamma`:

```bash
make verify-backward-cubic-reduction
```

The generated records are
[`macfarlane_g20_dimension_reduction_audit.json`](artifacts/generated-results/macfarlane_g20_dimension_reduction_audit.json),
[`macfarlane_f13_low_degree_invariants.json`](artifacts/generated-results/macfarlane_f13_low_degree_invariants.json),
and
[`backward_cubic_reduction_calibration.json`](artifacts/generated-results/backward_cubic_reduction_calibration.json).
The same target applies the two backward objectives and the pair-aware
collision policy to the retained restricted-minima archives, reconstructs
two current representatives exactly, and writes
[`backward_cubic_current_applications.json`](artifacts/generated-results/backward_cubic_current_applications.json).
The generic calibration also proves that the parent is isotrivial over
`t!=0`, its `t=0` fiber is triangular and injective, and every parent
collision can therefore be normalized to `t=1`; an exact `t=2` MacFarlane
collision is replayed as a regression.
It also runs the established `16 -> 24` rank-compressed BCW route with the
new reverse-companion regression enabled.
With the pinned external determinant certificate, the same audit updates the
external-certificate frontiers to `n_cub<=20` and, by homogeneous cotangent
lift, `n_HN,4<=40`; the internal dependency-free replay endpoints remain 21
and 42.
These commands calibrate the backward compiler and close stated direct
linear/degree-at-most-three routes; they do not construct a twelve-variable
map.

The next coordinate-pair reduction goes beyond pullback-fixed invariants.
It uses `s=F13_13=x13+x2^2` as a source coordinate and the target square
completion `y4-y8^2`.  The resulting exact relative form restricts at
`s=0` to a 12-variable degree-three Keller collision.  A direct sparse
determinant expansion and a separate standard-library implementation replay
the theorem.  Its cubic-output rank is six, giving a 19-variable
cubic-homogeneous parent and the updated bounds `n_cub<=19` and
`n_HN,4<=38`:

```bash
make verify-macfarlane-f12
```

The generated record is
[`macfarlane_f12_coordinate_pair_reduction.json`](artifacts/generated-results/macfarlane_f12_coordinate_pair_reduction.json).

The first exact continuation toward eleven variables classifies every linear
target coordinate whose pullback is a polynomial graph coordinate.  The raw
degree-three coefficient ideal is the unit ideal in all nine possible pivot
families.  At the literal triangular coordinates, every graph deletion has
degree four or five and at least one high-degree defect lies outside the
complete degree-at-most-three target-shear span in the other raw retained
outputs.  The two closest literal cases remain outside that span through
target degree four:

```bash
make verify-k12-coordinate-pair-frontier
```

The generated record is
[`k12_coordinate_pair_frontier.json`](artifacts/generated-results/k12_coordinate_pair_frontier.json).
This is a bounded obstruction, not a dimension-eleven lower bound; nonlinear
source coordinates and ordered multi-stage target automorphisms remain open.

The parameterized continuation then retains every linear target coordinate
whose pullback has a quadratic graph. Fixed full-column and augmented
minors, together with unit-ideal covers in the graph parameters, exclude
quadratic target completion for all six pivot families and cubic target
completion for all five single-defect families:

```bash
make verify-k12-parameterized-completion
```

The generated record is
[`k12_parameterized_completion_frontier.json`](artifacts/generated-results/k12_parameterized_completion_frontier.json).
This remains a bounded theorem. The multi-defect `z8` cubic completion is
handled by the next command; cubic graph corrections and ordered target
stages remain outside the combined scope.

The remaining multi-defect `z8` cubic system is assembled without a full
fraction-field expansion. Sparse modular elimination selects three minors,
which are then reconstructed exactly over the rational parameter ring.
Their determinant opens generate the unit ideal and every augmented
determinant is `9/7` times its column determinant:

```bash
make verify-k12-z8-cubic-completion
```

The generated record is
[`k12_z8_cubic_completion_frontier.json`](artifacts/generated-results/k12_z8_cubic_completion_frontier.json).
Together with the preceding command, this excludes one-stage cubic target
completion for all six quadratic graph-coordinate families.

The same sparse compiler extends through target degree four on every
single-defect family. Each parameter family has a nonzero constant
`990 x 990` column minor and a nonzero constant augmented minor:

```bash
make verify-k12-single-defect-quartic-completion
```

The generated record is
[`k12_single_defect_quartic_completion_frontier.json`](artifacts/generated-results/k12_single_defect_quartic_completion_frontier.json).
Only the much larger multi-defect `z8` quartic family remains in this graph
class.

The cross-construction audit compares the public dimension-38 route with the
independent `K12` route.  It checks the shared compressed cost `n+r=18`,
obstructs all seven source-affine linear/quadratic pivot completions of the
public eleven-variable lift.  On `K12`, it also obstructs quadratic target
completion of the nonlinear `z8` pivot and finds a fourteen-parameter
family of coordinated degree-preserving quadratic source shears.  A fixed
minor and an inconsistent exact Schur system of ranks `(5,6)` prove that
no member of that family lowers the cubic-output rank from six to five:

```bash
make verify-hvc38-cross-frontier
```

The generated record is
[`hvc38_cross_construction_frontier.json`](artifacts/generated-results/hvc38_cross_construction_frontier.json).
This is a bounded frontier computation, not a lower bound for quartic HVC.

The next gap-closure audit uses the square identities at the public `d`
pivot and local `z8` pivot to reduce nonlinear completion to a filtered
pullback calculation.  Good-prime ranks and matching exact kernels exclude
both pivots through target degree eight.  It then combines 140 quadratic
source columns with 792 elementary quadratic target columns.  The
36-dimensional exact high-degree kernel contains seventeen directions that
integrate to genuine triangular one-parameter source-target families.  On
their combined seventeen-parameter degree-three locus, the ideal obtained by
adjoining a selected cubic rank-six minor has Gröbner basis `[1]`:

```bash
make verify-hvc38-gap-closure
```

The generated record is
[`hvc38_gap_closure.json`](artifacts/generated-results/hvc38_gap_closure.json).
This excludes only the stated bounded pivot algebras and quadratic
source-target family; it does not prove minimality at dimension 38.

The maximal-block continuation enumerates all six maximal jointly affine
source blocks of `K12`.  For each block it combines every complementary
quadratic source shear with all 792 elementary quadratic target directions,
lifts the complete good-prime high-degree kernel over `QQ`, and verifies a
linearized rank-six Schur witness.  It then integrates every kernel
direction into one full triangular source-target family—including
source-only directions and directions that fail to preserve degree three
individually—and asks whether the exact degree-three locus can have
cubic-output rank at most five.  Pinned packets of at most 32 cubic minors
give unit ideals in Singular:

```bash
make verify-hvc38-maximal-block-closure
```

The generated record is
[`hvc38_maximal_block_closure.json`](artifacts/generated-results/hvc38_maximal_block_closure.json).
This closes the full quadratic left-right kernel class on all maximal
jointly affine blocks.  It remains a bounded theorem, not a dimension-38
minimality result.

The tensor continuation computes both natural coefficient flattenings of
the `K12` quadratic and cubic tensors and of the cubic-homogeneous `G19`
tensor.  Exact rational row reduction gives input-directional ranks `12`,
`12`, and `19`; hence all three common right kernels are zero.  The `G19`
output rank is `18`, with sole left annihilator the fixed `tau` output.  The
same checker replays both collisions and the companion scaling identity:

```bash
make verify-k12-tensor-module-frontier
```

The exact generated record is
[`k12_tensor_module_frontier.json`](artifacts/generated-results/k12_tensor_module_frontier.json).
This excludes constant linear tensor quotients of the displayed maps and
gives pure-cube decomposition lower bounds `12` and `19`; it does not
exclude nonlinear graphs, nonconstant modules, Schur elimination, or a
different tensor.

A separate finite-field scout enters the larger linear-coordinate families
whose graph corrections are genuinely cubic.  It searches all parameter
supports of size at most two with values in `{-2,-1,1,2}`, plus 250
deterministic random points per parameter count, over both `GF(101)` and
`GF(103)`.  Every bad retained output is tested against all 10 linear and 55
bilinear target monomials in the other raw outputs:

```bash
.venv/bin/python scripts/search_k12_cubic_graph_bilinear_completions.py \
  --support-max 2 --values=-2,-1,1,2 --random-samples 250 \
  --random-seed 20260804 --primes 101,103 --keep-closest 5 \
  --output artifacts/generated-results/k12_cubic_graph_bilinear_modular_search.json
```

The 15,688 modular evaluations contain no survivor at either prime.  The
generated
[`k12_cubic_graph_bilinear_modular_search.json`](artifacts/generated-results/k12_cubic_graph_bilinear_modular_search.json)
is a bounded discovery experiment, not a rational obstruction or a
dimension-eleven lower bound.

The exact continuation lifts the modular echelon rows across all nine
complete cubic graph-coordinate parameter spaces.  It combines constant
full-column/augmented minors with rank-stratified determinant covers over
`QQ`:

```bash
make verify-k12-cubic-graph-bilinear-obstruction
```

The generated record is
[`k12_cubic_graph_bilinear_obstruction.json`](artifacts/generated-results/k12_cubic_graph_bilinear_obstruction.json).
Five families have constant determinant ratios `-3,-9,-1,-1/2,-1/2`.
For `z5,z6,z7,z8`, the checker verifies exact determinant-open covers and
the residual closed strata, using explicit column relations for `z5,z6,z7`
and the complete BCR5 quadratic-family obstruction for the closed `z8`
stratum.  Together these certificates close all nine normalized linear graph
families through one-stage target degree two.

The cubic-target continuation first quotients by the certified bilinear
column space and adds the 220 cubic target monomials.  Its bounded two-prime
discovery pass uses all parameter supports of size at most one with values
in `{-2,-1,1,2}` plus 250 deterministic random points per parameter count:

```bash
make search-k12-schur-cubic-completions
```

All 2,902 modular evaluations have augmented-rank increment one and there is
no survivor.  The exact continuation is:

```bash
make verify-k12-graph-cubic-completion-obstruction
```

It reconstructs the complete 285-element target basis over the rational
parameter rings.  The pivots `z4,z8,z9,z10,z11,z12` have constant
full-column and augmented minors.  Determinant-open covers plus exact
closed-stratum relations handle `z5,z6,z7`.  Hence all nine normalized
linear graph families are obstructed through one-stage target degree three.
The generated records are
[`k12_schur_cubic_completion_modular_search.json`](artifacts/generated-results/k12_schur_cubic_completion_modular_search.json)
and
[`k12_graph_cubic_completion_obstruction.json`](artifacts/generated-results/k12_graph_cubic_completion_obstruction.json).

The restricted-minima continuation then changes the BCW circuit before
homogenization: it exposes two polynomial gates, cancels one complete
multi-term circuit block, and scores every partial trace by its Jacobian
power-rank profile.  The two frozen winners are a 22-variable
cubic-homogeneous collision of exact index 18 and a 24-variable collision of
exact generic rank 17 and index 18.  Singular certifies their generic kernel
dimensions, while independent standard-library audits multiply the full
polynomial Jacobians and verify `(JH)^17!=0`, `(JH)^18=0`.  A fifth grouped
atom cancels the first-coordinate circuit `x^2(3y+xz)`.  The expanded
32-family Pareto search finds a separate 22-variable cubic source whose
44-variable HN lift has exact generic Hessian rank 37:

These computations form the now-frozen upper-bound track.  The commands below
reproduce the recorded search; they are not an active broad-search queue.
If the program is reopened, it should begin with a theorem-directed
five-dimensional cubic classification or the invertibility-only question for
arbitrary cubic-homogeneous Keller maps with `(JH)^3=0`.

```bash
.venv/bin/python scripts/search_restricted_bcw_circuits.py \
  --width 64 --max-steps 24 --prebeam-factor 2 --partial-power-depth 8 \
  --skip-terminal-hessian-power \
  --enable-atom x2s --enable-atom v2r --enable-atom qb \
  --enable-atom v2h --enable-atom y2vb \
  --output artifacts/generated-results/restricted_bcw_circuit_search_v2_w64.json
.venv/bin/python scripts/search_rank37_gate_perturbations.py \
  --width 16 --max-steps 10 --prebeam-factor 3 \
  --partial-power-depth 8 \
  --output artifacts/generated-results/rank37_gate_perturbation_search.json
.venv/bin/python scripts/verify_index_reduced_bcw_22_route.py
python3 scripts/audit_index_reduced_bcw_22_independent.py
.venv/bin/python scripts/verify_rank_reduced_bcw_24_route.py
python3 scripts/audit_rank_reduced_bcw_24_independent.py
.venv/bin/python scripts/verify_hessian_rank_reduced_bcw_22_route.py
python3 scripts/audit_hessian_rank_reduced_bcw_22_independent.py
.venv/bin/python scripts/verify_hessian_rank_35_identity_slice.py
python3 scripts/audit_hessian_rank_35_identity_slice_independent.py
.venv/bin/python scripts/search_identity_slice_hessian_rank.py
.venv/bin/python scripts/search_identity_slice_local_perturbations.py
.venv/bin/python scripts/verify_hessian_rank_34_double_identity_slice.py
python3 scripts/audit_hessian_rank_34_double_identity_slice_independent.py
.venv/bin/python scripts/verify_index_three_inverse_model.py
.venv/bin/python scripts/verify_index_three_degree_bound_counterexample.py
.venv/bin/python scripts/derive_index_three_tree_obstruction.py
.venv/bin/python scripts/verify_restricted_minima_frontier.py
```

The second index-three command replays van den Essen's dimension-five
generic-rank-three automorphism, proves `(JH)^3=0`, verifies both inverse
compositions, and extracts the nonzero degree-eleven and degree-thirteen
terms.  The tree command independently evaluates the degree-eleven normal
form on the same tensor.  Together they disprove the proposed uniform
inverse-degree-nine bound while leaving the full-class invertibility-only
question open.

The resolved two-real theorem, the first classified three-real minimality
island, and the cross-conjecture minimum ledger have their own fast exact
target:

```bash
make verify-counterexample-scoreboard
```

The generated
`artifacts/generated-results/minimal_counterexample_scoreboard.json` records
the exact unrestricted-GVC failure dimension as three.  Its current
whole-file SHA-256 is
`09df0e398def5df799243c906066f0b469b17ccf63f7d9261e8944a96fe8f8b1`.

This proves GMC for every quadratic Gaussian polynomial in every dimension,
checks the two-weight and affine-circular-source obstructions in two real
variables, exactly excludes all 27 mixed-sign cubic three-weight supports on
their 72 nonvanishing charts using moments through order eight, excludes 29
of the 33 mixed-sign cubic four-weight supports on 97 charts using moments
through order six, and then excludes all four charts of the symmetric
exceptional support by three good-prime quotient-algebra certificates: in
each representative the tenth moment acts with rank 84 on the
84-dimensional order-eight quotient, and circular-coordinate reflection
supplies the fourth chart.  Seven further exact rational unit-ideal
calculations exclude the last three supports and 20 charts through moment
six.  Thus all 121 mixed-sign four-weight cubic charts are closed and a
cubic GMC(2) counterexample needs at least five rotational weights.  The
target also proves the Bessel--factorial moment formula for the three-level
family with support `{-1,0,1}` and computes 31 exact rational unit ideals:
6 charts in degree four through moment six, 10 charts in degree five through
moment eight, and 15 charts in degree six through moment nine.  A
prime-endpoint theorem now closes that family in every degree: at odd
prime \(p\), the orders \(p\) and \(2p\) isolate the \(C^p\) and \(D^p\)
endpoints according to the two possible \(U\)-adic order inequalities.
The bounded charts remain finite-cutoff regressions.  A companion exact
arithmetic check verifies the prime coefficient and factorial congruences;
an independent pure-Python audit reconstructs both polynomial endpoint
identities and both normalized factorial cases.  The unit-star regression
then checks the primitive invariants, the \(p,2p,3p\) endpoints, and all
three normalized-order cases for the smallest star; the theorem covers
every support `{0,1,-d_1,...,-d_q}` and its reflection.  Another regression derives
the finite radial-moment recurrence, constructs the four-dimensional
resolvent differential system for a centered degree-\((2,3)\) pair, and
checks it against the factorial series.  The same target now verifies the
all-degree first-cycle theorem for support `{-2,-1,1,2}`: it enumerates the
toric invariant moments, checks
`CT(P^(kp)) = CT(P^k)^p (mod p)`, and eliminates every unique, adjacent-tie,
four-way-tie, and boundary valuation face using invariant moments through
degree twelve.  The target also recomputes three
unit Groebner bases for the direct Long-style
collapse, checks the Dvorsky--Long five-variable GVC and five-pair SIC
identities by dependency-free exact sparse arithmetic, and writes the
dimension/rank/index/degree scoreboard.  It does not
use bounded support enumeration to settle GMC(2): the accompanying
lower-face theorem handles arbitrary rotational support.  Its audit checks
supporting-line minima, Frobenius constant-term dilation, and normalized
factorial isolation on representative stars, mixed semigroups, and cycles.
The same target now verifies the normalized rank-one three-real ansatz:
the first three moments cut out Long's family scheme-theoretically, a formal
square identity proves all-order vanishing, and deletion saturation proves
five-term and degree-four minimality inside that ansatz.  Global minimality
outside the ansatz remains open.

The credited factorially weighted multitorus regression is:

```bash
python3 scripts/audit_factorially_weighted_multitorus.py
```

It checks a rank-two exposed coefficient, the normalized congruence at two
prime dilations, strict torus separation with a mixed cutoff, and the
circular Gaussian embedding.  These finite exact identities are regressions,
not the all-order proof.  Long's theorem, the local proof audit, the
prime-separating arbitrary-torus synthesis, and the one-radial search sieve
are in
[`FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md`](extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md).
The checked-in Lean development still formalizes only the rank-one Gaussian
specialization.

The companion one-profile Hopf classification is checked directly by

```bash
.venv/bin/python scripts/verify_hopf_lift_classification.py
.venv/bin/python scripts/verify_hopf_lift_classification.py --require-singular
```

This exact regression expands the phase-integrated polynomial for
several windings and endpoint multiplicities, checks the full lower-jet
binomial ladder through order twenty, and includes a non-power endpoint
profile.  It also checks the quadratic profile
`R(z)=(1-z)(z-1/5)`, whose first adjacent detector vanishes, showing that
the mixed-moment nonvanishing hypothesis is essential.  The all-order result
and exact polynomiality criterion are proved in
`extended-geometry/HOPF_LIFT_CLASSIFICATION.md`; the command performs no
search in a general `V_d`.  Finally it verifies the three triangular pure
jets in the complete class
`p=x^(-1)(C(x)+D(x)t^2)`, `deg(C)<=1`, `deg(D)<=3`; they force
`D=-C^3`, after normalization, and hence reduce the all-order statement to
the endpoint theorem.  In the next complete rectangle
`deg(C)<=2`, `deg(D)<=4`, it verifies the four triangular coefficient
solutions, the fifth-jet exceptional branch `e=-a^2/12`, and the nonzero
sixth residual `1280*a^6/81081` that removes that branch.
For `deg(C)<=3`, `deg(D)<=5`, it verifies the five triangular solutions,
the residual jets `P_6,P_7,P_8`, their two exact quadratic resultants in
the final parameter, and the Euclidean certificate `gcd(Q_5,Q_6)=1`.
This forces the two remaining weighted parameters to vanish and proves
eight-jet uniqueness in that rectangle.
The same checker verifies through order twenty the uniform triangular
coefficient
`2^(m-1)*m!/(2m+1)!!` of the new `b_m` term in the `m`-th pure jet.
The written proof uses this to reconstruct `D` successively from `C` in
every fixed numerator degree and identifies the surviving eventual tangent
directions that obstruct a first-order uniform proof.
For `deg(C)<=4`, `deg(D)<=6`, it derives the first six triangular solutions
and uses exact SymPy arithmetic over `QQ` to verify both containments in the
displayed residual Groebner ideal.  That ideal is supported only at Long's
point but has the predicted two-dimensional tangent space.
For `deg(C)<=5`, `deg(D)<=7`, it derives seven triangular solutions and five
residual jets, verifies zero-dimensionality, performs exact FGLM conversion,
and checks both containments in the displayed lexicographic ideal.  The
support is again Long's point and the tangent dimension is three.
For `deg(C)<=6`, `deg(D)<=8`, the default command reconstructs all eight
coefficients of `D`, constructs the six exact residuals through jet
fourteen, checks their term counts `[19,28,37,51,64,83]`, and verifies
residual Jacobian rank one.  The second command additionally requires
Singular and computes the exact rational modular standard basis and its
FGLM lexicographic conversion.  It checks quotient length 32 and the
triangular support relations `E^8`, `F^5`, then nonzero pure powers of
`G`, `H`, and `I` after successive substitution.  Singular's modular
routine does not give a deterministic ideal-equality guarantee for this
nonhomogeneous input, so this remains computational evidence for
fourteen-jet uniqueness, not a proof.  The same command separately performs
two deterministic rational checks: the boundary
`84*E+54*F+5=0` gives the unit ideal, and exact lift matrices put
`E^6`, `48*F^3+7*E^5`, and
`972*G^2+864*F^2*E+29*E^5-108*E^4` in the residual specialization
`H=I=0`.  The remaining compact certificate is `H^3,I^3` in the full
residual ideal; equivalently, one may eliminate `I` on the certified
principal open and exclude the resulting four-variable `H!=0` chart.  The
command also computes the exact original `H=0` slice: its quotient has
length 17 and an eight-element lex basis successively forcing
`E=F=G=I=0`.  Therefore the sole possible extra locus is the principal
chart `H!=0`, equivalently the unresolved saturation `J_6:H^infinity`.
Finite-field unit lifts on that chart must first be normalized modulo the
syzygy module before CRT: the raw coefficients are prime-dependent.  At
primes 32003 and 32009 the normalized seven-multiplier support is identical,
with term counts `[778,834,814,732,646,487,692]`; ten normalized small
primes give a 150-bit modulus but do not yet suffice for balanced rational
reconstruction.  This longer reconstruction is not part of the command
above and should use a resumable checkpoint before being promoted to a
reproducer.

A final group of checks uses the first collision-coordinate values `0,1,-1`
to fix the multiplier, expands the homogeneous 42-variable quartic, descends
the contraction to `SIC(20)`, independently reconstructs the 628-term
40-variable Laplacian witness, and verifies an all-order inverse recurrence:

```bash
.venv/bin/python scripts/generate_image_vanishing_counterexamples.py
.venv/bin/python scripts/generate_identity_slice_counterexamples.py
python3 scripts/audit_identity_slice_counterexamples_independent.py
.venv/bin/python scripts/verify_inverse_coordinate_recurrence.py
```

The provenance-preserving compression audit is separate from the absolute
two-pair construction:

```bash
python3 scripts/audit_bcw_21_low_degree_invariants.py
python3 scripts/audit_bcw_21_sextic_defect_sectors.py
Singular -q scripts/audit_bcw_21_vertical_ideal.sing
python3 scripts/audit_bcw_21_septic_component_screen.py
.venv/bin/python scripts/audit_keller_near_invariant_backtrace.py
python3 scripts/audit_keller_observable_quotients.py
python3 scripts/audit_keller_provenance_compression.py
```

The first command gives characteristic-zero rank and Lie-image certificates
through degree five for the stored 21-variable map. It also records the near-invariant
`Q=X_18*X_20-X_6*X_8` and its one-term pullback defect. The second uses two
exact torus gradings to exclude both sextic correction channels and classify
all 220 sextic Lie sectors. Unique-row peeling completely certifies 25 of
the 28 dense sectors and reduces the remaining three to small exact cores.
The full sextic Lie kernel is generated by `X_20^6`, `X_20^4*Q`,
`X_20^2*Q^2`, and `Q^3`; their pullback defects extend the fixed-space
identity to `Q[X_20]_{<=6}`. It also verifies that the two degree-seven
correction sectors are exactly `X_20` times the already excluded sextic
sectors. For the remaining pure septic problem, reduction modulo `X_20`
has 657800 columns in 204 sectors; exact unique-row peeling removes 451891
and leaves 205909 columns in 79 sectors. This last statement is a support
reduction, not a classification of the residual kernel or its lifts. The
checker also verifies that the reduced derivation is constant vertical over
the fourteen-variable base and that `X_9^7` has a nonzero first lifting
obstruction because its required sextic correction sector is empty.
The Singular command proves that the six vertical coefficients generate a
height-two ideal with five displayed minimal primes, verifies their
intersection as the radical, computes
`dim_Q(B/(A_14,...,A_19))_8=158412`, and independently confirms that
`X_0^2*X_9^6` survives in the obstruction quotient. It requires Singular
with `primdec.lib`.
The component-screen command evaluates all 77520 degree-seven base
monomials on the five minimal components. It proves that 71588 support-one
base septics have nonzero first lifting obstruction, leaving 5932 monomials
for sectorwise cancellation and higher-order analysis. Exactly eight
monomials, `X_3^(7-j)*X_5^j`, have zero first obstruction.
Stacking all five restrictions in the 29 bidegree sectors gives modular rank
61060 on the full 77520-dimensional base-septic space. Hence over `Q` the
radical-level survivor space has dimension at most 16460. Embedded torsion
and higher `X_20`-adic lifts remain open.
The backtrace command reconstructs the frozen
17-step circuit and identifies `Q=c_4*s-v_3*v_5` as a determinantal
shared-factor gate residual whose stable-source restriction is `x^2*y*z`.
The observable command proves that any rational semiconjugate quotient carrying either
`X_0` or the restricted
quadratic observable has dimension at least 13; its longer rank plateau is
printed as experiment only, while a stacked rank-20 certificate excludes a
common constant translation direction behind that plateau. The final provenance command verifies that the normalized
three-variable canonical contraction fails at its first pure moment, checks
the full twenty-coordinate inverse-recurrence dependency closure after the
known identity slice, and reports the finite stored-circuit census.  Only
the degree-at-most-six invariant statement is a nonlinear quotient-class
obstruction; the circuit census is computation, not minimality.

The independent small-witness audit can also be run directly:

```bash
python3 scripts/audit_dvorsky_gvc5_counterexample.py
```

It verifies the two pre-\(\partial_t\) identities and the resulting GVC(5)
and SIC(5) failures through order eight.  The all-order binomial proof and
the separation from ordinary-Laplacian GVC are documented in
[`DVORSKY_GVC5_COUNTEREXAMPLE.md`](extended-geometry/DVORSKY_GVC5_COUNTEREXAMPLE.md).

The smaller three-pair Image-Mathieu witness is checked by

```bash
python3 scripts/verify_three_pair_image_mathieu_counterexample.py
```

For
\(f=\tau(t-y)(wz+vt)\) and \(g=y\), the dependency-free
checker verifies exact sparse contractions
\(\mathcal E(f^m)=0\) and
\([t]\mathcal E(gf^m)=(-1)^{m-1}(m+1)!m!\) through order ten, records the
four-term bidegree-\((2,2)\) artifact, and replays the two binomial identities used in
the all-order proof.  The same script independently reads the dehomogenized
seed as the four-term cubic Gaussian polynomial
\[
P=(1-Z_2)(W_1Z_1+W_2),\qquad Q=Z_2,
\]
and verifies by exact Wick contraction through order ten that
\(\mathbb E(P^m)=0\) and
\(\mathbb E(QP^m)=(-1)^{m-1}m!\).  The all-order proof and the comparison
with Long's displayed six-term four-real cubic are in
[`THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](extended-geometry/THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).

The sharp two-pair witness is checked independently by

```bash
python3 scripts/verify_two_pair_image_mathieu_counterexample.py
python3 scripts/audit_two_pair_image_mathieu_coefficient_extraction.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.SIC2C4FiniteSum
```

For
\[
\begin{gathered}
R=\xi _1z_1+\xi _2z_2,\quad Z=\xi _1z_2,\quad
W=2\xi _2z_1,\quad T=\xi _1z_1-\xi _2z_2,\\
F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right),\qquad Q=Z,
\end{gathered}
\]
the dependency-free checker verifies \(T^2=R^2-2ZW\), the sixteen-term
bidegree-\((4,4)\) expansion, the full-rank coefficient-matrix determinant
\(48\), and exact sparse contractions through order eight.  It also
replays the phase-extracted finite sums through order \(99\).  The written
Hopf-coordinate and beta-integral argument proves for every \(m\geq1\)
\[
\mathcal E_2(F^m)=0,\qquad
\mathcal E_2(QF^m)=\frac{(4m+2)!\,m!}{(2m+1)!!}.
\]
The second command audits a non-Gaussian all-order proof.  A formal
constant-term formula for contraction gives the two sums directly; the
pure sum is an \(m\)-th finite difference of a degree-\((m-1)\) polynomial,
while polynomial division reduces the mixed sum to
\[
B_m=\sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+1},\qquad
(2m+1)B_m=2mB_{m-1}.
\]
Thus \(B_m=2^m m!/(2m+1)!!\).  The audit separately checks the chart
identity, both chart constant-term expansions, the divisibility certificate,
finite differences, the general denominator-remainder invariance, and the
recurrence with exact rational arithmetic.  Its cutoff is a regression; the
displayed degree bounds and termwise identities in the written proof are
all-order.  The Lean command
formalizes the general finite-difference cancellation, the
denominator-remainder theorem for alternating quotient sums and its rank-one
endpoint-residue specialization, the specialized normalized products, and
the generalized repeated-pole beta recurrence, the finite
remainder-to-jet identity, and the order-one
factorial/double-factorial evaluation; it also formalizes the scalar chart
identity, the coefficient functional and algebraic beta identity, the
monomial and balanced-array contraction/coefficient-extraction equalities,
the selected chart coefficients, the formal integrals of the resulting
chart polynomials, the product-polynomial evaluation at natural numbers,
and both normalized displayed binomial-sum identities.  It now also
represents (4.3) literally in \(\mathbb Q[v][x,x^{-1}]\), proves the
all-order binomial expansion of its powers, identifies the pure and mixed
constant terms, and proves their final formal-integral values.  It also
defines the original four-variable \(F,Q\) as `MvPolynomial` objects and
proves that their displayed substitution gives the Laurent witness and its
mixed powers.  The remaining Lean integration seam is only a wrapper between
the generic balanced coefficient array in the contraction theorem and
`MvPolynomial.coeff`; the algebraic chart proof itself is formalized on both
sides of that representation boundary.

The positive-characteristic phase diagram is replayed separately by

```bash
.venv/bin/python scripts/verify_two_pair_sic_characteristic_p.py
```

It clears denominators with \(\widetilde F=2F\), checks that the quadric
chart is nondegenerate for every odd prime and that the coefficient tensor
has full rank away from \(2,3\), and verifies
\[
\mathcal E_{2,p}(\widetilde F^m)=0,\qquad
\mathcal E_{2,p}(Z\widetilde F^m)
=\overline{2^m(4m+2)!m!/(2m+1)!!}.
\]
Legendre floor sums, base-\(p\) digit sums, and the Lucas--Kummer carry
criterion all give the exact nonvanishing condition \(4m+2<p\).  The
same audit treats every \(R^k(2F)^r\) of degree \(d=4r+k\) through
degree twenty and verifies the uniform criterion \(dm+2<p\).
At prime-power level it checks the integral radial quotient
\[
A_{s+1}/A_s=16(4s+3)(4s+5)(s+1)^2
\]
and the resulting valuation monotonicity.  It also verifies the
non-radial re-entry for \(R(2F)\): the order-four moment is zero modulo
\(11^2\), but the order-five moment is \(22\) modulo \(11^2\).
The general signed consecutive-order valuation recurrence (4.14e) is
audited through degree sixteen, prime \(31\), and order forty.
The checker also computes the exact coefficient determinants and every
exceptional modular rank for \(R^k(2F)\), \(0\leq k\leq4\), through
degree eight; these are tabulated in (7.6) of the written proof.  The
four binomial-convolution diagonal symbols and their lower-Hessenberg
determinant recurrence are independently compared with expanded
polynomials for every \(0\leq k\leq20\).  The closed characteristic-two
rank formula
\[
\operatorname {rank}C_{R^k(2F)}
=2^{1+s_2(\lfloor(k+2)/2\rfloor)}
\]
is checked through \(k=128\).
For the non-power profiles it checks the universal necessary cutoff
\(4hm+2<p\), proves the closed height-two formula
\[
C_{2,m}=\frac{4^m m!}{\prod_{j=0}^m(4j+1)},
\]
and records the first higher-profile numerator-prime holes at
\((h,m,p)=(6,1,47)\) and \((4,5,89)\).
The checker also verifies the characteristic-two one-sided degeneration, the
characteristic-three Hilbert--Mumford unit ideal, and the naive Hasse
formulas
\(\mathcal H_2(\widetilde F^m)=16^m\) and
\(\mathcal H_2(Z\widetilde F^m)=2m16^m\).
It additionally checks the binomial intertwining and Lucas no-carry units
behind the complete Hasse Image-kernel theorem and its reduction to the
\(p\)-typical operator orders \(1,p,p^2,\ldots\).
For that modified Image it verifies the one-pair counterexample
\[
f=\xi z^p,\qquad g=z,
\]
whose pure moments vanish because \(\binom{pm}{m}=0\), while the mixed
moments are nonzero at every \(m=(p^e-1)/(p-1)\).
The written proof in
[`TWO_PAIR_SIC_CHARACTERISTIC_P.md`](extended-geometry/TWO_PAIR_SIC_CHARACTERISTIC_P.md)
shows division-freely that the Image-kernel identity survives in every
characteristic.  Frobenius gives
\(\mathcal E_{r,p}(f^p)=f(0,z)^p\), so the single \(p\)-th pure moment
forces the dual-degree-zero part to vanish; every fixed mixed contraction
then vanishes for all \(m\geq p\).  This proves ordinary
\(\operatorname{SIC}(r)\) for every \(r\) and \(p>0\), with sharp cutoff
witness \(f=\xi _1,\ g=z_1^{p-1}\).  The finite replay is not being used
as an all-order or periodicity argument.

The Frobenius/\(p\)-curvature bridge is tested by

```bash
.venv/bin/python scripts/research_two_pair_sic_frobenius_curvature.py
```

For
\[
M_{d,r}(m)=4^{rm}(dm+2)!((rm)!)^2/(2rm+1)!,
\]
it derives the coprime minimal order-one recurrence
\(A_{d,r}(m)M(m+1)=B_{d,r}(m)M(m)\) at nine radial rows.
At every good prime its recurrence-operator \(p\)-curvature is proved and
directly replayed as
\[
\prod_{i=0}^{p-1}\frac{B_{d,r}(m+i)}{A_{d,r}(m+i)}
=d^d(m^p-m)^d.
\]
Separately, the normalized angular beta period has Picard--Fuchs operator
\(\theta(2\theta+1)-x(\theta+1)^2\).  Its differential \(p\)-curvature is
computed at every odd prime through \(101\); it is always nonzero,
square-zero of rank one, with poles only at \(0,2\).  This bounded
differential calculation is not an all-prime proof.  The correlation audit
shows why neither curvature recovers the exact phase diagram: first radial
lifts have one common Picard--Fuchs operator, while the recurrence shift
norm cancels the separate zero/pole factors responsible for prime-power
re-entry.  The exact reusable mechanism is instead the local rule
\[
v_p(M(m+1))-v_p(M(m))=v_p(B_{d,r}(m))-v_p(A_{d,r}(m)).
\]
The status and the degree-eight same-curvature/different-phase control are
in
[`TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`](extended-geometry/TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md).
The resulting integral-lattice postprocessing stage is incorporated in
[`HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md`](extended-geometry/HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md)
and in the next-step protocols for the bidegree-\((3,3)\) and rank-two
bidegree-\((4,4)\) recurrence programmes.

The exact local geometry of this displayed \(F\) is checked by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_local_moduli.py
```

The all-order fourth-order continuation is checked separately by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_fourth_order.py
```

It derives the combined shifted beta-tail for numerator degrees
\(4,8,12,16\), rather than imposing an independent cutoff on each degree.
The exact \(3220\)-by-\(455\) universal polynomial-section system has rank
jump \(90\) to \(91\).  Restoring the eleven free cubic-lift parameters at
the reduced direction \((1,2,3,4,5)\) leaves one affine-linear equation
and one rank-one quadric, of dimension nine and degree two.  Its
discriminant square class is \(41\), so the fiber is two conjugate affine
\(9\)-planes over \(\mathbb Q(\sqrt {41})\) and has no rational point.
This certifies a nonradial geometric fourth-order lift, not a formal arc.

One explicit conjugate pair of fourth-order lifts is continued by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_fifth_order.py
```

The checker reconstructs the complete
\(\mathbb Q(\sqrt {41})\)-valued jet, verifies its four coefficients
against the all-order fourth tail, and derives the combined
degree-\(4,8,12,16,20\) fifth tail of rank \(56\).  After all thirteen
fourth-tangent corrections, the fifth rank jumps \(2\) to \(3\).  A
primitive obstruction supported only on residual rows \(12,13,14\) has
coefficients
\((2727113757934325760,-407042494824,17047)\).
This obstructs the selected conjugate jets, not the entire
nine-dimensional fourth-lift components.

The algebraization samples and their component-wide fifth obstructions are
computed by

```bash
.venv/bin/python scripts/research_two_pair_counterexample_algebraization.py
.venv/bin/python scripts/research_two_pair_counterexample_fifth_component.py
.venv/bin/python scripts/analyze_two_pair_counterexample_fifth_factor.py
```

For the generic rational direction \((2,-1,3,1,-2)\), the documented
direction \((1,2,3,4,5)\), and the pure apolar-odd direction
\((0,1,0,0,0)\), the exact fourth fiber again has dimension nine, degree
two, and discriminant square class \(41\).  At one point on each conjugate
pair, restoring the previously omitted eleven-dimensional cubic-tangent
kernel changes the fifth coefficient/augmented ranks from \(2/3\) to
\(4/5\), so the selected points remain obstructed.  The component-wide
commands then parameterize all nine coordinates on one component and
restore the full eleven-dimensional kernel.  After the rank-two
fourth-tangent image is eliminated, the remaining coefficient and
augmented ranks over the component function field are \(2/3\).  Every
coefficient \(3\)-by-\(3\) minor vanishes identically, while an augmented
\(3\)-by-\(3\) minor is a nonzero constant with nonzero quadratic norm.
Thus both conjugate affine \(9\)-planes are uniformly obstructed at fifth
order for all three directions.  The exact \(F_{1+s,1}\) control is
polynomial of parameter degree three, so its coefficients at orders
\(4,\ldots,12\) vanish.

The final command is a fast exact replay from nine stored component samples
on \(h_3=h_4=0\).  Within a quadratic projective ansatz, eight samples
reconstruct a selected constant augmented minor, up to a rational chart
unit, as
\(h_0A_1+\sqrt{41}B_2\), with \(A_1\) linear and \(B_2\) a nonsingular
quadratic form.  The ninth sample checks the reconstruction exactly.  Its
vanishing locus is a smooth conic over \(\mathbb Q(\sqrt{41})\); this
candidate exceptional locus has no rational projective point, by exact binary
discriminants.  This finite reconstruction does not prove the ansatz or a
universal identity, and it does not yet prove that the other augmented minors
have no common zero on the conic over \(\mathbb Q(\sqrt{41})\).
See
[`TWO_PAIR_COUNTEREXAMPLE_ALGEBRAIZATION_RESEARCH.md`](extended-geometry/TWO_PAIR_COUNTEREXAMPLE_ALGEBRAIZATION_RESEARCH.md).

The minimum-degree separating invariant and the low-degree invariant-ring
calculation are checked independently by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_missing_invariant.py
```

The degree-four moment-field continuation is replayed by

```bash
.venv/bin/python scripts/verify_degree_four_tau_even_parameters.py
.venv/bin/python scripts/research_degree_four_moment_field.py \
  --max-weight 16 --targets odd-square
.venv/bin/python scripts/verify_degree_four_diagonal_moment_field.py
.venv/bin/python scripts/verify_degree_four_single_phase_moment_fields.py
.venv/bin/python scripts/research_degree_four_phase_one_chart.py \
  --prime 101 --orders 1 2 3 4 5 6 7 8 9 10 \
  --threads 6 --timeout 300 --groebner-basis 2 \
  --compare-even-parameters 2 3 5 7 11 69 3 6 \
  --test-apolar-orbit 2 3 5 7 11 69 3 6 --certify-example
```

The first command constructs twenty-two algebraically independent
apolar-even trace invariants of degrees \(1,2^4,3^9,4^8\).  Their exact
modular Jacobian rank is \(22\), and their combined cotangent matrix with
\(\mu_1,\ldots,\mu_{22}\) still has rank \(22\).  The second command is a
bounded exact search: it proves that no relation
\(Q(\mu)+c_{234}^2P(\mu)=0\) of invariant weight at most sixteen exists
with that support.  It does not exclude a higher-weight denominator or a
higher even minimal polynomial.  The third command requires Singular.  On
the diagonal quartic slice it proves that the first five moments give a
finite parameter ring of quotient length \(120\), and that the complete
first-six-moment fiber through \((2,3,5,7,11)\) consists exactly of that
point and its reversal.  Finiteness then proves that the full diagonal
moment field is the reversal-fixed field, of exact generic degree two.
The reversal is also the \(\operatorname{SL}_2\) Weyl action on the
diagonal space, so these are two raw parameter points but one invariant
quotient point; this is a fixed-locus control, not a degree-two quotient
test.
The fourth command also requires Singular.  It repeats the finite-parameter
and exact-fiber calculation on all ten coordinate choices of a
\(\tau\)-even positive/negative direction pair in phases \(1,2,3,4\).
Every resulting six-dimensional parameter space has parameter quotient length
\(360\), and its first-seven-moment fiber through
\((2,3,5,7,11,221)\) is exactly the reduced reversal pair.  Openness
then gives exact degree two and fixed-field equality for a nonempty
Zariski-open family of raw direction-pair parameter spaces in every
phase.  The odd cubic is nonzero on exactly four coordinate cross-pairs
in phases one and two, so only those four are genuinely apolar-moving
quotient tests; the other six lie in the fixed locus.
The fifth command is an exact \(F_{101}\) experiment using `msolve`.  On
the eight-dimensional chart containing both positive and both negative
phase-one directions, the first ten moments have a reduced four-point
fiber.  All moments through order eleven and all twenty-two known even
parameters agree on the extra branch, while \(c_{234}\) changes from
\(11\) to \(-11\).  A four-equation orbit basis proves that this branch
is \(\operatorname{SL}_2\)-conjugate to the apolar reversal.  Thus the
four raw points form two candidate quotient points.  The same checker
verifies the exact rational branch \(u=5/3,w=6\): its first eleven moments
agree over \(\mathbb Q\), its odd cubic is \(-1728\) versus \(1728\), and
the matrix with rows \((0,-1/\sqrt3)\), \((\sqrt3,0)\) conjugates
\(\tau(p)\) to it.  It also reconstructs the four rational points
\(p,q,\tau(p),\tau(q)\) and verifies that the first-eight-moment Jacobian
is nonzero at each, so all four are reduced and isolated in
characteristic zero.  Fiber completeness remains proved only modulo
\(101\): additional characteristic-zero components have not been
excluded, and no characteristic-zero generic-degree conclusion is
claimed.
The full \(22\)-dimensional degree-four moment-field equality remains
open; see
[`DEGREE_FOUR_MOMENT_FIELD.md`](extended-geometry/DEGREE_FOUR_MOMENT_FIELD.md).

The completed-coordinate comparison in degrees three through five is
replayed by

```bash
.venv/bin/python scripts/research_completed_moment_algebra.py \
  --degrees 3 4 5 --max-weight 10
```

This exploratory checker constructs the quadratic Casimir decompositions,
verifies the coefficient rows
\(\binom{2d+1}{d-r}\) in \(\mu_2\), and gives exact modular moment-Jacobian
ranks \(13,22,33\).  It runs linear-denominator nonrelation searches for
the missing quadratics over the moments and over the moments with \(q_2\).
For the square of a first apolar-odd invariant it additionally tests the
proposed \((q_2,q_6)\) base and the full quadratic completion.  It also
records bounded Hilbert-series necessary tests for natural and minimally
corrected parameter-degree sequences, and checks the propagated
all-moment-zero witnesses in degrees four and five.  A zero relation
intersection is an exact bounded nonexistence certificate; Hilbert
compatibility is not a proof of a nullcone zero fiber.  See
[`COMPLETED_MOMENT_ALGEBRA_RESEARCH.md`](extended-geometry/COMPLETED_MOMENT_ALGEBRA_RESEARCH.md).

The automatic missing-invariant and \(d=6\) extension is replayed by

```bash
.venv/bin/python scripts/research_completed_moment_algebra.py \
  --degrees 3 4 5 6 --invariant-cutoff 6 \
  --skip-relation-tests --power-witness-cutoff 12 \
  --ladder-beta-check 32 \
  --output artifacts/generated-results/automatic_missing_invariants_d3_d6.json
```

The refined weight-zero-minus-weight-two calculation splits the invariant
spaces by the apolar involution through polynomial degree six and subtracts
the moment-monomial subspace.  It proves that the first missing degree is
two, with even multiplicity \(d-1\), and enumerates the first odd cubic
triples in degrees four through six.  It also certifies full modular
moment-Jacobian ranks \(13,22,33,46\), runs the Hilbert necessary tests for
candidate augmented parameter systems, and evaluates \(q_2\) on the
propagated all-order witnesses through \(d=6\).  The conclusion that \(q_2\)
removes the recorded witnesses is not a classification of every semistable
moment-zero component.  The same run verifies the radial Casimir recurrence
\[
q_{2r}^{(d+1)}(Rf)=(d-r+1)(d+r+2)q_{2r}^{(d)}(f),
\]
the resulting all-degree closed formulas for \(q_2,q_4\) on
\(R^{d-4}F_4\), and the exact power-witness pattern through \(F_4^{12}\).
A finite-difference proof using the chart expansion of \(F_4^m\) shows
for every \(m\geq1\) that all earlier quadratics vanish and
\(q_{2\lceil m/2\rceil}(F_4^m)\ne0\).  The run regresses the exact beta
sums used in that proof through \(m=32\).  The stronger formula listing
every surviving torus phase is recorded through \(m=12\) and remains
bounded evidence.

The completed-invariant comparison specialized to bidegree \((3,3)\) is
replayed by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_casimir_fiber.py
```

The default run does not launch the older direct boundary standard-basis
calculations.  It compares
\((\mu_1,\ldots,\mu_{12},\mu_{14})\),
\((\mu_1,\ldots,\mu_{12},q_2)\), and full-rank mixed
moment/Casimir systems of the same total invariant degree \(92\).  It
also evaluates the complete weight-\(14\) monomial spaces generated by
the lower moments with \(q_2\), and with \(q_2,q_4\).  The modular ranks
prove over characteristic zero that \(\mu_{14}\) is independent from
the pure Casimir span modulo the lower-moment span.  Hilbert compatibility
does not prove that any displayed zero fiber is the nullcone.

The new null-quadratic synchronization experiment is included explicitly by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_casimir_fiber.py \
  --run-normal-rank-locus \
  --run-residual-normal-probe \
  --run-complete-normal-fiber \
  --run-exceptional-normal-fibers \
  --timeout 300 \
  --power-bound 50
```

On the normalized chart \(F_2=X^2\), the synchronization ideal has the
seven normal coordinates \(s_3,s_4,s_5,s_6,t_2,t_3,t_4\).  The checker
computes the exact linear normal-symbol matrix of
\(\mu_2,\mu_3,\mu_4\), proves its generic rank is three, and proves that
the common divisor of its nonzero maximal minors is the displayed
irreducible cubic \(P\).  It then records two finite-field calculations
at \(p=32003\): the residual rank-drop locus after dividing by \(P\),
and the complete normal fiber of \(\mu_2,\ldots,\mu_{12}\) above the
fixed allowed-coordinate point \((20,27,36,47,60)\).  The latter has
dimension zero and quotient length \(195\).  Good reduction therefore
proves characteristic-zero transverse isolation at this point and on a
nonempty open subset of the synchronized chart.  The recorded coordinate
powers are finite-field certificates only.  The same run decomposes the
residual rank support exactly into two quadratic-field components and one
lower rational locus, proves that all three are disjoint from \(P=0\),
and tests full normal fibers at good reductions of exact algebraic points
on all four exceptional strata.  The three top-dimensional exceptional
fibers have quotient length \(195\); the lower locus has length \(197\).
Consequently transverse isolation holds on a nonempty open subset of
every linear-rank stratum.  Prefix computations through orders nine and
ten timed out and make no minimal-cutoff claim; proper closed subsets
inside the exceptional strata, the \(F_2=0\) chart, and global nullcone
equality remain open.

The exact diagonal fixed-field theorem in all three degrees is replayed
by

```bash
.venv/bin/python scripts/verify_completed_moment_diagonal_fields.py
```

For \(d=3,4,5\), the checker proves that the first \(d+1\) diagonal
moments form a parameter system of quotient length
\((d+1)!=24,120,720\).  Homogeneous finite-field standard bases at
\(32003\), followed by projective properness and the regular-sequence
Hilbert series, give the characteristic-zero finiteness statement.
An invertible midpoint/direction change then proves over \(\mathbb Q\)
that the first \(d+2\) moment fiber through the selected integral point is
exactly
\((y_0,\ldots,y_{d-1},s^2-1)\).  Hence the full diagonal moment field is
the reversal-fixed field and has exact generic degree two in each degree.
These are slice theorems, not statements about the full invariant
quotients.

The exact single-phase extension in degrees three and five is replayed
by

```bash
.venv/bin/python scripts/verify_completed_moment_single_phase_fields.py
```

For one matching apolar-eigendirection pair in every nonzero phase, the
checker proves that the first \(d+2\) moments have full Jacobian rank and
that adding \(\mu_{d+3}\) makes the moment-origin fiber finite.  Exact
standard bases over \(\mathbb F_{32003}\), weighted-projective
properness, and Nakayama lift the displayed two-point reversal fiber to
characteristic zero.  The quotient lengths at the special moment origin
are \(54\) for \(d=3\) and \(1934\) for \(d=5\).  Quintic cross-direction
slices in phases one and two have
\(c_{234}=-273686400/7\), so those reversal pairs are genuinely distinct
invariant-quotient points.  The remaining slices are raw parameter-space
fixed-field controls.  None of these slice certificates determines the
generic degree on the full invariant quotient.

The first branchwise global-\(d=4\) \(q_2\)-augmented nullcone attack is
replayed by

```bash
.venv/bin/python scripts/research_degree_four_q2_augmented_nullcone.py \
  --prime 32003 --max-jet 4 --composition native \
  --ordering dp --timeout 300
```

On the branch \(q_2=0,F_2\ne0\), normalize \(F_2\) to a highest-weight
square.  At one deterministic synchronized point, the checker expands
the first twenty-one moments in the twelve forbidden weight coordinates.
Moments two through five provide four formal pivots.  In the remaining
eight variables, the exact quadratic and cubic jet ideals over
\(\mathbb F_{32003}\) have dimensions six and four.  This is a bounded
normal-jet frontier, not a formal-isolation or global-nullcone theorem.
The native eight-variable quartic basis reaches its declared
\(300\)-second timeout, which supplies no mathematical evidence.
The cubic support decomposition and the cheaper dominant-sheet quartic
restrictions are replayed by

```bash
.venv/bin/python \
  scripts/research_degree_four_q2_cubic_decomposition.py \
  --prime 32003 --timeout 180
```

After the forced equation \(x_7=0\), the checker proves the radical of
the \(x_6=0\) cubic support, factors its binary cubic into a rational
sheet and an irreducible quadratic sheet over \(\mathbb F_{32003}\),
and records a distinct degree-nine off-\(x_6\) saturation with an
irreducible generic degree-nine fiber.  Quartic restriction collapses
both dominant sheets to the same three-plane.  The generic quartic
restriction on the off-axis saturation still times out, so formal
isolation, quintic normal terms, and all \(F_2=0\) boundary branches
remain open.  See
[`DEGREE_FOUR_Q2_AUGMENTED_NULLCONE.md`](extended-geometry/DEGREE_FOUR_Q2_AUGMENTED_NULLCONE.md).

The explicit first \(d=2\) moment relation is reconstructed and verified
by

```bash
.venv/bin/python scripts/verify_two_pair_d2_moment_field_relation.py
```

With factorial-normalized centered moments \(x_m\), this produces the
monic quintic relation for \(x_7\) over
\(\mathbb Q[x_2,\ldots,x_6]\).  Five good primes reconstruct its \(241\)
nonzero rational coefficients.  The checker then verifies the relation
exactly at \(418\) fixed integral coefficient matrices whose full ansatz
evaluation matrix has rank \(418\) modulo \(1000003\).  The generated
artifact stores the five sparse coefficient blocks.  The same checker
also proves, by five further modular rank jumps, that the centered
moments \(x_8,\ldots,x_{12}\) are each outside the algebra generated by
their predecessors.  Thus the first seven moments generate the moment
field but not the polynomial moment algebra.  Finally, it certifies that
the affine operator slice
\[
 \begin{pmatrix}0&1&0\\a&b&c\\d&e&f\end{pmatrix}
\]
has dense \(\mathrm{SL}_2\)-saturation and expands reconstructed formulas
for \(x_{13},\ldots,x_{18}\) identically to zero in
\(\mathbb Q[a,b,c,d,e,f]\).  Hence every moment through order \(18\)
belongs to the algebra generated by the first twelve.

The missing-invariant checker constructs the five Casimir projectors for
\(\operatorname{End}(\operatorname{Sym}^4)\), verifies the five quadratic
values and the exact expression of \(\mu_2\) in that basis, evaluates the
ten primitive cubic contractions, and certifies modulo \(1000003\) that
\(\mu_1,\ldots,\mu_{22}\) have Jacobian rank \(22\).  The nonzero modular
minor is an exact characteristic-zero algebraic-independence certificate,
not a numerical rank test.  It also constructs the equivariant apolar
adjoint, verifies its alternating signs on the five irreducible summands,
and checks that the first odd cubic is nonzero.  The written coefficient
reindexing proves that every moment is adjoint-even in all orders; hence
the moment and invariant fields differ and the algebra conductor is zero.
The same weight calculation finds the first odd invariant in bidegrees
two and three, and gives degree ten for the proved first-six-moment
parameter map on \(V_2\).  A nonzero fifteen-by-fifteen modular evaluation
determinant shows that \(\mu_7\) leaves the first-six moment field; the
prime degree-five tower then proves that the complete \(d=2\) moment field
is exactly the apolar-adjoint fixed field, of generic degree two in the
full invariant field.

The local-moduli checker computes the contraction-preserving stabilizer
and orbit, uses the
all-order Hopf coefficient identity to obtain the thirteen-dimensional
all-moment tangent space, records the seven independent quadratic lifting
obstructions, proves that their quotient radical is a five-plane, and
splits the quotient moment tangent as \(2+8\) and that five-plane as
\(1+4\) under the apolar involution.  It constructs a polynomial cubic
lift for every direction on that plane and verifies the defect-preserving
family
\[
F_{a,b}=\frac{aR+bZ}{2}
\left(2W(aR+bZ)^2-2abR^3-b^2R^2Z\right).
\]
The parity split reduces the certified cubic correction system from
\(490\times195\) to a consistent \(154\times73\) equivariant system and
reduces the proposed fourth-order continuation from \(980\times455\) to
at most \(644\times273\).
The generated artifact contains the exact tangent matrix, kernel basis,
obstruction quadrics, radical equations, and cubic correction.  The family
proves a positive-dimensional local quotient but is not claimed to exhaust
every reduced local branch; the first unresolved local equations are fourth
order.

The same checker verifies the propagated forms
\(F_d=R^{d-4}F\) for \(4\leq d\leq10\) through moment four.  The written
radial-shift identity proves, for every \(d\geq4\) and \(m\geq1\),
\[
\mathcal E_2(F_d^m)=0,\qquad
\mathcal E_2(QF_d^m)=\frac{(dm+2)!\,m!}{(2m+1)!!}.
\]
It also verifies through moment three and degrees \(4\leq d\leq15\) the
bounded-radial-order family
\[
G_{r,k}=R^kF^r,\qquad d=4r+k,\quad 0\leq k\leq3,
\]
whose all-order formula is
\[
\mathcal E_2(G_{r,k}^m)=0,\qquad
\mathcal E_2(QG_{r,k}^m)
=\frac{(dm+2)!\,(rm)!}{(2rm+1)!!}.
\]
The written proof shows that \(G_{r,k}\) has exact \(R\)-adic order \(k\).
Thus the family is \(R\)-primitive in degrees divisible by four and has
radial order at most three in every degree.
Finally, the checker verifies for \(1\leq h\leq4\) and \(m\leq3\) the
explicit non-power Hopf-profile family \(\Phi_h\) of degree \(4h\), whose
all-order mixed formula is
\[
\mathcal E_2(Q\Phi_h^m)
=(4hm+2)!\int_0^1
(1-v^2)^m(1+v^2)^{(h-1)m}\,dv.
\]
The written endpoint-contact proof shows that every pure moment vanishes,
the displayed detector is positive, and \(\Phi_h\) is \(R\)-primitive and
not a proper power.

The degree-five bilinear-multiplier obstruction is checked by

```bash
.venv/bin/python scripts/verify_two_pair_degree_five_multiplier_obstruction.py
```

For \(L=aR+bZ+cW+eT\), the checker derives the first four pure moments of
\(LF\).  On the nonradial branch, moments one and two eliminate \(c,a\);
moment three gives
\[
q(u)=8019u^4-623736u^2+3219760,
\]
and moment four gives
\[
p(u)=136323u^6-5359284u^4-174020976u^2-802761152.
\]
Their exact rational gcd is one.  Hence the first four moments force
\(L=aR\), excluding every nonradial bilinear lift of the quartic seed in
degree five.  This does not classify the full \(V_5\).

The sharp primitive finite-prefix family is replayed by

```bash
python3 scripts/verify_two_pair_primitive_prefix_obstruction.py
```

For \(G_{d,\lambda}=R^{d-4}F+\lambda Z^d\), the checker directly verifies
in degrees \(4\leq d\leq8\) that moments one through \(d\) vanish and
\[
\mathcal E_2(G_{d,\lambda}^{d+1})
=(d+1)\lambda(d(d+1)+1)!\frac{d!}{(2d+1)!!}.
\]
The written phase-support proof establishes the formula for every
\(d\geq4\).  For \(\lambda\ne0\), these are \(R\)-primitive sharp-prefix
points, not all-order counterexamples.
The same checker tests the stronger triangular statement for
\(4\leq d\leq7\): if
\(H=\sum c_jR^{d-j}Z^j\) has least nonzero phase \(s\), then
\(R^{d-4}F+H\) is detected exactly at moment \(s+1\).  The written proof is
all-order and excludes every nonzero positive-phase triangular correction.
It also checks the two-sided degree-five extreme ansatz
\(RF+aZ^5+bW^5\): moment two is \(921600ab\), and the remaining two
branches are excluded by explicit nonzero moments four and six.
Finally, for \(4\leq d\leq7\), it verifies the odd-height family
\[
J_{d,\lambda}=R^{d-4}F+\lambda Z^{d-1}T,
\]
whose all-order phase-parity proof gives zero moments below \(2d\) and
\[
\mathcal E_2(J_{d,\lambda}^{2d})
=\binom{2d}{2}\lambda^2(2d^2+1)!
\frac{(2d-2)!}{(4d-1)!!}.
\]
The same replay checks all opposite odd-height monomial pairs, all
opposite even-height pairs of phase \(s\geq3\), and the phase-two
exceptional branch in the displayed degree range against the written
all-degree formulas.  It also verifies the degree-five phase-one
elimination for
\(RF+aZT^4+bWT^4\): the normalized moments \(2,3,4\), their exact
lexicographic remainder, and the nonzero resultant
\(-418538718730248905250\).

The remaining uniform phase-one elimination is reproduced by

```bash
.venv/bin/python scripts/verify_two_pair_phase_one_uniform_obstruction.py
```

It derives the symbolic height-dependent moments, eliminates to a
quadratic and cubic in the remaining correction parameter, and factors
their resultant.  The only nonlinear height factor has degree \(31\)
and no root modulo \(29\).  Together with the preceding checker, this
excludes every correction
\(aZ^sT^{d-s}+bW^sT^{d-s}\) for every \(d>4\).

The first multi-pair local obstruction is checked by

```bash
.venv/bin/python scripts/verify_two_pair_degree_five_odd_height_quadratic_rigidity.py
```

For the degree-five odd-height correction
\(a_2Z^2T^3+a_4Z^4T+b_2W^2T^3+b_4W^4T\), it constructs the quadratic
coefficients of moments \(2,\ldots,10\).  Their exact Gröbner basis
contains the cube of every parameter, so their radical is the homogeneous
maximal ideal.  It also computes the Hilbert vector \((1,4,1)\), length
six, and a nondegenerate socle pairing, identifying the quotient as a
compressed quadratic Artin--Gorenstein algebra.  This is a local formal
obstruction.  The checker then treats the full ten-dimensional monomial
correction space: moments \(2,\ldots,7\) eliminate the six even-height
second-order variables, and projected moments \(8,\ldots,11\) give four
quadrics forming a length-sixteen complete intersection with Hilbert
vector \((1,4,6,4,1)\).  Thus \(RF\) is formally isolated in this full
space; this is not a global classification of finite corrections away
from \(RF\).

The higher-degree continuation is checked by

```bash
.venv/bin/python scripts/verify_two_pair_higher_degree_monomial_formal_rigidity.py
```

Using exact Hopf-angular integration, it treats the full monomial
correction spaces in degrees six and seven.  The consecutive linear
blocks have sizes six and eight, respectively.  In both degrees, the
projected odd-height obstruction consists of six quadrics forming a
complete intersection with Hilbert vector
\((1,6,15,20,15,6,1)\) and length \(64\).

The all-degree linear pivot theorem has the dependency-free regression

```bash
python3 scripts/verify_two_pair_linear_pivot.py
```

The written proof uses the endpoint factorization
`2p=-(1+x)(2+x)+(1-u)x^(-1)(1+x)^3` over
\(\mathbb Z_{(2)}\).  With
\(\delta_s=s+\nu_2(s!)\), the scaled local Smith exponents are two copies
of \(\delta_2,\delta_4,\ldots,\delta_d\) in even degree.  In degree
\(d=2h-1\), they are one \(\delta_1\), two copies of
\(\delta_3,\ldots,\delta_{2h-1}\), and one \(\delta_{2h+1}\).  The
checker independently reconstructs the terminating binomial entries and
verifies this pattern exactly in degrees 5 through 25.  This finite range
is a regression for the filtered proof, not the proof itself.

Degrees eight through eleven use the scalable good-prime replay

```bash
.venv/bin/python scripts/verify_two_pair_degree_eight_monomial_formal_rigidity.py
```

The checker constructs the exact rational pivot and projected systems,
reduces them at \(1000003\), and invokes Singular.  Degrees eight and nine
give eight-quadric quotients of dimension \(256\).  Degree ten gives a
ten-quadric quotient of dimension \(1024\), and degree eleven gives the
same ten-quadric dimension.  The corresponding ninth and eleventh
variable powers reduce to zero, proving the characteristic-zero complete
intersections.
See
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
This proves `not MN_d` for every `d>=4`, `not SIC(2)`, and, with the known
one-pair theorem, the exact minimum failing pair dimension two.  The finite
replay is not being used as the all-order proof.

The coefficient-rank frontier inside bidegree \((4,4)\) is replayed by

```bash
python3 scripts/verify_two_pair_sic_bidegree44_rank_frontier.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_invariants.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_all_order_audit.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_direct_chart.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_channel.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_boundaries.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_off_diagonal.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_off_diagonal_boundaries.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_rank_two_single_shear.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_rank_two_double_shear_real_prefix.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree44_rank_two_swap_slice.py
.venv/bin/python scripts/verify_two_variable_quartic_squarefree_pivot.py
.venv/bin/python scripts/verify_two_variable_quartic_two_root_finite.py
```

The dependency-free checker verifies the determinantal dimensions
\(r(10-r)\) for \(1\leq r\leq4\), exact representatives of ranks one
through four, and the coefficient formulas for pure contractions and all
four bilinear mixed multipliers through order four. It also checks the
nilpotent endomorphism trace screen on a fixed-flag one-sided chart and
replays the known rank-five determinant and mixed formula. The written
split-symbol argument excludes rank one for arbitrary SIC multipliers,
giving the rigorous interval \(2\leq r_{\min}\leq5\). Ranks two, three,
and four remain open; finite rank-two residuals are not treated as an
exact counterexample or a lower bound. See
[`TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md).

The swap-slice checker treats
\(F_P=\xi_1^4P(z_1,z_2)-\xi_2^4P(z_2,z_1)\).  It proves that odd moments
vanish by involution and that the even moments of orders
\(2,4,6,8,10\) have a length-twelve zero scheme whose radical consists
of \(P=(z_2\pm z_1)^4\).  Seven good primes reconstruct the triangular
lex ideal; exact rational reductions, the matching special-fiber length,
and four unit charts on the complete \(a_4=0\) projective boundary certify
the characteristic-zero result.  Both reduced points have coefficient
rank one.  The parity-even three-parameter core is already excluded by
orders \(2,4,6\).  This is an exact structured-slice theorem, not a global
rank-two exclusion.

The direct quotient checker fixes \(U=[I_2;A]\), verifies the pure and
mixed relative-period identities, and certifies the one-sided
function-field recurrence and mixed cutoff.  The two-row checker then
computes the exact degree-\(604\) seven-moment scheme and the
eighth-moment unit ideal on the dense support.  Its boundary companion
classifies all 135 proper support orbits: every mixed orbit has a sharp
cutoff at most seven, including separate rank-one and exact-rank-two
localizations on the three mixed-rank tori.
The off-diagonal two-row checker generalizes the dense calculation to all
ten row pairs.  Six simultaneous-reversal representatives suffice; after
fixing \(U=(e_r,e_s)\), every coefficient torus has exact rank two and its
exact \(\mathbb Q\)-moment ideal is a unit through \(\mu_8\).  Its boundary
companion classifies all 1174 proper support orbits of the ten charts,
separates the zero, fixed-flag, rank-one, mixed-rank, and exact-rank-two
strata, and checks 942 exact rank-two opens through \(\mu_{10}\).  It also
records exact degree-8 and degree-252 delayed fibres and their respective
\(\mu_9\)- and \(\mu_{10}\)-unit certificates.  Together the two checkers
prove all ten complete off-diagonal two-row coordinate subspaces SIC-safe.
The single-shear checker then moves outside the coordinate row planes.  On
all 60 direct quotient charts with one nonpivot entry \(U_{k\ell}=a\), it
verifies the exact-rank-two pivot minor and the localized identity
\(\mu_1=k!(4-k)!aB_{\ell k}\ne0\).  Hence every dense single-shear torus
already fails the first pure-moment premise.
The double-shear checker enumerates 150 labelled charts and 78 reversal
orbits, then treats one different-row/different-column representative.
Using only exact rational interval arithmetic, it proves a Krawczyk
inclusion for a unique real coefficient-torus zero of
\(\mu_1,\ldots,\mu_8\), verifies exact coefficient rank two, and proves
\(\mu_9/(37!)>0\) on the whole isolating box.  This certifies and kills one
real finite-prefix component; it does not classify the remaining complex
fibre or the other 77 double-shear orbits.

The cross-degree rank-stratified finite-prefix census is replayed by

```bash
python3 scripts/verify_rank_stratified_moment_census.py
```

For every determinantal rank stratum through rank four in balanced
degrees two through four, the checker computes the exact
diagonal-\(\mathrm{SL}_2\) invariant Hilbert coefficients through degree
\(85\) and a good-prime Jacobian of the dimension-sized consecutive
moment system. It recovers the certified degree-four rank-two coefficient
\(-5266\) and ambient degree-three coefficient \(-2186\). The new
rank-two cubic calculation has full nine-moment Jacobian and
\[
[t^{29}]H_{3,2}(t)\prod_{m=1}^{9}(1-t^m)=-58.
\]
Since the first four moments already cut the cubic rank-one Segre cone
down to the nullcone, the resulting semistable nine-moment point has
exact coefficient rank two. It is existential and finite-prefix only.
The first tested single replacement with full Jacobian and a nonnegative
candidate numerator through degree \(85\) is
\(\mu_1,\ldots,\mu_8,\mu_{12}\); orders ten and eleven do not pass. See
[`RANK_STRATIFIED_MOMENT_PROGRAM.md`](extended-geometry/RANK_STRATIFIED_MOMENT_PROGRAM.md).

The first holonomic probe on the degree-three rank-two row is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_holonomic_probe.py
```

It verifies the exact beta/constant-term period formula at two integral
exact-rank-two points, compiles the adjacent C++ moment engine in a
temporary directory, and computes 501 normalized moments at both points
over three primes.  In all six cases an order-\(27\), \(m\)-degree-\(11\)
recurrence ansatz uses 335 equations and passes 139 unused equations.  The
monic forward coefficient is common to every probe, reconstructs exactly,
and has no nonnegative integer root.  An order-\(27\), degree-\(10\)
ansatz fails at both points modulo \(1000003\).  These are exact bounded
computations and modular evidence for a universal recurrence shape, not a
creative-telescoping certificate.  The universal parameter denominator,
exceptional-locus stratification, and corrected-system bridge moments
remain open.  See
[`TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md).

The relative-cohomology refinement is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_relative_jacobian.py
```

It computes the exact logarithmic Jacobian length
\(18=2_{t=0}+2_{t=1}+14_{\rm interior}\), including the saturation
exponent six and an explicit eighteen-monomial basis.  At the same two
points and three primes it finds an order-\(18\), \(m\)-degree-\(18\)
recurrence with 83 unused equations; degree 17 fails at both points
modulo \(1000003\).  Its forward coefficient has eight common linear
factors and a point-dependent decic.  The leading \(m^{18}\) coefficient
has nonzero remainder in all eighteen relative-Jacobian coordinates, so
a naive degree-17 polynomial divergence certificate cannot prove the
recurrence.  The result remains exact quotient data plus modular
recurrence evidence, not a universal telescoping certificate.

The characteristic-zero cyclic splitting at the first integral rank-two
point is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_interior_cyclic_split.py
```

The checker proves the pairwise-comaximal critical-algebra decomposition
\(18=14_{\rm interior}+2_{t=0}+2_{t=1}\).  Eliminating
\(P=Q/u^3\) gives pairwise-coprime squarefree polynomials of degrees
\(14,2,2\), whose product is the degree-\(18\) relative eliminant up to
scalar.  Hence \(1,P,\ldots,P^{13}\) is an exact cyclic basis of the
interior algebra.  It also proves
\[
 ((uQ_u-3Q,tQ_t):(ut)^\infty)=I_{\rm interior},
\]
so the exact toric logarithmic critical rank at this fiber is \(14\).
Exact evaluation of the Chinese-remainder idempotents
against \(\nu_0,\ldots,\nu_{18}\) gives nonzero values on both endpoint
pairs.  Therefore the raw period does not descend through the ordinary
Jacobian quotient.  The same run computes the exact first divergence
seed
\[
u^{47}p_{\rm int}(P)=X(uQ_u-3Q)+YQ_t.
\]
The lift has 6750 and 6791 terms, its divergence has 6749 terms, and the
two nonzero endpoint restrictions have 45 and 85 terms.  An
\(m\)-dependent reduction of these retained terms is still required.
Their first ordinary normal forms occupy all \(14+2+2\) coordinates;
the discarded gradient parts must be lifted recursively rather than
silently set to zero.
The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_interior_cyclic_split.json`.

The exact fixed-fiber rational \(D\)-module seed is replayed by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs.py \
  --annihilator-only
```

This command requires Macaulay2, the `BernsteinSato` package, and
Normaliz.  It rewrites the generating integrand as
\(u^2/(u^3-zQ(u,t))\), computes 34 first-order annihilators of
\((u^3-zQ)^{-1}\), and then computes the 76-generator annihilator of the
specific numerator \(u^2\).  Both left ideals are checked exactly to be
holonomic of rank one.  The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs_research.json`.
This is an all-order integrand certificate, not yet the integrated
Picard--Fuchs operator.  Omitting `--annihilator-only` requests the long
sequential \(t,u\) pushforward; its output must still pass the relative
endpoint audit.

The exact shift-Ore comparison of the two sampled recurrence shapes is
replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_ore_gcd.py
```

At the same two points and three primes, the checker recomputes the
order-\(18\) and order-\(27\) operators through moment 500 and performs
left Euclidean division in
\(\mathbb F_p(m)[S;\,Sf(m)=f(m+1)S]\).  The order-\(27\) operator is not
a left multiple of the order-\(18\) operator.  Their greatest common
right divisor has order \(14\), primitive coefficient degree \(58\), and
left quotient orders \(4\) and \(13\), respectively.  The primitive
order-\(14\) operator is checked directly on all 487 available moment
rows at every sample.  Its order matches the exact interior length in
the \(2+2+14\) relative-Jacobian decomposition.  The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_gcd.json`.
This is an exact bounded modular factor calculation, not a universal
Picard--Fuchs certificate.

Three bounded research probes test shortcuts from the sampled factor to
a relative telescoping certificate:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_ore_reconstruct.py \
  --primes 1000003 1000033
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --steps 1
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_leading_syzygy.py
```

The first exposes all 885 modular coefficients of the fixed-fiber
order-\(14\), degree-\(58\) factor and tests balanced reconstruction at a
held-out prime.  The latter two work in the exact length-eighteen
relative quotient modulo \(1000003\).  They prove that the leading
\(m^{58}\) class is nonzero and that the divergence classes of all
leading Koszul corrections \(R(C,-A)\) have rank zero.  Therefore the
direct zero-boundary polynomial certificate cannot start; the full
\(14+2+2\) endpoint-extended connection is required.  These artifacts
record exact modular no-go calculations for those ansätze, not a
characteristic-zero Picard--Fuchs certificate.

The completed all-order certificate at the first point modulo \(1000003\)
is produced in three restartable Laurent-reduction chunks:

The chunks are intentionally local cache material (large checkpoints and raw
Singular certificates), not Git artifacts.  Run the following target; override
`LOCAL_CERTIFICATE_CACHE` to place the cache elsewhere:

```bash
make generate-rank-two-divergence-local
```

Its expanded commands are:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --mode interior --mapped-quotient --steps 20 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m38.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m58_m38.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_research.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --mode interior --mapped-quotient \
  --checkpoint-input \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m38.poly \
  --steps 20 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m18.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m38_m18.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_research_m38_m19.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --mode interior --mapped-quotient \
  --checkpoint-input \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m18.poly \
  --steps 18 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m0.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m18_m0.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_research_m18_m1.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_terminal_syzygy_block.py \
  --output artifacts/local/two_pair_sic_bidegree33_rank_two_terminal_syzygy_block_research.json \
  --R-output artifacts/local/two_pair_sic_bidegree33_rank_two_terminal_syzygy_R.sing
```

The uncorrected descent has a 307276-term terminal residual.  The final
command constructs the 298606-term Koszul correction \(R\) and the exact
identity \(T=QH(R)\).  The complete independent replay is

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_all_order_certificate.py
```

It replays all 58 coefficient identities in Singular, checks every
restart residual, independently verifies the terminal correction, and
proves that both corrected endpoint exponential-polynomials vanish.
The resulting theorem is the order-\(14\) recurrence for every
\(m\geq0\) at this fixed fiber over \(\mathbb F_{1000003}\).  The summary
artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_all_order_certificate.json`.
This does not reconstruct characteristic zero or a generic rank-two
parameter identity.

The finite characteristic-zero lift at the same fixed fiber is rebuilt in
three stages.  The first command is resumable but expensive: from an empty
cache, computing the 205 exact modular images takes roughly 45 minutes on
the reference machine.

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_ore_reconstruct.py \
  --point 0 --prime-count 205 --prime-start 1000000 \
  --holdout-count 5 \
  --image-cache \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json \
  --output \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_reconstruct_research.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_simultaneous_reconstruct.py \
  --kind common --prime-count 205 --holdout-count 5 \
  --cache \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json \
  --output \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift.json
```

The simultaneous LLL reconstruction uses the first 200 images and reserves
five as fresh holdouts.  It returns a primitive integer order-\(14\),
degree-\(58\) operator whose largest coefficient has 2397 bits.  Generate
the first exact characteristic-zero divergence level with

```bash
make generate-rank-two-char0-leading-local
```

Equivalently:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift.json \
  --mode interior --mapped-quotient --steps 1 --timeout 900 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_char0_interior_divergence_checkpoint_m57.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_char0_interior_divergence_certificate_m58_m57.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_char0_interior_divergence_research_m58_m57.json
```

On the reference machine this producer takes about 255 seconds and peaks
near 1.6 GB.  The combined independent verification is

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_characteristic_zero_lift.py
```

It replays all 205 modular images, 27 exact rational moment identities, the
known forward factor, and the \(m^{58}\) divergence identity over
\(\mathbb Q\).  Its Singular replay takes about 93 seconds and peaks near
1.2 GB.  The manifest is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_characteristic_zero_lift.json`.
This is exact finite verification, not an all-order characteristic-zero
certificate: levels \(m^{57}\) through \(m^1\), the terminal syzygy, and
both endpoint identities remain open.

An experimental reduction-based Picard--Fuchs route is retained in
`scripts/research_two_pair_sic_bidegree33_rank_two_picard_fuchs.jl`.  With
`MultivariateCreativeTelescoping.jl` 0.1.3 installed in an isolated Julia
environment, the generic projective CRT run is

```bash
timeout 900 julia --project=/tmp/sic33-mct-env \
  scripts/research_two_pair_sic_bidegree33_rank_two_picard_fuchs.jl \
  crt 1 original
```

The original-coordinate reference run reached the 900-second cap without an
operator.  The exact beta compression \(x=u, y=ut/(1-t)\) instead gives the
sixteen-term form

\[
 \frac{x+y}{(x+y)^3-z\Phi(x,y)}.
\]

The compact closed-cycle calculation finishes in roughly eight minutes and
peaks near 2 GB:

```bash
timeout 900 julia --project=/tmp/sic33-mct-env \
  scripts/research_two_pair_sic_bidegree33_rank_two_picard_fuchs.jl \
  crt 1 compact \
  artifacts/local/two_pair_sic_bidegree33_rank_two_compact_picard_fuchs.ore
```

It returns a differential order-eight closed-cycle operator.  The interval
period has a different inhomogeneous order-eight relation.  Generate 100 exact
modular images of that relation, the reference-prime structural comparison,
and its resumable local cache with

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_compact_relative_picard_fuchs.py \
  --prime-count 100 --prime-start 1000000 --jobs 4
```

From an empty cache this takes about six minutes on the reference machine.
It proves at \(p=1000003\) that the differential residual has degree 55,
converts its tail to an order-64, \(m\)-degree-eight shift operator \(R\), and
checks the exact modular factorization \(R=Q_{50}G_{14}\).  Reconstruct the
characteristic-zero differential operator from 95 images with five holdouts:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_simultaneous_reconstruct.py \
  --kind compact --prime-count 100 --holdout-count 5 \
  --dimensions 8 12 16 24 32 --offsets 8 \
  --output \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_characteristic_zero_lift.json
```

The exact bridge verifier replays all 100 images, exact rational forcing and
moment rows, 50 exact initial \(G_{14}\)-identities, positivity of the quotient
forward denominator, and the rational shift-Ore identity
\(R=Q_{50}G_{14}\):

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_compact_relative_pf_lift.py
```

The exact Ore division takes about eight minutes but stays below 150 MB.  Use
`--skip-exact-ore-division` for the two-second reconstruction, moment, and
forward-coefficient checks only.

The compact modular all-order certificate is generated in two divergence
chunks and one terminal solve:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json \
  --mode interior --mapped-quotient --steps 1 --timeout 600 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m7.json \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_certificate_m8_m7.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_divergence_m8_m7.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json \
  --mode interior --mapped-quotient \
  --checkpoint-input \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m7.json \
  --steps 7 --timeout 900 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m0.json \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_certificate_m7_m0.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_divergence_m7_m1.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_terminal_syzygy_block.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json \
  --certificate \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_certificate_m7_m0.sing \
  --terminal \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m0.poly \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_terminal_syzygy_research.json \
  --R-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_terminal_syzygy_R.sing
```

The independent combined replay is

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_compact_relative_modular_all_order.py
```

It verifies all eight coefficient identities, the 132615-term terminal
Koszul correction, and the zero endpoint trace, proving \(R\nu=0\) for every
\(m\geq0\) over \(\mathbb F_{1000003}\).  The compact characteristic-zero
operator and its exact factorization are proved, but its divergence and
endpoint identities remain open.  A first expanded rational descent level
reached 900 seconds and 8.6 GB; ordinary top-pole Griffiths reduction leaves
an 18-term remainder, so the next engine must use extended relative reduction
or reconstruct the eight modular certificate levels.

For the two-prime support scout, generate a second operator artifact with

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_compact_relative_picard_fuchs.py \
  --prime-count 1 --prime-start 1000003 --jobs 1 \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research_mod1000033.json
```

Repeat the two divergence commands and the terminal block solve above with
that operator, replacing `mod1000003` by `mod1000033` in every local output.
Then compare the exact Laurent supports with

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_compact_certificate_support.py
```

All \(Y_r\) supports, five of eight \(X_r\) supports, and the terminal support
agree; the other three \(X_r\) supports differ by one monomial each.  This is
an exact two-prime feasibility scout, not a rational reconstruction.

The exact border-basis calculation on the generic factor pencil is
replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_pencil_border.py
```

Over each of three finite rational-function fields it verifies
saturation exponent six, quotient length eighteen, and six reduced
nineteen-term border relations.  The three distinct monic coefficient
denominators have degrees \(74,88,94\), common gcd degree \(74\), and
coprime quotient degrees \(14,20\); their lcm is squarefree of degree
\(108\).  The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_pencil_border.json`.
The same run specializes all \(4,6,5\) base-field roots of that polynomial.
Every specialization remains coefficient-rank two.  Four roots of the
degree-\(74\) component preserve the full \(2+2+14\) relative profile,
while the ten accessible degree-\(14\) roots and one accessible
degree-\(20\) root lower one endpoint or interior length and have total
length seventeen.  This is an exact modular one-pencil classification, not
a universal parameter-space determinant or a classification of the
non-linear exceptional closed points.

The generic-pencil interpolation stress test is

```bash
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_rank_two_recurrence_line.py \
  --samples 256 --holdout 12 --jobs 4 --maximum-moment 390
```

The scaling-family control reconstructs the predicted homogeneous
degrees.  On the generic quadratic factor pencil, eight representative
coefficients have no rational interpolant in the sampled degree window
(combined degree at most 243).  This redirects the universal computation
to an eighteen-dimensional relative connection and determinant
representation instead of an expanded parameter interpolation.

The second checker works on the sixteen-dimensional rank-at-most-two
determinantal variety. At one exact rank-two factor chart it proves modulo
the good prime \(1000003\) that both
\(\mu_1,\ldots,\mu_{13}\) and
\(\mu_1,\ldots,\mu_{12},\mu_{14}\) have Jacobian rank thirteen. It then
computes the diagonal-\(\mathrm{SL}_2\) invariant Hilbert coefficients
from the two-row Cauchy decomposition. Degrees \(1,\ldots,13\) fail the
homogeneous-parameter test exactly:
\[
[t^{69}]H(t)\prod_{m=1}^{13}(1-t^m)=-5266.
\]
Thus their common rank-at-most-two zero fiber necessarily contains a
semistable point. The corrected degrees \(1,\ldots,12,14\) pass the
necessary Hilbert test through degree \(100\). Their numerator is
palindromic through degree \(82\), the top degree predicted by the
determinantal invariant ring's \(a\)-invariant \(-10\) and the candidate
parameter-degree sum \(92\), but their zero fiber remains open. The
rank-one finite-cutoff checks below close all collided-root strata, but
one squarefree uniform-specialization chart remains open. Thus exact rank
two is not yet forced for the semistable thirteen-moment point. No
all-order counterexample is claimed.

The third dependency-free command audits why recurrence work is currently
parked. It verifies
\[
\frac{\mu_m}{(4m+1)!}
=\operatorname{CT}_u\int_0^1
\Phi_C(1,u,t,(1-t)/u)^m\,dt
\]
on the rank-two factor chart and records the corresponding rational
generating function. It also proves that the displayed exact rank-two
Jacobian point has \(\mu_1=7414\), so it is not the existential
thirteen-moment survivor, and computes the generic rank-two Newton polygon
with normalized volume \(48\). The semistable fiber has no recorded closed
point or residue field, so no coefficient-specialized scalar recurrence
is claimed and \(\mu_{14}\) is not evaluated. See
[`TWO_PAIR_SIC_BIDEGREE44_RANK_TWO_ALL_ORDER_AUDIT.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_TWO_ALL_ORDER_AUDIT.md).

On the rank-one boundary, the same checker proves that moments one through
six have exact Jacobian rank six and computes their nonnegative Hilbert
numerator, of coefficient sum \(50\). Rank-one annihilator tensors are
Hilbert--Mumford unstable: for the one-parameter subgroup adapted to the
annihilating root, the missing extremal coefficient makes every surviving
tensor weight strictly positive.

The third checker refines that rank-one boundary on the squarefree
operator chart. After \(\mu_1\) and the normalization \(e=1\), it proves
that the coefficient pivot of the \(a\)-linear second moment has empty
common boundary with \(\mu_2,\ldots,\mu_6\) over
\(\mathbb Q(\lambda)\). On the pivot-open chart the three expected
annihilator sections have ideal
\[
(8c-3d^2,\ d(d-4)(d-4\lambda)),
\]
and the last four moments contain the cubic factor with multiplicities
\(1,1,2,2\) after imposing the first relation. The checker then computes
the gcd of all pairwise \(c\)-resultants. Its only extra quadratic factor
is supported on the pivot divisor because \(p^3\) reduces to zero there,
and the other three branches force
\((c,d)=(0,0),(6,4),(6\lambda^2,4\lambda)\).
Thus the generic squarefree fiber is exactly the annihilator sections.
The same checker closes the pivot-annihilator orbit
\(\lambda^2+4\lambda+1=0\) and the equianharmonic orbit
\(\lambda^2-\lambda+1=0\): each projective fiber has degree four and
the expected four-point radical, with eighth-power certificates.  The
harmonic orbit is closed by the separate exact \(\lambda=2\) anchor.
It also extracts the complete pivot-boundary exceptional gcd and the
three expected-branch gcds. Their sole new \(S_3\)-orbit is represented by
\[
22\lambda^4-54\lambda^3+\lambda^2-54\lambda+22=0;
\]
the checker closes that quartic-field fiber with the same degree and
radical and with eighth-power certificates. On the chart
\(8c-3d^2=0,\ d(d-4)(d-4\lambda)\ne0\), direct substitution and exact
division by the invertible cubic powers reduce the problem to a
three-variable unit ideal over \(\mathbb Q\). One Rabinowitsch membership
remains: on \(p\ne0\), the \(8c-3d^2\ne0\) chart should be supported only
at \(\lambda=0,1\).  Writing \(q=\lambda^4(\lambda-1)^4\) and
\(M=p(8c-3d^2)\), the checker reduces this to the target-only membership
\(qM^5\in(f_3,f_4,f_5,f_6)\).  Modulo \(101,103,107\), the least
saturation exponent is consistently \(5\) and the degree-order basis has
size \(87\).  At \(101\), the four lifted multipliers have degree/term
profiles `(34,5356)`, `(29,3679)`, `(27,3037)`, and `(22,1853)`.
These are exact finite-field identities, but the checker does not promote
them to a characteristic-zero certificate.

The target-only lift experiments can be run separately:

```bash
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py --prime 101
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --prime-lift 101 5
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py --direct 5
```

The first two commands reproduce the exponent, basis size, and finite-field
lift profile.  The third asks for the exact rational target-only lift; the
recorded run exceeded its 1,200-second bound.  A timeout is not evidence
against membership.

The resumable large-prime CRT experiment is recorded in
`artifacts/generated-results/two_variable_quartic_squarefree_crt.json`.
To rebuild it independently, use a new checkpoint path:

```bash
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --crt-lift 5 /tmp/two-variable-quartic-squarefree-crt.json \
  1000003 1000033 1000037 1000039 1000081 1000099 1000117 1000121 \
  1000133 1000151 1000159 1000171 1000183 1000187 1000193
```

The first thirteen images have a common 14,508-term support.  Twelve
build primes give a 240-bit modulus and `1000183` is the holdout; only
three balanced rational reconstructions agree at the holdout.  The final
two primes have different supports.  This rejects coefficientwise CRT of
these arbitrary lifts at the recorded bound, not membership over
\(\mathbb Q\).  The two attempted normalizations are:

```bash
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --compare-tracked-lifts 5 101 103 107
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --compare-syzygy-lifts 5 std component 101
```

The tracked transformations finish but have three different support
hashes.  The component-order syzygy normalization reached its 600-second
per-prime timeout at \(101\).

The fourth checker handles the remaining at-most-two-root normal forms
\(u^rv^{4-r}\). For \(r=0,4\), the first moment gives the one-sided
hyperplane. For \(r=1,2,3\), the first four moments have exactly the
expected two one-sided linear components, with eighth-power radical
certificates. Combined with the existing five-moment three-root theorem,
this closes every collided-root rank-one stratum. The single squarefree
uniform-specialization gate described above remains.

The dual-linear two-pair theorem is replayed by

```bash
python3 scripts/verify_dual_linear_sic2.py
```

This dependency-free audit accompanies the
[dual-linear `SIC(2)` theorem](extended-geometry/DUAL_LINEAR_SIC2.md).
For every \(p=w\mathbin{\cdot}H\), the first two contractions force
\(\operatorname{tr}JH=\det JH=0\), hence
\(H=c+(b,-a)f(ax+by)\).  If \(d=\deg f\) and
\(G=\deg_{x,y}g\), the proof gives
\(\mathcal E_2(gp^m)=0\) for \(m>(d+2)G\).  The normalized Keller case
needs only the first contraction and retains the sharper cutoff \(m>G\).
The checker verifies the second-contraction identity and replays both
cutoffs on exact integer examples.

The unrestricted bidegree-\((2,2)\) theorem is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_bidegree22_frontier.py
```

The checker proves that the natural four-parameter linear compression of
the three-pair witness is forced onto a strict one-sided weight branch by
its first two pure contractions.  Globally, it eliminates a four-parameter
finite-flag chart of the pair-linear one-sided nullcone and obtains twelve
generators for its prime ideal.  If \(I\) is generated by the first six pure
contractions and \(J\) is this nullcone ideal, exact reductions prove
\(I\subseteq J\), \(j_1\in I\), and \(j_r^7\in I\) for
\(2\le r\le12\).  Thus \(\sqrt I=J\), proving SIC(2) for every
bidegree-\((2,2)\) form, including dense eight- and nine-term forms.

As an independent sparse regression, the checker also enumerates all exact
supports of size at most seven, certifies the eight six-term and twelve
seven-term Laurent-curve charts, and verifies their hidden one-sided
factorizations.  The proof and exact claim boundary are in
[`TWO_PAIR_SIC_BIDEGREE22_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE22_FRONTIER.md).

The first mixed-bidegree ordinary-degree search is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_mixed_22_13_31.py
```

For
\[
F=F_{2,2}+F_{1,3}+F_{3,1},
\]
the checker enumerates all \(35\) nonconstant bidegrees of ordinary degree
at most seven and all \(256\) primitive positive/negative central circuits,
then treats this complete three-block stratum exactly.  On
the nonzero \(F_{1,3}\) branch, the dual-linear normal form gives
\(F_{1,3}=\xi _2z_1^3\).  Polynomial-valued contractions through order
three force the aligned triangular flag.  The remaining scalar
contractions through order four reduce to three polynomials whose exact
lexicographic basis contains \(h^3\), and whose \(h=0\) boundary contains
fourth powers of the other two central variables.  The surviving support
has a strict two-step exponent cone, giving an explicit all-multiplier
cutoff rather than a bounded moment prefix.  The \(F_{1,3}=0\) branch
reduces to the bidegree-\((2,2)\) theorem.  The proof and the still-open
ordinary-degree-\(<8\) collection classes are in
[`TWO_PAIR_SIC_ORDINARY_DEGREE_LT8_MIXED_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_ORDINARY_DEGREE_LT8_MIXED_FRONTIER.md).
The same replay checks the central factorial formula for
\(V_{1,d}\oplus V_{d,1}\), \(2\leq d\leq6\); the written exponent-cone
proof is uniform for every \(d\geq2\).

The next dual-linear mixed branch is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_mixed_33_14_41.py
```

For
\[
F=F_{3,3}+F_{1,4}+F_{4,1},\qquad F_{1,4}\ne0,
\]
the checker normalizes \(F_{1,4}=\xi _2z_1^4\), verifies the maximal-weight
linear ladder, and uses two full-polynomial output coefficients from
moments three and four to force the \(F_{3,3}\) coefficient matrix upper
triangular.  It constructs the four-variable diagonal/opposite-corner
core and its moments through order five.  Their exact 28-element rational
Gröbner basis contains tenth powers of all three diagonal parameters and
the fifth power of the corner parameter.  The residual exponent cone
then gives the explicit multiplier cutoff.  The checker does not close
the \(F_{1,4}=0\) boundary: the written proof identifies that boundary
exactly with the existing balanced bidegree-\((3,3)\) problem.

The first degree-eight diagonal-core probe is

```bash
.venv/bin/python scripts/explore_two_pair_sic_mixed_diagonal_core.py
```

The script uses a direct factorial composition formula for the core of
\(V_{4,4}\oplus V_{1,5}\oplus V_{5,1}\), avoiding expanded Wick
polynomials.  Moments two through six have exact rational Jacobian rank
five.  Singular gives basis size \(132\), quotient dimension \(360\), and
origin-support power reductions over both \(\mathbb F_{101}\) and
\(\mathbb F_{1009}\).  These are exact finite-field computations but only
evidence for the rational radical.  The script neither reconstructs a
characteristic-zero certificate nor checks the preceding \(d=4\)
triangularization equations.

The first non-dual-linear positive block is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_mixed_23_32_pure_summands.py
```

The checker constructs
\(V_{2,3}=\operatorname{Sym}^5\oplus\operatorname{Sym}^3
\oplus\operatorname{Sym}^1\) exactly.  It identifies the first contraction
with the linear projection, the pure cubic second-moment ideal with the
rational normal cubic, and the pure quintic ideal with the prime
tangential variety \(L^4M\), including an exact Singular primary
decomposition.  It then inserts all three nonzero orbit normal forms into
\(V_{2,3}\oplus V_{3,2}\).  Moments two through four eliminate every
negative-block coefficient capable of an unbounded central contraction,
and explicit exponent bounds prove eventual mixed vanishing.  The branch
with both positive irreducible summands nonzero remains open.  At one
explicit mixed-positive point, the same checker constructs the full
moment-one-through-four Jacobian.  Its two transverse excess directions
are obstructed at deformation orders two and four, respectively.  For the
second direction it retains all six parameters in both higher correction
spaces and verifies that the two order-four compatibility quadratics have
nonzero resultant
\(2283980165392458318151680000\).  This is an exact local jet exclusion,
not a global classification of the remaining branch.

The one-parameter generic strengthening is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_mixed_23_32_generic_local_gate.py
```

The raising operator removes the third incidence coefficient on
\(b(a-b)\ne0\), leaving \(B(u)=VZ^2(uVY+WZ)\).  Over
\(\mathbb Q(u)\), the checker retains all second- and third-order
corrections and factors the fourth-order resultant as
\[
648u^8(u+6)^4(4u+3)^2
(2u^3+10u^2+21u+9)^2S(u),
\]
where the displayed artifact records the irreducible sextic \(S\).
The linear subresultant also gives the unique fourth-order correction
\(\tau=-B(u)/A(u)\) at every root of \(S\), since the checker verifies
\(\operatorname{Res}(A,S)\ne0\).
It also recomputes \(u=0,-6,-3/4,\infty\) and one \(a=b,c\ne0\)
representative without the rational-function chart; every resulting
fourth-order obstruction ideal is the unit ideal.  This is a generic
local first-four-moment theorem, not a global SIC classification.

The balanced cubic two-variable GVC theorem is replayed by

```bash
.venv/bin/python scripts/verify_two_variable_cubic_gvc.py
```

For the three nonzero binary-cubic symbol orbits \(u^3\), \(u^2v\), and
\(uv(u+v)\), the checker derives the apolar moments directly from
coefficient expansion.  It verifies that moments through orders one,
three, and four, respectively, leave only the displayed one-sided normal
forms.  In the squarefree orbit it checks all three exact branch
factorizations and the constant annihilating direction of every surviving
pure cube.  The written degree argument gives
\(\Lambda^m(QP^m)=0\) for \(m>\deg Q\); this all-order conclusion, not a
bounded replay, proves the theorem.  See
[`TWO_VARIABLE_CUBIC_GVC_THEOREM.md`](extended-geometry/TWO_VARIABLE_CUBIC_GVC_THEOREM.md).

The quartic three-root continuation is replayed by

```bash
.venv/bin/python scripts/verify_two_variable_quartic_three_root_gvc.py
.venv/bin/python scripts/verify_two_variable_quartic_squarefree_generic.py
```

After normalizing the \((2,1,1)\) symbol to \(u^2v(u+v)\), the first
moment eliminates one polynomial coefficient.  The checker derives moments
two through five, proves containment in the three-component one-sided
ideal, and verifies that the fourth powers of all five generators of that
radical lie in the moment ideal.  The same
[`low-root note`](extended-geometry/TWO_VARIABLE_LOW_ROOT_GVC_THEOREMS.md)
proves the all-degree theorem for symbols with at most two roots via the
one-variable Duistermaat--van der Kallen constant-term theorem; that part
is a written proof rather than a bounded computation.  The squarefree
checker verifies the four symbolic annihilator sections and proves that at
cross-ratio \(2\), moments one through six have exactly the four expected
reduced projective zeros.  Proper-family upper semicontinuity then proves
the same equality on a nonempty Zariski-open set of cross-ratios.  This is
a generic theorem; finitely many exceptional squarefree orbits remain
possible.

The later
[`SPLIT_SYMBOL_GVC_THEOREM.md`](extended-geometry/SPLIT_SYMBOL_GVC_THEOREM.md)
requires no computer algebra.  It factors a homogeneous operator symbol
as directional derivatives and proves
\[
 \Lambda^m(P^m)
 =(m!)^d\operatorname{CT}
 \left(\frac{P(z+\sum_i t_i v_i)}{\prod_i t_i}\right)^m.
\]
Choose one generic \(z\) exposing the full finite \(t\)-support.
Duistermaat--van der Kallen gives an integral weight separating that
support from the origin, and the same weight works at every specialization
of \(z\).  A fixed translated multiplier cannot cross the linearly growing
gap.  Thus every homogeneous binary operator satisfies GVC for arbitrary
\(P\), including \(\deg P>\operatorname{ord}\Lambda\).  There is no
bounded replay to confuse with the proof.

The nonhomogeneous lowest-order extension and the rank obstruction to
natural conversions of the two-pair witness are checked by

```bash
python3 scripts/verify_separable_gvc_escape_obstructions.py
.venv/bin/python scripts/verify_binary_heat_quadratic_gvc.py
```

The first dependency-free command also rewrites
`artifacts/generated-results/separable_gvc_escape_obstructions.json`; its
current whole-file SHA-256 is
`3343e46cca1b9459f0a3f113278d1db610379e1c7083370290f30f32e420f226`.

The finite checker verifies that the witness matrix has determinant \(48\)
and rank five, that four separated rank-one channels cannot reach it, and
that coefficient extraction is not multiplicative.  The written proof in
[`SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md`](extended-geometry/SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
shows that if \(r\) is the lowest positive order of a nonhomogeneous binary
operator, then \(\deg P\leq r\) implies GVC.  Higher operator pieces can
enter a mixed value only a bounded number of times, and the split-symbol
Newton separator absorbs those bounded defects.  This is an all-order
argument, not a bounded search.  The same note records the arbitrary-degree
factor-unit extension \(\Lambda=\Lambda_0\Gamma\), with \(\Lambda_0\)
homogeneous split and \(\Gamma(0)\ne0\).  For a completely general
nonhomogeneous binary operator, every fixed number of leading homogeneous
mixed layers still vanishes eventually; any defect must move to unbounded
filtration depth.  The note also proves that fixed linear translation plus
one diagonal coefficient is possible exactly for split homogeneous
symbols.  Finally, it closes the binary linear-plus-quadratic heat class
for \(\deg P\le2\): after normalizing the linear part, the only irreducible
coefficient is \(C\partial_y^2\), the first equation gives
\(P=c(y^2-2Cx)+ey+f\), and the second is \(16C^2c^2\).  The surviving
affine transverse form has cutoff
\(m>2\deg_xQ+\deg_yQ\).  The second command derives the first-moment
normal form and checks the universal second-moment identity symbolically.
It also constructs the generic degree-six heat-harmonic polynomial and
checks
\(\Lambda^2(P^2)=4C^2(P_{yy})^2\).  The written product-rule proof gives
the corresponding highest-\(y\)-degree square in every degree and closes
every binary operator with nonzero linear part and no terms above order
two.  The checker also verifies the generic cubic drift--diffusion normal
form.  The cutoff is the written derivative count.
The same checker tests the degree-six family for
\(\partial_x+C\partial_y^2+E\partial_y^3\).  The written iterated
product-defect proof closes every
\(\partial_x+h(\partial_y)\), with cutoff
\(m>r\deg_xQ+\deg_yQ\) when \(r\) is the lowest order of \(h\).
The proof permits formal \(h\), since its differential action is locally
finite.  Formal Weierstrass division writes every binary symbol with
nonzero linear part as \(U(\xi,\eta)(\xi+q(\eta))\); the unit
\(U(\partial)\) and its inverse are locally finite differential
automorphisms.  This reduces the entire lowest-order-one class to the
separated theorem for arbitrary \(P\).
The symbolic checker also adds a completely general cubic operator piece
and verifies that the quadratic-\(P\) second moment remains
\(16C^2c^2\).  The written order argument extends this to arbitrary
higher pieces, proving safety for every \(\deg P\le2\) when the linear
operator part is nonzero.
For a general operator and cubic polynomial, the checker solves the first
equation and retains every order-four and order-five term that can occur
in the second moment.  It verifies the decisive coefficients
\[
144C^2p_9^2,\qquad 16C^2p_5^2,\qquad 648G^2p_9^2.
\]
Higher operator orders kill \(P^2\), and the retained terms cannot change
these branches.  This gives the written theorem for every binary operator
with nonzero linear part and every \(P\) of degree at most three.
The checker additionally retains the complete operator \(7\)-jet for a
quartic \(P\) and verifies the successive branch coefficients
\(2304C^4p_4^2\), \(15552G^2p_4^2\), and \(39168L^2p_4^2\).
These finite-jet calculations are regressions for the stronger formal
straightening theorem, rather than the source of its all-degree proof.
For lowest positive order two and cubic \(P\), it separately checks both
quadratic-symbol orbits.  In the double-line orbit it verifies the unique
second-moment cancellation and the decisive third moment
\(-4608C^3p_{xy^2}^3\); in the distinct-line orbit it verifies the
triangular second-moment branches.  The written strict weighted-degree
cutoffs turn the surviving branches into an all-order theorem.
For quartic \(P\), the checker retains every operator jet that can enter
the displayed moments.  It verifies the triangular distinct-root
second-moment closure, the double-line \(xy^3\) third-moment coefficient
\(-3604176H^3\), and the two residual polynomials \(S,T\) on the
double-line \(y^4\) branch.  Their exact resultant is a monomial times the
homogeneous sextic printed in Proposition 3.8, leaving only finitely many
nonweighted parameter ratios for the next moment.
It then evaluates moment four, eliminates the remaining coefficient, and
verifies that the resulting octavic has gcd one with the sextic.  The
written \((2,1)\)-weighted-face identity shows that higher operator jets
and lower polynomial terms cannot alter any pure moment, upgrading this
calculation to closure of the full \(r=2,\deg P\le4\) cell.
For \(r=3,\deg P=4\) with triple-root leading symbol, the checker verifies
the three stabilizer types.  In the \(x^2y^2\) type it checks the forced
weighted correction chain and terminal coefficient
\(3361505280U^4\); the \(xy^3\) and \(y^4\) types terminate at moments
three and two.  The written weighted-face argument supplies the all-order
mixed cutoffs.  It also checks all double-root branches and the squarefree
orbit.  For the latter it computes the leading-moment Gröbner basis,
reduces the three fourth-power tips by root permutation, retains the full
order-four and order-five jets at the \(x^4\) tip, and verifies the
terminal weighted-face value \(129392640T^3\).  Together with the written
orbit and weight arguments this proves GVC for every binary operator and
every \(P\) of degree at most four.  Finally, the checker verifies the
leading reduction in the first open \(r=2,\deg P=5\) row: the only three
top pairs are
\((\partial_x\partial_y,x^5)\),
\((\partial_x^2,xy^4)\), and
\((\partial_x^2,y^5)\).  It then retains the complete jets that can enter
the first two moments.  The distinct-root branch reduces to
\(P=f(x)+ay\) and
\(\Lambda=\partial_y\Gamma+H(\partial_x)\), with
\(\operatorname{ord}_{\min}H\ge6>\deg f\).  The written binomial
derivative count gives its all-order mixed cutoff.  For the \(xy^4\)
branch the checker verifies two nested third-moment faces, terminating in
\(-553153536H^3\) and \(-5430509568J^3\).  For the \(y^5\) branch it
checks the six weighted second-moment ratios and the four nonzero
third-moment residuals; the other two ratios are one-sided.  This closes
the full \(r=2,\deg P\le5\) cell.  The checker finally derives the next
\(r=3,\deg P=5\) leading reduction: four triple-root, three double-root,
and one squarefree top-form normal forms.  This checker supplies the
leading-face classification only; the separate degree-five frontier
checker below closes all eight nonhomogeneous correction systems.
The accompanying written no-go for formal umbral straightening is
proof-theoretic: conjugation by an algebra automorphism preserves the
Leibniz rule, while a locally finite formal constant-coefficient operator
is a derivation only when its symbol is linear.  No bounded computation is
used for that statement.

The eight cubic-leading quintic normal forms and the squarefree
quartic-leading cross-ratio row are closed by
[`BINARY_DEGREE_FIVE_GVC_FRONTIER.md`](extended-geometry/BINARY_DEGREE_FIVE_GVC_FRONTIER.md).
The default exact checker replays the triangular moment eliminations,
terminal one-sided faces, and explicit mixed-multiplier tails:

```bash
make verify-binary-degree-five-gvc
```

For an exploratory dump of all eight complete second-moment jets, with the
optional product-splitting and unit-pivot heuristic, run:

```bash
.venv/bin/python scripts/explore_binary_degree_five_gvc_frontier.py \
  --triangular-components
```

With Singular, replay both residual radicals and the uniform
squarefree-quartic projective saturation:

```bash
.venv/bin/python scripts/verify_binary_degree_five_gvc_frontier.py \
  --singular --singular-top
```

The separate modular screen enumerates \(6{,}696{,}142\) residual triples
over \(p=101,103,107\) and \(2{,}082{,}612\) squarefree-quartic projective
top forms over seven smaller primes:

```bash
make search-binary-degree-five-gvc
```

It regenerates
[`binary_degree_five_gvc_face_search.json`](artifacts/generated-results/binary_degree_five_gvc_face_search.json).
The modular record is an exhaustive bounded experiment on the displayed
faces; the characteristic-zero radical computations and weighted
face-separation argument are the proof.  The non-squarefree
quartic-leading nonhomogeneous rows remain outside this particular theorem,
so it alone does not give a universal binary degree-five corollary.

The quadruple-root partition \((4)\) is closed separately, with arbitrary
lower polynomial pieces and arbitrary higher Weierstrass operator jets, by:

```bash
.venv/bin/python scripts/verify_binary_quartic_quadruple_root_gvc.py
```

This exact checker replays the defect-one radical, all three minimal
branches, the terminal weight-eight equality chain through pure moment
five, and the final strict or one-sided weight separators.  The proof and
normalizations are documented in
[`BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md`](extended-geometry/BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md).

The other repeated-root quartic partitions are closed by:

```bash
.venv/bin/python scripts/verify_binary_quartic_triple_simple_root_gvc.py
.venv/bin/python scripts/verify_binary_quartic_double_root_gvc.py
```

The first command verifies the \((3+1)\) defect-one radical, every
projective branch through defect three, and its final separators.  The
second closes \((2+2)\) and \((2+1+1)\), including both two-parameter
weight-six equality systems and the isolated fifth-power components.
Together with the earlier degree-five checker, these exact calculations
prove binary GVC for every polynomial of degree at most five.  See
[`BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md`](extended-geometry/BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md)
and
[`BINARY_QUARTIC_DOUBLE_ROOT_GVC.md`](extended-geometry/BINARY_QUARTIC_DOUBLE_ROOT_GVC.md).

The first sextic frontier cell, with lowest symbol of order five and root
partition \((5)\), is closed by:

```bash
.venv/bin/python scripts/verify_binary_quintic_quintuple_root_gvc.py
```

This exact checker verifies the defect-one radical, all five projective
top-form branches, and the complete terminal equality chain through
operator order ten and pure moment six.  The proof includes arbitrary
lower pieces of \(P\) and arbitrary higher operator jets.  See
[`BINARY_QUINTIC_QUINTUPLE_ROOT_GVC.md`](extended-geometry/BINARY_QUINTIC_QUINTUPLE_ROOT_GVC.md).
The complete replay below incorporates this longest component and closes
the other quintic root partitions.

All six remaining quintic root partitions, and hence the complete
\((r,\deg P)=(5,6)\) row, are closed by:

```bash
.venv/bin/python scripts/verify_binary_quintic_all_root_partitions_gvc.py
```

This checker verifies the Hall-matching classification of the leading
pure-zero locus, the local correction systems for root multiplicities one
through five, generic strict quintic cofactors, and every final weighted
separator.  It invokes the quintuple-root replay for the multiplicity-five
cell.  See
[`BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md).
The complete \((r,\deg P)=(4,6)\) row is closed by:

```bash
.venv/bin/python scripts/verify_binary_quartic_all_root_partitions_gvc.py
```

This exact checker verifies the Hall leading-locus classification, every
repeated-root terminal face, both coupled pure-sixth-power endpoint
radicals, the triple-root terminal tail coefficients, and the simple-root
defect layers through defect four.  The finite-tail inequalities covering
arbitrary later operator jets are documented in
[`BINARY_QUARTIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_QUARTIC_ALL_ROOT_PARTITIONS_GVC.md).
The complete \((r,\deg P)=(3,6)\) row is closed by:

```bash
.venv/bin/python scripts/verify_binary_cubic_all_root_partitions_gvc.py
```

This exact checker verifies the cubic Hall locus, the triple-, double-,
and simple-root Newton ladders, and the two coupled weighted-face chart
covers.  It requires Singular and `msolve`; the latter is used over
characteristic zero to prove explicit affine saturations empty.  See
[`BINARY_CUBIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_CUBIC_ALL_ROOT_PARTITIONS_GVC.md).
The complete \((r,\deg P)=(2,6)\) row, and hence universal binary GVC
through polynomial degree six, is closed by:

```bash
.venv/bin/python scripts/verify_binary_quadratic_all_root_partitions_gvc.py
```

This exact checker verifies the quadratic Hall locus, the full
distinct-root first-equation reduction and second-moment ladder, all
half-integral and integral double-line Newton faces, and every primary and
secondary radical at the pure-sixth-power endpoint.  It requires Singular
over characteristic zero.  The arbitrary-jet weight-defect argument is in
[`BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md).
The quadratic-leading row on the next polynomial-degree-seven frontier is
closed by:

```bash
.venv/bin/python scripts/verify_binary_quadratic_septic_gvc.py
```

This exact checker verifies the seven-factor Hall locus, the complete
distinct-root first-equation reduction and second-moment ladder, every
half-integral and integral double-line face, and the full branch tree over
the pure seventh-power endpoint.  In particular it checks the extra
slope-three axis migration which does not occur in degree six.  The written
proof in
[`BINARY_QUADRATIC_SEPTIC_GVC.md`](extended-geometry/BINARY_QUADRATIC_SEPTIC_GVC.md)
uses the final common-threshold coordinate deficits to cover arbitrary
higher jets and a fixed multiplier.  A degree-seven counterexample must
therefore have lowest positive operator order three through six.

The cubic-leading row on the same frontier is closed by:

```bash
.venv/bin/python scripts/verify_binary_cubic_septic_gvc.py
```

This exact checker derives every crossing from the cubic two-wing normal
form, verifies the four exceptional monomial radicals and every origin
radical in their child intervals over \(\mathbb Q\), and audits all terminal
common thresholds.  It requires Singular.  The written proof in
[`BINARY_CUBIC_SEPTIC_GVC.md`](extended-geometry/BINARY_CUBIC_SEPTIC_GVC.md)
uses the uniform terminal-face theorem to absorb arbitrary strict jets and
the fixed multiplier.  A degree-seven counterexample must therefore have
lowest positive operator order four through six.

The remaining three septic rows, and hence binary GVC through polynomial
degree seven, are closed by:

```bash
.venv/bin/python scripts/verify_binary_high_order_septic_gvc.py
```

This exact checker constructs all 46 Hall charts for lowest orders four
through six, verifies 287 initial face radicals and 98 child face radicals
over (mathbb Q), checks the fifteen squarefree axis exceptions, and
audits strict marked-gap descent to every final common threshold.  It
requires Singular.  The proof and complete census are in
[`BINARY_DEGREE_SEVEN_GVC.md`](extended-geometry/BINARY_DEGREE_SEVEN_GVC.md).

Unrestricted binary GVC is proved without a degree census in
[`BINARY_GVC_ENVELOPE_CLOSURE.md`](extended-geometry/BINARY_GVC_ENVELOPE_CLOSURE.md).
The proof combines Hall localization, the shifted-ray endpoint theorem, and
the unequal common-threshold theorem.  It has no new computational premise:
the finite global lower and upper Newton envelopes cannot exchange horizontal
order while their weight gap is positive, and finite support forces that gap
to reach zero.  The first octic Ferrers face which exposed this argument has
the exact regression:

```bash
.venv/bin/python scripts/verify_binary_gvc_ferrers_regression.py
```

It uses `msolve` over characteristic zero to verify the radical
((A,S,T,BP,BQ,CQ)).  The first degree-nine gap-four staircase is a longer
optional replay:

```bash
.venv/bin/python scripts/verify_binary_gvc_ferrers_regression.py --gap-four
```

That optional command checks eight affine saturations and can take several
minutes.  These Ferrers calculations are regressions, not dependencies of
the unrestricted envelope proof.

The all-degree Hall localization and unequal-weight terminal-face theorem
have a dependency-free exact regression:

```bash
.venv/bin/python scripts/verify_binary_gvc_uniform_face_termination.py
```

It exhausts the Hall inequality through order twelve, checks the
coefficient-independent prime-valuation inequalities on small weighted
lattice segments, and verifies an exact weight-\((3,2)\) example whose
first moment cancels but whose fifth moment has the uniquely predicted
5-adic endpoint, together with the generic-translation multifactorial
identity underlying the Newton-intersection criterion.  It also checks the
homogeneous beta-integral identity, the failure of constant-term
extraction to commute with powers in the first two-channel example, and
the exponent-dependent factorial distortion under a toric blow-up.  It
also verifies the exact minimal Bernstein-circuit formula
\(\Phi(F^m)=(ac+bd)^m/(m+1)\), which identifies Long's rank-one
beta--torus circuit as a linear Hall annihilator rather than a GVC
counterexample.  Finally it replays the all-degree primitive cusp
parallelogram obstruction
\[
 E_{r,s}=T_r-T_s+\frac92(C_s-C_r)(C_r+C_s-2)
\]
and the moment-two closure of every sparse four-channel dilation.  The
same checker verifies the five-channel warning
\[
 A=X^2-\frac13Y^3,\qquad
 P=x^2+y^3+\frac{13}{30}+\frac{11}{2}x+xy^3:
\]
its first three scalar moments vanish, no four-channel restriction
inherits those three vanishings, and its fourth moment is \(1\,205\,760\).
This is a finite-prefix obstruction, not a pure-moment-zero pair.  The
checker also replays the all-even unit-line half-bridge theorem.  After
moments two and three determine \(u=de\) and \(v=ce^2\), it checks
\[
 M_4=\frac12\left(
 2Q_n-40C_nT_n+48T_n+81C_n^3-180C_n^2+132C_n-48
 \right),
\]
the exceptional value \(H_2=-480\), and the increasing multinomial ratio
\(Q_n/(C_nT_n)>20\) from \(n=4\) onward.  Finally it constructs the six
quadratic--cubic and two double-quadratic obstruction formulas from the
factorial weights \(W_{m,k}\), checks the double-quadratic determinants,
and finds all eight obstructions nonzero for every unequal
\(1\leq r,s\leq30\).  This last window is an exact regression, not the
all-\((r,s)\) inequality proof for the five expressions left open below.

The return classification and the first unbounded arithmetic certificate
are replayed by

```bash
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --prove-h00
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --prove-negative-corners
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --prove-three-more
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --verify-opposite-packet --limit 40
.venv/bin/python scripts/research_binary_gvc_frobenius_carry.py
.venv/bin/python scripts/research_binary_gvc_frobenius_carry.py \
  --radial-limit 3 --order-limit 2 --bridge-limit 29 --residue-limit 2
.venv/bin/python scripts/research_binary_gvc_ghost_shell.py
.venv/bin/python scripts/research_binary_gvc_quotient_graver.py
.venv/bin/python scripts/research_binary_gvc_witt_rees.py
.venv/bin/python scripts/research_binary_gvc_equal_radial_union.py
.venv/bin/python scripts/research_binary_gvc_reversal_union.py --prove-width-two
.venv/bin/python scripts/research_binary_gvc_torsion_torus_trace.py
.venv/bin/python scripts/verify_binary_gvc_torsion_torus_digit_separation.py
.venv/bin/python scripts/verify_binary_gvc_weighted_trace_obstruction.py
.venv/bin/python scripts/verify_binary_gvc_first_ghost_source_collapse_and_ray_rigidity.py
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --limit 100
```

The second Frobenius-carry command is the focused replay for the affine
singleton-localization theorem.  Its final blocks verify endpoint-fibre
signed-digit stabilization and factorial units at two successive primes,
replay the sharp interior family
\[
 ((p-1)/2,1,(p-1)/2)-(0,p,0)=((p-1)/2)(1,-2,1)
\]
through every odd prime at most \(29\), and verify the radial-carry Hasse
compression formula.  The finite checks are regressions for the written
all-prime Graver, exposed-vertex, and Wilson-unit proofs.  The same note's
global Hilbert-module theorem is a direct Dickson--Gordan argument: all
\((p,y)\) fibres are a finite module over the Hilbert basis of
\(B(p,y)=0\).  It justifies one fixed toric/Graver presentation per support;
the computation is not the source of that finiteness statement.

The ghost-shell command also constructs the complete primitive three-level
affine ghost
\[
 G_{p;u,v}(X)=\sum_{t=1}^{\lfloor p/(u+v)\rfloor}
 \frac{(p-1)!X^t}{(ut)!(vt)!(p-(u+v)t)!}.
\]
For the 12 coprime types \(1\leq u\leq v\leq6\), it tests the primes through
43, the reduced rational window with numerator of absolute value at most 40
and denominator at most 20, every cyclotomic order through 80, and all 2,139
primitive irreducible polynomials of degrees two and three with coefficient
height at most four.  The only non-support survivor is the centered value
\((u,v,X)=(1,1,1)\).  The note proves this for everywhere-good,
Galois-stable root-of-unity candidates at every width and torsion order by
the first-admissible-prime and Euler-\(\varphi\)-degree argument.  It does not
remove an arbitrary extra finite exceptional-prime set.  The finite search is still not an
all-prime classification of arbitrary non-torsion algebraic roots; the
constant-term, Wilson-coefficient, and cyclotomic-separation arguments printed
in the note are the exact theorems.

The first command checks the two closed determinant formulas and the exact
coefficient-positive proof that
\(\mathcal H_{0,0}(r,s)>0\) for every unequal endpoint pair.  Its forward
difference certificate has 932 nonnegative terms and its \(s=1,\ r\ge2\)
base certificate has 19.  This is an all-order symbolic certificate.  The
second command proves \(\mathcal H_{0,3}<0\) and
\(\mathcal D_{0,2}<0\) for every unequal endpoint pair.  It verifies the
three exact step-ratio bounds defining the coupled cone and four
coefficient-positive numerator expansions, of sizes 266, 361, 2236, and
2236.  The third command proves the uniform nonvanishing of
\(\mathcal H_{0,1},\mathcal H_{0,2},\mathcal H_{1,1}\), using monotone
ordered-tail cones, fixed-ray cones, and exact finite complements.  It
then closes \(\mathcal H_{1,0},\mathcal D_{0,1}\) on the last wedge
\(r>s\ge4\), using the increasing product \(L_nM_n/C_n^3\) and
coefficient-positive expansions of 408 and 1692 terms.  The fourth
command verifies the exact second-coefficient identity for the opposite
three-by-three packet and checks its central-binomial ratio dichotomy
through degree and endpoint order 40.  This is a bounded regression; the
unbounded proof is the strict Vandermonde supermultiplicativity argument
in Theorem 7.5 of the canonical note.  The fifth command checks the
binary Frobenius carry gap and its normalized unit formulas on 28,858
exact homogeneous return types, checks the nonhomogeneous jet--carry
score on 477 mixed-degree types, and constructs 111,930 one-sided triples and
5,787,067 two-sided pair-pairs through radial degree 40.  The unbounded
proof in Lemma 7.4 ter and Corollary 7.4 sexies is the corresponding
Legendre--Kummer and finite-field kernel calculation.  The sixth command
constructs the centered-triple and two-by-two ghost diagonal blocks.  It
verifies their universal factors \(X-1\) and \(X(X+1)\), respectively,
and displays their prime-dependent residual factors through prime 43.
It constructs the exact characteristic-zero beta diagonals and verifies
their common gcd \(X(X+1)(X^2+X+1)\) across the tested primes.  Its rational
cross-prime search of height 20 leaves only \(1\) in the centered block
and \(0,-1\) in the beta block, but does not see the algebraic cube-root
branch; that window is only a regression.  Proposition 7.4 quater proves
the universal and persistent factors for all primes at least five.  The
same command verifies
the terminal augmented blocks: the beta ordinary row \(1+X\) leaves only
the Hall value, and the centered Bessel endpoint rows \(U,U^2+2V\) have
Jacobian determinant \(2\) and force support loss.  Corollary 7.4
quinquies closes the isolated atom arithmetic; compatibility with the
common high-digit quotient remains unproved inside the parked Hall/carry
route.  The Hall-envelope proof does not need it.  The seventh command verifies
the first obstruction to circuit-only quotient peeling,
\[
 R_3B_1B_2=R_0B_3^2.
\]
On the projected support \(R=\{0,3\}\), \(B=\{1,2,3\}\), these are the
only two states of their color-count/level fiber, and their support-five
difference has no circuit move.  The checker verifies primitivity and
finds this as the first such two-color identity, up to reversal.  It also
checks the two terminal completions through \(R_2B_2^2\) and
\(R_3B_0B_3\).  Proposition 7.4 septies records the projected-scroll
obstruction and the finite Gröbner bound; Corollary 7.4 octies closes
this first block by the circuit-completion/radial-reversal dichotomy.
The checker also replays the explicit \(S(6)\) and \(S(5,4)\)
non-universal-Gröbner witnesses from
Bogart--Hemmecke--Petrović.  The projected \(S(5,4)\) fiber has exactly
four states, circuit-component sizes \(1+3\), and repeated-ray
factorial signatures \((2,2),(2,1,1),(1,1,1,1)\) with Stirling bases
\(16,4,1\).  Proposition 7.4 nonies proves that every whole exposed
scroll profile is terminal by a one-variable constant-term reduction.
The eighth command tests the first two normalized \(p\)-typical Witt
coordinates at \(p=13\) on two unit coefficient specializations of the
support-five, \(S(6)\), \(S(5,4)\), and first larger
reversal-symmetric packets.  In every case the first residual has
valuation exactly one, while both \(\mathcal G_2-\mathcal G_1\) and the
second Witt residual have valuation exactly two.  Proposition 7.4
decies proves the all-height result: the Laurent constant-term factor
and the signed normalized radial factorial are Gauss sequences, so
every exposed packet has an integral \(p\)-typical Witt recursion.
Ghost injectivity does not split a vanishing sum of several profile
Witt vectors; the remaining open step is a profile-separating Rees
initial idempotent or a separator/support-loss consequence of
least-profile cancellation.  The ninth command verifies the complete
equal-radial union identity directly through scale four.  It also finds
the smallest persistent failure of color-count saturation:
\(R=\{0,2\}\), \(B=\{1\}\), whose achievable red counts at scale \(N\)
are \(0,2,\ldots,2N\).  Despite those holes, the whole union is one
Laurent constant-term sequence; its first two symbolic rows reduce to
\(-14a^4\).  Proposition 7.4 undecies proves that color-count
saturation is unnecessary once all states at one oriented radial
vector are exposed.  It does not prove Hall/jet exposure of that
complete union or separate a coordinate-reversed pair.  The tenth
command saturates the four characteristic-zero endpoint charts for the
first coordinate-reversed Laurent width.  For support in \([-2,2]\)
and target slopes \(\pm1\), the charts
\((-1,1),(-1,2),(-2,1),(-2,2)\) close at rows \(2,4,4,8\).
It also exhausts the projective five-coefficient space modulo
\(5,7,11\), with no survivor through row eight.  The rational
saturations prove the width-two statement; the modular census is only
a regression.  Arbitrary reversal width is closed by the finite-trace
digit theorem below.  The eleventh
command verifies the exact regular-representation identity
\[
 \operatorname{CT}_{\mathbb Z\times C_q}(u^N)
 =q^{-1}\operatorname{CT}_{\mathbb Z}
   \operatorname{Tr}(\operatorname{Reg}_{C_q}(u)^N),
\]
its log-determinant generating series, and the compatible Frobenius
identity for \(p\equiv1\pmod q\), on exact \(C_2\) and \(C_3\)
examples.  It also diagonalizes the width-two reversal packet into its
two \(C_2\)-character components.  The twelfth command verifies signed
base-\(p\) digit uniqueness directly in one and two free Laurent
variables:
\[
 \operatorname {CT}(f^{n_0+n_1p+\cdots+n_sp^s})
 \equiv
 \prod_j\operatorname {CT}(f^{n_j})^{p^j}\pmod p.
\]
It also replays the Newton-identity endpoint for two, three, and five
trace components.  The characteristic-zero proof uses arbitrarily large
completely split primes: repeated equal digits recover every power sum
of the component moments, forcing componentwise vanishing.  Together
with character orthogonality and Duistermaat--van der Kallen, this closes
the identity-coefficient torsion--torus trace lemma and every
scale-compatible carry packet.  It does not prove that a Hall--jet shell
has that form.  The thirteenth command verifies the weighted mixed-digit
extension and two exact obstructions to that promotion.  The dilation
pair \(z+z^{-1}\), \(z^2+z^{-2}\) has equal constant-term power
sequences, so affine \(C_2\) character weights cancel at every pure row
while a fixed multiplier detects every odd row.  It also checks
\[
 v_{11}\!\left(\mathcal L((y^2+4xy+2x^2)^{12})\right)=3,
\]
where naive Laurent repeated-digit factorization predicts a
valuation-two nonzero residue.  The fourteenth command is
retained as an exact sign regression through 100; all eight all-order
obstruction conclusions are supplied independently by the first three
proof modes.  The accompanying
transverse-lattice proof shows that every order-four return in the 14
reduced types is generated by the primitive order-two and order-three
rows, so there is no separate primitive order-four branch to search.

The bounded counterexample probes which led to the repeated-digit theorem
can be replayed by

```bash
.venv/bin/python scripts/search_binary_gvc_torsion_torus_counterexample.py \
  pair --width 3 --height 1 --depth 14
.venv/bin/python scripts/search_binary_gvc_torsion_torus_counterexample.py \
  shared --width 4 --height 1 --depth 16
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --max-degree 12 --extra-depth 4
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --max-degree 3 --extra-depth 4 --rectangles
python3 scripts/verify_binary_gvc_translation_tangent_rigidity.py
python3 scripts/verify_positive_return_semigroup_jet_rigidity.py \
  --require-singular
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --character-order 3 --max-degree 8 --extra-depth 4
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --character-order 4 --max-degree 7 --extra-depth 4
```

They find no survivor.  Their longest zero prefixes are sparse binomial
lattice-delay examples, first obstructed at rows 10 and 14 respectively.
The third command tests character twists of one translated binomial row.
After quotienting scalar sign, the torus action \(t\mapsto-t\), and
coefficient reversal, it checks 8,188 twists and 81,924 exact
twist--slope rows through depth \(2d+4\), for every \(d\leq12\).  Every
collision is explained by one of those scale-compatible symmetries.
The fourth command adds all \(C_2\) twists on the \((2,2)\) and \((3,2)\)
binomial Taylor rectangles.  Its 2,304 twists and 4,352 exact moving rows
have no collision outside scalar, two-torus, reversal, and
coordinate-exchange symmetry through depths 12 and 14.

The fifth command replays the proved primitive translation-tangent theorem.
For every primitive \((d,r)\) through degree 12 it verifies that the matrix
\[
 \binom dj\binom{d(N-1)}{rN-j}
\]
has the one-dimensional kernel \(j-r\); it also checks the nonprimitive
power/subsequence rank jump, finite-field ranks away from the chosen integer
minors through prime 97, all 2,550 displayed spanning generators of the
universal blind tangent module on 225 rectangular slope cases through
bidegree \((6,6)\), the exact two-dimensional quadratic-Hessian kernel on
100 slope cases through bidegree \((5,5)\), and the factorially weighted
two-free-translate counterexample to module-only inheritance.  The theorems
and the two-direction no-go are unbounded and are proved in
[`BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md`](extended-geometry/BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md).
It implies flatness of every \(q^a\)-order character collision at primitive
one-direction slope once the underlying prime \(q\) is sufficiently large.
The bilinear ghost proves that two-dimensional separation cannot follow from
the first cyclotomic neighbourhood; its first potentially effective row is
quadratic or higher.  The positive return-word proof then shows that this
quadratic row separates every fixed nonflat integer label and does so modulo
all but finitely many primes.

The sixth command replays the broader
[positive return-semigroup jet theorem](extended-geometry/POSITIVE_RETURN_SEMIGROUP_JET_RIGIDITY.md).
Its dependency-free part checks the positive-return group-completion
mechanism on 27 bounded Cartesian configurations.  Singular then constructs
the centered \((2,2)\) derivative ideals from their return words: their
affine dimensions through jet orders one, two, three, and four are
\(4,2,1,0\), and the final quotient has dimension \(40\) with an
eleven-element standard basis.  This is an exact replay of the finite
four-jet certificate and was replayed with Singular 4.3.2.  The
all-configuration full-jet theorem and its
all-torsion-order corollary are proved in the note, not inferred from
the bounded computation.  The same note proves that independently marked
return polynomials at finitely many total degrees separate arbitrary paired
coefficient points modulo the coefficient torus.

The final two commands extend the earlier sign search using exact arithmetic
in \(\mathbb Z[\zeta_3]\) and \(\mathbb Z[i]\).  They test 63,972 and
123,792 moving rows and find no unexplained collision.  These \(C_3,C_4\)
searches, like the rectangular search, are bounded evidence for whether the
signed Hall shell inherits enough common marks before promotion.

The exact primitive translation-orbit census is

```bash
python3 scripts/research_binary_gvc_translation_observability.py \
  --radial-degree 4 \
  --output artifacts/generated-results/binary_gvc_translation_observability_span4.json
python3 scripts/research_binary_gvc_translation_observability.py \
  --radial-degree 5 \
  --modes operator,polynomial \
  --output artifacts/generated-results/binary_gvc_translation_observability_span5_one_colour.json
```

Normaliz 3.10.2 computes the projected Graver bases and Singular 4.3.2
computes the characteristic-zero derivative-orbit ideals after torus
saturation.  At span four, 48 of the 65 normalized mixed packets are already
factorially obstructed and all 17 survivors have empty one-colour and
independent-translation torus ideals.  Two survive only in the weaker
common-diagonal mode; they are the quartic Veronese identities displayed in
[`BINARY_GVC_PRIMITIVE_TRANSLATION_OBSERVABILITY.md`](extended-geometry/BINARY_GVC_PRIMITIVE_TRANSLATION_OBSERVABILITY.md).
The span-five one-colour run has 404 mixed packets, 125 factorial survivors,
and no torus survivor.  The written translation-degree proof gives the
zero-survivor conclusion for every projected span; these commands are exact
bounded regressions, not an extrapolated GVC theorem.  Both outputs also
replay the three-state span-two discriminant orbit identity and certify its
nonzero scale-two factorial obstruction under both fixed-external and
character-power phase laws, which shows why pairwise primitive separation
alone cannot promote a signed linear shell.

The exact small-shell all-scale-prefix eliminations are

```bash
.venv/bin/python scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --structural-certificate --maximum-wronskian-rank 7 \
  --maximum-affine-slope 5 --maximum-affine-offset 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_structural_certificate.json
python3 scripts/verify_binary_gvc_cobham_carry_obstruction.py \
  --primes 3,5,7,11,13 --maximum-period 512 \
  --output artifacts/generated-results/binary_gvc_cobham_carry_obstruction.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 3 \
  --maximum-operator-count 3 --maximum-polynomial-count 3 \
  --maximum-scale 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span2_counts3.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 3 --state-count 3 \
  --maximum-operator-count 2 --maximum-polynomial-count 2 \
  --maximum-scale 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span3_counts2.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 4 \
  --maximum-operator-count 2 --maximum-polynomial-count 2 \
  --maximum-scale 4 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span2_four_state_counts2.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 4 \
  --maximum-operator-count 3 --maximum-polynomial-count 3 \
  --maximum-scale 3 --require-factorial-pairing \
  --output artifacts/generated-results/binary_gvc_all_scale_factorial_pair_circuits_span2_counts3.json
```

The first command is the structural replay.  It checks the exact confluent
Vandermonde leading symbol for all 44 block partitions through rank seven,
uses Singular to prove that the coefficient-torus smooth-conic divisor branch
of the arbitrary four-state scale-1/2/3 ideal is empty, and verifies the
rank-two pair-block converse.  It also checks exact coefficient-span ranks in
the full-span and projectively constant common-base calibrations.  Normaliz
computes the 106 Graver moves of the 35-atom integer-affine factorial
universe with slopes at most five and offsets from -3 through 3; all are
certified sums of exact same-rational-boundary transfers.  The canonical note
also replays the residue-wise rational-minor criterion for an eventually
periodic additive factorial law and detects a one-entry perturbation.  The
canonical note proves the Wronskian identity at every rank, the common-base
differential-ideal collapse, the singular-conic classification, the
boundary-transfer presentation for every integer-affine factorial profile,
and the periodic additive reduction.  Cobham's theorem makes the same
reduction for one finite state sequence automatic in two multiplicatively
independent bases.

The second command proves why complete Hall carry states do not supply that
common sequence.  It checks the Kummer and two-state digit formulas for the
central-binomial carry indicator through \(N=10{,}000\), constructs
equal-residue/opposite-output witnesses for all proposed periods through
\(512\), and checks that every sparse ray \(q p^e\) is stationary.  The
arbitrary-period witness in the canonical note proves nonperiodicity for
every odd prime.  Hence a common two-base automatic refinement retaining
the carry state is impossible by Cobham; this is a proof-route obstruction,
not a GVC counterexample.
Consequently every fixed finite affine-ray family
\(h_j(t)f_j(t)^N\) splits by proportional bases for all large \(N\), with
arbitrary scale-dependent scalar coefficients and changing active support.
The residual identity inside one proportional-base correction space is not
split by this theorem.  The computation is an exact regression for the stated
cross-base proof, not the unrestricted GVC(2) certificate; that certificate
is the written Hall-envelope theorem.

For three states, the note proves at every span that scales one and two force
all orbit-function ratios to be constant; every transferring shell is
therefore torus-empty.  The commands audit the residual zero-transfer block.
The span-two/count-three run tests 11,988 signed triples and kills its 416
scale-one survivors by scale three.  The span-three/count-two run kills all
240 survivors among 8,408 candidates at scale two.  The four-state pilot
tests 928 signed quartets and kills all 40 scale-one survivors at scale two.
These are characteristic-zero Singular saturations.  The counts audit the
proved three-state theorem and the arbitrary-coefficient four-state
scale-1/2/3 theorem; they are not those proofs.  Unrestricted GVC(2) is
proved separately by Hall-envelope separation.

The final command uses the proved all-scale factorial-ray splitting to retain
the 2,882 fixed-sign quartet rows which admit an opposite-sign equal-factorial
pairing.  It finds 142 scale-one survivors.  Character-power scale two kills
all of them.  Fixed signs leave 60 exact all-scale pair cancellations; 12
survive the other one-colour tower and six survive both independent towers.
For those six the script checks both ideal containments and proves that the
saturated pair ideal equals
`(R1^2-4*R0*R2, B1^2-4*B0*B2)`.  This is the already-safe product Veronese
block.  To replay the larger unfiltered finite-prefix census of all 52,416
signed quartets (about ten minutes on the recorded machine), omit
`--require-factorial-pairing` and `--output`; it has 928 scale-one survivors,
868 fixed-sign scale-two obstructions, 60 fixed-sign scale-three survivors,
and no character-power scale-two survivor.

The fast final regression suite is

```bash
.venv/bin/python scripts/verify_binary_gvc_all.py
```

It runs the uniform Hall/weighted-face checker, the regular-trace checker,
the repeated-digit/Newton checker, the weighted affine/factorial obstruction
checker, the default degree-ten translation-twist search, and the structural
Wronskian/four-state certificate.  It is a historical-route regression, not
the written proof of unrestricted GVC(2).  The module-only version of the
former final lemma is
disproved, while distinct bases in every fixed finite affine-ray template are
now separated, and every common-base all-scale ideal has collapsed to a
finite coefficient span.  The exact hypotheses left inside the parked route
are projective carry-rank/safe-rank-drop classification of the bounded
correction circuits and uniform extraction of such a circuit from the growing
positive-density Cartesian face.

## Factorial trace independence

```bash
make verify-factorial-trace-independence
```

To classify two rays directly, encode gamma factors as
`slope:offset[:multiplicity]`:

```bash
python3 scripts/verify_factorial_trace_independence.py \
  --compare '2:1' '1:1/2,1:1'
python3 scripts/verify_factorial_trace_independence.py \
  --compare '2:1' '1:1:2'
```

The first comparison certifies Gauss duplication and reports exponential
base `4`; the second reports the separating signature `[1/2]-[0]`.

The optional independent symbolic replay is

```bash
make verify-factorial-trace-independence-sympy
```

It uses the repository SymPy environment to simplify 322 Gauss/shift and
1,384 census certificates directly as rational-function identities.  The
default target remains dependency-free.

The dependency-free checker reconstructs all 78,124 nonzero signed slope
vectors on slopes `1,...,7` with coefficients in `[-2,2]`.  It integrates
8,134 exact zero-sum translation-orbit divisors and classifies 67,524 products
of at most three rational-offset gamma atoms into 66,140 canonical
signatures, certifying all 1,384 collisions.  It also verifies 161 Gauss
refinements before and after integer shifts, 1,000 seeded signed
transformation cases, and 24 integer-affine reductions.  Its exact
successor-divisor census classifies 82,250 products of at most four
integer-offset atoms into 72,383 classes and decomposes all 9,867 collisions
into boundary transfers.  It also replays the one-scale and entropy
collisions, `m`-fold periodic/rational-slope symmetries, and 2,187 signed
Frobenius-dilation valuation profiles at `p=2,3,5`, and separates all 276 SIC
radial-moment families `(d,r)` with `d<=48`.  These are finite regressions for
the formulas; the complete characteristic-zero gamma-affine classifier,
integer-affine boundary presentation, and exact characteristic-`p` valuation
obstruction are proved in
`extended-geometry/FACTORIAL_TRACE_INDEPENDENCE.md`.
The same proof shows that a factorial ratio with finite value set is
constant, so finite nonzero carry, automaton, sign, or torsion alphabets
introduce no additional projective class.

## Binary GVC prime-power tomography

Requires Normaliz; the pinned run used Normaliz 3.10.2.

```bash
.venv/bin/python scripts/research_binary_gvc_prime_power_tomography.py \
  --output artifacts/generated-results/binary_gvc_prime_power_tomography.json.gz \
  --summary-output artifacts/generated-results/binary_gvc_prime_power_tomography_summary.json
```

The Lawrence lifting gives the exact two-colour Graver basis on levels
`0,...,6`: 8,559 raw relations, 1,584 after translation/dilation and finite
symmetries, and 1,490 mixed support-at-least-five packet candidates on 868
projected supports.  The probes use primes `2,3,5,7,11`, exponents
`m=q*p^e` with `1<=e<=2` and `1<=q<=3`, factorial units modulo `p^2`, and
the `C2,C3` marked-character traces.  Exactly two symmetry classes survive
the configured signature.  They represent

```text
R6*B2*B3 = R0*B5*B6
R6*B1*B2 = R0*B4*B5
```

and together are the reversal orbits of
`R6*B_a*B_(a+1)=R0*B_(a+3)*B_(a+4)`, `a=0,1,2`.  Each exact projected
fibre has two states and no support-at-most-four path.  Equality of the
partitions proves that all factorial, digit, and carry data agree at every
scale, not just in the finite prime window; `C2,C3` are also identically
blind.  A `C4` character separates both classes.  This is a projected
collision census, not a GVC(2) counterexample; the fixed-character Franel
theorem below subsequently closes the family after packet exposure.

The logical JSON result SHA-256 is
`685f60b5843bca33d32034a16a8b599dcfbf46c35e80e62f747f6d9715e285eb`.
The whole-file SHA-256 hashes of the compressed full artifact and compact
summary are, respectively,
`b1cf105d161b83a5e0c23dab0ed3b2cbc0a2726970e20db760ca5ae78eb5c09b`
and
`c9380d66bdad08cb30896c0c1b31d9c5211397d85e962b48d808d90baa832522`.
The precise model and caveats are in
[`BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md`](extended-geometry/BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md).

The strengthened primitive-relation census is:

```bash
.venv/bin/python scripts/research_binary_gvc_prime_power_tomography.py \
  --radial-degree 7 --primitive-only \
  --primes 5,7,11,13 --max-exponent 3 --max-quotient 3 \
  --unit-power 3 --torsion-orders 2,3 \
  --output artifacts/generated-results/binary_gvc_adelic_tomography_span7.json.gz \
  --summary-output artifacts/generated-results/binary_gvc_adelic_tomography_span7_summary.json
```

It computes 34,890 raw Graver relations, 6,601 normalized relations, 6,401
mixed primitive candidates, and the exact bases of 3,107 represented support
semigroups.  Of the candidates, 4,750 are separated by a configured total
valuation.  The remaining 1,651 are exactly the all-scale scalar factorial
collisions, so there are no accidental finite-window collisions and no row
first separated only by a unit residue.  Marked digit/carry data separate
221, `C2,C3` separate another 1,427, and three normalized relations remain.
The first valuation separators occur at `p=5` for 4,557 relations, at `p=7`
for 148 which survive every configured `p=5` probe, and at `p=11` for 45
which survive both earlier primes.  No first separator needs `p=13` or
`e>1` in this span.
They are precisely the span-seven orbits of

```text
R_(s+6)*B_a*B_(a+1) = R_s*B_(a+3)*B_(a+4).
```

Every member is an exact two-state primitive collision, and `C4` separates
it.  Among all 1,430 equal marked-partition relations, the first character
separator distribution is `C2:1244`, `C3:183`, `C4:3`.

The span-seven logical result SHA-256 is
`a8128639805f9e6e0047dc39e70b20f8e939b6a76213930ab44bd0b26a35dde3`.
The whole-file SHA-256 hashes of its compressed full artifact and summary
are
`dc64e57cac395b4f98cfcb8cc0ac1cdf03c598c5fba67253f67c18e379f7f035`
and
`b268bf8ed6cb7c564e154efdcd12b7e60659533fc8cd85e134b6c82870ccabd9`.
This is an exact bounded result in the projected two-colour model, not a Hall
promotion theorem or a counterexample.  It is not the unrestricted proof;
that is the Hall-envelope theorem.

The surviving family has an exact fixed-character termination theorem.
Replay its fibre and coefficient identities with:

```bash
python3 scripts/verify_binary_gvc_six_step_packet_termination.py
```

After row normalization the support is `R6,R0,B0,B1,B3,B4`.  At scale `N`,
the complete `C2,C3`-blind fibre is
`(N-t,t,N-t,N-t,t,t)`, `0<=t<=N`, and its normalized coefficient is

```text
binom(2N,N) * sum_t binom(N,t)^3 * U^(N-t) * V^t.
```

If a further fixed finite character has relative order `h`, the
endpoint-containing rows at scales `h` and `2h` are
`U^h+V^h` and
`U^(2h)+binom(2h,h)^3*U^h*V^h+V^(2h)`.  Their only common zero in
characteristic zero is `U=V=0`, because substitution from the first row
leaves `(2-binom(2h,h)^3)*U^(2h)=0`.  Thus every fixed finite-character
promotion of the six-step family is separated, terminal, or loses support.
The proof is general; the dependency-free script enumerates blind fibres
through scale 12 and checks character orders through 32 as regressions.
Prime-dependent affine-carry promotion to one fixed packet remains unproved
inside the parked route, but is no longer required for binary GVC.

## Binary GVC nonfree-factorization tomography

The all-span consecutive-residue theorem has a quick exact regression:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --verify-consecutive-residues
```

The written incidence-forest proof is uniform in the cyclic order.  The
regression checks orders \(2,\ldots,16\): the \(C_q,C_{q+1}\) histograms are
injective through span \(2q-1\), while at span \(2q\) their kernel is

```text
(1,...,1,0,-1,...,-1)
```

and decomposes into the \(q\) safe beta swaps
\(R_iB_{i+q+1}=R_{i+q+1}B_i\).  This proves that the fixed marked nonfree
factorization quotient is injective at every span; it does not prove that a
prime-dependent Hall shell inherits the required fixed markings.

Run the complete span-four Hilbert/factorial/Graver census with:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py
```

It computes 426 raw and 90 normalized projected Graver relations, giving 65
nonfree profiles.  The exact factorial map is injective on 52; one collision
lattice is reversal-only and 12 contain nonreversal same-vector relations.
For the 11 collision profiles with at most 20 atoms, the complete
factorial-compatible Graver basis has 308 primitive moves.  Packet partitions
separate 15, \(C_2\) separates 207, and \(C_3\) separates 86.  Atom labels
separate every profile, including the two larger deferred Graver lattices.
The first factorial-only square has
\(R=\{0,4\}\), \(B=\{0,1,2,4\}\), counts \((1,3)\), and radial vector
\((8,8)\); \(C_3\) separates it.

Run the larger atom-signature censuses with:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --radial-degree 5 --signature-only --torsion-orders 2,3,4 \
  --output artifacts/generated-results/binary_gvc_nonfree_factorization_span5_signature.json

.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --radial-degree 6 --signature-only --torsion-orders 2,3,4 \
  --output artifacts/generated-results/binary_gvc_nonfree_factorization_span6_signature.json
```

The span-five run has 400 profiles and first-injective distribution
\(169,80,143,8\) at partitions, \(C_2,C_3,C_4\); the span-six run has
1,469 profiles and distribution \(382,358,599,130\).  Neither has an
unresolved atom signature.  The logical result hashes for spans four, five,
and six are, respectively,
`97bcdb8049b34ef0fed2bd0c9f70102e7b06d828e57a27cb09d6036395402627`,
`567712c213029dc01ca749888d21fc838b49bf21604e0a92b73997575d2dc8fd`,
and
`b8f334669719ac829b70182c4c648b8a66a92f409dd65b8945c6309ea2a6ecde`.
Their whole-file hashes are
`888f1f465c67045a2a39a157d6d6cf4872f1cb1a76a79b248e31aa658fd21d2d`,
`991f531db27e707a2155080f84751f2b3c12f9f88c957775cf6e62f26776a18d`,
and
`052928f81694395d1369a4e2c1e9973ebd1fc4d28d10e9c2216f549d4c1b1e99`.
These are exact bounded projected-semigroup computations; the all-span
claim comes from the incidence proof, not extrapolation from the census.

The accompanying exact finite-moment search is

```bash
.venv/bin/python scripts/search_binary_gvc_five_channel_descent.py
.venv/bin/python scripts/search_binary_gvc_five_channel_descent.py --frontier-suite
.venv/bin/python scripts/search_binary_gvc_five_channel_descent.py \
  --frontier-suite \
  --json-output artifacts/generated-results/binary_gvc_five_channel_pivot_clusters.json
```

The first command runs the quicker \((1,2)\) and \((1,3)\) cases.  The
frontier suite adds \((1,4)\) and \((2,3)\).  It enumerates respectively
13,288, 41,728, 95,368, and 253,576 supports with zero through three added
channels in the order-four balance boxes.  Successive saturation finds
the first unit at moments two, three, and four with distributions
\[
\begin{array}{c|rrr}
(r,s)&M_2&M_3&M_4\\ \hline
(1,2)&13082&173&33\\
(1,3)&41299&394&35\\
(1,4)&94636&679&53\\
(2,3)&252442&1074&60
\end{array}
\]
and no torus survivor.  Points outside those boxes cannot occur in a
balanced return through order four, so the computations cover arbitrary
nonnegative support with at most five channels for all four endpoint
pairs.  They are exact rational finite-moment computations, not a proof
for arbitrary endpoint orders or three operator endpoints, and not the
Hall-envelope proof of unrestricted GVC(2).  The third command additionally records every fourth-pivot
support and canonicalizes its balanced-selection rows under endpoint
exchange and permutation of the three added channels.  The 181 supports
collapse to 14 return-matrix types, all already realized at `(1,2)`.  The
generated JSON has SHA-256
`59436a3617671c4ca47cd354b45cb74abc7b9787352e725c97ebce1a304ffa16`.

The
proof and the precise remaining restricted beta--torus
coupled-convolution obstruction are in
[`BINARY_GVC_UNIFORM_FACE_TERMINATION.md`](extended-geometry/BINARY_GVC_UNIFORM_FACE_TERMINATION.md).
The checker is a regression, not the proof.

The explicit homogeneous GVC(3) counterexample and its exact consequences
are replayed by

```bash
python3 scripts/verify_gvc3_homogeneous_counterexample.py
python3 scripts/verify_gvc3_homogeneous_spillovers.py
.venv/bin/python scripts/verify_gvc3_power_tail_and_minimum.py
.venv/bin/python scripts/verify_gvc3_independent_parity_quartic.py
.venv/bin/python scripts/verify_gvc3_isotropic_harmonic_channels.py
.venv/bin/python scripts/verify_gvc3_four_coherent_channels.py
.venv/bin/python \
  scripts/research_gvc3_degree10_distinct_four_channels.py \
  --workers 11 --timeout 900
.venv/bin/python \
  scripts/verify_gvc3_degree10_four_channel_collisions.py \
  --workers 8 --modular-timeout 180 --exact-timeout 300
.venv/bin/python \
  scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_compile5.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 4 --compile-only \
  --scan-linear-pivots --scan-pivot-max-order 4 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_pivot_scan4.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --linear-pivot a4 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_pivot_compile5.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --linear-pivot-boundary a4 --boundary-linear-pivot a2 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a2_pivot_compile5.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a3_pivot_compile5.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 3 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method quotient --report-quotient-dimension \
  --primes 101 103 107 --timeout 60 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a3_pivot_modular3.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 3 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method successive --report-quotient-dimension \
  --primes 101 103 107 --timeout 60 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a3_pivot_successive3.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 4 \
  --normalize-coefficient a1 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method successive --report-quotient-dimension \
  --primes 101 --timeout 5 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a1norm_a4_boundary_a3_pivot_successive4_p101.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 4 \
  --normalize-coefficient a1 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method successive --report-quotient-dimension \
  --primes 101 103 107 --timeout 90 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a1norm_a4_boundary_a3_pivot_component_cycle_modular3.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 3 \
  --linear-pivot-boundary a4 --saturation-method quotient \
  --primes 101 103 107 --timeout 60 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_modular3.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 \
  --linear-pivot-boundary a4 --saturation-method quotient \
  --primes 101 --timeout 300 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_modular5.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 4 \
  --normalize-coefficient a1 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method successive --report-quotient-dimension \
  --primes 101 --timeout 300 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a1norm_a4_boundary_a3_pivot_hbranch_p101.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 \
  --normalize-coefficient a1 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --coordinate-boundary-only --exact-coordinate-boundary-support \
  --primes 101 --timeout 600 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a1norm_a4_boundary_a3_pivot_exact_coordinate_collision5_p101.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 6 \
  --normalize-coefficient a1 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --compile-only --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a1norm_a4_boundary_a3_pivot_compile6.json
.venv/bin/python \
  scripts/research_gvc3_degree10_five_channel_slice.py \
  --lam 2 --mu 3 --max-order 8 \
  --modular-timeout 180 --exact-timeout 300
.venv/bin/python scripts/verify_gvc3_cusp_profile_suspension.py
.venv/bin/python scripts/research_gvc3_harmonic_cubic_profile.py \
  --cases alpha1_n3:8 alpha1_d_k:7 alpha1_dk:7 \
          alpha0:7 alpha0_radical:7 \
  --primes 101 103 107 --timeout 900 \
  --exact-all --msolve-threads 4
```

The first checker verifies polynomiality, homogeneity, primitivity, the
closed-form detector, and two exact finite replays through moment six.  The
all-order counterexample is proved in
[`THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md`](extended-geometry/THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md),
not inferred from the bounded replay.  The second checker replays the
homogeneous dimension and quadratic-rank spillovers.  The third checks the
maximal shifted-power formulas and the scoped one-profile minimum.  The
fourth verifies the exact characteristic-zero elimination for the
independent-linear quartic repair.  The fifth checks the full
two/three-channel isotropic coherent-state obstruction in balanced degrees
four through ten: it compiles invariant Reynolds moments, uses three modular
Gröbner runs to discover every chart cutoff, repeats all forty saturated
chart eliminations exactly over \(\mathbb Q\), and rejects the
\(\mathcal H_2\oplus\mathcal H_4\oplus\mathcal H_6\) near survivor at pure
moment five despite its surviving multiplier channel.  The sixth checks the
unique four-channel degree-eight coherent family.  It records the finite
proper-hypergeometric occupation/edge sum, proves all thirteen direction
collisions and the four-distinct `B=0` boundary exactly over \(\mathbb Q\),
and retains the `B!=0` unit bases at \(101,103,107\) as unpromoted modular
evidence.  The optional `--exact-open` flag attempts that remaining
characteristic-zero chart and fails rather than promoting a timeout.  The
next command checks all five pairwise-distinct four-of-five coherent profiles
in balanced degree ten.  Quotient saturation gives the unit ideal through
moment seven at \(101,103,107\) for every profile.  This is an exact statement
over the declared finite fields but remains unpromoted modular discovery;
direction collisions and the five-channel profile are not included.  The
following checker enumerates the thirteen nonterminal direction partitions
for each of the four new profiles containing \(\mathcal H_{10}\).  It uses
quotient saturation at (101,103,107) for cutoff discovery and replay, then
proves all 52 localized charts empty over \(\mathbb Q\) by literal msolve
bases \([1]\).  The exact cutoffs have distribution (21,3,5,22,1\) at
moments three, four, five, six, and nine.  This command requires Singular and
msolve.  The next compile-only command records the five-channel term growth
\(271,1142,3686\) through moments three, four, and five; it performs no
elimination.
The next two compile-only commands compare every linear \(\mu_3\) pivot and
record the selected \(a_4\) localization.  Their transformed term counts are
7,333 and 58,971 at moments four and five; they make no elimination claim.
On the exceptional \(A=B=0\) boundary, the next compile-only command solves
the irreducible four-channel equation \(A=0\) for each eligible nested pivot.
The \(a_3\) chart wins: its transformed \((B,\mu_4,\mu_5)\) counts are
648, 7,340, and 15,300, versus 958, 6,173, and 34,782 for \(a_2\).  The following
three-prime command records 60-second timeouts already on transformed \(B\)
with the complete localization product.  It makes no survivor claim and
motivates successive factor saturation.  The next command performs that
factor-by-factor replay at cutoff three: all ten irreducible factors complete
at 101, 103, and 107, showing that the earlier timeout was algorithmic.  The
\(a_1=1\), cutoff-four command is a bounded component diagnostic at 101.  It
certifies the \(a_0\)-colon at exponent two by two residual-intersection
minors and clears six further factors by unmixed boundary gcds.  The two
remaining linear boundaries and the final degree-17 factor form one coupled
component cycle; the artifact is explicitly partial and makes neither a unit
nor a characteristic-zero claim.  The next command solves the final factor
linearly for (a_0) and replays the reciprocal boundary gcd over 101, 103,
and 107.  All five cleared-restriction identities have zero remainder; after
the seven completed factors are removed, the degree-29 gcd has valuations
one on (143a_0+60a_2), twelve on (lambda-mu), and no residual factor.
Thus `factors_component_classified=10` and `support_closed=1`, while
`scheme_closed=0`: this is a support classification, not an explicit
generic saturated ideal.  The same artifact contains an exact generic normal
certificate over `Q(a2,a4,lam)`: the original two cleared restrictions have
collision valuations 12 and 20, their exact 728-term gcd has only the cleared
`a0`, `143*a0+60*a2`, and `(lam-mu)^12` factors, and a nonzero 959-term
normal-jet determinant proves generic local length 12.  Its already inverted
chart factors occur with valuations three, two, and two; the residual
exceptional factor has 115 terms and total degree 47.  On the deterministic
rational normal fiber
\((a_2,a_4,\lambda)=(2,3,5)\), the same command promotes the local calculation
over \(\mathbb Q\): the literal basis is `(dd,ee^12)` and the quotient length
is 12.  The dedicated `hbranch` replay additionally proves that the 115-term
exceptional factor is irreducible and squarefree, has generic normal length
18, and has no deeper stratum on the declared chart: its next degree-57
factor meets it only on the already inverted companion boundary.  The
cutoff-five exact coordinate-boundary command compiles the 38,015-term third
restriction.  Its three-generator gcd has 56 terms and is exactly the
product of the cleared `a0` and `143*a0+60*a2` factors.  After all chart
factors are removed, the third restriction on `mu=lam` is
`a2^3*lam^8*(lam-1)^36`; together with the second restriction's
`(lam-mu)^20` factor, this proves that the localized `D=0` fibre is empty.
The cutoff-six compile-only replay records original moment-six term count
9,263 and nested-pivot term count 59,418.  It makes no elimination claim.
The earlier full coordinate colon, component quotient, and larger-chart
quotient timeouts remain algorithmic diagnostics and make no survivor claim.
The following two commands retain the corresponding pre-nested 60- and
300-second time bounds; they also make no survivor claim.
The rational-slice checker then specializes the two cross-ratios to
\((2,3)\) before moment compilation.  On the exceptional \(a_4\)-pivot
boundary \(A=B=0\), quotient saturation has dimensions (2,1,0) through
moments three, four, and five and quotient length 36 at moment five; moment
six is the unit ideal at (101,103,107), and exact msolve elimination over
\(\mathbb Q\) returns \([1]\).  This is a theorem for that boundary slice,
not the generic five-channel chart.  The following command checks the full winding--profile--radial suspension,
its cusp identity,
complete phase ladder, top Reynolds--apolar contractions, exact trace depths,
and direct
shifted-power detectors for the non-power profile \(S=1+z\).  The final
command is an exact search with modular discovery replay in the complete
harmonic-cubic repair:
it compiles seven invariant weight channels, covers the nonzero-even-part
chart by three pivot strata, and audits the radical of the zero-even-part
boundary.  It requires Singular and msolve.  All nine projective and
boundary saturations are exact over \(\mathbb Q\); the three-prime runs are
retained as discovery replay.
These commands reproduce the corresponding `gvc3_*.json` artifacts.  The
degree-eight four-channel artifact deliberately has mixed exact/modular
status; the pairwise-distinct degree-ten artifact is modular only, its
direction-collision companion is exact over \(\mathbb Q\), and the
initial five-channel and pivot artifacts are compile-only.  The five-channel
rational-slice companion is exact over \(\mathbb Q\).  The final command writes the
separate exact harmonic-cubic
artifact in `artifacts/generated-results/`.  These results
disprove GVC in every dimension at least three but do not disprove the
ordinary-Laplacian/Hessian-nilpotent conjecture.

The integration replay on 2026-08-02 used system Python 3.12.3 for the first
two dependency-free commands and the locked `.venv` Python 3.13.5 for the
SymPy commands.  Regenerating the first two artifacts normalized JSON list
formatting without changing their parsed content.  The counterexample
artifact changed from
`sha256:f05e5bee5c9b9aab5e245026f99af30a4f379bd88b9f3163fbdd7859f56aba06`
to
`sha256:ef44b4d7390ca261c432d23bcfc7b262062d3027b4f46ca9ccef1b9c556ec04d`;
the spillover artifact changed from
`sha256:ef27d18337a34a527140a799e15a3a242b008bc1b6c22e403ea775d86df31b50`
to
`sha256:6810131f43f822c39d4abed682e97570254faed6b1e19080cba1559821f2a666`.
`jq` comparison of each old/new pair is exactly `true`.

The earlier exact three-variable tagged-lift reduction and its bounded
extensions are replayed by

```bash
.venv/bin/python scripts/research_three_variable_gvc_tagged_lift.py
```

The target first replays the coordinate-only detector extracted from the
two-pair witness through order six and uses the proved equal-channel
degree-tag identity.  Over
\(\mathbb Q\), it then computes the complete binary-cubic operator-jet
moment ideals for the literal Long tag: moments one through four remain
nonunit and moment five gives the unit ideal.  The same run performs three
explicitly experimental calculations over \(\mathbf F_{101}\): the
canonical rank-five auxiliary chart, the normalized factor-compatible
cubic-profile chart, and 200 deterministic general cubic-profile fibers.
Their output is written to
[`three_variable_gvc_tagged_lift.json`](artifacts/generated-results/three_variable_gvc_tagged_lift.json).
The exact formulas, scope distinctions, and remaining mixed-order target
are documented in
[`THREE_VARIABLE_GVC_TAGGED_LIFT.md`](extended-geometry/THREE_VARIABLE_GVC_TAGGED_LIFT.md).
These are retained architecture exclusions; the homogeneous counterexample
above now settles GVC(3) negatively.

The coupled order/degree-\((2,3,4)\) continuation is replayed by

```bash
.venv/bin/python scripts/research_three_channel_gvc_lift.py
```

It enumerates all \(56\) oriented rank-three parallelograms on the positive
weighted-quartic plane, computes their exact rational moment ideals, and
then checks the persistent five-term radical, the complete quartic repair
on the polynomial side, the sparse activated operator endpoints, and the
complete odd-quartic operator/polynomial jet.  The remaining radicals have
written all-order mixed cutoffs.  The generated
[`three_channel_gvc_lift.json`](artifacts/generated-results/three_channel_gvc_lift.json)
records these characteristic-zero results.  The odd chart closes at
moment six with radical \((A,S,RU)\); the simultaneous complete even-and-odd
quartic total space remains open.

The first repeated-root continuation is the migrating defect-one ansatz
\((\Lambda_4+\Lambda_5,P_5+P_4)\).  Run its faithful-characteristic samples,
followed by the conditioned defect-two \((\Lambda_6,P_3)\) search, with:

```bash
.venv/bin/python \
  scripts/search_binary_repeated_quartic_gvc_jets_mod_p.py
```

The pinned run takes roughly two minutes.  Use `--quick` for a small
non-pinning regression.  The generated
[`binary_repeated_quartic_gvc_jet_search.json`](artifacts/generated-results/binary_repeated_quartic_gvc_jet_search.json)
is bounded experimental evidence only.  Its scope, the support-separator
proof for the two fifth-moment survivors, and the unique defect-two
fourth-moment failure are documented in
[`BINARY_REPEATED_QUARTIC_GVC_JET_SEARCH.md`](extended-geometry/BINARY_REPEATED_QUARTIC_GVC_JET_SEARCH.md).

The first rank-efficient ordinary-Laplacian lift is excluded by

```bash
.venv/bin/python scripts/verify_dvorsky_one_pair_schur_obstruction.py
```

Pairing \(t\) with one new variable \(s\) makes
\[
 \partial_a\partial_d-\partial_b\partial_c+\partial_t\partial_s
\]
a nondegenerate quadratic operator in six variables.  The checker first
retains the homogeneous cubic regression, then parametrizes the
unrestricted transverse two-jet of an arbitrary polynomial or formal
harmonic lift.  It proves the exact axis identity
\[
 \widetilde\Delta^2(F^2)=12t^2-8\rho t
 \quad\text{modulo }(a,b,c,d,s).
\]
Thus no degree mixture can repair the canonical six-variable hyperplane
lift.  Different quadratic completions, additional blocks, and nonlinear
specializations remain open.  See
[`DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md`](extended-geometry/DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md).

The still-open bidegree-\((3,3)\) classification is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_sextic_slice.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_anchor_jacobians.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_boundary_family.py
```

The checker eliminates the full seven-parameter one-sided locus in the
sixteen-dimensional bidegree-\((3,3)\) coefficient space, giving its exact
dimension seven, projective degree twenty, and Gröbner-basis size \(148\).
It also evaluates the full moment Jacobian without expanding the
sixteen-variable moments and verifies a displayed nonzero
\(13\times13\) integer minor for \(\mu_1,\ldots,\mu_{13}\). Thus these
moments are algebraically independent and attain the invariant-quotient
dimension bound. It then computes the invariant Hilbert series from the
sixteen exact \(\mathrm{SL}_2\)-weights: the proposed numerator for
degrees \(1,\ldots,13\) has coefficient \(-2186\) in degree \(63\), so
those moments cannot be a homogeneous system of parameters and their zero
fiber has an extra semistable component. Replacing \(\mu_{13}\) by
\(\mu_{14}\) gives the least-total-degree corrected candidate; the checker
verifies another exact rank-thirteen minor and a nonnegative proposed
Hilbert numerator through degree \(100\), while making no zero-fiber
claim. It also inverts the full Clebsch--Gordan basis to recover the global
quadratic discriminant and proves that moments \(1,\ldots,4\) have only
the origin on the maximal-torus fixed diagonal slice, with four
seventh-power certificates.
It then constructs the highest \(\operatorname{Sym}^6\) summand by the
\(\mathfrak{sl}_2\) lowering chain.  On this binary-sextic slice, exact
elimination and power reductions prove that moments \(2,4,6,10\) have the
same radical as the \(L^4Q\) nullcone; the ten nullcone generators have
power certificate \((1,5,5,5,5,5,5,5,5,5)\).  Thus this slice contains no
SIC(2) counterexample.  The same checker proves direct ideal equality with
the \(L^3R\) quartic nullcone from moments \(2,3\), and with the \(L^2\)
quadratic nullcone from moment \(2\).  Thus all three pure irreducible
summands are closed.  On the normalized non-null quadratic branch, moments
through order six prove \(c^6\) lies in the
\(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\) moment ideal over
\(\mathbb Q\), excluding that branch.  For
\(\operatorname{Sym}^6\oplus\operatorname{Sym}^2\), the checker records
only the finite-field result over \(\mathbb F_{32003}\): even moments
through order fourteen give a basis of size \(7576\) and contain \(c^{25}\)
but not \(c^{24}\).  The exact characteristic-zero lift and the full
mixed-summand problem remain open; see
[`TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md).
The second checker normalizes a non-null quadratic to \(2XT\), covers the
full higher \(\operatorname{Sym}^6\oplus\operatorname{Sym}^4\) locus by
five residual-torus chart orbits, and proves by displayed nonzero exact
eleven-by-eleven Jacobian determinants that moments \(2,\ldots,12\) are
algebraically independent on every chart orbit.  This is a
dimension-sized coordinate theorem, not a zero-fiber exclusion. It also
verifies the exact normalized formula for \(\mu_2\), whose derivatives in
the five opposite-weight chart variables are respectively
\(-72,432,-1080,336,-1344\). Thus one variable is eliminated without
saturation on every chart. On \(s_0=1\), the checker then eliminates
\(s_6\) and proves that \(\mu_3\) is affine in \(s_5,t_4\), recording its
two explicit pivot coefficients and their common boundary. On the natural
two-parameter plane
\(t_0=a,t_3=3a,t_4=b,s_6=(14ab+70)/3\) in that boundary, it verifies
\(\mu_3=1866240a^3\), computes the fourth moment, and checks an explicit
unit certificate modulo \(a^3\). Thus this sparse boundary plane contains
no moment-zero point. Finally it evaluates the full chart Jacobian at
exact rational points in \(A\ne0\), \(A=0,B\ne0\), and \(A=B=0\).
Together with the independent gradients of \(A,B\), the nonzero
determinants prove maximal restricted differential ranks \(11,10,9\).
The checker also verifies the constant triangular pivots
\(\partial A/\partial t_3=1\) and
\(\partial B/\partial s_4=-3\), with both cross derivatives zero.
Thus \(A=0\) eliminates \(t_3\) globally.  The fully substituted
\(A=0,B\ne0\) export has nine effective variables after the \(\mu_3\)
pivot, and the \(A=B=0\) export eliminates \(t_3,s_4\) and has eight
effective variables.
The third checker enlarges the common-boundary plane to the exact
four-parameter family
\[
 s_4=-4q^2,\quad s_5=h,\quad
 (t_0,t_1,t_3,t_4)=(a,q,3a,b),\quad
 s_6=(14ab-168aq+70)/3.
\]
Here \(A=B=\mu_2=0\) identically.  Moments \(3,\ldots,6\) leave a
zero-dimensional quotient of length \(372\), and adjoining moment \(7\)
gives the unit ideal over \(\mathbb Q\).

The reduced common-boundary fiber calculation is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_boundary_fiber.py
```

After the constant \(t_3,s_4,s_6\) substitutions, it records a rational
\(\mu_3=0\) base point at which
\((\mu_4,\mu_5)\) cuts out a length-six quotient in the two fiber
variables \(s_5,t_4\).  Its six standard monomials are
\(1,t_4,t_4^2,t_4^3,s_5,s_5t_4\).  Openness promotes this to an exact
characteristic-zero rank-six theorem on a nonempty open of the
\(\mu_3=0\) base; it is not a zero-fiber exclusion.

The generic boundary quotient and its denominator strata are replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_boundary_generic_quotient.py
```

The checker constructs the sparse moments integrally and works over the
characteristic-zero rational-function field of the six base variables.
It obtains a three-element basis and quotient length six for
\((\mu_4,\mu_5)\).  Both \(\mu_6,\mu_7\) reduce to six fiber coordinates.
The \(\mu_6\) denominators are supported on
\(L=s_1t_0-t_1\) and
\(Q=s_1^2-s_2-(13/3)t_0^2\); \(\mu_7\) introduces the additional
quartic divisor \(J\) recorded in the artifact.  The checker also proves
that the \(t_4^3\)-coefficient of the \(\mu_6\) normal form is
affine-linear in \(t_2\), with constant derivative
\(-100078239744000\), so it globally eliminates \(t_2\) on the
generic \(LQJ\)-open.  It is also affine-linear in \(s_3\), with the
explicit alternate eight-term pivot \(H\).
In the adapted coordinates it verifies
\(J=(99Q+155t_0^2)^2+30420L^2\) and
\(H=32J+1179Q(99Q+155t_0^2)\), giving the exact localized exclusion
\((J,H):(LQ)^\infty=(1)\).
The checker also constructs an irreducible quartic number field and an
explicit point on \(J=\mu_3=0\).  At that point
\((\mu_4,\mu_5)\) has length five, with initial ideal
\((s_5^2,t_4^3,s_5t_4^2)\) and basis
\(1,t_4,t_4^2,s_5,s_5t_4\).  It then proves the same statement at the
generic point of \(J=0\).  Over
\(\mathbb Q(\alpha)(s_1,s_3,t_0,L,t_2)\), \(\alpha^2=-30420\), a
quadratic-pair fraction-free calculation constructs a three-element
Gröbner basis with supports \(6,7,6\) and leading monomials
\(s_5^2,s_5t_4^2,t_4^3\).  One pair is removed by the product criterion
and the final pair reduces exactly to zero in five steps.  Generic
rational-function-field calculations at both split roots modulo 47 and
101 independently reproduce the basis.  The same exact quadratic-pair
reducer sends \(\mu_6,\mu_7\) to normal forms supported on all five
standard monomials in respectively five and ten pseudo-reduction steps.
After solving the constant \(t_2\)-pivot, the checker forms the
642-term cubic \(P(s_3)\) coming from \(\mu_3\) and the cubic
\(s_5t_4\)-coefficient \(C(s_3)\) of the remaining \(\mu_6\) normal
form.  Their exact resultant factors as
\[
 \operatorname{Res}_{s_3}(P,C)=L^6Q^6\mathcal R_{63}.
\]
The residual factor has degree 63 and 6702 rational terms.  Its
degree-preserving reduction modulo 47 is irreducible, proving
\(\mathcal R_{63}\) irreducible over \(\mathbb Q\); reduction modulo
101 independently reproduces the degree and irreducibility.
For the next \(t_4^2\)-coefficient, the checker obtains
\(\operatorname{Res}_{s_3}(P,C_2)=L^9Q^6\mathcal T_{66}\).
Both reductions of \(\mathcal T_{66}\) are irreducible, and the checker
verifies modulo both primes that
\(\gcd(\mathcal R_{63},\mathcal T_{66})=1\).  Thus the principal
\(\mu_6\)-zero base has codimension at least two rather than a surviving
degree-63 divisorial component.
For the two cubics \(P=as_3^3+bs_3^2+cs_3+d\) and
\(C=es_3^3+fs_3^2+gs_3+h\), it also constructs the direct linear
pseudo-remainder \(V_1s_3+V_0\).  Over \(\mathbb Q\), \(V_1,V_0\)
have degrees \(60,63\) and respectively 2105 and 5170 terms.
Their reductions modulo both primes are coprime to
\(\mathcal R_{63}\).  Consequently the degree-63 incidence branch has
the dense rational pivot \(s_3=-V_0/V_1\); the checker does not assert
that this pivot covers every component of the residual codimension-two
intersection.
On the divisor strata \(L=0\), \(Q=0\), and \(L=Q=0\), the exact
characteristic-zero quotient lengths are respectively six, five, and
five.  Their changed standard monomial bases are recorded, and the
normal forms of both \(\mu_6,\mu_7\) occupy every basis coordinate.
Independent calculations over \(\mathbb F_{47}\) and
\(\mathbb F_{101}\) replay all quotient shapes.  This is a finite
quotient certificate, not a full boundary unit certificate.

The content-preserving corrected-boundary continuation is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py \
  --branch LQ --prime 0 --through 10 --export-only \
  --include-branch-table --deepest-ffnf --deepest-ffnf-through 7 \
  --t0-zero-branch-table --t0-open-rank-six --timeout 300 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_corrected_boundary_deepest.json
```

It reconstructs the generic, \(L=0\), \(Q=0\), and \(L=Q=0\)
\((\mu_4,\mu_5)\) algebras over \(\mathbb Q\), imports the separate exact
rank-five \(J=0\) result, and proves that on \(L=Q=0\) the algebra after
\((\mu_3,\mu_4,\mu_5)\) has length \(15\) over
\(\mathbb Q(s_1,s_3,t_0)\).  The checker deliberately does not call
`cleardenom` on these exports: that command removes base content and can
turn a nonconstant base equation into \(1\).  In the rank-fifteen algebra,
it pseudo-reduces \(\mu_6,\mu_7\) exactly without division by a base
polynomial; both primitive normal forms occupy all fifteen standard
monomials.  It also proves
\((\mu_3,\ldots,\mu_{10})=(1)\) over \(\mathbb Q\) on
\(L=Q=t_0=0\).  More generally it saturates the generic, \(L\), \(Q\),
\(J\), and \(L=Q\) strata by their specialized principal opens and returns
the unit ideal through \(\mu_{10}\) on every \(t_0=0\) stratum.  These five
opens partition the adapted \((L,Q)\)-plane, so this closes the entire
branchwise \(t_0=0\) divisor.  On the remaining open it normalizes
\(t_0=1\), sets \(u=s_0^{-1}\), eliminates \(t_4,t_3,s_4\), and verifies
that \(\mu_3\) is fiber-independent while \((\mu_4,\mu_5)\) cuts out an
exact rank-six algebra in \(s_6,s_5\).  The exact \(\mu_6\) normal form
occupies all six standard monomials; the three leading coefficients have
the explicit \(K,H,Q_*KJ_*H\) factorization recorded in the artifact,
with \(K=4A_*-15Q_*\) and \(H=4J_*-15A_*Q_*\).  On \(K=0\), a changed
basis retains length six away from \(\ell J_*=0\); on the reduced
\(K=H=0\) linear slice, another changed basis has leading ideal
\((s_5^2,s_6^3)\), length six, and a six-coordinate \(\mu_6\) normal
form.  The separate exact \(J_*=0\) calculation covers the remaining
\(K=J_*=0\) intersection.  Finally it parametrizes the rational conic
\(H=0\) by (5.12t) and obtains the exact generic leading ideal
\((s_6^2,s_5^3)\) and length six; the omitted parametrization point lies
on \(J_*=0\).
Exact and mod-\(47\) pseudo-reductions of
\(\mu_8\), and the reordered mod-\(47\) full deepest solve, reach the
recorded \(600\)-second bounds.  The \(t_0\)-open common-root equations,
lower-dimensional coefficient specializations, inherited \(Q_*,J_*\)
branch radicals, orders \(7,8,9,10,11,12,14\), and the rational radical
remain open.

The first \(t_0\)-open common-root step has a separate subsecond exact
certificate:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py \
  --branch generic --prime 0 --t0-open-fixed-fiber --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_open_fixed_fiber.json
```

At the rational base (5.12u), it verifies \(\mu_3=0\), nonvanishing of
\(Q_*J_*KH\), length six for
\(\mathbb Q[s_6,s_5]/(\mu_4,\mu_5)\), and the exact unit ideal
\((\mu_4,\mu_5,\mu_6)=(1)\).  This proves that the first norm is not
identically zero on the local base component, but does not expand or
classify its exceptional divisor.

The same fixed point extends to the exact rational \(\mu_3=0\) curve
(5.12v).  Its first norm and the next Fitting coefficient are replayed
by:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py \
  --branch generic --prime 0 --t0-open-curve-norm --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_open_curve_norm.json
```

The checker obtains an irreducible degree-\(198\) numerator for
\(\det M_{\mu_6}\), with denominator degree \(144\).  The coefficient of
\(z\) in \(\det(M_{\mu_6}+zM_{\mu_7})\) has numerator and denominator
degrees \(209,153\), and its numerator is coprime to the degree-\(198\)
norm.  Thus \((\mu_6,\mu_7)\) has no common root on the norm divisor
where the curve and border chart are defined.  The norm denominator
factors into the specialized \(Q_*=0\), curve-pole, and \(J_*=0\)
factors with degrees \(2,3,4\).  Exact degree-two and degree-four
number-field calculations give length five for \((\mu_4,\mu_5)\) and
the unit ideal after adjoining \(\mu_6,\mu_7\) on both \(Q_*=0\) and
\(J_*=0\).  The cubic factor is a genuine pole with coprime numerator,
not an affine point of the parametrized curve.  Thus every defined point
of this rational curve is excluded.

Directional modular interpolation shards for the first two Fitting
coefficients use paired roots of the quadratic \(\mu_3(s_3)\):

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 1019 --variable s1 --sample-count 450 \
  --training-count 400 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_fitting_s1_mod1019.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_t0_fitting_degree_scout.py
```

The generated scout aggregates all five directions through two base
points, at primes \(1019,2039\): 6750 paired samples fit fifteen rational
line reconstructions and 750 unused pairs verify them.  The common
observed denominator models are (5.12x).  This supplies stable degree
bounds only; it does not reconstruct the dense five-variable
numerators.

The same engine can evaluate the complete degree-six determinant pencil
on deterministic random paired-root shards and replay only its common-zero
candidates:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 47 --random-seed 102 --sample-count 450 \
  --max-attempts 2000 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_pencil_random_p47_seed102.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_t0_pencil_random_scout.py
```

The expanded aggregate also includes the smallest admissible prime and
the direct divisor mode:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 404 --sample-count 450 \
  --max-attempts 2000 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_pencil_random_p43_seed404.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 501 --stratum Q --sample-count 450 \
  --max-attempts 2000 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_random_p43_seed501.json
```

The aggregate contains forty-four generic shards at primes
\(43,47,59,71\): 19800 accepted paired bases, 39600 evaluated roots of
\(\mu_3\), and twenty direct length-one common roots through \(\mu_7\).
All twenty have block rank five through \(\mu_7\) and rank six after
adding \(\mu_8\); nineteen use the leading \(M_6\) pivot and one uses a
second pivot.  Four direct divisor scouts evaluate another 3600 roots.
The generic quotient lengths are five on \(Q,J\) and six on \(K,H\);
one reduced through-\(\mu_7\) point occurs on each of \(Q,J\), and both
are excluded by \(\mu_8\).  This is a bounded modular scout, not a
reconstruction or a characteristic-zero exclusion.

The specialization, rank-complement, leading-border, and projected-border
continuation uses only short guarded jobs.  Representative producers and
the aggregate verifier are:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 602 --stratum QJH \
  --sample-count 450 --max-attempts 10000 --retain-pairs --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_QJH_random_p43_seed602.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 901 --sample-count 225 \
  --max-attempts 6000 --pivot-scout --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_rank_complement_random_p43_seed901.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_leading.py \
  --prime 43 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_leading_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_border_resultant.py \
  --prime 43 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_border_resultant_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_leading.py \
  --prime 0 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_leading_exact.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_border_resultant.py \
  --prime 0 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_border_resultant_exact.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_t0_strata_rank_continuation.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_Q_residual_slice.py
```

The verifier pins 31 direct-stratum artifacts.  Across \(Q,J,K,H,KH\),
\(QJH,JH,JK\), \(a_2=0\), and the repeated-root stratum, seventeen
sampled reduced common roots through \(\mu_7\) are all excluded by direct
\(\mu_8\) evaluation.  Twelve rank-complement shards test another 6300
\(\mu_3\)-roots: four miss both selected pivots, but all four have full
rank in \([M_6\ M_7]\), and no joint-rank-at-most-four point occurs.

The `stratum_leading` command uses `liftstd`; the LCM of its leading
coefficients is the specialization border.  Modulo \(43\) it gives the
nine irreducible border profiles in (5.12z) and contains every retained
length-drop point.  On \(Q,J,JK\) its sampled zeros are exactly the
length drops.  The degree-at-most-four \(Q,J\) point-cloud kernels contain
only the ambient cubic equation and its five linear multiples.

The border-resultant command projects the quadratic \(\mu_3(s_3)\)
against this border and factors the result.  The residual factors have
the degree, term-count, and multiplicity profiles in (5.12aa).  The
linear pseudo-remainder \(A s_3+B\) is coprime to every residual factor,
so all have the dense pivot \(s_3=-B/A\).  These are exact finite-field
calculations, not characteristic-zero component certificates.

The two `--prime 0 --stratum Q` commands promote the \(Q\)-row.  Over
\(\mathbb Q\), its leading border is irreducible of degree \(36\) with
588 terms.  The exact degree-\(76\) resultant factors as
\(c\,u^{20}J_Q^4R_{20}^2\), where \(J_Q\) is the inherited four-term
quartic and \(R_{20}\) is irreducible of degree \(20\) with 200 terms.
The exact linear pivot is coprime to \(R_{20}\), and reduction modulo
43 matches the modular artifacts up to units.  The next step is custom
arithmetic in this degree-five extension and the remaining modular
degree-five/degree-six extensions; direct Gröbner recomputation is
deliberately not part of this command sequence.

The final verifier closes the exact \(s_1=\ell=0\) one-parameter slice
of that \(Q\)-residual component.  The degree-\(100\) factor of the
\(\mu_6\) norm is removed by the first
\(\det(M_{\mu_6}+zM_{\mu_7})\) coefficient.  On the only remaining
cubic factor, the Kummer modulus specializes to \(u^4\), so the whole
scheme-theoretic fiber is supported at \(u=0\) and is empty after the
chart localization \(u=s_0^{-1}\).  It writes
`artifacts/generated-results/two_pair_sic_bidegree33_Q_cubic_exceptional_factor.json`;
the whole-file SHA-256 is
`69caf2d4b83fc2d70e1ce46b945471b6d0666b0ebbe12de6d20ec660a6e7114a`.
This is an exact slice exclusion, not an exclusion of the full
\(Q\)-residual component.

Three corrected exact closed-fibre calibrations away from that slice use
\(\ell=s_1u-t_1\):

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual.py \
  --through 7 --timeout 20 \
  --specialize s1=5 --specialize ell=7 --specialize u=2 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_fiber_s1_5_ell_7_u_2_exact.json
```

The analogous artifacts at \((s_1,\ell,u)=(7,4,3)\) and
\((11,9,5)\) give the same exact degree-five extension, length-four
fibre, and unit ideal through \(\mu_7\).  They prove that the exceptional
Fitting locus is proper, not that it is empty.

The projective source of that length-four fibre and the smaller
quadratic-remainder elimination are checked with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual.py \
  --through 7 --timeout 60 --projective-probe \
  --specialize s1=5 --specialize ell=7 --specialize u=2 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_projective_probe_s1_5_ell_7_u_2_exact.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_Q_residual_infinity.py \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_Q_residual_infinity_exact.json
```

The second command proves symbolically that the highest fibre-degree
parts of \(\mu_4,\mu_5\) share
\(6s_1us_5-s_6\), that any second infinity point is supported on the
inherited factor \(J_Q=6084\ell^2+4805u^2\), and that the tangent
determinant is a nonzero scalar times \(u^5D\), where \(D\) is the
residual border factor.  The first command confirms full projective
coprimality and infinity length two at the exact closed fibre.  It also
reduces \(\mu_5,\mu_6,\mu_7\) modulo the quadratic \(\mu_4\): their
\((s_6,s_5)\)-degrees become \((1,3),(1,3),(1,4)\), and the four
polynomials generate the unit ideal at that fibre.  These facts justify
the smaller univariate subresultant continuation; they do not yet exclude
its exceptional parameter locus globally.

The closed-point determinant-pencil oracle is reproduced by:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --prime 43 --timeout 60 \
  --specialize s1=5 --specialize ell=7 --specialize u=2 \
  --write-moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_s1_5_ell_7_u_2_mod43.json \
  --singular-output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_s1_5_ell_7_u_2_mod43.singular.log
```

It verifies the pivot open, \(\mu_3\), and the leading border before
constructing the length-twenty quotient and all coefficients of
\(\det(M_{\mu_6}+zM_{\mu_7})\).  This is a modular interpolation oracle,
not a characteristic-zero certificate.

One batched line reconstruction and the two transverse degree scouts are
reproduced by:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --pivot-mode equation --prime 1009 \
  --scan-variable s1 --scan-values 0:500 \
  --specialize ell=7 --specialize u=2 \
  --moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --reconstruct-training-count 400 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_scan_s1_ell7_u2_500_mod1009.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --pivot-mode equation --prime 1009 \
  --scan-variable ell --scan-values 0:1009 \
  --specialize s1=5 --specialize u=2 \
  --moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --reconstruct-training-count 800 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_scan_ell_s1_5_u2_full_mod1009.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --pivot-mode equation --prime 1009 \
  --scan-variable u --scan-values 1:1009 \
  --specialize s1=5 --specialize ell=0 \
  --moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --reconstruct-training-count 800 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_scan_u_s1_5_ell0_full_mod1009.json
```

The first fit validates all 21 rational pencil coefficients on 99
held-out points.  Representative numerator/denominator degrees are
\(209/91\) in \(s_1\), \(420/212\) in \(\ell\), and \(270/100\) in
\(v=u^2\) on the \(\ell=0\) chart.  These are exact finite-field line
certificates and interpolation estimates, not characteristic-zero
multivariate reconstruction.

The corrected sparse unspecialized ratio chart
\(\lambda=(s_1u-t_1)/u,\ v=u^2\) is profiled with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage evaluated --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_evaluated_mod43.json
```

The pivot and moment stages finish quickly and the displayed evaluation
finishes in about thirty seconds.  The next unspecialized basis stage was
stopped after four minutes under a \(3\)-GB cap, so it is not a
reproduction command for a completed artifact.

The bounded pre-pivot quadratic elimination is reproduced by:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage prepivot5 --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_prepivot5_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage prepivot6 --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_prepivot6_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage prepivot7 --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_prepivot7_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage raw_equations --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_raw_equations_mod43.json
```

The first three commands pseudo-divide by the quadratic \(\mu_4\) before
substituting the dense \(s_3\)-pivot.  Their modular linear remainders
have respectively \(25,29,41\) raw monomials.  The last command forms
the three \(s_6\)-elimination equations in about five seconds; their
\((s_5,s_3)\)-degree pairs are \((5,12),(5,10),(6,12)\).
These are exact computations over the displayed finite-field residual
extension, not characteristic-zero component exclusions.  The guarded
`prepivot_cross6` stage confirms that expanding even the smallest
equation after the dense pivot exceeds 180 seconds; the intended next
consumer is therefore a resultant/subresultant implementation that keeps
the pivot linear.

The \(L=1\) trace/norm reconnaissance treats
\((s_3,s_5,t_4)\) as a rank-twelve finite fiber after
\(\mu_3,\mu_4,\mu_5\).  Export the rational-function-field generator
matrices, sample exact multiplication invariants at two good primes, and
replay every sampled joint-rank-drop point against the corrected later
moments with:

```bash
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_boundary_coefficients.py \
  --prime 47 --orders 2,3,4,5,6 --trace-norm --timeout 600
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_boundary_coefficients.py \
  --prime 47 --orders 2,3,4,5,6,7 \
  --trace-norm --trace-samples 1200 --timeout 600
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_boundary_coefficients.py \
  --prime 101 --orders 2,3,4,5,6,7 \
  --trace-norm --trace-samples 250 --timeout 600
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_boundary_trace_candidates.py \
  --prime 47 --timeout 600
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_boundary_trace_slice.py \
  --primes 47,101 --s1 1 --t0 1 --timeout 600
```

Every accepted specialization at both primes has leading ideal
\((s_3^2,s_5^2,s_5t_4^2,t_4^4)\).  The mod-\(47\) scan finds three
reduced length-one common \((\mu_6,\mu_7)\) fibers; corrected \(\mu_8\)
is nonzero at all three.  The commands write their exact finite-field
data under `artifacts/generated-results/`.  This is sampled modular
evidence only, not an exhaustion of the rank-drop locus or a
characteristic-zero nullcone certificate.  The final command is stronger:
on the complete slice \(L=s_1=t_0=1\), it finds length \(1128\) through
\(\mu_7\) modulo both primes and then uses Singular `modStd` with exactness
one to certify over \(\mathbb Q\) that adjoining corrected \(\mu_8\)
gives the unit ideal.  This is an exact slice exclusion, not a global
boundary certificate.  The same command keeps \(s_1\) free and checks the
larger \(L=t_0=1\) hyperslice with `msolve`, and separately keeps
\(t_0\) free on \(L=s_1=1\).  Corrected
\((\mu_3,\ldots,\mu_8)\) is the unit ideal modulo both primes and
directly over \(\mathbb Q\) on both hyperslices.  Each exact rational
`msolve` run uses deterministic sparse linear algebra (`-l 2`) and
outputs the one-element Gröbner basis \([1]\).  Singular's
separate verified modular reconstruction of the first larger hyperslice
hit the recorded \(600\)-second bound, but rational `msolve` supplies
both exact hyperslice certificates.  The corresponding unfixed
seven-variable \(L=1\) exact solve also hit its recorded \(600\)-second
bound, so the result is not a full \(L\)-open certificate.

Modular full-chart reconnaissance is available separately:

```bash
.venv/bin/python scripts/explore_two_pair_sic_bidegree33_full_anchor.py \
  --prime 43 \
  --orders 2,3,4,5,6,7,8,9,10,11,12,14 \
  --timeout 180 --backend msolve --charts s0 \
  --branch s0-A-open-sparse --msolve-linear-algebra 44
```

Recorded Singular and `msolve` runs reach the dense chart ideal quickly
but time out, including after the first pivot and on the common
third-moment pivot boundary. The fully substituted \(A\ne0\) branch
exports as eleven equations in ten variables, but the recorded `msolve`
run terminates inside the solver. Sparse principal-open encodings avoid
that expansion, but the full corrected \(A\)- and \(B\)-open systems still
exceed the recorded bounds. These are computational diagnostics, not
evidence for or against the anchor. Any later modular result remains
experimental until it has an exact characteristic-zero certificate.
The conceptual continuation introduces no new computed certificate:
[`TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md`](extended-geometry/TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md)
now records the bidegree-\((4,4)\) falsification of the formerly all-\(d\)
moment--nullcone conjecture, together with the surviving
Hilbert-series degree-selection layer, the quadratic-anchor target, and
the common-root synchronization induction. It also derives the all-\(d\)
nullcone dimension and the Krull-height lower bound of
\((d+1)^2-3\) global moments, and explains why the balanced higher-pair
analogue is already false for every \(n\geq3\) by the padded three-pair
counterexample.

The same witness has a separate factorial-functional translation and finite
prefix search:

```bash
make verify-factorial-moments
```

It extracts the torus diagonal of every checked power and verifies the
all-order binomial formula behind its zero factorial value.  The diagonal
sequence is not a power sequence: at order two its true translated middle
coefficient is `-8`, while squaring the first diagonal gives `-2`, with
factorial values zero and twelve respectively.  The target also reconstructs
cyclotomic linear forms in exact prime cyclotomic rings and verifies
witness-derived `2r`-term quartics whose first `2r-1` moments vanish for
`r=3,5,7,11`.  The displayed general identities and the homogeneous
linear-class minimality proof are in
[`FACTORIAL_MOMENT_WITNESSES.md`](extended-geometry/FACTORIAL_MOMENT_WITNESSES.md).
These are sharp finite-prefix examples, not counterexamples to the
Factorial Conjecture.

The first exact sparse frontier and fixed binary-form cutoffs are checked by

```bash
make verify-factorial-frontier
```

It exhausts all 3,276 three-monomial supports in two variables through total
degree six and finds no nonzero coefficient point with moments one through
three zero.  It also exhausts all 4,950 pairs of nontrivial monomial orbits
through degree six under the Dvorsky-aligned four-variable involution; odd
moments vanish automatically, but the exact second/fourth-moment gcd has no
nonzero root on any support.  Finally, projective rational Gröbner bases show
that for homogeneous binary forms of degrees one through four, the first
`d+1` moments force the zero form and the cutoff is sharp.  The formulas,
explicit sharp quadratic and cubic witnesses, and strict finite-search scope
are in
[`SPARSE_FACTORIAL_MOMENT_FRONTIER.md`](extended-geometry/SPARSE_FACTORIAL_MOMENT_FRONTIER.md).

The bounded four-variable descendant search can be replayed separately:

```bash
python3 scripts/search_dvorsky_gvc4_bounded.py
```

It exhausts the declared \(40\)-by-\(7{,}448\) normalized lattice slice,
checks pure contractions through order twelve, and screens the full space
of fixed linear multipliers on orders five through twelve.  It finds no
witness in that slice; this is a finite negative search, not a GVC(4) or
SIC(4) theorem.

One natural full-coefficient slice is now closed exactly:

```bash
.venv/bin/python scripts/verify_four_pair_dvorsky_slice_obstruction.py
```

For \(P=(t+a+b+d)(ad+bt)\) and the general ternary quadratic symbol
\(\Lambda=\partial_tR(\partial_a,\partial_b,\partial_d)\), the checker
constructs the first eight pure contractions, proves by exact Singular
reductions that their projective zero set consists of four rank-one square
directions, and verifies a strict weight gap for every direction.  The
written proof in
[`FOUR_PAIR_DVORSKY_SLICE_OBSTRUCTION.md`](extended-geometry/FOUR_PAIR_DVORSKY_SLICE_OBSTRUCTION.md)
upgrades this slice from a bounded lattice experiment to an all-order
no-counterexample theorem for arbitrary fixed multipliers.  It does not
prove GVC(4) or SIC(4).

The all-order nonvanishing proof is written in
[`IMAGE_VANISHING_COUNTEREXAMPLES.md`](extended-geometry/IMAGE_VANISHING_COUNTEREXAMPLES.md);
the generators check the finite artifacts and change-of-variable identities;
the dependency-free audit re-expands the 40-variable witness from scratch.
The dependency chain is `F1 -> LR1 -> IV1`: the foundational collision feeds
the essential cubic quotient, whose named inverse coordinate gives the direct
Image witness, identity slice, and the 40/42-variable Vanishing witnesses.
The parallel `LR1 -> GS1` branch is instead the nonexplicit route to
`not GMC(42)`.  The quantitative rank branch is
`LR1 -> LR2 -> LR3`, with `IV1 -> LR3` supplying the HN consequence
framework for the separate rank-37 realization.  These arrows record logical
or construction dependence.
The collision route retains witness sizes 20/40/42 and rank 37.
Historically, the independent Dvorsky--Long formulas lowered the certified
SIC and unrestricted GVC entries to 5/5.  The current repository witnesses
lower these to SIC pair-dimension 2 and unrestricted GVC dimension 3; binary
GVC makes the latter exact.  The ordinary-Laplacian 40 and homogeneous HN 42
entries are unchanged.  These are witness-ledger values, not literature-wide
minimality claims.
A local proof of the
fixed-dimensional DVEZ/Zhao implication, including Gaussian contraction, the
countable-union step, and formal inversion, completes the nonexplicit route to
`not GMC(42)`; `not GMC(158)` remains the exact conservative Long-route bound.
These high-dimensional bounds are retained as logical-transport regressions,
not active witness searches: Long's five-term three-real example already
settles all dimensions `n>=3`.
It also verifies the uniform weighted-seed Gaussian bridge:
first the standalone Gaussian--Lagrange identity for a nonlinear polynomial
map with nonzero constant terms, then the exact pencil branch, polynomial
determinant correction, and bounded Wick moments for canonical and split
seeds.  It also reverts the mixed-moment generating series to recover a
symbolic quartic and a concrete weighted quintic exactly, verifies the
optimal `N-3` normalized moment coordinates through degree eight and the
variable-scale `N-2` bound, and checks the determinantal reciprocal-series
equations, followed by a separate standard-library reconstruction.  The
all-order completed-ring and residue proof is
[`FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md`](extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md);
the bounded exact script is explicitly a regression rather than a substitute
for that proof.  These checks are part of
`verify-regressions`, not `verify-minimal`.

The first nontrivial exceptional partition complex has a separate exact
moment-coordinate certificate:

```bash
.venv/bin/python scripts/verify_degree_six_gaussian_moment_geometry.py
.venv/bin/python scripts/verify_moment_prony_determinantal_geometry.py
```

The first command derives the irreducible sextic equation of the all-double
component, parametrizes the all-triple curve, and verifies that their four
all-six collision points have scheme-theoretic intersection length two.  It
also transports both degree-six vertical Ritt hypersurfaces into optimal
moment coordinates, verifies the second displayed sextic, and proves that the
`2 o 3` Ritt surface is exactly the all-double exceptional component.  The
second command constructs the equal-multiplicity loci and both degree-eight
Ritt orders from log-Prony and Krylov minors in the optimal moments, compares
their ideals scheme-theoretically, retains the degree-six dual-number
intersection, proves the degree-eight Ritt intersection reduced, and exhibits
the cubic collision thickness in naive mixed-weight Fitting minors.  It then
replaces that marked scheme by the saturated unmarked
Christoffel--Hankel/subresultant ideal, covers the minimal degree-five
`3+2` case and one-node collision strata, and verifies the length-two
degree-eight mixed/all-double
intersection.

The general primitive-merger theorem and its first failure have a separate
exact certificate:

```bash
.venv/bin/python scripts/verify_omitted_intersection_algebra.py
```

It constructs the allocation hypergraph and merger-cycle spaces for the
degree-twelve, degree-eighteen, and first degree-twenty-four faces.  It also
derives the primitive dual-number block, finds the first nonminimal failure
`k[t]/(t^3)` when its root meets a common double atom in degree eight, and
proves that coalescing two pure transfer blocks gives
`k[X,Y]/(X^3,XY,Y^2)`, which has the same length and Hilbert vector as two
dual numbers but a two-dimensional socle.

The underlying Hessian-incidence and Ritt-intersection calculation is replayed by

```bash
.venv/bin/python scripts/verify_hessian_ritt_degree_six.py
.venv/bin/python scripts/verify_degree_six_ritt_atlas.py
```

The second checker refines the Hessian-incidence result on the normalized
seed chart: it computes the `2^3` and `3^2` omitted-value intersections, the
four doubled type-`(6)` collision points, factored affine-sheet boundary cuts,
and clean rational witnesses for all open pieces.

The complete degree-six boundary atlas requires both SymPy and Singular:

```bash
make verify-ritt-boundary
```

It proves that the two Ritt surfaces have respectively two and three exact
affine-boundary curves, supplies a rational Hessian-clean witness on every
curve, and computes the common-curve deletions: one reduced sextic
zero-cluster orbit, two rational plus four conjugate extra-root points, and
four Hessian/type-`(6)` collisions disjoint from the affine boundary.

The first genuine braid of complete decompositions has a separate
scheme-theoretic certificate:

```bash
make verify-ritt-2-complex
```

It builds the Ritt Coxeter 2-complex with commuting-square and braid
relations, verifies the Dickson coefficient map at all six degree-thirty
vertices, and compares the two path ideals around the `S_3` hexagon.  Both
paths have the same smooth `A^2` reduction and normalization.  One path is
reduced; the other has nilpotence index four, with one excess tangent
direction and normalization-defect annihilator `(z^2)`, supported on the
monomial divisor.  It also identifies dual-number and length-five
curvilinear slices of the defect, computes the latter's `K`-adic length
filtration `2,4,5,5`, and verifies that the path tangent dimensions are
unchanged when computed directly in the ambient polynomial and Hessian
coefficient spaces.  The full ideal and doubled-annihilator comparison is
then repeated independently on the opposite `5 o 3 o 2` endpoint chart.
Restoring the omitted linear-coefficient residual leaves every endpoint,
path, and boundary ideal unchanged on both charts, proving exact
scheme-theoretic Hessian transfer for this braid component.
The checker then audits the four remaining vertex charts.  The three
composite-omission sectors `10`, `15`, and `6` have respectively
nilpotence/annihilator data
`(4,z^2)`, `(3,z^2)`, and `(4,z^4)`; the complementary prime-omission path
is reduced in every sector, and opposite endpoint charts agree.
Their annihilator slices have
`(length, embedding dimension, Hilbert vector)` equal to
`(5,1,(1,1,1,1,1))`, `(4,2,(1,2,1))`, and
`(8,2,(1,2,2,2,1))`.  All three have one-dimensional socle; the latter two
are codimension-two Artin Gorenstein complete intersections.  Exact
elimination identifies the three slice algebras as
`Q[u]/(u^5)`, `Q[u,v]/(u^2,v^2)`, and
`Q[u,v]/(u^4,v^2)`.  Their conormal ranks are `1,2,2`, with residue-field
Koszul Tor ranks `(1,1)`, `(1,2,1)`, and `(1,2,1)`.

The coefficient-decorated cellular and Postnikov generalization is replayed
by

```bash
.venv/bin/python scripts/verify_hessian_ritt_cellular_cotangent_prototype.py
.venv/bin/python scripts/verify_degree42_ritt_conormal_transitivity.py
.venv/bin/python scripts/verify_degree42_ritt_postnikov_overlap.py
.venv/bin/python scripts/verify_cellular_postnikov_transitivity.py
.venv/bin/python scripts/verify_hessian_ritt_cotangent_descent.py
.venv/bin/python scripts/verify_ritt_cellular_prototype_completion.py
```

The first command verifies the vertex, move, commuting-cell, braid-cell, and
relative-path totalizations.  The second and third use Singular to prove
degree-forty-two conormal non-splitting, overlap vanishing, and the separation
of non-flat base-change Tor.  The last command is fast exact rational linear
algebra: it validates arbitrary finite equivariant module towers, replays the
degree-thirty one-layer degeneration and the actual degree-forty-two
base-square action matrices, and writes
`artifacts/generated-results/cellular_postnikov_transitivity.json`.
The final command verifies the skeletal descent boundary: the filled braid
is complete in dimension two, while the oriented permutohedron three-cell
kills the topological `H2` line of the four-factor Coxeter two-skeleton
without changing `H0` or `H1`.  It also constructs the normalized
face-poset bars and canonical subdivision maps for the relative half-braid
and filled braid.  Their mapping cones are exactly acyclic before and after
tensoring with coefficient blocks of dimensions `2`, `4`, and `6`; the
same run computes the actual degree-forty-two tangent images of all six
factor charts in Hessian coefficient space.  The vertex ranks are all nine,
the adjacent-move intersection ranks are `(8,5,6,6,5,8)`, and the common
intersection is the Dickson tangent plane plus the Hessian projection of
`(W+1)^36-1`.  Intersecting the four vertex images along each half-braid
then verifies the conormal flag `(5,6,6,7)` for all three opposite-pair
sectors, with composite omissions `6`, `14`, and `21`.  The result is
written to
`artifacts/generated-results/hessian_ritt_cotangent_descent.json`.

The two nonlinear rotated degree-forty-two sectors are expensive,
specialized computations:

```bash
.venv/bin/python scripts/explore_degree42_ritt_rotated_conormal_flags.py --word 237
.venv/bin/python scripts/explore_degree42_ritt_rotated_conormal_flags.py --word 327
.venv/bin/python scripts/verify_degree42_ritt_cut14_postnikov_overlap.py
.venv/bin/python scripts/verify_degree42_ritt_cut21_postnikov_overlap.py
.venv/bin/python scripts/verify_degree42_ritt_cut14_tensor_split_q4.py
.venv/bin/python scripts/verify_degree42_ritt_inverse_limit_sections.py
.venv/bin/python scripts/verify_degree42_ritt_completed_splits.py
```

They reconstruct the cut-`14` and cut-`21` residual ideals, change to seven
normal plus two Dickson-base coordinates, and compute the exact quotient
modulo the fourth maximal-ideal power.  Both have conormal flag
`(5,6,6,7)`, thin path equal to boundary through order four, and common
spectator dimensions `(1,3,6)`.  Their sector dimensions are respectively
`(1,4,9)` and `(1,4,10)`, rather than the existing cut-`6` profile
`(1,5,13)`.  The commands write
`artifacts/generated-results/degree42_ritt_rotated_conormal_jet_237.json`
and
`artifacts/generated-results/degree42_ritt_rotated_conormal_jet_327.json`.
The first two reconstructions also write reusable compressed ordinary-ideal
caches
`artifacts/generated-results/degree42_ritt_rotated_source_ideals_237.json.gz`
and
`artifacts/generated-results/degree42_ritt_rotated_source_ideals_327.json.gz`.
To reconstruct only one cache from scratch, add
`--rebuild-source --build-source-only`; this replaces the selected cache
and does not run the local Singular audit.  The last two commands consume
those caches and perform exact Nakayama and Artin--Rees tests.  They prove
completed thin-path/boundary equality and vanishing of the completed
quadratic overlaps for cuts `14` and `21`, and write the corresponding
`degree42_ritt_cut14_postnikov_overlap.json` and
`degree42_ritt_cut21_postnikov_overlap.json` artifacts.
The two rotated-conormal JSON files remain finite-jet computations; the
separate completed certificates cover both new sectors.
The order-four command consumes the compressed order-four action-matrix cache
and proves that the tensor-presented cut-`14` conormal extension splits
over `B/(tau,zeta)^4`, with dimensions `9 -> 13 -> 4` and cocycle ranks
`32=32`.

The next command consumes the two compressed order-seven action-matrix
caches and constructs orders five and six as quotients of those single
order-seven presentations.  Thus its sections commute with truncation by
construction.  For both cuts `14` and `21`, the dimensions at orders
`5,6,7` are respectively

```text
12 -> 17 -> 5
15 -> 21 -> 6
18 -> 25 -> 7
```

and the cocycle/coboundary ranks are `55=55`, `84=84`, and `119=119`.
The section-difference restriction maps have two-dimensional cokernels, so
these finite splits alone do not imply an inverse-limit split.

The last command consumes the completed two-variable presentation caches
and verifies explicit polynomial sections.  With `u=1+tau`, their generator
images are

```text
cut 14: e4 + (-3*u^2 + 2*zeta)*e6
cut 21: e4 + (-4*u^3 + 8*u*zeta)*e7.
```

Every spectator relation maps to zero in the total presentation, and each
section followed by projection is the identity.  Hence both extensions
split over `Q[[tau,zeta]]`; the completed extension and inverse-limit torsor
obstruction classes vanish.  The earlier extra cut-`21` fourth-jet
quotient-ring dimension is a non-flat base-change/Tor discrepancy, while
the different correction polynomials retain the genuine labelled-sector
asymmetry.  Restriction coherence between different factor charts remains
open.

The last command writes the fully explicit degree-forty-two and
degree-thirty factor/move/labelled-cell diagrams, the totalized complexes,
all certified filtration cohomology rows, and the exact first failed split
reduction.  It verifies uniform cellular `H2=0` for the prototype while
locating the genuine obstruction at degree-forty-two filtration order three:
the sector--spectator extension is non-split and its completed cotangent
connecting morphism is nonzero.  The result is
`artifacts/generated-results/ritt_cellular_prototype_completion.json`.

The direct conductor-first node/cusp ansatz is replayed by

```bash
.venv/bin/python scripts/verify_conductor_first_one_chart_obstruction.py
```

It constructs finite quadratic marked-root algebras over
`Q+t(t-1)Q[t]` and `Q+t^2Q[t]`, verifies discriminant descent through the
node and cusp conductors, and solves the complete reconstruction
polynomiality equations against a prescribed conductor pole.  The systems
are incompatible in every pole order by divisibility; exact Gröbner
regressions cover the first four orders.  The independent unit-group
obstruction shows that the sole-conductor localization is not affine space
after any polynomial stabilization.  The theorem is scoped to the direct
separated one-chart ansatz; multi-chart and ambient-coupled conductors remain
open.  The result is
`artifacts/generated-results/conductor_first_one_chart_obstruction.json`.

The smallest symmetric ambient-coupled escape is replayed by

```bash
.venv/bin/python scripts/verify_conductor_three_boundary_cox_fill.py
```

It replaces the conductor localization by the three-prime fill
`x*y*z=c(t)`, with `c=t(t-1)` for the node and `c=t^2` for the cusp.
Singular verifies the descended hypersurface equations, smoothness of the
nodal normalization, and the exact three-axis singular locus of the
cuspidal normalization.  The remaining exact checks show that the marked
root and reconstruction equations descend polynomially, but the descended
dualizing form pulls back as `Omega_norm/(x*y*z)`.  Affine-space recognition
then fails: the smooth nodal fill has Hodge--Deligne polynomial
`(uv)^3+2(uv)^2-uv`, while the normal cuspidal fill is singular.  This is an
obstruction for the symmetric Cox-product ansatz, not for asymmetric affine
modifications or a distributed target conductor ledger.  The result is
`artifacts/generated-results/conductor_three_boundary_cox_fill.json`.

The general affine boundary-obstruction regressions are:

```bash
make verify-boundary-obstruction-theory
```

The checker compiles one saturated module and one genuine boundary-torsion
module, including regular-element and distinguished-class certificates.  It
then computes the finite normal jets `Q[x]/(x^n)` for `n=1,...,6`: every
transition is surjective, but the least boundary-annihilation exponent grows
as `1,2,3,4,5,6`.  This is the exact control showing why bounded jet
saturation does not imply formal saturation without a uniform exponent.
Independent rational matrices verify the node and cusp conductor pullbacks,
rank-three finite-free tensor descent, strict bounded lifting, and a
non-strict degree-loss example.  The command writes
`artifacts/generated-results/boundary_obstruction_theory.json`.

The conductor-first existence certificate is:

```bash
.venv/bin/python scripts/verify_conductor_first_foundational_cusp_keller.py
```

It begins with `Q[u^2,u^3] subset Q[u]`, derives the translated cusp
`4*S^3+27*V^2=0`, and reconstructs the cubic seed
`H(W)=W^2(1-W)`.  The cubic marked-root discriminant descends through the
conductor and has a triple root at its conductor point.  The checker then
solves the full weighted source/target ledger: all apparent divisions cancel
to `(F1,F2/2,F3/2)`, the reconstruction is `x=-C/E_W`, and the exact
Jacobian is `-1/2`.  Three rational source points have common target
`(-1/4,0,0)`.  This meets the conductor-first existence criterion but
recovers the known foundational weighted mechanism rather than a new stable
class.  The result is
`artifacts/generated-results/conductor_first_foundational_cusp_keller.json`.

The characteristic-labelled Hessian--Ritt linear complex and its
positive-characteristic Frobenius summand are replayed together by

```bash
make verify-ritt-deformation-complex
```

The first checker verifies the characteristic-zero tree differential and
cellular baseline.  The second exhausts normalized one-sided composition
tangents over small fields and verifies the exact Hessian-cutoff trichotomy:
the full \(r-1\)-dimensional kernel for \(H'=0\), the one-dimensional
\(kx\) kernel for \(H'=a\ne0\), and no invisible tangent when
\(\deg H'\ge1\).  The exhaustion is a regression for the written
all-degree degree-additivity proof; it is not itself the proof.

The chart-independent missing-linear-coefficient test is:

```bash
make verify-hessian-synchronization
```

It constructs the canonical lift \(\lambda_{a,b}\) using only Hessian
coefficients.  Exact ambient and canonical-factor-chart ideal membership
proves every multiple intersection through degree `18` is synchronized
scheme-theoretically.  In degree `24`, fourteen pairs reduce directly on
canonical factor charts.  The final outer-cut pair `{2,3}` is certified after
transporting the degree-six Dickson collision through a generic quartic and
changing to `4 normal | 5 base` coordinates; its exact Groebner basis has
size `63`.  Thus every degree-`24` multiple intersection is synchronized, and
each ordinary polynomial intersection is exactly one graph over its Hessian
intersection.  The same target verifies the augmentation-ideal lengths,
point-cotangent homology, and intrinsic Tor ranks of the degree-thirty
transverse sector models.  Finally, five exact degree-`30` pair reductions
with basis sizes `11,6,95,6,11` form the cut spanning tree
`2-6-3-15-5-10`.  Therefore the global all-six degree-thirty intersection is
scheme-theoretically synchronized.  A `4 normal | 7 base` common-refinement
calculation also closes the nested pair `{2,10}` with basis size `4`; five
incomparable two-cut subintersections remain uncertified.

The four larger non-tree pair certificates are intentionally split from the
fast spanning-tree regression:

```bash
make audit-degree30-hessian-synchronization-pairs
```

They certify `{5,6}`, `{6,10}`, `{6,15}`, and `{10,15}` with exact basis
sizes `502,189,12,96`.  Together with the default target, ten of the fifteen
degree-thirty pairs are therefore certified.

The same target runs the rank-two Poisson pre-audit and the independent
completion certificate.  The first verifies that the
single displayed output `R=x(2-3xq)` is exactly the foundational third output
after a polynomial source automorphism, and proves that the naive choices
`S=F_1/2`, `T=F_2` have no polynomial `D`-completion.  The second derives the
pole-cancelling shear `Z -> Z-9Q^2`, constructs exact polynomial `T,D,S`,
checks all six brackets and determinant one, and transports the complete
three-point fiber.  A dependency-free sparse-polynomial implementation then
rebuilds the formulas and separately checks all six brackets, the determinant,
term counts, and collision.  This proves a repository rank-two Poisson
theorem; it does not assert that these are the unavailable manuscript's
formulas.

The all-degree rigidity step behind the transported Hessian cases is
replayed by

```bash
make verify-common-right-factor-synchronization
```

It verifies the triangular top-jet reconstruction for every common-right
degree occurring in degrees `30` and `42`, checks that the two degree
censuses each have exactly three decorated incomparable pairs, and verifies
the characteristic-two dual-number counterexample when the total outer
degree is not invertible.  The theorem itself works over every ring in which
that outer degree is a unit.

The first degree-`42` primary transport certificate is:

```bash
make verify-degree42-hessian-normal-jets
```

On the `{2,7}` pair it constructs the `5 normal | 6 base` common-cubic
power chart and proves over `QQ` that the synchronization defect belongs to
the Hessian residual ideal plus the fifth power of the normal ideal.  Thus
the full six-parameter component synchronizes through normal order four.
The exact basis has size `88`.

The conceptual all-order upgrade on the dense power chart is replayed by

```bash
make verify-degree42-conormal-rees-synchronization
```

The normal Jacobian has maximal-minor ideal `(w0^2)`.  Hence away from
`w0=0` the residual conormal map is onto, complete Nakayama identifies the
completed residual and normal ideals, and the synchronization defect
vanishes at every Rees order.  On `w0=0` the conormal rank is exactly
three, so the existing fourth-order certificate remains the correct global
statement and the all-order primary frontier is confined to that divisor.

The divisor itself has a two-normal-variable Rees reduction:

```bash
make verify-degree42-divisor-rees-reduction
```

Three unit residual pivots eliminate `x3,x4,x5`.  The remaining binary
quadrics have resultant
`(81/256)*w1^4*((t+e1*e2)^2-4*e1^3)`.  Off this resultant their Hilbert
vector is `(1,2,1)`, so the completed normal ideal has cube zero and the
fourth-order defect certificate becomes exact.  The unresolved all-order
locus is reduced to `V(w0,w1)` together with
`V(w0,(t+e1*e2)^2-4*e1^3)`.

Dense opens of both residual branches are closed by

```bash
make verify-degree42-kuranishi-branches
```

On the discriminant branch, normalization exposes one common quadratic
tangent; its cubic obstruction is nonzero on `D(w1*w2*t)`, giving the
initial ideal `(ell^2,ell*s,s^3)`.  On `w1=0`, the two binary cubic
Kuranishi forms have resultant
`-15625/262144*w2^6*A*B`, with `A,B` displayed in the canonical note.
On `D(w2*A*B)` their complete-intersection Hilbert vector is
`(1,2,3,2,1)`.  In both cases the existing membership modulo the fifth
normal power is therefore exact.  The same checker computes the first
subresultant on the exceptional divisors: generically the cubics share
exactly one explicit linear factor on each of `A=0` and `B=0`.  Their
next obstruction is consequently a one-variable quartic restriction,
whereas, on this `w0=w1=0` branch, the further equation `w2=0` is the
sevenfold monomial collision.

The degenerate part of the discriminant branch is closed by

```bash
make verify-degree42-discriminant-quartics
```

When `t=0`, the common-tangent cubic vanishes but the terminal quartic
coefficient is `5*w2/64`.  This remains true at the cusp `e1=t=0`,
where the quadratic ideal becomes a single square.  Consequently the
whole discriminant branch synchronizes on `D(w1*w2)`.  The remaining
support is only `V(w0,w2)` together with `V(w0,w1,A*B)`.

Geometrically, for `W(z)=z^3+w2*z^2+w1*z+w0` and
`U(z)=z*W(z)^2`, the generic `V(w0,w2)` core is the odd polynomial
`z^3*(z^2+w1)^2`, while `V(w0,w1)` has the contact-five core
`z^5*(z+w2)^2`; only their deepest intersection is the sevenfold monomial
`z^7`.  On the `A` and `B` divisors the unique common tangent line is
handled by the common-line residual-intersection theorem: blow up the
normal plane, eliminate the transverse coordinate, and read the first
nonzero coefficient of the resulting one-variable residual series.

The generic `A/B` quartic restrictions are closed by

```bash
make verify-degree42-ab-residual-quartics
```

At the exact characteristic-zero point
`(e1,e2,t,w2)=(1,1,3/5,1)` on `A=0`, the residual scalar is
`-4203/1280`.  At the good-prime point `(1,1,21,1)` on `B=0` modulo
`103`, it is `47`; here `A=1` and the subresultant coefficient
`alphaB=9`.  The latter nonzero reduction excludes an identically zero
characteristic-zero restriction on the irreducible divisor `B`.  Hence
both generic resultant divisors synchronize.  Only the proper quartic-zero
subloci, together with the odd core `V(w0,w2)`, can remain.

The part of the proper quartic-zero analysis supported on the higher-gcd
locus is resolved by

```bash
make verify-degree42-higher-gcd-strata
```

The reduced higher-gcd locus on `D(w2)` is the union of four weighted
curves.  Quartic envelopes close every punctured curve: their Hilbert
vectors are `(1,2,3,2)`, except for one `(1,2,3,3,1)` curve.  The only
point of this locus not closed by the nilpotence cutoff is their common
contact-five vertex `e1=e2=t=0`; there the cubics vanish and the quartics
retain a common cubic factor.  This does not yet exclude additional zero
divisors of the scalar quartic residual away from the higher-gcd locus.

Those scalar residual divisors are factored by

```bash
make verify-degree42-ab-residual-factors
```

On the rational normalization of `A=0`, the residual is
`-75/512*e1^2*w2*(4*e1-e2^2)*P_A/
(e2*(6*e1-e2^2)^3)`, where `P_A` is affine-linear in `w2`.  On `e2=1`,
the normalization of `B=0` over `q^2=-3` is
`e1=(1-r+q*(r-1)*(2*r-1))/2`,
`t=-(1+q)*(r-1)^2*(2*r-1)/2`; there the residual is
`75/1024*w2*(q-1)*(r-1)^2*(2*r-1)*P_B`, with `P_B` also affine-linear
in `w2`.  The `e1=0`, `4*e1=e2^2`, `r=1`, and `r=1/2` factors are
exactly the already-classified higher-gcd branches `P4` and `P3`.
Thus the new degree-one-gcd support consists only of the residual graphs
`P_A=0` and `P_B=0`.  The checker also certifies the solved graph
identities for `w2` and verifies that the apparent coefficient-zero
values `5*e1=e2^2` and `r=3/5` do not add vertical components.

Dense opens of both residual graphs are closed by

```bash
make verify-degree42-ab-residual-quintics
```

The checker follows one terminal equation through its quartic transverse
correction, re-solves the three pivot equations through fourth order, and
computes the invariant fifth residual.  On the `A` graph its exact value is
`-250011279/8192000` at
`(e1,e2,t,w2)=(1,1,3/5,567/100)`.  On the geometric `B` normalization it
is `14 mod 103` at `(q,r)=(10,0)`.  The resulting homogeneous envelope has
Hilbert vector `(1,2,3,2,1)`, so the pre-existing defect membership modulo
the fifth normal power becomes exact on both dense opens.

The complete cutoff chain—from the conormal open through the higher-gcd
quartic strata and residual-graph quintics—can be replayed with

```bash
make verify-degree42-kuranishi-cutoff-chain
```

The two pieces are combined by the single support ideal

```text
k = (w0, w1*w2, A*B*w2).
```

The global non-jet target is `I:k^infinity = I`.  By the shared
[support-saturation principle](verified/SUPPORT_SATURATION_PRINCIPLE.md),
this is equivalent to excluding associated primes of the residual algebra
over `V(k)` and is sufficient after normal completion.  It is strictly
weaker than proving the residual algebra flat over the full Ritt base.

The first exact compression of this target is:

```bash
make verify-degree42-depth-reduction
```

Residuals 5, 11, and 17 are global unit-triangular pivots, not only pivots
on `w0=0`.  They eliminate `x3,x4,x5` exactly and present the same residual
algebra using only `x1,x2` over the six-dimensional base.  The checker also
has exploratory `--method height` and `--method colon` modes for

```text
f = w0 + w1*w2 + A*B*w2.
```

The height mode searches for a codimension-two perfect reduced ideal and a
one-step dimension drop after adjoining `f`; the colon mode tests `I:f=I`
directly.  Neither mode is part of the verified target until its
characteristic-zero computation completes.

The next normal order has an exact good-prime certificate, and one rational
point on the remaining divisor has an exact untruncated characteristic-zero
certificate:

```bash
.venv/bin/python scripts/verify_degree42_transported_27_normal_jets.py \
  --prime 32003 --normal-order 5 --timeout 240
.venv/bin/python scripts/verify_degree42_transported_27_normal_jets.py \
  --base-values 1,2,3,0,5,6 --normal-order 0 --timeout 240
```

They return basis sizes `179` and `8`, respectively.  To attempt the
remaining generic characteristic-zero calculation directly on the
conormal divisor, use

```bash
.venv/bin/python scripts/verify_degree42_transported_27_normal_jets.py \
  --w0-zero --normal-order 5
```

That function-field calculation currently exceeds 300 seconds; the timeout
is a performance boundary, not a failed reduction.

The generated certificates are stored as the conservative
[`79-variable artifact`](artifacts/generated-results/long_bcw_79_counterexample.json)
and optimized
[`33-variable artifact`](artifacts/generated-results/shared_bcw_33_counterexample.json),
together with the
[`24-variable rank-compressed artifact`](artifacts/generated-results/rank_compressed_bcw_24_counterexample.json)
and final
[`22-variable constant-kernel quotient`](artifacts/generated-results/constant_kernel_bcw_22_counterexample.json),
together with the new
[`21-variable essential quotient`](artifacts/generated-results/essential_bcw_21_counterexample.json)
and its
[`20/40-dimensional identity-slice witnesses`](artifacts/generated-results/image_vanishing_counterexamples_20_40.json)
and
[`homogeneous 21/42-dimensional witnesses`](artifacts/generated-results/image_vanishing_counterexamples_21_42.json).
They record the sparse cubic maps, every reduction-step choice, and the three
exact collision points, together with the expanded contraction and quartic
polynomials; regeneration is deterministic.

## Complete active suite

```bash
make verify
```

## Degreewise theorem audit

The former standalone manuscript has been retired. The retained theorem
statement, proof dependencies, and reproduction commands are in
[`DEGREEWISE_MULTIPLICITY_AUDIT.md`](DEGREEWISE_MULTIPLICITY_AUDIT.md).

To compile every standalone paper with the same discovery rule used by CI,
run:

```bash
make verify-papers
```

To retain an environment record and complete log under `artifacts/`, run:

```bash
make verify-logged
```

Generated outputs, bounded scans, and exploratory search programs are not
part of the public proof navigation.  Existing generated artifacts live under
`artifacts/generated-results/`; historical search tools are preserved under
`archive/tooling/`.

## Cyclic and dihedral absolute inverse-Galois audit

The first Programme 3 checker replays the Dickson power-sum recurrence,
primitive derivative, reduced-branch pullback, and odd/even discriminant
formulas through degree twelve.  It also checks the low-degree \(D_3,D_4,D_5\)
cards and the determinant-minus-one derivative-unit suspension:

```bash
.venv/bin/python scripts/verify_cyclic_dihedral_keller_audit.py
```

The all-degree proof and the distinction between geometric and arithmetic
monodromy are in
[`extended-geometry/ABSOLUTE_INVERSE_GALOIS_CYCLIC_DIHEDRAL_AUDIT.md`](extended-geometry/ABSOLUTE_INVERSE_GALOIS_CYCLIC_DIHEDRAL_AUDIT.md).
The bounded replay is a regression certificate, not an exhaustive proof over
all degrees.

## \(\operatorname{PSL}_2(11)\) Keller monodromy action spectrum

The exact checker constructs
\(\operatorname{PSL}_2(\mathbb F_{11})\) as
\(\operatorname{SL}_2(\mathbb F_{11})/\{\pm I\}\), enumerates its natural
degree-twelve action and both exceptional degree-eleven \(A_5\)-coset
actions, and verifies Gassmann equivalence, cross-subdegrees \(5+6\), the two
rigid \((2,3,11)\) Nielsen orbits, and the genus-one/genus-zero
Riemann--Hurwitz split.  It identifies the normalized degree-five/six
components as the degree-\(55\) \(A_4\) quotient of genus one and the
degree-\(66\) \(D_{10}\) quotient of genus two, including both projection
passports.  Its symbolic half checks the corrected Klein Shabat factorization
over \(\mathbb Q(\sqrt{-11})\), the derivative, square discriminant, the
direct correspondence, irreducible boundary factors, the rank-six/rank-two
unit ledger, the exact two-step cubic reduction of the genus-one component,
its conductor-\(121\) model and \(j=-121\), the trace obstruction to an
isogeny with \(X_0(11)\), the canonical adjoint pencil and even
hyperelliptic model of the genus-two component, its two elliptic quotients
and \((2,2)\)-split Jacobian, the explicit \(X_0(11)\) equation and
degree-twelve \(j\)-map, and the exact positive-genus boundary-unit lattices
of ranks \(3,14,17\).  It also constructs both compact boundary pullback
matrices, computes their primitive rank-ten unit images and free cokernels
of ranks four/seven, and tests the minimal two-output derivative ledger.  It
also reduces projection exchange on the residual quotients, constructs
effective bases for all four/seven mask classes, exhausts the degree-five
simple-pole divisors, and certifies the first normal-support degrees two and
three:

```bash
.venv/bin/python scripts/verify_psl2_11_keller_action_spectrum.py
.venv/bin/python scripts/verify_psl2_11_keller_action_spectrum.py \
  --pari-mordell-weil
.venv/bin/python scripts/verify_psl2_11_keller_action_spectrum.py \
  --singular-normalization
.venv/bin/python scripts/verify_psl2_11_normalization_masks.py
Singular -q scripts/psl2_11_c5_boundary_pullbacks.sing
.venv/bin/python scripts/verify_psl2_11_c6_boundary_images.py
```

The second command requires PARI/GP and certifies the rational ranks and
Heegner generators used in the quadratic rank decompositions over
\(\mathbb Q(\sqrt{-11})\).  The third requires Singular \(4.4.1\).  It
exactly normalizes both affine correspondence curves, proves that their
conductor ideals equal their reduced singular ideals, verifies five/eight
ordinary nodes and all fifteen/twenty normalized boundary primes, and
replays the canonical elimination for the genus-two sextic and both
quadratic Cremona transformations from the nodal quintic to the sparse
cubic.  It then runs the exact normalization-module mask replay.  The fourth
is that robust standalone replay: it pins the Singular source hash, rejects
silent procedure failures, verifies the small \(C_5\) formulas, constructs
the \(C_6\) filtered Riemann--Roch spaces of dimensions \(23,29,41\), solves
all seven prime-power evaluation kernels, and imposes the two infinity-five
jets selecting the asymmetric mask.  The fifth labels the \(C_5\) boundary
primes by both projection colors and separates all three normalized node
branches.  The sixth requires Singular and SymPy; it replays the nineteen
direct \(C_6\) residue
traces through both elliptic quotients and recovers the remaining
\(q_2\)-node class from its principal fiber relation.

The action-level classification problem, determinant-one common-target
Keller charts, the natural prime-triangle genus formula, affine-completion
obstructions, and external
Jones--Zvonkin monodromy input are in
[`extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md`](extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md).

<!-- status-consumer: KAS1 f40f8588d37ade00 -->
<!-- status-consumer: KAS2 f56459cc921661ea -->
<!-- status-consumer: KAS3 b95d888270f98c59 -->
<!-- status-consumer: KAS4 45a513f714702919 -->
<!-- status-consumer: KAS5 2baa200b6712564f -->

## Absolute \(D_5\) affine-modification frontier

The degree-five precomputation checks the split derivative and branch
ledgers, the tangency of the two ramification colors, the singular product
and separated Cox fills, and the maximal-minor determinant identity for
affine-linear couplings with two and three new coordinates.  It also checks
the block determinant for arbitrary zero-section thickenings:

```bash
.venv/bin/python scripts/verify_absolute_dihedral_d5_modification_frontier.py
```

The all-degree affine-linear mask theorem, class-group calculations, and
the exact surviving nonlinear search conditions are in
[`extended-geometry/ABSOLUTE_DIHEDRAL_D5_MODIFICATION_FRONTIER.md`](extended-geometry/ABSOLUTE_DIHEDRAL_D5_MODIFICATION_FRONTIER.md).

## Nonlinear \(D_5\) obstruction classification

The nonlinear follow-up checks the complete branch-supported valuation
ledger, its primitive diagonal pole, the normalized target-cusp incidence
and tangent-rank drop, and the translated graph-section obstruction for two
and three auxiliary coordinates:

```bash
.venv/bin/python scripts/verify_d5_nonlinear_obstruction_classification.py
```

The eight-gate classification, proofs, and exact rank-five fibre
certificate requirements are in
[`extended-geometry/D5_NONLINEAR_MODIFICATION_OBSTRUCTION_CLASSIFICATION.md`](extended-geometry/D5_NONLINEAR_MODIFICATION_OBSTRUCTION_CLASSIFICATION.md).

## Canonical \(D_5\) two-mask blowdown

The first construction attempt verifies the determinant-\(\Delta\)
two-mask matrix and its adjugate inverse, proves the constant-linear
coefficient locks, checks the generic genus-two fibre used in the all-degree
automorphic rigidity theorem, and exhausts \(72\) coordinate assignments
from the unchanged and two ramification-incidence source charts.  It also
checks the first nonautomorphic contraction mismatch and solves the minimal
affine-normal tangential class:

```bash
.venv/bin/python scripts/verify_d5_two_mask_blowdown_obstructions.py
```

The chain-rule obstruction, all-degree genus-two proof, and precise
nonautomorphic continuation are in
[`extended-geometry/D5_TWO_MASK_BLOWDOWN_OBSTRUCTIONS.md`](extended-geometry/D5_TWO_MASK_BLOWDOWN_OBSTRUCTIONS.md).

## All-degree dihedral affine-completion obstructions

The uniform replay checks the Dickson branch identity, both components of
the even-degree branch, the determinant-\(\Delta_n\) two-mask blowdown, the
positive-genus automorphic-rigidity gate, the first nonautomorphic cusp
remainder, the affine-normal coefficient locks, and the nonlinear
normal-degree resonance, even factorization, and odd valuation-at-infinity
one-normal no-go gates for \(3\le n\le12\).  It also replays two
resonant-looking false positives and verifies that the bounded
one-normal search compiler has no open route:

```bash
.venv/bin/python scripts/verify_dihedral_all_degree_affine_completion_obstructions.py
```

The odd/even valuation ledgers and uniform proofs are in
[`extended-geometry/DIHEDRAL_ALL_DEGREE_AFFINE_COMPLETION_OBSTRUCTIONS.md`](extended-geometry/DIHEDRAL_ALL_DEGREE_AFFINE_COMPLETION_OBSTRUCTIONS.md).
The bounded replay is a regression certificate, not the proof of the
all-degree statements.

## \(S_4\) collision-frame Keller frontier

The collision-frame checker expands the decomposable absolute map
\(F\circ F\), verifies its determinant and degree-nine tower, and checks the
quartic factorization, six-edge action, discriminant, primitive conductor,
tiny normal form, and one-/two-normal determinant ledgers.  It also verifies
the rational determinant-one cotangent lift, the polynomial relative
logarithmic factorization, the natural resultant and Bezout obstruction
models, and finite-field point-count regressions for the displayed motivic
classes:

```bash
.venv/bin/python scripts/verify_s4_collision_frame_keller_frontier.py
```

The class-group calculations, the all-degree rank argument, and the
affine-space recognition requirements are written proofs in
[`extended-geometry/S4_COLLISION_FRAME_KELLER_FRONTIER.md`](extended-geometry/S4_COLLISION_FRAME_KELLER_FRONTIER.md).
The checker does not construct an ordinary polynomial Keller map for the
six-sheet collision cover.

## \(A_4\) Keller inverse-Galois frontier

The pure-target ledger, two-mask factorization, normalized-boundary
assembly, and root-incidence derivative-split checks are:

```bash
.venv/bin/python scripts/verify_a4_pure_target_ledger.py
.venv/bin/python scripts/verify_a4_two_mask_factorization.py
.venv/bin/python scripts/verify_a4_normalized_boundary_assembly.py
.venv/bin/python scripts/verify_a4_root_incidence_derivative_split.py
.venv/bin/python scripts/verify_a4_chart_unit_rank_four.py
.venv/bin/python scripts/verify_a4_two_mask_local_viability.py
.venv/bin/python scripts/verify_a4_affine_modification_obstruction.py
.venv/bin/python scripts/verify_a4_corrected_boundary_selector.py
.venv/bin/python scripts/verify_a4_boundary_coloring_surgery.py
Singular -q scripts/verify_a4_corrected_boundary_genus.sing
.venv/bin/python scripts/verify_a4_genus_zero_selector_search.py
Singular -q scripts/verify_a4_genus_zero_selector_search.sing
.venv/bin/python scripts/verify_a4_sharp_selector_plane.py
.venv/bin/python scripts/verify_a4_conic_principal_obstruction.py
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --conic-sieve
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --cubic-sieve
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --census-bound 6
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --census-bound 3 --include-q5
Singular -q scripts/verify_a4_degree_twenty_near_selector.sing
```

The normalized-boundary command verifies the determinant-one ambient
completion and the exact obstruction to the resulting automorphic assembly.
It also verifies the first birational nonautomorphic solution of the corrected
log-Jacobian equation, namely homogeneous radial scaling, and the abstract
Jacobian identity used to prove that no choice of two polynomial mask outputs
can repair that radial base after the pure-target lift.  It also verifies
the determinant-boundary matrix factorization organizing the two inverse
masks and the direct-adjugate Jacobian obstruction.  The accompanying note
identifies the actual descended coupling space as the ideal-contraction
quotient \(\iota^{-1}(\iota(B_\pi)\Gamma)/(B_\pi)\).  The checker also rejects
the quadratic polar map at an explicit smooth boundary point and proves that
multiplying the cubic by any one hyperplane produces no linear Jacobian
syzygy and hence no Saito-free quartic.  It now also verifies the
etale-component sieve for coupling: a nonzero contraction class requires
the pulled boundary to be supported on the critical divisor \(WKL\).  For
the normalized nonradial triple the pulled boundary specializes at \(W=0\)
to \(-133z_2^3\), so its contraction module is zero.  For the old boundary,
the same checker proves that \(S^2(Q^2+3QR+9R^2)\) is the unique new
degree-four contraction class.  Elimination and the \(L=0\) dominance
test prove the full formula
\((\mathcal B,S^2(Q^2+3QR+9R^2),S^2P^3)\), which closes every base-fixed
use.  The displayed two-mask pair gives exact polynomial inverse quotients,
but the resulting composite Jacobian is \(W^3K^3L\rho_Vz_1/2\).  The
differential argument on \(\mathcal B=S=0\) then closes every incidence
whose new base boundary remains proportional to \(\mathcal B\), including
nonradial feedback.  The final mod-101 identity sieve is exhaustive through
base degree six and finds only \(\mathcal B\) and \(\mathcal B^2\), so a
genuinely different reduced exceptional boundary has degree at least seven.
For the nonreduced hit, the checker also verifies that \(\mathcal B\) has
ordinary node tangent cone \(3x(x+3y)\) and replays the resolved intersection
ledger whose residual ramification class is \(-E_2\) with intersection
\(-1\) against the pulled-back line.  Together with the cited
totally-invariant-curve theorem, the written argument excludes every
degree-at-most-two triple with
\(\mathcal B(p,q,\rho)=v\mathcal B^2\), including triples using the two mask
variables.  A nonreduced realization must therefore start in degree at
least three and cancel its higher terms.
Finally, the same command verifies the dominant node-chord rechart
\((\mu,\lambda,t)\mapsto
(\mu+\lambda f_1(t),\mu+\lambda f_2(t),\lambda f_3(t))\), its boundary
\(27\mu\lambda^2t^3(t-1)^3\), and its Jacobian
\(3\lambda t^2(t-1)^2\).  The irreducible etale witnesses
\(M,H,N_2,WN_2-1\) force its contraction module to vanish, so a successful
coupling must replace at least one node-chord factor by a selector supported
only on \(W,K,L\) after pullback.  The five coordinate-zero witnesses and
the five irreducible shifted-coordinate witnesses reject all sixty
injective placements of \((\mu,\lambda,t)\), so that selector must be
non-coordinate.
The same command then replaces \(\mu\) by the old boundary and verifies the
first nonzero different-boundary coupling.  For \((\lambda,t)=(P,Q)\), its
boundary is \(27\mathcal B P^2Q^3(Q-1)^3\); the class
\(27P^2Q^3(Q-1)^3S^2C_K\) gives two exact polynomial inverse masks and an
explicit polynomial composite, but its Jacobian remains nonconstant.  The
stronger exterior-form calculation proves that every admissible numerator
for a coordinate pair contains
\(h=\lambda^2t^3(t-1)^3\), whereas the log equation contains only
\(9\lambda t(t-1)\).  All twenty ordered coordinate pairs are therefore
closed for arbitrary polynomial masks, not just for the displayed pair.
The root-incidence command verifies a localized two-coordinate
representation of `1/P'(T)`, generic root-field recovery, and the residual
orientation pole obstructing target-only polynomial pullback.  It also
checks the selected rational root and proves that the ordinary `(U,V)`
pullback retains Jacobian `H^6/(2*Theta*K^6*L^3)`.  Neither command verifies
an ordinary polynomial Keller map.  The rank-four command expands the
correct reciprocal `H^3/(4*K^3*L)` in the quartic root basis, verifies its
common denominator `B^2*rho*sigma`, and checks the resulting localized
two-mask determinant-one suspension together with its three genuine
coefficient-ring pole divisors.  It also resolves their common cluster by
four point blowups, computes the branchwise Newton orders of the full
numerator, and solves the resulting local divisor-allocation intervals.
Finally it verifies the forced selector `T+a^3`, its exceptional Cartier
transforms, and the all-degree obstruction to realizing the exact divisor
with masks in the original polynomial root algebra.  It does not construct
the new affine modification needed to adjoin those exceptional quotients.
The local-viability command verifies that the three components have one
common nontransverse cusp/tangency cluster and that the simple second mask
has local order deficit two.  The affine-modification command follows the
forced quotients through the `E3` and `F` charts, verifies their nonnormal
and singular loci, and checks that the smooth negative-definite full-chain
resolution is nonaffine.  The corrected-selector command factors the first
exceptional quartic, derives the higher-order simple-branch correction,
verifies the exact rational two-mask divisor on every resolved branch, and
checks the unimodular relatively ample coarse deletion
`{F,E3,B,Hhat}`.  It then verifies the simple/triple normalized splitting
above retained `E1,E2`, proves that strict `B,Hhat` give only two
horizontal deleted primes, and certifies that the normalized open has
relative class rank at least two.  Finally it computes the normalized
seven-curve intersection matrix, contracts it to the `A3` chain, verifies
the cyclic order-four discriminant group and an explicit odd curvette, and
proves that the forced strict `B` class has content two in the full dual
lattice.  Thus no boundary basis can contain it; the `B,Hhat` pair also
has determinantal divisor two.  The same command then verifies the
`B`-free split `T*Lchi/Hhat | (T+a^3)/Lchi`, its two irreducible curvette
boundaries, the unique unimodular exceptional completion
`{Fs,R2,Ft,R1,Hhat,Dchi,D14}`, and an effective divisor on that support
whose seven exceptional intersections are all positive.  The final
Singular command computes the geometric genus of the irreducible
degree-16 `Hhat` norm as 13, verifies
`Norm(Lchi)=chi*q14`, and computes genus 20 for the degree-14 component.
Since a smooth completion of `A2` has only rational boundary components,
either positive-genus divisor rules out the surviving `B`-free complement
as an affine plane.  The bounded genus-zero search computes the complete
ordinary-total-degree-three valuation spaces in selected-root degrees one
and two.  It verifies four new exact root-linear selectors and proves that
every exact selector in this ansatz has horizontal norm degree at least
sixteen; root-quadratic terms start at norm degree twenty-four.  The bound
is sharp.  The accompanying Singular command verifies irreducible
degree-sixteen norms of genera twelve and fourteen for two sharp
representatives, genus twelve for one displayed `a^3` perturbation, and
genus ten for the near-selector `(b+6)*T-81*rho`.
Those samples do not classify the remaining degree-sixteen parameter plane,
so no rational selector or affine reconstruction is claimed.  The final
Python command constructs its birational degree-ten strict model, proves
generic absolute irreducibility, and verifies the complete fixed-infinity
tangent hierarchy.  Its two terminal members are absolutely irreducible of
genera ten and nine.  The same command factors the moving critical
determinant as a fixed conic squared times one irreducible degree-23 curve,
analyzes the exceptional conic rank drops, and verifies the first two
rational selectors found on that curve.  The parameter `[77:-16:-8]` is
absolutely irreducible of genus twelve with two additional nodes.  The
parameter `[103:-16:8]` has a rational line component, and its coefficient
norm has a rational conic factor, but the remaining component is absolutely
irreducible of genus ten.  Finally, the command checks all 864 primitive
parameters through height six for absolute irreducibility and computes a
144-member height-three genus census.  Those two censuses are bounded
experiments; the full parameter discriminant is not implicitized, and no
affine-space Keller map is asserted.  The conic-principal command purifies
the rational component to the prime ideal \((q_2,\ell)\).  It verifies
\(\operatorname{Norm}(\ell)=-3q_2R_{10}\), with \(R_{10}\) absolutely
irreducible of genus two, and the birational strict factorization
\(3UK^3G_5/H^3\), with \(G_5\) absolutely irreducible of genus two.  Its
conormal computation finds two quadratic closed points, hence four
geometric points, where the conic prime needs two local generators.  This
proves that no coefficient-one principal divisor isolates the reduced
conic.  It then identifies the cluster curvette with class two in the
local \`Z/4\` group, hence local index two, and resolves the other quadratic
pair as an ordinary two-branch conductor node.  The two exact root-chart
preimages have coefficient Jacobians \`-32/27\` and \`-4\`; the conic is
divisorial on the first branch and only a transverse codimension-two
incidence on the second.  Conductor matching and the principal ideal
theorem therefore exclude a support-only principal divisor at every
positive multiplicity.  The full moving-discriminant image, alternative
rational selectors, and a Keller map remain open.  The degree-twenty line
checker constructs the six strict pullbacks for `q0,...,q5`, proves the
complete fixed `K`, `M`, and chart-`rho` divisibility kernels, and gives an
exact four-minor resultant-gcd certificate for every affine line in the full
six-dimensional space.  Its only rational line is `U=0`; the kernel there
is generated by the known `[103:-16:0:8:0:0]` direction and the nonzero-`q5`
direction `T*(a^2-4*rho)`, which is the old rational conic multiplied by
`T`.  The affine-line problem is closed; nonlinear rational components
remain open.  The conic-sieve command exhausts all 3,875 projective
degree-two forms over \(\mathbb F _5\).  It finds only the reductions of
\(K,M,\rho_V\) and one exceptional point; exhaustive lifting gives 5, 25,
and 0 incidence points modulo 25, 125, and 625.  Thus the exception has no
characteristic-zero lift in the good-reduction chart, while the nonreduced
\(K/M/\rho_V\) neighborhoods and degree-dropping reductions remain open.
The cubic-sieve command factors all 3,906 projective selector members over
\(\mathbb F _5\).  The only irreducible cubic factors are the fixed
\(A_*,L,H\) factors and the known cubic component of
\(T(a^2-4\rho)\).  All 38 non-\(H\) incidence points have no moving-factor
tangent.  The artificial \(H\) plane is the sole nontransverse residue; on
the slice \(\widehat R_3+x\widehat R_4+y\widehat R_5\), its eight bounded
lift counts are 5, 5, 25, 25, 125, 125, 625, and 625.  Those counts are an
experiment, not a characteristic-zero existence or nonexistence theorem.
The next two commands are explicitly bounded factorization
experiments: all 175,680
primitive height-at-most-six parameters on `q5=0`, and all 58,095 primitive
height-at-most-three parameters in the full six-dimensional space, are
reducible exactly on the projectivized `K/M` kernels.  The final Singular
checker proves that the one-jet near-selector has an absolutely irreducible
degree-16 genus-10 norm.  Its explicit exact jet correction has an
absolutely irreducible degree-18 strict curve of genus 31.  These results
narrow the selector and coupling searches; they do not construct a Keller
map.

### Davenport alternating coefficient pencils

The first coefficient class beyond AS14--AS17, and the stronger
all-length-two obstruction, are replayed by:

```bash
.venv/bin/python scripts/verify_davenport_independent_marked_line.py
```

The checker forms the general quadratic--quadratic alternating Jung
coordinate with all translations and lower coefficients present, saturates
the four highest fiber coefficients at the outer quadratic leading
coefficient, and obtains the unit ideal over
`Q[a]/(a^2+a+2)`.  It then proves the unbounded length-two statement from
the three exact Newton-leading coefficients.  These are exact symbolic
obstructions, not bounded searches.  No unit-gate solution survives, so no
candidate is promoted to the nonlinear-`U` translated-incidence equations.

## Plane degree-frontier audit

The universal logarithmic boundary-purity theorem and the two exact
matching-cokernel countermodels are checked by:

```bash
.venv/bin/python scripts/verify_universal_boundary_saturation.py
```

The written proof identifies the logarithmic cotangent cokernel on every
resolved SNC plane Keller compactification as a perfect codimension-one
module.  The checker verifies the finite nodal-tree presentations showing
that the same conclusion does not follow for an arbitrary conductor/gauge
matching matrix.

<!-- status-consumer: UCBS1 824720a8f727bdf8 -->

The log-conductor degree shift, determinant-insufficiency models, boundary
blowup factors, and normalized `(75,125)` terminal bracket are checked by:

```bash
.venv/bin/python scripts/verify_log_conductor_degree_shift.py
```

The written theorem places normalization mismatches in `H_Z^1` and proves
that the complete normalized logarithmic determinant has zero conductor
mismatch.  The checker exhibits two integrable Jacobian matrices with the
same determinant and generic branch data but different nodal `Fitt_1`
profiles.  It is an exact local regression for the new full-matrix/localized-
`c_2` target, not a proof of `JC(2)`.

<!-- status-consumer: LCDS1 5b4d92acd50d6c41 -->

The global logarithmic second-Chern-character identity, boundary-blowup law,
and current F2 budget specialization are checked by:

```bash
.venv/bin/python scripts/verify_logarithmic_ch2_budget.py
```

The command proves the bundle and logarithmic-ramification forms over exact
rationals, verifies that node blowups are log crepant while a smooth boundary
blowup lowers `(K+D)^2` by one, and replays both compiled F2 source and target
intersection calculations.  On the current lower-bound model it obtains
`B_sq(d)=(7*d-8)/2` and `B_dbl(d)=(7*d-13)/2`; after the exact cyclic root
class `27`, their degree-floor virtual residuals are `-10` and `17/2`.
These are diagnostics, not effective lengths or exclusions: an exact
divisorial/finite-length filtration and the target curve, pullback
factorization, and Chern module of the purity-forced affine component are
still required.

<!-- status-consumer: LCHB1 176bf85520516fa6 -->

The exact F2 terminal logarithmic node profiles and the corrected two-blowup
interior attachment count are checked by:

```bash
.venv/bin/python scripts/verify_f2_log_node_profiles.py
.venv/bin/python plane-jc/cas/verify_f2_75_125_global_attachment.py
.venv/bin/python plane-jc/cas/test_f2_75_125_global_attachment.py
```

The first command derives the transverse coefficient orders `0,-2,-2,+1`
directly from the terminal rational functions, verifies the tame exponent
matrices at both nodes of every interior chain, and applies Keller boundary
support to identify the smooth endpoint cokernel as `R/(w^3)`.  The second
and third replay
the corrected 19/31-component principal source lower bounds.  Intentional
regeneration of the pinned JSON uses
`plane-jc/cas/verify_f2_75_125_global_attachment.py --refresh`; its current
whole-file SHA-256 is
`419c970e322b16e1bfb6403dc36b1a38b95eb9a52403def6b7ee067c42fe8ddc`.

<!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->

The exact F2 carrier, principal-arm, and spectator logarithmic profiles are
checked by:

```bash
.venv/bin/python scripts/verify_f2_carrier_log_node_profiles.py
```

The command derives the two carrier target-node normal coordinates, verifies
the order pairs `(-1,1)`, `(-1,3)`, and `(-4,5)`, constructs the exact common
regular fans, and checks that every marked local exponent determinant is
`1`, `3`, or `5`.  It then inserts six missing components on each principal
arm, two squarefree spectator branches, or the five-component double-row
attachment, and rechecks unimodularity, Hodge inertia, and integral
adjunction.  The resulting `27/48`-component graphs are lower bounds; the
upstream carrier-extraction chain, outgoing terminal tail, affine purity row,
and uncompiled global centers are outside this theorem.

The checker source hash changed from
`f353a8a116306b43b5a11de02f59498a243ffe0023098b46429dbbe6a035eb1c` to
`537dac833c311ba8210dd84c90f352fdd3c73c318024f45a6821fa5d118925da`
when its returned graph record gained the exact `log_canonical_square` used by
the downstream Chern-budget replay.  The carrier classifications and printed
`27/48` result are unchanged.  Reproduce under `.python-version` and
`requirements.txt` with the command above.

<!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->

The upstream carrier-extraction logarithmic profile is checked by:

```bash
.venv/bin/python scripts/verify_f2_upstream_carrier_extraction.py
```

The command enumerates the certified degree/edge support, verifies that the
carrier-zero source ladder maps unimodularly to the target rays through
`(5,36)`, and evaluates the exact Keller two-form at the extraction-root
node.  It forces the primitive local model
`(T,z)=(W*U^5,W^3*U^18)` up to units and target shears, hence logarithmic
cokernel `R/(W^3*U^18)`.  It independently counts the branchwise matching
quotient `R/(W^3,U^18)` as length `54`.  This is a degree-one local class,
not an exclusion of `(75,125)`.

<!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->

The cyclic boundary matching/Chern charge and its blowup conservation are
checked by:

```bash
.venv/bin/python scripts/verify_log_cyclic_boundary_blowup_conservation.py
```

The command verifies the symbolic node and smooth-point blowup identities,
replays repeated exact graph blowups, and applies the conserved charge
`D^2/2` to `D_root=3E+18L`.  The raw F2 matching length changes from `54` to
`441` after a root blowup, while the self-intersection correction changes to
`-414`; their sum remains `27`.  This is the untwisted Cartier contribution,
not by itself the actual cokernel contribution.  The following kernel-line
check proves that the latter is `deg(K_root)+27`.

<!-- status-consumer: LCBBC1 b3eb4679f781c55f -->

The cyclic cokernel kernel-line/twist formula is checked by:

```bash
.venv/bin/python scripts/verify_log_cyclic_cokernel_twist.py
```

The command verifies `L=K tensor O_D(D)` from the rank-one restriction of a
rank-two bundle map and the GRR identity
`ch_2(i_*L)=deg_D(K)+D^2/2`.  For the F2 root it reduces the unknown global
contribution to the single formula `deg(K_root)+27`; it does not compute that
kernel-line degree.

<!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->

The contracted-packet kernel/Gauss-degree formula is checked by:

```bash
.venv/bin/python scripts/verify_log_kernel_gauss_degree.py
```

The command checks the basepoint-free pencil model and exact Chern arithmetic.
For a cyclic determinant packet contracted to a target point, the kernel line
is `gamma^*O(-1)` for a morphism to `P^1`; consequently its degree is the
negative of the nonnegative Gauss degree.  The F2 root contribution is reduced
to `27-e_root<=27`.  The checker does not calculate `e_root`.

<!-- status-consumer: LKGD1 8a357250b5005186 -->

Tangential-coordinate trivialization of the cyclic root kernel is checked by:

```bash
.venv/bin/python scripts/verify_log_tangential_kernel_trivialization.py
```

The command keeps a nonconstant unit `b` and verifies directly that
`d(W^3*U^18*b)` is divisible by `W^3*U^18` in the logarithmic basis.  Thus
the fixed target covector `dz` spans the kernel on the full nonreduced root
packet, not only on its reduced components.  It checks `e_root=0` and the
exact cyclic contribution `ch_2=27`.  Noncyclic attachments and the global
Chern cancellation ledger remain outside the claim.

<!-- status-consumer: LTKT1 32ac27318f16c20c -->

The F2 outgoing terminal tail is checked by:

```bash
.venv/bin/python scripts/verify_f2_outgoing_terminal_tail.py
```

The command derives the complete exponent map from the terminal support
halfspaces and verifies that the source rays
`(5,12),(2,5),(1,3),(0,1)` map to the target rays
`(5,2),(2,1),(1,1),(0,1)` with determinant one on every cone.  It also
checks exact middle-node coordinates and the invertible final endpoint after
the harmless target translation `P->P-c0`.  Thus the outgoing tail requires
no further fan refinement and has zero logarithmic cokernel.

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->

The F2 affine-purity frontier is checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_purity_frontier.py
```

The command verifies that none of the certified `27/48` source components
can be the purity-forced affine branch divisor, so the rigorous component
floors are `28/49`.  It checks the degree intervals `6..9375` and
`12..9375`, the parametrization bound `124`, and an explicit coarse purity
ledger at every remaining degree.  The latter is a numerical signature, not
a constructed cover; it proves that generic purity alone does not improve
the degree floors or determine the target curve.

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

The F2 affine target-curve atlas is checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_target_curve_atlas.py
```

The command enumerates the 24 normalization-degree rows `(3k,5k)`,
`1<=k<=24`, verifies the exact curve degrees `5k<=120` and the maximum 194
raw parametrization coefficients, and checks both diagonal cusp and
off-diagonal collision models for the divided-difference ideal.  The
degree-ratio, one-point-at-infinity, and affine-line obstructions used to
prove that every actual candidate lies on this atlas are the external
theorems cited in the canonical note.

<!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->

The first target chart's exact collision and conductor packet is checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_target_k1_collision.py
```

The command derives the collision quartic from the degree-`(3,5)` normal
form, verifies the diagonal criterion and target-value formulas, and checks
the exact witness `(a,b,c,d)=(1,0,0,0)`: quartic discriminant `-400`,
diagonal resultant `25`, tangent resultant `-10000`, and four distinct
collision images.  The genus ledger is `6=2+4`, from the `(2,5)` infinity
cusp and four ordinary affine nodes.

<!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->

The exact `k=1` implicit quintic and conductor-gradient bridge are checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_target_k1_implicit_conductor.py
```

The command derives the twelve-support resultant, verifies its degree-five
top form, and proves
`(F_P(p,q),F_Q(p,q))=C(t)*(q'(t),-p'(t))`, where `C` is the degree-eight
resultant packaging the four collision pairs.

<!-- status-consumer: PF2K1I1 a7582c1e36140840 -->

The fixed-coordinate `k=1` Keller-pullback interface is checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_k1_keller_pullback.py
```

The command undoes every target normalization in one triangular formula,
checks that the cleared expression is a genuine quintic, computes the
carrier residue and eight-jet interface, transports the four node fibers,
and verifies the invertible-Jacobian identity equating the source gradient
ideal with the pulled-back target Jacobian ideal.  Étale reducedness and
normalization base change then identify every affine source singularity as
an ordinary node with normalization defect one over one of those fibers.
The finite-flat fiber
split also checks the per-node bound `N_r<=d-1` and total bound `4(d-1)`.

<!-- status-consumer: PF2K1PB1 6f837229017243c4 -->

The normalized `k=1` carrier-jet factorization and fixed-coordinate terminal
factor audit are checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_k1_carrier_jet_factorization.py
```

The command proves the first four normalized carrier jets recover the four
target parameters, verifies the three forced weighted equations on jets
five through seven, compiles the exact weighted transport from raw fixed
carrier centers through `P0,Q0,Gamma,A,B`, separates the
translation-dependent `C^3` and `s^2`
factors from the invariant terminal leading cubic, and runs a bounded
625-specialization irreducibility regression.  It also factors the raw
seven-jet Jacobian as `3*Res(p',q')`, proving that free target translations
absorb the three normalized residuals on the immersed locus.  On the
`E_6+A_1` subfamily it verifies the closed Lagrange-inversion formula for
all carrier coefficients, the six monomial seven-jet equations, and rank
four after the three target transports; this is a codimension-three fixed
carrier gate.  At its `E_8` endpoint the checker verifies the four
scale-free raw-center equations forming a prime codimension-four complete
intersection.  The general target irreducibility statement uses the coprime
parametrization degrees `(3,5)`, not the bounded scan.

<!-- status-consumer: PF2K1JF1 7bc57f390f0531b5 -->

The primitive coprime carrier saturation and nonimmersion pattern is checked
by:

```bash
.venv/bin/python scripts/verify_coprime_carrier_jet_discriminant_pattern.py
```

The command verifies the universal three-parameter carrier transport and
normal-form count `N_*=m+n-1`, recomputes the exact determinant rows for
`(2,3)`, `(2,5)`, and `(3,4)` (with `(3,5)` supplied by the preceding
theorem), and proves that the generic `(3,5)` nonimmersion packet is one
ordinary cusp plus three nodes.  It factors the conductor as a double cusp
point times a squarefree sextic at an exact witness and checks raw seven-jet
corank one. The all-coprime-bidegree determinant law remains a conjecture.

<!-- status-consumer: PCJDP1 d4c16bb71dfc6b80 -->

The unibranch boundary-attachment `Fitt_1` class is checked by:

```bash
.venv/bin/python scripts/verify_log_unibranch_attachment_fitting.py
```

The command verifies that a minimal transverse SNC attachment over a target
branch of multiplicity `m_C` and local residue index `q_p` has logarithmic
matrix `diag(r,t^(q_p*m_C))`, determinant module
`R/(r*t^(q_p*m_C))`, and actual split cokernel
`R/(r) direct-sum R/(t^(q_p*m_C))`.  Their finite difference has length
`q_p*m_C`.  Over a complete residue-degree-`f` fiber the total is `m_C*f`;
the ordinary-cusp total is `2f`.

<!-- status-consumer: LUAF1 b0279670ffbd3fa5 -->

The complementary smooth-fold versus boundary-node cusp dichotomy is
checked by:

```bash
.venv/bin/python scripts/verify_log_cusp_attachment_dichotomy.py
```

The command verifies the universal monomial fold, including cancellation of
the linear transverse branch term, the factor
`J=const*r*t^((n-m)*q-1)`, and the exact Fitting ideal
`(r,t^(m*q-1+epsilon))`.  For the ordinary cusp it checks
`y^2-x^3=-r^2*(r+3*t^2/4)`, `J=3*r/2`, the reduced log matrix
`[[r,0],[t,r]]`, and the point quotient of length one.  It also replays the
complete-fiber ledger `m_C*f-h+c`; the ordinary-cusp range is `f..2f`, with
`2f` only at the node-saturated endpoint.

<!-- status-consumer: LCAD1 7b9c15d3dfae0337 -->

The `k=1` target-complement monodromy obstruction is checked by:

```bash
sage -python scripts/verify_f2_affine_k1_complement_monodromy.py
```

This optional SageMath/`sirocco` replay computes certified three-strand
braid monodromy on all five immersed collision partitions and on the generic
one-cusp and two-cusp strata.  Every van Kampen presentation reduces to
`Z`, so each stratum requires a second ramified affine component and has
conditional source floors `29/50`.  The replay separately certifies the
first escape: the noncyclic `E_6+A_1` complement admits a transitive
degree-six action whose meridians have cycle type `2+2+1+1`.  Thus the
remaining cusp attack must combine boundary attachment data with the
logarithmic Chern budget.

<!-- status-consumer: PF2K1M1 fafcbb3c2e6ceb2b -->

The concentrated `E_8` endpoint has a dependency-free icosahedral
permutation audit:

```bash
.venv/bin/python scripts/verify_f2_affine_k1_e8_monodromy.py
```

The command enumerates all transitive `S_6` images of the torus-knot group
`<a,b | a^3=b^5>` whose meridian `a^-1*b^2` has cycle type
`2+2+1+1`.  The `720` labeled actions form one conjugacy class, represented
by the exceptional six-point `A_5` action.  Thus E8 complement monodromy
allows, but rigidifies, the two-fixed-sheet packet.  The preferred longitude
has the same permutation as the meridian, so its peripheral orbits force
two distinct `(e,f)=(2,1)` source-boundary rows and exclude one `(2,2)` row.
The same command enumerates all `59` subgroups of `A_5`: the complete
fixed-sheet coset-action degrees in the F2 range are `6,10,15,30`, with
respectively `2,4,6,14` separate `(2,1)` rows.  It also checks the uniform
squarefree/double doubled residual ledgers
`7d-62+4N-4r(b-6)-s_X` and `7d-67+4N-4r(b-6)-s_X`.

<!-- status-consumer: PF2K1E8M1 bbb282c6bcfa62fc -->

The full simple-inertia closure, including nontrivial central normalization,
is replayed by:

```bash
.venv/bin/python scripts/verify_f2_affine_k1_e8_orbifold_atlas.py
```

This exact SymPy computation proves that the universal `M^2=1` cusp
quotient has order `240`, enumerates all `30` subgroup classes, and finds
exactly `13` fixed-sheet F2 actions in degrees
`6,10,12,15,20,24,30,40,60,120`.  The order-four center glues meridian
transpositions into peripheral rows `(2,f)` with `f=1,2,4`; the checker also
verifies the resulting Chern ledgers.  Runtime is about forty seconds.

<!-- status-consumer: PF2K1E8O1 4251750ed4e43c89 -->

The stable complete-chain cancellation is checked independently by:

```bash
.venv/bin/python scripts/verify_f2_affine_k1_complete_chain_budget.py
```

After the full Cartier determinant cycle and conormal kernel are subtracted,
the point budget is exactly `u-1`.  Every simple-inertia E8 row has cusp
lower `2R>u-1`; the checker records the required negative `Fitt_1` deficits
`3,7,5,10,13,17,23,27,33,53,113`.

<!-- status-consumer: PF2K1CB1 5cc386dba344a867 -->

The multi-component affine-ramification identity, complete degree-six
cubic-inertia cusp-group enumeration, terminal `A_6` passport match, and
generic split contracted local packet are checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_global_ramification_budget.py
```

<!-- status-consumer: PF2GRB1 aa3a0efd2e0ff277 -->

The geometric-degree-six terminal cubic classification, the exact
conic/line/constant high-contact factorization, the complete `r=6` cusp
surface, the four-point `r=8` cusp field, and the `r=10` SNC exclusion are
checked symbolically by:

```bash
.venv/bin/python scripts/verify_f2_geometric_degree_six_stein_reduction.py
```

The final `r=6` and `r=8` monodromy exclusions require SageMath with its
optional SIROCCO package.  They certify every continued root path before
enumerating all transposition-valued degree-six representations:

```bash
sage -python scripts/verify_f2_r6_cusp_braid.sage
sage -python scripts/verify_f2_r8_cusp_braid.sage
```

The first Sage command checks the generic `A_2` point and both exceptional
`E_6` scale classes.  The second checks all four complex embeddings of the
exact `r=8` quartic field.  Each row has seven equal/disjoint assignments and
no transitive assignment.  These commands exclude all three normal even
terminal rows; the odd normal and nonnormal conductor regimes remain.

<!-- status-consumer: PF2D6E1 d23d615295a1bf58 -->

The generic contracted-divisor Smith theorem, forced cubic E8 jets,
saturated normal form, length-four cyclic quotient, and global incidence
gate are checked by:

```bash
.venv/bin/python scripts/verify_log_contracted_divisor_smith_classification.py
```

<!-- status-consumer: LCDSC1 07dcd994b4faf092 -->

The remaining isolated-`Fitt_1` sign is closed by:

```bash
.venv/bin/python scripts/verify_log_cokernel_cyclic_submodule_positivity.py
```

The local theorem embeds the cyclic determinant module into every
generically cyclic `2 x 2` cokernel with an effective finite quotient.  Thus
isolated `Fitt_1` corrections cannot repair the E8 deficit, excluding every
one-component simple-inertia E8 completion.

<!-- status-consumer: LCSP1 8658eebeb1d65671 -->

The all-stratum `k=1` conductor conservation is checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_target_k1_conductor_conservation.py
```

The command proves that the collision resultant stays monic of degree eight,
replays the affine delta ledger `6-2=4`, and verifies both the concentrated
`(3,5)` cusp `C=t^8` and the cusp-plus-triple degeneration
`C=t^2*(t^3+1)^2`.  It also checks the étale source bounds `4(d-1)` for the
normalization quotient and `8(d-1)` for the conductor divisor.

<!-- status-consumer: PF2K1CC1 f152c82ef2d54c32 -->

The tame-node packet theorem is checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_k1_tame_node_packet.py
```

The command verifies the split and collided logarithmic matrices, their
Smith type `(1,e)`, and the complete regular fan of the cyclic model
`z^e=x*y`, including all `e-1` exceptional `(-2)` curves and their target
orders `(j,e-j)`.  Kato's chart criterion and log-étaleness of log blowups
then give zero logarithmic cokernel and localized `ch_2` for every fs tame
Kummer toroidal packet.  This is a conditional local classification: it
does not assert that the unresolved F2 source packet is toroidal.  The same
check verifies the general SNC residue-rank gate and the two first-jet
equations necessary for a rank-one packet to have singular determinant
support.

<!-- status-consumer: PF2K1TN1 521fb57f7e6abc1f -->

The affine strict-log-étale resolution theorem is checked by:

```bash
.venv/bin/python scripts/verify_affine_keller_strict_log_etale.py
```

The command verifies the Keller chain rule for node, cusp, tacnode, and
ordinary-triple target curves, checks both point-blowup charts under a
nonlinear determinant-one automorphism, and replays the identity logarithmic
matrix after strict étale base change.  The general theorem uses flat base
change of blowups: every affine conductor packet has zero relative log
cokernel, although its ordinary normalization length remains.

<!-- status-consumer: PAER1 60eb24b2232d159e -->

The affine-purity target puncture and terminal-attachment dichotomy are
checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_purity_puncture_attachment.py
```

The command verifies contact `k` with the target `(5,2)` divisor and the
leading residue formula, replays the terminal third branch fiber, and checks
the direct-SNC exponent matrix.  It then audits the resolved terminal node:
every exceptional valuation above it has positive orders on both target
parameters and therefore maps to one point.  The formal
`lambda=125/729,e=3` comparison is not an available affine-row slot, and no
conditional `29/50` component count follows.

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

The generic `k=1` affine-row logarithmic Chern packet is checked by:

```bash
.venv/bin/python scripts/verify_f2_affine_k1_log_ch2.py
```

The command verifies the cyclic local Smith profile `R/(u^e)`, computes the
kernel from the logarithmic conormal sequence, and tracks the target quintic
through the mandatory `(5,2)` center and up to eight carrier centers.  It
also checks `b=min(ord_u(w|_C),8)` in the carrier-normalized coordinate and
the explicit first-jet equation separating `b=1` from `b>=2`, then combines
the packet with the global budget for every signature
`(e,f,n,b)`.  The conditional degree-floor residuals are
`(e^2*n-2*e*f*(b-7)-20-s_X)/2` and
`(e^2*n-2*e*f*(b-7)+17-s_X)/2`, with exact opposite parity gates.  No
terminal `e=3` attachment is assumed.  On the generic cusp face, exact
minimal packets have incidence-sensitive boundary length `2f-h+c`, ranging
from `f` for smooth unramified folds to `2f` for an all-node fiber.  Booking
it subtracts `2(2f-h+c)` from the two doubled unidentified residual
numerators, without changing parity.

<!-- status-consumer: PF2K1L1 5221f5659fc19729 -->

The fixed-coordinate normalized sparse-support exclusions are replayed by:

```bash
.venv/bin/python plane-jc/cas/verify_sparse_support_exclusions.py
```

For `F=(x+P,y+Q)` with `P,Q` having no terms below degree two, the
arbitrary-degree proof classifies every exact support split `1+q` and
`q+1`, `q<=5`.  All charts are unit ideals except the quadratic, cubic, and
quartic
monomial shear chains; their Gröbner bases force compositions of two shears
with explicit polynomial inverses.  The same command checks explicit
Rabinowitsch unit identities for the balanced `2+2` class in arbitrary
degree: 256 determinant/divergence presence masks reduce to 15 linear
survivors and 20 canonical collision partitions, all exactly inconsistent.
The former `14,653,584`-support census through degree twelve remains as an
independent regression.  For `2+3`, the eleven Keller contributions give
2048 presence masks; 321 are compatible with the derivative/determinant
pattern, 98 survive the linear collision sieve, and the global exact
integer no-singleton formula is unsatisfiable.  Thus every `2+3` and `3+2`
chart has a singleton Rabinowitsch unit identity, and every normalized
support of cardinality at most five is invertible.  For support six, the
unbounded exponent formulas classify `1+5` as the quartic shear chain,
make every `2+4` and `4+2` chart a singleton unit ideal, and leave only the
common `{x^2,x*y,y^2}` exponent support in `3+3`; its saturated coefficient
ideal forces a directional quadratic shear with inverse `id-H`.  A separate
`5,290,000`-support census through degree six finds the same unique `3+3`
collision support.  Thus every normalized support of cardinality at most
six is invertible.  Minimizing support over
all tangent-to-identity affine normalizations shows that every noninvertible
plane Keller map has affine-normalized support complexity at least seven.
This is a support-cardinality theorem in fixed
normalized coordinates, not a new universal degree bound.  Support six is
the stopping point for sequential sparse-layer escalation; the next program
must connect affine support to Newton or boundary geometry.  The proof and
exact claim boundary are in
[`plane-jc/CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md`](plane-jc/CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md).
Intentional artifact regeneration uses `--refresh`.

The exact affine-support/Newton bridge audit is:

```bash
.venv/bin/python plane-jc/cas/verify_affine_support_newton_bridge.py
.venv/bin/python plane-jc/cas/classify_f2_75_125_layers.py
.venv/bin/python plane-jc/cas/reduce_f2_75_125_endpoint_system.py
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py
.venv/bin/python plane-jc/cas/test_f2_75_125_frontend.py
.venv/bin/python plane-jc/cas/generate_f2_modified_system.py --include-equations --output artifacts/generated-results/jc2_f2_modified_laurent_family.json
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py
.venv/bin/python plane-jc/cas/verify_f2_kummer_orbit_transfer.py
.venv/bin/python plane-jc/cas/verify_f2_terminal_residue_cover.py
.venv/bin/python plane-jc/cas/verify_f2_a6_simple_spectator_gluing.py
.venv/bin/python plane-jc/cas/verify_f2_75_125_global_attachment.py
.venv/bin/python plane-jc/cas/test_f2_75_125_global_attachment.py
.venv/bin/python plane-jc/cas/verify_f2_75_125_carrier_wronskian.py
.venv/bin/python plane-jc/cas/test_f2_75_125_carrier_wronskian.py
.venv/bin/python plane-jc/cas/verify_f2_75_125_carrier_specializations.py
.venv/bin/python plane-jc/cas/test_f2_75_125_carrier_specializations.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_nonlinear_forcing.py
.venv/bin/python plane-jc/cas/test_f2_75_125_nonlinear_forcing.py
.venv/bin/python plane-jc/cas/test_sparse_circuit_modp.py
.venv/bin/python plane-jc/cas/probe_f2_75_125_nonlinear_modular.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_tangent_obstruction.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py --regular-gauge --artifact artifacts/generated-results/jc2_f2_75_125_formal_homotopy_regular_gauge.json
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py --prime 61 --rho 19 --y 19 --maximum-order 8 --artifact artifacts/generated-results/jc2_f2_75_125_formal_homotopy_mod61.json
.venv/bin/python plane-jc/cas/verify_common_power_carrier_wronskian.py
.venv/bin/python plane-jc/cas/test_common_power_carrier_wronskian.py
```

The global-attachment command pair compiles the first global `(75,125)`
attachment layer:
the original nonmonomial valuation orders `(-25,5,12)`, the shared six-blowup
carrier and six-blowup principal arms, the complete minimal target boundary
ledger, five marked attachment slots per terminal packet, the three symbolic
degree ledgers, and strict candidate-mode gates for any proposed completion.
Target-valuation uniqueness removes the distinct-target double row.  The
pinned attachment output remains incomplete.  The carrier-Wronskian command
pair then
classify its carrier target row: they extract `(5,36)`, force the unique
squarefree cofactor `R=(v^2-3v+3)/25` or the two double-root parameters
`rho^2-3rho+1=0`, identify the cyclic degree-three and Belyi degree-six
carrier maps, and prove the squarefree simple spectators are unramified.
The lower-band realization, purity row, global degree, and global meridians
remain unknown.
<!-- status-consumer: PF2CW1 a7774b0fa736b64c -->
<!-- status-consumer: PF2GA1 57dea3062b1147fb -->

The carrier-specialization pair performs the exact next handoff.  It routes
the squarefree rational carrier outside the descent-eight double-root
component, specializes every exposed Schur and Hermite matrix at
`rho^2-3*rho+1=0`, and pins the zero-row kernel/cokernel bases through layer
`3`.  The double-carrier and descent-eight defect fields form a quartic
compositum with four embeddings.  The coupled zero rows have a raw successive
forcing-cokernel total of `53`, followed by the residual `7+6` target/Hermite
coordinates.  The nonlinear compiler then substitutes the endpoint circuits,
retains all `294+53` Laurent coordinates, appends the final thirteen
functionals and descent-eight incidence, and routes the squarefree carrier
through its later first-defect ledger.  The modular probe adds the explicit
`a!=0` Rabinowitsch coordinate over `F_31`.  Its localized seed has full
Jacobian rank `214`; the spacing-four staircase has `169` variables,
restricted rank `57`, and a consistent `112`-dimensional affine tangent
fiber.  Its full Jacobian cokernel has dimension `153`.  Of six sampled free
directions, none has a nonzero projected remainder after corrected
descending-column back-substitution.  The next compiler checks all `112`
coordinate lines and eight dense mixed lines through the exact degree-seven
bound; all are zero in the `153` cokernel.  The fixed-Jacobian homotopy then
lifts without obstruction through order `16` over `F_31` and order `8` over
the independent good split reduction `F_61`.  A regular seven-coordinate
gauge removes a repeated pivot-pole mode.  Direct `lambda=1` substitution in
the order-16 regular truncation still leaves `50` nonzero equations, so these
are formal jets rather than modular points or component exclusions.
<!-- status-consumer: PF2CS1 666da98d2d24669e -->
<!-- status-consumer: PF2NF1 cfd1da5136c0b6d0 -->

The last command pair abstracts the carrier calculation from F2.  For every
primitive degree-`k` common-power edge with coprime powers `(m,n)`, it proves
the forced first nonshear descent `k*(m+n-1)+1`, constructs the fixed-carrier
linear reconstruction matrix, and derives the universal three-point Hurwitz
passport.  The pinned partition census through `k=24` is a regression audit;
the theorem itself is symbolic and unbounded.  The resonant `k=2` and
imprimitive-multiplicity algorithms remain separate open continuations.
<!-- status-consumer: PCW1 94b10929118f151d -->

An optional exact characteristic-zero family sample requires Singular and
proves that the `r=5` P-only projected top-gap ideal is also `(1)`:

```bash
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py --extended-r5
```

The contact-only artifact was intentionally refreshed after its strategic
recommendation became historical, without changing the four-stratum census:

```bash
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py --refresh
```

Its current SHA-256 is
`77bffef9fed0ed9749f135a426de945dc27e226d43fce88bf6a45b79bb8a83e5`;
the refreshed checker SHA-256 recorded in `MATH_STATUS.json` is
`952d0955dee25eb96933d36f2510783cfe5610d04ac012bb0836749373ad6684`.
The refresh changes only the verdict context: the old recommended pivot is
now labelled historical and linked to `PF2KO1`, `PF2TR1`, and `PF2GC1`.

The first replay proves that coarse Newton vertices, geometric degree, and
nonproperness data cannot upper-bound affine-normalized support: a
Zariski-open family of triangular automorphisms has fixed coarse geometry
and support lower bound `d-2`.  It then verifies the Kummer chain-rule gate
and the live `(75,125)` F2 terminal characters `P={1,4}`,
`Q={0,1,3}` modulo five, which block constant-Jacobian descent.  The second
replay corrects the Laurent chart to `[t,z]=-z`.  Its upper window contains
all 35 zero layers (`39` through `5`), 665 band pairs, and 978 jet-reduced
linear parameters (973 after normalization).  The corrected post-jet
support-row upper bound is 5,344; exact reconstruction uses the jet equations
and all 165,980 compressed generators.  The same replay now carries the
corner bounds through the complete B0 tail: P bands `-75..15`, Q bands
`-125..25`, 2,418 jet-reduced parameters, 240 zero layers through `-200`,
13,741 band pairs, and 1,327,026 exact compressed generators.  This necessary
over-envelope does not by itself exclude F2.  The replay also proves that the
common-power top root is not an
arbitrary degree-18 polynomial: it has the exact two-parameter form
`H(t)=(1+u+...+u^4)^2*R(u^5)`, `u=1+t`, with `R` quadratic and
`R(1)=1/25`.  It also corrects the former upper-descent claim: the substitution
`p=C0^2*U`, `q=(-9/5)*C0^4*V` was an unproved divisibility restriction.  The
true exact source-band kernels at the first five descents have dimensions
`6,6,7,7,10`, because every P-band direction has the Q follower
`q=-3*C0^2*p`; the extra layer-35 direction is the commuting `C0^4` term.
The actual formal `lambda*C0^(-1)` resonance is at layer 10 and is not an
independent source-band kernel.  The intentionally refreshed layer artifact
also verifies the nonlinear resultant `1701*a^8`, forces root continuation
through descent 7, and leaves the first exact local residual
`27*y^2-9*y+1=0` at descent 8.  Its Q-band-one normalization excludes the
four fixed Kummer double-prime supports.  At this earliest spacing only the
nonzero double-root stratum of `R` remains; an exact normalized `P_3/Q_1` interpolation passes its first
local target jet.  At descent 40 the new 9- and 19-dimensional lower bands
show why `E5=0` is not a valid equation after the target: the correct object
is the layer-zero Fitting row.  The replay now reduces the complete target
operator to twelve local jets and two triangular residues.  It proves that
the local `P_3/Q_1` jet needs an off-grid lowest-`u` correction, derives the
edge equation `A'*D-B*C'=1/5`, and verifies an exact edge witness.  It also
proves that the all-`r` sparse witness completes only as an infinite formal
shear.  It classifies the polynomial order-two repair, proves that it never
terminates quadratically, and gives exact `r=3` unit eliminations for cubic
and quartic termination.  It now reconstructs the exact `v^5` Kummer packet:
the omitted fifth-binomial correction cancels the top conflict, its `8 x 10`
matrix has unit minor `3^5*5^16*e^13` over the rank-two branch, and a second
unit Bezout minor carries the edge recursion through `v^10`.  It also proves
the two all-`r` unit-minor formulas and traces their `18*r-1` pivots to exact
polynomial source combinations.  Their source-degree gaps are `12*r-8` on P
and `24*r-20` on Q, with minimum terminal-edge slack one; at `r=3` there are
`53` lifts with largest degrees `47<75` and `73<125`.  The two exact
confluent-CRT determinants at `0,1,w0` then
quotient the triangular `w=0` conditions and leave a rank-`24` global
Hermite coordinate module over the rank-two candidate algebra.  Its exact
fixed-endpoint cone has 22 target and 25 nonzero-weight layer-zero pairs.
The leading target coordinate is normalized, while the remaining ten
`w=1` rows have affine-linear determinant `75000`; degree-seven followers
divisible by `w^2` preserve the controlled `w=0` jets.  Eliminating those
ten variables leaves thirteen global Hermite coordinates.  It also
integrates layer zero into the length-15 algebra
`B[w]/(w^3*(w-1)^6*(w-w0)^6)`, whose quotient by constants is the exact
rank-14 Fitting residue.  Explicit target minors and a constructive
layer-zero span show why neither row alone contradicts the unrestricted old
B0 bands; their earlier triangular equations must be substituted first.
The new endpoint-reduction replay performs that substitution exactly.  Its
ten pivot circuits contain `1,489` terms and give a degree-eight upper bound
after forming brackets.  It derives `1,172` normalized zero-row coefficient
slots and retains all `1,061` active source coordinates, including the
P-minus-21/Q-minus-11 block which re-enters layer `3`.  Explicit tridiagonal
unit minors then eliminate `134` endpoint-disjoint new-Q coordinates on
layers `39..29`, leaving `219` upper Fitting slots and `927` active source
coordinates.  The surviving task is the coupled Schur/Fitting calculation
from descent `12` (layer `28`) through descent `37`, followed by the thirteen
`w=w0`/residue functionals, not another endpoint rank count.
The endpoint-reduction artifact has SHA-256
`9834ed2ba4e64a2b034a83cd0604140206f1a8192a561509c208d02e0a0ca189`;
its checker SHA-256 is
`762ff5e509abcf7701beea4a836b99853d33777f38cdf088eda356be23695858`.
<!-- status-consumer: PF2ER1 64378dad616fc3f2 -->
The artifact separately enumerates the still
open later first-defect spacings `9..90`.  The artifact
has SHA-256
`96e4fd2ff853fcba9d41a72973ea55a1acb2cd41eefad21e69efd6ae73df8b8b`.
The boundary-handoff replay then gives four exhaustive
contact partitions of 25.  It proves that their multiplicities do not
determine branch scales or finite-normalization rows; even the unsupported
strongest contact-to-row surrogate survives the degree-26 packet budget.
This stops the nonlinear triangular elimination, while the finite
character-resolved B0 system now includes every lower band.  Subsequent Kummer-orbit and
terminal-residue calculations, described below, bypass the contact surrogate
and reopen F2 at the global gluing stage; see
[`plane-jc/AFFINE_SUPPORT_NEWTON_BRIDGE.md`](plane-jc/AFFINE_SUPPORT_NEWTON_BRIDGE.md).
The frontend replay checks the forced chain and terminal normalization
independently.  The modified-system generator reconstructs the published
`r=2` systems and the 14- and 22-coefficient-function `r=3` windows under its
explicitly stated common-power ansatz.  It triangularly eliminates the power
rows, gives the remaining core in `B=A^r` coordinates, and presents the true
residue as an Artinian Fitting ideal: a `16 x 16` determinant for `d=2,h=2`
and the four-generator maximal-minor ideal for `d=3,h=2`.  It also checks the
certified `h=4` alternative, which has zero and two generic residue equations.
Every compact residue has an endpoint-binomial point.  More sharply, the
`d=2` residue is rational on a dense open, while the `d=3` residue has a
smooth fourfold meeting the coefficient torus and an exact cubic-invariant
subfamily.  The uniform congruence-support gate excludes all `d=2,3`
congruence sections under the certified `X^4` weight.  For the surviving
`d=3,h=2` section, the Laurent antiderivative constant is impossible and the
only monomial-bracket survivors lie on an exponent-seven, `lambda=0` ray;
its first member is stored as an exact eight-row formal residual point.
The calculation does not prove the ansatz, `d=2,3`, or the lower Laurent-`y`
support ledger, so it does not exclude the unrestricted family.
See
[`plane-jc/F2_MODIFIED_LAURENT_FAMILY.md`](plane-jc/F2_MODIFIED_LAURENT_FAMILY.md).
The generated JSON has SHA-256
`bca206498c153e41a2f31344015df2ce63890f8b12228b6d0ba2c0970eb87c85`;
its software assumptions are `.python-version`, `requirements.txt`, and exact
characteristic-zero SymPy arithmetic.
The modified-chart checker then derives `gamma=2`, the `d=3` monomial chart,
and every possible nonnegative-`xi` support position from the F2 corner
chain.  It preserves the translated binomial-jet relations instead of
treating the support box as independent: at `r=3` the exact P/Q image ranks
are `74/83` and `196/215`.  A primitive top-band relation excludes the formal
terminal point.  For the full projected top diagonal, the three P gap
equations define a length-27 Artinian algebra and the first Q gap has nonzero
resultant/multiplication determinant in it, so the combined ideal is `(1)`.
Thus every branch of the literal bracket-preserving polynomial projection is
excluded.  This also proves that naive deletion of the negative Laurent tail
cannot be the missing theorem.  The full seed has the explicit source lift
`x*(x*y^5-1)^2*R(x*y^5)`, and the corrected top tangent shows that its
first-five kernel dimensions are `6,6,7,7,10`; `lambda*C0^(-1)` is a
layer-10 formal resonance, not an independent source mode.  Deriving its
nonlinear `F`-tail cancellation remains open.  The
pinned chart artifact has SHA-256
`ac7dbc170cafbcf028079b9ccdb41afd78c333e3850abc0761f64cc056e7d7b8`.
<!-- status-consumer: PF2MCB1 6ff13314e0090f52 -->
The terminal checker first verifies the all-parameter degree-`2r` passport
`(2r-1,1)|(r,r)|(3,1^(2r-3))`, geometric monodromy `A_(2r)`, and discriminant
squareclass `(-1)^(r+1)*(2r-1)`.  It then specializes with the Kummer replay
to the terminal target row with transverse index one, residue degree six,
passport `(5,1)|(3,3)|(3,1,1,1)`, and geometric monodromy `A_6`.  Global
boundary gluing remains open.  The checker also proves that the natural
`A_6` action is four-transitive and primitive (so the residue cover has no
`2`-by-`3` factorization), that its target-fixed deck group is trivial, that
`e=1` gives zero transverse different, and that the residue formula is
parameter-free.  It emits residue-different packet `(4,2,2,2)` and verifies

```text
disc_s(125*s*(s+1)^5-r*(9*s^2+15*s+5)^3)
  = 5^17*r^4*(729*r-125)^2,
```

so the rational model has arithmetic `S_6` over `Q(r)`, geometric `A_6`, and
quadratic constant field `Q(sqrt(5))`.  After rescaling the third branch value
to one it is a Belyi map; its regular geometric `A_6` closure has signature
`(5,3,3)` and genus `25`.  The two target toric nodes have exactly three
preimages in the source-divisor interior, fixing three boundary-attachment
points with different contributions `(4,2,2)`.  The last contribution `2`
is at the source toric endpoint over the smooth third branch value.
The target valuation equality gives
geometric degree at least six, or at least twelve for two distinct packets
over the same target divisor.  Packets over distinct target divisors do not
add.  Since the certified target valuation is centered at infinity, the
affine-companion theorem supplies no `+1`; purity instead requires a separate
affine ramification row.  Global geometric monodromy has `A_6` as a
nonabelian simple composition factor.  Thus the remaining ledger has one
squarefree packet or two identical double-root packets attached to the same
versus distinct target components.  The final command exhausts the stated
two-transposition simple-spectator model: all six degree-seven genus-zero
classes survive and generate `S_7`, while paired-star witnesses show that the
coarse filters allow every larger remaining degree.  Adding the certified
endpoint/interior markings leaves three classes of signature `(5,3,1)` under
the strongest naive requirement that the connector anchor avoid both source
endpoints.  The same replay compares the two order-five structures exactly:
the terminal inertia normalizer is `AGL(1,5)`, its multiplier parity and the
arithmetic `S_6` sign both cut out `Q(sqrt(5))`, but `A_6` is perfect and the
full residue deck group is trivial.  Under the stronger assumption that one
Kummer orbit contributes five disjoint transpositions at one branch value,
Riemann--Hurwitz and a rational connected source boundary force degree eleven
and monodromy `S_11`; matching the five core anchors to the inertia support leaves one unoriented class, or four
Kummer-generator orientations.  Both gluing models are conditional because
sheet specialization is not a toroidal node-gluing theorem and the Kummer
spectators have not been assigned certified branch cycles.  These
calculations do not exclude `(75,125)`.

The expanded terminal checker is pinned in `MATH_STATUS.json` at SHA-256
`baa8fe7abdcf1652bc0a8636437b9505ebc6838bf2eebaf2492c03684fd63cbf`.
The software assumptions remain `.python-version` and `requirements.txt`;
the command above both recomputes the assertions and emits the final pass
marker.
<!-- status-consumer: PF2GC1 6ba3fd9eb6a0bcdf -->

The small deterministic regression of the published candidate tables is:

```bash
python3 plane-jc/cas/frontier_125_150.py
```

The direct unibranch finite-flat attack is replayed by:

```bash
.venv/bin/python plane-jc/cas/verify_unibranch_spectator_models.py
```

It verifies a universal rank-`n+1` family with a clean singular unibranch
boundary collision of length `n` and a separate étale spectator.  The
quartic member realizes both exact `3+1` and `2+2` frontier fibers.  These
are countermodels to a purely local exclusion, not Keller maps: deleting
the principal ramification curve gives `A1 x G_m` with a nonconstant unit,
not the distinguished `A2` open.  See
[`plane-jc/UNIBRANCH_SPECTATOR_COUNTERMODELS.md`](plane-jc/UNIBRANCH_SPECTATOR_COUNTERMODELS.md).

The global quartic Cox-lattice continuation is included in:

```bash
.venv/bin/python plane-jc/cas/test_plane_boundary_exclusion.py
.venv/bin/python plane-jc/cas/test_degree_zero_endpoint_pairing.py
Singular -q plane-jc/cas/quartic_completed_deletion.sing
```

It separates the one-boundary row and the same-target/different-target
versions of the two-boundary row.  Their target-pullback lattices have
respectively index-two, rank-one, and index-two defects.  In every case the
single primitive ramified-boundary character saturates the exponent
lattice; it is also the canonical/different class.  This is a reduction,
not a quartic exclusion.  In the rank-one row, `g=a*s_E^2` defines the
finite target-side normalization input.  Base change `B/A` to that
hypersurface to obtain a rank-four finite-free order and normalize it in
`k(x,y)(s_E)`.  The distinguished source `A2` supplies the complementary
normal Gorenstein hypersurface `g(P,Q)=a*s_E^2`; its coordinate ring
contains the finite normalization and gives the graded Zariski--Main open
immersion restricting over `s_E!=0` to
`A2 x G_m -> X x G_m`.  The completed calculation is now part of the
replay: at the `3+1` cusp and each branch of a `2+2` connector the order is
`a*s^2=r^2*ell`, its normalization adjoins `z=r*ell/s`, and its conductor
and canonical module are `(r,s)`.  The source deletion `D(r)` is locally
compatible in every chart.  Its transitions are compatible too:
`r_i=u_ij*r_j` gives `ell_i=u_ij^-2*ell_j` and
`z_i=u_ij^-1*z_j`.  The revised target is therefore the two-generated
degree-zero global-section algebra and its cusp/connector endpoint
pairing, not nonprincipality alone.  The replay also checks the sharper
module form: the degree `-1` square map into degree `-2` has affine
companion cokernel `k[x,y]/(h)`.  A descended unit or exceptional curve is
now a possible witness, not the statement being assumed.  The second
Python command proves the proposed four-filter semigroup search is not yet
finite: for every `n>=1` the reduced divisor
`(y^2-x^3)*(x*y-y^(n+1)-1)` has cusp pole pair `(2,3)` and connector pole
pairs `(1,0),(n,1)`, while the connector matrix is unimodular and all
degree-zero, rank-one-piece, and odd-square-cokernel conditions survive.
The bounded loop through `n=20` is only a regression for the uniform
algebraic proof.  It also checks the triangular coordinates
`(x,y+x^m)`: the connector pole pairs become `(1,m),(n,m*n)` without
changing the affine plane, graded bridge, or packet data.  Hence no bound
on raw coordinate-generator poles can be intrinsic.  A genuine finite
compiler first needs the marked multivaluation semigroup or an
automorphism-minimal coordinate pair, a degree-four bound on its minimized
pole height, and the conductor equivalence relation pairing connector
endpoints.  On the displayed connector the inverse shear
`u=x-y^n=t^-1, v=y=t` gives the exact minimum height two.
Finally, the two-ended family
`X=t+t^-1, Y_c=t^2+c*t^-2+t` has the same pole row `(1,2)` at both
endpoints for every nonzero `c`, but one quadratic shear cancels both
leading poles only when their residue ratios agree (`c=1`).  The compiler
must therefore retain simultaneous initial-residue data, not only numerical
semigroup values.  The imported
`plane-jc/cas/endpoint_valuation_compiler.py` enumerates every monomial
triangular shear that strictly lowers the two-endpoint pole height, in both
orientations, and terminates by integer descent.  It reduces every displayed
connector to height two and distinguishes the residue-matched and
residue-mismatched rows.  It now also exhausts complete lowering polynomial
shears by recursively retaining forced cancellation terms in strictly
descending degree, even when a leading prefix is height-neutral or
height-increasing.  The witness
`u=t^-1+t^2`,
`v=u^3-t^6+t^5+u^2`
has no lowering monomial shear, but the compiler finds
`P(u)=u^3+u^2`, which lowers total height from eleven to nine.  This closes
the one-polynomial-shear gap.  Proposition 6.6 now closes the reduced
alternating-direction peak gate as well.  At each marked valuation, the pole
change made by the second of two opposite degree-at-least-two factors is at
least the change made by the first.  The factor-height increments along a
reduced Jung word are therefore nondecreasing, so a globally lowering word
starts with a lowering complete factor.  As a bounded regression, the
checker exhausts all 49
ordered nonempty-support pairs on Laurent exponents `{-1,0,1}` and every
alternating two-step monomial shear of degrees one or two with coefficients
`+/-1`.  It finds 16 paths which lower height after a nondecreasing first
step; every initial pair already has a lowering complete polynomial shear,
so no terminal peak counterexample occurs in this grid.
The expanded marked multi-pole experiment is:

```bash
.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 2 --include-linear \
  --extended-seeds --max-seed-terms 2 --scan-all \
  --output artifacts/generated-results/marked_multipole_peak_search.json

.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 3 \
  --extended-seeds --max-seed-terms 2 \
  --output artifacts/generated-results/marked_multipole_peak_reduced_degree27.json
```

The first command tracks exact valuations, initial coefficients, conductor
pairing, and every factor height in two- and three-pole rational charts.  It
finds 166 delayed lowering paths, including three-pole peaks of shape
`4 -> 5 -> 3`; all 166 initial states already admit a complete lowering
shear.  The second command tests reduced alternating degree-two/three words
through length three and polydegree 27 on the extended seed basis.  It checks
93,440 words from complete-shear-terminal states and finds no counterexample.
These exact bounded experiments are regressions for the uniform pole-change
proof, not its logical basis.  The proof, generated-orbit continuation, and
the complementary signed and complete-factor commands are recorded in
[`plane-jc/MARKED_MULTIPOLE_PEAK_EXPERIMENT.md`](plane-jc/MARKED_MULTIPOLE_PEAK_EXPERIMENT.md).
See
[`plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md`](plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md).
The normalization and conductor formulas are written algebraic proofs.  The
Python checker replays their cusp factorization, determinantal identities,
and monomial conductor quotient; the independent Singular command computes
the cusp normalization and conductor and verifies normality of the
determinantal overring.

The conductor-decorated endpoint-semigroup continuation is:

```bash
.venv/bin/python plane-jc/cas/experiment_quartic_endpoint_semigroups.py \
  --max-connectors 4 \
  --max-pole 8 \
  --max-contact 8 \
  --cutoff 12 \
  --output artifacts/generated-results/quartic_endpoint_semigroups.json
```

For each bounded connector row it computes the cusp semigroup
`<2,3>`, the displayed connector polar-bound monoid, the exact signed
two-endpoint valuation semigroup
`{(u,v) in Z^2 : u+v<=0}`, the residue completion, the conductor endpoint
pairing, and the odd-square contact vector.  It also replays the 24 cusp
braid pairs and three connector sheet matchings and checks the rank-one
graded bridge.  The output is a bounded feasibility report, not a quartic
exclusion.  The uniform carrier
`(y^2-x^3)*product_i(x*y-y^(n_i+1)-lambda_i)` proves that the listed inputs,
even with conductor pairing restored, do not bound connector count,
displayed pole parameters, or completed contact.  See
[`plane-jc/QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md`](plane-jc/QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md).

The exact quartic linear-pencil calibration is:

```bash
.venv/bin/python plane-jc/cas/experiment_quartic_keller_pencil.py \
  --output artifacts/generated-results/quartic_keller_pencil_calibration.json
```

It compares the finite-free packet model
`(y,x^4-x^3+x*y)` with its target shear
`(y,x^4-x^3+x*y+y^2)`. The maps have the same `3+1` cusp, the same
`2+2` connector, and the same Jacobian determinant. The script resolves
the symbolic linear pencil on the strata `beta!=0` and
`[alpha:beta]=[1:0]`. For the unsheared map the generic fiber is `C*` and
the zeta function is `1`; after the target shear it has genus one, two
punctures, and zeta `1/(1+t^2)`. Both maps are non-Keller packet
countermodels. The output proves that the packet does not determine the
pencil and that the linear pencil is not invariant under nonlinear target
equivalence; it does not exclude either quartic Keller packet. See
[`plane-jc/KELLER_PENCIL_AT_INFINITY_EXPERIMENT.md`](plane-jc/KELLER_PENCIL_AT_INFINITY_EXPERIMENT.md).

The exact chart-aware boundary localization/Smith-normal-form prefilter is:

```bash
.venv/bin/python plane-jc/cas/boundary_lattice_prefilter.py
python3 scripts/verify_boundary_package_compiler.py
```

The first command now checks both localization orientations: boundary primes
mapping to the Picard lattice of a smooth completion, and units of a
certified affine UFD core mapping to all codimension-one fill primes.  Its
`N=3,5,6,7` balanced wild regressions have Smith torsion `2,4,5,6` and exact
named-class orders.  The second,
dependency-free implementation computes Smith factors from gcds of minors
inside the abstract finite-normalization package compiler and rejects the
stabilized `N=5` block by `Cl=Z/4` and `ord([L1])=4`.  These are exact integer
checks; normality, affine UFD structure and unit basis of the core, and completeness
of the listed boundary are theorem-bearing inputs.  When an extended torus
action is separately certified, the toric specialization is sufficient:
trivial class group and constant units force affine space.
The second command also imports the shared retained-root Euler evaluator.
Only an explicitly certified balanced chart with squarefree nonzero retained
roots and one omitted fierce boundary is eligible; `deg(A)>1` is then
rejected before stage-two reconstruction.

The new pre-coefficient front ends are:

```bash
.venv/bin/python plane-jc/cas/test_intrinsic_a2_boundary.py
.venv/bin/python plane-jc/cas/test_log_boundary_compiler.py
.venv/bin/python plane-jc/cas/test_poisson_square_rigidity.py
```

The first reconstructs the canonical class of a complete `A2` boundary,
checks the adjunction/Noether identities, and audits target pole vectors,
ramification, and intrinsic dicriticals.  It proves that a nonproper Keller
resolution needs canonical free depth at least three.  The second compiles
certified branch scales to toroidal proximity and complete
boundary data.  It extracts the local `(2,1),(3,1),(4,1)` rays from the
published `(72,108)` case tree, distinguishes them from the longer adapted
map-base ideals `(t,x^4),(t,x^6),(t,x^8)`, compiles the isolated source chains
of lengths `4,6,8`, and verifies that the common order-four step collapses all
three cases to the same eight-blowup translation graph.  The `F_4` transition
and affine-plane fill then give a unimodular `10 x 10` source boundary
passing adjunction.  The factor-residue tree is encoded symbolically.  The
unselected order-three factor avoids both the common order-four center and
the filled divisor.  The complete common-graph pole vector has no dicritical
component.  A smooth point of `E3` is the unique one-blowup zero-pole
extension; exact two-step witnesses over `Yinf`, `E4`, `E7`, and `E8` delimit
that numerical minimality claim.  The first weighted Wronskian instead
forces the actual `E3 intersect E4` cluster with ten simple children.  At the
plane-return corner, the Poisson-square edge produces a quartic common
factor; all five root-partition fans compile with complete matrices,
differents, conductors, and ramification vectors.  The primary split-factor
formula and alternate-factor chart then select the quadruple-root package
and control its transverse terms.  Both terminal cases have the same
23-component boundary with one degree-twelve dicritical, so the
chain-to-boundary gap is closed.
Its JSON report now includes `retained_root_euler_gate` and
`passes_search_gates`.  The regression checks nonlinear rejection, linear
passage, uncertified fallback, and absent/not-applicable input before
downstream boundary searches are launched.
The third classifies the entire
geometric reduced three-layer Poisson-square locus into the tangent closure
and the `C=0`, `A=0` components, with generic multiplicities `2,3,1`.  Its
exact principal-chart audit proves `I:d0^infinity=I`, excluding associated
primes on `d0=0`.  The normalized `d3,d2` colon filtration classifies the
complete associated-prime set: three minimal primes, three embedded
intersection surfaces, and two embedded core/intersection curves.
The four fast plane checks above run under:

```bash
make verify-plane-jc
```

The independent Singular scheme checks are:

```bash
make verify-plane-poisson-radical
make verify-plane-poisson-primary-charts
make verify-plane-poisson-separators
make verify-plane-poisson-primary-filtration
make verify-plane-poisson-filtered-modules
```

The exact 90 MB certificate archive and extracted replay source are pinned
under `plane-jc/external/zenodo-21479814/`.  Attachment hashes, environment
versions, the full replay command, and the independent hard-certificate
command are in
[`plane-jc/cas/README.md`](plane-jc/cas/README.md).  The independent checker
does not import the primary CAS or generation modules.

For the same archived exact replay with portable process-level CPU
parallelism and without changing the pinned external snapshot, run:

```bash
.venv/bin/python plane-jc/cas/verify_72_108_exact_fast.py --jobs 4
make verify-plane-72-108-exact-fast PYTHON=/absolute/path/to/venv/bin/python
```

GPU backends are intentionally deferred.  They may later be used for modular
or bounded discovery workloads only when their output has a separate portable
CPU verifier.

### Dessin-first no-vertical `(72,108)` closure

The no-vertical Laurent polygon has a separate exact reconstruction that does
not form a Gröbner basis of the original coefficient system.  It enumerates
the five dessins of passport
`(2^10,1)|(3^7)|(17,1,1,1,1)`, reconstructs their single `S_5`-transitive
quintic coefficient model, and verifies the normalized Belyi factorization.
It then compiles the successive linear maps for `(B,E)`, `(C,F)`, and `G'`.
Their shape/rank pairs are `(20x19,17)`, `(20x21,18)`, and `(20x12,12)`,
with complete kernel dimensions `2`, `3`, and `0`; the constant of `G` is a
one-dimensional target-translation kernel.  The terminal calculation uses
Singular only for 25 sparse degree-three/four equations in five deformation
parameters.  Saturation by the required `B_8!=0` open gives the unit ideal.

Generate the pinned artifact with:

```bash
.venv/bin/python plane-jc/cas/jc2_degree108_belyi_deformations.py
make refresh-plane-72-108-belyi-deformations
```

Replay the exact coefficient reconstruction, Galois audit, linear maps, and
terminal Singular unit ideal with:

```bash
.venv/bin/python scripts/verify_jc2_degree108_belyi_deformations.py
make verify-plane-72-108-belyi-deformations
```

Add `--quick` to the direct verifier command to skip Singular while replaying
all preceding exact stages against the pinned terminal result.  The canonical
claim boundary and formulas are in
[`plane-jc/JC2_72_108_BELYI_DEFORMATION_CLOSURE.md`](plane-jc/JC2_72_108_BELYI_DEFORMATION_CLOSURE.md).

### Case-1 full lower-band continuation

The archived Case-1 replay stops after bracket layer `-3`, at
`P:z^-5,Q:z^-4`, because its first thirteen compatibility equations already
give the contradiction.  The following quick audit pins the resulting
full-band ledger and exactly replays the first formerly omitted layer:

```bash
.venv/bin/python scripts/verify_case1_full_band_continuation.py
```

The specialized complete exact replay through `P:z^-8,Q:z^-12` is:

```bash
.venv/bin/python scripts/continue_case1_full_bands.py \
  --stop-layer -11 \
  --ledger-output \
  artifacts/generated-results/case1_full_band_continuation.json
```

It reads but never mutates the pinned external checkpoint.  Exact
`python-flint` arithmetic gives full column ranks
`11,9,7,5,4,3,2,1`, zero nullity in every layer, the same six parameters,
and 66 additional compatibility equations.  The lowest convolutions make
the full replay a long specialized calculation; it is intentionally outside
the routine check suite.  See
[`plane-jc/CASE1_FULL_BAND_CONTINUATION.md`](plane-jc/CASE1_FULL_BAND_CONTINUATION.md).

## Shared `JC_2`--`HC_4` isotropic boundary bridge

The combined programme is
[`JC2_HC4_SHARED_BOUNDARY_PROGRAM.md`](JC2_HC4_SHARED_BOUNDARY_PROGRAM.md).
Its exact first calculation identifies the cotangent determinant
\(\det\operatorname{Hess}(tP+mQ+H)=J(P,Q)^2\), the first isotropic Schur
remainder \(-\Phi_{mm}R(P)\), and the quartic packet's reduced conormal
residue \(2\ell\). It also verifies that every source-only Hessian direction
preserves the cotangent determinant and gives a square-zero relative
endomorphism. Together with `HC4MR1`, this is the exact restricted
equivalence `HC4MR2`: `JC2` is equivalent to pencil-admissible `HC4`.
The completed continuation computes the \(3+1\) cusp
and both \(2+2\) connector initials, proves that the relevant positive
associated-graded conductor maps are isomorphisms, and finds
\(\operatorname{Obs}_{\rm pair}=0\) for all 72 monodromy-compatible
labellings.  Run:

```bash
.venv/bin/python scripts/verify_jc2_hc4_isotropic_boundary_bridge.py
.venv/bin/python scripts/verify_jc2_hc4_global_jet_transport.py
.venv/bin/python scripts/verify_hc4_rank_one_pencil_recognition.py
.venv/bin/python scripts/verify_hc4_diagonal_rank_one_pencil_obstruction.py
```

The third command verifies `HC4MR3`: the coefficient quadrics of
`ell^T*adj(Hess(psi))*ell` define the constant-null-covector scheme; any
projective point supplies a square-zero rank-one pencil. It checks the
universal determinant identity, a nonlinear cotangent control, an oblique
linear rechart, and the four exact projective Groebner charts. This is a
sufficient pencil-recognition gate, not a proof that every HC4 potential
passes it.

The fourth command verifies `HC4MR4`, the first application to the direct
degree-five open problem.  It retains arbitrary lower quartic, cubic, and
quadratic coefficients on the diagonal rank-three Schur packet and extracts
six lower-independent metric coefficients.  Their three channel squares and
three immutable `t^2` squares prove that every nonaligned prolongation has an
empty constant-null-covector scheme.  The calculation rules out this
rank-one recognition method on the packet; it does not rule out the packet
or higher-rank/nonlinear pencil directions.

This does not verify isotropic-flag recognition for an arbitrary
four-variable potential.  It proves instead that the proposed local
paired initial-conormal cokernel cannot be nonzero without an additional
global transport between the two connector jet lines.  The second checker
shows that the marked affine-line normalization makes their projective jet
ratio intrinsic, with exact quartic value \([-1:1]\), but that the
dualizing-residue quotient is annihilated by the conductor.  The abstract
cusp-node family \(R_\lambda\subset k[T]\) realizes the varying ratio
\([-1:\lambda^2]\).  Its symmetric member is the exact two-generated plane
carrier \(k[T^2,T^3(T^2-1)]\), with equation
\(y^2=x^3(x-1)^2\), the anti-diagonal ratio \([-1:1]\), and the required
cusp/node conductor.  Its actual vertical braids are \(\sigma^3\) at the
cusp and \(\sigma^2\) at the node, so van Kampen imposes both the braid and
commutation relations.  These force equal meridian transpositions and
exclude a connected degree-four cover over this structured carrier.  Other
carriers still require their own global braid factorization.

## Direct Schur-descent audit for `HC_4`

The reusable six-to-five-variable construction is imported in
[`MENG_YANG_SCHUR_DESCENT_BRIDGE.md`](MENG_YANG_SCHUR_DESCENT_BRIDGE.md).
For \(\Phi=tA+B\), it separates the two exact hypotheses--constant bordered
Hessian determinant and identically singular reduced pencil
\(\operatorname{Hess}_w(B+sA)\)--from the automatic conclusion
\(\det\operatorname{Hess}(B+\lambda A^2/2+\mu A)=-\lambda c\).
The doubled-Keller row count proves the singular-pencil hypothesis directly,
so no homogeneous-face calculation is needed to reconstruct that descent.
The commands below concern the genuinely additional attempt to descend from
five variables to four.

Two algebraic strengthenings and the resulting research gates are in
[`SCHUR_DESCENT_CONTINUATIONS.md`](SCHUR_DESCENT_CONTINUATIONS.md).  The
exact scalar formula is
\[
\det\operatorname{Hess}\psi_{\lambda,\mu}
=\det M(\mu+\lambda A,w)-\lambda c,
\]
so a fixed descent needs only the specialized reduced determinant to be
constant.  The simultaneous \(r\)-pivot theorem uses the corank condition
\(\operatorname{rank}M\le m-r\); at \(n=3,r=2\) it gives a direct
six-to-four-variable bridge from any Keller collision admitting a jointly
affine two-dimensional source block.  The existing linear-direction audit
rules out such a block for the foundational gauge family, and the general
codimension-one block-affine theorem `SDX2` proves that every three-variable
Keller map with such a block is a polynomial automorphism, even after a
nonlinear source rechart.  Only mixed source--dual or coisotropic pivots
survive this route.

The parameter-uniform affine and low-degree graph audit of the v2
Meng--Yang family is:

```bash
.venv/bin/python scripts/verify_hc4_meng_yang_graph_obstructions.py
```

It verifies the three affine-normal chart coefficients
\(793152L^4,2160L^4,2160L^4\), the zero-dual-normal rank obstruction, and
the degree-at-most-three graph chain

\[
 \rho=-89/16,\qquad \sigma=-\tau,\qquad
 [t^5]D_R(0,t,t,0)=197LN^3/4.
\]

For degree four it retains every graph coefficient and computes the full
two-slope pencil \((0,t,ct,dt)\) only through degree eight.  The leading
square kills all ten \(x_1\)-free quartic jets containing \(y_2\), including
\([x_2^3y_2]R\).  The next coefficients reduce the remaining slice to

\[
 160\rho^2+1968\rho+6021=0,
 \qquad \operatorname{disc}=576\cdot34.
\]

At degree four the same truncated determinant then forces
\([y_1^3]R_3=[y_2^2]R_2=0\) and
\(8\rho^2+99\rho+279=0\).  The two quadratics have resultant
\(16959456\), so no branch survives over any characteristic-zero field.
These results are `HC4MYA1`, `HC4MYG3`, and `HC4MYG4`; degree five is the
first not fully classified single-graph degree.  See
[`HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md`](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md).

Verify the four exact all-degree frontend statements with

```bash
.venv/bin/python scripts/verify_hc4_all_degree_frontends.py
```

This checks the degree-free support and parity identities behind `HC4FSD1`,
the top-face and complementary-minor proof of `HC4FSD2`, the exact
all-lower-layer coefficient identities and unrestricted binary regressions
behind `HC4FSD3`, and the exact Taylor/cofactor proof of `HC4MYGJ2`.  The
original discovery regressions are
replayed with

```bash
.venv/bin/python scripts/research_hc4_all_degree_frontends.py \
  --minimum-degree 4 \
  --maximum-degree 8 \
  --maximum-normal-order 12 \
  --output artifacts/generated-results/hc4_all_degree_frontend_experiments.json
```

The finite Singular tables cover degrees four through eight and normal orders
through twelve.  They are regressions for the all-degree proofs, not their
logical basis.  `HC4FSD2` remains scoped to the minimal tower; arbitrary
lower layers can change its later faces, while `HC4FSD3` proves that they
cannot remove the exact binary obstruction.  See
[`HC4_ALL_DEGREE_FRONTEND_EXPERIMENTS.md`](HC4_ALL_DEGREE_FRONTEND_EXPERIMENTS.md).

<!-- status-consumer: HC4FSD3 1107bc6ff58456f5 -->

The first polynomial-termination and \(x=\infty\) frontends are checked by:

```bash
.venv/bin/python scripts/verify_hc4_meng_yang_polynomial_termination.py
```

This proves `HC4MYPT1`, the degree-free backward terminal equation

\[
 (d-1)U\det H_U-(d+1)\nabla U^{\mathsf T}
 \operatorname{adj}(H_U)\nabla U=0,
 \qquad d=2n+6,
\]

for a last normal coefficient \(R_n=U\ne0\).  Its highest tangential form
has zero ternary Hessian, so it is a cone after scalar extension.  It also
proves `HC4MYPT2`: the pure weight-five ansatz
\(R=(N/L)y^5f(xy)\) has no polynomial constant-determinant member, because
the dominant quartic coefficient is
\(-144a^4(m^2-m-21)\) for \(f\sim az^m\), and the resonance discriminant is
85.  The same checker proves `HC4MYPT3` for the coupled packet

\[
R=y^5f(xy)+y^3p\,g(xy)+y^3h(xy)+y^2q\,j(xy)
  +yp^2\ell(xy).
\]

If \(A,B,C,D,E\) are the five leading coefficients at a common maximal
degree, its terminal equation holds exactly when \(D=0\) or \(B=E=0\).
For `HC4MYPT4`, the checker extracts the complete 2,348-term weight-six
equation, verifies the cancellation of \(M,P,Q,h\), and compares all 396
Newton sectors.  Its upper envelope is

\[
\max\{5a+e+2d+22,\ 4a+2b+2d+22\}.
\]

The strict \(g^2\)-dominant chamber is empty by a negative-discriminant
calculation.  The other strict chamber first reduces to the explicit
resonance recorded in the note.  When \(j\) is constant, both strict
chambers first force \(j=3\); the checker substitutes this value, reduces
the remaining 512 terms to 121 Newton sectors, and excludes both chambers
at the next face.  For `HC4MYPT5`, reduction modulo \(a\) writes
\(30(d-2)=ka\); positivity leaves \(0\le k\le12\), and the thirteen exact
quadratic discriminants are nonsquares.  Hence there is no positive-\(a\)
integral resonance.  The endpoint \(a=0,d=2\) survives one further face
only on \(e=2b+2\), \(b\ge1\), with
\(360AE+B^2D(b+5)^2=0\).

For `HC4MYPT6`, the checker verifies that the full coupled equation is
affine in \(\ell\), contains no \(\ell'\) or \(\ell''\), and that its
coefficient factors as
\(2z(2Lz^5f-3Nz+2N)\mathscr Q\), where the nonzero 311-term factor
\(\mathscr Q\) is independent of \(g,h,\ell\).  It also checks the exact
axis remainder and the nonzero wall multiplier
\(-32L^4A^5D^2P_F(a,d)\).  Thus the balanced wall is unit-triangular in the
lower \(\ell\)-coefficients, while the resulting four-function divisibility
remainder and the exceptional ridge remain open.  Finally, exact truncated
recursion keeps the formal branch above the
collision-containing plane-flat near miss nonzero through normal order five,
and the complete terminal bracket is nonzero at every order one through
five.  The higher corrections are not asserted to retain the marks.  This
is a bounded computation, not an all-order nontermination theorem.  See
[`HC4_MENG_YANG_POLYNOMIAL_TERMINATION.md`](HC4_MENG_YANG_POLYNOMIAL_TERMINATION.md).

The reverse scalar-pivot classifier and simultaneous matrix-pivot equation
builder are:

```bash
.venv/bin/python scripts/verify_hc4_reverse_schur_descent.py
```

The command generates
`artifacts/generated-results/hc4_reverse_schur_descent.json`. It verifies
the exact identities supporting `HC4RSD1`. An identically singular scalar
pencil whose generic kernel line is constant in the four reduced variables
has generic corank one, and the
bordered-unit equation forces a common constant kernel direction of `A`
and `B`.  Every reduced collision fiber is then a three-variable
constant-Hessian gradient fiber, hence a singleton by `HC3`.  This closes
the homogeneous scalar cone-pencil stratum and intersects it trivially with
all 318/306 live affine-degree-two/three projective-polar rows.  The same
checker constructs the exact corank-minor, integrability, collision, and
Schur equations for matrix pivots.  Nonhomogeneous scalar pencils with an
`x`-moving kernel line, nonsingular scalar pencils with determinant-term
cancellation, and moving matrix-pivot kernel planes remain open.  See
[`HC4_REVERSE_SCHUR_DESCENT.md`](HC4_REVERSE_SCHUR_DESCENT.md).

Continue from constant kernel directions to the first affine moving-kernel
stratum with:

```bash
.venv/bin/python scripts/verify_hc4_affine_moving_kernel_pencils.py
```

Together with the classification argument in the accompanying note, these
exact calculations establish `HC4RSD2`.  The adjugate Piola identity
classifies a unimodular affine kernel vector as either constant or, up to
affine coordinates,
`(z,1,0,0)`.  The checker integrates the latter kernel without a degree
bound, obtains
`f=y*C(z)+(x-y*z)*C'(z)+G(z,w)`, and imposes the complete pencil and
bordered-unit equations.  They force the normal form recorded in
[`HC4_AFFINE_MOVING_KERNEL_PENCILS.md`](HC4_AFFINE_MOVING_KERNEL_PENCILS.md),
whose every scalar Schur descendant has an explicit triangular polynomial
inverse.  The generated artifact is
`artifacts/generated-results/hc4_affine_moving_kernel_pencils.json`.
The affine branch left by `HC4RSD2` is parameter-moving, with primitive
kernel form `v0(x)+s*v1(x)`.

Close that parameter-moving affine branch with:

```bash
.venv/bin/python scripts/verify_hc4_parameter_moving_affine_kernel_pencils.py
```

Together with the proof in
[`HC4_PARAMETER_MOVING_AFFINE_KERNEL_PENCILS.md`](HC4_PARAMETER_MOVING_AFFINE_KERNEL_PENCILS.md),
this establishes `HC4RSD3`. The bordered adjugate bounds a primitive
kernel's parameter degree by one. Piola reduces its affine linear part to
rank-one compression pencils; proportional and common-image pencils collapse
to a fixed line, the constant-at-infinity corner is inconsistent, and the
sole moving common-covector integral has Hessian rank at most two. The
generated ledger is
`artifacts/generated-results/hc4_parameter_moving_affine_kernel_pencils.json`.

The first degree-unbounded nonlinear continuation is:

```bash
.venv/bin/python scripts/verify_hc4_univariate_shear_kernel_pencils.py
```

This supports `HC4RSD4`. For an arbitrary nonconstant `P(z,w)`, the fixed
kernel `v=(P(z,w),1,0,0)` integrates completely to
`f=x*a(z,w)+y*b(z,w)+G(z,w)` with `db=-P*da`. The bordered unit kills the
curvature of the common transverse composite, forcing it onto one linear
form. The resulting univariate normal form has an explicit triangular
polynomial inverse for every Schur descendant. See
[`HC4_UNIVARIATE_SHEAR_KERNEL_PENCILS.md`](HC4_UNIVARIATE_SHEAR_KERNEL_PENCILS.md).
The generated ledger is
`artifacts/generated-results/hc4_univariate_shear_kernel_pencils.json`.

Continue from one transverse-polynomial shear component to a unimodular
pair with:

```bash
.venv/bin/python scripts/verify_hc4_two_component_quasitranslation_kernels.py
```

Together with
[`HC4_TWO_COMPONENT_QUASITRANSLATION_KERNELS.md`](HC4_TWO_COMPONENT_QUASITRANSLATION_KERNELS.md),
this supports `HC4RSD5`. For a fixed primitive kernel
`v=(P,Q,0,0)`, Piola and unimodularity first force `P,Q` to be independent
of the two active variables. Exact mixed-partial integration then gives
`f=x*a(z,w)+y*b(z,w)+G(z,w)` and `P*da+Q*db=0`. The bordered unit makes
`P,Q` algebraically dependent and supplies a polynomial frame of constant
determinant over their closed common composite. Its two differential
coefficients cannot both be nonzero by polynomial degree, so `(P,Q)` lies
on an affine line. A constant active-coordinate change reduces the kernel
to the `HC4RSD4` shear form, whose descendants have triangular polynomial
inverses. The generated ledger is
`artifacts/generated-results/hc4_two_component_quasitranslation_kernels.json`.

Test whether a direct HC4 candidate is covered by an affine singular
scalar pivot with:

```bash
.venv/bin/python scripts/verify_hc4_affine_pivot_coverage_gate.py
```

Together with
[`HC4_AFFINE_PIVOT_COVERAGE_GATE.md`](HC4_AFFINE_PIVOT_COVERAGE_GATE.md),
this supports `HC4RSD6`. For a constant-Hessian potential with Hessian
`H`, a nonzero constant covector `ell` gives an affine singular lift if and
only if `ell^T*adj(H)*ell` is a nonzero constant. On an essential-rank-three
quintic top, the pivot must annihilate the constant top kernel. The next
metric face and the existing Schur face then force a constant relation
`a^T*adj(C)*grad(s3)=0`. Equivalently, all 3-by-3 minors of the degree-eight
coefficient matrix of the Schur vector vanish. The diagonal exact
calibration satisfies the Schur equation but has coefficient rank three,
showing that the coverage gate is a genuine additional restriction. More
generally on that diagonal top, Schur divisibility leaves coefficients
`alpha,beta,gamma`, and affine coverage is confined exactly to
`alpha*beta*gamma=0`. The generated ledger is
`artifacts/generated-results/hc4_affine_pivot_coverage_gate.json`.

Impose the marked collision on every affine-pivot coverage component with:

```bash
.venv/bin/python scripts/verify_hc4_affine_pivot_collision_fibers.py
```

Together with
[`HC4_AFFINE_PIVOT_COLLISION_FIBERS.md`](HC4_AFFINE_PIVOT_COLLISION_FIBERS.md),
this supports `HC4RSD7`. In coordinates adapted to a constant pivot
covector `ell`, the metric numerator `ell^T*adj(Hess(psi))*ell` is, up to a
nonzero constant square, the Hessian determinant of the restriction to
`ell.x=c`. If that numerator is a nonzero constant, every pivot fiber is a
three-variable constant-Hessian potential. `HC3` makes its tangential
gradient injective, so equal full gradients and equal pivot values force
the two points to coincide. Thus an affine zero-corner scalar parent may
represent a four-variable potential but cannot inherit its marked collision
at a common parent pivot value, whether its reduced Hessian is singular or
lies in the nonsingular exact-remainder branch. The generated ledger is
`artifacts/generated-results/hc4_affine_pivot_collision_fibers.json`.

Begin the nonlinear scalar-pivot branch with:

```bash
.venv/bin/python scripts/verify_hc4_quadratic_pivot_rank_obstruction.py
```

Together with
[`HC4_QUADRATIC_PIVOT_RANK_OBSTRUCTION.md`](HC4_QUADRATIC_PIVOT_RANK_OBSTRUCTION.md),
this supports `HC4RSD8`. If `A` is quadratic and
`det Hess(B+s*A)=0`, then a constant nonzero five-variable parent Hessian
determinant forces `rank Hess(A)<=2`. Rank four is excluded by the leading
pencil coefficient. In rank three, splitting off the null direction first
kills `D_z^2 B`; the cleared bordered identity then says that a polynomial
square equals a nonzero constant times `det(s*Q3+H)`, which has degree three
in `s`. The bordered unit makes the affine entries of `grad(A)` generate the
unit ideal, forcing a nonzero linear slice on `ker Hess(A)` and hence the
normal form `A=w+u^T*Qr*u/2`. The checker also gives an exact rank-two
fixed-kernel calibration, so the bound is sharp. The remaining quadratic
frontier is the rank-one and rank-two moving-kernel locus. In rank two, the
passive binary Hessian of `B` is singular. Its rank-zero stratum has
`det(M)=det(D)^2` and reduces to the fixed-support two-component kernel
theorem `HC4RSD5`; only passive rank one is genuinely new. The generated
ledger is
`artifacts/generated-results/hc4_quadratic_pivot_rank_obstruction.json`.

Close the rank-two quadratic branch with:

```bash
.venv/bin/python scripts/verify_hc4_quadratic_rank_two_pivots.py
```

Together with
[`HC4_QUADRATIC_RANK_TWO_PIVOTS.md`](HC4_QUADRATIC_RANK_TWO_PIVOTS.md),
this supports `HC4RSD9`. In the hyperbolic normal form `A=x*y+w`, the
leading pencil and parent faces make `B` affine in the other passive
variable. The next pencil coefficient leaves one active channel, and the
bordered unit makes that channel affine nonconstant. Exact integration gives
`B=x*z+rho*(y+h(x)*A)^2/2+beta(x)*y+gamma(x)*A+delta(x)`.
The checker verifies parent determinant `rho`, descendant determinant
`-kappa*rho`, and the triangular recovery of `x`, the displayed square
coordinate, `A`, `y`, `w`, and `z` from the descendant gradient. Thus every
rank-two quadratic-pivot descendant is a polynomial automorphism. The
generated ledger is
`artifacts/generated-results/hc4_quadratic_rank_two_pivots.json`.

Close the final passive three-by-three branch with:

```bash
.venv/bin/python scripts/verify_hc4_quadratic_rank_one_pivots.py
```

Together with
[`HC4_QUADRATIC_RANK_ONE_PIVOTS.md`](HC4_QUADRATIC_RANK_ONE_PIVOTS.md),
this supports `HC4RSD10`. Normalize the rank-one pivot to
`A=x^2/2+w`. The leading reduced-pencil and parent faces make the passive
three-variable Hessian `E` singular and impose
`a^T*adj(E)*a=0`. Passive ranks zero and two contradict the generic
corank-one bordered unit, so `rank(E)=1`. The rank-one polynomial-Hessian
normal form and the exact identity
`det Hess(Phi)=rho*det(a,d,ell)^2` turn the surviving factor into a unit
frame. Its Wronskian equation fixes the projective direction and gives
`B=x*z+rho*(y+h(x)*w)^2/2+alpha(x)*y+gamma(x)*w+delta(x)`.
The checker verifies the universal block faces, the frame identity, parent
determinant `rho`, descendant determinant `-kappa*rho`, and the triangular
recovery of all four variables. Hence all quadratic scalar pivots in the
identically singular reduced-pencil programme are collision-free. The
generated ledger is
`artifacts/generated-results/hc4_quadratic_rank_one_pivots.json`.

Split and classify the scalar exact-cancellation branch with:

~~~bash
.venv/bin/python scripts/verify_hc4_scalar_cancellation_dichotomy.py
~~~

Together with
[`HC4_SCALAR_CANCELLATION_DICHOTOMY.md`](HC4_SCALAR_CANCELLATION_DICHOTOMY.md),
this supports HC4RSD11--HC4RSD16. For a nonzero pivot corner, completing
the square identifies the parent exactly with a four-variable
constant-Hessian pencil `psi+s*A`, including collisions. For a zero corner,
every graph-coordinate pivot `A=w+q(u)` factors the five-variable gradient
through ternary gradients `grad_u(C+tau*q)`; HC3 makes the parent injective.
The bordered unit forces every quadratic zero-corner pivot into this class.
For nonzero-corner quadratic pencils, rank four is impossible and rank three
has a complete triangular normal form. In rank two, the passive-rank-one
packet reduces to HC2, while passive rank zero is exactly the cotangent lift
of a plane Keller map and hence equivalent to JC2. In rank one, both
ternary singular-Hessian charts reduce to HC2 or the same JC2 packet. The
unit border freezes the active-variable-moving rational chart: in the
constant-kernel type it forces the kernel derivative back into the kernel;
in the exceptional type it fixes the distinguished passive covector or
collapses the form to the first type. Thus every quadratic scalar pencil
reduces to HC2 or exactly JC2. The generated ledger is
`artifacts/generated-results/hc4_scalar_cancellation_dichotomy.json`.

Continue to arbitrary rank-one and leading-rank-three cubic pencil
directions, including the residual tangent-ruling classification, with:

~~~bash
.venv/bin/python scripts/verify_hc4_higher_degree_pencil_obstructions.py
~~~

Together with
[`HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md`](HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md),
this supports HC4RSD17--HC4RSD23 and HC4RSD25--HC4RSD28. It verifies the global equivalence between
a constant-Hessian pencil and a polynomial nilpotent relative Hessian
endomorphism. It then observes that the rank-one proof is degree-independent,
closing every generic-Hessian-rank-one direction to HC2 or JC2. Finally it
checks the unique moving leading-rank-three cubic normal form
`A=w*z+y*b(z)+G(x,z)`: the last three pencil faces force first `C_r=C_x=0`,
then `D_rr=0`, and finally two polynomial units whose derivatives contradict
`b''!=0`. In the constant-kernel residue it verifies the binary null-cone
synchronization and the fixed-cylinder reduction to HC2 or the exact JC2
cotangent lift. It then checks the universal-field identity, every cubic
root chart, all four non-pure binary quartic Schur charts, and all
fourth-power-top correction charts. It continues through the quintic
simple-root square, all four repeated-root Schur ideals, both immutable
exceptional next-face coefficients, and the remaining lower transverse
squares. It also verifies the arbitrary-multiplicity root-valuation
resonance polynomial, enumerates its sextic root weights, and checks every
remaining weighted sextic Schur face including the positive-`z` equations.
Thus homogeneous border coefficients in every degree, arbitrary border
coefficients through degree four, degree-five coefficients with non-pure
leading quintic, and every degree-six coefficient with non-pure leading
sextic are closed. The generated ledger is
`artifacts/generated-results/hc4_higher_degree_pencil_obstructions.json`.

Close the remaining pure-fifth chart, and hence obtain HC4RSD24 and the
complete quintic bordered lemma HC4BL5, with:

~~~bash
.venv/bin/python scripts/verify_hc4_quintic_bordered_lemma.py
~~~

The checker verifies the passive quartic-Hessian first face, the two curved
full lower-tail unit ideals, all passive-affine projective charts, the
nonzero-transverse common-form radical, and the zero-transverse binary
Schur square. Together with HC4RSD23, every border coefficient through
degree five is fixed. The simple-root square in the higher-degree checker
also gives HC4RSD25, closing squarefree leading binary forms in every
degree. Its order-two/order-four double-root comparison gives HC4RSD26 and
closes the generic discriminant stratum. HC4RSD27--HC4RSD28 then close all
non-pure sextic tops.

Stabilize the first passive flag in the remaining pure-sixth chart with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_sextic_collision.py
~~~

This supports HC4RSD29. It verifies the passive-quintic Hessian first face,
the complete next-face factorization, the two discrete collision ratios,
and the three terminal coefficients `1`, `6`, and `-10/9` while retaining
every lower tail capable of entering those faces. Hence the quintic
correction to a pure-sixth top is binary after constant passive coordinates.
The generated ledger is
`artifacts/generated-results/hc4_pure_sextic_collision.json`. A lower
homogeneous component may still break this stabilized direction, so this is
a narrowing theorem rather than a complete degree-six bordered lemma.

Close every pure-sixth chart whose stabilized binary quintic correction has
nonzero passive curvature with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_sextic_lower_flag.py
~~~

This supports HC4RSD30. It verifies the factor equation
`H_yy*Q_1=(P_3,y)^2`, all five repeated-root cubic normal forms, the
terminal resonance coefficients `-1` and `-4`, and the five-step square
cascade eliminating every later transverse tail. The generated ledger is
`artifacts/generated-results/hc4_pure_sextic_lower_flag.json`. The remaining
degree-six scalar boundary has passive-affine quintic correction
`c_5=a*x^5+x^4*L(y,z)`.

Complete the quintic linear form and close every remaining curved quartic
with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_sextic_affine_quartic.py
~~~

This supports HC4RSD31. The checker verifies that
`c4hat=c_4-x^2*L^2/4` has singular passive Hessian, closes both misaligned
curvature ratios, checks the localized aligned resonance ideal `(D,s)`, and
removes the finite, infinite, and zero-Schur lower breaks. The generated
ledger is `artifacts/generated-results/hc4_pure_sextic_affine_quartic.json`.
The surviving degree-six scalar tower is
`c_5=a*x^5+x^4*L` and `c_4=b*x^4+x^3*M+x^2*L^2/4`.

Close that complete two-linear-form tower with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_sextic_two_linear_tower.py
~~~

This supports HC4RSD32. The checker computes the radical packet equations
and both rational charts for the remaining cubic face. For independent
`L,M`, it verifies the terminal coefficients `-4*p`, `-u/9`,
`-v^2/1296`, `-4*r^2`, and the base incompatibility `-1/46656`. It then
checks all dependent-rank normalizations and the global identities
`[J(h+v*D)]_v=D*h_ww*(D*D''-2*(D')^2)` and
`J(h+y*f+z*g)=-(f*g'-f'*g)^2`. The generated ledger is
`artifacts/generated-results/hc4_pure_sextic_two_linear_tower.json`.
Consequently every scalar degree-six leading direction is fixed; the next
scalar degree target begins at repeated-root leading forms of degree seven.

Close every non-pure septic leading form with:

~~~bash
.venv/bin/python scripts/verify_hc4_nonpure_septic.py
~~~

This supports HC4RSD33 and requires Singular. The root-valuation sieve leaves
exactly seventeen faces in transverse degrees four and five. The checker
keeps every same-weight term, notably the binary-linear `z^3` tail in degree
five, and computes exact characteristic-zero saturations on the zero-, one-,
and two-cross-ratio strata. Every coefficient of `g`, `q`, and the cubic
tail lies in the radical. The generated ledger is
`artifacts/generated-results/hc4_nonpure_septic.json`. Only the pure-seventh
chart remains in scalar degree seven; arbitrary repeated-root tops resume in
degree eight.

Open that pure-seventh chart through its first two complete faces with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_septic_opening.py
~~~

This supports the narrowing theorem HC4RSD34. With a fully generic sextic
correction it verifies
`[J]20=49*x^12*det Hess_(y,z)(c6)`. After the passive singular-Hessian
normal form `c6=H6(x,y)+k*x^5*z`, a fully generic quintic correction gives
`[J]19=(49/2)*x^12*(H6)_yy*(2*(c5)_zz-(8/7)*k^2*x^3)`.
Thus the curved chart has
`c5=R5(x,y)+z*P4(x,y)+(2/7)*k^2*x^3*z^2`. The generated ledger is
`artifacts/generated-results/hc4_pure_septic_opening.json`. This does not
close the pure-seventh chart; its degree-eighteen descendants and the
passive-affine sextic-correction boundary remain.

Extract and square-complete that degree-eighteen face with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_septic_degree18.py
~~~

This supports the narrowing theorem HC4RSD35. The checker retains every
quartic coefficient and verifies
`[J]18=x^8*(H_yy*(49*x^4*(c4)_zz-42*k*x^2*P+18*k^2*H-6*k^3*x^5*z)
-(7*x^2*P_y-4*k*H_y)^2)`. It then kills the quartic `z^4` coefficient,
fixes its `z^3` tail to `k^3*x*z^3/49`, and records the remaining binary
divisibility identity. On `k!=0`, its immutable coefficients force
`a6=a5=0`, factor the last curvature ratio into
`(-5*a4*k+14*p4)*(-a4*k+7*p4)`, and then decompose the last two equations
into two generic resonance packets, two double-root packets, and the
`x^3*L` and pure-`x^4` endpoints. The generated ledger is
`artifacts/generated-results/hc4_pure_septic_degree18.json`. The root charts
of the `k=0` identity, the six nonzero-`k` degree-seventeen descendants,
and the passive-affine sextic boundary remain open.

Close every nonzero-`k` descendant with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_septic_moving_closure.py
~~~

This supports HC4RSD36. The four nonzero-`a4` packets have immutable
degree-seventeen coefficients `72*a4^3*k^2/7` or `18*a4^3*k^2/7`; the
rank-drop packet has `24*a3^3*k^2/7`. The pure-`x^4` endpoint leaves two
ratios, killed in degree sixteen by `-24/49` and `256/1225`. The ledger is
`artifacts/generated-results/hc4_pure_septic_moving_closure.json`.

Close the `P_y=0` part of the zero-`k` curved chart with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_septic_kzero.py
~~~

This supports HC4RSD37. It verifies the complete recursive square (11.63),
classifies its nonzero-`p` cubic solutions, checks the normalized
obstructions `-648/49`, `16`, `12/7`, and `48/7`, and closes the zero-`p`
tails by `-4`, `-24`, and the two global linear-coordinate identities. The
ledger is `artifacts/generated-results/hc4_pure_septic_kzero.json`.

Close the `P_y!=0` part with the corrected coupled Wronskian:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_septic_kzero_wronskian.py
~~~

This supports HC4RSD38. The checker retains the cubic `z^3` tail, proves
the exact ordered two-linear-form classification (11.61), and closes all
five projective charts. The ledger is
`artifacts/generated-results/hc4_pure_septic_kzero_wronskian.json`.

Reduce the passive-affine septic boundary with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_septic_passive_affine.py
~~~

This supports HC4RSD39. It verifies the shifted passive-Hessian opening,
closes both misaligned ratios and all five aligned ordered-line charts, and
then reduces the two-linear-form tower to the eight explicit
degree-fifteen packets below (11.74). The ledger is
`artifacts/generated-results/hc4_pure_septic_passive_affine.json`.

Close all eight quartic packets simultaneously with:

~~~bash
.venv/bin/python scripts/verify_hc4_pure_septic_quartic_packets.py
~~~

This supports HC4RSD40. It verifies the common quartic-polar rank-one split,
whose nonzero locus is exactly the square-Hessian resonance
`3*A3^2=8*A2*A4`. The three transverse-direction packets die by `12/7` or
`-12/2401`; the five aligned resonances die in degrees twelve through eight.
The checker also verifies the global affine-transverse identity that makes
the only nonempty zero strata fixed cylinders. The ledger is
`artifacts/generated-results/hc4_pure_septic_quartic_packets.json`.

The same checker now continues into degree five.  On \(x_1=0\) it proves

\[
 D_R=\mathcal F(T,\partial T)-8LN^3S,
 \qquad T=R|_{x_1=0},\quad S=\partial_{x_1}R|_{x_1=0}.
\]

It follows that an \(x_1^2\)-divisible tail cannot repair any quartic graph
1-jet.  The leading quintic faces force
\(T_5=\kappa x_2^5\) and \(\partial_{y_2}T_4=0\).  For the complete v2 trace

\[
 T=\kappa x_2^5+d x_2^3y_1+\rho x_2^2y_2,
\]

exact resultants and first-transverse coefficients give a
characteristic-zero contradiction, even with every allowed off-plane term.
The checker also replays a rational graph that contains the marked collision
and has determinant \(17165601/25\) on all of \(x_1=0\), then detects its
nonzero coefficient \([x_1x_2^7]=22032/125\).  These are HC4MYGJ1,
HC4MYG5J, and HC4MYG5S; the general degree-five graph remains open.

The complementary rational top-cone and relative-linear checks are:

~~~bash
.venv/bin/python \
  scripts/verify_hc4_meng_yang_quintic_graph_normal_slice.py \
  --output \
  artifacts/generated-results/hc4_meng_yang_quintic_graph_normal_slice.json
.venv/bin/python \
  scripts/verify_hc4_meng_yang_quintic_q_kernel_slice.py \
  --output \
  artifacts/generated-results/hc4_meng_yang_quintic_q_kernel_slice.json
.venv/bin/python \
  scripts/verify_hc4_meng_yang_relative_linear_obstruction.py
~~~

The first forces the rational quintic top onto two constant-kernel charts;
the second substitutes the forced plane jet into the first
\(\partial_{y_2}\)-kernel slice with the complete degree-at-most-two lower
trace.  Its generic three immutable transverse equations generate the unit
ideal over \(\mathbb Q\); a square in the forced quartic normal jet kills the
\(y_2^2\) coefficient; and the exceptional nonzero \(x_2y_2\) branch ends in
two coprime transverse polynomials with resultant
\(986335129354383654912000\).  A separate resultant excludes the kernel
denominator chart.  Exact enumeration finds no admissible point modulo 101
or 103, hence no rational point to reconstruct.  As a diagnostic subfamily,
setting the \(y_1y_2\) coefficient to zero leaves two cubics with resultant
\(-108117004020524928=-2^7 3^{17}\cdot11\cdot13\cdot53\cdot863\).
The checker independently replays both transverse branches in the original
five-variable potential.  The third uses the external plane-Jacobian
degree-100 theorem to exclude collisions in every residual-linear correction
of degree at most 89.  The immediate graph target is the joint ideal obtained
from the forced plane normal jet, broader cubic and quartic traces, the
remaining top-cone charts, and first-transverse components above degree
three.  The generated JSON hashes for the projected normal slice and the
transverse kernel slice are,
respectively,
`cdf7fd3cb03dcaea616ce4e177ba87fddcffa423c5f29a0d5f6e4f5dc1e0fee5`
and
`dd28f7a44f6c813bdd422335133869ad7b4a513cc4426987930630c7cba859f9`.

The first bounded mixed canonical-pivot search is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_canonical_pivots.py \
  --output artifacts/generated-results/hc4_mixed_canonical_pivot_search.json
```

It searches 312 exact polynomial symplectic charts generated by quadratic
or cubic Hamiltonians in one mixed source--dual line, two commuting mixed
lines, or a cubic coisotropic-graph generator with one nonlinear constraint.
Pure source transformations are excluded by construction.  All 312 charts
have a scalar affine pivot and there are 258 jointly affine pairs.  Of 4320
specialized scalar-remainder trials, 240 are the exact inherited
`D=0` route and the other 4080 have modular nonconstancy witnesses.  Every
affine pair fails the simultaneous rank-at-most-two budget, and all 41796
complete descended determinants in the declared small repair box have
unequal values modulo `1000003`.  This is the finite-box result `HC4MCP1`,
not an exclusion of symbolic multi-parameter generators, mixed shear
compositions, coefficient-dependent repairs, or general coisotropic
embeddings.  See
[`HC4_MIXED_CANONICAL_PIVOT_SEARCH.md`](HC4_MIXED_CANONICAL_PIVOT_SEARCH.md).

The exact canonical primitives used by that search and by the DC2
Hessian-symbol optimizer now share one convention-explicit implementation.
Check its rank-two/rank-three Poisson signs, mixed-line alphabet, exact word
inverses, and symplectic identities with:

```bash
.venv/bin/python scripts/verify_canonical_transform_search.py
```

The two consumers retain separate Hessian/collision and
symplectic/Weyl/nonsurjectivity gates; see
[`HC4_DC2_CANONICAL_SEARCH_PROTOCOL.md`](HC4_DC2_CANONICAL_SEARCH_PROTOCOL.md).

The first genuinely compositional continuation is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_quadratic_words.py \
  --output artifacts/generated-results/hc4_mixed_quadratic_words.json
```

From 36 signed mixed quadratic Hamiltonian letters it forms all 1296
ordered words, keeps the 648 noncommuting words, removes 48 pure-source
cotangent maps, and deduplicates the remaining 600 exact linear symplectic
charts.  They have 1040 scalar affine pivots and 168 jointly affine pairs.
Exactly 864 scalar trials retain the inherited `D=0` mechanism.  Every pair
fails the rank-at-most-two budget and all 27216 complete determinants in the
same repair box are nonconstant by exact modular witnesses.  This is
`HC4MCP2`; words containing a cubic shear and general coefficient-dependent
repairs remain open.

The first word containing a cubic mixed shear is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_quadratic_words.py \
  --family quadratic-cubic \
  --output artifacts/generated-results/hc4_mixed_quadratic_cubic_words.json
```

It uses the 18 unit-time quadratic and 18 unit-time cubic mixed letters in
both orders.  Exact Poisson-bracket filtering and polynomial-map
deduplication leave 324 noncommuting nonlinear canonical words.  They have
576 scalar affine pivots and 108 jointly affine pairs, but no specialized
scalar remainder survives.  Every pair fails the corank budget and all
17496 complete repairs are nonconstant by exact modular witnesses.  The
post-gate audit also gives a nonconstant parent Hessian determinant in every
chart.  This is the normalized finite-box theorem `HC4MCP3`; signed flow
times, symbolic coefficient families, and coefficient-dependent repairs
remain open.

The fixed-order signed quadratic--cubic canonical box is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_quadratic_words.py \
  --family signed-quadratic-cubic \
  --output \
  artifacts/generated-results/hc4_canonical_signed_quadratic_cubic_words.json
```

It searches \(T_{H_2}\circ T_{H_1}\) with a signed quadratic \(H_1\) and
signed cubic \(H_2\).  Of 1296 raw words, 648 commute and are excluded.
The 648 noncommuting maps are distinct.  Exact support gives affine-pivot
dimension one for 432 words and dimension two for 216; the latter are
exactly the shared-dual words.  Their bracket-incidence census is 96 in
each one-sided type and 24 reciprocal words.  Every transformed reduced
Hessian pencil is generically rank four by an exact modular determinant
witness, every parent Hessian is nonconstant, and all 34992 complete
descended determinants are nonconstant in the declared repair box.  This
is the finite-box result `HC4MCP4`.  The box does not classify oblique
affine directions, symbolic coefficients, longer words, or
coefficient-dependent repairs.

The fixed-order symbolic coefficient closure is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/verify_hc4_symbolic_quadratic_cubic_words.py \
  --output \
  artifacts/generated-results/hc4_symbolic_quadratic_cubic_words.json
```

For each of the 54 noncommuting shared-dual support/sign patterns, it takes
\(H_1=aL_1^2\), \(H_2=bL_2^3\) over `Q[a,b]` and saturates by `a*b`.
Exact parent-Hessian determinant differences at integer probes give 14
localized monomial certificates and 40 unit standard bases in Singular.
Thus no pattern has a parent-constant specialization with `a*b != 0`.
This proves `HC4MCP5`, a coefficient-uniform parent obstruction for the
fixed degree order.  Zero coefficients, other Hamiltonian supports,
oblique directions, and longer words remain open.

The reverse-order symbolic coefficient family is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/verify_hc4_symbolic_quadratic_cubic_words.py \
  --order cubic-quadratic \
  --output \
  artifacts/generated-results/hc4_symbolic_cubic_quadratic_words.json
```

All 54 noncommuting shared-dual patterns now have an exact nonlinear
parent-preserving line.  The 48 one-sided cases force
`a = +/-1/2`; the six reciprocal cases force `a = +/-1/4`; in every case
`b` is arbitrary nonzero and the parent Hessian determinant is identically
`-16384`.  The checker verifies each complete six-variable identity in
Singular.  The retained coordinate affine-pivot dimensions are three in
24 patterns and two in 30 patterns.  All 102 constituent two-pivot pairs
have rank at least three and generic rank four for every `b != 0`, by exact
univariate minor and determinant certificates.  This is `HC4MCP6`.

The next unit quadratic--cubic commutator box is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_commutator_words.py \
  --output \
  artifacts/generated-results/hc4_mixed_quadratic_cubic_commutators.json
```

It searches
`T_-H2 o T_-H1 o T_H2 o T_H1` for the 18 unit quadratic and 18 unit cubic
mixed-line letters.  Of 324 pairs, 162 commute.  The 162 noncommuting
commutator maps are distinct and split into 72 cases in each one-sided
Poisson-incidence type and 18 reciprocal cases.  An exact modular
chain-rule evaluation of the transformed Hessian gives unequal determinant
values at integer points for every word.  Thus every parent determinant is
nonconstant in characteristic zero and there are no survivors requiring
potential expansion.  This finite parent obstruction is `HC4MCP7`.

The smallest reciprocal mixed-line coisotropic graph box is:

```bash
.venv/bin/python \
  scripts/search_hc4_noncoordinate_coisotropic_scalar_gate.py \
  --output \
  artifacts/generated-results/hc4_noncoordinate_coisotropic_scalar_gate.json
```

For ordered `i != j`, it takes
`K=q_i+rho*p_j`, `L=q_j+rho*p_i`, and `H=tau*K*L^2`.
The commuting mixed forms make the time-one flow polynomial and send
`p_i=0` to a nonlinear reciprocal mixed graph.  In the boxes
`rho,tau in {-2,-1,1,2}`, `lambda in {-1,1}`, and
`mu in {-1,0,1}`, the 96 charts have 128 affine scalar pivots.  All 768
graph-specialized scalar Schur remainders have exact nonconstancy witnesses
modulo `1000003`.  The post-gate audit also proves that all 96 parent
Hessian determinants are nonconstant.  Thus no collision-transfer or full
descended-determinant calculation is reached.  This is the bounded result
`HC4MCP8`; arbitrary rational parameters and general coisotropic embeddings
remain open.  See
[`HC4_NONCOORDINATE_COISOTROPIC_GATE.md`](HC4_NONCOORDINATE_COISOTROPIC_GATE.md).

Search the direct affine-unit scalar and unimodular polynomial matrix-pivot
objectives on the nonlinear parent-preserving `HC4MCP6` family with:

```bash
PYTHONHASHSEED=0 PYTHONPATH=scripts .venv/bin/python \
  scripts/search_hc4_nonlinear_unit_schur_blocks.py \
  --output \
  artifacts/generated-results/hc4_nonlinear_unit_schur_blocks.json
```

The 54 reverse-order cubic--quadratic support/sign patterns have exact
parent-constant resonances `a=+/-1/2` or `a=+/-1/4`; the search specializes
`b` to `{-2,-1,1,2}`, giving 216 genuinely nonlinear constant-Hessian
charts.  It exhausts the 364 projective directions in
`({-1,0,1}^6-{0})/{+/-1}`.  None of the 78624 scalar trials has a constant
nonzero second directional derivative.  Exact symbolic filtering leaves
4968 cubic-null directions and 32360 distinct candidate two-planes.  All
but 72 restricted Hessian determinants have modular nonconstancy witnesses;
the 72 reciprocal near misses are jointly quadratic exactly, but each full
determinant is a nonconstant degree-twelve, 28-term polynomial.  Thus the
box contains no scalar unit pivot and no unimodular polynomial `2x2` Schur
block.  This is `HC4MCP9`, not a classification of arbitrary coefficients,
directions, longer words, or polynomial symplectomorphisms.  The generated
artifact hash is
`7f235d427e0cf63e3aeddf198d6ade72c5478ae90774d526fb3a5610dae9286e`.

To remove the finite coefficient and direction boxes on the same HC4MCP6
resonance family, run the exact saturated scalar scheme and all 15 affine
Pluecker charts for each of the 54 patterns:

```bash
PYTHONHASHSEED=0 PYTHONPATH=scripts .venv/bin/python \
  scripts/search_hc4_nonlinear_unit_schur_blocks.py \
  --symbolic-classification \
  --jobs 4 \
  --output \
  artifacts/generated-results/hc4_symbolic_unit_schur_classification.json
```

This requires Singular.  The number of parallel pattern workers does not
change the sorted artifact.  Every sampled determinant certificate is exact:
a constant determinant would satisfy all sampled equalities, so a unit
sampled ideal excludes the complete constant locus without interpolation.
The generated artifact SHA-256 is
`e92465e4991e7635f07fcc70895995f5d0465a1c3b816a6c9a88643500865e30`.

The direct one-variable calculation for the `PC(2)` graph is:

```bash
.venv/bin/python scripts/verify_hc4_direct_schur_descent.py
```

It first verifies that the nonsymmetric `PC(2)` Jacobian cannot literally be
a Hessian Schur complement.  It then classifies every coordinate chart and
omitted graph coordinate for an arbitrary polynomial auxiliary function of
`(X,Y,W,D)`, under transfer of a pair from the certified fiber.  Component
ideals, explicit fixed points, and a two-variable Jacobian-mate lemma leave
only charts `0010` and `0011`; in both, every slice is `W` up to scaling and
retained-coordinate gauge.  Their canonical generating families retain the
rational collision.  Irreducibility reduces every polynomial quadratic pivot
to three cases, and exact Hessian evaluations exclude constant nonzero
determinant in all six.  Finally, the checker proves that the Meng--Yang
five-variable potential has no further polynomial partial Legendre transform
along a constant linear direction.  Non-coordinate graph embeddings and
mixed coefficient-dependent critical equations are not tested.  Pure
univariate higher-degree repairs are excluded abstractly by `SDX1`.

The double-Schur audit for the parameterized quadratic-gauge families is:

```bash
.venv/bin/python scripts/verify_hc4_double_schur_gauge_obstruction.py
```

It tests two exact all-parameter routes from the six-variable Meng doubling
to four variables.  For elimination of two repaired dual variables, the
coefficient of the square of the retained dual variable is the bordered
Hessian invariant `K(L)`.  Three cubic coefficient equations and an
all-degree leading-layer calculation prove `K(L)` is nonzero for every
nonzero constant target linear form throughout the root-engineered gauge
family.  For a source-first descent, the common first coordinate forces the
only possible affine source direction to be `z`; higher-degree decorations
destroy that direction, while the cubic coefficient row has no constant
nonzero second pivot.  The calculation does not test nonlinear symplectic
changes, nonlinear retained dual coefficients, or nonconstant pivots with
exceptional divisibility.

The first nonlinear continuation, using triangular target coordinates, is:

```bash
.venv/bin/python scripts/verify_hc4_triangular_target_shears.py
```

After the canonical linear normalization of the full cubic gauge family,
it writes `L=W+H(U,V)` for each target permutation and a general polynomial
`H` of total degree at most three.  Sparse exact coefficients of the
bordered invariant `K(L)` first force the cubic part of `H` to vanish and
then give the quadratic contradiction.  In the third orientation the
coefficient `[x^8]K=9` is independent of every shear parameter.  Thus no
candidate reaches the full four-variable Hessian check.  The length-two
quadratic composition is treated next; degree at least four and general
nonlinear target or source--dual symplectic changes are not tested.

The length-two quadratic continuation is:

```bash
.venv/bin/python scripts/verify_hc4_two_quadratic_target_shears.py
```

For each of the six ordered words it writes
`A=X_i+Q(X_j,X_k)` and `L=X_j+R(A,X_k)`, with general positive-degree
quadratics `Q,R`.  Exact axis and transverse-line coefficients exclude the
`A^2` and `A*X_k` coupling strata over the reals.  If both couplings vanish,
the word is one quadratic triangular shear and is covered by the preceding
cubic single-shear checker.  This does not test degree-four single shears,
words with higher-degree factors, length at least three, general target
automorphisms, or nonlinear source--dual symplectic changes.

The preliminary bounded modular screen is:

```bash
.venv/bin/python scripts/search_hc4_two_quadratic_target_shears.py --prime 3
```

That command is an experiment, not a proof; the exact checker above is the
status-bearing calculation.

The degree-four single-shear closure is:

```bash
.venv/bin/python scripts/verify_hc4_quartic_target_shears.py
```

It includes every lower-degree shear coefficient and uses five extreme
spatial coefficients to force the homogeneous quartic layer to vanish in
the `P`- and `B`-retained orientations.  The cubic checker then applies.  In
the `C`-retained orientation, `[x^8]K=9` remains independent of every
coefficient through degree four.  This exact characteristic-zero result
does not test degree at least five, words with cubic-or-higher factors,
length at least three, general target automorphisms, or nonlinear
source--dual symplectic changes.

The all-order binary-form symbol underlying the two \(X\)-caustic normal-jet
recursions is documented in
[`HC4_PC2_GRAPH_POLARIZATION_AUDIT.md`](HC4_PC2_GRAPH_POLARIZATION_AUDIT.md).
Its finite exact regression in both charts is:

```bash
.venv/bin/python scripts/verify_hc4_logarithmic_normal_symbol.py
```

The checker verifies through tensor order ten that the order-\(n\) symbol is
\(L\) times a unimodular coefficient-extraction matrix.  The all-order
statement is proved by the triangular determinant formula in the note.  This
establishes a divisor-supported graded jet module, not yet a finite-rank
logarithmic connection or a Bernstein--Sato polynomial.

The nonlinear toric descent of the Meng--Yang potential is:

```bash
.venv/bin/python scripts/verify_hc5_nonlinear_toric_descent.py
```

It constructs an explicit determinant-one dual-coordinate change with
`t=A`, verifies the polynomial unit-pivot critical solution, and computes the
natural four-variable determinant `16*J(x1*x2)^2`.  A relative toric `SL(2)`
correction cancels the nonconstant factor and gives determinant `64`, but the
descended gradient then has an explicit polynomial inverse.  The checker
also proves the all-degree toric radial obstruction: constant nonzero
determinant forces both radial factors to be units, so this class cannot
retain the Meng--Yang collision.  Non-toric changes and non-coordinate
coisotropic embeddings remain open.

The bounded non-toric relative correction is:

```bash
.venv/bin/python scripts/verify_hc4_nontoric_sl2_correction_degree4.py
```

For a general matrix `C(x,y) in SL(2)` whose four entries have total degree
at most four, it sets `G=beta*C` for the natural complementary coefficient
row.  The four-variable identity is
`det Hess(psi)=16*det(DG)^2`.  Singular proves that the 45 determinant-one
equations and 218 nonconstant-Jacobian equations generate the unit ideal
over `QQ`; collision equality is not needed.  A second exact calculation
excludes arbitrary affine `SL(2)` perturbations of the known degree-ten
toric correction after the collision equations are imposed.  Raw degree at
least five, quadratic-or-higher perturbations of the toric correction,
mixed base--dual changes, and non-coordinate coisotropic embeddings remain
open.  This checker requires `Singular` on `PATH`.

The shortest exact replay of the canonical Meng descent chain
`HC5T1 -> HC4MQ1 -> HC4MCK` is:

```bash
.venv/bin/python scripts/verify_hc5_nonlinear_toric_descent.py
.venv/bin/python scripts/verify_hc4_meng_sparse_quartic_obstruction.py
.venv/bin/python scripts/verify_hc4_meng_full_cubic_kernel.py
```

The final command requires `Singular` on `PATH`.  The support-three and
support-four commands below are exact historical checkpoints and targeted
regressions, not logical prerequisites for `HC4MCK`.

The first collision-first non-toric Hamiltonian screen is:

```bash
.venv/bin/python scripts/verify_hc4_meng_sparse_quartic_obstruction.py
```

After the unit-pivot descent and a polynomial base gauge, the transported
points are antipodal for `psi_0=2*y*r+4*x*s`.  The checker exhausts every
homogeneous quartic correction supported on at most four monomials, imposing
the collision before determinant work.  It rejects 42,953 isolated collision
solutions and 515 one-parameter families exactly.  The degree-eight
principal part then leaves only 232 isolated quartics and two exact family
members; none admits a constant-determinant correction by one cubic monomial.
For two cubic monomials, degrees seven and one give a linear rank gate;
rank-zero cases pass through degree-six conic linearization and a
degree-four/two lift.  Only four bivariate families reach the terminal
calculation, and their full determinant ideals are units modulo `1000003`.
The complete homogeneous-quartic reduction is:

```bash
.venv/bin/python scripts/verify_hc4_meng_full_quartic_reduction.py
```

It verifies an explicit complex congruence taking
`2*y*r+4*x*s` to one half the sum of four squares, checks the gradient
chain rule on all 35 quartic monomials, and transports the antipodal
collision.  The transformed gradient is `z+grad(h4)`, with homogeneous
cubic nonlinear part and symmetric Jacobian.  The external
de Bondt--van den Essen dimension-four theorem therefore makes every
constant-Jacobian member a polynomial automorphism.  This closes the
complete homogeneous quartic chart, and the same congruence argument closes
every nondegenerate quadratic renormalization plus one homogeneous quartic.
Mixed nonlinear degrees are not covered.

The complete dense cubic--quartic continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_dense_cubic_quartic_reduction.py
```

For `psi=q2+h3+h4`, homogeneous determinant layers and
Gordan--Noether give a constant direction satisfying
`D_v(h4)=D_v^2(h3)=0` in every rank of `Hess(h4)`.  The rank-one quartic
case uses an additional exact ternary Hessian-pencil calculation.  A
nonisotropic direction descends to `HC(3)`; the isotropic bordered form
reduces either to a degree-at-most-three plane Keller map or to an `HC(2)`
block after a binary quadratic invariant forces `s=constant+x`.  This
excludes all 20 cubic and all 35 quartic coefficients simultaneously,
without a support bound.  The external inputs are Gordan--Noether,
`HC(3)`, `HC(2)`, and Moh's plane degree bound.

The independent dense mixed-quartic coefficient regression is:

```bash
.venv/bin/python scripts/verify_hc4_meng_dense_mixed_quartic.py
```

It treats all 25 homogeneous quartic monomials of base--dual bidegrees
`(1,3)`, `(2,2)`, and `(3,1)`, then separately adjoins every pure-base
quartic and every pure-dual quartic.  Collision and the linear determinant
layer have rank 14, leaving 11, 16, and 16 parameters.  Singular proves that
the remaining exact coefficient ideals, with 262, 273, and 273 generators,
are unit ideals over `QQ`.  The complete homogeneous-quartic theorem
subsumes this one-sided coefficient calculation, and the dense
cubic--quartic theorem subsumes the later sparse cubic-kernel chain.
Quartic--sextic, simultaneous cubic--quartic--sextic, and non-coordinate
embeddings remain open.  This command requires `Singular` on `PATH`.

The finite-field continuation through exactly three cubic monomials is:

```bash
.venv/bin/python scripts/verify_hc4_meng_three_cubic_rank_gate.py
```

It checks all `234*binomial(20,3)=266760` quartic/triple pairs.  The combined
degree-seven/degree-one rank census is `5480, 53364, 130508, 77408` in ranks
zero through three.  After delegating support-at-most-two boundary loci to
the preceding certificate, unit full-determinant gcds or ideals exclude all
920 genuine rank-two lines, 2,952 genuine rank-one planes, and 5,480
rank-zero three-parameter spaces over `F_1000003`.  This is an exact
finite-field computation.

The characteristic-zero promotion is:

```bash
.venv/bin/python scripts/verify_hc4_meng_three_cubic_characteristic_zero.py
```

It reconstructs all 234 quartics over `QQ`, reproduces the same odd-layer
rank census over `QQ`, and obtains unit gcds or Gröbner ideals for every
genuine line, plane, and three-parameter space.  The maximum evaluation
prefix lengths are respectively 5, 7, and 9 points.  Thus cubic support at
most three is excluded in characteristic zero.  This bounded-support
checkpoint is subsumed by the full cubic-kernel theorem below.

The finite-field continuation through exactly four cubic monomials is:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic.py
```

Its two targeted stages are:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic_rank_gate.py
.venv/bin/python scripts/verify_hc4_meng_four_cubic_rank_zero.py
```

The first checker exhausts all
`234*binomial(20,4)=1133730` quartic/quadruple pairs.  The odd-layer ranks
zero through four are `5430, 79396, 353740, 504818, 190346`.  After
delegating support-at-most-three boundaries to `HC4MC3`, unit gcds or ideals
exclude all 466 genuine rank-three lines, 6,082 genuine rank-two planes, and
7,956 genuine rank-one three-spaces.  The second checker isolates the 5,430
rank-zero four-parameter spaces and proves that every full determinant ideal
is a unit modulo `1000003`, using at most twelve evaluation points.

The characteristic-zero promotion is:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic_characteristic_zero.py
```

It reconstructs the same rank and boundary census over `QQ`.  Exact rational
gcds or Gröbner ideals exclude all 466 lines, 6,082 planes, 7,956
three-spaces, and 5,430 four-spaces.  The maximum evaluation-prefix lengths
for ranks three through zero are `5, 7, 8, 11`.  Thus cubic support at most
four is excluded in characteristic zero.  This bounded-support checkpoint is
also subsumed by the full cubic-kernel theorem below.

The full cubic-kernel characteristic-zero checker is:

```bash
.venv/bin/python scripts/verify_hc4_meng_full_cubic_kernel.py
```

It requires `Singular` on `PATH`.  For all 234 rational quartics, it
constructs the complete odd determinant kernel in the 20-dimensional cubic
space, symbolically extracts every spatial coefficient of the Hessian
determinant, clears denominators, and adds the coefficient ideals in
descending spatial degree.  The exact odd ranks are 8 for 229 quartics, 7
for one, and 4 for four exceptional `u^3*L` quartics, giving kernel
dimensions 12, 13, and 16.  Singular reaches the unit ideal first at degree
six for 62 quartics, degree five for 16, and degree four for 156.  Therefore
arbitrary homogeneous cubic corrections are excluded in characteristic
zero whenever the collision quartic has support at most four.

The parallel sparse sextic collision-carrier checker is:

```bash
.venv/bin/python scripts/verify_hc4_meng_sparse_sextic_obstruction.py
```

It exhausts all supports of at most four among the 84 homogeneous sextic
monomials.  Collision-first linear algebra leaves 1,725,838 isolated points,
7,566 lines, and one plane.  Five leading-Hessian evaluations reduce the
isolated points to 748 candidates, all rejected by the first two full
determinant evaluations.  Exact rational gcds leave only two principal
lines, at parameters `-81/8` and `9/2`; determinant degree twelve rejects
both.  The unique binary-cubic plane in `y*r` and `x*s` has unit principal
coefficient ideal over `QQ`.  Thus no sextic-only correction supported on
at most four monomials works.  The later joint theorem `HC4JQS4` treats
mixed quartic--sextic corrections of combined support at most four, and
the later support-free theorem `HC4E46` subsumes both calculations.

The characteristic-zero continuation over the 234 quartic principal parts
is:

```bash
.venv/bin/python scripts/verify_hc4_meng_mixed_quartic_sextic.py
```

Exactly four of the 234 quartic principal parts have zero immutable
determinant-degree-two signature.  In their homogeneous zero-gradient
sextic kernel, the checker finds 976 three-support lines, 205,494
four-support lines, and 519 planes.  Principal cancellation leaves 121,146
lines.  For each quartic, 52,686 lines have projectively complete modular
certificates and the remaining 68,460 have unit exact rational gcds.
Separate exact affine-matrix expansion proves that all `519*4=2076`
three-evaluation plane ideals are units over `QQ`.  Thus no zero-gradient
sextic supported on at most four monomials repairs any of the 234 quartic
principal parts.  The later `HC4JQS4` theorem treats quartics outside that
sextic-free principal screen when the combined support is at most four;
the larger mixed chart and simultaneous cubic corrections remain open.

The joint total-support-four continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_joint_quartic_sextic_total_support_four.py
```

It exhausts all 6,133,820 genuinely mixed quartic--sextic supports of
combined size at most four.  Collision-first linear algebra leaves
5,225,684 isolated points, 44,300 lines, and 34 planes.  Descending
determinant-layer evaluations reject every isolated point modulo the
rank-preserving prime `1000003`.  Every line is reconstructed over `QQ`
and has unit exact evaluation gcd; every plane has unit exact
three-evaluation ideal.  Together with the pure quartic theorem `HC4MQ1`
and pure sextic theorem `HC4MS6`, this proves `HC4JQS4`: no
quartic-plus-sextic correction of total support at most four works when
there is no cubic term.  The later `HC4E46` theorem removes the support
bound completely; simultaneous cubic and sextic corrections remain open.

The support-free common-kernel continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_dense_rank_three_sextic_reduction.py
```

It checks the determinant identities in theorem `HC4DCK`.  The theorem
excludes every support-free correction having a constant direction \(v\)
with `D_v h6=0` and `D_v^2 h4=0`.  Gordan--Noether and the
degree-fourteen determinant layer make this automatic for a rank-three
sextic Hessian.  A non-isotropic quadratic pivot reduces to `HC_3`.  In the
isotropic bordered chart, the next two coefficients reduce either to a
degree-at-most-five plane Keller map or to a two-variable constant-Hessian
block.  Moh's degree bound and `HC_2` exclude both.  Only sextic Hessian
rank at most two with a variable quartic null direction is left by
`HC4DCK` itself; theorem `HC4E46` below closes it.  Simultaneous cubic and
sextic corrections remain open.

The source/dual reorganization of that rank-two boundary is:

```bash
.venv/bin/python scripts/verify_hc4_source_dual_bigrading.py
```

It checks the complete \(2+2\) bidegree ledger, the corrected weighted
Hessian-face identity, and the cotangent determinant
`det Hess(t*F+m*G+H)=Jac(F,G)^2`.  The canonical note proves theorem
`HC4SDW`: the dual-linear stratum is exactly a degree-at-most-five
`JC(2)` locus, and every vanishing sequence of rank-two residual Hessian
faces synchronizes by a Schur recursion to one rational projective cone.
The rotating example `(x*t+y*m)^2` is synchronized but nonconstant, so
the synchronization theorem alone leaves a nonlinear moving-cone
algebraization problem.

The support-free closure of the full even quartic--sextic chart is:

```bash
.venv/bin/python scripts/verify_hc4_even_quartic_sextic_closure.py
```

The primitive cone-degree lemma shows that bihomogeneity leaves only
`c*(X^T*M*U)^2` as a nonconstant rank-two residual cone.  A rank-one `M`
has a constant dual direction and is `HC4DCK`; for invertible `M`, the
dual-degree-four part of the spatial `z^4` determinant coefficient is

```text
48*c^4*det(M)^2*(X^T*M*U)^4,
```

and cannot be cancelled by the quadratic block, a dual-linear/source-only
quartic block, or the source-only sextic block.  The rank-one sextic
boundary becomes constant by one-base scaling, rank three is `HC4DCK`,
and rank zero is `HC4HQ1`.  This proves `HC4E46`: no support-free
potential `q2+h4+h6` with nondegenerate `q2` has both constant nonzero
Hessian determinant and an antipodal collision.  Only simultaneous cubic
and sextic interaction, and non-coordinate coisotropic embeddings, remain.

The rank-three part of the simultaneous cubic--quartic--sextic chart is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_three_reduction.py
```

For `psi=q2+h3+h4+h6`, the degree-fourteen and degree-thirteen determinant
layers on the one-dimensional sextic kernel force both
`D_t^2(h4)=0` and `D_t^2(h3)=0`.  The nonisotropic direction descends to
`HC(3)`.  In the isotropic chart, the checker reconstructs the complete
cubic bordered invariant.  The two rank-two binary-cubic orbits have
coefficient ideal `(qxm,qym,lm)^2`; the rank-one orbit has two explicit
radical branches, each with a constant missing direction.  The remaining
descent is either a plane Keller cotangent lift of degrees at most three
and five or an `HC(2)` block.  This proves `HC4T31` without a support
restriction.  Only simultaneous corrections with
`rank Hess(h6)<=2`, and non-coordinate coisotropic embeddings, remain.

The rank-two continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_two_reduction.py
```

The degree-twelve face makes the binary quartic Hessian on the constant
sextic kernel plane singular.  A constant cone immediately recovers the
common-direction reduction; if that binary Hessian is zero, the
degree-ten face makes the cubic binary Hessian singular, and total degree
three forbids a moving cone.  In the only moving quartic case, the
degree-eleven face forces the high-dual cubic to align as
`(X^T*M*U)*(alpha*t+beta*m)`.  The checker keeps this cubic and all
compatible lower blocks and proves that the later dual-degree-four face
is still the uncancellable
`48*c^4*det(M)^2*(X^T*M*U)^4`.  This proves `HC4T21`.

The rank-one continuation and complete coordinate-chart exhaustion are:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_one_reduction.py
```

The small-rank Hessian normal form makes the sextic a sixth power
`c*L^6`, with constant three-dimensional kernel `W`.  The degree-ten
face makes `Hess_W(h4)` singular.  Its rank-two branch aligns the cubic
at degree nine.  In rank one, the degree-eight face is exactly the binary
discriminant `a*d-b^2` of `Hess(h3)` on the constant quartic-kernel
plane; this is the `Sym^2` nullcone `L^2` from the SIC(2) binary-root
classification.  In rank zero, degree seven makes `Hess_W(h3)` singular.
One-base homogeneity makes every residual projective cone constant, so
each branch reaches the common-direction descent `HC4T31`.  This proves
`HC4T11`; ranks three, two, one, and zero then give `HC4TC1`, the complete
support-free obstruction for `q2+h3+h4+h6`.  Quintic and higher
homogeneous layers, as well as non-coordinate coisotropic embeddings,
remain outside this chart.

The first quintic extension of the common-direction descent is:

```bash
.venv/bin/python scripts/verify_hc4_quintic_common_direction.py
```

The checker proves the quartic bordered lemma `HC4BL4` by the binary-root
stratification used in SIC(2), including exact radical certificates for
the double-double and pure-fourth exceptional strata.  It then proves
`HC4CD5`: for `psi=q2+h3+h4+h5+h6`, a common direction satisfying
`D_v h6=0` and `D_v^2 h5=D_v^2 h4=D_v^2 h3=0` reduces by two Schur
steps to `HC(3)`, Moh's plane theorem, or `HC(2)`.  In sextic Hessian rank
three this closes the branch `D_v h5=0`.  The first remaining quintic
face is the exact divisibility

```text
det(Cbar) | grad(D_v h5)^T * adj(Cbar) * grad(D_v h5),
```

where `Cbar` is the nondegenerate ternary block of `Hess(h6)`.

The full Fermat-sextic diagonal Schur-norm stratum is checked by:

```bash
.venv/bin/python scripts/verify_hc4_quintic_diagonal_schur.py
cd formal/finite-etale-keller
lake env lean FiniteEtaleKeller/HC4QuinticDiagonal.lean
```

The first command retains an arbitrary base quintic and every lower
quartic, cubic, and quadratic coefficient.  Exact truncated determinant
expansion shows that the `lambda^13*t*x^4*y^4*z^4` face fixes the `t^3`
coefficient of `h3`, after which three `lambda^11*t^3` coefficients are
`1024*a^5`, `1024*b^5`, and `1024*c^5`.  Before that prolongation, an
exact 66-equation radical certificate proves that every quartic with
polynomial Fermat Schur norm is diagonal.  Thus `HC4QF1` closes the entire
Fermat-sextic quintic stratum.  Lean checks
the scalar Schur identity and the characteristic-zero fifth-power
conclusion; determinant coefficient extraction remains the Python
certificate.

The first full non-diagonal sextic pencil is:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_symmetric_sextic_schur.py
```

For `h6=(x^6+y^6+z^6)/30+mu*x^2*y^2*z^2`, the checker constructs the
complete 111-equation Schur-divisibility ideal for a generic 15-coefficient
quartic and six-coefficient quadratic quotient.  Saturation by `mu` has an
exact 261-element rational Gröbner basis.  The saturated ideal lies in the
21-variable coefficient origin, and the fourth powers of all 21 generators
reduce to zero.  Thus `HC4QS1` closes every `mu!=0` member; `HC4QF1` closes
`mu=0`.  The fourth-power reduced-ring endpoint is replayed in
`FiniteEtaleKeller/HC4QuinticDiagonal.lean`.

The two-parameter generic broadening is:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py
```

For
`h6=(x^6+y^6+z^6)/30+mu*x^2*y^2*z^2+nu*sum_(i!=j)x_i^4*x_j^2`,
six quotient pivots have determinant `4096*nu^12`.  After eliminating the
quadratic quotient and clearing `2*nu^2`, 114 intrinsic equations remain.
Over `Q(mu,nu)` their exact 117-element Gröbner basis contains the cubes of
all fifteen quartic coefficients.  This proves `HC4QSG2`, generic rigidity
on the two-parameter surface.  It does not exclude exceptional curves
inside `nu!=0`: exact uniform saturation timed out at 300 seconds, and
raw projective-chart runs timed out at 120 seconds.  Those timeouts are
diagnostics, not certificates.

The surface nevertheless contains an explicit exceptional Schur pair:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_radial_exceptional_schur.py
```

At `(mu,nu)=(1/5,1/10)`, the sextic is
`(x^2+y^2+z^2)^3/30`.  The checker proves that the nonzero quartic
`s4=(x^2+y^2+z^2)^2` has polynomial Schur quotient
`16*(x^2+y^2+z^2)` and that the sextic Hessian determinant is
`(x^2+y^2+z^2)^6/25`.  This verifies only the Schur face; it does not
by itself claim extension through the lower collision identities.  In
the invariant coordinates `R=x^2+y^2+z^2`,
`P2=x^2*y^2+x^2*z^2+y^2*z^2`, and `P3=x^2*y^2*z^2`, the surface is
`R^3/30+(nu-1/10)*R*P2+(mu-3*nu+1/10)*P3`.  The checker also proves that
the fixed radial quartic `R^2` has a polynomial Schur norm only when both
deformation coefficients vanish, so this pair occurs only at the stated
parameter point.

The lower-face prolongation is excluded by:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_radial_prolongation.py
```

The checker retains an arbitrary base quintic and every allowed quartic,
cubic, and quadratic coefficient.  Sparse exact determinant extraction
gives `2*(3*delta-32)/25` at `lambda^13*t*x^12`, forcing
`delta=32/3`, while `lambda^11*t^3*x^8` then equals `1024/25`.
Every arbitrary lower coefficient cancels.  Thus the radial Schur pair is
exceptional at the Schur face but cannot produce a constant-Hessian
collision.

The full even permutation-invariant quartic line is classified by:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_even_symmetric_schur.py
```

For `s4=a*R^2+b*P2`, the `a=1` chart has radical support at exactly two
points: the radial pair and the Fermat pair
`(mu,nu,s4)=(0,0,x^4+y^4+z^4)`.  The `a=0` chart is the unit ideal.
Consequently this complete symmetric-even slice contains no exceptional
parameter curve, and both isolated points are already excluded at lower
faces.

The component-directed exceptional-locus research transcript is:

```bash
.venv/bin/python scripts/verify_hc4_exceptional_schur_atlas.py
.venv/bin/python scripts/research_hc4_exceptional_schur_locus.py
.venv/bin/python scripts/research_hc4_exceptional_schur_locus.py \
  --exact-pure-chart --singular-timeout 900
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --extract-basis-denominators --basis-profile
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 0 --cube-torsion-stage generic
```

The first command is the exact full reduced atlas `HC4QSE4`.  It builds
the 120-equation projective incidence system, expresses the Hessian
determinant in the seven symmetric invariants, and proves that its
nonsquarefree locus consists only of Fermat and radial.  The apparent
third discriminant point `(5/7,-1/14)` has two distinct coprime Hessian
factors.  Hence the reduced incidence is the Fermat projective plane of
diagonal quartics plus the single radial quartic `R^2`.  The lower-face
checkers `HC4QF1` and `HC4QSE2` exclude both before additional antipodal
collision equations can contribute.

The later commands are the earlier modular reconstruction,
special-fiber geometry, and denominator research recorded in
[`HC4_EXCEPTIONAL_SCHUR_LOCUS.md`](HC4_EXCEPTIONAL_SCHUR_LOCUS.md) and
[`hc4_exceptional_schur_locus_modular.json`](artifacts/generated-results/hc4_exceptional_schur_locus_modular.json).
The even-quartic charts reconstruct the Fermat and radial parameter
points modulo \(47,101,103\).  The transformation-aware
`--extract-denominators` calculation timed out at its 900-second Singular
bound and supplies no certificate.  The basis-only denominator is the
constant \(2\), so the desired exceptional divisor must occur in the lift
certificates.  The displayed first sign-character-block cube lift also
timed out at 900 seconds; it has 191 target cubic monomials and 441
multiplication columns.  The direct characteristic-zero `a=1`
even-quartic elimination likewise reached its 900-second bound before
returning a standard basis.  Those interrupted routes are historical
diagnostics and are not used by `HC4QSE4`.

The degree-three cube-torsion research modes are:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 2 \
  --cube-torsion-stage finite-field --cube-prime 19
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 2 \
  --cube-torsion-stage specialize \
  --cube-mu-value=-5/3 --cube-nu-value=-1/6
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 2 \
  --cube-torsion-stage fiber \
  --cube-mu-value=-5/3 --cube-nu-value=-1/6
```

The construction uses the canonical multiplication presentation
`A^1710 -> A^680` and the fifteen cube targets, not a noncanonical
`A^114 -> A^15` compression.  Complete scans of the four
coefficient-monomial orbits modulo `11,13,17,19` reconstruct the radial
point and the additional point `(-5/3,-1/6)` on `nu!=0`.  Exact rational
specialization shows that precisely the three `x_i^2*x_j^2` coefficient
cubes survive at the new point, while every fourth power is zero; the
60-dimensional fiber is therefore supported at the coefficient origin.
This is a nilpotence-order jump, not a reduced exceptional Schur pair.
The even-block integral annihilator and function-field lift each reached a
900-second timeout; relation extraction and the zeroth Fitting ideal were
not reached.  See
[`HC4_FITTING_DENOMINATOR_EXTRACTION.md`](HC4_FITTING_DENOMINATOR_EXTRACTION.md)
and
[`hc4_fitting_denominator_extraction.json`](artifacts/generated-results/hc4_fitting_denominator_extraction.json).

The full-coefficient reduced-fiber scans are:

```bash
for prime in 7 11 13; do
  .venv/bin/python -u \
    scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
    --fourth-power-profile --fourth-prime "$prime"
done
```

The corresponding symbolic attempt is:

```bash
.venv/bin/python -u \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --fourth-power-profile --fourth-prime 7 \
  --fourth-stage annihilator --fourth-timeout 900
```

The zero sign-character block of degree four has 819 target monomials and
3474 equation-times-quadratic columns.  At all \(42+110+156\) parameter
points on `nu!=0`, membership of all fifteen coefficient fourth powers
certifies an empty reduced projective Schur fiber except at the radial
reduction.  The cube-torsion point `(-5/3,-1/6)` is certified
reduced-empty at every prime.  These scans cover \(\mathbb F_p\)-rational
parameters, not points over proper finite-field extensions and not the
characteristic-zero support ideal.  The symbolic \(\mathbb F_7[\mu,\nu]\)
annihilator command reached its 900-second bound before returning a
standard basis and supplies no factorization.  See
[`hc4_fourth_power_support.json`](artifacts/generated-results/hc4_fourth_power_support.json).

The direct collision-normalized finite-field experiment in degree bounds
five through eight is:

```bash
.venv/bin/python scripts/search_hc4_finite_field_potentials.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-bound 2 \
  --points 6 \
  --output artifacts/generated-results/hc4_finite_field_sparse_search.json
```

It exhausts one- and two-direction affine perturbations in the full linear
kernel of the normalized gradient-collision condition, totaling 45,181,194
coefficient choices.  No exact modular candidate survives.  This is a
bounded experiment, not evidence for unrestricted `HC_4`; its construction
and precise scope are in
[`HC4_FINITE_FIELD_SEARCH.md`](HC4_FINITE_FIELD_SEARCH.md).

The denser sampled-support coefficient solve is:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies uniform homogeneous mixed \
  --trials 2 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_dense_support_search.json
```

It forms the complete determinant coefficient ideal for 96 deterministic
supports of 6, 8, 10, or 12 collision-kernel directions.  Adjoining
`a_i^p-a_i` makes each Singular calculation an exact existence test over the
selected prime field.  All 192 support-prime ideals are unit ideals, with no
timeouts.  The support selection is sampled, so this remains a bounded
experiment.

The structurally guided companion forces every sampled support to contain as
many directions as possible that can alter both base determinant defects on
the normalized collision axis:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies axis \
  --trials 2 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_axis_support_search.json
```

All 64 ideals from these 32 additional supports are also unit ideals.

The principal-Hessian companion makes every degree-\(d\) correction a
three-variable cone, so its top homogeneous Hessian determinant vanishes
identically, then adds lower-degree monomials involving the omitted fourth
variable:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies cone2 cone3 \
  --trials 4 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_cone_bridge_search.json
```

All 256 full coefficient ideals from these 128 additional supports are unit
ideals.  Thus the failure occurs below the principal homogeneous Hessian
gate, but the support search remains bounded.

The non-coordinate cone search uses
`u=x2+lambda*x3`, `v=x3` for `lambda=-1,1,2`.  In the adapted coordinates
the quadratic base is `x0*x1+u*v-lambda*v^2`, the top correction depends only
on `x0,x1,u`, and the lower bridges involve `v`:

```bash
.venv/bin/python scripts/search_hc4_oblique_cone_bridges.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --slopes -1 1 2 \
  --trials 3 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_oblique_cone_bridge_search.json
```

All 288 exact ideals from the 144 oblique families are unit ideals, with no
timeout.  The adapted change has determinant one, so the Hessian determinant
identity is equivalent to the tied-monomial identity in the original
coordinates.

## Proof-carrying arithmetic compiler

The active common-arithmetic-fibers paper has a separate correspondence
compiler for its displayed Berend--Bilu example. One compact JSON
specification generates the Lean-readable coefficients and map, TeX macros
used directly by the manuscript, runnable SymPy and Sage inputs, and an
expanded sparse JSON certificate. The verifier checks that every generated
view is current, runs the SymPy input, and builds the Lean theorems equating
the generated paper polynomial, map, inverse polynomial, output scalings, and
distinguished targets with the existing formal definitions:

```bash
make verify-common-arithmetic-fibers-correspondence
```

Intentional regeneration is:

```bash
make refresh-common-arithmetic-fibers-example
```

The canonical source is
`papers/common-arithmetic-fibers/data/explicit-quintic-spec.json`; generated
results are recorded under `artifacts/generated-results/`, with the
Lean-readable module at
`formal/finite-etale-keller/FiniteEtaleKeller/GeneratedPaperExample.lean`.

The pinned ramified-quintic local-field specification compiles to a portable
minimal-gauge JSON proof object and a Lean specialization of the formal
algebra-to-Keller theorem.  Two further JSON objects certify the same
ramified quintic in the power-shift family at `m=2` and the connected cubic
`T^3-T-1` in the cubic family at `n=7`.  Replay all three maps with two
independent arithmetic implementations and check all three generated Lean
specializations with:

```bash
python3 scripts/verify_arithmetic_keller_certificate.py
gp -q -f scripts/verify_arithmetic_keller_certificate.gp
python3 scripts/verify_arithmetic_keller_certificate.py \
  artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json
ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json \
  gp -q -f scripts/verify_arithmetic_keller_certificate.gp
python3 scripts/verify_arithmetic_keller_certificate.py \
  artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json
ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json \
  gp -q -f scripts/verify_arithmetic_keller_certificate.gp
cd formal/finite-etale-keller
lake env lean FiniteEtaleKeller/GeneratedArithmeticQuintic.lean
lake env lean FiniteEtaleKeller/GeneratedArithmeticQuinticStableM2.lean
lake env lean FiniteEtaleKeller/GeneratedArithmeticCubicStableN7.lean
```

From the repository root, the combined non-mutating command is:

```bash
make verify-arithmetic-compilation
```

Intentional regeneration is separate:

```bash
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py \
  --stable-parameter 2 \
  --certificate artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json \
  --lean-module FiniteEtaleKeller.GeneratedArithmeticQuinticStableM2 \
  --lean formal/finite-etale-keller/FiniteEtaleKeller/GeneratedArithmeticQuinticStableM2.lean
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py \
  --spec arithmetic/specifications/connected_cubic_stable_n7.json \
  --certificate artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json \
  --lean formal/finite-etale-keller/FiniteEtaleKeller/GeneratedArithmeticCubicStableN7.lean
```

The format and exact claim boundary are documented in
[`arithmetic/PROOF_CARRYING_COMPILATION.md`](arithmetic/PROOF_CARRYING_COMPILATION.md).

## Support-saturation compiler

The reusable module compiler computes
\((N:I^\infty)/N=H_I^0(F/N)\), module associated primes, exact regular
elements, distinguished-class annihilators and radicals, finite normal-jet
transitions, and finite-tower uniform-exponent tests.  The shared JSON schema
also records the completion ideal, parameter/base variables, normal
variables, and exact-versus-modular assurance.  Compile the checked example
with:

```bash
.venv/bin/python scripts/compile_support_saturation.py \
  schemas/support_saturation_example.json \
  --output /tmp/support_saturation_certificate.json
```

Its fast exact calibration is:

```bash
make verify-support-saturation-compiler
```

The theorem-to-programme gate ledger has a dependency-free status check:

```bash
make verify-support-saturation-paths
```

This validates the common `G0`--`G4` contract, the cubic `C0`--`C3`,
plane-JC `P0`--`P4`, and restricted-Weyl `W0`--`W5` routes, and the
deliberately open status of all three application outcomes.  The broader
compiler target depends on this check, so a numerical certificate cannot
silently promote an unmet structural gate.

<!-- status-consumer: SST1 12c5cb15e8b6de26 -->

This target also runs the standard-library-only replay of the exact
characteristic-zero degree-42 \(c_6\) certificate.  It verifies
\[
c_6\notin J_6,\qquad w_0c_6,w_2c_6\in J_6
\]
on the specialized Ritt fiber, using explicit rational multipliers and a
finite-support Macaulay dual functional.  The narrower command is:

```bash
make verify-degree42-c6-macaulay
```

Intentional modular block-Wiedemann reconstruction over the pinned 31-bit
primes is:

```bash
make refresh-degree42-c6-macaulay
```

Intentional regeneration of the homogeneous cubic-symbol atlas, the six-row
cubic double-saturation stratification, the exact degree-42 finite-jet
computation over `GF(32003)`, and the normalized plane-JC cyclic `d3`
boundary layer is:

```bash
make refresh-support-saturation-cases
```

Regenerate only the cubic search stratification imported from the proved
formal-gauge cokernel atlas with:

```bash
.venv/bin/python scripts/compile_support_saturation_cases.py --case cubic-frontier
```

This writes
`artifacts/generated-results/support_saturation_cubic_annihilator_frontier.json`.
It closes further smooth-symbol quartic saturation searches, queues the six
singular squarefree cases by annihilator type, and places the
generically-étale/Keller gate before saturation for the double-line,
triple-line, and zero symbols.  It is a routing certificate, not a new
singular saturation computation.

Compile and replay the complete deterministic nongauge complements for the
six singular-squarefree rows with:

```bash
make refresh-cubic-double-saturation-stratification
make verify-cubic-double-saturation-stratification
```

The generated
`artifacts/generated-results/cubic_double_saturation_stratification.json`
stores all representative quartic tensors and complement ranks.  On each
full complement family, the replay proves exact cotangent saturation and a
parameter-independent multiplicity-six `Ext^2` obstruction.  It also builds
the weighted Rees presentations of `Omega` and of the cokernel controlling
`Ann_B(Omega)`, proves both are `t`-saturated with their literal central
initial presentations, and thereby promotes annihilator, support-module,
cotangent, and `Ext^2` base change to every geometric parameter fiber.  It
therefore closes `C2` and fails `C1` fiberwise for these quartic models.
Reduction of the intrinsic annihilator modulo the collision maximal ideal
also gives a parameter-independent six-dimensional Nakayama module, proving
that the Kähler different needs six local generators and is not Cartier on
any geometric fiber.  It does not promote either result to higher formal
orders or close the global Keller-compatible cubic programme.

<!-- status-consumer: KDSQ6 cd423f625f1f3cd2 -->

Compile intentionally and replay the nodal sextic different-persistence
certificate with:

```bash
make refresh-nodal-sextic-different-persistence
make verify-nodal-sextic-different-persistence
```

The pinned artifact proves that the homogeneous nodal quintic and sextic
compatible spaces decompose as `42=39+3` and `64=60+4`, with quotient bases
`y^2*eta,y*z*eta,z^2*eta` and
`y^3*eta,y^2*z*eta,y*z^2*eta,z^3*eta`.  On the exact family formed with the
two quartic quotient directions, strict Rees base change identifies the
intrinsic annihilator on every geometric fiber, and `J/nJ` is the scalar
extension of the six-dimensional central module.  Thus the Kähler different
remains non-Cartier through the complete sextic normal-form quotient.  The
command does not control order-seven corrections, normality, or Keller-open
compatibility.

<!-- status-consumer: NSDP6 c5f68253995b7b6a -->

Compile intentionally and replay the all-orders nodal different-persistence
certificate with:

```bash
make refresh-nodal-all-orders-different-persistence
make verify-nodal-all-orders-different-persistence
```

The checker independently replays the cyclic graded gauge quotient
`ker(C)/im(G_nod)=Q[y,z](-3)`, then works over the universal coefficient
ring `Q[u,x,y,z]`.  It verifies that the multiplication table of
`h_nod+u*eta` is affine-linear in `u`, that the intrinsic collision
Nakayama module is `Q[u]^6`, and that assigning `u` collision weight one
gives `t`-saturated Rees presentations with literal central initial modules
for both `Omega` and `coker(B -> Omega^3)`.  The formal proof then applies
the monic graph equation `u-f` to commute the annihilator with every
polynomial or power-series specialization in the collision ideal.  Thus the
nodal Kähler different remains six-generated and non-Cartier after every
compatible formal tail.  Normality, algebraization of the infinite formal
gauge, and Keller-open compatibility are outside the claim.

<!-- status-consumer: NADPALL 60218641ccdf6fac -->

Compile intentionally and replay the all-orders singular-squarefree
different-persistence certificate with:

```bash
make refresh-singular-squarefree-all-orders-different-persistence
make verify-singular-squarefree-all-orders-different-persistence
```

The checker first verifies the formally rigid smooth central row, including
its six-generator intrinsic different.  For the six singular squarefree
symbols it verifies exact minimal presentations of the graded gauge
cokernels and the generator counts `1,2,2,3,3,4`, then constructs each
universal normal-coefficient family and computes the intrinsic `J/nJ`
presentation as six copies of the
collision residue module, and verifies strict weight-one Rees packets for
`Omega` and `coker(B -> Omega^3)`.  The formal theorem applies the graph
equations successively; each has initial form monic in an unused coefficient,
so the annihilator commutes with every polynomial or formal compatible tail.
The command therefore closes all-orders non-Cartier persistence for every
squarefree cubic symbol.  It does not prove normality, algebraize an infinite
formal gauge, or construct a Keller open.

<!-- status-consumer: SSADPALL 584a6e05374612ee -->

The older degree-42 and plane cases require Singular; the plane primary
decomposition is the longer run.  Their exact scopes, especially the
remaining full characteristic-zero saturation and order-seven gaps, the
unproved generic all-order degree-42 statement, and the still-undefined
plane Case-1 conductor/residue module, are recorded in
[`extended-geometry/SUPPORT_SATURATION_COMPILER.md`](extended-geometry/SUPPORT_SATURATION_COMPILER.md).

## HC4 projective polar atlas

Regenerate and verify the low-degree projective-degree/Segre-signature
atlas with:

```bash
.venv/bin/python scripts/verify_hc4_projective_polar_atlas.py
```

Independently compute the graph-compactification and full-polar
multidegrees for the quadratic and cubic constant-Hessian calibrations
with Macaulay2:

```bash
M2 --script scripts/verify_projective_polar_calibrations.m2
```

The Python command writes
`artifacts/generated-results/hc4_projective_polar_atlas.json`.  The
formula, the graph-versus-polar distinction, the cotangent and
Meng--Yang controls, Wang's exclusion of all sixteen quadratic-gradient
affine-degree-two/three rows, the `HC4CQ1` exclusion of all 139
cubic-gradient rows, the 319 and 307 quartic-gradient numerical rows, the
rank-one/two/three leading-quintic determinant faces, the exact
rank-three cubic Schur gap, its squarefree Hessian-discriminant
obstruction and exact witness, the resulting potential-degree lower bound
five, the generic essential-rank top-gradient/Rees support sieve, its exact
intersection with the three atlas codimension columns, and the exact
nonexistence scope are documented in
[`HC4_PROJECTIVE_POLAR_GEOMETRY.md`](HC4_PROJECTIVE_POLAR_GEOMETRY.md).

Construct the universal 56-coefficient quintic top part, verify its
gradient/Hessian/Euler/Koszul and midpoint-collision identities, build the
generic essential Hessian-rank strata, and intersect their exact support
codimensions with the atlas using:

```bash
.venv/bin/python scripts/analyze_hc4_quintic_infinity_rees.py
```

Independently certify that the generic smooth rank-one/two/three top ideals
are equal-degree complete intersections of linear type, with the stated
pure-top projective degrees, using:

```bash
M2 --script scripts/verify_hc4_quintic_infinity_rees_strata.m2
```

The Python command writes
`artifacts/generated-results/hc4_quintic_infinity_rees_strata.json`.
The pure-top Segre vectors are degeneration calibrations, not completed
constant-Hessian Segre classes.  The exact restrictions on the actual
atlas are the support-codimension/Segre-vanishing filters; lower-layer
normal-cone multiplicities remain open at this stage.  The next checker
below closes the smooth rank-three vertex packet.
The recorded replay environment is the repository Python lock together with
Macaulay2 1.22 and its `Cremona` and `ReesAlgebra` packages over
\(\mathbb Q\).

Close the smooth essential rank-three, codimension-four packet by verifying
the \(\epsilon\)-flat length-\(256\) local complete intersection, the
socle bound \(\dim(Bs_3)\ge2\), and the resulting affine-degree bound
\(\delta\ge6\):

```bash
.venv/bin/python scripts/verify_hc4_rank3_vertex_colength.py
M2 --script scripts/verify_hc4_rank3_vertex_colength.m2
```

The Python command writes
`artifacts/generated-results/hc4_rank3_vertex_colength.json` and intersects
the theorem with the atlas, excluding the signatures
`(1,4,16,64,2)` and `(1,4,16,64,3)`.  The Macaulay2 command independently
checks the complete-intersection Hilbert function and exact Fermat and
deformed local calibrations.  The universal conclusion comes from the
flatness/socle proof in
[`HC4_PROJECTIVE_POLAR_GEOMETRY.md`](HC4_PROJECTIVE_POLAR_GEOMETRY.md), not
from extrapolating those representatives.

Refine the two codimension-three packets with the rank-two
constant-kernel/Schur calculation and the rank-three ordinary-singularity
incidence:

```bash
.venv/bin/python scripts/verify_hc4_codim3_gradient_strata.py
M2 --script scripts/verify_hc4_codim3_gradient_strata.m2
```

The Python command writes
`artifacts/generated-results/hc4_codim3_gradient_strata.json`.  It proves
that a nonzero rank-two kernel restriction \(h_4|_K\) synchronizes a
constant direction, closes the squarefree binary-Hessian branch through
`HC4CD5`, and forces \(\sigma_3=16\) on the nonsquarefree remainder.  It
also checks that the rank-three Schur cubic vanishes at every isolated
singular point where the top Hessian has rank two.  The Macaulay2 replay
checks the radical synchronization powers, transverse lengths \(64\to16\),
and a nodal Hessian calibration.  These are packet restrictions, not
unconditional deletions of codimension-three atlas rows.

Apply `PGS3` to the essential-rank-two singular binary-quintic packet with:

```bash
.venv/bin/python scripts/verify_hc4_binary_root_partition_segre.py
M2 --script scripts/verify_hc4_binary_root_partition_segre.m2
```

On the open stratum where a redundant active gradient has
\(X_0\)-order one with unit coefficient at every repeated root, a root of
multiplicity \(e\) contributes exactly \(e-1\) to \(\sigma_2\).  Hence a
binary quintic with \(q\) distinct roots has \(\sigma_2=5-q\); the generic
double-root packet retains only 51 and 50 atlas rows for affine degrees two
and three.  The Macaulay2 replay checks multiplicities \(2,3,4\).
The higher-\(X_0\)-torsion failure locus remains open, so this is not an
unconditional row deletion.

## All-dimensional projective-gradient Segre machinery

Verify the canonical
\((g_0,\ldots,g_n)\leftrightarrow(\sigma_1,\ldots,\sigma_n)\) transform,
the actual affine-gradient and full-polar constructors, the leading
integrability/Euler reconstruction, and regenerate the typed family registry
with:

```bash
.venv/bin/python scripts/verify_projective_gradient_segre_machinery.py
.venv/bin/python scripts/verify_projective_gradient_normal_slices.py
.venv/bin/python scripts/verify_projective_gradient_singular_slices.py
```

Independently compute the exact plane-cotangent and
quadratic-stabilization multidegrees, and replay representative
smooth-essential normal slices, with:

```bash
M2 --script scripts/verify_projective_gradient_segre_families.m2
M2 --script scripts/verify_projective_gradient_normal_slices.m2
M2 --script scripts/verify_projective_gradient_singular_slices.m2
```

The Python commands write
`artifacts/generated-results/projective_gradient_segre_registry.json` and
the smooth and singular normal-slice ledgers
`artifacts/generated-results/projective_gradient_normal_slices.json` and
`artifacts/generated-results/projective_gradient_singular_slices.json`.
Complete multidegree/Segre vectors, top-degree-only transport controls, and
explicit families with uncomputed vectors are distinct record types.  The
normal-slice artifact records the dimension-free complete-intersection
Hilbert series, filtered missing-generator bound, exact unit-penultimate
law, and the HC4 specializations `HC4PPG7` and `HC4PPG8`.  The canonical
singular ledger records the kernel-vertex/singularity join, the exact
truncated DVR-module formula, and a repeated-root binary quintic whose
lower quartics realize active lengths \(8,3,2\).  This proves that the
singular support alone does not determine its Segre multiplicity.  The
canonical scope and the resulting restrictions on cotangent, Schur, HN,
coefficient-scheme, and boundary-normalization consumers are documented in
[`PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`](PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md).

## Coefficient-space and Kuranishi calculations

Certify the exact full-box tangent ranks for `F_4,F_5,F_6`, their visible
seed ranks and source gauges, and one nonzero quadratic obstruction in each
degree with:

```bash
.venv/bin/python scripts/verify_all_degree_coefficient_tangents.py
```

Certify the complete characteristic-zero quartic quadratic Kuranishi rank and
the explicit reduced-family tangent ranks with:

```bash
.venv/bin/python scripts/verify_quartic_full_box_kuranishi.py
.venv/bin/python scripts/verify_generic_coefficient_family_tangents.py
```

Regenerate the modular first-order source filtration with:

```bash
.venv/bin/python scripts/research_filtered_source_tangent_profile.py \
  --prime 32003 \
  --json-output artifacts/generated-results/filtered_source_tangent_profiles_mod32003.json
```

Regenerate the modular quartic slices and the Singular/Macaulay2 input files
with:

```bash
.venv/bin/python scripts/research_quartic_coefficient_kuranishi.py \
  --prime 32003 --jet-order 8 \
  --json-output artifacts/generated-results/quartic_coefficient_kuranishi_mod32003.json \
  --singular-output artifacts/generated-results/quartic_coefficient_kuranishi_mod32003.sing \
  --macaulay2-output artifacts/generated-results/quartic_coefficient_kuranishi_mod32003.m2

.venv/bin/python scripts/research_quartic_generic_component.py \
  --prime 32003 --greedy-jet-order 6 \
  --json-output artifacts/generated-results/quartic_generic_component_mod32003.json \
  --singular-output artifacts/generated-results/quartic_generic_component_mod32003.sing \
  --singular-order3-output artifacts/generated-results/quartic_generic_component_order3_mod32003.sing \
  --macaulay2-output artifacts/generated-results/quartic_generic_component_mod32003.m2
```

The optional cubic-layer compilation takes several minutes.  The CAS files
are research inputs.  Their full primary decompositions have not completed
within the available memory and are not certificate artifacts.
The theorem/computation boundary and the all-degree formal-versus-algebraic
statement are documented in
[`extended-geometry/JELONEK_COEFFICIENT_COMPONENTS.md`](extended-geometry/JELONEK_COEFFICIENT_COMPONENTS.md).

## Free-discriminant and Saito-matrix experiment

Verify the four first marked-root discriminants, the three full-target
Saito bases, the fixed-\(P\) quadratic-gauge Saito basis, all regular
marked-root lifts, and the type-\(A_3\) reflection control with:

```bash
.venv/bin/python scripts/verify_free_discriminant_saito.py
```

The weighted and cancellation branch surfaces, and the full
quadratic-gauge branch and ledger divisors, have a separate exact Singular
nonfreeness certificate:

```bash
Singular -q scripts/verify_free_discriminant_saito_nonfree.sing
```

The Singular command takes about one minute.  It verifies codimension two
and minimal-resolution length three for both Jacobian ideals; it is not a
bounded search.  The formulas, the corrected Saito--incidence proposition,
and the external-candidate gates are in
[`cancellation/FREE_DISCRIMINANT_SAITO_EXPERIMENT.md`](cancellation/FREE_DISCRIMINANT_SAITO_EXPERIMENT.md).

## Bidegree-\((3,3)\) Rodrigues survivor and sparse census

Verify the full-rank all-order pure-moment survivor, its beta factorization,
the Rodrigues identity, and the arbitrary-multiplier SIC cutoff with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rodrigues_survivor.py
```

Reproduce its exact normalized null-quadratic local certificate and the
five-variable Singular slice with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_null_quadratic_s6.py \
  --orders 2,3,4,5,6,7,8,9,10,11 --skip-solver \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_null_quadratic_s6_local.json
```

The optional unsaturated full ten-variable modular solve is an experiment,
not part of the local theorem.  It can be run by omitting
`--skip-solver` and adding `--prime 43 --timeout 600`.

Reproduce the exact characteristic-zero anti-Weyl exclusion with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_anti_weyl.py \
  --prime 0 --through 14 --backend msolve \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_anti_weyl_normalized_msolve14_char0.json
```

Reproduce the exact isolated rank-two nine-moment component and its
nonlifting signs with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_finite_prefix.py
```

This exact calculation builds the primitive anti-Weyl
moments, performs the rational Krawczyk inclusion in the radius-\(10^{-10}\)
box, proves coefficient rank two, certifies tangent rank eight on the smooth
rank-two chart, and bounds the primitive moments with signs
\(\mu_{10}>0,\mu_{12}<0,\mu_{14}>0\).  It writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_finite_prefix.json`.
It also derives the five-variable anti-Weyl square-invariant quotient and,
using characteristic-zero `msolve`, proves the unit ideal for the corrected
rank-two system on this chart.

Reproduce the exact exclusions on the generic rank-two Hurwitz chart
with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz.py \
  --characteristic-zero --backend msolve --minor 01 \
  --lambda-value 0 --orders 2,3,4,5,6,7,8 \
  --timeout 60 --memory-gb 3 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_lambda0_char0.json

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz.py \
  --characteristic-zero --backend msolve --minor 01 \
  --orders 2,3,4,5,6,7,8 \
  --mu2-pivot-boundary-reduced secondary \
  --timeout 120 --memory-gb 3 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_secondary_boundary_char0.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_hurwitz_fixed_fibres.py

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_hurwitz_root_jet_slice.py

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz_linear_incidence.py \
  --characteristic-zero --orders 2,3,4,5,6,7 \
  --emit /tmp/hurwitz_p1_A0_Bopen_incidence_QQ_m7.ms \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_linear_incidence.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_hurwitz_module_descent.py

PYTHONPATH=scripts .venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz_base_fibres.py \
  --prime 29 --through 9 --workers 4 --timeout 60 \
  --output artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_base_fibres_p29.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz_base_fibres.py \
  --prime 31 --through 9 --workers 4 --timeout 60 \
  --output artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_base_fibres_p31.json
```

The first command proves the unit ideal on the full \(\lambda=0\) fibre
after localizing the quadratic discriminant and the \(01\) channel minor.
The second imposes both successive \(\mu_2\)-pivot equations and proves
that branch unit.  The third reconstructs the same exact moments once
and proves the complete localized fibres \(\lambda=-1,1,2\) unit over
\(\mathbb Q\), writing
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_fixed_fibres.json`.
The \(521\)-bit modular construction is an exact
integer recovery: the script checks the coefficient bound
\((3m)!\,52^m\) before invoking `msolve` over \(\mathbb Q\).
The fourth command imposes the coefficient slice \(b_0=0,b_1=1\).
Moments through \(\mu_6\) give a quotient of length 687; \(\mu_7\) leaves
one rational point of local length 26.  The checker verifies its exact
fixed-flag one-sided factorization, recurrence \(\nu_{m+1}=0\), degree-one
mixed value \(-2\), and degree-\(e\) cutoff \(m>e\).  It writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_root_jet_slice.json`.
The fifth command constructs the generic characteristic-zero incidence
system on the divisor \(P_1\ne0,A=0,\mathcal B\ne0\).  It first eliminates
\(b_2\), which is valid because \(A=0,M_{01}\ne0\) force \(b_0\ne0\),
then uses \(q=b_1/b_0,\ z=b_0a_2\) and retains the sparse equation
\(\mathcal Bz+\mathcal C=0\).  The ratio chart removes the invertible
\(b_0\)-contents \(1,3,4,5,6\) from \(\mu_3,\ldots,\mu_7\), reducing the
exact source to 593,624 bytes.  This remains a source-generation result,
not a claim that the full ideal is empty.

The sixth command certifies the first component layer over
\(\mathbb Q\).  It compares primitive exact polynomials coefficient by
coefficient with independently constructed degree-preserving reductions
modulo \(29\), then invokes Singular factorization.  On
\(\mathcal B\ne0\), the incidence-plus-\(\mu_4\) projection is the
irreducible degree-\(44\) hypersurface \(R_4=0\).  On
\(\mathcal B=\mathcal C=0\), the corrected quadratic-module norm contains
exactly the removable factor
\((3\lambda+3+4q)^2=(P_1/(3b_0))^2\); the residual degree-\(26\)
hypersurface is irreducible, as is the descended cubic-in-\(z\)
\(\mu_4\) equation of total degree \(27\).  On its leading-coefficient
open, the cubic--cubic pseudoremainder with descended \(\mu_5\) is an
irreducible degree-\(40\) quadratic in \(z\); the degree-\(9\) first
leading coefficient and degree-\(28\) new leading coefficient are also
certified irreducible over \(\mathbb Q\).  It then takes the
cubic--quadratic pseudo-remainder, whose norm has the exact residual
irreducible degree-\(118\) base factor \(K\).  Reducing the descended
\(\mu_6,\mu_7\) equations in the same quadratic module and pairing their
two coordinates with the linear remainder produces primitive irreducible
base factors \(J_6,J_7\) of degrees \(130,162\).  Every exact polynomial
is compared coefficient by coefficient with its independently
constructed degree-preserving reduction modulo \(29\).  The checker writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_module_descent.json`.
It reduces the generic finite prefix through \(\mu_7\) to
\(H=K=J_6=J_7=0\) in three variables; it does not prove that common zero
set empty and does not divide through the simultaneous linear-remainder
boundary; equivalence with the original component is stated only on the
birational open \(V\ne0\).

The last two commands are bounded modular routing calculations.  They
construct \(J_8,J_9\) by the same quotient-module recurrence and saturate
each rational \(q\)-fibre by \(b_0P_1L_4S_2VR_1\).  Through \(\mu_9\),
28 of 29 fibres are unit at \(p=29\), with one length-two survivor, and
27 of 31 are unit at \(p=31\), with four finite survivors.  The artifacts
record every fibre basis, timeout, source hash, polynomial profile, and
removed factor exponent.  These calculations do not cover extension-field
\(q\)-values and do not imply a characteristic-zero unit ideal.

For a numerical-algebraic complexity estimate only, run:

```bash
julia --project=. \
  scripts/research_two_pair_sic_bidegree33_rank_two_homotopy.jl
```

It reports mixed volume \(74\,144\) for the unreduced square system
\(\mu_2,\ldots,\mu_8\).  This is not an exclusion, a solution count, or
a characteristic-zero certificate.

The complete six-entry coefficient-torus census is split into four
independent exact characteristic-zero shards:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 0 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard0.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 1897 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard1.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 3794 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard2.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 5691 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard3.json

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard0.json \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard1.json \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard2.json \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard3.json \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_sparse_six_support_screen.json
```

For support size seven, replace the four starts by
`0,2800,5600,8400`, use `--limit 2800 --support-size 7`, name the
temporary files `two_pair_sic_bidegree33_sparse_support7_shard0.json`
through `shard3.json`, and combine them into
`artifacts/generated-results/two_pair_sic_bidegree33_sparse_support7_screen.json`.
The size-six census has exactly two normalized survivors, both on the
Rodrigues orbit; the size-seven census excludes all \(11{,}200\) mixed
coefficient tori.  Boundaries are covered by the separately verified
smaller-support results.

Certify that the two size-six nonunit systems each have exactly one
complex point, rather than merely one real box, with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_survivor_rur.py
```

For support size eight, use four contiguous shards of length 3,195 with
starts `0,3195,6390,9585`, `--support-size 8 --through 12`, and combine
them into
`artifacts/generated-results/two_pair_sic_bidegree33_sparse_support8_screen.json`.
The sole timeout in that census is reproduced exactly with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 8 --through 14 --start 8384 --limit 1 \
  --timeout 600 --threads 4 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_sparse_support8_parity_msolve14_char0.json
```

The recorded run takes about eight minutes.  Validate that unit
certificate, rerun full complex RURs for the fourteen nonunit systems,
and check their explicit one-sided normal forms with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_support8.py
```

Add `--rerun-parity` only when intentionally refreshing the pinned
long-running parity artifact.  The resulting theorem is a complete
coefficient-torus classification through support size eight: any actual
bidegree-\((3,3)\) SIC counterexample has at least nine nonzero standard
monomial coefficients.

The first complete nine-entry class, the sixteen \(3\times3\) coefficient
rectangles, is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_rectangle9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_two_row_fringe9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_cross_two9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_three_line9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_three_line4329.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_full_line43119.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_full_incidence32229.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_full_incidence_remaining9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_parity_channels.py
```

The checker uses six transpose/reversal orbits.  It forms exact finite
schemes over \(\mathbb Q\) through moment fourteen and proves that all 54
\(2\times2\)-minor localizations are unit ideals.  Thus every reduced
dense component has rank one; boundaries are covered by the size-eight
theorem.  This does not run or claim the full \(11{,}420\)-support
size-nine census.

The second checker closes 96 complete-two-row/column fringe supports
in 24 exact symmetry orbits.  It also verifies the full discrete census
of 11,420 mixed size-nine supports in 2,924 transpose/reversal orbits,
of which 30 are closed by the rectangle and fringe theorems.  The third
checker closes all 576 complete-row/complete-column cross-plus-two
supports in 156 exact symmetry orbits through \(\mu_{10}\); six of these
supports are outside the mixed census.  The fourth sparse checker closes
all 480 regular three-line supports in 120 exact symmetry orbits.  Its 114
unit ideals and six unique rational fixed-flag rank-two points close the
regular \(3+3+3\) class.  The fifth sparse checker closes all 1,148
\(4+3+2\) three-line supports in 287 exact symmetry orbits.  This finishes
every one of the 1,740 nine-entry supports with an empty row or column.
Overall, 2,310 mixed supports are closed in 591 orbits, with 9,110
supports in 2,333 orbits remaining.  The sixth sparse checker closes all
1,244 full-line \(4+3+1+1\) supports in 311 exact symmetry orbits through
\(\mu_{10}\).  Overall, 3,554 mixed supports are closed in 902 orbits,
with 7,866 supports in 2,022 orbits remaining.  The next checker closes
the \((3+2+2+2)^2\) class of 816 supports in 230 orbits.  It finds
228 exact unit ideals and two degree-one rational RURs; exact flag changes
certify both residual rank-three points as one-sided, including their
relative-period recurrence, initial vanishing, low-degree mixed values,
and cutoff.  The final sparse checker treats the remaining four incidence
types: 7,050 supports in 1,792 orbits.  All 1,792 dense coefficient-torus
ideals are exact characteristic-zero units through \(\mu_{10}\), completing
the 11,420-support, 2,924-orbit mixed size-nine census.  The last command
uses the reversal-centralizer orbit cover to classify the complete exact-rank-two
reversal-parity factor family.  Exact characteristic-zero msolve
calculations close both projective semistable charts through \(\mu_6\),
and a Singular minimal-prime decomposition finds exactly two components
on the invariant-zero boundary.  Fixed-flag factorization then certifies
their all-order recurrence, initial vanishing, nonzero degree-two mixed
values, and the mixed cutoff \(2m>e\).

## Degree-seven marked-root obstruction component

The exact characteristic-zero reconstruction and terminal order-seven check
for the degree-seven classical-symbol search are replayed by

```bash
.venv/bin/python \
  scripts/reconstruct_degree_seven_order_five_rational_chart.py \
  --holdout-prime 1103 \
  --output artifacts/generated-results/degree_seven_order_five_rational_chart.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_seven_order_five_survivor.py \
  --output artifacts/generated-results/degree_seven_order_five_survivor.json
```

The first command consumes fifteen checked-in modular interpolation artifacts
and proves stable rational reconstruction, exact agreement at the unused
prime, and an irreducible length-eight zero scheme over `Q`. The second
requires Singular. It recomputes the octic-field ranks, proves that the
genuine order-five scheme is a doubled affine four-space, and proves that the
104 complete order-seven equations generate the unit ideal.

## Degree-eight marked-root obstruction component

The degree-eight characteristic-zero reconstruction and terminal order-seven
check are replayed by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/reconstruct_degree_eight_order_five_rational_chart.py \
  --holdout-prime 1009 \
  --output artifacts/generated-results/degree_eight_order_five_rational_chart.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_eight_order_five_survivor.py \
  --output artifacts/generated-results/degree_eight_order_five_survivor.json
```

The first command consumes three 31-bit and three smaller checked-in modular
chart images. Three normalized generators are stable under removal of the
last build prime; exact Buchberger completion recovers the full basis, which
agrees at unused prime 1009. Singular proves that the zero scheme is one
irreducible degree-twelve point. The second command is expensive: it verifies
the exact residue-field ranks, the two-square thickening of affine
five-space, all 220 cubic interpolation nodes plus a holdout, and the
terminal order-seven unit ideal.

## Degree-nine marked-root order-five gate

The next classical row, its exact relative complex, four modular Fitting
planes, the nonlinear screens, and the pivot-chart degree probe are replayed
by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_dc2_marked_root_degree_ladder.py \
  --output artifacts/generated-results/dc2_marked_root_degree_ladder.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/derive_degree_nine_marked_root_shear.py --jobs 7 \
  --output artifacts/generated-results/degree_nine_marked_root_shear.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_nine_relative_quantization_obstruction.py \
  --output artifacts/generated-results/degree_nine_relative_quantization_obstruction.json

for prime in 23 29 31 37; do
  PYTHONPATH=scripts .venv/bin/python \
    scripts/search_degree_seven_order_five_fitting_locus.py \
    --degree 9 --prime "$prime" --jobs 8 \
    --output \
    "artifacts/generated-results/degree_nine_order_five_scan_gf${prime}.json"
done

PYTHONPATH=scripts .venv/bin/python \
  scripts/screen_degree_seven_order_five_survivors.py \
  artifacts/generated-results/degree_nine_order_five_scan_gf23.json \
  artifacts/generated-results/degree_nine_order_five_scan_gf31.json \
  --output artifacts/generated-results/degree_nine_order_five_nonlinear_screen.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/probe_degree_nine_order_five_chart.py --prime 1009 --jobs 12 \
  --output artifacts/generated-results/degree_nine_order_five_chart_degree_probe.json
```

The exact rank row is `(227,12,142,6,149,150,371)` in the order used by the
degree-ladder table. `p=19` is excluded because its order-three rank drops to
226. The good-prime Fitting scans find two points at 23, none at 29, one at
31, and none at 37; all three recorded points are finite-free rank-six
thickenings of affine six-space, with a perfect-cube `z11` relation. The last
command uses 166 samples and holdouts to prove modular chart
denominator degree 72 and the first two numerator degree staircases. These
commands do not reconstruct a characteristic-zero component and therefore do
not authorize degree-nine order-seven PBW computation.

The full modular chart driver is prepared as

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/interpolate_degree_nine_order_five_chart.py \
  --prime 1009 --jobs 12 \
  --output artifacts/generated-results/degree_nine_order_five_chart_gf1009.json
```

This `150 x 84` interpolation batch has not yet been completed or checked in.
Its output will still be modular; multiple images, rational reconstruction,
and an unused-prime reduction are required for a characteristic-zero claim.

## Optimized Hessian-symbol search for `DC_2`

The first cross-branch search from nilpotent Hessian pencils to rank-two
classical symplectic symbols is replayed by

```bash
.venv/bin/python scripts/search_dc2_hessian_symbol_candidates.py \
  --shortlist 8 \
  --output \
  artifacts/generated-results/dc2_hessian_symbol_optimization.json
```

The command exactly screens all 1,540 one- and two-monomial cubic/quartic
Hamiltonian Hessian pencils through `N^4=0`, polynomial Cayley integrability,
the symplectic identities, and determinant one.  All 184 integrable rows are
square-zero Moyal-flat shears.  It then scores 238 canonical noncommuting
two-pencil words and quantizes the best eight through the native-support
`hbar^5` equations.  Four top rows have a native order-three rank jump
`14 -> 15`, but the factorwise Weyl lift adds one correction monomial outside
that lattice and verifies the relations through order five.  These are exact
automorphism controls and a support-saturation warning, not `DC_2`
counterexamples or unrestricted quantization obstructions.

The exact regular-index-four `N^2 != 0` control and the reciprocal `R21`
stable-chart, affine-contact, and tame degree-four graph frontier are checked
by the command below.  It also proves the fiber-preserving stabilization
no-go from the nontrivial factorization `R=x*S` and verifies the surviving
stable-mixed `U_2` chart through sextic order, including the explicit nonzero
degree-seven remainder. Finally, it verifies on the degree-six and
degree-seven rows the Euler-homotopy recurrence that eliminates every
homogeneous defect to all orders in the completed local ring. This formal
trivialization is not a finite polynomial admission certificate. The same
run computes the exact straight-line Moser Pfaffian correction: it is nonzero,
has 124 terms in source degrees two through twenty-seven, and therefore gives
a genuinely nonterminating canonical geometric series. It then constructs
the constant-Pfaffian dilation path and its exact polynomial trivializing
field, proves that its time-one map is the nonpolynomial inverse of the
generic-degree-four graph, and excludes every polynomial Darboux normalizer
whose target-`b` component has a polynomial Poisson mate on its zero fiber.
The exact invariant curve has induced field
`((y^2-1)^2/144)*partial_y`; its time form has residues `-36,36`, proving
that the kernel has no polynomial slice. This excludes every elementary
target-`b` coordinate, including all displayed corrections and elementary
`b`-moving shears. Finally, it verifies the first three rows of the infinite
tame no-slice family `F_k=a+c^k*(b+a*d)`. The general residue formula proves
that this family works for every `k>=1`; its degree-three row makes the
target-coordinate obstruction sharp, while its tied signature `(k,k)`
isolates the R21 target's split pole/weight signature `(2,3)`. The bounded
semi-invariant census through degree seven is also checked, together with its
degree-eight counterexample. The decisive all-degree calculation identifies
the polynomial constant rings as Danielewski surfaces with exponents `k` and
three. Combined with the generic time-form pole order, this excludes every
`F_k`. Finally, the checker constructs the exact `k=2` conjugacy on `I!=0`
and verifies that its Jacobian is `-I/5184`. It then crosses this
affine-modification divisor with the degree-seven Bezout coordinate
`F=A*(Q+G^2/2)-B*P`: the displayed companion `W` makes `(F,G,C,W)` a
determinant-one polynomial coordinate system, and its zero hypersurface pulls
the complete R21 `b=0` two-form back exactly. Finally, it verifies the
polynomial cancellation `s=-2/C-I^2*L/144`, the full normal--tangent form
match, and the divergence-free first transverse jet. The jet's canonical
normal vector is not locally nilpotent, so the run does not claim a finite
four-dimensional completion. At the next order, it verifies that only three
normal form coefficients remain, constructs their polynomial Bezout
correction, and reduces volume compatibility to a one-variable coefficient
recurrence. The nonzero tail beginning `p_5=-17/5388768` proves that the
canonical minimal first jet has no polynomial second completion. It then
checks that the obstruction is affine-linear in an arbitrary shift, uses the
constant-ring grading to isolate `I^2*U*K[S]`, and verifies the symbolic
nonzero top-`J` coefficient for `I^2*U*S^m`, `m>=0`. Hence every polynomial
invariant shift is excluded and this degree-seven target-`b` coordinate has
no polynomial second completion. Finally, it checks the general identity
`{F_H,W_H}|_(F_H=0)=-H_G`, proving that compatible shifts are only the
ambient symplectic shears `H=G^2/2+k(C)`, and verifies that every other
Bezout pair changes the companion by `W->W-T*F`. Thus the no-go covers the
complete reciprocal-compatible Bezout ansatz. The same run allows independent
base translations `K(G,C),L(G,C)` and proves that their only full-form
difference is their fiber-admission PDE; every admitted affine-momentum pair
therefore reduces to the already-excluded canonical form.
The same certificate performs the global boundary audit. It computes the
fixed-`P^4` coefficient-degree profile `(16,13,11,16,11,11)`, pole order
nineteen, and rank-two leading form, then verifies the elementary
standard-orbit control `omega_std+a^N*da^dc` with arbitrary pole order `N+3`.
It checks that the Pfaffian is one and that `P*dQ+e*dR` is a polynomial
primitive, realizes the Danielewski exponent crossing `2 -> 3` as the
`I`-chart of the Rees algebra of `(I,U)`, and certifies
`x|X_P(x) != 0`, `x|X_Q(x) != 0`, `A|X_e(A) != 0`, and
`X_R=partial_e`. These last identities exclude the first three canonical
Hamiltonians from the locally nilpotent locus. The quartic-root cancellation
test further proves that the Hamiltonian-LND locus in `K[P,Q,R]` is exactly
`K[R]`; the same holds in `K[e,R]`. It does not classify arbitrary
Hamiltonian locally nilpotent derivations in the full source ring. The new
generic-fiber block then computes
`disc_z(S-c)=-x^3*(6*c*(1+x*y^2)^2-5*x*y^2-6)/18`. Its branch quadratic in
`u=x*y^2` has discriminant `24*c+25`, giving four simple branch points and
genus one. Together with factorial closure of LND kernels and the analogous
geometrically integral generic `R`-fiber, this proves the complete
stable-independent classification
`{f in K[x,y,z]: X_f is LND}=K[R]`. The final associated-graded block filters
`f=sum(e^j*f_j)` by `e`-degree. Its top layer either induces a forbidden
nonzero LND of `K[x,y,z]` fixing `R`, or puts `f_m` in `K[R]`. The next layer
uses `e|delta(e)` to make `f_m` constant; a nonzero constant would give an
LND with slice `R`, contradicting the verified factorization `R=x*S`.
Therefore every Hamiltonian LND has Hamiltonian in `K[R]`. This computes the
Hamiltonian LND-generator algebra as `K[R]`, computes its kernel intersection
as `K[x,y,z]`, and proves that `Omega_21` is not polynomially
symplectomorphic to the standard form. The checker certifies the
`R21`-specific identities; the written proof invokes the standard filtered
LND lemma, Rentschler's two-variable theorem, factorial closure, and the
slice theorem. The quantization block then identifies the commuting
inverse-Jacobian frame `(delta_P,delta_Q,delta_R)`, checks the graph-normal
Weyl pair `(delta_P,P-delta_Q)`, and verifies that
`R,delta_R,delta_Q,Q-delta_P` form two commuting Weyl pairs in its
centralizer. The written PBW recurrence, using the certified common constant
ring `K[R]`, proves that these four operators exhaust the centralizer. Hence
the natural ambient-`A_3` reduction is only a transported `A_2`, not a
non-surjective endomorphism. The Hamiltonian-LND classification separately
excludes any strict filtered PBW Weyl frame for the `R21` Poisson algebra.
Only an essentially filtration-collapsing or non-PBW identification remains
open; `(DC_2)` is not settled.

```bash
.venv/bin/python scripts/verify_dc2_higher_nilpotence_r21_frontier.py \
  --output \
  artifacts/generated-results/dc2_higher_nilpotence_r21_frontier.json
```

<!-- status-consumer: C1FBC1 0f14ef01fff25097 -->

## HC4 final regular `[4]` affine-flag packet

The exact first-order moving-frame audit and its second-order flatness
prolongation are replayed by

```bash
.venv/bin/python scripts/verify_hc4_affine_plane_bridge.py
.venv/bin/python scripts/verify_hc4_affine_plane_prolongation.py
# or: make verify-hc4-relative-nilpotent-final-packet
```

The first command certifies the affine middle foliation and its two transverse
Grassmann matrices.  The second eliminates all 68 derivative jets, excludes
rank-two source-kernel motion after imposing the constant HC4 motion
determinant, and certifies the lower-rank split `p*r=0` together with the exact
`E2` and `E3` flag tensors.  It regenerates
`artifacts/generated-results/hc4_affine_plane_prolongation.json`.  The final
`p=0` closure is the degree-one hyperplane-incidence proof in
`HC4_AFFINE_PLANE_SCHUBERT_BRIDGE.md`; it is not represented as a bounded
symbolic computation.

## HC4 nonreduced Hessian--Schur module

Replay the Fermat, radial, defect-degree, and quartic-denominator
calibrations for `HC4NHM1` with

```bash
.venv/bin/python scripts/verify_hc4_nonreduced_hessian_schur_module.py
# or: make verify-hc4-nonreduced-hessian-schur
```

The written normalization argument defines the minimal denominator
`P` clearing `C^(-1)*d` and the normalized corank-two defect divisor `B_Q`.
It proves
`deg(B_Q) >= 2*deg(Q)*(r+1-deg(P))`.  For the ternary-quintic packet
`r=3`, `P^2 | det(C)` and `deg(det(C))=9`; hence any survivor with one
defect-free essential component has `deg(P)=4` and
`det(C)=P^2*ell`.  The checker confirms the exact polynomial identities and
the Fermat/radial calibrations.  The duality, DVR, and line-bundle steps are
written proofs in
[`HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md`](HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md),
not bounded CAS claims.

Continue with the exact septuple-line pole/defect ladder and its extremal
kernel gate using

~~~bash
.venv/bin/python scripts/verify_hc4_direct_septuple_linear_hessian_gate.py
# or: make verify-hc4-direct-septuple-linear
~~~

For a squarefree residual cofactor, the written proof leaves exactly six
`(pole order, primitive kernel degree, defect length)` rows. On the extremal
`(3,2,2)` row it proves that a quadratic-pencil kernel is incompatible with
a repeated Hessian line: a transverse pencil drops the boundary Hessian to
rank at most one, while a tangent pencil makes the determinant line simple.
Thus only the nondegenerate-conic kernel remains there. The checker verifies
the pole/defect arithmetic, the general tangent-pencil determinant identity,
the conic adjugate calibration, the first-normal divisibility obstruction,
and the repeated-root binary-quintic sieve.  The latter leaves four automatic
root partitions and one exceptional `2+1+1+1` orbit with anharmonic invariant
`25/4`.  An exact saturated Singular elimination then proves that the
squarefree stratum is one projective orbit, represented by `u^5+v^5`.
The complete boundary-kernel solutions then have a common quadratic factor
in every case except the `4+1` family; its primitive locus is killed by an
exact first-normal remainder.  Hence the extremal `(3,2,2)` row is empty,
and the same moving-line-image argument closes both defect-four rows.  The
checker finally verifies the constant-kernel coefficient ladders.  A
transverse kernel kills every normal jet and forces determinant zero.  A
tangent kernel either makes the septuple coefficient vanish or makes its
residual quadratic a square.  Thus the exact septuple-line packet with
squarefree coprime cofactor is empty. See
[`HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md`](HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md).

Continue through the octuple and nonuple line packets with

~~~bash
.venv/bin/python scripts/verify_hc4_octuple_nonuple_linear_hessian_gate.py
# or: make verify-hc4-octuple-nonuple-linear
~~~

For either multiplicity, the pole/defect arithmetic gives ten rows.  The
degree-one and degree-two kernel rows are already excluded by the preceding
moving-line and conic gates.  On the sole new `(4,3,0)` row, the binary
quintic Hessian must be a perfect square.  The complete classification leaves
only root types `4+1` and `3+2`; the checker verifies that their boundary
kernel equations force a common factor.  Once these moving rows are removed,
the constant-kernel coefficient ladder makes the whole quintic independent
of one variable, so its Hessian determinant is zero.  Thus `HC4NHM3` closes
`det(C)=x^8*ell` and `det(C)=x^9`, equivalently the clean quartic-denominator
partition `P=x^4`.  The next clean partition is `P=x^3*y`, whose residual
line gives incidence types `x^7*y^2`, `x^6*y^3`, and `x^6*y^2*z`;
generic lower-Smith boundaries remain separate. See
[`HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md`](HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md).

Classify the first two-line incidence packet with

~~~bash
.venv/bin/python scripts/verify_hc4_two_line_quartic_denominator_packet.py
# or: make verify-hc4-two-line-quartic-denominator
~~~

`HC4NHM4` proves that every clean `P=x^3*y`, `det(C)=x^7*y^2` packet is,
after linear normalization,
`h5=A*x*y^4+x^4*(B*y+Gamma*z)/24+D*x^5/120`, with `A*Gamma!=0`.  The checker
verifies the determinant, the complete Schur space
`s3=a*x*y^2+b*x^3`, the primitive cleared vector, and the quotient
`d^T*adj(C)*d/det(C)=a^2*x/(3*A)`.  The minimal denominator is exactly
`x^3*y` when `a!=0`.  It also verifies the forced first prolongation term
`a^2*x*t^2/(6*A)` in `h3`. See
[`HC4_TWO_LINE_QUARTIC_DENOMINATOR_PACKET.md`](HC4_TWO_LINE_QUARTIC_DENOMINATOR_PACKET.md).

Exclude its full first prolongation with

~~~bash
.venv/bin/python scripts/verify_hc4_two_line_quartic_denominator_prolongation.py
# or: make verify-hc4-two-line-quartic-prolongation
~~~

The checker retains the complete ternary quartic `r4`, quadratic `Q2`,
cubic `r3`, and arbitrary four-variable quadratic `h2`.  The degree-ten
Schur face cancels, but the next determinant face has immutable coefficient
`[x^8*t][lambda^9]det(Hess(psi))=-Gamma^2*a^3/(54*A)`.  It is nonzero for
the genuinely two-component channel `A*Gamma*a!=0`, proving `HC4NHM5` and
closing this incidence before collision equations. See
[`HC4_TWO_LINE_QUARTIC_DENOMINATOR_PROLONGATION.md`](HC4_TWO_LINE_QUARTIC_DENOMINATOR_PROLONGATION.md).

Close the other two residual-line incidences in the `3+1` partition with

~~~bash
.venv/bin/python scripts/verify_hc4_remaining_three_one_quartic_denominator_gate.py
# or: make verify-hc4-remaining-three-one-quartic
~~~

After the moving-kernel exclusions, the essential multiplicity-six line has
a constant kernel.  The transverse case has zero determinant.  In the
tangent case, the exact residual cubic is forced to be a cube, immediately
excluding `x^6*y^2*z`.  For `x^6*y^3`, matching all higher coefficients
forces the induced boundary rank-two determinant to vanish.  The other
tangent branch has multiplicity at least seven.  Thus `HC4NHM6` closes both
packets before the Schur gradient or collision equations, completing the
clean `3+1` partition. See
[`HC4_REMAINING_THREE_ONE_QUARTIC_DENOMINATOR_GATE.md`](HC4_REMAINING_THREE_ONE_QUARTIC_DENOMINATOR_GATE.md).

Close the `2+2` and `2+1+1` quartic-denominator partitions with

~~~bash
.venv/bin/python scripts/verify_hc4_two_two_quartic_denominator_gate.py
.venv/bin/python scripts/verify_hc4_two_one_one_quartic_denominator_gate.py
# or: make verify-hc4-two-one-one-quartic
~~~

`HC4NHM7` excludes both residual-line incidences for `P=x^2*y^2`.  The
coincident incidence has immutable coefficient `j1^3*x^7*y*z/18`; after
normalizing the distinct residual line, the other incidence has immutable
coefficient `-9*h1^4*x^5*y^2*z^2/(2*beta)`.  `HC4NHM8` compares the local
residual root partitions on the double line and uses the same coefficients
to exclude every power and `4+1` concurrency degeneration.  See
[`HC4_TWO_TWO_QUARTIC_DENOMINATOR_GATE.md`](HC4_TWO_TWO_QUARTIC_DENOMINATOR_GATE.md)
and
[`HC4_TWO_ONE_ONE_QUARTIC_DENOMINATOR_GATE.md`](HC4_TWO_ONE_ONE_QUARTIC_DENOMINATOR_GATE.md).

Replay the exact frontend for the last squarefree clean partition with

~~~bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_denominator_frontend.py
# or: make verify-hc4-squarefree-quartic-frontend
~~~

`HC4NHM9` proves that pole order one forces a constant kernel on each of the
four denominator lines.  A tangent direction is equivalent to
`D_v(h5) in (L^2)` and automatically supplies the double determinant factor;
a transverse direction sharpens to `D_v(h5) in (L^3)` on the
generic-corank-one locus.  Thus the open squarefree partition is reduced to
three line arrangements and sixteen tangent/transverse flag patterns.  The
finite-field commands in
[`HC4_SQUAREFREE_QUARTIC_DENOMINATOR_FRONTEND.md`](HC4_SQUAREFREE_QUARTIC_DENOMINATOR_FRONTEND.md)
are explicitly experimental and do not close this partition.

Close the all-four-concurrent arrangement and the transverse-fourth half of
the exactly-three-concurrent arrangement with

~~~bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_concurrence_closure.py
# or: make verify-hc4-squarefree-quartic-concurrence
~~~

`HC4NHM10` uses the order-eight local Hessian ladder at a four-line pencil
point to show that every nonzero rank-at-most-one determinant has a single
repeated tangent direction; the only apparent rank-zero order-eight face is
`F4*Hess_2(F4)`, whose root multiplicities cannot be the square of a
squarefree quartic.  For `(x,y,x+y,z)`, exact quartic polar syzygies and
mixed-partial compatibility exclude all eight patterns in which the fourth
flag is transverse.  See
[`HC4_SQUAREFREE_QUARTIC_CONCURRENCE_CLOSURE.md`](HC4_SQUAREFREE_QUARTIC_CONCURRENCE_CLOSURE.md).

Close all sixteen no-three-concurrent flag patterns with

~~~bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_general_position_closure.py
# or: make verify-hc4-squarefree-quartic-general-position
~~~

`HC4NHM11` computes the quartic polar-syzygy ranks for zero/one tangent
flags, the four direction charts for exactly two tangent flags, and the eight
direction charts for three tangent flags.  Every rank-drop nullspace is a
binary cone or a pure fifth power.  Together with `HC4NHM10`, this leaves
exactly eight clean squarefree rows: the triple-concurrent arrangement with
a tangent fourth flag.  See
[`HC4_SQUAREFREE_QUARTIC_GENERAL_POSITION_CLOSURE.md`](HC4_SQUAREFREE_QUARTIC_GENERAL_POSITION_CLOSURE.md).

Close those final eight rows with

~~~bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_tangent_fourth_closure.py
# or: make verify-hc4-squarefree-quartic-tangent-fourth
~~~

`HC4NHM12` computes the four exact quartic-syzygy spaces for the symmetry
representatives `RRRT`, `TRRT`, `TTRT`, and `TTTT`.  For each, Singular
verifies that the mixed-partial ideal, saturated first by the quartic
rank-two minors and then by the direction rank-two minors, is the unit ideal.
Polar injectivity forces the two ranks to agree.  The only remaining
rank-one relations are repeated tangent pairs, and the unused flag gives an
elementary incompatible monomial polar.  This closes all forty-eight
split-squarefree flag rows and completes the clean generic-corank-one packet
whose quartic denominator splits into lines.  Clean nonlinear denominator
components remain separate.  See
[`HC4_SQUAREFREE_QUARTIC_TANGENT_FOURTH_CLOSURE.md`](HC4_SQUAREFREE_QUARTIC_TANGENT_FOURTH_CLOSURE.md).

Reduce the clean smooth-cubic component invariantly with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_cubic_orthogonal_normal_form.py
# or: make verify-hc4-smooth-cubic-orthogonal
~~~

`HC4NHM19` uses `K_Q(3)=O_Q` on the smooth cubic to make the restricted
rank-two quotient an orthogonal bundle.  Its symmetric form splits the
quotient into two isotropic line bundles whose degrees add to nine.  Global
generation leaves exactly `(0,9)`, `(2,7)`, `(3,6)`, and `(4,5)`; Hessian
integrability excludes `(0,9)`.  The checker verifies the universal
hyperbolic matrix, cross-product kernel, rank-one adjugate, trace-free
splitting identity, and degree ledger.  The elliptic bundle argument and
the exclusion of the trivial summand are written in
[`HC4_SMOOTH_CUBIC_ORTHOGONAL_NORMAL_FORM.md`](HC4_SMOOTH_CUBIC_ORTHOGONAL_NORMAL_FORM.md).
The other three degree packets remain open.

Start the next nonlinear clean denominator partition with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_conic_divisible_top_gate.py
# or: make verify-hc4-smooth-conic-divisible-top
~~~

`HC4NHM13` normalizes a smooth conic to `q=x*z-y^2` and the residual line
to its tangent or secant orbit. For a general conic-divisible quintic
`h5=q*G3`, Singular proves
`I(det(Hess(h5))-k*q^4*ell):(k)^infinity=(1)` in both orbits. Thus the
complete divisible-top subrow of the double-conic denominator `P=q^2` is
empty. The checker also verifies that `h5=q^2*L` has exact conic Hessian
multiplicity three. See
[`HC4_SMOOTH_CONIC_DIVISIBLE_TOP_GATE.md`](HC4_SMOOTH_CONIC_DIVISIBLE_TOP_GATE.md).

Replay the complementary clean smooth-quartic reciprocal frontend with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_reciprocal_frontend.py
# or: make verify-hc4-smooth-quartic-reciprocal
~~~

`HC4NHM14` verifies the reciprocal identities
`det(A)=Q*ell*mu` and `adj(A)=mu*C+lambda*d*d^T`, the ten residual-line
gradient representatives, their common-factor degrees, and the four
basepoint-free quadratic boundary-matrix families. The written binary-cubic
classification splits the clean irreducible-quartic packet into one
scalar-degenerate, nine simple residual-line, and ten doubled residual-line
rows. This is a finite frontend, not an exclusion. See
[`HC4_SMOOTH_QUARTIC_RECIPROCAL_FRONTEND.md`](HC4_SMOOTH_QUARTIC_RECIPROCAL_FRONTEND.md).

Replay the first generic basepoint-free simple-line gate with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_squarefree_line_generic.py
# or: make verify-hc4-smooth-quartic-squarefree-line-generic
~~~

`HC4NHM16` constructs the 81 exact reciprocal/Hessian coefficient equations
for the squarefree-line type `d0=(x^2,y^2,0)`.  Over the total parameter
function field, a staged Singular basis is the maximal ideal in the 18 active
deformation variables; the three surviving bottom-right coefficients leave
`det(A)=0`.  This excludes the generic incidence point only.  Exceptional
specializations, the complementary line chart, and the other three
basepoint-free types remain.  See
[`HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_GENERIC_GATE.md`](HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_GENERIC_GATE.md).

Replay the first exceptional-divisor slices with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_squarefree_line_exceptional_slices.py
# or one slice: add --group NAME
# or: make verify-hc4-smooth-quartic-squarefree-line-exceptional-slices
~~~

`HC4NHM17` verifies nine exact strata inside the first visible generic-basis
divisor.  They comprise the central `H2=0` locus and its first algebraic
pivot, both generic charts and the first algebraic `m^3=48` slice at
residual-line slope `tau=0`, and both
components at `tau=-1` together with their first secondary pivots.  Every
staged basis has only `det(A)=0` support.  The arbitrary-`tau` divisor,
further hidden denominator factors, and the complementary line chart remain.
See
[`HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_EXCEPTIONAL_SLICES.md`](HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_EXCEPTIONAL_SLICES.md).

Identify the visible pivot invariantly with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_pivot_polar_geometry.py
# or: make verify-hc4-smooth-quartic-pivot-polar
~~~

`HC4NHM20` verifies that the displayed pivot is the first variation of
`Res(s^3+t^3,H)` in the explicit quadratic direction `K_tau`.  In root-value
coordinates it is `k1*h2*h3+k2*h1*h3+k3*h1*h2`, hence a smooth conic when
`K_tau` is coprime to the cubic.  The degeneration resultant is a squarefree
degree-fifteen product of factors of degrees `1,2,4,8`; every exceptional
fiber is one resultant line plus one residual polar line.  This identifies
the divisor but does not eliminate the full reciprocal-Hessian ideal on it.
See
[`HC4_SMOOTH_QUARTIC_PIVOT_POLAR_GEOMETRY.md`](HC4_SMOOTH_QUARTIC_PIVOT_POLAR_GEOMETRY.md).

Exclude the generic point of that polar conic with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_polar_conic_generic_gate.py
# or: make verify-hc4-smooth-quartic-polar-conic-generic
~~~

`HC4NHM22` uses the universal polar point `[p:q:r]=[1:3:1]` to give a
rational parametrization over `Q(tau,m,c)`.  It specializes the same 81
reciprocal-Hessian equations, verifies the ten linear pivots, and computes a
standard basis from eight selected nonlinear coefficients.  The sixth power
of each of the 18 active deformation coordinates reduces to zero.  Hence the
generic set-theoretic support is the determinant-zero boundary matrix.  The
lower parametrization and basis-denominator strata and the non-generic
two-line fibers remain.  See
[`HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md`](HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md).

Reduce the fifteen exceptional line fibers by Fermat symmetry with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_fermat_symmetry_orbits.py
# or: make verify-hc4-smooth-quartic-fermat-symmetry
~~~

`HC4NHM23` verifies covariance of the complete normalized reciprocal packet
under `tau -> lambda*tau`, `lambda^2+lambda+1=0`, and the normalized
reflection `tau -> 1/tau`.  The degree-fifteen degeneration polynomial has
three geometric orbits of sizes `3,6,6`.  The first is `tau^3=-1`, so the
two generic line-component and first-secondary certificates at `tau=-1`
transport exactly to both roots of `tau^2-tau+1`.  Only the two six-point
line-fiber normal forms remain genuinely new.  See
[`HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md`](HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md).

Close the generic points of those final two normal forms with

~~~bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_final_line_orbits.py
# or: make verify-hc4-smooth-quartic-final-line-orbits
~~~

`HC4NHM24` fixes one Fermat root and works over the exact quartic field
`Q[tau]/(tau^4-4*tau^3+10*tau^2-4*tau+1)`, whose four embeddings meet both
six-point orbit types.  Modulo this quartic, the polar conic is the union of
the resultant and residual-polar lines.  On a generic parameterization of
each component, the checker rebuilds all 81 reciprocal-Hessian equations,
solves ten linear pivots, and reduces rows
`28,34,46,49,52,56,62,64,67,71,76`.  Their eleven-monomial coefficient
matrix has exact rank eleven on both components, so all active variables lie
in the ideal and the surviving boundary matrix has determinant zero.  After
denominator clearing, the witness determinant factors with profiles
`c^10*linear^17*cubic` on the resultant line and
`c^10*quadratic^17*septic` on the residual-polar line.  This is a
characteristic-zero function-field certificate, not a finite-field
inference.  Those finite lower linear-pivot and witness-determinant strata
remain open.
The checker hash is
`sha256:ce8b2830aa93bf7c96fd1ab3164890d0321e824c1941f6ca278daf2076381684`;
the imported 81-equation builder hash is
`sha256:99f3f94c9ee4bac0a489f25916ff290b076d33e7165e88b0a952754548c419ec`.
See
[`HC4_SMOOTH_QUARTIC_FINAL_LINE_ORBITS.md`](HC4_SMOOTH_QUARTIC_FINAL_LINE_ORBITS.md).

Continue the nonzero-restriction double-conic row with

~~~bash
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-18
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-14
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-10
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-6
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group support-one-two
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group support-three
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group harmonic-four
# long generic-cross-ratio replays
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-7111
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-6211
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-5311
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-5221
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-4411
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-4321
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-3331
# long one-chart polynomial replay and function-field replay
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-4222-root-chart
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group balanced-function-fields
# or: make verify-hc4-double-conic-normal-layers
.venv/bin/python scripts/verify_hc4_double_conic_balanced_four_root_closure.py
# or: make verify-hc4-double-conic-balanced-four-root
~~~

`HC4NHM15` uses the harmonic splitting
`h5=H5(f10)+q*H3(g6)+q^2*H1(k2)` and verifies the four binary-covariant
normal layers `Phi18`, `Phi14`, `Phi10`, and `Phi6`. Exact rational
Rabinowitsch charts exclude all fourteen decic partitions supported on at
most three points and all nine four-point partitions at harmonic
cross-ratio. Further characteristic-zero unit certificates close the seven
complete arbitrary-cross-ratio partitions `(7,1,1,1)`, `(6,2,1,1)`,
`(5,3,1,1)`, `(5,2,2,1)`, `(4,4,1,1)`, `(4,3,2,1)`, and `(3,3,3,1)`.
The exact double-root value chart and root permutation close `(4,2,2,2)`.
All three `(3,3,2,2)` charts are unit over the function field in the
cross-ratio. `HC4NHM18` then closes every exceptional fiber exactly. Its
endpoint coefficients force `A=C=0`; in the remaining coordinates the
middle coefficient is `B=16*u^3`, while two normal layers give
`3*u-v^2=0` and `v*(2*u-v^2)=0`. These contradict `B!=0`. Hence every
four-point partition is empty. The support-at-least-five families remain
separate; no Schur solution is claimed. See
[`HC4_DOUBLE_CONIC_NORMAL_LAYERS.md`](HC4_DOUBLE_CONIC_NORMAL_LAYERS.md) and
[`HC4_DOUBLE_CONIC_BALANCED_FOUR_ROOT_CLOSURE.md`](HC4_DOUBLE_CONIC_BALANCED_FOUR_ROOT_CLOSURE.md).

Audit the invariant-ring target and its required clean saturation with

~~~bash
.venv/bin/python scripts/verify_hc4_double_conic_invariant_saturation_gate.py
# or: make verify-hc4-double-conic-invariant-saturation
~~~

`HC4NHM21` defines the clean elimination ideal as the four-layer ideal
saturated by the three coefficients of `Phi2`.  The checker harmonically
decomposes `h5=x^5+z^5`, verifies all four normal layers vanish, and checks
that its restriction `s^10+t^10` has nonzero discriminant.  Thus no power of
the discriminant can belong to the unsaturated ideal.  It also records a
GIT-stable repeated-root decic with zero discriminant, showing that even a
correct post-saturation discriminant certificate closes only the squarefree
open.  The all-stable target is nullcone containment for the contraction to
the binary-decic invariant ring.  See
[`HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md`](HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md).

## HC4 direct repeated-linear Hessian-factor gates

Replay the direct homogeneous filtration and its all-degree squarefree
top-Hessian obstruction with

```bash
.venv/bin/python scripts/verify_hc4_direct_homogeneous_filtration.py
# or: make verify-hc4-direct-filtration
```

This regenerates
`artifacts/generated-results/hc4_direct_homogeneous_filtration.json` and
checks the determinant-face identities used by `HC4DIR2`.  Then replay the
exact all-degree normal-form and boundary-jet identities for `HC4-DIR3--27`
with

```bash
.venv/bin/python scripts/verify_hc4_direct_double_linear_hessian_gate.py
# or: make verify-hc4-direct-repeated-linear
```

The command regenerates
`artifacts/generated-results/hc4_direct_double_linear_hessian_gate.json`.  It
checks the degree gate forcing first transverse order `j=1` and the two
normal-form Hessian determinants with linear-factor multiplicities `D-2` and
at least `2(D-2)`.  This closes exact multiplicity two from degree five and
reduces exact multiplicity three to one degree-five additive split top.  The
`HC4-DIR3b` weighted channel fixes its `w^2*ell` coefficient and then gives
the nonzero terminal multiplier `-(81/400)*(alpha^4/C^2)`, closing exact
multiplicity three for every `D>=5`.  It also checks the lower-rank exact quadruple-linear
boundary coefficients and incompatible rank-one weights, closing that
boundary for every `D>=5`.  Finally it checks the quadratic-suspension
coefficient that eliminates the apparent rank-two order-two degree-six split.
It also verifies the rank-one fraction-field integration used by `HC4-DIR5`,
which collapses the remaining order-one system to
`D=6`, `f=C*ell^6+h_6(y,z)`, and `L=ell*partial_ell`.
The final weighted-binary channel check verifies `HC4-DIR6`: after the forced
`w^2*ell^2` and `w^3` coefficients are substituted, the next determinant
channel is `(128/3375)*(alpha^5/C^3)*det Hess_(y,z)(h_6)`, closing exact
linear multiplicity four.
For `HC4-DIR7`, it checks the two possible first-motion orders on the generic
corank-one exact quintuple boundary.  The order-two packet has terminal
multiplier `(32/21)*(alpha^3/C)`; the order-one packet forces
`B_1=25*alpha^2/(84*C)`, `gamma=125*alpha^3/(5292*C^2)`, and `eta=0`, then
ends in the nonzero channel `-9*gamma^2*det Hess_(y,z)(h_7)`.
For `HC4-DIR8`, it verifies the exact order-five coefficients of both
lower-rank boundary jets.  The field equations force respectively `m=2` and
`m=3`; the first is below `D=5`, while in the second the required coefficient
has degree one and hence zero second transverse derivative.  Thus exact
linear multiplicity five is closed in every boundary rank.
For `HC4-DIR9`, it checks the two high-order exact sextuple channels.  At
`j=3`, the forced constant `w^2` Hessian entry leaves terminal multiplier
`(24/7)*(alpha^3/C)`.  At `j=2`, the rank-one field integration reduces to
the degree-eight split packet and leaves terminal multiplier
`-(1875/3136)*(alpha^4/C^2)`.  Thus only the quadratic-vector `j=1` system
survives on the generic corank-one boundary; the Jacobian of that quadratic
field is a linear boundary matrix of rank at most two.  For `HC4-DIR10`, the checker
verifies the two exact order-six coefficients in the lower-rank jets; the
degree-five tangent packet and these two jets are reductions, not exclusions.
For `HC4-DIR11`, it verifies the five forced coefficients in the
quadratic-Jacobian rank-zero ladder and the terminal nonzero multiplier
`2187*alpha^7/(4302592*C^5)`.  The surviving generic `j=1` matrix therefore
has boundary rank one or two.  For `HC4-DIR12`, it checks the outer-product
normal form of the rank-one boundary Jacobian and the exact leading ratios in
the tangent field equations.  The axial branch reduces to the earlier
linear-field system, the normal orientation is eliminated by valuation, and
the tangent orientation is reduced to `t=3,4,5,6`.  For `t<6` the boundary
quadratic is a square; for `t=6`, a nonsquare boundary quadratic requires odd
`m`.  These four tangent rows are reductions, not exclusions.
For `HC4-DIR13`, it verifies the eliminated one-form identity for
`q/ell^2` and `F/ell^m`, the primitive monomial construction, and its exact
order `t=m+3-2k`.  The written quadratic-centralizer lemma splits each row
into a binary-composite pencil and a primitive conic packet.  Below `t=6`,
the composite field is invariant along its image and the top is linear in
the complementary tangent coordinate; `t=6` also has a transverse composite
orientation.  The primitive packet has parity `(3,even)`, `(4,odd)`,
`(5,even)`, or `(6,odd)` and, below `t=6`, normal form `q=y^2+ell*z`.
This remains a reduction rather than a closure.
For `HC4-DIR14`, it checks the Hessian determinant of the invariant composite
top `f=z*G(ell,y)+h(ell,y)`.  The valuation bound eliminates `t=5,6`; at
`t=3`, the order-four equation kills the pure boundary coefficient and
boundary rank two then makes the order-five coefficient nonzero.  The sole
invariant composite survivor is `t=4`, whose exact order-six coefficient is
`-16*a0*c0^2*(m+1)*(m+2)*y^(3m-6)`.
For `HC4-DIR15`, the `z` coefficient of the field equation makes the binary
part of `v` annihilate `G`.  If it vanishes, the logarithmic identity
`G*q^(m+1)=C*ell^(3m+3)` contradicts `q mod ell!=0`; otherwise `G` is a pure
linear power, and the only order-four possibility contradicts the nonzero
boundary derivative from `HC4-DIR14`.  Thus the invariant composite
orientation is closed.
For `HC4-DIR16`, the transverse composite top splits as
`f=P(ell,y)+h(ell,z)`.  The active part of `v` would lie in the kernel of the
invertible boundary Hessian of `h`, so it vanishes; absorbing the remaining
component reduces again to the impossible logarithmic identity.  Hence all
composite rank-one packets are closed.
For `HC4-DIR17`, `ord_ell(V)=t-2>=2` puts both constant vectors in the unique
boundary-Hessian kernel for `t>=4`, again forcing the pure-quadratic
contradiction.  At primitive `t=3`, the checker verifies
`[ell]Q(F)=c*(m/2)*y^m`, incompatible with `ell^3|Q(F)`.  Thus boundary
Jacobian rank one is completely closed and only rank two remains in the
generic sextuple system.
For `HC4-DIR18`, it checks the derivative-gcd count behind the rank-two
root gate.  A primitive left-kernel generator of a linear rank-two boundary
matrix has degree at most two; consequently a nonzero binary boundary form
has one, two, or three distinct roots.  The alternative normal packet has
`f mod ell=0` and normal component `Q_ell in (ell^2)`.  This is an exact
reduction, not a closure of those four packets.
For `HC4-DIR19--20`, the one-root order-two Hessian coefficient forces its
second jet to `C*y^m`; the `Q(f)` and matrix equations then require
incompatible values of the same constant field coefficient, closing that
profile.  `HC4-DIR21` verifies the shear-basis calculations that reduce the
two- and three-root first jets to finite monomial representatives.
`HC4-DIR22` checks their Schur complements.  It records the two explicit
outer second jets, verifies the residual quadratic denominator in 27
three-root cases, and proves in the written all-degree argument that
polynomiality kills the three-root parameter.  If the first jet vanishes,
exact sextuple order forces the first positive `ell`-jet to be `ell^8` and
therefore `D>=8`.
`HC4-DIR23` checks the normal-packet square boundary field, the scalar
coefficient that kills its diagonal tangent weight, and the immutable next
matrix coefficient.  `HC4-DIR24` checks the product-rule identity behind the
delayed-jet rank collapse; the written argument applies it to every
invertible binary boundary Hessian with two missing normal jets.
`HC4-DIR25` checks the polynomial boundary kernel of both outer jets, the
order-one scalar relations, and their shared tangent factor with nonzero
coefficient `-(2a+4b)*kappa^2`.  Thus every generic-corank-one exact-sextuple
packet is closed; the lower-rank boundary is handled next.
`HC4-DIR26` checks the lower-rank order-one recurrence weights and the sole
`m=3` resonance.  The first order-one jet is empty; the second is already a
subfamily of the order-two tangent top, with exact Hessian determinant
`-32*v^2*ell^6*(c*ell^3+10*y^3)`.  Hence the lower-rank handoff is one common
degree-five top geometry with two synchronized motion orders.
`HC4-DIR27` checks the bordered determinant split for the order-two family.
Its complete potential is affine in the kernel variable with quadratic pivot
`P=2*C*ell^2+linear`; a missing tangent linear part is polynomially
impossible, while a nonzero one invokes the registered fiber reduction
`HC4RSD12`.  Only the order-one degree-five resonance remains.

Close that last exact-sextuple resonance with

```bash
.venv/bin/python scripts/verify_hc4_exact_sextuple_pure_cube_scalar_parent.py
# or: make verify-hc4-exact-sextuple-pure-cube
```

`HC4DIR28` starts from the forced scalar parent
`H+w*P+eta*w^2/2`, where `P_3` is a pure cube.  For a nonzero corner, exact
Schur completion gives a ternary constant-Hessian pencil and `HC3` excludes
collisions.  For the zero corner, the bordered unit makes `grad(P)` nowhere
zero.  The written pure-cube gradient lemma then supplies a constant unit
direction of `P`; graph coordinates factor the determinant into a binary
constant-Hessian determinant, and `HC2` excludes collisions.  The checker
replays both universal block identities and the two scalar alternatives in
the critical-point lemma.  See
[`HC4_EXACT_SEXTUPLE_PURE_CUBE_SCALAR_PARENT.md`](HC4_EXACT_SEXTUPLE_PURE_CUBE_SCALAR_PARENT.md).
For that resonance the checker also verifies the forced pure-cube pivot top
`a_3=(4v/3)*ell^3`; its complete scalar-parent form is the next uniform
classification target.
The UFD/DVR argument showing that the required half radical divides the
adjugate motion vector is the written proof in
`HC4_DIRECT_DOUBLE_LINEAR_HESSIAN_GATE.md`, not a bounded search.

## Quotient-first rank-jump fingerprints and family replays

Rebuild the R17 mechanism/high-rank quotient fingerprints and the two
Fermigier quotient fingerprints with

```bash
PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_elkies_2026_rank_jump_fingerprints.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_fermigier_rank_jump_fingerprints.py --check
```

Replay all 60,815,684 primitive parameters in the frozen Fermigier global box
and then validate the laboratory registry with

```bash
PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_fermigier_rank_jump_replay.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_nagao_section7_rank_jump_fingerprint.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_nagao_section7_rank_jump_replay.py --check

.venv/bin/python elliptic-curves/scripts/run_rank_jump_laboratory.py
```

The complete replay is retrospective and assigns no negative labels to
unlabelled fibres. See
[`elliptic-curves/notes/FERMIGIER_RANK_JUMP_REPLAY.md`](elliptic-curves/notes/FERMIGIER_RANK_JUMP_REPLAY.md).
The Nagao replay likewise assigns no negative labels; its quotient and exact
18,244,819-row ranks are documented in
[`elliptic-curves/notes/NAGAO_SECTION7_RANK_JUMP_REPLAY.md`](elliptic-curves/notes/NAGAO_SECTION7_RANK_JUMP_REPLAY.md).

Run the prospective compact-R17 shell experiment, using the unchanged
weakest-block rule on every primitive parameter with `10000 < H <= 30000`,
with

```bash
.venv/bin/python elliptic-curves/scripts/run_r17_frozen_nagao_shell.py
```

The complete population has `972697152` rows. Exact evaluation of the
preexisting 39,120-bisection atlas on the frozen, ordinary-Nagao, and random
128-row lanes, followed by frozen-rank depth replay and compact summarization,
is documented in
[`elliptic-curves/notes/R17_FROZEN_NAGAO_SHELL_2026-09-02.md`](elliptic-curves/notes/R17_FROZEN_NAGAO_SHELL_2026-09-02.md).
It certifies seven separate quotient-rank-one gains beyond the generic 17 and
runs no unrestricted point search without the mandatory residual-Selmer gate.

## First-seventeen subgroup audit for record curves 273 and 302

Replay the exact coordinate, finite-Kummer, and bad-component codes together
with the 100-digit canonical-height/theta calculation using

```bash
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_record_first17_subgroups.py --check
```

The command checks
`artifacts/generated-results/elliptic-curves/record_first17_subgroups_v1.json`.
Its exact conclusions are relative saturation index one, intersection ranks
10 and 9 with the independently selected candidate cores, free quotient ranks
13 and 14, faithful quotient Kummer dimensions 13 and 14, and surjectivity of
the first seventeen onto each complete bad-component product.  Canonical
heights remain high-precision numerical data.  See
[`elliptic-curves/notes/RECORD_CURVES_273_302_FIRST17_SUBGROUPS.md`](elliptic-curves/notes/RECORD_CURVES_273_302_FIRST17_SUBGROUPS.md).

## Six-root Mestre two-section surface and dependence relation

<!-- status-consumer: EC-M2S-GERM c7707c1e464f037b -->

Certify the characteristic-zero root-coordinate component through
`(0,25,95,143,168,205)`, its two labelled affine sections, and the rational
square root of the leading invariant with

```bash
Singular -q elliptic-curves/cas/verify_mestre_two_section_root_surface.sing
```

The verifier reconstructs the Mestre remainder recursively, localizes the
two-equation root surface away from the non-seed resultant factor, and proves
all six ordinate-eliminated residuals vanish.  It also checks that the source
is smooth with the two labelled coordinates as local parameters.

Replay the accepted visible-dependence classification

```text
P2 = P1 + V(0,+) + V(r4,-) + V(r5,+) + V(r6,-)
```

on an irreducible degree-eight fibre, including a nonzero `L(6O)` cofactor,
with the modular command below.  The attempted characteristic-zero
determinant scripts were stopped during denominator-cleared covariant normal
forms and are preserved under `archive/elliptic-curves/cas/`; they are not
certificates.  The fully lucky finite-characteristic repetitions of the
local certificate are

```bash
Singular -q -u 17 elliptic-curves/cas/verify_mestre_two_section_root_surface.sing
Singular -q -u 31 elliptic-curves/cas/verify_mestre_two_section_root_surface.sing
```

At 29 the component and section identities pass, but the selected etale
Jacobian minor vanishes because its characteristic-zero numerator is
divisible by 29.  Prime 37 is good for the original incidence tangent but
unlucky for this particular component separator: the removed resultant factor
meets the reduced seed there.  See
[`elliptic-curves/notes/MESTRE_TWO_SECTION_INCIDENCE_GERM.md`](elliptic-curves/notes/MESTRE_TWO_SECTION_INCIDENCE_GERM.md).

The older finite-field audit

```bash
python3 elliptic-curves/cas/probe_mestre_two_section_split_infinity.py
```

is retained as an independent experiment.  Its pointwise square observations
at 17, 29, 31 and 37 are superseded by the characteristic-zero rational square
identity in the root-surface verifier.

An independent degree-eight modular fibre determinant and nonzero-cofactor
check is

```bash
MESTRE_RELATION_PRIME=17 MESTRE_COMPONENT_SQRT=1 \
MESTRE_SPECIALIZE_AB=3,6 \
MESTRE_TEST_RELATION_DETERMINANT=1 \
  sage elliptic-curves/cas/verify_mestre_two_section_component_relation.sage
```

## Rational-surface quadratic rank search

<!-- status-consumer: EC-K3-RES-QBC-E6A1-RHO19 7103fa2a1a4e7ba2 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-GENUINE-Q2-MW3 cd4314040bb028f7 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-EQUATION 827d75cb8d14d7f4 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-ARITHMETIC-RANK2 387d6237125637a3 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-SPECIALIZATION-RANK7 bf1d025228805b31 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT96-A7D7-GALOIS ba008502f0e5533f -->

Generate and byte-check the low-complexity rational-surface catalogue, the
Golay/NS0031 control import, the complete polynomial degree-`(2,2)` twist
section elimination, and the new one-modulus `E6+A1` Picard-19 family with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rational_surface_quadratic_rank_search.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rational_surface_quadratic_rank_search.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_k3_dissection.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_k3_dissection.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_genuine_q2_neighbors.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_genuine_q2_neighbors.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit103_rr_weierstrass.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit103_rr_weierstrass.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit96_rr_galois.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit96_rr_galois.sage --check
```

The replay proves the displayed generic rank decomposition `1+1`, K3 height
matrix `diag(1/3,3)`, generic Picard rank 19, exact MW saturation, generic
transcendental lattice `U(3)+<4>`, four singular-K3 boundary lattices, and the
complete zero-section obstruction in the nominal norm-four and norm-six
layers.  The final replay exhausts the genuine norm-eight quadratic layer:
`119` Weyl orbits, `116` primitive classes, `90` physical degree-two classes,
and `18` nef MW-rank-3 frames in four root types.  The ansatz completeness is
limited to nondegenerate polynomial twist sections with both `x` and `y` of
degree at most two, and the genuine neighbor census is limited to isotropic
classes `2*e+2*f-w` with `w^2=8`.  The orbit-103 replay additionally proves
the complete resolved basis for `P0+P1+A3_2`, its binary quartic with rational
origin, and the explicit `2I1*+2I3+4I1` Jacobian equation.  The final replay
proves orbit-103 arithmetic rank two and exhibits its anti-invariant third
geometric direction over `QQ(sqrt(-3))`.  The orbit-96 replay attaches the
physical E6 components, guards the quadratic coefficient parent, produces the
genuine `I8+I3*+7I1` equation, and proves arithmetic rank two with the same
`chi_-3` line.  See
[`elkies-k3/RATIONAL_SURFACE_QUADRATIC_RANK_SEARCH_2026-09-02.md`](elkies-k3/RATIONAL_SURFACE_QUADRATIC_RANK_SEARCH_2026-09-02.md)
and
[`elkies-k3/E6A1_RHO19_K3_DISSECTION_2026-09-02.md`](elkies-k3/E6A1_RHO19_K3_DISSECTION_2026-09-02.md)
and
[`elkies-k3/E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md`](elkies-k3/E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md)
and
[`elkies-k3/E6A1_RHO19_ORBIT103_WEIERSTRASS_2026-09-02.md`](elkies-k3/E6A1_RHO19_ORBIT103_WEIERSTRASS_2026-09-02.md)
and
[`elkies-k3/E6A1_RHO19_ORBIT96_WEIERSTRASS_GALOIS_2026-09-02.md`](elkies-k3/E6A1_RHO19_ORBIT96_WEIERSTRASS_GALOIS_2026-09-02.md).

Search for genuinely new rational specialization directions on orbit 103,
without using its known `sqrt(-3)` direction, with

```bash
python3 elkies-k3/scripts/search_e6a1_orbit103_specializations.py

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_e6a1_orbit103_specialization_rank.sage \
  --lane small_coefficient_finalists --candidate-limit 1000 \
  --height 0 --pari-effort 0 --max-rank 24 \
  --reduction-prime-bound 500

python3 elkies-k3/scripts/test_e6a1_orbit103_specializations.py
```

The trace pass exhausts `H(k)<=3000`, promotes eight `k` values through
`H(r)<=5000`, and runs a separate `k=1`, `H(r)<=500` coefficient lane.  The
exact promotion pass certifies the distribution
`2:33, 3:251, 4:366, 5:240, 6:103, 7:7` on all 1000 selected fibres.  A
system-GP fallback completes the 617 bundled-PARI precision errors.  Every lower bound
uses exact points and a combined mod-3 finite-quotient certificate against
`Q_plus,Q_minus`.  No fibre certifies rank at least eight, and neither the
trace scores nor failed point searches are rank upper bounds.  See
[`elkies-k3/E6A1_ORBIT103_SPECIALIZATION_SEARCH_2026-09-02.md`](elkies-k3/E6A1_ORBIT103_SPECIALIZATION_SEARCH_2026-09-02.md).

### E6 rank sum three and rationalized D6 frontier

<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->

<!-- status-consumer: EC-K3-RES-QBC-E6-II-Q2-MW4 3aa5084463780acc -->

<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART a94042dd2d76797c -->

Generate and byte-check the exact `E6` rank split `2+1`, generic Picard-rank
19 K3, saturated determinant-24 Neron--Severi lattice, and same-NS rootless
MW17 impossibility with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_ii_rank3_quadratic_base_change.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_ii_rank3_quadratic_base_change.sage --check
```

Generate and byte-check the complete first genuine `q=2` shell, including
the six nef `A6+D7/MW4` frames, with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_e6_ii_rank3_q2_neighbor_candidates.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_e6_ii_rank3_q2_neighbor_candidates.sage --check
```

This is a complete exact census of the first zero-neutral degree-two shell.
It proves geometric MW rank four for six neighbours, but it does not compile
their equations or prove that all four directions descend to `QQ(r)`.

Replay the rational `D6` equation, its polynomial marked-section chart, the
exact rank-zero correspondence obstruction, and the retained height-30
regression with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_rationalized_d6_rank2_section_chart.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_rationalized_d6_rank2_section_chart.sage --check
```

The obstruction is exact only inside the declared D6 polynomial chart; it is
not a rank-four nonexistence theorem for larger rational-function charts.

Run the complete good non-`j=0` `GF(11)` shared-simple-pole two-twist census
with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_e6_shared_pole_two_twist_sections_modp.sage \
  --prime 11 --enumerate --skip-msolve
```

This tests `9,663,060` section candidates.  Every survivor is classified as
either the dependent pair `S,-2S` or a constant-section component whose
rational rank-two compatibility conic has no nondegenerate `QQ` point.  It is
a complete modular result for the declared ansatz, not a global `2+2`
obstruction.  The proof and search boundaries are recorded in
[`elkies-k3/E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md`](elkies-k3/E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md)
and
[`elkies-k3/LOWER_ROOT_TWO_TWIST_SEARCH_2026-09-02.md`](elkies-k3/LOWER_ROOT_TWO_TWIST_SEARCH_2026-09-02.md).

<!-- status-consumer: EC-K3-RES-QBC-E6-RANK4-LINEAR-CHORD 3bcfe3534656b26f -->
<!-- status-consumer: EC-K3-E6-RANK4-ROOTLESS-Q2Q4-CENSUS 2351738f44774cfe -->
<!-- status-consumer: EC-K3-E6-RANK4-DET78-GLOBAL-ROOTFUL 648ec884ce7152bb -->

Replay the systematic E6 node-collision linear-chord incidence, its exact
genus-`0/2` unordered component decomposition, the `QQ(k)` parameterization
of the genus-zero quotient, the genus-one ordered-cover descent obstruction,
the independent `2+2` height witness over that ordered field, and the
saturated geometric determinant-78 `NS` with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_linear_chord_incidence.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_linear_chord_incidence.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_rootless_low_degree_search.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_rootless_low_degree_search.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_e6_rank4_det78_niemeier_frames.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_e6_rank4_det78_niemeier_frames.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_e6_rank4_det78_niemeier_frames.sage \
  --rootless-obstruction --check
```

The unordered quotient is `P1_QQ`, with smooth point `(S,M)=(2,16)` and

```text
S=-(k^2+1)*(k^4+2*k^2+13)/(4*k*(k^2+3)),
M=-2*(k^2+1)*(k^2+3)/k^3.
```

The ordered incidence needed for the four individual sections is instead
`r^2=k^4+6*k^2+13`, birational to the rank-zero curve `52a2`.  It has no
affine `QQ`-point; degree two is minimal for a nondegenerate affine point, with
one displayed over `QQ(sqrt(53))`.  Thus geometric rank four holds over the
ordered genus-one function field, while the descended `QQ(k)` family has
exact arithmetic rank two.  The first checker emits this descent certificate
and the full geometric determinant-78 NS marking.  The second exhausts the
zero-neutral old-degree two, three, and four shells and finds no rootless frame
among their `79,837` primitive classes.  The focused Niemeier checker then
covers 1,591 primitive `A3+A2+A1` anchors in all 23 rooted Niemeier lattices;
the residual root rank is always at least 14, which exceeds the final
auxiliary vector's `13/4` norm budget under the rootless hypothesis.  The full
Niemeier run enumerates `37,397` primitive embedding-cover points and
deduplicates them to exactly `1,549` J2 frame classes.  Their root ranks run
from 10 through 17, with distribution
`10:1, 11:45, 12:249, 13:543, 14:477, 15:200, 16:33, 17:1`; hence the maximum
MW rank is 7 and no rootless MW17 frame exists.  The reciprocal-automorphism
sum equals the exact genus mass
`1463420154787/4131952105881600`.  See
[`elkies-k3/E6_RANK4_DET78_NIEMEIER_CLASSIFICATION_2026-09-03.md`](elkies-k3/E6_RANK4_DET78_NIEMEIER_CLASSIFICATION_2026-09-03.md).

### Rootless genus first-moment certificate

<!-- status-consumer: EC-K3-ROOTLESS-GENUS-MASS 2f5b874c0c22133b -->

Generate and byte-check the exact Siegel weighted root averages for the
determinant-78, determinant-948, and determinant-950 rank-17 control genera:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rootless_genus_first_moment.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rootless_genus_first_moment.sage --check
```

The three mass-normalized average signed root counts are exactly

```text
det 78:  2913380886349/59299224796,
det 948: 7957563723128755857618/562456712956783562285,
det 950: 4967763637986279936/352882035745379473.
```

The determinant-78 value is independently reconstructed from all 1,549
classes and their automorphism orders.  Every value exceeds two, so the cheap
`average<2` sufficient criterion is inconclusive on all three controls; it
does not contradict the explicit rootless determinant-948 and 950 frames.
The exact higher-ADE mass inversion is proved but not yet implemented for
these genera.  See
[`elkies-k3/ROOTLESS_GENUS_THEORY_2026-09-03.md`](elkies-k3/ROOTLESS_GENUS_THEORY_2026-09-03.md).

### Arithmetic rank transfer and marking gate

<!-- status-consumer: EC-K3-ARITHMETIC-RANK-TRANSFER 3031dd2365a29cd5 -->
<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_arithmetic_rank_transfer.sage --check
```

This replays the finite Galois-module form of Shioda--Tate, including the
representation-ring rank-transfer identity.  It certifies H3 arithmetic rank
17 and the degree-two `norm12-orbit-11952` alternate-Q80 arithmetic rank 17,
the latter directly in the rational divisor basis `(F,O,Q1,...,Q17)` before
equation compilation.  It also checks the E6 unordered-incidence fixed rank
`2/4` and the orbit-103 `2+chi_-3` split.  It applies the same schema fail-closed to the current
NS0024 completed-core path and records that its ranks `4,12,12,17` remain
geometric-only until a rational source marking and field-defined target `U`
are constructed.  See
[`elkies-k3/ARITHMETIC_RANK_TRANSFER_2026-09-03.md`](elkies-k3/ARITHMETIC_RANK_TRANSFER_2026-09-03.md).
The alternate application proof is
[`elkies-k3/R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md`](elkies-k3/R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md).

### Integral rank-transfer and character-glue calculus

<!-- status-consumer: EC-K3-RELATIVE-U-BRIDGE-LIFTING 800e22abf69b91aa -->
<!-- status-consumer: EC-K3-LOCAL-BRIDGE-MUTATION-H1C 2db88fff92ef48b9 -->
<!-- status-consumer: EC-K3-NS0024-RELATIVE-U-FIRST-EDGE-OBSTRUCTION d57544697149506f -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-CALCULUS 7eeeeaa80d9b2bf3 -->
<!-- status-consumer: EC-K3-INTEGRAL-CHARACTER-GLUE 0b76d65366279037 -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-INTEGRAL-GLUE 52de13c8443f2b7d -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-BRIDGE-PREDICTOR-BENCHMARK 3127e24cc505f646 -->
<!-- status-consumer: EC-K3-E6-DET78-PROSPECTIVE-BRIDGE-NEGATIVE d23a0abd146c2ed9 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-THETA-CONVOLUTION 5ebbd3d242fdb3db -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-CORE-GENERATION d0d78c49b44f55ac -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-REVERSE-THETA eee16ce986ec0a1f -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-WEIL-COMPRESSION 34d2abea91a265f4 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MODULAR-DIMENSION-SIEVE 9622c6eb4d8522bd -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-GENERATION 9a7a1e01cb22f62e -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-CONTROLS 3cbde45fb2cb0f17 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-DIRECTED-Q80 80de8b6727cd3409 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-BIRTH-DEATH a755a3956c4c97cb -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-ROOT-SYSTEM-SIGNATURE d32b35b66a35627c -->
<!-- status-consumer: EC-K3-NS0024-INVERSE-ADE-MUTATION 5c56f07d14129837 -->
<!-- status-consumer: EC-K3-INVERSE-ADE-PROJECTIVE-BIRTH-STRATA b4a7edb452e6dcc7 -->

Generate and byte-check the equation-free census, then replay the local
bridge and involution graph-glue theorems with

```bash
python3 elkies-k3/scripts/build_integral_rank_transfer_claim_provenance.py
python3 elkies-k3/scripts/build_integral_rank_transfer_claim_provenance.py --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_bridge_reglue.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_bridge_reglue.sage \
  --relative-u-output artifacts/generated-results/elkies-k3-relative-u-bridge-lifting-regression-v1.json
sage -python elkies-k3/scripts/certify_r17_local_bridge_mutation.sage --check
sage -python elkies-k3/scripts/benchmark_integral_rank_transfer_bridge_predictor.sage --check
sage -python elkies-k3/scripts/benchmark_e6_det78_prospective_bridge_predictor.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_theta_convolution.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_core_generation.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_reverse_theta_masks.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_weil_compression.sage --check
sage -python elkies-k3/scripts/generate_integral_rank_transfer_masked_core_neighbors.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_masked_core_controls.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_q80_defect_completion.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_q80_defect_birth_death.sage --check
sage -python elkies-k3/scripts/analyze_small_genus_defect_graphs.sage --check
sage -python elkies-k3/scripts/certify_integral_rank_transfer_root_system_signature.sage --check
sage -python elkies-k3/scripts/certify_ns0024_inverse_ade_mutation.sage --check
sage -python elkies-k3/scripts/certify_inverse_ade_projective_birth_strata.sage --check
sage -python elkies-k3/scripts/certify_integral_character_glue_calculus.sage --check
sage -python elkies-k3/scripts/certify_r17_norm12_103b2_mw_glue.sage \
  --skip-specialization-saturation --check
sage -python elkies-k3/scripts/build_integral_rank_transfer_glue_census.sage
sage -python elkies-k3/scripts/build_integral_rank_transfer_glue_census.sage --check
```

The first two commands generate and byte-check the claim-level provenance
artifact from the canonical theorem headings and the literature/novelty map.
Coverage fails closed when a labelled statement is missing or stale; this is
a documentation integrity check, not a mathematical verifier.

The relative-`U` theorem and first NS0024 application are recorded in
[`elkies-k3/RELATIVE_U_BRIDGE_LIFTING_2026-09-03.md`](elkies-k3/RELATIVE_U_BRIDGE_LIFTING_2026-09-03.md).
The local-mutation checker proves the symbolic raw-Gram/parity identities,
replays the glue-support law on all 42 edges, and certifies the new published-
R17 degree-two `4A1/MW13` fibration with maximal non-cyclic
`ZZ/4+ZZ/8` bridge.  Exact theta coefficients distinguish its frame from
both stored H3 `4A1` frames; it is not a historical-route shortcut.
Rebuild the four completed frames, compare them with the known NS0024 route,
and export the root-adapted source with

```bash
sage -python elkies-k3/scripts/search_ns0024_relative_u_bridge_lifts.sage \
  --compare-known-only \
  --export-adapted-source artifacts/generated-results/elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt \
  --output artifacts/generated-results/elkies-k3-ns0024-completed-frame-comparison-v1.json
```

The exact first-vector obstruction for the prospective `D5+E8/MW4` to
`3A1+A2/MW12` edge is replayed in degree order by the following compact
classifications.  Here `q=d(d+t)`, so the degree-two command covers
`t=0,...,18`, and the degree-three and degree-four commands cover
`t=0,...,4`.  `--summary-only` omits individual negative witnesses without
changing the enumeration or root classification.

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt \
  --root-rank 13 --degree 2 \
  --q 4 --q 6 --q 8 --q 10 --q 12 --q 14 --q 16 --q 18 --q 20 --q 22 \
  --q 24 --q 26 --q 28 --q 30 --q 32 --q 34 --q 36 --q 38 --q 40 \
  --adapt-mw-at-least 12 --rank-growth-only --include-zero-mw --summary-only \
  --output artifacts/generated-results/elkies-k3-ns0024-relative-u-degree2-fibre-summary-v1.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt \
  --root-rank 13 --degree 3 --q 9 --q 12 --q 15 --q 18 --q 21 \
  --adapt-mw-at-least 12 --rank-growth-only --include-zero-mw --summary-only \
  --output artifacts/generated-results/elkies-k3-ns0024-relative-u-degree3-fibre-summary-v1.json
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt \
  --root-rank 13 --degree 4 --q 16 --q 20 --q 24 --q 28 --q 32 \
  --adapt-mw-at-least 12 --rank-growth-only --include-zero-mw --summary-only \
  --output artifacts/generated-results/elkies-k3-ns0024-relative-u-degree4-fibre-summary-v1.json
python3 elkies-k3/scripts/certify_ns0024_relative_u_first_edge.py
python3 elkies-k3/scripts/certify_ns0024_relative_u_first_edge.py --check
```

The three complete boxes have respective maximum child MW ranks `5,7,8`, so
the MW12 target is absent.  A future hit should be rerun without
`--summary-only` and then refined to a literal ordered target `U` with
`search_ns0024_relative_u_bridge_lifts.sage`; no abstract Kneser path is
accepted as a marking.

The bridge replay covers all 42 marked H3, Q80, NS0024, and Golay-720 edges,
including exact common cores, cyclic graph glue, and complete root transfer.
The retrospective H3 predictor replay then shows that the mandatory
`K+C_new` root-budget screen rejects only 178 of 2,892 historical candidates
(`1.066x` projected speedup).  Its complete terminal fixed-core census also
finds exploratory rank-two signal: maximizing the bridge minimum retains four
rootless classes among five candidates, versus five among all 14 classes
(`2.24x` precision enrichment, `80%` recall).  Because those cores and
determinants come from successful edges, neither retrospective test passes the
prospective new-algorithm gate.
The blind determinant-78 replay closes the first untouched-shell gate as a
negative control: all 277 primitive candidates have glue-coset minimum two,
despite child root ranks 12 through 15 and 31 distinct classes in the
mass-closed truth catalogue.  Hence bridge minimum alone is not a prospective
ranking statistic; the next design must retain the decorated root profile on
glue classes.  Since this held-out genus has no rootless class, positive
rootless recall remains untested.
The theta-convolution replay then implements the exact inverse gate on the
four complete terminal fixed-core censuses.  Starting from the four cyclic
bridge determinants, it independently enumerates all fourteen reduced
positive even binary bridges, derives all 28 oriented graph multipliers from
finite-form isotropy, and checks their theta convolutions without constructing
rank-17 children.  It reproduces every independent signed-root count and
selects exactly the five rootless classes.  This is a proved fixed-core
zero-support enumerator, not a speedup or automatic core-generation theorem.
The core-generation replay inverts one further layer for maximal graph glue.
It verifies on all 84 old/new presentations that
`q_K=q_W orthogonal_sum (-q_C)` and `det(K)=det(W)*det(C)`, so a proposed
binary bridge and target frame generate the rank-15 core genus before any
isotropic graph is listed.  Exhaustive even-genus enumeration gives
`48, 8, 8, 8` candidates at the four terminal determinants and the generated
finite form selects exactly one in each case.  The theta-decorated
discriminant form is then a
complete bounded-bridge completion signature.  Its first necessary gate
rejects all 277 primitive held-out E6 cores because each already contains a
root.  This is an exact bounded core-first procedure; an unbounded bridge
cutoff, uniform genus enumeration, and a complexity or speedup theorem remain
open.
The reverse-theta replay then derives the required zero cells of the core from
each bridge graph alone.  Its 28 terminal masks contain 18--87 cells.  After
sign symmetry they compile to 14 nonredundant masks of 10--44 cells.  Lazy
exact CVP queries reproduce every zero-support decision; complete core theta
tables are computed only afterward and reproduce every signed root count.
The four all-graph runs query respectively 13, 25, 94, and 13 distinct core
cells, with no rank-17 child construction.  The accompanying coefficient
bound makes the allowed rank-15 low-order signature set finite.  Direct
realization of only those signatures is still open.
The zero-orbit Weil replay removes the full discriminant-group coefficient
space from that next modular step.  Orthogonal averaging on the `q_W` factor
preserves every mask coordinate exactly and, after theta symmetry, reduces
the four terminal coefficient dimensions from
`16,560, 44,556, 181,450, 21,804` to
`864, 5,760, 24,960, 2,880`.  Good-prime `S,T` closure certifies that the
zero-generated cyclic submodule equals the full orthogonal-orbit quotient in
all four cases.  This is an exact modular-prefilter compression, not yet a
compatible-signature enumeration or core realization.
The same replay evaluates the invariant Riemann--Roch trace formula.  The
compressed weight-`15/2` modular dimensions are
`476, 3,121, 13,488, 1,563`, while the cusp dimensions are
`472, 3,120, 13,485, 1,562`.  Rank--nullity then proves that every 10--44-cell
terminal mask leaves a cusp kernel of dimension at least 461.  Thus the
linear modular zero-mask stage is not selective; affine normalization,
integrality, nonnegativity, and lattice realization remain separate.
The mask-aware neighbour replay then supplies the first positive constructive
core-generation control.  Starting from the canonical rootful representative
of the finite-form-forced Golay-720 core genus, its seven certified
good-prime neighbours produce a nonhistorical minimum-four core whose
class-2 reverse mask is empty.  The resulting order-23 completion is rootless
and isometric to the declared rank-17 target.  The pinned artifact also stores
the deterministic 34,571-neighbour discovery telemetry.  Rerun that bounded
search, rather than only the short path certificate, with

```bash
sage -python elkies-k3/scripts/generate_integral_rank_transfer_masked_core_neighbors.sage \
  --search --check
```

This control proves neither beam completeness nor a uniform speedup.
The prospective control replay then starts from the canonical representatives
of the H3, NS0024, and Q80 forced genera and follows the extracted exact
good-prime paths.  It proves new zero-mask H3 and NS0024 cores and their
rootless rank-17 completions.  The NS0024 child is a new isometry class with
the target discriminant form.  The Q80 path is certified only as a rootless
two-cell near miss.  To rerun the corrected support-diversity search itself,
use for example

```bash
sage -python elkies-k3/scripts/search_integral_rank_transfer_masked_core_controls.sage \
  --corridor NS0024 --root-descent --mask-cap 3 --support-diversity
```

That full driver is a bounded experiment; only the extracted paths and exact
completion checks are theorem fields.
The defect-directed Q80 replay then proves the exact neighbour-survival law
and closes the two-cell near miss.  Each stored line is nonorthogonal modulo
its good prime to every physical witness of the current defect.  The replay
checks that all old witnesses leave the new dual lattice, recomputes the
replacement witnesses, and verifies the sequence `4 -> 6 -> 4 -> 4 -> 0`.
The resulting class-2 completion is a new rootless determinant-948 rank-17
construction, but not a third rootless class: exact `qfisom` tests identify it
with alternate Q80 and exclude published R17.  It has 1,313 norm-four pairs,
automorphism-group order four, and the exact target local symbols at `2`, `3`,
and `79`.  Rerun the bounded discovery records with

```bash
sage -python elkies-k3/scripts/search_integral_rank_transfer_q80_defect_neighbors.sage --check
sage -python elkies-k3/scripts/search_integral_rank_transfer_q80_defect_beam.sage --check
```

The first checks 10,000 one-step lines; the second runs the four-generation
isometry-diverse directed beam.  Their sampling bounds are experiment fields,
not completeness claims.
The birth--death replay strengthens the survival law to a complete dual-layer
transition.  For every stored line it first adjusts the isotropic lift modulo
`p^2`, expands each reverse-mask cell into the `p` affine layers
`M+r+k_0+j*y/p`, and predicts the child witnesses without a child Gram
matrix.  It separately evaluates the full `Sigma_2` transition, covering
about ten thousand dual vectors and 5,377--5,485 occupied theta cells per
edge.  Only afterward does it materialize the neighbour; independent
child-dual enumeration gives the same complete theta profile and physical
vector set on all four edges and reproduces `4 -> 6 -> 4 -> 4 -> 0`.  This
proves the exact transition and zero-defect criterion.  It does not prove that
the counted abstract
`Sigma_2` data alone suffice, or that the layered implementation is uniformly
faster than direct child construction.
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-GRAPH-REACHABILITY e02f950eba79b32a -->

The small-genus graph replay then exhausts every projective isotropic line at
three declared good primes in each of three mass-closed even ternary genera.
It constructs all singleton, pair, and triple unions; labels every directed
SCC and condensation exit; stores prime/line witnesses for shortest zero
paths; and exhaustively solves the minimum sufficient prime-set problem.  The
determinant-112 and determinant-316 unrestricted 3-neighbour graphs are
strongly connected, while their defect-directed subgraphs have respectively
a closed two-state cycle and two closed singleton traps.  Prime 5 escapes all
of those states, every tested two/three-prime union is universally reachable,
and the exact minimum sets are recorded.  All three genera have one proper
spinor genus.  In determinant 126 the exact directed distance profile is
`1,2`, and the distance-two path has signed root defects `2 -> 2 -> 0`.  Thus
fixed-prime traps and nontrivial distance occur in complete finite graphs,
while those bounded graphs alone make no all-good-prime claim.
See
[`elkies-k3/DEFECT_GRAPH_SMALL_GENUS_DYNAMICS_2026-09-03.md`](elkies-k3/DEFECT_GRAPH_SMALL_GENUS_DYNAMICS_2026-09-03.md).

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MARKED-ROOTLESS-REACHABILITY 354cc7a9fc81f33e -->

The all-good-prime completion is theorem-level but has no computational
replay.  A finite discriminant/glue marking has compact-open stabilizer and
therefore defines a level class set in the sense of Chenevier, Theorem 5.9.
Within that set, a state reaches zero support through finitely many directed
good-prime edges exactly when its marked spinor/level component contains a
zero-support state.  In fact every sufficiently large prime with the required
spinor displacement gives a direct marked neighbour in the chosen rootless
class; the dual zero-layer survival law makes this edge automatically kill
every parent physical witness.  The resulting finite prime set is
non-effective.  This supplies neither a small-prime bound nor an equation or
elliptic-neighbour lift; see Theorem H0i.3 in
[`elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md`](elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md).
The root-system signature replay then enumerates every physical completion
root `k+c` on the four NS0024 stages and records all pairwise inner products.
It recovers `D5+E8`, `3A1+A2`, `3A1+A2`, and rootless, together with root
ranks, component discriminants, primitive closures, and torsion quotients.
At the third stage the core is rootless but the twelve completed roots occupy
five nonzero order-191 graph-glue labels.  The same metric classifier recovers
the Q80 `4A1` and `A1` controls.  The marked coordinates, rather than the
pairwise metric alone, certify primitive closure and exact Mordell--Weil
torsion.
The inverse-ADE replay composes this classifier with the affine dual-layer
transition before constructing the first NS0024 child.  It compiles 140
modular root-incidence forms on the `p=17` line: six vanish and 134 are
nonzero.  Exhaustive graph-glued affine CVP then predicts exactly those six
survivors, no births, no extra roots, and metric `3A1+A2`; the subsequently
materialized child has the identical physical root set.  This proves a finite
necessary-and-sufficient target predicate for a fixed line, not completeness
or a runtime bound for finding such a line.
The projective-birth-strata replay eliminates the affine variable by scaling
a born root `v` to `z=p*v`.  It exhausts 48 state/prime cases in three
mass-closed ternary genera and checks all 346 predicted root sets against
independent child enumerations.  The isotropic-quadric complement predicts
all 192 rootless lines without a marked target core.  A six-line index-two
graph-glue control additionally has 16 nonzero-coset births per line and exact
predicted/materialized equality, giving 352 exact comparisons in total.  See
[`elkies-k3/INVERSE_ADE_PROJECTIVE_BIRTH_STRATA_2026-09-03.md`](elkies-k3/INVERSE_ADE_PROJECTIVE_BIRTH_STRATA_2026-09-03.md).
The character replay exhausts the E6 `2+1` and `2+2` involution graphs after
the declared factor-12 integral scaling.  The norm-twelve byte check reuses
the pinned full saturation record; generating that artifact without the skip
flag runs every possible saturation prime in a separate memory-bounded
process.  The accompanying mass criterion is
a terminating genus-completeness certificate at J2 level; equation lifts and
J1 surface-automorphism orbits remain separate.

<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-TWO-TWIST-POLYNOMIAL ea0496c9566cfdc3 -->

<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-LOW-SLICE-ELIMINANTS 43d297285eb3655b -->

Replay the exact two-marked D5 seed and its complete `GF(11)` polynomial
two-twist census with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_d5_two_marked_two_twist_polynomial_modp.sage \
  --prime 11

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_d5_two_marked_two_twist_polynomial_modp.sage \
  --prime 11 --check
```

The exact part certifies an `I1*` seed with two independent invariant
sections.  The modular part checks 146,410 section incidences and finds one
quadratic twist with two distinct-`x` polynomial sections.  It is a modular
lifting candidate, not a characteristic-zero `2+2` or height-determinant
certificate.

Replay the unique `11`-adic and `13`-adic lifts in their regular low-section
slices, and compute the exact characteristic-zero eliminants, with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/lift_d5_two_marked_two_twist_low_section_slices.sage \
  --run-eliminants --eliminant-timeout 900 --threads 4

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/lift_d5_two_marked_two_twist_low_section_slices.sage \
  --run-eliminants --eliminant-timeout 900 --threads 4 --check
```

The saturated `A=-1` and `A=2` slices have `t`-eliminant degrees 142 and
128, with rational factor degrees `2,2,4,4,42,88` and `2,4,44,78`.
The modular points select the irreducible degree-88 and degree-78 factors.
Neither eliminant has a linear rational factor, so neither complete slice has
a `QQ` point.  This exact result does not close the full D5 polynomial chart.

<!-- status-consumer: EC-K3-RES-A4-TWO-POINT-TATE-SLICE-OBSTRUCTION b9729a0a8f2f17be -->

Replay the exact obstruction in the normalized two-point A4 Tate slice with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_a4_two_point_tate_slice_obstruction.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_a4_two_point_tate_slice_obstruction.sage --check
```

The `I5` jet has two branches.  Both force either an extra repeated
discriminant root (generically an `I2` fibre) or a linear dependence between
the two automatic points.  This is exact only for the declared `r=s=1` Tate
slice, not for a general two-marked A4 surface.

### Published-R17 barcode-targeted genus-one bisections

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-BISECTION-PILOT 80fa6e59107cc9e6 -->

<!-- status-consumer: EC-K3-R17-RANK28-INTEGRAL-CHARACTER-GLUE 617f1838d8581fcd -->

Generate and byte-check the norm-eight trace census and the exact genus-one
bisection pencil members through all eleven rank-28 exceptional targets with

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_genus_one_bisections.sage

sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_genus_one_bisections.sage \
  --check
```

The default pilot completely equation-ranks the 63,925 minimum-norm-eight
translation classes and executes the cheapest trace `-P2-P5`.  All eleven
target-fitted members have irreducible squarefree quartic branch polynomials,
exact Kummer barcode matches, and independent height-16 anti-invariant
sections.  Fitting one member of a genus-one pencil through each already known
point is not a point-discovery or rank-32 certificate.  See
[`elkies-k3/R17_RANK28_GENUS_ONE_BISECTIONS_2026-09-02.md`](elkies-k3/R17_RANK28_GENUS_ONE_BISECTIONS_2026-09-02.md).

### Published-R17 multisection-visibility filtration

<!-- status-consumer: EC-K3-R17-MULTISECTION-VISIBILITY-FILTRATION 2f41e9f4236f6c9e -->

Generate and byte-check the genus-qualified filtration at all four rank-25--28
controls with

```bash
sage -python \
  elkies-k3/scripts/certify_r17_multisection_visibility_filtration.sage

sage -python \
  elkies-k3/scripts/certify_r17_multisection_visibility_filtration.sage \
  --check
```

The checker fits the common norm-eight genus-one bisection pencil through all
38 displayed exceptional generators.  This proves that the literal all-genus
displayed filtration is full in degree two, with dimensions `8,9,10,11`.  It
separately imports the complete rational-bisection dimensions `5,3,2,1` and
retains the incomplete labels on the sampled degree-three and degree-four
rational-curve layers.  It does not claim a full Mordell--Weil quotient or an
exhaustive rational trisection/quadrisection calculation.

Extract the common integral character/glue carrier from that exact certificate
with

```bash
python3 elkies-k3/scripts/certify_rank28_integral_character_glue.py

python3 elkies-k3/scripts/certify_rank28_integral_character_glue.py --check
```

This proves that all eleven lifts fitted from the chosen trace pencil repeat
the same index-two graph glue `<16>_+ + <16>_- -> <8>+<8>` with trace
`-P2-P5`, while their eleven deck squareclasses remain distinct.  The
signature is relative to that pencil; it does not put eleven directions on
one cover or produce a new specialization.

### Published-R17 frozen-quartic simultaneous-splitting search

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-SIMULTANEOUS-SPLITTING-H10000 40fb0bc465e3e95c -->

Run the exhaustive compact projective scan through height `10,000` and the
canonical pointed-Jacobian subgroup search through `30P` with

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_simultaneous_splitting.sage

sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_simultaneous_splitting.sage \
  --check
```

The compact scan checks all `121,589,943` primitive parameters in the declared
box.  Exact replay finds only the fitted rank-28 control `-9529/5471`; there is
no new simultaneous split in the box or in the eleven cyclic subgroup ranges.
This is a bounded negative search, not a global rational-point obstruction.

### Published-R17 mixed-trace genus-one splitting search

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-MIXED-TRACE-SPLITTING-H10000 c7aa09836b842b60 -->

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_mixed_trace_splitting.sage

sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_mixed_trace_splitting.sage \
  --check
```

This constructs `77` exact target-fitted quartics from the seven cheapest
distinct finite-pole norm-eight traces and repeats the complete height-10,000
scan, requiring survivors from at least two trace pencils.  Exact replay finds
only the original control, where all 77 split.  The negative conclusion is
restricted to the selected trace prefix and compact parameter box.

### Published-R17 high-throughput genus-one splitting search

<!-- status-consumer: EC-K3-R17-GENUS1-HIGH-THROUGHPUT-SPLITTING cad3d98ce58c89e7 -->

The unbiased 100-plus-43 production profile replays with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_elkies_2026_genus_one_bisection_splitting.sage \
  --norm8-count 100 \
  --equation-pool-size 1000 \
  --output artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json \
  --local-directory artifacts/local/elkies-k3/r17-genus-one-bisection-splitting/production-v1
```

It compiles all 143 quartics into dynamic Legendre masks, evaluates 5,474,328
compact or deterministic large-coordinate parameters plus 900
pointed-Jacobian parameters, and exact-tests all 77,704 modular extremes.  The
single simultaneous split is at `t=1/25`; its norm-twelve point has one exact
finite-quotient direction beyond generic MW17.  The compact height-2000 box is
exhaustive, while the trace and other parameter populations remain bounded.
See
[`elkies-k3/R17_GENUS_ONE_BISECTION_SPLITTING_SEARCH_2026-09-02.md`](elkies-k3/R17_GENUS_ONE_BISECTION_SPLITTING_SEARCH_2026-09-02.md).

### Norm-twelve `0x103b2` Jacobian and targeted split sweeps

<!-- status-consumer: EC-K3-R17-NORM12-103B2-MW-LATTICE-SIEVE aa0d8718eb57de6f -->

Use the newly found point at `t=1/25` to construct the cover Jacobian, extract
17 independent points, and sieve the other 142 compiled covers over the
radius-two/support-five lattice shell:

```bash
sage -python \
  elkies-k3/scripts/search_r17_norm12_103b2_mw_lattice.sage
```

The complementary unit/support-eight shell is:

```bash
sage -python \
  elkies-k3/scripts/search_r17_norm12_103b2_mw_lattice.sage \
  --max-support 8 --coefficient-radius 1 \
  --output artifacts/generated-results/\
elkies-k3-r17-norm12-103b2-mw-lattice-unit-support8-v1.json
```

The first shell has 6,991,556 vectors and the second has 9,746,882; their
intersection has 242,114 vectors, so the union contains 16,496,324 signed
coefficient vectors.  Thirty-two exact good-prime reductions leave nine and
four local survivors respectively.  Exact replay recovers only the known
`t=1/25` overlap with norm-eight `0x0f6b1`; all other finite survivors are
nonsquares and one vector in each shell is exceptional for the affine inverse
map.  This proves rank at least 17, not exact rank 17: the PARI rank upper-bound
run and eclib saturation at `40251553` were stopped without results.  The
lattice shells are bounded and do not exclude other rational overlaps or
points outside the known subgroup.

The seven covers appearing as false all-prime survivors can be tested more
directly through their genus-three quotient curves:

<!-- status-consumer: EC-K3-R17-NORM12-103B2-HARD-FIBRE-PRODUCT-H300000 b4fef7ab54b922e0 -->

```bash
sage -python \
  elkies-k3/scripts/search_r17_norm12_103b2_hard_fibre_products.sage \
  --height 300000

sage -python \
  elkies-k3/scripts/search_r17_norm12_103b2_hard_fibre_products.sage \
  --height 300000 --check
```

All seven quotients have zero affine rational points in PARI's naive-height
`300000` search, so there are no simultaneous splits in that bounded range.
This does not determine the complete rational-point sets of the genus-three
curves.

### `0x103b2` hidden frame and pointed-cover rank control

<!-- status-consumer: EC-K3-R17-NORM12-103B2-ISOTROPIC-FRAME 47f3a0eb7ee50bcb -->
<!-- status-consumer: EC-K3-R17-POINTED-COVER-JACOBIAN-CONTROL-H10000 4bb087b3a1ebc684 -->

Classify the primitive isotropic bisection class and compare its orthogonal
frame with both certified determinant-948 rootless `J2` classes:

```bash
sage -python \
  elkies-k3/scripts/classify_r17_103b2_isotropic_frame.sage
sage -python \
  elkies-k3/scripts/classify_r17_103b2_isotropic_frame.sage --check
```

The expected identification is `published-R17-J2-class`: the frame is
rootless of rank 17 and determinant 948, isometric to published R17, and not
isometric to the alternate Q80 frame.

### Minimal `J2` accessibility through all norm-twelve pencils

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-MINIMAL-ACCESSIBILITY 631f50389e0a3283 -->

Classify all 43 exact norm-twelve genus-one bisection frames against both
mass-complete determinant-948 rootless controls:

```bash
sage -python \
  elkies-k3/scripts/classify_r17_norm12_isotropic_frames.sage
sage -python \
  elkies-k3/scripts/classify_r17_norm12_isotropic_frames.sage --check
```

The exact distribution is 33 published R17 and ten alternate Q80.  Every
alternate witness is a nef degree-two fibre sharing the published zero with
zero-section degree one.  Together with the general lower bound for distinct
`J2` classes, this proves elliptic incidence distance two.  See
[`elkies-k3/J2_GEOMETRIC_ACCESSIBILITY_2026-09-03.md`](elkies-k3/J2_GEOMETRIC_ACCESSIBILITY_2026-09-03.md).

### Direct `norm12-orbit-11952` alternate-Q80 equation

<!-- status-consumer: EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION 077c6409d76cbe63 -->

Compile the cheapest degree-two witness directly on the published R17 model:

```bash
sage -python \
  elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage

sage -python \
  elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage --check
```

The exact replay solves the eight-by-ten Riemann--Roch congruence for
`D=(3,2,w)`, constructs the binary-quartic pencil and pointed Jacobian, and
checks a polynomial `(8,12,24)` K3 model with irreducible squarefree
discriminant and `24I1` fibres.  It splits off `<D,D+O_old>=U`, proves the
rootless determinant-948 frame is alternate Q80 rather than published R17,
and transports sixteen old sections plus `orbit-0adf9` to a saturated
seventeen-section basis.  See
[`elkies-k3/R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md`](elkies-k3/R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).

<!-- status-consumer: EC-K3-R17-NORM12-11952-CONTROL-J-PREIMAGES 1ef38474a0d7f629 -->

Test exact rational-fibre accessibility of the four published rank-25--28
controls under the alternate `j`-map:

```bash
sage -python \
  elkies-k3/scripts/certify_r17_norm12_11952_control_j_preimages.sage

sage -python \
  elkies-k3/scripts/certify_r17_norm12_11952_control_j_preimages.sage --check
```

The four degree-24 preimage polynomials have no rational finite root and no
root at infinity.  Hence the published controls are not rational fibres of the
alternate family, even after allowing quadratic twists.

Run the exact `H=10000` point-map-relation control on the deterministic
seeded sample of ten other pointed covers:

```bash
sage -python \
  elkies-k3/scripts/control_pointed_cover_jacobian_ranks.sage
sage -python \
  elkies-k3/scripts/control_pointed_cover_jacobian_ranks.sage --check
```

Each control has only its two signed pointing points in the bounded PARI
search, hence no nonbase image; the positive regression has 58 nonbase images
of relation rank 17.  This is an exact bounded-search contrast, not an upper
bound for any control Jacobian.
