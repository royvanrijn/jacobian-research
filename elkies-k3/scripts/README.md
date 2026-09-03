# Elkies K3 script map

This directory is the executable history of the elliptic-neighbour reconstruction.
It contains current proof replays, reusable compiler code, active searches, regression
fixtures, and some historical attacks that remain in the root because other notes or
launchers still refer to them.

The base inventory audit was made against repository commit
`4eac04a442a132b696a85384f08b81569870a940`; later proof and search entry points
are documented in the topical sections below. Every file present at that audit
in `scripts/` and `scripts/archive/` was enumerated. Current and ambiguous entry
points were inspected against their code headers, outputs, certificates, and
the notes that consume them; the archive is classified by historical programme
rather than promoted file-by-file.

<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-MOD131-HORIZONTAL 2249c509c1217d7c -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-HORIZONTAL 688f0a5f6d989e9c -->
<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R18-COVER 6b4ee5bbc1afc01e -->
<!-- status-consumer: EC-K3-ELKIES-2026-R19-PAIRED f1e135d2ba803e80 -->
<!-- status-consumer: EC-K3-BISECT-BIQUADRATIC-R19 707bffd8b85f8f3e -->
<!-- status-consumer: EC-K3-H3-Q12O5867-POINT-FACTORY 9399c93ee42ee2a4 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-TWO-PRIMARY-BOUNDARY 783482d8f700105d -->

## Local Sage 10.9 installation

On this machine the pinned SageMath 10.9 build is installed from conda-forge at

```text
/home/royvanrijn/.local/share/jacobian-sage-10.9
```

It is isolated from the repository and system packages.  The launcher in this
conda build does not accept the historical `sage -python` spelling; invoke Sage
Python scripts with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python script.sage
```

or prepend its `bin` directory to `PATH` and use that environment's `python`.

## Status vocabulary

| status | meaning |
|---|---|
| **ACTIVE_PROOF** | Replays a claim used by a current note or pinned certificate. |
| **ACTIVE_COMPILER** | Reusable exact equation/Riemann--Roch/Jacobian infrastructure. |
| **ACTIVE_SEARCH** | Current bounded exploration. Its output is not a theorem until separately certified. |
| **UNPROMOTED_RESULT** | An exact-looking output exists, but the canonical note, locks, or `MATH_STATUS.json` have not admitted it yet. |
| **REGRESSION** | Reproduces a specialization or previously solved compiler case. |
| **HISTORICAL_DIAGNOSTIC** | Preserves a failed or superseded route and the evidence that rejected it. |
| **ARCHIVED_SNAPSHOT** | Historical source copy. Never use it as the current proof entry point. |

Location alone is not authority. A root-level script is authoritative only when a current
note or pinned certificate names it. Conversely, a failed script is worth retaining when
it records a bounded negative result, a normalization bug, or a useful local model.

## Quick decision guide

- For the story and lessons, read `../ELKIES_K3_PROCESS_ATLAS.md`.
- For pinned commands and hashes, use `success-path/ledger.json`.
- The physical equation proof continues from the P1229-pointed
  `4A1/MW13` q8/orbit376 child through q12/orbit5867 to the rootless `24I1`
  endpoint. Source identity, Picard rank 19, and full saturated R17 are exact.
- `lift_h92_q24_orbit42_resolved_rr_qq.sage` and
  `certify_h92_q24_orbit42_a11_equation_marking.sage` are the exact equation
  and marking proofs for that child.
- Work next in the compact published `t` chart on residual 2-descent and
  arithmetic specialization. The exact rank-25--28 positive controls and
  same-curve fail-closed Selmer gate precede every expensive point search.
  q12/orbit4484 remains a certified fallback, not an open requirement. Do not
  restart the withdrawn q6/orbit1307, q323, changed-zero, zero-pole,
  point-transport, or halving searches unless a specialization certificate
  requires one.

## Current proof and compiler entry points

### Fixed-u marked Q80 third-q12 search

- `certify_q80_fixed_u_marked_third_q12.sage` is the fail-closed
  **ACTIVE_SEARCH** for
  `D7+D5/MW5 --q12--> A5+A3+3A1/MW6`.  It specializes a predeclared
  low-height rational `u` list on the exact coefficient curve, constructs the
  source and both q4 children exactly over `QQ`, tests the forced first Q80
  marking, and independently searches the `D7+D5` child over several good
  primes.  Modular group-law words are retained only when they have
  `P.O=2`, height `8`, and identity-component fingerprints at both star
  fibres.  It can also export the direct denominator-two section scheme
  `h=W^2+h1*W+h0`, `x=N/h^2`, `y=M/h^3` in either a sparse auxiliary-`M`
  chart or a dense square-recursive chart and can run a bounded `msolve`
  probe.  Its current output is reconnaissance: the terminal marked-q12
  status remains unreachable until projective CRT/LLL, the complete connected
  `D7+D5` quotient, exact parent/child maps, and the
  `A5+A3+3A1/MW6` fibre gate all pass.

The affine audit of the forty declared low-height values finds no split forced
source marking; the separately tested projective point `u=infinity` is a good
exact q4/q4 specialization but has the same nonsquare marking obstruction.
The exact first-marking cover is simplified and checked by
[`simplify_q80_first_marked_cover_qq.sage`](simplify_q80_first_marked_cover_qq.sage)
to
`Y^2=2*s*(s^2+s+1/3)*(s^2+2*s+2)`.  A completed `hyperellratpoints` search
through height `10^8` finds only its visible branch points, so it supplies no
useful nonbranch rational marked specialization within that bound; this is
not a global rational-point theorem.
Among primes at most `19`, the audit selects `u=-2, p=19` as the first good
affine reduction.  Replay the corrected signed shell through length ten with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_q80_fixed_u_marked_third_q12.sage \
  --u=-2 --prime 19 --word-length 10 \
  --output artifacts/generated-results/q80-fixed-u-marked-third-q12-search.json
```

This bounded search has eight unsigned polynomial section pairs and no target
profile through length ten.  An optimized extension through length fifteen
visits 21,871 distinct subgroup points and still finds no target profile.  The
exact fourth/eighth-multiple Shioda replay now supersedes that word bound: the
eight pairs span only rank three, with saturated Gram
`[[11/4,-1/2,-5/4],[-1/2,3,-1/2],[-5/4,-1/2,1]]`, determinant `2`, and
trivial torsion residual.  Its complete height-eight shell is empty.  Thus the
required horizontal uses the two missing non-polynomial MW directions.

`scan_q80_po1_msolve_slices.sage` exhausts the direct `P.O=1` charts without
random slicing.  At `u=-2, p=19`, the 342 finite-pole `(z,l)` slices contain
exactly one sign pair, at `(15,4)` and `(15,15)`.  The decoded section is
literally replayed on the reduced Q80 equation and has height `19/4`, but its
coordinates in the rank-three polynomial lattice are `(1,1,0)`.  All 18
pole-at-infinity `l != 0` slices are empty.  Therefore neither `P.O=1` chart
supplies a missing direction; the next section search starts at `P.O=2`.

Reproduce the finite `P.O=1` export, exhaustive scan, and augmented height
replay with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_q80_fixed_u_marked_third_q12.sage \
  --u=-2 --prime 19 --prime-audit-only --direct-pole-order 1 \
  --direct-pole-location finite --direct-chart auxiliary \
  --direct-msolve-dir artifacts/generated-results/q80-fixed-u-minus2-p19-po1-msolve \
  --output artifacts/generated-results/q80-fixed-u-minus2-p19-po1-export.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_q80_po1_msolve_slices.sage \
  --input artifacts/generated-results/q80-fixed-u-minus2-p19-po1-msolve/q80-third-q12-um2d1-p19-po1-finite-auxiliary-sign+1.ms \
  --output-dir /tmp/q80-po1-p19-finite-slices-certified \
  --output artifacts/generated-results/q80-fixed-u-minus2-p19-po1-finite-slices.json \
  --workers 4 --msolve-threads 1 --timeout 30

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_q80_fixed_u_marked_third_q12.sage \
  --u=-2 --prime 19 --word-length 1 --height-shell \
  --po1-slice-certificate artifacts/generated-results/q80-fixed-u-minus2-p19-po1-finite-slices.json \
  --output artifacts/generated-results/q80-fixed-u-minus2-p19-height-shell-with-po1.json
```

The direct `P.O=2` exporter and its sliced variants are now superseded.  The
complete MW5 embedding identifies `Q=H-P4` as a polynomial section.  The
exact degree-four polynomial-section scheme has radical degree 12; its
quadratic factor `T^2+12T+3` supplies the unique Frobenius orbit with
`P.O=2`, height `8`, identity components at both star fibres, and the pinned
intersection fingerprint.  Two independent decoders certify the same orbit.

The connected quotient is also compiled.  Smith saturation of
`{1,x,h*z_H}` followed by shifted Popov reduction gives a seven-dimensional
ambient.  The resolved D7 complete ideal `(Y,U^2,ZU,Z^3)` contributes four
conditions and the complete D5 fibre contributes one; the rank-five matrix
leaves the required two-dimensional pencil.  Removing `h^2*x-Nx` gives an
irreducible moving equation of degrees `(2,9,3)` in `(new base,W,x)`.
Primitivity, square zero, old-fibre degree three, separability, and K3
adjunction certify generic genus one without invoking Singular's unsupported
normalization over an algebraic finite field.

Replay and verify the immutable 18-file p=19 pin with:

```bash
sage -python elkies-k3/scripts/certify_q80_po0_rur_third_q12_modp.sage
sage -python elkies-k3/scripts/compile_q80_third_q12_um2_p19_resolved_pencil.sage
sage -python elkies-k3/scripts/verify_q80_third_q12_um2_p19_resolved_genus.sage
python3 elkies-k3/scripts/pin_q80_third_q12_p19_pipeline.py
```

`sample_q80_third_q12_weierstrass_mod19_quadratic.sage` is the exact mapped
fibre primitive for the child Jacobian.  It uses the rational simple infinity
branch `xi=-6`, a finite integral basis, and exact local regularity at the
double branch `xi=-16` to recover `L(2P)` and `L(3P)`.  The deterministic
batch now retains 156 smooth mapped fibres: 148 interpolation fibres and
eight held out.  `interpolate_q80_third_q12_jacobian_mod19_quadratic.sage`
recovers and replays the generic long model, and
`interpolate_q80_third_q12_maps_mod19_quadratic.sage` jointly interpolates
compact bivariate map blocks.  The latter verifies literal generic
function-field identities in both directions; its inverse weighted bounds
are 4 for `W` and 10 for the old cubic coordinate.

`minimize_q80_third_q12_jacobian_mod19_quadratic.sage` clears the unique
Laurent-gauge pole by the exact admissible scaling and produces a short K3
model with coefficient degrees `(8,12)`.  Its degree-24 discriminant has
configuration `I6+I4+3I2+8I1`, hence root lattice `A5+A3+3A1`.  This proves
the exact minimal p=19 child equation and maps.  The component certifier then
transports eleven old D7/D5 curves and the old zero into the full
`A5+A3+3A1` graph, selects the simple-branch component `R5` as the new zero,
and sends the old zero to the `I6` base factor.  This is an exact finite-field
marking; it does not infer characteristic-zero MW rank from Shioda--Tate.
Replay the completed equation, map, marking, and invariant-encoding checkpoint
with:

```bash
python3 elkies-k3/scripts/batch_q80_third_q12_weierstrass_mod19_quadratic.py \
  --attempts 160 --workers 6
sage -python elkies-k3/scripts/interpolate_q80_third_q12_jacobian_mod19_quadratic.sage
sage -python elkies-k3/scripts/minimize_q80_third_q12_jacobian_mod19_quadratic.sage
sage -python elkies-k3/scripts/interpolate_q80_third_q12_maps_mod19_quadratic.sage
sage -python elkies-k3/scripts/certify_q80_third_q12_um2_p19_component_marking.sage
python3 elkies-k3/scripts/encode_q80_third_q12_p19_frobenius_invariants.py
```

The Frobenius encoder writes every parent, child, and map coefficient as
trace plus an anti-invariant multiple of `eta=2r+12`, with `eta^2=18`, and
keeps section sign, base `PGL2`, and Weierstrass scaling in separate ledgers.

The polynomial-closure producer is prime-independent as a finite-field
worker.  Complete children with the same resolved-pencil and minimal-fibre
profiles have been certified at `p=19,61,67,83,89,103,131`.  This local
agreement is not, by itself, a common quadratic lift.

Exact Hensel lifting of both `p=19` closure operands supplies the missing
global interpretation.  After the even substitution `q=l^2`, the five
closure equations have nonsingular Jacobian on both branches.  Newton--Hensel
lifting plus exact reconstruction gives two independent rational square
classes `q1,q2`, with literal replay in all six characteristic-zero closure
equations.  Composing the operands over
`QQ(a,b)`, `a^2=q1`, `b^2=q2`, proves an exact height-eight horizontal: its
`x` coordinate lies in `QQ(a*b)` and its `y` coordinate has `a` and `b`
parts.  All four sign conjugates pass a direct, formerly held-out `p=71`
surface test.  Replay these certificates with:

```bash
sage elkies-k3/scripts/lift_q80_third_q12_closure_operands_p19_qq.sage \
  --biquadratic-operands
sage elkies-k3/scripts/certify_q80_third_q12_biquadratic_horizontal_qq.sage
```

The exact splitting characters explain the earlier alignment puzzle.  At
`19,83,103`, `q1` is inert and `q2` splits; at `61,67,89,131`, the roles are
reversed; `q1*q2` is inert at every one of these primes.  The old local
generator therefore alternated between an `a`-presentation and a
`b`-presentation.  Its trace/norm CRT accumulation remains diagnostic and
must not be interpreted directly as rational reconstruction.

The exact connected pencil supplies the missing common quotient.  Every one
of its 63 coefficients lies in

```text
QQ(theta^2) = QQ(a*b),  theta=a+b,
omega=2*theta^2-2*(q1+q2)=4*a*b,  omega^2=16*q1*q2.
```

`align_q80_third_q12_exact_quadratic_pencil_primes.py` reduces this exact
quadratic field at all seven collected primes.  At each prime it finds the
unique signed scale from `omega` to the local anti-invariant generator and
the unique induced symmetric-square base-`PGL2` transformation for which all
63 moving-equation coefficients replay.  Thus the seven complete local
children are now canonically anchored to one global quadratic pencil; their
child and map coefficients still need to be transported through these base
gauges before a new CRT/LLL reconstruction is valid.  Replay the alignment
with:

```bash
python3 elkies-k3/scripts/align_q80_third_q12_exact_quadratic_pencil_primes.py --check
```

The explicit recovered `2x2` base matrices make the alignment executable.
`transport_q80_third_q12_long_jacobians_exact_quadratic.sage` applies their
inverse substitutions to the six generalized-pipeline long Jacobians and
rewrites every coefficient in the common basis `(1,omega)`.  It verifies the
transport by literal recomputation of the discriminant and `j` invariant.
The resulting seven-prime long-model coefficient schema has 292 slots.  Both
directions of the generic maps also transport with a common schema, producing
3,484 coefficient slots.  Compile the valid CRT ledgers with:

```bash
sage -python elkies-k3/scripts/transport_q80_third_q12_long_jacobians_exact_quadratic.sage --check
python3 elkies-k3/scripts/compile_q80_third_q12_long_jacobian_quadratic_crt.py --check
sage -python elkies-k3/scripts/transport_q80_third_q12_maps_exact_quadratic.sage --check
python3 elkies-k3/scripts/compile_q80_third_q12_maps_quadratic_crt.py --check
```

The CRT modulus is `7739891239523`.  It is not yet large enough for rational
reconstruction; additional complete primes remain to be transported.  The
mandatory legacy p=19 child and its two-way maps are included through a
separate exact generator/base-gauge alignment artifact.

The p-adic route avoids accumulating thousands of small primes.  The exact
pencil compiler writes all 63 coefficients modulo `19^64` directly in the
global `(1,omega)` basis.  The discriminant worker then lifts the intrinsic
`L^3 Q^2 D` factorization and its repeated root through five digits using a
rank-nine fixed finite-field Jacobian inverse.  It exports the normalized
candidate `(z^2+A*z+B)/(LQ)` for the cubic integral basis.  Replay with:

```bash
python3 elkies-k3/scripts/compile_q80_third_q12_exact_pencil_p19_adic.py --check
sage -python elkies-k3/scripts/lift_q80_third_q12_discriminant_factors_p19_adic.sage --check
python3 elkies-k3/scripts/verify_q80_third_q12_integral_basis_mod19_power.py --check
sage -python elkies-k3/scripts/compile_q80_third_q12_riemann_roch_p19_adic_sample.sage --check
python3 elkies-k3/scripts/interpolate_q80_third_q12_long_jacobian_p19_adic.py --check
```

The conductor/root lift and the candidate's generic integrality modulo `19^5`
are certified.  The dedicated verifier uses a single fixed global
`U`-denominator and checks the trace, second symmetric coefficient, and
determinant divisibilities by `LQ`, `(LQ)^2`, and `(LQ)^3` coefficientwise.
The first p-adic Riemann--Roch sample is also complete.  At the exact-gauge
base value `U=16+7*omega`, which reduces to the legacy `T=1` positive control,
the worker flattens the double-branch field to a quartic extension of
`QQ_19`, imposes regularity at both conjugate branches, and descends the
normalized pole-two and pole-three generators back to `QQ_19(omega)`.  It
recovers dimensions `1,1,2,3`, the long Weierstrass equation, and maps in both
directions through all five available digits.  The infinity roots are derived
from each base residue, so the same worker applies to arbitrary good
`U in QQ_19(omega)` rather than only the legacy sample.

Ninety-five accepted residue-distinct samples now support the generic long
equation.  Rational interpolation gives the exact-gauge coefficient degrees
`2/2,4/4,4/4,6/6,8/8` modulo `19^5`; two samples are held out, and every
interpolated coefficient reduces literally to the independently transported
generic p=19 model.  The maps are now generic as well.  Before interpolation,
the worker cancels the p-adically invisible common conductor factors by
solving `n*D-d*N=0` at the transported `W`-degree bounds.  Ninety-three maps
remain canonical (91 training, two held out), including the degree-`40/40`
inverse old-`x` map.  Every scalar map function replays the held-outs and
reduces to the transported generic p=19 map.  None of these modular
certificates is a characteristic-zero claim.

The conductor worker now compresses its rational functions after every Hensel
digit.  This removes the former valuation-four expression-growth wall:
factorization, repeated root, generic integrality, and the complete pinned
child with both maps have replayed through `19^16`.  The exact source currently
supports a next precision run near 60 digits before its `19^64` ceiling.

That ceiling has been removed: an exact `19^260` pencil now supports a fully
certified conductor, basis, and pinned child through `19^256`.  Rational
reconstruction remains underdetermined—successful-looking coordinates are at
the square-root modulus boundary and most coordinates fail—so no exact child
is claimed from that artifact.

The 1,024-digit gate is now complete.  The checkpoint launcher stores
canonical rational-function coefficients at the final target precision and
can resume both the factor and repeated-root stages.  Its pointwise Newton
mode solves the full 9-by-9 system at good `U` values and interpolates on the
certified supports, giving quadratic valuation growth.  The actual hybrid
execution is separately pinned because the immutable canonical worker's
embedded algorithm label describes only the fixed-digit baseline.  Replay
the provenance gate with:

```bash
python3 elkies-k3/scripts/certify_q80_third_q12_p19_precision1024_execution.py --check
```

The exact pencil through `19^1028`, factor/root lift, generic integral-basis
divisibility, and one complete child sample with maps both ways now pass
through `19^1024`.  A batch worker compiles residue-distinct samples with
per-output dependency/hash validation:

```bash
python3 elkies-k3/scripts/batch_q80_third_q12_riemann_roch_p19_adic_samples.py \
  --workers 8 --limit 20
```

Twenty samples are the support-minimal equation batch: seventeen determine
the degree-`8/8` coefficient and three are held out.  This default deliberately
does not consume the full seed inventory; pass a larger `--limit` when a later
calculation needs more residue classes.  The reconstruction worker then
performs projective LLL separately for each long coefficient and requires
reduction to every selected transported model:

```bash
sage -python elkies-k3/scripts/reconstruct_q80_third_q12_long_jacobian_p19_adic.sage
```

The preferred route now reconstructs invariants before the long model and
maps.  The degree-`24/24` `j`-map is composed symbolically modulo `19^1024`
from the five small-support long coefficients, so the existing 17 training
samples plus three held-outs suffice:

```bash
sage -python elkies-k3/scripts/reconstruct_q80_third_q12_j_map_p19_adic.sage
```

`reconstruct_q80_third_q12_j_map_p19_adic.sage` exploits
`numerator = scalar * degree8_polynomial^3` and the denominator multiplicity
shape `I6+I4+3I2+8I1`, and uses three p-adic held-outs.  It supports one joint
lattice across eight residue-distinct evaluations (`--intrinsic-basis
joint-evaluations`), any of the nine projective degree-eight coefficient charts
(`--c4-pivot 0` through `8`), and a coupled two-coordinate projective scale over
the quadratic coefficient field (`--reconstruction-granularity
quadratic-projective`).  Repeatable `--holdout-prime` arguments exclude those
primes from CRT/LLL and require literal replay before an artifact is written.
The output remains a reconstruction candidate, not a characteristic-zero
Jacobian theorem.  The active order is
`j -> (c4^3,Delta) -> minimal Jacobian -> maps`; the common invariant scaling
must be separated from base `PGL2` gauge and Weierstrass scaling.

The complete restructured `19^2048` audit is negative: all coefficient,
component, evaluation, base-normalized, pivot-chart, and quadratic-projective
formulations fail the untouched `p=199` replay.  Their shortest vectors remain
at the random-lattice boundary.  Exact closure operands already have heights
up to 36,335 bits, so `19^8192` would still be too small; the first justified
target is `19^12288`, with a random boundary near 46,400 bits.  The exact pencil
at precision 12291, checkpointed factor/root lift, five-digit generic basis
check, and all 20 residue-distinct Riemann--Roch samples are now complete.  The
sample worker uses compact conductor-power arithmetic for its exact relation
replay; this preserves the certified output while removing the former generic
rational-function normalization bottleneck.  The basis check is reproduced
with:

```bash
python3 elkies-k3/scripts/verify_q80_third_q12_integral_basis_mod19_power.py \
  --source artifacts/local/elkies-k3/q80-third-q12-exact-pencil-p19-adic-precision12291.json \
  --lift artifacts/local/elkies-k3/q80-third-q12-discriminant-factors-p19-adic-precision12288.json \
  --verification-digits 5 \
  --output artifacts/local/elkies-k3/q80-third-q12-integral-basis-mod19-power-lift12288-check5.json
```

Reconstruction at `19^12288` is still negative.  The separate-component and
coupled quadratic-projective degree-eight factors use 46,446-bit primitive
coordinates against a 46,448-bit random boundary, and independent `p=199`
replay rejects both with 25 of 25 numerator pairs and 24 of 25 denominator
pairs unequal.  Since `omega^2=16q1q2` itself has 33,886/33,890-bit
numerator/denominator height, the next route must retain the exact field and
denominator factors intrinsically instead of simply rerunning the same LLL at
larger precision.

The full hashes, failed holdouts, and height justification are recorded in
[`../Q80_THIRD_Q12_COMMON_PRODUCER_2026-09-01.md`](../Q80_THIRD_Q12_COMMON_PRODUCER_2026-09-01.md).

The independent eighth-prime replay is now complete at `p=163`.  The
horizontal, resolved pencil, genus-one gate, 72-sample Jacobian batch, generic
interpolation, quadratic-field alignment, and long-model transport were all
produced before reducing the candidate.  The alignment and transport workers
accept repeatable `--extra-prime PRIME PATH` arguments, avoiding a second
implementation.  In the legacy p=19 base gauge, replay with:

```bash
sage -python elkies-k3/scripts/certify_q80_fixed_u_marked_third_q12.sage \
  --u=-2 --prime 163 --prime-audit-only \
  --output artifacts/generated-results/q80-fixed-u-minus2-p163-good-reduction.json
sage -python elkies-k3/scripts/produce_q80_third_q12_polynomial_closure_modp.sage \
  --surface artifacts/generated-results/q80-fixed-u-minus2-p163-good-reduction.json \
  --system artifacts/local/elkies-k3/q80-third-q12-um2-p163-polynomial-closure.ms \
  --output artifacts/generated-results/q80-third-q12-um2-p163-polynomial-closure-scheme.json
msolve -t 16 \
  -f artifacts/local/elkies-k3/q80-third-q12-um2-p163-polynomial-closure.ms \
  -o artifacts/local/elkies-k3/q80-third-q12-um2-p163-polynomial-closure.solve
sage -python elkies-k3/scripts/certify_q80_third_q12_polynomial_closure_rur_modp.sage \
  --surface artifacts/generated-results/q80-fixed-u-minus2-p163-good-reduction.json \
  --scheme artifacts/generated-results/q80-third-q12-um2-p163-polynomial-closure-scheme.json \
  --solution artifacts/local/elkies-k3/q80-third-q12-um2-p163-polynomial-closure.solve \
  --output artifacts/generated-results/q80-third-q12-um2-p163-common-producer-horizontal.json
sage -python elkies-k3/scripts/compile_q80_third_q12_resolved_pencil_modp2.py \
  --surface artifacts/generated-results/q80-fixed-u-minus2-p163-good-reduction.json \
  --horizontal artifacts/generated-results/q80-third-q12-um2-p163-common-producer-horizontal.json \
  --output artifacts/generated-results/q80-third-q12-um2-p163-resolved-pencil.json
sage -python elkies-k3/scripts/verify_q80_third_q12_resolved_genus_modp2.sage \
  --input artifacts/generated-results/q80-third-q12-um2-p163-resolved-pencil.json \
  --output artifacts/generated-results/q80-third-q12-um2-p163-resolved-genus.json
python3 elkies-k3/scripts/batch_q80_third_q12_weierstrass_modp2.py \
  --input artifacts/generated-results/q80-third-q12-um2-p163-resolved-pencil.json \
  --attempts 72 --workers 4 \
  --sample-dir artifacts/local/elkies-k3/q80-third-q12-p163-weierstrass-samples \
  --output artifacts/generated-results/q80-third-q12-p163-weierstrass-sample-batch.json
sage -python elkies-k3/scripts/interpolate_q80_third_q12_jacobian_modp2.py \
  --input artifacts/generated-results/q80-third-q12-p163-weierstrass-sample-batch.json \
  --output artifacts/generated-results/q80-third-q12-p163-jacobian-interpolated.json
python3 elkies-k3/scripts/align_q80_third_q12_exact_quadratic_pencil_primes.py \
  --p19 artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json \
  --extra-prime 163 artifacts/generated-results/q80-third-q12-um2-p163-resolved-pencil.json \
  --output artifacts/generated-results/q80-third-q12-um2-exact-quadratic-pencils-p19-legacy-aligned-p163-heldout.json
sage -python elkies-k3/scripts/transport_q80_third_q12_long_jacobians_exact_quadratic.sage \
  --alignment artifacts/generated-results/q80-third-q12-um2-exact-quadratic-pencils-p19-legacy-aligned-p163-heldout.json \
  --extra-prime 163 artifacts/generated-results/q80-third-q12-p163-jacobian-interpolated.json \
  --output artifacts/generated-results/q80-third-q12-long-jacobians-exact-quadratic-gauge-p163-heldout.json
python3 elkies-k3/scripts/certify_q80_third_q12_j_map_p163_heldout.py
```

The last command returns
`FAIL_Q80_THIRD_Q12_J_MAP_BLIND_P163_REPLAY`: all 25 numerator pairs and 24
of 25 denominator pairs mismatch.  This rejects, rather than merely weakens,
the p19-adic/auxiliary-prime candidate.  The final rootless Q80 endpoint
`j`-map remains a separate open construction prerequisite for the E29 and
ICARM 398--400 historical-recognition test.

The continuation through `p=191` and `p=199`, together with the intrinsic
evaluation, Galois-component, base-normalization, and factored-`H` probes, is
recorded canonically in
[`../Q80_THIRD_Q12_COMMON_PRODUCER_2026-09-01.md`](../Q80_THIRD_Q12_COMMON_PRODUCER_2026-09-01.md).
Replay the two later blind comparisons with:

```bash
python3 elkies-k3/scripts/certify_q80_third_q12_j_map_p163_heldout.py \
  --prime 191 \
  --candidate artifacts/generated-results/q80-third-q12-j-map-p19-adic-reconstructed-eight-prime-qq.json \
  --output artifacts/generated-results/q80-third-q12-j-map-p191-heldout-replay.json
python3 elkies-k3/scripts/certify_q80_third_q12_j_map_p163_heldout.py \
  --prime 199 \
  --candidate artifacts/generated-results/q80-third-q12-j-map-p19-adic-reconstructed-nine-prime-qq.json \
  --output artifacts/generated-results/q80-third-q12-j-map-p199-heldout-replay.json
```

Both commands return `FAIL`, with 25 numerator-pair and 24 denominator-pair
mismatches.  The generic reconstructor accepts repeatable `--holdout-prime`,
`--reconstruction-granularity`, `--intrinsic-basis`, and
`--base-normalization` options; a held-out mismatch raises before an artifact
is written.  The long-model worker now CRT-combines the selected reconstruction
primes and reconstructs the monic linear factor in `H=(U+r)^2`, but its
nine-prime value of `r` is also rejected at `p=199`.

The long-equation worker remains useful after the invariant gate.  It too
produces only a candidate until literal characteristic-zero replay in the
63-term pencil.  Coefficientwise rational reconstruction remains
inadmissible: its successes still lie at the square-root modulus boundary.

The p=61 compiler is now complete.  The generic mapped-fibre adapter replaces
Sage's characteristic-`p` power-map normalization with Singular's
Grauert--Remmert module normalization and the same reversed-Hermite basis
reduction.  It retains all 72 deterministic samples, interpolates the generic
long Jacobian and maps both ways, minimizes to `I6+I4+3I2+8I1`, and transports
the `R5` zero and `A5+A3+3A1` marking.  Replay the completed child with:

```bash
python3 elkies-k3/scripts/batch_q80_third_q12_weierstrass_modp2.py \
  --input artifacts/generated-results/q80-third-q12-um2-p61-resolved-pencil.json \
  --attempts 72 --workers 6 \
  --sample-dir artifacts/local/elkies-k3/q80-third-q12-p61-weierstrass-samples \
  --output artifacts/generated-results/q80-third-q12-p61-weierstrass-sample-batch.json
sage -python elkies-k3/scripts/interpolate_q80_third_q12_jacobian_modp2.py \
  --input artifacts/generated-results/q80-third-q12-p61-weierstrass-sample-batch.json \
  --output artifacts/generated-results/q80-third-q12-p61-jacobian-interpolated.json
sage -python elkies-k3/scripts/minimize_q80_third_q12_jacobian_modp2.py \
  --input artifacts/generated-results/q80-third-q12-p61-jacobian-interpolated.json \
  --output artifacts/generated-results/q80-third-q12-p61-jacobian-minimal.json
python3 elkies-k3/scripts/align_q80_third_q12_full_children_primes.py
```

The full finite-field alignment certificate has 1,947 ordered coefficient
slots at each prime.  Compile and replay the untransported branch-mixed ledger
with:

```bash
python3 elkies-k3/scripts/compile_q80_third_q12_frobenius_crt_interface.py
```

This formally accumulates the seven local residue sets modulo
`7739891239523`, pins the alternating square-class diagnostic, and explicitly
rejects interpreting the centered CRT integers as common rational
coefficients.  The exact connected correction over the biquadratic field is:

```bash
sage -python elkies-k3/scripts/compile_q80_third_q12_biquadratic_resolved_pencil_qq.py
```

This now pins the exact two-dimensional resolved pencil: Smith degrees
`(0,0,6)`, ambient dimension seven, connected D7/D5 rank five, kernel
dimension two, and 63 moving-equation terms of degrees `(2,9,3)`.  The
117 MB artifact is
`artifacts/generated-results/q80-third-q12-um2-biquadratic-resolved-pencil-qq.json`
(SHA-256 `ac67210166cd414945e1fa373e8f0d5829ff8231daf83c764c376a32ff4b641e`).
The current **ACTIVE_COMPILER** gate is to transport each already certified
finite child and its two-way maps through the exact quadratic generator and
base-`PGL2` alignments, then rebuild the CRT/LLL interface in that common
gauge.  The characteristic-zero child Jacobian, maps, minimal fibres, and
transported marking remain subsequent gates; no characteristic-zero child
equation or Mordell--Weil rank is claimed here.

Certify the generic characteristic-zero genus with:

```bash
sage -python elkies-k3/scripts/verify_q80_third_q12_biquadratic_resolved_genus_qq.sage
```

This reduces the exact 63-term equation at the good prime `19`, proves that
reduction irreducible, and lifts irreducibility to characteristic zero.
Primitivity, completeness of the pencil, Bertini, and K3 adjunction then give
generic genus one.  It does not compute the Jacobian or its maps.

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-LINEAR-CONDUCTOR 957479f39bedd57b -->

Recover and certify the first generic discriminant conductor factor with:

```bash
sage -python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --base-value 0 --certify-generic-linear --check
```

The worker uses pair arithmetic in the exact quadratic descent field, strips
the reconstructed linear factor at `V=0`, and then certifies generically by an
integral expansion in `S=denominator(r)*W+numerator(r)`.  It proves exact
multiplicity three.  Its optional `--attempt-quartic-gcd` is deliberately
gated because the expanded exact degree-four gcd has not completed at the
current multi-million-bit coefficient heights.

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-SPECIALIZED-QUARTICS 725664f9e36ae8a7 -->

The dedicated exact Brown-PRS continuation is:

```bash
sage -python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --attempt-subresultant-prs \
  --output artifacts/generated-results/elkies-k3-q80-third-q12-exact-discriminant-specialization-v1.json
```

It atomically checkpoints every completed remainder to
`artifacts/local/elkies-k3/q80-third-q12-exact-quartic-subresultant-checkpoint-v1.pickle`
and resumes that trusted local checkpoint automatically.  Pass
`--restart-subresultant-prs` only to deliberately start a fresh sequence.  The
completed run returns `PASS_EXACT_SPECIALIZED_L3_Q2_D_FACTORIZATION`: the exact
degree-four gcd has maximum coordinate height 320,859 bits, and literal exact
division gives the remaining quartic at 1,735,258 bits.  Replay the completed
checkpoint and artifact with the same command plus `--check`.  This proves the
factorization at `V=0`, not the generic quartics over the full base.

Replay the strongest retained normalization for that quartic with:

```bash
sage -python \
  elkies-k3/scripts/audit_q80_third_q12_quartic_denominator_candidate.sage \
  --check
```

This reconstructs the common monic linear denominator `H(V)` from the
`19^12288` lift, rebuilds the exact-pencil discriminant independently at
163, 191, and 199, and verifies the predicted denominator on every
nonleading coefficient of the exponent-two quartic.  It remains candidate
data for the generic `V`-dependent quartic; exact `Q^2` divisibility is now
separately proved only at `V=0`, and no Jacobian is claimed.

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-GENERIC-QUARTICS aa704dc4685e4c9b -->

Recover the generic quartic without another LLL reconstruction by lifting the
exact repeated factor through the first `V`-adic jet:

```bash
sage -python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --attempt-generic-quartic-division \
  --output artifacts/generated-results/elkies-k3-q80-third-q12-exact-generic-quartic-factorization-v1.json
```

The worker determines `dQ/dV` uniquely from the exact `V=0` factorization,
replays all sixteen cleared numerator coordinates modulo `19^12288`, expands
the full characteristic-zero discriminant, and proves zero fraction-free
remainder after division by the reconstructed quartic square.  The common
`H` denominator cancels completely from the complementary quartic.  This
certifies the generic discriminant factors; it does not yet construct the
Jacobian or its birational maps.

### Rootless J2 classification controls

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-COMPLETE c6f054948b04b507 -->
<!-- status-consumer: EC-K3-H3-ROOTLESS-J1-UNIFORM-BOUND b71330a75ad2c9ad -->

- `audit_rootless_j2_completeness_track.sage` certifies the two non-isometric
  target-genus controls, rejects the mislabeled old neighbour corpus, and
  computes the exact genus mass obstruction.  It does not claim completeness.
- `certify_rootless_j2_niemeier_controls.sage` pins the correct rank-seven
  anti-discriminant auxiliary, constructs all cyclic discriminant gluings for
  both controls, identifies every ambient as `N(2A7+2D5)`, and exports exact
  primitive embedding and saturated-complement certificates.  Its norm-12
  cost probe explains why the completed all-Niemeier search quotients by
  stabilizer and Weyl actions.
- `build_rooted_niemeier_catalog.sage` certifies hash-pinned Gram models for
  all 23 rooted Niemeier lattices and applies the exact `D5` root-system gate;
  the Leech lattice is excluded separately because the auxiliary has roots.
- `enumerate_niemeier_d5_anchor_orbits.sage` reduces the thirteen admissible
  rooted classes to sixteen full-automorphism `D5` anchor orbits, including
  exact exceptional `E6`, `E7`, and `E8` subsystem enumeration in `--deep`
  mode.
- `enumerate_niemeier_auxiliary_sixth_dominant.sage` enumerates all 3,220
  primitive norm-12 sixth-vector representatives modulo the residual Weyl
  groups.
- `classify_rootless_j2_niemeier_first.sage` applies the strict-positive
  Dynkin-label rootlessness criterion, exactly enumerates 167 remaining
  rational fixed-space ellipsoids, and deduplicates twelve primitive rootless
  embedding representatives to exactly two frame classes: published R17 and
  alternate Q80.  This completes `J2` frame-isometry classification while
  deliberately not asserting that its retained cover counts are embedding-
  orbit or `J1` counts.
- `certify_rootless_j1_uniform_bound.py` is the dependency-free finite-form
  replay for Corollary H2a.  It locks the two complete rootless `J2` inputs,
  enumerates all eight isometries of the cyclic determinant-948 quadratic
  form, and verifies the four cosets modulo the rank-three Hodge image
  `{+1,-1}`.  It proves `2 <= #J1_rootless <= 8`, not the exact `J1` count.

The current proof boundary and replay commands are in
[`../ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md`](../ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md).

### Picard-19 lattice foundry

- `build_rank7_auxiliary_catalogue.sage` is the surface-first merge layer for
  the determinant-banded factory. It imports exact backend records, groups
  first by `(T,NS)`, then by partner auxiliary and frame isometry, retains all
  primitive ambient embeddings, and emits four determinant bands across 23
  rooted backends plus a separate Leech/`Co0` backend. Its current artifact has
  827 surface classes, 1,074 auxiliary classes, and 1,840 MW12--17 frame
  classes. All 96 backend-band
  shards remain explicitly open: this script is not itself an embedding
  enumerator and does not promote bounded foundry inputs to completeness.
- `build_rank7_t_arithmetic.sage` is the mandatory pre-solver arithmetic
  layer. For all 827 `(T,NS)` rows it computes primitive similarity data,
  rational isotropy, the even Clifford quaternion algebra, and the integral
  Clifford-order reduced discriminant. The current ledger has 550 isotropic
  and 277 anisotropic rows. Exactly 313 isotropic rows have a certified
  integral `U` splitting and therefore an exact `Gamma_0(N)` curve with cusps,
  genus, a genus-zero Hauptmodul lookup, and a cusp-width comparison against
  modular elliptic surfaces. The other 237 split orders retain typed-open
  signatures. On the anisotropic side, H3 is the exact `(D,N)=(6,79)` positive
  control with its genus-13 base curve, genus-two `w_474` quotient, and CM
  anchors of discriminants `-3` and `-24`; 171 further rows have only an
  Eichler-level candidate and 105 are certified non-Eichler at a ramified
  prime. Reduced discriminant alone is never promoted to an Eichler theorem.
- `build_rank7_all_niemeier_factory.py` is the supported orchestration entry
  point. It runs catalogue merge, T-arithmetic, and Pareto enrichment in that
  order and deliberately stops before equation solving.
- `build_rank7_surface_pareto.py` ranks all imported surfaces by the exact
  core metrics available for every row: maximum catalogued MW rank, easiest
  known source MW rank, source reducible-fibre support count, and determinant.
  It attaches exact pole, stabilizer, rootless short-vector, and bounded
  multisection evidence where existing artifacts provide it, and emits
  coverage-restricted enriched frontiers without imputing missing equations,
  fields, routes, or unresolved arithmetic geometry. Every row now carries
  the hash-matched T-arithmetic pre-solver gate. The current core
  frontier has fourteen of 827 surfaces; there are 39 exact-pole and 787
  nontrivial-symmetry rows, while the certified-route frontier is explicitly
  empty.
- `enumerate_24a1_octad_prefix_orbits.sage` induces the exact `M24` action on
  all 759 Golay octads and uses set-stabilizer augmentation plus exact GAP
  transporters. It proves the unordered-octad-subset orbit counts
  `1,3,16,206,10547` through size five, including the independent
  orbit-stabilizer mass identities. This is canonical-prefix infrastructure,
  not a rank-seven determinant-band census.
- `enumerate_24a1_octad_rank7_completion_shard.sage` completes a declared
  half-open range of those five-prefix orbits by two octads. It applies the
  union/MW gate, exact modular saturation prefilters, integral
  Construction-A saturation, prefix-stabilizer quotienting, residual-`M24`
  transporter deduplication, saturated complements, root counts, and ternary
  discriminant-form gates. The pinned 43 contiguous shards exhaust all 10,547
  five-prefix orbits and retain 3,051 local residual-M24 records; 3,015 have
  matching ternary genera. The separate multi-shard full-Weyl quotient closes
  the sign action and cross-shard duplicates for this complete positive
  seven-octad input; signed/non-octad generator languages remain open.
- `run_24a1_octad_rank7_completion_frontier.py` reproducibly runs or
  byte-checks that whole contiguous frontier, with an optional bounded number
  of parallel subprocesses. It orchestrates the exact shard enumerator; the
  manifest remains the gap/overlap and hash authority.
- `build_24a1_octad_completion_manifest.sage` discovers the exact completion
  shards, validates a gap-free and overlap-free prefix interval starting at
  zero, and pins every artifact hash and local accounting value. It is the
  shared input contract for the full-Weyl quotient and catalogue.
- `canonicalize_24a1_weyl_m24_shard.sage` applies the missing full
  `W(24A1) semidirect M24` quotient to that completion shard. It replaces the
  impossible enumeration of `2^24` signs by an exact comparison of the 24
  doubled physical coordinate covectors modulo sign, intrinsic auxiliary
  isometries, and GAP `M24` transporters. It accepts contiguous completion
  shards and preserves local provenance. The combined 3,051 records give five
  intrinsic auxiliary classes and 24 full embedding orbits; 18 pass the
  ternary-genus gate and are imported `(T,NS)` first into the catalogue. Every
  collapse has an explicit row-isometry,
  coordinate-permutation, and coordinate-sign witness, and every full orbit
  has a certified stabilizer order.
- `build_cross_niemeier_mod2_priority.sage` turns the umbral stabilizer idea
  into a 23-backend scheduling ledger. It gives first priority to the exact
  `A7^2 D5^2` classes `2B,2C,4A`, then to backends with order-four or swap
  component-permutation envelopes. Multiplicity is only a heuristic: a
  future embedding is retained by this experiment only after its full
  ambient stabilizer induces an action satisfying
  `rank_GF2(g_M-I)>0` on the complement modulo two. It imports the exact
  negative all-class fixed-coordinate control for `3E8` and keeps
  non-coordinate transposition- and three-cycle-invariant languages open. It
  also imports the positive all-class `3D8` control with 40 residual orbits
  and seven local surfaces, plus the exact negative 792-coordinate `2D12`
  swap control, 11,440-coordinate coupled `D10+2E7` involution control, and
  792-coordinate `2A12` order-two control. The positive `2A9+D6` control
  contributes 32 residual orbits and five local surfaces. The positive `3A8`
  control contributes 189 residual orbits and 25 local surfaces.
- `enumerate_2a7_2d5_4a_fixed_rank7.sage` closes the first exact
  symmetry-first family. It enumerates every primitive corank-one rank-seven
  sublattice of the common `4A` fixed lattice through determinant 5,000 using
  the dual-lattice determinant identity. All 336 embeddings have literal
  `2B,2C,4A` stabilizers and nontrivial complement actions modulo two; they
  reduce to three auxiliary and five frame isometry classes, including 304
  rootless/MW17 embeddings. The ternary gate rejects all of them because
  their discriminant groups have length seven. This is a complete negative
  result for the declared pointwise-fixed family, not for nontrivial
  seven-dimensional `4A` representations.
- `enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage` is the first positive
  symmetry-first factory shell. It tests all 11,440 rank-seven coordinate
  direct summands of a pinned integral LLL basis of `Fix(2C)`, filters by
  determinant 5,000, discriminant length three, and MW rank 12--17, and then
  closes the 97 survivors under the exact `Dih_4` section. All 97 section
  orbits pass the nontrivial mod-2 `2C` gate; 73 distinct discriminant forms
  each have one matching ternary genus. `(T,NS)`-first deduplication yields
  73 surfaces, 76 partner auxiliaries, and 86 frames. This is exact in the
  declared coordinate shell, not a complete fixed-lattice or full-Weyl
  enumeration.
- `enumerate_4d6_swap_fixed_high_mw_seed.sage` lifts an exact `S4` component
  section in `N(4D6)` and tests all 11,440 coordinate rank-seven summands of
  a pinned LLL basis of a component-transposition fixed lattice. Exactly 183
  fail discriminant length and all remaining 11,257 have MW rank below 12.
  This closes one bounded negative shell, not the `4D6` backend.
- `enumerate_6a4_double_swap_fixed_high_mw_seed.sage` enumerates all 240
  chamber-preserving component/diagram automorphisms of `N(6A4)` before
  searching the fixed lattice of a literal double component swap. The 11,440
  coordinate summands yield 161 MW12--13 seeds. Full residual-group
  canonicalization, exact ternary gates, and `(T,NS)`-first deduplication give
  42 source surfaces, 55 auxiliary classes, and 128 frames. Three surfaces
  overlap `N(2A7+2D5)`, so the global catalogue gains 39 surfaces. This is
  exact in the declared coordinate language and residual quotient, not a
  full fixed-lattice or Weyl census.
- `enumerate_4a5_d4_order4_fixed_high_mw_seed.sage` enumerates all 48
  chamber-preserving component/diagram automorphisms of `N(4A5+D4)` and
  selects a lifted order-four double component swap.  Of its 120 coordinate
  rank-seven summands, 39 survive with MW rank 13--17 and nontrivial mod-two
  complement action.  Full residual canonicalization and exact ternary gates
  give nine surfaces, nine auxiliary classes, and nine rootless frames, with
  post-dedup distribution `MW13:7, MW17:2`.  This is exact only in the
  declared coordinate language and residual quotient.
- `probe_4a6_4e6_residual_sections.sage` exhausts the 384 possible
  component/diagram maps in each of `N(4A6)` and `N(4E6)`.  The exact
  residual groups have orders 24 and 48, component images `A4` and `S4`, and
  seven and eight matrix conjugacy classes respectively.
- `probe_4a6_4e6_fixed_coordinate_shells.sage` tests all 26,064 coordinate
  summands from every nonidentity residual class of fixed rank at least
  seven.  It finds 86 `4A6` and 48 `4E6` high-MW/mod-two seeds, all from
  order-three classes; the involution, order-four, and order-six coordinate
  shells are negative at or before the MW12 gate.
- `canonicalize_4a6_4e6_fixed_coordinate_shells.sage` closes those seeds
  under the complete residual groups and applies the ternary and
  `(T,NS)`-first gates.  It gives nine `4A6` surfaces with ten frames
  (`MW12:2, MW13:8`) and one `4E6` MW12 surface/frame.  All ten are new in the
  global catalogue.  This is exact for the declared coordinate languages,
  not a full fixed-lattice or Weyl census.
- `probe_8a3_glue_code_residual_group.sage` recovers the order-256 glue code
  `N(8A3)/A3^8` in `(Z/4)^8`, exhausts all `2^8 8!` signed component maps,
  and certifies the complete order-2,688 integral residual group, its
  order-1,344 component image, and 16 matrix conjugacy classes.
- `probe_8a3_fixed_coordinate_shells.sage` exhausts 24,600 rank-seven
  coordinate summands from the seven nonidentity residual classes of fixed
  rank at least seven. Exactly 1,166 pass the determinant, length, MW12--17,
  and nontrivial mod-two gates.
- `canonicalize_8a3_fixed_coordinate_shells.sage` quotients those seeds under
  all 2,688 residual elements. A mod-251 RREF rejects impossible equalities;
  all candidates are rechecked by integral HNF. The 1,162 orbits give 435
  local `(T,NS)` surfaces, 523 auxiliaries, and 574 frames. Twenty-four
  surfaces overlap earlier backends, so the global catalogue gains 411.
- `probe_6d4_hexacode_residual_group.sage` recovers the order-64 glue code
  `N(6D4)/D4^6` in `((Z/2)^2)^6`, exhausts all `6^6 6!` local-triality and
  component-permutation maps, and certifies the complete order-2,160
  residual group, its kernel of order three, full `S6` component image, and
  16 matrix conjugacy classes.
- `probe_6d4_fixed_coordinate_shells.sage` scans all 25,416 coordinate
  summands attached to the eleven nonidentity residual classes of fixed rank
  at least seven. Exactly 472 pass the determinant, length, MW12--17, and
  nontrivial mod-two gates.
- `canonicalize_6d4_fixed_coordinate_shells.sage` closes those seeds under
  all 2,160 residual elements and applies the exact ternary/T-NS gates. The
  466 orbits give 218 local surfaces, 255 auxiliaries, and 289 frames. Fifty
  surfaces overlap prior backends, so the global catalogue gains 168.
- `probe_3d8_glue_residual_group.sage` recovers the order-eight
  `N(3D8)/D8^3` glue code, tests all 48 component/diagram maps, and certifies
  the natural order-six component-permutation `S3` with trivial diagram
  kernel.
- `probe_3d8_fixed_coordinate_shells.sage` scans all 11,448 coordinate
  summands for the transposition and three-cycle classes. The transposition
  shell gives 40 high-MW/mod-two seeds (`MW12:28, MW13:12`); the three-cycle
  shell is negative.
- `canonicalize_3d8_fixed_coordinate_shells.sage` closes those seeds under
  the residual `S3` and applies the exact ternary/T-NS gates. Its 40 orbits
  give seven local surfaces, auxiliaries, and frames. One surface overlaps
  the existing `2A7+2D5`/`6A4` class, so six surfaces but all seven
  auxiliaries and frames are new globally.
- `probe_2d12_glue_residual_group.sage` recovers the order-four
  `N(2D12)/D12^2` glue code, tests all eight component/diagram maps, and
  certifies the natural order-two component-swap residual group with trivial
  diagram kernel.
- `probe_2d12_fixed_coordinate_shells.sage` scans all 792 rank-seven
  coordinate summands of the rank-twelve swap-fixed lattice. None passes the
  determinant, length, MW12--17, and nontrivial mod-two gates, making residual
  canonicalization and catalogue import vacuous for this declared shell.
- `probe_d10_2e7_residual_group.sage` recovers the mixed root decomposition
  and tests all four `D10`-diagram/`E7`-component chamber maps. Only the
  simultaneous diagram involution and component swap lifts integrally,
  giving the complete order-two residual group.
- `probe_d10_2e7_fixed_coordinate_shells.sage` scans all 11,440 rank-seven
  coordinate summands of the rank-sixteen fixed lattice. None passes the
  determinant, length, MW12--17, and nontrivial mod-two gates, so the declared
  coordinate shell contributes no catalogue row.
- `probe_2a12_residual_group.sage` recovers the index-thirteen root quotient,
  tests all eight component/diagram maps, and certifies the complete cyclic
  order-four residual group with fixed ranks `24,12,6,6`.
- `probe_2a12_fixed_coordinate_shells.sage` scans all 792 coordinate summands
  for the sole nonidentity class of fixed rank at least seven. None passes the
  determinant, length, MW12--17, and nontrivial mod-two gates, so no residual
  canonicalization or catalogue import is needed for this declared shell.
- `probe_2a9_d6_residual_group.sage` recovers the index-twenty root quotient,
  tests all sixteen component/diagram maps, and certifies the complete cyclic
  order-four residual group with fixed ranks `24,16,10,10`.
- `probe_2a9_d6_fixed_coordinate_shells.sage` scans all 11,680 coordinate
  summands for the three nonidentity residual classes. The order-two shell is
  negative, while each inverse order-four class gives 32 qualified seeds.
- `canonicalize_2a9_d6_fixed_coordinate_shells.sage` quotients the 64 seeds
  by the exact residual `C4` and applies the ternary/T-NS gates. The 32 orbits
  give five local surfaces, auxiliaries, and frames (`MW13:4, MW17:1`). One
  surface overlaps `8A3`, so four surfaces but all five auxiliaries and frames
  are new globally.
- `probe_3a8_residual_group.sage` recovers the index-27 `A8^3` root quotient,
  tests all 48 component/diagram maps, and certifies the complete order-twelve
  residual group `{+/-1} x S3`, including its six matrix conjugacy classes.
- `probe_3a8_fixed_coordinate_shells.sage` scans all 13,032 coordinate
  summands for the four nonidentity classes of fixed rank at least seven. The
  transposition class gives 189 qualified seeds (`MW12:135, MW13:54`); the
  central involution, signed transposition, and three-cycle shells are
  negative.
- `canonicalize_3a8_fixed_coordinate_shells.sage` closes the seeds under the
  exact residual group and applies the ternary/T-NS gates. All 189 orbits pass
  and give 25 local surfaces, 30 auxiliaries, and 64 frames. Twenty surfaces
  overlap earlier backends, so five surfaces, twelve auxiliaries, and all 64
  frames are new globally.
- `probe_12a2_ternary_golay_residual_group.sage` intrinsically recovers the
  index-729 `A2^12` root quotient as the self-dual ternary Golay `[12,6,6]`
  code and certifies its full order-190,080 monomial group `2.M12`, its
  order-95,040 component image, and all 26 conjugacy classes.
- `probe_12a2_fixed_coordinate_shells.sage` scans all 13,968 coordinate
  summands in the nine nonidentity classes of fixed rank at least seven and
  retains 237 MW12--17 seeds.
- `canonicalize_12a2_fixed_coordinate_shells.sage` checks all 190,080 group
  images with an exact mod-two rejection filter and integral-HNF equality
  certificates. The 214 residual orbits give 210 K3-compatible orbits, 99
  local surfaces, 108 auxiliaries, and 151 frames. Globally, 52 surfaces, 86
  auxiliaries, and 143 frames are new.
- `probe_eta_only_niemeier_residual_groups.sage` intrinsically recovers the
  roots and quotients of the six remaining rooted systems and exhausts every
  product of component diagram maps. It proves trivial residual groups for
  `D24` and `D16+E8`, and order-two eta groups for `A24`, `A17+E7`,
  `A15+D9`, and `A11+D7+E6`.
- `probe_eta_only_niemeier_fixed_coordinate_shells.sage` scans all 35,112
  eligible coordinate summands for the four nontrivial eta classes. `A24`
  fails discriminant length uniformly; every length-admissible mixed-system
  frame has MW rank below 12. Hence all four declared coordinate languages
  are exact negative controls.
- `probe_3e8_residual_group.sage` recovers the three `E8` root components
  intrinsically and certifies the complete chamber residual `S3` as all six
  component permutations. Its transposition and three-cycle classes have
  fixed ranks 16 and 8 and fixed determinants 256 and 6,561.
- `probe_3e8_fixed_coordinate_shells.sage` scans all 11,448 rank-seven
  coordinate summands for those two nonidentity classes. None passes the
  determinant, length, MW12--17, and nontrivial mod-two gates, so no residual
  canonicalization or catalogue import is required for this declared shell.
- `build_leech_co0_backend.sage` derives the invariant integral Gram matrix of
  the AtlasRep 24-dimensional `2.Co1=Co0` representation. It certifies an even
  unimodular rank-24 lattice of minimum four with 196,560 minimal vectors,
  pins the two exact `Co0` generators, and marks rank-seven embedding-orbit
  enumeration open. This is the separate Leech ambient backend, not a Leech
  determinant-band census.
- `build_leech_minimal_line_action.sage` recovers all 98,280 antipodal
  norm-four pairs and the exact transitive `Co1` permutation action. Its line
  stabilizer has the `Co2` order and suborbits `1,4600,46575,47104`, identified
  intrinsically by absolute inner products `4,2,0,1`.
- `probe_leech_minimal_basis_coordinate_shell.sage` certifies a norm-four
  determinant-one ambient basis and exhausts all 346,104 of its rank-seven
  coordinate summands. Exact signed-basis canonicalization gives 221 types;
  207 pass discriminant length, 194 pass the ternary-genus gate, and
  `(T,NS)`-first deduplication gives 150 preliminary MW17 surface keys. This
  is an exact declared coordinate language before the `Co1` quotient, not an
  all-primitive determinant-band census. A 24-bit superset transform also
  records the complete norm-four-pair distribution of every complement,
  ranging from 931 to 2,160, for short-vector ranking and Co1 bucketing.
- `compare_leech_coordinate_targets_to_rooted_catalogue.py` compares those
  150 keys literally with the rooted catalogue. It finds 43 matches whose
  easiest currently catalogued rooted frames have MW-rank distribution
  `MW12:25, MW13:11, MW14:4, MW17:3`; 107 keys are absent from the current
  catalogue. It does not turn a rooted lattice frame into an equation or a
  certified physical neighbour corridor.
- `build_lattice_foundry.sage` consumes the complete H3 `J2` controls and the
  hash-pinned Niemeier catalogue, then runs a JSON-declared auxiliary mutation
  shell. It saturates every rank-seven auxiliary, retains full ambient
  markings, computes and isometry-deduplicates all rank-17 complements,
  classifies roots/MW ranks, solves the even ternary signature-`(2,1)`
  discriminant-form gate, and emits target, companion, route, and equation
  ledgers.
- The initial `one-root-control-shell-v1` run is complete only inside its 768
  declared mutations. It is a deterministic discovery shell, not a complete
  determinant-5000 auxiliary classification. Its separate locked regression
  still replays the complete H3 result: published R17 and alternate Q80, with
  no third rootless `J2` frame.
- `enumerate_lattice_foundry_prescribed_root_sources.sage` is the default
  source-discovery workflow. For each selected foundry NS class it fixes the
  stored rank-seven auxiliary, enumerates prescribed rank-15--17 root faces in
  the rooted Niemeier lattices, solves the remaining auxiliary embedding
  conditions exactly, and accepts a source only after primitive-embedding,
  saturated-complement, full-root-system, and genus gates. Exact repetitions
  with the same deterministic reduced Gram are merged, but distinct reduced
  Grams may still represent the same integral-isometry or `J2` class.
  Each source is attached to every catalogued MW15--17 target frame in the
  same NS class; choosing one target frame is not part of the source search.
  Because this construction works with the full auxiliary embedding rather
  than the cyclic discriminant gluing used by the older hunter, it also covers
  the foundry classes with noncyclic discriminant group.
- A prescribed-root run is exact only for its declared root-type, support, and
  Niemeier embedding shell. A retained complement certifies a complex
  lattice-level fibration and its exact root/MW rank; it does not construct a
  rational marking, a characteristic-zero equation, or a physical
  elliptic-neighbour corridor. A miss is not a classification outside the
  declared slice.
- The first complete production slice uses `N(3E8)`, all 48 foundry NS
  classes, and all-`A` root systems of rank 15--17 with two or three supports.
  It finds 97 reduced-Gram MW2 representatives in 23 NS classes: 64 of type
  `A2+A6+A7` and 33 of type `A1+2A7`. The all-ambient NS0001 control also
  recovers `E7+E8/MW2` and the pinned H3 binary height form. These are exact
  slice results, not rational source models or `J2`-deduplicated counts.
- `summarize_lattice_foundry_prescribed_root_shards.py` independently audits
  disjoint prescribed-root shards, checks their input hashes and declared
  MW/root-rank window, and records both shard-local source occurrences and
  exact repeated `(NS, reduced-Gram digest)` identities. It deliberately does
  not identify unequal reduced Grams up to lattice isometry or `J2`.
- The full-support rank-16/17 census covers all 48 foundry NS classes, all 13
  D5-admissible rooted Niemeier ambients, and all 16 stored D5 anchor orbits.
  It records 2,134 reduced-Gram representatives, all MW1, with every NS class
  covered and no MW0 row. There are no repeated exact `(NS, Gram digest)`
  identities across shards. The MW0 miss is exact only inside this declared
  embedding cover and is not a global non-existence theorem.
- `audit_lattice_foundry_rank1_section_poles.sage` solves the exact affine
  closest-vector problem for every primitive-root MW1 row and every multiple
  capable of improving the current norm. It certifies 1,342 rows with exact
  rational recomputation and independent 256-bit MPFR agreement, while 792
  nonprimitive-root rows remain explicitly open for torsion/glue analysis.
  This separates genuinely cheap generators from deceptively simple
  one-support fibre configurations.
- `probe_lattice_foundry_ns0011_source_ansatz_modp.sage` exhausts the
  normalized split semistable fibre ansatz for the selected
  `A2+A6+A8/MW1` source. `build_lattice_foundry_ns0011_pole2_section_modp.sage`
  compiles its exact I9-depth-two and I7-depth-one component data into
  `R(t)^2=H(t)` pole-two charts, and
  `scan_lattice_foundry_ns0011_pole2_sections_modp.sage` evaluates those
  charts exhaustively as polynomial functions. The complete `GF(5)` chart is
  empty after seven fibre models, 19 infinity charts, and 1,484,375 affine
  tuples. This is a local obstruction for the displayed normalized split
  chart, not a characteristic-zero nonexistence theorem; the `GF(7)` result is
  only a bounded fibre-sample pilot.
- `probe_lattice_foundry_ns0007_source_ansatz_modp.sage` targets the cheaper
  `A1+A3+2A6/MW1` alternative, whose exact generator has pole zero and only
  depth-one conditions at I2 and I4. It normalizes the four supports to
  `0,1,lambda,infinity` and imposes all 20 branch jets. The complete split
  `GF(5)` scan covers `3*5^8` normalized A polynomials, finds 966
  Hermite-compatible sign branches, and finds no squarefree
  `I2+I4+2I7+4I1` model. This is an exact obstruction only for that displayed
  characteristic-five chart. The stored `GF(7)` run is a bounded
  100,000-case negative pilot.
- The same probe now also accepts `--candidate ns0034` and optional fixed
  `lambda`/leading-`A` slices.  The primitive `NS0034-S008` source has
  `A2+A3+A4+A7/MW1`, height `19/8`, pole zero, and three rootless MW17
  endpoints.  A complete fixed `GF(7)` slice at `lambda=2,A8=1` checks all
  `7^7` remaining `A` polynomials and has a unique exact squarefree fibre
  model.  `scan_lattice_foundry_ns0034_pole0_sections_modp.sage` exhausts its
  `7^8` component-adapted section tuples and finds eight polynomial sections,
  but none with the required I4/I8 depths and I3/I5 identity components.
  This excludes only that fixed modular slice.
- `build_lattice_foundry_ns0034_fibre_hermite_modp.sage` writes the raw
  28-variable nodal-Hermite system.  The A-eliminated
  `build_lattice_foundry_ns0034_fibre_hermite_reduced_modp.sage` fixes an
  infinity square-root branch and reduces it to 20 cubic-or-lower equations.
  The pinned `p=7,lambda=2,A8=1,hi0=3` solve is positive-dimensional because
  the divisibility locus includes higher-order/collided fibre boundaries;
  exact residual-order testing is still required.
- `scan_lattice_foundry_ns0043_pole0_sections_modp.sage` reuses the complete
  seven-model `I9+I7+I3+5I1` fibre census for the stronger source-first
  `NS0043-S005` marking.  This source has MW1, height four, pole zero, zero
  component corrections, and four same-NS MW15 targets.  Exhausting both
  local quadratic-twist classes at `GF(5)` finds 54 polynomial sections in
  the nonsquare twist and none in the square twist, but no section meets all
  three smooth identity components.  The two stored `GF(7)` fibre models give
  the same bounded two-twist negative result.
- `probe_lattice_foundry_ns0030_source_ansatz_modp.sage` promotes a pole-zero
  semistable source attached to both MW15 and MW16 targets.  It imposes the
  ordered `I3+I2+I2+I7+I7` branch jets with exact corrections
  `0,0,1/2,6/7,10/7`.
  `combine_lattice_foundry_ns0030_source_ansatz_modp.py` certifies that the
  prefix and suffix are adjacent segments of the same coprime-stride
  permutation.  The combined `GF(5)` census therefore covers all
  `6*5^8=2,343,750` normalized `A` polynomials, finds 1,536
  Hermite-compatible signed branches, and no exact prescribed fibre orders.
  This is a complete obstruction only for the displayed normalized chart.
- `probe_lattice_foundry_ns0048_source_ansatz_modp.sage` puts the star fibre
  of the `A1+A4+A6+D5/MW1` source at infinity, reducing the short-model bounds
  to `(6,9)`.  The complete `GF(7)` census has six
  `I5+I7+I2+I1*+3I1` models.
  `scan_lattice_foundry_ns0048_pole0_sections_xonly_modp.sage` exhausts only
  the three remaining X parameters and finds one marked sign pair with height
  `37/14` and determinant 740.  The independent tensor scanner is retained as
  a cross-check.
- `audit_lattice_foundry_ns0048_marked_family_modp.sage` evaluates the full
  22-equation section-built family at the surviving points.  It proves
  Jacobian rank 18 in 19 variables at both `p=7` and a bounded `p=11` point.
  `lift_lattice_foundry_ns0048_marked_family_padic.sage` Hensel-lifts the
  latter through `11^80` after fixing a local parameter, but exact rational
  reconstruction fails.  `build_lattice_foundry_ns0048_tate_family_modp.sage`
  translates the section to `(0,0)` and removes the forced support square from
  the discriminant; this cuts the one-dimensional system from 21,606 to 9,071
  monomials.  None of these modular or p-adic certificates supplies a rational
  source equation or an MW16 corridor.
- `probe_lattice_foundry_ns0028_source_ansatz_modp.sage` tests the
  multisection-leading `A2+A6+A7/MW2` source.  Both exact generators have pole
  zero and only depth-one node contacts.  The complete normalized fibre
  censuses contain 25 models at `GF(5)` and 112 at `GF(7)`.
  `scan_lattice_foundry_ns0028_pole0_section_pairs_modp.sage` exhausts the
  three-parameter P and four-parameter Q X charts and both twist classes.  It
  finds no section at 5; at 7 it finds each individual section type, but never
  on the same fibre model.  This rejects the displayed two-prime charts, not
  the characteristic-zero lattice source.
- `build_lattice_foundry_ns0007_pole0_reduced_modp.sage` eliminates the
  global `a3,a4` coefficients from the pole-zero `NS0007-S025` chart and
  emits its 19-variable, 19-equation fixed-lambda system.
  `run_lattice_foundry_ns0007_p7_fixed_case_census.py` partitions the six
  displayed base-field coefficients into exact lexicographic ranges and uses
  isolated `msolve -g 1` processes to classify unit versus nonunit ideals.
  Only the fully expanded input is certifying: although
  `audit_lattice_foundry_ns0007_compact_msolve_encoding.sage` proves that the
  compact strings expand to the same Sage polynomials, msolve reports false
  nonunit ideals on the factored syntax.  The runner therefore rejects that
  syntax by default.  `combine_lattice_foundry_ns0007_p7_fixed_case_census.py`
  verifies contiguous expanded-shard coverage and fails on gaps, overlaps,
  timeouts, solver errors, or noncertifying syntax.
  `repair_lattice_foundry_ns0007_p7_fixed_case_census.py` replaces only
  declared timeout/error entries with checked terminal singleton replays and
  records both the raw artifact and every replay hash.
- `hunt_lattice_foundry_rootful_source.sage` and
  `run_lattice_foundry_mw3_broad_scout.py` are retained for bounded Kneser
  discovery provenance and byte-for-byte replay of their existing artifacts,
  not as the default source-search strategy. The hunter accepts any exact
  foundry frame, including MW15/MW16 starts, and `--allow-below-target`
  records the best exact bounded source when the requested root rank is
  missed. The original exact promoted `NS0024` sources are
  `5A1+A2+A5/MW5` in `N(A11+D7+E6)` and the more equation-friendly
  `A3+A4+A6/MW4` in `N(A15+D9)`; both are saturated complements of the same
  primitive rank-seven auxiliary.
- `score_lattice_foundry_sources.sage` ranks every stored same-NS source by
  the MW0--2-first equation objective, attaches all catalogued MW15--17 target
  frames, and leaves rational marking, Galois orbit, parametrization, and
  uncertified route costs explicitly unknown. Audited low-degree
  multisection richness is retained only as the final heuristic tie-break. It
  consumes both individual rootful-source certificates and the rows of the
  prescribed-root inventory, without widening the latter's declared finite
  search scope.  For primitive MW2 rows it ranks the exact cheapest complete
  basis maximum pole before the cheapest individual-section pole.
- `audit_lattice_foundry_rank2_section_basis_poles.sage` uses the exact MW
  height Schur complement to bound all tail classes capable of improving the
  displayed basis, solves every affine root-lattice CVP with double-double and
  MPFR-256 GSO arithmetic, recomputes the returned norms exactly, and compares
  all determinant-one tail pairs.  The pinned artifact covers all 97
  primitive MW2 rows in the declared all-A inventory; it is an exhaustive
  tail-class calculation with a stated numerical CVP branch boundary.
- `sample_lattice_foundry_multisection_spectrum.sage` computes complete
  degree-two low-height translation-coset spectra on selected rootless foundry
  frames and deterministic exact-CVP samples in degrees three and four. It
  replays the published R17 count of 39,120 rational bisection orbits; sampled
  higher-degree coordinates are discovery heuristics, not curve censuses or
  rank predictions.
- `analyze_lattice_foundry_umbral_orbits.sage` reconstructs the lambency-eight
  umbral group section for `N(2A7+2D5)`, computes full auxiliary-embedding
  stabilizers across the primitive gluing, and resolves their actions on
  norm-four vectors and rational degree-two cosets. Its degree-three data are
  exact only inside deterministic group-invariant samples; it claims no
  correspondence with an umbral module.
- `complete_lattice_foundry_degree3_spectrum.py` exhausts all `3^17`
  translation cosets for each selected rootless frame, using inversion to
  halve the CVP workload, exact integral norm recomputation for every returned
  candidate, deterministic 256-bit MPFR cross-precision audits, and resumable
  chunk checkpoints. The primary pinned batch is the five surfaces selected
  by the earlier cheapest-single-section MW2 ranking; a second artifact retains the five
  pre-prescribed-root route-aware leaders plus R17. Both lattice censuses are
  complete, but effectivity, nefness, irreducibility, arithmetic descent, and
  specialization rank gain remain open.
- `analyze_r17_multisection_diversity.py` refines the R17 count-only
  coordinate into a complete degree-two quotient metric and graph, a complete
  degree-three one-vertex spectrum, deterministic `Aut(M)`-closed
  degree-three/four graph samples, exact two-to-four torsion overlap,
  equation-complexity weights, and Gauss-local squareclass signatures.  It
  labels representative angles and all higher-degree graph fields by their
  gauge or sampling boundary.  Its `--comparison-only` mode also gives the
  complete common degree-two profile, full PARI automorphism action, and exact
  inherited degree-four mass for R17, `NS0032-F011`, and `NS0028-F005`.
  `--control-calibration-only` then intersects the R17 graph with the exact
  split-cover sets at the rank-21/25/26/27/28 controls and counts every induced
  component and clique together with lattice and exceptional-quotient spans;
  no new rank theorem is inferred.
- `enumerate_golay_det720_prescribed_root_sources.sage` applies the
  prescribed-root search directly to the rootless Golay-octad rank-seven
  auxiliary.  Its complete declared 23-Niemeier, root-rank-15/16,
  one-to-three-support run returns 4,823 distinct reduced-Gram MW1/MW2 source
  rows.  `audit_golay_det720_source_poles.sage` performs the complete
  Smith-quotient pole audit through frame norm eight; the compact
  `build_golay_det720_equation_first_shortlist.py` ledger retains all 177
  semistable rows with a complete basis of pole at most two.  Reduced-Gram
  distinctness is not asserted to be integral-isometry distinctness.
- `classify_golay_det720_ideal_source_isometries.sage` completely classifies
  the 48-row ideal cut (MW2, semistable, at most three supports, pole profile
  `[0,0]`) into three marked integral-isometry classes of sizes 35, 4, and 9.
  Their representatives are precisely the already-tested `3A5`, `A11+A4`,
  and `A3+A4+A8` charts, so the modular marking gates cover the whole cut.
- `probe_golay_det720_a11_a4_source_ansatz_modp.sage` and
  `scan_golay_det720_a11_a4_pole0_pairs_modp.sage` exhaust the normalized
  `I12+I5` chart and its two pole-zero generators at 5 and 7.  Individual
  sections occur, but no complete marked MW2 pair occurs in either twist
  class.  The generalized three-support fibre scanner together with
  `scan_golay_det720_three_support_pole0_pairs_modp.sage` similarly rejects
  the tested `A3+A4+A8` charts and finds 24 marked pairs in the nonsquare
  `GF(7)` `3A5` chart for `G720-S0128`.
- `certify_golay_det720_3a5_marked_gf7_lift.sage` certifies Jacobian rank 45
  in 46 variables and lifts all 55 marked equations through `7^8`.
  `certify_golay_det720_3a5_formal_smoothness.sage` then proves that the ten
  nonpivot section equations are forced by the discriminant/node identity and
  component depths.  The unit minor therefore gives a one-parameter formal
  `Z_7` marked family.  With `--free-parameter-integer 10`, the lift reaches
  a small rational point at precision `7^40`.
- `certify_golay_det720_3a5_source_qq.sage` rationally reconstructs that point
  and replays all 55 equations over `QQ`.  It certifies a split
  `3I6+6I1` model, rational section heights `5/6` and `4`, and a rank-19 NS
  sublattice of determinant `-720`.  The separate
  `certify_golay_det720_3a5_picard19.sage` uses exact `F_p` and `F_(p^2)`
  counts at 17 and 19; the rank-20 reductions have incompatible Artin--Tate
  square classes, proving geometric Picard rank 19.
- `certify_golay_det720_3a5_saturation_rejection.sage` detects an exact
  rational 3-torsion section and an exact rational half of the displayed
  height-four section on that `s6=10` model.  They give an index-six
  enlargement.  Exhausting the discriminant-form isotropic subgroups proves
  that it is maximal, so the full NS determinant is 20, the torsion is
  `Z/3`, and the free MW height Gram is `diag(5/6,1)`.  Consequently this
  exceptionally simple Picard-19 rational point is not the determinant-720
  K3; a replacement rational point must pass torsion and divisibility gates.
- `scan_golay_det720_3a5_rational_parameters.sage` repeats the exact
  rational-reconstruction and torsion/divisibility gates in each of the six
  etale marked `GF(7)` residue disks.  The pinned `|a|,b <= 40` scan tests
  1,478 reduced parameters, reconstructs exactly three points, and rejects
  all three as determinant-20 specializations.  Its bounded miss is not an
  irrationality proof for the remaining formal parameters.
- `compile_section_first_normal_forms.sage` is the reusable marked-equation
  frontend for source searches.  Its MW1 Tate-style chart fixes `P=(0,0)`.
  Its MW2 chart fixes `P=(0,0)`, `Q=(h*r,h^2*s)` and eliminates both section
  equations and the intersection divisor by one polynomial Bezout identity,
  leaving discriminant jets as the closed fibre-tuning equations.  The exact
  regression artifact translates both the rational Golay `3I6` pair and the
  NS0031 model-157 `GF(7)` pair, preserving their degree-two smooth
  intersection divisors and their `I_n` discriminant orders.  These are
  normal-form controls, not new saturation or characteristic-zero source
  certificates.
- `extract_rank7_catalogue_source_search_target.py` adapts any selected
  expanded-catalogue surface, partner auxiliary, and target frame to the
  ordered prescribed-root engine. It refuses a missing or stale T-arithmetic
  ledger and attaches the selected exact/typed-open arithmetic row before the
  target can reach an equation workflow; a typed-open curve identification
  keeps target extraction blocked. The generalized
  `enumerate_golay_det720_prescribed_root_sources.sage` retains its locked
  Golay defaults but also accepts these adapters and dynamic determinants.
  For `K3-04b86146cc6b284b`, the complete six-large-ambient search in the
  declared rank-14--17, one-to-three-support window finds 3,101 MW1--3
  reduced Grams, including three integrally isometric `A3+A4+A9/MW1`
  representatives with trivial torsion, height `5/2`, and pole profile `[1]`.
- `certify_k3_04b_equation_first_promotion.sage` explicitly merges those
  three representatives by integral isometries and promotes the resulting
  semistable-root-pattern MW1 source opposite the rootless MW17 target.  It
  records `T=U(5)+<20>` and keeps rational marking, full semistability,
  equation construction, and the neighbour corridor as open gates.  Its
  attached complete degree-three census has 12,095,162 rational and
  29,878,240 genus-one trisection translation cosets.
- `probe_lattice_foundry_ns0028_source_ansatz_modp.sage` now accepts
  three-support A-type sources of root rank 14--17 and correctly pads local
  A-jets beyond degree eight, as required by an `I10` support at infinity.
  Its exhaustive determinant-500 `GF(5)` chart finds three squarefree
  `I4+I5+I10+5I1` models among all `5^8` normalized A polynomials.
- `scan_k3_04b_a3_a4_a9_pole1_marking_modp.sage` exhausts the pole-one
  section chart on those models.  The seven depth `(2,0,5)` X-jets leave one
  numerator per denominator; the square twist has four signed MW1 sections
  and the nonsquare twist has none.  The two distinct positive seeds have
  rank-38 Jacobians in 40 variables and stop at `5^2` in
  `certify_k3_04b_a3_a4_a9_marked_gf5_hensel.sage`, so they are retained as
  finite-field feasibility evidence rather than characteristic-zero lifts.
  The same scanner exhausts all `7^8` normalized A polynomials at prime 7,
  finding six squarefree fibre models.  Both twist classes are marking-positive
  (two and four signed sections), and all three distinct marked seeds have the
  expected rank-39 Jacobian in 40 variables and lift through `7^8`.  This is a
  smooth local one-parameter precursor.  The exact node/discriminant identity
  in `certify_k3_04b_a3_a4_a9_formal_smoothness.sage` forces the 14 excess
  section equations from the fibre/component jets and the five middle
  residual coefficients.  The unit minor therefore proves a formally smooth
  one-dimensional `Z_7` marked family; rational algebraization remains open.
- `certify_k3_04b_a3_a4_a9_source_qq_rejection.sage` reconstructs the first
  small rational point (`m4=-20`) and rejects it at the saturation gate.  The
  displayed height-`5/2` pole-one section is five times a rational pole-zero
  section of height `1/10`.  Exact discriminant-form enumeration shows that
  the index-five enlargement is maximal, so the rank-19 primitive closure has
  determinant 20 rather than 500.
- `scan_k3_04b_a3_a4_a9_rational_parameters.sage` makes the bounded
  determinant-500 rational probe reproducible.  It tests 87 integral values
  of the free coordinate across all three smooth marked `GF(7)` disks at
  precision `7^40`; 86 have no full coefficientwise reconstruction and the
  sole exact point is matched literally to the determinant-20 rejection
  certificate.
- `classify_k3_04b_semistable_mw2_sources.sage` puts all nine retained
  semistable `A3+A4+A8/MW2` pole-`[0,0]` rows in one integral frame class.
  Their selected bases form two exact marking profiles of sizes three and
  six, but both require depths `(1,0,1)`, `(1,0,3)` and smooth intersection
  one.  `scan_k3_04b_a3_a4_a8_mw2_marking_modp.sage` exhausts those equation
  conditions on the 30 `GF(5)` and 114 `GF(7)` squarefree fibre models found
  by the generalized three-support scanner.  No twist gives a marked basis.
  The nonsquare `GF(7)` chart has 20 and 32 individual generator sections;
  all 28 component-matched pairs meet a reducible fibre.
- `certify_k3_10a_semistable_source_rejection.py` aggregates the complete
  ideal-source searches for both auxiliary classes of the determinant-750,
  genus-zero `Gamma_0(3)` surface `K3-10a14a46c14b3150`.  Its rootless MW17
  target remains exact, but neither auxiliary has a semistable MW0--2 source
  with at most three supports in the six-large-ambient cut.  This is a scoped
  negative result, not a complete fibration classification.
- `search_lattice_foundry_same_ns_compiler_routes.sage` now has a reproducible
  determinant-500 rank-first beam from the `A3+A4+A9/MW1` source to the named
  rootless MW17 frame.  With `q=4`, old-fibre degree two, pole order at most
  one, beam width eight, and depth twelve, it reaches root rank four (MW13)
  after four edges but never the target.  Its empty result is bounded by beam
  pruning and the declared compiler coordinates; it is not a complete
  neighbour-graph rejection.
- `extract_rank7_catalogue_source_search_target.py --lattice-only` separates
  pure same-surface source discovery from equation authorization.  It permits
  an exact auxiliary/target adapter when the T-arithmetic curve gate is still
  open, records `equation_work_authorized: false`, and does not promote the
  unresolved arithmetic row.
- `classify_k3_6ce_a5_a10_mw2_sources.sage` classifies the best determinant-384
  MW2-to-MW15 source cut.  The 24 semistable `A10+A5` rows with two supports
  and basis poles `[0,1]` form one integral frame class and one marking profile:
  depths `(0,4)` and `(0,2)` at `I6,I11`, with required smooth intersection
  one.  `probe_k3_6ce_a5_a10_fibre_ansatz_modp.sage` now reads any two-A-
  component rank-15 source orders from the selected lattice artifact and
  exhausts its normalized two-support fibre chart, while
  `scan_k3_6ce_a5_a10_mw2_marking_modp.sage` derives the local depths and
  required intersection and exhausts its pole-zero/pole-one section marking
  over each selected finite field and twist.  Their pinned defaults remain
  byte-compatible with the determinant-384 calculation.  The fibre scans
  find 152 squarefree models over `GF(5)` and 1,032 over `GF(7)`, but all four
  square/nonsquare marking charts are empty at the full basis gate.  The
  nonsquare `GF(7)` chart has 84 and 228 individual generator sections and 72
  component-matched pairs, all with the wrong smooth intersection.
- `classify_k3_6ce_a2_a5_a8_mw2_sources.sage` and
  `scan_k3_6ce_a2_a5_a8_mw2_marking_modp.sage` test the next determinant-384
  pole-`[0,1]` source class.  Its 17 rows give one frame class and two marking
  profiles.  The nonsquare `GF(5)`/`GF(7)` charts carry both individual
  generator types, but all 28/92 candidate pairs meet singular fibres; both
  profiles are therefore empty in the displayed normalized charts.
- `certify_k3_6ce_equation_first_candidate.py` aggregates that MW2-to-MW15
  discovery with its depth-eight same-NS corridor beam.  The beam reaches
  MW13 after five degree-two edges but not the exact `A2/MW15` target.  The
  certificate promotes determinant 384 only as the new lattice-source leader;
  rational marking, T-arithmetic curve identification, rootful-target
  multisections, and a complete corridor remain open.
- `scan_k3_04b_a3_a4_a8_mw2_marking_modp.sage` now also accepts any
  three-support semistable MW2 source with a pole-`[0,0]` physical basis,
  including repeated A-component orders, and derives the required pair
  intersection from the physical basis lattice.  Its determinant-500 default
  artifact is unchanged.
- `certify_k3_14ad_equation_first_candidate.py` aggregates the determinant-654
  source-first test opposite an `A1/MW16` target.  The six-large-ambient cut
  has 9 abstract MW1 and 411 MW2 rows, with 127 complete low-pole MW2 bases.
  Exhaustive normalized marking scans reject five cheapest profiles: three
  pole-`[0,0]` sources at primes 5 and 7 and two pole-`[0,1]` sources at prime
  5.  The selected pole-zero source reaches MW14 but not the exact MW16 target
  in the capped depth-eight degree-two beam.  These are scoped modular and
  bounded-route negatives, not characteristic-zero or graph obstructions.
- `scan_k3_04b_a3_a4_a9_pole1_marking_modp.sage` and
  `certify_k3_04b_a3_a4_a9_marked_gf5_hensel.sage` now derive arbitrary
  three-support root-rank-16 MW1 fibre orders and component depths while
  retaining byte-compatible determinant-500 defaults.  The marking scanner
  exhausts any affine numerator-X space left by fewer than seven component
  jets; the Hensel checker builds the corresponding dynamic fibre/component
  system.
- `certify_k3_cf7f_a4_2a6_mw1_formal_smoothness.sage` promotes the new
  determinant-714 `A4+2A6/MW1` source opposite an `A1/MW16` target.  Its
  nonsquare `GF(7)` chart has four marked sections on two models; both seeds
  have one-dimensional smooth tangents and lift through `7^8`.  The exact
  node/discriminant identity forces the eight omitted residual equations, so
  the branch is formally smooth over `Z_7`.
- `scan_k3_cf7f_a4_2a6_mw1_rational_parameters.sage` checks 23 small integral
  `m8` values across both smooth det714 residue disks through `7^40`; none has
  a full coefficientwise rational reconstruction.  The bounded miss is not
  an irrationality result.
- `certify_k3_cf7f_equation_first_candidate.sage` aggregates both auxiliary
  source censuses (463 and 548 rows), proves the two MW1 representatives are
  integrally isometric, and pins the modular, formal, rational-scan, and
  corridor gates.  The degree-two beam reaches an MW15 intermediate at depth
  six but not the exact MW16 target through depth eight.  Det714 is the active
  second formally smooth MW1 promotion after det500.
- `build_rank7_rational_moduli_source_optimizer.py` maintains the evaluated
  source-first queue for the six catalogue surfaces having rational moduli and
  a rootless MW17 frame.  Determinant 500 is the active formally smooth MW1
  promotion; all other five surfaces have exact negative results in the
  declared semistable MW0--2 six-large-ambient cut.  Each determinant-1296
  auxiliary has 402 terminal embeddings in the three A-containing ambients;
  determinants 1500 and 1728 have 120 and 768 respectively.  Every terminal
  is nonprimitive.  The completed MW17 queue now points explicitly to a
  rational-moduli MW15/MW16 expansion.
- `sample_lattice_foundry_multisection_spectrum.sage` and
  `complete_lattice_foundry_degree3_spectrum.py` also accept the expanded
  rank-seven catalogue's `surfaces` schema.  The determinant-500 target has
  29,040 rational and 63,895 genus-one bisection orbits; its complete
  degree-three census is checkpointed and byte-checkable separately.
- `build_golay_det720_foundry_adapter.py` exposes the certified rootless MW17
  Gram to the generic spectrum programs without pretending that it belonged
  to the original consolidated foundry search.  Its complete degree-three
  census has 15,717,830 rational and 31,988,690 genus-one trisection
  translation cosets.  These are below the NS0031 leaders, so equation cost
  and multisection richness remain separate optimization coordinates.
- `search_golay_det720_degree2_source_corridor.sage` exhausts all 64,515
  section-nonnegative integral genus-one degree-two classes on that rootless
  target: 64,355 have coset minimum eight and 160 have minimum twelve.
  Exactly 64,512 define primitive elliptic fibres.  Minimum-eight children
  have root system `rA1` with `1 <= r <= 11`, and every minimum-twelve child
  is rootless, so none is the marked `3A5/MW2` source.  This closes the full
  one-edge degree-two corridor, not higher-degree or multi-edge routes.
- `certify_lattice_foundry_route.sage` consumes an ordered route manifest,
  replays every primitive isotropic split, checks component/all-section and
  Proposition-C2 finite horizontal walls, composes determinant-one NS
  markings, and scores the terminal rootless frame. The two pinned routes have
  eleven and thirteen degree-two edges, use only `q=4,6`, and require zero
  physical Weyl repairs. The MW5 route discovers a fourth rootless `NS0024`
  class; the MW4 route lands on catalogue frame `NS0024-F005`.
- `probe_lattice_foundry_ns0024_source_ansatz_modp.sage` imposes the
  `I7+I5+I4+8I1` branch jets of the MW4 source over a finite field. Exact
  examples at 11, 13, and 17 prove fibre-stratum feasibility only; the four
  MW sections and `NS0024` marking remain open gates.
- `prepare_lattice_foundry_ns0024_edge1_compiler.sage` certifies the abstract
  q4/orbit1 handoff: horizontal `P3`, divisor
  `O+P3+2F-C2-2C3-C4`, four-dimensional chord ambient, and expected resolved
  rank two. `compile_lattice_foundry_ns0024_edge1_modp.sage` consumes a
  certified model over `GF(p)`, `GF(p^d)`, or a one-parameter function field
  and fails closed unless the exact resolved two-plane, quartic, and
  `A1+A2+A4+D5` child gates all pass.
- `recover_lattice_foundry_ns0024_mw4_family_resolved_modp.sage` now exports a
  joint resolved-depth13 MW3/MW4 ideal.  Surface and all marked-section
  coordinates remain algebraic, so a closed point may place the entire MW3
  marking over `GF(p^d)`.  The forced Q3 section, exact-depth saturations,
  identity-component open charts, generic surface hyperplanes, and MW3-only
  diagnostic mode are explicit.  `--explicit-formal-centers` retains the nine
  `I5/I4` center jets with sparse exact recurrences and avoids recursive
  inverse-power expansion; bounded msolve F4 batches are the reproducible
  low-memory route for the resulting joint system.
- `extract_lattice_foundry_ns0024_joint_gb_point.sage` decodes a
  zero-dimensional joint basis over an arbitrary irreducible residue-field
  modulus, scans the Frobenius-field points, and emits only after replaying all
  four equations, fibre orders, absolute profiles, and the full Gram.  Its
  arithmetic-realizability audit separates the surface, four-section, and
  auxiliary orientation orbits; on each surface-fixed Frobenius power it
  certifies the integral action and fixed rank in the marked MW4 lattice.
- `extract_lattice_foundry_ns0024_joint_rur_point.sage` is the primary
  arbitrary-degree path.  It requires the exporter's fixed full-coordinate
  RUR anchor, a squarefree degree-equals-quotient eliminant, exact substitution
  in the original joint system, and factor-local Frobenius decoding; each
  factor is then passed through the independent joint-GB source verifier and
  retained in a compact closed-point/Galois audit.  The first exact marking is
  sent through the adapter and edge-1 compiler automatically; `--no-edge1`
  exists only for extractor diagnostics.
- `adapt_lattice_foundry_ns0024_mw4_point_for_edge1.sage` is the lossless
  bridge from the residue-algebra recovery format. It accepts either a compact
  P4 point plus an indexed prime-field MW3 seed, or a joint point embedding the
  surface and all three MW3 sections over the same finite extension, and
  independently replays every source gate before granting compiler-input status.
  Both the original minimum-pole basis and the exact q4-containing resolved
  component-depth recommendation are supported; the latter binds `P2` as the
  same abstract q4/orbit1 horizontal class.
- `convert_lattice_foundry_ns0024_mw4_seed_to_point.py` losslessly converts a
  direct prime-field `MW4SEED` record into that same compact point schema.  It
  deliberately grants no geometric status beyond conversion; the adapter is
  still the independent four-section/profile/Gram gate.
- `run_lattice_foundry_ns0024_edge1_handoff.sage` is the single-command entry
  point. It dispatches a certified family directly to the compiler, or runs
  the marked-point adapter first when given the compact residue-algebra format;
  `--check` replays both generated artifacts without rewriting them.

The scope, counts, certified route, and open equation gates are recorded in
[`../LATTICE_FOUNDRY_REPORT_2026-09-01.md`](../LATTICE_FOUNDRY_REPORT_2026-09-01.md).
The source-first objective, bounded high-rank-frame trials, and multisection
score are recorded in
[`../LATTICE_FOUNDRY_SOURCE_FIRST_OBJECTIVE_2026-09-01.md`](../LATTICE_FOUNDRY_SOURCE_FIRST_OBJECTIVE_2026-09-01.md).

### Complete bisection pair-cover arithmetic

- `analyze_elkies_2026_bisection_pair_covers.sage` proves all 39,120
  individual conics rational, all branch quadratics irreducible and distinct,
  and all 765,167,640 distinct `V4` pair bases genus one with exact surface
  height matrix `diag(24,24)`.
- `catalogue_elkies_2026_immediate_point_pairs.sage` computes exact minimal
  Jacobians, conductors, and global root numbers for all 5,566 pair bases with
  a visible rational point at zero or infinity.
- `screen_elkies_2026_immediate_pair_ranks.sage` maintains a resumable bounded
  point ledger for either the immediate-point or control-selected catalogue.
  Only exact finite-quotient-certified independent points count toward its
  lower bounds; an empty bounded search is not a rank-zero result.
- `verify_elkies_2026_rank19_rank9_base.sage` is the short promoted replay for
  masks `42110:43109`: nine independent base-Jacobian points, the exact
  degree-two isogeny to the paired base, a birational pointed-quartic map to
  the `t`-line with nine explicit rational `(t,u,v)` points, and generic
  surface rank at least 19.
- `search_elkies_2026_rank9_paired_base.sage` is the first actual search on
  that promoted base. It records rank/saturation evidence, a 100-digit height
  Gram and LLL basis, enumerates by height while selecting separately for
  small `t(P)`, specializes all nineteen visible sections, tests the other
  39,118 covers, certifies displayed rank incrementally, and applies Nagao
  only after the explicit split count.
- `analyze_elkies_2026_high_rank_control_pair_bases.sage` constructs every
  pair in the four exact sets `S(t0)`. It materializes the pair points and
  computes their mod-2 finite-quotient incidence with the public exceptional
  complements, preventing an invalid automatic `published rank + 2` claim.
- `sieve_elkies_2026_rank9_paired_base.sage` exhausts complete canonical-height
  shells on the promoted base. A denominator-aware modular bitset sieve has no
  false negatives, and exact finite quotients certify every retained rank
  lower bound. The height-60 validation certifies all 1,640 fibres; the
  discovery shell `60 < h <= 150` tests another 99,200 parameters compactly.
- `catalogue_elkies_2026_control_pair_bases.sage` builds the 300 pair bases in
  `binom(S(3/8),2)`, including their pointed quartics, minimal Jacobians,
  conductors, root numbers, and all additional exact control incidences.
- `search_elkies_2026_control_pair_base_points.sage` LLL-reduces every
  certified positive-rank control pair, maps its bounded coefficient box
  through the exact degree-two isogeny and pointed-quartic inverse, performs
  the complete split sieve, materializes both signs, certifies the specialized
  ranks, and only then attaches Nagao tie-break scores. The certified
  radius-one box has 6,676 parameters; the split-only radius-two box has
  75,504.

The canonical theorem and full proof boundary are in
[`../BISECTION_PAIR_COVER_GEOMETRY_2026-08-31.md`](../BISECTION_PAIR_COVER_GEOMETRY_2026-08-31.md).

### Generic neighbour and compiler infrastructure

<!-- status-consumer: EC-K3-UNIVERSAL-DEGREE2-FIBRATION-COMPILER fd4b5d71c9497eaf -->

- `exact_neighbor_examples.sage`
- `run_exact_neighbor_engine.sage`
- `verify_exact_neighbor_engine.sage`
- `elliptic_neighbor_compiler.sage`
- `verify_elliptic_neighbor_compiler.sage`
- `elliptic_neighbor_compiler_field_generic.sage`
- `verify_elliptic_neighbor_compiler_field_generic.sage`

The core compiler includes the normalized vertical-padding calculation, the
exact Brandhorst--Elkies coefficient budget for a fully marked old-degree-two
pencil, and both the nonzero-trace chord and trace-zero `(1,x)` quartic
branches.  The resolved local quotient data remain explicit inputs; no
Kodaira/ADE label is silently promoted to a local module.  These are reusable
infrastructure rather than a claim that one particular route is best.

### H3 source-family recovery

The source is the H3 `E7+E8/MW2` family on the level-474 `H21 cap H92` component.
The following scripts are the current source-side proof chain:

- `classify_kumar_e7e8_anchors.sage`
- `deconstruct_x0679_quotients.sage`
- `verify_h21_h92_level474_branch.sage`
- `factor_h21_h92_level474_modp.sage`
- `reconstruct_h21_h92_level474_qq.sage`
- `export_h3_level474_source_family.sage`
- `prove_h3_level474_rational_points.m`
- `verify_h3_noncm_q6_source_anchor.sage`
- `verify_h21_q6_section_descent.sage`
- `lift_h21_p1_modular.sage`
- `verify_h92_section_descent.sage`

### Exact H3 q6 and q8 equation route

The characteristic-zero equation route is exact through `D13/MW4`:

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4.
```

Authoritative replays:

- `derive_h92_q6_child_q8_marking.sage`
- `derive_h92_q6_child_q8_physical_root_target.sage`
- `derive_h92_q6_child_q8_corrected2cover_qq.sage`
- `certify_h92_q6_child_jacobian.sage`

The corrected q8 route deliberately supersedes the degree-46, `true1600`, and
`corrected1278` experiments in `archive/`. See
[`../H3_Q8_REAUDIT_2026-08-22.md`](../H3_Q8_REAUDIT_2026-08-22.md).

### H3 selected degree-two corridor

The following scripts certify one complete lattice/chamber corridor from H3 to the
rootless rank-17 frame:

- `analyze_h3_first_q6_chamber.sage`
- `classify_h3_q6_child_q8_orbits.sage`
- `analyze_h3_d13_q4_chamber.sage`
- `analyze_h3_rank_growing_degree2_chain.sage`
- `analyze_h3_mw10_to_rootless_chambers.sage`
- `verify_h3_d13_to_mw17_path.sage`
- `verify_rank17_to_h3_reverse_transport.sage`

The corridor is

```text
H3 E7+E8/MW2
 --q6 --> E8+E6/MW3
 --q8 --> D13/MW4
 --q24--> D12/MW5
 --q6 --> A11/MW6
 --q8 --> 2A5/MW7
 --q4 --> 3A3/MW8
 --q4 --> A3+2A2/MW10
 --q4 --> 5A1/MW12
 --q4 --> 4A1/MW13
 --q4 --> 3A1/MW14
 --q4 --> 2A1/MW15
 --q4 --> A1/MW16
 --q6 --> rootless/MW17.
```

This is a **selected certified corridor**, not a shortest-path or equation-cost
optimality theorem. At D13, all stated proper presentations through `q=23` were checked
and the first rank growth occurs at `q=24`; three q24 orbits lead to `D12/MW5`, and
orbit 85 was selected. The later continuation follows deterministic first hits from
that one root-adapted frame. Lateral moves, larger-q exits, the other q24 children, and
alternative multi-step corridors can still be easier at equation level.

The last checker supplies the previously missing positive-frame isometry from
the corridor's rootless endpoint to pinned `rank17_gram.txt`, inverts the full
H3-to-R17 transport, and exports every stage basis in both H3 and pinned-R17
coordinates. It also retains the exact bridge between the dominant D13
lattice marking and the distinct component-nef D13 equation marking; that
bridge changes the embedded `U`.

The former q4/orbit230--q6/orbit1315 and q6/orbit1307 promotions are withdrawn
as equation-cost targets.  Exact physical-chamber audits show respectively a
pseudo-zero after the q4 return and a non-section component-10 zero after the
q6 Weyl repair.  The relevant route and correction scripts are:

- `search_h92_a5a5_zero_changing_loops.sage` — **ACTIVE_SEARCH** for the
  exhaustive 283-first-edge/558-zero q4/q6 loop ranking;
- `certify_h92_a5a5_q6o1307_promoted_route.sage` — **ACTIVE_PROOF** for the
  earlier abstract q6/orbit1307 splice, retained as an exact historical
  comparator only;
- `audit_h92_a5a5_q6o1307_physical_nef.sage` — exact physical-I6 Weyl repair;
  it preserves P1307, rejects component 10 as a zero, and exposes physical
  degree-one components 3, 5, and 9 with estimated RR profile `9 -> 3 -> 2`;
- `export_h92_a5a5_physical_source_marking.sage` — exports the exact
  component-9-zero `2A5` equation frame with all suffix/reverse target fibres;
- `export_h92_a5a5_zero_loop_returned_marking.sage` — exports any exact
  changed-zero `2A5/MW7` marking with its equation-A11 transport;
- `search_h92_d13_zero_changing_d12_presentations.sage --mode a5` — evaluates
  second zero changes using inherited-explicit curve degrees and per-edge
  horizontal-cost floors;
- `certify_h92_a5a5_q6o3372_q6o2052_promoted_route.sage --variant q230` —
  **ACTIVE_PROOF** for q4/orbit230, q4 return, q6/orbit1315, q4 return, q4
  exit, the exact current-`3A3` landing, and pinned R17 endpoint;
- `build_h92_a11_route_optimization_handoff.py` — the compact machine handoff
  generator consumed by
  [`../A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md`](../A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md).

The historical 4,199 and 10,334 totals are not equation-realizable as stated,
and the old q104/13,518 comparator is also non-nef in the physical I6 chamber.
`certify_h92_a5a5_direct_physical_q10.sage` reduces it exactly to q10 and
certifies the canonical current-3A3 landing.
`certify_h92_a5a5_direct_physical_q10_promoted_route.sage` composes that landing
to pinned R17 and promotes the 4,471 operational target (5,071 with the older
affine-omitting curve convention).  This is the current lifting target while
lower-score lateral searches continue.

The promoted replay target now begins one stage earlier:

- `search_h92_d13_zero_changing_d12_presentations.sage --mode d13` — exact
  search over the compact q4/q6/q8 equation-D13 frontier, all degree-one old
  D13 component zeros, exact return fibres, and the fixed current-D12 exit;
- `certify_h92_d13_q4o11_promoted_route.sage` — **ACTIVE_PROOF** for the
  q4/orbit11, q4-return, q24-exit splice, exact current-D12 landing, and full
  pinned-R17 endpoint;
- `search_h92_d13_zero_changing_d12_presentations.sage --mode a11` — the
  companion first-q8 audit.  It scans all 1,119 declared-nef candidates and
  ranks inherited-explicit costs; no presentation beats orbit12. Twenty-four
  nonprimitive-root re-zeroings are handled by saturated unimodular frames
  retaining the embedded simple-root lattice.
- `export_h92_first_q8_source_marking.sage` and
  `export_h92_first_q8_zero_loop_returned_marking.sage` — export the exact
  equation-explicit E8+E6 source and its q4/orbit11 changed-zero return state;
- `certify_h92_first_q8_q4o11_promoted_route.sage` — **ACTIVE_PROOF** for the
  q4/orbit11, q4-return, q4-exit replacement of the first q8, its exact
  equation-D13 landing, and the full pinned-R17 endpoint.

The D13 splice scores 25,323 against the measured-RR-calibrated direct q24
score 27,885. It rejoins the exact current D12 basis, keeps A11 q8/orbit12,
then uses the promoted q4/orbit230 and q6/orbit1315 double-zero suffix.

The earlier first q8 now has a cheaper exact replay: q4/orbit11, q4 return,
q4 exit scores 3,961 against 5,802 for direct q8. A complete second q4/q6/q8
zero-loop layer from the returned E8+E6 marking has 38 exact presentations and
no winner. Widened q10 degree-two, q6/q9/q12 degree-three, and
q8/q12/q16 degree-four searches likewise have no winner. The promoted full H3 route is recorded in
`elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json`.
The exact changed-zero D13 landing is exported separately and its 73-presentation
D12 crossover scan rejects carrying that zero into the canonical D13 splice.

### H3 q24 closeout and orbit42 frontier

The direct q24 equation edge is closed.  For a new replay, prefer the promoted
D13 q4/orbit11 zero-changing presentation above. Do not reopen
q32, native q24 suffixes, generic q24 section Hensel lifting, or easy orbit42
zero-pole searches as the active path.

Pinned q24 closeout scripts:

- `certify_h92_q24_equation_d13_to_pinned_r17.sage` — **ACTIVE_PROOF** for the
  current-equation D13 marking, q24/orbit85 D12 marking, pinned-R17 suffix, and
  exact backward-manifest transition/profile lock. This is still a lattice/NS
  boundary, not an equation-level child.
- `replay_h92_q24_d12_discovery_and_theta.sage` — **ACTIVE_PROOF** for the exact
  D12 fibre structure `D24 = Theta + sum m_i C_i`; use this Theta/component
  construction to build the q24 RR pencil.
- `probe_h92_q24_d12_component_valuation_rr_modp.sage` — **REGRESSION** for the
  modular q24/orbit85 component-valuation RR replay. It records the 56-dimensional
  geometric ambient, smooth-collision rank 48, component-cover resolved rank 6,
  `h0=2`, quartic degree 4, and D12 fibre data.
- `extract_h92_q24_d12_modp_signature.sage` — **REGRESSION** for the modular
  q24/orbit85 resolved-RR signature: ambient 56, collision rank 48, resolved rank
  6, kernel 2, quartic degree 4, D12 fibre data. New outputs retain the embedded
  component compiler status and remain `CANDIDATE_*` unless an actual
  equation-level terminal PASS exists; do not treat the modular signature alone
  as a D12 equation certificate.
- `run_h92_q24_orbit85_modp_stack.py` — **ACTIVE_PROOF_AID** to replay the
  fixed q24/orbit85 modular construction stack for canonical split-centre
  `C10` primes. It runs the preflight, affine component graph, effective D13
  transport, component-valuation RR compiler, and compact D12 signature
  extractor in order; it is for CRT/rational-reconstruction data collection,
  not Hensel lifting or neighbor search.
- `reconstruct_h92_q24_orbit85_rr_kernel_crt.sage` — **ACTIVE_PROOF_AID** to
  CRT/rational-reconstruct the canonical q24/orbit85 `2 x 56` resolved-RR
  kernel, binary quartic, and Jacobian data from modular signatures with a
  held-out prime check. A passing output is a QQ kernel candidate, not yet the
  full exact D12 equation certificate.
- `reconstruct_h92_q24_orbit85_compact_rows_lll.sage` — **ACTIVE_PROOF_AID** to
  recover the same q24/orbit85 kernel projectively in compact intrinsic
  `(A, BZ, C)` coordinates using CRT plus LLL and a held-out prime. This is the
  fallback when scalar RREF rational reconstruction is denominator-heavy.
- `extend_h92_q24_orbit85_crt_precision.py` — **ACTIVE_PROOF_AID** to add fresh
  q24/orbit85 modular signatures beyond the archived small-prime pool. It
  builds the q24 direct bridge, I9* resolution, canonical `C10` stack, and then
  retries scalar CRT plus compact-row LLL after each accepted prime.
- `build_h92_q24_orbit85_exact_construction_manifest.sage` — **ACTIVE_PROOF** for
  the q24/orbit85 exact-construction contract. It binds the exact Theta divisor,
  component thresholds, modular 2D RR kernel, quartic/Jacobian regression, and
  D12 gates into
  `artifacts/local/elkies-k3/q24-orbit85-exact-construction-manifest.json`.

Orbit42/D12 profile rule:

- `success-path/` — **ACTIVE_ROUTE_LEDGER** containing the version-locked
  launcher/index for every successful equation-lift stage, the exact artifact
  statuses and hashes, the complete pending route to pinned R17, and the
  separately labelled shortcut audits.  Canonical implementations stay in
  this directory and are not copied into the ledger directory.
- Do not use parity/minimum-`P.O` shortcuts for D12 correction classes. The
  R17-directed orbit42 class is the spinor class with correction `3`,
  `P.O=3`, denominator degree `3`, and no extra fibre twist. Live q24/D12
  profile helpers should use the exact D12 discriminant-class lookup table.
- `extract_h92_q24_orbit42_current_equation_bridge.sage` — **ACTIVE_PROOF** for
  the corrected orbit42 divisor in the current-equation D12 frame:
  `mw=(-1,0,-1,-1,0)`, height 7, correction 3, `P.O=3`, fibre twist 0.
- `recover_h92_q24_orbit42_current_equation_section_modp.sage`,
  `recover_h92_q24_d12_orbit42_section_modp.sage`,
  `archive/recover_h92_q24_pointed_zero_pole_sections.sage`, and
  `archive/recover_h92_q24_r17_a11_zero_pole_sections.sage` —
  **HISTORICAL_DIAGNOSTIC**.
  They prove the old easy zero-pole route does not generate the selected
  orbit42 target.
- `audit_h92_q24_orbit42_explicit_multisections.sage` — **REGRESSION** for the
  NS/Picard-level multisection span audit. The target is in the span, but the
  artifact does not compute the fibrewise Abel-Jacobi functions or the A11
  equation.
- `compile_h92_q24_orbit42_a11_chord_modp.sage` — **REGRESSION** only when fed an
  actual P42 section artifact. The preferred exact route is now the direct
  resolved-RR compiler for the degree-two divisor `D42`, not a zero-pole P42
  search.
- `preflight_h92_q24_orbit42_component_valuation_qq.sage` and
  `map_h92_q24_orbit42_i8star_physical_components_qq.sage` —
  **ACTIVE_EXACT_PREREQUISITES** for the resolved-RR compiler.  They pin the
  corrected divisor profile and the two exact physical spinor orientations.
- `recover_h92_q24_orbit42_zero_pole_smallprime.sage`,
  `scan_h92_q24_orbit42_zero_pole_model_modp.sage`, and
  `lift_h92_q24_orbit42_zero_pole_sections_qq.sage` — **EXACT_BOUNDARY_AUDIT**.
  Together they reconstruct eighteen exact rational identity-class zero-pole
  sections.  The modular scan is only a seed-selection aid; the final QQ
  identities are exact.
- `lift_h92_q24_orbit42_spinor_zero_pole_sections_qq.sage` —
  **EXACT_CONSTRUCTION_AID** for the remaining opposite spinor-class pair.
  It solves the cubic-`x`, constant-`y` cancellation branch directly over QQ
  and pins its mod-100003 residue.  Together with the preceding audit this
  represents the full twenty-point zero-pole shell exactly, but does not
  construct the resolved RR kernel or the A11 child.
- `construct_h92_q24_orbit42_exact_section_candidates_qq.sage` —
  **EXACT_POINT_AID_WITH_MODULAR_MARKING_BOUNDARY**.  Exact QQ(u) group law
  combines the 18+2 shell into four projective `(9,9,3)` points with exact
  Weierstrass identities.  Their four-way orbit42 identification uses the
  pinned mod-100003 shell isometry; resolved RR must still select the physical
  orientation.
- `lift_h92_q24_orbit42_resolved_rr_qq.sage` — **ACTIVE_PROOF**. It certifies
  the exact weighted resolved-RR two-plane, quartic, minimized Jacobian and
  A11 fibre classification.
- `certify_h92_q24_orbit42_a11_equation_marking.sage` —
  **ACTIVE_PROOF_WITH_GOOD_REDUCTION_MARKING_BOUNDARY**. The exact
  identity-shell degree fingerprint at p=100003 selects orbit64/mapping7 in
  the C10 orientation, with orbit65/mapping6 as spinor conjugate.
- `lift_h92_q24_a11_q8_residual_resolved_hensel.sage` and
  `derive_h92_q24_a11_q8_difference_qq.sage` — **EXACT_CONSTRUCTION**.  The
  first lifts the regular six-variable component-3 residual chart; the second
  uses fraction-free group law to derive the selected `(16,24,6)` q8
  horizontal without the former 36-variable lift.
- `lift_h92_q24_a11_q8_resolved_rr_qq.sage` — **ACTIVE_PROOF**.  It certifies
  the exact 14-to-2 recurrence RR plane, quartic, globally minimal
  `2I6+12I1` Jacobian, Euler number 24, and `2A5/MW7` classification without
  a Groebner basis.
- `certify_h92_q24_a11_q8_equation_marking_qq.sage` — **ACTIVE_PROOF**.  The
  nodal-cubic sign identity selects `old_A11_component_9` as the exact child
  zero, transports the affine component as a `(4,6,0)` section, attaches the
  two physical A5 chains, and verifies determinant-`-1` NS transports in both
  directions.
- `certify_h92_q4o208_physical_q4o1584_rr_qq.sage` and
  `certify_h92_q4o1584_equation_marking_qq.sage` — **ACTIVE_PROOF**. They
  certify the exact q4/orbit1584 RR plane and `D4+A3+3A1/MW7` Jacobian, then
  point the second-I6-affine zero and the old-`C0` branch required next.
- `certify_h92_q4o1584_physical_q4o164_rr_qq.sage` and
  `certify_h92_q4o164_c8_equation_marking_qq.sage` — **ACTIVE_PROOF**. They
  recover the q4/orbit164 pencil by rational branch-value interpolation,
  certify the `2A3+2A1/MW9` Jacobian, and identify old `C8` as its exact zero.
- `certify_h92_q4o323_horizontal_marking_qq.sage`,
  `lift_h92_q4o323_horizontal_via_summands_qq.sage`, and
  `compile_h92_q4o208_q4o1599_a3_2a2_qq.sage` — **CORRECTION AND EXACT
  FIXED-ADE REPLACEMENT**. The complete eleven-section graph excludes the
  former branch-33 q4/orbit323 claim, and exact summand reconstruction
  reproduces that same excluded branch. Exhaustive q4 enumeration instead
  finds lifted branch 16 at q4/orbit1599. Its direct `5 -> 3 -> 2` RR
  calculation gives `I4+2I3+14I1`, hence `A3+2A2/MW10`, without a Groebner
  basis. Branch 16 has sign `+1` and the same NS class in both complete graph
  solutions; after the certified raw-to-physical basis change its divisor is
  exactly `D=O+P+V` with zero MW tail in `V`. Thus its pinned-suffix marking is
  unambiguous for the q4/orbit1599 lattice edge; the remaining twofold choice
  concerns only an unused branch. Its child is not the stored canonical
  q4/orbit323 frame feeding q4/orbit207.
  <!-- status-consumer: EC-K3-H3-Q4O208-Q4O1599-QQ-A3-2A2 2018f08fd6b8e2a9 -->
- `construct_h92_q4o323_horizontal_by_halving_qq.sage` and
  `compile_h92_q4o208_q4o323_a3_2a2_qq.sage` — **EXACT CORRECTED
  FIXED-CORRIDOR EDGE**.  The first verifies the marked doubling relation
  `2*T=P8+2*P18+P33-2*C7`, resolves the inherited global inversion by
  `P.O=1`, and extracts the unique rational half from a linear factor of the
  duplication quartic.  The second identifies the vertical residual as
  components 2 and 3 of the second old `I4`, certifies the resolved
  `5 -> 3 -> 2` RR plane, and produces a minimal `I4+2I3+14I1` Jacobian,
  hence `A3+2A2/MW10`.  Neither script uses a Groebner basis or a new shell
  search.  Literal reuse of the stored q4/orbit207 class is non-nef, but
  `build_h92_q4o323_reflected_fixed_suffix_marking.sage` now carries the q323
  wall reflection through the whole suffix.  Nine further vertical-component
  reflections give its exact physical continuation: a q12 degree-two nef
  edge to `5A1/MW12`, with the next marked target still degree two and
  unimodular transport.  Its cheapest horizontal preflight is
  `D=O+P-4F`, `P.O=10`, height `65/3`, no vertical residual, and old
  q4/orbit208 degree 16.  The component-2 pointing is now exact; recovery of
  this horizontal remains open.
- `certify_h92_q4o323_component2_pointing_qq.sage` — **EXACT POINTING
  PREREQUISITE**.  A split-`I4` toric arc fixes the stored normalization as
  `W=+L0(u)` for `old_A11_component_2`, checks the global `81/729` pointed
  invariant identities, and gives the opposite branch as an exact rational
  section of degrees `x=(6,2)`, `y=(8,3)`.  This uses local series and
  univariate arithmetic only.
- `construct_h92_q4o323_p0_shell_modp.sage` — **MODULAR Q12 COMPILER
  FRONTIER**.  It applies the successful q12/orbit5867 five-functional
  polynomial-square method to the component-2-pointed q323 child.  The script
  has an explicit fibre-specialization gate: `p=31` is rejected because the
  two `I3` valuations become `(4,4)`, while `p=61` preserves `(3,3)` and
  produces the complete 602-section signed `P.O=0` shell.  Of these, 120 have
  full ordinary coefficient-Jacobian rank 12 and are regular Hensel
  candidates.  Resolved component signs, lattice naming, QQ lifts, and the
  physical q12 horizontal are not yet claimed.  No Groebner basis or
  elimination is used.
  <!-- status-consumer: EC-K3-H3-Q4O208-Q4O323-QQ-A3-2A2 a903147a9023d49f -->
- `lift_h92_q4o323_q207_deflated_qq.sage`,
  `compile_h92_q4o323_q207_smooth_rr_qq.sage`, and
  `scan_h92_q4o323_q207_candidates_by_smooth_rr_mod61.sage` — **NEGATIVE
  Q207 SEED AUDIT**.  Two deflations regularize modular candidate 5887 with
  rank sequence `69 -> 139 -> 280`; Newton lifting reconstructs an exact
  `QQ(u)` section at 1024 base-61 digits with maximum rational size 1604
  bits.  Its direct `22 -> 2` smooth RR compiler has rank 20 and `h0=2`, but
  the exact child is `6A1`, not the prescribed `5A1`, and its marked Abel
  trace also differs from q207.  Candidate 5903 fails the marked trace before
  lifting.  The modular RR scan compiles all 948 target-shape candidates in
  the stored shell without elimination; none has root rank five at `p=61`,
  and all eighteen candidates of minimum root rank six fail the independent
  marked trace.  These are construction rejections, not a non-existence
  theorem for the q207 horizontal.  Replay the exact negative lift and RR
  certificate with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/lift_h92_q4o323_q207_deflated_qq.sage \
    --candidate 5887 --maximum-precision 1024 --reconstruction-start 256
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/compile_h92_q4o323_q207_smooth_rr_qq.sage
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/scan_h92_q4o323_q207_candidates_by_smooth_rr_mod61.sage
  ```
- `recover_h92_q4o164_odd_i4_sections_modp.sage` and
  `lift_h92_q4o164_odd_i4_sections_qq.sage` — **EXACT/MODULAR FRONTIER
  DIAGNOSTIC**. They enumerate odd finite-`I4` polynomial sections without a
  Groebner basis, Newton-lift the regular branches, and verify every retained
  `QQ(t)` section by literal substitution. The optional bounded signed-sum
  audit replays with
  `--precision 180 --rational-sum-scan --output artifacts/generated-results/elkies-k3-h3-q4o164-odd-i4-rational-sum-scan-p41.json`;
  it filters 364 words at three fibres, interpolates 22 survivors in the
  marked q8 degree window, and finds no exact rational sum. This is a bounded
  construction audit, not a non-existence theorem for the missing ninth
  direction.
- `recover_h92_q4o164_zero_node_sections_modp.sage` and
  `lift_h92_q4o164_zero_node_sections_qq.sage` — **EXACT/MODULAR FRONTIER
  DIAGNOSTIC**. The exhaustive degree-four zero-node scan at `p=23` has one
  inverse pair; both lift over `QQ`, but specialization quotient tests keep
  them in the known rank-eight shell.
- `probe_h92_q4o164_inherited_p1_abel_trace_modp.sage` — **MODULAR q8
  CONSTRUCTION FRONTIER**. It transports the exact inherited `P1` curve
  through the three certified q4 models using their degree-one pointed
  quartic maps, obtaining degrees `3 -> 6 -> 7`. Fibrewise Abel reduction uses
  only the unique `7 x 8` kernel in `L(8O)`. At `p=131`, `--interpolate`
  reconstructs the full trace section with degrees `x=(32,28)` and
  `y=(48,42)` from 122 good fibres, with 61/31 holdout fibres and an exact
  modular Weierstrass identity. Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/probe_h92_q4o164_inherited_p1_abel_trace_modp.sage \
    --interpolate \
    --output artifacts/local/elkies-k3/q4o164-inherited-p1-abel-trace-section-mod131.json
  ```

  The trace is not the q8/orbit376 horizontal modulo the already named
  degree-one sections: their marked tails leave the exact correction
  `(0,0,0,0,-1,-1,2,0,0)`. This remains a modular construction, not an
  equation certificate.
- `audit_h92_q4o164_integral_basis_height_gram.sage` — **EXACT HEIGHT/MARKING
  AUDIT**. Raw node incidence had mislabeled the infinity-`I4` component of
  `B7`. Multiplying every basis section and pair sum by four clears all
  component groups; compact pole growth then gives the corrected height Gram
  with determinant `459/8`, rather than relying on unresolved node labels.
  A finite positive-definite enumeration finds 16 integral embeddings in the
  C8-pointed marked MW9 lattice, eight compatible with the valid
  `B0,...,B6` profiles. All eight contain the q8 residual but with different
  basis words. Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/audit_h92_q4o164_integral_basis_height_gram.sage
  ```

  The old coarse profile `[1,0,0,2]` for `B7` and the resulting height `3`
  for `N=2*B0+B5+B7` are withdrawn; exact pole growth gives `h(N)=13/4` and
  an odd infinity-`I4` label.
- `identify_h92_q4o164_q8_horizontal_mod131.sage` — **EXACT MODULAR q8
  HORIZONTAL**. It combines each of those eight embeddings with the complete
  inherited-`P1` Abel trace and the exact saturated relation
  `3*C8opp=-2*B0-3*B1-4*B2+3*B3-2*B4+2*B5-B6-2*B7`. Exactly one candidate
  has the certified q8 pole profile `x=(12,8)`, `y=(18,12)`, `P.O=4`, giving

  ```text
  H=T-C8opp-B0+2*B1+B2-3*B3-B4-2*B5+B7.
  ```

  Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/identify_h92_q4o164_q8_horizontal_mod131.sage
  ```

  Pass `--prime p` to replay any exact good-prime Abel-trace artifact.  For
  larger CRT primes, the trace probe accepts `--good-fibre-limit 100`; this
  retains 91 interpolation fibres and nine independent holdouts without
  scanning all of `GF(p)`.
- `reconstruct_h92_q4o164_q8_horizontal_crt_qq.sage` — **EXACT QQ q8
  HORIZONTAL**. It combines the monic-normalized compact coordinates of `H`
  over 22 good primes. Independent coefficient reconstruction exposes the
  large common scales; simultaneous projective LLL on the 22- and 32-entry
  coordinate vectors recovers primitive vectors of at most 363 and 526 bits
  from the 566-bit CRT modulus. It accepts the result only after exact
  substitution in the compact Weierstrass equation and reduction back to
  every input prime. This targets the much smaller `(12,8)/(18,12)` q8
  section directly rather than first reconstructing the `(32,28)/(47,42)`
  Abel trace. Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/reconstruct_h92_q4o164_q8_horizontal_crt_qq.sage
  ```

  Terminal status is `PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT`. The exact
  fourfold pole degree is 172, giving canonical height 11.  The horizontal is
  not yet the q8 child equation.
- `compile_h92_q4o164_q8o376_smooth_rr_qq.sage` — **EXACT QQ RESOLVED RR,
  4A1 JACOBIAN, AND P1229 POINTING**. It applies the earlier A11 q8 coefficient recurrence to the exact
  <!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-SMOOTH-RR c0bc752848961743 -->
  <!-- status-consumer: EC-K3-H3-Q12O5867-DEGREE1-SECTION-QQ 28056772646e9fc7 -->
  <!-- status-consumer: EC-K3-H3-Q12O5867-QQ-ROOTLESS 6a3dbee5942ddd0a -->
  inherited-`P1` horizontal. With `deg Z=4`, the saturated chord ambient is
  `a=AA/Z^2`, `deg AA<=8`, and `b=BB/Z`, `deg BB<=2`: dimension 12. The exact
  congruence `AA*X=BB*Y mod Z^2` has rank 8 and leaves a four-plane, verified
  independently by a full coefficient matrix and at all 22 CRT primes. The
  two full old-`I4` nonidentity `A3` chains give two independent toric quotient
  rows, so the resolved pencil has dimension two. Its quartic has Jacobian
  fibres `4I2+16I1`, Euler number 24, and root data `(4,8,16)`, hence
  `4A1/MW13` in the pinned rank-19 lattice. Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/compile_h92_q4o164_q8o376_smooth_rr_qq.sage
  ```

  The selected marked embedding and the exact `B0` tangent orient the finite
  `T=0` I4 chain as old components `6`, missing, `5`; thus the first pointed
  origin is old component 6. More economically, `P1229` is the nonidentity
  component of the finite I2 at `T=25281/168246841`. A rational arc on its
  exceptional conic selects quartic sign `-1`, and direct pointing verifies
  `81*A_pointed=A_child` and `729*B_pointed=B_child`. Terminal status is
  `PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO`. The downstream
  q12/orbit5867 equation is now complete. No Groebner basis is used.
- `construct_h92_q12o5867_degree1_branch_mod131.sage` and
  `lift_h92_q12o5867_degree1_branch_qq.sage` — **FIRST PROMOTED Q12 COMPILER
  SECTION**. The unique parent-degree-one branch has exact marked word
  `2*C8opp+B1+2*B2-2*B3+B4-2*B5+B6+B7`. Modulo 131, restricting that curve
  through the resolved q8 pencil and applying the P1229-pointed quartic map
  gives a polynomial child section of degrees `(4,6)` with component profile
  `(0,0,0,0)`. Its 13 coefficient equations have rank 12 in 12 variables,
  with selected minor determinant 50. Ordinary Hensel lifting and rational
  reconstruction recover the exact QQ section with maximum coefficient height
  800 bits. Replay with:

  ```bash
  sage -python elkies-k3/scripts/construct_h92_q12o5867_degree1_branch_mod131.sage
  sage -python elkies-k3/scripts/lift_h92_q12o5867_degree1_branch_qq.sage
  ```

  Terminal statuses are
  `PASS_MOD131_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD` and
  `PASS_EXACT_QQ_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD`. This remains a
  historical exact branch: the final equation compiler uses the replacement
  word below. No elimination or Groebner basis is used.

- `classify_h92_q12o5867_p0_shell_against_lattice_mod89.sage`,
  `lift_h92_q12o5867_replacement_word_seeds_qq.sage`,
  `construct_h92_q12o5867_target_horizontal_qq.sage`, and
  `compile_h92_q12o5867_smooth_rr_qq.sage` — **EXACT Q12/O5867 ROOTLESS
  EQUATION**. Complete shells at primes 83, 89, 137, and 151 do not realize
  the nominal Q1 branch. The equation-effective word is
  `499+500+69+511-489+933-913`, with actual parent degrees `(4,2,1,5)` and
  parent `a-b=(4,2,1,4)`. Exact Hensel lifts and group law construct its
  horizontal. The smooth chord calculation has ambient dimension 22, rank
  20, and `h0=2`; the resulting minimal Jacobian has degrees `(8,12,24)` and
  geometrically `24I1`. Replay with:

  ```bash
  sage -python elkies-k3/scripts/classify_h92_q12o5867_p0_shell_against_lattice_mod89.sage
  sage -python elkies-k3/scripts/lift_h92_q12o5867_replacement_word_seeds_qq.sage
  sage -python elkies-k3/scripts/construct_h92_q12o5867_target_horizontal_qq.sage
  sage -python elkies-k3/scripts/compile_h92_q12o5867_smooth_rr_qq.sage
  ```

  Terminal status is
  `PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN`. No elimination or
  Groebner basis is used.
- `construct_h92_q12o5867_p0_shell_modp.sage`,
  `identify_h92_q12o5867_two_primary_cosets_mod89.sage`,
  `attempt_h92_q12o5867_three_target_halvings_qq.sage`,
  `audit_h92_q12o5867_target_support_cross_prime.sage`, and
  `verify_h92_q12o5867_three_halvings_and_basis_qq.sage` — **CLOSED EXACT
  TWO-PRIMARY BOUNDARY**.  The shell constructor's
  `--all-component-profiles` mode exhausts all 16 equation incidence modes at
  `p=83,89,137`.  The historical current subgroup has Smith factors
  `1^10,2,2,2`.  Three minimum-`L1` independent targets have exact QQ
  polynomial halves (classes 146, 30, 22), literal doubling identities,
  profiles `0101,1101,1110`, and heights `3,5/2,5/2`; abstractly they reduce
  the index from eight to one.  Exact full-NS names are obtained from smooth
  intersection fingerprints, not from a global permutation of equation-mode
  bits.  The bounded equation-level pool nevertheless has rank 13 and index
  two, with height determinant 237 instead of `237/4`.  The declared boundary
  therefore supplies no saturated thirteen-section equation basis and no
  words for the 42 controls; q12/orbit5867 is closed as the rank-32
  point-production route while its exact arbitrary-point map is retained.
  The final certificate is
  `artifacts/generated-results/elkies-k3-q12o5867-two-primary-boundary.json`.
- `construct_h92_q12o5867_rootless_p0_shell_mod131.sage`,
  `select_h92_q12o5867_rootless_mod131_basis.py`,
  `lift_h92_q12o5867_rootless_selected_basis_qq.sage`, and
  `certify_h92_q12o5867_rootless_height_basis_qq.sage` — **EXACT ROOTLESS
  RANK-17 SECTION BASIS**. The exhaustive p=131 polynomial shell contains
  2,622 signed sections, all regular for Hensel lifting. A 17-section modular
  basis has height determinant 948 and determinant-minus-one transport to
  pinned R17. All seventeen reconstruct over QQ at 4,096 p-adic digits, with
  maximum coefficient sizes 7,895--7,933 bits. Exact QQ intersections recover
  the same Gram and prove rank at least 17; rootlessness proves torsion is
  trivial. Replay with:

  ```bash
  sage -python elkies-k3/scripts/construct_h92_q12o5867_rootless_p0_shell_mod131.sage
  python3 elkies-k3/scripts/select_h92_q12o5867_rootless_mod131_basis.py
  sage -python elkies-k3/scripts/lift_h92_q12o5867_rootless_selected_basis_qq.sage
  sage -python elkies-k3/scripts/certify_h92_q12o5867_rootless_height_basis_qq.sage
  ```

  Terminal status is
  `PASS_EXACT_QQ_Q12O5867_ROOTLESS_RANK17_HEIGHT_BASIS_PINNED`. No
  elimination or Groebner basis is used.

  <!-- status-consumer: EC-K3-H3-Q12O5867-QQ-R17-BASIS a2097150acf00645 -->
- `certify_h92_q12o5867_endpoint_qq.sage` — **EXACT ENDPOINT THEOREM**.
  The old finite `I2` support `v=0` gives an exact rational point on the q12
  binary quartic; its pointed generalized Weierstrass invariants satisfy the
  exact `81/729` identities with the stored endpoint equation. Exact counts
  over `F_p` and `F_{p^2}` at `p=131,137` give rank-20 reductions with
  incompatible Artin--Tate NS discriminant square classes, proving geometric
  Picard rank 19. The only possible proper integral enlargement of the
  determinant-948 section lattice has index two and odd new-vector norm 73,
  so evenness of the rootless MW height lattice excludes it. Replay with:

  ```bash
  sage -python elkies-k3/scripts/certify_h92_q12o5867_endpoint_qq.sage
  ```

  Terminal status is
  `PASS_EXACT_Q12O5867_SOURCE_IDENTITY_RHO19_FULL_MW_R17`. Thus the endpoint
  is exactly the certified H3 source K3 and its full geometric MW group is
  saturated R17 of rank 17, trivial torsion, and determinant 948. No
  Groebner basis or surface elimination is used.

  <!-- status-consumer: EC-K3-H3-Q12O5867-ENDPOINT-QQ a83b08acd921c32b -->
- `build_h92_q4o164_all_zero_q8_known_horizontal_audit.py` — **BOUNDED
  NEGATIVE ROUTING AUDIT**. For each of the eight certified effective q4/o164
  origins, a q8/degree-two shell was generated with
  `search_root_adapted_weyl_neighbors.sage --q 8 --degree 2
  --mw-vector-cap 10000` and filtered with
  `score_h92_marked_frontier_equation_cost.sage
  --require-known-horizontal`. The aggregate covers 119,220 primitive
  candidates and finds zero horizontals in the corresponding explicit-section
  subgroup. Replay the compact aggregate after the individual shell/score
  files are present with:

  ```bash
  python3 elkies-k3/scripts/build_h92_q4o164_all_zero_q8_known_horizontal_audit.py
  ```

  This is exact for the stored capped shells, not a global non-existence
  theorem. Together with the 364-word odd-I4 rational-sum scan, it supports
  the inherited-`P1` Abel-reduction construction target above.
- `certify_h92_marked_route_to_pinned_r17.sage` — **ACTIVE_ROUTE_PROOF**. It
  certifies the q323-free marked q4/q4/q8/q12 suffix and terminal integral
  isometry to pinned R17. It is not an equation certificate for q8/orbit376 or
  q12/orbit5867 (or fallback q12/orbit4484) and does not supply the endpoint
  section/rank package.
- `transport_h92_q24_a11_degree_one_shell_qq.sage` —
  **EXACT_CONSTRUCTION_AID**.  It transports the two degree-one identity
  curves and the degree-one spinor curve by Möbius inversion and exact binary
  quartic covariants, and verifies the pointed opposite relation over QQ.
- `audit_h92_q24_a11_missing_direction_alternatives.sage` —
  **EXACT_NEGATIVE_CONSTRUCTION_AUDIT**.  It proves that all three exact
  degree-one transports remain in the known fifth-coordinate-zero parent
  hyperplane, that q8 orbit2162 only reverses the required child coordinate,
  and that the smallest parent carrier is the D12 `P.O=4` vector
  `(0,0,0,0,1)`.
- `audit_h92_a11_quintic_bridge_zero_mismatch.sage` —
  **EXACT_REJECTION_CERTIFICATE**.  It replays the selected `R3`-zero frame
  and rejects the former quintic shortcut, which had mixed it with `A0`-zero
  coordinates.  The compatible degrees are 46 and 4, and the claimed word
  does not equal `M`.
- `audit_h92_a11_explicit_aj_carriers.sage` —
  **EXACT_ALTERNATIVE-CONSTRUCTION_AUDIT**.  It scans all stored explicit
  (-2)-classes in the selected marking.  No positive single carrier works;
  the first positive subset uses degree-40 and degree-44 traces, so it is a
  fallback rather than the active route.
- `score_h92_a11_equation_cost_neighbors.sage` and
  `certify_h92_a11_equation_cost_orbit849.sage` —
  **EXACT_ROUTE-PLANNING/AUDIT**.  The first scores all declared neighbours;
  the second applies the stronger nefness and marked-U gates.  Low-score
  orbits 849 and 591 fail nefness; passing lateral candidates do not yet have
  the required 2A5-to-pinned-R17 continuation.
- `probe_h92_q24_a11_close_p24_quintic_modp.sage` —
  **MODULAR_RESOLUTION_GATE**.  It composes both pointed quartics using only
  univariate arithmetic and rejects the naive unresolved restriction:
  q24 degree 14 and A11 degrees 39/41 replace the compatible
  strict-transform degree 46.  Both tangent signs vanish in the raw chord
  discriminant but are cancelled from the normalized quartic, locating the
  issue at the resolved base locus.
- `preflight_h92_q24_o12_p42_exact_q6_points.sage` and
  `run_h92_q24_orbit42_fast_parallel.py` — **EXACT_NEGATIVE_AUDIT** for the
  rejected q6-point transport.  The named equation-D13 coordinate conversion
  is exact, but `O12` and `P42` have q6 degrees `435` and `703`; the runner
  stops cleanly before the archived transport/orientation experiments.
- `analyze_h92_q24_orbit42_identity_halving.sage` and
  `recover_h92_q24_orbit42_by_identity_halving_qq.sage` —
  **SHORTCUT_BOUNDARY_AUDIT**.  The first proves the exact identity-shell
  doubling relation.  The second uses the exact QQ points but only a modular
  shell marking/halving census; its four rational degree-three candidates
  have degree-18 squarefree chord branches and no A11 fibre.  This is not a
  characteristic-zero non-existence theorem.
- `recover_h92_q24_exact_by_qq_trace_interpolation.sage` — **HISTORICAL_DIAGNOSTIC**
  unless it is refactored to serve the Theta/component construction directly.

Archived q24 dead ends:

- `archive/lift_h92_q24_direct_hensel.sage`
- `archive/exactify_h92_q24_from_padic_srr.sage`

These record why p-adic direct lifting/SRR was abandoned. Do not increase Hensel
precision or reconstruct the generic q24 section as the next step. The q24
`D12/MW5`, orbit42 `A11/MW6`, q8 orbit12 `2A5/MW7`, physical q4/orbit208
`3A3/MW8`, q4/orbit1584 `D4+A3+3A1/MW7`, and q4/orbit164
`2A3+2A1/MW9` children are closed. Continue from the C8-pointed child to
q8/orbit376.

### Q80 compiler and regression route

The generic Q80 lattice corridor and the exact CM24 characteristic-zero shadow are
separate claims. Current entry points include:

- `verify_q80_to_rootless_path.sage`
- `trace_q80_candidate1_marked_transport.sage`
- `score_h3_d12_q80_crossovers.sage` — exact transport of all eleven retained
  Q80 fibre classes into the current equation-side H3 D12 frame, with the
  negative crossover-cost audit documented in
  [`../H3_D12_Q80_CROSSOVER_AUDIT_2026-08-24.md`](../H3_D12_Q80_CROSSOVER_AUDIT_2026-08-24.md)
- `recover_q80_final_q6_via_basis_sections.sage`
- `certify_q80_final_q6_char0_rr_from_basis.sage`
- `compile_q80_final_q6_char0_child.sage`
- `verify_q80_lowq_cm24_equations.sage`

The CM24 child `4A2+A3+A5/MW2` is a regression/specialization shadow; it is not the
generic rootless `MW17` equation.

### Backtracks and comparison fibrations

These remain useful and reproducible, but they are not the H3 source construction:

- `verify_e6_neighbor_chain.sage`
- `verify_rank17_h8_split.sage`
- `analyze_rank17_h8_q9_fibers.sage`
- `verify_humbert8_d9e7_two_neighbor.sage`
- `verify_picard20_mw1_path.sage`
- `recover_mw1_*_glue.sage`

Use the names **Low-q MW2 Backtrack**, **E6 Backtrack**, **H2 Symmetry Comparison**,
and **H2 Minimal-MW Comparison** from
[`../CONSTRUCTION_ROUTES.md`](../CONSTRUCTION_ROUTES.md).

## Historical root-level launchers retained intentionally

The old E6/MW3 attack was run through scripts such as

- `start_e6_attack.sh`
- `run_e6_mw3_probe.py`
- `run_mw3_local_probe.py`
- `run_iv_mod_p2.py`

They are retained because they preserve exact bounded negative work and are referenced by
historical notes. They are **HISTORICAL_DIAGNOSTIC**, not source-construction entry
points. The guessed split E6 chart was never obtained by transporting the certified
neighbour chain, and its old survivors fail the missing off-diagonal height gate.

## Archive policy

[`archive/README.md`](archive/README.md) classifies the imported snapshots by failed or
superseded programme. Rules:

1. An archived script is never the authoritative proof source.
2. If the same filename exists in the root, use the root copy; the archived copy is an
   older snapshot unless a note explicitly says otherwise.
3. Keep unique failed scripts when they explain why an approach was abandoned.
4. Delete only demonstrated byte-identical copies or accidental empty files.
5. New failed experiments should be archived together with a short outcome note, not
   dropped into an unlabelled filename pile.

The 2026-08-23 audit removed five byte-identical ` (1)` copies and made no broad moves.
Moving historical root scripts without repairing every old command would make the record
less reproducible, not cleaner.

## Where to read the history

- [`../SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](../SCRIPT_ROUTE_AND_FAILURE_LEDGER.md) —
  route chronology, selection reasons, rejected assumptions, and remaining gaps.
- [`../CONSTRUCTION_ROUTES.md`](../CONSTRUCTION_ROUTES.md) — named geometric routes.
- [`../H3_Q8_REAUDIT_2026-08-22.md`](../H3_Q8_REAUDIT_2026-08-22.md) — exact q8 bug
  diagnosis and repair.
- [`../Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](../Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md)
  — Q80 terminal marking/RR closeout.

### q12/o4484 fallback section compiler

q12/orbit4484 remains fully lattice-certified as a fallback. Its horizontal
section has the following exact four-section group-law compilation plan:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_h92_q12o4484_p0_section_word.sage \
  --marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --q12-cost artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap10000-rootless-equation-cost.json \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o4484-four-p0-section-word.json
```

This exhausts the marked norm-at-most-four shell, imposes the physical simple
and affine component gate, rules out words of length at most three modulo the
explicit subgroup, and selects four genuine `P.O=0` chamber sections with
q4/o164 parent degrees `3,2,3,2`. It is an exact MW/group-law compiler target,
but q12/orbit5867 is preferred because its corresponding degrees are
`3,2,1,2` and its parent `a-b` values are `2,2,1,1`. The characteristic-zero
section equations are still to be constructed.

### Rootless-exit low-pole compiler frontier

The q12 shell was expanded from 10,000 to 50,000 MW representatives by:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt \
  --root-rank 4 --q 12 --degree 2 --mw-vector-cap 50000 \
  --adapt-mw-at-least 17 --rank-growth-only \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap50000-mw17-neighbors.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/score_h92_marked_frontier_equation_cost.sage \
  --neighbors artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap50000-mw17-neighbors.json \
  --marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --min-child-mw-rank 17 --retain 100 \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap50000-rootless-equation-cost.json
```

This retained 28 rootless candidates among 61,352 primitive candidates. The
full q12 MW shell has 114,347,416 vectors, so this is a bounded expansion, not
a global optimum theorem. Orbit indices change with the sample cap; subsequent
certificates identify candidates by the exact fibre vector and its SHA-256
fingerprint.

A doubled-cap stability run with the same command and
`--mw-vector-cap 100000` screened 133,331 primitive candidates and found the
identical set of 28 rootless fibres: vectors 50,001 through 100,000 contributed
no new rootless exit. Build the compact comparison after generating the larger
shell with:

```bash
python3 elkies-k3/scripts/build_h92_q12_cap_stability_audit.py \
  --smaller artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap50000-mw17-neighbors.json \
  --larger /tmp/q12d2-cap100000-mw17-neighbors.json \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12-cap50000-cap100000-rootless-stability-audit.json
```

This strengthens bounded stability only; it does not exhaust the full shell.

The low-pole construction is compared across all 28 retained q12 exits and all
four retained q20 exits by:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_h92_rootless_p0_section_word_frontier.sage \
  --marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --cost artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap50000-rootless-equation-cost.json \
  --cost artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q20d2-cap10000-rootless-equation-cost.json \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-rootless-p0-section-word-frontier.json
```

The enlarged exact marked quotient search retains fibre fingerprint
`246597bd...` (sample-local orbit43833) as a low-pole alternative. Its three
physical sections have current `P.O=(0,0,1)`, q4/o164 parent `a-b` values
`(0,1,2)`, and parent degrees `(3,2,3)`. It uses one fewer new branch and halves
the parent pole-proxy sum relative to orbit5867, but introduces a `P.O=1`
branch and four named corrections. A second exact fibre, fingerprint
`5d96a74a...` (sample-local orbit49112), uses the same three new sections and
differs only in its explicit correction. Both complete marked routes are
exactly pinned to R17.

The preferred optional endpoint remains q12/orbit5867. Its nominal lattice
word uses four polynomial `P.O=0` branches with q4/o164 parent degrees
`(3,2,1,2)` and parent `a-b` values `(2,2,1,1)`, outperforming orbit4484 on
both planning totals.
Its stable certificate identity is `q12o5867_after_q8o376` with exact fibre
fingerprint `d676cab5...`; the expanded-frontier index 30357 is sample-local.
Orbit5867 is now equation-explicit over characteristic zero. Its physical
compiler instead selected degrees `(4,2,1,5)` and `a-b=(4,2,1,4)` after the
nominal Q1 branch failed complete good-prime shell tests. The other two q12
choices remain compiler plans.

The retained low-pole alternative replays with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_degree2_candidate.sage \
  --source-marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --source-frame artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt \
  --fibre 6,2,-13,10,2,5,-12,2,9,7,-6,-22,6,2,-10,0,-1,-1,-2 \
  --candidate-label q12cap50k_o43833_after_q8o376 --target pinned_R17 \
  --frame-output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12-fp-6-2-m13-10-rootless-frame.txt \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12-fp-6-2-m13-10-rootless-certificate.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_route_to_pinned_r17.sage \
  --edge artifacts/generated-results/elkies-k3-h3-q4o208-physical-q4o1584-lateral-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12-fp-6-2-m13-10-rootless-certificate.json \
  --source q4o208_equation_physical_3A3_C5_zero \
  --status PASS_EXACT_Q323_FREE_Q4O1584_Q4O164_Q8O376_Q12_FP_M13_10_PINNED_R17_ROUTE \
  --output artifacts/generated-results/elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12-fp-m13-10-pinned-r17-route-certificate.json

python3 elkies-k3/scripts/build_h92_q12_expanded_compiler_pareto.py
```

The optional orbit5867 endpoint and its pinned-route composition replay with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_degree2_candidate.sage \
  --source-marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --source-frame artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt \
  --fibre 6,2,-16,18,2,5,-11,2,13,11,-6,-41,6,5,-14,0,-3,-1,-2 \
  --candidate-label q12o5867_after_q8o376 --target pinned_R17 \
  --frame-output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-frame.txt \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-certificate.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_route_to_pinned_r17.sage \
  --edge artifacts/generated-results/elkies-k3-h3-q4o208-physical-q4o1584-lateral-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-certificate.json \
  --source q4o208_equation_physical_3A3_C5_zero \
  --status PASS_EXACT_Q323_FREE_Q4O1584_Q4O164_Q8O376_Q12O5867_PINNED_R17_ROUTE \
  --output artifacts/generated-results/elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12o5867-pinned-r17-route-certificate.json
```

### Elkies 2026 compact-t positive controls and residual gate

The published compact chart is now the only active direct specialization
coordinate. Replay all seventeen sections and the exact coordinate match, then
certify the four public exceptional fibres with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/verify_elkies_2026_published_r17_target.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/match_h92_q12o5867_to_elkies_2026_qq.sage
python3 elliptic-curves/scripts/verify_elkies_2026_high_rank_calibrations.py
```

The last command imports public exact point sets of lengths 25--28 and
certifies one combined generic-plus-complement independence matrix at each
fibre, with quotient gains `8,9,10,11`.

The final neighbour's arbitrary-point map and its five-control backward
calibration are replayed with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/q12o5867_genus_one_point_factory.sage \
  --mode controls \
  --output artifacts/generated-results/elkies-k3-q12o5867-genus-one-point-factory-controls.json
```

The script exposes both the forward parent-to-published map and its inverse,
and checks all 42 public-complement points by exact round trip. See
[`../Q12O5867_GENUS_ONE_POINT_FACTORY_2026-08-31.md`](../Q12O5867_GENUS_ONE_POINT_FACTORY_2026-08-31.md)
for the formulas and the current MW13-coordinate boundary.

The accepted complete scoring calibration is:

```bash
python3 elkies-k3/scripts/calibrate_elkies_2026_positive_controls_nagao.py
```

It enumerates every primitive compact `t=a/b` through height 10,000, uses
three disjoint 34-prime ensembles, ranks first by weakest standardized block,
and fails unless every control lies in the top one percent. This is a
heuristic ranking gate, not rank evidence.

The next arithmetic gate is the actual residual 2-Selmer quotient:

```bash
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --timeout 300 --overwrite
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --backend eclib --timeout 300 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_eclib_v1.json \
  --overwrite
```

Both pinned attempts time out without a Selmer dimension and therefore forbid
search. Generate the unconditional basis-level Magma path with
`build_elkies_2026_rank28_relative_descent_magma.py`; it applies the dimension
gate before relative cover construction and contains no point search. Magma is
not installed on this host.

`probe_q12o5867_pari_two_cover.py`, `probe_q12o5867_ratpoints.py`,
`probe_q12o5867_section_charts.py`, `probe_q12o5867_mwrank.py`, and
`search_q12o5867_section_slope_slices.py` now require a passing
`--residual-selmer-gate` for the identical parameter and minimal curve. The
old low-complexity x-ansatz raw search is hard parked. BNF-free signatures,
norm-one candidates, incomplete class ledgers, and local candidates cannot
satisfy the gate.

### Fixed-corridor reverse lift from the q12/o5867 endpoint

<!-- status-consumer: EC-K3-H3-FIXED-REVERSE-A1-2A1-QQ 26bce707a77972c4 -->
<!-- status-consumer: EC-K3-H3-FIXED-REVERSE-4A1-QQ 8f3863b630d27e16 -->

The exact q12/o5867 rootless equation and its pinned 17-section basis now give
a cheaper construction direction for the historical fixed corridor.  The
terminal `A1 -> rootless` and `2A1 -> A1` arrows have been executed in reverse
over `QQ`, including their prescribed zeros and root components.  Replay with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/construct_h92_fixed_final_a1_horizontal_from_q12_endpoint_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_h92_fixed_final_a1_reverse_rr_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_fixed_final_a1_reverse_pointing_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/construct_h92_fixed_reverse_2a1_horizontal_from_a1_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_h92_fixed_final_a1_reverse_rr_qq.sage \
  --target artifacts/local/elkies-k3/fixed-reverse-2a1-horizontal-from-a1-qq.json \
  --surface artifacts/local/elkies-k3/fixed-final-a1-reverse-rr-qq.json \
  --output artifacts/local/elkies-k3/fixed-reverse-2a1-rr-qq.json \
  --target-status PASS_EXACT_QQ_FIXED_REVERSE_2A1_HORIZONTAL_ON_A1 \
  --surface-status PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_RR_JACOBIAN \
  --edge 'A1/MW16 reverse to 2A1/MW15' \
  --expected-i2-count 2 --expected-ade 2A1 --expected-mw-rank 15 \
  --result-status PASS_EXACT_QQ_FIXED_REVERSE_2A1_RR_JACOBIAN
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_fixed_reverse_2a1_pointing_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/construct_h92_fixed_reverse_3a1_horizontal_from_2a1_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_h92_fixed_final_a1_reverse_rr_qq.sage \
  --target artifacts/local/elkies-k3/fixed-reverse-3a1-horizontal-from-2a1-qq.json \
  --surface artifacts/local/elkies-k3/fixed-reverse-2a1-rr-qq.json \
  --output artifacts/local/elkies-k3/fixed-reverse-3a1-rr-qq.json \
  --target-status PASS_EXACT_QQ_FIXED_REVERSE_3A1_HORIZONTAL_ON_2A1 \
  --surface-status PASS_EXACT_QQ_FIXED_REVERSE_2A1_RR_JACOBIAN \
  --edge '2A1/MW15 reverse to 3A1/MW14' \
  --expected-i2-count 3 --expected-ade 3A1 --expected-mw-rank 14 \
  --result-status PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_fixed_reverse_2a1_pointing_qq.sage \
  --surface artifacts/local/elkies-k3/fixed-reverse-2a1-rr-qq.json \
  --curves artifacts/local/elkies-k3/fixed-reverse-3a1-horizontal-from-2a1-qq.json \
  --rr artifacts/local/elkies-k3/fixed-reverse-3a1-rr-qq.json \
  --output artifacts/local/elkies-k3/fixed-reverse-3a1-pointing-qq.json \
  --surface-status PASS_EXACT_QQ_FIXED_REVERSE_2A1_RR_JACOBIAN \
  --curves-status PASS_EXACT_QQ_FIXED_REVERSE_3A1_HORIZONTAL_ON_2A1 \
  --rr-status PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN \
  --result-status PASS_EXACT_QQ_FIXED_REVERSE_3A1_POINTING \
  --edge '3A1/MW14 --q4 orbit498--> 2A1/MW15' \
  --old-i2-support "$(jq -r '.effective_horizontal_components[0].child_I2_support' \
    artifacts/local/elkies-k3/fixed-reverse-2a1-pointing-qq.json)" \
  --horizontal-roots-key effective_3A1_horizontal_roots_on_2A1_source \
  --horizontal-root-class-key class_in_2A1_coordinates \
  --expected-child-i2-count 3 --remaining-vertical-root-count 1
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_h92_fixed_reverse_4a1_physical_nef.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_h92_fixed_reverse_4a1_physical_rr_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_fixed_reverse_2a1_pointing_qq.sage \
  --surface artifacts/local/elkies-k3/fixed-reverse-3a1-rr-qq.json \
  --curves artifacts/local/elkies-k3/fixed-reverse-4a1-horizontal-from-3a1-qq.json \
  --rr artifacts/local/elkies-k3/fixed-reverse-4a1-rr-qq.json \
  --output artifacts/local/elkies-k3/fixed-reverse-4a1-pointing-qq.json \
  --surface-status PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN \
  --curves-status PASS_EXACT_QQ_FIXED_REVERSE_4A1_HORIZONTAL_ON_3A1 \
  --rr-status PASS_EXACT_QQ_FIXED_REVERSE_4A1_RR_JACOBIAN \
  --result-status PASS_EXACT_QQ_FIXED_REVERSE_4A1_POINTING \
  --edge '4A1/MW13 --q4 orbit114--> 3A1/MW14' \
  --old-i2-support "$(jq -r '.effective_horizontal_components[1].child_I2_support' \
    artifacts/local/elkies-k3/fixed-reverse-3a1-pointing-qq.json)" \
  --horizontal-roots-key effective_4A1_horizontal_roots_on_3A1_source \
  --horizontal-root-class-key class_in_3A1_coordinates \
  --expected-child-i2-count 4 --remaining-vertical-root-count 2 \
  --fixed-zero-source 'effective nonidentity component of the third old 3A1 I2 fibre'
```

The first three smooth chord systems have dimensions `14 -> rank 12 -> h0 2`.
Their minimal Jacobians have fibres `I2+22I1`, `2I2+20I1`, and `3I2+18I1`.
The first horizontal is a small word in the exact endpoint basis.  For the
second, 21 endpoint height-four sections meet the A1 fibre once; an integral
Smith solve selects seven pointed images.  For the third, 78 short A1 sections
span the full rank-15 2A1 tail and integral words select the orbit498
horizontal and two horizontal roots.  Old reducible-fibre components select
the quartic origins by exhaustive mod-131 exceptional-conic sign gates
followed by exact `QQ` invariant identities.  No Groebner basis or surface
elimination is used.

The q4/orbit114 replay first performs four exact physical affine-Weyl
reflections.  Its complete component, all-section, and finite-horizontal-wall
gates pass.  The physical divisor `O+P-C1-9F` has the compact rank sequence
`45 -> 3 -> 2`; exact square stripping gives a quartic and a minimal
`4I2+16I1` Jacobian.  The prescribed old component points the zero, two roots
are horizontal and two are vertical, and the reflected full NS transport is
unimodular.  This closes the last four fixed-corridor arrows.  The earlier
57-to-15-to-2 genus-two module remains a certified negative construction for
the unreduced divisor.

The next q4/orbit52 gate currently replays as follows:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_h92_fixed_reverse_5a1_physical_nef.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/construct_h92_fixed_reverse_5a1_abel_word_mod131.sage \
  --prime 167 --interpolate
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compact_h92_fixed_reverse_4a1_crossratio_qq.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/reconstruct_h92_fixed_reverse_5a1_compact_sections_crt_qq.sage
```

The first command is an exact `QQ` chamber/nef/marking certificate.  The
second is a regular finite-field seed for the old `P.O=15` section and three
horizontal roots, using fibrewise Abel reduction before taking the integral
lattice words.  The compact-model command sends the four rational `I2`
supports to `0`, `1`, `923/3815`, and infinity and applies an exact rational
Weierstrass isomorphism, reducing the model to 215-bit coefficients.  The CRT
command consumes the good-prime seed files listed in its output and
reconstructs all four sections over `QQ(t)`, with compact degree fingerprints
`(34/30,50/45)`, `(6/2,9/3)`, `(6/2,8/3)`, and `(10/6,14/9)`.  Exact equation
substitution and a withheld-prime-167 replay pass; the largest reconstructed
rational coefficient is 816 bits.

This proves the q52 horizontal and its three horizontal roots over `QQ`, but
does not yet prove the fixed `5A1 -> 4A1` equation arrow.  The remaining gate
is the exact resolved `D=O+P-C2-6F` two-plane, binary quartic, minimized
`5A1` Jacobian, and prescribed-zero/component pointing.  The optional
`lift_h92_fixed_reverse_5a1_sections_qq.sage` is retained only as a
coefficient-growth diagnostic: direct lifting in the former million-bit
normalization remained unresolved even at precision `167^8192`.  None of
these constructions uses a Groebner basis or surface elimination.

## E6+A1 rational-surface quadratic-base-change family

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-EQUATION 827d75cb8d14d7f4 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-ARITHMETIC-RANK2 387d6237125637a3 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-SPECIALIZATION-RANK7 bf1d025228805b31 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT96-A7D7-GALOIS ba008502f0e5533f -->

The exact construction and its integral K3 dissection replay with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rational_surface_quadratic_rank_search.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_k3_dissection.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_genuine_q2_neighbors.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit103_rr_weierstrass.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit96_rr_galois.sage
```

The first script verifies the rational surface, twist section, quadratic
pullback, height matrix, saturation, and complete degree-`(2,2)` polynomial
ansatz.  The second pins `NS`, `T=U(3)+<4>`, four singular-K3 boundary
lattices, the determinant-36 neighbor frame, and the complete smallest
nominal degree-two obstruction.  The third performs the complete genuine
norm-eight quadratic-neighbour census, proves nefness for eighteen MW-rank-3
frames, and selects the sparse divisor `P0+P1+A3_2` for equation compilation.
The fourth certifies the complete resolved basis `<1,z>`, eliminates it to a
binary quartic with rational origin, and produces the explicit
`2I1*+2I3+4I1` Weierstrass model with two polynomial rational points.
The fifth certifies that those two points are independent and exhaust the
`QQ(k)(r)` rank and exhibits the anti-invariant third geometric direction over
`QQ(sqrt(-3))`.  The sixth corrects the former orbit-96 coefficient-parent
error, attaches the physical E6 components, compiles the genuine
`I8+I3*+7I1` model, and proves that its MW representation is also
`1+chi_-3+1`, hence has arithmetic rank two.

Use the orbit-103 descent obstruction as a rational-specialization baseline
with

```bash
python3 elkies-k3/scripts/search_e6a1_orbit103_specializations.py

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_e6a1_orbit103_specialization_rank.sage \
  --lane small_coefficient_finalists --candidate-limit 1000 \
  --height 0 --pari-effort 0 --max-rank 24 \
  --reduction-prime-bound 500

python3 elkies-k3/scripts/test_e6a1_orbit103_specializations.py
```

The first command performs the complete two-stage trace search in its declared
height boxes.  The second accepts a new point only after a combined exact
mod-3 finite-quotient test against `Q_plus,Q_minus`.  The bounded result has
seven rank-at-least-seven fibres and no certified rank-at-least-eight fibre.
All 1000 selected compact candidates complete; a system-GP fallback handles
the 617 reproducible bundled-PARI precision errors.  See
[`../E6A1_ORBIT103_SPECIALIZATION_SEARCH_2026-09-02.md`](../E6A1_ORBIT103_SPECIALIZATION_SEARCH_2026-09-02.md).
See
[`../E6A1_RHO19_K3_DISSECTION_2026-09-02.md`](../E6A1_RHO19_K3_DISSECTION_2026-09-02.md)
and
[`../E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md`](../E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md)
and
[`../E6A1_RHO19_ORBIT103_WEIERSTRASS_2026-09-02.md`](../E6A1_RHO19_ORBIT103_WEIERSTRASS_2026-09-02.md).
See also
[`../E6A1_RHO19_ORBIT96_WEIERSTRASS_GALOIS_2026-09-02.md`](../E6A1_RHO19_ORBIT96_WEIERSTRASS_GALOIS_2026-09-02.md).

## E6 rank-three and rationalized D6 continuation

<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->

<!-- status-consumer: EC-K3-RES-QBC-E6-II-Q2-MW4 3aa5084463780acc -->

<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART a94042dd2d76797c -->

The exact `E6` rank-sum-three construction replays with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_ii_rank3_quadratic_base_change.sage
```

It proves the rational-surface rank two, twist rank one, generic K3 profile
`2E6+A2/MW3`, Picard rank 19, saturated determinant-24 `NS`, and an exact
Blichfeldt obstruction to every rootless MW17 fibration in that `NS`.

The complete first genuine quadratic-neighbour shell replays with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_e6_ii_rank3_q2_neighbor_candidates.sage
```

It enumerates `268` norm-eight Weyl orbits and proves that the six
minimum-root classes are nef `A6+D7/MW4` Jacobian fibrations.  The shell has
no frame of MW rank above four; equation compilation and arithmetic descent
of the MW4 directions are not part of this replay.

The rationalized `D6` marked-section frontier, exact rank-zero correspondence
obstruction, and retained height-30 regression replay with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_rationalized_d6_rank2_section_chart.sage
```

This command proves that the declared polynomial chart has no nontrivial
rational section pair.  It does not exclude a larger D6 rational-function
chart or a rank-sum-four family.

The simultaneous E6 shared-simple-pole modular census replays with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_e6_shared_pole_two_twist_sections_modp.sage \
  --prime 11 --enumerate --skip-msolve
```

It is complete for the displayed `GF(11)` ansatz and classifies every
survivor into two rejected mechanisms; it is not a characteristic-zero
rank-four theorem.  See
[`../E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md`](../E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md)
and
[`../LOWER_ROOT_TWO_TWIST_SEARCH_2026-09-02.md`](../LOWER_ROOT_TWO_TWIST_SEARCH_2026-09-02.md).

The systematic E6 node-collision linear-chord incidence replays with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_linear_chord_incidence.sage
```

It eliminates both complete linear-chord pencils, decomposes the unordered
base into geometric genus-zero and genus-two components, and parametrizes the
genus-zero quotient over `QQ`.  Recovering the ordered section labels gives
the genus-one cover `r^2=k^4+6*k^2+13`, birational to the rank-zero curve
`52a2`, so there is no affine ordered `QQ`-point.  The geometric rank-four
height calculation therefore lives over that ordered function field; the
descended `QQ(k)` family has exact arithmetic rank two.  The same checker
certifies the anti-invariant height determinant `52`.  Saturating by the two
literal chord sections changes the pure-character index by four and gives
geometric MW determinant `13/3`, generic geometric `rho=19`, and
`abs(det NS)=78`.  The rootless-MW17 Hermite screen passes.  It also emits the
full saturated integral NS Gram and verifies its unimodular split against
`data/lattice/e6_rank4_det78_frame.txt`.

<!-- status-consumer: EC-K3-RES-QBC-E6-RANK4-LINEAR-CHORD 3bcfe3534656b26f -->
<!-- status-consumer: EC-K3-E6-RANK4-ROOTLESS-Q2Q4-CENSUS 2351738f44774cfe -->
<!-- status-consumer: EC-K3-E6-RANK4-DET78-GLOBAL-ROOTFUL 648ec884ce7152bb -->

The complete zero-neutral rootless search through old degree four replays with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_rootless_low_degree_search.sage
```

It exhausts `80,123` Weyl-dominant classes, of which `79,837` are primitive,
and finds no rootless child in degrees two, three, or four.  The complete
genus-wide J2 frame classification replays with

```bash
sage -python classify_e6_rank4_det78_niemeier_frames.sage --check
```

Its complete 1,591-anchor Niemeier census yields 37,397 primitive
embedding-cover points and exactly 1,549 integral-isometry classes.  The exact
root-rank distribution is `10:1, 11:45, 12:249, 13:543, 14:477, 15:200,
16:33, 17:1`, and the reciprocal-automorphism sum closes the exact genus mass.
Thus the determinant-78 frame genus is rootful at `O(NS)`/J2 level and has
maximum MW rank 7.  The faster `--rootless-obstruction --check` mode is an
independent residual-rank cross-check, not the primary class certificate.  See
[`../E6_RANK4_DET78_NIEMEIER_CLASSIFICATION_2026-09-03.md`](../E6_RANK4_DET78_NIEMEIER_CLASSIFICATION_2026-09-03.md).

## Rootless genus theory

<!-- status-consumer: EC-K3-ROOTLESS-GENUS-MASS 2f5b874c0c22133b -->

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  certify_rootless_genus_first_moment.sage --check
```

`certify_rootless_genus_first_moment.sage` is an **ACTIVE_PROOF** for the
degree-one Siegel gate in Theorem H6.  It computes the exact mass-normalized
signed-root averages of the determinant-78, 948, and 950 rank-17 genera from
their bad-prime local densities.  It independently recomputes the
determinant-78 average from the complete 1,549-class census and obtains exact
agreement.  All three averages exceed two, so this checker records the cheap
criterion as inconclusive; it does not implement the higher-degree ADE
representation-mass inversion.  See
[`../ROOTLESS_GENUS_THEORY_2026-09-03.md`](../ROOTLESS_GENUS_THEORY_2026-09-03.md).

## Arithmetic rank transfer

<!-- status-consumer: EC-K3-ARITHMETIC-RANK-TRANSFER 3031dd2365a29cd5 -->
<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  certify_arithmetic_rank_transfer.sage --check
```

The checker implements Theorem A2 of
[`../RANK_MUTATION_AND_LIFT_THEOREMS.md`](../RANK_MUTATION_AND_LIFT_THEOREMS.md).
From an integral geometric NS Gram, finite Galois generators, and the marked
`U` and fibre-root basis of each fibration, it computes the exact rational
fixed NS, root, and Mordell--Weil spaces.  For every declared edge it verifies
`[M_2]-[M_1]=[R_1]-[R_2]` by traces on every group element.  The input schema
is
[`../data/arithmetic/arithmetic-marking-v1.schema.json`](../data/arithmetic/arithmetic-marking-v1.schema.json),
and `--marking FILE` validates an additional self-contained record.

Pinned controls reproduce arithmetic rank 17 for H3/R17 and promote the
degree-two `norm12-orbit-11952` alternate-Q80 pencil to exact arithmetic rank
17 before equation compilation.  The artifact stores its full rational
divisor basis `(F,O,Q1,...,Q17)`, the columns of `<D,O+D>`, the identity
Galois action, and the rootless quotient calculation.  The controls also give
arithmetic rank two inside geometric rank four for the unordered E6 incidence
and the `2+chi_-3` orbit-103 decomposition.  The same certificate applies the gate to
the current NS0024 completed-core path and returns
`FAIL_CLOSED_GEOMETRIC_ONLY`: the stored Kneser path is not a field-defined
marked-`U` corridor.  See
[`../ARITHMETIC_RANK_TRANSFER_2026-09-03.md`](../ARITHMETIC_RANK_TRANSFER_2026-09-03.md).
The alternate application proof is
[`../R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md`](../R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md).

## Integral rank-transfer glue calculus

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

```bash
python3 build_integral_rank_transfer_claim_provenance.py --check
sage -python certify_integral_rank_transfer_bridge_reglue.sage --check
sage -python certify_integral_rank_transfer_bridge_reglue.sage \
  --relative-u-output ../../artifacts/generated-results/elkies-k3-relative-u-bridge-lifting-regression-v1.json
sage -python certify_r17_local_bridge_mutation.sage --check
sage -python benchmark_integral_rank_transfer_bridge_predictor.sage --check
sage -python benchmark_e6_det78_prospective_bridge_predictor.sage --check
sage -python certify_integral_rank_transfer_theta_convolution.sage --check
sage -python certify_integral_rank_transfer_core_generation.sage --check
sage -python certify_integral_rank_transfer_reverse_theta_masks.sage --check
sage -python certify_integral_rank_transfer_weil_compression.sage --check
sage -python generate_integral_rank_transfer_masked_core_neighbors.sage --check
sage -python certify_integral_rank_transfer_q80_defect_birth_death.sage --check
sage -python certify_integral_rank_transfer_root_system_signature.sage --check
sage -python certify_inverse_ade_projective_birth_strata.sage --check
sage -python certify_integral_character_glue_calculus.sage --check
sage -python certify_r17_norm12_103b2_mw_glue.sage \
  --skip-specialization-saturation --check
sage -python build_integral_rank_transfer_glue_census.sage --check
```

`build_integral_rank_transfer_claim_provenance.py` checks that every labelled
theorem, corollary, proposition, lemma, and exact finite control in the
canonical theorem note has one provenance row, then regenerates the
machine-readable claim ledger.  It does not change mathematical status or
checker outputs.

The bridge checker certifies the rank-15 common core, rank-two cyclic bridge
replacement, graph glue, and exact root transfer for all 42 selected marked
corridor edges.  Its opt-in relative-`U` replay additionally constructs the
cross matrix `A`, verifies `Gram(w_1,w_2)=A^t J A-J`, recovers both saturated
bridges and the old-fibre degree, and checks the determinant square-index law
in both orientations on all 42 edges.  The local-mutation checker strengthens
this to equal two-sided saturation and glue order, verifies that coprimality
forces maximal glue on 35 of the 42 edges, and certifies the new published-R17
degree-two `4A1/MW13` fibration with maximal non-cyclic `ZZ/4+ZZ/8` bridge.
Its exact theta coefficients distinguish the new frame from both stored H3
`4A1` frames.  The character checker exhausts the integral E6 `2+1` and
`2+2` involution graph glues.  The norm-twelve checker gives the exact
`0x103b2` cover-level visible lattice and preserves rank at least 18 at the
specialization.  The pinned artifact also proves the displayed specialized
rank-18 subgroup primitive by isolating every possible eclib saturation prime;
the skip-mode byte check reuses that record.  These lattice/J2 checks
intentionally separate equation and rank-upper-bound claims.

The first relative-`U` application to the new completed NS0024 route is
documented in
[`../RELATIVE_U_BRIDGE_LIFTING_2026-09-03.md`](../RELATIVE_U_BRIDGE_LIFTING_2026-09-03.md).
`search_ns0024_relative_u_bridge_lifts.sage` rebuilds the four completed
frames, compares them with the known route, and performs exact ordered
representations of `G_A` when a full intersection box is requested.  Before
that more expensive two-vector step,
`search_root_adapted_weyl_neighbors.sage --summary-only` can exhaust the
possible first vectors modulo the source Weyl group while retaining exact
counts and root histograms but omitting individual negative witnesses.  The
recorded degree-two, degree-three, and degree-four first-edge boxes contain no
MW12 child, so no prospective primitive target `U` exists in those boxes.
After generating the three compact shell summaries, replay their combined
boundary and counts with

```bash
python3 certify_ns0024_relative_u_first_edge.py --check
```

The core-generation checker proves and replays the next inversion.  In every
maximal graph presentation the core form is exactly
`q_W orthogonal_sum (-q_C)` and its determinant is `det(W)*det(C)`; all 84
old/new presentations satisfy these identities.  Hence a proposed binary
bridge generates a rank-15 core genus before any `U` or equation is known.
For the four terminal controls, exhaustive even-genus enumeration at the
relevant determinant gives `48, 8, 8, 8` choices and the generated finite form
selects exactly one each time.
For a bounded bridge universe, the core's discriminant-coset theta table
through norm two is a complete completion signature.  The root-at-zero gate
rejects all 277 held-out E6 cores before graph enumeration; four known
minimum-four terminal cores remain positive controls.  Genus-class
enumeration, a bridge-determinant cutoff, and a speedup theorem remain open.

The reverse-theta checker derives from each bridge graph the exact core cells
that must vanish.  It makes the accept/reject decision through lazy exact CVP
queries of only those cells, before independently constructing the full theta
table as a truth check.  All 28 terminal graph decisions and root counts agree.
Sign symmetry compiles them to 14 nonredundant masks of 10--44 cells, and no
rank-17 child is constructed.  The theorem
also bounds every nonconstant norm-at-most-two coefficient of a rootless
rank-15 core by 30, producing a finite allowed-signature sieve.  Realizing
only lattices with those signatures remains the open inverse-theta step.

The zero-orbit Weil checker compresses the modular feasibility stage without
changing any reverse-mask coordinate.  It reconstructs each terminal
`q_W`, computes `O(q_W)` and its exact orbit quotient, and certifies the
`S,T`-cyclic rank of the zero class at a good split prime.  Theta symmetry
reduces the four full coefficient dimensions
`16,560, 44,556, 181,450, 21,804` to
`864, 5,760, 24,960, 2,880`.  In all four controls the zero-generated cyclic
submodule equals the orthogonal-invariant quotient, so no additional cyclic
compression is available.  The checker does not yet enumerate modular forms
or realize new lattice cores.

The checker also evaluates the compressed invariant Riemann--Roch trace
formula.  It certifies modular dimensions `476, 3,121, 13,488, 1,563` and
cusp dimensions `472, 3,120, 13,485, 1,562`.  Since the fourteen masks have
only 10--44 cells, rank--nullity leaves cusp-kernel dimensions of at least
461 for every mask.  This is an exact negative result for linear modular
screening, not a decision about affine normalization or arithmetic theta
realizability.

`generate_integral_rank_transfer_masked_core_neighbors.sage` closes the first
bounded core-generation loop.  It starts from the canonical rootful
representative of the finite-form-forced Golay-720 rank-15 genus, not from a
historical core.  A stored seven-step good-prime neighbour path reaches a new
minimum-four core class with reverse-mask profile `(3,0)`; exact graph glue
then gives the declared rootless rank-17 target.  The normal `--check` mode
replays this short certificate.  Add `--search` to rerun the deterministic
34,571-neighbour beam whose exact score is occupied reverse-mask support.  It
is a positive bounded control, not a complete genus traversal or complexity
claim.

`search_integral_rank_transfer_masked_core_controls.sage` applies the core
generator to H3, NS0024, and Q80.  `--root-descent` removes the
Golay-specific four-root frontier cutoff; `--mask-cap 3` truncates only
nonzero ranking counts while leaving zero acceptance exact; and
`--support-diversity` reserves elite slots for distinct occupied masked
support signatures.  Rootful bridges are discarded before scoring because
their compulsory `(0,0)` mask cell can never vanish.

`certify_integral_rank_transfer_masked_core_controls.sage --check` is the
short deterministic theorem replay.  It reconstructs new H3 and NS0024 core
classes from paths of lengths eight and three, verifies their zero masks, and
constructs rootless rank-17 completions.  It also reconstructs Q80's best
two-cell near miss.  The long driver found NS0024 after 7,477 raw neighbours;
the corrected H3 and Q80 support-diversity controls each miss within 42,300
raw neighbours.  Those misses are bounded experiments, not genus-wide
obstructions.

`search_integral_rank_transfer_q80_defect_neighbors.sage` expands every
occupied Q80 mask cell into its exact physical dual vectors.  It samples only
isotropic lines nonorthogonal to all of them, which by the dual-neighbour
survival lemma removes every old witness.  Its 10,000-line control finds
1,397 rootless directed neighbours and no zero mask: all acquire replacement
defects.

`search_integral_rank_transfer_q80_defect_beam.sage` repeats this operation
and keeps at most one lowest-defect representative from each integral
isometry class.  It reaches zero support in generation four after 30,228
constructed neighbours.  `certify_integral_rank_transfer_q80_defect_completion.sage
--check` is the short exact replay.  It verifies the four witness-removal
steps, the replacement counts `4,6,4,4,0`, and the new rootless determinant-948
completion.  It also compares both rootless controls, counts 1,313 norm-four
pairs, computes automorphism-group order four, checks the exact local genus
symbols, and stores an explicit integral isometry to the alternate Q80 frame.
The historical `declared_target_frame` in this corridor is published R17.
Defect cardinality is not monotone; physical witness incidence is the
operative transition datum.

`certify_integral_rank_transfer_q80_defect_birth_death.sage --check`
implements the missing half of that datum.  It decomposes the child dual into
the affine layers `K_y^dual+j*y/p` and queries the equivalent masked cosets
`M+r+k_0+j*y/p` before constructing the child.  It also evaluates every
theta cell through norm two, giving a complete `Sigma_2` transition.  The
later materialized-child enumeration agrees profile-for-profile and
vector-for-vector and reproduces `4,6,4,4,0`.  This is a correctness and
regression certificate; it does not claim that the abstract
counted theta signature determines line pairings or that this implementation
has a uniform timing advantage.

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-GRAPH-REACHABILITY e02f950eba79b32a -->

`analyze_small_genus_defect_graphs.sage --check` is the complete finite-graph
control for reachability.  It mass-closes positive even ternary genera of
determinants 112, 126, and 316, enumerates every isotropic line at good primes
`{3,5,11}`, `{5,11,13}`, and `{3,5,7}`, stores every signed physical root and
birth/death count, and computes all singleton, pair, and triple unions.  Every
directed SCC is labelled with its condensation exits; every finite shortest
path stores a prime and isotropic-line witness.  Exhaustive subset selection
gives minimum sufficient sets `{5}`/`{11}`, `{5}`/`{11}`/`{13}`, and
`{5}`/`{7}`.  It exhibits fixed-`3` traps even though the unrestricted graph
is connected and each genus has one proper spinor genus, as well as the exact
distance-two defect sequence `2,2,0`.  It is a root-defect calibration, not a
complete rank-15 reverse-mask graph.

`certify_integral_rank_transfer_root_system_signature.sage --check` expands
the physical completion witnesses into complete root-line metrics.  For the
four NS0024 stages it verifies every identity
`<k+c,k'+c'>=<k,k'>+<c,c'>`, classifies the metric components as
`D5+E8`, `3A1+A2`, `3A1+A2`, and rootless, and computes the marked primitive
closures and torsion quotients.  The middle rootful completion with a
rootless core has twelve roots in five nonzero graph-glue labels.  Existing
Q80 `4A1` and `A1` frames are independent target-classification controls.
Pairwise root products determine ADE data, but the retained frame coordinates
are essential for primitive closure and exact Mordell--Weil torsion.

`certify_ns0024_inverse_ade_mutation.sage --check` composes that metric
signature with the good-prime dual-layer law to invert one ADE transition
before child construction.  For the first NS0024 `p=17` edge it records the
140 parent physical root lines and their modular linear forms on `y`: exactly
six must vanish and 134 must be nonzero.  It then exhausts the affine-CVP
layers joined through the fixed order-191 bridge, finds no births or extra
roots, and recognizes the six survivors as `3A1+A2`.  Materializing the child
afterward gives the identical physical root set.  This is an exact inverse
predicate and control, not yet an all-line solver or a complexity result.

`certify_inverse_ade_projective_birth_strata.sage --check` eliminates the
affine variable by scaling a born vector `v` to the parent dual-shell vector
`z=p*v` and retaining its projective reduction modulo `p`.  It exhausts all
346 isotropic lines in 48 state/prime cases from the three mass-closed
ternary controls; every predicted root set equals the independently
materialized child root set, and the forbidden-stratum complement gives all
192 rootless lines.  It also exhausts a six-line index-two graph-glue control
where every line has 16 births in the nonzero glue coset; the predicted and
materialized completed-child root sets again agree, for 352 exact comparisons
overall.  Expanding a rank-15 scaled shell has no claimed complexity bound.

`derive_r17_genus_one_bisection_twist_section.sage` descends the certified
`0x103b2` split bisection to an exact height-eight section on the quartic
twist. `export_elkies_2026_twist_polynomial_sections_modp.sage` accepts its
`--genus-one-label`; `run_twist_polynomial_sections_bruteforce.py` compiles
`bruteforce_twist_polynomial_sections_modp.cpp` and exhausts the exported
finite-field `P.O=0` blocks by coefficient enumeration, value interpolation,
or meet-in-the-middle bitsets. Finally,
`hensel_lift_r17_103b2_po0_sections.sage` audits the seven reduced `p=29`
branches through precision `29^800`. Only the known branch reconstructs
exactly. This is a short-shell result, not a full twist-rank upper bound; see
[`../R17_103B2_ANTI_INVARIANT_RANK_AUDIT_2026-09-03.md`](../R17_103B2_ANTI_INVARIANT_RANK_AUDIT_2026-09-03.md).

The predictor benchmark replays five H3 first-hit histories without
enumerating candidate-child roots.  The exact `K+C_new` lower-bound screen
rejects only 178 of 2,892 candidates, so it is preserved as a negative
experiment rather than promoted to a construction algorithm.  It also
exhausts 14 binary bridge classes on the four observed rootless terminal
cores.  Maximizing bridge minimum retains five classes, four rootless, versus
five rootless classes without screening (`2.24x` precision enrichment and
`80%` recall).  This second result is explicitly selected-core evidence, not
an out-of-sample q-neighbor benchmark.

`benchmark_e6_det78_prospective_bridge_predictor.sage` supplies the first
genuinely held-out shell.  It generates all 277 primitive zero-neutral
old-degree-two children of the determinant-78 E6 source before consulting the
mass-closed 1,549-class truth catalogue.  Every candidate has the same least
nonzero glue-coset minimum two, although root ranks range from 12 to 15.  This
is an exact negative control for that scalar score, not a rootless-positive
predictor benchmark.

`certify_integral_rank_transfer_theta_convolution.sage` implements the exact
replacement for that scalar score.  It caches all dual-coset theta
coefficients through norm two for the four terminal cores.  From their cyclic
bridge determinants it independently enumerates all fourteen reduced positive
even binary bridges, derives all 28 oriented graph multipliers from
finite-form isotropy, and convolves their tables.  It recovers every child
root count and all five rootless classes without constructing a rank-17 child
during prediction.  Its scope is the complete declared fixed-core universe;
it does not claim a speedup or automatic core generation.

<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-TWO-TWIST-POLYNOMIAL ea0496c9566cfdc3 -->

<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-LOW-SLICE-ELIMINANTS 43d297285eb3655b -->

The section-first D5 seed and its complete polynomial two-twist census replay
with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_d5_two_marked_two_twist_polynomial_modp.sage \
  --prime 11
```

The exact seed has an `I1*` fibre and two independent invariant sections,
certified by specialization to a rank-two basis over `QQ`.  The `GF(11)`
census finds one two-section quadratic twist in the declared degree-two
polynomial chart.  The regular low-section lift through that pair is excluded
over `QQ` by the exact eliminant below; the full polynomial chart remains
open.

`lift_d5_two_marked_two_twist_low_section_slices.sage` certifies the regular
local slices through the `p=11,13` survivors.  With `--run-eliminants` it
orders `t` last, eliminates the other seven variables exactly over `QQ`,
factors the resulting univariate polynomial, and selects the target factor by
the declared modular residue.  The two target fields have degrees 88 and 78;
neither complete saturated slice has a rational point.  Use an eliminant
timeout of at least 900 seconds on the current reference host.

<!-- status-consumer: EC-K3-RES-A4-TWO-POINT-TATE-SLICE-OBSTRUCTION b9729a0a8f2f17be -->

The simplest two-point A4 Tate slice is closed exactly by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_a4_two_point_tate_slice_obstruction.sage
```

On both `I5` branches its residual discriminant has a forced square linear
factor.  Removing the extra repeated fibre, generically `I2`, makes the two
marked points dependent.  The result is scoped to the normalized `r=s=1`
slice.

## Published-R17 genus-one bisection target pilot

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-BISECTION-PILOT 80fa6e59107cc9e6 -->

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_genus_one_bisections.sage \
  --check
```

This active proof entry point completely ranks the minimum-norm-eight trace
frontier, then uses the cheapest finite-pole trace to fit and verify the
genus-one pencil member through each of the eleven exact rank-28 targets.  It
certifies quartic branching, exact Kummer barcodes, lifted sections and
height-16 anti-invariant independence.  Its declared boundary is one trace
template and one known positive-control fibre, not a rank-32 search.

## Published-R17 multisection-visibility filtration

<!-- status-consumer: EC-K3-R17-MULTISECTION-VISIBILITY-FILTRATION 2f41e9f4236f6c9e -->

```bash
sage -python \
  elkies-k3/scripts/certify_r17_multisection_visibility_filtration.sage \
  --check
```

This reuses the exact common norm-eight trace and fits its genus-one
bisection pencil through all 38 displayed exceptional generators at the four
rank-25--28 controls.  It records the literal all-genus filtration, which is
already full in degree two, separately from the finite rational-curve
filtration.  The latter has a complete degree-two atlas but only bounded
degree-three and degree-four equation subatlases.

## Published-R17 frozen-quartic simultaneous splitting

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-SIMULTANEOUS-SPLITTING-H10000 40fb0bc465e3e95c -->

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_simultaneous_splitting.sage \
  --check
```

This active bounded search compiles
`scan_elkies_2026_rank28_genus_one_splitting.cpp`, exhaustively scans every
primitive `t=a/b` with `|a|,b <= 10000`, and exactly tests all simultaneous
square survivors.  It also enumerates `2P,...,30P` for the canonical
opposite-ordinate point on each of the eleven pointed quartic Jacobians.  The
only compact simultaneous split is the original positive control and there
are no subgroup cross-splits.  The result is not a global theorem about the
pairwise fibre products.

## Published-R17 mixed-trace genus-one splitting

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-MIXED-TRACE-SPLITTING-H10000 c7aa09836b842b60 -->

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_mixed_trace_splitting.sage \
  --check
```

This active bounded search constructs all eleven barcode-fitted quartics for
each of the seven equation-cheapest finite-pole norm-eight traces.  Its C++
hot loop retains only square survivors belonging to at least two distinct
trace pencils.  Among `121,589,943` primitive parameters through height
`10,000`, exact testing recovers only the original 77-fold positive control.
It does not search the rest of the norm-eight trace frontier.

## High-throughput genus-one bisection splitting

<!-- status-consumer: EC-K3-R17-GENUS1-HIGH-THROUGHPUT-SPLITTING cad3d98ce58c89e7 -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-MW-LATTICE-SIEVE aa0d8718eb57de6f -->

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_elkies_2026_genus_one_bisection_splitting.sage \
  --norm8-count 100 --equation-pool-size 1000
```

This active bounded search selects a structurally diverse equation-cheap
norm-eight sample, adds all 43 norm-twelve deep classes, constructs pointed
quartics without public exceptional targets, and drives
`scan_elkies_2026_genus_one_bisection_splitting.cpp`.  The C++ mask width is
dynamic and each parameter is ranked by weakest-prime-block cover count and
trace-mask rank.  Exact survivors are mapped back to the R17 fibre and tested
modulo the specialized generic MW17 subgroup.  The pinned production run
finds the exact `t=1/25` norm-eight/norm-twelve collision and one escaping
norm-twelve direction; its no-other-hit statement is bounded.

`search_r17_norm12_103b2_mw_lattice.sage` uses that new point to construct the
global minimal Jacobian of `0x103b2`.  A height-10000 quartic-point search gives
17 independent points, proving rank at least 17; the attempted exact PARI
upper bound and final `40251553` saturation were stopped and are not claimed.
The active checker instead sieves signed sparse coefficient vectors in this
known subgroup against the other 142 covers at 32 good primes and exactly
replays every survivor.  The radius-two/support-five and
unit/support-eight shells have 16,496,324 vectors in their union and recover
only the known `t=1/25`/`0x0f6b1` overlap.  This is a bounded cover-level
search, not a full Jacobian-group or K3 specialization-rank theorem.

```bash
sage -python search_r17_norm12_103b2_mw_lattice.sage
sage -python search_r17_norm12_103b2_mw_lattice.sage \
  --max-support 8 --coefficient-radius 1 \
  --output ../../artifacts/generated-results/\
elkies-k3-r17-norm12-103b2-mw-lattice-unit-support8-v1.json
```

`search_r17_norm12_103b2_hard_fibre_products.sage` isolates the seven
partners that survived all 32 local tests somewhere in those shells.  For
each partner it searches the necessary genus-three quotient
`z^2=f_103b2(t)f_partner(t)` exactly with PARI.  At naive height `300000` all
seven searches have zero affine quotient points, hence zero simultaneous
splits.  This remains a bounded point search, not a global rational-point
determination.

<!-- status-consumer: EC-K3-R17-NORM12-103B2-HARD-FIBRE-PRODUCT-H300000 b4fef7ab54b922e0 -->

```bash
sage -python search_r17_norm12_103b2_hard_fibre_products.sage --height 300000
```

`classify_r17_103b2_isotropic_frame.sage` performs the exact lattice test for
the genus-one bisection itself.  In `NS=U+R17(-1)` it forms
`D=(3,2,w_103b2)`, splits off an integral `U`, enumerates every norm-two root
of the orthogonal positive frame, and compares that frame with both certified
rootless determinant-948 `J2` classes.  The result is rootless/MW17 in the
published R17 frame class, not the alternate Q80 class.

`control_pointed_cover_jacobian_ranks.sage` runs the literal
`hyperellratpoints(H=10000)` to pointed-map to eclib-relation pipeline on a
seeded uniform sample of ten other pointed covers, with `0x103b2` replayed by
the same code.  All ten controls have only their signed pointing pair and no
nonbase mapped point; `0x103b2` has 58 nonbase images of relation rank 17.
The controls retain their separately certified rank-one generator, and the
zero visible nonbase ranks are not rank upper bounds.

<!-- status-consumer: EC-K3-R17-NORM12-103B2-ISOTROPIC-FRAME 47f3a0eb7ee50bcb -->
<!-- status-consumer: EC-K3-R17-POINTED-COVER-JACOBIAN-CONTROL-H10000 4bb087b3a1ebc684 -->

```bash
sage -python classify_r17_103b2_isotropic_frame.sage
sage -python control_pointed_cover_jacobian_ranks.sage
```

`classify_r17_norm12_isotropic_frames.sage` applies the same exact splitting
and two-control isometry test to all 43 norm-twelve genus-one bisections.  It
finds 33 published-frame copies and ten alternate-Q80 copies.  Every alternate
copy has old-fibre degree two, shares the old zero, and has zero-section degree
one.  Since distinct `J2` classes cannot have fibre intersection below two,
this proves the exact minimum accessibility distance.  The cheapest stored
alternate witness is `norm12-orbit-11952`.

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-MINIMAL-ACCESSIBILITY 631f50389e0a3283 -->

```bash
sage -python classify_r17_norm12_isotropic_frames.sage
sage -python classify_r17_norm12_isotropic_frames.sage --check
```

<!-- status-consumer: EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION 077c6409d76cbe63 -->

`compile_r17_norm12_orbit11952_qq.sage` performs the direct classical
two-neighbour hop for the cheapest alternate-Q80 witness.  It solves the
published-model `H^0(O(D))` kernel, constructs and points the quartic,
minimalizes its Jacobian to a `24I1` K3, identifies the rootless complement
with alternate Q80, and exports a saturated rank-17 section basis.  This is
the primary equation compiler for the alternate frame; the giant historical
Q80 transport below is a fallback.

```bash
sage -python compile_r17_norm12_orbit11952_qq.sage
sage -python compile_r17_norm12_orbit11952_qq.sage --check
```

<!-- status-consumer: EC-K3-R17-NORM12-11952-CONTROL-J-PREIMAGES 1ef38474a0d7f629 -->

`certify_r17_norm12_11952_control_j_preimages.sage` tests whether the four
published-R17 rank-25--28 control curves occur at rational parameters of the
alternate family.  It factors the exact degree-24 cross-multiplied `j`-preimage
polynomial for each control and checks both finite and infinite rational roots.
All four have no rational alternate preimage, so alternate-native calibration
fibres are required.

```bash
sage -python certify_r17_norm12_11952_control_j_preimages.sage
sage -python certify_r17_norm12_11952_control_j_preimages.sage --check
```

## Alternate-Q80 rootless equation handoff

`build_q80_alternate_final_divisor_handoff.sage` is the fail-closed
historical fallback input to the alternate-rootless equation compiler.  It
replays the physical nef final q6 fibre and zero, full determinant-one NS
transport, rootless determinant-948 MW17 child, and all 1,313 norm-four
pairs.  Its artifact deliberately records the immediate `A1/MW16`
characteristic-zero equation and the resolved function basis as open.  The
preferred construction now compiles `norm12-orbit-11952` directly on the
published R17 equation:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_q80_alternate_final_divisor_handoff.sage --check
```

`audit_q80_third_q12_descent_field_normalization.py` is an
**ACTIVE_COMPILER** arithmetic audit at the true upstream equation frontier.
It proves directly from the exact closure operands that the reduced numerator
of `q1*q2` is a square, and therefore replaces the rational-radicand generator
by `delta^2=denominator(q1*q2)` without attempting factorization:

```bash
python3 \
  elkies-k3/scripts/audit_q80_third_q12_descent_field_normalization.py --check
```

The normalization is not an equation reconstruction and does not assert that
the transformed coefficients are smaller.  See
[`../Q80_ALTERNATE_ROOTLESS_EQUATION_HANDOFF_2026-09-03.md`](../Q80_ALTERNATE_ROOTLESS_EQUATION_HANDOFF_2026-09-03.md).

`audit_q80_third_q12_pencil_basis_heights.py` performs that missing bounded
height test on all 63 exact moving-equation coefficients:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_q80_third_q12_pencil_basis_heights.py --check
```

The exact result rejects `delta` as a raw-coordinate simplification: it is
worse than `omega` on 54 terms and raises the maximum by 5,929 bits.  After
rational projective normalization it gives only a 10,888-bit (about 0.7
percent) improvement at a still 1,484,751-bit primitive maximum.  Integer
content is only 7--12 bits.  The script does not optimize multiplication by a
general quadratic-field element, base `PGL2`, or model transformations.

## Rule for future additions

Every new script should state near its header:

```text
status: ACTIVE_PROOF | ACTIVE_COMPILER | ACTIVE_SEARCH | REGRESSION | HISTORICAL_DIAGNOSTIC
claim: the exact claim or bounded search it supports
inputs: pinned files/certificates
outputs: generated artifact path
supersedes/superseded-by: optional script or note
```

A successful search result becomes a proof entry point only after an independent replay
checks its exact divisor, chamber/nefness, equation identity, and claimed fibre or
Mordell--Weil data.
