# q8/orbit376 RR branch audit

Date: 2026-08-26  
Scope: draft PR #35, branch `research/q8o376-rr-p2-scan`

## Audit boundary

This is a static audit of the mathematics, schemas, dependency chain, and all
scripts added by the branch. The ignored `artifacts/local/elkies-k3` inputs are
not available through the GitHub connector, so none of the Sage calculations
has been executed in this environment. A syntactically plausible exact check
is therefore classified separately from a replayed result.

The branch must remain draft. No q8/orbit376 equation status should be
promoted from this branch until the blocking items below are fixed and the
complete Sage replay succeeds on the repository machine.

## Executive result

The post-collision reduction is promising: conditional on the stated chord
module, the quotient by constants really is represented by a projective
quadratic `BB`, so the remaining modular discovery space is only
`P^2(GF(p))`. The CRT/projective-LLL and exact binary-quartic strategy is also
reasonable.

However, the current branch has one definite algebraic sign error in the
pointing stage, several provenance and cache weaknesses, and a larger proof
boundary: an exact semistable `4A1` pencil is not yet proved to be the marked
q8/orbit376 pencil. Even a successful unpointed reconstruction would currently
prove less than its PASS status says.

## Definite defects

### 1. Chord sign conventions are inconsistent

The collision congruence and fraction-free radicand used by the modular and
QQ compilers are consistent with

```text
m = (y_H - y)/(x - x_H).
```

For this convention the quadratic discriminant is

```text
m^4 - 6*x_H*m^2 - 8*y_H*m - 3*x_H^2 - 4*A.
```

Indeed, with `N=AA-u*Z^2` and `Db=-BB`, the compiler has
`m=N/(Z*Db)` and the stored `-8*Y*N*Db^3` term is correct.

The docstrings instead state

```text
m = (y - y_H)/(x - x_H),
```

whose discriminant has the opposite linear term,

```text
m^4 - 6*x_H*m^2 + 8*y_H*m - 3*x_H^2 - 4*A.
```

More importantly,
`point_h92_q4o164_q8o376_from_known_sections_qq.sage` currently uses

```python
chord = (y_signed - yH) / (x_value - xH)
```

while comparing against the quartic compiled with the first convention. Its
literal surface-to-quartic identity should therefore fail. The pointing code
must use

```python
chord = (yH - y_signed) / (x_value - xH)
```

or the entire collision/radicand convention must be changed coherently. The
probe and reconstruction docstrings should be corrected at the same time.

This is a blocking correctness issue, not merely notation.

### 2. Superseded v1 entry points were runnable and wrote the same artifact

The first P1229 certificate consumed nonexistent pointing keys
`degree_one`/`candidates`; the actual artifact emits
`degree_one_candidates`. Its companion runner called that broken script. Both
v1 and v2 certificates wrote
`q4o164-q8o376-p1229-origin-qq.json`, so accidentally running the old entry
point could overwrite a valid v2 result with an OPEN artifact.

The two obsolete files have been removed from the branch during this audit:

- `certify_h92_q4o164_q8o376_p1229_origin_qq.sage`
- `finish_h92_q4o164_q8o376_pointing_stack.py`

Only the `_v2` finishing path should remain until it is renamed canonically in
a later cleanup.

### 3. Cached modular scans are not tied to the code or current inputs

The parallel runner reuses a scan from only its schema, prime, and
`complete=true`. It does not compare the artifact's input hashes with the
current model, horizontal, and cost files, and it does not record/check the
probe script hash. A stale v2 artifact produced before a code correction can
therefore be silently reused.

The QQ reconstructor checks that all scan artifacts share one input-hash
mapping, but does not check that mapping against the current files it loads.
Thus a self-consistent stale scan set can be combined with a newer exact
parent.

Every consumer must recompute and compare current hashes. The modular schema
should also be bumped whenever the selector or branch formula changes.

## Major mathematical proof gaps

### 4. A `4A1` pencil is not yet the marked q8/orbit376 pencil

The modular selector asks for square vertical content and a semistable `4A1`
Jacobian. The QQ stage then proves that the reconstructed line has those
properties. It does not verify the target's exact vertical-root correction or
its class in the C8-pointed Neron--Severi lattice.

There can in principle be another line in the post-collision `P^2` with the
same ADE profile. Even finding `P1229` as a degree-one curve only identifies a
possible zero; it does not identify the embedded fibre class or full marked
`U`.

Until an exact class check exists, these claims are too strong:

- `...Q8O376_UNPOINTED_RR_AND_4A1_JACOBIAN`;
- `resolved_RR.divisor_target = q8/orbit376`;
- `route_lock.matches_promoted_route = true`.

The exact result should be described as an unmarked `4A1` pencil/Jacobian
candidate. Marked promotion requires either resolved component valuations or
an exact intersection fingerprint on enough explicit curves to solve for the
fibre class, followed by the certified unimodular NS transport.

### 5. The complete Riemann--Roch space is not proved

The branch constructs the constant section and one selected nonconstant
section. It does not explicitly construct two independent linear vertical
conditions or prove that the complete target `H^0(D)` is exactly this
2-plane. The values

```text
vertical_condition_rank = 2
kernel_dimension = 2
```

are inferred from the planning score plus the selected square-content line.
They are not yet an exact component-valuation rank calculation.

A found pencil can still be useful, but the artifact should not call it the
complete resolved RR space until those local conditions are certified.
Conversely, a failed `P^2` scan would not prove the marked edge impossible,
because the global chord-module ansatz itself has not been proved exhaustive
for the resolved divisor.

### 6. Pointing does not yet supply the full marked equation edge

The known-section pass verifies quartic identities and `c4/c6`, but it does
not serialize a complete birational map from the original K3 to the pointed
short child. The P1229 certificate identifies a zero class only. It does not
attach the four child `I2` component labels or the complete equation-to-lattice
frame.

The old v1 proof boundary stated this limitation more accurately than the v2
text. The v2 `matches_promoted_route=true` must be weakened to
`matches_route_selected_zero=true` until the fibre/components are attached.

## Additional implementation risks

### 7. Embedding 15 is hard-coded but not locked to the horizontal artifact

The marked-class matcher chooses `embedding_index=15` and checks that this
embedding is compatible with the old audit. It does not verify that the
currently loaded exact horizontal was reconstructed from the modular identity
that selected embedding 15. The horizontal artifact is loaded but only its
status is inspected.

Lock the selected identity

```text
H=T-C8opp-B0+2*B1+B2-3*B3-B4-2*B5+B7
```

and its embedding/index provenance explicitly.

### 8. Too few construction primes may make reconstruction unreliable

The runner permits two construction primes and one holdout. For primes near
40--60 this is only about eleven or twelve CRT bits. Even the default four
construction primes provide only a few dozen bits. This can cause rational
reconstruction to fail despite a correct direction.

The projective LLL fallback examines reduced basis rows only; it does not
enumerate all sufficiently short lattice vectors or combinations. This is a
false-negative risk rather than a false proof, but the runner should add
primes adaptively and use short-vector enumeration in the tiny 3-dimensional
lattice.

### 9. One bad or ambiguous prime aborts the entire parallel run

A denominator prime, degree-drop prime, or prime with an accidental second
`4A1` direction raises from its future and aborts the stack. The driver should
record rejected primes and continue until enough good primes in one stable
chart remain.

### 10. Direct reconstruction accepts stale/older scan shapes

The standalone reconstructor has fallback support for `scanned` and does not
require the v2 schema or a recognized candidate status. The runner is stricter,
but direct invocation through the default glob can consume older artifacts.
Require the exact schema, complete flag, prime validity, and expected input
paths/hashes.

### 11. Exact replay compares raw coefficient lists

Construction/holdout replay compares finite-field `A` and `B` coefficient
lists literally. At a prime where a leading coefficient vanishes, one side
may retain a trailing zero while Sage's polynomial `.list()` on the other side
trims it. This can reject an otherwise valid held-out prime. Canonicalize by
padding to the expected degree before comparing.

### 12. P1229 identification status is slightly over-eager

The marked-class artifact emits the special P1229 status whenever there is
exactly one direct match, even if a separate transported section matches
`-P1229`. The downstream certificate correctly checks the combined direct and
negative count, but the intermediate status should use that same total-count
gate.

### 13. The q12 shortcut should remain conditional

It may be possible to start q12/orbit5867 from only the P1229-pointed child,
but this has not been established. Its proposed four `P.O=0` section word and
marked target may require the child I2 labels to identify branches. For I2
fibres those labels should be relatively cheap: hitting the nodal point
already detects the nonidentity component. Skipping all component attachment
should be treated as an experiment, not as a certified shortcut.

## Parts that survive the static audit

Subject to the sign convention and exhaustiveness caveats above, the following
ideas are coherent:

1. `12 -> 4` by the congruence `AA*X=BB*Y mod Z^2`.
2. Quotienting the collision kernel by constants leaves the three `BB`
   coefficients, hence projective `P^2`.
3. The fraction-free chord expression and exact `Z^6` removal are consistent
   with `m=(y_H-y)/(x-x_H)`.
4. Dividing a verified square content and compiling binary-quartic invariants
   is a sensible way to obtain the child Jacobian.
5. The semistable finite/infinity fibre audit is structurally appropriate in
   characteristic different from two and three.
6. Multi-prime projective reconstruction plus an independent holdout is a
   good discovery-to-proof architecture once provenance is locked.
7. Fourfold pole growth and polarization against the rank-eight equation
   basis is a strong way to identify transported sections in marked MW9,
   provided the selected embedding is explicitly locked.

## Required promotion gates

Run these in order:

1. Fix the chord sign inconsistency and bump the modular artifact schema.
2. Add current-input and script-hash validation to scan reuse and every exact
   consumer.
3. Run smoke tests, then collect enough good primes adaptively. Keep pointing
   disabled until the sign fix is present.
4. Reconstruct over QQ and verify at several construction primes and at least
   two independently generated holdouts.
5. Verify the exact marked fibre class, not only ADE, using component
   valuations or a spanning explicit-curve intersection fingerprint.
6. Prove the complete resolved RR rank `12 -> 4 -> 2` for that marked divisor.
7. Point with the independently identified P1229 curve and serialize the exact
   birational/Weierstrass maps.
8. Attach enough child component data to import the certified NS transport.
9. Only then promote q8/orbit376 and begin the marked q12/orbit5867 equation
   lift.

## Safe current claim

Before replay, the branch provides an **experimental modular and exact
compiler pipeline** for discovering a low-dimensional `4A1` pencil from the
known q8 horizontal. It does not yet provide a certified q8/orbit376 equation
edge.
