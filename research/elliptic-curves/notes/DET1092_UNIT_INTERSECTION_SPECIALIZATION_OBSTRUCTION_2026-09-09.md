# Unit intersections predict inherited conic splits from generic data

## Result and classification

**New deduction/application of established intersection theory.** For a
minimum norm-ten bisection and an inherited section, a height-lattice
unit-intersection condition forces a rational conic point. If its original
elliptic fibre is smooth, both split branches are inherited. Thus this
condition constructs an exact specialization obstruction, not a new seed.

**Verified application.** One fixed generic conic, orbit8044, and578 signed
words supported on at most two displayed generators give14 unit
intersections. All14 are distinct, finite, unramified over the original
parameter line, and lie on smooth elliptic fibres. All28 branch points
are certified explicit generic words. No exceptional point, control
parameter or search artifact was a construction input.

The predicted list contains the old dependent control **`s=-528/3635`**.
That address was read only during a separate comparison after construction
and independent replay were sealed. This explains its failed seed incidence
as a forced intersection of two globally defined curves, independently of
the previous halving-cycle diagnosis. The other13 addresses are further
exact obstructions for this conic, not assertions that their elliptic
fibres have rank17 or lack other seeds.

The rule was developed retrospectively and uses an existing generic orbit.
Its execution is equation-only; this is not a claim of a previously frozen
blind discovery experiment. The578-word support restriction is coordinate
dependent. The intersection criterion itself is intrinsic. The number14
here has no asserted connection to the14 exceptional directions on302.

## 1. Audit and the intrinsic criterion

**Established literature.** The Shioda map identifies the negative
orthogonal intersection form with the Mordell--Weil height pairing; see
[Schuett--Shioda, *Elliptic Surfaces*, §§11.4--11.8](https://arxiv.org/pdf/0907.0298).
For this24-I1 K3 there are no reducible-fibre correction terms. Write

\[
 C_w=2O+4F+\phi(w),\qquad
 S_v=O+\frac{\lVert v\rVert^2}{2}F+\phi(v).
\]

Using `O^2=-2`, `O.F=1`, `F^2=0` gives

\[
 \boxed{C_w\cdot S_v=\lVert v\rVert^2-\langle w,v\rangle.}       \tag{1}
\]

The [index theorem](DET1092_RATIONAL_BISECTION_INDEX_AND_ODD_DIVISOR_CONSTRUCTION_2026-09-09.md)
already used this formula for displayed basis sections to construct
odd-degree divisors on the conics. The
[Euclidean formula](DET1092_EUCLIDEAN_BISECTION_FORMULA_2026-09-09.md)
already constructs their equations and polynomial maps. The
[halving-or-cycle theorem](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md)
already classified the specified old dependent split. Those are inputs,
not new results of this note.

**New specialization-obstruction deduction.** If (1) is1, the two distinct
smooth irreducible curves intersect in one geometric point, with
multiplicity1. Because both curves are defined over `Q`, that point is
fixed by Galois and rational. This also supplies a direct rational
basepoint on the conic, without odd-divisor descent.

Suppose its image `t_v` is a smooth original elliptic fibre. The branch
through the intersection is

\[
 Q_v=P_v(t_v),\qquad Q'_v=P_{w-v}(t_v),                     \tag{2}
\]

since the conic's trace is `P_w`. Hence both branches lie in the integral
specialized generic subgroup. The two roots may coincide at a ramified
value; the same inherited conclusion holds. In the actual14 applications
below they are distinct. The displayed18-section subgroup on the conic
base change acquires the explicit primitive specialization relation
`Q-P_v=0`; no assertion about its entire specialization kernel is needed.

This does not contradict generic rank18 on the conic: generic independence
and a relation at one rational base point are different statements.
Conversely, unit intersections need not account for every inherited
specialization, including odd-index saturation phenomena. Failure of this
criterion remains `UNKNOWN`.

## 2. A small, equation-only obstruction polynomial

The fixed trace word is

\[
 w=-e_2+e_{10}-e_{14},\qquad\lVert w\rVert^2=10.
\]

The complete bounded selector is all34 signed individual generators and
all544 signed sums on two distinct generators. It retains

\[
 \lVert v\rVert^2=4,\qquad\langle w,v\rangle=3.             \tag{3}
\]

These578 cheap lattice tests select14 words. Their exact conic
intersections give the following reduced parent parameters. Define the
primitive integral polynomial

\[
 \boxed{D(s)=\prod_{a/b\in\mathcal T}(bs-a),\qquad\deg D=14,} \tag{4}
\]

where every fraction is reduced and

```text
T = {
   8926758/9563789,        77335551/542052125,
  -3532573/2605690,         173994/239095,
      8412/5417,           7506249/16868462,
9448637025/44921394667,      -528/3635,
  19382781/296937400,     -62444544/466294325,
2663334983/381468770,         7999/19810,
   2397513/13274575,       4908471/12315172
}.
```

**Verified exact obstruction.** If `D(s)=0`, the original fibre is smooth,
this conic splits into two distinct rational points, and both belong to
the inherited generic subgroup by the words in the certificate. These
four properties are checked for every factor of (4). The complement of
this finite zero set receives no independence or incidence conclusion.

The small control is predicted particularly transparently by

\[
 v=-e_2+e_{15},\quad
 \lVert v\rVert^2=4,\quad\langle w,v\rangle=3.
\]

Its unique intersection has `s=-528/3635`. Formula (2) specializes to

\[
 Q=-P_2+P_{15},\qquad Q'=P_{10}-P_{14}-P_{15}.              \tag{5}
\]

This derives both the address and the dependence words from generic data;
it does not solve for a relation after importing a known split point.

## 3. Constructive extraction and independent proof

**Verified algorithmic application.** Use the Euclidean conic equation
and polynomial maps

\[
 W^2=q(t),\quad X=(b+hW)/2,\quad Y=-(hk+mW)/2.
\]

For each selected height-four section, `X_v,Y_v` are polynomials of degrees
at most4 and6. Choose `U h+V m=1` and form

\[
 L_v=U(2X_v-b)-V(2Y_v+hk).
\]

The finite intersection divisor is computed without a root search by

\[
 d_v=\gcd(2X_v-b-hL_v,\ 2Y_v+hk+mL_v,\ L_v^2-q).           \tag{6}
\]

For all14 selected words, (6) is linear. Its root gives `t_v`, and
`W_v=L_v(t_v)`. The fixed fractional-linear parent chart then gives `s_v`.
Equation (1) proves that this finite intersection is complete: there is
no second intersection hidden at infinity. The certificate verifies that
`q(t_v)` and the elliptic discriminant are nonzero in every case.

The producer constructs the conic afresh from generic sections. The
independent checker imports no producer: it uses a manual, reversed-order
rational-function group law for all section words and the negative trace,
verifies the Euclidean identities and full polynomial maps, checks (6),
and reconstructs the companion word `w-v` independently over `Q(t)` before
specializing it. All578 selected and rejected lattice exposures are saved.
The product (4) is verified primitive, squarefree and degree14.

The separate eleven-address comparison reads only the previous addresses,
after the generic prediction packet and independent proof have been sealed.
It matches the old dependent control. It does not match302, the productive
control `1926/2699`, or any of the eight unchanged null controls. Those ten
nonmatches mean only that they are outside this particular obstruction.
Their previously established split/nonsplit results are not reclassified.

## 4. What this explains, and what it does not

**New explanatory deduction.** A rational basepoint obtained by directly
intersecting the conic with an inherited unit-intersection section is
guaranteed to give an inherited specialization. Thus rationality proofs
and convenient parametrization basepoints can naturally supply perfectly
valid rational splits which are useless as seeds. This is a geometric
source of the dependent control, not a failure to search far enough or
evidence of a Sha obstruction. Moving elsewhere on that same conic can
produce independent specializations, as the old positive control proves.

The statement is deliberately about *direct unit intersections*. It does
not assert that every basepoint produced by an odd-divisor descent is
inherited, nor that the14-point list exhausts the dependent locus.

**Open.** This does not construct a positive302 seed or identify its useful
positive-genus carrier prospectively. The same conic remains nonsplit at302.
It does not explain amplification, add an exceptional direction, classify
the other fibres' full ranks, or authorize a full-atlas campaign.

## Reproduction

Construction took0.120 seconds internally. Both independent replay and
post-seal comparison completed within25-second process caps. No detached
job, point search, parameter scan, new CVP, class/unit computation, later
cascade point or production mutation was used.

- [Frozen generic-input protocol](../../artifacts/generated-results/elliptic-curves/det1092_unit_intersection_obstructions_v1/protocol.json)
- [All578 lattice exposures](../../artifacts/generated-results/elliptic-curves/det1092_unit_intersection_obstructions_v1/selection.json)
- [Fourteen exact obstructions](../../artifacts/generated-results/elliptic-curves/det1092_unit_intersection_obstructions_v1/summary.json)
- [Integral obstruction polynomial](../../artifacts/generated-results/elliptic-curves/det1092_unit_intersection_obstructions_v1/collision-polynomial.json)
- [Independent proof](../../artifacts/generated-results/elliptic-curves/det1092_unit_intersection_obstructions_v1/independent-replay.json)
- [Post-seal control comparison](../../artifacts/generated-results/elliptic-curves/det1092_unit_intersection_obstructions_v1/control-comparison.json)

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_unit_intersection_obstructions.sage check
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_unit_intersection_obstructions.sage compare
```

Checker SHA256: `ee5c6ddec6ee5ec8cbe5a456ecb7cedfbaf6e21087c1b9c69ee1fe7017a656bb`.
Replay SHA256: `cf52ff88d402488a4b9afa0f0365dfe3e8db7b6744f048ffffbf1427f8b916e7`.
The geometric argument is not formally verified or externally reviewed.
