# Three parents outside the retained Kihara path

The fixed parent ratios **p/q=1,2,3** give three additional K3 surfaces
over Q. Independent finite-field point counts distinguish each from every
previously retained parent and from one another. The proved tested
portfolio now contains **at least fourteen Q-distinct parents**.

Direct section equations and exact canonical-height matrices certify
rational generic section spans of ranks **7,9,12**, respectively. These
are lower bounds for the full generic groups. Full generic ranks and NS
lattices remain unknown. No specialized point-search boxes or new
near-record curves are part of this intake, and literature novelty is
not asserted.

Authority: `EC-KIHARA-THREE-POSITIVE-PARENTS-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).
The [expansion certificate](../../artifacts/generated-results/elliptic-curves/kihara_positive_parent_expansion_v1.json)
contains the models, sections, exact height matrices, all39 separation
witnesses and source bindings.

## Fixed population and construction

The [preceding coverage proof](KIHARA_QUADRATIC_SECTION_AND_PARENT_COVERAGE_2026-09-07.md)
shows that the existing rational rank14 path has only-1/2<p/q<0.
Its unrestricted six-centre formula is defined at positive ratios too.
Before arithmetic, this experiment fixes the first three positive integer
ratios1,2,3, with no replacement, ranking or score selection. A candidate
falling outside the usual20I1+I4 configuration is retained for geometric
classification, rather than counted as a failed rank search.

Set q=1 and normalize the six centres by their second entry. The degree12
product of the shifted roots admits the same square approximation, and
the remainder divided by T^2 is a quartic. Its twelve root points and the
thirteenth rational point from the unrestricted formula are checked
exactly. Taking the first root point as origin gives twelve Jacobian
points. The involution point K satisfies

```
6K = S1+...+S11.
```

Replace S1 by K before height calculations. This preserves their rational
span; an exact subgroup index is not inferred at these new special
parameters. In particular, listing twelve points does not prove rank12.

## Geometry changes at the small positive ratios

The exact minimal coefficient degrees are8 and12 in all three cases.
The discriminant is coprime to c4, so its finite repeated factors give
semistable fibres. All the I2 base places are rational. Infinity has a
rational split node.

| p/q | Finite fibres | Infinity | Certified section-span rank | Independent height determinant |
|---:|---|---|---:|---:|
|1|8I1+6I2|split I4|7|27/8|
|2|14I1+2I2|split I6|9|15/2|
|3|20I1|split I4|12|189|

All three Euler numbers are24. The determinants belong to the displayed
independent subgroups, not an asserted full MW or NS lattice.

The first producer left ratios1 and2 outside its fixed squarefree
discriminant check. The independent extension completed ratio1 and then
stopped at ratio2: the latter's discriminant degree is18, making infinity
I6 instead of I4. That completed first-parent checkpoint is preserved.
The corrected extension checks only ratios2 and3. It uses multiplier6
for I6 and4 for I4, each killing every finite I2 component as well.

For a nonzero multiple mP on identity components, the independent height
calculation uses

```
h(P) = (4 + max(deg denominator(x(mP)),
                deg numerator(x(mP))-4)) / m^2.
```

This avoids the producer's local node/component implementation. Every
pairing satisfies the exact parallelogram identity. The resulting Gram
rank gives the stated section-span rank; the selected principal minor
is positive definite. Zero multiples are handled separately, without
assuming absence of torsion on the exceptional parents.

## Parent separation

The independent verifier counts every smooth fibre with elliptic-curve
cardinality arithmetic. For nodal fibres it uses the exact split or
nonsplit node count and the resolution correction: p for I2, and the
split In cycle at infinity contributes n*p. Coefficient degrees, fibre
multiplicities and the smooth-reduction conditions are checked at every
retained prime.

| p/q | #X(F131) | #X(F239) | #X(F251) |
|---:|---:|---:|---:|
|1|19466|60938|67466|
|2|19340|60954|67475|
|3|19178|60929|67236|

Each new row differs at a common good prime from all twelve prior retained
presentations: five Kihara controls, six Mestre parents and X948. The new
rows also differ pairwise. Those39 comparisons prove Q-nonisomorphism.
The old twelve have only an eleven-parent lower bound under their retained
comparison certificates; adding these three therefore gives **at least14**.
The prior counts are reused from their already independently verified
certificates. No geometric nonisomorphism or new lattice type follows
solely from these Q-point counts.

## Replay, cost and next use

The [portable input](../../artifacts/generated-results/elliptic-curves/kihara_positive_parent_bundle_v1.json)
and [current standalone verifier](../cas/verify_kihara_positive_parents_v2.sage)
reconstruct the product, quartic, Jacobian, section map, heights and new
surface counts without repository imports. In an empty output directory:

```sh
sage -python verify_kihara_positive_parents_v2.sage \
  --input kihara_positive_parent_bundle_v1.json --output-directory OUTPUT
python3 elliptic-curves/cas/report_kihara_positive_parents.py --check
```

The retained execution used one worker with2GiB and a120-second cap per
stage. Total supervised cost is **29.316635372 seconds**, including a
0.640388375-second source-loading failure and the2.328151256-second
partially completed first standalone attempt. The remaining two-parent
replay took24.476874981 seconds. Completed arithmetic was not repeated
after the I6 correction.

These parents can now enter a separately frozen fibre comparison using
their own certified section subgroups. Generic rank alone does not predict
their exceptional-fibre yield. The [subsequent fixed pilot](KIHARA_POSITIVE_PARENT_POINT_PILOT_2026-09-07.md)
now completes294 boxes on six fibres, with three new directions and
strongest lower bound14. Its exact specialized spans and torsion handling
are recorded separately from this generic intake.
