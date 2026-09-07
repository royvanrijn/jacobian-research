# Curve 302: bounded maximal fixed-resolvent orbit slices

This is a narrowly bounded constructor experiment, not a class-group or rank
result.  The active bounded protocol is
[v3](../cas/CURVE302_FIXED_RESOLVENT_ORBIT_PROTOCOL_V3.json), with its
[v4 replay](../cas/curve302_fixed_resolvent_orbits_v4.sage) correcting the
Sage-10.9 API in the otherwise dormant squareclass-extraction branch.  Its
endpoint is still a new strict character outside M24 followed by a sealed
2-cover and `PointsQI` run.

## Mandatory maximal-order gate

The displayed descent cubic has order index

```text
9923258571576383661788160
```

in its maximal cubic order.  It was therefore never used as a fixed
resolvent form.  The independently recomputed normal maximal-order form is

\[
-9923258571576383661788160X^3
-7851180678455372984811918483571145X^2Y
+869829586504686453875887554455898189850XY^2
+76467148428008343737496240497823633059162876309Y^3,
\]

whose discriminant is the certified field discriminant of the cubic field.
The preflight is recorded in the
[v3 result](../../artifacts/generated-results/elliptic-curves/curve302_fixed_resolvent_orbits_v3.json)
and independently replayed in its
[verification](../../artifacts/generated-results/elliptic-curves/curve302_fixed_resolvent_orbits_v3_verification_v1.json).

## Bounded result

The original 49-pair `M11=0` affine chart was entirely the split
\(\mathbf Q\times K\) locus.  A corrected direct perturbation slice then
made 3969 declared determinant attempts and retained 56 distinct primitive
integral pairs satisfying \(4\det(xA-yB)=F_{\max}\) exactly.  Each of the 56
reconstructed associative quartic orders has an explicit nonzero zero
divisor (the same normalized basis element).  Thus none is a quartic field;
none reaches the maximality, `S4`, real-splitting, local-splitting,
squareclass, M24-independence, cover, or point-search gate.

The v2 artifact is retained as a software-failure record: its direct quotient
used \(F_{\max}\) instead of \(F_{\max}/4\), so it solved no pairs.  The v3
protocol and result correct precisely that scale and no mathematical outcome
is inferred from v2.

Consequently the desired chain remains open:

\[
\text{fixed-resolvent quartic field}
\longrightarrow\text{ new strict class outside }M_{24}
\longrightarrow2\text{-cover}\longrightarrow\text{PointsQI}.
\]

No cover job was created and no rank statement follows from these bounded
null slices.

## Replay

From `research/`:

```sh
sage -python elliptic-curves/cas/curve302_fixed_resolvent_orbits.sage check
sage -python elliptic-curves/cas/curve302_fixed_resolvent_orbits_v3.sage check
sage -python elliptic-curves/cas/verify_curve302_fixed_resolvent_orbits_v3.sage check
sage -python elliptic-curves/cas/curve302_fixed_resolvent_orbits_v4.sage check
```
