# The q12/orbit5867 genus-one point factory

<!-- status-consumer: EC-K3-H3-Q12O5867-POINT-FACTORY 9399c93ee42ee2a4 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-TWO-PRIMARY-BOUNDARY 783482d8f700105d -->

## Outcome

The final exact neighbour now has a reusable arbitrary-point map on a common
affine open:

```text
(4A1/MW13 base s, point P) <-> (published R17 base t, point Q).
```

The implementation is
[`scripts/q12o5867_genus_one_point_factory.sage`](scripts/q12o5867_genus_one_point_factory.sage).
It composes the stored q12 chord pencil, the pointed binary quartic at the old
finite `I2` support `s=0`, and the certified raw-q12 to published-`t`
isomorphism. Both directions are literal rational formulas. All 42 points in
the five public complement bases pass exact parent equation, child equation,
and forward/inverse round-trip checks.

This closes the birational engineering gate.  The bounded specialized-MW13
coordinate route has now also reached its declared stopping boundary.  Three
exact targeted halvings kill the abstract index-eight quotient, but the
equation-level characteristic-zero sections certified within the boundary
still span index two.  Consequently the artifact continues to record exact
parent points and 128-bit canonical heights without inventing MW13 coefficient
vectors, and q12/orbit5867 is closed as the proposed rank-32 point-production
route.  The exact birational map remains a theorem.

## Exact formulas

Write the parent horizontal used by q12 as

```text
H=(X/Z^2,Y/Z^3)
```

and the stored Riemann--Roch basis as `(AA_i,BB_i)`, `i=0,1`. For a parent
point `P=(x,y)` over `s`, put

```text
m   = (y+y(H))/(x-x(H)),
r_i = AA_i(s)+BB_i(s) Z(s) m,
u   = -r_0/r_1.
```

If `AA=AA_0+u AA_1`, `BB=BB_0+u BB_1`, and the stored chord radicand is

```text
quartic(s,u) * square(s,u)^2,
```

then its ordinate is

```text
W = BB(s,u)^2 * (2x+x(H)-m^2) / square(s,u).
```

The exact identity `W^2=quartic(s,u)` is checked. For

```text
quartic=e+d s+c s^2+b s^3+a s^4,  e=q0^2,
```

the standard nonbranch pointed-quartic formulas give

```text
xg = (2 q0 (W+q0)+d s)/s^2,
yg = (4 q0^2(W+q0)+2 q0(d s+c s^2)-d^2 s^2/(2q0))/s^3.
```

The raw rootless point is

```text
x_raw = 9  (xg+b2/12),
y_raw = 27 (yg+(a1*xg+a3)/2).
```

Finally, with the pinned constants from the coordinate-match certificate,

```text
t     = (alpha*u+beta)/(gamma*u+delta),
x_raw = kappa^2 (gamma*u+delta)^4 x_R17,
y_raw = kappa^3 (gamma*u+delta)^6 y_R17.
```

The script also implements the inverse pointed-quartic and chord formulas.
The affine implementation rejects the zero section, pointed origin, infinite
base, and denominator-zero loci explicitly; those loci need the usual
projective companion charts and are not silently discarded.

## Five-control calibration

The parent has finite `I2` supports `0,r1,r2`. To remove the unhelpful raw
base scale, the artifact also records

```text
z = s*(r1-r2)/(r1*(s-r2)),
```

which sends `(0,r1,r2)` to `(0,1,infinity)`. The following ranges are for the
complete public complement at each control. `z bits` is the bit length of its
projective height; the height column is the numerical canonical height on the
specialized parent elliptic curve.

| published `t` | points | atlas-invisible | `z` bits | parent height |
| --- | ---: | ---: | ---: | ---: |
| `-2/377` | 8 | 3 | 46--104 | 193.668--437.741 |
| `-308/251` | 9 | 6 | 55--99 | 227.376--404.745 |
| `2456/135` | 10 | 8 | 55--133 | 214.686--631.967 |
| `-9529/5471` | 11 | 10 | 66--116 | 319.355--558.684 |
| `3/8` | 4 | 0 | 25--70 | 90.390--302.711 |

For the ten rank-28 directions invisible to the rational-bisection atlas:

| direction | public source point | `z` bits | parent height |
| --- | ---: | ---: | ---: |
| Q1 | 1 | 105 | 431.854 |
| Q3 | 3 | 95 | 442.216 |
| Q4 | 4 | 113 | 472.844 |
| Q5 | 7 | 116 | 558.684 |
| Q6 | 8 | 85 | 387.604 |
| Q7 | 9 | 77 | 352.165 |
| Q8 | 11 | 101 | 457.818 |
| Q9 | 15 | 91 | 421.328 |
| Q10 | 19 | 66 | 319.355 |
| Q11 | 22 | 78 | 355.455 |

Thus the first simplicity test is negative in the coordinates and invariants
currently certified: none of the ten missing rank-28 directions becomes a
small parent base or a low-height parent point. This is a finite exact
calibration, not a theorem that no small MW13 word exists. That stronger test
must wait for the missing equation-level saturated MW13 basis.

Because the requested trigger did not fire, no new small-parent-fibre search
or quotient-escape claim is promoted here. Enumerating parent fibres before
pinning the thirteen generators would conflate the full specialized group
with the generic MW13 specialization.

## Degree-one MW13 recovery gate

The proposed shortcut through inverse-parent-degree-one shell sections fails
exactly.  The complete marked physical shell contains 23 such classes, and
their MW coordinate vectors have rank 12, not 13.  This is independent of the
finite-field enumeration: the same obstruction is computed in the full
938-class lattice shell.  Exhaustive shells at `p=83,89,137`, with component
profiles, Abel traces, and smooth pairwise intersections, give compatible
reductions of that rank bound.

The exact parent Shioda height lattice is nevertheless recovered abstractly
as the Schur complement of the four `A1` roots.  Its rank is 13 and its
determinant is `237/4`, as required by determinant mutation from the
determinant-948 rootless frame.  The full physical `P.O=0` shell spans this
saturated lattice, but the current two-profile `p=89` modular shell spans an
index-8 subgroup.  Therefore a successful equation-level recovery must both
admit inverse degree greater than one and either enumerate the missing `I2`
component profiles or certify three independent exact 2-divisions.

The fail-closed replay is:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/recover_p1229_mw13_basis_qq.sage
```

It writes
[`../artifacts/generated-results/elkies-k3-p1229-mw13-degree1-recovery-gate.json`](../artifacts/generated-results/elkies-k3-p1229-mw13-degree1-recovery-gate.json)
and returns status
`REJECTED_EXACT_DEGREE1_SHELL_CANNOT_SPAN_MW13`.

Consequently no MW13 words are assigned to the 42 controls, the ten invisible
rank-28 directions remain inconclusive, and parent-point enumeration remains
disallowed.  This rejects the proposed recovery workflow; it does not prove
that no saturated equation-level MW13 basis exists by a widened shell.

## Exhaustive two-primary boundary

The bounded follow-up enumerates every one of the 16 equation incidence modes
at the four `I2` fibres.  The complete signed polynomial shells have 890, 926,
and 938 sections at `p=83,89,137`; every profile is nonempty at every prime.
The profile histograms are retained machine-readably.  Equation-mode bits are
not used as a global linear relabeling of the marked component pairings: exact
Q4 and recovered-half controls show that such a relabeling is invalid.  Full
NS names are instead certified by exact smooth section intersections against
pinned named anchors, with the reductions at all three primes retained as
independent fingerprints.

The strict unique part of the historical `p=89` shell has row HNF

```text
e0
e1+e10+e11+e12
e2+e11+e12
e3
e4+e10+e11
e5+e11+e12
e6+e10+e12
e7+e11+e12
e8+e10+e11+e12
e9+e12
2e10
2e11
2e12
```

and Smith factors `1^10,2,2,2`.  Thus the quotient is canonically
`(Z/2Z)^3`.  Minimum-`L1` independent doubled targets can be represented by
physical classes 146, 30, and 22, with quotient keys `001`, `010`, and `101`
and current-section words

```text
2*146 = 144 + 170 - 873
2*30  = 4 + 14 + 145
2*22  = 43 - 500 + 726.
```

For each word, the duplication quartic over `QQ(u)` factors as `1+3`.  Its
linear factor gives a literal polynomial section.  Exact curve substitution,
doubling, full-NS intersection naming, Shioda height, and three-prime profile
fingerprints give

| class | quotient key | equation profile | height |
| ---: | :---: | :---: | ---: |
| 146 | `001` | `0101` | `3` |
| 30 | `010` | `1101` | `5/2` |
| 22 | `101` | `1110` | `5/2` |

Adding these three abstract MW vectors reduces the Smith index from eight to
one.  This is not yet an equation-level saturated basis: the 24 exact sections
available after the declared lifts have rank 13, Smith factors
`1^12,2`, and height determinant `237`, whereas the saturated determinant is
`237/4`.  The remaining exact-pool coset is represented in the marked shell by
classes `21,48,216,516,543,701`.  Four corresponding regular `p=89` branches
do not rationally reconstruct (`216,701,543,48`); the two reconstructed
branches acquire the exact intersection names `21 -> 26` and `516 -> 521`.
No fourth halving or wider shell is permitted by the stopping rule.

Therefore no thirteen-section saturated equation basis and no words for the
42 public controls are promoted.  In particular the ten atlas-invisible
rank-28 directions cannot be tested for short MW13 words.  They remain large
and unstructured under the exact invariants already available, so no bounded
parent-MW enumeration is authorized.  This closes q12/orbit5867 as the
rank-32 point-production route while retaining the arbitrary-point birational
theorem.

Replay the exact boundary with:

```bash
for p in 83 89 137; do
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/construct_h92_q12o5867_p0_shell_modp.sage \
    --prime "$p" --all-component-profiles --include-all-records \
    --output "artifacts/local/elkies-k3/q12o5867-p0-shell-all-profiles-mod${p}.json"
done

~/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/identify_h92_q12o5867_two_primary_cosets_mod89.sage
~/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/attempt_h92_q12o5867_three_target_halvings_qq.sage
~/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_h92_q12o5867_target_support_cross_prime.sage
~/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/verify_h92_q12o5867_three_halvings_and_basis_qq.sage
```

The final certificate is
[`../artifacts/generated-results/elkies-k3-q12o5867-two-primary-boundary.json`](../artifacts/generated-results/elkies-k3-q12o5867-two-primary-boundary.json).

## Replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/q12o5867_genus_one_point_factory.sage \
  --mode controls \
  --output artifacts/generated-results/elkies-k3-q12o5867-genus-one-point-factory-controls.json
```

Expected terminal status:

```text
Q12O5867POINTFACTORY|controls=42|roundtrips=exact|status=PASS_EXACT_Q12O5867_BIRATIONAL_POINT_MAP_AND_CONTROL_ROUNDTRIPS
```

The artifact is
[`../artifacts/generated-results/elkies-k3-q12o5867-genus-one-point-factory-controls.json`](../artifacts/generated-results/elkies-k3-q12o5867-genus-one-point-factory-controls.json).
Its exact point coordinates are intentionally kept machine-readable rather
than copied into this note.
