# Curve 302: retrospective residual visibility geometry

This is a post-search diagnostic on the sealed public rank-31 group. It does
not select a prospective centre, construct a new point, identify a canonical
splitting of `D/M24`, or give a rank upper bound. Mathematical status remains
in [`MATH_STATUS.json`](../../MATH_STATUS.json).

## Result

The anticipated uniform split does **not** occur in the declared diagnostic.
The four recovered locally visible directions have very small finite reduced
coordinates (and two have a projective-infinity representation), whereas the
strict sector is heterogeneous:

| deterministic quotient panel | finite coordinate decimal digits, by direction |
| --- | --- |
| recovered local (4) | `1, 1, 5, 2` |
| recovered strict (3) | `17, 56, 144` |
| residual strict (7) | `176, 67, 152, 147, 165, 4, 88` |

The small residual strict entry has finite height `2945`; the largest residual
entry has 176 decimal digits. Thus the actual result supports a strong
local-versus-strict contrast for the four locally visible directions, but it
does **not** support the proposed statement that all seven remaining strict
directions are uniformly invisible. The individual seven-direction bases are
deterministic integral complements, not canonical arithmetic directions, so
their individual coordinate minima should not be promoted to a quotient
invariant.

The exact data are in
[`curve302_residual_visibility_geometry_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_residual_visibility_geometry_v1.json).

## The common `M/2M` label space

There is an exact, useful correspondence behind the `2^17` equality. For the
full geometric parent lattice `M=M17`, degree-two divisor translation is

\[
 w\longmapsto w+2x,
\]

and a pointed-quartic centre changes by the same action

\[
 Q\longmapsto Q+2x.
\]

Consequently both are labelled by the same finite torsor `M/2M`, of size
`131072`. The retained LLL changes differ, so the audit records the explicit
basis conversion: a chart label `p` in the certified parent basis has
degree-two orbit mask `p B^{-1} mod 2`, where `B` is the degree-two row LLL
change to that basis.

The `40917` rational bisections are therefore exactly the norm-10 shell in
this shared label space. They are not automatically the deep empty-ball
centres chosen by the search. In the 448 exact chart rows here, 150 labels
are in that rational shell, 239 are genus-one lattice candidates, and 59 lie
outside the section-nonnegative degree-two survivors. This is label-level
geometry only: it neither constructs every bisection equation nor associates
a split special fibre to a label.

## Declared diagnostic and boundary

For each of a deterministic `4 local + 3 strict + 7 residual strict` direct
summand panel, the audit:

1. scans every `131071` nonzero `M17/2M17` class by deterministic Babai
   reduction in a 384-bit PARI specialized-height Gram rounded at `10^6`;
2. re-solves the best fixed 512 parity labels with DD and MPFR CVP, checking
   agreement in the integral rounded metric;
3. constructs and verifies the exact pointed quartic and translated target
   for the best 32 CVP labels; and
4. records its exact reduced projective coordinate, quartic square identity,
   and degree-two `M/2M` orbit label.

For a parity `p`, take the reduced centre `q_p=p+2u` and choose `p+2z`
closest to `2R`. The target represented in the chart is `R+(u-z)`, since

\[
 2(R+(u-z))-q_p=2R-(p+2z).
\]

This is the requested translation/CVP relation. Exact coordinates are minimum
only within the 32-chart vetted set. Because the final rational coordinate has
chart-dependent PGL2 distortion, the result does not claim a global
reduced-coordinate minimum across all `131072` labels. Projective infinity is
reported separately rather than assigned an artificial finite height.

## Replay

From `research/`:

```sh
sage -python elliptic-curves/cas/audit_curve302_residual_visibility_geometry.sage --check
python3 -m unittest elliptic-curves.tests.test_curve302_residual_visibility_geometry
```
