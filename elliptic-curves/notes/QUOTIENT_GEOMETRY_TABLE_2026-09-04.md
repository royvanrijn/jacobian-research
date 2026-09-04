# Specialization quotient geometry and the detector cutoff test

The complete cross-experiment table is
[`quotient_geometry_table_v1.json`](../../artifacts/generated-results/elliptic-curves/quotient_geometry_table_v1.json).
It contains five usable known R17 controls, all sixteen refreshed R17 ladder
fibres, and all nine A1/MW16 parent presentations: 30 presentation rows and
230 exact independent blind-recovery events.

The answer is mixed but decisive.  Small intrinsic quotient height is useful
fibre geometry, but it does not order first discovery.  All thirteen
exact-containment presentations with a nonempty, strictly partial initial
recovery fail the successive-minimum-prefix test.  Seven remain failures after
the adaptive/final stage.  The other three containment comparisons are
deliberately null: the existing verification does not place the full blind
subgroup integrally in the displayed subgroup for R17 curves 478 and 539, and
the curve-542 specialized MW16 subgroup is half-integral in the public basis.

The pointwise completing-the-square experiment is also recorded.  For all 230
recovered directions, the projection coefficients `a`, Schur complement
`lambda`, and nearest half-lattice phase are recomputed.  The nearest centre is
unchanged when the canonical-height Gram is rounded independently at scales
`10^5` and `10^6`.  For the 190 initial-chart events, the actual source phase
exceeds the pointwise best half-lattice phase by a median `9.33804`.  Thus even
the charts that found the points were generally not their nearest possible
midpoints; reduced-coordinate distortion is a separate material effect.

## What is recorded

For a specialized generic subgroup `M` inside each certified displayed
subgroup, the artifact records the full Schur-complement Neron--Tate Gram on
the displayed quotient, its determinant (labelled `regulator`), a complete
Fincke--Pohst enumeration sufficient to obtain the successive minima, and one
deterministic rank-increasing witness for every minimum.  Here

\[
  j_{\rm displayed}=\operatorname{rank}(G_{\rm displayed})-
  \operatorname{rank}(M);
\]

it is not the elliptic `j`-invariant.

Every exact independent recovery event has its intrinsic quotient energy
`q_t(P)`, first recovery stage, chart and half-lattice centre, phase, reduced
coordinate, and coordinate distortion.  If `Q` is the chart centre, the
replayed decomposition is

\[
 \log H(s(P))
 =2q_{t,Q}(P)+2\phi_{t,Q}(P)+\delta_s(P),
 \qquad q_{t,Q}(P)=\widehat h_{/M}(P-Q/2),
\]

with `q_{t,Q}(P)=q_t(P)` on the initial charts.  Hence the actual source-chart
window at search bound `B=100000` is

\[
 q_{t,Q}(P)\le {\log B-\delta_s(P)\over2}-\phi_{t,Q}(P).
\]

This identity holds for every one of the 230 recovered witnesses.  That is a
replay of visibility in the chart that found the point, not a prospective
assignment of phases and distortions to unseen directions.

For a single recovered point `R`, let `G` be the specialized generic height
Gram, `b_i=<P_i,R>`, and `a=G^{-1}b`.  The table separately checks

\[
 \lambda=\widehat h(R)-b^TG^{-1}b,
 \qquad
 \min_{m\in\mathbf Z^r}\widehat h\!\left(R-\frac12\sum_i m_iP_i\right)
 =\lambda+\operatorname{dist}_G\!\left(a,\frac12\mathbf Z^r\right)^2.
\]

The first term is intrinsic new height.  The distance is the best positional
phase available anywhere in the inherited half-lattice.  The actual-source
phase records the chart that first found `R`; for adaptive charts its centre
has a nonzero original quotient component, so the source-versus-optimum
comparison is intentionally left null.  Finally
`delta_s(R)=log H(s(R))-hhat(2R-Q)/2` records coordinate distortion.

## Human index of the complete table

`Gram` gives the dimension of the full matrix stored in the linked JSON.
Minima below are rounded only for this index; the artifact retains the full
PARI output.  `I`, `A`, `C43`, and `MWmax` mean initial, adaptive,
generic-deepest-43 control, and complete maximum-depth MW16 recovery stages.
The last column is equality with the deterministic successive-minimum prefix
at the final recovered rank; `?` is the fail-closed containment case.

| presentation | `j_displayed` | Gram | quotient regulator | successive minima | blind stages | prefix? |
|---|---:|---:|---:|---|---:|:---:|
| R17 control rank 21 | 4 | 4x4 | 3.97067e2 | 3.194, 5.655, 5.703, 5.715 | C43:4 | yes |
| R17 control rank 25 | 8 | 8x8 | 7.47287e6 | 8.345, 8.876, 9.081, 9.154, 10.19, 10.68, 10.74, 11.22 | C43:8 | yes |
| R17 control rank 26 | 9 | 9x9 | 1.92911e8 | 11.11, 11.49, 11.84, 12.06, 12.17, 12.68, 13.26, 13.54, 13.65 | C43:9 | yes |
| R17 control rank 27 | 10 | 10x10 | 2.20104e10 | 13.62, 14.48, 14.52, 16.06, 16.83, 17.34, 17.46, 17.60, 17.62, 17.72 | C43:10 | yes |
| R17 control rank 28 | 11 | 11x11 | 1.08767e12 | 16.53, 17.61, 18.87, 18.97, 19.21, 20.09, 20.17, 20.62, 20.70, 20.74, 21.95 | C43:9 | **no** |
| R17 ladder curve 478 | 4 | 4x4 | 4.66612e3 | 8.315, 9.671, 10.52, 10.53 | I:6 | ? |
| R17 ladder curve 498 | 6 | 6x6 | 7.81874e4 | 7.401, 7.757, 7.816, 8.485, 8.529, 9.598 | I:6 | yes |
| R17 ladder curve 531 | 11 | 11x11 | 4.59155e11 | 15.76, 16.01, 17.26, 17.41, 17.70, 18.02, 18.39, 18.67, 18.75, 19.48, 19.87 | I:9+A:2 | yes |
| R17 ladder curve 532 | 3 | 3x3 | 6.37954e1 | 4.380, 4.550, 5.334 | I:3 | yes |
| R17 ladder curve 534 | 11 | 11x11 | 5.98825e12 | 16.73, 18.02, 22.59, 22.82, 23.37, 23.52, 23.55, 23.73, 23.88, 24.18, 25.22 | I:4+A:7 | yes |
| R17 ladder curve 535 | 11 | 11x11 | 4.94543e12 | 18.42, 18.53, 20.50, 21.31, 21.39, 21.68, 22.41, 22.65, 23.12, 23.82, 23.92 | I:3+A:7 | **no** |
| R17 ladder curve 536 | 11 | 11x11 | 5.55650e11 | 13.96, 15.53, 17.61, 17.97, 18.17, 18.43, 18.52, 19.38, 20.15, 20.24, 20.95 | I:4+A:7 | yes |
| R17 ladder curve 537 | 10 | 10x10 | 7.39570e10 | 16.35, 16.49, 17.64, 17.76, 18.09, 18.20, 18.32, 18.78, 19.11, 19.86 | I:7+A:3 | yes |
| R17 ladder curve 538 | 5 | 5x5 | 9.90793e3 | 7.427, 8.053, 8.128, 8.182, 8.649 | I:5 | yes |
| R17 ladder curve 539 | 6 | 6x6 | 1.44261e5 | 6.808, 7.146, 8.059, 8.324, 9.324, 12.44 | I:6 | ? |
| R17 ladder curve 540 | 8 | 8x8 | 4.67914e7 | 11.43, 12.00, 12.23, 12.26, 12.78, 13.73, 14.06, 14.10 | I:8 | yes |
| R17 ladder curve 541 | 8 | 8x8 | 1.05896e7 | 8.735, 9.211, 9.447, 10.34, 11.06, 11.07, 11.42, 12.36 | I:8 | yes |
| R17 ladder curve 543 | 12 | 12x12 | 2.35304e13 | 19.35, 19.93, 20.13, 20.44, 20.90, 21.00, 21.18, 21.19, 21.43, 22.03, 22.34, 22.57 | I:3+A:9 | yes |
| R17 ladder curve 544 | 11 | 11x11 | 1.63692e13 | 21.64, 21.93, 24.45, 25.08, 25.19, 25.87, 25.95, 26.22, 26.95, 28.11, 28.42 | none | yes (empty) |
| R17 ladder curve 545 | 11 | 11x11 | 1.69594e12 | 15.39, 17.38, 17.84, 18.78, 19.09, 19.88, 19.93, 20.01, 20.15, 23.67, 24.65 | I:6+A:5 | yes |
| R17 ladder curve 546 | 8 | 8x8 | 7.98728e6 | 7.785, 8.309, 8.752, 9.210, 9.876, 9.939, 11.23, 11.86 | I:8 | yes |
| MW16 curve 398, parent 16875 | 14 | 14x14 | 3.08604e13 | 13.05, 13.93, 14.76, 14.94, 15.57, 15.92, 16.06, 16.10, 16.34, 16.50, 17.10, 17.11, 17.27, 17.68 | MWmax:5 | **no** |
| MW16 curve 398, parent 63669 | 14 | 14x14 | 3.08604e13 | 13.05, 13.93, 14.76, 14.94, 15.57, 15.92, 16.06, 16.10, 16.34, 16.50, 17.10, 17.11, 17.27, 17.68 | MWmax:5 | **no** |
| MW16 curve 400, parent 53042 | 12 | 12x12 | 2.65082e11 | 13.67, 14.00, 14.51, 14.66, 14.68, 14.74, 15.45, 15.56, 15.65, 16.23, 16.46, 16.66 | MWmax:5 | **no** |
| MW16 curve 400, parent 62992 | 12 | 12x12 | 2.65082e11 | 13.67, 14.00, 14.51, 14.66, 14.68, 14.74, 15.45, 15.56, 15.65, 16.23, 16.46, 16.66 | MWmax:5 | **no** |
| MW16 curve 401, parent 57487 | 11 | 11x11 | 3.89276e10 | 10.21, 13.59, 13.65, 13.85, 14.06, 14.59, 14.62, 14.76, 15.32, 15.41, 16.23 | MWmax:10 | **no** |
| MW16 curve 542, parent 30486 | 10 | 10x10 | 1.52375e8 | 7.146, 7.719, 7.802, 8.395, 9.239, 9.273, 9.535, 9.777, 10.59, 11.04 | MWmax:10 | ? |
| MW16 curve 548, parent 31627 | 8 | 8x8 | 1.22449e6 | 4.998, 6.232, 8.244, 8.251, 8.348, 8.566, 8.628, 9.043 | MWmax:8 | yes |
| MW16 curve 548, parent 54835 | 8 | 8x8 | 1.22449e6 | 4.998, 6.232, 8.244, 8.251, 8.348, 8.566, 8.628, 9.043 | MWmax:8 | yes |
| MW16 curve 548, parent 63647 | 8 | 8x8 | 1.22449e6 | 4.998, 6.232, 8.244, 8.251, 8.348, 8.566, 8.628, 9.043 | MWmax:8 | yes |

## The precise detector question

The necessary scalar-window test is:

> At every frozen detector stage, is its exact recovered rational quotient
> subspace the deterministic successive-minimum prefix of the same rank?
> Equivalently away from boundary ties, can the recovered subspace be exactly
> the span of the quotient directions below one intrinsic `q_t` cutoff?

No.  The rank-28 R17 control already recovers the fourth successive-minimum
witness while missing the third, and curve 535 recovers later minimum
directions while missing the fourth.  Five MW16 presentation rows also fail.
Thus quotient height is useful geometry but not a sufficient detector state:
the half-lattice phase and coordinate distortion select directions that are
not ordered by intrinsic `q_t` alone.

Curve 398 is the clearest separation.  Its displayed quotient has rank 14 and
successive minima `13.0465` through `17.6782`.  The first five recovered
directions have `lambda` from `15.5748` through `23.0175`; their subspace
contains only the fifth and eighth deterministic successive-minimum witnesses
and misses the first four.  Their actual source phases are `9.3671` through
`18.7229`, compared with pointwise optimal half-lattice phases `4.6039` through
`6.5612`, while their coordinate distortions are strongly negative,
`-63.5878` through `-51.5398`.  On this diagnostic, curve 398 is not explained
by a detector simply harvesting the intrinsically shortest new directions, nor
by selecting each point's nearest midpoint.  Chart-specific coordinate
compression is doing essential work.

Intrinsic height is not irrelevant at the fibre-selection scale.  The
zero-recovery ladder control, curve 544, has the table's highest quotient
minimum profile (`21.6425` through `28.4242`).  But nearby profiles do not give
a cutoff: curves 543 and 534 are fully recovered with successive minima up to
`22.5682` and `25.2237`.  The finite panel therefore supports using quotient
height as one family feature, not as a chart schedule or promotion rule.

This does not show that the richer chart-dependent mechanism is solved.  It
shows exactly why a scalar successive-minimum cutoff cannot solve it.  A
prospective version must predict phase and distortion for candidate charts
without seeing the hidden point.

## Replay and claim boundary

Run:

```sh
sage -python elliptic-curves/cas/build_quotient_geometry_table.sage --check
python3 -m unittest elliptic-curves/tests/test_quotient_geometry_table.py
```

Canonical heights, Grams, determinants, projection coefficients, CVP phases,
and distortion terms are 80-digit numerical diagnostics, not interval
certificates.  The two-scale CVP agreement is a stability check, not a rigorous
height bound.  Exact group identities, rational embeddings, Smith quotient
maps, and rational-subspace comparisons are replayed exactly.  Each quotient
lives in a certified displayed subgroup, which need not be the full
Mordell--Weil group.  The nine MW16 rows are nested in five target curves and
are not nine independent observations.  Original search time, chart, and
height bounds still apply to every miss.
