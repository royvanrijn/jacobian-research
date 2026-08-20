# E6 MW3 progress — 2026-08-20

This note records the current frontier of the alternate E6 MW-rank-3 fibration attack.

## Lattice/fibration target

The selected MW3 neighbor has

    ADE = E6 + A3^2 + A1^2
    root determinant = 192
    MW rank = 3

and exact reduced Mordell-Weil height Gram

    (1/12) * [ 23 -10  -8 ]
             [-10  23   1 ]
             [ -8   1  23 ]

with determinant 79/16.

The corresponding preferred Kodaira configuration is

    IV* + I4 + I4 + I2 + I2 + 4 I1.

We normalize IV* at infinity, I4 fibers at 0 and 1, and the two I2 fibers at lambda and mu. The IV* condition gives deg A <= 5 and deg B <= 8 in short Weierstrass form.

## Canonical P1 component data

From the component-label enumeration one canonical triple is

    P1 = (1,0,1,0,0), P1.O = 0
    P2 = (1,1,2,1,1), P2.O = 1
    P3 = (2,3,0,0,0), P3.O = 0

with all pairwise intersections equal to 2. P1 and P3 are polynomial sections in this normalization.

## High-coefficient triangular chain

For P1, the section-first ansatz simplifies to

    X = (s1-x0-x1)t^2 + x1*t + x0
    Y = t^4 + (-y0-y1-y2-1)t^3 + y2*t^2 + y1*t + y0.

The four highest section coefficients eliminate successively:

    P1_7 -> b4
    P1_6 -> b5
    P1_5 -> y2
    P1_4 -> a3.

The observed expressions remain comparatively small through these four steps. The resulting system has 13 variables and 12 equations: eight fiber equations plus P1_0..P1_3.

## Exact P1_0 parametrization

The constant section equation factors as

    P1_0 = y0^2 - (x0-s0)^2*(x0+2*s0).

It has the polynomial parametrization

    x0 = r0^2 - 2*s0,
    y0 = r0*(r0^2 - 3*s0).

Substitution kills P1_0 identically and replaces x0,y0 by one parameter r0. This leaves 12 variables and 11 displayed equations.

After the parametrization, each of P1_1,P1_2,P1_3 is only a seven-term polynomial and is affine-linear in a1, x1 and y1.

## Hidden dependency among P1_1,P1_2,P1_3

Solving P1_1 for y1 gives, on the generic branch,

    y1 = (a1 + 3*(r0^2-s0)*x1)/(2*r0).

After this exact substitution,

    P1_2 = 0
    P1_3 = 0

identically.

This is an important structural simplification: the final two low section coefficients are dependencies, not independent conditions. Thus P1 contributes only the expected conditions, and the genuine reduced locus consists of

    11 variables
    8 independent fiber equations
    expected dimension 3.

### Correction to an exploratory file

A quick temporary one-liner previously formed `/tmp/e6-param1.ms` by taking only the numerator of the rational y1 solution. Because the denominator is proportional to r0, discarding it is incorrect. Any solver result based on that temporary file should be ignored.

The committed scripts avoid this error. The coordinate-slice builder fixes r0 to a nonzero field element first, so the y1 denominator becomes a nonzero scalar; it then performs the exact substitution.

## Groebner/solver observations

The experiments also established several practical points:

1. msolve itself is functioning correctly. Earlier SIGSEGVs came from an exporter bug that wrote literal `\\n` characters, producing a malformed one-line input file.
2. Dense three-slice systems on the pre-parametrized representation can return `No solution` quickly because dependent section equations were still being counted as constraints.
3. A one-slice 13x13 system over GF(101) reaches large F4 matrices around degrees 18-19. Repeating over GF(31) gives nearly the same progression, showing the complexity is structural rather than prime-specific.
4. Dense affine elimination of the three final slice equations is computationally poor: substituting 9-term linear forms into degree-20+ fiber equations causes severe multinomial expansion.
5. Coordinate slices are much cheaper. Replacing selected variables by constants keeps the final 8x8 systems small enough for bounded parallel msolve exploration.

## Current search strategy

The current pipeline is:

1. Build the unsliced four-elimination E6/P1 system.
2. Parametrize P1_0 exactly.
3. For each search seed, choose a coordinate slice that includes r0 plus two other coordinates.
4. Set r0 to a nonzero constant.
5. Solve P1_1 -> y1 exactly after slicing.
6. Verify P1_2 and P1_3 vanish identically.
7. Export only the eight independent fiber equations in the eight remaining variables.
8. Run msolve with a bounded timeout.
9. Search many seeds/triples in parallel rather than spending unbounded time on one Groebner basis.

Recommended first run:

    sage elkies-k3/scripts/export_e6_p1_sliced.sage \
      --p 31 --seed 1 --slices 0 \
      --out artifacts/local/elkies-k3/e6-base.ms

    sage elkies-k3/scripts/parametrize_e6_p1_at_0.sage \
      --input artifacts/local/elkies-k3/e6-base.ms \
      --out artifacts/local/elkies-k3/e6-param0.ms

    python3 elkies-k3/scripts/run_e6_coordinate_slice_search.py \
      --input artifacts/local/elkies-k3/e6-param0.ms \
      --workers 8 --threads 2 --timeout 45 --seeds 32

A successful finite-field slice will give a concrete surface on the P1 locus. The next step is then to reconstruct its eliminated coefficients and impose P3, which should cut the dimension again toward the desired rank-3 K3 family.
