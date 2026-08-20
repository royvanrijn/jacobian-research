# E6 MW3 direct attack

Node:

    ADE = E6 + A3^2 + A1^2
    MW rank = 3

Exact reduced MW height Gram:

    (1/12) * [[23,-10,-8],[-10,23,1],[-8,1,23]]

Preferred Kodaira model:

    IV* + I4 + I4 + I2 + I2 + 4 I1

The IV* fiber is placed at infinity. In short Weierstrass form this directly forces

    deg A <= 5
    deg B <= 8

which is substantially simpler than the A10/I11 branch.

## Component/section data

The canonical component-label triple used for the first section attack is:

    P1 = (1,0,1,0,0), P1.O = 0
    P2 = (1,1,2,1,1), P2.O = 1
    P3 = (2,3,0,0,0), P3.O = 0

with all pairwise section intersections equal to 2. Thus P1 and P3 can be represented polynomially while P2 is the one-denominator section.

## P1 triangular reduction

The section-first P1 model has a strong exact triangular chain. Starting with 17 variables, the high section coefficients solve successively as

    P1_7 -> b4
    P1_6 -> b5
    P1_5 -> y2
    P1_4 -> a3

leaving 13 variables and 12 equations before using the low coefficients P1_0..P1_3.

The low constant equation factors as

    P1_0 = y0^2 - (x0-s0)^2 (x0+2*s0).

Use the rational parametrization

    x0 = r0^2 - 2*s0
    y0 = r0*(r0^2 - 3*s0).

This kills P1_0 identically and replaces x0,y0 by r0, giving 12 variables / 11 displayed equations.

After this parametrization P1_1 is affine-linear in y1. On the generic open branch its exact solution is

    y1 = (a1 + 3*(r0^2-s0)*x1) / (2*r0).

After substituting that exact relation, P1_2 and P1_3 vanish identically. Therefore they are dependent equations, not two additional constraints. The genuine P1 locus is consequently represented by 8 fiber equations in 11 variables, hence expected dimension 3.

Important correction: an exploratory /tmp one-liner briefly formed `e6-param1.ms` using only the numerator of this rational y1 expression. That file is invalid and must not be used. The committed coordinate-slice builder performs the exact elimination correctly.

## Solver lessons

Several useful negative experiments clarified the correct representation:

- Dense three-plane slicing before removing the dependent P1 equations gave instant `No solution` results after the msolve file newline bug was fixed.
- A single dense slice on the 13-variable system produced very large F4 matrices around degrees 18-19, essentially independent of using p=101 or p=31. The difficulty is structural, not a bad-prime effect.
- Exact P1_0 parametrization exposed the dependency P1_2=P1_3=0 after P1_1 elimination.
- Dense affine substitution of the three final slicing equations causes severe multinomial expansion in the high-degree fiber equations.
- Coordinate slices are dramatically cheaper. Fixing three coordinates to constants keeps the substituted fiber equations compact enough for practical msolve probes.

## Current reproducible pipeline

Build the unsliced four-elimination system:

    sage elkies-k3/scripts/export_e6_p1_sliced.sage \
      --p 31 --seed 1 --slices 0 \
      --out artifacts/local/elkies-k3/e6-base.ms

Parametrize P1_0 exactly:

    sage elkies-k3/scripts/parametrize_e6_p1_at_0.sage \
      --input artifacts/local/elkies-k3/e6-base.ms \
      --out artifacts/local/elkies-k3/e6-param0.ms

Then search many correct coordinate slices in parallel:

    python3 elkies-k3/scripts/run_e6_coordinate_slice_search.py \
      --input artifacts/local/elkies-k3/e6-param0.ms \
      --workers 8 --threads 2 --timeout 45 --seeds 32

Each job fixes r0 nonzero plus two additional coordinates, solves P1_1 -> y1 exactly, verifies P1_2=P1_3=0, exports an 8-variable / 8-equation system, then runs msolve with a bounded timeout.

The immediate goal is to obtain one finite-field hit, reconstruct the corresponding E6 surface and P1 section, and then impose P3 to descend from the three-dimensional P1 locus toward the desired rank-3 family.

See also `E6_MW3_PROGRESS_2026-08-20.md`.
