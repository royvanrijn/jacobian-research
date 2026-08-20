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
- Coordinate slices are dramatically cheaper. Fixing r0,s0,x1 is especially effective, reducing the final 8x8 probe to about 570 monomials with maximum degree 11, although direct F4 solving still grows rapidly around degrees 14-16.
- The best current method is not Groebner solving but a fast finite-field GCD scan for the symmetric lambda/mu fiber conditions after saturating the normalized-fiber and collision factors.

## Verified GF(31) point

The first concrete point on the intended E6/P1 locus has now been verified exactly over GF(31).

Coordinate slice:

    r0 = 4
    s0 = 18
    x1 = 27

Reduced core:

    a1 = 4
    a2 = 16
    a4 = 6
    s1 = 23

The saturated common polynomial for the two I2 positions is

    z^2 + 20*z + 29 = (z-24)(z-18).

Thus one labeling is

    lambda = 24
    mu     = 18
    sl     = 23
    sm     = 4.

All eight reduced fiber equations evaluate to zero. Swapping the two I2 fibers also works:

    lambda = 18, mu = 24, sl = 4, sm = 23.

This point is rediscovered independently by multiple scan seeds. Other split and irreducible squarefree quadratic common factors also occur, so the distinct-I2 locus is demonstrably nonempty mod 31.

## Current reproducible pipeline

Build the unsliced four-elimination system:

    sage elkies-k3/scripts/export_e6_p1_sliced.sage \
      --p 31 --seed 1 --slices 0 \
      --out artifacts/local/elkies-k3/e6-base.ms

Parametrize P1_0 exactly:

    sage elkies-k3/scripts/parametrize_e6_p1_at_0.sage \
      --input artifacts/local/elkies-k3/e6-base.ms \
      --out artifacts/local/elkies-k3/e6-param0.ms

For bounded coordinate-slice exploration:

    python3 elkies-k3/scripts/run_e6_coordinate_slice_search.py \
      --input artifacts/local/elkies-k3/e6-param0.ms \
      --workers 8 --threads 2 --timeout 45 --seeds 32

For the now-known modular point, reconstruct the full eliminated chain and verify the original section-first model:

    sage elkies-k3/scripts/reconstruct_e6_gf31_point.sage \
      --meta artifacts/local/elkies-k3/e6-base.meta.txt

The immediate goal is now end-to-end verification of the complete A(t),B(t),P1 model at this GF(31) point, followed by lifting/multi-prime reconstruction and then imposing P3.

See also `E6_MW3_PROGRESS_2026-08-20.md`.
