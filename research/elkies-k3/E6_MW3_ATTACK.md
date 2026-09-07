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

## Exact Neron--Severi transport recovered

The E6 frame is not merely another determinant-948 lattice.  The exact check

    sage -python elkies-k3/scripts/verify_e6_frame_ns_genus.sage

shows that

    U + (-M17)  and  U + (-E6frame)

have the same signature `(1,18)`, the same local genus, and cyclic
discriminant group `Z/948`.  Their rank is 19 while every primary
discriminant length is 1.  [Nikulin's Theorem
1.14.2](https://www.mathnet.ru/eng/im1677) therefore makes this indefinite
genus a single integral isometry class.  Thus the E6 frame belongs to the same
abstract Neron--Severi lattice as the recovered rank-17 fibration.

The previously missing integral transport has now also been recovered.  The
actual discovery path is

    rank17 --q=90--> MW7 --q=4--> MW4 --q=4--> E6/MW3.

[`scripts/verify_e6_neighbor_chain.sage`](scripts/verify_e6_neighbor_chain.sage)
reconstructs all three raw child frames from their pinned isotropic vectors and
checks the composite determinant-one isometry stored in
[`data/fibrations/e6_ns_transport_from_rank17.txt`](data/fibrations/e6_ns_transport_from_rank17.txt).
See [`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md) for the exact witnesses and
claim boundary.

This is an explicit Neron--Severi lattice transport, but not yet a
Weierstrass transport.  In particular it does not determine the rational
elliptic parameter, fiber fields of definition, or section formulas after the
three geometric neighbor operations.  The current reconstruction failure
therefore does not invalidate the E6 neighbor.

## Component/section data

The exact frame glue, not just the Shioda height equations, selects the
component-label orbit represented by:

    P1 = (1,0,1,0,0), P1.O = 0
    P2 = (1,1,2,1,1), P2.O = 1
    P3 = (2,3,0,0,0), P3.O = 0

with all pairwise section intersections equal to 2. Thus P1 and P3 can be represented polynomially while P2 is the one-denominator section.

[`scripts/recover_e6_mw3_component_glue.sage`](scripts/recover_e6_mw3_component_glue.sage)
recovers this map directly from the 17-dimensional frame.  The 32 triples
returned by the earlier height-only enumeration split into two 16-element
fiber-symmetry orbits.  Only the orbit above occurs in the frame glue.  The
reduced height lattice has automorphism group `{+I,-I}`, so the second orbit
is not another basis for this neighbor.

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

## First GF(31) reconstruction and correction

The first concrete point of the closed E6/P1 equation scheme was reconstructed
exactly over GF(31).

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

This point is rediscovered independently by multiple scan seeds.  However, an
exact discriminant audit gives `ord_1(Delta)=5`: it is an `I5` enhancement,
not a point of the required exact `I4` open stratum.  This distinction was
missed by the closed equations and is now an explicit promotion gate.

The corrected seven-core audit leaves two exact target-fiber P1 surfaces.  One
has no P2 numerator square.  The other has one signed pair of squares, but
doubling the section shows `I4` labels `(2,2)` instead of the target `(1,2)`.
Thus the audited slice contains no canonical P1+P2 seed.  See
`E6_P2_REDUCTION_2026-08-20.md`.

## Full-height correction and finite-field backtrack

The previous search called a section a P1+P2 seed after checking its fiber
components.  That is insufficient: the target Gram also requires

    <P1,P2> = -5/6.

For the selected component profiles, Shioda's formula makes this equivalent
to the exact intersection condition

    (P1+P2).O = 1.

The compiled scanner and exact checker now apply this gate before attempting
P3.  Complete scans of the declared rational chart give:

| field | core tests | exact-fiber records | valid P1 | component-valid P2 | target-height P1+P2 |
|---:|---:|---:|---:|---:|---:|
| `GF(5)` | 40,000 | 0 | 0 | 0 | 0 |
| `GF(7)` | 518,616 | 63 | 3 | 3 | 0 |
| `GF(11)` | 14,641,000 | 684 | 36 | 7 | 0 |
| `GF(13)` | 49,353,408 | 1,677 | 54 | 12 | 0 |
| `GF(17)` | 342,102,016 | 5,659 | 319 | 47 | 0 |

Thus 406,655,040 exact core tests produced 69 component-valid P2 sections
and no target two-generator height lattice.  Of those 69, 66 have pairing
`-11/6`, two have `1/6`, and one has `7/6`; none has `-5/6`.

The characteristic-11 surface that also has an independent polynomial P3 is
an exact diagnostic near miss, not a target seed.  Its height Gram is

    (1/12) * [[23,-22,16],[-22,23,-17],[16,-17,23]],

with determinant `13/48`, versus target determinant `79/16`.  This is checked
by
[`scripts/verify_e6_mw3_wrong_rank3_gf11.sage`](scripts/verify_e6_mw3_wrong_rank3_gf11.sage).

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

The next step is no longer to scan more primes in the old system.  The exact
lattice path is now pinned, and it shows what must be backtracked: execute the
three neighbor operations geometrically from a genuine explicit model, track
the fiber/component fields and sections, and derive the final Weierstrass
chart rather than assuming the split normalization.  Then rebuild the P1/P2
construction with `(P1+P2).O=1` imposed explicitly.  If that transported chart
is empty, reject this E6 neighbor as an explicit route and apply the same gates
to the A10 or A6/A4 neighbor.

See also `E6_MW3_PROGRESS_2026-08-20.md`.
