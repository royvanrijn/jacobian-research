# E6 MW3 progress — 2026-08-20

> **Historical snapshot.** This note records the 2026-08-20 frontier of the
> alternate E6/MW3 attack.  The split chart below was rejected, and the
> selected H3 equation corridor has since reached the certified rootless
> `24I1/MW17` endpoint.  Use [`README.md`](README.md) for the current residual
> 2-descent priority and the
> [`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md) for chronology.

## Neighbor-chain correction

The E6 frame's actual lattice provenance is now exact.  It was reached by

```text
rank17 --q=90--> MW7 --q=4--> MW4 --q=4--> E6/MW3,
```

not by a direct unrecorded change of basis.  All three primitive isotropic
vectors and the composite integral isometry are checked by
[`scripts/verify_e6_neighbor_chain.sage`](scripts/verify_e6_neighbor_chain.sage);
see [`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md).

This changes the interpretation of the finite-field work below.  The searched
split chart was inferred from the Kodaira symbols, not produced by executing
the recovered neighbor chain on the original K3.  Its exhaustive emptiness is
a valid bounded result about that chart, but it is not a rejection of the E6
neighbor.  At this snapshot the next step was a geometric backtrack through
the pinned chain, not another larger scan; that historical route decision has
since been superseded by the completed H3 corridor described above.

## Exact-open correction

The first reconstruction described below satisfies the displayed closed
equations, but an exact audit gives `ord_1(Delta)=5`, not `4`.  It is therefore
an `I5` boundary point, not a point of the intended exact-fiber open stratum.
The subsequent seven-core audit found two genuine target-fiber P1 surfaces,
but no canonical P2: the only numerator-square P2 has `I4` labels `(2,2)`
rather than the target `(1,2)`.  A later compiled all-chart audit over
`GF(5),GF(7),GF(11),GF(13),GF(17)` exposed a second missing promotion gate:
all 69 component-valid P2 sections have the wrong P1/P2 height pairing. See
[`E6_P2_REDUCTION_2026-08-20.md`](E6_P2_REDUCTION_2026-08-20.md) for the
corrected promotion gates and complete bounded results.

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

The searched chart normalized IV* at infinity, I4 fibers at 0 and 1, and the
two I2 fibers at lambda and mu.  This is a convenient split ansatz, not yet a
transported normalization of the target K3.  The IV* condition gives
`deg A <= 5` and `deg B <= 8` in short Weierstrass form.

## Canonical P1 component data

Exact discriminant-glue recovery from the 17-dimensional frame selects the
orbit represented by

    P1 = (1,0,1,0,0), P1.O = 0
    P2 = (1,1,2,1,1), P2.O = 1
    P3 = (2,3,0,0,0), P3.O = 0

with all pairwise intersections equal to 2. P1 and P3 are polynomial sections in this normalization.

This is certified by
[`scripts/recover_e6_mw3_component_glue.sage`](scripts/recover_e6_mw3_component_glue.sage).
The other 16-element orbit from the height-only component enumeration is
excluded by the frame glue and is not related by a Mordell-Weil basis change.

## Missing pair-intersection gate

Fiber components and section self-heights do not determine the off-diagonal
height entries.  For the certified profiles, the required equality

    <P1,P2> = -5/6

is equivalent to

    (P1+P2).O = 1.

The exact checker now computes the full Shioda height matrix before any P3
search.  Exhaustion of the rational chart over the five small fields tested
406,655,040 cores.  It found 69 P2 sections passing the old component gates
and zero passing this height gate.  The field-by-field counts are recorded in
[`E6_MW3_ATTACK.md`](E6_MW3_ATTACK.md).

## High-coefficient triangular chain

For P1, the section-first ansatz simplifies to

    X = (s1-x0-x1)t^2 + x1*t + x0
    Y = t^4 + (-y0-y1-y2-1)t^3 + y2*t^2 + y1*t + y0.

The four highest section coefficients eliminate successively:

    P1_7 -> b4
    P1_6 -> b5
    P1_5 -> y2
    P1_4 -> a3.

The resulting system has 13 variables and 12 equations: eight fiber equations plus P1_0..P1_3.

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

1. msolve itself is functioning correctly. Earlier SIGSEGVs came from an exporter bug that wrote literal `\\n` characters, producing a malformed one-line input file.
2. Dense three-slice systems on the pre-parametrized representation can return `No solution` quickly because dependent section equations were still being counted as constraints.
3. A one-slice 13x13 system over GF(101) reaches large F4 matrices around degrees 18-19. Repeating over GF(31) gives nearly the same progression, showing the complexity is structural rather than prime-specific.
4. Dense affine elimination of the three final slice equations is computationally poor: substituting 9-term linear forms into degree-20+ fiber equations causes severe multinomial expansion.
5. Coordinate slices are much cheaper. Replacing selected variables by constants keeps the final 8x8 systems small enough for bounded parallel msolve exploration.
6. The especially effective slice family fixes r0 and s0; fixing x1 as the third coordinate reduces the eight equations to only about 570 monomials total, with maximum degree 11. Even then direct F4 solving grows rapidly around degrees 14-16, so structural elimination is preferable.

## GCD reduction for the two I2 fibers

After fixing the cheap coordinate slice `(r0,s0,x1)` and eliminating `y1` exactly, the two I2-fiber triples are symmetric in lambda and mu. Eliminating the singular-root variables `sl,sm` produces two polynomial conditions in a generic fiber coordinate `z`.

The first apparent common factor is always

    z^3 (z-1),

coming from the already normalized fibers at 0 and 1. Additional repeated factors such as `(z-r)^2` are collision loci `lambda=mu` and must also be saturated away.

A fast finite-field line scan over GF(31), after saturating all powers of `z`, `z-1`, and rejecting repeated single-root factors, finds genuine squarefree common factors of degree at least two. This avoids symbolic subresultants and is several orders of magnitude faster than Groebner solving.

## First closed-equation GF(31) point

The first fully reconstructed distinct-root solution of the closed equations
is obtained on the seed-1 coordinate slice

    r0 = 4
    s0 = 18
    x1 = 27

with reduced core parameters

    a1 = 4
    a2 = 16
    a4 = 6
    s1 = 23.

The saturated common polynomial for the two I2 locations is

    z^2 + 20 z + 29 = (z-24)(z-18)  mod 31.

Taking

    lambda = 24
    mu     = 18

reconstructs

    sl = 23
    sm = 4,

and all eight reduced fiber equations vanish exactly over GF(31):

    [0,0,0,0,0,0,0,0].

Swapping the two I2 fibers also works exactly:

    lambda = 18, mu = 24,
    sl = 4, sm = 23.

The same core point is rediscovered independently by several random line seeds, confirming that it is not an artifact of a particular scan path.

Additional genuine squarefree degree-2 and degree-3 common factors were found over GF(31), including both split and irreducible quadratic cases. Thus the desired distinct-I2 locus is nonempty modulo 31 and appears repeatedly in the reduced search.

It is a useful algebraic reconstruction and smooth point of the closed
equation scheme, but the omitted exact-multiplicity open test places it on the
`I5` boundary.  It is not a promoted point of the intended E6/P1 fiber
configuration.

## Current search/reconstruction strategy

1. Build the unsliced four-elimination E6/P1 system.
2. Parametrize P1_0 exactly.
3. Use cheap coordinate slices fixing `r0,s0` and one additional coordinate.
4. Eliminate P1_1 -> y1 exactly and verify P1_2=P1_3=0.
5. Saturate normalized-fiber factors from the symmetric lambda/mu elimination.
6. Search for squarefree common factors of degree >=2 rather than running large Groebner bases.
7. Reconstruct `sl,sm` and verify all eight reduced fiber equations.
8. Reconstruct the four earlier triangular eliminations `a3,y2,b5,b4` and the complete A(t), B(t), X(t), Y(t).
9. Verify the full Weierstrass section identity, **exact** Kodaira
   multiplicities, squarefree residual discriminant, and deeper component
   labels.
10. Compute the full P1/P2 Shioda pairing and require `(P1+P2).O=1`.
11. Use only promoted modular points as seeds for lifting / multi-prime
    rational reconstruction, then impose P3 to cut toward the target
    one-dimensional rank-3 family.

Recommended reconstruction command after generating the base metadata:

    sage elkies-k3/scripts/reconstruct_e6_gf31_point.sage \
      --meta artifacts/local/elkies-k3/e6-base.meta.txt

The immediate milestone is a corrected P1/P2 equation system with the missing
intersection point parameterized from the start.  Enlarging the old
finite-field scan is not justified by the exhaustive small-field results.
