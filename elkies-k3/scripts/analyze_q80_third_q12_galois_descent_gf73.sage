#!/usr/bin/env sage
"""Recover the modular Galois action on the CM24 third-q12 base.

The characteristic-73 Jacobian was reconstructed using ``sqrt(-3)=17`` in
``reconstruct_q80_third_q12_jacobian_gf73.sage``.  Repeating the canonical
simple-infinity-place conversion with the conjugate embedding
``sqrt(-3)=-17=56`` gives the pinned samples below.  Forty-nine samples
interpolate the conjugate degree-(24,24) j-map and the remaining samples are
withheld.

The pole divisors then leave only a handful of possible fractional-linear
changes of the new base.  This script tests those candidates against the
complete rational functions.  It is an exact finite-characteristic descent
diagnostic, not a characteristic-zero lift.
"""

from sage.all import GF, Matrix, PolynomialRing, vector


finite = GF(73)
polynomial_ring = PolynomialRing(finite, "V")
V = polynomial_ring.gen()

# Output of analyze_q80_third_q12_cm24_weierstrass_gf73.sage with
# --j-root=56.  The failed values are singular fibers or local-presentation
# exceptions and are deliberately absent.
conjugate_samples = (
    (3, 6), (4, 8), (7, 67), (8, 34), (9, 19), (10, 9), (11, 27),
    (12, 45), (13, 41), (14, 51), (15, 47), (16, 72), (17, 51),
    (19, 63), (20, 34), (21, 28), (22, 20), (23, 51), (24, 7),
    (25, 13), (26, 42), (28, 6), (29, 31), (30, 48), (31, 30),
    (32, 53), (33, 19), (35, 31), (38, 9), (39, 51), (40, 53),
    (41, 19), (42, 51), (44, 44), (45, 60), (46, 30), (47, 46),
    (48, 8), (49, 7), (51, 48), (52, 41), (54, 53), (57, 19),
    (58, 34), (59, 0), (60, 0), (61, 34), (62, 60), (63, 67),
    # Withheld after interpolation.
    (64, 44), (65, 33), (1, 20), (66, 33), (68, 27),
)
training = conjugate_samples[:49]
withheld = conjugate_samples[49:]
assert len(training) == 49 and len(withheld) == 5


def interpolation_row(value, j_value):
    value = finite(value)
    j_value = finite(j_value)
    return [value**index for index in range(25)] + [
        -j_value*value**index for index in range(25)
    ]


interpolation = Matrix(
    finite, [interpolation_row(value, j_value) for value, j_value in training]
)
assert interpolation.rank() == 49
kernel = interpolation.right_kernel().basis()[0]
first_nonzero = next(value for value in kernel if value)
kernel /= first_nonzero
expected_kernel = vector(
    finite,
    [
        1, 19, 59, 37, 54, 33, 19, 14, 0, 2, 10, 15, 11, 63, 21,
        61, 35, 27, 52, 9, 67, 27, 42, 31, 43,
        63, 16, 63, 53, 61, 41, 37, 44, 71, 50, 7, 5, 53, 65, 4,
        22, 11, 21, 52, 16, 24, 46, 30, 3, 41,
    ],
)
assert kernel == expected_kernel
conjugate_numerator = polynomial_ring(list(kernel[:25]))
conjugate_denominator = polynomial_ring(list(kernel[25:]))
assert all(
    conjugate_denominator(value) != 0
    and conjugate_numerator(value)
    == finite(j_value)*conjugate_denominator(value)
    for value, j_value in withheld
)

# The original embedding's pinned j-map.
original_c4 = 4*(V+36)*(V+39)*(
    V**6+25*V**5+26*V**4+3*V**3+24*V**2+23*V+31
)
original_denominator = (
    10*(V+37)*(V+46)*(V+17)**2*(V+30)**2*(V+68)**2
    *(V+20)**7*(V+67)**7*(V**2+20*V+67)
)
original_numerator = original_c4**3


def rational_linear_roots_with_multiplicity(polynomial):
    roots = []
    for factor, multiplicity in polynomial.factor():
        if factor.degree() == 1:
            roots.append((-factor[0]/factor[1], multiplicity))
    return tuple(roots)


original_poles = rational_linear_roots_with_multiplicity(original_denominator)
conjugate_poles = rational_linear_roots_with_multiplicity(conjugate_denominator)
assert sorted(multiplicity for _, multiplicity in original_poles) == sorted(
    multiplicity for _, multiplicity in conjugate_poles
)


def points_of_multiplicity(poles, multiplicity):
    return tuple(point for point, order in poles if order == multiplicity)


def mobius_from_three_images(sources, targets):
    rows = []
    for source, target in zip(sources, targets):
        rows.append((source, 1, -target*source, -target))
    relation = Matrix(finite, rows).right_kernel()
    if relation.dimension() != 1:
        return None
    return vector(finite, relation.basis()[0])


def compose_degree_24(polynomial, numerator, denominator):
    return sum(
        polynomial[index]*numerator**index*denominator**(24-index)
        for index in range(25)
    )


old_sevens = points_of_multiplicity(original_poles, 7)
new_sevens = points_of_multiplicity(conjugate_poles, 7)
old_ones = points_of_multiplicity(original_poles, 1)
new_ones = points_of_multiplicity(conjugate_poles, 1)
assert tuple(map(len, (old_sevens, new_sevens, old_ones, new_ones))) == (2, 2, 2, 2)

candidates = []
for target_sevens in (old_sevens, tuple(reversed(old_sevens))):
    for source_one in new_ones:
        for target_one in old_ones:
            matrix = mobius_from_three_images(
                (new_sevens[0], new_sevens[1], source_one),
                (target_sevens[0], target_sevens[1], target_one),
            )
            if matrix is None:
                continue
            a, b, c, d = matrix
            if a*d-b*c == 0:
                continue
            phi_numerator = a*V+b
            phi_denominator = c*V+d
            composed_numerator = compose_degree_24(
                original_numerator, phi_numerator, phi_denominator
            )
            composed_denominator = compose_degree_24(
                original_denominator, phi_numerator, phi_denominator
            )
            if (
                conjugate_numerator*composed_denominator
                == conjugate_denominator*composed_numerator
            ):
                matrix /= next(value for value in matrix if value)
                candidates.append(tuple(matrix))

assert len(set(candidates)) == 0

# A convenient short Weierstrass normalization for the conjugate embedding.
# The interpolation kernel is projective.  Multiplication by 7 makes its
# numerator a cube and numerator-1728*denominator a square.
conjugate_scale = finite(7)
scaled_conjugate_numerator = conjugate_scale*conjugate_numerator
scaled_conjugate_denominator = conjugate_scale*conjugate_denominator
conjugate_c4 = 36*(V+13)*(V+14)*(
    V**3+49*V**2+32*V+69
)*(
    V**3+51*V**2+12*V+60
)
square_difference = (
    scaled_conjugate_numerator
    - finite(1728)*scaled_conjugate_denominator
)
square_factors = tuple(square_difference.factor())
assert all(multiplicity % 2 == 0 for _, multiplicity in square_factors)
conjugate_c6 = finite(20)
for factor, multiplicity in square_factors:
    if factor.degree() > 0:
        conjugate_c6 *= factor**(multiplicity//2)
conjugate_A = -conjugate_c4/finite(48)
conjugate_B = -conjugate_c6/finite(864)
assert conjugate_c4**3 == scaled_conjugate_numerator
assert conjugate_c6**2 == square_difference
assert (
    -finite(16)*(4*conjugate_A**3+27*conjugate_B**2)
    == scaled_conjugate_denominator
)

# A rational value of the normalized base has rational j only where the two
# split-prime embeddings agree.  This fixed-j divisor is therefore the
# modular shadow of the characteristic-zero specialization gate.
fixed_j_polynomial = (
    conjugate_numerator*original_denominator
    - conjugate_denominator*original_numerator
)
assert fixed_j_polynomial != 0
assert fixed_j_polynomial(finite(-27)) == 0
assert original_denominator(finite(-27)) != 0
assert conjugate_denominator(finite(-27)) != 0
assert fixed_j_polynomial(finite(-37)) == 0
assert original_denominator(finite(-37)) == 0

print(
    "Q80THIRDGALOISGF73|prime=73|sqrt_minus_3=17_to_56|"
    f"interpolation=49x50|rank=49|withheld={len(withheld)}|"
    f"conjugate_Delta={conjugate_denominator.factor()}",
    flush=True,
)
print(
    "Q80THIRDGALOISGF73|"
    "pgl2_descent_candidates=0|"
    f"conjugate_scale={conjugate_scale}|"
    f"conjugate_A={conjugate_A}|conjugate_B={conjugate_B}|"
    f"fixed_j_degree={fixed_j_polynomial.degree()}|"
    f"fixed_j_factor={fixed_j_polynomial.factor()}|"
    "smooth_fixed_j_residue=V:-27|singular_artifact=V:-37|"
    "status=PASS_NO_MODULAR_PGL2_DESCENT",
    flush=True,
)
