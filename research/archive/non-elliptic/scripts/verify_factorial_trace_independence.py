#!/usr/bin/env python3
"""Exact regressions for factorial/gamma-affine trace independence.

The all-scale theorem is proved in
``extended-geometry/FACTORIAL_TRACE_INDEPENDENCE.md``.  This dependency-free
checker audits signed shift-orbit reconstruction, the complete rational-offset
gamma signature, rational coboundary certificates, affine-shift
canonicalization, exact integer-affine boundary transfers, Stirling-invariant
examples, residue-class localization, and the exact Frobenius-dilation
collapse of p-adic factorial valuations.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import factorial
from random import Random


MAX_SLOPE = 7
COEFFICIENTS = range(-2, 3)
GAMMA_MAX_SLOPE = 6
GAMMA_OFFSET_DENOMINATOR = 12

GammaFactor = tuple[int, Fraction, int]
GammaRay = tuple[GammaFactor, ...]
RationalGammaFactor = tuple[Fraction, Fraction, int]
RationalGammaRay = tuple[RationalGammaFactor, ...]
AffineFactor = tuple[int, int]
AffineProfile = tuple[AffineFactor, ...]


def orbit_key(numerator: int, denominator: int) -> Fraction:
    """Translation orbit of -numerator/denominator in Q/Z.

    Negating every key is an automorphism of Q/Z, so the positive fractional
    representative is equivalent to the root convention in the proof.
    """

    return Fraction(numerator, denominator) % 1


def residue_class(value: Fraction) -> Fraction:
    """Canonical representative of a rational translation orbit in Q/Z."""

    return value % 1


SLOPE_ORBITS = {
    slope: tuple(orbit_key(numerator, slope) for numerator in range(1, slope + 1))
    for slope in range(1, MAX_SLOPE + 1)
}


def orbit_signature(vector: tuple[int, ...]) -> dict[Fraction, int]:
    signature: defaultdict[Fraction, int] = defaultdict(int)
    for slope, multiplicity in enumerate(vector, 1):
        for key in SLOPE_ORBITS[slope]:
            signature[key] += multiplicity
    return {key: value for key, value in signature.items() if value}


def recover_signed_vector(signature: dict[Fraction, int]) -> tuple[int, ...]:
    """Recover signed slope multiplicities by descending maximal slope."""

    remaining: defaultdict[Fraction, int] = defaultdict(int, signature)
    recovered = [0] * MAX_SLOPE
    for slope in range(MAX_SLOPE, 0, -1):
        pivot = orbit_key(1, slope)
        multiplicity = remaining[pivot]
        recovered[slope - 1] = multiplicity
        for key in SLOPE_ORBITS[slope]:
            remaining[key] -= multiplicity
            if not remaining[key]:
                del remaining[key]
    assert not remaining
    return tuple(recovered)


def verify_signed_orbit_reconstruction() -> int:
    checked = 0
    for vector in product(COEFFICIENTS, repeat=MAX_SLOPE):
        signature = orbit_signature(vector)
        assert recover_signed_vector(signature) == vector
        if any(vector):
            assert signature
            checked += 1
        else:
            assert not signature
    return checked


def gamma_shift_divisor(ray: GammaRay) -> dict[Fraction, int]:
    """Finite zero/pole divisor of G(n+1)/G(n), with constants removed."""

    divisor: defaultdict[Fraction, int] = defaultdict(int)
    for slope, offset, multiplicity in ray:
        assert slope > 0
        assert multiplicity
        for shift in range(slope):
            root = -(offset + shift) / slope
            divisor[root] += multiplicity
    return {root: value for root, value in divisor.items() if value}


def divisor_orbit_signature(
    divisor: dict[Fraction, int], translation_step: int = 1
) -> dict[Fraction, int]:
    assert translation_step > 0
    signature: defaultdict[Fraction, int] = defaultdict(int)
    for root, multiplicity in divisor.items():
        signature[root % translation_step] += multiplicity
    return {key: value for key, value in signature.items() if value}


def gamma_orbit_signature(ray: GammaRay) -> dict[Fraction, int]:
    """Signed unit-slope gamma offsets modulo integer translation.

    Negating the shift-quotient roots is an automorphism of Q/Z, so this is
    equivalent to ``divisor_orbit_signature(gamma_shift_divisor(ray))`` while
    matching the offsets in Gauss's multiplication formula directly.
    """

    signature: defaultdict[Fraction, int] = defaultdict(int)
    for slope, offset, multiplicity in ray:
        assert slope > 0
        assert multiplicity
        for shift in range(slope):
            signature[residue_class((offset + shift) / slope)] += multiplicity
    return {key: value for key, value in signature.items() if value}


def gamma_signature_key(ray: GammaRay) -> tuple[tuple[Fraction, int], ...]:
    return tuple(sorted(gamma_orbit_signature(ray).items()))


def subtract_divisors(
    left: dict[Fraction, int], right: dict[Fraction, int]
) -> dict[Fraction, int]:
    difference: defaultdict[Fraction, int] = defaultdict(int, left)
    for root, multiplicity in right.items():
        difference[root] -= multiplicity
    return {root: value for root, value in difference.items() if value}


def integrate_shift_divisor(
    divisor: dict[Fraction, int], translation_step: int = 1
) -> dict[Fraction, int]:
    """Construct E with D(x)=E(x+s)-E(x) on every s-translation orbit."""

    assert translation_step > 0
    orbits: defaultdict[Fraction, dict[int, int]] = defaultdict(dict)
    for root, multiplicity in divisor.items():
        residue = root % translation_step
        integer_part = (root - residue) / translation_step
        assert integer_part.denominator == 1
        step = int(integer_part)
        orbits[residue][step] = orbits[residue].get(step, 0) + multiplicity

    witness: dict[Fraction, int] = {}
    for residue, orbit_divisor in orbits.items():
        assert sum(orbit_divisor.values()) == 0
        first = min(orbit_divisor)
        last = max(orbit_divisor)
        running = 0
        for step in range(first, last + 1):
            running += orbit_divisor.get(step, 0)
            if running:
                witness[residue + (step + 1) * translation_step] = running
        assert running == 0

    audit_points = set(divisor) | set(witness)
    audit_points.update(root - translation_step for root in witness)
    for root in audit_points:
        assert (
            witness.get(root + translation_step, 0) - witness.get(root, 0)
            == divisor.get(root, 0)
        )
    return witness


def evaluate_factored_divisor(
    divisor: dict[Fraction, int], value: int
) -> Fraction:
    result = Fraction(1)
    for root, multiplicity in divisor.items():
        result *= (Fraction(value) - root) ** multiplicity
    return result


def gamma_shift_quotient(ray: GammaRay, value: int) -> Fraction:
    result = Fraction(1)
    for slope, offset, multiplicity in ray:
        for shift in range(slope):
            result *= (slope * value + offset + shift) ** multiplicity
    return result


def gamma_leading_constant(ray: GammaRay) -> Fraction:
    result = Fraction(1)
    for slope, _offset, multiplicity in ray:
        result *= Fraction(slope) ** (slope * multiplicity)
    return result


def certify_gamma_equivalence(
    left: GammaRay, right: GammaRay
) -> tuple[Fraction, dict[Fraction, int]]:
    """Return the constant and rational-coboundary divisor certificate."""

    assert gamma_orbit_signature(left) == gamma_orbit_signature(right)
    divisor = subtract_divisors(
        gamma_shift_divisor(left), gamma_shift_divisor(right)
    )
    assert not divisor_orbit_signature(divisor)
    witness = integrate_shift_divisor(divisor)
    constant = gamma_leading_constant(left) / gamma_leading_constant(right)

    checked = 0
    for value in range(20, 40):
        try:
            left_value = gamma_shift_quotient(left, value)
            right_value = gamma_shift_quotient(right, value)
            witness_value = evaluate_factored_divisor(witness, value)
            shifted_witness_value = evaluate_factored_divisor(witness, value + 1)
            assert left_value / right_value == (
                constant * shifted_witness_value / witness_value
            )
        except ZeroDivisionError:
            continue
        checked += 1
    assert checked >= 8
    return constant, witness


def verify_coboundary_integration() -> int:
    """Exhaust finite divisors on one nonintegral translation orbit."""

    positions = tuple(Fraction(step) + Fraction(1, 3) for step in range(-3, 4))
    checked = 0
    for coefficients in product(COEFFICIENTS, repeat=len(positions)):
        if not any(coefficients) or sum(coefficients):
            continue
        divisor = {
            root: multiplicity
            for root, multiplicity in zip(positions, coefficients)
            if multiplicity
        }
        assert not divisor_orbit_signature(divisor)
        integrate_shift_divisor(divisor)
        checked += 1
    return checked


def verify_gamma_examples() -> None:
    duplication = ((2, Fraction(1), 1),)
    duplication_refined = (
        (1, Fraction(1, 2), 1),
        (1, Fraction(1), 1),
    )
    constant, witness = certify_gamma_equivalence(
        duplication, duplication_refined
    )
    assert constant == 4
    assert not witness

    triplication = ((3, Fraction(1), 1),)
    triplication_refined = (
        (1, Fraction(1, 3), 1),
        (1, Fraction(2, 3), 1),
        (1, Fraction(1), 1),
    )
    constant, witness = certify_gamma_equivalence(
        triplication, triplication_refined
    )
    assert constant == 27
    assert not witness

    shifted = ((2, Fraction(2), 1),)
    unshifted = ((2, Fraction(1), 1),)
    constant, witness = certify_gamma_equivalence(shifted, unshifted)
    assert constant == 1
    assert witness == {Fraction(-1, 2): 1}

    central_numerator = ((2, Fraction(1), 1),)
    central_denominator = ((1, Fraction(1), 2),)
    assert gamma_orbit_signature(central_numerator) != gamma_orbit_signature(
        central_denominator
    )


def verify_gauss_refinements() -> tuple[int, int]:
    """Test all configured slope divisors, then add integer offset shifts."""

    refinements = 0
    nonconstant_coboundaries = 0
    for slope in range(2, 13):
        for divisor in range(2, slope + 1):
            if slope % divisor:
                continue
            for numerator in range(-3, 4):
                offset = Fraction(numerator, 7)
                left = ((slope, offset, 1),)
                right = tuple(
                    (
                        slope // divisor,
                        (offset + shift) / divisor,
                        1,
                    )
                    for shift in range(divisor)
                )
                constant, witness = certify_gamma_equivalence(left, right)
                assert constant == Fraction(divisor) ** slope
                assert not witness
                refinements += 1

                shifted_right = tuple(
                    (new_slope, new_offset + (index % 3) - 1, multiplicity)
                    for index, (new_slope, new_offset, multiplicity) in enumerate(
                        right
                    )
                )
                shifted_constant, shifted_witness = certify_gamma_equivalence(
                    left, shifted_right
                )
                assert shifted_constant == constant
                nonconstant_coboundaries += bool(shifted_witness)
    return refinements, nonconstant_coboundaries


def verify_gamma_signature_census() -> tuple[int, int, int, int, int]:
    """Classify products of at most three small rational-offset atoms."""

    atoms: tuple[GammaFactor, ...] = tuple(
        (slope, Fraction(numerator, GAMMA_OFFSET_DENOMINATOR), 1)
        for slope in range(1, GAMMA_MAX_SLOPE + 1)
        for numerator in range(GAMMA_OFFSET_DENOMINATOR)
    )
    representatives: dict[tuple[tuple[Fraction, int], ...], GammaRay] = {}
    class_sizes: defaultdict[tuple[tuple[Fraction, int], ...], int] = defaultdict(
        int
    )
    rays = 0
    collisions = 0
    for length in range(1, 4):
        for indices in combinations_with_replacement(range(len(atoms)), length):
            ray = tuple(atoms[index] for index in indices)
            signature = gamma_signature_key(ray)
            rays += 1
            class_sizes[signature] += 1
            if signature in representatives:
                _constant, witness = certify_gamma_equivalence(
                    ray, representatives[signature]
                )
                # Offsets were already normalized to [0,1), so Gauss
                # refinement leaves only an exponential constant.
                assert not witness
                collisions += 1
            else:
                representatives[signature] = ray

    nontrivial_classes = sum(size > 1 for size in class_sizes.values())
    largest_class = max(class_sizes.values())
    assert len(atoms) == 72
    assert rays == 67_524
    assert len(representatives) == 66_140
    assert collisions == 1_384
    assert nontrivial_classes == 1_219
    assert largest_class == 4
    return (
        rays,
        len(representatives),
        collisions,
        nontrivial_classes,
        largest_class,
    )


def verify_signed_gamma_property_cases() -> int:
    """Randomized signed rays under Gauss refinements and integer shifts."""

    generator = Random(0xFA_C7_0A)
    checked = 0
    for _case in range(1_000):
        factors: list[GammaFactor] = []
        transformed: list[GammaFactor] = []
        for _factor in range(generator.randint(1, 6)):
            slope = generator.randint(1, 8)
            offset = Fraction(generator.randint(-6, 6), generator.randint(1, 7))
            multiplicity = generator.choice((-2, -1, 1, 2))
            factors.append((slope, offset, multiplicity))

            divisors = [
                divisor
                for divisor in range(2, slope + 1)
                if slope % divisor == 0
            ]
            if divisors and generator.randrange(2):
                divisor = generator.choice(divisors)
                transformed.extend(
                    (
                        slope // divisor,
                        (offset + shift) / divisor,
                        multiplicity,
                    )
                    for shift in range(divisor)
                )
            else:
                integer_shift = generator.choice((-3, -2, -1, 1, 2, 3))
                transformed.append(
                    (slope, offset + integer_shift, multiplicity)
                )

        ray = tuple(factors)
        transformed_ray = tuple(transformed)
        assert gamma_orbit_signature(ray) == gamma_orbit_signature(
            transformed_ray
        )
        certify_gamma_equivalence(ray, transformed_ray)
        checked += 1

    for slope in range(1, 11):
        ray = ((slope, Fraction(2, 7), 1),)
        perturbed = ((slope, Fraction(2, 7) + Fraction(1, 101), 1),)
        assert gamma_orbit_signature(ray) != gamma_orbit_signature(perturbed)
    return checked


def canonical_factorial_ray(vector: tuple[int, ...], scale: int) -> Fraction:
    value = Fraction(1)
    for slope, multiplicity in enumerate(vector, 1):
        value *= Fraction(factorial(slope * scale)) ** multiplicity
    return value


def verify_integer_affine_shifts() -> int:
    """Check positive and negative shifts against their rational prefactors."""

    # H(n)=(2n+2)!/(n-1)! has canonical vector e_2-e_1 and
    # H/Phi=n(2n+1)(2n+2).
    vector = (-1, 1)
    checked = 0
    for scale in range(1, 13):
        shifted = Fraction(
            factorial(2 * scale + 2),
            factorial(scale - 1),
        )
        canonical = canonical_factorial_ray(vector, scale)
        prefactor = scale * (2 * scale + 1) * (2 * scale + 2)
        assert shifted == prefactor * canonical
        checked += 1

    # A second denominator shift checks more than one removed factor:
    # (3n-2)!/(3n)!=1/((3n-1)(3n)).
    for scale in range(1, 13):
        shifted = Fraction(factorial(3 * scale - 2), factorial(3 * scale))
        predicted = Fraction(1, (3 * scale - 1) * (3 * scale))
        assert shifted == predicted
        checked += 1
    return checked


def affine_successor_divisor(
    profile: dict[AffineFactor, int],
) -> dict[Fraction, int]:
    """Exact finite divisor of the consecutive factorial quotient.

    The factor ``(a*n+b)!`` contributes the roots
    ``-(b+1)/a,...,-(b+a)/a`` to its shift quotient.  Negating every root
    is harmless and makes the boundary-transfer presentation more legible.
    """

    divisor: defaultdict[Fraction, int] = defaultdict(int)
    for (slope, offset), multiplicity in profile.items():
        assert slope > 0
        if not multiplicity:
            continue
        for shift in range(1, slope + 1):
            divisor[Fraction(offset + shift, slope)] += multiplicity
    return {root: value for root, value in divisor.items() if value}


def affine_slope_vector(
    profile: dict[AffineFactor, int],
) -> dict[int, int]:
    vector: defaultdict[int, int] = defaultdict(int)
    for (slope, _offset), multiplicity in profile.items():
        vector[slope] += multiplicity
    return {slope: value for slope, value in vector.items() if value}


def affine_profile_difference(
    left: AffineProfile,
    right: AffineProfile,
) -> dict[AffineFactor, int]:
    difference: defaultdict[AffineFactor, int] = defaultdict(int)
    for factor in left:
        difference[factor] += 1
    for factor in right:
        difference[factor] -= 1
    return {factor: value for factor, value in difference.items() if value}


def boundary_transfer_decomposition(
    profile: dict[AffineFactor, int],
) -> tuple[tuple[int, AffineFactor, AffineFactor], ...]:
    """Decompose a zero exact divisor into elementary boundary transfers.

    A record ``(m,(a,k),(c,l))`` represents ``m`` times

        E(a,k)-E(a,k-1)-E(c,l)+E(c,l-1),

    where ``k/a=l/c``.  Its factorial product is the constant ``(a/c)^m``.
    The proof in the canonical note shows that this succeeds for every
    finitely supported integer profile with zero exact divisor.
    """

    if affine_successor_divisor(profile):
        raise ValueError("the affine factorial divisor is nonzero")

    baselines: defaultdict[int, int] = defaultdict(int)
    increments: defaultdict[AffineFactor, int] = defaultdict(int)
    for (slope, offset), multiplicity in profile.items():
        baselines[slope] += multiplicity
        if offset > 0:
            for boundary in range(1, offset + 1):
                increments[(slope, boundary)] += multiplicity
        elif offset < 0:
            for boundary in range(offset + 1, 1):
                increments[(slope, boundary)] -= multiplicity

    baselines = defaultdict(
        int,
        {slope: value for slope, value in baselines.items() if value},
    )
    if baselines:
        raise AssertionError(
            f"zero divisor retained nonzero slope baselines: {dict(baselines)}"
        )

    by_boundary: defaultdict[Fraction, list[tuple[AffineFactor, int]]] = (
        defaultdict(list)
    )
    for factor, multiplicity in increments.items():
        if multiplicity:
            slope, boundary = factor
            by_boundary[Fraction(boundary, slope)].append(
                (factor, multiplicity)
            )

    transfers: list[tuple[int, AffineFactor, AffineFactor]] = []
    for rational_boundary, terms in sorted(by_boundary.items()):
        if sum(multiplicity for _factor, multiplicity in terms):
            raise AssertionError(
                "zero divisor retained an unmatched translation edge at "
                f"{rational_boundary}"
            )
        terms.sort()
        reference = terms[0][0]
        for factor, multiplicity in terms[1:]:
            if multiplicity:
                transfers.append((multiplicity, factor, reference))

    reconstructed: defaultdict[AffineFactor, int] = defaultdict(int)
    for multiplicity, (slope, boundary), (
        reference_slope,
        reference_boundary,
    ) in transfers:
        assert Fraction(boundary, slope) == Fraction(
            reference_boundary, reference_slope
        )
        reconstructed[(slope, boundary)] += multiplicity
        reconstructed[(slope, boundary - 1)] -= multiplicity
        reconstructed[(reference_slope, reference_boundary)] -= multiplicity
        reconstructed[(reference_slope, reference_boundary - 1)] += multiplicity
    cleaned = {
        factor: value for factor, value in reconstructed.items() if value
    }
    if cleaned != profile:
        raise AssertionError(
            f"boundary transfers reconstruct {cleaned}, expected {profile}"
        )
    return tuple(transfers)


def boundary_transfer_constant(
    transfers: tuple[tuple[int, AffineFactor, AffineFactor], ...],
) -> Fraction:
    constant = Fraction(1)
    for multiplicity, (slope, _boundary), (
        reference_slope,
        _reference_boundary,
    ) in transfers:
        constant *= Fraction(slope, reference_slope) ** multiplicity
    return constant


def affine_factorial_value(
    profile: dict[AffineFactor, int],
    scale: int,
) -> Fraction:
    value = Fraction(1)
    for (slope, offset), multiplicity in profile.items():
        argument = slope * scale + offset
        if argument < 0:
            raise ValueError("factorial profile is not defined at this scale")
        value *= Fraction(factorial(argument)) ** multiplicity
    return value


def verify_exact_affine_boundary_census() -> tuple[int, int, int, int, int]:
    """Classify a bounded universe by the exact, not orbit, divisor.

    Every nontrivial collision is certified as a sum of the elementary
    boundary transfers from Theorem 5.2.  This is a regression for the
    unbounded telescoping proof, not a bounded extrapolation.
    """

    atoms: tuple[AffineFactor, ...] = tuple(
        (slope, offset)
        for slope in range(1, 6)
        for offset in range(-3, 4)
    )
    representatives: dict[
        tuple[tuple[Fraction, int], ...], AffineProfile
    ] = {}
    class_sizes: defaultdict[tuple[tuple[Fraction, int], ...], int] = (
        defaultdict(int)
    )
    profiles = 0
    collisions = 0
    for length in range(1, 5):
        for indices in combinations_with_replacement(range(len(atoms)), length):
            profile = tuple(atoms[index] for index in indices)
            signed_profile = affine_profile_difference(profile, ())
            signature = tuple(
                sorted(affine_successor_divisor(signed_profile).items())
            )
            class_sizes[signature] += 1
            profiles += 1
            representative = representatives.get(signature)
            if representative is None:
                representatives[signature] = profile
                continue

            difference = affine_profile_difference(profile, representative)
            if affine_slope_vector(difference):
                raise AssertionError(
                    "an exact affine-divisor collision changed its slope vector"
                )
            transfers = boundary_transfer_decomposition(difference)
            predicted = boundary_transfer_constant(transfers)
            if (
                affine_factorial_value(difference, 8) != predicted
                or affine_factorial_value(difference, 9) != predicted
            ):
                raise AssertionError(
                    "a boundary-transfer constant failed direct factorial replay"
                )
            collisions += 1

    nontrivial_classes = sum(size > 1 for size in class_sizes.values())
    largest_class = max(class_sizes.values())
    assert len(atoms) == 35
    assert profiles == 82_250
    assert len(representatives) == 72_383
    assert collisions == 9_867
    assert nontrivial_classes == 8_253
    assert largest_class == 6
    return (
        profiles,
        len(representatives),
        collisions,
        nontrivial_classes,
        largest_class,
    )


def factorial_product(parts: tuple[int, ...], scale: int) -> int:
    value = 1
    for part in parts:
        value *= factorial(part * scale)
    return value


def verify_one_scale_collision() -> None:
    left = (4, 1, 1, 1)
    right = (3, 2, 2)
    assert sum(left) == sum(right) == 7
    assert factorial_product(left, 1) == factorial_product(right, 1)
    assert factorial_product(left, 2) != factorial_product(right, 2)


def entropy_base(parts: tuple[int, ...]) -> int:
    value = 1
    for part in parts:
        value *= part**part
    return value


def inverse_moment(parts: tuple[int, ...], odd_power: int) -> Fraction:
    assert odd_power > 0 and odd_power % 2 == 1
    return sum((Fraction(1, part**odd_power) for part in parts), Fraction(0))


def verify_stirling_hierarchy() -> None:
    left = (12, 6, 4, 4, 4, 1)
    right = (9, 8, 8, 2, 2, 2)
    assert sum(left) == sum(right) == 31
    assert len(left) == len(right) == 6

    left_product = 1
    right_product = 1
    for part in left:
        left_product *= part
    for part in right:
        right_product *= part
    assert left_product == right_product == 4608
    assert entropy_base(left) == entropy_base(right)
    assert inverse_moment(left, 1) == 2
    assert inverse_moment(right, 1) == Fraction(67, 36)


def factorial_valuation(size: int, prime: int) -> int:
    assert size >= 0
    assert prime >= 2
    valuation = 0
    while size:
        size //= prime
        valuation += size
    return valuation


def primitive_valuation_profile(
    vector: tuple[int, ...], prime: int
) -> tuple[dict[int, int], int]:
    """Reduce p-power slope towers to primitive slopes plus linear drift."""

    primitive: defaultdict[int, int] = defaultdict(int)
    drift = 0
    for slope, multiplicity in enumerate(vector, 1):
        if not multiplicity:
            continue
        exponent = 0
        primitive_slope = slope
        while primitive_slope % prime == 0:
            primitive_slope //= prime
            exponent += 1
        primitive[primitive_slope] += multiplicity
        drift += (
            multiplicity
            * primitive_slope
            * (prime**exponent - 1)
            // (prime - 1)
        )
    return (
        {slope: value for slope, value in primitive.items() if value},
        drift,
    )


def valuation_ray(vector: tuple[int, ...], scale: int, prime: int) -> int:
    return sum(
        multiplicity * factorial_valuation(slope * scale, prime)
        for slope, multiplicity in enumerate(vector, 1)
    )


def verify_frobenius_dilation_valuations() -> int:
    """Audit the exact p-primitive reduction and its four-ray kernel."""

    checked = 0
    for prime in (2, 3, 5):
        for vector in product(range(-1, 2), repeat=6):
            primitive, drift = primitive_valuation_profile(vector, prime)
            for scale in (1, 2, 3, 5, 8, 13):
                reduced = drift * scale + sum(
                    multiplicity
                    * factorial_valuation(primitive_slope * scale, prime)
                    for primitive_slope, multiplicity in primitive.items()
                )
                assert valuation_ray(vector, scale, prime) == reduced
            checked += 1

        # b(F_(pa)-F_a)-a(F_(pb)-F_b)=0.  Choosing coprime a,b
        # avoids a degenerate overlap of the four displayed slopes.
        a = 1
        b = prime + 1
        relation = [0] * (prime * b)
        relation[prime * a - 1] += b
        relation[a - 1] -= b
        relation[prime * b - 1] -= a
        relation[b - 1] += a
        relation_vector = tuple(relation)
        assert any(relation_vector)
        for scale in range(1, 80):
            assert valuation_ray(relation_vector, scale, prime) == 0
    return checked


def verify_sic_moment_family_signatures() -> int:
    """Check the injective signatures of M_(d,r) through d=48."""

    signatures: dict[tuple[tuple[Fraction, int], ...], tuple[int, int]] = {}
    checked = 0
    for degree in range(4, 49):
        for seed_power in range(1, degree // 4 + 1):
            ray: GammaRay = (
                (degree, Fraction(3), 1),
                (seed_power, Fraction(1), 2),
                (2 * seed_power, Fraction(2), -1),
            )
            signature = gamma_signature_key(ray)
            assert signature not in signatures, (
                signatures.get(signature),
                (degree, seed_power),
            )
            signatures[signature] = (degree, seed_power)
            checked += 1
    assert checked == 276
    return checked


def verify_mfold_periodic_symmetry() -> None:
    """Check the periodic exponential-polynomial symmetry of m-fold terms."""

    for value in range(24):
        first = 1
        second = 1 if value % 2 == 0 else 2
        periodic_coefficient = Fraction(3 - (-1) ** value, 2)
        assert second == periodic_coefficient * first
        assert (1 if (value + 2) % 2 == 0 else 2) == second

    # A nonconstant 2-step rational coboundary certificate exercises C/(2Z).
    divisor = {Fraction(-1, 3): -1, Fraction(5, 3): 1}
    assert not divisor_orbit_signature(divisor, translation_step=2)
    witness = integrate_shift_divisor(divisor, translation_step=2)
    for value in range(10, 20):
        left = evaluate_factored_divisor(divisor, value)
        right = evaluate_factored_divisor(witness, value + 2) / (
            evaluate_factored_divisor(witness, value)
        )
        assert left == right


def rational_gamma_mshift_divisor(
    ray: RationalGammaRay, translation_step: int
) -> dict[Fraction, int]:
    divisor: defaultdict[Fraction, int] = defaultdict(int)
    for slope, offset, multiplicity in ray:
        increment = slope * translation_step
        assert increment.denominator == 1 and increment > 0
        for shift in range(int(increment)):
            divisor[-(offset + shift) / slope] += multiplicity
    return {root: value for root, value in divisor.items() if value}


def rational_gamma_mshift_quotient(
    ray: RationalGammaRay, translation_step: int, value: int
) -> Fraction:
    result = Fraction(1)
    for slope, offset, multiplicity in ray:
        increment = slope * translation_step
        assert increment.denominator == 1 and increment > 0
        for shift in range(int(increment)):
            result *= (slope * value + offset + shift) ** multiplicity
    return result


def rational_gamma_mshift_leading_constant(
    ray: RationalGammaRay, translation_step: int
) -> Fraction:
    result = Fraction(1)
    for slope, _offset, multiplicity in ray:
        increment = slope * translation_step
        assert increment.denominator == 1 and increment > 0
        result *= slope ** (int(increment) * multiplicity)
    return result


def verify_rational_slope_mfold_gamma() -> None:
    """Check duplication after passing to a common rational-slope period."""

    translation_step = 2
    left: RationalGammaRay = ((Fraction(1), Fraction(1), 1),)
    right: RationalGammaRay = (
        (Fraction(1, 2), Fraction(1, 2), 1),
        (Fraction(1, 2), Fraction(1), 1),
    )
    left_divisor = rational_gamma_mshift_divisor(left, translation_step)
    right_divisor = rational_gamma_mshift_divisor(right, translation_step)
    difference = subtract_divisors(left_divisor, right_divisor)
    assert not divisor_orbit_signature(difference, translation_step)
    witness = integrate_shift_divisor(difference, translation_step)
    assert not witness
    constant = rational_gamma_mshift_leading_constant(
        left, translation_step
    ) / rational_gamma_mshift_leading_constant(right, translation_step)
    assert constant == 4
    for value in range(1, 20):
        assert rational_gamma_mshift_quotient(
            left, translation_step, value
        ) == constant * rational_gamma_mshift_quotient(
            right, translation_step, value
        )


def parse_gamma_ray(specification: str) -> GammaRay:
    """Parse ``slope:offset[:multiplicity],...`` using exact fractions."""

    factors: list[GammaFactor] = []
    for raw_factor in specification.split(","):
        fields = raw_factor.strip().split(":")
        if len(fields) not in (2, 3):
            raise ValueError(
                "gamma factors must have slope:offset[:multiplicity] form"
            )
        slope = int(fields[0])
        offset = Fraction(fields[1])
        multiplicity = int(fields[2]) if len(fields) == 3 else 1
        if slope <= 0 or multiplicity == 0:
            raise ValueError("slopes must be positive and multiplicities nonzero")
        factors.append((slope, offset, multiplicity))
    if not factors:
        raise ValueError("a gamma ray must contain at least one factor")
    return tuple(factors)


def format_signature(signature: dict[Fraction, int]) -> str:
    if not signature:
        return "{}"
    entries = ", ".join(
        f"{residue}:{multiplicity:+d}"
        for residue, multiplicity in sorted(signature.items())
    )
    return "{" + entries + "}"


def format_coboundary(witness: dict[Fraction, int]) -> str:
    if not witness:
        return "1"
    return " * ".join(
        f"(n-({root}))^({multiplicity})"
        for root, multiplicity in sorted(witness.items())
    )


def compare_gamma_ray_specs(left_spec: str, right_spec: str) -> None:
    left = parse_gamma_ray(left_spec)
    right = parse_gamma_ray(right_spec)
    left_signature = gamma_orbit_signature(left)
    right_signature = gamma_orbit_signature(right)
    print(f"LEFT_SIGNATURE {format_signature(left_signature)}")
    print(f"RIGHT_SIGNATURE {format_signature(right_signature)}")
    if left_signature != right_signature:
        difference: defaultdict[Fraction, int] = defaultdict(int, left_signature)
        for residue, multiplicity in right_signature.items():
            difference[residue] -= multiplicity
        print(
            "DISTINCT_SIGNATURE "
            + format_signature(
                {key: value for key, value in difference.items() if value}
            )
        )
        return

    constant, witness = certify_gamma_equivalence(left, right)
    print("EQUIVALENT_SIGNATURE")
    print(f"EXPONENTIAL_BASE {constant}")
    print(f"RATIONAL_COBOUNDARY {format_coboundary(witness)}")
    print(
        "SEQUENCE_RATIO nonzero_constant * "
        f"({constant})^n * {format_coboundary(witness)}"
    )


def run_sympy_audit() -> None:
    """Replay every configured equivalence as a symbolic rational identity."""

    try:
        import sympy as sp
    except ImportError as error:
        raise SystemExit(
            "--sympy-audit requires the repository virtual environment"
        ) from error

    variable = sp.symbols("n")

    def sympy_fraction(value: Fraction):
        return sp.Rational(value.numerator, value.denominator)

    def symbolic_shift_quotient(ray: GammaRay):
        result = sp.S.One
        for slope, offset, multiplicity in ray:
            exact_offset = sympy_fraction(offset)
            for shift in range(slope):
                result *= (slope * variable + exact_offset + shift) ** multiplicity
        return sp.cancel(result)

    def symbolic_coboundary(witness: dict[Fraction, int]):
        result = sp.S.One
        for root, multiplicity in witness.items():
            result *= (variable - sympy_fraction(root)) ** multiplicity
        return sp.cancel(result)

    def audit_pair(left: GammaRay, right: GammaRay) -> None:
        constant, witness = certify_gamma_equivalence(left, right)
        rational_witness = symbolic_coboundary(witness)
        difference = sp.cancel(
            symbolic_shift_quotient(left) / symbolic_shift_quotient(right)
            - sympy_fraction(constant)
            * rational_witness.subs(variable, variable + 1)
            / rational_witness
        )
        assert difference == 0

    gauss_certificates = 0
    for slope in range(2, 13):
        for divisor in range(2, slope + 1):
            if slope % divisor:
                continue
            for numerator in range(-3, 4):
                offset = Fraction(numerator, 7)
                left = ((slope, offset, 1),)
                right = tuple(
                    (slope // divisor, (offset + shift) / divisor, 1)
                    for shift in range(divisor)
                )
                shifted_right = tuple(
                    (new_slope, new_offset + (index % 3) - 1, multiplicity)
                    for index, (new_slope, new_offset, multiplicity) in enumerate(
                        right
                    )
                )
                audit_pair(left, right)
                audit_pair(left, shifted_right)
                gauss_certificates += 2

    atoms: tuple[GammaFactor, ...] = tuple(
        (slope, Fraction(numerator, GAMMA_OFFSET_DENOMINATOR), 1)
        for slope in range(1, GAMMA_MAX_SLOPE + 1)
        for numerator in range(GAMMA_OFFSET_DENOMINATOR)
    )
    representatives: dict[tuple[tuple[Fraction, int], ...], GammaRay] = {}
    census_certificates = 0
    for length in range(1, 4):
        for indices in combinations_with_replacement(range(len(atoms)), length):
            ray = tuple(atoms[index] for index in indices)
            signature = gamma_signature_key(ray)
            if signature in representatives:
                audit_pair(ray, representatives[signature])
                census_certificates += 1
            else:
                representatives[signature] = ray

    assert gauss_certificates == 322
    assert census_certificates == 1_384
    print(
        "PASS SymPy rational-identity replay: "
        f"{gauss_certificates} Gauss/shift and "
        f"{census_certificates} census certificates"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify or compare factorial/gamma-affine trace rays."
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("LEFT", "RIGHT"),
        help=(
            "compare rays encoded as slope:offset[:multiplicity],...; "
            "offsets are exact fractions"
        ),
    )
    parser.add_argument(
        "--sympy-audit",
        action="store_true",
        help="replay every configured equivalence with SymPy rational algebra",
    )
    args = parser.parse_args()
    if args.compare and args.sympy_audit:
        parser.error("--compare and --sympy-audit are mutually exclusive")
    if args.compare:
        compare_gamma_ray_specs(*args.compare)
        return
    if args.sympy_audit:
        run_sympy_audit()
        return

    signed = verify_signed_orbit_reconstruction()
    coboundaries = verify_coboundary_integration()
    verify_gamma_examples()
    refinements, shifted_refinements = verify_gauss_refinements()
    (
        gamma_rays,
        gamma_classes,
        gamma_collisions,
        nontrivial_gamma_classes,
        largest_gamma_class,
    ) = verify_gamma_signature_census()
    signed_gamma_cases = verify_signed_gamma_property_cases()
    shifts = verify_integer_affine_shifts()
    (
        affine_profiles,
        affine_classes,
        affine_collisions,
        nontrivial_affine_classes,
        largest_affine_class,
    ) = verify_exact_affine_boundary_census()
    verify_one_scale_collision()
    verify_stirling_hierarchy()
    verify_mfold_periodic_symmetry()
    verify_rational_slope_mfold_gamma()
    valuation_profiles = verify_frobenius_dilation_valuations()
    sic_families = verify_sic_moment_family_signatures()
    print(
        "PASS signed factorial shift-orbit reconstruction: "
        f"{signed} nonzero vectors on slopes 1..{MAX_SLOPE} "
        "with coefficients -2..2"
    )
    print(
        "PASS rational shift-coboundary integration: "
        f"{coboundaries} zero-sum divisors on one seven-point orbit"
    )
    print(
        "PASS rational-offset gamma signature census: "
        f"{gamma_rays} rays -> {gamma_classes} classes, "
        f"{gamma_collisions} collisions in {nontrivial_gamma_classes} "
        f"nontrivial classes (largest {largest_gamma_class})"
    )
    print(
        "PASS Gauss refinement and integer-shift certificates: "
        f"{refinements} refinements, "
        f"{shifted_refinements} nonconstant rational coboundaries"
    )
    print(
        "PASS signed gamma property cases: "
        f"{signed_gamma_cases} deterministic randomized transformations"
    )
    print(f"PASS integer-affine factorial canonicalization: {shifts} exact values")
    print(
        "PASS exact affine-factorial boundary-transfer census: "
        f"{affine_profiles} profiles -> {affine_classes} classes, "
        f"{affine_collisions} collisions in {nontrivial_affine_classes} "
        f"nontrivial classes (largest {largest_affine_class}); every "
        "collision decomposes into elementary boundary squares"
    )
    print("PASS one-scale factorial collision separates on the second scale")
    print("PASS entropy collision separates at the first inverse-power moment")
    print(
        "PASS m-fold periodic symmetry, two-step coboundary, and rational-slope "
        "duplication"
    )
    print(
        "PASS Frobenius-dilation valuation reduction: "
        f"{valuation_profiles} signed profiles for p=2,3,5"
    )
    print(
        "PASS SIC radial-moment signature separation: "
        f"{sic_families} distinct (d,r) families through d=48"
    )
    print(
        "STATUS: complete rational-offset gamma signatures and exact "
        "integer-offset projective classes in characteristic zero; p-adic "
        "valuations require primitive-slope and unit/carry data"
    )


if __name__ == "__main__":
    main()
