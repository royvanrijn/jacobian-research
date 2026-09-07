#!/usr/bin/env sage
"""Compare the 074d9 record fibres at the local 2-descent/Kummer level.

The exact displayed quotient bases are imported from the certified 074d9
lineage and alternate-Q80 curve-12 artifacts.  For every quotient generator
this replay uses the standard cubic Kummer representative

    4*x(P) - zeta,

where zeta satisfies the monic integral completed-square 2-division cubic.
It records anonymous local-factor invariants rather than comparing power-basis
coordinates between different cubic etale algebras.

The finite places are 2, every bad prime of the curve, and the fixed common
good-prime block 53,67,71,79,83,97,101,113.  At each place the certificate
contains valuation parity, exact ambient and displayed-block Kummer
dimensions, component-image orders, and componentwise Hilbert symbols.  The
product of the component symbols is checked to be the trivial local Tate
pairing, as it must be for two points in the local Kummer image.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from run_fermigier_rank20_auxiliary_fingerprints import (  # noqa: E402
    f2_rank,
    qpari,
    two_adic_coords,
)

from sage.all import (  # noqa: E402
    EllipticCurve,
    GF,
    Matrix,
    PolynomialRing,
    QQ,
    ZZ,
    divisors,
    inverse_mod,
    pari,
)
from sage.version import version as sage_version  # noqa: E402


PROTOCOL = "R17074D9LOCALKUMMER"
SCHEMA = "elkies-k3.r17-074d9-quotient-arithmetic-blocks.v1"
LEGACY_SCHEMA = "elkies-k3.r17-074d9-local-kummer-meet.v1"
COMMON_GOOD_PRIMES = (53, 67, 71, 79, 83, 97, 101, 113)
TARGET_IDS = (351, 356, 376, 377, 385, 12)
RECORD_IDS = (356, 385)

PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
LINEAGE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
CURVE12 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
RIGID = ROOT / "artifacts/generated-results/elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
LEGACY_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-074d9-local-kummer-meet-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-074d9-quotient-arithmetic-blocks-v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def counted_multiset(values: list[Any]) -> list[dict[str, Any]]:
    counter = Counter(canonical_text(value) for value in values)
    return [
        {"count": counter[key], "value": json.loads(key)}
        for key in sorted(counter)
    ]


def signature_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def gf2_row_basis(rows: list[list[int]], width: int | None = None) -> list[list[int]]:
    """Return the deterministic RREF basis of a binary row space."""

    if width is None:
        width = len(rows[0]) if rows else 0
    if not rows:
        return []
    matrix = Matrix(GF(2), rows).echelon_form()
    return [
        [int(value) for value in row]
        for row in matrix.rows()
        if any(row)
    ]


def gf2_in_span(row: list[int], basis: list[list[int]]) -> bool:
    return f2_rank(basis + [row]) == f2_rank(basis)


def gf2_sum_rows(rows: list[list[int]], coefficients: list[int]) -> list[int]:
    if not rows:
        return []
    return [
        sum(int(coefficient) * int(row[index]) for coefficient, row in zip(coefficients, rows)) & 1
        for index in range(len(rows[0]))
    ]


def quotient_projection_data(
    labels: list[str], relation_rows: list[list[int]]
) -> dict[str, Any]:
    """Present F_2<labels>/<relation_rows> without choosing a complement."""

    relation_rref = gf2_row_basis(relation_rows, len(labels))
    matrix = Matrix(GF(2), relation_rref)
    pivots = list(matrix.pivots())
    nonpivots = [index for index in range(len(labels)) if index not in pivots]

    def project(row):
        reduced = [int(value) & 1 for value in row]
        for relation, pivot in zip(relation_rref, pivots):
            if reduced[pivot]:
                reduced = [left ^ right for left, right in zip(reduced, relation)]
        return [reduced[index] for index in nonpivots]

    named_images = []
    for index, label in enumerate(labels):
        row = [0] * len(labels)
        row[index] = 1
        named_images.append({"label": label, "quotient_coordinates": project(row)})
    return {
        "ambient_named_generators": labels,
        "relation_rref_over_f2": relation_rref,
        "relation_rank": len(relation_rref),
        "quotient_dimension": len(nonpivots),
        "canonical_presentation_basis_labels": [labels[index] for index in nonpivots],
        "canonical_presentation_basis_indices_zero_based": nonpivots,
        "named_generator_images": named_images,
        "basis_independence_statement": (
            "The quotient is the presented module on all twelve named generators. "
            "The ten nonpivot labels only give a deterministic coordinate presentation; "
            "all local targets kill the visible image and every Hilbert form is checked "
            "to annihilate the visible relation space."
        ),
    }


def quotient_local_module(
    image_rows: list[list[int]],
    quotient: dict[str, Any],
) -> dict[str, Any]:
    """Compute im(V -> A)/im(R -> A) for Q=V/R."""

    relation_rows = quotient["relation_rref_over_f2"]
    nonpivots = quotient["canonical_presentation_basis_indices_zero_based"]
    target_width = len(image_rows[0]) if image_rows else 0
    killed_rows = [gf2_sum_rows(image_rows, row) for row in relation_rows]
    killed_basis = gf2_row_basis(killed_rows, target_width)
    quotient_basis_target_rows = [image_rows[index] for index in nonpivots]

    image_basis = []
    coordinate_rows = []
    for row in quotient_basis_target_rows:
        found = None
        for mask in range(1 << len(image_basis)):
            trial = list(row)
            coordinates = []
            for index, basis_row in enumerate(image_basis):
                bit = (mask >> index) & 1
                coordinates.append(bit)
                if bit:
                    trial = [left ^ right for left, right in zip(trial, basis_row)]
            if gf2_in_span(trial, killed_basis):
                found = coordinates
                break
        if found is None:
            image_basis.append(list(row))
            for previous in coordinate_rows:
                previous.append(0)
            found = [0] * len(image_basis)
            found[-1] = 1
        coordinate_rows.append(found)

    quotient_dimension = len(image_basis)
    kernel_vectors = []
    for mask in range(1 << len(nonpivots)):
        row = [
            (mask >> index) & 1 for index in range(len(nonpivots))
        ]
        target = gf2_sum_rows(quotient_basis_target_rows, row)
        if gf2_in_span(target, killed_basis):
            kernel_vectors.append(row)
    kernel_basis = gf2_row_basis(kernel_vectors, len(nonpivots))

    named_rows = []
    for named in quotient["named_generator_images"]:
        named_rows.append(
            {
                "label": named["label"],
                "image_coordinates": gf2_sum_rows(
                    coordinate_rows, named["quotient_coordinates"]
                ),
            }
        )
    return {
        "construction": "image(V)/image(visible_relation_space)",
        "full_displayed_image_dimension": f2_rank(image_rows),
        "visible_image_dimension": len(killed_basis),
        "quotient_image_dimension": quotient_dimension,
        "canonical_quotient_basis_image_rows": coordinate_rows,
        "named_quotient_generator_images": named_rows,
        "kernel_in_canonical_quotient_coordinates": kernel_basis,
        "rank_nullity_verified": len(kernel_basis) + quotient_dimension
        == quotient["quotient_dimension"],
    }


def descended_pairing_matrix(
    matrix_rows: list[list[int]], quotient: dict[str, Any]
) -> dict[str, Any]:
    """Descend a binary bilinear form through the visible quotient."""

    matrix = Matrix(GF(2), matrix_rows)
    relation_matrix = Matrix(GF(2), quotient["relation_rref_over_f2"])
    radical_product = relation_matrix * matrix
    descends = radical_product.is_zero()
    indices = quotient["canonical_presentation_basis_indices_zero_based"]
    descended = [
        [int(matrix[left, right]) for right in indices]
        for left in indices
    ]
    return {
        "visible_relation_space_is_in_radical": descends,
        "visible_relation_pairing_obstruction_rows": [
            [int(value) for value in row] for row in radical_product.rows()
        ],
        "matrix_on_canonical_complement_presentation": descended,
        "quotient_pairing_matrix": descended if descends else None,
        "canonical_complement_matrix_rank": int(Matrix(GF(2), descended).rank()),
    }


def factor_discriminant(discriminant: ZZ, primes: list[ZZ]) -> tuple[list[dict[str, Any]], ZZ]:
    remainder = ZZ(discriminant)
    factors = []
    for prime in primes:
        exponent = 0
        while remainder % prime == 0:
            remainder //= prime
            exponent += 1
        if exponent == 0:
            raise ArithmeticError(f"declared bad prime {prime} misses the cubic discriminant")
        factors.append(
            {
                "prime": str(prime),
                "exponent": exponent,
                "primality_proved": bool(prime.is_prime(proof=True)),
            }
        )
    return factors, remainder


def kodaira_name(code: int) -> str:
    if code == 1:
        return "I0"
    names = {
        2: "II",
        3: "III",
        4: "IV",
        -1: "I0*",
        -2: "IV*",
        -3: "III*",
        -4: "II*",
    }
    if code in names:
        return names[code]
    if code >= 5:
        return f"I{code - 4}"
    if code <= -5:
        return f"I{-code - 4}*"
    raise ArithmeticError(f"unknown PARI Kodaira code {code}")


def reduction_kind(code: int) -> str:
    if code == 1:
        return "good"
    if code >= 5:
        return "multiplicative"
    return "additive"


def rational_valuation(value: QQ, prime: ZZ) -> int:
    value = QQ(value)
    return int(ZZ(value.numerator()).valuation(prime) - ZZ(value.denominator()).valuation(prime))


def point_in_E0(point, prime: ZZ, ainvs: tuple[QQ, ...]) -> bool:
    """Test membership in E_0(Q_p) on a minimal integral model."""

    if point.is_zero():
        return True
    x_coordinate, y_coordinate = map(QQ, point[:2])
    x_valuation = rational_valuation(x_coordinate, prime)
    y_valuation = rational_valuation(y_coordinate, prime)
    if x_valuation < 0 or y_valuation < 0:
        if x_valuation < 0 and y_valuation < 0:
            return True  # reduction to the smooth point at infinity
        raise ArithmeticError("an integral minimal-model point has inconsistent valuations")

    p = ZZ(prime)
    x_residue = ZZ(x_coordinate.numerator()) * inverse_mod(ZZ(x_coordinate.denominator()), p) % p
    y_residue = ZZ(y_coordinate.numerator()) * inverse_mod(ZZ(y_coordinate.denominator()), p) % p
    a1, a2, a3, a4, _a6 = [ZZ(value) % p for value in ainvs]
    derivative_x = (a1 * y_residue - 3 * x_residue**2 - 2 * a2 * x_residue - a4) % p
    derivative_y = (2 * y_residue + a1 * x_residue + a3) % p
    return derivative_x != 0 or derivative_y != 0


def component_image_order(point, prime: ZZ, tamagawa: int, ainvs: tuple[QQ, ...]) -> int:
    for candidate in divisors(ZZ(tamagawa)):
        if point_in_E0(candidate * point, prime, ainvs):
            return int(candidate)
    raise ArithmeticError("point component order does not divide the Tamagawa number")


def component_subgroup(generators, prime: ZZ, ainvs: tuple[QQ, ...], bound: int):
    """Enumerate the subgroup generated in E(Q_p)/E_0(Q_p)."""

    zero = generators[0].curve()(0) if generators else None
    representatives = [zero]
    changed = True
    while changed:
        changed = False
        for representative in list(representatives):
            for generator in generators:
                candidate = representative + generator
                if not any(
                    point_in_E0(candidate - known, prime, ainvs)
                    for known in representatives
                ):
                    representatives.append(candidate)
                    changed = True
                    if len(representatives) > bound:
                        raise ArithmeticError("displayed component subgroup exceeds Tamagawa number")
    return representatives


def component_member(point, subgroup, prime: ZZ, ainvs: tuple[QQ, ...]) -> bool:
    return any(point_in_E0(point - known, prime, ainvs) for known in subgroup)


def quotient_component_structure(
    *,
    points,
    prime: ZZ,
    tamagawa: int,
    ainvs: tuple[QQ, ...],
    quotient: dict[str, Any],
    visible_exact_rows: list[list[int]],
) -> dict[str, Any]:
    """Exact component quotient plus its mod-2 module induced from Q."""

    zero = points[0].curve()(0)

    def point_combination(coefficients):
        result = zero
        for coefficient, point in zip(coefficients, points):
            result += ZZ(coefficient) * point
        return result

    displayed_group = component_subgroup(points, prime, ainvs, tamagawa)
    visible_points = [point_combination(row) for row in visible_exact_rows]
    visible_group = component_subgroup(visible_points, prime, ainvs, tamagawa)
    if not all(component_member(point, displayed_group, prime, ainvs) for point in visible_group):
        raise ArithmeticError("visible component image escaped the displayed component subgroup")

    quotient_representatives = []
    for representative in displayed_group:
        if not any(
            component_member(
                representative - known, visible_group, prime, ainvs
            )
            for known in quotient_representatives
        ):
            quotient_representatives.append(representative)
    quotient_orders = []
    for representative in quotient_representatives:
        order = None
        for candidate in divisors(ZZ(tamagawa)):
            if component_member(candidate * representative, visible_group, prime, ainvs):
                order = int(candidate)
                break
        if order is None:
            raise ArithmeticError("component quotient order does not divide Tamagawa number")
        quotient_orders.append(order)
    quotient_group_order = len(displayed_group) // len(visible_group)
    if len(displayed_group) % len(visible_group):
        raise ArithmeticError("visible component image order does not divide displayed image order")
    if len(quotient_representatives) != quotient_group_order:
        raise ArithmeticError("component quotient coset enumeration has the wrong order")
    order_histogram = counted_multiset(quotient_orders)
    exponent = max(quotient_orders, default=1)
    filtration = [
        {
            "annihilator": int(candidate),
            "element_count": sum(order <= candidate and candidate % order == 0 for order in quotient_orders),
        }
        for candidate in divisors(ZZ(exponent))
    ]

    doubled_points = [2 * point for point in points]
    killed_mod2_group = component_subgroup(
        visible_points + doubled_points, prime, ainvs, tamagawa
    )
    nonpivots = quotient["canonical_presentation_basis_indices_zero_based"]
    basis_points = [points[index] for index in nonpivots]
    module_basis = []
    coordinate_rows = []
    for point in basis_points:
        found = None
        for mask in range(1 << len(module_basis)):
            trial = point
            coordinates = []
            for index, basis_point in enumerate(module_basis):
                bit = (mask >> index) & 1
                coordinates.append(bit)
                if bit:
                    trial -= basis_point
            if component_member(trial, killed_mod2_group, prime, ainvs):
                found = coordinates
                break
        if found is None:
            module_basis.append(point)
            for previous in coordinate_rows:
                previous.append(0)
            found = [0] * len(module_basis)
            found[-1] = 1
        coordinate_rows.append(found)

    kernel_vectors = []
    for mask in range(1 << len(basis_points)):
        coefficients = [(mask >> index) & 1 for index in range(len(basis_points))]
        if component_member(
            point_combination(
                [
                    coefficients[nonpivots.index(index)] if index in nonpivots else 0
                    for index in range(len(points))
                ]
            ),
            killed_mod2_group,
            prime,
            ainvs,
        ):
            kernel_vectors.append(coefficients)
    kernel_basis = gf2_row_basis(kernel_vectors, len(basis_points))
    named_rows = []
    for named in quotient["named_generator_images"]:
        named_rows.append(
            {
                "label": named["label"],
                "image_coordinates": gf2_sum_rows(
                    coordinate_rows, named["quotient_coordinates"]
                ),
            }
        )

    canonical_orders = []
    for point in basis_points:
        order = next(
            int(candidate)
            for candidate in divisors(ZZ(tamagawa))
            if component_member(candidate * point, visible_group, prime, ainvs)
        )
        canonical_orders.append(order)
    return {
        "exact_displayed_component_image_order": len(displayed_group),
        "exact_visible_component_image_order": len(visible_group),
        "exact_quotient_component_image_order": quotient_group_order,
        "exact_quotient_component_image_exponent": exponent,
        "exact_quotient_component_order_histogram": order_histogram,
        "exact_quotient_component_annihilator_filtration": filtration,
        "canonical_quotient_basis_component_orders_modulo_visible_image": canonical_orders,
        "mod2_component_module": {
            "construction": "image(displayed lattice)/(image(visible lattice)+2*image(displayed lattice))",
            "dimension": len(module_basis),
            "canonical_quotient_basis_image_rows": coordinate_rows,
            "named_quotient_generator_images": named_rows,
            "kernel_in_canonical_quotient_coordinates": kernel_basis,
            "rank_nullity_verified": len(kernel_basis) + len(module_basis)
            == quotient["quotient_dimension"],
        },
    }


def local_places(nf, prime: int):
    places = list(pari.idealprimedec(nf, prime))
    descriptors = []
    for place in places:
        ramification_index = int(place[2])
        residue_degree = int(place[3])
        descriptors.append(
            {
                "ramification_index": ramification_index,
                "residue_degree": residue_degree,
                "local_degree": ramification_index * residue_degree,
            }
        )
    if sum(item["local_degree"] for item in descriptors) != 3:
        raise ArithmeticError(f"local degrees above {prime} do not sum to three")
    return places, descriptors


def odd_squareclass_rows(nf, alphas, places):
    prepared = []
    for place in places:
        uniformizer_column = pari.idealappr(nf, place)
        uniformizer = pari.nfbasistoalg(nf, uniformizer_column)
        if int(pari.idealval(nf, uniformizer, place)) != 1:
            raise ArithmeticError("PARI returned a non-uniformizer")
        reduction = pari.nfmodprinit(nf, place)
        prepared.append((place, uniformizer, reduction))

    rows = []
    for alpha in alphas:
        row = []
        for place, uniformizer, reduction in prepared:
            valuation = int(pari.idealval(nf, alpha, place))
            unit = alpha / uniformizer**valuation
            residue = pari.nfmodpr(nf, unit, reduction)
            row.extend((valuation & 1, 0 if bool(pari.issquare(residue)) else 1))
        rows.append(row)
    return rows


def anonymous_histogram(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return counted_multiset(entries)


def place_symbol_multiset(symbols: list[int], descriptors: list[dict[str, Any]]) -> list[dict[str, Any]]:
    entries = []
    for symbol, descriptor in zip(symbols, descriptors):
        entries.append(
            {
                "ramification_index": descriptor["ramification_index"],
                "residue_degree": descriptor["residue_degree"],
                "hilbert_symbol": int(symbol),
            }
        )
    return anonymous_histogram(entries)


def direction_squareclass_multiset(row: list[int], descriptors: list[dict[str, Any]]) -> list[dict[str, Any]]:
    entries = []
    for index, descriptor in enumerate(descriptors):
        entries.append(
            {
                "ramification_index": descriptor["ramification_index"],
                "residue_degree": descriptor["residue_degree"],
                "valuation_parity": int(row[2 * index]),
                "unit_nonsquare": int(row[2 * index + 1]),
            }
        )
    return anonymous_histogram(entries)


def valuation_support(alpha, nf, places, descriptors):
    entries = []
    for place, descriptor in zip(places, descriptors):
        if int(pari.idealval(nf, alpha, place)) & 1:
            entries.append(
                {
                    "ramification_index": descriptor["ramification_index"],
                    "residue_degree": descriptor["residue_degree"],
                }
            )
    return anonymous_histogram(entries)


def selected_signature(local_record: dict[str, Any], selected_labels: list[str]) -> dict[str, Any]:
    selected = set(selected_labels)
    directions = [
        row for row in local_record["directions"] if row["label"] in selected
    ]
    pairs = [
        row
        for row in local_record["pairwise_hilbert_symbols"]
        if row["left"] in selected and row["right"] in selected
    ]
    direction_values = []
    for row in directions:
        direction_values.append(
            {
                "component_image_order": row["component_group_image_order"],
                "valuation_parity_support": row["valuation_parity_support"],
                "odd_local_squareclass_multiset": row.get("odd_local_squareclass_multiset"),
                "self_component_hilbert_symbol_multiset": row[
                    "self_component_hilbert_symbol_multiset"
                ],
            }
        )
    signature = {
        "reduction_kind": local_record["reduction_kind"],
        "kodaira_symbol": local_record["kodaira_symbol"],
        "tamagawa_number": local_record["tamagawa_number"],
        "ambient_local_kummer_dimension": local_record[
            "ambient_local_kummer_dimension"
        ],
        "selected_block_image_dimension": local_record[
            "selected_block_image_dimension"
        ],
        "direction_invariant_multiset": counted_multiset(direction_values),
        "pairwise_component_hilbert_multiset": counted_multiset(
            [row["component_hilbert_symbol_multiset"] for row in pairs]
        ),
    }
    signature["sha256"] = signature_hash(signature)
    return signature


def audit_local_place(
    *,
    curve,
    pari_curve,
    ainvs,
    nf,
    alphas,
    labels,
    points,
    prime,
    is_bad,
    selected_labels,
    quotient=None,
    visible_exact_rows=None,
):
    places, descriptors = local_places(nf, prime)
    local_reduction = pari.elllocalred(pari_curve, prime)
    conductor_exponent = int(local_reduction[0])
    kodaira_code = int(local_reduction[1])
    minimal_change = [str(value) for value in local_reduction[2]]
    tamagawa = int(local_reduction[3])
    if minimal_change != ["1", "0", "0", "0"]:
        raise ArithmeticError("the pinned public model is not locally minimal")
    if bool(is_bad) != (kodaira_code != 1):
        raise ArithmeticError(f"bad-prime classification changed at p={prime}")

    odd_rows = None
    if prime == 2:
        _basis, _origins, two_rows = two_adic_coords(pari, nf, places, alphas)
        local_image_rows = two_rows
        block_image_dimension = f2_rank(two_rows)
        selected_indices = [index for index, label in enumerate(labels) if label in selected_labels]
        _sb, _so, selected_local_image_rows = two_adic_coords(
            pari, nf, places, [alphas[index] for index in selected_indices]
        )
        selected_image_dimension = f2_rank(selected_local_image_rows)
    else:
        odd_rows = odd_squareclass_rows(nf, alphas, places)
        local_image_rows = odd_rows
        block_image_dimension = f2_rank(odd_rows)
        selected_local_image_rows = [
            row for row, label in zip(odd_rows, labels) if label in selected_labels
        ]
        selected_image_dimension = f2_rank(selected_local_image_rows)

    ambient_dimension = len(places) - 1 + int(prime == 2)
    if block_image_dimension > ambient_dimension or selected_image_dimension > ambient_dimension:
        raise ArithmeticError("known point images exceed E(Q_p)/2E(Q_p)")

    direction_records = []
    component_pairing_matrices = [
        [[0] * len(labels) for _ in labels] for _place in places
    ]
    for index, (label, alpha, point) in enumerate(zip(labels, alphas, points)):
        self_symbols = [int(pari.nfhilbert(nf, alpha, alpha, place)) for place in places]
        if any(symbol not in (-1, 1) for symbol in self_symbols):
            raise ArithmeticError("a local Hilbert symbol is not +/-1")
        if (-1) ** sum(symbol == -1 for symbol in self_symbols) != 1:
            raise ArithmeticError("self local Tate pairing is nontrivial")
        for place_index, symbol in enumerate(self_symbols):
            component_pairing_matrices[place_index][index][index] = int(symbol == -1)
        self_multiset = place_symbol_multiset(self_symbols, descriptors)
        direction = {
            "label": label,
            "valuation_parity_support": valuation_support(
                alpha, nf, places, descriptors
            ),
            "component_group_image_order": component_image_order(
                point, ZZ(prime), tamagawa, ainvs
            ),
            "self_component_hilbert_symbol_multiset": self_multiset,
            "self_corestricted_local_tate_symbol": 1,
        }
        if odd_rows is not None:
            direction["odd_local_squareclass_multiset"] = direction_squareclass_multiset(
                odd_rows[index], descriptors
            )
        direction_records.append(direction)

    pairing_records = []
    for left_index in range(len(labels)):
        for right_index in range(left_index + 1, len(labels)):
            symbols = [
                int(
                    pari.nfhilbert(
                        nf,
                        alphas[left_index],
                        alphas[right_index],
                        place,
                    )
                )
                for place in places
            ]
            negative_count = sum(symbol == -1 for symbol in symbols)
            if any(symbol not in (-1, 1) for symbol in symbols):
                raise ArithmeticError("a local Hilbert symbol is not +/-1")
            corestricted = (-1) ** negative_count
            if corestricted != 1:
                raise ArithmeticError("point Kummer images are not locally isotropic")
            for place_index, symbol in enumerate(symbols):
                bit = int(symbol == -1)
                component_pairing_matrices[place_index][left_index][right_index] = bit
                component_pairing_matrices[place_index][right_index][left_index] = bit
            pairing_records.append(
                {
                    "left": labels[left_index],
                    "right": labels[right_index],
                    "component_hilbert_symbol_multiset": place_symbol_multiset(
                        symbols, descriptors
                    ),
                    "corestricted_local_tate_symbol": corestricted,
                }
            )

    record = {
        "rational_prime": str(prime),
        "place_kind": "bad" if is_bad else "fixed_common_good",
        "reduction_kind": reduction_kind(kodaira_code),
        "kodaira_symbol": kodaira_name(kodaira_code),
        "conductor_exponent": conductor_exponent,
        "minimal_discriminant_valuation": rational_valuation(
            QQ(curve.discriminant()), ZZ(prime)
        ),
        "tamagawa_number": tamagawa,
        "local_factor_descriptors": descriptors,
        "ambient_local_kummer_dimension": ambient_dimension,
        "quotient_basis_image_dimension": block_image_dimension,
        "selected_block_labels": selected_labels,
        "selected_block_image_dimension": selected_image_dimension,
        "known_residual_localization": quotient_local_module(
            local_image_rows, quotient_projection_data(labels, [])
        ),
        "selected_block_localization": quotient_local_module(
            selected_local_image_rows,
            quotient_projection_data(selected_labels, []),
        ),
        "directions": direction_records,
        "pairwise_hilbert_symbols": pairing_records,
        "all_pairwise_corestricted_local_tate_symbols_trivial": True,
    }
    if quotient is not None:
        if visible_exact_rows is None:
            raise ArithmeticError("quotient audit is missing exact visible rows")
        if prime == 2:
            factor_image_rows = []
            for place in places:
                _fb, _fo, rows = two_adic_coords(pari, nf, [place], alphas)
                factor_image_rows.append(rows)
        else:
            factor_image_rows = [
                [row[2 * index : 2 * index + 2] for row in odd_rows]
                for index in range(len(places))
            ]

        factor_records = []
        for descriptor, rows, matrix_rows in zip(
            descriptors, factor_image_rows, component_pairing_matrices
        ):
            descended = descended_pairing_matrix(matrix_rows, quotient)
            factor_records.append(
                {
                    "ramification_index": descriptor["ramification_index"],
                    "residue_degree": descriptor["residue_degree"],
                    "local_degree": descriptor["local_degree"],
                    "quotient_local_kummer_module": quotient_local_module(rows, quotient),
                    "hilbert_pairing": {
                        "encoding": "1 means Hilbert symbol -1; 0 means +1",
                        "matrix_on_twelve_named_generators": matrix_rows,
                        **descended,
                    },
                }
            )
        factor_records.sort(key=canonical_text)
        for index, factor_record in enumerate(factor_records, start=1):
            factor_record["canonical_anonymous_factor_id"] = f"F{index}"

        total_pairing = Matrix(GF(2), component_pairing_matrices[0])
        for matrix_rows in component_pairing_matrices[1:]:
            total_pairing += Matrix(GF(2), matrix_rows)
        if not total_pairing.is_zero():
            raise ArithmeticError("componentwise Hilbert matrices do not corestrict to zero")
        record["visible_quotient_arithmetic"] = {
            "quotient_local_kummer_module": quotient_local_module(
                local_image_rows, quotient
            ),
            "component_group_filtration": quotient_component_structure(
                points=points,
                prime=ZZ(prime),
                tamagawa=tamagawa,
                ainvs=ainvs,
                quotient=quotient,
                visible_exact_rows=visible_exact_rows,
            ),
            "anonymous_local_factors": factor_records,
            "all_componentwise_hilbert_forms_descend": all(
                factor["hilbert_pairing"][
                    "visible_relation_space_is_in_radical"
                ]
                for factor in factor_records
            ),
            "corestricted_hilbert_form_is_zero": True,
        }
    record["selected_block_canonical_signature"] = selected_signature(
        record, selected_labels
    )
    return record


def load_inputs():
    public = json.loads(PUBLIC.read_text())
    lineage = json.loads(LINEAGE.read_text())
    curve12 = json.loads(CURVE12.read_text())
    rigid = json.loads(RIGID.read_text())
    if public.get("status") != "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES":
        raise ArithmeticError("public ICARM projection is not certified")
    if lineage.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ArithmeticError("074d9 quotient source is not certified")
    if curve12.get("status") != "PROVED_CURVE12_NATIVE_ALTERNATE_Q80_AND_DISPLAYED_QUOTIENT":
        raise ArithmeticError("curve-12 quotient source is not certified")
    if rigid.get("status") != "PASS_EXACT_COMPLETE_074D9_CROSS_FIBRE_BISECTION_TRANSFER":
        raise ArithmeticError("rigid-transfer source is not certified")
    return public, lineage, curve12, rigid


def quotient_bases(lineage, curve12):
    bases = {
        int(row["curve_id"]): list(row["free_basis_modulo_generic_group"])
        for row in lineage["exceptional_quotients"]
    }
    bases[12] = list(
        curve12["displayed_exceptional_quotient"][
            "free_basis_modulo_specialized_generic"
        ]
    )
    expected = {
        351: [f"P{index}" for index in range(18, 26)],
        356: [f"P{index}" for index in range(18, 30)],
        376: [f"P{index}" for index in range(18, 23)],
        377: [f"P{index}" for index in range(18, 24)],
        385: [f"P{index}" for index in range(18, 30)],
        12: ["P2", "P11", "P4", "P3", "P6", "P8", "P17", "P10", "P28", "P24", "P19", "P15"],
    }
    if bases != expected:
        raise ArithmeticError("a certified displayed quotient basis changed")
    return bases


def rigid_complements(rigid, bases):
    complements = {}
    provenance = {}
    quotient_presentations = {}
    records = {int(row["curve_id"]): row for row in rigid["fibres"]}
    for curve_id in RECORD_IDS:
        span_rows = [
            row["finite_quotient_class_modulo_generic_17"][
                "displayed_quotient_coordinates_over_f2"
            ]
            for row in records[curve_id]["records"]
        ]
        exact_rows = [
            row["exact_displayed_free_quotient_class"]["coordinates"]
            for row in records[curve_id]["records"]
        ]
        matrix = Matrix(GF(2), span_rows)
        if matrix.rank() != 2:
            raise ArithmeticError("record rigid span no longer has rank two")
        pivots = list(matrix.echelon_form().pivots())
        nonpivots = [index for index in range(12) if index not in pivots]
        complements[curve_id] = [bases[curve_id][index] for index in nonpivots]
        quotient_presentations[curve_id] = quotient_projection_data(
            bases[curve_id], span_rows
        )
        quotient_presentations[curve_id]["visible_exact_lattice_generators"] = exact_rows
        provenance[curve_id] = {
            "quotient_basis": bases[curve_id],
            "rigid_span_vectors_over_f2": span_rows,
            "rigid_span_exact_integral_vectors": exact_rows,
            "canonical_rref_pivot_labels": [bases[curve_id][index] for index in pivots],
            "canonical_coordinate_complement_labels": complements[curve_id],
            "rigid_span_rank": 2,
            "complement_dimension": 10,
        }
    return complements, provenance, quotient_presentations


def audit_curve(
    public_record,
    labels,
    selected_labels,
    parameter,
    chart,
    quotient_presentation=None,
):
    curve_id = int(public_record["id"])
    ainvs = tuple(QQ(value) for value in public_record["ainvs"])
    curve = EllipticCurve(QQ, list(ainvs))
    pari_curve = pari.ellinit([qpari(pari, value) for value in ainvs])
    points_by_label = {
        f"P{index + 1}": curve(QQ(point[0]), QQ(point[1]))
        for index, point in enumerate(public_record["points"])
    }
    points = [points_by_label[label] for label in labels]

    a1, a2, a3, a4, a6 = ainvs
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    if not all(value.denominator() == 1 for value in (b2, b4, b6)):
        raise ArithmeticError("the pinned global minimal model is not integral")
    polynomial_ring = PolynomialRing(ZZ, "z")
    z = polynomial_ring.gen()
    polynomial = z**3 + ZZ(b2) * z**2 + ZZ(8 * b4) * z + ZZ(16 * b6)
    if not polynomial.is_irreducible():
        raise ArithmeticError("the 2-division cubic is reducible over Q")
    discriminant = abs(ZZ(polynomial.discriminant()))
    if discriminant != 256 * abs(ZZ(curve.discriminant())):
        raise ArithmeticError("completed-square cubic discriminant identity failed")

    bad_primes = [ZZ(value) for value in public_record["bad_primes"]]
    factorization, remainder = factor_discriminant(discriminant, bad_primes)
    if remainder != 1 or not all(row["primality_proved"] for row in factorization):
        raise ArithmeticError("the bad-prime list is not a proved complete factorization")
    if any(ZZ(prime) in bad_primes for prime in COMMON_GOOD_PRIMES):
        raise ArithmeticError("the fixed common-good block contains a bad prime")

    pari.addprimes(bad_primes)
    nf = pari.nfinit([pari(polynomial), bad_primes])
    if list(pari.nfcertify(nf)):
        raise ArithmeticError("the factor-supplied cubic number field failed certification")
    theta = pari(f"Mod(z,{polynomial})")

    alphas = []
    kummer_records = []
    for label, point in zip(labels, points):
        x_coordinate, y_coordinate = map(QQ, point[:2])
        alpha = qpari(pari, 4 * x_coordinate) - theta
        norm = QQ(str(pari.nfeltnorm(nf, alpha)))
        square_root = 4 * (2 * y_coordinate + a1 * x_coordinate + a3)
        if norm != square_root**2:
            raise ArithmeticError("Kummer representative norm is not the expected square")
        alphas.append(alpha)
        kummer_records.append(
            {
                "label": label,
                "representative_coefficients_in_1_zeta_zeta2": [
                    str(4 * x_coordinate),
                    "-1",
                    "0",
                ],
                "norm": str(norm),
                "canonical_norm_square_root": str(square_root),
            }
        )

    finite_primes = sorted(set(int(prime) for prime in bad_primes) | set(COMMON_GOOD_PRIMES))
    local_records = []
    for prime in finite_primes:
        local_records.append(
            audit_local_place(
                curve=curve,
                pari_curve=pari_curve,
                ainvs=ainvs,
                nf=nf,
                alphas=alphas,
                labels=labels,
                points=points,
                prime=prime,
                is_bad=ZZ(prime) in bad_primes,
                selected_labels=selected_labels,
                quotient=quotient_presentation,
                visible_exact_rows=(
                    quotient_presentation["visible_exact_lattice_generators"]
                    if quotient_presentation is not None
                    else None
                ),
            )
        )
        print(
            f"{PROTOCOL}|curve={curve_id}|prime={prime}|"
            f"basis_dim={local_records[-1]['quotient_basis_image_dimension']}|"
            f"selected_dim={local_records[-1]['selected_block_image_dimension']}|status=PASS",
            flush=True,
        )

    local_index = {int(row["rational_prime"]): row for row in local_records}
    good_dimension_vector = [
        local_index[prime]["ambient_local_kummer_dimension"]
        for prime in COMMON_GOOD_PRIMES
    ]
    result = {
        "curve_id": curve_id,
        "role": "record" if curve_id in RECORD_IDS else "control",
        "native_chart": chart,
        "family_parameter": parameter,
        "global_minimal_model": [str(value) for value in ainvs],
        "bad_primes": [str(prime) for prime in bad_primes],
        "two_division_etale_algebra": {
            "generator": "zeta",
            "defining_polynomial": str(polynomial),
            "defining_polynomial_coefficients_ascending": [
                str(polynomial[index]) for index in range(4)
            ],
            "discriminant": str(discriminant),
            "proved_prime_factorization": factorization,
            "globally_irreducible_cubic_field": True,
            "kummer_map": "P maps to [4*x(P)-zeta]; the factor 4 is a global square",
        },
        "certified_displayed_exceptional_quotient_basis": labels,
        "selected_comparison_block": selected_labels,
        "kummer_images": kummer_records,
        "fixed_common_good_ambient_dimension_vector": good_dimension_vector,
        "local_places": local_records,
    }
    if quotient_presentation is not None:
        result["visible_quotient_presentation"] = quotient_presentation
    return result


def quotient_arithmetic_block_decomposition(curve_record: dict[str, Any]) -> dict[str, Any]:
    """Assemble the bad-place quotient tensor and its labelled direction graph."""

    quotient = curve_record["visible_quotient_presentation"]
    basis_labels = quotient["canonical_presentation_basis_labels"]
    dimension = quotient["quotient_dimension"]
    bad_places = sorted(
        [row for row in curve_record["local_places"] if row["place_kind"] == "bad"],
        key=lambda row: ZZ(row["rational_prime"]),
    )
    factor_coordinates = []
    factor_matrices = []
    pairing_descent_obstructions = []
    for place in bad_places:
        for factor in place["visible_quotient_arithmetic"]["anonymous_local_factors"]:
            factor_coordinates.append(
                {
                    "rational_prime": place["rational_prime"],
                    "local_factor": factor["canonical_anonymous_factor_id"],
                    "ramification_index": factor["ramification_index"],
                    "residue_degree": factor["residue_degree"],
                }
            )
            pairing = factor["hilbert_pairing"]
            if pairing["visible_relation_space_is_in_radical"]:
                factor_matrices.append(pairing["quotient_pairing_matrix"])
            else:
                pairing_descent_obstructions.append(
                    {
                        "rational_prime": place["rational_prime"],
                        "local_factor": factor["canonical_anonymous_factor_id"],
                        "ramification_index": factor["ramification_index"],
                        "residue_degree": factor["residue_degree"],
                        "visible_relation_pairing_obstruction_rows": pairing[
                            "visible_relation_pairing_obstruction_rows"
                        ],
                    }
                )

    bad_place_local_structures = [
        {
            "rational_prime": place["rational_prime"],
            "reduction_kind": place["reduction_kind"],
            "kodaira_symbol": place["kodaira_symbol"],
            "tamagawa_number": place["tamagawa_number"],
            **place["visible_quotient_arithmetic"],
        }
        for place in bad_places
    ]
    if pairing_descent_obstructions:
        quotient_compatible_signatures = []
        for place in bad_places:
            arithmetic = place["visible_quotient_arithmetic"]
            component = arithmetic["component_group_filtration"]
            quotient_compatible_signatures.append(
                [
                    {
                        "local_kummer_image": arithmetic[
                            "quotient_local_kummer_module"
                        ]["canonical_quotient_basis_image_rows"][direction],
                        "component_mod2_image": component["mod2_component_module"][
                            "canonical_quotient_basis_image_rows"
                        ][direction],
                        "component_order_modulo_visible": component[
                            "canonical_quotient_basis_component_orders_modulo_visible_image"
                        ][direction],
                    }
                    for direction in range(dimension)
                ]
            )

        def quotient_compatible_distinguished(place_indices):
            signatures = [
                canonical_text(
                    [
                        quotient_compatible_signatures[index][direction]
                        for index in place_indices
                    ]
                )
                for direction in range(dimension)
            ]
            return len(set(signatures)) == dimension

        compatible_minimum_sets = []
        for cardinality in range(len(bad_places) + 1):
            for place_indices in combinations(range(len(bad_places)), cardinality):
                if quotient_compatible_distinguished(place_indices):
                    compatible_minimum_sets.append(
                        [bad_places[index]["rational_prime"] for index in place_indices]
                    )
            if compatible_minimum_sets:
                break

        compatible_full_classes = {}
        for direction in range(dimension):
            key = canonical_text(
                [
                    quotient_compatible_signatures[index][direction]
                    for index in range(len(bad_places))
                ]
            )
            compatible_full_classes.setdefault(key, []).append(basis_labels[direction])

        obstruction_rows = [
            row
            for obstruction in pairing_descent_obstructions
            for row in obstruction["visible_relation_pairing_obstruction_rows"]
            if any(row)
        ]
        return {
            "status": "NOT_DEFINED_COMPONENT_HILBERT_FORMS_FAIL_TO_DESCEND",
            "scope": (
                "The F_2 quotient of the twelve displayed exceptional directions by "
                "the two exact rigid-visible directions; no full Mordell-Weil quotient is asserted."
            ),
            "quotient_presentation": quotient,
            "bad_place_local_structures": bad_place_local_structures,
            "pairing_descent_obstructions": pairing_descent_obstructions,
            "pairing_descent_obstruction_span": {
                "ambient_dual_coordinate_order": quotient["ambient_named_generators"],
                "dimension": f2_rank(obstruction_rows),
                "rref_basis": gf2_row_basis(
                    obstruction_rows,
                    len(quotient["ambient_named_generators"]),
                ),
            },
            "nontrivial_pairing_graph": {
                "status": "NOT_DEFINED_ON_QUOTIENT",
                "reason": (
                    "Changing a quotient lift by a visible relation changes at least "
                    "one componentwise Hilbert pairing."
                ),
            },
            "corestricted_local_tate_pairing_control": {
                "all_rational_bad_place_matrices_zero": True,
                "graph_edges": [],
                "formal_coordinate_component_partition": [1] * dimension,
                "boundary": (
                    "The genuine corestricted local Tate pairing is identically zero "
                    "on images of local points. Its edgeless graph is basis-independent "
                    "but contains no coupling information and selects no canonical "
                    "one-dimensional summands."
                ),
            },
            "pair_placewise_pairing_vector_span": {
                "status": "NOT_DEFINED_ON_QUOTIENT"
            },
            "minimum_bad_places_distinguishing_canonical_directions_without_hilbert_data": {
                "uses": [
                    "quotient local Kummer images",
                    "quotient mod-2 component images",
                    "exact component orders modulo the visible component image",
                ],
                "exists": bool(compatible_minimum_sets),
                "minimum_cardinality": (
                    len(compatible_minimum_sets[0]) if compatible_minimum_sets else None
                ),
                "lexicographically_first_set": (
                    compatible_minimum_sets[0] if compatible_minimum_sets else None
                ),
                "number_of_minimum_sets": len(compatible_minimum_sets),
                "all_minimum_sets": compatible_minimum_sets,
                "full_bad_place_indistinguishability_classes": sorted(
                    compatible_full_classes.values()
                ),
                "boundary": (
                    "This distinguishes the fixed canonical presentation directions "
                    "using only quotient-compatible data; it does not repair the "
                    "non-descending Hilbert tensor."
                ),
            },
            "minimum_bad_places_distinguishing_directions_with_full_requested_structure": {
                "status": "NOT_DEFINED_FOR_THE_REQUESTED_GLOBAL_STRUCTURE"
            },
            "indecomposable_components": {
                "status": "NOT_DEFINED_ON_QUOTIENT",
                "reason": (
                    "The vector-valued Hilbert tensor required to define coupling "
                    "does not descend through the two visible directions."
                ),
            },
            "basis_dependent_complement_partition_deliberately_not_reported": True,
        }

    pair_factor_vectors = []
    pair_place_vectors = []
    edge_evidence = {}
    for left in range(dimension):
        for right in range(left + 1, dimension):
            factor_vector = [matrix[left][right] for matrix in factor_matrices]
            place_vector = []
            offset = 0
            nontrivial_places = []
            for place in bad_places:
                factor_count = len(
                    place["visible_quotient_arithmetic"]["anonymous_local_factors"]
                )
                bit = int(any(factor_vector[offset : offset + factor_count]))
                place_vector.append(bit)
                if bit:
                    nontrivial_places.append(place["rational_prime"])
                offset += factor_count
            pair_factor_vectors.append(
                {
                    "left": basis_labels[left],
                    "right": basis_labels[right],
                    "vector": factor_vector,
                }
            )
            pair_place_vectors.append(
                {
                    "left": basis_labels[left],
                    "right": basis_labels[right],
                    "vector": place_vector,
                }
            )
            if nontrivial_places:
                edge_evidence[(left, right)] = nontrivial_places

    adjacency = [set() for _ in range(dimension)]
    edges = []
    for (left, right), primes in sorted(edge_evidence.items()):
        adjacency[left].add(right)
        adjacency[right].add(left)
        edges.append(
            {
                "left": basis_labels[left],
                "right": basis_labels[right],
                "nontrivial_rational_places": primes,
            }
        )
    components = []
    unseen = set(range(dimension))
    while unseen:
        start = min(unseen)
        stack = [start]
        component = []
        unseen.remove(start)
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbour in sorted(adjacency[vertex], reverse=True):
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))
    components.sort(key=lambda values: (-len(values), values))

    direction_vectors = []
    for direction in range(dimension):
        vector = []
        for matrix in factor_matrices:
            vector.extend(matrix[direction])
        direction_vectors.append(
            {"label": basis_labels[direction], "vector": vector}
        )

    def direction_signature(place, direction):
        arithmetic = place["visible_quotient_arithmetic"]
        component = arithmetic["component_group_filtration"]
        return {
            "local_kummer_image": arithmetic["quotient_local_kummer_module"][
                "canonical_quotient_basis_image_rows"
            ][direction],
            "component_mod2_image": component["mod2_component_module"][
                "canonical_quotient_basis_image_rows"
            ][direction],
            "component_order_modulo_visible": component[
                "canonical_quotient_basis_component_orders_modulo_visible_image"
            ][direction],
            "componentwise_hilbert_rows": [
                factor["hilbert_pairing"]["quotient_pairing_matrix"]
                [direction]
                for factor in arithmetic["anonymous_local_factors"]
            ],
        }

    local_direction_signatures = [
        [direction_signature(place, direction) for direction in range(dimension)]
        for place in bad_places
    ]

    def distinguished(place_indices):
        signatures = [
            canonical_text(
                [local_direction_signatures[index][direction] for index in place_indices]
            )
            for direction in range(dimension)
        ]
        return len(set(signatures)) == dimension

    minimum_sets = []
    for cardinality in range(len(bad_places) + 1):
        for place_indices in combinations(range(len(bad_places)), cardinality):
            if distinguished(place_indices):
                minimum_sets.append(
                    [bad_places[index]["rational_prime"] for index in place_indices]
                )
        if minimum_sets:
            break

    full_signature_classes = {}
    for direction in range(dimension):
        key = canonical_text(
            [
                local_direction_signatures[index][direction]
                for index in range(len(bad_places))
            ]
        )
        full_signature_classes.setdefault(key, []).append(basis_labels[direction])

    factor_vector_rows = [row["vector"] for row in pair_factor_vectors]
    place_vector_rows = [row["vector"] for row in pair_place_vectors]
    direction_vector_rows = [row["vector"] for row in direction_vectors]
    return {
        "scope": (
            "The F_2 quotient of the twelve displayed exceptional directions by "
            "the two exact rigid-visible directions; no full Mordell-Weil quotient is asserted."
        ),
        "quotient_presentation": quotient,
        "status": "PASS_QUOTIENT_ARITHMETIC_BLOCK_DECOMPOSITION",
        "bad_place_local_structures": bad_place_local_structures,
        "pairing_coordinate_order": factor_coordinates,
        "direction_pairing_vectors": direction_vectors,
        "direction_pairing_vector_span": {
            "dimension": f2_rank(direction_vector_rows),
            "rref_basis": gf2_row_basis(
                direction_vector_rows,
                len(direction_vector_rows[0]) if direction_vector_rows else 0,
            ),
        },
        "pair_placewise_pairing_vectors": pair_place_vectors,
        "pair_placewise_pairing_vector_span": {
            "coordinate_order": [place["rational_prime"] for place in bad_places],
            "dimension": f2_rank(place_vector_rows),
            "rref_basis": gf2_row_basis(
                place_vector_rows,
                len(place_vector_rows[0]) if place_vector_rows else 0,
            ),
        },
        "pair_local_factor_pairing_vectors": pair_factor_vectors,
        "pair_local_factor_pairing_vector_span": {
            "dimension": f2_rank(factor_vector_rows),
            "rref_basis": gf2_row_basis(
                factor_vector_rows,
                len(factor_vector_rows[0]) if factor_vector_rows else 0,
            ),
        },
        "nontrivial_pairing_graph": {
            "vertices": basis_labels,
            "edges": edges,
            "connected_components": [
                [basis_labels[index] for index in component] for component in components
            ],
            "dimension_partition": [len(component) for component in components],
            "lift_independence": (
                "Edges are unchanged when a displayed lift is altered by a visible "
                "relation, because every component Hilbert form annihilates that relation space."
            ),
            "coordinate_caveat": (
                "Connectivity is for the fixed named RREF presentation basis; it is not "
                "claimed invariant under an arbitrary GL(10,F_2) change of quotient basis."
            ),
        },
        "minimum_bad_places_distinguishing_canonical_directions": {
            "exists": bool(minimum_sets),
            "minimum_cardinality": len(minimum_sets[0]) if minimum_sets else None,
            "lexicographically_first_set": minimum_sets[0] if minimum_sets else None,
            "number_of_minimum_sets": len(minimum_sets),
            "all_minimum_sets": minimum_sets,
            "full_bad_place_indistinguishability_classes": sorted(
                full_signature_classes.values()
            ),
        },
        "canonical_direction_graph_components_are_not_intrinsic_module_summands": True,
    }


def target_comparison(curves: dict[int, dict[str, Any]]):
    local = {
        curve_id: {
            int(row["rational_prime"]): row for row in curves[curve_id]["local_places"]
        }
        for curve_id in RECORD_IDS
    }
    two_signatures = {
        str(curve_id): local[curve_id][2]["selected_block_canonical_signature"]
        for curve_id in RECORD_IDS
    }
    common_good = []
    for prime in COMMON_GOOD_PRIMES:
        left = local[356][prime]
        right = local[385][prime]
        left_signature = left["selected_block_canonical_signature"]
        right_signature = right["selected_block_canonical_signature"]
        common_good.append(
            {
                "prime": prime,
                "curve356_ambient_local_kummer_dimension": left[
                    "ambient_local_kummer_dimension"
                ],
                "curve385_ambient_local_kummer_dimension": right[
                    "ambient_local_kummer_dimension"
                ],
                "ambient_dimensions_equal": left["ambient_local_kummer_dimension"]
                == right["ambient_local_kummer_dimension"],
                "full_selected_block_signature_equal": left_signature
                == right_signature,
                "curve356_signature_sha256": left_signature["sha256"],
                "curve385_signature_sha256": right_signature["sha256"],
            }
        )

    def bad_category_multiset(curve_id: int):
        values = []
        for row in curves[curve_id]["local_places"]:
            if row["place_kind"] != "bad":
                continue
            signature = dict(row["selected_block_canonical_signature"])
            signature.pop("sha256")
            values.append(signature)
        return counted_multiset(values)

    left_bad = bad_category_multiset(356)
    right_bad = bad_category_multiset(385)
    return {
        "record_pair": [356, 385],
        "selected_blocks": {
            "356": curves[356]["selected_comparison_block"],
            "385": curves[385]["selected_comparison_block"],
        },
        "two_adic": {
            "curve356_signature": two_signatures["356"],
            "curve385_signature": two_signatures["385"],
            "signatures_equal": two_signatures["356"] == two_signatures["385"],
            "curve356_component_image_order_multiset": counted_multiset(
                [
                    row["component_group_image_order"]
                    for row in local[356][2]["directions"]
                    if row["label"] in curves[356]["selected_comparison_block"]
                ]
            ),
            "curve385_component_image_order_multiset": counted_multiset(
                [
                    row["component_group_image_order"]
                    for row in local[385][2]["directions"]
                    if row["label"] in curves[385]["selected_comparison_block"]
                ]
            ),
        },
        "fixed_common_good_places": common_good,
        "common_good_ambient_dimension_agreement_count": sum(
            row["ambient_dimensions_equal"] for row in common_good
        ),
        "common_good_full_block_signature_agreement_count": sum(
            row["full_selected_block_signature_equal"] for row in common_good
        ),
        "bad_place_anonymous_category_multisets": {
            "356": left_bad,
            "385": right_bad,
            "equal": left_bad == right_bad,
        },
    }


def control_comparison(curves: dict[int, dict[str, Any]]):
    rows = []
    ids = sorted(curves)
    for left_index, left in enumerate(ids):
        for right in ids[left_index + 1 :]:
            left_vector = curves[left]["fixed_common_good_ambient_dimension_vector"]
            right_vector = curves[right]["fixed_common_good_ambient_dimension_vector"]
            rows.append(
                {
                    "left_curve": left,
                    "right_curve": right,
                    "agreement_count_out_of_8": sum(
                        a == b for a, b in zip(left_vector, right_vector)
                    ),
                    "vectors_equal": left_vector == right_vector,
                }
            )
    record_agreement = next(
        row["agreement_count_out_of_8"]
        for row in rows
        if (row["left_curve"], row["right_curve"]) == (356, 385)
    )
    maximum_other_agreement = max(
        row["agreement_count_out_of_8"]
        for row in rows
        if (row["left_curve"], row["right_curve"]) != (356, 385)
    )
    return {
        "fixed_common_good_dimension_vectors": {
            str(curve_id): curves[curve_id]["fixed_common_good_ambient_dimension_vector"]
            for curve_id in ids
        },
        "pairwise_agreements": rows,
        "record_pair_agreement_count_out_of_8": record_agreement,
        "maximum_other_pair_agreement_count_out_of_8": maximum_other_agreement,
        "record_pair_is_uniquely_more_similar_than_controls": (
            record_agreement > maximum_other_agreement
        ),
    }


def build():
    public, lineage, curve12, rigid = load_inputs()
    bases = quotient_bases(lineage, curve12)
    complements, complement_provenance, quotient_presentations = rigid_complements(
        rigid, bases
    )
    public_by_id = {int(row["id"]): row for row in public["records"]}
    parameters = {
        int(row["curve_id"]): row["parameter"]
        for row in lineage["target_isomorphisms"]
        if row["chart"] == "norm12-orbit-074d9"
    }
    parameters[12] = curve12["native_fibre"]["parameter"]
    charts = {curve_id: "norm12-orbit-074d9" for curve_id in TARGET_IDS}
    charts[12] = "norm12-orbit-11952"

    curves = {}
    for curve_id in TARGET_IDS:
        selected_labels = complements.get(curve_id, bases[curve_id])
        curves[curve_id] = audit_curve(
            public_by_id[curve_id],
            bases[curve_id],
            selected_labels,
            parameters[curve_id],
            charts[curve_id],
            quotient_presentations.get(curve_id),
        )

    arithmetic_blocks = {
        str(curve_id): quotient_arithmetic_block_decomposition(curves[curve_id])
        for curve_id in RECORD_IDS
    }
    expected_obstructions = {
        "356": ({"13", "23", "37", "139"}, 4, 5, 22),
        "385": ({"5", "29", "37", "41", "73", "109", "127"}, 10, 4, 1),
    }
    for curve_id, (primes, span_dimension, minimum_size, minimum_count) in expected_obstructions.items():
        block = arithmetic_blocks[curve_id]
        actual_primes = {
            row["rational_prime"] for row in block["pairing_descent_obstructions"]
        }
        minimum = block[
            "minimum_bad_places_distinguishing_canonical_directions_without_hilbert_data"
        ]
        if block["status"] != "NOT_DEFINED_COMPONENT_HILBERT_FORMS_FAIL_TO_DESCEND":
            raise ArithmeticError("the expected quotient Hilbert descent obstruction disappeared")
        if actual_primes != primes:
            raise ArithmeticError("the quotient Hilbert obstruction-prime set changed")
        if block["pairing_descent_obstruction_span"]["dimension"] != span_dimension:
            raise ArithmeticError("the quotient Hilbert obstruction-span dimension changed")
        if (
            minimum["minimum_cardinality"] != minimum_size
            or minimum["number_of_minimum_sets"] != minimum_count
        ):
            raise ArithmeticError("the quotient-compatible distinguishing-place minimum changed")

    record = target_comparison(curves)
    controls = control_comparison(curves)
    if record["two_adic"]["signatures_equal"]:
        raise ArithmeticError("the expected exact two-adic separation disappeared")
    if record["common_good_ambient_dimension_agreement_count"] != 3:
        raise ArithmeticError("the fixed common-good dimension comparison changed")
    if record["common_good_full_block_signature_agreement_count"] != 0:
        raise ArithmeticError("a fixed common-good block signature unexpectedly matches")
    if record["bad_place_anonymous_category_multisets"]["equal"]:
        raise ArithmeticError("the anonymous bad-place profiles unexpectedly match")
    if controls["record_pair_is_uniquely_more_similar_than_controls"]:
        raise ArithmeticError("the record pair became uniquely similar among controls")

    output = {
        "schema": SCHEMA,
        "status": "PASS_EXACT_QUOTIENT_DATA_AND_HILBERT_DESCENT_OBSTRUCTION",
        "summary": {
            "record_curves": list(RECORD_IDS),
            "control_curves": [351, 376, 377, 12],
            "record_exceptional_quotient_dimensions": [12, 12],
            "record_rigid_span_dimensions": [2, 2],
            "record_selected_complement_dimensions": [10, 10],
            "fixed_common_good_primes": list(COMMON_GOOD_PRIMES),
            "all_bad_primes_covered": True,
            "all_global_cubic_discriminant_factorizations_proved_complete": True,
            "all_pairwise_local_tate_products_trivial": True,
            "record_two_adic_signatures_equal": False,
            "record_common_good_full_block_signature_agreement_count": 0,
            "record_common_good_ambient_dimension_agreement_count": 3,
            "record_quotient_arithmetic_block_statuses": {
                curve_id: arithmetic_blocks[curve_id]["status"]
                for curve_id in sorted(arithmetic_blocks)
            },
            "component_hilbert_descent_obstruction_primes": {
                curve_id: sorted(
                    {
                        row["rational_prime"]
                        for row in arithmetic_blocks[curve_id][
                            "pairing_descent_obstructions"
                        ]
                    },
                    key=ZZ,
                )
                for curve_id in sorted(arithmetic_blocks)
            },
            "conclusion": (
                "The ten-dimensional rigid-invisible blocks at curves 356 and 385 "
                "do not carry the proposed quotient Hilbert-pairing structure: "
                "componentwise Hilbert forms fail to annihilate the visible relation "
                "space at split bad primes, so no intrinsic 10=... partition is defined."
            ),
        },
        "method": {
            "ambient_local_dimension": (
                "If r_p is the number of field factors of the cubic etale "
                "algebra over Q_p, then dim E(Q_p)/2E(Q_p)=r_p-1 for odd p "
                "and r_p for p=2."
            ),
            "odd_local_squareclasses": (
                "Each local field factor uses valuation parity and the "
                "squareclass of the residue of its unit part."
            ),
            "two_adic_block_rank": (
                "Exact PARI nfislocalpower tests greedily row-reduce the "
                "displayed Kummer elements in the product of completions."
            ),
            "component_image": (
                "The automorphism-invariant image order is the least divisor "
                "m of the Tamagawa number for which mP lies in E_0(Q_p)."
            ),
            "pairing": (
                "PARI nfhilbert computes each component symbol; their product "
                "is the corestricted local Tate pairing in the standard "
                "cubic 2-descent realization."
            ),
        },
        "rigid_invisible_complement_provenance": {
            str(curve_id): row for curve_id, row in complement_provenance.items()
        },
        "record_quotient_arithmetic_blocks": arithmetic_blocks,
        "curves": [curves[curve_id] for curve_id in TARGET_IDS],
        "record_pair_comparison": record,
        "control_comparison": controls,
        "crt_and_inward_search_gate": {
            "status": "NOT_RUN_NEGATIVE_LOCAL_FINGERPRINT_GATE",
            "reason": (
                "No invariant signature is shared by the two ten-direction blocks; "
                "there are therefore no justified common local congruences to combine."
            ),
        },
        "claim_boundary": [
            "The exact Kummer classes concern the certified displayed subgroups and do not assert that those subgroups are the full Mordell-Weil groups.",
            "The local Tate products are universally trivial on local point images; discrimination comes from anonymous componentwise Hilbert data, valuation supports, local dimensions, and component-image orders.",
            "The negative conclusion rejects a shared mechanism characterized by this declared local fingerprint. It does not prove that no more abstract common rank-jump mechanism can exist.",
            "The common-good block is fixed in advance. No post-hoc prime selection or CRT parameter search is performed.",
            "The earlier RREF nonpivot blocks are deterministic coordinate complements, not quotient-invariant Hilbert modules. The new descent ledger reports no indecomposable partition when a visible relation is outside a component-form radical.",
        ],
        "inputs": {
            relative(path): digest(path)
            for path in (PUBLIC, LINEAGE, CURVE12, RIGID)
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_074d9_local_kummer_meet.sage --check"
        ),
        "software_assumptions": {
            "sage": str(sage_version),
            "pari": ".".join(str(part) for part in pari.version()),
        },
    }
    return output


def legacy_projection(output: dict[str, Any]) -> dict[str, Any]:
    """Retain byte-replay of the superseded coordinate-complement certificate."""

    legacy = json.loads(json.dumps(output))
    legacy["schema"] = LEGACY_SCHEMA
    legacy["status"] = "PASS_EXACT_NEGATIVE_LOCAL_KUMMER_MEET"
    legacy["summary"].pop("record_quotient_arithmetic_block_statuses")
    legacy["summary"].pop("component_hilbert_descent_obstruction_primes")
    legacy["summary"]["conclusion"] = (
        "The ten-dimensional rigid-invisible blocks at curves 356 and 385 "
        "do not share the declared invariant local Kummer fingerprint."
    )
    legacy.pop("record_quotient_arithmetic_blocks")
    for provenance in legacy["rigid_invisible_complement_provenance"].values():
        provenance.pop("rigid_span_exact_integral_vectors")
    for curve in legacy["curves"]:
        curve.pop("visible_quotient_presentation", None)
        for place in curve["local_places"]:
            place.pop("visible_quotient_arithmetic", None)
    legacy["claim_boundary"] = legacy["claim_boundary"][:-1]
    return legacy


def quotient_certificate_projection(output: dict[str, Any]) -> dict[str, Any]:
    """Keep the new certificate narrow while the full v1 replay stays preserved."""

    return {
        key: output[key]
        for key in (
            "schema",
            "status",
            "summary",
            "method",
            "rigid_invisible_complement_provenance",
            "record_quotient_arithmetic_blocks",
            "claim_boundary",
            "inputs",
            "reproducing_command",
            "software_assumptions",
        )
    } | {
        "preserved_coordinate_complement_certificate": {
            "path": relative(LEGACY_OUTPUT),
            "sha256": digest(LEGACY_OUTPUT),
            "replaying_command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_074d9_local_kummer_meet.sage "
                "--legacy-coordinate-complement --check"
            ),
        }
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--legacy-coordinate-complement",
        action="store_true",
        help="replay the preserved v1 coordinate-complement artifact",
    )
    args = parser.parse_args()
    record = build()
    if args.legacy_coordinate_complement:
        record = legacy_projection(record)
    else:
        record = quotient_certificate_projection(record)
    serialized = json.dumps(record, indent=2, sort_keys=True) + "\n"
    default_output = LEGACY_OUTPUT if args.legacy_coordinate_complement else OUTPUT
    output = (args.output or default_output).resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored local-Kummer meet differs from exact replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        f"{PROTOCOL}|records=356,385|controls=351,376,377,12|"
        f"common_good=8|result={'LEGACY_COORDINATE_COMPLEMENT' if args.legacy_coordinate_complement else 'QUOTIENT_HILBERT_DESCENT_OBSTRUCTION'}|status=PASS|"
        f"output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
