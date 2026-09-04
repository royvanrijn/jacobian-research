#!/usr/bin/env sage-python
"""Certify the NS0031 rational-marking modular obstruction.

The exact ternary transcendental lattice has split even Clifford algebra.  This
checker identifies its integral norm-one group as

    Gamma_ns(4) intersection Gamma_0(37),

computes the congruence signature, and excludes the two noncuspidal rational
points of X_0(37) by a mod-4 Frobenius trace test.  The cited theorem of Velu
that determines X_0(37)(QQ), and the standard marked-K3/ternary-spin moduli
correspondence, are theorem inputs recorded in the canonical proof note.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from math import gcd
from pathlib import Path

from sage.all import CliffordAlgebra, EllipticCurve_from_j, QQ, ZZ, QuadraticForm, matrix


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CATALOGUE = GENERATED / "elkies-k3-rank7-auxiliary-catalogue-v1.json"
T_ARITHMETIC = GENERATED / "elkies-k3-rank7-t-arithmetic-v1.json"
SOURCE = GENERATED / (
    "elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
OUTPUT = GENERATED / "elkies-k3-ns0031-qq-marking-obstruction-v1.json"
SURFACE_ID = "K3-d1b1381f87d69f1c"
SOURCE_ID = "NS0031-S001"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def tuple_matrix(value, modulus):
    return tuple(int(entry % modulus) for entry in value.list())


def multiply(left, right, modulus):
    a, b, c, d = left
    e, f, g, h = right
    return (
        (a * e + b * g) % modulus,
        (a * f + b * h) % modulus,
        (c * e + d * g) % modulus,
        (c * f + d * h) % modulus,
    )


def determinant(value, modulus):
    a, b, c, d = value
    return (a * d - b * c) % modulus


def trace(value, modulus):
    return (value[0] + value[3]) % modulus


def sl2(modulus):
    return [
        value
        for value in product(range(modulus), repeat=4)
        if determinant(value, modulus) == 1 % modulus
    ]


def left_cosets(group, subgroup, modulus):
    unseen = set(group)
    representatives = []
    lookup = {}
    while unseen:
        representative = min(unseen)
        coset = {multiply(element, representative, modulus) for element in subgroup}
        index = len(representatives)
        representatives.append(representative)
        for element in coset:
            lookup[element] = index
        unseen -= coset
    return representatives, lookup


def projective_line_prime(prime):
    return [(0, 1)] + [(1, value) for value in range(prime)]


def canonical_projective_row(value, prime):
    left, right = value
    if left % prime:
        return (1, (right * pow(left, -1, prime)) % prime)
    return (0, 1)


def cycle_lengths(permutation):
    seen = set()
    result = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        result.append(length)
    return sorted(result)


def build_payload():
    catalogue = json.loads(CATALOGUE.read_text())
    t_payload = json.loads(T_ARITHMETIC.read_text())
    source_payload = json.loads(SOURCE.read_text())

    surface = next(
        row for row in catalogue["surfaces"] if row["surface_id"] == SURFACE_ID
    )
    t_row = next(
        row for row in t_payload["surfaces"] if row["surface_id"] == SURFACE_ID
    )
    source = next(
        row for row in source_payload["sources"] if row["source_id"] == SOURCE_ID
    )

    transcendental = matrix(ZZ, surface["surface_key"]["transcendental_gram"])
    expected = matrix(ZZ, [[0, 0, 4], [0, 74, 1], [4, 1, -2]])
    assert transcendental == expected
    assert transcendental.det() == -1184
    assert t_row["literal_transcendental_gram"] == rows(transcendental)
    assert t_row["rational_isotropy"]["primitive_isotropic_vector"] == [1, 0, 0]
    assert t_row["rational_isotropy"]["integral_u_split"]["isotropic_divisibility"] == 4
    assert t_row["clifford"]["quaternion_discriminant"] == 1
    assert t_row["clifford"]["integral_even_clifford_order"]["reduced_discriminant"] == 592

    # Exact split representation of the even Clifford order.  The basis is
    # 1, e0e1, e0e2, e1e2.  Conjugation by diag(4,1) makes it integral.
    quadratic = QuadraticForm(QQ, transcendental)
    clifford = CliffordAlgebra(quadratic)
    e0, e1, e2 = clifford.gens()
    even_basis = [clifford.one(), e0 * e1, e0 * e2, e1 * e2]
    split_basis = [
        matrix(QQ, [[1, 0], [0, 1]]),
        matrix(QQ, [[0, 1], [0, 0]]),
        matrix(QQ, [[4, 0], [0, 0]]),
        matrix(QQ, [[1, QQ(1) / 4], [148, 0]]),
    ]
    clifford_keys = list(clifford.basis().keys())

    def coordinates(value):
        coefficients = value.monomial_coefficients()
        return [coefficients.get(clifford_keys[index], 0) for index in (0, 4, 5, 6)]

    def split_image(value):
        coordinate = coordinates(value)
        return sum(
            (coordinate[index] * split_basis[index] for index in range(4)),
            matrix(QQ, 2, 2),
        )

    for left in even_basis:
        for right in even_basis:
            assert split_image(left * right) == split_image(left) * split_image(right)

    conjugator = matrix(QQ, [[4, 0], [0, 1]])
    integral_basis = [
        conjugator * value * conjugator.inverse() for value in split_basis
    ]
    expected_integral_basis = [
        matrix(ZZ, [[1, 0], [0, 1]]),
        matrix(ZZ, [[0, 4], [0, 0]]),
        matrix(ZZ, [[4, 0], [0, 0]]),
        matrix(ZZ, [[1, 1], [37, 0]]),
    ]
    assert integral_basis == expected_integral_basis

    # O = {[[A,B],[C,D]] in M_2(ZZ): 37|C,
    #       B=C/37 (mod 4), A-D=C/37 (mod 4)}.
    # Hence O^1 reduces to the nonsplit Cartan norm-one group at 4 and the
    # upper triangular Borel at 37.
    group4 = sl2(4)
    h4 = sorted(
        value
        for value in group4
        if (value[2] - value[1]) % 4 == 0
        and (value[2] - value[0] + value[3]) % 4 == 0
    )
    cartan4 = sorted(
        (
            (a + b) % 4,
            b % 4,
            b % 4,
            a % 4,
        )
        for a, b in product(range(4), repeat=2)
        if ((a + b) * a - b * b) % 2 == 1
    )
    cartan4 = sorted(set(cartan4))
    assert len(cartan4) == 12
    assert h4 == sorted(value for value in cartan4 if determinant(value, 4) == 1)
    assert len(h4) == 6 and len(group4) == 48
    # theta^2-theta-1 is irreducible modulo 2, so this is the unramified
    # nonsplit Cartan, not a split Cartan in disguised coordinates.
    assert all((root * root - root - 1) % 2 for root in range(2))

    representatives4, lookup4 = left_cosets(group4, h4, 4)
    assert len(representatives4) == 8
    prime = 37
    projective = projective_line_prime(prime)
    points = [(index, value) for index in range(8) for value in projective]
    point_lookup = {value: index for index, value in enumerate(points)}

    def permutation(operator):
        operator4 = tuple(value % 4 for value in operator)
        operator37 = tuple(value % prime for value in operator)
        result = []
        for index4, projective_row in points:
            image4 = lookup4[multiply(representatives4[index4], operator4, 4)]
            a, b, c, d = operator37
            x, y = projective_row
            image37 = canonical_projective_row(
                ((x * a + y * c) % prime, (x * b + y * d) % prime), prime
            )
            result.append(point_lookup[(image4, image37)])
        return result

    s_operator = (0, -1, 1, 0)
    r_operator = (0, -1, 1, 1)
    t_operator = (1, 1, 0, 1)
    s_cycles = cycle_lengths(permutation(s_operator))
    r_cycles = cycle_lengths(permutation(r_operator))
    t_cycles = cycle_lengths(permutation(t_operator))
    modular_index = len(points)
    elliptic2 = s_cycles.count(1)
    elliptic3 = r_cycles.count(1)
    cusps = len(t_cycles)
    genus = QQ(1) + QQ(modular_index) / 12 - QQ(elliptic2) / 4 - QQ(elliptic3) / 3 - QQ(cusps) / 2
    assert modular_index == 304
    assert elliptic2 == 0 and elliptic3 == 4
    assert t_cycles == [4, 4, 148, 148]
    assert genus == 23

    # Velu's exact determination of X_0(37)(QQ) leaves two noncuspidal
    # rational points.  Their j-values both fail the necessary Cartan trace
    # condition at the good prime 19.  Quadratic twisting changes a_19 by a
    # sign, which is still 2 modulo 4.
    rational_j_values = [-7 * 11**3, -7 * 137**3 * 2083**3]
    cartan_trace_determinants = sorted(
        {(trace(value, 4), determinant(value, 4)) for value in cartan4}
    )
    assert (2, 3) not in cartan_trace_determinants
    exclusions = []
    for j_value in rational_j_values:
        curve = EllipticCurve_from_j(QQ(j_value)).global_minimal_model()
        assert curve.j_invariant() == j_value
        assert gcd(int(curve.discriminant()), 19) == 1
        isogenies = curve.isogenies_prime_degree(ZZ(37))
        assert [int(isogeny.degree()) for isogeny in isogenies] == [37]
        frobenius_trace = int(curve.ap(19))
        assert frobenius_trace == -6
        trace_determinant = (frobenius_trace % 4, 19 % 4)
        assert trace_determinant == (2, 3)
        assert trace_determinant not in cartan_trace_determinants
        exclusions.append(
            {
                "j": j_value,
                "minimal_ainvariants": list(map(int, curve.ainvs())),
                "rational_isogeny_degree": 37,
                "good_prime": 19,
                "frobenius_trace": frobenius_trace,
                "trace_determinant_mod_4": list(trace_determinant),
                "in_nonsplit_cartan_trace_determinant_set": False,
                "quadratic_twist_invariance": "-a19 is also 2 modulo 4",
            }
        )

    assert source["determinant"] == 1184
    assert source["source"]["root_type"] == "A1+2A7"
    assert source["source"]["root_lattice_primitive"] is True
    assert source["source"]["torsion"] == 1
    assert source["source"]["mw_height_gram"] == [["2", "1"], ["1", "41/8"]]

    return {
        "schema": "elkies-k3.ns0031-qq-marking-obstruction.v1",
        "status": "PASS_NS0031_QQ_RATIONAL_MARKING_OBSTRUCTION",
        "surface_id": SURFACE_ID,
        "ns_id": "NS0031",
        "ns_rank": 19,
        "ns_determinant": 1184,
        "transcendental_lattice": {
            "gram": rows(transcendental),
            "determinant": int(transcendental.det()),
            "primitive_isotropic_vector": [1, 0, 0],
            "isotropic_divisibility": 4,
        },
        "even_clifford_order": {
            "split_basis_before_conjugation": [
                [[str(entry) for entry in row] for row in value.rows()]
                for value in split_basis
            ],
            "integral_basis_after_diag_4_1_conjugation": [
                rows(value) for value in integral_basis
            ],
            "integral_matrix_conditions": [
                "37 divides C",
                "B = C/37 modulo 4",
                "A-D = C/37 modulo 4",
            ],
            "reduced_discriminant": 592,
        },
        "norm_one_modular_curve": {
            "label": "X_ns(4) fiber_product_X(1) X_0(37)",
            "group": "Gamma_ns(4) intersection Gamma_0(37)",
            "congruence_level": 148,
            "index_in_PSL2Z": modular_index,
            "elliptic_orbits_order_2": elliptic2,
            "elliptic_orbits_order_3": elliptic3,
            "cusp_widths": t_cycles,
            "genus": int(genus),
            "mod_4_full_nonsplit_cartan_order": len(cartan4),
            "mod_4_norm_one_order": len(h4),
        },
        "x0_37_rational_point_gate": {
            "external_theorem": "Velu's determination of X_0(37)(QQ)",
            "noncuspidal_rational_j_values": rational_j_values,
            "nonsplit_cartan_trace_determinants_mod_4": [
                list(value) for value in cartan_trace_determinants
            ],
            "excluded_lifts": exclusions,
            "conclusion": (
                "The fibre product X_ns(4) x_{X(1)} X_0(37) has no "
                "noncuspidal QQ-rational point."
            ),
        },
        "source_frame": {
            "source_id": SOURCE_ID,
            "root_type": "A1+2A7",
            "mw_rank": 2,
            "mw_height_gram": [["2", "1"], ["1", "41/8"]],
            "root_lattice_primitive": True,
            "torsion": 1,
        },
        "arithmetic_conclusion": (
            "No characteristic-zero K3 over QQ with geometric NS=NS0031 can "
            "have all nineteen Neron-Severi divisor classes rational. In "
            "particular a rootless NS0031 fibration cannot have a saturated "
            "rational MW17 basis over QQ(t)."
        ),
        "external_references": [
            {
                "name": "J. Velu, Les points rationnels de X_0(37)",
                "citation": "Bull. Soc. Math. France, Memoire 37 (1974), 169-179",
                "doi": "10.24033/msmf.145",
                "url": "https://www.numdam.org/articles/10.24033/msmf.145/",
            }
        ],
        "input_hashes": {
            relative(CATALOGUE): digest(CATALOGUE),
            relative(T_ARITHMETIC): digest(T_ARITHMETIC),
            relative(SOURCE): digest(SOURCE),
        },
        "proof_boundary": {
            "proved": (
                "The exact split Clifford-order model, congruence subgroup, "
                "signature, and mod-4 Frobenius exclusions of both "
                "noncuspidal rational X_0(37) points are replayed here. With "
                "Velu's theorem and the standard marked-K3 ternary-spin "
                "moduli correspondence, a full QQ-rational NS0031 marking is "
                "impossible."
            ),
            "not_proved": (
                "The checker does not reprove Velu's global determination of "
                "X_0(37)(QQ), or the general period/Clifford moduli theorem. "
                "It does not obstruct geometric NS0031 surfaces, models over "
                "larger number fields, or QQ models with only a proper "
                "Galois-invariant sublattice."
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    rendered = json.dumps(build_payload(), indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != rendered:
            raise SystemExit(f"stale artifact: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    print("NS0031QQMARKING|group=Gamma_ns(4)&Gamma0(37)|genus=23|status=PASS")


if __name__ == "__main__":
    main()
