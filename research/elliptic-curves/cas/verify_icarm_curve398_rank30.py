#!/usr/bin/env python3
"""Certify ICARM curve 398 and record its public construction boundary.

The rank lower bound uses exact point membership and exhaustive finite
good-reduction quotients.  PARI/GP is used independently for local reduction,
primality, the root number, and a numerical height diagnostic.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve398 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    POINTS,
    SHORT_POINTS,
    c4,
    c6,
    discriminant,
    on_curve,
    short_coefficients,
)
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    finite_curve_points,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


PROTOCOL = "R30ICARM398"
ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_curve398_rank30_and_construction_v1.json"
)
A1_REDUCTION = ROOT / "elliptic-curves/data/icarm_curve398_known_a1_mod179.json"
REPRODUCING_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/verify_icarm_curve398_rank30.py"
)

PUBLIC_CURVE_SOURCE = "https://elliptic-rank.icarm.cloud/curve/398"
PUBLIC_JSON_SOURCE = "https://elliptic-rank.icarm.cloud/curve/398.json"
PUBLIC_JSON_SHA256_2026_09_04 = (
    "0682131699cdfe07283aeb09b19a448bb375f8f1765a01eb5b14c0aa7105994c"
)
PRIMARY_ANNOUNCEMENT = (
    "https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;d593eaaa.2608&S="
)

PUBLIC_CONDUCTOR = int(
    "2835647668537242470520670933169702570972603599817629325681756626410520535791722336193301955569324372206950385073317702302276416849610"
)
PUBLIC_DISCRIMINANT = int(
    "1311021171349839305076869751181166178696551783990116503463598327754948038424504439658681699321482006153838737418363257993512584919966486940354703657840019569754443427840000000"
)

# prime -> discriminant valuation.  All reductions are multiplicative I_n.
DISCRIMINANT_FACTORIZATION = (
    (2, 18),
    (3, 10),
    (5, 7),
    (7, 4),
    (19, 3),
    (29, 4),
    (31, 2),
    (37, 2),
    (41, 2),
    (59, 2),
    (73, 2),
    (79, 3),
    (101, 2),
    (151, 2),
    (197, 2),
    (20084070565614383, 1),
    (277131689105980733414153, 1),
    (
        91586120381369539248864736998169886452607476058975046580346704317125237,
        1,
    ),
)

# prime -> (conductor exponent, PARI Kodaira code, Tamagawa number, local sign)
EXPECTED_LOCAL = {
    2: (1, 22, 18, -1),
    3: (1, 14, 10, -1),
    5: (1, 11, 7, -1),
    7: (1, 8, 4, -1),
    19: (1, 7, 3, -1),
    29: (1, 8, 4, -1),
    31: (1, 6, 2, 1),
    37: (1, 6, 2, 1),
    41: (1, 6, 2, 1),
    59: (1, 6, 2, 1),
    73: (1, 6, 2, 1),
    79: (1, 7, 3, -1),
    101: (1, 6, 2, -1),
    151: (1, 6, 2, 1),
    197: (1, 6, 2, -1),
    20084070565614383: (1, 5, 1, 1),
    277131689105980733414153: (1, 5, 1, -1),
    91586120381369539248864736998169886452607476058975046580346704317125237: (
        1,
        5,
        1,
        -1,
    ),
}

EXPECTED_CERTIFICATE_PRIMES = (
    11,
    13,
    17,
    43,
    47,
    53,
    61,
    67,
    83,
    103,
    107,
    109,
    113,
    131,
    139,
    149,
    157,
    173,
    179,
    193,
    199,
    211,
    227,
    229,
    239,
)

# Counts of displayed points reducing to the node of the singular cubic.
# A missing count means the point denominator is divisible by that prime.
EXPECTED_NODE_INCIDENCE = {
    2: (29, (3,)),
    3: (28, (15,)),
    5: (29, ()),
    7: (24, ()),
    19: (25, ()),
    29: (28, ()),
    31: (22, ()),
    37: (24, ()),
    41: (21, ()),
    59: (23, ()),
    73: (24, ()),
    79: (26, ()),
    101: (15, ()),
    151: (24, ()),
    197: (21, ()),
    20084070565614383: (0, ()),
    277131689105980733414153: (0, ()),
    91586120381369539248864736998169886452607476058975046580346704317125237: (
        0,
        (),
    ),
}

LOCAL_PATTERN = re.compile(
    r"^LOCAL\|(\d+)\|(\d+)\|(-?\d+)\|\[([^]]+)\]\|(\d+)\|(-?\d+)$",
    re.MULTILINE,
)


def factor_product(factors: tuple[tuple[int, int], ...]) -> int:
    answer = 1
    for prime, exponent in factors:
        answer *= prime**exponent
    return answer


def gp_rational(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def finite_group_order(prime: int) -> int:
    coefficients = short_coefficients()
    return len(
        finite_curve_points(
            int(coefficients[3]) % prime,
            int(coefficients[4]) % prime,
            prime,
        )
    )


def _trim(polynomial: list[int], prime: int) -> list[int]:
    polynomial = [coefficient % prime for coefficient in polynomial]
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def _remainder(left: list[int], right: list[int], prime: int) -> list[int]:
    left = _trim(left, prime)
    right = _trim(right, prime)
    inverse = pow(right[-1], -1, prime)
    while len(left) >= len(right) and left != [0]:
        multiplier = left[-1] * inverse % prime
        offset = len(left) - len(right)
        for index, coefficient in enumerate(right):
            left[index + offset] -= multiplier * coefficient
        left = _trim(left, prime)
    return left


def _polynomial_gcd(
    left: list[int], right: list[int], prime: int
) -> list[int]:
    left = _trim(left, prime)
    right = _trim(right, prime)
    while right != [0]:
        left, right = right, _remainder(left, right, prime)
    inverse = pow(left[-1], -1, prime)
    return _trim([coefficient * inverse for coefficient in left], prime)


def singular_node(prime: int) -> tuple[int, int]:
    a4 = int(GENERAL_WEIERSTRASS_COEFFICIENTS[3])
    a6 = int(GENERAL_WEIERSTRASS_COEFFICIENTS[4])
    if prime <= 3:
        nodes = []
        for x_value in range(prime):
            for y_value in range(prime):
                equation = y_value**2 + x_value * y_value - x_value**3 - a4 * x_value - a6
                derivative_x = y_value - 3 * x_value**2 - a4
                derivative_y = 2 * y_value + x_value
                if equation % prime == derivative_x % prime == derivative_y % prime == 0:
                    nodes.append((x_value, y_value))
        if len(nodes) != 1:
            raise AssertionError(f"expected one node modulo {prime}, got {nodes}")
        return nodes[0]

    # Completing the square gives z^2=g(x), with z=2*y+x.  The node has
    # z=0 and x equal to the common root of g and g'.
    common = _polynomial_gcd(
        [4 * a6, 4 * a4, 1, 4],
        [4 * a4, 2, 12],
        prime,
    )
    if len(common) != 2:
        raise AssertionError(f"singular gcd modulo {prime} is not linear")
    x_value = -common[0] % prime
    y_value = -x_value * pow(2, -1, prime) % prime
    return x_value, y_value


def node_incidence() -> dict[int, tuple[int, tuple[int, ...]]]:
    answer = {}
    for prime, _valuation in DISCRIMINANT_FACTORIZATION:
        node = singular_node(prime)
        hits = 0
        unavailable = []
        for index, (x_value, y_value) in enumerate(POINTS, 1):
            if x_value.denominator % prime == 0 or y_value.denominator % prime == 0:
                unavailable.append(index)
                continue
            reduced = (
                x_value.numerator * pow(x_value.denominator, -1, prime) % prime,
                y_value.numerator * pow(y_value.denominator, -1, prime) % prime,
            )
            hits += reduced == node
        answer[prime] = (hits, tuple(unavailable))
    return answer


def _convolution(left: list[int], right: list[int], prime: int) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return [value % prime for value in answer]


def verify_known_a1_exclusion() -> dict[str, object]:
    document = json.loads(A1_REDUCTION.read_text(encoding="utf-8"))
    prime = int(document["prime"])
    family = document["family"]
    coefficient_a = list(family["A_coefficients_low_to_high_mod_prime"])
    coefficient_b = list(family["B_coefficients_low_to_high_mod_prime"])
    a_cubed = _convolution(_convolution(coefficient_a, coefficient_a, prime), coefficient_a, prime)
    b_squared = _convolution(coefficient_b, coefficient_b, prime)
    family_discriminant = [
        -16
        * (
            4 * (a_cubed[index] if index < len(a_cubed) else 0)
            + 27 * (b_squared[index] if index < len(b_squared) else 0)
        )
        % prime
        for index in range(max(len(a_cubed), len(b_squared)))
    ]
    family_c4 = [(-48 * value) % prime for value in coefficient_a]
    family_c4_cubed = _convolution(
        _convolution(family_c4, family_c4, prime), family_c4, prime
    )
    degree = max(len(family_c4_cubed), len(family_discriminant))
    recognition = [
        (
            (family_c4_cubed[index] if index < len(family_c4_cubed) else 0)
            * (PUBLIC_DISCRIMINANT % prime)
            - pow(c4(), 3, prime)
            * (
                family_discriminant[index]
                if index < len(family_discriminant)
                else 0
            )
        )
        % prime
        for index in range(degree)
    ]
    recognition = _trim(recognition, prime)
    expected = document["recognition"]
    if recognition != expected["coefficients_low_to_high_mod_prime"]:
        raise AssertionError("known A1 recognition polynomial changed")
    roots = [
        value
        for value in range(prime)
        if sum(
            coefficient * pow(value, index, prime)
            for index, coefficient in enumerate(recognition)
        )
        % prime
        == 0
    ]
    if len(recognition) - 1 != 24 or recognition[-1] == 0 or roots:
        raise AssertionError("known A1 no-projective-root witness failed")
    return {
        "prime": prime,
        "degree": len(recognition) - 1,
        "finite_roots": roots,
        "root_at_infinity": False,
        "coefficients_sha256": sha256(
            ",".join(map(str, recognition)).encode()
        ).hexdigest(),
        "conclusion": expected["conclusion"],
    }


def pari_diagnostics() -> dict[str, object]:
    ainvs = ",".join(gp_rational(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS)
    points = ",".join(
        f"[{gp_rational(x_value)},{gp_rational(y_value)}]"
        for x_value, y_value in POINTS
    )
    commands = [
        "default(realprecision,100);",
        f"E=ellinit([{ainvs}]);",
    ]
    for prime in EXPECTED_LOCAL:
        commands.extend(
            [
                f"L=elllocalred(E,{prime});",
                (
                    f'print("LOCAL|{prime}|",L[1],"|",L[2],"|",L[3],"|",'
                    f'L[4],"|",ellrootno(E,{prime}));'
                ),
                f'print("PRIME|{prime}|",isprime({prime}));',
            ]
        )
    commands.extend(
        [
            "ISOMAT=ellisomat(E);",
            'print("ISOGENY_CLASS_SIZE|",#ISOMAT[1]);',
            'print("ISOGENY_DEGREE_MATRIX|",ISOMAT[2]);',
            f"P=[{points}];",
            "H=ellheightmatrix(E,P);",
            'print("REGULATOR|",matdet(H));',
            'print("HEIGHT_FIRST|",H[1,1]);',
            'print("HEIGHT_LAST|",H[30,30]);',
            'print("PARI|",version());',
            "quit",
        ]
    )
    completed = subprocess.run(
        ["gp", "-q"],
        input="\n".join(commands) + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=120,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"PARI/GP stderr: {completed.stderr.strip()}")

    rows = {}
    for match in LOCAL_PATTERN.finditer(completed.stdout):
        prime = int(match.group(1))
        change = tuple(item.strip() for item in match.group(4).split(","))
        if change != ("1", "0", "0", "0"):
            raise AssertionError(f"the public model is not minimal at {prime}")
        rows[prime] = (
            int(match.group(2)),
            int(match.group(3)),
            int(match.group(5)),
            int(match.group(6)),
        )
    if rows != EXPECTED_LOCAL:
        raise AssertionError(f"local-reduction fingerprint changed: {rows}")
    primality = {
        int(prime): int(status)
        for prime, status in re.findall(r"^PRIME\|(\d+)\|(\d+)$", completed.stdout, re.MULTILINE)
    }
    if primality != {prime: 1 for prime in EXPECTED_LOCAL}:
        raise AssertionError("a pinned bad-prime factor failed PARI isprime")

    def field(label: str) -> str:
        match = re.search(rf"^{label}\|(.+)$", completed.stdout, re.MULTILINE)
        if match is None:
            raise AssertionError(f"PARI/GP omitted {label}")
        return match.group(1).strip()

    root_number = -1
    for _conductor_exponent, _kodaira, _tamagawa, local_sign in rows.values():
        root_number *= local_sign
    isogeny_class_size = int(field("ISOGENY_CLASS_SIZE"))
    isogeny_degree_matrix = field("ISOGENY_DEGREE_MATRIX")
    if isogeny_class_size != 1 or isogeny_degree_matrix != "Mat(1)":
        raise AssertionError(
            "curve 398 no longer has a singleton rational isogeny class"
        )
    return {
        "version": field("PARI"),
        "local": rows,
        "root_number": root_number,
        "regulator_numeric_100_digits": field("REGULATOR"),
        "first_basis_height_numeric_100_digits": field("HEIGHT_FIRST"),
        "last_basis_height_numeric_100_digits": field("HEIGHT_LAST"),
        "rational_isogeny_class_size": isogeny_class_size,
        "rational_isogeny_degree_matrix": isogeny_degree_matrix,
    }


def build_certificate() -> dict[str, object]:
    if len(POINTS) != 30 or len(set(POINTS)) != 30:
        raise AssertionError("expected thirty distinct public points")
    if not all(on_curve(point) for point in POINTS):
        raise AssertionError("a public point is off curve 398")
    short_a = short_coefficients()[3]
    short_b = short_coefficients()[4]
    for x_value, y_value in SHORT_POINTS:
        if y_value**2 != x_value**3 + short_a * x_value + short_b:
            raise AssertionError("the integral short-model transport failed")
    print(f"{PROTOCOL}|stage=membership|checked=30|status=PASS", flush=True)

    if discriminant() != PUBLIC_DISCRIMINANT:
        raise AssertionError("the public discriminant changed")
    if factor_product(DISCRIMINANT_FACTORIZATION) != PUBLIC_DISCRIMINANT:
        raise AssertionError("the discriminant factorization changed")
    if factor_product(tuple((prime, 1) for prime in EXPECTED_LOCAL)) != PUBLIC_CONDUCTOR:
        raise AssertionError("the conductor factorization changed")

    short = short_coefficients()
    two_torsion_prime = find_two_torsion_certificate_prime(short, prime_bound=500)
    if two_torsion_prime != 23:
        raise AssertionError(f"least no-2-torsion witness changed to {two_torsion_prime}")
    torsion_counts = {11: finite_group_order(11), 23: finite_group_order(23)}
    if torsion_counts != {11: 18, 23: 31} or gcd(*torsion_counts.values()) != 1:
        raise AssertionError(f"finite torsion witnesses changed: {torsion_counts}")

    signatures = find_mod2_reduction_certificate(short, SHORT_POINTS, prime_bound=2000)
    certificate_primes = tuple(signature.prime for signature in signatures)
    rank = combined_mod2_rank(signatures, 30)
    row_count = sum(len(signature.rows) for signature in signatures)
    if certificate_primes != EXPECTED_CERTIFICATE_PRIMES or rank != 30 or row_count != 31:
        raise AssertionError(
            f"finite-reduction certificate changed: primes={certificate_primes}, "
            f"rank={rank}, rows={row_count}"
        )
    print(
        f"{PROTOCOL}|stage=mod2|rank=30|rows=31|primes=25|status=PASS",
        flush=True,
    )

    incidence = node_incidence()
    if incidence != EXPECTED_NODE_INCIDENCE:
        raise AssertionError(f"bad-prime node-incidence fingerprint changed: {incidence}")
    a1_exclusion = verify_known_a1_exclusion()
    print(
        f"{PROTOCOL}|stage=construction|known_a1_modp=179|projective_roots=0|status=PASS",
        flush=True,
    )

    pari = pari_diagnostics()
    if pari["root_number"] != 1:
        raise AssertionError(f"global root number changed to {pari['root_number']}")
    if pari["rational_isogeny_class_size"] != 1:
        raise AssertionError("curve 398 acquired a nontrivial rational isogeny")
    print(
        f"{PROTOCOL}|stage=local|bad_primes=18|root_number=1|status=PASS",
        flush=True,
    )

    integral_pairs = sum(
        x_value.denominator == y_value.denominator == 1 for x_value, y_value in POINTS
    )
    return {
        "schema_version": 1,
        "artifact_kind": "exact_elliptic_curve_rank_lower_bound_and_construction_dissection",
        "curve_id": "icarm_curve_398",
        "claim": "rank E(Q) >= 30",
        "claim_status": "exact unconditional lower bound; no unconditional exact-rank claim",
        "public_sources": {
            "curve": PUBLIC_CURVE_SOURCE,
            "json": PUBLIC_JSON_SOURCE,
            "json_sha256_retrieved_2026_09_04": PUBLIC_JSON_SHA256_2026_09_04,
            "primary_construction_announcement": PRIMARY_ANNOUNCEMENT,
        },
        "curve": {
            "ainvs": [str(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS],
            "c4": str(c4()),
            "c6": str(c6()),
            "discriminant": str(PUBLIC_DISCRIMINANT),
            "discriminant_factorization": [list(item) for item in DISCRIMINANT_FACTORIZATION],
            "conductor": str(PUBLIC_CONDUCTOR),
            "conductor_factorization": [[str(prime), 1] for prime in EXPECTED_LOCAL],
            "global_minimal_model": True,
            "root_number": pari["root_number"],
            "torsion_order": 1,
        },
        "points": [[str(x_value), str(y_value)] for x_value, y_value in POINTS],
        "point_membership_checks": 30,
        "displayed_integral_point_pairs": integral_pairs,
        "independence_certificate": {
            "method": "finite good-reduction quotients E(F_p)/2E(F_p)",
            "relation_prime": 2,
            "no_rational_2_torsion_witness_prime": two_torsion_prime,
            "full_torsion_witness_group_orders": torsion_counts,
            "combined_binary_rank": rank,
            "rows": [
                {
                    "prime": signature.prime,
                    "group_order": signature.group_order,
                    "doubled_subgroup_order": signature.doubled_subgroup_order,
                    "quotient_dimension": signature.quotient_dimension,
                    "matrix_rows": [list(row) for row in signature.rows],
                }
                for signature in signatures
            ],
        },
        "local_reduction": [
            {
                "prime": str(prime),
                "discriminant_valuation": dict(DISCRIMINANT_FACTORIZATION)[prime],
                "conductor_exponent": data[0],
                "pari_kodaira_code": data[1],
                "kodaira": f"I{data[1] - 4}",
                "tamagawa_number": data[2],
                "local_root_number": data[3],
                "public_basis_points_reducing_to_node": incidence[prime][0],
                "point_indices_unavailable_at_prime": list(incidence[prime][1]),
            }
            for prime, data in EXPECTED_LOCAL.items()
        ],
        "height_diagnostic": {
            "regulator_numeric_100_digits": pari["regulator_numeric_100_digits"],
            "first_basis_height_numeric_100_digits": pari[
                "first_basis_height_numeric_100_digits"
            ],
            "last_basis_height_numeric_100_digits": pari[
                "last_basis_height_numeric_100_digits"
            ],
            "used_for_rank_claim": False,
        },
        "rational_isogeny_certificate": {
            "method": "PARI ellisomat over Q",
            "isomorphism_class_count": pari["rational_isogeny_class_size"],
            "minimal_isogeny_degree_matrix": pari[
                "rational_isogeny_degree_matrix"
            ],
            "conclusion": (
                "The Q-isogeny class contains only curve 398's Q-isomorphism "
                "class; there is no nontrivial rational isogeny to another curve."
            ),
        },
        "construction_provenance": {
            "publicly_resolved": (
                "Elkies and Klagsbrun found the curve in September 2025 by "
                "searching elliptic fibrations of the determinant-948 X948 K3 "
                "with one I2 (or III) fibre and generic Mordell--Weil rank 16."
            ),
            "search_description": (
                "The announcement says hundreds of such fibrations were searched "
                "with an improved codebase using Drew Sutherland's smalljac library."
            ),
            "exact_fibration_parameter_and_section_map": "NOT_PUBLIC_UNKNOWN",
            "previous_exact_exclusions": (
                "The repository already excludes the published R17 chart and all six "
                "known norm-twelve chart classes by exact j-recognition."
            ),
            "known_equation_explicit_a1_test": a1_exclusion,
            "claim_boundary": (
                "The modulo-179 test excludes only the repository's one equation-explicit "
                "fixed-corridor A1/MW16 fibration. It does not identify or exhaust the "
                "hundreds of A1 fibrations used by Elkies and Klagsbrun."
            ),
        },
        "generation": {
            "command": REPRODUCING_COMMAND,
            "python": sys.version.split()[0],
            "pari_gp": pari["version"],
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "model_data_sha256": sha256(
                Path(__file__).with_name("icarm_curve398.py").read_bytes()
            ).hexdigest(),
            "known_a1_reduction_sha256": sha256(A1_REDUCTION.read_bytes()).hexdigest(),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale pinned certificate: rerun {REPRODUCING_COMMAND}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(
        f"{PROTOCOL}|stage=done|rank_lower_bound=30|exact_rank=false|"
        f"mode={'check' if args.check else 'write'}|output={args.output}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
