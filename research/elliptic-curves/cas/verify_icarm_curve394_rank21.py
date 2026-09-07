#!/usr/bin/env python3
"""Replay ICARM curve 394 as an exact rank-at-least-21 conductor anchor.

The compact Elkies R17 family is specialized at t=3/8 and globally minimized
to the public model.  The seventeen exact generic sections plus four public
points receive one finite-reduction independence certificate of rank 21.
PARI is used only for exact minimalization, primality, local Tate data, and the
root number; it is not used for the rank decision.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS)]

from analyze_icarm_7fff_zip_sequence import kodaira_symbol  # noqa: E402
from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    source_point_to_target,
    verify_finite_quotient_certificate,
    weierstrass_invariants,
)
from pari_bridge import pari_version  # noqa: E402
from search_extra_points import gp_rational, run_gp  # noqa: E402


Q = Fraction
PROTOCOL = "R21ICARM394"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve394_rank21_v1.json"
)
REPRODUCING_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/verify_icarm_curve394_rank21.py"
)
PUBLIC_SOURCE = "https://elliptic-rank.icarm.cloud/curve/394"


PUBLIC_MODEL = (
    Q(1),
    Q(0),
    Q(0),
    Q(-354803089674674467206754048738095),
    Q(2558194545892203175112161719607326368645810580537),
)
PUBLIC_POINTS = tuple(
    (Q(x), Q(y))
    for x, y in (
        ("-21518735422155486", "-478256579137186088822007"),
        ("-16162051516824546", "-2017627317695268505559667"),
        ("1038924778080097686/49", "-733011743810123946757418901/343"),
        ("-14545633423847586", "-2154421321999041799859307"),
        ("-10943447218123569/4", "14984597983246077098605569/8"),
        ("10125440811318864", "61344439585521052023093"),
        ("789701034488124", "-1509469661698186222233897"),
        ("228385634227786349/16", "-40580603327174753093901823/64"),
        ("-21627184418104416", "340283723703987370480203"),
        ("8007779902166814", "-480109043182118745601707"),
        ("12800689418722964", "337581497100471707837693"),
        ("9953363946011664", "113068979322481106243643"),
        ("11790929866630764", "118246915756400740772493"),
        ("945820846806594", "-1491127257554750266981227"),
        ("16992684482716974", "-1198247797994777274230547"),
        ("20859855244894382", "-2057638595340369190119955"),
        ("13765446221801214", "531556029601305747519693"),
        ("-2970309185321538", "-1893637547361814114859787"),
        ("-21725671106393538", "109128648895082432309325"),
        ("11620482431468214", "-66260339865015976889907"),
        ("67111185729163791/4", "9219933694407166168945929/8"),
    )
)
PUBLIC_COMPLEMENT_INDICES_ZERO_BASED = (0, 1, 2, 3)
EXPECTED_CERTIFICATE_PRIMES = (
    11,
    17,
    31,
    37,
    41,
    43,
    53,
    61,
    67,
    73,
    83,
    97,
    101,
    107,
    109,
    131,
    137,
    149,
    157,
)


PUBLIC_DISCRIMINANT = int(
    "31362810004008297563332830974145357962340308661810559854716023966067971553607881810277592152000000"
)
DISCRIMINANT_FACTORIZATION = (
    (2, 9),
    (3, 8),
    (5, 6),
    (7, 4),
    (13, 5),
    (23, 3),
    (29, 4),
    (89, 2),
    (43207, 1),
    (226549, 1),
    (22823593909227592035983291, 1),
    (44013936637595415741483513793, 1),
)
PUBLIC_CONDUCTOR = int(
    "1593562111507190066539814084004447718921281851572777685020200143306222910"
)
CONDUCTOR_FACTORIZATION = tuple((prime, 1) for prime, _ in DISCRIMINANT_FACTORIZATION)

# prime -> (f, PARI Kodaira code, Tamagawa number, local root number)
EXPECTED_LOCAL = {
    2: (1, 13, 9, -1),
    3: (1, 12, 8, -1),
    5: (1, 10, 6, -1),
    7: (1, 8, 4, -1),
    13: (1, 9, 5, -1),
    23: (1, 7, 3, -1),
    29: (1, 8, 4, -1),
    89: (1, 6, 2, 1),
    43207: (1, 5, 1, -1),
    226549: (1, 5, 1, -1),
    22823593909227592035983291: (1, 5, 1, 1),
    44013936637595415741483513793: (1, 5, 1, -1),
}
LOCAL_PATTERN = re.compile(
    r"^LOCAL\|(\d+)\|(\d+)\|(-?\d+)\|\[([^]]+)\]\|(\d+)\|(-?\d+)$",
    re.MULTILINE,
)


def factor_product(factors: tuple[tuple[int, int], ...]) -> int:
    return prod(prime**exponent for prime, exponent in factors)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def local_reduction_replay() -> dict[str, object]:
    ainvs = ",".join(gp_rational(value) for value in PUBLIC_MODEL)
    commands = ["default(realprecision,80);", f"E=ellinit([{ainvs}]);"]
    for prime in EXPECTED_LOCAL:
        commands.extend(
            [
                f'print("PRIME|{prime}|",isprime({prime}));',
                f"L=elllocalred(E,{prime});",
                (
                    f'print("LOCAL|{prime}|",L[1],"|",L[2],"|",L[3],"|",'
                    f'L[4],"|",ellrootno(E,{prime}));'
                ),
            ]
        )
    commands.extend(
        [
            'print("MINIMAL|",Vec(ellminimalmodel(E))[1..5]);',
            'print("ROOT|",ellrootno(E));',
            f'print("LOGN|",log({PUBLIC_CONDUCTOR}));',
            "quit",
        ]
    )
    output, _wall = run_gp(
        "\n".join(commands) + "\n", timeout=60.0, stack_bytes=500_000_000
    )
    prime_rows = re.findall(r"^PRIME\|(\d+)\|(\d+)$", output, re.MULTILINE)
    if [(int(prime), int(flag)) for prime, flag in prime_rows] != [
        (prime, 1) for prime in EXPECTED_LOCAL
    ]:
        raise AssertionError("a declared discriminant factor failed exact primality")

    rows: dict[int, tuple[int, int, int, int]] = {}
    local_records = []
    valuation_by_prime = dict(DISCRIMINANT_FACTORIZATION)
    for match in LOCAL_PATTERN.finditer(output):
        prime = int(match.group(1))
        change = tuple(item.strip() for item in match.group(4).split(","))
        if change != ("1", "0", "0", "0"):
            raise AssertionError(f"public model is not minimal at {prime}: {change}")
        row = (
            int(match.group(2)),
            int(match.group(3)),
            int(match.group(5)),
            int(match.group(6)),
        )
        rows[prime] = row
        local_records.append(
            {
                "prime": str(prime),
                "discriminant_valuation": valuation_by_prime[prime],
                "conductor_exponent": row[0],
                "kodaira_code": row[1],
                "kodaira_symbol": kodaira_symbol(row[1]),
                "minimal_change": list(change),
                "tamagawa_number": row[2],
                "local_root_number": row[3],
            }
        )
    if rows != EXPECTED_LOCAL:
        raise AssertionError(f"local-reduction fingerprint changed: {rows}")

    minimal_match = re.search(r"^MINIMAL\|\[(.*?)\]$", output, re.MULTILINE)
    root_match = re.search(r"^ROOT\|(-?\d+)$", output, re.MULTILINE)
    log_match = re.search(r"^LOGN\|(\S+)$", output, re.MULTILINE)
    if minimal_match is None or root_match is None or log_match is None:
        raise AssertionError("PARI omitted a global diagnostic")
    minimal = tuple(item.strip() for item in minimal_match.group(1).split(","))
    if minimal != tuple(str(value) for value in PUBLIC_MODEL):
        raise AssertionError("PARI changed the public global minimal model")
    if int(root_match.group(1)) != -1:
        raise AssertionError("the global root number changed")
    return {
        "global_minimal_model": list(minimal),
        "local_reductions": local_records,
        "root_number": -1,
        "tamagawa_product": prod(row[2] for row in EXPECTED_LOCAL.values()),
        "log_conductor_numeric_80_digits": log_match.group(1),
    }


def exact_log_bound() -> dict[str, object]:
    decimal_digits = len(str(PUBLIC_CONDUCTOR))
    upper = Q(decimal_digits * 231, 100)
    target = Q(17325, 100)
    if not PUBLIC_CONDUCTOR < 10**decimal_digits or not upper < target:
        raise AssertionError("the elementary exact log-conductor bound failed")
    return {
        "conductor_less_than_power_of_ten": f"10^{decimal_digits}",
        "deduced_log_10_upper_bound": "231/100",
        "deduced_log_conductor_upper_bound": str(upper),
        "strict_target": str(target),
        "strict_target_proved_exactly": True,
    }


def build_certificate() -> dict[str, object]:
    sys.set_int_max_str_digits(0)
    invariants = weierstrass_invariants(PUBLIC_MODEL)
    if invariants["discriminant"] != PUBLIC_DISCRIMINANT:
        raise AssertionError("public discriminant does not match the model")
    if factor_product(DISCRIMINANT_FACTORIZATION) != PUBLIC_DISCRIMINANT:
        raise AssertionError("discriminant factorization is incomplete")
    if factor_product(CONDUCTOR_FACTORIZATION) != PUBLIC_CONDUCTOR:
        raise AssertionError("conductor factorization is incomplete")
    if len(PUBLIC_POINTS) != 21 or any(
        not is_on_weierstrass_curve(PUBLIC_MODEL, point) for point in PUBLIC_POINTS
    ):
        raise ArithmeticError("the public 21-point list failed exact replay")

    data = load_q12o5867_data(MODEL, SECTIONS)
    specialization = evaluate_projective_specialization(data, 3, 8)
    minimal_model, minimal_change, minimal_metadata = global_minimal_model_with_change(
        specialization.model
    )
    if minimal_model != PUBLIC_MODEL:
        raise ArithmeticError("the compact t=3/8 fibre missed ICARM curve 394")
    generic_points = tuple(
        source_point_to_target(point, minimal_change)
        for point in specialization.points
    )
    if len(generic_points) != 17 or any(
        not is_on_weierstrass_curve(PUBLIC_MODEL, point) for point in generic_points
    ):
        raise ArithmeticError("a specialized generic section missed the public model")

    combined_points = generic_points + tuple(
        PUBLIC_POINTS[index] for index in PUBLIC_COMPLEMENT_INDICES_ZERO_BASED
    )
    short_model, short_change = short_certificate_model(PUBLIC_MODEL)
    short_points = tuple(
        source_point_to_target(point, short_change) for point in combined_points
    )
    certificate = build_finite_quotient_certificate(
        short_model, short_points, relation_prime=2, prime_bound=1000
    )
    verify_finite_quotient_certificate(short_model, short_points, certificate)
    if (
        not certificate["certified_independent"]
        or certificate["combined_rank_over_relation_field"] != 21
        or tuple(certificate["certificate_primes"]) != EXPECTED_CERTIFICATE_PRIMES
    ):
        raise ArithmeticError("the rank-21 finite-reduction certificate changed")
    if certificate["torsion_witness"] != {"prime": 19, "group_order": 25}:
        raise ArithmeticError("the trivial-torsion witness changed")
    print(
        f"{PROTOCOL}|stage=rank|generic=17|public_complement=4|rank=21|status=PASS",
        flush=True,
    )

    local = local_reduction_replay()
    print(
        f"{PROTOCOL}|stage=conductor|bad_primes={len(EXPECTED_LOCAL)}"
        f"|logN={local['log_conductor_numeric_80_digits']}|status=PASS",
        flush=True,
    )
    script_path = Path(__file__)
    return {
        "schema": "elliptic-curves.icarm-curve394-rank21.v1",
        "status": "PASS_EXACT_ICARM_CURVE394_RANK21_CONDUCTOR_REPLAY",
        "claim": "rank E(Q) >= 21 and log conductor < 173.25",
        "claim_boundary": (
            "unconditional rank lower bound and exact conductor only; no "
            "Selmer upper bound or exact-rank claim"
        ),
        "public_source": PUBLIC_SOURCE,
        "compact_family_source": "arXiv:2608.25406 and pinned repository R17 model",
        "parameter": "3/8",
        "curve": {
            "global_minimal_model": [str(value) for value in PUBLIC_MODEL],
            "discriminant": str(PUBLIC_DISCRIMINANT),
            "discriminant_factorization": [
                [str(prime), exponent]
                for prime, exponent in DISCRIMINANT_FACTORIZATION
            ],
            "conductor": str(PUBLIC_CONDUCTOR),
            "conductor_factorization": [
                [str(prime), exponent] for prime, exponent in CONDUCTOR_FACTORIZATION
            ],
            "exact_log_bound": exact_log_bound(),
            **local,
        },
        "compact_specialization": {
            "exact_public_minimal_model_equality": True,
            "minimalization": minimal_metadata,
            "generic_section_count": 17,
            "generic_points_on_public_model": True,
        },
        "public_point_replay": {
            "count": len(PUBLIC_POINTS),
            "all_exactly_on_curve": True,
            "points": [[str(x), str(y)] for x, y in PUBLIC_POINTS],
        },
        "rank_lower_bound": {
            "certified_rank": 21,
            "basis": "17 compact generic sections plus public points 1,2,3,4",
            "public_complement_indices_one_based": [1, 2, 3, 4],
            "short_certificate_model": [str(value) for value in short_model],
            "combined_points": [[str(x), str(y)] for x, y in combined_points],
            "finite_quotient_independence": certificate,
        },
        "generation": {
            "command": REPRODUCING_COMMAND,
            "checker_sha256": file_sha256(script_path),
            "model_sha256": file_sha256(MODEL),
            "sections_sha256": file_sha256(SECTIONS),
            "pari_gp": pari_version(),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit(f"stale pinned certificate: rerun {REPRODUCING_COMMAND}")
        print(f"{PROTOCOL}|artifact={args.output}|status=PASS")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(f"{PROTOCOL}|artifact={args.output}|status=WROTE")


if __name__ == "__main__":
    main()
