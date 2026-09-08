#!/usr/bin/env python3
"""Replay Fermigier's 22 published points on his record curve.

The point coordinates below are transcribed from Theorem 1 on pages 362--363
of Fermigier's 1997 paper and cross-checked against Dujella's record table.
This verifier performs exact curve-membership checks over ``Fraction`` and
asks PARI/GP to repeat exact membership plus a high-precision numerical
Neron--Tate height-pairing calculation.  It also transports the printed points
to the normalized short model and proves their independence exactly: their
images in a finite product of good-reduction quotients
``E(F_p)/2E(F_p)`` have full binary column rank, while a separate good prime
proves ``E(Q)[2]=0``.  Infinite descent then gives an unconditional,
portable rank-at-least-22 certificate independent of the numerical regulator.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import platform
import shutil
import subprocess
from typing import Any

from ek_k3 import rational_to_string
from fermigier_mestre import (
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
)
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from pari_bridge import pari_version
from verify_fermigier_benchmark import PUBLISHED_CONDUCTOR, PUBLISHED_MODEL


PRIMARY_PAPER = {
    "author": "Stephane Fermigier",
    "title": "Une courbe elliptique definie sur Q de rang >= 22",
    "journal": "Acta Arithmetica 82.4 (1997), 359--363",
    "doi": "10.4064/aa-82-4-359-363",
    "url": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
    "theorem_and_point_pages": "362--363",
    "downloaded_pdf_sha256": (
        "9e0455228382c74b0e558b80b28346d5440531eb797b2539d3d379f1c86d77e4"
    ),
}
SECONDARY_CROSS_CHECK = (
    "https://web.math.pmf.unizg.hr/~duje/tors/rk22.html"
)
PAPER_DETERMINANT_APPROX = Decimal("1.299202e22")
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_fermigier_rank22_points.py"
)

# PARI's exact ``ellminimalmodel`` change from the normalized short model to
# Fermigier's printed minimal model.  For ``v=(u,r,s,t)``, the inverse point
# change is
#
#   x_short = u^2*x_min + r,
#   y_short = u^3*y_min + s*u^2*x_min + t.
SHORT_TO_MINIMAL_CHANGE = (
    Fraction(14, 507),
    Fraction(49, 771147),
    Fraction(7, 507),
    Fraction(1372, 130323843),
)
EXPECTED_CERTIFICATE_PRIMES = (
    29,
    43,
    67,
    73,
    79,
    83,
    89,
    101,
    103,
    107,
    109,
    127,
    131,
    137,
    149,
    191,
    223,
)
EXPECTED_TWO_TORSION_CERTIFICATE_PRIME = 31
STRICT_LOG_CONDUCTOR_TARGET = Fraction(4568, 25)
E_RATIONAL_UPPER_BOUND = Fraction(1359141, 500000)


def q(value: str) -> Fraction:
    """Parse a source coordinate as an exact rational number."""

    return Fraction(value)


# Keep the source order P_1,...,P_22.  The strings make visual comparison
# with the typeset theorem possible and avoid any decimal interpretation.
PUBLISHED_POINT_STRINGS = (
    (
        "32741153161482344264/3025",
        "-223089674587110979578532169697/166375",
    ),
    (
        "215521674613198983365/24649",
        "-6872949155061353554235704378947/3869893",
    ),
    (
        "637312541911044643/81",
        "-1420356190129296832193564087/729",
    ),
    (
        "-11906250919327880080/361",
        "-16580788535875788634285886853/6859",
    ),
    ("-136152345735493381/4", "-14482270545045735913281693/8"),
    (
        "-27830298157016213012252/7134241",
        "72099692861364392796183359497454267/19055557711",
    ),
    ("4127671322151440", "2626107692045613116291646"),
    ("6175679781777296", "2266254335997033124678449"),
    ("12047255022287093", "1061993236525943920980477"),
    (
        "416685837455186583191/32761",
        "5321268222786709669160311587369/5929741",
    ),
    (
        "149915813139075767108024/10220809",
        "8704326838108646949177663157917117/32675926373",
    ),
    ("58759417448623559/4", "2030968553150713398654657/8"),
    (
        "237195157887349854919517/16024009",
        "-11477798111611307979707215505421441/64144108027",
    ),
    (
        "9568474434078537574436/687241",
        "319520556343135681977874272805086/569722789",
    ),
    (
        "1725892668710258675291/177241",
        "117378050663464845770966453025039/74618461",
    ),
    (
        "-35277008506980340471/1024",
        "48766027143946934186731674507/32768",
    ),
    (
        "-2752742763529705669/121",
        "6000532252185982381233585699/1331",
    ),
    ("-18552633109178014", "-4665466215824339436717966"),
    (
        "-113251707338691187737649969/3304065361",
        "310152527894831470820009872373229341739/189920981015641",
    ),
    (
        "-7572001778163591251/729",
        "-86590661426506799357663502953/19683",
    ),
    (
        "-380526048554032285152211/11242609",
        "73081235744931307684790623068490233/37696467977",
    ),
    (
        "-1503889497722021588110681/42784681",
        "-160705885170116750151534640924719585/279854598421",
    ),
)
PUBLISHED_POINTS = tuple(
    (q(x_value), q(y_value)) for x_value, y_value in PUBLISHED_POINT_STRINGS
)


def curve_residual(point: tuple[Fraction, Fraction]) -> Fraction:
    """Return LHS minus RHS of the generalized Weierstrass equation."""

    a1, a2, a3, a4, a6 = PUBLISHED_MODEL
    x_value, y_value = point
    return (
        y_value**2
        + a1 * x_value * y_value
        + a3 * y_value
        - x_value**3
        - a2 * x_value**2
        - a4 * x_value
        - a6
    )


def minimal_point_to_short(
    point: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    """Apply the exact inverse minimal-model point change."""

    u, r, s, t = SHORT_TO_MINIMAL_CHANGE
    x_minimal, y_minimal = point
    x_short = u**2 * x_minimal + r
    y_short = u**3 * y_minimal + s * u**2 * x_minimal + t
    return x_short, y_short


def exact_independence_certificate() -> dict[str, Any]:
    """Certify the 22 printed points by finite good reductions."""

    coefficients = FermigierMestreFamily.coefficients(
        NORMALIZED_RECORD_PARAMETER
    )
    short_points = tuple(minimal_point_to_short(point) for point in PUBLISHED_POINTS)
    _, _, _, coefficient_a, coefficient_b = coefficients
    if any(
        y_value**2 != x_value**3 + coefficient_a * x_value + coefficient_b
        for x_value, y_value in short_points
    ):
        raise AssertionError("the exact minimal-to-short point change failed")

    signatures = find_mod2_reduction_certificate(
        coefficients, short_points, prime_bound=500
    )
    exact_rank = combined_mod2_rank(signatures, len(short_points))
    certificate_primes = tuple(signature.prime for signature in signatures)
    if exact_rank != 22:
        raise AssertionError("finite reductions did not certify all 22 points")
    if certificate_primes != EXPECTED_CERTIFICATE_PRIMES:
        raise AssertionError("the deterministic reduction certificate changed")
    two_torsion_prime = find_two_torsion_certificate_prime(coefficients)
    if two_torsion_prime != EXPECTED_TWO_TORSION_CERTIFICATE_PRIME:
        raise AssertionError("the rational 2-torsion certificate prime changed")

    return {
        "normalized_record_parameter": rational_to_string(
            NORMALIZED_RECORD_PARAMETER
        ),
        "short_weierstrass_coefficients": [
            rational_to_string(value) for value in coefficients
        ],
        "minimal_to_short_change": [
            rational_to_string(value) for value in SHORT_TO_MINIMAL_CHANGE
        ],
        "transported_points": [
            {
                "label": f"P{index}",
                "x": rational_to_string(point[0]),
                "y": rational_to_string(point[1]),
                "exact_short_curve_membership_checked": True,
            }
            for index, point in enumerate(short_points, start=1)
        ],
        "two_torsion_certificate_prime": two_torsion_prime,
        "certificate_primes": list(certificate_primes),
        "finite_reduction_signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
        "combined_exact_rank_over_F2": exact_rank,
        "certified_algebraic_rank_lower_bound": exact_rank,
        "height_pairing_not_used_in_certificate": True,
        "argument": (
            "full binary column rank forces every integral relation coefficient "
            "to be even; E(Q)[2]=0 permits infinite descent, so every coefficient "
            "vanishes"
        ),
    }


def points_sha256() -> str:
    """Hash the canonical source-order coordinate transcription."""

    payload = json.dumps(PUBLISHED_POINT_STRINGS, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def exact_strict_conductor_comparison(conductor: int) -> dict[str, Any]:
    """Prove exactly that ``log(conductor)`` exceeds ``4568/25``.

    The exponential series through ``1/7!`` and the geometric tail bound

    ``sum(k>=8,1/k!) <= (1/8!)/(1-1/9) = 9/(8*8!)``

    give a rational upper bound for ``e``.  A single integer-power comparison
    then proves ``e^4568 < conductor^25``.  No floating-point logarithm enters
    the strict target decision.
    """

    if conductor <= 1:
        raise ValueError("the conductor must exceed one")
    series_bound = sum(Fraction(1, factorial(k)) for k in range(8))
    series_bound += Fraction(9, 8 * factorial(8))
    if not series_bound < E_RATIONAL_UPPER_BOUND:
        raise AssertionError("the pinned rational upper bound for e is invalid")
    power_inequality = (
        E_RATIONAL_UPPER_BOUND.numerator**4568
        < E_RATIONAL_UPPER_BOUND.denominator**4568 * conductor**25
    )
    if not power_inequality:
        raise AssertionError("the exact strict conductor comparison failed")
    return {
        "strict_log_target": rational_to_string(STRICT_LOG_CONDUCTOR_TARGET),
        "series_upper_bound_for_e": rational_to_string(series_bound),
        "rational_upper_bound_for_e": rational_to_string(E_RATIONAL_UPPER_BOUND),
        "series_bound_is_strictly_smaller": True,
        "integer_power_inequality": (
            "1359141^4568 < 500000^4568 * conductor^25"
        ),
        "integer_power_inequality_holds": power_inequality,
        "exact_conclusion": "log(conductor) > 4568/25 = 182.72",
        "meets_strict_log_conductor_target": False,
    }


def gp_rational(value: Fraction) -> str:
    return f"({rational_to_string(value)})"


def parse_precisions(value: str) -> tuple[int, ...]:
    try:
        precisions = tuple(int(item) for item in value.split(",") if item)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "precisions must be comma-separated integers"
        ) from error
    if len(precisions) < 2 or any(precision < 38 for precision in precisions):
        raise argparse.ArgumentTypeError(
            "provide at least two decimal precisions, each at least 38 digits"
        )
    if tuple(sorted(set(precisions))) != precisions:
        raise argparse.ArgumentTypeError("precisions must be unique and increasing")
    return precisions


def block(lines: list[str], name: str) -> list[str]:
    start = lines.index(f"{name}_BEGIN") + 1
    end = lines.index(f"{name}_END")
    return lines[start:end]


def pari_height_replay(
    *, precisions: tuple[int, ...], timeout: float
) -> dict[str, Any]:
    """Check the points and recompute the height matrix at each precision."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    model = ",".join(str(value) for value in PUBLISHED_MODEL)
    points = ",".join(
        f"[{gp_rational(x_value)},{gp_rational(y_value)}]"
        for x_value, y_value in PUBLISHED_POINTS
    )
    commands = [
        f"default(realprecision,{precisions[0]});",
        f"E=ellinit([{model}]);",
        f"P=[{points}];",
        "G=ellglobalred(E);",
        'print("CONDUCTOR_BEGIN");',
        "print(G[1]);",
        'print("CONDUCTOR_END");',
        'print("POINTS_BEGIN");',
        "print(#P);",
        "print(vecsum(vector(#P,i,ellisoncurve(E,P[i]))));",
        "print(vector(#P,i,ellisoncurve(E,P[i])));",
        'print("POINTS_END");',
    ]
    for precision in precisions:
        commands.extend(
            [
                f"default(realprecision,{precision});",
                "H=ellheightmatrix(E,P);",
                "D=matdet(H);",
                "EV=mateigen(H,1)[1];",
                f'print("HEIGHT_{precision}_BEGIN");',
                "print(D);",
                "print(matrank(H));",
                "print(vecmin(EV));",
                "print(vecmax(EV));",
                f'print("HEIGHT_{precision}_END");',
            ]
        )
    commands.append("quit")
    result = subprocess.run(
        [executable, "-q"],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP replay failed: {result.stderr.strip()}")
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    conductor = int(block(lines, "CONDUCTOR")[0])
    point_data = block(lines, "POINTS")
    count = int(point_data[0])
    on_curve_count = int(point_data[1])
    membership_vector = [
        int(value) for value in point_data[2].strip("[]").split(",")
    ]
    height_runs: list[dict[str, Any]] = []
    for precision in precisions:
        height_data = block(lines, f"HEIGHT_{precision}")
        height_runs.append(
            {
                "decimal_precision": precision,
                "determinant": height_data[0],
                "numerical_matrix_rank": int(height_data[1]),
                "smallest_eigenvalue": height_data[2],
                "largest_eigenvalue": height_data[3],
            }
        )
    return {
        "point_count": count,
        "on_curve_count": on_curve_count,
        "all_on_curve": count == on_curve_count == len(PUBLISHED_POINTS),
        "membership_vector": membership_vector,
        "conductor": conductor,
        "height_matrix_runs": height_runs,
    }


def decimal_comparison(height_runs: list[dict[str, Any]]) -> dict[str, str]:
    """Measure determinant stability and agreement with the paper's rounding."""

    determinants = [Decimal(run["determinant"]) for run in height_runs]
    with localcontext() as context:
        context.prec = max(run["decimal_precision"] for run in height_runs)
        final = determinants[-1]
        relative_stability = abs(determinants[-2] - final) / abs(final)
        relative_paper_difference = abs(final - PAPER_DETERMINANT_APPROX) / abs(final)
        return {
            "highest_precision_determinant": str(final),
            "highest_precision_scientific_7_significant_digits": f"{final:.6E}",
            "relative_change_between_last_two_precisions": str(relative_stability),
            "paper_approximation": str(PAPER_DETERMINANT_APPROX),
            "relative_difference_from_paper_rounded_value": str(relative_paper_difference),
        }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--precisions",
        type=parse_precisions,
        default=(96, 192),
        help="increasing comma-separated PARI decimal precisions",
    )
    parser.add_argument("--timeout", type=float, default=180.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic-curves"
        / "elliptic_fermigier_rank22_points.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.timeout <= 0:
        raise SystemExit("--timeout must be positive")
    if len(PUBLISHED_POINTS) != 22:
        raise AssertionError("the primary-source transcription must contain 22 points")

    residuals = tuple(curve_residual(point) for point in PUBLISHED_POINTS)
    if any(residual != 0 for residual in residuals):
        failures = [index for index, residual in enumerate(residuals, start=1) if residual]
        raise AssertionError(f"published points failed exact membership: {failures}")

    pari = pari_height_replay(precisions=args.precisions, timeout=args.timeout)
    if pari["conductor"] != PUBLISHED_CONDUCTOR:
        raise AssertionError("PARI did not reproduce Fermigier's conductor")
    if not pari["all_on_curve"] or pari["membership_vector"] != [1] * 22:
        raise AssertionError("PARI rejected at least one published point")
    for run in pari["height_matrix_runs"]:
        if run["numerical_matrix_rank"] != 22:
            raise AssertionError("PARI height matrix was not numerically full rank")
        if Decimal(run["determinant"]) <= 0 or Decimal(run["smallest_eigenvalue"]) <= 0:
            raise AssertionError("PARI height matrix was not numerically positive definite")

    comparison = decimal_comparison(pari["height_matrix_runs"])
    if comparison["highest_precision_scientific_7_significant_digits"] != "1.299202E+22":
        raise AssertionError("height determinant did not reproduce the paper's rounding")

    exact_certificate = exact_independence_certificate()
    if exact_certificate["certified_algebraic_rank_lower_bound"] != 22:
        raise AssertionError("the exact rank lower-bound certificate changed")
    strict_conductor_comparison = exact_strict_conductor_comparison(
        pari["conductor"]
    )

    point_records = [
        {
            "label": f"P{index}",
            "x": rational_to_string(point[0]),
            "y": rational_to_string(point[1]),
            "exact_equation_residual": rational_to_string(residual),
        }
        for index, (point, residual) in enumerate(
            zip(PUBLISHED_POINTS, residuals, strict=True), start=1
        )
    ]
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": {
            "published_theorem": (
                "Fermigier proves that P1,...,P22 are independent; hence "
                "rank E(Q) >= 22"
            ),
            "verified_exactly_here": (
                "all 22 transcribed rational coordinates satisfy the published "
                "curve equation, and exact finite good-reduction quotients "
                "independently certify their Z-independence; PARI replays the "
                "conductor and an integer-power inequality proves the strict "
                "conductor target is missed"
            ),
            "verified_numerically_here": (
                "the Neron--Tate height matrix is positive definite and has "
                "numerical rank 22 at two recorded precisions; its determinant "
                "reproduces Fermigier's published approximation"
            ),
            "not_claimed": (
                "the exact Mordell--Weil rank is not bounded above here; the "
                "certificate proves only rank at least 22"
            ),
        },
        "curve": {
            "model": list(PUBLISHED_MODEL),
            "equation": "y^2 + x*y + y = x^3 + a4*x + a6",
        },
        "conductor_replay_and_exact_target_comparison": {
            "conductor": pari["conductor"],
            **strict_conductor_comparison,
        },
        "published_rank_lower_bound": 22,
        "exact_rank_claim": None,
        "certified_rank_lower_bound_in_this_artifact": 22,
        "points": point_records,
        "point_transcription_sha256": points_sha256(),
        "exact_python_checks": {
            "point_count": len(PUBLISHED_POINTS),
            "all_on_curve": all(residual == 0 for residual in residuals),
            "arithmetic": "fractions.Fraction; no floating-point operations",
        },
        "pari_replay": pari,
        "determinant_comparison": comparison,
        "exact_independence_certificate": exact_certificate,
        "independence_classification": {
            "lower_bound_source": "published theorem in Fermigier 1997",
            "independent_replay": "high-precision numerical height pairing",
            "exact_certificate_generated_here": True,
            "exact_certificate_method": (
                "full-column-rank images in a product of E(F_p)/2E(F_p), "
                "plus a separate proof that E(Q)[2] is trivial"
            ),
            "ellsaturation_run": False,
            "ellsaturation_reason": (
                "not needed: the printed points themselves have a direct exact "
                "finite-reduction independence certificate"
            ),
        },
        "sources": {
            "primary": PRIMARY_PAPER,
            "secondary_coordinate_cross_check": SECONDARY_CROSS_CHECK,
            "pari_documentation": (
                "https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html"
            ),
        },
        "parameters": {
            "precisions": list(args.precisions),
            "timeout_seconds": args.timeout,
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

    final_run = pari["height_matrix_runs"][-1]
    print(f"wrote {args.output}")
    print("exact Python membership: 22/22")
    print("exact PARI membership: 22/22")
    print(
        f"height determinant ({final_run['decimal_precision']} digits): "
        f"{final_run['determinant']}"
    )
    print(f"numerical height-matrix rank: {final_run['numerical_matrix_rank']}")
    print("exact finite-reduction rank lower bound: 22")
    print("independence status: exact certificate plus published/numerical replay")


if __name__ == "__main__":
    main()
