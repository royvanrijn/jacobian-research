#!/usr/bin/env sage-python
"""Exactly specialize and certify the anchored-MW18 Nagao finalists.

Every H<=1000 survivor is reconstructed from the certified conic map.  The
seventeen direct R17 sections and the cover section are evaluated exactly on
the raw short fibre.  Integral independence is then certified, when possible,
by full column rank in a product of finite quotients
``E(F_p)/ell E(F_p)`` together with a good-reduction proof that
``E(Q)[ell]=0``.  Failure to find full rank within the declared prime bound is
recorded as UNKNOWN, never as dependence.

This script does not compute an upper rank bound, minimize the fibres, perform
a Selmer descent, or search for further points.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import QQ


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from mod_l_reduction_independence import (  # noqa: E402
    combined_mod_l_rank,
    find_mod_l_reduction_certificate,
    find_no_rational_l_torsion_prime,
    mod_l_reduction_signature,
)


CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json"
)
CAMPAIGN = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-nagao-h1000-summary-v1.json"
)
PUBLISHED_SECTIONS = (
    ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
)
MOD_L_HELPER = CAS / "mod_l_reduction_independence.py"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json"
)
CHECKPOINT_DIRECTORY = (
    ROOT / "artifacts/local/elkies-k3/r17-extreme-anchored-mw18-specializations-h1000"
)
EXPECTED_CAMPAIGN_STATUS = (
    "PASS_BOUNDED_HEURISTIC_EXTREME_ANCHORED_MW18_NAGAO_CAMPAIGN"
)
EXPECTED_COVER_STATUS = "PASS_EXACT_EXTREME_ANCHORED_MW18_COVERS"
MODULI = (2, 3, 5)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def qtext(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def fraction(value) -> Fraction:
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def polynomial_value(coefficients, value):
    result = QQ(0)
    for coefficient in reversed(coefficients):
        result = result * value + QQ(coefficient)
    return result


def rational_function_value(record, value):
    numerator = polynomial_value(
        record["numerator_coefficients_low_to_high"], value
    )
    denominator = polynomial_value(
        record["denominator_coefficients_low_to_high"], value
    )
    if denominator == 0:
        raise ZeroDivisionError("a rational section or conic map met its pole")
    return numerator / denominator


def rational_function_projective_value(record, numerator, denominator):
    if denominator:
        return rational_function_value(
            record, QQ(numerator) / QQ(denominator)
        )
    numerator_coefficients = record["numerator_coefficients_low_to_high"]
    denominator_coefficients = record["denominator_coefficients_low_to_high"]
    degree = max(len(numerator_coefficients), len(denominator_coefficients)) - 1
    top_numerator = (
        QQ(numerator_coefficients[degree])
        if degree < len(numerator_coefficients)
        else QQ(0)
    )
    top_denominator = (
        QQ(denominator_coefficients[degree])
        if degree < len(denominator_coefficients)
        else QQ(0)
    )
    if top_denominator == 0:
        raise ZeroDivisionError("a conic map sends r=infinity to base infinity")
    return top_numerator / top_denominator


def maximum_rational_bits(values) -> int:
    result = 0
    for value in values:
        value = QQ(value)
        result = max(
            result,
            abs(int(value.numerator())).bit_length(),
            int(value.denominator()).bit_length(),
        )
    return result


def point_record(point) -> dict[str, str]:
    return {"x": qtext(point[0]), "y": qtext(point[1])}


def signature_record(signature) -> dict[str, object]:
    return {
        "prime": signature.prime,
        "group_order": signature.group_order,
        "multiple_subgroup_order": signature.multiple_subgroup_order,
        "quotient_dimension": signature.quotient_dimension,
        "rows": [list(row) for row in signature.rows],
    }


def published_generic_values(parameter):
    payload = json.loads(PUBLISHED_SECTIONS.read_text())
    if payload.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_SECTIONS_AND_CHORDS":
        raise ArithmeticError("the published R17 section source is not exact")
    values = []
    for expected_index, record in enumerate(payload["sections"]):
        if int(record["basis_index"]) != expected_index:
            raise ArithmeticError("the published R17 sections changed order")
        x_value = polynomial_value(record["x_coefficients_low_to_high"], parameter)
        if expected_index == 0:
            y_value = polynomial_value(
                record["y_coefficients_low_to_high"], parameter
            )
        else:
            chord = record["chord"]
            reference = values[int(chord["reference_basis_index"])]
            slope = polynomial_value(
                chord["slope_coefficients_low_to_high"], parameter
            )
            y_value = reference[1] + slope * (x_value - reference[0])
        values.append((x_value, y_value))
    if len(values) != 17:
        raise ArithmeticError("the published R17 basis is incomplete")
    return values


def direct_generic_values(model, parameter):
    records = model.get("sections", {}).get("records", [])
    if len(records) != 17:
        raise ArithmeticError("a direct extreme chart has no exact R17 basis")
    values = []
    for expected_index, record in enumerate(records):
        if int(record["basis_index"]) != expected_index:
            raise ArithmeticError("a direct extreme-chart basis changed order")
        values.append(
            (
                rational_function_value(record["X"], parameter),
                rational_function_value(record["Y"], parameter),
            )
        )
    return values


def family_coefficients(model, parameter):
    source = model.get("weierstrass_model", model)
    return (
        polynomial_value(source["A_coefficients_low_to_high"], parameter),
        polynomial_value(source["B_coefficients_low_to_high"], parameter),
    )


def cover_by_label(certificate, label):
    historical = certificate["historical_rank28_anchor"]
    if historical["label"] == label:
        return historical, ROOT / historical["direct_model"], "historical-rank28"
    matches = []
    for chart in certificate["charts"]:
        for fibre in chart["fibres"]:
            for cover in fibre["covers"]:
                if cover["label"] == label:
                    matches.append(
                        (cover, ROOT / chart["direct_model"], f"curve-{fibre['curve_id']}")
                    )
    if len(matches) != 1:
        raise ArithmeticError(f"cover label {label} is absent or ambiguous")
    return matches[0]


def candidate_id(anchor_id: str, label: str, numerator: int, denominator: int) -> str:
    sign = "m" if numerator < 0 else "p"
    compact_label = label.replace("-orbit-", "-").replace("orbit-", "")
    return f"{anchor_id}-{compact_label}-{sign}{abs(numerator)}d{denominator}"


def load_context(args):
    certificate_bytes = args.certificate.read_bytes()
    campaign_bytes = args.campaign.read_bytes()
    certificate = json.loads(certificate_bytes)
    campaign = json.loads(campaign_bytes)
    if certificate.get("status") != EXPECTED_COVER_STATUS:
        raise ArithmeticError("the anchored cover certificate is not passing")
    if campaign.get("status") != EXPECTED_CAMPAIGN_STATUS:
        raise ArithmeticError("the H<=1000 Nagao campaign is not passing")
    if campaign.get("cover_certificate_sha256") != sha256(certificate_bytes).hexdigest():
        raise ArithmeticError("the Nagao campaign points to another cover certificate")
    ledger_paths = [ROOT / row["source"] for row in campaign["priority_order_by_top_nagao_score"]]
    if len(ledger_paths) != 9 or len(set(ledger_paths)) != 9:
        raise ArithmeticError("the campaign does not bind nine distinct cover ledgers")
    ledgers = []
    for path, summary_row in zip(ledger_paths, campaign["priority_order_by_top_nagao_score"]):
        if digest(path) != summary_row["source_sha256"]:
            raise ArithmeticError(f"Nagao ledger hash drift: {relative(path)}")
        ledger = json.loads(path.read_text())
        if ledger.get("status") != "PASS_BOUNDED_HEURISTIC_EXTREME_ANCHORED_MW18_NAGAO_SIEVE":
            raise ArithmeticError("an input Nagao ledger is not passing")
        ledgers.append((path, ledger))
    requested = sum(len(ledger["finalists"]) for _path, ledger in ledgers)
    if requested != int(campaign["total_final_survivor_count"]) or requested != 178:
        raise ArithmeticError("the frozen campaign no longer has 178 finalists")
    return certificate, campaign, ledgers


def input_provenance(args, certificate, ledgers):
    paths = [
        args.certificate,
        args.campaign,
        PUBLISHED_SECTIONS,
        MOD_L_HELPER,
        Path(__file__).resolve(),
    ]
    paths.extend(path for path, _ledger in ledgers)
    paths.extend(
        ROOT / chart["direct_model"] for chart in certificate["charts"]
    )
    paths.append(ROOT / certificate["historical_rank28_anchor"]["direct_model"])
    unique = []
    for path in paths:
        path = path.resolve()
        if path not in unique:
            unique.append(path)
    return {relative(path): digest(path) for path in unique}


def specialize_candidate(*, cover, model_path, anchor_id, ledger, finalist, prime_bound):
    numerator, denominator = map(int, finalist["projective_pair"])
    r_value = None if denominator == 0 else QQ(numerator) / QQ(denominator)
    parameterization = cover["anchor_line_parameterization"]
    parameter = rational_function_projective_value(
        parameterization["t_of_r"], numerator, denominator
    )
    cover_value = rational_function_projective_value(
        parameterization["u_of_r"], numerator, denominator
    )
    if parameter != QQ(finalist["base_t"]):
        raise ArithmeticError("the exact conic map disagrees with the Nagao ledger")
    branch_value = polynomial_value(
        cover["branch_quadratic_coefficients_low_to_high"], parameter
    )
    if cover_value**2 != branch_value:
        raise ArithmeticError("the exact candidate missed its cover conic")

    model = json.loads(model_path.read_text())
    coefficient_a, coefficient_b = family_coefficients(model, parameter)
    if 4 * coefficient_a**3 + 27 * coefficient_b**2 == 0:
        raise ArithmeticError("the finalist specializes to a singular fibre")
    if model.get("sections", {}).get("records"):
        generic_points = direct_generic_values(model, parameter)
    else:
        generic_points = published_generic_values(parameter)
    new = cover["eighteenth_section"]
    new_point = (
        polynomial_value(new["x0_coefficients_low_to_high"], parameter)
        + polynomial_value(new["x1_coefficients_low_to_high"], parameter) * cover_value,
        polynomial_value(new["y0_coefficients_low_to_high"], parameter)
        + polynomial_value(new["y1_coefficients_low_to_high"], parameter) * cover_value,
    )
    points = tuple(generic_points) + (new_point,)
    for point in points:
        if point[1] ** 2 != point[0] ** 3 + coefficient_a * point[0] + coefficient_b:
            raise ArithmeticError("an exact specialized section missed the raw fibre")

    short_model = (Fraction(0), Fraction(0), Fraction(0), fraction(coefficient_a), fraction(coefficient_b))
    rational_points = tuple((fraction(x), fraction(y)) for x, y in points)
    attempts = []
    selected_modulus = None
    best_rank = 0
    for modulus in MODULI:
        signatures = find_mod_l_reduction_certificate(
            short_model,
            rational_points,
            modulus=modulus,
            prime_bound=prime_bound,
        )
        finite_rank = combined_mod_l_rank(
            signatures, len(rational_points), modulus
        )
        best_rank = max(best_rank, finite_rank)
        try:
            torsion_prime = find_no_rational_l_torsion_prime(
                short_model, modulus=modulus, prime_bound=min(prime_bound, 200)
            )
            torsion_group_order = mod_l_reduction_signature(
                short_model, (), torsion_prime, modulus
            ).group_order
            torsion_record = {
                "certified": True,
                "prime": torsion_prime,
                "group_order": torsion_group_order,
                "group_order_not_divisible_by_modulus": True,
            }
        except ValueError:
            torsion_record = {
                "certified": False,
                "prime_bound": min(prime_bound, 200),
            }
        attempts.append(
            {
                "modulus": modulus,
                "prime_bound": prime_bound,
                "finite_quotient_rank": finite_rank,
                "certificate_primes": [signature.prime for signature in signatures],
                "signatures": [signature_record(signature) for signature in signatures],
                "no_rational_modulus_torsion": torsion_record,
                "integral_independence_certified": (
                    finite_rank == len(points) and torsion_record["certified"]
                ),
            }
        )
        if attempts[-1]["integral_independence_certified"]:
            selected_modulus = modulus
            break

    identifier = candidate_id(anchor_id, cover["label"], numerator, denominator)
    certified = selected_modulus is not None
    all_coordinates = [coordinate for point in points for coordinate in point]
    return {
        "candidate_id": identifier,
        "anchor_id": anchor_id,
        "cover_label": cover["label"],
        "r": "infinity" if r_value is None else qtext(r_value),
        "r_projective_pair": [numerator, denominator],
        "base_t": qtext(parameter),
        "cover_u": qtext(cover_value),
        "nagao": finalist,
        "raw_short_model": ["0", "0", "0", qtext(coefficient_a), qtext(coefficient_b)],
        "raw_model_maximum_coefficient_bits": maximum_rational_bits(
            (coefficient_a, coefficient_b)
        ),
        "specialized_points": {
            "generic_R17": [point_record(point) for point in generic_points],
            "cover_section": point_record(new_point),
            "all_section_identities_verified_exactly": True,
            "maximum_coordinate_bits": maximum_rational_bits(all_coordinates),
        },
        "independence": {
            "method": (
                "full column rank in a product of E(F_p)/ell E(F_p), plus a "
                "good-reduction group order prime to ell proving E(Q)[ell]=0"
            ),
            "attempts": attempts,
            "selected_modulus": selected_modulus,
            "finite_mod_l_certified_rank_lower_bound": 18 if certified else best_rank,
            "status": (
                "CERTIFIED_INTEGRALLY_INDEPENDENT_RANK_AT_LEAST_18"
                if certified
                else "UNKNOWN_FULL_INDEPENDENCE_WITHIN_DECLARED_PRIME_BOUNDS"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--checkpoint-directory", type=Path, default=CHECKPOINT_DIRECTORY)
    parser.add_argument("--prime-bound", type=int, default=500)
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.prime_bound < 50:
        parser.error("--prime-bound must be at least 50")

    certificate, campaign, ledgers = load_context(args)
    provenance = input_provenance(args, certificate, ledgers)
    input_key = sha256(
        json.dumps(provenance, separators=(",", ":"), sort_keys=True).encode()
    ).hexdigest()
    records = []
    failures = []
    completed = 0
    total = sum(len(ledger["finalists"]) for _path, ledger in ledgers)
    args.checkpoint_directory.mkdir(parents=True, exist_ok=True)
    for _ledger_path, ledger in ledgers:
        cover, model_path, expected_anchor = cover_by_label(
            certificate, ledger["cover"]["label"]
        )
        if ledger["cover"]["anchor_id"] != expected_anchor:
            raise ArithmeticError("a Nagao ledger moved to another anchor")
        for finalist in ledger["finalists"]:
            numerator, denominator = map(int, finalist["projective_pair"])
            identifier = candidate_id(
                expected_anchor, cover["label"], numerator, denominator
            )
            checkpoint = args.checkpoint_directory / f"{identifier}.json"
            record = None
            if not args.no_resume and checkpoint.is_file():
                saved = json.loads(checkpoint.read_text())
                if saved.get("input_key") == input_key and saved.get("prime_bound") == args.prime_bound:
                    record = saved.get("record")
            if record is None:
                try:
                    record = specialize_candidate(
                        cover=cover,
                        model_path=model_path,
                        anchor_id=expected_anchor,
                        ledger=ledger,
                        finalist=finalist,
                        prime_bound=args.prime_bound,
                    )
                    checkpoint.write_text(
                        json.dumps(
                            {
                                "input_key": input_key,
                                "prime_bound": args.prime_bound,
                                "record": record,
                            },
                            separators=(",", ":"),
                            sort_keys=True,
                        )
                        + "\n"
                    )
                except (ArithmeticError, ValueError, ZeroDivisionError) as error:
                    failures.append(
                        {
                            "candidate_id": identifier,
                            "reason": type(error).__name__,
                            "detail": str(error),
                        }
                    )
                    completed += 1
                    print(
                        f"R17ANCHORMW18SPECIALIZE|completed={completed}/{total}|"
                        f"id={identifier}|status=FAIL|reason={type(error).__name__}",
                        flush=True,
                    )
                    continue
            records.append(record)
            completed += 1
            print(
                f"R17ANCHORMW18SPECIALIZE|completed={completed}/{total}|"
                f"id={identifier}|modulus={record['independence']['selected_modulus']}|"
                f"rank={record['independence']['finite_mod_l_certified_rank_lower_bound']}|"
                "status=COMPLETE",
                flush=True,
            )

    certified = sum(
        record["independence"]["status"]
        == "CERTIFIED_INTEGRALLY_INDEPENDENT_RANK_AT_LEAST_18"
        for record in records
    )
    payload = {
        "schema": "elkies-k3.r17-extreme-anchored-mw18-specializations.v1",
        "status": "COMPLETE_EXACT_MW18_FINALIST_SPECIALIZATION_AUDIT",
        "requested_finalist_count": total,
        "successful_specialization_count": len(records),
        "structural_failure_count": len(failures),
        "certified_rank_at_least_18_count": certified,
        "independence_unknown_count": len(records) - certified,
        "prime_bound": args.prime_bound,
        "moduli_tried_in_order": list(MODULI),
        "candidates": records,
        "structural_failures": failures,
        "inputs": provenance,
        "input_key": input_key,
        "next_gate": {
            "stage": "residual_2_selmer_filter",
            "eligible_candidate_count": certified,
            "rank32_required_residual_dimension_over_certified_MW18": 14,
            "point_search_before_residual_gate": False,
        },
        "claim_boundary": [
            "Every successful row contains eighteen exact rational points on its raw specialized short curve.",
            "A full finite mod-ell certificate proves those eighteen points integrally independent and hence rank at least 18.",
            "A rank-deficient bounded finite-reduction attempt is UNKNOWN, not a dependence or rank upper bound.",
            "No row proves that the generic cover rank or specialized Mordell-Weil rank is exactly 18.",
            "No minimization, complete Selmer computation, or further-point search occurs here.",
        ],
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "specialize_r17_extreme_anchored_mw18_finalists.sage --check --no-resume"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != serialized:
            raise ArithmeticError("stored specialization audit differs from fresh replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"R17ANCHORMW18SPECIALIZE|requested={total}|successful={len(records)}|"
        f"certified_rank18={certified}|unknown={len(records)-certified}|"
        f"failures={len(failures)}|output={relative(args.output)}|status=COMPLETE",
        flush=True,
    )


if __name__ == "__main__":
    main()
