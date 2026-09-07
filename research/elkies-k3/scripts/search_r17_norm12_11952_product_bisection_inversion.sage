#!/usr/bin/env sage-python
"""Invert alternate-Q80 product characters through norm-eight bisection pencils.

The input layer is the complete table of section-nonnegative degree-two
isotropic translation classes ``D_w=(2,2,w)``, ``(w,w)=8``.  For the trace
``tau=(Nx/h^2,Ny/h^3)`` attached to a representative, regular residual chord
slopes are

    M=M0+lambda*h^2.

The chord discriminant divided by ``h^6`` is a quartic in the old base and a
homogeneous quartic in the projective pencil parameter.  This script compares
that five-coefficient vector with each of the seventeen product quartics from
the exact rank-one V4 shortlist.

For speed, a class is rejected at a good prime only when *every* point of
``P1(F_p)`` fails the projective coefficient comparison.  This is an exact
characteristic-zero obstruction: a rational pencil parameter has a reduction
in ``P1(F_p)``, including at infinity and at a possible modular base point.
Any modular survivors are resolved by exact univariate coefficient gcds over
QQ.  A synthetic control quartic is first constructed from the first exact
norm-eight pencil and must be recovered by the same inversion path.

The search proves statements about the image of the bisection coboundary map
``P |-> P-sigma(P)``.  It does not silently assume that every integral twist
section is a coboundary; that separate 2-primary descent condition is recorded
explicitly in the output.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
PRIORITY_CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json"
PRIORITY_TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
SHORTLIST = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
RANKS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-bisection-inversion-v1.json"
DEFAULT_LEDGER = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-bisection-inversion-v1.tsv"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def coefficient_text(polynomial) -> list[str]:
    if not polynomial:
        return ["0"]
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def parse_vector(text: str):
    result = vector(ZZ, [ZZ(value) for value in text.split()])
    if len(result) != 17:
        raise ValueError("expected seventeen section coordinates")
    return result


def parse_polynomial(record, ring, coefficient_field=QQ):
    return ring([coefficient_field(QQ(value)) for value in record])


def parse_rational_function(record, ring, field, coefficient_field=QQ):
    numerator = parse_polynomial(
        record["numerator_coefficients_low_to_high"], ring, coefficient_field
    )
    denominator = parse_polynomial(
        record["denominator_coefficients_low_to_high"], ring, coefficient_field
    )
    return field(numerator / denominator)


def q_lambda_family(h, Nx, Ny, M0, A, ring):
    """Return q_0,...,q_4 for q(lambda)=sum q_j lambda^j."""

    h2 = h**2
    h6 = h**6
    numerators = [
        M0**4 - 6 * M0**2 * Nx - 8 * M0 * Ny - 3 * Nx**2 - 4 * A * h**4,
        4 * M0**3 * h2 - 12 * M0 * h2 * Nx - 8 * h2 * Ny,
        6 * M0**2 * h2**2 - 6 * h2**2 * Nx,
        4 * M0 * h6,
        h2 * h6,
    ]
    result = []
    for numerator in numerators:
        quotient, remainder = ring(numerator).quo_rem(h6)
        if remainder:
            raise ArithmeticError("pencil discriminant coefficient is not divisible by h^6")
        quotient = ring(quotient)
        if quotient.degree() > 4:
            raise ArithmeticError("pencil branch coefficient exceeds old-base degree four")
        result.append(quotient)
    if result[3] != 4 * M0 or result[4] != h2:
        raise ArithmeticError("closed-form leading pencil coefficients changed")
    return result


def polynomial_vector(polynomial, length=5):
    return tuple(polynomial[index] for index in range(length))


def projective_key(values):
    values = tuple(values)
    pivot = next((value for value in values if value), None)
    if pivot is None:
        return None
    inverse = pivot**-1
    return tuple(value * inverse for value in values)


def reciprocal_polynomial(polynomial, bound, ring):
    if polynomial.degree() > bound:
        raise ArithmeticError("reciprocal polynomial exceeds declared degree")
    u = ring.gen()
    return ring(
        sum(polynomial[index] * u ** (bound - index) for index in range(polynomial.degree() + 1))
    )


def exact_model_and_basis(direct, helper):
    ring = PolynomialRing(QQ, "u")
    field = ring.fraction_field()
    model = direct["weierstrass_model"]
    A = parse_polynomial(model["A_coefficients_low_to_high"], ring)
    B = parse_polynomial(model["B_coefficients_low_to_high"], ring)
    Delta = parse_polynomial(model["discriminant_coefficients_low_to_high"], ring)
    if Delta != -16 * (4 * A**3 + 27 * B**2):
        raise ArithmeticError("alternate-Q80 discriminant identity changed")
    curve = EllipticCurve(field, [A, B])
    basis = []
    for index, record in enumerate(direct["sections"]["records"]):
        if record["basis_index"] != index or record["equation_verified"] is not True:
            raise ArithmeticError("alternate section basis order/status changed")
        X = parse_rational_function(record["X"], ring, field)
        Y = parse_rational_function(record["Y"], ring, field)
        basis.append(curve(X, Y))
    if len(basis) != 17:
        raise ArithmeticError("alternate equation basis has wrong rank")
    return {
        "ring": ring,
        "field": field,
        "A": A,
        "B": B,
        "Delta": Delta,
        "curve": curve,
        "basis": basis,
        "helper": helper,
    }


def trace_from_word(context, word, multiples=None):
    curve = context["curve"]
    if multiples is None:
        return sum(
            (int(coefficient) * point for coefficient, point in zip(word, context["basis"])),
            curve(0),
        )
    return sum(
        (multiples[index][int(coefficient)] for index, coefficient in enumerate(word)),
        curve(0),
    )


def trace_pencil_family(context, trace):
    """Compile one trace in a chart with both trace poles finite."""

    ring = context["ring"]
    helper = context["helper"]
    X, Y = context["field"](trace[0]), context["field"](trace[1])
    frame = helper["trace_chord_frame"](X, Y, ring)
    chart = "finite"
    A, B, Delta = context["A"], context["B"], context["Delta"]
    if frame["h"].degree() != 2:
        chart = "inverted_at_infinity"
        A = reciprocal_polynomial(A, 8, ring)
        B = reciprocal_polynomial(B, 12, ring)
        Delta = reciprocal_polynomial(Delta, 24, ring)
        X = helper["invert_rational"](X, 4, ring, context["field"])
        Y = helper["invert_rational"](Y, 6, ring, context["field"])
        frame = helper["trace_chord_frame"](X, Y, ring)
    if frame["h"].degree() != 2:
        raise ArithmeticError("trace does not retain two finite poles in either base chart")
    family = q_lambda_family(
        frame["h"], frame["Nx"], frame["Ny"], frame["M0"], A, ring
    )
    return {
        "chart": chart,
        "A": A,
        "B": B,
        "Delta": Delta,
        **frame,
        "q_lambda_family": family,
    }


def build_positive_control(exact_context, first_row):
    word = parse_vector(first_row["section_basis_w"])
    trace = trace_from_word(exact_context, word)
    pencil = trace_pencil_family(exact_context, trace)
    helper = exact_context["helper"]
    for parameter in map(QQ, (0, 1, -1, 2, -2, 3, -3)):
        M = pencil["M0"] + parameter * pencil["h"]**2
        try:
            data = helper["chord_data_from_slope_numerator"](
                pencil["h"],
                pencil["Nx"],
                pencil["Ny"],
                M,
                pencil["A"],
                pencil["B"],
                pencil["Delta"],
                exact_context["ring"],
                exact_context["field"],
                expected_q_degree=4,
            )
        except ArithmeticError:
            continue
        q = data["q"]
        if (
            not q.is_irreducible()
            or q.gcd(pencil["Delta"]).degree()
            or q.gcd(pencil["h"]).degree()
        ):
            continue
        return {
            "label": "synthetic-known-norm8-bisection",
            "priority_rank": int(first_row["priority_rank"]),
            "orbit_mask": int(first_row["orbit_mask"]),
            "orbit_hex": first_row["orbit_hex"],
            "section_basis_w": list(map(int, word)),
            "construction_chart": pencil["chart"],
            "pencil_parameter_lambda": rational_text(parameter),
            "trace": {
                "h_coefficients_low_to_high": coefficient_text(pencil["h"]),
                "Nx_coefficients_low_to_high": coefficient_text(pencil["Nx"]),
                "Ny_coefficients_low_to_high": coefficient_text(pencil["Ny"]),
                "M0_coefficients_low_to_high": coefficient_text(pencil["M0"]),
            },
            "branch_quartic_coefficients_low_to_high": coefficient_text(q),
            "branch_quartic_irreducible_over_QQ": True,
            "branch_quartic_squarefree": True,
            "branch_quartic_coprime_to_surface_discriminant": True,
            "lifted_bisection": {
                "cover": "s^2=q(u)",
                "x0_coefficients_low_to_high": coefficient_text(data["x0"]),
                "x1_coefficients_low_to_high": coefficient_text(data["x1"]),
                "y0_coefficients_low_to_high": coefficient_text(data["y0"]),
                "y1_coefficients_low_to_high": coefficient_text(data["y1"]),
            },
        }
    raise ArithmeticError("failed to construct a squarefree synthetic control")


def selected_rank_one_targets(shortlist, ranks):
    exact_keys = {
        row["pair_key"]
        for row in ranks["results"]
        if row.get("status") == "completed"
        and int(row.get("rank_lower_bound", -1)) == 1
        and int(row.get("rank_upper_bound", -1)) == 1
    }
    by_key = {row["pair_key"]: row for row in shortlist["pairs"]}
    ordered_keys = [
        key for key in ranks["completed_pairs_ranked_by_upper_then_lower_bound"]
        if key in exact_keys
    ]
    if len(exact_keys) != 17 or set(ordered_keys) != exact_keys:
        raise ArithmeticError("expected exactly seventeen certified rank-one pair bases")
    targets = []
    for index, key in enumerate(ordered_keys):
        row = by_key[key]
        targets.append(
            {
                "index": index,
                "pair_key": key,
                "shortlist_rank": int(row["shortlist_rank"]),
                "labels": row["labels"],
                "product_quartic_coefficients_low_to_high": row[
                    "product_quartic_coefficients_low_to_high"
                ],
                "base_jacobian_rank": 1,
            }
        )
    return targets


def modular_context(prime, direct, helper, maximum_coefficient):
    coefficient_field = GF(prime)
    ring = PolynomialRing(coefficient_field, "u")
    field = ring.fraction_field()
    model = direct["weierstrass_model"]
    A = parse_polynomial(
        model["A_coefficients_low_to_high"], ring, coefficient_field
    )
    B = parse_polynomial(
        model["B_coefficients_low_to_high"], ring, coefficient_field
    )
    Delta = parse_polynomial(
        model["discriminant_coefficients_low_to_high"], ring, coefficient_field
    )
    if Delta != -16 * (4 * A**3 + 27 * B**2) or Delta.gcd(Delta.derivative()).degree():
        raise ArithmeticError(f"p={prime} is not a good alternate-Q80 model prime")
    curve = EllipticCurve(field, [A, B])
    basis = []
    for record in direct["sections"]["records"]:
        X = parse_rational_function(record["X"], ring, field, coefficient_field)
        Y = parse_rational_function(record["Y"], ring, field, coefficient_field)
        basis.append(curve(X, Y))
    context = {
        "prime": int(prime),
        "coefficient_field": coefficient_field,
        "ring": ring,
        "field": field,
        "A": A,
        "B": B,
        "Delta": Delta,
        "curve": curve,
        "basis": basis,
        "helper": helper,
    }
    context["multiples"] = [
        {
            coefficient: coefficient * point
            for coefficient in range(-maximum_coefficient, maximum_coefficient + 1)
        }
        for point in basis
    ]
    return context


def target_lookup(context, targets, control):
    ring = context["ring"]
    coefficient_field = context["coefficient_field"]
    rows = []
    for target in targets:
        rows.append(
            (
                target["pair_key"],
                parse_polynomial(
                    target["product_quartic_coefficients_low_to_high"],
                    ring,
                    coefficient_field,
                ),
            )
        )
    rows.append(
        (
            control["label"],
            parse_polynomial(
                control["branch_quartic_coefficients_low_to_high"],
                ring,
                coefficient_field,
            ),
        )
    )
    lookups = {}
    for chart in ("finite", "inverted_at_infinity"):
        lookup = {}
        for index, (label, polynomial) in enumerate(rows):
            if chart == "inverted_at_infinity":
                polynomial = reciprocal_polynomial(polynomial, 4, ring)
            key = projective_key(polynomial_vector(polynomial))
            if key is None:
                raise ArithmeticError(
                    f"target {label} vanishes identically modulo {context['prime']}"
                )
            lookup.setdefault(key, []).append(index)
        lookups[chart] = lookup
    return lookups


def modular_matches(context, word, lookups, target_count):
    trace = trace_from_word(context, word, context["multiples"])
    if trace.is_zero():
        raise ArithmeticError("norm-eight trace reduced to zero")
    pencil = trace_pencil_family(context, trace)
    family = pencil["q_lambda_family"]
    lookup = lookups[pencil["chart"]]
    field = context["coefficient_field"]
    matches = {index: [] for index in range(target_count + 1)}
    for parameter in field:
        q = family[4]
        for coefficient in reversed(family[:4]):
            q = q * parameter + coefficient
        key = projective_key(polynomial_vector(q))
        indices = range(target_count + 1) if key is None else lookup.get(key, ())
        for index in indices:
            matches[index].append(str(int(parameter)))
    # lambda=infinity is the leading homogeneous coefficient q_4=h^2.
    key = projective_key(polynomial_vector(family[4]))
    indices = range(target_count + 1) if key is None else lookup.get(key, ())
    for index in indices:
        matches[index].append("infinity")
    return pencil["chart"], {index: values for index, values in matches.items() if values}


def exact_parameters(family, target, ring):
    parameter_ring = PolynomialRing(QQ, "lambda")
    coefficient_polynomials = [
        parameter_ring([family[power][base_degree] for power in range(5)])
        for base_degree in range(5)
    ]
    target_vector = [QQ(target[index]) for index in range(5)]
    pivot = next(index for index, value in enumerate(target_vector) if value)
    equations = [
        target_vector[pivot] * coefficient_polynomials[index]
        - target_vector[index] * coefficient_polynomials[pivot]
        for index in range(5)
        if index != pivot
    ]
    nonzero = [equation for equation in equations if equation]
    if not nonzero:
        raise ArithmeticError("branch family is projectively constant")
    common = nonzero[0]
    for equation in nonzero[1:]:
        common = common.gcd(equation)
    roots = [root for root, _multiplicity in common.roots(QQ)] if common.degree() else []
    infinity_q = family[4]
    infinity_key = [infinity_q[index] for index in range(5)]
    infinity_matches = all(
        target_vector[pivot] * infinity_key[index]
        == target_vector[index] * infinity_key[pivot]
        for index in range(5)
    )
    return roots, infinity_matches, common


def squareclass_multiplier(polynomial, target):
    target_vector = polynomial_vector(target)
    polynomial_vector_value = polynomial_vector(polynomial)
    pivot = next(index for index, value in enumerate(target_vector) if value)
    multiplier = QQ(polynomial_vector_value[pivot] / target_vector[pivot])
    if any(
        polynomial_vector_value[index] != multiplier * target_vector[index]
        for index in range(5)
    ):
        raise ArithmeticError("candidate coefficient vectors are not proportional")
    return multiplier


def exact_resolve(row, target_indices, targets, exact_context):
    word = parse_vector(row["section_basis_w"])
    trace = trace_from_word(exact_context, word)
    pencil = trace_pencil_family(exact_context, trace)
    ring = exact_context["ring"]
    records = []
    hits = []
    for target_index in target_indices:
        target_record = targets[target_index]
        target = parse_polynomial(
            target_record["product_quartic_coefficients_low_to_high"], ring
        )
        if pencil["chart"] == "inverted_at_infinity":
            target = reciprocal_polynomial(target, 4, ring)
        roots, infinity_matches, common = exact_parameters(
            pencil["q_lambda_family"], target, ring
        )
        target_result = {
            "target_index": int(target_index),
            "pair_key": target_record["pair_key"],
            "coefficient_gcd_degree": int(common.degree()),
            "rational_finite_parameters": [rational_text(root) for root in roots],
            "infinity_matches": bool(infinity_matches),
            "squareclass_hits": [],
        }
        for parameter in roots:
            q = sum(
                (parameter**power * pencil["q_lambda_family"][power] for power in range(5)),
                ring.zero(),
            )
            multiplier = squareclass_multiplier(q, target)
            if multiplier.is_square():
                hit = {
                    "target_index": int(target_index),
                    "pair_key": target_record["pair_key"],
                    "parameter": rational_text(parameter),
                    "construction_chart": pencil["chart"],
                    "square_multiplier": rational_text(multiplier.sqrt()),
                    "branch_quartic_coefficients_low_to_high": coefficient_text(q),
                }
                target_result["squareclass_hits"].append(hit)
                hits.append(hit)
        if infinity_matches:
            multiplier = squareclass_multiplier(pencil["q_lambda_family"][4], target)
            if multiplier.is_square():
                hit = {
                    "target_index": int(target_index),
                    "pair_key": target_record["pair_key"],
                    "parameter": "infinity",
                    "construction_chart": pencil["chart"],
                    "square_multiplier": rational_text(multiplier.sqrt()),
                    "branch_quartic_coefficients_low_to_high": coefficient_text(
                        pencil["q_lambda_family"][4]
                    ),
                }
                target_result["squareclass_hits"].append(hit)
                hits.append(hit)
        records.append(target_result)
    return records, hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", default="131,137")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ledger-output", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    primes = [int(value) for value in arguments.primes.split(",") if value]
    if not primes or arguments.start < 0 or arguments.limit is not None and arguments.limit <= 0:
        parser.error("require primes, nonnegative --start, and positive --limit")

    direct = json.loads(DIRECT.read_text())
    priority_certificate = json.loads(PRIORITY_CERTIFICATE.read_text())
    shortlist = json.loads(SHORTLIST.read_text())
    ranks = json.loads(RANKS.read_text())
    if priority_certificate["status"] != "PASS_EXACT_COMPLETE_ALTERNATE_NORM8_PENCIL_PRIORITY":
        raise ArithmeticError("norm-eight priority certificate is not complete")
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("alternate-Q80 section basis is not saturated")
    with PRIORITY_TABLE.open(newline="") as stream:
        all_rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(all_rows) != 63917:
        raise ArithmeticError("norm-eight priority table is not complete")
    stop = len(all_rows) if arguments.limit is None else min(
        len(all_rows), arguments.start + arguments.limit
    )
    rows = all_rows[arguments.start:stop]
    if not rows:
        raise ArithmeticError("selected norm-eight interval is empty")

    helper = runpy.run_path(str(HELPER))
    exact_context = exact_model_and_basis(direct, helper)
    control = build_positive_control(exact_context, all_rows[0])
    targets = selected_rank_one_targets(shortlist, ranks)
    target_count = len(targets)
    maximum_coefficient = max(
        abs(int(value))
        for row in rows
        for value in row["section_basis_w"].split()
    )
    modular_contexts = {
        prime: modular_context(prime, direct, helper, maximum_coefficient)
        for prime in primes
    }
    lookups = {
        prime: target_lookup(context, targets, control)
        for prime, context in modular_contexts.items()
    }

    ledger_rows = []
    unresolved = []
    bad_reductions = Counter()
    prime_obstructions = Counter()
    control_hits = []
    for offset, row in enumerate(rows):
        word = parse_vector(row["section_basis_w"])
        survivor_indices = set(range(target_count))
        residues = {}
        good_primes = []
        for prime in primes:
            try:
                chart, matches = modular_matches(
                    modular_contexts[prime], word, lookups[prime], target_count
                )
            except (ArithmeticError, ValueError, ZeroDivisionError) as error:
                bad_reductions[(prime, type(error).__name__)] += 1
                continue
            good_primes.append(prime)
            present = set(matches) & set(range(target_count))
            survivor_indices &= present
            residues[str(prime)] = {
                "chart": chart,
                "target_residues": {
                    str(index): matches[index]
                    for index in sorted(present & survivor_indices)
                },
            }
            if (
                int(row["orbit_mask"]) == int(control["orbit_mask"])
                and target_count in matches
            ):
                expected_parameter = control["pencil_parameter_lambda"]
                expected_residue = str(
                    int(modular_contexts[prime]["coefficient_field"](QQ(expected_parameter)))
                )
                if expected_residue not in matches[target_count]:
                    raise ArithmeticError("synthetic control recovered at wrong parameter")
                control_hits.append(
                    {
                        "prime": prime,
                        "orbit_mask": int(row["orbit_mask"]),
                        "parameter_residue": expected_residue,
                    }
                )
            if not survivor_indices:
                prime_obstructions[prime] += 1
                break
        result = {
            "priority_rank": int(row["priority_rank"]),
            "orbit_mask": int(row["orbit_mask"]),
            "orbit_hex": row["orbit_hex"],
            "good_primes": good_primes,
            "surviving_target_indices": sorted(survivor_indices),
            "modular_residues": residues if survivor_indices else {},
        }
        ledger_rows.append(result)
        if survivor_indices:
            unresolved.append((row, sorted(survivor_indices), result))
        if arguments.progress_every and (offset + 1) % arguments.progress_every == 0:
            print(
                f"BISECTIONINVERSIONPROGRESS|done={offset + 1}|total={len(rows)}"
                f"|unresolved={len(unresolved)}",
                flush=True,
            )

    exact_resolutions = []
    hits = []
    for row, indices, modular_result in unresolved:
        resolutions, row_hits = exact_resolve(row, indices, targets, exact_context)
        exact_resolutions.append(
            {
                "priority_rank": int(row["priority_rank"]),
                "orbit_mask": int(row["orbit_mask"]),
                "orbit_hex": row["orbit_hex"],
                "modular_surviving_target_indices": indices,
                "targets": resolutions,
            }
        )
        for hit in row_hits:
            hits.append(
                {
                    "priority_rank": int(row["priority_rank"]),
                    "orbit_mask": int(row["orbit_mask"]),
                    "orbit_hex": row["orbit_hex"],
                    "section_basis_w": list(map(int, parse_vector(row["section_basis_w"]))),
                    **hit,
                }
            )

    if arguments.start == 0 and stop == len(all_rows):
        expected_control_primes = {
            prime
            for prime in primes
            if any(item["prime"] == prime for item in control_hits)
        }
        if not expected_control_primes:
            raise ArithmeticError("synthetic positive control was not recovered")
        coverage_status = "PASS_EXACT_COMPLETE_NORM8_BISECTION_INVERSION"
    else:
        coverage_status = "PASS_EXACT_BOUNDED_NORM8_BISECTION_INVERSION_INTERVAL"

    ledger_fields = [
        "priority_rank",
        "orbit_mask",
        "orbit_hex",
        "good_primes",
        "surviving_target_indices",
    ]
    ledger_lines = ["\t".join(ledger_fields)]
    for row in ledger_rows:
        ledger_lines.append(
            "\t".join(
                [
                    str(row["priority_rank"]),
                    str(row["orbit_mask"]),
                    row["orbit_hex"],
                    ",".join(map(str, row["good_primes"])),
                    ",".join(map(str, row["surviving_target_indices"])),
                ]
            )
        )
    ledger_text = "\n".join(ledger_lines) + "\n"

    payload = {
        "schema": "elkies-k3.r17-norm12-11952-product-bisection-inversion.v1",
        "status": coverage_status,
        "dictionary": {
            "forward_map": (
                "A QQ-rational genus-one bisection C with branch d pulls back to "
                "sections P,sigma(P); T=P-sigma(P) is anti-invariant and is a section "
                "of E^(d). Section translation changes P and sigma(P) together and "
                "therefore leaves T and d unchanged."
            ),
            "inverse_map": (
                "An anti-invariant twist section T comes from such a bisection exactly "
                "when T lies in (1-sigma)E(QQ(u)(sqrt(d))); equivalently there are "
                "P and tau in E(QQ(u)) with 2P=T+tau. The translation class is the "
                "class of tau modulo 2E(QQ(u))."
            ),
            "minimal_layer": (
                "For a rootless K3, a section-nonnegative degree-two isotropic class "
                "disjoint from O is D_tau=O+tau=(2,2,w), with height(tau)=w.M.w=8. "
                "Thus all minimal bisection coboundaries occur in the enumerated "
                "minimum-norm-eight parity cosets."
            ),
            "integral_warning": (
                "Surjectivity after tensoring with QQ does not imply integral "
                "surjectivity: the cokernel of 1-sigma is 2-primary. A negative "
                "bisection inversion therefore excludes minimal product-character "
                "sections only after a separate saturation/descent proof that the "
                "relevant minimal twist section is a coboundary."
            ),
        },
        "positive_control": {
            **control,
            "recovered_modular_hits": control_hits,
            "status": "PASS_SYNTHETIC_KNOWN_BISECTION_RECOVERED",
        },
        "targets": targets,
        "search": {
            "primes": primes,
            "interval": {"start": arguments.start, "stop": stop},
            "complete_class_count": len(all_rows),
            "searched_class_count": len(rows),
            "coefficient_comparison": (
                "projective equality of the five old-base quartic coefficients over "
                "every lambda in P1(F_p), followed by exact QQ coefficient gcds"
            ),
            "lambda_infinity_included": True,
            "modular_zero_vector_treated_as_survivor": True,
            "bad_reduction_histogram": {
                f"p{prime}:{name}": count
                for (prime, name), count in sorted(bad_reductions.items())
            },
            "first_obstructing_prime_histogram": {
                str(prime): count for prime, count in sorted(prime_obstructions.items())
            },
            "modular_survivor_class_count": len(unresolved),
            "exact_resolution_count": len(exact_resolutions),
            "squareclass_hit_count": len(hits),
        },
        "ledger": relative(arguments.ledger_output),
        "ledger_sha256": hashlib.sha256(ledger_text.encode()).hexdigest(),
        "exact_resolutions": exact_resolutions,
        "hits": hits,
        "conclusion": (
            (
                "At least one selected product character is realized by a member of "
                "the complete norm-eight bisection-pencil layer."
                if hits
                else "No selected product character is realized by a member of the "
                "searched norm-eight bisection-pencil layer."
            )
        ),
        "claim_boundary": {
            "exact_if_complete": (
                "When the full 63,917-class interval is present, the bisection-image "
                "negative is exhaustive for the norm-eight/pole-order-zero layer."
            ),
            "not_implied_without_integral_descent": (
                "A no-hit result alone does not prove that the product-twist MW group "
                "has no disjoint height-eight section; a non-coboundary class in the "
                "2-primary Tate-cohomology quotient is not represented by a bisection pencil."
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                HELPER,
                DIRECT,
                PRIORITY_CERTIFICATE,
                PRIORITY_TABLE,
                SHORTLIST,
                RANKS,
            )
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact elliptic-curve group law over rational function fields",
                "finite-field polynomial arithmetic",
                "exact univariate gcd and rational roots",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "search_r17_norm12_11952_product_bisection_inversion.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.ledger_output.exists() or arguments.ledger_output.read_text() != ledger_text:
            raise ArithmeticError("stored bisection-inversion ledger differs from replay")
        if not arguments.output.exists() or arguments.output.read_text() != serialized:
            raise ArithmeticError("stored bisection-inversion certificate differs from replay")
    else:
        arguments.ledger_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.ledger_output.write_text(ledger_text)
        arguments.output.write_text(serialized)
    print(
        f"BISECTIONINVERSION|classes={len(rows)}|targets={target_count}"
        f"|modular_survivors={len(unresolved)}|hits={len(hits)}"
        f"|status={coverage_status}|output={relative(arguments.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
