#!/usr/bin/env sage-python
"""Invert the 17 product quartics through the 49 deep trace parities.

The completed norm-eight inversion covers the disjoint genus-one bisection
carriers.  For any pulled-back height-16 anti-invariant point with zero Tate
class, integral character glue leaves one further possibility: an invariant
trace parity of minimum norm twelve.  The alternate-Q80 Mordell--Weil lattice
has exactly 49 such deep parities.

For a norm-twelve trace ``tau=(Nx/h^2,Ny/h^3)``, ``deg(h)=4``.  Every regular
height-ten half-point carrier has chord slope

    M = M0 + lambda*h^2.

The chord discriminant ``q_lambda`` has degree at most eight.  Its normalized
double cover has target quartic character ``d`` exactly when

    q_lambda = d*r^2

for a polynomial ``r`` of degree at most two, including the possible constant
square factor.  This script tests every lambda in P1(F_p) and rejects a
trace/target pair only at a good prime with no such factorization.  Survivors
are resolved over QQ.  A quartic produced by the first deep trace at lambda=0
is recovered through the same modular and exact path as a positive control.

Together with the norm-eight inversion, a no-hit result excludes the zero
Tate class for any height-eight product-twist section satisfying the stated
rootless/direct local-height hypotheses.  The existence of such a section
and the size of the rest of the anti-invariant lattice remain open.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import runpy

from sage.all import GF, PolynomialRing, QQ, ZZ, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
NORM8_INVERSION_SCRIPT = (
    ROOT / "elkies-k3/scripts/search_r17_norm12_11952_product_bisection_inversion.sage"
)
DIRECT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
PARITY = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json"
)
NORM8_INVERSION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-bisection-inversion-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-deep-trace-inversion-v1.json"
)
DEFAULT_LEDGER = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-deep-trace-inversion-v1.tsv"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def coefficient_text(polynomial, length=None) -> list[str]:
    if length is None:
        length = max(1, polynomial.degree() + 1)
    return [rational_text(polynomial[index]) for index in range(length)]


def deep_trace_family(context, word, inverse, helper):
    """Return the degree-at-most-eight chord family for one deep trace."""

    trace = inverse["trace_from_word"](context, word, context.get("multiples"))
    if trace.is_zero():
        raise ArithmeticError("deep trace reduced to zero")
    ring = context["ring"]
    field = context["field"]
    X, Y = field(trace[0]), field(trace[1])
    A, B, Delta = context["A"], context["B"], context["Delta"]
    frame = helper["trace_chord_frame"](X, Y, ring)
    chart = "finite"
    if frame["h"].degree() != 4:
        chart = "inverted_at_infinity"
        A = inverse["reciprocal_polynomial"](A, 8, ring)
        B = inverse["reciprocal_polynomial"](B, 12, ring)
        Delta = inverse["reciprocal_polynomial"](Delta, 24, ring)
        X = helper["invert_rational"](X, 4, ring, field)
        Y = helper["invert_rational"](Y, 6, ring, field)
        frame = helper["trace_chord_frame"](X, Y, ring)
    if frame["h"].degree() != 4:
        raise ArithmeticError("deep trace does not have four finite poles in either chart")

    h = frame["h"]
    Nx = frame["Nx"]
    Ny = frame["Ny"]
    M0 = frame["M0"]
    h2 = h**2
    h6 = h**6
    numerators = [
        M0**4 - 6 * M0**2 * Nx - 8 * M0 * Ny - 3 * Nx**2 - 4 * A * h**4,
        4 * M0**3 * h2 - 12 * M0 * h2 * Nx - 8 * h2 * Ny,
        6 * M0**2 * h2**2 - 6 * h2**2 * Nx,
        4 * M0 * h6,
        h2 * h6,
    ]
    family = []
    for numerator in numerators:
        quotient, remainder = ring(numerator).quo_rem(h6)
        if remainder:
            raise ArithmeticError("deep chord coefficient is not divisible by h^6")
        quotient = ring(quotient)
        if quotient.degree() > 8:
            raise ArithmeticError("deep chord coefficient exceeds degree eight")
        family.append(quotient)
    if family[3] != 4 * M0 or family[4] != h2:
        raise ArithmeticError("deep chord leading coefficients changed")
    return {
        "chart": chart,
        "A": A,
        "B": B,
        "Delta": Delta,
        **frame,
        "q_lambda_family": family,
    }


def evaluated_family(family, parameter, ring):
    if parameter == "infinity":
        return ring(family[4])
    value = family[4]
    for coefficient in reversed(family[:4]):
        value = value * parameter + coefficient
    return ring(value)


def squarefactor_match(polynomial, target):
    """Test q=target*r^2, preserving the constant squareclass."""

    if not polynomial:
        return None
    quotient, remainder = polynomial.quo_rem(target)
    if remainder or not quotient.is_square():
        return False
    return quotient.sqrt()


def build_positive_control(exact_context, deep_record, inverse, helper):
    word = vector(ZZ, deep_record["section_basis_w"])
    pencil = deep_trace_family(exact_context, word, inverse, helper)
    for parameter in map(QQ, (0, 1, -1, 2, -2, 3, -3)):
        q = evaluated_family(
            pencil["q_lambda_family"], parameter, exact_context["ring"]
        )
        if (
            q.degree() != 4
            or not q.is_squarefree()
            or not q.is_irreducible()
            or q.gcd(pencil["Delta"]).degree()
        ):
            continue
        M = pencil["M0"] + parameter * pencil["h"]**2
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
        if data["q"] != q or not data["branch_fibres_smooth"]:
            raise ArithmeticError("deep positive-control chord verification failed")
        return {
            "label": "synthetic-known-deep-trace-quartic",
            "orbit_mask": int(deep_record["orbit_mask"]),
            "orbit_hex": deep_record["orbit_hex"],
            "section_basis_w": list(map(int, word)),
            "construction_chart": pencil["chart"],
            "pencil_parameter_lambda": rational_text(parameter),
            "branch_quartic_coefficients_low_to_high": coefficient_text(q, 5),
            "branch_quartic_irreducible_over_QQ": True,
            "branch_quartic_squarefree": True,
            "branch_quartic_coprime_to_surface_discriminant": True,
            "exact_lifted_bisection_coefficients": {
                key: coefficient_text(data[key]) for key in ("x0", "x1", "y0", "y1")
            },
        }
    raise ArithmeticError("failed to construct the deep-trace positive control")


def target_polynomial(context, target, chart, inverse):
    polynomial = inverse["parse_polynomial"](
        target["product_quartic_coefficients_low_to_high"],
        context["ring"],
        context.get("coefficient_field", QQ),
    )
    if chart == "inverted_at_infinity":
        polynomial = inverse["reciprocal_polynomial"](
            polynomial, 4, context["ring"]
        )
    return polynomial


def modular_matches(context, word, targets, control, inverse, helper):
    pencil = deep_trace_family(context, word, inverse, helper)
    polynomials = [
        target_polynomial(context, target, pencil["chart"], inverse)
        for target in targets
    ]
    polynomials.append(
        target_polynomial(
            context,
            {
                "product_quartic_coefficients_low_to_high": control[
                    "branch_quartic_coefficients_low_to_high"
                ]
            },
            pencil["chart"],
            inverse,
        )
    )
    matches = {index: [] for index in range(len(polynomials))}
    parameters = list(context["coefficient_field"]) + ["infinity"]
    for parameter in parameters:
        q = evaluated_family(
            pencil["q_lambda_family"], parameter, context["ring"]
        )
        for index, target in enumerate(polynomials):
            root = squarefactor_match(q, target)
            if root is None or root is not False:
                matches[index].append(
                    "infinity" if parameter == "infinity" else str(int(parameter))
                )
    return pencil["chart"], {
        index: residues for index, residues in matches.items() if residues
    }


def exact_parameters(pencil, target, context):
    """Return all rational lambda for which target divides q_lambda."""

    parameter_ring = PolynomialRing(QQ, "lambda")
    parameter = parameter_ring.gen()
    base_ring = PolynomialRing(parameter_ring, "u")
    q = sum(
        (
            base_ring(pencil["q_lambda_family"][power]) * parameter**power
            for power in range(5)
        ),
        base_ring.zero(),
    )
    divisor = base_ring(target)
    quotient, remainder = q.quo_rem(divisor)
    equations = [parameter_ring(remainder[index]) for index in range(4)]
    nonzero = [equation for equation in equations if equation]
    if not nonzero:
        return None, quotient
    common = nonzero[0]
    for equation in nonzero[1:]:
        common = common.gcd(equation)
    return [root for root, unused in common.roots(QQ)], quotient


def exact_resolve(deep_record, indices, targets, exact_context, inverse, helper):
    word = vector(ZZ, deep_record["section_basis_w"])
    pencil = deep_trace_family(exact_context, word, inverse, helper)
    records = []
    hits = []
    for index in indices:
        target_record = targets[index]
        target = target_polynomial(
            exact_context, target_record, pencil["chart"], inverse
        )
        parameters, unused_symbolic_quotient = exact_parameters(
            pencil, target, exact_context
        )
        record = {
            "target_index": int(index),
            "pair_key": target_record["pair_key"],
            "divisibility_identically_in_lambda": parameters is None,
            "rational_divisibility_parameters": (
                [] if parameters is None else [rational_text(value) for value in parameters]
            ),
            "squareclass_hits": [],
        }
        # An identically divisible family requires a separate symbolic square
        # analysis.  Retain it unresolved rather than claiming a negative.
        if parameters is not None:
            for parameter in parameters:
                q = evaluated_family(
                    pencil["q_lambda_family"], parameter, exact_context["ring"]
                )
                square_root = squarefactor_match(q, target)
                if square_root is None:
                    record.setdefault("zero_polynomial_parameters", []).append(
                        rational_text(parameter)
                    )
                elif square_root is not False:
                    hit = {
                        "target_index": int(index),
                        "pair_key": target_record["pair_key"],
                        "parameter": rational_text(parameter),
                        "construction_chart": pencil["chart"],
                        "square_factor_coefficients_low_to_high": coefficient_text(
                            square_root
                        ),
                    }
                    record["squareclass_hits"].append(hit)
                    hits.append(hit)
        infinity_q = pencil["q_lambda_family"][4]
        infinity_root = squarefactor_match(infinity_q, target)
        if infinity_root is None:
            record["zero_polynomial_at_infinity"] = True
        elif infinity_root is not False:
            hit = {
                "target_index": int(index),
                "pair_key": target_record["pair_key"],
                "parameter": "infinity",
                "construction_chart": pencil["chart"],
                "square_factor_coefficients_low_to_high": coefficient_text(
                    infinity_root
                ),
            }
            record["squareclass_hits"].append(hit)
            hits.append(hit)
        records.append(record)
    return records, hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", default="131,137")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ledger-output", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    primes = [int(value) for value in arguments.primes.split(",") if value]
    if not primes or any(prime < 5 or not ZZ(prime).is_prime() for prime in primes):
        parser.error("--primes must be a nonempty list of odd primes at least five")

    direct = json.loads(DIRECT.read_text())
    parity = json.loads(PARITY.read_text())
    norm8_inversion = json.loads(NORM8_INVERSION.read_text())
    if parity.get("status") != "PASS_EXACT_PRODUCT_TATE_PARITY_REDUCTION":
        raise ArithmeticError("Tate-parity certificate is not exact")
    if parity.get("mathematical_status") != "UNKNOWN_PRODUCT_TWIST_TATE_COHOMOLOGY":
        raise ArithmeticError("input parity certificate crossed its proof boundary")
    if norm8_inversion.get("status") != "PASS_EXACT_COMPLETE_NORM8_BISECTION_INVERSION":
        raise ArithmeticError("norm-eight inversion is not complete")
    if norm8_inversion["search"]["squareclass_hit_count"] != 0:
        raise ArithmeticError("norm-eight inversion now has a product-character hit")
    deep_records = parity["invariant_trace_parity"]["deep_norm12_classes"]
    targets = norm8_inversion["targets"]
    if len(deep_records) != 49 or len(targets) != 17:
        raise ArithmeticError("expected exactly 49 deep parities and 17 targets")

    inverse = runpy.run_path(str(NORM8_INVERSION_SCRIPT))
    helper = runpy.run_path(str(HELPER))
    exact_context = inverse["exact_model_and_basis"](direct, helper)
    control = build_positive_control(exact_context, deep_records[0], inverse, helper)
    maximum_coefficient = max(
        abs(int(entry))
        for record in deep_records
        for entry in record["section_basis_w"]
    )
    modular_contexts = {
        prime: inverse["modular_context"](
            prime, direct, helper, maximum_coefficient
        )
        for prime in primes
    }

    ledger_rows = []
    unresolved = []
    bad_reductions = Counter()
    first_obstructions = Counter()
    control_hits = []
    for deep_record in deep_records:
        word = vector(ZZ, deep_record["section_basis_w"])
        survivors = set(range(len(targets)))
        good_primes = []
        modular_residues = {}
        for prime in primes:
            try:
                chart, matches = modular_matches(
                    modular_contexts[prime],
                    word,
                    targets,
                    control,
                    inverse,
                    helper,
                )
            except (ArithmeticError, ValueError, ZeroDivisionError) as error:
                bad_reductions[(prime, type(error).__name__)] += 1
                continue
            good_primes.append(prime)
            present = set(matches) & set(range(len(targets)))
            survivors &= present
            modular_residues[str(prime)] = {
                "chart": chart,
                "target_residues": {
                    str(index): matches[index]
                    for index in sorted(survivors)
                    if index in matches
                },
            }
            if int(deep_record["orbit_mask"]) == int(control["orbit_mask"]):
                control_index = len(targets)
                expected = control["pencil_parameter_lambda"]
                expected_residue = str(
                    int(modular_contexts[prime]["coefficient_field"](QQ(expected)))
                )
                if expected_residue not in matches.get(control_index, []):
                    raise ArithmeticError("deep positive control was not recovered")
                control_hits.append(
                    {
                        "prime": prime,
                        "orbit_mask": int(deep_record["orbit_mask"]),
                        "parameter_residue": expected_residue,
                    }
                )
            if not survivors:
                first_obstructions[prime] += 1
                break
        row = {
            "orbit_mask": int(deep_record["orbit_mask"]),
            "orbit_hex": deep_record["orbit_hex"],
            "minimum_norm": 12,
            "good_primes": good_primes,
            "surviving_target_indices": sorted(survivors),
            "modular_residues": modular_residues if survivors else {},
        }
        ledger_rows.append(row)
        if survivors:
            unresolved.append((deep_record, sorted(survivors)))

    exact_resolutions = []
    hits = []
    unresolved_symbolic_families = []
    for deep_record, indices in unresolved:
        records, row_hits = exact_resolve(
            deep_record, indices, targets, exact_context, inverse, helper
        )
        exact_resolutions.append(
            {
                "orbit_mask": int(deep_record["orbit_mask"]),
                "orbit_hex": deep_record["orbit_hex"],
                "targets": records,
            }
        )
        for record in records:
            if record["divisibility_identically_in_lambda"]:
                unresolved_symbolic_families.append(
                    {
                        "orbit_mask": int(deep_record["orbit_mask"]),
                        "target_index": record["target_index"],
                        "pair_key": record["pair_key"],
                    }
                )
        for hit in row_hits:
            hits.append(
                {
                    "orbit_mask": int(deep_record["orbit_mask"]),
                    "orbit_hex": deep_record["orbit_hex"],
                    "section_basis_w": deep_record["section_basis_w"],
                    **hit,
                }
            )

    if not control_hits:
        raise ArithmeticError("no good prime recovered the deep positive control")
    complete = not unresolved_symbolic_families
    status = (
        "PASS_EXACT_COMPLETE_DEEP_TRACE_PRODUCT_INVERSION"
        if complete
        else "INCOMPLETE_DEEP_TRACE_PRODUCT_INVERSION"
    )

    fields = [
        "orbit_mask",
        "orbit_hex",
        "minimum_norm",
        "good_primes",
        "surviving_target_indices",
    ]
    lines = ["\t".join(fields)]
    for row in ledger_rows:
        lines.append(
            "\t".join(
                (
                    str(row["orbit_mask"]),
                    row["orbit_hex"],
                    str(row["minimum_norm"]),
                    ",".join(map(str, row["good_primes"])),
                    ",".join(map(str, row["surviving_target_indices"])),
                )
            )
        )
    ledger_text = "\n".join(lines) + "\n"

    payload = {
        "schema": "elkies-k3.r17-norm12-11952-product-deep-trace-inversion.v1",
        "status": status,
        "mathematical_status": "UNKNOWN_PRODUCT_TWIST_MW_GROUPS",
        "dictionary": {
            "residual_trace_layer": (
                "The 49 invariant parity classes of minimum norm twelve are the "
                "only possible zero-class carriers left for any height-eight T after "
                "the rootless height bound and complete norm-eight inversion."
            ),
            "regular_slope": "M=M0+lambda*h^2 with deg(h)=4",
            "branch_squareclass_test": (
                "q_lambda has the target quartic squareclass iff "
                "q_lambda=target*r^2 for a polynomial r of degree at most two."
            ),
            "height": (
                "For pulled heights tau^2=24 and T^2=16, the half point "
                "R=(tau+T)/2 has height 10 and R.O=1."
            ),
        },
        "positive_control": {
            **control,
            "recovered_modular_hits": control_hits,
            "status": "PASS_SYNTHETIC_DEEP_TRACE_QUARTIC_RECOVERED",
        },
        "targets": targets,
        "search": {
            "primes": primes,
            "deep_trace_class_count": len(deep_records),
            "product_target_count": len(targets),
            "trace_target_comparisons": len(deep_records) * len(targets),
            "lambda_space": "every point of P1(F_p), including infinity",
            "zero_modular_q_treated_as_survivor": True,
            "constant_squareclass_preserved": True,
            "bad_reduction_histogram": {
                f"p{prime}:{name}": count
                for (prime, name), count in sorted(bad_reductions.items())
            },
            "first_obstructing_prime_histogram": {
                str(prime): count for prime, count in sorted(first_obstructions.items())
            },
            "modular_surviving_trace_count": len(unresolved),
            "modular_surviving_trace_target_pairs": sum(
                len(indices) for unused, indices in unresolved
            ),
            "exact_resolution_count": len(exact_resolutions),
            "unresolved_symbolic_family_count": len(unresolved_symbolic_families),
            "squareclass_hit_count": len(hits),
        },
        "ledger": relative(arguments.ledger_output),
        "ledger_sha256": sha256(ledger_text.encode()).hexdigest(),
        "exact_resolutions": exact_resolutions,
        "unresolved_symbolic_families": unresolved_symbolic_families,
        "hits": hits,
        "conclusion": (
            "No selected product quartic occurs in any norm-eight or deep norm-twelve "
            "trace carrier. Therefore every height-eight section T on a selected "
            "product twist, if one exists and satisfies the stated direct/local "
            "height gates, has nonzero class in Hhat^-1. This conclusion does not "
            "assume anti-invariant rank one."
            if complete and not hits
            else "The deep trace inversion did not close the zero-class carrier cases."
        ),
        "conditional_rank_one_corollary": (
            "If in addition the anti-invariant lattice is primitive rank one <16> "
            "generated by T, then Hhat^-1 is Z/2 and [T] is its unique nonzero class."
            if complete and not hits
            else "No rank-one corollary is available before the carrier cases close."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (
                HELPER,
                NORM8_INVERSION_SCRIPT,
                DIRECT,
                PARITY,
                NORM8_INVERSION,
            )
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact elliptic-curve group law over rational function fields",
                "finite-field polynomial square testing",
                "exact polynomial division and rational roots",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "search_r17_norm12_11952_product_deep_trace_inversion.sage --check"
        ),
        "proof_boundary": {
            "proved": (
                "All 49x17 deep trace/target cases are covered. A no-hit "
                "result, combined with the complete norm-eight inversion and the "
                "height/covering-spectrum identity, excludes the zero Tate class for "
                "any height-eight product-twist section under the stated direct/local "
                "height hypotheses."
            ),
            "not_proved": (
                "The existence, rank, saturation, or full lattice of any product-twist "
                "anti-invariant Mordell--Weil group is not computed. The full quotient "
                "may have other nonzero classes, especially in anti-rank above one. No "
                "height-eight product section or characteristic-zero rank-20 direction "
                "is found."
            ),
        },
    }

    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    ledger_output = arguments.ledger_output.resolve()
    if arguments.check:
        if not ledger_output.exists() or ledger_output.read_text() != ledger_text:
            raise SystemExit(f"stale or missing ledger: {ledger_output}")
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
    else:
        ledger_output.parent.mkdir(parents=True, exist_ok=True)
        ledger_output.write_text(ledger_text)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "PRODUCTDEEPTRACE|"
        f"traces={len(deep_records)}|targets={len(targets)}|"
        f"modular_survivor_pairs={payload['search']['modular_surviving_trace_target_pairs']}|"
        f"hits={len(hits)}|status={status}|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
