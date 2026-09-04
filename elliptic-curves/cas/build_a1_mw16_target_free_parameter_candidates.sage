#!/usr/bin/env sage-python
"""Sample and specialize anonymous A1/MW16 parameters without record targets.

The only family input is the sanitized five-fibration/nine-presentation
template.  A fixed height box and fixed local-prime blocks select parameters;
the exact saturated MW16 basis is then specialized on each selected fibre.
No known rank-jump parameter, curve model, point, rank, or local control is
loaded or evaluated.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
import sys
from time import perf_counter

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, vector


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/a1_mw16_family_template_v1.json"
MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
COMMON = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/a1_mw16_target_free_parameter_candidates_h300_v1.json"
PRIME_BLOCKS = (
    (19, 41, 43, 61, 71, 73, 79, 83),
    (89, 107, 113, 127, 131, 137, 139, 151),
    (157, 163, 167, 173, 179, 181, 191, 193, 197),
)


def load_source(name: str, path: Path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_sha256(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def qtext(value) -> str:
    value = QQ(value)
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def point_record(point):
    if point.is_zero():
        return {"infinity": True}
    return {"x": qtext(point[0]), "y": qtext(point[1])}


def polynomial_from_record(record, ring):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return ring.fraction_field()(numerator) / ring.fraction_field()(denominator)


def binary_quartic_invariants(quartic, coefficient_ring):
    e, d, c, b, a = [coefficient_ring(quartic[index]) for index in range(5)]
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    return invariant_i, invariant_j


def child_geometry(trace, old_a, old_b, old_ring, chord):
    frame = chord.trace_chord_frame(trace[0], trace[1], old_ring)
    h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("anonymous presentation left the finite-pole chart")
    parameter_ring = PolynomialRing(QQ, "lambda")
    parameter_variable = parameter_ring.gen()
    bivariate_ring = PolynomialRing(parameter_ring, "t")
    hh, nnx, nny, mm0 = map(bivariate_ring, (h, nx, ny, m0))
    slope_numerator = mm0 + parameter_variable * hh**2
    numerator = (
        slope_numerator**4
        - 6 * slope_numerator**2 * nnx
        - 8 * slope_numerator * nny
        - 3 * nnx**2
        - 4 * bivariate_ring(old_a) * hh**4
    )
    quartic, remainder = numerator.quo_rem(hh**6)
    if remainder or quartic.degree() != 4:
        raise ArithmeticError("residual chord did not produce a binary quartic")
    invariant_i, invariant_j = binary_quartic_invariants(quartic, parameter_ring)
    return h, nx, m0, quartic, -27 * invariant_i, -27 * invariant_j


def invert_mobius(function, value, ring):
    numerator = ring(function.numerator())
    denominator = ring(function.denominator())
    if numerator.degree() > 1 or denominator.degree() > 1:
        raise ArithmeticError("stored old-section base map is not Mobius")
    n0, n1 = numerator[0], numerator[1]
    d0, d1 = denominator[0], denominator[1]
    bottom = value * d1 - n1
    if not bottom:
        raise ZeroDivisionError("selected section meets quartic infinity")
    answer = (n0 - value * d0) / bottom
    if function(answer) != value:
        raise ArithmeticError("stored Mobius inversion failed")
    return QQ(answer)


def specialize(
    *, parameter, source_points, base_maps, old_ring, h, nx, m0,
    quartic, child_a, child_b
):
    specialized_a = QQ(child_a(parameter))
    specialized_b = QQ(child_b(parameter))
    if 4 * specialized_a**3 + 27 * specialized_b**2 == 0:
        return None, "SINGULAR_FIBRE"
    fixed_m = m0 + parameter * h**2
    fixed_quartic = old_ring([QQ(quartic[index](parameter)) for index in range(5)])
    sum_x = old_ring((fixed_m**2 - nx) // h**2)
    quartic_points = []
    try:
        for source_point, base_map in zip(source_points, base_maps):
            old_parameter = invert_mobius(base_map, parameter, old_ring)
            x_value = QQ(source_point[0](old_parameter))
            y_value = QQ(source_point[1](old_parameter))
            h_value = QQ(h(old_parameter))
            if not h_value:
                return None, "SELECTED_SECTION_MEETS_QUARTIC_INFINITY"
            w_value = (2 * x_value - QQ(sum_x(old_parameter))) / h_value
            if w_value**2 != fixed_quartic(old_parameter):
                raise ArithmeticError("selected old section missed the sampled quartic")
            quartic_points.append((old_parameter, w_value))
    except ZeroDivisionError:
        return None, "SELECTED_SECTION_SPECIALIZATION_UNDEFINED"

    t0, w0 = quartic_points[0]
    if not w0:
        return None, "POINTED_QUARTIC_ORIGIN_IS_BRANCH_POINT"
    shift_ring = PolynomialRing(QQ, "z")
    z = shift_ring.gen()
    shifted = shift_ring(fixed_quartic(t0 + z))
    ee, dd, cc, bb, aa = [QQ(shifted[index]) for index in range(5)]
    if ee != w0**2:
        raise ArithmeticError("pointed quartic constant term changed")
    a1g = dd / w0
    a2g = cc - dd**2 / (4 * w0**2)
    a3g = 2 * w0 * bb
    a4g = -4 * w0**2 * aa
    a6g = a2g * a4g
    b2g = a1g**2 + 4 * a2g
    b4g = a1g * a3g + 2 * a4g
    b6g = a3g**2 + 4 * a6g
    c4g = b2g**2 - 24 * b4g
    c6g = -b2g**3 + 36 * b2g * b4g - 216 * b6g
    if 81 * (-c4g / 48) != specialized_a or 729 * (-c6g / 864) != specialized_b:
        raise ArithmeticError("pointed quartic normalization missed the raw short model")

    raw_curve = EllipticCurve(QQ, [specialized_a, specialized_b])
    raw_points = []
    for old_parameter, w_value in quartic_points[1:]:
        zz = old_parameter - t0
        if not zz:
            return None, "GENERIC_SECTION_COLLIDES_WITH_SELECTED_ZERO"
        x_general = (2 * w0 * (w_value + w0) + dd * zz) / zz**2
        y_general = (
            4 * w0**2 * (w_value + w0)
            + 2 * w0 * dd * zz
            + (2 * w0 * cc - dd**2 / (2 * w0)) * zz**2
        ) / zz**3
        raw_points.append(
            raw_curve(
                9 * (x_general + b2g / 12),
                27 * (y_general + (a1g * x_general + a3g) / 2),
            )
        )
    if len(raw_points) != 16:
        raise ArithmeticError("sampled fibre MW16 point list is incomplete")
    return {
        "raw_short_model": ["0", "0", "0", qtext(specialized_a), qtext(specialized_b)],
        "raw_generic_points": [point_record(point) for point in raw_points],
        "j_invariant": qtext(raw_curve.j_invariant()),
    }, None


def family_model(common, presentation, input_path: Path):
    pencil = presentation["pencil"]
    source_hash = canonical_sha256(
        {
            "presentation_id": presentation["presentation_id"],
            "A": pencil["A_coefficients_low_to_high"],
            "B": pencil["B_coefficients_low_to_high"],
        }
    )
    return common.FamilyModel(
        source=input_path.resolve(),
        source_sha256=source_hash,
        a_coefficients=tuple(Fraction(value) for value in pencil["A_coefficients_low_to_high"]),
        b_coefficients=tuple(Fraction(value) for value in pencil["B_coefficients_low_to_high"]),
        a_degree=8,
        b_degree=12,
        coordinate=f"{presentation['presentation_id']}:lambda",
        coefficient_source_keys=(
            "pencil.A_coefficients_low_to_high",
            "pencil.B_coefficients_low_to_high",
        ),
    )


def candidate_id(presentation_id: str, numerator: int, denominator: int) -> str:
    sign = "m" if numerator < 0 else "p"
    return f"{presentation_id}-{sign}{abs(numerator)}d{denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=300)
    parser.add_argument("--height-bucket-width", type=int, default=50)
    parser.add_argument("--keep-per-bucket", default="16,8,4")
    parser.add_argument("--finalists-per-presentation", type=int, default=12)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    keep = tuple(int(value) for value in args.keep_per_bucket.split(","))
    if args.height_bound < 1 or args.height_bucket_width < 1:
        raise SystemExit("height bounds must be positive")
    if len(keep) != len(PRIME_BLOCKS) or any(value < 1 for value in keep):
        raise SystemExit("--keep-per-bucket requires three positive entries")

    template = json.loads(args.input.read_text())
    if template.get("status") != "PASS_TARGET_FREE_A1_MW16_FAMILY_PRESENTATIONS":
        raise ArithmeticError("anonymous family template is not passing")
    if template.get("presentation_count") != 9 or template.get("fibration_count") != 5:
        raise ArithmeticError("anonymous A1/MW16 family partition changed")
    common = load_source("a1_mw16_nagao_common", COMMON)
    chord = load_source("a1_mw16_chord", CHORD)

    source = json.loads(MODEL.read_text())
    old_ring = PolynomialRing(QQ, "t")
    old_field = old_ring.fraction_field()
    old_a = old_ring(source["weierstrass_model"]["A_coefficients_low_to_high"])
    old_b = old_ring(source["weierstrass_model"]["B_coefficients_low_to_high"])
    old_curve = EllipticCurve(old_field, [old_a, old_b])
    old_basis = tuple(
        old_curve(
            polynomial_from_record(record["X"], old_ring),
            polynomial_from_record(record["Y"], old_ring),
        )
        for record in source["sections"]["records"]
    )

    records = []
    failures = []
    presentation_rows = []
    started = perf_counter()
    for presentation in template["presentations"]:
        row_started = perf_counter()
        model = family_model(common, presentation, args.input)
        table_blocks, rejected = common.build_residue_tables(model, PRIME_BLOCKS)
        survivors, stages = common.run_staged_sieve(
            numerator_bound=args.height_bound,
            denominator_bound=args.height_bound,
            table_blocks=table_blocks,
            keep_per_bucket=keep,
            bucket_width=args.height_bucket_width,
        )
        finalists = sorted(survivors, key=common.candidate_sort_key)[
            : args.finalists_per_presentation
        ]

        marking = presentation["source_marking"]
        trace_vector = vector(ZZ, marking["trace_section_basis_w"])
        trace = sum(
            (coefficient * point for coefficient, point in zip(trace_vector, old_basis) if coefficient),
            old_curve(0),
        )
        h, nx, m0, quartic, child_a, child_b = child_geometry(
            trace, old_a, old_b, old_ring, chord
        )
        section_vectors = (
            vector(ZZ, marking["new_zero_source_section_basis_coordinates"]),
        ) + tuple(
            vector(ZZ, row)
            for row in marking["generic_source_section_basis_coordinates"]
        )
        source_points = tuple(
            sum(
                (coefficient * point for coefficient, point in zip(section_vector, old_basis) if coefficient),
                old_curve(0),
            )
            for section_vector in section_vectors
        )
        base_maps = tuple(
            polynomial_from_record(record, old_ring)
            for record in marking["base_maps_lambda_of_old_t"]
        )
        if len(source_points) != 17 or len(base_maps) != 17:
            raise ArithmeticError("anonymous section marking is incomplete")

        successful_ids = []
        for finalist in finalists:
            finalist_record = common.candidate_record(finalist)
            numerator, denominator = map(int, finalist_record["projective_pair"])
            if denominator == 0:
                failures.append({
                    "presentation_id": presentation["presentation_id"],
                    "parameter": "infinity",
                    "reason": "INFINITY_CHART_NOT_IMPLEMENTED",
                })
                continue
            parameter = QQ(numerator) / QQ(denominator)
            specialization, failure = specialize(
                parameter=parameter,
                source_points=source_points,
                base_maps=base_maps,
                old_ring=old_ring,
                h=h,
                nx=nx,
                m0=m0,
                quartic=quartic,
                child_a=child_a,
                child_b=child_b,
            )
            if failure is not None:
                failures.append({
                    "presentation_id": presentation["presentation_id"],
                    "parameter": finalist_record["parameter"],
                    "reason": failure,
                })
                continue
            identifier = candidate_id(presentation["presentation_id"], numerator, denominator)
            successful_ids.append(identifier)
            records.append({
                "candidate_id": identifier,
                "presentation_id": presentation["presentation_id"],
                "fibration_id": presentation["fibration_id"],
                "parameter": finalist_record["parameter"],
                "projective_pair": [numerator, denominator],
                "nagao": finalist_record,
                "generic_height_gram": presentation["generic_height_gram"],
                "generic_rank": 16,
                **specialization,
            })
        presentation_rows.append({
            "presentation_id": presentation["presentation_id"],
            "fibration_id": presentation["fibration_id"],
            "pencil_sha256": model.source_sha256,
            "usable_prime_blocks": [list(block.keys()) for block in table_blocks],
            "rejected_primes": list(rejected),
            "stages": stages,
            "finalists": [common.candidate_record(candidate) for candidate in finalists],
            "successful_candidate_ids": successful_ids,
            "runtime_seconds": perf_counter() - row_started,
        })
        print(
            f"A1MW16SAMPLE|presentation={presentation['presentation_id']}|"
            f"finalists={len(finalists)}|specialized={len(successful_ids)}|status=PASS",
            flush=True,
        )

    groups = []
    for record in records:
        curve = EllipticCurve(QQ, [QQ(value) for value in record["raw_short_model"]])
        for group in groups:
            if record["j_invariant"] == group["j_invariant"] and curve.is_isomorphic(group["curve"]):
                group["candidate_ids"].append(record["candidate_id"])
                break
        else:
            groups.append({
                "j_invariant": record["j_invariant"],
                "curve": curve,
                "candidate_ids": [record["candidate_id"]],
            })
    exact_classes = []
    for index, group in enumerate(groups, 1):
        class_id = f"QISO-{index:04d}"
        candidate_ids = sorted(group["candidate_ids"])
        exact_classes.append({
            "class_id": class_id,
            "j_invariant": group["j_invariant"],
            "candidate_ids": candidate_ids,
        })
        for record in records:
            if record["candidate_id"] in candidate_ids:
                record["q_isomorphism_class_id"] = class_id

    payload = {
        "schema": "elliptic-curves.a1-mw16-target-free-parameter-candidates.v1",
        "status": "PASS_TARGET_FREE_A1_MW16_PARAMETER_CANDIDATES",
        "geometry": {
            "exact_fibration_count": 5,
            "coordinate_presentation_count": 9,
            "presentation_role": "nested coordinate-height search charts",
        },
        "search": {
            "projective_height_box": args.height_bound,
            "numerator_interval": [-args.height_bound, args.height_bound],
            "denominator_interval": [1, args.height_bound],
            "primitive_pairs_only": True,
            "height_bucket_width": args.height_bucket_width,
            "requested_prime_blocks": [list(block) for block in PRIME_BLOCKS],
            "keep_per_bucket": list(keep),
            "finalists_per_presentation": args.finalists_per_presentation,
        },
        "presentations": presentation_rows,
        "requested_finalist_count": sum(len(row["finalists"]) for row in presentation_rows),
        "successful_specialization_count": len(records),
        "structural_failure_count": len(failures),
        "structural_failures": failures,
        "exact_q_isomorphism_class_count": len(groups),
        "q_isomorphism_classes": exact_classes,
        "duplicate_q_isomorphism_groups": [
            row for row in exact_classes if len(row["candidate_ids"]) > 1
        ],
        "candidates": records,
        "next_gate": {
            "stage": "bounded_half_lattice_jump_recovery",
            "authorized_search": (
                "complete exact maximum-depth MW16 stratum at reduced-coordinate "
                "height at most 100000 and at most 15 seconds per chart"
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (args.input, MODEL, COMMON, CHORD, Path(__file__))
        },
        "runtime_seconds": perf_counter() - started,
        "claim_boundary": [
            "Candidate selection reads only anonymous family equations and fixed local-prime blocks.",
            "No known-record parameter, curve model, point, rank, jump size, or local control is loaded.",
            "Nagao scores schedule a bounded detector and provide no rank evidence.",
            "Each successful row has sixteen exact specialized generic points on its raw short model.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/build_a1_mw16_target_free_parameter_candidates.sage --check"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file():
            raise ArithmeticError("stored target-free candidate ledger is absent")
        stored = json.loads(args.output.read_text())
        replay = json.loads(encoded)
        for document in (stored, replay):
            document.pop("runtime_seconds", None)
            for row in document["presentations"]:
                row.pop("runtime_seconds", None)
                for stage in row["stages"]:
                    stage.pop("runtime_seconds", None)
                    stage.pop("parameters_per_second", None)
        if stored != replay:
            raise ArithmeticError("stored target-free candidate ledger differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(
        f"A1MW16SAMPLE|requested={payload['requested_finalist_count']}|"
        f"specialized={len(records)}|q_isomorphism_classes={len(groups)}|"
        f"failures={len(failures)}|output={relative(args.output)}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
