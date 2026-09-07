#!/usr/bin/env python3
"""Build source-pinned relative 2-Selmer jobs for the compact R17 family.

The generated Magma programs deliberately separate discovery from control
labelling.  They compute the complete 2-Selmer group unconditionally, map the
specialized generic sections into it, construct quotient representatives and
search their 2-covers *before* declaring any held-out exceptional points.
Only after the blind phase are the public control points mapped into the same
quotient basis.

This builder does not claim that Magma completed.  Its manifest records exact
inputs and generated program hashes; completed transcripts are interpreted by
``parse_elkies_2026_relative_2selmer_suite.py``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import importlib
import json
from pathlib import Path
import re
import sys
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
HIGH_RANK_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
RANK21_CONTROL = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve394_rank21_v1.json"
)
NAGAO_CANDIDATES = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_compact_t_nagao_positive_control_h10000_v1.json"
)
RECORD_PUBLIC_FIBRES = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
RECORD_LINEAGE_FIBRES = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
RECORD_RIGID_DIRECTIONS = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
)
SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-input.v1"
STATUS = "EXACT_INPUTS_MAGMA_COMPLETION_REQUIRED"
PROTOCOL = "ELKIESR17REL2"
GENERIC_RANK = 17

sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    is_on_weierstrass_curve,
    source_point_to_target,
    verify_finite_quotient_certificate,
)


Q = Fraction
PARAMETER_RE = re.compile(r"^(-?\d+)/(\d+)$")
CONTROL_MODULES = {
    "-2/377": "elkies_rank25",
    "-308/251": "elkies_rank26",
    "2456/135": "elkies_rank27",
    "-9529/5471": "elkies_rank28",
}


@dataclass(frozen=True)
class RelativeCase:
    case_id: str
    role: str
    parameter: str
    parameter_pair: tuple[int, int]
    model: tuple[Q, ...]
    generic_points: tuple[tuple[Q, Q], ...]
    exceptional_points: tuple[tuple[Q, Q], ...]
    certified_rank_lower_bound: int | None
    nagao_record: dict[str, Any] | None
    provenance: dict[str, str]
    rigid_quotient_rows: tuple[tuple[int, ...], ...] = ()
    rigid_direction_labels: tuple[str, ...] = ()
    rigid_complement_point_labels: tuple[str, ...] = ()


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_parameter(value: str) -> tuple[int, int]:
    match = PARAMETER_RE.fullmatch(value)
    if match is None:
        raise ValueError(f"parameter is not a normalized rational: {value}")
    numerator, denominator = (int(item) for item in match.groups())
    value_q = Q(numerator, denominator)
    if (value_q.numerator, value_q.denominator) != (numerator, denominator):
        raise ValueError(f"parameter is not reduced with positive denominator: {value}")
    return numerator, denominator


def rational_pair(row: Sequence[object]) -> tuple[Q, Q]:
    if len(row) != 2:
        raise ValueError("an elliptic point needs two affine coordinates")
    return Q(str(row[0])), Q(str(row[1]))


def _f2_rank(rows: Iterable[Sequence[int]]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def _leftmost_rref_pivots(rows: Iterable[Sequence[int]], width: int) -> tuple[int, ...]:
    matrix = [[int(bit) & 1 for bit in row] for row in rows]
    if any(len(row) != width for row in matrix):
        raise ValueError("a rigid quotient row has the wrong width")
    pivot_row = 0
    pivots = []
    for column in range(width):
        source = next(
            (index for index in range(pivot_row, len(matrix)) if matrix[index][column]),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        for index in range(len(matrix)):
            if index != pivot_row and matrix[index][column]:
                matrix[index] = [
                    left ^ right
                    for left, right in zip(matrix[index], matrix[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return tuple(pivots)


def load_record_pair_cases() -> list[RelativeCase]:
    """Load the two rank-29 ``074d9`` fibres with their rigid quotient plane."""

    public = json.loads(RECORD_PUBLIC_FIBRES.read_text())
    lineage = json.loads(RECORD_LINEAGE_FIBRES.read_text())
    rigid = json.loads(RECORD_RIGID_DIRECTIONS.read_text())
    if public.get("status") != "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES":
        raise ValueError("the pinned ICARM public-fibre corpus is not passing")
    if lineage.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ValueError("the exact R17 lineage certificate is not passing")
    if rigid.get("status") != "PASS_EXACT_COMPLETE_074D9_CROSS_FIBRE_BISECTION_TRANSFER":
        raise ValueError("the rigid-bisection transfer certificate is not passing")

    public_by_id = {int(row["id"]): row for row in public["records"]}
    independence_by_id = {
        int(row["curve_id"]): row for row in lineage["displayed_point_independence"]
    }
    quotient_by_id = {
        int(row["curve_id"]): row for row in lineage["exceptional_quotients"]
    }
    rigid_by_id = {int(row["curve_id"]): row for row in rigid["fibres"]}
    provenance = {
        str(path.relative_to(ROOT)): file_sha256(path)
        for path in (
            RECORD_PUBLIC_FIBRES,
            RECORD_LINEAGE_FIBRES,
            RECORD_RIGID_DIRECTIONS,
        )
    }

    cases = []
    for curve_id in (356, 385):
        source = public_by_id[curve_id]
        independence = independence_by_id[curve_id]
        quotient = quotient_by_id[curve_id]
        rigid_fibre = rigid_by_id[curve_id]
        if (
            source.get("representative") != "norm12-orbit-074d9"
            or int(source.get("snapshot_rank_lower_bound", -1)) != 29
            or int(independence.get("proved_displayed_subgroup_rank", -1)) != 29
            or quotient.get("specialized_generic_group_equals_first_seventeen_displayed_points")
            is not True
            or quotient.get("free_basis_modulo_generic_group")
            != [f"P{index}" for index in range(18, 30)]
            or int(quotient.get("free_rank", -1)) != 12
            or int(rigid_fibre.get("displayed_quotient_rank", -1)) != 12
            or int(rigid_fibre.get("class_span_rank", -1)) != 2
        ):
            raise ArithmeticError(f"record-fibre certificate mismatch for curve {curve_id}")
        model = tuple(Q(value) for value in source["ainvs"])
        points = tuple(rational_pair(point) for point in source["points"])
        if len(points) != 29 or any(
            not is_on_weierstrass_curve(model, point) for point in points
        ):
            raise ArithmeticError(f"invalid public point list for curve {curve_id}")
        signature_rows = [
            [int(bit) for bit in row]
            for signature in independence["mod2_reduction_signatures"]
            for row in signature["rows"]
        ]
        columns = [
            [row[index] for row in signature_rows] for index in range(len(points))
        ]
        if _f2_rank(columns) != 29:
            raise ArithmeticError(f"public rank-29 certificate failed for curve {curve_id}")
        rigid_rows = tuple(
            tuple(
                int(bit)
                for bit in record["finite_quotient_class_modulo_generic_17"][
                    "displayed_quotient_coordinates_over_f2"
                ]
            )
            for record in rigid_fibre["records"]
        )
        if len(rigid_rows) != 2 or _f2_rank(rigid_rows) != 2:
            raise ArithmeticError(f"rigid quotient plane changed for curve {curve_id}")
        pivots = set(_leftmost_rref_pivots(rigid_rows, 12))
        complement = tuple(
            f"P{18 + index}" for index in range(12) if index not in pivots
        )
        parameter = str(source["representative_parameter"]["affine_parameter"])
        pair = parse_parameter(parameter)
        cases.append(
            RelativeCase(
                case_id=f"record-r29-{curve_id}",
                role="record-rank29-residual-selmer-target",
                parameter=parameter,
                parameter_pair=pair,
                model=model,
                generic_points=points[:GENERIC_RANK],
                exceptional_points=points[GENERIC_RANK:],
                certified_rank_lower_bound=29,
                nagao_record=None,
                provenance=dict(provenance),
                rigid_quotient_rows=rigid_rows,
                rigid_direction_labels=tuple(
                    str(record["label"]) for record in rigid_fibre["records"]
                ),
                rigid_complement_point_labels=complement,
            )
        )
    return cases


def specialize(parameter: tuple[int, int]) -> tuple[tuple[Q, ...], tuple[tuple[Q, Q], ...]]:
    data = load_q12o5867_data(MODEL, SECTIONS)
    raw = evaluate_projective_specialization(data, *parameter)
    minimal_model, change, _metadata = global_minimal_model_with_change(raw.model)
    points = tuple(source_point_to_target(point, change) for point in raw.points)
    if len(points) != GENERIC_RANK:
        raise ArithmeticError("the compact specialization did not produce seventeen sections")
    if any(not is_on_weierstrass_curve(minimal_model, point) for point in points):
        raise ArithmeticError("a specialized generic section is not on the minimal model")
    return tuple(minimal_model), points


def _verify_combined_certificate(
    model: tuple[Q, ...],
    generic: tuple[tuple[Q, Q], ...],
    exceptional: tuple[tuple[Q, Q], ...],
    certificate: dict[str, Any],
) -> None:
    short_model, change = short_certificate_model(model)
    short_points = tuple(
        source_point_to_target(point, change) for point in generic + exceptional
    )
    verify_finite_quotient_certificate(short_model, short_points, certificate)
    expected = len(generic) + len(exceptional)
    if (
        certificate.get("relation_prime") != 2
        or certificate.get("combined_rank_over_relation_field") != expected
        or certificate.get("certified_independent") is not True
    ):
        raise ArithmeticError("the held-out control subgroup certificate is not full rank")


def load_rank21_case() -> RelativeCase:
    document = json.loads(RANK21_CONTROL.read_text())
    if document.get("status") != "PASS_EXACT_ICARM_CURVE394_RANK21_CONDUCTOR_REPLAY":
        raise ValueError("the rank-21 mechanism control is not passing")
    parameter = str(document["parameter"])
    pair = parse_parameter(parameter)
    model, generic = specialize(pair)
    expected_model = tuple(Q(value) for value in document["curve"]["global_minimal_model"])
    if model != expected_model:
        raise ArithmeticError("the compact t=3/8 fibre missed the pinned rank-21 model")
    public = tuple(rational_pair(point) for point in document["public_point_replay"]["points"])
    indices = tuple(
        int(index) - 1
        for index in document["rank_lower_bound"]["public_complement_indices_one_based"]
    )
    exceptional = tuple(public[index] for index in indices)
    _verify_combined_certificate(
        model,
        generic,
        exceptional,
        document["rank_lower_bound"]["finite_quotient_independence"],
    )
    return RelativeCase(
        case_id="control-r21-t3_8",
        role="held-out-positive-control",
        parameter=parameter,
        parameter_pair=pair,
        model=model,
        generic_points=generic,
        exceptional_points=exceptional,
        certified_rank_lower_bound=21,
        nagao_record=None,
        provenance={
            str(RANK21_CONTROL.relative_to(ROOT)): file_sha256(RANK21_CONTROL),
            str(MODEL.relative_to(ROOT)): file_sha256(MODEL),
            str(SECTIONS.relative_to(ROOT)): file_sha256(SECTIONS),
        },
    )


def load_high_rank_cases() -> list[RelativeCase]:
    document = json.loads(HIGH_RANK_CONTROLS.read_text())
    if document.get("status") != "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS":
        raise ValueError("the rank-25--28 control suite is not passing")
    cases: list[RelativeCase] = []
    for row in document["fibres"]:
        parameter = str(row["parameter"])
        pair = parse_parameter(parameter)
        model, generic = specialize(pair)
        expected_model = tuple(Q(value) for value in row["minimal_model"])
        if model != expected_model:
            raise ArithmeticError(f"compact specialization missed control {parameter}")
        module = importlib.import_module(CONTROL_MODULES[parameter])
        public = tuple(module.POINTS)
        indices = tuple(
            int(index) - 1
            for index in row["public_positive_control"][
                "selected_public_point_indices_one_based"
            ]
        )
        exceptional = tuple(public[index] for index in indices)
        if any(not is_on_weierstrass_curve(model, point) for point in exceptional):
            raise ArithmeticError(f"a public complement point misses control {parameter}")
        _verify_combined_certificate(
            model,
            generic,
            exceptional,
            row["public_positive_control"][
                "combined_generic_plus_complement_independence"
            ],
        )
        rank = int(row["locally_certified_rank_lower_bound"])
        cases.append(
            RelativeCase(
                case_id=f"control-r{rank}",
                role="held-out-positive-control",
                parameter=parameter,
                parameter_pair=pair,
                model=model,
                generic_points=generic,
                exceptional_points=exceptional,
                certified_rank_lower_bound=rank,
                nagao_record=None,
                provenance={
                    str(HIGH_RANK_CONTROLS.relative_to(ROOT)): file_sha256(
                        HIGH_RANK_CONTROLS
                    ),
                    str(MODEL.relative_to(ROOT)): file_sha256(MODEL),
                    str(SECTIONS.relative_to(ROOT)): file_sha256(SECTIONS),
                    str((CAS / f"{CONTROL_MODULES[parameter]}.py").relative_to(ROOT)): file_sha256(
                        CAS / f"{CONTROL_MODULES[parameter]}.py"
                    ),
                },
            )
        )
    return cases


def load_nagao_cases(limit: int) -> list[RelativeCase]:
    if limit < 0:
        raise ValueError("candidate limit must be nonnegative")
    document = json.loads(NAGAO_CANDIDATES.read_text())
    if document.get("status") != "PASS_POSITIVE_CONTROL_SCORING_GATE":
        raise ValueError("the frozen high-Nagao candidate list is not passing")
    cases: list[RelativeCase] = []
    for position, row in enumerate(document["finalists"][:limit], start=1):
        parameter = str(row["parameter"])
        pair = parse_parameter(parameter)
        model, generic = specialize(pair)
        cases.append(
            RelativeCase(
                case_id=f"nagao-{position:04d}",
                role="prospective-high-Nagao-candidate",
                parameter=parameter,
                parameter_pair=pair,
                model=model,
                generic_points=generic,
                exceptional_points=(),
                certified_rank_lower_bound=None,
                nagao_record=dict(row),
                provenance={
                    str(NAGAO_CANDIDATES.relative_to(ROOT)): file_sha256(NAGAO_CANDIDATES),
                    str(MODEL.relative_to(ROOT)): file_sha256(MODEL),
                    str(SECTIONS.relative_to(ROOT)): file_sha256(SECTIONS),
                },
            )
        )
    return cases


def magma_q(value: Q | int | str) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"({value.numerator}/{value.denominator})"


def magma_points(name: str, points: Sequence[tuple[Q, Q]]) -> str:
    if not points:
        return f"{name} := [ E | ];"
    rows = ",\n    ".join(
        f"E![{magma_q(x)}, {magma_q(y)}, 1]" for x, y in points
    )
    return f"{name} := [\n    {rows}\n];"


def magma_string(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def build_magma(
    case: RelativeCase,
    *,
    search_bound: int,
    enumerate_class_limit: int,
) -> str:
    if search_bound < 0 or enumerate_class_limit < 1:
        raise ValueError("search bound must be nonnegative and class limit positive")
    model = ", ".join(magma_q(value) for value in case.model)
    generic = magma_points("generic", case.generic_points)
    exceptional = magma_points("exceptional", case.exceptional_points)
    provenance_hash = sha256(
        json.dumps(case.provenance, sort_keys=True).encode()
    ).hexdigest()
    parameter = magma_string(case.parameter)
    case_id = magma_string(case.case_id)
    role = magma_string(case.role)
    return f'''// Generated by build_elkies_2026_relative_2selmer_suite.py
// Protocol {PROTOCOL} v1. Public exceptional points occur only after blind_end.
// Bound=-1 requests unconditional class groups; no GRH class-group bound is set.

Q := Rationals();
E := EllipticCurve([{model}]);
{generic}
assert #generic eq {GENERIC_RANK};
assert &and[P in E : P in generic];
T2, T2map := TwoTorsionSubgroup(E);
two_torsion_dim := #[n : n in Invariants(T2) | IsEven(n)];
printf "{PROTOCOL}|version=1|stage=input|case={case_id}|role={role}|parameter={parameter}|generic_count={GENERIC_RANK}|held_out_exceptional_count={len(case.exceptional_points)}|search_bound={search_bound}|enumerate_class_limit={enumerate_class_limit}|provenance_sha256={provenance_hash}|magma=%o\\n", GetVersion();

// Only the generic subgroup is supplied as a local-image hint.
generic_hints := {{ P[1] : P in generic }};
SetVerbose("TwoDescent", 1);
printf "{PROTOCOL}|stage=two_selmer|status=start|bound=-1|raw=true|hints=generic_only\\n";
two := MultiplicationByMMap(E, 2);
mu, cover_map := DescentMaps(two);
selmer_started := Realtime();
S2, AtoS, toVec, factor_base, returned_hints := SelmerGroup(
    two : Bound := -1, Raw := true, Hints := generic_hints
);
selmer_seconds := Realtime(selmer_started);
selmer_invariants := Invariants(S2);
assert &and[n eq 2 : n in selmer_invariants];
selmer_dim := #selmer_invariants;
F2 := GF(2);

function ClassBits(s)
    return Vector(F2, [ Integers() ! c : c in Eltseq(s) ]);
end function;

generic_classes := [ AtoS(mu(P)) : P in generic ];
generic_matrix := Matrix(F2, [ ClassBits(s) : s in generic_classes ]);
generic_kummer_rank := Rank(generic_matrix);
assert generic_kummer_rank eq {GENERIC_RANK};
residual_dim := selmer_dim - generic_kummer_rank;
assert residual_dim ge 0;
printf "{PROTOCOL}|stage=two_selmer|status=complete|seconds=%o|total_dim=%o|two_torsion_dim=%o|generic_kummer_rank=%o|residual_dim=%o|factor_base_size=%o\\n",
    selmer_seconds, selmer_dim, two_torsion_dim, generic_kummer_rank, residual_dim, #factor_base;
for i in [1..#generic_classes] do
    printf "{PROTOCOL}|stage=generic_class|index=%o|selmer_bits=%o\\n", i, Eltseq(generic_classes[i]);
end for;

// Extend the certified generic rows to a basis of the full Selmer space.
full_rows := [ generic_matrix[i] : i in [1..Nrows(generic_matrix)] ];
quotient_selmer_rows := [];
current_rank := Rank(Matrix(F2, full_rows));
for i in [1..selmer_dim] do
    candidate := Vector(F2, [ F2 | j eq i select 1 else 0 : j in [1..selmer_dim] ]);
    trial := Matrix(F2, full_rows cat [candidate]);
    if Rank(trial) gt current_rank then
        Append(~full_rows, candidate);
        Append(~quotient_selmer_rows, candidate);
        current_rank +:= 1;
    end if;
end for;
assert current_rank eq selmer_dim;
assert #quotient_selmer_rows eq residual_dim;
change_basis := Matrix(F2, full_rows);
change_basis_inverse := change_basis^-1;

function SelmerElementFromBits(bits)
    answer := S2 ! 0;
    for i in [1..selmer_dim] do
        if bits[i] eq 1 then answer +:= S2.i; end if;
    end for;
    return answer;
end function;

function QuotientBits(s)
    coordinates := ClassBits(s) * change_basis_inverse;
    return Vector(F2, [ coordinates[i] : i in [generic_kummer_rank+1..selmer_dim] ]);
end function;

quotient_basis := [ SelmerElementFromBits(row) : row in quotient_selmer_rows ];
for i in [1..#quotient_basis] do
    printf "{PROTOCOL}|stage=quotient_basis|index=%o|selmer_bits=%o|quotient_bits=%o|alpha=%o\\n",
        i, Eltseq(quotient_basis[i]), Eltseq(QuotientBits(quotient_basis[i])), quotient_basis[i] @@ AtoS;
end for;

nonzero_class_count := 2^residual_dim - 1;
enumerate_all := nonzero_class_count le {enumerate_class_limit};
blind_targets := [];
blind_target_labels := [];
if enumerate_all then
    for mask in [1..nonzero_class_count] do
        s := S2 ! 0;
        for i in [1..residual_dim] do
            if ((mask div 2^(i-1)) mod 2) eq 1 then s +:= quotient_basis[i]; end if;
        end for;
        Append(~blind_targets, s);
        Append(~blind_target_labels, Sprintf("all-%o", mask));
    end for;
else
    blind_targets := quotient_basis;
    blind_target_labels := [ Sprintf("basis-%o", i) : i in [1..residual_dim] ];
end if;
printf "{PROTOCOL}|stage=blind_plan|residual_dim=%o|nonzero_class_count=%o|enumerate_all=%o|target_count=%o\\n",
    residual_dim, nonzero_class_count, enumerate_all, #blind_targets;

// Construct and search quotient covers without the exceptional control points.
recovered_rows := [];
blind_started := Realtime();
for i in [1..#blind_targets] do
    s := blind_targets[i];
    alpha := s @@ AtoS;
    cover_started := Realtime();
    C, CtoE := TwoCover(alpha : E := E);
    construction_seconds := Realtime(cover_started);
    f, h := HyperellipticPolynomials(C);
    search_started := Realtime();
    points := Points(C : Bound := {search_bound});
    search_seconds := Realtime(search_started);
    if #points gt 0 then
        recovered := CtoE(points[1]);
        recovered_class := AtoS(mu(recovered));
        recovered_quotient := QuotientBits(recovered_class);
        Append(~recovered_rows, recovered_quotient);
        printf "{PROTOCOL}|stage=blind_cover|index=%o|label=%o|quotient_bits=%o|alpha=%o|quartic_f=%o|quartic_h=%o|construction_seconds=%o|search_seconds=%o|search_status=point_found|cover_point=%o|elliptic_point=%o|recovered_quotient_bits=%o\\n",
            i, blind_target_labels[i], Eltseq(QuotientBits(s)), alpha, f, h,
            construction_seconds, search_seconds, points[1], recovered,
            Eltseq(recovered_quotient);
    else
        printf "{PROTOCOL}|stage=blind_cover|index=%o|label=%o|quotient_bits=%o|alpha=%o|quartic_f=%o|quartic_h=%o|construction_seconds=%o|search_seconds=%o|search_status=no_point_within_bound\\n",
            i, blind_target_labels[i], Eltseq(QuotientBits(s)), alpha, f, h,
            construction_seconds, search_seconds;
    end if;
end for;
recovered_rank := #recovered_rows eq 0 select 0 else Rank(Matrix(F2, recovered_rows));
blind_seconds := Realtime(blind_started);
printf "{PROTOCOL}|stage=blind_end|seconds=%o|target_count=%o|recovered_class_count=%o|recovered_quotient_rank=%o\\n",
    blind_seconds, #blind_targets, #recovered_rows, recovered_rank;

// Held-out control labelling starts here. These points were not hints or seeds.
{exceptional}
assert &and[P in E : P in exceptional];
exceptional_classes := [ AtoS(mu(P)) : P in exceptional ];
exceptional_rows := [ QuotientBits(s) : s in exceptional_classes ];
exceptional_rank := #exceptional_rows eq 0 select 0 else Rank(Matrix(F2, exceptional_rows));
for i in [1..#exceptional_classes] do
    printf "{PROTOCOL}|stage=exceptional_class|index=%o|selmer_bits=%o|quotient_bits=%o\\n",
        i, Eltseq(exceptional_classes[i]), Eltseq(exceptional_rows[i]);
end for;
unexplained_dim := residual_dim - exceptional_rank;
assert unexplained_dim ge 0;
printf "{PROTOCOL}|stage=classification|status=complete|exceptional_count=%o|exceptional_quotient_rank=%o|known_realized_class_count=%o|unexplained_dim=%o|unrealized_class_count=%o|blind_recovered_rank=%o\\n",
    #exceptional, exceptional_rank, 2^exceptional_rank, unexplained_dim,
    2^residual_dim - 2^exceptional_rank, recovered_rank;
'''


def case_record(case: RelativeCase, program_path: Path, program: str) -> dict[str, Any]:
    return {
        "case_id": case.case_id,
        "role": case.role,
        "parameter": case.parameter,
        "projective_parameter": list(case.parameter_pair),
        "global_minimal_model": [str(value) for value in case.model],
        "generic_section_count": len(case.generic_points),
        "held_out_exceptional_point_count": len(case.exceptional_points),
        "certified_rank_lower_bound": case.certified_rank_lower_bound,
        "nagao_record": case.nagao_record,
        "rigid_quotient_rows": [list(row) for row in case.rigid_quotient_rows],
        "rigid_direction_labels": list(case.rigid_direction_labels),
        "rigid_complement_point_labels": list(case.rigid_complement_point_labels),
        "provenance": case.provenance,
        "program": str(program_path),
        "program_sha256": sha256(program.encode()).hexdigest(),
    }


def write_suite(
    cases: Iterable[RelativeCase],
    *,
    output_dir: Path,
    manifest: Path,
    search_bound: int,
    enumerate_class_limit: int,
    overwrite: bool,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for case in cases:
        program = build_magma(
            case,
            search_bound=search_bound,
            enumerate_class_limit=enumerate_class_limit,
        )
        path = output_dir / f"{case.case_id}.m"
        if path.exists() and not overwrite:
            raise FileExistsError(path)
        path.write_text(program)
        records.append(case_record(case, path, program))
        print(f"{PROTOCOL}|stage=wrote|case={case.case_id}|program={path}")
    document = {
        "schema": SCHEMA,
        "status": STATUS,
        "method": {
            "backend": "Magma SelmerGroup([2])/DescentMaps/TwoCover",
            "unconditional_class_group_request": "Bound=-1",
            "selmer_raw_basis": True,
            "generic_embedding": "AtoS(mu(P)) for each specialized generic section",
            "blind_phase_uses_public_exceptional_points": False,
            "cover_search_bound": search_bound,
            "enumerate_nonzero_quotient_class_limit": enumerate_class_limit,
            "large_quotient_policy": "construct and search a canonical quotient basis only",
        },
        "cases": records,
        "case_count": len(records),
        "claim_boundary": (
            "This manifest certifies exact generated inputs only. It contains no completed "
            "2-Selmer group, rank upper bound, or successful cover-search claim."
        ),
    }
    manifest.parent.mkdir(parents=True, exist_ok=True)
    if manifest.exists() and not overwrite:
        raise FileExistsError(manifest)
    manifest.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--candidate-count", type=int, default=10)
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--candidates-only", action="store_true")
    parser.add_argument("--record-pair-only", action="store_true")
    parser.add_argument("--search-bound", type=int, default=1000)
    parser.add_argument("--enumerate-class-limit", type=int, default=255)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if sum((args.controls_only, args.candidates_only, args.record_pair_only)) > 1:
        parser.error(
            "--controls-only, --candidates-only, and --record-pair-only are mutually exclusive"
        )
    cases: list[RelativeCase] = []
    if args.record_pair_only:
        cases = load_record_pair_cases()
    elif not args.candidates_only:
        cases = [load_rank21_case(), *load_high_rank_cases()]
    if not (args.controls_only or args.record_pair_only):
        cases.extend(load_nagao_cases(args.candidate_count))
    document = write_suite(
        cases,
        output_dir=args.output_dir,
        manifest=args.manifest,
        search_bound=args.search_bound,
        enumerate_class_limit=args.enumerate_class_limit,
        overwrite=args.overwrite,
    )
    print(
        f"{PROTOCOL}|stage=suite_complete|cases={document['case_count']}"
        f"|status={document['status']}|manifest={args.manifest}"
    )


if __name__ == "__main__":
    main()
