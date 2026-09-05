#!/usr/bin/env sage-python
"""Export exact finite-field section systems at a declared ``P.O``.

For the short singleton twist of arithmetic genus three, a section with
``P.O=k`` has coprime coordinates

    x=X/H^2, y=Y/H^3,

with ``deg(H)=k``, ``deg(X)<=6+2k``, and ``deg(Y)<=9+3k``.  On any chart at
which the section is affine, take H monic and fix the leading affine point.
The upper half of the equation recursively determines Y; the lower half is
exported to msolve.  The union of k+1 distinct smooth chart fibres is
exhaustive, since a degree-k H cannot vanish at all of them.

This is an exact mod-p search space.  Solver completion and lifting to
characteristic zero are separate steps.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import sys
from dataclasses import asdict
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"elliptic-curves/cas"))
from research_runtime.regulator import Surface
from production_search_gates import function_field_gate_record

COVERS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
)
MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
LABELS = (
    "074d9-orbit-04b07",
    "074d9-orbit-11a44",
    "074d9-orbit-11279",
    "074d9-orbit-080fa",
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--label", choices=LABELS, required=True)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--intersection", type=int, required=True)
parser.add_argument(
    "--allow-singular-reduction",
    action="store_true",
    help="permit a coefficient-wise valid but globally nongood discovery prime",
)
parser.add_argument("--reduction-only", action="store_true", help="export finite-field proof work despite a rational section exclusion")
parser.add_argument("--covers", type=Path, default=COVERS)
parser.add_argument("--model", type=Path, default=MODEL)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=ROOT / "artifacts/local/elkies-k3/r17-074d9-twist-section-ladder",
)
args = parser.parse_args()

prime = int(args.prime)
intersection = int(args.intersection)
if prime < 5 or not GF(prime).is_prime_field():
    raise ValueError("--prime must be an odd prime at least five")
if intersection < 0:
    raise ValueError("--intersection must be nonnegative")
field = GF(prime)

covers = json.loads(args.covers.read_text())
if covers.get("status") != "PASS_EXACT_COMPLETE_074D9_CROSS_FIBRE_BISECTION_TRANSFER":
    raise ValueError("unexpected cover-certificate status")
source = next(
    row
    for fibre in covers["fibres"]
    for row in fibre["records"]
    if row["label"] == args.label
)
model = json.loads(args.model.read_text())
if model.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
    raise ValueError("unexpected lineage-model status")
representative = model["representative"]
surface = Surface(representative["A_coefficients_low_to_high"],
                  representative["B_coefficients_low_to_high"],
                  source["branch_quadratic_coefficients_low_to_high"])
section_gate = function_field_gate_record(surface=surface, target_rank=1,
    search_limits={"finite_field_prime": prime, "chart_count": intersection+1})
if not args.reduction_only and not section_gate["search_budget_gate"]["bounded_search_authorized"]:
    raise SystemExit("EXCLUDED_BY_THEOREM before section equations: " +
                     ", ".join(section_gate["theorem_pruning"]["theorems"]))

base_ring = PolynomialRing(field, "t")
t = base_ring.gen()
q = base_ring(
    [
        field(QQ(value))
        for value in source["branch_quadratic_coefficients_low_to_high"]
    ]
)
A0 = base_ring(
    [field(QQ(value)) for value in representative["A_coefficients_low_to_high"]]
)
B0 = base_ring(
    [field(QQ(value)) for value in representative["B_coefficients_low_to_high"]]
)
base_discriminant = -field(16) * (field(4) * A0**3 + field(27) * B0**2)
coefficient_model_valid = (
    q.degree() == 2
    and q.is_squarefree()
    and A0.degree() == 8
    and B0.degree() == 12
)
good_reduction = (
    coefficient_model_valid
    and base_discriminant.degree() == 24
    and base_discriminant.is_squarefree()
    and q.gcd(base_discriminant).degree() == 0
)
if not coefficient_model_valid or (not good_reduction and not args.allow_singular_reduction):
    raise ArithmeticError("declared prime is not permitted for the singleton twist")
A = q**2 * A0
B = q**3 * B0
if A.degree() != 12 or B.degree() != 18:
    raise ArithmeticError("short twist lost the arithmetic-genus-three degrees")


def fibre_data(coefficient_a, coefficient_b):
    discriminant = -field(16) * (
        field(4) * coefficient_a**3 + field(27) * coefficient_b**2
    )
    points = []
    two_torsion = []
    for leading_x in field:
        rhs = leading_x**3 + coefficient_a * leading_x + coefficient_b
        if rhs == 0:
            two_torsion.append(int(leading_x))
        if rhs.is_square():
            leading_y = rhs.sqrt()
            points.append((leading_x, leading_y))
            if leading_y:
                points.append((leading_x, -leading_y))
    return discriminant, points, two_torsion


def transformed_model(chart_parameter):
    if chart_parameter is None:
        return A, B, field(A[12]), field(B[18]), "original_infinity"
    c = field(chart_parameter)
    transformed_a = base_ring(
        sum(A[index] * (c * t + 1) ** index * t ** (12 - index)
            for index in range(13))
    )
    transformed_b = base_ring(
        sum(B[index] * (c * t + 1) ** index * t ** (18 - index)
            for index in range(19))
    )
    return transformed_a, transformed_b, field(A(c)), field(B(c)), f"t={int(c)}"


charts = []
for chart_parameter in (None, *tuple(field)):
    chart_A, chart_B, leading_a, leading_b, chart_label = transformed_model(
        chart_parameter
    )
    discriminant, points, two_torsion = fibre_data(leading_a, leading_b)
    if discriminant and not two_torsion:
        charts.append(
            (chart_parameter, chart_label, chart_A, chart_B, leading_a, leading_b,
             discriminant, points)
        )
    if len(charts) == intersection + 1:
        break
if len(charts) != intersection + 1:
    raise ArithmeticError("not enough smooth rational chart fibres without 2-torsion")

tag = args.label.removeprefix("074d9-orbit-")
output_dir = (
    args.output_dir.resolve() / tag / f"p{prime}" / f"intersection-{intersection}"
)
output_dir.mkdir(parents=True, exist_ok=True)
x_degree = 6 + 2 * intersection
y_degree = 9 + 3 * intersection
systems = []
for chart_index, chart in enumerate(charts):
    (
        chart_parameter,
        chart_label,
        chart_A,
        chart_B,
        leading_a,
        leading_b,
        discriminant,
        points,
    ) = chart
    for point_index, (leading_x, leading_y) in enumerate(points):
        names = (
            tuple(f"h{index}" for index in range(intersection - 1, -1, -1))
            + tuple(f"x{index}" for index in range(x_degree - 1, -1, -1))
        )
        coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
        variables = coefficient_ring.gens_dict()
        polynomial_ring = PolynomialRing(coefficient_ring, "T")
        T = polynomial_ring.gen()
        H = T**intersection + sum(
            variables[f"h{index}"] * T**index for index in range(intersection)
        )
        X = coefficient_ring(leading_x) * T**x_degree + sum(
            variables[f"x{index}"] * T**index for index in range(x_degree)
        )
        AS = polynomial_ring([coefficient_ring(value) for value in chart_A])
        BS = polynomial_ring([coefficient_ring(value) for value in chart_B])
        rhs = X**3 + AS * X * H**4 + BS * H**6
        y_coefficients = {y_degree: coefficient_ring(leading_y)}
        for degree in range(2 * y_degree - 1, y_degree - 1, -1):
            index = degree - y_degree
            known = sum(
                y_coefficients[left] * y_coefficients[degree - left]
                for left in y_coefficients
                if degree - left in y_coefficients
            )
            y_coefficients[index] = (
                rhs[degree] - known
            ) / (2 * coefficient_ring(leading_y))
        Y = sum(y_coefficients[index] * T**index for index in range(y_degree + 1))
        residual = Y**2 - rhs
        equations = [coefficient_ring(residual[degree]) for degree in range(y_degree)]
        if not all(equation != 0 for equation in equations):
            raise ArithmeticError("recursive section system contains a zero equation")
        system_path = output_dir / f"chart-{chart_index:02d}-block-{point_index:03d}.ms"
        text = ",".join(names) + f"\n{prime}\n"
        text += ",\n".join(str(equation).replace("**", "^") for equation in equations)
        text += "\n"
        system_path.write_text(text)
        systems.append(
            {
                "chart_index": chart_index,
                "block_index": point_index,
                "leading_x_y": [int(leading_x), int(leading_y)],
                "path": str(system_path.relative_to(ROOT)),
                "sha256": digest(system_path),
            }
        )

record = {
    "schema": "elkies-k3.r17-074d9-twist-section-ladder-msolve-export.v1",
    "status": (
        "PASS_EXACT_COMPLETE_MODP_SECTION_SYSTEM_EXPORT"
        if good_reduction
        else "PASS_EXACT_COMPLETE_MODP_DISCOVERY_SYSTEM_EXPORT"
    ),
    "proof_boundary": (
        "The union of the displayed chart systems is the exact coefficient-wise "
        "mod-p section scheme at the declared P.O. Solver completion and "
        "characteristic-zero lifting are separate. A nongood prime is discovery-only."
    ),
    "exact_surface": asdict(surface),
    "section_gate": section_gate,
    "purpose": "finite_field_proof" if args.reduction_only else "rational_section_search",
    "label": args.label,
    "prime": prime,
    "good_reduction": bool(good_reduction),
    "search_role": "rank_certificate_sieve" if good_reduction else "discovery_only",
    "intersection_P_dot_O": intersection,
    "coordinate_degree_bounds_H_X_Y": [intersection, x_degree, y_degree],
    "variable_count_per_block": 6 + 3 * intersection,
    "equation_count_per_block": y_degree,
    "chart_count": len(charts),
    "chart_cover_reason": "degree-k H cannot vanish at k+1 distinct chart fibres",
    "charts": [
        {
            "chart_index": index,
            "original_fibre": chart[1],
            "a_b": [int(chart[4]), int(chart[5])],
            "discriminant": int(chart[6]),
            "affine_point_count": len(chart[7]),
            "rational_two_torsion_x": [],
        }
        for index, chart in enumerate(charts)
    ],
    "systems": systems,
    "inputs": {
        str(path.resolve().relative_to(ROOT)): digest(path)
        for path in (args.covers, args.model)
    },
}
record_path = output_dir / "export.json"
record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(
    f"R17074D9SECTIONLADDER|label={args.label}|p={prime}|P.O={intersection}"
    f"|charts={len(charts)}|blocks={len(systems)}|variables={6 + 3 * intersection}"
    f"|equations={y_degree}|output={record_path}|status=PASS_EXPORTED",
    flush=True,
)
