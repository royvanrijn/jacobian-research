#!/usr/bin/env sage-python
"""Search two polynomial twist sections on an exact two-marked D5 seed.

status: ACTIVE_SEARCH
claim: exact seed certificate plus modular experiment only

The rational elliptic surface is

    y^2 = x^3 + f*x + g,
    f = u^2*(-3*u^2-13*u-3),
    g = u^3*(2*u^3-17*u^2-12*u-2).

It has an I1* fibre at u=0 and two independent sections

    P=(-u-2*u^2, u^2),
    Q=(-u-u^2, u^2*(1-2*u)).

For every monic squarefree quadratic d over GF(p), the modular search is
complete for polynomial twist sections

    d*y_i^2 = x_i^3+f*x_i+g,  deg(x_i),deg(y_i) <= 2.

The script also exports the three inverse-variable modular Groebner charts
covering x_1 != x_2.  Use --run-msolve to attempt them.  A modular pair is
only a lifting candidate; the bounded relation scan is not an independence
certificate over the function field.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

from sage.all import EllipticCurve, GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL_ROOT = ROOT / "artifacts/local/elkies-k3/d5-two-marked-two-twist-polynomial"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def polynomial_coefficients(poly, length):
    return [int(poly[index]) for index in range(length)]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=11)
parser.add_argument("--relation-bound", type=int, default=12)
parser.add_argument("--run-msolve", action="store_true")
parser.add_argument("--timeout", type=float, default=60.0)
parser.add_argument("--threads", type=int, default=4)
parser.add_argument(
    "--msolve", type=Path, default=Path(shutil.which("msolve") or "msolve")
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

prime = int(arguments.prime)
field = GF(prime)
if prime < 5 or not field.is_prime_field():
    raise ValueError("--prime must be an odd prime at least five")
if arguments.relation_bound < 1:
    raise ValueError("--relation-bound must be positive")


# Exact characteristic-zero certificate for the section-first seed.
exact_ring = PolynomialRing(QQ, "u")
u_exact = exact_ring.gen()
f_exact = u_exact**2*(-3*u_exact**2-13*u_exact-3)
g_exact = u_exact**3*(2*u_exact**3-17*u_exact**2-12*u_exact-2)
P_exact = (-u_exact-2*u_exact**2, u_exact**2)
Q_exact = (-u_exact-u_exact**2, u_exact**2*(1-2*u_exact))
for x_value, y_value in (P_exact, Q_exact):
    if y_value**2 != x_value**3+f_exact*x_value+g_exact:
        raise ArithmeticError("exact marked-section identity failed")

discriminant_exact = -16*(4*f_exact**3+27*g_exact**2)
if discriminant_exact.valuation(u_exact) != 7:
    raise ArithmeticError("the declared D5 seed lost ord_u(Delta)=7")
if not (discriminant_exact//u_exact**7).is_squarefree():
    raise ArithmeticError("the residual finite discriminant is not squarefree")

# Independence specializes injectively for relations.  At u=1 the two points
# are a signed basis of an exactly rank-two elliptic curve over QQ.
specialized_curve = EllipticCurve(QQ, [0, 0, 0, -19, -29])
specialized_P = specialized_curve(-3, 1)
specialized_Q = specialized_curve(-2, -1)
if specialized_curve.rank(proof=True) != 2:
    raise ArithmeticError("specialized rank certificate changed")
specialized_generators = specialized_curve.gens(proof=True)
if not (
    specialized_P in specialized_generators
    and -specialized_Q in specialized_generators
):
    raise ArithmeticError("the specialized marked points are no longer the basis")


# Modular exhaustive fibrewise search.
scalar_ring = PolynomialRing(field, "u")
u = scalar_ring.gen()
f = u**2*(-field(3)*u**2-field(13)*u-field(3))
g = u**3*(field(2)*u**3-field(17)*u**2-field(12)*u-field(2))
field_values = list(field)
twists = []
for d1_value in field_values:
    for d0_value in field_values:
        d_value = u**2+d1_value*u+d0_value
        if d_value.is_squarefree():
            twists.append(d_value)

sections_by_twist = defaultdict(dict)
candidate_x_tested = 0
incidences_tested = 0
for x2_value in field_values:
    for x1_value in field_values:
        for x0_value in field_values:
            candidate_x_tested += 1
            x_value = x2_value*u**2+x1_value*u+x0_value
            numerator = x_value**3+f*x_value+g
            for d_value in twists:
                incidences_tested += 1
                quotient, remainder = numerator.quo_rem(d_value)
                if remainder or not quotient.is_square():
                    continue
                y_value = quotient.sqrt()
                if y_value.degree() > 2:
                    continue
                d_key = tuple(polynomial_coefficients(d_value, 3))
                x_key = tuple(polynomial_coefficients(x_value, 3))
                sections_by_twist[d_key][x_key] = polynomial_coefficients(y_value, 3)


relation_bound = int(arguments.relation_bound)
multiple_twists = []
function_field = scalar_ring.fraction_field()
for d_key, section_map in sorted(sections_by_twist.items()):
    if len(section_map) < 2:
        continue
    d_value = scalar_ring(list(d_key))
    twist_curve = EllipticCurve(
        function_field,
        [0, 0, 0, d_value**2*f, d_value**3*g],
    )
    section_rows = []
    points = []
    for x_key, y_coefficients in sorted(section_map.items()):
        x_value = scalar_ring(list(x_key))
        y_value = scalar_ring(y_coefficients)
        point = twist_curve(
            function_field(d_value*x_value),
            function_field(d_value**2*y_value),
        )
        points.append(point)
        section_rows.append(
            {"x": list(x_key), "y": y_coefficients}
        )
    relations = []
    first, second = points[:2]
    for first_coefficient in range(-relation_bound, relation_bound+1):
        for second_coefficient in range(-relation_bound, relation_bound+1):
            if first_coefficient == 0 and second_coefficient == 0:
                continue
            if first_coefficient*first+second_coefficient*second == twist_curve(0):
                relations.append([first_coefficient, second_coefficient])
    multiple_twists.append(
        {
            "d": list(d_key),
            "distinct_x_section_count": len(section_map),
            "sections": section_rows,
            "relations_in_box": relations,
        }
    )


# Export the direct modular Groebner systems.  These charts are deliberately
# retained even when fibrewise enumeration is faster, so larger sliced systems
# can reuse exactly the same coefficient ideal.
names = (
    "d1", "d0",
    "x12", "x11", "x10", "y12", "y11", "y10",
    "x22", "x21", "x20", "y22", "y21", "y20",
    "inverse",
)
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
variables = coefficient_ring.gens_dict()
u_ring = PolynomialRing(coefficient_ring, "u")
u_symbolic = u_ring.gen()
d_symbolic = u_symbolic**2+variables["d1"]*u_symbolic+variables["d0"]
f_symbolic = u_symbolic**2*(-field(3)*u_symbolic**2-field(13)*u_symbolic-field(3))
g_symbolic = u_symbolic**3*(
    field(2)*u_symbolic**3-field(17)*u_symbolic**2
    -field(12)*u_symbolic-field(2)
)
equations = []
for section_index in (1, 2):
    x_symbolic = sum(
        variables[f"x{section_index}{degree}"]*u_symbolic**degree
        for degree in range(3)
    )
    y_symbolic = sum(
        variables[f"y{section_index}{degree}"]*u_symbolic**degree
        for degree in range(3)
    )
    residual = x_symbolic**3+f_symbolic*x_symbolic+g_symbolic-d_symbolic*y_symbolic**2
    equations.extend(coefficient_ring(residual[degree]) for degree in range(7))

output_dir = LOCAL_ROOT/f"p{prime}"
output_dir.mkdir(parents=True, exist_ok=True)
groebner_records = []
for pivot in range(3):
    difference = variables[f"x1{pivot}"]-variables[f"x2{pivot}"]
    open_factor = variables["d1"]**2-field(4)*variables["d0"]
    saturation = coefficient_ring(variables["inverse"]*open_factor*difference-field(1))
    system_path = output_dir/f"pivot-x{pivot}.ms"
    solution_path = output_dir/f"pivot-x{pivot}.solve"
    log_path = output_dir/f"pivot-x{pivot}.log"
    with system_path.open("w") as stream:
        stream.write(",".join(names)+f"\n{prime}\n")
        exported = equations+[saturation]
        for equation_index, equation in enumerate(exported):
            stream.write(str(equation).replace("**", "^"))
            stream.write(",\n" if equation_index+1 < len(exported) else "\n")
    classification = "not_run"
    returncode = None
    timed_out = False
    if arguments.run_msolve:
        try:
            completed = subprocess.run(
                [
                    str(arguments.msolve.resolve()), "-f", str(system_path),
                    "-o", str(solution_path), "-t", str(arguments.threads),
                    "-v", "1",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=arguments.timeout,
                check=False,
            )
            log_path.write_text(completed.stdout)
            returncode = int(completed.returncode)
            if solution_path.exists():
                content = solution_path.read_text(errors="replace").strip()
                if content == "[-1]:":
                    classification = "empty_over_algebraic_closure"
                elif content.startswith("[0,"):
                    classification = "zero_dimensional"
                elif content.startswith("[1,"):
                    classification = "positive_dimensional"
                else:
                    classification = "unparsed"
            else:
                classification = "no_solution_file"
        except subprocess.TimeoutExpired as error:
            output = error.stdout or ""
            if isinstance(output, bytes):
                output = output.decode(errors="replace")
            log_path.write_text(output)
            timed_out = True
            classification = "timeout"
    groebner_records.append(
        {
            "pivot": f"x1{pivot}-x2{pivot}",
            "classification": classification,
            "timed_out": timed_out,
            "returncode": returncode,
            "system": str(system_path.relative_to(ROOT)),
            "system_sha256": digest(system_path),
            "solution": (
                str(solution_path.relative_to(ROOT)) if solution_path.exists() else None
            ),
            "log": str(log_path.relative_to(ROOT)) if log_path.exists() else None,
        }
    )


summary = {
    "status": "MODULAR_EXPERIMENT",
    "claim_boundary": (
        "The D5 seed and its two invariant sections are exact over QQ(u). "
        "The twist census is complete only over the displayed finite field "
        "for monic squarefree quadratic d and polynomial x,y of degree at most two."
    ),
    "prime": prime,
    "exact_seed": {
        "a_b_c": [-13, -17, -12],
        "discriminant_order_at_u0": 7,
        "residual_finite_discriminant_squarefree": True,
        "marked_sections": [
            {"x": [0, -1, -2], "y": [0, 0, 1, 0]},
            {"x": [0, -1, -1], "y": [0, 0, 1, -2]},
        ],
        "specialization_u": 1,
        "specialized_curve": "y^2=x^3-19*x-29",
        "specialized_rank": 2,
        "specialized_points_form_signed_basis": True,
    },
    "modular_ansatz": {
        "d": "u^2+d1*u+d0, squarefree",
        "section_equation": "d*y_i^2=x_i^3+f*x_i+g",
        "degree_x": 2,
        "degree_y": 2,
        "complete_over_GF_p": True,
    },
    "counts": {
        "squarefree_monic_quadratic_twists": len(twists),
        "candidate_x_polynomials": candidate_x_tested,
        "section_incidence_tests": incidences_tested,
        "section_incidences": sum(len(rows) for rows in sections_by_twist.values()),
        "twists_with_a_section": len(sections_by_twist),
        "twists_with_at_least_two_distinct_x_sections": len(multiple_twists),
    },
    "relation_scan_bound": relation_bound,
    "warning": (
        "No relation in the displayed coefficient box is not a proof of "
        "function-field independence.  A QQ lift and height determinant are required."
    ),
    "multiple_twists": multiple_twists,
    "modular_groebner_charts": groebner_records,
}

GENERATED.mkdir(parents=True, exist_ok=True)
output_path = GENERATED/f"elkies-k3-d5-two-marked-two-twist-polynomial-p{prime}-v1.json"
serialized = json.dumps(summary, indent=2, sort_keys=True)+"\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale or missing generated summary: {output_path}")
else:
    output_path.write_text(serialized)

print(
    "D5_TWO_MARKED_TWO_TWIST|"
    f"p={prime}|twists={len(twists)}|x={candidate_x_tested}|"
    f"incidences={sum(len(rows) for rows in sections_by_twist.values())}|"
    f"multi={len(multiple_twists)}|status=MODULAR_EXPERIMENT"
)
print(output_path)
