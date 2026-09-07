#!/usr/bin/env sage-python
"""Search a shared-simple-pole E6 two-twist-section chart over GF(p).

status: ACTIVE_SEARCH
claim: modular Groebner experiment only
inputs: --prime, --timeout, --threads
outputs: one generated summary and local msolve systems/solutions

The rank-two E6 rational surface is kept in the polynomial marked chart

    y^2 = x^3 + a*u*x + u*(u+c),
    c=2-a,  a*(r+1)=2*(r^2+r+1).

For a monic quadratic twist d and a shared monic simple pole H=u-h, two
candidate twist sections are imposed simultaneously as

    x_i=X_i/H^2,  y_i=Y_i/H^3,
    d*Y_i^2=X_i^3+a*u*X_i*H^4+u*(u+c)*H^6,

with deg(X_i)<=2 and deg(Y_i)<=3.  Fixing the sign at infinity by making
Y_i monic lets the coefficients of u^2,u,1 in Y_i be recovered recursively
from degrees 7,6,5.  Five residual equations remain per section.  Together
with the two marked-surface equations this is a square system in 12
variables.  Three saturation charts cover X_1 != X_2.

An empty modular chart is only an obstruction for solutions integral in that
normalization at the displayed prime.  A modular point is only a discovery
candidate until it is lifted and its height determinant is certified.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL_ROOT = ROOT / "artifacts/local/elkies-k3/e6-shared-pole-two-twist-sections"
GENERATED = ROOT / "artifacts/generated-results"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def classify_solution(path):
    if not path.exists():
        return "missing"
    content = path.read_text(errors="replace").strip()
    if content == "[-1]:":
        return "empty_over_algebraic_closure"
    if content.startswith("[1,"):
        return "positive_dimensional"
    if content.startswith("[0,"):
        return "zero_dimensional"
    return "unparsed"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=7)
parser.add_argument("--timeout", type=float, default=180.0)
parser.add_argument("--threads", type=int, default=4)
parser.add_argument("--msolve", type=Path, default=Path(shutil.which("msolve") or "msolve"))
parser.add_argument("--enumerate", action="store_true")
parser.add_argument("--skip-msolve", action="store_true")
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

prime = int(arguments.prime)
field = GF(prime)
if prime < 5 or not field.is_prime_field():
    raise ValueError("--prime must be an odd prime at least five")
if arguments.timeout <= 0 or arguments.threads <= 0:
    raise ValueError("--timeout and --threads must be positive")

names = (
    "a", "r", "d1", "d0", "h",
    "x12", "x11", "x10", "x22", "x21", "x20",
    "inverse",
)
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
variables = coefficient_ring.gens_dict()
a = variables["a"]
r = variables["r"]
d1 = variables["d1"]
d0 = variables["d0"]
h = variables["h"]
inverse = variables["inverse"]

u_ring = PolynomialRing(coefficient_ring, "u")
u = u_ring.gen()
d = u**2 + d1*u + d0
H = u-h
B = u*(u+field(2)-a)


def reduced_section_equations(prefix):
    X = (
        variables[f"x{prefix}2"]*u**2
        + variables[f"x{prefix}1"]*u
        + variables[f"x{prefix}0"]
    )
    higher_y = u**3
    recovered = {}
    # Degree 8 is the automatic monic identity.  Degrees 7,6,5 recover
    # y2,y1,y0 successively; their coefficients are all exactly 2.
    for degree, label in ((7, "y2"), (6, "y1"), (5, "y0")):
        residual = d*higher_y**2 - (
            X**3 + a*u*X*H**4 + B*H**6
        )
        value = -residual[degree] / field(2)
        recovered[label] = coefficient_ring(value)
        higher_y += value*u**(degree-5)
    residual = d*higher_y**2 - (X**3 + a*u*X*H**4 + B*H**6)
    if any(residual[degree] != 0 for degree in range(5, 9)):
        raise ArithmeticError("high-degree recursion did not close")
    equations = [coefficient_ring(residual[degree]) for degree in range(5)]
    return X, higher_y, recovered, equations


X1, Y1, recovered1, section1 = reduced_section_equations("1")
X2, Y2, recovered2, section2 = reduced_section_equations("2")
surface_equations = [a*(r+field(1))-field(2)*(r**2+r+field(1))]
equations = surface_equations + section1 + section2
if len(equations)+1 != len(names):
    raise ArithmeticError("shared-pole system is not square")

base_open = (
    r*(r-field(1))*(r+field(1))*(r+field(2))*(field(2)*r+field(1))
    * a*(d1**2-field(4)*d0)
)


enumeration = None
if arguments.enumerate:
    # Once r,d,H are fixed there are only p^3 possible X-polynomials.  The
    # high-degree recursion determines Y uniquely with the chosen sign.  This
    # is a complete finite-field enumeration of the same normalized ansatz and
    # supplies parameter fibres on which later Groebner lifting is much smaller.
    scalar_ring = PolynomialRing(field, "T")
    T = scalar_ring.gen()
    field_values = list(field)
    excluded_r = {
        field(0), field(1), -field(1), -field(2), -field(1)/field(2)
    }
    parameter_fibres = []
    sections_by_surface_twist = {}
    globals_tested = 0
    section_candidates_tested = 0
    for r_value in field_values:
        if r_value in excluded_r:
            continue
        a_value = field(2)*(r_value**2+r_value+field(1))/(r_value+field(1))
        if a_value == 0:
            continue
        c_value = field(2)-a_value
        for d1_value in field_values:
            for d0_value in field_values:
                if d1_value**2-field(4)*d0_value == 0:
                    continue
                d_value = T**2+d1_value*T+d0_value
                for h_value in field_values:
                    globals_tested += 1
                    H_value = T-h_value
                    B_value = T*(T+c_value)
                    sections = []
                    for x2_value in field_values:
                        for x1_value in field_values:
                            for x0_value in field_values:
                                section_candidates_tested += 1
                                X_value = (
                                    x2_value*T**2+x1_value*T+x0_value
                                )
                                Y_value = T**3
                                for degree in (7, 6, 5):
                                    residual_value = d_value*Y_value**2-(
                                        X_value**3
                                        + a_value*T*X_value*H_value**4
                                        + B_value*H_value**6
                                    )
                                    coefficient = -residual_value[degree]/field(2)
                                    Y_value += coefficient*T**(degree-5)
                                residual_value = d_value*Y_value**2-(
                                    X_value**3
                                    + a_value*T*X_value*H_value**4
                                    + B_value*H_value**6
                                )
                                if residual_value == 0:
                                    x_function = scalar_ring(X_value)/scalar_ring(H_value**2)
                                    y_function = scalar_ring(Y_value)/scalar_ring(H_value**3)
                                    canonical_key = (
                                        tuple(int(value) for value in x_function.numerator().list()),
                                        tuple(int(value) for value in x_function.denominator().list()),
                                        tuple(int(value) for value in y_function.numerator().list()),
                                        tuple(int(value) for value in y_function.denominator().list()),
                                    )
                                    sections.append(
                                        {
                                            "X": [int(x0_value), int(x1_value), int(x2_value)],
                                            "Y": [int(value) for value in Y_value.list()],
                                            "canonical_function_key": [list(part) for part in canonical_key],
                                        }
                                    )
                                    parameter_key = (
                                        int(r_value), int(a_value), int(c_value),
                                        int(d0_value), int(d1_value),
                                    )
                                    aggregate = sections_by_surface_twist.setdefault(
                                        parameter_key, {}
                                    )
                                    record = aggregate.setdefault(
                                        canonical_key,
                                        {
                                            "x_numerator": list(canonical_key[0]),
                                            "x_denominator": list(canonical_key[1]),
                                            "y_numerator": list(canonical_key[2]),
                                            "y_denominator": list(canonical_key[3]),
                                            "pole_presentations": [],
                                        },
                                    )
                                    record["pole_presentations"].append(int(h_value))
                    if len(sections) >= 2:
                        first, second = sections[:2]
                        substitutions = {
                            a: a_value,
                            r: r_value,
                            d1: d1_value,
                            d0: d0_value,
                            h: h_value,
                        }
                        for prefix, section in (("1", first), ("2", second)):
                            for index, value in enumerate(section["X"]):
                                substitutions[variables[f"x{prefix}{index}"]] = field(value)
                        core_variables = tuple(
                            variable for variable in coefficient_ring.gens()
                            if variable != inverse
                        )
                        if any(equation.subs(substitutions) != 0 for equation in equations):
                            raise ArithmeticError("enumerated section pair missed the square system")
                        jacobian = matrix(
                            field,
                            [
                                [
                                    equation.derivative(variable).subs(substitutions)
                                    for variable in core_variables
                                ]
                                for equation in equations
                            ],
                        )
                        parameter_fibres.append(
                            {
                                "r": int(r_value),
                                "a": int(a_value),
                                "c": int(c_value),
                                "d": [int(d0_value), int(d1_value), 1],
                                "h": int(h_value),
                                "section_count_with_monic_Y": len(sections),
                                "first_pair_jacobian_rank": int(jacobian.rank()),
                                "sections": sections,
                            }
                        )
    independent_pole_fibres = []
    for parameter_key, section_map in sorted(sections_by_surface_twist.items()):
        if len(section_map) < 2:
            continue
        r_value, a_value, c_value, d0_value, d1_value = parameter_key
        mechanism = "unclassified"
        mechanism_data = {}
        if d0_value == 0 and d1_value == c_value:
            mechanism = "dependent_d_equals_B"
            function_field = scalar_ring.fraction_field()
            d_polynomial = scalar_ring([d0_value, d1_value, 1])
            base_B = T*(T+field(c_value))
            twist_curve = EllipticCurve(
                function_field,
                [d_polynomial**2*field(a_value)*T, d_polynomial**3*base_B],
            )
            points = []
            for section in section_map.values():
                x_value = function_field(scalar_ring(section["x_numerator"])) / function_field(
                    scalar_ring(section["x_denominator"])
                )
                y_value = function_field(scalar_ring(section["y_numerator"])) / function_field(
                    scalar_ring(section["y_denominator"])
                )
                points.append(
                    twist_curve(d_polynomial*x_value, d_polynomial**2*y_value)
                )
            base_index = next(
                index for index, section in enumerate(section_map.values())
                if section["x_numerator"] in ([], [0])
                and section["x_denominator"] == [1]
                and section["y_numerator"] == [1]
                and section["y_denominator"] == [1]
            )
            base_point = points[base_index]
            relations = []
            for point in points:
                matches = [n for n in range(-6, 7) if n*base_point == point]
                relations.append(matches[0] if matches else None)
            mechanism_data = {"multiples_of_S_0_1": relations}
        elif d1_value == 0:
            constant_sections = []
            for section in section_map.values():
                if (
                    len(section["x_numerator"]) == 1
                    and section["x_denominator"] == [1]
                    and section["y_numerator"] == [1]
                    and section["y_denominator"] == [1]
                ):
                    k_value = field(section["x_numerator"][0])
                    if (
                        field(a_value)*k_value+field(c_value) == 0
                        and k_value**3 == field(d0_value)
                    ):
                        constant_sections.append(int(k_value))
            if constant_sections:
                mechanism = "constant_section_component"
                mechanism_data = {"constant_x_values": constant_sections}
        independent_pole_fibres.append(
            {
                "r": r_value,
                "a": a_value,
                "c": c_value,
                "d": [d0_value, d1_value, 1],
                "distinct_canonical_section_count": len(section_map),
                "mechanism": mechanism,
                "mechanism_data": mechanism_data,
                "sections": list(section_map.values()),
            }
        )
    if prime == 11 and any(
        row["mechanism"] == "unclassified" for row in independent_pole_fibres
    ):
        raise ArithmeticError("GF(11) shared-pole survivor mechanism changed")
    enumeration = {
        "complete_for_declared_modular_ansatz": True,
        "global_parameter_fibres_tested": globals_tested,
        "section_candidates_tested": section_candidates_tested,
        "fibres_with_at_least_two_sections": len(parameter_fibres),
        "parameter_fibres": parameter_fibres,
        "surface_twist_fibres_with_two_distinct_sections_across_poles": len(
            independent_pole_fibres
        ),
        "independent_pole_parameter_fibres": independent_pole_fibres,
    }

output_dir = LOCAL_ROOT / f"p{prime}"
output_dir.mkdir(parents=True, exist_ok=True)
msolve_path = arguments.msolve.resolve()
records = []
for pivot in range(3):
    difference = variables[f"x1{pivot}"]-variables[f"x2{pivot}"]
    saturation_equation = coefficient_ring(inverse*base_open*difference-field(1))
    system_path = output_dir / f"pivot-x{pivot}.ms"
    solution_path = output_dir / f"pivot-x{pivot}.solve"
    log_path = output_dir / f"pivot-x{pivot}.log"
    with system_path.open("w") as stream:
        stream.write(",".join(names)+f"\n{prime}\n")
        exported = equations+[saturation_equation]
        for index, equation in enumerate(exported):
            stream.write(str(equation).replace("**", "^"))
            stream.write(",\n" if index+1 < len(exported) else "\n")
    if arguments.skip_msolve:
        classification = "not_run"
        returncode = None
        timed_out = False
    elif arguments.check:
        classification = classify_solution(solution_path)
        returncode = 0 if solution_path.exists() else None
        timed_out = False
    else:
        try:
            completed = subprocess.run(
                [
                    str(msolve_path), "-f", str(system_path),
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
            timed_out = False
        except subprocess.TimeoutExpired as error:
            output = error.stdout or ""
            if isinstance(output, bytes):
                output = output.decode(errors="replace")
            log_path.write_text(output)
            returncode = None
            timed_out = True
        classification = "timeout" if timed_out else classify_solution(solution_path)
    records.append(
        {
            "pivot": f"x1{pivot}-x2{pivot}",
            "classification": classification,
            "timed_out": timed_out,
            "returncode": returncode,
            "system": str(system_path.relative_to(ROOT)),
            "system_sha256": digest(system_path),
            "solution": str(solution_path.relative_to(ROOT)) if solution_path.exists() else None,
            "solution_sha256": digest(solution_path) if solution_path.exists() else None,
            "log": str(log_path.relative_to(ROOT)) if log_path.exists() else None,
        }
    )

# Characteristic-zero rejection of the constant-section survivor mechanism.
# A constant twist point (k,1) on d=u^2-h^2 forces
# a*k+c=0, k^3=-h^2, a+c=2.  Writing k=-s^2 and h=s^3 gives the displayed
# rational parametrization.  Compatibility with the marked E6 rank-two
# relation is a quadratic in r with a strictly negative discriminant over QQ.
proof_ring = PolynomialRing(QQ, names=("proof_s", "proof_r"))
proof_s, proof_r = proof_ring.gens()
rank_two_compatibility = (
    (1+proof_s**2)*proof_r**2+proof_s**2*proof_r+proof_s**2
)
compatibility_discriminant = proof_s**4-4*(1+proof_s**2)*proof_s**2
expected_discriminant = -proof_s**2*(4+3*proof_s**2)
if compatibility_discriminant != expected_discriminant:
    raise ArithmeticError("constant-section compatibility discriminant changed")

payload = {
    "schema": "elkies-k3.e6-shared-pole-two-twist-sections-modp.v1",
    "status": "MODULAR_EXPERIMENT",
    "prime": prime,
    "surface": {
        "equation": "y^2=x^3+a*u*x+u*(u+2-a)",
        "marked_rank_two_relations": ["a*(r+1)=2*(r^2+r+1)"],
        "generic_fibres": "IV*+II+2I1",
    },
    "twist": "d=u^2+d1*u+d0",
    "shared_pole": "H=u-h",
    "section_ansatz": {
        "x_i": "X_i/H^2, deg(X_i)<=2",
        "y_i": "Y_i/H^3, deg(Y_i)<=3, Y_i monic",
        "recovered_y1": {key: str(value) for key, value in recovered1.items()},
        "recovered_y2": {key: str(value) for key, value in recovered2.items()},
    },
    "variable_count": len(names),
    "equation_count_before_saturation": len(equations),
    "saturation_open": (
        "a*r*(r-1)*(r+1)*(r+2)*(2*r+1)*(d1^2-4*d0)"
        "*(selected X1-X2 coefficient)"
    ),
    "charts": records,
    "finite_field_enumeration": enumeration,
    "exact_mechanism_rejections": {
        "dependent_d_equals_B": "the GF(11) canonical points are S and -2*S",
        "constant_section_parameterization": {
            "k": "-s^2",
            "h": "s^3",
            "a": "2/(1+s^2)",
            "c": "2*s^2/(1+s^2)",
            "rank_two_compatibility": str(rank_two_compatibility),
            "discriminant_in_r": str(compatibility_discriminant.factor()),
            "rational_nonzero_solutions": 0,
            "reason": "-s^2*(4+3*s^2) is negative for every nonzero rational s",
        },
    },
    "proof_boundary": (
        "Each completed chart is an exact msolve saturation computation over the displayed "
        "finite field for the declared shared-simple-pole ansatz. Empty charts do not exclude "
        "rational solutions nonintegral at this prime, other pole structures, or higher degree. "
        "Modular points do not prove two independent characteristic-zero twist directions."
    ),
}
summary_path = GENERATED / f"elkies-k3-e6-shared-pole-two-twist-sections-p{prime}-v1.json"
encoded = json.dumps(payload, indent=2, sort_keys=True)+"\n"
if arguments.check:
    if not summary_path.exists() or summary_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {summary_path}")
else:
    summary_path.write_text(encoded)

print(
    f"E6TWIST2|p={prime}|variables={len(names)}|equations={len(equations)}|"
    f"charts={','.join(row['classification'] for row in records)}|status=MODULAR_EXPERIMENT",
    flush=True,
)
print(f"OUTPUT|{summary_path}", flush=True)
