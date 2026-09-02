#!/usr/bin/env sage-python
"""Lift the regular low-section slices through the modular D5 survivors.

status: ACTIVE_SEARCH
claim: exact slice reduction and p-adic lifting experiment only

For p=11 and p=13 the first twist section is restricted to

    x=u*(A+B*u),  y=C*u^2.

Its coefficient equations force A=-1 or A=2 (away from the double branch
intersection), and are rationalized by

    B=t^2-2,  C=t*(t^2-3).

After substituting the resulting monic quadratic twist, the second section
has six coefficients.  The seven section equations, including t, form a
square system.  This script checks the simple roots at p=11 and p=13, lifts
them uniquely p-adically, attempts coefficientwise rational reconstruction,
and exports the exact characteristic-zero systems for msolve.  The p=7
survivor is tested directly for lifting through p^3.  The two p=17 survivors
are classified on the forced d0=0 bad-fibre boundary and tested after fixing
their displayed k values.

A p-adic branch is not a rational lift.  Only literal substitution of a
reconstructed rational tuple changes that status.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ, ZZ, vector
from sage.arith.misc import rational_reconstruction


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3/d5-two-marked-two-twist-low-slices"


SLICE_RECORDS = (
    {
        "label": "Aminus1-p11",
        "branch": -1,
        "prime": 11,
        # t,x20,x21,x22,y20,y21,y22
        "point": (8, 1, 10, 3, 7, 8, 3),
        "d": (9, 9, 1),
        "sections": (
            {"x": [0, 10, 7], "y": [0, 0, 4]},
            {"x": [1, 10, 3], "y": [7, 8, 3]},
        ),
    },
    {
        "label": "A2-p13",
        "branch": 2,
        "prime": 13,
        "point": (2, 11, 2, 8, 5, 4, 3),
        "d": (8, 10, 1),
        "sections": (
            {"x": [0, 2, 2], "y": [0, 0, 2]},
            {"x": [11, 2, 8], "y": [5, 4, 3]},
        ),
    },
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def build_system(base_ring, branch, saturate_distinct=False):
    names = ["t", "x20", "x21", "x22", "y20", "y21", "y22"]
    if saturate_distinct:
        names.append("inverse")
    ring = PolynomialRing(
        base_ring,
        names=names,
        order="degrevlex",
    )
    variables = ring.gens()
    t, x20, x21, x22, y20, y21, y22 = variables[:7]
    u_ring = PolynomialRing(ring, "u")
    u = u_ring.gen()
    f = u**2*(-3*u**2-13*u-3)
    g = u**3*(2*u**3-17*u**2-12*u-2)
    x_value = x20+x21*u+x22*u**2
    y_value = y20+y21*u+y22*u**2
    denominator = t**2*(t**2-3)**2
    if branch == -1:
        scaled_twist = (
            denominator*u**2-t**2*(3*t**2+1)*u+1
        )
        d0_numerator = ring.one()
        d1_numerator = -t**2*(3*t**2+1)
    elif branch == 2:
        d0_numerator = 9*t**2-56
        d1_numerator = 6*t**4-37*t**2+27
        scaled_twist = (
            denominator*u**2+d1_numerator*u+d0_numerator
        )
    else:
        raise ValueError("branch must be -1 or 2")
    residual = scaled_twist*y_value**2-denominator*(x_value**3+f*x_value+g)
    equations = [ring(residual[degree]) for degree in range(7)]
    if saturate_distinct:
        # x20 != 0 separates the target from the tautological first-section
        # component.  t*(t^2-3) != 0 removes the denominator-clearing boundary.
        equations.append(variables[7]*x20*t*(t**2-3)-1)
    jacobian = Matrix(
        ring,
        [[equation.derivative(variable) for variable in variables]
         for equation in equations],
    )
    return {
        "ring": ring,
        "variables": variables,
        "equations": equations,
        "jacobian": jacobian,
        "twist_denominator": denominator,
        "d0_numerator": d0_numerator,
        "d1_numerator": d1_numerator,
    }


def evaluate_integer(polynomial, variables, values):
    return ZZ(polynomial.subs(dict(zip(variables, values))))


def reconstruct_tuple(residues, modulus):
    reconstructed = []
    for residue in residues:
        try:
            reconstructed.append(rational_reconstruction(residue % modulus, modulus))
        except (ArithmeticError, ValueError):
            return None
    return reconstructed


def hensel_square(record, exact, digits, checkpoints):
    prime = record["prime"]
    variables = exact["variables"]
    equations = exact["equations"]
    jacobian = exact["jacobian"]
    size = len(variables)
    if len(equations) != size:
        raise ValueError("Hensel system must be square")
    field = GF(prime)
    values = [ZZ(value % prime) for value in record["point"]]
    substitutions = dict(zip(variables, values))
    jacobian_mod_p = Matrix(
        field,
        [[jacobian[row, column].subs(substitutions)
          for column in range(size)] for row in range(size)],
    )
    determinant = int(jacobian_mod_p.det())
    if determinant > prime//2:
        determinant -= prime
    if determinant == 0:
        raise ArithmeticError(f"singular declared slice root: {record['label']}")
    inverse = jacobian_mod_p.inverse()
    modulus = ZZ(prime)
    reconstruction_rows = []
    for exponent in range(1, digits):
        residuals = [
            evaluate_integer(equation, variables, values)
            for equation in equations
        ]
        if any(residual % modulus for residual in residuals):
            raise ArithmeticError("Hensel divisibility invariant failed")
        correction = inverse*vector(
            field, [-ZZ(residual//modulus) for residual in residuals]
        )
        values = [
            values[index]+modulus*ZZ(correction[index])
            for index in range(size)
        ]
        modulus *= prime
        if exponent+1 not in checkpoints:
            continue
        candidate = reconstruct_tuple(values, modulus)
        exact_solution = False
        coefficient_bits = None
        if candidate is not None:
            candidate_substitutions = dict(zip(variables, candidate))
            exact_solution = all(
                equation.subs(candidate_substitutions) == 0
                for equation in equations
            )
            coefficient_bits = [
                max(abs(value.numerator()).nbits(),
                    abs(value.denominator()).nbits())
                for value in candidate
            ]
        reconstruction_rows.append(
            {
                "digits": exponent+1,
                "all_coordinates_reconstructed": candidate is not None,
                "coefficient_bits": coefficient_bits,
                "exact_rational_solution": exact_solution,
                "solution": [str(value) for value in candidate]
                    if exact_solution else None,
            }
        )
        if exact_solution:
            break
    return {
        "jacobian_rank": int(jacobian_mod_p.rank()),
        "jacobian_determinant_mod_p": determinant,
        "hensel_digits_requested": digits,
        "reconstruction_checkpoints": reconstruction_rows,
        "exact_rational_solution_found": any(
            row["exact_rational_solution"] for row in reconstruction_rows
        ),
    }


def hensel_lift(record, digits, checkpoints):
    return hensel_square(
        record, build_system(ZZ, record["branch"]), digits, checkpoints
    )


def msolve_text(system):
    names = [str(variable) for variable in system["variables"]]
    equations = [str(equation).replace("**", "^")
                 for equation in system["equations"]]
    return ",".join(names)+"\n0\n"+",\n".join(equations)+"\n"


def build_full_coefficient_system(base_ring):
    names = (
        "d1", "d0",
        "x12", "x11", "x10", "y12", "y11", "y10",
        "x22", "x21", "x20", "y22", "y21", "y20",
    )
    ring = PolynomialRing(base_ring, names=names, order="degrevlex")
    variables = ring.gens()
    named = ring.gens_dict()
    u_ring = PolynomialRing(ring, "u")
    u = u_ring.gen()
    d_value = u**2+named["d1"]*u+named["d0"]
    f = u**2*(-3*u**2-13*u-3)
    g = u**3*(2*u**3-17*u**2-12*u-2)
    equations = []
    for section_index in (1, 2):
        x_value = sum(
            named[f"x{section_index}{degree}"]*u**degree
            for degree in range(3)
        )
        y_value = sum(
            named[f"y{section_index}{degree}"]*u**degree
            for degree in range(3)
        )
        residual = x_value**3+f*x_value+g-d_value*y_value**2
        equations.extend(ring(residual[degree]) for degree in range(7))
    jacobian = Matrix(
        ring,
        [[equation.derivative(variable) for variable in variables]
         for equation in equations],
    )
    return {
        "ring": ring,
        "variables": variables,
        "equations": equations,
        "jacobian": jacobian,
    }


def finite_field_dependence_flags(prime, d_coefficients, sections):
    field = GF(prime)
    scalar_ring = PolynomialRing(field, "u")
    u = scalar_ring.gen()
    d_value = scalar_ring(d_coefficients)
    f = u**2*(-field(3)*u**2-field(13)*u-field(3))
    g = u**3*(field(2)*u**3-field(17)*u**2-field(12)*u-field(2))
    function_field = scalar_ring.fraction_field()
    curve = EllipticCurve(
        function_field, [0, 0, 0, d_value**2*f, d_value**3*g]
    )
    points = []
    for section in sections:
        x_value = scalar_ring(section["x"])
        y_value = scalar_ring(section["y"])
        points.append(curve(d_value*x_value, d_value**2*y_value))
    first, second = points
    return {
        "equal_or_opposite": bool(second == first or second == -first),
        "second_is_plus_or_minus_twice_first": bool(
            second == 2*first or second == -2*first
        ),
        "first_is_plus_or_minus_twice_second": bool(
            first == 2*second or first == -2*second
        ),
    }


def certify_p7_vertical_obstruction():
    prime = 7
    field = GF(prime)
    exact = build_full_coefficient_system(ZZ)
    variables = exact["variables"]
    equations = exact["equations"]
    point = [3, 6, 0, 3, 3, 3, 4, 6, 0, 4, 3, 3, 0, 1]
    substitutions = dict(zip(variables, point))
    jacobian = Matrix(
        field,
        [[exact["jacobian"][row, column].subs(substitutions)
          for column in range(14)] for row in range(14)],
    )
    if jacobian.rank() != 13:
        raise ArithmeticError("p=7 full coefficient Jacobian rank changed")
    right_kernel = jacobian.right_kernel().basis()[0]
    left_kernel = jacobian.left_kernel().basis()[0]
    values = vector(ZZ, point)
    residuals = vector(
        ZZ, [equation.subs(substitutions) for equation in equations]
    )
    if any(value % prime for value in residuals):
        raise ArithmeticError("declared p=7 point is not a solution")
    first_rhs = vector(field, [-ZZ(value//prime) for value in residuals])
    particular = jacobian.solve_right(first_rhs)
    first_lifts = []
    obstructions = []
    for parameter in field:
        correction = particular+parameter*right_kernel
        lifted = vector(
            ZZ,
            [values[index]+prime*ZZ(correction[index])
             for index in range(14)],
        )
        lifted_residuals = vector(
            ZZ,
            [equation.subs(dict(zip(variables, lifted)))
             for equation in equations],
        )
        if any(value % (prime**2) for value in lifted_residuals):
            raise ArithmeticError("invalid p=7 first Hensel lift")
        second_rhs = vector(
            field, [-ZZ(value//(prime**2)) for value in lifted_residuals]
        )
        obstruction = left_kernel.dot_product(second_rhs)
        first_lifts.append([int(value) for value in lifted])
        obstructions.append(int(obstruction))
    if obstructions != [1]*prime:
        raise ArithmeticError("p=7 second-order obstruction changed")
    dependence = finite_field_dependence_flags(
        prime,
        [6, 3, 1],
        (
            {"x": [3, 3, 0], "y": [6, 4, 3]},
            {"x": [3, 4, 0], "y": [1, 0, 3]},
        ),
    )
    if any(dependence.values()):
        raise ArithmeticError("p=7 survivor entered an obvious dependence component")
    return {
        "status": "VERTICAL_MODULAR_POINT_NO_LIFT_MOD_343",
        "prime": prime,
        "point_variable_order": [str(variable) for variable in variables],
        "point": point,
        "jacobian_rank": 13,
        "right_kernel": [int(value) for value in right_kernel],
        "left_kernel": [int(value) for value in left_kernel],
        "lifts_mod_49": len(first_lifts),
        "left_kernel_obstructions_mod_7": obstructions,
        "lifts_mod_343": 0,
        "squarefree_twist_discriminant_mod_7": 6,
        "distinct_section_pivot": "x11-x21=-1",
        "obvious_dependence_components": dependence,
        "constant_section_component": False,
        "claim": (
            "No coefficient tuple over Z_7 reducing to this ordered survivor "
            "solves the shared fourteen-equation ideal."
        ),
    }


def build_p17_boundary_system(base_ring, fixed_k=None):
    names = (
        ("a1", "b1", "c1", "e1", "a2", "b2", "c2", "e2")
        if fixed_k is not None else
        ("k", "a1", "b1", "c1", "e1", "a2", "b2", "c2", "e2")
    )
    ring = PolynomialRing(base_ring, names=names, order="degrevlex")
    variables = ring.gens()
    if fixed_k is None:
        k_value = variables[0]
        offset_values = (1, 5)
    else:
        k_value = ring(fixed_k)
        offset_values = (0, 4)
    equations = []
    for offset in offset_values:
        a_value, b_value, c_value, e_value = variables[offset:offset+4]
        equations.extend(
            (
                a_value**3-3*a_value-2-k_value*c_value**2,
                3*a_value**2*b_value-2*k_value*c_value*e_value
                    -c_value**2-13*a_value-3*b_value-12,
                3*a_value*b_value**2-k_value*e_value**2
                    -2*c_value*e_value-3*a_value-13*b_value-17,
                b_value**3-e_value**2-3*b_value+2,
            )
        )
    jacobian = Matrix(
        ring,
        [[equation.derivative(variable) for variable in variables]
         for equation in equations],
    )
    return {
        "ring": ring,
        "variables": variables,
        "equations": list(equations),
        "jacobian": jacobian,
    }


def certify_p17_boundary(digits, checkpoints):
    prime = 17
    field = GF(prime)
    declared = (
        {
            "label": "d0-boundary-k5-p17",
            "prime": prime,
            "k": 5,
            "curve_point": [5, 8, 16, 11, 2, 14, 7, 8, 1],
            "point": [8, 16, 11, 2, 14, 7, 8, 1],
            "sections": (
                {"x": [0, 8, 16], "y": [0, 11, 2]},
                {"x": [0, 14, 7], "y": [0, 8, 1]},
            ),
        },
        {
            "label": "d0-boundary-k8-p17",
            "prime": prime,
            "k": 8,
            "curve_point": [8, 4, 11, 6, 5, 11, 14, 3, 1],
            "point": [4, 11, 6, 5, 11, 14, 3, 1],
            "sections": (
                {"x": [0, 4, 11], "y": [0, 6, 5]},
                {"x": [0, 11, 14], "y": [0, 3, 1]},
            ),
        },
    )
    curve_system = build_p17_boundary_system(ZZ)
    curve_variables = curve_system["variables"]
    records = []
    scalar_ring = PolynomialRing(field, "u")
    u_mod = scalar_ring.gen()
    f_mod = u_mod**2*(-field(3)*u_mod**2-field(13)*u_mod-field(3))
    g_mod = u_mod**3*(field(2)*u_mod**3-field(17)*u_mod**2
                      -field(12)*u_mod-field(2))
    original_discriminant_mod = -field(16)*(4*f_mod**3+27*g_mod**2)
    residual_discriminant_mod = original_discriminant_mod//u_mod**7
    for source in declared:
        substitutions = dict(zip(curve_variables, source["curve_point"]))
        curve_jacobian = Matrix(
            field,
            [[curve_system["jacobian"][row, column].subs(substitutions)
              for column in range(9)] for row in range(8)],
        )
        if curve_jacobian.rank() != 8:
            raise ArithmeticError("p=17 boundary curve lost regularity")
        fixed_system = build_p17_boundary_system(ZZ, fixed_k=source["k"])
        lift_record = {key: source[key] for key in ("label", "prime", "point")}
        fixed_checkpoints = [value for value in checkpoints if value <= 400]
        lift = hensel_square(
            lift_record, fixed_system, min(digits, 400), fixed_checkpoints
        )
        dependence = finite_field_dependence_flags(
            prime, [0, source["k"], 1], source["sections"]
        )
        if any(dependence.values()):
            raise ArithmeticError("p=17 survivor entered an obvious dependence component")
        good_other_branch = residual_discriminant_mod(-field(source["k"]))
        if not good_other_branch:
            raise ArithmeticError("p=17 second branch point is not a good old fibre")
        records.append(
            {
                "label": source["label"],
                "k": source["k"],
                "boundary_curve_jacobian_rank": 8,
                "boundary_curve_tangent_dimension": 1,
                "fixed_k_lift": lift,
                "old_discriminant_unit_at_u_minus_k": int(good_other_branch),
                "obvious_dependence_components": dependence,
                "constant_section_component": False,
            }
        )

    coefficient_ring = PolynomialRing(QQ, "k")
    k = coefficient_ring.gen()
    u_ring = PolynomialRing(coefficient_ring, "u")
    u = u_ring.gen()
    f = u**2*(-3*u**2-13*u-3)
    g = u**3*(2*u**3-17*u**2-12*u-2)
    original_discriminant = -16*(4*f**3+27*g**2)
    minimal_a = (u+k)**2*(-3*u**2-13*u-3)
    minimal_b = (u+k)**3*(2*u**3-17*u**2-12*u-2)
    minimal_discriminant = -16*(4*minimal_a**3+27*minimal_b**2)
    if minimal_discriminant != (u+k)**6*original_discriminant//u**6:
        raise ArithmeticError("p=17 boundary minimal discriminant identity failed")
    if minimal_discriminant.valuation(u) != 1:
        raise ArithmeticError("p=17 boundary fibre at zero is not generically I1")
    return {
        "status": "BAD_FIBRE_BOUNDARY_EXCLUDED_FROM_K3_SOURCE_LOCUS",
        "forced_equation": "d0*c_i^2=0 for each section",
        "local_unit_condition": "c_i!=0 forces d0=0",
        "boundary_chart": {
            "d": "u*(u+k)",
            "x_i": "u*(a_i+b_i*u)",
            "y_i": "u*(c_i+e_i*u)",
            "equations_per_section": [
                "a^3-3*a-2-k*c^2=0",
                "3*a^2*b-2*k*c*e-c^2-13*a-3*b-12=0",
                "3*a*b^2-k*e^2-2*c*e-3*a-13*b-17=0",
                "b^3-e^2-3*b+2=0",
            ],
        },
        "minimal_twist_model": {
            "A": "(u+k)^2*(-3*u^2-13*u-3)",
            "B": "(u+k)^3*(2*u^3-17*u^2-12*u-2)",
            "degree_A": int(minimal_a.degree()),
            "degree_B": int(minimal_b.degree()),
            "affine_degree_discriminant": int(minimal_discriminant.degree()),
            "total_discriminant_degree": 12,
            "generic_fibre_at_u0": "I1",
            "generic_fibre_at_u_minus_k": "I0*",
            "generic_fibre_at_infinity": "I1",
            "chi": 1,
            "surface_type": "rational elliptic surface, not K3",
        },
        "records": records,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--digits", type=int, default=800)
parser.add_argument("--run-msolve", action="store_true")
parser.add_argument("--timeout", type=float, default=300.0)
parser.add_argument("--threads", type=int, default=4)
parser.add_argument(
    "--msolve", type=Path, default=Path(shutil.which("msolve") or "msolve")
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.digits < 2:
    raise ValueError("--digits must be at least two")

checkpoints = sorted(set(
    value for value in (20, 40, 80, 120, 200, 400, 800, arguments.digits)
    if value <= arguments.digits
))
LOCAL.mkdir(parents=True, exist_ok=True)
records = []
for declared in SLICE_RECORDS:
    record = dict(declared)
    exact_system = build_system(ZZ, record["branch"], saturate_distinct=True)
    system_path = LOCAL/f"{record['label']}.ms"
    solution_path = LOCAL/f"{record['label']}.solve"
    log_path = LOCAL/f"{record['label']}.log"
    content = msolve_text(exact_system)
    if arguments.check:
        if not system_path.exists() or system_path.read_text() != content:
            raise SystemExit(f"stale or missing msolve input: {system_path}")
    else:
        system_path.write_text(content)
    lift = hensel_lift(record, arguments.digits, checkpoints)
    dependence = finite_field_dependence_flags(
        record["prime"], record["d"], record["sections"]
    )
    if any(dependence.values()):
        raise ArithmeticError("low-section survivor entered an obvious dependence component")
    record["obvious_dependence_components"] = dependence
    record["constant_section_component"] = False
    msolve_record = {
        "classification": "not_run",
        "returncode": None,
        "timed_out": False,
        "system": str(system_path.relative_to(ROOT)),
        "system_sha256": digest(system_path),
        "solution": None,
        "log": None,
    }
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
            msolve_record["returncode"] = int(completed.returncode)
            msolve_record["classification"] = (
                "completed" if completed.returncode == 0 else "failed"
            )
        except subprocess.TimeoutExpired as error:
            output = error.stdout or ""
            if isinstance(output, bytes):
                output = output.decode(errors="replace")
            log_path.write_text(output)
            msolve_record["timed_out"] = True
            msolve_record["classification"] = "timeout"
        if solution_path.exists():
            msolve_record["solution"] = str(solution_path.relative_to(ROOT))
        if log_path.exists():
            msolve_record["log"] = str(log_path.relative_to(ROOT))
    record.update(lift)
    record["msolve"] = msolve_record
    records.append(record)

summary = {
    "schema": "elkies-k3.d5-two-marked-two-twist-low-slices.v1",
    "status": "EXACT_LOCAL_CLASSIFICATION_WITH_PADIC_EXPERIMENTS",
    "claim_boundary": (
        "The slice reductions, modular Jacobian ranks, p=7 obstruction "
        "through 7^3, and p=17 bad-fibre minimalization are exact. "
        "Each simple modular root has a unique p-adic lift in its slice. "
        "A p-adic lift is not a characteristic-zero rational point; only "
        "literal substitution of a reconstructed QQ tuple would prove one."
    ),
    "first_section_slice": {
        "x": "u*(A+B*u)",
        "y": "C*u^2",
        "forced_A_equation": "(A-2)*(A+1)^2=0",
        "parameterization": "B=t^2-2, C=t*(t^2-3)",
        "branches": {
            "A=-1": {
                "d0": "1/(t^2*(t^2-3)^2)",
                "d1": "-(3*t^2+1)/(t^2-3)^2",
            },
            "A=2": {
                "d0": "(9*t^2-56)/(t^2*(t^2-3)^2)",
                "d1": "(6*t^4-37*t^2+27)/(t^2*(t^2-3)^2)",
            },
        },
    },
    "records": records,
    "p7_full_coefficient_local_obstruction": certify_p7_vertical_obstruction(),
    "p17_bad_fibre_boundary": certify_p17_boundary(
        arguments.digits, checkpoints
    ),
}

GENERATED.mkdir(parents=True, exist_ok=True)
output_path = GENERATED/"elkies-k3-d5-two-marked-two-twist-low-slices-v1.json"
serialized = json.dumps(summary, indent=2, sort_keys=True)+"\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale or missing generated summary: {output_path}")
else:
    output_path.write_text(serialized)

for record in records:
    print(
        "D5_LOW_SLICE|"
        f"label={record['label']}|rank={record['jacobian_rank']}|"
        f"det={record['jacobian_determinant_mod_p']}|"
        f"qq={record['exact_rational_solution_found']}|"
        f"msolve={record['msolve']['classification']}"
    )
print(
    "D5_P7_LOCAL|lifts_mod_49=7|lifts_mod_343=0|"
    "status=VERTICAL_MODULAR_POINT_NO_LIFT"
)
for record in summary["p17_bad_fibre_boundary"]["records"]:
    print(
        "D5_P17_BOUNDARY|"
        f"k={record['k']}|rank={record['boundary_curve_jacobian_rank']}|"
        f"qq={record['fixed_k_lift']['exact_rational_solution_found']}|"
        "surface=RATIONAL_NOT_K3"
    )
print(output_path)
