#!/usr/bin/env sage-python
"""Decode arbitrary-degree NS0024 joint RUR factors and certify one point.

This is the primary closed-point extractor for the resolved-depth13 joint
system.  It requires the fixed ``rur_anchor`` separator, a reduced
zero-dimensional slice, and a squarefree degree-D eliminant.  Each irreducible
factor is decoded in its own residue field ``GF(p^d)``; no degree is guessed
and no common splitting field is constructed.

Every coordinate is substituted into the original exported equations modulo
the factor.  The resulting maximal ideal is then handed to the independent
joint-GB source verifier, which checks the four sections, fibre profile,
absolute components, and full Gram before any compact point is emitted.
"""

import argparse
import ast
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
JOINT_VERIFIER = HERE / "extract_lattice_foundry_ns0024_joint_gb_point.sage"
EDGE1_HANDOFF = HERE / "run_lattice_foundry_ns0024_edge1_handoff.sage"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(polynomial, value):
    answer = value.parent().zero()
    for coefficient in reversed(polynomial.list()):
        answer = answer * value + value.parent()(coefficient)
    return answer


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--system", type=Path, required=True)
parser.add_argument("--solution", type=Path, required=True, help="msolve RUR output")
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--source-output", type=Path)
parser.add_argument("--edge-output", type=Path)
parser.add_argument(
    "--no-edge1",
    action="store_true",
    help="diagnostic only: stop after the exact marked-point artifact",
)
args = parser.parse_args()
system_path = args.system.resolve()
solution_path = args.solution.resolve()
output_path = args.output.resolve()
if (args.source_output is None) != (args.edge_output is None):
    raise SystemExit("--source-output and --edge-output must be supplied together")
if args.no_edge1 and args.source_output is not None:
    raise SystemExit("explicit edge outputs cannot be combined with --no-edge1")

system_lines = system_path.read_text().splitlines()
names = tuple(system_lines[0].split(","))
prime = ZZ(system_lines[1])
if "rur_anchor" not in names:
    raise SystemExit("system was not exported with --fixed-rur-anchor")
anchor_index = names.index("rur_anchor")

solution_text = solution_path.read_text().strip()
parsed = ast.literal_eval(solution_text[:-1] if solution_text.endswith(":") else solution_text)
if parsed[0] != 0:
    raise ArithmeticError("msolve did not return a zero-dimensional RUR")
payload = parsed[1]
output_prime, variable_count, quotient_degree, output_names = payload[:4]
separating_vector = payload[4]
if int(output_prime) != prime or int(variable_count) != len(names):
    raise ArithmeticError("RUR header disagrees with the exported system")
if tuple(output_names) != names:
    raise ArithmeticError("RUR variable order changed")
if len(separating_vector) != len(names):
    raise ArithmeticError("RUR separating vector has the wrong length")
if any(
    int(value) % prime != (1 if index == anchor_index else 0)
    for index, value in enumerate(separating_vector)
):
    raise ArithmeticError("msolve did not use the fixed rur_anchor as separator")

parametrization = payload[5]
if parametrization[0] != 1:
    raise ArithmeticError("RUR has multiple parametrization blocks")
elimination_data, denominator_data, coordinate_data = parametrization[1]
base = GF(prime)
parameter_ring = PolynomialRing(base, "T")
T = parameter_ring.gen()
elimination = parameter_ring(elimination_data[1]).monic()
denominator = parameter_ring(denominator_data[1])
if elimination.degree() != int(quotient_degree):
    raise ArithmeticError("eliminant degree is not the quotient dimension")
if not elimination.is_squarefree():
    raise ArithmeticError("joint slice is not certified reduced/separated")
if denominator.gcd(elimination).degree() != 0:
    raise ArithmeticError("RUR denominator vanishes on the eliminant")

coordinate_polynomials = []
for block in coordinate_data:
    if len(block) != 1 or block[0][0] not in (-1, 0):
        raise ArithmeticError("unsupported multiblock coordinate numerator")
    coordinate_polynomials.append(
        None if block[0][0] == -1 else parameter_ring(block[0][1])
    )
if len(coordinate_polynomials) == len(names):
    coordinate_names = list(names)
elif len(coordinate_polynomials) == len(names) - 1:
    coordinate_names = [name for name in names if name != "rur_anchor"]
else:
    raise ArithmeticError("RUR coordinate count does not match the system")

factorization = tuple(elimination.factor())
if any(int(exponent) != 1 for unused, exponent in factorization):
    raise ArithmeticError("eliminant factorization is not squarefree")
print(
    "NS0024JOINTRUR|p={}|quotient={}|factors={}|degrees={}".format(
        prime,
        quotient_degree,
        len(factorization),
        ",".join(str(factor.degree()) for factor, unused in factorization),
    ),
    flush=True,
)

accepted = None
factor_failures = []
with tempfile.TemporaryDirectory(prefix="ns0024-joint-rur-") as temporary:
    temporary_path = Path(temporary)
    for factor_index, (factor, unused) in enumerate(factorization):
        factor = factor.monic()
        degree = int(factor.degree())
        if degree == 1:
            field = base
            root = -field(factor[0]) / field(factor[1])
            modulus_text = str(PolynomialRing(base, "z")([factor[0], factor[1]]))
        else:
            modulus_ring = PolynomialRing(base, "z")
            modulus = modulus_ring(factor.list()).monic()
            field = GF(prime**degree, "z", modulus=modulus)
            root = field.gen()
            modulus_text = str(modulus)
        if root ** (prime**degree) != root:
            raise ArithmeticError("residue-field generator fails Frobenius closure")
        if degree > 1 and len({root ** (prime**index) for index in range(degree)}) != degree:
            raise ArithmeticError("irreducible factor has a short Frobenius orbit")
        denominator_value = evaluate(denominator, root)
        if not denominator_value:
            raise ArithmeticError("RUR denominator vanished after factor restriction")
        assignments = {}
        for name, polynomial in zip(coordinate_names, coordinate_polynomials):
            # msolve records x_i*denominator + numerator = 0.
            assignments[name] = (
                field.zero() if polynomial is None
                else -evaluate(polynomial, root) / denominator_value
            )
        if "rur_anchor" not in assignments:
            assignments["rur_anchor"] = root
        separator_value = sum(
            field(separating_vector[index]) * assignments[name]
            for index, name in enumerate(names)
        )
        if separator_value != root:
            raise ArithmeticError("decoded coordinates disagree with the fixed separator")

        coordinate_ring = PolynomialRing(field, names=names, order="degrevlex")
        generators = coordinate_ring.gens_dict()
        equations = []
        for raw_line in system_lines[2:]:
            raw_line = raw_line.strip()
            if raw_line.endswith(","):
                raw_line = raw_line[:-1]
            if raw_line:
                equations.append(coordinate_ring(raw_line.replace("^", "**")))
        substitution = {generators[name]: assignments[name] for name in names}
        if any(equation.subs(substitution) for equation in equations):
            raise ArithmeticError("RUR factor fails exact substitution in the joint system")

        maximal_basis = [generators[name] - assignments[name] for name in names]
        gb_path = temporary_path / "factor-{}.gb".format(factor_index)
        gb_path.write_text(
            "[\n" + ",\n".join(str(polynomial) for polynomial in maximal_basis) + "\n]:\n"
        )
        candidate_path = temporary_path / "factor-{}.json".format(factor_index)
        command = [
            sys.executable,
            str(JOINT_VERIFIER),
            "--system", str(system_path),
            "--gb", str(gb_path),
            "--prime", str(prime),
            "--generator", "z",
            "--modulus", modulus_text,
            "--output", str(candidate_path),
        ]
        completed = subprocess.run(
            command, cwd=ROOT, text=True, capture_output=True, check=False
        )
        if completed.returncode:
            factor_failures.append(
                {
                    "factor": str(factor),
                    "degree": degree,
                    "verifier_tail": (completed.stdout + completed.stderr)[-1000:],
                }
            )
            continue
        accepted = json.loads(candidate_path.read_text())
        accepted["rur"] = {
            "quotient_dimension": int(quotient_degree),
            "elimination_polynomial": str(elimination),
            "selected_factor": str(factor),
            "selected_factor_degree": degree,
            "fixed_separator": "rur_anchor",
            "squarefree_degree_equals_quotient_dimension": True,
            "original_joint_equations_replayed": True,
            "frobenius_orbit_size": degree,
        }
        accepted["inputs"] = {
            "system": str(system_path),
            "system_sha256": digest(system_path),
            "rur_solution": str(solution_path),
            "rur_solution_sha256": digest(solution_path),
        }
        break

if accepted is None:
    print("NS0024JOINTRURFAILURES|{}".format(json.dumps(factor_failures)), flush=True)
    raise SystemExit("no RUR factor passes the exact resolved MW4 source marking")
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(accepted, indent=2, sort_keys=True) + "\n")
edge_text = "skipped"
if not args.no_edge1:
    source_output = (
        args.source_output.resolve()
        if args.source_output is not None
        else output_path.with_name(output_path.stem + "-source.json")
    )
    edge_output = (
        args.edge_output.resolve()
        if args.edge_output is not None
        else output_path.with_name(output_path.stem + "-edge1.json")
    )
    subprocess.run(
        [
            sys.executable,
            str(EDGE1_HANDOFF),
            "--input", str(output_path),
            "--source-output", str(source_output),
            "--edge-output", str(edge_output),
        ],
        cwd=ROOT,
        check=True,
    )
    edge_text = str(edge_output)
print(
    "NS0024JOINTRURPOINT|p={}|degree={}|output={}|edge1={}|status=PASS".format(
        prime, accepted["rur"]["selected_factor_degree"], output_path, edge_text
    ),
    flush=True,
)
