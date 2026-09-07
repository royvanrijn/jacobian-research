#!/usr/bin/env sage-python
"""Export the trace-zero genus-five R17 normalization scheme to msolve.

For a certified R17 model ``y^2=x^3+A(t)x+B(t)``, put

    x=x0+x1*t+...+x4*t^4.

The associated trace-zero bisection has rational normalization exactly when

    x^3+A*x+B = S(t)^2 * Q(t),

with monic ``deg(S)=5`` and squarefree ``deg(Q)=2``.  The ten residual
coefficient equations, saturated by the leading coefficient and optionally
the discriminant of Q, are the finite polynomial-x subchart. A common Q
attached to three independent x-polynomials would give the requested
twist-rank-three candidate; independence and characteristic-zero lifting
remain separate gates. General trace-zero bisections can have x=X/Q rather
than polynomial x, so this is not the complete ``P.O=0`` twist-section chart.
"""

# status: ACTIVE_COMPILER
# claim: export the finite polynomial-x trace-zero normalization subchart
# inputs: exact native 074d9 lineage certificate or explicit published model
# outputs: caller-selected msolve systems and compact export manifest

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--saturate-discriminant",
        action="store_true",
        help="exclude repeated Q inside msolve instead of filtering it afterward",
    )
    args = parser.parse_args()
    if not ZZ(args.prime).is_prime() or args.prime < 5:
        parser.error("--prime must be a prime at least five")

    model_path = args.model.resolve()
    document = json.loads(model_path.read_text())
    if document.get("status") == "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        equation = document["representative"]
        if equation.get("chart") != "norm12-orbit-074d9":
            raise ValueError("the exact lineage representative changed")
        coordinate = "norm12-orbit-074d9 native u"
    elif document.get("status") == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
        equation = document
        coordinate = document["coordinate"]
    else:
        raise ValueError("the input is not a supported certified R17 model")
    prime = int(args.prime)
    field = GF(prime)
    names = tuple(f"x{index}" for index in range(5)) + tuple(
        f"s{index}" for index in range(5)
    ) + (("vq", "vd") if args.saturate_discriminant else ("vq",))
    coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
    generators = coefficient_ring.gens()
    x_coefficients = generators[:5]
    s_coefficients = generators[5:10]
    vq = generators[10]
    vd = generators[11] if args.saturate_discriminant else None
    t_ring = PolynomialRing(coefficient_ring, "t")
    t = t_ring.gen()
    def coefficient_mod(value):
        value = QQ(value)
        denominator = int(value.denominator()) % prime
        if denominator == 0:
            raise ZeroDivisionError(f"model coefficient denominator vanishes at p={prime}")
        return field(int(value.numerator()) % prime) / field(denominator)

    A = t_ring([coefficient_mod(value) for value in equation["A_coefficients_low_to_high"]])
    B = t_ring([coefficient_mod(value) for value in equation["B_coefficients_low_to_high"]])
    x = sum(x_coefficients[index] * t**index for index in range(5))
    square_factor = t**5 + sum(
        s_coefficients[index] * t**index for index in range(5)
    )
    branch = x**3 + A * x + B
    q2 = coefficient_ring(branch[12])
    q1 = coefficient_ring(branch[11]) - 2 * s_coefficients[4] * q2
    q0 = (
        coefficient_ring(branch[10])
        - (s_coefficients[4] ** 2 + 2 * s_coefficients[3]) * q2
        - 2 * s_coefficients[4] * q1
    )
    quadratic = q0 + q1 * t + q2 * t**2
    residual = branch - square_factor**2 * quadratic
    equations = [coefficient_ring(residual[index]) for index in range(10)]
    discriminant = q1**2 - 4 * q0 * q2
    equations.append(vq * q2 - 1)
    if vd is not None:
        equations.append(vd * discriminant - 1)

    text = ",".join(names) + f"\n{prime}\n"
    text += ",\n".join(str(equation).replace("^", "^") for equation in equations)
    text += "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    metadata = args.metadata if args.metadata.is_absolute() else ROOT / args.metadata
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text)
    payload = {
        "schema": "elkies-k3.r17-tracezero-genus5-normalization-msolve-export.v1",
        "status": "PASS_FINITE_POLYNOMIAL_X_SUBCHART_EXPORTED",
        "prime": prime,
        "coordinate": coordinate,
        "variable_names": list(names),
        "equation_count": len(equations),
        "normalization": "S monic of degree five",
        "saturations": (
            ["q2 != 0", "q1^2-4*q0*q2 != 0"]
            if args.saturate_discriminant
            else ["q2 != 0"]
        ),
        "repeated_quadratic_filter": (
            "inside_msolve"
            if args.saturate_discriminant
            else "required_after_solution_extraction"
        ),
        "msolve_input": relative(output),
        "msolve_input_sha256": sha256(text.encode()).hexdigest(),
        "inputs": {relative(model_path): digest(model_path)},
        "proof_boundary": (
            "This is the polynomial-x trace-zero subchart modulo the displayed "
            "prime. It excludes degree drops and repeated quadratic branch divisors. "
            "It does not include the general x=X/Q form arising from a polynomial "
            "section on the short twist. Solving the modular system, treating all "
            "denominator and infinity charts, lifting to characteristic zero, hashing "
            "covers, and proving section independence are separate gates."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "export_r17_tracezero_genus5_normalizations_msolve.sage "
            f"--model {relative(model_path)} --prime {prime} --output {relative(output)} "
            f"--metadata {relative(metadata)}"
            + (" --saturate-discriminant" if args.saturate_discriminant else "")
        ),
    }
    metadata.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "R17TRACEZEROEXPORT"
        f"|p={prime}|variables={len(names)}|equations={len(equations)}"
        f"|output={relative(output)}|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
