#!/usr/bin/env python3
"""Generate modular inverse-pencil and Keller-target witnesses by cycle type."""

import argparse
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.chebotarev import (  # noqa: E402
    find_pencil_factorization_witness,
    integer_partitions,
    rational_good_reduction_certificate,
)


def parse_cycle_type(text):
    return tuple(int(part.strip()) for part in text.split(",") if part.strip())


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Find squarefree factorization types for H(w)-s*w+t over F_p and "
            "return balanced integer lifts and the C=1 Keller targets."
        )
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--degree",
        type=int,
        help="use the canonical inverse polynomial H=w^(N-1)*(1-w)",
    )
    source.add_argument(
        "--polynomial",
        help="a rational polynomial in w, for example 'w**4-w**5'",
    )
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument(
        "--cycle-type",
        action="append",
        default=[],
        help="comma-separated factor degrees; repeat the option as needed",
    )
    parser.add_argument(
        "--c",
        default="1",
        help="weighted-model constant c; the lifted C=1 target is (t/c,s,1)",
    )
    args = parser.parse_args()
    if not sp.isprime(args.prime):
        parser.error("--prime must be prime")

    w = sp.symbols("w")
    if args.degree is not None:
        if args.degree < 3:
            parser.error("--degree must be at least 3")
        H = sp.expand(w ** (args.degree - 1) * (1 - w))
    else:
        H = sp.expand(sp.sympify(args.polynomial, locals={"w": w}))
    c = sp.Rational(args.c)
    if c == 0:
        parser.error("--c must be nonzero")
    seed = sp.diff(H, w)
    admissibility_values = (
        H.subs(w, 0),
        H.subs(w, 1),
        seed.subs(w, 0),
        seed.subs(w, 1) + c,
    )
    if any(sp.cancel(value) != 0 for value in admissibility_values):
        parser.error(
            "H must satisfy H(0)=H(1)=H'(0)=0 and H'(1)=-c"
        )
    kappa = sp.cancel(sp.diff(seed, w).subs(w, 1) / c)
    if kappa == -2:
        parser.error("the weighted model is undefined when H''(1)/c=-2")
    a0 = sp.cancel(-(1 + kappa) / (2 + kappa))
    degree = sp.Poly(H, w, domain=sp.QQ).degree()
    requested = (
        [parse_cycle_type(text) for text in args.cycle_type]
        if args.cycle_type
        else list(integer_partitions(degree))
    )

    certificate = rational_good_reduction_certificate(H, w, c=c, a0=a0)
    if certificate["bad_integer"] % args.prime == 0:
        parser.error(
            f"prime {args.prime} divides the conservative bad-prime certificate; "
            "choose a certified good prime"
        )

    output = {
        "inverse_polynomial": str(H),
        "degree": degree,
        "prime": args.prime,
        "bad_integer": certificate["bad_integer"],
        "weighted_parameter_a0": str(a0),
        "witnesses": [],
    }
    for cycle_type in requested:
        witness = find_pencil_factorization_witness(
            H, w, args.prime, cycle_type
        )
        if witness is None:
            output["witnesses"].append(
                {"cycle_type": list(cycle_type), "found": False}
            )
            continue
        target_A = sp.cancel(sp.Rational(witness["lifted_intercept"]) / c)
        entry = {
            "cycle_type": list(witness["cycle_type"]),
            "found": True,
            "finite_field_pencil": {
                "s": witness["slope"],
                "t": witness["intercept"],
                "unit": witness["factorization_unit"],
                "factors": [str(factor) for factor in witness["factors"]],
            },
            "balanced_integer_lift": {
                "s": witness["lifted_slope"],
                "t": witness["lifted_intercept"],
                "inverse_polynomial": str(witness["lifted_polynomial"]),
            },
            "keller_target_C_equals_1": {
                "A": str(target_A),
                "B": witness["lifted_slope"],
                "C": 1,
            },
        }
        if witness["cycle_type"] == (degree,):
            entry["rational_certificate"] = (
                "irreducible over Q by irreducibility modulo p"
            )
        output["witnesses"].append(entry)

    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
