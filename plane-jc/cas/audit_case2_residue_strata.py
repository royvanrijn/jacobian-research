#!/usr/bin/env python3
"""Exclude the decomposable Case-2 dicritical residues exactly.

For the published transformed ``(72,108)`` Case-2 pair, solve ``(J3),(J2)``
over the certified first-block field and solve the determined part of
``(J1)`` for ``G`` exactly as in the archived replay.  A hypothetical
dicritical normalization cover of degree ``delta`` means

    C(t) = Cbar(h(t)),  G(t) = Gbar(h(t))

for a polynomial right component ``h`` of degree ``delta``.  Since
``deg(C),deg(G)=(8,12)``, only ``delta=2,4`` are nontrivial.

This checker reconstructs the unique general monic right component from the
top coefficients of ``C``, computes the exact polynomial-decomposition
remainders for both ``C`` and ``G``, and asks Singular for a characteristic-
zero standard basis over ``K0=Q[u]/(H)``.  Both remainder ideals are the unit
ideal.  No ``(J1)`` compatibility equation and no ``(J0)`` equation is used.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from types import ModuleType


EXPECTED = {
    2: {
        "constraint_count": 9,
        "constraint_degrees": (14, 10, 6, 25, 21, 17, 13, 9, 5),
        "singular_input_sha256": (
            "8bd56f701554cfe909f5c5dbdc78139e"
            "68583043dc1783ce1bffafd48da0297a"
        ),
    },
    4: {
        "constraint_count": 12,
        "constraint_degrees": (
            14,
            12,
            10,
            25,
            23,
            21,
            17,
            15,
            13,
            9,
            7,
            5,
        ),
        "singular_input_sha256": (
            "7c4aac7c464441872edc525d8b337bf3"
            "53f427763cfcbcdcf17c4fd5254bc93a"
        ),
    },
}


@dataclass(frozen=True)
class ResidueStratumAudit:
    cover_degree: int
    outer_degrees: tuple[int, int]
    constraint_count: int
    constraint_degrees: tuple[int, ...]
    uses_j1_compatibility: bool
    uses_j0: bool
    singular_input_sha256: str
    unit_ideal: bool


def default_replay_root() -> Path:
    external = Path(__file__).resolve().parents[1] / "external" / "zenodo-21479814"
    matches = sorted(external.glob("*/release_bundle/exact_replay"))
    if len(matches) != 1:
        raise RuntimeError("pass the exact_replay directory explicitly")
    return matches[0]


def load_exact_core(replay_root: Path) -> ModuleType:
    path = replay_root / "exact_core.py"
    spec = importlib.util.spec_from_file_location("case2_residue_exact_core", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def solve_case2_through_j1(ec: ModuleType):
    """Return ``B,E,C,F,G`` and the unused J1 compatibility equations."""

    B, E, C, F, parameter_count = ec.high_layers()
    if parameter_count != 3:
        raise AssertionError("Case-2 high layers must leave three parameters")

    rows = [[ec.K(0) for _ in range(12)] for _ in range(20)]
    for j in range(1, 13):
        for i, coefficient in enumerate(ec.Avec):
            degree = i + j - 1
            if degree < 20:
                rows[degree][j - 1] += coefficient * (2 * j)

    known = ec.tadd(
        ec.tmul(B, ec.tder(F)),
        ec.tscale(ec.tmul(ec.tder(B), F), ec.K(-1)),
    )
    known = ec.tadd(
        known,
        ec.tscale(ec.tmul(ec.tder(C), E), ec.K(-2)),
    )
    kept_rows = []
    right_hand_sides = []
    for degree, row in enumerate(rows):
        value = known[degree] if degree < len(known) else ec.PP()
        if any(row) or value:
            kept_rows.append(row)
            right_hand_sides.append(-value)

    solved, free, compatibility, final_count = ec.linear_solve(
        kept_rows,
        right_hand_sides,
        parameter_count,
    )
    if free or final_count != 3:
        raise AssertionError("J1 must determine all twelve G coefficients")
    G = [ec.PP() for _ in range(13)]
    for j in range(1, 13):
        G[j] = solved[j - 1]
    return B, E, C, F, G, tuple(compatibility)


def subtract_scaled(ec: ModuleType, left, right, scale):
    return ec.tadd(left, ec.tscale(right, -scale))


def composition_remainders(
    ec: ModuleType,
    polynomial,
    right_component,
    outer_degree: int,
):
    """Divide by successive powers of a monic right component."""

    powers = {1: right_component}
    for exponent in range(2, outer_degree + 1):
        powers[exponent] = ec.tmul(
            powers[exponent - 1],
            right_component,
        )
    remainder = list(polynomial)
    right_degree = len(right_component) - 1
    for exponent in range(outer_degree, 0, -1):
        target_degree = exponent * right_degree
        coefficient = (
            remainder[target_degree]
            if target_degree < len(remainder)
            else ec.PP()
        )
        remainder = subtract_scaled(
            ec,
            remainder,
            powers[exponent],
            coefficient,
        )
    return tuple(
        remainder[degree]
        for degree in range(1, len(remainder))
        if remainder[degree]
    )


def right_component_constraints(ec: ModuleType, C, G, cover_degree: int):
    """Construct all positive-degree decomposition remainders."""

    if cover_degree == 2:
        # After making h monic and absorbing its constant in the outer
        # polynomials, h=t^2+a*t.  Since C is monic of degree eight,
        # C_7=4a.
        a = C[7] * ec.K(4).inv()
        right_component = [ec.PP(), a, ec.PP.const(1)]
        c_remainders = composition_remainders(
            ec,
            C,
            right_component,
            4,
        )
        g_remainders = composition_remainders(
            ec,
            G,
            right_component,
            6,
        )
    elif cover_degree == 4:
        # Write h=t^4+a*t^3+b*t^2+c*t.  The leading part
        # C=h^2+lambda*h+constant reconstructs a,b,c successively.
        inverse_two = ec.K(2).inv()
        a = C[7] * inverse_two
        b = (C[6] - a * a) * inverse_two
        c = (C[5] - 2 * a * b) * inverse_two
        right_component = [
            ec.PP(),
            c,
            b,
            a,
            ec.PP.const(1),
        ]
        c_remainders = composition_remainders(
            ec,
            C,
            right_component,
            2,
        )
        g_remainders = composition_remainders(
            ec,
            G,
            right_component,
            3,
        )
    else:
        raise ValueError("only the nontrivial cover degrees 2 and 4 apply")

    unique = []
    seen = set()
    for polynomial in c_remainders + g_remainders:
        serialized = polynomial.sing(3)
        if serialized not in seen:
            seen.add(serialized)
            unique.append(polynomial)
    return tuple(unique)


def singular_minpoly(ec: ModuleType) -> str:
    terms = []
    for exponent in range(len(ec.MOD)):
        coefficient = ec.MOD[exponent]
        if not coefficient:
            continue
        scalar = (
            str(int(coefficient.p))
            if coefficient.q == 1
            else f"({int(coefficient.p)}/{int(coefficient.q)})"
        )
        monomial = (
            ""
            if exponent == 0
            else ("u" if exponent == 1 else f"u^{exponent}")
        )
        if not monomial:
            term = scalar
        elif coefficient == 1:
            term = monomial
        elif coefficient == -1:
            term = f"-{monomial}"
        else:
            term = f"{scalar}*{monomial}"
        terms.append(term)
    return "+".join(terms).replace("+-", "-")


def singular_input(
    ec: ModuleType,
    constraints,
    cover_degree: int,
) -> str:
    marker = f"CASE2_COVER_{cover_degree}_UNIT"
    return "\n".join(
        (
            'LIB "resources.lib";',
            "Resources::setcores(1);",
            'LIB "nfmodstd.lib";',
            "ring R=(0,u),(r,s,h),dp;",
            f"minpoly={singular_minpoly(ec)};",
            "option(redSB);",
            "ideal I=",
            ",\n".join(polynomial.sing(3) for polynomial in constraints)
            + ";",
            "ideal J=nfmodStd(I);",
            f'if(size(J)>0 && J[1]==1){{print("{marker}");}}',
            "quit;",
            "",
        )
    )


def run_stratum(
    ec: ModuleType,
    C,
    G,
    cover_degree: int,
    singular: str,
    timeout: int,
) -> ResidueStratumAudit:
    constraints = right_component_constraints(ec, C, G, cover_degree)
    degrees = tuple(max(map(sum, polynomial.d)) for polynomial in constraints)
    expected = EXPECTED[cover_degree]
    if (
        len(constraints) != expected["constraint_count"]
        or degrees != expected["constraint_degrees"]
    ):
        raise AssertionError(
            f"unexpected cover-{cover_degree} decomposition constraints"
        )

    source = singular_input(ec, constraints, cover_degree)
    digest = hashlib.sha256(source.encode()).hexdigest()
    if digest != expected["singular_input_sha256"]:
        raise AssertionError(
            f"cover-{cover_degree} Singular input hash drifted"
        )
    marker = f"CASE2_COVER_{cover_degree}_UNIT"
    with tempfile.TemporaryDirectory(prefix="jc2-residue-") as directory:
        path = Path(directory) / f"case2_cover_{cover_degree}.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if completed.returncode != 0 or marker not in completed.stdout:
        raise RuntimeError(
            f"cover-{cover_degree} exact elimination failed\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )

    return ResidueStratumAudit(
        cover_degree=cover_degree,
        outer_degrees=(8 // cover_degree, 12 // cover_degree),
        constraint_count=len(constraints),
        constraint_degrees=degrees,
        uses_j1_compatibility=False,
        uses_j0=False,
        singular_input_sha256=digest,
        unit_ideal=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-root", type=Path)
    parser.add_argument("--singular", default=shutil.which("Singular"))
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    if not args.singular:
        raise RuntimeError("Singular is required for the exact residue audit")

    replay_root = (args.replay_root or default_replay_root()).resolve()
    ec = load_exact_core(replay_root)
    _, _, C, _, G, compatibility = solve_case2_through_j1(ec)
    if not compatibility:
        raise AssertionError("the unused J1 compatibility system is missing")

    audits = tuple(
        run_stratum(
            ec,
            C,
            G,
            cover_degree,
            args.singular,
            args.timeout,
        )
        for cover_degree in (2, 4)
    )
    for audit in audits:
        print(json.dumps(asdict(audit), sort_keys=True))
    print("CASE2_DECOMPOSABLE_RESIDUE_STRATA_PASS")


if __name__ == "__main__":
    main()
