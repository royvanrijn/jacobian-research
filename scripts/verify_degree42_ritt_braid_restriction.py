#!/usr/bin/env python3
"""Verify completed section transport on the degree-42 237--327 Ritt edge.

The exact adjacent Ritt transition is reconstructed from the generic
``2 o 3 o 7`` polynomial by canonical ``3 o 14`` reconstruction followed
by canonical ``2 o 7`` reconstruction of the inner factor.  Linearizing
the resulting chart map in the seven graph-normal variables gives the
completed first-conormal transition.  The verifier then checks, for both
composite omissions 14 and 21, that this transition preserves the total
and spectator presentations and intertwines the explicit completed
sections.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import (  # noqa: E402
    W,
    Z,
    canonical_reconstruction,
)
from explore_degree42_ritt_rotated_conormal_flags import (  # noqa: E402
    build_chart,
    chart_coordinates,
    graph_normal_map,
)
from research_degree42_ritt_completed_presentations import (  # noqa: E402
    completed_presentation_output,
    presentation_cache,
)
from research_degree42_tensor_extension import parse_sections  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_braid_restriction.json"
)
SOURCE_WORD = (2, 3, 7)
TARGET_WORD = (3, 2, 7)
LAMBDA = {
    SOURCE_WORD: (
        "6*u^15",
        "6*u^8",
        "3*u",
        "1",
        "0",
        "0",
        "0",
    ),
    TARGET_WORD: (
        "12*u^22",
        "6*u^8",
        "2*u",
        "1",
        "0",
        "0",
        "0",
    ),
}
SECTIONS = {
    14: {
        "correction_component": 6,
        "correction": "-3*u^2+2*zeta",
    },
    21: {
        "correction_component": 7,
        "correction": "-4*u^3+8*u*zeta",
    },
}


def reconstructed_factors(
    polynomial: sp.Expr,
    outer_degree: int,
    inner_degree: int,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the canonically reconstructed normalized factor pair."""

    _, _, reconstruction, _ = canonical_reconstruction(
        polynomial,
        outer_degree,
        inner_degree,
    )
    inner = W**inner_degree + sum(
        reconstruction[sp.Symbol(f"b{outer_degree}{inner_degree}_{power}")]
        * W**power
        for power in range(1, inner_degree)
    )
    outer = Z**outer_degree + sum(
        reconstruction[sp.Symbol(f"a{outer_degree}{inner_degree}_{power}")]
        * Z**power
        for power in range(1, outer_degree)
    )
    return sp.expand(outer), sp.expand(inner)


def adjacent_transition() -> tuple[
    sp.Matrix,
    tuple[sp.Expr, sp.Expr],
]:
    """Return the first-conormal pullback and the two base corrections."""

    source_parameters, source_factor_variables, polynomial = build_chart(
        SOURCE_WORD
    )
    target_parameters, target_factor_variables = chart_coordinates(
        TARGET_WORD
    )
    target_outer, degree_fourteen_inner = reconstructed_factors(
        polynomial,
        3,
        14,
    )
    target_middle, target_inner = reconstructed_factors(
        degree_fourteen_inner,
        2,
        7,
    )
    target_factors = (
        sp.Poly(target_outer, Z),
        sp.Poly(target_middle, Z),
        sp.Poly(target_inner, W),
    )
    target_values = {
        variable: factor.nth(power)
        for factor, variables in zip(
            target_factors,
            target_factor_variables,
        )
        for power, variable in enumerate(variables, 1)
    }
    assert set(target_values) == set(target_parameters)

    source_normals, (tau, zeta), source_graph = graph_normal_map(
        SOURCE_WORD,
        source_factor_variables,
    )
    pulled_values = {
        variable: sp.cancel(
            value.subs(source_graph, simultaneous=True)
        )
        for variable, value in target_values.items()
    }
    inner_variables = target_factor_variables[-1]
    target_u = sp.cancel(pulled_values[inner_variables[-1]] / 7)
    target_tau = sp.cancel(target_u - 1)
    target_zeta = sp.cancel(
        (21 * target_u**2 - pulled_values[inner_variables[-2]]) / 7
    )

    target_normals, target_base, target_graph = graph_normal_map(
        TARGET_WORD,
        target_factor_variables,
    )
    zero_target_normals = {
        normal: sp.Integer(0) for normal in target_normals
    }
    target_base_graph = {
        variable: value.subs(zero_target_normals).subs(
            {
                target_base[0]: target_tau,
                target_base[1]: target_zeta,
            },
            simultaneous=True,
        )
        for variable, value in target_graph.items()
    }
    dependent_target_variables = tuple(
        variable
        for variables in target_factor_variables
        for variable in variables
        if variable not in {
            inner_variables[-1],
            inner_variables[-2],
        }
    )
    target_normal_values = tuple(
        sp.cancel(pulled_values[variable] - target_base_graph[variable])
        for variable in dependent_target_variables
    )
    zero_source_normals = {
        normal: sp.Integer(0) for normal in source_normals
    }
    assert all(
        sp.cancel(value.subs(zero_source_normals)) == 0
        for value in target_normal_values
    )
    base_corrections = (
        sp.cancel(target_tau - tau),
        sp.cancel(target_zeta - zeta),
    )
    conormal_jacobian = sp.Matrix(
        [
            [
                sp.cancel(
                    sp.diff(value, normal).subs(zero_source_normals)
                )
                for normal in source_normals
            ]
            for value in target_normal_values
        ]
    )
    # A target conormal covector pulls back by the transpose Jacobian.
    return conormal_jacobian.T, base_corrections


def singular(expression: sp.Expr) -> str:
    """Serialize a rational expression for Singular."""

    return str(sp.cancel(expression)).replace("**", "^")


def section_assignments(
    name: str,
    word: tuple[int, int, int],
    omission: int,
) -> str:
    """Return Singular assignments for one rank-one section matrix."""

    correction_row = SECTIONS[omission]["correction_component"]
    correction = SECTIONS[omission]["correction"]
    assignments = []
    for column, coefficient in enumerate(LAMBDA[word], 1):
        assignments.append(f"{name}[4,{column}]={coefficient};")
        assignments.append(
            f"{name}[{correction_row},{column}]="
            f"({correction})*({coefficient});"
        )
    return "".join(assignments)


def presentation(
    word: tuple[int, int, int],
    omission: int,
    section: str,
) -> list[str]:
    """Return one parsed total or spectator completed presentation."""

    output = completed_presentation_output(word, omission)
    sections = parse_sections(
        output,
        ("PRESENTATION_TOTAL", "PRESENTATION_SPECTATOR"),
    )
    result = sections[section]
    assert result
    return result


def transport_audit(pullback: sp.Matrix) -> dict[int, dict[str, int]]:
    """Check module and section transport for both composite omissions."""

    matrix_assignments = "".join(
        f"P[{row + 1},{column + 1}]={singular(pullback[row, column])};"
        for row in range(7)
        for column in range(7)
    )
    blocks = []
    for omission in SECTIONS:
        source_total = presentation(
            SOURCE_WORD,
            omission,
            "PRESENTATION_TOTAL",
        )
        source_spectator = presentation(
            SOURCE_WORD,
            omission,
            "PRESENTATION_SPECTATOR",
        )
        target_total = presentation(
            TARGET_WORD,
            omission,
            "PRESENTATION_TOTAL",
        )
        target_spectator = presentation(
            TARGET_WORD,
            omission,
            "PRESENTATION_SPECTATOR",
        )
        blocks.append(
            f"""
module PTS{omission}={",".join(source_total)};
module PSS{omission}={",".join(source_spectator)};
module PTT{omission}={",".join(target_total)};
module PST{omission}={",".join(target_spectator)};
matrix AS{omission}[7][7];
matrix AT{omission}[7][7];
{section_assignments(f"AS{omission}", SOURCE_WORD, omission)}
{section_assignments(f"AT{omission}", TARGET_WORD, omission)}
module totalRemainder{omission}=
  simplify(reduce(P*PTT{omission},std(PTS{omission})),2);
module spectatorRemainder{omission}=
  simplify(reduce(P*PST{omission},std(PSS{omission})),2);
module coherenceRemainder{omission}=simplify(
  reduce(
    P*AT{omission}-AS{omission}*P,
    std(PTS{omission})
  ),
  2
);
print("BRAID_RESTRICTION_{omission}");
print(size(totalRemainder{omission}));
print(size(spectatorRemainder{omission}));
print(size(coherenceRemainder{omission}));
"""
        )
    program = f"""
ring b=0,(tau,zeta),dp;
poly u=1+tau;
matrix P[7][7];
{matrix_assignments}
{"".join(blocks)}
"""
    singular_executable = shutil.which("Singular")
    assert singular_executable is not None, "Singular is required"
    result = subprocess.run(
        [singular_executable, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=600,
    )
    if "? error" in result.stdout:
        raise RuntimeError(result.stdout)
    compact = " ".join(
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip() and not line.strip().startswith("//")
    )
    audits = {}
    for omission in SECTIONS:
        match = re.search(
            rf"BRAID_RESTRICTION_{omission} "
            r"([0-9]+) ([0-9]+) ([0-9]+)",
            compact,
        )
        assert match is not None, result.stdout
        values = tuple(map(int, match.groups()))
        assert values == (0, 0, 0)
        audits[omission] = {
            "total_relation_remainder_generators": values[0],
            "spectator_relation_remainder_generators": values[1],
            "section_coherence_remainder_generators": values[2],
        }
    return audits


def main() -> None:
    pullback, base_corrections = adjacent_transition()
    tau, zeta = sp.symbols("tau zeta")
    constant_matrix = pullback.subs({tau: 0, zeta: 0})
    determinant_constant = sp.factor(constant_matrix.det())
    assert determinant_constant != 0
    assert all(correction == 0 for correction in base_corrections)
    assert pullback[3:, :] == sp.eye(7)[3:, :]
    audits = transport_audit(pullback)
    matrix_strings = [
        [singular(pullback[row, column]) for column in range(7)]
        for row in range(7)
    ]
    matrix_bytes = json.dumps(
        matrix_strings,
        separators=(",", ":"),
    ).encode()
    caches = [
        presentation_cache(word, omission)
        for omission in SECTIONS
        for word in (SOURCE_WORD, TARGET_WORD)
    ]
    output = {
        "schema": "degree42-ritt-braid-restriction.v1",
        "status": "exact completed adjacent-edge transport theorem",
        "edge": {
            "source_word": SOURCE_WORD,
            "target_word": TARGET_WORD,
            "move": "swap the adjacent coprime outer degrees 2 and 3",
            "construction": (
                "canonical 3-o-14 reconstruction followed by canonical "
                "2-o-7 reconstruction, then first-conormal linearization"
            ),
            "base_corrections": [
                singular(correction) for correction in base_corrections
            ],
            "pullback_rank": pullback.rank(),
            "constant_determinant": str(determinant_constant),
            "last_four_rows_identity": True,
            "matrix_sha256": hashlib.sha256(matrix_bytes).hexdigest(),
        },
        "presentation_caches": [
            {
                "path": str(cache.relative_to(ROOT)),
                "sha256": hashlib.sha256(cache.read_bytes()).hexdigest(),
            }
            for cache in caches
        ],
        "sectors": [
            {
                "composite_omission": omission,
                "source_section": (
                    f"e4+({SECTIONS[omission]['correction']})*"
                    f"e{SECTIONS[omission]['correction_component']}"
                ),
                "target_section": (
                    f"e4+({SECTIONS[omission]['correction']})*"
                    f"e{SECTIONS[omission]['correction_component']}"
                ),
                **audits[omission],
                "adjacent_restriction_class": "0",
            }
            for omission in SECTIONS
        ],
        "consequence": (
            "For omissions 14 and 21, the completed splitting is intrinsic "
            "across the 237--327 Ritt overlap: the total module, spectator "
            "quotient, and chosen section form a commuting restriction "
            "square."
        ),
        "sector_asymmetry": (
            "The chart transition does not exchange or remove the labelled "
            "correction terms. It fixes the last four normal covectors, so "
            "the cut-14 e6 correction and cut-21 e7 correction transport "
            "unchanged. The asymmetry is therefore labelled-sector data, "
            "not a mismatch between the 237 and 327 normalizations."
        ),
        "theorem_boundary": (
            "This is one adjacent edge of the six-vertex filled braid. "
            "Transport around the remaining five edges and the resulting "
            "closed braid/commuting-cell coherence classes are not computed "
            "here."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_braid_restriction.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: exact 237-to-327 completed conormal transition is invertible")
    print("PASS: cut-14 section has zero adjacent restriction class")
    print("PASS: cut-21 section has zero adjacent restriction class")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
