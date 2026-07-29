#!/usr/bin/env python3
"""Test degree-42 ideal and conormal extensions by tensor presentations.

Unlike image kernels after truncating the three quotient rings, this script
first presents ``K/I_6`` and ``K/I_boundary`` as modules and then tensors
those presentations with

    R / ((tau,zeta)^a + (n0,...,n6)^b).

Consequently a splitting of the completed R-module surjection would induce a
splitting in every displayed finite calculation.

The specialization ``b=1`` kills the full normal ideal ``K``.  It therefore
computes finite base quotients of the first cotangent homology modules

    K/(I_6 + K^2) -> K/(I_boundary + K^2),

rather than only the underlying ideal quotients.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from research_degree42_cellular_extension import (  # noqa: E402
    dickson_normal_map,
    extension_cocycle_certificate,
    serialize_polynomial,
    source_ideal_data,
    splitting_solution,
)


def parse_sections(
    output: str,
    markers: tuple[str, ...],
) -> dict[str, list[str]]:
    sections = {marker: [] for marker in markers}
    active: str | None = None
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if line in sections:
            active = line
            continue
        if active is None or not line or line.startswith("//"):
            continue
        sections[active].append(line)
    return sections


def vector_expression(
    text: str,
    symbols: tuple[sp.Symbol, ...],
    module_rank: int,
) -> tuple[sp.Expr, ...]:
    assert text.startswith("[") and text.endswith("]"), text
    entries = text[1:-1].split(",")
    names = {str(symbol): symbol for symbol in symbols}
    result = tuple(
        sp.expand(
            sp.sympify(
                re.sub(r"\^", "**", entry),
                locals=names,
            )
        )
        for entry in entries
    )
    assert len(result) <= module_rank
    return result + (sp.Integer(0),) * (module_rank - len(result))


def module_coordinate_matrix(
    expressions: list[str],
    basis_text: list[str],
    symbols: tuple[sp.Symbol, ...],
    module_rank: int,
) -> sp.Matrix:
    basis = tuple(
        vector_expression(text, symbols, module_rank)
        for text in basis_text
    )
    rank = module_rank
    basis_support: dict[tuple[int, tuple[int, ...]], int] = {}
    for index, vector in enumerate(basis):
        supports = []
        for component, polynomial in enumerate(vector):
            poly = sp.Poly(polynomial, *symbols)
            supports.extend(
                (component, monomial)
                for monomial, coefficient in poly.terms()
                if coefficient
            )
        assert len(supports) == 1, supports
        basis_support[supports[0]] = index

    matrix = sp.zeros(len(basis), len(expressions))
    for column, text in enumerate(expressions):
        vector = vector_expression(text, symbols, module_rank)
        assert len(vector) == rank
        for component, polynomial in enumerate(vector):
            for monomial, coefficient in sp.Poly(
                polynomial, *symbols
            ).terms():
                if not coefficient:
                    continue
                matrix[basis_support[(component, monomial)], column] = (
                    coefficient
                )
    return matrix


def singular_output(
    base_order: int,
    normal_order: int,
) -> tuple[str, tuple[sp.Symbol, ...]]:
    parameters, factor_variables, thick, boundary = source_ideal_data()
    normals, base_coordinates, images = dickson_normal_map(factor_variables)
    local_variables = normals + base_coordinates
    map_images = ",".join(
        serialize_polynomial(images[parameter]) for parameter in parameters
    )
    normal_module = ",".join(f"[{normal}]" for normal in normals)

    action_blocks = []
    markers = []
    for target, basis, presentation in (
        ("TOTAL", "BT", "PT"),
        ("SPECTATOR", "BS", "PS"),
    ):
        for variable in local_variables:
            marker = f"ACTION_{target}_{variable}"
            markers.append(marker)
            loop = f"index_{target.lower()}_{variable}"
            action_blocks.append(
                f"""
print("{marker}");
for (int {loop}=1; {loop}<=ncols({basis}); {loop}++)
{{
  print(reduce(({variable})*{basis}[{loop}],{presentation}));
}}
"""
            )

    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={",".join(thick)};
ideal IBsource={",".join(boundary)};
ring q=0,({",".join(map(str, local_variables))}),(dp({len(normals)}),dp(2));
map phi=source,{map_images};
option(redSB);
ideal IT=phi(ITsource);
ideal IB=phi(IBsource);
module MIT=IT;
module MIB=IB;
module MK={normal_module};
module RT=modulo(MK,MIT);
module RS=modulo(MK,MIB);
ideal baseIdeal={",".join(map(str, base_coordinates))};
ideal normalIdeal={",".join(map(str, normals))};
ideal F=baseIdeal^{base_order}+normalIdeal^{normal_order};
module FT=F*freemodule({len(normals)});
module PT=std(RT+FT);
module PS=std(RS+FT);
module BT=kbase(PT);
module BS=kbase(PS);
print("BASIS_TOTAL");
for (int index_total=1; index_total<=ncols(BT); index_total++)
{{
  print(BT[index_total]);
}}
print("BASIS_SPECTATOR");
for (int index_spectator=1; index_spectator<=ncols(BS); index_spectator++)
{{
  print(BS[index_spectator]);
}}
print("PROJECTION");
for (int index_projection=1; index_projection<=ncols(BT); index_projection++)
{{
  print(reduce(BT[index_projection],PS));
}}
{"".join(action_blocks)}
"""
    singular = shutil.which("Singular")
    assert singular is not None
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=7200,
    )
    return result.stdout, local_variables


def audit(base_order: int, normal_order: int) -> dict[str, object]:
    output, symbols = singular_output(base_order, normal_order)
    action_markers = tuple(
        f"ACTION_{target}_{variable}"
        for target in ("TOTAL", "SPECTATOR")
        for variable in symbols
    )
    markers = (
        "BASIS_TOTAL",
        "BASIS_SPECTATOR",
        "PROJECTION",
    ) + action_markers
    sections = parse_sections(output, markers)
    total_basis = sections["BASIS_TOTAL"]
    spectator_basis = sections["BASIS_SPECTATOR"]
    projection = module_coordinate_matrix(
        sections["PROJECTION"], spectator_basis, symbols, 7
    )
    total_actions = tuple(
        module_coordinate_matrix(
            sections[f"ACTION_TOTAL_{variable}"],
            total_basis,
            symbols,
            7,
        )
        for variable in symbols
    )
    spectator_actions = tuple(
        module_coordinate_matrix(
            sections[f"ACTION_SPECTATOR_{variable}"],
            spectator_basis,
            symbols,
            7,
        )
        for variable in symbols
    )

    assert projection.rank() == len(spectator_basis)
    assert all(
        projection * total_action == spectator_action * projection
        for total_action, spectator_action in zip(
            total_actions, spectator_actions
        )
    )
    assert all(
        left * right == right * left
        for actions in (total_actions, spectator_actions)
        for index, left in enumerate(actions)
        for right in actions[index + 1 :]
    )
    kernel = sp.Matrix.hstack(*projection.nullspace())
    kernel_actions = tuple(
        sp.Matrix.hstack(
            *(
                kernel.gauss_jordan_solve(
                    action * kernel[:, column]
                )[0]
                for column in range(kernel.cols)
            )
        )
        for action in total_actions
    )
    splits, section = splitting_solution(
        projection, total_actions, spectator_actions
    )
    certificate = extension_cocycle_certificate(
        kernel,
        projection,
        kernel_actions,
        total_actions,
        spectator_actions,
        tuple(map(str, symbols)),
    )
    zero_action_variables = [
        str(symbol)
        for symbol, total_action, spectator_action in zip(
            symbols, total_actions, spectator_actions
        )
        if total_action.is_zero_matrix and spectator_action.is_zero_matrix
    ]
    return {
        "base_order": base_order,
        "normal_order": normal_order,
        "dimensions": {
            "kernel": kernel.cols,
            "total": len(total_basis),
            "spectator": len(spectator_basis),
        },
        "splits_as_R_module": splits,
        "one_section": section.tolist() if section is not None else None,
        "zero_action_variables": zero_action_variables,
        "certificate": certificate,
    }


def conormal_audit(base_order: int = 2) -> dict[str, object]:
    """Audit the first Postnikov conormal projection modulo a base jet.

    Tensoring the presentations of ``K/I`` with

        R / (K + (tau,zeta)^base_order)

    gives ``K/(I+K^2)`` modulo the displayed base power.  These are the
    first homology modules of the two relative cotangent complexes for the
    quotient maps to ``R/K``.
    """

    result = audit(base_order=base_order, normal_order=1)
    result["interpretation"] = {
        "total": (
            "(K/(I_6+K^2)) tensor_B "
            "B/(tau,zeta)^base_order"
        ),
        "spectator": (
            "(K/(I_boundary+K^2)) tensor_B "
            "B/(tau,zeta)^base_order"
        ),
        "base_ring": "B=Q[[tau,zeta]]=completed R/K",
    }
    return result


def postnikov_overlap_audit() -> dict[str, object]:
    """Separate the completed conormal kernel from finite-base-change Tor.

    Put ``I=I_6`` and ``J=I_boundary``.  The first homology of the left
    transitivity term is

        S = J/(I+KJ).

    Its kernel in ``K/(I+K^2)`` is the quadratic overlap

        (J intersect (I+K^2))/(I+KJ).

    Both finite module lengths below are certified by an explicit
    Artin--Rees cutoff containment, rather than inferred from stabilization.
    """

    parameters, factor_variables, thick, boundary = source_ideal_data()
    normals, base_coordinates, images = dickson_normal_map(factor_variables)
    local_variables = normals + base_coordinates
    map_images = ",".join(
        serialize_polynomial(images[parameter]) for parameter in parameters
    )
    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={",".join(thick)};
ideal IBsource={",".join(boundary)};
ring q=0,({",".join(map(str, local_variables))}),(dp({len(normals)}),dp(2));
map phi=source,{map_images};
option(redSB);
ideal IT=phi(ITsource);
ideal IB=phi(IBsource);
ideal K={",".join(map(str, normals))};
ideal baseIdeal={",".join(map(str, base_coordinates))};
ideal maximalIdeal={",".join(map(str, local_variables))};

ideal sectorDenominator=IT+K*IB+baseIdeal^2*IB;
ideal sectorDenominatorStd=std(sectorDenominator);
ideal maximalIdeal4=maximalIdeal^4;
ideal sectorCutoffIntersection=intersect(IB,maximalIdeal4);
ideal sectorCutoffRemainder=simplify(
  reduce(sectorCutoffIntersection,sectorDenominatorStd),2
);
print("SECTOR_SOURCE");
print(size(sectorCutoffIntersection));
print(size(sectorCutoffRemainder));
print(vdim(std(sectorDenominator+maximalIdeal4)));
print(vdim(std(IB+maximalIdeal4)));

ideal overlapNumerator=intersect(IB,IT+K^2);
ideal overlapDenominator=IT+K*IB+baseIdeal^2*overlapNumerator;
ideal overlapDenominatorStd=std(overlapDenominator);
ideal maximalIdeal5=maximalIdeal^5;
ideal overlapCutoffIntersection=intersect(
  overlapNumerator,maximalIdeal5
);
ideal overlapCutoffRemainder=simplify(
  reduce(overlapCutoffIntersection,overlapDenominatorStd),2
);
print("QUADRATIC_OVERLAP");
print(size(overlapNumerator));
print(size(overlapCutoffIntersection));
print(size(overlapCutoffRemainder));
print(vdim(std(overlapDenominator+maximalIdeal5)));
print(vdim(std(overlapNumerator+maximalIdeal5)));
"""
    singular = shutil.which("Singular")
    assert singular is not None
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=7200,
    )
    compact = " ".join(
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip() and not line.strip().startswith("//")
    )
    sector_match = re.search(
        r"SECTOR_SOURCE ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)",
        compact,
    )
    overlap_match = re.search(
        r"QUADRATIC_OVERLAP "
        r"([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)",
        compact,
    )
    assert sector_match is not None, result.stdout
    assert overlap_match is not None, result.stdout
    sector_values = tuple(map(int, sector_match.groups()))
    overlap_values = tuple(map(int, overlap_match.groups()))
    sector_dimension = sector_values[2] - sector_values[3]
    overlap_dimension = overlap_values[3] - overlap_values[4]
    return {
        "sector_source_mod_base_square": {
            "artin_rees_cutoff": 4,
            "cutoff_intersection_generators": sector_values[0],
            "cutoff_remainder_generators": sector_values[1],
            "denominator_quotient_length": sector_values[2],
            "numerator_quotient_length": sector_values[3],
            "dimension": sector_dimension,
        },
        "quadratic_overlap_mod_base_square": {
            "numerator_generators": overlap_values[0],
            "artin_rees_cutoff": 5,
            "cutoff_intersection_generators": overlap_values[1],
            "cutoff_remainder_generators": overlap_values[2],
            "denominator_quotient_length": overlap_values[3],
            "numerator_quotient_length": overlap_values[4],
            "dimension": overlap_dimension,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-order", type=int, default=2)
    parser.add_argument("--normal-order", type=int, default=2)
    arguments = parser.parse_args()
    print(audit(arguments.base_order, arguments.normal_order))


if __name__ == "__main__":
    main()
