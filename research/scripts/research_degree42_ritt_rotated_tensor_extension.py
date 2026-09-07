#!/usr/bin/env python3
"""Research rotated degree-42 conormal extensions by tensor presentations.

For a selected chart, this presents

    K/I_thick -> K/I_boundary

as modules before tensoring with

    R / ((tau,zeta)^a + (n0,...,n6)^b).

Thus a splitting of the completed module projection would induce a
splitting in every finite calculation below.  With ``b=1`` this is the
base-changed first-conormal projection, not merely a quotient-ring image.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree42_ritt_rotated_conormal_flags import (  # noqa: E402
    graph_normal_map,
    rotated_source_ideal_data,
    serialize_ideal,
)
from research_degree42_cellular_extension import (  # noqa: E402
    vector_space_section,
    splitting_solution,
)
from research_degree42_tensor_extension import (  # noqa: E402
    module_coordinate_matrix,
    parse_sections,
)


def matrix_cache(
    word: tuple[int, int, int],
    base_order: int,
    normal_order: int,
) -> Path:
    """Return the compressed cache for Singular bases and reductions."""

    label = "".join(map(str, word))
    return (
        ROOT
        / "artifacts"
        / "generated-results"
        / (
            f"degree42_ritt_rotated_tensor_matrices_{label}"
            f"_b{base_order}_n{normal_order}.txt.gz"
        )
    )


def exact_rank_certificate(
    inclusion: sp.Matrix,
    projection: sp.Matrix,
    sector_actions: tuple[sp.Matrix, ...],
    total_actions: tuple[sp.Matrix, ...],
    spectator_actions: tuple[sp.Matrix, ...],
    action_names: tuple[str, ...],
) -> dict[str, object]:
    """Return the adapted cocycle with fraction-free exact rank tests."""

    section = vector_space_section(projection)
    adapted_basis = inclusion.row_join(section)
    adapted_inverse = adapted_basis.inv()
    sector_dimension = inclusion.cols
    spectator_dimension = section.cols
    couplings = []
    for sector_action, total_action, spectator_action in zip(
        sector_actions, total_actions, spectator_actions
    ):
        adapted_action = adapted_inverse * total_action * adapted_basis
        assert (
            adapted_action[:sector_dimension, :sector_dimension]
            == sector_action
        )
        assert (
            adapted_action[sector_dimension:, sector_dimension:]
            == spectator_action
        )
        assert (
            adapted_action[sector_dimension:, :sector_dimension]
            == sp.zeros(spectator_dimension, sector_dimension)
        )
        couplings.append(
            adapted_action[:sector_dimension, sector_dimension:]
        )

    h_variables = sp.symbols(
        f"h0:{sector_dimension * spectator_dimension}"
    )
    correction = sp.Matrix(
        sector_dimension, spectator_dimension, h_variables
    )
    expressions = []
    for sector_action, spectator_action in zip(
        sector_actions, spectator_actions
    ):
        expressions.extend(
            sector_action * correction - correction * spectator_action
        )
    coboundary, _ = sp.linear_eq_to_matrix(expressions, h_variables)
    cocycle = sp.Matrix(
        [value for coupling in couplings for value in coupling]
    )
    coboundary_rank = coboundary.to_DM().rank()
    augmented_rank = coboundary.row_join(-cocycle).to_DM().rank()
    return {
        "adapted_coupling": {
            name: coupling.tolist()
            for name, coupling in zip(action_names, couplings)
        },
        "coboundary_rank": coboundary_rank,
        "augmented_rank": augmented_rank,
        "obstruction_detected": augmented_rank > coboundary_rank,
    }


def rotated_singular_output(
    word: tuple[int, int, int],
    base_order: int,
    normal_order: int,
) -> tuple[str, tuple[sp.Symbol, ...]]:
    """Return bases, projection, and coordinate actions from Singular."""

    (
        parameters,
        factor_variables,
        thick,
        _thin,
        boundary,
        _thick_omission,
        _thin_omission,
    ) = rotated_source_ideal_data(word)
    normals, base_coordinates, images = graph_normal_map(
        word, factor_variables
    )
    local_variables = normals + base_coordinates
    cache = matrix_cache(word, base_order, normal_order)
    if cache.is_file():
        with gzip.open(cache, "rt") as source:
            return source.read(), local_variables
    map_images = ",".join(
        str(images[parameter]).replace("**", "^")
        for parameter in parameters
    )
    normal_module = ",".join(f"[{normal}]" for normal in normals)
    action_blocks = []
    for target, basis, presentation in (
        ("TOTAL", "BT", "PT"),
        ("SPECTATOR", "BS", "PS"),
    ):
        for variable in local_variables:
            loop = f"index_{target.lower()}_{variable}"
            action_blocks.append(
                f"""
print("ACTION_{target}_{variable}");
for (int {loop}=1; {loop}<=ncols({basis}); {loop}++)
{{
  print(reduce(({variable})*{basis}[{loop}],{presentation}));
}}
"""
            )

    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={serialize_ideal(thick)};
ideal IBsource={serialize_ideal(boundary)};
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
    assert singular is not None, "Singular is required"
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=7200,
    )
    required_markers = (
        "BASIS_TOTAL",
        "BASIS_SPECTATOR",
        "PROJECTION",
    )
    missing_markers = tuple(
        marker for marker in required_markers if marker not in result.stdout
    )
    if missing_markers:
        raise RuntimeError(
            "Singular returned an incomplete rotated tensor transcript; "
            f"missing {missing_markers}. stderr: {result.stderr.strip()}"
        )
    with gzip.open(cache, "wt", compresslevel=9) as target:
        target.write(result.stdout)
    return result.stdout, local_variables


def rotated_tensor_audit(
    word: tuple[int, int, int],
    base_order: int = 2,
    normal_order: int = 1,
) -> dict[str, object]:
    """Return the finite module extension and its exact splitting test."""

    output, symbols = rotated_singular_output(
        word, base_order, normal_order
    )
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
    module_rank = 7
    projection = module_coordinate_matrix(
        sections["PROJECTION"], spectator_basis, symbols, module_rank
    )
    total_actions = tuple(
        module_coordinate_matrix(
            sections[f"ACTION_TOTAL_{variable}"],
            total_basis,
            symbols,
            module_rank,
        )
        for variable in symbols
    )
    spectator_actions = tuple(
        module_coordinate_matrix(
            sections[f"ACTION_SPECTATOR_{variable}"],
            spectator_basis,
            symbols,
            module_rank,
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
    certificate = exact_rank_certificate(
        kernel,
        projection,
        kernel_actions,
        total_actions,
        spectator_actions,
        tuple(map(str, symbols)),
    )

    def serialized_matrix(matrix: sp.Matrix) -> list[list[str]]:
        return [[str(entry) for entry in row] for row in matrix.tolist()]

    (
        _parameters,
        _factor_variables,
        _thick,
        _thin,
        _boundary,
        thick_omission,
        thin_omission,
    ) = rotated_source_ideal_data(word)
    return {
        "word": word,
        "thick_composite_omission": thick_omission,
        "thin_prime_omission": thin_omission,
        "base_order": base_order,
        "normal_order": normal_order,
        "dimensions": {
            "kernel": kernel.cols,
            "total": len(total_basis),
            "spectator": len(spectator_basis),
        },
        "splits_as_R_module": splits,
        "one_section": (
            serialized_matrix(section) if section is not None else None
        ),
        "zero_action_variables": [
            str(symbol)
            for symbol, total_action, spectator_action in zip(
                symbols, total_actions, spectator_actions
            )
            if total_action.is_zero_matrix
            and spectator_action.is_zero_matrix
        ],
        "finite_module_presentation": {
            "variables": [str(symbol) for symbol in symbols],
            "projection": serialized_matrix(projection),
            "total_actions": [
                serialized_matrix(action) for action in total_actions
            ],
            "spectator_actions": [
                serialized_matrix(action) for action in spectator_actions
            ],
        },
        "certificate": certificate,
        "interpretation": {
            "total": (
                "(K/(I_thick+K^2)) tensor_B "
                "B/(tau,zeta)^base_order"
            ),
            "spectator": (
                "(K/(I_boundary+K^2)) tensor_B "
                "B/(tau,zeta)^base_order"
            ),
            "base_ring": "B=Q[[tau,zeta]]=completed R/K",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--word", choices=("237", "327"), required=True)
    parser.add_argument("--base-order", type=int, default=2)
    parser.add_argument("--normal-order", type=int, default=1)
    arguments = parser.parse_args()
    word = tuple(int(character) for character in arguments.word)
    print(
        rotated_tensor_audit(
            word,
            base_order=arguments.base_order,
            normal_order=arguments.normal_order,
        )
    )


if __name__ == "__main__":
    main()
