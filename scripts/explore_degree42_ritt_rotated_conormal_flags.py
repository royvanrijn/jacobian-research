#!/usr/bin/env python3
"""Explore the two uncomputed degree-42 completed braid sectors.

The certified chart ``273`` represents the sector omitting composite cut
``6``.  The opposite-pair representatives ``237`` and ``327`` omit
respectively composite cuts ``14`` and ``21`` on one half-braid.  This
script constructs either chart directly, derives its two paths from the
filled Ritt braid, changes to seven graph-normal and two Dickson-base
coordinates, and asks Singular for the first local flag invariants.

This is exploratory.  In particular, a finite jet or an annihilation result
is not promoted to a completed statement until a separate Artin--Rees or
Nakayama cutoff is proved.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import (  # noqa: E402
    canonical_reconstruction,
)
from jcsearch.ritt_complex import (  # noqa: E402
    MoveType,
    compose_factors,
    dickson_vertex_factors,
    symmetric_braid_complex,
)


W = sp.Symbol("W")
DEGREE = 42
ALL_CUTS = frozenset((2, 3, 6, 7, 14, 21))
PRIME_CUTS = frozenset((2, 3, 7))
COMPOSITE_CUTS = frozenset((6, 14, 21))


def serialize_ideal(equations: list[sp.Expr] | list[str]) -> str:
    """Serialize a nonempty SymPy ideal for Singular."""

    assert equations
    return ",".join(
        str(equation).replace("**", "^")
        for equation in equations
    )


def source_ideal_cache(word: tuple[int, int, int]) -> Path:
    """Return the one-time cache for a rotated chart's ordinary equations."""

    label = "".join(map(str, word))
    return (
        ROOT
        / "artifacts"
        / "generated-results"
        / f"degree42_ritt_rotated_source_ideals_{label}.json.gz"
    )


def chart_coordinates(
    word: tuple[int, int, int],
) -> tuple[tuple[sp.Symbol, ...], tuple[tuple[sp.Symbol, ...], ...]]:
    """Reconstruct the stable parameter order without composing the chart."""

    label = "".join(map(str, word))
    factor_variables = tuple(
        sp.symbols(f"x{label}_{position}_1:{degree}")
        for position, degree in enumerate(word)
    )
    parameters = tuple(
        variable
        for variables in factor_variables
        for variable in variables
    )
    return parameters, factor_variables


def build_chart(
    word: tuple[int, int, int],
) -> tuple[
    tuple[sp.Symbol, ...],
    tuple[tuple[sp.Symbol, ...], ...],
    sp.Expr,
]:
    """Return generic monic-original factors and their composition."""

    parameters, stable_factor_variables = chart_coordinates(word)
    factor_variables = []
    factors = []
    for variables, degree in zip(stable_factor_variables, word):
        factor_variables.append(variables)
        factors.append(
            W**degree
            + sum(
                variables[power - 1] * W**power
                for power in range(1, degree)
            )
        )
    return (
        parameters,
        tuple(factor_variables),
        compose_factors(tuple(factors), W),
    )


def rotated_source_ideal_data(
    word: tuple[int, int, int],
) -> tuple[
    tuple[sp.Symbol, ...],
    tuple[tuple[sp.Symbol, ...], ...],
    list[sp.Expr] | list[str],
    list[sp.Expr] | list[str],
    list[sp.Expr] | list[str],
    int,
    int,
]:
    """Load or construct the ordinary thick, thin, and boundary equations."""

    (
        base_cuts,
        thick_cuts,
        thick_omission,
        thin_cuts,
        thin_omission,
    ) = path_flag(word)
    cache = source_ideal_cache(word)
    parameters, factor_variables = chart_coordinates(word)
    if cache.is_file():
        with gzip.open(cache, "rt") as source:
            cached = json.load(source)
        assert cached["parameters"] == [str(parameter) for parameter in parameters]
        assert tuple(cached["word"]) == word
        return (
            parameters,
            factor_variables,
            cached["thick"],
            cached["thin"],
            cached["boundary"],
            thick_omission,
            thin_omission,
        )

    parameters, factor_variables, polynomial = build_chart(word)
    requested_cuts = tuple(sorted(ALL_CUTS - base_cuts))
    residuals = {}
    for cut in requested_cuts:
        print("BUILD_RESIDUAL", cut, flush=True)
        residuals[cut] = raw_canonical_residuals(
            polynomial,
            cut,
            DEGREE // cut,
            minimum_coefficient_degree=1,
        )

    def equations_for(cuts: frozenset[int]) -> list[sp.Expr]:
        return [
            equation
            for cut in sorted(cuts - base_cuts)
            for equation in residuals[cut]
        ]

    thick = equations_for(thick_cuts)
    thin = equations_for(thin_cuts)
    boundary = equations_for(ALL_CUTS)
    cached = {
        "schema": "degree42-ritt-rotated-source-ideals.v1",
        "word": word,
        "parameters": [str(parameter) for parameter in parameters],
        "thick_composite_omission": thick_omission,
        "thin_prime_omission": thin_omission,
        "thick": [str(equation).replace("**", "^") for equation in thick],
        "thin": [str(equation).replace("**", "^") for equation in thin],
        "boundary": [
            str(equation).replace("**", "^") for equation in boundary
        ],
        "construction": (
            "raw canonical residual numerators including coefficient degree "
            "one; nonzero rational rescaling does not change the ideals"
        ),
    }
    with gzip.open(cache, "wt", compresslevel=9) as target:
        json.dump(cached, target, indent=2)
        target.write("\n")
    print(f"WROTE_SOURCE_CACHE {cache.relative_to(ROOT)}", flush=True)
    return (
        parameters,
        factor_variables,
        thick,
        thin,
        boundary,
        thick_omission,
        thin_omission,
    )


def raw_canonical_residuals(
    polynomial: sp.Expr,
    a: int,
    b: int,
    *,
    minimum_coefficient_degree: int = 1,
) -> list[sp.Expr]:
    """Return canonical residual numerators without global primitive parts.

    Taking ``primitive(Poly(...))`` is useful for canonical serialized
    certificates but dominates runtime on the rotated degree-42 charts.
    Multiplication by a nonzero rational scalar does not change the ideal,
    so the local exploration can retain the raw numerator.
    """

    degree = a * b
    source, composition, reconstruction, used_degrees = (
        canonical_reconstruction(polynomial, a, b)
    )
    residuals = []
    for coefficient_degree in range(minimum_coefficient_degree, degree):
        if coefficient_degree in used_degrees:
            continue
        residual = sp.together(
            composition.nth(coefficient_degree).subs(reconstruction)
            - source.nth(coefficient_degree)
        ).as_numer_denom()[0]
        if residual != 0:
            residuals.append(residual)
    return residuals


def path_flag(
    word: tuple[int, int, int],
) -> tuple[
    frozenset[int],
    frozenset[int],
    int,
    frozenset[int],
    int,
]:
    """Return base, thick-path, thick omission, thin-path, thin omission."""

    braid = symmetric_braid_complex((2, 3, 7), MoveType.CHEBYSHEV)
    vertex_by_word = {vertex.word: vertex for vertex in braid.vertices}
    adjacency = {vertex.word: [] for vertex in braid.vertices}
    for edge in braid.edges:
        left, right = edge.endpoints
        adjacency[left].append(right)
        adjacency[right].append(left)
    endpoint = tuple(reversed(word))
    paths: list[tuple[tuple[int, int, int], ...]] = []

    def extend(path: tuple[tuple[int, int, int], ...]) -> None:
        if len(path) == 4:
            if path[-1] == endpoint:
                paths.append(path)
            return
        for neighbor in adjacency[path[-1]]:
            if neighbor not in path:
                extend(path + (neighbor,))

    extend((word,))
    assert len(paths) == 2
    base_cuts = frozenset(vertex_by_word[word].cuts)
    data = []
    for path in paths:
        cuts = frozenset().union(
            *(frozenset(vertex_by_word[vertex].cuts) for vertex in path)
        )
        omitted = ALL_CUTS - cuts
        assert len(omitted) == 1
        data.append((cuts, next(iter(omitted))))
    thick_cuts, thick_omission = next(
        item for item in data if item[1] in COMPOSITE_CUTS
    )
    thin_cuts, thin_omission = next(
        item for item in data if item[1] in PRIME_CUTS
    )
    return (
        base_cuts,
        thick_cuts,
        thick_omission,
        thin_cuts,
        thin_omission,
    )


def graph_normal_map(
    word: tuple[int, int, int],
    factor_variables: tuple[tuple[sp.Symbol, ...], ...],
) -> tuple[
    tuple[sp.Symbol, ...],
    tuple[sp.Symbol, sp.Symbol],
    dict[sp.Symbol, sp.Expr],
]:
    """Return the exact ``7 normal | 2 base`` Dickson chart map."""

    inner_degree = word[-1]
    assert inner_degree >= 3
    tau, zeta = sp.symbols("tau zeta")
    translation = 1 + tau
    factors = dickson_vertex_factors(word, W, translation, zeta)
    inner_top = factor_variables[-1][-1]
    inner_next = factor_variables[-1][-2]
    free_variables = {inner_top, inner_next}
    dependent_variables = tuple(
        variable
        for variables in factor_variables
        for variable in variables
        if variable not in free_variables
    )
    normals = sp.symbols(f"n0:{len(dependent_variables)}")
    normal_by_variable = dict(zip(dependent_variables, normals))
    images = {}
    for factor, variables in zip(factors, factor_variables):
        polynomial = sp.Poly(factor, W)
        for power, variable in enumerate(variables, 1):
            images[variable] = polynomial.nth(power)
            if variable in normal_by_variable:
                images[variable] += normal_by_variable[variable]
    assert sp.expand(images[inner_top] - inner_degree * translation) == 0
    assert sp.expand(
        images[inner_next]
        - (
            comb(inner_degree, 2) * translation**2
            - inner_degree * zeta
        )
    ) == 0
    return normals, (tau, zeta), images


def singular_local_audit(
    word: tuple[int, int, int],
    parameters: tuple[sp.Symbol, ...],
    factor_variables: tuple[tuple[sp.Symbol, ...], ...],
    thick: list[sp.Expr],
    thin: list[sp.Expr],
    boundary: list[sp.Expr],
) -> str:
    """Run the local flag audit after the graph-normal coordinate change."""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    normals, base_coordinates, images = graph_normal_map(
        word, factor_variables
    )
    local_variables = normals + base_coordinates
    map_images = ",".join(
        str(sp.expand(images[parameter])).replace("**", "^")
        for parameter in parameters
    )
    zero_substitutions = "\n".join(
        f"JT=subst(JT,{variable},0);"
        f"JTHIN=subst(JTHIN,{variable},0);"
        f"JB=subst(JB,{variable},0);"
        f"JK=subst(JK,{variable},0);"
        for variable in local_variables
    )
    zeta = base_coordinates[1]
    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={serialize_ideal(thick)};
ideal ITHINsource={serialize_ideal(thin)};
ideal IBsource={serialize_ideal(boundary)};
ring q=0,({",".join(map(str, local_variables))}),(dp({len(normals)}),dp(2));
map phi=source,{map_images};
option(redSB);
proc iszeroideal(ideal J)
{{
  ideal reducedGenerators=simplify(J,2);
  return(size(reducedGenerators)==0);
}}
proc issubsetstd(ideal A, ideal G)
{{
  ideal remainder=reduce(A,G);
  return(iszeroideal(remainder));
}}
proc firstannihilatingpower(ideal A, ideal G, poly z, int maximumPower)
{{
  ideal current=A;
  int exponent;
  for(exponent=0; exponent<=maximumPower; exponent++)
  {{
    if(issubsetstd(current,G))
    {{
      return(exponent);
    }}
    current=z*current;
  }}
  return(-1);
}}
ideal IT=phi(ITsource);
ideal ITHIN=phi(ITHINsource);
ideal IB=phi(IBsource);
ideal K={",".join(map(str, normals))};
ideal maximalIdeal={",".join(map(str, local_variables))};
ideal maximalIdeal2=maximalIdeal*maximalIdeal;
ideal maximalIdeal3=maximalIdeal2*maximalIdeal;
ideal maximalIdeal4=maximalIdeal3*maximalIdeal;
ideal GIT=std(IT+maximalIdeal4);
ideal GTHIN=std(ITHIN+maximalIdeal4);
ideal GIB=std(IB+maximalIdeal4);
ideal GK=std(K+maximalIdeal4);
matrix JT=jacob(IT);
matrix JTHIN=jacob(ITHIN);
matrix JB=jacob(IB);
matrix JK=jacob(K);
{zero_substitutions}
print("IDEAL_FLAG_Q4");
print(issubsetstd(IT,GIB));
print(issubsetstd(ITHIN,GIB));
print(issubsetstd(IB,GTHIN));
print(issubsetstd(IB,GK));
print("CONORMAL_RANKS");
print(rank(JT));
print(rank(JTHIN));
print(rank(JB));
print(rank(JK));
print("BASE_ANNIHILATION_EXPONENTS_Q4");
print(firstannihilatingpower(IB,GIT,{zeta},50));
print(firstannihilatingpower(K,GIB,{zeta},50));
print("JET_LENGTHS_Q1");
print(vdim(std(IT+maximalIdeal)));
print(vdim(std(ITHIN+maximalIdeal)));
print(vdim(std(IB+maximalIdeal)));
print(vdim(std(K+maximalIdeal)));
print("JET_LENGTHS_Q2");
print(vdim(std(IT+maximalIdeal2)));
print(vdim(std(ITHIN+maximalIdeal2)));
print(vdim(std(IB+maximalIdeal2)));
print(vdim(std(K+maximalIdeal2)));
print("JET_LENGTHS_Q3");
print(vdim(std(IT+maximalIdeal3)));
print(vdim(std(ITHIN+maximalIdeal3)));
print(vdim(std(IB+maximalIdeal3)));
print(vdim(std(K+maximalIdeal3)));
print("JET_LENGTHS_Q4");
print(vdim(GIT));
print(vdim(GTHIN));
print(vdim(GIB));
print(vdim(GK));
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=7200,
    )
    return result.stdout


def record_exact_output(
    word: tuple[int, int, int],
    thick_omission: int,
    thin_omission: int,
    equation_counts: tuple[int, int, int],
    output: str,
) -> None:
    """Assert and record the exact fourth-jet output."""

    tokens = output.split()

    def values_after(label: str, count: int) -> tuple[int, ...]:
        start = tokens.index(label) + 1
        return tuple(int(value) for value in tokens[start : start + count])

    observed = {
        "ideal_flag_q4": values_after("IDEAL_FLAG_Q4", 4),
        "conormal_ranks": values_after("CONORMAL_RANKS", 4),
        "base_annihilation_exponents_q4": values_after(
            "BASE_ANNIHILATION_EXPONENTS_Q4", 2
        ),
        "jet_lengths_q1": values_after("JET_LENGTHS_Q1", 4),
        "jet_lengths_q2": values_after("JET_LENGTHS_Q2", 4),
        "jet_lengths_q3": values_after("JET_LENGTHS_Q3", 4),
        "jet_lengths_q4": values_after("JET_LENGTHS_Q4", 4),
    }
    expected_by_word = {
        (2, 3, 7): {
            "ideal_flag_q4": (1, 1, 1, 1),
            "conormal_ranks": (5, 6, 6, 7),
            "base_annihilation_exponents_q4": (3, 1),
            "jet_lengths_q1": (1, 1, 1, 1),
            "jet_lengths_q2": (5, 4, 4, 3),
            "jet_lengths_q3": (13, 9, 9, 6),
            "jet_lengths_q4": (25, 16, 16, 10),
        },
        (3, 2, 7): {
            "ideal_flag_q4": (1, 1, 1, 1),
            "conormal_ranks": (5, 6, 6, 7),
            "base_annihilation_exponents_q4": (3, 1),
            "jet_lengths_q1": (1, 1, 1, 1),
            "jet_lengths_q2": (5, 4, 4, 3),
            "jet_lengths_q3": (13, 9, 9, 6),
            "jet_lengths_q4": (26, 16, 16, 10),
        },
    }
    assert observed == expected_by_word[word]
    label = "".join(map(str, word))
    artifact = (
        ROOT
        / "artifacts"
        / "generated-results"
        / f"degree42_ritt_rotated_conormal_jet_{label}.json"
    )
    data = {
        "schema": "degree42-ritt-rotated-conormal-jet.v1",
        "status": "exact fourth-maximal-adic-jet computation",
        "word": word,
        "opposite_word": tuple(reversed(word)),
        "thick_composite_omission": thick_omission,
        "thin_prime_omission": thin_omission,
        "equation_counts": {
            "thick": equation_counts[0],
            "thin": equation_counts[1],
            "boundary": equation_counts[2],
        },
        **observed,
        "sector_layer_dimensions_q2_q3_q4": tuple(
            observed[f"jet_lengths_q{order}"][0]
            - observed[f"jet_lengths_q{order}"][2]
            for order in (2, 3, 4)
        ),
        "spectator_layer_dimensions_q2_q3_q4": tuple(
            observed[f"jet_lengths_q{order}"][2]
            - observed[f"jet_lengths_q{order}"][3]
            for order in (2, 3, 4)
        ),
        "theorem_boundary": (
            "These are exact quotient calculations modulo the fourth "
            "maximal-ideal power. They do not prove equality, annihilation, "
            "or non-splitting in the completed local ring."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/explore_degree42_ritt_rotated_conormal_flags.py "
            f"--word {label}"
        ),
    }
    artifact.write_text(json.dumps(data, indent=2) + "\n")
    print(f"PASS: exact fourth jet matches the pinned {label} profile")
    print(f"PASS: wrote {artifact.relative_to(ROOT)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--word",
        choices=("237", "327"),
        required=True,
        help="opposite-pair representative to audit",
    )
    parser.add_argument(
        "--build-source-only",
        action="store_true",
        help="write/reuse the ordinary-equation cache and stop before Singular",
    )
    parser.add_argument(
        "--rebuild-source",
        action="store_true",
        help="replace the selected compressed source cache from residuals",
    )
    arguments = parser.parse_args()
    word = tuple(int(character) for character in arguments.word)
    cache = source_ideal_cache(word)
    if arguments.rebuild_source and cache.is_file():
        cache.unlink()
    (
        base_cuts,
        thick_cuts,
        thick_omission,
        thin_cuts,
        thin_omission,
    ) = path_flag(word)
    print(
        "PATH_DATA",
        word,
        "base",
        sorted(base_cuts),
        "thick_omission",
        thick_omission,
        "thin_omission",
        thin_omission,
        flush=True,
    )
    (
        parameters,
        factor_variables,
        thick,
        thin,
        boundary,
        cached_thick_omission,
        cached_thin_omission,
    ) = rotated_source_ideal_data(word)
    assert cached_thick_omission == thick_omission
    assert cached_thin_omission == thin_omission
    equation_counts = (len(thick), len(thin), len(boundary))
    print(
        "EQUATION_COUNTS",
        *equation_counts,
        flush=True,
    )
    if arguments.build_source_only:
        print("PASS: source ideals are cached")
        return
    output = singular_local_audit(
        word,
        parameters,
        factor_variables,
        thick,
        thin,
        boundary,
    )
    print(output)
    record_exact_output(
        word,
        thick_omission,
        thin_omission,
        equation_counts,
        output,
    )


if __name__ == "__main__":
    main()
