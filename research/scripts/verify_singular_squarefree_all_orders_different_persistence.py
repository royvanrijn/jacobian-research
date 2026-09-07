#!/usr/bin/env python3
"""All-orders different persistence for all squarefree ternary cubics.

For each singular squarefree ternary-cubic symbol, this checker chooses a
minimal set of degree-three compatible tensors whose classes generate the
exact graded gauge cokernel.  It verifies the full quotient presentation,
then studies the universal coefficient family

    h + u_1 eta_1 + ... + u_r eta_r.

The intrinsic collision Nakayama module of the Kahler different is the
constant module Q[u_1,...,u_r]^6.  With every u_i assigned collision weight
one, strict Rees packets for Omega and coker(B -> Omega^3) have their literal
central initial presentations.

Successive formal gauge elimination writes every compatible higher tail as
sum f_i eta_i.  The graph equations u_i-f_i form a filtered regular sequence:
at each stage their initial forms are monic in the next u_i.  Thus the
annihilator commutes with arbitrary polynomial and formal graph
specialization, and the specialized Kahler different has six minimal local
generators.  The smooth row has no normal coefficients: the independent
formal-rigidity theorem reduces it to the six-generator central calculation.
Normality and Keller-open compatibility are not asserted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import verify_cubic_formal_gauge_cokernel_atlas as atlas  # noqa: E402
import verify_cubic_symbol_double_saturation as cubic  # noqa: E402
import verify_universal_cubic_cotangent_saturation as smooth  # noqa: E402


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "singular_squarefree_all_orders_different_persistence.json"
)

SOURCES: dict[str, tuple[sp.Expr, ...]] = {
    "nodal": (cubic.Z**3,),
    "cuspidal": (cubic.Z**3, cubic.X * cubic.Z**2),
    "line-transverse-conic": (cubic.Y**3, cubic.X**3),
    "line-tangent-conic": (
        cubic.Y**2 * cubic.Z,
        cubic.Y**3,
        cubic.X * cubic.Y**2,
    ),
    "triangle": (cubic.Z**3, cubic.Y**3, cubic.X**3),
    "concurrent-lines": (
        cubic.Z**3,
        cubic.Y * cubic.Z**2,
        cubic.Y**2 * cubic.Z,
        cubic.X * cubic.Z**2,
    ),
}

# Columns of the exact pruned presentation of ker(C)/im(G_h), after the
# quotient basis changes performed by exact pruning.  The source tensors
# below are independently verified to generate the same quotient.
EXPECTED_GAUGE_PRESENTATIONS: dict[str, list[list[sp.Expr]]] = {
    "nodal": [[cubic.x]],
    "cuspidal": [
        [3 * cubic.x, -cubic.z],
        [0, cubic.x],
    ],
    "line-transverse-conic": [
        [0, cubic.z],
        [cubic.y, 0],
    ],
    "line-tangent-conic": [
        [0, cubic.y, cubic.z],
        [0, 0, cubic.y],
        [-6 * cubic.y, -2 * cubic.z, cubic.x],
        [3 * cubic.y**2, 0, -cubic.z**2],
    ],
    "triangle": [
        [0, 0, cubic.z],
        [0, cubic.y, 0],
        [cubic.x, 0, 0],
    ],
    "concurrent-lines": [
        [3 * cubic.x, -cubic.y, -cubic.z, 0],
        [0, 2 * cubic.x, 0, -2 * cubic.y + cubic.z],
        [0, 0, 2 * cubic.x, cubic.y - 2 * cubic.z],
        [0, 0, 0, cubic.x],
    ],
}

EXPECTED_ANNIHILATORS = {
    "nodal": "x",
    "cuspidal": "x2",
    "line-transverse-conic": "yz",
    "line-tangent-conic": "y3",
    "triangle": "xyz",
    "concurrent-lines": "x3",
}


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def tensor_dictionary(source: sp.Expr) -> dict[tuple[int, int, int], sp.Expr]:
    tensor = atlas.symbol_tensor(source)
    return {
        triple: sp.expand(tensor[row])
        for row, triple in enumerate(smooth.TRIPLES)
    }


def gauge_quotient_certificate(
    name: str,
    tensors: tuple[dict[tuple[int, int, int], sp.Expr], ...],
) -> dict[str, Any]:
    """Verify generation and the exact minimal gauge-cokernel matrix."""

    compatibility = smooth.compatibility_matrix()
    leading_tensor = atlas.symbol_tensor(cubic.CUBIC_STRATA[name])
    gauge = smooth.gauge_matrix(leading_tensor)
    eta_modules = ",".join(
        atlas.singular_module(
            sp.Matrix([tensor[triple] for triple in smooth.TRIPLES])
        )
        for tensor in tensors
    )
    expected = ",".join(
        cubic.singular_vector(column)
        for column in EXPECTED_GAUGE_PRESENTATIONS[name]
    )
    program = f"""
ring quotient_ring=0,(x,y,z),dp;
module C={atlas.singular_module(compatibility)};
module G={atlas.singular_module(gauge)};
module K=syz(C);
module ETAS={eta_modules};
module H=G,ETAS;
H=std(H);
module generation_difference=simplify(reduce(K,H),2);
module Q=std(prune(std(modulo(K,G))));
module expected={expected};
expected=std(expected);
module quotient_minus_expected=simplify(reduce(Q,expected),2);
module expected_minus_quotient=simplify(reduce(expected,Q),2);
ideal ANN=std(quotient(G,K));
print("GENERATION_DIFFERENCE="+string(size(generation_difference)));
print("PRESENTATION_DIFFERENCE="+string(
  size(quotient_minus_expected)+size(expected_minus_quotient)
));
print("PRUNED_RANK="+string(nrows(Q)));
print("ANNIHILATOR="+string(ANN));
quit;
"""
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=1800,
    )
    values: dict[str, str] = {}
    for line in completed.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in {
                "GENERATION_DIFFERENCE",
                "PRESENTATION_DIFFERENCE",
                "PRUNED_RANK",
                "ANNIHILATOR",
            }:
                values[key] = value
    assert set(values) == {
        "GENERATION_DIFFERENCE",
        "PRESENTATION_DIFFERENCE",
        "PRUNED_RANK",
        "ANNIHILATOR",
    }, completed.stdout + completed.stderr
    result: dict[str, Any] = {
        "generation_difference": int(values["GENERATION_DIFFERENCE"]),
        "presentation_difference": int(values["PRESENTATION_DIFFERENCE"]),
        "minimal_generator_count": int(values["PRUNED_RANK"]),
        "annihilator_singular_syntax": values["ANNIHILATOR"],
    }
    assert result == {
        "generation_difference": 0,
        "presentation_difference": 0,
        "minimal_generator_count": len(tensors),
        "annihilator_singular_syntax": EXPECTED_ANNIHILATORS[name],
    }, (name, result)
    return result


def expected_nakayama_presentation_difference(
    name: str,
    tensors: tuple[dict[tuple[int, int, int], sp.Expr], ...],
) -> int:
    """Compare the universal pruned J/nJ presentation with (x,y,z)^6."""

    program, _parameters = cubic.singular_subspace_program(
        cubic.CUBIC_STRATA[name],
        tensors,
    )
    expected_generators: list[list[sp.Expr]] = []
    for component in range(6):
        for variable in cubic.BASE_VARIABLES:
            generator = [sp.Integer(0)] * 6
            generator[component] = variable
            expected_generators.append(generator)
    program = program.replace(
        "quit;",
        "module expected_different_generator_presentation="
        + ",".join(map(cubic.singular_vector, expected_generators))
        + ";expected_different_generator_presentation=std("
        + "expected_different_generator_presentation);"
        + "module universal_minus_expected=simplify(reduce("
        + "pruned_different_generator_presentation,"
        + "expected_different_generator_presentation),2);"
        + "module expected_minus_universal=simplify(reduce("
        + "expected_different_generator_presentation,"
        + "pruned_different_generator_presentation),2);"
        + 'print("EXPECTED_DIFFERENT_PRESENTATION_DIFFERENCE="'
        + "+string(size(universal_minus_expected)"
        + "+size(expected_minus_universal)));"
        + "quit;",
    )
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=1800,
    )
    prefix = "EXPECTED_DIFFERENT_PRESENTATION_DIFFERENCE="
    values = [
        int(line[len(prefix) :])
        for line in completed.stdout.splitlines()
        if line.startswith(prefix)
    ]
    assert len(values) == 1, completed.stdout + completed.stderr
    return values[0]


def expected_universal_family(
    parameter_count: int,
    ext_difference: int,
) -> dict[str, int]:
    return {
        "parameter_count": parameter_count,
        "cotangent_saturation_generators": 0,
        "support_module_dimension": parameter_count + 2,
        "support_ext3_vector_dimension": 0,
        "support_ext2_dimension": parameter_count,
        "support_ext2_multiplicity": 6,
        "support_ext2_parameter_axis_radical_difference": 0,
        "support_ext2_central_pruned_presentation_difference": ext_difference,
        "support_ext2_pruned_presentation_rank": 3,
        "support_ext2_collision_square_action_generators": 0,
        "different_generator_module_dimension": parameter_count,
        "different_generator_module_multiplicity": 6,
        "different_generator_parameter_axis_radical_difference": 0,
        "different_generator_central_pruned_presentation_difference": 0,
        "different_generator_pruned_presentation_rank": 6,
    }


def smooth_central_certificate() -> dict[str, int]:
    """Verify the zero-parameter smooth row used with formal rigidity."""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    completed = subprocess.run(
        [singular, "-q"],
        input=cubic.singular_program(cubic.CUBIC_STRATA["smooth"]),
        text=True,
        capture_output=True,
        check=True,
        timeout=1800,
    )
    wanted = {
        "SATURATION_GENERATORS",
        "SUPPORT_DIMENSION",
        "DIFFERENT_GENERATOR_DIMENSION",
        "DIFFERENT_GENERATOR_MULTIPLICITY",
        "DIFFERENT_MINIMAL_GENERATORS",
        "EXT3_VECTOR_DIMENSION",
        "EXT2_DIMENSION",
        "EXT2_VECTOR_DIMENSION",
        "EXT2_TOP_DIMENSION",
        "EXT2_SQUARE_ACTION_GENERATORS",
    }
    values: dict[str, int] = {}
    for line in completed.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in wanted:
                values[key] = int(value)
    assert set(values) == wanted, completed.stdout + completed.stderr
    expected = {
        "SATURATION_GENERATORS": 0,
        "SUPPORT_DIMENSION": 2,
        "DIFFERENT_GENERATOR_DIMENSION": 0,
        "DIFFERENT_GENERATOR_MULTIPLICITY": 6,
        "DIFFERENT_MINIMAL_GENERATORS": 6,
        "EXT3_VECTOR_DIMENSION": 0,
        "EXT2_DIMENSION": 0,
        "EXT2_VECTOR_DIMENSION": 6,
        "EXT2_TOP_DIMENSION": 3,
        "EXT2_SQUARE_ACTION_GENERATORS": 0,
    }
    assert values == expected
    return {
        "cotangent_saturation_generators": values[
            "SATURATION_GENERATORS"
        ],
        "support_dimension": values["SUPPORT_DIMENSION"],
        "different_dimension": values["DIFFERENT_GENERATOR_DIMENSION"],
        "different_multiplicity": values[
            "DIFFERENT_GENERATOR_MULTIPLICITY"
        ],
        "different_minimal_generators": values[
            "DIFFERENT_MINIMAL_GENERATORS"
        ],
        "support_ext2_vector_dimension": values[
            "EXT2_VECTOR_DIMENSION"
        ],
        "support_ext2_top_dimension": values["EXT2_TOP_DIMENSION"],
        "support_ext2_collision_square_action_generators": values[
            "EXT2_SQUARE_ACTION_GENERATORS"
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="intentionally replace the pinned generated artifact",
    )
    args = parser.parse_args()
    cubic.FACTOR_SINGULAR_EXPRESSIONS = False
    compatibility = smooth.compatibility_matrix()
    cases: dict[str, Any] = {}
    smooth_central = smooth_central_certificate()

    for name, sources in SOURCES.items():
        tensors = tuple(tensor_dictionary(source) for source in sources)
        for tensor in tensors:
            column = sp.Matrix(
                [tensor[triple] for triple in smooth.TRIPLES]
            )
            assert (
                compatibility * column
            ).applyfunc(sp.expand) == sp.zeros(6, 1)

        quotient = gauge_quotient_certificate(name, tensors)
        family = cubic.run_singular_subspace_certificate(
            cubic.CUBIC_STRATA[name],
            tensors,
            timeout=1800,
        )
        ext_difference = 4 if name == "line-tangent-conic" else 0
        assert family == expected_universal_family(
            len(tensors), ext_difference
        )
        nakayama_difference = expected_nakayama_presentation_difference(
            name, tensors
        )
        assert nakayama_difference == 0
        rees = cubic.run_singular_rees_base_change_certificate(
            cubic.CUBIC_STRATA[name],
            tensors,
            timeout=1800,
            parameter_collision_weights=(1,) * len(tensors),
        )
        assert rees == {
            "parameter_count": len(tensors),
            "cotangent_rees_torsion_generators": 0,
            "cotangent_initial_presentation_difference": 0,
            "annihilator_cokernel_rees_torsion_generators": 0,
            "annihilator_cokernel_initial_presentation_difference": 0,
        }
        cases[name] = {
            "source_cubics": [sp.sstr(source) for source in sources],
            "tensor_components": [
                [sp.sstr(tensor[triple]) for triple in smooth.TRIPLES]
                for tensor in tensors
            ],
            "gauge_quotient": quotient,
            "gauge_presentation_columns": [
                [sp.sstr(entry) for entry in column]
                for column in EXPECTED_GAUGE_PRESENTATIONS[name]
            ],
            "universal_coefficient_family": family,
            "expected_nakayama_presentation_difference": (
                nakayama_difference
            ),
            "weighted_rees_certificate": rees,
        }

    exact_data = {
        "smooth_formally_rigid_central": smooth_central,
        "singular_cases": cases,
    }
    artifact = {
        "schema": "singular-squarefree-all-orders-different-persistence.v1",
        "case": "seven-complete-formal-squarefree-cubic-slices",
        "mathematical_scope": (
            "Exact characteristic-zero theorem for every compatible formal "
            "tensor with any squarefree ternary-cubic leading symbol. The "
            "smooth row is formally rigid and its central intrinsic "
            "different has six generators. In the six singular rows, "
            "minimal degree-three generators give a "
            "recursive formal normal form sum(f_i*eta_i). On every universal "
            "coefficient family the intrinsic different Nakayama module is "
            "Q[u_1,...,u_r]^6, and weight-one Rees packets are strict with "
            "central associated graded. The graph equations u_i-f_i have "
            "successively monic initial forms, so the intrinsic annihilator "
            "commutes with every polynomial or formal graph specialization. "
            "Every specialized Kahler different is six-generated and "
            "non-Cartier. Normality, algebraization of infinite gauges, and "
            "Keller-open compatibility are not asserted."
        ),
        "formal_graph_lemma": {
            "coefficient_counts": {
                "smooth": 0,
                **{
                    name: len(sources)
                    for name, sources in SOURCES.items()
                },
            },
            "coefficient_weights": 1,
            "normal_form": "h+sum_i f_i*eta_i with ord(f_i)>=1",
            "regular_sequence": "u_1-f_1,...,u_r-f_r",
            "regularity_reason": (
                "after each preceding graph quotient, the next initial "
                "equation is monic in the unused coefficient u_i on a "
                "central associated-graded packet polynomial in u_i"
            ),
        },
        "exact_data": exact_data,
        "exact_data_sha256": canonical_sha256(exact_data),
        "proved": [
            "the smooth formally rigid central different has six minimal generators",
            "the displayed tensors generate each complete graded gauge cokernel",
            "every compatible higher tail has a recursive multi-coefficient formal normal form",
            "C2 passes after every polynomial or formal graph specialization",
            "the intrinsic Kahler different has six minimal generators after every such specialization",
            "the Kahler different is non-Cartier to all formal orders for every squarefree cubic symbol",
        ],
        "not_proved": [
            "normality of the normal-form algebras",
            "algebraization of arbitrary infinite formal gauges",
            "existence or compatibility of a Keller open",
        ],
        "reproduce": (
            ".venv/bin/python scripts/"
            "verify_singular_squarefree_all_orders_different_persistence.py"
        ),
    }
    if args.refresh:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(artifact, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        pinned = json.loads(OUTPUT.read_text(encoding="utf-8"))
        assert artifact == pinned, "stale squarefree all-orders artifact"

    print(
        "PASS: smooth has coefficient count 0 and six-generator "
        "all-orders different by formal rigidity"
    )
    for name, sources in SOURCES.items():
        print(
            "PASS: "
            f"{name} has normal-form coefficient count {len(sources)} and "
            "six-generator all-orders different"
        )
    action = "wrote" if args.refresh else "replayed"
    print(f"PASS: {action} {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
