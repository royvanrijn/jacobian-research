#!/usr/bin/env python3
"""Compile the repository support-saturation case studies and search queues.

The cases have intentionally different theorem scopes:

* ``cubic`` runs the homogeneous cubic-symbol orbit atlas.  This is an exact
  leading-model certificate, not a singular-squarefree nonhomogeneous
  higher-lift theorem.
* ``cubic-frontier`` imports the proved formal-gauge cokernel atlas and
  compiles the remaining cubic work by annihilator type.  It is a routing
  certificate, not a new saturation computation.
* ``cubic-double-strata`` chooses exact complements to the quartic formal
  gauge images on the six singular squarefree symbols and computes both the
  cotangent-saturation and support-hull layers on the full complement
  families, together with the collision Nakayama quotient of the Kähler
  different.  Strict weighted-Rees packets promote both saturation layers,
  the annihilator quotient, and the six-generator non-Cartier conclusion to
  every geometric parameter fiber.  The case does not control higher formal
  orders or Keller compatibility.
* ``degree42`` uses the full unit-pivot-reduced core on the base fiber
  ``(e1,e2,t)=(1,2,3)``.  It certifies an untruncated characteristic-zero
  support class and computes the complete order-six/order-seven support
  saturations at the established good prime 32003.  It does not itself claim
  a characteristic-zero finite-jet lift or generic all-order synchronization;
  the separate sparse Macaulay certificate closes the fixed ``c6`` statement
  at order six only.
* ``plane`` runs the normalized cyclic ``d3`` multiplication-kernel layer of
  the Poisson-square coefficient scheme.  This is the currently defined
  boundary-residue proxy; the Case-1 conductor/residue matching module from
  OPEN_PROBLEMS_FOR_MAP_EXTENSIONS.md has not yet been constructed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.support_saturation import (  # noqa: E402
    CompilerOptions,
    ModulePresentation,
    NormalFiltration,
    PolynomialRing,
    SupportSaturationCompiler,
    SupportSaturationProblem,
    certificate_json,
)
from verify_cubic_symbol_double_saturation import (  # noqa: E402
    CUBIC_STRATA,
    differential_relations,
    homogeneous_monomials,
    quartic_kernel_basis_tensors,
    run_singular_rees_base_change_certificate,
    run_singular_subspace_certificate,
    singular_polynomial,
)
from verify_cubic_formal_gauge_cokernel_atlas import (  # noqa: E402
    symbol_tensor,
)
import verify_universal_cubic_cotangent_saturation as cubic_formal  # noqa: E402
from verify_degree42_order7_known_witness import C6  # noqa: E402
from verify_degree42_transported_27_normal_jets import (  # noqa: E402
    serialize,
    transformed_problem,
)


GENERATED = ROOT / "artifacts" / "generated-results"
FORMAL_GAUGE_ATLAS = GENERATED / "cubic_formal_gauge_cokernel_atlas.json"

SINGULAR_SQUAREFREE_SYMBOLS = (
    "nodal",
    "cuspidal",
    "line-transverse-conic",
    "line-tangent-conic",
    "triangle",
    "concurrent-lines",
)
DEGENERATE_SYMBOLS = ("double-line", "triple-line", "zero")

EXPECTED_SINGULAR_ANNIHILATORS = {
    "nodal": "(x)",
    "cuspidal": "(x^2)",
    "line-transverse-conic": "(yz)",
    "line-tangent-conic": "(y^3)",
    "triangle": "(xyz)",
    "concurrent-lines": "(x^3)",
}

EXPECTED_NONGAUGE_COMPLEMENT_INDICES = {
    "nodal": (0, 1),
    "cuspidal": (0, 1, 2, 5),
    "line-transverse-conic": (15, 19, 20, 23),
    "line-tangent-conic": (8, 12, 15, 16, 20, 21),
    "triangle": (0, 1, 15, 19, 20, 23),
    "concurrent-lines": (0, 1, 2, 3, 4, 5, 6, 7),
}


def cubic_case() -> dict[str, Any]:
    """Compile every homogeneous ternary-cubic orbit representative."""

    compiler = SupportSaturationCompiler(
        CompilerOptions(
            associated_primes="decompose",
            saturation_strategy="compute",
            timeout_seconds=300,
        )
    )
    certificates: dict[str, Any] = {}
    for name, cubic in CUBIC_STRATA.items():
        relations = differential_relations(cubic)
        generators = tuple(
            tuple(singular_polynomial(entry) for entry in relation)
            for relation in relations
        )
        presentation = ModulePresentation(
            ring=PolynomialRing(("x", "y", "z")),
            rank=12,
            generators=generators,
            label=f"homogeneous-cubic-cotangent-{name}",
        )
        certificates[name] = compiler.compile_problem(
            SupportSaturationProblem(
                presentation=presentation,
                support_ideal=("x", "y", "z"),
                parameter_base_variables=(),
                normal_variables=("x", "y", "z"),
            )
        )
    assert all(
        item["saturation"]["equal_to_presentation"]
        for item in certificates.values()
    )
    return {
        "schema": "support-saturation-case-family.v1",
        "case": "universal-cubic-homogeneous-symbol-atlas",
        "mathematical_scope": (
            "Exact cotangent saturation and associated primes for the ten "
            "homogeneous ternary-cubic orbit representatives; this does not "
            "prove saturation for arbitrary nonhomogeneous higher lifts."
        ),
        "certificates": certificates,
    }


def cubic_frontier_case() -> dict[str, Any]:
    """Route cubic normalization work using the formal-gauge theorem."""

    atlas = json.loads(FORMAL_GAUGE_ATLAS.read_text())
    rows = atlas["rows"]
    assert rows["smooth"]["formal_rigidity_above_degree_three"]
    assert rows["smooth"]["quartic_nongauge_dimension"] == 0
    assert {
        name: rows[name]["annihilator"]
        for name in SINGULAR_SQUAREFREE_SYMBOLS
    } == EXPECTED_SINGULAR_ANNIHILATORS
    assert all(
        not rows[name]["formal_rigidity_above_degree_three"]
        and rows[name]["support_dimension"] == 2
        and rows[name]["quartic_nongauge_dimension"] > 0
        for name in SINGULAR_SQUAREFREE_SYMBOLS
    )
    assert all(
        rows[name]["annihilator"] == "(0)"
        and rows[name]["support_dimension"] == 3
        for name in DEGENERATE_SYMBOLS
    )

    singular_queue = [
        {
            "annihilator_type": rows[name]["annihilator"],
            "symbol": name,
            "quartic_nongauge_dimension": rows[name][
                "quartic_nongauge_dimension"
            ],
            "support_dimension": rows[name]["support_dimension"],
            "support_multiplicity": rows[name]["support_multiplicity"],
            "next_adapter_input": (
                "deformation-dependent cotangent presentation restricted "
                "to a representative of this nongauge quotient"
            ),
        }
        for name in SINGULAR_SQUAREFREE_SYMBOLS
    ]
    degenerate_queue = [
        {
            "symbol": name,
            "annihilator_type": rows[name]["annihilator"],
            "quartic_nongauge_dimension": rows[name][
                "quartic_nongauge_dimension"
            ],
            "generic_rank_over_Q[x,y,z]": rows[name][
                "generic_rank_over_Q[x,y,z]"
            ],
            "gate_before_saturation": (
                "generically etale and Keller compatibility"
            ),
        }
        for name in DEGENERATE_SYMBOLS
    ]
    return {
        "schema": "support-saturation-search-frontier.v1",
        "case": "cubic-formal-gauge-annihilator-frontier",
        "source_artifact": str(FORMAL_GAUGE_ATLAS.relative_to(ROOT)),
        "source_exact_matrix_sha256": atlas["exact_matrix_sha256"],
        "mathematical_scope": (
            "Workflow compiled from the proved formal-gauge cokernel atlas. "
            "It closes further quartic saturation searches for the smooth "
            "symbol, queues the six singular squarefree cases by exact "
            "annihilator type on their finite-dimensional quartic nongauge "
            "quotients, and places the generically-etale/Keller gate before "
            "saturation for double-line, triple-line, and zero symbols. "
            "This artifact does not prove singular-symbol saturation or "
            "Keller compatibility."
        ),
        "stratification": {
            "smooth": {
                "symbol": "smooth",
                "quartic_nongauge_dimension": 0,
                "universal_24_parameter_cotangent_saturation": "proved",
                "next_tasks": [
                    "global algebraization of the formal gauge",
                    "boundary and Keller-open compatibility",
                ],
            },
            "singular_squarefree": {
                "search_key": "annihilator_type",
                "queue": singular_queue,
            },
            "double_triple_line_and_zero": {
                "search_key": "compatibility_before_saturation",
                "queue": degenerate_queue,
            },
        },
    }


def _quartic_tensor_column(
    tensor: dict[tuple[int, int, int], sp.Expr],
) -> sp.Matrix:
    monomials = homogeneous_monomials(4)
    return sp.Matrix(
        [
            sp.Poly(
                sp.expand(tensor[triple]),
                *cubic_formal.cubic_audit.BASE_VARIABLES,
            ).coeff_monomial(monomial)
            for triple in cubic_formal.TRIPLES
            for monomial in monomials
        ]
    )


def _tensor_record(
    tensor: dict[tuple[int, int, int], sp.Expr],
) -> list[str]:
    return [
        sp.sstr(sp.expand(tensor[triple]))
        for triple in cubic_formal.TRIPLES
    ]


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def _nongauge_complement(
    name: str,
) -> tuple[
    tuple[int, ...],
    tuple[dict[tuple[int, int, int], sp.Expr], ...],
    dict[str, int],
]:
    """Choose a deterministic complement of the quartic gauge image."""

    basis = quartic_kernel_basis_tensors()
    compatible_matrix = sp.Matrix.hstack(
        *(_quartic_tensor_column(tensor) for tensor in basis)
    )
    gauge = cubic_formal.gauge_matrix(symbol_tensor(CUBIC_STRATA[name]))
    variables = cubic_formal.cubic_audit.BASE_VARIABLES
    action_matrix = sp.Matrix.hstack(
        *(
            _quartic_tensor_column(
                {
                    triple: sp.expand(gauge[row, column] * variable)
                    for row, triple in enumerate(cubic_formal.TRIPLES)
                }
            )
            for column in range(gauge.cols)
            for variable in variables
        )
    )
    compatible_rank = compatible_matrix.rank()
    action_rank = action_matrix.rank()
    selected: list[int] = []
    span = action_matrix
    span_rank = action_rank
    for index in range(compatible_matrix.cols):
        enlarged = span.row_join(compatible_matrix[:, index])
        enlarged_rank = enlarged.rank()
        if enlarged_rank > span_rank:
            selected.append(index)
            span = enlarged
            span_rank = enlarged_rank
    selected_indices = tuple(selected)
    assert compatible_rank == 24
    assert span_rank == compatible_rank
    assert selected_indices == EXPECTED_NONGAUGE_COMPLEMENT_INDICES[name]
    return (
        selected_indices,
        tuple(basis[index] for index in selected_indices),
        {
            "compatible_quartic_dimension": compatible_rank,
            "linear_gauge_image_dimension": action_rank,
            "nongauge_quotient_dimension": compatible_rank - action_rank,
        },
    )


def cubic_double_stratification_case() -> dict[str, Any]:
    """Compute both cubic saturation layers on every nongauge quartic row."""

    # Expanded serialization is dramatically faster for the eight-parameter
    # concurrent-lines row and leaves the exact polynomial input unchanged.
    cubic_formal.cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    atlas = json.loads(FORMAL_GAUGE_ATLAS.read_text())
    basis = quartic_kernel_basis_tensors()
    basis_record = [_tensor_record(tensor) for tensor in basis]
    rows: dict[str, Any] = {}
    for name in SINGULAR_SQUAREFREE_SYMBOLS:
        indices, representatives, quotient = _nongauge_complement(name)
        assert quotient["nongauge_quotient_dimension"] == atlas["rows"][name][
            "quartic_nongauge_dimension"
        ]
        computation = run_singular_subspace_certificate(
            CUBIC_STRATA[name],
            representatives,
            timeout=1800,
        )
        base_change = run_singular_rees_base_change_certificate(
            CUBIC_STRATA[name],
            representatives,
            timeout=1800,
        )
        parameter_count = len(indices)
        expected = {
            "parameter_count": parameter_count,
            "cotangent_saturation_generators": 0,
            "support_module_dimension": parameter_count + 2,
            "support_ext3_vector_dimension": 0,
            "support_ext2_dimension": parameter_count,
            "support_ext2_multiplicity": 6,
            "support_ext2_parameter_axis_radical_difference": 0,
            "support_ext2_central_pruned_presentation_difference": 0,
            "support_ext2_pruned_presentation_rank": 3,
            "support_ext2_collision_square_action_generators": 0,
            "different_generator_module_dimension": parameter_count,
            "different_generator_module_multiplicity": 6,
            "different_generator_parameter_axis_radical_difference": 0,
            "different_generator_central_pruned_presentation_difference": 0,
            "different_generator_pruned_presentation_rank": 6,
        }
        assert computation == expected, (name, computation)
        expected_base_change = {
            "parameter_count": parameter_count,
            "cotangent_rees_torsion_generators": 0,
            "cotangent_initial_presentation_difference": 0,
            "annihilator_cokernel_rees_torsion_generators": 0,
            "annihilator_cokernel_initial_presentation_difference": 0,
        }
        assert base_change == expected_base_change, (name, base_change)
        representative_record = [
            _tensor_record(tensor) for tensor in representatives
        ]
        rows[name] = {
            "formal_gauge_cokernel_annihilator": atlas["rows"][name][
                "annihilator"
            ],
            **quotient,
            "compatible_basis_indices": list(indices),
            "representative_tensor_components": representative_record,
            "representative_sha256": _canonical_sha256(
                representative_record
            ),
            "exact_computation": computation,
            "fiberwise_base_change_certificate": base_change,
            "quartic_model_gate_results": {
                "C1_support_hull": (
                    "fails_on_every_geometric_parameter_fiber: strict Rees "
                    "base change identifies the intrinsic support module, "
                    "and Ext^2 is the parameter-independent nonzero "
                    "multiplicity-six central module"
                ),
                "C2_cotangent_torsion": (
                    "passes_on_every_geometric_parameter_fiber: the strict "
                    "Rees initial module is the saturated central cotangent "
                    "module"
                ),
                "Cartier_Kahler_different": (
                    "fails_on_every_geometric_parameter_fiber: the local "
                    "Kahler different has six minimal generators at the "
                    "collision"
                ),
            },
            "base_change_status": (
                "certified for every geometric specialization by the "
                "t-saturated cotangent and annihilator-cokernel Rees modules"
            ),
            "programme_status": (
                "C1 remains open for Keller-compatible cubic "
                "normalizations; this quartic model family supplies an "
                "obstruction that the geometric hypotheses must exclude"
            ),
        }
    return {
        "schema": "cubic-double-saturation-stratification.v3",
        "case": "singular-squarefree-cubic-quartic-double-saturation",
        "source_artifact": str(FORMAL_GAUGE_ATLAS.relative_to(ROOT)),
        "source_exact_matrix_sha256": atlas["exact_matrix_sha256"],
        "compatible_quartic_basis_sha256": _canonical_sha256(basis_record),
        "tensor_component_order": [
            list(triple) for triple in cubic_formal.TRIPLES
        ],
        "mathematical_scope": (
            "Exact characteristic-zero calculation on a deterministic "
            "linear complement of the formal gauge image in the complete "
            "quartic compatible-tensor space for each singular squarefree "
            "cubic symbol. It proves universal cotangent saturation on "
            "each displayed quartic nongauge family and a nonzero, "
            "parameter-independent multiplicity-six Ext^2 support-hull "
            "obstruction. The same exact presentations compute the "
            "collision Nakayama module J/nJ for the Kahler different "
            "J=Ann_B(Omega): it is the scalar extension of a six-dimensional "
            "central module, so J is not locally principal on any geometric "
            "parameter fiber. Strict Rees certificates for Omega and for the "
            "cokernel of B -> Omega^3 prove flatness over the parameter "
            "ring, commute the annihilator and support module with every "
            "geometric base change, and promote both gate results to all "
            "parameter fibers. It does not classify higher-order formal "
            "deformations or prove Keller-open compatibility."
        ),
        "rows": rows,
        "uniform_conclusion": {
            "C2_on_quartic_models": (
                "proved on every geometric fiber of all six families"
            ),
            "C1_on_quartic_models": (
                "fails on every geometric fiber of all six families"
            ),
            "fiberwise_base_change": "proved for all six quartic families",
            "Cartier_Kahler_different_on_quartic_models": (
                "fails on every geometric fiber of all six families"
            ),
            "global_cubic_programme": "open",
        },
    }


def degree42_case(characteristic: int = 32003) -> dict[str, Any]:
    """Compile the degree-42 support witness and finite-jet fiber."""

    normals, bases, residuals, _defect = transformed_problem()
    source_variables = normals + bases
    e1, e2, translation, w0, w1, w2 = bases
    setup = f"""
ring source={characteristic},({",".join(map(str,source_variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
poly p3=subst(I[5],{normals[2]},0);
ideal I3=subst(I,{normals[2]},p3);
poly p4=subst(I3[11],{normals[3]},0);
ideal I4=subst(I3,{normals[3]},p4);
poly p5=subst(I4[17],{normals[4]},0);
ideal I5=subst(I4,{normals[4]},p5);
ideal ReducedCoreSource=I5;

ring q={characteristic},(u,v,{w0},{w1},{w2}),(dp(2),dp(3));
map phi=source,u,v,0,0,0,1,2,3,{w0},{w1},{w2};
ideal ReducedCore=phi(ReducedCoreSource);
module InputPresentation=module(ReducedCore);
"""
    # At (e1,e2,t)=(1,2,3), A=-12 and B=1359, so AB=-16308.
    presentation = ModulePresentation(
        ring=PolynomialRing(
            ("u", "v", str(w0), str(w1), str(w2)),
            characteristic=characteristic,
            ordering="(dp(2),dp(3))",
        ),
        rank=1,
        generators=(),
        label="degree42-ritt-fiber-e1_1-e2_2-t_3",
        singular_setup=setup,
    )
    compiler = SupportSaturationCompiler(
        CompilerOptions(
            associated_primes="regularity",
            saturation_strategy="compute",
            timeout_seconds=3600 if characteristic == 0 else 1200,
            basis_algorithm="slimgb",
            torsion_exponent_bound=4,
        )
    )
    certificate = compiler.compile_distinguished_jets(
        presentation,
        (
            str(w0),
            f"{w1}*{w2}",
            f"-16308*{w2}",
        ),
        (C6,),
        NormalFiltration(
            ("u", "v"),
            (6, 7),
            strategy="full_saturation",
        ),
    )
    support_setup = f"""
ring support_source=0,({",".join(map(str,source_variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
poly p3=subst(I[5],{normals[2]},0);
ideal I3=subst(I,{normals[2]},p3);
poly p4=subst(I3[11],{normals[3]},0);
ideal I4=subst(I3,{normals[3]},p4);
poly p5=subst(I4[17],{normals[4]},0);
ideal I5=subst(I4,{normals[4]},p5);

ring support_homogeneous=0,(u,v,{w0},{w1},{w2},H),dp;
map support_phi=support_source,u,v,0,0,0,1,2,3,{w0},{w1},{w2};
ideal SupportCore=support_phi(I5);
ideal HomogeneousCore=homog(SupportCore,H);
ideal HomogeneousBasis=groebner(HomogeneousCore);
ideal DehomogeneousBasis=std(subst(HomogeneousBasis,H,1));

ring support_ring=0,(u,v,{w0},{w1},{w2}),dp;
map dehomogenize=support_homogeneous,u,v,{w0},{w1},{w2},1;
ideal SupportBasis=std(dehomogenize(DehomogeneousBasis));
module InputPresentation=module(SupportBasis);
attrib(InputPresentation,"isSB",1);
"""
    support_presentation = ModulePresentation(
        ring=PolynomialRing(("u", "v", str(w0), str(w1), str(w2))),
        rank=1,
        generators=(),
        label="degree42-exact-core-support-witness-e1_1-e2_2-t_3",
        singular_setup=support_setup,
    )
    support_witness = SupportSaturationCompiler(
        CompilerOptions(
            timeout_seconds=600,
            torsion_exponent_bound=4,
        )
    ).compile_distinguished_support_witness(
        support_presentation,
        (str(w0), str(w2)),
        (f"v*{w0}^2*(2*u-5*v)",),
    )
    return {
        "schema": "support-saturation-case-family.v1",
        "case": "degree42-ritt-synchronization-finite-jet-fiber",
        "characteristic": characteristic,
        "base_specialization": {
            str(e1): 1,
            str(e2): 2,
            str(translation): 3,
        },
        "untruncated_base": {
            "status": "nonzero_support_witness_certified",
            "full_saturation_module_computed": False,
            "support_witness": support_witness,
            "observed_attempts": [
                {
                    "characteristic": 0,
                    "presentation": "terminal-two-residual",
                    "outcome": "timeout",
                    "seconds": 3600,
                },
                {
                    "characteristic": 32003,
                    "presentation": "full-reduced-core",
                    "method": "direct-colon",
                    "outcome": "timeout",
                    "seconds": 600,
                },
                {
                    "characteristic": 0,
                    "presentation": "full-reduced-core",
                    "method": "full-order-six/order-seven-saturation",
                    "outcome": "timeout",
                    "seconds": 3600,
                },
            ],
        },
        "mathematical_scope": (
            "An exact characteristic-zero nonzero support class for the "
            "untruncated full reduced core, plus an exact distinguished-class "
            f"finite-jet computation in characteristic {characteristic} at "
            "orders six and seven. The full saturation module and associated "
            "prime list are not computed. The finite-jet computation does "
            "not lift itself to characteristic zero or imply generic "
            "completed synchronization."
        ),
        "certificate": certificate,
    }


PLANE_VARIABLES = (
    "a0",
    "a1",
    "a2",
    "a3",
    "b0",
    "b1",
    "b2",
    "b3",
    "b4",
    "c0",
    "c1",
    "c2",
    "d1",
    "d2",
    "d3",
)

PLANE_RELATIONS = (
    "2*a2*c2-3*a3*c1",
    "4*a1*c2-a2*c1-6*a3*c0",
    "6*a0*c2+a1*c1-4*a2*c0",
    "3*a0*c1-2*a1*c0",
    "-4*b4*c2+6*a3*d3",
    "-2*b3*c2-6*b4*c1+7*a2*d3+3*a3*d2",
    "-4*b3*c1-8*b4*c0+8*a1*d3+4*a2*d2",
    (
        "2*b1*c2-2*b2*c1-6*b3*c0+9*a0*d3+5*a1*d2"
        "+a2*d1-3*a3"
    ),
    "4*b0*c2-4*b2*c0+6*a0*d2+2*a1*d1-2*a2",
    "2*b0*c1-2*b1*c0+3*a0*d1-a1",
    "2*b4*d3",
    "3*b3*d3",
    "4*b2*d3+b3*d2-2*b4*d1",
    "5*b1*d3+2*b2*d2-b3*d1-4*b4",
    "6*b0*d3+3*b1*d2-3*b3",
    "4*b0*d2+b1*d1-2*b2",
    "2*b0*d1-b1-1",
)

PLANE_D3_SUPPORT = (
    "d3",
    "d2",
    "c2",
    "c1",
    "b4",
    "b3",
    "a3",
    "a2",
    "a1",
    "b1*d1+4*b0*d2-2*b2",
    "2*b0*d1-b1-1",
    "a0*d1",
    "b1*c0",
    "a0*c0",
)


def plane_case() -> dict[str, Any]:
    """Compile the normalized cyclic d3 boundary-torsion layer."""

    setup = f"""
ring plane_ring=0,({",".join(PLANE_VARIABLES)}),dp;
ideal CoefficientIdeal={",".join(PLANE_RELATIONS)};
ideal D3Annihilator=std(quotient(CoefficientIdeal,ideal(d3)));
module InputPresentation=module(D3Annihilator);
"""
    presentation = ModulePresentation(
        ring=PolynomialRing(PLANE_VARIABLES),
        rank=1,
        generators=(),
        label="plane-jc-normalized-d3-multiplication-kernel",
        singular_setup=setup,
    )
    compiler = SupportSaturationCompiler(
        CompilerOptions(
            associated_primes="decompose",
            saturation_strategy="compute",
            timeout_seconds=3600,
            regular_search_bound=0,
        )
    )
    certificate = compiler.compile_problem(
        SupportSaturationProblem(
            presentation=presentation,
            support_ideal=PLANE_D3_SUPPORT,
            parameter_base_variables=PLANE_VARIABLES,
            normal_variables=(),
            distinguished_class=("1",),
        )
    )
    assert not certificate["local_cohomology"]["zero"]
    assert certificate["distinguished_class"][
        "belongs_to_local_cohomology"
    ]
    return {
        "schema": "support-saturation-case-family.v1",
        "case": "plane-jc-poisson-square-d3-boundary-layer",
        "mathematical_scope": (
            "Exact normalized cyclic d3 multiplication-kernel layer of the "
            "Poisson-square coefficient algebra. This is not the undefined "
            "Case-1 conductor/residue matching module."
        ),
        "certificate": certificate,
    }


CASES = {
    "cubic": (
        cubic_case,
        "support_saturation_universal_cubic_symbols.json",
    ),
    "cubic-frontier": (
        cubic_frontier_case,
        "support_saturation_cubic_annihilator_frontier.json",
    ),
    "cubic-double-strata": (
        cubic_double_stratification_case,
        "cubic_double_saturation_stratification.json",
    ),
    "degree42": (
        degree42_case,
        "support_saturation_degree42_ritt_fiber_mod32003.json",
    ),
    "plane": (
        plane_case,
        "support_saturation_plane_jc_boundary_layer.json",
    ),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=("all", *CASES),
        default="all",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=GENERATED,
    )
    parser.add_argument(
        "--degree42-characteristic",
        type=int,
        default=32003,
        help="use 0 for the long characteristic-zero degree-42 run",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected = CASES if args.case == "all" else {args.case: CASES[args.case]}
    for name, (builder, filename) in selected.items():
        result = (
            builder(args.degree42_characteristic)
            if name == "degree42"
            else builder()
        )
        if name == "degree42" and args.degree42_characteristic != 32003:
            suffix = (
                "char0"
                if args.degree42_characteristic == 0
                else f"mod{args.degree42_characteristic}"
            )
            filename = f"support_saturation_degree42_ritt_fiber_{suffix}.json"
        output = args.output_dir / filename
        output.write_text(certificate_json(result))
        try:
            displayed_output = output.relative_to(ROOT)
        except ValueError:
            displayed_output = output
        compact = {
            "case": result["case"],
            "output": str(displayed_output),
        }
        print(json.dumps(compact, sort_keys=True))
        print(f"PASS: compiled support-saturation case {name}")


if __name__ == "__main__":
    main()
