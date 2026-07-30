#!/usr/bin/env python3
"""Compile the three repository support-saturation case studies.

The cases have intentionally different theorem scopes:

* ``cubic`` runs the homogeneous cubic-symbol orbit atlas.  This is an exact
  leading-model certificate, not the still-open arbitrary higher-order
  universal cubic theorem.
* ``degree42`` uses the full unit-pivot-reduced core on the base fiber
  ``(e1,e2,t)=(1,2,3)``.  It certifies an untruncated characteristic-zero
  support class and computes the complete order-six/order-seven support
  saturations at the established good prime 32003.  It does not claim a
  characteristic-zero finite-jet lift or generic all-order synchronization.
* ``plane`` runs the normalized cyclic ``d3`` multiplication-kernel layer of
  the Poisson-square coefficient scheme.  This is the currently defined
  boundary-residue proxy; the Case-1 conductor/residue matching module from
  OPEN_PROBLEMS_FOR_MAP_EXTENSIONS.md has not yet been constructed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

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
    certificate_json,
)
from verify_cubic_symbol_double_saturation import (  # noqa: E402
    CUBIC_STRATA,
    differential_relations,
    singular_polynomial,
)
from verify_degree42_order7_known_witness import C6  # noqa: E402
from verify_degree42_transported_27_normal_jets import (  # noqa: E402
    serialize,
    transformed_problem,
)


GENERATED = ROOT / "artifacts" / "generated-results"


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
        certificates[name] = compiler.compile(
            presentation,
            ("x", "y", "z"),
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
    certificate = compiler.compile(
        presentation,
        PLANE_D3_SUPPORT,
        distinguished_class=("1",),
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
        compact = {
            "case": result["case"],
            "output": str(output.relative_to(ROOT)),
        }
        print(json.dumps(compact, sort_keys=True))
        print(f"PASS: compiled support-saturation case {name}")


if __name__ == "__main__":
    main()
