#!/usr/bin/env sage
"""Analyze a split-prime q=80 slope-1/12 sextic relation space.

The input is a finite formal relation space, so the resulting ideal geometry
is a modular candidate only.  This deliberately keeps to fast invariants:
affine dimension, projective Hilbert polynomial, and an optional two-coordinate
elimination.
"""

import argparse
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--plane-pair", nargs=2, metavar=("X", "Y"))
parser.add_argument("--plane-delta", action="store_true")
parser.add_argument("--plane-normalize", action="store_true")
arguments = parser.parse_args()
payload = json.loads(Path(arguments.input).read_text())
if (
    payload.get("schema") != "q80-cm24-split-prime-formal-branch-v1"
    or payload.get("kind") != "canonical_centered_relation_space"
):
    raise ValueError("unexpected split-prime relation artifact")

field = GF(payload["prime"])
names = tuple(payload["centered_variables"])
ring = PolynomialRing(field, names=names, order="degrevlex")
variables = ring.gens()
monomials = tuple(
    ring.prod(variable**exponent for variable, exponent in zip(variables, exponents))
    for exponents in payload["monomial_exponents"]
)
generators = tuple(
    sum(
        (field(coefficient)*monomial for coefficient, monomial in zip(row, monomials)),
        ring.zero(),
    )
    for row in payload["rref_basis"]
)
ideal = ring.ideal(generators)

started = time.monotonic()
dimension = ideal.dimension()
print(
    f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage=affine|"
    f"generators={len(generators)}|dimension={dimension}|"
    f"seconds={time.monotonic()-started:.3f}",
    flush=True,
)

started = time.monotonic()
# Homogenizing the affine ideal (rather than each displayed generator)
# computes the actual projective closure via a Groebner basis.
projective_ideal = ideal.homogenize()
hilbert = projective_ideal.hilbert_polynomial()
print(
    f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage=projective|"
    f"dimension={projective_ideal.dimension()}|hilbert={hilbert}|"
    f"degree={hilbert[1] if hilbert.degree()==1 else 'NA'}|"
    f"arithmetic_genus={1-hilbert[0] if hilbert.degree()==1 else 'NA'}|"
    f"seconds={time.monotonic()-started:.3f}",
    flush=True,
)

if arguments.plane_pair:
    keep = tuple(arguments.plane_pair)
    if len(set(keep)) != 2 or not set(keep) <= set(names):
        raise ValueError("plane coordinates must be two distinct centered variables")
    eliminated = tuple(variable for variable in variables if str(variable) not in keep)
    started = time.monotonic()
    elimination = ideal.elimination_ideal(eliminated)
    plane_ring = PolynomialRing(field, names=keep, order="degrevlex")
    plane_generators = tuple(plane_ring(str(value)) for value in elimination.gens())
    plane_ideal = plane_ring.ideal(plane_generators)
    factors = ()
    if len(plane_generators) == 1:
        factors = tuple(
            (factor.total_degree(), multiplicity)
            for factor, multiplicity in plane_generators[0].factor()
        )
    print(
        f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage=plane|"
        f"coordinates={','.join(keep)}|generators={len(plane_generators)}|"
        f"degrees={tuple(value.total_degree() for value in plane_generators)}|"
        f"bidegrees={tuple((value.degree(plane_ring.gen(0)), value.degree(plane_ring.gen(1))) for value in plane_generators)}|"
        f"factors={factors}|seconds={time.monotonic()-started:.3f}",
        flush=True,
    )
    if arguments.plane_delta:
        if len(plane_generators) != 1 or len(factors) != 1 or factors[0][1] != 1:
            raise ValueError("delta analysis needs one irreducible plane equation")
        from sage.libs.singular.function_factory import ff

        delta_loc = ff.normal__lib.deltaLoc

        def chart_delta(polynomial, extra_equations, label):
            chart_ring = polynomial.parent()
            singular = chart_ring.ideal(
                [polynomial]
                + [polynomial.derivative(variable) for variable in chart_ring.gens()]
                + list(extra_equations)
            )
            if singular == chart_ring.ideal(1):
                print(
                    f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage={label}|"
                    "components=0|delta=0",
                    flush=True,
                )
                return {"delta": 0, "rational_points": 0, "rational_branches": 0}
            components = singular.minimal_associated_primes()
            values = tuple(
                tuple(delta_loc(polynomial, component)) for component in components
            )
            degrees = tuple(component.vector_space_dimension() for component in components)
            delta = sum(int(value[0]) for value in values)
            rational_points = sum(degree == 1 for degree in degrees)
            rational_branches = sum(
                int(value[2])
                for value, degree in zip(values, degrees)
                if degree == 1
            )
            print(
                f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage={label}|"
                f"components={len(components)}|component_degrees={degrees}|"
                f"invariants={values}|delta={delta}|"
                f"rational_points={rational_points}|rational_branches={rational_branches}",
                flush=True,
            )
            return {
                "delta": delta,
                "rational_points": rational_points,
                "rational_branches": rational_branches,
            }

        delta_started = time.monotonic()
        affine_data = chart_delta(plane_generators[0], (), "delta_affine")
        projective_ring = PolynomialRing(field, names=("X", "Y", "Z"), order="degrevlex")
        X, Y, Z = projective_ring.gens()
        degree = plane_generators[0].total_degree()
        homogeneous_polynomial = projective_ring.zero()
        for exponents, coefficient in plane_generators[0].dict().items():
            homogeneous_polynomial += (
                field(coefficient)*X**exponents[0]*Y**exponents[1]
                *Z**(degree-sum(exponents))
            )
        y_chart_ring = PolynomialRing(field, names=("y", "z"), order="degrevlex")
        y, z = y_chart_ring.gens()
        y_chart = y_chart_ring(homogeneous_polynomial(1, y, z))
        infinity_data = chart_delta(y_chart, (z,), "delta_infinity")
        x_chart_ring = PolynomialRing(field, names=("x", "z"), order="degrevlex")
        x, z2 = x_chart_ring.gens()
        x_chart = x_chart_ring(homogeneous_polynomial(x, 1, z2))
        x_origin = x_chart_ring.ideal(x, z2)
        if all(value in x_origin for value in (
            x_chart, x_chart.derivative(x), x_chart.derivative(z2)
        )):
            extra_data = chart_delta(x_chart, (x, z2), "delta_infinity_extra")
        else:
            extra_data = {"delta": 0, "rational_points": 0, "rational_branches": 0}
        total_delta = sum(data["delta"] for data in (affine_data, infinity_data, extra_data))
        geometric_genus = (degree-1)*(degree-2)//2-total_delta
        print(
            f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage=delta_total|"
            f"plane_degree={degree}|delta={total_delta}|"
            f"geometric_genus={geometric_genus}|"
            f"seconds={time.monotonic()-delta_started:.3f}",
            flush=True,
        )
        projective_points = sum(
            not homogeneous_polynomial(1, value, 0)
            for value in field
        ) + int(not homogeneous_polynomial(0, 1, 0))
        projective_points += sum(
            not plane_generators[0](left, right)
            for left in field for right in field
        )
        rational_singular_points = sum(
            data["rational_points"] for data in (affine_data, infinity_data, extra_data)
        )
        rational_branches = sum(
            data["rational_branches"] for data in (affine_data, infinity_data, extra_data)
        )
        normalization_points = (
            projective_points-rational_singular_points+rational_branches
        )
        trace = int(payload["prime"])+1-normalization_points
        print(
            f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage=point_count|"
            f"singular_model_points={projective_points}|"
            f"rational_singular_points={rational_singular_points}|"
            f"rational_branches={rational_branches}|"
            f"normalization_points={normalization_points}|trace={trace}",
            flush=True,
        )
        if int(payload["prime"]) != 79:
            # The two elliptic quotients of
            # u^2=16*t^6-19*t^4+88*t^2-48.  The first uses x=t^2;
            # the second uses x=1/t^2.  Any elliptic quotient of the source
            # Jacobian at a good prime must have one of these Frobenius traces.
            source_quotients = (
                EllipticCurve(QQ, [0, -19, 0, 1408, -12288]),
                EllipticCurve(QQ, [0, 88, 0, 912, 36864]),
            )
            source_traces = tuple(
                curve.change_ring(field).trace_of_frobenius()
                for curve in source_quotients
            )
            print(
                f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage=source_compare|"
                f"candidate_trace={trace}|source_quotient_traces={source_traces}|"
                f"direct_match={int(trace in source_traces)}|"
                f"match_up_to_quadratic_twist={int(trace in source_traces or -trace in source_traces)}",
                flush=True,
            )
    if arguments.plane_normalize:
        if len(plane_generators) != 1 or len(factors) != 1 or factors[0][1] != 1:
            raise ValueError("normalization needs one irreducible plane equation")
        from sage.libs.singular.function_factory import ff

        normal_started = time.monotonic()
        normalization = ff.normal__lib.normalP(plane_ideal)
        module = tuple(normalization[0][0])
        print(
            f"Q80CM24SPLITIDEAL|prime={payload['prime']}|stage=normalization|"
            f"components={len(normalization[0])}|module_generators={len(module)}|"
            f"module_degrees={tuple(value.total_degree() for value in module)}|"
            f"delta_data={normalization[1]}|seconds={time.monotonic()-normal_started:.3f}",
            flush=True,
        )
