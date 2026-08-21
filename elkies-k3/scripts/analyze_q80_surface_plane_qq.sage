#!/usr/bin/env sage
"""Certify and parametrize the exact characteristic-zero q80 D,Q plane.

The original projective model has a unique singular point at infinity, so the
adjoint algorithm must first be applied in a coordinate frame with smooth
intersection at infinity.  We use

    z_old = z_new - y,

which moves that point into the affine chart.  Singular's bounded local
delta computation in that frame certifies genus zero.  Singular's fast
choice-four adjoints are used only as an auxiliary linear system: their CM24
osculating ratio has degree three and produces a degree-eight bridge curve.
We parametrize that bridge, recover Q from a linear subresultant, and certify
the resulting birational map by exact substitution in both directions.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument(
    "--plane",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-DQ-plane.json",
)
parser.add_argument(
    "--series",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-order28-series.json",
)
parser.add_argument(
    "--adjoint-cache",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-DQ-adjoint-zminusy.json",
)
parser.add_argument(
    "--output",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-DQ-parameter.json",
)
parser.add_argument(
    "--bridge-cache",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-auxiliary-bridge-parameter.json",
)
parser.add_argument(
    "--composition-cache",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-reduced-composition.json",
)
parser.add_argument(
    "--reduced-q",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-reduced-Q.json",
)
parser.add_argument(
    "--timeout",
    type=int,
    default=35,
    help="hard timeout in seconds for Singular's local ideal-quotient adjoint call",
)
arguments = parser.parse_args()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def primitive_polynomial(polynomial):
    if not polynomial:
        return polynomial
    denominator = lcm(c.denominator() for c in polynomial.coefficients())
    integral = polynomial.parent()(denominator*polynomial)
    content = gcd(tuple(ZZ(c) for c in integral.coefficients()))
    integral = integral//content
    if integral.leading_coefficient() < 0:
        integral = -integral
    return integral


plane_path = Path(arguments.plane)
series_path = Path(arguments.series)
plane_payload = json.loads(plane_path.read_text())
series_payload = json.loads(series_path.read_text())
plane_hash = sha256(plane_path)
if plane_payload.get("slope") != "8/87":
    raise ValueError("plane artifact is not the slope 8/87 branch")
if series_payload.get("order") != len(series_payload["series"]["D"]):
    raise ValueError("series artifact has inconsistent order metadata")

affine_ring = PolynomialRing(QQ, names=("D", "Q"))
D, Q = affine_ring.gens()
plane_polynomial = affine_ring(plane_payload["polynomial"])
plane_degree = plane_polynomial.total_degree()
if plane_degree != 10:
    raise ValueError(f"expected a degree-ten plane, got degree {plane_degree}")

projective_ring = PolynomialRing(QQ, names=("x", "y", "z"))
x, y, z = projective_ring.gens()
projective_polynomial = projective_ring(
    sum(
        coefficient*x**exponents[0]*y**exponents[1]
        * z**(plane_degree-sum(exponents))
        for exponents, coefficient in plane_polynomial.dict().items()
    )
)

# New coordinates map to old coordinates by (x,y,z) -> (x,y,z-y).
adjoint_curve = projective_ring(projective_polynomial(z=z-y))
infinity_singularity_ideal = projective_ring.ideal(
    [adjoint_curve]
    + [adjoint_curve.derivative(variable) for variable in projective_ring.gens()]
    + [z]
)
if infinity_singularity_ideal.dimension() != 0:
    raise ArithmeticError("chosen coordinate change leaves a singularity at infinity")
if not adjoint_curve(x=0, y=1, z=0):
    raise ArithmeticError("chosen coordinate change leaves (0:1:0) on the plane")
integral_adjoint_curve = primitive_polynomial(adjoint_curve)


def parse_singular_polynomial(value):
    for variable in ("x", "y", "z"):
        value = re.sub(rf"{variable}([0-9]+)", rf"{variable}^\1", value)
    value = re.sub(r"([0-9])([xyz])", r"\1*\2", value)
    value = re.sub(r"([xyz])(?=[xyz])", r"\1*", value)
    return projective_ring(value)


adjoint_cache_path = Path(arguments.adjoint_cache)
adjoint_generators = None
if adjoint_cache_path.is_file():
    cache = json.loads(adjoint_cache_path.read_text())
    if (
        cache.get("schema") == "q80-cm24-qq-DQ-adjoint-v2"
        and cache.get("source_plane_sha256") == plane_hash
        and cache.get("coordinate_change") == "z_old=z_new-y"
        and cache.get("adjoint_choice") == 4
    ):
        adjoint_generators = tuple(
            projective_ring(value) for value in cache["generators"]
        )

sage_local = Path(os.environ.get("SAGE_LOCAL", ""))
singular_path = sage_local/"bin/Singular"
if adjoint_generators is None:
    if not singular_path.is_file():
        raise FileNotFoundError("cannot locate Sage's bundled Singular")
    singular_program = f'''
LIB "paraplanecurves.lib";
ring R=0,(x,y,z),dp;
poly f={integral_adjoint_curve};
ideal AI=adjointIdeal(f,4);
int adjoint_index;
for(adjoint_index=1;adjoint_index<=size(AI);adjoint_index++)
{{
  "Q80ADJOINT|"+string(adjoint_index)+"|"+string(AI[adjoint_index]);
}}
'''
    try:
        completed = subprocess.run(
            [str(singular_path), "-q"],
            input=singular_program,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=arguments.timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        partial = error.stdout or ""
        if isinstance(partial, bytes):
            partial = partial.decode(errors="replace")
        print(partial, end="", flush=True)
        raise TimeoutError(
            f"Singular adjoint computation exceeded the {arguments.timeout}s hard limit"
        ) from error
    if completed.returncode:
        print(completed.stdout, end="", flush=True)
        raise RuntimeError(
            f"Singular adjoint computation failed with exit code {completed.returncode}"
        )
    adjoint_strings = {}
    for line in completed.stdout.splitlines():
        if line.startswith("Q80ADJOINT|"):
            _, index, polynomial = line.split("|", 2)
            adjoint_strings[int(index)] = polynomial
    if not adjoint_strings:
        raise RuntimeError("Singular did not return the adjoint ideal")
    adjoint_generators = tuple(
        parse_singular_polynomial(adjoint_strings[index])
        for index in sorted(adjoint_strings)
    )
    adjoint_cache_path.parent.mkdir(parents=True, exist_ok=True)
    adjoint_cache_path.write_text(
        json.dumps(
            {
                "schema": "q80-cm24-qq-DQ-adjoint-v2",
                "source_plane": str(plane_path),
                "source_plane_sha256": plane_hash,
                "coordinate_change": "z_old=z_new-y",
                "adjoint_choice": int(4),
                "adjoint_method": "local_ideal_quotient",
                "generators": [str(value) for value in adjoint_generators],
            },
            indent=2,
            sort_keys=True,
        )+"\n"
    )
    print(
        "Q80SURFACEQQ|stage=adjoint|cache=written|"
        f"generators={len(adjoint_generators)}|path={adjoint_cache_path}",
        flush=True,
    )
else:
    print(
        "Q80SURFACEQQ|stage=adjoint|cache=hit|"
        f"generators={len(adjoint_generators)}|path={adjoint_cache_path}",
        flush=True,
    )

# The local-genus routine assumes the line at infinity is nonsingular.  Our
# transformed model satisfies that hypothesis, so this is an independent
# exact check on the degree-seven adjoint dimension.
library_path = sage_local/"share/singular/LIB/paraplanecurves.lib"
with tempfile.TemporaryDirectory(prefix="q80-plane-genus-") as temporary:
    wrapper_path = Path(temporary)/"q80_plane_genus.lib"
    wrapper_path.write_text(
        library_path.read_text()
        + '''
proc q80ExactGeomGenus(poly f)
{
  list geometry=geomGenusLA(f);
  return(geometry[1]);
}
'''
    )
    genus_program = f'''
LIB "{wrapper_path}";
ring R=0,(x,y,z),dp;
poly f={integral_adjoint_curve};
"Q80EXACTGENUS|"+string(q80ExactGeomGenus(f));
'''
    genus_completed = subprocess.run(
        [str(singular_path), "-q"],
        input=genus_program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=arguments.timeout,
        check=False,
    )
genus_lines = tuple(
    line for line in genus_completed.stdout.splitlines()
    if line.startswith("Q80EXACTGENUS|")
)
if genus_completed.returncode or len(genus_lines) != 1:
    print(genus_completed.stdout, end="", flush=True)
    raise RuntimeError("exact local-genus computation failed")
exact_geometric_genus = ZZ(genus_lines[0].split("|", 1)[1])
print(
    f"Q80SURFACEQQ|stage=exact_local_delta|geometric_genus={exact_geometric_genus}",
    flush=True,
)


def homogeneous_piece(target_degree):
    monomials = tuple(projective_ring.monomials_of_degree(target_degree))
    rows = []
    for generator in adjoint_generators:
        generator_degree = generator.total_degree()
        if generator_degree > target_degree:
            continue
        for multiplier in projective_ring.monomials_of_degree(
            target_degree-generator_degree
        ):
            polynomial = generator*multiplier
            rows.append(
                vector(
                    QQ,
                    [
                        polynomial.monomial_coefficient(monomial)
                        for monomial in monomials
                    ],
                )
            )
    if not rows:
        return ()
    space = span(QQ, rows)
    return tuple(
        primitive_polynomial(
            sum(coefficient*monomial for coefficient, monomial in zip(row, monomials))
        )
        for row in space.basis()
    )


choice4_degree_seven_forms = homogeneous_piece(plane_degree-3)
auxiliary_forms = homogeneous_piece(plane_degree-2)
if exact_geometric_genus != 0:
    raise ArithmeticError(
        f"transformed exact local-delta genus is {exact_geometric_genus}, expected zero"
    )
# Choice 4 is two conditions larger than the true conductor on this example:
# it has dimensions 2 and 11 in degrees 7 and 8.  We do not use those counts
# as a genus certificate.  The 11 degree-eight restrictions instead form an
# auxiliary O(10) system; its osculating ratio is accepted only after exact
# two-sided birational substitution below.
if len(choice4_degree_seven_forms) != 2:
    raise ArithmeticError(
        f"unexpected choice-4 degree-seven dimension {len(choice4_degree_seven_forms)}"
    )
if len(auxiliary_forms) != plane_degree+1:
    raise ArithmeticError(
        f"choice-4 degree-eight dimension is {len(auxiliary_forms)}, expected 11"
    )
auxiliary_degree = len(auxiliary_forms)-1

precision = int(series_payload["order"])
series_ring = PowerSeriesRing(QQ, "h", default_prec=precision)


def artifact_series(name, offset=QQ.zero()):
    values = [QQ(value) for value in series_payload["series"][name]]
    values[0] += offset
    return series_ring(values, precision)


D_series = artifact_series("D", QQ(1)/2)
Q_series = artifact_series("Q", QQ(9)/4)


def evaluate_adjoint_series(polynomial):
    # New homogeneous coordinates of an old affine point (D,Q,1) are
    # (D,Q,1+Q), since z_old=z_new-y.
    return sum(
        coefficient*D_series**exponents[0]*Q_series**exponents[1]
        * (1+Q_series)**exponents[2]
        for exponents, coefficient in polynomial.dict().items()
    )


adjoint_series = tuple(evaluate_adjoint_series(form) for form in auxiliary_forms)
jet_matrix = Matrix(
    QQ,
    [
        [coordinate[order] for coordinate in adjoint_series]
        for order in range(auxiliary_degree+1)
    ],
)
if jet_matrix.rank() != auxiliary_degree+1:
    raise ArithmeticError("auxiliary restrictions do not have a full CM24 jet")
high_kernel = jet_matrix.matrix_from_rows(range(auxiliary_degree)).right_kernel_matrix()
low_kernel = jet_matrix.matrix_from_rows(range(auxiliary_degree-1)).right_kernel_matrix()
if high_kernel.nrows() != 1 or low_kernel.nrows() != 2:
    raise ArithmeticError("unexpected osculating-flag dimensions")
high_vector = vector(QQ, high_kernel.row(0))
low_vector = next(
    vector(QQ, row)
    for row in low_kernel.rows()
    if Matrix(QQ, [high_vector, row]).rank() == 2
)


def normalized_vector(row):
    pivot = next(value for value in row if value)
    return row/pivot


high_vector = normalized_vector(high_vector)
low_vector = normalized_vector(low_vector)
high_form_new = primitive_polynomial(
    sum(coefficient*form for coefficient, form in zip(high_vector, auxiliary_forms))
)
low_form_new = primitive_polynomial(
    sum(coefficient*form for coefficient, form in zip(low_vector, auxiliary_forms))
)
high_series = evaluate_adjoint_series(high_form_new)
low_series = evaluate_adjoint_series(low_form_new)
if (high_series.valuation(), low_series.valuation()) != (
    auxiliary_degree, auxiliary_degree-1
):
    raise ArithmeticError("osculating auxiliary forms have the wrong CM24 orders")

# Express the parameter forms back in the original homogeneous coordinates.
high_form_old = primitive_polynomial(high_form_new(z=z+y))
low_form_old = primitive_polynomial(low_form_new(z=z+y))
parameter_precision = precision-low_series.valuation()
parameter_series = (high_series/low_series).add_bigoh(parameter_precision)
if parameter_series.valuation() != 1:
    raise ArithmeticError("osculating ratio is not a local parameter")
inverse_series = parameter_series.reverse().add_bigoh(parameter_precision)
D_in_parameter = D_series(inverse_series).add_bigoh(parameter_precision)
Q_in_parameter = Q_series(inverse_series).add_bigoh(parameter_precision)


def pade_candidate(series, numerator_degree, denominator_degree):
    coefficient_count = series.prec()
    coefficients = tuple(series[index] for index in range(coefficient_count))
    if coefficient_count < numerator_degree+denominator_degree+1:
        return None
    if denominator_degree == 0:
        if any(coefficients[numerator_degree+1:]):
            return None
        denominator = (QQ.one(),)
    else:
        rows = tuple(range(numerator_degree+1, coefficient_count))
        matrix_values = Matrix(
            QQ,
            [
                [
                    coefficients[index-offset] if index >= offset else QQ.zero()
                    for offset in range(1, denominator_degree+1)
                ]
                for index in rows
            ],
        )
        target = vector(QQ, [-coefficients[index] for index in rows])
        if matrix_values.rank() != denominator_degree:
            return None
        if matrix_values.augment(target).rank() != denominator_degree:
            return None
        denominator = (QQ.one(),)+tuple(matrix_values.solve_right(target))
    numerator = tuple(
        sum(
            denominator[offset]*coefficients[index-offset]
            for offset in range(min(denominator_degree, index)+1)
        )
        for index in range(numerator_degree+1)
    )
    return numerator, denominator


def minimal_pade_candidates(series, maximum_degree):
    for total_degree in range(2*maximum_degree+1):
        candidates = []
        for denominator_degree in range(maximum_degree+1):
            numerator_degree = total_degree-denominator_degree
            if not 0 <= numerator_degree <= maximum_degree:
                continue
            candidate = pade_candidate(
                series, numerator_degree, denominator_degree
            )
            if candidate is not None:
                candidates.append(
                    (numerator_degree, denominator_degree, candidate)
                )
        if candidates:
            return tuple(candidates)
    return ()


parameter_polynomial_ring = PolynomialRing(QQ, "t")
t = parameter_polynomial_ring.gen()
parameter_field = parameter_polynomial_ring.fraction_field()


def rational_function(result):
    numerator_degree, denominator_degree, (numerator, denominator) = result
    numerator_polynomial = sum(
        coefficient*t**index for index, coefficient in enumerate(numerator)
    )
    denominator_polynomial = sum(
        coefficient*t**index for index, coefficient in enumerate(denominator)
    )
    return parameter_field(numerator_polynomial/denominator_polynomial)


affine_high_form = affine_ring(high_form_old(x=D, y=Q, z=1))
affine_low_form = affine_ring(low_form_old(x=D, y=Q, z=1))


def live_coordinate_eliminant(keep_name, keep_series):
    """Eliminate the other plane coordinate and select the live factor."""
    coefficient_ring = PolynomialRing(QQ, names=("K", "T"))
    K, T = coefficient_ring.gens()
    elimination_ring = PolynomialRing(coefficient_ring, "E")
    E = elimination_ring.gen()

    def embed(polynomial):
        result = elimination_ring.zero()
        for (D_exponent, Q_exponent), coefficient in polynomial.dict().items():
            if keep_name == "D":
                result += coefficient*K**D_exponent*E**Q_exponent
            else:
                result += coefficient*E**D_exponent*K**Q_exponent
        return result

    embedded_plane = embed(plane_polynomial)
    embedded_parameter = embed(affine_high_form)-T*embed(affine_low_form)
    subresultants = tuple(embedded_plane.subresultants(embedded_parameter))
    linear_subresultants = tuple(
        value for value in subresultants if value.degree() == 1
    )
    if not linear_subresultants:
        raise ArithmeticError(f"{keep_name} elimination has no linear subresultant")
    resultant = coefficient_ring(
        embedded_plane.resultant(embedded_parameter)
    )
    factors = tuple(
        (primitive_polynomial(factor), int(exponent))
        for factor, exponent in resultant.factor()
    )

    def evaluate_factor(factor):
        return sum(
            coefficient*keep_series**exponents[0]
            * parameter_series**exponents[1]
            for exponents, coefficient in factor.dict().items()
        )

    live = tuple(
        factor
        for factor, exponent in factors
        if factor.degree(K) > 0
        and factor.degree(T) > 0
        and not evaluate_factor(factor)
    )
    if len(live) != 1:
        raise ArithmeticError(
            f"{keep_name} resultant has {len(live)} live factors; "
            f"degrees={tuple((factor.degree(K), factor.degree(T), exponent) for factor, exponent in factors)}"
        )
    return live[0], factors, linear_subresultants[-1]


D_eliminant, D_resultant_factors, Q_linear_subresultant = live_coordinate_eliminant(
    "D", D_series
)
bridge_ring = D_eliminant.parent()
bridge_K, bridge_T = bridge_ring.gens()
bridge_degree = D_eliminant.total_degree()
if (D_eliminant.degree(bridge_K), D_eliminant.degree(bridge_T), bridge_degree) != (
    3, 5, 8
):
    raise ArithmeticError(
        "unexpected auxiliary bridge degrees "
        f"{D_eliminant.degree(bridge_K)},{D_eliminant.degree(bridge_T)},{bridge_degree}"
    )

# The bidegree-(3,5), total-degree-eight bridge is cheaper to parametrize than the
# original degree-ten plane.  Singular returns homogeneous coordinates for
# (K,T,Z) as polynomials in (s,t).
bridge_projective_ring = PolynomialRing(QQ, names=("x", "y", "z"))
bridge_x, bridge_y, bridge_z = bridge_projective_ring.gens()
bridge_projective_polynomial = bridge_projective_ring(
    sum(
        coefficient*bridge_x**exponents[0]*bridge_y**exponents[1]
        * bridge_z**(bridge_degree-sum(exponents))
        for exponents, coefficient in D_eliminant.dict().items()
    )
)
integral_bridge = primitive_polynomial(bridge_projective_polynomial)
bridge_origin_multiplicity = min(
    sum(exponents) for exponents in D_eliminant.dict()
)
bridge_origin_tangent = sum(
    coefficient*bridge_K**exponents[0]*bridge_T**exponents[1]
    for exponents, coefficient in D_eliminant.dict().items()
    if sum(exponents) == bridge_origin_multiplicity
)
bridge_infinity = factor(
    bridge_projective_polynomial(z=0)
)
from sage.libs.singular.function_factory import ff
bridge_local_ring = PolynomialRing(
    QQ, names=("u", "v"), order="negdegrevlex"
)
bridge_u, bridge_v = bridge_local_ring.gens()
bridge_local_y = bridge_local_ring(
    bridge_projective_polynomial(x=bridge_u, y=1, z=bridge_v)
)
bridge_local_x = bridge_local_ring(
    bridge_projective_polynomial(x=1, y=bridge_u, z=bridge_v)
)
bridge_hne_y = ff.hnoether__lib.hnexpansion(bridge_local_y, "ess")
bridge_hne_x = ff.hnoether__lib.hnexpansion(bridge_local_x, "ess")
bridge_delta_y = ZZ(ff.hnoether__lib.delta(bridge_local_y))
bridge_delta_x = ZZ(ff.hnoether__lib.delta(bridge_local_x))
print(
    "Q80SURFACEQQ|stage=bridge_geometry|"
    f"bidegree=3,5|total_degree={bridge_degree}|"
    f"origin_multiplicity={bridge_origin_multiplicity}|"
    f"origin_tangent={factor(bridge_origin_tangent)}|"
    f"infinity={bridge_infinity}|"
    f"infinity_y_conjugacy={bridge_hne_y[1]}|infinity_y_delta={bridge_delta_y}|"
    f"infinity_x_conjugacy={bridge_hne_x[1]}|infinity_x_delta={bridge_delta_x}",
    flush=True,
)
bridge_parameter_homogeneous_ring = PolynomialRing(QQ, names=("s", "r"))
bridge_s, bridge_r = bridge_parameter_homogeneous_ring.gens()


def parse_bridge_polynomial(value):
    value = re.sub(r"s([0-9]+)", r"s^\1", value)
    value = re.sub(r"t([0-9]+)", r"r^\1", value)
    value = re.sub(r"(?<![A-Za-z])t(?![A-Za-z])", "r", value)
    value = re.sub(r"([0-9])([sr])", r"\1*\2", value)
    value = re.sub(r"([sr])(?=[sr])", r"\1*", value)
    return bridge_parameter_homogeneous_ring(value)


bridge_cache_path = Path(arguments.bridge_cache)
bridge_map = None
if bridge_cache_path.is_file():
    bridge_cache = json.loads(bridge_cache_path.read_text())
    if (
        bridge_cache.get("schema")
        == "q80-cm24-qq-auxiliary-bridge-parameter-v1"
        and bridge_cache.get("bridge_equation") == str(D_eliminant)
        and bridge_cache.get("method")
        == "singular_choice4_rational_normal_inversion"
    ):
        bridge_map = tuple(
            bridge_parameter_homogeneous_ring(value)
            for value in bridge_cache["bridge_map"]
        )
        print(
            "Q80SURFACEQQ|stage=bridge_parameter|cache=hit|"
            f"path={bridge_cache_path}",
            flush=True,
        )

if bridge_map is None:
    bridge_library_source = library_path.read_text()
    bridge_old_call = (
        'ideal AI = adjointIdeal(f,list(choice,"rattestyes/firstchecksdone"));'
    )
    bridge_new_call = (
        'ideal AI = adjointIdeal(f,list(4,"rattestyes/firstchecksdone"));'
    )
    if bridge_library_source.count(bridge_old_call) != 1:
        raise RuntimeError("unexpected paraplanecurves.lib parameter pipeline")
    bridge_library_source = bridge_library_source.replace(
        bridge_old_call, bridge_new_call
    )
    with tempfile.TemporaryDirectory(
        prefix="q80-bridge-parameter-"
    ) as bridge_temporary_name:
        bridge_library_path = (
            Path(bridge_temporary_name)/"q80_bridge_parameter.lib"
        )
        bridge_library_path.write_text(bridge_library_source)
        bridge_program = f'''
LIB "{bridge_library_path}";
ring R=0,(x,y,z),dp;
poly f={integral_bridge};
def RP1=paraPlaneCurve(f,"local");
setring RP1;
int parameter_index;
for(parameter_index=1;parameter_index<=size(PARA);parameter_index++)
{{
  "Q80BRIDGE|"+string(parameter_index)+"|"+string(PARA[parameter_index]);
}}
'''
        try:
            bridge_completed = subprocess.run(
                [str(singular_path), "-q"],
                input=bridge_program,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=min(arguments.timeout, 20),
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            raise TimeoutError(
                "Singular bridge parametrization exceeded its 20s hard limit"
            ) from error
    if bridge_completed.returncode:
        print(bridge_completed.stdout, end="", flush=True)
        raise RuntimeError("Singular failed to parametrize the degree-eight bridge")
    bridge_strings = {}
    for line in bridge_completed.stdout.splitlines():
        if line.startswith("Q80BRIDGE|"):
            _, index, polynomial = line.split("|", 2)
            bridge_strings[int(index)] = polynomial
    if tuple(sorted(bridge_strings)) != (1, 2, 3):
        print(bridge_completed.stdout, end="", flush=True)
        raise RuntimeError("Singular returned no complete bridge parameterization")
    bridge_map = tuple(
        parse_bridge_polynomial(bridge_strings[index])
        for index in (1, 2, 3)
    )
    bridge_cache_path.parent.mkdir(parents=True, exist_ok=True)
    bridge_cache_path.write_text(
        json.dumps(
            {
                "schema": "q80-cm24-qq-auxiliary-bridge-parameter-v1",
                "bridge_equation": str(D_eliminant),
                "bridge_map": [str(value) for value in bridge_map],
                "method": "singular_choice4_rational_normal_inversion",
            },
            indent=2,
            sort_keys=True,
        )+"\n"
    )
bridge_common_factor = gcd(bridge_map)
if not bridge_common_factor.is_constant():
    bridge_map = tuple(value//bridge_common_factor for value in bridge_map)
bridge_map_degree = max(value.total_degree() for value in bridge_map)
if bridge_map_degree != bridge_degree:
    raise ArithmeticError(
        f"bridge map degree is {bridge_map_degree}, expected {bridge_degree}"
    )
bridge_projective_residual = bridge_parameter_homogeneous_ring(
    bridge_projective_polynomial(
        x=bridge_map[0], y=bridge_map[1], z=bridge_map[2]
    )
)
if bridge_projective_residual:
    raise ArithmeticError("bridge parameterization misses its plane")

final_parameter_ring = PolynomialRing(QQ, "u")
u = final_parameter_ring.gen()


def primitive_pair(numerator, denominator):
    coefficient_denominator = lcm(
        coefficient.denominator()
        for polynomial in (numerator, denominator)
        for coefficient in polynomial.coefficients()
    )
    numerator = numerator.parent()(coefficient_denominator*numerator)
    denominator = denominator.parent()(coefficient_denominator*denominator)
    content = gcd(
        tuple(
            ZZ(coefficient)
            for polynomial in (numerator, denominator)
            for coefficient in polynomial.coefficients()
        )
    )
    numerator = numerator//content
    denominator = denominator//content
    if denominator.leading_coefficient() < 0:
        numerator = -numerator
        denominator = -denominator
    return numerator, denominator


D_bridge_gcd = gcd(bridge_map[0], bridge_map[2])
T_bridge_gcd = gcd(bridge_map[1], bridge_map[2])
D_bridge_pair = primitive_pair(
    bridge_map[0]//D_bridge_gcd, bridge_map[2]//D_bridge_gcd
)
T_bridge_pair = primitive_pair(
    bridge_map[1]//T_bridge_gcd, bridge_map[2]//T_bridge_gcd
)
D_bridge_affine = tuple(
    final_parameter_ring(value(s=1, r=u)) for value in D_bridge_pair
)
T_bridge_affine = tuple(
    final_parameter_ring(value(s=1, r=u)) for value in T_bridge_pair
)
print(
    "Q80SURFACEQQ|stage=bridge_function_reduction|"
    f"D_degree={max(value.degree() for value in D_bridge_affine)}|"
    f"T_degree={max(value.degree() for value in T_bridge_affine)}|"
    f"D_cancel_degree={D_bridge_gcd.total_degree()}|"
    f"T_cancel_degree={T_bridge_gcd.total_degree()}",
    flush=True,
)

composition_cache_path = Path(arguments.composition_cache)
composition_payload = None
if composition_cache_path.is_file():
    candidate = json.loads(composition_cache_path.read_text())
    if (
        candidate.get("schema") == "q80-cm24-qq-reduced-composition-v1"
        and candidate.get("bridge_equation") == str(D_eliminant)
        and candidate.get("D_numerator") == str(D_bridge_affine[0])
        and candidate.get("D_denominator") == str(D_bridge_affine[1])
        and candidate.get("T_numerator") == str(T_bridge_affine[0])
        and candidate.get("T_denominator") == str(T_bridge_affine[1])
    ):
        composition_payload = candidate

if composition_payload is None:
    # Recover Q from the degree-one subresultant.  Pairwise cancellation of
    # the bridge map first keeps D and T at degrees five and three.
    def evaluate_bridge_numerator(polynomial):
        denominator_degrees = (
            polynomial.degree(polynomial.parent().gen(0)),
            polynomial.degree(polynomial.parent().gen(1)),
        )
        numerator = final_parameter_ring.zero()
        for exponents, coefficient in polynomial.dict().items():
            numerator += (
                coefficient
                * D_bridge_affine[0]**exponents[0]
                * D_bridge_affine[1]**(
                    denominator_degrees[0]-exponents[0]
                )
                * T_bridge_affine[0]**exponents[1]
                * T_bridge_affine[1]**(
                    denominator_degrees[1]-exponents[1]
                )
            )
        return numerator, denominator_degrees

    Q_constant_numerator, Q_constant_denominator_degree = (
        evaluate_bridge_numerator(Q_linear_subresultant[0])
    )
    Q_coefficient_numerator, Q_coefficient_denominator_degree = (
        evaluate_bridge_numerator(Q_linear_subresultant[1])
    )
    if not Q_coefficient_numerator:
        raise ArithmeticError("linear Q subresultant degenerates on the bridge")
    Q_unreduced_numerator = -Q_constant_numerator
    Q_unreduced_denominator = Q_coefficient_numerator
    for coordinate, denominator in enumerate(
        (D_bridge_affine[1], T_bridge_affine[1])
    ):
        degree_difference = (
            Q_coefficient_denominator_degree[coordinate]
            - Q_constant_denominator_degree[coordinate]
        )
        if degree_difference >= 0:
            Q_unreduced_numerator *= denominator**degree_difference
        else:
            Q_unreduced_denominator *= denominator**(-degree_difference)
    composition_payload = {
        "schema": "q80-cm24-qq-reduced-composition-v1",
        "bridge_equation": str(D_eliminant),
        "D_numerator": str(D_bridge_affine[0]),
        "D_denominator": str(D_bridge_affine[1]),
        "T_numerator": str(T_bridge_affine[0]),
        "T_denominator": str(T_bridge_affine[1]),
        "Q_unreduced_numerator": str(Q_unreduced_numerator),
        "Q_unreduced_denominator": str(Q_unreduced_denominator),
    }
    composition_cache_path.parent.mkdir(parents=True, exist_ok=True)
    composition_cache_path.write_text(
        json.dumps(composition_payload, indent=2, sort_keys=True)+"\n"
    )
    composition_status = "created"
else:
    composition_status = "hit"
print(
    "Q80SURFACEQQ|stage=reduced_composition|"
    f"Q_degrees={final_parameter_ring(composition_payload['Q_unreduced_numerator']).degree()},"
    f"{final_parameter_ring(composition_payload['Q_unreduced_denominator']).degree()}|"
    f"cache={composition_status}|path={composition_cache_path}",
    flush=True,
)

reduced_q_path = Path(arguments.reduced_q)


def load_reduced_q():
    if not reduced_q_path.is_file():
        return None
    candidate = json.loads(reduced_q_path.read_text())
    if (
        candidate.get("schema") != "q80-cm24-qq-reduced-Q-v1"
        or candidate.get("input", {}).get("sha256")
        != sha256(composition_cache_path)
        or candidate.get("degrees") != [10, 10]
    ):
        return None
    return (
        final_parameter_ring(candidate["Q_numerator"]),
        final_parameter_ring(candidate["Q_denominator"]),
    )


reduced_q_pair = load_reduced_q()
if reduced_q_pair is None:
    reducer_path = Path(__file__).with_name(
        "reconstruct_q80_reduced_q_modular.sage"
    )
    try:
        reduced_completed = subprocess.run(
            [
                sys.executable,
                str(reducer_path),
                "--input", str(composition_cache_path),
                "--output", str(reduced_q_path),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise TimeoutError("modular Q reconstruction exceeded 30s") from error
    print(reduced_completed.stdout, end="", flush=True)
    if reduced_completed.returncode:
        raise RuntimeError("modular Q reconstruction failed")
    reduced_q_pair = load_reduced_q()
    if reduced_q_pair is None:
        raise ArithmeticError("modular Q artifact failed its provenance check")
    reduced_q_status = "created"
else:
    reduced_q_status = "hit"
Q_reduced_numerator, Q_reduced_denominator = reduced_q_pair
print(
    "Q80SURFACEQQ|stage=reduced_Q|"
    f"cache={reduced_q_status}|path={reduced_q_path}|degrees=10,10",
    flush=True,
)

# A direct degree-fifteen representative has a degree-five common factor.
raw_projective_map = (
    D_bridge_affine[0]*Q_reduced_denominator,
    Q_reduced_numerator*D_bridge_affine[1],
    D_bridge_affine[1]*Q_reduced_denominator,
)
raw_map_common_factor = gcd(raw_projective_map)
X_affine, Y_affine, Z_affine = tuple(
    coordinate//raw_map_common_factor for coordinate in raw_projective_map
)
raw_map_common_factor_degree = raw_map_common_factor.degree()
print(
    "Q80SURFACEQQ|stage=projective_reduction|"
    f"raw_common_factor_degree={raw_map_common_factor_degree}|"
    f"reduced_degrees={X_affine.degree()},{Y_affine.degree()},{Z_affine.degree()}",
    flush=True,
)

D_degrees = tuple(value.degree() for value in D_bridge_affine)
Q_degrees = tuple(value.degree() for value in reduced_q_pair)
map_degree = max(X_affine.degree(), Y_affine.degree(), Z_affine.degree())
if map_degree != plane_degree:
    raise ArithmeticError(
        f"parameter map has degree {map_degree}, expected {plane_degree}"
    )


def evaluate_projective_form_on_affine_map(polynomial):
    result = final_parameter_ring.zero()
    for exponents, coefficient in polynomial.dict().items():
        result += (
            coefficient
            * X_affine**exponents[0]
            * Y_affine**exponents[1]
            * Z_affine**exponents[2]
        )
    return result


high_parameter_numerator = evaluate_projective_form_on_affine_map(
    high_form_old
)
low_parameter_numerator = evaluate_projective_form_on_affine_map(
    low_form_old
)
auxiliary_parameter_residual = (
    T_bridge_affine[1]*high_parameter_numerator
    - T_bridge_affine[0]*low_parameter_numerator
)
if auxiliary_parameter_residual:
    raise ArithmeticError("reduced map misses the auxiliary bridge ratio")

homogeneous_parameter_ring = PolynomialRing(QQ, names=("s", "t"))
s, homogeneous_t = homogeneous_parameter_ring.gens()


def homogenize_parameter(polynomial):
    result = homogeneous_parameter_ring.zero()
    for exponent, coefficient in polynomial.dict().items():
        if isinstance(exponent, tuple):
            exponent = exponent[0]
        result += (
            coefficient*homogeneous_t**exponent*s**(map_degree-exponent)
        )
    return result


homogeneous_map = tuple(
    homogenize_parameter(polynomial)
    for polynomial in (X_affine, Y_affine, Z_affine)
)
common_factor = gcd(homogeneous_map)
if not common_factor.is_constant():
    raise ArithmeticError("homogeneous parameter map has a common factor")
projective_residual = homogeneous_parameter_ring.zero()
for exponents, coefficient in projective_polynomial.dict().items():
    projective_residual += (
        coefficient
        * homogeneous_map[0]**exponents[0]
        * homogeneous_map[1]**exponents[1]
        * homogeneous_map[2]**exponents[2]
    )
if projective_residual:
    raise ArithmeticError("homogeneous map misses the exact projective plane")
plane_residual = projective_residual


def function_record(numerator, denominator, degrees):
    return {
        "value": f"({numerator})/({denominator})",
        "numerator": str(numerator),
        "denominator": str(denominator),
        "degrees": [int(value) for value in degrees],
    }


output_payload = {
    "schema": "q80-cm24-qq-DQ-parameter-v1",
    "scope": "exact_global_birational_parameterization_of_the_reconstructed_DQ_plane",
    "status": "exact_forward_identity_and_birational_degree_formula",
    "slope": "8/87",
    "plane_degree": int(plane_degree),
    "arithmetic_genus": int((plane_degree-1)*(plane_degree-2)//2),
    "geometric_genus": int(0),
    "adjoint_coordinate_change": "z_old=z_new-y",
    "choice4_auxiliary_dimensions": {
        "degree_7": len(choice4_degree_seven_forms),
        "degree_8": len(auxiliary_forms),
        "caveat": (
            "These dimensions are not used as a conductor/genus certificate; "
            "choice 4 is two conditions larger than the exact conductor here."
        ),
    },
    "auxiliary_parameter": {
        "value": f"({high_form_old})/({low_form_old})",
        "high_form_new_coordinates": str(high_form_new),
        "low_form_new_coordinates": str(low_form_new),
        "high_order_at_cm24": int(high_series.valuation()),
        "low_order_at_cm24": int(low_series.valuation()),
        "parameter_at_cm24": "t=0",
    },
    "degree_eight_bridge": {
        "equation": str(D_eliminant),
        "degrees": [
            int(D_eliminant.degree(bridge_K)),
            int(D_eliminant.degree(bridge_T)),
            int(bridge_degree),
        ],
        "homogeneous_map": [str(value) for value in bridge_map],
        "map_degree": int(bridge_map_degree),
        "projective_residual": str(bridge_projective_residual),
        "final_parameter": "u=t/s",
        "auxiliary_T_of_u": (
            f"({T_bridge_affine[0]})/({T_bridge_affine[1]})"
        ),
    },
    "centered_coordinates": plane_payload["centered_coordinates"],
    "functions": {
        "D": function_record(
            D_bridge_affine[0], D_bridge_affine[1], D_degrees
        ),
        "Q": function_record(
            Q_reduced_numerator, Q_reduced_denominator, Q_degrees
        ),
    },
    "homogeneous_map": {
        name: str(value)
        for name, value in zip(("D", "Q", "Z"), homogeneous_map)
    },
    "checks": {
        "plane_residual": str(plane_residual),
        "auxiliary_parameter_residual": str(auxiliary_parameter_residual),
        "projective_residual": str(projective_residual),
        "map_degree": int(map_degree),
        "common_factor": str(common_factor),
        "birational_degree_from_degree_formula": 1,
        "series_coefficients": int(precision),
        "parameter_series_precision": int(parameter_precision),
    },
    "inputs": [
        {"path": str(plane_path), "sha256": plane_hash},
        {"path": str(series_path), "sha256": sha256(series_path)},
        {"path": str(adjoint_cache_path), "sha256": sha256(adjoint_cache_path)},
    ],
}

output_path = Path(arguments.output)
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(output_payload, indent=2, sort_keys=True)+"\n")

print(
    "Q80SURFACEQQ|stage=geometry|plane_degree=10|arithmetic_genus=36|"
    "exact_local_delta_genus=0|choice4_degree7=2|choice4_degree8=11|"
    "auxiliary_osculating_degree=10|"
    f"D_degrees={D_degrees[0]},{D_degrees[1]}|"
    f"Q_degrees={Q_degrees[0]},{Q_degrees[1]}|map_degree={map_degree}|"
    "plane_residual=0|auxiliary_parameter_residual=0|projective_residual=0|"
    f"output={output_path}|status=PASS_EXACT_GENUS0_PARAMETERIZATION",
    flush=True,
)
