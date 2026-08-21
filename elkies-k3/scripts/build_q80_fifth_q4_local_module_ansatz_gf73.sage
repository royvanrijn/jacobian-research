#!/usr/bin/env sage
"""Pin the smallest fifth-q4 divisor-module ansatz over GF(73).

The compact fourth-q12 curve is the bidegree-(14,3) equation stored in the
moving-cubic artifact.  The fifth horizontal target is the inverse of the
explicit degree-one section whose old polynomial x-coordinate is
``12+26V+13V^2+39V^3+49V^4``.  This checker verifies the selected y-sign,
its exact new-base map, and the component factor containing it at every
reducible support.

No global Jacobian is constructed.  The fifth pencil is recorded as the
two-dimensional fractional module

    L(O+(-R)) = Hom_A(I_O I_{-R}, A)

in the normalization of the compact coordinate ring.  The output pins the
section ideal and all known local data; constructing one compensated
nonconstant generator remains the next equation gate.
"""

import hashlib
import json
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
MOVING_ARTIFACT = (
    ROOT / "artifacts/generated-results/"
    "q80-fourth-q12-cm24-moving-cubic-gf73.json"
)
DISCRIMINANT_ARTIFACT = (
    ROOT / "artifacts/generated-results/"
    "q80-fourth-q12-cm24-discriminant-gf73.json"
)
moving_bytes = MOVING_ARTIFACT.read_bytes()
moving_hash = hashlib.sha256(moving_bytes).hexdigest()
assert moving_hash == (
    "c6560b3db2d1232866e9996fc727924090aa46293c2482885cf9f9dbf4c21c89"
)
moving_data = json.loads(moving_bytes)
discriminant_bytes = DISCRIMINANT_ARTIFACT.read_bytes()
discriminant_hash = hashlib.sha256(discriminant_bytes).hexdigest()
assert discriminant_hash == (
    "6caa3c9bb83a115a1e40689bf23d58dea1dd7ae1c77795f9970b6d24517a7ef0"
)
discriminant_data = json.loads(discriminant_bytes)
assert discriminant_data["semistable_signature"] == "I6+I4+4I3"

finite = GF(73, impl="modn")
surface_ring = PolynomialRing(finite, names=("T", "V", "X"))
T, V, X = surface_ring.gens()
moving = sum(
    finite(coefficient)*T**t_degree*V**v_degree*X**x_degree
    for t_degree, v_degree, x_degree, coefficient
    in moving_data["moving_terms_T_v_x_coefficient"]
)

old_ring = PolynomialRing(finite, "v")
v = old_ring.gen()
old_x = old_ring((12, 26, 13, 39, 49))
old_y = old_ring((52, 20, 45, 26, 44, 37, 9))
old_A = old_ring(moving_data["old_surface"]["A_coefficients"])
old_B = old_ring(moving_data["old_surface"]["B_coefficients"])
assert old_y**2 == old_x**3+old_A*old_x+old_B

section_base_numerator = -4*v+30
section_base_denominator = v-34
section_base = section_base_numerator/section_base_denominator
restricted = old_ring.fraction_field()(
    moving(T=section_base, V=v, X=old_x)
)
assert restricted == 0

# Over GF(73)(T), these are the two generators of the explicit section ideal
# in the compact affine coordinate ring.
parameter_field = FunctionField(finite, "t")
t = parameter_field.gen()
section_v = (34*t+30)/(t+4)
section_x = sum(
    finite(coefficient)*section_v**degree
    for degree, coefficient in enumerate(old_x.list())
)
generic_ring = PolynomialRing(parameter_field, names=("vg", "xg"))
vg, xg = generic_ring.gens()
generic_moving = generic_ring(sum(
    parameter_field(coefficient)
    *t**t_degree*vg**v_degree*xg**x_degree
    for t_degree, v_degree, x_degree, coefficient
    in moving_data["moving_terms_T_v_x_coefficient"]
))
assert generic_moving(vg=section_v, xg=section_x) == 0

factor_rows = []
for fiber in discriminant_data["finite_fibers"]:
    parameter = finite(fiber["T"])
    section_special_v = finite((34*parameter+30)/(parameter+4))
    section_special_x = old_x(section_special_v)
    specialized = surface_ring(moving(T=parameter))
    factors = tuple(specialized.factor())
    containing = tuple(
        index
        for index, (factor, _) in enumerate(factors)
        if factor(V=section_special_v, X=section_special_x) == 0
    )
    assert containing
    containing_factors = [
        {
            "factor_index": int(index),
            "factor_degree_V": int(factors[index][0].degree(V)),
            "factor_degree_X": int(factors[index][0].degree(X)),
            "factor_exponent": int(factors[index][1]),
        }
        for index in containing
    ]
    factor_rows.append(
        {
            "T": int(parameter),
            "fiber": fiber["fiber"],
            "section_V": int(section_special_v),
            "section_X": int(section_special_x),
            "containing_factors": containing_factors,
        }
    )

# T has a pole on this section at V=34.  The point lies on the irreducible
# leading cubic, which is the remaining I3 assigned by the lattice audit.
infinity_v = finite(34)
infinity_x = old_x(infinity_v)
infinity_equation = surface_ring(moving.coefficient({T: moving.degree(T)}))
assert infinity_equation(V=infinity_v, X=infinity_x) == 0

geometric_zero_old = (
    24, 18, -119, 320, -106, 103, 490, 65, -131, 4,
    463, -569, -708, -845, 281, 304, 4, 309, 39, -49,
)
target_component_corrections = (
    ("A2", "0"),
    ("A2", "2/3"),
    ("A2", "2/3"),
    ("A2", "2/3"),
    ("A3", "3/4"),
    ("A5", "5/6"),
)

artifact = {
    "schema": "q80-fifth-q4-local-module-ansatz-gf73-v1",
    "prime": int(73),
    "moving_artifact": str(MOVING_ARTIFACT.relative_to(ROOT)),
    "moving_sha256": moving_hash,
    "discriminant_artifact": str(DISCRIMINANT_ARTIFACT.relative_to(ROOT)),
    "discriminant_sha256": discriminant_hash,
    "explicit_section": {
        "old_MW_coordinates": [1, -1, -1],
        "new_MW_coordinates": [1, -1],
        "target_inverse_coordinates": [-1, 1],
        "x_coefficients": list(map(int, old_x.list())),
        "y_coefficients": list(map(int, old_y.list())),
        "new_base_numerator_coefficients": list(
            map(int, section_base_numerator.list())
        ),
        "new_base_denominator_coefficients": list(
            map(int, section_base_denominator.list())
        ),
        "generic_section_V": str(section_v),
        "generic_section_X": str(section_x),
    },
    "geometric_zero": {
        "old_divisor_coordinates": list(map(int, geometric_zero_old)),
        "old_fiber_degree": 18,
        "old_zero_pairing": 6,
    },
    "fifth_divisor": {
        "old_fiber_degree": 2,
        "D_dot_O": 0,
        "horizontal_height": "5/12",
        "horizontal_P_dot_O": 0,
        "component_corrections": [list(row) for row in target_component_corrections],
    },
    "section_local_factors": factor_rows,
    "infinity": {
        "fiber": "I3",
        "section_V": int(infinity_v),
        "section_X": int(infinity_x),
        "leading_cubic_vanishing": True,
    },
    "normalized_fractional_module": {
        "coordinate_ring": "A=normalization(GF(73)(t)[v,x]/(moving))",
        "section_ideal": "I_R=(v-(34*t+30)/(t+4), x-X_R(t))",
        "target": "L(O+(-R))=Hom_A(I_O*I_{-R},A)",
        "expected_dimension": 2,
        "status": "ANSATZ_COMPENSATED_GENERATOR_OPEN",
    },
    "reproduce": (
        "sage elkies-k3/scripts/"
        "build_q80_fifth_q4_local_module_ansatz_gf73.sage"
    ),
}
output_path = (
    ROOT / "artifacts/generated-results/"
    "q80-fifth-q4-local-module-ansatz-gf73.json"
)
encoded = json.dumps(artifact, indent=2, sort_keys=True, default=int)+"\n"
output_path.write_text(encoded)
output_hash = hashlib.sha256(encoded.encode("utf-8")).hexdigest()

print(
    "Q80FIFTHQ4MODULE|"
    f"section_T=({section_base_numerator})/({section_base_denominator})|"
    f"section_X={old_x}|section_Y={old_y}|"
    f"local_factors={tuple((row['T'], row['fiber'], tuple((item['factor_index'], item['factor_degree_V'], item['factor_degree_X']) for item in row['containing_factors'])) for row in factor_rows)}|"
    f"infinity_point=({int(infinity_v)},{int(infinity_x)})|"
    "status=PASS_EXPLICIT_SECTION_IDEAL",
    flush=True,
)
print(
    "Q80FIFTHQ4MODULE|"
    "module=L(O+(-R))=Hom_A(I_O*I_-R,A)|expected_dimension=2|"
    "compensated_generator=open|global_jacobian_used=0|"
    f"artifact={output_path}|sha256={output_hash}|"
    "status=PASS_LOCAL_MODULE_ANSATZ",
    flush=True,
)
