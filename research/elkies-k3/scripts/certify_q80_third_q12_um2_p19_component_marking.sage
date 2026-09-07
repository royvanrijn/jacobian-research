#!/usr/bin/env sage -python
"""Certify the transported p=19 A5+A3+3A1 component and zero marking."""

import contextlib
import hashlib
import io
import json
import runpy
import sys
import tempfile
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, ZZ, vector


ROOT = Path(__file__).resolve().parents[2]
CHAMBER = ROOT / "elkies-k3/scripts/analyze_q80_second_neighbor_chamber.sage"
COMPILER = ROOT / "elkies-k3/scripts/compile_q80_third_q12_um2_p19_resolved_pencil.sage"
MINIMAL = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-minimal.json"
MAPS = ROOT / "artifacts/generated-results/q80-third-q12-p19-birational-maps.json"
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-component-marking.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


minimal = json.loads(MINIMAL.read_text())
maps = json.loads(MAPS.read_text())
if minimal.get("status") != "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC":
    raise ValueError("minimal child is not certified")
if maps.get("status") != "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_MOD19_QUADRATIC":
    raise ValueError("generic maps are not certified")

# Reuse the canonical lattice/chamber calculation and the immutable resolved
# compiler without copying either implementation.  The compiler writes only
# to a temporary path during this replay.
with contextlib.redirect_stdout(io.StringIO()):
    chamber = runpy.run_path(str(CHAMBER))
with tempfile.TemporaryDirectory() as temporary_directory:
    saved_argv = sys.argv
    sys.argv = [str(COMPILER), "--output", str(Path(temporary_directory) / "resolved.json")]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            compiler = runpy.run_path(str(COMPILER))
    finally:
        sys.argv = saved_argv

D = vector(ZZ, chamber["third_reduced"])
NS = chamber["third_ns"]
old_curves = [(name, vector(ZZ, curve)) for name, curve in chamber["third_curves"]]


def intersection(left, right):
    return ZZ(vector(ZZ, left) * NS * vector(ZZ, right))


curve_by_name = dict(old_curves)
degrees = {name: intersection(D, curve) for name, curve in old_curves}
if intersection(D, D) != 0 or degrees["O"] != 0:
    raise ArithmeticError("wrong new fibre class or old-zero degree")
vertical = [(name, curve) for name, curve in old_curves if degrees[name] == 0]
if len(vertical) != 11:
    raise ArithmeticError("transported vertical old-curve count is not eleven")

vertical_gram = Matrix(ZZ, [[intersection(left, right) for unused, right in vertical] for unused, left in vertical])
unseen = set(range(len(vertical)))
connected = []
while unseen:
    todo = [min(unseen)]
    unseen.remove(todo[0])
    component = []
    while todo:
        index = todo.pop()
        component.append(index)
        for other in list(unseen):
            if vertical_gram[index, other]:
                unseen.remove(other)
                todo.append(other)
    connected.append(tuple(sorted(component)))
connected.sort(key=lambda component: (len(component), component))
if [len(component) for component in connected] != [1, 1, 1, 3, 5]:
    raise ArithmeticError("transported old curves do not have A5+A3+3A1 graph")
if sorted(abs(vertical_gram.matrix_from_rows_and_columns(component, component).det()) for component in connected) != [2, 2, 2, 4, 6]:
    raise ArithmeticError("transported root determinants are wrong")

components_by_size = {}
for component in connected:
    components_by_size.setdefault(len(component), []).append(
        [vertical[index][0] for index in component]
    )
if components_by_size[5] != [["O", "R4", "R6", "R12", "Theta0_2"]]:
    raise ArithmeticError("unexpected transported A5 support")
if components_by_size[3] != [["R8", "R9", "R10"]]:
    raise ArithmeticError("unexpected transported A3 support")
if sorted(label[0] for label in components_by_size[1]) != ["R2", "R3", "R7"]:
    raise ArithmeticError("unexpected transported A1 supports")

# The old D7 fibre has component multiplicities recorded by the affine-root
# calculation.  Among its degree-one components, R5 has multiplicity one and
# R1 multiplicity two.  These are exactly the simple xi=-6 and double xi=3
# infinity branches of the resolved cubic, so the canonical Laurent origin is
# the old component R5.
source_components = chamber["third_components"]
source_affine = chamber["third_affine_data"]
d7_index = next(index for index, component in enumerate(source_components) if len(component) == 7)
d7_component_indices = source_components[d7_index]
d7_multiplicities = source_affine[d7_index][2]
d7_multiplicity = {
    f"R{root_index + 1}": ZZ(d7_multiplicities[position])
    for position, root_index in enumerate(d7_component_indices)
}
degree_one_d7 = {
    name: d7_multiplicity[name]
    for name in d7_multiplicity
    if degrees[name] == 1
}
if degree_one_d7 != {"R1": 2, "R5": 1}:
    raise ArithmeticError("D7 infinity branches do not orient as R1 double/R5 simple")
new_zero_name = "R5"
new_zero = curve_by_name[new_zero_name]

# Complete each displayed finite-A chain to its affine I_n cycle.  For I2,
# the two components meet twice; for n>2 the graph is an ordinary cycle.
root_groups = [components_by_size[5][0], components_by_size[3][0]]
root_groups += [group for group in components_by_size[1]]
cycles = []
for labels in root_groups:
    old_components = [curve_by_name[label] for label in labels]
    missing = D - sum(old_components, vector(ZZ, [0] * len(D)))
    full_components = old_components + [missing]
    size = len(full_components)
    if any(intersection(component, component) != -2 for component in full_components):
        raise ArithmeticError("completed fibre has a non-root component")
    if any(intersection(D, component) != 0 for component in full_components):
        raise ArithmeticError("completed fibre component is not vertical")
    if sum(full_components, vector(ZZ, [0] * len(D))) != D:
        raise ArithmeticError("completed fibre components do not sum to the fibre")
    gram = Matrix(ZZ, [[intersection(left, right) for right in full_components] for left in full_components])
    expected = Matrix(ZZ, size, size)
    for index in range(size):
        expected[index, index] = -2
        expected[index, (index + 1) % size] += 1
        expected[(index + 1) % size, index] += 1
    # The input old-chain order is not necessarily path order, so compare
    # valencies and edge multiplicities rather than literal matrices.
    if sorted(sum(max(ZZ(0), gram[i, j]) for j in range(size)) for i in range(size)) != sorted(
        sum(max(ZZ(0), expected[i, j]) for j in range(size)) for i in range(size)
    ):
        raise ArithmeticError("completed fibre does not have affine-A incidence")
    zero_hits = [intersection(new_zero, component) for component in full_components]
    if sorted(zero_hits) != [0] * (size - 1) + [1]:
        raise ArithmeticError("selected zero does not meet one component of completed fibre")
    identity_index = zero_hits.index(1)
    nonidentity = [component for index, component in enumerate(full_components) if index != identity_index]
    root_gram = -Matrix(ZZ, [[intersection(left, right) for right in nonidentity] for left in nonidentity])
    if root_gram.det() != size:
        raise ArithmeticError("zero-oriented nonidentity root determinant is wrong")
    cycles.append(
        {
            "kodaira": f"I{size}",
            "root": f"A{size - 1}",
            "transported_old_components": labels,
            "missing_component_in_parent_NS": [int(value) for value in missing],
            "identity_component": (
                labels[identity_index] if identity_index < len(labels) else "missing_component"
            ),
            "zero_intersections_in_displayed_order_plus_missing": [int(value) for value in zero_hits],
            "nonidentity_component_classes_in_parent_NS": [
                [int(value) for value in component] for component in nonidentity
            ],
        }
    )

if sorted(cycle["kodaira"] for cycle in cycles) != ["I2", "I2", "I2", "I4", "I6"]:
    raise ArithmeticError("completed transported fibre multiset is wrong")

# Pin the base gauge on the I6 fibre.  At the old zero O, x has pole order 2,
# so x^2 is the unique leading term of both kernel numerators.  Its coefficient
# ratio is the constant value of V=N1/N0 on O.
kernel_numerators = compiler["kernel_numerators"]
old_zero_base = compiler["base_ring"].fraction_field()(
    kernel_numerators[1][2] / kernel_numerators[0][2]
)
if old_zero_base.numerator().degree() != 0 or old_zero_base.denominator().degree() != 0:
    raise ArithmeticError("old zero is not vertical in the resolved base ratio")
old_zero_base = old_zero_base.numerator()[0] / old_zero_base.denominator()[0]

base_finite = GF(19)
modulus = PolynomialRing(base_finite, "m")
m = modulus.gen()
finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)
r = finite.gen()
old_zero_coordinates = list(old_zero_base.list()) + [base_finite.zero(), base_finite.zero()]
if finite(old_zero_coordinates[0]) + finite(old_zero_coordinates[1]) * r != r + 5:
    raise ArithmeticError("unexpected new-base value on the old zero")
V_ring = PolynomialRing(finite, "V")
V = V_ring.gen()


def finite_element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * r


delta = V_ring(
    [
        finite_element(value)
        for value in minimal["minimal_short_weierstrass"]["discriminant_coefficients_low_to_high_1_r"]
    ]
)
i6_factor = next(factor.monic() for factor, exponent in delta.factor() if exponent == 6)
i6_root = -i6_factor[0] / i6_factor[1]
if i6_root != r + 5:
    raise ArithmeticError("transported old zero does not land on the I6 factor")

output = {
    "schema": "elkies-k3.q80-third-q12-component-marking-modp2.v1",
    "status": "PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_MOD19_QUADRATIC",
    "specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"},
    "new_fibre_class_in_parent_D7D5_NS": [int(value) for value in D],
    "old_curve_new_fibre_degrees": {name: int(degree) for name, degree in degrees.items()},
    "transported_vertical_old_curves": [name for name, unused in vertical],
    "transported_root_graph": {
        "components": components_by_size,
        "type": "A5+A3+3A1",
        "rank": 11,
        "determinant": 192,
    },
    "zero_orientation": {
        "selected_new_zero": new_zero_name,
        "equation_branch": "W=infinity, xi=-6 (simple)",
        "double_branch": "R1, xi=3, old-D7 multiplicity 2",
        "simple_branch": "R5, xi=-6, old-D7 multiplicity 1",
        "selected_zero_class_in_parent_D7D5_NS": [int(value) for value in new_zero],
    },
    "completed_fibre_cycles": cycles,
    "base_alignment": {
        "old_zero_section": "O",
        "old_zero_new_base_value": "r + 5",
        "minimal_discriminant_I6_factor": str(i6_factor),
        "match": True,
        "three_I2_fibres": "transported singleton A1 curves are canonical up to permutation of the three I2 base roots",
    },
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in (CHAMBER, COMPILER, MINIMAL, MAPS)
    ],
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "transport of eleven old D7/D5 components and the old zero into the complete A5+A3+3A1 root graph",
            "completion to three I2, one I4, and one I6 component cycles in the parent NS",
            "orientation by the selected simple-branch zero R5",
            "old-zero base alignment with the order-six minimal discriminant factor",
        ],
        "not_proved": [
            "a characteristic-zero Mordell--Weil rank inference from finite-field Shioda--Tate",
            "an ordering of the three A1 components among the three I2 base roots",
            "a second-prime or characteristic-zero alignment",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/certify_q80_third_q12_um2_p19_component_marking.sage",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80THIRDQ12MARKING|zero=R5|old_O_base=r+5|"
    "fibres=I6+I4+3I2|roots=A5+A3+3A1|"
    "status=PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_MOD19_QUADRATIC",
    flush=True,
)
