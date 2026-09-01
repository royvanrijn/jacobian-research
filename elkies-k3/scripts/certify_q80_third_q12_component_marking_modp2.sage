#!/usr/bin/env sage -python
"""Transport the certified Q80 component marking to a common-producer prime."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--integral-marking", type=Path,
    default=RESULTS / "q80-third-q12-p19-component-marking.json",
)
parser.add_argument(
    "--alignment", type=Path,
    default=RESULTS / "q80-third-q12-um2-p19-p61-common-producer-alignment.json",
)
parser.add_argument("--pencil", type=Path, required=True)
parser.add_argument("--minimal", type=Path, required=True)
parser.add_argument("--maps", type=Path, required=True)
parser.add_argument("--batch", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
for name in ("integral_marking", "alignment", "pencil", "minimal", "maps", "batch", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


integral = load(args.integral_marking)
alignment = load(args.alignment)
pencil = load(args.pencil)
minimal = load(args.minimal)
maps = load(args.maps)
batch = load(args.batch)
expected = (
    (integral, "PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_MOD19_QUADRATIC"),
    (alignment, "PASS_EXACT_SECOND_PRIME_COMMON_PRODUCER_ALIGNMENT"),
    (pencil, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"),
    (minimal, "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_COMMON_PRODUCER"),
    (maps, "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_COMMON_PRODUCER"),
    (batch, "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER"),
)
if any(payload.get("status") != status for payload, status in expected):
    raise ValueError("one or more marking inputs are uncertified")
specialization = pencil["specialization"]
if any(payload["specialization"] != specialization for payload in (minimal, maps, batch)):
    raise ValueError("marking input specializations disagree")
prime = int(specialization["prime"])
if prime == 19:
    raise ValueError("use the immutable p=19 marking certifier for the control prime")
common_resolved = alignment["common_signature"]["resolved_divisor"]
observed_resolved = {
    "smith_degrees": pencil["saturated_ambient"]["smith_degrees"],
    "ambient_dimension": pencil["saturated_ambient"]["dimension"],
    "generator_weights": pencil["saturated_ambient"]["generator_weights"],
    "column_labels_generator_degree": pencil["saturated_ambient"]["column_labels_generator_degree"],
    "D7_gate_rank": pencil["resolved_gates"]["D7"]["rank"],
    "combined_gate_rank": pencil["resolved_gates"]["combined_rank"],
    "kernel_dimension": pencil["resolved_gates"]["kernel_dimension"],
    "moving_degrees_T_W_x": pencil["moving_equation"]["degrees_T_W_x"],
}
if observed_resolved != common_resolved:
    raise ValueError("target prime does not have the certified common resolved signature")

modulus_text = specialization["extension_modulus"]
pieces = modulus_text.replace("r^2 + ", "").replace("*r + ", " ").split()
if len(pieces) != 2:
    raise ValueError("cannot parse quadratic extension modulus")
linear, constant = map(int, pieces)
base_finite = GF(prime)
modulus_ring = PolynomialRing(base_finite, "m")
m = modulus_ring.gen()
finite = GF(prime**2, "r", modulus=m**2 + linear * m + constant)
r = finite.gen()


def element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * r


# In the saturated ambient, column (generator 0, degree 2) is the x^2
# leading term at the old zero.  The ratio of its two kernel coefficients is
# therefore the new-base value on O, exactly as in the p=19 integral replay.
labels = [tuple(label) for label in pencil["saturated_ambient"]["column_labels_generator_degree"]]
old_zero_column = labels.index((0, 2))
kernel = pencil["resolved_gates"]["kernel"]
old_zero_base = element(kernel[1][old_zero_column]) / element(kernel[0][old_zero_column])

V_ring = PolynomialRing(finite, "V")
V = V_ring.gen()
leading_old_x = sum(
    element(coordinates) * V**t_degree
    for t_degree, w_degree, x_degree, coordinates
    in pencil["moving_equation"]["terms_T_W_x_coefficient_1_r"]
    if x_degree == 3
)
if leading_old_x.degree() != 2 or leading_old_x(old_zero_base) != 0:
    raise ArithmeticError("old zero does not lie on the resolved moving cubic")
if tuple((factor.degree(), int(exponent)) for factor, exponent in leading_old_x.factor()) != ((1, 2),):
    raise ArithmeticError("old-zero leading coefficient is not the expected double linear factor")

delta = V_ring([
    element(value)
    for value in minimal["minimal_short_weierstrass"][
        "discriminant_coefficients_low_to_high_1_r"
    ]
])
i6_factors = [factor.monic() for factor, exponent in delta.factor() if int(exponent) == 6]
if len(i6_factors) != 1 or i6_factors[0].degree() != 1:
    raise ArithmeticError("minimal child has no unique rational I6 factor")
i6_factor = i6_factors[0]
i6_root = -i6_factor[0] / i6_factor[1]
if i6_root != old_zero_base:
    raise ArithmeticError("transported old zero does not land on the I6 fibre")

sample_summary = batch["training_samples"][0]
sample_path = ROOT / sample_summary["path"]
if sha256(sample_path) != sample_summary["sha256"]:
    raise ArithmeticError("mapped sample hash mismatch")
sample = load(sample_path)
if sample["infinity"]["simple_branch"] != "xi=-6" or sample["infinity"]["double_branch"] != "xi=3":
    raise ArithmeticError("mapped sample has the wrong oriented infinity branches")
if sample["infinity"]["finite_integral_basis_valuations_simple_double"] != [
    [0, 0], [-3, -6], [-2, -4]
]:
    raise ArithmeticError("mapped sample has the wrong infinity normalization")

graph = integral["transported_root_graph"]
cycles = integral["completed_fibre_cycles"]
zero = integral["zero_orientation"]
if graph["type"] != "A5+A3+3A1" or graph["rank"] != 11:
    raise ArithmeticError("integral transported root graph changed")
if sorted(cycle["kodaira"] for cycle in cycles) != ["I2", "I2", "I2", "I4", "I6"]:
    raise ArithmeticError("integral completed fibre cycles changed")
if zero["selected_new_zero"] != "R5" or zero["simple_branch"] != "R5, xi=-6, old-D7 multiplicity 1":
    raise ArithmeticError("integral zero orientation changed")


def coordinates(value):
    values = list(finite(value).list()) + [base_finite.zero(), base_finite.zero()]
    return [int(values[0]), int(values[1])]


output = {
    "schema": "elkies-k3.q80-third-q12-component-marking-modp2.v2",
    "status": "PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_COMMON_PRODUCER",
    "specialization": specialization,
    "transported_root_graph": graph,
    "completed_fibre_cycles": cycles,
    "zero_orientation": zero,
    "base_alignment": {
        "old_zero_section": "O",
        "old_zero_new_base_coefficients_1_r": coordinates(old_zero_base),
        "moving_cubic_old_x3_factorization": [
            [str(factor.monic()), int(exponent)]
            for factor, exponent in leading_old_x.factor()
        ],
        "minimal_discriminant_I6_factor": str(i6_factor),
        "minimal_discriminant_I6_root_coefficients_1_r": coordinates(i6_root),
        "match": True,
        "three_I2_fibres": "canonical up to permutation of the three transported A1 curves",
    },
    "mapped_branch_replay": {
        "simple_branch": "R5, xi=-6, old-D7 multiplicity 1",
        "double_branch": "R1, xi=3, old-D7 multiplicity 2",
        "integral_basis_valuations_simple_double": [[0, 0], [-3, -6], [-2, -4]],
        "sample": {"path": str(sample_path.relative_to(ROOT)), "sha256": sha256(sample_path)},
    },
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in (
            args.integral_marking, args.alignment, args.pencil,
            args.minimal, args.maps, args.batch,
        )
    ],
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "transport of the integral A5+A3+3A1 component graph to the independently aligned prime",
            "orientation by the old R5 simple branch as the new zero",
            "old-zero base alignment with the unique order-six minimal discriminant factor",
            "explicit mapped-branch replay on the target-prime curve",
        ],
        "not_proved": [
            "a characteristic-zero Mordell--Weil rank inference from finite-field Shioda--Tate",
            "an ordering of the three A1 components among the three I2 base roots",
            "a characteristic-zero coefficient reconstruction",
        ],
    },
    "reproduce": (
        "sage -python elkies-k3/scripts/certify_q80_third_q12_component_marking_modp2.sage "
        f"--pencil {args.pencil} --minimal {args.minimal} --maps {args.maps} "
        f"--batch {args.batch} --output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONMARKING|prime={prime}|zero=R5|"
    f"old_O_base={coordinates(old_zero_base)}|fibres=I6+I4+3I2|"
    "roots=A5+A3+3A1|"
    "status=PASS_EXACT_TRANSPORTED_THIRD_Q12_COMPONENT_MARKING_COMMON_PRODUCER"
)
