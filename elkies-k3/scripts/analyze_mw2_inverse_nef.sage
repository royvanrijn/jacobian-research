from sage.all import *
from pathlib import Path
import argparse


parser = argparse.ArgumentParser(
    description="Audit transported inverse-neighbor fibers against MW2 section curves."
)
parser.add_argument("--bound", type=int, default=20)
args = parser.parse_args()

ROOT = Path(__file__).resolve().parents[2]
load(str(ROOT / "elkies-k3" / "scripts" / "verify_mw2_rank17_transport.sage"))

# Save these before the section-class construction reuses Sage's global
# preparser temporaries.
backward_reductions = {
    stage_name: reduce_against_explicit_curves(old_fiber)
    for stage_name, old_fiber, _, _ in backward_classes
}

component_minimum_cache = {}


def cached_component_minimum(component, point, gram, expected):
    key = (component, tuple(point), expected)
    if key not in component_minimum_cache:
        component_minimum_cache[key] = minimized_component_coordinates(point, gram, expected)
    return component_minimum_cache[key]


def section_in_explicit(m, n):
    target_vector = m * target_basis.row(0) + n * target_basis.row(1)
    lift = None
    for projected, representative in projected_cosets:
        difference = target_vector - projected
        if all(QQ(value).denominator() == 1 for value in difference):
            lift = vector(QQ, representative) + difference * C
            break
    assert lift is not None and all(QQ(value).denominator() == 1 for value in lift)
    root_part = vector(QQ, lift) - target_vector * C
    root_coordinates = (root_part * terminal_frame * R.transpose()) * GR_inverse

    minimized = []
    for component, gram in enumerate(
        basis * terminal_frame * basis.transpose() for basis in component_simple_bases
    ):
        left, right = component_bounds[component]
        point = fractional_class(root_coordinates[left:right])
        order = class_order(point, 6)
        expected = QQ(0) if order == 1 else expected_nonzero_norm[component]
        minimized.extend(cached_component_minimum(component, point, gram, expected))
    frame_lift = target_vector * C + vector(QQ, minimized) * R
    assert all(QQ(value).denominator() == 1 for value in frame_lift)
    frame_lift = vector(ZZ, frame_lift)
    zero_intersection = ZZ((frame_lift * terminal_frame * frame_lift - 4) / 2)
    standard_section = vector(
        ZZ, [zero_intersection + 1, 1] + list(frame_lift)
    )
    explicit_section = standard_section * explicit_terminal_basis.inverse()
    assert all(QQ(value).denominator() == 1 for value in explicit_section)
    explicit_section = vector(ZZ, explicit_section)
    assert explicit_section * explicit_gram * explicit_section == -2
    assert explicit_section * explicit_gram * terminal_fiber == 1
    return vector(ZZ, explicit_section), zero_intersection


sections = {}
for m in range(-args.bound, args.bound + 1):
    for n in range(-args.bound, args.bound + 1):
        section, zero_intersection = section_in_explicit(m, n)
        sections[(m, n)] = (section, zero_intersection)

assert sections[(0, 0)][0] == vector(ZZ, [0, 1] + [0] * 17)
assert sections[(1, 0)][0] == vector(ZZ, [0] * 17 + [1, 0])
assert sections[(0, 1)][0] == vector(ZZ, [0] * 18 + [1])

# The first inverse pencil has horizontal group sum 3*P1-4*P2.  Record the
# exact vertical remainder of the natural O+Q pole member; these coordinates
# are the input for a valuation-driven Riemann--Roch construction.
first_divisor = backward_reductions["a5_d4_2a2_a1_mw3"][0]
first_horizontal = sections[(3, -4)][0]
first_zero = sections[(0, 0)][0]
first_vertical = first_divisor - first_zero - first_horizontal
print(
    f"MW2NEF|stage=a5_d4_2a2_a1_mw3|pole_member=O+(3P1-4P2)+V"
    f"|Q={tuple(first_horizontal)}|V={tuple(first_vertical)}"
    f"|V_intersections={tuple(first_vertical * explicit_gram)}",
    flush=True,
)

for stage_name, old_fiber, _, _ in backward_classes:
    component_reduced, component_reflections = backward_reductions[stage_name]
    for variant, divisor in (("raw", old_fiber), ("component_reduced", component_reduced)):
        negative = []
        for (m, n), (section, zero_intersection) in sections.items():
            intersection = ZZ(divisor * explicit_gram * section)
            if intersection < 0:
                negative.append((intersection, m, n, zero_intersection, section))
        negative.sort(key=lambda item: (item[0], item[1], item[2]))
        print(
            f"MW2NEF|stage={stage_name}|variant={variant}|bound={args.bound}"
            f"|component_reflections={len(component_reflections)}"
            f"|degree={divisor * explicit_gram * terminal_fiber}"
            f"|negative_sections={len(negative)}"
            f"|minimum={None if not negative else negative[0][:4]}",
            flush=True,
        )
        for intersection, m, n, zero_intersection, section in negative[:20]:
            print(
                f"MW2NEF|stage={stage_name}|variant={variant}|section=({m},{n})"
                f"|D.S={intersection}|S.O={zero_intersection}"
                f"|class={tuple(section)}",
                flush=True,
            )

print("MW2NEF|status=PASS", flush=True)
