#!/usr/bin/env sage -python
"""
Certify the literal horizontal support of the third H3->q6 transported section.

The historical one-point generic-fibre representative is

    43*O + P,   MW(P)=22*(-P1)-P2,

and its exact source-vertical correction is

    -22 E7_1 -33 E7_2 -44 E7_3 -66 E7_4
    -55 E7_5 -44 E7_6 -33 E7_7 -2389 F.

That representative is correct but computationally poor: collapsing the
22 repeated -P1 points into the high-height group-law point P creates both
the large E7 correction and the huge fibre twist.

This audit keeps the literal effective horizontal divisor instead:

    H_lit = 21*O + 22*(-P1) + (-P2),

where (-P2) is the inverse of the pinned frame P2 class and is the sign
represented by the reconstructed Hensel coordinates used elsewhere in the
repo.

It proves in the pinned Neron--Severi lattice that the SAME transported q6
section is

    S3 = H_lit - 46*F.

Thus all E7 and E8 root corrections disappear and the fibre twist drops from
-2389 to -46.

It also records a 44-element generic-fibre basis strategy:
  * the standard balanced basis of L(21 O + 21(-P1));
  * one extra function m1^22-x^11, which has pole 22 at -P1 but its order-22
    pole at O cancels;
  * one chord m2 with simple poles at O and -P2.

The last basis statement is a generic-fibre pole/Riemann--Roch certificate;
resolved endpoint and smooth-collision lattices remain the next step.

Run:
  sage -python ~/Downloads/audit_h92_q6_third_literal_support.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
TRANSPORT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-h3-q6-weyl-section-transport.json"
)
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
OUTPUT = (
    args.output.resolve()
    if args.output
    else ROOT
    / "artifacts/local/elkies-k3"
    / "q6-third-literal-support.json"
)

for path in (FRAME, TRANSPORT, CORE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

transport = json.loads(TRANSPORT.read_text())
assert transport["status"] == "PASS_EXACT_Q6_WEYL_SECTION_TRANSPORT"

frame = load_gram(FRAME)
assert frame.nrows() == frame.ncols() == 17
ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)

F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
roots = matrix(ZZ, [list(root) for root in simple])
root_gram = roots * ns * roots.transpose()
assert abs(root_gram.det()) == 2  # E7 + E8.

# These are exactly the section classes used by the q6 transport certifier.
minus_p1 = vector(
    ZZ,
    [5, 1]
    + [-value for value in (2, 3, 4, 6, 5, 4, 3)]
    + [0] * 8
    + [1, 0],
)
p2_frame = vector(ZZ, [22, 1] + [0] * 16 + [1])

for name, point in (("minus_P1", minus_p1), ("P2_frame", p2_frame)):
    assert point * ns * point == -2, name
    assert point * ns * F == 1, name

# P2 has zero reducible-fibre correction and P2.O=21.  Its inverse is
# therefore obtained by negating its MW projection while preserving O+23F.
minus_p2 = vector(ZZ, [22, 1] + [0] * 16 + [-1])
assert minus_p2 * ns * minus_p2 == -2
assert minus_p2 * ns * F == 1
assert minus_p2 * ns * O == 21
assert all(minus_p2 * ns * root == 0 for root in simple)


def source_shioda(section):
    horizontal = section - O - (section * ns * O + 2) * F
    assert horizontal * ns * F == 0
    assert horizontal * ns * O == 0
    # Explicit identity matrix avoids relying on a helper loaded from CORE.
    from sage.all import identity_matrix
    projection = (
        identity_matrix(QQ, 19)
        - ns * roots.transpose() * root_gram.inverse() * roots
    )
    return vector(QQ, horizontal) * projection


phi1 = source_shioda(minus_p1)
phi2 = source_shioda(p2_frame)
height = matrix(QQ, [
    [-left * ns * right for right in (phi1, phi2)]
    for left in (phi1, phi2)
])
assert height == matrix(QQ, ((QQ(21) / 2, 3), (3, 46)))
assert source_shioda(minus_p2) == -phi2

# Reconstruct exactly the high point used by the old one-point description,
# following certify_h3_q6_weyl_section_transport.sage itself.
third_shioda = 22 * phi1 - phi2
third_height = -QQ(third_shioda * ns * third_shioda)
assert third_height.denominator() == 1
third_generic_q = (
    vector(QQ, O)
    + (third_height / 2) * vector(QQ, F)
    + third_shioda
)
assert all(value in ZZ for value in third_generic_q)
third_generic = vector(ZZ, third_generic_q)
assert third_generic * ns * third_generic == -2
assert third_generic * ns * F == 1

# Rebuild the old vertical correction from its certified source basis.
third = transport["third_vertical_correction"]
expected_basis = [
    "old_E7_1", "old_E7_2", "old_E7_3", "old_E7_4",
    "old_E7_5", "old_E7_6", "old_E7_7",
    "old_E8_1", "old_E8_2", "old_E8_3", "old_E8_4",
    "old_E8_5", "old_E8_6", "old_E8_7", "old_E8_8",
    "old_F",
]
assert third["basis"] == expected_basis
vertical_coordinates = vector(ZZ, third["coordinates"])
vertical_basis = roots.stack(matrix(ZZ, [list(F)]))
assert vertical_basis.nrows() == 16
old_vertical = vertical_coordinates * vertical_basis

# Old representation.
old_horizontal = 43 * O + third_generic
target = old_horizontal + old_vertical
assert target * ns * target == -2
assert target * ns * F == 44

# New literal support.  This is the key identity.
literal_horizontal = 21 * O + 22 * minus_p1 + minus_p2
assert literal_horizontal * ns * F == 44
assert target == literal_horizontal - 46 * F

# All non-fibre vertical correction has disappeared.
literal_vertical_coordinates = vector(ZZ, [0] * 15 + [-46])
assert target == literal_horizontal + literal_vertical_coordinates * vertical_basis

# The actual q6 nef fibre is the first marked divisor O+(-P1)-F.
q6_fibre = O + minus_p1 - F
assert q6_fibre * ns * q6_fibre == 0
assert q6_fibre * ns * F == 2
assert literal_horizontal * ns * q6_fibre == 93
assert (-46 * F) * ns * q6_fibre == -92
assert target * ns * q6_fibre == 1

# Bind the literal support to the reusable compiler contract.
scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
support_certificate = scope["certify_generic_fibre_horizontal_support"](
    ns,
    target,
    F,
    (
        ("old_O", 21, O),
        ("minus_P1", 22, minus_p1),
        ("minus_P2", 1, minus_p2),
    ),
    vertical_support=(),
    fiber_twist=-46,
    expected_old_fiber_degree=44,
)
assert support_certificate["old_fiber_degree"] == 44
assert support_certificate["fiber_twist"] == -46

# Generic-fibre basis:
# L(21O+21A) has the standard 42-element balanced chord basis.
balanced = scope["balanced_marked_chord_power_basis"](
    21, chord_symbol="m1", x_symbol="x"
)
assert len(balanced) == 42

# Two extra functions enlarge
#   21O+21A -> 21O+22A -> 21O+22A+C.
#
# extra_A = m1^22-x^11:
# at A, m1 has a simple pole and x is regular, so the pole order is 22;
# at O, m1^2-x is regular because the order-two leading poles cancel.
# Factoring
#   m1^22-x^11=(m1^2-x)*sum_{j=0}^{10} m1^(20-2j)*x^j
# bounds its O-pole by 20, hence it lies in L(21O+22A).
#
# extra_C = m2 is the marked chord whose poles are O and C=-P2.
generic_basis = list(balanced) + [
    {
        "kind": "extra_minus_P1_pole",
        "formula": "m1^22-x^11",
        "pole_order_at_O_upper_bound": 20,
        "pole_order_at_minus_P1": 22,
        "pole_order_at_minus_P2": 0,
        "proof": (
            "m1^2-x is O-regular; factor a^11-b^11. "
            "At minus_P1, m1 has a simple pole and x is regular."
        ),
    },
    {
        "kind": "minus_P2_chord",
        "formula": "m2",
        "pole_order_at_O": 1,
        "pole_order_at_minus_P1": 0,
        "pole_order_at_minus_P2": 1,
        "proof": "m2 is the chord through +P2, hence its uncancelled pole is -P2.",
    },
]
assert len(generic_basis) == 44

# Independence is certified by nested pole profiles:
# 42 balanced functions span L(21O+21A);
# extra_A has a new order-22 pole at A;
# m2 has a new pole at C.
# Degree 44 on a genus-one curve has h0=44, so the 44 independent functions
# are complete.

print(
    "Q6THIRDLITERAL|"
    "old_horizontal=43O+P|"
    "old_E7=-22,-33,-44,-66,-55,-44,-33|"
    "old_E8=0,0,0,0,0,0,0,0|old_F=-2389|"
    "literal=21O+22(-P1)+(-P2)-46F|"
    "literal_roots=0|status=PASS",
    flush=True,
)
print(
    "Q6THIRDLITERAL_DEGREE|"
    "old_fibre=44|q6_horizontal=93|q6_fibre_twist=-92|"
    "q6_target=1|status=PASS",
    flush=True,
)
print(
    "Q6THIRDLITERAL_BASIS|"
    "balanced=42|extra=m1^22-x^11,m2|dimension=44|"
    "status=PASS_GENERIC_RR",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q6-third-literal-support.v1",
    "status": "PASS_EXACT_Q6_THIRD_LITERAL_SUPPORT",
    "old_one_point_representation": {
        "horizontal": "43*O + P, MW(P)=22*(-P1)-P2",
        "vertical_coordinates": list(map(int, vertical_coordinates)),
        "fiber_twist": int(vertical_coordinates[-1]),
    },
    "literal_representation": {
        "horizontal_support": [
            {"name": "old_O", "multiplicity": 21},
            {"name": "minus_P1", "multiplicity": 22},
            {"name": "minus_P2", "multiplicity": 1},
        ],
        "vertical_root_coordinates": [0] * 15,
        "fiber_twist": -46,
        "identity": "S3=21*O+22*(-P1)+(-P2)-46*F",
        "compiler_certificate": {
            "old_fiber_degree": support_certificate["old_fiber_degree"],
            "fiber_twist": support_certificate["fiber_twist"],
        },
    },
    "q6_degree_balance": {
        "literal_horizontal": 93,
        "fiber_twist": -92,
        "target": 1,
    },
    "generic_rr_basis": {
        "dimension": 44,
        "balanced_subspace": {
            "divisor": "21*O+21*(-P1)",
            "dimension": 42,
            "basis": list(balanced),
        },
        "extra_functions": generic_basis[-2:],
        "completeness": (
            "Nested pole profiles give 44 independent functions; "
            "genus-one Riemann-Roch for degree 44 gives h0=44."
        ),
    },
    "boundary": (
        "This removes the artificial E7/E8 root correction and reduces the "
        "pure fibre twist to -46. It certifies the generic-fibre basis strategy "
        "but does not yet assemble the endpoint/smooth-collision coefficient "
        "lattices or recover explicit q6-child coordinates for the section."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}")
