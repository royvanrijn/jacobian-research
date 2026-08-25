#!/usr/bin/env sage -python
"""List exact A11 MW markings compatible with the pointed-opposite section.

The companion QQ construction proves the section has height 4/3, P.O=0 and
I12 correction 8/3.  This checker exhaustively enumerates the selected A11
MW lattice at that height and retains exactly those vectors with the same
closest-root correction.  It also records their cosets modulo the eighteen
known equation-shell sections.  No equation-coordinate identification is
asserted unless the profile is unique.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import IntegralLattice, QQ, ZZ, lcm, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-pointed-opposite-mw-candidates.json",
)
args = parser.parse_args()

POINT = LOCAL / "q24-a11-pointed-opposite-section-qq.json"
FRAME = LOCAL / (
    "q24-downstream-lift/d12-c10a-zero-q6-frames/"
    "q6-o0064-r11-n132-d12-ad4a027cb197.txt"
)
BRIDGE = LOCAL / "q24-a11-target-coset-bridge.json"
MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
INPUTS = (POINT, FRAME, BRIDGE, MISMATCH)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


point = json.loads(POINT.read_text())
bridge = json.loads(BRIDGE.read_text())
mismatch = json.loads(MISMATCH.read_text())
assert point["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_SECTION_QQ"
assert bridge["status"] == "PASS_EXACT_A11_TARGET_COSET_BRIDGE"
assert mismatch["status"] == "REJECT_A11_QUINTIC_BRIDGE_ZERO_MISMATCH"

target_height = QQ(point["height_profile"]["height"])
target_correction = QQ(point["height_profile"]["local_correction"])
target_pole = ZZ(point["height_profile"]["P_dot_O"])
frame = load_matrix(FRAME)
root_rank = 11
root = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
tail = frame[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling
denominator = lcm(value.denominator() for value in height.list())
height_integral = matrix(ZZ, denominator * height)
enumeration = pari(height_integral).qfminim(ZZ(target_height * denominator), flag=2)
vectors = {tuple([0] * 6)}
for column in matrix(ZZ, enumeration[2]).columns():
    vectors.add(tuple(column))
    vectors.add(tuple(-column))

root_lattice = IntegralLattice(root)


def profile(values):
    z = vector(ZZ, values)
    h = QQ(z * height * z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * 6)
        norm = QQ(lifted * frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
    correction = minimum - h
    pole = (h + correction - 4) / 2
    return h, correction, pole


shell = [
    vector(ZZ, row)
    for row in bridge["exact_identity_shell"]["MW_vectors_in_equation_order"]
]
shell_module = matrix(ZZ, [list(row) for row in shell]).row_module()
target = vector(ZZ, bridge["selected_bridge"]["mw"])
close_mw = vector(
    ZZ,
    mismatch["correct_selected_R3_transport"]["close_P24"]["A11_MW_Abel_Jacobi"],
)
assert close_mw[-1] == target[-1] == 1


def integral_word(generators, wanted):
    generator_matrix = matrix(ZZ, [list(row) for row in generators])
    hnf, transform = generator_matrix.hermite_form(
        include_zero_rows=False, transformation=True
    )
    if hnf.nrows() != hnf.ncols():
        return None
    basis_coefficients = hnf.solve_left(wanted)
    if not all(value in ZZ for value in basis_coefficients):
        return None
    coefficients = vector(ZZ, basis_coefficients) * transform
    assert coefficients * generator_matrix == wanted
    return coefficients


def shell_word(wanted):
    wanted = vector(ZZ, wanted)
    if wanted[-1] != 0:
        return None
    best = None
    target5 = vector(ZZ, wanted[:5])
    for indices in itertools.combinations(range(len(shell)), 5):
        basis = matrix(ZZ, [list(shell[index][:5]) for index in indices])
        if not basis.det():
            continue
        local = basis.solve_left(target5)
        if not all(value in ZZ for value in local):
            continue
        coefficients = vector(ZZ, [0] * len(shell))
        for index, value in zip(indices, local):
            coefficients[index] = ZZ(value)
        assert coefficients * matrix(ZZ, [list(row) for row in shell]) == wanted
        score = (
            max(abs(value) for value in coefficients),
            sum(abs(value) for value in coefficients),
            tuple(coefficients),
        )
        if best is None or score < best[0]:
            best = (score, coefficients)
    return None if best is None else best[1]

candidates = []
for values in sorted(vectors):
    h, correction, pole = profile(values)
    if (h, correction, pole) != (target_height, target_correction, target_pole):
        continue
    z = vector(ZZ, values)
    saturated = matrix(ZZ, [list(row[:5]) for row in shell] + [list(z[:5])]).row_module()
    saturation_index = (
        abs(ZZ(saturated.basis_matrix().det())) if saturated.rank() == 5 else None
    )
    short_word = None
    for pointed_coefficient in sorted(range(-5, 6), key=lambda value: (abs(value), value)):
        coefficients = shell_word(target - close_mw - pointed_coefficient * z)
        if coefficients is not None:
            short_word = (pointed_coefficient, coefficients)
            break
    candidates.append(
        {
            "mw": [int(value) for value in z],
            "in_exact_identity_shell_lattice": z in shell_module,
            "sixth_coordinate": int(z[-1]),
            "first_five_index_after_adjoining": (
                None if saturation_index is None else int(saturation_index)
            ),
            "bridge_with_close_P24": (
                None
                if short_word is None
                else {
                    "shell_coefficients": [int(value) for value in short_word[1]],
                    "pointed_opposite_coefficient": int(short_word[0]),
                    "close_P24_coefficient": 1,
                    "verified_MW_sum": [int(value) for value in target],
                }
            ),
        }
    )

payload = {
    "schema": "elkies-k3.h3-a11-pointed-opposite-mw-candidates.v1",
    "status": "PASS_EXACT_A11_POINTED_OPPOSITE_MW_PROFILE_ENUMERATION",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "profile": {
        "height": str(target_height),
        "local_correction": str(target_correction),
        "P_dot_O": int(target_pole),
    },
    "candidate_count": len(candidates),
    "candidates": candidates,
    "unique_profile_identification": len(candidates) == 1,
    "close_P24": {
        "A11_degree": int(
            mismatch["correct_selected_R3_transport"]["close_P24"]["A11_degree"]
        ),
        "A11_MW_Abel_Jacobi": [int(value) for value in close_mw],
    },
    "proof_boundary": (
        "Complete selected-lattice enumeration at the exact height/correction "
        "profile. If multiple candidates remain, an equation-side good-reduction "
        "fingerprint is still required to identify the constructed QQ point."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11POINTEDMW|height={}|correction={}|PO={}|candidates={}|sixth={}|status={}".format(
        target_height,
        target_correction,
        target_pole,
        len(candidates),
        ",".join(str(row["sixth_coordinate"]) for row in candidates),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
