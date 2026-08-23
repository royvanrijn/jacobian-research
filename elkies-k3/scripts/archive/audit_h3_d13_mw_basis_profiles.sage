#!/usr/bin/env sage -python
"""
Certify effective section profiles for the four pinned D13/MW4 generators
and measure their degrees on the preceding q6 fibration.

This uses only the already-certified dominant D13 coordinate chain.  It does
not identify the dominant and component-nef D13 presentations; q6 fibre degree
is invariant under later q6-component Weyl reflections.

For each MW basis vector z=e_i:
  1. compute its exact D13 discriminant class;
  2. derive the minimal D13 local correction (0, 1, or 13/4);
  3. derive P.O from Shioda's height formula;
  4. solve the corresponding D13 root-lattice CVP for an integral effective
     section lift w with norm height+correction;
  5. construct the (-2)-section [a,1,w];
  6. pull it into the q6 simple frame and source H3 NS.

The q24 target (0,-1,1,1) is independently checked against its known witness.

Run:
  sage -python ~/Downloads/audit_h3_d13_mw_basis_profiles.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix,
    lcm, matrix, vector
)

Q24_WITNESS = vector(ZZ, (
    0, 5, 0, 1, 2, 1, 2, 2, 2, 2, 4, 8, 2, 0, -1, 1, 1,
))


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
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def class_order(dual):
    order = ZZ(1)
    for value in dual:
        order = lcm(order, ZZ(QQ(value).denominator()))
    return order


def d13_minimal_correction(dual, root):
    """
    D13 has discriminant group Z/4:
      order 1 -> 0
      order 2 -> vector class, correction 1
      order 4 -> spinor class, correction 13/4.
    Check the quadratic class modulo 2 as an independent guard.
    """
    order = class_order(dual)
    raw = QQ(dual * root * dual)

    def mod_two(value):
        value = QQ(value)
        return value - 2 * (value / 2).floor()

    expected = {
        ZZ(1): QQ(0),
        ZZ(2): QQ(1),
        ZZ(4): QQ(13) / 4,
    }
    assert order in expected
    correction = expected[order]
    assert mod_two(raw) == mod_two(correction)
    return order, correction


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--cvp-cap", type=int, default=4096)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
PINNED = (
    ROOT / "elkies-k3/data/fibrations/"
           "h3_q6_q8_d13_mw4_root_adapted_frame.txt"
)
ORBITS = (
    ROOT / "artifacts/generated-results/"
           "elkies-k3-h3-q6-q8-orbits.json"
)
OUTPUT = (
    args.output.resolve()
    if args.output
    else ROOT / "artifacts/local/elkies-k3/d13-mw-basis-profiles.json"
)

source_frame = load_gram(FRAME)
G = load_gram(PINNED)
data = json.loads(ORBITS.read_text())
assert data["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert G.det() == source_frame.det() == 948

U = matrix(ZZ, ((0, 1), (1, 0)))
source_ns = block_diagonal_matrix(U, -source_frame)
d13_ns = block_diagonal_matrix(U, -G)

# ---------------------------------------------------------------------------
# Certified dominant coordinate chain, identical to the passing invariant
# audit.
# ---------------------------------------------------------------------------

B6 = matrix(ZZ, data["q6"]["neighbor_basis_in_source_ns"])
assert abs(B6.det()) == 1
q6_raw_ns = B6 * source_ns * B6.transpose()
assert q6_raw_ns[:2, :2] == U
q6_raw = -q6_raw_ns[2:, 2:]

root_mw = matrix(ZZ, data["q6"]["root_mw_basis_in_child"])
assert abs(root_mw.det()) == 1
assert (
    root_mw * q6_raw * root_mw.transpose()
    == matrix(ZZ, data["q6"]["root_adapted_gram"])
)

simple14 = matrix(ZZ, data["q8"]["simple_root_change_in_root_block"])
simple_change = block_diagonal_matrix(simple14, identity_matrix(ZZ, 3))
simple_to_raw = simple_change * root_mw
simple_frame = matrix(ZZ, data["q8"]["simple_frame_gram"])
assert simple_to_raw * q6_raw * simple_to_raw.transpose() == simple_frame

Bsimple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), simple_to_raw) * B6
)
q6_simple_ns = block_diagonal_matrix(U, -simple_frame)
assert Bsimple * source_ns * Bsimple.transpose() == q6_simple_ns

hits = data["q8"]["d13_mw4_hits"]
dominant_hit = hits[0]
assert matrix(ZZ, dominant_hit["d13_root_adapted_gram"]) == G

B8 = matrix(ZZ, dominant_hit["neighbor_basis_in_q6_ns"])
raw_d13 = matrix(ZZ, dominant_hit["child_frame"])
assert (
    B8 * q6_simple_ns * B8.transpose()
    == block_diagonal_matrix(U, -raw_d13)
)

A13 = matrix(ZZ, dominant_hit["d13_root_adapted_basis_in_child"])
assert A13 * raw_d13 * A13.transpose() == G

Bpinned_to_simple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * B8
)
assert (
    Bpinned_to_simple * q6_simple_ns * Bpinned_to_simple.transpose()
    == d13_ns
)

Fq6_simple = vector(ZZ, [1, 0] + [0] * 17)
Oq6_simple = vector(ZZ, [-1, 1] + [0] * 17)
Fh3 = vector(ZZ, [1, 0] + [0] * 17)
Oh3 = vector(ZZ, [-1, 1] + [0] * 17)

print(
    "D13MWPROFILE_COORD|q6=PASS|simple=PASS|dominant=PASS|pinned=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# D13 MW lattice data.
# ---------------------------------------------------------------------------

root = G[:13, :13]
coupling = G[:13, 13:]
tail = G[13:, 13:]
H = tail - coupling.transpose() * root.inverse() * coupling

EXPECTED_H = matrix(QQ, [
    [QQ(3)/4, QQ(1)/4, -QQ(1)/4, 0],
    [QQ(1)/4, QQ(11)/4, QQ(1)/4, 1],
    [-QQ(1)/4, QQ(1)/4, QQ(11)/4, -1],
    [0, 1, -1, 46],
])
assert H == EXPECTED_H
assert H.det() == 237
assert root.det() == 4

root_lattice = IntegralLattice(root)

def effective_lift_for_mw(z):
    z = vector(ZZ, z)
    base = vector(ZZ, [0] * 13 + list(z))

    # Orthogonal projection of the quotient lift onto D13 roots.
    pairings = base * G[:, :13]
    dual = vector(QQ, pairings) * root.inverse()
    order, correction = d13_minimal_correction(dual, root)
    height = QQ(z * H * z)
    target_norm = height + correction
    assert target_norm in ZZ and target_norm >= 4 and target_norm % 2 == 0

    # Find the closest D13 root shift. enumerate_close_vectors() is ordered
    # by distance; the first accepted target-norm lift is therefore minimal.
    iterator = root_lattice.enumerate_close_vectors(-dual)
    seen_smaller = []
    chosen = None
    for unused in range(args.cvp_cap):
        shift = vector(ZZ, next(iterator))
        candidate = base + vector(ZZ, list(shift) + [0] * 4)
        norm = ZZ(candidate * G * candidate)
        if norm < target_norm:
            seen_smaller.append((tuple(map(int, shift)), int(norm)))
            continue
        if norm == target_norm:
            chosen = candidate
            break
        # enumeration has passed target distance without a solution
        if norm > target_norm:
            break

    assert not seen_smaller, (
        "found an integral representative below the ADE minimal correction",
        seen_smaller[:10],
    )
    assert chosen is not None, (
        "CVP did not find the predicted minimal representative",
        tuple(z), target_norm,
    )

    norm = ZZ(chosen * G * chosen)
    pole = (norm - 4) // 2
    assert QQ(pole) == (height + correction - 4) / 2
    a = (norm - 2) // 2
    section = vector(ZZ, [a, 1] + list(chosen))
    assert section * d13_ns * section == -2
    assert section * d13_ns * vector(ZZ, [1,0]+[0]*17) == 1
    assert section * d13_ns * vector(ZZ, [-1,1]+[0]*17) == pole

    return {
        "mw": z,
        "height": height,
        "class_order": order,
        "correction": correction,
        "pole": ZZ(pole),
        "lift": chosen,
        "section": section,
    }


profiles = []
for i in range(4):
    z = vector(ZZ, [ZZ(j == i) for j in range(4)])
    profile = effective_lift_for_mw(z)

    section_simple = profile["section"] * Bpinned_to_simple
    section_source = section_simple * Bsimple

    assert section_simple * q6_simple_ns * section_simple == -2
    assert section_source * source_ns * section_source == -2

    q6_degree = ZZ(section_simple * q6_simple_ns * Fq6_simple)
    q6_O = ZZ(section_simple * q6_simple_ns * Oq6_simple)
    h3_degree = ZZ(section_source * source_ns * Fh3)
    h3_O = ZZ(section_source * source_ns * Oh3)

    profile.update({
        "section_simple": section_simple,
        "section_source": section_source,
        "q6_degree": q6_degree,
        "q6_O": q6_O,
        "h3_degree": h3_degree,
        "h3_O": h3_O,
    })
    profiles.append(profile)

    print(
        "D13MWPROFILE|"
        f"generator={i+1}|height={profile['height']}|"
        f"class_order={profile['class_order']}|"
        f"correction={profile['correction']}|PdotO={profile['pole']}|"
        f"q6_degree={q6_degree}|q6_O={q6_O}|"
        f"h3_degree={h3_degree}|h3_O={h3_O}|"
        f"lift={','.join(map(str,profile['lift']))}",
        flush=True,
    )

assert [p["pole"] for p in profiles] == [0, 1, 1, 21]

# ---------------------------------------------------------------------------
# Independent q24 target check.
# ---------------------------------------------------------------------------

z24 = vector(ZZ, (0, -1, 1, 1))
h24 = QQ(z24 * H * z24)
assert h24 == 47

base24 = vector(ZZ, Q24_WITNESS)
pair24 = vector(QQ, base24 * G[:, :13])
dual24 = pair24 * root.inverse()
order24, correction24 = d13_minimal_correction(dual24, root)
assert order24 == 2 and correction24 == 1
assert base24 * G * base24 == 48
pole24 = (48 - 4) // 2
assert pole24 == 22

P24 = vector(ZZ, [23, 1] + list(base24))
P24_simple = P24 * Bpinned_to_simple
P24_source = P24_simple * Bsimple
assert P24 * d13_ns * P24 == -2
assert P24_simple * q6_simple_ns * P24_simple == -2
assert P24_source * source_ns * P24_source == -2

P24_q6_degree = ZZ(P24_simple * q6_simple_ns * Fq6_simple)
assert P24_q6_degree == 49

print(
    "D13MWPROFILE_Q24|mw=0,-1,1,1|height=47|"
    "class_order=2|correction=1|PdotO=22|"
    f"q6_degree={P24_q6_degree}|"
    "relation=-G2+G3+G4",
    flush=True,
)

# The q24 MW vector uses the high fourth generator.  This assertion guards
# against mistakenly searching only in the span of the three low-pole points.
assert z24 == -vector(ZZ, (0,1,0,0)) + vector(ZZ, (0,0,1,0)) + vector(ZZ, (0,0,0,1))

payload = {
    "schema": "elkies-k3.h3-d13-mw-basis-profiles.v1",
    "status": "PASS_EXACT_D13_MW_BASIS_PROFILES",
    "height_gram": [[str(v) for v in row] for row in H.rows()],
    "generators": [
        {
            "index": i + 1,
            "mw_coordinates": list(map(int, p["mw"])),
            "height": str(p["height"]),
            "D13_discriminant_class_order": int(p["class_order"]),
            "D13_local_correction": str(p["correction"]),
            "P_dot_O": int(p["pole"]),
            "effective_frame_lift": list(map(int, p["lift"])),
            "section_pinned_D13_NS": list(map(int, p["section"])),
            "q6_degree": int(p["q6_degree"]),
            "q6_abstract_zero_intersection": int(p["q6_O"]),
            "source_H3_degree": int(p["h3_degree"]),
            "source_H3_zero_intersection": int(p["h3_O"]),
        }
        for i, p in enumerate(profiles)
    ],
    "q24_target": {
        "mw_coordinates": [0,-1,1,1],
        "height": "47",
        "D13_local_correction": "1",
        "P_dot_O": 22,
        "q6_degree": int(P24_q6_degree),
        "MW_relation": "-G2+G3+G4",
    },
    "boundary": (
        "This certifies lattice-effective section classes and pullback degrees. "
        "It does not yet give rational Weierstrass coordinates for the D13 MW generators."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("D13MWPROFILE_RESULT|status=PASS_EXACT_D13_MW_BASIS_PROFILES")
print(f"OUTPUT|{OUTPUT}")
